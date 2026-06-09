/* test_phase4_predict.c
 * Phase 4 verification: run cm_FastCalibrate on each CM in a file and print
 * the predicted ECM parameters (lambda, mu_extrap, mu_orig) per mode.
 *
 * Usage: ./test_phase4_predict <cm_file>
 *
 * Output format (one block per CM):
 *   # CM N: <name>  clen=<clen>
 *   ECMLC  <lambda>  <mu_extrap>  <mu_orig>
 *   ECMGC  <lambda>  <mu_extrap>  <mu_orig>
 *   ECMLI  <lambda>  <mu_extrap>  <mu_orig>
 *   ECMGI  <lambda>  <mu_extrap>  <mu_orig>
 *
 * Matches the output of: python3 fast_cmcalibrate.py <cm_file>
 */
#include "esl_config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_alphabet.h"

#include "infernal.h"
#include "cm_fast_calibrate.h"

/* Infernal mode index → mode name */
static const char *mode_names[EXP_NMODES] = {
  "ECMGC",   /* EXP_CM_GC = 0 */
  "ECMGI",   /* EXP_CM_GI = 1 */
  "ECMLC",   /* EXP_CM_LC = 2 */
  "ECMLI",   /* EXP_CM_LI = 3 */
};

int
main(int argc, char **argv)
{
  char          errbuf[eslERRBUFSIZE];
  CM_FILE      *cmfp = NULL;
  CM_t         *cm   = NULL;
  ESL_ALPHABET *abc  = NULL;
  int           ncm  = 0;
  int           hstatus, status;

  if (argc != 2) {
    fprintf(stderr, "usage: test_phase4_predict <cm_file>\n");
    return 1;
  }

  if ((status = cm_file_Open(argv[1], NULL, FALSE, &cmfp, errbuf)) != eslOK) {
    fprintf(stderr, "Error opening CM file: %s\n", errbuf);
    return 1;
  }

  while ((hstatus = cm_file_Read(cmfp, TRUE, &abc, &cm)) == eslOK)
    {
      ncm++;
      printf("# CM %d: %s  clen=%d\n", ncm,
             cm->name ? cm->name : "(unnamed)", cm->clen);

      if ((status = cm_FastCalibrate(cm)) != eslOK) {
        fprintf(stderr, "cm_FastCalibrate failed for %s (status=%d)\n",
                cm->name ? cm->name : "(unnamed)", status);
        FreeCM(cm); cm = NULL;
        continue;
      }

      int i;
      for (i = 0; i < EXP_NMODES; i++)
        printf("%-6s  %.10f  %.10f  %.10f  %.10g\n",
               mode_names[i],
               cm->expA[i]->lambda,
               cm->expA[i]->mu_extrap,
               cm->expA[i]->mu_orig,
               cm->expA[i]->dbsize > 0
                 ? (double) cm->expA[i]->nrandhits / (double) cm->expA[i]->dbsize
                 : 0.0);

      FreeCM(cm);
      cm = NULL;
    }

  if (hstatus != eslEOF) {
    fprintf(stderr, "Error reading CM file at CM %d\n", ncm + 1);
    cm_file_Close(cmfp);
    esl_alphabet_Destroy(abc);
    return 1;
  }

  cm_file_Close(cmfp);
  esl_alphabet_Destroy(abc);
  return 0;
}

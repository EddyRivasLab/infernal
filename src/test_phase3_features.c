/* test_phase3_features.c
 * Phase 3 verification: extract all 27 features from a CM file and print them.
 *
 * Usage: ./test_phase3_features <cm_file>
 *
 * For each CM in the file, prints one line per feature:
 *   <feature_name>  <value_with_10_decimal_places>
 * prefixed with a CM header line showing the CM name.
 *
 * Output format is identical to test_phase2_features but covers all
 * FAST_CAL_NFEAT (27) features including Phase 3 topo features.
 */
#include "esl_config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_alphabet.h"

#include "infernal.h"
#include "cm_fast_calibrate.h"

int
main(int argc, char **argv)
{
  char          errbuf[eslERRBUFSIZE];
  CM_FILE      *cmfp = NULL;
  CM_t         *cm   = NULL;
  ESL_ALPHABET *abc  = NULL;
  double        feats[FAST_CAL_NFEAT];
  int           i;
  int           ncm  = 0;
  int           hstatus;
  int           status;

  if (argc != 2) {
    fprintf(stderr, "usage: test_phase3_features <cm_file>\n");
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

      if ((status = cm_FastCalibrate_ExtractFeatures(cm, feats)) != eslOK) {
        fprintf(stderr, "Error extracting features from CM %s\n",
                cm->name ? cm->name : "(unnamed)");
        FreeCM(cm);
        cm = NULL;
        continue;
      }

      for (i = 0; i < FAST_CAL_NFEAT; i++)
        printf("%-30s %.10f\n", cm_FastCalibrate_FeatureName(i), feats[i]);

      FreeCM(cm);
      cm = NULL;
    }

  if (hstatus != eslEOF) {
    fprintf(stderr, "Error reading CM file at CM %d: unexpected status %d\n",
            ncm + 1, hstatus);
    cm_file_Close(cmfp);
    esl_alphabet_Destroy(abc);
    return 1;
  }

  cm_file_Close(cmfp);
  esl_alphabet_Destroy(abc);
  return 0;
}

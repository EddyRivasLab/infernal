/* test_phase1_print_models.c
 * Phase 1 verification: call cm_FastCalibrate_PrintModels() and exit.
 * Build with: see Makefile (or manual build command in REPORT)
 * Run with:  ./test_phase1_print_models 2> phase1_dump.txt
 */
#include "esl_config.h"
#include <stdio.h>
#include "cm_fast_calibrate.h"

int main(void)
{
  int status;
  status = cm_FastCalibrate_PrintModels(stderr);
  cm_FastCalibrateCleanup();
  return (status == 0 /* eslOK */ ) ? 0 : 1;
}

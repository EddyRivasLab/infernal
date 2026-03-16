/* test-p7band.c
 * EPN, Mon Mar 16 2026
 * 
 * Test program for MSV-derived P7 HMM banding.
 * 
 * Tests the full pipeline:
 *   1. Run p7_GMSV() to get full MSV matrix
 *   2. Traceback with my_p7_GTraceMSV() to extract alignment
 *   3. Prune alignment with prune_i2k()
 *   4. Derive bands with p7_pins2bands()
 *   5. Report statistics
 * 
 * Eventually will add:
 *   6. Run banded p7_GForward/Backward and compare to unbanded
 */

#include <esl_config.h>
#include <p7_config.h>
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "infernal.h"

static ESL_OPTIONS options[] = {
  /* name           type         default   env  range   toggles   reqs   incomp     help                                                      docgroup*/
  { "-h",           eslARG_NONE,   FALSE,  NULL, NULL,      NULL,  NULL,  NULL,      "show brief help on version and usage",                        0 },
  { "--pad",        eslARG_INT,    "3",    NULL, "n>=0",    NULL,  NULL,  NULL,      "set band padding to <n> (creates bands of width 2*<n>+1)",    0 },
  { "--minscore",   eslARG_REAL,   "0.0",  NULL, NULL,      NULL,  NULL,  NULL,      "minimum match emission score for alignment position",        0 },
  { "--minlen",     eslARG_INT,    "1",    NULL, "n>=1",    NULL,  NULL,  NULL,      "minimum n-mer length",                                        0 },
  { "--minend",     eslARG_INT,    "0",    NULL, "n>=0",    NULL,  NULL,  NULL,      "minimum distance from end",                                   0 },
  { "--minmprob",   eslARG_REAL,   "0.0",  NULL, "0<=x<=1", NULL,  NULL,  NULL,      "minimum match state probability",                             0 },
  { "--minmcprob",  eslARG_REAL,   "0.0",  NULL, "0<=x<=1", NULL,  NULL,  NULL,      "minimum cumulative match probability for n-mer",              0 },
  { "--maxiprob",   eslARG_REAL,   "1.0",  NULL, "0<=x<=1", NULL,  NULL,  NULL,      "maximum insert state probability",                            0 },
  { "--maxilprob",  eslARG_REAL,   "1.0",  NULL, "0<=x<=1", NULL,  NULL,  NULL,      "maximum insert state probability to left",                    0 },
  { "--phi",        eslARG_NONE,   FALSE,  NULL, NULL,      NULL,  NULL,  NULL,      "calculate and use phi (occupancy) probabilities for pruning", 0 },
  { "--dump",       eslARG_NONE,   FALSE,  NULL, NULL,      NULL,  NULL,  NULL,      "dump the bands (i, i2k, kmin, kmax) to stdout",              0 },
  { "--verbose",    eslARG_NONE,   FALSE,  NULL, NULL,      NULL,  NULL,  NULL,      "verbose output",                                              0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <cmfile> <seqfile>";
static char banner[] = "test MSV-derived P7 HMM banding";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *cmfile  = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  int             status;
  char            errbuf[eslERRBUFSIZE];
  ESL_ALPHABET   *abc     = NULL;
  CM_FILE        *cmfp    = NULL;
  CM_t           *cm      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  ESL_SQ         *sq      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_GMX         *gx      = NULL;
  P7_BG          *bg      = NULL;
  P7_TRACE       *p7tr    = NULL;
  double        **phi     = NULL;
  int            *i2k     = NULL;
  int            *kmin    = NULL;
  int            *kmax    = NULL;
  int             ncells;
  int             i, k;
  int             n_aligned, n_gaps;
  float           bandwidth_avg, bandwidth_min, bandwidth_max;
  
  /* Get command line options */
  int    pad       = esl_opt_GetInteger(go, "--pad");
  float  minscore  = esl_opt_GetReal   (go, "--minscore");
  int    minlen    = esl_opt_GetInteger(go, "--minlen");
  int    minend    = esl_opt_GetInteger(go, "--minend");
  float  minmprob  = esl_opt_GetReal   (go, "--minmprob");
  float  minmcprob = esl_opt_GetReal   (go, "--minmcprob");
  float  maxiprob  = esl_opt_GetReal   (go, "--maxiprob");
  float  maxilprob = esl_opt_GetReal   (go, "--maxilprob");
  int    do_phi    = esl_opt_GetBoolean(go, "--phi");
  int    do_dump   = esl_opt_GetBoolean(go, "--dump");
  int    be_verbose= esl_opt_GetBoolean(go, "--verbose");

  /*********************************************** 
   * Read CM file
   ***********************************************/
  if (be_verbose) printf("# Reading CM from %s...\n", cmfile);
  
  if ((status = cm_file_Open(cmfile, NULL, FALSE, &cmfp, errbuf)) != eslOK) 
    cm_Fail("Failed to open covariance model save file %s\n%s\n", cmfile, errbuf);
  if ((status = cm_file_Read(cmfp, TRUE, &abc, &cm)) != eslOK)
    cm_Fail("Failed to read CM from %s\n", cmfile);
  cm_file_Close(cmfp);
  
  if (be_verbose) printf("# Read CM: %s (%d consensus positions)\n", cm->name, cm->clen);

  /* Configure CM - this will create P7 HMMs if they don't exist */
  if ((status = cm_Configure(cm, errbuf, -1)) != eslOK)
    cm_Fail("Failed to configure CM\n%s\n", errbuf);

  /*********************************************** 
   * Read sequence file
   ***********************************************/
  if (be_verbose) printf("# Reading sequence from %s...\n", seqfile);
  
  status = esl_sqfile_Open(seqfile, eslSQFILE_UNKNOWN, NULL, &sqfp);
  if (status != eslOK) 
    cm_Fail("Failed to open sequence file %s\n", seqfile);
  
  sq = esl_sq_CreateDigital(abc);
  status = esl_sqio_Read(sqfp, sq);
  if (status != eslOK) 
    cm_Fail("Failed to read sequence from %s\n", seqfile);
  
  esl_sqfile_Close(sqfp);
  
  if (be_verbose) printf("# Read sequence: %s (length %d)\n", sq->name, (int)sq->n);

  /*********************************************** 
   * Setup for banding
   ***********************************************/
  
  /* Create P7 objects */
  gm   = p7_profile_Create(cm->mlp7->M, abc);
  gx   = p7_gmx_Create(cm->mlp7->M, sq->n);
  bg   = p7_bg_Create(abc);
  p7tr = p7_trace_Create();
  
  /* Convert CM's P7 HMM to a profile */
  if ((status = p7_ProfileConfig(cm->mlp7, bg, gm, sq->n, p7_LOCAL)) != eslOK)
    cm_Fail("Failed to configure P7 profile\n");
  
  /* Allocate phi (occupancy probabilities) - required by p7_Seq2Bands() */
  /* Set to defaults that won't cause any pruning unless --phi is used */
  if (be_verbose && do_phi) printf("# Calculating phi (occupancy probabilities)...\n");
  
  ESL_ALLOC(phi, sizeof(double *) * (cm->mlp7->M + 1));
  for (k = 0; k <= cm->mlp7->M; k++) {
    ESL_ALLOC(phi[k], sizeof(double) * 3);
    if (do_phi) {
      /* TODO: Calculate actual phi values from Forward/Backward */
      phi[k][HMMMATCH]  = 1.0;  /* dummy values for now */
      phi[k][HMMINSERT] = 0.0;
      phi[k][HMMDELETE] = 0.0;
    } else {
      /* Safe defaults: won't trigger pruning in prune_i2k() */
      phi[k][HMMMATCH]  = 1.0;  /* high match prob = won't prune */
      phi[k][HMMINSERT] = 0.0;  /* low insert prob = won't prune */
      phi[k][HMMDELETE] = 0.0;
    }
  }

  /*********************************************** 
   * Derive bands from MSV
   ***********************************************/
  
  if (be_verbose) {
    printf("#\n");
    printf("# Deriving bands from MSV trace:\n");
    printf("#   Padding:                %d (band width = %d)\n", pad, 2*pad + 1);
    printf("#   Min match score:        %.3f\n", minscore);
    printf("#   Min n-mer length:       %d\n", minlen);
    printf("#   Min distance from end:  %d\n", minend);
    if (do_phi) {
      printf("#   Min match probability:  %.3f\n", minmprob);
      printf("#   Min cumul match prob:   %.3f\n", minmcprob);
      printf("#   Max insert probability: %.3f\n", maxiprob);
      printf("#   Max insert prob (left): %.3f\n", maxilprob);
    }
    printf("#\n");
  }
  
  /* Call the main banding function */
  status = p7_Seq2Bands(cm, errbuf, gm, gx, bg, p7tr, sq->dsq, sq->n,
                        phi, minscore, minlen, minend, 
                        minmprob, minmcprob, maxiprob, maxilprob, pad,
                        &i2k, &kmin, &kmax, &ncells);
  
  if (status == eslEINCOMPAT) {
    printf("# WARNING: MSV trace was discontiguous - all alignments removed\n");
  } else if (status != eslOK) {
    cm_Fail("p7_Seq2Bands() failed\n%s\n", errbuf);
  }

  /*********************************************** 
   * Calculate and display statistics
   ***********************************************/
  
  /* Count aligned positions and gaps */
  n_aligned = 0;
  n_gaps = 0;
  for (i = 1; i <= sq->n; i++) {
    if (i2k[i] > 0) n_aligned++;
    else n_gaps++;
  }
  
  /* Calculate average bandwidth */
  bandwidth_avg = 0.0;
  bandwidth_min = sq->n;
  bandwidth_max = 0;
  for (i = 1; i <= sq->n; i++) {
    float bw = kmax[i] - kmin[i] + 1;
    bandwidth_avg += bw;
    if (bw < bandwidth_min) bandwidth_min = bw;
    if (bw > bandwidth_max) bandwidth_max = bw;
  }
  bandwidth_avg /= sq->n;
  
  printf("#\n");
  printf("# Results:\n");
  printf("#\n");
  printf("# Sequence:            %s\n", sq->name);
  printf("# Length:              %d\n", (int)sq->n);
  printf("# Model:               %s\n", cm->name);
  printf("# Model length:        %d\n", cm->mlp7->M);
  printf("#\n");
  printf("# MSV trace alignment:\n");
  printf("#   Aligned positions: %d (%.1f%%)\n", n_aligned, 100.0 * n_aligned / sq->n);
  printf("#   Gap positions:     %d (%.1f%%)\n", n_gaps, 100.0 * n_gaps / sq->n);
  printf("#\n");
  printf("# Band statistics:\n");
  printf("#   Padding:           %d\n", pad);
  printf("#   Cells in bands:    %d\n", ncells);
  printf("#   Full matrix cells: %d\n", (int)sq->n * cm->mlp7->M);
  printf("#   Fraction:          %.4f\n", (float)ncells / ((float)sq->n * cm->mlp7->M));
  printf("#   Speedup potential: %.2fx\n", ((float)sq->n * cm->mlp7->M) / ncells);
  printf("#\n");
  printf("# Bandwidth per position:\n");
  printf("#   Average:           %.1f\n", bandwidth_avg);
  printf("#   Minimum:           %.0f\n", bandwidth_min);
  printf("#   Maximum:           %.0f\n", bandwidth_max);
  printf("#\n");

  /*********************************************** 
   * Optionally dump bands
   ***********************************************/
  
  if (do_dump) {
    printf("# Band dump (i=position, i2k=aligned_node, kmin=band_min, kmax=band_max):\n");
    printf("#%8s %9s %9s %9s %9s\n", "i", "i2k[i]", "kmin[i]", "kmax[i]", "width");
    for (i = 1; i <= sq->n; i++) {
      printf(" %9d %9d %9d %9d %9d\n", 
             i, i2k[i], kmin[i], kmax[i], kmax[i] - kmin[i] + 1);
    }
  }

  /*********************************************** 
   * Cleanup and exit
   ***********************************************/
  
  if (i2k)  free(i2k);
  if (kmin) free(kmin);
  if (kmax) free(kmax);
  if (phi) {
    for (k = 0; k <= cm->mlp7->M; k++) free(phi[k]);
    free(phi);
  }
  
  p7_trace_Destroy(p7tr);
  p7_bg_Destroy(bg);
  p7_gmx_Destroy(gx);
  p7_profile_Destroy(gm);
  esl_sq_Destroy(sq);
  FreeCM(cm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  
  return 0;

 ERROR:
  cm_Fail("Memory allocation error\n");
  return 1;
}

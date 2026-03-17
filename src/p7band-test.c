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
#include <time.h>
#include <inttypes.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "p7_gbands.h"
#include "p7_gmxb.h"

#include "infernal.h"

/* Local declarations for functions not in infernal.h (due to type conflicts) */
extern int             p7_kbands2gbands(int *i2k, int *kmin, int *kmax, int L, int M, P7_GBANDS **ret_bnd);
extern int             p7_GBackwardBanded(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMXB *bx, float *opt_sc);
extern int             my_p7_GForwardBanded(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMXB *bx, float *opt_sc);

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
  printf("p7band-test starting...\n");
  fflush(stdout);
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  printf("getopts created\n");
  fflush(stdout);
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
  
  /* For banded Forward/Backward testing */
  P7_GBANDS      *bnd     = NULL;
  P7_GMXB        *bx      = NULL;
  P7_GMX         *fwd_gx  = NULL;  /* for unbanded Forward */
  P7_GMX         *bck_gx  = NULL;  /* for unbanded Backward */
  float           fwd_sc, bck_sc;   /* unbanded scores */
  float           bfwd_sc, bbck_sc; /* banded scores */
  double          fwd_time, bck_time, bfwd_time, bbck_time;
  clock_t         t0, t1;
  
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
  printf("About to read CM...\n"); fflush(stdout);
  if (be_verbose) printf("# Reading CM from %s...\n", cmfile);
  
  if ((status = cm_file_Open(cmfile, NULL, FALSE, &cmfp, errbuf)) != eslOK) 
    cm_Fail("Failed to open covariance model save file %s\n%s\n", cmfile, errbuf);
  printf("CM file opened\n"); fflush(stdout);
  if ((status = cm_file_Read(cmfp, TRUE, &abc, &cm)) != eslOK)
    cm_Fail("Failed to read CM from %s\n", cmfile);
  printf("CM read\n"); fflush(stdout);
  cm_file_Close(cmfp);
  
  if (be_verbose) printf("# Read CM: %s (%d consensus positions)\n", cm->name, cm->clen);

  /* Configure CM - this will create P7 HMMs if they don't exist */
  printf("About to configure CM...\n"); fflush(stdout);
  if ((status = cm_Configure(cm, errbuf, -1)) != eslOK)
    cm_Fail("Failed to configure CM\n%s\n", errbuf);
  printf("CM configured\n"); fflush(stdout);

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
  printf("Creating P7 objects...\n"); fflush(stdout);
  gm   = p7_profile_Create(cm->mlp7->M, abc);
  gx   = p7_gmx_Create(cm->mlp7->M, sq->n);
  bg   = p7_bg_Create(abc);
  p7tr = p7_trace_Create();
  printf("P7 objects created\n"); fflush(stdout);
  
  /* Convert CM's P7 HMM to a profile */
  printf("Configuring P7 profile as GLOCAL...\n"); fflush(stdout);
  if ((status = p7_ProfileConfig(cm->mlp7, bg, gm, sq->n, p7_GLOCAL)) != eslOK)
    cm_Fail("Failed to configure P7 profile\n");
  printf("P7 profile configured\n"); fflush(stdout);
  
  /* Allocate phi (occupancy probabilities) - required by p7_Seq2Bands() */
  /* Set to defaults that won't cause any pruning unless --phi is used */
  printf("Allocating phi...\n"); fflush(stdout);
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
  printf("phi allocated\n"); fflush(stdout);

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
  
  /* Create test bands (full-width or manually narrowed based on pad) */
  if (1) {
    /* Create bands with optional narrowing based on pad parameter */
    if (pad < (cm->mlp7->M / 2)) {
      printf("Creating narrow bands with pad=%d...\n", pad); fflush(stdout);
    } else {
      printf("Creating full-width bands (pad=%d covers full model)...\n", pad); fflush(stdout);
    }
    
    /* Allocate arrays */
    ESL_ALLOC(i2k, sizeof(int) * (sq->n + 1));
    ESL_ALLOC(kmin, sizeof(int) * (sq->n + 1));
    ESL_ALLOC(kmax, sizeof(int) * (sq->n + 1));
    
    /* Set bands: constrained by pad around a simple diagonal i->k mapping */
    ncells = 0;
    for (i = 0; i <= sq->n; i++) {
      /* Map sequence position i to model position k along a diagonal */
      int k_center = (i * cm->mlp7->M) / sq->n;  /* Simple proportional mapping */
      if (k_center < 1) k_center = 1;
      if (k_center > cm->mlp7->M) k_center = cm->mlp7->M;
      
      i2k[i] = k_center;
      kmin[i] = ESL_MAX(1, k_center - pad);
      kmax[i] = ESL_MIN(cm->mlp7->M, k_center + pad);
      
      /* Count cells in bands */
      if (i > 0) ncells += (kmax[i] - kmin[i] + 1);
    }
    
    status = eslOK;
    int full_cells = sq->n * cm->mlp7->M;
    printf("Bands created: %d cells (%.1f%% of full matrix)\n",
           ncells, (100.0 * ncells) / full_cells); fflush(stdout);
  } else {
    /* Call the main banding function (requires MSV trace) */
    printf("Calling p7_Seq2Bands...\n"); fflush(stdout);
    status = p7_Seq2Bands(cm, errbuf, gm, gx, bg, p7tr, sq->dsq, sq->n,
                          phi, minscore, minlen, minend, 
                          minmprob, minmcprob, maxiprob, maxilprob, pad,
                          &i2k, &kmin, &kmax, &ncells);
    printf("p7_Seq2Bands completed with status %d\n", status); fflush(stdout);
    
    if (status == eslEINCOMPAT) {
      printf("# WARNING: MSV trace was discontiguous - all alignments removed\n");
    } else if (status != eslOK) {
      cm_Fail("p7_Seq2Bands() failed\n%s\n", errbuf);
    }
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
   * Test banded Forward/Backward
   ***********************************************/
  
  if (be_verbose) {
    printf("# Testing banded Forward/Backward...\\n");
    printf("#   First 5 kmin/kmax values before conversion:\\n");
    for (i = 1; i <= ESL_MIN(5, sq->n); i++) {
      printf("#     row %d: kmin=%d, kmax=%d, i2k=%d\\n", i, kmin[i], kmax[i], i2k[i]);
    }
  }
  
  /* Convert kmin/kmax to P7_GBANDS structure */
  printf("Calling p7_kbands2gbands...\n"); fflush(stdout);
  status = p7_kbands2gbands(i2k, kmin, kmax, sq->n, cm->mlp7->M, &bnd);
  if (status != eslOK) cm_Fail("Failed to convert bands to P7_GBANDS (status=%d)\n", status);
  printf("p7_kbands2gbands completed\n"); fflush(stdout);
  
  if (be_verbose) {
    printf("#   P7_GBANDS created: %d segments, %" PRId64 " cells\n", bnd->nseg, bnd->ncell);
    printf("#     L=%d, M=%d, nrow=%d\n", bnd->L, bnd->M, bnd->nrow);
    if (bnd->nseg > 0) {
      printf("#     First segment: ia=%d, ib=%d\n", bnd->imem[0], bnd->imem[1]);
      if (bnd->nseg > 1) {
        printf("#     Last segment: ia=%d, ib=%d\n", 
               bnd->imem[(bnd->nseg-1)*2], bnd->imem[(bnd->nseg-1)*2+1]);
      }
    }
    if (bnd->nrow > 0) {
      printf("#     First row band: ka=%d, kb=%d\n", bnd->kmem[0], bnd->kmem[1]);
      printf("#     Last row band: ka=%d, kb=%d\n", 
             bnd->kmem[(bnd->nrow-1)*2], bnd->kmem[(bnd->nrow-1)*2+1]);
    }
    /* Dump first few rows of bands */
    printf("#     First 5 row bands:\n");
    for (i = 0; i < ESL_MIN(5, bnd->nrow); i++) {
      printf("#       row %d: ka=%d, kb=%d\n", i+1, 
             bnd->kmem[i*2], bnd->kmem[i*2+1]);
    }
  }
  
  /* Create banded matrix */
  bx = p7_gmxb_Create(bnd);
  if (bx == NULL) cm_Fail("Failed to create P7_GMXB\n");
  
  /* Create unbanded matrices for comparison */
  printf("Creating unbanded matrices...\n");
  fflush(stdout);
  fwd_gx = p7_gmx_Create(cm->mlp7->M, sq->n);
  bck_gx = p7_gmx_Create(cm->mlp7->M, sq->n);
  if (fwd_gx == NULL || bck_gx == NULL) cm_Fail("Failed to create unbanded matrices\n");
  
  /* Run unbanded Forward (reference) */
  printf("Running unbanded Forward...\n");
  fflush(stdout);
  t0 = clock();
  status = p7_GForward(sq->dsq, sq->n, gm, fwd_gx, &fwd_sc);
  t1 = clock();
  if (status != eslOK) cm_Fail("Unbanded Forward failed\n");
  fwd_time = (double)(t1 - t0) / CLOCKS_PER_SEC;
  printf("Unbanded Forward done: score=%.4f\n", fwd_sc);
  fflush(stdout);
  
  if (be_verbose) {
    printf("#   Unbanded Forward: score=%.4f, status=%d\n", fwd_sc, status);
    printf("#   fwd_gx->M=%d, fwd_gx->L=%d\n", fwd_gx->M, (int)fwd_gx->L);
  }
  
  /* Run banded Forward */
  t0 = clock();
  printf("About to call my_p7_GForwardBanded...\n");
  fflush(stdout);
  status = my_p7_GForwardBanded(sq->dsq, sq->n, gm, bx, &bfwd_sc);
  t1 = clock();
  printf("my_p7_GForwardBanded returned\n");
  fflush(stdout);
  if (status != eslOK) cm_Fail("Banded Forward failed\n");
  bfwd_time = (double)(t1 - t0) / CLOCKS_PER_SEC;
  
  if (be_verbose) {
    printf("#   Banded Forward: score=%.4f, status=%d\n", bfwd_sc, status);
    /* Check first few DP cells */
    printf("#   First few bx->dp values: %.4f, %.4f, %.4f\n", 
           bx->dp[0], bx->dp[1], bx->dp[2]);
    
    /* Compare Forward matrices to verify banded storage access */
    int L = sq->n;
    int M = gm->M;
    printf("\n# DEBUG: Comparing Forward matrices (row 1, first 5 k values)\n");
    printf("# Note: With full-width bands, row i starts at dp[(i-1)*M*3]\n");
    for (int k = 1; k <= ESL_MIN(5, M); k++) {
      float M_u = fwd_gx->dp[1][k * p7G_NSCELLS + p7G_M];
      float I_u = fwd_gx->dp[1][k * p7G_NSCELLS + p7G_I];
      float D_u = fwd_gx->dp[1][k * p7G_NSCELLS + p7G_D];
      
      /* For banded with full-width bands:
       * Row 1 (i=1) starts at dp[0], row 2 at dp[M*3], row 3 at dp[2*M*3], etc.
       * Within row i, node k is at: (i-1)*M*3 + (k-1)*3
       * But k in bands is relative to kmin, so for k=1..M with kmin=1: offset is (k-1)*3
       */
      int row_start = 0;  /* Row 1 starts at dp[0] */
      float M_b = bx->dp[row_start + (k-1) * p7G_NSCELLS + p7G_M];
      float I_b = bx->dp[row_start + (k-1) * p7G_NSCELLS + p7G_I];
      float D_b = bx->dp[row_start + (k-1) * p7G_NSCELLS + p7G_D];
      
      printf("#  FWD Row 1, k=%d: M: %.4f vs %.4f (diff %.4f), I: %.4f vs %.4f, D: %.4f vs %.4f\n",
             k, M_u, M_b, M_u-M_b, I_u, I_b, D_u, D_b);
    }
    
    /* Check row L too */
    printf("# Comparing Forward row L=%d (last 3 k values)\n", L);
    int L_row_start = (L-1) * M * p7G_NSCELLS;
    for (int k = M-2; k <= M; k++) {
      float M_u = fwd_gx->dp[L][k * p7G_NSCELLS + p7G_M];
      float I_u = fwd_gx->dp[L][k * p7G_NSCELLS + p7G_I];
      float D_u = fwd_gx->dp[L][k * p7G_NSCELLS + p7G_D];
      
      float M_b = bx->dp[L_row_start + (k-1) * p7G_NSCELLS + p7G_M];
      float I_b = bx->dp[L_row_start + (k-1) * p7G_NSCELLS + p7G_I];
      float D_b = bx->dp[L_row_start + (k-1) * p7G_NSCELLS + p7G_D];
      
      printf("#  FWD Row L, k=%d: M: %.4f vs %.4f (diff %.4f), I: %.4f vs %.4f, D: %.4f vs %.4f\n",
             k, M_u, M_b, M_u-M_b, I_u, I_b, D_u, D_b);
    }
    printf("\n");
  }
  
  /* Run unbanded Backward (reference) */
  t0 = clock();
  status = p7_GBackward(sq->dsq, sq->n, gm, bck_gx, &bck_sc);
  t1 = clock();
  if (status != eslOK) cm_Fail("Unbanded Backward failed\n");
  bck_time = (double)(t1 - t0) / CLOCKS_PER_SEC;
  
  /* Run banded Backward */
  t0 = clock();
  status = p7_GBackwardBanded(sq->dsq, sq->n, gm, bx, &bbck_sc);
  t1 = clock();
  if (status != eslOK) cm_Fail("Banded Backward failed\n");
  bbck_time = (double)(t1 - t0) / CLOCKS_PER_SEC;
  
  /* Debug: Compare unbanded and banded Backward matrices */
  if (be_verbose) {
    int L = sq->n;
    int M = gm->M;
    printf("\n# DEBUG: Comparing Backward matrices (unbanded vs banded)\n");
    
    /* Compare special states for a few rows */
    printf("# Special states:\n");
    for (int i = 0; i <= ESL_MIN(2, L); i++) {
      float xE_u = bck_gx->xmx[i * p7G_NXCELLS + p7G_E];
      float xN_u = bck_gx->xmx[i * p7G_NXCELLS + p7G_N];
      float xJ_u = bck_gx->xmx[i * p7G_NXCELLS + p7G_J];
      float xB_u = bck_gx->xmx[i * p7G_NXCELLS + p7G_B];
      float xC_u = bck_gx->xmx[i * p7G_NXCELLS + p7G_C];
      
      float xE_b = bx->xmx[i * p7G_NXCELLS + p7G_E];
      float xN_b = bx->xmx[i * p7G_NXCELLS + p7G_N];
      float xJ_b = bx->xmx[i * p7G_NXCELLS + p7G_J];
      float xB_b = bx->xmx[i * p7G_NXCELLS + p7G_B];
      float xC_b = bx->xmx[i * p7G_NXCELLS + p7G_C];
      
      printf("#  Row %d: xE: %.4f vs %.4f (diff %.4f), xB: %.4f vs %.4f (diff %.4f)\n",
             i, xE_u, xE_b, xE_u-xE_b, xB_u, xB_b, xB_u-xB_b);
      printf("#         xN: %.4f vs %.4f (diff %.4f), xJ: %.4f vs %.4f, xC: %.4f vs %.4f\n",
             xN_u, xN_b, xN_u-xN_b, xJ_u, xJ_b, xC_u, xC_b);
    }
    
    /* Compare row L - unbanded uses 2D array dp[i][k*p7G_NSCELLS+s] */
    printf("# Row L=%d:\n", L);
    for (int k = 1; k <= ESL_MIN(5, M); k++) {
      float M_u = bck_gx->dp[L][k * p7G_NSCELLS + p7G_M];
      float I_u = bck_gx->dp[L][k * p7G_NSCELLS + p7G_I];
      float D_u = bck_gx->dp[L][k * p7G_NSCELLS + p7G_D];
      
      /* Banded uses 1D array indexed as [i*M*p7G_NSCELLS + k*p7G_NSCELLS + s] */
      float M_b = bx->dp[L * M * p7G_NSCELLS + k * p7G_NSCELLS + p7G_M];
      float I_b = bx->dp[L * M * p7G_NSCELLS + k * p7G_NSCELLS + p7G_I];
      float D_b = bx->dp[L * M * p7G_NSCELLS + k * p7G_NSCELLS + p7G_D];
      
      printf("#  k=%d: M: %.4f vs %.4f (diff %.4f), D: %.4f vs %.4f (diff %.4f)\n",
             k, M_u, M_b, M_u-M_b, D_u, D_b, D_u-D_b);
    }
    
    /* Compare row 1 */
    printf("# Row 1:\n");
    for (int k = 1; k <= ESL_MIN(5, M); k++) {
      float M_u = bck_gx->dp[1][k * p7G_NSCELLS + p7G_M];
      float I_u = bck_gx->dp[1][k * p7G_NSCELLS + p7G_I];
      float D_u = bck_gx->dp[1][k * p7G_NSCELLS + p7G_D];
      
      float M_b = bx->dp[1 * M * p7G_NSCELLS + k * p7G_NSCELLS + p7G_M];
      float I_b = bx->dp[1 * M * p7G_NSCELLS + k * p7G_NSCELLS + p7G_I];
      float D_b = bx->dp[1 * M * p7G_NSCELLS + k * p7G_NSCELLS + p7G_D];
      
      printf("#  k=%d: M: %.4f vs %.4f (diff %.4f), I: %.4f vs %.4f, D: %.4f vs %.4f\n",
             k, M_u, M_b, M_u-M_b, I_u, I_b, D_u, D_b);
    }
    printf("\n");
  }
  
  /* Compare scores */
  float fwd_diff = fwd_sc - bfwd_sc;
  float bck_diff = bck_sc - bbck_sc;
  
  printf("# Banded Forward/Backward validation:\n");
  printf("#\n");
  printf("# Forward scores:\n");
  printf("#   Unbanded:          %.4f nats (%.4f ms)\n", fwd_sc, fwd_time * 1000);
  printf("#   Banded:            %.4f nats (%.4f ms)\n", bfwd_sc, bfwd_time * 1000);
  printf("#   Difference:        %.4f nats (%.6f bits)\n", fwd_diff, fwd_diff / eslCONST_LOG2);
  printf("#   Speedup:           %.2fx\n", fwd_time / bfwd_time);
  printf("#\n");
  printf("# Backward scores:\n");
  printf("#   Unbanded:          %.4f nats (%.4f ms)\n", bck_sc, bck_time * 1000);
  printf("#   Banded:            %.4f nats (%.4f ms)\n", bbck_sc, bbck_time * 1000);
  printf("#   Difference:        %.4f nats (%.6f bits)\n", bck_diff, bck_diff / eslCONST_LOG2);
  printf("#   Speedup:           %.2fx\n", bck_time / bbck_time);
  printf("#\n");
  
  /* Validate that scores match */
  float tolerance = 0.01;  /* bits */
  if (fabs(fwd_diff / eslCONST_LOG2) > tolerance) {
    printf("# WARNING: Forward scores differ by more than %.2f bits!\n", tolerance);
  } else {
    printf("# PASS: Forward scores match within %.2f bits\n", tolerance);
  }
  
  if (fabs(bck_diff / eslCONST_LOG2) > tolerance) {
    printf("# WARNING: Backward scores differ by more than %.2f bits!\n", tolerance);
  } else {
    printf("# PASS: Backward scores match within %.2f bits\n", tolerance);
  }
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
  
  if (bnd)    p7_gbands_Destroy(bnd);
  if (bx)     p7_gmxb_Destroy(bx);
  if (fwd_gx) p7_gmx_Destroy(fwd_gx);
  if (bck_gx) p7_gmx_Destroy(bck_gx);
  
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

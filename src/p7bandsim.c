/* p7bandsim.c
 * EPN, Wed Apr 9 2026
 *
 * Testbed for empirically deriving per-state p7 band pads.
 *
 * Algorithm:
 *   1. Read CM, configure (load fp7/mlp7).
 *   2. For i in 1..N:
 *        a. EmitParsetree(cm, ...) -> (cm_tr, esq)
 *        b. p7_Seq2BandsVit(esq, gm, pad=0) -> kmin/kmax (per-i bands derived
 *           from Viterbi trace of seq alone, no knowledge of CM truth).
 *        c. For each emitting CM state v in cm_tr, get true HMM node k via
 *           cm->cp9map->cs2hn[v][0/1] and the residue position(s) i (and j
 *           for MATPs). Compute band deficit:
 *             if kmin[i] <= k <= kmax[i]: deficit = 0
 *             else: deficit = max(kmin[i] - k, k - kmax[i])
 *        d. Bin per HMM node k: record deficit.
 *   3. After N samples, for each HMM node k:
 *        pad[k] = X-percentile of observed deficits (e.g., 95th).
 *   4. Output: M lines "k pad[k] n_observations".
 *
 * Usage: p7bandsim [options] <cmfile>
 */

#include <esl_config.h>
#include <p7_config.h>
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
#include "esl_vectorops.h"

#include "hmmer.h"
#include "infernal.h"

static ESL_OPTIONS options[] = {
  /* name           type       default  env  range   toggles reqs incomp     help                                 docgroup */
  { "-h",           eslARG_NONE,  FALSE, NULL, NULL,    NULL, NULL, NULL,    "show brief help",                     0 },
  { "-N",           eslARG_INT,   "1000",NULL, "n>0",   NULL, NULL, NULL,    "number of parsetrees to emit",        0 },
  { "-X",           eslARG_REAL,  "0.95",NULL, "0<x<=1",NULL, NULL, NULL,    "target percentile (e.g. 0.95)",       0 },
  { "-s",           eslARG_INT,   "0",   NULL, "n>=0",  NULL, NULL, NULL,    "RNG seed (0 = pick random)",          0 },
  { "--vitlocal",   eslARG_NONE,  FALSE, NULL, NULL,    NULL, NULL, NULL,    "use LOCAL p7 profile for Viterbi (default GLOCAL)", 0 },
  { "--dump",       eslARG_NONE,  FALSE, NULL, NULL,    NULL, NULL, NULL,    "dump per-sample, per-state deficits", 0 },
  { "--pad-out",    eslARG_OUTFILE,NULL, NULL, NULL,    NULL, NULL, NULL,    "write per-state pad vector to <f>",    0 },
  { "--flank",      eslARG_INT,   "500", NULL, "n>=0",  NULL, NULL, NULL,    "embed emit in <n> random residues each side",   0 },
  { "--wscale",     eslARG_REAL,  NULL,  NULL, "x>0",   NULL, NULL, NULL,    "set total embed length to <x>*cm->W (overrides --flank)", 0 },
  {  0,0,0,0,0,0,0,0,0,0 },
};

static char usage[]  = "[-options] <cmfile>";
static char banner[] = "empirical per-state p7 band pad simulation";

static int cmpint(const void *a, const void *b) { int x = *(const int *)a, y = *(const int *)b; return (x<y)?-1:(x>y); }

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go        = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *cmfile    = esl_opt_GetArg(go, 1);
  int             N         = esl_opt_GetInteger(go, "-N");
  double          X         = esl_opt_GetReal(go, "-X");
  uint32_t        seed      = (uint32_t) esl_opt_GetInteger(go, "-s");
  int             do_local  = esl_opt_GetBoolean(go, "--vitlocal");
  int             do_dump   = esl_opt_GetBoolean(go, "--dump");
  char           *pad_outfile = esl_opt_IsOn(go, "--pad-out") ? esl_opt_GetString(go, "--pad-out") : NULL;
  int             flank     = esl_opt_GetInteger(go, "--flank");
  int             use_wscale = esl_opt_IsOn(go, "--wscale");
  double          wscale    = use_wscale ? esl_opt_GetReal(go, "--wscale") : 0.0;

  ESL_RANDOMNESS *r         = esl_randomness_Create(seed);
  ESL_ALPHABET   *abc       = NULL;
  CM_FILE        *cmfp      = NULL;
  CM_t           *cm        = NULL;
  char            errbuf[eslERRBUFSIZE];
  int             status;

  if(cm_file_Open(cmfile, NULL, FALSE, &cmfp, errbuf) != eslOK) cm_Fail("Failed to open CM: %s", errbuf);
  if(cm_file_Read(cmfp, TRUE, &abc, &cm) != eslOK)              cm_Fail("Failed to read CM");
  if(cm_Configure(cm, errbuf, -1) != eslOK)                     cm_Fail("cm_Configure: %s", errbuf);

  if(cm->fp7 == NULL) cm_Fail("CM has no fp7 filter HMM (need cmpress'd model with mlp7)");
  if(cm->cp9map == NULL) cm_Fail("CM has no cp9map (cm_Configure should have built it)");

  P7_HMM     *hmm = cm->fp7;
  P7_BG      *bg  = p7_bg_Create(abc);
  P7_PROFILE *gm  = p7_profile_Create(hmm->M, abc);
  /* Configure profile for length L=400 (will reconfig per emitted seq) */
  if(do_local) p7_ProfileConfig(hmm, bg, gm, 400, p7_GLOCAL); /* will overwrite */
  P7_GMX     *gx  = p7_gmx_Create(hmm->M, 400);
  P7_TRACE   *p7tr = p7_trace_Create();

  int M = hmm->M;

  /* Per-HMM-node deficit storage. Use dynamic arrays. */
  int **deficits   = malloc(sizeof(int *) * (M + 1));
  int  *def_n      = calloc(M + 1, sizeof(int));
  int  *def_alloc  = calloc(M + 1, sizeof(int));
  for(int k = 0; k <= M; k++) { deficits[k] = NULL; def_alloc[k] = 0; }

  int total_emit_obs = 0;

  for(int s = 0; s < N; s++) {
    Parsetree_t *cm_tr = NULL;
    ESL_SQ      *esq   = NULL;
    int          L;
    char         name[32];
    snprintf(name, sizeof(name), "sim%d", s);
    if((status = EmitParsetree(cm, errbuf, r, name, TRUE, &cm_tr, &esq, &L)) != eslOK)
      cm_Fail("EmitParsetree: %s", errbuf);

    /* Embed emitted seq in random flanking residues, mimicking real cmsearch
     * use where hits live inside long sequences. If --wscale is set, total
     * embed length = round(wscale * cm->W); flanks fill to reach that length.
     * Otherwise use fixed --flank residues on each side. */
    int sample_flank = flank;
    if (use_wscale) {
      int L_target = (int)(wscale * cm->W + 0.5);
      sample_flank = (L_target - L) / 2;
      if (sample_flank < 0) sample_flank = 0;
    }
    int L_emb = 2 * sample_flank + L;
    ESL_DSQ *emb = malloc(sizeof(ESL_DSQ) * (L_emb + 2));
    emb[0] = emb[L_emb + 1] = eslDSQ_SENTINEL;
    for(int p = 1; p <= sample_flank;        p++) emb[p] = esl_rnd_FChoose(r, bg->f, abc->K);
    for(int p = 1; p <= L;                   p++) emb[sample_flank + p] = esq->dsq[p];
    for(int p = sample_flank + L + 1; p <= L_emb; p++) emb[p] = esl_rnd_FChoose(r, bg->f, abc->K);
    int true_offset = sample_flank; /* emit position p maps to embedded position p+sample_flank */

    /* Reconfigure p7 profile to embedded seq length and run Viterbi banding (pad=0) */
    p7_ProfileConfig(hmm, bg, gm, L_emb, do_local ? p7_LOCAL : p7_GLOCAL);
    p7_gmx_GrowTo(gx, M, L_emb);
    int *i2k = NULL, *kmin = NULL, *kmax = NULL, ncells = 0;
    if(p7_Seq2BandsVit(errbuf, gm, gx, bg, p7tr, emb, L_emb, /*pad=*/0, /*nodepad=*/NULL, /*hopback=*/0,
                       &i2k, &kmin, &kmax, &ncells) != eslOK) {
      free(i2k); free(kmin); free(kmax); free(emb);
      FreeParsetree(cm_tr); esl_sq_Destroy(esq);
      continue;
    }
    /* Sanity-check the Viterbi trace: it must contain enough M states to
     * be a real alignment of the embedded emit. Count distinct k values
     * in i2k; if too few, skip this sample (Viterbi found only spurious
     * matches in the flanks, e.g., emitted seq scored lower than noise). */
    {
      int distinct_k = 0;
      static int seen[10000];
      for(int kk = 0; kk <= M && kk < 10000; kk++) seen[kk] = 0;
      for(int p = 1; p <= L_emb; p++) {
        if(i2k[p] >= 1 && i2k[p] <= M && i2k[p] < 10000 && !seen[i2k[p]]) {
          seen[i2k[p]] = 1;
          distinct_k++;
        }
      }
      if(do_dump) fprintf(stderr, "# sample=%d distinct_k=%d M=%d filt_thresh=%d\n", s, distinct_k, M, M/2);
      if(distinct_k < M / 2) {
        /* Vit found <50% of HMM states; alignment is degenerate, skip. */
        free(i2k); free(kmin); free(kmax); free(emb);
        FreeParsetree(cm_tr); esl_sq_Destroy(esq);
        continue;
      }
      /* Check that Vit found the embedded emit, not a random flanking match.
       * Count pins in the emit region [sample_flank+1..sample_flank+L]. */
      int emit_pins = 0;
      for(int p = sample_flank + 1; p <= sample_flank + L; p++) {
        if(i2k[p] >= 1 && i2k[p] <= M) emit_pins++;
      }
      if(emit_pins < L / 4) {
        if(do_dump) fprintf(stderr, "# sample=%d skipped: only %d emit pins (need %d)\n", s, emit_pins, L/4);
        free(i2k); free(kmin); free(kmax); free(emb);
        FreeParsetree(cm_tr); esl_sq_Destroy(esq);
        continue;
      }
    }
    if(ncells == 0) {
      free(i2k); free(kmin); free(kmax); free(emb);
      FreeParsetree(cm_tr); esl_sq_Destroy(esq);
      continue;
    }

    /* Walk parsetree, for each emitting CM state record deficit at the
     * true HMM node k vs the band at residue position i (or j). */
    for(int t = 0; t < cm_tr->n; t++) {
      int v = cm_tr->state[t];
      int ipos = cm_tr->emitl[t];
      int jpos = cm_tr->emitr[t];
      int hn1  = cm->cp9map->cs2hn[v][0];
      int hn2  = cm->cp9map->cs2hn[v][1];
      int hs1  = cm->cp9map->cs2hs[v][0];
      int hs2  = cm->cp9map->cs2hs[v][1];

      /* Helper macro: record deficit for (true_k, residue_pos) into bin true_k */
      #define RECORD(true_k, pos_local) do {                                        \
        int pos = (pos_local) + true_offset;                                        \
        if((true_k) >= 1 && (true_k) <= M && pos >= 1 && pos <= L_emb) {            \
          int def;                                                                  \
          if(kmin[pos] == -1 || kmax[pos] == -1) {                                  \
            def = M; /* no band at this pos: maximum deficit */                     \
          } else if((true_k) >= kmin[pos] && (true_k) <= kmax[pos]) {               \
            def = 0;                                                                \
          } else if((true_k) < kmin[pos]) {                                         \
            def = kmin[pos] - (true_k);                                             \
          } else {                                                                  \
            def = (true_k) - kmax[pos];                                             \
          }                                                                         \
          int bin = (true_k);                                                       \
          if(def_n[bin] >= def_alloc[bin]) {                                        \
            int newsz = def_alloc[bin] ? def_alloc[bin] * 2 : 16;                   \
            deficits[bin] = realloc(deficits[bin], sizeof(int) * newsz);            \
            def_alloc[bin] = newsz;                                                 \
          }                                                                         \
          deficits[bin][def_n[bin]++] = def;                                        \
          total_emit_obs++;                                                         \
          if(do_dump) printf("# sample=%d v=%d pos=%d true_k=%d kmin=%d kmax=%d def=%d\n", \
                             s, v, pos, (true_k), kmin[pos], kmax[pos], def);       \
        }                                                                           \
      } while(0)

      /* Determine residue position for hn1/hn2 based on state type.
       * For most match states the consensus position is at i (left);
       * MATR_MR consensus is at j (right). For MATP_MP, hn1 is left
       * (residue i), hn2 is right (residue j). For inserts: IL emits i,
       * IR emits j. */
      if(hn1 >= 0) {
        int pos1;
        if(hs1 == 0) { /* match */
          pos1 = (cm->stid[v] == MATR_MR) ? jpos : ipos;
          RECORD(hn1, pos1);
        } else if(hs1 == 1) { /* insert */
          pos1 = (cm->sttype[v] == IR_st) ? jpos : ipos;
          RECORD(hn1, pos1);
        }
        /* deletes: no residue, skip */
      }
      if(hn2 >= 0) {
        if(hs2 == 0) { RECORD(hn2, jpos); }
        else if(hs2 == 1) { RECORD(hn2, jpos); }
      }

      #undef RECORD
    }

    free(i2k); free(kmin); free(kmax); free(emb);
    FreeParsetree(cm_tr); esl_sq_Destroy(esq);
  }

  /* Aggregate: for each HMM node, sort deficits and pick the X-percentile */
  printf("# CM: %s   M=%d   N=%d  X=%.3f  total_obs=%d\n", cm->name, M, N, X, total_emit_obs);
  printf("# k\tpad\tn_obs\tmean_def\tmax_def\n");
  int *pad_per_k = calloc(M + 1, sizeof(int));
  int max_pad = 0;
  for(int k = 1; k <= M; k++) {
    int n = def_n[k];
    if(n == 0) {
      pad_per_k[k] = 0;
      printf("%d\t0\t0\t0.0\t0\n", k);
      continue;
    }
    qsort(deficits[k], n, sizeof(int), cmpint);
    int idx = (int)(X * (n - 1) + 0.5);
    if(idx >= n) idx = n - 1;
    int pad = deficits[k][idx];
    int max = deficits[k][n - 1];
    double mean = 0; for(int i = 0; i < n; i++) mean += deficits[k][i]; mean /= n;
    printf("%d\t%d\t%d\t%.2f\t%d\n", k, pad, n, mean, max);
    pad_per_k[k] = pad;
    if(pad > max_pad) max_pad = pad;
  }
  printf("# max_pad_over_all_nodes=%d\n", max_pad);

  /* Histogram: count states with each pad value */
  int *pad_hist = calloc(max_pad + 2, sizeof(int));
  int  unset_n = 0;
  for(int k = 1; k <= M; k++) {
    if(def_n[k] == 0) { unset_n++; continue; }
    int idx = (int)(X * (def_n[k] - 1) + 0.5);
    if(idx >= def_n[k]) idx = def_n[k] - 1;
    pad_hist[deficits[k][idx]]++;
  }
  printf("# Pad histogram (pad value : count : cumulative %% of M=%d)\n", M);
  int cum = unset_n;
  if(unset_n > 0) printf("# unset\t%d\t%5.1f\n", unset_n, 100.0 * cum / M);
  for(int p = 0; p <= max_pad; p++) {
    if(pad_hist[p] == 0) continue;
    cum += pad_hist[p];
    printf("# pad=%d\t%d\t%5.1f\n", p, pad_hist[p], 100.0 * cum / M);
  }
  free(pad_hist);

  /* Write per-state pad vector to file (--pad-out).
   * Format: header line "# p7bandsim pad-vector <CM_name> M=<M> X=<X> N=<N>"
   * followed by M lines "<k> <pad>". */
  if(pad_outfile != NULL) {
    FILE *pf = fopen(pad_outfile, "w");
    if(pf == NULL) cm_Fail("Failed to open --pad-out file for writing: %s", pad_outfile);
    fprintf(pf, "# p7bandsim pad-vector  CM=%s  M=%d  X=%.3f  N=%d  total_obs=%d\n",
            cm->name, M, X, N, total_emit_obs);
    for(int k = 0; k <= M; k++) fprintf(pf, "%d\t%d\n", k, pad_per_k[k]);
    fclose(pf);
  }
  free(pad_per_k);

  /* Cleanup */
  for(int k = 0; k <= M; k++) if(deficits[k]) free(deficits[k]);
  free(deficits); free(def_n); free(def_alloc);
  p7_trace_Destroy(p7tr);
  p7_gmx_Destroy(gx);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  FreeCM(cm);
  cm_file_Close(cmfp);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}

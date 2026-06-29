/* test_pins2trace.c -- brief 137 Phase 3a regression test.
 *
 * For each sequence, compare:
 *   reference : p7_GViterbi -> p7_GTrace on the full (unbanded) profile
 *   test      : p7_Seq2BandsIBV_dnc(delta=0) -> p7_IBVPins2Trace
 * on the core M/I residue->model-position assignment (alignment-equivalence).
 *
 * The IBV DP is glocal-global, so the meaningful comparison is p7_UNIGLOCAL
 * (-g).  Use -l to force p7_UNILOCAL instead.
 *
 *   gcc ... -o test_pins2trace test_pins2trace.c -linfernal -lhmmer -leasel -lm
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
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "infernal.h"

static ESL_OPTIONS options[] = {
  { "-h", eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "show help",                       0 },
  { "-g", eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, "-l", "glocal: p7_UNIGLOCAL [default]",  0 },
  { "-l", eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, "-g", "local:  p7_UNILOCAL",             0 },
  { "-v", eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "verbose: dump per-seq mismatch",  0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <cmfile> <seqfile>";
static char banner[] = "brief137 pins-to-trace regression test";

/* Build per-residue (state,k) arrays for emitted residues 1..L from a trace.
 * state[i] in {p7T_M, p7T_I, 0}; 0 means residue i emitted by N/C (flanking). */
static void
trace_to_resmap(const P7_TRACE *tr, int L, int *state, int *kof)
{
  int z;
  for (int i = 0; i <= L; i++) { state[i] = 0; kof[i] = 0; }
  for (z = 0; z < tr->N; z++) {
    if ((tr->st[z] == p7T_M || tr->st[z] == p7T_I) && tr->i[z] >= 1 && tr->i[z] <= L) {
      state[tr->i[z]] = tr->st[z];
      kof[tr->i[z]]   = tr->k[z];
    }
  }
}

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char         *cmfile  = esl_opt_GetArg(go, 1);
  char         *seqfile = esl_opt_GetArg(go, 2);
  int           p7mode  = esl_opt_GetBoolean(go, "-l") ? p7_UNILOCAL : p7_UNIGLOCAL;
  int           verbose = esl_opt_GetBoolean(go, "-v");
  char          errbuf[eslERRBUFSIZE];

  ESL_ALPHABET *abc  = NULL;
  CM_FILE      *cmfp = NULL;
  CM_t         *cm   = NULL;
  P7_HMM       *hmm  = NULL;
  P7_BG        *bg   = NULL;
  P7_PROFILE   *gm   = NULL;
  P7_GMX       *gx   = NULL;
  ESL_SQFILE   *sqfp = NULL;
  ESL_SQ       *sq   = NULL;
  int           status;

  /* Read + configure CM (populates cm->fp7). */
  if (cm_file_Open(cmfile, NULL, FALSE, &cmfp, errbuf) != eslOK) p7_Fail("cm_file_Open: %s", errbuf);
  if (cm_file_Read(cmfp, TRUE, &abc, &cm)              != eslOK) p7_Fail("cm_file_Read failed");
  cm_file_Close(cmfp);
  if (cm_Configure(cm, errbuf, -1) != eslOK) p7_Fail("cm_Configure: %s", errbuf);
  if (cm->fp7 == NULL) p7_Fail("CM has no fp7 filter HMM");
  hmm = cm->fp7;

  bg = p7_bg_Create(abc);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, 400, p7mode);
  gx = p7_gmx_Create(hmm->M, 400);

  if (esl_sqfile_Open(seqfile, eslSQFILE_UNKNOWN, NULL, &sqfp) != eslOK) p7_Fail("open seqfile failed");
  sq = esl_sq_CreateDigital(abc);

  printf("# CM=%s  M=%d  mode=%s\n", cmfile, hmm->M, (p7mode==p7_UNIGLOCAL)?"UNIGLOCAL":"UNILOCAL");

  /* Emission-objective audit: IBV match lod = log2(mat/ins); gm MSC = log(mat/bg).
   * Equal iff ins==bg.  Also compare in bits directly (gm MSC is in nats). */
  { double max_insbg=0, max_lod=0;
    for (int k=1;k<hmm->M;k++) for (int x=0;x<abc->K;x++){
      double d = fabs((double)hmm->ins[k][x]-(double)bg->f[x]);
      if (d>max_insbg) max_insbg=d;
    }
    for (int k=1;k<=hmm->M;k++) for (int x=0;x<abc->K;x++){
      if (hmm->mat[k][x]<=0||hmm->ins[k][x]<=0) continue;
      double ibv_bits = log((double)hmm->mat[k][x]/(double)hmm->ins[k][x])/M_LN2;
      double gm_bits  = (double)p7P_MSC(gm,k,x)/M_LN2;   /* MSC stored in nats */
      double d = fabs(ibv_bits-gm_bits);
      if (d>max_lod) max_lod=d;
    }
    printf("# EMISSION AUDIT: max|ins-bg|=%.6g  max|IBVmatch_bits - gmMSC_bits|=%.6g (0 => emissions identical)\n",
           max_insbg, max_lod);
  }
  printf("# %-22s %6s %8s %7s %8s  %5s %5s %5s  %3s %3s %5s  %s\n",
         "seq", "L", "vit_sc", "scgap", "pins_eq", "ncore", "mism", "flank", "rNM", "rIN", "valid", "verdict");

  int tot_seq=0, tot_clean=0, tot_drift=0, tot_fail=0;

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK) {
    int   L = sq->n;
    float vsc;
    int  *i2k=NULL,*kmin=NULL,*kmax=NULL,ncells=0;
    P7_TRACE *rtr=NULL, *ttr=NULL;
    int  *rstate=malloc(sizeof(int)*(L+1)), *rk=malloc(sizeof(int)*(L+1));
    int  *tstate=malloc(sizeof(int)*(L+1)), *tk=malloc(sizeof(int)*(L+1));
    int   ncore=0, match=0, mism=0, flank=0, walker_ok=1, valid_ok=1;

    p7_ReconfigLength(gm, L);
    p7_bg_SetLength(bg, L);

    /* Reference. */
    p7_gmx_GrowTo(gx, hmm->M, L);
    p7_GViterbi(sq->dsq, L, gm, gx, &vsc);
    rtr = p7_trace_Create();
    if (p7_GTrace(sq->dsq, L, gm, gx, rtr) != eslOK) p7_Fail("GTrace failed on %s", sq->name);

    /* Test (IBV pins -> trace), delta=0, no boundary widen. */
    if ((status = p7_Seq2BandsIBV_dnc(cm, errbuf, sq->dsq, L, 0, 0, FALSE, FALSE,
                                      P7IBV_MODE_DELTA, 0, /* brief 140 */
                                      &i2k, &kmin, &kmax, &ncells)) != eslOK)
      p7_Fail("p7_Seq2BandsIBV_dnc failed on %s: %s", sq->name, errbuf);
    /* RAW-PIN validity (mode-independent): is the per-row argmax sequence a
     * realizable monotone Plan7 path?  rawNM = non-monotone steps;
     * rawINF = infeasible insert->match-with-gap steps (I_k then match >k+1,
     * which would require deletes after an insert). */
    int rawNM=0, rawINF=0; {
      int kprev=0, prevdelta=99;   /* prevdelta: i2k[i-1]-i2k[i-2]; 0 => prev was insert */
      for (int i=1;i<=L;i++){
        int d = i2k[i]-kprev;
        if (d<0) rawNM++;
        if (d>1 && prevdelta==0) rawINF++;   /* insert followed by match needing deletes */
        prevdelta = (d<0)?99:d;
        kprev=i2k[i];
      }
    }

    status = p7_IBVPins2Trace(gm, sq->dsq, L, i2k, kmin, kmax, ncells, &ttr);
    if (status != eslOK) { walker_ok = 0; }

    float rsc=0, tsc=0;   /* trace scores under the SAME gm */
    if (walker_ok) {
      p7_trace_Score(rtr, sq->dsq, gm, &rsc);
      p7_trace_Score(ttr, sq->dsq, gm, &tsc);
      if (p7_trace_Validate(ttr, abc, sq->dsq, errbuf) != eslOK) {
        valid_ok = 0; printf("#   %s validate: %s\n", sq->name, errbuf);
        /* find first illegal D-entry or I->E and dump its neighborhood */
        int bad=-1;
        for (int z=1; z<ttr->N; z++) {
          int p=ttr->st[z-1], c=ttr->st[z];
          if (c==p7T_D && p!=p7T_M && p!=p7T_D) { bad=z; break; }
          if (c==p7T_E && p!=p7T_M && p!=p7T_D) { bad=z; break; }
        }
        printf("#     bad@z=%d region:", bad);
        for (int z=(bad>4?bad-4:0); z<ttr->N && z<bad+3; z++)
          printf(" %c%d/%d", "?MDISNBECTJ"[ttr->st[z]], ttr->k[z], ttr->i[z]);
        printf("\n");
      }
      trace_to_resmap(rtr, L, rstate, rk);
      trace_to_resmap(ttr, L, tstate, tk);
      int first_mism=-1, last_mism=-1;
      for (int i = 1; i <= L; i++) {
        int rin = (rstate[i] != 0), tin = (tstate[i] != 0);
        if (!rin && !tin) continue;
        if (rin && tin) {
          ncore++;
          if (rstate[i]==tstate[i] && rk[i]==tk[i]) match++;
          else { mism++; if(first_mism<0)first_mism=i; last_mism=i;
                 if (verbose && mism<=8) printf("#   %s i=%d ref=(%c,%d) test=(%c,%d)\n",
                          sq->name, i, "?MDI"[rstate[i]], rk[i], "?MDI"[tstate[i]], tk[i]); }
        } else flank++;   /* one side put it in core, other in N/C */
      }
      /* #1 mechanism test: is the gm-Viterbi optimal path itself INSIDE the IBV
       * Delta=0 co-optimal band?  If a ref-trace core cell falls outside
       * [kmin,kmax], the IBV DP does not consider gm's optimum co-optimal ->
       * the two DPs optimize different score functions (not a walker artifact). */
      int ref_oob=0, ref_oob_first=-1, ref_oob_k=-1, ref_band_lo=-1, ref_band_hi=-1;
      for (int i=1;i<=L;i++){
        if (rstate[i]!=0) {  /* ref put residue i in core at model pos rk[i] */
          if (rk[i] < kmin[i] || rk[i] > kmax[i]) {
            ref_oob++;
            if (ref_oob_first<0){ ref_oob_first=i; ref_oob_k=rk[i]; ref_band_lo=kmin[i]; ref_band_hi=kmax[i]; }
          }
        }
      }
      if ((rsc-tsc) > 0.01)
        printf("#   %s GAP=%.2f mism=%d flank=%d | refPathOutsideIBVband=%d cells%s\n",
               sq->name, rsc-tsc, mism, flank, ref_oob,
               ref_oob? "" : " (gm-optimum IS in band -> tie/quantization)");
      if (ref_oob)
        printf("#       first ref-OOB @i=%d: gm wants M%d but IBV band=[%d,%d]\n",
               ref_oob_first, ref_oob_k, ref_band_lo, ref_band_hi);
    }

    const char *verdict;
    if (!walker_ok)            { verdict="WALKER_ERR"; tot_fail++; }
    else if (!valid_ok)       { verdict="INVALID_TR"; tot_fail++; }
    else if (mism==0 && flank==0) { verdict="clean";   tot_clean++; }
    else                      { verdict="drift";       tot_drift++; }
    tot_seq++;

    double eq = (ncore+flank>0) ? (100.0*match/(ncore+flank)) : 0.0;
    /* score gap: GViterbi-trace score minus IBV-pins-trace score, both under gm.
     * 0 => pins are a co-optimal Viterbi alignment (cell diffs are ties). */
    float scgap = walker_ok ? (rsc - tsc) : 0.0;
    printf("  %-22s %6d %8.2f %7.2f %7.1f%%  %5d %5d %5d  %3d %3d %5s  %s\n",
           sq->name, L, vsc, scgap, eq, ncore, mism, flank, rawNM, rawINF,
           walker_ok ? (valid_ok?"yes":"NO") : "-", verdict);
    if (!walker_ok) printf("#   walker returned status %d on %s\n", status, sq->name);

    free(i2k); free(kmin); free(kmax);
    free(rstate); free(rk); free(tstate); free(tk);
    p7_trace_Destroy(rtr); if (ttr) p7_trace_Destroy(ttr);
    esl_sq_Reuse(sq);
  }

  printf("# ---- %d seqs: %d clean, %d drift, %d fail ----\n",
         tot_seq, tot_clean, tot_drift, tot_fail);

  esl_sqfile_Close(sqfp);
  esl_sq_Destroy(sq);
  p7_gmx_Destroy(gx);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  FreeCM(cm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}

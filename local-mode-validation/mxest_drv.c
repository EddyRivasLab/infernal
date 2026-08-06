/* mxest_drv.c : brief 26_0430-225 self-check driver.
 *
 * For each sequence: derives CP9 bands, prints PREDICTED peak Mb from the
 * new pre-alignment estimators (cm_CheckptAlignSizeNeededHB /
 * cm_CheckptTrAlignSizeNeededHB / cm_DnCAlignSizeNeededHB /
 * cm_TrDnCAlignSizeNeededHB), then calls the REAL engines directly
 * (test-only, bypassing their _Qualifies() gates, mirroring the existing
 * hbdnc_drv.c / rl3_ckpttralign_drv.c pattern) with INFERNAL_CKPT_VERBOSE
 * and DNC_MEM_VERBOSE set so the MEASURED ground-truth peak Mb prints to
 * stderr right after each PREDICTED line, for by-hand comparison.
 *
 * kpin/bkind are passed as NULL to the bps>0 Post/OptAcc passes: confirmed
 * (cm_dpalign.c reading, brief 26_0430-225) that kpin only narrows a B_st's
 * k-search RANGE, never affects any deck alloc/free call -- so NULL (full,
 * unpinned combine) gives byte-identical peak_bytes to a real pinned run,
 * while sidestepping the need to thread kpin from a pass-1 parse tree.
 *
 * Usage: mxest_drv [-g] [--tau x] [--notrunc] <cmfile> <seqfile>
 */
#include <esl_config.h>
#include <p7_config.h>
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
  { "-h",        eslARG_NONE, FALSE,  NULL, NULL, NULL, NULL, NULL, "help", 0 },
  { "-g",        eslARG_NONE, FALSE,  NULL, NULL, NULL, NULL, NULL, "global config (default local)", 0 },
  { "--tau",     eslARG_REAL, "1e-7", NULL, NULL, NULL, NULL, NULL, "HMM band tail loss tau", 0 },
  { "--mxsize",  eslARG_REAL, "40000",NULL, NULL, NULL, NULL, NULL, "size limit (Mb) passed to the real engines", 0 },
  { "--notrunc", eslARG_NONE, FALSE,  NULL, NULL, NULL, NULL, NULL, "skip the truncated (rung-4) checks", 0 },
  { "--nodnc",   eslARG_NONE, FALSE,  NULL, NULL, NULL, NULL, NULL, "skip D&C checks (genome-scale: D&C ground truth can be slow)", 0 },
  { 0,0,0,0,0,0,0,0,0,0 },
};

int main(int argc, char **argv)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  char errbuf[eslERRBUFSIZE];
  int status;
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK || esl_opt_VerifyConfig(go) != eslOK) {
    printf("bad opts: %s\n", go->errbuf); exit(1);
  }
  if (esl_opt_GetBoolean(go, "-h")) { printf("usage: mxest_drv [-g] [--tau x] [--notrunc] [--nodnc] <cmfile> <seqfile>\n"); exit(0); }
  int   do_global = esl_opt_GetBoolean(go, "-g");
  float tau       = esl_opt_GetReal(go, "--tau");
  int   do_trunc  = ! esl_opt_GetBoolean(go, "--notrunc");
  int   do_dnc    = ! esl_opt_GetBoolean(go, "--nodnc");
  float size_limit= esl_opt_GetReal(go, "--mxsize");
  if (esl_opt_ArgNumber(go) != 2) { printf("usage: mxest_drv [-g] [--tau x] [--notrunc] [--nodnc] <cmfile> <seqfile>\n"); exit(1); }
  char *cmfile  = esl_opt_GetArg(go, 1);
  char *seqfile = esl_opt_GetArg(go, 2);

  setenv("INFERNAL_CKPT_VERBOSE", "1", 1);
  setenv("DNC_MEM_VERBOSE", "1", 1);

  CM_FILE *cmfp = NULL;
  ESL_ALPHABET *abc = NULL;
  CM_t *cm = NULL;
  if (cm_file_Open(cmfile, NULL, TRUE, &cmfp, errbuf) != eslOK) cm_Fail("open cm: %s", errbuf);
  if (cm_file_Read(cmfp, TRUE, &abc, &cm) != eslOK) cm_Fail("read cm: %s", cmfp->errbuf);
  cm_file_Close(cmfp);

  cm->align_opts |= CM_ALIGN_HBANDED;
  if (! do_global) {
    cm->config_opts |= CM_CONFIG_LOCAL;
    cm->config_opts |= CM_CONFIG_HMMLOCAL;
    cm->config_opts |= CM_CONFIG_HMMEL;
  }
  if (do_trunc) cm->config_opts |= CM_CONFIG_TRUNC; /* allocates cm->trhb_emx etc (cm_modelconfig.c:218-223) */
  cm->tau = tau;
  if (cm_Configure(cm, errbuf, -1) != eslOK) cm_Fail("configure: %s", errbuf);

  init_ilogsum();
  FLogsumInit();

  int has_bif = (CMCountStatetype(cm, B_st) > 0) ? TRUE : FALSE;
  printf("# CM: %s  M=%d clen=%d  mode=%s  tau=%g  has_bif=%d\n",
         cm->name, cm->M, cm->clen, do_global?"global":"local", tau, has_bif);

  ESL_SQFILE *sqfp = NULL;
  if (esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_UNKNOWN, NULL, &sqfp) != eslOK) cm_Fail("open seqfile %s", seqfile);
  ESL_SQ *sq = esl_sq_CreateDigital(abc);

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK) {
    int L = sq->n;
    printf("# --- seq %s/%d ---\n", sq->name, L);

    if ((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx,
                                sq->dsq, 1, L, cm->cp9b, FALSE, PLI_PASS_STD_ANY, 0)) != eslOK)
      cm_Fail("cp9_Seq2Bands: %s", errbuf);

    /* ---- non-truncated --ckpt: predicted vs measured ---- */
    float ckptdpmb=0., emxmb=0., cp9mxmb=0., totmb=0.;
    /* this driver derives bands via unbanded cp9_Seq2Bands() above -- pass
     * NULL/NULL (brief 26_0430-226) so cp9mxmb reflects that unbanded cost. */
    if ((status = cm_CheckptAlignSizeNeededHB(cm, errbuf, L, NULL, NULL, &ckptdpmb, &emxmb, &cp9mxmb, &totmb)) != eslOK)
      cm_Fail("cm_CheckptAlignSizeNeededHB: %s", errbuf);
    printf("PREDICT ckpt notrunc: ckptdpMb=%.4f emxMb=%.4f cp9mxMb=%.4f totMb=%.4f\n", ckptdpmb, emxmb, cp9mxmb, totmb);
    fflush(stdout);
    if (! has_bif) {
      char *ppstr=NULL; Parsetree_t *tr=NULL; float avgpp=0., sc=0.;
      status = cm_CheckptAlignHB(cm, errbuf, sq->dsq, L, size_limit, cm->hb_emx, &ppstr, &tr, &avgpp, &sc);
      if (status != eslOK) fprintf(stderr, "# MEASURE ckpt notrunc bps0 FAILED: %s\n", errbuf);
      if (ppstr) free(ppstr); if (tr) FreeParsetree(tr);
    } else {
      Parsetree_t *tr_cyk=NULL; float cyk_sc=0., post_Z=0.;
      char *ppstr=NULL; Parsetree_t *tr=NULL; float avgpp=0., pp=0.;
      status = cm_CheckptCYKAlignHB(cm, errbuf, sq->dsq, L, size_limit, &tr_cyk, &cyk_sc);
      if (status != eslOK) fprintf(stderr, "# MEASURE ckpt notrunc CYK-pass FAILED: %s\n", errbuf);
      if (tr_cyk) FreeParsetree(tr_cyk);
      status = cm_CheckptPostAlignHB(cm, errbuf, sq->dsq, L, size_limit, cm->hb_emx, NULL, &post_Z);
      if (status != eslOK) fprintf(stderr, "# MEASURE ckpt notrunc Post-pass FAILED: %s\n", errbuf);
      /* brief 26_0430-225 driver note: kpin=NULL (full/unpinned combine) has
       * been observed to produce a garbage Z in GLOBAL mode for this bps>0
       * pass (root cause not chased -- out of scope for this brief, an
       * existing engine property of the unpinned combine, not something the
       * new estimator functions depend on: kpin never affects any alloc/
       * free call, only VALUES). Skip the OA-pass call rather than crash on
       * downstream garbage. */
      if (status == eslOK && isfinite(post_Z) && fabs(post_Z) < 1e6) {
        status = cm_CheckptOptAccAlignHB(cm, errbuf, sq->dsq, L, size_limit, cm->hb_emx, NULL, &ppstr, &tr, &avgpp, &pp);
        if (status != eslOK) fprintf(stderr, "# MEASURE ckpt notrunc OA-pass FAILED: %s\n", errbuf);
      } else {
        fprintf(stderr, "# MEASURE ckpt notrunc OA-pass SKIPPED (Post-pass Z=%.4f looked unsafe with kpin=NULL)\n", post_Z);
      }
      if (ppstr) free(ppstr); if (tr) FreeParsetree(tr);
    }

    /* ---- non-truncated D&C: predicted vs measured ---- */
    if (do_dnc) {
      float vjdmb=0., shmb=0., dnctotmb=0.;
      if ((status = cm_DnCAlignSizeNeededHB(cm, errbuf, L, &vjdmb, &shmb, &dnctotmb)) != eslOK)
        cm_Fail("cm_DnCAlignSizeNeededHB: %s", errbuf);
      printf("PREDICT dnc notrunc: vjdMb=%.4f shMb=%.4f totMb=%.4f\n", vjdmb, shmb, dnctotmb);
      fflush(stdout);
      Parsetree_t *tr = NULL;
      float sc = CYKDivideAndConquerHB(cm, sq->dsq, L, 0, 1, L, &tr, cm->cp9b);
      (void) sc;
      if (tr) FreeParsetree(tr);
    }

    /* ---- truncated (rung-4), mode J only for this spot-check ---- */
    if (do_trunc) {
      int fill_L, fill_R, fill_T;
      cm_TrFillFromMode(TRMODE_J, &fill_L, &fill_R, &fill_T);
      float tckptdpmb=0., temxmb=0., tcp9mxmb=0., ttotmb=0.;
      if ((status = cm_CheckptTrAlignSizeNeededHB(cm, errbuf, L, TRMODE_J, NULL, NULL, &tckptdpmb, &temxmb, &tcp9mxmb, &ttotmb)) != eslOK)
        cm_Fail("cm_CheckptTrAlignSizeNeededHB: %s", errbuf);
      printf("PREDICT ckpt trunc(J): ckptdpMb=%.4f emxMb=%.4f cp9mxMb=%.4f totMb=%.4f\n", tckptdpmb, temxmb, tcp9mxmb, ttotmb);
      fflush(stdout);
      if (! has_bif) {
        char *ppstr=NULL; Parsetree_t *tr=NULL; char rmode=0; float avgpp=0., sc=0.;
        status = cm_CheckptTrAlignHB(cm, errbuf, sq->dsq, L, size_limit, TRMODE_J, PLI_PASS_5P_AND_3P_FORCE,
                                      cm->trhb_emx, &ppstr, &tr, &rmode, &avgpp, &sc);
        if (status != eslOK) fprintf(stderr, "# MEASURE ckpt trunc bps0 FAILED: %s\n", errbuf);
        if (ppstr) free(ppstr); if (tr) FreeParsetree(tr);
      } else {
        Parsetree_t *tr_cyk=NULL; char rmode=0; float cyk_sc=0., post_sc=0.;
        char *ppstr=NULL; Parsetree_t *tr=NULL; float avgpp=0., pp=0.;
        status = cm_CheckptTrCYKAlignHB(cm, errbuf, sq->dsq, L, size_limit, PLI_PASS_5P_AND_3P_FORCE, &tr_cyk, &rmode, &cyk_sc);
        if (status != eslOK) fprintf(stderr, "# MEASURE ckpt trunc CYK-pass FAILED: %s\n", errbuf);
        if (tr_cyk) FreeParsetree(tr_cyk);
        status = cm_CheckptTrPostAlignHB(cm, errbuf, sq->dsq, L, size_limit, TRMODE_J, PLI_PASS_5P_AND_3P_FORCE,
                                          cm->trhb_emx, NULL, NULL, NULL, NULL, NULL, &post_sc, &rmode);
        if (status != eslOK) fprintf(stderr, "# MEASURE ckpt trunc Post-pass FAILED: %s\n", errbuf);
        if (status == eslOK && isfinite(post_sc) && fabs(post_sc) < 1e6) {
          status = cm_CheckptTrOptAccAlignHB(cm, errbuf, sq->dsq, L, size_limit, TRMODE_J, PLI_PASS_5P_AND_3P_FORCE,
                                              cm->trhb_emx, NULL, NULL, NULL, NULL, NULL, &ppstr, &tr, &rmode, &avgpp, &pp);
          if (status != eslOK) fprintf(stderr, "# MEASURE ckpt trunc OA-pass FAILED: %s\n", errbuf);
        } else {
          fprintf(stderr, "# MEASURE ckpt trunc OA-pass SKIPPED (Post-pass sc=%.4f looked unsafe with kpin=NULL)\n", post_sc);
        }
        if (ppstr) free(ppstr); if (tr) FreeParsetree(tr);
      }

      if (do_dnc) {
        float tvjdmb=0., tshmb=0., tdnctotmb=0.;
        if ((status = cm_TrDnCAlignSizeNeededHB(cm, errbuf, L, TRMODE_J, &tvjdmb, &tshmb, &tdnctotmb)) != eslOK)
          cm_Fail("cm_TrDnCAlignSizeNeededHB: %s", errbuf);
        printf("PREDICT dnc trunc(J): vjdMb=%.4f shMb=%.4f totMb=%.4f\n", tvjdmb, tshmb, tdnctotmb);
        fflush(stdout);
        Parsetree_t *tr = NULL; char rmode = 0;
        float sc = TrCYKDivideAndConquerHB(cm, sq->dsq, L, 0, 1, L, PLI_PASS_5P_AND_3P_FORCE, TRMODE_J, &rmode, &tr, cm->cp9b);
        (void) sc;
        if (tr) FreeParsetree(tr);
      }
    }

    esl_sq_Reuse(sq);
  }

  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  FreeCM(cm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}

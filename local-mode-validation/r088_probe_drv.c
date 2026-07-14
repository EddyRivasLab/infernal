/* r088_probe_drv.c : brief 26_0610-088 investigation.
 *
 * For each seq, FORCE L mode and dump BOTH parses:
 *   (a) the D&C TrCYKDivideAndConquerHB(preset=L) parse (the buggy 4-node one), and
 *   (b) the trusted monolithic oracle cm_TrAlignHB(preset=L, do_optacc=FALSE) parse.
 * Prints scores and full ParsetreeDump for eyeball comparison of the structural collapse.
 *
 * Usage: r088_probe_drv [-g] [--tau x] [--mode J|L|R] <cmfile> <seqfile>
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
  { "-h",    eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "help", 0 },
  { "-g",    eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "global config (default local)", 0 },
  { "--tau", eslARG_REAL,  "1e-7",NULL, NULL, NULL, NULL, NULL, "HMM band tail loss tau", 0 },
  { "--mode",eslARG_STRING,"L",   NULL, NULL, NULL, NULL, NULL, "force mode J/L/R", 0 },
  { "--mxsize", eslARG_REAL,"40000",NULL,NULL,NULL, NULL, NULL, "size limit (Mb)", 0 },
  { 0,0,0,0,0,0,0,0,0,0 },
};

int main(int argc, char **argv)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  char errbuf[eslERRBUFSIZE];
  int status;
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK || esl_opt_VerifyConfig(go) != eslOK)
    { printf("bad opts: %s\n", go->errbuf); exit(1); }
  if (esl_opt_GetBoolean(go, "-h") || esl_opt_ArgNumber(go) != 2)
    { printf("usage: r088_probe_drv [-g] [--tau x] [--mode J|L|R] <cmfile> <seqfile>\n"); exit(0); }
  int    do_global = esl_opt_GetBoolean(go, "-g");
  float  tau       = esl_opt_GetReal(go, "--tau");
  float  size_limit= esl_opt_GetReal(go, "--mxsize");
  char  *modestr   = esl_opt_GetString(go, "--mode");
  char  *cmfile    = esl_opt_GetArg(go, 1);
  char  *seqfile   = esl_opt_GetArg(go, 2);
  int    pass_idx  = PLI_PASS_5P_AND_3P_ANY;
  char   fmode = (modestr[0]=='J')?TRMODE_J:(modestr[0]=='R')?TRMODE_R:TRMODE_L;

  CM_FILE *cmfp = NULL; ESL_ALPHABET *abc = NULL; CM_t *cm = NULL;
  if (cm_file_Open(cmfile, NULL, TRUE, &cmfp, errbuf) != eslOK) cm_Fail("open cm: %s", errbuf);
  if (cm_file_Read(cmfp, TRUE, &abc, &cm) != eslOK) cm_Fail("read cm: %s", cmfp->errbuf);
  cm_file_Close(cmfp);
  cm->config_opts |= CM_CONFIG_TRUNC;
  cm->align_opts  |= CM_ALIGN_HBANDED;
  if (! do_global) { cm->config_opts |= CM_CONFIG_LOCAL; cm->config_opts |= CM_CONFIG_HMMLOCAL; cm->config_opts |= CM_CONFIG_HMMEL; }
  cm->tau = tau;
  if (cm_Configure(cm, errbuf, -1) != eslOK) cm_Fail("configure: %s", errbuf);
  init_ilogsum(); FLogsumInit();

  printf("# CM=%s M=%d mode=%s FORCE=%c tau=%g LB=%d LE=%d\n", cm->name, cm->M,
         do_global?"global":"local", MarginalMode(fmode), tau,
         (cm->flags&CMH_LOCAL_BEGIN)?1:0, (cm->flags&CMH_LOCAL_END)?1:0);

  ESL_SQFILE *sqfp = NULL;
  ESL_SQ *sq = esl_sq_CreateDigital(abc);
  if (esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_UNKNOWN, NULL, &sqfp) != eslOK) cm_Fail("open seqfile %s", seqfile);
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK) {
    int L = sq->n;
    ESL_DSQ *dsq = sq->dsq;
    if ((status = cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, dsq, 1, L, cm->cp9b, FALSE, pass_idx, 0)) != eslOK)
      cm_Fail("cp9_Seq2Bands: %s", errbuf);

    printf("\n========================================================\n");
    printf("SEQ=%s L=%d  Jvalid[0]=%d Lvalid[0]=%d Rvalid[0]=%d Tvalid[0]=%d\n",
           sq->name, L, cm->cp9b->Jvalid[0], cm->cp9b->Lvalid[0], cm->cp9b->Rvalid[0], cm->cp9b->Tvalid[0]);

    /* (a) D&C forced-mode */
    Parsetree_t *tr_dnc = NULL; char m_dnc = TRMODE_UNKNOWN;
    float sc_dnc = TrCYKDivideAndConquerHB(cm, dsq, L, 0, 1, L, pass_idx, fmode, &m_dnc, &tr_dnc, cm->cp9b);
    float ptsc_dnc = 0.; char pterr[eslERRBUFSIZE];
    ParsetreeScore(cm, NULL, NULL, tr_dnc, dsq, FALSE, &ptsc_dnc, NULL, NULL, NULL, NULL);

    /* (b) oracle forced-mode */
    CM_TR_HB_MX *mx = cm_tr_hb_mx_Create(cm);
    CM_TR_HB_SHADOW_MX *shmx = cm_tr_hb_shadow_mx_Create(cm);
    Parsetree_t *tr_or = NULL; char m_or = TRMODE_UNKNOWN; float sc_or = IMPOSSIBLE;
    status = cm_TrAlignHB(cm, errbuf, dsq, L, size_limit, fmode, pass_idx, FALSE, FALSE,
                          mx, shmx, NULL, NULL, NULL, NULL, &tr_or, &m_or, NULL, &sc_or);
    if (status != eslOK) cm_Fail("oracle fail %s: %s", sq->name, errbuf);
    float ptsc_or = 0.;
    ParsetreeScore(cm, NULL, NULL, tr_or, dsq, FALSE, &ptsc_or, NULL, NULL, NULL, NULL);

    printf("--- D&C forced-%c: sc=%.6f ParsetreeScore=%.6f nodes=%d mode_ret=%s\n",
           MarginalMode(fmode), sc_dnc, ptsc_dnc, tr_dnc?tr_dnc->n:-1, MarginalMode(m_dnc));
    ParsetreeDump(stdout, tr_dnc, cm, dsq);
    printf("--- ORACLE forced-%c: sc=%.6f ParsetreeScore=%.6f nodes=%d mode_ret=%s\n",
           MarginalMode(fmode), sc_or, ptsc_or, tr_or?tr_or->n:-1, MarginalMode(m_or));
    ParsetreeDump(stdout, tr_or, cm, dsq);

    if (tr_dnc) FreeParsetree(tr_dnc);
    if (tr_or)  FreeParsetree(tr_or);
    cm_tr_hb_mx_Destroy(mx); cm_tr_hb_shadow_mx_Destroy(shmx);
    esl_sq_Reuse(sq);
  }
  esl_sq_Destroy(sq); esl_sqfile_Close(sqfp);
  FreeCM(cm); esl_alphabet_Destroy(abc); esl_getopts_Destroy(go);
  return 0;
}

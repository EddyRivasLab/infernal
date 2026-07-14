/* r093_dump_drv.c : brief 26_0610-093 -- dump D&C vs oracle parsetrees for a forced mode.
 * Usage: r093_dump_drv [-g] [--tau x] [--mode J|L|R|T] <cmfile> <seqfile>
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
  { "--mode",eslARG_STRING, "L",  NULL, NULL, NULL, NULL, NULL, "force mode J/L/R/T", 0 },
  { "--mxsize", eslARG_REAL,"40000",NULL,NULL,NULL, NULL, NULL, "size limit (Mb)", 0 },
  { 0,0,0,0,0,0,0,0,0,0 },
};

int main(int argc, char **argv)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  char errbuf[eslERRBUFSIZE];
  int status;
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) { printf("bad opts: %s\n", go->errbuf); exit(1); }
  if (esl_opt_GetBoolean(go, "-h") || esl_opt_ArgNumber(go) != 2) { printf("usage: r093_dump_drv [-g] [--tau x] [--mode L] <cmfile> <seqfile>\n"); exit(0); }
  int    do_global = esl_opt_GetBoolean(go, "-g");
  float  tau       = esl_opt_GetReal(go, "--tau");
  float  size_limit= esl_opt_GetReal(go, "--mxsize");
  char  *modestr   = esl_opt_GetString(go, "--mode");
  char  *cmfile    = esl_opt_GetArg(go, 1);
  char  *seqfile   = esl_opt_GetArg(go, 2);
  int    pass_idx  = PLI_PASS_5P_AND_3P_ANY;
  char   mode = TRMODE_J;
  if      (modestr[0]=='L') mode = TRMODE_L;
  else if (modestr[0]=='R') mode = TRMODE_R;
  else if (modestr[0]=='T') mode = TRMODE_T;
  else                      mode = TRMODE_J;

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

  ESL_SQFILE *sqfp = NULL;
  ESL_SQ *sq = esl_sq_CreateDigital(abc);
  if (esl_sqfile_OpenDigital(abc, seqfile, eslSQFILE_UNKNOWN, NULL, &sqfp) != eslOK) cm_Fail("open seqfile %s", seqfile);
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK) {
    int L = sq->n; ESL_DSQ *dsq = sq->dsq;
    if (cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, dsq, 1, L, cm->cp9b, FALSE, pass_idx, 0) != eslOK)
      cm_Fail("cp9_Seq2Bands: %s", errbuf);

    char dnc_mode = TRMODE_UNKNOWN; Parsetree_t *tr_dnc = NULL;
    float sc_dnc = TrCYKDivideAndConquerHB(cm, dsq, L, 0, 1, L, pass_idx, mode, &dnc_mode, &tr_dnc, cm->cp9b);

    CM_TR_HB_MX *mx = cm_tr_hb_mx_Create(cm);
    CM_TR_HB_SHADOW_MX *shmx = cm_tr_hb_shadow_mx_Create(cm);
    Parsetree_t *tr_or = NULL; char or_mode = TRMODE_UNKNOWN; float sc_or = IMPOSSIBLE;
    status = cm_TrAlignHB(cm, errbuf, dsq, L, size_limit, mode, pass_idx, FALSE, FALSE,
                          mx, shmx, NULL, NULL, NULL, NULL, &tr_or, &or_mode, NULL, &sc_or);
    if (status != eslOK) cm_Fail("oracle forced mode failed seq %s: %s", sq->name, errbuf);

    int peq = ParsetreeCompare(tr_or, tr_dnc);
    printf("### seq=%s L=%d mode=%c sc_dnc=%.6f sc_or=%.6f dev=%.3e ParsetreeCompare=%s nodes_dnc=%d nodes_or=%d\n",
           sq->name, L, "JLRT?"[(int)mode], sc_dnc, sc_or, fabsf(sc_dnc-sc_or), peq?"IDENTICAL":"DIFFER",
           tr_dnc?tr_dnc->n:-1, tr_or?tr_or->n:-1);
    printf("===== D&C parsetree =====\n");
    if (tr_dnc) ParsetreeDump(stdout, tr_dnc, cm, dsq);
    printf("===== ORACLE parsetree =====\n");
    if (tr_or) ParsetreeDump(stdout, tr_or, cm, dsq);

    if (tr_dnc) FreeParsetree(tr_dnc);
    if (tr_or)  FreeParsetree(tr_or);
    cm_tr_hb_mx_Destroy(mx); cm_tr_hb_shadow_mx_Destroy(shmx);
    esl_sq_Reuse(sq);
  }
  esl_sq_Destroy(sq); esl_sqfile_Close(sqfp);
  FreeCM(cm); esl_alphabet_Destroy(abc); esl_getopts_Destroy(go);
  return 0;
}

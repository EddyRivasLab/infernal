#include <esl_config.h>
#include <p7_config.h>
#include <stdio.h>
#include <stdlib.h>
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "hmmer.h"
#include "infernal.h"

int main(int argc, char **argv) {
  char errbuf[eslERRBUFSIZE];
  CM_FILE *cmfp = NULL; ESL_ALPHABET *abc = NULL; CM_t *cm = NULL;
  cm_file_Open(argv[1], NULL, TRUE, &cmfp, errbuf);
  cm_file_Read(cmfp, TRUE, &abc, &cm);
  cm_file_Close(cmfp);
  cm->config_opts |= CM_CONFIG_TRUNC;
  cm->align_opts  |= CM_ALIGN_HBANDED;
  cm->config_opts |= CM_CONFIG_LOCAL; cm->config_opts |= CM_CONFIG_HMMLOCAL; cm->config_opts |= CM_CONFIG_HMMEL;
  cm->tau = 1e-7;
  cm_Configure(cm, errbuf, -1);
  init_ilogsum(); FLogsumInit();

  ESL_SQFILE *sqfp = NULL;
  ESL_SQ *sq = esl_sq_CreateDigital(abc);
  esl_sqfile_OpenDigital(abc, argv[2], eslSQFILE_UNKNOWN, NULL, &sqfp);
  esl_sqio_Read(sqfp, sq);
  int L = sq->n;
  ESL_DSQ *dsq = sq->dsq;
  int pass_idx = PLI_PASS_5P_AND_3P_ANY;
  cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, dsq, 1, L, cm->cp9b, FALSE, pass_idx, 0);

  CM_TR_HB_MX *mx = cm_tr_hb_mx_Create(cm);
  CM_TR_HB_SHADOW_MX *shmx = cm_tr_hb_shadow_mx_Create(cm);
  Parsetree_t *tr = NULL; char mode = TRMODE_UNKNOWN; float sc = IMPOSSIBLE;
  int status = cm_TrAlignHB(cm, errbuf, dsq, L, 40000., TRMODE_J, pass_idx,
                         FALSE, FALSE, mx, shmx, NULL, NULL, NULL, NULL, &tr, &mode, NULL, &sc);
  printf("FORCED J-mode mono: status=%d sc=%f mode=%d (nodes=%d) trpenalty=%f\n", status, sc, mode, tr->n, tr->trpenalty);
  ParsetreeDump(stdout, tr, cm, dsq);
  return 0;
}

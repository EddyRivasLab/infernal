/* r092_oracle_row.c : dump oracle forced-L inside-HB Ldp[v][L][d] for ALL d at j=L,
 * for chosen states. Usage: r092_oracle_row <cmfile> <seqfile> v1 [v2 ...] */
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
  setbuf(stdout, NULL);
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
  int b; char mode; float sc;
  cm_TrCYKInsideAlignHB(cm, errbuf, dsq, L, 1e9, TRMODE_L, pass_idx, mx, shmx, &b, &mode, &sc);
  CP9Bands_t *cp9b = cm->cp9b;
  int j = L;
  for (int a = 3; a < argc; a++) {
    int v = atoi(argv[a]);
    if (j < cp9b->jmin[v] || j > cp9b->jmax[v]) { printf("v=%d j=%d OUT-OF-JBAND\n", v, j); continue; }
    int jp_v = j - cp9b->jmin[v];
    int hdmn = cp9b->hdmin[v][jp_v], hdmx = cp9b->hdmax[v][jp_v];
    printf("# v=%d j=%d dband[%d..%d]\n", v, j, hdmn, hdmx);
    for (int d = hdmn; d <= hdmx; d++) {
      int dp = d - hdmn;
      float lv = (mx->Ldp && mx->Ldp[v]) ? mx->Ldp[v][jp_v][dp] : -9999;
      printf("  v=%d d=%d Ldp=%f\n", v, d, lv);
    }
  }
  return 0;
}

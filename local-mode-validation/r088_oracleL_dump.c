/* r088_oracleL_dump.c : dump oracle forced-L inside-HB Ldp[v][L][L] (+ Jdp) for
 * chosen states, to compare against the D&C's begin-at-471 inside (claimed 9.403).
 * Usage: r088_oracleL_dump <cmfile> <seqfile> [v1 v2 ...] */
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

static void dumpcell(CM_t *cm, CM_TR_HB_MX *mx, int v, int j, int d) {
  CP9Bands_t *cp9b = cm->cp9b;
  if (j < cp9b->jmin[v] || j > cp9b->jmax[v]) { printf("  state %d [j=%d d=%d]: OUT-OF-JBAND (jmin=%d jmax=%d)\n", v,j,d, cp9b->jmin[v], cp9b->jmax[v]); return; }
  int jp_v = j - cp9b->jmin[v];
  int hdmn = cp9b->hdmin[v][jp_v], hdmx = cp9b->hdmax[v][jp_v];
  if (d < hdmn || d > hdmx) { printf("  state %d [j=%d d=%d]: OUT-OF-DBAND (dband@j[%d..%d])\n", v,j,d, hdmn, hdmx); return; }
  int dp_v = d - hdmn;
  float lv = (mx->Ldp && mx->Ldp[v]) ? mx->Ldp[v][jp_v][dp_v] : -9999;
  float jv = (mx->Jdp && mx->Jdp[v]) ? mx->Jdp[v][jp_v][dp_v] : -9999;
  printf("  state %d [j=%d d=%d]: Jdp=%f  Ldp=%f  (jband[%d..%d] dband@j[%d..%d] sttype=%d)\n",
         v, j, d, jv, lv, cp9b->jmin[v],cp9b->jmax[v], hdmn, hdmx, cm->sttype[v]);
}

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
  int status = cm_TrCYKInsideAlignHB(cm, errbuf, dsq, L, 1e9, TRMODE_L, pass_idx, mx, shmx, &b, &mode, &sc);
  printf("oracle inside-HB forced-L: status=%d sc=%f b=%d mode=%d L=%d\n", status, sc, b, mode, L);
  if (status != eslOK) { printf("ERR: %s\n", errbuf); return 1; }

  int i;
  printf("--- full-span cells (j=L=%d, d=L=%d) ---\n", L, L);
  for (i = 3; i < argc; i++) { int v = atoi(argv[i]); dumpcell(cm, mx, v, L, L); }
  return 0;
}

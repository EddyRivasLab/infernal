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

/* dump oracle inside-HB Jalpha[v][j][d] for a requested (v,j,d) */
static void dumpcell(CM_t *cm, CM_TR_HB_MX *mx, const char *tag, int v, int j, int d) {
  CP9Bands_t *cp9b = cm->cp9b;
  if (j < cp9b->jmin[v] || j > cp9b->jmax[v]) { printf("%s Jalpha[%d][%d][%d] = OUT-OF-JBAND (jmin=%d jmax=%d)\n", tag, v,j,d, cp9b->jmin[v], cp9b->jmax[v]); return; }
  int jp_v = j - cp9b->jmin[v];
  int hdmn = cp9b->hdmin[v][jp_v], hdmx = cp9b->hdmax[v][jp_v];
  if (d < hdmn || d > hdmx) { printf("%s Jalpha[%d][%d][%d] = OUT-OF-DBAND (jband[%d..%d] dband@j[%d..%d])\n", tag, v,j,d, cp9b->jmin[v],cp9b->jmax[v], hdmn, hdmx); return; }
  int dp_v = d - hdmn;
  printf("%s Jalpha[%d][%d][%d] = %f  (jband[%d..%d] dband@j[%d..%d])\n", tag, v, j, d, mx->Jdp[v][jp_v][dp_v], cp9b->jmin[v],cp9b->jmax[v], hdmn, hdmx);
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
  int status = cm_TrCYKInsideAlignHB(cm, errbuf, dsq, L, 1e9, TRMODE_J, pass_idx, mx, shmx, &b, &mode, &sc);
  printf("oracle inside-HB forced-J: status=%d sc=%f b=%d mode=%d\n", status, sc, b, mode);
  if (status != eslOK) { printf("ERR: %s\n", errbuf); return 1; }

  /* the D&C winning split cells */
  dumpcell(cm, mx, "[bif   ]", 54, 69, 60);
  dumpcell(cm, mx, "[Lchild]", 121, 41, 32);
  dumpcell(cm, mx, "[Rchild]", 55, 69, 28);
  /* the claimed entry state 3, full span */
  dumpcell(cm, mx, "[entry3]", 3, L, L);
  dumpcell(cm, mx, "[st0   ]", 0, L, L);
  printf("--- inside wedge path (state:j:d) ---\n");
  dumpcell(cm, mx, "w", 3,  77, 77);
  dumpcell(cm, mx, "w", 6,  76, 76);
  dumpcell(cm, mx, "w", 12, 75, 74);
  dumpcell(cm, mx, "w", 18, 74, 72);
  dumpcell(cm, mx, "w", 24, 73, 70);
  dumpcell(cm, mx, "w", 30, 72, 68);
  dumpcell(cm, mx, "w", 36, 71, 66);
  dumpcell(cm, mx, "w", 42, 70, 64);
  dumpcell(cm, mx, "w", 48, 69, 62);
  dumpcell(cm, mx, "w", 51, 69, 61);
  dumpcell(cm, mx, "w", 54, 69, 60);
  return 0;
}

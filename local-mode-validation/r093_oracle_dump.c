/* r093_oracle_dump.c : dump oracle forced-R inside (Jdp/Ldp/Rdp) and outside (Rdp)
 * cells, to compare against the D&C generic-splitter term decomposition.
 * Usage: r093_oracle_dump <cmfile> <seqfile> <mode J|L|R|T>
 */
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

static float getcell(CM_t *cm, float ***dp, int v, int j, int d) {
  CP9Bands_t *cp9b = cm->cp9b;
  if (dp==NULL || dp[v]==NULL) return 111111.0;
  if (j < cp9b->jmin[v] || j > cp9b->jmax[v]) return 222222.0;
  int jp_v = j - cp9b->jmin[v];
  int hdmn = cp9b->hdmin[v][jp_v], hdmx = cp9b->hdmax[v][jp_v];
  if (d < hdmn || d > hdmx) return 333333.0;
  return dp[v][jp_v][d - hdmn];
}

int main(int argc, char **argv) {
  char errbuf[eslERRBUFSIZE];
  setbuf(stdout, NULL);
  CM_FILE *cmfp = NULL; ESL_ALPHABET *abc = NULL; CM_t *cm = NULL;
  cm_file_Open(argv[1], NULL, TRUE, &cmfp, errbuf);
  cm_file_Read(cmfp, TRUE, &abc, &cm);
  cm_file_Close(cmfp);
  char mode = TRMODE_R;
  if (argc>3) { if(argv[3][0]=='L')mode=TRMODE_L; else if(argv[3][0]=='J')mode=TRMODE_J; else if(argv[3][0]=='T')mode=TRMODE_T; else mode=TRMODE_R; }
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
  int L = sq->n; ESL_DSQ *dsq = sq->dsq;
  int pass_idx = PLI_PASS_5P_AND_3P_ANY;
  cp9_Seq2Bands(cm, errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, dsq, 1, L, cm->cp9b, FALSE, pass_idx, 0);

  CM_TR_HB_MX *imx  = cm_tr_hb_mx_Create(cm);
  CM_TR_HB_MX *omx  = cm_tr_hb_mx_Create(cm);
  CM_TR_HB_SHADOW_MX *shmx = cm_tr_hb_shadow_mx_Create(cm);
  int b; char rmode; float sc;
  int status = cm_TrCYKInsideAlignHB(cm, errbuf, dsq, L, 1e9, mode, pass_idx, imx, shmx, &b, &rmode, &sc);
  printf("oracle inside-HB forced=%c: status=%d sc=%f b=%d rmode=%c\n", "TRLJ?"[(int)mode], status, sc, b, "TRLJ?"[(int)rmode]);
  status = cm_TrCYKOutsideAlignHB(cm, errbuf, dsq, L, 1e9, mode, pass_idx, TRUE, omx, imx);
  printf("oracle outside-HB (do_check): status=%d errbuf='%s'\n", status, errbuf);

  /* cells: either "-f <file>" (lines: plane v j d) or v,j,d quads on cmdline */
  if (argc>=6 && argv[4][0]=='-' && argv[4][1]=='f') {
    FILE *fp = fopen(argv[5],"r"); char pl; int v,j,d;
    while (fscanf(fp," %c %d %d %d", &pl,&v,&j,&d)==4) {
      float ***idp = (pl=='J')?imx->Jdp:(pl=='L')?imx->Ldp:(pl=='R')?imx->Rdp:imx->Tdp;
      float iv = getcell(cm, idp, v, j, d);
      printf("cell plane=%c v=%d j=%d d=%d  INSIDE=%.5f\n", pl, v, j, d, iv);
    }
    fclose(fp);
    return 0;
  }
  int a;
  for (a = 4; a+3 < argc; a += 4) {
    char pl = argv[a][0];
    int v = atoi(argv[a+1]), j = atoi(argv[a+2]), d = atoi(argv[a+3]);
    float iv, ov;
    float ***idp = (pl=='J')?imx->Jdp:(pl=='L')?imx->Ldp:(pl=='R')?imx->Rdp:imx->Tdp;
    float ***odp = (pl=='J')?omx->Jdp:(pl=='L')?omx->Ldp:(pl=='R')?omx->Rdp:omx->Tdp;
    iv = getcell(cm, idp, v, j, d);
    ov = getcell(cm, odp, v, j, d);
    printf("cell plane=%c v=%d j=%d d=%d  INSIDE=%.5f  OUTSIDE=%.5f\n", pl, v, j, d, iv, ov);
  }
  return 0;
}

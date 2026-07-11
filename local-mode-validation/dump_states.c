#include <esl_config.h>
#include <p7_config.h>
#include <stdio.h>
#include "easel.h"
#include "esl_alphabet.h"
#include "hmmer.h"
#include "infernal.h"
int main(int argc, char **argv) {
  char errbuf[eslERRBUFSIZE];
  CM_FILE *cmfp = NULL; ESL_ALPHABET *abc = NULL; CM_t *cm = NULL;
  cm_file_Open(argv[1], NULL, TRUE, &cmfp, errbuf);
  cm_file_Read(cmfp, TRUE, &abc, &cm);
  cm_file_Close(cmfp);
  for (int v = 0; v <= 60; v++) {
    printf("v=%d nd=%d ndtype=%d sttype=%d stid=%d plast=%d pnum=%d cfirst=%d cnum=%d\n",
           v, cm->ndidx[v], cm->ndtype[cm->ndidx[v]], cm->sttype[v], cm->stid[v], cm->plast[v], cm->pnum[v], cm->cfirst[v], cm->cnum[v]);
  }
  return 0;
}

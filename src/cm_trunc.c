/* CM_TR_OPTS and CM_TR_PENALTIES: data structures with information
 * relevant to truncated DP functions.
 * 
 * Contents:
 *    1. CM_TR_SINFO data structure functions,
 *       per-scan info for truncated CM scanners
 *    2. CM_TR_PENALTIES data structure functions,
 *       info on truncated alignment score penalties.
 *
 * EPN, Sat Jan 21 14:02:26 2012
 */
#include "esl_config.h"
#include "p7_config.h"
#include "config.h"

#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "easel.h"
#include "hmmer.h"
#include "funcs.h"
#include "structs.h"

/*****************************************************************
 *   1. CM_TR_OPTS data structure functions,
 *      truncated alignment options
 *****************************************************************/

/* Function: cm_tr_opts_Create()
 * Date:     EPN, Tue Nov  8 08:27:16 2011
 *
 * Purpose:  Allocate and initialize a CM_TR_OPTS object.
 * 
 * Returns:  Newly allocated CM_TR_OPTS object. NULL if out
 *           of memory. Needs a <cm> to create truncation
 *           penalties.
 */
CM_TR_OPTS *
cm_tr_opts_Create()
{
  int status;
  CM_TR_OPTS *tro = NULL;
  ESL_ALLOC(tro, sizeof(CM_TR_OPTS));

  /* set default values */
  tro->allow_L     = TRUE; /* L mode alignments are allowed */
  tro->allow_R     = TRUE; /* R mode alignments are allowed */
  tro->force_i0_RT = TRUE; /* in R and T alignments i0 must be included */
  tro->force_j0_LT = TRUE; /* in L and T alignments j0 must be included */

  return tro;

 ERROR:
  if(tro != NULL) free(tro);
  return NULL;
}

/* Function: cm_tr_opts_PenaltyIdx()
 * Date:     EPN, Sun Jan 22 16:03:21 2012
 *
 * Purpose:  Return the appropriate truncation 
 *           penalty index given truncation options.
 *           in <tro>.
 * 
 * Returns:  truncation penalty index, either
 *           TRPENALTY_5P_OR_3P or TRPENALTY_5P_AND_3P
 */
int
cm_tr_opts_PenaltyIdx(CM_TR_OPTS *tro)
{
  if(tro->allow_R && tro->allow_L) return TRPENALTY_5P_AND_3P;
  else                             return TRPENALTY_5P_OR_3P;
}

/* Function: cm_tr_opts_Dump()
 * Date:     EPN, Sat Jan 21 13:47:49 2012
 *
 * Purpose:  Print contents of a CM_TR_OPTS to
 *           stream <fp> for inspection.
 * 
 * Returns:  void
 */
void
cm_tr_opts_Dump(FILE *fp, CM_TR_OPTS *tro)
{
  fprintf(fp, "CM_TR_OPTS dump\n");
  fprintf(fp, "------------------\n");
  fprintf(fp, "allow_R     = %s\n",  tro->allow_R     ? "TRUE" : "FALSE");
  fprintf(fp, "allow_L     = %s\n",  tro->allow_L     ? "TRUE" : "FALSE");
  fprintf(fp, "force_i0_RT = %s\n",  tro->force_i0_RT ? "TRUE" : "FALSE");
  fprintf(fp, "force_j0_LT = %s\n",  tro->force_j0_LT ? "TRUE" : "FALSE");

  return;
}

/*****************************************************************
 *   2. CM_TR_PENALTIES data structure functions,
 *       info on truncated alignment score penalties.
 *****************************************************************/

/* Function: cm_tr_penalties_Create()
 * Date:     EPN, Sat Jan 21 12:03:52 2012
 *
 * Purpose:  Allocate and initialize a CM_TR_PENALTIES object.
 *           A CM and its emit map are required to determine
 *           truncation penalty scores.
 * 
 * Returns:  Newly allocated CM_TR_PENALTIES object. NULL if out
 *           of memory.
 */
CM_TR_PENALTIES *
cm_tr_penalties_Create(CM_t *cm)
{
  int status;
  int v, nd;
  int subtree_clen;

  if(cm == NULL || cm->emap == NULL) goto ERROR;

  CM_TR_PENALTIES *trp = NULL;
  ESL_ALLOC(trp, sizeof(CM_TR_PENALTIES));

  trp->M = cm->M;

  /* allocate and initialized the penalty arrays */
  ESL_ALLOC(trp->ptyAA,  sizeof(float *) * NTRPENALTY); 
  ESL_ALLOC(trp->iptyAA, sizeof(int *)   * NTRPENALTY);

  trp->ptyAA[TRPENALTY_5P_OR_3P]   = NULL;
  trp->ptyAA[TRPENALTY_5P_AND_3P]  = NULL;
  trp->iptyAA[TRPENALTY_5P_OR_3P]  = NULL;
  trp->iptyAA[TRPENALTY_5P_AND_3P] = NULL;
  ESL_ALLOC(trp->ptyAA[TRPENALTY_5P_OR_3P],   sizeof(float) * cm->M);
  ESL_ALLOC(trp->ptyAA[TRPENALTY_5P_AND_3P],  sizeof(float) * cm->M);
  ESL_ALLOC(trp->iptyAA[TRPENALTY_5P_OR_3P],  sizeof(int)   * cm->M);
  ESL_ALLOC(trp->iptyAA[TRPENALTY_5P_AND_3P], sizeof(int)   * cm->M);

  for(v = 0; v < cm->M; v++) { 
    if(StateDelta(cm->sttype[v]) > 0 || cm->stid[v] == BIF_B || cm->stid[v] == ROOT_S) { 
      /* truncated alignments can only begin at any emitter, B_st or ROOT_S */
      nd = cm->ndidx[v];
      subtree_clen = cm->emap->rpos[nd]-cm->emap->lpos[nd]+1;
      if(cm->ndtype[nd] != MATP_nd && cm->ndtype[nd] != MATL_nd) subtree_clen--; /* lpos was one less than what we want */
      if(cm->ndtype[nd] != MATP_nd && cm->ndtype[nd] != MATR_nd) subtree_clen--; /* rpos was one more than what we want */
      trp->ptyAA[TRPENALTY_5P_AND_3P][v] = sreLOG2(2. / (subtree_clen * (subtree_clen+1))); /* penalty for 5' *and* 3' truncation */
      trp->ptyAA[TRPENALTY_5P_OR_3P][v]  = sreLOG2(1. / subtree_clen);                      /* penalty for 5' *or*  3' truncation */
      trp->iptyAA[TRPENALTY_5P_AND_3P][v] = Prob2Score(2. / (subtree_clen * (subtree_clen+1)), 1.0); /* penalty for 5' *and* 3' truncation */
      trp->iptyAA[TRPENALTY_5P_OR_3P][v]  = Prob2Score(1. / subtree_clen, 1.0);                      /* penalty for 5' *or*  3' truncation */
    }
    else { 
      trp->ptyAA[TRPENALTY_5P_AND_3P][v] = IMPOSSIBLE;
      trp->ptyAA[TRPENALTY_5P_OR_3P][v]  = IMPOSSIBLE;
      trp->iptyAA[TRPENALTY_5P_AND_3P][v] = -INFTY;
      trp->iptyAA[TRPENALTY_5P_OR_3P][v]  = -INFTY;
    }
  }
  /*cm_tr_penalties_Dump(stdout, cm, trp);*/

  return trp;

 ERROR:
  if(trp != NULL) cm_tr_penalties_Destroy(trp);
  return NULL;
}

/* Function: cm_tr_penalties_Sizeof()
 * Date:     EPN, Sat Jan 21 15:37:58 2012
 *
 * Purpose:  Calculate and return the size of a CM_TR_PENALTIES
 *           object in Mb.
 */

float 
cm_tr_penalties_Sizeof(CM_TR_PENALTIES *trp)
{
  float bytes = 0.;

  if(trp == NULL) return 0.;
  
  bytes = sizeof(CM_TR_PENALTIES);
  bytes += sizeof(float *) * NTRPENALTY; /* ptyAA, 1st dim */
  bytes += sizeof(int *)   * NTRPENALTY; /* iptyAA, 1st dim */
  bytes += sizeof(float) * NTRPENALTY * trp->M; /* ptyAA, 2nd dim */
  bytes += sizeof(int)   * NTRPENALTY * trp->M; /* ptyAA, 2nd dim */

  return bytes / 1000000.;
}

/* Function:  cm_tr_penalties_Dump()
 *
 * Purpose:   Print contents of the <CM_TR_PENALTIES> <trp> to
 *            stream <fp> for inspection. 
 *
 * Returns:   void
 */
void
cm_tr_penalties_Dump(FILE *fp, const CM_t *cm, const CM_TR_PENALTIES *trp)
{
  int v, nd, subtree_clen;

  fprintf(fp, "CM_TR_PENALTIES dump\n");
  fprintf(fp, "--------------------\n");
  fprintf(fp, "M = %d\n",  trp->M);
  fprintf(fp, "\n");
  fprintf(fp, "%5s  %5s  %7s  %7s  %10s  %10s  %10s  %10s\n", "stidx", "ndidx", "stid", "st_clen", "f5P_OR_3P", "f5P_AND_3P", "i5P_OR_3P", "i5P_AND_3P");
  fprintf(fp, "%5s  %5s  %7s  %7s  %10s  %10s  %10s  %10s\n", "-----", "-----", "-------", "-------", "----------", "----------", "----------", "----------");
  for(v = 0; v < cm->M; v++) { 
    nd = cm->ndidx[v];
    subtree_clen = cm->emap->rpos[nd]-cm->emap->lpos[nd]+1;
    if(cm->ndtype[nd] != MATP_nd && cm->ndtype[nd] != MATL_nd) subtree_clen--; /* lpos was one less than what we want */
    if(cm->ndtype[nd] != MATP_nd && cm->ndtype[nd] != MATR_nd) subtree_clen--; /* rpos was one more than what we want */
    if(NOT_IMPOSSIBLE(trp->ptyAA[TRPENALTY_5P_OR_3P][v]) && NOT_IMPOSSIBLE(trp->ptyAA[TRPENALTY_5P_AND_3P][v])) { 
      fprintf(fp, "%5d  %5d  %-7s  %7d  %10.3f  %10.3f  %10d  %10d\n", v, cm->ndidx[v], CMStateid(cm->stid[v]), subtree_clen, 
	      trp->ptyAA[TRPENALTY_5P_AND_3P][v], 
	      trp->ptyAA[TRPENALTY_5P_OR_3P][v],
	      trp->iptyAA[TRPENALTY_5P_AND_3P][v], 
	      trp->iptyAA[TRPENALTY_5P_OR_3P][v]);
    }
    else { 
      fprintf(fp, "%5d  %5d  %-7s  %7d  %10s  %10s  %10s  %10s\n", v, cm->ndidx[v], CMStateid(cm->stid[v]), subtree_clen, 
	      "IMPOSSBLE", "IMPOSSBLE", "-INFTY", "-INFTY");
    }
  }

  return;
}

/* Function: cm_tr_penalties_Destroy()
 * Date:     EPN, Sat Jan 21 10:30:53 2012
 *
 * Purpose:  Destroy a CM_TR_PENALTIES object.
 * 
 * Returns:  void
 */
void
cm_tr_penalties_Destroy(CM_TR_PENALTIES *trp)
{
  int i;

  if(trp == NULL) return;
  if(trp->ptyAA != NULL) { 
    for(i = 0; i < NTRPENALTY; i++) { 
      if(trp->ptyAA[i] != NULL) free(trp->ptyAA[i]);
    }
    free(trp->ptyAA);
  }
  if(trp->iptyAA != NULL) { 
    for(i = 0; i < NTRPENALTY; i++) { 
      if(trp->iptyAA[i] != NULL) free(trp->iptyAA[i]);
    }
    free(trp->iptyAA);
  }

  free(trp);
  trp = NULL;
  
  return;
}

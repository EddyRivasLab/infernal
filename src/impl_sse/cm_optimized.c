#include "impl_sse.h"

#define SCALE_W 500.0	/* set a default scaling factor for 16-bit int scores */

/* External API 
 * int cm_optimized_Convert(const CM_t *cm, CM_OPTIMIZED *ocm);
 * void cm_optimized_Free(CM_OPTIMIZED *ocm);
 */

/* Internal funcs */
static int16_t wordify(float scale_w, float sc);

/* Fuction:  cm_optimized_Convert()
 * Author:   DLK, Tue May 05 2009
 */
CM_OPTIMIZED*
cm_optimized_Convert(const CM_t *cm)
{
  int status;
  int v, y, yoffset;
  int npairs = 0;
  int nsinglets = 0;
  int16_t *oesc_start;
  int cur_cell;
  CM_OPTIMIZED *ocm;

  ESL_ALLOC(ocm, sizeof(CM_OPTIMIZED));

  /* Copy basic information */
  ocm->M = cm->M;
  ocm->sttype = cm->sttype;	/* careful, these point at the external resource */
  ocm->cfirst = cm->cfirst;
  ocm->cnum   = cm->cnum;
  ocm->abc    = cm->abc;

  /* count pairs, singlets */
  for(v = 0; v < cm->M; v++) {
    switch(ocm->sttype[v]) {
    case IL_st:
    case ML_st:
    case IR_st:
    case MR_st:
      nsinglets++;
      break;
    case MP_st:
      npairs++;
      break;
    }
  }
  
  /* Allocate memory */
  ESL_ALLOC(ocm->tsc,     sizeof(int16_t *) * ocm->M);
  ESL_ALLOC(ocm->oesc,    sizeof(int16_t *) * ocm->M);
  ESL_ALLOC(ocm->beginsc, sizeof(int16_t  ) * ocm->M);
  ESL_ALLOC(ocm->endsc,   sizeof(int16_t  ) * ocm->M);

  ESL_ALLOC(ocm->tsc[0],  sizeof(int16_t  ) * ocm->M * MAXCONNECT);
  ESL_ALLOC(oesc_start,   sizeof(int16_t  ) * (nsinglets * ocm->abc->Kp + npairs * ocm->abc->Kp * ocm->abc->Kp));
  for (v = 1; v<ocm->M; v++) {
    //ocm->tsc[v]  = ocm->tsc[v-1]  + sizeof(int16_t) * MAXCONNECT;
    ocm->tsc[v]  = ocm->tsc[v-1]  + MAXCONNECT;
  }
  cur_cell = 0;
  for (v = 0; v<ocm->M; v++) {
    switch(ocm->sttype[v]) {
      case ML_st:
      case IL_st:
      case MR_st:
      case IR_st:
        ocm->oesc[v] = oesc_start + cur_cell;
        cur_cell += ocm->abc->Kp;
        break;
      case MP_st:
        ocm->oesc[v] = oesc_start + cur_cell;
        cur_cell += ocm->abc->Kp * ocm->abc->Kp;
        break;
      default:
        ocm->oesc[v] = NULL;
        break;
    }
  }

  /* Scale score values */
  ocm->scale_w = SCALE_W;
  for (v = 0; v<ocm->M; v++) {
    y = ocm->cfirst[v];
    for (yoffset = 0; y < ocm->cnum[v]; y++) {
      ocm->tsc[v][y] = wordify(ocm->scale_w, cm->tsc[v][y]);
    }

    if (ocm->sttype[v] == MP_st) {
      for (y = 0; y < (ocm->abc->Kp * ocm->abc->Kp); y++) {
        ocm->oesc[v][y] = wordify(ocm->scale_w, cm->oesc[v][y]);
      }
    }
    else if (ocm->sttype[v] == ML_st || ocm->sttype[v] == IL_st ||
             ocm->sttype[v] == MR_st || ocm->sttype[v] == IR_st ) {
      for (y = 0; y < ocm->abc->Kp; y++) {
        ocm->oesc[v][y] = wordify(ocm->scale_w, cm->oesc[v][y]);
      }
    }

    ocm->beginsc[v] = wordify(ocm->scale_w, cm->beginsc[v]);
    ocm->endsc[v]   = wordify(ocm->scale_w, cm->endsc[v]);
  }
  ocm->el_selfsc = wordify(ocm->scale_w, cm->el_selfsc);

  return ocm;

ERROR: 
  cm_Fail("Memory allocation error.\n");
  return 0; /* never reached */
}

/* Function:  cm_optimized_Free()
 * Author:    DLK, Tue May 05 2009
 */
void
cm_optimized_Free(CM_OPTIMIZED *ocm)
{
  int v = 0;

  if (ocm->tsc[0] != NULL) free(ocm->tsc[0]);
  while (ocm->oesc[v] == NULL && v < ocm->M) { v++; }
  if (ocm->oesc[v] != NULL) free(ocm->oesc[v]);

  if (ocm->endsc != NULL) free(ocm->endsc);
  if (ocm->beginsc != NULL) free(ocm->beginsc);
  if (ocm->oesc != NULL) free(ocm->oesc);
  if (ocm->tsc != NULL) free(ocm->tsc);

  return;
}

/* wordify()
 * Converts log probability score to a rounded signed 16-bit integer cost.
 * Both emissions and transitions for ViterbiFilter get this treatment.
 * No bias term needed, because we use signed words. 
 *   e.g. a score of +3.2, with scale 500.0, becomes +1600.
 */
static int16_t
wordify(float scale_w, float sc)
{
  sc  = roundf(scale_w * sc);
  if      (sc >=  32767.0) return  32767;
  else if (sc <= -32768.0) return -32768;
  else return (int16_t) sc;
}


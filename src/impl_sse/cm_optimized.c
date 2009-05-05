#include "impl_sse.h"

#define SCALE_W 500.0	/* set a default scaling factor for 16-bit int scores */

/* External API 
 * int cm_optimized_Convert(const CM_t *cm, CM_OPTIMIZED *ocm);
 */

/* Internal funcs */
static int16_t wordify(float scale_w, float sc);

/* Fuction:  cm_optimized_Convert()
 * Author:   DLK, Tue May 05 2009
 */
int
cm_optimized_Convert(const CM_t *cm, CM_OPTIMIZED *ocm)
{
  int status;
  int v, y, yoffset;

  if (ocm == NULL)
    ESL_ALLOC(ocm, sizeof(ocm));
  else
    cm_Fail("Not set up to re-use CM_OPTIMIZED models!\n");

  /* Copy basic information */
  ocm->M = cm->M;
  ocm->sttype = cm->sttype;	/* careful, these point at the external resource */
  ocm->cfirst = cm->cfirst;
  ocm->cnum   = cm->cnum;
  ocm->abc    = cm->abc;
  
  /* Allocate memory */
  ESL_ALLOC(ocm->tsc,     sizeof(int16_t *) * ocm->M);
  ESL_ALLOC(ocm->oesc,    sizeof(int16_t *) * ocm->M);
  ESL_ALLOC(ocm->beginsc, sizeof(int16_t  ) * ocm->M);
  ESL_ALLOC(ocm->endsc,   sizeof(int16_t  ) * ocm->M);

  ESL_ALLOC(ocm->tsc[0],  sizeof(int16_t  ) * ocm->M * MAXCONNECT);
  ESL_ALLOC(ocm->oesc[0], sizeof(int16_t  ) * ocm->M * ocm->abc->Kp);
  for (v = 1; v<ocm->M; v++) {
    ocm->tsc[v]  = ocm->tsc[v-1]  + sizeof(int16_t) * MAXCONNECT;
    ocm->oesc[v] = ocm->oesc[v-1] + sizeof(int16_t) * ocm->abc->Kp;
  }

  /* Scale score values */
  ocm->scale_w = SCALE_W;
  for (v = 0; v<ocm->M; v++) {
    y = ocm->cfirst[v];
    for (yoffset = 0; y < ocm->cnum[v]; y++) {
      ocm->tsc[v][y] = wordify(ocm->scale_w, cm->tsc[v][y]);
    }

    for (y = 0; y < ocm->abc->Kp; y++) {
      ocm->oesc[v][y] = wordify(ocm->scale_w, cm->oesc[v][y]);
    }

    ocm->beginsc[v] = wordify(ocm->scale_w, cm->beginsc[v]);
    ocm->endsc[v]   = wordify(ocm->scale_w, cm->endsc[v]);
  }
  ocm->el_selfsc = wordify(ocm->scale_w, cm->el_selfsc);

  return 0;

ERROR: 
  cm_Fail("Memory allocation error.\n");
  return 0; /* never reached */
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


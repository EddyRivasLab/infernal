#include "impl_sse.h"
#include "esl_stack.h"
#include "esl_vectorops.h"

#define SCALE_W 500.0	/* set a default scaling factor for 16-bit int scores */

/* External API 
 * int cm_optimized_Convert(const CM_t *cm, CM_OPTIMIZED *ocm);
 * void cm_optimized_Free(CM_OPTIMIZED *ocm);
 * int16_t wordify(float scale_w, float sc);
 */

/* Internal funcs */

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
  ocm->plast  = cm->plast;
  ocm->pnum   = cm->pnum;
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
    for (yoffset = 0; yoffset < ocm->cnum[v]; yoffset++) {
      ocm->tsc[v][yoffset] = wordify(ocm->scale_w, cm->tsc[v][yoffset]);
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
  return NULL; /* never reached */
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

/* Fuction:  cm_consensus_Convert()
 * Author:   DLK, Thu May 14 2009
 */
CM_CONSENSUS*
cm_consensus_Convert(CM_t *cm)
{
  int status;
  int v, w, x, y, z;
  int offset;
  int pairs, singlets;
  int nstates, nfrags, total_frags;
  int ebases, q;
  uint8_t *oesc_next;
  float max = 0.0;
  CM_CONSENSUS *ccm = NULL;

  ESL_STACK *oldstate;
  ESL_STACK *newstate;

  oldstate = esl_stack_ICreate();
  newstate = esl_stack_ICreate();

  ESL_ALLOC(ccm, sizeof(CM_CONSENSUS));

  /* Set alphabet */
  if (cm->abc->type == eslNONSTANDARD) { cm_Fail("Failed to handle nonstandard alphabet\n"); }
  ccm->abc = esl_alphabet_Create(cm->abc->type);

  /* Set number of consensus states */
  pairs   = CMCountNodetype(cm, MATP_nd);
  singlets= CMCountNodetype(cm, MATL_nd) + CMCountNodetype(cm, MATR_nd);
  /*
  pairs   = CMCountStatetype(cm, MP_st);
  singlets= CMCountStatetype(cm, ML_st) + CMCountStatetype(cm, MR_st);
  */
  /*
  ccm->M  = 0;
  ccm->M += CMCountStatetype(cm, MP_st);
  ccm->M += CMCountStatetype(cm, ML_st);
  ccm->M += CMCountStatetype(cm, MR_st);
  */
  ccm->M  = pairs + singlets;
  ccm->M += CMCountStatetype(cm,  S_st);
  ccm->M += CMCountStatetype(cm,  B_st);
  ccm->M += CMCountStatetype(cm,  E_st);

  //FIXME - taking these defaults directly from HMMER - need to be checked!!!
  for (v = 0; v < cm->M; v++) {
    switch(cm->sttype[v]) {
      case ML_st:
      case MR_st:
        max = ESL_MAX(max, esl_vec_FMax(cm->oesc[v], cm->abc->Kp));
        break;
      case MP_st:
        max = ESL_MAX(max, esl_vec_FMax(cm->oesc[v], cm->abc->Kp*cm->abc->Kp));
        break;
      default:
        break;
    }
  }
  ccm->scale_b = 3.0 / eslCONST_LOG2;
  ccm->base_b  = 195;
  ccm->bias_b  = -1 * unbiased_byteify(ccm, max);

  /* Set state connection */
  ESL_ALLOC(ccm->next,   sizeof(int ) * ccm->M);
  ESL_ALLOC(ccm->sttype, sizeof(char) * ccm->M);
  x = v = 0;

  /* Setup emissions */
  ESL_ALLOC(ccm->oesc,    sizeof(uint8_t *) * ccm->M);
  ESL_ALLOC(ccm->oesc[0], sizeof(uint8_t  ) * (singlets * ccm->abc->Kp + pairs * ccm->abc->Kp * ccm->abc->Kp));
  oesc_next = ccm->oesc[0];

  nstates = nfrags = total_frags = 0;
  ebases = 0;
  while (v != -1) {
    ccm->sttype[x] = cm->sttype[v];
    if (cm->sttype[v] == E_st) {
      ccm->next[x] = -1;
      ccm->oesc[x] = NULL;

      nfrags = nstates*(nstates+1)/2;
      total_frags += nfrags;
      for (q = 1; q <= nstates; q++) {
        ebases += q*(nstates-q+1)*StateDelta(ccm->sttype[x-q]);
      }
      nstates = 0;

      /* pop next x and v off stack */
      if (esl_stack_IPop(oldstate, &v) == eslEOD) { v = -1; }
      if (esl_stack_IPop(newstate, &x) == eslEOD) { x = -1; }
    }
    else if (cm->sttype[v] == B_st) {
      w = cm->cfirst[v];
      y = cm->cnum[v];

      offset  = CMSubtreeCountStatetype(cm, w, MP_st);
      offset += CMSubtreeCountNodetype(cm, w, MATL_ML);
      offset += CMSubtreeCountNodetype(cm, w, MATR_MR);
      offset += CMSubtreeCountStatetype(cm, w,  S_st);
      offset += CMSubtreeCountStatetype(cm, w,  B_st);
      offset += CMSubtreeCountStatetype(cm, w,  E_st);

      ccm->next[x] = x+offset+1;
      ccm->oesc[x] = NULL;

      nfrags = nstates*(nstates+1)/2;
      total_frags += nfrags;
      for (q = 1; q <= nstates; q++) {
        ebases += q*(nstates-q+1)*StateDelta(ccm->sttype[x-q]);
      }
      nstates = 0;

      /* push y and x+offset on stacks */
      esl_stack_IPush(oldstate, y);
      esl_stack_IPush(newstate, x+offset+1);

      x++;
      v = w;
    }
    else {
      y = cm->cfirst[v];

      while (y<cm->M && (cm->sttype[y] == D_st || cm->sttype[y] == IL_st || cm->sttype[y] == IR_st || cm->stid[y] == MATP_ML || cm->stid[y] == MATP_MR)) { y++; }
      //if (y == cm->M) y = -1;

      ccm->next[x] = x+1;
      ccm->oesc[x] = oesc_next;
      if (cm->sttype[v] == MP_st) {
        nstates++;
        for (z = 0; z < (ccm->abc->Kp * ccm->abc->Kp); z++) {
          ccm->oesc[x][z] = biased_byteify(ccm, cm->oesc[v][z]);
        }
        oesc_next += ccm->abc->Kp * ccm->abc->Kp;
      }
      else if (cm->stid[v] == MATL_ML || cm->stid[v] == MATR_MR) {
        nstates++;
        for (z = 0; z < ccm->abc->Kp; z++) {
          ccm->oesc[x][z] = biased_byteify(ccm, cm->oesc[v][z]);
        }
        oesc_next += ccm->abc->Kp;
      }
      x++;
      v = y;
    }

  }

  ccm->sc_frag = sreLOG2(1./total_frags);
  ccm->e_fraglen = ((float) ebases)/total_frags;

  esl_stack_Destroy(oldstate);
  esl_stack_Destroy(newstate);

  return ccm;

ERROR:
  cm_Fail("Memory allocation error.\n");
  return NULL;
}

/* Function:  cm_consensus_Free()
   Author:    DLK, Thu May 14 2009
 */
void
cm_consensus_Free(CM_CONSENSUS *ccm)
{
  if (ccm->next    != NULL) free(ccm->next);
  if (ccm->sttype  != NULL) free(ccm->sttype);
  if (ccm->oesc[0] != NULL) free(ccm->oesc[0]);
  if (ccm->oesc    != NULL) free(ccm->oesc);

  esl_alphabet_Destroy(ccm->abc);

  free(ccm); ccm = NULL;

  return;
}

/* wordify()
 * Converts log probability score to a rounded signed 16-bit integer cost.
 * Both emissions and transitions for ViterbiFilter get this treatment.
 * No bias term needed, because we use signed words. 
 *   e.g. a score of +3.2, with scale 500.0, becomes +1600.
 */
int16_t
wordify(float scale_w, float sc)
{
  sc  = roundf(scale_w * sc);
  if      (sc >=  32767.0) return  32767;
  else if (sc <= -32768.0) return -32768;
  else return (int16_t) sc;
}

/* biased_byteify()
 * Converts original log-odds residue score to a rounded biased uchar cost.
 * e.g. a score of +3.2, with scale 3.0 and bias 12, becomes 2.
 *    3.2*3 = 9.6; rounded = 10; bias-10 = 2.
 * When used, we add the bias, then subtract this cost.
 * (A cost of +255 is our -infinity "prohibited event")
 */
uint8_t
biased_byteify(CM_CONSENSUS *ccm, float sc)
{
  uint8_t b;

  sc  = -1.0f * roundf(ccm->scale_b * sc);                          /* ugh. sc is now an integer cost represented in a float...
    */
  b   = (sc > 255 - ccm->bias_b) ? 255 : (uint8_t) sc + ccm->bias_b; /* and now we cast, saturate, and bias it to an unsigned char cost
... */
  return b;
}

/* unbiased_byteify()
 * Convert original transition score to a rounded uchar cost
 * e.g. a score of -2.1, with scale 3.0, becomes a cost of 6.
 * (A cost of +255 is our -infinity "prohibited event")
 */
uint8_t
unbiased_byteify(CM_CONSENSUS *ccm, float sc)
{
  uint8_t b;

  sc  = -1.0f * roundf(ccm->scale_b * sc);       /* ugh. sc is now an integer cost represented in a float...    */
  b   = (sc > 255.) ? 255 : (uint8_t) sc;       /* and now we cast and saturate it to an unsigned char cost... */
  return b;
}

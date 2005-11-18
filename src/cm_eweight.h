/************************************************************
 * @LICENSE@
 ************************************************************/

/* cm_eweight.h
 * based on:
 */
/* lsjfuncs.h
 * Declarations of external functions used in lsj_eweight.c
 * (Entropy-based sequence weighting)
 *
 * Steve Johnson
 * CVS $Id: lsjfuncs.h 944 2004-05-24 15:49:07Z eddy $
 */


#include "config.h"
#include "structs.h"
#include "squid.h"
#include "msa.h"


/*
extern float CM_Eweight(struct CM_t *cm,  struct Prior_t *pri, 
		     float numb_seqs, float targetent);
*/

extern double CM_Eweight(CM_t *cm,  Prior_t *pri, 
		     float numb_seqs, float targetent);

extern void ModelContent(float *ent1, float *ent2, int M);
extern void CMRescale(CM_t *hmm, float scale);

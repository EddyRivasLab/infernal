/************************************************************
 * @LICENSE@
 ************************************************************/

/* lsjfuncs.h
 * Declarations of external functions used in lsj_eweight.c
 * (Entropy-based sequence weighting)
 *
 * Steve Johnson
 * CVS $Id: lsjfuncs.h 944 2004-05-24 15:49:07Z eddy $
 */


#include "hmmer_config.h"
#include "hmmer_structs.h"
#include "squid.h"
#include "msa.h"

extern float Eweight(struct plan7_s *hmm,  struct p7prior_s *pri, 
		     float numb_seqs, float entwgt);
extern void ModelContent(float *ent1, float *ent2, int M);

/* p7_modelmaker.c
 * EPN, Tue Aug  5 15:32:34 2008
 * SVN $Id: cm_modelmaker.c 2327 2008-02-13 22:09:06Z nawrockie $
 *
 * Construct a p7 model from CM and it's CP9 HMM.
 *
 *****************************************************************
 * @LICENSE@
 *****************************************************************
 */


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

#include "easel.h"		
#include "esl_msa.h"		
#include "esl_stack.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "hmmer.h"

#include "funcs.h"
#include "structs.h"

/* Function: BuildP7HMM_MatchEmitsOnly()
 * Incept:   EPN, Tue Aug  5 15:33:00 2008
 * 
 * Purpose:  Create and fill a P7_HMM object from a CM and it's CP9 HMM.
 *           Copy only the match emissions of the CP9 HMM, the rest of 
 *           the p7 model parameters are irrelevant. 
 * 
 * Args:     cm        - the cm
 *           ret_p7    - RETURN: new p7 model 
 *           ret_gm    - RETURN: new p7 generic profile
 *           ret_om    - RETURN: new p7 optimized profile
 *           
 * Return:   eslOK   on success
 *
 * Throws:   eslEINCOMPAT on contract violation
 *           eslEMEM on memory error
 */
int
BuildP7HMM_MatchEmitsOnly(CM_t *cm, P7_HMM **ret_p7, P7_PROFILE **ret_gm, P7_OPROFILE **ret_om)
{
  int        status;
  P7_HMM     *hmm = NULL;        /* RETURN: new hmm */
  P7_PROFILE  *gm = NULL;        /* RETURN: new generic profile */
  P7_OPROFILE *om = NULL;        /* RETURN: new optimized profile */
  int        k;

  if(cm->cp9 == NULL)         return eslEINCOMPAT; 
  if(! (cm->flags & CMH_CP9)) return eslEINCOMPAT; 

  if ((hmm    = p7_hmm_Create(cm->clen, cm->abc)) == NULL)  return eslEMEM;
  if ((status = p7_hmm_Zero(hmm))                 != eslOK) return status;

  /* copy only match emissions */
  for (k = 1; k <= cm->clen; k++) esl_vec_FCopy(cm->cp9->mat[k], cm->abc->K, hmm->mat[k]);

  /* parameterize, hacked from hmmer/p7_prior.c::p7_ParameterEstimation() */
  /* match transitions */
  for (k = 1; k <= hmm->M; k++) esl_vec_FNorm(hmm->t[k],   3); 

  /* insert transitions */
  for (k = 1; k <= hmm->M; k++) esl_vec_FNorm(hmm->t[k]+3, 2); 

  /* delete transitions */
  for (k = 1; k < hmm->M; k++) esl_vec_FNorm(hmm->t[k]+5, 2); 
  /* For k=0, which is unused; convention sets TMM=1.0, TMD=0.0
   * For k=M, TMM = 1.0 (to the E state) and TMD=0.0 (no next D; must go to E).
   */
  hmm->t[0][p7H_DM] = hmm->t[hmm->M][p7H_DM] = 1.0;
  hmm->t[0][p7H_DD] = hmm->t[hmm->M][p7H_DD] = 0.0;

  /* insert emissions */
  for (k = 0; k <= hmm->M; k++) esl_vec_FNorm(hmm->ins[k], hmm->abc->K); /* normalize inserts (0.25 each) */

  p7_hmm_SetName(hmm, cm->name);
  p7_hmm_SetAccession(hmm, cm->acc);
  p7_hmm_SetDescription(hmm, cm->desc);
  p7_hmm_SetCtime(hmm);
  ESL_ALLOC(hmm->comlog, sizeof(char)* (strlen(cm->comlog->bcom)+1));
  *(hmm->comlog) = '\0'; /* need this to make strcat work */
  strcat(hmm->comlog, cm->comlog->bcom);

  hmm->eff_nseq = cm->eff_nseq;
  hmm->nseq     = cm->nseq;
  hmm->checksum = 0;

  /* make the profiles */
  gm = p7_profile_Create (hmm->M, hmm->abc);
  om = p7_oprofile_Create(hmm->M, hmm->abc);

  *ret_p7 = hmm;
  *ret_gm = gm;
  *ret_om = om;

  return eslOK;

 ERROR: 
  if(hmm != NULL) p7_hmm_Destroy(hmm);
  return status;
}

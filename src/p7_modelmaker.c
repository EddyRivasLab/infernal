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
#include "esl_config.h"
#include "p7_config.h"
#include "config.h"

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

/* EPN, Mon Aug 25 09:00:10 2008
 * Had difficulty compiling infernal with hmmer as a subdir and with
 * impl_sse. Because I currently don't need OPROFILEs I #if 0ed out the code
 * here and made a new version of this func which doesn't return an omx (BELOW).
 */
#if 0
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
  if(cm->comlog != NULL && cm->comlog->bcom != NULL) { 
    ESL_ALLOC(hmm->comlog, sizeof(char)* (strlen(cm->comlog->bcom)+1));
    *(hmm->comlog) = '\0'; /* need this to make strcat work */
    strcat(hmm->comlog, cm->comlog->bcom);
  }
  else hmm->comlog = NULL;

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
#endif

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
 *           
 * Return:   eslOK   on success
 *
 * Throws:   eslEINCOMPAT on contract violation
 *           eslEMEM on memory error
 */
int
BuildP7HMM_MatchEmitsOnly(CM_t *cm, P7_HMM **ret_p7, P7_PROFILE **ret_gm)
{
  int        status;
  P7_HMM     *hmm = NULL;        /* RETURN: new hmm */
  P7_PROFILE  *gm = NULL;        /* RETURN: new generic profile */
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
  if(cm->comlog != NULL && cm->comlog->bcom != NULL) { 
    ESL_ALLOC(hmm->comlog, sizeof(char)* (strlen(cm->comlog->bcom)+1));
    *(hmm->comlog) = '\0'; /* need this to make strcat work */
    strcat(hmm->comlog, cm->comlog->bcom);
  }
  else hmm->comlog = NULL;

  hmm->eff_nseq = cm->eff_nseq;
  hmm->nseq     = cm->nseq;
  hmm->checksum = 0;

  /* make the profiles */
  gm = p7_profile_Create (hmm->M, hmm->abc);

  *ret_p7 = hmm;
  *ret_gm = gm;


  return eslOK;

 ERROR: 
  if(hmm != NULL) p7_hmm_Destroy(hmm);
  return status;
}

/* Function: CP9_to_P7()
 * Incept:   EPN, Fri Sep 24 13:46:37 2010
 * 
 * Purpose:  Create and fill a P7_HMM object from a CM and it's CP9 HMM.
 * 
 * Args:     cm        - the cm, must have a cp9 model in it.
 *           ret_p7    - RETURN: new p7 model 
 *           
 * Return:   eslOK   on success
 *
 * Throws:   eslEINCOMPAT on contract violation
 *           eslEMEM on memory error
 */
int
CP9_to_P7(CM_t *cm, P7_HMM **ret_p7)
{
  int        status;
  P7_HMM     *hmm = NULL;        /* RETURN: new hmm */
  int        k;

  if(cm->cp9 == NULL)         return eslEINCOMPAT; 
  if(! (cm->flags & CMH_CP9)) return eslEINCOMPAT; 

  if ((hmm    = p7_hmm_Create(cm->clen, cm->abc)) == NULL)  return eslEMEM;
  if ((status = p7_hmm_Zero(hmm))                 != eslOK) return status;

  /* copy transitions */
  for (k = 0; k <= hmm->M; k++) { 
    hmm->t[k][p7H_MM] = cm->cp9->t[k][CTMM];
    hmm->t[k][p7H_MI] = cm->cp9->t[k][CTMI];
    hmm->t[k][p7H_MD] = cm->cp9->t[k][CTMD];
    hmm->t[k][p7H_IM] = cm->cp9->t[k][CTIM];
    hmm->t[k][p7H_II] = cm->cp9->t[k][CTII];
    hmm->t[k][p7H_DM] = cm->cp9->t[k][CTDM];
    hmm->t[k][p7H_DD] = cm->cp9->t[k][CTDD];
    /* note: the cp9 CTDI and CTID transitions do not exist the p7 model */
  }
  /* normalize match transitions */
  for (k = 1; k <= hmm->M; k++) esl_vec_FNorm(hmm->t[k],  3); 
  /* normalize insert transitions */
  for (k = 0; k < hmm->M; k++) esl_vec_FNorm(hmm->t[k]+3, 2); 
  /* normalize delete transitions */
  for (k = 1; k < hmm->M; k++) esl_vec_FNorm(hmm->t[k]+5, 2); 

  /* enforce HMMER conventions */
  hmm->t[hmm->M][p7H_MD] = 0.0;
  esl_vec_FNorm(hmm->t[hmm->M], 3);
  hmm->t[0][p7H_DM] = hmm->t[hmm->M][p7H_DM] = 1.0;
  hmm->t[0][p7H_DD] = hmm->t[hmm->M][p7H_DD] = 0.0;

  /* enforce INFERNAL CP9 convention, the 0'th node's MM transition is really begin[0] */
  hmm->t[0][p7H_MM] = cm->cp9->begin[1];
  esl_vec_FNorm(hmm->t[0], 3);

  /* match emissions: copy, then normalize (should be unnec actually) */
  for (k = 1; k <= cm->clen; k++) esl_vec_FCopy(cm->cp9->mat[k], cm->abc->K, hmm->mat[k]);
  for (k = 1; k <= cm->clen; k++) esl_vec_FNorm(cm->cp9->mat[k], cm->abc->K);
  /* special case */
  esl_vec_FSet(hmm->mat[0], hmm->abc->K, 0.);
  hmm->mat[0][0] = 1.0;

  /* insert emissions: copy, then normalize (should be unnec actually) */
  for (k = 0; k <= cm->clen; k++) esl_vec_FCopy(cm->cp9->ins[k], cm->abc->K, hmm->ins[k]);
  for (k = 0; k <= cm->clen; k++) esl_vec_FNorm(cm->cp9->ins[k], cm->abc->K);

  /* copy the all important max length parameter */
  hmm->max_length = cm->W;

  p7_hmm_SetName       (hmm, cm->name);
  p7_hmm_SetAccession  (hmm, cm->acc);
  p7_hmm_SetDescription(hmm, cm->desc);
  p7_hmm_SetCtime      (hmm);
  if(cm->comlog != NULL && cm->comlog->bcom != NULL) { 
    ESL_ALLOC(hmm->comlog, sizeof(char)* (strlen(cm->comlog->bcom)+1));
    *(hmm->comlog) = '\0'; /* need this to make strcat work */
    strcat(hmm->comlog, cm->comlog->bcom);
  }
  else hmm->comlog = NULL;

  hmm->eff_nseq = cm->eff_nseq;
  hmm->nseq     = cm->nseq;
  hmm->checksum = 0;

  /* set the model composition */
  if ((status = p7_hmm_SetComposition(hmm)) != eslOK) goto ERROR;

  /* finally, calibrate the model */
  if((status = p7_Calibrate(hmm, NULL, NULL, NULL, NULL, NULL)) != eslOK) goto ERROR; 

  *ret_p7 = hmm;

  return eslOK;

 ERROR: 
  if(hmm != NULL) p7_hmm_Destroy(hmm);
  return status;
}


/* Function: dump_p7()
 * Incept:   EPN, Fri Sep 24 14:22:49 2010
 * 
 * Purpose:  Dump parameters of a p7 HMM to a file.
 * 
 * Args:     hmm       - the p7 HMM
 *           fp        - the file to print to
 *           
 * Return:   eslOK   on success
 *
 * Throws:   eslEINCOMPAT on contract violation
 *           eslEMEM on memory error
 */
int
dump_p7(P7_HMM *hmm, FILE *fp)
{
  return p7_hmmfile_WriteASCII(fp, -1, hmm);
}

/************************************************************
 * @LICENSE@
 ************************************************************/
/* cp9.c based on HMMER 2.x's plan7.c
 * EPN 02.27.06
 * 
 * Support for CM-Plan 9 HMM data structure, CP9_t.
 * 
 * All of the CP9 code is based on analogous HMMER 2.x Plan 7 HMM
 * code.  
 * 
 */

#include "esl_config.h"
#include "p7_config.h"
#include "config.h"

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>

#include "easel.h"
#include "esl_msa.h"
#include "esl_vectorops.h"

#include "funcs.h"
#include "structs.h"

/* Functions: AllocCPlan9(), AllocCPlan9Shell(), AllocCPlan9Body(), FreeCPlan9()
 * 
 * Purpose:   Allocate or free a CPlan9 HMM structure.
 *            Can either allocate all at once (AllocCPlan9()) or
 *            in two steps (AllocCPlan9Shell(), AllocCPlan9Body()).
 *            The two step method is used in CP9_hmmio.c where we start
 *            parsing the header of an HMM file but don't 
 *            see the size of the model 'til partway thru the header.
 */
CP9_t *
AllocCPlan9(int M, const ESL_ALPHABET *abc) 
{
  CP9_t *hmm;

  hmm = AllocCPlan9Shell();
  AllocCPlan9Body(hmm, M, abc);
  return hmm;
}  
CP9_t *
AllocCPlan9Shell(void) 
{
  int    status;
  CP9_t *hmm;

  ESL_ALLOC(hmm, sizeof(CP9_t));
  hmm->abc = NULL;

  hmm->M = 0;

  hmm->t      = NULL;
  hmm->mat    = NULL;
  hmm->ins    = NULL;
  
  hmm->tsc     = hmm->msc     = hmm->isc     = NULL;
  hmm->tsc_mem = hmm->msc_mem = hmm->msc_mem = NULL;

  hmm->begin  = NULL;
  hmm->end    = NULL;

  hmm->bsc = hmm->bsc_mem = NULL;
  hmm->esc = hmm->esc_mem = NULL;

  hmm->otsc = NULL;

  hmm->has_el      = NULL;
  hmm->el_from_ct  = NULL;
  hmm->el_from_idx = NULL;
  hmm->el_from_cmnd= NULL;

  hmm->flags = 0;
  return hmm;

 ERROR:
  cm_Fail("Memory allocation error.\n");
  return NULL; /* never reached */
}  

void
AllocCPlan9Body(CP9_t *hmm, int M, const ESL_ALPHABET *abc) 
{
  int status;
  int k, x;

  hmm->abc = abc;

  hmm->M = M;

  ESL_ALLOC(hmm->t,   (M+1) *           sizeof(float *));
  ESL_ALLOC(hmm->mat, (M+1) *           sizeof(float *));
  ESL_ALLOC(hmm->ins, (M+1) *           sizeof(float *));
  ESL_ALLOC(hmm->t[0],(cp9_NTRANS*(M+1))     *  sizeof(float));
  ESL_ALLOC(hmm->mat[0],(abc->K*(M+1)) * sizeof(float));
  ESL_ALLOC(hmm->ins[0],(abc->K*(M+1)) * sizeof(float));

  ESL_ALLOC(hmm->tsc, cp9_NTRANS *       sizeof(int *));
  ESL_ALLOC(hmm->msc, hmm->abc->Kp   *   sizeof(int *));
  ESL_ALLOC(hmm->isc, hmm->abc->Kp   *   sizeof(int *)); 
  ESL_ALLOC(hmm->tsc_mem,(cp9_NTRANS*(M+1))     *       sizeof(int));
  ESL_ALLOC(hmm->msc_mem,(hmm->abc->Kp*(M+1)) * sizeof(int));
  ESL_ALLOC(hmm->isc_mem,(hmm->abc->Kp*(M+1)) * sizeof(int));

  hmm->tsc[0] = hmm->tsc_mem;
  hmm->msc[0] = hmm->msc_mem;
  hmm->isc[0] = hmm->isc_mem;

  /* transition scores reordered */
  ESL_ALLOC(hmm->otsc, sizeof(int)   * (M+1)  * cp9O_NTRANS);

  /* note allocation strategy for important 2D arrays -- trying
   * to keep locality as much as possible, cache efficiency etc.
   */
  for (k = 1; k <= M; k++) {
    hmm->mat[k] = hmm->mat[0] + k * abc->K;
    hmm->ins[k] = hmm->ins[0] + k * abc->K;
    hmm->t[k]   = hmm->t[0]   + k * cp9_NTRANS;
  }
  for (x = 1; x < hmm->abc->Kp; x++) {
    hmm->msc[x] = hmm->msc[0] + x * (M+1);
    hmm->isc[x] = hmm->isc[0] + x * (M+1);
  }
  for (x = 0; x < cp9_NTRANS; x++)
    hmm->tsc[x] = hmm->tsc[0] + x * (M+1);

  /* tsc[x][0] is used as a boundary condition sometimes [Viterbi()],
   * so set to -inf always.
   */
  for (x = 0; x < cp9_NTRANS; x++)
    hmm->tsc[x][0] = -INFTY;

  ESL_ALLOC(hmm->begin, (M+1) * sizeof(float));
  ESL_ALLOC(hmm->end,   (M+1) * sizeof(float));

  ESL_ALLOC(hmm->bsc_mem, (M+1) * sizeof(int));
  ESL_ALLOC(hmm->esc_mem, (M+1) * sizeof(int));

  ESL_ALLOC(hmm->null, (abc->K) * sizeof(float));

  hmm->bsc = hmm->bsc_mem;
  hmm->esc = hmm->esc_mem;

  /* end[0], begin[0], esc[0] and bsc[0] are never
   * used, set them to 0. and -INFTY */
  hmm->end[0] = hmm->begin[0] = -INFTY;
  hmm->esc[0] = hmm->bsc[0] = -INFTY;
  
  ESL_ALLOC(hmm->has_el,     (M+1) * sizeof(int));
  ESL_ALLOC(hmm->el_from_ct, (M+2) * sizeof(int));
  ESL_ALLOC(hmm->el_from_idx,(M+2) * sizeof(int *));
  ESL_ALLOC(hmm->el_from_cmnd,(M+2) * sizeof(int *));

  return;
 ERROR:
  cm_Fail("Memory allocation error.");
}  


void
FreeCPlan9(CP9_t *hmm)
{
  int k;
  if (hmm->null       != NULL) free(hmm->null);
  if (hmm->bsc_mem    != NULL) free(hmm->bsc_mem);
  if (hmm->begin      != NULL) free(hmm->begin);
  if (hmm->esc_mem    != NULL) free(hmm->esc_mem);
  if (hmm->end        != NULL) free(hmm->end);
  if (hmm->msc_mem    != NULL) free(hmm->msc_mem);
  if (hmm->isc_mem    != NULL) free(hmm->isc_mem);
  if (hmm->tsc_mem    != NULL) free(hmm->tsc_mem);
  if (hmm->mat        != NULL) free(hmm->mat[0]);
  if (hmm->ins        != NULL) free(hmm->ins[0]);
  if (hmm->t          != NULL) free(hmm->t[0]);
  if (hmm->msc        != NULL) free(hmm->msc);
  if (hmm->isc        != NULL) free(hmm->isc);
  if (hmm->tsc        != NULL) free(hmm->tsc);
  if (hmm->otsc       != NULL) free(hmm->otsc);
  if (hmm->mat        != NULL) free(hmm->mat);
  if (hmm->ins        != NULL) free(hmm->ins);
  if (hmm->t          != NULL) free(hmm->t);
  if (hmm->has_el     != NULL) free(hmm->has_el);
  if (hmm->el_from_ct != NULL) free(hmm->el_from_ct);
  if(hmm->el_from_idx != NULL)
    {
      for(k = 0; k <= hmm->M+1; k++)
	if(hmm->el_from_idx[k] != NULL)
	  free(hmm->el_from_idx[k]);
      free(hmm->el_from_idx);
    }
  if(hmm->el_from_cmnd != NULL)
    {
      for(k = 0; k <= hmm->M+1; k++)
	if(hmm->el_from_cmnd[k] != NULL)
	  free(hmm->el_from_cmnd[k]);
      free(hmm->el_from_cmnd);
    }

  free(hmm);
}

/* Function: ZeroCPlan9()
 * 
 * Purpose:  Zeros the counts/probabilities fields in a model.  
 *           Leaves null model untouched. 
 */
void
ZeroCPlan9(CP9_t *hmm)
{
  int k;
  esl_vec_FSet(hmm->ins[0], hmm->abc->K, 0.);
  esl_vec_FSet(hmm->t[0], cp9_NTRANS, 0.);
  for (k = 1; k <= hmm->M; k++)
    {
      esl_vec_FSet(hmm->t[k], cp9_NTRANS, 0.);
      esl_vec_FSet(hmm->mat[k], hmm->abc->K, 0.);
      esl_vec_FSet(hmm->ins[k], hmm->abc->K, 0.);
    }
  esl_vec_FSet(hmm->begin+1, hmm->M, 0.);
  esl_vec_FSet(hmm->end+1, hmm->M, 0.);

  /* initialize the el_* data structures, these
   * depend on the CM guide tree and will be set
   * when the CP9 is constructed from the CM.
   */
  for (k = 0; k <= (hmm->M); k++)
    {
      hmm->has_el[k]      = FALSE;         
      hmm->el_from_ct[k]  = 0;
      hmm->el_from_idx[k] = NULL; 
      hmm->el_from_cmnd[k] = NULL; 
    }
  /* special case hmm->M+1 corresponds to the E state here */
  hmm->el_from_ct[(hmm->M+1)]  = 0;
  hmm->el_from_idx[(hmm->M+1)] = NULL; 
  hmm->el_from_cmnd[(hmm->M+1)] = NULL; 

  hmm->flags &= ~CPLAN9_HASBITS;	/* invalidates scores */
  hmm->flags &= ~CPLAN9_HASPROB;	/* invalidates probabilities */
  hmm->el_self = 0.; /* EL self transition probability */
}


/* Function: CPlan9SetNullModel()
 * 
 * Purpose:  Set the null model section of an HMM.
 *           Convenience function.
 *
 *            Assumes null* is allocated to hmm->abc->K
 */
void
CPlan9SetNullModel(CP9_t *hmm, float *null, float p1)
{
  int x;
  for (x = 0; x < hmm->abc->K; x++)
    hmm->null[x] = null[x];
  hmm->p1 = p1;
}

/*
 * Function: DuplicateCP9
 * Date:     EPN, Thu Jun 28 13:37:22 2007
 * Purpose:  Given a template cm 'src_cm' copy it's CP9 
 *           to the cm 'dest_cm'. dest_cm->cp9 and dest_cm->cp9map
 *           are alloc'ed here.
 *
 * Args:
 *           src_cm         the source CM, must have valid cp9
 *           dest_cm        the destination CM we're copying src_cm->cp9
 *                          to
 */
void
DuplicateCP9(CM_t *src_cm, CM_t *dest_cm)
{
  int       k,x;	          /* counter over nodes */

  /* Contract checks */
  if(!(src_cm->flags & CMH_CP9))
    cm_Fail("ERROR in DuplicateCP9() src_cm CMH_CP9 flag down.\n");

  CPlan9Renormalize(src_cm->cp9);
  CP9Logoddsify(src_cm->cp9);
  
  /* We can fill the map before we copy the CP9 */
  /* Allocate and initialize the cp9map */
  dest_cm->cp9map = AllocCP9Map(dest_cm);
  /* Map the CM states to CP9 states and nodes and vice versa */
  CP9_map_cm2hmm(dest_cm, dest_cm->cp9map, 0);

  /* Create the new model and copy everything over */
  dest_cm->cp9 = AllocCPlan9(dest_cm->cp9map->hmm_M, src_cm->abc);
  ZeroCPlan9(dest_cm->cp9);
  CPlan9SetNullModel(dest_cm->cp9, dest_cm->null, 1.0); /* set p1 = 1.0 which corresponds to the CM */
  CPlan9InitEL(dest_cm, dest_cm->cp9); /* set up hmm->el_from_ct and hmm->el_from_idx data, which
					* explains how the EL states are connected in the HMM. */

  /* Copy the transitions and emission probs and scores */
  for(k = 0; k <= dest_cm->cp9->M; k++)
    {
      dest_cm->cp9->t[k][CTMM] = src_cm->cp9->t[k][CTMM];
      dest_cm->cp9->t[k][CTMI] = src_cm->cp9->t[k][CTMI];
      dest_cm->cp9->t[k][CTMD] = src_cm->cp9->t[k][CTMD];
      dest_cm->cp9->t[k][CTMEL] = src_cm->cp9->t[k][CTMEL];


      dest_cm->cp9->t[k][CTIM] = src_cm->cp9->t[k][CTIM];
      dest_cm->cp9->t[k][CTII] = src_cm->cp9->t[k][CTII];
      dest_cm->cp9->t[k][CTID] = src_cm->cp9->t[k][CTID];

      dest_cm->cp9->t[k][CTDM] = src_cm->cp9->t[k][CTDM];
      dest_cm->cp9->t[k][CTDI] = src_cm->cp9->t[k][CTDI];
      dest_cm->cp9->t[k][CTDD] = src_cm->cp9->t[k][CTDD];

      dest_cm->cp9->begin[k] = src_cm->cp9->begin[k];
      dest_cm->cp9->end[k] = src_cm->cp9->end[k];

      for(x = 0; x < src_cm->cp9->abc->K; x++)
	{
	  dest_cm->cp9->mat[k][x] = src_cm->cp9->mat[k][x];
	  dest_cm->cp9->ins[k][x] = src_cm->cp9->ins[k][x];
	}
    }

  dest_cm->cp9->p1        = src_cm->cp9->p1;
  dest_cm->cp9->el_self   = src_cm->cp9->el_self;
  dest_cm->cp9->el_selfsc = src_cm->cp9->el_selfsc;

  CPlan9Renormalize(dest_cm->cp9);/* shouldn't be nec */
  CP9Logoddsify(dest_cm->cp9); /* fill in all the integer log odds scores:
				* msc, isc, bsc, esc, tsc, the *sc_mem
				* pointers were set up in AllocCPlan9() */
    
  /*
    FILE *fp;
    fp = fopen("destcp9" ,"w");
    debug_print_cp9_params(fp, dest_cm->cp9, TRUE);
    fclose(fp);
    fp = fopen("srccp9" ,"w");
    debug_print_cp9_params(fp, src_cm->cp9, TRUE);
    fclose(fp);
  */

  /* finally create the CP9Bands */
  dest_cm->cp9b = AllocCP9Bands(dest_cm, dest_cm->cp9);

  dest_cm->cp9->flags = src_cm->cp9->flags;
}

/* Function: cp9_GetNCalcsPerResidue()
 * Date:     EPN, Thu Jan 17 06:12:37 2008
 * 
 * Returns: eslOK on success, eslEINCOMPAT on contract violation.
 *          <ret_cp9_ncalcs_per_res> set as millions of DP calculations 
 *          per residue for the CP9 HMM.
 */
int
cp9_GetNCalcsPerResidue(CP9_t *cp9, char *errbuf, float *ret_cp9_ncalcs_per_res)
{
  int cp9_ntrans;
  float cp9_ncalcs_per_res;
  
  if(cp9 == NULL)                    ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_GetNCalcsPerRes(), cp9 == NULL.");
  if(ret_cp9_ncalcs_per_res == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_GetNCalcsPerRes(), ret_cp9_ncalcs_per_res == NULL.");

  /* determine millions of CP9 DP calcs per residue */
  cp9_ntrans = NHMMSTATETYPES * NHMMSTATETYPES; /* 3*3 = 9 transitions in global mode */
  if(cp9->flags & CPLAN9_LOCAL_BEGIN) cp9_ntrans++; 
  if(cp9->flags & CPLAN9_LOCAL_END)   cp9_ntrans++; 
  if(cp9->flags & CPLAN9_EL)          cp9_ntrans++; 
  cp9_ncalcs_per_res = (cp9_ntrans * cp9->M) / 1000000.; /* convert to millions of calcs per residue */

  *ret_cp9_ncalcs_per_res = cp9_ncalcs_per_res;
  return eslOK;
}

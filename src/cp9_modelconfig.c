/************************************************************
 * @LICENSE@
 ************************************************************/
/* cp9_modelconfig.c
 * EPN, Wed Dec  5 12:56:42 2007
 *
 * Note: all of these functions originated in cp9.c [EPN 02.27.06]
 * 
 * Configuring a CP9 HMM into different global or local modes.
 * Included in this file are functions for configuring HMMs that were
 * built for 'sub CMs'.
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

/* Function: CP9Logoddsify()
 * 
 * Purpose:  Take an HMM with valid probabilities, and
 *           fill in the integer log-odds score section of the model.
 *           
 *    Notes on log-odds scores (simplified from plan7.c):
 *         type of parameter       probability        score
 *         -----------------       -----------        ------
 *         any emission             p_x           log_2 p_x/null_x
 *         any transition           t_x           log_2 t_x
 *             
 * Args:      hmm          - the hmm to calculate scores in.
 *                  
 * Return:    (void)
 *            hmm scores are filled in.
 */  
void
CP9Logoddsify(CP9_t *hmm)
{
  int k;			/* counter for model position */
  int x;			/* counter for symbols        */
  int *sc;
  int status;

  if (hmm->flags & CPLAN9_HASBITS) return;

  ESL_ALLOC(sc, hmm->abc->Kp * sizeof(int));

  /* Symbol emission scores
   */

  sc[hmm->abc->K]     = -INFTY; /* gap character */
  sc[hmm->abc->Kp-1]  = -INFTY; /* missing data character */
  sc[hmm->abc->Kp-2]  = -INFTY; /* non-residue data character */

  /* Insert emission scores, relies on sc[K, Kp-1] initialization to -inf above */
  for (k = 0; k <= hmm->M; k++) {
    for (x = 0; x < hmm->abc->K; x++) 
      sc[x] = Prob2Score(hmm->ins[k][x], hmm->null[x]);
    esl_abc_IExpectScVec(hmm->abc, sc, hmm->null); 
    for (x = 0; x < hmm->abc->Kp; x++) {
      hmm->isc[x][k] = sc[x];
    }
  }

  /* Match emission scores, relies on sc[K, Kp-1] initialization to -inf above */
  for (k = 1; k <= hmm->M; k++) {
    for (x = 0; x < hmm->abc->K; x++) 
      sc[x] = Prob2Score(hmm->mat[k][x], hmm->null[x]);
    esl_abc_IExpectScVec(hmm->abc, sc, hmm->null); 
    for (x = 0; x < hmm->abc->Kp; x++) {
      hmm->msc[x][k] = sc[x];
    }
  }
  
  for (k = 0; k <= hmm->M; k++)
    {
      hmm->tsc[CTMM][k] = Prob2Score(hmm->t[k][CTMM], 1.0);
      hmm->tsc[CTMI][k] = Prob2Score(hmm->t[k][CTMI], 1.0);
      hmm->tsc[CTMD][k] = Prob2Score(hmm->t[k][CTMD], 1.0);
      hmm->tsc[CTMEL][k] = Prob2Score(hmm->t[k][CTMEL], 1.0);
      hmm->tsc[CTIM][k] = Prob2Score(hmm->t[k][CTIM], 1.0);
      hmm->tsc[CTII][k] = Prob2Score(hmm->t[k][CTII], 1.0);
      hmm->tsc[CTID][k] = Prob2Score(hmm->t[k][CTID], 1.0);
      if(k != 0)
	{
	  hmm->tsc[CTDM][k] = Prob2Score(hmm->t[k][CTDM], 1.0);
	  hmm->tsc[CTDI][k] = Prob2Score(hmm->t[k][CTDI], 1.0);
	  hmm->tsc[CTDD][k] = Prob2Score(hmm->t[k][CTDD], 1.0);
	}
      else
	{
	  hmm->tsc[CTDM][k] = -INFTY;
	  hmm->tsc[CTDD][k] = -INFTY; /*D_0 doesn't exist*/
	  hmm->tsc[CTDI][k] = -INFTY;
	}
      if(k != 0)
	{
	  hmm->bsc[k]   = Prob2Score(hmm->begin[k], 1.0);
	  //if(hmm->flags & CPLAN9_LOCAL_END) hmm->esc[k]   = 0;
	  //else hmm->esc[k]   = -INFTY;
	  hmm->esc[k] = Prob2Score(hmm->end[k], 1.0);
	}
    }
  hmm->el_selfsc = Prob2Score(hmm->el_self, 1.0);

  /* Finally, fill the efficiently reordered transition scores for this HMM. */
  for (k = 0 ; k <= hmm->M; k++) {
    int *otsc_k = hmm->otsc + k*cp9O_NTRANS;
    otsc_k[cp9O_MM] = hmm->tsc[CTMM][k];
    otsc_k[cp9O_MI] = hmm->tsc[CTMI][k];
    otsc_k[cp9O_MD] = hmm->tsc[CTMD][k];
    otsc_k[cp9O_IM] = hmm->tsc[CTIM][k];
    otsc_k[cp9O_II] = hmm->tsc[CTII][k];
    otsc_k[cp9O_DM] = hmm->tsc[CTDM][k];
    otsc_k[cp9O_DD] = hmm->tsc[CTDD][k];
    otsc_k[cp9O_ID] = hmm->tsc[CTID][k];
    otsc_k[cp9O_DI] = hmm->tsc[CTDI][k];
    otsc_k[cp9O_BM] = hmm->bsc[k];
    otsc_k[cp9O_MEL]= hmm->tsc[CTMEL][k];
    otsc_k[cp9O_ME] = hmm->esc[k];
  }

  hmm->flags |= CPLAN9_HASBITS;	/* raise the log-odds ready flag */

  free(sc);

  return;

 ERROR:
  cm_Fail("Memory allocation error.\n");
  return; /* never reached */
}

/* Function: CPlan9Renormalize()
 * 
 * Purpose:  Take an HMM in counts form, and renormalize
 *           all of its probability vectors. Also enforces
 *           CM Plan9 restrictions on nonexistent transitions.
 *           
 * Args:     hmm - the model to renormalize.
 *                 
 * Return:   (void)
 *           hmm is changed.
 */                          
void
CPlan9Renormalize(CP9_t *hmm)
{
  int   k;			/* counter for model position */
  float d;			/* denominator */

				/* match emissions */
  esl_vec_FSet(hmm->mat[0], hmm->abc->K, 0.);   /*M_0 is B state, non-emitter*/
  for (k = 1; k <= hmm->M; k++) 
    esl_vec_FNorm(hmm->mat[k], hmm->abc->K);
				/* insert emissions */
  for (k = 0; k <= hmm->M; k++)
    esl_vec_FNorm(hmm->ins[k], hmm->abc->K);

				/* begin transitions */
  d = esl_vec_FSum(hmm->begin+1, hmm->M) + hmm->t[0][CTMI] + hmm->t[0][CTMD] + hmm->t[0][CTMEL]; 
  /* hmm->t[0][CTMEL] should always be 0., can't local end from the M_0 == B state */
  esl_vec_FScale(hmm->begin+1, hmm->M, 1./d);
  hmm->t[0][CTMI] /= d;
  hmm->t[0][CTMD] /= d;
  hmm->t[0][CTMEL] /= d;

  esl_vec_FNorm(hmm->t[0] + cp9_TRANS_INSERT_OFFSET, cp9_TRANS_NINSERT);	        /* transitions out of insert for node 0 (state N)*/
  esl_vec_FSet (hmm->t[0] + cp9_TRANS_DELETE_OFFSET, cp9_TRANS_NDELETE, 0.);    
				/* main model transitions */
  for (k = 1; k <= hmm->M; k++) /* safe for node M too, hmm->t[hmm->M][CTMM] should be 0.*/
    {
      d = esl_vec_FSum(hmm->t[k], cp9_TRANS_NMATCH) + hmm->end[k]; 
      esl_vec_FScale(hmm->t[k], cp9_TRANS_NMATCH, 1./d);
      hmm->end[k] /= d;

      esl_vec_FNorm(hmm->t[k] + cp9_TRANS_INSERT_OFFSET, cp9_TRANS_NINSERT);	/* insert */
      esl_vec_FNorm(hmm->t[k] + cp9_TRANS_DELETE_OFFSET, cp9_TRANS_NDELETE);	/* delete */
    }
                                 /* null model emissions */
  esl_vec_FNorm(hmm->null, hmm->abc->K);

  hmm->flags &= ~CPLAN9_HASBITS;	/* clear the log-odds ready flag */
  hmm->flags |= CPLAN9_HASPROB;	/* set the probabilities OK flag */
}

/* Function: CPlan9SWConfig()
 * EPN 05.30.06
 * based on SRE's Plan7SWConfig() from HMMER's plan7.c
 * 
 * Purpose:  Set the alignment independent parameters of
 *           a CM Plan 9 model to hmmsw (Smith/Waterman) configuration.
 *           
 * Notes:    The desideratum for begin/end probs is that all fragments ij
 *           (starting at match i, ending at match j) are
 *           equiprobable -- there is no information in the choice of
 *           entry/exit. There are M(M+1)/2 possible choices of ij, so
 *           each must get a probability of 2/M(M+1). This prob is the
 *           product of a begin, an end, and all the not-end probs in
 *           the path between i,j. 
 *            
 *           Thus: entry/exit is asymmetric because of the left/right
 *           nature of the HMM/profile. Entry probability is distributed
 *           simply by assigning p_x = pentry / (M-1) to M-1 
 *           internal match states. However, the same approach doesn't
 *           lead to a flat distribution over exit points. Exit p's
 *           must be corrected for the probability of a previous exit
 *           from the model. Requiring a flat distribution over exit
 *           points leads to an easily solved piece of algebra, giving:
 *                      p_1 = pexit / (M-1)
 *                      p_x = p_1 / (1 - (x-1) p_1)
 *
 *           Modified EPN, Thu Feb  7 15:54:16 2008, as follows:
 *           To better match a locally configured CM, if <do_match_local_cm>
 *           we disallow insertions before the first (emitting) match state, 
 *           (from I_0), and after the final (emitting) match state,
 *           (from I_M). I_0 maps to ROOT_IL and I_M maps to ROOT_IR
 *           which can never be entered in a locally configured CM
 *           (b/c the ROOT_S state MUST jump into a local begin state, which
 *            are always match states>). Also we disallow a M_0->D_1 transition
 *           because these would be impossible in a locally configured CM.
 *
 *           <do_match_local_cm> is usually TRUE, unless we're configuring
 *           the CP9 specifically for eventual sub CM alignment, where
 *           the goal is simply find the most likely start/end point
 *           of the alignment with this CP9 (in that case we want
 *           I_0 and I_M reachable).
 *           
 * Args:     hmm    - the CM Plan 9 model w/ data-dep prob's valid
 *           pentry - probability of an internal entry somewhere;
 *                    will be evenly distributed over M-1 match states
 *           pexit  - probability of an internal exit somewhere; 
 *                    will be distributed over M-1 match states.
 *           do_match_local_cm - TRUE to make I_0, D_1 and I_M unreachable
 *                    to better match a locally configured CM.
 *           first_cm_ndtype - only used if do_match_local_cm is TRUE
 *                             if it's MATL or MATP then D_1 should be unreachable (it is in the CM)
 *                             if it's MATR or MATP then D_M should be unreachable (it is in the CM)
 *                    
 * Return:   (void)
 *           HMM probabilities are modified.
 */
void
CPlan9SWConfig(CP9_t *hmm, float pentry, float pexit, int do_match_local_cm, int first_cm_ndtype)
{
  float basep;			/* p1 for exits: the base p */
  int   k;			/* counter over states      */
  float d;
  
  /* No special (*x* states in Plan 7) states in CM Plan 9 */

  /*for (k = 1; k <= hmm->M; k++) printf("before anything: end[%d]: %f\n", k, hmm->end[k]);*/
  /* Configure entry.
   */
  if(do_match_local_cm) { 
    hmm->t[0][CTMI] = 0.;
    hmm->t[0][CTMM] = 0.;  /* already was 0.0, transition from M_0 to M_1 is begin[1] */
    hmm->t[0][CTMEL] = 0.; /* already was 0.0, can never do a local end from M_0 */
    if((first_cm_ndtype == MATL_nd) || (first_cm_ndtype == MATP_nd)) { /* CM can't possibly reach the CM delete state that maps to D_1, make D_1 unreachable too */
      hmm->t[0][CTMD] = 0.;
    }

    hmm->t[hmm->M][CTMI] = 0.;
    hmm->t[hmm->M][CTDI] = 0.;
    if((first_cm_ndtype == MATR_nd) || (first_cm_ndtype == MATP_nd)) { /* CM can't possibly reach the CM delete state that maps to D_M, make D_M unreachable too */
      hmm->t[hmm->M][CTMD] = 0.;
    }

    /* renormalize transitions out of M_M */
    d = esl_vec_FSum(hmm->t[hmm->M], cp9_TRANS_NMATCH) + hmm->end[hmm->M]; 
    esl_vec_FScale(hmm->t[hmm->M], cp9_TRANS_NMATCH, 1./d);
    hmm->end[hmm->M] /= d;
    
    /* renormalize transitions out of D_M */
    esl_vec_FNorm(hmm->t[hmm->M] + cp9_TRANS_DELETE_OFFSET, cp9_TRANS_NDELETE);	/* delete */
  }

  hmm->begin[1] = (1. - pentry) * (1. - (hmm->t[0][CTMI] + hmm->t[0][CTMD] + hmm->t[0][CTMEL]));
  esl_vec_FSet(hmm->begin+2, hmm->M-1, (pentry * (1.- (hmm->t[0][CTMI] + hmm->t[0][CTMD] + hmm->t[0][CTMEL]))) / (float)(hmm->M-1));
  /* note: hmm->t[0][CTMEL] == 0. (can't locally end from begin) 
   *       and if do_match_local_cm, hmm->t[0][CTMI] and hmm->t[0][CTMD] were just set to 0. 
   */
  
  /* Configure exit.
   * Don't touch hmm->end[hmm->M]
   */

  basep = pexit / (float) (hmm->M-1);
  for (k = 1; k < hmm->M; k++)
    hmm->end[k] = basep / (1. - basep * (float) (k-1));
  CPlan9RenormalizeExits(hmm, 1);
  /*for (k = 1; k <= hmm->M; k++) printf("after renormalizing: end[%d]: %f\n", k, hmm->end[k]);*/

  hmm->flags       &= ~CPLAN9_HASBITS; /* reconfig invalidates log-odds scores */
  hmm->flags       |= CPLAN9_LOCAL_BEGIN; /* local begins now on */
  hmm->flags       |= CPLAN9_LOCAL_END;   /* local ends now on */

  CP9Logoddsify(hmm);
}

/* Function: CPlan9SWConfigEnforce()
 * EPN, Fri Feb  9 05:47:37 2007
 * based on SRE's Plan7SWConfig() from HMMER's plan7.c
 * 
 * Purpose:  Set the alignment independent parameters of
 *           a CM Plan 9 model to hmmsw (Smith/Waterman) configuration.
 *           Same as CPlan9SWConfig but enforces a contiguous subset of
 *           nodes start at x, ending at y must be entered by forbidding
 *           local entries after x and local exits before y.
 *           
 * Args:     hmm    - the CM Plan 9 model w/ data-dep prob's valid
 *           pentry - probability of an internal entry somewhere;
 *                    will be evenly distributed over (enf_start) 
 *                    match states
 *           pexit  - probability of an internal exit somewhere; 
 *                    will be distributed over (M-enf_end) match 
 *                    states.
 *           enf_start_pos - HMM node where enforced node subset begins         
 *           enf_end_pos   - HMM node where enforced node subset ends
 * Return:   (void)
 *           HMM probabilities are modified.
 */
void
CPlan9SWConfigEnforce(CP9_t *hmm, float pentry, float pexit,
		      int enf_start_pos, int enf_end_pos)
{
  float basep;			/* p1 for exits: the base p */
  int   k;			/* counter over states      */

  /* No special (*x* states in Plan 7) states in CM Plan 9 */

  /* Configure entry.
   * To match CM, we enforce the only way out of the B state (M_0)
   * is through a local begin into a match state 
   */
  hmm->t[0][CTMI] = hmm->t[0][CTMD] = hmm->t[0][CTMEL] = 0.;
  hmm->begin[1] = 1. - pentry;
  for (k = 2; k <= enf_start_pos; k++)
    hmm->begin[k] = pentry / (float)(enf_start_pos-1);

  /* OLD WAY (more smith-waterman-like, less CM-like) EPN, Thu Jun 21 15:30:46 2007
     hmm->begin[1] = (1. - pentry) * (1. - (hmm->t[0][CTMI] + hmm->t[0][CTMD])); 
  for (k = 2; k <= enf_start_pos; k++)
    hmm->begin[k] = (pentry * (1.- (hmm->t[0][CTMI] + hmm->t[0][CTMD]))) / (float)(enf_start_pos-1);
  for (k = (enf_start_pos+1); k <= hmm->M; k++)
    hmm->begin[k] = 0.;
  */
    
  /* Configure exit.
   * Don't touch hmm->end[hmm->M]
   */
  if(enf_end_pos == hmm->M) /* no local exit possible */
    basep = 0.0;
  else
    basep = pexit / (float) (hmm->M-enf_end_pos);

  for (k = 0; k < enf_end_pos; k++)
    hmm->end[k] = 0.;
  for (k = enf_end_pos; k < hmm->M; k++)
    hmm->end[k] = basep / (1. - basep * (float) ((k-enf_end_pos)-1));
  CPlan9RenormalizeExits(hmm, 1);
  hmm->flags       &= ~CPLAN9_HASBITS;     /* reconfig invalidates log-odds scores */
  hmm->flags       |= CPLAN9_LOCAL_BEGIN; /* local begins now on */
  hmm->flags       |= CPLAN9_LOCAL_END;   /* local ends now on */
}


/* Function: CPlan9ELConfig()
 * Incept:   EPN, Tue Jun 19 09:50:52 2007
 * 
 * Purpose:  Turn EL local ends in a CM Plan 9 HMM on based on 
 *           the local end probs in the CM. 
 *           
 * Args:     cp9 - the cp9 HMM
 *           cm  - the CM the cp9 was built from
 *                    
 * Return:   (void)
 *           HMM probabilities are modified.
 */
void
CPlan9ELConfig(CP9_t *cp9, CM_t *cm)
{
  /*printf("IN CPlan9ELConfig\n");*/
  /* Contract checks */
  if(cp9->M != cm->clen)  cm_Fail("CPlan9ELConfig(), cp9 and cm model length do not match");
  if(cm->cp9map == NULL)  cm_Fail("CPlan9ELConfig(), cm->cp9map is NULL.\n");
  if(cp9->flags & CPLAN9_EL) cm_Fail("CPlan9ELConfig(), cp9's CPLAN9_EL flag is already up.");
  
  int v;
  int k;                     /* counter over HMM nodes */
  int nd;
  int seen_exit;
  float to_el_prob;
  float norm_factor;
  int   nexits;

  /* If the CM has local ends on, check to make sure all non-zero 
   * local end probabilities in the CM are identical (within reasonable 
   * precision), use that probability to set all HMM transitions to 
   * EL states.
   */
  if(cm->flags & CMH_LOCAL_END) { 
    seen_exit  = FALSE;
    to_el_prob = 0.;
    for(v = 0; v < cm->M; v++) {
      nd = cm->ndidx[v];
      if (((cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
	  cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BEGL_nd ||
	  cm->ndtype[nd] == BEGR_nd) && 
	 cm->ndtype[nd+1] != END_nd) && cm->nodemap[nd] == v) {
      /* this should have a non-zero local end probability */
      if(fabs(cm->end[v] - 0.) < 0.00001) /* non-zero */
	cm_Fail("In CPlan9ELConfig(), CM state %d should have non-zero local end prob, but it doesn't.\n", v);
      if(!seen_exit) {
	to_el_prob = cm->end[v];
	seen_exit  = TRUE;
      }
      else if(fabs(to_el_prob - cm->end[v]) > 0.00001)
	cm_Fail("In CPlan9ELConfig(), not all CM states EL probs are identical.\n");
      }
    }
    if(! seen_exit && cm->nodes != 3) cm_Fail("In CPlan9ELConfig(), CM_LOCAL_END flag up, cm->nodes != 3, but all CM local end probs are zero."); 
  }
  else {
    /* CM_LOCAL_END flag is down, local ends are off in the CM 
     * We figure out what the local end prob would be given cm->pend
     * and set the HMM local end probs based on that. 
     * First, count internal nodes MATP, MATL, MATR, BEGL, BEGR that aren't
     * adjacent to END nodes.
     */
    nexits = 0;
    for (nd = 1; nd < cm->nodes; nd++) {
      if ((cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
	   cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BEGL_nd ||
	   cm->ndtype[nd] == BEGR_nd) && 
	  cm->ndtype[nd+1] != END_nd)
	nexits++;
    } 
    to_el_prob = cm->pend / (float) nexits;
  }

  /* transitions from HMM node 0 to EL is impossible */
  cp9->t[0][CTMEL] = 0.;
  for(k = 1; k <= cp9->M; k++) 
    {
      if(cp9->has_el[k])
	{
	  cp9->t[k][CTMEL] = to_el_prob;
	  norm_factor = 1. - (cp9->t[k][CTMEL] / (1. - cp9->end[k]));
	  cp9->t[k][CTMM] *= norm_factor;
	  cp9->t[k][CTMI] *= norm_factor;
	  cp9->t[k][CTMD] *= norm_factor;
	  /* cp9->end[k] untouched */
	}
    }
  cp9->flags &= ~CPLAN9_HASBITS;	/* clear the log-odds ready flag */

  CP9Logoddsify(cp9);

  cp9->flags |= CPLAN9_EL;          /* EL end locals now on */
  /*debug_print_cp9_params(cp9);*/
  return;
}

/* Function: CPlan9NoEL()
 * Incept:   EPN, Tue Jun 19 09:50:52 2007
 * 
 * Purpose:  Turn EL local ends off in a CM Plan 9 HMM
 *           
 * Args:     cp9 - the HMM
 *                    
 * Return:   (void)
 *           HMM probabilities are modified.
 */
void
CPlan9NoEL(CP9_t *cp9)
{
  /* Contract checks */
  if(! (cp9->flags & CPLAN9_EL)) cm_Fail("CPlan9NoEL, CPLAN9_EL flag is already down.");
  
  int k;                     /* counter over HMM nodes */

  for(k = 0; k <= cp9->M; k++) cp9->t[k][CTMEL] = 0.;
  CPlan9RenormalizeExits(cp9, 1);

  cp9->flags &= ~CPLAN9_HASBITS; /* clear the log-odds ready flag */
  CP9Logoddsify(cp9);

  cp9->flags &= ~CPLAN9_EL;      /* EL local ends now off */

  return;
}

/* Function: CPlan9InitEL()
 * Incept:   EPN, Tue Jun 19 13:10:56 2007
 * 
 * Purpose:  Initialize a CP9 HMM for possible EL local ends
 *           by determining how the EL states should be connected
 *           based on the CM node topology.
 *           
 * Args:     cp9 - the CP9 HMM, built from cm
 *           cm  - the CM the cp9 was built from 
 *
 * Return:   (void)
 */
void
CPlan9InitEL(CP9_t *cp9, CM_t *cm)
{
  int status;
  int k;                     /* counter over HMM nodes */
  int nd;
  int *tmp_el_from_ct;

  if(cm->emap == NULL) cm_Fail("CPlan9InitEL() cm->emap is NULL");

  /* First copy the CM el self transition score/probability: */
  cp9->el_self   = sreEXP2(cm->el_selfsc);
  cp9->el_selfsc = Prob2Score(cp9->el_self, 1.0);

  /* For each HMM node k, we can transit FROM >= 0 EL states from 
   * HMM nodes kp. Determine how many such valid transitions exist
   * from each node, then allocate and fill cp9->el_from_idx[k] and 
   * cp9->el_from_cmnd arrays based on that.
   * This two-pass method saves memory b/c we only allocate for
   * what we'll need.
   */

  /* Initialize to 0 */
  for(k = 0; k <= cp9->M; k++) 
    {
      cp9->el_from_ct[k] = 0;
      cp9->has_el[k] = FALSE;
    }
  cp9->el_from_ct[(cp9->M+1)] = 0; /* special case, we can get to E state from EL states */
    
  /* first pass to get number of valid transitions */
  for(nd = 0; nd < cm->nodes; nd++)
    {
      if ((cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
	   cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BEGL_nd ||
	   cm->ndtype[nd] == BEGR_nd) && 
	  cm->ndtype[nd+1] != END_nd)
	{
	  /*printf("HMM node %d can be reached from HMM node %d's EL state\n", cm->emap->rpos[nd], cm->emap->lpos[nd]);*/
	  cp9->el_from_ct[cm->emap->rpos[nd]]++;
	  cp9->has_el[cm->emap->lpos[nd]] = TRUE;
	}
    }

  /* allocate cp9->el_from_idx[k], cp9->el_from_cmnd for all k */
  for(k = 0; k <= (cp9->M+1); k++) 
    {
      if(cp9->el_from_idx[k] != NULL) /* if !NULL we already filled it, shouldn't happen */
	cm_Fail("ERROR in CPlan9InitEL() el_from_idx has already been initialized\n");
      if(cp9->el_from_ct[k] > 0)
	{
	  ESL_ALLOC(cp9->el_from_idx[k], sizeof(int) * cp9->el_from_ct[k]);
	  ESL_ALLOC(cp9->el_from_cmnd[k],sizeof(int) * cp9->el_from_ct[k]);
	}
      /* else it remains NULL */
    }

  /* now fill in cp9->el_from_idx, we need a new counter array */
  ESL_ALLOC(tmp_el_from_ct, sizeof(int) * (cp9->M+2));
  for(k = 0; k <= (cp9->M+1); k++) 
    tmp_el_from_ct[k] = 0;
  for(nd = 0; nd < cm->nodes; nd++)
    {
      if ((cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
	   cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BEGL_nd ||
	   cm->ndtype[nd] == BEGR_nd) && 
	  cm->ndtype[nd+1] != END_nd)
	{
	  k = cm->emap->rpos[nd];
	  cp9->el_from_idx[k][tmp_el_from_ct[k]] = cm->emap->lpos[nd];
	  cp9->el_from_cmnd[k][tmp_el_from_ct[k]] = nd;
	  tmp_el_from_ct[k]++;
	}
    }

  /* Debugging printfs */
  /*  for(k = 0; k <= (cp9->M+1); k++) 
    {
      for(c = 0; c < cp9->el_from_ct[k]; c++)
	printf("cp9->el_from_idx[%3d][%2d]: %4d\n", k, c, cp9->el_from_idx[k][c]);
      if(cp9->has_el[k])
      printf("node k:%3d HAS an EL!\n", k);
      }*/

  /* Free memory and exit */
  free(tmp_el_from_ct);
  return;

 ERROR:
  cm_Fail("Memory allocation error.");
}


/* Function: CPlan9RenormalizeExits()
 * EPN 05.30.06 based on SRE's Plan7RenormalizeExits() from
 *                       HMMER's plan7.c.
 *
 * Date:     SRE, Fri Aug 14 11:22:19 1998 [St. Louis]
 *
 * Purpose:  Renormalize just the match state transitions;
 *           for instance, after a Config() function has
 *           modified the exit distribution.
 *
 * Args:     hmm - hmm to renormalize
 *           spos   - first consensus column modelled by original
 *                    CP9 HMM the sub CP9 HMM models. Often 1.
 * Returns:  void
 */
void
CPlan9RenormalizeExits(CP9_t *hmm, int spos)
{
  int   k;
  float d;

  /* We can't exit from node 0 so we start renormalizing at node 1 */
  for (k = 1; k < hmm->M; k++)
    {
      if(k != (spos-1)) /* we can't exit from the M_spos-1 */
	{
	  d = esl_vec_FSum(hmm->t[k], 4);
	  /* esl_vec_FScale(hmm->t[k], 4, 1./(d + d*hmm->end[k])); */
	  esl_vec_FScale(hmm->t[k], 4, (1.-hmm->end[k])/d);
	}
    }
  /* Take care of hmm->M node, which is special */
  d = hmm->t[hmm->M][CTMI] + hmm->t[hmm->M][CTMEL]; /* CTMD is IMPOSSIBLE, CTMM is hmm->end[hmm-M] */
  if(! (fabs(d-0.) < eslSMALLX1)) { /* don't divide by d if it's zero */
    hmm->t[hmm->M][CTMI] *= (1.-hmm->end[hmm->M])/d;
    hmm->t[hmm->M][CTMEL] *= (1.-hmm->end[hmm->M])/d;
  }
  return;
}

/* Function:  CP9EnforceHackMatchScores()
 * Incept:    EPN, Fri Feb  9 11:06:31 2007
 *
 * Purpose:   Make all match emissions 0, except those enforce
 *            a specified subsequence (it's assumed the CP9 
 *            is already set up for this enforcement). 
 *
 * Args:      cp9           - the CP9 HMM 
 *            enf_start_pos - first posn of enforced subseq
 *            enf_end_pos   - last  posn of enforced subseq
 * Returns:   (void)
 */
void
CP9EnforceHackMatchScores(CP9_t *cp9, int enf_start_pos, int enf_end_pos)
{
  int k, x;
  for (k = 1; k < enf_start_pos; k++) /* M_0 is the begin state, it's silent */
    for (x = 0; x < cp9->abc->Kp; x++)
      cp9->msc[x][k] = 0.;
  for (k = enf_end_pos+1; k <= cp9->M; k++)
    for (x = 0; x < cp9->abc->Kp; x++)
      cp9->msc[x][k] = 0.;
}

/************************************************************************
 * Functions stolen from HMMER 2.4 for use with CM plan 9 HMMs.
 * Eventually, these should go away, replaced with Easel funcs. 
 * These first 3 were stolen from HMMER:mathsupport.c
 * 
 * Score2Prob()
 * Prob2Score()
 * Scorify()
 * 
 * NOTE: ILogSum() (and auxiliary funcs associated with it) used to be here
 * but moved to logsum.c (EPN, Sat Sep  8 15:49:47 2007)
 */

/* Function: Prob2Score()
 * 
 * Purpose:  Convert a probability to a scaled integer log_2 odds score. 
 *           Round to nearest integer (i.e. note use of +0.5 and floor())
 *           Return the score. 
 */
int
Prob2Score(float p, float null)
{
  if(p == 0.0) return -INFTY;
  else         return (int) floor(0.5 + INTSCALE * sreLOG2(p/null));
}

/* Function: Score2Prob()
 * 
 * Purpose:  Convert an integer log_2 odds score back to a probability;
 *           needs the null model probability, if any, to do the conversion.
 */
float 
Score2Prob(int sc, float null)
{
  if (sc == -INFTY) return 0.;
  else              return (null * sreEXP2((float) sc / INTSCALE));
}

/* Function: Scorify()
 * 
 * Purpose:  Convert a scaled integer log-odds score to a floating
 *           point score for output. (could be a macro but who cares.)
 */
float 
Scorify(int sc)
{
  return ((float) sc / INTSCALE);
}

/* Function: CP9_reconfig2sub()
 * EPN 10.16.06
 * 
 * Purpose:  Given a CM Plan 9 HMM and a start position
 *           (spos) and end position (epos) that a sub CM models, 
 *           reconfigure the HMM so that it can only start in the 
 *           node that models spos (spos_nd) end in the node that 
 *           models epos (epos_nd).
 *
 *           If we're reconfiguring a CP9 HMM that ONLY models the
 *           consensus columns spos to epos, then spos_nd == 1 
 *           and epos_nd == hmm->M, but this is not necessarily true.
 *           We may be reconfiguring a CP9 HMM that models the
 *           full alignment including positions before and/or after
 *           spos and epos. In this case spos_nd == spos and
 *           epos_nd == epos;
 *           
 * Args:     hmm         - the CP9 model w/ data-dep prob's valid
 *           spos        - first consensus column modelled by some original
 *                         full length, template CP9 HMM that 'hmm' models.
 *           epos        - final consensus column modelled by some original
 *                         CP9 HMM that 'hmm' models.
 *           spos_nd     - the node of 'hmm' that models spos.
 *                         (1 if 'hmm' only has (epos-spos+1) nodes 
 *                         (spos if 'hmm' has a node for each column of original aln)
 *           epos_nd     - the node of the 'hmm' in that models epos.
 *                         (hmm->M if 'hmm' only has (epos-spos+1) nodes 
 *                         (epos if 'hmm' has a node for each column of original aln)
 *           orig_phi    - the 2D phi array for the original CP9 HMM.         
 * Return:   (void)
 *           HMM probabilities are modified.
 */
void
CP9_reconfig2sub(CP9_t *hmm, int spos, int epos, int spos_nd,
		 int epos_nd, double **orig_phi)
{
  /* Make the necessary modifications. Since in cmalign --sub mode this
   * function will be called potentially once for each sequence, we 
   * don't want to call CP9Logoddsify(), but rather only logoddsify
   * the parameters that are different.
   */

  /* Configure entry.
   * Exactly 3 ways to start, B->M_1 (hmm->begin[1]), B->I_0 (hmm->t[0][CTMI]),
   *                      and B->D_1 (hmm->t[0][CTMD])
   */
  /* prob of starting in M_spos is (1. - prob of starting in I_spos-1) as there is no D_spos-1 -> M_spos trans */
      
  if(spos > 1)
    {
      hmm->begin[spos_nd] = 1.-((orig_phi[spos-1][HMMINSERT] * (1. - hmm->t[spos-1][CTII])) + 
			        (orig_phi[spos  ][HMMDELETE] - (orig_phi[spos-1][HMMINSERT] * hmm->t[spos-1][CTID])));
      hmm->t[spos_nd-1][CTMI] =   (orig_phi[spos-1][HMMINSERT] * (1. - hmm->t[spos-1][CTII]));
      hmm->t[spos_nd-1][CTMD] =    orig_phi[spos  ][HMMDELETE] - (orig_phi[spos-1][HMMINSERT] * hmm->t[spos-1][CTID]);
      hmm->t[spos_nd-1][CTMM] = 0.; /* probability of going from B(M_0) to M_1 is begin[1] */
      hmm->t[spos_nd-1][CTMEL] = 0.; /* can't go to EL from B(M_0) */
      hmm->t[spos_nd-1][CTDM] = 0.; /* D_0 doesn't exist */
      hmm->t[spos_nd-1][CTDI] = 0.; /* D_0 doesn't exist */
      hmm->t[spos_nd-1][CTDD] = 0.; /* D_0 doesn't exist */
      
      hmm->bsc[spos_nd]       = Prob2Score(hmm->begin[1], 1.0);

      hmm->tsc[CTMM][spos_nd-1] = -INFTY; /* probability of going from B(M_0) to M_1 is begin[1] */
      hmm->tsc[CTMEL][spos_nd-1] = -INFTY; 
      hmm->tsc[CTDM][spos_nd-1] = -INFTY; /* D_0 doesn't exist */
      hmm->tsc[CTDI][spos_nd-1] = -INFTY; /* D_0 doesn't exist */
      hmm->tsc[CTDD][spos_nd-1] = -INFTY; /* D_0 doesn't exist */
      
      hmm->tsc[CTMI][spos_nd-1] = Prob2Score(hmm->t[spos_nd-1][CTMI], 1.0);
      hmm->tsc[CTMD][spos_nd-1] = Prob2Score(hmm->t[spos_nd-1][CTMD], 1.0);
    }

  if(epos < hmm->M)
    {
      hmm->end[epos_nd]      = hmm->t[epos][CTMM] + hmm->t[epos][CTMD];
      hmm->t[epos_nd][CTDM] += hmm->t[epos][CTDD];
      hmm->t[epos_nd][CTIM] += hmm->t[epos][CTID];
      hmm->t[epos_nd][CTMM]  = 0.; /* M->E is actually end[M] */
      hmm->t[epos_nd][CTMEL]  = 0.; 
      hmm->t[epos_nd][CTMD]  = 0.; /* D_M+1 doesn't exist */
      hmm->t[epos_nd][CTDD]  = 0.; /* D_M+1 doesn't exist */
      hmm->t[epos_nd][CTID]  = 0.; /* D_M+1 doesn't exist */
      
      hmm->esc[epos_nd]       = Prob2Score(hmm->end[epos_nd], 1.0);
      hmm->tsc[CTDM][epos_nd] = Prob2Score(hmm->t[epos_nd][CTDM], 1.0);
      hmm->tsc[CTIM][epos_nd] = Prob2Score(hmm->t[epos_nd][CTIM], 1.0);
      hmm->tsc[CTMM][epos_nd] = -INFTY; /* M->E is actually end[M] */
      hmm->tsc[CTMEL][epos_nd] = -INFTY; 
      hmm->tsc[CTMD][epos_nd] = -INFTY; /* D_M+1 doesn't exist */
      hmm->tsc[CTDD][epos_nd] = -INFTY; /* D_M+1 doesn't exist */
      hmm->tsc[CTID][epos_nd] = -INFTY; /* D_M+1 doesn't exist */
    }
  hmm->flags |= CPLAN9_HASBITS;	/* raise the log-odds ready flag */

  return;
}

/* cm_modelconfig.c
 * SRE, Wed May  8 14:30:38 2002 [St. Louis]
 * SVN $Id$
 * 
 * Configuring a model into different global or local modes.
 * 
 ******************************************************************
 * @LICENSE@
 ******************************************************************
 */

#include "esl_config.h"
#include "config.h"

#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_stack.h"
#include "esl_vectorops.h"

#include "funcs.h"
#include "structs.h"

/*
 * Function: ConfigCM()
 * Date:     EPN, Thu Jan  4 06:36:09 2007
 * Purpose:  Configure a CM for alignment or search based on cm->config_opts,
 *           cm->align_opts and cm->search_opts. 
 *           Always builds CP9 HMM (it's fast).
 *           Calculates query dependent bands (QDBs) if nec.
 *           QDBs can also be passed in. 
 * 
 * Args:     CM           - the covariance model
 *           preset_dmin  - supplied dmin values, NULL if none
 *           preset_dmax  - supplied dmax values, NULL if none
 *
 * Returns:   <eslOK> on success.
 *            <eslEINVAL> on contract violation.
 *            <eslEMEM> on memory allocation error.
 */
int 
ConfigCM(CM_t *cm, int *preset_dmin, int *preset_dmax)
{
  int status;
  float swentry, swexit;
  int do_calc_qdb   = FALSE;
  int do_preset_qdb = FALSE;
  int v;
  
  /* Check if we need to calculate QDBs and/or build a CP9 HMM. */
  if(cm->config_opts & CM_CONFIG_QDB) {
    if(preset_dmin == NULL && preset_dmax == NULL) 
      do_calc_qdb   = TRUE;
    else 
      do_preset_qdb = TRUE;
  }

  /* Build the CP9 HMM and associated data */
  /* IMPORTANT: do this before setting up CM for local mode
   * if we already have these, free them (wasteful but safe, 
   * and not a big deal b/c we 'should' only call ConfigCM() once
   * per CM.  
   */
  if(cm->cp9map     != NULL) FreeCP9Map(cm->cp9map);
  if(cm->cp9b       != NULL) FreeCP9Bands(cm->cp9b);
  if(cm->cp9        != NULL) FreeCPlan9(cm->cp9);
  if(cm->cp9_mx     != NULL) FreeCP9Matrix(cm->cp9_mx);
  if(cm->cp9_bmx    != NULL) FreeCP9Matrix(cm->cp9_bmx);

  if(!(build_cp9_hmm(cm, &(cm->cp9), &(cm->cp9map), FALSE, 0.0001, 0))) cm_Fail("Couldn't build a CP9 HMM from the CM\n");
  cm->cp9b = AllocCP9Bands(cm, cm->cp9);
  /* create the CP9 matrices, we init to 1 row, which is tiny so it's okay
   * that we have two of them, we only grow them as needed, cp9_bmx is 
   * only needed if we're doing Forward -> Backward -> Posteriors.
   */
  cm->cp9_mx  = CreateCP9Matrix(1, cm->cp9->M);
  cm->cp9_bmx = CreateCP9Matrix(1, cm->cp9->M);
  cm->flags |= CMH_CP9; /* raise the CP9 flag */
  
  /* Possibly configure the CM for local alignment. */

  if (cm->config_opts & CM_CONFIG_LOCAL)
    { 
      ConfigLocal(cm, cm->pbegin, cm->pend);
      CMLogoddsify(cm);
    }
  /* Possibly configure the CP9 for local alignment
   * Note: CP9 local/glocal config does not have to necessarily match CM config,
   *       although all the executables are currently setup so that when CM goes
   *       local HMM goes local.
   */
  if((cm->config_opts & CM_CONFIG_HMMLOCAL) || (cm->align_opts  & CM_ALIGN_SUB)) {
    if(cm->align_opts & CM_ALIGN_SUB) {
      /* To get spos and epos for the sub_cm, 
       * we config the HMM to local mode with equiprobable start/end points.
       * and we DO NOT make I_0 and I_M unreachable. These states map to ROOT_IL and ROOT_IR,
       * which are unreachable in a locally configured CM. 
       */
      swentry= ((cm->cp9->M)-1.)/cm->cp9->M; /* all start pts equiprobable, including 1 */
      swexit = ((cm->cp9->M)-1.)/cm->cp9->M; /* all end   pts equiprobable, including M */
      CPlan9SWConfig(cm->cp9, swentry, swexit, FALSE); /* FALSE means don't make I_0, D_1, I_M unreachable */
    }
    else { 
      /* if we're setting up the HMM for local search/alignment but NOT for a sub CM alignment, 
       * we DO make I_0 and I_M unreachable. These states map to ROOT_IL and ROOT_IR,
       * which are unreachable in a locally configured CM. 
       */
      /* CPlan9CMLocalBeginConfig(cm); */
      swentry = cm->pbegin;
      swexit  = cm->pbegin;
      /* swentry= ((cm->cp9->M)-1.)/cm->cp9->M; *//* all start pts equiprobable, including 1 */
      /* swexit = ((cm->cp9->M)-1.)/cm->cp9->M; *//* all end   pts equiprobable, including M */
      CPlan9SWConfig(cm->cp9, swentry, swexit, TRUE); /* TRUE means do make I_0, D_1, I_M unreachable to match the CM */
    }
    CP9Logoddsify(cm->cp9);
  }
  if(cm->config_opts & CM_CONFIG_HMMEL)
    CPlan9ELConfig(cm);

  /*
    FILE *fp;
    fp = fopen("temphmm1" ,"w");
    debug_print_cp9_params(fp, cm->cp9, TRUE);
    fclose(fp);
    
    fp = fopen("tempcm1" ,"w");
    debug_print_cm_params(fp, cm);
    fclose(fp);
  */

  /* If nec, set up the query dependent bands, this has to be done after 
   * local is set up because we want to consider local begins, but NOT local
   * ends. This is inefficient, local ends are set up twice currently, once
   * before QDB calc, then turned off in BandCalculationEngine() then 
   * back on. This should be fixed. */
  if (do_calc_qdb)
    {
      if(cm->flags & CMH_QDB) cm_Fail("ERROR in ConfigCM() CM already has QDBs\n");
      ConfigQDB(cm);
    }
  else if(do_preset_qdb)
    {
      if(cm->flags & CMH_QDB) cm_Fail("ERROR in ConfigCM() CM already has QDBs\n");
      ESL_ALLOC(cm->dmin, sizeof(int) * cm->M);
      ESL_ALLOC(cm->dmax, sizeof(int) * cm->M);
      for(v = 0; v < cm->M; v++)
	{
	  cm->dmin[v] = preset_dmin[v];
	  cm->dmax[v] = preset_dmax[v];
	}
      /* Set W as dmax[0], we're wasting time otherwise, looking at
       * hits that are bigger than we're allowing with QDB. */
      cm->W = cm->dmax[0];
      cm->flags |= CMH_QDB; /* raise the QDB flag */
    }
  else { /* we didn't set up QDBs, but we still need to set cm->W, we 
	  * set it as dmax[0] from the QDB calculation using cm->beta, 
	  * but don't save dmin/dmax.
	  * (This if, else if, esle could definitely be better organized) 
	  */
    int safe_windowlen = cm->clen * 2;
    int *dmin, *dmax;
    while(!(BandCalculationEngine(cm, safe_windowlen, cm->beta, FALSE, &(dmin), &(dmax), NULL, NULL)))
      {
	free(dmin);
	free(dmax);
	safe_windowlen *= 2;
	if(safe_windowlen > (cm->clen * 1000))
	  cm_Fail("ConfigCM(), safe_windowlen big: %d\n", safe_windowlen);
      }
    cm->W = dmax[0];
    free(dmin);
    free(dmax);
  }

  /*
  fp = fopen("temphmm2" ,"w");
  debug_print_cp9_params(fp, cm->cp9, TRUE);
  fclose(fp);

  fp = fopen("tempcm2" ,"w");
  debug_print_cm_params(fp, cm);
  fclose(fp);
  */

  /*cm_Fail("done.\n");*/
  /* We need to ensure that cm->el_selfsc * W >= IMPOSSIBLE
   * (cm->el_selfsc is the score for an EL self transition) This is
   * done because we potentially multiply cm->el_selfsc * W, and add
   * that to IMPOSSIBLE. 
   */
  if((cm->el_selfsc * cm->W) < IMPOSSIBLE)
    { 
      cm->el_selfsc = (IMPOSSIBLE / (cm->W+1));
      cm->iel_selfsc = -INFTY;
    }
  /* Potentially, overwrite transitions with non-probabilistic 
   * RSEARCH transitions, we do this after setting up QDBs and
   * the CP9 HMM, so they'll correspond to the probabilistic 
   * transitions that existed prior to overwriting with RSEARCH
   * transitions. Transitions scores are overwritten in CMLogoddsify() 
   */
  CMLogoddsify(cm);

  /*debug_print_cm_params(stdout, cm);
    debug_print_cp9_params(stdout, cm->cp9, TRUE);*/
  return eslOK;

 ERROR:
  cm_Fail("Memory allocation error.");
  return status; /* NOTREACHED */
}

/*
 * Function: ConfigLocal
 * Purpose:  Configure a CM for local alignment by spreading 
 *           p_internal_start local entry probability evenly
 *           across all internal nodes, and by spreading
 *           p_internal_exit local exit probability evenly
 *           across all internal nodes.
 * 
 * Args:     cm               - the covariance model
 *           p_internal_start - prob mass to spread for local begins
 *           p_internal_exit  - prob mass to spread for local ends
 */        

void
ConfigLocal(CM_t *cm, float p_internal_start, float p_internal_exit)
{
  int status;
  int v;			/* counter over states */
  int nd;			/* counter over nodes */
  int nstarts;			/* number of possible internal starts */
  int had_scanmatrix;           /* true if CM had a scan matrix when function was entered */

  /* contract check */
  if(cm->flags & CMH_LOCAL_BEGIN)
    cm_Fail("ERROR in ConfigLocal(), CMH_LOCAL_BEGIN flag already up.\n");
  if(cm->flags & CMH_LOCAL_END)
    cm_Fail("ERROR in ConfigLocal(), CMH_LOCAL_END flag already up.\n");

  if(cm->flags & CMH_SCANMATRIX) had_scanmatrix = TRUE;
  else had_scanmatrix = FALSE;

  /*****************************************************************
   * Internal entry.
   *****************************************************************/
  /* Count "internal" nodes: MATP, MATL, MATR, and BIF nodes.
   * Ignore all start nodes, and also node 1 (which is always the
   * "first" node and gets an entry prob of 1-p_internal_start).
   */
  nstarts = 0;
  for (nd = 2; nd < cm->nodes; nd++) {
    if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
    	cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BIF_nd) 
      nstarts++;
  }

  /* Zero everything.
   */
  for (v = 0; v < cm->M; v++)  cm->begin[v] = 0.;

  /* Erase the previous transition p's from node 0. The only
   * way out of node 0 is going to be local begin transitions
   * from the root v=0 directly to MATP_MP, MATR_MR, MATL_ML,
   * and BIF_B states.
   *
   * EPN, Wed Feb 14 13:56:02 2007: First we want to save
   * the transition probs we're about to zero, so we can
   * globalize this model later if we have to.
   * cm->root_trans is NULL prior to this, so 
   * if we call ConfigGlobal() w/o executing this code,
   * we'll die (with an EASEL contract exception).
   */
  if(cm->root_trans == NULL) /* otherwise they've already been set */
    {
      ESL_ALLOC(cm->root_trans, sizeof(float) * cm->cnum[0]);
      for (v = 0; v < cm->cnum[0]; v++)
	cm->root_trans[v] = cm->t[0][v];
    }
  for (v = 0; v < cm->cnum[0]; v++) cm->t[0][v] = 0.;

  /* Node 1 gets prob 1-p_internal_start.
   */
  cm->begin[cm->nodemap[1]] = 1.-p_internal_start;

  /* Remaining nodes share p_internal_start.
   */
  for (nd = 2; nd < cm->nodes; nd++) {
    if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
    	cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BIF_nd)  
      cm->begin[cm->nodemap[nd]] = p_internal_start/(float)nstarts;
  }
  cm->flags |= CMH_LOCAL_BEGIN;
  
  /*****************************************************************
   * Internal exit.
   *****************************************************************/
  ConfigLocalEnds(cm, p_internal_exit);

  /* new local probs invalidate log odds scores and QDBs */
  cm->flags &= ~CMH_BITS;
  /* recalc QDBs if they exist */
  if(cm->flags & CMH_QDB) {
    free(cm->dmin);
    free(cm->dmax);
    cm->dmin = NULL;
    cm->dmax = NULL;
    cm->flags &= ~CMH_QDB;
    ConfigQDB(cm);
  }      
  /* ConfigQDB should rebuild scan matrix, if it existed */
  if((! (cm->flags & CMH_SCANMATRIX)) && had_scanmatrix)
     cm_Fail("ConfigLocal(), CM had a scan matrix, but ConfigQDB didn't rebuild it.");
     
  CMLogoddsify(cm);
  return;

 ERROR:
  cm_Fail("Memory allocation error.");
}

/*
 * Function: ConfigGlobal
 * Purpose:  Configure a CM in local alignment mode to global
 *           alignment mode.
 *
 * Args:     CM               - the covariance model
 */        

void
ConfigGlobal(CM_t *cm)
{
  int v;			/* counter over states */

  /*printf("in configGlobal\n");*/
  /* Contract check: local begins MUST be active, if not then cm->root_trans (the 
   * transition probs from state 0 before local configuration) will be NULL, 
   * so we can't copy them back into cm->t[0], which is a problem. This is fragile. */
  if(!(cm->flags & CMH_LOCAL_BEGIN))
    cm_Fail("ERROR in ConfigGlobal() trying to globally configure a CM that has no local begins.");
  if(!(cm->flags & CMH_LOCAL_END))
    cm_Fail("ERROR in ConfigGlobal() trying to globally configure a CM that has no local ends.");
  if(cm->root_trans == NULL)
    cm_Fail("ERROR in ConfigGlobal() cm->root_trans NULL. CM must have been configured with local begins before we can configure it back to global");
  
  /*****************************************************************
   * Make local begins impossible
   *****************************************************************/
  for (v = 0; v < cm->M; v++)  cm->begin[v] = 0.;
  /* Now reset transitions out of ROOT_S to their initial state, 
   * ConfigLocal() zeroes these guys, but they've been saved in 
   * the CM data structure in (C
   * in the CM data structure. */
  for (v = 0; v < cm->cnum[0]; v++)  cm->t[0][v] = cm->root_trans[v];
  /*printf("in ConfigGlobal, printing transitions from the root.\n");
    for (v = 0; v < cm->cnum[0]; v++)  printf("cm->t[0][v:%d]: %f\n", v, cm->t[0][v]);*/

  cm->flags &= ~CMH_LOCAL_BEGIN; /* drop the local begin flag */
  
  /*****************************************************************
   * Make local ends impossible
   *****************************************************************/
  ConfigNoLocalEnds(cm);

  /* new probs invalidate log odds scores and QDB */
  cm->flags &= ~CMH_BITS;
  /* Recalc QDBs if they exist */
  if(cm->flags & CMH_QDB) {
    free(cm->dmin);
    free(cm->dmax);
    cm->dmin = NULL;
    cm->dmax = NULL;
    cm->flags &= ~CMH_QDB;
    ConfigQDB(cm);
  }      
  /* free and rebuild scan matrix to correspond to new QDBs, if it exists */
  if(cm->flags & CMH_SCANMATRIX) {
    int do_float = cm->smx->flags & cmSMX_HAS_FLOAT;
    int do_int   = cm->smx->flags & cmSMX_HAS_INT;
    cm_FreeScanMatrixForCM(cm);
    cm_CreateScanMatrixForCM(cm, do_float, do_int);
  }

  CMLogoddsify(cm);
  return;
}

/**************************************************************************
 * EPN 10.02.06
 * Function: ConfigNoLocalEnds()
 *
 * Purpose:  Set the probability of local ends to 0.0 for all states.
 *           This function was introduced for use in BandCalculationEngine()
 *           because allowing local ends when calculating bands dramatically
 *           widens all bands and decreases search acceleration. So this
 *           is the ad-hoc fix.                    
 * 
 * Args:    
 * CM_t *cm               - the CM
 * Returns: (void) 
 */
void
ConfigNoLocalEnds(CM_t *cm)
{
  int v;			/* counter over states */
  int nd;                       /* counter over nodes  */

  /* Contract check */
  if(!(cm->flags & CMH_LOCAL_END))
    cm_Fail("ERROR in ConfigNoLocalEnds() CMH_LOCAL_END flag already down.\n");

  for (v = 0; v < cm->M; v++) cm->end[v] = 0.;
  /* Now, renormalize transitions */
  for (nd = 1; nd < cm->nodes; nd++) {
    if ((cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
	 cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BEGL_nd ||
	 cm->ndtype[nd] == BEGR_nd) && 
	cm->ndtype[nd+1] != END_nd)
      {
	v = cm->nodemap[nd];
	esl_vec_FNorm(cm->t[v], cm->cnum[v]);
      }
  }
  /* Disable the local end probs in the CP9 */
  if((cm->flags |= CMH_CP9) && (cm->cp9->flags |= CPLAN9_EL))
    CPlan9NoEL(cm);

  cm->flags &= ~CMH_LOCAL_END; /* turn off local ends flag */
  /* new probs invalidate log odds scores */
  cm->flags &= ~CMH_BITS;
  /* local end changes don't invalidate QDBs, which means
   * they don't affect the ScanMatrix either. 
   */

  return;
}

/*
 * Function: ConfigLocalEnds()
 * Date:     EPN 10.02.06
 * Purpose:  Given a probability of local ends, spread the probability of
 *           local ends evenly across all states from which local ends are
 *           permitted (see code).
 */

void
ConfigLocalEnds(CM_t *cm, float p_internal_exit)
{
  int v;			/* counter over states */
  int nd;			/* counter over nodes */
  int nexits;			/* number of possible internal ends */
  float denom;

  /* Contract check */
  if(cm->flags & CMH_LOCAL_END)
    cm_Fail("ERROR in ConfigLocalEnds() CMH_LOCAL_END flag already up.\n");

  /* Count internal nodes MATP, MATL, MATR, BEGL, BEGR that aren't
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
  /* Spread the exit probability across internal nodes.
   * Currently does not compensate for the decreasing probability
   * of reaching a node, the way HMMER does: therefore the probability
   * of exiting at later nodes is actually lower than the probability 
   * of exiting at earlier nodes. This should be a small effect.
   */
  for (v = 0; v < cm->M; v++) cm->end[v] = 0.;
  for (nd = 1; nd < cm->nodes; nd++) {
    if ((cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
	 cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BEGL_nd ||
	 cm->ndtype[nd] == BEGR_nd) && 
	cm->ndtype[nd+1] != END_nd)
      {
	v = cm->nodemap[nd];
	cm->end[v] = p_internal_exit / (float) nexits;
				/* renormalize the main model transition distribution */
	denom = esl_vec_FSum(cm->t[v], cm->cnum[v]);
	denom += cm->end[v];
	esl_vec_FScale(cm->t[v], cm->cnum[v], 1./denom);
	/* cm->t[v] vector will purposefully no longer sum to 1., 
	 * if we were to append cm->end[v] as a new number in the vector, it would sum to 1. 
	 */
      }
  }

  cm->flags |= CMH_LOCAL_END;

  /* new probs invalidate log odds scores */
  cm->flags &= ~CMH_BITS;
  /* local end changes don't invalidate QDBs, which means
   * they don't affect the ScanMatrix either. 
   */

  return;
}

/*
 * Function: ConfigLocal_DisallowELEmissions()
 * Purpose:  Silence the EL state.
 */
void
ConfigLocal_DisallowELEmissions(CM_t *cm)
{
  /* Set the EL self transition score to as close to IMPOSSIBLE 
   * as we can while still guaranteeing we won't get underflow errors.
   * we need cm->el_selfsc * W >= IMPOSSIBLE 
   * because we will potentially multiply cm->el_selfsc * W, and add that to 
   * 2 * IMPOSSIBLE, and IMPOSSIBLE must be > -FLT_MAX/3 so we can add it together 3 
   * times (see structs.h). 
   */
  cm_Fail("ConfigLocal_DisallowELEmissions is deprecated.");
  cm->el_selfsc = (IMPOSSIBLE / (cm->W+1));
  cm->iel_selfsc = -INFTY; 
  cm->iel_selfsc = -100 * INTSCALE; 
  return;
}

/*
 * Function: ConfigQDB
 * Date:     EPN, Thu May  3 14:37:09 2007
 * Purpose:  Configure a CM's query dependent bands (QDBs).
 * Args:
 *           CM           - the covariance model
 */
int
ConfigQDB(CM_t *cm)
{
  int safe_windowlen;

  /* Contract check */
  if(cm->flags & CMH_QDB)
    cm_Fail("ERROR in ConfigQDB() CMH_QDB flag already up.\n");

  safe_windowlen = cm->W * 2;
  if(cm->dmin != NULL) {
    free(cm->dmin);
    cm->dmin = NULL;
  }
  if(cm->dmax != NULL) {
    free(cm->dmax);
    cm->dmax = NULL;
  }
  /*debug_print_cm_params(cm);*/
  while(!(BandCalculationEngine(cm, safe_windowlen, cm->beta, 0, &(cm->dmin), &(cm->dmax), NULL, NULL)))
    {
      free(cm->dmin);
      free(cm->dmax);
      cm->dmin = NULL;
      cm->dmax = NULL;
      safe_windowlen *= 2;
      if(safe_windowlen > (cm->clen * 1000))
	cm_Fail("ERROR safe_windowlen big: %d\n", safe_windowlen);
    }
  /* Set W as dmax[0], we're wasting time otherwise, looking at
   * hits that are bigger than we're allowing with QDB. */
  cm->W = cm->dmax[0];
  cm->flags |= CMH_QDB; /* raise the QDB flag */

  /* free and rebuild scan matrix to correspond to new QDBs, if it exists */
  if(cm->flags & CMH_SCANMATRIX) {
    int do_float = cm->smx->flags & cmSMX_HAS_FLOAT;
    int do_int   = cm->smx->flags & cmSMX_HAS_INT;
    cm_FreeScanMatrixForCM(cm);
    cm_CreateScanMatrixForCM(cm, do_float, do_int);
  }

  CMLogoddsify(cm); /* QDB calculation invalidates log odds scores */
  return eslOK;
}

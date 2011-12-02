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
#include "p7_config.h"
#include "config.h"

#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_stack.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "funcs.h"
#include "structs.h"

/* Function: ConfigCM()
 * Date:     EPN, Thu Jan  4 06:36:09 2007
 * Purpose:  Configure a CM for alignment or search based on cm->config_opts,
 *           cm->align_opts and cm->search_opts. 
 *           Always builds CP9 HMM(s) (it's fast).
 *           Calculates query dependent bands (QDBs) if nec.
 *           QDBs can also be passed in. 
 * 
 *           If <sub_mother_cm> and <sub_mother_map> are non-NULL, <cm> is 
 *           a sub CM constructed from <sub_mother_cm>. In this case, we're
 *           doing alignment and are constructing a new sub CM for each target
 *           sequence so running time should be minimized. Special functions
 *           for building the CP9 HMMs and for logoddsifying the model are called
 *           that are faster than the normal versions b/c they can just copy some
 *           of the parameters of the mother model instead of calc'ing them.
 *          
 * Args:     cm            - the covariance model
 *           errbuf        - for error messages
 *           always_calc_W - TRUE to always calculate W even if we're not calcing
 *                           QDBs, FALSE to only calc W if we're calcing QDBs
 *           mother_cm     - if non-NULL, <cm> is a sub CM construced from 
 *                           <mother_cm>. In this case we use <mother_map>
 *                           to help streamline the two steps that dominate the 
 *                           running time of this function (b/c speed is an issue):
 *                           building a cp9 HMM, and logoddsifying the model.
 *           mother_map    - must be non-NULL iff <mother_cm> is non-NULL, the 
 *                           map from <cm> to <mother_cm>.
 *
 * Returns:   <eslOK> on success.
 *            <eslEINVAL> on contract violation.
 *            <eslEMEM> on memory allocation error.
 */
int 
ConfigCM(CM_t *cm, char *errbuf, int always_calc_W, CM_t *mother_cm, CMSubMap_t *mother_map)
{
  int status;
  float swentry, swexit;
  int have_mother;
  have_mother = (mother_cm != NULL && mother_map != NULL)  ? TRUE : FALSE;

  /* contract check */
  if(mother_cm != NULL && mother_map == NULL)  ESL_FAIL(eslEINCOMPAT, errbuf, "ConfigCM(), mother_cm != NULL but mother_map == NULL (both must be NULL or both non-NULL).");
  if(mother_cm == NULL && mother_map != NULL)  ESL_FAIL(eslEINCOMPAT, errbuf, "ConfigCM(), mother_cm == NULL but mother_map != NULL (both must be NULL or both non-NULL).");
  if(have_mother && cm->config_opts & CM_CONFIG_LOCAL) ESL_FAIL(eslEINCOMPAT, errbuf, "ConfigCM(), configuring a sub CM, but CM_CONFIG_LOCAL config flag up.");

  /* TEMP */
  ESL_STOPWATCH *w;
  w = esl_stopwatch_Create();
  char          time_buf[128];  /* string for printing timings (safely holds up to 10^14 years) */
  /* TEMP */
  
  /* Build the emitmap, if necessary */
  if(cm->emap == NULL) cm->emap = CreateEmitMap(cm);

  /* Allocate the truncated HMM banded matrices.
   * The non-truncated versions (cm->hbmx and cm->ohbmx)
   * are allocated in CreateCMBody(), but these require
   * knowledge of B states so are allocated here.
   */
  if(cm->trhbmx  != NULL) cm_tr_hb_mx_Destroy(cm->trhbmx);
  if(cm->trohbmx != NULL) cm_tr_hb_mx_Destroy(cm->trohbmx);
  cm->trhbmx  = cm_tr_hb_mx_Create(cm);
  cm->trohbmx = cm_tr_hb_mx_Create(cm);

  /* Allocate the emit matrices for alignment traceback, 
   * Initially these are small, and only grown as needed.
   */
  if(cm->ehbmx   != NULL) cm_hb_emit_mx_Destroy(cm->ehbmx);
  if(cm->trehbmx != NULL) cm_tr_hb_emit_mx_Destroy(cm->trehbmx);
  cm->ehbmx   = cm_hb_emit_mx_Create(cm);
  cm->trehbmx = cm_tr_hb_emit_mx_Create(cm);

  /* Allocate the shadow matrices for alignment traceback, 
   * Initially these are small, and only grown as needed.
   */
  if(cm->shhbmx   != NULL) cm_hb_shadow_mx_Destroy(cm->shhbmx);
  if(cm->trshhbmx != NULL) cm_tr_hb_shadow_mx_Destroy(cm->trshhbmx);
  cm->shhbmx   = cm_hb_shadow_mx_Create(cm);
  cm->trshhbmx = cm_tr_hb_shadow_mx_Create(cm);

  /* Build the CP9 HMM and associated data */
  /* IMPORTANT: do this before setting up CM for local mode
   * if we already have these, free them (wasteful but safe, 
   * and not a big deal b/c we 'should' only call ConfigCM() once
   * per CM.  
   */
  if(cm->cp9map     != NULL) { FreeCP9Map(cm->cp9map);     cm->cp9map = NULL; }
  if(cm->cp9b       != NULL) { FreeCP9Bands(cm->cp9b);     cm->cp9b   = NULL; }
  if(cm->cp9loc     != NULL) { FreeCPlan9(cm->cp9loc);     cm->cp9loc = NULL; }
  if(cm->cp9glb     != NULL) { FreeCPlan9(cm->cp9glb);     cm->cp9glb = NULL; }
  if(cm->cp9_mx     != NULL) { FreeCP9Matrix(cm->cp9_mx);  cm->cp9_mx = NULL; }
  if(cm->cp9_bmx    != NULL) { FreeCP9Matrix(cm->cp9_bmx); cm->cp9_bmx = NULL; }
  cm->flags &= ~CMH_CP9GLB;
  cm->flags &= ~CMH_CP9LOC;
  
  esl_stopwatch_Start(w);  
  if(have_mother) { 
    if(!(sub_build_cp9_hmm_from_mother(cm, errbuf, mother_cm, mother_map, &(cm->cp9loc), &(cm->cp9map), FALSE, 0.0001, 0))) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "Couldn't build a CP9 HMM from the sub CM and it's mother\n");
  }
  else { 
    if(!(build_cp9_hmm(cm, &(cm->cp9loc), &(cm->cp9map), FALSE, 0.0001, 0))) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "Couldn't build a CP9 HMM from the CM\n");
  }

  esl_stopwatch_Stop(w); 
  FormatTimeString(time_buf, w->user, TRUE);
#if PRINTNOW
  fprintf(stdout, "\tcp9 build time        %11s\n", time_buf);
#endif

  cm->cp9b = AllocCP9Bands(cm->M, cm->cp9loc->M);
  /* create the CP9 matrices, we init to 1 row, which is tiny so it's okay
   * that we have two of them, we only grow them as needed, cp9_bmx is 
   * only needed if we're doing Forward -> Backward -> Posteriors.
   */
  cm->cp9_mx  = CreateCP9Matrix(1, cm->cp9loc->M);
  cm->cp9_bmx = CreateCP9Matrix(1, cm->cp9loc->M);
  cm->flags |= CMH_CP9LOC; /* raise the CP9 flag */
  
  /* build the p7 and profiles from the cp9 */
  if(cm->mlp7    != NULL) p7_hmm_Destroy(cm->mlp7);
#if 0 
  if((status = BuildP7HMM_MatchEmitsOnly(cm, &(cm->mlp7)) != eslOK)) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "Couldn't build a p7 HMM from the CM\n");
#endif
  /* EPN, Tue Nov  9 09:46:48 2010 */
  if((status = cm_cp9_to_p7(cm, cm->cp9loc)) != eslOK) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "Couldn't build a p7 from the CM");

  /* Possibly configure the CM for local alignment. */
  if (cm->config_opts & CM_CONFIG_LOCAL) { 
    ConfigLocal(cm, cm->pbegin, cm->pend);
    CMLogoddsify(cm);
    cm->cp9glb = NULL;        /* actually, this will already be NULL, but just to be clear */
    cm->flags &= ~CMH_CP9GLB; 
  }
  else { /* CM will stay configured in global mode */
    /* Clone the globally configured CP9 HMM cm->cp9 before its put into local mode,
     * the global version will be cm->cp9glb and will be used for standard
     * (non-truncated) alignment. cm->cp9 is still needed for truncated alignment.
     */
    if((cm->cp9glb = cp9_Clone(cm->cp9loc)) == NULL) ESL_FAIL(eslFAIL, errbuf, "couldn't clone cm->cp9loc into cm->cp9glb");
    cm->flags |= CMH_CP9GLB; 
  }

  /* Always configure one of the CP9 HMMs for local begins/ends alignment (cm->cp9) 
   * we'll use this for truncated alignment regardless of whether the CM is
   * in local or glocal mode. Optionally turn on ELs for cm->cp9. Always leave
   * cm->cp9glb untouched, it is either NULL (if the CM is not in local mode) 
   * or it remains in glocal mode and is used for standard (non-truncated) 
   * alignment.
   */
  swentry = cm->pbegin;
  swexit  = cm->pbegin; 
  CPlan9SWConfig(cm->cp9loc, swentry, swexit, 
		 FALSE, /* allow I_0, D_1, I_M to be reachable */
		 cm->ndtype[1]); 
  CP9Logoddsify(cm->cp9loc);
  if(cm->config_opts & CM_CONFIG_HMMEL) { 
    CPlan9ELConfig(cm);
  }

  /* If nec, set up the query dependent bands, this has to be done after 
   * local is set up because we want to consider local begins, but NOT local
   * ends. This is inefficient, local ends are set up twice currently, once
   * before QDB calc, then turned off in BandCalculationEngine() then 
   * back on. */
  if(cm->config_opts & CM_CONFIG_QDB) {
    if(cm->flags & CMH_QDB) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "ERROR in ConfigCM() CM already has QDBs\n");
    ConfigQDBAndW(cm, TRUE);
  }

  else if(always_calc_W) ConfigQDBAndW(cm, FALSE); /* FALSE says: don't calculate QDBs, W will still be calc'ed and set */
  
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
  esl_stopwatch_Start(w);  
  if(have_mother) { if((status = SubCMLogoddsify(cm, errbuf, mother_cm, mother_map)) != eslOK) return status; }
  else            CMLogoddsify(cm);
  esl_stopwatch_Stop(w); 
  FormatTimeString(time_buf, w->user, TRUE);
#if PRINTNOW
  fprintf(stdout, "\tCM logoddsify       %11s\n", time_buf);
#endif
  esl_stopwatch_Destroy(w);

  /*debug_print_cm_params(stdout, cm);
    debug_print_cp9_params(stdout, cm->cp9loc, TRUE);
    debug_print_cp9_params(stdout, cm->cp9glb, TRUE);*/
  return eslOK;
}

/* Function: ConfigLocal
 * Incept:   EPN, Tue Nov 29 13:46:29 2011 [updated] 
 *
 * Purpose:  Configure a CM for local alignment by spreading
 *           <p_internal_start> local entry probability evenly across
 *           all internal nodes, and by spreading <p_internal_exit>
 *           local exit probability evenly across all internal nodes.
 *
 *           Local entry probabilities for truncated alignment
 *           (cm->trbegin and cm->trbeginsc) are configured slightly
 *           differently, to allow inserts at the ends of the
 *           alignment. As with standard local begins,
 *           <p_internal_start> is spread evenly across all <nstarts>
 *           internal nodes, but entries into inserts are also
 *           allowed. The <p_internal_start>/<nstarts> probability
 *           for each node is divided between match and insert states
 *           of that node using a <psi> vector, psi[v] is the expected
 *           number of times state v is entered.
 *
 *           Note: to reproduce how Diana scored truncated alignments
 *           in the Kolbe and Eddy 2009 paper, set trbegin[v] to 
 *           (2 / (cm->clen * (cm->clen + 1))) for all v.
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
  float p;                      /* p_internal_start / nstarts */

  /* variables used for setting truncated begin probs (trbegin) */
  char  ***tmap   = NULL;       /* transition map, used to calc psi */
  int      i, j;                /* counters for transition maps */
  double  *psi    = NULL;       /* psi[v] = expected number of times state v entered (occupancy) */
  float    tmpvec[3];           /* temporary vector for normalizing trbegin */

  /* contract check */
  if(cm->flags & CMH_LOCAL_BEGIN)
    cm_Fail("ERROR in ConfigLocal(), CMH_LOCAL_BEGIN flag already up.\n");
  if(cm->flags & CMH_LOCAL_END)
    cm_Fail("ERROR in ConfigLocal(), CMH_LOCAL_END flag already up.\n");
  if(cm->flags & CM_IS_SUB)
    cm_Fail("ERROR in ConfigLocal(), CM is a sub CM, we can't localize it, we could if we first checked if its invalid (see cm_from_guide()).");

  had_scanmatrix = (cm->flags & CMH_SCANMATRIX) ? TRUE : FALSE;

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
  esl_vec_FSet(cm->begin,   cm->M, 0.);
  esl_vec_FSet(cm->trbegin, cm->M, 0.);

  /* Fill psi, used for setting trbegin */
  ESL_ALLOC(psi, sizeof(double) * cm->M);
  make_tmap(&tmap);
  fill_psi(cm, psi, tmap);

  /* Erase the previous transition p's from node 0. The only way out
   * of node 0 in standard scanners/aligners is going to be local
   * begin transitions from the root v=0 directly to MATP_MP, MATR_MR,
   * MATL_ML, and BIF_B states. In truncated scanners/aligners we also
   * allow transitions into inserts, see below. 
   *
   * First we want to save the transition probs we're about to zero,
   * so we can globalize this model later if we have to.
   * cm->root_trans is NULL prior to this, so if we call
   * ConfigGlobal() w/o executing this code, we'll die (with an Easel
   * contract exception).
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

  /* for trbegin, we allow transitions into ROOT_IL and ROOT_IR */
  esl_vec_FSet(tmpvec, 3, 0.);
  assert(cm->nodemap[1] == 3);
  tmpvec[0] = psi[1]; /* ROOT_IL */
  tmpvec[1] = psi[2]; /* ROOT_IR */
  tmpvec[2] = psi[3]; /* main state of node 1 (MATP_MP, MATR_MR, MATL_ML or BIF_B) */
  esl_vec_FNorm(tmpvec, 3);
  cm->trbegin[1] = tmpvec[0] * (1.-p_internal_start);
  cm->trbegin[2] = tmpvec[1] * (1.-p_internal_start);
  cm->trbegin[3] = tmpvec[2] * (1.-p_internal_start);

  /* Remaining nodes share p_internal_start.
   */
  p = p_internal_start / (float) nstarts;
  for (nd = 2; nd < cm->nodes; nd++) {
    if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
    	cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BIF_nd)  
      cm->begin[cm->nodemap[nd]] = p;
  }
  /* set trbegin */
  esl_vec_FSet(tmpvec, 3, 0.);
  esl_vec_FNorm(tmpvec, 3);
  for (nd = 2; nd < cm->nodes; nd++) {
    v = cm->nodemap[nd];
    tmpvec[0] = psi[v];
    if (cm->ndtype[nd] == MATP_nd) { 
      tmpvec[1] = psi[v+4]; /* MATP_IL */
      tmpvec[2] = psi[v+5]; /* MATP_IR */
      esl_vec_FNorm(tmpvec, 3);
      cm->trbegin[v]   = tmpvec[0] * p;
      cm->trbegin[v+4] = tmpvec[1] * p;
      cm->trbegin[v+5] = tmpvec[2] * p;
    }
    else if(cm->ndtype[nd] == MATL_nd) { 
      tmpvec[1] = psi[v+2]; /* MATL_IL */
      esl_vec_FNorm(tmpvec, 2);
      cm->trbegin[v]   = tmpvec[0] * p;
      cm->trbegin[v+2] = tmpvec[1] * p;
    }
    else if(cm->ndtype[nd] == MATR_nd) { 
      tmpvec[1] = psi[v+2]; /* MATR_IR */
      esl_vec_FNorm(tmpvec, 2);
      cm->trbegin[v]   = tmpvec[0] * p;
      cm->trbegin[v+2] = tmpvec[1] * p;
    }
    else if(cm->ndtype[nd] == BIF_nd) { 
      cm->trbegin[v] = p;
    }      
  }

  cm->flags |= CMH_LOCAL_BEGIN;
  
  /* New local begin probs invalidate QDBs if:
   * (A) they were calc'ed in glocal mode OR
   * (B) p_internal_start is non-default (not what was used 
   *     by cmbuild to calc the QDBs used in the CM file) 
   * If either is true, we have to recalc QDBs:
   * 
   * Note: we do this before configuring local ends because
   *       doing to calc QDBs we would have to turn local ends off.
   */
  if(((cm->flags & CMH_QDB) && (cm->flags & CMH_QDB_GLOBAL)) ||
     (fabs(cm->pbegin - DEFAULT_PBEGIN) > eslSMALLX1)) { 
    free(cm->dmin);
    free(cm->dmax);
    cm->dmin = NULL;
    cm->dmax = NULL;
    cm->flags &= ~CMH_QDB;
    cm->flags &= ~CMH_QDB_GLOBAL;
    cm->flags &= ~CMH_QDB_LOCAL;
    ConfigQDBAndW(cm, TRUE); /* TRUE says: calc QDBs */
  }      
  /* ConfigQDBAndW should rebuild scan matrix, if it existed */
  if((! (cm->flags & CMH_SCANMATRIX)) && had_scanmatrix)
    cm_Fail("ConfigLocal(), CM had a scan matrix, but ConfigQDBAndW didn't rebuild it.");

  /*****************************************************************
   * Internal exit.
   *****************************************************************/
  ConfigLocalEnds(cm, p_internal_exit);

  /* new local probs invalidate bit scores */
  cm->flags &= ~CMH_BITS;
     
  CMLogoddsify(cm);

  /* free tmap */
  if(tmap != NULL) { 
    for(i = 0; i < UNIQUESTATES; i++) { 
      for(j = 0; j < NODETYPES; j++) { 
	free(tmap[i][j]);
      }
      free(tmap[i]);
    }
    free(tmap);
  }
  if(psi    != NULL) free(psi);
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
  int had_scanmatrix;           /* true if CM had a scan matrix when function was entered */

  had_scanmatrix = (cm->flags & CMH_SCANMATRIX) ? TRUE : FALSE;

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
  for (v = 0; v < cm->M; v++)  { cm->begin[v] = cm->trbegin[v] = 0.; }
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

  /* Recalc QDBs if they exist and were calc'ed in local mode */
  if((cm->flags & CMH_QDB) && (cm->flags & CMH_QDB_LOCAL)) {
    free(cm->dmin);
    free(cm->dmax);
    cm->dmin = NULL;
    cm->dmax = NULL;
    cm->flags &= ~CMH_QDB;
    cm->flags &= ~CMH_QDB_LOCAL;
    cm->flags &= ~CMH_QDB_GLOBAL;
    ConfigQDBAndW(cm, TRUE); /* TRUE says: calc QDBs */
    /* ConfigQDBAndW should rebuild scan matrix, if it existed */
    if((! (cm->flags & CMH_SCANMATRIX)) && had_scanmatrix)
      cm_Fail("ConfigLocal(), CM had a scan matrix, but ConfigQDBAndW didn't rebuild it.");
  }
  /* new probs invalidate log odds scores */
  cm->flags &= ~CMH_BITS;

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
  if((cm->flags |= CMH_CP9LOC) && (cm->cp9loc->flags |= CPLAN9_EL))
    CPlan9NoEL(cm);

  esl_vec_FSet(cm->end, cm->M, 0.);

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
  if(cm->flags & CM_IS_SUB)
    cm_Fail("ERROR in ConfigLocalEnds(), CM is a sub CM, we can't localize it, we could if we first checked if its invalid (see cm_from_guide()).");

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

/* Function: ConfigQDBAndW
 * Date:     EPN, Wed Feb 13 17:44:02 2008
 *
 * Purpose:  Configure a CM's query dependent bands (QDBs) and/or 
 *           window length (W).
 * 
 * Args:     CM           - the covariance model
 *           do_calc_qdb  - TRUE to calculate QDBs and set cm->dmin, cm->dmax
 *                          using cm->beta_qdb
 */
int
ConfigQDBAndW(CM_t *cm, int do_calc_qdb)
{
  int mode;
  int v;
  int safe_windowlen;
  int *dmin, *dmax;
  int local_begins_on = (cm->flags & CMH_LOCAL_BEGIN) ? TRUE : FALSE;

  /* Three possible modes, depending on input args. 
   * We'll have to do the band calculation either:
   *
   * 1. one time with beta == cm->beta_W to calculate cm->W.
   * 2. one time with beta == cm->beta_qdb to calculate cm->dmin and
   *    cm->dmax bands, and implicitly W.
   * 3. twice, once with beta == cm->beta_W, to calculate cm->W,
   *    and again with beta == cm->beta_qdb to calculate cm->dmin
   *    cm->dmax bands. 
   */

  if(!do_calc_qdb) mode = 1;
  else { /* do_calc_qdb */
    if((cm->beta_W - cm->beta_qdb) > eslSMALLX1) { 
      mode = 3;
      /* TRUE if cm->beta_W used to calc W is greater than cm->beta_qdb,
       * in this case, we're in mode 3, cm->W will be less than
       * cm->dmax[0] (and cm->dmax[v] for other v as well). 
       * That's okay, we'll truncate those bands on d to 
       * never exceed W, but they'll be wider for some v than
       * they would have if we used cm->beta_W to calc qdbs.
       * (Imagine a CM with node 1 == BIF_nd, the BEGL and BEGR
       * subtrees can have wide bands, they just can't both 
       * combine to have a BIF subtree that exceeds cm->W, this
       * allows insertions in BEGL or BEGR subtrees, but not both.
       */
    }
    else mode = 2; 
    /* only calculate bands with cm->beta_qdb, then set W
     * as cm->dmax[0]. This may give W less than we would've
     * got with cm->beta_W. But if we're using the QDBs anyway,
     * the biggest hit we'll possibly get is cm->dmax[0] residues.
     */
  }

  /* run band calculation(s) */
  if(mode == 1 || mode == 3) { /* calculate cm->W */
    safe_windowlen = cm->clen * 3;
    while(!(BandCalculationEngine(cm, safe_windowlen, cm->beta_W, FALSE, &(dmin), &(dmax), NULL, NULL))) { 
      free(dmin);
      free(dmax);
      safe_windowlen *= 2;
      if(safe_windowlen > (cm->clen * 1000)) cm_Fail("ConfigQDBAndW(), mode 2, safe_windowlen big: %d\n", safe_windowlen);
    }
    cm->W = dmax[0];
    free(dmin);
    free(dmax);
  }
  if(mode == 2 || mode == 3) { /* calculate QDBs */
    safe_windowlen = cm->clen * 3;
    /* Contract check */
    if(cm->flags & CMH_QDB) cm_Fail("ConfigQDBAndW(): about to calculate QDBs, but CMH_QDB flag is already up.\n");
    if(cm->dmin != NULL) { free(cm->dmin); cm->dmin = NULL; }
    if(cm->dmax != NULL) { free(cm->dmax); cm->dmax = NULL; } 
    while(!(BandCalculationEngine(cm, safe_windowlen, cm->beta_qdb, FALSE, &(cm->dmin), &(cm->dmax), NULL, NULL))) { 
      free(cm->dmin); cm->dmin = NULL;
      free(cm->dmax); cm->dmax = NULL;
      safe_windowlen *= 2;
      if(safe_windowlen > (cm->clen * 1000)) cm_Fail("ConfigQDBAndW(), mode 2, safe_windowlen big: %d\n", safe_windowlen);
    }
    if(mode == 2) { /* set W as dmax[0], we're wasting time otherwise, looking at
		     * hits that are bigger than we're allowing with QDB. */
      cm->W = cm->dmax[0];
    } /* else, mode == 3, we set cm->W in loop above, it will be less than dmax[0] */

    else { /* mode == 3 */
      /* Quick check to make sure that cm->W >= dmin[v] for all v. If it's not, something went wrong */
      for(v = 0; v < cm->M; v++) { 
	if(cm->W < cm->dmin[v]) cm_Fail("ConfigQDBAndW(), mode 3, cm->W set as %d with beta: %g, but dmin[v:%d] (%d) exceeds it. QDBs calc'ed with beta: %g. This shouldn't happen.\n", cm->W, cm->beta_W, v, dmin[v], cm->beta_qdb); }
    }
    cm->flags |= CMH_QDB; /* raise the QDB flag */
    if(local_begins_on) { cm->flags |= CMH_QDB_LOCAL;  cm->flags &= ~CMH_QDB_GLOBAL; }
    else                { cm->flags |= CMH_QDB_GLOBAL; cm->flags &= ~CMH_QDB_LOCAL;  }
  }
  /* free and rebuild scan matrix to correspond to new QDBs and/or W, 
   * if it exists, this is where QDBs are potentially truncated 
   * in mode 3, that is: for all v, dmax[v] is reassigned as 
   * min(cm->dmax[v], cm->W) (note: we've ensured that 
   * dmin[v] < cm->W for all v above) 
   */
  if(cm->flags & CMH_SCANMATRIX) {
    int do_float = cm->smx->flags & cmSMX_HAS_FLOAT;
    int do_int   = cm->smx->flags & cmSMX_HAS_INT;
    cm_FreeScanMatrixForCM(cm);
    cm_CreateScanMatrixForCM(cm, do_float, do_int);
  }
  /* set cm->mlp7->max_length as cm->W */
  if(cm->mlp7 != NULL) cm->mlp7->max_length = cm->W;

  /*if(mode == 1 || mode == 3) printf("TEMP leaving ConfigQDBAndW(), mode: %d, set cm->W as:     %d with beta_W:   %g\n", mode, cm->W, cm->beta_W);
    if(mode == 2)              printf("TEMP leaving ConfigQDBAndW(), mode: %d, set cm->W as:     %d with beta_W:   %g\n", mode, cm->W, cm->beta_qdb);
    if(mode == 2 || mode == 3) printf("TEMP leaving ConfigQDBAndW(), mode: %d, set qdbs dmax[0]: %d with beta_qdb: %g\n", mode, cm->dmax[0], cm->beta_qdb);*/

  CMLogoddsify(cm); /* QDB calculation invalidates log odds scores */
  return eslOK;
}


#if 0
/* Function: ConfigQDB
 * Date:     EPN, Thu May  3 14:37:09 2007
 * Purpose:  Configure a CM's query dependent bands (QDBs).
 * Args:
 *           CM           - the covariance model
 */
int
ConfigQDB(CM_t *cm)
{
  int safe_windowlen;
  int local_begins_on = (cm->flags & CMH_LOCAL_BEGIN) ? TRUE : FALSE;

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
  if(local_begins_on) cm->flags |= CMH_QDB_LOCAL;
  else                cm->flags |= CMH_QDB_GLOBAL;

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
#endif

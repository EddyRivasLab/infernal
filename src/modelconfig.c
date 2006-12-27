/* modelconfig.c
 * SRE, Wed May  8 14:30:38 2002 [St. Louis]
 * SVN $Id$
 * 
 * Configuring a model into different global or local search modes.
 * 
 ******************************************************************
 * @LICENSE@
 ******************************************************************
 */

#include "squid.h"
#include "vectorops.h"

#include "structs.h"
#include "funcs.h"

void
ConfigLocal(CM_t *cm, float p_internal_start, float p_internal_exit)
{
  int v;			/* counter over states */
  int nd;			/* counter over nodes */
  int nstarts;			/* number of possible internal starts */

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
   */
  for (v = 0; v < cm->cnum[0]; v++)  cm->t[0][v] = 0.;

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
  cm->flags |= CM_LOCAL_BEGIN;
  
  /*****************************************************************
   * Internal exit.
   *****************************************************************/
  ConfigLocalEnds(cm, p_internal_exit);
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
  for (v = 0; v < cm->M; v++) cm->end[v] = 0.;
  /* Now, renormalize transitions */
  for (nd = 1; nd < cm->nodes; nd++) {
    if ((cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
	 cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BEGL_nd ||
	 cm->ndtype[nd] == BEGR_nd) && 
	cm->ndtype[nd+1] != END_nd)
      {
	v = cm->nodemap[nd];
	FNorm(cm->t[v], cm->cnum[v]);
      }
  }
  cm->flags &= ~CM_LOCAL_END; /* turn off local ends flag */
  return;
}

/**************************************************************************
 * EPN 10.02.06
 * Function: ConfigLocalEnds()
 *
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
	denom = FSum(cm->t[v], cm->cnum[v]);
	denom += cm->end[v];
	FScale(cm->t[v], cm->cnum[v], 1./denom);
      }
  }
  cm->flags |= CM_LOCAL_END;
  return;
}

void
ConfigLocal_fullsub(CM_t *cm, float p_internal_start, 
		    float p_internal_exit, int sstruct_nd,
		    int estruct_nd)
{
  int v;			/* counter over states */
  int nd;			/* counter over nodes */
  int nstarts;			/* number of possible internal starts */

  printf("in ConfigLocal_fullsub(), sstruct_nd: %d | estruct_nd: %d\n", sstruct_nd, estruct_nd);

  /* Currently, EL emissions in fullsub mode are disallowed.
   * To achieve, this set the EL self transition score to as close to IMPOSSIBLE 
   * as we can while still guaranteeing we won't get underflow errors.
   * we need cm->el_selfsc * W * 3 >= IMPOSSIBLE 
   * because we will potentially multiply cm->el_selfsc * W, and add that to 
   * 2 * IMPOSSIBLE, and IMPOSSIBLE must be > -FLT_MAX/3 so we can add it together 3 
   * times (see structs.h). 
   */
  ConfigLocal_DisallowELEmissions(cm);

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
   */
  for (v = 0; v < cm->cnum[0]; v++)  cm->t[0][v] = 0.;

  /* Node submap->sstruct gets prob 1-p_internal_start.
   */
  cm->begin[cm->nodemap[sstruct_nd]] = 1.-p_internal_start;
  printf("set cm->begin[%d]: %f\n", cm->nodemap[sstruct_nd], cm->begin[cm->nodemap[sstruct_nd]]);


  /* Remaining nodes share p_internal_start.
   */
  for (nd = 1; nd < cm->nodes; nd++) {
    if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
    	cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BIF_nd)  
      if(nd != sstruct_nd)
	cm->begin[cm->nodemap[nd]] = p_internal_start/(float)nstarts;
  }
  cm->flags |= CM_LOCAL_BEGIN;
  
  /*****************************************************************
   * Internal exit.
   *****************************************************************/
  ConfigLocalEnds(cm, p_internal_exit);
  return;
}

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
  cm->el_selfsc = (IMPOSSIBLE / (cm->W+1));
  return;
}

 /**************************************************************
  * Function: ConfigLocal_fullsub_post()
  * EPN, Mon Nov 13 13:32:27 2006
  * 
  * Purpose:  Configure a CM for local alignment in fullsub mode 
  *           using posterior probabilites from a CP9 HMM posterior
  *           decode of a sequence. Originally written for 'fullsub' 
  *           sub CM construction. Allow local begins into sub CM
  *           nodes that emit left before or at submap->sstruct,
  *           and....
  *
  * Args:     sub_cm      - the sub cm built from orig_cm
  *           orig_cm     - the original cm
  *           orig_cp9map - map from orig CM to orig CP9 HMM and vice versa
  *           submap      - the map from orig CM to sub CM and vice versa
  *           post        - posterior matrix, already filled
  *           L           - length of sequence to align
  * Returns:  
  * VOID
  */

 void 
 ConfigLocal_fullsub_post(CM_t *sub_cm, CM_t *orig_cm, CP9Map_t *orig_cp9map, CMSubMap_t *submap,
			  struct cp9_dpmatrix_s *post, int L)
 {
   int k;                       /* counter over HMM nodes */
   int v;			/* counter over states */
   int nd;			/* counter over nodes */
   int nstarts;			/* number of possible internal starts */
   int nexits;			/* number of possible internal ends */
   float denom;
   float sum_beg, sum_end;
   int lpos, rpos, lins, rins;
   int orig_v1, orig_v2;
   int orig_il, orig_ir;
   int orig_nd;
   CMEmitMap_t *sub_emap;           /* consensus emit map for the sub CM */
   int orig_v;

   char **nodetypes;
   char **sttypes;

   nodetypes = malloc(sizeof(char *) * 8);
   nodetypes[0] = "BIF";
   nodetypes[1] = "MATP";
   nodetypes[2] = "MATL";
   nodetypes[3] = "MATR";
   nodetypes[4] = "BEGL";
   nodetypes[5] = "BEGR";
   nodetypes[6] = "ROOT";
   nodetypes[7] = "END";

   sttypes = malloc(sizeof(char *) * 10);
   sttypes[0] = "D";
   sttypes[1] = "MP";
   sttypes[2] = "ML";
   sttypes[3] = "MR";
   sttypes[4] = "IL";
   sttypes[5] = "IR";
   sttypes[6] = "S";
   sttypes[7] = "E";
   sttypes[8] = "B";
   sttypes[9] = "EL";

  /* Currently, EL emissions in fullsub mode are disallowed.
   * To achieve, this set the EL self transition score to as close to IMPOSSIBLE 
   * as we can while still guaranteeing we won't get underflow errors.
   * we need cm->el_selfsc * W >= IMPOSSIBLE 
   * because we will potentially multiply cm->el_selfsc * W, and add that to 
   * 2 * IMPOSSIBLE, and IMPOSSIBLE must be > -FLT_MAX/3 so we can add it together 3 
   * times (see structs.h). 
   */
   ConfigLocal_DisallowELEmissions(sub_cm);

   sum_beg= sum_end = 0.;

   /* Zero all begin probs */
   for (v = 0; v < sub_cm->M; v++)  sub_cm->begin[v] = 0.;
   sub_emap = CreateEmitMap(sub_cm);

   printf("sstruct: %d\n", submap->sstruct);
   /* Fill in local begins up to sstruct */
   for(nd = 1; nd < submap->sstruct; nd++)
     {
       if(sub_cm->ndtype[nd] != MATL_nd)
	 Die("ERROR in ConfigLocal_fullsub_post: sstruct: %d but sub_cm node %d not MATL\n", submap->sstruct, nd);
       if(sub_emap->lpos[nd] != nd)
	 Die("ERROR in ConfigLocal_fullsub_post: sstruct: %d but sub_cm node %d lpos is %d (should be %d)\n", submap->sstruct, sub_emap->lpos[nd], nd);

       v = sub_cm->nodemap[nd];
       sub_cm->begin[v]  = Score2Prob(post->mmx[1][nd], 1.);
       sub_cm->begin[v] += Score2Prob(post->imx[1][nd], 1.);
       sum_beg += sub_cm->begin[v];
     }       
   /* Now fill in local begin for posn sstruct */
   orig_nd = orig_cp9map->pos2nd[submap->sstruct];
   orig_v  = orig_cm->nodemap[orig_nd];
   v       = submap->o2s_smap[orig_v][0];
   if(sub_cm->stid[v] != MATP_MP && sub_cm->stid[v] != MATL_ML)
     Die("ERROR in ConfigLocal_fullsub_post: 1st state of orig_cm node that maps to sstruct does not map to a sub CM MATP_MP or MATL_ML\n");
   nd = submap->sstruct;
   sub_cm->begin[v]  = Score2Prob(post->mmx[1][submap->sstruct], 1.);
   sub_cm->begin[v] += Score2Prob(post->imx[1][submap->sstruct], 1.);

   sum_beg += sub_cm->begin[v];
   printf("sum beg: %f\n", sum_beg);

   FNorm(sub_cm->begin, sub_cm->M);
   sub_cm->flags |= CM_LOCAL_BEGIN;

   free(sttypes);
   free(nodetypes);
   return;
 }

void
ConfigLocalEnforce(CM_t *cm, float p_internal_start, float p_internal_exit,
		   int enf_start, int enf_end)
{
  int v;			/* counter over states */
  int nd;			/* counter over nodes */
  int nstarts;			/* number of possible internal starts */
  int enf_start_pos;            /* consensus left position node enf_start emits to */
  int enf_end_pos;              /* consensus left position node enf_end   emits to */
  int nexits;			/* number of possible internal ends */
  float denom;
  CMEmitMap_t *emap;           /* consensus emit map for the CM */

  /* We want every parse to go through the MATL stretch from enf_start
   * to enf_end. To enforce this we disallow local begin and ends that
   * would allow parses to miss these nodes. */
  for(nd = enf_start; nd <= enf_end; nd++)
    {
      if(cm->ndtype[nd] != MATL_nd)
	Die("ERROR, trying to enforce a non-MATL stretch (node: %d not MATL).\n", nd);
    }
  emap = CreateEmitMap(cm); /* diff from ConfigLocalEnds() */
  enf_start_pos = emap->lpos[enf_start];
  enf_end_pos   = emap->lpos[enf_end];

  /* The following code is copied from ConfigLocal() and ConfigLocalEnds()
   * with modification to disallow local begins before enf_start and 
   * disallow exits from between closest_start and enf_end. This 
   * implementation sets local entry to the first node as 1-p_internal_start,
   * and end from final node as 1-p_internal_exit. */

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
      if(emap->lpos[nd] <= enf_start_pos &&
	 emap->rpos[nd] >= enf_end_pos) /* diff from ConfigLocalEnds() */
	nstarts++;
  }

  /* Zero everything.
   */
  for (v = 0; v < cm->M; v++)  cm->begin[v] = 0.;

  /* Erase the previous transition p's from node 0. The only
   * way out of node 0 is going to be local begin transitions
   * from the root v=0 directly to MATP_MP, MATR_MR, MATL_ML,
   * and BIF_B states.
   */
  for (v = 0; v < cm->cnum[0]; v++)  cm->t[0][v] = 0.;

  /* Node 1 gets prob 1-p_internal_start.
   */
  cm->begin[cm->nodemap[1]] = 1.-p_internal_start;

  /* Remaining nodes share p_internal_start.
   */
  for (nd = 2; nd < cm->nodes; nd++) {
    if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
    	cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BIF_nd)  
      if(emap->lpos[nd] <= enf_start_pos &&
	 emap->rpos[nd] >= enf_end_pos) /* diff from ConfigLocalEnds() */
	{
	  /*printf("enabling local begin into nd: %d lpos: %d rpos: %d s: %d e: %d\n", nd, emap->lpos[nd], emap->rpos[nd], enf_start_pos, enf_end_pos);*/
	  cm->begin[cm->nodemap[nd]] = p_internal_start/(float)nstarts;
	}
      else
	;/*printf("NOT enabling local begin into nd: %d lpos: %d rpos: %d s: %d e: %d\n", nd, emap->lpos[nd], emap->rpos[nd], enf_start_pos, enf_end_pos);*/
	
  }
  cm->flags |= CM_LOCAL_BEGIN;

  /*****************************************************************
   * Internal exit.
   *****************************************************************/
  /* Count internal nodes MATP, MATL, MATR, BEGL, BEGR that aren't
   * adjacent to END nodes.
   */
  nexits = 0;
  for (nd = 1; nd < cm->nodes; nd++) {
    if ((cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
	 cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BEGL_nd ||
	 cm->ndtype[nd] == BEGR_nd) && 
	cm->ndtype[nd+1] != END_nd)
      if(emap->lpos[nd] >= enf_end_pos || 
	 emap->rpos[nd] <  enf_start_pos) /* diff from ConfigLocalEnds() */
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
	if(emap->lpos[nd] >= enf_end_pos || 
	   emap->rpos[nd] <  enf_start_pos) /* diff from ConfigLocalEnds() */
	  {
	    /*printf("enabling local end from nd: %d lpos: %d rpos: %d s: %d e: %d\n", nd, emap->lpos[nd], emap->rpos[nd], enf_start_pos, enf_end_pos);*/
	    v = cm->nodemap[nd];
	    cm->end[v] = p_internal_exit / (float) nexits;
	    /* renormalize the main model transition distribution */
	    denom = FSum(cm->t[v], cm->cnum[v]);
	    denom += cm->end[v];
	    FScale(cm->t[v], cm->cnum[v], 1./denom);
	  }
	else
	  {
	    ;/*printf("NOT enabling local end from nd: %d lpos: %d rpos: %d s: %d e: %d\n", nd, emap->lpos[nd], emap->rpos[nd], enf_start_pos, enf_end_pos);*/
	  }
      }
  }
  cm->flags |= CM_LOCAL_END;
  FreeEmitMap(emap);
  return;
}


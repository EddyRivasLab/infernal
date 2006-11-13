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
#include "hmmer_funcs.h"


void
ConfigLocal(CM_t *cm, float p_internal_start, float p_internal_exit)
{
  int v;			/* counter over states */
  int nd;			/* counter over nodes */
  int nstarts;			/* number of possible internal starts */
  int nexits;			/* number of possible internal ends */
  float denom;

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
  int nexits;			/* number of possible internal ends */
  float denom;

  printf("in ConfigLocal_fullsub(), sstruct_nd: %d | estruct_nd: %d\n", sstruct_nd, estruct_nd);

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


 /**************************************************************
  * Function: ConfigLocal_fullsub_post()
  * EPN, Mon Nov 13 13:32:27 2006
  * 
  * Purpose:  Configure a CM for local alignment in fullsub mode 
  *           using posterior probabilites from a CP9 HMM posterior
  *           decode of a sequence. Originally written for 'fullsub' 
  *           sub CM construction. 
  *
  * Args:     cm        - the covariance model
  *           cp9map    - map from CM to CP9 HMM and vice versa
  *           post      - posterior matrix, already filled
  *           L         - length of sequence to align
  * Returns:  
  * VOID
  */

 void 
 ConfigLocal_fullsub_post(CM_t *cm, CP9Map_t *cp9map, struct cp9_dpmatrix_s *post, int L)
 {
   int k;                       /* counter over HMM nodes */
   int v;			/* counter over states */
   int nd;			/* counter over nodes */
   int nstarts;			/* number of possible internal starts */
   int nexits;			/* number of possible internal ends */
   float denom;
   float sum_beg, sum_end;
   int lpos, rpos, lins, rins;
   sum_beg= sum_end = 0.;

   /* Zero all begin probs */
   for (v = 0; v < cm->M; v++)  cm->begin[v] = 0.;

   /* Go through each CM node, filling in begin and end probs using the CP9 posteriors */
   for(nd = 0; nd < cm->nodes; nd++)
     {
       lpos = rpos = lins = rins = -1;
       if(cm->ndtype[nd] == MATP_nd)
	 {
	   lpos = cp9map->nd2lpos[k]; /* HMM node whose match state emits same left  pos as MATP_MP */
	   rpos = cp9map->nd2rpos[k]; /* HMM node whose match state emits same right pos as MATP_MP */
	   lins = cp9map->cs2hn[(cm->nodemap[nd] + 4)][0]; 
	   /* HMM node whose insert state emits same left pos as MATP_IL */
	   rins = cp9map->cs2hn[(cm->nodemap[nd] + 5)][0]; 
	   /* HMM node whose insert state emits same right pos as MATP_IR */
	   printf("cp9map->cs2hn[%3d][0]: %d\n", (cm->nodemap[nd]+4), cp9map->cs2hn[(cm->nodemap[nd]+4)][0]);
	   printf("cp9map->cs2hn[%3d][1]: %d\n", (cm->nodemap[nd]+4), cp9map->cs2hn[(cm->nodemap[nd]+4)][1]);
	   assert(cp9map->cs2hn[(cm->nodemap[nd] + 4)][1] == -1);
	   assert(cp9map->cs2hn[(cm->nodemap[nd] + 5)][1] == -1);
	 }
       else if(cm->ndtype[nd] == MATL_nd)
	 {
	   lpos = cp9map->nd2lpos[k]; /* HMM node whose match state emits same pos as MATL_ML */
	   lins = cp9map->cs2hn[(cm->nodemap[nd] + 2)][0]; 
	   /* HMM node whose insert state emits same pos as MATL_IL */
	 }
       else if(cm->ndtype[nd] == MATR_nd)
	 {
	   rpos = cp9map->nd2lpos[k]; /* HMM node whose match state emits same pos as MATR_MR */
	   rins = cp9map->cs2hn[(cm->nodemap[nd] + 2)][0]; 
	   /* HMM node whose insert state emits same pos as MATR_IR */
	 }
       else
	 continue;
       /* Set begin and end probs as sums of posterior probabilities of
	* HMM match and insert states emitting position 1 and L respectively.
	* NOTE: The way I handle the contribution of HMM insert posteriors 
	*       to end probs is WRONG because local ends are only allowed
	*       from consensus match states (to guarantee D&C will work).
	*       so you can never get to an insert of the same node before 
	*       doing a local end from that node. This is kind-of okay b/c
	*       the EL state is like an insert state. But I think a more robust 
	*       solution exists, although it may be more complex.
	*/
       v = cm->nodemap[nd];
       if(lpos != -1 && lins != -1)
	 {
	   cm->begin[v] += Score2Prob(post->mmx[1][lpos], 1.);
	   cm->begin[v] += Score2Prob(post->imx[1][lins], 1.);
	   cm->end[v]   += Score2Prob(post->mmx[L][lpos], 1.);
	   cm->end[v]   += Score2Prob(post->imx[L][lins], 1.);
	 }
       if(rpos != -1 && rins != -1)
	 {
	   cm->begin[v] += Score2Prob(post->mmx[1][rpos], 1.);
	   cm->begin[v] += Score2Prob(post->imx[1][rins], 1.);
	   cm->end[v]   += Score2Prob(post->mmx[L][rpos], 1.);
	   cm->end[v]   += Score2Prob(post->imx[L][rins], 1.);
	 }
     }
   for (v = 0; v < cm->M; v++)
     {
       sum_beg += cm->begin[v];
       sum_end += cm->end[v];
     }
   printf("sum beg: %f\n", sum_beg);
   printf("sum end: %f\n", sum_end);
   
   exit(1);
   return;
 }

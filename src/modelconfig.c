/* modelconfig.c
 * SRE, Wed May  8 14:30:38 2002 [St. Louis]
 * CVS $Id$
 * 
 * Configuring a model into different global or local search modes.
 * 
 ******************************************************************
 * @LICENSE@
 ******************************************************************
 */

#include "squid.h"
#include "structs.h"
#include "funcs.h"


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
  /* Spread the entry probability across internal nodes.
   */
  for (v = 0; v < cm->M; v++) cm->begin[v] = 0.;
  cm->begin[cm->nodemap[1]] = 1.-p_internal_exit;  /*v=3, the first state.*/
  for (nd = 2; nd < cm->nodes; nd++) {
    if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
    	cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BIF_nd)  
      cm->begin[cm->nodemap[nd]] = p_internal_start/(float)nstarts;
  }
  /* Make the begin transition the only way out of the root.
   */
  for (v = 0; v < cm->cnum[0]; v++)
    cm->t[0][v] = 0.;

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

/* EPN 11.28.05
 * CP9_cm2wrhmm.c 
 * 
 * Functions to build a Weinberg/Ruzzo maximum likelihood HMM from a CM. 
 * Uses the "CM plan 9" HMM architecture from cplan9.h and cplan.c. 
 * These functions use ideas/equations from Zasha Weinberg's thesis 
 * (notably p.122-124), but a CM plan 9 HMM is not exactly a Weinberg-Ruzzo
 * maximum likelihood heuristic HMM (though its close). The most
 * important difference is that a CM plan 9 HMM sometimes has a single
 * insert state that maps to 2 CM insert states. Due to self loops on
 * all insert states, it's impossible for a single HMM insert state to 
 * model two CM insert states, but the code below approximates it 
 * pretty well (it gets as close as possible (I think)).
 *
 * CP9_cm2wrhmm() is the main function, it sets the probabilities
 * of the HMM and checks that the expected number of times each
 * HMM state is entered is within a given threshold (0.0001 by default)
 * of the expected number of times each corresponding CM state is
 * entered.
 */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "squid.h"
#include "sre_stack.h"
#include "dirichlet.h"

#include "funcs.h"
#include "hmmer_structs.h"
#include "structs.h"

#include "esl_vectorops.h"

#include "stopwatch.h"          /* squid's process timing module        */

#include "cplan9.h"

static float
cm2hmm_emit_prob(CM_t *cm, int x, int i, int k, int *node_cc_left, int *node_cc_right, int *cc_node_map);

static void
cm2hmm_special_trans_cp9(CM_t *cm, struct cplan9_s *hmm, double *psi, char ***tmap, int **cs2hn_map,
			 int **cs2hs_map, int ***hns2cs_map);

static void
cm2hmm_trans_probs_cp9(CM_t *cm, struct cplan9_s *hmm, int k, double *psi, char ***tmap, 
		       int *cc_node_map, int **cs2hn_map, int **cs2hs_map, int ***hns2cs_map);

static void
fill_phi_cp9(struct cplan9_s *hmm, double **phi);

static void
hmm_add_single_trans_cp9(CM_t *cm, struct cplan9_s *hmm, int a, int b, int k, int hmm_trans_idx, 
			 double *psi, char ***tmap, int **cs2hn_map, int **cs2hs_map, int ***hns2cs_map);

static float
cm_sum_subpaths_cp9(CM_t *cm, int start, int end, char ***tmap, int **cs2hn_map, int **cs2hs_map, 
		    int ***hns2cs_map, int k, double *psi);

static int
check_psi_vs_phi_cp9(CM_t *cm, double *psi, double **phi, int ***hns2cs_map, int hmm_M,
		     double threshold, int debug_level);

static int 
check_cm_adj_bp(CM_t *cm, int *cc_node_map, int hmm_M);


static void 
CP9_fake_tracebacks(char **aseq, int nseq, int alen, int *matassign, struct cp9trace_s ***ret_tr);

static void 
CP9TraceCount(struct cplan9_s *hmm, char *dsq, float wt, struct cp9trace_s *tr);

static char *
CP9Statetype(char st);

static int 
CP9_node_chi_squared(struct cplan9_s *ahmm, struct cplan9_s *shmm, int nd, float thresh, 
		     int dual_mapping_insert);

/**************************************************************************
 * EPN 03.12.06
 * CP9_cm2wrhmm()
 *
 * Purpose:  Given a CM, an HMM and maps from nodes of HMM to CM and vice versa
 *           rewrite the HMM as a Weinberg/Ruzzo HMM with probabilities that
 *           are as close as possible to the CM's.
 * 
 * Args:    
 * CM_t *cm          - the CM
 * cplan9_s *hmm     - counts form CM plan 9 HMM
 * int ncc           - number of consensus columns
 * int *node_cc_left - consensus column each node's left emission maps to
 *                     [0..(cm->nodes-1)], -1 if maps to no consensus column
 * int *node_cc_right- consensus column each node's right emission corresponds to
 *                     [0..(cm->nodes-1)], -1 if maps to no consensus column
 * int *cc_node_map  - node that each consensus column maps to (is modelled by)
 *                     [1..hmm_ncc]
 * int **cs2hn_map   - 2D CM state to HMM node map, 1st D - CM state index
 *                     2nd D - 0 or 1 (up to 2 matching HMM states), value: HMM node
 *                     that contains state that maps to CM state, -1 if none.
 * int **cs2hs_map   - 2D CM state to HMM node map, 1st D - CM state index
 *                     2nd D - 2 elements for up to 2 matching HMM states, 
 *                     value: HMM STATE (0(M), 1(I), 2(D) that maps to CM state, -1 if none.
 * 
 *                     For example: HMM node cs2hn_map[v][0], state cs2hs_map[v][0]
 *                                  maps to CM state v.
 * 
 * int ***hns2cs_map  - 3D HMM node-state to CM state map, 1st D - HMM node index, 2nd D - 
 *                      HMM state (0(M), 1(I), 2(D)), 3rd D - 2 elements for up to 
 *                      2 matching CM states, value: CM states that map, -1 if none.
 *              
 *                     For example: CM states hns2cs_map[k][0][0] and hns2cs_map[k][0][1]
 *                                  map to HMM node k's match state.
 * int debug_level   - [0..3] tells the function what level of debugging print
 *                     statements to print.
 * Returns: TRUE if CP9 is constructed, FALSE if we get some error. 
 *          Its also possible one of the functions called within this function
 *          will print an error statement and exit.
 */
int
CP9_cm2wrhmm(CM_t *cm, struct cplan9_s *hmm, int *node_cc_left, int *node_cc_right, 
	     int *cc_node_map, int **cs2hn_map, int **cs2hs_map, int ***hns2cs_map, 
	     int debug_level)
{
  int       k;                 /* counter of consensus columns (HMM nodes)*/
  int       i,j;
  double    *psi;              /* expected num times each state visited in CM */
  double   **phi;              /* expected num times each state visited in HMM*/
  char     ***tmap;
  int *ap;                     /* CM state(s) (1 or 2) that maps to HMM state in node k*/
  int k_state;                 /* 0, 1 or 2, state in hmm node k*/

  double psi_vs_phi_threshold; /* the threshold that mapping (potentially summed) psi and 
				* phi values are allowed to be different by, without throwing
				* an error.*/
  int ret_val;                 /* return value */

  ap = malloc(sizeof(int) * 2);

  if(debug_level > 1)
    {
      printf("-------------------------------------------------\n");
      printf("In cm2wrhmm_cp9()\n");
    }

  psi = malloc(sizeof(double) * cm->M);
  make_tmap(&tmap);
  fill_psi(cm, psi, tmap);
  
  ZeroCPlan9(hmm);
  CPlan9SetNullModel(hmm, cm->null, 1.0); /* set p1 = 1.0 which corresponds to the CM */
  /* Special case 1st insert state maps to state 1 in the CM */
  for(i = 0; i < MAXABET; i++)
    {
      hmm->ins[0][i] = cm->e[1][i];
    }
  for(k = 1; k <= hmm->M; k++)
    {      
      for(i = 0; i < MAXABET; i++)
	{
	  hmm->mat[k][i] = 0.0;
	  hmm->ins[k][i] = 0.0;
	}
      /* First, take care of the match state. */
      k_state = HMMMATCH;
      ap[0] = hns2cs_map[k][k_state][0];
      ap[1] = hns2cs_map[k][k_state][1];
      /* ap[0] is a CM state that maps to HMM node k's match state */
      /* ap[1] is potentially another CM state that maps to HMM node k's match state
         (ex. if node k maps to the left half of a MATP node), and potentially = -1
         if no other state maps to hmm node k's match state.*/
      /* psi[ap[0]] is the expected number of times cm state ap[0] is entered. */
      for(i = 0; i < MAXABET; i++)
	{
	  hmm->mat[k][i] += psi[ap[0]] * 
	    cm2hmm_emit_prob(cm, ap[0], i, k, node_cc_left, node_cc_right, cc_node_map);
	  if(ap[1] != -1)
	    hmm->mat[k][i] += psi[ap[1]] *
	      cm2hmm_emit_prob(cm, ap[1], i, k, node_cc_left, node_cc_right, cc_node_map);
	}
      
      /* Now, do the insert state. */
      k_state = HMMINSERT;
      ap[0] = hns2cs_map[k][k_state][0];
      ap[1] = hns2cs_map[k][k_state][1];
      /* ap[0] is the only CM state that maps to HMM node k's insert state */
      /* ap[1] should be -1 unless k = hmm->M. */
      /* psi[ap[0]] is the expected number of times cm state ap[0] is entered. */
      for(i = 0; i < MAXABET; i++)
	{
	  hmm->ins[k][i] += psi[ap[0]] *
	    cm2hmm_emit_prob(cm, ap[0], i, k, node_cc_left, node_cc_right, cc_node_map);
	  if(ap[1] != -1)
	    hmm->ins[k][i] += psi[ap[1]] *
	      cm2hmm_emit_prob(cm, ap[1], i, k, node_cc_left, node_cc_right, cc_node_map);
	}
    }
  
  /* Done with emissions, fill in transitions of HMM (significantly more complex) */

  /* Step 1. Fill 'special' transitions, those INTO node 1, the N->N and N->M_1 transitions,
   * as well as transitions OUT of node M.
   */
  cm2hmm_special_trans_cp9(cm, hmm, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);

  for(k = 1; k < hmm->M; k++)
    {
      cm2hmm_trans_probs_cp9(cm, hmm, k, psi, tmap, cc_node_map, cs2hn_map, cs2hs_map, hns2cs_map);
    }

  CPlan9Renormalize(hmm);
  /* Create and fill phi to check to make sure our HMM is "close enough" to our CM.
   * phi[k][0..2] is the expected number of times HMM node k state 0 (match), 1(insert),
   * or 2(delete) is entered. These should be *very close* (within 0.0001) to the psi 
   * values for the CM states that they map to (psi[v] is the expected number of times
   * state v is entered in the CM). The reason we can't demand closer precision than
   * 0.0001 is that some HMM insert states map to 2 CM insert states, and its impossible
   * to perfectly mirror two CM insert states with 1 HMM insert states (due to self-loop
   * issues).
   */
  phi = malloc(sizeof(double *) * (hmm->M+1));
  for(k = 0; k <= hmm->M; k++)
    {
      phi[k] = malloc(sizeof(double) * 3);
    }
  fill_phi_cp9(hmm, phi);

  /*debug_print_cp9_params(hmm);*/
  psi_vs_phi_threshold = 0.0001;
  if(check_cm_adj_bp(cm, cc_node_map, hmm->M))
    {
      psi_vs_phi_threshold = 0.01;
    /* if check_cm_adj_bp() returns TRUE, then the following rare (but possible) situation
     * is true. Two adjacent consensus columns are modelled by the same MATP node. In this
     * case it is difficult for the CM plan 9 architecture to mirror the insert state
     * that maps to both the MATP_IL and MATP_IR of this node. It is more difficult than
     * for any other possible CM topology situation, so we relax the threshold when
     * checking psi and phi.
     */
    }    
  ret_val = check_psi_vs_phi_cp9(cm, psi, phi, hns2cs_map, hmm->M, psi_vs_phi_threshold, debug_level);
  CP9Logoddsify(hmm);

  free(ap);
  for(k = 0; k <= hmm->M; k++)
    {
      free(phi[k]);
    }
  free(phi);
  free(psi);
  for(i = 0; i < UNIQUESTATES; i++)
    {
      for(j = 0; j < NODETYPES; j++)
	free(tmap[i][j]);
      free(tmap[i]);
    }
  free(tmap);
  return ret_val;
}

/**************************************************************************
 * EPN 12.02.05
 * fill_psi()
 *
 * Purpose:  Fill psi matrix. Psi[v] is the expected number of times
 *           state v is entered.
 * 
 * Args:    
 * CM_t *cm          - the CM
 * double *psi       - psi[v] is expected number of times v is entered
 * char ***tmap      - eases coding transition use, hard-coded
 * 
 * Returns: (void) 
 */
void
fill_psi(CM_t *cm, double *psi, char ***tmap)
{
  int v; /*first state in cm node n*/
  int y;
  int x;
  char tmap_val;
  int n;
  double summed_psi;
  int nstates;
  int is_insert;

  /*psi[v] is the 'expected number of times state v is entered'.*/
  for (v = 0; v <= cm->M-1; v++)
    {
      psi[v] = 0.;
      if(cm->sttype[v] == IL_st || cm->sttype[v] == IR_st)
	is_insert = 1;
      else
	is_insert = 0;

      if(cm->sttype[v] == S_st)
	{
	  /* no transitions into start states - they're necessarily
	   * visited in every parse.
	   */
	  psi[v] = 1.0;
	}
      else if(is_insert)
	{
	  for (y = cm->pnum[v]-1; y >= 1; y--)
	    {
	      x = cm->plast[v] - y;
	      /* x is a parent of v, we're adding contribution 
	       * of transition from x to v. */
	      tmap_val = tmap[(int) cm->stid[x]][(int) cm->ndtype[cm->ndidx[v]+is_insert]][(int) cm->stid[v]];
	      if(tmap_val == -1)
		{
		  printf("tmap ERROR 1\n");
		  printf("v: %d | pnum[v]: %d | plast[v]: %d | y: %d | x: %d | d1: %d | d2: %d | d3: %d\n", v, cm->pnum[v], cm->plast[v], y, x, cm->stid[(int) x], (cm->ndtype[(int) cm->ndidx[v]+is_insert]), cm->stid[(int) v]);
		  exit(1);
		}
	      /*printf("before: psi[%d]: %f\n", v, psi[v]);
		printf("x: %d | tmap_val: %d | cm->t[x][tmap_val] : %f\n", x, tmap_val, cm->t[x][tmap_val]);*/
	      psi[v] += psi[x] * cm->t[x][(int) tmap_val];
	      /*printf("after: psi[%d]: %f\n", v, psi[v]);*/
	    }
	  /*printf("added self loop contribution of %f\n", (psi[v] * cm->t[v][0] / (1-cm->t[v][0])));*/
	  psi[v] += psi[v] * (cm->t[v][0] / (1-cm->t[v][0])); /*the contribution of the self insertion loops*/
	  /*printf("SL after: psi[%d]: %f\n", v, psi[v]);*/
	}
      else
	{
	  for (y = cm->pnum[v]-1; y >= 0; y--)
	    /*ERROR If t[y][v] is invalid, should be some number dependent on type of state y is 
	     *and type of state v is. I need a transition map.*/
	    {
	      x = cm->plast[v] - y;
	      /* x is a parent of v, we're adding contribution 
	       * of transition from x to v. */
	      tmap_val = tmap[(int) cm->stid[x]][(int) cm->ndtype[cm->ndidx[v]]][(int) cm->stid[v]];
	      
	      if(tmap_val == -1)
	      {
		printf("tmap ERROR 2\n");
		printf("v: %d | pnum[v]: %d | plast[v]: %d | y: %d | x: %d | d1: %d | d2: %d | d3: %d\n", v, cm->pnum[v], cm->plast[v], y, x, cm->stid[x], cm->ndtype[cm->ndidx[x]+1], cm->stid[v]);
		exit(1);
	      }
	      /*printf("before: psi[%d]: %f\n", v, psi[v]);
		printf("x: %d | y: %d | tmap_val: %d | cm->t[x][tmap_val] : %f\n", x, y, tmap_val, cm->t[x][tmap_val]);
	      */
	      psi[v] += psi[x] * cm->t[x][(int) tmap_val];
	      /*printf("after: psi[%d]: %f\n", v, psi[v]);*/
	    }
	}
      /*printf("psi[%d]: %15f\n", v, psi[v]);*/
    }  
  /* Sanity check. For any node the sum of psi values over
   * all split set states should be 1.0. 
   */
  for(n = 0; n < cm->nodes; n++)
    {
      summed_psi = 0.;
      if(cm->ndtype[n] == ROOT_nd)
	nstates = 3;
      else if(cm->ndtype[n] == BEGL_nd)
	nstates = 1;
      else if(cm->ndtype[n] == BEGR_nd)
	nstates = 2;
      else if(cm->ndtype[n] == BIF_nd)
	nstates = 1;
      else if(cm->ndtype[n] == MATP_nd)
	nstates = 6;
      else if(cm->ndtype[n] == MATL_nd)
	nstates = 3;
      else if(cm->ndtype[n] == MATR_nd)
	nstates = 3;
      else if(cm->ndtype[n] == END_nd)
	nstates = 1;
      else 
	Die("ERROR: bogus node type: %d\n", n);
      for(v = cm->nodemap[n]; v < cm->nodemap[n] + nstates; v++)
	if(cm->sttype[v] != IL_st && cm->sttype[v] != IR_st)
	  summed_psi += psi[v];
      if((summed_psi < 0.999) || (summed_psi > 1.001))
	{
	  printf("ERROR: summed psi of split states in node %d not 1.0 but : %f\n", n, summed_psi);
	  exit(1);
	}
      /* printf("split summed psi[%d]: %f\n", n, summed_psi);*/
    }
}

/**************************************************************************
 * EPN 12.02.05
 * make_tmap()
 *
 * Purpose:  Make the predefined transition map which tells you
 *           the index of a given transition from any of the 74
 *           transition sets.
 * 
 * char ***tmap;  A 3D char array. 
 *               1st D: statetype of v
 *               2nd D: type of downstream node.
 *               3rd D: statetype of y, that we're transitioning to.
 *               value: the index of v->y in cm->t[v]
 * Returns: (void) 
 */
void
make_tmap(char ****ret_tmap)
{
  int i,j,k;
  char ***tmap;

  tmap = malloc(sizeof(char **) * UNIQUESTATES);
  for(i = 0; i < UNIQUESTATES; i++)
    {
      tmap[i] = malloc(sizeof(char *) * NODETYPES);
      for(j = 0; j < NODETYPES; j++)
	{
	  tmap[i][j] = malloc(sizeof(char) * UNIQUESTATES);
	  for(k = 0; k < UNIQUESTATES; k++)
	    {
	      tmap[i][j][k] = -1;
	    }
	}
    }

  /*following code block generated by: 
   *perl ~nawrocki/notebook/5_1128_hmnl_ml_hmm/scripts/gen_tmap.pl
   */
  tmap[ROOT_S][BIF_nd][ROOT_IL] = 0;
  tmap[ROOT_S][BIF_nd][ROOT_IR] = 1;
  tmap[ROOT_S][BIF_nd][BIF_B] = 2;

  tmap[ROOT_S][MATP_nd][ROOT_IL] = 0;
  tmap[ROOT_S][MATP_nd][ROOT_IR] = 1;
  tmap[ROOT_S][MATP_nd][MATP_MP] = 2;
  tmap[ROOT_S][MATP_nd][MATP_ML] = 3;
  tmap[ROOT_S][MATP_nd][MATP_MR] = 4;
  tmap[ROOT_S][MATP_nd][MATP_D] = 5;

  tmap[ROOT_S][MATL_nd][ROOT_IL] = 0;
  tmap[ROOT_S][MATL_nd][ROOT_IR] = 1;
  tmap[ROOT_S][MATL_nd][MATL_ML] = 2;
  tmap[ROOT_S][MATL_nd][MATL_D] = 3;

  tmap[ROOT_S][MATR_nd][ROOT_IL] = 0;
  tmap[ROOT_S][MATR_nd][ROOT_IR] = 1;
  tmap[ROOT_S][MATR_nd][MATR_MR] = 2;
  tmap[ROOT_S][MATR_nd][MATR_D] = 3;

  tmap[ROOT_IL][BIF_nd][ROOT_IL] = 0;
  tmap[ROOT_IL][BIF_nd][ROOT_IR] = 1;
  tmap[ROOT_IL][BIF_nd][BIF_B] = 2;

  tmap[ROOT_IL][MATP_nd][ROOT_IL] = 0;
  tmap[ROOT_IL][MATP_nd][ROOT_IR] = 1;
  tmap[ROOT_IL][MATP_nd][MATP_MP] = 2;
  tmap[ROOT_IL][MATP_nd][MATP_ML] = 3;
  tmap[ROOT_IL][MATP_nd][MATP_MR] = 4;
  tmap[ROOT_IL][MATP_nd][MATP_D] = 5;

  tmap[ROOT_IL][MATL_nd][ROOT_IL] = 0;
  tmap[ROOT_IL][MATL_nd][ROOT_IR] = 1;
  tmap[ROOT_IL][MATL_nd][MATL_ML] = 2;
  tmap[ROOT_IL][MATL_nd][MATL_D] = 3;

  tmap[ROOT_IL][MATR_nd][ROOT_IL] = 0;
  tmap[ROOT_IL][MATR_nd][ROOT_IR] = 1;
  tmap[ROOT_IL][MATR_nd][MATR_MR] = 2;
  tmap[ROOT_IL][MATR_nd][MATR_D] = 3;

  tmap[ROOT_IR][BIF_nd][ROOT_IR] = 0;
  tmap[ROOT_IR][BIF_nd][BIF_B] = 1;

  tmap[ROOT_IR][MATP_nd][ROOT_IR] = 0;
  tmap[ROOT_IR][MATP_nd][MATP_MP] = 1;
  tmap[ROOT_IR][MATP_nd][MATP_ML] = 2;
  tmap[ROOT_IR][MATP_nd][MATP_MR] = 3;
  tmap[ROOT_IR][MATP_nd][MATP_D] = 4;

  tmap[ROOT_IR][MATL_nd][ROOT_IR] = 0;
  tmap[ROOT_IR][MATL_nd][MATL_ML] = 1;
  tmap[ROOT_IR][MATL_nd][MATL_D] = 2;

  tmap[ROOT_IR][MATR_nd][ROOT_IR] = 0;
  tmap[ROOT_IR][MATR_nd][MATR_MR] = 1;
  tmap[ROOT_IR][MATR_nd][MATR_D] = 2;

  tmap[BEGL_S][BIF_nd][BIF_B] = 0;

  tmap[BEGL_S][MATP_nd][MATP_MP] = 0;
  tmap[BEGL_S][MATP_nd][MATP_ML] = 1;
  tmap[BEGL_S][MATP_nd][MATP_MR] = 2;
  tmap[BEGL_S][MATP_nd][MATP_D] = 3;

  tmap[BEGR_S][BIF_nd][BEGR_IL] = 0;
  tmap[BEGR_S][BIF_nd][BIF_B] = 1;

  tmap[BEGR_S][MATP_nd][BEGR_IL] = 0;
  tmap[BEGR_S][MATP_nd][MATP_MP] = 1;
  tmap[BEGR_S][MATP_nd][MATP_ML] = 2;
  tmap[BEGR_S][MATP_nd][MATP_MR] = 3;
  tmap[BEGR_S][MATP_nd][MATP_D] = 4;

  tmap[BEGR_S][MATL_nd][BEGR_IL] = 0;
  tmap[BEGR_S][MATL_nd][MATL_ML] = 1;
  tmap[BEGR_S][MATL_nd][MATL_D] = 2;

  tmap[BEGR_IL][BIF_nd][BEGR_IL] = 0;
  tmap[BEGR_IL][BIF_nd][BIF_B] = 1;

  tmap[BEGR_IL][MATP_nd][BEGR_IL] = 0;
  tmap[BEGR_IL][MATP_nd][MATP_MP] = 1;
  tmap[BEGR_IL][MATP_nd][MATP_ML] = 2;
  tmap[BEGR_IL][MATP_nd][MATP_MR] = 3;
  tmap[BEGR_IL][MATP_nd][MATP_D] = 4;

  tmap[BEGR_IL][MATL_nd][BEGR_IL] = 0;
  tmap[BEGR_IL][MATL_nd][MATL_ML] = 1;
  tmap[BEGR_IL][MATL_nd][MATL_D] = 2;

  tmap[MATP_MP][BIF_nd][MATP_IL] = 0;
  tmap[MATP_MP][BIF_nd][MATP_IR] = 1;
  tmap[MATP_MP][BIF_nd][BIF_B] = 2;

  tmap[MATP_MP][MATP_nd][MATP_IL] = 0;
  tmap[MATP_MP][MATP_nd][MATP_IR] = 1;
  tmap[MATP_MP][MATP_nd][MATP_MP] = 2;
  tmap[MATP_MP][MATP_nd][MATP_ML] = 3;
  tmap[MATP_MP][MATP_nd][MATP_MR] = 4;
  tmap[MATP_MP][MATP_nd][MATP_D] = 5;

  tmap[MATP_MP][MATL_nd][MATP_IL] = 0;
  tmap[MATP_MP][MATL_nd][MATP_IR] = 1;
  tmap[MATP_MP][MATL_nd][MATL_ML] = 2;
  tmap[MATP_MP][MATL_nd][MATL_D] = 3;

  tmap[MATP_MP][MATR_nd][MATP_IL] = 0;
  tmap[MATP_MP][MATR_nd][MATP_IR] = 1;
  tmap[MATP_MP][MATR_nd][MATR_MR] = 2;
  tmap[MATP_MP][MATR_nd][MATR_D] = 3;

  tmap[MATP_MP][END_nd][MATP_IL] = 0;
  tmap[MATP_MP][END_nd][MATP_IR] = 1;
  tmap[MATP_MP][END_nd][END_E] = 2;

  tmap[MATP_ML][BIF_nd][MATP_IL] = 0;
  tmap[MATP_ML][BIF_nd][MATP_IR] = 1;
  tmap[MATP_ML][BIF_nd][BIF_B] = 2;

  tmap[MATP_ML][MATP_nd][MATP_IL] = 0;
  tmap[MATP_ML][MATP_nd][MATP_IR] = 1;
  tmap[MATP_ML][MATP_nd][MATP_MP] = 2;
  tmap[MATP_ML][MATP_nd][MATP_ML] = 3;
  tmap[MATP_ML][MATP_nd][MATP_MR] = 4;
  tmap[MATP_ML][MATP_nd][MATP_D] = 5;

  tmap[MATP_ML][MATL_nd][MATP_IL] = 0;
  tmap[MATP_ML][MATL_nd][MATP_IR] = 1;
  tmap[MATP_ML][MATL_nd][MATL_ML] = 2;
  tmap[MATP_ML][MATL_nd][MATL_D] = 3;

  tmap[MATP_ML][MATR_nd][MATP_IL] = 0;
  tmap[MATP_ML][MATR_nd][MATP_IR] = 1;
  tmap[MATP_ML][MATR_nd][MATR_MR] = 2;
  tmap[MATP_ML][MATR_nd][MATR_D] = 3;

  tmap[MATP_ML][END_nd][MATP_IL] = 0;
  tmap[MATP_ML][END_nd][MATP_IR] = 1;
  tmap[MATP_ML][END_nd][END_E] = 2;

  tmap[MATP_MR][BIF_nd][MATP_IL] = 0;
  tmap[MATP_MR][BIF_nd][MATP_IR] = 1;
  tmap[MATP_MR][BIF_nd][BIF_B] = 2;

  tmap[MATP_MR][MATP_nd][MATP_IL] = 0;
  tmap[MATP_MR][MATP_nd][MATP_IR] = 1;
  tmap[MATP_MR][MATP_nd][MATP_MP] = 2;
  tmap[MATP_MR][MATP_nd][MATP_ML] = 3;
  tmap[MATP_MR][MATP_nd][MATP_MR] = 4;
  tmap[MATP_MR][MATP_nd][MATP_D] = 5;

  tmap[MATP_MR][MATL_nd][MATP_IL] = 0;
  tmap[MATP_MR][MATL_nd][MATP_IR] = 1;
  tmap[MATP_MR][MATL_nd][MATL_ML] = 2;
  tmap[MATP_MR][MATL_nd][MATL_D] = 3;

  tmap[MATP_MR][MATR_nd][MATP_IL] = 0;
  tmap[MATP_MR][MATR_nd][MATP_IR] = 1;
  tmap[MATP_MR][MATR_nd][MATR_MR] = 2;
  tmap[MATP_MR][MATR_nd][MATR_D] = 3;

  tmap[MATP_MR][END_nd][MATP_IL] = 0;
  tmap[MATP_MR][END_nd][MATP_IR] = 1;
  tmap[MATP_MR][END_nd][END_E] = 2;

  tmap[MATP_D][BIF_nd][MATP_IL] = 0;
  tmap[MATP_D][BIF_nd][MATP_IR] = 1;
  tmap[MATP_D][BIF_nd][BIF_B] = 2;

  tmap[MATP_D][MATP_nd][MATP_IL] = 0;
  tmap[MATP_D][MATP_nd][MATP_IR] = 1;
  tmap[MATP_D][MATP_nd][MATP_MP] = 2;
  tmap[MATP_D][MATP_nd][MATP_ML] = 3;
  tmap[MATP_D][MATP_nd][MATP_MR] = 4;
  tmap[MATP_D][MATP_nd][MATP_D] = 5;

  tmap[MATP_D][MATL_nd][MATP_IL] = 0;
  tmap[MATP_D][MATL_nd][MATP_IR] = 1;
  tmap[MATP_D][MATL_nd][MATL_ML] = 2;
  tmap[MATP_D][MATL_nd][MATL_D] = 3;

  tmap[MATP_D][MATR_nd][MATP_IL] = 0;
  tmap[MATP_D][MATR_nd][MATP_IR] = 1;
  tmap[MATP_D][MATR_nd][MATR_MR] = 2;
  tmap[MATP_D][MATR_nd][MATR_D] = 3;

  tmap[MATP_D][END_nd][MATP_IL] = 0;
  tmap[MATP_D][END_nd][MATP_IR] = 1;
  tmap[MATP_D][END_nd][END_E] = 2;

  tmap[MATP_IL][BIF_nd][MATP_IL] = 0;
  tmap[MATP_IL][BIF_nd][MATP_IR] = 1;
  tmap[MATP_IL][BIF_nd][BIF_B] = 2;

  tmap[MATP_IL][MATP_nd][MATP_IL] = 0;
  tmap[MATP_IL][MATP_nd][MATP_IR] = 1;
  tmap[MATP_IL][MATP_nd][MATP_MP] = 2;
  tmap[MATP_IL][MATP_nd][MATP_ML] = 3;
  tmap[MATP_IL][MATP_nd][MATP_MR] = 4;
  tmap[MATP_IL][MATP_nd][MATP_D] = 5;

  tmap[MATP_IL][MATL_nd][MATP_IL] = 0;
  tmap[MATP_IL][MATL_nd][MATP_IR] = 1;
  tmap[MATP_IL][MATL_nd][MATL_ML] = 2;
  tmap[MATP_IL][MATL_nd][MATL_D] = 3;

  tmap[MATP_IL][MATR_nd][MATP_IL] = 0;
  tmap[MATP_IL][MATR_nd][MATP_IR] = 1;
  tmap[MATP_IL][MATR_nd][MATR_MR] = 2;
  tmap[MATP_IL][MATR_nd][MATR_D] = 3;

  tmap[MATP_IL][END_nd][MATP_IL] = 0;
  tmap[MATP_IL][END_nd][MATP_IR] = 1;
  tmap[MATP_IL][END_nd][END_E] = 2;

  tmap[MATP_IR][BIF_nd][MATP_IR] = 0;
  tmap[MATP_IR][BIF_nd][BIF_B] = 1;

  tmap[MATP_IR][MATP_nd][MATP_IR] = 0;
  tmap[MATP_IR][MATP_nd][MATP_MP] = 1;
  tmap[MATP_IR][MATP_nd][MATP_ML] = 2;
  tmap[MATP_IR][MATP_nd][MATP_MR] = 3;
  tmap[MATP_IR][MATP_nd][MATP_D] = 4;

  tmap[MATP_IR][MATL_nd][MATP_IR] = 0;
  tmap[MATP_IR][MATL_nd][MATL_ML] = 1;
  tmap[MATP_IR][MATL_nd][MATL_D] = 2;

  tmap[MATP_IR][MATR_nd][MATP_IR] = 0;
  tmap[MATP_IR][MATR_nd][MATR_MR] = 1;
  tmap[MATP_IR][MATR_nd][MATR_D] = 2;

  tmap[MATP_IR][END_nd][MATP_IR] = 0;
  tmap[MATP_IR][END_nd][END_E] = 1;

  tmap[MATL_ML][BIF_nd][MATL_IL] = 0;
  tmap[MATL_ML][BIF_nd][BIF_B] = 1;

  tmap[MATL_ML][MATP_nd][MATL_IL] = 0;
  tmap[MATL_ML][MATP_nd][MATP_MP] = 1;
  tmap[MATL_ML][MATP_nd][MATP_ML] = 2;
  tmap[MATL_ML][MATP_nd][MATP_MR] = 3;
  tmap[MATL_ML][MATP_nd][MATP_D] = 4;

  tmap[MATL_ML][MATL_nd][MATL_IL] = 0;
  tmap[MATL_ML][MATL_nd][MATL_ML] = 1;
  tmap[MATL_ML][MATL_nd][MATL_D] = 2;

  tmap[MATL_ML][MATR_nd][MATL_IL] = 0;
  tmap[MATL_ML][MATR_nd][MATR_MR] = 1;
  tmap[MATL_ML][MATR_nd][MATR_D] = 2;

  tmap[MATL_ML][END_nd][MATL_IL] = 0;
  tmap[MATL_ML][END_nd][END_E] = 1;

  tmap[MATL_D][BIF_nd][MATL_IL] = 0;
  tmap[MATL_D][BIF_nd][BIF_B] = 1;

  tmap[MATL_D][MATP_nd][MATL_IL] = 0;
  tmap[MATL_D][MATP_nd][MATP_MP] = 1;
  tmap[MATL_D][MATP_nd][MATP_ML] = 2;
  tmap[MATL_D][MATP_nd][MATP_MR] = 3;
  tmap[MATL_D][MATP_nd][MATP_D] = 4;

  tmap[MATL_D][MATL_nd][MATL_IL] = 0;
  tmap[MATL_D][MATL_nd][MATL_ML] = 1;
  tmap[MATL_D][MATL_nd][MATL_D] = 2;

  tmap[MATL_D][MATR_nd][MATL_IL] = 0;
  tmap[MATL_D][MATR_nd][MATR_MR] = 1;
  tmap[MATL_D][MATR_nd][MATR_D] = 2;

  tmap[MATL_D][END_nd][MATL_IL] = 0;
  tmap[MATL_D][END_nd][END_E] = 1;

  tmap[MATL_IL][BIF_nd][MATL_IL] = 0;
  tmap[MATL_IL][BIF_nd][BIF_B] = 1;

  tmap[MATL_IL][MATP_nd][MATL_IL] = 0;
  tmap[MATL_IL][MATP_nd][MATP_MP] = 1;
  tmap[MATL_IL][MATP_nd][MATP_ML] = 2;
  tmap[MATL_IL][MATP_nd][MATP_MR] = 3;
  tmap[MATL_IL][MATP_nd][MATP_D] = 4;

  tmap[MATL_IL][MATL_nd][MATL_IL] = 0;
  tmap[MATL_IL][MATL_nd][MATL_ML] = 1;
  tmap[MATL_IL][MATL_nd][MATL_D] = 2;

  tmap[MATL_IL][MATR_nd][MATL_IL] = 0;
  tmap[MATL_IL][MATR_nd][MATR_MR] = 1;
  tmap[MATL_IL][MATR_nd][MATR_D] = 2;

  tmap[MATL_IL][END_nd][MATL_IL] = 0;
  tmap[MATL_IL][END_nd][END_E] = 1;

  tmap[MATR_MR][BIF_nd][MATR_IR] = 0;
  tmap[MATR_MR][BIF_nd][BIF_B] = 1;

  tmap[MATR_MR][MATP_nd][MATR_IR] = 0;
  tmap[MATR_MR][MATP_nd][MATP_MP] = 1;
  tmap[MATR_MR][MATP_nd][MATP_ML] = 2;
  tmap[MATR_MR][MATP_nd][MATP_MR] = 3;
  tmap[MATR_MR][MATP_nd][MATP_D] = 4;

  tmap[MATR_MR][MATR_nd][MATR_IR] = 0;
  tmap[MATR_MR][MATR_nd][MATR_MR] = 1;
  tmap[MATR_MR][MATR_nd][MATR_D] = 2;

  tmap[MATR_D][BIF_nd][MATR_IR] = 0;
  tmap[MATR_D][BIF_nd][BIF_B] = 1;

  tmap[MATR_D][MATP_nd][MATR_IR] = 0;
  tmap[MATR_D][MATP_nd][MATP_MP] = 1;
  tmap[MATR_D][MATP_nd][MATP_ML] = 2;
  tmap[MATR_D][MATP_nd][MATP_MR] = 3;
  tmap[MATR_D][MATP_nd][MATP_D] = 4;

  tmap[MATR_D][MATR_nd][MATR_IR] = 0;
  tmap[MATR_D][MATR_nd][MATR_MR] = 1;
  tmap[MATR_D][MATR_nd][MATR_D] = 2;

  tmap[MATR_IR][BIF_nd][MATR_IR] = 0;
  tmap[MATR_IR][BIF_nd][BIF_B] = 1;

  tmap[MATR_IR][MATP_nd][MATR_IR] = 0;
  tmap[MATR_IR][MATP_nd][MATP_MP] = 1;
  tmap[MATR_IR][MATP_nd][MATP_ML] = 2;
  tmap[MATR_IR][MATP_nd][MATP_MR] = 3;
  tmap[MATR_IR][MATP_nd][MATP_D] = 4;

  tmap[MATR_IR][MATR_nd][MATR_IR] = 0;
  tmap[MATR_IR][MATR_nd][MATR_MR] = 1;
  tmap[MATR_IR][MATR_nd][MATR_D] = 2;

  tmap[BIF_B][BEGL_nd][BEGL_S] = 0;

  tmap[BIF_B][BEGR_nd][BEGR_S] = 0;

  *ret_tmap = tmap;
  return;
}


/**************************************************************************
 * EPN 03.13.06
 * cm2hmm_emit_prob()
 *
 * Purpose:  For a specific CM state, determine the probability of emitting
 *           residue i of the 4-letter RNA alphabet {A,C,G,U}. 
 * 
 *
 * Args:    
 * CM_t *cm          - the CM
 * int x             - the CM state 
 * int i             - the residue index in A=0, C=1, G=2, U=3
 * int k             - HMM node CM state x maps to
 * int *node_cc_left - consensus column each node's left emission maps to
 *                     [0..(cm->nodes-1)], -1 if maps to no consensus column
 * int *node_cc_right- consensus column each node's right emission corresponds to
 *                     [0..(cm->nodes-1)], -1 if maps to no consensus column
 * int *cc_node_map  - node that each consensus column maps to (is modelled by)
 *                     [1..hmm_ncc]
 * Returns: (float) probability of emitting letter i from correct half of CM state x.
 */
static float
cm2hmm_emit_prob(CM_t *cm, int x, int i, int k, int *node_cc_left, int *node_cc_right, 
		 int *cc_node_map)
{
  float ret_eprob;
  int is_left;
  int j;

  if(node_cc_left[cc_node_map[k]] == k)
    is_left = TRUE;
  else 
    is_left = FALSE;
  
  ret_eprob = 0.;

  /* trivial for non MATP_MP */
  if(cm->stid[x] != MATP_MP)
    ret_eprob = cm->e[x][i];
  
  else if(cm->stid[x] == MATP_MP)
    {
      /* determine which of the 16 indices to use */
      if(is_left)
	for(j = (i*MAXABET); j < ((i+1)*MAXABET); j++)
	  ret_eprob += cm->e[x][j];
      else
	for(j = i; j < (MAXABET*MAXABET); j+=MAXABET)
	  ret_eprob += cm->e[x][j];
    }
  return ret_eprob;
}



/**************************************************************************
 * EPN 03.12.06
 * cm2hmm_special_trans_cp9
 *
 * Purpose:  Fill the special transition probabilities (those INTO HMM node 1)
 *           for a CM plan 9 HMM given a CM.
 * 
 * Args:    
 * CM_t *cm          - the CM
 * cplan9_s *hmm      - the HMM
 * int *node_cc_left - consensus column each node's left emission maps to
 *                     [0..(cm->nodes-1)], -1 if maps to no consensus column
 * int *node_cc_right- consensus column each node's right emission corresponds to
 *                     [0..(cm->nodes-1)], -1 if maps to no consensus column
 * int *cc_node_map  - node that each consensus column maps to (is modelled by)
 *                     [1..hmm_ncc]
 * double *psi       - psi[v] is the expected number of times state v is entered
 *                     in a CM parse
 * char ***tmap;      - the hard-coded transition map
 * int **cs2hn_map   - 2D CM state to HMM node map, 1st D - CM state index
 *                     2nd D - 0 or 1 (up to 2 matching HMM states), value: HMM node
 *                     that contains state that maps to CM state, -1 if none.
 * int **cs2hs_map   - 2D CM state to HMM node map, 1st D - CM state index
 *                     2nd D - 2 elements for up to 2 matching HMM states, 
 *                     value: HMM STATE (0(M), 1(I), 2(D) that maps to CM state, -1 if none.
 * 
 *                     For example: HMM node cs2hn_map[v][0], state cs2hs_map[v][0]
 *                                  maps to CM state v.
 * 
 * int ***hns2cs_map  - 3D HMM node-state to CM state map, 1st D - HMM node index, 2nd D - 
 *                      HMM state (0(M), 1(I), 2(D)), 3rd D - 2 elements for up to 
 *                      2 matching CM states, value: CM states that map, -1 if none.
 *              
 *                     For example: CM states hns2cs_map[k][0][0] and hns2cs_map[k][0][1]
 *                                  map to HMM node k's match state.
 * Returns: (void) 
 */
static void
cm2hmm_special_trans_cp9(CM_t *cm, struct cplan9_s *hmm, double *psi, char ***tmap, int **cs2hn_map,
			 int **cs2hs_map, int ***hns2cs_map)
{
  int *ap; /* CM states a' that map to HMM state a, 
	      * ap[1] is -1 if only 1 CM state maps to a*/
  int *bp; /* CM states b' that map to HMM state b, 
	      * bp[1] is -1 if only 1 CM state maps to b*/
  int k_state; /*either HMMMATCH, HMMINSERT, or HMMDELETE*/
  int k;      /* HMM node index */
  int hmm_trans_idx; /*0-8;  CTMM, CTMI, CTMD, CTIM, CTII, CTID, CTDM, CTDI, or CTDD*/
  float d;

  ap = malloc(sizeof (int) * 2);
  bp = malloc(sizeof (int) * 2);
  
  /* Fill all special transitions with virtual counts, later normalize these into 
   * probabilities. CM p9 Special transitions are: transitions into HMM node 1, and 
   * transitions out of HMM node hmm->M, and the N->N self transition.
   *
   * Strategy:
   *
   * for each HMM transition from special state B (M_0) (maps to CM's ROOT_S) and N (I_0)
   *       (maps to CM's ROOT_IL) to state b in HMM node 1:
   *     - determine CM state(s) b' (can be 1 or 2) that map to HMM state b.
   *       (we must be moving down the CM (i.e. state num (ROOT_S) = 0 and (ROOT_IL) = 1)
   *     - determine sum(B->b') or sum(N->b') the summed probability of all paths from state(s)
   *       B to b' or N to b'.
   * 
   *       set virtual counts of transitions as follows:
   *       hmm->t[0][CTMI] (B->N)   = sum(over state(s) b' that map to N) psi[0] * ROOT_S -> b' 
   *                                 (just ROOT_S->ROOT_IL)
   *       hmm->t[0][CTMD] (B->D_1) = sum(over state(s) b' that map to D_1) psi[0] * ROOT_S -> b'
   *       begin[1] (B->M_1) = sum(over state(s) b' that map to M_1) psi[0] * ROOT_S -> b'
   *       hmm->t[0][CTIM]   = sum(over state(s) b' that map to M_1) psi[1] * ROOT_IL -> b' 
   *       hmm->t[0][CTII]   = sum(over state(s) b' that map to N)   psi[1] * ROOT_IL -> b' 
   *                      (just ROOT_IL->ROOT_IL)
   *       hmm->t[0][CTID]   = sum(over state(s) b' that map to D_1) psi[1] * ROOT_IL -> b' 
   *
   *       note: psi[0] = 1.0         
   *  
   * fill transitions out of node hmm->M. 
   * for each HMM transition from state a in HMM node M to either state END_E or state b 
   *       in node M:
   *     - determine CM state(s) a' and b' (can be 1 or 2 for each a' and b') that map 
   *       to HMM states a and b respectively
   *
   *     - IF a' > b' (we're moving DOWN the CM)
   *       determine sum(a'->b') the summed probability of all paths from state(s)
   *       a' to b'.
   * 
   *       set virtual counts of HMM transition a->b as 
   *       sum(over states a' and b') psi[a'] * a' -> b'
   *
   *     - ELSE IF a' < b' (we're moving UP the CM)
   *       determine sum(b'->a') the summed probability of all paths from state(s)
   *       b' to a'
   * 
   *       set virtual counts of HMM transition a->b as 
   *       sum(over states a' and b') psi[b'] * b' -> a'
   * 
   * Where psi[x] is the expected number of times CM state x is entered, as calc'ed
   * in fill_psi. psi[0] = 1.0
   * 
   */

  /* Set up the probability of transitioning into the first HMM node's states. */

  /* Transition 1 CTMM*/
  /* Special case, transition prob into M_1 from B is actually 
   * hmm->begin[1]; first we calc what it should be, then 
   * we switch hmm->t[0][CTMM] and hmm->begin[1].
   */
  k = 0;
  k_state = HMMMATCH;
  ap[0] = hns2cs_map[k][k_state][0];
  ap[1] = hns2cs_map[k][k_state][1];

  k_state = HMMMATCH;
  bp[0] = hns2cs_map[k+1][k_state][0];
  bp[1] = hns2cs_map[k+1][k_state][1];
  /*printf("0 CTMM: k: %4d | n: %4d | ap[0]: %4d ap[1]: %4d | bp[0]: %4d bp[1]: %4d\n", k, n, ap[0], ap[1], bp[0], bp[1]);*/
  hmm_trans_idx = CTMM;
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  /*hmm_set_single_trans_cp9(cm, hmm, ap, bp, k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);*/
  /* switch 'em */
  hmm->begin[1] = hmm->t[0][CTMM];
  hmm->t[0][CTMM] = 0.;
  
  /* Transition 2: CTMI; B -> N
   */
  k = 0;
  k_state = HMMMATCH;
  ap[0] = hns2cs_map[k][k_state][0];
  ap[1] = hns2cs_map[k][k_state][1];

  k_state = HMMINSERT;
  bp[0] = hns2cs_map[k][k_state][0];
  bp[1] = hns2cs_map[k][k_state][1];
  hmm_trans_idx = CTMI;
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);

  /* Transition 3: CTMD; B -> D_1 */
  k = 0;
  k_state = HMMMATCH;
  ap[0] = hns2cs_map[k][k_state][0];
  ap[1] = hns2cs_map[k][k_state][1];

  k_state = HMMDELETE;
  bp[0] = hns2cs_map[k+1][k_state][0];
  bp[1] = hns2cs_map[k+1][k_state][1];
  hmm_trans_idx = CTMD;
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);

  /* Transition 4: CTIM; N -> M_1*/
  k = 0;
  k_state = HMMINSERT;
  ap[0] = hns2cs_map[k][k_state][0];
  ap[1] = hns2cs_map[k][k_state][1];

  k_state = HMMMATCH;
  bp[0] = hns2cs_map[k+1][k_state][0];
  bp[1] = hns2cs_map[k+1][k_state][1];
  hmm_trans_idx = CTIM;
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);

  /* Transition 5: CTII; N -> N */
  k = 0;
  k_state = HMMINSERT;
  ap[0] = hns2cs_map[k][k_state][0];
  ap[1] = hns2cs_map[k][k_state][1];

  k_state = HMMINSERT;
  bp[0] = hns2cs_map[k][k_state][0];
  bp[1] = hns2cs_map[k][k_state][1];
  hmm_trans_idx = CTII;
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);

  /* Transition 6: CTID; N -> D_1 */
  k = 0;
  k_state = HMMINSERT;
  ap[0] = hns2cs_map[k][k_state][0];
  ap[1] = hns2cs_map[k][k_state][1];

  k_state = HMMDELETE;
  bp[0] = hns2cs_map[k+1][k_state][0];
  bp[1] = hns2cs_map[k+1][k_state][1];
  hmm_trans_idx = CTID;
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);

  /* Transitions 7-9, CTDM, CTDI, CTDD, all 0.0, there's no D_0 state */
  hmm->t[0][CTDM] = 0.;
  hmm->t[0][CTDI] = 0.;
  hmm->t[0][CTDD] = 0.;

  /* Finally, normalize the transition probabilities
   * Not strictly necessary, a CPlan9Renormalize() call will do this*/

  d = FSum(hmm->begin+1, hmm->M) + hmm->t[0][CTMI] + hmm->t[0][CTMD];
  FScale(hmm->begin+1, hmm->M, 1./d);
  hmm->t[0][CTMI] /= d;
  hmm->t[0][CTMD] /= d;

  FNorm(hmm->t[0]+3, 3);	/* transitions out of insert for node 0 (state N)*/

  k = 0;
  /*
  printf("S hmm->t[%d][CTMM]: %f\n", k, hmm->t[k][CTMM]);
  printf("S hmm->begin[%3d]: %f\n", 1, hmm->begin[1]);
  printf("S hmm->t[%d][CTMI]: %f\n", k, hmm->t[k][CTMI]);
  printf("S hmm->t[%d][CTMD]: %f\n", k, hmm->t[k][CTMD]);
  printf("S hmm->t[%d][CTIM]: %f\n", k, hmm->t[k][CTIM]);
  printf("S hmm->t[%d][CTII]: %f\n", k, hmm->t[k][CTII]);
  printf("S hmm->t[%d][CTID]: %f\n", k, hmm->t[k][CTID]);
  printf("S hmm->t[%d][CTDM]: %f\n", k, hmm->t[k][CTDM]);
  printf("S hmm->t[%d][CTDI]: %f\n", k, hmm->t[k][CTDI]);
  printf("S hmm->t[%d][CTDD]: %f\n\n", k, hmm->t[k][CTDD]);  
  printf("\n");
  */
  /*********************************************************************/
  /*********************************************************************/
  /* Now handle transitions OUT of HMM node M */
  /* Transition 1: CTMM; this is M_M -> E
   * a = node M, match state.
   * b = E state
   */
  k = hmm->M;
  k_state = HMMMATCH;
  ap[0] = hns2cs_map[k][k_state][0];
  ap[1] = hns2cs_map[k][k_state][1];

  bp[0] = (cm->M) - 1;
  bp[1] = -1;
  hmm_trans_idx = CTMM;
  /* special case, this transition is hmm->end[hmm->M] NOT hmm->t[M][CTMM]. */
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  /* switch 'em */
  hmm->end[hmm->M] = hmm->t[hmm->M][CTMM];
  hmm->t[hmm->M][CTMM] = 0.;
  
  /* Transition 2: CTMI; this is M_M -> I_M
   * a = node M, match state
   * b = node M, insert state
   */
  k = hmm->M;
  k_state = HMMMATCH;
  ap[0] = hns2cs_map[k][k_state][0];
  ap[1] = hns2cs_map[k][k_state][1];

  k_state = HMMINSERT;
  bp[0] = hns2cs_map[k][k_state][0];
  bp[1] = hns2cs_map[k][k_state][1];
  /*printf("CTMI: k: %4d | n: %4d | ap[0]: %4d ap[1]: %4d | bp[0]: %4d bp[1]: %4d\n", k, cc_node_map[k], ap[0], ap[1], bp[0], bp[1]);*/
  hmm_trans_idx = CTMI;
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);

  /* Transition 3: CTMD, no corresponding transitiona
   * This is an illegal HMM transition.
   */
  hmm->t[hmm->M][CTMD] = 0.;

  /* Transition 4: CTIM; this is I_M -> E
   * a = node k, insert state.
   * b = E state
   */
  k = hmm->M;
  k_state = HMMINSERT;
  ap[0] = hns2cs_map[k][k_state][0];
  ap[1] = hns2cs_map[k][k_state][1];

  bp[0] = (cm->M) - 1;
  bp[1] = -1;
  hmm_trans_idx = CTIM;
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);

  /* Transition 5: CTII; this is I_M -> I_M
   * a = node M, insert state.
   * b = node M, insert state.
   */
  k = hmm->M;
  k_state = HMMINSERT;
  ap[0] = hns2cs_map[k][k_state][0];
  ap[1] = hns2cs_map[k][k_state][1];

  k_state = HMMINSERT;
  bp[0] = hns2cs_map[k][k_state][0];
  bp[1] = hns2cs_map[k][k_state][1];
  /*printf("CTII: k: %4d | n: %4d | ap[0]: %4d ap[1]: %4d | bp[0]: %4d bp[1]: %4d\n", k, cc_node_map[k], ap[0], ap[1], bp[0], bp[1]);*/
  hmm_trans_idx = CTII;
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);

  /* Transition 6: CTID, no corresponding transition.
   * This is an illegal HMM transition.
   */
  hmm->t[hmm->M][CTID] = 0.;

  /* Transition 7: CTDM; this is D_M -> E
   * a = node M, delete state.
   * b = E state
   */
  k = hmm->M;
  k_state = HMMDELETE;
  ap[0] = hns2cs_map[k][k_state][0];
  ap[1] = hns2cs_map[k][k_state][1];

  bp[0] = (cm->M) - 1;
  bp[1] = -1;
  hmm_trans_idx = CTDM;
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);

  /* Transition 8: CTDI - this is D_M -> I_M
   * a = node M, delete state.
   * b = node M, insert state.
   */
  k = hmm->M;
  k_state = HMMDELETE;
  ap[0] = hns2cs_map[k][k_state][0];
  ap[1] = hns2cs_map[k][k_state][1];

  k_state = HMMINSERT;
  bp[0] = hns2cs_map[k][k_state][0];
  bp[1] = hns2cs_map[k][k_state][1];
  /*printf("CTDI: k: %4d | n: %4d | ap[0]: %4d ap[1]: %4d | bp[0]: %4d bp[1]: %4d\n", k, cc_node_map[k], ap[0], ap[1], bp[0], bp[1]);*/
  hmm_trans_idx = CTDI;
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);

  /* Transition 9: CTDD, no corresponding transition.
   * This is an illegal HMM transition.
   */
  hmm->t[hmm->M][CTDD] = 0.;

  /* Finally, normalize the transition probabilities
   * Not strictly necessary, a CPlan9Renormalize() call will do this*/
  k = hmm->M;
  d = FSum(hmm->t[k], 3) + hmm->end[k]; 
  FScale(hmm->t[k], 3, 1./d);
  hmm->end[k] /= d;
  FNorm(hmm->t[k]+3, 3);
  FNorm(hmm->t[k]+6, 3);

  /*
  printf("S hmm->t[%d][CTMM]: %f\n", k, hmm->t[k][CTMM]);
  printf("S hmm->end[%3d]   : %f\n", k, hmm->end[k]);
  printf("S hmm->t[%d][CTMI]: %f\n", k, hmm->t[k][CTMI]);
  printf("S hmm->t[%d][CTMD]: %f\n", k, hmm->t[k][CTMD]);
  printf("S hmm->t[%d][CTIM]: %f\n", k, hmm->t[k][CTIM]);
  printf("S hmm->t[%d][CTII]: %f\n", k, hmm->t[k][CTII]);
  printf("S hmm->t[%d][CTID]: %f\n", k, hmm->t[k][CTID]);
  printf("S hmm->t[%d][CTDM]: %f\n", k, hmm->t[k][CTDM]);
  printf("S hmm->t[%d][CTDI]: %f\n", k, hmm->t[k][CTDI]);
  printf("S hmm->t[%d][CTDD]: %f\n\n", k, hmm->t[k][CTDD]);  
  printf("\n");
  */

  free(ap);
  free(bp);
  return;
}

/**************************************************************************
 * EPN 03.12.06
 * cm2hmm_trans_probs_cp9()
 *
 * Purpose:  Fill transition "virtual counts" from HMM node k to 
 *           the following HMM node (k+1) given a CM, then normalize to
 *           probabilities.
 *  
 * Reference: Zasha Weinberg thesis p.123-124         
 *
 * Args:    
 * CM_t *cm          - the CM
 * cplan9_s *hmm     - the CM plan 9 HMM
 * int k             - the HMM node we're filling transitions for
 * double *psi       - psi[v] is the expected number of times state v is entered
 *                     in a CM parse
 * char ***tmap;     - the hard-coded transition map
 * int *cc_node_map  - node that each consensus column maps to (is modelled by)
 *                     [1..hmm_ncc]
 * int **cs2hn_map   - 2D CM state to HMM node map, 1st D - CM state index
 *                     2nd D - 0 or 1 (up to 2 matching HMM states), value: HMM node
 *                     that contains state that maps to CM state, -1 if none.
 * int **cs2hs_map   - 2D CM state to HMM node map, 1st D - CM state index
 *                     2nd D - 2 elements for up to 2 matching HMM states, 
 *                     value: HMM STATE (0(M), 1(I), 2(D) that maps to CM state, -1 if none.
 * 
 *                     For example: HMM node cs2hn_map[v][0], state cs2hs_map[v][0]
 *                                  maps to CM state v.
 * 
 * int ***hns2cs_map  - 3D HMM node-state to CM state map, 1st D - HMM node index, 2nd D - 
 *                      HMM state (0(M), 1(I), 2(D)), 3rd D - 2 elements for up to 
 *                      2 matching CM states, value: CM states that map, -1 if none.
 *              
 *                     For example: CM states hns2cs_map[k][0][0] and hns2cs_map[k][0][1]
 *                                  map to HMM node k's match state.
 * Returns: (void)    
 */
static void
cm2hmm_trans_probs_cp9(CM_t *cm, struct cplan9_s *hmm, int k, double *psi, char ***tmap, 
		       int *cc_node_map, int **cs2hn_map, int **cs2hs_map, int ***hns2cs_map)
{
  int *ap; /* CM states a' that map to HMM state a, 
	      * ap[1] is -1 if only 1 CM state maps to a*/
  int *bp; /* CM states b' that map to HMM state b, 
	      * bp[1] is -1 if only 1 CM state maps to b*/
  
  int n;   /* CM node that maps to HMM node k */
  int k_state; /*either HMMMATCH, HMMINSERT, or HMMDELETE*/

  int hmm_trans_idx; /*0-8;  CTMM, CTMI, CTMD, CTIM, CTII, CTID, CTDM, CTDI, or CTDD*/
  float d;

  ap = malloc(sizeof (int) * 2);
  bp = malloc(sizeof (int) * 2);
  n = cc_node_map[k];
  /* Fill all 9 transitions with virtual counts, later normalize these into 
   * probabilities.
   *
   * Strategy:
   *
   * for each HMM transition from state a in node k to state b in either node k or k+1:
   *     - determine CM state(s) a' and b' (can be 1 or 2 for each a' and b') that map 
   *       to HMM states a and b respectively
   *
   *     - IF a' > b' (we're moving DOWN the CM)
   *       determine sum(a'->b') the summed probability of all paths from state(s)
   *       a' to b'.
   * 
   *       set virtual counts of HMM transition a->b as 
   *       sum(over states a' and b') psi[a'] * a' -> b'
   *
   *     - ELSE IF a' < b' (we're moving UP the CM)
   *       determine sum(b'->a') the summed probability of all paths from state(s)
   *       b' to a'
   * 
   *       set virtual counts of HMM transition a->b as 
   *       sum(over states a' and b') psi[b'] * b' -> a'
   * 
   * Where psi[x] is the expected number of times CM state x is entered, as calc'ed
   * in fill_psi.
   * 
   */

  /* Transition 1: CTMM
   * a = node k, match state.
   * b = node k+1, match state.
   */
  k_state = HMMMATCH;
  ap[0] = hns2cs_map[k][k_state][0];
  ap[1] = hns2cs_map[k][k_state][1];

  k_state = HMMMATCH;
  bp[0] = hns2cs_map[k+1][k_state][0];
  bp[1] = hns2cs_map[k+1][k_state][1];
  hmm_trans_idx = CTMM;
  ////printf("CTMM: k: %4d | n: %4d | ap[0]: %4d ap[1]: %4d | bp[0]: %4d bp[1]: %4d\n", k, cc_node_map[k], ap[0], ap[1], bp[0], bp[1]);
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  /*hmm_set_single_trans_cp9(cm, hmm, ap, bp, k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);*/

  /* Transition 2: CTMI
   * a = node k, match state.
   * b = node k, insert state.
   */
  k_state = HMMMATCH;
  ap[0] = hns2cs_map[k][k_state][0];
  ap[1] = hns2cs_map[k][k_state][1];

  k_state = HMMINSERT;
  bp[0] = hns2cs_map[k][k_state][0];
  bp[1] = hns2cs_map[k][k_state][1];
  ////printf("CTMI: k: %4d | n: %4d | ap[0]: %4d ap[1]: %4d | bp[0]: %4d bp[1]: %4d\n", k, cc_node_map[k], ap[0], ap[1], bp[0], bp[1]);
  hmm_trans_idx = CTMI;
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);

  /* Transition 3: CTMD
   * a = node k, match state.
   * b = node k+1, delete state.
   */
  k_state = HMMMATCH;
  ap[0] = hns2cs_map[k][k_state][0];
  ap[1] = hns2cs_map[k][k_state][1];

  k_state = HMMDELETE;
  bp[0] = hns2cs_map[k+1][k_state][0];
  bp[1] = hns2cs_map[k+1][k_state][1];
  ////printf("CTMD: k: %4d | n: %4d | ap[0]: %4d ap[1]: %4d | bp[0]: %4d bp[1]: %4d\n", k, cc_node_map[k], ap[0], ap[1], bp[0], bp[1]);
  hmm_trans_idx = TMD;
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);

  /* Transition 4: CTIM
   * a = node k, insert state.
   * b = node k+1, match state.
   */
  k_state = HMMINSERT;
  ap[0] = hns2cs_map[k][k_state][0];
  ap[1] = hns2cs_map[k][k_state][1];

  k_state = HMMMATCH;
  bp[0] = hns2cs_map[k+1][k_state][0];
  bp[1] = hns2cs_map[k+1][k_state][1];
  ////printf("CTIM: k: %4d | n: %4d | ap[0]: %4d ap[1]: %4d | bp[0]: %4d bp[1]: %4d\n", k, cc_node_map[k], ap[0], ap[1], bp[0], bp[1]);
  hmm_trans_idx = CTIM;
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);

  /* Transition 5: CTII
   * a = node k, insert state.
   * b = node k, insert state.
   */
  k_state = HMMINSERT;
  ap[0] = hns2cs_map[k][k_state][0];
  ap[1] = hns2cs_map[k][k_state][1];

  k_state = HMMINSERT;
  bp[0] = hns2cs_map[k][k_state][0];
  bp[1] = hns2cs_map[k][k_state][1];
  ////printf("CTII: k: %4d | n: %4d | ap[0]: %4d ap[1]: %4d | bp[0]: %4d bp[1]: %4d\n", k, cc_node_map[k], ap[0], ap[1], bp[0], bp[1]);

  hmm_trans_idx = CTII;
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);

  /* Transition 6: CTID - a CM plan 9 transition not in plan 7
   * a = node k, insert state.
   * b = node k+1, delete state.
   */
  k_state = HMMINSERT;
  ap[0] = hns2cs_map[k][k_state][0];
  ap[1] = hns2cs_map[k][k_state][1];

  k_state = HMMDELETE;
  bp[0] = hns2cs_map[k+1][k_state][0];
  bp[1] = hns2cs_map[k+1][k_state][1];
  ////printf("CTID: k: %4d | n: %4d | ap[0]: %4d ap[1]: %4d | bp[0]: %4d bp[1]: %4d\n", k, cc_node_map[k], ap[0], ap[1], bp[0], bp[1]);
  hmm_trans_idx = CTID;
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);

  /* Transition 7: CTDM
   * a = node k, delete state.
   * b = node k+1, match state.
   */
  k_state = HMMDELETE;
  ap[0] = hns2cs_map[k][k_state][0];
  ap[1] = hns2cs_map[k][k_state][1];

  k_state = HMMMATCH;
  bp[0] = hns2cs_map[k+1][k_state][0];
  bp[1] = hns2cs_map[k+1][k_state][1];
  ////printf("CTDM: k: %4d | n: %4d | ap[0]: %4d ap[1]: %4d | bp[0]: %4d bp[1]: %4d\n", k, cc_node_map[k], ap[0], ap[1], bp[0], bp[1]);
  hmm_trans_idx = CTDM;
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);

  /* Transition 8: CTDI - a CM plan 9 transition not in plan 9 
   * a = node k, delete state.
   * b = node k, insert state.
   */
  k_state = HMMDELETE;
  ap[0] = hns2cs_map[k][k_state][0];
  ap[1] = hns2cs_map[k][k_state][1];

  k_state = HMMINSERT;
  bp[0] = hns2cs_map[k][k_state][0];
  bp[1] = hns2cs_map[k][k_state][1];
  ////printf("CTDI: k: %4d | n: %4d | ap[0]: %4d ap[1]: %4d | bp[0]: %4d bp[1]: %4d\n", k, cc_node_map[k], ap[0], ap[1], bp[0], bp[1]);
  hmm_trans_idx = CTDI;
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);

  /* Transition 9: CTDD
   * a = node k, delete state.
   * b = node k+1, delete state.
   */
  k_state = HMMDELETE;
  ap[0] = hns2cs_map[k][k_state][0];
  ap[1] = hns2cs_map[k][k_state][1];

  k_state = HMMDELETE;
  bp[0] = hns2cs_map[k+1][k_state][0];
  bp[1] = hns2cs_map[k+1][k_state][1];
  ////printf("CTDD: k: %4d | n: %4d | ap[0]: %4d ap[1]: %4d | bp[0]: %4d bp[1]: %4d\n", k, cc_node_map[k], ap[0], ap[1], bp[0], bp[1]);
  hmm_trans_idx = CTDD;
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[0], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[0], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);
  hmm_add_single_trans_cp9(cm, hmm, ap[1], bp[1], k, hmm_trans_idx, psi, tmap, cs2hn_map, cs2hs_map, hns2cs_map);

  /* Finally, normalize the transition probabilities
   * Not strictly necessary, a CPlan9Renormalize() call will do this*/
  d = FSum(hmm->t[k], 3) + hmm->end[k]; 
  FScale(hmm->t[k], 3, 1./d);
  hmm->end[k] /= d;

  FNorm(hmm->t[k]+3, 3);
  FNorm(hmm->t[k]+6, 3);
  /* Printing transition probs for HMM HERE!!*/
  /*
    printf("hmm->t[%d][CTMM]: %f\n", k, hmm->t[k][CTMM]);
    printf("hmm->t[%d][CTMI]: %f\n", k, hmm->t[k][CTMI]);
    printf("hmm->t[%d][CTMD]: %f\n", k, hmm->t[k][CTMD]);
    printf("hmm->t[%d][CTIM]: %f\n", k, hmm->t[k][CTIM]);
    printf("hmm->t[%d][CTII]: %f\n", k, hmm->t[k][CTII]);
    printf("hmm->t[%d][CTID]: %f\n", k, hmm->t[k][CTID]);
    printf("hmm->t[%d][CTDM]: %f\n", k, hmm->t[k][CTDM]);
    printf("hmm->t[%d][CTDI]: %f\n", k, hmm->t[k][CTDI]);
    printf("hmm->t[%d][CTDD]: %f\n\n", k, hmm->t[k][CTDD]);  
    printf("\n");
  */
  free(ap);
  free(bp);
  return;
}

/**************************************************************************
 * EPN 03.14.06
 * fill_phi_cp9()
 *
 * Purpose:  Fill phi matrix for a CM plan 9 HMM. 
 *           phi[k][v] is the expected number of times 
 *           HMM node k's state v (either 0 (match), 1 (insert) or
 *           (2) delete) is entered in the HMM.
 * 
 * Args:    
 * cplan9_s *hmm      - the HMM
 * double **phi       - phi array, phi[k][v] is expected number of times
 *                     state v (0 = match, 1 insert, 2 = delete) in 
 *                     node k is visited. Node 0 is special, 
 *                     state 0 = B state, state 1 = N_state, state 2 = NULL
 * 
 * Returns: (void) 
 */
static void
fill_phi_cp9(struct cplan9_s *hmm, double **phi)
{
  int k;

  /* Take care of special states, case by case. */
  phi[0][0] =  1.0; /*B state, necessary start*/
  phi[0][1] =  phi[0][0] * hmm->t[0][CTMI];
  phi[0][1] += phi[0][1] * (hmm->t[0][CTII] / (1.-hmm->t[0][CTII])); /* self insert contribution */
  phi[0][2] = 0.;
  /*
  printf("phi[0][0]: %f\n", phi[0][0]);
  printf("phi[0][1]: %f\n", phi[0][1]);
  printf("phi[0][2]: %f\n", phi[0][2]);
  */

  /* Handle all other nodes (including M) */
  for (k = 1; k <= hmm->M; k++)
    {
      /* match could've come from k-1 match, k-1 insert or k-1 delete */
      phi[k][0] =  phi[k-1][0] * hmm->t[k-1][CTMM];
      phi[k][0] += phi[k-1][2] * hmm->t[k-1][CTDM];
      phi[k][0] += phi[k-1][1] * hmm->t[k-1][CTIM];
      phi[k][0] += hmm->begin[k];

      /* again, we have to do deletes prior to inserts */
      /* deletes could've come from k-1 match, k-1 delete, k-1 insert */
      phi[k][2] =  phi[k-1][0] * hmm->t[k-1][CTMD];
      phi[k][2] += phi[k-1][1] * hmm->t[k-1][CTID];
      phi[k][2] += phi[k-1][2] * hmm->t[k-1][CTDD];
 
      /* inserts could've come from k match, k delete, or k insert */
      phi[k][1] =  phi[k][0] * hmm->t[k][CTMI];
      phi[k][1] += phi[k][2] * hmm->t[k][CTDI];
      /* self loops are special */
      phi[k][1] += (phi[k][1] * (hmm->t[k][CTII] / (1-hmm->t[k][CTII])));

      /*printf("phi[%d][0]: %f\n", k, phi[k][0]);
	printf("phi[%d][1]: %f\n", k, phi[k][1]);
	printf("phi[%d][2]: %f\n", k, phi[k][2]);
      */
    }
}


/**************************************************************************
 * EPN 03.15.06
 * Function: hmm_add_single_trans_cp9()
 *
 * Purpose:  Add a virtual counts contribution to a single CM plan 9 HMM transition. 
 *  
 * Args:    
 * CM_t *cm          - the CM
 * cplan9_s *hmm     - the HMM
 * int a             - a CM state that maps to HMM state we're transitioning out of
 * int b             - a CM state that maps to HMM state we're transitioning into
 * int k             - the HMM node we're setting a single transition for
 * int hmm_trans_idx - 0-8, the HMM transition index we're setting
 * double *psi       - psi[v] is the expected number of times state v is entered
 *                     in a CM parse
 * char ***tmap;      - the hard-coded transition map
 * int **cs2hn_map   - 2D CM state to HMM node map, 1st D - CM state index
 *                     2nd D - 0 or 1 (up to 2 matching HMM states), value: HMM node
 *                     that contains state that maps to CM state, -1 if none.
 * int **cs2hs_map   - 2D CM state to HMM node map, 1st D - CM state index
 *                     2nd D - 2 elements for up to 2 matching HMM states, 
 *                     value: HMM STATE (0(M), 1(I), 2(D) that maps to CM state, -1 if none.
 * 
 *                     For example: HMM node cs2hn_map[v][0], state cs2hs_map[v][0]
 *                                  maps to CM state v.
 * 
 * int ***hns2cs_map  - 3D HMM node-state to CM state map, 1st D - HMM node index, 2nd D - 
 *                      HMM state (0(M), 1(I), 2(D)), 3rd D - 2 elements for up to 
 *                      2 matching CM states, value: CM states that map, -1 if none.
 *              
 *                     For example: CM states hns2cs_map[k][0][0] and hns2cs_map[k][0][1]
 *                                  map to HMM node k's match state.
 * Returns: (void) 
 */
static void
hmm_add_single_trans_cp9(CM_t *cm, struct cplan9_s *hmm, int a, int b, int k, int hmm_trans_idx, 
			double *psi, char ***tmap, int **cs2hn_map, int **cs2hs_map, int ***hns2cs_map)
{
  /* check if we've got real CM state ids */
  if(a == -1 || b == -1)
    return;

  if(a <= b) /* going DOWN the CM */
    hmm->t[k][hmm_trans_idx] += psi[a] * cm_sum_subpaths_cp9(cm, a, b, tmap, cs2hn_map, cs2hs_map, hns2cs_map, k, psi);
  else if (a > b) /* going UP the CM */
    hmm->t[k][hmm_trans_idx] += psi[b] * cm_sum_subpaths_cp9(cm, b, a, tmap, cs2hn_map, cs2hs_map, hns2cs_map, k, psi);
}

/**************************************************************************
 * EPN 02.24.06
 * cm_sum_subpaths_cp9()
 *
 * Purpose:  Calculated probability of getting from one state (start) to 
 *           another (end) in a CM, taking special considerations. 
 * 
 *           Sum the probability of all subpaths that start 
 *           at "start" and end at "end" (ignoring "end"->"end" and
 *           "start" -> "start" transitions if they exist)
 * 
 *           This function is used to help determine CM plan 9 HMM
 *           transition probabilities. If we're trying to set a particular
 *           transition (1 of 9) out of HMM node k, we ignore the contribution 
 *           of subparses that correspond to other transitions out of node k.
 *           For example, we don't want to include the probability of an 
 *           insert(node a) ->match (node a+1) sub parse when calculating
 *           the transition probability for match(node a) -> match (node a+1), 
 *           (CTMM) because CTIM maps to that transition. 
 *          
 * Args:    
 * CM_t *cm          - the CM
 * int start         - state index we're starting at
 * int end           - state index we're ending in
 * int **cs2hn_map   - 2D CM state to HMM node map, 1st D - CM state index
 *                     2nd D - 0 or 1 (up to 2 matching HMM states), value: HMM node
 *                     that contains state that maps to CM state, -1 if none.
 * int **cs2hs_map   - 2D CM state to HMM node map, 1st D - CM state index
 *                     2nd D - 2 elements for up to 2 matching HMM states, 
 *                     value: HMM STATE (0(M), 1(I), 2(D) that maps to CM state, -1 if none.
 * 
 *                     For example: HMM node cs2hn_map[v][0], state cs2hs_map[v][0]
 *                                  maps to CM state v.
 * 
 * int ***hns2cs_map  - 3D HMM node-state to CM state map, 1st D - HMM node index, 2nd D - 
 *                      HMM state (0(M), 1(I), 2(D)), 3rd D - 2 elements for up to 
 *                      2 matching CM states, value: CM states that map, -1 if none.
 *              
 *                     For example: CM states hns2cs_map[k][0][0] and hns2cs_map[k][0][1]
 *                                  map to HMM node k's match state.
 * int k             - HMM node we're calc'ing transition (out of node k) for
 * double *psi       - psi[v] is the expected number of times state v is entered
 *                     in a CM parse
 * Returns: Float, the summed probability of all subpaths through the CM
 *          starting at "start" and ending at "end".
 */
static float
cm_sum_subpaths_cp9(CM_t *cm, int start, int end, char ***tmap, int **cs2hn_map, 
		 int **cs2hs_map, int ***hns2cs_map, int k, double *psi)
{
  int s_n; /* CM node that maps to HMM node with start state */
  int e_n; /* CM node that maps to HMM node with end state */
  int v; /* state index in CM */
  
  double *sub_psi; /*sub_psi[v] is the expected number of times state v is
		    * entered given we started at state "start", this is
		    * the summed probability of all paths starting at "start"
		    * and ending at v.
		    */
  float to_return;
  int y;
  int x;
  char tmap_val;
  int n_v; /* CM node containing state v*/
  int is_insert; /* 1 if v is insert, 0 if not */
  float insert_to_start; /* is start is not an insert and the insert state of 
			  * HMM node k < start, this is the
			  * contribution of insert -> start (which should be ignored
			  * because the TMI and TDI transitions solely map to
			  * this transition probability.
			  */

  /*printf("\nin cm_sum_subpaths_cp9, start: %d | end: %d\n", start, end);*/
  if(start > end)
    {
      printf("ERROR in cm_sum_subpaths_cp9: start: %d > end: %d\n", start, end);
      exit(1);
    }
  if(start == end)
    {
      if(cm->sttype[start] != IL_st && cm->sttype[start] != IR_st)
	{
	  /* This is possible, though unlikely. It only happens (I think and hope) if
	   * we've got two adjacent (k and k+1) columns modelled by the same MATP node n, 
	   * so when setting a transition from node k ->  node k+1, we enter this function
	   * three times with start == end. Once with each MATP_MP, MATP_ML, MATP_MR, and MATP_D
	   * because these states map to the beginning and ending states of the M_k->M_k+1, 
	   * M_k->D_k+1, D_k->M_k+1, and D_k->D_k+1 transitions respectively. 
	   * The solution is to return 1.0, 
	   * so the contribution is simply psi[MATP_M*] (in the function 
	   * hmm_add_single_trans_cp9 that called this function).
	   */
	  if((cm->stid[start] != MATP_MP && cm->stid[start] != MATP_D) && 
	     (cm->stid[start] != MATP_ML && cm->stid[start] != MATP_MR))
	    {
	      printf("ERROR asking for self transition of non-insert, non-MATP state: %d\n", start);
	      exit(1);
	    }
	  return 1.0;
	}
      return cm->t[start][0];
    }
  to_return = 0.;
  s_n = cm->ndidx[start];
  e_n = cm->ndidx[end];
  
  sub_psi = malloc(sizeof(double) * (end - start + 1));
  /* Initialize sub_psi[0]. Need to check if we need to ignore the probability
   * mass from the CM insert state(s) that maps to the HMM insert state of this node 
   * (these insert states are hns2cs_map[k][1][0] and (potentially) hns2cs_map[k][1][1]) 
   * that goes through "start" (which must map to either the M or 
   * D state of this HMM node).
   */
  sub_psi[0] = 1.; /* have to start in "start" */
  
  if((cm->sttype[start] != IL_st && cm->sttype[start] != IR_st) &&
     (cm->sttype[end] != IL_st && cm->sttype[end] != IR_st))
    {
      insert_to_start = 0.;
      if(hns2cs_map[k][1][0] < start) 
	insert_to_start = psi[hns2cs_map[k][1][0]] * cm_sum_subpaths_cp9(cm, hns2cs_map[k][1][0], start, tmap, cs2hn_map, cs2hs_map, hns2cs_map, k, psi);
      if((hns2cs_map[k][1][1] != -1) && (hns2cs_map[k][1][1] < start))
	insert_to_start += psi[hns2cs_map[k][1][1]] * cm_sum_subpaths_cp9(cm, hns2cs_map[k][1][1], start, tmap, cs2hn_map, cs2hs_map, hns2cs_map, k, psi);
      sub_psi[0] -= insert_to_start / psi[start];
    }
  /* note: when cm_sum_subpaths_cp9 is called recursively above
   * it will never result in another recursive call, 
   * because its "start" is an insert state.  
   */
  
  for (v = (start+1); v <= end; v++) 
    {
      sub_psi[v-start] = 0.;
      n_v = cm->ndidx[v];
      if(cm->sttype[v] == IL_st || cm->sttype[v] == IR_st)
	is_insert = 1;
      else
	is_insert = 0;
      
      if(cm->sttype[v] == S_st)
	{
	  /* previous state is necessarily either a BIF_B or a END_E, either
	   * way, there's no transitions FROM previous state to this state, so
	   * we handle this in a special way.*/
	  sub_psi[v-start] = sub_psi[(v-1)-start] * 1.;
	}
      if((v != end && is_insert) && ((cs2hn_map[v][0] == k) || cs2hn_map[v][1] == k))
	{
	  /*skip the contribution*/
	  /*printf("v: %d | skipping the contribution\n", v);
	    printf("\tcs2hn_map[%d][0] : %d\n", v, cs2hn_map[v][0]);
	    printf("\tcs2hn_map[%d][1] : %d\n", v, cs2hn_map[v][1]);
	    printf("\tcs2hn_map[start][0] : %d\n", v, cs2hn_map[start][0]);
	    printf("\tcs2hn_map[start][1] : %d\n", v, cs2hn_map[start][1]);
	  */
	}
      else 
	{
	  for (y = cm->pnum[v]-1; y >= is_insert; y--) 
	    {
	      x = cm->plast[v] - y;
	      /* x is a parent of v, we're adding contribution 
	       * of transition from x to v. */
	      tmap_val = tmap[(int) cm->stid[x]][(int) cm->ndtype[cm->ndidx[v]+is_insert]][(int) cm->stid[v]];
	      if(tmap_val == -1)
		{
		  printf("tmap ERROR 1\n");
		  printf("v: %d | pnum[v]: %d | plast[v]: %d | y: %d | x: %d | d1: %d | d2: %d | d3: %d\n", v, cm->pnum[v], cm->plast[v], y, x, ((int) cm->stid[x]), ((int) (cm->ndtype[cm->ndidx[v]+is_insert])), ((int) cm->stid[v]));
		  exit(1);
		}
	      if((x - start) < 0)
		sub_psi[v-start] += 0.;
	      else
		sub_psi[v-start] += sub_psi[x-start] * cm->t[x][(int) tmap_val];
	    }
	  if(v != end && is_insert) /* we don't want to include the probability of an 
				       insert self-loop to end itself */
	    {	  
	      sub_psi[v-start] += sub_psi[v-start] * (cm->t[v][0] / (1-cm->t[v][0])); 
	      /*the contribution of the self insertion loops*/
	    }
	}
    }
  to_return = (float) sub_psi[end-start];
  free(sub_psi);
  return to_return;
}

/**************************************************************************
 * EPN 03.19.06
 * check_psi_vs_phi_cp9()
 *
 * Purpose:  Check that psi and phi values for CM states that map to HMM states
 *           are within a certain leeway threshold. 
 *
 * Args:    
 * CM_t *cm          - the CM
 * double *psi       - psi[v] is expected number of times CM state v is entered
 * double **phi      - phi array, phi[k][v] is expected number of times
 *                     HMM state v (0 = match, 1 insert, 2 = delete) in 
 *                     node k is visited.
 * int ***hns2cs_map - 3D HMM node-state to CM state map, 1st D - HMM node index, 2nd D - 
 *                     HMM state (0(M), 1(I), 2(D)), 3rd D - 2 elements for up to 
 *                     2 matching CM states, value: CM states that map, -1 if none.
 * int hmm_M         - number of nodes in HMM
 * double threshold  - the threshold that mapping (potentially summed) psi and 
 *                     phi values are allowed to be different by, without throwing an error.
 * int print_flag    - TRUE to print out the values, FALSE not to 
 * Returns: TRUE: if CM and HMM are "close enough" (see code)
 *          FALSE: otherwise
 */
static int
check_psi_vs_phi_cp9(CM_t *cm, double *psi, double **phi, int ***hns2cs_map, int hmm_M,
                     double threshold, int debug_level)
{
  int v; /* CM state index*/ 
  int k;
  int y;
  double summed_psi;
  int *ap; /*CM state that maps to HMM state in node k*/
  int k_state; /*0, 1 or 2, state in hmm node k*/
  int dual_mapping_insert; /* TRUE if the HMM insert state maps to 2 CM insert states, 
			    * these situations are impossible for the HMM to mirror exactly
			    * (with a single insert state) due to the self insert transitions. 
			    * For example there's only one way for the HMM to emit 3 residues
			    * from this insert state, but 4 ways for the CM to emit 3 residues
			    * from the combo of the 2 CM inserts (3 from A, 3 from B, 1 from A and
			    * 2 from B, or 2 from A and 1 from B.) This is impossible to model
			    * with the same probabilities by the single HMM insert self transition. 
			    */
  int violation;
  int v_ct; /* Number of violations not involving dual insert states. */
  double diff;
  int ret_val; /* return value */

  if (cm->flags & CM_LOCAL_BEGIN) Die("internal error: we're in CM local mode while trying to build a CP9 HMM");

  ret_val = TRUE;
  v_ct = 0;
  ap = malloc(sizeof(int) * 2);

  if(debug_level >= 0)
    printf("\n");
  for(v = 0; v < cm->M; v++)
    {
      if(cm->stid[v] != BIF_B)
	{      
	  for(y = 0; y < cm->cnum[v]; y++)
	  {
	    if(debug_level > 1)
	      {
		printf("cm->t[%d][%d]: %f\n", v, y, cm->t[v][y]);
	      }
	  }
	}
      if(debug_level > 1)
	{
	  printf("\n");
	}
    }		       
  
  for (k = 0; k <= hmm_M; k++)
    {
      k_state = HMMMATCH;
      ap[0] = hns2cs_map[k][k_state][0];
      ap[1] = hns2cs_map[k][k_state][1];
      summed_psi = psi[ap[0]];
      if(ap[1] != -1)
	summed_psi += psi[ap[1]];
      violation = FALSE;
      diff = phi[k][0] - summed_psi;
      if((diff > threshold) || ((-1. * diff) > threshold))
	{
	  violation = TRUE;
	  v_ct++;
	}
      if(violation)
	printf("M k: %4d | phi: %f | psi: %f VIOLATION (%f)\n", k, phi[k][0], summed_psi, diff);
      else if(debug_level > 1)
	printf("M k: %4d | phi: %f | psi: %f\n", k, phi[k][0], summed_psi);

      k_state = HMMINSERT;
      ap[0] = hns2cs_map[k][k_state][0];
      ap[1] = hns2cs_map[k][k_state][1];
      summed_psi = psi[ap[0]];
      dual_mapping_insert = FALSE;
      if(ap[1] != -1)
	{
	  dual_mapping_insert = TRUE;
	  summed_psi += psi[ap[1]]; 
	  /* This is technically incorrect. Summing the psi of ap[0] and
	   * ap[1] here does not give the expected number of times either is
	   * visited. This is because they are not independent, some
	   * paths can go through both ap[0] and ap[1]. For a single HMM match
	   * and delete state, if two CM states map, they are from
	   * the same node's split set (ex: MATP_MP and MATP_ML) so
	   * they are independent (no path can go through both) and
	   * adding psi[ap[0]] and psi[ap[1]] is appropriate.  Here
	   * we don't bother to figure out what the correct psi value
	   * is (b/c its complex, and this is just a function for help 
	   * with debugging anyway).
	   */
	}
      violation = FALSE;
      diff = phi[k][1] - summed_psi;
      if((diff > threshold) || ((-1. * diff) > threshold))
	{
	  if(!(dual_mapping_insert)) 
	    {
	      violation = TRUE;
	      v_ct++;
	    }
	}
      if(violation)
	printf("I k: %4d | phi: %f | psi: %f VIOLATION (%f)\n", k, phi[k][1], summed_psi, diff);
      else if(dual_mapping_insert && debug_level > 1)
	printf("I k: %4d | phi: %f | psi: %f (DUAL INSERT so psi != phi)\n", k, phi[k][1], summed_psi);
      else if(debug_level > 1)
	printf("I k: %4d | phi: %f | psi: %f\n", k, phi[k][1], summed_psi);
      
      k_state = HMMDELETE;
      ap[0] = hns2cs_map[k][k_state][0];
      ap[1] = hns2cs_map[k][k_state][1];
     if(k == 0)
	summed_psi = 0.;      /*no such state in HMM or CM*/
      else
	{
	  summed_psi = psi[ap[0]];
	  if(ap[1] != -1)
	    summed_psi += psi[ap[1]];
	}
      violation = FALSE;
      diff = phi[k][2] - summed_psi;
      if((diff > threshold) || ((-1. * diff) > threshold))
	{
	  violation = TRUE;
	  v_ct++;
	}
      if(violation)
	printf("D k: %4d | phi: %f | psi: %f VIOLATION (%f)\n\n", k, phi[k][2], summed_psi, diff);
      else if(debug_level > 1)
	{
	  printf("D k: %4d | phi: %f | psi: %f\n\n", k, phi[k][2], summed_psi);
	}
    }
  free(ap);
  
  if(v_ct > 0)
    {
      printf("ERROR, %d HMM states violate the %f threshold b/t psi and phi.\n", v_ct, threshold);
      /*exit(1);*/
      ret_val = FALSE;
    }
  return ret_val;
}

/**************************************************************************
 * EPN 03.13.06
 * debug_print_cp9_params()
 *
 * Purpose:  Print out emission and transition probabilities and scores
 *           for a CM plan 9 HMM.
 *
 * Args:    
 * cplan9_s *hmm     - counts form CM plan 9 HMm
 * Returns: (void) 
 */
void
debug_print_cp9_params(struct cplan9_s *hmm)
{
  int k, i;

  printf("Printing CP9 HMM parameters in debug_print_cp9_params:\n\n");

  for(i = 0; i < MAXABET; i++)
    {
      printf("\tins[%d][%d] = %f | %d\n", 0, i, hmm->ins[0][i], hmm->isc[i][0]);
    }  
  printf("\n");

  k=0;
  printf("\tCTMM[%d] = %f | %d\n", k, hmm->t[0][CTMM], hmm->tsc[CTMM][0]);
  printf("\tCTMI[%d] = %f | %d\n", k, hmm->t[0][CTMI], hmm->tsc[CTMI][0]);
  printf("\tCTMD[%d] = %f | %d\n", k, hmm->t[0][CTMD], hmm->tsc[CTMD][0]);
  printf("\tCTIM[%d] = %f | %d\n", k, hmm->t[0][CTIM], hmm->tsc[CTIM][0]);
  printf("\tCTII[%d] = %f | %d\n", k, hmm->t[0][CTII], hmm->tsc[CTII][0]);
  printf("\tCTID[%d] = %f | %d\n", k, hmm->t[0][CTID], hmm->tsc[CTID][0]);
  printf("\tCTDM[%d] = %f | %d\n", k, hmm->t[0][CTDM], hmm->tsc[CTDM][0]);
  printf("\tCTDI[%d] = %f | %d\n", k, hmm->t[0][CTDI], hmm->tsc[CTDI][0]);
  printf("\tCTDD[%d] = %f | %d\n", k, hmm->t[0][CTDD], hmm->tsc[CTDD][0]);
  
  for(k = 1; k <= hmm->M; k++)
    {      
      printf("Node: %d\n", k);
      for(i = 0; i < MAXABET; i++)
	printf("mat[%3d][%3d] = %.3f | %d\n", k, i, hmm->mat[k][i], hmm->msc[i][k]);

      for(i = 0; i < MAXABET; i++)
	printf("ins[%3d][%3d] = %.3f | %d\n", k, i, hmm->ins[k][i], hmm->isc[i][k]);

      printf("\n");
      printf("\tCTMM[%d] = %f | %d\n", k, hmm->t[k][CTMM], hmm->tsc[CTMM][k]);
      printf("\tCTMI[%d] = %f | %d\n", k, hmm->t[k][CTMI], hmm->tsc[CTMI][k]);
      printf("\tCTMD[%d] = %f | %d\n", k, hmm->t[k][CTMD], hmm->tsc[CTMD][k]);
      printf("\tCTIM[%d] = %f | %d\n", k, hmm->t[k][CTIM], hmm->tsc[CTIM][k]);
      printf("\tCTII[%d] = %f | %d\n", k, hmm->t[k][CTII], hmm->tsc[CTII][k]);
      printf("\tCTID[%d] = %f | %d\n", k, hmm->t[k][CTID], hmm->tsc[CTID][k]);
      printf("\tCTDM[%d] = %f | %d\n", k, hmm->t[k][CTDM], hmm->tsc[CTDM][k]);
      printf("\tCTDI[%d] = %f | %d\n", k, hmm->t[k][CTDI], hmm->tsc[CTDI][k]);
      printf("\tCTDD[%d] = %f | %d\n", k, hmm->t[k][CTDD], hmm->tsc[CTDD][k]);
      printf("\t beg[%d] = %f | %d\n", k, hmm->begin[k], hmm->bsc[k]);
      printf("\t end[%d] = %f | %d\n", k, hmm->end[k], hmm->esc[k]);
      printf("\n");
    }
}

/**************************************************************************
 * EPN 03.13.06
 * check_cm_adj_bp()
 *
 * Purpose:  Check if two consensus columns (HMM nodes) are modelled by the
 *           same MATP node.
 * 
 * Args:    
 * CM_t *cm          - the CM
 * int *cc_node_map  - node that each consensus column maps to (is modelled by)
 *                     [1..hmm_ncc]
 * int hmm_M         - number of consensus columns (= num of HMM nodes)
 * Returns: TRUE if two adjacent consensus columns are modelled by the same MATP_nd
 *          FALSE if not
 */
static int
check_cm_adj_bp(CM_t *cm, int *cc_node_map, int hmm_M)
{
  int k, prev_k;
  prev_k = cc_node_map[1];
  for(k = 2; k <= hmm_M; k++)
    {
      if(cc_node_map[k] == prev_k)
	return TRUE;
      prev_k = cc_node_map[k];
    }
  return FALSE;
}

/**************************************************************************
 * EPN 03.13.06
 * CP9_check_wrhmm()
 *
 * Purpose:  Given a CM and a CM plan 9 hmm that is supposed to mirror 
 *           the CM as closely as possible (Weinberg-Ruzzo style, with
 *           differences due to differences in the CM Plan 9 architecture
 *           and the architecture that they use.
 * 
 *           Current strategy for checking is to build psi and phi arrays
 *           (see descriptions below) and check that they are "consistent",
 *           i.e. corresponding values are within a small threshold (0.0001).
 * 
 * Args:    
 * CM_t *cm          - the CM
 * cplan9_s *hmm     - the HMM
 * int ***hns2cs_map - 3D HMM node-state to CM state map, 1st D - HMM node index, 2nd D - 
 *                     HMM state (0(M), 1(I), 2(D)), 3rd D - 2 elements for up to 
 *                     2 matching CM states, value: CM states that map, -1 if none.
 *                     For example: CM states hns2cs_map[k][0][0] and hns2cs_map[k][0][1]
 *                                  map to HMM node k's match state.
 * int *cc_node_map  - node that each consensus column maps to (is modelled by)
 *                     [1..hmm_ncc]
 * int debug_level   - debugging level for verbosity of debugging printf statements
 *
 * Returns: TRUE: if CM and HMM are "close enough" (see code)
 *          FALSE: otherwise
 */
int 
CP9_check_wrhmm(CM_t *cm, struct cplan9_s *hmm, int ***hns2cs_map, int *cc_node_map,
		int debug_level)
{
  char     ***tmap;    /* hard-coded transition map */
  double    *psi;      /* psi[v] is the expected num times state v visited in CM */
  double   **phi;      /* phi[k][v] is expected number of times
			* state v (0 = match, 1 insert, 2 = delete) visited in CP9 HMM*/
  int k, i, j;         /* counters */
  double psi_vs_phi_threshold; /* the threshold that mapping (potentially summed) psi and 
				* phi values are allowed to be different by, without throwing
				* an error.*/
  int ret_val;         /* return value */
  psi = malloc(sizeof(double) * cm->M);
  make_tmap(&tmap);
  fill_psi(cm, psi, tmap);

  CPlan9Renormalize(hmm);
  /* Create and fill phi to check to make sure our HMM is "close enough" to our CM.
   * phi[k][0..2] is the expected number of times HMM node k state 0 (match), 1(insert),
   * or 2(delete) is entered. These should be *very close* (within 0.0001) to the psi 
   * values for the CM states that they map to (psi[v] is the expected number of times
   * state v is entered in the CM). The reason we can't demand closer precision than
   * 0.0001 is that some HMM insert states map to 2 CM insert states, and its impossible
   * to perfectly mirror two CM insert states with 1 HMM insert states (due to self-loop
   * issues).
   */
  phi = malloc(sizeof(double *) * (hmm->M+1));
  for(k = 0; k <= hmm->M; k++)
    {
      phi[k] = malloc(sizeof(double) * 3);
    }
  fill_phi_cp9(hmm, phi);

  /*debug_print_cp9_params(hmm);*/
  psi_vs_phi_threshold = 0.0001;
  if(check_cm_adj_bp(cm, cc_node_map, hmm->M))
    {
      psi_vs_phi_threshold = 0.01;
    /* if check_cm_adj_bp() returns TRUE, then the following rare (but possible) situation
     * is true. Two adjacent consensus columns are modelled by the same MATP node. In this
     * case it is difficult for the CM plan 9 architecture to mirror the insert state
     * that maps to both the MATP_IL and MATP_IR of this node. It is more difficult than
     * for any other possible CM topology situation, so we relax the threshold when
     * checking psi and phi.
     */
    }    
  ret_val = check_psi_vs_phi_cp9(cm, psi, phi, hns2cs_map, hmm->M, psi_vs_phi_threshold, debug_level);
  for(k = 0; k <= hmm->M; k++)
    {
      free(phi[k]);
    }
  free(phi);
  free(psi);
  for(i = 0; i < UNIQUESTATES; i++)
    {
      for(j = 0; j < NODETYPES; j++)
	free(tmap[i][j]);
      free(tmap[i]);
    }
  free(tmap);
  return ret_val;
}

/**************************************************************************
 * EPN 09.01.06
 * Function: CP9_check_wrhmm_by_sampling()
 *
 * Purpose:  Given a CM and a CM plan 9 hmm that is supposed to mirror 
 *           the CM as closely as possible (Weinberg-Ruzzo style, with
 *           differences due to differences in the CM Plan 9 architecture
 *           and the architecture that they use), check if the CM and CP9
 *           actually do correspond as closely as possible by generating
 *           an alignment from the CM using it to build a new CP9 HMM
 *           without using pseudocounts (the infinite MSA idea introduced
 *           by Zasha Weinberg in the ML-HMM Bioinformatics paper), and
 *           comparing the parameters of that CP9 HMM to the one we've
 *           built.
 *           
 *           NOTE: code to sample an alignment from the CM taken from
 *                 cmemit.c (which was ported from HMMER's hmmemit.c).
 * Args:    
 * CM_t *cm          - the CM
 * cplan9_s *hmm     - the HMM
 * int ***hns2cs_map - 3D HMM node-state to CM state map, 1st D - HMM node index, 2nd D - 
 *                     HMM state (0(M), 1(I), 2(D)), 3rd D - 2 elements for up to 
 *                     2 matching CM states, value: CM states that map, -1 if none.
 *                     For example: CM states hns2cs_map[k][0][0] and hns2cs_map[k][0][1]
 *                                  map to HMM node k's match state.
 * float thresh      - probability threshold for chi squared test
 * int nseq          - number of sequences to sample to build the new HMM.
 *
 * Returns: TRUE: if CM and HMM are "close enough" (see code)
 *          FALSE: otherwise
 */
int 
CP9_check_wrhmm_by_sampling(CM_t *cm, struct cplan9_s *hmm, int ***hns2cs_map, float thresh,
			    int nseq)
{
  int ret_val;         /* return value */
  Parsetree_t **tr;             /* Parsetrees of emitted aligned sequences */
  char    **dsq;                /* digitized sequences                     */
  SQINFO            *sqinfo;    /* info about sequences (name/desc)        */
  MSA               *msa;       /* alignment */
  float             *wgt;
  int i, idx, nd;
  int L;
  int apos;
  int *matassign;
  struct cp9trace_s **cp9_tr;   /* fake tracebacks for each seq            */
  struct cplan9_s  *shmm;       /* the new, CM plan9 HMM; built by sampling*/
  int msa_nseq;                 /* this is the number of sequences per MSA,
				 * current strategy is to sample (nseq/nseq_per_msa)
				 * alignments from the CM, and add counts from
				 * each to the shmm in counts form (to limit memory)
				 */
  int nsampled;                 /* number of sequences sampled thus far */
  int *dual_mapping_insert; /* dual_mapping_insert[nd] is TRUE if the HMM insert state 
			     * of node nd maps to 2 CM insert states, 
			     * these situations are impossible for the HMM to mirror exactly
			     * (with a single insert state) due to the self insert transitions. 
			     * For example there's only one way for the HMM to emit 3 residues
			     * from this insert state, but 4 ways for the CM to emit 3 residues
			     * from the combo of the 2 CM inserts (3 from A, 3 from B, 1 from A and
			     * 2 from B, or 2 from A and 1 from B.) This is impossible to model
			     * with the same probabilities by the single HMM insert self transition. 
			     */

  ret_val = TRUE;
  /* Determine which nodes of the HMM have dual mapping inserts */
  dual_mapping_insert = MallocOrDie(sizeof(int) * hmm->M);
  for(nd = 0; nd <= hmm->M; nd++)
    if(hns2cs_map[nd][HMMINSERT][1] != -1)
      dual_mapping_insert[nd] = 1;
    else
      dual_mapping_insert[nd] = 0;
      
  msa_nseq = 1000;
  /* Allocate and normalize the new HMM we're going to build by sampling from
   * the CM.
   */
  shmm = AllocCPlan9(hmm->M);
  ZeroCPlan9(shmm);

  CPlan9Renormalize(hmm);
  CMRenormalize(cm);

  /* sample MSA(s) from the CM */
  nsampled = 0;
  dsq    = MallocOrDie(sizeof(char *)             * msa_nseq);
  tr     = MallocOrDie(sizeof(Parsetree_t)        * msa_nseq);
  sqinfo = MallocOrDie(sizeof(SQINFO)             * msa_nseq);
  wgt    = MallocOrDie(sizeof(float)              * msa_nseq);
  FSet(wgt, msa_nseq, 1.0);

  while(nsampled < nseq)
    {
      /*printf("nsampled: %d\n", nsampled);*/
      if(nsampled != 0)
	{
	  MSAFree(msa);
	  free(matassign);
	  for (i = 0; i < msa_nseq; i++)
	    {
	      CP9FreeTrace(cp9_tr[i]);
	      FreeParsetree(tr[i]);
	      free(dsq[i]);
	    }
	  free(cp9_tr);
	}
      /* Emit msa_nseq parsetrees from the CM */
      if(nsampled + msa_nseq > nseq)
	msa_nseq = nseq - nsampled;
      for (i = 0; i < msa_nseq; i++)
	{
	  EmitParsetree(cm, &(tr[i]), NULL, &(dsq[i]), &L);
	  sprintf(sqinfo[i].name, "seq%d", i+1);
	  sqinfo[i].len   = L;
	  sqinfo[i].flags = SQINFO_NAME | SQINFO_LEN;
	}
      /* Build a new MSA from these parsetrees */
      msa = Parsetrees2Alignment(cm, dsq, sqinfo, NULL, tr, msa_nseq, TRUE);

      /* Determine match assignment from RF annotation
       */
      matassign = (int *) MallocOrDie (sizeof(int) * (msa->alen+1));
      matassign[0] = 0;
      for (apos = 0; apos < msa->alen; apos++)
	{
	  matassign[apos+1] = 0;
	  if (!isgap(msa->rf[apos])) 
	    matassign[apos+1] = 1;
	}
      /* make fake tracebacks for each seq */
      CP9_fake_tracebacks(msa->aseq, msa->nseq, msa->alen, matassign, &cp9_tr);
      
      /* build model from tracebacks (code from HMMER's modelmakers.c::matassign2hmm() */
      for (idx = 0; idx < msa->nseq; idx++) {
	CP9TraceCount(shmm, dsq[idx], msa->wgt[idx], cp9_tr[idx]);
      }
      nsampled += msa_nseq;
    }

  /* The new shmm is in counts form, filled with observations from MSAs sampled
   * from the CM. 
   * We want to do a series of chi-squared tests to determine the probability that
   * the observed samples from the CM were not taken from the corresponding CM Plan 9
   * HMM probability distributions (the CM Plan 9 is supposed to exactly mirror the 
   * CM in this way).
   */

  /* Node 0 is special */
  for(nd = 0; nd <= shmm->M; nd++)
    {
      if(!(CP9_node_chi_squared(hmm, shmm, nd, thresh, dual_mapping_insert[nd])))
	{
	  /*Die("ERROR: chi_squared test failed for node: %d\n", nd);*/
	  printf("ERROR: chi_squared test failed for node: %d\n", nd);
	  ret_val = FALSE;
	}
    }
  /*Next, renormalize shmm and logoddisfy it */
  CPlan9Renormalize(shmm);
  CP9Logoddsify(shmm);

  /*
  printf("PRINTING BUILT HMM PARAMS:\n");
  debug_print_cp9_params(hmm);
  printf("DONE PRINTING BUILT HMM PARAMS:\n");


  printf("PRINTING SAMPLED HMM PARAMS:\n");
  debug_print_cp9_params(shmm);
  printf("DONE PRINTING SAMPLED HMM PARAMS:\n");
  */

  /* Output the alignment */
  /*WriteStockholm(stdout, msa);*/

  free(dual_mapping_insert);
  FreeCPlan9(shmm);
  /* compare new CP9 to existing CP9 */

  return ret_val;
}

/* Function: CP9_fake_tracebacks()
 * 
 * Purpose:  From a consensus assignment of columns to MAT/INS, construct fake
 *           tracebacks for each individual sequence.
 *           
 * Args:     aseqs     - alignment [0..nseq-1][0..alen-1]
 *           nseq      - number of seqs in alignment
 *           alen      - length of alignment in columns
 *           matassign - assignment of column 1 if MAT, 0 if INS; 
 *                       [1..alen] (off one from aseqs)
 *           ret_tr    - RETURN: array of tracebacks
 *           
 * Return:   (void)
 *           ret_tr is alloc'ed here. Caller must free.
 */          
static void
CP9_fake_tracebacks(char **aseq, int nseq, int alen, int *matassign,
		struct cp9trace_s ***ret_tr)
{
  struct cp9trace_s **tr;
  int  idx;                     /* counter over sequences          */
  int  i;                       /* position in raw sequence (1..L) */
  int  k;                       /* position in HMM                 */
  int  apos;                    /* position in alignment columns   */
  int  tpos;			/* position in traceback           */
  int  first_match;             /* first match column */
  int  last_match;              /* last match column */

  tr = (struct cp9trace_s **) MallocOrDie (sizeof(struct cp9trace_s *) * nseq);
  
  first_match = -1;
  last_match  = -1;
  for (apos = 0; apos < alen; apos++)
    {
      if(matassign[apos+1] && first_match == -1) first_match = apos;
      if(matassign[apos+1]) last_match = apos;
    }

  for (idx = 0; idx < nseq; idx++)
    {
      CP9AllocTrace(alen+2, &tr[idx]);  /* allow room for B & E */
      
				/* all traces start with M_0 state (the B state)... */
      tr[idx]->statetype[0] = CSTB;
      tr[idx]->nodeidx[0]   = 0;
      tr[idx]->pos[0]       = 0;

      i = 1;
      k = 0;
      tpos = 1;

      for (apos = 0; apos < alen; apos++)
        {
	  tr[idx]->statetype[tpos] = CSTBOGUS; /* bogus, deliberately, to debug */

	  if (matassign[apos+1] && ! isgap(aseq[idx][apos]))
	    {			/* MATCH */
	      k++;		/* move to next model pos */
	      tr[idx]->statetype[tpos] = CSTM;
	      tr[idx]->nodeidx[tpos]   = k;
	      tr[idx]->pos[tpos]       = i;
	      i++;
	      tpos++;
	    }	      
          else if (matassign[apos+1])
            {                   /* DELETE */
	      /* We should be careful about S/W transitions; but we have 
	       * an ambiguity, based on the MSA, we can't tell if we
	       * did a local begin (some M->E transition) or if we
	       * went through a bunch of D state's before the first match 
	       * B->D_1 -> D_2 .... -> M_x. For now, we assume we're not in
	       * S/W mode, and treat it as the latter case, see
	       * HMMER's modelmaker.c:fake_tracebacks() for code
	       * on one *would* implement the S/W consideration IF
	       * there wasn't a B->D_1 transition allowed.
	       */
	      k++;		/* *always* move on model when match column seen */
	      tr[idx]->statetype[tpos] = CSTD;
	      tr[idx]->nodeidx[tpos]   = k;
	      tr[idx]->pos[tpos]       = 0;
	      tpos++;
            }
	  else if (! isgap(aseq[idx][apos]))
	    {			/* INSERT */
	      tr[idx]->statetype[tpos] = CSTI;
              tr[idx]->nodeidx[tpos]   = k;
              tr[idx]->pos[tpos]       = i;
	      i++;
	      tpos++;
	    }
	}
       /* all traces end with E state */
	      /* We should be careful about S/W transitions; but we have 
	       * an ambiguity, based on the MSA, we can't tell if we
	       * did a local end (some M->E transition) or if we
	       * went through a bunch of D state's before the final 
	       * D_M -> E transition. For now, we assume we're not in
	       * S/W mode, and treat it as the latter case, see
	       * HMMER's modelmaker.c:fake_tracebacks() for code
	       * on one *would* implement the S/W consideration IF
	       * there wasn't a D_M -> E transition allowed.
	       */
      tr[idx]->statetype[tpos] = CSTE;
      tr[idx]->nodeidx[tpos]   = 0;
      tr[idx]->pos[tpos]       = 0;
      tpos++;
      tr[idx]->tlen = tpos;
    }    /* end for sequence # idx */

  *ret_tr = tr;
  return;
}

/* Function: CP9TraceCount() 
 * EPN 09.04.06 based on Eddy's P7TraceCount() from HMMER's trace.c
 * 
 * Purpose:  Count a traceback into a count-based HMM structure.
 *           (Usually as part of a model parameter re-estimation.)
 *           
 * Args:     hmm   - counts-based CM Plan 9 HMM
 *           dsq   - digitized sequence that traceback aligns to the HMM (1..L)
 *           wt    - weight on the sequence
 *           tr    - alignment of seq to HMM
 *           
 * Return:   (void)
 */
void
CP9TraceCount(struct cplan9_s *hmm, char *dsq, float wt, struct cp9trace_s *tr)
{
  int tpos;                     /* position in tr */
  int i;			/* symbol position in seq */
  
  for (tpos = 0; tpos < tr->tlen; tpos++)
    {
      i = tr->pos[tpos];

      /* Emission counts. 
       */
      if (tr->statetype[tpos] == CSTM) 
	SingletCount(hmm->mat[tr->nodeidx[tpos]], dsq[i], wt);
      else if (tr->statetype[tpos] == CSTI) 
	SingletCount(hmm->ins[tr->nodeidx[tpos]], dsq[i], wt);

      /* State transition counts
       */
      switch (tr->statetype[tpos]) {
      case CSTB:
	switch (tr->statetype[tpos+1]) {
	case CSTM: hmm->begin[tr->nodeidx[tpos+1]] += wt; break;
	case CSTI: hmm->t[0][CTMI]                 += wt; break;
	case CSTD: hmm->t[0][CTMD]                 += wt; break;
	default:      
	  Die("illegal state transition %s->%s in traceback", 
	      CP9Statetype(tr->statetype[tpos]), 
	      CP9Statetype(tr->statetype[tpos+1]));
	}
	break;
      case CSTM:
	switch (tr->statetype[tpos+1]) {
	case CSTM: hmm->t[tr->nodeidx[tpos]][CTMM] += wt; break;
	case CSTI: hmm->t[tr->nodeidx[tpos]][CTMI] += wt; break;
	case CSTD: hmm->t[tr->nodeidx[tpos]][CTMD] += wt; break;
	case CSTE: hmm->end[tr->nodeidx[tpos]]     += wt; break;
	default:    
	  Die("illegal state transition %s->%s in traceback", 
	      CP9Statetype(tr->statetype[tpos]), 
	      CP9Statetype(tr->statetype[tpos+1]));
	}
	break;
      case CSTI:
	switch (tr->statetype[tpos+1]) {
	case CSTM: hmm->t[tr->nodeidx[tpos]][CTIM] += wt; break;
	case CSTI: hmm->t[tr->nodeidx[tpos]][CTII] += wt; break;
	case CSTD: hmm->t[tr->nodeidx[tpos]][CTID] += wt; break;
	case CSTE: 
	  /* This should only happen from the final insert (I_M) state */
	  if((tpos+1) != (tr->tlen-1))
	    Die("illegal state transition %s->%s (I is not final insert) in traceback", 
		CP9Statetype(tr->statetype[tpos]), 
		CP9Statetype(tr->statetype[tpos+1]));
	  hmm->t[tr->nodeidx[tpos]][CTIM] += wt; break;
	  break;
	default:    
	  Die("illegal state transition %s->%s in traceback", 
	      CP9Statetype(tr->statetype[tpos]), 
	      CP9Statetype(tr->statetype[tpos+1]));
	}
	break;
      case CSTD:
	switch (tr->statetype[tpos+1]) {
	case CSTM: hmm->t[tr->nodeidx[tpos]][CTDM] += wt; break;
	case CSTI: hmm->t[tr->nodeidx[tpos]][CTDI] += wt; break;
	case CSTD: hmm->t[tr->nodeidx[tpos]][CTDD] += wt; break;
	case CSTE: 
	  /* This should only happen from the final delete (D_M) state */
	  if((tpos+1) != (tr->tlen-1))
	    Die("illegal state transition %s->%s (D is not final delete) in traceback", 
		CP9Statetype(tr->statetype[tpos]), 
		CP9Statetype(tr->statetype[tpos+1]));
	  hmm->t[tr->nodeidx[tpos]][CTDM] += wt; break;
	  break;
	default:    
	  Die("illegal state transition %s->%s in traceback", 
	      CP9Statetype(tr->statetype[tpos]), 
	      CP9Statetype(tr->statetype[tpos+1]));
	}
	break;
      case CSTE:
	break; /* E is the last. It makes no transitions. */

      default:
	Die("illegal state %s in traceback", 
	    CP9Statetype(tr->statetype[tpos]));
      }
    }
}


/* Function: CP9Statetype()
 * 
 * Purpose:  Returns the state type in text.
 * Example:  CP9Statetype(M) = "M"
 */
char *
CP9Statetype(char st)
{
  switch (st) {
  case CSTM: return "M";
  case CSTD: return "D";
  case CSTI: return "I";
  case CSTB: return "B";
  case CSTE: return "E";
  default: return "BOGUS";
  }
}

/* Function: CP9_node_chi_squared()
 * 
 * Purpose : Given two CM Plan 9 HMMs, one in normalized form (ahmm), and one in
 *           raw counts form (shmm), determine for a specific node nd,
 *           the probability the samples implicit in the shmm were taken 
 *           from the corresponding probability distribution in ahmm using 
 *           the chi-squared test. Return FALSE IFF the probability that some 
 *           set of counts was not taken from the corresponding ahmm distribution
 *           (probability null hypothesis is rejected) is above the probability 
 *           threshold thresh.
 *             
 * Args:    
 * cplan9_s *ahmm    - an HMM with the probability distributions that define
 *                     the null hypothesis for the chi-squared tests (we
 *                     think the counts in shmm were taken from these distributions).
 * cplan9_s *shmm    - an HMM in counts form
 * int         nd    - node of the HMM to check. 
 * float   thresh    - probability of rejecting null hypotheses must be below
 *                     this value for all cases.
 * int dual_mapping_insert - TRUE if HMM node nd maps to 2 insert states in the CM,
 *                           for such states, we don't require that we pass the
 *                           chi squared test.
 */
int
CP9_node_chi_squared(struct cplan9_s *ahmm, struct cplan9_s *shmm, int nd, float thresh,
		     int dual_mapping_insert)
{
  double p;
  int x;
  float chi_sq;
  float m_nseq, i_nseq, d_nseq;
  float check_m_nseq, check_i_nseq;

  if(nd > shmm->M || nd >  ahmm->M)
    Die("ERROR CP9_node_chi_squared() is being grossly misused.\n");

  CPlan9Renormalize(ahmm);
  CP9Logoddsify(ahmm);

  /* First determine the sum of the counts for each state. This
   * is the number of samples that visited match, insert and 
   * delete state of this node.
   */
  
  if(nd == 0)
    m_nseq  = shmm->begin[1]; /* this is begin -> M_1 transition */
  else if(nd == shmm->M)
    m_nseq  = shmm->end[nd];   /* this is M_M -> end transition */
  else
    m_nseq  = shmm->t[nd][CTMM];
  m_nseq += shmm->t[nd][CTMI];
  m_nseq += shmm->t[nd][CTMD];
  
  i_nseq  = shmm->t[nd][CTIM];
  i_nseq += shmm->t[nd][CTII];
  i_nseq += shmm->t[nd][CTID];

  if(nd != 0) /* node 0 has no delete state */
    {
      d_nseq  = shmm->t[nd][CTDM];
      d_nseq += shmm->t[nd][CTDI];
      d_nseq += shmm->t[nd][CTDD];
    }

  if(nd != 0)
    {
      check_m_nseq = 0.;
      for (x = 0; x < MAXABET; x++) check_m_nseq += shmm->mat[nd][x];
      if((check_m_nseq >= m_nseq && ((check_m_nseq - m_nseq) > 0.0001)) ||
	 (check_m_nseq  < m_nseq && ((m_nseq - check_m_nseq) > 0.0001)))     
	{
	  Die("ERROR: node: %d has different number of sampled match emissions and transitions.\n");
	}
    }
  check_i_nseq = 0.;
  for (x = 0; x < MAXABET; x++) check_i_nseq += shmm->ins[nd][x];
  if((check_i_nseq >= i_nseq && ((check_i_nseq - i_nseq) > 0.0001)) ||
     (check_i_nseq  < i_nseq && ((i_nseq - check_i_nseq) > 0.0001)))     
    {
      Die("ERROR: node: %d has different number of sampled insert emissions and transitions.\n");
    }

  /* Check emissions */
  /* MATCH */
  if(nd != 0)
    {
      chi_sq = 0.;
      for (x = 0; x < MAXABET; x++) 
	chi_sq += (pow((ahmm->mat[nd][x] - (shmm->mat[nd][x] / m_nseq)), 2) / 
		   ((shmm->mat[nd][x] / m_nseq)));

      p = IncompleteGamma((MAXABET-1)/2., chi_sq/2.);
      /*printf("%4d E M P: %f m_nseq %f\n", nd, p, m_nseq);*/
      if((1. - p) > thresh)
	{
	  printf("%4d E M P: %f m_nseq %f\n", nd, p, m_nseq);
	  return FALSE;
	}
    }

  /* INSERT */
  chi_sq = 0.;
  for (x = 0; x < MAXABET; x++) 
    chi_sq += (pow((ahmm->ins[nd][x] - (shmm->ins[nd][x] / i_nseq)), 2) / 
	       ((shmm->ins[nd][x] / i_nseq)));

  p = IncompleteGamma((MAXABET-1)/2., chi_sq/2.);
  /*printf("%4d E I P: %f\n", nd, p);*/
  if((1. - p) > thresh)
    {
      printf("%4d E I P: %f i_nseq: %f\n", nd, p, i_nseq);
      return FALSE;
    }
  /* check transitions */
  /* out of match */
  if(nd == 0)
    chi_sq =  (pow((ahmm->begin[1] - (shmm->begin[1] / m_nseq)), 2) / 
	       ((shmm->begin[1] / m_nseq)));
  else if(nd == shmm->M)
    chi_sq =  (pow((ahmm->end[nd] - (shmm->end[nd] / m_nseq)), 2) / 
	       ((shmm->end[nd] / m_nseq)));
  else
    chi_sq   =  (pow((ahmm->t[nd][CTMM] - (shmm->t[nd][CTMM] / m_nseq)), 2) / 
	        ((shmm->t[nd][CTMM] / m_nseq)));
  chi_sq += (pow((ahmm->t[nd][CTMI] - (shmm->t[nd][CTMI] / m_nseq)), 2) / 
	        ((shmm->t[nd][CTMI] / m_nseq)));
  if(nd != shmm->M)
    chi_sq += (pow((ahmm->t[nd][CTMD] - (shmm->t[nd][CTMD] / m_nseq)), 2) / 
	          ((shmm->t[nd][CTMD] / m_nseq)));
  /* In CP9TraceCount(), we don't allow possibility of M->E transitions (local ends)
   * except if M = shmm->M, so we ignore them here also 
   */
  /*chi_sq += (pow((ahmm->end[nd]     - (shmm->end[nd]     / m_nseq)), 2) / 
	        ((shmm->end[nd]     / m_nseq)));
		
  */
  p = IncompleteGamma(1., chi_sq/2.);
  /*printf("%4d T M P: %f m_nseq: %f\n", nd, p, m_nseq);*/
  if((1. - p) > thresh)
    {
      printf("%4d T M P: %f m_nseq %f\n", nd, p, m_nseq);
      return FALSE;
    }

  if(!dual_mapping_insert)
    {  /* out of insert */
      chi_sq =  (pow((ahmm->t[nd][CTIM] - (shmm->t[nd][CTIM] / i_nseq)), 2) / 
		    ((shmm->t[nd][CTIM] / i_nseq)));
      chi_sq += (pow((ahmm->t[nd][CTII] - (shmm->t[nd][CTII] / i_nseq)), 2) / 
		    ((shmm->t[nd][CTII] / i_nseq)));
      if(nd != shmm->M)
	chi_sq += (pow((ahmm->t[nd][CTID] - (shmm->t[nd][CTID] / i_nseq)), 2) / 
	            ((shmm->t[nd][CTID] / i_nseq)));
      p = IncompleteGamma(1., chi_sq/2.);
      /*printf("%4d T I P: %f i_nseq: %f\n", nd, p, i_nseq);*/
      if((1. - p) > thresh)
	{
	  printf("%4d T I P: %f i_nseq: %f\n", nd, p, i_nseq);
	  return FALSE;
	}
    }
  else
    {
      /*printf("nd: %d is a dual mapping insert\n", nd);*/
    }      

  if(nd != 0)
    {
      /* out of delete */
      chi_sq =  (pow((ahmm->t[nd][CTDM] - (shmm->t[nd][CTDM] / d_nseq)), 2) / 
	            ((shmm->t[nd][CTDM] / d_nseq)));
      chi_sq += (pow((ahmm->t[nd][CTDI] - (shmm->t[nd][CTDI] / d_nseq)), 2) / 
	            ((shmm->t[nd][CTDI] / d_nseq)));
      if(nd != shmm->M)
	chi_sq += (pow((ahmm->t[nd][CTDD] - (shmm->t[nd][CTDD] / d_nseq)), 2) / 
		      ((shmm->t[nd][CTDD] / d_nseq)));
      p = IncompleteGamma(1., chi_sq/2.);
      /*printf("%4d T D P: %f d_nseq: %f\n", nd, p, d_nseq);*/
      if((1. - p) > thresh)
	{
	  printf("%4d T D P: %f d_nseq: %f\n", nd, p, d_nseq);
	  return FALSE;
	}
    }
  return TRUE;
}
  

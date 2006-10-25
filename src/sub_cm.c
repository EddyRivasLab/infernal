/*****************************************************************
 *  sub_cm.c
 * EPN 07.25.06
 * 
 * Building submodels (sub CMs) from a template CM, that represent
 * a contiguous subset of the consensus columns that were modelled 
 * by the template CM. Options exist to build a sub CM that
 * models the same consensus columns of the template, but without
 * structure outside a contiguous subset of the columns.
 *
 ***************************************************************** 
 * @LICENSE@
 ***************************************************************** 
 */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "squid.h"
#include "sre_stack.h"

#include "structs.h"
#include "funcs.h"
#include "hmmband.h"
#include "hmmer_funcs.h"

static void map_orig2sub_cm(CM_t *orig_cm, CM_t *sub_cm, int ***ret_orig2sub_smap, int ***ret_sub2orig_smap,
			    int sub_start, int sub_end);
static void map_orig2sub_cm2(CM_t *orig_cm, CM_t *sub_cm, int ***ret_orig2sub_smap, int ***ret_sub2orig_smap,
			     int sub_start, int sub_end);
static void map_orig2sub_cm3(CM_t *orig_cm, CM_t *sub_cm, int ***ret_orig2sub_smap, int ***ret_sub2orig_smap,
			     int sub_start, int sub_end, int **ret_sub2orig_id, int *imp_cc, int print_flag);
static int  map_orig2sub_cm_helper(CM_t *orig_cm, CM_t *sub_cm, int **orig2sub_smap, int **sub2orig_smap, int v_o, int v_s);
static void cm2sub_cm_check_id_next_node(CM_t *orig_cm, CM_t *sub_cm, int orig_nd, int sub_nd,
					 int *sub2orig_id, int *node_cc_left_o, int *node_cc_right_o,
					 int *node_cc_left_s, int *node_cc_right_s, int sub_start, int print_flag);
static void cm2sub_cm_emit_probs(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, int v_s, int v_o1, int v_o2,
				 int *sub2orig_id);
static void cm2sub_cm_trans_probs(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, char ***tmap, int v_s, 
				  int **orig2sub_smap, int **sub2orig_smap, int *sub2orig_id);
static void cm2sub_cm_trans_probs_S(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, char ***tmap, int v_start, 
				    int **orig2sub_smap, int **sub2orig_smap);
static void cm2sub_cm_trans_probs_B_E(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, char ***tmap, int v_end, 
				    int **orig2sub_smap, int **sub2orig_smap);
static void cm2sub_cm_trans_probs_B_E_dual_helper(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, char ***tmap, int v_s, 
						  int v_s_insert, int bif_end_yoffset, int **orig2sub_smap, 
						  int **sub2orig_smap);
static void cm2sub_cm_trans_probs_B_E2(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, char ***tmap, int v_end, 
				       int **orig2sub_smap, int **sub2orig_smap, int print_flag);
static void cm_add_single_trans(CM_t *orig_cm, CM_t *sub_cm, int **orig2sub_smap, int **sub2orig_smap, int orig_v, int orig_y, 
				int sub_v, int yoffset, double *orig_psi, char ***tmap);
static float cm_sum_subpaths(CM_t *orig_cm, CM_t *sub_cm, int **orig2sub_smap, int start, int end, int sub_start,
			     char ***tmap, double *orig_psi, int orig_insert1, int orig_insert2, int orig_insert3,
			     int orig_insert4);
static void debug_print_cm_params(CM_t *cm);
static int check_orig_psi_vs_sub_psi(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, double *sub_psi, 
				     int **orig2sub_smap, int **sub2orig_smap, double threshold, 
				     int print_flag);
static int orig2sub_state_check(int **orig2sub_smap, int orig_a, int orig_b);
static int sub_trans_check(CM_t *sub_cm, int **orig2sub_smap, int orig_a, int orig_b, int min_sub_v);
static int trans_check_helper(CM_t *cm, int a, int b);
static int cm2sub_cm_subtract_root_subpaths(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, char ***tmap, int sub_nd, 
					    int **orig2sub_smap, int **sub2orig_smap);
static void cm2sub_cm_subtract_root_subpaths2(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, char ***tmap, 
					     int **orig2sub_smap, int **sub2orig_smap, int print_flag);
static int cm2sub_cm_subtract_root_subpaths_helper(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, char ***tmap, 
						   int **orig2sub_smap, int **sub2orig_smap, int orig_il_or_s, int orig_ir);
static void cm2sub_cm_subtract_root_subpaths_helper2(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, char ***tmap, 
						     int **orig2sub_smap, int **sub2orig_smap, int orig_il, int orig_ir,
						     int il_iteration, int ir_iteration, int print_flag);
static void cm_subtract_subpaths(CM_t *orig_cm, CM_t *sub_cm, int **orig2sub_smap, int **sub2orig_smap, int orig_v, int orig_y, 
				 int orig_i, int sub_v, int yoffset, double *orig_psi, char ***tmap);
static void cm2sub_cm_find_impossible_misc_cases(CM_t *orig_cm, CM_t *sub_cm, int **orig2sub_smap, int **sub2orig_smap,
						 int *cc_node_map_o, int *node_cc_left_o,
						  int *node_cc_right_o, int *cc_node_map_s, int *node_cc_left_s, int *node_cc_right_s,
						  int *imp_cc, int spos, int epos, int print_flag);
static void cm2sub_cm_find_impossible_matr_cases(CM_t *orig_cm, CM_t *sub_cm, int *cc_node_map_o, int *node_cc_left_o, 
						 int *node_cc_right_o, int *cc_node_map_s, int *node_cc_left_s, int *node_cc_right_s,
						 int *imp_cc, int spos, int epos, int print_flag);
static void debug_print_misc_sub_cm_info(CM_t *orig_cm, CM_t *sub_cm, int **orig2sub_smap, int **sub2orig_smap, int spos, int epos,
					 int *cc_node_map_o, int *node_cc_left_o, int *node_cc_right_o);


/**************************************************************************** 
 * Function:  build_sub_cm()
 * EPN 08.28.06 
 *
 * Purpose:   Given a template CM (orig_cm) built from an MSA of consensus columns 
 *            1..M, build a new, "sub CM" (sub_cm), that only models well-nested base pairs 
 *            that occur (have both left and right halves) between consensus columns 
 *            [spos..epos]. If do_fullsub == TRUE, the sub CM will still model 
 *            all the consensus columns orig_cm models, with only structure removed.
 *            If do_fullsub == FALSE, the sub_cm will only model the consensus
 *            columns between [spos..epos].
 *
 *            We attempt to construct the sub CM such that: 
 *  
 *            1. All 'mapped' states in the orig_cm and sub_cm have identical 
 *               'psi' values, i.e. the expected number of times entered in a single parse
 *               is identical.
 *
 *            2. Alignments sampled from it should 'follow the same probability distributions'
 *               as alignments sampled from the template CM and then (if do_fullsub == TRUE)
 *               truncated, removing columns outside [spos..epos] Specifically, 'follow
 *               the same probability distributions' means if we built a ML HMM (Weinberg 
 *               style) from the resulting infinite alignments, we would have identical 
 *               HMMs for the sub-CM alignments and for the truncated template CMs. It turns
 *               out this is impossible if do_fullsub == FALSE for certain situations
 *               that cause differences in the topology of the orig_cm and sub_cm. These
 *               situations, so-called 'impossible cases', cause the distributions out of a
 *               single ML HMM node to be different between the orig_cm, but are relatively 
 *               rare (a rough, rough estimate is 1% of ML HMM nodes are affected).  We can 
 *               predict when all of these situations will occur, but sometimes when we 
 *               predict a ML HMM node will be off in this way, it is in fact identical 
 *               between the sub_cm and orig_cm for reasons I can't quite grasp.
 *            
 *            This function builds a sub_cm and then potentially performs up to 3 checks 
 *            to see if 1 and 2 above are satisfied (it only checks if 2 holds for 
 *            non-impossible cases, and that we can accurately predict the impossible 
 *            cases). The 3 potential checks are: 
 * 
 *            A. Calculate the expected number of times that each state in the orig_cm
 *               and sub_cm is entered (psi[v] = exp # times state v is entered), and 
 *               make sure these values are within 'threshold' of each other for all
 *               states that map between the orig_cm and sub_cm. This is performed
 *               by the 'check_orig_psi_vs_sub_psi' function, which is called within
 *               this function if either the 'do_acheck' and/or 'do_scheck' parameters
 *               passed in are TRUE.
 *
 *            B. Analytically (not by sampling) building 2 ML HMMs, one from the orig_cm 
 *               and one from the sub_cm, predicting the impossible cases and testing t
 *               o make sure all the corresponding parameters for non-impossible cases 
 *               are within 'threshold' of each other. This is done within the
 *               check_sub_cm() function which is called from this function if the
 *               'do_acheck' parameter passed in is set to TRUE. 
 *            
 *            C. Build a CP9 ML HMM analytically from the sub_cm. Sample a deep MSA
 *               from the orig_cm and potentially truncate outside [spos..epos]. Use
 *               counts from the MSA to build a ML HMM for the orig_cm. Then perform
 *               chi-squared tests to see if the counts in the orig_cm's HMM could
 *               have been generated by the distributions in the sub_cm's HMM. This
 *               check is performed if the 'do_scheck' parameter passed in is set
 *               to TRUE.
 *
 * Args:      orig_cm      - the original model, which we're going to (potentially)
 *                           remove some structure from
 *            ret_sub_cm   - the new sub_cm built from orig_cm with some structure removed.
 *            spos         - the first position (consensus column) we want to keep structure
 *                           for, and if do_subfull == FALSE the first column we want 
 *                           the sub_cm to model at all 
 *            epos         - the last position we want to keep structure for
 *                           and if do_subfull == FALSE the last column we want 
 *                           the sub_cm to model at all 
 *     ***ret_orig2sub_smap- 2D state map from orig_cm (template) to sub_cm.
 *                           1st dimension - state index in orig_cm 
 *                           2nd D - 2 elements for up to 2 matching sub_cm states, 
 *                           SHOULD NOT BE PRE-ALLOCATED
 *     ***ret_sub2orig_smap- 2D state map from orig_cm (template) to sub_cm.
 *                           1st dimension - state index in sub_cm (0..sub_cm->M-1)
 *                           2nd D - 2 elements for up to 2 matching orig_cm states, 
 *                           SHOULD NOT BE PRE-ALLOCATED
 *            **ret_imp_cc - [0..(epos-spos+1)] ret_imp_cc[k] != 0 if it is 
 *                           impossible for CP9 node (consensus column) k to be
 *                           calculated in the sub_cm to have distros to match the
 *                           corresponding CP9 node in the original CM - due to
 *                           topological differences in the architecture of the
 *                           sub_cm and orig_cm.
 *        **ret_apredict_ct- For 1 'atest's: the  number of times we predict we'll 
 *			     fail the test for an HMM node for
 *			     each of 7 cases - 5 different reasons we predict we'll fail.
 *  **ret_apredict_total_ct- For ALL 'atest's: the  number of times we predict we'll 
 *			     fail the test for an HMM node for
 *			     each of 7 cases - 5 different reasons we predict we'll fail.
 *        **ret-spredict_ct- For 1 sampling test: the  number of times we predict we'll 
 *			     fail the test for an HMM node for
 *			     each of 7 cases - 5 different reasons we predict we'll fail.
 *  **ret_spredict_total_ct- For ALL sampling tests: the  number of times we predict we'll 
 *			     fail the test for an HMM node for
 *			     each of 7 cases - 5 different reasons we predict we'll fail.
 *            threshold    - threshold for psi and phi tests
 *            do_fullsub   - TRUE to build a sub CM that models ALL consensus columns
 *                           in orig_cm, removing structure outside [spos..epos]
 *            do_acheck    - TRUE to perform psi and analytical HMM tests (A & B above)
 *            do_scheck    - TRUE to perform psi and sampling  HMM tests (A & C above)
 *            chi_thresh   - if a chi-square test during the sampling check is below
 *                           this threshold, fail.
 *            nsamples     - number of samples to use for sampling check (do_scheck)
 *                           NULL is acceptable if do_scheck == FALSE
 *            print_flag   - TRUE to print debugging statements    
 *
 * Returns:   int            TRUE if all checks passed or no checks were performed,
 *                           FALSE otherwise
 */
int 
build_sub_cm(CM_t *orig_cm, CM_t **ret_cm, int spos, int epos, int ***ret_orig2sub_smap, 
	     int ***ret_sub2orig_smap, int **ret_imp_cc, int **ret_apredict_ct, 
	     int **ret_awrong_predict_ct, int **ret_spredict_ct, int **ret_swrong_predict_ct, 
	     float threshold, int do_fullsub, int do_acheck, int do_scheck, float chi_thresh, 
	     int nsamples, int print_flag)
{
  CM_t            *sub_cm;      /* new covariance model, a submodel of the template */
  CMConsensus_t   *con;         /* growing consensus info for orig_cm               */
  Parsetree_t     *mtr;         /* master structure tree from the alignment         */
  char            *sub_cstr;    /* consensus substructure display string            */
  int             *sub_ct;	/* 0..con->clen-1 base pair partners array          */
  int              cpos;        /* position counter within orig_cm                  */
  int              sub_cpos;    /* position counter within sub_cm                   */
  CMEmitMap_t     *orig_emap;   /* consensus emit map for the original, template CM */
  CMEmitMap_t     *sub_emap;    /* consensus emit map for the sub_cm                */
  char          ***tmap;        /* hard-coded transition map, for convenience       */
  double          *orig_psi;    /* expected num times each state visited in orig_cm */
  double          *sub_psi;     /* expected num times each state visited in sub_cm  */
  int              v_o;         /* state counter for orig_cm                        */             
  int              v_o1;        /* state counter for orig_cm                        */             
  int              v_o2;        /* state counter for orig_cm                        */             
  int              v_s;         /* state counter for sub_cm                         */             
  int              v_s1;        /* state counter for sub_cm                         */             
  int              v_s2;        /* state counter for sub_cm                         */             
  int              n_s;         /* node  counter for sub_cm                         */             
  int              i,j;         /* counters                                         */
  FILE            *ofp;         /* an open output file                              */
  int              cc;          /* counter over consensus columns                   */
  int              ret_val;     /* keeps track of whether checks were passed        */
  int              model_start; /* first consensus (match) column of the orig_cm to 
				 * model with the sub_cm (if do_fullsub, this is 1) */
  int              model_end;   /* last consensus (match) column of the orig_cm to 
				 * model with the sub_cm (if do_fullsub, this is  
				 * total number of columns modelled by the orig_cm  */
  int             *imp_cc;      /* imp_cc[cc] != 0 if we predict that HMM node cc     
                                 * will be impossible to accurately predict (see 
				 * above), and imp_cc[cc] == 0 otherwise; 
				 * [0..(epos-spos+1)] if do_fullsub == FALSE
				 * [0..ncols] otherwise, ncols = num match cols in
				 * orig_cm                                          */
  int             *sub2orig_id; /* sub2orig_id[v_s] == TRUE if sub_cm state v_s maps
                                 * to an identical state (same emissions/transitions)
				 * in the orig_cm, FALSE otherwise. [0..sub_cm->M]  */
  int             *apredict_ct; /* Array of length 6 that holds the number of 
				 * times we predict an HMM node will fail */
  int             *spredict_ct; /* Array of length 6 that holds the number of 
				 * times we predict an HMM node will fail */
  int        *awrong_predict_ct;/* Array of length 6 that holds the number of 
			         * times we predict an HMM node will fail 
				 * it doesn't for the analytical test.          */
  int        *swrong_predict_ct;/* Array of length 6 that holds the number of 
				 * times we predict an HMM node will fail but
				 * it doesn't for the sampling test.          */
  int              **orig2sub_smap;/* 2D state map from orig_cm (template) to sub_cm.
				    * 1st dimension - state index in orig_cm 
				    * 2nd D - 2 elements for up to 2 matching sub_cm states */
  int              **sub2orig_smap;/* 2D state map from orig_cm (template) to sub_cm.
				    * 1st dimension - state index in sub_cm (0..sub_cm->M-1)
				    * 2nd D - 2 elements for up to 2 matching orig_cm states */
  int              **orig2sub_smap_copy;/* 2D state map from orig_cm (template) to sub_cm.
				    * 1st dimension - state index in orig_cm 
				    * 2nd D - 2 elements for up to 2 matching sub_cm states */
  int              **sub2orig_smap_copy;/* 2D state map from orig_cm (template) to sub_cm.
				    * 1st dimension - state index in sub_cm (0..sub_cm->M-1)
				    * 2nd D - 2 elements for up to 2 matching orig_cm states */

  apredict_ct = MallocOrDie(sizeof(int) * 7);
  spredict_ct = MallocOrDie(sizeof(int) * 7);
  awrong_predict_ct = MallocOrDie(sizeof(int) * 7);
  swrong_predict_ct = MallocOrDie(sizeof(int) * 7);
  for(i = 0; i <= 6; i++)
    {
      apredict_ct[i] = 0;
      spredict_ct[i] = 0;
      awrong_predict_ct[i] = 0;
      swrong_predict_ct[i] = 0;
    }

  /* Get the consensus sequence and consensus structure information from the original CM */
  con = CreateCMConsensus(orig_cm, 3.0, 1.0);
  if(print_flag)
    {
      printf("con->cseq    : %s\n", con->cseq);
      printf("con->cstr    : %s\n", con->cstr);
      printf("spos         : %d\n", spos);
      printf("epos         : %d\n", epos);
      printf("clen         : %d\n", con->clen);
    }

  if(do_fullsub)
    {
      model_start = 1;
      model_end   = con->clen;
    }
  else
    {
      model_start = spos;
      model_end   = epos;
    }

  /* Fill a new ct array for the sub_cm. The sub_cm will only model the consensus columns
   * between model_start and model_end, and only the structure between spos
   * and epos. First copy the template (original) CMs ct array but only for the 
   * appropriate consensus columns that lie in between both structure and model boundarIesx
   * Next, eliminate any structure that lies outside the structure boundaries.
   */

  sub_ct = MallocOrDie(sizeof(int) * (model_end - model_start + 1));
  /* First just copy ct array for model boundaries from con->ct */
  for (cpos = (model_start-1); cpos < model_end; cpos++)
    {
      sub_cpos = cpos - (model_start-1);
      if(con->ct[cpos] != -1 && 
	 (con->ct[cpos] <  (model_start-1) ||
	  con->ct[cpos] >=  model_end))
	sub_ct[sub_cpos] = -1;
      else
	sub_ct[sub_cpos] = con->ct[cpos];
    }
  /* Second remove structure outside structural boundaries */
  for (cpos = (model_start-1); cpos < model_end; cpos++)
    {
      sub_cpos = cpos - (model_start-1);
      if ((cpos+1) < spos || (cpos+1) > epos) /* cpos goes 1..clen, but ct is indexed
					       * 0..clen-1.*/
	{ 
	  /* CreateCMConsensus() uses -1 in ct[] to indicate single 
	   * stranded (different convention than WUSS2ct()). */
	  if (sub_ct[sub_cpos] != -1) 
	    sub_ct[sub_ct[sub_cpos]] = -1; 
	  sub_ct[sub_cpos] = -1;
	}
    }

  /* Construct the new structure ss_cons based on the template CM 
   * ct array.
   * We could do this similar to how display.c::CreateCMConsensus()
   * does it to get the fully formatted WUSS ([{<>}]) string but 
   * lazily we just do <> bps here.
   */
  sub_cstr = MallocOrDie(sizeof(char) * (model_end - model_start + 2));
  for (cpos = (model_start-1); cpos < model_end; cpos++)
    {
      sub_cpos = cpos - (model_start-1);
      if(sub_ct[sub_cpos] == -1)         sub_cstr[sub_cpos] = '.'; 
      else if (sub_ct[sub_cpos]  > cpos) sub_cstr[sub_cpos] = '<';
      else if (sub_ct[sub_cpos]  < cpos) sub_cstr[sub_cpos] = '>';
      else Die("ERROR: weird error in build_sub_cm()\n");
    }
  sub_cstr[(model_end-model_start+1)] = '\0';

  /* Build the new sub_cm given the new consensus structure. But don't
   * parameterize it yet.
   */
  ConsensusModelmaker(sub_cstr, (model_end-model_start+1), &sub_cm, &mtr);
  /* Rebalance the CM for optimization of D&C */
  CM_t *new;
  new = CMRebalance(sub_cm);
  FreeCM(sub_cm);
  sub_cm = new;

  if(print_flag)
    {
      printf("\n\norig struct: %s\n", con->cstr);
      printf("\n\nnew struct : %s\n", sub_cstr);
    }
  /* Parameterize the new sub_cm based on the template CM. 
   * First we need the emit maps for the template CM and
   * the new sub_cm.
   */
  orig_emap = CreateEmitMap(orig_cm);
  sub_emap = CreateEmitMap(sub_cm);

  if(print_flag)
    {
      printf("\n\n\nDumping original CM emitmap\n");
      DumpEmitMap(stdout, orig_emap, orig_cm);
      printf("\n\n\nDumping sub_cm emitmap\n");
      DumpEmitMap(stdout, sub_emap, sub_cm);
    }
  /* Allocate and initialize the imp_cc[] array */
  imp_cc = MallocOrDie(sizeof(int) * (model_end-model_start+2));
  for(cc = 0; cc <= (model_end-model_start+1); cc++)
    {
      imp_cc[cc] = FALSE;
    }
  /* Map states from orig_cm to sub_cm and vice versa. */
  map_orig2sub_cm3(orig_cm, sub_cm, &orig2sub_smap, &sub2orig_smap, model_start, model_end,
		   &sub2orig_id, imp_cc, print_flag);

  /* Fill orig_psi, which we need to determine the sub_cm parameters. */
  make_tmap(&tmap);
  orig_psi = malloc(sizeof(double) * orig_cm->M);
  fill_psi(orig_cm, orig_psi, tmap);

  CMZero(sub_cm);
  CMSetNullModel(sub_cm, orig_cm->null);
  sub_cm->el_selfsc = orig_cm->el_selfsc;

  /* Fill in emission probabilities */
  for(v_s = 0; v_s < sub_cm->M; v_s++)
    {
      if(sub_cm->sttype[v_s] != S_st &&
	 sub_cm->sttype[v_s] != D_st &&
	 sub_cm->sttype[v_s] != B_st &&
	 sub_cm->sttype[v_s] != E_st &&
	 sub_cm->sttype[v_s] != EL_st)
	cm2sub_cm_emit_probs(orig_cm, sub_cm, orig_psi, v_s, sub2orig_smap[v_s][0], sub2orig_smap[v_s][1], sub2orig_id);
    }
  /* Fill in transition virtual counts.
   * First handle non-B,S,E states, we'll deal with these later.
   * The reason we have to wait is that we can't (I don't think at least) 
   * unambiguously map the sub_cm B, S, or E states to orig_cm states.
   */
  for(v_s = 0; v_s < sub_cm->M; v_s++)
    {
      if(v_s == 0 || 
	 (sub_cm->sttype[v_s] != S_st &&
	  sub_cm->sttype[v_s] != B_st &&
	  sub_cm->sttype[v_s] != E_st))
	cm2sub_cm_trans_probs(orig_cm, sub_cm, orig_psi, tmap, v_s, orig2sub_smap, sub2orig_smap, sub2orig_id);
    }

  //printf("\n\n\nAddressing 090806\n");
  /* Address problem 090806 (in the 00LOG of ~/notebook/6_0725_inf_sub_cm/), by
   * retraversing the structure and subtracting out subpaths that have been counted twice
   * for a special situation involving the two inserts of ROOT and MATP states 
   */
    for(n_s = 0; n_s < sub_cm->nodes; n_s++)
    {
      if(sub_cm->ndtype[n_s] == MATP_nd)
	{
	  
	  if((sub2orig_smap[sub_cm->nodemap[n_s] + 4][1] != -1) ||
	     (sub2orig_smap[sub_cm->nodemap[n_s] + 5][1] != -1))
	    Die("ERROR, MATP_IL or MATP_IR node: %d map to 2 cm states\n", n_s);
	  if(sub2orig_smap[sub_cm->nodemap[n_s] + 4][0] != (sub2orig_smap[sub_cm->nodemap[n_s] + 5][0] - 1))
	    Die("ERROR, MATP_IL or MATP_IR node: %d don't map to adjacent orig_cm states\n");
	}
      if(sub_cm->ndtype[n_s] == ROOT_nd && sub_cm->ndtype[n_s+1] != BIF_nd) /* ROOT->BIFs are handled special
									     * (see next loop) */
	//cm2sub_cm_subtract_root_subpaths(orig_cm, sub_cm, orig_psi, tmap, n_s, orig2sub_smap, sub2orig_smap);
	cm2sub_cm_subtract_root_subpaths2(orig_cm, sub_cm, orig_psi, tmap, orig2sub_smap, sub2orig_smap, 
					  print_flag);
    }

  /* Go back through and fill in the transitions into E and B states and out of S states */
  for(v_s = 0; v_s < sub_cm->M; v_s++)
    {
      if(sub_cm->sttype[v_s] == S_st)
	cm2sub_cm_trans_probs_S(orig_cm, sub_cm, orig_psi, tmap, v_s, orig2sub_smap, sub2orig_smap);

      if(sub_cm->sttype[v_s] == E_st)
	//cm2sub_cm_trans_probs_B_E(orig_cm, sub_cm, orig_psi, tmap, v_s, orig2sub_smap, sub2orig_smap);
	cm2sub_cm_trans_probs_B_E2(orig_cm, sub_cm, orig_psi, tmap, v_s, orig2sub_smap, sub2orig_smap,
				   print_flag);

      if(sub_cm->sttype[v_s] == B_st)
	{
	  //cm2sub_cm_trans_probs_B_E(orig_cm, sub_cm, orig_psi, tmap, v_s, orig2sub_smap, sub2orig_smap);
	  cm2sub_cm_trans_probs_B_E2(orig_cm, sub_cm, orig_psi, tmap, v_s, orig2sub_smap, sub2orig_smap,
				     print_flag);
	  /* convention is to leave transitions out of BIF_B as 0.0, all the code knows they're obligate */
	}
    }

  /* Find and detach original CM insert states that map to 2 sub CM states. 
   * Since we only need 1 sub CM state for each of these cases, we make one of 
   * them 'dead', by setting transitions into it as 0.0.
   */
  
  if(print_flag) printf("\n\nLooking for dead inserts in sub CM...\n");
  for(v_s1 = 0; v_s1 < sub_cm->M; v_s1++)
    {
      if(sub_cm->sttype[v_s1] == IL_st || sub_cm->sttype[v_s1] == IR_st)
	{
	  v_o1 = sub2orig_smap[v_s1][0];
	  v_o2 = sub2orig_smap[v_s1][1];
	  if(orig2sub_smap[v_o1][1] != -1)
	    {
	      if(orig2sub_smap[v_o1][0] == v_s1)
		v_s2 = orig2sub_smap[v_o1][1];
	      else
		v_s2 = orig2sub_smap[v_o1][0];

	      /* v_s1 maps to an orig_cm state that maps to 2 sub_cm insert states */
	      if(print_flag) printf("Found what should be a dead state: v_s1: %d maps to v_o1: %d which maps to sub CM states %d and %d\n", v_s1, v_o1, orig2sub_smap[v_o1][0], orig2sub_smap[v_o1][1]);
	      if(v_o2 != -1)
		{
		  /* v_s1 maps to 2 orig_cm inserts, we have to make sure that v_s2 also maps to 
		   * both of these inserts */
		  if(!((v_o1 == sub2orig_smap[v_s2][0] && v_o2 == sub2orig_smap[v_s2][1]) ||
		       (v_o1 == sub2orig_smap[v_s2][1] && v_o2 == sub2orig_smap[v_s2][0])))
		    Die("ERROR: found sub_cm dead state but sub_cm state %d maps to 2 orig_cm states (%d and %d), and sub_cm state %d does not map to these same 2 states (but rather: %d and %d). I thought this was impossible (case 1).\n", v_s1, v_s2, v_o1, v_o2, sub2orig_smap[v_s2][0], sub2orig_smap[v_s2][1]);;
		}
	      if(!(cm_detach_state(sub_cm, v_s1, v_s2)))
		Die("ERROR: problem with finding a state to kill: neither v_s1: %d and v_s2: %d map to END_E-1!\n", v_s1, v_s2);
	    }
	}
    }

  /* Finally renormalize the CM */
  CMRenormalize(sub_cm);
  CMLogoddsify(sub_cm);
  if(print_flag)  printf("done renormalizing\n");

  if(print_flag)
    {
      ofp = fopen("sub.cm", "w");
      if(print_flag)  printf("%-40s ... ", "Saving model to file"); fflush(stdout);
      if(print_flag) CMFileWrite(ofp, sub_cm, FALSE);
      if(print_flag)  printf("done.\n");
    }

  if(do_acheck || do_scheck)
    {
      sub_psi = malloc(sizeof(double) * sub_cm->M);
      fill_psi(sub_cm, sub_psi, tmap);
    }

  if(print_flag)
  {
    printf("\nDEBUG PRINT OF ORIG_CM PARAMETERS:\n");
    debug_print_cm_params(orig_cm);
    printf("\nDEBUG PRINT OF SUB_CM PARAMETERS:\n");
    debug_print_cm_params(sub_cm);
 }    

  /* Perform checks if desired */
  ret_val = 0;
  if(do_acheck || do_scheck)
    ret_val += (!check_orig_psi_vs_sub_psi(orig_cm, sub_cm, orig_psi, sub_psi, orig2sub_smap, 
					   sub2orig_smap, threshold, print_flag));
  if(do_acheck)
    ret_val += (!check_sub_cm(orig_cm, sub_cm, model_start, model_end, threshold, imp_cc, 
			      apredict_ct, awrong_predict_ct, print_flag));
  if(do_scheck)
    ret_val += (!check_sub_cm_by_sampling(orig_cm, sub_cm, model_start, model_end, chi_thresh, 
					  nsamples, imp_cc, spredict_ct, swrong_predict_ct, print_flag));

  /* Cleanup and exit. */
  for(i = 0; i < UNIQUESTATES; i++)
    {
      for(j = 0; j < NODETYPES; j++)
	free(tmap[i][j]);
      free(tmap[i]);
    }
  free(tmap);

  free(sub_cstr);
  free(sub_ct);
  FreeEmitMap(orig_emap);
  FreeEmitMap(sub_emap);

  free(sub2orig_id);

  FreeCMConsensus(con);

  FreeParsetree(mtr);

  free(orig_psi);
  if(do_acheck || do_scheck)
    free(sub_psi);

  *ret_cm = sub_cm;
  if(ret_imp_cc != NULL)
    *ret_imp_cc = imp_cc;
  if(ret_apredict_ct != NULL)
    *ret_apredict_ct = apredict_ct;
  if(ret_spredict_ct != NULL)
    *ret_spredict_ct = spredict_ct;
  if(ret_awrong_predict_ct != NULL)
    *ret_awrong_predict_ct = awrong_predict_ct;
  if(ret_swrong_predict_ct != NULL)
    *ret_swrong_predict_ct = swrong_predict_ct;

  /* temporary solution to weird memory errors, copy the smaps */
  if(ret_orig2sub_smap != NULL)
    {
      orig2sub_smap_copy = MallocOrDie(sizeof(int *) * (orig_cm->M));
      for(v_o = 0; v_o < orig_cm->M; v_o++)
	{
	  orig2sub_smap_copy[v_o] = MallocOrDie(sizeof(int) * 2);
	  orig2sub_smap_copy[v_o][0] = orig2sub_smap[v_o][0];
	  orig2sub_smap_copy[v_o][1] = orig2sub_smap[v_o][1];
	}
      *ret_orig2sub_smap = orig2sub_smap_copy;
    }
  if(ret_sub2orig_smap != NULL)
    {
      sub2orig_smap_copy = MallocOrDie(sizeof(int *) * (sub_cm->M));
      for(v_s = 0; v_s < sub_cm->M; v_s++)
	{
	  sub2orig_smap_copy[v_s] = MallocOrDie(sizeof(int) * 2);
	  sub2orig_smap_copy[v_s][0] = sub2orig_smap[v_s][0];
	  sub2orig_smap_copy[v_s][1] = sub2orig_smap[v_s][1];
	}
      *ret_sub2orig_smap = sub2orig_smap_copy;
    }

  for(v_s = 0; v_s < orig_cm->M; v_s++)
    free(orig2sub_smap[v_s]);
  free(orig2sub_smap);
  for(v_s = 0; v_s < sub_cm->M; v_s++)
    free(sub2orig_smap[v_s]);
  free(sub2orig_smap);

  if(ret_val == 0)
    return TRUE; /* Passed all checks (or no checks were performed) */
  else
    return FALSE; /* Failed at least 1 check */
}

/**************************************************************
 * Function: CP9NodeForPosn()
 * EPN 07.25.06 Benasque, Spain
 * 
 * Purpose:  Determine the node of the CP9 HMM that is most likely to 
 *           have emitted (from either its Match or Insert state)
 *           a given posn in the target sequence.
 *
 * Args:     hmm       - the CM plan 9 HMM
 *           i0        - first posn of target subseq with info in posterior matrix
 *           j0        - last posn of target subseq with info in posterior matrix
 *           x         - posn of target subsequence we're interested in
 *           L         - last position of target sequence 
 *           post      - the posterior matrix for the hmm
 *           ret_node  - RETURN: index of node with highest probability of emitting x
 *           ret_type  - RETURN: type of state in ret_node with highest probability 
 *
 */
void
CP9NodeForPosn(struct cplan9_s *hmm, int i0, int j0, int x, struct cp9_dpmatrix_s *post, 
	       int *ret_node, int *ret_type)
{
  /* post->mmx[i][k]: posterior probability that posn i was emitted from node k's 
     match state */  
  int  max_k;    /* node index with highest posterior probability of emitting posn x */
  int  max_type; /* type of state in max_k node with max probability '0' for match, 
		    '1' for insert */
  int  max_sc;   /* score (log probability) from post matrix for max_k node max_type state type */
  int  k;        /* counter over nodes */
  if(x > j0 || x < i0)
    Die("ERROR in CP9NodeForPosn(), asking for position x: %d outside subseq bounds i0: %d j0: %d\n", x, i0, j0);

  if(post->mmx[x][0] > post->imx[x][0])
    {
      max_sc     = post->mmx[x][0];
      max_type   = 0; /* match */
    }
  else
    {
      max_sc     = post->imx[x][0];
      max_type   = 1; /* insert */
    }
  max_k    = 0; 

  for(k = 1; k <= hmm->M; k++)
    {
      if(post->mmx[x][k] > max_sc)
	{
	  max_k  = k;
	  max_sc = post->mmx[x][k];
	  max_type = 0; /* match */
	}
      if(post->imx[x][k] > max_sc)
	{
	  max_k  = k;
	  max_sc = post->imx[x][k];
	  max_type = 1; /* insert */
	}
    }
  if(max_type == 0)
    printf("MATCH | mx->mmx[%3d][%3d]: %9d | %8f\n", x, max_k, post->mmx[x][max_k], Score2Prob(post->mmx[x][max_k], 1.));
  else
    printf("INSERT | mx->imx[%3d][%3d]: %9d | %8f\n", x, max_k, post->imx[x][max_k], Score2Prob(post->imx[x][max_k], 1.));

  *ret_node = max_k;
  *ret_type = max_type;
  return;
}


/**********************************************************
 * Function:  StripWUSSGivenCC()
 * EPN 09.07.05
 *
 * Purpose:   Strips a secondary structure string in WUSS notation 
 *            of base pair information for specific match (consensus) columns.
 *            namely those before the first match column given by first_match,
 *            and after the last match column, given by last_match
 *            The msa->ss_cons secondary structure string is modified.
 *            
 *            Characters <([{  are converted to :   (left base of base pairs)
 *            Characters >)]}  are converted to :   (right base of base pairs)
 *            Characters _-,   are converted to :   (unpaired bases)
 *            Characters  .:~  are untouched        
 *            Pseudoknot characters are converted to : as well.
 *
 * Args:      msa         - the multiple sequence alignment
 *            dsq         - the sequences in the msa
 *            gapthresh   - the gap threshold for calling a match column
 *            first_match - first match column to keep structure for
 *            last_match  - last match column to keep structure for
 * Returns:   (void)
 */
void
StripWUSSGivenCC(MSA *msa, char **dsq, float gapthresh, int first_match, int last_match)
{
  int            *matassign;	/* 0..alen-1 array; 0=insert col, 1=match col */
  int gaps;
  char *s;
  int apos;
  int idx;
  int cc;
  int            *ct;		/* 0..alen-1 base pair partners array         */

  /* 1. Determine match/insert assignments
   *    matassign is 1..alen. Values are 1 if a match column, 0 if insert column.
   */
  matassign = MallocOrDie(sizeof(int) * (msa->alen+1));
  for (apos = 1; apos <= msa->alen; apos++)
    {
      for (gaps = 0, idx = 0; idx < msa->nseq; idx++)
	if (dsq[idx][apos] == DIGITAL_GAP) gaps++;
      matassign[apos] = ((double) gaps / (double) msa->nseq > gapthresh) ? 0 : 1;
    }

  /* 2. Determine a "ct" array, base-pairing partners for each position.
   *    Disallow/ignore pseudoknots. (That's what the FALSE flag does.)
   *    ct[] values give the index of a base pairing partner, or 0 for unpaired positions.
   *    Even though msa->ss_cons is in the 0..alen-1 coord system of msa, ct[]
   *    comes back in the 1..alen coord system of dsq.
   */
  if (! WUSS2ct(msa->ss_cons, msa->alen, FALSE, &ct))  
    Die("Consensus structure string is inconsistent"); 

  /* 3. Make sure the consensus structure "ct" is consistent with the match assignments.
   *    Wipe out all structure in insert columns; including the base-paired 
   *    partner of insert-assigned columns. 
   *    Also, remove structure outside of the consensus columns that 
   *    map to the HMM nodes first_match and last_match.
   */
  cc = 0;
  for (apos = 1; apos <= msa->alen; apos++)
    {
      if (! matassign[apos])
	{ 
	  if (ct[apos] != 0)  ct[ct[apos]] = 0;
	  ct[apos] = 0;
	}
      else /* matassign[apos] == 1 */
	{
	  cc++; 
	  if(cc < first_match || cc > last_match)
	  {
	    if (ct[apos] != 0)  ct[ct[apos]] = 0;
	    ct[apos] = 0;
	  }
	}
    }

  /* Next construct the new msa->ss_cons based on the ct array.
   * We should do this similar to display.c::CreateCMConsensus()
   * does it to get the fully formatted WUSS ([{<>}]) string but 
   * lazily we just do <> bps here.
   */
  for (apos = 1; apos <= msa->alen; apos++)
    {
      if      (ct[apos] == 0   ) msa->ss_cons[apos-1] = '.';
      else if (ct[apos]  > apos) msa->ss_cons[apos-1] = '<';
      else if (ct[apos]  < apos) msa->ss_cons[apos-1] = '>';
      else Die("ERROR: weird error in StripWUSSGivenCC\n");
    }

  /*
  apos = 1;
  cc   = 0;
  printf("first_match: %d | last_match: %d\n");
  for (s = msa->ss_cons; *s != '\0'; s++)
    {
      if(cc < first_match || cc > last_match)
	if ((*s != '~') && (*s != '.')) 
	  *s = ':';
    }
  */
  return;
}


/**************************************************************************
 * EPN 08.30.06
 * Function: map_orig2sub_cm()
 * derived loosely from hmmband.c::CP9_map_cm2hmm_and_hmm2cm()
 * Purpose:  Determine maps between a template CM and a sub CM, which probably
 *           has less structure (less MATP n_odes) than the template.
 * Args:    
 * CM_t *orig_cm       - the original template CM
 * CM_t *sub_cm        - the sub CM
 * int ***ret_orig2sub_smap - 2D state map from orig_cm (template) to sub_cm.
 *                       1st dimension - state index in orig_cm 
 *                       2nd D - 2 elements for up to 2 matching sub_cm states, .
 *                       Allocated here, responsibility of caller to free.
 * int ***ret_sub2orig_smap - 2D state map from orig_cm (template) to sub_cm.
 *                       1st dimension - state index in sub_cm (0..sub_cm->M-1)
 *                       2nd D - 2 elements for up to 2 matching orig_cm states, 
 *                       Allocated here, responsibility of caller to free.
 * int sub_start       - first consensus column in orig_cm that sub_cm models
 * int sub_end         -  last consensus column in orig_cm that sub_cm models

 * Returns: (void) 
 */
void
map_orig2sub_cm(CM_t *orig_cm, CM_t *sub_cm, int ***ret_orig2sub_smap, int ***ret_sub2orig_smap,
		int sub_start, int sub_end)
{

  int v_o; /* state counter for original CM */
  int v_s; /* state counter for sub CM */
  int is_left; /* TRUE if HMM node k maps to left half of CM node n */
  CMEmitMap_t *orig_emap;         /* consensus emit map for the original, template CM */
  CMEmitMap_t *sub_emap;         /* consensus emit map for the sub_cm */
  int         *o_node_cc_left; /* consensus column each CM node's left emission maps to
	                      * [0..(cm->nodes-1)], -1 if maps to n_o consensus column*/
  int               *o_node_cc_right;/* consensus column each CM node's right emission maps to
			            * [0..(cm->nodes-1)], -1 if maps to n_o consensus column*/
  int               *o_cc_node_map;  /* node that each consensus column maps to (is modelled by)
			            * [1..hmm_nmc] */
  int         *s_node_cc_left; /* consensus column each CM node's left emission maps to
	                      * [0..(cm->nodes-1)], -1 if maps to n_o consensus column*/
  int               *s_node_cc_right;/* consensus column each CM node's right emission maps to
			            * [0..(cm->nodes-1)], -1 if maps to n_o consensus column*/
  int               *s_cc_node_map;  /* sub_cm node that each consensus column maps to (is modelled by)
			            * [1..hmm_nmc] */
  int cc;
  int n_o;
  int n_s;
  char **nodetypes;
  char **sttypes;
  int *orig_cc2lins_map;
  int *orig_cc2rins_map;
  int **orig2sub_smap;
  int **sub2orig_smap;

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

  /* sanity check */
  if(sub_cm->M > orig_cm->M)
    Die("ERROR: sub_cm has more states than orig_cm in map_orig2sub_cm()\n");

  /* Get emitmap's for each CM */
  orig_emap = CreateEmitMap(orig_cm);
  sub_emap = CreateEmitMap(sub_cm);

  /* map the nodes of each CM to consensus column indices and vice versa */
  map_consensus_columns(orig_cm, orig_emap->clen, &o_node_cc_left, &o_node_cc_right,
			&o_cc_node_map, 0);
  map_consensus_columns(sub_cm,  sub_emap->clen,  &s_node_cc_left, &s_node_cc_right,
			&s_cc_node_map, 0);

  /* Allocate the maps */
  orig2sub_smap = MallocOrDie(sizeof(int *) * (orig_cm->M));
  for(v_o = 0; v_o < orig_cm->M; v_o++)
    orig2sub_smap[v_o] = MallocOrDie(sizeof(int) * 2);

  sub2orig_smap = MallocOrDie(sizeof(int *) * (sub_cm->M));
  for(v_s = 0; v_s < sub_cm->M; v_s++)
    sub2orig_smap[v_s] = MallocOrDie(sizeof(int) * 2);
    
  /* Initialize the maps */
  for(v_o = 0; v_o < orig_cm->M; v_o++)
    {
      orig2sub_smap[v_o][0] = -1;
      orig2sub_smap[v_o][1] = -1;
    }
  for(v_s = 0; v_s < sub_cm->M; v_s++)
    {
      sub2orig_smap[v_s][0] = -1;
      sub2orig_smap[v_s][1] = -1;
    }

  /* To make things easy, and avoid going traversing the emap multiple times
   * wastefully, we make maps of which nodes have an insert state that
   * inserts AFTER (in case of *lmap) and BEFORE (in case of *rmap)
   * each consensus node. 
   * orig_cc2lins_map[cc] = v, where orig_cm state v is an IL_st that
   *                           emits AFTER consensus column cc.
   * orig_cc2rins_map[cc] = v, where orig_cm state v is an IR_st that
   *                           emits BEFORE consensus column cc.
   * if no such state exists the value will be -1.
   */

  /* Allocate and initialize */
  orig_cc2lins_map = MallocOrDie(sizeof(int) * (orig_emap->clen + 2));
  orig_cc2rins_map = MallocOrDie(sizeof(int) * (orig_emap->clen + 2));
  for(cc = 0; cc <= orig_emap->clen+1; cc++)
    {
      orig_cc2lins_map[cc] = -1;
      orig_cc2rins_map[cc] = -1;
    }

  /* ROOT is special */
  orig_cc2lins_map[0] = 1; /* ROOT_IL */
  orig_cc2rins_map[orig_emap->clen+1] = 2; /* ROOT_IR */
  
  for(n_o = 0; n_o < orig_cm->nodes; n_o++)
    {
      switch (orig_cm->ndtype[n_o]) {
      case MATP_nd:
	orig_cc2lins_map[orig_emap->lpos[n_o]] = orig_cm->nodemap[n_o] + 4; /* MATP_IL */
	orig_cc2rins_map[orig_emap->rpos[n_o]] = orig_cm->nodemap[n_o] + 5; /* MATP_IR */
	break;
	
      case MATL_nd:
	orig_cc2lins_map[orig_emap->lpos[n_o]] = orig_cm->nodemap[n_o] + 2; /* MATL_IL */
	break;
	
      case MATR_nd:
	orig_cc2rins_map[orig_emap->rpos[n_o]] = orig_cm->nodemap[n_o] + 2; /* MATR_IR */
	break;
	
      case BEGR_nd:
	orig_cc2lins_map[orig_emap->lpos[n_o]] = orig_cm->nodemap[n_o] + 1; /* BEGR_IL */
	break;
	
      default: {} /*do nothing*/
      }
    }
  for(cc = 0; cc <= orig_emap->clen+1; cc++)
    {
      //printf("!cc: %d lins: %d rins: %d\n", cc, orig_cc2lins_map[cc], orig_cc2rins_map[cc]);
    }
  
  for(cc = sub_start; cc <= sub_end; cc++)
    {
      n_o = o_cc_node_map[cc];
      //printf("n_o: %d\n", n_o);
    }
  /* Using the *cc_node_maps, traverse the consensus columns left 
   * to right and determine state mappings for MATP, MATL and MATR nodes
   * first. For orig_cm MATL_nd's that map to sub_cm MATR_nd's and 
   * orig_cm MATL_nd's that map to sub_cm MATL_nd's we have to be
   * careful when mapping inserts (because IR's emit on right and
   * IL's emit on left), to make this easier we use the
   * orig_cc2lins_map and orig_cc2rins_map's.
   */
  for(cc = sub_start; cc <= sub_end; cc++)
    {
      printf("\n\nCC: %d\n\n", cc);
      n_o = o_cc_node_map[cc];
      printf("n_o: %d\n", n_o);
      n_s = s_cc_node_map[(cc-sub_start+1)]; /* note offset */
      if(orig_cm->ndtype[n_o] != MATP_nd &&
	 orig_cm->ndtype[n_o] != MATL_nd &&
	 orig_cm->ndtype[n_o] != MATR_nd)
	  Die("ERROR 1o in map_orig2sub_cm(), n_o: %d\n", n_o);
      if(sub_cm->ndtype[n_s] != MATP_nd &&
	 sub_cm->ndtype[n_s] != MATL_nd &&
	 sub_cm->ndtype[n_s] != MATR_nd)
	  Die("ERROR 1s in map_orig2sub_cm()\n");
      if(orig_cm->ndtype[n_o] == sub_cm->ndtype[n_s])
	{
	  /* same node type maps to same consensus column */
	  /* map all corresponding states */
	  v_o = orig_cm->nodemap[n_o];
	  v_s =  sub_cm->nodemap[n_s];
	  while(orig_cm->ndidx[v_o] == n_o)
	    {

	      if(sub_cm->ndidx[v_s] != n_s)
		Die("ERROR 2 in map_orig2sub_cm()\n");
	      map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);

	      /* if we're mapping insert states, we may have another
	       * insert that should map to this state */
	      if(sub_cm->sttype[v_s] == IL_st || sub_cm->sttype[v_s] == IR_st)
		{
		  if(orig_cc2lins_map[cc-1] != -1 && orig_cc2lins_map[cc-1] != v_o) 
		    {
		      /* orig_cc2lins_map[cc-1] is orig_cm left insert state (if any) that inserts BEFORE cc */
		      map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, orig_cc2lins_map[cc-1], v_s);
		      printf("DUAL INS mapping v_s:%d to v_o:%d\n", v_s, orig_cc2lins_map[cc-1]);
		    }
		  if(orig_cc2rins_map[cc] != -1 && orig_cc2rins_map[cc] != v_o) 
		    {
		      /* orig_cc2rins_map[cc] is orig_cm right insert state (if any) that inserts BEFORE cc */
		      map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, orig_cc2rins_map[cc], v_s);
		      printf("DUAL INS mapping v_s:%d to v_o:%d\n", v_s, orig_cc2rins_map[cc]);
		    }		
		}  
	      v_o++; 
	      v_s++;
	    }
	}
      /* n_o and n_s are of different node types. 
       * case 1: n_o is a MATP, n_s is a MATL 
       * case 2: n_o is a MATP, n_s is a MATR 
       * case 3: n_o is a MATL, n_s is a MATR
       * case 4: n_o is a MATR, n_s is a MATL
       * any other cases are impossible */
      else if(orig_cm->ndtype[n_o] == MATP_nd) /*case 1 or 2*/
	{
	  if(o_node_cc_left[o_cc_node_map[cc]] == cc)
	    is_left = TRUE;
	  else 
	    is_left = FALSE;

	  v_o = orig_cm->nodemap[n_o];
	  v_s =  sub_cm->nodemap[n_s];
	  if(sub_cm->ndtype[n_s] != MATL_nd && sub_cm->ndtype[n_s] != MATR_nd)
	    {
	      printf("orig_cm %4d %4d %4s %2s | sub_cm %4d %4d %4s %2s\n", n_o, v_o, nodetypes[orig_cm->ndtype[n_o]], sttypes[orig_cm->sttype[v_o]], n_s, v_s, nodetypes[sub_cm->ndtype[n_s]], sttypes[sub_cm->sttype[v_s]]);
	      Die("ERROR 3 in map_orig2sub_cm()\n");
	    }
	  
	  if(sub_cm->ndtype[n_s] == MATL_nd) /* case 1 */
	    {
	      if(is_left)
		{
		  printf("n_o MATP n_s MATL\n\n");
		  /* Case by case (uchh..) 
		   * v_o is MATP_MP 
		   * v_s is MATL_ML */
		  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);
		  v_o++; /* v_o is MATP_ML */
		  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);
		  
		  v_s++; /* v_s is MATL_D */
		  v_o++; /* v_o is MATP_MR */
		  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);
		  v_o++; /* v_o is MATP_D */
		  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);
		  
		  v_s++; /* v_s is MATL_IL */
		  v_o++; /* v_o is MATP_IL */
		  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);
		}
	      else 
		{
		  printf("n_o MATP n_s MATL\n\n");
		  /* Case by case (uchh..) 
		   * v_o is MATP_MP 
		   * v_s is MATL_ML */
		  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);
		  v_o++; /* v_o is MATP_ML */
		  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);
		  
		  v_s++; /* v_s is MATL_D */
		  v_o++; /* v_o is MATP_MR */
		  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);
		  v_o++; /* v_o is MATP_D */
		  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);
		  
		  v_s++; /* v_s is MATL_IL */
		  v_o++; /* v_o is MATP_IL */
		  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);
		}
		
	    }
	  else if(sub_cm->ndtype[n_s] == MATR_nd) /* case 2 */
	    {
	      /* Case by case (uchh..) */
	      /* v_o is MATP_MP */
	      /* v_s is MATR_MR */
	      if(orig_cm->sttype[v_o] != MP_st)
		Die("WHOA, state v_o: %d is not a MATP_MP\n");
	      if(sub_cm->sttype[v_s] != MR_st)
		Die("WHOA, state v_s: %d is not a MATR_MR\n");

	      printf("mapping MATP_MP orig v_o: %d to MATR_R sub v_s: %d\n", v_o, v_s);
	      map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);
	      v_o += 2; /* v_o is MATP_MR */
	      if(orig_cm->sttype[v_o] != MR_st)
		Die("WHOA, state v_o: %d is not a MATP_MR\n");
	      map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);

	      v_s++; /* v_s is MATR_D */
	      if(sub_cm->sttype[v_s] != D_st)
		Die("WHOA, state v_s: %d is not a MATR_D\n");
	      v_o--; /* v_o is MATP_ML */
	      if(orig_cm->sttype[v_o] != ML_st)
		Die("WHOA, state v_o: %d is not a MATP_ML\n");
	      map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);
	      v_o += 2; /* v_o is MATP_D */
	      if(orig_cm->sttype[v_o] != D_st)
		Die("WHOA, state v_o: %d is not a MATP_D\n");
	      map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);
	      
	      v_s++; /* v_s is MATR_IR */
	      if(sub_cm->sttype[v_s] != IR_st)
		Die("WHOA, state v_s: %d is not a MATR_IR\n");
	      v_o += 2; /* v_o is MATP_IR */
	      if(orig_cm->sttype[v_o] != IR_st)
		Die("WHOA, state v_o: %d is not a MATP_IR\n");
	      map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);
	    }
	  else 
	    {
	      Die("ERROR for consensus column %d in map_orig2sub_cm() orig node is MATP, but sub node is not MATP, MATL or MATR\n", cc);
	    }
	}
      else if(orig_cm->ndtype[n_o] == MATL_nd) /*case 3: n_o is MATL n_s is MATR */
	{
	  v_o = orig_cm->nodemap[n_o]; /* v_o is MATL_ML */
	  v_s =  sub_cm->nodemap[n_s]; /* v_s is MATR_MR */
	  /* sanity check */
	  if(sub_cm->ndtype[n_s] != MATR_nd)
	    Die("ERROR 4 in map_orig2sub_cm()\n");
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);

	  v_o++; /* v_o is MATL_D */
	  v_s++; /* v_s is MATR_D */
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);

	  v_o = orig_cc2lins_map[cc-1]; /* orig_cm left insert state (if any) that inserts BEFORE cc */
	  v_s++; /* v_s is MATR_IR, potentially two orig_cm states map to v_s */
	  if(v_o != -1)
	    map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);
	  v_o = orig_cc2rins_map[cc]; /* orig_cm right insert state (if any) that inserts BEFORE cc */
	  if(v_o != -1)
	    map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);

	  /* sanity check */
	  if(orig_cc2lins_map[cc-1] == -1 && orig_cc2rins_map[cc] == -1)
	    Die("ERROR 5 in map_orig2sub_cm() no insert state in orig_cm emits before cc:%d\n", cc);	    
	}
      else if(orig_cm->ndtype[n_o] == MATR_nd) /*case 4: n_o is MATR n_s is MATL */
	{
	  v_o = orig_cm->nodemap[n_o]; /* v_o is MATR_MR */
	  v_s =  sub_cm->nodemap[n_s]; /* v_s is MATL_ML */
	  /* sanity check */
	  if(sub_cm->ndtype[n_s] != MATL_nd)
	    Die("ERROR 6 in map_orig2sub_cm()\n");
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);

	  v_o++; /* v_o is MATR_D */
	  v_s++; /* v_s is MATL_D */
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);

	  v_o = orig_cc2lins_map[cc]; /* orig_cm left insert state (if any) that inserts AFTER cc */
	  v_s++; /* v_s is MATL_IL, potentially two orig_cm states map to v_s */
	  if(v_o != -1)
	    map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);
	  v_o = orig_cc2rins_map[cc+1]; /* orig_cm right insert state (if any) that inserts AFTER cc */
	  if(v_o != -1)
	    map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);

	  /* sanity check */
	  if(orig_cc2lins_map[cc] == -1 && orig_cc2rins_map[cc+1] == -1)
	    Die("ERROR 6 in map_orig2sub_cm() no insert state in orig_cm emits after cc:%d\n", cc);	    
	}
      else
	Die("ERROR 7 in map_orig2sub_cm() n_o: %d and n_s: %d, cases 1-4 don't apply.\n");
    }     
  
  /* SPECIAL CASE: Map the ROOT_S's to each other */
  v_o = 0; /* ROOT_S */
  v_s = 0; /* ROOT_S */
  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);

  /* Map ROOT_IL, ROOT_IR, and BEGR_IL's in an ugly manner, by looking for the
   * insert states in the emap that map to them. 
   */

  /* Now map ROOT_IL, ROOT_IR and all BEGR_IL's, we could have either IL or IR
   * states in orig_cm that map to these guys, so we need to use
   * both the lmap and rmap. In some cases there may be two orig_cm insert states that
   * map to a single insert state in sub_cm (due to an ambiguity in the Infernal
   * grammar).
   */

  /* ROOT_IL */
  v_s = 1;
  v_o = orig_cc2lins_map[sub_start-1]; /* we could have a left emitter that inserts before 
					* consensus column sub_start */
  if(v_o != -1)
    map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);

  v_o = orig_cc2rins_map[sub_start]; /* we could have a right emitter that inserts before 
					    * consensus column sub_start */
  if(v_o != -1)
    map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);

  /* ROOT_IR */
  v_s = 2;
  v_o = orig_cc2lins_map[sub_end]; /* we could have a left emitter that inserts after
					* consensus column sub_end */
  if(v_o != -1)
    map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);

  v_o = orig_cc2rins_map[sub_end+1]; /* we could have a right emitter that inserts after
					  * consensus column sub_end */
  if(v_o != -1)
    map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);

  /* BEGR_ILs */
  for(n_s = 0; n_s < sub_cm->nodes; n_s++)
    {
      if(sub_cm->ndtype[n_s] == BEGR_nd)
	{
	  printf("trying to fill in BEGR_IL for node: %d\n", n_s);
	  v_s = sub_cm->nodemap[n_s] + 1;
	  printf("checking lpos (%d + %d -1): %d\n", sub_emap->lpos[n_s], sub_start, (sub_emap->lpos[n_s] + sub_start-1));
	  v_o = orig_cc2lins_map[(sub_emap->lpos[n_s] + sub_start - 1)];
	  if(v_o != -1)
	    map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);

	  printf("checking rpos (%d + %d): %d\n", sub_emap->lpos[n_s], sub_start, (sub_emap->lpos[n_s] + sub_start));
	  v_o = orig_cc2rins_map[(sub_emap->lpos[n_s] + sub_start)];
	  printf("\t\tv_o: %d | v_s: %d\n", v_o, v_s);
	  if(v_o != -1)
	    map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);
	}
    }
  
  /* NOTE: We ignore mapping B states, E states and S states (except ROOT_S). We'll handle these guys 
   * specially when we fill in transition probabilities. The reason we don't map them is that
   * I think its impossible to do it robustly. I saw cases with SSU where there were BIF nodes in the
   * sub_cm for which I couldn't figure out which BIF nodes (and correspondingly BEGL and BEGR nodes)
   * in the orig_cm they should map to. Regardless, even if there is a pretty, nice way I'm abandoning
   * an attempt to find it for a crude way - we will handle the transitions involving these guys
   * specially, without a need for a mapping between CMs.
   */

  /* print sub2orig_smap, checking consistency with orig2sub_smap along the way. */
  for(v_s = 0; v_s < sub_cm->M; v_s++)
    {
      n_s = sub_cm->ndidx[v_s];
      v_o = sub2orig_smap[v_s][0];

      if(sub_cm->sttype[v_s] == E_st) 	
	{
	  printf("sub v:%4d   END\n", v_s);
	  continue;
	}
      if(sub_cm->sttype[v_s] == S_st) 	
	{
	  printf("sub v:%4d   START\n", v_s);
	  continue;
	}
      if(sub_cm->sttype[v_s] == B_st) 	
	{
	  printf("sub v:%4d   BIF_B\n", v_s);
	  continue;
	}
      if(v_o == -1 && sub_cm->sttype[v_s] != E_st)
	Die("ERROR sub_cm state: %d type: %s node type: %s doesn't map to any state in orig_cm\n", v_s, sttypes[sub_cm->sttype[v_s]], nodetypes[sub_cm->ndtype[n_s]]);

      n_o = orig_cm->ndidx[v_o];
      
      printf("sub v:%4d(%4d) %6s%6s | orig v:%4d(%4d) %6s%6s\n", v_s, n_s, nodetypes[sub_cm->ndtype[n_s]], sttypes[sub_cm->sttype[v_s]], v_o, n_o, nodetypes[orig_cm->ndtype[n_o]], sttypes[orig_cm->sttype[v_o]]);

      /* check to make sure orig2sub_smap is consistent */
      if(orig2sub_smap[v_o][0] != v_s && orig2sub_smap[v_o][1] != v_s)
	{
	  Die("ERROR inconsistency; neither orig2sub_smap[%d][0] and [1] is %d\n", v_o, v_s);
	}

      v_o = sub2orig_smap[v_s][1];
      if(v_o != -1)
	{
	  n_o = orig_cm->ndidx[v_o];
	  printf("                            | orig v:%4d(%4d) %6s%6s\n", v_o, n_o, nodetypes[orig_cm->ndtype[n_o]], sttypes[orig_cm->sttype[v_o]]);
	  /* check to make sure orig2sub_smap is consistent */
	  if(orig2sub_smap[v_o][0] != v_s && orig2sub_smap[v_o][1] != v_s)
	    {
	      Die("ERROR inconsistency; neither orig2sub_smap[%d][0] and [1] is %d\n", v_o, v_s);
	    }
	}
    }

  *ret_orig2sub_smap = orig2sub_smap;
  *ret_sub2orig_smap = sub2orig_smap;

  return;
}


/**************************************************************************
 * EPN 08.30.06
 * Function: map_orig2sub_cm2()
 * 
 * Purpose:  Determine maps between a template CM and a sub CM, which probably
 *           has less structure (less MATP nodes) than the template. Do this
 *           a bit indirectly, first get maps from each CM to a CP9 HMM,
 *           then use these maps to map the CMs to each other.
 * Args:    
 * CM_t *orig_cm       - the original template CM
 * CM_t *sub_cm        - the sub CM
 * int ***ret_orig2sub_smap - 2D state map from orig_cm (template) to sub_cm.
 *                       1st dimension - state index in orig_cm 
 *                       2nd D - 2 elements for up to 2 matching sub_cm states, .
 *                       Allocated here, responsibility of caller to free.
 * int ***ret_sub2orig_smap - 2D state map from orig_cm (template) to sub_cm.
 *                       1st dimension - state index in sub_cm (0..sub_cm->M-1)
 *                       2nd D - 2 elements for up to 2 matching orig_cm states, 
 *                       Allocated here, responsibility of caller to free.
 * int sub_start       - first consensus column in orig_cm that sub_cm models
 * int sub_end         -  last consensus column in orig_cm that sub_cm models

 * Returns: (void) 
 */
void
map_orig2sub_cm2(CM_t *orig_cm, CM_t *sub_cm, int ***ret_orig2sub_smap, int ***ret_sub2orig_smap,
		 int sub_start, int sub_end)
{
  int **cs2hn_map_o;
  int **cs2hs_map_o;
  int ***hns2cs_map_o;
  int **cs2hn_map_s;
  int **cs2hs_map_s;
  int ***hns2cs_map_s;
  struct cplan9_s  *cp9hmm_s;        /* constructed CM p9 HMM from the sub_cm */
  struct cplan9_s  *cp9hmm_o;        /* constructed CM p9 HMM from the orig_cm */

  int         *node_cc_left_o; /* consensus column each CM node's left emission maps to
	                      * [0..(cm->nodes-1)], -1 if maps to n_o consensus column*/
  int         *node_cc_right_o;/* consensus column each CM node's right emission maps to
			            * [0..(cm->nodes-1)], -1 if maps to n_o consensus column*/
  int         *cc_node_map_o;  /* node that each consensus column maps to (is modelled by)
			            * [1..hmm_nmc] */
  int         *node_cc_left_s; /* consensus column each CM node's left emission maps to
	                      * [0..(cm->nodes-1)], -1 if maps to n_o consensus column*/
  int          *node_cc_right_s;/* consensus column each CM node's right emission maps to
			            * [0..(cm->nodes-1)], -1 if maps to n_o consensus column*/
  int         *cc_node_map_s;  /* sub_cm node that each consensus column maps to (is modelled by)
			            * [1..hmm_nmc] */
  int k_s;  /* HMM node counter */
  int ks_s; /* HMM state counter (0(Match) 1(insert) or 2(delete)*/
  int k_o;  /* HMM node counter */
  int v_o;
  int v_s;
  int v; /* state index in CM */
  CMEmitMap_t *orig_emap;         /* consensus emit map for the original, template CM */
  CMEmitMap_t *sub_emap;         /* consensus emit map for the sub_cm */
  int cc;
  int n_o;
  int n_s;
  char **nodetypes;
  char **sttypes;
  int *orig_cc2lins_map;
  int *orig_cc2rins_map;
  int **orig2sub_smap;
  int **sub2orig_smap;
  int cp9_M;
  int sub_k;
  int orig_k;
  int sub_nd;
  int orig_nd;
  int is_right;
  int easy_stretch;

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

  nodetypes = malloc(sizeof(char *) * 8);
  nodetypes[0] = "BIF";
  nodetypes[1] = "MATP";
  nodetypes[2] = "MATL";
  nodetypes[3] = "MATR";
  nodetypes[4] = "BEGL";
  nodetypes[5] = "BEGR";
  nodetypes[6] = "ROOT";
  nodetypes[7] = "END";


  /* sanity check */
  if(sub_cm->M > orig_cm->M)
    Die("ERROR: sub_cm has more states than orig_cm in map_orig2sub_cm()\n");

  /* Get emitmap's for each CM */
  orig_emap = CreateEmitMap(orig_cm);
  sub_emap = CreateEmitMap(sub_cm);

  /* To make things easy, and avoid going traversing the emap multiple times
   * wastefully, we make maps of which nodes have an insert state that
   * inserts AFTER (in case of *lmap) and BEFORE (in case of *rmap)
   * each consensus node. 
   * orig_cc2lins_map[cc] = v, where orig_cm state v is an IL_st that
   *                           emits AFTER consensus column cc.
   * orig_cc2rins_map[cc] = v, where orig_cm state v is an IR_st that
   *                           emits BEFORE consensus column cc.
   * if no such state exists the value will be -1.
   */

  /* Allocate and initialize */
  orig_cc2lins_map = MallocOrDie(sizeof(int) * (orig_emap->clen + 2));
  orig_cc2rins_map = MallocOrDie(sizeof(int) * (orig_emap->clen + 2));
  for(cc = 0; cc <= orig_emap->clen+1; cc++)
    {
      orig_cc2lins_map[cc] = -1;
      orig_cc2rins_map[cc] = -1;
    }
  /* fill in the map */
  /* ROOT is special */
  orig_cc2lins_map[0] = 1; /* ROOT_IL */
  orig_cc2rins_map[orig_emap->clen+1] = 2; /* ROOT_IR */
  for(n_o = 0; n_o < orig_cm->nodes; n_o++)
    {
      switch (orig_cm->ndtype[n_o]) {
      case MATP_nd:
	orig_cc2lins_map[orig_emap->lpos[n_o]] = orig_cm->nodemap[n_o] + 4; /* MATP_IL */
	orig_cc2rins_map[orig_emap->rpos[n_o]] = orig_cm->nodemap[n_o] + 5; /* MATP_IR */
	break;
	
      case MATL_nd:
	orig_cc2lins_map[orig_emap->lpos[n_o]] = orig_cm->nodemap[n_o] + 2; /* MATL_IL */
	break;
	
      case MATR_nd:
	orig_cc2rins_map[orig_emap->rpos[n_o]] = orig_cm->nodemap[n_o] + 2; /* MATR_IR */
	break;
	
      case BEGR_nd:
	orig_cc2lins_map[orig_emap->lpos[n_o]] = orig_cm->nodemap[n_o] + 1; /* BEGR_IL */
	break;
	
      default: {} /*do nothing*/
      }
    }
  for(cc = 0; cc <= orig_emap->clen+1; cc++)
    {
      printf("!cc: %d lins: %d rins: %d\n", cc, orig_cc2lins_map[cc], orig_cc2rins_map[cc]);
    }

  cp9_M = 0;
  for(v = 0; v <= orig_cm->M; v++)
    {
      if(orig_cm->stid[v] ==  MATP_MP)
	cp9_M += 2;
      else if(orig_cm->stid[v] == MATL_ML || orig_cm->stid[v] == MATR_MR)
	cp9_M++;
    }
  cp9hmm_o = AllocCPlan9(cp9_M);
  ZeroCPlan9(cp9hmm_o);
  CPlan9SetNullModel(cp9hmm_o, orig_cm->null, 1.0); /* set p1 = 1.0 which corresponds to the CM */

  cp9hmm_s = AllocCPlan9((sub_end-sub_start+1));
  ZeroCPlan9(cp9hmm_s);
  CPlan9SetNullModel(cp9hmm_s, sub_cm->null, 1.0); /* set p1 = 1.0 which corresponds to the CM */

  /* map the nodes of each CM to consensus column indices and vice versa */
  map_consensus_columns(orig_cm, cp9hmm_o->M, &node_cc_left_o, &node_cc_right_o,
			&cc_node_map_o, 0);
  map_consensus_columns(sub_cm,  cp9hmm_s->M, &node_cc_left_s, &node_cc_right_s,
			&cc_node_map_s, 0);
  /* map the CM states to CP9 states and nodes and vice versa */
  CP9_map_cm2hmm_and_hmm2cm(orig_cm, cp9hmm_o, node_cc_left_o, node_cc_right_o, 
			    cc_node_map_o, &cs2hn_map_o, &cs2hs_map_o, 
			    &hns2cs_map_o, 0);
  CP9_map_cm2hmm_and_hmm2cm(sub_cm, cp9hmm_s, node_cc_left_s, node_cc_right_s, 
			    cc_node_map_s, &cs2hn_map_s, &cs2hs_map_s, 
			    &hns2cs_map_s, 0);

  /* Allocate the maps */
  orig2sub_smap = MallocOrDie(sizeof(int *) * (orig_cm->M));
  for(v_o = 0; v_o < orig_cm->M; v_o++)
    orig2sub_smap[v_o] = MallocOrDie(sizeof(int) * 2);

  sub2orig_smap = MallocOrDie(sizeof(int *) * (sub_cm->M));
  for(v_s = 0; v_s < sub_cm->M; v_s++)
    sub2orig_smap[v_s] = MallocOrDie(sizeof(int) * 2);
    
  /* Initialize the maps */
  for(v_o = 0; v_o < orig_cm->M; v_o++)
    {
      orig2sub_smap[v_o][0] = -1;
      orig2sub_smap[v_o][1] = -1;
    }
  for(v_s = 0; v_s < sub_cm->M; v_s++)
    {
      sub2orig_smap[v_s][0] = -1;
      sub2orig_smap[v_s][1] = -1;
    }

  /* step through the states of sub_cm, filling in maps */
  /* ROOT is special: */
  v_s = 0; /* ROOT_S */
  v_o = 0; /* ROOT_S */
  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);

  /* ROOT_IL inserts before orig cc sub_start */
  k_s = 1; /* insert */
  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
			 hns2cs_map_o[sub_start-1][k_s][0],
			 hns2cs_map_s[0][k_s][0]);     
  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
			 hns2cs_map_o[sub_start-1][k_s][0],
			 hns2cs_map_s[0][k_s][1]);     
  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
			 hns2cs_map_o[sub_start-1][k_s][1],
			 hns2cs_map_s[0][k_s][0]);     
  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
			 hns2cs_map_o[sub_start-1][k_s][1],
			 hns2cs_map_s[0][k_s][1]);     

  /* ROOT_IR inserts after orig cc sub_end */
  k_s = 1; /* insert */
  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
			 hns2cs_map_o[sub_end][k_s][0],
			 hns2cs_map_s[sub_end-sub_start+1][k_s][0]);     
  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
			 hns2cs_map_o[sub_end][k_s][0],
			 hns2cs_map_s[sub_end-sub_start+1][k_s][1]);     
  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
			 hns2cs_map_o[sub_end][k_s][1],
			 hns2cs_map_s[sub_end-sub_start+1][k_s][0]);     
  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
			 hns2cs_map_o[sub_end][k_s][1],
			 hns2cs_map_s[sub_end-sub_start+1][k_s][1]);     

  /* BEGR_ILs */

  for(n_s = 0; n_s < sub_cm->nodes; n_s++)
    {
      if(sub_cm->ndtype[n_s] == BEGR_nd)
	{
	  printf("trying to fill in BEGR_IL for node: %d\n", n_s);
	  v_s = sub_cm->nodemap[n_s] + 1;
	  printf("checking lpos (%d + %d -1): %d\n", sub_emap->lpos[n_s], sub_start, (sub_emap->lpos[n_s] + sub_start-1));
	  v_o = orig_cc2lins_map[(sub_emap->lpos[n_s] + sub_start - 1)];
	  if(v_o != -1)
	    map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);

	  printf("checking rpos (%d + %d): %d\n", sub_emap->lpos[n_s], sub_start, (sub_emap->lpos[n_s] + sub_start));
	  v_o = orig_cc2rins_map[(sub_emap->lpos[n_s] + sub_start)];
	  printf("\t\tv_o: %d | v_s: %d\n", v_o, v_s);
	  if(v_o != -1)
	    map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);
	}
    }

  for(sub_k = 1; sub_k <= cp9hmm_s->M; sub_k++)
    {
      orig_k = sub_k + sub_start - 1;
      
      printf("sub_k: %d | orig_k: %d\n", sub_k, orig_k);
      sub_nd  = cc_node_map_s[sub_k];
      orig_nd = cc_node_map_o[orig_k];
      printf("sub_nd: %d (%s) orig_nd: %d (%s)\n", sub_nd, nodetypes[sub_cm->ndtype[sub_nd]],
	     orig_nd, nodetypes[orig_cm->ndtype[orig_nd]]);

      /* first check for the easy case: when same node type
       * maps to the same consensus column.
       */
      if(orig_cm->ndtype[orig_nd] == sub_cm->ndtype[sub_nd])
	{
	  easy_stretch = TRUE;

	  if(node_cc_left_o[orig_nd] == orig_k)
	    is_right = FALSE;
	  else if(node_cc_right_o[orig_nd] == orig_k)
	    is_right = TRUE;
	  else
	    Die("You don't understand CMs 0\n");

	  printf("EASY CASE\n");
	  /* map all corresponding states */
	  v_o = orig_cm->nodemap[orig_nd];
	  v_s =  sub_cm->nodemap[sub_nd];
	  while(orig_cm->ndidx[v_o] == orig_nd)
	    {

	      if(sub_cm->ndidx[v_s] != sub_nd)
		Die("ERROR 2 in map_orig2sub_cm()\n");
	      map_orig2sub_cm_helper(orig_cm, sub_cm,
				     orig2sub_smap, sub2orig_smap, v_o, v_s);

	      /* if we're mapping insert states in NON-MATP nodes, 
	       * we may have another orig_cm insert that should map to this state */
	      if((sub_cm->ndtype[sub_cm->ndidx[v_s]] != MATP_nd) &&
		 (sub_cm->sttype[v_s] == IL_st || sub_cm->sttype[v_s] == IR_st))
		{
		  k_s = 1; /* insert */		 
		  printf("!*! sub_k: %d | orig_k: %d\n", sub_k, orig_k);
		  if(hns2cs_map_o[orig_k][k_s][0] != -1 && hns2cs_map_o[(orig_k-is_right)][k_s][1] != v_o)
		    map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
					   hns2cs_map_o[orig_k-is_right][k_s][0], v_s);
		  if(hns2cs_map_o[orig_k][k_s][1] != -1 && hns2cs_map_o[(orig_k-is_right)][k_s][1] != v_o)
		    map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
					   hns2cs_map_o[orig_k-is_right][k_s][1], v_s);
		}
	      v_o++; 
	      v_s++;
	    }
	}
      /* end easy case check */

      else
	{
	  if(easy_stretch && (is_right && sub_cm->ndtype[sub_nd] == MATR_nd))
	    {
	      /* previous column was an 'easy case' where we had an identical node
	       * in the sub_cm to the orig_cm, in this case we need to check for
	       * a special situation here: when prev column was modelled by a right
	       * half of an orig_cm node, and the current sub_cm node is a MATR,
	       * in this case we may have missed a mapping from the MATR_IR of this
	       * state to the orig_cm state that inserts AFTER the prev_cc.
	       */
	      printf("SPECIAL EASY CASE STRETCH NOW FALSE, MATR CASE\n");
	      k_s = 1; /*insert*/
	      map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				     hns2cs_map_o[orig_k-1][k_s][0],
				     hns2cs_map_s[sub_k-1][k_s][0]);     
	      map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				     hns2cs_map_o[orig_k-1][k_s][0],
				     hns2cs_map_s[sub_k-1][k_s][1]);     
	      map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				     hns2cs_map_o[orig_k-1][k_s][1],
				     hns2cs_map_s[sub_k-1][k_s][0]);     
	      map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				     hns2cs_map_o[orig_k-1][k_s][1],
				     hns2cs_map_s[sub_k-1][k_s][1]);     
	    }

	  easy_stretch = FALSE;

	  printf("sub_k: %d | orig_k: %d\n", sub_k, orig_k);
	  k_s = 0; /* match */
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				 hns2cs_map_o[orig_k][k_s][0],
				 hns2cs_map_s[sub_k][k_s][0]);     
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				 hns2cs_map_o[orig_k][k_s][0],
				 hns2cs_map_s[sub_k][k_s][1]);     
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				 hns2cs_map_o[orig_k][k_s][1],
				 hns2cs_map_s[sub_k][k_s][0]);     
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				 hns2cs_map_o[orig_k][k_s][1],
				 hns2cs_map_s[sub_k][k_s][1]);     
	  
	  k_s = 1; /* insert */
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				 hns2cs_map_o[orig_k][k_s][0],
				 hns2cs_map_s[sub_k][k_s][0]);     
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				 hns2cs_map_o[orig_k][k_s][0],
				 hns2cs_map_s[sub_k][k_s][1]);     
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				 hns2cs_map_o[orig_k][k_s][1],
				 hns2cs_map_s[sub_k][k_s][0]);     
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				 hns2cs_map_o[orig_k][k_s][1],
				 hns2cs_map_s[sub_k][k_s][1]);     
	  
	  
	  k_s = 2; /* delete */
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				 hns2cs_map_o[orig_k][k_s][0],
				 hns2cs_map_s[sub_k][k_s][0]);     
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				 hns2cs_map_o[orig_k][k_s][0],
				 hns2cs_map_s[sub_k][k_s][1]);     
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				 hns2cs_map_o[orig_k][k_s][1],
				 hns2cs_map_s[sub_k][k_s][0]);     
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				 hns2cs_map_o[orig_k][k_s][1],
				 hns2cs_map_s[sub_k][k_s][1]);    
	}
    }
  /* NOTE: We ignore mapping B states, E states and S states (except ROOT_S). We'll handle these guys 
   * specially when we fill in transition probabilities. The reason we don't map them is that
   * I think its impossible to do it robustly. I saw cases with SSU where there were BIF nodes in the
   * sub_cm for which I couldn't figure out which BIF nodes (and correspondingly BEGL and BEGR nodes)
   * in the orig_cm they should map to. Regardless, even if there is a pretty, nice way I'm abandoning
   * an attempt to find it for a crude way - we will handle the transitions involving these guys
   * specially, without a need for a mapping between CMs.
   */
  
  /* print sub2orig_smap, checking consistency with orig2sub_smap along the way. */
  printf("\n\n\n KACHOW! MAP\n\n\n");
  for(v_s = 0; v_s < sub_cm->M; v_s++)
    {
      n_s = sub_cm->ndidx[v_s];
      v_o = sub2orig_smap[v_s][0];
      /*      printf("v_s: %d\n", v_s);
	      printf("n_s: %d\n", n_s);
	      printf("v_o: %d\n", v_o);*/
      if(sub_cm->sttype[v_s] == E_st) 	
	{
	  printf("sub v:%4d   END\n", v_s);
	  continue;
	}
      if(sub_cm->sttype[v_s] == S_st) 	
	{
	  printf("sub v:%4d   START\n", v_s);
	  continue;
	}
      if(sub_cm->sttype[v_s] == B_st) 	
	{
	  printf("sub v:%4d   BIF_B\n", v_s);
	  continue;
	}
      if(v_o == -1 && sub_cm->sttype[v_s] != E_st)
	Die("ERROR sub_cm state: %d type: %s node type: %s doesn't map to any state in orig_cm\n", v_s, sttypes[sub_cm->sttype[v_s]], nodetypes[sub_cm->ndtype[n_s]]);

      n_o = orig_cm->ndidx[v_o];
      
      printf("sub v:%4d(%4d) %6s%6s | orig v:%4d(%4d) %6s%6s\n", v_s, n_s, nodetypes[sub_cm->ndtype[n_s]], sttypes[sub_cm->sttype[v_s]], v_o, n_o, nodetypes[orig_cm->ndtype[n_o]], sttypes[orig_cm->sttype[v_o]]);

      /* check to make sure orig2sub_smap is consistent */
      if(orig2sub_smap[v_o][0] != v_s && orig2sub_smap[v_o][1] != v_s)
	{
	  Die("ERROR inconsistency; neither orig2sub_smap[%d][0] and [1] is %d\n", v_o, v_s);
	}

      v_o = sub2orig_smap[v_s][1];
      if(v_o != -1)
	{
	  n_o = orig_cm->ndidx[v_o];
	  printf("                              | orig v:%4d(%4d) %6s%6s\n", v_o, n_o, nodetypes[orig_cm->ndtype[n_o]], sttypes[orig_cm->sttype[v_o]]);
	  /* check to make sure orig2sub_smap is consistent */
	  if(orig2sub_smap[v_o][0] != v_s && orig2sub_smap[v_o][1] != v_s)
	    {
	      Die("ERROR inconsistency; neither orig2sub_smap[%d][0] and [1] is %d\n", v_o, v_s);
	    }
	}
    }


  /* Clean up and return */
  free(node_cc_left_o);
  free(node_cc_right_o);
  free(cc_node_map_o);
  free(node_cc_left_s);
  free(node_cc_right_s);
  free(cc_node_map_s);
  for(v = 0; v <= orig_cm->M; v++)
    {
      free(cs2hn_map_o[v]);
      free(cs2hs_map_o[v]);
    }
  free(cs2hn_map_o);
  for(v = 0; v <= sub_cm->M; v++)
    {
      free(cs2hn_map_s[v]);
      free(cs2hs_map_s[v]);
    }
  free(cs2hn_map_s);

  for(k_o = 0; k_o <= cp9hmm_o->M; k_o++)
    {
      for(ks_s = 0; ks_s < 3; ks_s++)
	free(hns2cs_map_o[k_o][ks_s]);
      free(hns2cs_map_o[k_o]);
    }
  free(hns2cs_map_o);
  FreeCPlan9(cp9hmm_o);

  for(k_s = 0; k_s <= cp9hmm_s->M; k_s++)
    {
      for(ks_s = 0; ks_s < 3; ks_s++)
	free(hns2cs_map_s[k_s][ks_s]);
      free(hns2cs_map_s[k_s]);
    }
  free(hns2cs_map_s);
  FreeCPlan9(cp9hmm_s);



  *ret_orig2sub_smap = orig2sub_smap;
  *ret_sub2orig_smap = sub2orig_smap;

  return;
}


/**************************************************************************
 * EPN 08.30.06
 * map_orig2sub_cm_helper()
 *
 * helper function for map_orig2sub_cm(), 
 *
 * Purpose:  Fill in specific part of the map, given orig_v (orig_cm state),
 *           sub_v (sub_cm state)
 * Args:    
 * CM_t *orig_cm
 * CM_t *sub_cm
 * int **orig2sub_smap
 * int **sub2orig_smap
 * int   orig_v; 
 * int   sub_v; 
 * Returns: 1 if a new mapping was made, 0 if not 
 */
static int
map_orig2sub_cm_helper(CM_t *orig_cm, CM_t *sub_cm, int **orig2sub_smap, int **sub2orig_smap, int orig_v, int sub_v)
{
  int sub_nd;
  int orig_nd;
  int is_insert;
  int temp;

  //printf("\nin helper: orig_v: %d sub_v: %d\n", orig_v, sub_v);

  if(orig_v == -1 || sub_v == -1)
    return 0;

  /* check to see if we already have this mapping */
  if(orig2sub_smap[orig_v][0] == sub_v || orig2sub_smap[orig_v][1] == sub_v)
    {
      //printf("already have mapping, returning\n");
      return 0;
    }

  orig_nd = orig_cm->ndidx[orig_v];
  sub_nd  =  sub_cm->ndidx[sub_v];

  is_insert = FALSE;
  if(sub_cm->sttype[sub_v] == IL_st || sub_cm->sttype[sub_v] == IR_st)
    {
      if(!(orig_cm->sttype[orig_v] == IL_st || orig_cm->sttype[orig_v] == IR_st))
	Die("INSERT error in map_orig2sub_cm_helper()\n");
      is_insert = TRUE;
      
      /* Make sure that neither orig_v or sub_v are detached insert states,
       * if they are we return b/c they're irrelevant and we don't want that 
       * information in the maps. */
      if(orig_cm->sttype[(orig_v+1)] == E_st)
	return 0;
      //if(sub_cm->sttype[(sub_v+1)] == E_st)
      //return 0;
    }      
  /* 09.14.06  I think that some of the code that checks for cases where we wnat to avoid mapping is unnecessary!,
   * but I'm not sure what...*/
      
  /* check for cases we want to avoid mapping, we exploit the fact that we know a MATP_nd in a sub_cm MUST be mapped 
   * to a corresponding MATP_nd in the original template CM: */

  if(sub_cm->ndtype[sub_nd] == MATP_nd)
    {
      if(orig_cm->ndtype[orig_nd] == MATP_nd && sub_cm->sttype[sub_v] != orig_cm->sttype[orig_v])
	{
	  //printf("return case 1\n");
	  return 0;
	}
      if(orig_cm->ndtype[orig_nd] != MATP_nd)
	{
	  if(!(is_insert))
	    {
	      //printf("return case 2\n");
	      return 0;
	    }
	  else
	    {
	      /* we have a sub MATP insert we're mapping to an orig non-MATP insert, we don't want to do this
	       * but first we want to check to make sure we don't already have a mapping for these two guys,
	       * and if we do, obliterate it.
	       */
  	      //printf("in else\n");
	      if(orig2sub_smap[orig_v][0] != -1 && sub2orig_smap[sub_v][0] != -1)
		{
		  if((orig2sub_smap[orig_v][1] != -1 || (sub2orig_smap[sub_v][1] != -1)))
		    Die("ERROR with sub MATP to orig non-MATP insert 1\n");
		  if(orig2sub_smap[orig_v][0] == sub_v)
		    {
		      printf("obliterate 1\n");
  	      /* obliterate */ 
		      orig2sub_smap[orig_v][0] = -1;
		      sub2orig_smap[sub_v][0]  = -1;
		    }
		}
	      //printf("return case 3\n");
	      return 0;
	    }
	}
    }

  /* check for case that we're mapping a orig MATP_I* to a sub MATP_I*, but we already have a mapping for the
   * sub MATP_I*, in which case we obliterate it. 
   */

  if(is_insert && (orig_cm->ndtype[orig_nd] == MATP_nd && sub_cm->ndtype[sub_nd] == MATP_nd))
    {
      if(orig2sub_smap[orig_v][1] != -1)
	Die("Obliteration error 0\n");
      
      //printf("orig2sub_smap[orig_v:%d][0]: %d\n", orig_v, orig2sub_smap[orig_v][0]);
      if((orig2sub_smap[orig_v][0] != sub_v) && (orig2sub_smap[orig_v][0] != -1))
	{
	  if(sub2orig_smap[orig2sub_smap[orig_v][0]][0] == sub_v)
	    {
	      sub2orig_smap[orig2sub_smap[orig_v][0]][0] = -1;
	      orig2sub_smap[orig_v][0] = -1;
	    }
	  else if(sub2orig_smap[orig2sub_smap[orig_v][0]][1] == sub_v)
	    {
	      sub2orig_smap[orig2sub_smap[orig_v][0]][1] = -1;
	      orig2sub_smap[orig_v][0] = -1;
	    }
	  else
	    {
	      //printf("\n\nobliteration error\n");
	      temp = orig2sub_smap[orig_v][0];
	      //printf("\twe thought either sub2orig_smap[%d][0]:%d == sub_v:%d\n\tor sub2orig_smap[%d][1]:%d == sub_v: %d\n", temp, sub2orig_smap[temp][0], sub_v, temp, sub2orig_smap[temp][1], sub_v);
	    }
	}
    }

  /* check for case that we've already mapped this orig_cm MATP_I* state to a sub_cm MATP_I* state,
   * in which case we want to skip the current mapping.
   */

  /* 091806 09.18.06, I commented out below line b/c comment above implies it had an error, I fixed the
   * error in the line below it, but I'm not sure if this won't break something else.
   */
  //if(is_insert && (orig_cm->ndtype[orig_nd] == MATP_nd && sub_cm->ndtype[sub_nd] != MATP_nd))
  if(is_insert && (orig_cm->ndtype[orig_nd] == MATP_nd && sub_cm->ndtype[sub_nd] == MATP_nd))
    {
      if(orig2sub_smap[orig_v][0] != -1)
	{
	  if((orig2sub_smap[orig_v][1] != -1 || (sub2orig_smap[sub_v][1] != -1)) || (sub2orig_smap[sub_v][0] == -1))
	    {
	      printf("orig2sub_smap[orig_v:%d][0]       (%d)\n", orig_v, orig2sub_smap[orig_v][0]);
	      printf("orig2sub_smap[orig_v:%d][1] != -1 (%d)\n", orig_v, orig2sub_smap[orig_v][1]);
	      printf("sub2orig_smap[sub_v:%d][1]  != -1 (%d)\n", sub_v,  sub2orig_smap[sub_v][1]);
	      printf("sub2orig_smap[sub_v:%d][0]  == -1 (%d)\n", sub_v,  sub2orig_smap[sub_v][0]);
	      Die("ERROR with sub non-MATP to orig MATP insert 2\n");
	    }
	  if(sub_cm->ndtype[sub_cm->ndidx[orig2sub_smap[orig_v][0]]] == MATP_nd)
	    {
	      //printf("return case 4\n");
	      return 0;
	    }
	}
    }


  /* fill in orig2sub_smap */
  if(orig2sub_smap[orig_v][0] == -1)
    {
      if (orig2sub_smap[orig_v][1] != -1) 
	Die("ERROR in map_orig2sub_cm_helper, orig2sub_smap[%d][0] is -1 but orig2sub_smap[%d][1] is not, this shouldn't happen.\n", orig_v, orig_v);
      else
	{
	  orig2sub_smap[orig_v][0] = sub_v;
	  //printf("orig2sub filling 0: orig2sub_smap[%d][0]: %d\n", orig_v, orig2sub_smap[orig_v][0]);
	}
    }
  else if (orig2sub_smap[orig_v][1] != -1)
    {
      if(orig2sub_smap[orig_v][0] == sub_v || orig2sub_smap[orig_v][1] == sub_v)
	/* abort!, we already have this mapping */
	return 0;  
      else
	Die("ERROR in map_orig2sub_cm_helper, orig2sub_smap[%d][0] is not -1 and orig2sub_smap[%d][1] is not -1, this shouldn't happen.\n", orig_v, orig_v);
    }
  else /* orig2sub_smap[orig_v][0] != -1 && orig2sub_smap[orig_v][1] == -1 */
    {
      if(orig2sub_smap[orig_v][0] == sub_v || orig2sub_smap[orig_v][1] == sub_v)
	/* abort!, we already have this mapping */
	return 0; 
      orig2sub_smap[orig_v][1] = sub_v;
      //printf("orig2sub filling 1: orig2sub_smap[%d][1]: %d\n", orig_v, orig2sub_smap[orig_v][1]);
    }

  /* now fill in sub2orig_smap */
  //printf("sub2orig_smap[%d][0]: %d | sub2orig_smap[%d][1]: %d\n", sub_v,  sub2orig_smap[sub_v][0],  sub_v,  sub2orig_smap[sub_v][1]);
  if(sub2orig_smap[sub_v][0] == -1)
    {
      if (sub2orig_smap[sub_v][1] != -1)
	Die("ERROR in map_sub2orig_cm_helper, sub2orig_smap[%d][0] is -1 but sub2orig_smap[%d][1] is not, this shouldn't happen.\n", sub_v, sub_v);
      else
	{
	  sub2orig_smap[sub_v][0] = orig_v;
	  //printf("sub2orig filling 0: sub2orig_smap[%d][0]: %d\n", sub_v, sub2orig_smap[sub_v][0]);
	}
    }
  else if (sub2orig_smap[sub_v][1] != -1)
    Die("ERROR in map_sub2orig_cm_helper, sub2orig_smap[%d][0] is not -1 and sub2orig_smap[%d][1] is not -1, this shouldn't happen.\n", sub_v, sub_v);
  else /* sub2orig_smap[sub_v][0] != -1 && sub2orig_smap[sub_v][1] == -1 */
    {
      sub2orig_smap[sub_v][1] = orig_v;
      //printf("sub2orig filling 1: sub2orig_smap[%d][1]: %d\n", sub_v, sub2orig_smap[sub_v][1]);
    }
  return 1;
}

/**************************************************************************
 * EPN 08.31.06
 * cm2sub_cm_emit_probs()
 *
 * Purpose:  For a specific sub CM state v_s, determine the 
 *           emission probabilities using the emission probs 
 *           of state(s) (up to 2) in the original CM that
 *           map to it. We weight the contribution of each
 *           state to the emission probability by it's psi
 *           value in orig_cm (psi[v] is expected number of 
 *           times state v is entered in a parse).
 *
 * Args:    
 * CM_t *orig_cm     - the original, template CM
 * CM_t *sub_cm      - the sub CM 
 * double *orig_psi  - orig_psi[v] is the expected number of times state v is entered
 *                     in a CM parse
 * int v_s           - the sub_cm state we're filling emissions for
 * int v_o1          - orig_cm state v_s maps to
 * int v_o2          - orig_cm state v_s maps to (-1 if v_s maps to only 1 state)
 * Returns: void
 */
static void
cm2sub_cm_emit_probs(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, int v_s, int v_o1, int v_o2,
		     int *sub2orig_id)
{
  int is_left;
  int i, j;

  //printf("\nin cm2sub_cm_emit_probs v_s: %d, v_o1: %d, v_o2: %d\n", v_s, v_o1, v_o2);

  if(v_o1 == -1)
    Die("ERROR in cm2sub_cm_emit_probs, sub_cm state %d maps to 0 states in orig_cm\n", v_s);

  if(sub_cm->sttype[v_s] == MP_st)
    {
      for(i = 0; i < (MAXABET*MAXABET); i++)
	sub_cm->e[v_s][i] = orig_cm->e[v_o1][i];
      return;
    }

  if(sub2orig_id[v_s] == TRUE)
    {
      /* must be a singlet emitter */
	for(i = 0; i < MAXABET; i++)
	  sub_cm->e[v_s][i] = orig_cm->e[v_o1][i];
      /* No FNorm's necessary (assuming the orig_cm is normalized), since we're
       * building a new CM for each sequence in --sub mode, we skip it for speed.
       */
      return;
    }

  /* If we get here, v_s is a singlet emitter. */
  
  /* There are two cases when two states can map to v_s.
   * Case 1: one of them is an MP_st,
   * Case 2: one is an IL_st and one is an IR_st (ambiguity in CM architecture)
   * These are the only cases where we need to weight emission probs by orig_psi values,
   * and subsequently only cases we need to call FNorm() for */
  if(orig_cm->sttype[v_o1] == MP_st)
    {
      if(orig_cm->sttype[v_o2] == ML_st)
	{
	  is_left = TRUE;
	  //printf("v_o1 MATP, v_s LEFT case\n");
	}
      else if(orig_cm->sttype[v_o2] == MR_st)
	{
	  is_left = FALSE;
	  //printf("v_o1 MATP, v_s RIGHT case\n");
	}
      else
	Die("ERROR v_s: %d maps to a MP_st and another non-ML and non-MR state\n");

      if(v_o2 == -1)
	Die("ERROR sub_cm state: %d maps to a single state orig MATP_MP state: %d\n", v_s, v_o1);
      
      for(i = 0; i < MAXABET; i++)
	if(is_left)
	  for(j = (i*MAXABET); j < ((i+1)*MAXABET); j++)
	    sub_cm->e[v_s][i] += orig_psi[v_o1] * orig_cm->e[v_o1][j];
	else
	  for(j = i; j < (MAXABET*MAXABET); j+=MAXABET)
	    sub_cm->e[v_s][i] += orig_psi[v_o1] * orig_cm->e[v_o1][j];
      if(orig_cm->sttype[v_o2] == MP_st)
	Die("ERROR sub_cm state: %d maps to two MATP_MP states\n", v_s);
      
      /*v_o2 must be ML, MR, IL, or IR, which can all be handled identically */
      //      printf("v_o2 ML, MR, IL or IR\n");
      for(i = 0; i < MAXABET; i++)
	sub_cm->e[v_s][i] += orig_psi[v_o2] * orig_cm->e[v_o2][i];
      FNorm(sub_cm->e[v_s], MAXABET);
      return;
    }
  else if(v_o2 != -1)
    {
      //printf("dual insert case\n");
      if(!((orig_cm->sttype[v_o1] == IL_st && orig_cm->sttype[v_o2] == IR_st) ||
	   (orig_cm->sttype[v_o1] == IR_st && orig_cm->sttype[v_o2] == IL_st)))
	Die("ERROR sub_cm state: %d maps to two states (%d and %d), but not case 1: one is MATP_MP OR case 2: one is IL and one is IR\n", v_s, v_o1, v_o2);
      for(i = 0; i < MAXABET; i++)
	sub_cm->e[v_s][i] = orig_psi[v_o1] * orig_cm->e[v_o1][i];
      for(i = 0; i < MAXABET; i++)
	sub_cm->e[v_s][i] = orig_psi[v_o2] * orig_cm->e[v_o2][i];
      FNorm(sub_cm->e[v_s], MAXABET);
      return;
    }

  if(v_o2 != -1)
    Die("ERROR: dumbass, v_o2 isn't -1\n");

  /* If we get here, v_s maps to a single singlet emitter in orig_cm, v_o1 */
  //printf("singlet emitter case\n");
  for(i = 0; i < MAXABET; i++)
    sub_cm->e[v_s][i] = orig_cm->e[v_o1][i];

  return;
}

/**************************************************************************
 * EPN 08.31.06
 * cm2sub_cm_trans_probs()
 *
 * Purpose:  For a specific sub CM state v_s, determine the 
 *           transition probabilities going out of v_s, 
 *           using the psi values for the original template
 *           CM (orig_cm) (psi[v] is expected number of 
 *           times state v is entered in a parse).
 *           We fill the sub_cm->t[v_s] with 'virtual
 *           counts', then we'll normalize to probabilities
 *           later (outside this function).
 * 
 *           This is based on CP9_cm2wrhmm.c::cm2hmm_trans_probs_cp9()
 *           which was based on formulas/ideas in Zasha Weinberg's
 *           thesis (p.123).
 *
 * Args:    
 * CM_t *orig_cm     - the original, template CM
 * CM_t *sub_cm      - the sub CM 
 * double *orig_psi  - orig_psi[v] is the expected number of times state v is entered
 *                     in a CM parse
 * char ***tmap      - the hard-coded transition map
 * int v_s           - the sub_cm state we're filling transitions for
 * orig2sub_smap     - 2D state map from orig_cm (template) to sub_cm.
 *                     1st dimension - state index in orig_cm 
 *                     2nd D - 2 elements for up to 2 matching sub_cm states, .
 * sub2orig_smap     - 2D state map from sub_cm to orig_cm.
 *                     1st dimension - state index in sub_cm (0..sub_cm->M-1)
 *                     2nd D - 2 elements for up to 2 matching orig_cm states, 
 * Returns: void
 */
static void
cm2sub_cm_trans_probs(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, char ***tmap, int v_s, int **orig2sub_smap,
		      int **sub2orig_smap, int *sub2orig_id)
{
  int v_o;
  int yoffset;
  int y_s;
  int y_o;

  //printf("in cm2sub_cm_trans_probs: v_s: %d\n", v_s);

  if(sub2orig_id[v_s] == TRUE) /* v_s is identical to sub2orig_smap[v_s][0] */
    {
      v_o = sub2orig_smap[v_s][0]; 
      for(yoffset = 0; yoffset < sub_cm->cnum[v_s]; yoffset++)
	{
	  
	  /* we can just copy the transitions */
	  sub_cm->t[v_s][yoffset] = orig_psi[v_o] * orig_cm->t[v_o][yoffset];
	}
      return;
    }

  /* start with the first orig_cm state that maps to v_s */
  v_o = sub2orig_smap[v_s][0];
  //printf("\tv_o: %d\n", v_o);
  if(v_o == -1)
    {
      if(sub_cm->sttype[v_s] != S_st ||
	 sub_cm->sttype[v_s] != E_st ||
	 sub_cm->sttype[v_s] != B_st)
	/* special cases, S_st, E_st, B_st */
	Die("ERROR, sub_cm state v_s: %d maps to no state in sub_cm, but it's not a B, E or S state\n", v_s);
    }
  else
    {
      if(sub_cm->sttype[v_s] == S_st ||
	 sub_cm->sttype[v_s] == E_st ||
	 sub_cm->sttype[v_s] == B_st)
	if(v_s != 0)
	  Die("ERROR, sub_cm state v_s: %d is S, E or B but maps to a orig_cm state: v_o:%d\n", v_o);

      for(yoffset = 0; yoffset < sub_cm->cnum[v_s]; yoffset++)
	{
	  y_s = sub_cm->cfirst[v_s] + yoffset;
	  y_o = sub2orig_smap[y_s][0];
	  //printf("\ty_s: %d | y_o: %d\n", y_s, y_o);
	  
	  if(y_o == -1)
	    {
	      if(sub_cm->sttype[y_s] != B_st &&
		 sub_cm->sttype[y_s] != E_st)
		Die("ERROR, sub_cm state y_s: %d (v_s: %d) maps to no state in sub_cm, but it's not a B or E state\n", y_s, v_s);
	      /* We're transitioning to a BIF_B or END_E; We don't have a map to the orig_cm state
	       * for these guys (its impossible to robustly map these guys), so we handle this
	       * after we've filled in virtual counts for all the other transitions.
	       */
	    }
	  else
	    cm_add_single_trans(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, y_o, v_s, yoffset, orig_psi, tmap);

	  y_o = sub2orig_smap[y_s][1];
	  cm_add_single_trans(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, y_o, v_s, yoffset, orig_psi, tmap);
	}
    }

  /* move on to the second orig_cm state that maps to v_s */
  v_o = sub2orig_smap[v_s][1];
  if(v_o != -1)
    {
      for(yoffset = 0; yoffset < sub_cm->cnum[v_s]; yoffset++)
	{
	  y_s = sub_cm->cfirst[v_s] + yoffset;
	  y_o = sub2orig_smap[y_s][0];
	  if(y_o == -1)
	    {
	      if(sub_cm->sttype[y_s] != B_st &&
		 sub_cm->sttype[y_s] != E_st)
		Die("ERROR, sub_cm state y_s: %d (v_s: %d) maps to no state in sub_cm, but it's not a B or E state\n", y_s, v_s);
	      /* We're transitioning to a BIF_B or END_E; We don't have a map to the orig_cm state
	       * for these guys (its impossible to robustly map these guys), so we handle this
	       * after we've filled in virtual counts for all the other transitions.
	       */
	    }
	  else
	    cm_add_single_trans(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, y_o, v_s, yoffset, orig_psi, tmap);

	  y_o = sub2orig_smap[y_s][1];
	  cm_add_single_trans(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, y_o, v_s, yoffset, orig_psi, tmap);
	}
    }
  return;
}

/**************************************************************************
 * EPN 09.15.06
 * cm2sub_cm_trans_probs_S()
 *
 * Purpose:  For a specific sub CM S state v_start, fill in virtual counts
 *           for transitions out of v_s. We do this in its own seperate function
 *           because we can't robustly map S states in a sub_cm to S states
 *           in an orig CM (if its possible - I can't figure out how to do it).
 * 
 * Args:    
 * CM_t *orig_cm     - the original, template CM
 * CM_t *sub_cm      - the sub CM 
 * double *orig_psi  - orig_psi[v] is the expected number of times state v is entered
 *                     in a CM parse
 * char ***tmap      - the hard-coded transition map
 * int v_start       - the sub_cm start state we're filling virtual counts of transitions into
 * orig2sub_smap     - 2D state map from orig_cm (template) to sub_cm.
 *                     1st dimension - state index in orig_cm 
 *                     2nd D - 2 elements for up to 2 matching sub_cm states, .
 * sub2orig_smap     - 2D state map from sub_cm to orig_cm.
 *                     1st dimension - state index in sub_cm (0..sub_cm->M-1)
 *                     2nd D - 2 elements for up to 2 matching orig_cm states, 
 * Returns: void
 */
static void
cm2sub_cm_trans_probs_S(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, char ***tmap, int v_start, int **orig2sub_smap,
			int **sub2orig_smap)
{
  int yoffset;
  int y_s;

  int sub_nd;

  int v_s_insert;
  int v_o_insert;

  float sum;
  float il_psi;
  float temp_psi;
  float temp_psi_sum;
  float diff;

  int orig_il1, orig_il2, orig_ir1, orig_ir2;

  //printf("in cm2sub_cm_trans_probs_S: v_start: %d\n", v_start);

  sub_nd     = sub_cm->ndidx[v_start];

  if(sub_cm->ndtype[sub_nd] == BEGL_nd)
    {
      /* This is the easy case, we have transitions into each of the states in the split-set
       * of the next node. The only way to reach each of these states is from the BEGL
       * so we just weight each by the psi values for the matching orig_cm states.
       */
      if(sub_cm->ndtype[sub_nd + 1] == BIF_nd)
	sub_cm->t[v_start][0] = 1.0; /* BEGL_S -> BIF_B */
      else
	for(yoffset = 0; yoffset < sub_cm->cnum[v_start]; yoffset++)
	  {
	    y_s = sub_cm->cfirst[v_start] + yoffset;
	    //printf("updating sub_cm->t[%d][%d]\n", v_start, yoffset);
	    sub_cm->t[v_start][yoffset] = orig_psi[sub2orig_smap[y_s][0]];
	    if(sub2orig_smap[y_s][1] != -1)
	      sub_cm->t[v_start][yoffset] += orig_psi[sub2orig_smap[y_s][1]];
	  }
    }
  
  else if(sub_cm->ndtype[sub_nd] == BEGR_nd)
    {
      //printf("!!!SETTING TRANSITIONS OUT OF BEGR_S, NOT SURE IF THIS IS IMPLEMENTED CORRECTLY!!!!\n\n");
      /* More complicated than the BEGL case b/c we need to handle the
       * BEGR_S -> BEGR_IL transition as well as BEGR_S -> next node 
       * split set transitions.
       * We know the BEGR_IL -> BEGR_IL self transition though because
       * we were able to map BEGR_IL to 1 or 2 orig_cm states.
       */

      v_s_insert = v_start + 1;
      if(sub_cm->ndtype[sub_nd + 1] == BIF_nd)
	{
	  //printf("!!!SPECIAL CASE BEGR -> BIF! v_ct\n");

	  v_o_insert = sub2orig_smap[v_s_insert][0];
	  if(sub2orig_smap[v_s_insert][1] == -1)
	    {
	      diff = sub_cm->t[v_s_insert][0] - (orig_psi[v_o_insert] * orig_cm->t[v_o_insert][0]);
	      if(diff >= 0. && diff > 0.0001)
		Die("ERROR, code for calc'ing BEGR_IL -> BIF_B is wrong, sub_cm->t[v_s_insert:%d][0] should be %f (based on your understanding) but its really %f\n", v_s_insert, (orig_psi[v_o_insert] * orig_cm->t[v_o_insert][0]), (sub_cm->t[v_s_insert][0]));
	      if(diff <= 0. && diff < -0.0001)
		Die("ERROR, code for calc'ing BEGR_IL -> BIF_B is wrong, sub_cm->t[v_s_insert:%d][0] should be %f (based on your understanding) but its really %f\n", v_s_insert, (orig_psi[v_o_insert] * orig_cm->t[v_o_insert][0]), (sub_cm->t[v_s_insert][0]));
	    }
	  
	  /* BEGR_IL -> BEGR_IL will be the sum over possibly two orig_cm states v_o_insert that 
	   * map to sub_cm state BEGR_IL of: 
	   *    orig_psi[v_o_insert] * orig_cm->t[v_o_insert][0],
	   * so BEGR_IL -> BIF_B is calc'ed as follows */
	  sub_cm->t[v_s_insert][1] = orig_psi[v_o_insert] * (1. - orig_cm->t[v_o_insert][0]);
	  if(sub2orig_smap[v_s_insert][1] != -1)
	    {
	      /* we have to factor in the other state as well. */
	      Die("ERROR, BEGR_IL: %d maps to 2 orig_cm states, was hoping this was impossible - need to implement.\n", v_s_insert);
	      v_o_insert = sub2orig_smap[v_s_insert][1];
	      sub_cm->t[v_s_insert][1] += orig_psi[v_o_insert] * (1. - orig_cm->t[v_o_insert][0]);
	    }
	  /* Now we normalize the transitions out of BEGR_IL, so we can calculate BEGR_S -> BEGR_IL */
	  FNorm(sub_cm->t[v_s_insert], sub_cm->cnum[v_s_insert]);
	  il_psi   = orig_psi[sub2orig_smap[v_s_insert][0]];
	  if(sub2orig_smap[v_s_insert][1] != -1)
	    {
	      Die("ERROR, BEGR_IL: %d maps to 2 orig_cm states, was hoping this was impossible - need to implement.\n", v_s_insert);
	      il_psi+= orig_psi[sub2orig_smap[v_s_insert][1]];
	    }
	  sub_cm->t[v_start][0] = (1. - sub_cm->t[v_s_insert][0]) * il_psi; /* set BEGR_S -> BEGR_IL */
	  sub_cm->t[v_start][1] = 1. - sub_cm->t[v_start][0]; /* set BEGR_S -> BIF_B */
	}

      else
	{ /* next node is not a BIF node */
	  /* First, normalize the probabilities out of the BEGR_IL, so
	   * we can calculate what BEGR_S -> BEGR_IL should be (we need
	   * to know BEGR_IL -> BEGR_IL probability and orig_psi of
	   * the states that map to BEGR_IL to do this).
	   */
	  sum = 0.;
	  FNorm(sub_cm->t[v_s_insert], sub_cm->cnum[v_s_insert]);
	  il_psi   = orig_psi[sub2orig_smap[v_s_insert][0]];
	  if(sub2orig_smap[v_s_insert][1] != -1)
	    {
	      Die("ERROR, BEGR_IL: %d maps to 2 orig_cm states, was hoping this was impossible - need to implement.\n", v_s_insert);
	      il_psi += orig_psi[sub2orig_smap[v_s_insert][1]];
	    }
	  sub_cm->t[v_start][0] = (1. - sub_cm->t[v_s_insert][0]) * il_psi; /* set BEGR_S -> BEGR_IL */
	  sum = sub_cm->t[v_start][0];
	  
	  for(yoffset = 1; yoffset < sub_cm->cnum[v_start]; yoffset++) /* note we start at yoffset = 1 
									* BEGR_S -> first state of next
									* node's split set. */
	    {
	      y_s = sub_cm->cfirst[v_start] + yoffset;
	      temp_psi   = orig_psi[sub2orig_smap[y_s][0]];
	      if(sub2orig_smap[y_s][1] != -1)
		temp_psi += orig_psi[sub2orig_smap[y_s][1]];
	      
	      sub_cm->t[v_start][yoffset] = temp_psi - il_psi * sub_cm->t[v_s_insert][yoffset];
	      sum += sub_cm->t[v_start][yoffset];
	    }
	  //printf("BEGR->NON BIF  SUM: %f\n", sum);
	  if(sum < 1.0 && ((1.0 - sum) > 0.001))
	    Die("ERROR calculating transitions out of BEGR_S incorrectly\n");
	  if(sum > 1.0 && ((sum - 1.0) > 0.001))
	    Die("ERROR calculating transitions out of BEGR_S incorrectly\n");
	}      
    }
  else if(sub_cm->ndtype[sub_nd] == ROOT_nd)
    {
      //printf("KACHOW! in cm2sub_cm_trans_probs_S(), ROOT_nd\n");
      /* the only case we have to worry about is if the next node is BIF node,
       * otherwise the transitions out of ROOT_S have already been set.
       */
      if(sub_cm->ndtype[sub_nd + 1] == BIF_nd)
	{
	  //printf("!!!SPECIAL CASE ROOT -> BIF!\n");
	  /* Before we do anything we have to check to see if we need to subtract
	   * any subpaths from ROOT_S -> ROOT_IR that have been double counted:
	   */
	  /* if orig_ir < orig_il, we've counted paths from
	   * ir -> il correctly for ROOT_IL -> ROOT_IR and
	   *        incorrectly for ROOT_S  -> ROOT_IR
	   */

	  orig_il1 = sub2orig_smap[1][0];
	  orig_il2 = sub2orig_smap[1][1];
	  
	  orig_ir1 = sub2orig_smap[2][0];
	  orig_ir2 = sub2orig_smap[2][1];
	  
	  if(orig_ir1 < orig_il1)
	    {
	      sub_cm->t[0][1] -= orig_psi[orig_ir1] * 
		cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ir1, orig_il1, 0, tmap, orig_psi, -1, -2, -3, -4);
	    }
	  if(orig_ir2 != -1 && orig_ir2 < orig_il1)
	    {
	      sub_cm->t[0][1] -= orig_psi[orig_ir2] * 
		cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ir2, orig_il1, 0, tmap, orig_psi, -1, -2, -3, -4);
	    }
	  if(orig_il2 != -1 && orig_ir1 < orig_il2)
	    {
	      sub_cm->t[0][1] -= orig_psi[orig_ir1] * 
		cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ir1, orig_il2, 0, tmap, orig_psi, -1, -2, -3, -4);
	    }
	  if((orig_ir2 != -1 && orig_il2 != -1) && (orig_ir2 < orig_il2))
	    {
	      sub_cm->t[0][1] -= orig_psi[orig_ir2] * 
		cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ir2, orig_il2, 0, tmap, orig_psi, -1, -2, -3, -4);
	    }
	  
	  /* Next, set transition from ROOT_S -> BIF_B, we know that
	   * ROOT_S -> ROOT_IL and ROOT_S -> ROOT_IR were set using 
	   * orig_cm subpaths originating at ROOT_S, so we know that
	   * the virtual counts out of ROOT_S should NOT be scaled, in
	   * other words, they can be treated as probabilities and should
	   * sum to 1.0.
	   */
	  sub_cm->t[0][2] = 1.0 - (sub_cm->t[0][0] + sub_cm->t[0][1]); /* set ROOT_S->BIF_B */

	  v_s_insert = v_start + 1; /* ROOT_IL */
	  v_o_insert = sub2orig_smap[v_s_insert][0];
	  if(sub2orig_smap[v_s_insert][1] == -1)
	    {
	      diff = sub_cm->t[v_s_insert][0] - (orig_psi[v_o_insert] * orig_cm->t[v_o_insert][0]);
	      if(diff >= 0. && diff >= 0.0001)
		Die("ERROR, code for calc'ing ROOT_IL -> BIF_B is wrong, sub_cm->t[v_s_insert:%d][0] should be %f (based on your understanding) but its really %f\n", v_s_insert, (orig_psi[v_o_insert] * orig_cm->t[v_o_insert][0]), (sub_cm->t[v_s_insert][0]));
	      if(diff <= 0. && diff <= -0.0001)
		Die("ERROR, code for calc'ing ROOT_IL -> BIF_B is wrong, sub_cm->t[v_s_insert:%d][0] should be %f (based on your understanding) but its really %f\n", v_s_insert, (orig_psi[v_o_insert] * orig_cm->t[v_o_insert][0]), (sub_cm->t[v_s_insert][0]));
	    }
	  /* We want to figure out what the virtual counts for the traansition ROOT_IL -> BIF_B should
	   * be, but this depends on how we filled the virtual counts for the other two 
	   * transitions (-> ROOT_IL and -> ROOT_IR) out of ROOT_IL, so we SHOULD revisit
	   * their calculation. Here, I'm attempting a trick that SHOULD work: figuring out
	   * what the self_insert probability IL->IL should be based on the self insert probabilities
	   * of the up to 2 orig_cm states that sub_cm ROOT_IL maps to, then using that scaling
	   * factor (the actual probability / the virtual counts currently in sub_cm->t[ROOT_IL][0])
	   * to scale both ROOT_IL->ROOT_IL and ROOT_IL->ROOT_IR, then we can just set ROOT_IL->BIF
	   * as 1.0 - (ROOT_IL->ROOT_IL + ROOT_IL->ROOT_IR).
	   */

	  temp_psi_sum = orig_psi[v_o_insert];
	  
	  if(sub2orig_smap[v_s_insert][1] != -1)
	    {
	      v_o_insert    = sub2orig_smap[v_s_insert][1];
	      temp_psi_sum += orig_psi[v_o_insert];
	    }	      
	  sub_cm->t[v_s_insert][0] /= temp_psi_sum; /* ROOT_IL -> ROOT_IL */
	  sub_cm->t[v_s_insert][1] /= temp_psi_sum; /* ROOT_IL -> ROOT_IR */
	  sub_cm->t[v_s_insert][2]  = 1. - (sub_cm->t[v_s_insert][0] + sub_cm->t[v_s_insert][1]);
 	                                            /* ROOT_IL -> BIF_B */

	  /* move on to calc'ing ROOT_IR -> BIF */
	  v_s_insert = v_start + 2; /* ROOT_IR */
	  v_o_insert = sub2orig_smap[v_s_insert][0];
	  if(sub2orig_smap[v_s_insert][1] == -1)
	    {
	      diff = sub_cm->t[v_s_insert][0] - (orig_psi[v_o_insert] * orig_cm->t[v_o_insert][0]);
	      if(diff >= 0. && diff > 0.0001)
		Die("ERROR, code for calc'ing ROOT_IR -> BIF_B is wrong, sub_cm->t[v_s_insert:%d][0] should be %f (based on your understanding) but its really %f\n", v_s_insert, (orig_psi[v_o_insert] * orig_cm->t[v_o_insert][0]), (sub_cm->t[v_s_insert][0]));
	      if(diff <= 0. && diff < -0.0001)
		Die("ERROR, code for calc'ing ROOT_IR -> BIF_B is wrong, sub_cm->t[v_s_insert:%d][0] should be %f (based on your understanding) but its really %f\n", v_s_insert, (orig_psi[v_o_insert] * orig_cm->t[v_o_insert][0]), (sub_cm->t[v_s_insert][0]));
	    }
	  /* We can calculate what the ROOT_IR -> ROOT_IR probability should be, and then
	   * just take 1.0 minus that probability to set ROOT_IR -> BIF.
	   */
	  temp_psi_sum = orig_psi[v_o_insert];
	  
	  if(sub2orig_smap[v_s_insert][1] != -1)
	    {
	      v_o_insert    = sub2orig_smap[v_s_insert][1];
	      temp_psi_sum += orig_psi[v_o_insert];
	    }	      
	  sub_cm->t[v_s_insert][0] /= temp_psi_sum; /* ROOT_IR -> ROOT_IR */
	  sub_cm->t[v_s_insert][1]  = 1. - (sub_cm->t[v_s_insert][0]);
 	                                            /* ROOT_IR -> BIF_B */
	}
      else
	{
	  /* ROOT -> non-BIF node, we have already handled this, so we return. */
	}
    }
  //printf("leaving cm2sub_cm_trans_probs_S\n\n");
  return;
}

/**************************************************************************
 * EPN 09.13.06
 * cm2sub_cm_trans_probs_B_E()
 *
 * Purpose:  For a specific sub CM B or E state v_be fill in virtual counts
 *           for transitions into v_be. We do this in its own seperate function,
 *           because we can't robustly map E states in a sub_cm to E states
 *           in an orig CM (if its possible - I can't figure out how to do it).
 * 
 * Args:    
 * CM_t *orig_cm     - the original, template CM
 * CM_t *sub_cm      - the sub CM 
 * double *orig_psi  - orig_psi[v] is the expected number of times state v is entered
 *                     in a CM parse
 * char ***tmap      - the hard-coded transition map
 * int v_be         - the sub_cm END state we're filling virtual counts of transitions into
 * orig2sub_smap     - 2D state map from orig_cm (template) to sub_cm.
 *                     1st dimension - state index in orig_cm 
 *                     2nd D - 2 elements for up to 2 matching sub_cm states, .
 * sub2orig_smap     - 2D state map from sub_cm to orig_cm.
 *                     1st dimension - state index in sub_cm (0..sub_cm->M-1)
 *                     2nd D - 2 elements for up to 2 matching orig_cm states, 
 * Returns: void
 */
static void
cm2sub_cm_trans_probs_B_E(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, char ***tmap, int v_be, int **orig2sub_smap,
			int **sub2orig_smap)
{
  int v_o;
  int v_s;

  int sub_nd;
  int psub_nd;
  int v_s_insert;
  int bif_end_yoffset;

  //printf("in cm2sub_cm_trans_probs_B_E: v_s: %d\n", v_be);

  sub_nd =  sub_cm->ndidx[v_be];
  psub_nd = sub_nd - 1;

  switch (sub_cm->ndtype[psub_nd]) {
  case MATP_nd:
    //printf("prev node type MATP\n");
    /* psub_nd is a MATP: each state in psub_nd transits to 
     * exactly 3 states, the MATP_IL, MATP_IR and the BIF or END state v_be, for which
     * we are trying to fill virtual transition counts into.
     *
     * The case for MATP_nd's is actually simpler than for other nodes because we
     * can exploit the fact that if a MATP node exists in the sub_cm it necessarily
     * must map to a MATP node in the original CM (we're not adding or changing any 
     * base pairs, just deleting some possibly), so every sub MATP state must correspond
     * to exactly 1 original MATP state.
       */
    bif_end_yoffset = 2; /* cm->t[][2] goes to BIF_B or END_E */
    v_s = sub_cm->nodemap[psub_nd]; /* v_s is MATP_MP */
    v_o = sub2orig_smap[v_s][0];
    if(sub2orig_smap[v_s][1] == -1)
      {
	/* virtual counts for the MATP_MP -> BIF or END is just the remainder of the 
	 * prob mass (orig_psi[MATP_MP]) that doesn't go into MATP_IL or MATP_IR.
	 * that doesn't go into MAT*_I* */
	sub_cm->t[v_s][2] = orig_psi[v_o] * (1. - ((sub_cm->t[v_s][0] + sub_cm->t[v_s][1]) / 
						   orig_psi[v_o]));
      }
    else 
      Die("ERROR MAT_MP state: %d maps to 2 states in the orig_cm\n", v_s);
    
    /* move on to transitions out of next state MATP_ML, we do this guy analogously */
    v_s++; /* v_s is MATP_ML */
    v_o = sub2orig_smap[v_s][0];
    if(sub2orig_smap[v_s][1] == -1)
	  sub_cm->t[v_s][2] = orig_psi[v_o] * (1. - ((sub_cm->t[v_s][0] + sub_cm->t[v_s][1]) / 
						     orig_psi[v_o]));
    else
      Die("ERROR MATP_ML state: %d maps to 2 states in the orig_cm\n", v_s);
    
    /* MATP_MR */
    v_s++; /* v_s is MATP_MR */
    v_o = sub2orig_smap[v_s][0];
    if(sub2orig_smap[v_s][1] == -1)
      sub_cm->t[v_s][2] = orig_psi[v_o] * (1. - ((sub_cm->t[v_s][0] + sub_cm->t[v_s][1]) / 
						 orig_psi[v_o]));
    else
      Die("ERROR MATP_MR state: %d maps to 2 states in the orig_cm\n", v_s);
    
    /* MATP_D */
    v_s++; /* v_s is MATP_D */
    v_o = sub2orig_smap[v_s][0];
    if(sub2orig_smap[v_s][1] == -1)
      sub_cm->t[v_s][2] = orig_psi[v_o] * (1. - ((sub_cm->t[v_s][0] + sub_cm->t[v_s][1]) / 
						 orig_psi[v_o]));
    else
      Die("ERROR MATP_D state: %d maps to 2 states in the orig_cm\n", v_s);
    
    /* MATP_IR */
    v_s++; /* v_s is MATP_IR */
    v_o = sub2orig_smap[v_s][0];
    if(sub2orig_smap[v_s][1] == -1)
      sub_cm->t[v_s][2] = orig_psi[v_o] * (1. - ((sub_cm->t[v_s][0] + sub_cm->t[v_s][1]) / 
						 orig_psi[v_o]));
    else
      Die("ERROR MATP_IR state: %d maps to 2 states in the orig_cm\n", v_s);
    
    /* MATP_IR */
    v_s++; /* v_s is MATP_IR */
    v_o = sub2orig_smap[v_s][0];
    if(sub2orig_smap[v_s][1] == -1)
      sub_cm->t[v_s][1] = orig_psi[v_o] * (1. - (sub_cm->t[v_s][0] / 
						 orig_psi[v_o]));
    else
      Die("ERROR MATP_IL state: %d maps to 2 states in the orig_cm\n", v_s);
    break;

  case MATL_nd:
  case MATR_nd:
    //printf("prev node type MATL or MATR\n");
    /* psub_nd is a MATL or MATR and we know each state in psub_nd transits to 
     * exactly 2 states, the insert state of psub_nd and the BIF or END state v_be, which
     * we are trying to fill in transitions to.
     */
    bif_end_yoffset = 1; /* cm->t[][1] goes to BIF_B or END_E */
    v_s = sub_cm->nodemap[psub_nd]; /* v_s is MATL_ML or MATR_MR */
    v_s_insert = sub_cm->nodemap[psub_nd] + 2; /* v_s_insert is MATL_IL or MATR_IR */
    v_o = sub2orig_smap[v_s][0];
    if(sub2orig_smap[v_s][1] == -1)
      {
	/* virtual counts for the MAT*_* -> BIF or END is just the remainder of the orig_psi[MAT*_*] 
	 * that doesn't go into MAT*_I* */
	sub_cm->t[v_s][1] = orig_psi[v_o] * (1. - (sub_cm->t[v_s][0] / orig_psi[v_o]));
      }
    else /* v_s maps to 2 states in orig_cm, this makes things more complicated. */
      {
	/* all the info we have is orig_psi[v_o] and the transition virtual counts into MATL_IL, but
	 * if two orig_cm states map to v_s, then we need more info: specifically the contribution
	 * of each of these two orig_cm states to the transition virtual counts, so we need to wastefully
	 * recall cm_sum_subpaths().
	 */
	cm2sub_cm_trans_probs_B_E_dual_helper(orig_cm, sub_cm, orig_psi, tmap, v_s, v_s_insert, bif_end_yoffset,
					      orig2sub_smap, sub2orig_smap);
      }
    /* move on to transitions out of next state MAT*_D, we do this guy analogously */
    v_s++; /* v_s is MATL_D or MATR_D */
    v_o = sub2orig_smap[v_s][0];
    if(sub2orig_smap[v_s][1] == -1)
      sub_cm->t[v_s][1] = orig_psi[v_o] * (1. - (sub_cm->t[v_s][0] / orig_psi[v_o]));
    else
      cm2sub_cm_trans_probs_B_E_dual_helper(orig_cm, sub_cm, orig_psi, tmap, v_s, v_s_insert, bif_end_yoffset,
					    orig2sub_smap, sub2orig_smap);
    
    /* move on to transitions out of next state MAT*_I*, we do this guy analogously */
    v_s++; /* v_s is MATL_IL or MATR_IR */
    v_o = sub2orig_smap[v_s][0];
    if(sub2orig_smap[v_s][1] == -1)
      sub_cm->t[v_s][1] = orig_psi[v_o] * (1. - (sub_cm->t[v_s][0] / orig_psi[v_o]));
    else
      cm2sub_cm_trans_probs_B_E_dual_helper(orig_cm, sub_cm, orig_psi, tmap, v_s, v_s_insert, bif_end_yoffset,
					    orig2sub_smap, sub2orig_smap);
    break;
    
  case ROOT_nd:
  case BEGL_nd:
  case BEGR_nd:
    { /* we handle these cases in cm2sub_cm_trans_probs_S() */ }
    break;

  default: 
    Die("ERROR bogus node type transiting to END or BIF\n");
    break;
  }

  //printf("Returning from cm2sub_cm_trans_probs_B_E\n");
  return;
}


/**************************************************************************
 * EPN 09.13.06
 * Function: cm2sub_cm_trans_probs_B_E_dual_helper()
 *
 * Purpose:  Add virtual counts contribution into a B or E state for a sub_cm state
 *           that maps to 2 orig_cm states (a special case).
 *
 * See related functions for explanation of parameters. 
 * Returns: (void) 
 */
static void
cm2sub_cm_trans_probs_B_E_dual_helper(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, char ***tmap, int v_s, 
				      int v_s_insert, int bif_end_yoffset, int **orig2sub_smap, int **sub2orig_smap)
{
  int v_o;
  int v_o_insert;
  float contribution;
  int start;
  int end;
  int going_down;
  int is_insert; 
  
  is_insert = FALSE;
  if(sub_cm->sttype[v_s] == IL_st || sub_cm->sttype[v_s] == IR_st)
    is_insert = TRUE;

  /* deal with the first orig_cm state that maps to v_s */
  v_o = sub2orig_smap[v_s][0];
  v_o_insert = sub2orig_smap[v_s_insert][0]; 

  printf("in cm2sub_cm_trans_prob_B_E_dual_helper(), v_s: %d v_s_insert: %d\n\n", v_s, v_s_insert);
  printf("\ttrans to insert: sub_cm->t[v_s:%d][0]: %f\n", v_s, sub_cm->t[v_s][0]);
  printf("\tbeg sub_cm->t[v_s:%d][%d]: %f\n", v_s, bif_end_yoffset, sub_cm->t[v_s][bif_end_yoffset]);
  printf("first v_o: %d v_o_insert: %d\n", v_o, v_o_insert);
  if(v_o <= v_o_insert) /* going DOWN the CM */
    {
      start = v_o;
      end   = v_o_insert;
      going_down = TRUE;
    }
  else /* going UP the CM */
    {
      end = v_o;
      start = v_o_insert;
      going_down = FALSE;
    }
  contribution = cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, start, end, v_s, tmap, orig_psi, -1, -2, -3, -4); 

  v_o_insert = sub2orig_smap[v_s_insert][1]; 
  if(v_o_insert != -1)
    {
      if(!is_insert && ((!(going_down) && (v_o <= v_o_insert)) ||
			(going_down && (v_o > v_o_insert))))
	 /* shouldn't happen, a non-insert maps to 
	  * to 2 orig_cm states, one of which is before 
	  * v_o_insert and one of which is after it */

	Die("ERROR, cm2sub_cm_trans_probs_B_E_dual_helper() you don't understand CM architecture 0, v_s: %d, v_s_insert: %d, v_o: %d, v_o_insert: %d\n", v_s, v_s_insert, v_o, v_o_insert);
      if(going_down)
	{
	  start = v_o;
	  end   = v_o_insert;
	}
      else
	{
	  end   = v_o;
	  start = v_o_insert;
	}
      contribution += cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, start, end, v_s, tmap, orig_psi, -1, -2, -3, -4); 
    }
  printf("\t\tcontribution: %f | orig_psi[%d]: %f\n", contribution, start, orig_psi[start]);
  sub_cm->t[v_s][bif_end_yoffset] = orig_psi[start] * (1. - contribution);
  printf("\tmid sub_cm->t[v_s:%d][%d]: %f\n", v_s, bif_end_yoffset, sub_cm->t[v_s][bif_end_yoffset]);

  /* now for the second orig_cm state that maps to v_s */
  v_o = sub2orig_smap[v_s][1];
  if(v_o == -1)
    Die("ERROR in cm2sub_cm_trans_probs_B_E_dual_helper\n");

  v_o_insert = sub2orig_smap[v_s_insert][0]; 
  printf("second v_o: %d v_o_insert: %d\n", v_o, v_o_insert);
  if(v_o <= v_o_insert) /* going DOWN the CM */
    {
      start = v_o;
      end = v_o_insert;
      going_down = TRUE;
    }
  else /* going UP the CM */
    {
      end = v_o;
      start = v_o_insert;
      going_down = FALSE;
    }
  contribution = cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, start, end, v_s, tmap, orig_psi, -1, -2, -3, -4); 
  v_o_insert = sub2orig_smap[v_s_insert][1]; 
  if(v_o_insert != -1)
    {
      if((!(going_down) && (v_o <= v_o_insert)) ||
	  (going_down && (v_o > v_o_insert))) /* shouldn't happen, a non-insert maps to 
					       * to 2 orig_cm states, one of which is before 
					       * v_o_insert and one of which is after it */
	 Die("ERROR, cm2sub_cm_trans_probs_B_E_dual_helper() you don't understand CM architecture 0\n");
      if(going_down)
	{
	  start = v_o;
	  end   = v_o_insert;
	}
      else
	{
	  end   = v_o;
	  start = v_o_insert;
	}
      contribution += cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, start, end, v_s, tmap, orig_psi, -1, -2, -3, -4); 
    }

  sub_cm->t[v_s][bif_end_yoffset] += orig_psi[start] * (1. - contribution);
  
  printf("\t\tcontribution: %f | orig_psi[%d]: %f\n", contribution, start, orig_psi[start]);
  printf("\tend sub_cm->t[v_s:%d][%d]: %f\n", v_s, bif_end_yoffset, sub_cm->t[v_s][bif_end_yoffset]);
  return;
}

/**************************************************************************
 * EPN 09.21.06
 * cm2sub_cm_trans_probs_B_E2()
 *
 * Purpose:  For a specific sub CM B or E state v_be fill in virtual counts
 *           for transitions into v_be. We do this in its own seperate function,
 *           because we can't robustly map E states in a sub_cm to E states
 *           in an orig CM (if its possible - I can't figure out how to do it).
 * 
 * Args:    
 * CM_t *orig_cm     - the original, template CM
 * CM_t *sub_cm      - the sub CM 
 * double *orig_psi  - orig_psi[v] is the expected number of times state v is entered
 *                     in a CM parse
 * char ***tmap      - the hard-coded transition map
 * int v_be         - the sub_cm END state we're filling virtual counts of transitions into
 * orig2sub_smap     - 2D state map from orig_cm (template) to sub_cm.
 *                     1st dimension - state index in orig_cm 
 *                     2nd D - 2 elements for up to 2 matching sub_cm states, .
 * sub2orig_smap     - 2D state map from sub_cm to orig_cm.
 *                     1st dimension - state index in sub_cm (0..sub_cm->M-1)
 *                     2nd D - 2 elements for up to 2 matching orig_cm states, 
 * print_flag        - TRUE to print useful debugging info
 * Returns: void
 */
static void
cm2sub_cm_trans_probs_B_E2(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, char ***tmap, int v_be, int **orig2sub_smap,
			   int **sub2orig_smap, int print_flag)
{
  int orig_v;
  int sub_v;

  int sub_nd;
  int psub_nd;
  int sub_i;
  int bif_end_yoffset;
  int orig_v1, orig_v2, orig_i1, orig_i2;
  float contribution;
  int orig_il, orig_ir, sub_il, sub_ir;

  if(print_flag) printf("in cm2sub_cm_trans_probs_B_E: sub_v: %d\n", v_be);

  sub_nd =  sub_cm->ndidx[v_be];
  psub_nd = sub_nd - 1;

  switch (sub_cm->ndtype[psub_nd]) {
  case MATP_nd:
    if(print_flag) printf("prev node type MATP\n");
    /* psub_nd is a MATP: each state in psub_nd transits to 
     * exactly 3 states, the MATP_IL, MATP_IR and the BIF or END state v_be, for which
     * we are trying to fill virtual transition counts into.
     *
     * The case for MATP_nd's is actually simpler than for other nodes because we
     * can exploit the fact that if a MATP node exists in the sub_cm it necessarily
     * must map to a MATP node in the original CM (we're not adding or changing any 
     * base pairs, just deleting some possibly), so every sub MATP state must correspond
     * to exactly 1 original MATP state.
       */

    bif_end_yoffset = 2; /* cm->t[][2] goes to BIF_B or END_E */
    sub_il = sub_cm->nodemap[psub_nd] + 4; 
    sub_ir = sub_cm->nodemap[psub_nd] + 5; 

    orig_il = sub2orig_smap[sub_il][0];
    orig_ir = sub2orig_smap[sub_ir][0];

    for(sub_v = sub_cm->nodemap[psub_nd]; sub_v < sub_il; sub_v++)
      {
	orig_v = sub2orig_smap[sub_v][0];
	sub_cm->t[sub_v][bif_end_yoffset] = orig_psi[orig_v] - 
	  (sub_cm->t[sub_v][0] + sub_cm->t[sub_v][1]);
      }
    /* now do MATP_IL and MATP_IR */
    sub_cm->t[sub_il][bif_end_yoffset] = orig_psi[orig_il] - 
      (sub_cm->t[sub_il][0] + sub_cm->t[sub_il][1]);

    bif_end_yoffset = 1;
    sub_cm->t[sub_ir][bif_end_yoffset] = orig_psi[orig_ir] - 
      sub_cm->t[sub_ir][0];
    
    break;

  case MATL_nd:
  case MATR_nd:
    if(print_flag) printf("prev node type MATL or MATR\n");
    /* psub_nd is a MATL or MATR and we know each state in psub_nd transits to 
     * exactly 2 states, the insert state of psub_nd and the BIF or END state v_be, which
     * we are trying to fill in transitions to.
     */
    bif_end_yoffset = 1; /* cm->t[][1] goes to BIF_B or END_E */
    
    sub_i = sub_cm->nodemap[psub_nd] + 2; /* sub_i is MATL_IL or MATR_IR */
    orig_i1 = sub2orig_smap[sub_i][0];
    orig_i2 = sub2orig_smap[sub_i][1];

    for(sub_v = sub_cm->nodemap[psub_nd]; sub_v < sub_i; sub_v++)
      {
	orig_v1 = sub2orig_smap[sub_v][0];
	orig_v2 = sub2orig_smap[sub_v][1];
	
	if(orig_v1 < orig_i1)
	  {
	    contribution = cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, 
					   orig_v1, orig_i1, sub_v, tmap, orig_psi, -1, -2, -3, -4); 
	    sub_cm->t[sub_v][bif_end_yoffset] = orig_psi[orig_v1] * (1. - contribution);
	    if(print_flag) printf("curr v1 < i1 sub_cm->t[sub_v:%d][1] now: %f (added: psi:%f * 1-cont: %f (%f))\n", sub_v, sub_cm->t[sub_v][1], orig_psi[orig_v1], (1.-contribution), (orig_psi[orig_v1] * (1. - contribution)));
	  }
	else
	  {
	    contribution = orig_psi[orig_i1] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, 
							       orig_i1, orig_v1, sub_v, tmap, orig_psi, -1, -2, -3, -4); 
	    sub_cm->t[sub_v][bif_end_yoffset] = orig_psi[orig_v1] * (1. - (contribution / orig_psi[orig_v1]));
	    if(print_flag) printf("curr i1 < v1 sub_cm->t[sub_v:%d][1] now: %f (added: psi:%f * 1-cont: %f (%f))\n", sub_v, sub_cm->t[sub_v][1], orig_psi[orig_v1], (1.-contribution), (orig_psi[orig_v1] * (1. - contribution)));
	  }
	if(orig_v2 != -1)
	  {
	    if(orig_v2 < orig_i1)
	      {
		contribution = cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, 
					       orig_v2, orig_i1, sub_v, tmap, orig_psi, -1, -2, -3, -4); 
		sub_cm->t[sub_v][bif_end_yoffset] += orig_psi[orig_v2] * (1. - contribution);
		if(print_flag) printf("curr v2 < i1 sub_cm->t[sub_v:%d][1] now: %f (added: psi:%f * 1-cont: %f (%f))\n", sub_v, sub_cm->t[sub_v][1], orig_psi[orig_v2], (1.-contribution), (orig_psi[orig_v2] * (1. - contribution)));
	      }
	    else
	      {
		contribution = orig_psi[orig_i1] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, 
								   orig_i1, orig_v2, sub_v, tmap, orig_psi, -1, -2, -3, -4); 
		sub_cm->t[sub_v][bif_end_yoffset] += orig_psi[orig_v2] * (1. - (contribution / orig_psi[orig_v2]));
		if(print_flag) printf("curr i1 < v2 sub_cm->t[sub_v:%d][1] now: %f (added: psi:%f * 1-cont: %f (%f))\n", sub_v, sub_cm->t[sub_v][1], orig_psi[orig_v2], (1.-contribution), (orig_psi[orig_v2] * (1. - contribution)));
	      }
	  }
	
	if(orig_i2 != -1)
	  {
	    if(orig_v1 < orig_i2)
	      {
		contribution = cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, 
					       orig_v1, orig_i2, sub_v, tmap, orig_psi, -1, -2, -3, -4); 
		sub_cm->t[sub_v][bif_end_yoffset] += orig_psi[orig_v1] * (1. - contribution);
	      }
	    else
	      {
		contribution = orig_psi[orig_i2] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, 
								   orig_i2, orig_v1, sub_v, tmap, orig_psi, -1, -2, -3, -4); 
		sub_cm->t[sub_v][bif_end_yoffset] += orig_psi[orig_v1] * (1. - (contribution / orig_psi[orig_v1]));
	      }
	    if(orig_v2 != -1)
	      {
		if(orig_v2 < orig_i2)
		  {
		    contribution = cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, 
						   orig_v2, orig_i2, sub_v, tmap, orig_psi, -1, -2, -3, -4); 
		    sub_cm->t[sub_v][bif_end_yoffset] += orig_psi[orig_v2] * (1. - contribution);
		  }
		else
		  {
		    contribution = orig_psi[orig_i2] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, 
								       orig_i2, orig_v2, sub_v, tmap, orig_psi, -1, -2, -3, -4); 
		    sub_cm->t[sub_v][bif_end_yoffset] += orig_psi[orig_v2] * (1. - (contribution / orig_psi[orig_v2]));
		  }
	      }
	  }
      }
    /* now set the sub_i -> B or E transition prob */
    contribution = cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, 
				   orig_i1, orig_i1, sub_i, tmap, orig_psi, -1, -2, -3, -4); 
    sub_cm->t[sub_i][1] = orig_psi[orig_i1] * (1. - contribution);
    if(orig_i2 != -1)
      {
	contribution = cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, 
				       orig_i2, orig_i2, sub_i, tmap, orig_psi, -1, -2, -3, -4); 
	sub_cm->t[sub_i][1] += orig_psi[orig_i2] * (1. - contribution);
      }
    break;
	
  case ROOT_nd:
  case BEGL_nd:
  case BEGR_nd:
    { /* we handle these cases in cm2sub_cm_trans_probs_S() */ }
    break;

  default: 
    Die("ERROR bogus node type transiting to END or BIF\n");
    break;
  }

  if(print_flag) printf("Returning from cm2sub_cm_trans_probs_B_E\n");
  return;
}

/**************************************************************************
 * EPN 08.31.06
 * Function: cm_add_single_trans()
 *
 * Purpose:  Add a virtual counts contribution to a single CM transition,
 *           and subtract subpaths if necessary.
 * 
 * See related functions for explanation of parameters. 
 * Returns: (void) 
 */
static void
cm_add_single_trans(CM_t *orig_cm, CM_t *sub_cm, int **orig2sub_smap, int **sub2orig_smap, int orig_v, int orig_y, 
		    int sub_v, int yoffset, double *orig_psi, char ***tmap)
{
  int y_s;

  y_s = sub_cm->cfirst[sub_v] + yoffset;
  //printf("\t\tin cm_add_single_trans, orig_v: %d | orig_y: %d | sub_v: %d | y_s: %d\n", orig_v, orig_y, sub_v, y_s);
  //printf("\t\tbeg sub_cm->t[sub_v:%d][yoffset:%d] : %.9f\n", sub_v, yoffset, sub_cm->t[sub_v][yoffset]);

  /* check if we've got real CM state ids */
  if(orig_v == -1 || orig_y == -1)
    return;

  int start;
  int end;
  int orig_insert1, orig_insert2, orig_insert3, orig_insert4;
  int sub_insert1, sub_insert2;
  int n_s;
  int sub_i;
  int orig_i1;
  int orig_i2;

  if(orig_v <= orig_y) /* going DOWN the CM */
    {
      start = orig_v;
      end = orig_y;
    }
  else /* going UP the CM */
    {
      start = orig_y;
      end = orig_v;
    }
  
  /* find the insert states in orig_cm that map to the sub_cm insert states in the same node as sub_v,
   * cm_sum_subpaths needs to know these states */
  orig_insert1 = orig_insert2 = orig_insert3 = orig_insert4 = -1;
  sub_insert1  = sub_insert2 = -1;
  n_s          = sub_cm->ndidx[sub_v];
  /* find the inserts in n_s, we know each state in a node has all the inserts in that node
   * as children (if there are any), and further these must be the top 2 lowest numbered children.
   * We're exploiting our knowledge of the CM architecture here -- not good. */
  if(sub_cm->sttype[sub_v] != B_st)
    if(sub_cm->sttype[sub_cm->cfirst[sub_v]] == IL_st || 
       sub_cm->sttype[sub_cm->cfirst[sub_v]] == IR_st)
      {
	sub_insert1 = sub_cm->cfirst[sub_v];
	if(sub_cm->sttype[sub_cm->cfirst[sub_v]+1] == IR_st)
	  sub_insert2 = sub_cm->cfirst[sub_v] + 1;
      }
  if(sub_insert1 != -1)
    {
      orig_insert1 = sub2orig_smap[sub_insert1][0];
      orig_insert2 = sub2orig_smap[sub_insert1][1];
      if(sub_insert2 != -1)
	{
	  orig_insert3 = sub2orig_smap[sub_insert2][0];
	  orig_insert4 = sub2orig_smap[sub_insert2][1];
	}
    }

  //  printf("\t\t\toi1: %d oi2: %d oi3: %d oi4: %d\n", orig_insert1, orig_insert2, orig_insert3, orig_insert4);
  
  sub_cm->t[sub_v][yoffset] += orig_psi[start] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, start, end, sub_v, tmap, orig_psi, orig_insert1, orig_insert2, orig_insert3, orig_insert4);

  //printf("\t\tmid sub_cm->t[sub_v:%d][yoffset:%d] : %.9f\n", sub_v, yoffset, sub_cm->t[sub_v][yoffset]);

  /* 09.21.06 situation check: see if we have a situation where we've counted
   * some subpaths that we didn't want to - i.e. that will be counted again
   * out of other states. This only applies if we're in a MATL or MATR node:
   * Six potential cases, 1 of which is impossible due to CM construction rules:
   * 
   * Let orig_v is the state we're transitioning out of 
   *     orig_y is the state we're going into, 
   *     orig_i is the orig_cm state that maps to the insert state in the same node as orig_v, 
   *
   * Cases 1A, 1B, 1C all have orig_v < orig_y
   * Case  1A: orig_v < orig_i < orig_y
   * Case  1B: orig_i < orig_v < orig_y
   * Case  1C: orig_v < orig_y < orig_i IMPOSSIBLE
   *
   * Cases 2A, 2B, 2C all have orig_y < orig_v
   * Case  2A: orig_y < orig_i  < orig_v
   * Case  2B: orig_i < orig_y < orig_v
   * Case  2C: orig_y < orig_v  < orig_i
   */
  if(sub_cm->ndtype[sub_cm->ndidx[sub_v]] == MATL_nd || 
     sub_cm->ndtype[sub_cm->ndidx[sub_v]] == MATR_nd)
    {
      /* First determine insert state of this node */
      sub_i = sub_cm->nodemap[sub_cm->ndidx[sub_v]] + 2;
      if(sub_i != y_s) /* don't subtract paths when we're going to sub_i */
	{
	  //printf("sub_v: %d sub_i: %d\n", sub_v, sub_i);
	  orig_i1 = sub2orig_smap[sub_i][0];
	  orig_i2 = sub2orig_smap[sub_i][1];
	  
	  /* call our helper function */
	  cm_subtract_subpaths(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, orig_v, orig_y,
			       orig_i1, sub_v, yoffset, orig_psi, tmap);
	  
	  if(orig_i2 != -1)
	    cm_subtract_subpaths(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, orig_v, orig_y,
				 orig_i2, sub_v, yoffset, orig_psi, tmap);
	}
    }
}	  

/**************************************************************************
 * EPN 09.21.06
 * Function: cm_subtract_subpaths()
 *
 * Purpose:  Subtract subpaths from a given state. There's
 *           6 possible cases.
 *
 *           Let orig_v is the state we're transitioning out of 
 *           orig_y is the state we're going into, 
 *           orig_i is the orig_cm state that maps to the insert state in the same node as orig_v, 
 *
 *           Cases 1A, 1B, 1C all have orig_v < orig_y
 *           Case  1A: orig_v < orig_i < orig_y
 *           Case  1B: orig_i < orig_v < orig_y
 *           Case  1C: orig_v < orig_y < orig_i
 *
 *           Cases 2A, 2B, 2C all have orig_y < orig_v
 *           Case  2A: orig_y < orig_i < orig_v
 *           Case  2B: orig_i < orig_y < orig_v
 *           Case  2C: orig_y < orig_v < orig_i
 *
 *           See notes in cm_add_single_trans().
 * 
 * See related functions for explanation of parameters. 
 * Returns: (void) 
 */
static void
cm_subtract_subpaths(CM_t *orig_cm, CM_t *sub_cm, int **orig2sub_smap, int **sub2orig_smap, int orig_v, int orig_y,
		     int orig_i, int sub_v, int yoffset, double *orig_psi, char ***tmap)
{

  if(orig_i == -1)
    return;

  /* Determine which case we have for i1 */
  if(orig_v < orig_i && orig_i < orig_y)
    {
      //      printf("sub_cm->t[v:%d] [yoff:%d] case non-root 1A\n", sub_v, yoffset);
      /* this is handled correctly by cm_sum_subpaths() */
    }
  if(orig_i < orig_v && orig_v < orig_y)
    {
      //printf("sub_cm->t[v:%d] [yoff:%d] case non-root 1B\n", sub_v, yoffset);
      /* We've counted paths from
       * i -> v -> y correctly for sub_i -> sub_y 
       *           incorrectly for sub_v -> sub_y
       */
      sub_cm->t[sub_v][yoffset] -= orig_psi[orig_i] * 
	cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_i, orig_v, sub_v, tmap, orig_psi, -1, -2, -3 , -4) *
	cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_v, orig_y, sub_v, tmap, orig_psi, -1, -2, -3 , -4);
    }

  if(orig_v < orig_y && orig_y < orig_i)
    {
      Die("sub_cm->t[v:%d] [yoff:%d] case non-root 1C\n", sub_v, yoffset);
      /* impossible case */
    }


  if(orig_y < orig_i && orig_i < orig_v)
    {
      //printf("sub_cm->t[v:%d] [yoff:%d] case non-root 2A\n", sub_v, yoffset);      
      /* this is handled correctly by cm_sum_subpaths() */
    }
  if(orig_i < orig_y && orig_y < orig_v)
    {
      //printf("sub_cm->t[v:%d] [yoff:%d] case non-root 2B\n", sub_v, yoffset);
      /* We've counted paths from
       * i -> y -> v correctly for sub_i -> sub_y 
       *           incorrectly for sub_v -> sub_y
       */
      sub_cm->t[sub_v][yoffset] -= orig_psi[orig_i] * 
	cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_i, orig_y, sub_v, tmap, orig_psi, -1, -2, -3 , -4) *
	cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_y, orig_v, sub_v, tmap, orig_psi, -1, -2, -3 , -4);
    }
  if(orig_y < orig_v && orig_v < orig_i)
    {
      //printf("sub_cm->t[v:%d] [yoff:%d] case non-root 2C\n", sub_v, yoffset);
      /* We've counted paths from
       * y -> v -> i correctly for sub_i -> sub_y 
       *           incorrectly for sub_v -> sub_y
       */
      sub_cm->t[sub_v][yoffset] -= orig_psi[orig_y] * 
	cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_y, orig_v, sub_v, tmap, orig_psi, -1, -2, -3 , -4) *
	cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_v, orig_i, sub_v, tmap, orig_psi, -1, -2, -3 , -4);
    }
  return;
}

/**************************************************************************
 * EPN 08.31.06
 * Function: cm_sum_subpaths()
 *
 * Purpose:  Calculate probability of getting from one state (start) to 
 *           another (end) in a CM, taking special considerations. 
 * 
 *           This function is similar to CP9_cm2wrhmm::cm_sum_subpaths_cp9()
 *           but was written for mapping transitions in a template CM (orig_cm)
 *           to those in a sub CM built from that template.
 * 
 *           Sum the probability of all subpaths that start 
 *           at "start" and end at "end" (ignoring "end"->"end" and
 *           "start" -> "start" transitions if they exist)
 * 
 *           NOT SURE HOW TO PORT THIS:!
 *           When getting a transition probability for the sub_cm state v, we ignore the 
 *           contribution of subparses that correspond to other transitions out of 
 *           the same node that state v is in.
 *           For example, we don't want to include the probability of an 
 *           insert(node a) ->match (node a+1) sub parse when calculating
 *           the transition probability for match(node a) -> match (node a+1), 
 *           (CTMM) because CTIM maps to that transition. 
 *          
 * Args:    
 * CM_t *orig_cm     - the original, template CM
 * CM_t *sub_cm      - the sub CM
 * int start         - state index we're starting at
 * int end           - state index we're ending in
 * int sub_start     - sub_cm state index of start we're starting at
 * double *orig_psi       - for orig_cm: orig_psi[v] is the expected number of times state v is entered
 *                     in a CM parse
 * int orig_insert1  - index of an orig_cm insert state that maps to the same sub_cm node as
 *                     1 of the potentially 2 sub_cm states that orig_cm state 'start' maps to. 
 *                     -1 if none. There are up to 4 of these, up to 2 for up to 2 different
 *                     sub_cm nodes.
 * int orig_insert2
 * int orig_insert3
 * int orig_insert4
 * Returns: Float, the summed probability of all subpaths through the CM
 *          starting at "start" and ending at "end".
 */
static float
cm_sum_subpaths(CM_t *orig_cm, CM_t *sub_cm, int **orig2sub_smap, int start, int end, int sub_start, char ***tmap, double *orig_psi,
		int orig_insert1, int orig_insert2, int orig_insert3, int orig_insert4)
{
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
  int sub_end1; /* up to 2 sub_cm states that contain a state mapping */
  int sub_end2; /* to orig_cm state 'end'; -1 if none.*/
  int sub_v1;   /* up to 2 sub_cm states that contain a state mapping */
  int sub_v2;   /* to orig_cm state v; -1 if none.*/
  int skip;     /* TRUE to skip state v's contribution */

  int temp;

  //printf("\t\t\tin cm_sum_subpaths, start: %d | end: %d\n", start, end);
  
  if(start > end)
    {
      printf("ERROR in cm_sum_subpaths2: start: %d > end: %d\n", start, end);
      temp = end;
      end = start;
      start = temp;
    }
  if(start == end)
    {
      if(orig_cm->sttype[start] != IL_st && orig_cm->sttype[start] != IR_st)
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
	  if((orig_cm->stid[start] != MATP_MP && orig_cm->stid[start] != MATP_D) && 
	     (orig_cm->stid[start] != MATP_ML && orig_cm->stid[start] != MATP_MR))
	    {
	      printf("ERROR asking for self transition of non-insert, non-MATP state: %d\n", start);
	      exit(1);
	    }
	  return 1.0;
	}
      /* else we just return the self-insert probability */
      //printf("\t\t\tReturning self insert prob: %f\n", orig_cm->t[start][0]);
      return orig_cm->t[start][0];
    }
  to_return = 0.;
  
  sub_psi = malloc(sizeof(double) * (end - start + 1));

  /* Initialize sub_psi[0].*/
  sub_psi[0] = 1.; /* have to start in "start" */
  
  /* Check if we need to ignore the probability
   * mass from the ORIG_CM insert state(s) that maps to the SUB CM insert state(s) of this node .
   * These states are found in cm_add_single_trans() which calls this function, and are passed in 
   */
  if((orig_cm->sttype[start] != IL_st && orig_cm->sttype[start] != IR_st) &&
     (orig_cm->sttype[end]   != IL_st && orig_cm->sttype[end]   != IR_st))
    {
      insert_to_start = 0.;
      /*
      if(orig_insert1 != -1 && orig_insert1 < start) 
	insert_to_start += orig_psi[orig_insert1] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_insert1, start, tmap, orig_psi, -1, -1, -1, -1);
      if(orig_insert2 != -1 && orig_insert2 < start) 
	insert_to_start += orig_psi[orig_insert2] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_insert2, start, tmap, orig_psi, -1, -1, -1, -1);
      if(orig_insert3 != -1 && orig_insert3 < start) 
	insert_to_start += orig_psi[orig_insert3] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_insert3, start, tmap, orig_psi, -1, -1, -1, -1);
      if(orig_insert4 != -1 && orig_insert4 < start) 
	insert_to_start += orig_psi[orig_insert4] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_insert4, start, tmap, orig_psi, -1, -1, -1, -1);
      sub_psi[0] -= insert_to_start / orig_psi[start];
      */
      //printf("\t\tinsert_to_start: %f sub_psi[0]: %f\n", insert_to_start, sub_psi[0]);
    }
  /* note: when cm_sum_subpaths is called recursively above
   * it will never result in another recursive call, 
   * because its "start" is an insert state. This is why we can get away
   * with passing -1 for all: orig_insert1, orig_insert2, orig_insert3, orig_insert4
   */

  if(orig2sub_smap[end][0] != -1) sub_end1 = orig2sub_smap[end][0];
  else sub_end1 = -1;
  if(orig2sub_smap[end][1] != -1) sub_end2 = orig2sub_smap[end][1];
  else sub_end2 = -1;
  
  for (v = (start+1); v <= end; v++) 
    {
      sub_psi[v-start] = 0.;
      n_v = orig_cm->ndidx[v];
      if(orig_cm->sttype[v] == IL_st || orig_cm->sttype[v] == IR_st)
	is_insert = 1;
      else
	is_insert = 0;
      
      if(orig_cm->sttype[v] == S_st)
	{
	  /* previous state is necessarily either a BIF_B or a END_E, either
	   * way, there's no transitions FROM previous state to this state, so
	   * we handle this in a special way.*/
	  sub_psi[v-start] = sub_psi[(v-1)-start] * 1.;
	}
      
      if(orig2sub_smap[v][0] != -1) sub_v1 = orig2sub_smap[v][0];
      else sub_v1 = -1;
      if(orig2sub_smap[v][1] != -1) sub_v2 = orig2sub_smap[v][1];
      else sub_v2 = -1;
      
      //printf("\t\t\tv: %d, se1: %d se2: %d sv1: %d sv2: %d\n", v, sub_end1, sub_end2, sub_v1, sub_v2);

      /* We want to check if there is a sub_cm transition from the sub_cm state that 
       * maps to v (v_s) to end exists, if so then we want to skip the contribution here.
       * Unless v = end.
       */
      /*if(v != end && is_insert)*/
      skip = FALSE;
      if(v != end)
	{
	  if((is_insert) && ((sub_trans_check(sub_cm, orig2sub_smap, v, end,   sub_start)) ||
			     (sub_trans_check(sub_cm, orig2sub_smap, v, start, sub_start))))
	    {
	      /* if v is an insert state and there exists a transition from 
	       * sub_v -> sub_end or sub_end -> sub_v (where sub_v is a sub_cm state
	       * that maps to orig_cm state v, and sub_end is a sub_cm state that 
	       * maps to orig_cm state end), we skip the contribution because
	       * it will be counted when compiling counts for sub_v -> sub_end or 
	       * sub_end -> sub_v.
	       */
	      skip = TRUE;
	      //printf("\t\t\tskip true case 911\n");
	    }
	}
      if(!skip)
	{
	  for (y = orig_cm->pnum[v]-1; y >= is_insert; y--) 
	    {
	      x = orig_cm->plast[v] - y;
	      /* x is a parent of v, we're adding contribution 
	       * of transition from x to v. */
	      tmap_val = tmap[(int) orig_cm->stid[x]][(int) orig_cm->ndtype[orig_cm->ndidx[v]+is_insert]][(int) orig_cm->stid[v]];
	      if(tmap_val == -1)
		{
		  printf("tmap ERROR 1\n");
		  printf("v: %d | pnum[v]: %d | plast[v]: %d | y: %d | x: %d | d1: %d | d2: %d | d3: %d\n", v, orig_cm->pnum[v], orig_cm->plast[v], y, x, ((int) orig_cm->stid[x]), ((int) (orig_cm->ndtype[orig_cm->ndidx[v]+is_insert])), ((int) orig_cm->stid[v]));
		  exit(1);
		}
	      if((x - start) < 0)
		sub_psi[v-start] += 0.;
	      else
		sub_psi[v-start] += sub_psi[x-start] * orig_cm->t[x][(int) tmap_val];
	    }
	  if(v != end && is_insert) /* we don't want to include the probability of an 
				       insert self-loop to end itself */
	    {	  
	      sub_psi[v-start] += sub_psi[v-start] * (orig_cm->t[v][0] / (1-orig_cm->t[v][0])); 
	      /*the contribution of the self insertion loops*/
	    }
	}
      else
	{	
	  /*skip the contribution*/
	  //printf("v: %d | skipping the contribution\n", v);
	}
    }
  //printf("***returning from cm_sum_subpaths: %f\n", sub_psi[(end-start)]);

  to_return = (float) sub_psi[end-start];
  free(sub_psi);
  return to_return;
}
/**************************************************************************
 * EPN 09.01.06
 * debug_print_cm_params()
 *
 * Purpose:  Print out emission and transition probabilities and scores
 *           for a CM.
 *
 * Args:    
 * CM_t *cm     
 * Returns: (void) 
 */
void
debug_print_cm_params(CM_t *cm)
{
  int v, i;
  int yoffset;

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

  printf("cm->nodes: %d\n", cm->nodes);
  printf("cm->M:     %d\n", cm->M);
  for(v = 0; v < cm->M; v++)
    {
      printf("v:%4d:%4d %4s %2s\n", v, cm->ndidx[v], nodetypes[cm->ndtype[cm->ndidx[v]]], sttypes[cm->sttype[v]]);
      if(cm->sttype[v] == MP_st)
	{
	  printf("\tE: ");
	  for(i = 0; i < MAXABET*MAXABET; i++)
	    printf("%0.3f ", cm->e[v][i]);
	  printf("\n");
	}
      else if(cm->sttype[v] == ML_st ||
	      cm->sttype[v] == MR_st ||
	      cm->sttype[v] == IL_st ||
	      cm->sttype[v] == IR_st)
	{
	  printf("\tE: ");
	  for(i = 0; i < MAXABET; i++)
	    printf("%0.3f ", cm->e[v][i]);
	  printf("\n");
	}
      if(cm->sttype[v] != B_st && cm->sttype[v] != E_st)
	{
	  printf("\tT: ");
	  for(yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	    printf("%0.3f ", cm->t[v][yoffset]);
	  printf("\n");
	}	    
      else if(cm->sttype[v] == B_st)
	{
	  printf("\tL: %d | R: %d\n", cm->cfirst[v], cm->cnum[v]);
	}
      else if(cm->sttype[v] == E_st)
	printf("\n\n");
    }
  free(nodetypes);
  free(sttypes);
  printf("\n\n");
  return;
}

/**************************************************************************
 * EPN 09.01.06
 * Function: check_sub_cm_by_sampling()
 *
 * Purpose:  Given a CM and a sub CM that is supposed to mirror 
 *           the CM as closely as possible between two given consensus
 *           columns (spos and epos), check that the sub_cm was correctly 
 *           constructed. 
 *           
 *           The current approach is to build a CM Plan 9 HMM from the
 *           sub CM, then sample from the CM and see if the samples 
 *           were likely drawn from the CM Plan 9 distributions. 
 *           This is done inside CP9_cm2wrhmm::CP9_check_wrhmm_by_sampling().

 * Args:    
 * CM_t *orig_cm     - the original, template CM
 * CM_t  *sub_cm     - the sub CM built from the orig_cm
 * int spos          - first consensus column in cm that hmm models (often 1)
 * int epos          -  last consensus column in cm that hmm models 
 * float thresh      - P value threshold for chi-squared tests (often 0.01)
 * int nseq          - number of samples to draw from orig_cm
 * int *imp_cc;      - imp_cc[k] = 1 if CP9 node k is an impossible case to get 
 *		       the correct transition distros for the sub_cm. 
 * int *predict_ct   - number of times we predict we'll fail the test for
 *                     an HMM node, for each of the 7 cases.
 * int *wrong_predict_ct - number of times we predict we'll fail the test for
 *                         an HMM node and we get it right, for each of the 7 cases.
 * int print_flag    - TRUE to print useful info for debugging
 * Returns: TRUE: if CM and sub CM are "close enough" (see code)
 *          FALSE: otherwise
 */
int 
check_sub_cm_by_sampling(CM_t *orig_cm, CM_t *sub_cm, int spos, int epos, float thresh, int nseq,
			 int *imp_cc, int *predict_ct, int *wrong_predict_ct, int print_flag)
{
  struct cplan9_s       *sub_hmm; /* constructed CP9 HMM from the sub_cm */
  int ret_val;         /* return value */


  int *node_cc_left;    /* consensus column each CM node's left emission maps to
			 * [0..(cm->nodes-1)], -1 if maps to no consensus column*/
  int *node_cc_right;   /* consensus column each CM node's right emission maps to
			 * [0..(cm->nodes-1)], -1 if maps to no consensus column*/
  int *cc_node_map;     /* node that each consensus column maps to (is modelled by)
			 * [1..hmm_nmc] */
  int **cs2hn_map;      /* 2D CM state to HMM node map, 1st D - CM state index
		         * 2nd D - 0 or 1 (up to 2 matching HMM states), value: HMM node
		         * that contains state that maps to CM state, -1 if none.*/
  int **cs2hs_map;      /* 2D CM state to HMM node map, 1st D - CM state index
		         * 2nd D - 2 elements for up to 2 matching HMM states, 
		         * value: HMM STATE (0(M), 1(I), 2(D) that maps to CM state,
		         * -1 if none.
		         * For example: HMM node cs2hn_map[v][0], state cs2hs_map[v][0]
                         * maps to CM state v.*/
  int ***hns2cs_map;    /* 3D HMM node-state to CM state map, 1st D - HMM node index, 2nd D - 
		         * HMM state (0(M), 1(I), 2(D)), 3rd D - 2 elements for up to 
		         * 2 matching CM states, value: CM states that map, -1 if none.
		         * For example: CM states hsn2cs_map[k][0][0] and hsn2cs_map[k][0][1]
		         * map to HMM node k's match state.*/
  int debug_level;
  int v, k, ks;
  int cc;

  debug_level = 0;
  ret_val = TRUE;

  /* Build a CP9 HMM from the sub_cm */
  sub_hmm = AllocCPlan9((epos-spos+1));
  ZeroCPlan9(sub_hmm);
  CPlan9SetNullModel(sub_hmm, sub_cm->null, 1.0); /* set p1 = 1.0 which corresponds to the CM */

  map_consensus_columns(sub_cm, sub_hmm->M, &node_cc_left, &node_cc_right,
			&cc_node_map, debug_level);
  
  if(print_flag) printf("sub_hmm->M: %d\n", sub_hmm->M);
  CP9_map_cm2hmm_and_hmm2cm(sub_cm, sub_hmm, node_cc_left, node_cc_right, 
			    cc_node_map, &cs2hn_map, &cs2hs_map, 
			    &hns2cs_map, debug_level);
  
  if(!(CP9_cm2wrhmm(sub_cm, sub_hmm, node_cc_left, node_cc_right, cc_node_map, cs2hn_map,
		    cs2hs_map, hns2cs_map, debug_level)))
    Die("Couldn't build a CM Plan 9 HMM from the sub CM.\n");
  
  if(!(CP9_check_wrhmm_by_sampling(orig_cm, sub_hmm, spos, epos, hns2cs_map, thresh, nseq, imp_cc, 
				   predict_ct, wrong_predict_ct, print_flag)))
    Die("CM Plan 9 built from sub_cm fails sampling check using orig_cm; sub_cm was built incorrectly.!\n");
  else
    if(print_flag) printf("CM Plan 9 built from sub_cm passed sampling check; sub_cm was built correctly.\n");

  free(node_cc_left);
  free(node_cc_right);
  free(cc_node_map);
  for(v = 0; v <= sub_cm->M; v++)
    {
      free(cs2hn_map[v]);
      free(cs2hs_map[v]);
    }
  free(cs2hn_map);
  free(cs2hs_map);
  for(k = 0; k <= sub_hmm->M; k++)
    {
      for(ks = 0; ks < 3; ks++)
	free(hns2cs_map[k][ks]);
      free(hns2cs_map[k]);
    }
  free(hns2cs_map);
  FreeCPlan9(sub_hmm);
}

/**************************************************************************
 * EPN 09.06.06
 * Function: check_sub_cm_by_sampling2()
 *
 * Purpose:  Given a CM and a sub CM that is supposed to mirror 
 *           the CM as closely as possible between two given consensus
 *           columns (spos and epos), check that the sub_cm was correctly 
 *           constructed. 
 *           
 *           The approach is to sample from the CM and the sub_cm 
 *           and use those samples to build two CP9 HMMs, then 
 *           compare those two CP9 HMMs.
 *
 * Args:    
 * CM_t *orig_cm     - the original, template CM
 * CM_t  *sub_cm     - the sub CM built from the orig_cm
 * int spos          - first consensus column in cm that hmm models (often 1)
 * int epos          -  last consensus column in cm that hmm models 
 * int nseq          - number of sequences to sample to build the new HMMs.
 *
 * Returns: TRUE: if CM and sub CM are "close enough" (see code)
 *          FALSE: otherwise
 */
int 
check_sub_cm_by_sampling2(CM_t *orig_cm, CM_t *sub_cm, int spos, int epos, int nseq)
{
  struct cplan9_s       *orig_hmm; /* constructed CP9 HMM from the sub_cm */
  struct cplan9_s       *sub_hmm; /* constructed CP9 HMM from the sub_cm */
  int ret_val;         /* return value */
  Parsetree_t **tr;             /* Parsetrees of emitted aligned sequences */
  char    **dsq;                /* digitized sequences                     */
  char    **seq;                /* actual sequences (real letters)         */
  SQINFO            *sqinfo;    /* info about sequences (name/desc)        */
  MSA               *msa;       /* alignment */
  float             *wgt;
  int i;
  int idx;
  int L;
  int apos;
  int *matassign;
  int *useme;
  struct cp9trace_s **cp9_tr;   /* fake tracebacks for each seq            */
  int msa_nseq;                 /* this is the number of sequences per MSA,
				 * current strategy is to sample (nseq/nseq_per_msa)
				 * alignments from the CM, and add counts from
				 * each to the shmm in counts form (to limit memory)
				 */
  int nsampled;                 /* number of sequences sampled thus far */
  int debug_level;
  int cc;

  printf("\n\n*****************\nin check_sub_cm_by_sampling2()\n\n");

  debug_level = 0;
  ret_val = TRUE;
  msa_nseq = 1000;

  /* Build two CP9 HMMs */
  /* the orig_hmm only models consensus positions spos to epos of the orig_cm */
  orig_hmm = AllocCPlan9((epos-spos+1));
  ZeroCPlan9(orig_hmm);
  CPlan9SetNullModel(orig_hmm, orig_cm->null, 1.0); /* set p1 = 1.0 which corresponds to the CM */

  sub_hmm = AllocCPlan9((epos-spos+1));
  ZeroCPlan9(sub_hmm);
  CPlan9SetNullModel(sub_hmm, sub_cm->null, 1.0); /* set p1 = 1.0 which corresponds to the CM */

  /**************************************************/
  /* first sample from the orig_cm and use the samples to fill in orig_hmm
  /* sample MSA(s) from the CM */
  nsampled = 0;
  dsq    = MallocOrDie(sizeof(char *)             * msa_nseq);
  seq    = MallocOrDie(sizeof(char *)             * msa_nseq);
  tr     = MallocOrDie(sizeof(Parsetree_t)        * msa_nseq);
  sqinfo = MallocOrDie(sizeof(SQINFO)             * msa_nseq);
  wgt    = MallocOrDie(sizeof(float)              * msa_nseq);
  FSet(wgt, msa_nseq, 1.0);

  while(nsampled < nseq)
    {
      /*printf("nsampled: %d\n", nsampled);*/
      if(nsampled != 0)
	{
	  /* clean up from previous MSA */
	  MSAFree(msa);
	  free(matassign);
	  free(useme);
	  for (i = 0; i < msa_nseq; i++)
	    {
	      CP9FreeTrace(cp9_tr[i]);
	      FreeParsetree(tr[i]);
	      free(dsq[i]);
	      free(seq[i]);
	    }
	  free(cp9_tr);
	}
      /* Emit msa_nseq parsetrees from the CM */
      if(nsampled + msa_nseq > nseq)
	msa_nseq = nseq - nsampled;
      for (i = 0; i < msa_nseq; i++)
	{
	  EmitParsetree(orig_cm, &(tr[i]), &(seq[i]), &(dsq[i]), &L);
	  sprintf(sqinfo[i].name, "seq%d", i+1);
	  sqinfo[i].len   = L;
	  sqinfo[i].flags = SQINFO_NAME | SQINFO_LEN;
	}
      /* Build a new MSA from these parsetrees */
      msa = Parsetrees2Alignment(orig_cm, dsq, sqinfo, wgt, tr, msa_nseq, TRUE);
      
      /* Truncate the alignment prior to consensus column spos and after 
	 consensus column epos */
      useme = (int *) MallocOrDie (sizeof(int) * (msa->alen+1));
      for (apos = 0, cc = 0; apos < msa->alen; apos++)
	{
	  /* Careful here, placement of cc++ increment is impt, 
	   * we want all inserts between cc=spos-1 and cc=spos,
	   * and between cc=epos and cc=epos+1.
	   */
	  if(cc < (spos-1) || cc > epos)
	    useme[apos] = 0;
	  else
	    useme[apos] = 1;
	  if (!isgap(msa->rf[apos])) 
	    { 
	      cc++; 
	      if(cc == (epos+1))
		 useme[apos] = 0; 
		 /* we misassigned this guy, overwrite */ 
	    }
	}

      MSAShorterAlignment(msa, useme);

      /* Shorten the dsq's */
      for (i = 0; i < msa_nseq; i++)
	{
	  MakeDealignedString(msa->aseq[i], msa->alen, msa->aseq[i], &(seq[i])); 
	  free(dsq[i]);
	  dsq[i] = DigitizeSequence(seq[i], strlen(seq[i]));
	}
	  
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
	CP9TraceCount(orig_hmm, dsq[idx], msa->wgt[idx], cp9_tr[idx]);
      }
      nsampled += msa_nseq;
    }

  /*Next, renormalize the orig_hmm and logoddisfy it */
  CPlan9Renormalize(orig_hmm);
  CP9Logoddsify(orig_hmm);

  /* clean up from previous MSA */
  MSAFree(msa);
  free(matassign);
  free(useme);
  for (i = 0; i < msa_nseq; i++)
    {
      CP9FreeTrace(cp9_tr[i]);
      FreeParsetree(tr[i]);
      free(dsq[i]);
      free(seq[i]);
    }
  free(cp9_tr);
  
  /**************************************************/
  /* Now for the sub_hmm */
  /* first sample from the sub_cm and use the samples to fill in sub_hmm
  /* sample MSA(s) from the CM */
  nsampled = 0;
  dsq    = MallocOrDie(sizeof(char *)             * msa_nseq);
  seq    = MallocOrDie(sizeof(char *)             * msa_nseq);
  tr     = MallocOrDie(sizeof(Parsetree_t)        * msa_nseq);
  sqinfo = MallocOrDie(sizeof(SQINFO)             * msa_nseq);
  wgt    = MallocOrDie(sizeof(float)              * msa_nseq);
  FSet(wgt, msa_nseq, 1.0);

  while(nsampled < nseq)
    {
      /*printf("nsampled: %d\n", nsampled);*/
      if(nsampled != 0)
	{
	  /* clean up from previous MSA */
	  MSAFree(msa);
	  free(matassign);
	  for (i = 0; i < msa_nseq; i++)
	    {
	      CP9FreeTrace(cp9_tr[i]);
	      FreeParsetree(tr[i]);
	      free(dsq[i]);
	      free(seq[i]);
	    }
	  free(cp9_tr);
	}
      /* Emit msa_nseq parsetrees from the CM */
      if(nsampled + msa_nseq > nseq)
	msa_nseq = nseq - nsampled;
      for (i = 0; i < msa_nseq; i++)
	{
	  EmitParsetree(sub_cm, &(tr[i]), &(seq[i]), &(dsq[i]), &L);
	  sprintf(sqinfo[i].name, "seq%d", i+1);
	  sqinfo[i].len   = L;
	  sqinfo[i].flags = SQINFO_NAME | SQINFO_LEN;
	}
      /* Build a new MSA from these parsetrees */
      msa = Parsetrees2Alignment(sub_cm, dsq, sqinfo, wgt, tr, msa_nseq, TRUE);
      
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
	CP9TraceCount(sub_hmm, dsq[idx], msa->wgt[idx], cp9_tr[idx]);
      }
      nsampled += msa_nseq;
    }

  /*Next, renormalize the sub_hmm and logoddisfy it */
  CPlan9Renormalize(sub_hmm);
  CP9Logoddsify(sub_hmm);
  /**************************************************/

  printf("PRINTING SAMPLED ORIG HMM PARAMS:\n");
  debug_print_cp9_params(orig_hmm);
  printf("DONE PRINTING SAMPLED ORIG HMM PARAMS:\n");


  printf("PRINTING SAMPLED SUB HMM PARAMS:\n");
  debug_print_cp9_params(sub_hmm);
  printf("DONE PRINTING SAMPLED SUB HMM PARAMS:\n");

  FreeCPlan9(orig_hmm);
  FreeCPlan9(sub_hmm);
  return TRUE;
}


/**************************************************************************
 * EPN 09.07.06
 * check_orig_psi_vs_sub_psi()
 *
 * Purpose:  Check that the psi values for an original, template CM and 
 *           a sub CM are withing a given leeway threshold, given maps
 *           from the states of the original to the sub and vice versa.
 *
 * Args:    
 * CM_t *orig_cm     - the original, template CM
 * CM_t  *sub_cm     - the sub CM
 * double *orig_psi  - psi[v] is expected number of times CM state v is entered
 * double  *sub_psi
 * int **orig2sub_smap - 2D state map from orig_cm (template) to sub_cm.
 *                       1st dimension - state index in orig_cm 
 *                       2nd D - 2 elements for up to 2 matching sub_cm states, 
 * int **sub2orig_smap - 2D state map from orig_cm (template) to sub_cm.
 *                       1st dimension - state index in sub_cm (0..sub_cm->M-1)
 *                       2nd D - 2 elements for up to 2 matching orig_cm states, 
 * double threshold  - the threshold that mapping (potentially summed) psi 
 *                     values are allowed to be different by, without throwing an error.
 * int print_flag    - TRUE to print out the values, FALSE not to 
 * Returns: TRUE: if orig_cm and sub_cm are "close enough" (see code)
 *          FALSE: otherwise
 */
static int
check_orig_psi_vs_sub_psi(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, double *sub_psi, 
			  int **orig2sub_smap, int **sub2orig_smap, double threshold, 
			  int print_flag)
{
  int v_s; /* sub_cm state index*/ 
  int v_o; /* orig_cm state index*/ 

  double temp_psi;
  int violation;
  int v_ct; /* Number of violations */
  int ret_val; /* return value */
  int dual_mapping_insert;
  int dead_insert;
  int is_insert;
  int root_v_ct;
  float diff;

  if(print_flag)
    {
      printf("Printing psi in check_orig_psi_vs_sub_psi():\n");
      
      for(v_o = 0; v_o < orig_cm->M; v_o++)
	printf("orig_psi[%4d]: %.6f\n", v_o, orig_psi[v_o]);
      
    }

  ret_val = TRUE;
  v_ct = 0;
  root_v_ct = 0;
  if(print_flag == TRUE)
    printf("\n");
  for(v_s = 0; v_s < sub_cm->M; v_s++)
    {
      dual_mapping_insert = FALSE;
      dead_insert         = FALSE;

      if(sub_cm->sttype[v_s] == IL_st || sub_cm->sttype[v_s] == IR_st)
	is_insert = TRUE;
      else
	is_insert = FALSE;

      if(print_flag)
	printf("\tv_s: %4d (%.6f) ", v_s, sub_psi[v_s]);
      
      v_o = sub2orig_smap[v_s][0];

      if(sub_cm->sttype[v_s+1] == E_st)
	{
	  if(orig2sub_smap[v_o][1] != -1)
	    dead_insert = TRUE;
	}

      if(v_o != -1)
	{
	  if(print_flag)
	    printf("v_o1: %4d (%.6f) ", v_o, orig_psi[v_o]);
	  
	  temp_psi = orig_psi[v_o];
	  v_o = sub2orig_smap[v_s][1];
	  if(v_o != -1)
	    {
	      if(is_insert) /* this insert state maps to 2 orig_cm inserts */
		dual_mapping_insert = TRUE;

	      temp_psi += orig_psi[v_o];
	      /* NOTE: if is_insert, this is technically incorrect. 
	       * Summing the psi of the two orig_cm states that sub_cm v_s
	       * maps to does not give the expected number of times either
	       * is visited. This is because they are not independent, some
	       * paths can go through both of these insert states.
	       * However, if v_s is not an insert, then this is correct b/c
	       * only other cases involve MATP nodes (where a MATL_ML 
	       * maps to MATP_MP AND MATP_ML for example), in which
	       * the two orig_cm states are in the same split set and
	       * thus no one path can visit them both.
	       */

	      if(print_flag)
		printf("v_o2: %4d (%.6f)\n", v_o, orig_psi[v_o]);
	    }
	  else
	    if(print_flag) printf("\n");
	  
	  violation = FALSE;
	  if(dead_insert)
	    temp_psi = 0.0;
	  
	  /* 10.20.06 Found an exceedingly rare case (2 cases in all possible sub CM 
	   * models of RMARK) where psi test fails with diff < 0.00002 but > 0.00001
	   * (our default). Both cases involve insert self loops with p > 0.9, 
	   * the reason (I'm pretty sure) these guys fail is because when the
	   * contribution of the self insertion loop is included in fill_psi, even
	   * if the self-insert probs for a sub_cm and orig_cm state are equal, if
	   * the psi values for that state BEFORE the contribution of the self insert
	   * is added are > 0.0000001 or so, the self insert contributions amplifies
	   * that difference above our 0.00001 threshold. This only happens if the
	   * self insert probs are really high and explains the rareness of this case.
	   * The approach to fixing it is to subtract out the self loop contribution
	   * prior to checking if we exceed our threshold. 
	   */
	  if(is_insert && !dead_insert)
	    {
	      diff = (sub_psi[v_s] - (sub_psi[v_s] * sub_cm->t[v_s][0])) -
		(temp_psi - (temp_psi * orig_cm->t[sub2orig_smap[v_s][0]][0]));
	    }
	  else
	    diff = sub_psi[v_s] - temp_psi;
	  if((diff > threshold) || ((-1. * diff) > threshold))
	    {
	      if(!dual_mapping_insert)
		{
		  violation = TRUE;
		  v_ct++;
		  if((sub_cm->ndidx[v_s] == 0) || (sub_cm->ndidx[v_s] == 1))
		    {
		      root_v_ct++;
		    }
		}
	    }
	  if(violation)
	    printf("sub: %.6f | orig: %.6f | diff: %.6f VIOLATION\n\n", sub_psi[v_s], temp_psi, (sub_psi[v_s]-temp_psi));
	  else if(dead_insert && print_flag)
	    printf("sub: %.6f | orig: %.6f | diff: %.6f (DEAD INSERT)\n\n", sub_psi[v_s], temp_psi, (sub_psi[v_s]-temp_psi));
	  else if(dual_mapping_insert && print_flag)
	    printf("sub: %.6f | orig: %.6f | diff: %.6f (DUAL INSERT)\n\n", sub_psi[v_s], temp_psi, (sub_psi[v_s]-temp_psi));
	  else if(print_flag)
	    printf("sub: %.6f | orig: %.6f | diff: %.6f\n\n", sub_psi[v_s], temp_psi, (sub_psi[v_s]-temp_psi));
	}
      else
	{
	  if(sub_cm->sttype[v_s] != E_st &&
	     sub_cm->sttype[v_s] != B_st &&
	     sub_cm->sttype[v_s] != S_st &&
	     sub_cm->sttype[v_s] != EL_st)
	    Die("ERROR state v_s:%d maps to nothing and its not E,B,S,EL\n", v_s);
	  if(print_flag) printf("E B S or EL\n");
	}
    }

  if(v_ct > 0)
    {
      printf("ERROR, v_ct: %d sub_cm states violate the %f threshold b/t psi.\n", v_ct, threshold);
      /*exit(1);*/
      ret_val = FALSE;
    }
  else
    if(print_flag) printf("KACHOW! v_ct is 0 with thresh: %f!\n", threshold);

  if(root_v_ct > 0)
    printf("ROOT v_ct is %d with thresh: %f!\n", root_v_ct, threshold);

  return ret_val;
}

/**************************************************************************
 * EPN 09.11.06
 * orig2sub_state_check()
 *
 * Return true if either of the two sub_cm states that potentially map to 
 * orig_cm state orig_a are identical to either of the two
 * sub_cm states that map to sub_cm state orig_b.
 *
 * Args:    
 * int **orig2sub_smap
 * int   orig_a; 
 * int   orig_b; 
 * Returns: (v_oid) 
 */
static int
orig2sub_state_check(int **orig2sub_smap, int orig_a, int orig_b)
{
  int sub_a1, sub_a2, sub_b1, sub_b2;
  
  printf("in orig2sub_state_check() orig_a: %d orig_b: %d\n", orig_a, orig_b);

  sub_a1 = orig2sub_smap[orig_a][0];
  sub_a2 = (orig2sub_smap[orig_a][1] != -1) ? orig2sub_smap[orig_a][1] : -2;

  sub_b1 = (orig2sub_smap[orig_b][0] != -1) ? orig2sub_smap[orig_b][0] : -3;
  sub_b2 = (orig2sub_smap[orig_b][1] != -1) ? orig2sub_smap[orig_b][1] : -4;

  if((sub_a1 == sub_b1 || sub_a1 == sub_b2) ||
     (sub_a2 == sub_b1 || sub_a2 == sub_b2))
    return TRUE;
  
  /* else */
  return FALSE;
}


/**************************************************************************
 * EPN 09.11.06
 * sub_trans_check()
 *
 * Return TRUE if (1) either of the two sub_cm states that map to 
 * orig_cm state orig_a are either children OR parents of either of the two
 * sub_cm states that map to sub_cm state orig_b and (2) the sub_cm state
 * of interest (return_sub_v) is >= min_sub_v.
 *
 * Args:    
 * CM_t  sub_cm,
 * int **orig2sub_smap
 * int   orig_a; 
 * int   orig_b; 
 * int   min_sub_v;
 * Returns: TRUE or FALSE;
 */
static int
sub_trans_check(CM_t *sub_cm, int **orig2sub_smap, int orig_a, int orig_b,
		int min_sub_v)
{
  int sub_a1, sub_a2, sub_b1, sub_b2;
  float diff;
  
  //printf("in sub_trans_check() orig_a: %d orig_b: %d min_sub_v: %d\n", orig_a, orig_b, min_sub_v);

  sub_a1 = orig2sub_smap[orig_a][0];
  sub_a2 = orig2sub_smap[orig_a][1];
  sub_b1 = orig2sub_smap[orig_b][0];
  sub_b2 = orig2sub_smap[orig_b][1];

  if(trans_check_helper(sub_cm, sub_a1, sub_b1) >= min_sub_v) return TRUE;
  //if(trans_check_helper(sub_cm, sub_b1, sub_a1) >= min_sub_v) return TRUE;
  if(trans_check_helper(sub_cm, sub_a1, sub_b2) >= min_sub_v) return TRUE;
  //if(trans_check_helper(sub_cm, sub_b2, sub_a1) >= min_sub_v) return TRUE;
  if(trans_check_helper(sub_cm, sub_a2, sub_b1) >= min_sub_v) return TRUE;
  //if(trans_check_helper(sub_cm, sub_b1, sub_a2) >= min_sub_v) return TRUE;
  if(trans_check_helper(sub_cm, sub_a2, sub_b2) >= min_sub_v) return TRUE;
  //if(trans_check_helper(sub_cm, sub_b2, sub_a2) >= min_sub_v) return TRUE;

  return FALSE;
}  

/**************************************************************************
 * EPN 09.11.06
 * trans_check_helper()
 *
 * Return TRUE if there's a transition in the CM from state a to state b.
 *
 * Args:    
 * CM_t  cm,
 * int   a; 
 * int   b; 
 * Returns: a if b is a child of a (a->b exists)
 *          -1 otherwise
 */
int
trans_check_helper(CM_t *cm, int a, int b)
{
  //  printf("\t**in trans_check_helper a: %d | b: %d\n", a, b);

  if((a == -1 || b == -1) || (a > b))
    return -1;

  if((b - cm->cfirst[a]) < cm->cnum[a])
    { /*printf("returning a: %d\n", a)*/; return a; }

  return -1;

}


/**************************************************************************
 * EPN 09.14.06
 * Function: cm2sub_cm_subtract_root_subpaths
 *
 * Purpose:  When building a sub CM (sub_cm) from an original, template
 *           CM (orig_cm), there's special considerations that must be
 *           taken involving sub_cm nodes with 2 insert states where there
 *           isn't an exactly identical node in the original CM. The only
 *           node that potentially meets this criteria is the ROOT_nd because
 *           all MATP nodes in sub_cm must necessarily also exist in the
 *           orig_cm. 
 *           
 *           The situation we have to fix is when the sub_cm ROOT_IL and
 *           sub_cm ROOT_IR map to orig_cm states in two different nodes.
 *           In this case, the probabilities out of the sub_cm ROOT_IL and       
 *           ROOT_IR will be incorrectly calculated using the standard
 *           method (cm_sum_subpaths() based). We have to go back
 *           and subtract out some subpaths that will have been counted 
 *           twice using that method.

 * Args:    
 * CM_t *orig_cm     - the original, template CM
 * CM_t *sub_cm      - the sub CM
 * int start         - state index we're starting at
 * int end           - state index we're ending in
 * int sub_start     - sub_cm state index of start we're starting at
 * double *orig_psi  - for orig_cm: orig_psi[v] is the expected number of times state v is entered
 *                     in a CM parse
 * int orig_insert1  - index of an orig_cm insert state that maps to the same sub_cm node as
 *                     1 of the potentially 2 sub_cm states that orig_cm state 'start' maps to. 
 *                     -1 if none. There are up to 4 of these, up to 2 for up to 2 different
 *                     sub_cm nodes.
 * int orig_insert2
 * int orig_insert3
 * int orig_insert4
 * Returns: Float, the summed probability of all subpaths through the CM
 *          starting at "start" and ending at "end".
 */
static int
cm2sub_cm_subtract_root_subpaths(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, char ***tmap, int sub_nd, int **orig2sub_smap, int **sub2orig_smap)

{
  int orig_il1, orig_il2, orig_ir1, orig_ir2;
  double temp_psi1;
  double temp_psi2;
  int ir1_case;
  int ir2_case;
  int is_matp;
  int sub_il;
  int sub_ir;
  int yoffset;
  int sub_y;
  int curr_case;
  int orig_root_s;
  int orig_y1;
  int orig_y2;
  float diff;
  float contribution1;
  float contribution2;
  int   start, end;
  float sum;

  /* Most of the work here is done in a helper function, which checks to see which of the 
   * 3 possible cases we have. There are 4 calls to this helper function because their
   * be 1 to 4 instances of each case, 1 if sub_il and sub_ir both map to exactly 1 orig_cm state,
   * up to 4, if they both map to 2 orig_cm states. 
   * Case 1 occurs if orig_il* < orig_ir*
   * Case 2 occurs if orig_il* > orig_ir*
   * If there are less than 4 calls (at least one of
   * sub_cm ROOT_IL and sub_cm ROOT_IR map to a single orig_cm state), then it is possible that
   * their can be instances of both cases (orig_ir1 < orig_il1 < orig_ir2 for example.)
   * 
   * There are however 2 guarantees we can exploit: 
   *
   * 1. If there  are 4 orig_cm states (both sub_il and sub_ir map to 2 orig_cm states), then it
   * MUST be the case that only one case exists (both orig_il's come before both orig_ir's (case 1)
   * or both orig_il's come after both orig_ir's (case 2)).
   * This is because the 4 state situation only occurs if both sub_il and sub_ir are dual mapping
   * inserts that map to an ambiguous insert column in the orig_cm. Because there is only one 
   * of these per 'stem' (more exactly, one per END_E state), there cannot be a case where 
   * both cases 1 and 2 exist in such a situation. If there is, then I misunderstand CM
   * architecture. We rely on this guarantee in the code below, but check for a violation
   * to it, and die if one occurs.
   * 
   * 2. If sub_il maps to exactly 2 orig_il states (orig_il1 and orig_il2) and sub_ir maps
   * to exactly 1 orig_ir state (orig_ir1) then either both orig_il states come before orig_ir1 (case 1),
   * or both orig_il states come after orig_ir1. This is because in this case the sub_il state
   * is a dual mapping insert corresponding to an ambiguous insert column in the orig_cm
   * and it's impossible for the orig_il1 and orig_il2 to 'straddle' orig_ir - it just is.
   * We rely on this guarantee in the code, but we check for it and die if its violated.
   */


  sub_il = 1; /* sub ROOT_IL index */
  sub_ir = 2; /* sub ROOT_IR index */

  orig_il1 = sub2orig_smap[sub_il][0];
  orig_il2 = sub2orig_smap[sub_il][1];
  
  orig_ir1 = sub2orig_smap[sub_ir][0];
  orig_ir2 = sub2orig_smap[sub_ir][1];

  printf("in cm2sub_cm_subtract_root_subpaths orig_il1: %d | orig_il2: %d | orig_ir1: %d | orig_ir2: %d\n", orig_il1, orig_il2, orig_ir1, orig_ir2);

  /* check to make sure our guarantee (1) when we have 4 states mentioned above is met */
  if((orig_il1 != -1 && orig_il2 != -1) &&
     (orig_ir1 != -1 && orig_ir1 != -1))
    {
      if((orig_il1 < orig_ir1 && orig_il1 < orig_ir2) &&
	 (orig_il2 < orig_ir1 && orig_il2 < orig_ir2))
	{ /* okay */ }
      if((orig_il1 > orig_ir1 && orig_il1 > orig_ir2) &&
	 (orig_il2 > orig_ir1 && orig_il2 > orig_ir2))
	{ /* okay */ }
      else
	{ /* NOT okay */ 
	  Die("ERROR in cm2sub_cm_subtract_root_subpaths, both sub_il and sub_ir are dual mapping inserts, but both cases 1 and 2 exist! You odn't understand CM architecure\n");
	}
    }

  /* check to make sure our guarantee (2) mentioned above when we have 2 orig_il states and 1 orig_ir state
   * mentioned above is met */
  if((orig_il1 != -1 && orig_il2 != -1) &&
     (orig_ir1 != -1 && orig_ir2 == -1))
    {
      if(orig_il1 < orig_ir1 && orig_il2 < orig_ir1)
	{ /* okay */ }
      if(orig_il1 > orig_ir1 && orig_il2 > orig_ir1)
	{ /* okay */ }
      else
	{ /* NOT okay */ 
	  Die("ERROR in cm2sub_cm_subtract_root_subpaths, sub_il is a dual mapping inserts, sub_ir is not, but both cases 1 and 2 exist! You odn't understand CM architecure\n");
	}
    }

  ir1_case = cm2sub_cm_subtract_root_subpaths_helper(orig_cm, sub_cm, orig_psi, tmap, 
						     orig2sub_smap, sub2orig_smap, orig_il1, orig_ir1);
  if(ir1_case == 0)
    Die("ERROR in cm2sub_cm_subtract_root_subpaths, neither case 1 or case 2\n");
  curr_case = cm2sub_cm_subtract_root_subpaths_helper(orig_cm, sub_cm, orig_psi, tmap, 
						      orig2sub_smap, sub2orig_smap, orig_il2, orig_ir1);
  if((ir1_case == 1 && curr_case == 2) ||
     (ir1_case == 2 && curr_case == 1))
    {
      /* This should never happen it violates guarantee 2 above. */
      Die("ERROR in cm2sub_cm_subtract_root_subpaths, sub_il is a dual mapping inserts, sub_ir is not, but both cases 1 and 2 exist! You odn't understand CM architecure\n Further, your guarantee 2 check is flawed.");
    }

  if(ir1_case == 1)
    {
      /* we have to adjust counts out of ROOT_S, because we're double counting those paths
       * that start at ROOT_S and go through sub_cm ROOT_IR to split states in node 1, we 
       * have counted them once in the virtual counts out of ROOT_S and again in the virtual
       * counts out of sub_cm ROOT_IR.
       */
      orig_root_s = 0; /* index of orig_cm ROOT_S and sub_cm ROOT_S */
      cm2sub_cm_subtract_root_subpaths_helper(orig_cm, sub_cm, orig_psi, tmap,  
					      orig2sub_smap, sub2orig_smap, orig_root_s, orig_ir1);
    }

  ir2_case = cm2sub_cm_subtract_root_subpaths_helper(orig_cm, sub_cm, orig_psi, tmap,  
						     orig2sub_smap, sub2orig_smap, orig_il1, orig_ir2);
  curr_case = cm2sub_cm_subtract_root_subpaths_helper(orig_cm, sub_cm, orig_psi, tmap,  
						      orig2sub_smap, sub2orig_smap, orig_il2, orig_ir2);
  if((ir2_case == 1 && curr_case == 2) ||
     (ir2_case == 2 && curr_case == 1))
    {
      /* This should never happen it violates guarantee 1 above. */
      Die("ERROR in cm2sub_cm_subtract_root_subpaths, both sub_il and sub_ir are dual mapping inserts, but both cases 1 and 2 exist! You don't understand CM architecure. Further your guarantee 1 check is flawed.\n");
    }
  if(ir2_case == 1)
    {
      /* we have to adjust counts out of ROOT_S, because we're double counting those paths
       * that start at ROOT_S and go through sub_cm ROOT_IR to split states in node 1, we 
       * have counted them once in the virtual counts out of ROOT_S and again in the virtual
       * counts out of sub_cm ROOT_IR.
       */
      orig_root_s = 0; /* index of orig_cm ROOT_S and sub_cm ROOT_S */
      cm2sub_cm_subtract_root_subpaths_helper(orig_cm, sub_cm, orig_psi, tmap,  
					      orig2sub_smap, sub2orig_smap, orig_root_s, orig_ir2);
    }

  if(ir1_case == 2 || ir2_case == 2)
    {
      /* adjust the virtual counts for transitions out of ROOT_IL */
      printf("KACHOW! case2 FLAG\n");
      
      printf("after t[0][1] = %f\n", sub_cm->t[0][1]);

      /* First determine temp_psi1 and temp_psi2, where: 
       * temp_psi1 = what orig_psi[orig_il1] would be w/o the subpath of orig_ir1 -> orig_il1 and orig_ir2 -> orig_il1
       * temp_psi2 = what orig_psi[orig_il2] would be w/o the subpath of orig_ir1 -> orig_il2 and orig_ir2 -> orig_il2
       */
      
      /* This takes several steps: */
      /* first remove contribution of self-insert */
      temp_psi1   = orig_psi[orig_il1] / (1 + (orig_cm->t[orig_il1][0]/(1-(orig_cm->t[orig_il1][0]))));
      if(orig_il2 != -1)
	temp_psi2 = orig_psi[orig_il2] / (1 + (orig_cm->t[orig_il2][0]/(1-(orig_cm->t[orig_il2][0]))));
      
      printf("0 temp_psi1: %f\n", temp_psi1);
      if(orig_il2 != -1) printf("0 temp_psi2: %f\n", temp_psi2);
      
      /* subtract contribution of orig_ir1 and orig_ir2 to orig_il1 and orig_i12 */
      if(ir1_case == 2) 
	temp_psi1   -= orig_psi[orig_ir1] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ir1, orig_il1, sub_ir, tmap, orig_psi, -1, -2, -3, -4); 
      if(ir2_case == 2)
	temp_psi1 -= orig_psi[orig_ir2] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ir2, orig_il1, sub_ir, tmap, orig_psi, -1, -2, -3, -4); 
      printf("1 temp_psi1: %f\n", temp_psi1);
      
      if(orig_il2 != -1)
	{
	  if(ir1_case == 2)
	    temp_psi2   -= orig_psi[orig_ir1] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ir1, orig_il2, sub_ir, tmap, orig_psi, -1, -2, -3, -4); 
	  if(ir2_case == 2)
	    temp_psi2 -= orig_psi[orig_ir2] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ir2, orig_il2, sub_ir, tmap, orig_psi, -1, -2, -3, -4); 
	  printf("1 temp_psi2: %f\n", temp_psi2);
	}
      
      /* add contribution of self-insert back */
      temp_psi1 += temp_psi1 * (orig_cm->t[orig_il1][0] / (1 - orig_cm->t[orig_il1][0]));
      if(orig_il2 != -1)
	temp_psi2 += temp_psi2 * (orig_cm->t[orig_il2][0] / (1 - orig_cm->t[orig_il2][0]));
      
      printf("2 temp_psi1: %f\n", temp_psi1);
      if(orig_il2 != -1)
	printf("2 temp_psi2: %f\n", temp_psi2);
      
      /* Finally scale the transitions out of ROOT_IL to the split states of node 1 by the appropriate scaling factor.
       * if orig_il2 == -1 this is easy, the scaling factor is simply: (temp_psi / orig_psi[orig_il1])
       * if orig_il2 != -1 its more complicated because we have to recalculate the relative contribution of orig_il1 and
       *                   orig_il2 towards the virtual counts out of sub_il. 
       */
      if(orig_il2 == -1) /* easier case */
	{
	  for(yoffset = 0; yoffset < sub_cm->cnum[sub_il]; yoffset++)
	    {
	      sub_y = sub_cm->cfirst[sub_il] + yoffset; 
	      if(sub_cm->ndidx[sub_y] != sub_cm->ndidx[sub_il])
		sub_cm->t[sub_il][yoffset] *= (temp_psi1/orig_psi[orig_il1]);
	      /* we only want to adjust the virtual counts into split states of the next node */
	    }
	}
      else /* harder case */
	{
	  printf("KACHOW!!!! HARDER ROOT_S CASE 2 CASE!\n");
	  for(yoffset = 0; yoffset < sub_cm->cnum[sub_il]; yoffset++)
	    {
	      sub_y = sub_cm->cfirst[sub_il] + yoffset; 

	      if(sub_cm->ndidx[sub_y] != sub_cm->ndidx[sub_il])
		{
		  /* first we need to find the one or two orig_cm states that map to sub_y */
		  orig_y1 = sub2orig_smap[sub_y][0];	      
		  orig_y2 = sub2orig_smap[sub_y][1];
		  
		  /* now recalculate subpath prob from orig_il1 -> orig_y1 and orig_il1 -> orig_y2*/
		  if(orig_il1 < orig_y1)
		    {
		      start = orig_il1;
		      end   = orig_y1;
		    }
		  else
		    {
		      start = orig_y1;
		      end   = orig_il1;
		    }
		  contribution1 = orig_psi[start] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, start, end, sub_il, tmap, orig_psi, -1, -2, -3, -4); 
		  if(orig_y2 != -1)
		    {
		      if(orig_il1 < orig_y2)
			{
			  start = orig_il1;
			  end   = orig_y2;
			}
		      else
			{
			  start = orig_y2;
			  end   = orig_il1;
			}
		      contribution1 += orig_psi[start] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, start, end, sub_il, tmap, orig_psi, -1, -2, -3, -4); 
		    }
		  
		  /* recalculate subpath prob from orig_il2 -> orig_y1 and orig_il2 -> orig_y2 */
		  if(orig_il2 < orig_y1)
			{
			  start = orig_il2;
			  end   = orig_y1;
			}
		      else
			{
			  start = orig_y1;
			  end   = orig_il2;
			}
		  contribution2 = orig_psi[start] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, start, end, sub_il, tmap, orig_psi, -1, -2, -3, -4); 

		  if(orig_y2 != -1)
		    {
		      if(orig_il2 < orig_y2)
			{
			  start = orig_il2;
			  end   = orig_y2;
			}
		      else
			{
			  start = orig_y2;
			  end   = orig_il2;
			}
		      contribution2 += orig_psi[start] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, start, end, sub_il, tmap, orig_psi, -1, -2, -3, -4); 
		    }
		  
		  diff = (contribution1 + contribution2) - sub_cm->t[sub_il][yoffset];
		  if((diff >= 0 && diff >  0.0001) ||
		     (diff <= 0 && diff < -0.0001))
		    Die("ERROR in cm2sub_cm_subtract_root_subpaths() contribution1:%f and contribution2%f don't add up to sub_cm->t[1][yoffset:%d]\n", contribution1, contribution2, yoffset, sub_cm->t[sub_il][yoffset]);
		  
		  sum = contribution1 + contribution2;
		  contribution1 /= sum;
		  contribution2 /= sum;
		  
		  /* finally adjust the counts */
		  sub_cm->t[sub_il][yoffset] *= ((contribution1 * (temp_psi1/orig_psi[orig_il1])) + 
						 (contribution2 * (temp_psi2/orig_psi[orig_il2])));

		  /* we only want to adjust the virtual counts into split states of the next node */
		}
	    }
	}
    }
  return;
}

/**************************************************************************
 * EPN 09.15.06
 * Function: cm2sub_cm_subtract_root_subpaths_helper()
 *
 * Purpose:  Helper function that actually does the subpath subtraction
 *           for the ROOT_nd of sub CMs. It first checks which of 3 cases
 *           we have. 
 * 
 *           Case 1 means we need to subtract out subpaths from the 
 *           sub_cm ROOT_S and ROOT_IL to split states of node 1. This occurs if 
 *           ROOT_IR maps to orig_cm states that come AFTER the orig_cm
 *           states that map to sub_cm ROOT_IL. 
 * 
 *           Case 2 means we need to subtract out subpaths from the 
 *           sub_cm ROOT_S to sub_cm ROOT_IL, and then adjust counts out
 *           of ROOT_IL. This occurs if ROOT_IR maps to orig_cm states 
 *           that come BEFORE the orig_cm states that map to sub_cm 
 *           ROOT_IL (and case 3 does not apply).
 *
 *           Case 3 means we DO NOT need to subtract out any subpaths 
 *           out of ROOT_S or ROOT_IL. This occurs if ROOT_IR maps to 
 *           orig_cm states that come BEFORE the orig_cm states that 
 *           map to sub_cm ROOT_IL AND all the split states of node 1
 *           come BEFORE BOTH the orig_cm states that map to ROOT_IL and
 *           ROOT_IR.
 *
 *           This function is called 4 times, once for each possible
 *           pair of orig_cm states that map to sub_cm ROOT_IL and 
 *           ROOT_IR.
 *
 *
 * Args:    
 * CM_t *orig_cm     - the original, template CM
 * CM_t *sub_cm      - the sub CM
 * double *orig_psi  - for orig_cm: orig_psi[v] is the expected number of times state v is entered
 *                     in a CM parse
 * char ***tmap      - the hard-coded transition map
 * orig2sub_smap     - 2D state map from orig_cm (template) to sub_cm.
 *                     1st dimension - state index in orig_cm 
 *                     2nd D - 2 elements for up to 2 matching sub_cm states, .
 * sub2orig_smap     - 2D state map from sub_cm to orig_cm.
 *                     1st dimension - state index in sub_cm (0..sub_cm->M-1)
 *                     2nd D - 2 elements for up to 2 matching orig_cm states, 
 * int orig_il_or_s  - an orig_cm state that maps to sub_cm ROOT_IL OR 0 for sub_cm state ROOT_S
 * int orig_ir       - an orig_cm state that maps to sub_cm ROOT_IR
 * 
 * Returns: 0, if either orig_il_or_s or orig_ir is -1 
 *          1, if orig_il_or_s < orig_ir (case 1)
 *          2, if orig_il_or_s > orig_ir (case 2)
 *          3, if orig_ss > orig_il_or_s > orig_ir for all orig_ss (split set states of node 1)
 */
static int
cm2sub_cm_subtract_root_subpaths_helper(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, char ***tmap, 
					int **orig2sub_smap, int **sub2orig_smap, int orig_il_or_s, int orig_ir)

{
  int sub_il_or_s;
  int sub_ir;
  int yoffset;
  int sub_y;
  int orig_y;
  int orig_ss;
  int orig_ss1;
  int orig_ss2;
  int skip;

  if(orig_il_or_s == 0)
    {
      sub_il_or_s = 0; /* ROOT_S sub_cm index */
    }
  else
    sub_il_or_s = 1;   /* ROOT_IL sub_cm index */
  sub_ir = 2;          /* ROOT_IR sub_cm index */

  printf("\n\nin cm2sum_cm_subtract_root_subpaths_helper: orig_il_or_s: %d orig_ir: %d\n", orig_il_or_s, orig_ir);
  if(orig_il_or_s == -1 || orig_ir == -1)
    {
      return 0;
    }

  /* check for case 1 or case 2 (should be 1 or the other): */

  /* first check for case 1 */
  /* There could be 1 to 4 instances of each case, 1 if sub_il_or_s and sub_ir both map to exactly 1 orig_cm state,
   * up to 4, if they both map to 2 orig_cm states. We check for each explicitly:                       
   */
  if(orig_il_or_s < orig_ir) /* case 1 if orig_il_or_s == 0 (ROOT_S) this will always be true */
    {
      /* for each split set state of the next node we can reach, we want to subtract the probability
       * that we will go through ROOT_IR to get there, because we've already counted these in the
       * virtual counts out of ROOT_IR.
       */
      for(yoffset = 0; yoffset < sub_cm->cnum[sub_il_or_s]; yoffset++)
	{
	  sub_y = sub_cm->cfirst[sub_il_or_s] + yoffset; 

	  skip = FALSE;
	  if(sub_cm->ndidx[sub_y] != sub_cm->ndidx[sub_il_or_s]) /* if sub_y is in the split set state of the next node */
	    {
	      printf("beg subtract helper loop: sub_cm->t[sub_il_or_s:%d][yoffset:%d]: %f\n", sub_il_or_s, yoffset, sub_cm->t[sub_il_or_s][yoffset]);

	      orig_ss = sub2orig_smap[sub_y][0];
	      if(orig_ss == -1)
		Die("ERROR in cm2sub_cm_subtract_root_subpaths_helper(), case 1 but orig_ss is -1 (orig_ss is orig_cm state that maps to sub_cm state: %d)\n", sub_y);
	      if(orig_ss < orig_il_or_s)
		Die("ERROR in cm2sub_cm_subtract_root_subpaths_helper(), case 1 but orig_ss: %d < orig_il_or_s: %d I thought this was impossible...\n", orig_ss, orig_il_or_s);

	      if(orig_ss > orig_ir)
		{
		  /* We've already correctly handled this situation by skipping the
		   * contribution of orig_ir -> orig_ss when calculating transition
		   * virtual counts out of orig_il_or_s in cm_sum_subpaths().
		   */
		  skip = TRUE;
		}

	      if(!skip)
		sub_cm->t[sub_il_or_s][yoffset] -= orig_psi[orig_il_or_s] * 
		  cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_il_or_s, orig_ss, sub_il_or_s, 
				  tmap, orig_psi, -1, -2, -3, -4) * 
		  cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ss, orig_ir, sub_y, 
				  tmap, orig_psi, -1, -2, -3, -4);
	      
	      skip = FALSE;
	      if(sub2orig_smap[sub_y][1] != -1)
		{
		  orig_ss = sub2orig_smap[sub_y][1];
		  if(orig_ss > orig_ir)
		    {
		      /* We've already correctly handled this situation by skipping the
		       * contribution of orig_ir -> orig_ss when calculating transition
		       * virtual counts out of orig_il_or_s in cm_sum_subpaths().
		       */
		      skip = TRUE;
		    }
		  if(orig_ss < orig_il_or_s)
		    Die("ERROR in cm2sub_cm_subtract_root_subpaths_helper(), case 1 but orig_ss: %d < orig_il_or_s: %d I thought this was impossible...\n", orig_ss, orig_il_or_s);


		  if(!skip)
		    sub_cm->t[sub_il_or_s][yoffset] -= orig_psi[orig_il_or_s] * 
		      cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_il_or_s, orig_ss, sub_il_or_s, 
				      tmap, orig_psi, -1, -2, -3, -4) * 
		      cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ss, orig_ir, sub_y, 
				      tmap, orig_psi, -1, -2, -3, -4);
		}
	      printf("end subtract helper loop: sub_cm->t[sub_il_or_s:%d][yoffset:%d]: %f\n", sub_il_or_s, yoffset, sub_cm->t[sub_il_or_s][yoffset]);
	    }
	}	    
      return 1;
    }

  /* check for case 2, case 3, and case 4 */
  if(orig_il_or_s > orig_ir && orig_ir > sub2orig_smap[sub_cm->nodemap[1]][0])
    {
      /* CASE 3: orig_ss(node1) < orig_ir < orig_il (09.21.06 case 2C) */
      /* Check to make sure that all orig_cm states that map to a split state in node 1 also
       * follow this (I think this is guaranteed, because split states are either matches or deletes,
       * which necessarily map to the same orig_cm NODE, (even in the case of dual mappings MATR_MR to MATP_MP
       * MATP_MR for example)), but we check in case I'm wrong.
       */
      
      /* perform check */
      for(yoffset = 2; yoffset < sub_cm->cnum[0]; yoffset++) /* note we start at yoffset = 2
							      * ROOT_S -> first state of next
							      * node's split set. */
	{
	  sub_y = yoffset + sub_cm->cfirst[0];
	  orig_y = sub2orig_smap[sub_y][0];
	  if(!((orig_il_or_s > orig_y) && (orig_ir > orig_y)))
	    Die("ERROR in cm2sub_cm_subtract_root_subpaths_helper() CASE 3 BUT GUARANTEE VIOLATED 1\n");
	  orig_y = sub2orig_smap[sub_y][1];
	  if(orig_y != -1)
	    if(!((orig_il_or_s > orig_y) && (orig_ir > orig_y)))
	      Die("ERROR in cm2sub_cm_subtract_root_subpaths_helper() CASE 3 BUT GUARANTEE VIOLATED 1\n");
	}

      for(yoffset = 0; yoffset < sub_cm->cnum[0]; yoffset++)
	printf("CASE 3 before t[0][%d] = %f\n", yoffset, sub_cm->t[0][yoffset]);

      /* Adjust the virtual counts for transitions out of ROOT_S. 
       * When compiling virtual counts out of ROOT_S, we have double counted 
       * paths from ROOT_S to all children of ROOT_S that go 
       * through the orig_cm states that map to sub_cm ROOT_IL AND sub_cm ROOT_IR. 
       * These paths have also been counted when compiling counts out of ROOT_IL ROOT_IR.
       */
      for(yoffset = 0; yoffset < sub_cm->cnum[0]; yoffset++)
	{
	  sub_y = sub_cm->cfirst[0] + yoffset; 
	  if(sub_cm->ndidx[sub_y] != 0)
	  {
	    orig_ss1 = sub2orig_smap[sub_y][0];
	    orig_ss2 = sub2orig_smap[sub_y][1];

	    /* subtract paths from root_s -> orig_ss1 -> orig_il */
	    sub_cm->t[0][yoffset] -= 
	      cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, 0, orig_ss1, 
			      0, tmap, orig_psi, -1, -2, -3, -4) * 
	      cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ss1, orig_il_or_s, 
			      0, tmap, orig_psi, -1, -2, -3, -4);

	    printf("\t*!* removed 0 -> ss1: %d -> il: %d sub_cm->t[0][%d]: %f (removed: %f)\n\n", orig_ss1, orig_il_or_s, yoffset, sub_cm->t[0][yoffset], (cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, 0, orig_ss1, 
			      0, tmap, orig_psi, -1, -2, -3, -4) * 
	      cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ss1, orig_il_or_s, 
			      0, tmap, orig_psi, -1, -2, -3, -4)));

	    /* subtract paths from root_s -> orig_ss1 -> orig_ir */
	    sub_cm->t[0][yoffset] -= 
	      cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, 0, orig_ss1, 
			      0, tmap, orig_psi, -1, -2, -3, -4) * 
	      cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ss1, orig_ir, 
			      0, tmap, orig_psi, -1, -2, -3, -4);
	    printf("\t*!* removed 1 -> ss1: %d -> ir: %d sub_cm->t[0][%d]: %f (removed: %f)\n\n", orig_ss1, orig_ir, yoffset, sub_cm->t[0][yoffset], (cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, 0, orig_ss1, 
			      0, tmap, orig_psi, -1, -2, -3, -4) * 										    cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ss1, orig_ir, 
			      0, tmap, orig_psi, -1, -2, -3, -4)));


	    if(orig_ss2 != -1)
	      {
		/* subtract paths from root_s -> orig_ss2 -> orig_il */
		sub_cm->t[0][yoffset] -=
		  cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, 0, orig_ss2, 
				  0, tmap, orig_psi, -1, -2, -3, -4) * 
		  cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ss2, orig_il_or_s, 
				  0, tmap, orig_psi, -1, -2, -3, -4);

		printf("\t*!* removed 3 -> ss2: %d -> il: %d sub_cm->t[0][%d]: %f (removed: %f)\n\n", orig_ss2, orig_il_or_s, yoffset, sub_cm->t[0][yoffset], (		  cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, 0, orig_ss2, 
				  0, tmap, orig_psi, -1, -2, -3, -4) * 
		  cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ss2, orig_il_or_s, 
				  0, tmap, orig_psi, -1, -2, -3, -4)));
 
		/* subtract paths from root_s -> orig_ss2 -> orig_ir */
		sub_cm->t[0][yoffset] -= 
		  cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, 0, orig_ss2, 
				  0, tmap, orig_psi, -1, -2, -3, -4) * 
		  cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ss2, orig_ir, 
				  0, tmap, orig_psi, -1, -2, -3, -4);

		printf("\t*!* removed 4 -> ss2: %d -> ir: %d sub_cm->t[0][%d]: %f (removed: %f)\n\n", orig_ss2, orig_ir, yoffset, sub_cm->t[0][yoffset], (		  cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, 0, orig_ss2, 
				  0, tmap, orig_psi, -1, -2, -3, -4) * 
		  cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ss2, orig_ir, 
				  0, tmap, orig_psi, -1, -2, -3, -4)));

	      }
	  }

	}
      /* Finally, we've over counted paths from orig_cm states ROOT_IR -> ROOT_IL, once when filling virtual counts 
       * from sub_cm states ROOT_IL -> ROOT_IR and once when filling sub_cm states ROOT_S -> ROOT_IR, they
       * belong in the ROOT_IL -> ROOT_IR transition so we subtract them here from the ROOT_S -> ROOT_IR */ 
      sub_cm->t[0][1] -= orig_psi[orig_ir] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ir, orig_il_or_s, sub_ir, tmap, orig_psi, -1, -2, -3, -4);
      
      for(yoffset = 0; yoffset < sub_cm->cnum[0]; yoffset++)
	printf("CASE 3 after t[0][%d] = %f\n", yoffset, sub_cm->t[0][yoffset]);

	printf("!!!RETURNING CASE 3 for orig_il_or_s: %d orig_ir: %d\n", orig_il_or_s, orig_ir);
	return 3;
    }
  else if(orig_il_or_s > sub2orig_smap[sub_cm->nodemap[1]][0] && sub2orig_smap[sub_cm->nodemap[1]][0] > orig_ir)
  {
      /* CASE 4: orig_ir < orig_ss(node1) < orig_il */
      /* Check to make sure that all orig_cm states that map to a split state in node 1 also
       * follow this (I think this is guaranteed, because split states are either matches or deletes,
       * which necessarily map to the same orig_cm NODE, (even in the case of dual mappings MATR_MR to MATP_MP
       * MATP_MR for example)), but we check in case I'm wrong.
       */
      
      /* perform check */
      for(yoffset = 2; yoffset < sub_cm->cnum[0]; yoffset++) /* note we start at yoffset = 2
							      * ROOT_S -> first state of next
							      * node's split set. */
	{
	  sub_y = yoffset + sub_cm->cfirst[0];
	  orig_y = sub2orig_smap[sub_y][0];
	  if(!((orig_il_or_s > orig_y) && (orig_y > orig_ir)))
	    Die("ERROR in cm2sub_cm_subtract_root_subpaths_helper() CASE 4 BUT GUARANTEE VIOLATED 1\n");
	  orig_y = sub2orig_smap[sub_y][1];
	  if(orig_y != -1)
	    if(!((orig_il_or_s > orig_y) && (orig_y > orig_ir)))
	      Die("ERROR in cm2sub_cm_subtract_root_subpaths_helper() CASE 4 BUT GUARANTEE VIOLATED 1\n");
	}

      for(yoffset = 0; yoffset < sub_cm->cnum[0]; yoffset++)
	printf("CASE 4 before t[0][%d] = %f\n", yoffset, sub_cm->t[0][yoffset]);


      /* Adjust the virtual counts for transitions out of ROOT_S. 
       * When compiling virtual counts out of ROOT_S and ROOT_IL, we have double counted 
       * paths from ROOT_S to all children of ROOT_S that go 
       * through the orig_cm states that map to sub_cm ROOT_IL (but not ROOT_IR)
       * These paths have also been counted when compiling counts out of ROOT_IL.
       *
       * We also need to subtract subpaths out of ROOT_IL to split states of
       * node 1 that go through ROOT_IR, as these have also been counted
       * in the transition out of ROOT_IR. This is necessary because
       * orig_ir < orig_il and so cm_sum_subpaths can't ignore the
       * probability that originated at orig_ir.
       */
      for(yoffset = 0; yoffset < sub_cm->cnum[0]; yoffset++)
	{
	  sub_y = sub_cm->cfirst[0] + yoffset; 
	  if(sub_cm->ndidx[sub_y] != 0)
	  {
	    orig_ss1 = sub2orig_smap[sub_y][0];
	    orig_ss2 = sub2orig_smap[sub_y][1];

	    /* First adjust counts out of ROOT_S */
	    /* subtract paths from root_s -> orig_ss1 -> orig_il */
	    sub_cm->t[0][yoffset] -= 
	      cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, 0, orig_ss1, 
			      0, tmap, orig_psi, -1, -2, -3, -4) * 
	      cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ss1, orig_il_or_s, 
			      0, tmap, orig_psi, -1, -2, -3, -4);

	    if(orig_ss2 != -1)
	      {
		/* subtract paths from root_s -> orig_ss2 -> orig_il */
		sub_cm->t[0][yoffset] -=
		  cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, 0, orig_ss2, 
				  0, tmap, orig_psi, -1, -2, -3, -4) * 
		  cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ss2, orig_il_or_s, 
				  0, tmap, orig_psi, -1, -2, -3, -4);
	      }
	    
	    /* Adjust counts out of ROOT_IL */
	    /* subtract paths from orig_ir -> orig_ss1 -> orig_il */
	    sub_cm->t[1][yoffset] -= orig_psi[orig_ir] * 
	      cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ir, orig_ss1, 
			      1, tmap, orig_psi, -1, -2, -3, -4) * 
	      cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ss1, orig_il_or_s, 
			      1, tmap, orig_psi, -1, -2, -3, -4);

	    if(orig_ss2 != -1)
	      {
		/* subtract paths from orig_ir -> orig_ss2 -> orig_il */
		sub_cm->t[1][yoffset] -= orig_psi[orig_ir] * 
		  cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ir, orig_ss2, 
				  1, tmap, orig_psi, -1, -2, -3, -4) * 
		  cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ss2, orig_il_or_s, 
				  1, tmap, orig_psi, -1, -2, -3, -4);
	      }
	  }

	}
      /* Finally, we've over counted paths from orig_cm states ROOT_IR -> ROOT_IL, once when filling virtual counts 
       * from sub_cm states ROOT_IL -> ROOT_IR and once when filling sub_cm states ROOT_S -> ROOT_IR, they
       * belong in the ROOT_IL -> ROOT_IR transition so we subtract them here from the ROOT_S -> ROOT_IR */ 
      sub_cm->t[0][1] -= orig_psi[orig_ir] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ir, orig_il_or_s, sub_ir, tmap, orig_psi, -1, -2, -3, -4);
      
      for(yoffset = 0; yoffset < sub_cm->cnum[0]; yoffset++)
	printf("CASE 4 after t[0][%d] = %f\n", yoffset, sub_cm->t[0][yoffset]);

	printf("!!!RETURNING CASE 4 for orig_il_or_s: %d orig_ir: %d\n", orig_il_or_s, orig_ir);
	return 4;
  }

  else if(orig_il_or_s > orig_ir)
    {
      /* Case 2 (so-named for historical reasons) */
      /* We've over counted paths from orig_cm states ROOT_IR -> ROOT_IL, once when filling virtual counts 
       * from sub_cm states ROOT_IL -> ROOT_IR and once when filling sub_cm states ROOT_S -> ROOT_IR, they
       * belong in the ROOT_IL -> ROOT_IR transition so we subtract them here from the ROOT_S -> ROOT_IR */ 
     sub_cm->t[0][1] -= orig_psi[orig_ir] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ir, orig_il_or_s, sub_ir, tmap, orig_psi, -1, -2, -3, -4);
     return 2;
  }
}


/**************************************************************************
 * EPN 09.21.06
 * Function: cm2sub_cm_subtract_root_subpaths2()
 *
 * Purpose:  When building a sub CM (sub_cm) from an original, template
 *           CM (orig_cm), there's special considerations that must be
 *           taken involving sub_cm nodes with 2 insert states where there
 *           isn't an exactly identical node in the original CM. The only
 *           node that potentially meets this criteria is the ROOT_nd because
 *           all MATP nodes in sub_cm must necessarily also exist in the
 *           orig_cm. 
 *           
 *           The problem is that we have overcounted certain subpaths when
 *           transitioning out of the sub_cm ROOT_S and ROOT_IL states.
 *       
 *           For each sub_cm ROOT_IL state, there are up to 2 orig_cm inserts
 *           that map to it (orig_il1 and orig_il2), and analagously for sub_cm
 *           ROOT_IR (orig_ir1 and orig_ir2). Then there are up to 2 sub_cm states
 *           that map to each of the split set states in sub_cm node 1 (orig_ss1 
 *           and orig_ss2). The subpaths that we have overcounted is dependent
 *           on the relationship between orig_il*, orig_ir*, and orig_ss*.
 *           There are six cases of this relationship:
 *
 *             cases 1A, 1B, 1C apply when il < ir.
 *             case 1A: il < ir < ss (this is correctly handled by cm_sum_subpaths())
 *             case 1B: il < ss < ir
 *             case 1C: ss < il < ir
 *
 *             cases 2A, 2B, 2C apply when ir < il.
 *             case 2A: ir < il < ss
 *             case 2B: ir < ss < il
 *             case 2C: ss < ir < il
 *           
 *           These cases are not explicitly checked for in the code but were
 *           useful for determining the correct strategy to use here, and
 *           are mentioned in comments throughout the code. 
 *
 *           This function calls a helper function for each pair of orig_il* and 
 *           orig_ir* (up to 4 possible pairs), which actually subtracts
 *           the path based on the relationships of orig_il*, orig_ir* and 
 *           orig_ss*.  
 * 
 *           Importantly, though there are potentially many orig_ss* states 
 *           (up to 2 that map to each sub_cm node 1 split set state), 
 *           they all must be close enough in state indices that we'll
 *           never have a case where for a specific il* and ir* we get more
 *           than 1 of the 6 possible cases for all orig_ss*. This is
 *           because either the orig_ss are in a sub_cm MATP node - in
 *           which case they MUST map to exactly 1 state each in the orig_cm
 *           (to the corresponding MATP node) and are contiguous state indices,
 *           or they are a MATL or MATR which either map to a corresponding
 *           MATL, MATR (in which case they're contiguous) or each sub_cm 
 *           split state maps to exactly 2 states in the same orig_cm MATP
 *           node, (for example MATL_ML might map to MATP_MP and MATP_ML)
 *           which are all in the same node and thus have no insert states
 *           with a state index between the indices of the 2 states they
 *           map to (ouch).
 *
 *
 * Args:    
 * CM_t *orig_cm     - the original, template CM
 * CM_t *sub_cm      - the sub CM
 * double *orig_psi  - for orig_cm: orig_psi[v] is the expected number of times state v is entered
 *                     in a CM parse
 * int print_flag    - print useful debugging info
 * Returns: Float, the summed probability of all subpaths through the CM
 *          starting at "start" and ending at "end".
 */
static void
cm2sub_cm_subtract_root_subpaths2(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, char ***tmap, 
				  int **orig2sub_smap, int **sub2orig_smap, int print_flag)

{
  int orig_il1, orig_il2, orig_ir1, orig_ir2;
  int sub_il;
  int sub_ir;

  sub_il = 1; /* sub ROOT_IL index */
  sub_ir = 2; /* sub ROOT_IR index */

  orig_il1 = sub2orig_smap[sub_il][0];
  orig_il2 = sub2orig_smap[sub_il][1];
  
  orig_ir1 = sub2orig_smap[sub_ir][0];
  orig_ir2 = sub2orig_smap[sub_ir][1];

  if(print_flag) 
    printf("in cm2sub_cm_subtract_root_subpaths2 orig_il1: %d | orig_il2: %d | orig_ir1: %d | orig_ir2: %d\n", orig_il1, orig_il2, orig_ir1, orig_ir2);

  cm2sub_cm_subtract_root_subpaths_helper2(orig_cm, sub_cm, orig_psi, tmap, 
					   orig2sub_smap, sub2orig_smap, orig_il1, orig_ir1,
					   1, 1, print_flag);
  cm2sub_cm_subtract_root_subpaths_helper2(orig_cm, sub_cm, orig_psi, tmap, 
					   orig2sub_smap, sub2orig_smap, orig_il1, orig_ir2,
					   2, 1, print_flag);
  cm2sub_cm_subtract_root_subpaths_helper2(orig_cm, sub_cm, orig_psi, tmap, 
					   orig2sub_smap, sub2orig_smap, orig_il2, orig_ir1,
					   1, 2, print_flag);
  cm2sub_cm_subtract_root_subpaths_helper2(orig_cm, sub_cm, orig_psi, tmap, 
					   orig2sub_smap, sub2orig_smap, orig_il2, orig_ir2,
					   2, 2, print_flag);

  /* the last two ints passed to cm2sub_cm_subtract_root_subpaths_helper2() tell the function
   * the how many times that function has been called with the passed in orig_il* and orig_ir*
   * respectively, this is important so we don't over-subtract.
   */
  return;
}

/**************************************************************************
 * EPN 09.21.06
 * Function: cm2sub_cm_subtract_root_subpaths_helper2()
 *
 * Purpose:  Helper function that actually does the subpath subtraction
 *           for the ROOT_nd of sub CMs. See additional notes in 
 *           comments of cm2sub_cm_subtract_root_subpaths2().
 *
 *           This function is called 4 times, once for each possible
 *           pair of orig_cm states that map to sub_cm ROOT_IL and 
 *           ROOT_IR.
 *
 * Args:    
 * CM_t *orig_cm     - the original, template CM
 * CM_t *sub_cm      - the sub CM
 * double *orig_psi  - for orig_cm: orig_psi[v] is the expected number of times state v is entered
 *                     in a CM parse
 * char ***tmap      - the hard-coded transition map
 * orig2sub_smap     - 2D state map from orig_cm (template) to sub_cm.
 *                     1st dimension - state index in orig_cm 
 *                     2nd D - 2 elements for up to 2 matching sub_cm states, .
 * sub2orig_smap     - 2D state map from sub_cm to orig_cm.
 *                     1st dimension - state index in sub_cm (0..sub_cm->M-1)
 *                     2nd D - 2 elements for up to 2 matching orig_cm states, 
 * int orig_il       - an orig_cm state that maps to sub_cm ROOT_IL
 * int orig_ir       - an orig_cm state that maps to sub_cm ROOT_IR
 * int il_iteration  - either 1 or 2, 1 if this is the first call to this function
 *                     using the orig_il passed in, 2 if its the second call with this orig_il
 * int ir_iteration  - either 1 or 2, 1 if this is the first call to this function
 *                     using the orig_ir passed in, 2 if its the second call with this orig_ir
 * int print_flag    - print useful debugging info
 * 
 * Returns: (void)
 */
static void
cm2sub_cm_subtract_root_subpaths_helper2(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, char ***tmap, 
					 int **orig2sub_smap, int **sub2orig_smap, int orig_il, int orig_ir,
					 int il_iteration, int ir_iteration, int print_flag)

{
  int sub_root_s;
  int sub_il;
  int sub_ir;
  int yoffset;
  int sub_y;
  int orig_y;
  int orig_ss;
  int orig_ss1;
  int orig_ss2;

  if(print_flag) printf("\n\nin cm2sum_cm_subtract_root_subpaths_helper: orig_il: %d orig_ir: %d\n", orig_il, orig_ir);
  if(orig_il == -1 || orig_ir == -1)
    {
      if(print_flag) printf("il or ir is -1, returning\n");
      return;
    }

  for(yoffset = 0; yoffset < sub_cm->cnum[0]; yoffset++)
    if(print_flag) printf("Before t[0][%d] = %f\n", yoffset, sub_cm->t[0][yoffset]);

  sub_root_s = 0; /* sub_cm ROOT_S index */
  sub_il     = 1; /* sub_cm ROOT_IL index */
  sub_ir     = 2; /* sub_cm ROOT_IR index */

  orig_ss = sub2orig_smap[3][0]; /* orig_ss is the 1 (of possibly 2) orig_cm states that map to the first
				  * state in sub_cm node 1 (sub_cm state 3)
				  */
  /* Check to make sure that all the split states meet our guarantee (not necessary) */
  for(yoffset = 2; yoffset < sub_cm->cnum[0]; yoffset++) /* note we start at yoffset = 2
							  * ROOT_S -> first state of next
							  * node's split set. */
    {
	  sub_y = yoffset + sub_cm->cfirst[0];
	  orig_y = sub2orig_smap[sub_y][0];
	  if((orig_il > orig_ss && orig_il < orig_y) ||
	     (orig_il < orig_ss && orig_il > orig_y) ||
	     (orig_ir > orig_ss && orig_ir < orig_y) ||
	     (orig_ir < orig_ss && orig_ir > orig_y))
	    Die("ERROR in cm2sub_cm_subtract_root_subpaths_helper2() split set state guarantee violated!\n");
    }
  

  /* Check for which of the 6 cases we have (not actually necessary) */
  if((orig_il < orig_ir) && (orig_ir < orig_ss))
     if(print_flag) printf("ROOT NODE KACHOW case 1A\n");
  if((orig_il < orig_ss) && (orig_ss < orig_ir))
     if(print_flag) printf("ROOT NODE KACHOW case 1B\n");
  if((orig_ss < orig_il) && (orig_il < orig_ir))
     if(print_flag) printf("ROOT NODE KACHOW case 1C\n");

  if((orig_ir < orig_il) && (orig_il < orig_ss))
     if(print_flag) printf("ROOT NODE KACHOW case 2A\n");
  if((orig_ir < orig_ss) && (orig_ss < orig_il))
     if(print_flag) printf("ROOT NODE KACHOW case 2B\n");
  if((orig_ss < orig_ir) && (orig_ir < orig_il))
     if(print_flag) printf("ROOT NODE KACHOW case 2C\n");


  /* First adjust counts out of sub_cm ROOT_S */

  /* if orig_ir < orig_il, we've counted paths from
   * ir -> il correctly for ROOT_IL -> ROOT_IR and
   *        incorrectly for ROOT_S  -> ROOT_IR
   */
  if(orig_ir < orig_il)
    {
      if(print_flag) printf("0 sub from S->IR, IR -> IL\n"); /* cases 2A, 2B, 2C */
      sub_cm->t[sub_root_s][1] -= orig_psi[orig_ir] * 
	cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ir, orig_il, 0, tmap, orig_psi, -1, -2, -3, -4);
    }
  
  /* All the remaining paths to be subtracted involve the split set state of 
   * node 1. */
  for(yoffset = 0; yoffset < sub_cm->cnum[0]; yoffset++)
    {
      sub_y = sub_cm->cfirst[0] + yoffset; 
      orig_ss1 = sub2orig_smap[sub_y][0];
      orig_ss2 = sub2orig_smap[sub_y][1];
      if(sub_cm->ndidx[sub_y] != 0)
	{
	  /* if orig_ss < orig_il, we've counted paths from
	   * root_s -> ss -> il correctly for ROOT_IL -> ROOT_SS and
	   *                  incorrectly for ROOT_S  -> ROOT_SS
	   *
	   * if il_iteration != 1, then we've already subtracted this
	   * quantity in a previous call.
	   */
	  if(il_iteration == 1 && orig_ss < orig_il) /* cases 1C, 2B, 2C */
	    {
	      if(print_flag) printf("1 (1c, 2b, 2c)  sub from S, S -> SS -> IL\n");
	      /* subtract paths from root_s -> orig_ss1 -> orig_il */
	      sub_cm->t[sub_root_s][yoffset] -= 
		cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, 0, orig_ss1, 
				0, tmap, orig_psi, -1, -2, -3, -4) * 
		cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ss1, orig_il, 
				0, tmap, orig_psi, -1, -2, -3, -4);
	      
	      if(orig_ss2 != -1)
		/* subtract paths from root_s -> orig_ss2 -> orig_il */
		sub_cm->t[sub_root_s][yoffset] -= 
		  cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, 0, orig_ss2, 
				  0, tmap, orig_psi, -1, -2, -3, -4) * 
		  cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ss2, orig_il, 
				  0, tmap, orig_psi, -1, -2, -3, -4);
	    }

	  /* if orig_ss < orig_ir, we've counted paths from
	   * root_s -> ss -> ir correctly for ROOT_IR -> ROOT_SS and
	   *                  incorrectly for ROOT_S  -> ROOT_SS
	   *
	   * if ir_iteration != 1, then we've already subtracted this
	   * quantity in a previous call.
	   */
	  if(ir_iteration == 1 && orig_ss < orig_ir) /* cases 1B, 1C, 2C */
	    {
	      if(print_flag) printf("2 (1b, 1c, 2c) sub from S, S -> SS -> IR\n");
	      /* subtract paths from root_s -> orig_ss1 -> orig_ir */
	      sub_cm->t[sub_root_s][yoffset] -= 
		cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, 0, orig_ss1, 
				0, tmap, orig_psi, -1, -2, -3, -4) * 
		cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ss1, orig_ir, 
				0, tmap, orig_psi, -1, -2, -3, -4);
	      
	      if(orig_ss2 != -1)
	      /* subtract paths from root_s -> orig_ss2 -> orig_ir */
		sub_cm->t[sub_root_s][yoffset] -= 
		  cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, 0, orig_ss2, 
				  0, tmap, orig_psi, -1, -2, -3, -4) * 
		  cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ss2, orig_ir, 
				  0, tmap, orig_psi, -1, -2, -3, -4);
	      
	    }

	  /* Adjust counts out of ROOT_IL if necessary */

	  /* if orig_ir < orig_ss < orig_il, we've counted paths from
	   * ir -> ss -> il correctly for ROOT_IL -> ROOT_IR and
	   *              incorrectly for ROOT_IL -> ROOT_SS
	   */
	  if(orig_ir < orig_ss && orig_ss < orig_il) /* case 2B only */
	    {
	      if(print_flag) printf("3 (2b) sub from IL, IR -> SS -> IL\n");

	      /* subtract paths from orig_ir -> orig_ss1 -> orig_il */
	      sub_cm->t[sub_il][yoffset] -= orig_psi[orig_ir] * 
		cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ir, orig_ss1, 
				sub_il, tmap, orig_psi, -1, -2, -3, -4) * 
		cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ss1, orig_il, 
				sub_il, tmap, orig_psi, -1, -2, -3, -4);

	      if(orig_ss2 != -1)
		/* subtract paths from orig_ir -> orig_ss2 -> orig_il */
		sub_cm->t[1][yoffset] -= orig_psi[orig_ir] * 
		  cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ir, orig_ss2, 
				  sub_il, tmap, orig_psi, -1, -2, -3, -4) * 
		  cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ss2, orig_il, 
				  sub_il, tmap, orig_psi, -1, -2, -3, -4);
	      }

	  /* if orig_ss < orig_il < orig_ir, we've counted paths from
	   * ss -> il -> ir correctly for ROOT_IR -> ROOT_SS and
	   *              incorrectly for ROOT_IL -> ROOT_SS
	   */
	  if(orig_ss < orig_il && orig_il < orig_ir) /* case 1C only */
	    {
	      if(print_flag) printf("4 (1c) sub from IL, SS -> IL -> IR\n");

	      /* subtract paths from orig_ss1 -> orig_il (add self insert) -> orig_ir */
	      sub_cm->t[sub_il][yoffset] -= orig_psi[orig_ss1] * 
		cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ss1, orig_il, 
				sub_il, tmap, orig_psi, -1, -2, -3, -4) *
		(1. + (orig_cm->t[orig_il][0] / (1 - orig_cm->t[orig_il][0]))) *
		cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_il, orig_ir, 
				sub_il, tmap, orig_psi, -1, -2, -3, -4);

	      if(orig_ss2 != -1)
		/* subtract paths from orig_ss2 -> orig_il (add self insert) -> orig_ir */
		sub_cm->t[sub_il][yoffset] -= orig_psi[orig_ss2] * 
		  cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ss2, orig_il, 
				  sub_il, tmap, orig_psi, -1, -2, -3, -4) *
		  (1. + (orig_cm->t[orig_il][0] / (1 - orig_cm->t[orig_il][0]))) *
		  cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_il, orig_ir, 
				  sub_il, tmap, orig_psi, -1, -2, -3, -4);
	    }

	  /* if orig_ir < orig_il < orig_ss, we've counted paths from
	   * ir -> il -> ss correctly for ROOT_IR -> ROOT_SS and
	   *              incorrectly for ROOT_IL -> ROOT_SS
	   */
	  if(orig_ir < orig_il && orig_il < orig_ss) /* case 2A only */
	    {
	      if(print_flag) printf("5 (2a) sub from IL, IR -> IL -> SS\n");

	      /* subtract paths from orig_ir -> orig_il (add self insert) -> orig_ss1 */
	      sub_cm->t[sub_il][yoffset] -= orig_psi[orig_ir] * 
		cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ir, orig_il, 
				sub_il, tmap, orig_psi, -1, -2, -3, -4) *
		(1. + (orig_cm->t[orig_il][0] / (1 - orig_cm->t[orig_il][0]))) *
		cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_il, orig_ss1, 
				sub_il, tmap, orig_psi, -1, -2, -3, -4);

	      if(orig_ss2 != -1)
		/* subtract paths from orig_ir -> orig_il (add self insert) -> orig_ss2 */
	      sub_cm->t[sub_il][yoffset] -= orig_psi[orig_ir] * 
		cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ir, orig_il, 
				sub_il, tmap, orig_psi, -1, -2, -3, -4) *
		(1. + (orig_cm->t[orig_il][0] / (1 - orig_cm->t[orig_il][0]))) *
		cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_il, orig_ss2, 
				sub_il, tmap, orig_psi, -1, -2, -3, -4);

	    }

	  /* if orig_il < orig_ss < orig_ir, we've counted paths from
	   * il -> ss -> ir correctly for ROOT_IL -> ROOT_IR and
	   *              incorrectly for ROOT_IL -> ROOT_SS
	   */
	  if(orig_il < orig_ss && orig_ss < orig_ir) /* case 1B only */
	    {
	      if(print_flag) printf("6 (1b) sub from IL, IL -> SS -> IR\n");

	      if(print_flag) printf("1B before sub 1: sub_cm->t[sub_il:%d][yoffset:%d]: %f\n", sub_il, yoffset, sub_cm->t[sub_il][yoffset]);
	      /* subtract paths from orig_il -> orig_ss1 -> orig_ir */
	      sub_cm->t[sub_il][yoffset] -= orig_psi[orig_il] * 
		cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_il, orig_ss1, 
				sub_il, tmap, orig_psi, -1, -2, -3, -4) *
		cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ss1, orig_ir, 
				sub_il, tmap, orig_psi, -1, -2, -3, -4);

	      if(print_flag) printf("1B before sub 2: sub_cm->t[sub_il:%d][yoffset:%d]: %f\n", sub_il, yoffset, sub_cm->t[sub_il][yoffset]);
	      if(orig_ss2 != -1)
		/* subtract paths from orig_il -> orig_ss2 -> orig_ir */
		sub_cm->t[sub_il][yoffset] -= orig_psi[orig_il] * 
		  cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_il, orig_ss2, 
				  sub_il, tmap, orig_psi, -1, -2, -3, -4) *
		  cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ss2, orig_ir, 
				  sub_il, tmap, orig_psi, -1, -2, -3, -4);
	      if(print_flag) printf("1B after sub 2: sub_cm->t[sub_il:%d][yoffset:%d]: %f\n", sub_il, yoffset, sub_cm->t[sub_il][yoffset]);

	    }
	}
    }
  for(yoffset = 0; yoffset < sub_cm->cnum[0]; yoffset++)
    if(print_flag) printf("After t[0][%d] = %f\n", yoffset, sub_cm->t[0][yoffset]);
  
  return;
}

/**************************************************************************
 * EPN 09.22.06
 * Function: map_orig2sub_cm3()
 * 
 * Purpose:  Determine maps between a template CM and a sub CM, which probably
 *           has less structure (less MATP nodes) than the template. Do this
 *           a bit indirectly, first get maps from each CM to a CP9 HMM,
 *           then use these maps to map the CMs to each other.
 * Args:    
 * CM_t *orig_cm       - the original template CM
 * CM_t *sub_cm        - the sub CM
 * int ***ret_orig2sub_smap - 2D state map from orig_cm (template) to sub_cm.
 *                       1st dimension - state index in orig_cm 
 *                       2nd D - 2 elements for up to 2 matching sub_cm states, .
 *                       Allocated here, responsibility of caller to free.
 * int ***ret_sub2orig_smap - 2D state map from orig_cm (template) to sub_cm.
 *                       1st dimension - state index in sub_cm (0..sub_cm->M-1)
 *                       2nd D - 2 elements for up to 2 matching orig_cm states, 
 *                       Allocated here, responsibility of caller to free.
 * int sub_start       - first consensus column in orig_cm that sub_cm models
 * int sub_end         -  last consensus column in orig_cm that sub_cm models
 * int **ret_sub2orig_id - [0..sub_cm->M], 1 if sub2orig_id_orig[sub_v] == 1
 *                       IFF there is a state in the template, original CM
 *                       that is identical to sub_v, i.e. will have the
 *                       exact same emission and transition distribution.
 * int  *imp_cc        - [0..(sub_end-sub_start+1)] ret_imp_cc[k] = 1 if it is 
 *                       impossible for CP9 node (consensus column) k to be
 *                       calculated in the sub_cm to have distros to match the
 *                       corresponding CP9 node in the original CM - due to
 *                       topological differences in the architecture of the
 *                       sub_cm and orig_cm.
 * int print_flag      - TRUE to print out useful debugging info
 * Returns: (void) 
 */
void
map_orig2sub_cm3(CM_t *orig_cm, CM_t *sub_cm, int ***ret_orig2sub_smap, int ***ret_sub2orig_smap,
		 int sub_start, int sub_end, int **ret_sub2orig_id, int *imp_cc, int print_flag)
{
  int **cs2hn_map_o;
  int **cs2hs_map_o;
  int ***hns2cs_map_o;
  int **cs2hn_map_s;
  int **cs2hs_map_s;
  int ***hns2cs_map_s;
  struct cplan9_s  *cp9hmm_s;        /* constructed CM p9 HMM from the sub_cm */
  struct cplan9_s  *cp9hmm_o;        /* constructed CM p9 HMM from the orig_cm */

  int         *node_cc_left_o; /* consensus column each CM node's left emission maps to
	                      * [0..(cm->nodes-1)], -1 if maps to n_o consensus column*/
  int         *node_cc_right_o;/* consensus column each CM node's right emission maps to
			            * [0..(cm->nodes-1)], -1 if maps to n_o consensus column*/
  int         *cc_node_map_o;  /* node that each consensus column maps to (is modelled by)
			            * [1..hmm_nmc] */
  int         *node_cc_left_s; /* consensus column each CM node's left emission maps to
	                      * [0..(cm->nodes-1)], -1 if maps to n_o consensus column*/
  int          *node_cc_right_s;/* consensus column each CM node's right emission maps to
			            * [0..(cm->nodes-1)], -1 if maps to n_o consensus column*/
  int         *cc_node_map_s;  /* sub_cm node that each consensus column maps to (is modelled by)
			            * [1..hmm_nmc] */
  int k_s;  /* HMM node counter */
  int k_o;  /* HMM node counter */
  int ks_s; /* HMM state counter (0(Match) 1(insert) or 2(delete)*/
  int ks_o; /* HMM state counter (0(Match) 1(insert) or 2(delete)*/
  int v_o;
  int v_s;
  int v; /* state index in CM */
  CMEmitMap_t *orig_emap;         /* consensus emit map for the original, template CM */
  CMEmitMap_t *sub_emap;         /* consensus emit map for the sub_cm */
  int cc;
  int n_o;
  int n_s;
  char **nodetypes;
  char **sttypes;
  int *orig_cc2lins_map;
  int *orig_cc2rins_map;
  int **orig2sub_smap;
  int **sub2orig_smap;
  int *sub2orig_id;
  int cp9_M;
  int sub_k;
  int orig_k;
  int sub_nd;
  int orig_nd;
  int is_right;
  int id_stretch; 

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

  nodetypes = malloc(sizeof(char *) * 8);
  nodetypes[0] = "BIF";
  nodetypes[1] = "MATP";
  nodetypes[2] = "MATL";
  nodetypes[3] = "MATR";
  nodetypes[4] = "BEGL";
  nodetypes[5] = "BEGR";
  nodetypes[6] = "ROOT";
  nodetypes[7] = "END";

  /* sanity check */
  if(sub_cm->M > orig_cm->M)
    Die("ERROR: sub_cm has more states than orig_cm in map_orig2sub_cm()\n");

  /* Get emitmap's for each CM */
  orig_emap = CreateEmitMap(orig_cm);
  sub_emap = CreateEmitMap(sub_cm);

  /* To make things easy, and avoid going traversing the emap multiple times
   * wastefully, we make maps of which nodes have an insert state that
   * inserts AFTER (in case of *lmap) and BEFORE (in case of *rmap)
   * each consensus node. 
   * orig_cc2lins_map[cc] = v, where orig_cm state v is an IL_st that
   *                           emits AFTER consensus column cc.
   * orig_cc2rins_map[cc] = v, where orig_cm state v is an IR_st that
   *                           emits BEFORE consensus column cc.
   * if no such state exists the value will be -1.
   */

  /* Allocate and initialize */
  sub2orig_id = MallocOrDie(sizeof(int) * (sub_cm->M + 1));
  for(v = 0; v <= sub_cm->M; v++)
    sub2orig_id[v] = FALSE;

  orig_cc2lins_map = MallocOrDie(sizeof(int) * (orig_emap->clen + 2));
  orig_cc2rins_map = MallocOrDie(sizeof(int) * (orig_emap->clen + 2));
  for(cc = 0; cc <= orig_emap->clen+1; cc++)
    {
      orig_cc2lins_map[cc] = -1;
      orig_cc2rins_map[cc] = -1;
    }
  /* fill in the map */
  /* ROOT is special */
  orig_cc2lins_map[0] = 1; /* ROOT_IL */
  orig_cc2rins_map[orig_emap->clen+1] = 2; /* ROOT_IR */
  for(n_o = 0; n_o < orig_cm->nodes; n_o++)
    {
      switch (orig_cm->ndtype[n_o]) {
      case MATP_nd:
	orig_cc2lins_map[orig_emap->lpos[n_o]] = orig_cm->nodemap[n_o] + 4; /* MATP_IL */
	orig_cc2rins_map[orig_emap->rpos[n_o]] = orig_cm->nodemap[n_o] + 5; /* MATP_IR */
	break;
	
      case MATL_nd:
	orig_cc2lins_map[orig_emap->lpos[n_o]] = orig_cm->nodemap[n_o] + 2; /* MATL_IL */
	break;
	
      case MATR_nd:
	orig_cc2rins_map[orig_emap->rpos[n_o]] = orig_cm->nodemap[n_o] + 2; /* MATR_IR */
	break;
	
      case BEGR_nd:
	orig_cc2lins_map[orig_emap->lpos[n_o]] = orig_cm->nodemap[n_o] + 1; /* BEGR_IL */
	break;
	
      default: {} /*do nothing*/
      }
    }

  if(print_flag)
    {
      for(cc = 0; cc <= orig_emap->clen+1; cc++)
	{
	  printf("!cc: %d lins: %d rins: %d\n", cc, orig_cc2lins_map[cc], orig_cc2rins_map[cc]);
	}
    }

  cp9_M = 0;
  for(v = 0; v <= orig_cm->M; v++)
    {
      if(orig_cm->stid[v] ==  MATP_MP)
	cp9_M += 2;
      else if(orig_cm->stid[v] == MATL_ML || orig_cm->stid[v] == MATR_MR)
	cp9_M++;
    }
  cp9hmm_o = AllocCPlan9(cp9_M);
  ZeroCPlan9(cp9hmm_o);
  CPlan9SetNullModel(cp9hmm_o, orig_cm->null, 1.0); /* set p1 = 1.0 which corresponds to the CM */

  cp9hmm_s = AllocCPlan9((sub_end-sub_start+1));
  ZeroCPlan9(cp9hmm_s);
  CPlan9SetNullModel(cp9hmm_s, sub_cm->null, 1.0); /* set p1 = 1.0 which corresponds to the CM */

  /* map the nodes of each CM to consensus column indices and vice versa */
  map_consensus_columns(orig_cm, cp9hmm_o->M, &node_cc_left_o, &node_cc_right_o,
			&cc_node_map_o, 0);
  map_consensus_columns(sub_cm,  cp9hmm_s->M, &node_cc_left_s, &node_cc_right_s,
			&cc_node_map_s, 0);
  /* map the CM states to CP9 states and nodes and vice versa */
  CP9_map_cm2hmm_and_hmm2cm(orig_cm, cp9hmm_o, node_cc_left_o, node_cc_right_o, 
			    cc_node_map_o, &cs2hn_map_o, &cs2hs_map_o, 
			    &hns2cs_map_o, 0);
  CP9_map_cm2hmm_and_hmm2cm(sub_cm, cp9hmm_s, node_cc_left_s, node_cc_right_s, 
			    cc_node_map_s, &cs2hn_map_s, &cs2hs_map_s, 
			    &hns2cs_map_s, 0);

  if(print_flag)
    {
      printf("PRINTING CC NODE MAP ORIG\n");
      for(cc = 1; cc <= orig_emap->clen; cc++)
	printf("cc_node_map_o[%d]: %d\n", cc, cc_node_map_o[cc]);
      printf("DONE PRINTING CC NODE MAP ORIG\n");
      
      printf("PRINTING CC NODE MAP SUB\n");
      for(cc = 1; cc <= sub_emap->clen; cc++)
	printf("cc_node_map_s[%d]: %d\n", cc, cc_node_map_s[cc]);
      printf("DONE PRINTING CC NODE MAP SUB\n");
    }

  /* Allocate the maps */
  orig2sub_smap = MallocOrDie(sizeof(int *) * (orig_cm->M));
  for(v_o = 0; v_o < orig_cm->M; v_o++)
    orig2sub_smap[v_o] = MallocOrDie(sizeof(int) * 2);

  sub2orig_smap = MallocOrDie(sizeof(int *) * (sub_cm->M));
  for(v_s = 0; v_s < sub_cm->M; v_s++)
    sub2orig_smap[v_s] = MallocOrDie(sizeof(int) * 2);
    
  /* Initialize the maps */
  for(v_o = 0; v_o < orig_cm->M; v_o++)
    {
      orig2sub_smap[v_o][0] = -1;
      orig2sub_smap[v_o][1] = -1;
    }
  for(v_s = 0; v_s < sub_cm->M; v_s++)
    {
      sub2orig_smap[v_s][0] = -1;
      sub2orig_smap[v_s][1] = -1;
    }

  /* step through the states of sub_cm, filling in maps */
  /* ROOT is special: */
  v_s = 0; /* ROOT_S */
  v_o = 0; /* ROOT_S */
  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);

  /* ROOT_IL inserts before orig cc sub_start */
  k_s = 1; /* insert */
  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
			 hns2cs_map_o[sub_start-1][k_s][0],
			 hns2cs_map_s[0][k_s][0]);     
  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
			 hns2cs_map_o[sub_start-1][k_s][0],
			 hns2cs_map_s[0][k_s][1]);     
  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
			 hns2cs_map_o[sub_start-1][k_s][1],
			 hns2cs_map_s[0][k_s][0]);     
  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
			 hns2cs_map_o[sub_start-1][k_s][1],
			 hns2cs_map_s[0][k_s][1]);     

  /* ROOT_IR inserts after orig cc sub_end */
  k_s = 1; /* insert */
  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
			 hns2cs_map_o[sub_end][k_s][0],
			 hns2cs_map_s[sub_end-sub_start+1][k_s][0]);     
  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
			 hns2cs_map_o[sub_end][k_s][0],
			 hns2cs_map_s[sub_end-sub_start+1][k_s][1]);     
  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
			 hns2cs_map_o[sub_end][k_s][1],
			 hns2cs_map_s[sub_end-sub_start+1][k_s][0]);     
  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
			 hns2cs_map_o[sub_end][k_s][1],
			 hns2cs_map_s[sub_end-sub_start+1][k_s][1]);     

  /* BEGR_ILs */

  for(n_s = 0; n_s < sub_cm->nodes; n_s++)
    {
      if(sub_cm->ndtype[n_s] == BEGR_nd)
	{
	  v_s = sub_cm->nodemap[n_s] + 1;
	  v_o = orig_cc2lins_map[(sub_emap->lpos[n_s] + sub_start - 1)];
	  if(v_o != -1)
	    map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);

	  v_o = orig_cc2rins_map[(sub_emap->lpos[n_s] + sub_start)];
	  if(v_o != -1)
	    map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, v_o, v_s);
	}
    }

  id_stretch = FALSE;
  for(sub_k = 1; sub_k <= cp9hmm_s->M; sub_k++)
    {
      orig_k = sub_k + sub_start - 1;
      
      sub_nd  = cc_node_map_s[sub_k];
      orig_nd = cc_node_map_o[orig_k];

      /* First check for the easy case: when same node type
       * maps to the same consensus column. If this is 
       * true than we can set sub2orig_id[sub_v] = TRUE
       * for all states in this node.
       */
      if(orig_cm->ndtype[orig_nd] == sub_cm->ndtype[sub_nd])
	{
	  cm2sub_cm_check_id_next_node(orig_cm, sub_cm, orig_nd, sub_nd, sub2orig_id, node_cc_left_o,
				       node_cc_right_o, node_cc_left_s, node_cc_right_s, sub_start, 
				       print_flag);

	  if(node_cc_left_o[orig_nd] == orig_k)
	    is_right = FALSE;
	  else if(node_cc_right_o[orig_nd] == orig_k)
	    is_right = TRUE;
	  else
	    Die("You don't understand CMs 0\n");

	  //printf("EASY CASE\n");
	  /* map all corresponding states */
	  v_o = orig_cm->nodemap[orig_nd];
	  v_s =  sub_cm->nodemap[sub_nd];
	  while(orig_cm->ndidx[v_o] == orig_nd)
	    {
	      if(sub_cm->ndidx[v_s] != sub_nd)
		Die("ERROR 2 in map_orig2sub_cm()\n");
	      map_orig2sub_cm_helper(orig_cm, sub_cm,
				     orig2sub_smap, sub2orig_smap, v_o, v_s);
	      /* If we're mapping insert states in NON-MATP nodes, 
	       * we may have another orig_cm insert that should map to this state */
	      if((sub_cm->ndtype[sub_cm->ndidx[v_s]] != MATP_nd) &&
		 (sub_cm->sttype[v_s] == IL_st || sub_cm->sttype[v_s] == IR_st))
		{
		  k_s = 1; /* insert */		 
		  //printf("!*! sub_k: %d | orig_k: %d\n", sub_k, orig_k);
		  if(hns2cs_map_o[orig_k][k_s][0] != -1 && hns2cs_map_o[(orig_k-is_right)][k_s][1] != v_o)
		    {
		      if(map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
						hns2cs_map_o[orig_k-is_right][k_s][0], v_s))
			sub2orig_id[v_s] = FALSE; /* this state maps to 2 guys, its not an id case */
		    }
		  if(hns2cs_map_o[orig_k][k_s][1] != -1 && hns2cs_map_o[(orig_k-is_right)][k_s][1] != v_o)
		    {
		      if(map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
						hns2cs_map_o[orig_k-is_right][k_s][1], v_s))
			sub2orig_id[v_s] = FALSE; /* this state maps to 2 guys, its not an id case */
		    }
		}
	      v_o++; 
	      v_s++;
	      id_stretch = TRUE;
	    }
	}
      /* end easy case check */

      else
	{
	  if(id_stretch && (is_right && sub_cm->ndtype[sub_nd] == MATR_nd))
	    {
	      /* previous column was an 'id case' where we had an identical node
	       * in the sub_cm to the orig_cm, in this case we need to check for
	       * a special situation here: when prev column was modelled by a right
	       * half of an orig_cm node, and the current sub_cm node is a MATR,
	       * in this case we may have missed a mapping from the MATR_IR of this
	       * state to the orig_cm state that inserts AFTER the prev_cc.
	       */
	      //printf("SPECIAL EASY CASE STRETCH NOW FALSE, MATR CASE\n");
	      k_s = 1; /*insert*/
	      map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				     hns2cs_map_o[orig_k-1][k_s][0],
				     hns2cs_map_s[sub_k-1][k_s][0]);     
	      map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				     hns2cs_map_o[orig_k-1][k_s][0],
				     hns2cs_map_s[sub_k-1][k_s][1]);     
	      map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				     hns2cs_map_o[orig_k-1][k_s][1],
				     hns2cs_map_s[sub_k-1][k_s][0]);     
	      map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				     hns2cs_map_o[orig_k-1][k_s][1],
				     hns2cs_map_s[sub_k-1][k_s][1]);     
	    }

	  id_stretch = FALSE;

	  //printf("sub_k: %d | orig_k: %d\n", sub_k, orig_k);
	  k_s = 0; /* match */
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				 hns2cs_map_o[orig_k][k_s][0],
				 hns2cs_map_s[sub_k][k_s][0]);     
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				 hns2cs_map_o[orig_k][k_s][0],
				 hns2cs_map_s[sub_k][k_s][1]);     
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				 hns2cs_map_o[orig_k][k_s][1],
				 hns2cs_map_s[sub_k][k_s][0]);     
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				 hns2cs_map_o[orig_k][k_s][1],
				 hns2cs_map_s[sub_k][k_s][1]);     
	  
	  k_s = 1; /* insert */
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				 hns2cs_map_o[orig_k][k_s][0],
				 hns2cs_map_s[sub_k][k_s][0]);     
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				 hns2cs_map_o[orig_k][k_s][0],
				 hns2cs_map_s[sub_k][k_s][1]);     
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				 hns2cs_map_o[orig_k][k_s][1],
				 hns2cs_map_s[sub_k][k_s][0]);     
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				 hns2cs_map_o[orig_k][k_s][1],
				 hns2cs_map_s[sub_k][k_s][1]);     
	  
	  
	  k_s = 2; /* delete */
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				 hns2cs_map_o[orig_k][k_s][0],
				 hns2cs_map_s[sub_k][k_s][0]);     
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				 hns2cs_map_o[orig_k][k_s][0],
				 hns2cs_map_s[sub_k][k_s][1]);     
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				 hns2cs_map_o[orig_k][k_s][1],
				 hns2cs_map_s[sub_k][k_s][0]);     
	  map_orig2sub_cm_helper(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				 hns2cs_map_o[orig_k][k_s][1],
				 hns2cs_map_s[sub_k][k_s][1]);    
	}
    }
  /* NOTE: We ignore mapping B states, E states and S states (except ROOT_S). We'll handle these guys 
   * specially when we fill in transition probabilities. The reason we don't map them is that
   * I think its impossible to do it robustly. I saw cases with SSU where there were BIF nodes in the
   * sub_cm for which I couldn't figure out which BIF nodes (and correspondingly BEGL and BEGR nodes)
   * in the orig_cm they should map to. Regardless, even if there is a pretty, nice way I'm abandoning
   * an attempt to find it for a crude way - we will handle the transitions involving these guys
   * specially, without a need for a mapping between CMs.
   */
  
  /* Check to make sure sub2orig_id makes sense */
  for(v = 0; v <= sub_cm->M; v++)
    {
      if(sub2orig_id[v] == TRUE)
	{
	  if(print_flag) printf("sub2orig_id[%d] TRUE\n", v);
	  if(sub2orig_smap[v][1] != -1)
	    {
	      /*Die("ERROR, sub2orig_id[v:%d] but sub2orig_smap[v:%d][0]: %d and sub2orig_smap[v:%d][1]: %d\n", v, v, sub2orig_smap[v][0], v, sub2orig_smap[v][1]);*/
	    }
	}
      else
	{
	  if(print_flag) printf("sub2orig_id[%d] FALSE\n", v);
	}
    }
  
  /* print sub2orig_smap, checking consistency with orig2sub_smap along the way. */
  if(print_flag) printf("\n\n\n KACHOW! MAP\n\n\n");
  for(v_s = 0; v_s < sub_cm->M; v_s++)
    {
      n_s = sub_cm->ndidx[v_s];
      v_o = sub2orig_smap[v_s][0];
      /*      printf("v_s: %d\n", v_s);
	      printf("n_s: %d\n", n_s);
	      printf("v_o: %d\n", v_o);*/
      if(print_flag)
	{
	  if(sub_cm->sttype[v_s] == E_st) 	
	    {
	      printf("sub v:%4d   END\n", v_s);
	    }
	  if(sub_cm->sttype[v_s] == S_st) 	
	    {
	      printf("sub v:%4d   START\n", v_s);
	    }
	  if(sub_cm->sttype[v_s] == B_st) 	
	    {
	      printf("sub v:%4d   BIF_B\n", v_s);
	    }
	}

      if((sub_cm->sttype[v_s] != B_st && 
	  sub_cm->sttype[v_s] != S_st) && 
	 sub_cm->sttype[v_s] != E_st)
	{
	  if(v_o == -1 && sub_cm->sttype[(v_s+1)] == E_st) /* v_s is a dead insert */
	    continue;
	  if(v_o == -1 && sub_cm->sttype[v_s] != E_st)
	    Die("ERROR sub_cm state: %d type: %s node type: %s doesn't map to any state in orig_cm\n", v_s, sttypes[sub_cm->sttype[v_s]], nodetypes[sub_cm->ndtype[n_s]]);

	  n_o = orig_cm->ndidx[v_o];
	  if(print_flag) printf("sub v:%4d(%4d) %6s%6s | orig v:%4d(%4d) %6s%6s\n", v_s, n_s, nodetypes[sub_cm->ndtype[n_s]], sttypes[sub_cm->sttype[v_s]], v_o, n_o, nodetypes[orig_cm->ndtype[n_o]], sttypes[orig_cm->sttype[v_o]]);
	  /* check to make sure orig2sub_smap is consistent */
	  if(orig2sub_smap[v_o][0] != v_s && orig2sub_smap[v_o][1] != v_s)
	    {
	      Die("ERROR inconsistency; neither orig2sub_smap[%d][0] and [1] is %d\n", v_o, v_s);
	    }

	  v_o = sub2orig_smap[v_s][1];
	  if(v_o != -1)
	    {
	      n_o = orig_cm->ndidx[v_o];
	      if(print_flag) printf("                              | orig v:%4d(%4d) %6s%6s\n", v_o, n_o, nodetypes[orig_cm->ndtype[n_o]], sttypes[orig_cm->sttype[v_o]]);
	      /* check to make sure orig2sub_smap is consistent */
	      if(orig2sub_smap[v_o][0] != v_s && orig2sub_smap[v_o][1] != v_s)
		{
		  Die("ERROR inconsistency; neither orig2sub_smap[%d][0] and [1] is %d\n", v_o, v_s);
		}
	    }
	}
    }

  /* 10.05.06
   * Look for 'impossible' cases where we know the sub_cm construction procedure fails, in that
   * the distribution of transitions out of CP9 nodes built from the sub_cm will be the same distros
   * out of corresponding CP9 nodes built from the full CM. */
  cm2sub_cm_find_impossible_misc_cases(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, 
				       cc_node_map_o, node_cc_left_o, node_cc_right_o, 
				       cc_node_map_s, node_cc_left_s, node_cc_right_s, imp_cc, 
				       sub_start, sub_end, print_flag);
  cm2sub_cm_find_impossible_matr_cases(orig_cm, sub_cm, cc_node_map_o, node_cc_left_o, node_cc_right_o,
				       cc_node_map_s, node_cc_left_s, node_cc_right_s, imp_cc, 
				       sub_start, sub_end, print_flag);
  /*debug_print_misc_sub_cm_info(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, sub_start, sub_end, cc_node_map_o, 
    node_cc_left_o, node_cc_right_o);*/

  /* Clean up and return */
  free(node_cc_left_o);
  free(node_cc_right_o);
  free(cc_node_map_o);
  free(node_cc_left_s);
  free(node_cc_right_s);
  free(cc_node_map_s);
  free(orig_cc2lins_map);
  free(orig_cc2rins_map);
  FreeEmitMap(orig_emap);
  FreeEmitMap(sub_emap);
  for(v = 0; v <= orig_cm->M; v++)
    {
      free(cs2hn_map_o[v]);
      free(cs2hs_map_o[v]);
    }
  free(cs2hn_map_o);
  free(cs2hs_map_o);
  for(v = 0; v <= sub_cm->M; v++)
    {
      free(cs2hn_map_s[v]);
      free(cs2hs_map_s[v]);
    }
  free(cs2hn_map_s);
  free(cs2hs_map_s);

  for(k_o = 0; k_o <= cp9hmm_o->M; k_o++)
    {
      for(ks_s = 0; ks_s < 3; ks_s++)
	free(hns2cs_map_o[k_o][ks_s]);
      free(hns2cs_map_o[k_o]);
    }
  free(hns2cs_map_o);
  FreeCPlan9(cp9hmm_o);

  for(k_s = 0; k_s <= cp9hmm_s->M; k_s++)
    {
      for(ks_s = 0; ks_s < 3; ks_s++)
	free(hns2cs_map_s[k_s][ks_s]);
      free(hns2cs_map_s[k_s]);
    }
  free(hns2cs_map_s);
  FreeCPlan9(cp9hmm_s);

  free(sttypes);
  free(nodetypes);

  *ret_sub2orig_id   = sub2orig_id;
  *ret_orig2sub_smap = orig2sub_smap;
  *ret_sub2orig_smap = sub2orig_smap;

  return;
}

/**************************************************************************
 * EPN 09.11.06
 * cm2sub_cm_check_id_next_node()
 *
 * It's known that sub_nd and orig_nd are identical nodes, in that they 
 * are of the same type and model the same column of the seed alignment.
 * In this function, we check if the next node of the sub_cm and orig_nd
 * have the same relationship and if so, we update information in 
 * the sub2orig_id array, setting the sub_cm states of this node to 
 * TRUE.
 * 
 * Returns: void
 *
 * Args:    
 * CM_t  orig_cm
 * CM_t  sub_cm
 * int   orig_nd
 * int   sub_nd
 * int *sub2orig_id
 * int* *node_cc_left_o
 * int* *node_cc_right_o
 * int* *node_cc_left_s
 * int* *node_cc_right_s
 * int   sub_start
 */

static void
cm2sub_cm_check_id_next_node(CM_t *orig_cm, CM_t *sub_cm, int orig_nd, int sub_nd,
			     int *sub2orig_id, int *node_cc_left_o, int *node_cc_right_o,
			     int *node_cc_left_s, int *node_cc_right_s, int sub_start,
			     int print_flag)
{
  int v_s;

  if(print_flag) printf("in cm2sub_cm_check_id_next_node()\n");

  if((orig_nd+1) > (orig_cm->nodes-1))
    return;
  if((sub_nd+1)  > (sub_cm->nodes-1))
    return;

  if(orig_cm->ndtype[orig_nd+1] != sub_cm->ndtype[sub_nd+1])
    {
      if(print_flag) printf("cm2sub_cm_check_id_next_node: %d return 1\n", sub_nd);
      return;
    }
  if(node_cc_left_o[orig_nd+1] == -1 & node_cc_left_s[sub_nd+1] != -1)
    {
      if(print_flag) printf("cm2sub_cm_check_id_next_node: %d return 2\n", sub_nd);
      return;
    }
  if(node_cc_right_o[orig_nd+1] == -1 & node_cc_right_s[sub_nd+1] != -1)
    {
      if(print_flag) printf("cm2sub_cm_check_id_next_node: %d return 3\n", sub_nd);
      return;
    }
  if(node_cc_left_o[orig_nd+1] != (node_cc_left_s[sub_nd+1] + sub_start - 1))
    {
      if(print_flag) printf("cm2sub_cm_check_id_next_node: %d return 4\n", sub_nd);
      return;
    }
  if(node_cc_right_o[orig_nd+1] != (node_cc_right_s[sub_nd+1] + sub_start - 1))
    {
      if(print_flag) printf("cm2sub_cm_check_id_next_node: %d return 5\n", sub_nd);
      return;
    }
  /* if we get here, the next node is identical */
  v_s = sub_cm->nodemap[sub_nd];
  while(sub_cm->ndidx[v_s] == sub_nd)
    {
      if(print_flag) printf("setting sub2orig_id[v_s:%d] to TRUE\n", v_s);
      sub2orig_id[v_s++] = TRUE;
    }
}  
  
/**************************************************************************
 * EPN 10.05.06
 * cm2sub_cm_find_impossible_misc_cases
 *
 * For certain situations, the conversion of an orig_cm to a sub_cm loses
 * some information that makes it impossible for a CP9 trained from the sub_cm
 * to exactly match a CP9 trained from the orig_cm for the corresponding
 * columns. One case where it is impossible involves start states as follows:
 *
 * if for any k k=spos..epos-1
 * X >= 0 start states exists in the orig_cm in a node between: 
 *          cc_node_map_o[k] -> cc_node_map_o[k+1]
 * AND Y >= 1 start states exist in the sub_cm in a node between:
 *          cc_node_map_s[k-spos+1] -> cc_node_map_s[k-spos+1+1]
 * where X != Y 
 *
 * AND further one or both of the two nodes in the sub_cm (cc_node_map_s[k-spos+1] 
 * OR cc_node_map_s[k-spos+1+1] must be a MATP. 
 * 
 * This is because for the sub_cm paths that would go from CP9 node k to k+1
 * were forced to go through a start state where as for the orig_cm there
 * was not a requirement to go through a start. Therefore when the
 * sub_cm was constructed it lost some information about the original
 * transitions.
 *
 * Also, a special situation of this case occurs with transitions out of
 * CP9 node k=spos-1 to k=spos, and out of CP9 node k=epos to k=epos+1. 
 * The sub_cm node that models columns spos-1 and epos-1 is the ROOT node, 
 * which only really models inserts in those columns. It turns out that
 * the transition distributions will nearly always be screwy out of these
 * two nodes, except in the case when the sub_cm ROOT_IL and ROOT_IR map
 * to an IL and IR state respectivley in the original *within the same 
 * orig_cm node*. In which case the transitions out of node 0 will be 
 * identical to those out of node spos-1 in the orig_cm. 
 * There are a few other rare cases where the transitions will be 
 * identical, but an exhaustive understanding of them eludes me and
 * it would take too long to explain those that I do know of (see
 * ~nawrockie/notebook/6_0725_inf_sub_cm/00LOG for more).
 *
 * Returns: void
 *
 * Args:    
 * CM_t  orig_cm
 * CM_t  sub_cm
 * int *cc_node_map_o
 * int *cc_node_map_s
 * int *imp_cc
 * int spos;
 * int epos;
 */

static void
cm2sub_cm_find_impossible_misc_cases(CM_t *orig_cm, CM_t *sub_cm, int **orig2sub_smap, int **sub2orig_smap, 
				     int *cc_node_map_o, int *node_cc_left_o, int *node_cc_right_o,
				     int *cc_node_map_s, int *node_cc_left_s, int *node_cc_right_s, int *imp_cc, int spos, int epos,
				     int print_flag)
{
  int k;
  int sub_starts; 
  int orig_starts;
  int orig_nd1;
  int orig_nd2;
  int sub_nd1;
  int sub_nd2;
  int temp;
  int nd;
  int orig_cm_M;
  int orig_special_matps;
  int sub_special_matl;
  int sub_both_matps;
  CMEmitMap_t *orig_emap;         /* consensus emit map for the original, template CM */
  CMEmitMap_t *sub_emap;          /* consensus emit map for the sub CM */
  int case7_flag;
  int curr_sub_nd;
  int curr_orig_nd;
  int orig_next_start;

  int orig_il1;
  int orig_il2;
  int orig_ir1;
  int orig_ir2;
  int other_il;
  int other_ir;

  char **nodetypes;
  nodetypes = malloc(sizeof(char *) * 8);
  nodetypes[0] = "BIF";
  nodetypes[1] = "MATP";
  nodetypes[2] = "MATL";
  nodetypes[3] = "MATR";
  nodetypes[4] = "BEGL";
  nodetypes[5] = "BEGR";
  nodetypes[6] = "ROOT";
  nodetypes[7] = "END";

  if(print_flag) printf("in cm2sub_find_impossible_misc_cases()\n");

  orig_emap = CreateEmitMap(orig_cm);
  sub_emap = CreateEmitMap(sub_cm);

  /* We know that the sub CM node 0 and node epos-spos+1 are nearly always going to 
   * get the transition distributions wrong (see comments above in the function
   * explanation). There are a couple of cases they should get right 
   * (see ~nawrockie/notebook/6_0725_inf_sub_cm/00LOG for details), we 
   * check for these here:
   */
  imp_cc[0] = 1;
  imp_cc[(epos-spos+1)] = 1;
  /* check if the orig_cm states that model sub_cm ROOT_IL and ROOT_IR are
   * from the same node. */
  orig_il1 = sub2orig_smap[1][0]; /* 1st of up to 2 states that maps to sub_cm's ROOT_IL */
  orig_il2 = sub2orig_smap[1][1]; /* 2nd state that maps to sub_cm's ROOT_IL or -1 if only 1 maps*/
  orig_ir1 = sub2orig_smap[2][0]; /* 1st of up to 2 states that maps to sub_cm's ROOT_IR */
  orig_ir2 = sub2orig_smap[2][1]; /* 2nd state that maps to sub_cm's ROOT_IR or -1 if only 1 maps*/

  /* We ASSUME that ambiguities have been removed, i.e. if two insert states map to either ROOT_IL
   * or ROOT_IR, one of them has been detached. We exploit this knowledge.
   */
  if(orig_il2 != -1)
    {
      if(orig_cm->sttype[orig_il1+1] == E_st)
	orig_il1 = orig_il2; /* orig_il1 was detached */
      else if(orig_cm->sttype[orig_il2+1] == E_st)
	{
	  /* do nothing */
	}
      else
	Die("ERROR, can't determine which state was detached\n");
    }
  if(orig_ir2 != -1)
    {
      if(orig_cm->sttype[orig_ir1+1] == E_st)
	orig_ir1 = orig_ir2; /* orig_ir1 was detached */
      else if(orig_cm->sttype[orig_ir2+1] == E_st)
	{
	  /* do nothing */
	}
      else
	Die("ERROR, can't determine which state was detached\n");
    }

  /* Now orig_il1 and orig_ir1 map to the ONLY insert states that map to sub_cm 
   * ROOT_IL and ROOT_IR respectively. */
  if(orig_cm->ndidx[orig_il1] == orig_cm->ndidx[orig_ir1])
    {
      /* we can get these guys right in this case */
      imp_cc[0] = FALSE; 
      imp_cc[(epos-spos+1)] = FALSE;
    }
  /* Check for the other case we know about where transition
   * distribution out of epos-spos+1 should be okay */
  else if((orig_cm->ndtype[orig_cm->ndidx[orig_il1]] == ROOT_nd ||
	   orig_cm->ndtype[orig_cm->ndidx[orig_il1]] == MATP_nd) && 
	  (orig_cm->ndtype[orig_cm->ndidx[orig_ir1]] == ROOT_nd ||
	   orig_cm->ndtype[orig_cm->ndidx[orig_ir1]] == MATP_nd))
    {
      if(orig_cm->sttype[orig_il1] == IL_st && 
	 orig_cm->sttype[orig_ir1] == IR_st)
	{
	  other_il = orig_il1 + 1;
	  other_ir = orig_ir1 - 1;
	  if(other_il > orig_ir1 && other_ir < orig_il1)
	    {
	      imp_cc[(epos-spos+1)] == FALSE;
	    }
	}
    }

  for(k = 0; k <= (epos-spos); k++)
    {
      if(k == 0 && imp_cc[0])
	continue;

      if(print_flag) printf("k: %d\n", k);
      if((k+spos-1) == 0)
	orig_nd1 = 0;
      else
	orig_nd1 = cc_node_map_o[k+spos-1];
      orig_nd2 = cc_node_map_o[k+spos-1+1];

      orig_special_matps = FALSE;
      if((orig_cm->ndtype[orig_nd1] == MATP_nd && 
	  orig_cm->ndtype[orig_nd2] == MATP_nd) &&
	 (node_cc_right_o[orig_nd1] == (k+spos-1) && 
	  node_cc_left_o[orig_nd2]  == (k+spos)))
	{
	  if((node_cc_left_o[orig_nd1] < spos) ||
	     (node_cc_right_o[orig_nd2] > epos))
	  {	
	    /* This is a special case */
	    orig_special_matps = TRUE;
	  }
	}

      if(orig_nd2 < orig_nd1)
      {
	temp = orig_nd1;
	orig_nd1 = orig_nd2;
	orig_nd2 = temp;
      }
      orig_starts = 0;

      if(print_flag) printf("orig_nd1: %d | orig_nd2: %d\n", orig_nd1, orig_nd2);
      for(nd = orig_nd1; nd <= orig_nd2; nd++)
	{
	  if(print_flag) printf("orig_cm->ndtype[%d]: %s L: %4d R: %4d (spos: %4d) (epos: %4d)\n", nd, nodetypes[orig_cm->ndtype[nd]], orig_emap->lpos[nd], orig_emap->rpos[nd], spos, epos);
	  if(orig_cm->ndtype[nd] == BEGL_nd || 
	     orig_cm->ndtype[nd] == BEGR_nd)
	    { orig_starts++; }
	}
      if(k == 0)
	sub_nd1 = 0;
      else
	sub_nd1 = cc_node_map_s[k];
      sub_nd2 = cc_node_map_s[k+1];

      sub_special_matl = FALSE;
      if(sub_cm->ndtype[sub_nd1] == MATL_nd)
	sub_special_matl = TRUE;

      if(sub_nd2 < sub_nd1)
	{
	  temp = sub_nd1;
	  sub_nd1 = sub_nd2;
	  sub_nd2 = temp;
	}

      sub_both_matps = FALSE;
      if((sub_cm->ndtype[sub_nd1] == MATP_nd && sub_cm->ndtype[sub_nd2] == MATP_nd) &&
	 (!(node_cc_left_s[sub_nd1] < node_cc_left_s[sub_nd2] && node_cc_right_s[sub_nd1] > node_cc_right_s[sub_nd2])))
	sub_both_matps = TRUE;

      sub_starts = 0;
      if(print_flag) printf("sub_nd1: %d | sub_nd2: %d\n", sub_nd1, sub_nd2);
      for(nd = sub_nd1; nd <= sub_nd2; nd++)
	{
	  if(print_flag) printf("sub_cm->ndtype[%d]: %s L: %4d R: %4d (spos: %4d) (epos: %4d)\n", nd, nodetypes[sub_cm->ndtype[nd]], (sub_emap->lpos[nd]+spos-1), (sub_emap->rpos[nd]+spos-1), spos, epos);
	  if(sub_cm->ndtype[nd] == BEGL_nd || sub_cm->ndtype[nd] == BEGR_nd)
	    {
	      sub_starts++;
	      if(sub_cm->ndtype[sub_nd1] != MATP_nd && sub_cm->ndtype[sub_nd2] != MATP_nd)
		Die("ERROR in cm2sub_cm_find_impossible_misc_cases() found impossible case not involving any MATP in the sub_cm, k: %d epos-spos+1: %d\n", k, (epos-spos+1));
	      if(orig_cm->ndtype[orig_nd1] != MATP_nd && orig_cm->ndtype[orig_nd2] != MATP_nd)
		Die("ERROR in cm2sub_cm_find_impossible_misc_cases() found impossible case not involving any MATP in the orig_cm\n, k: %d | epos-spos+1: %d", k, (epos-spos+1));
	    }	    
	}

      if(print_flag) printf("sub_starts: %d | orig_starts: %d\n", sub_starts, orig_starts);

      if(sub_starts > 0 && orig_starts == 0)
	imp_cc[k] = 3;
      else if((sub_starts > 0 && sub_starts > orig_starts) &&
	      sub_both_matps == TRUE)
	imp_cc[k] = 4;
      else if(sub_starts > 0 && orig_special_matps && sub_special_matl)
	imp_cc[k] = 5;
      else if(sub_starts == 0 && orig_starts > 0 && orig_special_matps && sub_special_matl)
	imp_cc[k] = 6;

    }  
  FreeEmitMap(orig_emap);
  FreeEmitMap(sub_emap);
  free(nodetypes);
}  

/**************************************************************************
 * EPN 10.18.06
 * cm2sub_cm_find_impossible_matr_cases
 *
 * For certain situations, the conversion of an orig_cm to a sub_cm loses
 * some information that makes it impossible for a CP9 trained from the sub_cm
 * to exactly match a CP9 trained from the orig_cm for the corresponding
 * columns. One case where it is impossible involves MATR nodes as follows:
 *
 * if for any k k=spos..epos-1:
 * cc_node_map_o[k] + 1  = cc_node_map_o[k+1] AND  (NOT TRUE!)
 * cc_node_map_s[k] + 1 != cc_node_map_s[k+1] AND
 * cc_node_map_o[k]     == cc_node_map_s[k]   == MATL AND
 * cc_node_map_o[k+1]   == cc_node_map_s[k+1] == MATP AND
 * all nodes between cc_node_map_s[k] .. cc_node_map[k+1] are BIF, BEG*, or MATR nodes AND
 * the orig_cm nodes that the stretch of MATR nodes map to are all MATP nodes
 *     for which the other half of the nodes map to positions outside the sub_cm.
 * 
 * Here we explicitly check for these situations and set imp_cc[k] = TRUE for
 * any k that satisfy all the criteria.
 * Returns: void
 *
 * Args:    
 * CM_t  orig_cm
 * CM_t  sub_cm
 * int *cc_node_map_o
 * int *cc_node_map_s
 * int *imp_cc
 * int spos;
 * int epos;
 */

static void
cm2sub_cm_find_impossible_matr_cases(CM_t *orig_cm, CM_t *sub_cm, int *cc_node_map_o, int *node_cc_left_o,
				     int *node_cc_right_o, int *cc_node_map_s, int *node_cc_left_s, 
				     int *node_cc_right_s, int *imp_cc, int spos, int epos, int print_flag)
{
  int sub_k;
  int orig_k;
  int orig_nd;
  int sub_nd;
  int next_sub_nd;
  int next_orig_nd;
  int sub_nd1;
  int sub_nd2;
  int orig_nd1;
  int orig_nd2;
  int tmp_k;
  int tmp_sub_nd;
  int tmp_orig_nd;
  int current_case;
  int correct_node_types_flag;
  int orig_matr_stretch_flag;
  int sub_matr_stretch_flag;

  char **nodetypes;
  nodetypes = malloc(sizeof(char *) * 8);
  nodetypes[0] = "BIF";
  nodetypes[1] = "MATP";
  nodetypes[2] = "MATL";
  nodetypes[3] = "MATR";
  nodetypes[4] = "BEGL";
  nodetypes[5] = "BEGR";
  nodetypes[6] = "ROOT";
  nodetypes[7] = "END";

  for(sub_k = 1; sub_k <= (epos-spos); sub_k++)
    {
      orig_k = sub_k + spos - 1;
      if(print_flag) printf("CASE 2 k: %d\n", sub_k);
      /* Strategy: we check each of our criteria independently */
      
      /* Find out what the node that models the next column is */
      sub_nd = cc_node_map_s[sub_k];
      next_sub_nd = cc_node_map_s[sub_k+1];
      orig_nd = cc_node_map_o[orig_k];
      next_orig_nd = cc_node_map_o[orig_k+1];
      if(sub_nd < next_sub_nd)
	{
	  sub_nd1 = sub_nd;
	  sub_nd2 = next_sub_nd;
	}
      else
	{
	  sub_nd1 = next_sub_nd;
	  sub_nd2 = sub_nd;
	}
      if(orig_nd < next_orig_nd)
	{
	  orig_nd1 = orig_nd;
	  orig_nd2 = next_orig_nd;
	}
      else
	{
	  orig_nd1 = next_orig_nd;
	  orig_nd2 = orig_nd;
	}
      
      if(print_flag) printf("CASE 2 k: %d  sub_nd: %d  next_sub_nd: %d\n", sub_k, sub_nd, next_sub_nd);
      if(print_flag) printf("CASE 2 k: %d orig_nd: %d next_orig_nd: %d\n", sub_k, orig_nd, next_orig_nd);

      /* Make sure that next_sub_nd is and sub_nd are not consecutive */
      if(sub_nd1 != sub_nd2-1)
	{
	  /* Check if min(sub_nd, next_sub_nd) and min(orig_nd, next_orig_nd) are both MATLs and
	     if max(sub_nd, next_sub_nd) and max(orig_nd, next_orig_nd) are both MATPs */
	  if((orig_cm->ndtype[orig_nd1] == MATL_nd && sub_cm->ndtype[sub_nd1] == MATL_nd) && 
	     (orig_cm->ndtype[orig_nd2] == MATP_nd && sub_cm->ndtype[sub_nd2] == MATP_nd))
	    correct_node_types_flag = TRUE;
	  else
	    correct_node_types_flag = FALSE;
	  if(print_flag) printf("CASE 2 k: %d correct_node_types_flag: %d\n", sub_k, correct_node_types_flag);
	  
	  if(correct_node_types_flag == TRUE)
	    {
	      /* Determine if the next_orig_nd is the next left emitting
	       * node of the orig_cm after orig_nd */
	      orig_matr_stretch_flag = TRUE;
	      for(tmp_orig_nd = orig_nd1+1; tmp_orig_nd < orig_nd2; tmp_orig_nd++)
		{
		  if(orig_cm->ndtype[tmp_orig_nd] != MATR_nd)
		    {
		      orig_matr_stretch_flag = FALSE;
		      break;
		    }
		}
	      if(print_flag) printf("CASE 2 k: %d orig_matr_stretch_flag: %d\n", sub_k, orig_matr_stretch_flag);
	      if(orig_matr_stretch_flag == TRUE)
	      {
		  /* Check if all the sub_cm nodes between sub_nd and 
		   * next_sub_nd are MATRs. */
		  sub_matr_stretch_flag = TRUE;
		  for(tmp_sub_nd = sub_nd1+1; tmp_sub_nd < sub_nd2; tmp_sub_nd++)
		    if(sub_cm->ndtype[tmp_sub_nd] != MATR_nd)
		      {
			sub_matr_stretch_flag = FALSE;
			break;
		      }
		  if(print_flag) printf("CASE 2 k: %d sub_matr_stretch_flag: %d\n", sub_k, sub_matr_stretch_flag);
		  if(sub_matr_stretch_flag == TRUE)
		    {
		      /* This should be a MATR impossible case, 
		       * Check all the criteria we *think* are always true in this situation,
		       * Die if what we think is wrong. 
		       *
		       * The orig_cm nodes that map to the same consensus columns
		       * as the MATR nodes in the sub_cm must either be:

		       * as the MATR nodes in the sub_cm must either be:
		       * 
		       * 1 orig_cm MATPs or MATRs for which the RIGHT half maps to the same column
		       *   modelled by the corresponding sub_cm MATR, and in the case of the 
		       *   MATPs the LEFT half maps to consensus columns before spos.
		       * 2 orig_cm nodes that map to the same consensus columns
		       *   as the MATR nodes in the sub_cm are ALL either orig_cm MATPs
		       *   or MATLs for which the LEFT half maps to the same 
		       *   column modelled by the corresponding sub_cm MATR, and in the
		       *   case of the MATPs, the RIGHT half maps to consensus columns 
		       *   after epos.
		       *
		       * Usually the entire set of orig_cm nodes that map to the sub_cm 
		       * MATRs are either one or the other type, but I've seen rare cases
		       * where there's a stretch of one type, and then a stretch of the
		       * other.
		       */
		      tmp_sub_nd = cc_node_map_s[sub_k] + 1;
		      tmp_k       = node_cc_right_s[tmp_sub_nd] + spos - 1;
		      tmp_orig_nd = cc_node_map_o[tmp_k];
		      
		      for(tmp_sub_nd = sub_nd1+1; tmp_sub_nd < sub_nd2; tmp_sub_nd++)
			{
			  tmp_k       = node_cc_right_s[tmp_sub_nd] + spos - 1;
			  tmp_orig_nd = cc_node_map_o[tmp_k];
			  
			  if(print_flag) printf("10.18.06: %s %3d %3d | %s %3d %3d (xcc: %3d)\n", nodetypes[sub_cm->ndtype[tmp_sub_nd]], tmp_sub_nd, (node_cc_right_s[tmp_sub_nd]+spos-1), nodetypes[orig_cm->ndtype[tmp_orig_nd]], node_cc_left_o[tmp_orig_nd], node_cc_right_o[tmp_orig_nd], tmp_k);
			  if(node_cc_right_o[tmp_orig_nd] == (node_cc_right_s[tmp_sub_nd]+spos-1)) /* Case 1 above */
			    {
			      if(orig_cm->ndtype[tmp_orig_nd] != MATP_nd && orig_cm->ndtype[tmp_orig_nd] != MATR_nd)
				Die("ERROR 2 in cm2sub_cm_find_impossible_matr_cases() found impossible MATR case that can't be classified as case 1 or case 2, k: %d | spos: %d epos: %d", sub_k, spos, epos);
			      if(orig_cm->ndtype[tmp_orig_nd] == MATP_nd && node_cc_left_o[tmp_orig_nd] >= spos)
				Die("ERROR 3 in cm2sub_cm_find_impossible_matr_cases() found impossible MATR case that can't be classified as case 1 or case 2, k: %d | spos: %d epos: %d", sub_k, spos, epos);
			    }
			  else if(node_cc_left_o[tmp_orig_nd] == (node_cc_right_s[tmp_sub_nd]+spos-1)) /* Case 2 above */
			    {
			      if(orig_cm->ndtype[tmp_orig_nd] != MATP_nd && orig_cm->ndtype[tmp_orig_nd] != MATL_nd)
				Die("ERROR 4 in cm2sub_cm_find_impossible_matr_cases() found impossible MATR case that can't be classified as case 1 or case 2, k: %d | spos: %d epos: %d", sub_k, spos, epos);
			      if(orig_cm->ndtype[tmp_orig_nd] == MATP_nd && node_cc_right_o[tmp_orig_nd] <= epos)
				Die("ERROR 5 in cm2sub_cm_find_impossible_matr_cases() found impossible MATR case that can't be classified as case 1 or case 2, k: %d | spos: %d epos: %d | node_cc_right_o[%d]: %d", sub_k, spos, epos, tmp_orig_nd, node_cc_right_o[tmp_orig_nd]);
			    }
			}
		      /* if we get here, we've satisfied all of our criteria */
		      imp_cc[sub_k] = 2;
		    }
		}
	    }
	}
    }
  free(nodetypes);
}  

/**************************************************************************
 * EPN 10.06.06
 * Function: debug_print_misc_sub_cm_info()
 **************************************************************************/
static void
debug_print_misc_sub_cm_info(CM_t *orig_cm, CM_t *sub_cm, int **orig2sub_smap, int **sub2orig_smap, int spos, int epos, 
			     int *cc_node_map_o, int *node_cc_left_o, int *node_cc_right_o)
{
  int orig_il1;
  int orig_il2;
  int orig_ir1;
  int orig_ir2;
  int orig_ss;

  char **nodetypes;
  char **sttypes;
  char **sides;

  int orig_n1_type;
  int side_idx;

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

  sides = malloc(sizeof(char *) * 3);
  sides[0] = "L";
  sides[1] = "R";
  sides[2] = "N";


  orig_il1 = sub2orig_smap[1][0]; /* 1st of up to 2 states that maps to sub_cm's ROOT_IL */
  orig_il2 = sub2orig_smap[1][1]; /* 2nd state that maps to sub_cm's ROOT_IL or -1 if only 1 maps*/
  orig_ir1 = sub2orig_smap[2][0]; /* 1st of up to 2 states that maps to sub_cm's ROOT_IR */
  orig_ir2 = sub2orig_smap[2][1]; /* 2nd state that maps to sub_cm's ROOT_IR or -1 if only 1 maps*/

  /* We ASSUME that ambiguities have been removed, i.e. if two insert states map to either ROOT_IL
   * or ROOT_IR, one of them has been detached. We exploit this knowledge.
   */
  if(orig_il2 != -1)
    {
      if(orig_cm->sttype[orig_il1+1] == E_st)
	orig_il1 = orig_il2; /* orig_il1 was detached */
      else if(orig_cm->sttype[orig_il2+1] == E_st)
	{
	  /* do nothing */
	}
      else
	Die("ERROR, can't determine which state was detached in debug_print_misc_sub_cm_info\n");
    }
  if(orig_ir2 != -1)
    {
      if(orig_cm->sttype[orig_ir1+1] == E_st)
	orig_ir1 = orig_ir2; /* orig_ir1 was detached */
      else if(orig_cm->sttype[orig_ir2+1] == E_st)
	{
	  /* do nothing */
	}
      else
	Die("ERROR, can't determine which state was detached in debug_print_misc_sub_cm_info\n");
    }

  /* Now orig_il1 and orig_ir1 map to the ONLY insert states that map to sub_cm 
   * ROOT_IL and ROOT_IR respectively.
   */
  printf("10.16.06 IL1: %3d %4s %2s | start:   %3d | end:   %3d\n", orig_il1, nodetypes[orig_cm->ndtype[orig_cm->ndidx[orig_il1]]], sttypes[orig_cm->sttype[orig_il1]], spos, epos);

  if(sub_cm->ndtype[1] == BIF_nd)
    {
      orig_n1_type = 0;
    }
  else
    {
      orig_n1_type = orig_cm->ndtype[orig_cm->ndidx[(sub2orig_smap[3][0])]];
    }

  if(orig_n1_type == MATP_nd && sub_cm->ndtype[1] != MATR_nd)
    {
      if(node_cc_left_o[orig_cm->ndidx[(sub2orig_smap[3][0])]] == spos)
	side_idx = 0;
      else if(node_cc_right_o[orig_cm->ndidx[(sub2orig_smap[3][0])]] == spos)
	side_idx = 1;
      else
	Die("ERROR MATP confusion! orig_cm node: %d | left: %d | right: %d | spos: %d\n", (orig_cm->ndidx[(sub2orig_smap[3][0])]), (node_cc_left_o[orig_cm->ndidx[(sub2orig_smap[3][0])]]), (node_cc_right_o[orig_cm->ndidx[(sub2orig_smap[3][0])]]), spos);
    }
  else if (orig_n1_type == MATP_nd && sub_cm->ndtype[1] == MATR_nd)
    {
      if(node_cc_left_o[orig_cm->ndidx[(sub2orig_smap[3][0])]] == epos)
	side_idx = 0;
      else if(node_cc_right_o[orig_cm->ndidx[(sub2orig_smap[3][0])]] == epos)
	side_idx = 1;
      else
	Die("ERROR MATP confusion! orig_cm node: %d | left: %d | right: %d | spos: %d\n", (orig_cm->ndidx[(sub2orig_smap[3][0])]), (node_cc_left_o[orig_cm->ndidx[(sub2orig_smap[3][0])]]), (node_cc_right_o[orig_cm->ndidx[(sub2orig_smap[3][0])]]), spos);
    }
  else
    {
      side_idx = 2;
    }
  printf("10.16.06 IR1: %3d %4s %2s | sub-n1: %4s | orig-n1: %4s%1s | case:   ", orig_ir1, nodetypes[orig_cm->ndtype[orig_cm->ndidx[orig_ir1]]], sttypes[orig_cm->sttype[orig_ir1]], nodetypes[sub_cm->ndtype[1]], nodetypes[orig_n1_type], sides[side_idx]);

  /* figure out 'case' of ROOT transitions */
  orig_ss = sub2orig_smap[3][0]; /* orig_ss is the 1 (of possibly 2) orig_cm states that map to the first
				  * state in sub_cm node 1 (sub_cm state 3)
				  */
  if((orig_il1 < orig_ir1) && (orig_ir1 < orig_ss))
     printf("1A\n");
  if((orig_il1 < orig_ss) && (orig_ss < orig_ir1))
     printf("1B\n");
  if((orig_ss < orig_il1) && (orig_il1 < orig_ir1))
     printf("1C\n");

  if((orig_ir1 < orig_il1) && (orig_il1 < orig_ss))
     printf("2A\n");
  if((orig_ir1 < orig_ss) && (orig_ss < orig_il1))
     printf("2B\n");
  if((orig_ss < orig_ir1) && (orig_ir1 < orig_il1))
     printf("2C\n");

  printf("\n");

  /* Begin 10.17.06 info */
  int ilmap, irmap;
  int ildual, irdual;
  int iloff, iroff;
  int other_insert_il;
  int other_insert_ir;
  int ilother_case;
  int irother_case;

  if(orig_cm->sttype[orig_il1] == IL_st)
    ilmap = 0; /* sides[0] = "L" */
  else
    ilmap = 1; /* sides[1] = "R" */
  if(orig_cm->sttype[orig_ir1] == IL_st)
    irmap = 0; /* sides[0] = "L" */
  else
    irmap = 1; /* sides[1] = "R" */

  if(orig_cm->ndtype[orig_cm->ndidx[orig_il1]] == MATP_nd ||
     orig_cm->ndtype[orig_cm->ndidx[orig_il1]] == ROOT_nd)
    ildual = TRUE;
  else
    ildual = FALSE;
  if(orig_cm->ndtype[orig_cm->ndidx[orig_ir1]] == MATP_nd ||
     orig_cm->ndtype[orig_cm->ndidx[orig_ir1]] == ROOT_nd)
    irdual = TRUE;
  else
    irdual = FALSE;

  iloff = -1;
  if(ildual == TRUE)
    {
      /* check if other insert state in orig_cm node that has insert that
       * maps to sub_cm ROOT_IL maps to a state in the sub_cm 
       */
      if(ilmap == 0) /* ROOT_IL maps to a IL */
	other_insert_il = orig_il1 + 1;
      else           /* ROOT_IL maps to a IR */
	other_insert_il = orig_il1 - 1;
      if(orig2sub_smap[other_insert_il][0] == -1 && 
	 orig2sub_smap[other_insert_il][1] == -1)
	iloff = TRUE;
      else
	iloff = FALSE;
    }

  iroff = -1;
  if(irdual == TRUE)
    {
      if(irmap == 0) /* ROOT_IR maps to a IL */
	other_insert_ir = orig_ir1 + 1;
      else           /* ROOT_IR maps to a IR */
	other_insert_ir = orig_ir1 - 1;
      if(orig2sub_smap[other_insert_ir][0] == -1 && 
	 orig2sub_smap[other_insert_ir][1] == -1)
	iroff = TRUE;
      else
	iroff = FALSE;
    }

  CMEmitMap_t *orig_emap;         /* consensus emit map for the original, template CM */
  int other_cc_il, other_cc_ir;
  orig_emap = CreateEmitMap(orig_cm);
  other_cc_il = -1;
  other_cc_ir = -1;
  if(ildual == TRUE)
    {
      if(orig_cm->sttype[other_insert_il] == IL_st) /* sub ROOT_IL maps to IR, other maps to IL */
 	{
	  other_cc_il = orig_emap->lpos[orig_cm->ndidx[other_insert_il]] + 1;
	  if(other_cc_il > spos)
	    Die("ERROR FUNKY\n");
	  ildual = 4;
	}
      else /* ROOT_IL maps to IL, other maps to IR */
	{
	  other_cc_il = orig_emap->rpos[orig_cm->ndidx[other_insert_il]] - 1;
	  if(other_cc_il < epos)
	    ildual = 1;
	  if(other_cc_il == epos)
	    ildual = 2;
	  if(other_cc_il > epos)
	    ildual = 3;
	}	    
    }
  /*printf("10.17.06 other_insert_il: %d other_cc_il: %d spos: %d epos: %d ildual: %d\n", other_insert_il, other_cc_il, spos, epos, ildual);*/
  if(irdual == TRUE)
    {
      if(orig_cm->sttype[other_insert_ir] == IL_st) /* sub ROOT_IR maps to IR, other maps to IL */
	{
	  other_cc_ir = orig_emap->lpos[orig_cm->ndidx[other_insert_ir]] + 1;
	  if(other_cc_ir > spos)
	    irdual = 1;
	  if(other_cc_ir == spos)
	    irdual = 2;
	  if(other_cc_ir < spos)
	    irdual = 3;
	}	  
      else /* ROOT_IR maps to IL, other maps to IR */
	{
	  other_cc_ir = orig_emap->rpos[orig_cm->ndidx[other_insert_ir]] - 1;
	  if(other_cc_ir < epos)
	    Die("ERROR FUNKY\n");
	  irdual = 4;
	}
    }	    
  /*printf("10.17.06 other_insert_ir: %d other_cc_ir: %d spos: %d epos: %d irdual: %d\n", other_insert_ir, other_cc_ir, spos, epos, irdual);*/
  
  printf("10.17.06 ilmap: %5s ildual: %2d iloff: %2d irmap: %5s irdual: %2d iroff: %2d subn1: %4s orign1: %4s%1s\n", sides[ilmap], ildual, iloff, sides[irmap], irdual, iroff, nodetypes[sub_cm->ndtype[1]], nodetypes[orig_n1_type], sides[side_idx]);

  int start_flag;
  int v;
  int vend;

  if(orig_il1 < orig_ir1)
    {
      v = orig_il1;
      vend = orig_ir1;
    }
  else
    {
      v = orig_ir1;
      vend = orig_il1;
    }
  start_flag = 0;
  for(; v <= vend; v++)
    if(orig_cm->sttype[v] == S_st)
      start_flag = 1;

  /*  if(start_flag)
    printf("10.17.06 STARTFLAG!\n");
  else
    printf("10.17.06 no startflag\n");
  */
}


/**************************************************************************
 * EPN 10.16.06
 * Function: check_sub_cm()
 *
 * Purpose:  Given a CM and a sub CM that is supposed to mirror 
 *           the CM as closely as possible between two given consensus
 *           columns (spos and epos), check that the sub_cm was correctly 
 *           constructed. 
 *           
 *	     1. Build a CP9 HMM (cp9_1) from the sub_cm.
 *	     2. Build a CP9 HMM (cp9_2) from the full cm.
 *	     3. Reconfig cp9_2 so start node is spos and end node is epos.
 *	     4. Check corresponding parameters of cp9_1 and cp9_2 to make
 *	        sure they're within pthresh.
 *
 * Args:    
 * CM_t *orig_cm     - the original, template CM
 * CM_t  *sub_cm     - the sub CM built from the orig_cm
 * int spos          - first consensus column in cm that hmm models (often 1)
 * int epos          -  last consensus column in cm that hmm models (often ncols)
 * float thresh      - probability threshold for parameters of cp9_1 and cp9_2
 * int *imp_cc;      - imp_cc[k] = 1 if CP9 node k is an impossible case to get 
 *		       the correct transition distros for the sub_cm. 
 * int *imp_cc;      - imp_cc[k] = 1 if CP9 node k is an impossible case to get 
 *		       the correct transition distros for the sub_cm. 
 * int *predict_ct   - array storing the number of HMM nodes we predicted 
 *                     would be wrong, for each of the 7 cases.
 * int *wrong_predict_ct - array storing the number of HMM nodes we predicted 
 *                     would be wrong, but we got right, for each of the 7 cases.
 * int print_flag    - TRUE to print useful debugging info
 * 
 * Returns: TRUE: if CM and sub CM are "close enough" (see code)
 *          FALSE: otherwise
 */
int 
check_sub_cm(CM_t *orig_cm, CM_t *sub_cm, int spos, int epos, float pthresh, int *imp_cc,
	     int *predict_ct, int *wrong_predict_ct, int print_flag)
{
  struct cplan9_s       *sub_hmm; /* constructed CP9 HMM from the sub_cm */
  struct cplan9_s       *orig_hmm; /* constructed CP9 HMM from the original cm 
				    * this will be reconfiged to match the sub_hmm */
  int ret_val;         /* return value */


  /* original CM paramters *_o */
  int *node_cc_left_o;    /* consensus column each CM node's left emission maps to
			   * [0..(cm->nodes-1)], -1 if maps to no consensus column*/
  int *node_cc_right_o;   /* consensus column each CM node's right emission maps to
			   * [0..(cm->nodes-1)], -1 if maps to no consensus column*/
  int *cc_node_map_o;     /* node that each consensus column maps to (is modelled by)
			   * [1..hmm_nmc] */
  int **cs2hn_map_o;      /* 2D CM state to HMM node map, 1st D - CM state index
			   * 2nd D - 0 or 1 (up to 2 matching HMM states), value: HMM node
			   * that contains state that maps to CM state, -1 if none.*/
  int **cs2hs_map_o;      /* 2D CM state to HMM node map, 1st D - CM state index
			   * 2nd D - 2 elements for up to 2 matching HMM states, 
			   * value: HMM STATE (0(M), 1(I), 2(D) that maps to CM state,
			   * -1 if none.
			   * For example: HMM node cs2hn_map[v][0], state cs2hs_map[v][0]
                           * maps to CM state v.*/
  int ***hns2cs_map_o;    /* 3D HMM node-state to CM state map, 1st D - HMM node index, 2nd D - 
			   * HMM state (0(M), 1(I), 2(D)), 3rd D - 2 elements for up to 
			   * 2 matching CM states, value: CM states that map, -1 if none.
			   * For example: CM states hsn2cs_map[k][0][0] and hsn2cs_map[k][0][1]
			   * map to HMM node k's match state.*/
  int orig_cp9_M;         /* Number of consensus columns modelled by the orig_CM */

  /* sub CM paramters *_s */
  int *node_cc_left_s;
  int *node_cc_right_s; 
  int *cc_node_map_s;   
  int **cs2hn_map_s;    
  int **cs2hs_map_s;    
  int ***hns2cs_map_s;  
			
  int debug_level;
  int v, k, ks;
  int cc;
  double **orig_phi;
  int *violation;
  int v_ct;
  int predict_total_ct;
  int wrong_predict_total_ct;
  int i;
  float diff;
  int nd;
  int orig_nd;

  char **nodetypes;
  nodetypes = malloc(sizeof(char *) * 8);
  nodetypes[0] = "BIF";
  nodetypes[1] = "MATP";
  nodetypes[2] = "MATL";
  nodetypes[3] = "MATR";
  nodetypes[4] = "BEGL";
  nodetypes[5] = "BEGR";
  nodetypes[6] = "ROOT";
  nodetypes[7] = "END";

  char **sttypes;
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

  ret_val = TRUE;
  v_ct = 0;
  wrong_predict_total_ct = 0;

  /* Build a CP9 HMM (cp9_1) from the sub_cm */
  sub_hmm = AllocCPlan9((epos-spos+1));
  ZeroCPlan9(sub_hmm);
  CPlan9SetNullModel(sub_hmm, sub_cm->null, 1.0); /* set p1 = 1.0 which corresponds to the CM */

  map_consensus_columns(sub_cm, sub_hmm->M, &node_cc_left_s, &node_cc_right_s,
			&cc_node_map_s, print_flag);
  
  if(print_flag) printf("sub_hmm->M: %d\n", sub_hmm->M);
  CP9_map_cm2hmm_and_hmm2cm(sub_cm, sub_hmm, node_cc_left_s, node_cc_right_s, 
			    cc_node_map_s, &cs2hn_map_s, &cs2hs_map_s, 
			    &hns2cs_map_s, print_flag);
  
  if(!(CP9_cm2wrhmm(sub_cm, sub_hmm, node_cc_left_s, node_cc_right_s, cc_node_map_s, cs2hn_map_s,
		    cs2hs_map_s, hns2cs_map_s, print_flag)))
    Die("Couldn't build a CM Plan 9 HMM from the sub CM.\n");
  

  /* Build a CP9 HMM from the full, original cm */
  orig_cp9_M = 0;
  for(v = 0; v <= orig_cm->M; v++)
    {
      if(orig_cm->stid[v] ==  MATP_MP)
	orig_cp9_M += 2;
      else if(orig_cm->stid[v] == MATL_ML || orig_cm->stid[v] == MATR_MR)
	orig_cp9_M++;
    }
  orig_hmm = AllocCPlan9(orig_cp9_M);
  ZeroCPlan9(orig_hmm);
  CPlan9SetNullModel(orig_hmm, orig_cm->null, 1.0); /* set p1 = 1.0 which corresponds to the CM */
  map_consensus_columns(orig_cm, orig_hmm->M, &node_cc_left_o, &node_cc_right_o,
			&cc_node_map_o, print_flag);
  
  if(print_flag) printf("orig_hmm->M: %d\n", orig_hmm->M);
  CP9_map_cm2hmm_and_hmm2cm(orig_cm, orig_hmm, node_cc_left_o, node_cc_right_o, 
			    cc_node_map_o, &cs2hn_map_o, &cs2hs_map_o, 
			    &hns2cs_map_o, print_flag);
  
  if(!(CP9_cm2wrhmm(orig_cm, orig_hmm, node_cc_left_o, node_cc_right_o, cc_node_map_o, cs2hn_map_o,
		    cs2hs_map_o, hns2cs_map_o, print_flag)))
    Die("Couldn't build a CM Plan 9 HMM from the sub CM.\n");
  
  /* Reconfig the orig_hmm so that it can only start in the spos node, and end from the epos node */
  /* Build the sub CP9 HMM by copying as much of the original cp9_hmm as possible */
  fill_phi_cp9(orig_hmm, &orig_phi, 1);
  CP9_reconfig2sub(orig_hmm, spos, epos, spos, epos, orig_phi);

  if(print_flag)
    {
      printf("PRINTING BUILT SUB HMM PARAMS:\n");
      debug_print_cp9_params(sub_hmm);
      printf("DONE PRINTING BUILT SUB HMM PARAMS:\n");
      
      printf("PRINTING BUILT & RECONFIGED ORIG HMM PARAMS:\n");
      debug_print_cp9_params(orig_hmm);
      printf("DONE PRINTING BUILT & RECONFIGEDSAMPLED HMM PARAMS:\n");
    }

  /* Check the parameters of the two CP9 HMMs */
  if(print_flag)
    {
      printf("COMPARING CP9 HMM parameters in check_sub_cm()\n");
      printf("orig | sub\n");
    }
  violation = MallocOrDie(sizeof(int) * (epos-spos+2));
  for(k = 0; k <= sub_hmm->M; k++)
    {      
      violation[k] = FALSE;
      if(print_flag) printf("Node: %d\n", k);
      if(k > 0)
	{
	  for(i = 0; i < MAXABET; i++)
	    {
	      diff = orig_hmm->mat[(spos+k-1)][i] - sub_hmm->mat[k][i];
	      if(print_flag) printf("mat[%d][%d] = %8.5f | %8.5f | (%8.5f)\n", 0, i, orig_hmm->mat[(spos+k-1)][i], sub_hmm->mat[k][i], diff);
	      if((diff > 0 && diff > pthresh) || (diff < 0 && diff < (-1. * pthresh)))
		{
		  Die("EMISSION PROBABILITY INCORRECT!\n");
		}
	    }
	}
      for(i = 0; i < MAXABET; i++)
	{
	  diff = orig_hmm->ins[(spos+k-1)][i] - sub_hmm->ins[k][i];
	  if(print_flag) printf("ins[%d][%d] = %8.5f | %8.5f | (%8.5f)\n", 0, i, orig_hmm->ins[(spos+k-1)][i], sub_hmm->ins[k][i], diff);
	  if((diff > 0 && diff > pthresh) || (diff < 0 && diff < (-1. * pthresh)))
	    {
	      Die("EMISSION PROBABILITY INCORRECT!\n");
	    }
	}

      /* Transitions */
      if(print_flag) printf("\n");
      diff = orig_hmm->t[(spos+k-1)][CTMM] - sub_hmm->t[k][CTMM];
      if((diff > 0 && diff > pthresh) || (diff < 0 && diff < (-1. * pthresh)))
	{
	  violation[k] = TRUE;
	  if(print_flag) printf("\tCTMM[%d] = %8.5f | %8.5f | %8.5f VIOLATION\n", k, orig_hmm->t[(spos+k-1)][CTMM], sub_hmm->t[k][CTMM], diff);
	}
      else
	if(print_flag) printf("\tCTMM[%d] = %8.5f | %8.5f | %8.5f\n", k, orig_hmm->t[(spos+k-1)][CTMM], sub_hmm->t[k][CTMM], diff);
	
      diff = orig_hmm->t[(spos+k-1)][CTMI] - sub_hmm->t[k][CTMI];
      if((diff > 0 && diff > pthresh) || (diff < 0 && diff < (-1. * pthresh)))
	{
	  violation[k] = TRUE;
	  if(print_flag) printf("\tCTMI[%d] = %8.5f | %8.5f | %8.5f VIOLATION\n", k, orig_hmm->t[(spos+k-1)][CTMI], sub_hmm->t[k][CTMI], diff);
	}
      else
	if(print_flag) printf("\tCTMI[%d] = %8.5f | %8.5f | %8.5f\n", k, orig_hmm->t[(spos+k-1)][CTMI], sub_hmm->t[k][CTMI], diff);

      diff = orig_hmm->t[(spos+k-1)][CTMD] - sub_hmm->t[k][CTMD];
      if((diff > 0 && diff > pthresh) || (diff < 0 && diff < (-1. * pthresh)))
	{
	  violation[k] = TRUE;
	  if(print_flag) printf("\tCTMD[%d] = %8.5f | %8.5f | %8.5f VIOLATION\n", k, orig_hmm->t[(spos+k-1)][CTMD], sub_hmm->t[k][CTMD], diff);
	}
      else
	if(print_flag) printf("\tCTMD[%d] = %8.5f | %8.5f | %8.5f\n", k, orig_hmm->t[(spos+k-1)][CTMD], sub_hmm->t[k][CTMD], diff);

      diff = orig_hmm->t[(spos+k-1)][CTIM] - sub_hmm->t[k][CTIM];
      if((diff > 0 && diff > pthresh) || (diff < 0 && diff < (-1. * pthresh)))
	{
	  violation[k] = TRUE;
	  if(print_flag) printf("\tCTIM[%d] = %8.5f | %8.5f | %8.5f VIOLATION\n", k, orig_hmm->t[(spos+k-1)][CTIM], sub_hmm->t[k][CTIM], diff);
	}
      else
	if(print_flag) printf("\tCTIM[%d] = %8.5f | %8.5f | %8.5f\n", k, orig_hmm->t[(spos+k-1)][CTIM], sub_hmm->t[k][CTIM], diff);

      diff = orig_hmm->t[(spos+k-1)][CTII] - sub_hmm->t[k][CTII];
      if((diff > 0 && diff > pthresh) || (diff < 0 && diff < (-1. * pthresh)))
	{
	  violation[k] = TRUE;
	  if(print_flag) printf("\tCTII[%d] = %8.5f | %8.5f | %8.5f VIOLATION\n", k, orig_hmm->t[(spos+k-1)][CTII], sub_hmm->t[k][CTII], diff);
	}
      else
	if(print_flag) printf("\tCTII[%d] = %8.5f | %8.5f | %8.5f\n", k, orig_hmm->t[(spos+k-1)][CTII], sub_hmm->t[k][CTII], diff);

      diff = orig_hmm->t[(spos+k-1)][CTID] - sub_hmm->t[k][CTID];
      if((diff > 0 && diff > pthresh) || (diff < 0 && diff < (-1. * pthresh)))
	{
	  violation[k] = TRUE;
	  if(print_flag) printf("\tCTID[%d] = %8.5f | %8.5f | %8.5f VIOLATION\n", k, orig_hmm->t[(spos+k-1)][CTID], sub_hmm->t[k][CTID], diff);
	}
      else
	if(print_flag) printf("\tCTID[%d] = %8.5f | %8.5f | %8.5f\n", k, orig_hmm->t[(spos+k-1)][CTID], sub_hmm->t[k][CTID], diff);

      diff = orig_hmm->t[(spos+k-1)][CTDM] - sub_hmm->t[k][CTDM];
      if((diff > 0 && diff > pthresh) || (diff < 0 && diff < (-1. * pthresh)))
	{
	  violation[k] = TRUE;
	  if(print_flag) printf("\tCTDM[%d] = %8.5f | %8.5f | %8.5f VIOLATION\n", k, orig_hmm->t[(spos+k-1)][CTDM], sub_hmm->t[k][CTDM], diff);
	}
      else
	if(print_flag) printf("\tCTDM[%d] = %8.5f | %8.5f | %8.5f\n", k, orig_hmm->t[(spos+k-1)][CTDM], sub_hmm->t[k][CTDM], diff);

      diff = orig_hmm->t[(spos+k-1)][CTDI] - sub_hmm->t[k][CTDI];
      if((diff > 0 && diff > pthresh) || (diff < 0 && diff < (-1. * pthresh)))
	{
	  violation[k] = TRUE;
	  if(print_flag) printf("\tCTDI[%d] = %8.5f | %8.5f | %8.5f VIOLATION\n", k, orig_hmm->t[(spos+k-1)][CTDI], sub_hmm->t[k][CTDI], diff);
	}
      else
	if(print_flag) printf("\tCTDI[%d] = %8.5f | %8.5f | %8.5f\n", k, orig_hmm->t[(spos+k-1)][CTDI], sub_hmm->t[k][CTDI], diff);

      diff = orig_hmm->t[(spos+k-1)][CTDD] - sub_hmm->t[k][CTDD];
      if((diff > 0 && diff > pthresh) || (diff < 0 && diff < (-1. * pthresh)))
	{
	  violation[k] = TRUE;
	  if(print_flag) printf("\tCTDD[%d] = %8.5f | %8.5f | %8.5f VIOLATION\n", k, orig_hmm->t[(spos+k-1)][CTDD], sub_hmm->t[k][CTDD], diff);
	}
      else
	if(print_flag) printf("\tCTDD[%d] = %8.5f | %8.5f | %8.5f\n", k, orig_hmm->t[(spos+k-1)][CTDD], sub_hmm->t[k][CTDD], diff);

      if(k > 0)
	{
	  diff = orig_hmm->begin[(spos+k-1)] - sub_hmm->begin[k];
	  if((diff > 0 && diff > pthresh) || (diff < 0 && diff < (-1. * pthresh)))
	    {
	      violation[0] = TRUE; /* begin actually has to do with the transition distro out of node 0 */
	      if(print_flag) printf("\t beg[%d] = %8.5f | %8.5f | %8.5f VIOLATION\n", k, orig_hmm->begin[(spos+k-1)], sub_hmm->begin[k], diff);
	    }
	  else
	    if(print_flag) printf("\t beg[%d] = %8.5f | %8.5f | %8.5f\n", k, orig_hmm->begin[(spos+k-1)], sub_hmm->begin[k], diff);

	  diff = orig_hmm->end[(spos+k-1)] - sub_hmm->end[k];
	  if((diff > 0 && diff > pthresh) || (diff < 0 && diff < (-1. * pthresh)))
	    {
	      violation[k] = TRUE;
	      if(print_flag) printf("\t end[%d] = %8.5f | %8.5f | %8.5f VIOLATION\n", k, orig_hmm->end[(spos+k-1)], sub_hmm->end[k], diff);
	    }
	  else
	    if(print_flag) printf("\t end[%d] = %8.5f | %8.5f | %8.5f\n", k, orig_hmm->end[(spos+k-1)], sub_hmm->end[k], diff);
	}
    }

  /* Add-up the violations */
  for(nd = 0; nd <= (epos-spos+1); nd++)
    {
      if(violation[nd] && imp_cc[nd] == 0)
	{
	  v_ct++;
	  printf("VIOLATION[%3d]: TRUE | spos: %3d | epos: %3d | imp_cc: %d\n", nd, spos, epos, imp_cc[nd]);
	}
      else if(violation[nd] && imp_cc[nd] != 0)
	{
	  predict_ct[imp_cc[nd]]++;
	  predict_total_ct++;
	  if(print_flag)
	    printf("PREDICTED VIOLATION[%3d]: TRUE | spos: %3d | epos: %3d | imp_cc: %d\n", nd, spos, epos, imp_cc[nd]);
	}
      else if(!violation[nd] && imp_cc[nd] != 0)
	{
	  predict_ct[imp_cc[nd]]++;
	  predict_total_ct++;
	  wrong_predict_ct[imp_cc[nd]]++;
	  wrong_predict_total_ct++;
	  if(print_flag) printf("NON-VIOLATION[%3d] %3d : spos: %3d | epos: %3d | imp_cc: %d\n", nd, wrong_predict_total_ct, spos, epos, imp_cc[nd]);
	}
    }

  /* Clean up and return */
  for(k = 0; k <= orig_hmm->M; k++)
    free(orig_phi[k]);
  free(orig_phi);

  free(violation);
  free(node_cc_left_o);
  free(node_cc_right_o);
  free(cc_node_map_o);
  for(v = 0; v <= orig_cm->M; v++)
    {
      free(cs2hn_map_o[v]);
      free(cs2hs_map_o[v]);
    }
  free(cs2hn_map_o);
  free(cs2hs_map_o);
  for(k = 0; k <= orig_hmm->M; k++)
    {
      for(ks = 0; ks < 3; ks++)
	free(hns2cs_map_o[k][ks]);
      free(hns2cs_map_o[k]);
    }
  free(hns2cs_map_o);
  FreeCPlan9(orig_hmm);

  free(node_cc_left_s);
  free(node_cc_right_s);
  free(cc_node_map_s);
  for(v = 0; v <= sub_cm->M; v++)
    {
      free(cs2hn_map_s[v]);
      free(cs2hs_map_s[v]);
    }
  free(cs2hn_map_s);
  free(cs2hs_map_s);
  for(k = 0; k <= sub_hmm->M; k++)
    {
      for(ks = 0; ks < 3; ks++)
	free(hns2cs_map_s[k][ks]);
      free(hns2cs_map_s[k]);
    }
  free(hns2cs_map_s);
  FreeCPlan9(sub_hmm);

  free(nodetypes);
  free(sttypes);

  if(v_ct > 0)
    return FALSE;
  else
    return TRUE;
}

/**************************************************************************
 * EPN 10.23.06
 * Function: sub_cm2cm_parsetree()
 *
 *
 * Returns: void
 *

 */

int
sub_cm2cm_parsetree(CM_t *orig_cm, CM_t *sub_cm, Parsetree_t **ret_orig_tr, Parsetree_t *sub_tr, int spos, 
		    int epos, int **orig2sub_smap, int **sub2orig_smap, int do_fullsub, int print_flag)
{
  int *ss_used; 
  int *ss_emitl; 
  int *ss_emitr; 
  int *il_ct; 
  int *ir_ct; 
  int *il_used; 
  int *ir_used; 
  int *tr_nd_for_bifs;
  /*int *il_emitl;
    int *ir_emitr;*/

  int x;
  int nd;
  int sub_v;
  int orig_v1; 
  int orig_v2;
  int orig_nd1;
  int orig_nd2;
  int nodes_used;
  Parsetree_t *orig_tr;
  int cm_nd;
  int emitl;
  int emitr;
  int i;
  int parent_tr_nd;
  int is_delete;
  int is_left;
  int is_right;

  Nstack_t    *pda;
  int          pos;
  int          ss;
  int          on_right;

  int emitl_flag;
  int emitr_flag;

  char **nodetypes;
  nodetypes = malloc(sizeof(char *) * 8);
  nodetypes[0] = "BIF";
  nodetypes[1] = "MATP";
  nodetypes[2] = "MATL";
  nodetypes[3] = "MATR";
  nodetypes[4] = "BEGL";
  nodetypes[5] = "BEGR";
  nodetypes[6] = "ROOT";
  nodetypes[7] = "END";

  char **sttypes;
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

  if(print_flag) printf("orig_cm nodes: %d\n", orig_cm->nodes);

  ss_used        = MallocOrDie(sizeof(int) * orig_cm->nodes + 1);
  ss_emitl       = MallocOrDie(sizeof(int) * orig_cm->nodes + 1);
  ss_emitr       = MallocOrDie(sizeof(int) * orig_cm->nodes + 1);
  il_used        = MallocOrDie(sizeof(int) * orig_cm->nodes + 1);
  ir_used        = MallocOrDie(sizeof(int) * orig_cm->nodes + 1);
  il_ct          = MallocOrDie(sizeof(int) * orig_cm->nodes + 1);
  ir_ct          = MallocOrDie(sizeof(int) * orig_cm->nodes + 1);
  /*il_emitl       = MallocOrDie(sizeof(int) * orig_cm->nodes + 1); 
    ir_emitr       = MallocOrDie(sizeof(int) * orig_cm->nodes + 1);*/
  tr_nd_for_bifs = MallocOrDie(sizeof(int) * orig_cm->nodes + 1);
  /* i*_emitl[nd] is the last residue emitted by the i* state 
   * of node nd, the first is (il_emitl[nd] - il_ct[nd] + 1)
   * or (ir_emitr[nd] + ir_ct[nd] - 1)
   */

  for(nd = 0; nd < orig_cm->nodes; nd++)
    {
      ss_used[nd]   = -1;
      ss_emitl[nd]  = -1;
      ss_emitr[nd]  = -1;
      il_ct[nd]     =  0; /* the number of times the IL state was used in the sub_cm parse */
      ir_ct[nd]     =  0; /* the number of times the IR state was used in the sub_cm parse */
      il_used[nd]   = -1;
      ir_used[nd]   = -1;
      /*il_emitl[nd]   = -1;
	ir_emitr[nd]   = -1;*/
      tr_nd_for_bifs[nd] = -1; /* this will remain -1 except for bif nodes */
    }

  for(x = 0; x < sub_tr->n; x++)
    {
      sub_v    = sub_tr->state[x];
      if(print_flag) printf("x: %d sub_v: %d\n", x, sub_v);
      orig_v1  = sub2orig_smap[sub_v][0];
      orig_v2  = sub2orig_smap[sub_v][1];
      if(print_flag) printf("orig_v1: %d | orig_v2: %d\n", orig_v1, orig_v2);
      if(orig_v1 == -1)
	{
	  if(sub_cm->sttype[sub_v] != S_st &&
	     sub_cm->sttype[sub_v] != E_st &&
	     sub_cm->sttype[sub_v] != B_st &&
	     sub_cm->sttype[sub_v] != EL_st)
	    Die("ERROR 0 in sub_cm2cm_parstree()\n");
	  continue;
	}
      orig_nd1 = orig_cm->ndidx[orig_v1];
      if(orig_v2 != -1)
	orig_nd2 = orig_cm->ndidx[orig_v2];
      else
	orig_nd2 = -1;
      
      /* We assume that no sub_cm insert states can map to 2 orig_cm inserts */
      if(orig_cm->sttype[orig_v1] == IL_st)
	{
	  il_used[orig_nd1] = orig_v1;
	  /*il_emitl[orig_nd1] = sub_tr->emitl[x];*/
	  il_ct[orig_nd1]++;
	  if(orig_v2 != -1)
	    Die("ERROR 1 in sub_cm2cm_parstree()\n");
	}
      else if(orig_cm->sttype[orig_v1] == IR_st)
	{
	  ir_used[orig_nd1] = orig_v1;
	  /*ir_emitr[orig_nd1] = sub_tr->emitr[x];*/
	  ir_ct[orig_nd1]++;
	  if(orig_v2 != -1)
	    Die("ERROR 2 in sub_cm2cm_parstree()\n");
	}
      else if(sub_cm->ndtype[sub_cm->ndidx[sub_v]] == MATP_nd)
	{
	  ss_used[orig_nd1] = orig_v1;
	  if(orig_v2 != -1)
	    Die("ERROR 3 in sub_cm2cm_parsetree()\n");
	}
      else if(orig_cm->ndtype[orig_nd1] == MATP_nd)
	{
	  if(orig_v2 == -1)
	    Die("ERROR 4 in sub_cm2cm_parsetree()\n");
	  /* We have to figure out which MATP split state sub_v corresponds to. */
	  if(sub_cm->ndtype[sub_cm->ndidx[sub_v]] != MATL_nd && 
	     sub_cm->ndtype[sub_cm->ndidx[sub_v]] != MATR_nd)
	    Die("ERROR 5 in sub_cm2cm_parsetree()\n");
	  if(orig_cm->ndtype[orig_nd2] != MATP_nd)
	    Die("ERROR 6 in sub_cm2cm_parsetree()\n");
	  
	  if(sub_cm->sttype[sub_v] == D_st)
	    {
	      if(ss_used[orig_nd1] == -1 || 
		 orig_cm->sttype[ss_used[orig_nd1]] == D_st)
		ss_used[orig_nd1] = orig_cm->nodemap[orig_nd1] + 3; /* MATP_D */
	      
	      /*else we do nothing, orig_cm->sttype[ss_used[orig_nd1]] is already
	       * either a ML_st or an MR_st */
	    }
	  else /* sub_cm->sttype[sub_v] != D_st */
	    {
	      if(ss_used[orig_nd1] == -1 || 
		 orig_cm->sttype[ss_used[orig_nd1]] == D_st)
		{
		  /* Figure out if sub_v maps to the left or right half of the MATP node */
		  if(orig_cm->sttype[orig_v1] == ML_st || 
		     orig_cm->sttype[orig_v2] == ML_st)
		    ss_used[orig_nd1] = orig_cm->nodemap[orig_nd1] + 1; /* MATP_ML */
		  else if(orig_cm->sttype[orig_v1] == MR_st || 
			  orig_cm->sttype[orig_v2] == MR_st)
		      ss_used[orig_nd1] = orig_cm->nodemap[orig_nd1] + 2; /* MATP_MR */
		  else
		    Die("ERROR 7 in sub_cm2cm_parsetree()\n");
		}
	      /* below is the only line we really need: 
		 else
		 ss_used[orig_nd1] = orig_cm->nodemap[orig_nd1]; */
	      
	      else if(orig_cm->sttype[ss_used[orig_nd1]] == ML_st) /* just for safety; should erase eventually */
		{
		  if(orig_cm->sttype[orig_v1] == MR_st || 
		     orig_cm->sttype[orig_v2] == MR_st) /* just for safety; should erase eventually */
		    {
		      ss_used[orig_nd1] = orig_cm->nodemap[orig_nd1]; /* MATP_MP */
		    }
		  else 
		    Die("ERROR 8 in sub_cm2cm_parsetree()\n");
		}
	      else if(orig_cm->sttype[ss_used[orig_nd1]] == MR_st) /* just for safety; should erase eventually */
		{
		  if(orig_cm->sttype[orig_v1] == ML_st || 
		     orig_cm->sttype[orig_v2] == ML_st) /* just for safety; should erase eventually */
		  {
		    ss_used[orig_nd1] = orig_cm->nodemap[orig_nd1]; /* MATP_MP */
		  }
		  else 
		    Die("ERROR 9 in sub_cm2cm_parsetree()\n");
		}
	    }
	}
      else
	{
	  if(orig_v2 != -1)
	    Die("ERROR 5 in sub_cm2cm_parsetree()\n");
	  ss_used[orig_nd1] = orig_v1;
	}
    }

  /* Some MATL, MATR, and MATP nodes in the orig_cm might have 0 states that map to any state 
   * used in the sub_cm parse. If we're not allowing local begins and ends, our strategy is 
   * to claim that in the orig_cm parse the D state of these nodes was used.
   */
  for(nd = 0; nd < orig_cm->nodes; nd++)
    {
      if(ss_used[nd] == -1)
	{
	  if(orig_cm->ndtype[nd] == MATP_nd)
	    ss_used[nd] = orig_cm->nodemap[nd] + 3; /* MATP_D */
	  if(orig_cm->ndtype[nd] == MATL_nd ||
	     orig_cm->ndtype[nd] == MATR_nd)
	    ss_used[nd] = orig_cm->nodemap[nd] + 1; /* MAT{L,R}_D */
	  if(orig_cm->ndtype[nd] == BIF_nd  ||
	     orig_cm->ndtype[nd] == BEGL_nd || 
	     orig_cm->ndtype[nd] == BEGR_nd ||
	     orig_cm->ndtype[nd] == END_nd)
	    ss_used[nd] = orig_cm->nodemap[nd];     /* BIF_B, END_E or BEG{L,R}_S */
	}
    }      

  /* Determine emitl and emitr for each state in the orig_cm parse */
  /* This code is noticeably cleaner than the rest - it's Sean's, from 
   * CreateEmitMap() adapted for our purposes here.
   */

  pos   = 1;
  ss    = 0;
  pda   = CreateNstack();
  PushNstack(pda, 0);		/* 0 = left side. 1 would = right side. */
  PushNstack(pda, ss);
  while (PopNstack(pda, &ss))
    {
      PopNstack(pda, &on_right);

      if (on_right) 
	{
	  pos += ir_ct[orig_cm->ndidx[ss]]; /* account for right inserts */
	  if (orig_cm->sttype[ss] == MP_st || orig_cm->sttype[ss] == MR_st) 
	    pos++;
	  ss_emitr[orig_cm->ndidx[ss]] = pos - 1;
	}
      else
	{
	  ss_emitl[orig_cm->ndidx[ss]] = pos;
	  if (orig_cm->sttype[ss] == MP_st || orig_cm->sttype[ss] == ML_st) 
	    pos++;

	  if (orig_cm->sttype[ss] == B_st)
	    {
				/* push the BIF back on for its right side  */
	      PushNstack(pda, 1);
	      PushNstack(pda, ss);
                            /* push node index for right child */
	      PushNstack(pda, 0);
	      PushNstack(pda, orig_cm->cnum[ss]);
                            /* push node index for left child */
	      PushNstack(pda, 0);
	      PushNstack(pda, orig_cm->cfirst[ss]);
	    }
	  else
	    {
				/* push the node back on for right side */
	      PushNstack(pda, 1);
	      PushNstack(pda, ss);
				/* push split state of child node on */
	      if (orig_cm->sttype[ss] != E_st) {
		PushNstack(pda, 0);
		PushNstack(pda, ss_used[orig_cm->ndidx[ss]+1]);
	      }
	    }
	  pos += il_ct[orig_cm->ndidx[ss]]; /* account for left inserts */
	}
    }      
  
  if(print_flag)
    {
      for(nd = 0; nd < orig_cm->nodes; nd++)
	{
	  printf("ss_used[%4d] (%4s) first state(%4d) | ", nd, nodetypes[orig_cm->ndtype[nd]], orig_cm->nodemap[nd]);
	  if(ss_used[nd] != -1)
	    printf("%4d (%2s) | L: %3d R: %3d\n", ss_used[nd], sttypes[orig_cm->sttype[ss_used[nd]]], ss_emitl[nd], ss_emitr[nd]);
	  else
	    printf("%4d\n", -1);
	}
      for(nd = 0; nd < orig_cm->nodes; nd++)
	printf("il_ct[%4d] (st used: %4d) ct: %4d | ir_ct[%4d] ct: %4d (st used: %4d)\n", nd, il_used[nd], il_ct[nd], nd, ir_used[nd], ir_ct[nd]);

    }

  orig_tr = CreateParsetree();
  nodes_used = 0;
  for(cm_nd = 0; cm_nd < orig_cm->nodes; cm_nd++)
    {
      emitl_flag = 0;
      emitr_flag = 0;
      if(orig_cm->sttype[ss_used[cm_nd]] == MP_st || 
	 orig_cm->sttype[ss_used[cm_nd]] == ML_st)
	emitl_flag = 1;
      if(orig_cm->sttype[ss_used[cm_nd]] == MP_st || 
	 orig_cm->sttype[ss_used[cm_nd]] == MR_st)
	emitr_flag = 1;

      /* At least 1 state in each node must be visited (if we're not in local mode) */
      if(orig_cm->ndtype[cm_nd] == BEGR_nd)
	{
	  parent_tr_nd =  tr_nd_for_bifs[orig_cm->ndidx[orig_cm->plast[orig_cm->nodemap[cm_nd]]]];
	  if(print_flag) printf("tr_nd_for_bifs[%d]\n", (orig_cm->ndidx[orig_cm->plast[orig_cm->nodemap[cm_nd]]]));
	  if(print_flag) printf("parent_tr_nd for cm_nd %d: %d\n", cm_nd, parent_tr_nd);
	  InsertTraceNode(orig_tr, parent_tr_nd, TRACE_RIGHT_CHILD, ss_emitl[cm_nd], ss_emitr[cm_nd], ss_used[cm_nd]);
	  orig_tr->nxtr[parent_tr_nd]  = orig_tr->n; /* Go back and fix nxtr for the BIF parent of this BEGR */
	}
      else
	{
	  InsertTraceNode(orig_tr, orig_tr->n-1, TRACE_LEFT_CHILD, ss_emitl[cm_nd], ss_emitr[cm_nd], ss_used[cm_nd]);
	  if(print_flag) printf("inserted trace node for orig_cm st %4s | emitl: %d | emitr: %d\n", sttypes[orig_cm->sttype[ss_used[cm_nd]]], ss_emitl[cm_nd], ss_emitr[cm_nd]);
	}

      /* Note: if we've just added a trace node for a BIF state, it's incomplete, in that it 
       * doesn't have the nextr correctly set. We'll go back and set this when we get to the 
       * right child (BEGR) of this BIF */
      if(orig_cm->ndtype[cm_nd] == BIF_nd)
	{
	  tr_nd_for_bifs[cm_nd] = orig_tr->n;
	  if(print_flag) printf("set tr_nd_for_bifs[%d]: %d\n", cm_nd, orig_tr->n);
	}

      /* Add left inserts, if any */
      for(i = 0; i < il_ct[cm_nd]; i++)
	{
	  InsertTraceNode(orig_tr, orig_tr->n-1, TRACE_LEFT_CHILD, (ss_emitl[cm_nd] + emitl_flag + i), ss_emitr[cm_nd], il_used[cm_nd]);
	  if(print_flag) printf("inserted trace node for orig_cm st %4s | emitl: %d | emitr: %d\n", sttypes[orig_cm->sttype[il_used[cm_nd]]], orig_tr->emitl[orig_tr->n-1], orig_tr->emitr[orig_tr->n-1]);
	}
      /* Add right inserts, if any */
      for(i = 0; i < ir_ct[cm_nd]; i++)
	{
	  InsertTraceNode(orig_tr, orig_tr->n-1, TRACE_LEFT_CHILD, ss_emitl[cm_nd], (ss_emitr[cm_nd] - emitr_flag - i), ir_used[cm_nd]);
	  if(print_flag) printf("inserted trace node for orig_cm st %4s | emitl: %d | emitr: %d\n", sttypes[orig_cm->sttype[ir_used[cm_nd]]], orig_tr->emitl[orig_tr->n-1], orig_tr->emitr[orig_tr->n-1]);
	}
      if(print_flag) printf("END nd: %4d | emitl: %4d | emitr: %4d\n", cm_nd, emitl, emitr);
    }      
  *ret_orig_tr = orig_tr;

  free(nodetypes);
  free(sttypes);
  return 1;
}

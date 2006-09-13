/* subCM.c
 * EPN 07.25.06
 * 
 * Building submodels (subCMs) from a template CM, with potentially less
 * structure. 
 ***************************************************************** 
 * @LICENSE@
 ***************************************************************** 
 */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>

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
static void map_orig2sub_cm_helper(CM_t *orig_cm, CM_t *sub_cm, int **orig2sub_smap, int **sub2orig_smap, int v_o, int v_s);
static void cm2sub_cm_emit_probs(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, int v_s, int v_o1, int v_o2);
static void cm2sub_cm_trans_probs(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, char ***tmap, int v_s, 
				  int **orig2sub_smap, int **sub2orig_smap);
static void cm_add_single_trans(CM_t *orig_cm, CM_t *sub_cm, int **orig2sub_smap, int **sub2orig_smap, int v_o, 
				int y_o, int v_s, int yoffset, double *orig_psi, char ***tmap);
static float cm_sum_subpaths(CM_t *orig_cm, CM_t *sub_cm, int **orig2sub_smap, int start, int end, int sub_start,
			     char ***tmap, double *orig_psi, int orig_insert1, int orig_insert2, int orig_insert3,
			     int orig_insert4);
static void debug_print_cm_params(CM_t *cm);
static int check_orig_psi_vs_sub_psi(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, double *sub_psi, 
				     int **orig2sub_smap, int **sub2orig_smap, double threshold, 
				     int print_flag);
static int orig2sub_state_check(int **orig2sub_smap, int orig_a, int orig_b);
static int sub_trans_check(CM_t *sub_cm, int **orig2sub_smap, int orig_a, int orig_b, int min_sub_v);

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
 * Function:  StripWUSS()
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

/**************************************************************************** 
 * Function:  BuildSubCM()
 * EPN 08.28.06 
 *
 * ***************WRITE DESCRIPTION!!!!!!!!!
 *
 * Args:      orig_cm    - the original model, which we're going to (potentially)
 *                         remove some structure from
 *            ret_cm     - the new subCM built from orig_cm with some structure removed.
 *            struct_start - the first position (consensus column) we want to keep structure
 *                         for.
 *            struct_end  - the last position we want to keep structure for
 *            model_start - the first position we want to model with the subCM
 *            model_end   - the last position we want to model with the subCM 
 *            orig2sub_smap - 2D state map from orig_cm (template) to sub_cm.
 *                            1st dimension - state index in orig_cm 
 *                            2nd D - 2 elements for up to 2 matching sub_cm states, 
 *            sub2orig_smap - 2D state map from orig_cm (template) to sub_cm.
 *                            1st dimension - state index in sub_cm (0..sub_cm->M-1)
 *                            2nd D - 2 elements for up to 2 matching orig_cm states, 
 * 
 * Returns:   (void)
 */
void
BuildSubCM(CM_t *orig_cm, CM_t **ret_cm, int struct_start, int struct_end, int model_start,
	   int model_end, int **orig2sub_smap, int **sub2orig_smap)
{
  CM_t            *sub_cm;       /* new covariance model, a submodel of the template */
  CMConsensus_t   *con;         /* growing consensus info for orig_cm*/
  Parsetree_t     *mtr;         /* master structure tree from the alignment*/
  char     *sub_cstr;            /* consensus substructure display string  */
  int      *sub_ct;		/* 0..con->clen-1 base pair partners array         */
  int n;
  int cpos;
  CMEmitMap_t *orig_emap;         /* consensus emit map for the original, template CM */
  CMEmitMap_t *sub_emap;         /* consensus emit map for the subCM */
  char     ***tmap;
  double    *orig_psi;              /* expected num times each state visited in template CM */
  double    *sub_psi;              /* expected num times each state visited in template CM */
  int sub_cpos;
  int v_s;
  int v_o1;
  int v_o2;
  int n_s;
  int orig_il1, orig_il2, orig_ir1, orig_ir2;
  double temp_psi;
  double temp_sum;
  int adjust_flag;

  /* Get the consensus sequence and consensus structure information from the original CM */
  con = CreateCMConsensus(orig_cm, 3.0, 1.0);
  printf("con->cseq    : %s\n", con->cseq);
  printf("con->cstr    : %s\n", con->cstr);
  printf("struct_start : %d\n", struct_start);
  printf("struct_end   : %d\n", struct_end);
  printf("model_start  : %d\n", model_start);
  printf("model_end    : %d\n", model_end);

  /* Fill a new ct array for the subCM. The subCM will only model the consensus columns
   * between model_start and model_end, and only the structure between struct_start
   * and struct_end. First copy the template (original) CMs ct array but only for the 
   * appropriate consensus columns that lie in between both structure and model boundaries
   * Next, eliminate any structure that lies outside the structure boundaries 
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
      /*printf("B sub_ct[cpos=%3d]: %3d\n", sub_cpos, sub_ct[sub_cpos]);*/
      if ((cpos+1) < struct_start || (cpos+1) > struct_end) /* cpos goes 1..clen, but ct is indexed
							     * 0..clen-1.*/
	{ 

	  if (sub_ct[sub_cpos] != -1)  sub_ct[sub_ct[sub_cpos]] = -1; /* CreateCMConsensus() uses
								       * -1 in ct[] to indicate single
								       * stranded (different convention
								       * than WUSS2ct()). */
	  sub_ct[sub_cpos] = -1;
	}
      /*printf("A sub_ct[cpos=%3d]: %3d\n\n", sub_cpos, sub_ct[sub_cpos]);*/
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
      /*printf("sub_ct[cpos=%3d]: %3d\n", cpos, sub_ct[cpos]);*/
      if(sub_ct[sub_cpos] == -1)         sub_cstr[sub_cpos] = '.'; 
      else if (sub_ct[sub_cpos]  > cpos) sub_cstr[sub_cpos] = '<';
      else if (sub_ct[sub_cpos]  < cpos) sub_cstr[sub_cpos] = '>';
      else Die("ERROR: weird error in BuildSubCM()\n");
    }
  sub_cstr[(model_end-model_start+1)] = '\0';

  /* Build the new subCM given the new consensus structure. But don't
   * parameterize it yet.
   */
  ConsensusModelmaker(sub_cstr, (model_end-model_start+1), &sub_cm, &mtr);

  printf("\n\norig struct: %s\n", con->cstr);
  printf("\n\nnew struct : %s\n", sub_cstr);

  /* Parameterize the new subCM based on the template CM. 
   * First we need the emit maps for the template CM and
   * the new subCM.
   */
  orig_emap = CreateEmitMap(orig_cm);
  sub_emap = CreateEmitMap(sub_cm);
  printf("\n\n\nDumping original CM emitmap\n");
  DumpEmitMap(stdout, orig_emap, orig_cm);
  printf("\n\n\nDumping subCM emitmap\n");
  DumpEmitMap(stdout, sub_emap, sub_cm);

  /* Map states from orig_cm to sub_cm and vice versa. */
  /*  map_orig2sub_cm(orig_cm, sub_cm, &orig2sub_smap, &sub2orig_smap, model_start,
      model_end);*/
  map_orig2sub_cm2(orig_cm, sub_cm, &orig2sub_smap, &sub2orig_smap, model_start,
		   model_end);

  /*
  printf("con->clen    : %d\n", con->clen);
  printf("con->cseq    : %s\n", con->cseq);
  printf("con->cstr    : %s\n", con->cstr);
  printf("struct_start : %d\n", struct_start);
  printf("struct_end   : %d\n", struct_end);
  printf("model_start  : %d\n", model_start);
  printf("model_end    : %d\n", model_end);
  printf("\n\norig struct: %s\n", con->cstr);
  printf("\n\nnew struct : %s\n", sub_cstr);
  */

  orig_psi = malloc(sizeof(double) * orig_cm->M);
  make_tmap(&tmap);
  fill_psi(orig_cm, orig_psi, tmap);
  CMZero(sub_cm);
  CMSetNullModel(sub_cm, orig_cm->null);

  /* Fill in emission probabilities */
  for(v_s = 0; v_s < sub_cm->M; v_s++)
    {
      if(sub_cm->sttype[v_s] != S_st &&
	 sub_cm->sttype[v_s] != D_st &&
	 sub_cm->sttype[v_s] != B_st &&
	 sub_cm->sttype[v_s] != E_st &&
	 sub_cm->sttype[v_s] != EL_st)
	cm2sub_cm_emit_probs(orig_cm, sub_cm, orig_psi, v_s, sub2orig_smap[v_s][0], sub2orig_smap[v_s][1]);
    }
  /* Fill in transition probabilities */
  for(v_s = 0; v_s < sub_cm->M; v_s++)
    {
      if(v_s == 0 || 
	 (sub_cm->sttype[v_s] != S_st &&
	  sub_cm->sttype[v_s] != B_st &&
	  sub_cm->sttype[v_s] != E_st))
	cm2sub_cm_trans_probs(orig_cm, sub_cm, orig_psi, tmap, v_s, orig2sub_smap, sub2orig_smap);
    }


  printf("\n\n\nAddressing 090806\n");
  /* Address problem 090806 (in the 00LOG of ~/notebook/6_0725_inf_destruct), by
   * retraversing the structure and subtracting out subpaths that have been counted twice
   * for a special situation involving the two inserts of ROOT and MATP states */
  for(n_s = 0; n_s < sub_cm->nodes; n_s++)
    {
      adjust_flag = FALSE;
      if(sub_cm->ndtype[n_s] == ROOT_nd)
	{
	  v_s = sub_cm->nodemap[n_s]; /* v_s == ROOT_s */
	  orig_il1 = sub2orig_smap[(v_s+1)][0];
	  orig_il2 = sub2orig_smap[(v_s+1)][1];
	  
	  orig_ir1 = sub2orig_smap[(v_s+2)][0];
	  orig_ir2 = sub2orig_smap[(v_s+2)][1];
	  
	  printf("\n\norig_il1: %d | orig_il2: %d | orig_ir1: %d | orig_ir2: %d\n", orig_il1, orig_il2, orig_ir1, orig_ir2);
	  printf("before t[0][1] = %f\n", sub_cm->t[v_s][1]);

	  if((orig_il1 != -1 && orig_ir1 != -1) && 
	     (orig_il1 > orig_ir1))
	    {
	      sub_cm->t[v_s][1] -= orig_psi[orig_ir1] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ir1, orig_il1, 2, tmap, orig_psi, -1, -2, -3, -4);
	      adjust_flag = TRUE;
	    }
	  if((orig_il1 != -1 && orig_ir2 != -1) && 
	     (orig_il1 > orig_ir2))
	    {
	      sub_cm->t[v_s][1] -= orig_psi[orig_ir2] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ir2, orig_il1, 2, tmap, orig_psi, -1, -2, -3, -4);
	      adjust_flag = TRUE;
	    }
	  if((orig_il2 != -1 && orig_ir1 != -1) && 
	     (orig_il2 > orig_ir1))
	    {
	      sub_cm->t[0][1] -= orig_psi[orig_ir1] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ir1, orig_il2, 2, tmap, orig_psi, -1, -2, -3, -4);
	      adjust_flag = TRUE;
	    }
	  if((orig_il2 != -1 && orig_ir1 != -1) && 
	     (orig_il2 > orig_ir2))
	    {
	      sub_cm->t[0][1] -= orig_psi[orig_ir2] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ir2, orig_il2, 2, tmap, orig_psi, -1, -2, -3, -4);
	      adjust_flag = TRUE;
	    }
	  printf("after t[0][1] = %f\n", sub_cm->t[v_s][1]);
	  
	  /* now adjust the virtual counts for transitions out of ROOT_IL */
	  if(adjust_flag)
	    {	  if(orig_ir1 != -1)
	      {
		/* determine what orig_psi[orig_il1] would be without the subpath of 
		 * orig_z -> orig_v (orig_z = orig_cm state that maps to sub_z=ROOT_IR)
		 *                  (orig_v = orig_y = orig_cm state that maps to ROOT_IL)
		 */
		/* first remove contribution of self-insert */
		temp_psi = orig_psi[orig_il1] / (1 + (orig_cm->t[orig_il1][0]/(1-(orig_cm->t[orig_il1][0]))));
		printf("0 temp_psi: %f\n", temp_psi);
		/* now subtract contribution of orig_z to orig_v */
		temp_psi -= orig_psi[orig_ir1] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, orig_ir1, orig_il1, 2, tmap, orig_psi, -1, -2, -3, -4); /* HEREHEREHERE */
		printf("1 temp_psi: %f\n", temp_psi);
		/* now add contribution of self-insert back */
		temp_psi += temp_psi * (orig_cm->t[orig_il1][0] / (1 - orig_cm->t[orig_il1][0]));
		printf("2 temp_psi: %f\n", temp_psi);
		/* finally scale the transitions out of ROOT_IL to non ROOT_IR states by temp_psi / orig_psi[orig_il1] */
		/* IF NEXT NODE IS MATL or MATR: */
		sub_cm->t[1][2] *= (temp_psi/orig_psi[orig_il1]);
		sub_cm->t[1][3] *= (temp_psi/orig_psi[orig_il1]);
		/* THIS CODE IS SET UP TO HANDLE FOUR.CM ONLY!!!*/
		printf("sub_cm->t[1][0]: %f\n", sub_cm->t[1][0]);
		printf("sub_cm->t[1][1]: %f\n", sub_cm->t[1][1]);
		printf("sub_cm->t[1][2]: %f\n", sub_cm->t[1][2]);
		printf("sub_cm->t[1][3]: %f\n", sub_cm->t[1][3]);
		temp_sum = sub_cm->t[1][0] + sub_cm->t[1][1] + sub_cm->t[1][2] + sub_cm->t[1][3];
		printf("SUM: %f\n", temp_sum);
	      }
	    }
	  /* else if MATP_nd TO DO IF ROOT_ND CASE ACTUALLY WORKS */
	}
    }
  /* Finally renormalize the CM */
  CMRenormalize(sub_cm);

  sub_psi = malloc(sizeof(double) * sub_cm->M);
  fill_psi(sub_cm, sub_psi, tmap);

  printf("\nDEBUG PRINT OF ORIG_CM PARAMETERS:\n");
  debug_print_cm_params(orig_cm);
  printf("\nDEBUG PRINT OF SUB_CM PARAMETERS:\n");
  debug_print_cm_params(sub_cm);

  printf("\n\ncalling check_orig_psi_vs_sub_psi()\n");
  check_orig_psi_vs_sub_psi(orig_cm, sub_cm, orig_psi, sub_psi, orig2sub_smap, sub2orig_smap, 0.001, TRUE);
  printf("\n\ndone check_orig_psi_vs_sub_psi()\n");
  exit(1);
  

  free(sub_cstr);
  free(sub_ct);
  FreeEmitMap(orig_emap);
  FreeEmitMap(sub_emap);

  *ret_cm = sub_cm;
  return;
}

/* Function: ConsensusModelmaker()
 * EPN 08.29.06 based closely on HandModelMaker:
 *              SRE 29 Feb 2000 [Seattle]; from COVE 2.0 code
 * 
 * Purpose:  Construct a model given a stated structure. The structure
 *           is provided via a "ss_cons" (consensus structure) line, as would
 *           occur in an annotated SELEX or Stockholm file. Only > and < characters
 *           in this line are interpreted (as base pairs). Pseudoknots, 
 *           if annotated, are ignored.
 *           
 *           All positions/columns of the given structure are considered 
 *           consensus, and this is the difference b/t this function and
 *           HandModelmaker. Also, this function does not take in a MSA 
 *           data structure. It was originally written for building a 
 *           new CM with less structure (MATPs) than a template CM.
 *           
 * Args:     ss_cons   - input consensus structure string 
 *           len       - length of ss_cons, number of consensus columns
 *           ret_cm    - RETURN: new model                      (maybe NULL)
 *           ret_gtr   - RETURN: guide tree for alignment (maybe NULL)
 *           
 * Return:   void
 *           cm is allocated here. FreeCM(*ret_cm).
 *           gtr is allocated here. FreeTrace().
 */
void
ConsensusModelmaker(char *ss_cons, int clen, CM_t **ret_cm, Parsetree_t **ret_gtr)
{
  CM_t           *cm;		/* new covariance model                       */
  Parsetree_t    *gtr;		/* guide tree for alignment                   */
  Nstack_t       *pda;		/* pushdown stack used in building gtr        */
  int            *matassign;	/* 0..alen-1 array; 0=insert col, 1=match col */
  int            *ct;		/* 0..alen-1 base pair partners array         */
  int             apos;		/* counter over columns of alignment          */
  int             idx;		/* counter over sequences in the alignment    */
  int             v;		/* index of current node                      */
  int             i,j,k;	/* subsequence indices                        */
  int  type;			/* type of node we're working on              */
  int  diff, bestdiff, bestk;   /* used while finding optimal split points    */   
  int  nnodes;			/* number of nodes in CM                      */
  int  nstates;			/* number of states in CM                     */

  if (ss_cons == NULL)
    Die("No consensus structure annotation available in ConsensusModelmaker().");

  /* 1. Determine a "ct" array, base-pairing partners for each position.
   *    Disallow/ignore pseudoknots. (That's what the FALSE flag does.)
   *    ct[] values give the index of a base pairing partner, or 0 for unpaired positions.
   *    Even though ss_cons is in the 0..clen-1 coord system of msa, ct[]
   *    comes back in the 1..alen coord system of the sequence.
   */
  printf("in ConsensusModelmaker: clen: %d,  %s\n", clen, ss_cons);

  if (! WUSS2ct(ss_cons, clen, FALSE, &ct))  
    Die("Consensus structure string is inconsistent"); 

  /* 2. Construct a guide tree. 
   *    This codes is borrowed from HandModelmaker(), where it
   *    was originally borrowed from yarn's KHS2Trace().
   *    
   *    We also keep track of how many states we'll need in the final CM,
   *    so we'll know how much to allocate -- and the number of nodes,
   *    for informational purposes.
   */
  nstates = nnodes = 0;
  gtr = CreateParsetree();	/* the parse tree we'll grow        */
  pda = CreateNstack();		/* a pushdown stack for our indices */

  /* Construction strategy has to make sure we number the nodes in
   * preorder traversal: for bifurcations, we can't attach the right 
   * child until we've fully traversed the left side. Therefore, we have
   * to push what we intend to attach, and pop it later. And since we
   * don't know an index for the node until we attach it, we have no
   * place to put the node's data except the stack -- so we have to
   * push several numbers onto the stack: what type of node, what
   * subseq it's responsible for (emitl...emitr), and what node
   * index it attaches to.
   * 
   * Note that we have to deal with the fact that ct is off-by-one
   * in both indices and values: e.g. the base pairing partner 
   * j of residue i is ct[i-1]-1. 
   */
  PushNstack(pda, -1);		/* what node it's attached to */
  PushNstack(pda, 1);		/* emitl */
  PushNstack(pda, clen);	/* emitr */
  PushNstack(pda, ROOT_nd);	/* "state" (e.g. node type) */

  while (PopNstack(pda, &type))	/* pop a node type to attach */
    {
      PopNstack(pda, &j);
      PopNstack(pda, &i);	/* i..j == subseq we're responsible for */
      PopNstack(pda, &v);	/* v = index of parent node in gtr */

      /* This node accounts for i..j, but we usually don't know how yet.
       * Six possibilities:
       *    i > j; this is an END state; do nothing.
       *    this is already assigned as a BEGIN; push i,j
       *    i is unpaired; this is a MATL state; push i+1, j
       *    j is unpaired; this is a MATR state; push i,j-1
       *    i,j pair to each other; this is a MATP state; push i+1,j-1
       *    i,j pair but not to each other; this is a BIFURC state;
       *        pick mid ip <= mid < jp; push BEGIN i,mid and working i,mid,
       *        and push BEGIN mid+1,j and working mid+1,j
       */
      if (i > j) {
	v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, END_nd);
	nstates += 1;		/* END_nd -> E_st */
	nnodes++;
      }

      else if (type == ROOT_nd) { /* try to push i,j; but deal with INSL and INSR */
	v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, ROOT_nd);
	PushNstack(pda, v);	/* here v==0 always. */
	PushNstack(pda, i);
	PushNstack(pda, j);
	PushNstack(pda, DUMMY_nd); /* we don't know yet what the next node will be */
	nstates += 3;		/* ROOT_nd -> S_st, IL_st, IR_st */
	nnodes++;
      }

      else if (type == BEGL_nd) {    /* no inserts */
	v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, BEGL_nd);
	PushNstack(pda, v);	
	PushNstack(pda, i);
	PushNstack(pda, j);
	PushNstack(pda, DUMMY_nd); /* we don't know yet what the next node will be */
	nstates += 1;		/* BEGL_nd -> S_st */
	nnodes++;
      }

      else if (type == BEGR_nd)  { /* look for INSL */
	v = InsertTraceNode(gtr, v, TRACE_RIGHT_CHILD, i, j, BEGR_nd);
	PushNstack(pda, v);	
	PushNstack(pda, i);
	PushNstack(pda, j);
	PushNstack(pda, DUMMY_nd); /* we don't know yet what the next node will be */
	nstates += 2;		/* BEGR_nd -> S_st IL_st */
	nnodes++;
      }

      else if (ct[i] == 0) {
	 	/* i unpaired. This is a MATL node; allow INSL */
	v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, MATL_nd);
	i++;
	PushNstack(pda, v);
	PushNstack(pda, i);
	PushNstack(pda, j);
	PushNstack(pda, DUMMY_nd); /* we don't know yet what the next node will be */
	nstates += 3;		/* MATL_nd -> ML_st, D_st, IL_st */
	nnodes++;
      }

      else if (ct[j] == 0) { 	/* j unpaired. MATR node. Deal with INSR */
	v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, MATR_nd);
	j--;
	PushNstack(pda, v);
	PushNstack(pda, i);
	PushNstack(pda, j);
	PushNstack(pda, DUMMY_nd); /* we don't know yet what the next node will be */
	nstates += 3;		/* MATR_nd -> MR_st, D_st, IL_st */
	nnodes++;
      }

      else if (ct[i] == j) { /* i,j paired to each other. MATP. deal with INSL, INSR */
	v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, MATP_nd);
	i++;
	j--;
	PushNstack(pda, v);
	PushNstack(pda, i);
	PushNstack(pda, j);
	PushNstack(pda, DUMMY_nd); /* we don't know yet what the next node will be */
	nstates += 6;		/* MATP_nd -> MP_st, ML_st, MR_st, D_st, IL_st, IR_st */
	nnodes++;
      }

      else /* i,j paired but not to each other. BIFURC. no INS. */
	{
	  /* Here's the first of two places where we can optimize the topology 
           * of a CM. (The other comes from choosing a state traversal order when 
           * building the CM.) Imagine a multifurcation of four domains:
           *   [1]..[2]..[3]..[4]
           * The "default leftwise" rule means that BEGL/INSL generates the 
           * intervening sequences, so we must model the four stems as:
           *   [1],    ..[2],    ..[3],   ..[4]
           * but we have a choice of how we bifurcate:
           *   (1,2)(3,4)    1,(2,(3,4))   ((1,2),3),4 
           * Our choice affects the time and memory requirements of a divide
           * conquer alignment algorithm. (1,(2,(3,4)) is most efficient for
           * memory; (1,2)(3,4) is most efficient for time.
           * 
           * So we may want to choose carefully from several possible split 
           * points k (3, in the above example). A priori we only know one 
           * possible midpoint precisely: ct[i]+1, the next base after closing 
           * domain 1. We can find the others by scanning for them, and we can 
           * be reasonably efficient about scanning by using ct[] to instantly 
           * skip subdomains.
	   */
	  /* One possible rule: optimize by finding most balanced split.
           * Each stop of the following loop gives a possible midpoint k, which is
           * then evaluated, keeping track of the best split so far.
           */
	  v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, BIF_nd);

	  bestk    = ct[i]+1;
	  bestdiff = clen;
	  for (k = ct[i] + 1; k < ct[j]; k = ct[k] + 1) 
	    {
	      diff = abs(i+j-2*k); /* = len2-len1-1, where len2 = j-k+1, len1= k-i */
	      if (diff < bestdiff) {
		bestdiff = diff; 
		bestk    = k;
	      }
	      while (ct[k] == 0) k++;
	    }
				/* push the right BEGIN node first */
	  PushNstack(pda, v);	
	  PushNstack(pda, bestk);
	  PushNstack(pda, j);
	  PushNstack(pda, BEGR_nd);
				/* then push the left BEGIN node */
	  PushNstack(pda, v);	
	  PushNstack(pda, i);
	  PushNstack(pda, bestk-1);
	  PushNstack(pda, BEGL_nd);
	  nstates += 1;		/* BIF_nd -> B_st */
	  nnodes++;
	}
    }	/* while something's on the stack */
  FreeNstack(pda);
  free(ct);

  /* OK, we've converted ct into gtr -- gtr is a tree structure telling us the
   * arrangement of consensus nodes. Now do the drill for constructing a full model 
   * using this guide tree.
   */
  cm = CreateCM(nnodes, nstates);
  cm_from_guide(cm, gtr);
  CMZero(cm);

  if (ret_cm  != NULL) *ret_cm  = cm;  else FreeCM(cm);
  if (ret_gtr != NULL) *ret_gtr = gtr; else FreeParsetree(gtr);
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
  int k; /* HMM node counter */
  int ks; /* HMM state counter (0(Match) 1(insert) or 2(delete)*/
  int n; /* CM node that maps to HMM node k */
  int nn; /* CM node index */
  int n_begr; /* CM node index */
  int is_left; /* TRUE if HMM node k maps to left half of CM node n */
  int is_right; /* TRUE if HMM node k maps to right half of CM node n */
  int v; /* state index in CM */
  int v1, v2;
  CMEmitMap_t *orig_emap;         /* consensus emit map for the original, template CM */
  CMEmitMap_t *sub_emap;         /* consensus emit map for the subCM */
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
  int found_bif;
  char **nodetypes;
  char **sttypes;
  int orig2sub_mapped_lpos;
  int orig2sub_mapped_rpos;
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
      printf("!cc: %d lins: %d rins: %d\n", cc, orig_cc2lins_map[cc], orig_cc2rins_map[cc]);
    }
  
  for(cc = sub_start; cc <= sub_end; cc++)
    {
      n_o = o_cc_node_map[cc];
      printf("n_o: %d\n", n_o);
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
  printf("\n\n\n");
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
	  printf("                          | orig v:%4d(%4d) %6s%6s\n", v_o, n_o, nodetypes[orig_cm->ndtype[n_o]], sttypes[orig_cm->sttype[v_o]]);
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
  int ks_o; /* HMM state counter (0(Match) 1(insert) or 2(delete)*/
  int v_o;
  int v_s;
  int n; /* CM node that maps to HMM node k */
  int nn; /* CM node index */
  int n_begr; /* CM node index */
  int is_left; /* TRUE if HMM node k maps to left half of CM node n */
  int is_right; /* TRUE if HMM node k maps to right half of CM node n */
  int v; /* state index in CM */
  int v1, v2;
  CMEmitMap_t *orig_emap;         /* consensus emit map for the original, template CM */
  CMEmitMap_t *sub_emap;         /* consensus emit map for the subCM */
  int cc;
  int n_o;
  int n_s;
  int found_bif;
  char **nodetypes;
  char **sttypes;
  int orig2sub_mapped_lpos;
  int orig2sub_mapped_rpos;
  int *orig_cc2lins_map;
  int *orig_cc2rins_map;
  int **orig2sub_smap;
  int **sub2orig_smap;
  int cp9_M;
  int sub_k;
  int orig_k;
  int sub_nd;
  int orig_nd;

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
  /*
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
  */
  for(sub_k = 1; sub_k <= cp9hmm_s->M; sub_k++)
    {
      orig_k = sub_k + sub_start - 1;
      
      printf("sub_k: %d | orig_k: %d\n", sub_k, orig_k);
      sub_nd  = cc_node_map_s[sub_k];
      orig_nd = cc_node_map_o[orig_k];
      printf("sub_nd: %d (%s) orig_nd: %d (%s)\n", sub_nd, nodetypes[sub_cm->ndtype[sub_nd]],
	     orig_nd, nodetypes[orig_cm->ndtype[orig_nd]]);

      printf("sub_k: %d | orig_k: %d\n", sub_k, orig_k);
      /* We have to be careful about sub_cm MP nodes:
       * they necessarily map to a orig_cm MP node, and
       * we only want to map MATP_MP to MATP_MP, MATP_ML
       * to MATP_ML, MATP_MR to MATP_MR here and MATP_D
       * to MATP_D, so we have to explicitly check here,
       * otherwise we'd have problems: for example we'd
       * try to map MATP_MP to MATP_ML (because they 
       * both map to the same consensus column).
       * This is checked for in map_orig2sub_cm_helper()
       */

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

      /* We have to be careful here, as we 
       * don't want to map other inserts to the MATP_IL and
       * MATP_IR, even if they emit in the same
       * place, we want the MATP_I* of the sub_cm
       * to exactly mirror the MATP_I* of the orig_cm,
       * ambiguous as it is. We check for this in map_orig2sub_cm_helper().
       */      
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

  /* NOTE: We ignore mapping B states, E states and S states (except ROOT_S). We'll handle these guys 
   * specially when we fill in transition probabilities. The reason we don't map them is that
   * I think its impossible to do it robustly. I saw cases with SSU where there were BIF nodes in the
   * sub_cm for which I couldn't figure out which BIF nodes (and correspondingly BEGL and BEGR nodes)
   * in the orig_cm they should map to. Regardless, even if there is a pretty, nice way I'm abandoning
   * an attempt to find it for a crude way - we will handle the transitions involving these guys
   * specially, without a need for a mapping between CMs.
   */

  /* print sub2orig_smap, checking consistency with orig2sub_smap along the way. */
  printf("\n\n\n");
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
	  printf("                             | orig v:%4d(%4d) %6s%6s\n", v_o, n_o, nodetypes[orig_cm->ndtype[n_o]], sttypes[orig_cm->sttype[v_o]]);
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
 * Returns: (orig_vid) 
 */
static void
map_orig2sub_cm_helper(CM_t *orig_cm, CM_t *sub_cm, int **orig2sub_smap, int **sub2orig_smap, int orig_v, int sub_v)
{
  int sub_nd;
  int orig_nd;
  int is_insert;
  int temp_sub_v;
  int temp;

  printf("\nin helper: orig_v: %d sub_v: %d\n", orig_v, sub_v);

  if(orig_v == -1 || sub_v == -1)
    return;

  orig_nd = orig_cm->ndidx[orig_v];
  sub_nd  =  sub_cm->ndidx[sub_v];

  is_insert = FALSE;
  if(sub_cm->sttype[sub_v] == IL_st || sub_cm->sttype[sub_v] == IR_st)
    {
      if(!(orig_cm->sttype[orig_v] == IL_st || orig_cm->sttype[orig_v] == IR_st))
	Die("INSERT error in map_orig2sub_cm_helper()\n");
      is_insert = TRUE;
    }
      
  /* check for cases we want to avoid mapping, we exploit the fact that we know a MATP_nd in a sub_cm MUST be mapped 
   * to a corresponding MATP_nd in the original template CM: */
  if(sub_cm->ndtype[sub_nd] == MATP_nd)
    {
      if(orig_cm->ndtype[orig_nd] == MATP_nd && sub_cm->sttype[sub_v] != orig_cm->sttype[orig_v])
	{
	  printf("return case 1\n");
	  return;
	}
      if(orig_cm->ndtype[orig_nd] != MATP_nd)
	{
	  if(!(is_insert))
	    {
	      printf("return case 2\n");
	      return;
	    }
	  else
	    {
	      /* we have a sub MATP insert we're mapping to an orig non-MATP insert, we don't want to do this
	       * but first we want to check to make sure we don't already have a mapping for these two guys,
	       * and if we do, obliterate it.
	       */
	      printf("in else\n");
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
	      printf("return case 3\n");
	      return;
	    }
	}
    }
  /* check for case that we're mapping a orig MATP_I* to a sub MATP_I*, but we already have a mapping for the
   * sub MATP_I*, in which case we obliterate it. 
   */
  if(is_insert && (orig_cm->ndtype[orig_nd] == MATP_nd && sub_cm->ndtype[sub_nd] == MATP_nd))
    {
      printf("kachow\n");

      if(orig2sub_smap[orig_v][1] != -1)
	Die("Obliteration error 0\n");
      
      printf("orig2sub_smap[orig_v:%d][0]: %d\n", orig_v, orig2sub_smap[orig_v][0]);
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
	      printf("\n\nobliteration error\n");
	      temp = orig2sub_smap[orig_v][0];
	      printf("\twe thought either sub2orig_smap[%d][0]:%d == sub_v:%d\n\tor sub2orig_smap[%d][1]:%d == sub_v: %d\n", temp, sub2orig_smap[temp][0], sub_v, temp, sub2orig_smap[temp][1], sub_v);
	    }
	}
    }

  /* check for case that we've already mapped this orig_cm MATP_I* state to a sub_cm MATP_I* state,
   * in which case we want to skip the current mapping.
   */
  if(is_insert && (orig_cm->ndtype[orig_nd] == MATP_nd && sub_cm->ndtype[sub_nd] != MATP_nd))
    {
      if(orig2sub_smap[orig_v][0] != -1)
	{
	  if((orig2sub_smap[orig_v][1] != -1 || (sub2orig_smap[sub_v][1] != -1)) || (sub2orig_smap[sub_v][0] == -1))
	    Die("ERROR with sub non-MATP to orig MATP insert 2\n");
	  if(sub_cm->ndtype[sub_cm->ndidx[orig2sub_smap[orig_v][0]]] == MATP_nd)
	    {
	      printf("return case 4\n");
	      return;
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
	  printf("orig2sub filling 0: orig2sub_smap[%d][0]: %d\n", orig_v, orig2sub_smap[orig_v][0]);
	}
    }
  else if (orig2sub_smap[orig_v][1] != -1)
    {
      if(orig2sub_smap[orig_v][0] == sub_v || orig2sub_smap[orig_v][1] == sub_v)
	/* abort!, we already have this mapping */
	return; 
      else
	Die("ERROR in map_orig2sub_cm_helper, orig2sub_smap[%d][0] is not -1 and orig2sub_smap[%d][1] is not -1, this shouldn't happen.\n", orig_v, orig_v);
    }
  else /* orig2sub_smap[orig_v][0] != -1 && orig2sub_smap[orig_v][1] == -1 */
    {
      if(orig2sub_smap[orig_v][0] == sub_v || orig2sub_smap[orig_v][1] == sub_v)
	/* abort!, we already have this mapping */
	return; 
      orig2sub_smap[orig_v][1] = sub_v;
      printf("orig2sub filling 1: orig2sub_smap[%d][1]: %d\n", orig_v, orig2sub_smap[orig_v][1]);
    }

  /* now fill in sub2orig_smap */
  printf("sub2orig_smap[%d][0]: %d | sub2orig_smap[%d][1]: %d\n", sub_v,  sub2orig_smap[sub_v][0],  sub_v,  sub2orig_smap[sub_v][1]);
  if(sub2orig_smap[sub_v][0] == -1)
    {
      if (sub2orig_smap[sub_v][1] != -1)
	Die("ERROR in map_sub2orig_cm_helper, sub2orig_smap[%d][0] is -1 but sub2orig_smap[%d][1] is not, this shouldn't happen.\n", sub_v, sub_v);
      else
	{
	  sub2orig_smap[sub_v][0] = orig_v;
	  printf("sub2orig filling 0: sub2orig_smap[%d][0]: %d\n", sub_v, sub2orig_smap[sub_v][0]);
	}
    }
  else if (sub2orig_smap[sub_v][1] != -1)
    Die("ERROR in map_sub2orig_cm_helper, sub2orig_smap[%d][0] is not -1 and sub2orig_smap[%d][1] is not -1, this shouldn't happen.\n", sub_v, sub_v);
  else /* sub2orig_smap[sub_v][0] != -1 && sub2orig_smap[sub_v][1] == -1 */
    {
      sub2orig_smap[sub_v][1] = orig_v;
      printf("sub2orig filling 1: sub2orig_smap[%d][1]: %d\n", sub_v, sub2orig_smap[sub_v][1]);
    }
  return;
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
cm2sub_cm_emit_probs(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, int v_s, int v_o1, int v_o2)
{
  int is_left;
  int i, j;

  printf("\nin cm2sub_cm_emit_probs v_s: %d, v_o1: %d, v_o2: %d\n", v_s, v_o1, v_o2);

  if(v_o1 == -1)
    Die("ERROR in cm2sub_cm_emit_probs, sub_cm state %d maps to 0 states in orig_cm\n", v_s);

  /* We handle MP's special */
  if(sub_cm->sttype[v_s] == MP_st)
    {
      printf("MATP case\n");
      /* sanity check */
      if(v_o2 != -1 || orig_cm->sttype[v_o1] != MP_st)
	Die("ERROR filling emission parameters for sub_cm state: %d\n", v_s);
      for(i = 0; i < (MAXABET*MAXABET); i++)
	sub_cm->e[v_s][i] = orig_cm->e[v_o1][i];
      /* The FNorm shouldn't be necessary, and since we build a new CM for each sequence in --sub mode, 
	 we skip it */
      /* FNorm(sub_cm->e[v_s], MAXABET*MAXABET)*/
      return;
    }

  /* if we get here, v_s is a singlet emitter */
  
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
	  printf("v_o1 MATP, v_s LEFT case\n");
	}
      else if(orig_cm->sttype[v_o2] == MR_st)
	{
	  is_left = FALSE;
	  printf("v_o1 MATP, v_s RIGHT case\n");
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
      printf("v_o2 ML, MR, IL or IR\n");
      for(i = 0; i < MAXABET; i++)
	sub_cm->e[v_s][i] += orig_psi[v_o2] * orig_cm->e[v_o2][i];
      FNorm(sub_cm->e[v_s], MAXABET);
      return;
    }
  else if(v_o2 != -1)
    {
      printf("dual insert case\n");
      if(!((orig_cm->sttype[v_o1] == IL_st && orig_cm->sttype[v_o2] == IR_st) ||
	   (orig_cm->sttype[v_o1] == IR_st && orig_cm->sttype[v_o2] == IL_st)))
	Die("ERROR sub_cm state: %d maps to two states, but not case 1: one is MATP_MP OR case 2: one is IL and one is IR\n", v_o1, v_o2);
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
  printf("singlet emitter case\n");
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
 *           First, fill the sub_cm->t[v_s] with 'virtual
 *           counts', then normalize to probabilities.
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
 * int v_s           - the sub_cm state we're filling emissions for
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
		      int **sub2orig_smap)
{
  int v_o;
  int yoffset;
  int y_s;
  int y_o;

  printf("in cm2sub_cm_trans_probs: v_s: %d\n", v_s);

  /* start with the first orig_cm state that maps to v_s */
  v_o = sub2orig_smap[v_s][0];
  printf("\tv_o: %d\n", v_o);
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
	  printf("\ty_s: %d | y_o: %d\n", y_s, y_o);
	  
	  if(y_o == -1)
	    {
	      if(sub_cm->sttype[y_s] != B_st &&
		 sub_cm->sttype[y_s] != E_st)
		Die("ERROR, sub_cm state y_s: %d (v_s: %d) maps to no state in sub_cm, but it's not a B or E state\n", y_s, v_s);
	      /* We're transitioning to a BIF_B or END_E; not sure how to handle. */
	      /* TEMPORARY SOLUTION: */
	      sub_cm->t[v_s][yoffset] += 1.0;
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
	      /* We're transitioning to a BIF_B or END_E; not sure how to handle. */
	      /* TEMPORARY SOLUTION: */
	      sub_cm->t[v_s][yoffset] += 1.0;
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
 * EPN 08.31.06
 * Function: cm_add_single_trans()
 *
 * Purpose:  Add a virtual counts contribution to a single CM transition.
 * 
 * See related functions for explanation of parameters. 
 * Returns: (void) 
 */
static void
cm_add_single_trans(CM_t *orig_cm, CM_t *sub_cm, int **orig2sub_smap, int **sub2orig_smap, int v_o, int y_o, 
		    int v_s, int yoffset, double *orig_psi, char ***tmap)
{
  int y_s;
  y_s = sub_cm->cfirst[v_s] + yoffset;
  printf("\t\tin cm_add_single_trans, v_o: %d | y_o: %d | v_s: %d | y_s: %d\n", v_o, y_o, v_s, y_s);
  printf("\t\tbeg sub_cm->t[v_s:%d][yoffset:%d] : %.9f\n", v_s, yoffset, sub_cm->t[v_s][yoffset]);

  /* check if we've got real CM state ids */
  if(v_o == -1 || y_o == -1)
    return;

  int start;
  int end;
  int orig_insert1, orig_insert2, orig_insert3, orig_insert4;
  int sub_insert1, sub_insert2;
  int n_s;
  

  if(v_o <= y_o) /* going DOWN the CM */
    {
      start = v_o;
      end = y_o;
    }
  else /* going UP the CM */
    {
      start = y_o;
      end = v_o;
    }
  
  /* find the insert states in orig_cm that map to the sub_cm insert states in the same node as v_s,
   * cm_sum_subpaths needs to know these states */
  orig_insert1 = orig_insert2 = orig_insert3 = orig_insert4 = -1;
  sub_insert1  = sub_insert2 = -1;
  n_s          = sub_cm->ndidx[v_s];
  /* find the inserts in n_s, we know each state in a node has all the inserts in that node
   * as children (if there are any), and further these must be the top 2 lowest numbered children.
   * We're exploiting our knowledge of the CM architecture here -- not good. */
  if(sub_cm->sttype[v_s] != B_st)
    if(sub_cm->sttype[sub_cm->cfirst[v_s]] == IL_st || 
       sub_cm->sttype[sub_cm->cfirst[v_s]] == IR_st)
      {
	sub_insert1 = sub_cm->cfirst[v_s];
	if(sub_cm->sttype[sub_cm->cfirst[v_s]+1] == IR_st)
	  sub_insert2 = sub_cm->cfirst[v_s] + 1;
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

  printf("\t\t\toi1: %d oi2: %d oi3: %d oi4: %d\n", orig_insert1, orig_insert2, orig_insert3, orig_insert4);
  
  sub_cm->t[v_s][yoffset] += orig_psi[start] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, start, end, v_s, tmap, orig_psi, orig_insert1, orig_insert2, orig_insert3, orig_insert4);

  printf("\t\tend sub_cm->t[v_s:%d][yoffset:%d] : %.9f\n", v_s, yoffset, sub_cm->t[v_s][yoffset]);
  
}

/**************************************************************************
 * EPN 08.31.06
 * Function: cm_sum_subpaths()
 *
 * Purpose:  Calculated probability of getting from one state (start) to 
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
 * int sub_start     - state index in start we're starting at
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
  int s_n;
  int e_n;
  int sub_end1; /* up to 2 sub_cm states that contain a state mapping */
  int sub_end2; /* to orig_cm state 'end'; -1 if none.*/
  int sub_v1;   /* up to 2 sub_cm states that contain a state mapping */
  int sub_v2;   /* to orig_cm state v; -1 if none.*/
  int skip;     /* TRUE to skip state v's contribution */

  int temp;

  printf("\t\t\tin cm_sum_subpaths, start: %d | end: %d\n", start, end);
  
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
      printf("\t\t\tReturning self insert prob: %f\n", orig_cm->t[start][0]);
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
      printf("\t\tinsert_to_start: %f sub_psi[0]: %f\n", insert_to_start, sub_psi[0]);
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
      
      printf("\t\t\tv: %d, se1: %d se2: %d sv1: %d sv2: %d\n", v, sub_end1, sub_end2, sub_v1, sub_v2);

      /* We want to check if there is a sub_cm transition from the sub_cm state that 
       * maps to v (v_s) to end exists, if so then we want to skip the contribution here.
       * Unless v = end.
       */
      /*if(v != end && is_insert)*/
      skip = FALSE;
      if(v != end)
	{
	  if((is_insert) && ((sub_trans_check(sub_cm, orig2sub_smap, v, end, sub_start)) ||
			     (sub_trans_check(sub_cm, orig2sub_smap, v, start, sub_start))))
	    {
	      /* if v is an insert state and there exists a transition from 
	       * sub_v -> sub_end or sub_end -> sub_v (where sub_v is a sub_cm state
	       * that maps to orig_cm state v, and sub_end is a sub_cm state that 
	       * maps to orig_cm state end), we skip the contribution because
	       * it will counted when compiling counts for sub_v -> sub_end or 
	       * sub_end -> sub_v.
	       */
	      skip = TRUE;
	      printf("\t\t\tskip true case 911\n");
	    }

	  if((is_insert) && (orig2sub_state_check(orig2sub_smap, v, end)))
	    {
	      /* if v is an insert state that maps to the same orig_cm state that end 
	       * maps to, then we don't want to double count its contribution (it will be 
	       * counted in a subsequent cm_sum_subpaths() call), so we skip it here.
	       */
	      /*skip = TRUE;*/
	      printf("\t\t\tskip true case CP9\n");
	    }

	  if(sub_v1 > sub_start) /* otherwise we don't want to skip */
	    {
	      /* check for possible sub_cm transitions */
	      if(sub_v1 != -1 && sub_cm->sttype[sub_v1] != B_st)
		{
		  if(sub_end1 != -1 && (sub_end1 > sub_cm->cfirst[sub_v1] && 
					((sub_end1 - sub_cm->cfirst[sub_v1]) < sub_cm->cnum[sub_v1])))
		    {
		      /*skip = TRUE;*/ /* a transition from sub_v1 -> sub_end1 exists in sub_cm */
		      printf("\t\t\t\tskip true case 1\n");
		    }
		  if(sub_end1 != -1 && (sub_end2 > sub_cm->cfirst[sub_v1] && 
					((sub_end2 - sub_cm->cfirst[sub_v1]) < sub_cm->cnum[sub_v1])))
		    {
		      /*skip = TRUE;*/ /* a transition from sub_v1 -> sub_end2 exists in sub_cm */
		      printf("\t\t\t\tskip true case 2\n");
		    }
		}
	    }
	  
	  if(sub_v2 > sub_start) /* otherwise we don't want to skip */
	    {
	      if(sub_v2 != -1 && sub_cm->sttype[sub_v2] != B_st)
		{
		  if(sub_end1 != -1 && (sub_end1 > sub_cm->cfirst[sub_v2] && 
					((sub_end1 - sub_cm->cfirst[sub_v1]) < sub_cm->cnum[sub_v2])))
		    {
		      /*skip = TRUE;*/ /* a transition from sub_v2 -> sub_end1 exists in sub_cm */
		      printf("\t\t\t\tskip true case 2\n");
		    }
		  
		  if(sub_end2 != -1 && (sub_end2 > sub_cm->cfirst[sub_v2] 
					&& ((sub_end2 - sub_cm->cfirst[sub_v2]) < sub_cm->cnum[sub_v2])))
		    {
		      /*skip = TRUE;*/ /* a transition from sub_v1 -> sub_end2 exists in sub_cm */
			  printf("\t\t\t\tskip true case 2\n");
		    }
		}
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
	  printf("v: %d | skipping the contribution\n", v);
	}
    }
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
 * Function: check_subCM_by_sampling()
 *
 * Purpose:  Given a CM and a sub CM that is supposed to mirror 
 *           the CM as closely as possible between two given consensus
 *           columns (spos and epos), check that the subCM was correctly 
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
 *
 * Returns: TRUE: if CM and sub CM are "close enough" (see code)
 *          FALSE: otherwise
 */
int 
check_subCM_by_sampling(CM_t *orig_cm, CM_t *sub_cm, int spos, int epos)
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

  debug_level = 0;
  ret_val = TRUE;

  /* Build a CP9 HMM from the sub_cm */
  sub_hmm = AllocCPlan9((epos-spos+1));
  ZeroCPlan9(sub_hmm);
  CPlan9SetNullModel(sub_hmm, sub_cm->null, 1.0); /* set p1 = 1.0 which corresponds to the CM */

  map_consensus_columns(sub_cm, sub_hmm->M, &node_cc_left, &node_cc_right,
			&cc_node_map, debug_level);
  
  printf("sub_hmm->M: %d\n", sub_hmm->M);
  CP9_map_cm2hmm_and_hmm2cm(sub_cm, sub_hmm, node_cc_left, node_cc_right, 
			    cc_node_map, &cs2hn_map, &cs2hs_map, 
			    &hns2cs_map, debug_level);
  
  if(!(CP9_cm2wrhmm(sub_cm, sub_hmm, node_cc_left, node_cc_right, cc_node_map, cs2hn_map,
		    cs2hs_map, hns2cs_map, debug_level)))
    Die("Couldn't build a CM Plan 9 HMM from the CM.\n");
  
  if(!(CP9_check_wrhmm_by_sampling(orig_cm, sub_hmm, spos, epos, hns2cs_map, 0.05, 100000)))
    Die("CM Plan 9 fails sampling check!\n");
  else
    printf("CM Plan 9 passed sampling check.\n");

}

/**************************************************************************
 * EPN 09.06.06
 * Function: check_subCM_by_sampling2()
 *
 * Purpose:  Given a CM and a sub CM that is supposed to mirror 
 *           the CM as closely as possible between two given consensus
 *           columns (spos and epos), check that the subCM was correctly 
 *           constructed. 
 *           
 *           The approach is to sample from the CM and the subCM 
 *           and use those samples to build two CP9 HMMs, then 
 *           compare those two CP9 HMMs.

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
check_subCM_by_sampling2(CM_t *orig_cm, CM_t *sub_cm, int spos, int epos, int nseq)
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
  int i, idx, nd;
  int L;
  int apos;
  int *matassign;
  int *useme;
  struct cp9trace_s **cp9_tr;   /* fake tracebacks for each seq            */
  struct cplan9_s  *shmm;       /* the new, CM plan9 HMM; built by sampling*/
  int msa_nseq;                 /* this is the number of sequences per MSA,
				 * current strategy is to sample (nseq/nseq_per_msa)
				 * alignments from the CM, and add counts from
				 * each to the shmm in counts form (to limit memory)
				 */
  int nsampled;                 /* number of sequences sampled thus far */
  int debug_level;
  int cc, cp9_M, v;

  printf("\n\n*****************\nin check_subCM_by_sampling2()\n\n");

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
  int k;
  int y;
  double temp_psi;
  int violation;
  int v_ct; /* Number of violations */
  double diff;
  int ret_val; /* return value */

  if(print_flag)
    {
      printf("Printing orig_psi:\n");
      
      for(v_o = 0; v_o < orig_cm->M; v_o++)
	printf("orig_psi[%4d]: %.4f\n", v_o, orig_psi[v_o]);
      
    }

  ret_val = TRUE;
  v_ct = 0;

  if(print_flag == TRUE)
    printf("\n");
  for(v_s = 0; v_s < sub_cm->M; v_s++)
    {
      if(print_flag)
	printf("\tv_s: %4d (%.4f) ", v_s, sub_psi[v_s]);
      
      v_o = sub2orig_smap[v_s][0];
      if(v_o != -1)
	{
	  if(print_flag)
	    printf("v_o1: %4d (%.4f) ", v_o, orig_psi[v_o]);
	  
	  temp_psi = orig_psi[v_o];
	  v_o = sub2orig_smap[v_s][1];
	  if(v_o != -1)
	    {
	      temp_psi += orig_psi[v_o];
	      if(print_flag)
		printf("v_o2: %4d (%.4f)\n", v_o, orig_psi[v_o]);
	    }
	  else
	    if(print_flag) printf("\n");
	  
	  printf("sub: %.4f | orig: %.4f | diff: %.4f\n\n", sub_psi[v_s], temp_psi, (sub_psi[v_s]-temp_psi));
	  
	}
      else
	{
	  if(sub_cm->sttype[v_s] != E_st &&
	     sub_cm->sttype[v_s] != B_st &&
	     sub_cm->sttype[v_s] != S_st &&
	     sub_cm->sttype[v_s] != EL_st)
	    Die("ERROR state v_s:%d maps to nothing and its not E,B,S,EL\n", v_s);
	  printf("E B S or EL\n");
	}
    }
  

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
  
  printf("in sub_trans_check() orig_a: %d orig_b: %d min_sub_v: %d\n", orig_a, orig_b, min_sub_v);

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
static int
trans_check_helper(CM_t *cm, int a, int b)
{
  printf("\t**in trans_check_helper a: %d | b: %d\n");

  if((a == -1 || b == -1) || (a > b))
    return -1;

  if((b - cm->cfirst[a]) < cm->cnum[a])
    { printf("returning a: %d\n", a); return a; }

  return -1;

}


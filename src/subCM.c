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
#include "hmmer_funcs.h"

static void map_orig2sub_cm(CM_t *orig_cm, CM_t *sub_cm, int ***ret_orig2sub_smap, int ***ret_sub2orig_smap,
			    int sub_start, int sub_end);
static void map_orig2sub_cm_helper(int **orig2sub_smap, int **sub2orig_smap, int v_o, int v_s);
static void cm2sub_cm_emit_probs(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, int v_s, int v_o1, int v_o2);
static void cm2sub_cm_trans_probs(CM_t *orig_cm, CM_t *sub_cm, double *orig_psi, char ***tmap, int v_s, 
				  int **orig2sub_smap, int **sub2orig_smap);
static void cm_add_single_trans(CM_t *orig_cm, CM_t *sub_cm, int **orig2sub_smap, int v_o, int y_o,
				int v_s, int yoffset, double *orig_psi, char ***tmap);
static float cm_sum_subpaths(CM_t *orig_cm, CM_t *sub_cm, int **orig2sub_smap, int start, int start_s, int end,
			     char ***tmap, double *psi);
static void debug_print_cm_params(CM_t *cm);
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


/* Function:  StripWUSS()
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
  int sub_cpos;
  int v_s;
  int v_o1;
  int v_o2;

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
      printf("B sub_ct[cpos=%3d]: %3d\n", sub_cpos, sub_ct[sub_cpos]);
      if ((cpos+1) < struct_start || (cpos+1) > struct_end) /* cpos goes 1..clen, but ct is indexed
							     * 0..clen-1.*/
	{ 

	  if (sub_ct[sub_cpos] != -1)  sub_ct[sub_ct[sub_cpos]] = -1; /* CreateCMConsensus() uses
								       * -1 in ct[] to indicate single
								       * stranded (different convention
								       * than WUSS2ct(). */
	  sub_ct[sub_cpos] = -1;
	}
      printf("A sub_ct[cpos=%3d]: %3d\n\n", sub_cpos, sub_ct[sub_cpos]);
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
   * the new, subCM.
   */
  orig_emap = CreateEmitMap(orig_cm);
  sub_emap = CreateEmitMap(sub_cm);
  printf("\n\n\nDumping original CM emitmap\n");
  DumpEmitMap(stdout, orig_emap, orig_cm);
  printf("\n\n\nDumping subCM emitmap\n");
  DumpEmitMap(stdout, sub_emap, sub_cm);

  /* Map states from orig_cm to sub_cm and vice versa. */
  map_orig2sub_cm(orig_cm, sub_cm, &orig2sub_smap, &sub2orig_smap, model_start,
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
      if(sub_cm->sttype[v_s] != S_st &&
	 sub_cm->sttype[v_s] != B_st &&
	 sub_cm->sttype[v_s] != E_st)
	cm2sub_cm_trans_probs(orig_cm, sub_cm, orig_psi, tmap, v_s, orig2sub_smap, sub2orig_smap);
    }

  /* Finally renormalize the CM */
  CMRenormalize(sub_cm);

  printf("\nDEBUG PRINT OF ORIG_CM PARAMETERS:\n");
  debug_print_cm_params(orig_cm);
  printf("\nDEBUG PRINT OF SUB_CM PARAMETERS:\n");
  debug_print_cm_params(sub_cm);

  free(sub_cstr);
  free(sub_ct);
  FreeEmitMap(orig_emap);
  FreeEmitMap(sub_emap);
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
 * int sub_end         - lasst consensus column in orig_cm that sub_cm models
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
  o_node_cc_left  = malloc(sizeof(int) * orig_cm->nodes);
  o_node_cc_right = malloc(sizeof(int) * orig_cm->nodes);
  o_cc_node_map   = malloc(sizeof(int) * (orig_emap->clen + 1));
  map_consensus_columns(orig_cm, orig_emap->clen, o_node_cc_left, o_node_cc_right,
			o_cc_node_map, 0);
  s_node_cc_left  = malloc(sizeof(int) * sub_cm->nodes);
  s_node_cc_right = malloc(sizeof(int) * sub_cm->nodes);
  s_cc_node_map   = malloc(sizeof(int) * (sub_emap->clen + 1));
  map_consensus_columns(sub_cm,  sub_emap->clen,  s_node_cc_left, s_node_cc_right,
			s_cc_node_map, 0);

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
  orig_cc2lins_map = MallocOrDie(sizeof(int) * orig_emap->clen +2);
  orig_cc2rins_map = MallocOrDie(sizeof(int) * orig_emap->clen +2);
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
      n_s = s_cc_node_map[(cc-sub_start+1)]; /* note offset */
      if(orig_cm->ndtype[n_o] != MATP_nd &&
	 orig_cm->ndtype[n_o] != MATL_nd &&
	 orig_cm->ndtype[n_o] != MATR_nd)
	  Die("ERROR 1o in map_orig2sub_cm()\n");
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
	      printf("mapping v_s:%d to v_o:%d\n", v_s, v_o);
	      if(sub_cm->ndidx[v_s] != n_s)
		Die("ERROR 2 in map_orig2sub_cm()\n");
	      map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);
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
	  v_o = orig_cm->nodemap[n_o];
	  v_s =  sub_cm->nodemap[n_s];
	  if(sub_cm->ndtype[n_s] != MATL_nd && sub_cm->ndtype[n_s] != MATR_nd)
	    {
	      printf("orig_cm %4d %4d %4s %2s | sub_cm %4d %4d %4s %2s\n", n_o, v_o, nodetypes[orig_cm->ndtype[n_o]], sttypes[orig_cm->sttype[v_o]], n_s, v_s, nodetypes[sub_cm->ndtype[n_s]], sttypes[sub_cm->sttype[v_s]]);
	      Die("ERROR 3 in map_orig2sub_cm()\n");
	    }
	  
	  if(sub_cm->ndtype[n_s] == MATL_nd) /* case 1 */
	    {
	      /* Case by case (uchh..) 
	       * v_o is MATP_MP 
	       * v_s is MATL_ML */
	      map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);
	      v_o++; /* v_o is MATP_ML */
	      map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);

	      v_s++; /* v_s is MATL_D */
	      v_o++; /* v_o is MATP_MR */
	      map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);
	      v_o++; /* v_o is MATP_D */
	      map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);
	      
	      v_s++; /* v_s is MATL_IL */
	      v_o++; /* v_o is MATP_IL */
	      map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);
	    }
	  else if(sub_cm->ndtype[n_s] == MATR_nd) /* case 2 */
	    {
	      /* Case by case (uchh..)
	      /* v_o is MATP_MP */
	      /* v_s is MATR_MR */
	      if(orig_cm->sttype[v_o] != MP_st)
		Die("WHOA, state v_o: %d is not a MATP_MP\n");
	      if(sub_cm->sttype[v_s] != MR_st)
		Die("WHOA, state v_s: %d is not a MATR_MR\n");

	      printf("mapping MATP_MP orig v_o: %d to MATR_R sub v_s: %d\n", v_o, v_s);
	      map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);
	      v_o += 2; /* v_o is MATP_MR */
	      if(orig_cm->sttype[v_o] != MR_st)
		Die("WHOA, state v_o: %d is not a MATP_MR\n");
	      map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);

	      v_s++; /* v_s is MATR_D */
	      if(sub_cm->sttype[v_s] != D_st)
		Die("WHOA, state v_s: %d is not a MATR_D\n");
	      v_o--; /* v_o is MATP_ML */
	      if(orig_cm->sttype[v_o] != ML_st)
		Die("WHOA, state v_o: %d is not a MATP_ML\n");
	      map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);
	      v_o += 2; /* v_o is MATP_D */
	      if(orig_cm->sttype[v_o] != D_st)
		Die("WHOA, state v_o: %d is not a MATP_D\n");
	      map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);
	      
	      v_s++; /* v_s is MATR_IR */
	      if(sub_cm->sttype[v_s] != IR_st)
		Die("WHOA, state v_s: %d is not a MATR_IR\n");
	      v_o += 2; /* v_o is MATP_IR */
	      if(orig_cm->sttype[v_o] != IR_st)
		Die("WHOA, state v_o: %d is not a MATP_IR\n");
	      map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);
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
	  map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);

	  v_o++; /* v_o is MATL_D */
	  v_s++; /* v_s is MATR_D */
	  map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);

	  v_o = orig_cc2lins_map[cc-1]; /* orig_cm left insert state (if any) that inserts BEFORE cc */
	  v_s++; /* v_s is MATR_IR, potentially two orig_cm states map to v_s */
	  if(v_o != -1)
	    map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);
	  v_o = orig_cc2rins_map[cc]; /* orig_cm right insert state (if any) that inserts BEFORE cc */
	  if(v_o != -1)
	    map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);

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
	  map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);

	  v_o++; /* v_o is MATR_D */
	  v_s++; /* v_s is MATL_D */
	  map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);

	  v_o = orig_cc2lins_map[cc]; /* orig_cm left insert state (if any) that inserts AFTER cc */
	  v_s++; /* v_s is MATL_IL, potentially two orig_cm states map to v_s */
	  if(v_o != -1)
	    map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);
	  v_o = orig_cc2rins_map[cc+1]; /* orig_cm right insert state (if any) that inserts AFTER cc */
	  if(v_o != -1)
	    map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);

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
  map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);

  /* Map ROOT_IL, ROOT_IR, and BEGR_IL's in an ugly manner, by looking for the
   * insert states in the emap that map to them. 
   */

  /* Now map ROOT_IL, ROOT_IR and all BEGR_IL's, we could have either IL or IR
   * states in orig_cm that map to these guys, so we need to use
   * both the lmap and rmap. In some cases there may be two orig_cm insert states that
   * map to a single insert state in sub_cm (due to an ambiguity in the Infernal
   * grammar).

  /* ROOT_IL */
  v_s = 1;
  v_o = orig_cc2lins_map[sub_start-1]; /* we could have a left emitter that inserts before 
					* consensus column sub_start */
  if(v_o != -1)
    map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);

  v_o = orig_cc2rins_map[sub_start]; /* we could have a right emitter that inserts before 
					    * consensus column sub_start */
  if(v_o != -1)
    map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);

  /* ROOT_IR */
  v_s = 2;
  v_o = orig_cc2lins_map[sub_end]; /* we could have a left emitter that inserts after
					* consensus column sub_end */
  if(v_o != -1)
    map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);

  v_o = orig_cc2rins_map[sub_end+1]; /* we could have a right emitter that inserts after
					  * consensus column sub_end */
  if(v_o != -1)
    map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);

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
	    map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);

	  printf("checking rpos (%d + %d): %d\n", sub_emap->lpos[n_s], sub_start, (sub_emap->lpos[n_s] + sub_start));
	  v_o = orig_cc2rins_map[(sub_emap->lpos[n_s] + sub_start)];
	  printf("\t\tv_o: %d | v_s: %d\n", v_o, v_s);
	  if(v_o != -1)
	    map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);
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
	  printf("                        | orig v:%4d(%4d) %6s%6s\n", v_o, n_o, nodetypes[orig_cm->ndtype[n_o]], sttypes[orig_cm->sttype[v_o]]);
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
 * map_orig2sub_cm_helper()
 *
 * helper function for map_orig2sub_cm(), 
 *
 * Purpose:  Fill in specific part of the map, given v_o (orig_cm state),
 *           v_s (sub_cm state)
 * Args:    
 * int **orig2sub_smap
 * int **sub2orig_smap
 * int   v_o; 
 * int   v_s; 
 * Returns: (v_oid) 
 */
static void
map_orig2sub_cm_helper(int **orig2sub_smap, int **sub2orig_smap, int v_o, int v_s)
{
  /* fill in orig2sub_smap */
  if(orig2sub_smap[v_o][0] == -1)
    if (orig2sub_smap[v_o][1] != -1) 
      Die("ERROR in map_orig2sub_cm_helper, orig2sub_smap[%d][0] is -1 but orig2sub_smap[%d][1] is not, this shouldn't happen.\n", v_o, v_o);
    else
      orig2sub_smap[v_o][0] = v_s;
  else if (orig2sub_smap[v_o][1] != -1)
    Die("ERROR in map_orig2sub_cm_helper, orig2sub_smap[%d][0] is not -1 and orig2sub_smap[%d][1] is not -1, this shouldn't happen.\n", v_o, v_o);
  else if (orig2sub_smap[v_o][0] == v_s)
    {
      /* abort!, we already have this mapping */
      return; 
    }
  else /* orig2sub_smap[v_o][0] != -1 && orig2sub_smap[v_o][1] != -1 */
    orig2sub_smap[v_o][1] = v_s;

  /* now fill in sub2orig_smap */
  if(sub2orig_smap[v_s][0] == -1)
    if (sub2orig_smap[v_s][1] != -1)
      Die("ERROR in map_sub2orig_cm_helper, sub2orig_smap[%d][0] is -1 but sub2orig_smap[%d][1] is not, this shouldn't happen.\n", v_s, v_s);
    else
      sub2orig_smap[v_s][0] = v_o;
  else if (sub2orig_smap[v_s][1] != -1)
    Die("ERROR in map_sub2orig_cm_helper, sub2orig_smap[%d][0] is not -1 and sub2orig_smap[%d][1] is not -1, this shouldn't happen.\n", v_s, v_s);
  else /* sub2orig_smap[v_s][0] != -1 && sub2orig_smap[v_s][1] != -1 */
    sub2orig_smap[v_s][1] = v_o;

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

  if(v_o1 == -1)
    Die("ERROR in cm2sub_cm_emit_probs, sub_cm state %d maps to 0 states in orig_cm\n", v_s);

  /* We handle MP's special */
  if(sub_cm->sttype[v_s] == MP_st)
    {
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
  if(sub_cm->sttype[v_s] == ML_st ||
     sub_cm->sttype[v_s] == IL_st)
    is_left = TRUE;
  else 
    is_left = FALSE;
  
  /* Only way two states can map to v_s is if one of them is MP_st,
   * These are the only states we need to weight emission probs by orig_psi values,
   * and subsequently only states we need to call FNorm() for */
  if(orig_cm->sttype[v_o1] == MP_st)
    {
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
      for(i = 0; i < MAXABET; i++)
	sub_cm->e[v_s][i] += orig_psi[v_o2] * orig_cm->e[v_o2][i];
      FNorm(sub_cm->e[v_s], MAXABET);
      return;
    }
      
  /* If we get here, v_s maps to a single singlet emitter in orig_cm, v_o1 */
  if(v_o2 != -1)
    Die("ERROR sub_cm state: %d maps to two states, neither is MATP_MP\n", v_o1, v_o2);

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
	Die("ERROR, sub_cm state v_s: %d is S, E or B but maps to a orig_cm state: v_o:%d\n", v_o);

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
	    }
	  else
	    cm_add_single_trans(orig_cm, sub_cm, orig2sub_smap, v_o, y_o, v_s, yoffset, orig_psi, tmap);

	  y_o = sub2orig_smap[y_s][1];
	  cm_add_single_trans(orig_cm, sub_cm, orig2sub_smap, v_o, y_o, v_s, yoffset, orig_psi, tmap);
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
	    }
	  else
	    cm_add_single_trans(orig_cm, sub_cm, orig2sub_smap, v_o, y_o, v_s, yoffset, orig_psi, tmap);

	  y_o = sub2orig_smap[y_s][1];
	  cm_add_single_trans(orig_cm, sub_cm, orig2sub_smap, v_o, y_o, v_s, yoffset, orig_psi, tmap);
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
cm_add_single_trans(CM_t *orig_cm, CM_t *sub_cm, int **orig2sub_smap, int v_o, int y_o, int v_s, int yoffset, double *orig_psi, char ***tmap)
{
  /* check if we've got real CM state ids */
  if(v_o == -1 || y_o == -1)
    return;

  if(v_o <= y_o) /* going DOWN the CM */
    sub_cm->t[v_s][yoffset] += orig_psi[v_o] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, v_o, y_o, v_s, tmap, orig_psi);
  else if (v_o > y_o) /* going UP the CM */
    sub_cm->t[v_s][yoffset] += orig_psi[y_o] * cm_sum_subpaths(orig_cm, sub_cm, orig2sub_smap, y_o, v_o, v_s, tmap, orig_psi);
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
 * int start_s       - state index in sub_cm that start maps to
 * double *psi       - for orig_cm: psi[v] is the expected number of times state v is entered
 *                     in a CM parse
 * Returns: Float, the summed probability of all subpaths through the CM
 *          starting at "start" and ending at "end".
 */
static float
cm_sum_subpaths(CM_t *orig_cm, CM_t *sub_cm, int **orig2sub_smap, int start, int end, int start_s, char ***tmap, double *psi)
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

  /*printf("\nin cm_sum_subpaths2, start: %d | end: %d\n", start, end);*/
  if(start > end)
    {
      printf("ERROR in cm_sum_subpaths2: start: %d > end: %d\n", start, end);
      exit(1);
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
      return orig_cm->t[start][0];
    }
  to_return = 0.;
  
  sub_psi = malloc(sizeof(double) * (end - start + 1));
  /* Initialize sub_psi[0]. Need to check if we need to ignore the probability
   * mass from the ORIG_CM insert state(s) that maps to the HMM insert state of this node 
   * (these insert states are hns2cs_map[k][1][0] and (potentially) hns2cs_map[k][1][1]) 
   * that goes through "start" (which must map to either the M or 
   * D state of this HMM node).
   */
  sub_psi[0] = 1.; /* have to start in "start" */
  
  if((orig_cm->sttype[start] != IL_st && orig_cm->sttype[start] != IR_st) &&
     (orig_cm->sttype[end] != IL_st && orig_cm->sttype[end] != IR_st))
    {
      /* 08.31.06 TO DO: NOT SURE HOW TO DO THE EQUIVALENT WITH CM TO CM 
       * I THINK I MIGHT NEED TO KNOW THE INDEX OF THE INSERTS IN THE sub_cm
       * NODE THAT CONTAINS v_s

      insert_to_start = 0.;
      if(hns2cs_map[k][1][0] < start) 
	insert_to_start = psi[hns2cs_map[k][1][0]] * cm_sum_subpaths(cm, hns2cs_map[k][1][0], start, tmap, cs2hn_map, cs2hs_map, hns2cs_map, k, psi);
      if((hns2cs_map[k][1][1] != -1) && (hns2cs_map[k][1][1] < start))
	insert_to_start += psi[hns2cs_map[k][1][1]] * cm_sum_subpaths(cm, hns2cs_map[k][1][1], start, tmap, cs2hn_map, cs2hs_map, hns2cs_map, k, psi);
      sub_psi[0] -= insert_to_start / psi[start];
      */
    }
  /* note: when cm_sum_subpaths is called recursively above
   * it will never result in another recursive call, 
   * because its "start" is an insert state.  
   */
  
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
      /*if((v != end && is_insert) && ((cs2hn_map[v][0] == k) || cs2hn_map[v][1] == k))*/
      /* HEREHEREHERE */
      /* DO I HAVE TO CHECK 4 POSSIBILITIES? HERE */

      if((v != end && is_insert) && ((sub_cm->nodemap[orig2sub_smap[v][0]] == sub_cm->nodemap[orig2sub_smap[start][0]] ||
				      sub_cm->nodemap[orig2sub_smap[v][1]] == sub_cm->nodemap[orig2sub_smap[start][0]]) ||
				     (sub_cm->nodemap[orig2sub_smap[v][0]] == sub_cm->nodemap[orig2sub_smap[start][1]] ||
				      sub_cm->nodemap[orig2sub_smap[v][1]] == sub_cm->nodemap[orig2sub_smap[start][1]])))
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

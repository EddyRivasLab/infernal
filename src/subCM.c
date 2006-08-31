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

static void map_orig2sub_cm(CM_t *orig_cm, CM_t *sub_cm, int **orig2sub_smap, int **sub2orig_smap,
			    int sub_start, int sub_end);
static void map_orig2sub_cm_helper(int **orig2sub_smap, int **sub2orig_smap, int v_o, int v_s);

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
  /* first just copy ct array for model boundaries from con->ct */
  for (cpos = (model_start-1); cpos < model_end; cpos++)
    {
      sub_cpos = cpos - (model_start-1);
      sub_ct[sub_cpos] = con->ct[cpos];
    }
  /* second remove structure outside structural boundaries */
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
	      sub_ct[cpos] = -1;
	    }
      printf("A sub_ct[cpos=%3d]: %3d\n\n", cpos, sub_ct[cpos]);
    }

  /* Construct the new structure ss_cons based on the template CM 
   * ct array.
   * We could do this similar to how display.c::CreateCMConsensus()
   * does it to get the fully formatted WUSS ([{<>}]) string but 
   * lazily we just do <> bps here.
   */
  sub_cstr = MallocOrDie(sizeof(char) * (model_end - model_start + 1));
  for (cpos = (model_start-1); cpos < model_end; cpos++)
    {
      sub_cpos = cpos - (model_start-1);
      /*printf("sub_ct[cpos=%3d]: %3d\n", cpos, sub_ct[cpos]);*/
      if(sub_ct[sub_cpos] == -1)         sub_cstr[sub_cpos] = '.'; 
      else if (sub_ct[sub_cpos]  > cpos) sub_cstr[sub_cpos] = '<';
      else if (sub_ct[sub_cpos]  < cpos) sub_cstr[sub_cpos] = '>';
      else Die("ERROR: weird error in BuildSubCM()\n");
    }

  /* Build the new subCM given the new consensus structure. But don't
   * parameterize it yet.
   */
  ConsensusModelmaker(sub_cstr, con->clen, &sub_cm, &mtr);

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
  map_orig2sub_cm(orig_cm, sub_cm, orig2sub_smap, sub2orig_smap, model_start,
		  model_end);

  printf("con->clen    : %d\n", con->clen);
  printf("con->cseq    : %s\n", con->cseq);
  printf("con->cstr    : %s\n", con->cstr);
  printf("struct_start : %d\n", struct_start);
  printf("struct_end   : %d\n", struct_end);
  printf("model_start  : %d\n", model_start);
  printf("model_end    : %d\n", model_end);
  printf("\n\norig struct: %s\n", con->cstr);
  printf("\n\nnew struct : %s\n", sub_cstr);

  exit(1);

  orig_psi = malloc(sizeof(double) * orig_cm->M);
  make_tmap(&tmap);
  CMZero(sub_cm);
  CMSetNullModel(sub_cm, orig_cm->null);

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
  printf("in ConsensusModelmaker: %s\n", ss_cons);

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
 * map_orig2sub_cm()
 * derived loosely from hmmband.c::CP9_map_cm2hmm_and_hmm2cm()
 * Purpose:  Determine maps between a template CM and a sub CM, which probably
 *           has less structure (less MATP n_odes) than the template.
 * Args:    
 * CM_t *orig_cm       - the original template CM
 * CM_t *sub_cm        - the sub CM
 * int **orig2sub_smap - 2D state map from orig_cm (template) to sub_cm.
 *                       1st dimension - state index in orig_cm 
 *                       2nd D - 2 elements for up to 2 matching sub_cm states, .
 *                       Allocated here, responsibility of caller to free.
 * int **sub2orig_smap - 2D state map from orig_cm (template) to sub_cm.
 *                       1st dimension - state index in sub_cm (0..sub_cm->M-1)
 *                       2nd D - 2 elements for up to 2 matching orig_cm states, 
 *                       Allocated here, responsibility of caller to free.
 * int sub_start       - first consensus column in orig_cm that sub_cm models
 * int sub_end         - lasst consensus column in orig_cm that sub_cm models
 * Returns: (void) 
 */
void
map_orig2sub_cm(CM_t *orig_cm, CM_t *sub_cm, int **orig2sub_smap, int **sub2orig_smap,
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
      if(v_o < sub_cm->M)
	{
	  sub2orig_smap[v_o][0] = -1;
	  sub2orig_smap[v_o][1] = -1;
	}
    }

  /* Using the *cc_node_maps, traverse the consensus columns left 
   * to right and determine state mappings for MATP, MATL and MATR nodes
   * first.
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
      else /* n_o and n_s are of different node types. */
	/* case 1: n_o is a MATP, n_s is a MATL 
	 * case 2: n_o is a MATP, n_s is a MATR 
	 * any other cases are impossible */
	{
	  v_o = orig_cm->nodemap[n_o];
	  v_s =  sub_cm->nodemap[n_s];
	  if(orig_cm->ndtype[n_o] != MATP_nd)
	    Die("ERROR 3 in map_orig2sub_cm()\n");
	  if(sub_cm->ndtype[n_s] == MATL_nd)
	    {
	      /* Case by case (uchh..)
	      /* v_o is MATP_MP */
	      /* v_s is MATL_ML */
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
	  else if(sub_cm->ndtype[n_s] == MATR_nd)
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
	} /* end of else (n_o and n_s are different node types */
    }     
  
  /* TO FIX: for now just map the corresponding ROOT states to each other,
   * this isn't correct because a ROOT_IL in orig_cm will emit before the 
   * consensus column 1 while a ROOT_IL in sub_cm will emit before the
   * consensus column sub_start (which may be > 1)
   */
  v_o = 0; /* ROOT_S */
  v_s = 0; /* ROOT_S */
  map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);
  v_o++; /* v_o is ROOT_IL */
  v_s++; /* v_o is ROOT_IL */
  map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);
  v_o++; /* v_o is ROOT_IR */
  v_s++; /* v_o is ROOT_IR */
  map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);

  /* NOTE we ignore END_E states because they won't ever affect emission
   * OR transition probabilities.
   */

  /* map remaining node types: BIF, BEGL, BEGR */
  n_s = sub_cm->nodes - 1;
  n_o = orig_cm->nodes - 1;

  for(n_s = sub_cm->nodes - 1; n_s > 0; n_s--)
    {
      if(sub_cm->ndtype[n_s] == BIF_nd)
	{
	  found_bif = FALSE;
	  while(!found_bif)
	    {
	      /* Find the BIF node in orig_cm that maps to n_s. 
	       * The corresponding orig_cm node is the highest
	       * numbered node n_o that satisfies:
	       *    sub_emap[lpos] >= orig_emap[lpos] 
	       * && sub_emap[rpos] <= orig_emap[rpos] 
	       */
	      if(orig_cm->ndtype[n_o] != BIF_nd ||
		 (orig_cm->ndtype[n_o] == BIF_nd &&
		  (sub_emap->lpos[n_s] < orig_emap->lpos[n_o] ||
		   sub_emap->rpos[n_s] < orig_emap->rpos[n_o])))
		{ n_o--; } 
	      else
		{
		  printf("BIF match sub n: %d orig n: %d\n", n_o, n_s);
		  found_bif = TRUE;
		  /* We have a match */
		  v_o = orig_cm->nodemap[n_o]; /* v_o is BIF_B */
		  v_s =  sub_cm->nodemap[n_s]; /* v_s is BIF_B */
		  map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);
		  /* we also know the BEGL and BEGR's map */
		  /* BEGL first */
		  v_o = orig_cm->cfirst[orig_cm->nodemap[n_o]]; /* v_o is BEGL_S */
		  v_s =  sub_cm->cfirst[sub_cm->nodemap[n_s]];       /* v_s is BEGL_S */
		  map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);
		  /* now BEGR */
		  v_o = orig_cm->cnum[orig_cm->nodemap[n_o]]; /* v_o is BEGR_S */
		  v_s =  sub_cm->cnum[sub_cm->nodemap[n_s]];       /* v_s is BEGR_S */
		  map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);
		  /* and finally BEGR_IL */
		  /* TO FIX! this orig_cm BEGR_IL doesn't necessarily map 
		   * to the sub_cm BEGR_IL because these states may emit 
		   * in front of different consensus columns. */
		  v_o++; /* v_o is BEGR_IL */
		  v_s++; /* v_o is BEGR_IL */
		  map_orig2sub_cm_helper(orig2sub_smap, sub2orig_smap, v_o, v_s);

		  n_o--;
		}
	      if(n_o == 0)
		Die("ERROR no orig_cm node maps to BIF node %d of sub_cm\n", n_s);
	    }
	}
    }	  

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
      if(v_o == -1 && sub_cm->sttype[v_s] != E_st)
	Die("ERROR sub_cm state: %d type: %s node type: %s doesn't map to any state in orig_cm\n", v_s, sttypes[sub_cm->sttype[v_s]], nodetypes[sub_cm->ndtype[n_s]]);

      n_o = orig_cm->ndidx[v_o];
      
      printf("sub v:%4d %6s%6s | orig v:%4d %6s%6s\n", v_s, nodetypes[sub_cm->ndtype[n_s]], sttypes[sub_cm->sttype[v_s]], v_o, nodetypes[orig_cm->ndtype[n_o]], sttypes[orig_cm->sttype[v_o]]);

      /* check to make sure orig2sub_smap is consistent */
      if(orig2sub_smap[v_o][0] != v_s && orig2sub_smap[v_o][1] != v_s)
	{
	  Die("ERROR inconsistency; neither orig2sub_smap[%d][0] and [1] is %d\n", v_o, v_s);
	}

      v_o = sub2orig_smap[v_s][1];
      if(v_o != -1)
	{
	  n_o = orig_cm->ndidx[v_o];
	  printf("                        | orig v:%4d %6s%6s\n", v_o, nodetypes[orig_cm->ndtype[n_o]], sttypes[orig_cm->sttype[v_o]]);
	  /* check to make sure orig2sub_smap is consistent */
	  if(orig2sub_smap[v_o][0] != v_s && orig2sub_smap[v_o][1] != v_s)
	    {
	      Die("ERROR inconsistency; neither orig2sub_smap[%d][0] and [1] is %d\n", v_o, v_s);
	    }
	}
    }
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


/* destruct.c
 * EPN 07.25.06
 * 
 * Removing structure from the CM when aligning non-full length seqs.
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

/* Function:  DestructCM()
 * EPN 08.28.06 
 *
 * Args:      orig_cm    - the original model, which we're going to (potentially)
 *                         remove some structure from
 *            ret_cm     - the new CM built from orig_cm with some structure removed.
 *            start_pos  - the first position (consensus column) we want to keep structure
 *                         for.
 *            end_pos    - the last position we want to keep structure for
 *            
 * Returns:   (void)
 */
void
DestructCM(CM_t *orig_cm, CM_t **ret_cm, int start_pos, int end_pos)
{
  CM_t            *ds_cm;       /* new covariance model, the destructed template */
  CMConsensus_t   *con;         /* growing consensus info for orig_cm*/
  Parsetree_t     *mtr;         /* master structure tree from the alignment*/
  char     *ds_cstr;            /* destructed consensus structure display string  */
  int      *ds_ct;		/* 0..con->clen-1 base pair partners array         */
  char **nodetypes;
  int n;
  int cpos;
  CMEmitMap_t *orig_emap;         /* consensus emit map for the original, template CM */
  CMEmitMap_t *ds_emap;         /* consensus emit map for the destructed CM */
  
  nodetypes = malloc(sizeof(char *) * 8);
  nodetypes[0] = "BIF";
  nodetypes[1] = "MATP";
  nodetypes[2] = "MATL";
  nodetypes[3] = "MATR";
  nodetypes[4] = "BEGL";
  nodetypes[5] = "BEGR";
  nodetypes[6] = "ROOT";
  nodetypes[7] = "END";

  /* Get the consensus sequence and consensus structure information from the original CM */
  con = CreateCMConsensus(orig_cm, 3.0, 1.0);
  printf("con->cseq: %s\n", con->cseq);
  printf("con->cstr: %s\n", con->cstr);
  printf("start_pos: %d\n", start_pos);
  printf("end_pos  : %d\n", end_pos);

  /* Fill a new ct array for the destructed CM. First copy the template (original) CMs
   * ct array. Next, eliminate any structure that lies outside the region from
   * start_pos and end_pos, the positions the HMM told us are likely to be the end points
   * of the alignment.
   */

  ds_ct = MallocOrDie(sizeof(int) * (con->clen));
  for (cpos = 0; cpos < con->clen; cpos++)
    ds_ct[cpos] = con->ct[cpos];
  for (cpos = 0; cpos < con->clen; cpos++)
    {
      printf("B ds_ct[cpos=%3d]: %3d\n", cpos, ds_ct[cpos]);
      if ((cpos+1) < start_pos || (cpos+1) > end_pos) /* cpos goes 1..clen, but ct is indexed
						       * 0..clen-1.*/
	{ 
	  if (ds_ct[cpos] != -1)  ds_ct[ds_ct[cpos]] = -1; /* CreateCMConsensus() uses
							    * -1 in ct[] to indicate single
							    * stranded (different convention
							    * than WUSS2ct(). */
	  ds_ct[cpos] = -1;
	}
      printf("A ds_ct[cpos=%3d]: %3d\n\n", cpos, ds_ct[cpos]);
    }

  /* Construct the new structure ss_cons based on the template CM 
   * ct array.
   * We should do this similar to how display.c::CreateCMConsensus()
   * does it to get the fully formatted WUSS ([{<>}]) string but 
   * lazily we just do <> bps here.
   */
  ds_cstr = MallocOrDie(sizeof(char) * (con->clen + 1));
  for (cpos = 0; cpos < con->clen; cpos++)
    {
      /*printf("ds_ct[cpos=%3d]: %3d\n", cpos, ds_ct[cpos]);*/
      if(ds_ct[cpos] == -1)         ds_cstr[cpos] = '.'; 
      else if (ds_ct[cpos]  > cpos) ds_cstr[cpos] = '<';
      else if (ds_ct[cpos]  < cpos) ds_cstr[cpos] = '>';
      else Die("ERROR: weird error in DestructCM()\n");
    }

  /* Build the new destructed CM given the new consensus structure. But don't
   * parameterize it yet.
   */
  ConsensusModelmaker(ds_cstr, con->clen, &ds_cm, &mtr);

  printf("\n\norig struct: %s\n", con->cstr);
  printf("\n\nnew struct : %s\n", ds_cstr);

  /* Parameterize the new CM based on the template CM. 
   * First we need the emit maps for the template CM and
   * the new, destructed CM.
   */
  orig_emap = CreateEmitMap(orig_cm);
  ds_emap = CreateEmitMap(ds_cm);
  printf("\n\n\nDumping original CM emitmap\n");
  DumpEmitMap(stdout, orig_emap, orig_cm);
  printf("\n\n\nDumping destructed CM emitmap\n");
  DumpEmitMap(stdout, ds_emap, ds_cm);

  

  free(ds_cstr);
  free(ds_ct);
  free(nodetypes);
  FreeEmitMap(orig_emap);
  FreeEmitMap(ds_emap);
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

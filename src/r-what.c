/* r-what.c
 * DLK
 *
 * *********************************************
 * @LICENSES
 * *********************************************
 */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>

#include "squid.h"

#include "structs.h"
#include "funcs.h"

#include "dlk_priorityqueue.h"
#include "r-what.h"

/* Function: PA_Copy()
 * Author:   DLK
 *
 * Purpose:  Copy a partial alignment struct
 */
PA_t*
PA_Copy(PA_t *orig)
{
  PA_t *dup;
  dup = malloc(sizeof(PA_t));

  dup->init_v = orig->init_v;
  dup->init_j = orig->init_j;
  dup->init_d = orig->init_d;
	 
  dup->cur_v = orig->cur_v;
  dup->cur_j = orig->cur_j;
  dup->cur_d = orig->cur_d;
       
  dup->current_sc = orig->current_sc;
  dup->upper_bound_sc = orig->upper_bound_sc;

  dup->terminated = orig->terminated;
		   
  return dup;
}

/* Function: BPA_Copy()
 * Author:   DLK
 *
 * Purpose:  Copy a branched partial alignment struct
 */
BPA_t*
BPA_Copy(BPA_t *orig)
{
  BPA_t *dup;
  dup = malloc(sizeof(BPA_t));

  if (orig->chunk  != NULL)
    dup->chunk = PA_Copy(orig->chunk);
  else
    dup->chunk = NULL;
  if (orig->left_child != NULL)
    dup->left_child  = BPA_Copy(orig->left_child);
  else
    dup->left_child  = NULL;
  if (orig->right_child != NULL)
    dup->right_child = BPA_Copy(orig->right_child);
  else
    dup->right_child = NULL;

  dup->terminated = orig->terminated;

  return dup;
}

/* Function: BPA_Free()
 * Author:   DLK
 *
 * Purpose:  Free memory from BPA
 */
void
BPA_Free(BPA_t *root)
{
  if (root != NULL) {
    if (root->right_child != NULL) {
      BPA_Free(root->right_child);
    }
    if (root->left_child != NULL) {
      BPA_Free(root->left_child);
    }
    if (root->chunk != NULL) {
      free(root->chunk);
    }
    free(root);
  }
  
  return;
}

/* Function: BPA_Update()
 * Author:   DLK
 *
 * Purpose:  Re-evaluate and update overall BPA properties
 */
void
BPA_Update(BPA_t *root)
{
  if (root != NULL) {
    BPA_Update(root->left_child);
    BPA_Update(root->right_child);
    root->current_sc = root->chunk->current_sc;
    root->upper_bound_sc = 0.0;
    if (root->left_child != NULL) {
      root->current_sc += root->left_child->current_sc;
      root->upper_bound_sc += root->left_child->upper_bound_sc;
    }
    if (root->right_child != NULL) {
      root->current_sc += root->right_child->current_sc;
      root->upper_bound_sc += root->right_child->upper_bound_sc;
    }
    if (root->left_child == NULL && root->right_child == NULL) {
      root->upper_bound_sc = root->chunk->upper_bound_sc;
      if (root->chunk->terminated) { root->terminated = 1; }
    }
    if (root->left_child->terminated && root->right_child->terminated) {
      root->terminated = 1;
    }
  }

  return;
}

/* Function: BPA_Current_Score()
 * Author:   DLK
 *
 * Purpose:  Calculate current score for a branched partial alignment
 */
float
BPA_Current_Score(BPA_t *root)
{
  float score = 0.0;

  if (root != NULL) {
    score = root->chunk->current_sc;
    score += BPA_Current_Score(root->left_child);
    score += BPA_Current_Score(root->right_child);
  }

  return score;
}

/* Function: BPA_Upper_Bound()
 * Author:   DLK
 *
 * Purpose:  Calculate the upper bound for a branched partial alignment
 *           Note that this is an incremental score:
 *           current score + upper bound = bound on total alignment score
 */
float
BPA_Upper_Bound(BPA_t *root)
{
  float score = 0.0;

  if (root != NULL) {
    if (root->left_child == NULL & root->right_child == NULL)
      score += root->chunk->upper_bound_sc;
    else {
      score += BPA_Upper_Bound(root->left_child);
      score += BPA_Upper_Bound(root->right_child);
    }
  }

  return score;
}

/* Function: MaxSubsequenceScore()
 * Author:   DLK
 *
 * Purpose:  For a cm, compute a matrix of maximum possible
 *           subsequence scores for each state v and length
 *           d, up to a window length W
 *
 * Notes:    In the matrix, a value max_sc[v][d] represents
 *           the max sc for a sequence of length d accounted
 *           for by the states BELOW v in the model.  That
 *           is, any bases emitted by v are NOT included in
 *           the size of d.
 *
 * Args:     cm		- the covariance model
 *           W		- max d: max size of hit
 *           ret_max_sc - RETURN: matrix of max scores
 *
 * Returns:  ret_max_sc
 *           ret_max_sc allocated here; caller free's w/ free()
 */
void
MaxSubsequenceScore(CM_t *cm, int W, float ***ret_max_sc)
{
  int v;		/* state index */
  int d;		/* subsequence length */
  float **max_sc;	/* matrix of max scores */
  int i,j;		/* counter variables */
  float sc;		/* max score for a state */

  int w,y;		/* child state indices */
  int yoffset;		/* offset to a child state */

  char seq[] = "ACGU";	/* sequence enumerating the bases */
  int L = 4;		/* length of sequence */
  char *dsq;

  dsq = DigitizeSequence(seq,L);

  max_sc = MallocOrDie(sizeof(float *) * cm->M);
  max_sc[0] = MallocOrDie(sizeof(float) * cm->M * W);
  for (v=1; v<cm->M; v++)
    max_sc[v] = max_sc[0] + W*v;

  /* Initialize for d=0 */
  for (v=cm->M-1; v>=0; v--)
  {
    max_sc[v][0] = IMPOSSIBLE;

    if (cm->sttype[v] == E_st) max_sc[v][0] = 0;
    else if (cm->sttype[v] == B_st)
    {
      w = cm->cfirst[v];
      y = cm->cnum[v];
      max_sc[v][0] = max_sc[w][0] + max_sc[y][0];
    }
    else
    {
      y = cm->cfirst[v];
      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
      {
	if (cm->sttype[y+yoffset] == D_st ||
	    cm->sttype[y+yoffset] == S_st ||
	    cm->sttype[y+yoffset] == E_st ||
	    cm->sttype[y+yoffset] == B_st)
	{ sc = cm->tsc[v][yoffset] + max_sc[y+yoffset][0]; }
	else
	{ sc = IMPOSSIBLE; }
	if (sc > max_sc[v][0]) { max_sc[v][0] = sc; }
      }
    }
    if (max_sc[v][0] < IMPROBABLE) { max_sc[v][0] = IMPOSSIBLE; }
  }

  /* Initialize for d=1 */
  for (v=cm->M-1; v>=0; v--)
  {
    max_sc[v][1] = IMPOSSIBLE;
    if      (cm->sttype[v] == E_st) { max_sc[v][1] = IMPOSSIBLE; }
    else if (cm->sttype[v] == B_st) 
    {
      w = cm->cfirst[v];
      y = cm->cnum[v];
      max_sc[v][1] = max_sc[w][1] + max_sc[y][0];
      sc           = max_sc[w][0] + max_sc[y][1];
      if (sc > max_sc[v][1]) { max_sc[v][1] = sc; }
    }
    else
    {
      y = cm->cfirst[v];
      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
      {
	if      (cm->sttype[y+yoffset] == D_st ||
	         cm->sttype[y+yoffset] == S_st ||
	         cm->sttype[y+yoffset] == E_st ||
	         cm->sttype[y+yoffset] == B_st)
        { sc = cm->tsc[v][yoffset] + max_sc[y+yoffset][1]; }
	else if (cm->sttype[y+yoffset] == ML_st ||
	         cm->sttype[y+yoffset] == MR_st ||
		 cm->sttype[y+yoffset] == IL_st ||
		 cm->sttype[y+yoffset] == IR_st)
	{
	  sc = IMPOSSIBLE;
          for (i=0; i<L; i++)
	  { if (cm->esc[y+yoffset][(int) dsq[i]] > sc) { sc = cm->esc[y+yoffset][(int) dsq[i]]; } }
	  sc += cm->tsc[v][yoffset] + max_sc[y+yoffset][0];
	}
	else
	{ sc = IMPOSSIBLE; }
	if (sc > max_sc[v][1]) { max_sc[v][1] = sc; }
      }
    }
    if (max_sc[v][1] < IMPROBABLE) { max_sc[v][1] = IMPOSSIBLE; }
  }

  /* Fill matrix for d >= 2 */
  for (d=2; d<W; d++)
  {
    for (v=cm->M-1; v>=0; v--)
    {
      max_sc[v][d] = IMPOSSIBLE;
      if      (cm->sttype[v] == E_st) { max_sc[v][d] = IMPOSSIBLE; }
      else if (cm->sttype[v] == B_st) 
      {
        w = cm->cfirst[v];
        y = cm->cnum[v];
	sc = IMPOSSIBLE;
	for (i=0; i<=d; i++)
	{
          sc = max_sc[w][i] + max_sc[y][d-i];
          if (sc > max_sc[v][d]) { max_sc[v][d] = sc; }
	}
      }
      else
      {
        y = cm->cfirst[v];
        for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
        {
          if      (cm->sttype[y+yoffset] == D_st ||
                   cm->sttype[y+yoffset] == S_st ||
                   cm->sttype[y+yoffset] == E_st ||
                   cm->sttype[y+yoffset] == B_st)
          { sc = cm->tsc[v][yoffset] + max_sc[y+yoffset][d]; }
          else if (cm->sttype[y+yoffset] == ML_st ||
                   cm->sttype[y+yoffset] == MR_st ||
		   cm->sttype[y+yoffset] == IL_st ||
		   cm->sttype[y+yoffset] == IR_st)
          {
            sc = IMPOSSIBLE;
            for (i=0; i<L; i++)
              if (cm->esc[y+yoffset][(int) dsq[i]] > sc)
	        sc = cm->esc[y+yoffset][(int) dsq[i]];  
            sc += cm->tsc[v][yoffset] + max_sc[y+yoffset][d-1];
          }
          else if (cm->sttype[y+yoffset] == MP_st)
          {
	    sc = IMPOSSIBLE;
	    for (i=0; i<L; i++)
	      for (j=0; j<L; j++)
		if (cm->esc[y+yoffset][(int) (dsq[i]*Alphabet_size+dsq[j])] > sc)
		  sc = cm->esc[y+yoffset][(int) (dsq[i]*Alphabet_size+dsq[j])];
	    sc += cm->tsc[v][yoffset] + max_sc[y+yoffset][d-2];
	  }

          if (sc > max_sc[v][d]) { max_sc[v][d] = sc; }
        }
      }
      if (max_sc[v][d] < IMPROBABLE) { max_sc[v][d] = IMPOSSIBLE; }
    }
  }

  *ret_max_sc = max_sc;

  return;
}

/* Function: AstarExtension()
 * Author:   DLK
 *
 * Purpose:  An SCFG alignment algorithm equivalent to CYK,
 *           but using A* (Astar) to operate in the outside-
 *           to-inside direction, extending from a given 
 *           initial state
 *
 * Args:     cm      - the covariance model
 *           dsq     - the digitized sequence
 *           init_v  - initial v state in the model
 *           init_j  - initial j position in the sequence
 *           lower_d - 
 *           upper_d -
 *                   - specify the range of d values to be
 *                     considered.  lower_d == upper_d if
 *                     extending from known pair of positions
 *           init_sc - initial score
 *           max_sc  - a matrix of maximum subsequence scores
 *           cutoff  - minimum score threshold - used for pruning
 *
 * Returns:  PA_t    - partial alignment
 *
 *
 */
PA_t*
AstarExtension(CM_t *cm, char *dsq, int init_v, int init_j, int lower_d, int upper_d,
        float init_sc, float **max_sc, float cutoff)
{
  int d;
  int i;
  int yoffset;
  float tsc, esc, total_sc;
  PA_t *pa;
  PA_t *child_pa;

  PriorityQueue_t *alignPQ;

  alignPQ = CreatePQ();

  for (d=lower_d; d<=upper_d; d++) {
    pa = malloc(sizeof(PA_t));

    pa->init_v = pa->cur_v = init_v;
    pa->init_j = pa->cur_j = init_j;
    pa->init_d = pa->cur_d = d;
    pa->current_sc = init_sc;
    pa->upper_bound_sc = max_sc[init_v][d-2];	/* i and j already accounted for */

    total_sc = pa->current_sc + pa->upper_bound_sc;
    if (total_sc > cutoff) { EnqueuePQ(alignPQ, pa, total_sc); }
    else { free(pa); }
  }

  while ((pa = DequeuePQ(alignPQ)) != NULL) {
    if (cm->sttype[pa->cur_v] == E_st || cm->sttype[pa->cur_v] == B_st) {
      break;
    }

    /* For each child state */
    for (yoffset = 0; yoffset < cm->cnum[pa->cur_v]; yoffset++) {
      /* Create partial alignment copy */
      child_pa = PA_Copy(pa);
      /* Calculate tsc and new v */
      tsc = cm->tsc[pa->cur_v][yoffset];
      child_pa->cur_v = cm->cfirst[pa->cur_v] + yoffset;

      /* Get new j,d and calculate esc */
      if      (cm->sttype[child_pa->cur_v] == D_st) {
        esc = 0.0;
      }
      else if (cm->sttype[child_pa->cur_v] == MP_st) {
        child_pa->cur_j--;
        child_pa->cur_d -= 2;
        i = child_pa->cur_j-child_pa->cur_d+1;
        if (dsq[i] < Alphabet_size && dsq[child_pa->cur_j] < Alphabet_size)
          esc = cm->esc[child_pa->cur_v][(int) (dsq[i]*Alphabet_size+dsq[child_pa->cur_j])];
        else
          esc = DegeneratePairScore(cm->esc[child_pa->cur_v],dsq[i],dsq[child_pa->cur_j]);
	//if (child_pa->cur_d < 2)
	  //esc = IMPOSSIBLE;
      }
      else if (cm->sttype[child_pa->cur_v] == ML_st || cm->sttype[child_pa->cur_v] == IL_st) {
        child_pa->cur_d--;
        i = child_pa->cur_j-child_pa->cur_d+1;
        if (dsq[i] < Alphabet_size) 
          esc = cm->esc[child_pa->cur_v][(int) dsq[i]];
        else
          esc = DegenerateSingletScore(cm->esc[child_pa->cur_v],dsq[i]);
	//if (child_pa->cur_d < 1)
	  //esc = IMPOSSIBLE;
      }
      else if (cm->sttype[child_pa->cur_v] == MR_st || cm->sttype[child_pa->cur_v] == IR_st) {
        child_pa->cur_j--;
        child_pa->cur_d--;
        if (dsq[child_pa->cur_j] < Alphabet_size)
          esc = cm->esc[child_pa->cur_v][(int) dsq[child_pa->cur_j]];
        else
          esc = DegenerateSingletScore(cm->esc[child_pa->cur_v],dsq[child_pa->cur_j]);
	//if (child_pa->cur_d < 1)
	  //esc = IMPOSSIBLE;
      }
      else if (cm->sttype[child_pa->cur_v] == E_st || cm->sttype[child_pa->cur_v] == B_st) {
	esc = 0.0;
      }

      child_pa->current_sc = child_pa->current_sc + tsc + esc;
      child_pa->upper_bound_sc = max_sc[child_pa->cur_v][child_pa->cur_d-2];

      total_sc = child_pa->current_sc + child_pa->upper_bound_sc;
      if (total_sc > cutoff) { EnqueuePQ(alignPQ, child_pa, total_sc); }
      else { free(child_pa); }
    }

    free(pa);
  }
  
  if (pa != NULL) {
    /* Top alignment info is in pa.  Store it */
    child_pa = PA_Copy(pa);
    free(pa);
  }
  else {
    child_pa = NULL;
  }

  /* Empty PQ and release memory */
  while ((pa = DequeuePQ(alignPQ)) != NULL) {
    free(pa);
  }
  FreePQ(alignPQ);

  return child_pa;
}

/* Function: AstarBAlign()
 * Author:   DLK
 *
 * Purpose:  An SCFG alignment algorithm equivalent to CYK,
 *           but using A* to operate in the outside-to-inside
 *           direction, extending from a given initial state
 *           and position
 *
 *           Includes bifurcation
 *
 *           Slow and not intended for primary scanning use -
 *           tends to recalculate values because alpha[v][j][d]
 *           is not stored
 *
 * Args:     cm      - the covariance model
 *           dsq     - the digitized sequence
 *           align   - initial partial alignment
 *           max_sc  - a matrix of maximu subsequence scores
 *           cutoff  - minimum score threshold - used for pruning
 *
 * Returns:  align   - finished partial alignment
 *           return value - 0 if alignment found, 1 otherwise
 */

int
AstarBAlign(CM_t *cm, char *dsq, BPA_t *init_align, float **max_sc, float cutoff)
{
  int retval;
  int i,d;
  int yoffset;
  float tsc, esc, total_sc;
  BPA_t *align;
  BPA_t *working;
  BPA_t *child_align;
  BPA_t *child_working;
  PriorityQueue_t *alignPQ;

  alignPQ = CreatePQ();
  align = BPA_Copy(init_align);

  BPA_Update(align);
  total_sc = align->current_sc + align->upper_bound_sc;
  if (total_sc > cutoff) { EnqueuePQ(alignPQ, align, total_sc); }
  else { BPA_Free(align); }

  while (((align = DequeuePQ(alignPQ)) != NULL) && !(align->terminated)) {
    /* Move working pointer to next unterminated chunk */
    working = align;
    while (working->chunk->terminated) {
      if      ((working->left_child != NULL) && !working->left_child->terminated)
        { working = working->left_child; }
      else if ((working->right_child != NULL) && !working->right_child->terminated)
        { working = working->right_child; }
      else {
	BPA_Update(align);
	total_sc = align->current_sc + align->upper_bound_sc;
	if (total_sc > cutoff) { EnqueuePQ(alignPQ, align, total_sc); }
        else { BPA_Free(align); }
      }
    }

    /* Extend chunk */
    if (cm->sttype[working->chunk->cur_v] == E_st) {
      /* terminate chunk */
      working->chunk->terminated = 1;
      working->terminated = 0;
      BPA_Update(align);
      total_sc = align->current_sc + align->upper_bound_sc;
      if (total_sc > cutoff) { EnqueuePQ(alignPQ, align, total_sc); }
      else { BPA_Free(align); }
    }
    else if (cm->sttype[working->chunk->cur_v] == B_st) {
      /* bifurcate */
      for (d=0; d<=working->chunk->cur_v; d++) {
        child_align = BPA_Copy(align);
        child_working = child_align;
        while (child_working->chunk->terminated) {
          if      ((child_working->left_child != NULL) && !child_working->left_child->terminated)
            { child_working = child_working->left_child; }
          else if ((child_working->right_child != NULL) && !child_working->right_child->terminated)
            { child_working = child_working->right_child; }
          else {
            fprintf(stderr,"WARNING: copy doesn't match original!\n");
            exit(1);
          }
        }

        child_working->terminated = 1;
        /* Allocate new memory */
        child_working->left_child = malloc(sizeof(BPA_t));
        child_working->left_child->chunk = malloc(sizeof(PA_t));
        child_working->left_child->left_child = NULL;
        child_working->left_child->right_child = NULL;
        child_working->right_child = malloc(sizeof(BPA_t));
        child_working->right_child->chunk = malloc(sizeof(PA_t));
        child_working->right_child->left_child = NULL;
        child_working->right_child->right_child = NULL;

        /* Initialize left and right children */
        /* Note that the initial j,d coordinates are outside the sequence chunk   *
         * actually covered by v and the states beneath it on the parse tree      *
         * When reporting, use init_i+1 and init_j-1 (with care!)                 */
        child_working->left_child->chunk->init_v = cm->cfirst[child_working->chunk->cur_v];
        child_working->left_child->chunk->cur_v  = cm->cfirst[child_working->chunk->cur_v];
        child_working->left_child->chunk->init_j = child_working->chunk->cur_j - d + 2;
        child_working->left_child->chunk->cur_j  = child_working->chunk->cur_j - d + 2;
        child_working->left_child->chunk->init_d = child_working->chunk->cur_d - d + 2;
        child_working->left_child->chunk->cur_d  = child_working->chunk->cur_d - d + 2;

        child_working->right_child->chunk->init_v = cm->cnum[child_working->chunk->cur_v];
        child_working->right_child->chunk->cur_v  = cm->cnum[child_working->chunk->cur_v];
        child_working->right_child->chunk->init_j = child_working->chunk->cur_j;
        child_working->right_child->chunk->cur_j  = child_working->chunk->cur_j;
        child_working->right_child->chunk->init_d = d;
        child_working->right_child->chunk->cur_d  = d;

        child_working->left_child->chunk->current_sc = 0.0;
        child_working->right_child->chunk->current_sc = 0.0;
        child_working->left_child->chunk->upper_bound_sc =
          max_sc[child_working->left_child->chunk->cur_v][child_working->left_child->chunk->cur_d - 2];
        child_working->right_child->chunk->upper_bound_sc =
          max_sc[child_working->right_child->chunk->cur_v][child_working->right_child->chunk->cur_d - 2];

        BPA_Update(child_align);
        total_sc = child_align->chunk->current_sc + child_align->chunk->upper_bound_sc;
        if (total_sc > cutoff) { EnqueuePQ(alignPQ, child_align, total_sc); }
        else { BPA_Free(child_align); }
      }
      BPA_Free(align);
    }
    else {
      /* simple extension */
      for (yoffset = 0; yoffset < cm->cnum[working->chunk->cur_v]; yoffset++) {
        child_align = BPA_Copy(align);
        child_working = child_align;
        while (child_working->chunk->terminated) {
          if      ((child_working->left_child != NULL) && !child_working->left_child->terminated)
            { child_working = child_working->left_child; }
          else if ((child_working->right_child != NULL) && !child_working->right_child->terminated)
            { child_working = child_working->right_child; }
          else {
            fprintf(stderr,"WARNING: copy doesn't match orignal!\n");
            exit(1);
          }
        }
      
        /* Calculate tsc and new v */
        tsc = cm->tsc[working->chunk->cur_v][yoffset];
        child_working->chunk->cur_v = cm->cfirst[working->chunk->cur_v] + yoffset;

        /* Get new j,d and calculate esc */
        if      (cm->sttype[child_working->chunk->cur_v] == D_st ||
                 cm->sttype[child_working->chunk->cur_v] == E_st ||
                 cm->sttype[child_working->chunk->cur_v] == B_st) {
          esc = 0.0;
        } 
        else if (cm->sttype[child_working->chunk->cur_v] == MP_st) {
          child_working->chunk->cur_j--;
          child_working->chunk->cur_d -= 2;
          i = child_working->chunk->cur_j - child_working->chunk->cur_d + 1;
          if (dsq[i] < Alphabet_size && dsq[child_working->chunk->cur_j] < Alphabet_size)
            esc = cm->esc[child_working->chunk->cur_v][(int) (dsq[i]*Alphabet_size+dsq[child_working->chunk->cur_j])];
          else
            esc = DegeneratePairScore(cm->esc[child_working->chunk->cur_v],dsq[i],dsq[child_working->chunk->cur_j]);
        }
        else if (cm->sttype[child_working->chunk->cur_v] == ML_st || cm->sttype[child_working->chunk->cur_v] == IL_st) {
          child_working->chunk->cur_d--;
          i = child_working->chunk->cur_j - child_working->chunk->cur_d + 1;
          if (dsq[i] < Alphabet_size)
            esc = cm->esc[child_working->chunk->cur_v][(int) dsq[i]];
          else
            esc = DegenerateSingletScore(cm->esc[child_working->chunk->cur_v],dsq[i]);
        }
        else if (cm->sttype[child_working->chunk->cur_v] == MR_st || cm->sttype[child_working->chunk->cur_v] == IR_st) {
          child_working->chunk->cur_j--;
          child_working->chunk->cur_d--;
          if (dsq[child_working->chunk->cur_j] < Alphabet_size)
            esc = cm->esc[child_working->chunk->cur_v][(int) dsq[child_working->chunk->cur_j]];
          else
            esc = DegenerateSingletScore(cm->esc[child_working->chunk->cur_v],dsq[child_working->chunk->cur_j]);
        }

        child_working->chunk->current_sc = child_working->chunk->current_sc + tsc + esc;
        child_working->chunk->upper_bound_sc = max_sc[child_working->chunk->cur_v][child_working->chunk->cur_d-2];

        BPA_Update(align);
        total_sc = child_align->chunk->current_sc + child_align->chunk->upper_bound_sc;
        if (total_sc > cutoff) { EnqueuePQ(alignPQ, child_align, total_sc); }
        else { BPA_Free(child_align); }
      }
      BPA_Free(align);
    }
  }

  if (align != NULL) {
    /* Top alignment is in align - store for return */
    BPA_Free(init_align);
    init_align = align;
    retval = 0;
  }
  else {
    retval = 1;
  }

  while ((align = DequeuePQ(alignPQ)) != NULL) {
    BPA_Free(align);
  }
  FreePQ(alignPQ);

  return retval;
}

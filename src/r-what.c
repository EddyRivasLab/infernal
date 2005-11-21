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
		   
  return dup;
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
 *
 * Returns:  PA_t    - partial alignment
 *
 *
 */
PA_t*
AstarExtension(CM_t *cm, char *dsq, int init_v, int init_j, int lower_d, int upper_d,
        float init_sc, float **max_sc)
{
  int d;
  int i;
  int yoffset;
  float tsc, esc;
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

    EnqueuePQ(alignPQ, pa, pa->current_sc + pa->upper_bound_sc);
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
	if (child_pa->cur_d < 2)
	  esc = IMPOSSIBLE;
      }
      else if (cm->sttype[child_pa->cur_v] == ML_st || cm->sttype[child_pa->cur_v] == IL_st) {
        child_pa->cur_d--;
        i = child_pa->cur_j-child_pa->cur_d+1;
        if (dsq[i] < Alphabet_size) 
          esc = cm->esc[child_pa->cur_v][(int) dsq[i]];
        else
          esc = DegenerateSingletScore(cm->esc[child_pa->cur_v],dsq[i]);
	if (child_pa->cur_d < 1)
	  esc = IMPOSSIBLE;
      }
      else if (cm->sttype[child_pa->cur_v] == MR_st || cm->sttype[child_pa->cur_v] == IR_st) {
        child_pa->cur_j--;
        child_pa->cur_d--;
        if (dsq[child_pa->cur_j] < Alphabet_size)
          esc = cm->esc[child_pa->cur_v][(int) dsq[child_pa->cur_j]];
        else
          esc = DegenerateSingletScore(cm->esc[child_pa->cur_v],dsq[child_pa->cur_j]);
	if (child_pa->cur_d < 1)
	  esc = IMPOSSIBLE;
      }
      else if (cm->sttype[child_pa->cur_v] == E_st || cm->sttype[child_pa->cur_v] == B_st) {
	esc = 0.0;
      }

      child_pa->current_sc = child_pa->current_sc + tsc + esc;
      child_pa->upper_bound_sc = max_sc[child_pa->cur_v][child_pa->cur_d-2];

      EnqueuePQ(alignPQ, child_pa, child_pa->current_sc + child_pa->upper_bound_sc);
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

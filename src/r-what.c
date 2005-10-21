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
      for (yoffset = 0; yoffset < cm->cnum[v]-1; yoffset++)
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
      for (yoffset = 0; yoffset < cm->cnum[v]-1; yoffset++)
      {
	if      (cm->sttype[y+yoffset] == D_st ||
	         cm->sttype[y+yoffset] == S_st ||
	         cm->sttype[y+yoffset] == E_st ||
	         cm->sttype[y+yoffset] == B_st)
        { sc = cm->tsc[v][yoffset] + max_sc[y+yoffset][1]; }
	else if (cm->sttype[y+yoffset] == ML_st ||
	         cm->sttype[y+yoffset] == MR_st)
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
      if      (cm->sttype[v] == E_st) { max_sc[v][d] = IMPOSSIBLE; }
      else if (cm->sttype[v] == B_st) 
      {
        w = cm->cfirst[v];
        y = cm->cnum[v];
	sc = IMPOSSIBLE;
	for (i=0; i<=d; i++)
	{
          sc = max_sc[w][i] + max_sc[y][d-i];
          if (sc > max_sc[v][1]) { max_sc[v][1] = sc; }
	}
      }
      else
      {
        y = cm->cfirst[v];
        for (yoffset = 0; yoffset < cm->cnum[v]-1; yoffset++)
        {
          if      (cm->sttype[y+yoffset] == D_st ||
                   cm->sttype[y+yoffset] == S_st ||
                   cm->sttype[y+yoffset] == E_st ||
                   cm->sttype[y+yoffset] == B_st)
          { sc = cm->tsc[v][yoffset] + max_sc[y+yoffset][d]; }
          else if (cm->sttype[y+yoffset] == ML_st ||
                   cm->sttype[y+yoffset] == MR_st)
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

/* scancyk.c
 * SRE, Thu May  2 11:50:48 2002 [AA 3050 SFO->STL]
 * SVN $Id$
 * 
 * CYK alignment: multihit, local, database scanning mode.
 * [xref STL6 p47]
 * 
 ***************************************************************** 
 * @LICENSE@
 ***************************************************************** 
 */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>

#include "squid.h"

#include "structs.h"
#include "funcs.h"

/*
 * Function: CreateResults ()
 * Date:     RJK, Mon Apr 1 2002 [St. Louis]
 * Purpose:  Creates a results type of specified size
 */
scan_results_t *CreateResults (int size) {
  scan_results_t *results;
  results = MallocOrDie (sizeof(scan_results_t));
  results->num_results = 0;
  results->num_allocated = size;
  results->data = MallocOrDie(sizeof(scan_result_node_t)*size);
  return (results);
}

/* Function: ExpandResults ()
 * Date:     RJK, Mon Apr 1 2002 [St. Louis]
 * Purpose:  Expans a results structure by specified amount
 */
void ExpandResults (scan_results_t *results, int additional) {
  results->data = ReallocOrDie (results->data, 
				sizeof(scan_result_node_t)*
				(results->num_allocated+additional));
  results->num_allocated+=additional;
}

/*
 * Function: FreeResults()
 * Date:     RJK, Mon Apr 1 2002 [St. Louis]
 * Purpose:  Frees a results structure
 */
void FreeResults (scan_results_t *r) {
  int i;
  if (r != NULL) {
    for (i=0; i<r->num_results; i++) {
      if (r->data[i].tr != NULL) {
	FreeParsetree(r->data[i].tr);
      }
    }
    free (r->data);
    free(r);
  }
}


/*
 * Function: compare_results()
 * Date:     RJK, Wed Apr 10, 2002 [St. Louis]
 * Purpose:  Compares two scan_result_node_ts based on score and returns -1
 *           if first is higher score than second, 0 if equal, 1 if first
 *           score is lower.  This results in sorting by score, highest
 *           first.
 */
int compare_results (const void *a_void, const void *b_void) {
  scan_result_node_t *a, *b;
 
  a = (scan_result_node_t *)a_void;
  b = (scan_result_node_t *)b_void;

  if (a->score < b->score)
    return (1);
  else if (a->score > b->score)
    return (-1);
  else if (a->start < b->start)
    return (1);
  else if (a->start > b->start)
    return (-1);
  else
    return (0);
}

/*
 * Function: remove_overlapping_hits ()
 * Date:     EPN, Tue Apr  3 14:36:38 2007
 * Plucked verbatim out of RSEARCH
 * 
 * RSEARCH date: RJK, Sun Mar 31, 2002 [LGA Gate D7]
 *
 * Purpose:  Given a list of hits, removes overlapping hits to produce
 * a list consisting of at most one hit covering each nucleotide in the
 * sequence.  Works as follows:
 * 1.  quicksort hits 
 * 2.  For each hit, sees if any nucleotide covered yet
 *     If yes, remove hit
 *     If no, mark each nt as covered
 */
void remove_overlapping_hits (scan_results_t *results, int L) 
{
  char *covered_yet;
  int x,y;
  int covered;
  scan_result_node_t swap;

  if (results == NULL)
    return;

  if (results->num_results == 0)
    return;

  covered_yet = MallocOrDie (sizeof(char)*(L+1));
  for (x=0; x<=L; x++)
    covered_yet[x] = 0;

  sort_results (results);

  for (x=0; x<results->num_results; x++) {
    covered = 0;
    for (y=results->data[x].start; y<=results->data[x].stop && !covered; y++) {
      if (covered_yet[y] != 0) {
	covered = 1;
      } 
    }
    if (covered == 1) {
      results->data[x].start = -1;        /* Flag -- remove later to keep sorted */
    } else {
      for (y=results->data[x].start; y<=results->data[x].stop; y++) {
	covered_yet[y] = 1;
      }
    }
  }
  free (covered_yet);

  for (x=0; x<results->num_results; x++) {
    while (results->num_results > 0 &&
	   results->data[results->num_results-1].start == -1)
      results->num_results--;
    if (x<results->num_results && results->data[x].start == -1) {
      swap = results->data[x];
      results->data[x] = results->data[results->num_results-1];
      results->data[results->num_results-1] = swap;
      results->num_results--;
    }
  }
  while (results->num_results > 0 &&
	 results->data[results->num_results-1].start == -1)
    results->num_results--;

  sort_results(results);
}

/*
 * Function: sort_results()
 * Date:    RJK,  Sun Mar 31, 2002 [AA Flight 2869 LGA->STL]
 * Purpose: Given a results array, sorts it with a call to qsort
 *
 */
void sort_results (scan_results_t *results) 
{
  qsort (results->data, results->num_results, sizeof(scan_result_node_t), compare_results);
}

/*
 * Function: report_hit()
 * Date:     RJK, Sun Mar 31, 2002 [LGA Gate D7]
 *
 * Given j,d, coordinates, a score, and a scan_results_t data type,
 * adds result into the set of reportable results.  Naively adds hit.
 *
 * Non-overlap algorithm is now done in the scanning routine by Sean's
 * Semi-HMM code.  I've just kept the hit report structure for convenience.
 */
void report_hit (int i, int j, int bestr, float score, scan_results_t *results) 
{

  if(results == NULL) 
    Die("in report_hit, but results is NULL\n");
  if (results->num_results == results->num_allocated) 
    ExpandResults (results, INIT_RESULTS);

  results->data[results->num_results].score = score;
  results->data[results->num_results].start = i;
  results->data[results->num_results].stop = j;
  results->data[results->num_results].bestr = bestr;
  results->data[results->num_results].tr = NULL;
  results->num_results++;
}

/* Function: CYKScan()
 * Date:     SRE, Thu May  2 11:56:11 2002 [AA 3050 SFO->STL]
 *
 * Purpose:  Scan a sequence for matches to a covariance model.
 *           Multiple nonoverlapping hits and local alignment.
 *
 * Args:     cm          - the covariance model
 *           dsq         - digitized sequence to search; i0..j0
 *           i0          - start of target subsequence (1 for full seq)
 *           j0          - end of target subsequence (L for full seq)
 *           W           - max d: max size of a hit
 *           cutoff      - minimum score to report 
 *           results     - scan_results_t to add to; if NULL, don't add to it
 *
 * Returns:  score of best overall hit
 */
float 
CYKScan(CM_t *cm, char *dsq, int i0, int j0, int W, 
	float cutoff, scan_results_t *results)
{
  float  ***alpha;              /* CYK DP score matrix, [v][j][d] */
  int      *bestr;              /* auxil info: best root state at alpha[0][cur][d] */
  float    *gamma;              /* SHMM DP matrix for optimum nonoverlap resolution */
  int      *gback;              /* traceback pointers for SHMM */ 
  float    *savesc;             /* saves score of hit added to best parse at j */
  int      *saver;		/* saves initial non-ROOT state of best parse ended at j */
  int       v;			/* a state index, 0..M-1 */
  int       w, y;		/* child state indices */
  int       yoffset;		/* offset to a child state */
  int       i,j;		/* index of start/end positions in sequence, 0..L */
  int       d;			/* a subsequence length, 0..W */
  int       k;			/* used in bifurc calculations: length of right subseq */
  int       prv, cur;		/* previous, current j row (0 or 1) */
  float     sc;			/* tmp variable for holding a score */
  int       jp;			/* index into BEGL_S decks: jp=j%(W+1) */
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       gamma_j;            /* j index in the gamma matrix, which is indexed 0..j0-i0+1, 
				 * while j runs from i0..j0 */
  int       gamma_i;            /* i index in the gamma* data structures */
  float     best_score;         /* Best overall score from semi-HMM to return */
  float     best_neg_score;     /* Best score overall score to return, used if all scores < 0 */
  int       bestd;              /* d value of best hit thus far seen for j (used if greedy strategy) */

  /*****************************************************************
   * alpha allocations.
   * The scanning matrix is indexed [v][j][d]. 
   *    v ranges from 0..M-1 over states in the model.
   *    j takes values 0 or 1: only the previous (prv) or current (cur) row
   *      with the exception of BEGL_S, where we have to have a whole W+1xW+1
   *      deck in memory, and j ranges from 0..W, and yes it must be square
   *      because we'll use a rolling pointer trick thru it
   *    d ranges from 0..W over subsequence lengths.
   * Note that E memory is shared: all E decks point at M-1 deck.
   *****************************************************************/

  best_score     = IMPOSSIBLE;
  best_neg_score = IMPOSSIBLE;
  L = j0-i0+1;
  if (W > L) W = L; 

  if(dsq == NULL)
    Die("in CYKScan, dsq is NULL\n");

  alpha = MallocOrDie (sizeof(float **) * cm->M);
  for (v = cm->M-1; v >= 0; v--) {	/* reverse, because we allocate E_M-1 first */
    if (cm->stid[v] == BEGL_S)
      {
	alpha[v] = MallocOrDie(sizeof(float *) * (W+1));
	for (j = 0; j <= W; j++)
	  alpha[v][j] = MallocOrDie(sizeof(float) * (W+1));
      }
    else if (cm->sttype[v] == E_st && v < cm->M-1) 
      alpha[v] = alpha[cm->M-1];
    else 
      {
	alpha[v] = MallocOrDie(sizeof(float *) * 2);
	for (j = 0; j < 2; j++) 
	  alpha[v][j] = MallocOrDie(sizeof(float) * (W+1));
      }
  }
  bestr = MallocOrDie(sizeof(int) * (W+1));

  /*****************************************************************
   * alpha initializations.
   * We initialize on d=0, subsequences of length 0; these are
   * j-independent. Any generating state (P,L,R) is impossible on d=0.
   * E=0 for d=0. B,S,D must be calculated. 
   * Also, for MP, d=1 is impossible.
   * Also, for E, all d>0 are impossible.
   *****************************************************************/ 
  for (v = cm->M-1; v >= 0; v--)
    {
      alpha[v][0][0] = IMPOSSIBLE;

      if      (cm->sttype[v] == E_st)  alpha[v][0][0] = 0;
      else if (cm->sttype[v] == MP_st) alpha[v][0][1] = alpha[v][1][1] = IMPOSSIBLE;
      else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) 
	{
	  y = cm->cfirst[v];
	  alpha[v][0][0] = cm->endsc[v];
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	    if ((sc = alpha[y+yoffset][0][0] + cm->tsc[v][yoffset]) > alpha[v][0][0]) 
	      alpha[v][0][0] = sc;
          /* ...we don't bother to look at local alignment starts here... */
	  bestr[0] = -1;
	  if (alpha[v][0][0] < IMPOSSIBLE) alpha[v][0][0] = IMPOSSIBLE;	
	}
      else if (cm->sttype[v] == B_st) 
	{
	  w = cm->cfirst[v];
	  y = cm->cnum[v];
	  alpha[v][0][0] = alpha[w][0][0] + alpha[y][0][0]; 
	}

      alpha[v][1][0] = alpha[v][0][0];
      if (cm->stid[v] == BEGL_S) 
	for (j = 2; j <= W; j++) 
	  alpha[v][j][0] = alpha[v][0][0];
    }
  for (d = 1; d <= W; d++)
    alpha[cm->M-1][0][d] = alpha[cm->M-1][1][d] = IMPOSSIBLE;

  /*****************************************************************
   * gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits.
   *****************************************************************/ 
  gamma    = MallocOrDie(sizeof(float) * (L+1));
  gamma[0] = 0; 
  gback    = MallocOrDie(sizeof(int)   * (L+1));
  gback[0] = -1;
  savesc   = MallocOrDie(sizeof(float) * (L+1));
  saver    = MallocOrDie(sizeof(int)   * (L+1));

  /*****************************************************************
   * The main loop: scan the sequence from position 1 to L.
   *****************************************************************/
  for (j = i0; j <= j0; j++) 
    {
      gamma_j = j-i0+1; /* j is actual index in j, gamma_j is offset j index in gamma* data structures */
      cur = j%2;
      prv = (j-1)%2;
      for (v = cm->M-1; v > 0; v--) /* ...almost to ROOT; we handle ROOT specially... */
	{
	  if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	    {
	      if (cm->stid[v] == BEGL_S) jp = j % (W+1); else jp = cur;
	      for (d = 1; d <= W && d <= gamma_j; d++) 
		{
		  y = cm->cfirst[v];
		  alpha[v][jp][d] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][cur][d] + cm->tsc[v][yoffset]) > alpha[v][jp][d]) 
		      alpha[v][jp][d] = sc;
		  if (alpha[v][jp][d] < IMPOSSIBLE) alpha[v][jp][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == MP_st) 
	    {
	      for (d = 2; d <= W && d <= gamma_j; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][prv][d-2] + cm->tsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  
		  i = j-d+1;
		  if (dsq[i] < Alphabet_size && dsq[j] < Alphabet_size)
		    alpha[v][cur][d] += cm->esc[v][(int) (dsq[i]*Alphabet_size+dsq[j])];
		  else
		    alpha[v][cur][d] += DegeneratePairScore(cm->esc[v], dsq[i], dsq[j]);
		  
		  if (alpha[v][cur][d] < IMPOSSIBLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) 
	    {
	      for (d = 1; d <= W && d <= gamma_j; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v]))); 
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][cur][d-1] + cm->tsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  
		  i = j-d+1;
		  if (dsq[i] < Alphabet_size)
		    alpha[v][cur][d] += cm->esc[v][(int) dsq[i]];
		  else
		    alpha[v][cur][d] += DegenerateSingletScore(cm->esc[v], dsq[i]);
		  
		  if (alpha[v][cur][d] < IMPOSSIBLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) 
	    {
	      for (d = 1; d <= W && d <= gamma_j; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][prv][d-1] + cm->tsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  
		  if (dsq[j] < Alphabet_size)
		    alpha[v][cur][d] += cm->esc[v][(int) dsq[j]];
		  else
		    alpha[v][cur][d] += DegenerateSingletScore(cm->esc[v], dsq[j]);
		  
		  if (alpha[v][cur][d] < IMPOSSIBLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == B_st) 
	    {
	      w = cm->cfirst[v];
	      y = cm->cnum[v];
	      i = j-d+1;
	      for (d = 1; d <= W && d <= gamma_j; d++) 
		{
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));

		  for (k = 0; k <= d; k++) /* k is length of right fragment */
		    {
		      jp = (j-k)%(W+1);	   /* jp is rolling index into BEGL_S deck j dimension */
		      if ((sc = alpha[w][jp][d-k] + alpha[y][cur][k]) > alpha[v][cur][d])
			alpha[v][cur][d] = sc;
		    }
		  if (alpha[v][cur][d] < IMPOSSIBLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	} /* end loop over decks v>0 */

      /* Finish up with the ROOT_S, state v=0; and deal w/ local begins.
       * 
       * If local begins are off, the hit must be rooted at v=0.
       * With local begins on, the hit is rooted at the second state in
       * the traceback (e.g. after 0), the internal entry point. Divide & conquer
       * can only handle this if it's a non-insert state; this is guaranteed
       * by the way local alignment is parameterized (other transitions are
       * -INFTY), which is probably a little too fragile of a method. 
       */
      for (d = 1; d <= W && d <= gamma_j; d++) 
	{
	  y = cm->cfirst[0];
	  alpha[0][cur][d] = alpha[y][cur][d] + cm->tsc[0][0];
	  bestr[d]         = 0;	/* root of the traceback = root state 0 */
	  for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++)
	    if ((sc = alpha[y+yoffset][cur][d] + cm->tsc[0][yoffset]) > alpha[0][cur][d]) 
	      alpha[0][cur][d] = sc;

	  if (cm->flags & CM_LOCAL_BEGIN) {
	    for (y = 1; y < cm->M; y++) {
	      if (cm->stid[y] == BEGL_S) sc = alpha[y][j%(W+1)][d] + cm->beginsc[y];
	      else                       sc = alpha[y][cur][d]     + cm->beginsc[y];
	      if (sc > alpha[0][cur][d]) {
		alpha[0][cur][d] = sc;
		bestr[d]         = y;
	      }
	    }
	  }
	  if (alpha[0][cur][d] < IMPOSSIBLE) alpha[0][cur][d] = IMPOSSIBLE;
	  if (alpha[0][cur][d] > best_neg_score) best_neg_score = alpha[0][cur][d];
	}

      if(!(cm->search_opts & CM_SEARCH_GREEDY)) /* resolve overlaps optimally */
	{
	  /* The little semi-Markov model that deals with multihit parsing:
	   */
	  gamma[gamma_j]  = gamma[gamma_j-1] + 0; 
	  gback[gamma_j]  = -1;
	  savesc[gamma_j] = IMPOSSIBLE;
	  saver[gamma_j]  = -1;
	  for (d = 1; d <= W && d <= gamma_j; d++) 
	    {
	      i = j-d+1;
	      gamma_i = j-d+1-i0+1;
	      sc = gamma[gamma_i-1] + alpha[0][cur][d] + cm->sc_boost; 
	      /* sc_boost is experimental technique for finding hits < 0 bits. 
	       * value is 0.0 if technique not used. */
	      if (sc > gamma[gamma_j])
		{
		  gamma[gamma_j]  = sc;
		  gback[gamma_j]  = i;
		  savesc[gamma_j] = alpha[0][cur][d]; 
		  saver[gamma_j]  = bestr[d];
		}
	    }
	}
      else
	{
	  /* Resolving overlaps greedily (RSEARCH style),  
	   * At least one hit is sent back for each j here.
	   * However, some hits can already be removed for the greedy overlap
	   * resolution algorithm.  Specifically, at the given j, any hit with a
	   * d of d1 is guaranteed to mask any hit of lesser score with a d > d1 */
	  /* First, report hit with d of 1 if > cutoff */
	  if (alpha[0][cur][1] >= cutoff) 
	    if(results != NULL) 
	      {
		report_hit (j, j, bestr[1], alpha[0][cur][1], results);
	      }
	  bestd = 1;
	  if (alpha[0][cur][1] > best_score) 
	    best_score = alpha[0][cur][1];
	  
	  /* Now, if current score is greater than maximum seen previous, report
	     it if >= cutoff and set new max */
	  for (d = 2; d <= W && d <= gamma_j; d++) 
	    {
	      if (alpha[0][cur][d] > best_score) 
		best_score = alpha[0][cur][d];
	      if (alpha[0][cur][d] > alpha[0][cur][bestd]) 
		{
		  if (alpha[0][cur][d] >= cutoff)
		    if(results != NULL) 
		      {
			report_hit (j-d+1, j, bestr[d], alpha[0][cur][d], results);
		      }
		  bestd = d;
		}
	    }
	}
    } /* end loop over end positions j */
  
  /*****************************************************************
   * we're done with alpha, free it; everything we need is in gamma.
   *****************************************************************/ 
  for (v = 0; v < cm->M; v++) 
    {
      if (cm->stid[v] == BEGL_S) {                     /* big BEGL_S decks */
	for (j = 0; j <= W; j++) free(alpha[v][j]);
	free(alpha[v]);
      } else if (cm->sttype[v] == E_st && v < cm->M-1) { /* avoid shared E decks */
	continue;
      } else {
	free(alpha[v][0]);
	free(alpha[v][1]);
	free(alpha[v]);
      }
    }
  free(alpha);
  free(bestr);

  /*****************************************************************
   * Traceback stage.
   * Recover all hits: an (i,j,sc) triple for each one.
   *****************************************************************/ 
  if(!(cm->search_opts & CM_SEARCH_GREEDY)) /* resolve overlaps optimally */
    {
      j     = j0;
      while (j >= i0) 
	{
	  gamma_j = j-i0+1;
	  if (gback[gamma_j] == -1) /* no hit */
	    j--; 
	  else                /* a hit, a palpable hit */
	    {
	      if(savesc[gamma_j] > best_score) 
		best_score = savesc[gamma_j];
	      if(savesc[gamma_j] >= cutoff && results != NULL) /* report the hit */
		report_hit(gback[gamma_j], j, saver[gamma_j], savesc[gamma_j], results);
	      j = gback[gamma_j]-1;
	    }
	}
    }
  free(gback);
  free(gamma);
  free(savesc);
  free(saver);
  if(best_score <= 0.) /* there were no hits found by the semi-HMM, no hits above 0 bits */
    best_score = best_neg_score;

  return best_score;
}

/* Function: CYKScanRequires()
 * Date:     SRE, Mon May  6 18:48:18 2002 [St. Louis]
 *
 * Purpose:  Return the memory required by CYKScan(), in megabytes.
 */
float
CYKScanRequires(CM_t *cm, int L, int W)
{
  float ram;
  int   v,j;

				/* alpha allocations */
  ram = (float) (sizeof(float **) * cm->M);
  for (v = cm->M-1; v >= 0; v--)
    {
      if (cm->stid[v] == BEGL_S)
	{
	  ram += (float) (sizeof(float *) * (W+1));
	  for (j = 0; j <= W; j++)
	    ram += (float) (sizeof(float) * (W+1));
	}
      else if (cm->sttype[v] == E_st && v < cm->M-1) 
	continue;
      else 
	{
	  ram += (float) (sizeof(float *) * 2);
	  for (j = 0; j < 2; j++) 
	    ram += (float) (sizeof(float) * (W+1));
	}
    }
  ram += (float) (sizeof(float) * (L+1)); /* gamma allocation */
  ram += (float) (sizeof(int)   * (L+1)); /* gback allocation */
  ram += (float) (sizeof(float) * (L+1)); /* savesc allocation */
  return (ram / 1000000.);
}

/* Function: CYKScanRFWrapper()
 * Date:     EPN, Thu May  3 15:48:03 2007
 *
 * Purpose:  Scan a sequence for matches to a covariance model.
 *           Multiple nonoverlapping hits and local alignment.
 *
 *           Called by Zasha Weinberg's rigorous filter code.
 *           Created solely to maintain compatibility with rigfilters.
 *           This function has the same prototype as the old
 *           (0.72 and before) CYKScan() that was called
 *           by the original rigfilters code. We call
 *           the new function, and convert the results to 
 *           the old format.
 *
 * Args:     cm        - the covariance model
 *           dsq       - digitized sequence to search; i0..j0
 *           i0        - start of target subsequence (1 for full seq)
 *           j0        - end of target subsequence (L for full seq)
 *           W         - max d: max size of a hit
 *           ret_nhits - RETURN: number of hits
 *           ret_hitr  - RETURN: start states of hits, 0..nhits-1
 *           ret_hiti  - RETURN: start positions of hits, 0..nhits-1
 *           ret_hitj  - RETURN: end positions of hits, 0..nhits-1
 *           ret_hitsc - RETURN: scores of hits, 0..nhits-1            
 *           min_thresh- minimum score to report (EPN via Alex Coventry 03.11.06)
 *
 * Returns:  
 *           hiti, hitj, hitsc are allocated here; caller free's w/ free().
 */
void
CYKScanRFWrapper(CM_t *cm, char *dsq, int i0, int j0, int W, 
		 int *ret_nhits, int **ret_hitr, int **ret_hiti, int **ret_hitj, float **ret_hitsc, 
		 float min_thresh)
{
  int       nhits;		/* # of hits in optimal parse */
  int      *hitr;		/* initial state indices of hits in optimal parse */
  int      *hiti;               /* start positions of hits in optimal parse */
  int      *hitj;               /* end positions of hits in optimal parse */
  float    *hitsc;              /* scores of hits in optimal parse */
  scan_results_t *results;      /* the results data structure, we extract hits from */
  int       h;                  /* counter over hits */

  hitr = NULL;
  hiti = NULL;
  hitj = NULL;
  hitsc = NULL;
  results = CreateResults(INIT_RESULTS);
  CYKScan(cm, dsq, i0, j0, W, 0.0, /* cutoff, in 0.72 and before, implicitly 0.0 */
	  results);
  nhits = results->num_results;
  if(nhits > 0)
    {
      hitr  = MallocOrDie(sizeof(int)   * nhits);
      hitj  = MallocOrDie(sizeof(int)   * nhits);
      hiti  = MallocOrDie(sizeof(int)   * nhits);
      hitsc = MallocOrDie(sizeof(float) * nhits);
      for(h = 0; h < nhits; h++)
	{
	  hitr[h]  = results->data[h].bestr;
	  hiti[h]  = results->data[h].start;
	  hitj[h]  = results->data[h].stop;
	  hitsc[h] = results->data[h].score;
	}
    }  
  *ret_nhits = nhits;
  *ret_hitr  = hitr;
  *ret_hiti  = hiti;
  *ret_hitj  = hitj;
  *ret_hitsc = hitsc;
  return;
}

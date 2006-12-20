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
 * Function: sort_results()
 * Date:    RJK,  Sun Mar 31, 2002 [AA Flight 2869 LGA->STL]
 * Purpose: Given a results array, sorts it with a call to qsort
 *
 */
void sort_results (scan_results_t *results) {
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
void report_hit (int i, int j, int bestr, float score, scan_results_t *results) {
  if (results->num_results == results->num_allocated) {
    ExpandResults (results, INIT_RESULTS);
  }
  results->data[results->num_results].score = score;
  results->data[results->num_results].start = i;
  results->data[results->num_results].stop = j;
  results->data[results->num_results].bestr = bestr;
  results->data[results->num_results].tr = NULL;
  results->num_results++;
}

/*
 * Function: remove_overlapping_hits ()
 * Date:     RJK, Sun Mar 31, 2002 [LGA Gate D7]
 *
 * Purpose:  Given a list of hits, removes overlapping hits to produce
 * a list consisting of at most one hit covering each nucleotide in the
 * sequence.  Works as follows:
 * 1.  Bubble sort hits (I know this is ridiculously slow, but I want
 *     to get something working.  I can replace this later with a faster
 *     algorithm.)  (O(N^2))
 * 2.  For each hit, sees if any nucleotide covered yet
 *     If yes, remove hit
 *     If no, mark each nt as covered
 */
void remove_overlapping_hits (scan_results_t *results, int L) {
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
 *           score_boost - boost in bits to temporarily add to all scores, 
 *                         experimental technique for finding significant 
 *                         hits < 0 bits. 0.0 if technique not used.
 *           results     - scan_results_t to add to; if NULL, don't add to it
 *
 * Returns:  score of best overall hit
 */
float 
CYKScan(CM_t *cm, char *dsq, int i0, int j0, int W, 
	float cutoff, float score_boost, scan_results_t *results)
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
  float     best_score;         /* Best overall score to return */
  /*int     updated_flag;*/         /* strategy 2 */
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

  best_score = IMPOSSIBLE;
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
      gamma_j = j-i0+1; /* j is actual index in j, gamma_j is offeset j index in gamma* data structures */
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
	  if (alpha[0][cur][d] > best_score) best_score = alpha[0][cur][d];
	}

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
	  sc = gamma[gamma_i-1] + alpha[0][cur][d] + score_boost; 
	  if (sc > gamma[gamma_j])
	    {
	      gamma[gamma_j]  = sc;
	      gback[gamma_j]  = i;
	      savesc[gamma_j] = alpha[0][cur][d]; 
	      saver[gamma_j]  = bestr[d];
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
  if(results != NULL)
    {
      j     = j0;
      while (j >= i0) 
	{
	  gamma_j = j-i0+1;
	  if (gback[gamma_j] == -1) /* no hit */
	    j--; 
	  else                /* a hit, a palpable hit */
	    {
	      if(savesc[gamma_j] >= cutoff) /* report the hit */
		report_hit(gback[gamma_j], j, saver[gamma_j], savesc[gamma_j], results);
	      j = gback[gamma_j]-1;
	    }
	}
    }
  free(gback);
  free(gamma);
  free(savesc);
  free(saver);

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

/* Function: iCYKScan()
 * Date:     EPN, Tue Dec 19 13:54:26 2006
 *           based on CYKScan() SRE, Thu May  2 11:56:11 2002 [AA 3050 SFO->STL]
 *
 * Purpose:  Scan a sequence for matches to a covariance model.
 *           Multiple nonoverlapping hits and local alignment.
 *           Diff with CYKScan(): here scaled int log odds scores are used
 *           instead of float log odds scores.
 *
 * Args:     cm          - the covariance model
 *           dsq         - digitized sequence to search; i0..j0
 *           i0          - start of target subsequence (1 for full seq)
 *           j0          - end of target subsequence (L for full seq)
 *           W           - max d: max size of a hit
 *           cutoff      - minimum score to report 
 *           score_boost - boost in bits to temporarily add to all scores, 
 *                         experimental technique for finding significant 
 *                         hits < 0 bits. 0.0 if technique not used.
 *           results     - scan_results_t to add to; if NULL, don't add to it
 *
 * Returns:  score of best overall hit
 */
float 
iCYKScan(CM_t *cm, char *dsq, int i0, int j0, int W, 
	 float cutoff, float score_boost, scan_results_t *results)
{
  int    ***alpha;              /* CYK DP score matrix, [v][j][d] */
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
  float     best_score;         /* Best overall score to return */
  /*int     updated_flag;*/         /* strategy 2 */
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

  best_score = IMPOSSIBLE;
  L = j0-i0+1;
  if (W > L) W = L; 

  printf("in iCYKScan\n");
  if(dsq == NULL)
    Die("in CYKScan, dsq is NULL\n");

  alpha = MallocOrDie (sizeof(int **) * cm->M);
  for (v = cm->M-1; v >= 0; v--) {	/* reverse, because we allocate E_M-1 first */
    if (cm->stid[v] == BEGL_S)
      {
	alpha[v] = MallocOrDie(sizeof(int *) * (W+1));
	for (j = 0; j <= W; j++)
	  alpha[v][j] = MallocOrDie(sizeof(int) * (W+1));
      }
    else if (cm->sttype[v] == E_st && v < cm->M-1) 
      alpha[v] = alpha[cm->M-1];
    else 
      {
	alpha[v] = MallocOrDie(sizeof(int *) * 2);
	for (j = 0; j < 2; j++) 
	  alpha[v][j] = MallocOrDie(sizeof(int) * (W+1));
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
      alpha[v][0][0] = -INFTY;

      if      (cm->sttype[v] == E_st)  alpha[v][0][0] = 0;
      else if (cm->sttype[v] == MP_st) alpha[v][0][1] = alpha[v][1][1] = -INFTY;
      else if (cm->sttype[v] == S_st || cm->sttype[v] == D_st) 
	{
	  y = cm->cfirst[v];
	  alpha[v][0][0] = cm->iendsc[v];
	  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
	    if ((sc = alpha[y+yoffset][0][0] + cm->itsc[v][yoffset]) > alpha[v][0][0]) 
	      alpha[v][0][0] = sc;
          /* ...we don't bother to look at local alignment starts here... */
	  bestr[0] = -1;
	  if (alpha[v][0][0] < -INFTY) alpha[v][0][0] = -INFTY;	
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
    alpha[cm->M-1][0][d] = alpha[cm->M-1][1][d] = -INFTY;

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
      gamma_j = j-i0+1; /* j is actual index in j, gamma_j is offeset j index in gamma* data structures */
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
		  alpha[v][jp][d] = cm->iendsc[v] + (cm->iel_selfsc * (d - StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][cur][d] + cm->itsc[v][yoffset]) > alpha[v][jp][d]) 
		      alpha[v][jp][d] = sc;
		  if (alpha[v][jp][d] < -INFTY) alpha[v][jp][d] = -INFTY;
		}
	    }
	  else if (cm->sttype[v] == MP_st) 
	    {
	      for (d = 2; d <= W && d <= gamma_j; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->iendsc[v] + (cm->iel_selfsc * (d - StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][prv][d-2] + cm->itsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  
		  i = j-d+1;
		  if (dsq[i] < Alphabet_size && dsq[j] < Alphabet_size)
		    alpha[v][cur][d] += cm->iesc[v][(int) (dsq[i]*Alphabet_size+dsq[j])];
		  else
		    alpha[v][cur][d] += iDegeneratePairScore(cm->iesc[v], dsq[i], dsq[j]);
		  
		  if (alpha[v][cur][d] < -INFTY) alpha[v][cur][d] = -INFTY;
		}
	    }
	  else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) 
	    {
	      for (d = 1; d <= W && d <= gamma_j; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->iendsc[v] + (cm->iel_selfsc * (d - StateDelta(cm->sttype[v]))); 
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][cur][d-1] + cm->itsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  
		  i = j-d+1;
		  if (dsq[i] < Alphabet_size)
		    alpha[v][cur][d] += cm->iesc[v][(int) dsq[i]];
		  else
		    alpha[v][cur][d] += iDegenerateSingletScore(cm->iesc[v], dsq[i]);
		  
		  if (alpha[v][cur][d] < -INFTY) alpha[v][cur][d] = -INFTY;
		}
	    }
	  else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) 
	    {
	      for (d = 1; d <= W && d <= gamma_j; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->iendsc[v] + (cm->iel_selfsc * (d - StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][prv][d-1] + cm->itsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  
		  if (dsq[j] < Alphabet_size)
		    alpha[v][cur][d] += cm->iesc[v][(int) dsq[j]];
		  else
		    alpha[v][cur][d] += iDegenerateSingletScore(cm->iesc[v], dsq[j]);
		  
		  if (alpha[v][cur][d] < -INFTY) alpha[v][cur][d] = -INFTY;
		}
	    }
	  else if (cm->sttype[v] == B_st) 
	    {
	      w = cm->cfirst[v];
	      y = cm->cnum[v];
	      i = j-d+1;
	      for (d = 1; d <= W && d <= gamma_j; d++) 
		{
		  alpha[v][cur][d] = cm->iendsc[v] + (cm->iel_selfsc * (d - StateDelta(cm->sttype[v])));

		  for (k = 0; k <= d; k++) /* k is length of right fragment */
		    {
		      jp = (j-k)%(W+1);	   /* jp is rolling index into BEGL_S deck j dimension */
		      if ((sc = alpha[w][jp][d-k] + alpha[y][cur][k]) > alpha[v][cur][d])
			alpha[v][cur][d] = sc;
		    }
		  if (alpha[v][cur][d] < -INFTY) alpha[v][cur][d] = -INFTY;
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
	  alpha[0][cur][d] = alpha[y][cur][d] + cm->itsc[0][0];
	  bestr[d]         = 0;	/* root of the traceback = root state 0 */
	  for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++)
	    if ((sc = alpha[y+yoffset][cur][d] + cm->itsc[0][yoffset]) > alpha[0][cur][d]) 
	      alpha[0][cur][d] = sc;

	  if (cm->flags & CM_LOCAL_BEGIN) {
	    for (y = 1; y < cm->M; y++) {
	      if (cm->stid[y] == BEGL_S) sc = alpha[y][j%(W+1)][d] + cm->ibeginsc[y];
	      else                       sc = alpha[y][cur][d]     + cm->ibeginsc[y];
	      if (sc > alpha[0][cur][d]) {
		alpha[0][cur][d] = sc;
		bestr[d]         = y;
	      }
	    }
	  }
	  if (alpha[0][cur][d] < -INFTY) alpha[0][cur][d] = -INFTY;
	  if (Scorify(alpha[0][cur][d]) > best_score) best_score = Scorify(alpha[0][cur][d]);
	}

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
	  sc = gamma[gamma_i-1] + Scorify(alpha[0][cur][d]) + score_boost; 
	  if (sc > gamma[gamma_j])
	    {
	      gamma[gamma_j]  = sc;
	      gback[gamma_j]  = i;
	      savesc[gamma_j] = Scorify(alpha[0][cur][d]); 
	      saver[gamma_j]  = bestr[d];
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
  if(results != NULL)
    {
      j     = j0;
      while (j >= i0) 
	{
	  gamma_j = j-i0+1;
	  if (gback[gamma_j] == -1) /* no hit */
	    j--; 
	  else                /* a hit, a palpable hit */
	    {
	      if(savesc[gamma_j] >= cutoff) /* report the hit */
		report_hit(gback[gamma_j], j, saver[gamma_j], savesc[gamma_j], results);
	      j = gback[gamma_j]-1;
	    }
	}
    }
  free(gback);
  free(gamma);
  free(savesc);
  free(saver);

  return best_score;
}

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

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_sqio.h"
#include "esl_vectorops.h"

#include "funcs.h"
#include "structs.h"

/*
 * Function: CreateResults ()
 * Date:     RJK, Mon Apr 1 2002 [St. Louis]
 * Purpose:  Creates a results type of specified size
 */
search_results_t *CreateResults (int size) {
  int status;
  search_results_t *results;
  int i;

  if(size == 0) return NULL;

  ESL_ALLOC(results, sizeof(search_results_t));
  results->num_results = 0;
  results->num_allocated = size;
  ESL_ALLOC(results->data, sizeof(search_result_node_t)*size);
  for(i = 0; i < size; i++) {
    results->data[i].start  = -1;    
    results->data[i].stop   = -1;    
    results->data[i].bestr  = -1;    
    results->data[i].score  = IMPOSSIBLE;    
    results->data[i].tr     = NULL;
    results->data[i].pcode1 = NULL;
    results->data[i].pcode2 = NULL;
  }
  return (results);
 ERROR:
  esl_fatal("Memory allocation error.");
  return NULL; /* never reached */
}

/* Function: ExpandResults ()
 * Date:     RJK, Mon Apr 1 2002 [St. Louis]
 * Purpose:  Expands a results structure by specified amount
 */
void ExpandResults (search_results_t *results, int additional) {
  int status;
  void *tmp;
  int i;
  ESL_RALLOC(results->data, tmp, sizeof(search_result_node_t) * (results->num_allocated+additional));

  for(i = results->num_allocated; i < (results->num_allocated + additional); i++) {
    results->data[i].start  = -1;    
    results->data[i].stop   = -1;    
    results->data[i].bestr  = -1;    
    results->data[i].score  = IMPOSSIBLE;    
    results->data[i].tr     = NULL;
    results->data[i].pcode1 = NULL;
    results->data[i].pcode2 = NULL;
  }

  results->num_allocated+=additional;
  return;
 ERROR:
  esl_fatal("Memory reallocation error.");
}

/* Function: AppendResults()
 * Date:     EPN, Wed Aug 29 08:58:28 2007
 *
 * Purpose:  Add result nodes from one results structure onto
 *           another by copying data and manipulating pointers. 
 *           Originally written to add results returned from 
 *           an MPI worker to a growing 'master' results structure 
 *           in the MPI master. 

 *           The search_results_node_t's (results->data)
 *           must have their start, stop, bestr, and score data
 *           copied because the results->data is not a set 
 *           of pointers (the whole reason for this is so 
 *           we can call quicksort() on the hits, which requires
 *           we have a fixed distance between them in memory).
 *           The parsetree, however, can have it's pointer 
 *           simply switched.
 *
 *           i0 is an offset in the sequence, (i0-1) is added to
 *           data[i].start and data[i].stop for hits in src_results. 
 *           i0 == 1 means no offset. i0 != 1 usually useful for hits 
 *           returned by MPI workers who were searching a database
 *           subsequence.
 *
 * Note:     Because the dest_results->data now points to some of the
 *           src_results->data, be careful not to free src_results with
 *           FreeResults() if you don't want to lose dest_results.
 */
void AppendResults (search_results_t *src_results, search_results_t *dest_results, int i0) {
  int i, ip;
  for(i = 0; i < src_results->num_results; i++) 
    {
      ip = dest_results->num_results;
      report_hit (src_results->data[i].start+i0-1, src_results->data[i].stop+i0-1, 
		  src_results->data[i].bestr,      src_results->data[i].score,
		  dest_results);
      if(src_results->data[i].tr != NULL)
	(*dest_results).data[ip].tr = (*src_results).data[i].tr;
      if(src_results->data[i].pcode1 != NULL)
	(*dest_results).data[ip].pcode1 = (*src_results).data[i].pcode1;
      if(src_results->data[i].pcode2 != NULL)
	(*dest_results).data[ip].pcode2 = (*src_results).data[i].pcode2;
    }
  return;
}

/*
 * Function: FreeResults()
 * Date:     RJK, Mon Apr 1 2002 [St. Louis]
 * Purpose:  Frees a results structure
 */
void FreeResults (search_results_t *r) {
  int i;
  if (r != NULL) {
    for (i=0; i < r->num_allocated; i++) {
      if (r->data[i].tr != NULL)     FreeParsetree(r->data[i].tr);
      if (r->data[i].pcode1 != NULL) free(r->data[i].pcode1);
      if (r->data[i].pcode2 != NULL) free(r->data[i].pcode2);
    }
    free (r->data);
    free(r);
  }
}


/*
 * Function: compare_results()
 * Date:     RJK, Wed Apr 10, 2002 [St. Louis]
 * Purpose:  Compares two search_result_node_ts based on score and returns -1
 *           if first is higher score than second, 0 if equal, 1 if first
 *           score is lower.  This results in sorting by score, highest
 *           first.
 */
int compare_results (const void *a_void, const void *b_void) {
  search_result_node_t *a, *b;
 
  a = (search_result_node_t *)a_void;
  b = (search_result_node_t *)b_void;

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
 * Plucked verbatim out of RSEARCH.
 * RSEARCH date: RJK, Sun Mar 31, 2002 [LGA Gate D7]
 *
 * Purpose:  Given a list of hits, removes overlapping hits to produce
 * a list consisting of at most one hit covering each nucleotide in the
 * sequence.  Works as follows:
 * 1.  quicksort hits 
 * 2.  For each hit, sees if any nucleotide covered yet
 *     If yes, remove hit
 *     If no, mark each nt as covered
 * 
 * Args: 
 *        i0    - first position of subsequence results work for (1 for full seq)
 *        j0    - last  position of subsequence results work for (L for full seq)
 */
void remove_overlapping_hits (search_results_t *results, int i0, int j0)
{
  int status;
  char *covered_yet;
  int x,y;
  int covered;
  int L;
  int yp;          /* offset position, yp = y-i0+1 */
  search_result_node_t swap;

  if (results == NULL)
    return;

  if (results->num_results == 0)
    return;

  L = j0-i0+1;
  ESL_ALLOC(covered_yet, sizeof(char)*(L+1));
  for (x=0; x<=L; x++)
    covered_yet[x] = 0;

  sort_results (results);

  for (x=0; x<results->num_results; x++) {
    covered = 0;
    for (y=results->data[x].start; y<=results->data[x].stop && !covered; y++) {
      {
	yp = y-i0+1; 
	assert(yp > 0 && yp <= L);
	if (covered_yet[yp] != 0) {
	  covered = 1;
	} 
      }
    }
    if (covered == 1) {
      results->data[x].start = -1;        /* Flag -- remove later to keep sorted */
    } else {
      for (y=results->data[x].start; y<=results->data[x].stop; y++) {
	yp = y-i0+1; 
	covered_yet[yp] = 1;
      }
    }
  }
  free (covered_yet);

  for (x=0; x < results->num_results; x++) {
    while (results->num_results > 0 &&
	   results->data[results->num_results-1].start == -1)
      results->num_results--;
    if (x < results->num_results && results->data[x].start == -1) {
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
  return;

 ERROR:
  esl_fatal("Memory allocation error.");
}

/*
 * Function: sort_results()
 * Date:    RJK,  Sun Mar 31, 2002 [AA Flight 2869 LGA->STL]
 * Purpose: Given a results array, sorts it with a call to qsort
 *
 */
void sort_results (search_results_t *results) 
{
  qsort (results->data, results->num_results, sizeof(search_result_node_t), compare_results);
}

/*
 * Function: report_hit()
 * Date:     RJK, Sun Mar 31, 2002 [LGA Gate D7]
 *
 * Given j,d, coordinates, a score, and a search_results_t data type,
 * adds result into the set of reportable results.  Naively adds hit.
 *
 * Non-overlap algorithm is now done in the scanning routine by Sean's
 * Semi-HMM code.  I've just kept the hit report structure for convenience.
 */
void report_hit (int i, int j, int bestr, float score, search_results_t *results) 
{

  if(results == NULL) 
    esl_fatal("in report_hit, but results is NULL\n");
  if (results->num_results == results->num_allocated) 
    ExpandResults (results, INIT_RESULTS);

  results->data[results->num_results].score = score;
  results->data[results->num_results].start = i;
  results->data[results->num_results].stop = j;
  results->data[results->num_results].bestr = bestr;
  results->data[results->num_results].tr = NULL;
  results->data[results->num_results].pcode1 = NULL;
  results->data[results->num_results].pcode2 = NULL;
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
 *           results     - search_results_t to add to; if NULL, don't add to it
 *
 * Returns:  score of best overall hit
 */
float 
CYKScan(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, 
	float cutoff, search_results_t *results)
{
  /* Contract check */
  if(j0 < i0)
    esl_fatal("ERROR in CYKScan, i0: %d j0: %d\n", i0, j0);
  if(dsq == NULL)
    esl_fatal("ERROR in CYKScan, dsq is NULL\n");

  int       status;
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

  best_score     = IMPOSSIBLE;
  best_neg_score = IMPOSSIBLE;
  L = j0-i0+1;
  if (W > L) W = L; 

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

  ESL_ALLOC(alpha, (sizeof(float **) * cm->M));
  for (v = cm->M-1; v >= 0; v--) {	/* reverse, because we allocate E_M-1 first */
    if (cm->stid[v] == BEGL_S)
      {
	ESL_ALLOC(alpha[v], (sizeof(float *) * (W+1)));
	for (j = 0; j <= W; j++)
	  ESL_ALLOC(alpha[v][j], (sizeof(float) * (W+1)));
      }
    else if (cm->sttype[v] == E_st && v < cm->M-1) 
      alpha[v] = alpha[cm->M-1];
    else 
      {
	ESL_ALLOC(alpha[v], sizeof(float *) * 2);
	for (j = 0; j < 2; j++) 
	  ESL_ALLOC(alpha[v][j], (sizeof(float) * (W+1)));
      }
  }
  ESL_ALLOC(bestr, (sizeof(int) * (W+1)));

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
  ESL_ALLOC(gamma,  sizeof(float) * (L+1));
  gamma[0] = 0;
  ESL_ALLOC(gback,  sizeof(int)   * (L+1));
  gback[0] = -1;
  ESL_ALLOC(savesc, sizeof(float) * (L+1));
  ESL_ALLOC(saver,  sizeof(int)   * (L+1));

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
		  if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][(dsq[i]*cm->abc->K+dsq[j])];
		  else
		    alpha[v][cur][d] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
		  
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
		  if (dsq[i] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][dsq[i]];
		  else
		    alpha[v][cur][d] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
		  
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
		  if (dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][dsq[j]];
		  else
		    alpha[v][cur][d] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
		  
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
	  if (cm->flags & CMH_LOCAL_BEGIN) {
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

      if(!(cm->search_opts & CM_SEARCH_CMGREEDY)) /* resolve overlaps optimally */
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
	      sc = gamma[gamma_i-1] + alpha[0][cur][d];
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
  if(!(cm->search_opts & CM_SEARCH_CMGREEDY)) /* resolve overlaps optimally */
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
  /* printf("CYKScan() return score: %10.4f\n", best_score); */
  return best_score;

 ERROR:
  esl_fatal("Memory allocation error.\n");
  return 0.; /* never reached */
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


/* Function: CountScanDPCalcs()
 * Date:     EPN, Wed Aug 22 09:08:03 2007
 *
 * Purpose:  Count all DP calcs for a CM scan against a 
 *           sequence of length L. Similar to smallcyk.c's
 *           CYKDemands() but takes into account number of
 *           transitions from each state, and is concerned
 *           with a scanning dp matrix, not an alignment matrix.
 *
 * Args:     cm     - the model
 *           L      - length of sequence
 *           use_qdb- TRUE to enforce cm->dmin and cm->dmax for calculation
 *
 * Returns: (float) the total number of DP calculations, either using QDB or not.
 */
float
CountScanDPCalcs(CM_t *cm, int L, int use_qdb)
{
  int v, j;
  float dpcalcs = 0.;
  float dpcalcs_bif = 0.;
  
  /* Contract check */
  if(cm->W > L) esl_fatal("ERROR in CountScanDPCalcs(), cm->W: %d exceeds L: %d\n", cm->W, L);

  float  dpcells     = 0.;
  float  dpcells_bif = 0.;
  int d,w,y,kmin,kmax, bw;

  if(! use_qdb) 
    {
      dpcells = (cm->W * L) - (cm->W * (cm->W-1) * .5); /* fillable dp cells per state (deck) */
      for (j = 1; j < cm->W; j++)
	dpcells_bif += ((j    +2) * (j    +1) * .5) - 1;
      for (j = cm->W; j <= L; j++)
	dpcells_bif += ((cm->W+2) * (cm->W+1) * .5) - 1; 
      dpcalcs_bif = CMCountStatetype(cm, B_st) * dpcells_bif; /* no choice of transition */
      for(v = 0; v < cm->M; v++)
	{
	  if(cm->sttype[v] != B_st && cm->sttype[v] != E_st)
	    {
	      dpcalcs += dpcells * cm->cnum[v]; /* cnum choices of transitions */

	      /* non-obvious subtractions that match implementation in scancyk.c::CYKScan() */
	      if(v == 0) dpcalcs  -= dpcells; 
	      if(cm->sttype[v] == MP_st) dpcalcs  -= L * cm->cnum[v];
	    }
	}
    }
  else /* use_qdb */
    {
      for(v = 0; v < cm->M; v++)
	{
	  if(cm->sttype[v] != B_st && cm->sttype[v] != E_st)
	    {
	      bw = cm->dmax[v] - cm->dmin[v] + 1; /* band width */
	      if(cm->dmin[v] == 0) bw--;
	      dpcalcs += ((bw * L) - (bw * (bw-1) * 0.5)) * cm->cnum[v];

	      /* non-obvious subtractions that match implementation in bandcyk.c::CYKBandedScan() */
	      if(v == 0) dpcalcs  -= ((bw * L) - (bw * (bw-1) * 0.5)); 
	      if(cm->sttype[v] == MP_st) dpcalcs  -= bw * cm->cnum[v];
	    }

	  else if(cm->sttype[v] == B_st)
	    {
	      w = cm->cfirst[v];
	      y = cm->cnum[v];
	      for (j = 1; j <= L; j++)
		{
		  d = (cm->dmin[v] > 0) ? cm->dmin[v] : 1;
		  for (; d <= cm->dmax[v] && d <= j; d++)
		    {
		      if(cm->dmin[y] > (d-cm->dmax[w])) kmin = cm->dmin[y];
		      else kmin = d-cm->dmax[w];
		      if(kmin < 0) kmin = 0;
		      if(cm->dmax[y] < (d-cm->dmin[w])) kmax = cm->dmax[y];
		      else kmax = d-cm->dmin[w];
		      if(kmin <= kmax)
			{
			  bw = (kmax - kmin + 1);
			  dpcalcs_bif += bw;
			}
		    }
		}
	    }
	}
    }
  /*printf("%d CountScanDPCalcs dpc     %.0f\n", use_qdb, dpcalcs);
    printf("%d CountScanDPCalcs dpc_bif %.0f\n", use_qdb, dpcalcs_bif);
    printf("%d CountScanDPCalcs total   %.0f\n", use_qdb, dpcalcs + dpcalcs_bif);*/
  return dpcalcs + dpcalcs_bif;
}

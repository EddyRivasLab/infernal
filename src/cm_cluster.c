/************************************************************
 * @LICENSE@
 ************************************************************/

/* cm_cluster.c
 * EPN, Wed Mar 21 17:25:21 2007
 * 
 * Functions to support building multiple CMs from 
 * a single input MSA in cmbuild, with each CM built from
 * a cluster of seqs in the MSA.
 * 
 */

#include "config.h"
#include "squidconf.h"

#include <stdio.h>
#include <math.h>
#include <float.h>

#include "easel.h"		
#include "esl_tree.h"
#include "esl_dmatrix.h"
#include "esl_stack.h"
#include "esl_distance.h"
#include "esl_vectorops.h"
#include "squid.h"
#include "structs.h"
#include "funcs.h"

static int select_node(ESL_TREE *T, double *diff, double mindiff, 
		       int **ret_clust, int *ret_nc);
static float find_mindiff(ESL_TREE *T, double *diff, int target_nc, 
			  int **ret_clust, int *ret_nc);

/* Function: MSADivide()
 * EPN, Wed Mar 21 17:26:39 2007
 * 
 * Purpose:  Given an MSA, divide it into multiple MSAs, each with
 *           a different cluster of the original sequences. Each
 *           MSA will be used to construct a separate CM.
 *
 *           Different modes:
 *           
 *        1. if(do_all): each seq is its own cluster, so
 *           the number of new MSAs is number of seqs in input 
 *           master MSA. 
 *
 *        2. if(!do_all) && target_nc == 0: define clusters
 *           such that we maximize the number of clusters while
 *           satisfying: minimum fractional difference b/t any 
 *           2 seqs in different clusters >= 'mindiff'. 
 *           The contract states that mindiff > 0. in this case.
 *           
 *        3. if(!do_all) && mindiff == 0.: define clusters 
 *           such that we have exactly 'target_nc' clusters by
 *           searching for the 'mindiff' that gives exactly
 *           'target_nc' clusters. (We guarantee we can do this
 *           by rounding 'diff' fractional difference values b/t
 *           seqs to nearest 0.001). 
 *
 *        *. if(do_pickone): in mode 2 or 3, we select a single
 *           sequence from each cluster to represent that cluster. 
 *           The sequence is chosen that has the minimum average
 *           fractional difference with all other seqs in the cluster.
 *           (NOT YET IMPLEMENTED)
 *
 * Args:    
 * MSA    *mmsa         - the master MSA, we cluster the seqs in this guy
 *                        and build a new MSA from each cluster
 * int     do_all       - TRUE (mode 1): each seq is its own cluster
 * int     target_nc    - number of clusters to define (0 indicates mode 2)
 * float   mindiff      - the minimum fractional difference allowed
 *                        between 2 seqs of different clusters
 *                        (0. indicates mode 3) 
 * int     do_pickone   - TRUE to pick a single seq representative
 *                        from each cluster (see * above) NOT YET IMPLEMENTED
 * int     do_orig      - TRUE to include the master MSA as one of the new MSAs
 * int    *ret_num_msa  - number of MSAs in ret_MSA
 * MSA  ***ret_cmsa     - new MSAs, one for each cluster
 *           
 * Return: ret_cmsa (alloc'ed here) and ret_num_msa
 */
int 
MSADivide(MSA *mmsa, int do_all, int target_nc, float mindiff, int do_pickone,
	  int do_orig, int *ret_num_msa, MSA ***ret_cmsa)
{
  int   status;        /* Easel status code */
  MSA **cmsa;          /* the new MSAs we're creating from clusters of seqs in mmsa */
  int   i;             /* counter over sequences */
  int   m;             /* counter over new MSAs */
  int   n;             /* counter over tree nodes */
  char  buffer[50];    /* for naming the new MSAs */
  ESL_TREE    *T;      /* the tree, created by Single-Linkage Clustering */
  ESL_DMATRIX *D;      /* the distance matrix */
  double *diff;        /* [0..T->N-2], diff[n]= min distance between any leaf in right and
		        * left subtree of node n of tree T */
  double *minld;       /* [0..T->N-2], min dist from node to any taxa in left  subtree */
  double *minrd;       /* [0..T->N-2], min dist from node to any taxa in right subtree */
  int     nc;          /* number of clusters/MSAs  */
  int    *clust;       /* [0..T->N-1], cluster number (0..nc-1) this seq is in */
  int    *c_nseq;      /* [0..nc-1], current number of seqs in the MSA for this cluster */
  int    *csize;       /* [0..nc-1], size of each cluster */

  /* Contract check */
  if((!do_all) && (target_nc == 0 && mindiff == 0.)) 
    ESL_EXCEPTION(eslEINCOMPAT, "target_nc is 0 and mindiff is 0.0, exactly one must be non-zero");
  if((!do_all) && (target_nc != 0 && mindiff != 0.)) 
    ESL_EXCEPTION(eslEINCOMPAT, "target_nc is not 0 and mindiff is not 0.0, exactly one must be non-zero");

  if(do_pickone) 
    ESL_EXCEPTION(eslEINCOMPAT, "do_pickone behavior not yet implemented.");

  /* Mode 1: Each seq becomes own MSA. Easy. */
  if(do_all)
    {
      ESL_ALLOC(clust, sizeof(int) * (mmsa->nseq));
      ESL_ALLOC(csize, sizeof(int) * (mmsa->nseq));
      nc = 0;
      /* each seq is its own cluster */
      for(i = 0; i < mmsa->nseq; i++)
	{
	  clust[i] = nc++;
	  csize[i] = 1;
	}
   }
  else /* Mode 2 or Mode 3 */ 
    {
      /* Create distance matrix and infer tree by single linkage clustering */
      esl_dst_CDiffMx(mmsa->aseq, mmsa->nseq, &D);
      esl_tree_SingleLinkage(D, &T);
      esl_tree_SetTaxaParents(T);
      /*esl_tree_WriteNewick(stdout, T);*/
      esl_tree_Validate(T, NULL);

      /* HEREHEREHERE */
      /* determine the diff values: 
       * (use: n_child > n, unless n's children are taxa)
       * diff[n] is minimum distance between any taxa (leaf) in left subtree of 
       * n to any taxa in right subtree of n. 
       */
      ESL_ALLOC(diff,  (sizeof(double) * (T->N - 1)));  /* one for each node */
      ESL_ALLOC(minld, (sizeof(double) * (T->N - 1))); 
      ESL_ALLOC(minrd, (sizeof(double) * (T->N - 1))); 
      for (n = (T->N-2); n >= 0; n--)
	{
	  minld[n] = T->ld[n] + ((T->left[n]  > 0) ? (minld[T->left[n]])  : 0);
	  minrd[n] = T->rd[n] + ((T->right[n] > 0) ? (minrd[T->right[n]]) : 0);
	  diff[n] = minld[n] + minrd[n];
	  diff[n] *= 1000.; 
	  diff[n] = (float) ((int) diff[n]);
	  diff[n] /= 1000.; 
	  printf("diff[n:%d]: %f\n", n, diff[n]);
	}
      free(minld);
      free(minrd);
      /*for (n = 0; n < (T->N-1); n++)
	printf("diff[n:%d]: %f\n", n, diff[n]);
	for (n = 0; n < (T->N-1); n++)
	printf("left[n:%d]: %d right[n:%d]: %d\n", n, T->left[n], n, T->right[n]);*/
      
      if(target_nc == 0) /* Mode 2 */
	{
	  /* Define clusters that are at least mindiff different
	   * from each other. */
	  select_node(T, diff, mindiff, &clust, &nc);
	}
      else /* Mode 3, mindiff == 0.0 (it's in the contract) */
	{
	  /* Find the minimum fractional difference (mindiff) that 
	   * gives exactly target_nc clusters, also define clusters
	   * based on that mindiff, this is all done with find_mindiff(),
	   * which does a binary search for mindiff, we're guaranteed to 
	   * find exactly target_nc clusters b/c diff values are rounded
	   * to nearest 0.001. */
	  if(target_nc > (T->N)) target_nc = T->N; /* max num clusters is num seqs */
	  mindiff = find_mindiff(T, diff, target_nc, &clust, &nc);
	  printf("nc: %d target_nc: %d\n", nc, target_nc);
	}
      /* Determine the size of each cluster */
      ESL_ALLOC(csize, (sizeof(int) * (nc)));
      esl_vec_ISet(csize, nc, 0);
      for(i = 0; i < mmsa->nseq; i++)
	csize[clust[i]]++;

      /*printf("Distance matrix:\n");
	esl_dmatrix_Dump(stdout, D, NULL, NULL);*/
    }

  /* Create one new MSA for each cluster,
   * if(do_orig): keep the original MSA as cmsa[nc] */
  if(do_orig) 
    {
      ESL_ALLOC(cmsa, (sizeof(MSA *) * (nc+1)));
      cmsa[nc] = mmsa;
    }
  else
    ESL_ALLOC(cmsa, (sizeof(MSA *) * (nc)));

   for(m = 0; m < nc; m++)
    {
      cmsa[m]             = MSAAlloc(csize[m], mmsa->alen);
      if(mmsa->desc != NULL) cmsa[m]->desc = sre_strdup(mmsa->desc, -1);
      cmsa[m]->ss_cons    = sre_strdup(mmsa->ss_cons, -1);
      cmsa[m]->nseq       = csize[m];
      n = sprintf (buffer, ".%d", (m+1));
      cmsa[m]->name       = sre_strdup(mmsa->name, -1);
      sre_strcat(&cmsa[m]->name, -1, buffer, (n+1));
      printf("cmsa[m:%d]->name: %s\n", m, cmsa[m]->name);
      printf("mmsa->name: %s\n", mmsa->name);
    }


  /* add each seq to the appropriate MSA */
  ESL_ALLOC(c_nseq, (sizeof(int) * (nc)));
  esl_vec_ISet(c_nseq, nc, 0);
  for(i = 0; i < mmsa->nseq; i++)
    {
      m = clust[i];
      if(m != -1) 
	{
	  cmsa[m]->aseq[c_nseq[m]]     = sre_strdup(mmsa->aseq[i], -1);
	  cmsa[m]->sqname[c_nseq[m]++] = sre_strdup(mmsa->sqname[i], -1);
	}
    }

  if(do_orig) *ret_num_msa = nc+1;
  else        *ret_num_msa = nc;
  *ret_cmsa = cmsa;

  esl_tree_Destroy(T);
  esl_dmatrix_Destroy(D);
  free(diff);
  free(clust);
  free(c_nseq);
  free(csize);
  return eslOK;

 ERROR: 
  if(diff  != NULL) free(diff);
  if(minld != NULL) free(minld);
  if(minrd != NULL) free(minrd);
  if(clust != NULL) free(clust);
  if(c_nseq!= NULL) free(c_nseq);
  if(csize != NULL) free(csize);
  if(cmsa  != NULL) 
    {
      for(m = 0; m < nc; m++)
	if(cmsa[m] != NULL) MSAFree(cmsa[m]);
      free(cmsa);
    }
  return status;
}

/* Function: select_node()
 * EPN, Fri Mar 23 08:48:37 2007 
 * Adapted from SRE's select_node() in maketestset.c originally written
 * for the PROFMARK HMMER benchmark.
 * 
 * Purpose:  Define clusters of the taxa (seqs) in the tree such
 *           that minimum disparity b/t any 2 seqs in different 
 *           clusters is greater than <mindiff> and the number of
 *           clusters is maximized. Return the index of the node
 *           of the tree under which the largest cluster belongs.
 *           
 *           For high disparities, this cluster may contain all
 *           the sequences, and we'll return the root node (0).
 *
 * Args:    
 * ESL_TREE *T        - the tree
 * double   *diff     - [0..T->N-2]: for each node of the tree, the minimum
 *                      distance (sum of branch lengths) from any taxa (leaf)
 *                      in left subtree to any taxa in right subtree.
 * double    mindiff  - (see description above)
 * int     **ret_clust- [0..T->N-1] cluster number this seq is in, alloc'ed, filled here
 * int      *ret_nc   - number of clusters
 *
 * Returns: node index (as explained in Purpose)
 */
static int
select_node(ESL_TREE *T, double *diff, double mindiff, int **ret_clust, 
	    int *ret_nc)
{
  int status;     /* Easel status code */
  ESL_STACK *ns1; /* stack for traversing tree */
  ESL_STACK *ns2; /* another stack for traversing tree */
  int c;	  /* counter for clusters */
  int best;       /* index of current best node */
  int maxsize;    /* size of cluster for best node */
  int n, np;      /* counters over tree nodes */
  int *clust;     /* [1..T->N-1] cluster number this seq is in */

  /*printf("in selec_node mindiff: %f T->N: %d\n", mindiff, T->N);*/
  /* set tree cladesizes if not already set */
  if(T->cladesize == NULL) 
    esl_tree_SetCladesizes(T);

  ESL_ALLOC(clust, (sizeof(int) * T->N));
  esl_vec_ISet(clust, T->N, 0);

  ns1 = esl_stack_ICreate();
  ns2 = esl_stack_ICreate();

  esl_stack_IPush(ns1, 0);	/* push root on stack to start */
  maxsize  = 0;
  best     = 0;
  c        = 0;
  while (esl_stack_IPop(ns1, &n) != eslEOD)
    {
      if ((n == 0 || diff[T->parent[n]] > mindiff) &&
	  diff[n] <= mindiff)
	{			/* we're at a cluster */
	  if (T->cladesize[n] > maxsize) 
	    {
	      maxsize = T->cladesize[n];
	      best = n;
	    }
	  /* determine all taxa in the clade rooted at n*/
	  esl_stack_IPush(ns2, n);	
	  while (esl_stack_IPop(ns2, &np) != eslEOD)
	  {
	    /*printf("np: %d T->left[np]: %d\n", np, T->left[np]);*/
	    if(T->left[np]  <= 0) clust[(-1*T->left[np])]  = c;
	    else esl_stack_IPush(ns2, T->left[np]);
	    if(T->right[np] <= 0) clust[(-1*T->right[np])]  = c;
	    else esl_stack_IPush(ns2, T->right[np]);
	  }
	  c++;
	}
      else			/* we're not a cluster, keep traversing */
	{
	  /*printf("n: %d T->left[n]: %d\n", n, T->left[n]);*/
	  if(T->left[n]  <= 0) clust[(-1*T->left[n])]  = c++; /* single seq with its own cluster */
	  else esl_stack_IPush(ns1, T->left[n]);
	  if(T->right[n] <= 0) clust[(-1*T->right[n])] = c++; /* single seq with its own cluster */
	  else esl_stack_IPush(ns1, T->right[n]);
	}
    }
  esl_stack_Destroy(ns1);
  esl_stack_Destroy(ns2);
  *ret_nc = c;
  *ret_clust = clust;
  printf("nc: %d(%d) best: %d maxsize: %d nc: %d\n\n", *ret_nc, c, best, maxsize, c);
  for(n = 0; n < T->N; n++)
    {
      printf("clust[%d]: %d\n", n, clust[n]);
    }
  return best;
 ERROR: 
  if(clust != NULL) free(clust);
  return status;
}


/* Function: find_mindiff()
 * EPN, Fri Mar 23 18:59:42 2007
 * 
 * Purpose:  Given a tree resulting from single linkage clustering,
 *           find the min fractional difference (mindiff) that when used to
 *           define clusters (such that no seq in cluster A is less
 *           than mindiff different than any seq in cluster B), 
 *           gives >= target_nc.
 *
 * Args:    
 * ESL_TREE *T        - the tree
 * double   *diff     - [0..T->N-2]: for each node of the tree, the minimum
 *                      distance (sum of branch lengths) from any taxa (leaf)
 *                      in left subtree to any taxa in right subtree.
 * int      target_nc - number of clusters we want
 * int     **ret_clust- [0..T->N-1] cluster number this seq is in, alloc'ed, filled here
 * int      *ret_nc   - number of clusters
 *
 * Returns: fractional difference (as explained in Purpose)
 */
static float
find_mindiff(ESL_TREE *T, double *diff, int target_nc, int **ret_clust, int *ret_nc)
{
  float high = 1.0;
  float low  = 0.0;
  int   high_nc = 0;
  int   low_nc = 0;
  float mindiff = 0.5;
  int curr_nc = -1;
  int keep_going = TRUE;
  float thresh = 0.001;
  int *clust = NULL;

  /* Contract check */
  if(target_nc > T->N) 
    ESL_EXCEPTION(eslEINCOMPAT, "desired number of clusters is greater than number of seqs in the tree");

  while(keep_going)
    {
      if(clust != NULL) free(clust);
      select_node(T, diff, mindiff, &clust, &curr_nc);
      if(curr_nc < target_nc)
	{
	  high       = mindiff;
	  high_nc    = curr_nc;
	  mindiff   -= (mindiff - low) / 2.;
	  if((fabs(high-0.) < thresh) && (fabs(low-0.) < thresh))  keep_going = FALSE; 
	  /* stop, high and low have converged at 0. */
	  printf("LOWER   nc: %d mindiff: %f low: %f high: %f\n", curr_nc, mindiff, low, high);
	}
      else /* curr_nc >= target_nc */
	{
	  low        = mindiff;
	  low_nc     = curr_nc;
	  mindiff   += (high - mindiff) / 2.;
	  if(fabs(high-low) < thresh)  keep_going = FALSE; /* stop, high and low have converged */
	  printf("GREATER nc: %d mindiff: %f low: %f high: %f\n", curr_nc, mindiff, low, high);
	}
    }
  /* it's possible we can't reach our target, if so, set mindiff as minimum value that gives 
   * less than target_nc clusters. */
  if(curr_nc != target_nc)
    {
      printf("targ: %d curr: %d low: %d (%f) high: %d (%f)\n", target_nc, curr_nc, low_nc, low, high_nc, high);
      if(high_nc < target_nc)
	{
	  mindiff = high;
	  select_node(T, diff, mindiff, &clust, &curr_nc);
	}
      else
	while(high_nc > target_nc)
	  {
	    high += thresh;
	    if(high > 1.0)  ESL_EXCEPTION(eslEINCONCEIVABLE, "mindiff has risen above 1.0");
	    mindiff = high;
	    select_node(T, diff, mindiff, &clust, &curr_nc);
	    high_nc = curr_nc;
	  }
    }
  printf("FINAL mindiff: %f\n", mindiff);  
  *ret_nc    = curr_nc;
  *ret_clust = clust;

  return mindiff;
}


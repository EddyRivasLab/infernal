/* bandcyk.c
 * SRE, Wed Nov 20 07:46:56 2002 [flight home from Airlie mtg]
 * SVN $Id$
 * 
 * Banded CYK implementation.
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************  
 */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "squid.h"
#include "vectorops.h"
#include "sre_stack.h"

#include "structs.h"
#include "funcs.h"


void
BandExperiment(CM_t *cm)
{
  int  W;
  int *dmin, *dmax;

  W = 1000;
  while (! BandCalculationEngine(cm, W, 0.00001, FALSE, &dmin, &dmax, NULL, FALSE))
    {
      W += 1000;
      SQD_DPRINTF1(("increasing W to %d, redoing band calculation...\n", W));
    }
}

/* Function:  BandCalculationEngine()
 * Incept:    SRE, Sat Oct 11 14:17:40 2003 [St. Louis]
 *
 * Purpose:   Given a CM and a maximum length W;
 *            calculate probability densities gamma_v(n), probability
 *            of a parse subtree rooted at state v emitting a sequence
 *            of length n.
 *            Then use these to return bounds dmin[v] and dmax[v] which
 *            include a probability mass of >= 1-2(p_thresh).
 *            Each truncated tail (left and right) contains <= p_thresh
 *            probability mass.
 *             
 *            Let L_v(n) be the cumulative probability distribution,
 *            P(length <= n), for state v:
 *                L_v(n) = \sum_{i=0}^{n}  \gamma_v(i)
 *                
 *            For each state v, find dmin such that the probability
 *            of missing a hit is <= p on the low side:
 *                dmin = max_dmin L_v(dmin-1) <= p
 *                
 *            On the high side, let H_v(n) be 1-L_v(n): e.g.
 *            P(length > n) for state v. But it is important not 
 *            to calculate this as 1-L_v(n); because of numerical
 *            roundoff issues, it must be done as:
 *                H_v(n) = \sum_{i=n+1}{\infty} \gamma_v(i)
 *                
 *            Then for each state v, find a dmax such that the 
 *            probability of missing a hit is <= p on the high
 *            side:
 *               dmax = min_dmax  H_v(dmax) <= p
 *               
 *            Note on truncation error:
 *            Of course we can't calculate the sum to \infty; we have
 *            to truncate somewhere. Truncation error must be negligible,
 *            else our choice of bands will depend on the choice of W.
 *            See BandTruncationNegligible() for the test.
 *            
 *
 * Args:      cm        - model to build the bands for
 *            W         - maximum subsequence length W.
 *            p_thresh  - tail probability mass; bounds will be set
 *                        so that we miss <= p_thresh of the mass in
 *                        the left tail and the right tail.
 *            save_densities - TRUE if we want to keep all of the
 *                        gamma matrix. Probably only useful if you're
 *                        also asking for gamma to be returned. Memory usage
 *                        is O(WM). If FALSE, uses a memory-efficient O(W lnM) 
 *                        algorithm instead; if you get the gamma matrix
 *                        back, you're guaranteed the root gamma[0]
 *                        is valid, but don't count on anything else.
 *            ret_dmin  - RETURN: dmin[v] is the minimum subsequence length
 *                        (inclusive) that satisfies p_thresh for left tail.
 *                        Pass NULL if you don't want dmin back.
 *            ret_dmax  - RETURN: dmax[v] is the maximum subsequence length
 *                        (inclusive) that satisfies p_thresh for right tail.
 *                        Pass NULL if you don't want dmax back.
 *            ret_gamma - RETURN: gamma[v][n], [0..M-1][0..W], is the prob
 *                        density Prob(length=n | parse subtree rooted at v).
 *            do_local  - TRUE to factor in possibility of jumping from root
 *                        to any consensus state (see EPN for changes)
 *
 * Returns:   1 on success.
 *            0 if W was too small; caller needs to increase W and
 *            call the engine again.
 *            
 *            The dependency on a sensible W a priori is an annoyance;
 *            it may be possible to write an algorithm that calculates
 *            W on the fly, but I don't see it.
 *
 * Xref:      STL7 p.127 - rearranged calculations of previous BandBounds()
 *                         imp., for better numerical precision
 *            STL7 p.128 - justification of truncation error calculations. 
 *            STL7 p.130 - tests/evaluations.
 */
int
BandCalculationEngine(CM_t *cm, int W, double p_thresh, int save_densities,
		      int **ret_dmin, int **ret_dmax, double ***ret_gamma,
		      int do_local)
{
  double **gamma;               /* P(length = n) for each state v            */
  double  *tmp;
  int     *dmin;                /* lower bound for band. */
  int     *dmax;                /* upper bound for band. */
  int      v;			/* counter over states, 0..M-1               */
  int      y;			/* counter over connected states             */
  int      n;			/* counter over lengths, 0..W */
  int      dv;			/* Delta for state v */
  int      leftn;		/* length of left subsequence under a bifurc */
  double   pdf; 		/* P(<=n) or P(>=n) for this state v         */
  int     *touch;               /* touch[y] = # higher states depending on y */
  Mstack_t *beamstack;          /* pool of beams we can reuse  */
  int      status;		/* return status. */
  int      nd;                  /* counter over nodes */

  /* gamma[v][n] is Prob(state v generates subseq of length n)
   */
  gamma = MallocOrDie(sizeof(double *) * cm->M);        
  for (v = 0; v < cm->M; v++) gamma[v] = NULL;

  /* dmin[v] and dmax[v] are the determined bounds that we return.
   */
  dmin = MallocOrDie(sizeof(int) * cm->M);
  dmax = MallocOrDie(sizeof(int) * cm->M);  
  if (ret_dmin != NULL) *ret_dmin = NULL;
  if (ret_dmax != NULL) *ret_dmax = NULL;

  /* beamstack is a trick for reusing memory: a pushdown stack of 
   * "beams" (gamma[v] rows) we can reuse.
   */
  beamstack = CreateMstack();

  /* The second component of memory saving is the "touch" array.
   * touch[y] is the number of states above state [y] that will
   * depend on y but haven't been calculated yet. When we're done
   * calculating a new state v, we decrement touch[y] for all
   * y \in C_v. Any time touch[y] reaches 0, we put that beam
   * back into the pool for reuse.
   */
  touch  = MallocOrDie(sizeof(int) * cm->M);
  for (v = 0; v < cm->M; v++) touch[v] = cm->pnum[v];

  /* Allocate and initialize the shared end beam.
   */
  gamma[cm->M-1] = MallocOrDie(sizeof(double) * (W+1));
  DSet(gamma[cm->M-1], W+1, 0.);
  gamma[cm->M-1][0] = 1.0;

  for (v = cm->M-1; v >= 0; v--)
    {
      /* Get a beam of memory from somewhere.
       *   1. If we're an E_st, we're sharing the end beam, and
       *      it's already initialized for us; don't do anything
       *      else to it. Bounds are 0..0 by definition.
       *   2. If there's a beam in the pool we can reuse, take it
       *      and set it back to 0's. 
       *   3. Else, allocate and initialize to 0's.
       */
      if (cm->sttype[v] == E_st) {
	gamma[v] = gamma[cm->M-1];
	dmin[v]  = dmax[v] = 0;
	continue;
      }

      if ((gamma[v] = PopMstack(beamstack)) == NULL)
	gamma[v] = MallocOrDie(sizeof(double) * (W+1));		    
      DSet(gamma[v], W+1, 0.);


      /* Recursively calculate prob density P(length=n) for this state v.
       * (The heart of the algorithm is right here.)
       */
      if (cm->sttype[v] == B_st) 
	{			/* a bifurcation state: */
	  pdf = 0.;
	  for (n = 0; n <= W; n++)
	    {
	      for (leftn = 0; leftn <= n; leftn++) 
		gamma[v][n] += gamma[cm->cfirst[v]][leftn]*gamma[cm->cnum[v]][n-leftn];
	      pdf += gamma[v][n];
	    }
	}
      /*EPN 11.11.05 adding following else if () to handle local begins*/
      else if (do_local && v == 0) /*state 0 is the one and only ROOT_S state*/
	{
	  pdf = 0.;
	  for (n = 0; n <= W; n++)
	    {
	      /* Step through by nodes, not states. Only one local begin transition 
	       * is possible into each MATP, MATL, MATR, and BIF nodes, specifically
	       * to the first state of that node.
	       */
	      for (nd = 1; nd < cm->nodes; nd++)
		if (cm->ndtype[nd] == MATP_nd || 
		    cm->ndtype[nd] == MATL_nd ||
		    cm->ndtype[nd] == MATR_nd || 
		    cm->ndtype[nd] == BIF_nd) 
		  {
		    gamma[v][n] += cm->begin[cm->nodemap[nd]] * gamma[cm->nodemap[nd]][n];
		    /* cm->begin[y] is probability we transition to state y from root */
		  }
	      pdf += gamma[v][n];
	    }	      
	}
      /*end EPN block*/
      else 
	{
	  pdf = 0.;
	  dv = StateDelta(cm->sttype[v]);
	  for (n = dv; n <= W; n++)
	    {
	      for (y = 0; y < cm->cnum[v]; y++)
		{
		  gamma[v][n] += cm->t[v][y] * gamma[cm->cfirst[v] + y][n-dv];
		  /* EPN 11.11.05 
		   * Factor in local exit transition probability where appropriate.
		   */
		  if((do_local) && (y == 0)
		     && ((cm->ndtype[cm->ndidx[v]] == MATP_nd) ||
			 (cm->ndtype[cm->ndidx[v]] == MATL_nd) ||
			 (cm->ndtype[cm->ndidx[v]] == MATR_nd) ||
			 (cm->ndtype[cm->ndidx[v]] == BEGL_nd) ||
			 (cm->ndtype[cm->ndidx[v]] == BEGR_nd)))
		    {
		      /* That horrendous if() is survived if
		       * (1) do_local is TRUE
		       * (2) state y is the first state of it's node
		       * (3) state y belongs to a node that allows local
		       *     exits from its first state.
		       * For such states, we want to factor in the local
		       * exit transition probability.
		       */
		      gamma[v][n] += cm->end[v] * gamma[cm->cfirst[v]][n-dv];
		    }
		  /*end EPN block*/
		}
	      pdf += gamma[v][n];
	    }
	}

      /* Make sure we've captured "enough" of the distribution (e.g.,
       * we have captured the right tail; our truncation error is
       * negligible).
       *   Of our 3 criteria, we apply two to every state:
       *     1. we're on the right side of the density (pdf is > 0.5
       *        would be enough, but we use .999)
       *     2. gamma_v(W) < p * DBL_EPSILON
       *        Must be true if \sum_{i=W+1...\infty} g(i) < p*DBL_EPSILON,
       *        which is really what we're trying to prove
       */
      if (pdf <= 0.999 || gamma[v][W] > p_thresh * DBL_EPSILON)
	{
	  /* fail; truncation error is unacceptable; 
	   * caller is supposed to increase W and rerun. 
	   */
	  /*
	    printf("pdf : %f\n", pdf);
	    printf("gamma[v][W] : %f\n", gamma[v][W]);
	  */
	  status = 0; 
	  goto CLEANUP;
	}
      
      /* Renormalize this beam. (Should we really be doing this?)
       */
      if (pdf > 1.0) DNorm(gamma[v], W+1);

      /* Determine our left bound, dmin.
       */
      pdf = 0.;
      for (n = 0; n <= W; n++)
	{
	  pdf += gamma[v][n];
	  if (pdf > p_thresh) { dmin[v] = n; break; }
	}

      /* And our right bound, dmax.
       */
      pdf = 0.;
      for (n = W; n >= 0; n--)
	{
	  pdf += gamma[v][n];
	  if (pdf > p_thresh) { dmax[v] = n; break; }
	}

      /* Reuse memory where possible, using the "touch" trick:
       *   look at all children y \in C_v.
       *   decrement touch[y]
       *   if touch[y] reaches 0, no higher state v depends on this
       *     state's numbers; release the memory.
       *   we're reusing the end state for every E, so don't free it
       *   'til we're done.
       * But it if the save_densities flag is up, don't do this - the
       * caller is telling us to keep the whole gamma matrix around,
       * prob because it's going to be returned and examined.
       * 
       * EPN 11.11.05
       * If we're doing banded local, we don't want to free any match states,
       * to enforce this, we (hackishly) just don't reuse beams. Although
       * we could just save the match state beams...
       */
      if ((! save_densities) && (! do_local)) {
	if (cm->sttype[v] == B_st)
	  {  /* connected children of a B st are handled specially, remember */
	    y = cm->cfirst[v]; PushMstack(beamstack, gamma[y]); gamma[y] = NULL;
	    y = cm->cnum[v];   PushMstack(beamstack, gamma[y]); gamma[y] = NULL;
	  }
	else
	  {
	    for (y = cm->cfirst[v]; y < cm->cfirst[v]+cm->cnum[v]; y++)
	      {
		touch[y]--;
		if (touch[y] == 0 && cm->sttype[y] != E_st) {
		  PushMstack(beamstack, gamma[y]); 
		  gamma[y] = NULL;
		}
	      }
	  }
      }

    } /*end loop up through all states v*/

  /* EPN 11.13.05
   * Step through band for each step and ensure that dmax[v] <= dmax[0]
   * for all v. Because all hits must be rooted at state 0, it doesn't
   * make sense to allow the maximum subseq length for any state to 
   * be greater than the maximum allowable length for state 0. Important
   * because the maximum length of a hit in a banded scan or alignment
   * can be reset to dmax[0] after this step.
   */
  for (v = 1; v <= cm->M-1; v++)
    if (dmax[v] > dmax[0]) dmax[v] = dmax[0];

  if (! BandTruncationNegligible(gamma[0], dmax[0], W, NULL)) 
    { status = 0; goto CLEANUP; }

  status = 1;

 CLEANUP:
  free(touch);

  if (ret_dmin  != NULL) *ret_dmin = dmin;   else free(dmin);
  if (ret_dmax  != NULL) *ret_dmax = dmax;   else free(dmax);
  if (ret_gamma != NULL) *ret_gamma = gamma; else FreeBandDensities(cm, gamma);
  
  /* Free the reused stack of beams.
   */
  while ((tmp = PopMstack(beamstack)) != NULL) free(tmp);
  FreeMstack(beamstack);

  return status;
}



/* Function:  BandTruncationNegligible()
 * Incept:    SRE, Sun Oct 19 10:43:21 2003 [St. Louis]
 *
 * Purpose:   Verifies that for this choice of W, the truncation error,
 *                 D = \sum_{n=W+1}{\infty} \gamma_v(n), 
 *            will not affect our calculation of the right bound, b (dmax),
 *            based on the probability mass we can observe,
 *                 C = \sum_{n=b+1}{W} \gamma_v(n).     
 *                 
 *            Specifically, we want D such that C + D = C; that is,
 *                 D < C * DBL_EPSILON.
 *                 
 *            We assume that the tail of \gamma is decreasing
 *            geometrically. This lets us predict the tail is
 *                  \gamma_v(n) = \gamma_v(W+1) \beta^{n-W-1} 
 *                  
 *            and from the sum of an infinite geometric series, combined
 *            with \gamma_v(W+1) = \beta \gamma_v(W), we obtain:
 *                            \gamma_v(W) \beta
 *                    D'  =     -----------------
 *                                1 - \beta
 *                                    
 *            How well D' approximates the true truncated tail mass D
 *            depends on how valid the assumption of geometrically
 *            decreasing tails is. For a single insert state, one
 *            obtains a geometrically decreasing tail.  For a mixture
 *            of geometric distributions, if one fits the low side,
 *            one overestimates the rate of convergence to 0, and so
 *            underestimates D; D' <= D. But if anything, we want D' >= D;
 *            an upper bound on D lets us prove D < C * epsilon.
 *            Puzzlingly, this does not seem to be a problem. For
 *            a variety of models and states, D' is indeed an
 *            overestimate of D; empirically, the tail density converges to zero
 *            supergeometrically, which I can't explain.
 *            
 *            Using ret_beta to verify:
 *            
 *               D' = (beta / (1.-beta)) * density[W];
 *               D  = \sum_{n=b+1}{\infty} density[W];
 *               D' >= D.
 *               
 *               test for an even stronger criterion:
 *               let estimated density g(n) = gamma[W] * \beta^(n-W);
 *               g(n) >= gamma[n] for all n > W.
 *            
 *            In the testsuite, "check_bandtruncation" empirically verifies 
 *            D' > D. 
 *
 * Args:      density    - one density \gamma_v() calculated by BandCalculationEngine();
 *                         usually the root (if W is big enough for the root, it's
 *                         big enough for every state).
 *            b          - the left bound dmax[v]
 *            W          - the maximum length; gamma_v[] runs [0..W]
 *            ret_beta   - RETURN (optional): the geometric decay constant \beta,
 *                         obtained by simple linear fit to log gamma().
 *                         
 *
 * Returns:   1 if truncation error is negligible (D' < C * DBL_EPSILON)
 *            0 if truncation error is not negigible, and caller will
 *              have to worry about increasing W.
 *
 * Xref: STL7 p.128.  
 */
int
BandTruncationNegligible(double *density, int b, int W, double *ret_beta)
{
  double logbeta;
  double beta;		/* geometric decay parameter                  */
  double C;		/* area under density from b+1..W inclusive   */
  double D;		/* area under unseen density from W+1..\infty */
  int    i;
  
  /* Sum up how much probability mass we do see,
   * in the truncated tail from b+1..W.
   */
  C = 0.;
  for (i = b+1; i <= W; i++) C += density[i];

  /* If density is falling off as a geometric, log(beta) is 
   * the slope of log(p). Estimate slope quickly and crudely, by a
   * simple 2-point fit at our boundaries b+1 and W.  
   */
  logbeta = (log(density[W]) - log(density[b+1])) / (W - b - 1);
  beta = exp(logbeta);
	     
  /* We can now guess at the missing probability mass from W+1...\infty,
   * because a finite geometric series converges to 1/(1+\beta).
   */
  D = (beta / (1.-beta)) * density[W];

  if (ret_beta != NULL) *ret_beta = beta;

  if (D < C * DBL_EPSILON) return 1;
  else                     return 0;
}  
  
/* Function:  BandMonteCarlo()
 * Incept:    SRE, Fri Oct 17 08:01:42 2003 [St. Louis]
 *
 * Purpose:   Calculate the gamma_v densities by Monte Carlo simulation,
 *            by sampling parsetrees from the CM. gamma[v][i] is the
 *            observed number of subsequences of length i rooted at v
 *            in the sample. These counts are left unnormalized; different
 *            states are reached with different probabilities, and the 
 *            caller may want to test the Monte Carlo observed counts
 *            against predicted counts (see bandcyk-montecarlo-test in
 *            the testsuite).
 *
 * Args:      cm         - the model to sample from
 *            nsample    - number of Monte Carlo sampled parsetrees    
 *            W          - maximum subsequence length in densities
 *            ret_gamma  - RETURN: gamma[v][n], [0..M-1][0..W] as
 *                         unnormalized observed counts.
 *
 * Returns:   1 on success.
 *            0 on failure: one or more samples had a length too great to be
 *              captured by W.
 *
 *            Caller frees the returned gamma with FreeBandDensities().
 */
int
BandMonteCarlo(CM_t *cm, int nsample, int W, double ***ret_gamma)
{
  Parsetree_t  *tr;             /* sampled parsetree */
  double      **gamma;          /* RETURN: the densities */
  int           i;		/* counter over samples */
  int           seqlen;		/* length of sampled sequence */
  int           tidx;		/* index on parsetree */
  int           v;		/* state used at a parsetree node */
  int           n;		/* subseq length at a parsetree node */
  int           status;		/* return status. */
  
  /* Allocate gamma, completely; and initialize to zeros. 
   * For consistency w/ BandCalculationEngine(), allocate a single
   * shared end deck at M-1, and point other ends at it - even
   * though P(n=0) = 1.0 by definition at the E's and we don't
   * really need to calculate it. Then we can use FreeBandDensities()
   * for gamma matrices alloc'ed in either function.
   */                                                     
  gamma          = MallocOrDie(sizeof(double *) * cm->M);
  gamma[cm->M-1] = MallocOrDie(sizeof(double) * (W+1)); 
  DSet(gamma[cm->M-1], W+1, 0.0);
  for (v = 0; v < cm->M-1; v++)
    {
      if (cm->sttype[v] != E_st)
	{
	  gamma[v] = MallocOrDie(sizeof(double) * (W+1));
	  DSet(gamma[v], W+1, 0.0);
	}
      else
	gamma[v] = gamma[cm->M-1];
    }

  /* Count Monte Carlo samples of subsequence lengths for
   * all nodes of sampled parsetrees.
   */
  status = 1;			
  for (i = 0; i < nsample; i++)
    {
      EmitParsetree(cm, &tr, NULL, NULL, &seqlen);
      if (seqlen > W) {
	FreeParsetree(tr);
	status = 0;		/* set status to FAILED */
	continue;
      }

      /* The parsetree, though it's a tree, is stored as an
       * array - so traversing it in preorder is trivial. It's
       * already arranged in preorder.
       */        
      for (tidx = 0; tidx < tr->n; tidx++) 
	{
	  v = tr->state[tidx];
	  n = (cm->sttype[v] == E_st) ? 0 : tr->emitr[tidx] - tr->emitl[tidx] + 1;
	  gamma[v][n] += 1.;
	}
      FreeParsetree(tr);
    }

  /* Return gamma, the observed counts (unnormalized densities).
   */
  *ret_gamma = gamma;
  return status;
}


/* Function:  FreeBandDensities()
 * Incept:    SRE, Thu Oct 16 08:30:47 2003 [St. Louis]
 *
 * Purpose:   Free a gamma[] array that was returned by BandCalculationEngine().
 *            Best to handle this with a special function because of the reuse
 *            of the END rows - only cm->M-1 is actually allocated, and other
 *            ENDs just point at that one. Too easy to double free() if we
 *            leave this tricky business to the caller.                 
 *            
 * Args:      cm    - the model we build the band densities, gamma[], for.
 *            gamma - the band densities. Doesn't matter if this is a full
 *                    matrix (save_densities = TRUE) or a partial matrix
 *                    (save_densities = FALSE).
 *
 * Returns:   (void)
 *
 * Xref:      STL7 p130.
 */
void
FreeBandDensities(CM_t *cm, double **gamma)
{
  int v;
  for (v = 0; v < cm->M; v++) 
    if (cm->sttype[v] != E_st && gamma[v] != NULL) 
      { free(gamma[v]); gamma[v] = NULL; }
  free(gamma[cm->M-1]);		/* free the end state */
  free(gamma);
}  


/* Function: BandDistribution()
 * Date:     SRE, Wed Nov 20 07:15:31 2002 [flight home from Washington DC]
 *
 * Purpose:  Given a model, and a DP maximum window size W; calculate
 *           cumulative distribution g_v(n) = P(len <= n | subgraph rooted at v).
 *           Return as matrix g[v][n].
 *
 * Args:     cm   - CM to calculate length distribution for.
 *           W    - maximum DP window size.       
 *           do_local - TRUE to factor in possibility of jumping from root
 *                      to any consensus state (see EPN for changes)
 *
 * Returns:  gamma[v][n] (0..M-1, 0..W).
 *           Caller free's w/ DMX2Free(gamma).
 *
 * Deprecated 11.11.05 EPN
 * Use BandCalculationEngine() instead. 
 * This function does not properly handle potential truncation error
 * problems. No longer supported.
 * 
 */
double **
BandDistribution(CM_t *cm, int W, int do_local)
{
  double **gamma;            /* gamma[v][n] = log P(length n | state v); [0..W][0..M-1] */
  int      n,x;
  int      v,y;

  n     = MAX(MAXCONNECT, W+1);
  gamma = DMX2Alloc(cm->M, W+1);
  
  printf("BandDistribution() is deprecated.\nUse BandCalculationEngine() instead (its better).\n");
  exit(1);

  for (n = 0; n <= W; n++)
    for (v = cm->M-1; v >= 0; v--)
      {
	gamma[v][n] = 0.;

	switch (cm->sttype[v]) {
	case S_st:
	  /*EPN Handle local begins.*/
	  if(do_local && v == 0) 
	    {
	      for (y = 0; y < cm->M; y++) {
		if(cm->sttype[y] == MP_st ||
		   cm->sttype[y] == ML_st ||
		   cm->sttype[y] == MR_st ||
		   cm->sttype[y] == B_st) {
		  gamma[v][n] += cm->begin[y] * gamma[y][n];
		  /* cm->begin[y] is probability we transition to state y from root */
		}
	      }
	    }
	  else
	    {
	      for (y = 0; y < cm->cnum[v]; y++)
		{
		  gamma[v][n] += cm->t[v][y] * gamma[cm->cfirst[v] + y][n];
		  /* EPN 11.11.05 - factor in local exit transition probability*/
		  if((do_local) && (y == 0))
		    {
		      /* That if() is survived if
		       * (1) do_local is TRUE
		       * (2) state y is the first state of it's node
		       *     (this means its node is a MATL or MATR node)
		       * (3) state v isn't the ROOT_S state (would have entered
		       *     earlier if(do_local && v==0) if that were true)
		       *     This means v is either a BEGL_S or BEGR_S
		       * For such states, we want to factor in the local
		       * exit transition probability.
		       */
		      gamma[v][n] += cm->end[v] * gamma[cm->cfirst[v] + y][n];
		    }
		}
	      /*end EPN block*/	  
	    }
	  break;
	  
	case D_st:
	  for (y = 0; y < cm->cnum[v]; y++)
	    gamma[v][n] += cm->t[v][y] * gamma[cm->cfirst[v] + y][n];
	  break;
	  
	case ML_st:
	case MR_st:
	  if (n >= 1) 
	    {
	      for (y = 0; y < cm->cnum[v]; y++)
		{
		  gamma[v][n] += cm->t[v][y] * gamma[cm->cfirst[v] + y][n-1];
		  /* EPN 11.11.05 - factor in local exit transition probability*/
		  if((do_local) && (y == 0))
		    {
		      /* That if() is survived if
		       * (1) do_local is TRUE
		       * (2) state y is the first state of it's node
		       *     (this means its node is a MATL or MATR node)
		       * For such states, we want to factor in the local
		       * exit transition probability.
		       */
		      gamma[v][n] += cm->end[v] * gamma[cm->cfirst[v] + y][n-1];
		    }
		}
	    }
	    /*end EPN block*/	  
	  break;

	case IL_st:
	case IR_st:
	  if (n >= 1) {
	    for (y = 0; y < cm->cnum[v]; y++)
	      gamma[v][n] += cm->t[v][y] * gamma[cm->cfirst[v] + y][n-1];
	  }
	  break;

	case MP_st:
	  if (n >= 2) 
	    {
	      for (y = 0; y < cm->cnum[v]; y++)
		gamma[v][n] += cm->t[v][y] * gamma[cm->cfirst[v] + y][n-2];
	      /* EPN 11.11.05 - factor in local exit transition probability*/
	      if(do_local)
		{
		  /* If we get here, state v is a MATP_MP state,
		   * so y = 0. 
		   * For such states, we want to factor in the local
		   * exit transition probability.
		   */
		  gamma[v][n] += cm->end[v] * gamma[cm->cfirst[v] + y][n-2];
		}
	    }
	  break;
  
	case B_st:
	  for (x = 0; x <= n; x++) /* x = length of left child */
	    gamma[v][n] += gamma[cm->cfirst[v]][x] * gamma[cm->cnum[v]][n-x];
	  break;

	case E_st:
	  if (n == 0) gamma[v][n] = 1.;
	  break;
	
	default: 
	  Die("gamma on fire");
	}
      }

  /* Reduce numerical imprecision issues: renormalize the distributions.
   * Don't do this if you're debugging; it makes all other problems go away too.
   */
  for (v = 0; v < cm->M; v++)
    DNorm(gamma[v], W+1);

  /* Convert to cumulative distribution, P(len <= n) 
   * Squash any numerical imprecision that overflows past 1.0.
   */
  for (v = 0; v < cm->M; v++)
    for (n = 1; n <= W; n++)
      {
	gamma[v][n] = gamma[v][n] + gamma[v][n-1];
	if (gamma[v][n] > 1.) gamma[v][n] = 1.;
      }
				
  return gamma;
}

/* Function: BandBounds()
 * Date:     SRE, Fri Sep 26 10:01:41 2003 [AA 886 from Salk]
 *
 * Purpose: Given a cumulative probability distribution gamma[n] =
 *          P(length <= n) for each state v and a tail
 *          cutoff probability p.
 * 
 *          Find dmin and dmax for each v, such that the probability
 *          of missing a hit is <= p on either the low side or high side:
 *             gamma[dmin-1]   <= p
 *             1 - gamma[dmax] <= p
 * 
 *          The total probability mass these bounds encompass is:
 *             gamma[dmax] - gamma[dmin-1] >= 1-2p   
 *
 * Args:    gamma    - cumulative probability distribution P(length <= n) for state v;
 *                     [0..v..M-1][0..W] 
 *          M        - # of states in CM
 *          W        - maximum subsequence length in DP         
 *          p        - tail cutoff probability
 *          ret_dmin - RETURN: minimum d bound for each state v; [0..v..M-1]
 *          ret_dmax - RETURN: maximum d bound for each state v; [0..v..M-1]
 *
 * Returns: (void) 
 *          dmin, dmax are allocated here. Caller must free.
 */
void
BandBounds(double **gamma, int M, int W, double p, int **ret_dmin, int **ret_dmax)
{
  int     *dmin, *dmax;
  int      v;

#if (SRE_CONLEVEL >= 1)
  for (v = 0; v < M; v++) 
    {
      assert(gamma[v][W] <= 1.0 && gamma[v][0] >= 0.0);
      assert(gamma[v][W] >= gamma[v][0]);
      assert(p >= 0. && p <= 0.5);
    }
#endif

  dmin = MallocOrDie(sizeof(int) * M);
  dmax = MallocOrDie(sizeof(int) * M);

  for (v = 0; v < M; v++)
    {
      dmin[v] = 0; while (dmin[v] <= W && gamma[v][dmin[v]]   <= p)    dmin[v]++;
      dmax[v] = W; while (dmax[v] > 0  && gamma[v][dmax[v]-1] >= 1.-p) dmax[v]--;
    }

#if (SRE_CONLEVEL >= 1)
  for (v = 0; v < M; v++)
    {
      if (dmin[v] > 0) 
	assert(gamma[v][dmax[v]] - gamma[v][dmin[v]-1] >= 1.-2*p);
      else
	assert(gamma[v][dmax[v]] >= 1.-2*p);
      assert(dmin[v] >= 0 && dmin[v] <= W);
      assert(dmax[v] >= 0 && dmax[v] <= W);
      assert(dmax[v] >= dmin[v]);
    }
#endif  
  
  *ret_dmin = dmin;
  *ret_dmax = dmax;
}
      

/* A couple of quick hacks. ...
 * Print an XMGRACE xy file for a specified v, showing the
 * cumulative distribution. Needed this for the R01 renewal.
 * SRE, Wed Feb 19 08:35:32 2003
 */
void
PrintBandGraph(FILE *fp, double **gamma, int *min, int *max, int v, int W)
{
  int n;

  for (n = 0; n <= W; n++)
    fprintf(fp, "%d %.6f\n", n, gamma[v][n]);
  fprintf(fp, "&\n");
  fprintf(fp, "%d  0\n",   min[v]);
  fprintf(fp, "%d  1.0\n", min[v]);
  fprintf(fp, "&\n");
  fprintf(fp, "%d  0\n",   max[v]);
  fprintf(fp, "%d  1.0\n", max[v]);
  fprintf(fp, "&\n");
}
/* ... and, estimate the total savings in DP cells filled.
 */
void
PrintDPCellsSaved(CM_t *cm, int *min, int *max, int W)
{
  int v;
  int after, before;

  before = after = 0;
  for (v = 0; v < cm->M; v++) 
    {
      if (cm->sttype[v] != E_st) {
	after  += max[v] - min[v] + 1;
	before += W;
      }
    }
  printf("Before:  something like %d\n", before);
  printf("After:   something like %d\n", after);
  printf("Speedup: maybe %.2f fold\n", (float) before / (float) after);
}


/* Function: CYKBandedScan()
 * Date:     SRE, Fri Sep 26 09:53:42 2003 [AA 886, returning from Salk Institute]
 *
 * Purpose:  Scan a sequence for matches to a covariance model, using the
 *           banded algorithm.
 *           Allows multiple nonoverlapping hits and local alignment.
 *           Derived from scancyk.c.
 *
 *           dmin,dmax set the bounds. Both are arrays, 0..v..cm->M. 
 *           The band for v is dmin[v]..dmax[v], inclusive; that is,
 *           dmin[v] is the minimum allowed d for state v (inclusive);
 *           dmax[v] is the maximum; that is, v can only account for
 *           subsequences of length >= dmin[v] and <= dmax[v].
 *
 * Args:     cm        - the covariance model
 *           dsq       - digitized sequence to search; 1..L
 *           dmin      - minimum bound on d for state v; 0..M
 *           dmax      - maximum bound on d for state v; 0..M          
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
CYKBandedScan(CM_t *cm, char *dsq, int *dmin, int *dmax, int i0, int j0, int W, 
	      int *ret_nhits, int **ret_hitr, int **ret_hiti, int **ret_hitj, float **ret_hitsc, float min_thresh)
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
  int       jp;			/* rolling index into BEGL_S decks: jp=j%(W+1) */
  int       nhits;		/* # of hits in optimal parse */
  int      *hitr;		/* initial state indices of hits in optimal parse */
  int      *hiti;               /* start positions of hits in optimal parse */
  int      *hitj;               /* end positions of hits in optimal parse */
  float    *hitsc;              /* scores of hits in optimal parse */
  int       alloc_nhits;	/* used to grow the hit arrays */
  int       jmax;               /* when imposing bands, maximum j value in alpha matrix */
  int       kmax;               /* for B_st's, maximum k value consistent with bands*/
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       gamma_j;            /* j index in the gamma matrix, which is indexed 0..j0-i0+1, 
				 * while j runs from i0..j0 */
  int       gamma_i;            /* i index in the gamma* data structures */

  /* EPN 08.11.05 Next line prevents wasteful computations when imposing
   * bands before the main recursion.  There is no need to worry about
   * alpha cells corresponding to subsequence distances within the windowlen
   * (W) but LONGER than the full sequence (L).  Saves a significant amount 
   * of time if W is much larger than necessary, and the search sequences 
   * are short (as in a possible benchmark).
   */
  L = j0-i0+1;
  if (W > L) W = L; 

  /*PrintDPCellsSaved(cm, dmin, dmax, W);*/

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
   *
   * and, for banding: any cell outside our bands is impossible.
   * These inits are never changed in the recursion, so even with the
   * rolling, matrix face reuse strategy, this works.
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
	  /* treat EL as emitting only on self transition */
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
	for (j = 2; j < W; j++) 
	  alpha[v][j][0] = alpha[v][0][0];
    }
  /* Impose the bands.
   *   (note: E states have all their probability on d=0, so dmin[E] = dmax[E] = 0;
   *    the first loop will be skipped, the second initializes the E states.)
   */
  for (v = 0; v < cm->M; v++)
    {
      if(cm->stid[v] == BEGL_S) jmax = W; 
      else jmax = 1;

      for (d = 0; d < dmin[v] && d <=W; d++) 
	for(j = 0; j <= jmax; j++)
	  alpha[v][j][d] = IMPOSSIBLE;
      
      for (d = dmax[v]+1; d <= W;      d++) 
	for(j = 0; j <= jmax; j++)
	  alpha[v][j][d] = IMPOSSIBLE;
    }

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
	      for (d = dmin[v]; d <= dmax[v] && d <= gamma_j; d++) 
		{
		  y = cm->cfirst[v];
		  alpha[v][jp][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][cur][d] + cm->tsc[v][yoffset]) > alpha[v][jp][d]) 
		      alpha[v][jp][d] = sc;
		  if (alpha[v][jp][d] < IMPROBABLE) alpha[v][jp][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == MP_st) 
	    {
	      for (d = dmin[v]; d <= dmax[v] && d <= gamma_j; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][prv][d-2] + cm->tsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  
		  i = j-d+1;
		  if (dsq[i] < Alphabet_size && dsq[j] < Alphabet_size)
		    alpha[v][cur][d] += cm->esc[v][(int) (dsq[i]*Alphabet_size+dsq[j])];
		  else
		    alpha[v][cur][d] += DegeneratePairScore(cm->esc[v], dsq[i], dsq[j]);
		  
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) 
	    {
	      for (d = dmin[v]; d <= dmax[v] && d <= gamma_j; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][cur][d-1] + cm->tsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  
		  i = j-d+1;
		  if (dsq[i] < Alphabet_size)
		    alpha[v][cur][d] += cm->esc[v][(int) dsq[i]];
		  else
		    alpha[v][cur][d] += DegenerateSingletScore(cm->esc[v], dsq[i]);
		  
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) 
	    {
	      for (d = dmin[v]; d <= dmax[v] && d <= gamma_j; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][prv][d-1] + cm->tsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  
		  if (dsq[j] < Alphabet_size)
		    alpha[v][cur][d] += cm->esc[v][(int) dsq[j]];
		  else
		    alpha[v][cur][d] += DegenerateSingletScore(cm->esc[v], dsq[j]);
		  
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == B_st) 
	    {
	      w = cm->cfirst[v];
	      y = cm->cnum[v];
	      i = j-d+1;
	      for (d = dmin[v]; d <= dmax[v] && d <= gamma_j; d++) 
		{
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));

		  /*EPN : Make sure k is consistent with bands in state w and state y.
		    Not sure if this is necessary because the speed-up will be 
		    very small (if its even a speed-up due to extra computations), 
		    but it does make the banded approach more consistent.   
		    original line : 
		    for (k = 0; k <= d; k++)
		  */

		  /* k is the length of the right fragment */
		  if(dmin[y] > (d-dmax[w])) k = dmin[y];
		  else k = d-dmax[w];
		  if(k < 0) k = 0;

		  if(dmax[y] < (d-dmin[w])) kmax = dmax[y];
		  else kmax = d-dmin[w];
		  
		  for (; k <= kmax; k++)
		    {
		      jp = (j-k)%(W+1);	   /* jp is rolling index into BEGL_S deck j dimension */
		      if ((sc = alpha[w][jp][d-k] + alpha[y][cur][k]) > alpha[v][cur][d])
			alpha[v][cur][d] = sc;
		    }
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
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
      for (d = dmin[0]; d <= dmax[0] && d <= gamma_j; d++)
	{
	  y = cm->cfirst[0];
	  alpha[0][cur][d] = alpha[y][cur][d] + cm->tsc[0][0];
	  bestr[d]         = 0;	/* root of the traceback = root state 0 */
	  for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++)
	    if ((sc = alpha[y+yoffset][cur][d] + cm->tsc[0][yoffset]) > alpha[0][cur][d]) 
	      alpha[0][cur][d] = sc;
	}
      
      /* EPN 11.09.05 
       * The following loop that deals with local begins was modified
       * to enforce bands on all states y that are possible internal entry
       * points. Old code block went from [0..d] in the d dimension
       * for state y.
       * ref: ~nawrocki/notebook/5_1109_inf_local_banded_spd/00LOG
       */

      if (cm->flags & CM_LOCAL_BEGIN) {
	for (y = 1; y < cm->M; y++) {
	  d = (dmin[y] > dmin[0]) ? dmin[y]:dmin[0];
	  /*if (dmin[y] > dmin[0]) d = dmin[y];
	    else d = dmin[0];*/
	  for (; d <= dmax[y] && d <= gamma_j; d++)
	    {
	      if (cm->stid[y] == BEGL_S) sc = alpha[y][j%(W+1)][d] + cm->beginsc[y];
	      else                       sc = alpha[y][cur][d]     + cm->beginsc[y];
	      if (sc > alpha[0][cur][d]) {
		alpha[0][cur][d] = sc;
		bestr[d]         = y;
	      }
	    }
	}
	if (alpha[0][cur][d] < IMPROBABLE) alpha[0][cur][d] = IMPOSSIBLE;
      }
      
      /* The little semi-Markov model that deals with multihit parsing:
       */
      gamma[gamma_j]  = gamma[gamma_j-1] + 0; /* extend without adding a new hit */
      gback[gamma_j]  = -1;
      savesc[gamma_j] = IMPOSSIBLE;
      saver[gamma_j]  = -1;
      for (d = dmin[0]; d <= dmax[0] && d <= gamma_j; d++) 
	{
	  i = j-d+1;
	  gamma_i = j-d+1-i0+1;
	  sc = gamma[gamma_i-1] + alpha[0][cur][d]  - min_thresh; 
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
  alloc_nhits = 10;
  hitr  = MallocOrDie(sizeof(int)   * alloc_nhits);
  hitj  = MallocOrDie(sizeof(int)   * alloc_nhits);
  hiti  = MallocOrDie(sizeof(int)   * alloc_nhits);
  hitsc = MallocOrDie(sizeof(float) * alloc_nhits);
  
  j     = j0;
  nhits = 0;
  while (j >= i0) {
    gamma_j = j-i0+1;
    if (gback[gamma_j] == -1) /* no hit */
      j--; 
    else                /* a hit, a palpable hit */
      {
	hitr[nhits]   = saver[gamma_j];
	hitj[nhits]   = j;
	hiti[nhits]   = gback[gamma_j];
	hitsc[nhits]  = savesc[gamma_j];
	nhits++;
	j = gback[gamma_j]-1;
	
	if (nhits == alloc_nhits) {
	  hitr  = ReallocOrDie(hitr,  sizeof(int)   * (alloc_nhits + 10));
	  hitj  = ReallocOrDie(hitj,  sizeof(int)   * (alloc_nhits + 10));
	  hiti  = ReallocOrDie(hiti,  sizeof(int)   * (alloc_nhits + 10));
	  hitsc = ReallocOrDie(hitsc, sizeof(float) * (alloc_nhits + 10));
	  alloc_nhits += 10;
	}
      }
  }
  free(gback);
  free(gamma);
  free(savesc);
  free(saver);

  *ret_nhits = nhits;
  *ret_hitr  = hitr;
  *ret_hiti  = hiti;
  *ret_hitj  = hitj;
  *ret_hitsc = hitsc;
  return;
}


/* Function: BandedParsetreeDump()
 * Date:     SRE, Sat Sep 27 15:44:22 2003 [St. Louis]
 * Xref:     STL7 pg. 119
 *
 * Purpose:  Generate a detailed picture of a parsetree data structure,
 *           annotated with relevant information from the sequence
 *           and the model -- and show agreement (or disagreement)
 *           with calculated bands.
 *           
 *           Modified from ParsetreeDump().
 *
 * Args:    fp    - FILE to write output to.
 *          tr    - parsetree to examine.
 *          cm    - model that was aligned to dsq to generate the parsetree
 *          dsq   - sequence that was aligned to cm to generate the parsetree
 *          gamma - cumulative subsequence length probability distributions
 *                  used to generate the bands; from BandDistribution(); [0..v..M-1][0..W]
 *          W     - maximum window length W (gamma distributions range up to this)        
 *          dmin  - minimum subseq length for each state; from BandBounds(); [0..v..M-1]
 *          dmax  - maximum ""
 *
 * Returns:  (void)
 */
void
BandedParsetreeDump(FILE *fp, Parsetree_t *tr, CM_t *cm, char *dsq, 
		    double **gamma, int W, int *dmin, int *dmax)
{
  int   x;
  char  syml, symr;
  float tsc;
  float esc;
  int   v,y;
  int   L;

  fprintf(fp, "%5s %6s %6s %7s %5s %5s %5s %5s %5s %5s %5s %5s\n",
	  " idx ", "emitl", "emitr", "state", " nxtl", " nxtr", " prv ", " tsc ", " esc ", 
	  " L   ", " dmin", " dmax");
  fprintf(fp, "%5s %6s %6s %7s %5s %5s %5s %5s %5s %5s %5s %5s\n",
	  "-----", "------", "------", "-------", "-----","-----", "-----","-----", "-----",
	  "-----", "-----", "-----");
  for (x = 0; x < tr->n; x++)
    {
      v = tr->state[x];

      /* Set syml, symr: one char representation of what we emit, or ' '.
       * Set esc:        emission score, or 0.
       * Only P, L, R states have emissions.
       */
      syml = symr = ' ';
      esc = 0.;
      if (cm->sttype[v] == MP_st) {
	syml = Alphabet[(int)dsq[tr->emitl[x]]]; 
	symr = Alphabet[(int)dsq[tr->emitr[x]]];
	esc  = DegeneratePairScore(cm->esc[v], dsq[tr->emitl[x]], dsq[tr->emitr[x]]);
      } else if (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st) {
	syml = Alphabet[(int)dsq[tr->emitl[x]]];
	esc  = DegenerateSingletScore(cm->esc[v], dsq[tr->emitl[x]]);
      } else if (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st) {
	symr = Alphabet[(int)dsq[tr->emitr[x]]];
	esc  = DegenerateSingletScore(cm->esc[v], dsq[tr->emitr[x]]);
      }

      /* Set tsc: transition score, or 0.
       * B, E, and the special EL state (M, local end) have no transitions.
       */
      tsc = 0.;
      if (v != cm->M && cm->sttype[v] != B_st && cm->sttype[v] != E_st) {
	y = tr->state[tr->nxtl[x]];

	if (v == 0 && (cm->flags & CM_LOCAL_BEGIN))
	  tsc = cm->beginsc[y];
	else if (y == cm->M) /* CM_LOCAL_END is presumably set, else this wouldn't happen */
	  tsc = cm->endsc[v] + (cm->el_selfsc * (tr->emitr[x] - tr->emitl[x] + 1 - StateDelta(cm->sttype[v])));
	else 		/* y - cm->first[v] gives us the offset in the transition vector */
	  tsc = cm->tsc[v][y - cm->cfirst[v]];
      }

      /* Print the info line for this state
       */
      L = tr->emitr[x]-tr->emitl[x]+1;
      fprintf(fp, "%5d %5d%c %5d%c %5d%-2s %5d %5d %5d %5.2f %5.2f %5d %5d %5d %2s\n",
	      x, tr->emitl[x], syml, tr->emitr[x], symr, tr->state[x], 
	      Statetype(cm->sttype[v]), tr->nxtl[x], tr->nxtr[x], tr->prv[x], tsc, esc,
	      L, dmin[v], dmax[v],
	      (L >= dmin[v] && L <= dmax[v]) ? "" : "!!");
    }
  fprintf(fp, "%5s %6s %6s %7s %5s %5s %5s %5s %5s %5s %5s %5s\n",
	  "-----", "------", "------", "-------", "-----","-----", "-----","-----", "-----",
	  "-----", "-----", "-----");
  fflush(fp);
} 

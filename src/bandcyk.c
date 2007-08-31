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

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_random.h"
#include "esl_vectorops.h"
#include "esl_stack.h"
#include "esl_sqio.h"

#include "funcs.h"
#include "structs.h"

void
BandExperiment(CM_t *cm)
{
  int  W;
  int *dmin, *dmax;

  W = 1000;
  while (! BandCalculationEngine(cm, W, 0.00001, FALSE, &dmin, &dmax, NULL))
    {
      W += 1000;
      ESL_DPRINTF1(("increasing W to %d, redoing band calculation...\n", W));
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
		      int **ret_dmin, int **ret_dmax, double ***ret_gamma)
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
  ESL_STACK *beamstack;         /* pool of beams we can reuse  */
  int      status;		/* return status. */
  int      nd;                  /* counter over nodes */
  int      yoffset;             /* counter over children */
  int      reset_local_ends;    /* TRUE if we erased them and need to reset */
  int      reset_cp9_local_ends;/* TRUE if we erased them and need to reset */
  /* If we're in local to avoid extremely wide bands due to 
   * the permissive nature of local ends, we make local ends
   * impossible for the band calculation than make them
   * possible again before exiting this function.
   */
  reset_local_ends = reset_cp9_local_ends = FALSE;
  if(cm->flags & CM_LOCAL_END)
    {
      reset_local_ends = TRUE;
      if((cm->flags & CM_CP9) && cm->cp9->flags & CPLAN9_EL)
	reset_cp9_local_ends = TRUE;
      ConfigNoLocalEnds(cm);
    }
  /* gamma[v][n] is Prob(state v generates subseq of length n)
   */
  ESL_ALLOC(gamma, sizeof(double *) * cm->M);        
  for (v = 0; v < cm->M; v++) gamma[v] = NULL;

  /* dmin[v] and dmax[v] are the determined bounds that we return.
   */
  ESL_ALLOC(dmin, sizeof(int) * cm->M);
  ESL_ALLOC(dmax, sizeof(int) * cm->M);  
  if (ret_dmin != NULL) *ret_dmin = NULL;
  if (ret_dmax != NULL) *ret_dmax = NULL;

  /* beamstack is a trick for reusing memory: a pushdown stack of 
   * "beams" (gamma[v] rows) we can reuse.
   */
  beamstack = esl_stack_PCreate();

  /* The second component of memory saving is the "touch" array.
   * touch[y] is the number of states above state [y] that will
   * depend on y but haven't been calculated yet. When we're done
   * calculating a new state v, we decrement touch[y] for all
   * y \in C_v. Any time touch[y] reaches 0, we put that beam
   * back into the pool for reuse.
   */
  ESL_ALLOC(touch, sizeof(int) * cm->M);
  for (v = 0; v < cm->M; v++) touch[v] = cm->pnum[v];

  /* Allocate and initialize the shared end beam.
   */
  ESL_ALLOC(gamma[cm->M-1], (sizeof(double) * (W+1)));
  esl_vec_DSet(gamma[cm->M-1], W+1, 0.);
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

      if (esl_stack_PPop(beamstack, (void **) &gamma[v]) == eslEOD)
	ESL_ALLOC(gamma[v], sizeof(double) * (W+1));		    
      esl_vec_DSet(gamma[v], W+1, 0.);

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
      else if ((cm->flags & CM_LOCAL_BEGIN) && v == 0) /*state 0 is the one and only ROOT_S state*/
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
	      for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		{
		  y = cm->cfirst[v] + yoffset;
		  gamma[v][n] += cm->t[v][yoffset] * gamma[y][n-dv];
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
	  /*printf("truncation error unacceptable, failing.\n");
	    printf("p_thresh: %g\n", p_thresh);
	    printf("pdf : %g\n", pdf);
	    printf("gamma[v][W] : %g\n", gamma[v][W]);*/
	  status = 0; 
	  goto CLEANUP;
	}
      
      /* Renormalize this beam. (Should we really be doing this?)
       */
      if (pdf > 1.0) esl_vec_DNorm(gamma[v], W+1);

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
      if ((! save_densities) && (! (cm->flags & CM_LOCAL_BEGIN))) {
	if (cm->sttype[v] == B_st)
	  {  /* connected children of a B st are handled specially, remember */
	    y = cm->cfirst[v]; esl_stack_PPush(beamstack, gamma[y]); gamma[y] = NULL;
	    y = cm->cnum[v];   esl_stack_PPush(beamstack, gamma[y]); gamma[y] = NULL;
	  }
	else
	  {
	    for (y = cm->cfirst[v]; y < cm->cfirst[v]+cm->cnum[v]; y++)
	      {
		touch[y]--;
		if (touch[y] == 0 && cm->sttype[y] != E_st) {
		  esl_stack_PPush(beamstack, gamma[y]); 
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

  /* If we're in local mode, we set all local ends to impossible at
   * the beginning of this function, we set them back here.
   * Careful we can only call this once */
  if(reset_local_ends)
    {
      ConfigLocalEnds(cm, cm->pend);
      reset_local_ends = FALSE;
      if(reset_cp9_local_ends)
	{
	  CPlan9ELConfig(cm);
	  reset_cp9_local_ends = FALSE;
	}
    }
  if (! BandTruncationNegligible(gamma[0], dmax[0], W, NULL)) 
    { status = 0; goto CLEANUP; }

  status = 1;

 CLEANUP:
  free(touch);
  /* If we're in local mode, we set all local ends to impossible at
   * the beginning of this function, we set them back here.
   * Careful we can only call this once */
  if(reset_local_ends)
    {
      ConfigLocalEnds(cm, cm->pend);
      reset_local_ends = FALSE;
      if(reset_cp9_local_ends)
	{
	  CPlan9ELConfig(cm);
	  reset_cp9_local_ends = FALSE;
	}
    }
  if (ret_dmin  != NULL) *ret_dmin = dmin;   else free(dmin);
  if (ret_dmax  != NULL) *ret_dmax = dmax;   else free(dmax);
  if (ret_gamma != NULL) *ret_gamma = gamma; else FreeBandDensities(cm, gamma);
  
  /* Free the reused stack of beams.
   */
  while (esl_stack_PPop(beamstack, (void **) &tmp) != eslEOD) free(tmp);
  esl_stack_Destroy(beamstack);
  return status;

 ERROR:
  esl_fatal("Memory allocation error.\n");
  return eslFAIL; /* never reached */
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
 *
 * Note:     Dies (esl_fatal()) from memory allocation error, without
 *           cleanup.
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
  char         *name;           /* name for the seq we've emitted */
  ESL_RANDOMNESS  *r = NULL;    /* source of randomness */

  /* Create and seed RNG */
  if ((r = esl_randomness_CreateTimeseeded()) == NULL) 
    esl_fatal("Failed to create random number generator: probably out of memory");

  /* Allocate gamma, completely; and initialize to zeros. 
   * For consistency w/ BandCalculationEngine(), allocate a single
   * shared end deck at M-1, and point other ends at it - even
   * though P(n=0) = 1.0 by definition at the E's and we don't
   * really need to calculate it. Then we can use FreeBandDensities()
   * for gamma matrices alloc'ed in either function.
   */                                                     
  ESL_ALLOC(gamma, (sizeof(double *) * cm->M));
  ESL_ALLOC(gamma[cm->M-1], (sizeof(double) * (W+1))); 
  esl_vec_DSet(gamma[cm->M-1], W+1, 0.0);
  for (v = 0; v < cm->M-1; v++)
    {
      if (cm->sttype[v] != E_st)
	{
	  ESL_ALLOC(gamma[v], sizeof(double) * (W+1));
	  esl_vec_DSet(gamma[v], W+1, 0.0);
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
      sprintf(name, "seq%d", i+1);
      EmitParsetree(cm, r, NULL, FALSE, &tr, NULL, &seqlen);
      free(name);
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
  esl_randomness_Destroy(r);
  return status;

 ERROR:
  esl_fatal("Memory allocation error.\n");
  return eslFAIL; /* never reached */
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
 *           dsq       - the digitized sequence
 *           dmin      - minimum bound on d for state v; 0..M
 *           dmax      - maximum bound on d for state v; 0..M          
 *           i0        - start of target subsequence (1 for full seq)
 *           j0        - end of target subsequence (L for full seq)
 *           W         - max d: max size of a hit
 *           cutoff    - minimum score to report 
 *           results   - search_results_t to add to; if NULL, don't add to it
 *
 * Returns:  score of best overall hit
 *
 * Note:     Dies (esl_fatal()) from memory allocation error, without
 *           cleanup.
 */
float
CYKBandedScan(CM_t *cm, ESL_DSQ *dsq, int *dmin, int *dmax, int i0, int j0, int W, 
	      float cutoff, search_results_t *results)
{
  printf("in CYKBandedScan\n");
  /* Contract check */
  if(j0 < i0)
    esl_fatal("ERROR in CYKBandedScan, i0: %d j0: %d\n", i0, j0);
  if(!(cm->flags & CM_QDB))
    esl_fatal("ERROR in CYKBandedScan, QDBs invalid\n");
  if(dsq == NULL)
    esl_fatal("ERROR in CYKBandedScan, dsq is NULL.\n");

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
  int       jp;			/* rolling index into BEGL_S decks: jp=j%(W+1) */
  int       jmax;               /* when imposing bands, maximum j value in alpha matrix */
  int       kmax;               /* for B_st's, maximum k value consistent with bands*/
  int       L;                  /* length of the subsequence (j0-i0+1) */
  int       gamma_j;            /* j index in the gamma matrix, which is indexed 0..j0-i0+1, 
				 * while j runs from i0..j0 */
  int       gamma_i;            /* i index in the gamma* data structures */
  int       curr_dmin;          /* temporary value for min d in for loops */
  int       curr_dmax;          /* temporary value for max d in for loops */
  float     best_score;         /* Best overall score from semi-HMM to return */
  float     best_neg_score;     /* Best score overall score to return, used if all scores < 0 */
  int       bestd;              /* d value of best hit thus far seen for j (used if greedy strategy) */

  /*PrintDPCellsSaved(cm, dmin, dmax, W);*/
  best_score     = IMPOSSIBLE;
  best_neg_score = IMPOSSIBLE;
  L = j0-i0+1;
  if (W > L) W = L; 

  /*printf("in CYKBandedScan i0: %d j0: %d\n", i0, j0);*/

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
	for (j = 2; j <= W; j++) 
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
	  /* determine max d we're allowing for this state v and this position j */
	  curr_dmax = (gamma_j < dmax[v]) ? gamma_j : dmax[v];
	  if(curr_dmax > W) curr_dmax = W;
	  curr_dmin = (dmin[v] > 0) ? dmin[v] : 1;
	  if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	    {
	      if (cm->stid[v] == BEGL_S) jp = j % (W+1); else jp = cur;
	      for (d = curr_dmin; d <= curr_dmax; d++) 
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
	      for (d = curr_dmin; d <= curr_dmax; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][prv][d-2] + cm->tsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  i = j-d+1;
		  if (dsq[i] < cm->abc->K && dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][(int) (dsq[i]*cm->abc->K+dsq[j])];
		  else
		    alpha[v][cur][d] += DegeneratePairScore(cm->abc, cm->esc[v], dsq[i], dsq[j]);
		  
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) 
	    {
	      for (d = curr_dmin; d <= curr_dmax; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][cur][d-1] + cm->tsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  i = j-d+1;
		  if (dsq[i] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][(int) dsq[i]];
		  else
		    alpha[v][cur][d] += esl_abc_FAvgScore(cm->abc, dsq[i], cm->esc[v]);
		  
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) 
	    {
	      for (d = curr_dmin; d <= curr_dmax; d++)
		{
		  y = cm->cfirst[v];
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d-StateDelta(cm->sttype[v])));
		  for (yoffset = 0; yoffset < cm->cnum[v]; yoffset++)
		    if ((sc = alpha[y+yoffset][prv][d-1] + cm->tsc[v][yoffset]) > alpha[v][cur][d])
		      alpha[v][cur][d] = sc;
		  if (dsq[j] < cm->abc->K)
		    alpha[v][cur][d] += cm->esc[v][(int) dsq[j]];
		  else
		    alpha[v][cur][d] += esl_abc_FAvgScore(cm->abc, dsq[j], cm->esc[v]);
		  
		  if (alpha[v][cur][d] < IMPROBABLE) alpha[v][cur][d] = IMPOSSIBLE;
		}
	    }
	  else if (cm->sttype[v] == B_st) 
	    {
	      w = cm->cfirst[v];
	      y = cm->cnum[v];
	      i = j-d+1;
	      for (d = curr_dmin; d <= curr_dmax; d++) 
		{
		  alpha[v][cur][d] = cm->endsc[v] + (cm->el_selfsc * (d - StateDelta(cm->sttype[v])));
		  /* k is the length of the right fragment */
		  /* Careful, make sure k is consistent with bands in state w and state y. */
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

      /* determine max d we're allowing for the root state and this position j */
      curr_dmax = (gamma_j < dmax[0]) ? gamma_j : dmax[0];
      if(curr_dmax > W) curr_dmax = W;
      curr_dmin = (dmin[0] > 0) ? dmin[0] : 1;

      for (d = curr_dmin; d <= curr_dmax; d++)
	{
	  y = cm->cfirst[0];
	  alpha[0][cur][d] = alpha[y][cur][d] + cm->tsc[0][0];
	  bestr[d]         = 0;	/* root of the traceback = root state 0 */
	  for (yoffset = 1; yoffset < cm->cnum[0]; yoffset++)
	    if ((sc = alpha[y+yoffset][cur][d] + cm->tsc[0][yoffset]) > alpha[0][cur][d]) 
	      alpha[0][cur][d] = sc;
	  if (alpha[0][cur][d] < IMPROBABLE) alpha[0][cur][d] = IMPOSSIBLE;
	  if (alpha[0][cur][d] > best_neg_score) best_neg_score = alpha[0][cur][d];
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
	  for (; (d <= dmax[y] && d <= gamma_j) && d <= W; d++)
	    {
	      if (cm->stid[y] == BEGL_S) sc = alpha[y][j%(W+1)][d] + cm->beginsc[y];
	      else                       sc = alpha[y][cur][d]     + cm->beginsc[y];
	      if (sc > alpha[0][cur][d]) {
		alpha[0][cur][d] = sc;
		bestr[d]         = y;
	      }
	      if (alpha[0][cur][d] < IMPROBABLE) alpha[0][cur][d] = IMPOSSIBLE;
	      if (alpha[0][cur][d] > best_neg_score) best_neg_score = alpha[0][cur][d];
	    }
	}
      }
      
      if(!(cm->search_opts & CM_SEARCH_CMGREEDY)) /* resolve overlaps optimally */
	{
	  /* The little semi-Markov model that deals with multihit parsing:
	   */
	  gamma[gamma_j]  = gamma[gamma_j-1] + 0; /* extend without adding a new hit */
	  gback[gamma_j]  = -1;
	  savesc[gamma_j] = IMPOSSIBLE;
	  saver[gamma_j]  = -1;
	  for (d = dmin[0]; d <= curr_dmax; d++) 
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
	  if (alpha[0][cur][dmin[v]] >= cutoff) 
	    if(results != NULL) 
	      report_hit (j, j, bestr[1], alpha[0][cur][dmin[v]], results);
	  bestd = 1;
	  if (alpha[0][cur][dmin[v]] > best_score) 
	    best_score = alpha[0][cur][dmin[v]];
	  
	  /* Now, if current score is greater than maximum seen previous, report
	     it if >= cutoff and set new max */
	  for (d = dmin[v]; d <= curr_dmax; d++) 
	    {
	      if (alpha[0][cur][d] > best_score) 
		best_score = alpha[0][cur][d];
	      if (alpha[0][cur][d] > alpha[0][cur][bestd]) 
		{
		  if (alpha[0][cur][d] >= cutoff)
		    if(results != NULL) 
		      report_hit (j-d+1, j, bestr[d], alpha[0][cur][d], results);
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
  return best_score;

 ERROR:
  esl_fatal("Memory allocation error.\n");
  return 0.; /* never reached */
}

/* EPN 07.22.05
 * ExpandBands()
 * Function: ExpandBands
 *
 * Purpose:  Called when the sequence we are about to align 
 *           using bands is either shorter in length than
 *           the dmin on the root state, or longer in length
 *           than the dmax on the root state.
 *            
 *           This function expands the bands on ALL states
 *           v=1..cm->M-1 in the following manner :
 *           
 *           case 1 : target len < dmin[0]
 *                    subtract (dmin[0]-target len) from
 *                    dmin of all states, and ensure
 *                    dmin[v]>=0 for all v.
 *                    Further :
 *                    if cm->sttype[v] == MP_st ensure dmin[v]>=2;
 *                    if cm->sttype[v] == IL_st || ML_st ensure dmin[v]>=1;
 *                    if cm->sttype[v] == IR_st || MR_st ensure dmin[v]>=1;
 *                        
 *           case 2 : target len > dmax[0]
 *                    add (target len-dmax[0] to dmax
 *                    of all states.
 *
 *           Prior to handling such situtations with this
 *           hack, the program would choke and die.  This
 *           hacky approach is used as a simple, inefficient
 *           not well thought out, but effective way to 
 *           solve this problem.
 * 
 * Args:    cm       - the CM
 *          tlen     - length of target sequence about to be aligned
 *          dmin     - minimum d bound for each state v; [0..v..M-1]
 *                     may be modified in this function
 *          dmax     - maximum d bound for each state v; [0..v..M-1]
 *                     may be modified in this function
 *
 * Returns: (void) 
 */

void
ExpandBands(CM_t *cm, int tlen, int *dmin, int *dmax)
{
  int v;
  int diff;
  int root_min;
  int root_max;
  int M = cm->M;
  root_min = dmin[0];
  root_max = dmax[0];

  if(tlen < root_min)
    {
      diff = root_min - tlen;
      for(v=0; v<M; v++)
	{
	  dmin[v] -= diff;
	  if((cm->sttype[v] == MP_st) && (dmin[v] < 2)) 
	    dmin[v] = 2;
	  else if(((cm->sttype[v] == IL_st) || (cm->sttype[v] == ML_st)) 
		  && (dmin[v] < 1)) 
	    dmin[v] = 1;
	  else if(((cm->sttype[v] == IR_st) || (cm->sttype[v] == MR_st)) 
		  && (dmin[v] < 1)) 
	    dmin[v] = 1;
	  else if(dmin[v] < 0) 
	    dmin[v] = 0;
	}
    }
  else if(tlen > root_max)
    {
      diff = tlen - root_min;
      for(v=0; v<M; v++)
	{
	  dmax[v] += diff;
	}
    }
}

/* EPN 08.15.05
 * qdb_trace_info_dump()
 * Function: qdb_trace_info_dump
 *
 * Purpose:  Called when the user has enabled the --banddump
 *           options.  This function determines how close the
 *           trace was to the bands at each state in the trace,
 *           and prints out that information in differing levels
 *           of verbosity depending on an input parameter 
 *           (bdump_level).
 * 
 * Args:    tr       - the parsetree (trace)
 *          dmin     - minimum d bound for each state v; [0..v..M-1]
 *                     may be modified in this function
 *          dmax     - maximum d bound for each state v; [0..v..M-1]
 *                     may be modified in this function
 *          bdump_level - level of verbosity
 * Returns: (void) 
 */

void
qdb_trace_info_dump(CM_t *cm, Parsetree_t *tr, int *dmin, int *dmax, int bdump_level)
{
  int    status;
  char **sttypes;
  char **nodetypes;
  int v, i, j, d, tpos;
  int mindiff;            /* d - dmin[v] */
  int maxdiff;            /* dmax[v] - d */

  ESL_ALLOC(sttypes, sizeof(char *) * 10);
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

  ESL_ALLOC(nodetypes, sizeof(char *) * 8);
  nodetypes[0] = "BIF";
  nodetypes[1] = "MATP";
  nodetypes[2] = "MATL";
  nodetypes[3] = "MATR";
  nodetypes[4] = "BEGL";
  nodetypes[5] = "BEGR";
  nodetypes[6] = "ROOT";
  nodetypes[7] = "END";

  for (tpos = 0; tpos < tr->n; tpos++)
    {
      v  = tr->state[tpos];
      i = tr->emitl[tpos];
      j = tr->emitr[tpos];
      d = j-i+1;

      if(cm->sttype[v] != EL_st)
	{
	  mindiff = d-dmin[v];
	  maxdiff = dmax[v]-d;
	  if(bdump_level > 1 || ((mindiff < 0) || (maxdiff < 0)))
	    printf("%-4s %-3s v: %4d | d: %4d | dmin: %4d | dmax: %4d | %3d | %3d |\n", nodetypes[(int) cm->ndtype[(int) cm->ndidx[v]]], sttypes[(int) cm->sttype[v]], v, d, dmin[v], dmax[v], mindiff, maxdiff);
	}
      else
	{
	  if(bdump_level > 1)
	    printf("%-8s v: %4d | d: %4d |\n", sttypes[(int) cm->sttype[v]], v, d);
	}
    }
  free(sttypes);
  free(nodetypes);
  return;

 ERROR:
  esl_fatal("Memory allocation error.");
}


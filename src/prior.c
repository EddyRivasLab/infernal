/* Support for Dirichlet-mixture priors. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "squidconf.h"
#include "config.h"

#include "squid.h"
#include "vectorops.h"

#include "sre_stack.h"
#include "structs.h"
#include "funcs.h"

/* Function: LogNorm()
 * 
 * Purpose:  Normalize a vector of log likelihoods, changing it
 *           to a probability vector. Be careful of overflowing exp().
 *           Implementation adapted from Graeme Mitchison.
 *
 * Args:     vec - vector destined to become log probabilities
 *           n   - length of vec 
 *
 * NOTE : EPN - this function probably exists in squid under a different
 *              name, but I couldn't find it.  This is copied from the
 *              mathsupport.c file in HMMER 2.3.2
 */
void
LogNorm(float *vec, int n)
{
  int   x;
  float max   = -1.0e30;
  float denom = 0.;

  for (x = 0; x < n; x++)
    if (vec[x] > max) max = vec[x];
  for (x = 0; x < n; x++)
    if (vec[x] > max - 50.)
      denom += exp(vec[x] - max);
  for (x = 0; x < n; x++)
    if (vec[x] > max - 50.)
      vec[x] = exp(vec[x] - max) / denom;
    else
      vec[x] = 0.0;
}
 
/* Function: AllocPrior(), FreePrior()
 * 
 * Purpose:  Allocation and free'ing of a prior structure.
 *           Very simple, but might get more complex someday.
 */
struct prior_s *
AllocPrior(void)
{ return (struct prior_s *) MallocOrDie (sizeof(struct prior_s)); }
void
FreePrior(struct prior_s *pri)
{ free(pri); }


/* Function: ReadPrior()
 * 
 * Purpose:  Input a transition prior from disk file.
 */
struct prior_s *
ReadPrior(char *prifile) 
{
  FILE             *fp;
  struct prior_s *pri;
  char             *buf;
  char             *tok;
  char             *s;
  int               n;
  int               toklen;
  int               q, x;
  int               i;       /*counter over transition sets*/
  int               j;       /*counter over components in a mixture*/
  int               k;       /*counter over alphas*/
  int               curr_state_id; 
  int               curr_next_node_id; 
  int               a, b;

  if ((fp = fopen(prifile, "r")) == NULL)
    Die("Failed to open INFERNAL prior file %s\n", prifile);
  pri = AllocPrior();

  /*set some defaults*/
  pri->mbpasize = Alphabet_size * Alphabet_size;
  pri->mntasize = Alphabet_size;
  pri->iasize = Alphabet_size;
  
  for(a = 0; a < UNIQUESTATES; a++)
    {
      for(b = 0; b < NODETYPES; b++)
	{
	  pri->tsetmap[a][b] = -1;
	}
    }

  
  /***********************************************************
   * Format of Dirichlet file : 
   * 
   * line 1 : "DIRICHLET"
   * line 2 : number of transition sets T
   * 
   * next section repeated T times (once for each transition set)
   *   line 1 : <node state ID number> <next_node_ID number>
   *   line 2 : <alphabet size for current transition set>
   *   line 3 : <number of components> 
   *            for each component : 
   *            <mixture coefficients> (separated by a space)
   *            <alphas> (sep by space)
   *   
   * emission parameters : 
   * first, match base pair priors : 
   *   line 1 : <number of components> 
   *            for each component : 
   *            <mixture coefficients> (separated by a space)
   *            <alphas> (sep by space)
   *
   * second, match singlet priors :
   *   line 1 : <number of components> 
   *            for each component : 
   *            <mixture coefficients> (separated by a space)
   *            <alphas> (sep by space)
   *
   * finally, insertion priors :
   *   line 1 : <number of components> 
   *            for each component : 
   *            <mixture coefficients> (separated by a space)
   *            <alphas> (sep by space)
   *
   *     


  /* First entry is the strategy: 
   * Only standard Dirichlet prior (simple or mixture) is supported so far
   */
  //  sptr = Getword(fp, sqdARG_STRING);
  n = 0;
  buf = NULL;
  
  sre_fgets(&buf, &n, fp);
  if      (strncmp(buf, "DIRICHLET", 9) == 0) pri->strategy = PRI_DCHLET;
  else Die("No such prior strategy %s; failed to parse file %s", s, prifile);
 
  pri->tsetnum = atoi(sre_fgets(&buf, &n, fp));
  if (pri->tsetnum < 0)
    Die("%d is bad; need at least one transition set", pri->tsetnum);
  if (pri->tsetnum > MAXTRANSSETS)
    Die("%d is bad, too many transition sets (MAXTRANSSET = %d)\n", MAXTRANSSETS);
  
  /*
   * transition parameters : 
   * next section repeated pri->tsetnum times (once for each transition set)
   *   line 1 : <node state ID number> <next_node_ID number>
   *   line 2 : <alphabet size for current transition set>
   *   line 3 : <number of components> 
   *            for each component : 
   *            <mixture coefficients> (separated by a space)
   *            <alphas> (sep by space)
   */
  
  for(i = 0; i < pri->tsetnum; i++)
    {
      sre_fgets(&buf, &n, fp);
      s = buf;
      if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) Die("ERROR : (A1) reading in transitions in prior file\n");
      curr_state_id = atoi(tok);
      //      printf("current state id is %d\n", curr_state_id);
      if(curr_state_id > UNIQUESTATES) Die("ERROR : (A2) reading in transitions in prior file\n");
      if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) Die("Wrong format of prior file\n");

      curr_next_node_id = atoi(tok);
      if(curr_next_node_id > NODETYPES) Die("ERROR : (A3) reading in transitions in prior file\n");
      //printf("current next node id is %d\n", curr_next_node_id);
      
      // add information to setmap
      pri->tsetmap[curr_state_id][curr_next_node_id] = i;

      //get alphabet size for current transition set
      sre_fgets(&buf, &n, fp);
      s = buf;
      if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) Die("Wrong format of prior file\n");
      //printf("A4 tok is %s\n", tok);
      pri->tasize[i] = atoi(tok);
      if(pri->tasize[i] > MAXABET) Die("ERROR : (A4) reading in transitions in prior file\nalph size is %d\n", pri->tasize[i]);

      //get number of components for current transition set
      sre_fgets(&buf, &n, fp);
      //printf("A5 buf is %s\n", buf);
      pri->tnum[i] = atoi(buf);
      if(pri->tnum[i] > MAXDCHLET) Die("ERROR : (A5) reading in transitions in prior file\n");

      for(j = 0; j < pri->tnum[i]; j++)
	{
	  //get mixture coefficient for current transition set i
	  //and current component j
	  sre_fgets(&buf, &n, fp);
	  //printf("A6 buf is %s\n", buf);
	  pri->tq[i][j] = (float) atof(buf);
	  if(pri->tq[i][j] > 1.0) Die("ERROR : (A6) reading in transitions in prior file\npri->tq[%d][%d] is %d\n", i, j, pri->tq[i][j]);

	  //get alphas
	  sre_fgets(&buf, &n, fp);
	  s = buf;
	  for(k = 0; k < pri->tasize[i]; k++)
	    {
	      if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) Die("ERROR : (A7) reading in transitions in prior file\nalph size is %d\n", pri->tasize[i]);
	      pri->t[i][j][k] = (float) atof(tok);
	      //	      printf("reading in transition priors\n");
	      //printf("tok is %s\n", tok);
	      //printf("pri->[%d][%d][%d] is %f\n", i, j, k, pri->t[i][j][k]);
	    }
	}
    }

  /*
   * emission parameters : 
   * first, match base pair priors : 
   *            <number of components> 
   *            for each component : 
   *            <mixture coefficients> (separated by a space)
   *            <alphas> (sep by space)
   */
  
  //get number of components for match base pair prior
  sre_fgets(&buf, &n, fp);
  pri->mbpnum = atoi(buf);
  if(pri->mbpnum > MAXDCHLET) Die("ERROR : (B1) reading in emission bps\n");
  for(j = 0; j < pri->mbpnum; j++)
    {
      //get mixture coefficient for current transition set i
      //and current component j
      sre_fgets(&buf, &n, fp);
      pri->mbpq[j] = (float) atof(buf);
      if(pri->mbpq[j] > 1.0) Die("ERROR : (B2) reading in emission bps\n");
      
      //get alphas
      sre_fgets(&buf, &n, fp);
      s = buf;
      for(k = 0; k < pri->mbpasize; k++)
	{
	  if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) Die("ERROR : (B3) reading in emission bps\n");
	  pri->mbp[j][k] = (float) atof(tok);
	  printf("read bp prior mbp[%d][%d] %f\n", j, k, pri->mbp[j][k]);
	}
    }
  
  /*
   * second, match singlet priors :
   *            <number of components> 
   *            for each component : 
   *            <mixture coefficients> (separated by a space)
   *            <alphas> (sep by space)
   */

  //get number of components for match singlet prior
  sre_fgets(&buf, &n, fp);
  pri->mntnum = atoi(buf);
  if(pri->mntnum > MAXDCHLET) Die("ERROR : (C1) reading in emission nts\n");
  for(j = 0; j < pri->mntnum; j++)
    {
      //get mixture coefficient for current transition set i
      //and current component j
      sre_fgets(&buf, &n, fp);
      pri->mntq[j] = (float) atof(buf);
      if(pri->mntq[j] > 1.0) Die("ERROR : (C2) reading in emission nts\n");
      
      //get alphas
      sre_fgets(&buf, &n, fp);
      s = buf;
      for(k = 0; k < pri->mntasize; k++)
	{
	  if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) Die("ERROR : (C3) reading in emission nts\n");
	  pri->mnt[j][k] = (float) atof(tok);
	}
    }
  
  /*
   * finally, insertion priors :
   *   line 1 : <number of components> 
   *            for each component : 
   *            <mixture coefficients> (separated by a space)
   *            <alphas> (sep by space)
   */

  //get number of components for match singlet prior
  sre_fgets(&buf, &n, fp);
  pri->inum = atoi(buf);
  if(pri->inum > MAXDCHLET) Die("ERROR : (E1) reading in emission inserts\n");
  
  for(j = 0; j < pri->inum; j++)
    {
      //get mixture coefficient for current transition set i
      //and current component j
      sre_fgets(&buf, &n, fp);
      pri->iq[j] = (float) atof(buf);
      if(pri->iq[j] > 1.0) Die("ERROR : (E2) reading in emission inserts\n");
      
      //get alphas
      sre_fgets(&buf, &n, fp);
      s = buf;
      for(k = 0; k < pri->iasize; k++)
	{
	  if ((tok = sre_strtok(&s, " \t\n", &toklen)) == NULL) Die("ERROR : (E3) reading in emission inserts\n");
	  pri->i[j][k] = (float) atof(tok);
	}
    }

  printf("LEAVING READPRIOR\n");

  fclose(fp);
  return pri;

 FAILURE:
  Die("Error with the prior file.\n");

}

/* Function: PriorifyCM()
 * 
 * Purpose:  Add pseudocounts to a CM using Dirichlet priors,
 *           and renormalize the CM.
 * 
 * Args:     CM -- the CM to add counts to (counts form)
 *           pri -- the Dirichlet prior to use
 *           
 * Return:   (void)
 *           CM returns in probability form.
 */          
void
PriorifyCM(struct cm_s *cm, struct prior_s *pri)
{
  int v;			/* counter for model position   */
  float d;			/* a denominator */
  int setnum;                   /* number of set to use */
  int nxtndtype;                /* type of next node */
  int nselfndtrans;             /* number of 'self-node' transitions,
				   defined by number of transitions
				   from state v to another state that
				   is a member of the same node as v. */
                                 

    /* EPN - No analog for below commented out code in CMs b/c 
     * begin and end states have no transition priors 
     * (I think this is right.)
     */

  /* Model-dependent transitions are handled simply; Laplace.
   */
    //  FSet(hmm->begin+2, hmm->M-1, 0.);     /* wipe internal BM entries */
    //FSet(hmm->end+1, hmm->M-1, 0.);	/* wipe internal ME exits   */
    //d = hmm->tbd1 + hmm->begin[1] + 2.;
    //hmm->tbd1        = (hmm->tbd1 + 1.)/ d;
    //hmm->begin[1]    = (hmm->begin[1] + 1.)/ d;
    //hmm->end[hmm->M] = 1.0;

  /* Main model transitions and emissions
   */
  for (v = 0; v < cm->M; v++)
    {
      /* Priorify transition vector if not a BIF or E state */
      if(cm->sttype[v] != B_st && cm->sttype[v] != E_st)
	{
	  /* Determine which transition set to use. 
	   * We already know the current state id
	   * but we need to know the type of the next node. First
	   * we need to know how many transitions from the current state
	   * v to other states in the same node there are, which
	   * depends on the current node type.  Then we 
	   * can figure out the type of the next node.
	   */

	  if(cm->ndtype[cm->ndidx[v]] == MATP_nd || 
	     cm->ndtype[cm->ndidx[v]] == ROOT_nd)
	    {
	      if(cm->sttype[v] != IR_st)
		{
		  nselfndtrans = 2;
		}
	      else
		{
		  nselfndtrans = 1;
		}
	    } 
	  else if(cm->ndtype[cm->ndidx[v]] == MATL_nd || 
		  cm->ndtype[cm->ndidx[v]] == MATR_nd || 
		  cm->ndtype[cm->ndidx[v]] == BEGR_nd)
	    {
	      nselfndtrans = 1;
	    }
	  else if(cm->ndtype[cm->ndidx[v]] == BEGL_nd)
	    {
	      nselfndtrans = 0;
	    }
	  
	  nxtndtype = cm->ndtype[cm->ndidx[cm->cfirst[v] + nselfndtrans]];
	  setnum = pri->tsetmap[cm->stid[v]][nxtndtype];

	  /* Some debugging print statements 
	  printf("v is %d\nstid is %d\nnext node type is %d\n\n", v, cm->stid[v], nxtndtype);
	  printf("nselfndtrans is %d\n", nselfndtrans);
	  printf("cm->cfirst[%d] is %d\n", v, cm->cfirst[v]);
	  printf("current node index is %d\n", cm->ndidx[v]);
	  printf("current node type is %d\n", cm->ndtype[cm->ndidx[v]]);
	  printf("cm->ndidx[%d] is %d\n", (cm->cfirst[v] + nselfndtrans), cm->ndidx[(cm->cfirst[v] + nselfndtrans)]);
	  printf("setnum is %d\n", setnum);
	  */

	  PriorifyTransitionVector(cm->t[v], pri, pri->tq[setnum], setnum);
	}
      
      /* Now add in the emission priors */
      if(cm->sttype[v] == MP_st)
	{
	  PriorifyEmissionVector(cm->e[v], pri, pri->mbpnum, pri->mbpasize, pri->mbpq, pri->mbp, NULL);
	}
      else if (cm->sttype[v] == ML_st || cm->sttype[v] == MR_st)
	{
	  PriorifyEmissionVector(cm->e[v], pri, pri->mntnum, pri->mntasize, pri->mntq, pri->mnt, NULL);
	}
      else if (cm->sttype[v] == IL_st || cm->sttype[v] == IR_st)
	{
	  PriorifyEmissionVector(cm->e[v], pri, pri->inum, pri->iasize, pri->iq, pri->i, NULL);
	}

    }
}

/* Function: PriorifyEmissionVector()
 * 
 * Purpose:  Add prior pseudocounts to an observed 
 *           emission count vector and renormalize. 
 *
 *           Can return the posterior mixture probabilities
 *           P(q | counts) if ret_mix[MAXDCHLET] is passed.
 *           Else, pass NULL.  
 * 
 * Args:     vec     - the 4 or 20-long vector of counts to modify
 *           pri     - prior data structure
 *           num     - pri->mbpnum, pri->mntnum or pri->inum; # of mixtures
 *           asize   - pri->mbpasize (16), pri->mntasize (4), or pri->iasize (4)
 *           eq      - pri->mbpq, pri->mntq or pri->iq; prior mixture probabilities
 *           e       - pri->i or pri->mnt, or pri->mbp; Dirichlet components          
 *           ret_mix - filled with posterior mixture probabilities, or NULL
 *                   
 * Return:   (void)
 *           The counts in vec are changed and normalized to probabilities.
 *
 * NOTE : Copied and morphed from HMMER 2.3.2
 */                  
void
PriorifyEmissionVector(float *vec, struct prior_s *pri, 
		       int num, int asize, float eq[MAXDCHLET], float e[MAXDCHLET][(MAXABET * MAXABET)],
		       float *ret_mix)
{
  int   x;                      /* counter over vec                     */
  int   q;                      /* counter over mixtures                */
  float mix[MAXDCHLET];         /* posterior distribution over mixtures */
  float totc;                   /* total counts                         */
  float tota;                   /* total alpha terms                    */
  float xi;                     /* X_i term, Sjolander eq. 41           */

  /* Calculate mix[], which is the posterior probability
   * P(q | n) of mixture component q given the count vector n
   *
   * (side effect note: note that an insert vector in a PAM prior
   * is passed with num = 1, bypassing pam prior code; this means
   * that inserts cannot be mixture Dirichlets...)
   * [SRE, 12/24/00: the above comment is cryptic! what the hell does that
   *  mean, inserts can't be mixtures? doesn't seem to be true. it 
   *  may mean that in a PAM prior, you can't have a mixture for inserts,
   *  but I don't even understand that. The insert vectors aren't passed
   *  with num=1!!]
   */
  printf("in PriorifyEmissionVector num is %d\n", num);
  for(q = 0; q < num; q++)
    {
      printf("vec %d is %f\n", q, vec[q]);
      printf("eq %d is %f\n", q, eq[q]);
    }
  printf("Alphabet size is %d\n\n", asize);


  mix[0] = 1.0;
  if (pri->strategy == PRI_DCHLET && num > 1) 
    {
      for (q = 0; q < num; q++) 
	{
	  mix[q] =  eq[q] > 0.0 ? log(eq[q]) : -999.;
	  mix[q] += Dchlet_logp_counts(vec, asize, e[q]);
	}
      LogNorm(mix, num);      /* now mix[q] is P(component_q | n) */
    }
  else if (pri->strategy == PRI_PAM && num > 1) 
    {		/* pam prior uses aa frequencies as `P(q|n)' */
      for (q = 0; q < asize; q++) 
	mix[q] = vec[q];
      FNorm(mix, asize);
    }

  /* Convert the counts to probabilities, following Sjolander (1996) 
   */
  totc = FSum(vec, asize);
  for (x = 0; x < asize; x++) {
    xi = 0.0;
    for (q = 0; q < num; q++) {
      tota = FSum(e[q], asize);
      xi += mix[q] * (vec[x] + e[q][x]) / (totc + tota);
    }
    vec[x] = xi;
  }
  FNorm(vec, asize);

  if (ret_mix != NULL)
    for (q = 0; q < num; q++)
      ret_mix[q] = mix[q];
}



/* Function: PriorifyTransitionVector()
 * 
 * Purpose:  Add prior pseudocounts to transition vector,
 *           
 * Args:     t     - state transitions, counts   
 *           pri   - Dirichlet prior information
 *           tq    - prior distribution over Dirichlet components.
 *                   (overrides pri->tq[]; used for alternative
 *                   methods of conditioning prior on structural data)  
 *           setnum - the transition set number
 * Return:   (void)
 *           t is changed, and renormalized -- comes back as
 *           probability vectors.
 * NOTE : Copied and morphed from HMMER 2.3.2
 */          
void
PriorifyTransitionVector(float *t, struct prior_s *pri, 
			   float tq[MAXDCHLET], int setnum)
{
  int   ts;
  int   q;
  float mix[MAXDCHLET];
  float totc;                   /* total counts */
  float tota;                   /* alpha terms */
  float xi;                     /* Sjolander's X_i term */
  int   x;                      /* counter over t */

  mix[0] = 1.0;			/* default is simple one component */
  if ((pri->strategy == PRI_DCHLET || pri->strategy == PRI_PAM) && pri->tnum[setnum] > 1)
    {
      for (q = 0; q < pri->tnum[setnum]; q++)
        {
          mix[q] =  tq[q] > 0.0 ? log(tq[q]) : -999.;
	  mix[q] += Dchlet_logp_counts(t, pri->tasize[setnum], pri->t[setnum][q]);
	  //LINES BELOW ARE FROM HMMER CODE
	  //mix[q] += Dchlet_logp_counts(t,   3, pri->t[setnum][q]);   /* 3 match  */
          //mix[q] += Dchlet_logp_counts(t+3, 2, pri->t[setnum][q]+3); /* 2 insert */
	  //mix[q] += Dchlet_logp_counts(t+5, 2, pri->t[q]+5); /* 2 delete */
        }
      LogNorm(mix, pri->tnum[setnum]); /* mix[q] is now P(q | counts) */
    }
				/* precalc some denominators */

  /*Code chunk below from HMMER 2.3.2 because a transition
   *vector passed into this function is in fact 3 transition vectors
   *for the 3 states in each node - Match, Insert, and Delete
   *This is not the current strategy here for INFERNAL, but might be useful if
   *that strategy changes 
   */
  /*  totm = FSum(t,3);		
  toti = t[TIM] + t[TII];
  totd = t[TDM] + t[TDD];

  for (ts = 0; ts < 7; ts++)  
    {
      xi = 0.0;
      for (q = 0; q < pri->tnum; q++)
        {
	  switch (ts) {
	  case TMM: case TMI: case TMD: 
	    xi += mix[q] * (t[ts] + pri->t[q][ts]) / 
	      (totm + FSum(pri->t[q], 3)); 
	    break;
	  case TIM: case TII: 
	    xi += mix[q] * (t[ts] + pri->t[q][ts]) / 
	      (toti + pri->t[q][TIM] + pri->t[q][TII]);
	    break;
	  case TDM: case TDD: 
	    xi += mix[q] * (t[ts] + pri->t[q][ts]) / 
	      (totd + pri->t[q][TDM] + pri->t[q][TDD]);
	    break;
	  }
        }
      t[ts] = xi;
      }*/
  //  FNorm(t,   3);		/* match  */
  //FNorm(t+3, 2);		/* insert */
  //FNorm(t+5, 2);		/* delete */

  //Chunk below copied and modified from ProbifyEmissionVector()
  /* Convert the counts to probabilities, following Sjolander (1996) 
   */

  /* Debugging print statements */
  //printf("\nin priorify transition vec setnum is %d\n", setnum);
  //printf("before incorporating priors\n");
  //  for(x = 0; x < pri->tasize[setnum]; x++) 
  //  {
  //    printf("counts t[%d] is %f\n", x, t[x]);
  //  }
  //printf("\nafter incorporating priors\n");
  
  totc = FSum(t, pri->tasize[setnum]);
  float total = 0;
  for (x = 0; x < pri->tasize[setnum]; x++) {
    xi = 0.0;
    for (q = 0; q < pri->tnum[setnum]; q++) {
      tota = FSum(pri->t[setnum][q], pri->tasize[setnum]);
      //      printf("pri->t[%d][%d][%d] is %f\n", setnum, q, x, pri->t[setnum][q][x]);
      xi += mix[q] * (t[x] + pri->t[setnum][q][x]) / (totc + tota);
      assert(!isnan(xi));
    }
    t[x] = xi;
    //printf("new counts t[%d] is %f\n", x, t[x]);

    total += xi;
    //    printf("%1.3f\n", xi);    
  }
  //  printf("\n");
  assert(abs(total-1) < 1e-6); // Check that we actually have probabilities.
  FNorm(t, pri->tasize[setnum]);
  for (x = 0; x < pri->tasize[setnum]; x++) {
    //printf("normalized t[%d] is %f\n", x, t[x]);
    //printf("about to assert t[%d] = %f is > 0\n", x, t[x]);
    assert(t[x] > 0);
  }  
}


/* Below are Alex's default prior structures, but they are not consistent
 * with the new prior structure defined in this file, and need to be 
 * morphed into the correct structure.  I've kept them here to save the
 * numbers. EPN 01.31.05
 */

struct prior_s *default_single_prior(void) {
  struct prior_s *pri;
  int i, j;
#define numcomponents 8
  assert(numcomponents < MAXDCHLET);
  float q[numcomponents] = {
    0.085091850427, 0.015935406086, 0.102013232739, 0.415954530541,
    0.074470557341, 0.055442639402, 0.118379098369, 0.132712685095
  };
  float m[numcomponents][MAXABET] = {
    {0.575686380127, 0.756214632926, 0.340269621276, 13.774558068728, },
    {153.865583955384, 0.235000107300, 0.356622653787, 0.006812718667, },
    {176.440373997567, 0.935905951648, 1.292808081312, 1.617069444109, },
    {1.696250324914, 1.128033754503, 0.955462899400, 1.676465850057, },
    {0.074365531036, 0.039185613484, 0.063868972113, 0.042432587902, },
    {0.615068901818, 14.630712353118, 0.298404817403, 0.864718655041, },
    {1.163176461349, 0.408090165233, 11.188793743319, 0.699118301558, },
    {16.417200192194, 0.980503286582, 1.132071515554, 1.376129445524, },
  };
  pri = MallocOrDie(sizeof(struct prior_s));
  pri->mntnum = numcomponents;
  pri->mntasize = MAXABET;
  for (i = 0; i < numcomponents; i++) {
    pri->mntq[i] = q[i];
    for (j = 0; j < pri->mntasize; j++) {
      pri->mnt[i][j] = m[i][j];
    }
  }
#undef numcomponents
  return pri;
}

struct prior_s *default_basepair_prior(void) {
  struct prior_s *pri;
  int i, j;
#define numcomponents 9
  assert(numcomponents < MAXDCHLET);
  float q[numcomponents] = {
    0.030512242264, 0.070312169889, 0.118499696300, 0.181025557995,
    0.188791659665, 0.157630937531, 0.041708924031, 0.095930656547,
    0.115588155778    
  };
  float m[numcomponents][MAXABET*MAXABET] = {
    {0.571860339721, 0.605642194896, 0.548004739487, 1.570353271532, 0.591611867703, 0.469713257214, 1.447411319683, 0.600381079228, 0.520096937350, 1.867142019076, 0.470428282443, 1.165356324744, 1.528348208160, 0.686072963473, 1.072148274499, 0.659833749087},
    {0.116757286812, 0.052661180881, 0.067541712113, 0.258482314714, 0.152527972588, 0.034460232010, 0.416430364713, 0.051541326273, 0.079542103337, 0.162883420833, 0.042615616796, 0.123363759874, 0.922897266376, 0.078567729294, 0.315242459757, 0.116457644231},
    {0.028961414077, 0.022849036260, 0.120089637379, 0.509884713979, 0.142464495045, 0.079507804767, 21.835608089779, 0.070200164694, 0.005189494879, 0.540651647339, 0.117833357497, 0.128182594376, 1.766866842025, 0.016341625779, 0.832665494899, 0.058379188171},
    {0.000926960236, 0.008100076237, 0.001794303710, 0.114209483231, 0.001459159085, 0.000053878201, 0.072605927746, 0.005533021345, 0.003941720307, 0.095421675098, 0.004844990769, 0.072393572779, 0.099144450569, 0.002561491533, 0.043103588084, 0.008080970629},
    {0.002163861165, 0.007785521817, 0.003483930554, 0.625515668281, 0.018621932500, 0.001352139642, 1.371471086809, 0.007920737783, 0.000946403264, 0.688821384972, 0.002203762108, 0.192533693864, 0.979473608513, 0.000916007398, 0.347662973488, 0.020677924150},
    {0.083035113547, 0.166815168558, 0.042669979127, 3.415107328082, 0.023530116520, 0.047677945396, 1.183956650707, 0.059920099115, 0.076614058723, 5.434261851985, 0.095284240991, 0.889915882997, 1.201576769946, 0.074453244946, 0.397879304331, 0.130525904952},
    {0.217001113139, 0.388746098242, 0.134680826556, 24.923110155367, 0.102582693868, 0.131678864943, 1.150978162882, 0.256720461728, 0.150993730345, 3.200824712363, 0.077595421397, 1.025428618792, 1.228870901327, 0.143610901605, 0.406308970402, 0.322809888354},
    {0.129043208355, 0.112308496092, 0.116841517642, 2.878927926806, 0.306789207829, 0.078411064993, 6.377836578660, 0.114524370807, 0.094192610036, 2.566493997218, 0.096694574300, 0.791295335090, 6.907854285192, 0.132657156809, 1.225349985791, 0.296596767798},
    {0.005830777296, 0.153807106950, 0.003131256711, 1.340589241710, 0.006802639527, 0.135277067812, 0.487492640368, 0.009160116179, 0.068942867388, 29.409376576276, 0.099733235653, 0.722700985558, 0.500134122079, 0.124671165331, 0.105694456385, 0.025741311658},
  };
  pri = MallocOrDie(sizeof(struct prior_s));
  pri->mbpnum = numcomponents;
  pri->mbpasize = MAXABET*MAXABET;
  for (i = 0; i < pri->mbpnum; i++) {
    pri->mbpq[i] = q[i];
    for (j = 0; j < pri->mbpasize; j++) {
      pri->mbp[i][j] = m[i][j];
    }
  }
#undef numcomponents
  return pri;
}

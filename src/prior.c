/* prior.c
 * Dirichlet priors for parameterizing a new model. 
 *
 * Original code from Eric Nawrocki. Adapted by SRE.
 * SRE, Thu Apr  7 10:44:13 2005
 * SVN $Id$
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <easel.h>
#include <esl_vectorops.h>
#include <esl_fileparser.h>

#include "prior.h"


/* Function: Prior_Read()
 * 
 * Purpose:  Input a transition prior from an open stream
 *           (probably an open file).
 */
Prior_t *
Prior_Read(FILE *fp) 
{
  Prior_t         *pri;
  ESL_FILEPARSER  *efp;
  char            *tok;
  int              toklen;
  int              status;

  int              i;       /*counter over transition sets*/
  int              j;       /*counter over components in a mixture*/
  int              k;       /*counter over alphas*/
  int              curr_state_id; 
  int              curr_next_node_id; 
  int              a, b;

  pri = MallocOrDie (sizeof(Prior_t));
  pri->maxnq     = 0;
  pri->maxnalpha = 0;

  for(a = 0; a < UNIQUESTATES; a++)
    for(b = 0; b < NODETYPES; b++)
      pri->tsetmap[a][b] = -1;

  if ((efp = esl_fileparser_Create(fp)) == NULL)
    Die("Failed to associate open prior file stream with fileparser");
  esl_fileparser_SetCommentChar(efp, '#');

  /* First entry is the strategy: "Dirichlet" is the only possibility now. */
  if ((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) 
    Die("%s\nPrior file parse failed, on first (Dirichlet) field", efp->errbuf);
  if (strcasecmp(tok, "Dirichlet") != 0)
    Die("No such prior strategy %s\n", tok);
 
  /* Second entry is NTRANSSETS, which ought to be 74 */
  if ((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK)
    Die("%s\nPrior file parse failed reading NTRANSSETS", efp->errbuf);
  pri->tsetnum = atoi(tok);
  pri->t = MallocOrDie(sizeof(ESL_MIXDCHLET *) * pri->tsetnum);
  
  /* Transition section: a whole bunch of mixture Dirichlets.
   */
  for (i = 0; i < pri->tsetnum; i++)
    {
      if ((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK)
	Die("%s\nPrior file parse failed at line %d reading unique statetype", 
	    efp->errbuf, efp->linenumber);
      if ((curr_state_id = UniqueStateCode(tok)) == -1)
	Die("%s is not a uniq state;\nPrior file parse failed, line %d\n",
	    tok, efp->linenumber);

      if ((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK)
	Die("%s\nPrior file parse failed reading node code", efp->errbuf);
      if ((curr_next_node_id = NodeCode(tok)) == -1)
	Die("%s is not a node code;\nPrior file parse failed, line %d\n",
	    tok, efp->linenumber);

      pri->tsetmap[curr_state_id][curr_next_node_id] = i;

      if (esl_mixdchlet_Read(efp, &(pri->t[i])) != eslOK)
	Die("%s\nPrior file parse failed, reading transition prior %d at line %d.",
	    efp->errbuf, i, efp->linenumber);
      if (pri->t[i]->N > pri->maxnq)     pri->maxnq     = pri->t[i]->N;
      if (pri->t[i]->K > pri->maxnalpha) pri->maxnalpha = pri->t[i]->K;
    }
  
  /* Consensus base pair emission prior section.
   */
  if (esl_mixdchlet_Read(efp, &(pri->mbp)) != eslOK) 
    Die("%s\nPrior file parse failed in base pair priors at line %d\n", 
	efp->errbuf, efp->linenumber);
  if (pri->mbp->N > pri->maxnq)     pri->maxnq     = pri->mbp->N;
  if (pri->mbp->K > pri->maxnalpha) pri->maxnalpha = pri->mbp->K;

  /* Consensus singlet emission prior section.
   */
  if (esl_mixdchlet_Read(efp, &(pri->mnt)) != eslOK) 
    Die("%s\nPrior file parse failed in consensus singlet priors at line %d\n", 
	efp->errbuf, efp->linenumber);
  if (pri->mnt->N > pri->maxnq)     pri->maxnq     = pri->mnt->N;
  if (pri->mnt->K > pri->maxnalpha) pri->maxnalpha = pri->mnt->K;

  /* Nonconsensus singlet emission prior section.
   */
  if (esl_mixdchlet_Read(efp, &(pri->i)) != eslOK)  
    Die("%s\nPrior file parse failed in nonconsensus singlet priors at line %d\n", 
	efp->errbuf, efp->linenumber);
  if (pri->i->N > pri->maxnq)     pri->maxnq     = pri->i->N;
  if (pri->i->K > pri->maxnalpha) pri->maxnalpha = pri->i->K;

  esl_fileparser_Destroy(efp);
  return pri;
}



/* Function:  Prior_Destroy()
 * Incept:    SRE, Mon Apr 11 10:06:34 2005 [St. Louis]
 *
 * Purpose:   Free's a prior.
 */
void
Prior_Destroy(Prior_t *pri)
{
  int i;
  if (pri == NULL) return;
  if (pri->t != NULL) 
    {
      for (i = 0; i < pri->tsetnum; i++)
	esl_mixdchlet_Destroy(pri->t[i]);
      free(pri->t);
    }
  esl_mixdchlet_Destroy(pri->mbp);
  esl_mixdchlet_Destroy(pri->mnt);
  esl_mixdchlet_Destroy(pri->i);
  free(pri);
}
      
/* Function: PriorifyCM()
 * 
 * Purpose:  Given a CM containing counts; add pseudocounts to a CM 
 *           using Dirichlet priors, and renormalize the CM.
 *           
 *           The Easel Dirichlet routines are in double-precision, 
 *           whereas the CM is in floats, so we have some internal vector
 *           conversion going on here.
 * 
 * Args:     CM -- the CM to add counts to (counts form)
 *           pri -- the Dirichlet prior to use
 *           
 * Return:   (void)
 *           CM is changed from counts to probability form.
 */          
void
PriorifyCM(CM_t *cm, Prior_t *pri)
{
  int v;		/* counter for model position   */
  int setnum;           /* number of set to use */
  int nxtndtype;        /* type of next node */
  double *counts;	/* double copy of floating-pt counts in the CM */
  double *probs;	/* double copy of new probability parameters */
  double *mixq;		/* posterior probs of mixture components, P(q | c) */
  int     i;

  /* Create our temporary buffers; counts, probs, and mixq.
   */
  counts = MallocOrDie(sizeof(double) * pri->maxnalpha);
  probs  = MallocOrDie(sizeof(double) * pri->maxnalpha);
  mixq   = MallocOrDie(sizeof(double) * pri->maxnq);
                                 
  for (v = 0; v < cm->M; v++)
    {
      /* Priorify transition vector if not a BIF or E state */
      if (cm->sttype[v] != B_st && cm->sttype[v] != E_st)
	{
	  /* Determine which transition set to use. 
	   * Current unique state id is easy (cm->stid[v]). 
           * Type of next node is a little trickier. The trick is 
           * to use the ndidx of the *last* state this state v
           * connects to. This is guaranteed to be in the next node,
           * cannot be an insert state of the current node.
	   */
	  nxtndtype = cm->ndtype[cm->ndidx[cm->cfirst[v] + cm->cnum[v] - 1]];
	  setnum = pri->tsetmap[cm->stid[v]][nxtndtype];
	  for (i = 0; i < cm->cnum[v]; i++)
	    counts[i] = (double) cm->t[v][i];

	  esl_mixdchlet_MPParameters(counts, cm->cnum[v], 
				     pri->t[setnum],
				     mixq, probs);

	  for (i = 0; i < cm->cnum[v]; i++)
	    cm->t[v][i] = (float) probs[i];
	}
      

      /* Emission priors
       */
      if (cm->sttype[v] == MP_st)
	{       /* Consensus base pairs */
	  for (i = 0; i < MAXABET*MAXABET; i++)
	    counts[i] = (double) cm->e[v][i];

	  esl_mixdchlet_MPParameters(counts, MAXABET*MAXABET,
				     pri->mbp,
				     mixq, probs);

	  for (i = 0; i < MAXABET*MAXABET; i++)
	    cm->e[v][i] = (float) probs[i];
	}
      else if (cm->stid[v] == MATL_ML || cm->stid[v] == MATR_MR)
	{      /* Consensus singlets */
	  for (i = 0; i < MAXABET; i++)
	    counts[i] = (double) cm->e[v][i];

	  esl_mixdchlet_MPParameters(counts, MAXABET,
				     pri->mnt,
				     mixq, probs);

	  for (i = 0; i < MAXABET; i++)
	    cm->e[v][i] = (float) probs[i];
	}
      else if (cm->sttype[v] == IL_st || cm->sttype[v] == IR_st ||
	       cm->stid[v] == MATP_ML || cm->stid[v] == MATP_MR)
	{	/* nonconsensus singlets */
	  for (i = 0; i < MAXABET; i++)
	    counts[i] = (double) cm->e[v][i];

	  esl_mixdchlet_MPParameters(counts, MAXABET,
				     pri->i,
				     mixq, probs);

	  for (i = 0; i < MAXABET; i++)
	    cm->e[v][i] = (float) probs[i];
	}
    }/* end loop over states v */

  free(mixq);
  free(counts);
  free(probs);
}



/* Below are Alex's default prior structures, but they are not consistent
 * with the new prior structure defined in this file, and need to be 
 * morphed into the correct structure.  I've kept them here to save the
 * numbers. EPN 01.31.05
 */
/*
struct prior_s *default_single_prior(void) {
  struct prior_s *pri;
  int i, j;
#define numcomponents 8
  assert(numcomponents < MAXDCHLET);
  double q[numcomponents] = {
    0.085091850427, 0.015935406086, 0.102013232739, 0.415954530541,
    0.074470557341, 0.055442639402, 0.118379098369, 0.132712685095
  };
  double m[numcomponents][MAXABET] = {
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
  double q[numcomponents] = {
    0.030512242264, 0.070312169889, 0.118499696300, 0.181025557995,
    0.188791659665, 0.157630937531, 0.041708924031, 0.095930656547,
    0.115588155778    
  };
  double m[numcomponents][MAXABET*MAXABET] = {
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
*/

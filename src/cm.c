/* cm.c
 * SRE, Sat Jul 29 09:01:20 2000 [St. Louis]
 * SVN $Id$
 * 
 * Routines for dealing with the CM data structure.
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************  
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "easel.h"
#include "esl_vectorops.h"
#include "esl_alphabet.h"
#include "esl_stack.h"

#include "funcs.h"
#include "structs.h"

/* Function: CreateCM(); CreateCMShell(); CreateCMBody()
 * Date:     SRE, Sat Jul 29 09:02:16 2000 [St. Louis]
 *
 * Purpose:  Create a covariance model, given the number of states 
 *           and nodes that should be in it.
 *
 *           Allocation is usually one step: CreateCM(N, M).
 *           
 *           In cmio.c, allocation is two-step: CreateCMShell()
 *           then CreateCMBody(cm, N, M). This way we can create
 *           a model, start storing header info in it (including M
 *           and N themselves), and only after we read M and N
 *           do we allocate the bulk of the model.
 *
 * Args:     nnodes  =  number of nodes in the model
 *           nstates = number of states in the model
 *
 * Returns:  ptr to allocated cm. 
 *           Caller is responsible for free'ing the cm.
 */
CM_t *
CreateCM(int nnodes, int nstates, const ESL_ALPHABET *abc)
{
  CM_t *cm;

  cm = CreateCMShell();
  if (cm != NULL) CreateCMBody(cm, nnodes, nstates, abc);
  return cm;
}
CM_t *
CreateCMShell(void)
{
  int status;
  CM_t *cm;

  ESL_ALLOC(cm, sizeof(CM_t));
				/* general information: added later */
  cm->abc    = NULL;

  cm->name   = NULL;
  cm->acc    = NULL;
  cm->desc   = NULL;
  cm->annote = NULL;

				/* null model information */
  cm->null   = NULL;

				/* structural information */
  cm->M      = 0;
  cm->clen   = 0;
  cm->sttype = NULL;
  cm->ndidx  = NULL;
  cm->stid   = NULL;
  cm->cfirst = NULL;
  cm->cnum   = NULL;
  cm->plast  = NULL;
  cm->pnum   = NULL;
				/* node->state map information */
  cm->nodes  = 0;
  cm->nodemap= NULL;
  cm->ndtype = NULL;
				/* parameter information */
  cm->t      = NULL;
  cm->e      = NULL;
  cm->begin  = NULL;
  cm->end    = NULL;
  cm->tsc    = NULL;
  cm->esc    = NULL;
  cm->beginsc= NULL;
  cm->endsc  = NULL;

  cm->flags         = 0;

  cm->W      = 200;           /* for backwards compatibility */
  cm->el_selfsc = 0.;         /* this is backwards compatible also */
  
  cm->dmin   = NULL;
  cm->dmax   = NULL;
  cm->beta   = DEFAULT_BETA;     /* 1E-7 the default beta (tail loss for QDB) */
  cm->tau    = DEFAULT_TAU;      /* 1E-7 the default tau  (tail loss for HMM banding) */
  cm->cp9    = NULL;          
  cm->cp9map = NULL;
  cm->enf_start   = 0;
  cm->enf_seq     = NULL;
  cm->enf_scdiff  = 0.;
  cm->sc_boost    = 0.;
  cm->ffract      = 0.;
  cm->cutoff_type = DEFAULT_CM_CUTOFF_TYPE; /* score cutoff */
  cm->cutoff      = DEFAULT_CM_CUTOFF;      /* 0.0 bits */
  cm->cp9_cutoff_type = DEFAULT_CP9_CUTOFF_TYPE; /* score cutoff */
  cm->cp9_cutoff      = DEFAULT_CP9_CUTOFF;      /* 0.0 bits */
  cm->cp9_sc_boost = 0.;
  cm->root_trans   = NULL;
  cm->hmmpad       = DEFAULT_HMMPAD; /* 0 residues */
  cm->stats        = NULL;
  cm->pbegin       = DEFAULT_PBEGIN; /* summed probability of internal local begin */
  cm->pend         = DEFAULT_PEND;   /* summed probability of internal local end */

  cm->ga     = 0.;  /* only valid if cm->flags & CMH_GA */
  cm->tc     = 0.;  /* only valid if cm->flags & CMH_TC */
  cm->nc     = 0.;  /* only valid if cm->flags & CMH_NC */
  cm->eff_nseq = 0.;  
  cm->nseq = 0;
  cm->clen = 0;
  return cm;

 ERROR:
  esl_fatal("Memory allocation error.\n");
  return NULL; /* never reached */
}
void
CreateCMBody(CM_t *cm, int nnodes, int nstates, const ESL_ALPHABET *abc)
{
  int status;
  int v;
                                /* alphabet, only a reference */
  cm->abc    = abc; 
				/* structural information */
  cm->M      = nstates;

				/* null model information */
  CMAllocNullModel(cm);         

  ESL_ALLOC(cm->sttype, (nstates+1) * sizeof(char));
  ESL_ALLOC(cm->ndidx,   nstates    * sizeof(int));
  ESL_ALLOC(cm->stid,   (nstates+1) * sizeof(char));
  ESL_ALLOC(cm->cfirst,  nstates    * sizeof(int));
  ESL_ALLOC(cm->cnum,    nstates    * sizeof(int));
  ESL_ALLOC(cm->plast,   nstates    * sizeof(int));
  ESL_ALLOC(cm->pnum,    nstates    * sizeof(int));
				/* node->state map information */
  cm->nodes  = nnodes;
  ESL_ALLOC(cm->nodemap, nnodes  * sizeof(int));
  ESL_ALLOC(cm->ndtype,  nnodes  * sizeof(char));
  
  /* parameter information */
  /* level 1 */
  ESL_ALLOC(cm->t,    (nstates) * sizeof(float *));
  ESL_ALLOC(cm->e,    (nstates) * sizeof(float *));
  ESL_ALLOC(cm->tsc,  (nstates) * sizeof(float *));
  ESL_ALLOC(cm->esc,  (nstates) * sizeof(float *));
  ESL_ALLOC(cm->itsc, (nstates) * sizeof(int *));
  ESL_ALLOC(cm->iesc, (nstates) * sizeof(int *));
  cm->t[0]   = NULL;
  cm->e[0]   = NULL;
  cm->tsc[0] = NULL;
  cm->esc[0] = NULL;
  cm->itsc[0]= NULL;
  cm->iesc[0]= NULL;
  ESL_ALLOC(cm->begin,   (nstates) * sizeof(float));
  ESL_ALLOC(cm->end,     (nstates) * sizeof(float));
  ESL_ALLOC(cm->beginsc, (nstates) * sizeof(float));
  ESL_ALLOC(cm->endsc,   (nstates) * sizeof(float));
  ESL_ALLOC(cm->ibeginsc,(nstates) * sizeof(int));
  ESL_ALLOC(cm->iendsc,  (nstates) * sizeof(int));
  
  /* level 2 */
  ESL_ALLOC(cm->t[0],    MAXCONNECT * nstates * sizeof(float));
  ESL_ALLOC(cm->e[0],    cm->abc->K * cm->abc->K * nstates * sizeof(float));
  ESL_ALLOC(cm->tsc[0],  MAXCONNECT * nstates * sizeof(float));
  ESL_ALLOC(cm->esc[0],  cm->abc->K * cm->abc->K * nstates * sizeof(float));
  ESL_ALLOC(cm->itsc[0], MAXCONNECT * nstates * sizeof(int));
  ESL_ALLOC(cm->iesc[0], cm->abc->K * cm->abc->K * nstates * sizeof(int));
  for (v = 0; v < nstates; v++) 
    {
      cm->e[v]    = cm->e[0]    + v * (cm->abc->K * cm->abc->K);
      cm->t[v]    = cm->t[0]    + v * MAXCONNECT;
      cm->esc[v]  = cm->esc[0]  + v * (cm->abc->K * cm->abc->K);
      cm->tsc[v]  = cm->tsc[0]  + v * MAXCONNECT;
      cm->iesc[v] = cm->iesc[0] + v * (cm->abc->K * cm->abc->K);
      cm->itsc[v] = cm->itsc[0] + v * MAXCONNECT;
    }

  /* the EL state at M is special: we only need state
   * type info recorded, so functions looking at parsetrees  
   * can interpret what an "M" index means.
   */
  cm->sttype[cm->M] = EL_st;
  cm->stid[cm->M]   = END_EL;

  cm->flags         = 0;
  cm->config_opts   = 0;
  cm->align_opts    = 0;
  cm->search_opts   = 0;
  cm->dmin          = NULL;
  cm->dmax          = NULL;
  cm->cp9           = NULL;
  cm->cp9map        = NULL;
  /* we'll allocate the cp9 and cp9map only if nec inside ConfigCM() */
  return;

 ERROR:
  return; /* never reached */
}


/* Function: CMZero()
 * Date:     SRE, Mon Jul 31 19:14:31 2000 [St. Louis]
 *
 * Purpose:  Initialize the probability parameters and scores of a CM to zero.
 *
 * Returns:  (void)
 */
void 
CMZero(CM_t *cm)
{
  int v;			/* counter over states                 */

  for (v = 0; v < cm->M; v++) {
    esl_vec_FSet(cm->e[v],    (cm->abc->K * cm->abc->K), 0.);
    esl_vec_FSet(cm->t[v],    MAXCONNECT,                0.);
    esl_vec_FSet(cm->esc[v],  (cm->abc->K * cm->abc->K), 0.);
    esl_vec_FSet(cm->tsc[v],  MAXCONNECT,                0.);
    esl_vec_ISet(cm->iesc[v], (cm->abc->K * cm->abc->K), 0);
    esl_vec_ISet(cm->itsc[v], MAXCONNECT,                0);
  }
  esl_vec_FSet(cm->begin,    cm->M, 0.);
  esl_vec_FSet(cm->end,      cm->M, 0.);
  esl_vec_FSet(cm->beginsc,  cm->M, 0.);
  esl_vec_FSet(cm->endsc,    cm->M, 0.);
  esl_vec_ISet(cm->ibeginsc, cm->M, 0);
  esl_vec_ISet(cm->iendsc,   cm->M, 0);
}

/* Function:  CMRenormalize()
 * Incept:    SRE, Wed Aug 14 14:16:55 2002 [St. Louis]
 *
 * Purpose:   Renormalize all the probability distributions in a CM.
 *            Used by cmio.c's flatfile parser, for example.
 *
 * Xref:      STL6 p.108
 */
void
CMRenormalize(CM_t *cm)
{
  int v;

  esl_vec_FNorm(cm->null, cm->abc->K);
  for (v = 0; v < cm->M; v++)
    {
      if (cm->cnum[v] > 0 && cm->sttype[v] != B_st)
	esl_vec_FNorm(cm->t[v], cm->cnum[v]);
      
      if (cm->sttype[v] == ML_st || cm->sttype[v] == MR_st || cm->sttype[v] == IL_st || cm->sttype[v] == IR_st)
	esl_vec_FNorm(cm->e[v], cm->abc->K);
      if (cm->sttype[v] == MP_st)
	esl_vec_FNorm(cm->e[v], cm->abc->K * cm->abc->K);
    }
  if (cm->flags & CM_LOCAL_BEGIN) esl_vec_FNorm(cm->begin, cm->M);
  if (cm->flags & CM_LOCAL_END)   esl_fatal("Renormalization of models in local end mode not supported yet");
}


/* Function: FreeCM()
 * Date:     SRE, Sat Jul 29 11:22:32 2000 [St. Louis]
 *
 * Purpose:  Free a CM data structure.
 *
 * Args:     cm - the model to free. (duh).
 *
 * Note:      Remember, leave reference pointer to abc alone.
 *            This is under the application's control not ours.
 *
 * Returns:  (void)
 */
void
FreeCM(CM_t *cm)
{
  if (cm->name   != NULL) free(cm->name);
  if (cm->acc    != NULL) free(cm->acc);
  if (cm->desc   != NULL) free(cm->desc);
  if (cm->annote != NULL) free(cm->annote);
  /* This gives a memory error, b/c it's copied from optarg
   * in cmsearch.c (I can't find where optarg is alloc'ed,
   * worst case scenario: small memory leak here */
  /*if (cm->enf_seq != NULL) free(cm->enf_seq);*/

  free(cm->null);
  free(cm->sttype);
  free(cm->ndidx);
  free(cm->stid);
  free(cm->cfirst);
  free(cm->cnum);
  free(cm->plast);
  free(cm->pnum);
  free(cm->nodemap);
  free(cm->ndtype);

  free(cm->t[0]);
  free(cm->t);
  free(cm->e[0]);
  free(cm->e);
  free(cm->tsc[0]);
  free(cm->tsc);
  free(cm->esc[0]);
  free(cm->esc);
  free(cm->itsc[0]);
  free(cm->itsc);
  free(cm->iesc[0]);
  free(cm->iesc);
  free(cm->begin);
  free(cm->end);
  free(cm->beginsc);
  free(cm->endsc);
  free(cm->ibeginsc);
  free(cm->iendsc);
  free(cm->dmin);
  free(cm->dmax);
  if(cm->cp9map     != NULL) FreeCP9Map(cm->cp9map);
  if(cm->cp9        != NULL) FreeCPlan9(cm->cp9);
  if(cm->root_trans != NULL) free(cm->root_trans);
  if(cm->stats      != NULL) FreeCMStats(cm->stats);
  free(cm);
}

/* Function: DefaultNullModel()
 * Date:     SRE, Tue Aug  1 15:31:52 2000 [St. Louis]
 *
 * Purpose:  Allocate and initialize a float vector
 *           that will be a template null model to 
 *           equiprobable (e.g. 0.25)
 */
int
DefaultNullModel(const ESL_ALPHABET *abc, float **ret_null)
{
  /* Contract check */
  if(abc      == NULL) esl_fatal("ERROR in CMCreateNullModel, cm->abc is NULL.\n");

  int status;
  float *null = NULL;
  ESL_ALLOC(null, sizeof(float) * abc->K);
  int x;
  for (x = 0; x < abc->K; x++)
    null[x] = 1./(float) abc->K;
  *ret_null = null;
  return eslOK;

 ERROR:
  if(null != NULL) free(null);
  return status;
}

/* Function: CMAllocNullModel()
 *
 * Purpose:  Allocate the null model section of a CM
 *           and fill it with default, equiprobable 
 *           null distro.
 */
int
CMAllocNullModel(CM_t *cm)
{
  int status;

  /* Contract check */
  if(cm->abc  == NULL) esl_fatal("ERROR in CMAllocNullModel, cm->abc is NULL.\n");
  if(cm->null != NULL) esl_fatal("ERROR in CMAllocNullModel, cm->null is not NULL.\n");

  status = DefaultNullModel(cm->abc, &(cm->null));
  return status;
}

/* Function: CMSetNullModel()
 *
 * Purpose:  Set the null model section of a CM.
 */
void
CMSetNullModel(CM_t *cm, float *null)
{
  /* Contract check */
  if(cm->abc  == NULL) esl_fatal("ERROR in CMCreateNullModel, cm->abc is NULL.\n");
  if(cm->null == NULL) esl_fatal("ERROR in CMSetNullModel, cm->null is NULL.\n");

  int x;
  for (x = 0; x < cm->abc->K; x++)
    cm->null[x] = null[x];
}


/* Function: CMReadNullModel()
 * EPN 10.19.05
 * based on SRE's HMMER's cm.c's P7ReadNullModel() 
 *
 * Purpose:  Read a CM null model from a file.
 *           ret_null is filled with a newly allocated
 *           float vector that is the null model.
 *
 * Returns:  eslOK on success.
 */
int
CMReadNullModel(const ESL_ALPHABET *abc, char *nullfile, float **ret_null)
{
  /* Contract check */
  if(abc  == NULL) esl_fatal("ERROR in CMReadNullModel, abc is NULL.\n");

  int status;
  float *null = NULL;
  FILE *fp;
  char *buf;
  char *s;
  int   n;			/* length of buf */
  int   x;
  char *tok;
  int   toklen;
  float sum;

  ESL_ALLOC(null, sizeof(float) * abc->K);
  buf = NULL;
  n   = 0;
  sum = 0.;
  /* Expects a file with cm->abc->K lines that don't begin with "# ".
   * The first token of each of these 4 lines is read as 
   * the background probability of A, C, G, and U (in that order)
   * Then does a check to make sure the 4 read in values
   * sum to 1.0 exactly.
   */
  if ((fp = fopen(nullfile, "r")) == NULL)
    esl_fatal("Failed to open null model file %s\n", nullfile);
  
  /* parse the file */
  x = 0;
  while(x < abc->K) {
    if((status = esl_fgets(&buf, &n, fp)) != eslOK) goto ERROR;
    s   = buf;
    if((status = esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto ERROR;
    if(strcmp(tok, "#") != 0)
      {      
	null[x] = atof(tok);
	sum += null[x];
	x++;
      }
  }
  /*fragile*/
  if(sum > 1.00001 || sum < 0.99999)
    esl_fatal("%s is not in CM null model file format.\nThere are not %d background probabilities that sum to exactly 1.0", nullfile, abc->K);
  esl_vec_FNorm(null, abc->K);
    
  *ret_null = null;
  fclose(fp);
  return eslOK;

 ERROR:
  fclose(fp);
  if(null != NULL) free(null);
  return status;
}

/* Function: CMSimpleProbify()
 * Date:     SRE, Tue Aug  1 11:07:17 2000 [St. Louis]
 *
 * Purpose:  Convert a counts-based CM to probability form, using
 *           a plus-one Laplace prior.
 */
void
CMSimpleProbify(CM_t *cm)
{
  int v,x;

  for (v = 0; v < cm->M; v++) 
    {
      /* Transitions. B, E have no transition probabilities.
       */
      if (cm->sttype[v] != B_st && cm->sttype[v] != E_st) 
	{
	  for (x = 0; x < cm->cnum[v]; x++) cm->t[v][x] += 1.0; /* Laplace prior */
	  esl_vec_FNorm(cm->t[v], cm->cnum[v]);	                        /* normalize to a probability */
	}

      /* Emissions.
       */
      if (cm->sttype[v] == MP_st) 
	{
	  for (x = 0; x < cm->abc->K*cm->abc->K; x++) cm->e[v][x] += 1.0;
	  esl_vec_FNorm(cm->e[v], cm->abc->K*cm->abc->K);
	}
      else if (cm->sttype[v] == ML_st || cm->sttype[v] == MR_st || 
	       cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) 
	{
	  for (x = 0; x < cm->abc->K; x++) cm->e[v][x] += 1.0;
	  esl_vec_FNorm(cm->e[v], cm->abc->K);
	}
    }
}

/* Function: rsearch_CMProbifyEmissions()
 * Date:     EPN, Wed Mar 14 06:14:51 2007
 *
 * Purpose:  Convert emissions in a counts-based CM built from a single sequence, 
 *           we expect 1 count in each vector to be 1.0, all others 0.0,
 *           to probability form, using a RIBOSUM matrix with background
 *           and target frequencies.
 * 
 *           The code that does this in RSEARCH is buildcm.c::SingleSequenceLogoddsify(), 
 *           but that's different in that it fills in log odds scores, here we
 *           fill in probabilities, derived from the RIBOSUM log odds scores.
 * 
 * Returns:   <eslOK> on success.           
 *
 * Throws:    <eslEINVAL> if an emission vector does not have exactly 1 non-zero
 *                        count that is exactly 1.0 (or within 0.000001 of it)
 */
int
rsearch_CMProbifyEmissions(CM_t *cm, fullmat_t *fullmat)
{
  int v,x,y;
  int cur_emission;
  float thresh;
  int found_ct_flag;
  thresh = 0.000001;



  /* Check the contract. */
  if(fullmat->scores_flag) ESL_EXCEPTION(eslEINVAL, "in rsearch_CMProbifyEmissions(), matrix is in log odds mode, it should be in probs mode");
  if(!(cm->flags & CM_RSEARCHEMIT)) ESL_EXCEPTION(eslEINVAL, "in rsearch_CMProbifyEmissions(), CM_RSEARCHEMIT flag is down");
  
  for (v = 0; v < cm->M; v++) 
    {
      found_ct_flag = FALSE;
      if (cm->stid[v] == MATP_MP) 
	{
	  /* First, figure out which letter was in the query */
	  
	  for (x=0; x<cm->abc->K; x++) 
	    for (y=0; y<cm->abc->K; y++) 
	      if (fabs(cm->e[v][x*cm->abc->K+y] - 0.) > thresh) 
		{
		  if(found_ct_flag)
		    {
		      for (x=0; x<cm->abc->K; x++) 
			for (y=0; y<cm->abc->K; y++) 
			  printf("cm->e[v:%d][%d]: %f\n", v, (x*cm->abc->K+y), cm->e[v][(x*cm->abc->K+y)]);
		      ESL_EXCEPTION(eslEINVAL, "cm->e[v:%d] a MATP_MP has > 1 non-zero count"); 
		    }
		  cur_emission = numbered_basepair(cm->abc->sym[x], cm->abc->sym[y]);
		  found_ct_flag = TRUE;
		}
	  /* Now, set emission probs as target probs in correct cells of score matrix */
	  for (x=0; x<cm->abc->K*cm->abc->K; x++) 
	    cm->e[v][x] = fullmat->paired->matrix[matrix_index(cur_emission, x)];
	  esl_vec_FNorm(cm->e[v], cm->abc->K*cm->abc->K);
	}
      else if (cm->stid[v] == MATL_ML || cm->stid[v] == MATR_MR)
	{
	  for (x=0; x<cm->abc->K; x++) 
	    if (fabs(cm->e[v][x] - 0.) > thresh) 
		{
		  if(found_ct_flag) ESL_EXCEPTION(eslEINVAL, "cm->e[v:%d] a MAT{L,R}_M{L,R} has > 1 non-zero count"); 
		  cur_emission = numbered_nucleotide(cm->abc->sym[x]);
		  found_ct_flag = TRUE;
		}
	  /* Now, set emission probs as target probs in correct cells of score matrix */
	  for (x=0; x<cm->abc->K; x++) 
	    cm->e[v][x] = fullmat->unpaired->matrix[matrix_index(cur_emission, x)];
	  esl_vec_FNorm(cm->e[v], cm->abc->K);
	}
      else if (cm->stid[v] == MATP_ML || cm->stid[v] == MATP_MR)
	{
	  /* RSEARCH technique: determine residue emitted to left and right, 
	   * use target freqs from unpaired matrix for this residue. 
	   * Alternative technique: determine residue emitted to left and right, 
	   * marginalize target freqs from paired matrix for this residue. 
	   * RSEARCH technique currently implemented */
	  for (x=0; x<cm->abc->K; x++) 
	    for (y=0; y<cm->abc->K; y++) 
	      {
		if (cm->stid[v] == MATP_ML && (fabs(cm->e[(v-1)][x*cm->abc->K + y] - 0.) > thresh))
		  cur_emission = numbered_nucleotide(cm->abc->sym[x]);
		else if (cm->stid[v] == MATP_MR && (fabs(cm->e[(v-2)][x*cm->abc->K + y] - 0.) > thresh))
		  cur_emission = numbered_nucleotide(cm->abc->sym[y]);
		/* We don't have to check we have only 1 non-zero count, we've already
		 * done so when we filled e for the MATP_MP in the same node as v */
	      }		
	  /* fill emission probs */
	  for (x=0; x<cm->abc->K; x++) 
	    cm->e[v][x] = fullmat->unpaired->matrix[matrix_index(cur_emission, x)];
	  esl_vec_FNorm(cm->e[v], cm->abc->K);
	}
      else if (cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) 
	{
	  /* Don't give any score for emissions matching to an Insert state,
	   * but make sure we don't have any counts in any of these guys */
	  for (x = 0; x < cm->abc->K; x++) 
	    if(fabs(cm->e[v][x] - 0.) > thresh) ESL_EXCEPTION(eslEINVAL, "cm->e[v:%d] an I{L,R} has > 0 non-zero count", v); 
	  esl_vec_FNorm(cm->e[v], cm->abc->K); /* these will have all been zero */
	}
    }
  return eslOK;
}

/* Function: CMLogoddsify()
 * Date:     SRE, Tue Aug  1 15:18:26 2000 [St. Louis]
 *
 * Purpose:  Convert the probabilities in a CM to log-odds 
 *           EPN 12.19.06: also fill in integer log-odds scores.
 */
void
CMLogoddsify(CM_t *cm)
{
  int v, x, y;

  for (v = 0; v < cm->M; v++)
    {
      if (cm->sttype[v] != B_st && cm->sttype[v] != E_st)
	for (x = 0; x < cm->cnum[v]; x++)
	  {
	    cm->tsc[v][x]  = sreLOG2(cm->t[v][x]);
	    cm->itsc[v][x] = Prob2Score(cm->t[v][x], 1.0);
	    /*printf("cm->t[%4d][%2d]: %f itsc->e: %f itsc: %d\n", v, x, cm->t[v][x], Score2Prob(cm->itsc[v][x], 1.0), cm->itsc[v][x]);*/
	  }	    
      if (cm->sttype[v] == MP_st)
	for (x = 0; x < cm->abc->K; x++)
	  for (y = 0; y < cm->abc->K; y++)
	    {
	      cm->esc[v][x*cm->abc->K+y]  = sreLOG2(cm->e[v][x*cm->abc->K+y] / (cm->null[x]*cm->null[y]));
	      cm->iesc[v][x*cm->abc->K+y] = Prob2Score(cm->e[v][x*cm->abc->K+y], (cm->null[x]*cm->null[y]));
	      /*printf("cm->e[%4d][%2d]: %f iesc->e: %f iesc: %d\n", v, (x*cm->abc->K+y), cm->e[v][(x*cm->abc->K+y)], Score2Prob(cm->iesc[v][x*cm->abc->K+y], (cm->null[x]*cm->null[y])), cm->iesc[v][(x*cm->abc->K+y)]);*/
	    }
      if (cm->sttype[v] == ML_st || cm->sttype[v] == MR_st ||
	  cm->sttype[v] == IL_st || cm->sttype[v] == IR_st)
	for (x = 0; x < cm->abc->K; x++)
	  {
	    cm->esc[v][x]  = sreLOG2(cm->e[v][x] / cm->null[x]);
	    cm->iesc[v][x] = Prob2Score(cm->e[v][x], cm->null[x]);
	    /*printf("cm->e[%4d][%2d]: %f esc: %f null[%d]: %f\n", v, x, cm->e[v][x], cm->esc[v][x], x, cm->null[x]);*/
	    /*printf("cm->e[%4d][%2d]: %f iesc->e: %f iesc: %d\n", v, x, cm->e[v][x], Score2Prob(cm->iesc[v][x], (cm->null[x])), cm->iesc[v][x]);*/
	  }
      /* These work even if begin/end distributions are inactive 0's,
       * sreLOG2 will set beginsc, endsc to -infinity.
       */
      cm->beginsc[v]  = sreLOG2(cm->begin[v]);
      cm->ibeginsc[v] = Prob2Score(cm->begin[v], 1.0);
      /*printf("cm->begin[%4d]: %f ibeginsc->e: %f ibeginsc: %d\n", v, cm->begin[v], Score2Prob(cm->ibeginsc[v], 1.0), cm->ibeginsc[v]);*/

      cm->endsc[v]    = sreLOG2(cm->end[v]);
      cm->iendsc[v]   = Prob2Score(cm->end[v], 1.0);
      /*printf("cm->end[%4d]: %f iendsc->e: %f iendsc: %d\n\n", v, cm->end[v], Score2Prob(cm->iendsc[v], 1.0), cm->iendsc[v]);*/
    }

  cm->iel_selfsc = Prob2Score(sreEXP2(cm->el_selfsc), 1.0);
  /*printf("cm->el_selfsc: %f prob: %f cm->iel_selfsc: %d prob: %f\n", cm->el_selfsc, 
	 (sreEXP2(cm->el_selfsc)), cm->iel_selfsc, (Score2Prob(cm->iel_selfsc, 1.0)));
	 printf("-INFTY: %d prob: %f 2^: %f\n", -INFTY, (Score2Prob(-INFTY, 1.0)), sreEXP2(-INFTY));*/

  /* Potentially, overwrite transitions with non-probabilistic 
   * RSEARCH transitions. Currently only default transition
   * parameters are allowed, these are defined as DEFAULT_R*
   * in structs.h */
  if(cm->flags & CM_RSEARCHTRANS)
    {
      float           alpha =   DEFAULT_RALPHA; 
      float           beta =    DEFAULT_RBETA;
      float           alphap =  DEFAULT_RALPHAP;
      float           betap =   DEFAULT_RBETAP;
      float           beginsc = DEFAULT_RBEGINSC;
      float           endsc =   DEFAULT_RENDSC;
      int             nd;
      /* First do the normal transitions */
      for (v=0; v<cm->M; v++) 
	{
	  if (cm->sttype[v] != B_st && cm->sttype[v] != E_st) 
	    {
	      for (x=0; x<cm->cnum[v]; x++) 
		{
		  cm->tsc[v][x] = -1. * rsearch_calculate_gap_penalty 
		    (cm->stid[v], cm->stid[cm->cfirst[v]+x], 
		     cm->ndtype[cm->ndidx[v]], cm->ndtype[cm->ndidx[cm->cfirst[v]+x]],
		     alpha, beta, alphap, betap);
		  /* alphas and rbetas were positive -- gap score is a penalty, so
		     multiply by -1 */
		  cm->itsc[v][x] = INTSCALE * cm->tsc[v][x];
		}
	    }
	}
      /* Overwrite local begin and end scores */
      for (v=cm->M - 1; v>=0; v--) {
	cm->beginsc[v] = IMPOSSIBLE;
	cm->endsc[v] = IMPOSSIBLE;
      }
      
      /* beginsc states */
      for (nd = 2; nd < cm->nodes; nd++) {
	if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
	    cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BIF_nd)
	  {	 
	    cm->beginsc[cm->nodemap[nd]] = beginsc;
	    cm->ibeginsc[cm->nodemap[nd]] = INTSCALE * beginsc;
	  }
      }
      
      /* endsc states */
      for (nd = 1; nd < cm->nodes; nd++) {
	if ((cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
	     cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BEGL_nd ||
	     cm->ndtype[nd] == BEGR_nd) &&
	    cm->ndtype[nd+1] != END_nd)
	  {
	  cm->endsc[cm->nodemap[nd]] = endsc;
	  cm->iendsc[cm->nodemap[nd]] = INTSCALE * endsc;
	  }
      }
      
      cm->flags |= CM_LOCAL_BEGIN;
      cm->flags |= CM_LOCAL_END;
    }
  /* raise flag saying we have valid log odds scores */
  cm->flags |= CMH_BITS;
}

/* Function:  CMHackInsertScores()
 * Incept:    SRE, Wed Jul 24 09:48:22 2002 [St. Louis]
 *
 * Purpose:   Temporary (I hope): make all insert scores 0.
 *            If you let inserts train on the data, you can get
 *            positive insert emission scores. Local alignments,
 *            in particular, can then consist of just a couple of
 *            consensus states and a long string of insert 
 *            states, hitting base-composition-biased sequence
 *            with very high score. This is a Bad Thing.
 *            
 *            The long term solution for this problem will
 *            go in with mixture Dirichlet priors, but for now
 *            (with only Laplace coded), this'll appease the
 *            pitchfork and torches mob at Cambridge.
 *
 * Args:      cm - the model 
 *
 * Returns:   (void)
 *
 * Xref:      STL6 p.93.
 */
void
CMHackInsertScores(CM_t *cm)
{
  int v, x;
  for (v = 0; v < cm->M; v++)
    {
      if (cm->sttype[v] == IL_st || cm->sttype[v] == IR_st)
	for (x = 0; x < cm->abc->K; x++)
	  {
	    cm->esc[v][x]  = 0.;
	    cm->iesc[v][x] = 0.;
	  }
    }
  if(cm->cp9 != NULL)
    CP9HackInsertScores(cm->cp9);
}


/* Function: CountStatetype(), CMSubtreeCountStatetype(), CMSegmentCountStatetype
 * Date:     SRE, Wed Aug  2 09:15:00 2000 [St. Louis]
 *
 * Purpose:  Conveniences for counting the # of occurrences
 *           of a particular state type in a CM. Useful for
 *           "how many bifurcations does this model have", etc.
 *          
 *           CMSubtreeCountStatetype() only counts underneath     
 *           a particular subtree rooted at state v
 *
 * Args:     cm   - the model
 *           r    - the root of the subtree to start from (inclusive)
 *           z    - end of the subtree to stop at (inclusive) 
 *           type - a state type (e.g. E_st or MP_st)    
 *
 * Returns:  how many states of that type are in the model
 */
int
CMSegmentCountStatetype(CM_t *cm, int r, int z, char type)
{
  int count = 0;
  int v;
  for (v = r; v <= z; v++) 
    if (cm->sttype[v] == type) count++;
  return count;
}
int
CMSubtreeCountStatetype(CM_t *cm, int v, char type)
{
  int unsatisfied_starts = 1;
  int count = 0;

  while (unsatisfied_starts) {
    if (cm->sttype[v] == B_st) unsatisfied_starts++;
    if (cm->sttype[v] == E_st) unsatisfied_starts--; 
    if (cm->sttype[v] == type) count++;
    v++;
  }
  return count;
}
int
CMCountStatetype(CM_t *cm, char type)
{
  return CMSubtreeCountStatetype(cm, 0, type);
}
int 
CMSubtreeFindEnd(CM_t *cm, int r)
{
  int unsatisfied_starts = 1;

  while (unsatisfied_starts) {
    if (cm->sttype[r] == B_st) unsatisfied_starts++;
    if (cm->sttype[r] == E_st) unsatisfied_starts--; 
    r++;
  }
  return (r-1);
}

/* Function: CalculateStateIndex()
 * Date:     SRE, Mon Jul 31 15:37:55 2000 [St. Louis]
 *
 * Purpose:  Given a node index and a unique state type, use the CM's
 *           nodemap to calculate and return a state index in the CM.
 *
 *           Doesn't check that the node type matches what's implied
 *           by the utype! (e.g., if you pass utype==MATP_MP, the node
 *           had better be a MATP.)
 *
 * Args:     cm     - the covariance model
 *           node   - node index, 0..cm->nodes-1
 *           utype  - unique statetype, e.g. MATP_MP
 *
 * Returns:  a state index, 0..cm->M-1
 *
 * Used in:  modelmaker.c:transmogrify() 
 */
int
CalculateStateIndex(CM_t *cm, int node, char utype)
{
  int base;

  base = cm->nodemap[node];
  switch (utype) {
  case ROOT_S:  return base;
  case ROOT_IL: return base+1;
  case ROOT_IR: return base+2;
  case BEGL_S:  return base;
  case BEGR_S:  return base;
  case BEGR_IL: return base+1;
  case MATP_MP: return base;
  case MATP_ML: return base+1;
  case MATP_MR: return base+2;
  case MATP_D:  return base+3;  
  case MATP_IL: return base+4;
  case MATP_IR: return base+5; 
  case MATL_ML: return base;
  case MATL_D:  return base+1;
  case MATL_IL: return base+2;
  case MATR_MR: return base;
  case MATR_D:  return base+1;
  case MATR_IR: return base+2;
  case END_E:   return base;
  case BIF_B:   return base;
  default: esl_fatal("bogus utype %d in CalculateStateIndex()", utype);
  }
  return base;			/* not used */
}

/* Function:  TotalStatesInNode(), SplitStatesInNode(), InsertStatesInNode()
 * Incept:    SRE, Thu Aug  8 09:57:59 2002 [St. Louis]
 *
 * Purpose:   Returns the number of states in a node type.
 *
 * Args:      ndtype  - type of node (cm->ndtype[])
 */
int
TotalStatesInNode(int ndtype)
{
  switch (ndtype) {
  case BIF_nd:   return 1;
  case MATP_nd:  return 6;
  case MATL_nd:  return 3;
  case MATR_nd:  return 3;
  case BEGL_nd:  return 1;
  case BEGR_nd:  return 2;
  case ROOT_nd:  return 3;
  case END_nd:   return 1;
  default:       esl_fatal("Bogus node type %d", ndtype);
  }
  return 0;/*NOTREACHED*/
}
int
SplitStatesInNode(int ndtype)
{
  switch (ndtype) {
  case BIF_nd:   return 1;
  case MATP_nd:  return 4;
  case MATL_nd:  return 2;
  case MATR_nd:  return 2;
  case BEGL_nd:  return 1;
  case BEGR_nd:  return 1;
  case ROOT_nd:  return 1;
  case END_nd:   return 1;
  default:       esl_fatal("Bogus node type %d", ndtype);
  }
  return 0;/*NOTREACHED*/
}
int
InsertStatesInNode(int ndtype)
{
  switch (ndtype) {
  case BIF_nd:   return 0;
  case MATP_nd:  return 2;
  case MATL_nd:  return 1;
  case MATR_nd:  return 1;
  case BEGL_nd:  return 0;
  case BEGR_nd:  return 1;
  case ROOT_nd:  return 2;
  case END_nd:   return 0;
  default:       esl_fatal("Bogus node type %d", ndtype);
  }
  return 0;/*NOTREACHED*/
}




/* Function:  StateDelta(), StateLeftDelta(), StateRightDelta()
 * Incept:    SRE, Thu Oct  9 11:23:13 2003 [St. Louis]
 *
 * Purpose:   Convenience functions, mirroring some notation in Durbin et al.
 *            and elsewhere. \Delta notation simplifies some expositions
 *            of dynamic programming code.
 *            
 *            \Delta^R_v = 1 if the state emits right; else 0
 *            \Delta^L_v = 1 if the state emits left;  else 0
 *            \Delta_v   = 2 for pairwise, 1 for singlet, 0 for mute states.
 *            
 *            B_st, EL_st are special cases - Delta is returned as zero,
 *            but can't be used the same way.                                
 *
 * Args:      sttype   - state type code, e.g. MP_st
 *
 * Returns:   (see above)
 */
int
StateDelta(int sttype)
{
  switch (sttype) {
  case D_st:  return 0;
  case MP_st: return 2;
  case ML_st: return 1;
  case MR_st: return 1;
  case IL_st: return 1;
  case IR_st: return 1;
  case S_st:  return 0;
  case E_st:  return 0;
  case B_st:  return 0;
  case EL_st: return 0;
  default: esl_fatal("bogus state type %d\n", sttype);
  }
  /*NOTREACHED*/
  return 0;
}
int
StateLeftDelta(int sttype)
{
  switch (sttype) {
  case D_st:  return 0;
  case MP_st: return 1;
  case ML_st: return 1;
  case MR_st: return 0;
  case IL_st: return 1;
  case IR_st: return 0;
  case S_st:  return 0;
  case E_st:  return 0;
  case B_st:  return 0;
  case EL_st: return 0;
  default: esl_fatal("bogus state type %d\n", sttype);
  }
  /*NOTREACHED*/
  return 0;
}
int
StateRightDelta(int sttype)
{
  switch (sttype) {
  case D_st:  return 0;
  case MP_st: return 1;
  case ML_st: return 0;
  case MR_st: return 1;
  case IL_st: return 0;
  case IR_st: return 1;
  case S_st:  return 0;
  case E_st:  return 0;
  case B_st:  return 0;
  case EL_st: return 0;
  default: esl_fatal("bogus state type %d\n", sttype);
  }
  /*NOTREACHED*/
  return 0;
}

/* Function: PrintCM()
 * Date:     SRE, Sat Jul 29 10:55:16 2000 [St. Louis]
 *
 * Purpose:  Debugging: show a tabular representation of a CM structure.
 *
 * Args:     fp - output stream (e.g. stdout)
 *           cm - the CM to show
 *
 * Returns:  (void)
 */
void
PrintCM(FILE *fp, CM_t *cm)
{
  int x;

  fprintf(fp, "%5s %6s %5s %6s %7s %6s %5s %5s %5s\n",
	  " idx ","sttype", "ndidx", "ndtype", "  stid ", "cfirst", " cnum", "plast", " pnum");
  fprintf(fp, "%5s %6s %5s %6s %7s %5s %5s %5s %5s\n",
	  "-----", "------", "-----", "------","-------","------","-----", "-----", "-----");
  
  for (x = 0; x < cm->M; x++)
    {
      fprintf(fp, "%5d %-6s %5d %6s %-7s %6d %5d %5d %5d\n",
	      x, Statetype(cm->sttype[x]), cm->ndidx[x], 
	      Nodetype(cm->ndtype[cm->ndidx[x]]), UniqueStatetype(cm->stid[x]),
	      cm->cfirst[x], cm->cnum[x],
	      cm->plast[x], cm->pnum[x]);
    }
}

/* Function: SummarizeCM()
 * Date:     SRE, Sat Jul 29 12:19:31 2000 [St. Louis]
 *
 * Purpose:  Print some summary information about a new CM;
 *           called by cmbuild after each new model construction.
 *
 * Args:     fp - output stream (e.g. stdout)
 *           cm - cm to summarize
 *
 * Returns:  (void)
 */
void
SummarizeCM(FILE *fp, CM_t *cm)
{
  int x;
  int count[UNIQUESTATES];

  for (x = 0; x < UNIQUESTATES; x++) count[x] = 0;

  for (x = 0; x < cm->M; x++)
    count[(int) cm->stid[x]]++;
  
  fprintf(fp, "Summary report for CM structure:\n");
  fprintf(fp, "--------------------------------------\n");
  fprintf(fp, "Total states:       %4d\n", cm->M);
  fprintf(fp, "Total nodes:        %4d\n", cm->nodes);
  fprintf(fp, "Bifurcations:       %4d\n", count[BIF_B]);
  fprintf(fp, "MATP nodes:         %4d\n", count[MATP_MP]);
  fprintf(fp, "MATL nodes:         %4d\n", count[MATL_ML]);
  fprintf(fp, "MATR nodes:         %4d\n", count[MATR_MR]);
  fprintf(fp, "Consensus columns:  %4d    (2*MATP+MATL+MATR)\n",
	  count[MATP_MP]*2+count[MATL_ML]+count[MATR_MR]);
  fprintf(fp, "Base pairs:         %4d    (MATP)\n", count[MATP_MP]);
  fprintf(fp, "Single stranded:    %4d    (MATL+MATR)\n", count[MATL_ML]+count[MATR_MR]);
  /*fprintf(fp, "W: max hit size:    %d\n", cm->W);*/

}

/* Functions: Statetype(), Nodetype(), UniqueStatetype();
 *            StateCode(), NodeCode(), UniqueStateCode()
 * Date:      SRE, Sat Jul 29 11:07:47 2000 [St. Louis]
 *
 * Purpose:   Translate internal flags into human-readable strings, 
 *            for clearer debugging output (*type functions);
 *            or vice versa (*Code functions)
 * 
 * Args:      type - a state type, node type, or unique statetype code
 *            s    - string representing a code
 *
 * Returns:   the appropriate string; or the appropriate code.
 */
char *
Statetype(int type) 
{
  switch (type) {
  case D_st:  return "D";
  case MP_st: return "MP";
  case ML_st: return "ML";
  case MR_st: return "MR";
  case IL_st: return "IL";
  case IR_st: return "IR";
  case S_st:  return "S";
  case E_st:  return "E";
  case B_st:  return "B";
  case EL_st: return "EL";
  default: esl_fatal("bogus state type %d\n", type);
  }
  return ""; /*NOTREACHED*/
}
int
StateCode(char *s)
{
  if      (strcmp(s, "D")  == 0) return D_st;
  else if (strcmp(s, "MP") == 0) return MP_st;
  else if (strcmp(s, "ML") == 0) return ML_st;
  else if (strcmp(s, "MR") == 0) return MR_st;
  else if (strcmp(s, "IL") == 0) return IL_st;
  else if (strcmp(s, "IR") == 0) return IR_st;
  else if (strcmp(s, "S")  == 0) return S_st;
  else if (strcmp(s, "E")  == 0) return E_st;
  else if (strcmp(s, "B")  == 0) return B_st;
  else if (strcmp(s, "EL") == 0) return EL_st;
  return -1;
}
char *
Nodetype(int type) 
{
  switch (type) {
  case DUMMY_nd: return "-";
  case BIF_nd:   return "BIF";
  case MATP_nd:  return "MATP";
  case MATL_nd:  return "MATL";
  case MATR_nd:  return "MATR";
  case BEGL_nd:  return "BEGL";
  case BEGR_nd:  return "BEGR";
  case ROOT_nd:  return "ROOT";
  case END_nd:   return "END";
  default: esl_fatal("bogus node type %d\n", type);
  }
  return "";
}
int
NodeCode(char *s)
{
  if      (strcmp(s, "BIF")  == 0) return BIF_nd;
  else if (strcmp(s, "MATP") == 0) return MATP_nd;
  else if (strcmp(s, "MATL") == 0) return MATL_nd;
  else if (strcmp(s, "MATR") == 0) return MATR_nd;
  else if (strcmp(s, "BEGL") == 0) return BEGL_nd;
  else if (strcmp(s, "BEGR") == 0) return BEGR_nd;
  else if (strcmp(s, "ROOT") == 0) return ROOT_nd;
  else if (strcmp(s, "END")  == 0) return END_nd;
  return -1;
}
char *
UniqueStatetype(int type)
{
  switch (type) {
  case DUMMY:   return "DUMMY";   
  case ROOT_S:  return "ROOT_S";
  case ROOT_IL: return "ROOT_IL";
  case ROOT_IR: return "ROOT_IR";
  case BEGL_S : return "BEGL_S";
  case BEGR_S : return "BEGR_S";
  case BEGR_IL: return "BEGR_IL";
  case MATP_MP: return "MATP_MP";
  case MATP_ML: return "MATP_ML";
  case MATP_MR: return "MATP_MR";
  case MATP_D : return "MATP_D";
  case MATP_IL: return "MATP_IL";
  case MATP_IR: return "MATP_IR";
  case MATL_ML: return "MATL_ML";
  case MATL_D : return "MATL_D";
  case MATL_IL: return "MATL_IL";
  case MATR_MR: return "MATR_MR";
  case MATR_D : return "MATR_D";
  case MATR_IR: return "MATR_IR";
  case END_E  : return "END_E";
  case BIF_B  : return "BIF_B";
  case END_EL : return "END_EL";
  default: esl_fatal("bogus unique state type %d\n", type);
  }
  return "";
}
int
UniqueStateCode(char *s)
{
  if      (strcmp(s, "ROOT_S")  == 0) return ROOT_S;
  else if (strcmp(s, "ROOT_IL") == 0) return ROOT_IL;  
  else if (strcmp(s, "ROOT_IR") == 0) return ROOT_IR;  
  else if (strcmp(s, "BEGL_S")  == 0) return BEGL_S;  
  else if (strcmp(s, "BEGR_S")  == 0) return BEGR_S;  
  else if (strcmp(s, "BEGR_IL") == 0) return BEGR_IL;  
  else if (strcmp(s, "MATP_MP") == 0) return MATP_MP;  
  else if (strcmp(s, "MATP_ML") == 0) return MATP_ML;  
  else if (strcmp(s, "MATP_MR") == 0) return MATP_MR;  
  else if (strcmp(s, "MATP_D")  == 0) return MATP_D;  
  else if (strcmp(s, "MATP_IL") == 0) return MATP_IL;  
  else if (strcmp(s, "MATP_IR") == 0) return MATP_IR;  
  else if (strcmp(s, "MATL_ML") == 0) return MATL_ML;  
  else if (strcmp(s, "MATL_D")  == 0) return MATL_D;  
  else if (strcmp(s, "MATL_IL") == 0) return MATL_IL;  
  else if (strcmp(s, "MATR_MR") == 0) return MATR_MR;  
  else if (strcmp(s, "MATR_D")  == 0) return MATR_D;  
  else if (strcmp(s, "MATR_IR") == 0) return MATR_IR;  
  else if (strcmp(s, "BIF_B")   == 0) return BIF_B;
  else if (strcmp(s, "END_E")   == 0) return END_E;
  else if (strcmp(s, "END_EL")  == 0) return END_EL;
  else esl_fatal("bogus unique statetype %s\n", s);
  return 0; /*NOTREACHED*/
}
int
DeriveUniqueStateCode(int ndtype, int sttype)
{
  switch (ndtype) {
  case BIF_nd:   
    switch (sttype) {
    case B_st:  return BIF_B;
    default:    return -1;
    }
  case MATP_nd:  
    switch (sttype) {
    case D_st:  return MATP_D;
    case MP_st: return MATP_MP;
    case ML_st: return MATP_ML;
    case MR_st: return MATP_MR;
    case IL_st: return MATP_IL;
    case IR_st: return MATP_IR;
    default:    return -1;
    }
  case MATL_nd:  
    switch (sttype) {
    case D_st:  return MATL_D;
    case ML_st: return MATL_ML;
    case IL_st: return MATL_IL;
    default:    return -1;
    }
  case MATR_nd:  
    switch (sttype) {
    case D_st:  return MATR_D;
    case MR_st: return MATR_MR;
    case IR_st: return MATR_IR;
    default:    return -1;
    }
  case BEGL_nd:  
    switch (sttype) {
    case S_st:  return BEGL_S;
    default:    return -1;
    }
  case BEGR_nd:  
    switch (sttype) {
    case S_st:  return BEGR_S;
    case IL_st: return BEGR_IL;
    default:    return -1;
    }
  case ROOT_nd:  
    switch (sttype) {
    case S_st:  return ROOT_S;
    case IL_st: return ROOT_IL;
    case IR_st: return ROOT_IR;
    default:    return -1;
    }
  case END_nd:   
    switch (sttype) {
    case E_st:  return END_E;
    default:    return -1;
    }
  default: 
    return -1;
  }
}

/* Function: CMRebalance()
 * Date:     SRE, Mon Apr  8 11:40:46 2002 [St. Louis]
 *
 * Purpose:  Rebalance a CM tree to guarantee O(N^2 log N) memory in
 *           smallcyk.c's divide and conquer algorithm.
 * 
 *           Input: a CM that's numbered in preorder traversal: 
 *           visit root, visit left, visit right. (e.g., left
 *           child S always visited before right child S, 
 *           cfirst[w] < cnum[y], as produced by modelmaker.c).
 *           
 *           Output: a renumbered CM, in a modified preorder traversal:
 *           visit root, visit min weight child, visit max weight child,
 *           where weight is the # of extra CYK decks that'll need to
 *           be held in memory to calculate this subgraph.
 *           
 * Args:     cm - the old CM
 *
 * Returns:  A new CM. 
 *           Caller is responsible for free'ing this with FreeCM().
 */
CM_t *
CMRebalance(CM_t *cm)
{
  int       status;
  ESL_STACK *pda = NULL;  /* stack used for traversing old CM */
  CM_t     *new;          /* new CM we're creating */
  int      *wgt;          /* # of extra CYK decks required to calc subgraphs */
  int      *newidx;       /* newidx[v] = old CM state v's new index in new CM */
  int       v, w, y,z;	  /* state indices in old CM */
  int       nv;		  /* state index in new CM */
  int       x;		  /* counter over transitions, residues, nodes */

  /* Create the new model. Copy information that's unchanged by
   * renumbering the CM.
   */
  new = CreateCM(cm->nodes, cm->M, cm->abc);
  esl_strdup(cm->name, -1, &(new->name));
  esl_strdup(cm->acc,  -1, &(new->acc));
  esl_strdup(cm->desc, -1, &(new->desc));
  new->flags = cm->flags;
  new->clen  = cm->clen;
  new->nseq     = cm->nseq;
  new->eff_nseq = cm->eff_nseq;
  if(cm->flags & CMH_GA) new->ga = cm->ga;
  if(cm->flags & CMH_TC) new->tc = cm->tc;
  if(cm->flags & CMH_NC) new->nc = cm->nc;

  for (x = 0; x < cm->abc->K; x++) new->null[x] = cm->null[x];

  /* Calculate "weights" (# of required extra decks) on every B and S state.
   * Recursive rule here is: 1 + min(wgt[left], wgt[right]).
   */
  ESL_ALLOC(wgt, sizeof(int) * cm->M);
  for (v = cm->M-1; v >= 0; v--) 
    {
      if      (cm->sttype[v] == E_st) /* initialize unbifurcated segments with 1 */
	wgt[v] = 1; 
      else if (cm->sttype[v] == B_st) /* "cfirst"=left S child. "cnum"=right S child. */
	wgt[v] = 1 + ESL_MIN(wgt[cm->cfirst[v]], wgt[cm->cnum[v]]);
      else 
	wgt[v] = wgt[v+1];            /* all other states propagate up to S */
    }

  /* Now, preorder traverse the new CM. At each bifurcation, we want
   * to visit the S with minimum weight first. v is an index on the
   * old CM, and we hop it around using this traversal order and a
   * pushdown stack. nv is an index on the new CM, which just moves
   * in preorder traversal 0..cm->M-1.
   * 
   */
  v = 0;
  z = cm->M-1;
  pda = esl_stack_ICreate();
  ESL_ALLOC(newidx, sizeof(int) * cm->M);
  for (nv = 0; nv < cm->M; nv++)
    {    
      /* Keep a map of where the old states are going in new CM 
       * old state v becomes newidx[v] in the new model.
       * This is guaranteed to be a one to one map.
       */
      newidx[v] = nv;		

      /* Copy old v to new nv. 
       * First, the easy stuff, that's unaffected by renumbering.
       */
      new->sttype[nv] = cm->sttype[v];
      new->ndidx[nv]  = cm->ndidx[v];
      new->stid[nv]   = cm->stid[v];
      new->pnum[nv]   = cm->pnum[v];
      for (x = 0; x < MAXCONNECT; x++) {
	new->t[nv][x]   = cm->t[v][x];
	new->tsc[nv][x] = cm->t[v][x];
      }
      for (x = 0; x < cm->abc->K*cm->abc->K; x++) {
	new->e[nv][x] = cm->e[v][x];
	new->esc[nv][x] = cm->esc[v][x];
      }

      /* Slightly harder - the plast connection for nv, to the last
       * of 1-6 parent states. We use the newidx map to get it from plast[v].
       */
      if (nv != 0) new->plast[nv] = newidx[cm->plast[v]];
      else         new->plast[nv] = -1;	/* ROOT. */

      /* Now, figure out next v, and make cfirst, cnum connections.
       * 
       * If we're a B, then traverse to the lighter child S state first.
       * Remember the overload in CM struct: cfirst = idx of left child; 
       * cnum = idx of right child. So if we visit left w first, cfirst=nv+1; 
       * if we visit right y first, cnum=nv+1. Getting the second child
       * index is a little tricky: we rely on knowing that 
       * the # of states in the first subgraph we visit is y-w,
       * so we know the second child index is nv+y-w+1.
       * 
       * If we're an E, pop the next v off the stack. cfirst=-1,cnum=0, because
       * it has no children.
       * 
       * Else, the next v is just v++. cfirst for new nv can be calculated by using the
       * offset in the old model: e.g. nv + (cfirst[v] - v). cnum is unchanged.
       * 
       */
      if (cm->sttype[v] == B_st) 
	{
	  w = cm->cfirst[v];	/* left child of v*/
	  y = cm->cnum[v];	/* right child of v*/

	  if (wgt[w] <= wgt[y])	/* left (w) lighter or same weight? visit w first, defer y */
	    { 
	      esl_stack_IPush(pda, y); 
	      esl_stack_IPush(pda, z);
	      v = w; 
	      z = y-1;
	      new->cfirst[nv] = nv+1;     /* left child is nv+1 */
	      new->cnum[nv]   = nv+y-w+1; 
	    }  
	  else			/* right (y) lighter? visit y first, defer w */
	    { 
	      esl_stack_IPush(pda, w); 
	      esl_stack_IPush(pda, y-1);
	      v = y;		/* z unchanged. */
	      new->cfirst[nv] = nv+z-y+2; 
	      new->cnum[nv]   = nv+1;     /* right child is nv+1 */
	    }
	}
      else if (cm->sttype[v] == E_st) 
	{
	  new->cfirst[nv] = -1;
	  new->cnum[nv]   = 0;
	  esl_stack_IPop(pda, &z);
	  esl_stack_IPop(pda, &v);
	}
      else	
	{
	  new->cfirst[nv] = nv + (cm->cfirst[v]-v); /* use offset in old model */
	  new->cnum[nv]   = cm->cnum[v];            /* cnum unchanged. */
	  v++;
	}
    }

  /* Deal with the renumbered begin and end transition distributions,
   * using the newidx[v] map.
   */
  for (v = 0; v < cm->M; v++)
    {
      new->begin[newidx[v]] = cm->begin[v];
      new->end[newidx[v]]   = cm->end[v];
    }

  /* Guide tree numbering is unchanged - still in preorder.
   * Associate nodes with new state numbering.
   */
  for (x = 0; x < new->nodes; x++) 
    {
      new->nodemap[x] = newidx[cm->nodemap[x]];
      new->ndtype[x]  = cm->ndtype[x];
    }

  free(wgt);
  free(newidx);
  esl_stack_Destroy(pda);
  return new;

 ERROR:
  esl_fatal("Memory allocation error.\n");
  return NULL; /* never reached */
}

/*
 * EPN, Wed Mar 21 09:29:55 2007
 * 
 * rsearch_calculate_gap_penalty (FROM RSEARCH::buildcm.c)
 *
 * Given the from state, the to state, and the gap parameters, returns
 * the gap penalty.
 */
float rsearch_calculate_gap_penalty (char from_state, char to_state, 
				     int from_node, int to_node, 
				     float input_alpha, float input_beta, 
				     float input_alphap, float input_betap) {
  int from_class, to_class;
  double alpha, beta;          /* Alpha or beta values to use */

  /* There are potentially 400 different combinations of state pairs here.
     To make it manageable, break down into the 6 classes on p. 8 of lab
     book 7, numbered as follows
     0    M    ROOT_S, BEGL_S, BEGR_S, MATP_MP, MATL_ML, MATR_MR, END_E, BIF_B
     1    IL   ROOT_IL, BEGR_IL, MATP_IL, MATL_IL
     2    DL   MATP_MR, MATL_D
     3    IR   ROOT_IR, MATP_IR, MATR_IR
     4    DR   MATP_ML, MATR_D
     5    DB   MATP_D
  */
  switch (from_state) {
  case MATP_D:
    from_class = DB_cl;
    break;
  case MATP_ML:
  case MATR_D:
    from_class = DR_cl;
    break;
  case ROOT_IR:
  case MATP_IR:
  case MATR_IR:
    from_class = IR_cl;
    break;
  case MATP_MR:
  case MATL_D:
    from_class = DL_cl;
    break;
  case ROOT_IL:
  case BEGR_IL:
  case MATP_IL:
  case MATL_IL:
    from_class = IL_cl;
    break;
  default:
    from_class = M_cl;
  }

  switch (to_state) {
  case MATP_D:
    to_class = DB_cl;
    break;
  case MATP_ML:
  case MATR_D:
    to_class = DR_cl;
    break;
  case ROOT_IR:
  case MATP_IR:
  case MATR_IR:
    to_class = IR_cl;
    break;
  case MATP_MR:
  case MATL_D:
    to_class = DL_cl;
    break;
  case ROOT_IL:
  case BEGR_IL:
  case MATP_IL:
  case MATL_IL:
    to_class = IL_cl;
    break;
  default:
    to_class = M_cl;
  }

  /* Now set alpha and beta according to state classes and nodes */
  /* Alpha is alpha' for MATP->MATP, alpha otherwise */
  if (from_node == MATP_nd && to_node == MATP_nd)
    alpha = input_alphap;
  else
    alpha = input_alpha;
  /* Beta is beta' iff from_cl is DB and MATP->MATP */
  if (from_class == DB_cl && from_node == MATP_nd && to_node == MATP_nd)
    beta = input_betap;
  else 
    beta = input_beta;

  /* Now that we have the proper class, return the appropriate gap penalty */
  if (from_class == M_cl) {
    if (to_class == M_cl) {
      return (0.);
    } else if (to_class == DB_cl) {
      return (alpha);
    } else {
      return (0.5*alpha);
    }
  } else if (from_class == IL_cl) {
    if (to_class == M_cl) {
      return (beta + 0.5*alpha);
    } else if (to_class == IL_cl) {
      return (beta);
    } else if (to_class == DB_cl) {
      return (beta+1.5*alpha);
    } else {
      return (beta + alpha);
    }
  } else if (from_class == DL_cl) {
    if (to_class == M_cl) {
      return (beta + 0.5*alpha);
    } else if (to_class == DL_cl) {
      return (beta);
    } else if (to_class == DB_cl) {
      return (beta + 0.5*alpha);
    } else {
      return (beta + alpha);
    }
  } else if (from_class == IR_cl) {
    if (to_class == M_cl) {
      return (beta + 0.5*alpha);
    } else if (to_class == IR_cl) {
      return (beta);
    } else if (to_class == DB_cl) {
      return (beta+1.5*alpha);
    } else {
      return (beta + alpha);
    }
  } else if (from_class == DR_cl) {
    if (to_class == M_cl) {
      return (beta + 0.5*alpha);
    } else if (to_class == DR_cl) {
      return (beta);
    } else if (to_class == DB_cl) {
      return (beta + 0.5*alpha);
    } else {
      return (beta + alpha);
    }
  } else {                /* DB_cl */
    if (to_class == IL_cl || to_class == IR_cl) {
      return (2*beta + 1.5*alpha);
    } else if (to_class == M_cl) {
      return (2*beta + alpha);
    } else if (to_class == DB_cl) {
      return (2*beta);
    } else {
      return (2*beta + 0.5*alpha);
    }
  }
  return (0);
}

/*
 * Function: ExponentiateCM
 * Date:     EPN, Sun May 20 13:10:06 2007
 * Purpose:  Exponentiate the emission and transition probabilities 
 *           of a CM by z. If CM is in local mode, put it in global mode,
 *           exponentiate it and put it back in local mode, otherwise
 *           the cm->end probabilities would change, and we don't want
 *           that.
 * 
 * Args:
 *           CM           - the covariance model
 *           z            - factor to exponentiate by
 */
int
ExponentiateCM(CM_t *cm, double z)
{
  /*printf("in ExponentiateCM, z: %f\n", z);*/
  int v;
  int x,y;
  int local_flag = FALSE;

  /* If in local mode, configure to global first. */
  if(cm->flags & CM_LOCAL_BEGIN || cm->flags & CM_LOCAL_END) 
    {
      ConfigGlobal(cm);
      local_flag = TRUE;
    }

  for(v = 0; v < cm->M; v++)
    {
      if (cm->sttype[v] != B_st && cm->sttype[v] != E_st)
	for (x = 0; x < cm->cnum[v]; x++)
	  cm->t[v][x]  = pow(cm->t[v][x], z);
      if (cm->sttype[v] == MP_st)
	{
	  for (x = 0; x < cm->abc->K; x++)
	    for (y = 0; y < cm->abc->K; y++)
	      cm->e[v][x*cm->abc->K+y]  = pow(cm->e[v][x*cm->abc->K+y], z);
	}
      if (cm->sttype[v] == ML_st || cm->sttype[v] == MR_st ||
	  cm->sttype[v] == IL_st || cm->sttype[v] == IR_st)
	{
	  for (x = 0; x < cm->abc->K; x++)
	    cm->e[v][x]  = pow(cm->e[v][x], z);
	}
    }

  CMRenormalize(cm);

  if(local_flag) ConfigLocal(cm, cm->pbegin, cm->pend);
  /* new probs invalidate log odds scores */
  cm->flags &= ~CMH_BITS;
  return eslOK;
}


/*
 * Function: DuplicateCM
 * Date:     EPN, Thu May 24 09:57:12 2007
 * Purpose:  Given a template CM 'cm', copy it's params into a
 *           new CM which is allocated here and must be
 *           freed by the caller. 
 * 
 * Args:
 *           src          - the template covariance model
 */
CM_t *
DuplicateCM(CM_t *cm)
{
  int       status;
  int       v;	          /* counter over states */
  int       x;		  /* counter over transitions, residues, nodes */
  CM_t     *new;
  ESL_ALPHABET *abc;
  abc = esl_alphabet_Create(cm->abc->type);

  /* Create the new model and copy everything over except the cp9 and stats */
  new = CreateCM(cm->nodes, cm->M, cm->abc);
  esl_strdup(cm->name, -1, &(new->name));
  esl_strdup(cm->acc,  -1, &(new->acc));
  esl_strdup(cm->desc, -1, &(new->desc));
  new->flags       = cm->flags;
  new->search_opts = cm->search_opts;
  new->align_opts  = cm->align_opts;
  new->config_opts = cm->config_opts;

  new->nodes      = cm->nodes;
  for(x = 0; x < cm->nodes; x++)
    {
      new->nodemap[x]   = cm->nodemap[x];
      new->ndtype[x]    = cm->ndtype[x];
    }
  for(v = 0; v < cm->M; v++)
    {
      new->sttype[v]  = cm->sttype[v];
      new->ndidx[v]   = cm->ndidx[v];
      new->stid[v]    = cm->stid[v];

      new->cfirst[v]  = cm->cfirst[v];
      new->cnum[v]    = cm->cnum[v];

      new->pnum[v]    = cm->pnum[v];
      new->plast[v]   = cm->plast[v];

      new->begin[v]   = cm->begin[v];
      new->beginsc[v] = cm->beginsc[v];
      new->ibeginsc[v]= cm->ibeginsc[v];
      new->end[v]     = cm->end[v];
      new->endsc[v]   = cm->endsc[v];
      new->iendsc[v]  = cm->iendsc[v];

      /* copy transitions and emissions*/
      for (x = 0; x < MAXCONNECT; x++)
	{
	  new->t[v][x]     = cm->t[v][x];
	  new->tsc[v][x]   = cm->tsc[v][x];
	  new->itsc[v][x]  = cm->itsc[v][x];
	}
      for (x = 0; x < cm->abc->K * cm->abc->K; x++)
	{
	  new->e[v][x]     = cm->e[v][x];
	  new->esc[v][x]   = cm->esc[v][x];
	  new->iesc[v][x]  = cm->iesc[v][x];
	}
    }      
  if(cm->dmin != NULL && cm->dmax != NULL)
    {
      ESL_ALLOC(new->dmin, sizeof(int) * cm->M);
      ESL_ALLOC(new->dmax, sizeof(int) * cm->M);
      for(v = 0; v < cm->M; v++)
	{
	  new->dmin[v] = cm->dmin[v];
	  new->dmax[v] = cm->dmax[v];
	}
    }
  else 
    {
      new->dmin = new->dmax = NULL;
      new->flags &= ~CM_QDB;
    }
  new->W      = cm->W;
  new->el_selfsc  = cm->el_selfsc;
  new->iel_selfsc = cm->iel_selfsc;
  new->beta  = cm->beta;
  new->tau   = cm->tau;
  new->enf_start = cm->enf_start;
  if(cm->enf_seq != NULL)
    esl_strdup(cm->enf_seq, -1, &(new->enf_seq));
  else new->enf_seq = NULL;
  new->enf_scdiff = cm->enf_scdiff;
  new->sc_boost   = cm->sc_boost;
  new->ffract     = cm->ffract;
  new->cutoff_type= cm->cutoff_type;
  new->cutoff     = cm->cutoff;
  new->cp9_cutoff_type = cm->cp9_cutoff_type;
  new->cp9_cutoff = cm->cp9_cutoff;
  new->cp9_sc_boost = cm->cp9_sc_boost;
  if(cm->root_trans == NULL)
    new->root_trans = NULL;
  else
    {
      ESL_ALLOC(new->root_trans, sizeof(float) * cm->cnum[0]);
      for (v = 0; v < cm->cnum[0]; v++)
	new->root_trans[v] = cm->root_trans[v];
    }
  new->hmmpad = cm->hmmpad;

  new->cp9  = NULL;

  /* Copy the CM stats if they exist */
  if(cm->flags & CM_GUMBEL_STATS)
    {
      new->stats = AllocCMStats(cm->stats->np);
      CopyCMStats(cm->stats, new->stats);
    }

  /* Copy the CP9 if it exists */
  if(cm->flags & CM_CP9)
    {
      DuplicateCP9(cm, new);
      new->flags |= CM_CP9; /* raise the CP9 flag */
    }

  return new;

 ERROR:
  esl_fatal("Memory allocation error.\n");
  return NULL; /* never reached */
}

/* Function:  cm_banner()
 * Synopsis:  print standard INFERNAL application output header
 *            Based on p7_banner from HMMER3 dev code.
 * Incept:    EPN, Fri May 25 15:05:42 2007
 *
 * Purpose:   Print the standard INFERNAL command line application banner
 *            to <fp>, constructing it from <progname> (the name of the
 *            program) and a short one-line description <banner>.
 *              
 *            <progname> would typically be an application's
 *            <argv[0]>, rather than a fixed string. This allows the
 *            program to be renamed, or called under different names
 *            via symlinks. Any path in the <progname> is discarded;
 *            for instance, if <progname> is "/usr/local/bin/cmcalibrate",
 *            "cmcalibrate" is used as the program name.
 *            
 * Note:    
 *    Needs to pick up preprocessor #define's from config.h,
 *    as set by ./configure.
 *
 * Returns:   (void)
 */
void
cm_banner(FILE *fp, char *progname, char *banner)
{
  char *appname = NULL;

  if (esl_FileTail(progname, FALSE, &appname) != eslOK) appname = progname;

  fprintf(fp, "# %s :: %s\n", appname, banner);
  fprintf(fp, "# INFERNAL %s (%s)\n", PACKAGE_VERSION, PACKAGE_DATE);
  fprintf(fp, "# %s\n", PACKAGE_COPYRIGHT);
  fprintf(fp, "# %s\n", PACKAGE_LICENSE);
  fprintf(fp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");

  if (appname != NULL) free(appname);
  return;
}

/* Function:  cm_Validate()
 * Incept:    EPN, Fri Jul 27 14:59:55 2007 [Janelia]
 *
 * Purpose:   Validates some internals of the CM structure <cm>.
 * 
 *            Probability vectors are validated to sum up to
 *            within a fractional tolerance <tol> of 1.0.
 *
 *            Probably only useful for debugging and development,
 *            not production code.
 *
 * Returns:   <eslOK> if <cm> internals look fine.
 *            Returns <eslFAIL> if something is wrong.
 */
int
cm_Validate(CM_t *cm, float tol, char *errbuf)
{
  int status;
  int v;
  int clen = 0;

  if (cm             == NULL)       ESL_XFAIL(eslFAIL, errbuf, "CM is a null pointer");
  if (cm->M          <  1)          ESL_XFAIL(eslFAIL, errbuf, "CM has M < 1");
  if (cm->abc        == NULL)       ESL_XFAIL(eslFAIL, errbuf, "CM has no alphabet reference");
  if (cm->abc->type  == eslUNKNOWN) ESL_XFAIL(eslFAIL, errbuf, "CM's alphabet is set to unknown");
  
  for (v = 0; v < cm->M; v++)
    {
      if(StateDelta(cm->sttype[v]) == 2)
	{
	  if (esl_vec_FValidate(cm->e[v], (cm->abc->K * cm->abc->K), tol, NULL) != eslOK)
	    ESL_XFAIL(eslFAIL, errbuf, "e[%d] fails pvector validation", v); 
	}
      else if(StateDelta(cm->sttype[v]) > 0)
	{ 
	  if (esl_vec_FValidate(cm->e[v], cm->abc->K, tol, NULL) != eslOK) 
	    ESL_XFAIL(eslFAIL, errbuf, "e[%d] fails pvector validation", v); 
	}
      if (cm->sttype[v] != B_st && cm->sttype[v] != E_st)
	if(esl_vec_FValidate(cm->t[v], cm->cnum[v], tol, NULL) != eslOK) 
	  ESL_XFAIL(eslFAIL, errbuf, "t[%d] fails pvector validation", v);
      if(cm->stid[v] == MATL_ML) clen++;
      if(cm->stid[v] == MATR_MR) clen++;
      if(cm->stid[v] == MATP_MP) clen+=2;
    }
  if(cm->clen != clen) ESL_XFAIL(eslFAIL, errbuf, "consensus length %d not correctly stored in CM, should be %d", cm->clen, clen);
  
  return eslOK;

 ERROR:
  return status;
}

/* Function: CMStatetype()
 * 
 * Purpose:  Returns the CM state type in text.
 * Example:  CP9Statetype(MP_st) = "MP"
 */
char *
CMStatetype(char st)
{
  switch (st) {
  case D_st:  return "D";
  case MP_st: return "MP";
  case ML_st: return "ML";
  case MR_st: return "MR";
  case IL_st: return "IL";
  case IR_st: return "IR";
  case S_st:  return "S";
  case E_st:  return "E";
  case B_st:  return "B";
  case EL_st: return "EL";
  default: return "BOGUS";
  }
}

/* Function: CMNodetype()
 * 
 * Purpose:  Returns the CM state type in text.
 * Example:  CP9Statetype(MATP_nd) = "MATP"
 */
char *
CMNodetype(char nd)
{
  switch (nd) {
  case DUMMY_nd:  return "DUMMY";
  case BIF_nd:    return "BIF";
  case MATP_nd:   return "MATP";
  case MATL_nd:   return "MATL";
  case MATR_nd:   return "MATR";
  case BEGL_nd:   return "BEGL";
  case BEGR_nd:   return "BEGR";
  case ROOT_nd:   return "ROOT";
  case END_nd:    return "END";
  default: return "BOGUS";
  }
}

/* Function: CMStateid()
 * 
 * Purpose:  Returns the CM state id in text.
 * Example:  CP9Statetype(MATP_MP) = "MATP_MP"
 */
char *
CMStateid(char st)
{
  switch (st) {
  case DUMMY:   return "DUMMY";
  case ROOT_S:  return "ROOT_S";
  case ROOT_IL: return "ROOT_IL";
  case ROOT_IR: return "ROOT_IR";
  case BEGL_S:  return "BEGL_S";
  case BEGR_S:  return "BEGR_S";
  case BEGR_IL: return "BEGR_IL";
  case MATP_MP: return "MATP_MP";
  case MATP_ML: return "MATP_ML";
  case MATP_MR: return "MATP_MR";
  case MATP_D:  return "MATP_D";
  case MATP_IL: return "MATP_IL";
  case MATP_IR: return "MATP_IR";
  case MATL_ML: return "MATL_ML";
  case MATL_D:  return "MATL_D";
  case MATL_IL: return "MATL_IL";
  case MATR_MR: return "MATR_MR";
  case MATR_D:  return "MATR_D";
  case MATR_IR: return "MATR_IR";
  case END_E:   return "END_E";
  case BIF_B:   return "BIF_B";
  case END_EL:  return "END_EL";
  default: return "BOGUS";
  }
}


/*****************************************************************
 * Convenience routines for setting fields in an CM. (from p7_cm.c)
 *****************************************************************/ 
/* Function: cm_SetName()
 * Incept:   EPN, Fri Jul 27 16:49:49 2007 [Janelia]
 * 
 * Purpose:  Set or change the name of a CM to <name>.
 *           Any trailing whitespace (including newline) is chopped off.     
 *      
 * Returns:  <eslOK> on success.
 *
 * Throws:   <eslEMEM> on allocation error, and original name (if any) 
 *           remains.
 */
int
cm_SetName(CM_t *cm, char *name)
{
  int   status;
  void *tmp;
  int   n;

  if (name == NULL) {
    if (cm->name != NULL) free(cm->name); 
    cm->name = NULL;
  } else {
    n = strlen(name);
    ESL_RALLOC(cm->name, tmp, sizeof(char)*(n+1));
    strcpy(cm->name, name);
    if ((status = esl_strchop(cm->name, n)) != eslOK) goto ERROR;
  }
  return eslOK;

 ERROR:
  return status;
}

/* Function: cm_SetAccession()
 * Incept:   SRE, Mon Jan  1 16:53:53 2007 [Casa de Gatos]
 * 
 * Purpose:  Set or change the accession number of a Plan7 CM to <acc>,
 *           and raise the <ACC> flag. Trailing whitespace (including newline) 
 *           is chopped.  
 *           
 *           If <acc> is <NULL>, unset the CM's accession (if any) and drop 
 *           the <ACC> flag.
 *
 * Returns:  <eslOK> on success.
 *
 * Throws:   <eslEMEM> on allocation error, and original name (if any) 
 *           remains.
 */
int
cm_SetAccession(CM_t *cm, char *acc)
{
  int   status;
  void *tmp;
  int   n;

  if (acc == NULL) {
    if (cm->acc != NULL) free(cm->acc); 
    cm->acc = NULL;
    cm->flags &= ~CMH_ACC;
  } else {
    n = strlen(acc);
    ESL_RALLOC(cm->acc, tmp, sizeof(char)*(n+1));
    strcpy(cm->acc, acc);
    if ((status = esl_strchop(cm->acc, n)) != eslOK) goto ERROR;
    cm->flags |= CMH_ACC;
  }
  return eslOK;

 ERROR:
  return status;
}

/* Function: cm_SetDescription()
 * Incept:   SRE, Mon Jan  1 16:59:28 2007 [Casa de Gatos]
 * 
 * Purpose:  Set or change the description line of a Plan7 CM. 
 *           Trailing whitespace (including newline) is chopped.
 */
int
cm_SetDescription(CM_t *cm, char *desc)
{
  int   status;
  void *tmp;
  int   n;

  if (desc == NULL) 
    {
      if (cm->desc != NULL) free(cm->desc); 
      cm->desc   = NULL;
      cm->flags &= ~CMH_DESC;
    }
  else
    {
      n = strlen(desc);
      ESL_RALLOC(cm->desc, tmp, sizeof(char)*(n+1));
      strcpy(cm->desc, desc);
      if ((status = esl_strchop(cm->desc, n)) != eslOK) goto ERROR;
      cm->flags |= CMH_DESC;
    }
  return eslOK;

 ERROR:
  return status;
}

/* Function: cm_AppendComlog()
 * Incept:   SRE, Mon Jan  1 18:23:42 2007 [Casa de Gatos]
 * 
 * Purpose:  Concatenate command line options and append as a new line in the
 *           command line log. Command line log is multiline, with each line
 *           ending in newline char, except for last line.
 *           
 * Returns:  <eslOK> on success.
 * 
 * Throws:   <eslEMEM> on allocation failure.          
 */
int
cm_AppendComlog(CM_t *cm, int argc, char **argv)
{
  int   status;
  void *tmp;
  int   n;
  int   i;

  /* figure out length of added command line, and (re)allocate comlog */
  n = argc-1;	/* account for 1 space per arg, except last one */
  for (i = 0; i < argc; i++)
    n += strlen(argv[i]);

  if (cm->comlog != NULL) {
    n += strlen(cm->comlog) + 1; /* +1 for the \n we're going to add to the old comlog */
    ESL_RALLOC(cm->comlog, tmp, sizeof(char)* (n+1));
    strcat(cm->comlog, "\n");
  } else {
    ESL_ALLOC(cm->comlog, sizeof(char)* (n+1));
    *(cm->comlog) = '\0'; /* need this to make strcat work */
  }

  for (i = 0; i < argc-1; i++)
    {
      strcat(cm->comlog, argv[i]);
      strcat(cm->comlog, " ");
    }
  strcat(cm->comlog, argv[argc-1]);
  return eslOK;

 ERROR:
  return status;
}

/* Function: cm_SetCtime()
 * Date:     SRE, Wed Oct 29 11:53:19 1997 [TWA 721 over the Atlantic]
 * 
 * Purpose:  Set the <ctime> field in a new CM to the current time.
 *
 *           This function is not reentrant and not threadsafe, because
 *           it calls the nonreentrant ANSI C ctime() function.
 * 
 * Returns:  <eslOK> on success.
 * 
 * Throws:   <eslEMEM> on allocation failure. <eslESYS> if the <time()>
 *           system call fails to obtain the calendar time.
 */
int
cm_SetCtime(CM_t *cm)
{
  int    status;
  char  *s = NULL;
  time_t date;

  if ((date   = time(NULL))                       == -1) { status = eslESYS; goto ERROR; }
  if ((status = esl_strdup(ctime(&date), -1, &s)) != eslOK) goto ERROR;
  if ((status = esl_strchop(s, -1))               != eslOK) goto ERROR;
  
  if (cm->ctime != NULL) free(cm->ctime);
  cm->ctime = s;
  return eslOK;

 ERROR:
  if (s != NULL) free(s);
  return status;
}
/*---------------- end, internal-setting routines ---------------*/

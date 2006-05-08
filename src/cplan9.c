/************************************************************
 * @LICENSE@
 ************************************************************/
/* cplan9.c based on plan7.c
 * EPN 02.27.06
 * 
 * Support for CM-Plan 9 HMM data structure, cplan9_s.
 * Including support for a dp matrix structure for cplan_s,
 * cp9_dpmatrix_s.
 */

#include "hmmer_config.h"
#include "squidconf.h"
#include "cplan9.h"

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include "squid.h"
#include "funcs.h"
#include "structs.h"
#include "hmmer_funcs.h"
#include "hmmer_structs.h"


/* Functions: AllocCPlan9(), AllocCPlan9Shell(), AllocCPlan9Body(), FreeCPlan9()
 * 
 * Purpose:   Allocate or free a CPlan9 HMM structure.
 *            Can either allocate all at one (AllocCPlan9()) or
 *            in two steps (AllocPlan7Shell(), AllocPlan7Body()).
 *            The two step method is used in hmmio.c where we start
 *            parsing the header of an HMM file but don't 
 *            see the size of the model 'til partway thru the header.
 */
struct cplan9_s *
AllocCPlan9(int M) 
{
  struct cplan9_s *hmm;

  hmm = AllocCPlan9Shell();
  AllocCPlan9Body(hmm, M);
  return hmm;
}  
struct cplan9_s *
AllocCPlan9Shell(void) 
{
  struct cplan9_s *hmm;

  hmm    = (struct cplan9_s *) MallocOrDie (sizeof(struct cplan9_s));
  hmm->M = 0;

  hmm->name     = NULL;
  hmm->acc      = NULL;
  hmm->desc     = NULL;
  hmm->rf       = NULL;
  hmm->cs       = NULL;
  hmm->ca       = NULL;
  hmm->comlog   = NULL; 
  hmm->nseq     = 0;
  hmm->ctime    = NULL;
  hmm->map      = NULL;
  hmm->checksum = 0;

  hmm->tpri = NULL;
  hmm->mpri = NULL;
  hmm->ipri = NULL;

  hmm->ga1 = hmm->ga2 = 0.0;
  hmm->tc1 = hmm->tc2 = 0.0;
  hmm->nc1 = hmm->nc2 = 0.0;

  hmm->t      = NULL;
  hmm->mat    = NULL;
  hmm->ins    = NULL;
  
  hmm->tsc     = hmm->msc     = hmm->isc     = NULL;
  hmm->tsc_mem = hmm->msc_mem = hmm->msc_mem = NULL;

  hmm->begin  = NULL;
  hmm->end    = NULL;

  hmm->bsc = hmm->bsc_mem = NULL;
  hmm->esc = hmm->esc_mem = NULL;

				/* DNA translation is not enabled by default */
  hmm->dnam   = NULL;
  hmm->dnai   = NULL;
  hmm->dna2   = -INFTY;
  hmm->dna4   = -INFTY;
			/* statistical parameters set to innocuous empty values */
  hmm->mu     = 0.; 
  hmm->lambda = 0.;
  
  hmm->flags = 0;
  return hmm;
}  
#ifndef ALTIVEC  /* in Altivec port, replaced in fast_algorithms.c */
void
AllocCPlan9Body(struct cplan9_s *hmm, int M) 
{
  int k, x;

  hmm->M = M;

  hmm->rf     = MallocOrDie ((M+2) * sizeof(char));
  hmm->cs     = MallocOrDie ((M+2) * sizeof(char));
  hmm->ca     = MallocOrDie ((M+2) * sizeof(char));
  hmm->map    = MallocOrDie ((M+1) * sizeof(int));

  hmm->t      = MallocOrDie ((M+1) *           sizeof(float *));
  hmm->mat    = MallocOrDie ((M+1) *           sizeof(float *));
  hmm->ins    = MallocOrDie ((M+2) *           sizeof(float *));
  hmm->t[0]   = MallocOrDie ((9*(M+1))     *       sizeof(float));
  hmm->mat[0] = MallocOrDie ((HMMERMAXABET*(M+1)) * sizeof(float));
  hmm->ins[0] = MallocOrDie ((HMMERMAXABET*(M+1)) *     sizeof(float));

  hmm->tsc     = MallocOrDie (9     *           sizeof(int *));
  hmm->msc     = MallocOrDie (MAXCODE   *       sizeof(int *));
  hmm->isc     = MallocOrDie (MAXCODE   *       sizeof(int *)); 
  hmm->tsc_mem = MallocOrDie ((9*(M+2))     *       sizeof(int));
  hmm->msc_mem = MallocOrDie ((MAXCODE*(M+1)) * sizeof(int));
  hmm->isc_mem = MallocOrDie ((MAXCODE*(M+2)) *     sizeof(int));

  hmm->tsc[0] = hmm->tsc_mem;
  hmm->msc[0] = hmm->msc_mem;
  hmm->isc[0] = hmm->isc_mem;

  /* note allocation strategy for important 2D arrays -- trying
   * to keep locality as much as possible, cache efficiency etc.
   */
  for (k = 1; k <= M; k++) {
    hmm->mat[k] = hmm->mat[0] + k * HMMERMAXABET;
    hmm->ins[k] = hmm->ins[0] + k * HMMERMAXABET;
    hmm->t[k]   = hmm->t[0]   + k * 9;
  }
  for (x = 1; x < MAXCODE; x++) {
    hmm->msc[x] = hmm->msc[0] + x * (M+1);
    hmm->isc[x] = hmm->isc[0] + x * (M+1);
  }
  for (x = 0; x < 9; x++)
    hmm->tsc[x] = hmm->tsc[0] + x * (M+1);

  /* tsc[x][0] is used as a boundary condition sometimes [Viterbi()],
   * so set to -inf always.
   */
  for (x = 0; x < 9; x++)
    hmm->tsc[x][0] = -INFTY;

  hmm->begin  = MallocOrDie  ((M+1) * sizeof(float));
  hmm->end    = MallocOrDie  ((M+1) * sizeof(float));

  hmm->bsc_mem  = MallocOrDie  ((M+1) * sizeof(int));
  hmm->esc_mem  = MallocOrDie  ((M+1) * sizeof(int));

  hmm->bsc = hmm->bsc_mem;
  hmm->esc = hmm->esc_mem;

  return;
}  
#endif /* ALTIVEC */

void
FreeCPlan9(struct cplan9_s *hmm)
{
  if (hmm->name    != NULL) free(hmm->name);
  if (hmm->acc     != NULL) free(hmm->acc);
  if (hmm->desc    != NULL) free(hmm->desc);
  if (hmm->rf      != NULL) free(hmm->rf);
  if (hmm->cs      != NULL) free(hmm->cs);
  if (hmm->ca      != NULL) free(hmm->ca);
  if (hmm->comlog  != NULL) free(hmm->comlog);
  if (hmm->ctime   != NULL) free(hmm->ctime);
  if (hmm->map     != NULL) free(hmm->map);
  if (hmm->tpri    != NULL) free(hmm->tpri);
  if (hmm->mpri    != NULL) free(hmm->mpri);
  if (hmm->ipri    != NULL) free(hmm->ipri);
  if (hmm->bsc_mem != NULL) free(hmm->bsc_mem);
  if (hmm->begin   != NULL) free(hmm->begin);
  if (hmm->esc_mem != NULL) free(hmm->esc_mem);
  if (hmm->end     != NULL) free(hmm->end);
  if (hmm->msc_mem != NULL) free(hmm->msc_mem);
  if (hmm->isc_mem != NULL) free(hmm->isc_mem);
  if (hmm->tsc_mem != NULL) free(hmm->tsc_mem);
  if (hmm->mat     != NULL) free(hmm->mat[0]);
  if (hmm->ins     != NULL) free(hmm->ins[0]);
  if (hmm->t       != NULL) free(hmm->t[0]);
  if (hmm->msc     != NULL) free(hmm->msc);
  if (hmm->isc     != NULL) free(hmm->isc);
  if (hmm->tsc     != NULL) free(hmm->tsc);
  if (hmm->mat     != NULL) free(hmm->mat);
  if (hmm->ins     != NULL) free(hmm->ins);
  if (hmm->t       != NULL) free(hmm->t);
  if (hmm->dnam    != NULL) free(hmm->dnam);
  if (hmm->dnai    != NULL) free(hmm->dnai);
  free(hmm);
}

/* Function: ZeroCPlan9()
 * 
 * Purpose:  Zeros the counts/probabilities fields in a model.  
 *           Leaves null model untouched. 
 */
void
ZeroCPlan9(struct cplan9_s *hmm)
{
  int k;
  FSet(hmm->ins[0], Alphabet_size, 0.);
  FSet(hmm->t[0], 9, 0.);
  for (k = 1; k <= hmm->M; k++)
    {
      FSet(hmm->t[k], 9, 0.);
      FSet(hmm->mat[k], Alphabet_size, 0.);
      FSet(hmm->ins[k], Alphabet_size, 0.);
    }
  FSet(hmm->begin+1, hmm->M, 0.);
  FSet(hmm->end+1, hmm->M, 0.);
  hmm->flags &= ~CPLAN9_HASBITS;	/* invalidates scores */
  hmm->flags &= ~CPLAN9_HASPROB;	/* invalidates probabilities */
}


/* Function: CPlan9SetName()
 * 
 * Purpose:  Change the name of a CPlan9 HMM. Convenience function.
 *      
 * Note:     Trailing whitespace and \n's are chopped.     
 */
void
CPlan9SetName(struct cplan9_s *hmm, char *name)
{
  if (hmm->name != NULL) free(hmm->name);
  hmm->name = Strdup(name);
  StringChop(hmm->name);
}
/* Function: Cplan9SetAccession()
 * 
 * Purpose:  Change the accession number of a Cplan9 HMM. Convenience function.
 *      
 * Note:     Trailing whitespace and \n's are chopped.     
 */
void
CPlan9SetAccession(struct cplan9_s *hmm, char *acc)
{
  if (hmm->acc != NULL) free(hmm->acc);
  hmm->acc = Strdup(acc);
  StringChop(hmm->acc);
  hmm->flags |= CPLAN9_ACC;
}

/* Function: CPlan9SetDescription()
 * 
 * Purpose:  Change the description line of a Cplan9 HMM. Convenience function.
 * 
 * Note:     Trailing whitespace and \n's are chopped.
 */
void
CPlan9SetDescription(struct cplan9_s *hmm, char *desc)
{
  if (hmm->desc != NULL) free(hmm->desc);
  hmm->desc = Strdup(desc);
  StringChop(hmm->desc); 
  hmm->flags |= CPLAN9_DESC;
}

/* Function: CPlan9ComlogAppend()
 * Date:     SRE, Wed Oct 29 09:57:30 1997 [TWA 721 over Greenland] 
 * 
 * Purpose:  Concatenate command line options and append to the
 *           command line log.
 */
void
CPlan9ComlogAppend(struct cplan9_s *hmm, int argc, char **argv)
{
  int len;
  int i;

  /* figure out length of command line, w/ spaces and \n */
  len = argc;
  for (i = 0; i < argc; i++)
    len += strlen(argv[i]);

  /* allocate */
  if (hmm->comlog != NULL)
    {
      len += strlen(hmm->comlog);
      hmm->comlog = ReallocOrDie(hmm->comlog, sizeof(char)* (len+1));
    }
  else
    {
      hmm->comlog = MallocOrDie(sizeof(char)* (len+1));
      *(hmm->comlog) = '\0'; /* need this to make strcat work */
    }

  /* append */
  strcat(hmm->comlog, "\n");
  for (i = 0; i < argc; i++)
    {
      strcat(hmm->comlog, argv[i]);
      if (i < argc-1) strcat(hmm->comlog, " ");
    }
}

/* Function: CPlan9SetCtime()
 * Date:     SRE, Wed Oct 29 11:53:19 1997 [TWA 721 over the Atlantic]
 * 
 * Purpose:  Set the ctime field in a new HMM to the current time.
 */
void
CPlan9SetCtime(struct cplan9_s *hmm)
{
  time_t date = time(NULL);
  if (hmm->ctime != NULL) free(hmm->ctime);
  hmm->ctime = Strdup(ctime(&date));
  StringChop(hmm->ctime);
}


/* Function: CPlan9SetNullModel()
 * 
 * Purpose:  Set the null model section of an HMM.
 *           Convenience function.
 */
void
CPlan9SetNullModel(struct cplan9_s *hmm, float null[HMMERMAXABET], float p1)
{
  int x;
  for (x = 0; x < Alphabet_size; x++)
    hmm->null[x] = null[x];
  hmm->p1 = p1;
}


/* Function: CP9Logoddsify()
 * 
 * Purpose:  Take an HMM with valid probabilities, and
 *           fill in the integer log-odds score section of the model.
 *           
 *    Notes on log-odds scores (simplified from plan7.c):
 *         type of parameter       probability        score
 *         -----------------       -----------        ------
 *         any emission             p_x           log_2 p_x/null_x
 *         any transition           t_x           log_2 t_x
 *             
 * Args:      hmm          - the hmm to calculate scores in.
 *                  
 * Return:    (void)
 *            hmm scores are filled in.
 */  
void
CP9Logoddsify(struct cplan9_s *hmm)
{
  int k;			/* counter for model position */
  int x;			/* counter for symbols        */

  /*float accum;
    float tbm, tme;
  */
  if (hmm->flags & CPLAN9_HASBITS) return;

  /* Symbol emission scores
   */

  /* emission scores from state N, ins[0] */
  for (x = 0; x < Alphabet_size; x++) 
    hmm->isc[x][0] =  Prob2Score(hmm->ins[0][x], hmm->null[x]); 
  for (x = Alphabet_size; x < Alphabet_iupac; x++) 
    hmm->isc[x][0] = DegenerateSymbolScore(hmm->ins[0], hmm->null, x);

  for (k = 1; k <= hmm->M; k++) 
    {
      for (x = 0; x < Alphabet_size; x++) 
	{
	  hmm->msc[x][k] = Prob2Score(hmm->mat[k][x], hmm->null[x]);
	  hmm->isc[x][k] =  Prob2Score(hmm->ins[k][x], hmm->null[x]); 
	}
      /* degenerate match/insert emissions */
      for (x = Alphabet_size; x < Alphabet_iupac; x++) 
	{
	  hmm->msc[x][k] = DegenerateSymbolScore(hmm->mat[k], hmm->null, x);
	  hmm->isc[x][k] = DegenerateSymbolScore(hmm->ins[k], hmm->null, x);
	}
    }

  for (k = 0; k <= hmm->M; k++)
    {
      hmm->tsc[CTMM][k] = Prob2Score(hmm->t[k][CTMM], 1.0);
      hmm->tsc[CTMI][k] = Prob2Score(hmm->t[k][CTMI], 1.0);
      hmm->tsc[CTMD][k] = Prob2Score(hmm->t[k][CTMD], 1.0);
      hmm->tsc[CTIM][k] = Prob2Score(hmm->t[k][CTIM], 1.0);
      hmm->tsc[CTII][k] = Prob2Score(hmm->t[k][CTII], 1.0);
      hmm->tsc[CTID][k] = Prob2Score(hmm->t[k][CTID], 1.0);
      if(k != 0)
	{
	  hmm->tsc[CTDM][k] = Prob2Score(hmm->t[k][CTDM], 1.0);
	  hmm->tsc[CTDI][k] = Prob2Score(hmm->t[k][CTDI], 1.0);
	  hmm->tsc[CTDD][k] = Prob2Score(hmm->t[k][CTDD], 1.0);
	}
      else
	{
	  hmm->tsc[CTDM][k] = -INFTY;
	  hmm->tsc[CTDD][k] = -INFTY; /*D_0 doesn't exist*/
	  hmm->tsc[CTDI][k] = -INFTY;
	}
      if(k != 0)
	{
	  hmm->bsc[k]   = Prob2Score(hmm->begin[k], 1.0);
	  hmm->esc[k]   = Prob2Score(hmm->end[k], 1.0);
	}
    }
  hmm->flags |= CPLAN9_HASBITS;	/* raise the log-odds ready flag */
}


/* Function:  CPlan9Rescale()
 * Incept:    Steve Johnson, 3 May 2004
 *            eweights code incorp: SRE, Thu May 20 10:34:03 2004 [St. Louis]
 *
 * Purpose:   Scale a counts-based HMM by some factor, for
 *            adjusting counts to a new effective sequence number.
 *
 * Args:      hmm        - counts based HMM.
 *            scale      - scaling factor (e.g. eff_nseq/nseq); 1.0= no scaling.
 *
 * Returns:   (void)
 */
void 
CPlan9Rescale(struct cplan9_s *hmm, float scale)
{
  int k;

  /* emissions and transitions in the main model.
   * Note that match states are 1..M, insert states are 1..M-1,
   * and only nodes 1..M-1 have a valid array of transitions.
   */
  for(k = 1; k <= hmm->M; k++) 
    FScale(hmm->mat[k], Alphabet_size, scale);
  for(k = 0; k <=  hmm->M; k++) 
    FScale(hmm->ins[k], Alphabet_size, scale);
  for(k = 0; k <  hmm->M; k++) 
    FScale(hmm->t[k],   9,             scale);

  /* begin, end transitions; only valid [1..M] */
  FScale(hmm->begin+1, hmm->M, scale);
  FScale(hmm->end+1,   hmm->M, scale);
  
  return;
}


/* Function: CPlan9Renormalize()
 * 
 * Purpose:  Take an HMM in counts form, and renormalize
 *           all of its probability vectors. Also enforces
 *           CM Plan9 restrictions on nonexistent transitions.
 *           
 * Args:     hmm - the model to renormalize.
 *                 
 * Return:   (void)
 *           hmm is changed.
 */                          
void
CPlan9Renormalize(struct cplan9_s *hmm)
{
  int   k;			/* counter for model position */
  float d;			/* denominator */

				/* match emissions */
  FSet(hmm->mat[0], Alphabet_size, 0.);   /*M_0 is B state, non-emitter*/
  for (k = 1; k <= hmm->M; k++) 
    FNorm(hmm->mat[k], Alphabet_size);
				/* insert emissions */
  for (k = 0; k <= hmm->M; k++)
    FNorm(hmm->ins[k], Alphabet_size);

				/* begin transitions */
  d = FSum(hmm->begin+1, hmm->M) + hmm->t[0][CTMI] + hmm->t[0][CTMD];
  FScale(hmm->begin+1, hmm->M, 1./d);
  hmm->t[0][CTMI] /= d;
  hmm->t[0][CTMD] /= d;

  FNorm(hmm->t[0]+3, 3);	        /* transitions out of insert for node 0 (state N)*/
  FSet( hmm->t[0]+6, 3, 0.);    
				/* main model transitions */
  for (k = 1; k <= hmm->M; k++) /* safe for node M too, hmm->t[hmm->M][CTMM] should be 0.*/
    {
      d = FSum(hmm->t[k], 3) + hmm->end[k]; 
      FScale(hmm->t[k], 3, 1./d);
      hmm->end[k] /= d;

      FNorm(hmm->t[k]+3, 3);	/* insert */
      FNorm(hmm->t[k]+6, 3);	/* delete */
    }
				/* null model emissions */
  FNorm(hmm->null, Alphabet_size);

				/* special transitions, none?*/

  hmm->flags &= ~CPLAN9_HASBITS;	/* clear the log-odds ready flag */
  hmm->flags |= CPLAN9_HASPROB;	/* set the probabilities OK flag */
}

/* NO CPlan9RenormalizeExits(). 
   I can't see a need for it yet.
*/

/* NO CPlan9 configuration functions.
 */


/* Function: AllocCPlan9Matrix()
 * based on  AllocPlan7Matrix() <-- this function's comments below 
 * Date:     SRE, Tue Nov 19 07:14:47 2002 [St. Louis]
 *
 * Purpose:  Used to be the main allocator for dp matrices; we used to
 *           allocate, calculate, free. But this spent a lot of time
 *           in malloc(). Replaced with Create..() and Resize..() to
 *           allow matrix reuse in P7Viterbi(), the main alignment 
 *           engine. But matrices are alloc'ed by other alignment engines
 *           too, ones that are less frequently called and less 
 *           important to optimization of cpu performance. Instead of
 *           tracking changes through them, for now, provide
 *           an Alloc...() call with the same API that's just a wrapper.
 *
 * Args:     rows  - generally L+1, or 2; # of DP rows in seq dimension to alloc
 *           M     - size of model, in nodes
 *           xmx, mmx, imx, dmx 
 *                 - RETURN: ptrs to four mx components as a convenience
 *
 * Returns:  mx
 *           Caller free's w/ FreeCPlan9Matrix()
 */
struct cp9_dpmatrix_s *
AllocCPlan9Matrix(int rows, int M, int ***mmx, int ***imx, int ***dmx, int ***emx)
{
  struct cp9_dpmatrix_s *mx;
  mx = CreateCPlan9Matrix(rows-1, M, 0, 0);
  if (mmx != NULL) *mmx = mx->mmx;
  if (imx != NULL) *imx = mx->imx;
  if (dmx != NULL) *dmx = mx->dmx;
  if (emx != NULL) *emx = mx->emx;
  return mx;
}

/* Function: SizeCPlan9Matrix()
 * Date:     EPN, 04.21.06
 *
 * Purpose:  Return the size of a cp9_dpmatrix_s data structure with
 *           specified dimensions.
 *
 * Args:     rows  - generally L+1, or 2; # of DP rows in seq dimension to alloc
 *           M     - size of model, in nodes
 *
 * Returns:  size of cp9_dpmatrix_s in megabytes
 */
float
SizeCPlan9Matrix(int rows, int M)
{
  float ram;

  ram  = (float) (sizeof(struct cp9_dpmatrix_s));
  ram += (float) (sizeof(int *) * (rows+1) * 4);         /* mx->*mx */
  ram += (float) (sizeof(int)   * (rows+1) * (M+2) * 3); /* mx->*mx_mem */
  ram += (float) (sizeof(int)   * (rows+1));             /* mx->emx_mem */
  return (ram / 1000000.);
}


/* Function: FreeCPlan9Matrix()
 * based on  FreePlan7Matrix() <-- this function's comments below  
 * Purpose:  Free a dynamic programming matrix allocated by CreatePlan7Matrix().
 * 
 * Return:   (void)
 */
void
FreeCPlan9Matrix(struct cp9_dpmatrix_s *mx)
{
  free (mx->mmx_mem);
  free (mx->imx_mem);
  free (mx->dmx_mem);
  free (mx->emx_mem);
  free (mx->mmx);
  free (mx->imx);
  free (mx->dmx);
  free (mx->emx);
  free (mx);
}

#ifndef ALTIVEC
/* Function: CreateCPlan9Matrix()
 * based on  CreatePlan7Matrix() <-- this function's comments below  
 * Purpose:  Create a dynamic programming matrix for standard Forward,
 *           Backward, or Viterbi, with scores kept as scaled log-odds
 *           integers. Keeps 2D arrays compact in RAM in an attempt 
 *           to maximize cache hits. 
 *           
 *           The mx structure can be dynamically grown, if a new
 *           HMM or seq exceeds the currently allocated size. Dynamic
 *           growing is more efficient than an alloc/free of a whole
 *           matrix for every new target. The ResizePlan7Matrix()
 *           call does this reallocation, if needed. Here, in the
 *           creation step, we set up some pads - to inform the resizing
 *           call how much to overallocate when it realloc's. 
 *           
 * Args:     N     - N+1 rows are allocated, for sequence.  
 *           M     - size of model in nodes
 *           padN  - over-realloc in seq/row dimension, or 0
 *           padM  - over-realloc in HMM/column dimension, or 0
 *                 
 * Return:   mx
 *           mx is allocated here. Caller frees with FreeCPlan9Matrix(mx).
 */
struct cp9_dpmatrix_s *
CreateCPlan9Matrix(int N, int M, int padN, int padM)
{
  struct cp9_dpmatrix_s *mx;
  int i;

  mx          = (struct cp9_dpmatrix_s *) MallocOrDie (sizeof(struct cp9_dpmatrix_s));
  mx->mmx     = (int **) MallocOrDie (sizeof(int *) * (N+1));
  mx->imx     = (int **) MallocOrDie (sizeof(int *) * (N+1));
  mx->dmx     = (int **) MallocOrDie (sizeof(int *) * (N+1));
  mx->emx     = (int **) MallocOrDie (sizeof(int *) * (N+1));
  mx->mmx_mem = (void *) MallocOrDie (sizeof(int) * ((N+1)*(M+2)));
  mx->imx_mem = (void *) MallocOrDie (sizeof(int) * ((N+1)*(M+2)));
  mx->dmx_mem = (void *) MallocOrDie (sizeof(int) * ((N+1)*(M+2)));
  mx->emx_mem = (void *) MallocOrDie (sizeof(int) * (N+1) * (1));

  /* The indirect assignment below looks wasteful; it's actually
   * used for aligning data on 16-byte boundaries as a cache 
   * optimization in the fast altivec implementation
   */
  mx->mmx[0] = (int *) mx->mmx_mem;
  mx->imx[0] = (int *) mx->imx_mem;
  mx->dmx[0] = (int *) mx->dmx_mem;
  mx->emx[0] = (int *) mx->emx_mem;
  for (i = 1; i <= N; i++)
    {
      mx->mmx[i] = mx->mmx[0] + (i*(M+2));
      mx->imx[i] = mx->imx[0] + (i*(M+2));
      mx->dmx[i] = mx->dmx[0] + (i*(M+2));
      mx->emx[i] = mx->emx[0] + (i*1);
    }

  mx->maxN = N;
  mx->maxM = M;
  mx->padN = padN;
  mx->padM = padM;
  
  return mx;
}
#endif /*ALTIVEC*/

#ifndef ALTIVEC
/* Function: ResizeCPlan9Matrix()
 * based on  ResizePlan7Matrix() <-- this function's comments below  
 * Purpose:  Reallocate a dynamic programming matrix, if necessary,
 *           for a problem of NxM: sequence length N, model size M.
 *           (N=1 for small memory score-only variants; we allocate
 *           N+1 rows in the DP matrix.) 
 *           
 *           We know (because of the way hmmsearch and hmmpfam are coded)
 *           that only one of the two dimensions is going to change
 *           in size after the first call to ResizePlan7Matrix();
 *           that is, for hmmsearch, we have one HMM of fixed size M
 *           and our target sequences may grow in N; for hmmpfam,
 *           we have one sequence of fixed size N and our target models
 *           may grow in M. What we have to watch out for is P7SmallViterbi()
 *           working on a divide and conquer problem and passing us N < maxN,
 *           M > maxM; we should definitely *not* reallocate a smaller N.
 *           Since we know that only one dimension is going to grow,
 *           we aren't scared of reallocating to maxN,maxM. (If both
 *           M and N could grow, we would be more worried.)
 *
 *           Returns individual ptrs to the four matrix components
 *           as a convenience.
 *           
 * Args:     mx    - an already allocated model to grow.
 *           N     - seq length to allocate for; N+1 rows
 *           M     - size of model
 *           xmx, mmx, imx, dmx 
 *                 - RETURN: ptrs to four mx components as a convenience
 *                   
 * Return:   (void)
 *           mx is (re)allocated here.
 */
void
ResizeCPlan9Matrix(struct cp9_dpmatrix_s *mx, int N, int M, 
		   int ***mmx, int ***imx, int ***dmx, int ***emx)
{
  int i;
  printf("N: %d | maxN: %d | M: %d | maxM: %d\n", N, mx->maxN, M, mx->maxM);


  if (N <= mx->maxN && M <= mx->maxM) goto DONE;
  
  if (N > mx->maxN) {
    N          += mx->padN; 
    mx->maxN    = N; 
    mx->mmx     = (int **) ReallocOrDie (mx->mmx, sizeof(int *) * (mx->maxN+1));
    mx->imx     = (int **) ReallocOrDie (mx->imx, sizeof(int *) * (mx->maxN+1));
    mx->dmx     = (int **) ReallocOrDie (mx->dmx, sizeof(int *) * (mx->maxN+1));
    mx->emx     = (int **) ReallocOrDie (mx->emx, sizeof(int *) * (mx->maxN+1));
  }

  if (M > mx->maxM) {
    M += mx->padM; 
    mx->maxM = M; 
  }

  mx->mmx_mem = (void *) ReallocOrDie (mx->mmx_mem, sizeof(int) * ((mx->maxN+1)*(mx->maxM+2)));
  mx->imx_mem = (void *) ReallocOrDie (mx->imx_mem, sizeof(int) * ((mx->maxN+1)*(mx->maxM+2)));
  mx->dmx_mem = (void *) ReallocOrDie (mx->dmx_mem, sizeof(int) * ((mx->maxN+1)*(mx->maxM+2)));
  mx->emx_mem = (void *) ReallocOrDie (mx->emx_mem, sizeof(int) * ((mx->maxN+1)*(1)));

  mx->mmx[0] = (int *) mx->mmx_mem;
  mx->imx[0] = (int *) mx->imx_mem;
  mx->dmx[0] = (int *) mx->dmx_mem;
  mx->emx[0] = (int *) mx->emx_mem;

  for (i = 1; i <= mx->maxN; i++)
    {
      mx->mmx[i] = mx->mmx[0] + (i*(mx->maxM+2));
      mx->imx[i] = mx->imx[0] + (i*(mx->maxM+2));
      mx->dmx[i] = mx->dmx[0] + (i*(mx->maxM+2));
      mx->emx[i] = mx->emx[0] + (i*(1));
    }

 DONE:
  if (mmx != NULL) *mmx = mx->mmx;
  if (imx != NULL) *imx = mx->imx;
  if (dmx != NULL) *dmx = mx->dmx;
  if (emx != NULL) *emx = mx->emx;
}
#endif /*ALTIVEC*/




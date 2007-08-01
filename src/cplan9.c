/************************************************************
 * @LICENSE@
 ************************************************************/
/* cplan9.c based on plan7.c
 * EPN 02.27.06
 * 
 * Support for CM-Plan 9 HMM data structure, cplan9_s.
 * Including support for a dp matrix structure for cplan_s,
 * cp9_dpmatrix_s.
 * 
 * All of the CP9 code is based on analogous Plan 7 HMM code. There
 * are aspects of the plan 7 HMM data structure that I've kept in the
 * CM Plan 9 structre, but that are not yet used.
 * 
 * Included in this file are functions for configuring HMMs that were
 * built for 'sub CMs'.
 * 
 * At the end of this file are some functions that were stolen from
 * HMMER 2.4 and placed here without modification.
 */

#include "squidconf.h"
#include "cplan9.h"

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include "squid.h"
#include "funcs.h"
#include "structs.h"

static void rightjustify(char *s, int n);

/* Functions: AllocCPlan9(), AllocCPlan9Shell(), AllocCPlan9Body(), FreeCPlan9()
 * 
 * Purpose:   Allocate or free a CPlan9 HMM structure.
 *            Can either allocate all at once (AllocCPlan9()) or
 *            in two steps (AllocCPlan9Shell(), AllocCPlan9Body()).
 *            The two step method is used in CP9_hmmio.c where we start
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

  hmm->t      = NULL;
  hmm->mat    = NULL;
  hmm->ins    = NULL;
  
  hmm->tsc     = hmm->msc     = hmm->isc     = NULL;
  hmm->tsc_mem = hmm->msc_mem = hmm->msc_mem = NULL;

  hmm->begin  = NULL;
  hmm->end    = NULL;

  hmm->bsc = hmm->bsc_mem = NULL;
  hmm->esc = hmm->esc_mem = NULL;

  hmm->has_el      = NULL;
  hmm->el_from_ct  = NULL;
  hmm->el_from_idx = NULL;
  hmm->el_from_cmnd= NULL;

  hmm->flags = 0;
  return hmm;
}  

void
AllocCPlan9Body(struct cplan9_s *hmm, int M) 
{
  int k, x;

  hmm->M = M;

  hmm->t      = MallocOrDie ((M+1) *           sizeof(float *));
  hmm->mat    = MallocOrDie ((M+1) *           sizeof(float *));
  hmm->ins    = MallocOrDie ((M+1) *           sizeof(float *));
  hmm->t[0]   = MallocOrDie ((10*(M+1))     *  sizeof(float));
  hmm->mat[0] = MallocOrDie ((MAXABET*(M+1)) * sizeof(float));
  hmm->ins[0] = MallocOrDie ((MAXABET*(M+1)) *     sizeof(float));

  hmm->tsc     = MallocOrDie (10     *           sizeof(int *));
  hmm->msc     = MallocOrDie (MAXDEGEN   *       sizeof(int *));
  hmm->isc     = MallocOrDie (MAXDEGEN   *       sizeof(int *)); 
  hmm->tsc_mem = MallocOrDie ((10*(M+1))     *       sizeof(int));
  hmm->msc_mem = MallocOrDie ((MAXDEGEN*(M+1)) * sizeof(int));
  hmm->isc_mem = MallocOrDie ((MAXDEGEN*(M+1)) *     sizeof(int));

  hmm->tsc[0] = hmm->tsc_mem;
  hmm->msc[0] = hmm->msc_mem;
  hmm->isc[0] = hmm->isc_mem;

  /* note allocation strategy for important 2D arrays -- trying
   * to keep locality as much as possible, cache efficiency etc.
   */
  for (k = 1; k <= M; k++) {
    hmm->mat[k] = hmm->mat[0] + k * MAXABET;
    hmm->ins[k] = hmm->ins[0] + k * MAXABET;
    hmm->t[k]   = hmm->t[0]   + k * 10;
  }
  for (x = 1; x < MAXDEGEN; x++) {
    hmm->msc[x] = hmm->msc[0] + x * (M+1);
    hmm->isc[x] = hmm->isc[0] + x * (M+1);
  }
  for (x = 0; x < 10; x++)
    hmm->tsc[x] = hmm->tsc[0] + x * (M+1);

  /* tsc[x][0] is used as a boundary condition sometimes [Viterbi()],
   * so set to -inf always.
   */
  for (x = 0; x < 10; x++)
    hmm->tsc[x][0] = -INFTY;

  hmm->begin  = MallocOrDie  ((M+1) * sizeof(float));
  hmm->end    = MallocOrDie  ((M+1) * sizeof(float));

  hmm->bsc_mem  = MallocOrDie  ((M+1) * sizeof(int));
  hmm->esc_mem  = MallocOrDie  ((M+1) * sizeof(int));

  hmm->bsc = hmm->bsc_mem;
  hmm->esc = hmm->esc_mem;

  hmm->has_el      = MallocOrDie  ((M+1) * sizeof(int));
  hmm->el_from_ct  = MallocOrDie  ((M+2) * sizeof(int));
  hmm->el_from_idx = MallocOrDie  ((M+2) * sizeof(int *));
  hmm->el_from_cmnd = MallocOrDie  ((M+2) * sizeof(int *));

  return;
}  


void
FreeCPlan9(CP9_t *hmm)
{
  int k;
  if (hmm->bsc_mem    != NULL) free(hmm->bsc_mem);
  if (hmm->begin      != NULL) free(hmm->begin);
  if (hmm->esc_mem    != NULL) free(hmm->esc_mem);
  if (hmm->end        != NULL) free(hmm->end);
  if (hmm->msc_mem    != NULL) free(hmm->msc_mem);
  if (hmm->isc_mem    != NULL) free(hmm->isc_mem);
  if (hmm->tsc_mem    != NULL) free(hmm->tsc_mem);
  if (hmm->mat        != NULL) free(hmm->mat[0]);
  if (hmm->ins        != NULL) free(hmm->ins[0]);
  if (hmm->t          != NULL) free(hmm->t[0]);
  if (hmm->msc        != NULL) free(hmm->msc);
  if (hmm->isc        != NULL) free(hmm->isc);
  if (hmm->tsc        != NULL) free(hmm->tsc);
  if (hmm->mat        != NULL) free(hmm->mat);
  if (hmm->ins        != NULL) free(hmm->ins);
  if (hmm->t          != NULL) free(hmm->t);
  if (hmm->has_el     != NULL) free(hmm->has_el);
  if (hmm->el_from_ct != NULL) free(hmm->el_from_ct);
  if(hmm->el_from_idx != NULL)
    {
      for(k = 0; k <= hmm->M+1; k++)
	if(hmm->el_from_idx[k] != NULL)
	  free(hmm->el_from_idx[k]);
      free(hmm->el_from_idx);
    }
  if(hmm->el_from_cmnd != NULL)
    {
      for(k = 0; k <= hmm->M+1; k++)
	if(hmm->el_from_cmnd[k] != NULL)
	  free(hmm->el_from_cmnd[k]);
      free(hmm->el_from_cmnd);
    }

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
  FSet(hmm->t[0], 10, 0.);
  for (k = 1; k <= hmm->M; k++)
    {
      FSet(hmm->t[k], 10, 0.);
      FSet(hmm->mat[k], Alphabet_size, 0.);
      FSet(hmm->ins[k], Alphabet_size, 0.);
    }
  FSet(hmm->begin+1, hmm->M, 0.);
  FSet(hmm->end+1, hmm->M, 0.);

  /* initialize the el_* data structures, these
   * depend on the CM guide tree and will be set
   * when the CP9 is constructed from the CM.
   */
  for (k = 0; k <= (hmm->M); k++)
    {
      hmm->has_el[k]      = FALSE;         
      hmm->el_from_ct[k]  = 0;
      hmm->el_from_idx[k] = NULL; 
      hmm->el_from_cmnd[k] = NULL; 
    }
  /* special case hmm->M+1 corresponds to the E state here */
  hmm->el_from_ct[(hmm->M+1)]  = 0;
  hmm->el_from_idx[(hmm->M+1)] = NULL; 
  hmm->el_from_cmnd[(hmm->M+1)] = NULL; 

  hmm->flags &= ~CPLAN9_HASBITS;	/* invalidates scores */
  hmm->flags &= ~CPLAN9_HASPROB;	/* invalidates probabilities */
  hmm->el_self = 0.; /* EL self transition probability */
}


/* Function: CPlan9SetNullModel()
 * 
 * Purpose:  Set the null model section of an HMM.
 *           Convenience function.
 */
void
CPlan9SetNullModel(struct cplan9_s *hmm, float null[MAXABET], float p1)
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
	  hmm->msc[x][k] =  Prob2Score(hmm->mat[k][x], hmm->null[x]);
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
      hmm->tsc[CTME][k] = Prob2Score(hmm->t[k][CTME], 1.0);
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
  hmm->el_selfsc = Prob2Score(hmm->el_self, 1.0);
  hmm->flags |= CPLAN9_HASBITS;	/* raise the log-odds ready flag */
}


/* Function:  CPlan9Rescale() 
 *            EPN based on Steve Johnsons plan 7 version
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
    FScale(hmm->t[k],   10,             scale);

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
  d = FSum(hmm->begin+1, hmm->M) + hmm->t[0][CTMI] + hmm->t[0][CTMD] + hmm->t[0][CTME]; 
  /* hmm->t[0][CTME] should always be 0., can't local end from the M_0 == B state */
  FScale(hmm->begin+1, hmm->M, 1./d);
  hmm->t[0][CTMI] /= d;
  hmm->t[0][CTMD] /= d;
  hmm->t[0][CTME] /= d;

  FNorm(hmm->t[0]+4, 3);	        /* transitions out of insert for node 0 (state N)*/
  FSet( hmm->t[0]+7, 3, 0.);    
				/* main model transitions */
  for (k = 1; k <= hmm->M; k++) /* safe for node M too, hmm->t[hmm->M][CTMM] should be 0.*/
    {
      d = FSum(hmm->t[k], 4) + hmm->end[k]; 
      FScale(hmm->t[k], 4, 1./d);
      hmm->end[k] /= d;

      FNorm(hmm->t[k]+4, 3);	/* insert */
      FNorm(hmm->t[k]+7, 3);	/* delete */
    }
				/* null model emissions */
  FNorm(hmm->null, Alphabet_size);

				/* special transitions, none?*/

  hmm->flags &= ~CPLAN9_HASBITS;	/* clear the log-odds ready flag */
  hmm->flags |= CPLAN9_HASPROB;	/* set the probabilities OK flag */
}

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
AllocCPlan9Matrix(int rows, int M, int ***mmx, int ***imx, int ***dmx, int ***elmx, int **erow)
{
  struct cp9_dpmatrix_s *mx;
  mx = CreateCPlan9Matrix(rows-1, M, 0, 0);
  if (mmx != NULL) *mmx = mx->mmx;
  if (imx != NULL) *imx = mx->imx;
  if (dmx != NULL) *dmx = mx->dmx;
  if (elmx!= NULL) *elmx= mx->elmx;
  if (erow!= NULL) *erow= mx->erow;
  
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
  ram += (float) (sizeof(int *) * (rows+1) * 5);         /* mx->*mx */
  ram += (float) (sizeof(int)   * (rows+1) * (M+1) * 4); /* mx->*mx_mem */
  ram += (float) (sizeof(int)   * (rows+1));             /* mx->erow */
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
  free (mx->elmx_mem);
  free (mx->mmx);
  free (mx->imx);
  free (mx->dmx);
  free (mx->elmx);
  free (mx->erow);
  free (mx);
}

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
  mx->elmx    = (int **) MallocOrDie (sizeof(int *) * (N+1)); 
  /* slightly wasteful, some nodes can't go to EL (ex: right half of MATPs) */
  mx->erow    = (int *) MallocOrDie (sizeof(int) * (N+1));
  mx->mmx_mem = (void *) MallocOrDie (sizeof(int) * ((N+1)*(M+1)));
  mx->imx_mem = (void *) MallocOrDie (sizeof(int) * ((N+1)*(M+1)));
  mx->dmx_mem = (void *) MallocOrDie (sizeof(int) * ((N+1)*(M+1)));
  mx->elmx_mem= (void *) MallocOrDie (sizeof(int) * ((N+1)*(M+1)));

  /* The indirect assignment below looks wasteful; it's actually
   * used for aligning data on 16-byte boundaries as a cache 
   * optimization in the fast altivec implementation
   */
  mx->mmx[0] = (int *) mx->mmx_mem;
  mx->imx[0] = (int *) mx->imx_mem;
  mx->dmx[0] = (int *) mx->dmx_mem;
  mx->elmx[0]= (int *) mx->elmx_mem;
  for (i = 1; i <= N; i++)
    {
      mx->mmx[i] = mx->mmx[0] + (i*(M+1));
      mx->imx[i] = mx->imx[0] + (i*(M+1));
      mx->dmx[i] = mx->dmx[0] + (i*(M+1));
      mx->elmx[i]= mx->elmx[0]+ (i*(M+1));
    }

  mx->maxN = N;
  mx->maxM = M;
  mx->padN = padN;
  mx->padM = padM;

  return mx;
}

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
 *           xmx, mmx, imx, dmx, elmx, erow 
 *                 - RETURN: ptrs to four mx components as a convenience
 *                   
 * Return:   (void)
 *           mx is (re)allocated here.
 */
void
ResizeCPlan9Matrix(struct cp9_dpmatrix_s *mx, int N, int M, 
		   int ***mmx, int ***imx, int ***dmx, int ***elmx, int **erow)
{
  int i;
  /*printf("N: %d | maxN: %d | M: %d | maxM: %d\n", N, mx->maxN, M, mx->maxM);*/

  if (N <= mx->maxN && M <= mx->maxM) goto DONE;
  
  if (N > mx->maxN) {
    N          += mx->padN; 
    mx->maxN    = N; 
    mx->mmx     = (int **) ReallocOrDie (mx->mmx, sizeof(int *) * (mx->maxN+1));
    mx->imx     = (int **) ReallocOrDie (mx->imx, sizeof(int *) * (mx->maxN+1));
    mx->dmx     = (int **) ReallocOrDie (mx->dmx, sizeof(int *) * (mx->maxN+1));
    mx->elmx    = (int **) ReallocOrDie (mx->elmx,sizeof(int *) * (mx->maxN+1));
    mx->erow    = (void *) ReallocOrDie (mx->erow,sizeof(int)   * (mx->maxN+1));
  }

  if (M > mx->maxM) {
    M += mx->padM; 
    mx->maxM = M; 
  }

  mx->mmx_mem = (void *) ReallocOrDie (mx->mmx_mem, sizeof(int) * ((mx->maxN+1)*(mx->maxM+1)));
  mx->imx_mem = (void *) ReallocOrDie (mx->imx_mem, sizeof(int) * ((mx->maxN+1)*(mx->maxM+1)));
  mx->dmx_mem = (void *) ReallocOrDie (mx->dmx_mem, sizeof(int) * ((mx->maxN+1)*(mx->maxM+1)));
  mx->elmx_mem= (void *) ReallocOrDie (mx->elmx_mem,sizeof(int) * ((mx->maxN+1)*(mx->maxM+1)));

  mx->mmx[0] = (int *) mx->mmx_mem;
  mx->imx[0] = (int *) mx->imx_mem;
  mx->dmx[0] = (int *) mx->dmx_mem;
  mx->elmx[0]= (int *) mx->elmx_mem;

  for (i = 1; i <= mx->maxN; i++)
    {
      mx->mmx[i] = mx->mmx[0] + (i*(mx->maxM+1));
      mx->imx[i] = mx->imx[0] + (i*(mx->maxM+1));
      mx->dmx[i] = mx->dmx[0] + (i*(mx->maxM+1));
      mx->elmx[i]= mx->elmx[0]+ (i*(mx->maxM+1));
    }

 DONE:
  if (mmx != NULL) *mmx = mx->mmx;
  if (imx != NULL) *imx = mx->imx;
  if (dmx != NULL) *dmx = mx->dmx;
  if (elmx!= NULL) *elmx= mx->elmx;
  if (erow != NULL) *erow = mx->erow;
}

/* Function: CPlan9SWConfig()
 * EPN 05.30.06
 * based on SRE's Plan7SWConfig() from HMMER's plan7.c
 * 
 * Purpose:  Set the alignment independent parameters of
 *           a CM Plan 9 model to hmmsw (Smith/Waterman) configuration.
 *           
 * Notes:    The desideratum for begin/end probs is that all fragments ij
 *           (starting at match i, ending at match j) are
 *           equiprobable -- there is no information in the choice of
 *           entry/exit. There are M(M+1)/2 possible choices of ij, so
 *           each must get a probability of 2/M(M+1). This prob is the
 *           product of a begin, an end, and all the not-end probs in
 *           the path between i,j. 
 *            
 *           Thus: entry/exit is asymmetric because of the left/right
 *           nature of the HMM/profile. Entry probability is distributed
 *           simply by assigning p_x = pentry / (M-1) to M-1 
 *           internal match states. However, the same approach doesn't
 *           lead to a flat distribution over exit points. Exit p's
 *           must be corrected for the probability of a previous exit
 *           from the model. Requiring a flat distribution over exit
 *           points leads to an easily solved piece of algebra, giving:
 *                      p_1 = pexit / (M-1)
 *                      p_x = p_1 / (1 - (x-1) p_1)
 *           
 * Args:     hmm    - the CM Plan 9 model w/ data-dep prob's valid
 *           pentry - probability of an internal entry somewhere;
 *                    will be evenly distributed over M-1 match states
 *           pexit  - probability of an internal exit somewhere; 
 *                    will be distributed over M-1 match states.
 *                    
 * Return:   (void)
 *           HMM probabilities are modified.
 */
void
CPlan9SWConfig(struct cplan9_s *hmm, float pentry, float pexit)
{
  float basep;			/* p1 for exits: the base p */
  int   k;			/* counter over states      */

  /* No special (*x* states in Plan 7) states in CM Plan 9 */

  /*for (k = 1; k <= hmm->M; k++)
    printf("before anything: end[%d]: %f\n", k, hmm->end[k]);*/
  /* Configure entry.
   */
  hmm->begin[1] = (1. - pentry) * (1. - (hmm->t[0][CTMI] + hmm->t[0][CTMD] + hmm->t[0][CTME]));
  FSet(hmm->begin+2, hmm->M-1, (pentry * (1.- (hmm->t[0][CTMI] + hmm->t[0][CTMD] + hmm->t[0][CTME]))) 
       / (float)(hmm->M-1));
  
  /* Configure exit.
   * Don't touch hmm->end[hmm->M]
   */

  basep = pexit / (float) (hmm->M-1);
  for (k = 1; k < hmm->M; k++)
    hmm->end[k] = basep / (1. - basep * (float) (k-1));
  CPlan9RenormalizeExits(hmm, 1);
  /*for (k = 1; k <= hmm->M; k++)
    printf("after renormalizing: end[%d]: %f\n", k, hmm->end[k]);*/

  hmm->flags       &= ~CPLAN9_HASBITS; /* reconfig invalidates log-odds scores */
  hmm->flags       |= CPLAN9_LOCAL_BEGIN; /* local begins now on */
  hmm->flags       |= CPLAN9_LOCAL_END;   /* local ends now on */

  CP9Logoddsify(hmm);
}

/* Function: CPlan9CMLocalBeginConfig()
 * Incept:   EPN, Thu Jun 21 15:43:29 2007
 * based on SRE's Plan7SWConfig() from HMMER's plan7.c
 * 
 * Purpose:  Set up a CM Plan 9 HMM to mimic CM local begins as closely
 *           as it can. We can't enforce that a begin/end point are chosen
 *           the same way a CM's are, as the choice of a CM local begin
 *           (in non-truncated CYK mode) defines both a start and end point,
 *           and some start/end combinations are impossible. For the CP9
 *           we allow all possible start/end combos.
 *           
 * Args:     cm    - the CM, must have valid cm->cp9, we'll use
 *                   the CM local begin probs to set the cm->cp9s
 *                   begin/end probs.
 *                    
 * Return:   (void)
 *           HMM probabilities are modified.
 */
void
CPlan9CMLocalBeginConfig(CM_t *cm)
{
  /* Contract checks */
  if(cm->cp9 == NULL)
    esl_fatal("ERROR in CPlan9CMLocalBeginConfig, cm->cp9 is NULL.\n");
  if(cm->cp9map == NULL)
    esl_fatal("ERROR in CPlan9CMLocalBeginConfig, cm->cp9map is NULL.\n");
  if(!(cm->flags & CM_CP9))
     esl_fatal("ERROR in CPlan9CMLocalBeginConfig, CM_CP9 flag is down.");
  if(!(cm->flags & CM_LOCAL_BEGIN))
     esl_fatal("ERROR in CPlan9CMLocalBeginConfig, CM_LOCAL_BEGIN flag is down.");
  if(!(cm->flags & CM_LOCAL_END))
     esl_fatal("ERROR in CPlan9CMLocalBeginConfig, CP9_LOCAL_BEGIN flag is already up.");
  if(cm->cp9->flags & CPLAN9_LOCAL_END)
     esl_fatal("ERROR in CPlan9CMLocalBeginConfig, CP9_LOCAL_END flag is already up.");

  /* Configure entry.
   * To match CM, we enforce the only way out of the B state (M_0)
   * is through a local begin into a match state 
   */
  esl_fatal("In CPlan9CMLocalBeginConfig(), function not yet finished.\n");

  /* To do: determine which nodes we can begin into (for example can't begin
   * into a node modeled by right half of a MATP). And those nodes we can
   * end out of (ex: can't end out of left half of MATP). I think we can't
   * get the local begin/end scores match the cm->begin scores exactly 
   * though b/c hmm->begin[1] should match cm->begin[cm->nodemap[1]], but
   * then the begin[2]-->[M] should only be half what they are in the 
   * CM (because the end[2]->end[M] will incur a penalty also.
   */

  float hmm_begin_end_prob;
  hmm_begin_end_prob = cm->begin[cm->nodemap[1]];
  /***************************************************/
  /* BELOW IS INCOMPLETE! */
  /*cm->cp9->t[0][CTMI] = cm->cp9->t[0][CTMD] = cm->cp9->t[0][CTME] = 0.;
  cm->cp9->begin[1] = hmm_begin_end_prob;
  FSet(cm->cp9->begin+2, cm->cp9->M-1, ((hmm_begin_end_prob/2.) / (float)(cm->cp9->M-1)));*/

  /* Configure hmm->ends, there is no equivalent in the CM, as cm->end's are local
   * ends (these correspond to EL states in the HMM see CPlan9ELConfig). A 
   * CM local begin defines a start and end, we can't enforce that here, what
   * we do is set end probabilities as 
   * 
   * Don't touch cm->cp9->end[cm->cp9->M]
   */

  /*basep = pexit / (float) (cm->cp9->M-1);
  for (k = 1; k < cm->cp9->M; k++)
    cm->cp9->end[k] = basep / (1. - basep * (float) (k-1));
    CPlan9RenormalizeExits(cm->cp9, 1);*/

  cm->cp9->flags       &= ~CPLAN9_HASBITS; /* reconfig invalidates log-odds scores */
  cm->cp9->flags       |= CPLAN9_LOCAL_BEGIN; /* local begins now on */
  cm->cp9->flags       |= CPLAN9_LOCAL_END;   /* local ends now on */

  CP9Logoddsify(cm->cp9);
}

/* Function: CPlan9ELConfig()
 * Incept:   EPN, Tue Jun 19 09:50:52 2007
 * 
 * Purpose:  Turn EL local ends in a CM Plan 9 HMM on based on 
 *           the local end probs in the CM. 
 *           
 * Args:     cm     - the CM, must have valid CP9 HMM
 *                    
 * Return:   (void)
 *           HMM probabilities are modified.
 */
void
CPlan9ELConfig(CM_t *cm)
{
  /*printf("IN CPlan9ELConfig\n");*/
  /* Contract checks */
  if(cm->cp9 == NULL)
    esl_fatal("ERROR in CPlan9ELConfig, cm->cp9 is NULL.\n");
  if(cm->cp9map == NULL)
    esl_fatal("ERROR in CPlan9ELConfig, cm->cp9map is NULL.\n");
  if(!(cm->flags & CM_CP9))
     esl_fatal("ERROR in CPlan9ELConfig, CM_CP9 flag is down.");
  if(!(cm->flags & CM_LOCAL_END))
     esl_fatal("ERROR in CPlan9ELConfig, CM_LOCAL_END flag is down.");
  if(cm->cp9->flags & CPLAN9_EL)
     esl_fatal("ERROR in CPlan9ELConfig, CP9_EL flag is already up.");
  
  int v;
  int k;                     /* counter over HMM nodes */
  int nd;
  int seen_exit;
  float to_el_prob;
  float norm_factor;

  /* Check to make sure all non-zero local end probabilities 
   * in the CM are identical (within reasonable precision), 
   * use that probability to set all HMM transitions to EL states.
   */
  seen_exit  = FALSE;
  to_el_prob = 0.;
  for(v = 0; v < cm->M; v++)
    {
      nd = cm->ndidx[v];
      if (((cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
	    cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BEGL_nd ||
	    cm->ndtype[nd] == BEGR_nd) && 
	   cm->ndtype[nd+1] != END_nd) && cm->nodemap[nd] == v)
	{
	  /* this should have a non-zero local end probability */
	  if(fabs(cm->end[v] - 0.) < 0.00001) /* non-zero */
	    esl_fatal("In CPlan9ELConfig(), CM state %d should have non-zero local end prob, but it doesn't.\n", v);
	  if(!seen_exit)	  
	    {
	      to_el_prob = cm->end[v];
	      seen_exit  = TRUE;
	    }
	  else if(fabs(to_el_prob - cm->end[v]) > 0.00001)
	    esl_fatal("In CPlan9ELConfig(), not all CM states EL probs are identical.\n");
	}
    } 

  /* transitions from HMM node 0 to EL is impossible */
  cm->cp9->t[0][CTME] = 0.;
  for(k = 1; k <= cm->cp9->M; k++) 
    {
      if(cm->cp9->has_el[k])
	{
	  cm->cp9->t[k][CTME] = to_el_prob;
	  norm_factor = 1. - (cm->cp9->t[k][CTME] / (1. - cm->cp9->end[k]));
	  cm->cp9->t[k][CTMM] *= norm_factor;
	  cm->cp9->t[k][CTMI] *= norm_factor;
	  cm->cp9->t[k][CTMD] *= norm_factor;
	  /* cm->cp9->end[k] untouched */
	}
    }
  cm->cp9->flags &= ~CPLAN9_HASBITS;	/* clear the log-odds ready flag */

  CP9Logoddsify(cm->cp9);

  cm->cp9->flags |= CPLAN9_EL;          /* EL end locals now on */
  /*debug_print_cp9_params(cm->cp9);*/
  return;
}

/* Function: CPlan9NoEL()
 * Incept:   EPN, Tue Jun 19 09:50:52 2007
 * 
 * Purpose:  Turn EL local ends off in a CM Plan 9 HMM
 *           
 * Args:     cm     - the CM, must have valid CP9 HMM
 *                    
 * Return:   (void)
 *           HMM probabilities are modified.
 */
void
CPlan9NoEL(CM_t *cm)
{
  /* Contract checks */
  if(cm->cp9 == NULL)
    esl_fatal("ERROR in CPlan9ELConfig, cm->cp9 is NULL.\n");
  if(cm->cp9map == NULL)
    esl_fatal("ERROR in CPlan9ELConfig, cm->cp9map is NULL.\n");
  if(!(cm->flags & CM_CP9))
     esl_fatal("ERROR in CPlan9ELConfig, CM_CP9 flag is down.");
  if(!(cm->cp9->flags & CPLAN9_EL))
     esl_fatal("ERROR in CPlan9ELConfig, CP9_EL flag is already down.");
  
  int k;                     /* counter over HMM nodes */

  for(k = 0; k <= cm->cp9->M; k++) 
    cm->cp9->t[k][CTME] = 0.;
  CPlan9RenormalizeExits(cm->cp9, 1);

  cm->cp9->flags &= ~CPLAN9_HASBITS;	/* clear the log-odds ready flag */
  CP9Logoddsify(cm->cp9);

  cm->cp9->flags &= ~CPLAN9_EL;          /* EL end locals now off */

  return;
}

/* Function: CPlan9InitEL()
 * Incept:   EPN, Tue Jun 19 13:10:56 2007
 * 
 * Purpose:  Initialize a CP9 HMM for possible EL local ends
 *           by determining how the EL states should be connected
 *           based on the CM node topology.
 *           
 * Args:     cm     - the CM
 *           cp9    - the CP9 HMM, built from cm
 *
 * Return:   (void)
 */
void
CPlan9InitEL(CM_t *cm, CP9_t *cp9)
{
  CMEmitMap_t *emap;         /* consensus emit map for the CM */
  int k;                     /* counter over HMM nodes */
  int nd;
  int *tmp_el_from_ct;

  /* First copy the CM el self transition score/probability: */
  cp9->el_self   = sreEXP2(cm->el_selfsc);
  cp9->el_selfsc = Prob2Score(cp9->el_self, 1.0);

  /* For each HMM node k, we can transit FROM >= 0 EL states from 
   * HMM nodes kp. Determine how many such valid transitions exist
   * from each node, then allocate and fill cp9->el_from_idx[k] and 
   * cp9->el_from_cmnd arrays based on that.
   * This two-pass method saves memory b/c we only allocate for
   * what we'll need.
   */
  emap = CreateEmitMap(cm); 

  /* Initialize to 0 */
  for(k = 0; k <= cp9->M; k++) 
    {
      cp9->el_from_ct[k] = 0;
      cp9->has_el[k] = FALSE;
    }
  cp9->el_from_ct[(cp9->M+1)] = 0; /* special case, we can get to E state from EL states */
    
  /* first pass to get number of valid transitions */
  for(nd = 0; nd < cm->nodes; nd++)
    {
      if ((cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
	   cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BEGL_nd ||
	   cm->ndtype[nd] == BEGR_nd) && 
	  cm->ndtype[nd+1] != END_nd)
	{
	  /*printf("HMM node %d can be reached from HMM node %d's EL state\n", emap->rpos[nd], emap->lpos[nd]);*/
	  cp9->el_from_ct[emap->rpos[nd]]++;
	  cp9->has_el[emap->lpos[nd]] = TRUE;
	}
    }

  /* allocate cp9->el_from_idx[k], cp9->el_from_cmnd for all k */
  for(k = 0; k <= (cp9->M+1); k++) 
    {
      if(cp9->el_from_idx[k] != NULL) /* if !NULL we already filled it, shouldn't happen */
	esl_fatal("ERROR in CPlan9InitEL() el_from_idx has already been initialized\n");
      if(cp9->el_from_ct[k] > 0)
	{
	  cp9->el_from_idx[k] = MallocOrDie(sizeof(int) * cp9->el_from_ct[k]);
	  cp9->el_from_cmnd[k] = MallocOrDie(sizeof(int) * cp9->el_from_ct[k]);
	}
      /* else it remains NULL */
    }

  /* now fill in cp9->el_from_idx, we need a new counter array */
  tmp_el_from_ct = MallocOrDie(sizeof(int) * (cp9->M+2));
  for(k = 0; k <= (cp9->M+1); k++) 
    tmp_el_from_ct[k] = 0;
  for(nd = 0; nd < cm->nodes; nd++)
    {
      if ((cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
	   cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BEGL_nd ||
	   cm->ndtype[nd] == BEGR_nd) && 
	  cm->ndtype[nd+1] != END_nd)
	{
	  k = emap->rpos[nd];
	  cp9->el_from_idx[k][tmp_el_from_ct[k]] = emap->lpos[nd];
	  cp9->el_from_cmnd[k][tmp_el_from_ct[k]] = nd;
	  tmp_el_from_ct[k]++;
	}
    }

  /* Debugging printfs */
  /*  for(k = 0; k <= (cp9->M+1); k++) 
    {
      for(c = 0; c < cp9->el_from_ct[k]; c++)
	printf("cp9->el_from_idx[%3d][%2d]: %4d\n", k, c, cp9->el_from_idx[k][c]);
      if(cp9->has_el[k])
      printf("node k:%3d HAS an EL!\n", k);
      }*/

  /* Free memory and exit */
  free(tmp_el_from_ct);
  FreeEmitMap(emap);
}

/* Function: CPlan9GlobalConfig()
 * EPN 09.24.06
 * based on SRE's Plan7GlobalConfig() from HMMER's plan7.c
 * 
 * Purpose:  Set the alignment independent parameters of
 *           a CM Plan 9 model to (Needleman/Wunsch) configuration.
 *           Make all transitions to EL states impossible.
 *           
 * Args:     hmm    - the CM Plan 9 model w/ data-dep prob's valid
 *           pentry - probability of an internal entry somewhere;
 *                    will be evenly distributed over M-1 match states
 *           pexit  - probability of an internal exit somewhere; 
 *                    will be distributed over M-1 match states.
 *                    
 * Return:   (void)
 *           HMM probabilities are modified.
 */
void
CPlan9GlobalConfig(struct cplan9_s *hmm)
{
  int k;
  /* No special (*x* states in Plan 7) states in CM Plan 9 */

  /* Configure entry.
   * Exactly 3 ways to start, B->M_1 (hmm->begin[1]), B->I_0 (hmm->t[0][CTMI]),
   *                      and B->D_1 (hmm->t[0][CTMD])
   */
  hmm->begin[1] = 1. - (hmm->t[0][CTMI] + hmm->t[0][CTMD] + hmm->t[0][CTME]); 
  /* this is okay, hmm->t[0] is never changed, even during local
   * configuration */
  FSet(hmm->begin+2, hmm->M-1, 0.);
  
  hmm->end[hmm->M] = 1. - hmm->t[hmm->M][CTMI];
  FSet(hmm->end+1, hmm->M-1, 0.);
  CPlan9RenormalizeExits(hmm, 1);

  /* Make all transitions to EL impossible, node 0, M special and should 
   * always have CTME transition as impossible. */
  for(k = 1; k < hmm->M; k++)
    {
      hmm->t[k][CTME] = 0.;
      FNorm(hmm->t[k], 4); /* renormalize transitions out of node k */
    }
  hmm->flags       &= ~CPLAN9_HASBITS; /* reconfig invalidates log-odds scores */
  hmm->flags       &= ~CPLAN9_LOCAL_BEGIN; /* local begins now off */
  hmm->flags       &= ~CPLAN9_LOCAL_END;   /* local ends now off */
  hmm->flags       &= ~CPLAN9_EL;          /* EL end locals now off */

  CP9Logoddsify(hmm);
}


/* Function: CPlan9ConfigEnforce()
 * EPN, Fri Feb  9 05:47:37 2007
 * based on SRE's Plan7SWConfig() from HMMER's plan7.c
 * 
 * Purpose:  Set the alignment independent parameters of
 *           a CM Plan 9 model to hmmsw (Smith/Waterman) configuration.
 *           Same as CPlan9SWConfig but enforces a contiguous subset of
 *           nodes start at x, ending at y must be entered by forbidding
 *           local entries after x and local exits before y.
 *           
 * Args:     hmm    - the CM Plan 9 model w/ data-dep prob's valid
 *           pentry - probability of an internal entry somewhere;
 *                    will be evenly distributed over (enf_start) 
 *                    match states
 *           pexit  - probability of an internal exit somewhere; 
 *                    will be distributed over (M-enf_end) match 
 *                    states.
 *           enf_start_pos - HMM node where enforced node subset begins         
 *           enf_end_pos   - HMM node where enforced node subset ends
 * Return:   (void)
 *           HMM probabilities are modified.
 */
void
CPlan9SWConfigEnforce(struct cplan9_s *hmm, float pentry, float pexit,
		      int enf_start_pos, int enf_end_pos)
{
  float basep;			/* p1 for exits: the base p */
  int   k;			/* counter over states      */

  /* No special (*x* states in Plan 7) states in CM Plan 9 */

  /* Configure entry.
   * To match CM, we enforce the only way out of the B state (M_0)
   * is through a local begin into a match state 
   */
  hmm->t[0][CTMI] = hmm->t[0][CTMD] = hmm->t[0][CTME] = 0.;
  hmm->begin[1] = 1. - pentry;
  for (k = 2; k <= enf_start_pos; k++)
    hmm->begin[k] = pentry / (float)(enf_start_pos-1);


  /* OLD WAY (more smith-waterman-like, less CM-like) EPN, Thu Jun 21 15:30:46 2007
     hmm->begin[1] = (1. - pentry) * (1. - (hmm->t[0][CTMI] + hmm->t[0][CTMD])); 
  for (k = 2; k <= enf_start_pos; k++)
    hmm->begin[k] = (pentry * (1.- (hmm->t[0][CTMI] + hmm->t[0][CTMD]))) / (float)(enf_start_pos-1);
  for (k = (enf_start_pos+1); k <= hmm->M; k++)
    hmm->begin[k] = 0.;
  */
    
  /* Configure exit.
   * Don't touch hmm->end[hmm->M]
   */
  if(enf_end_pos == hmm->M) /* no local exit possible */
    basep = 0.0;
  else
    basep = pexit / (float) (hmm->M-enf_end_pos);

  for (k = 0; k < enf_end_pos; k++)
    hmm->end[k] = 0.;
  for (k = enf_end_pos; k < hmm->M; k++)
    hmm->end[k] = basep / (1. - basep * (float) ((k-enf_end_pos)-1));
  CPlan9RenormalizeExits(hmm, 1);
  hmm->flags       &= ~CPLAN9_HASBITS;     /* reconfig invalidates log-odds scores */
  hmm->flags       |= CPLAN9_LOCAL_BEGIN; /* local begins now on */
  hmm->flags       |= CPLAN9_LOCAL_END;   /* local ends now on */
}

/* Function: CPlan9RenormalizeExits()
 * EPN 05.30.06 based on SRE's Plan7RenormalizeExits() from
 *                       HMMER's plan7.c.
 *
 * Date:     SRE, Fri Aug 14 11:22:19 1998 [St. Louis]
 *
 * Purpose:  Renormalize just the match state transitions;
 *           for instance, after a Config() function has
 *           modified the exit distribution.
 *
 * Args:     hmm - hmm to renormalize
 *           spos   - first consensus column modelled by original
 *                    CP9 HMM the sub CP9 HMM models. Often 1.
 * Returns:  void
 */
void
CPlan9RenormalizeExits(struct cplan9_s *hmm, int spos)
{
  int   k;
  float d;

  /* We can't exit from node 0 so we start renormalizing at node 1 */
  for (k = 1; k < hmm->M; k++)
    {
      if(k != (spos-1)) /* we can't exit from the M_spos-1 */
	{
	  d = FSum(hmm->t[k], 4);
	  /* FScale(hmm->t[k], 4, 1./(d + d*hmm->end[k])); */
	  FScale(hmm->t[k], 4, (1.-hmm->end[k])/d);
	}
    }
  /* Take care of hmm->M node, which is special */
  d = hmm->t[hmm->M][CTMI] + hmm->t[hmm->M][CTME]; /* CTMD is IMPOSSIBLE, CTMM is hmm->end[hmm-M] */
  hmm->t[hmm->M][CTMI] *= (1.-hmm->end[hmm->M])/d;
  hmm->t[hmm->M][CTME] *= (1.-hmm->end[hmm->M])/d;
  return;
}

/* Function: CP9AllocTrace(), CP9ReallocTrace(), CP9FreeTrace()
 * 
 * Purpose:  allocation and freeing of traceback structures
 */
void
CP9AllocTrace(int tlen, CP9trace_t **ret_tr)
{
  CP9trace_t *tr;
  
  tr =            MallocOrDie (sizeof(CP9trace_t));
  tr->statetype = MallocOrDie (sizeof(char) * tlen);
  tr->nodeidx   = MallocOrDie (sizeof(int)  * tlen);
  tr->pos       = MallocOrDie (sizeof(int)  * tlen);
  *ret_tr = tr;
}
void
CP9ReallocTrace(CP9trace_t *tr, int tlen)
{
  tr->statetype = ReallocOrDie (tr->statetype, tlen * sizeof(char));
  tr->nodeidx   = ReallocOrDie (tr->nodeidx,   tlen * sizeof(int));
  tr->pos       = ReallocOrDie (tr->pos,       tlen * sizeof(int));
}
void 
CP9FreeTrace(CP9trace_t *tr)
{
  if (tr == NULL) return;
  free(tr->pos);
  free(tr->nodeidx);
  free(tr->statetype);
  free(tr);
}


/* Function: CP9_2sub_cp9()
 * EPN 09.24.06
 * 
 * Purpose:  Given a template CM Plan 9 HMM, build a sub-model that
 *           models only a subset of the consensus columns of the
 *           original alignment. This requires a bit of care for
 *           the initial and final node of the sub CP9, and 
 *           straightforward copying of parameters for the rest.
 *          
 *           The new CP9 is constructed in Global Needleman/Wunsch
 *           mode. The orig_hmm MUST be in global mode. THIS IS
 *           CHECKED FOR IN A VERY FRAGILE MANNER!
 *       
 *           The approach here is to allocate and fill the new 
 *           sub CP9. There might be a better way - transforming
 *           the original CP9 into the new sub CP9 using a method
 *           involving pointer rearrangement, but I'm not sure
 *           how to do this.
 *             
 * Args:     orig_hmm    - the CP9 model w/ data-dep prob's valid
 *           ret_sub_hmm - the new sub CP9 hmm, allocated here, must
 *                         be freed by caller.
 *           spos        - first consensus column modelled by original
 *                         CP9 HMM the sub CP9 HMM models.
 *           epos        - final consensus column modelled by original
 *                         CP9 HMM the sub CP9 HMM models.
 *           orig_phi    - the 2D phi array for the original CP9 HMM.         
 * Return:   (void)
 *           HMM probabilities are modified.
 */
void
CP9_2sub_cp9(struct cplan9_s *orig_hmm, struct cplan9_s **ret_sub_hmm, int spos, int epos, double **orig_phi)
{
  struct cplan9_s       *sub_hmm;       
  int i, x;
  int orig_pos;

  sub_hmm = AllocCPlan9((epos-spos+1));

  for(x = 0; x < MAXABET; x++)
    {
      sub_hmm->null[x] = orig_hmm->null[x];
    }
  /* No special (*x* states in Plan 7) states in CM Plan 9 */

  /* First we just copy the parameters for spos..epos from the template HMM.
   * This is *slightly* wasteful, as we'll overwrite a few of these later.
   */
  for(i = 0; i <= (epos-spos+1); i++)
    {
      orig_pos = i + spos - 1;

      if(i > 0)
	{
	  for(x = 0; x < MAXABET; x++)
	    {
	      sub_hmm->mat[i][x] = orig_hmm->mat[orig_pos][x];
	      sub_hmm->msc[x][i] = orig_hmm->msc[x][orig_pos];
	    }

	  sub_hmm->begin[i]   = orig_hmm->begin[orig_pos];
	  sub_hmm->end[i]     = orig_hmm->end[orig_pos];
	  sub_hmm->bsc[i]     = orig_hmm->bsc[orig_pos];
	  sub_hmm->esc[i]     = orig_hmm->esc[orig_pos];
	  if((i > 1) && ((0. - sub_hmm->begin[i] > 0.00000001) ||
			 (sub_hmm->begin[i] - 0. > 0.00000001)))
	    {
	      Die("ERROR in cp9_2sub_cp9() is original CP9 HMM not in global (NW) mode? i: %d\n", i);
	    }
	}
      for(x = 0; x < MAXABET; x++)
	{
	  sub_hmm->ins[i][x] = orig_hmm->ins[orig_pos][x];
	  sub_hmm->isc[x][i] = orig_hmm->isc[x][orig_pos];
	}

      for(x = 0; x < 10; x++)
	{
	  sub_hmm->t[i][x]   = orig_hmm->t[orig_pos][x];
	  sub_hmm->tsc[x][i] = orig_hmm->tsc[x][orig_pos];
	}
      
    }

  /* Make the necessary modifications. */
  CP9_reconfig2sub(sub_hmm, spos, epos, 1, sub_hmm->M, orig_phi);

  sub_hmm->el_self   = orig_hmm->el_self;
  sub_hmm->el_selfsc = orig_hmm->el_selfsc;

  sub_hmm->flags |= CPLAN9_HASBITS;	/* raise the log-odds ready flag */
  *ret_sub_hmm = sub_hmm;
  return;
}

/* Function: CP9_reconfig2sub()
 * EPN 10.16.06
 * 
 * Purpose:  Given a CM Plan 9 HMM and a start position
 *           (spos) and end position (epos) that a sub CM models, 
 *           reconfigure the HMM so that it can only start in the 
 *           node that models spos (spos_nd) end in the node that 
 *           models epos (epos_nd).
 *
 *           If we're reconfiguring a CP9 HMM that ONLY models the
 *           consensus columns spos to epos, then spos_nd == 1 
 *           and epos_nd == hmm->M, but this is not necessarily true.
 *           We may be reconfiguring a CP9 HMM that models the
 *           full alignment including positions before and/or after
 *           spos and epos. In this case spos_nd == spos and
 *           epos_nd == epos;
 *           
 * Args:     hmm         - the CP9 model w/ data-dep prob's valid
 *           spos        - first consensus column modelled by some original
 *                         full length, template CP9 HMM that 'hmm' models.
 *           epos        - final consensus column modelled by some original
 *                         CP9 HMM that 'hmm' models.
 *           spos_nd     - the node of 'hmm' that models spos.
 *                         (1 if 'hmm' only has (epos-spos+1) nodes 
 *                         (spos if 'hmm' has a node for each column of original aln)
 *           epos_nd     - the node of the 'hmm' in that models epos.
 *                         (hmm->M if 'hmm' only has (epos-spos+1) nodes 
 *                         (epos if 'hmm' has a node for each column of original aln)
 *           orig_phi    - the 2D phi array for the original CP9 HMM.         
 * Return:   (void)
 *           HMM probabilities are modified.
 */
void
CP9_reconfig2sub(struct cplan9_s *hmm, int spos, int epos, int spos_nd,
		 int epos_nd, double **orig_phi)
{
  /* Make the necessary modifications. Since in cmalign --sub mode this
   * function will be called potentially once for each sequence, we 
   * don't want to call CP9Logoddsify(), but rather only logoddsify
   * the parameters that are different.
   */

  /* Configure entry.
   * Exactly 3 ways to start, B->M_1 (hmm->begin[1]), B->I_0 (hmm->t[0][CTMI]),
   *                      and B->D_1 (hmm->t[0][CTMD])
   */
  /* prob of starting in M_spos is (1. - prob of starting in I_spos-1) as there is no D_spos-1 -> M_spos trans */
      
  if(spos > 1)
    {
      hmm->begin[spos_nd] = 1.-((orig_phi[spos-1][HMMINSERT] * (1. - hmm->t[spos-1][CTII])) + 
			        (orig_phi[spos  ][HMMDELETE] - (orig_phi[spos-1][HMMINSERT] * hmm->t[spos-1][CTID])));
      hmm->t[spos_nd-1][CTMI] =   (orig_phi[spos-1][HMMINSERT] * (1. - hmm->t[spos-1][CTII]));
      hmm->t[spos_nd-1][CTMD] =    orig_phi[spos  ][HMMDELETE] - (orig_phi[spos-1][HMMINSERT] * hmm->t[spos-1][CTID]);
      hmm->t[spos_nd-1][CTMM] = 0.; /* probability of going from B(M_0) to M_1 is begin[1] */
      hmm->t[spos_nd-1][CTME] = 0.; /* can't go to EL from B(M_0) */
      hmm->t[spos_nd-1][CTDM] = 0.; /* D_0 doesn't exist */
      hmm->t[spos_nd-1][CTDI] = 0.; /* D_0 doesn't exist */
      hmm->t[spos_nd-1][CTDD] = 0.; /* D_0 doesn't exist */
      
      hmm->bsc[spos_nd]       = Prob2Score(hmm->begin[1], 1.0);

      hmm->tsc[CTMM][spos_nd-1] = -INFTY; /* probability of going from B(M_0) to M_1 is begin[1] */
      hmm->tsc[CTME][spos_nd-1] = -INFTY; 
      hmm->tsc[CTDM][spos_nd-1] = -INFTY; /* D_0 doesn't exist */
      hmm->tsc[CTDI][spos_nd-1] = -INFTY; /* D_0 doesn't exist */
      hmm->tsc[CTDD][spos_nd-1] = -INFTY; /* D_0 doesn't exist */
      
      hmm->tsc[CTMI][spos_nd-1] = Prob2Score(hmm->t[spos_nd-1][CTMI], 1.0);
      hmm->tsc[CTMD][spos_nd-1] = Prob2Score(hmm->t[spos_nd-1][CTMD], 1.0);
    }

  if(epos < hmm->M)
    {
      hmm->end[epos_nd]      = hmm->t[epos][CTMM] + hmm->t[epos][CTMD];
      hmm->t[epos_nd][CTDM] += hmm->t[epos][CTDD];
      hmm->t[epos_nd][CTIM] += hmm->t[epos][CTID];
      hmm->t[epos_nd][CTMM]  = 0.; /* M->E is actually end[M] */
      hmm->t[epos_nd][CTME]  = 0.; 
      hmm->t[epos_nd][CTMD]  = 0.; /* D_M+1 doesn't exist */
      hmm->t[epos_nd][CTDD]  = 0.; /* D_M+1 doesn't exist */
      hmm->t[epos_nd][CTID]  = 0.; /* D_M+1 doesn't exist */
      
      hmm->esc[epos_nd]       = Prob2Score(hmm->end[epos_nd], 1.0);
      hmm->tsc[CTDM][epos_nd] = Prob2Score(hmm->t[epos_nd][CTDM], 1.0);
      hmm->tsc[CTIM][epos_nd] = Prob2Score(hmm->t[epos_nd][CTIM], 1.0);
      hmm->tsc[CTMM][epos_nd] = -INFTY; /* M->E is actually end[M] */
      hmm->tsc[CTME][epos_nd] = -INFTY; 
      hmm->tsc[CTMD][epos_nd] = -INFTY; /* D_M+1 doesn't exist */
      hmm->tsc[CTDD][epos_nd] = -INFTY; /* D_M+1 doesn't exist */
      hmm->tsc[CTID][epos_nd] = -INFTY; /* D_M+1 doesn't exist */
    }
  hmm->flags |= CPLAN9_HASBITS;	/* raise the log-odds ready flag */

  return;
}

/***********************************************************
 * NOTE: INCOMPLETE! DO NOT USE WITHOUT FINISHING!
 *
 * Function: sub_CPlan9GlobalConfig()
 * EPN 09.24.06
 * based on SRE's Plan7GlobalConfig() from HMMER's plan7.c
 * 
 * Purpose:  Set the alignment independent parameters of
 *           a CM Plan 9 model to (Needleman/Wunsch) configuration
 *           for a 'sub' model. Used for sub_cm alignment, in which
 *           the model only aligns to a contiguous subset of 
 *           consensus columns from spos..epos inclusively. The 
 *           local entry is set to 1.0 for spos and exit
 *           is set to 1.0 for epos. Also transitions into 
 *           inserts and deletes that correspond to spos, and
 *           out of inserts and deletes that correspond to 
 *           epos are made impossible (see code).
 *           
 * Args:     hmm    - the CM Plan 9 model w/ data-dep prob's valid
 *           spos   - first consensus column modelled by original
 *                    CP9 HMM the sub CP9 HMM models.
 *           epos   - final consensus column modelled by original
 *                    CP9 HMM the sub CP9 HMM models.
 *           phi    - the 2D phi array for the original CP9 HMM.         
 * Return:   (void)
 *           HMM probabilities are modified.
 */
void
sub_CPlan9GlobalConfig(struct cplan9_s *hmm, int spos, int epos, double **phi)
{
  int i;
  Die("ERROR, sub_CPlan9GlobalConfig() was abandoned on 09.27.06 by EPN. It is incomplete.\n");
  /* No special (*x* states in Plan 7) states in CM Plan 9 */

  /* Configure entry.
   * Exactly 3 ways to start, B->M_1 (hmm->begin[1]), B->I_0 (hmm->t[0][CTMI]),
   *                      and B->D_1 (hmm->t[0][CTMD])
   */
  if(spos > 1)
    {
      /* prob of starting in M_spos is (1. - prob of starting in I_spos-1) as there is no D_spos-1 -> M_spos trans */
      hmm->begin[spos] = 1. - ((phi[spos-1][HMMINSERT] * (1. - hmm->t[spos-1][CTII])) + 
			       (phi[spos  ][HMMDELETE] - (phi[spos-1][HMMINSERT] * hmm->t[spos-1][CTID])));
      hmm->t[spos-1][CTMI] =   (phi[spos-1][HMMINSERT] * (1. - hmm->t[spos-1][CTII]));
      hmm->t[spos-1][CTMD] =    phi[spos  ][HMMDELETE] - (phi[spos-1][HMMINSERT] * hmm->t[spos-1][CTID]);
      hmm->t[spos-1][CTMM] = 0.; /* probability of going from B to M_1 is begin[1] */

      /* eliminate the possibility of going from B -> I_0, the only other way to 
       * start a parse besides B -> M_1 */
      hmm->t[0][CTMI] = 0.;
    }

  for(i = 1; i < spos; i++)
    hmm->begin[i] = 0.;
  for(i = spos+1; i <= hmm->M; i++)
    hmm->begin[i] = 0.;

  if(epos < hmm->M)
    {
      //      hmm->end[epos]      = hmm->t[epos][CTMM] + (phi[epos] * hmm->t[epos][CTMD]) + (;
      hmm->end[hmm->M] = hmm->t[epos][CTMM] + hmm->t[epos][CTMD];
      hmm->t[epos][CTDM] += hmm->t[epos][CTDD];
      hmm->t[epos][CTIM] += hmm->t[epos][CTID];
      hmm->t[epos][CTDD]  = 0.; /* no D state in final node M */
      hmm->t[epos][CTID]  = 0.; /* no D state in final node M */
    }
  /* EPN 09.27.06 NEED SOME WAY OF ENSURING THAT NODE epos+1 IS NEVER REACHED. */

  for(i = 1; i < epos; i++)
    hmm->end[i] = 0.;
  for(i = epos+1; i <= hmm->M; i++)
    hmm->end[i] = 0.;

  CPlan9RenormalizeExits(hmm, spos);
  hmm->flags       &= ~CPLAN9_HASBITS; /* reconfig invalidates log-odds scores */
}

/************************************************************************
 * Functions stolen from HMMER 2.4 for use with CM plan 9 HMMs.
 * Eventually, these should go away, replaced with Easel funcs. 
 * These first 4 were stolen from HMMER:mathsupport.c
 * 
 * ILogSum() (and auxiliary funcs associated with it)
 * Score2Prob()
 * Prob2Score()
 * Scorify()
 * 
 * And 1 was stolen from HMMER:plan7.c:
 * DegenerateSymbolScore()
 *                           
 ************************************************************************/
/* Function: ILogsum()
 * 
 * Purpose:  Return the scaled integer log probability of
 *           the sum of two probabilities p1 and p2, where
 *           p1 and p2 are also given as scaled log probabilities.
 *         
 *           log(exp(p1)+exp(p2)) = p1 + log(1 + exp(p2-p1)) for p1 > p2
 *           
 *           For speed, builds a lookup table the first time it's called.
 *           LOGSUM_TBL is set to 20000 by default, in config.h.
 *
 *           Because of the one-time initialization, we have to
 *           be careful in a multithreaded implementation... hence
 *           the use of pthread_once(), which forces us to put
 *           the initialization routine and the lookup table outside
 *           ILogsum(). (Thanks to Henry Gabb at Intel for pointing
 *           out this problem.)
 *           
 * Args:     p1,p2 -- scaled integer log_2 probabilities to be summed
 *                    in probability space.
 *                    
 * Return:   scaled integer log_2 probability of the sum.
 */

static int ilogsum_lookup[LOGSUM_TBL];
static void 
init_ilogsum(void)
{
  int i;
  for (i = 0; i < LOGSUM_TBL; i++) 
    ilogsum_lookup[i] = (int) (INTSCALE * 1.44269504 * 
	   (log(1.+exp(0.69314718 * (float) -i/INTSCALE))));
}
int 
ILogsum(int p1, int p2)
{
  int    diff;
#ifdef HMMER_THREADS
  static pthread_once_t firsttime = PTHREAD_ONCE_INIT;
  pthread_once(&firsttime, init_ilogsum);
#else
  static int firsttime = 1;
  if (firsttime) { init_ilogsum(); firsttime = 0; }
#endif
  if(p1 == -INFTY) return p2; /* EPN */
  if(p2 == -INFTY) return p1; /* EPN */

  diff = p1-p2;
  if      (diff >=  LOGSUM_TBL) return p1;
  else if (diff <= -LOGSUM_TBL) return p2;
  else if (diff > 0)            return p1 + ilogsum_lookup[diff];
  else                          return p2 + ilogsum_lookup[-diff];
} 

/* Function: Prob2Score()
 * 
 * Purpose:  Convert a probability to a scaled integer log_2 odds score. 
 *           Round to nearest integer (i.e. note use of +0.5 and floor())
 *           Return the score. 
 */
int
Prob2Score(float p, float null)
{
  if   (p == 0.0) return -INFTY;
  else            return (int) floor(0.5 + INTSCALE * sreLOG2(p/null));
}

/* Function: Score2Prob()
 * 
 * Purpose:  Convert an integer log_2 odds score back to a probability;
 *           needs the null model probability, if any, to do the conversion.
 */
float 
Score2Prob(int sc, float null)
{
  if (sc == -INFTY) return 0.;
  else              return (null * sreEXP2((float) sc / INTSCALE));
}


/* Function: Scorify()
 * 
 * Purpose:  Convert a scaled integer log-odds score to a floating
 *           point score for output. (could be a macro but who cares.)
 */
float 
Scorify(int sc)
{
  return ((float) sc / INTSCALE);
}

/* Stolen from HMMER-2.4::cplan9.c()*/
/* Function: DegenerateSymbolScore()
 * 
 * Purpose:  Given a sequence character x and an hmm emission probability
 *           vector, calculate the log-odds (base 2) score of
 *           the symbol.
 *          
 *           Easy if x is in the emission alphabet, but not so easy
 *           is x is a degenerate symbol. The "correct" Bayesian
 *           philosophy is to calculate score(X) by summing over
 *           p(x) for all x in the degenerate symbol X to get P(X),
 *           doing the same sum over the prior to get F(X), and
 *           doing log_2 (P(X)/F(X)). This gives an X a zero score,
 *           for instance.
 *           
 *           Though this is correct in a formal Bayesian sense --
 *           we have no information on the sequence, so we can't
 *           say if it's random or model, so it scores zero --
 *           it sucks, big time, for scoring biological sequences.
 *           Sequences with lots of X's score near zero, while
 *           real sequences have average scores that are negative --
 *           so the X-laden sequences appear to be lifted out
 *           of the noise of a full histogram of a database search.
 *           Correct or not, this is highly undesirable.
 *           
 *           So therefore we calculated the expected score of
 *           the degenerate symbol by summing over all x in X:
 *                 e_x log_2 (p(x)/f(x))
 *           where the expectation of x, e_x, is calculated from
 *           the random model.
 *
 *           Empirically, this works; it also has a wooly hand-waving
 *           probabilistic justification that I'm happy enough about.
 *           
 * Args:     p      - probabilities of normal symbols
 *           null   - null emission model
 *           ambig  - index of the degenerate character in Alphabet[]
 *                    
 * Return:   the integer log odds score of x given the emission
 *           vector and the null model, scaled up by INTSCALE.              
 */
int 
DegenerateSymbolScore(float *p, float *null, int ambig)
{
  int x;
  float numer = 0.;
  float denom = 0.;

  for (x = 0; x < Alphabet_size; x++) {
    if (Degenerate[ambig][x]) {
      numer += null[x] * sreLOG2(p[x] / null[x]);
      denom += null[x];
    }
  }
  return (int) (INTSCALE * numer / denom);
}

/* Function:  CP9HackInsertScores()
 * Incept:    EPN, Fri Feb  9 10:59:12 2007
 *
 * Purpose:   Make all inserts 0. Usually called from CMHackInsertScores()
 *            to make the HMM inserts match the CM inserts.
 *
 * Args:      cp9 - the CP9 HMM 
 *
 * Returns:   (void)
 */
void
CP9HackInsertScores(CP9_t *cp9)
{
  int k, x;
  for (k = 0; k <= cp9->M; k++)
    /* CP9 HMMs have insert states in nodes 0 and M */
    for (x = 0; x < MAXDEGEN; x++)
      cp9->isc[x][k] = 0.;
}
/* Function:  CP9EnforceHackMatchScores()
 * Incept:    EPN, Fri Feb  9 11:06:31 2007
 *
 * Purpose:   Make all match emissions 0, except those enforce
 *            a specified subsequence (it's assumed the CP9 
 *            is already set up for this enforcement). 
 *
 * Args:      cp9           - the CP9 HMM 
 *            enf_start_pos - first posn of enforced subseq
 *            enf_end_pos   - last  posn of enforced subseq
 * Returns:   (void)
 */
void
CP9EnforceHackMatchScores(CP9_t *cp9, int enf_start_pos, int enf_end_pos)
{
  int k, x;
  for (k = 1; k < enf_start_pos; k++) /* M_0 is the begin state, it's silent */
    for (x = 0; x < MAXDEGEN; x++)
      cp9->msc[x][k] = 0.;
  for (k = enf_end_pos+1; k <= cp9->M; k++)
    for (x = 0; x < MAXDEGEN; x++)
      cp9->msc[x][k] = 0.;
}

/* Function: CP9_fake_tracebacks()
 * 
 * Purpose:  From a consensus assignment of columns to MAT/INS, construct fake
 *           tracebacks for each individual sequence.
 *           
 * Args:     aseqs     - alignment [0..nseq-1][0..alen-1]
 *           nseq      - number of seqs in alignment
 *           alen      - length of alignment in columns
 *           matassign - assignment of column 1 if MAT, 0 if INS; 
 *                       [1..alen] (off one from aseqs)
 *           ret_tr    - RETURN: array of tracebacks
 *           
 * Return:   (void)
 *           ret_tr is alloc'ed here. Caller must free.
 */          
void
CP9_fake_tracebacks(char **aseq, int nseq, int alen, int *matassign,
		CP9trace_t ***ret_tr)
{
  CP9trace_t **tr;
  int  idx;                     /* counter over sequences          */
  int  i;                       /* position in raw sequence (1..L) */
  int  k;                       /* position in HMM                 */
  int  apos;                    /* position in alignment columns   */
  int  tpos;			/* position in traceback           */
  int  first_match;             /* first match column */
  int  last_match;              /* last match column */

  tr = (CP9trace_t **) MallocOrDie (sizeof(CP9trace_t *) * nseq);
  
  first_match = -1;
  last_match  = -1;
  for (apos = 0; apos < alen; apos++)
    {
      if(matassign[apos+1] && first_match == -1) first_match = apos;
      if(matassign[apos+1]) last_match = apos;
    }

  for (idx = 0; idx < nseq; idx++)
    {
      CP9AllocTrace(alen+2, &tr[idx]);  /* allow room for B & E */
      
				/* all traces start with M_0 state (the B state)... */
      tr[idx]->statetype[0] = CSTB;
      tr[idx]->nodeidx[0]   = 0;
      tr[idx]->pos[0]       = 0;

      i = 1;
      k = 0;
      tpos = 1;

      for (apos = 0; apos < alen; apos++)
        {
	  tr[idx]->statetype[tpos] = CSTBOGUS; /* bogus, deliberately, to debug */

	  if (matassign[apos+1] && ! isgap(aseq[idx][apos]))
	    {			/* MATCH */
	      k++;		/* move to next model pos */
	      tr[idx]->statetype[tpos] = CSTM;
	      tr[idx]->nodeidx[tpos]   = k;
	      tr[idx]->pos[tpos]       = i;
	      i++;
	      tpos++;
	    }	      
          else if (matassign[apos+1])
            {                   /* DELETE */
	      /* We should be careful about S/W transitions; but we have 
	       * an ambiguity, based on the MSA, we can't tell if we
	       * did a local begin (some M->E transition) or if we
	       * went through a bunch of D state's before the first match 
	       * B->D_1 -> D_2 .... -> M_x. For now, we assume we're not in
	       * S/W mode, and treat it as the latter case, see
	       * HMMER's modelmaker.c:fake_tracebacks() for code
	       * on one *would* implement the S/W consideration IF
	       * there wasn't a B->D_1 transition allowed.
	       */
	      k++;		/* *always* move on model when match column seen */
	      tr[idx]->statetype[tpos] = CSTD;
	      tr[idx]->nodeidx[tpos]   = k;
	      tr[idx]->pos[tpos]       = 0;
	      tpos++;
            }
	  else if (! isgap(aseq[idx][apos]))
	    {			/* INSERT */
	      tr[idx]->statetype[tpos] = CSTI;
              tr[idx]->nodeidx[tpos]   = k;
              tr[idx]->pos[tpos]       = i;
	      i++;
	      tpos++;
	    }
	}
       /* all traces end with E state */
      /* We should be careful about S/W transitions; but we have 
       * an ambiguity, based on the MSA, we can't tell if we
       * did a local end (some M->E transition) or if we
       * went through a bunch of D state's before the final 
       * D_M -> E transition. For now, we assume we're not in
       * S/W mode, and treat it as the latter case, see
       * HMMER's modelmaker.c:fake_tracebacks() for code
       * on one *would* implement the S/W consideration IF
       * there wasn't a D_M -> E transition allowed.
       */
      tr[idx]->statetype[tpos] = CSTE;
      tr[idx]->nodeidx[tpos]   = 0;
      tr[idx]->pos[tpos]       = 0;
      tpos++;
      tr[idx]->tlen = tpos;
    }    /* end for sequence # idx */

  *ret_tr = tr;
  return;
}

/* Function: CP9TraceCount() 
 * EPN 09.04.06 based on Eddy's P7TraceCount() from HMMER's trace.c
 * 
 * Purpose:  Count a traceback into a count-based HMM structure.
 *           (Usually as part of a model parameter re-estimation.)
 *           
 * Args:     hmm   - counts-based CM Plan 9 HMM
 *           dsq   - digitized sequence that traceback aligns to the HMM (1..L)
 *           wt    - weight on the sequence
 *           tr    - alignment of seq to HMM
 *           
 * Return:   (void)
 */
void
CP9TraceCount(struct cplan9_s *hmm, char *dsq, float wt, CP9trace_t *tr)
{
  int tpos;                     /* position in tr */
  int i;			/* symbol position in seq */
  
  for (tpos = 0; tpos < tr->tlen; tpos++)
    {
      i = tr->pos[tpos];

      /* Emission counts. 
       */
      if (tr->statetype[tpos] == CSTM) 
	SingletCount(hmm->mat[tr->nodeidx[tpos]], dsq[i], wt);
      else if (tr->statetype[tpos] == CSTI) 
	SingletCount(hmm->ins[tr->nodeidx[tpos]], dsq[i], wt);

      /* State transition counts
       */
      switch (tr->statetype[tpos]) {
      case CSTB:
	switch (tr->statetype[tpos+1]) {
	case CSTM: hmm->begin[tr->nodeidx[tpos+1]] += wt; break;
	case CSTI: hmm->t[0][CTMI]                 += wt; break;
	case CSTD: hmm->t[0][CTMD]                 += wt; break;
	default:      
	  Die("illegal state transition %s->%s in traceback", 
	      CP9Statetype(tr->statetype[tpos]), 
	      CP9Statetype(tr->statetype[tpos+1]));
	}
	break;
      case CSTM:
	switch (tr->statetype[tpos+1]) {
	case CSTM: hmm->t[tr->nodeidx[tpos]][CTMM] += wt; break;
	case CSTI: hmm->t[tr->nodeidx[tpos]][CTMI] += wt; break;
	case CSTD: hmm->t[tr->nodeidx[tpos]][CTMD] += wt; break;
	case CSTE: hmm->end[tr->nodeidx[tpos]]     += wt; break;
	default:    
	  Die("illegal state transition %s->%s in traceback", 
	      CP9Statetype(tr->statetype[tpos]), 
	      CP9Statetype(tr->statetype[tpos+1]));
	}
	break;
      case CSTI:
	switch (tr->statetype[tpos+1]) {
	case CSTM: hmm->t[tr->nodeidx[tpos]][CTIM] += wt; break;
	case CSTI: hmm->t[tr->nodeidx[tpos]][CTII] += wt; break;
	case CSTD: hmm->t[tr->nodeidx[tpos]][CTID] += wt; break;
	case CSTE: 
	  /* This should only happen from the final insert (I_M) state */
	  if((tpos+1) != (tr->tlen-1))
	    Die("illegal state transition %s->%s (I is not final insert) in traceback", 
		CP9Statetype(tr->statetype[tpos]), 
		CP9Statetype(tr->statetype[tpos+1]));
	  hmm->t[tr->nodeidx[tpos]][CTIM] += wt; break;
	  break;
	default:    
	  Die("illegal state transition %s->%s in traceback", 
	      CP9Statetype(tr->statetype[tpos]), 
	      CP9Statetype(tr->statetype[tpos+1]));
	}
	break;
      case CSTD:
	switch (tr->statetype[tpos+1]) {
	case CSTM: hmm->t[tr->nodeidx[tpos]][CTDM] += wt; break;
	case CSTI: hmm->t[tr->nodeidx[tpos]][CTDI] += wt; break;
	case CSTD: hmm->t[tr->nodeidx[tpos]][CTDD] += wt; break;
	case CSTE: 
	  /* This should only happen from the final delete (D_M) state */
	  if((tpos+1) != (tr->tlen-1))
	    Die("illegal state transition %s->%s (D is not final delete) in traceback", 
		CP9Statetype(tr->statetype[tpos]), 
		CP9Statetype(tr->statetype[tpos+1]));
	  hmm->t[tr->nodeidx[tpos]][CTDM] += wt; break;
	  break;
	default:    
	  Die("illegal state transition %s->%s in traceback", 
	      CP9Statetype(tr->statetype[tpos]), 
	      CP9Statetype(tr->statetype[tpos+1]));
	}
	break;
      case CSTE:
	break; /* E is the last. It makes no transitions. */

      default:
	Die("illegal state %s in traceback", 
	    CP9Statetype(tr->statetype[tpos]));
      }
    }
}


/* Function: CP9Statetype()
 * 
 * Purpose:  Returns the state type in text.
 * Example:  CP9Statetype(M) = "M"
 */
char *
CP9Statetype(char st)
{
  switch (st) {
  case CSTM: return "M";
  case CSTD: return "D";
  case CSTI: return "I";
  case CSTB: return "B";
  case CSTE: return "E";
  case CSTEL: return "L";
  default: return "BOGUS";
  }
}

/* Function: CP9TraceScore()
 *           based on HMMER 2.3.2's P7TraceScore by SRE
 *
 * Purpose:  Score a traceback and return the score in scaled bits.
 * Incept:   EPN, Wed May 30 06:07:14 2007
 *           
 * Args:     hmm   - HMM with valid log odds scores.
 *           dsq   - digitized sequence that traceback aligns to the HMM (1..L)
 *           tr    - alignment of seq to HMM
 *           
 * Return:   (void)
 */
float
CP9TraceScore(CP9_t *hmm, char *dsq, CP9trace_t *tr)
{
  int score;			/* total score as a scaled integer */
  int tpos;                     /* position in tr */
  char sym;		        /* digitized symbol in dsq */
  
  /*CP9PrintTrace(stdout, tr, hmm, dsq); */
  score = 0;
  for (tpos = 0; tpos < tr->tlen-1; tpos++)
    {
      sym = dsq[tr->pos[tpos]];

      /* Emissions from M and I states.
       */
      if (tr->statetype[tpos] == CSTM) 
	score += hmm->msc[(int) sym][tr->nodeidx[tpos]];
      else if (tr->statetype[tpos] == CSTI) 
	score += hmm->isc[(int) sym][tr->nodeidx[tpos]];

      /* State transitions. Including EL emissions, EL emits on transition 
       */
      score += CP9TransitionScoreLookup(hmm, 
				     tr->statetype[tpos], tr->nodeidx[tpos],
				     tr->statetype[tpos+1], tr->nodeidx[tpos+1]);
    }
  return Scorify(score);
}

/* Function: CP9PrintTrace()
 *           based on HMMER's 2.3.2 P7PrintTrace()
 *
 * Purpose:  Print out a traceback structure.
 *           If hmm is non-NULL, also print transition and emission scores.
 * Incept:   EPN, Wed May 30 06:07:57 2007
 *           
 * Args:     fp  - stderr or stdout, often
 *           tr  - trace structure to print
 *           hmm - NULL or hmm containing scores to print
 *           dsq - NULL or digitized sequence trace refers to.                
 */
void
CP9PrintTrace(FILE *fp, CP9trace_t *tr, CP9_t *hmm, char *dsq)
{
  int          tpos;		/* counter for trace position */
  unsigned int sym;
  int          sc; 

  if (tr == NULL) {
    fprintf(fp, " [ trace is NULL ]\n");
    return;
  }

  if (hmm == NULL) {
    fprintf(fp, "st  node   rpos  - traceback len %d\n", tr->tlen);
    fprintf(fp, "--  ---- ------\n");
    for (tpos = 0; tpos < tr->tlen; tpos++) {
      fprintf(fp, "%1s  %4d %6d\n", 
	      CP9Statetype(tr->statetype[tpos]),
	      tr->nodeidx[tpos],
	      tr->pos[tpos]);
    } 
  } else {
    if (!(hmm->flags & CPLAN9_HASBITS))
      Die("oi, you can't print scores from that hmm, it's not ready.");

    sc = 0;
    fprintf(fp, "st  node   rpos  transit emission - traceback len %d\n", tr->tlen);
    fprintf(fp, "--  ---- ------  ------- --------\n");
    for (tpos = 0; tpos < tr->tlen; tpos++) {
      if (dsq != NULL) sym = dsq[tr->pos[tpos]];

      fprintf(fp, "%1s  %4d %6d  %7d", 
	      CP9Statetype(tr->statetype[tpos]),
	      tr->nodeidx[tpos],
	      tr->pos[tpos],
	      (tpos < tr->tlen-1) ? 
	      CP9TransitionScoreLookup(hmm, tr->statetype[tpos], tr->nodeidx[tpos],
				    tr->statetype[tpos+1], tr->nodeidx[tpos+1]) : 0);

      if (tpos < tr->tlen-1)
	sc += CP9TransitionScoreLookup(hmm, tr->statetype[tpos], tr->nodeidx[tpos],
				       tr->statetype[tpos+1], tr->nodeidx[tpos+1]);

      if (dsq != NULL) {
	if (tr->statetype[tpos] == CSTM)  
	  {
	    fprintf(fp, " %8d %c", hmm->msc[(int) sym][tr->nodeidx[tpos]], 
		    Alphabet[(int) sym]);
	    sc += hmm->msc[(int) sym][tr->nodeidx[tpos]];
	  }
	else if (tr->statetype[tpos] == CSTI) 
	  {
	    fprintf(fp, " %8d %c", hmm->isc[(int) sym][tr->nodeidx[tpos]], 
		    (char) tolower((int) Alphabet[(int) sym]));
	    sc += hmm->isc[(int) sym][tr->nodeidx[tpos]];
	  }
	else if (tr->statetype[tpos] == CSTEL) 
	  {
	    if(tr->statetype[tpos-1] == CSTEL) /* we will emit on self transit */
	      {
		fprintf(fp, " %8d %c", 0,
			(char) tolower((int) Alphabet[(int) sym]));
	      }
	    else /* we just entered EL, no emission */
	      {
		fprintf(fp, " %8s %c", "-", '-');
	      }
	  }
	/*	else if ((tr->statetype[tpos] == STN && tr->statetype[tpos-1] == STN) ||
		 (tr->statetype[tpos] == STC && tr->statetype[tpos-1] == STC) ||
		 (tr->statetype[tpos] == STJ && tr->statetype[tpos-1] == STJ))
	  {
	    fprintf(fp, " %8d %c", 0, (char) tolower((int) hmmer_Alphabet[(int) sym]));
	    }*/
      } else {
	fprintf(fp, " %8s %c", "-", '-');
      }


      fputs("\n", fp);
    }
    fprintf(fp, "                 ------- --------\n");
    fprintf(fp, "           total: %6d\n\n", sc);
  }
}

/* Function: CP9TransitionScoreLookup()
 *           based on HMMER's 2.3.2 function of same name
 *
 * Incept:   EPN, Wed May 30 06:09:04 2007
 * Purpose:  Convenience function used in CP9PrintTrace() and CP9TraceScore();
 *           given state types and node indices for a transition,
 *           return the integer score for that transition. 
 */
int
CP9TransitionScoreLookup(struct cplan9_s *hmm, char st1, int k1, 
			 char st2, int k2)
{
  switch (st1) {
  case CSTB:
    switch (st2) {
    case CSTM: return hmm->bsc[k2]; 
    case CSTI: return hmm->tsc[CTMI][0];
    case CSTD: return hmm->tsc[CTMD][0];
    default:      Die("illegal %s->%s transition", CP9Statetype(st1), CP9Statetype(st2));
    }
    break;
  case CSTM:
    switch (st2) {
    case CSTM: return hmm->tsc[CTMM][k1];
    case CSTI: return hmm->tsc[CTMI][k1];
    case CSTD: return hmm->tsc[CTMD][k1];
    case CSTE: return hmm->esc[k1];
    case CSTEL: return hmm->tsc[CTME][k1];
    default:      Die("illegal %s->%s transition", CP9Statetype(st1), CP9Statetype(st2));
    }
    break;
  case CSTI:
    switch (st2) {
    case CSTM: return hmm->tsc[CTIM][k1];
    case CSTI: return hmm->tsc[CTII][k1];
    case CSTD: return hmm->tsc[CTID][k1];
    case CSTE: return hmm->tsc[CTIM][k1]; /* This should only happen from the final insert (I_M) state */
    default:      Die("illegal %s->%s transition", CP9Statetype(st1), CP9Statetype(st2));
    }
    break;
  case CSTD:
    switch (st2) {
    case CSTM: return hmm->tsc[CTDM][k1]; 
    case CSTI: return hmm->tsc[CTDI][k1];
    case CSTD: return hmm->tsc[CTDD][k1];
    case CSTE: return hmm->tsc[CTDM][k1]; /* This should only happen from the final delete (D_M) state */
    default:      Die("illegal %s->%s transition", CP9Statetype(st1), CP9Statetype(st2));
    }
    break;
  case CSTEL:
    switch (st2) {
    case CSTM: return 0; /* transition to EL penalty incurred when M->EL transition takes place */
    case CSTE: return 0; /* transition to EL penalty incurred when M->EL transition takes place */
    case CSTEL: return hmm->el_selfsc; /* penalty for EL->EL self transition loop */
    default:      Die("illegal %s->%s transition", CP9Statetype(st1), CP9Statetype(st2));
    }
    break;
  case CSTE: /* this should never happen, it means we transitioned from E, which is not
	      * allowed. */
    Die("illegal %s->%s transition", CP9Statetype(st1), CP9Statetype(st2));
    break;
  default:        Die("illegal state %s in traceback", CP9Statetype(st1));
  }
  /*NOTREACHED*/
  return 0;
}


/* Function: CP9ViterbiTrace()
 * Date:     EPN, Wed May 30 17:32:05 2007
 *           based on HMMER 2.3.2's P7ViterbiTrace()
 *
 * Purpose:  Traceback of a Viterbi matrix: i.e. retrieval 
 *           of optimum alignment.
 *           
 * Args:     hmm    - hmm, log odds form, used to make mx
 *           dsq    - sequence aligned to (digital form) 1..N  
 *           i0     - first residue of sequence, often 1
 *           j0     - last residue of sequence, often L
 *           mx     - the matrix to trace back in, N x hmm->M
 *           ret_tr - RETURN: traceback.
 *           
 * Return:   (void)
 *           ret_tr is allocated here. Free using CP9FreeTrace().
 */
void
CP9ViterbiTrace(struct cplan9_s *hmm, char *dsq, int i0, int j0,
		struct cp9_dpmatrix_s *mx, CP9trace_t **ret_tr)
{
  CP9trace_t *tr;
  int curralloc;		/* current allocated length of trace */
  int tpos;			/* position in trace */
  int i;			/* position in seq (1..N) */
  int k;			/* position in model (1..M) */
  int *erow, **mmx, **imx, **dmx, **elmx;
  int sc;			/* temp var for pre-emission score */
  int error_flag; 
  int c;                        /* counter over possible EL states */

  /* Overallocate for the trace.
   * B- ... - E  : 2 states + N is minimum trace;
   * add N more as buffer.             
   */
  curralloc = (j0-i0+1) * 2 + 2; 
  CP9AllocTrace(curralloc, &tr);

  mmx = mx->mmx;
  imx = mx->imx;
  dmx = mx->dmx;
  elmx= mx->elmx;
  erow= mx->erow;

  /* Initialization of trace
   * We do it back to front; ReverseTrace() is called later.
   */
  tr->statetype[0] = CSTE;
  tr->nodeidx[0]   = 0;
  tr->pos[0]       = 0;
  tpos = 1;
  i    = j0;			/* current i (seq pos) we're trying to assign */

  /* Traceback
   */
  while (tr->statetype[tpos-1] != CSTB) {
    error_flag = FALSE;
    switch (tr->statetype[tpos-1]) {
    case CSTM:			/* M connects from i-1,k-1, B or an EL*/
      /*printf("CSTM k: %d i:%d \n", k, i);*/
      sc = mmx[i+1][k+1] - hmm->msc[(int) dsq[i+1]][k+1];
      if (sc <= -INFTY) { CP9FreeTrace(tr); *ret_tr = NULL; return; }
      else if (sc == mmx[i][k] + hmm->tsc[CTMM][k])
	{
	  tr->statetype[tpos] = CSTM;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = i--;
	}
      else if (sc == imx[i][k] + hmm->tsc[CTIM][k])
	{
	  tr->statetype[tpos] = CSTI;
	  tr->nodeidx[tpos]   = k;
	  tr->pos[tpos]       = i--;
	}
      else if (sc == dmx[i][k] + hmm->tsc[CTDM][k])
	{
	  tr->statetype[tpos] = CSTD;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = 0;
	}
      else /* Check if we came from an EL state (could be more than 1 choice) */ 
	{
	  error_flag = TRUE;
	  /* note we look at el_from_ct[k+1] not k for same reason we look
	   * at bsc[k+1] above, we're going backwards, this is a tricky off-by-one */
	  for(c = 0; c < hmm->el_from_ct[k+1]; c++) /* el_from_ct[k+1] is >= 0 */
	    {
	      /* transition penalty to EL incurred when EL was entered */
	      if(sc == elmx[i][hmm->el_from_idx[k+1][c]])
		{
		  tr->statetype[tpos] = CSTEL;
		  k = hmm->el_from_idx[(k+1)][c];
		  tr->nodeidx[tpos]   = k;
		  tr->pos[tpos]       = i--; 
		  error_flag = FALSE;
		  break;
		}
	    }
	}
      if(error_flag)
	{
	  /* one last possibility, we came from B, check this last, in
	   * case hmm->bsc[k+1] happens to be identical to sc but
	   * we're not done the parse yet (i.e. one of the cases
	   * above equaled sc). */
	  if (sc == hmm->bsc[k+1])
	    {
	      tr->statetype[tpos] = CSTB;
	      tr->nodeidx[tpos]   = 0;
	      tr->pos[tpos]       = 0;
	      if(tr->pos[tpos-1] != 1)
		Die("traceback failed: premature begin");
	      error_flag = FALSE;
	    }
	}
      if(error_flag)
	Die("traceback failed");
      break;

    case CSTD:			/* D connects from M,D,I, (D_1 also connects from B (M_0) */
      /*printf("CSTD k: %d i:%d \n", k, i);*/
      if (dmx[i][k+1] <= -INFTY) { CP9FreeTrace(tr); *ret_tr = NULL; return; }
      else if(k == 0) /* D_0 connects from B(M_0), and I_0 */
	{
	  if(dmx[i][k+1] == mmx[i][k] + hmm->tsc[CTMD][k])
	    {
	      tr->statetype[tpos] = CSTB;
	      tr->nodeidx[tpos]   = 0;
	      tr->pos[tpos]       = 0;
	    }
	  if (dmx[i][k+1] == imx[i][k] + hmm->tsc[CTID][k])
	    {
	      tr->statetype[tpos] = CSTI;
	      tr->nodeidx[tpos]   = k--;
	      tr->pos[tpos]       = i--;
	    }
	} /* else k != 0 */
      else if (dmx[i][k+1] == mmx[i][k] + hmm->tsc[CTMD][k])
	{
	  tr->statetype[tpos] = CSTM;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = i--;
	}
      else if (dmx[i][k+1] == imx[i][k] + hmm->tsc[CTID][k]) 
	{
	  tr->statetype[tpos] = CSTI;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = i--;
	}
      else if (dmx[i][k+1] == dmx[i][k] + hmm->tsc[CTDD][k]) 
	{
	  tr->statetype[tpos] = CSTD;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = 0;
	}
      else Die("traceback failed");
      break;

    case CSTI:			/* I connects from M,I,D, (I_0 connects from B also(*/
      /*printf("CSTI k: %d i:%d \n", k, i);*/
      sc = imx[i+1][k] - hmm->isc[(int) dsq[i+1]][k];
      if (sc <= -INFTY) { CP9FreeTrace(tr); *ret_tr = NULL; return; }
      else if(k == 0) /* I_0 connects from B(M_0), and I_0 */
	{
	  if(sc == mmx[i][k] + hmm->tsc[CTMI][k])
	    {
	      tr->statetype[tpos] = CSTB;
	      tr->nodeidx[tpos]   = 0;
	      tr->pos[tpos]       = 0;
	    }
	  if (sc == imx[i][k] + hmm->tsc[CTII][k])
	    {
	      tr->statetype[tpos] = CSTI;
	      tr->nodeidx[tpos]   = k;
	      tr->pos[tpos]       = i--;
	    }
	}
      /* else k != 0 */
      else if (sc == mmx[i][k] + hmm->tsc[CTMI][k])
	{
	  tr->statetype[tpos] = CSTM;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = i--;
	}

      else if (sc == imx[i][k] + hmm->tsc[CTII][k])
	{
	  tr->statetype[tpos] = CSTI;
	  tr->nodeidx[tpos]   = k;
	  tr->pos[tpos]       = i--;
	}
      else if (sc == dmx[i][k] + hmm->tsc[CTDI][k])
	{
	  tr->statetype[tpos] = CSTD;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = 0;
	}
      else Die("traceback failed");
      break;

    case CSTE:			/* E connects from any M state. k set here 
				 * also can connect from I_M or D_M (diff from p7) 
				 * or even EL_M if it exists */
      if (erow[i] <= -INFTY) { CP9FreeTrace(tr); *ret_tr = NULL; return; }
      if (erow[i] == imx[i][hmm->M] + hmm->tsc[CTIM][hmm->M])
	{
	  k = hmm->M;
	  tr->statetype[tpos] = CSTI;
	  tr->nodeidx[tpos]   = k;
	  tr->pos[tpos]       = i--;
	}
      else if (erow[i] == dmx[i][hmm->M] + hmm->tsc[CTDM][hmm->M])
	{
	  k = hmm->M;
	  tr->statetype[tpos] = CSTD;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = 0;
	}
      else
	{
	  error_flag = TRUE;
	  for (k = hmm->M; k >= 1; k--)
	    if (erow[i] == mmx[i][k] + hmm->esc[k])
	      {
		tr->statetype[tpos] = CSTM;
		tr->nodeidx[tpos]   = k--;
		tr->pos[tpos]       = i--;
		error_flag = FALSE;
		break;
	      }
	  if(error_flag)
	    {
	      /* Check if we came from an EL state (could be more than 1 choice) */ 
	      /* hmm->el_from-ct[hmm->M+1] is # of ELs that can transit to E (END) */
	      for(c = hmm->el_from_ct[hmm->M+1]-1; c >= 0; c--) /* el_from_ct[] is >= 0 */
		{
		  /* transition penalty to EL incurred when EL was entered */
		  if(erow[i] == elmx[i][hmm->el_from_idx[hmm->M+1][c]])
		    {
		      tr->statetype[tpos] = CSTEL;
		      k = hmm->el_from_idx[(hmm->M+1)][c];
		      tr->nodeidx[tpos]   = k;
		      tr->pos[tpos]       = i;
		      error_flag = FALSE;
		      break;
		    }
		}
	    }
	}
      if (k < 0 || error_flag) Die("traceback failed");
      break;

    case CSTEL:			/* EL connects from certain M states and itself */
      /*printf("CSTEL k: %d i:%d \n", k, i);*/
      /* check if we are staying in the EL */
      sc = elmx[i+1][k];
      if (sc == elmx[i][k] + hmm->el_selfsc) /* i >= 2, first residue must be emitted by a match, not an EL */
	{
	  tr->statetype[tpos] = CSTEL;
	  tr->nodeidx[tpos]   = k;
	  tr->pos[tpos]       = i--;
	}
      else if(sc  == mmx[i+1][k]   + hmm->tsc[CTME][k])    /* M->EL->M with 0 self loops in EL */
	{
	  tr->statetype[tpos] = CSTM;
	  tr->nodeidx[tpos]   = k--;
	  tr->pos[tpos]       = i+1; /* special case, we decremented i prematurely b/c we 
				      * had no way of knowing it was our last visit to EL, before
				      * we went to M (since we're working backwards this is actually
				      * the first visit to EL). 
				      */
	}
      else Die("traceback failed");
      break;

    default:
      Die("traceback failed");

    } /* end switch over statetype[tpos-1] */
    
    tpos++;
    if (tpos == curralloc) 
      {				/* grow trace if necessary  */
	curralloc += (j0-i0+1);
	CP9ReallocTrace(tr, curralloc);
      }

  } /* end traceback, at S state; tpos == tlen now */
  tr->tlen = tpos;
  CP9ReverseTrace(tr);

  *ret_tr = tr;
}

/* Function: CP9ReverseTrace()
 * Date:     EPN, Wed May 30 17:52:18 2007
 *           identical to SRE's P7ReverseTrace() from HMMER 2.3.2
 *
 * Purpose:  Reverse the arrays in a traceback structure.
 *           Tracebacks from Forward() and Viterbi() are
 *           collected backwards, and call this function
 *           when they're done.
 *           
 *           It's possible to reverse the arrays in place
 *           more efficiently; but the realloc/copy strategy
 *           has the advantage of reallocating the trace
 *           into the right size of memory. (Tracebacks
 *           overallocate.)
 *           
 * Args:     tr - the traceback to reverse. tr->tlen must be set.
 *                
 * Return:   (void)
 *           tr is modified.
 */                
void
CP9ReverseTrace(CP9trace_t *tr)
{
  char  *statetype;
  int   *nodeidx;
  int   *pos;
  int    opos, npos;

  /* Allocate
   */
  statetype = MallocOrDie (sizeof(char)* tr->tlen);
  nodeidx   = MallocOrDie (sizeof(int) * tr->tlen);
  pos       = MallocOrDie (sizeof(int) * tr->tlen);
  
  /* Reverse the trace.
   */
  for (opos = tr->tlen-1, npos = 0; npos < tr->tlen; npos++, opos--)
    {
      statetype[npos] = tr->statetype[opos];
      nodeidx[npos]   = tr->nodeidx[opos];
      pos[npos]       = tr->pos[opos];
    }

  /* Swap old, new arrays.
   */
  free(tr->statetype);
  free(tr->nodeidx);
  free(tr->pos);
  tr->statetype = statetype;
  tr->nodeidx   = nodeidx;
  tr->pos       = pos;
}



/* Function: CP9Traces2Alignment()
 *           based on SRE's P7Traces2Alignment() from HMMER 2.3.2
 *
 * Purpose:  Convert an array of traceback structures for a set
 *           of sequences into a new multiple alignment. Modified
 *           from HMMER to account for possible EL local-end 
 *           insertions (which don't exist in P7). Including EL
 *           insertions requires an emit map from the CM.
 *           
 *           Insertions/ELs are put into lower case and 
 *           are not aligned; instead, Nterm is right-justified,
 *           Cterm is left-justified, and internal insertions
 *           are split in half and the halves are justified in
 *           each direction (the objective being to increase
 *           the chances of getting insertions aligned well enough
 *           for them to become a match). SAM gap char conventions
 *           are used: - in match columns, . in insert columns
 * 
 * Args:     cm         - the CM the CP9 was built from, needed to get emitmap,
 *                        so we know where to put EL transitions
 *           sq         - sequences 
 *           nseq       - number of sequences
 *           tr         - array of tracebacks
 *           do_full    - TRUE to always include all match columns in alignment
 *
 * Return:   MSA structure; NULL on failure.
 *           Caller responsible for freeing msa with MSAFree(msa);
 */          
MSA *
CP9Traces2Alignment(CM_t *cm, ESL_SQ **sq, float *wgt, int nseq, 
		    CP9trace_t **tr, int do_full)
{
  /* Contract checks */
  if(cm->cp9 == NULL)
    esl_fatal("ERROR in CP9Traces2Alignment, cm->cp9 is NULL.\n");
  if(cm->cp9map == NULL)
    esl_fatal("ERROR in CP9Traces2Alignment, cm->cp9map is NULL.\n");
  if(!(cm->flags & CM_CP9))
     esl_fatal("ERROR in CP9Traces2Alignment, CM_CP9 flag is down.");

  MSA   *msa;                   /* RETURN: new alignment */
  int    idx;                   /* counter for sequences */
  int    alen;                  /* width of alignment */
  int   *maxins;                /* array of max inserts between aligned columns */
  int   *maxels;                /* array of max ELs emissions between aligned columns */
  int   *matmap;                /* matmap[k] = apos of match k [1..M] */
  int    nins;                  /* counter for inserts */
  int    cpos;                  /* HMM node, consensus position */
  int    apos;                  /* position in aligned sequence (0..alen-1)*/
  int    rpos;                  /* position in raw digital sequence (1..L)*/
  int    tpos;                  /* position counter in traceback */
  int    epos;                  /* position ctr for EL insertions */
  int    statetype;		/* type of current state, e.g. STM */
  CMEmitMap_t *emap;            /* consensus emit map for the CM */
  int         *imap;            /* first apos for an insert following a cpos */
  int         *elmap;           /* first apos for an EL following a cpos */
  int         *matuse;          /* TRUE if we need a cpos in mult alignment */
  int         *eluse;           /* TRUE if we have an EL after cpos in alignment */
  int        **eposmap;         /* [seq idx][CP9 node idx] where each EL should emit to */
  int         *iuse;            /* TRUE if we have an I after cpos in alignment */
  CMConsensus_t *con;           /* consensus information for the CM */
  int          next_match;      /* used for filling eposmap */
  int          c;               /* counter over possible EL froms */
  emap = CreateEmitMap(cm);

  /* Here's the problem. We want to align the match states in columns,
   * but some sequences have inserted symbols in them; we need some
   * sort of overall knowledge of where the inserts are and how long
   * they are in order to create the alignment. 
   * 
   * Here's our trick. maxins[] and maxels[] are 0..hmm->M arrays; 
   * maxins[i] stores the maximum number of times insert substate 
   * i was used. maxels[i] stores the max number of times an EL insertion
   * occurs after insert substate i. maxins[i] + maxels[i] 
   * is the maximum number of gaps to insert between canonical 
   * column i and i+1.  maxins[0], maxels[0] is the N-term tail; 
   * maxins[M], maxels[0] is the C-term tail.
   */
  matuse  = (int *) MallocOrDie(sizeof(int)  * (emap->clen+1));   
  eluse   = (int *) MallocOrDie(sizeof(int)  * (emap->clen+1));   
  iuse    = (int *) MallocOrDie(sizeof(int)  * (emap->clen+1));   
  maxins  = (int *) MallocOrDie(sizeof(int) * (emap->clen+1));
  maxels  = (int *) MallocOrDie(sizeof(int) * (emap->clen+1));
  matmap  = (int *) MallocOrDie(sizeof(int) * (emap->clen+1));
  imap    = (int *) MallocOrDie(sizeof(int) * (emap->clen+1));   
  elmap   = (int *) MallocOrDie(sizeof(int) * (emap->clen+1));   
  eposmap   = (int **) MallocOrDie(sizeof(int *) * nseq); 
  /* eposmap is 2D b/c different traces can have different epos
   * (position where EL inserts) for the same EL state, for example:
   * an EL state for node 9 may reconnect at node 25 in one parse
   * and node 50 in another if there's a CM MATL node with subtree
   * lpos=9 rpos=50, and a CM BEGL node with subtree lpos=9 rpos=25,
   * i.e. there are 2 CM EL states being mirrored by 1 HMM EL state. 
   */

  for (cpos = 0; cpos <= emap->clen; cpos++)
    {
      if(!do_full || cpos == 0)
	matuse[cpos] = 0;
      else
	matuse[cpos] = 1;
      maxins[cpos] = maxels[cpos] = 0;
      iuse[cpos] = eluse[cpos] = imap[cpos] = elmap[cpos] = 0;
    }
  
  /* Look at all the traces; find maximum length of
   * insert needed at each of the clen+1 possible gap
   * points. (There are two types of insert, I/EL)
   * Also find whether we don't need some of the match
   * (consensus) columns.
   */
  for (idx = 0; idx < nseq; idx++) 
    {
      eposmap[idx] = (int *) MallocOrDie(sizeof(int) * (emap->clen+1));   
      for (cpos = 0; cpos <= emap->clen; cpos++) 
	{
	  iuse[cpos] = eluse[cpos] = 0;
	  eposmap[idx][cpos] = -1;
	}      

      /* Determine the eposmap, the cpos EL's go into for each cpos for each seq.
       * This depends on the first match state entered after each EL, so we go bottom up */
      next_match = -1;
      for (tpos = tr[idx]->tlen - 1; tpos >= 0; tpos--) 
	{
	  statetype = tr[idx]->statetype[tpos]; /* just for clarity */
	  cpos      = tr[idx]->nodeidx[tpos];      
	  if(statetype == CSTM) next_match = cpos;
	  if(statetype == CSTE) next_match = cm->cp9->M+1;
	  if(statetype == CSTEL) eposmap[idx][cpos] = next_match; /* this will be overwritten below */
	}
      for (cpos = 0; cpos <= emap->clen; cpos++) 
	{
	  if(eposmap[idx][cpos] != -1)
	    {
	      /*printf("cpos: %d eposmap[idx][cpos]: %d ct: %d\n", cpos, eposmap[idx][cpos], cm->cp9->el_from_ct[eposmap[idx][cpos]]);*/
	      /* determine the epos based on the CM emit map and cm->cp9->el* data structures */
	      for(c = 0; c < cm->cp9->el_from_ct[eposmap[idx][cpos]]; c++)
		{
		  if(cm->cp9->el_from_idx[eposmap[idx][cpos]][c] == cpos)
		    {
		      eposmap[idx][cpos] = emap->epos[cm->cp9->el_from_cmnd[eposmap[idx][cpos]][c]];
		      break;
		    }
		  if(c == (cm->cp9->el_from_ct[eposmap[idx][cpos]] - 1))
		    esl_fatal("Couldn't determine epos for cpos: %d\n", cpos);
		}
	    }
	}

      for (tpos = 0; tpos < tr[idx]->tlen; tpos++)
	{
	  cpos = tr[idx]->nodeidx[tpos];
	  switch (tr[idx]->statetype[tpos]) {
	    case CSTI: iuse[cpos]++; break;
	  case CSTM: matuse[tr[idx]->nodeidx[tpos]] = 1; break;
	  case CSTEL: 
	    eluse[eposmap[idx][cpos]]++; 
	    break;
	  case CSTD:
	  case CSTE:
	  case CSTB:
	    break;
	  default:
	    Die("Traces2Alignment reports unrecognized statetype %c", 
		CP9Statetype(tr[idx]->statetype[tpos]));
	  }
	} /* end looking at trace i */
      for (cpos = 0; cpos <= emap->clen; cpos++) 
	{
	  if (iuse[cpos]  > maxins[cpos]) maxins[cpos]  = iuse[cpos];
	  if (eluse[cpos] > maxels[cpos]) maxels[cpos]  = eluse[cpos]-1; /* EL only emits on self loops */
	}
    } /* end calculating lengths used by all traces */

  /***********************************************
   * Construct the alignment
   ***********************************************/
  
  /* Now we can calculate the total length of the multiple alignment, alen;
   * and the maps imap and  elmap that turn a cpos into an apos
   * in the multiple alignment: e.g. for an insert that follows consensus 
   * position cpos, put it at or after apos = imap[cpos] in aseq[][].
   */
  
  matmap[0] = -1; /* M_0 is B state, non-emitter */
  alen = 0;
  for (cpos = 0; cpos <= emap->clen; cpos++)
    {
      if (matuse[cpos]) 
	{
	  matmap[cpos] = alen; 
	  alen++;
	} 
      else 
	matmap[cpos] = -1;
      
      imap[cpos]  = alen;
      alen += maxins[cpos];
      elmap[cpos] = alen; 
      alen += maxels[cpos];
    }
                                /* allocation for new alignment */
  msa = MSAAlloc(nseq, alen);

  for (idx = 0; idx < nseq; idx++) 
    {
      for (cpos = 0; cpos <= emap->clen; cpos++)
	{
	  iuse[cpos] = eluse[cpos] = 0;
	}
      /* blank an aseq */
      for (apos = 0; apos < alen; apos++)
	msa->aseq[idx][apos] = '.';
      for (cpos = 0; cpos <= emap->clen; cpos++)
	if (matmap[cpos] != -1) msa->aseq[idx][matmap[cpos]] = '-';
      msa->aseq[idx][alen] = '\0';

      /* align the sequence */
      apos = 0;
      for (tpos = 0; tpos < tr[idx]->tlen; tpos++) 
	{
	  statetype = tr[idx]->statetype[tpos]; /* just for clarity */
	  rpos      = tr[idx]->pos[tpos]; 
	  cpos      = tr[idx]->nodeidx[tpos];
	  
	  if (statetype == CSTM) 
	    {
	      apos = matmap[cpos];
	      msa->aseq[idx][apos] = Alphabet[(int) sq[idx]->dsq[rpos]];
	    }
	  else if (statetype == CSTD) 
	    apos = matmap[cpos]+1;	/* need for handling D->I; xref STL6/p.117 */
	  else if (statetype == CSTI) 
	    {
	      apos = imap[cpos] + iuse[cpos];
	      msa->aseq[idx][apos] = (char) tolower((int) Alphabet[(int) sq[idx]->dsq[rpos]]);
	      iuse[cpos]++;
	    }
	  else if (statetype == CSTEL) 
	    {
	      /*printf("CSTEL cpos: %d rpos: %d epos: %d\n", cpos, rpos);*/
	      epos = eposmap[idx][cpos];
	      if(tr[idx]->statetype[tpos-1] == CSTEL) /* we don't emit on first EL visit */
		{
		  apos = elmap[epos] + eluse[epos];
		  msa->aseq[idx][apos] = (char) tolower((int) Alphabet[(int) sq[idx]->dsq[rpos]]);
		  eluse[epos]++;
		}
	    }
	  else if (statetype == CSTE)
	    apos = matmap[emap->clen]+1;	/* set position for C-term tail */
	}
      /* N-terminal extension is right-justified.
       * Internal inserts are split in half, and C-term is right-justified.
       * C-terminal extension remains left-justified.
       */
      rightjustify(msa->aseq[idx], maxins[0]);
      
      for (cpos = 1; cpos < emap->clen; cpos++) 
	{
	  if (maxins[cpos] > 1) 
	    {
	      for (nins = 0, apos = matmap[cpos]+1; islower((int) (msa->aseq[idx][apos])); apos++)
		nins++;
	      nins /= 2;		/* split the insertion in half */
	      rightjustify(msa->aseq[idx]+matmap[cpos]+1+nins, maxins[cpos]-nins);
	    }
	}
    }
  /***********************************************
   * Build the rest of the MSA annotation.
   ***********************************************/
        
  msa->nseq = nseq;
  msa->alen = alen;
  msa->au   = MallocOrDie(sizeof(char) * (strlen(PACKAGE_VERSION)+10));
  sprintf(msa->au, "Infernal %s", PACKAGE_VERSION);
  /* copy sqinfo array and weights */
  for (idx = 0; idx < nseq; idx++)
    {
      msa->sqname[idx] = sre_strdup(sq[idx]->name, -1);
      msa->sqlen[idx]  = sq[idx]->n;
      /*if (sqinfo[idx].flags & SQINFO_ACC)
	MSASetSeqAccession(msa, idx, sqinfo[idx].acc);
	if (sqinfo[idx].flags & SQINFO_DESC)
	MSASetSeqDescription(msa, idx, sqinfo[idx].desc);

      if (sqinfo[idx].flags & SQINFO_SS) {
	if (msa->ss == NULL) msa->ss = MallocOrDie(sizeof(char *) * nseq);
	MakeAlignedString(msa->aseq[idx], alen, 
			  sqinfo[idx].ss, &(msa->ss[idx]));
      }
      if (sqinfo[idx].flags & SQINFO_SA) {
	if (msa->sa == NULL) msa->sa = MallocOrDie(sizeof(char *) * nseq);
	MakeAlignedString(msa->aseq[idx], alen, 
			  sqinfo[idx].sa, &(msa->sa[idx]));
			  }*/
      if (wgt == NULL) msa->wgt[idx] = 1.0;
      else             msa->wgt[idx] = wgt[idx];
    }

  /* Construct the secondary structure consensus line, msa->ss_cons:
   *       IL, IR are annotated as .
   *       EL is annotated as ~
   *       and match columns use the structure code.
   * Also the primary sequence consensus/reference coordinate system line,
   * msa->rf.
   */
  msa->ss_cons = MallocOrDie(sizeof(char) * (alen+1));
  msa->rf = MallocOrDie(sizeof(char) * (alen+1));
  con = CreateCMConsensus(cm, 3.0, 1.0);

  for (cpos = 0; cpos <= emap->clen; cpos++) 
    {
      if (matuse[cpos]) 
	{ /* CMConsensus is off-by-one right now, 0..clen-1 relative to cpos's 1..clen */
	  if (con->ct[cpos-1] != -1 && matuse[con->ct[cpos-1]+1] == 0) {
	    msa->ss_cons[matmap[cpos]] = '.';
	    msa->rf[matmap[cpos]]      = con->cseq[cpos-1];
	  } else {
	    msa->ss_cons[matmap[cpos]] = con->cstr[cpos-1];	
	    msa->rf[matmap[cpos]]      = con->cseq[cpos-1];
	  }
	}
      if (maxins[cpos] > 0) 
	for (apos = imap[cpos]; apos < imap[cpos] + maxins[cpos]; apos++)
	  {
	    msa->ss_cons[apos] = '.';
	    msa->rf[apos] = '.';
	  }
      if (maxels[cpos] > 0)
	{
	  for (apos = elmap[cpos]; apos < elmap[cpos] + maxels[cpos]; apos++)
	  {
	    msa->ss_cons[apos] = '~';
	    msa->rf[apos] = '~';
	  }
	}
    }
  msa->ss_cons[alen] = '\0';
  msa->rf[alen] = '\0';

  /* Free and return */
  FreeCMConsensus(con);
  FreeEmitMap(emap);
  free(matuse);
  free(eluse);
  free(iuse);
  free(maxins);
  free(maxels);
  free(matmap);
  free(imap);
  free(elmap);
  for(idx = 0; idx < nseq; idx++)
    free(eposmap[idx]);
  free(eposmap);
  return msa;
}
/* Function: rightjustify()
 * 
 * Purpose:  Given a gap-containing string of length n,
 *           pull all the non-gap characters as far as
 *           possible to the right, leaving gaps on the
 *           left side. Used to rearrange the positions
 *           of insertions in HMMER alignments.
 */
static void
rightjustify(char *s, int n)
{
  int npos;
  int opos;

  npos = n-1;
  opos = n-1;
  while (opos >= 0) {
    if (isgap(s[opos])) opos--;
    else                s[npos--]=s[opos--];  
  }
  while (npos >= 0) 
    s[npos--] = '.';
}


/*
 * Function: DuplicateCP9
 * Date:     EPN, Thu Jun 28 13:37:22 2007
 * Purpose:  Given a template cm 'src_cm' copy it's CP9 
 *           to the cm 'dest_cm'. dest_cm->cp9 and dest_cm->cp9map
 *           are alloc'ed here.
 *
 * Args:
 *           src_cm         the source CM, must have valid cp9
 *           dest_cm        the destination CM we're copying src_cm->cp9
 *                          to
 */
void
DuplicateCP9(CM_t *src_cm, CM_t *dest_cm)
{
  int       k,x;	          /* counter over nodes */

  /* Contract checks */
  if(!(src_cm->flags & CM_CP9))
    esl_fatal("ERROR in DuplicateCP9() src_cm CM_CP9 flag down.\n");

  CPlan9Renormalize(src_cm->cp9);
  CP9Logoddsify(src_cm->cp9);
  
  /* We can fill the map before we copy the CP9 */
  /* Allocate and initialize the cp9map */
  dest_cm->cp9map = AllocCP9Map(dest_cm);
  /* Map the CM states to CP9 states and nodes and vice versa */
  CP9_map_cm2hmm(dest_cm, dest_cm->cp9map, 0);

  /* Create the new model and copy everything over */
  dest_cm->cp9 = AllocCPlan9(dest_cm->cp9map->hmm_M);
  ZeroCPlan9(dest_cm->cp9);
  CPlan9SetNullModel(dest_cm->cp9, dest_cm->bg->f, 1.0); /* set p1 = 1.0 which corresponds to the CM */
  CPlan9InitEL(dest_cm, dest_cm->cp9); /* set up hmm->el_from_ct and hmm->el_from_idx data, which
					* explains how the EL states are connected in the HMM. */

  /* Copy the transitions and emission probs and scores */
  for(k = 0; k <= dest_cm->cp9->M; k++)
    {
      dest_cm->cp9->t[k][CTMM] = src_cm->cp9->t[k][CTMM];
      dest_cm->cp9->t[k][CTMI] = src_cm->cp9->t[k][CTMI];
      dest_cm->cp9->t[k][CTMD] = src_cm->cp9->t[k][CTMD];
      dest_cm->cp9->t[k][CTME] = src_cm->cp9->t[k][CTME];


      dest_cm->cp9->t[k][CTIM] = src_cm->cp9->t[k][CTIM];
      dest_cm->cp9->t[k][CTII] = src_cm->cp9->t[k][CTII];
      dest_cm->cp9->t[k][CTID] = src_cm->cp9->t[k][CTID];

      dest_cm->cp9->t[k][CTDM] = src_cm->cp9->t[k][CTDM];
      dest_cm->cp9->t[k][CTDI] = src_cm->cp9->t[k][CTDI];
      dest_cm->cp9->t[k][CTDD] = src_cm->cp9->t[k][CTDD];

      dest_cm->cp9->begin[k] = src_cm->cp9->begin[k];
      dest_cm->cp9->end[k] = src_cm->cp9->end[k];

      for(x = 0; x < Alphabet_size; x++)
	{
	  dest_cm->cp9->mat[k][x] = src_cm->cp9->mat[k][x];
	  dest_cm->cp9->ins[k][x] = src_cm->cp9->ins[k][x];
	}
    }

  dest_cm->cp9->p1        = src_cm->cp9->p1;
  dest_cm->cp9->el_self   = src_cm->cp9->el_self;
  dest_cm->cp9->el_selfsc = src_cm->cp9->el_selfsc;

  CPlan9Renormalize(dest_cm->cp9);/* shouldn't be nec */
  CP9Logoddsify(dest_cm->cp9); /* fill in all the integer log odds scores:
				* msc, isc, bsc, esc, tsc, the *sc_mem
				* pointers were set up in AllocCPlan9() */
  if(src_cm->config_opts & CM_CONFIG_ZEROINSERTS)
    CP9HackInsertScores(dest_cm->cp9);
    
  /*
    FILE *fp;
    fp = fopen("destcp9" ,"w");
    debug_print_cp9_params(fp, dest_cm->cp9);
    fclose(fp);
    fp = fopen("srccp9" ,"w");
    debug_print_cp9_params(fp, src_cm->cp9);
    fclose(fp);
  */

  dest_cm->cp9->flags = src_cm->cp9->flags;
}


/* Following functions for CPlan9 HMMs were deprecated 01.04.07,
 * we never use these aspects of a CP9 HMM.
 */
#if 0
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
#endif


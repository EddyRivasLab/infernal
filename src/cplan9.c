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
 * CM Plan 9 structre, but is not yet used.
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
			/* statistical parameters set to innocuous empty values */
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
  hmm->ins    = MallocOrDie ((M+2) *           sizeof(float *));
  hmm->t[0]   = MallocOrDie ((9*(M+1))     *       sizeof(float));
  hmm->mat[0] = MallocOrDie ((MAXABET*(M+1)) * sizeof(float));
  hmm->ins[0] = MallocOrDie ((MAXABET*(M+1)) *     sizeof(float));

  hmm->tsc     = MallocOrDie (9     *           sizeof(int *));
  hmm->msc     = MallocOrDie (MAXDEGEN   *       sizeof(int *));
  hmm->isc     = MallocOrDie (MAXDEGEN   *       sizeof(int *)); 
  hmm->tsc_mem = MallocOrDie ((9*(M+2))     *       sizeof(int));
  hmm->msc_mem = MallocOrDie ((MAXDEGEN*(M+1)) * sizeof(int));
  hmm->isc_mem = MallocOrDie ((MAXDEGEN*(M+2)) *     sizeof(int));

  hmm->tsc[0] = hmm->tsc_mem;
  hmm->msc[0] = hmm->msc_mem;
  hmm->isc[0] = hmm->isc_mem;

  /* note allocation strategy for important 2D arrays -- trying
   * to keep locality as much as possible, cache efficiency etc.
   */
  for (k = 1; k <= M; k++) {
    hmm->mat[k] = hmm->mat[0] + k * MAXABET;
    hmm->ins[k] = hmm->ins[0] + k * MAXABET;
    hmm->t[k]   = hmm->t[0]   + k * 9;
  }
  for (x = 1; x < MAXDEGEN; x++) {
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


void
FreeCPlan9(CP9_t *hmm)
{
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
  /*printf("N: %d | maxN: %d | M: %d | maxM: %d\n", N, mx->maxN, M, mx->maxM);*/

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

  /* Configure entry.
   */
  for (k = 1; k <= hmm->M; k++)
    hmm->begin[k] = 2. * (float) (hmm->M-k+1) / (float) hmm->M / (float) (hmm->M+1);

  hmm->begin[1] = (1. - pentry) * (1. - (hmm->t[0][CTMI] + hmm->t[0][CTMD]));
  FSet(hmm->begin+2, hmm->M-1, (pentry * (1.- (hmm->t[0][CTMI] + hmm->t[0][CTMD]))) / (float)(hmm->M-1));
  
  /* Configure exit.
   * hmm->end[hmm->M] = 1. - hmm->t[hmm->M][XTMI] 
   *
   * hmm->t[hmm->M][CTIM] is M_I->E transition
   * hmm->t[hmm->M][CTDM] is M_D->E transition
   * and we don't want to ignore these.
   */

  basep = pexit / (float) (hmm->M-1);
  for (k = 1; k < hmm->M; k++)
    hmm->end[k] = basep / (1. - basep * (float) (k-1));
  CPlan9RenormalizeExits(hmm, 1);
  hmm->flags       &= ~CPLAN9_HASBITS; /* reconfig invalidates log-odds scores */
}

/* Function: CPlan9GlobalConfig()
 * EPN 09.24.06
 * based on SRE's Plan7GlobalConfig() from HMMER's plan7.c
 * 
 * Purpose:  Set the alignment independent parameters of
 *           a CM Plan 9 model to (Needleman/Wunsch) configuration.
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
  /* No special (*x* states in Plan 7) states in CM Plan 9 */

  /* Configure entry.
   * Exactly 3 ways to start, B->M_1 (hmm->begin[1]), B->I_0 (hmm->t[0][CTMI]),
   *                      and B->D_1 (hmm->t[0][CTMD])
   */
  hmm->begin[1] = 1. - (hmm->t[0][CTMI] + hmm->t[0][CTMD]);
  FSet(hmm->begin+2, hmm->M-1, 0.);
  
  hmm->end[hmm->M] = 1. - hmm->t[hmm->M][CTMI];
  FSet(hmm->end+1, hmm->M-1, 0.);

  CPlan9RenormalizeExits(hmm, 1);
  hmm->flags       &= ~CPLAN9_HASBITS; /* reconfig invalidates log-odds scores */
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
	  d = FSum(hmm->t[k], 3);
	  /* FScale(hmm->t[k], 3, 1./(d + d*hmm->end[k])); */
	  FScale(hmm->t[k], 3, (1.-hmm->end[k])/d);
	}
    }
  return;
}

/* Function: CP9AllocTrace(), CP9ReallocTrace(), CP9FreeTrace()
 * 
 * Purpose:  allocation and freeing of traceback structures
 */
void
CP9AllocTrace(int tlen, struct cp9trace_s **ret_tr)
{
  struct cp9trace_s *tr;
  
  tr =            MallocOrDie (sizeof(struct cp9trace_s));
  tr->statetype = MallocOrDie (sizeof(char) * tlen);
  tr->nodeidx   = MallocOrDie (sizeof(int)  * tlen);
  tr->pos       = MallocOrDie (sizeof(int)  * tlen);
  *ret_tr = tr;
}
void
CP9ReallocTrace(struct cp9trace_s *tr, int tlen)
{
  tr->statetype = ReallocOrDie (tr->statetype, tlen * sizeof(char));
  tr->nodeidx   = ReallocOrDie (tr->nodeidx,   tlen * sizeof(int));
  tr->pos       = ReallocOrDie (tr->pos,       tlen * sizeof(int));
}
void 
CP9FreeTrace(struct cp9trace_s *tr)
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

      for(x = 0; x < 9; x++)
	{
	  sub_hmm->t[i][x]   = orig_hmm->t[orig_pos][x];
	  sub_hmm->tsc[x][i] = orig_hmm->tsc[x][orig_pos];
	}
      
    }

  /* Make the necessary modifications. */
  CP9_reconfig2sub(sub_hmm, spos, epos, 1, sub_hmm->M, orig_phi);

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
      hmm->t[spos_nd-1][CTDM] = 0.; /* D_0 doesn't exist */
      hmm->t[spos_nd-1][CTDI] = 0.; /* D_0 doesn't exist */
      hmm->t[spos_nd-1][CTDD] = 0.; /* D_0 doesn't exist */
      
      hmm->bsc[spos_nd]       = Prob2Score(hmm->begin[1], 1.0);

      hmm->tsc[CTMM][spos_nd-1] = -INFTY; /* probability of going from B(M_0) to M_1 is begin[1] */
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
      hmm->t[epos_nd][CTMD]  = 0.; /* D_M+1 doesn't exist */
      hmm->t[epos_nd][CTDD]  = 0.; /* D_M+1 doesn't exist */
      hmm->t[epos_nd][CTID]  = 0.; /* D_M+1 doesn't exist */
      
      hmm->esc[epos_nd]       = Prob2Score(hmm->end[epos_nd], 1.0);
      hmm->tsc[CTDM][epos_nd] = Prob2Score(hmm->t[epos_nd][CTDM], 1.0);
      hmm->tsc[CTIM][epos_nd] = Prob2Score(hmm->t[epos_nd][CTIM], 1.0);
      hmm->tsc[CTMM][epos_nd] = -INFTY; /* M->E is actually end[M] */
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

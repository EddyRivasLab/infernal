/* CP9_MX functions: dynamic programming matrix for CP9 HMMs
 * 
 * EPN, Wed Nov 28 05:11:51 2007
 * SVN $Id$
 */

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"

#include "funcs.h"
#include "structs.h"

/*****************************************************************
 *   1. CP9_MX data structure functions,
 *      matrix of integer log odd scores for CP9 HMM alignment/search
 *****************************************************************/

/* Function: CreateCP9Matrix()
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
 * Args:     N     - N+1 rows are allocated, usually N == 1 for 
 *                   scanning in memory efficient mode, or N == L, length of sequence.  
 *           M     - size of model in nodes
 *                 
 * Return:   mx
 *           mx is allocated here. Caller frees with FreeCP9Matrix(mx).
 */
CP9_MX *
CreateCP9Matrix(int N, int M)
{
  int status;
  CP9_MX *mx;
  int i;

  ESL_ALLOC(mx,      sizeof(CP9_MX));
  ESL_ALLOC(mx->mmx, sizeof(int *) * (N+1));
  ESL_ALLOC(mx->imx, sizeof(int *) * (N+1));
  ESL_ALLOC(mx->dmx, sizeof(int *) * (N+1));
  ESL_ALLOC(mx->elmx,sizeof(int *) * (N+1)); 
  /* slightly wasteful, some nodes can't go to EL (for ex: right half of MATPs) */
  ESL_ALLOC(mx->erow,    sizeof(int) * (N+1));
  ESL_ALLOC(mx->mmx_mem, sizeof(int) * ((N+1)*(M+1)));
  ESL_ALLOC(mx->imx_mem, sizeof(int) * ((N+1)*(M+1)));
  ESL_ALLOC(mx->dmx_mem, sizeof(int) * ((N+1)*(M+1)));
  ESL_ALLOC(mx->elmx_mem,sizeof(int) * ((N+1)*(M+1)));

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

  mx->M = M;
  mx->rows = N;
  mx->size_Mb =  (float) sizeof(CP9_MX);
  mx->size_Mb += (float) (sizeof(int *) * (mx->rows+1) * 4); /* mx->*mx ptrs */
  mx->size_Mb += (float) (sizeof(int)   * (mx->rows+1) * (M+1) * 4); /* mx->*mx_mem */
  mx->size_Mb += (float) (sizeof(int)   * (mx->rows+1));             /* mx->erow */
  mx->size_Mb /= 1000000.;

  return mx;

 ERROR:
  cm_Fail("Memory allocation error.");
  return NULL; /* never reached */
}


/* Function: FreeCP9Matrix()
 * based on  FreePlan7Matrix() <-- this function's comments below  
 * Purpose:  Free a dynamic programming matrix allocated by CreatePlan7Matrix().
 * 
 * Return:   (void)
 */
void
FreeCP9Matrix(CP9_MX *mx)
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

/* 
 * Function: GrowCP9Matrix()
 *
 * Purpose:  Reallocate a CP9 dp matrix, if necessary, for seq for
 *           length N, or 2 rows (if we're scanning in memory 
 *           efficient mode, in this case N == 1, nrows = N+1).
 * 
 *           Note: unlike HMMER, M never changes, so we only have
 *           to worry about increasing the number of rows if nec.
 *           
 *           Returns individual ptrs to the four matrix components
 *           as a convenience.
 *           
 *           This function allocates the requested matrix regardless
 *           of it's size.
 * 
 * Args:     mx    - an already allocated matrix to grow.
 *           N     - seq length to allocate for; N+1 rows
 *           M     - size of model, contract enforces this must == mx->M
 *           mmx, imx, dmx, elmx, erow 
 *                 - RETURN: ptrs to four mx components as a convenience
 *                   
 * Return:   eslOK on success, eslEINCOMPAT if contract is violated,
 *           mx is (re)allocated here.
 */
int
GrowCP9Matrix(CP9_MX *mx, char *errbuf, int N, int M, int ***mmx, int ***imx, int ***dmx, int ***elmx, int **erow)
{
  int status;
  void *p;
  int i;

  if(mx->M != M) ESL_FAIL(eslEINCOMPAT, errbuf, "GrowCP9Matrix(), mx->M: %d != M passed in: %d\n", mx->M, M);
  if(N < 0)      ESL_FAIL(eslEINCOMPAT, errbuf, "GrowCP9Matrix(), N: %d < 0\n", N);

  if (N <= mx->rows) goto DONE;

  /* If we didn't goto DONE, then N > mx->rows, we need more space */
  ESL_RALLOC(mx->mmx,  p, sizeof(int *) * (N+1));
  ESL_RALLOC(mx->imx,  p, sizeof(int *) * (N+1));
  ESL_RALLOC(mx->dmx,  p, sizeof(int *) * (N+1));
  ESL_RALLOC(mx->elmx, p, sizeof(int *) * (N+1)); 
  ESL_RALLOC(mx->erow, p, sizeof(int)   * (N+1));
  ESL_RALLOC(mx->mmx_mem,  p, sizeof(int) * ((N+1)*(M+1)));
  ESL_RALLOC(mx->imx_mem,  p, sizeof(int) * ((N+1)*(M+1)));
  ESL_RALLOC(mx->dmx_mem,  p, sizeof(int) * ((N+1)*(M+1)));
  ESL_RALLOC(mx->elmx_mem, p, sizeof(int) * ((N+1)*(M+1)));
  
  mx->mmx[0]  = mx->mmx_mem;
  mx->imx[0]  = mx->imx_mem;
  mx->dmx[0]  = mx->dmx_mem;
  mx->elmx[0] = mx->elmx_mem;
  for (i = 1; i <= N; i++) {
    mx->mmx[i] = mx->mmx[0] + (i*(M+1));
    mx->imx[i] = mx->imx[0] + (i*(M+1));
    mx->dmx[i] = mx->dmx[0] + (i*(M+1));
    mx->elmx[i]= mx->elmx[0]+ (i*(M+1));
  }
  mx->size_Mb =  (float) sizeof(CP9_MX);
  mx->size_Mb += (float) (sizeof(int *) * (N+1) * 4); /* mx->*mx ptrs */
  mx->size_Mb += (float) (sizeof(int)   * (N+1) * (M+1) * 4); /* mx->*mx_mem */
  mx->size_Mb += (float) (sizeof(int)   * (N+1));             /* mx->erow */
  mx->size_Mb /= 1000000.;

 DONE:
  mx->rows = N;
  if (mmx != NULL) *mmx = mx->mmx;
  if (imx != NULL) *imx = mx->imx;
  if (dmx != NULL) *dmx = mx->dmx;
  if (elmx!= NULL) *elmx= mx->elmx;
  if (erow != NULL) *erow = mx->erow;
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, ("ResizeCP9Matrix(), memory reallocation error."));
}


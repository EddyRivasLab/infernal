/* The MSV filter implementation; SSE version.
 * 
 * A "filter" is a one-row, O(M), DP implementation that calculates
 * an approximated nat score (i.e. in limited precision - uchar, for 
 * example) and may have limited numeric range. It will return 
 * <eslERANGE> if its numeric range is exceeded, in which case the caller
 * will have to obtain the score by another (probably slower) method.
 * 
 * Contents:
 *   1. p7_MSVFilter() implementation.
 *   2. Benchmark driver.
 *   3. Unit tests
 *   4. Test driver
 *   5. Example
 *   6. Copyright and license information
 * 
 * SRE, Sun Nov 25 11:26:48 2007 [Casa de Gatos]
 * SVN $Id: impl_sse.c 2509 2008-07-30 14:45:52Z eddys $
 */
#include "config.h"
#include "p7_config.h"

#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"
#include "esl_sse.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_sse.h"

#include "funcs.h"
#include "structs.h"

/*****************************************************************
 * 1. The p7_MSVFilter() DP implementation.
 *****************************************************************/

/* Function:  p7_MSVFilter()
 * Synopsis:  Calculates MSV score, vewy vewy fast, in limited precision.
 * Incept:    SRE, Wed Dec 26 15:12:25 2007 [Janelia]
 *
 * Purpose:   Calculates an approximation of the MSV score for sequence
 *            <dsq> of length <L> residues, using optimized profile <om>,
 *            and a preallocated one-row DP matrix <ox>. Return the 
 *            estimated MSV score (in nats) in <ret_sc>.
 *            
 *            Score may overflow (and will, on high-scoring
 *            sequences), but will not underflow. 
 *            
 *            The model may be in any mode, because only its match
 *            emission scores will be used. The MSV filter inherently
 *            assumes a multihit local mode, and uses its own special
 *            state transition scores, not the scores in the profile.
 *
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of dsq in residues          
 *            om      - optimized profile
 *            ox      - DP matrix
 *            ret_sc  - RETURN: MSV score (in nats)          
 *                      
 * Note:      We misuse the matrix <ox> here, using only a third of the
 *            first dp row, accessing it as <dp[0..Q-1]> rather than
 *            in triplets via <{MDI}MX(q)> macros, since we only need
 *            to store M state values. We know that if <ox> was big
 *            enough for normal DP calculations, it must be big enough
 *            to hold the MSVFilter calculation.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if the score overflows the limited range; in
 *            this case, this is a high-scoring hit.
 *
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
my_p7_MSVFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, P7_GMX *gx, float *ret_sc)
{
  register __m128i mpv;            /* previous row values                                       */
  register __m128i xEv;		   /* E state: keeps max for Mk->E as we go                     */
  register __m128i xBv;		   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register __m128i sv;		   /* temp storage of 1 curr row value in progress              */
  register __m128i biasv;	   /* emission bias in a vector                                 */
  uint8_t  xE, xB, xC;             /* special states' scores                                    */
  int i;			   /* counter over sequence positions 1..L                      */
  int q;			   /* counter over vectors 0..nq-1                              */
  int Q        = p7O_NQU(om->M);   /* segment length: # of vectors                              */
  __m128i *dp  = ox->dpu[0];	   /* we're going to use dp[0][0..q..Q-1], not {MDI}MX(q) macros*/
  __m128i *rsc;			   /* will point at om->ru[x] for residue x[i]                  */

  /* Check that the DP matrix is ok for us. */
  if (Q > ox->allocQ16)  ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
  ox->M   = om->M;

  /* Initialization. In offset unsigned arithmetic, -infinity is 0, and 0 is om->base.
   */
  biasv = _mm_set1_epi8((int8_t) om->bias); /* yes, you can set1() an unsigned char vector this way */
  for (q = 0; q < Q; q++)
    dp[q] = _mm_setzero_si128();
  xB   = om->base - om->tjb;                /* remember, all values are costs to be subtracted. */
  xC   = 0;

  p7_omx_CopyMSVRow2gmx(ox, om, gx, 0, xE, 0, xC, xB, xC);
#if p7_DEBUGGING
  if (ox->debugging) p7_omx_DumpMSVRow(ox, 0, 0, 0, xC, xB, xC);
#endif

  for (i = 1; i <= L; i++)
    {
      rsc = om->rm[dsq[i]];
      xEv = _mm_setzero_si128();      
      xBv = _mm_set1_epi8((int8_t) (xB - om->tbm));

      /* Right shifts by 1 byte. 4,8,12,x becomes x,4,8,12. 
       * Because ia32 is littlendian, this means a left bit shift.
       * Zeros shift on automatically, which is our -infinity.
       */
      mpv = _mm_slli_si128(dp[Q-1], 1);   
      for (q = 0; q < Q; q++)
	{
	  /* Calculate new MMXo(i,q); don't store it yet, hold it in sv. */
	  sv   = _mm_max_epu8(mpv, xBv);
	  sv   = _mm_adds_epu8(sv, biasv);     
	  sv   = _mm_subs_epu8(sv, *rsc);   rsc++;
	  xEv  = _mm_max_epu8(xEv, sv);	

	  mpv   = dp[q];   	  /* Load {MDI}(i-1,q) into mpv */
	  dp[q] = sv;       	  /* Do delayed store of M(i,q) now that memory is usable */
	}	  

      /* Now the "special" states, which start from Mk->E (->C, ->J->B) */
      xE = esl_sse_hmax_epu8(xEv);
      if (xE >= 255 - om->bias) { *ret_sc = eslINFINITY; return eslERANGE; }	/* immediately detect overflow */

      xC = ESL_MAX(xC,        xE  - om->tec);
      xB = ESL_MAX(om->base,  xC) - om->tjb;
#if 1
      p7_omx_CopyMSVRow2gmx(ox, om, gx, i, xE, 0, xC, xB, xC); 
#endif
#if p7_DEBUGGING
      if (ox->debugging) p7_omx_DumpMSVRow(ox, i, xE, 0, xC, xB, xC);   
#endif
    } /* end loop over sequence residues 1..L */

  if (ox->debugging) { 
    p7_gmx_Dump(stdout, gx);
    ESL_DMATRIX *D;
    double min;
    double max;
    FILE *hfp;
    p7_gmx_Match2DMatrix(gx, FALSE, &D, &min, &max);
    hfp = fopen("cur.ps", "w");
    dmx_Visualize(hfp, D, min, max);
    fclose(hfp);
    esl_dmatrix_Destroy(D);
  }

  /* finally C->T, and add our missing precision on the NN,CC,JJ back */
  *ret_sc = ((float) (xC - om->tjb) - (float) om->base);
  *ret_sc /= om->scale;
  *ret_sc -= 3.0; /* that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ */

  return eslOK;
}
/*------------------ end, p7_MSVFilter() ------------------------*/

           
/* Function:  p7_omx_CopyMSVRow2gmx()
 * Based on p7_omx_DumpMSVRow()
 * Synopsis:  Copy one row of the MSV uchar version of a DP matrix to 
 *            a generic matrix in the correct position.
 * 
 * Incept:    EPN, Wed Aug  6 10:16:25 2008
 *
 * Purpose:   Copy current row of uchar part of DP matrix <ox> to generic
 *            matrix.
 * 
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure. * Args:      
 */
int
p7_omx_CopyMSVRow2gmx(P7_OMX *ox, const P7_OPROFILE *om, P7_GMX *gx, int rowi, uint8_t xE, uint8_t xN, uint8_t xJ, uint8_t xB, uint8_t xC)
{
  __m128i *dp = ox->dpu[0];	
  int      M  = ox->M;
  int      Q  = p7O_NQU(M);
  uint8_t *v  = NULL;		/* array of unstriped scores  */
  int      q,z,k;
  union { __m128i v; uint8_t i[16]; } tmp;
  int      status;

  ESL_ALLOC(v, sizeof(unsigned char) * ((Q*16)+1));
  v[0] = 0;

  /* Unpack and unstripe, then copy M's. */
  for (q = 0; q < Q; q++) {
    tmp.v = dp[q];
    for (z = 0; z < 16; z++) v[q+Q*z+1] = tmp.i[z];
  }
  for (k = 0; k <= M; k++) { 
    gx->dp[rowi][k * p7G_NSCELLS + p7G_M] = ((float) v[k] - (float) om->base) / om->scale;
    gx->dp[rowi][k * p7G_NSCELLS + p7G_I] = 0.;
    gx->dp[rowi][k * p7G_NSCELLS + p7G_D] = 0.;
  }
  /* The specials */
  gx->xmx[rowi * p7G_NXCELLS + 0] = ((float) xE - (float) om->base) / om->scale;
  gx->xmx[rowi * p7G_NXCELLS + 1] = ((float) xN - (float) om->base) / om->scale;
  gx->xmx[rowi * p7G_NXCELLS + 2] = ((float) xJ - (float) om->base) / om->scale;
  gx->xmx[rowi * p7G_NXCELLS + 3] = ((float) xB - (float) om->base) / om->scale;
  gx->xmx[rowi * p7G_NXCELLS + 4] = ((float) xC - (float) om->base) / om->scale;

  free(v);
  return eslOK;

ERROR:
  free(v);
  return status;
}


/* Function:  p7_gmx_Match2DMatrix()
 * Synopsis:  Copy the dp match cells of a generic matrix 
 *            to a ESL_DMATRIX, for visualization with esl_dmx_Visualize()
 * 
 * Incept:    SRE, Fri Jul 13 09:56:04 2007 [Janelia]
 *
 * Purpose:   Dump matrix <gx> to stream <fp> for diagnostics.
 */
int
p7_gmx_Match2DMatrix(P7_GMX *gx, int do_diff, ESL_DMATRIX **ret_D, double *ret_min, double *ret_max)
{
  int i, k;
  ESL_DMATRIX *D;
  double min =  eslINFINITY;
  double max = -eslINFINITY;
  float sc;

  D = esl_dmatrix_Create(gx->M+1, gx->L);
  /* fill k == 0 row, the X matrix E state (logically, the begin state) scores */
  for (i = 1; i <= gx->L; i++) { 
    D->mx[0][(i-1)] = gx->xmx[i * p7G_NXCELLS + 0];
  }

  for (i = 1; i <= gx->L; i++) { 
    for (k = 1; k <= gx->M; k++) { 
      sc = gx->dp[i][k * p7G_NSCELLS + p7G_M];
      if(do_diff) { 
	if((i < 2) || (k < 2)) D->mx[k][(i-1)] = 0.;
	else                   D->mx[k][(i-1)] = sc - ESL_MAX(D->mx[k][(i-2)], D->mx[0][(i-2)]);
      }
      else { 
	D->mx[k][(i-1)] = sc;
      }
      min = ESL_MIN(min, D->mx[k][(i-1)]);
      max = ESL_MAX(max, D->mx[k][(i-1)]);
    }
  }

  *ret_D = D;
  *ret_min = min;
  *ret_max = max;
  return eslOK;
}

/****************************************************************
 * Stolen from hmmer/h3/heatmap.c SVN revision 2171
 * as dmx_Visualize. Then modified so that the full 
 * matrix is printed (not half split diagonally).
 */
/* my_dmx_Visualize()
 * Incept:    SRE, Wed Jan 24 11:58:21 2007 [Janelia]
 *
 * Purpose:   
 *            
 *            Color scheme roughly follows Tufte, Envisioning
 *            Information, p.91, where he shows a beautiful
 *            bathymetric chart. The CMYK values conjoin two
 *            recommendations from ColorBrewer (Cindy Brewer
 *            and Mark Harrower) 
 *            [http://www.personal.psu.edu/cab38/ColorBrewer/ColorBrewer.html],
 *            specifically the 9-class sequential2 Blues and
 *            9-class sequential YlOrBr.
 * 
 *            Might eventually become part of Easel, once mature?
 *           
 * Note:      Binning rules basically follow same convention as
 *            esl_histogram. nb = xmax-xmin/w, so w = xmax-xmin/nb; 
 *            picking bin is (int) ceil((x - xmin)/w) - 1. (xref
 *            esl_histogram_Score2Bin()). This makes bin b contain
 *            values bw+min < x <= (b+1)w+min. (Which means that 
 *            min itself falls in bin -1, whoops - but we catch
 *            all bin<0 and bin>=nshades and put them in the extremes.
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */
int
my_dmx_Visualize(FILE *fp, ESL_DMATRIX *D, double min, double max, double min2fill)
{
   int    nshades   = 18;
   double cyan[]    = { 1.00, 1.00, 0.90, 0.75, 0.57, 0.38, 0.24, 0.13, 0.03,
			0.00, 0.00, 0.00, 0.00, 0.00, 0.07, 0.20, 0.40, 0.60};
   double magenta[] = { 0.55, 0.45, 0.34, 0.22, 0.14, 0.08, 0.06, 0.03, 0.01,
			0.00, 0.03, 0.11, 0.23, 0.40, 0.55, 0.67, 0.75, 0.80};
   double yellow[]  = { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
			0.10, 0.25, 0.40, 0.65, 0.80, 0.90, 1.00, 1.00, 1.00};
   double black[]   = { 0.30, 0.07, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
			0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};
   double w;			
   int    i,j;
   int    bin;
   int    boxsize;		/* box size in points */
   int    xcoord, ycoord;	/* postscript coords in points */
   int    leftmargin, rightmargin;
   int    bottommargin, topmargin;
   float  fboxsize;		/* box size in fractional points */

   /* Set some defaults that might become arguments later.
    */
   leftmargin   = rightmargin = 20;
   bottommargin = topmargin   = 20;

   /* Determine some working parameters 
    */
   w = (max-min) / (double) nshades; /* w = bin size for assigning values->colors*/
   boxsize = ESL_MAX(1, (ESL_MIN((792 - bottommargin) / D->n, 
				 (612 - leftmargin)   / D->m)));
   fboxsize= ESL_MIN( (792. - ((float) bottommargin + topmargin))   / (float) D->n, 
		      (612. - ((float) leftmargin   + rightmargin)) / (float) D->m);


   fprintf(fp, "%.4f %.4f scale\n", (fboxsize/(float) boxsize), (fboxsize/(float) boxsize));
   /* printf("n: %d\nm: %d\n", D->n, D->m); */
   for (i = 0; i < D->n; i++) {
     /* printf("\n"); */
     /* for (j = i; j < D->n; j++) */
     for (j = 0; j < D->m; j++)
       {
	 /* printf("i: %4d j: %4d %5.1f\n", i, j, D->mx[i][j]); */
	 xcoord = j * boxsize + leftmargin;
	 ycoord = (D->m-(i+1)) * boxsize + bottommargin; /* difference w/heatmap.c: (D->m-i+1) */
	 
	 if      (D->mx[i][j]  <  min2fill)    continue;
	 /* if      ((i > 0) && (j > 0) && (D->mx[i][j] <=  D->mx[i-1][j-1]) && (D->mx[i][j] < min2fill))    continue;*/
	 else if (D->mx[i][j] == -eslINFINITY) bin = 0;
	 else if (D->mx[i][j] ==  eslINFINITY) bin = nshades-1;
	 else {
	   bin    = (int) ceil((D->mx[i][j] - min) / w) - 1;
	   if (bin < 0)        bin = 0;
	   if (bin >= nshades) bin = nshades-1;
	 }
	 
	 printf("%4d %4d %10.3f\n", i, j, D->mx[i][j]);
	 fprintf(fp, "newpath\n");
	 fprintf(fp, "  %d %d moveto\n", xcoord, ycoord);
	 fprintf(fp, "  0  %d rlineto\n", boxsize);
	 fprintf(fp, "  %d 0  rlineto\n", boxsize);
	 fprintf(fp, "  0 -%d rlineto\n", boxsize);
	 fprintf(fp, "  closepath\n");
	 fprintf(fp, " %.2f %.2f %.2f %.2f setcmykcolor\n",
		 cyan[bin], magenta[bin], yellow[bin], black[bin]);
	 fprintf(fp, "  fill\n");
       }
   }
  fprintf(fp, "showpage\n");
  return eslOK;
}


/* Function: my_p7_GTraceMSV()
 * Incept:   EPN, Mon Aug 11 10:01:23 2008
 * 
 * Purpose:  Traceback of a MSV matrix: retrieval of optimum alignment.
 *           
 *           Based on p7_GTrace().
 *
 *           This function is currently implemented as a
 *           reconstruction traceback, rather than using a shadow
 *           matrix. Because H3 uses floating point scores, and we
 *           can't compare floats for equality, we have to compare
 *           floats for near-equality and therefore, formally, we can
 *           only guarantee a near-optimal traceback. However, even in
 *           the unlikely event that a suboptimal is returned, the
 *           score difference from true optimal will be negligible.
 *           
 * Args:     dsq    - digital sequence aligned to, 1..L 
 *           L      - length of <dsq>
 *           gm     - profile
 *           mx     - MSV matrix to trace, L x M
 *           tr     - storage for the recovered traceback.
 *           
 * Return:   <eslOK> on success.
 *           <eslFAIL> if even the optimal path has zero probability;
 *           in this case, the trace is set blank (<tr->N = 0>).
 *           <eslEINCOMPAT> if the optimal trace is discontiguous wrt
 *           the sequence, node k > kp emits residue i < ip, while kp
 *           emits ip. In this case, trace is invalid - caller must 
 *           know this.
 *
 * Note:     Care is taken to evaluate the prev+tsc+emission
 *           calculations in exactly the same order that Viterbi did
 *           them, lest you get numerical problems with
 *           a+b+c = d; d-c != a+b because d,c are nearly equal.
 *           (This bug appeared in dev: xref J1/121.)
 */

/* for HMMER3 P7 HMMs */
#define P7MMX(i,k) (dp[(i)][(k) * p7G_NSCELLS + p7G_M])
#define P7XMX(i,s) (xmx[(i) * p7G_NXCELLS + (s)])
#define P7MSC(k)   (rsc[(k) * p7P_NR     + p7P_MSC])

/* for CP9 HMMs */
#define TSC(s,k) (tsc[(k) * cp9O_NTRANS + (s)])

int
my_p7_GTraceMSV(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_GMX *gx, P7_TRACE *tr, int **ret_i2k, int **ret_k2i, float **ret_isc)
{
  int     status = eslOK;
  int     i;			/* position in seq (1..L) */
  int     k;			/* position in model (1..M) */
  int     M   = gm->M;
  float **dp  = gx->dp;
  float  *xmx = gx->xmx;
  float   tol = 1e-5;
  float   esc = p7_profile_IsLocal(gm) ? 0 : -eslINFINITY;
  /* new vars */
  float        tloop = logf((float) L / (float) (L+3));
  float        tmove = logf(     3.0f / (float) (L+3));
  float        tbmk  = logf(     2.0f / ((float) gm->M * (float) (gm->M+1)));
  float        tec   = logf(0.5f);

  int *k2i; /* [0.1..k..gm->M] = i, residue i emitted from node k's match state in MSV trace */
  int *i2k; /* [0.1..i..L]     = k, residue i emitted from node k's match state in MSV trace */
  float *isc; /*[0.1..i..L]    = sc, match emission score for residue i is sc */

  ESL_ALLOC(k2i, sizeof(int)   * (gm->M+1));
  ESL_ALLOC(i2k, sizeof(int)   * (L+1));
  ESL_ALLOC(isc, sizeof(float) * (L+1));
  esl_vec_ISet(k2i, (gm->M+1), -1);
  esl_vec_ISet(i2k, (L+1),     -1);
  esl_vec_FSet(isc, (L+1),     -eslINFINITY);

  if ((status = p7_trace_Reuse(tr)) != eslOK) goto ERROR;

  /* Initialization.
   * (back to front. ReverseTrace() called later.)
   */
  if ((status = p7_trace_Append(tr, p7T_T, 0, 0)) != eslOK) goto ERROR;
  if ((status = p7_trace_Append(tr, p7T_C, 0, 0)) != eslOK) goto ERROR;
  i    = L;			/* next position to explain in seq */

  /* Traceback
   */
  while (tr->st[tr->N-1] != p7T_S) {
    float const *rsc = gm->rsc[dsq[i]];

    switch (tr->st[tr->N-1]) {
    case p7T_C:		/* C(i) comes from E(i) */
      if   (P7XMX(i,p7G_C) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible C reached at i=%d", i);

      if (esl_FCompare(P7XMX(i, p7G_C), P7XMX(i-1, p7G_C) + tloop, tol) == eslOK) {
	tr->i[tr->N-1]    = i--;  /* first C doesn't emit: subsequent ones do */
	status = p7_trace_Append(tr, p7T_C, 0, 0);
      } else if (esl_FCompare(P7XMX(i, p7G_C), P7XMX(i, p7G_E) + tec, tol) == eslOK) 
	status = p7_trace_Append(tr, p7T_E, 0, 0);
      else ESL_XEXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);
      break;

    case p7T_E:		/* E connects from any M state. k set here */
      if (P7XMX(i, p7G_E) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible E reached at i=%d", i);

      if (esl_FCompare(P7XMX(i, p7G_E), P7MMX(i,M), tol) == eslOK) { k = M; status = p7_trace_Append(tr, p7T_M, M, i); }
      else {
	for (k = M-1; k >= 1; k--)
	  if (esl_FCompare(P7XMX(i, p7G_E), P7MMX(i,k) + esc, tol) == eslOK)
	    { status = p7_trace_Append(tr, p7T_M, k, i); break; }
	if (k < 0) ESL_XEXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
      }
      break;

    case p7T_M:			/* M connects from i-1,k-1, or B */
      if (P7MMX(i,k) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible M reached at k=%d,i=%d", k,i);
      if      (esl_FCompare(P7MMX(i,k), P7XMX(i-1,p7G_B) + tbmk  + P7MSC(k), tol) == eslOK) status = p7_trace_Append(tr, p7T_B, 0,   0);
      else if (esl_FCompare(P7MMX(i,k), P7MMX(i-1,k-1)           + P7MSC(k), tol) == eslOK) { 
	status = p7_trace_Append(tr, p7T_M, k-1, i-1);
	if(k2i[(k-1)] != -1) { status = eslEINCOMPAT; printf("! discontiguous trace k2i[k-1=%d] != -1 (%d) i-1 = %d\n", k-1, k2i[(k-1)], i-1); goto ERROR;} 
	if(i2k[(i-1)] != -1) { status = eslEINCOMPAT; printf("! discontiguous trace i2k[i-1=%d] != -1 (%d) k-1 = %d\n", i-1, i2k[(i-1)], k-1); goto ERROR;} 
	k2i[(k-1)] = i-1;
	i2k[(i-1)] = k-1;
	isc[(i-1)] = P7MSC(k);
      }
      else ESL_XEXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);

      if (status != eslOK) goto ERROR;
      k--; 
      i--;
      break;

    case p7T_N:			/* N connects from S, N */
      if (P7XMX(i, p7G_N) <= p7_IMPOSSIBLE) ESL_XEXCEPTION(eslFAIL, "impossible N reached at i=%d", i);

      if (i == 0) status = p7_trace_Append(tr, p7T_S, 0, 0);
      else if (esl_FCompare(P7XMX(i,p7G_N), P7XMX(i-1, p7G_N) + tloop, tol) == eslOK)
	{
	  tr->i[tr->N-1] = i--;
	  status = p7_trace_Append(tr, p7T_N, 0, 0);
	} 
      else ESL_XEXCEPTION(eslFAIL, "N at i=%d couldn't be traced", i);
      break;

    case p7T_B:			/* B connects from N, J */
      if (P7XMX(i,p7G_B) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible B reached at i=%d", i);

      if (esl_FCompare(P7XMX(i,p7G_B), P7XMX(i, p7G_N) + tmove, tol)  == eslOK)
	status = p7_trace_Append(tr, p7T_N, 0, 0);
      else if (esl_FCompare(P7XMX(i,p7G_B),  P7XMX(i, p7G_J) + tmove, tol) == eslOK)
	status = p7_trace_Append(tr, p7T_J, 0, 0);
      else  ESL_XEXCEPTION(eslFAIL, "B at i=%d couldn't be traced", i);
      break;

    case p7T_J:			/* J connects from E(i) or J(i-1) */
      if (P7XMX(i,p7G_J) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible J reached at i=%d", i);

      if (esl_FCompare(P7XMX(i,p7G_J), P7XMX(i-1,p7G_J) + tloop, tol) == eslOK) {
	tr->i[tr->N-1] = i--;
	status = p7_trace_Append(tr, p7T_J, 0, 0);
      } else if (esl_FCompare(P7XMX(i,p7G_J), P7XMX(i-1,p7G_E) + tec, tol) == eslOK) /* note: P7XMX(i-1,p7G_E) differs from Viterbi traceback, where it's P7XMX(i,p7G_E), not sure why */
	status = p7_trace_Append(tr, p7T_E, 0, 0);
      else  ESL_XEXCEPTION(eslFAIL, "J at i=%d couldn't be traced", i);
      break;

    default: ESL_XEXCEPTION(eslFAIL, "bogus state in traceback");
    } /* end switch over statetype[tpos-1] */

    if (status != eslOK) goto ERROR;
  } /* end traceback, at S state */

  if ((status = p7_trace_Reverse(tr)) != eslOK) goto ERROR;
  if (ret_i2k != NULL) { *ret_i2k = i2k; } else free(i2k);
  if (ret_k2i != NULL) { *ret_k2i = k2i; } else free(k2i);
  if (ret_isc != NULL) { *ret_isc = isc; } else free(isc);
  return eslOK;

 ERROR:
  if (ret_i2k != NULL) { *ret_i2k = i2k; } else free(i2k);
  if (ret_k2i != NULL) { *ret_k2i = k2i; } else free(k2i);
  if (ret_isc != NULL) { *ret_isc = isc; } else free(isc);
  return status;
}


    
/* Function: Parsetree2i_to_k()
 * Date:     EPN, Mon Aug 11 13:41:38 2008
 *
 * Purpose:  Given a parsetree, fill a vector i2k[0.1..i..L] = k, saying that
 *           residue i is emitted into consensus column k. If k <= 0, this implies
 *           residue i was inserted after consensus column (-1 * k).
 *
 * Returns:  eslOK on success
 *           eslEINCOMPAT on contract violation.
 */
int 
Parsetree2i_to_k(CM_t *cm, CMEmitMap_t *emap, int L, char *errbuf, Parsetree_t *tr, int **ret_i2k)
{
  int status;                   /* Easel status code */
  int tidx;			/* counter through positions in the parsetree        */
  int v;			/* state index in CM */
  int nd;
  int *i2k;

  ESL_ALLOC(i2k, sizeof(int) * (L+1));
  esl_vec_ISet(i2k, (L+1), -1 * (cm->clen + 1));
  
  /* contract check */
  if(emap   == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "Parsetree2i_to_k(): emap == NULL.");

		/* trivial preorder traverse, since we're already numbered that way */
  for (tidx = 0; tidx < tr->n; tidx++) {
    v = tr->state[tidx];        	/* index of parent state in CM */
    nd = cm->ndidx[v];
    if (v == cm->M) continue;      	/* special case: v is EL, local alignment end */
    switch (cm->sttype[v]) { 
    case MP_st: 
      i2k[tr->emitl[tidx]] = emap->lpos[nd];
      i2k[tr->emitr[tidx]] = emap->rpos[nd];
      break;

    case ML_st:
      i2k[tr->emitl[tidx]] = emap->lpos[nd];
      break;
      
    case MR_st: 
      i2k[tr->emitr[tidx]] = emap->rpos[nd];
      break;

    case IL_st: 
      i2k[tr->emitl[tidx]] = (-1 * emap->lpos[nd]);
      break;

    case IR_st: 
      i2k[tr->emitr[tidx]] = (-1 * (emap->rpos[nd] - 1)); /* IR's emit to before consensus column rpos */
      break;

    }
  }

  *ret_i2k = i2k;
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "Parsetree2i_to_k(), memory allocation error.");
  return status; /* NEVERREACHED */
}

/* Function: prune_i2k()
 * Incept:   EPN, Mon Aug 11 15:02:35 2008
 * 
 * Purpose:  Prune an i2k array of pins. Optionally prune by the following
 *           criteria (in this order):
 *           o match score of the pin 
 *           o length (n) of n-mer each pin exists in
 *           o proximity to end of n-mer 
 *           
 * Args:     i2k        - the input pin array, modified (pruned) in place
 *           isc        - score (match emission) for each pin
 *           L          - length of current sequence
 *           phi        - phi array, phi[k][v] is expected number of times (probability)
 *                        state v (0 = match, 1 insert, 2 = delete) in 
 *                        node k is *entered*. Node 0 is special, state 0 = B state, state 1 = N_state, state 2 = NULL
 *                        Calculated *without* taking insert->insert transitions into account.
 *           min_sc     - minimum score to allow as a pin, 0. to allow any score
 *           min_len    - min n-mer size, 1 to allow any size
 *           min_end    - min distance from end to allow, prune away any others, 0 to not prune based on end proximity
 *           min_mprob  - min match phi probability to allow in a pin 
 *           min_mcprob - min cumulative match phi probability to allow in a nmer pin
 *           max_iprob  - max insert phi probability to allow in a pin 
 *           max_ilprob - max insert phi probability to allow in a state to the left of a pin 
 *           
 * Return:   <eslOK> on success.
 *
 */
int
prune_i2k(int *i2k, float *isc, int L, double **phi, float min_sc, int min_len, int min_end, float min_mprob, float min_mcprob, float max_iprob, float max_ilprob)
{
  int     i, j;

  int do_end;
  int do_sc;
  int do_len;
  float   tol = 1e-5;
  float mcprob = 1.; /* cumulative match probability */

  do_sc   = (esl_FCompare(0., min_sc, tol) == eslOK) ? FALSE : TRUE;
  do_end  = (min_end == 0)  ? FALSE : TRUE;
  do_len  = (min_len == 1) ? FALSE : TRUE;

  if(do_sc) { 
    for(i = 1; i <= L; i++) if((i2k[i] != -1) && (isc[i] < min_sc)) i2k[i] = -1;
  }

  int n = 0;

  /* first pass, prune on phi values */
  for(i = 1; i <= L; i++) { 
    if(i2k[i] != -1) { /* position i is currently pinned */
       if ((phi[i2k[i]][HMMMATCH]      < min_mprob)   || /* match  prob too low */
	   (phi[i2k[i]][HMMINSERT]     > max_iprob)   || /* insert to the right prob too high */
	   (phi[(i2k[i]-1)][HMMINSERT] > max_ilprob)) {  /* insert to the left  prob too high */
	 /* remove pin */
	 /*printf("removing pin %4d %.4f\n", i, phi[i2k[i]][HMMMATCH]);*/
	 i2k[i] = -1;
       }
    }
  }

  /* second pass, prune on nmer length, match score, match phi probability, and distance from end */
  if(do_len || do_end) { /* determine the size n of the n-mer each residue i is in */
    for(i = 1; i <= L; i++) { 
      if((i2k[i] == -1) || /* position i is not pinned */
	 ((i2k[(i-1)] != -1) && !(i2k[(i-1)] == (i2k[i]-1)))) { /* position i is pinned to k, position i-1 is pinned to kp, but k and kp are not consecutive */
	/* we just ended an n-mer (n >= 0) */
	if(n > 0 && n < min_len) { 
	  for(j = (i - n); j < i; j++) i2k[j] = -1; /* remove the n-mer */
	  mcprob = 1.;
	}
	else if (do_end && n >= min_len) { /* n >= min_len, remove those within do_end of end */
	  for(j = (i - n);        j < (i - (n+1)) + min_end; j++) i2k[j] = -1; /* remove the part of the n-mer within nend residues of the beginning edge */
	  for(j = (i - min_end);  j < i;                     j++) i2k[j] = -1; /* remove the part of the n-mer within nend residues of the end edge */
	} 
	if(i2k[i] == -1) { /* position i is not pinned */
	  mcprob = 1.;  
	  n = 0; 
	}
	else { /* position i is pinned to k, position i-1 is pinned to kp, but k and kp are not consecutive */ 
	  mcprob = phi[(i2k[i])][HMMMATCH]; 
	  n = 1;
	  if(mcprob < min_mcprob) { i2k[i] = -1; n = 0; mcprob = 1.; }
	}
      }
      else { 
	n++; 
	mcprob *= phi[(i2k[i])][HMMMATCH]; 
	if(mcprob < min_mcprob) { /* remove the n-mer up to this point */
	  for(j = (i - n); j < i; j++) i2k[j] = -1; /* remove the n-mer */
	  mcprob = 1.;
	  n = 0;
	}
      }
    }
    /* deal with possibility that last n residues were an n-mer, with n < min_len */
    if(n > 0 && n < min_len) { for(j = (i - n); j < i; j++) i2k[j] = -1; /* remove the n-mer */}
  }

  return eslOK;
}

/* Function: p7_pins2bands()
 * Incept:   EPN, Mon Aug 11 15:41:44 2008
 * 
 * Purpose:  Given an i2k pins array, determine the imin and imax bands.
 *           
 * Args:     i2k      - the input pin array, modified (pruned) in place
 *           errbuf   - for error messages
 *           L        - length of current sequence
 *           M        - number of nodes in the HMM
 *           pad      - pad on each side of pin, if pad = 3, we allow +/- 3 residues from pin
 *           ret_kmin - [0.i..L] = k, min node k for residue i
 *           ret_kmax - [0.i..L] = k, max node k for residue i
 *           ret_ncells - number of cells within bands, to return
#if 0
 *           ret_imin - [0.k..M] = i, min residue i for node k
 *           ret_imax - [0.k..M] = i, max residue i for node k
#endif
 *
 * Return:   <eslOK> on success.
 *
 */
int
p7_pins2bands(int *i2k, char *errbuf, int L, int M, int pad, int **ret_kmin, int **ret_kmax, int *ret_ncells)
{
  int     status;

  int i;
  int kn = 1;
  int kx = M;
  int *kmin, *kmax;

  ESL_ALLOC(kmin, sizeof(int) * (L+1));
  ESL_ALLOC(kmax, sizeof(int) * (L+1));

  /* traverse residues left to right to get kmins */
  for(i = 0; i <= L; i++) { 
    if(i2k[i] != -1) { 
      if(kn >= i2k[i] && kn > 1) ESL_FAIL(eslFAIL, errbuf, "p7_pins2bands() error i: %d, i2k[i]: %d but current kn: %d\n", i, i2k[i], kn); 
      kn = ESL_MAX(1, i2k[i] - pad);
    }
    kmin[i] = kn;
  }

  /* traverse nodes right to left to get imaxs */
  for(i = L; i >= 0; i--) { 
    if(i2k[i] != -1) { 
      if(kx <= i2k[i] && kx < L) ESL_FAIL(eslFAIL, errbuf, "p7_pins2bands() error: i: %d, i2k[i]: %d but current kx: %d\n", i, i2k[i], kx); 
      kx = ESL_MIN(M, i2k[i] + pad);
    }
    kmax[i] = kx;
  }

/* M_0 == B state, which must start the parse with i == 0 */
  kmin[0] = 0; 
  kmax[0] = 0; 

  /* get number of cells if wanted */
  int ncells;
  if(ret_ncells != NULL) { 
    ncells = 0;
    for(i = 1; i <= L; i++) ncells += kmax[i] - kmin[i] + 1;
    *ret_ncells = ncells;
  }

  if(ret_kmin != NULL) { *ret_kmin = kmin; } else free(kmin);
  if(ret_kmax != NULL) { *ret_kmax = kmax; } else free(kmax);
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "p7_pins2bands() memory error.");
  return status; /* NEVERREACHED */

#if 0
  /* if we want to get imin/imax instead of kmin/kmax */
  int in = 1;
  int ix = L;
  int *imin;
  int *imax;
  int *k2i;

  ESL_ALLOC(imin, sizeof(int) * (M+1));
  ESL_ALLOC(imax, sizeof(int) * (M+1));
  ESL_ALLOC(k2i,  sizeof(int) * (M+1));

  imin[0] = imax[0] = -1;
  esl_vec_ISet(k2i, (M+1), -1);

  for(i = 1; i <= L; i++) if(i2k[i] != -1) k2i[i2k[i]] = i;

  /* traverse nodes left to right to get imins */
  for(k = 1; k <= M; k++) { 
    if(k2i[k] != -1) { 
      if(k2i[k] != 1 && in >= k2i[k]) ESL_FAIL(eslFAIL, errbuf, "p7_pins2bands() error k: %d, k2i[k]: %d but current in: %d\n", k, k2i[k], in); 
      in = ESL_MAX(1, k2i[k] - pad);
    }
    imin[k] = in;
  }

  /* traverse nodes right to left to get imaxs */
  for(k = M; k >= 1; k--) { 
    if(k2i[k] != L && k2i[k] != -1) { 
      if(ix <= k2i[k]) ESL_FAIL(eslFAIL, errbuf, "p7_pins2bands() error: k: %d, k2i[k]: %d but current ix: %d\n", k, k2i[k], ix); 
      ix = ESL_MIN(L, k2i[k] + pad);
    }
    imax[k] = ix;
  }

  free(k2i);

  /* get number of cells if wanted */
  int ncells;
  if(ret_ncells != NULL) { 
    ncells = 0;
    for(k = 1; k <= M; k++) ncells += imax[k] - imin[k] + 1;
    *ret_ncells = ncells;
  }
  if(ret_imin != NULL) { *ret_imin = imin; } else free(imin);
  if(ret_imax != NULL) { *ret_imax = imax; } else free(imax);
#endif
}

/* Function: DumpP7Bands()
 * Incept:   EPN, Thu Aug 14 08:45:54 2008
 * 
 * Purpose:  Given i2k and kmin, kmax arrays, print them.
 *           
 * Args:     i2k      - [0.k..M] = i, node k is pinned to residue i
 *           kmin     - [0.i..L] = k, min node k for residue i
 *           kmax     - [0.i..L] = k, max node k for residue i
 *           L        - length of current sequence
 *
 * Return:   <eslOK> on success.
 *
 */
int
DumpP7Bands(FILE *fp, int *i2k, int *kmin, int *kmax, int L)
{
  int i;

  fprintf(fp, "# %4s  %4s  %4s ... %4s\n", "i", "i2k", "kmin", "kmax");
  fprintf(fp, "# %4s  %4s  %13s\n", "----", "----", "-------------");
  for(i = 0; i <= L; i++) { 
    fprintf(fp, "  %4d  %4d  %4d ... %4d\n", i, i2k[i], kmin[i], kmax[i]);
  }
  return eslOK;

}

/* Function: cp9_ForwardP7B()
 * 
 * Purpose:  Runs the banded Forward dynamic programming algorithm on an
 *           input subsequence (i0-j0). Complements cp9_BackwardP7B().  
 *           The 'P7B' suffix indicates plan 7 HMM derived bands
 *           in the kmin and kmax arrays are applied.  This function
 *           was derived from cp9_Forward(), differences from that
 *           function were introduced solely to impose bands on the 
 *           matrix. 
 *         
 *           Because of the bands, some options that exist to cp9_Forward()
 *           (like be_efficient and do_scan and reporting hits in scan mode) 
 *           are not available here. This function is meant to be used
 *           solely for first stage of a Forward, Backward, Posterior type
 *           calculation. 
 *
 *           Also due to bands, only L is passed as seq length, instead
 *           of i0 and j0 (seq start/stop).  This simplifies the application
 *           of bands kmin[i]/kmax[i] refers to residue i.
 *
 *           See additional notes in cp9_Forward() "Purpose" section.
 *
 * Args:
 *           cm        - the covariance model, includes cm->cp9: a CP9 HMM
 *           errbuf    - char buffer for error messages
 *           mx        - the matrix, expanded to correct size (if nec), and filled here
 *           dsq       - sequence in digitized form
 *           L         - start of target subsequence (1 for beginning of dsq)
 *           kmin      - [0.1..i..N] minimum k for residue i
 *           kmax      - [0.1..i..N] maximum k for residue i
 *           ret_sc    - RETURN: log P(S|M,bands)/P(S|R), as a bit score
 *
 * Returns:  eslOK on success; 
 *           eslEINCOMPAT on contract violation;
 */
#define INBAND(i,k) ((k >= kmin[i]) && (k <= kmax[i]))

int
cp9_ForwardP7B(CM_t *cm, char *errbuf, CP9_MX *mx, ESL_DSQ *dsq, int L, int *kmin, int *kmax, float *ret_sc)
{
  int          status;
  int          i;           /* j-W: position in the subsequence                             */
  int          cur, prv;    /* rows in DP matrix 0 or 1                                     */
  int          k;           /* CP9 HMM node position                                        */
  int        **mmx;         /* DP matrix for match  state scores [0..1][0..cm->cp9->M]      */
  int        **imx;         /* DP matrix for insert state scores [0..1][0..cm->cp9->M]      */
  int        **dmx;         /* DP matrix for delete state scores [0..1][0..cm->cp9->M]      */
  int        **elmx;        /* DP matrix for EL state scores [0..1][0..cm->cp9->M]          */
  int         *erow;        /* end score for each position [0..1]                           */
  int          c;           /* counter for EL states                                        */
  int          M;           /* cm->cp9->M, query length, number of consensus nodes of model */

  int          kp, kn, kx, kpcur, kpprv, kpprv_el;

  /* Contract checks */
  if(cm->cp9 == NULL)                  ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_Forward, cm->cp9 is NULL.\n");
  if(dsq == NULL)                      ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_Forward, dsq is NULL.");
  if(mx == NULL)                       ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_Forward, mx is NULL.\n");
  if(mx->M != cm->clen)                ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_Forward, mx->M != cm->clen.\n");
  if(cm->clen != cm->cp9->M)           ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_Forward, cm->clen != cm->cp9->M.\n");
  if(kmin == NULL || kmax == NULL)     ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_Forward, kmin and/or kmax == NULL.\n");
    
  M = cm->cp9->M;

  int const *tsc = cm->cp9->otsc; /* ptr to efficiently ordered transition scores           */

  /* gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits. */

  /* Grow DP matrix if nec, to either 2 rows or L+1 rows (depending on be_efficient), 
   * stays M+1 columns */
  if((status = GrowCP9Matrix(mx, errbuf, L, M, kmin, kmax, &mmx, &imx, &dmx, &elmx, &erow)) != eslOK) return status;
  ESL_DPRINTF2(("cp9_ForwardP7B(): CP9 matrix size: %.8f Mb rows: %d.\n", mx->size_Mb, mx->rows));
  
  /* Initialization of the zero row. */
  mmx[0][0] = 0;      /* M_0 is state B, and everything starts in B */
  imx[0][0] = -INFTY; /* I_0 is state N, can't get here without emitting*/
  dmx[0][0] = -INFTY; /* D_0 doesn't exist. */
  elmx[0][0]= -INFTY; /* can't go from B to EL state */
  erow[0]   = -INFTY;   
  
  /* Because there's a D state for every node 1..M, 
     dmx[0][k] is possible for all k 1..M */
  int sc;
  i = 0;
  kn = ESL_MAX(1, kmin[0]);
  kp = kmin[0] - kn;
  for (k = kn; k <= kmax[0]; k++, kp++) { 
    assert(kp >= 0);
    mmx[0][kp] = imx[0][kp] = elmx[0][kp] = -INFTY;      /* need seq to get here */
    sc = -INFTY;
    if(kp > 0) { /* if kp == 0, kp-1 is outside the bands */
      sc = ILogsum(ILogsum(mmx[0][kp-1] + TSC(cp9O_MD,k-1),
			   imx[0][kp-1] + TSC(cp9O_ID,k-1)),
		   dmx[0][kp-1] + TSC(cp9O_DD,k-1));
    }
    dmx[0][k] = sc;
  }
  /* We can do a full parse through all delete states. */
  erow[0] = -INFTY;
  if(INBAND(0, M)) { erow[0] = dmx[0][M] + TSC(cp9O_DM,M); }
     
  /*****************************************************************
   * The main loop: scan the sequence from position 1 to L.
   *****************************************************************/
  /* Recursion. */
  
  for (i = 1; i <= L; i++) 
    { 
      prv = i-1;
      cur = i;
      int const *isc = cm->cp9->isc[dsq[i]];
      int const *msc = cm->cp9->msc[dsq[i]];
      int endsc     = -INFTY;
      int el_selfsc = cm->cp9->el_selfsc;
      int sc;

      if(kmin[i] == 0) { 
	mmx[i][0]  = -INFTY;
	dmx[i][0]  = -INFTY;  /*D_0 is non-existent*/
	elmx[i][0] = -INFTY;  /*no EL state for node 0 */
	sc = ILogsum(ILogsum(mmx[i-1][0] + TSC(cp9O_MI,0),
			     imx[i-1][0] + TSC(cp9O_II,0)),
		     dmx[i-1][0] + TSC(cp9O_DI,0));
	imx[i][0] = sc + isc[0];
	kn = 1; /* kmin[i] + 1 */
      }
      else { 
	kn = kmin[i];
      }

      /*match state*/
      kn = ESL_MAX(kn, (kmin[i-1]+1)); /* start at first cell from which we can look back to a valid cell at *mx[i-1][k-1] */
      kx = ESL_MIN(kmax[i], kmax[i-1]+1);

      /* NOT SURE ABOUT THIS AND HOW IT COUPLES WITH BLOCK ABOVE kmin[i] == 0 */ 
      for (kpcur = 0;            kpcur < (kn - kmin[i]);   kpcur++) mmx[i][kpcur] = -INFTY; /* impossible to reach these guys */
      for (kpcur = kx-kmin[i]+1; kpcur <= kmax[i]-kmin[i]; kpcur++) mmx[i][kpcur] = -INFTY; /* impossible to reach these guys */

      for (kpcur = 0;            kpcur < (kn - kmin[i]);   kpcur++) elmx[i][kpcur] = -INFTY; /* impossible to reach these guys */
      for (kpcur = kx-kmin[i]+1; kpcur <= kmax[i]-kmin[i]; kpcur++) elmx[i][kpcur] = -INFTY; /* impossible to reach these guys */

      kpcur = kn - kmin[i]; /* unnec, loop above ends with this */
      kpprv = kn - kmin[i-1];
      for (k = kn; k <= kx; k++, kpcur++, kpprv++) {
	/*printf("M i: %d kpprv: %d\n", i, kpprv);*/
	assert((kpprv-1) >= 0);

	sc = ILogsum(ILogsum(mmx[i-1][kpprv-1] + TSC(cp9O_MM,k-1),
			     imx[i-1][kpprv-1] + TSC(cp9O_IM,k-1)),
		     dmx[i-1][kpprv-1] + TSC(cp9O_DM,k-1));

	/* FIX ME: inefficient! check B->M_K transition */
	if(INBAND(i-1, 0)) { /* if i-1 is in k == 0's band */
	  assert(kmin[(i-1)] == 0);
	  if(mmx[i-1][0] != -INFTY) 
	   sc = ILogsum(sc, mmx[i-1][0] + TSC(cp9O_BM,k));
	}

	/* check possibility we came from an EL, if they're valid */
	for(c = 0; c < cm->cp9->el_from_ct[k]; c++) { /* el_from_ct[k] is >= 0 */
	  if(INBAND(i-1, cm->cp9->el_from_idx[k][c])) { 
	    kpprv_el = cm->cp9->el_from_idx[k][c] - kmin[(i-1)];
	    sc = ILogsum(sc, elmx[i-1][kpprv_el]);
	  }
	} /* transition penalty to EL incurred when EL was entered */
	if(sc != -INFTY) { 
	  mmx[i][kpcur] = sc + msc[k];
	  /* E state update */
	  endsc = ILogsum(endsc, mmx[i][kpcur] + TSC(cp9O_ME,k));
	}
	else { 
	  mmx[i][kpcur] = -INFTY;
	  /* don't update E state */
	}
	/*printf("k: %4d mmx[i:%4d][kpcur:%4d]: %d\n", k, i, kpcur, mmx[i][kpcur]);*/

	/* el state */
	sc = -INFTY;
	if((cm->cp9->flags & CPLAN9_EL) && cm->cp9->has_el[k]) /* not all HMM nodes have an EL state (for ex: 
								  HMM nodes that map to right half of a MATP_MP) */
	  {
	    if(mmx[i][kpcur] != -INFTY) { 
	      kpprv_el = k - kmin[(i-1)];
	      sc = ILogsum(mmx[i][kpcur]  + TSC(cp9O_MEL,k), /* transitioned from cur node's match state */
			   elmx[i-1][kpprv_el] + el_selfsc);      /* transitioned from cur node's EL state emitted ip on transition */
	    }
	  }
	elmx[i][kpcur] = sc;
      }

      /* insert state*/
      kn = ESL_MAX(kmin[i], kmin[i-1]);
      kx = ESL_MIN(kmax[i], kmax[i-1]);

      for (kpcur = 0;            kpcur < (kn - kmin[i]);   kpcur++) imx[i][kpcur] = -INFTY; /* impossible to reach these guys */
      for (kpcur = kx-kmin[i]+1; kpcur <= kmax[i]-kmin[i]; kpcur++) imx[i][kpcur] = -INFTY; /* impossible to reach these guys */

      kpcur = kn - kmin[i]; /* unnec, loop above ends with this */
      kpprv = kn - kmin[i-1];
      for (k = kn; k <= kx; k++, kpcur++, kpprv++) { 
	/*insert state*/
	assert(kpprv >= 0);
	/* HERE, EVENTUALLY IF kmin/kmax differ b/t Match and Inserts: 
	 * only look at match states from k that have i-1 within band */
	/* all insert states from k should have i-1 within band */
	/*printf("I i: %d kpprv: %d\n", i, kpprv);*/
	sc = ILogsum(ILogsum(mmx[i-1][kpprv] + TSC(cp9O_MI,k),
			     imx[i-1][kpprv] + TSC(cp9O_II,k)),
		     dmx[i-1][kpprv] + TSC(cp9O_DI,k));
	if(sc != -INFTY) imx[i][kpcur] = sc + isc[k];
	else             imx[i][kpcur] = -INFTY;
	/*printf("k: %4d imx[i:%4d][kpcur:%4d]: %d\n", k, i, kpcur, imx[i][kpcur]);*/
      }

      /*delete state*/
      kn = kmin[i]+1;

      for (kpcur = 0; kpcur < (kn - kmin[i]); kpcur++) dmx[i][kpcur] = -INFTY; /* impossible to reach these guys */
      kpcur = kn - kmin[i]; /* unnec, loop above ends with this */
      for (k = kn; k <= kmax[i]; k++, kpcur++) { /* should I be adding one for delete off-by-one?? */
	sc = ILogsum(ILogsum(mmx[i][kpcur-1] + TSC(cp9O_MD,k-1),
			     imx[i][kpcur-1] + TSC(cp9O_ID,k-1)),
		     dmx[i][kpcur-1] + TSC(cp9O_DD,k-1));
	dmx[i][kpcur] = sc;
	/*printf("k: %4d dmx[i:%4d][kpcur:%4d]: %d\n", k, i, kpcur, dmx[i][kpcur]);*/
      }
	  /*printf("mmx [jp:%d][%d]: %d\n", jp, k, mmx[j][k]);
	    printf("imx [jp:%d][%d]: %d\n", jp, k, imx[j][k]);
	    printf("dmx [jp:%d][%d]: %d\n", jp, k, dmx[j][k]);
	    printf("elmx[jp:%d][%d]: %d\n", jp, k, elmx[j][k]);*/

      if(INBAND(i, M)) { 
	endsc = ILogsum(ILogsum(endsc, dmx[i][M-kmin[i]] + TSC(cp9O_DM,M)), /* transition from D_M -> end */
			imx[i][M-kmin[i]] + TSC(cp9O_IM,M)); /* transition from I_M -> end */
	for(c = 0; c < cm->cp9->el_from_ct[M+1]; c++) { /* el_from_ct[k] is >= 0 */
	  if(INBAND(i, cm->cp9->el_from_idx[M+1][c])) { 
	    kpprv_el = cm->cp9->el_from_idx[M+1][c] - kmin[i];
	    endsc = ILogsum(endsc, elmx[i][kpprv_el]);
	  }
	}
      }
	/* transition penalty to EL incurred when EL was entered */
      /*printf("endsc: %d\n", endsc);*/

      erow[i] = endsc;
    } /* end loop over end positions i */
  
  *ret_sc = Scorify(erow[L]);
  printf("cp9_ForwardP7B() return score: %10.4f\n", Scorify(erow[L]));
  ESL_DPRINTF1(("cp9_ForwardP7B() return score: %10.4f\n", Scorify(erow[L])));

  return eslOK;
}

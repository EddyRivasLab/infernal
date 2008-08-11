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
 *
 * Note:     Care is taken to evaluate the prev+tsc+emission
 *           calculations in exactly the same order that Viterbi did
 *           them, lest you get numerical problems with
 *           a+b+c = d; d-c != a+b because d,c are nearly equal.
 *           (This bug appeared in dev: xref J1/121.)
 */

#define MMX(i,k) (dp[(i)][(k) * p7G_NSCELLS + p7G_M])
#define IMX(i,k) (dp[(i)][(k) * p7G_NSCELLS + p7G_I])
#define DMX(i,k) (dp[(i)][(k) * p7G_NSCELLS + p7G_D])
#define XMX(i,s) (xmx[(i) * p7G_NXCELLS + (s)])

#define TSC(s,k) (tsc[(k) * p7P_NTRANS + (s)])
#define MSC(k)   (rsc[(k) * p7P_NR     + p7P_MSC])
#define ISC(k)   (rsc[(k) * p7P_NR     + p7P_ISC])

int
my_p7_GTraceMSV(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_GMX *gx, P7_TRACE *tr, int **ret_i2k, int **ret_k2i)
{
  int     status;
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

  ESL_ALLOC(k2i, sizeof(int) * (gm->M+1));
  ESL_ALLOC(i2k, sizeof(int) * (L+1));
  esl_vec_ISet(k2i, (gm->M+1), -1);
  esl_vec_ISet(i2k, (L+1),     -1);

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
      if   (XMX(i,p7G_C) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible C reached at i=%d", i);

      if (esl_FCompare(XMX(i, p7G_C), XMX(i-1, p7G_C) + tloop, tol) == eslOK) {
	tr->i[tr->N-1]    = i--;  /* first C doesn't emit: subsequent ones do */
	status = p7_trace_Append(tr, p7T_C, 0, 0);
      } else if (esl_FCompare(XMX(i, p7G_C), XMX(i, p7G_E) + tec, tol) == eslOK) 
	status = p7_trace_Append(tr, p7T_E, 0, 0);
      else ESL_XEXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);
      break;

    case p7T_E:		/* E connects from any M state. k set here */
      if (XMX(i, p7G_E) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible E reached at i=%d", i);

      if (esl_FCompare(XMX(i, p7G_E), MMX(i,M), tol) == eslOK) { k = M; status = p7_trace_Append(tr, p7T_M, M, i); }
      else {
	for (k = M-1; k >= 1; k--)
	  if (esl_FCompare(XMX(i, p7G_E), MMX(i,k) + esc, tol) == eslOK)
	    { status = p7_trace_Append(tr, p7T_M, k, i); break; }
	if (k < 0) ESL_XEXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
      }
      break;

    case p7T_M:			/* M connects from i-1,k-1, or B */
      if (MMX(i,k) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible M reached at k=%d,i=%d", k,i);
      if      (esl_FCompare(MMX(i,k), XMX(i-1,p7G_B) + tbmk  + MSC(k), tol) == eslOK) status = p7_trace_Append(tr, p7T_B, 0,   0);
      else if (esl_FCompare(MMX(i,k), MMX(i-1,k-1)           + MSC(k), tol) == eslOK) { 
	status = p7_trace_Append(tr, p7T_M, k-1, i-1);
	if(k2i[(k-1)] != -1) ESL_XEXCEPTION(eslFAIL, "k2i[k-1=%d] != -1 (%d) i-1 = %d\n", k-1, k2i[(k-1)], i-1);
	if(i2k[(i-1)] != -1) ESL_XEXCEPTION(eslFAIL, "i2k[i-1=%d] != -1 (%d) k-1 = %d\n", i-1, i2k[(i-1)], k-1);
	k2i[(k-1)] = i-1;
	i2k[(i-1)] = k-1;
      }
      else ESL_XEXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);

      if (status != eslOK) goto ERROR;
      k--; 
      i--;
      break;

    case p7T_N:			/* N connects from S, N */
      if (XMX(i, p7G_N) <= p7_IMPOSSIBLE) ESL_XEXCEPTION(eslFAIL, "impossible N reached at i=%d", i);

      if (i == 0) status = p7_trace_Append(tr, p7T_S, 0, 0);
      else if (esl_FCompare(XMX(i,p7G_N), XMX(i-1, p7G_N) + tloop, tol) == eslOK)
	{
	  tr->i[tr->N-1] = i--;
	  status = p7_trace_Append(tr, p7T_N, 0, 0);
	} 
      else ESL_XEXCEPTION(eslFAIL, "N at i=%d couldn't be traced", i);
      break;

    case p7T_B:			/* B connects from N, J */
      if (XMX(i,p7G_B) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible B reached at i=%d", i);

      if (esl_FCompare(XMX(i,p7G_B), XMX(i, p7G_N) + tmove, tol)  == eslOK)
	status = p7_trace_Append(tr, p7T_N, 0, 0);
      else if (esl_FCompare(XMX(i,p7G_B),  XMX(i, p7G_J) + tmove, tol) == eslOK)
	status = p7_trace_Append(tr, p7T_J, 0, 0);
      else  ESL_XEXCEPTION(eslFAIL, "B at i=%d couldn't be traced", i);
      break;

    case p7T_J:			/* J connects from E(i) or J(i-1) */
      if (XMX(i,p7G_J) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible J reached at i=%d", i);

      if (esl_FCompare(XMX(i,p7G_J), XMX(i-1,p7G_J) + tloop, tol) == eslOK) {
	tr->i[tr->N-1] = i--;
	status = p7_trace_Append(tr, p7T_J, 0, 0);
      } else if (esl_FCompare(XMX(i,p7G_J), XMX(i-1,p7G_E) + tec, tol) == eslOK) /* note: XMX(i-1,p7G_E) differs from Viterbi traceback, where it's XMX(i,p7G_E), not sure why */
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
  return eslOK;

 ERROR:
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

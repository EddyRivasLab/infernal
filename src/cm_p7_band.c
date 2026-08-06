/* cm_p7_band.c
 * 
 * Functions for p7 HMM banding.
 * BEWARE: only partially implemented.
 */

#include <esl_config.h>
#include <p7_config.h>
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdint.h>
#include <assert.h>

#include <xmmintrin.h>
#include <emmintrin.h>
#include <smmintrin.h>

#include "easel.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sse.h"
#ifdef HMMER_THREADS
#include "esl_threads.h"
#include "esl_workqueue.h"
#endif
#include "esl_vectorops.h"

#include "hmmer.h"
#include "p7_gbands.h"
#include "p7_gmxb.h"

#include "infernal.h"

/* Function:  pn_match_bands_enforce_monotone()
 * Incept:    EPN, 2026-04-24
 *
 * Purpose:   Enforce a monotone reachability sweep on the match-state
 *            HMM bands pn_min_m[k] and pn_max_m[k], k=1..M.
 *
 *            After: pn_min_m[k] <= min(pn_min_m[k..M])
 *                   pn_max_m[k] >= max(pn_max_m[1..k])
 *
 *            Rationale: when pn_min/max_m come from p7-banded F/B
 *            posteriors, cells outside the p7 Viterbi diagonal band
 *            have zero mass, so pn_min_m[k] at downstream nodes can
 *            exclude envelope positions that are reachable via valid
 *            HMM paths whose match/delete layout differs from the p7
 *            Viterbi trace. The sweep relaxes only the match bands and
 *            only where a predecessor/successor bound was strictly
 *            tighter. Unset entries (pn_min_m[k]==-1) are skipped.
 *
 *            If <do_print> is TRUE, emits one line to stderr:
 *              pnmono[<ctx>] M=.. L=.. sum_width_before=.. sum_width_after=.. delta=..
 *
 * Returns:   (void). pn_min_m/pn_max_m are updated in place.
 */
void
pn_match_bands_enforce_monotone(int *pn_min_m, int *pn_max_m, int M, int L,
                                int do_print, const char *ctx)
{
  int k;
  long sum_before = 0, sum_after = 0;

  if(do_print) {
    for(k = 1; k <= M; k++) {
      if(pn_min_m[k] != -1) sum_before += (pn_max_m[k] - pn_min_m[k] + 1);
    }
  }

  /* Reachability: M_k must precede M_{k+1}. So:
   *   pn_min_m[k] <= pn_min_m[k'] for all k' > k  (earliest-start non-decreasing in k)
   *   pn_max_m[k] <= pn_max_m[k'] for all k' > k  (latest-start also non-decreasing)
   * If p7-banded posteriors leave a downstream node pn_min_m[k']=p0 but an
   * upstream node pn_min_m[k]=p1 with p1 > p0, that violates reachability.
   * Relax: sweep right-to-left, widen pn_min_m[k] down to min(pn_min_m[k..M]).
   * Symmetrically: sweep left-to-right, widen pn_max_m[k] up to max(pn_max_m[1..k]). */
  {
    int running_min = L + 2;
    for(k = M; k >= 1; k--) {
      if(pn_min_m[k] != -1 && pn_min_m[k] < running_min) running_min = pn_min_m[k];
      if(pn_min_m[k] != -1 && running_min < pn_min_m[k])  pn_min_m[k] = running_min;
    }
  }
  {
    int running_max = -1;
    for(k = 1; k <= M; k++) {
      if(pn_max_m[k] > running_max) running_max = pn_max_m[k];
      if(pn_max_m[k] != -1 && running_max > pn_max_m[k]) pn_max_m[k] = running_max;
    }
  }

  if(do_print) {
    for(k = 1; k <= M; k++) {
      if(pn_min_m[k] != -1) sum_after += (pn_max_m[k] - pn_min_m[k] + 1);
    }
    fprintf(stderr, "pnmono[%s] M=%d L=%d sum_width_before=%ld sum_width_after=%ld delta=%ld\n",
            ctx ? ctx : "?", M, L, sum_before, sum_after, sum_after - sum_before);
  }
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
#define CP9TSC(s,k) (tsc[(k) * cp9O_NTRANS + (s)])

int
my_p7_GTraceMSV(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_GMX *gx, P7_TRACE *tr, int **ret_i2k, int **ret_k2i, float **ret_isc, int **ret_iconflict)
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
  int *iconflict; /* [0.1..i..L] = TRUE | FALSE. TRUE to eventually remove the kmer with pin at i */
  float *isc; /*[0.1..i..L]    = sc, match emission score for residue i is sc */
  
  ESL_ALLOC(k2i, sizeof(int)   * (gm->M+1));
  ESL_ALLOC(i2k, sizeof(int)   * (L+1));
  ESL_ALLOC(iconflict, sizeof(int)   * (L+1));
  ESL_ALLOC(isc, sizeof(float) * (L+1));
  esl_vec_ISet(k2i, (gm->M+1), -1);
  esl_vec_ISet(i2k, (L+1),     -1);
  esl_vec_ISet(iconflict, (L+1), FALSE);
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

      if (esl_FCompare_old(P7XMX(i, p7G_C), P7XMX(i-1, p7G_C) + tloop, tol) == eslOK) {
	tr->i[tr->N-1]    = i--;  /* first C doesn't emit: subsequent ones do */
	status = p7_trace_Append(tr, p7T_C, 0, 0);
      } else if (esl_FCompare_old(P7XMX(i, p7G_C), P7XMX(i, p7G_E) + tec, tol) == eslOK) 
	status = p7_trace_Append(tr, p7T_E, 0, 0);
      else ESL_XEXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);
      break;

    case p7T_E:		/* E connects from any M state. k set here */
      if (P7XMX(i, p7G_E) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible E reached at i=%d", i);

      if (esl_FCompare_old(P7XMX(i, p7G_E), P7MMX(i,M), tol) == eslOK) { k = M; status = p7_trace_Append(tr, p7T_M, M, i); }
      else {
	for (k = M-1; k >= 1; k--)
	  if (esl_FCompare_old(P7XMX(i, p7G_E), P7MMX(i,k) + esc, tol) == eslOK)
	    { status = p7_trace_Append(tr, p7T_M, k, i); break; }
	if (k < 0) ESL_XEXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
      }
      break;

    case p7T_M:			/* M connects from i-1,k-1, or B */
      if (P7MMX(i,k) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible M reached at k=%d,i=%d", k,i);
      if      (esl_FCompare_old(P7MMX(i,k), P7XMX(i-1,p7G_B) + tbmk  + P7MSC(k), tol) == eslOK) status = p7_trace_Append(tr, p7T_B, 0,   0);
      else if (esl_FCompare_old(P7MMX(i,k), P7MMX(i-1,k-1)           + P7MSC(k), tol) == eslOK) { 
	status = p7_trace_Append(tr, p7T_M, k-1, i-1);
	/*if(k2i[(k-1)] != -1) { status = eslEINCOMPAT; printf("! discontiguous trace k2i[k-1=%d] != -1 (%d) i-1 = %d\n", k-1, k2i[(k-1)], i-1); goto ERROR;} */
	/*if(i2k[(i-1)] != -1) { status = eslEINCOMPAT; printf("! discontiguous trace i2k[i-1=%d] != -1 (%d) k-1 = %d\n", i-1, i2k[(i-1)], k-1); goto ERROR;} */
	if(i2k[(i-1)] != -1 || k2i[(k-1)] != -1) {
	  if (getenv("MSVBAND_DBG")) {
	    fprintf(stderr, "#MSVBAND_DBG: conflict at i=%d k=%d: i2k[i-1=%d]=%d k2i[k-1=%d]=%d\n",
	            i, k, i-1, i2k[(i-1)], k-1, k2i[(k-1)]);
	  }
	  iconflict[(k2i[(k-1)])] = TRUE; /* eventually remove pin kmer that included k we *thought* pinned i-1 */
	  iconflict[(i-1)]        = TRUE; /* eventually remove pin kmer that includes k we currently think pins i-1 */
	}
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
      if (P7XMX(i, p7G_N) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible N reached at i=%d", i);

      if (i == 0) status = p7_trace_Append(tr, p7T_S, 0, 0);
      else if (esl_FCompare_old(P7XMX(i,p7G_N), P7XMX(i-1, p7G_N) + tloop, tol) == eslOK)
	{
	  tr->i[tr->N-1] = i--;
	  status = p7_trace_Append(tr, p7T_N, 0, 0);
	} 
      else ESL_XEXCEPTION(eslFAIL, "N at i=%d couldn't be traced", i);
      break;

    case p7T_B:			/* B connects from N, J */
      if (P7XMX(i,p7G_B) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible B reached at i=%d", i);

      if (esl_FCompare_old(P7XMX(i,p7G_B), P7XMX(i, p7G_N) + tmove, tol)  == eslOK)
	status = p7_trace_Append(tr, p7T_N, 0, 0);
      else if (esl_FCompare_old(P7XMX(i,p7G_B),  P7XMX(i, p7G_J) + tmove, tol) == eslOK)
	status = p7_trace_Append(tr, p7T_J, 0, 0);
      else  ESL_XEXCEPTION(eslFAIL, "B at i=%d couldn't be traced", i);
      break;

    case p7T_J:			/* J connects from E(i) or J(i-1) */
      if (P7XMX(i,p7G_J) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible J reached at i=%d", i);

      if (esl_FCompare_old(P7XMX(i,p7G_J), P7XMX(i-1,p7G_J) + tloop, tol) == eslOK) {
	tr->i[tr->N-1] = i--;
	status = p7_trace_Append(tr, p7T_J, 0, 0);
      } else if (esl_FCompare_old(P7XMX(i,p7G_J), P7XMX(i,p7G_E) + tec, tol) == eslOK) 
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
  if (ret_iconflict != NULL) { *ret_iconflict = iconflict; } else free(iconflict);
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
 *           iconflict  - [0.1..i..L] TRUE if i:i2k[i] pin is involved in a conflict, remove the whole nmer pin
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
prune_i2k(int *i2k, int *iconflict, float *isc, int L, double **phi, float min_sc, int min_len, int min_end, float min_mprob, float min_mcprob, float max_iprob, float max_ilprob)
{
  int     i, j;

  int do_end;
  int do_sc;
  int do_len;
  float   tol = 1e-5;
  float mcprob = 1.; /* cumulative match probability */
  int n = 0;
  int ip;

  do_sc   = (esl_FCompare_old(0., min_sc, tol) == eslOK) ? FALSE : TRUE;
  do_end  = (min_end == 0)  ? FALSE : TRUE;
  do_len  = (min_len == 1) ? FALSE : TRUE;

  /* pass 1, prune on conflict from trace, this must be first pass */
  for(i = 1; i <= L; i++) { 
    if(iconflict[i] == TRUE) {
      /* remove entire nmer's pins */
      ip = i;
      while((ip > 0)  && (i2k[ip] != -1)) i2k[ip--] = -1;
      ip = i;
      while((ip <= L) && (i2k[ip] != -1)) i2k[ip++] = -1;
    }
  }

  /* pass 2, prune on scores */
  if(do_sc) { 
    for(i = 1; i <= L; i++) if((i2k[i] != -1) && (isc[i] < min_sc)) i2k[i] = -1;
  }

  /* pass 3, prune on phi values */
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

  /* pass 4, prune on nmer length, match score, match phi probability, and distance from end */
  if(do_len || do_end) { /* determine the size n of the n-mer each residue i is in */
    for(i = 1; i <= L; i++) { 
      if((i2k[i] == -1) || /* position i is not pinned */
	 ((i2k[(i-1)] != -1) && !(i2k[(i-1)] == (i2k[i]-1)))) { /* position i is pinned to k, position i-1 is pinned to kp, but k and kp are not consecutive */
	/* we just ended an n-mer (n >= 0) */
	if(n > 0 && n < min_len) { 
	  for(j = (i - n); j < i; j++) i2k[j] = -1; /* remove the n-mer */
	  mcprob = 1.;
	}
	else if (do_end && n > 0 && n >= min_len) { /* n >= min_len, remove those within do_end of end */
	  for(j = (i - n); j < ESL_MIN(i, (i - (n+1)) + min_end); j++) i2k[j] = -1; /* remove the part of the n-mer within min_end residues of the beginning edge */
	  for(j = ESL_MAX(1, i - min_end); j < i; j++) i2k[j] = -1;                   /* remove the part of the n-mer within min_end residues of the end edge */
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
    /* deal with possibility that last n residues were an n-mer ending at position L */
    if(n > 0 && n < min_len) { for(j = (i - n); j < i; j++) i2k[j] = -1; }
    else if(do_end && n > 0 && n >= min_len) {
      for(j = (i - n); j < ESL_MIN(i, (i - (n+1)) + min_end); j++) i2k[j] = -1;
      for(j = ESL_MAX(1, i - min_end); j < i; j++) i2k[j] = -1;
    }
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
  int kn = 0;
  int kx = M;
  int *kmin, *kmax;

  ESL_ALLOC(kmin, sizeof(int) * (L+1));
  ESL_ALLOC(kmax, sizeof(int) * (L+1));

  /* traverse residues left to right to get kmins */
  for(i = 0; i <= L; i++) {
    if(i2k[i] != -1) {
      if(kn >= i2k[i] && kn > 1) {
	i2k[i] = -1; /* non-monotone pin from multi-segment MSV trace; remove and continue */
      }
      else {
	kn = ESL_MAX(1, i2k[i] - pad);
      }
    }
    kmin[i] = kn;
  }

  /* traverse nodes right to left to get imaxs */
  for(i = L; i >= 0; i--) {
    if(i2k[i] != -1) {
      if(kx <= i2k[i] && kx < M) { i2k[i] = -1; } /* non-monotone pin; remove and continue */
      else                        { kx = ESL_MIN(M, i2k[i] + pad); }
    }
    kmax[i] = kx;
  }

  /* M_0 == B state, which must start the parse with i == 0 */
  kmin[0] = 0;

  /* D-state bridge: ensure the band at each emitting position i
   * extends far enough in k to reach the next emitting position's k
   * via delete states. Without this, deletion runs in the trace
   * (M_k(i) → D_{k+1} → ... → M_{k'}(i+1)) can't propagate through
   * the banded Forward DP because the D-state columns k+1..k'-1
   * aren't in the band at row i.
   *
   * Fix: for each pair of consecutive pins (i with k) and (i' with k'),
   * if k' > kmax[i], widen kmax[i] to k' (and also kmax at any
   * intermediate unpinned rows between i and i'). Symmetrically,
   * if k < kmin[i'], widen kmin[i'] down to k.
   */
  {
    int prev_pin_i = -1, prev_pin_k = -1;
    for(i = 0; i <= L; i++) {
      if(i2k[i] != -1) {
        if(prev_pin_i >= 0) {
          /* Bridge from previous pin to current pin */
          int pk = prev_pin_k;
          int ck = i2k[i];
          if(ck > pk) {
            /* Need D-state columns pk+1..ck at rows prev_pin_i..i-1 */
            int j;
            for(j = prev_pin_i; j < i; j++) {
              if(kmax[j] < ck) kmax[j] = ESL_MIN(M, ck);
            }
          }
          if(pk > ck) {
            /* Reverse: need D-state columns ck..pk-1 at rows prev_pin_i+1..i */
            int j;
            for(j = prev_pin_i + 1; j <= i; j++) {
              if(kmin[j] > ck) kmin[j] = ESL_MAX(1, ck);
            }
          }
        }
        prev_pin_i = i;
        prev_pin_k = i2k[i];
      }
    }
  }

  /* get number of cells if wanted.
   * Brief 26_0628-027: accumulate in int64. At genome scale (M~L~1.5e5, and sparse
   * anchors -> wide bands) the true cell count reaches ~1e10, overflowing the
   * int32 accumulator and wrapping to <=0, which would spuriously trip the
   * caller's ncells==0 fallback. ret_ncells stays int (every caller uses it
   * only as a >0/==0 discriminator + debug avg_bw; it never sizes memory), so
   * we saturate the returned value at INT_MAX and emit the true int64 count on
   * a debug line for genome-scale band-size measurement. */
  if(ret_ncells != NULL) {
    int64_t ncells64 = 0;
    for(i = 1; i <= L; i++) ncells64 += (int64_t)(kmax[i] - kmin[i] + 1);
    if(ncells64 > (int64_t) INT_MAX) {
      fprintf(stderr, "#NCELLS64 p7_pins2bands L=%d M=%d ncells=%lld (saturated to INT_MAX in int32 ret_ncells)\n",
              L, M, (long long) ncells64);
      *ret_ncells = INT_MAX;
    } else {
      *ret_ncells = (int) ncells64;
    }
  }

  if(ret_kmin != NULL) { *ret_kmin = kmin; } else free(kmin);
  if(ret_kmax != NULL) { *ret_kmax = kmax; } else free(kmax);
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "p7_pins2bands() memory error.");
  return status; /* NEVERREACHED */
}

/* Function: p7_pins2bands_nodepad()
 * Date:     EPN*, Sat Apr  5 2026
 *
 * Purpose:  Same as p7_pins2bands() but uses a per-node pad array
 *           <nodepad> instead of a single pad value. nodepad[k] is
 *           the pad for node k.
 *
 * Args:     i2k      - the input pin array, modified (pruned) in place
 *           errbuf   - for error messages
 *           L        - length of current sequence
 *           M        - number of nodes in the HMM
 *           nodepad  - [0..M] per-node pad array
 *           hopback  - if >0, dilate band per residue with min/max of i2k
 *                      across the 2*hopback+1-pin window (D1 strategy).
 *                      Vit pins are monotone in (i,k); the dilation is
 *                      computed in O(L) with a sliding window over the
 *                      pinned positions in trace order. 0 = off (no-op).
 *           alpha    - brief 26_0628-043: distance-scaled slack coefficient for the
 *                      interpolated-ramp inter-pin band (brief 26_0628-042); see
 *                      RAMP_SLACK_ALPHA comment below. Callers that don't
 *                      expose a tunable knob should pass 0.75 (brief 26_0628-042's
 *                      validated default) to preserve prior behavior.
 *           ret_kmin - [0.i..L] = k, min node k for residue i
 *           ret_kmax - [0.i..L] = k, max node k for residue i
 *           ret_ncells - number of cells within bands, to return
 *
 * Return:   <eslOK> on success.
 *
 */
int
p7_pins2bands_nodepad(int *i2k, char *errbuf, int L, int M, int *nodepad,
                      int hopback, double alpha,
                      int **ret_kmin, int **ret_kmax, int *ret_ncells)
{
  int     status;

  int i;
  int kn = 0;
  int kx = M;
  int *kmin, *kmax;

  ESL_ALLOC(kmin, sizeof(int) * (L+1));
  ESL_ALLOC(kmax, sizeof(int) * (L+1));

  /* traverse residues left to right to get kmins */
  for(i = 0; i <= L; i++) {
    if(i2k[i] != -1) {
      if(kn >= i2k[i] && kn > 1) {
	i2k[i] = -1; /* non-monotone pin from multi-segment MSV trace; remove and continue */
      }
      else {
	kn = ESL_MAX(1, i2k[i] - nodepad[i2k[i]]);
      }
    }
    kmin[i] = kn;
  }

  /* traverse nodes right to left to get imaxs */
  for(i = L; i >= 0; i--) {
    if(i2k[i] != -1) {
      if(kx <= i2k[i] && kx < M) { i2k[i] = -1; } /* non-monotone pin; remove and continue */
      else                        { kx = ESL_MIN(M, i2k[i] + nodepad[i2k[i]]); }
    }
    kmax[i] = kx;
  }

  /* M_0 == B state, which must start the parse with i == 0 */
  kmin[0] = 0;

  /* D-state bridge (same as p7_pins2bands): widen bands between consecutive
   * pins to allow delete-state transitions through the banded Forward DP. */
  {
    int prev_pin_i = -1, prev_pin_k = -1;
    for(i = 0; i <= L; i++) {
      if(i2k[i] != -1) {
        if(prev_pin_i >= 0) {
          int pk = prev_pin_k;
          int ck = i2k[i];
          if(ck > pk) {
            int j;
            for(j = prev_pin_i; j < i; j++) {
              if(kmax[j] < ck) kmax[j] = ESL_MIN(M, ck);
            }
          }
          if(pk > ck) {
            int j;
            for(j = prev_pin_i + 1; j <= i; j++) {
              if(kmin[j] > ck) kmin[j] = ESL_MAX(1, ck);
            }
          }
        }
        prev_pin_i = i;
        prev_pin_k = i2k[i];
      }
    }
  }

  /* Interpolated ramp (brief 26_0628-042). The base scans + bridge block above give
   * every residue in a gap between two consecutive pins the SAME flat band,
   * [prev_pin_k - pad, next_pin_k + pad] -- cost O(gap_len * gap_k) per gap,
   * quadratic when both are large (kmerchain's typical sparse-anchor output;
   * see brief 26_0628-041). Replace that flat band on interior rows of long gaps
   * with a band centered on the linear interpolation between the two
   * bracketing pins, widened only by nodepad (the same local-wiggle margin
   * used at a pin itself) -- not by the full model distance spanned by the
   * gap. This assumes the true path's k-trajectory across a long gap tracks
   * close to the diagonal implied by the two pins, which is reasonable for
   * kmerchain anchors specifically because they were chained on near-zero
   * diagonal drift in the first place (p7_Seq2BandsKmerChain's gapc chain
   * cost penalizes exactly the drift this ramp assumes is small).
   *
   * Short gaps are deliberately left flat (untouched): with few interior
   * rows the flat cost is cheap anyway, and a short gap carrying a large
   * model jump is exactly the "one real deletion run concentrated in a
   * single row" case the D-state bridge above exists to protect -- the ramp
   * does not attempt to prove correctness for that case, so it stays out of
   * its way. The gap_len > 2*pad_here threshold ties "short" to the same
   * nodepad tunable already used for the pin margin itself, rather than an
   * arbitrary constant.
   */
  {
    int prev_pin_i = -1, prev_pin_k = -1;
    for(i = 0; i <= L; i++) {
      if(i2k[i] != -1) {
        if(prev_pin_i >= 0) {
          int pi = prev_pin_i, pk = prev_pin_k;
          int ci = i,          ck = i2k[i];
          int gap_len  = ci - pi;
          int pad_here = ESL_MAX(nodepad[pk], nodepad[ck]);
          if (ck != pk && gap_len > 2 * pad_here) {
            /* RAMP_SLACK_ALPHA: empirically bracketed (brief 26_0628-042), not derived
             * from a correctness proof -- pure nodepad-width margin (alpha=0)
             * measurably WORSENED one real panel case (norovirus MT372469.1,
             * --p7band --p7kmerchain: bit score 1358.60 -> 1168.35, a true
             * regression, not just a placement difference), because a long
             * on-diagonal gap can still contain a real, non-uniformly placed
             * indel that pure linear interpolation clips. Growing the margin
             * by alpha per residue of distance from the nearer bracketing pin
             * (capped, since d is bounded by gap_len/2) recovers that case
             * (1358.60 -> 1352.56 at alpha=0.75, within this thread's
             * established ~0.2%-noise tolerance) while still cutting ncells
             * substantially versus the flat default (2.6M vs 4.7M on that same
             * gap -- alpha=0 gave 1.4M but wasn't safe). This is a tunable
             * correctness/speed tradeoff, not a proven-safe bound either --
             * see brief 26_0628-042's summary "Design reasoning" section. brief 26_0628-043:
             * alpha is now a caller-supplied parameter (default 0.75 preserved
             * at every call site; only --p7kmerchain's cmalign call site
             * exposes it as a runtime option, --p7kmerchain-alpha).
             */
            int j;
            for(j = pi + 1; j < ci; j++) {
              double frac = (double) (j - pi) / (double) gap_len;
              double kexp = (double) pk + frac * (double) (ck - pk);
              int d = (j - pi); if (ci - j < d) d = ci - j; /* distance to nearer anchor */
              int slack = (int) (alpha * (double) d);
              int klo = (int) floor(kexp);
              int khi = (int) ceil(kexp);
              if (klo < 0) klo = 0; if (klo > M) klo = M;
              if (khi < 0) khi = 0; if (khi > M) khi = M;
              kmin[j] = ESL_MAX(1, klo - nodepad[klo] - slack);
              kmax[j] = ESL_MIN(M, khi + nodepad[khi] + slack);
            }
          }
        }
        prev_pin_i = i;
        prev_pin_k = i2k[i];
      }
    }
  }

  /* D1 hop-back dilation pass.
   * For each pinned residue i, take min/max of i2k across the 2*hopback+1
   * window of pinned positions centred on i (clipped at trace ends), and
   * widen kmin[i]/kmax[i] to that range. Vit pins are monotone in (i,k),
   * so the window-min is i2k of the (hopback)th preceding pinned position,
   * and the window-max is i2k of the (hopback)th following pinned position.
   * O(L) total. hopback==0 is a no-op.
   */
  if (hopback > 0) {
    int *pin_pos = NULL;   /* [0..npin-1] residue index of each pinned position, in trace order */
    int  npin = 0;
    int  p;
    ESL_ALLOC(pin_pos, sizeof(int) * (L+2));
    for (i = 0; i <= L; i++) {
      if (i2k[i] != -1) { pin_pos[npin++] = i; }
    }
    for (p = 0; p < npin; p++) {
      int lo = p - hopback; if (lo < 0)       lo = 0;
      int hi = p + hopback; if (hi > npin-1)  hi = npin - 1;
      int kw_min = i2k[pin_pos[lo]];
      int kw_max = i2k[pin_pos[hi]];
      int ii = pin_pos[p];
      int new_kmin = ESL_MAX(1, kw_min);
      int new_kmax = ESL_MIN(M, kw_max);
      if (kmin[ii] > new_kmin) kmin[ii] = new_kmin;
      if (kmax[ii] < new_kmax) kmax[ii] = new_kmax;
    }
    free(pin_pos);
  }

  /* get number of cells if wanted.
   * Brief 26_0628-027: accumulate in int64. At genome scale (M~L~1.5e5, and sparse
   * anchors -> wide bands) the true cell count reaches ~1e10, overflowing the
   * int32 accumulator and wrapping to <=0, which would spuriously trip the
   * caller's ncells==0 fallback. ret_ncells stays int (every caller uses it
   * only as a >0/==0 discriminator + debug avg_bw; it never sizes memory), so
   * we saturate the returned value at INT_MAX and emit the true int64 count on
   * a debug line for genome-scale band-size measurement. */
  if(ret_ncells != NULL) {
    int64_t ncells64 = 0;
    for(i = 1; i <= L; i++) ncells64 += (int64_t)(kmax[i] - kmin[i] + 1);
    if(ncells64 > (int64_t) INT_MAX) {
      fprintf(stderr, "#NCELLS64 p7_pins2bands_nodepad L=%d M=%d ncells=%lld (saturated to INT_MAX in int32 ret_ncells)\n",
              L, M, (long long) ncells64);
      *ret_ncells = INT_MAX;
    } else {
      *ret_ncells = (int) ncells64;
    }
  }

  if(ret_kmin != NULL) { *ret_kmin = kmin; } else free(kmin);
  if(ret_kmax != NULL) { *ret_kmax = kmax; } else free(kmax);
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "p7_pins2bands_nodepad() memory error.");
  return status; /* NEVERREACHED */
}

/* brief 26_0628-045/047: small-M gate for kmerchain (also gated the earlier
 * best-window-anchor deriver before its removal by brief 26_0628-066).
 * Root cause (rmark4 MIR2655 and 4 other catastrophic-loss families, all
 * M=84-400): at this M range these divergent structural-RNA test families
 * have too little exact-match identity to the model's argmax consensus for
 * k-mer seeding to work -- often 0-1 true k=10 matches exist anywhere in the
 * whole model, outnumbered by coincidental short-k-mer collisions that
 * chain/cluster together and out-compete the lone true anchor. This is the
 * SAME mechanism brief 26_0628-022 root-caused and gated for the cmsearch fast-anchor
 * path (norovirus's non-repeat coincidental-collision residual); brief 26_0628-045
 * ported brief 26_0628-022's validated M<4,000 threshold here as a hardcoded,
 * default-on constant. Brief 26_0628-047 replaced that with a real, user-facing,
 * OFF-BY-DEFAULT cmalign option (--p7kmerchain-mgate <M>, cm->p7_kmerchain_mgate,
 * 0=off) -- fires (M < threshold) only when the user opts in, same gating
 * logic as before. Below the threshold, the deriver immediately reports
 * "no anchor" (ret_ncells=0), and the caller falls back per brief 26_0628-047's new
 * shared fallback-selection mechanism (default --p7ibv; see cmalign.c). */

#define KMW_BIN         200   /* model-window (bin) width, matches brief 26_0628-023/26_0628-025 B */
#define KMW_TOL         15    /* diagonal-cluster tolerance (brief 26_0628-023/26_0628-025 TOL)    */
#define KMW_MIN_CORRECT 3     /* floor: min on-diagonal hits for a real window     */
#define KMW_DRIFT       25    /* diagonal drift allowed when expanding outward     */
#define KMW_EXP_FLOOR   2     /* min on-diagonal hits to accept an expansion bin   */

static const int kmw_kvals[] = { 10, 15, 20, 25, 30 };  /* all <=31 => uint64-encodable */
#define KMW_NK ((int)(sizeof(kmw_kvals)/sizeof(kmw_kvals[0])))

/* brief 26_0628-046: per-query k>=mink zero-hits signal gate for kmerchain (also
 * gated the earlier best-window-anchor deriver before its removal by brief 26_0628-066).
 * `nrawk` is the raw hit count per k-tier (kmw_kvals order); fires
 * (returns TRUE) iff mink>0 and every tier with k>=mink has zero hits anywhere
 * in the model for this query -- i.e. there is no exact-match content long
 * enough to carry real signal, independent of M. Disabled (mink<=0, default)
 * by construction returns FALSE. See brief 26_0628-046 for the mechanism/rationale. */
static int
kmer_ngate_fires(const int *nrawk, int mink)
{
  int ki, sum = 0;
  if (mink <= 0) return FALSE;
  for (ki = 0; ki < KMW_NK; ki++) if (kmw_kvals[ki] >= mink) sum += nrawk[ki];
  return (sum == 0);
}

typedef struct { uint64_t code; int pos; } kmw_kmer_t; /* a target k-mer occurrence */
typedef struct { int j, t, k; }            kmw_hit_t;  /* model j / target t / k-length */

static int kmw_kmer_cmp(const void *a, const void *b) {
  uint64_t ca = ((const kmw_kmer_t *)a)->code, cb = ((const kmw_kmer_t *)b)->code;
  return (ca > cb) - (ca < cb);
}
/* first index in sorted arr[0..n-1] whose code >= key */
static int kmw_lower_bound(const kmw_kmer_t *arr, int n, uint64_t key) {
  int lo = 0, hi = n;
  while (lo < hi) { int mid = (lo+hi)/2; if (arr[mid].code < key) lo = mid+1; else hi = mid; }
  return lo;
}

/* encode digital k-mer of length k (<=31) at v[pos] (2 bits/residue).
 * *ok=0 if any residue is degenerate (>3, i.e. not a canonical ACGU). */
static uint64_t kmw_encode(const ESL_DSQ *v, int pos, int k, int *ok) {
  uint64_t code = 0; int o;
  for (o = 0; o < k; o++) {
    ESL_DSQ r = v[pos+o];
    if (r > 3) { *ok = 0; return 0; }
    code |= ((uint64_t) r) << (2*o);
  }
  *ok = 1;
  return code;
}

/*****************************************************************
 * Brief 26_0628-027: k-mer seed-and-chain (minimap2/BLAST-style) guide-deriver for
 * genome-scale (>40kb) anchoring.
 *
 * An earlier deriver (brief 26_0628-026, removed by brief 26_0628-066) picked the
 * single best 200nt model window and seeded only from it -- which at genome scale
 * (a) is fooled by repeat-driven false diagonals (HSV locks onto an internal
 * repeat block) and (b) cannot span a genome from one anchor (HSV 41% / MPXV 88%
 * coverage). This deriver instead collects ALL k-mer seeds genome-wide and chains
 * them by GLOBAL diagonal/colinearity consistency. A repeat block's seeds are
 * internally consistent with each other but NOT colinear with the sequence-
 * spanning chain elsewhere, so a standard colinear chaining DP structurally
 * out-competes an isolated repeat block with a chain that covers far more of the
 * genome.
 *
 * Reuses the k-mer index infrastructure from that earlier deriver unchanged
 * (kmw_encode / kmw_kmer_cmp / kmw_lower_bound, k in {10,15,20,25,30}) and feeds
 * the winning chain's pins through the SAME unmodified p7_pins2bands_nodepad
 * consumer.
 *****************************************************************/

/* chaining-DP tunables (first-pass, minimap2-style; deliberately un-tuned) */
#define KMC_AVGK      15.0    /* representative k for gap-cost scaling            */
#define KMC_MAX_QGAP  30000   /* max query gap between chained anchors (bounds DP)*/
#define KMC_MAX_DGAP  10000   /* max implied-indel (|diagonal diff|) on one link  */
#define KMC_MAX_ITER  5000    /* max predecessors examined per anchor (DP cap)    */
#define KMC_MIN_ANCHOR 10     /* min merged exact-match length to keep an anchor  */
                               /* brief 26_0628-029: lowered 20->10 (=smallest indexed k) to */
                               /* admit far more anchors on sparse-anchor genomes    */
                               /* (MPXV: 56->79 chain anchors, avg_bw 1223->279.6,   */
                               /* peak RSS 12.95G->4.73G), zero accuracy cost.       */
#define KMC_GAP_LIN   0.01    /* linear gap-cost coefficient                      */
#define KMC_GAP_LOG   0.5     /* log2 gap-cost coefficient                        */
#define KMC_HSV_RLO   9401    /* HSV known false-repeat block, model lo (brief026)*/
#define KMC_HSV_RHI   9600    /*  ... model hi -- for the repeat-trap check       */

/* a raw k-mer seed as a query interval on a diagonal: d = model_j - query_t,
 * query positions [tlo,thi]. Point seeds are merged along a diagonal into
 * maximal exact-match anchors (BLAST/minimap2-style) to collapse the multi-k
 * and per-position seed redundancy before chaining. */
typedef struct { int d, tlo, thi; } kmc_ival_t;

/* sort intervals by diagonal d ascending, then query start tlo ascending. */
static int kmc_ival_cmp(const void *a, const void *b) {
  const kmc_ival_t *x = (const kmc_ival_t *)a, *y = (const kmc_ival_t *)b;
  if (x->d   != y->d)   return (x->d   > y->d)   - (x->d   < y->d);
  return (x->tlo > y->tlo) - (x->tlo < y->tlo);
}

/* sort anchors by query position t ascending, tiebreak model j ascending. */
static int kmc_seed_cmp(const void *a, const void *b) {
  const kmw_hit_t *x = (const kmw_hit_t *)a, *y = (const kmw_hit_t *)b;
  if (x->t != y->t) return (x->t > y->t) - (x->t < y->t);
  return (x->j > y->j) - (x->j < y->j);
}

/* brief 26_0628-035: temporary peak-RSS attribution instrumentation, gated by
 * BRIEF035_MEMPOINT (mirrors P135B_FB_INSTRUMENT convention). Reads
 * /proc/self/status VmRSS. Revert before finishing if not worth keeping. */
static long
brief035_rss_kb(void)
{
  FILE *fp = fopen("/proc/self/status", "r");
  char line[256];
  long rss = -1;
  if (fp == NULL) return -1;
  while (fgets(line, sizeof(line), fp) != NULL) {
    if (strncmp(line, "VmRSS:", 6) == 0) { sscanf(line+6, "%ld", &rss); break; }
  }
  fclose(fp);
  return rss;
}

/* Function: p7_Seq2BandsKmerChain()
 * Date:     Brief 26_0628-027, 2026-07-03
 *
 * Purpose:  Derive p7 bands from a genome-wide k-mer seed-and-chain instead of
 *           a single best window. Collect all k-mer seeds (model vs target),
 *           chain them by global colinearity/diagonal consistency (O(N^2)-with-
 *           bounded-lookback DP), backtrack the single best chain, emit its
 *           seeds as multi-segment sparse pins, then the UNMODIFIED
 *           p7_pins2bands_nodepad.
 *
 * Args:     cm         - covariance model (uses cm->fp7 for consensus, cm->p7bpad)
 *           errbuf     - for error messages
 *           dsq        - digital target sequence, 1..L
 *           L          - length of dsq
 *           nodepad    - [0..M] per-node pad array, or NULL for uniform cm->p7bpad
 *           do_trunc   - brief 26_0628-033: CM_ALIGN_TRUNC flag, mirrors p7_Seq2BandsWV's
 *                        do_trunc argument for signature-shape consistency. Unused
 *                        internally: this deriver chains exact k-mer seeds against
 *                        the raw consensus, with no begin/end-anywhere (Tgm) score
 *                        DP and no glocal-only boundary assumption for do_trunc to
 *                        relax.
 *           ret_i2k/ret_kmin/ret_kmax    - RETURN band arrays (caller frees)
 *           ret_ncells - RETURN total banded cells (saturated int; 0 => no
 *                        usable chain, caller falls back to unbanded Viterbi).
 *
 * Return:   eslOK on success (including the ncells=0 "no chain" case).
 */
int
p7_Seq2BandsKmerChain(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, int *nodepad,
                      int do_trunc,
                      int **ret_i2k, int **ret_kmin, int **ret_kmax, int *ret_ncells,
                      double *ret_a_s, double *ret_b_s, double *ret_bd_s)
{
  int status = eslOK;
  (void) do_trunc; /* brief 26_0628-033: no-op, see function header comment */
  int M = cm->fp7->M;
  int K = cm->abc->K;
  ESL_DSQ    *cons  = NULL;         /* model consensus, digital, cons[1..M]   */
  kmc_ival_t *ival  = NULL;         /* raw seed intervals (all k), pre-merge  */
  int         nival = 0, ivalcap = 0;
  kmw_hit_t  *seeds = NULL;         /* merged maximal-exact-match anchors      */
  int         nseed = 0, seedcap = 0;
  int         nrawk[KMW_NK];        /* raw seed count per k (instrumentation)  */
  kmw_kmer_t *idx   = NULL;         /* per-k target k-mer index               */
  float      *f     = NULL;         /* best chain score ending at anchor i     */
  int        *pre   = NULL;         /* chain predecessor                      */
  int        *chain = NULL;         /* backtracked chain (anchor indices)      */
  int        *i2k   = NULL, *kmin = NULL, *kmax = NULL;
  int        *local_nodepad = NULL;
  int ki, j, x, i;
  int min_anchor;                   /* brief 26_0628-078: KMC_MIN_ANCHOR or env override */
  /* brief 26_0628-059 (extended 2026-07-11, brief 190 follow-up): optional
   * stage a/b/bd timing -- a = seed finding (raw hits + merge), b = colinear
   * chaining DP + backtrack + pin emission, bd = converting the winning
   * chain's pins into HMM bands (the p7_pins2bands_nodepad call). If we bail
   * before reaching the chaining DP (M-gate/N-gate/no anchors), all elapsed
   * time is attributed to stage a. If we reach the chaining DP but bail
   * before pins2bands (alloc failure), the remainder is attributed to b and
   * bd is 0. */
  struct timespec _stagea_t0, _stageb_t0, _stagebd_t0, _stage_texit;
  int _stageb_t0_set = FALSE, _stagebd_t0_set = FALSE;
  if (ret_a_s != NULL) clock_gettime(CLOCK_MONOTONIC, &_stagea_t0);

  *ret_i2k = NULL; *ret_kmin = NULL; *ret_kmax = NULL; *ret_ncells = 0;
  if (cm->p7_kmerchain_mgate > 0 && M < cm->p7_kmerchain_mgate) {
    fprintf(stderr, "#KMERCHAIN L=%d M=%d gated=small-M (M<%d): falling back to unbanded\n", L, M, cm->p7_kmerchain_mgate);
    status = eslOK;   /* brief 26_0628-045/047 opt-in small-M gate; caller falls back per brief 26_0628-047's fallback mechanism */
    goto CLEANUP;     /* brief 26_0628-059: route through CLEANUP so stage a/b timing (if requested) still gets filled in */
  }
  for (ki = 0; ki < KMW_NK; ki++) nrawk[ki] = 0;
  if (getenv("BRIEF035_MEMPOINT") != NULL)
    fprintf(stderr, "#MEMPOINT before_kmerchain L=%d M=%d rss_kb=%ld\n", L, M, brief035_rss_kb());

  /* 1. model consensus (argmax match emission per node) */
  ESL_ALLOC(cons, sizeof(ESL_DSQ) * (M+2));
  cons[0] = eslDSQ_SENTINEL; cons[M+1] = eslDSQ_SENTINEL;
  for (j = 1; j <= M; j++) {
    int argmax = 0; float best = cm->fp7->mat[j][0];
    for (x = 1; x < K; x++) if (cm->fp7->mat[j][x] > best) { best = cm->fp7->mat[j][x]; argmax = x; }
    cons[j] = (ESL_DSQ) argmax;
  }

  /* 2. collect ALL k-mer seeds genome-wide as (diagonal, query-interval) records.
   *    A raw genome-scale seed count is ~1e6 (measured: HSV 833k), dominated by
   *    per-position and multi-k redundancy -- far too many for an O(N^2) chain.
   *    We collect each seed as a query interval on its diagonal, then merge
   *    contiguous/overlapping intervals per diagonal into maximal exact-match
   *    anchors, which collapses the redundancy by ~1-2 orders of magnitude and
   *    yields chain-anchors that are colinear match blocks. */
  for (ki = 0; ki < KMW_NK; ki++) {
    int k = kmw_kvals[ki], ntgt = 0, t, p;
    if (k > M || k > L) continue;
    ESL_ALLOC(idx, sizeof(kmw_kmer_t) * (L - k + 1));
    for (t = 1; t <= L - k + 1; t++) {
      int ok; uint64_t code = kmw_encode(dsq, t, k, &ok);
      if (ok) { idx[ntgt].code = code; idx[ntgt].pos = t; ntgt++; }
    }
    qsort(idx, ntgt, sizeof(kmw_kmer_t), kmw_kmer_cmp);
    for (j = 1; j <= M - k + 1; j++) {
      int ok; uint64_t code = kmw_encode(cons, j, k, &ok);
      if (!ok) continue;
      for (p = kmw_lower_bound(idx, ntgt, code); p < ntgt && idx[p].code == code; p++) {
        int tt = idx[p].pos;
        if (nival == ivalcap) {
          int newcap = ivalcap ? ivalcap*2 : 8192;
          void *tmp = realloc(ival, sizeof(kmc_ival_t) * newcap);
          if (tmp == NULL) { status = eslEMEM; goto ERROR; }
          ival = tmp; ivalcap = newcap;
        }
        ival[nival].d = j - tt; ival[nival].tlo = tt; ival[nival].thi = tt + k - 1;
        nival++; nrawk[ki]++;
      }
    }
    free(idx); idx = NULL;
  }

  /* merge intervals per diagonal into maximal exact-match anchors. Sort by
   * (diagonal, tlo); an anchor extends while the next interval on the same
   * diagonal overlaps or abuts (tlo <= cur_thi + 1). Keep anchors of length
   * >= KMC_MIN_ANCHOR to drop isolated short (mostly spurious) matches. */
  /* brief 26_0628-078: optional source-level raise of the anchor-admission
   * floor (env KMC_ANCHOR_FLOOR, unset => compiled KMC_MIN_ANCHOR=10). Tests
   * whether requiring a longer exact match to be an anchor AT ALL beats the
   * per-link KMC_MIN_OPEN floor: it suppresses the same chance-collision
   * anchors, but unlike the per-link floor it also removes them as
   * same-diagonal extension pins, which is where 26_0628-029's MPXV band-width
   * win came from. */
  {
    char *afl_env = getenv("KMC_ANCHOR_FLOOR");
    min_anchor = (afl_env != NULL && atoi(afl_env) > 0) ? atoi(afl_env) : KMC_MIN_ANCHOR;
  }
  if (nival > 0) {
    qsort(ival, nival, sizeof(kmc_ival_t), kmc_ival_cmp);
    int r = 0;
    while (r < nival) {
      int d = ival[r].d, tlo = ival[r].tlo, thi = ival[r].thi;
      int s = r + 1;
      while (s < nival && ival[s].d == d && ival[s].tlo <= thi + 1) {
        if (ival[s].thi > thi) thi = ival[s].thi;
        s++;
      }
      int len = thi - tlo + 1;
      if (len >= min_anchor) {
        if (nseed == seedcap) {
          int newcap = seedcap ? seedcap*2 : 4096;
          void *tmp = realloc(seeds, sizeof(kmw_hit_t) * newcap);
          if (tmp == NULL) { status = eslEMEM; goto ERROR; }
          seeds = tmp; seedcap = newcap;
        }
        seeds[nseed].j = tlo + d; seeds[nseed].t = tlo; seeds[nseed].k = len;
        nseed++;
      }
      r = s;
    }
  }
  fprintf(stderr, "#KMERCHAIN L=%d M=%d nraw=%d perk=[k10:%d k15:%d k20:%d k25:%d k30:%d] nanchor=%d (minlen=%d)\n",
          L, M, nival, nrawk[0], nrawk[1], nrawk[2], nrawk[3], nrawk[4], nseed, min_anchor);

  /* brief 26_0628-046: per-query k>=mink zero-hits signal gate (unvalidated, default
   * off -- see kmer_ngate_fires() header comment). Checked right after the raw
   * hit collection/perk breakdown above, independent of the M-gate above. */
  if (kmer_ngate_fires(nrawk, cm->p7_kmerchain_mink)) {
    fprintf(stderr, "#KMERCHAIN L=%d M=%d gated=zero-hits (k>=%d finds no exact match anywhere in model): "
                    "falling back to unbanded\n", L, M, cm->p7_kmerchain_mink);
    status = eslOK; goto CLEANUP;
  }

  /* brief 26_0628-045: dump every pre-chain anchor (not just the winning chain), gated
   * by BRIEF045_SEEDDUMP (silent no-op by default, same convention as brief
   * 041's BRIEF041_CHAINDUMP). Lets us see whether a given model region has
   * ANY candidate anchor at all, vs has one that the chaining DP rejected. */
  if (getenv("BRIEF045_SEEDDUMP") != NULL) {
    int si;
    for (si = 0; si < nseed; si++)
      fprintf(stderr, "#KMERCHAIN_SEED L=%d M=%d idx=%d/%d model=[%d,%d] query=[%d,%d] len=%d\n",
              L, M, si, nseed, seeds[si].j, seeds[si].j + seeds[si].k - 1,
              seeds[si].t, seeds[si].t + seeds[si].k - 1, seeds[si].k);
  }
  if (getenv("BRIEF035_MEMPOINT") != NULL)
    fprintf(stderr, "#MEMPOINT after_seedmerge L=%d nraw=%d nanchor=%d rss_kb=%ld\n", L, nival, nseed, brief035_rss_kb());
  free(ival); ival = NULL;

  if (nseed == 0) {   /* no anchors: caller falls back to unbanded Viterbi */
    fprintf(stderr, "#KMERCHAIN L=%d M=%d chain=NONE (no anchors)\n", L, M);
    status = eslOK; goto CLEANUP;
  }

  /* 3. colinear chaining DP. Sort anchors by query t (tiebreak model j). For
   *    anchor i, the best predecessor j has t_j<t_i and model_j<model_i
   *    (strictly colinear); link score = f[j] + overlap-adjusted match(i) -
   *    gap_cost, where gap_cost penalizes the implied indel |diag_i - diag_j|
   *    (minimap2-style). Lookback is bounded by KMC_MAX_QGAP (anchors sorted by
   *    t => break once exceeded) and by KMC_MAX_ITER predecessors examined, so
   *    the DP is O(N * min(W, MAX_ITER)). */
  if (ret_a_s != NULL) { clock_gettime(CLOCK_MONOTONIC, &_stageb_t0); _stageb_t0_set = TRUE; }
  qsort(seeds, nseed, sizeof(kmw_hit_t), kmc_seed_cmp);
  ESL_ALLOC(f,   sizeof(float) * nseed);
  ESL_ALLOC(pre, sizeof(int)   * nseed);
  int   best_end = 0;
  float best_sc  = -1.0;
  long   nlink    = 0;   /* DP work counter (instrumentation) */
  /* brief 26_0628-078: soft, env-gated gap-cost knobs on the colinear-chaining
   * DP link cost below. All no-ops when unset (gapc arithmetic is then
   * bit-identical to the pre-078 expression).
   *
   *   KMC_GAP_OPEN=<x>   affine "cost to open a new diagonal": a constant x
   *                      charged once per dgap>0 link, so a chance anchor of
   *                      length k only pays for itself when k > x -- but a
   *                      short anchor may still open a diagonal when the
   *                      alternative is far worse, so legitimate excursions
   *                      (and their pins) are never lost outright. (An affine
   *                      gap cost; the existing terms are the extension part.)
   *   KMC_SPAN_LOG=<c>   span-aware term c*log2(max(qgap,rgap)) per dgap>0
   *                      link -- the gap cost otherwise ignores absolute span
   *                      entirely, which is what makes two cheap hops beat one
   *                      real jump (26_0430's peer message, section 2).
   *   KMC_AVGK_OVERRIDE / KMC_GAP_LIN_OVERRIDE / KMC_GAP_LOG_OVERRIDE
   *                      re-tune the three existing gap-cost constants; the
   *                      compiled defaults were never tuned for wrong-diagonal
   *                      suppression (KMC_AVGK=15 no longer "represents" the
   *                      effective seed length, which is 10 since 26_0628-029). */
  char  *gopen_env = getenv("KMC_GAP_OPEN");
  char  *spanl_env = getenv("KMC_SPAN_LOG");
  char  *avgk_env  = getenv("KMC_AVGK_OVERRIDE");
  char  *glin_env  = getenv("KMC_GAP_LIN_OVERRIDE");
  char  *glog_env  = getenv("KMC_GAP_LOG_OVERRIDE");
  double gap_open  = (gopen_env != NULL) ? atof(gopen_env) : 0.0;
  double span_log  = (spanl_env != NULL) ? atof(spanl_env) : 0.0;
  double avgk      = (avgk_env  != NULL) ? atof(avgk_env)  : KMC_AVGK;
  double glin      = (glin_env  != NULL) ? atof(glin_env)  : KMC_GAP_LIN;
  double glog      = (glog_env  != NULL) ? atof(glog_env)  : KMC_GAP_LOG;
  int    span_on   = (spanl_env != NULL && span_log != 0.0);
  for (i = 0; i < nseed; i++) {
    int   ti = seeds[i].t, ji = seeds[i].j, di = ji - ti;
    float fi = (float) seeds[i].k;   /* base weight = chain-start contribution */
    int   pi = -1;
    int   jj, iter = 0;
    for (jj = i-1; jj >= 0; jj--) {
      int qgap = ti - seeds[jj].t;
      if (qgap > KMC_MAX_QGAP) break;   /* sorted by t: nothing earlier is closer */
      if (++iter > KMC_MAX_ITER) break; /* bound predecessors examined per anchor */
      if (qgap <= 0) continue;          /* need strictly increasing query pos */
      int rgap = ji - seeds[jj].j;
      if (rgap <= 0) continue;          /* need strictly increasing model pos (colinear) */
      int dgap = di - (seeds[jj].j - seeds[jj].t); if (dgap < 0) dgap = -dgap;
      if (dgap > KMC_MAX_DGAP) continue;
      nlink++;
      int   mn    = qgap < rgap ? qgap : rgap;
      int   match = mn < seeds[i].k ? mn : seeds[i].k;   /* overlap-adjusted */
      float gapc  = (dgap > 0) ? (gap_open + glin*avgk*(double)dgap + glog*log2((double)dgap)
                                  + (span_on ? span_log*log2((double)(qgap > rgap ? qgap : rgap)) : 0.0)) : 0.0;
      float sc    = f[jj] + (float) match - gapc;
      if (sc > fi) { fi = sc; pi = jj; }
    }
    f[i] = fi; pre[i] = pi;
    if (fi > best_sc) { best_sc = fi; best_end = i; }
  }

  /* 4. backtrack the single best chain */
  ESL_ALLOC(chain, sizeof(int) * nseed);
  int nc = 0;
  for (i = best_end; i != -1; i = pre[i]) chain[nc++] = i;   /* reversed (end->start) */
  /* chain[0..nc-1] is end..start; reverse to start..end */
  for (i = 0; i < nc/2; i++) { int tmp = chain[i]; chain[i] = chain[nc-1-i]; chain[nc-1-i] = tmp; }

  /* 5. chain diagnostics + HSV repeat-trap check (brief 26_0628-027 central question) */
  {
    int cs = chain[0], ce = chain[nc-1];
    int mlo = seeds[cs].j, mhi = seeds[ce].j + seeds[ce].k - 1;
    int qlo = seeds[cs].t, qhi = seeds[ce].t + seeds[ce].k - 1;
    int nrep = 0, c;   /* chain anchors OVERLAPPING the HSV repeat block */
    for (c = 0; c < nc; c++) {
      int aj = seeds[chain[c]].j, ajx = seeds[chain[c]].j + seeds[chain[c]].k - 1;
      if (aj <= KMC_HSV_RHI && ajx >= KMC_HSV_RLO) nrep++;
    }
    fprintf(stderr, "#KMERCHAIN L=%d M=%d chain_nanchor=%d score=%.1f model_span=[%d,%d] query_span=[%d,%d] "
                    "repeat_anchors_in[%d,%d]=%d nlink=%ld\n",
            L, M, nc, best_sc, mlo, mhi, qlo, qhi, KMC_HSV_RLO, KMC_HSV_RHI, nrep, nlink);
    /* brief 26_0628-041 (ported from infernal-brief040-stagetime): per-anchor chain dump,
     * gated by BRIEF041_CHAINDUMP (silent no-op by default, mirrors
     * BRIEF035_MEMPOINT convention). */
    if (getenv("BRIEF041_CHAINDUMP") != NULL) {
      for (c = 0; c < nc; c++) {
        int aj  = seeds[chain[c]].j,  ajx = seeds[chain[c]].j + seeds[chain[c]].k - 1;
        int at  = seeds[chain[c]].t,  atx = seeds[chain[c]].t + seeds[chain[c]].k - 1;
        int qgap = -1, mgap = -1, dgap = -1;
        if (c+1 < nc) {
          int nj = seeds[chain[c+1]].j, nt = seeds[chain[c+1]].t;
          qgap = nt - atx - 1;
          mgap = nj - ajx - 1;
          dgap = (nj - nt) - (aj - at);
        }
        fprintf(stderr, "#KMERCHAIN_ANCHOR L=%d M=%d idx=%d/%d model=[%d,%d] query=[%d,%d] len=%d "
                        "gap_to_next: query=%d model=%d diagdelta=%d\n",
                L, M, c, nc, aj, ajx, at, atx, seeds[chain[c]].k, qgap, mgap, dgap);
      }
    }
    if (getenv("BRIEF035_MEMPOINT") != NULL)
      fprintf(stderr, "#MEMPOINT after_chainDP L=%d nseed=%d nlink=%ld rss_kb=%ld\n", L, nseed, nlink, brief035_rss_kb());
  }

  /* 6. emit multi-segment pins from every seed in the winning chain */
  ESL_ALLOC(i2k, sizeof(int) * (L+1));
  esl_vec_ISet(i2k, L+1, -1);
  {
    int c, o;
    for (c = 0; c < nc; c++) {
      int jj = seeds[chain[c]].j, tt = seeds[chain[c]].t, kk = seeds[chain[c]].k;
      for (o = 0; o < kk; o++) {
        int T = tt + o, J = jj + o;
        if (T >= 1 && T <= L && J >= 1 && J <= M && i2k[T] == -1) i2k[T] = J;
      }
    }
  }

  /* report pin span (reach) -- directly comparable to brief 26_0628-026's cover= */
  {
    int npin = 0, tmin = L+1, tmax = 0;
    for (i = 1; i <= L; i++) if (i2k[i] != -1) { npin++; if (i < tmin) tmin = i; if (i > tmax) tmax = i; }
    fprintf(stderr, "#KMERCHAIN L=%d M=%d npins=%d pin_tspan=[%d,%d] cover=%.3f\n",
            L, M, npin, (npin? tmin:0), tmax, (double) npin / (double) L);
    if (getenv("BRIEF035_MEMPOINT") != NULL)
      fprintf(stderr, "#MEMPOINT after_pinemission L=%d npins=%d rss_kb=%ld\n", L, npin, brief035_rss_kb());
  }

  /* 7. pins -> bands via the UNMODIFIED consumer */
  if (ret_bd_s != NULL) { clock_gettime(CLOCK_MONOTONIC, &_stagebd_t0); _stagebd_t0_set = TRUE; }
  if (nodepad == NULL) {
    int k2;
    ESL_ALLOC(local_nodepad, sizeof(int) * (M+1));
    for (k2 = 0; k2 <= M; k2++) local_nodepad[k2] = cm->p7bpad;
    nodepad = local_nodepad;
  }
  /* brief 26_0628-043: --p7kmerchain-alpha overrides brief 26_0628-042's ramp-slack alpha
   * (default 0.75, cm->p7_kmerchain_ramp_alpha initialized in cm.c). */
  if ((status = p7_pins2bands_nodepad(i2k, errbuf, L, M, nodepad, 0, cm->p7_kmerchain_ramp_alpha, &kmin, &kmax, ret_ncells)) != eslOK) goto ERROR;

  *ret_i2k = i2k; *ret_kmin = kmin; *ret_kmax = kmax;
  i2k = kmin = kmax = NULL;   /* handed off to caller */

 CLEANUP:
  if (idx)           free(idx);
  if (ival)          free(ival);
  if (seeds)         free(seeds);
  if (f)             free(f);
  if (pre)           free(pre);
  if (chain)         free(chain);
  if (cons)          free(cons);
  if (local_nodepad) free(local_nodepad);
  if (i2k)           free(i2k);
  if (kmin)          free(kmin);
  if (kmax)          free(kmax);
  if (getenv("BRIEF035_MEMPOINT") != NULL)
    fprintf(stderr, "#MEMPOINT after_kmerchain_return L=%d rss_kb=%ld\n", L, brief035_rss_kb());
  if (ret_a_s != NULL) {
    clock_gettime(CLOCK_MONOTONIC, &_stage_texit);
    if (_stagebd_t0_set) {
      *ret_a_s  = (_stageb_t0.tv_sec  - _stagea_t0.tv_sec)  + (_stageb_t0.tv_nsec  - _stagea_t0.tv_nsec)  / 1e9;
      *ret_b_s  = (_stagebd_t0.tv_sec - _stageb_t0.tv_sec)  + (_stagebd_t0.tv_nsec - _stageb_t0.tv_nsec)  / 1e9;
      if (ret_bd_s != NULL)
        *ret_bd_s = (_stage_texit.tv_sec - _stagebd_t0.tv_sec) + (_stage_texit.tv_nsec - _stagebd_t0.tv_nsec) / 1e9;
    } else if (_stageb_t0_set) {
      *ret_a_s = (_stageb_t0.tv_sec - _stagea_t0.tv_sec) + (_stageb_t0.tv_nsec - _stagea_t0.tv_nsec) / 1e9;
      *ret_b_s = (_stage_texit.tv_sec - _stageb_t0.tv_sec) + (_stage_texit.tv_nsec - _stageb_t0.tv_nsec) / 1e9;
      if (ret_bd_s != NULL) *ret_bd_s = 0.;
    } else {
      *ret_a_s = (_stage_texit.tv_sec - _stagea_t0.tv_sec) + (_stage_texit.tv_nsec - _stagea_t0.tv_nsec) / 1e9;
      *ret_b_s = 0.;
      if (ret_bd_s != NULL) *ret_bd_s = 0.;
    }
  }
  return status;

 ERROR:
  if (errbuf != NULL) snprintf(errbuf, eslERRBUFSIZE, "p7_Seq2BandsKmerChain() memory error");
  goto CLEANUP;
}

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
 *           input sequence (1..L). Complements cp9_BackwardP7B().  
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
 *           This function requires that local EL states are turned off, if 
 *           they're not it returns. A separate function should exist to align
 *           with EL states on, I don't think it's worth the extra computations to
 *           make 1 function that handles both.
 *
 *           See additional notes in cp9_Forward() "Purpose" section.
 *
 * Args:
 *           cp9       - the CP9 HMM
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
cp9_ForwardP7B(CP9_t *cp9, char *errbuf, CP9_MX *mx, ESL_DSQ *dsq, int L, int *kmin, int *kmax, float *ret_sc)
{
  int          status;
  int          i;           /* j-W: position in the subsequence                             */
  int          k;           /* CP9 HMM node position                                        */
  int        **mmx;         /* DP matrix for match  state scores [0..1][0..cp9->M]          */
  int        **imx;         /* DP matrix for insert state scores [0..1][0..cp9->M]          */
  int        **dmx;         /* DP matrix for delete state scores [0..1][0..cp9->M]          */
  int        **elmx;        /* DP matrix for EL state scores [0..1][0..cp9->M]              */
  int         *erow;        /* end score for each position [0..1]                           */
  int          M;           /* cp9->M, query length, number of consensus nodes of model */
  int          kp, kn, kx, kpcur, kpprv;

  /* Contract checks */
  if(cp9 == NULL)                      ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_ForwardP7B, cp9 is NULL.\n");
  if(dsq == NULL)                      ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_ForwardP7B, dsq is NULL.");
  if(mx == NULL)                       ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_ForwardP7B, mx is NULL.\n");
  if(mx->M != cp9->M)                 ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_ForwardP7B, mx->M != cp9->M.\n");
  if(kmin == NULL || kmax == NULL)     ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_ForwardP7B, kmin and/or kmax == NULL.\n");
  if(cp9->flags & CPLAN9_EL)           ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_ForwardP7B, cp9 EL flag up.\n");
    
  M = cp9->M;

  int const *tsc = cp9->otsc; /* ptr to efficiently ordered transition scores           */

  /* gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits. */

  /* Rearrange DP matrix for this seq */
  if((status = GrowCP9Matrix(mx, errbuf, L, M, kmin, kmax, &mmx, &imx, &dmx, &elmx, &erow)) != eslOK) return status;
  ESL_DPRINTF2(("#DEBUG: cp9_ForwardP7B(): CP9 matrix size: %.8f Mb rows: %d.\n", mx->size_Mb, mx->rows));
  
  /* Initialization of the zero row. */
  mmx[0][0] = 0;      /* M_0 is state B, and everything starts in B */
  imx[0][0] = -INFTY; /* I_0 is state N, can't get here without emitting*/
  dmx[0][0] = -INFTY; /* D_0 doesn't exist. */
  elmx[0][0]= -INFTY; /* can't go from B to EL state */
  erow[0]   = -INFTY;   
  
  /* Because there's a D state for every node 1..M, 
     dmx[0][k] is possible for all k 1..M (if it's within kmin[0]..kmax[0]) */
  int sc;
  i = 0;
  kn = ESL_MAX(1, kmin[0]);
  kp = kn - kmin[0];
  for (k = kn; k <= kmax[0]; k++, kp++) { 
    assert(kp >= 0);
    mmx[0][kp] = imx[0][kp] = elmx[0][kp] = -INFTY;      /* need seq to get here */
    sc = -INFTY;
    if(kp > 0) { /* if kp == 0, kp-1 is outside the bands */
      sc = ILogsum(ILogsum(mmx[0][kp-1] + CP9TSC(cp9O_MD,k-1),
			   imx[0][kp-1] + CP9TSC(cp9O_ID,k-1)),
		   dmx[0][kp-1] + CP9TSC(cp9O_DD,k-1));
    }
    dmx[0][kp] = sc;
  }
  /* We can do a full parse through all delete states. */
  erow[0] = -INFTY;
  if(INBAND(0, M)) { erow[0] = dmx[0][M] + CP9TSC(cp9O_DM,M); }
     
  /*****************************************************************
   * The main loop: scan the sequence from position 1 to L.
   *****************************************************************/
  /* Recursion. */
  
  for (i = 1; i <= L; i++) 
    { 
      int const *isc = cp9->isc[dsq[i]];
      int const *msc = cp9->msc[dsq[i]];
      int endsc     = -INFTY;
      int sc;

      if(kmin[i] == 0) { 
	mmx[i][0]  = -INFTY;
	dmx[i][0]  = -INFTY;  /*D_0 is non-existent*/
	elmx[i][0] = -INFTY;  /*no EL state for node 0 */
	sc = ILogsum(ILogsum(mmx[i-1][0] + CP9TSC(cp9O_MI,0),
			     imx[i-1][0] + CP9TSC(cp9O_II,0)),
		     dmx[i-1][0] + CP9TSC(cp9O_DI,0));
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
      for (kpcur = 0; kpcur < ESL_MIN(kn-kmin[i], kmax[i]-kmin[i]+1); kpcur++) mmx[i][kpcur] = -INFTY; /* impossible to reach these guys */
      for (kpcur = ESL_MAX(0,kx-kmin[i]+1); kpcur <= kmax[i]-kmin[i]; kpcur++) mmx[i][kpcur] = -INFTY; /* impossible to reach these guys */

      kpcur = kn - kmin[i]; /* unnec, loop above ends with this */
      kpprv = kn - kmin[i-1];
      for (k = kn; k <= kx; k++, kpcur++, kpprv++) {
	/*printf("M i: %d kpprv: %d\n", i, kpprv);*/
	assert((kpprv-1) >= 0);

	sc = ILogsum(ILogsum(mmx[i-1][kpprv-1] + CP9TSC(cp9O_MM,k-1),
			     imx[i-1][kpprv-1] + CP9TSC(cp9O_IM,k-1)),
		     dmx[i-1][kpprv-1] + CP9TSC(cp9O_DM,k-1));

	/* FIX ME: inefficient! check B->M_K transition */
	if(INBAND(i-1, 0)) { /* if i-1 is in k == 0's band */
	  assert(kmin[(i-1)] == 0);
	  if(mmx[i-1][0] != -INFTY)
	   sc = ILogsum(sc, mmx[i-1][0] + CP9TSC(cp9O_BM,k));
	}

	if(sc != -INFTY) {
	  mmx[i][kpcur] = sc + msc[k];
	  /* E state update */
	  endsc = ILogsum(endsc, mmx[i][kpcur] + CP9TSC(cp9O_ME,k));
	}
	else { 
	  mmx[i][kpcur] = -INFTY;
	  /* don't update E state */
	}
	///printf("mmx[i:%4d][k:%4d(%4d)]: %d\n", i, k, kpcur, mmx[i][kpcur]);
      }

      /* insert state*/
      kn = ESL_MAX(kmin[i], kmin[i-1]);
      kx = ESL_MIN(kmax[i], kmax[i-1]);

      for (kpcur = 0; kpcur < ESL_MIN(kn-kmin[i], kmax[i]-kmin[i]+1); kpcur++) imx[i][kpcur] = -INFTY; /* impossible to reach these guys */
      for (kpcur = ESL_MAX(0,kx-kmin[i]+1); kpcur <= kmax[i]-kmin[i]; kpcur++) imx[i][kpcur] = -INFTY; /* impossible to reach these guys */

      kpcur = kn - kmin[i]; /* unnec, loop above ends with this */
      kpprv = kn - kmin[i-1];
      for (k = kn; k <= kx; k++, kpcur++, kpprv++) { 
	/*insert state*/
	assert(kpprv >= 0);
	/* HERE, EVENTUALLY IF kmin/kmax differ b/t Match and Inserts: 
	 * only look at match states from k that have i-1 within band */
	/* all insert states from k should have i-1 within band */
	/*printf("I i: %d kpprv: %d\n", i, kpprv);*/
	sc = ILogsum(ILogsum(mmx[i-1][kpprv] + CP9TSC(cp9O_MI,k),
			     imx[i-1][kpprv] + CP9TSC(cp9O_II,k)),
		     dmx[i-1][kpprv] + CP9TSC(cp9O_DI,k));
	if(sc != -INFTY) imx[i][kpcur] = sc + isc[k];
	else             imx[i][kpcur] = -INFTY;
	///printf("imx[i:%4d][k:%4d(%4d)]: %d\n", i, k, kpcur, imx[i][kpcur]);
      }

      /*delete state*/
      kn = kmin[i]+1;

      for (kpcur = 0; kpcur < (kn - kmin[i]); kpcur++) dmx[i][kpcur] = -INFTY; /* impossible to reach these guys */
      kpcur = kn - kmin[i]; /* unnec, loop above ends with this */
      for (k = kn; k <= kmax[i]; k++, kpcur++) { /* should I be adding one for delete off-by-one?? */
	sc = ILogsum(ILogsum(mmx[i][kpcur-1] + CP9TSC(cp9O_MD,k-1),
			     imx[i][kpcur-1] + CP9TSC(cp9O_ID,k-1)),
		     dmx[i][kpcur-1] + CP9TSC(cp9O_DD,k-1));
	dmx[i][kpcur] = sc;
	///printf("dmx[i:%4d][k:%4d(%4d)]: %d\n", i, k, kpcur, dmx[i][kpcur]);
      }
	  /*printf("mmx [jp:%d][%d]: %d\n", jp, k, mmx[j][k]);
	    printf("imx [jp:%d][%d]: %d\n", jp, k, imx[j][k]);
	    printf("dmx [jp:%d][%d]: %d\n", jp, k, dmx[j][k]);
	    printf("elmx[jp:%d][%d]: %d\n", jp, k, elmx[j][k]);*/

      if(INBAND(i, M)) { 
	endsc = ILogsum(ILogsum(endsc, dmx[i][M-kmin[i]] + CP9TSC(cp9O_DM,M)), /* transition from D_M -> end */
			imx[i][M-kmin[i]] + CP9TSC(cp9O_IM,M)); /* transition from I_M -> end */
      }
      /* transition penalty to EL incurred when EL was entered */
      /*printf("endsc: %d\n", endsc);*/

      erow[i] = endsc;
    } /* end loop over end positions i */
  
  *ret_sc = Scorify(erow[L]);
  ESL_DPRINTF1(("#DEBUG: cp9_ForwardP7B() return score: %10.4f\n", Scorify(erow[L])));

  return eslOK;
}



/* Function: cp9_ForwardP7B_OLD_WITH_EL()
 * 
 * Purpose:  Runs the banded Forward dynamic programming algorithm on an
 *           input sequence (1..L). Complements cp9_BackwardP7B().  
 *           The 'P7B' suffix indicates plan 7 HMM derived bands
 *           in the kmin and kmax arrays are applied.  This function
 *           was derived from cp9_Forward(), differences from that
 *           function were introduced solely to impose bands on the 
 *           matrix. 
 *         
 *           Because of the bands, some options that exist to cp9_Forward()
 *           (like be_efficient and do_scan and reporting hits in scan mode) 
 *           are not available here. This function is meant to be used
 *           solely for first stage of a ForwardP7B, BackwardP7B, Posterior type
 *           calculation. 
 *
 *           Also due to bands, only L is passed as seq length, instead
 *           of i0 and j0 (seq start/stop).  This simplifies the application
 *           of bands kmin[i]/kmax[i] refers to residue i.
 *
 *           See additional notes in cp9_Forward() "Purpose" section.
 *
 *           NOTE: EPN, Thu Aug 14 09:54:36 2008
 *           This is the original version written to potentially handle EL
 *           states, but I think omitting them when they're off is significantly
 *           more efficient, so I wrote a competing function called
 *           cp9_ForwardP7B(). Note the included EL dp steps HAVE NOT BEEN
 *           TESTED YET! I just left them here as a starting point if I 
 *           ever want to implement them.
 *
 * Args:
 *           cp9       - the CP9 HMM
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

int
cp9_ForwardP7B_OLD_WITH_EL(CP9_t *cp9, char *errbuf, CP9_MX *mx, ESL_DSQ *dsq, int L, int *kmin, int *kmax, float *ret_sc)
{
  int          status;
  int          i;           /* j-W: position in the subsequence                             */
  int          k;           /* CP9 HMM node position                                        */
  int        **mmx;         /* DP matrix for match  state scores [0..1][0..cp9->M]      */
  int        **imx;         /* DP matrix for insert state scores [0..1][0..cp9->M]      */
  int        **dmx;         /* DP matrix for delete state scores [0..1][0..cp9->M]      */
  int        **elmx;        /* DP matrix for EL state scores [0..1][0..cp9->M]          */
  int         *erow;        /* end score for each position [0..1]                           */
  int          c;           /* counter for EL states                                        */
  int          M;           /* cp9->M, query length, number of consensus nodes of model */
  int          kn, kx, kpcur, kpprv, kpprv_el;

  /* Contract checks */
  if(cp9 == NULL)                      ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_ForwardP7B_OLD_WITH_EL, cm->cp9 is NULL.\n");
  if(dsq == NULL)                      ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_ForwardP7B_OLD_WITH_EL, dsq is NULL.");
  if(mx == NULL)                       ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_ForwardP7B_OLD_WITH_EL, mx is NULL.\n");
  if(mx->M != cp9->M)                  ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_ForwardP7B_OLD_WITH_EL, mx->M != cp9->M.\n");
  if(kmin == NULL || kmax == NULL)     ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_ForwardP7B_OLD_WITH_EL, kmin and/or kmax == NULL.\n");
    
  M = cp9->M;

  int const *tsc = cp9->otsc; /* ptr to efficiently ordered transition scores           */

  /* gamma allocation and initialization.
   * This is a little SHMM that finds an optimal scoring parse
   * of multiple nonoverlapping hits. */

  /* Grow DP matrix if nec, to either 2 rows or L+1 rows (depending on be_efficient), 
   * stays M+1 columns */
  if((status = GrowCP9Matrix(mx, errbuf, L, M, kmin, kmax, &mmx, &imx, &dmx, &elmx, &erow)) != eslOK) return status;
  ESL_DPRINTF2(("#DEBUG: cp9_ForwardP7B_OLD_WITH_EL(): CP9 matrix size: %.8f Mb rows: %d.\n", mx->size_Mb, mx->rows));
  
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
  kn    = ESL_MAX(1, kmin[0]+1);
  kpcur = kn - kmin[0];
  for (k = kn; k <= kmax[0]; k++, kpcur++) { 
    assert(kpcur >= 1);
    mmx[0][kpcur] = imx[0][kpcur] = elmx[0][kpcur] = -INFTY;      /* need seq to get here */
    sc = ILogsum(ILogsum(mmx[0][kpcur-1] + CP9TSC(cp9O_MD,k-1),
			 imx[0][kpcur-1] + CP9TSC(cp9O_ID,k-1)),
		 dmx[0][kpcur-1] + CP9TSC(cp9O_DD,k-1));
    dmx[0][kpcur] = sc;
  }
  /* We can do a full parse through all delete states. */
  erow[0] = -INFTY;
  if(INBAND(0, M)) { erow[0] = dmx[0][M] + CP9TSC(cp9O_DM,M); }
     
  /*****************************************************************
   * The main loop: scan the sequence from position 1 to L.
   *****************************************************************/
  /* Recursion. */
  
  for (i = 1; i <= L; i++) 
    { 
      int const *isc = cp9->isc[dsq[i]];
      int const *msc = cp9->msc[dsq[i]];
      int endsc     = -INFTY;
      int el_selfsc = cp9->el_selfsc;
      int sc;

      if(kmin[i] == 0) { 
	mmx[i][0]  = -INFTY;
	dmx[i][0]  = -INFTY;  /*D_0 is non-existent*/
	elmx[i][0] = -INFTY;  /*no EL state for node 0 */
	sc = ILogsum(ILogsum(mmx[i-1][0] + CP9TSC(cp9O_MI,0),
			     imx[i-1][0] + CP9TSC(cp9O_II,0)),
		     dmx[i-1][0] + CP9TSC(cp9O_DI,0));
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
      for (kpcur = 0; kpcur < ESL_MIN(kn-kmin[i], kmax[i]-kmin[i]+1); kpcur++) mmx[i][kpcur] = -INFTY; /* impossible to reach these guys */
      for (kpcur = ESL_MAX(0,kx-kmin[i]+1); kpcur <= kmax[i]-kmin[i]; kpcur++) mmx[i][kpcur] = -INFTY; /* impossible to reach these guys */

      for (kpcur = 0; kpcur < ESL_MIN(kn-kmin[i], kmax[i]-kmin[i]+1); kpcur++) elmx[i][kpcur] = -INFTY; /* impossible to reach these guys */
      for (kpcur = ESL_MAX(0,kx-kmin[i]+1); kpcur <= kmax[i]-kmin[i]; kpcur++) elmx[i][kpcur] = -INFTY; /* impossible to reach these guys */

      kpcur = kn - kmin[i]; /* unnec, loop above ends with this */
      kpprv = kn - kmin[i-1];
      for (k = kn; k <= kx; k++, kpcur++, kpprv++) {
	/*printf("M i: %d kpprv: %d\n", i, kpprv);*/
	assert((kpprv-1) >= 0);

	sc = ILogsum(ILogsum(mmx[i-1][kpprv-1] + CP9TSC(cp9O_MM,k-1),
			     imx[i-1][kpprv-1] + CP9TSC(cp9O_IM,k-1)),
		     dmx[i-1][kpprv-1] + CP9TSC(cp9O_DM,k-1));

	/* FIX ME: inefficient! check B->M_K transition */
	if(INBAND(i-1, 0)) { /* if i-1 is in k == 0's band */
	  assert(kmin[(i-1)] == 0);
	  if(mmx[i-1][0] != -INFTY) 
	   sc = ILogsum(sc, mmx[i-1][0] + CP9TSC(cp9O_BM,k));
	}

	/* check possibility we came from an EL, if they're valid */
	for(c = 0; c < cp9->el_from_ct[k]; c++) { /* el_from_ct[k] is >= 0 */
	  if(INBAND(i-1, cp9->el_from_idx[k][c])) { 
	    kpprv_el = cp9->el_from_idx[k][c] - kmin[(i-1)];
	    sc = ILogsum(sc, elmx[i-1][kpprv_el]);
	  }
	} /* transition penalty to EL incurred when EL was entered */
	if(sc != -INFTY) { 
	  mmx[i][kpcur] = sc + msc[k];
	  /* E state update */
	  endsc = ILogsum(endsc, mmx[i][kpcur] + CP9TSC(cp9O_ME,k));
	}
	else { 
	  mmx[i][kpcur] = -INFTY;
	  /* don't update E state */
	}
	/*printf("k: %4d mmx[i:%4d][kpcur:%4d]: %d\n", k, i, kpcur, mmx[i][kpcur]);*/

	/* el state */
	sc = -INFTY;
	if((cp9->flags & CPLAN9_EL) && cp9->has_el[k]) /* not all HMM nodes have an EL state (for ex:
								  HMM nodes that map to right half of a MATP_MP) */
	  {
	    sc = mmx[i][kpcur] + CP9TSC(cp9O_MEL,k); /* M_k -> EL_k transition */
	    if(INBAND(i-1, k)) { /* EL_k self-loop: only if k was in band at previous row */
	      kpprv_el = k - kmin[(i-1)];
	      sc = ILogsum(sc, elmx[i-1][kpprv_el] + el_selfsc);
	    }
	  }
	elmx[i][kpcur] = sc;
      }

      /* insert state*/
      kn = ESL_MAX(kmin[i], kmin[i-1]);
      kx = ESL_MIN(kmax[i], kmax[i-1]);

      for (kpcur = 0; kpcur < ESL_MIN(kn-kmin[i], kmax[i]-kmin[i]+1); kpcur++) imx[i][kpcur] = -INFTY; /* impossible to reach these guys */
      for (kpcur = ESL_MAX(0,kx-kmin[i]+1); kpcur <= kmax[i]-kmin[i]; kpcur++) imx[i][kpcur] = -INFTY; /* impossible to reach these guys */

      kpcur = kn - kmin[i]; /* unnec, loop above ends with this */
      kpprv = kn - kmin[i-1];
      for (k = kn; k <= kx; k++, kpcur++, kpprv++) { 
	/*insert state*/
	assert(kpprv >= 0);
	/* HERE, EVENTUALLY IF kmin/kmax differ b/t Match and Inserts: 
	 * only look at match states from k that have i-1 within band */
	/* all insert states from k should have i-1 within band */
	/*printf("I i: %d kpprv: %d\n", i, kpprv);*/
	sc = ILogsum(ILogsum(mmx[i-1][kpprv] + CP9TSC(cp9O_MI,k),
			     imx[i-1][kpprv] + CP9TSC(cp9O_II,k)),
		     dmx[i-1][kpprv] + CP9TSC(cp9O_DI,k));
	if(sc != -INFTY) imx[i][kpcur] = sc + isc[k];
	else             imx[i][kpcur] = -INFTY;
	/*printf("k: %4d imx[i:%4d][kpcur:%4d]: %d\n", k, i, kpcur, imx[i][kpcur]);*/
      }

      /*delete state*/
      kn = kmin[i]+1;

      for (kpcur = 0; kpcur < (kn - kmin[i]); kpcur++) dmx[i][kpcur] = -INFTY; /* impossible to reach these guys */
      kpcur = kn - kmin[i]; /* unnec, loop above ends with this */
      for (k = kn; k <= kmax[i]; k++, kpcur++) { /* should I be adding one for delete off-by-one?? */
	sc = ILogsum(ILogsum(mmx[i][kpcur-1] + CP9TSC(cp9O_MD,k-1),
			     imx[i][kpcur-1] + CP9TSC(cp9O_ID,k-1)),
		     dmx[i][kpcur-1] + CP9TSC(cp9O_DD,k-1));
	dmx[i][kpcur] = sc;
	/*printf("k: %4d dmx[i:%4d][kpcur:%4d]: %d\n", k, i, kpcur, dmx[i][kpcur]);*/
      }
	  /*printf("mmx [jp:%d][%d]: %d\n", jp, k, mmx[j][k]);
	    printf("imx [jp:%d][%d]: %d\n", jp, k, imx[j][k]);
	    printf("dmx [jp:%d][%d]: %d\n", jp, k, dmx[j][k]);
	    printf("elmx[jp:%d][%d]: %d\n", jp, k, elmx[j][k]);*/

      if(INBAND(i, M)) { 
	endsc = ILogsum(ILogsum(endsc, dmx[i][M-kmin[i]] + CP9TSC(cp9O_DM,M)), /* transition from D_M -> end */
			imx[i][M-kmin[i]] + CP9TSC(cp9O_IM,M)); /* transition from I_M -> end */
	for(c = 0; c < cp9->el_from_ct[M+1]; c++) { /* el_from_ct[k] is >= 0 */
	  if(INBAND(i, cp9->el_from_idx[M+1][c])) { 
	    kpprv_el = cp9->el_from_idx[M+1][c] - kmin[i];
	    endsc = ILogsum(endsc, elmx[i][kpprv_el]);
	  }
	}
      }
	/* transition penalty to EL incurred when EL was entered */
      /*printf("endsc: %d\n", endsc);*/

      erow[i] = endsc;
    } /* end loop over end positions i */
  
  *ret_sc = Scorify(erow[L]);
  ESL_DPRINTF1(("#DEBUG: cp9_ForwardP7B_OLD_WITH_EL() return score: %10.4f\n", Scorify(erow[L])));

  return eslOK;
}

/* Function: cp9_BackwardP7B()
 * 
 * Purpose:  Runs the banded Backward dynamic programming algorithm on an
 *           input sequence (1..L). Complements cp9_ForwardP7B().  
 *           The 'P7B' suffix indicates plan 7 HMM derived bands
 *           in the kmin and kmax arrays are applied.  This function
 *           was derived from cp9_Backward(), differences from that
 *           function were introduced solely to impose bands on the 
 *           matrix. 
 *
 *           Because of the bands, some options that exist to cp9_Forward()
 *           (like be_efficient and do_scan and reporting hits in scan mode) 
 *           are not available here. This function is meant to be used
 *           solely for second stage of a Forward, Backward, Posterior type
 *           calculation. 
 *
 *           Also due to bands, only L is passed as seq length, instead
 *           of i0 and j0 (seq start/stop).  This simplifies the application
 *           of bands kmin[i]/kmax[i] refers to residue i.
 *
 *           This function requires that local EL states are turned off, if 
 *           they're not it returns. A separate function should exist to align
 *           with EL states on, I don't think it's worth the extra computations to
 *           make 1 function that handles both.
 *
 *           See additional notes in cp9_Backward() "Purpose" section.
 *
 * Args:     
 *           cp9       - the CP9 HMM
 *           errbuf    - char buffer for error messages
 *           mx        - the matrix, expanded to correct size (if nec), and filled here
 *           dsq       - sequence in digitized form
 *           L         - start of target subsequence (1 for beginning of dsq)
 *           kmin      - [0.1..i..N] minimum k for residue i
 *           kmax      - [0.1..i..N] maximum k for residue i
 *           ret_sc    - RETURN: log P(S|M,bands)/P(S|R), as a bit score, this is B->M[0][0]
 *
 * Returns:  eslOK on success; 
 *           eslEINCOMPAT on contract violation;
 */
int
cp9_BackwardP7B(CP9_t *cp9, char *errbuf, CP9_MX *mx, ESL_DSQ *dsq, int L, int *kmin, int *kmax, float *ret_sc)
{
  int          status;
  int          i;           /*     j-W: position in the subsequence                         */
  int          k;           /* CP9 HMM node position                                        */
  int        **mmx;         /* DP matrix for match  state scores [0..1][0..cp9->M]      */
  int        **imx;         /* DP matrix for insert state scores [0..1][0..cp9->M]      */
  int        **dmx;         /* DP matrix for delete state scores [0..1][0..cp9->M]      */
  int        **elmx;        /* DP matrix for EL state scores [0..1][0..cp9->M]          */
  int         *erow;        /* end score for each position [0..1]                           */
  int          c;           /* counter for EL states */
  int          M;           /* cp9->M */
  int          kpcur, kpprv, kpcur_el;
  int          kprv, kprvn, kprvx;
  int          kn, kx;

  /* Contract checks */
  if(cp9 == NULL)                      ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_BackwardP7B, cp9 is NULL.\n");
  if(dsq == NULL)                      ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_BackwardP7B, dsq is NULL.");
  if(mx == NULL)                       ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_BackwardP7B, mx is NULL.\n");
  if(mx->M != cp9->M)                  ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_BackwardP7B, mx->M != cm->clen.\n");
  /* EL states are handled in the body of this function (unlike cp9_ForwardP7B) */
    
  M = cp9->M;

  int const *tsc = cp9->otsc; /* ptr to efficiently ordered transition scores           */

  /* Rearrange DP matrix for this seq */
  if((status = GrowCP9Matrix(mx, errbuf, L, M, kmin, kmax, &mmx, &imx, &dmx, &elmx, &erow)) != eslOK) return status;
  ESL_DPRINTF2(("#DEBUG: cp9_BackwardP7B(): CP9 matrix size: %.8f Mb rows: %d.\n", mx->size_Mb, mx->rows));

  /* Initialization of the L row. */
  i = L;

  /*******************************************************************
   * 0 Handle EL, looking at EL_k->E for all valid k.
   * we're going backwards so we have to work out of order, we could get 
   * around this by storing the nodes each EL goes TO in an el_to_ct[] vec. */
  /* init to -INFTY */
  kpcur = 0;
  for (k = kmin[i]; k <= kmax[i]; k++, kpcur++) elmx[i][kpcur] = -INFTY;
  if(cp9->flags & CPLAN9_EL)
    {
      for(c = 0; c < cp9->el_from_ct[cp9->M+1]; c++) /* el_from_ct[cp9->M+1] holds # ELs that can go to END */
	if(INBAND(i, cp9->el_from_idx[M+1][c])) { 
	  kpcur_el = cp9->el_from_idx[M+1][c] - kmin[i];
	  elmx[i][kpcur_el] = 0.;
	}
    }
  /*******************************************************************/

  /* elmx[cur][cp9->M] is either 0 (if EL_M exists (it would nec be in el_from_idx[cp9->M+1] array if it does, so
   * it would be filled with 0 in above loop), or -INFTY if it doesn't exist. We don't add possibility of EL_M -> EL_M
   * self loop b/c it's impossible to do that without emitting, and we've already seen our last res emitted. 
   * either way we don't have to modify it */

  if(INBAND(i, M)) { /* if i is in k == M's band */
    assert(M == kmax[i]);
    kpcur = M-kmin[i];
    mmx[i][kpcur]  = 0. + 
      ILogsum(elmx[i][kpcur] + CP9TSC(cp9O_MEL, M),/* M_M<-EL_M<-E, with 0 self loops in EL_M */
	      CP9TSC(cp9O_ME,M));                      /* M_M<-E ... everything ends in E (the 0; 2^0=1.0) */
    mmx[i][kpcur] += cp9->msc[dsq[i]][M];  /* ... + emitted match symbol */
    imx[i][kpcur]  = 0. + CP9TSC(cp9O_IM,M);     /* I_M<-E ... everything ends in E (the 0; 2^0=1.0) */
    imx[i][kpcur] += cp9->isc[dsq[i]][M];  /* ... + emitted insert symbol */
    dmx[i][kpcur]  = CP9TSC(cp9O_DM,M);          /* D_M<-E */
    kx = M-1;
    kpcur--;
  }
  else { kx = kmax[i]; kpcur = kmax[i]-kmin[i]; }
  /*******************************************************************
   * No need to look at EL_k->M_M b/c elmx[i] with i == L means last emitted residue was L+1 
   * and this is impossible if we've come from M_M (only would be valid if we were coming from
   * E which is handled above with the EL_k->E code). 
   *******************************************************************/

  for (k = kx; k >= kmin[i]; k--, kpcur--)
    {
      mmx[i][kpcur]  = 0 + CP9TSC(cp9O_ME,k);  /*M_k<- E */

      if(INBAND(i, k+1)) { 
	mmx[i][kpcur]  = ILogsum(mmx[i][kpcur], dmx[i][kpcur+1] + CP9TSC(cp9O_MD,k));
      }
      if(cp9->flags & CPLAN9_EL)
	mmx[i][kpcur]  = ILogsum(mmx[i][kpcur], elmx[i][kpcur] + CP9TSC(cp9O_MEL,k));
      
      mmx[i][kpcur] += cp9->msc[dsq[i]][k];
      
      /*******************************************************************
       * No need to look at EL_k->M_M b/c elmx[i] with i == L means last emitted residue was L+1 
       * and this is impossible if we've come from M_M (only would be valid if we were coming from
       * E which is handled above with the EL_k->E code). 
       *******************************************************************/

      if(INBAND(i, k+1)) { 
	imx[i][kpcur]  = dmx[i][kpcur+1] + CP9TSC(cp9O_ID,k);
	imx[i][kpcur] += cp9->isc[dsq[i]][k];

	dmx[i][kpcur]  = dmx[i][kpcur+1] + CP9TSC(cp9O_DD,k);
      }
      else { 
	imx[i][kpcur] = -INFTY;
	dmx[i][kpcur] = -INFTY;
      }
      /* elmx[i][k] was set above, out of order */

      ////printf("mmx[i:%4d][k:%4d(kp:%4d)] %10d\n", i, k, kpcur, mmx[i][kpcur]);
      ////printf("imx[i:%4d][k:%4d(kp:%4d)] %10d\n", i, k, kpcur, imx[i][kpcur]);
      ////printf("dmx[i:%4d][k:%4d(kp:%4d)] %10d\n", i, k, kpcur, dmx[i][kpcur]);
    }
  
  if(INBAND(i, 0)) { 
    /* remember M_0 is special, the B state, a non-emitter */
    mmx[i][0]  = dmx[i][1] + CP9TSC(cp9O_MD,0); /* M_0(B)->D_1, no seq emitted, all deletes */
    /* above line is diff from CPBackwardOLD() which has mmx[i][0] = -INFTY; */
    imx[i][0]  = dmx[i][1] + CP9TSC(cp9O_ID,0);
    imx[i][0] += cp9->isc[dsq[i]][0];
    
    dmx[i][0]   = -INFTY; /*D_0 doesn't exist*/
    elmx[i][0]  = -INFTY; /*EL_0 doesn't exist*/
  }
     
  /*****************************************************************
   * The main loop: scan the sequence from position j0-1 to i0.
   *****************************************************************/
  /* Reision */
  for (i = L-1; i >= 1; i--) 
    {
      /* init EL mx to -INFTY */
      kpcur = 0;
      for (k = kmin[i]; k <= kmax[i]; k++, kpcur++) elmx[i][kpcur] = -INFTY;
      
      /* Deal with node k == M first */
      /* elmx[i][k] could have come from self (EL_k), we 
       * can't have come from END b/c we haven't emitted the last res of the seq yet.
       */
      if(INBAND(i, M)) { 
	kpcur = M-kmin[i];
	if((cp9->flags & CPLAN9_EL) && (cp9->has_el[M]))
	  elmx[i][kpcur] = elmx[i][kpcur] + cp9->el_selfsc;

	if(INBAND(i+1, M)) { 
	  kpprv = M-kmin[i+1];
	  mmx[i][kpcur]  = imx[i+1][kpprv] + CP9TSC(cp9O_MI,M);
	  mmx[i][kpcur] += cp9->msc[dsq[i]][M];

	  imx[i][kpcur]  = imx[i+1][kpprv] + CP9TSC(cp9O_II,M);
	  imx[i][kpcur] += cp9->isc[dsq[i]][M];
      
	  dmx[i][kpcur]  = imx[i+1][kpprv] + CP9TSC(cp9O_DI,M);
	}
	else { 
	  mmx[i][kpcur] = imx[i][kpcur] = dmx[i][kpcur] = -INFTY;
	}

	if((cp9->flags & CPLAN9_EL) && (cp9->has_el[M]))
	  mmx[i][kpcur] = ILogsum(mmx[i][kpcur], elmx[i][kpcur] + CP9TSC(cp9O_MEL,M));
	 	 
	/*******************************************************************
	 * 1b Handle EL, looking at EL_k->M_M for all valid k.
	 * EL_k->M_M transition, which has no transition penalty */
	if(INBAND(i+1, M)) { 
	  if(cp9->flags & CPLAN9_EL)
	    {
	      for(c = 0; c < cp9->el_from_ct[M]; c++) /* el_from_ct[M] holds # ELs that can go to M_M */
		if(INBAND(i, cp9->el_from_idx[M][c])) { 
		  kpcur_el = cp9->el_from_idx[M][c] - kmin[i];
		  elmx[i][kpcur_el] = ILogsum(elmx[i][kpcur_el], mmx[i+1][kpprv]);
		}
	    }
	}
      }

      /*********************************************************/
      /* MATCH: *_k <- M_k+1 transitions FROM a match state*/
      kn = ESL_MAX(kmin[i], kmin[i+1]-1); /* start at first cell from which we can look ahead to a valid cell at mmx[i-1][k-1] */
      kn = ESL_MAX(kn, 1); /* kn can't go all the way down to 0, that's a special case, handled outside the main loop */
      kx = ESL_MIN(kmax[i], kmax[i+1]-1);

      for (kpcur = 0; kpcur < ESL_MIN(kn-kmin[i], kmax[i]-kmin[i]+1); kpcur++) mmx[i][kpcur] = imx[i][kpcur] = dmx[i][kpcur] = elmx[i][kpcur] = -INFTY; /* impossible to reach these guys */
      for (kpcur = ESL_MAX(0,kx-kmin[i]+1); kpcur <= kmax[i]-kmin[i]; kpcur++) mmx[i][kpcur] = imx[i][kpcur] = dmx[i][kpcur] = elmx[i][kpcur] = -INFTY; /* impossible to reach these guys */

      kpcur = kx - kmin[i]; /* unnec, loop above ends with this */
      kpprv = kx - kmin[i+1];
      for (k = kx; k >= kn; k--, kpcur--, kpprv--)
	{
	  /* Handle EL, looking at EL_k->M_k for all valid k and EL_k->EL_k
	   * we're going backwards so we have to work out of order
	   * we could get around this by storing the nodes each EL goes TO
	   * in an el_to_ct[] vector. */
	  if(cp9->flags & CPLAN9_EL) {
	    for(c = 0; c < cp9->el_from_ct[k]; c++) { /* el_from_ct[k] holds # ELs that can go to M_k */
	      if(INBAND(i, cp9->el_from_idx[k][c])) { 
		kpcur_el = cp9->el_from_idx[k][c] - kmin[i];
		elmx[i][kpcur_el] = ILogsum(elmx[i][kpcur_el], mmx[i+1][kpprv]);
		/* EL<-M, penalty incurred when we enter EL (i.e. leave going backwards) */
	      }
	    }
	  }
	  
	  /* Finish off elmx[i][k] with possibility of coming from self (EL_k), 
	   * elmx[i][k] will have been filled by block above for ks > current k,
	   * no M_k -> EL_k' with k' > k */
	  if(INBAND(i+1, k)) {
	    if((cp9->flags & CPLAN9_EL) && (cp9->has_el[k]))
	      elmx[i][kpcur] = ILogsum(elmx[i][kpcur], elmx[i+1][kpprv] + cp9->el_selfsc);
	  }
	  mmx[i][kpcur] = mmx[i+1][kpprv+1] + CP9TSC(cp9O_MM,k);
	  imx[i][kpcur] = mmx[i+1][kpprv+1] + CP9TSC(cp9O_IM,k);
	  dmx[i][kpcur] = mmx[i+1][kpprv+1] + CP9TSC(cp9O_DM,k);
	}

      /*********************************************************/
      /* INSERTIONS: *_k <- I_k+1 transitions FROM a insert state*/
      kn = ESL_MAX(kmin[i], kmin[i+1]); /* start at first cell from which we can look ahead to a valid cell at imx[i-1][k] */
      kn = ESL_MAX(kn, 1); /* kn can't go all the way down to 0, that's a special case, handled outside the main loop */
      kx = ESL_MIN(kmax[i], kmax[i+1]);

      kpcur = kx - kmin[i]; /* unnec, loop above ends with this */
      kpprv = kx - kmin[i+1];
      for (k = kx; k >= kn; k--, kpcur--, kpprv--)
	{
	  mmx[i][kpcur] = ILogsum(mmx[i][kpcur], imx[i+1][kpprv] + CP9TSC(cp9O_MI,k));
	  imx[i][kpcur] = ILogsum(imx[i][kpcur], imx[i+1][kpprv] + CP9TSC(cp9O_II,k));
	  dmx[i][kpcur] = ILogsum(dmx[i][kpcur], imx[i+1][kpprv] + CP9TSC(cp9O_DI,k));
	}

      /*********************************************************/
      /* DELETIONS: *_k <- D_k+1 transitions FROM a delete state*/
      kn = ESL_MAX(kmin[i], kmin[i]-1); /* start at first cell from which we can look ahead to a valid cell at imx[i-1][k] */
      kn = ESL_MAX(kn, 1); /* kn can't go all the way down to 0, that's a special case, handled outside the main loop */
      kx = ESL_MIN(kmax[i], kmax[i]-1);

      kpcur = kx - kmin[i]; /* unnec, loop above ends with this */
      for (k = kx; k >= kn; k--, kpcur--)
	{
	  mmx[i][kpcur] = ILogsum(mmx[i][kpcur], dmx[i][kpcur+1] + CP9TSC(cp9O_MD,k));
	  imx[i][kpcur] = ILogsum(imx[i][kpcur], dmx[i][kpcur+1] + CP9TSC(cp9O_ID,k));
	  dmx[i][kpcur] = ILogsum(dmx[i][kpcur], dmx[i][kpcur+1] + CP9TSC(cp9O_DD,k));

	  /* now add in the emission score */
	  mmx[i][kpcur] += cp9->msc[dsq[i]][k];
	  imx[i][kpcur] += cp9->isc[dsq[i]][k];
	}
      /* there's one valid cell that we didn't consider in above loop to add emission score
       * b/c D_k <- D_k+1 (so k+1 is out of bounds for Delete) */
      for(k = kx+1; k <= kmax[i]; k++) { 
	kpcur = k - kmin[i]; /* unnec, loop above ends with this */
	mmx[i][kpcur] += cp9->msc[dsq[i]][k];
	imx[i][kpcur] += cp9->isc[dsq[i]][k];
      }
      /*********************************************************/
      kpcur = kmax[i] - kmin[i];
      for(k = kmax[i]; k >= kmin[i] && k > 0; k--, kpcur--) { 
	////printf("mmx[i:%4d][k:%4d(kp:%4d)] %10d\n", i, k, kpcur, mmx[i][kpcur]);
	////printf("imx[i:%4d][k:%4d(kp:%4d)] %10d\n", i, k, kpcur, imx[i][kpcur]);
	////printf("dmx[i:%4d][k:%4d(kp:%4d)] %10d\n", i, k, kpcur, dmx[i][kpcur]);
      }
      /* special case when k == 0 */
      kpcur = 0;
      kpprv = 0 - kmin[i+1];
      if(INBAND(i, 0)) { 
	assert(kmin[i] == 0);
	dmx[i][kpcur]  = -INFTY; /* D_0 does not exist */
	elmx[i][kpcur] = -INFTY; /* EL_0 does not exist */

	/* INSERT k=0 */
	imx[i][kpcur] = -INFTY;
	/* imx[i][0] is filled same as imx[i][1..k] in the loop above */
	if(INBAND(i+1, 1)) { 
	  if(mmx[i+1][kpprv+1] != -INFTY) 
	    imx[i][kpcur] = ILogsum(imx[i][kpcur], mmx[i+1][kpprv+1] + CP9TSC(cp9O_IM,0));
	}
	if(INBAND(i+1, 0)) { 
	  if(imx[i+1][kpprv] != -INFTY) 
	    imx[i][kpcur] = ILogsum(imx[i][kpcur], imx[i+1][kpprv] + CP9TSC(cp9O_II,0));
	}
	if(INBAND(i, 1)) { 
	  if(dmx[i][kpcur+1] != -INFTY) 
	    imx[i][kpcur] = ILogsum(imx[i][kpcur], dmx[i][kpcur+1] + CP9TSC(cp9O_ID,0));
	}

	/*M_0 is the B state, it doesn't emit, and can be reached from any match via a begin transition */
	kprvn = ESL_MAX(1, kmin[i+1]); 
	kprvx = kmax[i+1]; 
	/* careful, don't change kpcur - this is a special case the M_0 state, we loop over all M_k children, only kpprv changes */
	kpprv = kprvx - kmin[i+1];
	mmx[i][kpcur] = -INFTY;
	for(kprv = kprvx; kprv >= kprvn; kprv--, kpprv--) { 
	  if(mmx[i+1][kpprv] != -INFTY)
	    mmx[i][kpcur] = ILogsum(mmx[i][kpcur], (mmx[i+1][kpprv] + CP9TSC(cp9O_BM,kprv)));
	}
	k = 0;
	if(INBAND(i+1, 0)) { 
	  kpprv = k - kmin[i+1];
	  if(imx[i+1][kpprv] != -INFTY) { 
	    mmx[i][kpcur] = ILogsum(mmx[i][kpcur], (imx[i+1][kpprv] + CP9TSC(cp9O_MI,0)));
	  }
	}
	if(INBAND(i, 1)) { 
	  if(dmx[i][kpcur+1] != -INFTY) { 
	    mmx[i][kpcur] = ILogsum(mmx[i][kpcur], (dmx[i][kpcur+1] + CP9TSC(cp9O_MD,0)));     /* B->D_1 */
	  }
	}
      }
      ////printf("mmx[i:%4d][k:%4d(kp:%4d)] %10d\n", i, 0, kpcur, mmx[i][kpcur]);
      ////printf("imx[i:%4d][k:%4d(kp:%4d)] %10d\n", i, 0, kpcur, imx[i][kpcur]);
      ////printf("dmx[i:%4d][k:%4d(kp:%4d)] %10d\n", i, 0, kpcur, dmx[i][kpcur]);
    }

  /*******************************************************************/
  /* Special case: i == 0 */
     
  /* initialize all match, inserts, deletes and ELs to -INFTY 
   * deletes and M_0 cells MAY get changed later, inserts and ELs always stay -INFTY
   * b/c we need at least 1 residue to get to those cells */

  i = 0;
  kpcur = k - kmin[0];
  for (kpcur = 0; kpcur <= kmax[0] - kmin[0]; kpcur++) mmx[i][kpcur] = imx[i][kpcur] = dmx[i][kpcur] = elmx[i][kpcur] = -INFTY;

  /* D_M(i == 0) <- I_M(i = 1) */
  if(INBAND(i, M)) { 
    kpcur = M - kmin[i];
    kpprv = M - kmin[i+1];
    if(INBAND(i+1, M)) { 
      dmx[i][kpcur]  = imx[i+1][kpprv] + CP9TSC(cp9O_DI,M); 
    }
  }

  /* update D_k(i == 0) cells */
  /* D_k(i == 0) <- M_k+1(i == 1) */
  kn = ESL_MAX(kmin[i], kmin[i+1]-1); /* start at first cell from which we can look back to a valid cell at *mx[i-1][k-1] */
  kn = ESL_MAX(kn, 1); /* kn can't go all the way down to 0, that's a special case, handled outside the main loop */
  kx = ESL_MIN(kmax[i], kmax[i+1]-1);
  kpcur = kx - kmin[i]; 
  kpprv = kx - kmin[i+1];
  for (k = kx; k >= kn; k--, kpcur--, kpprv--)
    dmx[i][kpcur]  = mmx[i+1][kpprv+1] + CP9TSC(cp9O_DM,k);
  
  /* D_k(i == 0) <- I_k(i == 1) */
  kn = ESL_MAX(kmin[i], kmin[i+1]); /* start at first cell from which we can look ahead to a valid cell at imx[i-1][k] */
  kn = ESL_MAX(kn, 1); /* kn can't go all the way down to 0, that's a special case, handled outside the main loop */
  kx = ESL_MIN(kmax[i], kmax[i+1]);
  kpcur = kx - kmin[i]; /* unnec, loop above ends with this */
  kpprv = kx - kmin[i+1];
  for (k = kx; k >= kn; k--, kpcur--, kpprv--)
    dmx[i][kpcur] = ILogsum(dmx[i][kpcur], imx[i+1][kpprv] + CP9TSC(cp9O_DI,k));
  
  /* D_k(i == 0) <- D_k+1(i == 0) */
  kn = ESL_MAX(kmin[i], kmin[i]-1); /* start at first cell from which we can look ahead to a valid cell at imx[i-1][k] */
  kn = ESL_MAX(kn, 1); /* kn can't go all the way down to 0, that's a special case, handled outside the main loop */
  kx = ESL_MIN(kmax[i], kmax[i]-1);
  kpcur = kx - kmin[i]; /* unnec, loop above ends with this */
  for (k = kx; k >= kn; k--, kpcur--)
    dmx[i][kpcur] = ILogsum(dmx[i][kpcur], dmx[i][kpcur+1] + CP9TSC(cp9O_DD,k));

    ////printf("mmx[i:%4d][k:%4d(kp:%4d)] %10d\n", i, k, kpcur, mmx[i][kpcur]);
    ////printf("imx[i:%4d][k:%4d(kp:%4d)] %10d\n", i, k, kpcur, imx[i][kpcur]);
    ////printf("dmx[i:%4d][k:%4d(kp:%4d)] %10d\n", i, k, kpcur, dmx[i][kpcur]);

  /* Case when k == 0 (i is still 0) */
  k = 0;
  if(INBAND(i, 0)) { 
    assert(kmin[i]  == 0);
    imx[i][0] = -INFTY; /* need seq to get here */
    dmx[i][0]   = -INFTY; /* D_0 does not exist */
    elmx[i][0]  = -INFTY; /* EL_0 does not exist */
    mmx[i][0] = -INFTY;

    /*M_0 is the B state, it doesn't emit, and can be reached from any match via a begin transition */
    kprvn = ESL_MAX(1, kmin[i+1]); 
    kprvx = kmax[i+1]; 

    kpcur = 0; /* special case, M_0 = B state */
    /* careful, this is a different assignment to kpcur b/c it's the M_0 state, we loop over all M_k children, only kpprv changes */
    kpprv = kprvx - kmin[i+1];
    mmx[i][kpcur] = -INFTY;
    for(kprv = kprvx; kprv >= kprvn; kprv--, kpprv--) { 
      if(mmx[i+1][kpprv] != -INFTY)
	mmx[i][kpcur] = ILogsum(mmx[i][kpcur], (mmx[i+1][kpprv] + CP9TSC(cp9O_BM,kprv)));
    }
    k = 0;
    if(INBAND(i+1, 0)) { 
      kpprv = k - kmin[i+1];
      if(imx[i+1][kpprv] != -INFTY) { 
	mmx[i][kpcur] = ILogsum(mmx[i][kpcur], (imx[i+1][kpprv] + CP9TSC(cp9O_MI,0)));
      }
    }
    if(INBAND(i, 1)) { 
      if(dmx[i][kpcur+1] != -INFTY) { 
	mmx[i][kpcur] = ILogsum(mmx[i][kpcur], (dmx[i][kpcur+1] + CP9TSC(cp9O_MD,0)));     /* B->D_1 */
      }
    }
  }       
  /* No EL contribution here, can't go B->EL_* */
      
  /* final score (fsc) is Scorify(mmx[i][0]); */
  /**********************************************************************************/
  /* End of Backward recursion */
  
  if(ret_sc != NULL) *ret_sc = Scorify(mmx[i][0]);
  ESL_DPRINTF1(("#DEBUG: cp9_BackwardP7B() return score: %10.4f\n", Scorify(mmx[i][0])));
  return eslOK;
}


/*****************************************************************
 * Float-precision P7-banded CP9 F/B/Posterior triple.
 *
 * Mechanical mirrors of cp9_ForwardP7B, cp9_BackwardP7B, cp9_PosteriorP7B
 * above. ILogsum -> p7_FLogsum, -INFTY -> -eslINFINITY, int sentinels ->
 * float sentinels, all model scores read on the fly via Scorify().
 *
 * Used only by the truncated cmalign --p7band band derivation path
 * (cp9_IterateSeq2BandsP7B do_trunc=TRUE branch).
 *****************************************************************/

int
cp9_ForwardP7BF(CP9_t *cp9, char *errbuf, CP9_FMX *mx, ESL_DSQ *dsq, int L, int *kmin, int *kmax, float *ret_sc)
{
  int          status;
  int          i;
  int          k;
  float      **mmx;
  float      **imx;
  float      **dmx;
  float      **elmx;
  float       *erow;
  int          M;
  int          kp, kn, kx, kpcur, kpprv;

  if(cp9 == NULL)                      ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_ForwardP7BF, cp9 is NULL.\n");
  if(dsq == NULL)                      ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_ForwardP7BF, dsq is NULL.");
  if(mx == NULL)                       ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_ForwardP7BF, mx is NULL.\n");
  if(mx->M != cp9->M)                  ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_ForwardP7BF, mx->M != cp9->M.\n");
  if(kmin == NULL || kmax == NULL)     ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_ForwardP7BF, kmin/kmax NULL.\n");
  /* EL states handled below via CPLAN9_EL guard, mirroring cp9_ForwardP7B_OLD_WITH_EL. */

  M = cp9->M;
  int const *tsc = cp9->otsc;

  if((status = GrowCP9FMatrix(mx, errbuf, L, M, kmin, kmax, &mmx, &imx, &dmx, &elmx, &erow)) != eslOK) return status;

  /* Init zero row */
  mmx[0][0] = 0.;
  imx[0][0] = -eslINFINITY;
  dmx[0][0] = -eslINFINITY;
  elmx[0][0]= -eslINFINITY;
  erow[0]   = -eslINFINITY;

  float sc;
  i = 0;
  kn = ESL_MAX(1, kmin[0]);
  kp = kn - kmin[0];
  for (k = kn; k <= kmax[0]; k++, kp++) {
    assert(kp >= 0);
    mmx[0][kp] = imx[0][kp] = elmx[0][kp] = -eslINFINITY;
    sc = -eslINFINITY;
    if(kp > 0) {
      sc = p7_FLogsum(p7_FLogsum(mmx[0][kp-1] + Scorify(CP9TSC(cp9O_MD,k-1)),
				 imx[0][kp-1] + Scorify(CP9TSC(cp9O_ID,k-1))),
		      dmx[0][kp-1] + Scorify(CP9TSC(cp9O_DD,k-1)));
    }
    dmx[0][kp] = sc;
  }
  erow[0] = -eslINFINITY;
  if(INBAND(0, M)) { erow[0] = dmx[0][M] + Scorify(CP9TSC(cp9O_DM,M)); }

  for (i = 1; i <= L; i++) {
      int const *isc_i = cp9->isc[dsq[i]];
      int const *msc_i = cp9->msc[dsq[i]];
      float endsc = -eslINFINITY;
      float sc;

      if(kmin[i] == 0) {
	mmx[i][0]  = -eslINFINITY;
	dmx[i][0]  = -eslINFINITY;
	elmx[i][0] = -eslINFINITY;
	sc = p7_FLogsum(p7_FLogsum(mmx[i-1][0] + Scorify(CP9TSC(cp9O_MI,0)),
				   imx[i-1][0] + Scorify(CP9TSC(cp9O_II,0))),
			dmx[i-1][0] + Scorify(CP9TSC(cp9O_DI,0)));
	imx[i][0] = sc + Scorify(isc_i[0]);
	kn = 1;
      }
      else {
	kn = kmin[i];
      }

      /* match */
      kn = ESL_MAX(kn, (kmin[i-1]+1));
      kx = ESL_MIN(kmax[i], kmax[i-1]+1);

      for (kpcur = 0; kpcur < ESL_MIN(kn-kmin[i], kmax[i]-kmin[i]+1); kpcur++) mmx[i][kpcur] = -eslINFINITY;
      for (kpcur = ESL_MAX(0,kx-kmin[i]+1); kpcur <= kmax[i]-kmin[i]; kpcur++) mmx[i][kpcur] = -eslINFINITY;

      for (kpcur = 0; kpcur < ESL_MIN(kn-kmin[i], kmax[i]-kmin[i]+1); kpcur++) elmx[i][kpcur] = -eslINFINITY;
      for (kpcur = ESL_MAX(0,kx-kmin[i]+1); kpcur <= kmax[i]-kmin[i]; kpcur++) elmx[i][kpcur] = -eslINFINITY;

      kpcur = kn - kmin[i];
      kpprv = kn - kmin[i-1];
      for (k = kn; k <= kx; k++, kpcur++, kpprv++) {
	assert((kpprv-1) >= 0);

	sc = p7_FLogsum(p7_FLogsum(mmx[i-1][kpprv-1] + Scorify(CP9TSC(cp9O_MM,k-1)),
				   imx[i-1][kpprv-1] + Scorify(CP9TSC(cp9O_IM,k-1))),
			dmx[i-1][kpprv-1] + Scorify(CP9TSC(cp9O_DM,k-1)));

	if(INBAND(i-1, 0)) {
	  assert(kmin[(i-1)] == 0);
	  if(mmx[i-1][0] != -eslINFINITY)
	    sc = p7_FLogsum(sc, mmx[i-1][0] + Scorify(CP9TSC(cp9O_BM,k)));
	}

	if (cp9->flags & CPLAN9_EL) {
	  int c_el, kpprv_el;
	  for (c_el = 0; c_el < cp9->el_from_ct[k]; c_el++) {
	    if (INBAND(i-1, cp9->el_from_idx[k][c_el])) {
	      kpprv_el = cp9->el_from_idx[k][c_el] - kmin[i-1];
	      sc = p7_FLogsum(sc, elmx[i-1][kpprv_el]);
	    }
	  }
	}

	if(sc != -eslINFINITY) {
	  mmx[i][kpcur] = sc + Scorify(msc_i[k]);
	  endsc = p7_FLogsum(endsc, mmx[i][kpcur] + Scorify(CP9TSC(cp9O_ME,k)));
	}
	else {
	  mmx[i][kpcur] = -eslINFINITY;
	}

	{
	  float el_sc = -eslINFINITY;
	  if ((cp9->flags & CPLAN9_EL) && cp9->has_el[k]) {
	    el_sc = mmx[i][kpcur] + Scorify(CP9TSC(cp9O_MEL, k)); /* M_k -> EL_k */
	    if (INBAND(i-1, k)) {                                  /* EL self-loop */
	      int kpprv_el = k - kmin[i-1];
	      el_sc = p7_FLogsum(el_sc, elmx[i-1][kpprv_el] + Scorify(cp9->el_selfsc));
	    }
	  }
	  elmx[i][kpcur] = el_sc;
	}
      }

      /* brief 26_0430-134: BM (B->M_k) coverage fix per brief 26_0430-130 phase 3.
       * The match loop above fills only k in [kn,kx], clipped by the diagonal
       * MM-predecessor bound (kx=ESL_MIN(kmax[i],kmax[i-1]+1) at :1979).  The
       * begin transition B=(i-1,0)->M_k and the EL-from-into-M contributions do
       * not need that diagonal predecessor, but folding them inside the match
       * loop gates them on it.  When row 0 is pinned by IBV (kmax[0]=0), row 1
       * enters only M_1, dropping begin mass into M_2..M_M that cp9_BackwardP7BF
       * counts over the full child band -> F<B asymmetry (brief 26_0430-130).
       *
       * Fix: a dedicated pass over the full row band [max(1,kmin[i]),kmax[i]]
       * MINUS the [kn,kx] range already handled by the match loop (skip to avoid
       * double-count).  These supplementary cells have no in-band diagonal
       * predecessor (low range k<=kmin[i-1], high range k>kmax[i-1]+1), so sc
       * starts at -inf; we add begin + EL-from-into-M and fold the result into
       * endsc (M->E) and elmx (M->EL) exactly as the match loop does.  Runs
       * BEFORE the delete recursion and erow accumulation so they observe the
       * corrected mmx[i][k]/endsc/elmx.  Gated on INBAND(i-1,0) && mmx[i-1][0] !=
       * -inf (begin mass present).  kn/kx still hold the match-loop values here
       * (the insert loop reassigns them below).  :1979 stays UNCHANGED.
       */
      if(INBAND(i-1, 0) && mmx[i-1][0] != -eslINFINITY) {
	int k_lo = ESL_MAX(1, kmin[i]);
	int k_hi = kmax[i];
	for (k = k_lo; k <= k_hi; k++) {
	  if(k >= kn && k <= kx) continue;   /* already filled by the match loop */
	  kpcur = k - kmin[i];

	  sc = mmx[i-1][0] + Scorify(CP9TSC(cp9O_BM,k));  /* begin; diagonal predecessor out of band here */

	  if (cp9->flags & CPLAN9_EL) {
	    int c_el, kpprv_el;
	    for (c_el = 0; c_el < cp9->el_from_ct[k]; c_el++) {
	      if (INBAND(i-1, cp9->el_from_idx[k][c_el])) {
		kpprv_el = cp9->el_from_idx[k][c_el] - kmin[i-1];
		sc = p7_FLogsum(sc, elmx[i-1][kpprv_el]);
	      }
	    }
	  }

	  if(sc != -eslINFINITY) {
	    mmx[i][kpcur] = sc + Scorify(msc_i[k]);
	    endsc = p7_FLogsum(endsc, mmx[i][kpcur] + Scorify(CP9TSC(cp9O_ME,k)));
	  }
	  else {
	    mmx[i][kpcur] = -eslINFINITY;
	  }

	  {
	    float el_sc = -eslINFINITY;
	    if ((cp9->flags & CPLAN9_EL) && cp9->has_el[k]) {
	      el_sc = mmx[i][kpcur] + Scorify(CP9TSC(cp9O_MEL, k)); /* M_k -> EL_k */
	      if (INBAND(i-1, k)) {                                  /* EL self-loop */
		int kpprv_el = k - kmin[i-1];
		el_sc = p7_FLogsum(el_sc, elmx[i-1][kpprv_el] + Scorify(cp9->el_selfsc));
	      }
	    }
	    elmx[i][kpcur] = el_sc;
	  }
	}
      }

      /* insert */
      kn = ESL_MAX(kmin[i], kmin[i-1]);
      kx = ESL_MIN(kmax[i], kmax[i-1]);

      for (kpcur = 0; kpcur < ESL_MIN(kn-kmin[i], kmax[i]-kmin[i]+1); kpcur++) imx[i][kpcur] = -eslINFINITY;
      for (kpcur = ESL_MAX(0,kx-kmin[i]+1); kpcur <= kmax[i]-kmin[i]; kpcur++) imx[i][kpcur] = -eslINFINITY;

      kpcur = kn - kmin[i];
      kpprv = kn - kmin[i-1];
      for (k = kn; k <= kx; k++, kpcur++, kpprv++) {
	assert(kpprv >= 0);
	sc = p7_FLogsum(p7_FLogsum(mmx[i-1][kpprv] + Scorify(CP9TSC(cp9O_MI,k)),
				   imx[i-1][kpprv] + Scorify(CP9TSC(cp9O_II,k))),
			dmx[i-1][kpprv] + Scorify(CP9TSC(cp9O_DI,k)));
	if(sc != -eslINFINITY) imx[i][kpcur] = sc + Scorify(isc_i[k]);
	else                   imx[i][kpcur] = -eslINFINITY;
      }

      /* delete */
      kn = kmin[i]+1;

      for (kpcur = 0; kpcur < (kn - kmin[i]); kpcur++) dmx[i][kpcur] = -eslINFINITY;
      kpcur = kn - kmin[i];
      for (k = kn; k <= kmax[i]; k++, kpcur++) {
	sc = p7_FLogsum(p7_FLogsum(mmx[i][kpcur-1] + Scorify(CP9TSC(cp9O_MD,k-1)),
				   imx[i][kpcur-1] + Scorify(CP9TSC(cp9O_ID,k-1))),
			dmx[i][kpcur-1] + Scorify(CP9TSC(cp9O_DD,k-1)));
	dmx[i][kpcur] = sc;
      }

      if(INBAND(i, M)) {
	endsc = p7_FLogsum(p7_FLogsum(endsc, dmx[i][M-kmin[i]] + Scorify(CP9TSC(cp9O_DM,M))),
			   imx[i][M-kmin[i]] + Scorify(CP9TSC(cp9O_IM,M)));
      }

      erow[i] = endsc;

      if (cp9->flags & CPLAN9_EL) {
	int c_el;
	for (c_el = 0; c_el < cp9->el_from_ct[M+1]; c_el++) {
	  if (INBAND(i, cp9->el_from_idx[M+1][c_el])) {
	    int kpel = cp9->el_from_idx[M+1][c_el] - kmin[i];
	    erow[i] = p7_FLogsum(erow[i], elmx[i][kpel]);
	  }
	}
      }
  }

  *ret_sc = erow[L];
  return eslOK;
}


int
cp9_BackwardP7BF(CP9_t *cp9, char *errbuf, CP9_FMX *mx, ESL_DSQ *dsq, int L, int *kmin, int *kmax, float *ret_sc)
{
  int          status;
  int          i;
  int          k;
  float      **mmx;
  float      **imx;
  float      **dmx;
  float      **elmx;
  float       *erow;
  int          c;
  int          M;
  int          kpcur, kpprv, kpcur_el;
  int          kprv, kprvn, kprvx;
  int          kn, kx;

  if(cp9 == NULL)                      ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_BackwardP7BF, cp9 is NULL.\n");
  if(dsq == NULL)                      ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_BackwardP7BF, dsq is NULL.");
  if(mx == NULL)                       ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_BackwardP7BF, mx is NULL.\n");
  if(mx->M != cp9->M)                  ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_BackwardP7BF, mx->M != cp9->M.\n");

  M = cp9->M;
  int const *tsc = cp9->otsc;

  if((status = GrowCP9FMatrix(mx, errbuf, L, M, kmin, kmax, &mmx, &imx, &dmx, &elmx, &erow)) != eslOK) return status;

  i = L;
  kpcur = 0;
  for (k = kmin[i]; k <= kmax[i]; k++, kpcur++) elmx[i][kpcur] = -eslINFINITY;
  if(cp9->flags & CPLAN9_EL) {
    for(c = 0; c < cp9->el_from_ct[cp9->M+1]; c++)
      if(INBAND(i, cp9->el_from_idx[M+1][c])) {
	kpcur_el = cp9->el_from_idx[M+1][c] - kmin[i];
	elmx[i][kpcur_el] = 0.;
      }
  }

  if(INBAND(i, M)) {
    assert(M == kmax[i]);
    kpcur = M-kmin[i];
    mmx[i][kpcur]  = 0. +
      p7_FLogsum(elmx[i][kpcur] + Scorify(CP9TSC(cp9O_MEL, M)),
		 Scorify(CP9TSC(cp9O_ME,M)));
    mmx[i][kpcur] += Scorify(cp9->msc[dsq[i]][M]);
    imx[i][kpcur]  = 0. + Scorify(CP9TSC(cp9O_IM,M));
    imx[i][kpcur] += Scorify(cp9->isc[dsq[i]][M]);
    dmx[i][kpcur]  = Scorify(CP9TSC(cp9O_DM,M));
    kx = M-1;
    kpcur--;
  }
  else { kx = kmax[i]; kpcur = kmax[i]-kmin[i]; }

  for (k = kx; k >= kmin[i]; k--, kpcur--) {
      mmx[i][kpcur]  = 0 + Scorify(CP9TSC(cp9O_ME,k));
      if(INBAND(i, k+1)) {
	mmx[i][kpcur]  = p7_FLogsum(mmx[i][kpcur], dmx[i][kpcur+1] + Scorify(CP9TSC(cp9O_MD,k)));
      }
      if(cp9->flags & CPLAN9_EL)
	mmx[i][kpcur]  = p7_FLogsum(mmx[i][kpcur], elmx[i][kpcur] + Scorify(CP9TSC(cp9O_MEL,k)));
      mmx[i][kpcur] += Scorify(cp9->msc[dsq[i]][k]);

      if(INBAND(i, k+1)) {
	imx[i][kpcur]  = dmx[i][kpcur+1] + Scorify(CP9TSC(cp9O_ID,k));
	imx[i][kpcur] += Scorify(cp9->isc[dsq[i]][k]);

	dmx[i][kpcur]  = dmx[i][kpcur+1] + Scorify(CP9TSC(cp9O_DD,k));
      }
      else {
	imx[i][kpcur] = -eslINFINITY;
	dmx[i][kpcur] = -eslINFINITY;
      }
  }

  if(INBAND(i, 0)) {
    mmx[i][0]  = dmx[i][1] + Scorify(CP9TSC(cp9O_MD,0));
    imx[i][0]  = dmx[i][1] + Scorify(CP9TSC(cp9O_ID,0));
    imx[i][0] += Scorify(cp9->isc[dsq[i]][0]);

    dmx[i][0]   = -eslINFINITY;
    elmx[i][0]  = -eslINFINITY;
  }

  for (i = L-1; i >= 1; i--) {
      kpcur = 0;
      for (k = kmin[i]; k <= kmax[i]; k++, kpcur++) elmx[i][kpcur] = -eslINFINITY;

      if(INBAND(i, M)) {
	kpcur = M-kmin[i];
	if((cp9->flags & CPLAN9_EL) && (cp9->has_el[M]))
	  elmx[i][kpcur] = elmx[i][kpcur] + Scorify(cp9->el_selfsc);

	if(INBAND(i+1, M)) {
	  kpprv = M-kmin[i+1];
	  mmx[i][kpcur]  = imx[i+1][kpprv] + Scorify(CP9TSC(cp9O_MI,M));
	  mmx[i][kpcur] += Scorify(cp9->msc[dsq[i]][M]);

	  imx[i][kpcur]  = imx[i+1][kpprv] + Scorify(CP9TSC(cp9O_II,M));
	  imx[i][kpcur] += Scorify(cp9->isc[dsq[i]][M]);

	  dmx[i][kpcur]  = imx[i+1][kpprv] + Scorify(CP9TSC(cp9O_DI,M));
	}
	else {
	  mmx[i][kpcur] = imx[i][kpcur] = dmx[i][kpcur] = -eslINFINITY;
	}

	if((cp9->flags & CPLAN9_EL) && (cp9->has_el[M]))
	  mmx[i][kpcur] = p7_FLogsum(mmx[i][kpcur], elmx[i][kpcur] + Scorify(CP9TSC(cp9O_MEL,M)));

	if(INBAND(i+1, M)) {
	  if(cp9->flags & CPLAN9_EL) {
	    for(c = 0; c < cp9->el_from_ct[M]; c++)
	      if(INBAND(i, cp9->el_from_idx[M][c])) {
		kpcur_el = cp9->el_from_idx[M][c] - kmin[i];
		elmx[i][kpcur_el] = p7_FLogsum(elmx[i][kpcur_el], mmx[i+1][kpprv]);
	      }
	  }
	}
      }

      /* MATCH transitions */
      kn = ESL_MAX(kmin[i], kmin[i+1]-1);
      kn = ESL_MAX(kn, 1);
      kx = ESL_MIN(kmax[i], kmax[i+1]-1);

      for (kpcur = 0; kpcur < ESL_MIN(kn-kmin[i], kmax[i]-kmin[i]+1); kpcur++) mmx[i][kpcur] = imx[i][kpcur] = dmx[i][kpcur] = elmx[i][kpcur] = -eslINFINITY;
      for (kpcur = ESL_MAX(0,kx-kmin[i]+1); kpcur <= kmax[i]-kmin[i]; kpcur++) mmx[i][kpcur] = imx[i][kpcur] = dmx[i][kpcur] = elmx[i][kpcur] = -eslINFINITY;

      kpcur = kx - kmin[i];
      kpprv = kx - kmin[i+1];
      for (k = kx; k >= kn; k--, kpcur--, kpprv--) {
	  if(cp9->flags & CPLAN9_EL) {
	    for(c = 0; c < cp9->el_from_ct[k]; c++) {
	      if(INBAND(i, cp9->el_from_idx[k][c])) {
		kpcur_el = cp9->el_from_idx[k][c] - kmin[i];
		elmx[i][kpcur_el] = p7_FLogsum(elmx[i][kpcur_el], mmx[i+1][kpprv]);
	      }
	    }
	  }

	  if(INBAND(i+1, k)) {
	    if((cp9->flags & CPLAN9_EL) && (cp9->has_el[k]))
	      elmx[i][kpcur] = p7_FLogsum(elmx[i][kpcur], elmx[i+1][kpprv] + Scorify(cp9->el_selfsc));
	  }
	  mmx[i][kpcur] = mmx[i+1][kpprv+1] + Scorify(CP9TSC(cp9O_MM,k));
	  imx[i][kpcur] = mmx[i+1][kpprv+1] + Scorify(CP9TSC(cp9O_IM,k));
	  dmx[i][kpcur] = mmx[i+1][kpprv+1] + Scorify(CP9TSC(cp9O_DM,k));
      }

      /* INSERT transitions */
      kn = ESL_MAX(kmin[i], kmin[i+1]);
      kn = ESL_MAX(kn, 1);
      kx = ESL_MIN(kmax[i], kmax[i+1]);

      kpcur = kx - kmin[i];
      kpprv = kx - kmin[i+1];
      for (k = kx; k >= kn; k--, kpcur--, kpprv--) {
	  mmx[i][kpcur] = p7_FLogsum(mmx[i][kpcur], imx[i+1][kpprv] + Scorify(CP9TSC(cp9O_MI,k)));
	  imx[i][kpcur] = p7_FLogsum(imx[i][kpcur], imx[i+1][kpprv] + Scorify(CP9TSC(cp9O_II,k)));
	  dmx[i][kpcur] = p7_FLogsum(dmx[i][kpcur], imx[i+1][kpprv] + Scorify(CP9TSC(cp9O_DI,k)));
      }

      /* DELETE transitions */
      kn = ESL_MAX(kmin[i], kmin[i]-1);
      kn = ESL_MAX(kn, 1);
      kx = ESL_MIN(kmax[i], kmax[i]-1);

      kpcur = kx - kmin[i];
      for (k = kx; k >= kn; k--, kpcur--) {
	  mmx[i][kpcur] = p7_FLogsum(mmx[i][kpcur], dmx[i][kpcur+1] + Scorify(CP9TSC(cp9O_MD,k)));
	  imx[i][kpcur] = p7_FLogsum(imx[i][kpcur], dmx[i][kpcur+1] + Scorify(CP9TSC(cp9O_ID,k)));
	  dmx[i][kpcur] = p7_FLogsum(dmx[i][kpcur], dmx[i][kpcur+1] + Scorify(CP9TSC(cp9O_DD,k)));

	  mmx[i][kpcur] += Scorify(cp9->msc[dsq[i]][k]);
	  imx[i][kpcur] += Scorify(cp9->isc[dsq[i]][k]);
      }
      for(k = kx+1; k <= kmax[i]; k++) {
	kpcur = k - kmin[i];
	mmx[i][kpcur] += Scorify(cp9->msc[dsq[i]][k]);
	imx[i][kpcur] += Scorify(cp9->isc[dsq[i]][k]);
      }

      /* k == 0 special case */
      kpcur = 0;
      kpprv = 0 - kmin[i+1];
      if(INBAND(i, 0)) {
	assert(kmin[i] == 0);
	dmx[i][kpcur]  = -eslINFINITY;
	elmx[i][kpcur] = -eslINFINITY;

	imx[i][kpcur] = -eslINFINITY;
	if(INBAND(i+1, 1)) {
	  if(mmx[i+1][kpprv+1] != -eslINFINITY)
	    imx[i][kpcur] = p7_FLogsum(imx[i][kpcur], mmx[i+1][kpprv+1] + Scorify(CP9TSC(cp9O_IM,0)));
	}
	if(INBAND(i+1, 0)) {
	  if(imx[i+1][kpprv] != -eslINFINITY)
	    imx[i][kpcur] = p7_FLogsum(imx[i][kpcur], imx[i+1][kpprv] + Scorify(CP9TSC(cp9O_II,0)));
	}
	if(INBAND(i, 1)) {
	  if(dmx[i][kpcur+1] != -eslINFINITY)
	    imx[i][kpcur] = p7_FLogsum(imx[i][kpcur], dmx[i][kpcur+1] + Scorify(CP9TSC(cp9O_ID,0)));
	}

	kprvn = ESL_MAX(1, kmin[i+1]);
	kprvx = kmax[i+1];
	kpprv = kprvx - kmin[i+1];
	mmx[i][kpcur] = -eslINFINITY;
	for(kprv = kprvx; kprv >= kprvn; kprv--, kpprv--) {
	  if(mmx[i+1][kpprv] != -eslINFINITY)
	    mmx[i][kpcur] = p7_FLogsum(mmx[i][kpcur], (mmx[i+1][kpprv] + Scorify(CP9TSC(cp9O_BM,kprv))));
	}
	k = 0;
	if(INBAND(i+1, 0)) {
	  kpprv = k - kmin[i+1];
	  if(imx[i+1][kpprv] != -eslINFINITY) {
	    mmx[i][kpcur] = p7_FLogsum(mmx[i][kpcur], (imx[i+1][kpprv] + Scorify(CP9TSC(cp9O_MI,0))));
	  }
	}
	if(INBAND(i, 1)) {
	  if(dmx[i][kpcur+1] != -eslINFINITY) {
	    mmx[i][kpcur] = p7_FLogsum(mmx[i][kpcur], (dmx[i][kpcur+1] + Scorify(CP9TSC(cp9O_MD,0))));
	  }
	}
      }
  }

  /* Special case: i == 0 */
  i = 0;
  kpcur = k - kmin[0];
  for (kpcur = 0; kpcur <= kmax[0] - kmin[0]; kpcur++) mmx[i][kpcur] = imx[i][kpcur] = dmx[i][kpcur] = elmx[i][kpcur] = -eslINFINITY;

  if(INBAND(i, M)) {
    kpcur = M - kmin[i];
    kpprv = M - kmin[i+1];
    if(INBAND(i+1, M)) {
      dmx[i][kpcur]  = imx[i+1][kpprv] + Scorify(CP9TSC(cp9O_DI,M));
    }
  }

  kn = ESL_MAX(kmin[i], kmin[i+1]-1);
  kn = ESL_MAX(kn, 1);
  kx = ESL_MIN(kmax[i], kmax[i+1]-1);
  kpcur = kx - kmin[i];
  kpprv = kx - kmin[i+1];
  for (k = kx; k >= kn; k--, kpcur--, kpprv--)
    dmx[i][kpcur]  = mmx[i+1][kpprv+1] + Scorify(CP9TSC(cp9O_DM,k));

  kn = ESL_MAX(kmin[i], kmin[i+1]);
  kn = ESL_MAX(kn, 1);
  kx = ESL_MIN(kmax[i], kmax[i+1]);
  kpcur = kx - kmin[i];
  kpprv = kx - kmin[i+1];
  for (k = kx; k >= kn; k--, kpcur--, kpprv--)
    dmx[i][kpcur] = p7_FLogsum(dmx[i][kpcur], imx[i+1][kpprv] + Scorify(CP9TSC(cp9O_DI,k)));

  kn = ESL_MAX(kmin[i], kmin[i]-1);
  kn = ESL_MAX(kn, 1);
  kx = ESL_MIN(kmax[i], kmax[i]-1);
  kpcur = kx - kmin[i];
  for (k = kx; k >= kn; k--, kpcur--)
    dmx[i][kpcur] = p7_FLogsum(dmx[i][kpcur], dmx[i][kpcur+1] + Scorify(CP9TSC(cp9O_DD,k)));

  k = 0;
  if(INBAND(i, 0)) {
    assert(kmin[i] == 0);
    imx[i][0] = -eslINFINITY;
    dmx[i][0]   = -eslINFINITY;
    elmx[i][0]  = -eslINFINITY;
    mmx[i][0] = -eslINFINITY;

    kprvn = ESL_MAX(1, kmin[i+1]);
    kprvx = kmax[i+1];

    kpcur = 0;
    kpprv = kprvx - kmin[i+1];
    mmx[i][kpcur] = -eslINFINITY;
    for(kprv = kprvx; kprv >= kprvn; kprv--, kpprv--) {
      if(mmx[i+1][kpprv] != -eslINFINITY)
	mmx[i][kpcur] = p7_FLogsum(mmx[i][kpcur], (mmx[i+1][kpprv] + Scorify(CP9TSC(cp9O_BM,kprv))));
    }
    k = 0;
    if(INBAND(i+1, 0)) {
      kpprv = k - kmin[i+1];
      if(imx[i+1][kpprv] != -eslINFINITY) {
	mmx[i][kpcur] = p7_FLogsum(mmx[i][kpcur], (imx[i+1][kpprv] + Scorify(CP9TSC(cp9O_MI,0))));
      }
    }
    if(INBAND(i, 1)) {
      if(dmx[i][kpcur+1] != -eslINFINITY) {
	mmx[i][kpcur] = p7_FLogsum(mmx[i][kpcur], (dmx[i][kpcur+1] + Scorify(CP9TSC(cp9O_MD,0))));
      }
    }
  }

  if(ret_sc != NULL) *ret_sc = mmx[i][0];
  return eslOK;
}


int
cp9_PosteriorP7BF(ESL_DSQ *dsq, char *errbuf, int L, CP9_t *hmm, CP9_FMX *fmx, CP9_FMX *bmx, CP9_FMX *pmx, int *kmin, int *kmax)
{
  int i;
  int k;
  float sc;
  int M = hmm->M;
  int kp, kn, kx;

  if(bmx != pmx) GrowCP9FMatrix(pmx, errbuf, L, M, kmin, kmax, NULL, NULL, NULL, NULL, NULL);

  sc = bmx->mmx[0][0];

  assert(kmin[0] == 0);
  pmx->mmx[0][0] = fmx->mmx[0][0] + bmx->mmx[0][0] - sc;
  pmx->imx[0][0] = -eslINFINITY;
  pmx->dmx[0][0] = -eslINFINITY;
  i  = 0;
  kn = ESL_MAX(kmin[i], 1);
  kx = kmax[i];
  kp = kn - kmin[i];
  for (k = kn; k <= kx; k++, kp++) {
    pmx->mmx[0][kp] = -eslINFINITY;
    pmx->imx[0][kp] = -eslINFINITY;
    pmx->dmx[0][kp] = fmx->dmx[0][kp] + bmx->dmx[0][kp] - sc;
  }

  for (i = 1; i <= L; i++) {
    k = 0;
    if(INBAND(i,0)) {
      kp = k - kmin[i];
      assert(kp == 0);
      pmx->mmx[i][kp] = ESL_MAX(fmx->mmx[i][kp] + bmx->mmx[i][kp] - sc, -eslINFINITY);
      pmx->imx[i][kp] = ESL_MAX(fmx->imx[i][kp] + bmx->imx[i][kp] - Scorify(hmm->isc[dsq[i]][0]) - sc, -eslINFINITY);
      pmx->dmx[i][kp] = -eslINFINITY;
    }

    kn = ESL_MAX(kmin[i], 1);
    kx = kmax[i];
    kp = kn - kmin[i];
    for(k = kn; k <= kx; k++, kp++) {
      pmx->mmx[i][kp] = ESL_MAX(fmx->mmx[i][kp] + bmx->mmx[i][kp] - Scorify(hmm->msc[dsq[i]][k]) - sc, -eslINFINITY);
      pmx->imx[i][kp] = ESL_MAX(fmx->imx[i][kp] + bmx->imx[i][kp] - Scorify(hmm->isc[dsq[i]][k]) - sc, -eslINFINITY);
      pmx->dmx[i][kp] = ESL_MAX(fmx->dmx[i][kp] + bmx->dmx[i][kp] - sc, -eslINFINITY);
    }
  }
  return eslOK;
}


/* Function: cp9_CheckFBP7B()
 * 
 * Purpose:  Debugging function to make sure the P7 banded
 *           DP functions cp9_ForwardP7B() and cp9_BackwardP7B()
 *            are working by checking:
 *           For all positions i, and states k within kmin[i]..kmax[i]:
 *             sum_k f[i][k] * b[i][k] = P(x|hmm)
 *           
 * Args:     fmx    - p7 banded forward dp matrix, already filled
 *           bmx    - p7 banded backward dp matrix, already filled
 *           hmm    - the model
 *           sc     - P(x|hmm, p7 bands) the probability of the entire
 *                    seq given the model
 *           i0     - start of target subsequence (often 1, beginning of dsq)
 *           j0     - end of target subsequence (often L, end of dsq)
 *           dsq    - the digitized sequence
 *           
 * Note about sequence position indexing: although this function
 * works on a subsequence from i0 to j0, fmx and bmx have offset indices,
 * from 1 to L, with L = j0-i0+1.
 * 
 * Return:   eslOK on success;
 *           eslFAIL if any residue fails check
 */
int
cp9_CheckFBP7B(CP9_MX *fmx, CP9_MX *bmx, CP9_t *hmm, char *errbuf, float sc, int i0, int j0, ESL_DSQ *dsq, int *kmin, int *kmax)
{
  if(fmx == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_CheckFBP7B(), fmx is NULL.\n");
  if(bmx == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_CheckFBP7B(), bmx is NULL.\n");
  if(dsq == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_CheckFBP7B(), dsq is NULL.");

  int k, i;
  float max_diff;  /* maximum allowed difference between sc and 
		    * sum_k f[i][k] * b[i][k] for any i */
  float diff;
  int fb_sum;
  float fb_sc;
  int   L;		/* subsequence length */
  int   kp;             /* k', relative k within band; k-kmin[i] */
  int to_add;

  L  = j0-i0+1;		/* the length of the subsequence */
  max_diff = 0.1;       /* tolerance, must be within .1 bits of original score */

  /* In all possible paths through the model, each residue of the sequence must have 
   * been emitted by exactly 1 insert, match or EL state. */
  for (i = 1; i <= L; i++) {
    fb_sum = -INFTY;
    kp = 0;
    for (k = kmin[i]; k <= kmax[i]; k++, kp++) {
      if     (fmx->mmx[i][kp] == -INFTY) to_add = -INFTY;
      else if(bmx->mmx[i][kp] == -INFTY) to_add = -INFTY;
      else {
	to_add = fmx->mmx[i][kp] + bmx->mmx[i][kp];
	if(k > 0) to_add -= hmm->msc[dsq[i]][k];
      }
      /* hmm->msc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx
       * unless, we're talking about M_0, the B state, it doesn't emit */
      fb_sum = ILogsum(fb_sum, to_add);
      
      /*printf("fmx->mmx[i:%4d][k:%4d(%4d)]: %d\n", i, k, kp, fmx->mmx[i][kp]);
	printf("bmx->mmx[i:%4d][k:%4d(%4d)]: %d sum: %d\n", i, k, kp, (bmx->mmx[i][kp]-hmm->msc[dsq[i]][k]), fb_sum);*/
      
      if     (fmx->imx[i][kp] == -INFTY) to_add = -INFTY;
      else if(bmx->imx[i][kp] == -INFTY) to_add = -INFTY;
      else  {
	to_add  = fmx->imx[i][kp] + bmx->imx[i][kp]; 
	to_add -= hmm->isc[dsq[i]][k];
      }
      /*hmm->isc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
      fb_sum = ILogsum(fb_sum, to_add);

      /*printf("fmx->imx[i:%4d][k:%4d(%4d)]: %d\n", i, k, kp, fmx->imx[i][kp]);
	printf("bmx->imx[i:%4d][k:%4d(%4d)]: %d sum: %d\n", i, k, kp, (bmx->imx[i][kp]-hmm->isc[dsq[i]][k]), fb_sum);*/

      if     (fmx->elmx[i][kp] == -INFTY) to_add = -INFTY;
      else if(bmx->elmx[i][kp] == -INFTY) to_add = -INFTY;
      else  {
	to_add  = fmx->elmx[i][kp] + bmx->elmx[i][kp]; 
	/* EL emissions are by definition zero scoring */
      }
      fb_sum = ILogsum(fb_sum, to_add);
      
      /*printf("fmx->elmx[i:%4d][k:%4d(%4d)]: %d\n", i, k, kp, fmx->elmx[i][kp]);
	printf("bmx->elmx[i:%4d][k:%4d(%4d)]: %d sum: %d\n", i, k, kp, bmx->elmx[i][kp], fb_sum);*/
    }
    fb_sc  = Scorify(fb_sum);
    diff = fabs(fb_sc - sc);
    /*printf("FB CHECK: i: %4d %10.4f %10.4f (%10.4f)\n", i, fb_sc, sc, diff);*/
    if((fabs(diff) > max_diff)) 
      ESL_FAIL(eslFAIL, errbuf, "cp9_CheckFB(), residue at posn i:%d violates sum_k f[i][k]*b[i][k]=P(x|hmm), sum_k = %.4f bits (should be %.4f)\n", i, fb_sc, sc);
  }
  ESL_DPRINTF1(("#DEBUG: cp9_CheckFB() passed, Forward/Backward matrices pass check.\n"));
  /*printf("cp9_CheckFB() passed, Forward/Backward matrices pass check.\n");*/
  return eslOK;
}



/* Function: cp9_Seq2BandsP7B
 * Date    : EPN, Fri Aug 15 13:43:21 2008
 *
 * Purpose:  Given a CM with precalc'ed CP9 HMM, CP9Map, and HMMER3 plan 7
 *           HMM bands for the CP9 HMM DP matrices, a sequence and
 *           a CP9Bands_t structure, calculate the CP9 HMM bands and store them
 *           in the CP9Bands_t structure.
 *
 * Args:     cm          - the covariance model
 *           errbuf      - char buffer for reporting errors
 *           fmx         - CP9 dp matrix for Forward()
 *           bmx         - CP9 dp matrix for Backward()
 *           pmx         - CP9 dp matrix to fill with posteriors, can == bmx
 *           dsq         - sequence in digitized form (1..L, offset so dsq[1] is first residue)
 *           L           - length of sequence we're aligning (1..L)
 *           cp9b        - PRE-ALLOCATED, the HMM bands for this sequence, filled here.
 *           kmin        - P7 dervied band to enforce: [0.i..L] = k, min node k for residue i
 *           kmax        - P7 derived band to enforce: [0.i..L] = k, min node k for residue i
 *           i0          - first position in original sequence coords (for cp9_HMM2ijBands)
 *           j0          - final position in original sequence coords (for cp9_HMM2ijBands)
 *           pass_idx    - pipeline pass index, determines truncation mode
 *           debug_level - verbosity level for debugging printf()s
 * Return:  eslOK on success;
 *
 */
int
cp9_Seq2BandsP7B(CM_t *cm, char *errbuf, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *pmx, ESL_DSQ *dsq, int L, CP9Bands_t *cp9b, int *kmin, int *kmax, int i0, int j0, int pass_idx, int debug_level, int do_pnmono, int do_pnmono_print)
{
  int   status;
  float sc;
  CP9_t *cp9 = NULL;

  /* Contract checks */
  if(cm->cp9map == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_Seq2BandsP7B, but cm->cp9map is NULL.\n");
  if(dsq == NULL)        ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_Seq2BandsP7B, dsq is NULL.");
  if(!((cm->align_opts & CM_ALIGN_HBANDED) || (cm->search_opts & CM_SEARCH_HBANDED)))        ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_Seq2BandsP7B, CM_ALIGN_HBANDED and CM_SEARCH_HBANDED flags both down, exactly 1 must be up.\n");
  if((cm->search_opts & CM_SEARCH_HMMALNBANDS) && (!(cm->search_opts & CM_SEARCH_HBANDED))) ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_Seq2BandsP7B, CM_SEARCH_HMMALNBANDS flag raised, but not CM_SEARCH_HBANDED flag, this doesn't make sense\n");
  if(cm->tau > 0.5)      ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_Seq2BandsP7B, cm->tau (%f) > 0.5, we can't deal.", cm->tau);

  switch(pass_idx) {
    case PLI_PASS_5P_ONLY_FORCE:                cp9 = cm->Rcp9; break;
    case PLI_PASS_3P_ONLY_FORCE:                cp9 = cm->Lcp9; break;
    case PLI_PASS_5P_AND_3P_FORCE:
    case PLI_PASS_5P_AND_3P_ANY:                cp9 = cm->Tcp9; break;
    default:                                    cp9 = cm->cp9;  break;
  }
  if(cp9 == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_Seq2BandsP7B, cp9 is NULL (pass_idx %d).\n", pass_idx);

  /* Phase 1: P7-banded CP9 Forward + Backward (tau-independent) */
  if((status = cp9_ForwardP7B_OLD_WITH_EL(cp9, errbuf, fmx, dsq, L, kmin, kmax, &sc)) != eslOK) return status;
  if((status = cp9_BackwardP7B(cp9, errbuf, bmx, dsq, L, kmin, kmax, NULL)) != eslOK) return status;

  if(cm->align_opts & CM_ALIGN_CHECKFB) {
    if((status = cp9_CheckFBP7B(fmx, bmx, cp9, errbuf, sc, 1, L, dsq, kmin, kmax)) != eslOK) return status;
    printf("Forward/Backward matrices checked.\n");
  }

  /* Phase 2: F/B -> HMM bands -> CM bands (tau-dependent) */
  if((status = cp9_FBMatrices2BandsP7B(cm, errbuf, cp9, fmx, bmx, pmx, dsq, cp9b, kmin, kmax,
				       L, i0, j0, pass_idx, debug_level, do_pnmono, do_pnmono_print)) != eslOK) return status;
  return eslOK;
}


/* Function: cp9_Seq2BandsP7BF()
 *
 * Float-precision mirror of cp9_Seq2BandsP7B(). Allocates CP9_FMX matrices,
 * runs P7-banded float CP9 Forward + Backward, and dispatches to
 * cp9_FBMatrices2BandsF for posterior + band derivation.
 *
 * Currently called only by cp9_IterateSeq2BandsP7B's do_trunc=TRUE branch.
 * Don't add new callers without auditing — the float path is intentionally
 * scoped to truncated cmalign --p7band only.
 */
int
cp9_Seq2BandsP7BF(CM_t *cm, char *errbuf, CP9_FMX *fmx, CP9_FMX *bmx, CP9_FMX *pmx, ESL_DSQ *dsq, int L, CP9Bands_t *cp9b, int *kmin, int *kmax, int i0, int j0, int pass_idx, int debug_level, int do_pnmono, int do_pnmono_print)
{
  int   status;
  float sc;
  CP9_t *cp9 = NULL;

  if(cm->cp9map == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_Seq2BandsP7BF, but cm->cp9map is NULL.\n");
  if(dsq == NULL)        ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_Seq2BandsP7BF, dsq is NULL.");
  if(!((cm->align_opts & CM_ALIGN_HBANDED) || (cm->search_opts & CM_SEARCH_HBANDED)))        ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_Seq2BandsP7BF, CM_ALIGN_HBANDED and CM_SEARCH_HBANDED flags both down.\n");
  if(cm->tau > 0.5)      ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_Seq2BandsP7BF, cm->tau (%f) > 0.5.", cm->tau);

  switch(pass_idx) {
    case PLI_PASS_5P_ONLY_FORCE:                cp9 = cm->Rcp9; break;
    case PLI_PASS_3P_ONLY_FORCE:                cp9 = cm->Lcp9; break;
    case PLI_PASS_5P_AND_3P_FORCE:
    case PLI_PASS_5P_AND_3P_ANY:                cp9 = cm->Tcp9; break;
    default:                                    cp9 = cm->cp9;  break;
  }
  if(cp9 == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_Seq2BandsP7BF, cp9 is NULL (pass_idx %d).\n", pass_idx);

  if((status = cp9_ForwardP7BF (cp9, errbuf, fmx, dsq, L, kmin, kmax, &sc))   != eslOK) return status;
  if((status = cp9_BackwardP7BF(cp9, errbuf, bmx, dsq, L, kmin, kmax, NULL))  != eslOK) return status;

  if((status = cp9_FBMatrices2BandsF(cm, errbuf, cp9, fmx, bmx, pmx, dsq, cp9b, kmin, kmax,
				     L, i0, j0, pass_idx, debug_level, do_pnmono, do_pnmono_print)) != eslOK) return status;
  return eslOK;
}


/* Function: cp9_FBMatrices2BandsP7B()
 * Date:     EPN, 2026-04-30
 *
 * Purpose:  Phase 2 of p7-banded CP9 band derivation. Given filled P7-banded
 *           CP9 Forward (fmx) and Backward (bmx) matrices, derive HMM bands
 *           using cm->tau, then convert to CM bands. This is the tau-dependent
 *           half that can be iterated by cp9_IterateSeq2BandsP7B().
 *
 * Args:     cm          - the CM
 *           errbuf      - for error messages
 *           cp9         - the CP9 HMM (cm->cp9, Lcp9, Rcp9, or Tcp9 depending on pass_idx)
 *           fmx         - filled P7-banded CP9 Forward matrix (read-only)
 *           bmx         - filled P7-banded CP9 Backward matrix (read-only)
 *           pmx         - CP9 matrix for posteriors (filled here; must NOT alias bmx if iterating)
 *           dsq         - digitized sequence (1..L)
 *           cp9b        - CP9 bands structure (filled here)
 *           kmin, kmax  - P7-derived per-residue node bands
 *           L           - sequence length
 *           i0, j0      - subsequence bounds in original coords
 *           pass_idx    - pipeline pass index
 *           debug_level - verbosity
 *           do_pnmono, do_pnmono_print - pnmono flags
 *
 * Returns:  eslOK on success.
 */
int
cp9_FBMatrices2BandsP7B(CM_t *cm, char *errbuf, CP9_t *cp9, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *pmx,
			ESL_DSQ *dsq, CP9Bands_t *cp9b, int *kmin, int *kmax,
			int L, int i0, int j0, int pass_idx, int debug_level,
			int do_pnmono, int do_pnmono_print)
{
  int status;
  int use_sums      = ((cm->align_opts & CM_ALIGN_SUMS) || (cm->search_opts & CM_SEARCH_SUMS)) ? TRUE : FALSE;
  int do_old_hmm2ij = ((cm->align_opts & CM_ALIGN_HMM2IJOLD) || (cm->search_opts & CM_SEARCH_HMM2IJOLD)) ? TRUE : FALSE;
  int do_trunc      = cm_pli_PassAllowsTruncation(pass_idx);

  /* Step 2: F/B -> HMM bands. */
  if(use_sums) {
    printf("USE SUMS!\n");
    exit(1);
  }
  else {
    if((status = cp9_FB2HMMBandsP7B(cp9, errbuf, dsq, fmx, bmx, pmx, cp9b, L, cp9b->hmm_M,
				    (1.-cm->tau), do_old_hmm2ij, kmin, kmax, debug_level,
				    do_pnmono, do_pnmono_print)) != eslOK) return status;
    cp9b->tau = cm->tau;
  }
  if(debug_level > 0) cp9_DebugPrintHMMBands(stdout, L, cp9b, cm->tau, 1);

  /* Step 2b: Shift HMM bands from 1..L to i0..j0 coordinate system. */
  if(i0 != 1) {
    int offset = i0 - 1;
    int k;
    for(k = 0; k <= cp9b->hmm_M; k++) {
      if(cp9b->pn_min_m[k] != -1) { cp9b->pn_min_m[k] += offset; cp9b->pn_max_m[k] += offset; }
      if(cp9b->pn_min_i[k] != -1) { cp9b->pn_min_i[k] += offset; cp9b->pn_max_i[k] += offset; }
      if(cp9b->pn_min_d[k] != -1) { cp9b->pn_min_d[k] += offset; cp9b->pn_max_d[k] += offset; }
    }
  }

  /* Step 2c: Set truncation candidate valid arrays. */
  if(do_trunc) {
    cp9_PredictStartAndEndPositionsP7B(pmx, cp9b, kmin, kmax, i0, j0);
    if((status = cp9_MarginalCandidatesFromStartEndPositions(cm, cp9b, pass_idx, errbuf)) != eslOK) return status;
  }
  else {
    esl_vec_ISet(cp9b->Jvalid, cm->M+1, TRUE);
    esl_vec_ISet(cp9b->Lvalid, cm->M+1, FALSE);
    esl_vec_ISet(cp9b->Rvalid, cm->M+1, FALSE);
    esl_vec_ISet(cp9b->Tvalid, cm->M+1, FALSE);
  }

  /* Step 3: HMM bands -> CM bands. */
  if(do_old_hmm2ij) {
    if((status = cp9_HMM2ijBands_OLD(cm, errbuf, cm->cp9b, cm->cp9map, i0, j0, TRUE, debug_level)) != eslOK) return status;
  }
  else {
    if((status = cp9_HMM2ijBands(cm, errbuf, cp9, cm->cp9b, cm->cp9map, i0, j0, TRUE, do_trunc, debug_level)) != eslOK) return status;
  }
  if((status = cp9_GrowHDBands(cp9b, errbuf)) != eslOK) return status;
  ij2d_bands(cm, cp9b, do_trunc, debug_level);

#if eslDEBUGLEVEL >= 1
  if((status = cp9_ValidateBands(cm, errbuf, cp9b, i0, j0, do_trunc)) != eslOK) return status;
  ESL_DPRINTF1(("#DEBUG: bands validated.\n"));
#endif
  if(debug_level > 0) debug_print_ij_bands(cm);
  if(debug_level > 0) PrintDPCellsSaved_jd(cm, cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax, L);

  return eslOK;
}


/* Function: cp9_IterateSeq2BandsP7B()
 * Date:     EPN, 2026-04-30
 *
 * Purpose:  Like cp9_IterateSeq2Bands(), but uses P7-banded CP9 Forward/Backward.
 *           Runs P7-banded CP9 F/B once (tau-independent), then iteratively
 *           tightens tau/thresh until the resulting CM DP matrix fits within
 *           <size_limit> Mb, or tau reaches <maxtau>.
 *
 * Args:     cm          - the CM
 *           errbuf      - for error messages
 *           dsq         - digitized sequence (1..L)
 *           L           - sequence length
 *           kmin, kmax  - P7-derived per-residue node bands (1..L indexed)
 *           i0, j0      - subsequence bounds in original coords
 *           pass_idx    - pipeline pass index
 *           size_limit  - max allowed CM DP matrix size in Mb
 *           doing_search - TRUE if bands for search, FALSE for alignment
 *           do_sample   - TRUE if we'll sample a parsetree
 *           do_post     - TRUE if we'll do posterior alignment
 *           maxtau      - max allowed cm->tau value
 *           do_pnmono, do_pnmono_print - pnmono flags
 *           ret_Mb      - RETURN: matrix Mb for final bands (can be NULL)
 *
 * Returns:  eslOK on success.
 *           eslERANGE if matrix still exceeds size_limit at maxtau.
 */
int
cp9_IterateSeq2BandsP7B(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, int *kmin, int *kmax,
			int i0, int j0, int pass_idx, float size_limit,
			int doing_search, int do_sample, int do_post,
			double maxtau, int do_pnmono, int do_pnmono_print, float *ret_Mb)
{
  int      status;
  int      do_trunc = cm_pli_PassAllowsTruncation(pass_idx);
  float    cp9mx_Mb = 0.;
  float    hbmx_Mb  = 0.;
  float    tot_Mb;
  int      tau_at_limit     = FALSE;
  int      thresh1_at_limit = (do_trunc) ? FALSE : TRUE;
  int      thresh2_at_limit = (do_trunc) ? FALSE : TRUE;
  CP9_t   *cp9 = NULL;
  CP9_FMX *fmx_f = NULL, *bmx_f = NULL, *pmx_f = NULL; /* float matrices, used by do_trunc path */
  float    sc;

  /* Caller audit: cp9_Seq2BandsP7B and cp9_Seq2BandsP7BF must remain narrow.
   * cp9_Seq2BandsP7B currently has 2 callers: this function (only on the
   * !do_trunc branch below, indirectly via cp9_FBMatrices2BandsP7B) and
   * cm_pipeline.c (search/scan, untouched). cp9_Seq2BandsP7BF has 1 caller:
   * this function (only on the do_trunc branch below, indirectly via
   * cp9_FBMatrices2BandsF). The float entry point is intentionally
   * unreachable from cmscan / cmsearch / non-truncated cmalign.
   */

  /* Contract checks */
  if(cm->cp9map == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_IterateSeq2BandsP7B, cm->cp9map is NULL.");
  if(dsq == NULL)        ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_IterateSeq2BandsP7B, dsq is NULL.");
  if(cm->tau > 0.5)      ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_IterateSeq2BandsP7B, cm->tau (%f) > 0.5.", cm->tau);

  switch(pass_idx) {
    case PLI_PASS_5P_ONLY_FORCE:                cp9 = cm->Rcp9; break;
    case PLI_PASS_3P_ONLY_FORCE:                cp9 = cm->Lcp9; break;
    case PLI_PASS_5P_AND_3P_FORCE:
    case PLI_PASS_5P_AND_3P_ANY:                cp9 = cm->Tcp9; break;
    default:                                    cp9 = cm->cp9;  break;
  }
  if(cp9 == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_IterateSeq2BandsP7B, cp9 is NULL (pass_idx %d).", pass_idx);

  /* Brief 26_0430-215: env-gated DIRECT p7->cp9 band inversion, bypassing the
   * numerically-fragile CP9 F/B posterior (which collapses ALL per-node bands
   * to -1 on a tight band -- see brief 187 comment below). Mirrors the robust
   * HMM-path band conversion (p7_kbands2gbands). Uses p7bands_to_cp9bands with
   * pocc=NULL (extent-based). Test path for whether a tightened band can drive
   * a valid CM alignment without the F/B collapse. */
  if(getenv("P215_DIRECT") != NULL) {
    status = p7bands_to_cp9bands(cm, errbuf, kmin, kmax, L, cm->cp9b, i0, j0, pass_idx, NULL, 0);
    if(ret_Mb != NULL) *ret_Mb = 0.;
    return status;
  }

  if(do_trunc) {
    /* brief 26_0430-150 (144-B Phase 2): checkpointed banded CP9 float F/B for the
     * truncated path (--p7ibv-ckpt, or CP9_CKPT env). Byte-identical bands in
     * O(sqrt(L)*avg_bw) memory. Wired BEFORE the float matrices below are
     * allocated; the tau-ratchet recomputes the checkpointed float F/B each
     * bump (no cached pmx), instrumented via #CP9_CKPTF_TAU. */
    /* brief 26_0430-187 (2026-07-05): temporary, pragmatic reroute. The float non-ckpt
     * F/B below (cp9_ForwardP7BF/cp9_BackwardP7BF/cp9_FBMatrices2BandsF) was
     * shown (briefs 26_0430-185/186) to lose Fwd=Bwd consistency and collapse
     * numerically on genuinely-truncated genome-scale sequences (L~147K,
     * ~15x past its documented L=10K validation bound), producing garbage
     * per-node bands. The checkpointed double-precision twin (brief 26_0430-154) is
     * already validated correct at both sub-genome and genome scale. Rather
     * than write a second, parallel double-precision non-checkpointed kernel
     * right now, we force ALL do_trunc traffic through the checkpointed
     * double path unconditionally, regardless of --p7ibv-ckpt/CP9_CKPT. This
     * CP9-level F/B is already p7-banded (small footprint either way), so
     * this costs at most some recomputation overhead, not memory. The float
     * non-ckpt code below is intentionally left in place, unmodified, but is
     * now dead code for do_trunc pending a real double-precision
     * non-checkpointed rewrite (a wanted future item, out of scope here). */
    if(TRUE) {
      /* brief 26_0430-167 (tau-ratchet single-pass): the old while(1) loop recomputed the
       * WHOLE checkpointed float F/B on every bump (up to ~26 passes). The F/B is
       * tau/thresh-independent (brief 26_0430-166 Q3), so the driver runs it ONCE: step 0
       * via the single-call path (zero overhead for the common 0-bump case), and
       * if that doesn't fit, ONE multi-threshold F/B sweep over the whole ratchet
       * grid + scan-and-pick. Output is byte-identical to the old loop. */
      int nbump = 0;
      status = cp9_IterateSeq2BandsP7BF_chk_multi(cm, errbuf, cp9, dsq, L, kmin, kmax, i0, j0,
                                                  pass_idx, size_limit, doing_search, do_sample,
                                                  do_post, maxtau, do_pnmono, do_pnmono_print,
                                                  &nbump, &hbmx_Mb);
      if(status != eslOK && status != eslERANGE) goto ERROR;
      if(getenv("CP9_CKPT_VERBOSE") != NULL) fprintf(stderr, "#CP9_CKPTF_TAU L=%d tau_bumps=%d (single-pass tabulation; one checkpointed float F/B)\n", L, nbump);
      if(ret_Mb != NULL) *ret_Mb = hbmx_Mb;
      return status; /* eslOK or eslERANGE, exactly as the old loop returned */
    }

    /* Float path: use float-DP CP9 F/B/Posterior to avoid the ~3% per-cell
     * precision drift that destabilizes pocc-based sp/ep prediction in
     * truncated mode. Allocate local float matrices; we do NOT reuse
     * cm->cp9_mx / cm->cp9_bmx (those are int CP9_MX). */
    if((fmx_f = CreateCP9FMatrix(1, cp9->M)) == NULL) ESL_XFAIL(eslEMEM, errbuf, "cp9_IterateSeq2BandsP7B: OOM allocating fmx_f");
    if((bmx_f = CreateCP9FMatrix(1, cp9->M)) == NULL) ESL_XFAIL(eslEMEM, errbuf, "cp9_IterateSeq2BandsP7B: OOM allocating bmx_f");
    if((pmx_f = CreateCP9FMatrix(1, cp9->M)) == NULL) ESL_XFAIL(eslEMEM, errbuf, "cp9_IterateSeq2BandsP7B: OOM allocating pmx_f");

    /* Phase 1: float P7-banded F/B (tau-independent), run once. */
    if((status = cp9_ForwardP7BF (cp9, errbuf, fmx_f, dsq, L, kmin, kmax, &sc))   != eslOK) goto ERROR;
    if((status = cp9_BackwardP7BF(cp9, errbuf, bmx_f, dsq, L, kmin, kmax, NULL))  != eslOK) goto ERROR;

    /* Phase 2: iterate tau/thresh until matrix fits. */
    while(1) {
      if((status = cp9_FBMatrices2BandsF(cm, errbuf, cp9, fmx_f, bmx_f, pmx_f, dsq, cm->cp9b,
					 kmin, kmax, L, i0, j0, pass_idx, 0, do_pnmono, do_pnmono_print)) != eslOK) goto ERROR;
      if(doing_search) {
	if(do_trunc) { if((status = cm_tr_hb_mx_SizeNeeded(cm, errbuf, cm->cp9b, j0-i0+1, NULL, NULL, NULL, NULL, &hbmx_Mb)) != eslOK) goto ERROR; }
	else         { if((status = cm_hb_mx_SizeNeeded   (cm, errbuf, cm->cp9b, j0-i0+1, NULL, &hbmx_Mb)) != eslOK) goto ERROR; }
      }
      else {
	if(do_trunc) { status = cm_TrAlignSizeNeededHB(cm, errbuf, j0-i0+1, size_limit, do_sample, do_post, NULL, NULL, NULL, &cp9mx_Mb, &hbmx_Mb, &tot_Mb); }
	else         { status = cm_AlignSizeNeededHB  (cm, errbuf, j0-i0+1, size_limit, do_sample, do_post, NULL, NULL, NULL, &cp9mx_Mb, &hbmx_Mb, &tot_Mb); }
	if(status != eslOK && status != eslERANGE) goto ERROR;
      }
      if(hbmx_Mb < size_limit)                                  break;
      if(tau_at_limit && thresh1_at_limit && thresh2_at_limit)  break;
      if(! tau_at_limit) { cm->tau *= TAU_MULTIPLIER; if(cm->tau >= maxtau) { cm->tau = maxtau; tau_at_limit = TRUE; } }
      if(! thresh1_at_limit) { cm->cp9b->thresh1 += DELTA_CP9BANDS_THRESH1; if(cm->cp9b->thresh1 >= MAX_CP9BANDS_THRESH1) { cm->cp9b->thresh1 = MAX_CP9BANDS_THRESH1; thresh1_at_limit = TRUE; } }
      if(! thresh2_at_limit) { cm->cp9b->thresh2 -= DELTA_CP9BANDS_THRESH2; if(cm->cp9b->thresh2 <= MIN_CP9BANDS_THRESH2) { cm->cp9b->thresh2 = MIN_CP9BANDS_THRESH2; thresh2_at_limit = TRUE; } }
    }

    FreeCP9FMatrix(fmx_f); fmx_f = NULL;
    FreeCP9FMatrix(bmx_f); bmx_f = NULL;
    FreeCP9FMatrix(pmx_f); pmx_f = NULL;
    if(ret_Mb != NULL) *ret_Mb = hbmx_Mb;
    if(hbmx_Mb > size_limit) return eslERANGE;
    return eslOK;
  }

  /* Non-truncated path (brief 26_0628-073 diagnosis, brief 26_0628-074 fix):
   * route band derivation through the SAME double-precision checkpointed CP9
   * F/B that the do_trunc path uses (cp9_IterateSeq2BandsP7BF_chk_multi),
   * regardless of do_trunc. The legacy int CP9 F/B (cp9_ForwardP7B_OLD_WITH_EL
   * / cp9_BackwardP7B / cp9_FBMatrices2BandsP7B) computed diffuse posterior
   * occupancies for low-bp-density models, producing cp9b bands up to ~7500x
   * wider than this kernel on the identical p7 input (norovirus: 170M vs 22.5K
   * hd_needed) and a 20-3800x --notrunc slowdown (briefs 26_0628-072/26_0628-073).
   * The double-ckpt kernel honors tau exactly as --trunc mode already does;
   * banded parses then score slightly lower (brief 26_0628-073 measured -0.16 to
   * -0.89 bits, alignment span preserved, no clipping) -- the correct, expected
   * consequence of banding, signed off in brief 26_0628-074.
   *
   * Both the int-ckpt escape hatch (--p7ibv-ckpt / CP9_CKPT) and the default int
   * path are removed here; like do_trunc, this path now ignores --p7ibv-ckpt /
   * CP9_CKPT (the double-ckpt kernel is itself checkpointed, so their
   * O(sqrt(L)*avg_bw) memory goal is still met). The int helper functions stay
   * defined -- cp9_IterateSeq2BandsP7BF_chk_multi calls cp9_ForwardP7B_OLD_WITH_EL
   * / cp9_BackwardP7B internally. */
  {
    int nbump = 0;
    status = cp9_IterateSeq2BandsP7BF_chk_multi(cm, errbuf, cp9, dsq, L, kmin, kmax, i0, j0,
                                                pass_idx, size_limit, doing_search, do_sample,
                                                do_post, maxtau, do_pnmono, do_pnmono_print,
                                                &nbump, &hbmx_Mb);
    if(status != eslOK && status != eslERANGE) goto ERROR;
    if(getenv("CP9_CKPT_VERBOSE") != NULL) fprintf(stderr, "#CP9_CKPTF_TAU L=%d tau_bumps=%d (non-truncated; one checkpointed double F/B)\n", L, nbump);
    if(ret_Mb != NULL) *ret_Mb = hbmx_Mb;
    return status;
  }

 ERROR:
  if(fmx_f) FreeCP9FMatrix(fmx_f);
  if(bmx_f) FreeCP9FMatrix(bmx_f);
  if(pmx_f) FreeCP9FMatrix(pmx_f);
  if(ret_Mb != NULL) *ret_Mb = 0.;
  return status;
}


/* Function: p7bands_to_cp9bands
 * Date    : EPN, 2026-03-26
 *
 * Purpose:  Convert p7 HMM bands (kmin[i]/kmax[i] per residue) directly to
 *           CP9 HMM bands (pn_min/pn_max per node), bypassing the O(L*M)
 *           CP9 Forward/Backward DP.
 *
 *           For each CP9 node k, the band covers all positions i where k is
 *           within the p7 band: kmin[i] <= k <= kmax[i].  This is guaranteed
 *           to produce bands at least as permissive (wide) as the p7-banded
 *           CP9 F/B approach, because CP9 F/B constrained by p7 bands cannot
 *           assign non-zero posterior to positions outside those bands.
 *
 *           NOTE: This function is CURRENTLY UNUSED (no active callers on the
 *           defppp dispatch path as of 2026-04-21). It is inferior to the
 *           two-phase p7banded_post_to_pn_bands + p7pn_bands_to_cp9cm_bands
 *           approach used by defppp, because it has no access to F/B posteriors
 *           and therefore cannot compute per-node occupancy for the cp9-style
 *           thresh1/thresh2 sp/ep derivation (it falls back to the legacy
 *           extent-based sp1=sp2/ep1=ep2 collapse when pocc==NULL).
 *           Consider removing this function if no future caller materializes.
 *
 * Args:     cm          - the covariance model
 *           errbuf      - char buffer for reporting errors
 *           kmin        - p7 bands: kmin[i] = min node k for residue i [0..L]
 *           kmax        - p7 bands: kmax[i] = max node k for residue i [0..L]
 *           L           - length of target subsequence (1..L)
 *           cp9b        - PRE-ALLOCATED CP9 bands, filled here
 *           i0          - first position in original sequence coords
 *           j0          - final position in original sequence coords
 *           pass_idx    - pipeline pass index, determines truncation mode
 *           pocc        - [0..M] per-node match-posterior occupancy OR NULL.
 *                         When non-NULL, enables cp9-style thresh1/thresh2 sweep
 *                         for sp1/sp2/ep1/ep2 (same as p7pn_bands_to_cp9cm_bands).
 *                         When NULL, falls back to legacy extent-based collapse
 *                         (sp1=sp2, ep1=ep2). This function cannot compute pocc
 *                         internally (no F/B matrices available); callers must
 *                         supply it from an upstream posterior computation or
 *                         pass NULL to get legacy behavior.
 *           debug_level - verbosity level for debugging printf()s
 *
 * Returns:  eslOK on success
 */
int
p7bands_to_cp9bands(CM_t *cm, char *errbuf, int *kmin, int *kmax, int L,
		    CP9Bands_t *cp9b, int i0, int j0, int pass_idx,
		    const float *pocc, int debug_level)
{
  int   status;
  int   i, k;
  int   M            = cp9b->hmm_M;
  int   do_old_hmm2ij;
  int   do_trunc;
  CP9_t *cp9;

  /* Contract checks */
  if(cm->cp9map == NULL)
    ESL_FAIL(eslEINCOMPAT, errbuf, "p7bands_to_cp9bands, cm->cp9map is NULL.\n");
  if(!((cm->align_opts & CM_ALIGN_HBANDED) || (cm->search_opts & CM_SEARCH_HBANDED)))
    ESL_FAIL(eslEINCOMPAT, errbuf, "p7bands_to_cp9bands, neither CM_ALIGN_HBANDED nor CM_SEARCH_HBANDED is set.\n");
  if(cm->tau > 0.5)
    ESL_FAIL(eslEINCOMPAT, errbuf, "p7bands_to_cp9bands, cm->tau (%f) > 0.5.\n", cm->tau);

  do_old_hmm2ij = ((cm->align_opts & CM_ALIGN_HMM2IJOLD) || (cm->search_opts & CM_SEARCH_HMM2IJOLD)) ? TRUE : FALSE;
  do_trunc      = cm_pli_PassAllowsTruncation(pass_idx);

  switch(pass_idx) {
    case PLI_PASS_5P_ONLY_FORCE:                cp9 = cm->Rcp9; break;
    case PLI_PASS_3P_ONLY_FORCE:                cp9 = cm->Lcp9; break;
    case PLI_PASS_5P_AND_3P_FORCE:
    case PLI_PASS_5P_AND_3P_ANY:                cp9 = cm->Tcp9; break;
    default:                                    cp9 = cm->cp9;  break;
  }
  if(cp9 == NULL)
    ESL_FAIL(eslEINCOMPAT, errbuf, "p7bands_to_cp9bands, cp9 is NULL (pass_idx %d).\n", pass_idx);

  /* Initialize all bands to "unset" sentinel values:
   * pn_min_* = L+2  (larger than any valid position, signals "not yet seen from left")
   * pn_max_* = -1   (smaller than any valid position, signals "not yet seen from right") */
  for(k = 0; k <= M; k++) {
    cp9b->pn_min_m[k] = L + 2;
    cp9b->pn_max_m[k] = -1;
    cp9b->pn_min_i[k] = L + 2;
    cp9b->pn_max_i[k] = -1;
    cp9b->pn_min_d[k] = L + 2;
    cp9b->pn_max_d[k] = -1;
  }

  /* Direct inversion: sweep i=1..L for match/insert/delete; handle i=0 specially.
   * At i=0 (the "begin" position before any residues):
   *   - Only M_0 (begin state) has non-zero posterior => set pn_min/max_m[0] = 0
   *   - D_k for k>0 can have non-zero posterior (entering deletes from begin) => update d bands
   *   - M_k and I_k for k>0 are -INFTY at i=0 => must NOT update m/i bands for k>0
   *   Setting pn_min_m[k] = 0 for k>0 would cause imin[v] = 0 in cp9_HMM2ijBands,
   *   which violates the requirement imin[v] >= i0 >= 1. */
  cp9b->pn_min_m[0] = cp9b->pn_max_m[0] = 0;  /* M_0 (begin) is always at position 0 */
  for(k = 1; k <= kmax[0]; k++) {               /* D_k for k>0 may be active at i=0 */
    if(0 < cp9b->pn_min_d[k]) cp9b->pn_min_d[k] = 0;
    if(0 > cp9b->pn_max_d[k]) cp9b->pn_max_d[k] = 0;
  }

  for(i = 1; i <= L; i++) {
    for(k = kmin[i]; k <= kmax[i]; k++) {
      if(i < cp9b->pn_min_m[k]) cp9b->pn_min_m[k] = i;
      if(i > cp9b->pn_max_m[k]) cp9b->pn_max_m[k] = i;
      if(i < cp9b->pn_min_i[k]) cp9b->pn_min_i[k] = i;
      if(i > cp9b->pn_max_i[k]) cp9b->pn_max_i[k] = i;
      if(i < cp9b->pn_min_d[k]) cp9b->pn_min_d[k] = i;
      if(i > cp9b->pn_max_d[k]) cp9b->pn_max_d[k] = i;
    }
  }

  /* Convert unset states (min still > max) to -1 sentinel used by downstream code */
  for(k = 0; k <= M; k++) {
    if(cp9b->pn_min_m[k] > cp9b->pn_max_m[k]) cp9b->pn_min_m[k] = cp9b->pn_max_m[k] = -1;
    if(cp9b->pn_min_i[k] > cp9b->pn_max_i[k]) cp9b->pn_min_i[k] = cp9b->pn_max_i[k] = -1;
    if(cp9b->pn_min_d[k] > cp9b->pn_max_d[k]) cp9b->pn_min_d[k] = cp9b->pn_max_d[k] = -1;
  }
  cp9b->pn_min_d[0] = cp9b->pn_max_d[0] = -1; /* D_0 does not exist */

  cp9b->tau = cm->tau;

  /* Shift HMM bands from 1..L to i0..j0 coordinate system if needed */
  if(i0 != 1) {
    int offset = i0 - 1;
    for(k = 0; k <= M; k++) {
      if(cp9b->pn_min_m[k] != -1) { cp9b->pn_min_m[k] += offset; cp9b->pn_max_m[k] += offset; }
      if(cp9b->pn_min_i[k] != -1) { cp9b->pn_min_i[k] += offset; cp9b->pn_max_i[k] += offset; }
      if(cp9b->pn_min_d[k] != -1) { cp9b->pn_min_d[k] += offset; cp9b->pn_max_d[k] += offset; }
    }
  }

  /* Set truncation candidate valid arrays.
   * Two modes depending on whether the caller supplies per-node posterior occupancy:
   * - pocc != NULL: cp9-style thresh1/thresh2 sweep (sp1 != sp2, ep1 != ep2 possible),
   *   Rmarg/Lmarg from union over both. Mirrors p7pn_bands_to_cp9cm_bands.
   * - pocc == NULL: legacy extent-based collapse (sp1=sp2, ep1=ep2). */
  if(do_trunc) {
    if(pocc != NULL) {
      cp9b->sp1 = cp9b->sp2 = -1;
      for(k = 1; k <= M; k++) {
        if(cp9b->pn_min_m[k] == -1 && cp9b->pn_min_i[k] == -1 && cp9b->pn_min_d[k] == -1) continue;
        if(cp9b->sp1 == -1 && pocc[k] > cp9b->thresh1) cp9b->sp1 = k;
        if(cp9b->sp2 == -1 && pocc[k] > cp9b->thresh2) cp9b->sp2 = k;
        if(cp9b->sp1 != -1 && cp9b->sp2 != -1) break;
      }
      if(cp9b->sp1 == -1) cp9b->sp1 = M + 1;
      if(cp9b->sp2 == -1) cp9b->sp2 = M + 1;

      if(cp9b->sp1 == M + 1 && cp9b->sp2 == M + 1) {
        cp9b->ep1 = cp9b->ep2 = 0;
      } else {
        cp9b->ep1 = cp9b->ep2 = -1;
        for(k = M; k >= 1; k--) {
          if(cp9b->pn_min_m[k] == -1 && cp9b->pn_min_i[k] == -1 && cp9b->pn_min_d[k] == -1) continue;
          if(cp9b->ep1 == -1 && pocc[k] > cp9b->thresh1) cp9b->ep1 = k;
          if(cp9b->ep2 == -1 && pocc[k] > cp9b->thresh2) cp9b->ep2 = k;
          if(cp9b->ep1 != -1 && cp9b->ep2 != -1) break;
        }
        if(cp9b->ep1 == -1) cp9b->ep1 = 0;
        if(cp9b->ep2 == -1) cp9b->ep2 = 0;
      }

      if(cp9b->sp1 == M + 1) {
        cp9b->Rmarg_imin = i0; cp9b->Rmarg_imax = j0;
      } else {
        int rmin = INT_MAX, rmax = INT_MIN;
        if(cp9b->sp1 != (M+1) && cp9b->pn_min_m[cp9b->sp1] >= 0) rmin = ESL_MIN(rmin, cp9b->pn_min_m[cp9b->sp1]);
        if(cp9b->sp1 != (M+1) && cp9b->pn_min_i[cp9b->sp1] >= 0) rmin = ESL_MIN(rmin, cp9b->pn_min_i[cp9b->sp1]);
        if(cp9b->sp1 != (M+1) && cp9b->pn_min_d[cp9b->sp1] >= 0) rmin = ESL_MIN(rmin, cp9b->pn_min_d[cp9b->sp1]);
        if(cp9b->sp2 != (M+1) && cp9b->pn_min_m[cp9b->sp2] >= 0) rmin = ESL_MIN(rmin, cp9b->pn_min_m[cp9b->sp2]);
        if(cp9b->sp2 != (M+1) && cp9b->pn_min_i[cp9b->sp2] >= 0) rmin = ESL_MIN(rmin, cp9b->pn_min_i[cp9b->sp2]);
        if(cp9b->sp2 != (M+1) && cp9b->pn_min_d[cp9b->sp2] >= 0) rmin = ESL_MIN(rmin, cp9b->pn_min_d[cp9b->sp2]);
        if(rmin == INT_MAX || cp9b->sp1 == (M+1) || cp9b->sp2 == (M+1)) rmin = i0;
        cp9b->Rmarg_imin = ESL_MAX(i0, ESL_MIN(j0 + 1, rmin));

        if(cp9b->sp1 != (M+1) && cp9b->pn_max_m[cp9b->sp1] >= 0) rmax = ESL_MAX(rmax, cp9b->pn_max_m[cp9b->sp1]);
        if(cp9b->sp1 != (M+1) && cp9b->pn_max_i[cp9b->sp1] >= 0) rmax = ESL_MAX(rmax, cp9b->pn_max_i[cp9b->sp1]);
        if(cp9b->sp1 != (M+1) && cp9b->pn_max_d[cp9b->sp1] >= 0) rmax = ESL_MAX(rmax, cp9b->pn_max_d[cp9b->sp1]);
        if(cp9b->sp2 != (M+1) && cp9b->pn_max_m[cp9b->sp2] >= 0) rmax = ESL_MAX(rmax, cp9b->pn_max_m[cp9b->sp2]);
        if(cp9b->sp2 != (M+1) && cp9b->pn_max_i[cp9b->sp2] >= 0) rmax = ESL_MAX(rmax, cp9b->pn_max_i[cp9b->sp2]);
        if(cp9b->sp2 != (M+1) && cp9b->pn_max_d[cp9b->sp2] >= 0) rmax = ESL_MAX(rmax, cp9b->pn_max_d[cp9b->sp2]);
        if(rmax == INT_MIN || cp9b->sp1 == (M+1) || cp9b->sp2 == (M+1)) rmax = j0 + 1;
        cp9b->Rmarg_imax = ESL_MAX(i0, ESL_MIN(j0 + 1, rmax));
      }

      if(cp9b->ep1 == 0) {
        cp9b->Lmarg_jmin = i0 - 1; cp9b->Lmarg_jmax = j0;
      } else {
        int lmin = INT_MAX, lmax = INT_MIN;
        if(cp9b->ep1 != 0 && cp9b->pn_min_m[cp9b->ep1] >= 0) lmin = ESL_MIN(lmin, cp9b->pn_min_m[cp9b->ep1]);
        if(cp9b->ep1 != 0 && cp9b->pn_min_i[cp9b->ep1] >= 0) lmin = ESL_MIN(lmin, cp9b->pn_min_i[cp9b->ep1]);
        if(cp9b->ep1 != 0 && cp9b->pn_min_d[cp9b->ep1] >= 0) lmin = ESL_MIN(lmin, cp9b->pn_min_d[cp9b->ep1] - 1);
        if(cp9b->ep2 != 0 && cp9b->pn_min_m[cp9b->ep2] >= 0) lmin = ESL_MIN(lmin, cp9b->pn_min_m[cp9b->ep2]);
        if(cp9b->ep2 != 0 && cp9b->pn_min_i[cp9b->ep2] >= 0) lmin = ESL_MIN(lmin, cp9b->pn_min_i[cp9b->ep2]);
        if(cp9b->ep2 != 0 && cp9b->pn_min_d[cp9b->ep2] >= 0) lmin = ESL_MIN(lmin, cp9b->pn_min_d[cp9b->ep2] - 1);
        if(lmin == INT_MAX || cp9b->ep1 == 0 || cp9b->ep2 == 0) lmin = i0 - 1;
        cp9b->Lmarg_jmin = ESL_MAX(i0 - 1, ESL_MIN(j0, lmin));

        if(cp9b->ep1 != 0 && cp9b->pn_max_m[cp9b->ep1] >= 0) lmax = ESL_MAX(lmax, cp9b->pn_max_m[cp9b->ep1]);
        if(cp9b->ep1 != 0 && cp9b->pn_max_i[cp9b->ep1] >= 0) lmax = ESL_MAX(lmax, cp9b->pn_max_i[cp9b->ep1]);
        if(cp9b->ep1 != 0 && cp9b->pn_max_d[cp9b->ep1] >= 0) lmax = ESL_MAX(lmax, cp9b->pn_max_d[cp9b->ep1] - 1);
        if(cp9b->ep2 != 0 && cp9b->pn_max_m[cp9b->ep2] >= 0) lmax = ESL_MAX(lmax, cp9b->pn_max_m[cp9b->ep2]);
        if(cp9b->ep2 != 0 && cp9b->pn_max_i[cp9b->ep2] >= 0) lmax = ESL_MAX(lmax, cp9b->pn_max_i[cp9b->ep2]);
        if(cp9b->ep2 != 0 && cp9b->pn_max_d[cp9b->ep2] >= 0) lmax = ESL_MAX(lmax, cp9b->pn_max_d[cp9b->ep2] - 1);
        if(lmax == INT_MIN || cp9b->ep1 == 0 || cp9b->ep2 == 0) lmax = j0;
        cp9b->Lmarg_jmax = ESL_MAX(i0 - 1, ESL_MIN(j0, lmax));
      }
    } else {
      /* Legacy: no posterior available; extent-based collapse (sp1=sp2, ep1=ep2). */
      int sp = M + 1, ep = 0;
      for(i = 1; i <= L; i++) {
        if(kmin[i] < sp) sp = kmin[i];
        if(kmax[i] > ep) ep = kmax[i];
      }
      if(sp < 1)     sp = 1;
      if(sp > M + 1) sp = M + 1;
      if(ep < 0)     ep = 0;
      if(ep > M)     ep = M;
      cp9b->sp1 = cp9b->sp2 = sp;
      cp9b->ep1 = cp9b->ep2 = ep;

      if(cp9b->sp1 == M + 1) {
        cp9b->Rmarg_imin = i0; cp9b->Rmarg_imax = j0;
      } else {
        int rmarg_imin = INT_MAX, rmarg_imax = INT_MIN;
        if(cp9b->pn_min_m[sp] >= 0) rmarg_imin = ESL_MIN(rmarg_imin, cp9b->pn_min_m[sp]);
        if(cp9b->pn_min_i[sp] >= 0) rmarg_imin = ESL_MIN(rmarg_imin, cp9b->pn_min_i[sp]);
        if(cp9b->pn_min_d[sp] >= 0) rmarg_imin = ESL_MIN(rmarg_imin, cp9b->pn_min_d[sp]);
        if(rmarg_imin == INT_MAX) rmarg_imin = i0;
        cp9b->Rmarg_imin = ESL_MAX(i0, ESL_MIN(j0 + 1, rmarg_imin));

        if(cp9b->pn_max_m[sp] >= 0) rmarg_imax = ESL_MAX(rmarg_imax, cp9b->pn_max_m[sp]);
        if(cp9b->pn_max_i[sp] >= 0) rmarg_imax = ESL_MAX(rmarg_imax, cp9b->pn_max_i[sp]);
        if(cp9b->pn_max_d[sp] >= 0) rmarg_imax = ESL_MAX(rmarg_imax, cp9b->pn_max_d[sp]);
        if(rmarg_imax == INT_MIN) rmarg_imax = j0 + 1;
        cp9b->Rmarg_imax = ESL_MAX(i0, ESL_MIN(j0 + 1, rmarg_imax));
      }

      if(cp9b->ep1 == 0) {
        cp9b->Lmarg_jmin = i0 - 1; cp9b->Lmarg_jmax = j0;
      } else {
        int lmarg_jmin = INT_MAX, lmarg_jmax = INT_MIN;
        if(cp9b->pn_min_m[ep] >= 0) lmarg_jmin = ESL_MIN(lmarg_jmin, cp9b->pn_min_m[ep]);
        if(cp9b->pn_min_i[ep] >= 0) lmarg_jmin = ESL_MIN(lmarg_jmin, cp9b->pn_min_i[ep]);
        if(cp9b->pn_min_d[ep] >= 0) lmarg_jmin = ESL_MIN(lmarg_jmin, cp9b->pn_min_d[ep] - 1);
        if(lmarg_jmin == INT_MAX) lmarg_jmin = i0 - 1;
        cp9b->Lmarg_jmin = ESL_MAX(i0 - 1, ESL_MIN(j0, lmarg_jmin));

        if(cp9b->pn_max_m[ep] >= 0) lmarg_jmax = ESL_MAX(lmarg_jmax, cp9b->pn_max_m[ep]);
        if(cp9b->pn_max_i[ep] >= 0) lmarg_jmax = ESL_MAX(lmarg_jmax, cp9b->pn_max_i[ep]);
        if(cp9b->pn_max_d[ep] >= 0) lmarg_jmax = ESL_MAX(lmarg_jmax, cp9b->pn_max_d[ep] - 1);
        if(lmarg_jmax == INT_MIN) lmarg_jmax = j0;
        cp9b->Lmarg_jmax = ESL_MAX(i0 - 1, ESL_MIN(j0, lmarg_jmax));
      }
    }

    if((status = cp9_MarginalCandidatesFromStartEndPositions(cm, cp9b, pass_idx, errbuf)) != eslOK) return status;
  } else {
    esl_vec_ISet(cp9b->Jvalid, cm->M + 1, TRUE);
    esl_vec_ISet(cp9b->Lvalid, cm->M + 1, FALSE);
    esl_vec_ISet(cp9b->Rvalid, cm->M + 1, FALSE);
    esl_vec_ISet(cp9b->Tvalid, cm->M + 1, FALSE);
  }

  /* HMM bands -> CM ij bands */
  if(do_old_hmm2ij) {
    if((status = cp9_HMM2ijBands_OLD(cm, errbuf, cm->cp9b, cm->cp9map, i0, j0, TRUE, debug_level)) != eslOK) return status;
  } else {
    if((status = cp9_HMM2ijBands(cm, errbuf, cp9, cm->cp9b, cm->cp9map, i0, j0, TRUE, do_trunc, debug_level)) != eslOK) return status;
  }

  /* CM ij bands -> CM d bands */
  if((status = cp9_GrowHDBands(cp9b, errbuf)) != eslOK) return status;
  ij2d_bands(cm, cp9b, do_trunc, debug_level);

#if eslDEBUGLEVEL >= 1
  if((status = cp9_ValidateBands(cm, errbuf, cp9b, i0, j0, do_trunc)) != eslOK) return status;
  ESL_DPRINTF1(("#DEBUG: p7bands_to_cp9bands bands validated.\n"));
#endif
  if(debug_level > 0) debug_print_ij_bands(cm);

  return eslOK;
}


/* Function: p7banded_post_to_cp9bands
 * Date    : EPN, 2026-03-26
 *
 * Purpose:  Derive CP9 HMM bands from p7 glocal banded Forward and Backward
 *           posterior probabilities, bypassing CP9 Forward/Backward entirely.
 *
 *           For each CP9 node k and state type (M/I/D), scans the banded
 *           (i,k) cells and includes position i in the band for node k if
 *           the posterior probability exceeds <thresh>:
 *
 *             post_M(i,k) = expf(gxfb->dp[...M...] + gxbb->dp[...M...] - fwdsc)
 *
 *           The p7 glocal F and B matrices are computed during F5 (envelope
 *           definition) and are already available in pli->gxfb / pli->gxbb.
 *           This approach gives tighter bands than the Viterbi-inversion in
 *           p7bands_to_cp9bands() because posteriors concentrate probability
 *           mass near the most likely alignment path, while still bypassing
 *           the O(L*M) CP9 F/B DP.
 *
 *           Delete bands are set equal to match bands for each node k (D
 *           states do not emit, so we mirror the match coverage).
 *
 *           The i=0 special case is handled identically to p7bands_to_cp9bands:
 *           M_0 (begin state) is always at position 0; D_k for k in the first
 *           band row may be active at position 0.
 *
 * Args:     cm          - the covariance model
 *           errbuf      - char buffer for reporting errors
 *           gxfb        - p7 banded glocal Forward matrix (from F5)
 *           gxbb        - p7 banded glocal Backward matrix (from F5)
 *           fwdsc       - Forward score in nats (from p7_GForwardBanded)
 *           bnd         - band structure for gxfb/gxbb (passed explicitly to avoid use-after-free;
 *                         caller keeps it alive via pli->p7bnd)
 *           ws          - absolute start of the window that gxfb/gxbb were computed over
 *                         (1-indexed in original sequence); used to map window-relative
 *                         gxfb rows to absolute/envelope positions
 *           L           - length of target subsequence (= j0 - i0 + 1)
 *           cp9b        - PRE-ALLOCATED CP9 bands, filled here
 *           i0          - first position (absolute) of the envelope
 *           j0          - final position (absolute) of the envelope
 *           pass_idx    - pipeline pass index, determines truncation mode
 *           thresh      - posterior probability threshold (e.g. 1e-5):
 *                         include position i for node k if post(i,k) >= thresh
 *           debug_level - verbosity level for debugging printf()s
 *
 * Returns:  eslOK on success
 */
int
p7banded_post_to_cp9bands(CM_t *cm, char *errbuf,
			  P7_GMXB *gxfb, P7_GMXB *gxbb, float fwdsc,
			  P7_GBANDS *bnd, int ws, int L, CP9Bands_t *cp9b, int i0, int j0,
			  int pass_idx, float thresh, int debug_level)
{
  int          status;
  int          g, i, k;
  int          M            = cp9b->hmm_M;
  int          do_old_hmm2ij;
  int          do_trunc;
  CP9_t       *cp9;
  int          ia, ib;
  int          kac, kbc;
  int         *bnd_ip;
  int         *bnd_kp;
  float const *fwd_dp;
  float const *bck_dp;
  float        log_thresh;
  float        fM, bM, fI, bI;
  int          nk;
  int          first_kac, first_kbc; /* band of first row, for i=0 D_k handling */
  float       *pocc = NULL;          /* [0..M] per-node match-posterior occupancy */

  /* Contract checks */
  if(cm->cp9map == NULL)
    ESL_FAIL(eslEINCOMPAT, errbuf, "p7banded_post_to_cp9bands: cm->cp9map is NULL.\n");
  if(!((cm->align_opts & CM_ALIGN_HBANDED) || (cm->search_opts & CM_SEARCH_HBANDED)))
    ESL_FAIL(eslEINCOMPAT, errbuf, "p7banded_post_to_cp9bands: neither CM_ALIGN_HBANDED nor CM_SEARCH_HBANDED is set.\n");
  if(cm->tau > 0.5)
    ESL_FAIL(eslEINCOMPAT, errbuf, "p7banded_post_to_cp9bands: cm->tau (%f) > 0.5.\n", cm->tau);
  if(gxfb == NULL || gxbb == NULL)
    ESL_FAIL(eslEINCOMPAT, errbuf, "p7banded_post_to_cp9bands: gxfb or gxbb is NULL.\n");
  if(bnd == NULL)
    ESL_FAIL(eslEINCOMPAT, errbuf, "p7banded_post_to_cp9bands: bnd is NULL.\n");
  if(bnd->nrow == 0)
    ESL_FAIL(eslEINCOMPAT, errbuf, "p7banded_post_to_cp9bands: bnd has no banded rows.\n");

  ESL_ALLOC(pocc, sizeof(float) * (M + 1));
  for(k = 0; k <= M; k++) pocc[k] = 0.0f;

  do_old_hmm2ij = ((cm->align_opts & CM_ALIGN_HMM2IJOLD) || (cm->search_opts & CM_SEARCH_HMM2IJOLD)) ? TRUE : FALSE;
  do_trunc      = cm_pli_PassAllowsTruncation(pass_idx);

  switch(pass_idx) {
    case PLI_PASS_5P_ONLY_FORCE:                cp9 = cm->Rcp9; break;
    case PLI_PASS_3P_ONLY_FORCE:                cp9 = cm->Lcp9; break;
    case PLI_PASS_5P_AND_3P_FORCE:
    case PLI_PASS_5P_AND_3P_ANY:                cp9 = cm->Tcp9; break;
    default:                                    cp9 = cm->cp9;  break;
  }
  if(cp9 == NULL)
    ESL_FAIL(eslEINCOMPAT, errbuf, "p7banded_post_to_cp9bands: cp9 is NULL (pass_idx %d).\n", pass_idx);

  log_thresh = logf(thresh);

  /* Initialize all bands to "unset" sentinel values */
  for(k = 0; k <= M; k++) {
    cp9b->pn_min_m[k] = L + 2;
    cp9b->pn_max_m[k] = -1;
    cp9b->pn_min_i[k] = L + 2;
    cp9b->pn_max_i[k] = -1;
    cp9b->pn_min_d[k] = L + 2;
    cp9b->pn_max_d[k] = -1;
  }

  /* Main posterior sweep over banded (i,k) cells.
   * Walk fwd->dp and bck->dp in parallel (same P7_GBANDS layout).
   * bnd is passed explicitly (not via gxfb->bnd) to avoid use-after-free
   * when the caller transfers bnd ownership to pli->p7bnd. */
  bnd_ip = bnd->imem;
  bnd_kp = bnd->kmem;
  fwd_dp = gxfb->dp;
  bck_dp = gxbb->dp;

  /* Record first row's band for i=0 D_k handling below */
  first_kac = bnd->kmem[0];
  first_kbc = bnd->kmem[1];

  for(g = 0; g < bnd->nseg; g++) {
    ia = *bnd_ip++;
    ib = *bnd_ip++;

    for(i = ia; i <= ib; i++) {
      kac = *bnd_kp++;
      kbc = *bnd_kp++;
      nk  = kbc - kac + 1;

      /* Convert window-relative row i to absolute position, then to
       * envelope-relative position ep (1-indexed within [i0..j0]).
       * Skip rows outside the envelope. */
      int abs_i = ws + i - 1;
      int ep    = abs_i - i0 + 1;   /* 1-indexed within envelope */
      if(abs_i >= i0 && abs_i <= j0) {
	for(k = kac; k <= kbc; k++) {
	  fM = fwd_dp[(k-kac)*p7G_NSCELLS + p7G_M];
	  bM = bck_dp[(k-kac)*p7G_NSCELLS + p7G_M];
	  pocc[k] += expf(fM + bM - fwdsc);

	  if(fM + bM - fwdsc > log_thresh) {
	    if(ep < cp9b->pn_min_m[k]) cp9b->pn_min_m[k] = ep;
	    if(ep > cp9b->pn_max_m[k]) cp9b->pn_max_m[k] = ep;
	    /* Delete bands mirror match bands */
	    if(k > 0) {
	      if(ep < cp9b->pn_min_d[k]) cp9b->pn_min_d[k] = ep;
	      if(ep > cp9b->pn_max_d[k]) cp9b->pn_max_d[k] = ep;
	    }
	  }

	  if(k < M) {
	    fI = fwd_dp[(k-kac)*p7G_NSCELLS + p7G_I];
	    bI = bck_dp[(k-kac)*p7G_NSCELLS + p7G_I];
	    if(fI + bI - fwdsc > log_thresh) {
	      if(ep < cp9b->pn_min_i[k]) cp9b->pn_min_i[k] = ep;
	      if(ep > cp9b->pn_max_i[k]) cp9b->pn_max_i[k] = ep;
	    }
	  }
	}
      }
      fwd_dp += nk * p7G_NSCELLS;
      bck_dp += nk * p7G_NSCELLS;
    }
  }

  /* Handle i=0 special case (same logic as p7bands_to_cp9bands):
   *   M_0 (begin state) is always at position 0.
   *   D_k for k in first_kac..first_kbc may be active at position 0
   *   (entering the model via begin -> delete transitions). */
  cp9b->pn_min_m[0] = cp9b->pn_max_m[0] = 0;
  for(k = first_kac; k <= first_kbc; k++) {
    if(k > 0) {
      if(0 < cp9b->pn_min_d[k]) cp9b->pn_min_d[k] = 0;
      if(0 > cp9b->pn_max_d[k]) cp9b->pn_max_d[k] = 0;
    }
  }

  /* Convert unset states to -1 sentinel */
  for(k = 0; k <= M; k++) {
    if(cp9b->pn_min_m[k] > cp9b->pn_max_m[k]) cp9b->pn_min_m[k] = cp9b->pn_max_m[k] = -1;
    if(cp9b->pn_min_i[k] > cp9b->pn_max_i[k]) cp9b->pn_min_i[k] = cp9b->pn_max_i[k] = -1;
    if(cp9b->pn_min_d[k] > cp9b->pn_max_d[k]) cp9b->pn_min_d[k] = cp9b->pn_max_d[k] = -1;
  }
  cp9b->pn_min_d[0] = cp9b->pn_max_d[0] = -1; /* D_0 does not exist */

  cp9b->tau = cm->tau;

  /* Shift HMM bands from 1..L to i0..j0 coordinate system if needed */
  if(i0 != 1) {
    int offset = i0 - 1;
    for(k = 0; k <= M; k++) {
      if(cp9b->pn_min_m[k] != -1) { cp9b->pn_min_m[k] += offset; cp9b->pn_max_m[k] += offset; }
      if(cp9b->pn_min_i[k] != -1) { cp9b->pn_min_i[k] += offset; cp9b->pn_max_i[k] += offset; }
      if(cp9b->pn_min_d[k] != -1) { cp9b->pn_min_d[k] += offset; cp9b->pn_max_d[k] += offset; }
    }
  }

  /* Set truncation candidate valid arrays using cp9-style thresh1/thresh2 sweep.
   * pocc[k] was accumulated above (unconditionally) as Σ_i exp(fM+bM-fwdsc).
   * Mirrors p7pn_bands_to_cp9cm_bands and hmmband.c:cp9_PredictStartAndEndPositions. */
  if(do_trunc) {
    cp9b->sp1 = cp9b->sp2 = -1;
    for(k = 1; k <= M; k++) {
      if(cp9b->pn_min_m[k] == -1 && cp9b->pn_min_i[k] == -1 && cp9b->pn_min_d[k] == -1) continue;
      if(cp9b->sp1 == -1 && pocc[k] > cp9b->thresh1) cp9b->sp1 = k;
      if(cp9b->sp2 == -1 && pocc[k] > cp9b->thresh2) cp9b->sp2 = k;
      if(cp9b->sp1 != -1 && cp9b->sp2 != -1) break;
    }
    if(cp9b->sp1 == -1) cp9b->sp1 = M + 1;
    if(cp9b->sp2 == -1) cp9b->sp2 = M + 1;

    if(cp9b->sp1 == M + 1 && cp9b->sp2 == M + 1) {
      cp9b->ep1 = cp9b->ep2 = 0;
    } else {
      cp9b->ep1 = cp9b->ep2 = -1;
      for(k = M; k >= 1; k--) {
        if(cp9b->pn_min_m[k] == -1 && cp9b->pn_min_i[k] == -1 && cp9b->pn_min_d[k] == -1) continue;
        if(cp9b->ep1 == -1 && pocc[k] > cp9b->thresh1) cp9b->ep1 = k;
        if(cp9b->ep2 == -1 && pocc[k] > cp9b->thresh2) cp9b->ep2 = k;
        if(cp9b->ep1 != -1 && cp9b->ep2 != -1) break;
      }
      if(cp9b->ep1 == -1) cp9b->ep1 = 0;
      if(cp9b->ep2 == -1) cp9b->ep2 = 0;
    }

    if(cp9b->sp1 == M + 1) {
      cp9b->Rmarg_imin = i0; cp9b->Rmarg_imax = j0;
    } else {
      int rmin = INT_MAX, rmax = INT_MIN;
      if(cp9b->sp1 != (M+1) && cp9b->pn_min_m[cp9b->sp1] >= 0) rmin = ESL_MIN(rmin, cp9b->pn_min_m[cp9b->sp1]);
      if(cp9b->sp1 != (M+1) && cp9b->pn_min_i[cp9b->sp1] >= 0) rmin = ESL_MIN(rmin, cp9b->pn_min_i[cp9b->sp1]);
      if(cp9b->sp1 != (M+1) && cp9b->pn_min_d[cp9b->sp1] >= 0) rmin = ESL_MIN(rmin, cp9b->pn_min_d[cp9b->sp1]);
      if(cp9b->sp2 != (M+1) && cp9b->pn_min_m[cp9b->sp2] >= 0) rmin = ESL_MIN(rmin, cp9b->pn_min_m[cp9b->sp2]);
      if(cp9b->sp2 != (M+1) && cp9b->pn_min_i[cp9b->sp2] >= 0) rmin = ESL_MIN(rmin, cp9b->pn_min_i[cp9b->sp2]);
      if(cp9b->sp2 != (M+1) && cp9b->pn_min_d[cp9b->sp2] >= 0) rmin = ESL_MIN(rmin, cp9b->pn_min_d[cp9b->sp2]);
      if(rmin == INT_MAX || cp9b->sp1 == (M+1) || cp9b->sp2 == (M+1)) rmin = i0;
      cp9b->Rmarg_imin = ESL_MAX(i0, ESL_MIN(j0 + 1, rmin));

      if(cp9b->sp1 != (M+1) && cp9b->pn_max_m[cp9b->sp1] >= 0) rmax = ESL_MAX(rmax, cp9b->pn_max_m[cp9b->sp1]);
      if(cp9b->sp1 != (M+1) && cp9b->pn_max_i[cp9b->sp1] >= 0) rmax = ESL_MAX(rmax, cp9b->pn_max_i[cp9b->sp1]);
      if(cp9b->sp1 != (M+1) && cp9b->pn_max_d[cp9b->sp1] >= 0) rmax = ESL_MAX(rmax, cp9b->pn_max_d[cp9b->sp1]);
      if(cp9b->sp2 != (M+1) && cp9b->pn_max_m[cp9b->sp2] >= 0) rmax = ESL_MAX(rmax, cp9b->pn_max_m[cp9b->sp2]);
      if(cp9b->sp2 != (M+1) && cp9b->pn_max_i[cp9b->sp2] >= 0) rmax = ESL_MAX(rmax, cp9b->pn_max_i[cp9b->sp2]);
      if(cp9b->sp2 != (M+1) && cp9b->pn_max_d[cp9b->sp2] >= 0) rmax = ESL_MAX(rmax, cp9b->pn_max_d[cp9b->sp2]);
      if(rmax == INT_MIN || cp9b->sp1 == (M+1) || cp9b->sp2 == (M+1)) rmax = j0 + 1;
      cp9b->Rmarg_imax = ESL_MAX(i0, ESL_MIN(j0 + 1, rmax));
    }

    if(cp9b->ep1 == 0) {
      cp9b->Lmarg_jmin = i0 - 1; cp9b->Lmarg_jmax = j0;
    } else {
      int lmin = INT_MAX, lmax = INT_MIN;
      if(cp9b->ep1 != 0 && cp9b->pn_min_m[cp9b->ep1] >= 0) lmin = ESL_MIN(lmin, cp9b->pn_min_m[cp9b->ep1]);
      if(cp9b->ep1 != 0 && cp9b->pn_min_i[cp9b->ep1] >= 0) lmin = ESL_MIN(lmin, cp9b->pn_min_i[cp9b->ep1]);
      if(cp9b->ep1 != 0 && cp9b->pn_min_d[cp9b->ep1] >= 0) lmin = ESL_MIN(lmin, cp9b->pn_min_d[cp9b->ep1] - 1);
      if(cp9b->ep2 != 0 && cp9b->pn_min_m[cp9b->ep2] >= 0) lmin = ESL_MIN(lmin, cp9b->pn_min_m[cp9b->ep2]);
      if(cp9b->ep2 != 0 && cp9b->pn_min_i[cp9b->ep2] >= 0) lmin = ESL_MIN(lmin, cp9b->pn_min_i[cp9b->ep2]);
      if(cp9b->ep2 != 0 && cp9b->pn_min_d[cp9b->ep2] >= 0) lmin = ESL_MIN(lmin, cp9b->pn_min_d[cp9b->ep2] - 1);
      if(lmin == INT_MAX || cp9b->ep1 == 0 || cp9b->ep2 == 0) lmin = i0 - 1;
      cp9b->Lmarg_jmin = ESL_MAX(i0 - 1, ESL_MIN(j0, lmin));

      if(cp9b->ep1 != 0 && cp9b->pn_max_m[cp9b->ep1] >= 0) lmax = ESL_MAX(lmax, cp9b->pn_max_m[cp9b->ep1]);
      if(cp9b->ep1 != 0 && cp9b->pn_max_i[cp9b->ep1] >= 0) lmax = ESL_MAX(lmax, cp9b->pn_max_i[cp9b->ep1]);
      if(cp9b->ep1 != 0 && cp9b->pn_max_d[cp9b->ep1] >= 0) lmax = ESL_MAX(lmax, cp9b->pn_max_d[cp9b->ep1] - 1);
      if(cp9b->ep2 != 0 && cp9b->pn_max_m[cp9b->ep2] >= 0) lmax = ESL_MAX(lmax, cp9b->pn_max_m[cp9b->ep2]);
      if(cp9b->ep2 != 0 && cp9b->pn_max_i[cp9b->ep2] >= 0) lmax = ESL_MAX(lmax, cp9b->pn_max_i[cp9b->ep2]);
      if(cp9b->ep2 != 0 && cp9b->pn_max_d[cp9b->ep2] >= 0) lmax = ESL_MAX(lmax, cp9b->pn_max_d[cp9b->ep2] - 1);
      if(lmax == INT_MIN || cp9b->ep1 == 0 || cp9b->ep2 == 0) lmax = j0;
      cp9b->Lmarg_jmax = ESL_MAX(i0 - 1, ESL_MIN(j0, lmax));
    }

    if((status = cp9_MarginalCandidatesFromStartEndPositions(cm, cp9b, pass_idx, errbuf)) != eslOK) goto ERROR;
  } else {
    esl_vec_ISet(cp9b->Jvalid, cm->M + 1, TRUE);
    esl_vec_ISet(cp9b->Lvalid, cm->M + 1, FALSE);
    esl_vec_ISet(cp9b->Rvalid, cm->M + 1, FALSE);
    esl_vec_ISet(cp9b->Tvalid, cm->M + 1, FALSE);
  }

  /* HMM bands -> CM ij bands */
  if(do_old_hmm2ij) {
    if((status = cp9_HMM2ijBands_OLD(cm, errbuf, cm->cp9b, cm->cp9map, i0, j0, TRUE, debug_level)) != eslOK) goto ERROR;
  } else {
    if((status = cp9_HMM2ijBands(cm, errbuf, cp9, cm->cp9b, cm->cp9map, i0, j0, TRUE, do_trunc, debug_level)) != eslOK) goto ERROR;
  }

  /* CM ij bands -> CM d bands */
  if((status = cp9_GrowHDBands(cp9b, errbuf)) != eslOK) goto ERROR;
  ij2d_bands(cm, cp9b, do_trunc, debug_level);

#if eslDEBUGLEVEL >= 1
  if((status = cp9_ValidateBands(cm, errbuf, cp9b, i0, j0, do_trunc)) != eslOK) goto ERROR;
  ESL_DPRINTF1(("#DEBUG: p7banded_post_to_cp9bands bands validated.\n"));
#endif
  if(debug_level > 0) debug_print_ij_bands(cm);

  free(pocc);
  return eslOK;

 ERROR:
  if(pocc) free(pocc);
  return status;
}


/* Function: p7banded_post_to_pn_bands
 * Date    : EPN, 2026-03-26
 *
 * Purpose:  Phase 1 of p7banded_post_to_cp9bands: sweep the banded Forward/
 *           Backward matrices and fill per-node pn_min/max arrays.
 *
 *           Extracted so that the posterior sweep can be run at F5 time
 *           (while gxfb/gxbb are valid) and the results stored per envelope,
 *           decoupling them from the later CYK dispatch where gxfb/gxbb may
 *           already point to a different window.
 *
 * Args:     gxfb      - p7 banded glocal Forward matrix (log-space)
 *           gxbb      - p7 banded glocal Backward matrix (log-space)
 *           fwdsc     - Forward score in nats
 *           bnd       - band structure for gxfb/gxbb
 *           ws        - absolute start of the window (1-indexed)
 *           M         - number of HMM nodes (= length of pn_* arrays - 1)
 *           i0        - first absolute position of the envelope
 *           j0        - final absolute position of the envelope
 *           thresh    - posterior probability threshold
 *           L         - envelope length (j0 - i0 + 1)
 *           pn_min_m  - [0..M] pre-allocated output: min match position
 *           pn_max_m  - [0..M] pre-allocated output: max match position
 *           pn_min_i  - [0..M] pre-allocated output: min insert position
 *           pn_max_i  - [0..M] pre-allocated output: max insert position
 *           pn_min_d  - [0..M] pre-allocated output: min delete position
 *           pn_max_d  - [0..M] pre-allocated output: max delete position
 *           pocc      - [0..M] pre-allocated OR NULL: per-node match-posterior
 *                       occupancy sum (pocc[k] = Σ_i exp(fM+bM-fwdsc)).  Used
 *                       downstream by p7pn_bands_to_cp9cm_bands together with
 *                       cp9b->thresh1/thresh2 to compute sp1/sp2/ep1/ep2
 *                       (mimics the cp9 path's cp9_PredictStartAndEndPositions).
 *                       Pass NULL if the caller does not want it.
 *
 *           Positions are stored in 1..L (envelope-relative) coordinates.
 *           Unset nodes have pn_min[k] = pn_max[k] = -1.
 *           Node 0 match band is always set to 0 (begin state special case).
 *
 * Returns:  eslOK on success.
 */
int
p7banded_post_to_pn_bands(P7_GMXB *gxfb, P7_GMXB *gxbb, float fwdsc,
                           P7_GBANDS *bnd, int ws, int M, int i0, int j0,
                           float thresh, int L,
                           int *pn_min_m, int *pn_max_m,
                           int *pn_min_i, int *pn_max_i,
                           int *pn_min_d, int *pn_max_d,
                           float *pocc,
                           int do_pnmono, int do_pnmono_print)
{
  int          g, i, k;
  int          ia, ib;
  int          kac, kbc;
  int         *bnd_ip;
  int         *bnd_kp;
  float const *fwd_dp;
  float const *bck_dp;
  float        log_thresh;
  float        fM, bM, fI, bI;
  int          nk;
  int          first_kac, first_kbc;
  int          abs_i, ep;

  log_thresh = logf(thresh);

  for(k = 0; k <= M; k++) {
    pn_min_m[k] = L + 2;
    pn_max_m[k] = -1;
    pn_min_i[k] = L + 2;
    pn_max_i[k] = -1;
    pn_min_d[k] = L + 2;
    pn_max_d[k] = -1;
  }
  if(pocc != NULL) { for(k = 0; k <= M; k++) pocc[k] = 0.0f; }

  bnd_ip    = bnd->imem;
  bnd_kp    = bnd->kmem;
  fwd_dp    = gxfb->dp;
  bck_dp    = gxbb->dp;
  first_kac = bnd->kmem[0];
  first_kbc = bnd->kmem[1];

  for(g = 0; g < bnd->nseg; g++) {
    ia = *bnd_ip++;
    ib = *bnd_ip++;
    for(i = ia; i <= ib; i++) {
      kac = *bnd_kp++;
      kbc = *bnd_kp++;
      nk  = kbc - kac + 1;
      abs_i = ws + i - 1;
      ep    = abs_i - i0 + 1;
      if(abs_i >= i0 && abs_i <= j0) {
        for(k = kac; k <= kbc; k++) {
          fM = fwd_dp[(k-kac)*p7G_NSCELLS + p7G_M];
          bM = bck_dp[(k-kac)*p7G_NSCELLS + p7G_M];
          if(fM + bM - fwdsc > log_thresh) {
            if(ep < pn_min_m[k]) pn_min_m[k] = ep;
            if(ep > pn_max_m[k]) pn_max_m[k] = ep;
            if(k > 0) {
              if(ep < pn_min_d[k]) pn_min_d[k] = ep;
              if(ep > pn_max_d[k]) pn_max_d[k] = ep;
            }
          }
          if(pocc != NULL) pocc[k] += expf(fM + bM - fwdsc);
          if(k < M) {
            fI = fwd_dp[(k-kac)*p7G_NSCELLS + p7G_I];
            bI = bck_dp[(k-kac)*p7G_NSCELLS + p7G_I];
            if(fI + bI - fwdsc > log_thresh) {
              if(ep < pn_min_i[k]) pn_min_i[k] = ep;
              if(ep > pn_max_i[k]) pn_max_i[k] = ep;
            }
          }
        }
      }
      fwd_dp += nk * p7G_NSCELLS;
      bck_dp += nk * p7G_NSCELLS;
    }
  }

  /* i=0 special case: begin state M_0 always at position 0;
   * D_k for k in first band row may be active at position 0. */
  pn_min_m[0] = pn_max_m[0] = 0;
  for(k = first_kac; k <= first_kbc; k++) {
    if(k > 0) {
      if(0 < pn_min_d[k]) pn_min_d[k] = 0;
      if(0 > pn_max_d[k]) pn_max_d[k] = 0;
    }
  }

  /* I_0 and I_M (5'/3' flanking CP9 inserts) have no p7 counterparts
   * (p7 uses N and C states), so the posterior sweep never populates
   * their bands. Derive tight bounds from the adjacent match bands:
   * I_0 emits before M_1, I_M emits after M_M. Without this,
   * cp9_HMM2ijBands pins ROOT_IL/IR to envelope edges and the CM
   * pays mis-match penalty for flanking residues.
   */
  if (pn_max_m[1] >= 1) {
    pn_min_i[0] = 1;
    pn_max_i[0] = pn_max_m[1] - 1;
  }
  if (pn_min_m[M] >= 1 && pn_min_m[M] <= L) {
    pn_min_i[M] = pn_min_m[M] + 1;
    pn_max_i[M] = L;
  }

  /* Convert unset entries to -1 sentinel */
  for(k = 0; k <= M; k++) {
    if(pn_min_m[k] > pn_max_m[k]) pn_min_m[k] = pn_max_m[k] = -1;
    if(pn_min_i[k] > pn_max_i[k]) pn_min_i[k] = pn_max_i[k] = -1;
    if(pn_min_d[k] > pn_max_d[k]) pn_min_d[k] = pn_max_d[k] = -1;
  }
  pn_min_d[0] = pn_max_d[0] = -1; /* D_0 does not exist */

  if(do_pnmono) pn_match_bands_enforce_monotone(pn_min_m, pn_max_m, M, L, do_pnmono_print, "post_thresh");

  return eslOK;
}


/* Function: p7banded_post_to_pn_bands_tau
 * Date    : EPN*, Sat Apr  5 2026
 *
 * Purpose:  Cumulative-tau variant of p7banded_post_to_pn_bands().
 *           Instead of a hard posterior probability threshold, trims
 *           per-node bands using a cumulative probability fraction
 *           <tau>.  For each node k, the band is set so that at most
 *           tau/2 of the total posterior mass is excluded from each
 *           side (left and right).
 *
 *           Three passes over the banded F/B matrices:
 *             Pass 1 (forward): compute total_m[k] and total_i[k]
 *                    via log-sum-exp, and save per-row metadata
 *                    (dp pointer offsets, kac/kbc, abs_i) for reverse.
 *             Pass 2 (forward): left-trim — accumulate per-node
 *                    probability from left, set pn_min when cumsum
 *                    reaches tau/2.
 *             Pass 3 (backward): right-trim — walk rows in reverse
 *                    using saved pointers, set pn_max when cumsum
 *                    from right reaches tau/2.
 *
 * Args:     gxfb      - p7 banded glocal Forward matrix (log-space)
 *           gxbb      - p7 banded glocal Backward matrix (log-space)
 *           fwdsc     - Forward score in nats
 *           bnd       - band structure for gxfb/gxbb
 *           ws        - absolute start of the window (1-indexed)
 *           M         - number of HMM nodes
 *           i0        - first absolute position of the envelope
 *           j0        - final absolute position of the envelope
 *           tau       - cumulative probability fraction to trim
 *           L         - envelope length (j0 - i0 + 1)
 *           pn_min_m  - [0..M] pre-allocated output: min match position
 *           pn_max_m  - [0..M] pre-allocated output: max match position
 *           pn_min_i  - [0..M] pre-allocated output: min insert position
 *           pn_max_i  - [0..M] pre-allocated output: max insert position
 *           pn_min_d  - [0..M] pre-allocated output: min delete position
 *           pn_max_d  - [0..M] pre-allocated output: max delete position
 *
 *           Positions are stored in 1..L (envelope-relative) coordinates.
 *           Unset nodes have pn_min[k] = pn_max[k] = -1.
 *           Node 0 match band is always set to 0 (begin state special case).
 *
 * Returns:  eslOK on success.
 */
int
p7banded_post_to_pn_bands_tau(P7_GMXB *gxfb, P7_GMXB *gxbb, float fwdsc,
                               P7_GBANDS *bnd, int ws, int M, int i0, int j0,
                               float tau, int L,
                               int *pn_min_m, int *pn_max_m,
                               int *pn_min_i, int *pn_max_i,
                               int *pn_min_d, int *pn_max_d,
                               int do_pnmono, int do_pnmono_print)
{
  int          status;
  int          g, i, k;
  int          ia, ib;
  int          kac, kbc;
  int         *bnd_ip;
  int         *bnd_kp;
  float const *fwd_dp;
  float const *bck_dp;
  float        fM, bM, fI, bI;
  float        log_post, prob;
  int          nk;
  int          first_kac, first_kbc;
  int          abs_i, ep;
  float        half_tau = tau / 2.0f;

  /* per-node totals (log space) */
  float *total_m = NULL;
  float *total_i = NULL;
  /* per-node cumulative sums (probability space) for left and right trim */
  float *left_cumsum_m  = NULL;
  float *left_cumsum_i  = NULL;
  float *right_cumsum_m = NULL;
  float *right_cumsum_i = NULL;

  /* Per-row saved data for reverse pass.
   * We count total rows in the banded matrix first. */
  int    nrows = 0;
  int   *row_kac   = NULL;
  int   *row_kbc   = NULL;
  int   *row_ep    = NULL;
  int   *row_abs_i = NULL;
  int   *row_fwd_off = NULL;  /* offset into gxfb->dp for this row's first cell */
  int   *row_bck_off = NULL;  /* offset into gxbb->dp for this row's first cell */

  /* Count total rows */
  bnd_ip = bnd->imem;
  for(g = 0; g < bnd->nseg; g++) {
    ia = *bnd_ip++;
    ib = *bnd_ip++;
    nrows += (ib - ia + 1);
  }

  /* Allocations */
  ESL_ALLOC(total_m,       sizeof(float) * (M+1));
  ESL_ALLOC(total_i,       sizeof(float) * (M+1));
  ESL_ALLOC(left_cumsum_m, sizeof(float) * (M+1));
  ESL_ALLOC(left_cumsum_i, sizeof(float) * (M+1));
  ESL_ALLOC(right_cumsum_m,sizeof(float) * (M+1));
  ESL_ALLOC(right_cumsum_i,sizeof(float) * (M+1));
  ESL_ALLOC(row_kac,       sizeof(int) * nrows);
  ESL_ALLOC(row_kbc,       sizeof(int) * nrows);
  ESL_ALLOC(row_ep,        sizeof(int) * nrows);
  ESL_ALLOC(row_abs_i,     sizeof(int) * nrows);
  ESL_ALLOC(row_fwd_off,   sizeof(int) * nrows);
  ESL_ALLOC(row_bck_off,   sizeof(int) * nrows);

  for(k = 0; k <= M; k++) {
    total_m[k]       = -eslINFINITY;
    total_i[k]       = -eslINFINITY;
    left_cumsum_m[k] = 0.0f;
    left_cumsum_i[k] = 0.0f;
    right_cumsum_m[k]= 0.0f;
    right_cumsum_i[k]= 0.0f;
  }

  for(k = 0; k <= M; k++) {
    pn_min_m[k] = L + 2;
    pn_max_m[k] = -1;
    pn_min_i[k] = L + 2;
    pn_max_i[k] = -1;
    pn_min_d[k] = L + 2;
    pn_max_d[k] = -1;
  }

  first_kac = bnd->kmem[0];
  first_kbc = bnd->kmem[1];

  /* ---- Pass 1: compute total_m[k], total_i[k] via log-sum-exp;
   *              save per-row metadata for reverse pass ---- */
  bnd_ip = bnd->imem;
  bnd_kp = bnd->kmem;
  fwd_dp = gxfb->dp;
  bck_dp = gxbb->dp;
  int r = 0;

  for(g = 0; g < bnd->nseg; g++) {
    ia = *bnd_ip++;
    ib = *bnd_ip++;
    for(i = ia; i <= ib; i++) {
      kac = *bnd_kp++;
      kbc = *bnd_kp++;
      nk  = kbc - kac + 1;
      abs_i = ws + i - 1;
      ep    = abs_i - i0 + 1;

      /* save row metadata */
      row_kac[r]     = kac;
      row_kbc[r]     = kbc;
      row_ep[r]      = ep;
      row_abs_i[r]   = abs_i;
      row_fwd_off[r] = (int)(fwd_dp - gxfb->dp);
      row_bck_off[r] = (int)(bck_dp - gxbb->dp);

      if(abs_i >= i0 && abs_i <= j0) {
        for(k = kac; k <= kbc; k++) {
          fM = fwd_dp[(k-kac)*p7G_NSCELLS + p7G_M];
          bM = bck_dp[(k-kac)*p7G_NSCELLS + p7G_M];
          log_post = fM + bM - fwdsc;
          total_m[k] = p7_FLogsum(total_m[k], log_post);
          if(k > 0) {
            /* D states share the match posterior for banding purposes */
          }
          if(k < M) {
            fI = fwd_dp[(k-kac)*p7G_NSCELLS + p7G_I];
            bI = bck_dp[(k-kac)*p7G_NSCELLS + p7G_I];
            log_post = fI + bI - fwdsc;
            total_i[k] = p7_FLogsum(total_i[k], log_post);
          }
        }
      }
      fwd_dp += nk * p7G_NSCELLS;
      bck_dp += nk * p7G_NSCELLS;
      r++;
    }
  }

  /* ---- Pass 2 (forward): left-trim ---- */
  for(r = 0; r < nrows; r++) {
    kac   = row_kac[r];
    kbc   = row_kbc[r];
    abs_i = row_abs_i[r];
    ep    = row_ep[r];
    nk    = kbc - kac + 1;
    fwd_dp = gxfb->dp + row_fwd_off[r];
    bck_dp = gxbb->dp + row_bck_off[r];

    if(abs_i >= i0 && abs_i <= j0) {
      for(k = kac; k <= kbc; k++) {
        /* match */
        fM = fwd_dp[(k-kac)*p7G_NSCELLS + p7G_M];
        bM = bck_dp[(k-kac)*p7G_NSCELLS + p7G_M];
        log_post = fM + bM - fwdsc;
        if(total_m[k] > -eslINFINITY) {
          prob = expf(log_post - total_m[k]);
          left_cumsum_m[k] += prob;
          if(left_cumsum_m[k] >= half_tau && pn_min_m[k] > L + 1) {
            pn_min_m[k] = ep;
          }
          /* delete bands track match positions */
          if(k > 0 && left_cumsum_m[k] >= half_tau && pn_min_d[k] > L + 1) {
            pn_min_d[k] = ep;
          }
        }
        /* insert */
        if(k < M) {
          fI = fwd_dp[(k-kac)*p7G_NSCELLS + p7G_I];
          bI = bck_dp[(k-kac)*p7G_NSCELLS + p7G_I];
          log_post = fI + bI - fwdsc;
          if(total_i[k] > -eslINFINITY) {
            prob = expf(log_post - total_i[k]);
            left_cumsum_i[k] += prob;
            if(left_cumsum_i[k] >= half_tau && pn_min_i[k] > L + 1) {
              pn_min_i[k] = ep;
            }
          }
        }
      }
    }
  }

  /* ---- Pass 3 (backward): right-trim ---- */
  for(r = nrows - 1; r >= 0; r--) {
    kac   = row_kac[r];
    kbc   = row_kbc[r];
    abs_i = row_abs_i[r];
    ep    = row_ep[r];
    nk    = kbc - kac + 1;
    fwd_dp = gxfb->dp + row_fwd_off[r];
    bck_dp = gxbb->dp + row_bck_off[r];

    if(abs_i >= i0 && abs_i <= j0) {
      for(k = kac; k <= kbc; k++) {
        /* match */
        fM = fwd_dp[(k-kac)*p7G_NSCELLS + p7G_M];
        bM = bck_dp[(k-kac)*p7G_NSCELLS + p7G_M];
        log_post = fM + bM - fwdsc;
        if(total_m[k] > -eslINFINITY) {
          prob = expf(log_post - total_m[k]);
          right_cumsum_m[k] += prob;
          if(right_cumsum_m[k] >= half_tau && pn_max_m[k] == -1) {
            pn_max_m[k] = ep;
          }
          if(k > 0 && right_cumsum_m[k] >= half_tau && pn_max_d[k] == -1) {
            pn_max_d[k] = ep;
          }
        }
        /* insert */
        if(k < M) {
          fI = fwd_dp[(k-kac)*p7G_NSCELLS + p7G_I];
          bI = bck_dp[(k-kac)*p7G_NSCELLS + p7G_I];
          log_post = fI + bI - fwdsc;
          if(total_i[k] > -eslINFINITY) {
            prob = expf(log_post - total_i[k]);
            right_cumsum_i[k] += prob;
            if(right_cumsum_i[k] >= half_tau && pn_max_i[k] == -1) {
              pn_max_i[k] = ep;
            }
          }
        }
      }
    }
  }

  /* i=0 special case: begin state M_0 always at position 0;
   * D_k for k in first band row may be active at position 0. */
  pn_min_m[0] = pn_max_m[0] = 0;
  for(k = first_kac; k <= first_kbc; k++) {
    if(k > 0) {
      if(0 < pn_min_d[k]) pn_min_d[k] = 0;
      if(0 > pn_max_d[k]) pn_max_d[k] = 0;
    }
  }

  /* I_0 and I_M (5'/3' flanking CP9 inserts) have no p7 counterparts;
   * see p7banded_post_to_pn_bands() for full explanation. */
  if (pn_max_m[1] >= 1) {
    pn_min_i[0] = 1;
    pn_max_i[0] = pn_max_m[1] - 1;
  }
  if (pn_min_m[M] >= 1 && pn_min_m[M] <= L) {
    pn_min_i[M] = pn_min_m[M] + 1;
    pn_max_i[M] = L;
  }

  /* Convert unset entries to -1 sentinel */
  for(k = 0; k <= M; k++) {
    if(pn_min_m[k] > pn_max_m[k]) pn_min_m[k] = pn_max_m[k] = -1;
    if(pn_min_i[k] > pn_max_i[k]) pn_min_i[k] = pn_max_i[k] = -1;
    if(pn_min_d[k] > pn_max_d[k]) pn_min_d[k] = pn_max_d[k] = -1;
  }
  pn_min_d[0] = pn_max_d[0] = -1; /* D_0 does not exist */

  if(do_pnmono) pn_match_bands_enforce_monotone(pn_min_m, pn_max_m, M, L, do_pnmono_print, "post_tau");

  free(total_m);
  free(total_i);
  free(left_cumsum_m);
  free(left_cumsum_i);
  free(right_cumsum_m);
  free(right_cumsum_i);
  free(row_kac);
  free(row_kbc);
  free(row_ep);
  free(row_abs_i);
  free(row_fwd_off);
  free(row_bck_off);
  return eslOK;

 ERROR:
  if(total_m)        free(total_m);
  if(total_i)        free(total_i);
  if(left_cumsum_m)  free(left_cumsum_m);
  if(left_cumsum_i)  free(left_cumsum_i);
  if(right_cumsum_m) free(right_cumsum_m);
  if(right_cumsum_i) free(right_cumsum_i);
  if(row_kac)        free(row_kac);
  if(row_kbc)        free(row_kbc);
  if(row_ep)         free(row_ep);
  if(row_abs_i)      free(row_abs_i);
  if(row_fwd_off)    free(row_fwd_off);
  if(row_bck_off)    free(row_bck_off);
  return status;
}


/* Function: p7pn_bands_to_cp9cm_bands
 * Date    : EPN, 2026-03-26
 *
 * Purpose:  Phase 2 of p7banded_post_to_cp9bands: convert pre-filled pn_min/max
 *           arrays (from p7banded_post_to_pn_bands) into CP9 HMM bands and then
 *           CM ij/d bands.
 *
 *           Copies pn_min/max into cp9b->pn_min/max, applies the i0..j0 coordinate
 *           shift, handles truncation candidates, then calls cp9_HMM2ijBands,
 *           cp9_GrowHDBands, and ij2d_bands.
 *
 * Args:     cm          - covariance model
 *           errbuf      - error buffer
 *           pn_min_m..pn_max_d - pre-filled pn arrays [0..M] (1..L envelope-relative)
 *           pocc        - [0..M] per-node match-posterior occupancy OR NULL.
 *                         When non-NULL, enables cp9-style computation of sp1/sp2/ep1/ep2
 *                         via cp9b->thresh1/thresh2 comparisons (mimics
 *                         hmmband.c:cp9_PredictStartAndEndPositions). When NULL, falls back
 *                         to the legacy collapse sp1=sp2=min_k(pn_min_m!=-1), ep1=ep2=max_k.
 *           cp9b        - pre-allocated CP9 bands, filled here
 *           i0          - first absolute position of the envelope
 *           j0          - final absolute position of the envelope
 *           L           - envelope length (j0 - i0 + 1)
 *           pass_idx    - pipeline pass index
 *           debug_level - verbosity
 *
 * Returns:  eslOK on success, error code on failure.
 */
int
p7pn_bands_to_cp9cm_bands(CM_t *cm, char *errbuf,
                           int *pn_min_m, int *pn_max_m,
                           int *pn_min_i, int *pn_max_i,
                           int *pn_min_d, int *pn_max_d,
                           const float *pocc,
                           CP9Bands_t *cp9b, int i0, int j0, int L,
                           int pass_idx, int debug_level)
{
  int    status;
  int    k;
  int    M            = cp9b->hmm_M;
  int    do_old_hmm2ij;
  int    do_trunc;
  CP9_t *cp9;

  do_old_hmm2ij = ((cm->align_opts & CM_ALIGN_HMM2IJOLD) || (cm->search_opts & CM_SEARCH_HMM2IJOLD)) ? TRUE : FALSE;
  do_trunc      = cm_pli_PassAllowsTruncation(pass_idx);
  switch(pass_idx) {
    case PLI_PASS_5P_ONLY_FORCE:                cp9 = cm->Rcp9; break;
    case PLI_PASS_3P_ONLY_FORCE:                cp9 = cm->Lcp9; break;
    case PLI_PASS_5P_AND_3P_FORCE:
    case PLI_PASS_5P_AND_3P_ANY:                cp9 = cm->Tcp9; break;
    default:                                    cp9 = cm->cp9;  break;
  }

  /* Copy pn_min/max into cp9b */
  for(k = 0; k <= M; k++) {
    cp9b->pn_min_m[k] = pn_min_m[k];
    cp9b->pn_max_m[k] = pn_max_m[k];
    cp9b->pn_min_i[k] = pn_min_i[k];
    cp9b->pn_max_i[k] = pn_max_i[k];
    cp9b->pn_min_d[k] = pn_min_d[k];
    cp9b->pn_max_d[k] = pn_max_d[k];
  }
  cp9b->tau = cm->tau;

  /* Shift HMM bands from 1..L to i0..j0 coordinate system if needed */
  if(i0 != 1) {
    int offset = i0 - 1;
    for(k = 0; k <= M; k++) {
      if(cp9b->pn_min_m[k] != -1) { cp9b->pn_min_m[k] += offset; cp9b->pn_max_m[k] += offset; }
      if(cp9b->pn_min_i[k] != -1) { cp9b->pn_min_i[k] += offset; cp9b->pn_max_i[k] += offset; }
      if(cp9b->pn_min_d[k] != -1) { cp9b->pn_min_d[k] += offset; cp9b->pn_max_d[k] += offset; }
    }
  }

  /* Set truncation candidate valid arrays.
   *
   * Two modes:
   * - pocc != NULL: cp9-style — sp1/sp2 are leftmost k where pocc[k] crosses
   *   thresh1/thresh2; ep1/ep2 symmetric. Rmarg/Lmarg take union over sp1,sp2
   *   and ep1,ep2 (mirrors hmmband.c:cp9_PredictStartAndEndPositions).
   * - pocc == NULL: legacy fallback — collapses sp1=sp2=leftmost_k_with_coverage,
   *   ep1=ep2=rightmost_k. Used by paths that don't cache pocc (e.g. the
   *   --p7post_tau variant and the final-stage band-shift call at
   *   cm_pipeline.c:5807).
   *
   */
  if(do_trunc) {
    if(pocc != NULL) {
      /* Part 1: sp1/sp2 by threshold sweep (cp9_PredictStartAndEndPositions analogue) */
      cp9b->sp1 = cp9b->sp2 = -1;
      for(k = 1; k <= M; k++) {
        if(cp9b->pn_min_m[k] == -1 && cp9b->pn_min_i[k] == -1 && cp9b->pn_min_d[k] == -1) continue;
        if(cp9b->sp1 == -1 && pocc[k] > cp9b->thresh1) cp9b->sp1 = k;
        if(cp9b->sp2 == -1 && pocc[k] > cp9b->thresh2) cp9b->sp2 = k;
        if(cp9b->sp1 != -1 && cp9b->sp2 != -1) break;
      }
      if(cp9b->sp1 == -1) cp9b->sp1 = M + 1;
      if(cp9b->sp2 == -1) cp9b->sp2 = M + 1;

      /* Part 2: ep1/ep2 symmetric sweep from the right */
      if(cp9b->sp1 == M + 1 && cp9b->sp2 == M + 1) {
        cp9b->ep1 = cp9b->ep2 = 0;
      } else {
        cp9b->ep1 = cp9b->ep2 = -1;
        for(k = M; k >= 1; k--) {
          if(cp9b->pn_min_m[k] == -1 && cp9b->pn_min_i[k] == -1 && cp9b->pn_min_d[k] == -1) continue;
          if(cp9b->ep1 == -1 && pocc[k] > cp9b->thresh1) cp9b->ep1 = k;
          if(cp9b->ep2 == -1 && pocc[k] > cp9b->thresh2) cp9b->ep2 = k;
          if(cp9b->ep1 != -1 && cp9b->ep2 != -1) break;
        }
        if(cp9b->ep1 == -1) cp9b->ep1 = 0;
        if(cp9b->ep2 == -1) cp9b->ep2 = 0;
      }

      /* Part 3: Rmarg_i{min,max} from UNION of sp1 and sp2 bands (hmmband.c:4857-4884) */
      if(cp9b->sp1 == M + 1) {
        cp9b->Rmarg_imin = i0;
        cp9b->Rmarg_imax = j0;
      } else {
        int rmin = INT_MAX, rmax = INT_MIN;
        if(cp9b->sp1 != (M+1) && cp9b->pn_min_m[cp9b->sp1] >= 0) rmin = ESL_MIN(rmin, cp9b->pn_min_m[cp9b->sp1]);
        if(cp9b->sp1 != (M+1) && cp9b->pn_min_i[cp9b->sp1] >= 0) rmin = ESL_MIN(rmin, cp9b->pn_min_i[cp9b->sp1]);
        if(cp9b->sp1 != (M+1) && cp9b->pn_min_d[cp9b->sp1] >= 0) rmin = ESL_MIN(rmin, cp9b->pn_min_d[cp9b->sp1]);
        if(cp9b->sp2 != (M+1) && cp9b->pn_min_m[cp9b->sp2] >= 0) rmin = ESL_MIN(rmin, cp9b->pn_min_m[cp9b->sp2]);
        if(cp9b->sp2 != (M+1) && cp9b->pn_min_i[cp9b->sp2] >= 0) rmin = ESL_MIN(rmin, cp9b->pn_min_i[cp9b->sp2]);
        if(cp9b->sp2 != (M+1) && cp9b->pn_min_d[cp9b->sp2] >= 0) rmin = ESL_MIN(rmin, cp9b->pn_min_d[cp9b->sp2]);
        if(rmin == INT_MAX || cp9b->sp1 == (M+1) || cp9b->sp2 == (M+1)) rmin = i0;
        cp9b->Rmarg_imin = ESL_MAX(i0, ESL_MIN(j0 + 1, rmin));

        if(cp9b->sp1 != (M+1) && cp9b->pn_max_m[cp9b->sp1] >= 0) rmax = ESL_MAX(rmax, cp9b->pn_max_m[cp9b->sp1]);
        if(cp9b->sp1 != (M+1) && cp9b->pn_max_i[cp9b->sp1] >= 0) rmax = ESL_MAX(rmax, cp9b->pn_max_i[cp9b->sp1]);
        if(cp9b->sp1 != (M+1) && cp9b->pn_max_d[cp9b->sp1] >= 0) rmax = ESL_MAX(rmax, cp9b->pn_max_d[cp9b->sp1]);
        if(cp9b->sp2 != (M+1) && cp9b->pn_max_m[cp9b->sp2] >= 0) rmax = ESL_MAX(rmax, cp9b->pn_max_m[cp9b->sp2]);
        if(cp9b->sp2 != (M+1) && cp9b->pn_max_i[cp9b->sp2] >= 0) rmax = ESL_MAX(rmax, cp9b->pn_max_i[cp9b->sp2]);
        if(cp9b->sp2 != (M+1) && cp9b->pn_max_d[cp9b->sp2] >= 0) rmax = ESL_MAX(rmax, cp9b->pn_max_d[cp9b->sp2]);
        if(rmax == INT_MIN || cp9b->sp1 == (M+1) || cp9b->sp2 == (M+1)) rmax = j0 + 1;
        cp9b->Rmarg_imax = ESL_MAX(i0, ESL_MIN(j0 + 1, rmax));
      }

      /* Part 4: Lmarg_j{min,max} from UNION of ep1 and ep2 bands (hmmband.c:4886-4914) */
      if(cp9b->ep1 == 0) {
        cp9b->Lmarg_jmin = i0 - 1;
        cp9b->Lmarg_jmax = j0;
      } else {
        int lmin = INT_MAX, lmax = INT_MIN;
        if(cp9b->ep1 != 0 && cp9b->pn_min_m[cp9b->ep1] >= 0) lmin = ESL_MIN(lmin, cp9b->pn_min_m[cp9b->ep1]);
        if(cp9b->ep1 != 0 && cp9b->pn_min_i[cp9b->ep1] >= 0) lmin = ESL_MIN(lmin, cp9b->pn_min_i[cp9b->ep1]);
        if(cp9b->ep1 != 0 && cp9b->pn_min_d[cp9b->ep1] >= 0) lmin = ESL_MIN(lmin, cp9b->pn_min_d[cp9b->ep1] - 1);
        if(cp9b->ep2 != 0 && cp9b->pn_min_m[cp9b->ep2] >= 0) lmin = ESL_MIN(lmin, cp9b->pn_min_m[cp9b->ep2]);
        if(cp9b->ep2 != 0 && cp9b->pn_min_i[cp9b->ep2] >= 0) lmin = ESL_MIN(lmin, cp9b->pn_min_i[cp9b->ep2]);
        if(cp9b->ep2 != 0 && cp9b->pn_min_d[cp9b->ep2] >= 0) lmin = ESL_MIN(lmin, cp9b->pn_min_d[cp9b->ep2] - 1);
        if(lmin == INT_MAX || cp9b->ep1 == 0 || cp9b->ep2 == 0) lmin = i0 - 1;
        cp9b->Lmarg_jmin = ESL_MAX(i0 - 1, ESL_MIN(j0, lmin));

        if(cp9b->ep1 != 0 && cp9b->pn_max_m[cp9b->ep1] >= 0) lmax = ESL_MAX(lmax, cp9b->pn_max_m[cp9b->ep1]);
        if(cp9b->ep1 != 0 && cp9b->pn_max_i[cp9b->ep1] >= 0) lmax = ESL_MAX(lmax, cp9b->pn_max_i[cp9b->ep1]);
        if(cp9b->ep1 != 0 && cp9b->pn_max_d[cp9b->ep1] >= 0) lmax = ESL_MAX(lmax, cp9b->pn_max_d[cp9b->ep1] - 1);
        if(cp9b->ep2 != 0 && cp9b->pn_max_m[cp9b->ep2] >= 0) lmax = ESL_MAX(lmax, cp9b->pn_max_m[cp9b->ep2]);
        if(cp9b->ep2 != 0 && cp9b->pn_max_i[cp9b->ep2] >= 0) lmax = ESL_MAX(lmax, cp9b->pn_max_i[cp9b->ep2]);
        if(cp9b->ep2 != 0 && cp9b->pn_max_d[cp9b->ep2] >= 0) lmax = ESL_MAX(lmax, cp9b->pn_max_d[cp9b->ep2] - 1);
        if(lmax == INT_MIN || cp9b->ep1 == 0 || cp9b->ep2 == 0) lmax = j0;
        cp9b->Lmarg_jmax = ESL_MAX(i0 - 1, ESL_MIN(j0, lmax));
      }
    } else {
      /* Legacy fallback: extent-based collapse (sp1=sp2, ep1=ep2). */
      int sp = M + 1, ep = 0;
      for(k = 1; k <= M; k++) {
        if(cp9b->pn_min_m[k] != -1 && k < sp) sp = k;
        if(cp9b->pn_max_m[k] != -1 && k > ep) ep = k;
      }
      if(sp < 1)     sp = 1;
      if(sp > M + 1) sp = M + 1;
      if(ep < 0)     ep = 0;
      if(ep > M)     ep = M;
      cp9b->sp1 = cp9b->sp2 = sp;
      cp9b->ep1 = cp9b->ep2 = ep;

      if(cp9b->sp1 == M + 1) {
        cp9b->Rmarg_imin = i0;
        cp9b->Rmarg_imax = j0;
      } else {
        int rmarg_imin = INT_MAX, rmarg_imax = INT_MIN;
        if(cp9b->pn_min_m[sp] >= 0) rmarg_imin = ESL_MIN(rmarg_imin, cp9b->pn_min_m[sp]);
        if(cp9b->pn_min_i[sp] >= 0) rmarg_imin = ESL_MIN(rmarg_imin, cp9b->pn_min_i[sp]);
        if(cp9b->pn_min_d[sp] >= 0) rmarg_imin = ESL_MIN(rmarg_imin, cp9b->pn_min_d[sp]);
        if(rmarg_imin == INT_MAX)    rmarg_imin = i0;
        cp9b->Rmarg_imin = ESL_MAX(i0, ESL_MIN(j0 + 1, rmarg_imin));

        if(cp9b->pn_max_m[sp] >= 0) rmarg_imax = ESL_MAX(rmarg_imax, cp9b->pn_max_m[sp]);
        if(cp9b->pn_max_i[sp] >= 0) rmarg_imax = ESL_MAX(rmarg_imax, cp9b->pn_max_i[sp]);
        if(cp9b->pn_max_d[sp] >= 0) rmarg_imax = ESL_MAX(rmarg_imax, cp9b->pn_max_d[sp]);
        if(rmarg_imax == INT_MIN)    rmarg_imax = j0 + 1;
        cp9b->Rmarg_imax = ESL_MAX(i0, ESL_MIN(j0 + 1, rmarg_imax));
      }

      if(cp9b->ep1 == 0) {
        cp9b->Lmarg_jmin = i0 - 1;
        cp9b->Lmarg_jmax = j0;
      } else {
        int lmarg_jmin = INT_MAX, lmarg_jmax = INT_MIN;
        if(cp9b->pn_min_m[ep] >= 0) lmarg_jmin = ESL_MIN(lmarg_jmin, cp9b->pn_min_m[ep]);
        if(cp9b->pn_min_i[ep] >= 0) lmarg_jmin = ESL_MIN(lmarg_jmin, cp9b->pn_min_i[ep]);
        if(cp9b->pn_min_d[ep] >= 0) lmarg_jmin = ESL_MIN(lmarg_jmin, cp9b->pn_min_d[ep] - 1);
        if(lmarg_jmin == INT_MAX)    lmarg_jmin = i0 - 1;
        cp9b->Lmarg_jmin = ESL_MAX(i0 - 1, ESL_MIN(j0, lmarg_jmin));

        if(cp9b->pn_max_m[ep] >= 0) lmarg_jmax = ESL_MAX(lmarg_jmax, cp9b->pn_max_m[ep]);
        if(cp9b->pn_max_i[ep] >= 0) lmarg_jmax = ESL_MAX(lmarg_jmax, cp9b->pn_max_i[ep]);
        if(cp9b->pn_max_d[ep] >= 0) lmarg_jmax = ESL_MAX(lmarg_jmax, cp9b->pn_max_d[ep] - 1);
        if(lmarg_jmax == INT_MIN)    lmarg_jmax = j0;
        cp9b->Lmarg_jmax = ESL_MAX(i0 - 1, ESL_MIN(j0, lmarg_jmax));
      }
    }
    if((status = cp9_MarginalCandidatesFromStartEndPositions(cm, cp9b, pass_idx, errbuf)) != eslOK) return status;
  } else {
    esl_vec_ISet(cp9b->Jvalid, cm->M + 1, TRUE);
    esl_vec_ISet(cp9b->Lvalid, cm->M + 1, FALSE);
    esl_vec_ISet(cp9b->Rvalid, cm->M + 1, FALSE);
    esl_vec_ISet(cp9b->Tvalid, cm->M + 1, FALSE);
  }

  /* HMM bands -> CM ij bands */
  if(do_old_hmm2ij) {
    if((status = cp9_HMM2ijBands_OLD(cm, errbuf, cm->cp9b, cm->cp9map, i0, j0, TRUE, debug_level)) != eslOK) return status;
  } else {
    if((status = cp9_HMM2ijBands(cm, errbuf, cp9, cm->cp9b, cm->cp9map, i0, j0, TRUE, do_trunc, debug_level)) != eslOK) return status;
  }

  /* CM ij bands -> CM d bands */
  if((status = cp9_GrowHDBands(cp9b, errbuf)) != eslOK) return status;
  ij2d_bands(cm, cp9b, do_trunc, debug_level);

#if eslDEBUGLEVEL >= 1
  if((status = cp9_ValidateBands(cm, errbuf, cp9b, i0, j0, do_trunc)) != eslOK) return status;
  ESL_DPRINTF1(("#DEBUG: p7pn_bands_to_cp9cm_bands bands validated.\n"));
#endif
  if(debug_level > 0) debug_print_ij_bands(cm);

  return eslOK;
}


/* Function: cm_BandsFromParsetree()
 * Date    : EPN, Mon Apr  7 2026
 *
 * Purpose:  Derive CP9 HMM bands (pn_min_m/max_m/min_i/max_i/min_d/max_d)
 *           from a CYK parsetree, with a fixed half-width pad. Designed
 *           to be the CM analog of p7_pins2bands().
 *
 *           Walk the parsetree, recording the (i,j) positions where each
 *           CM node's match/insert/delete states emit. Map CM states to
 *           HMM consensus positions via cm->cp9map->cs2hn/cs2hs. Then fill
 *           the pn_min/max arrays with [pos-pad, pos+pad] clamped to
 *           the envelope, and do left/right sweeps to fill positions
 *           that don't have a direct contribution from the parsetree.
 *
 *           Coordinates: parsetree emitl/emitr are in 1..L (subseq local)
 *           coordinates. The output pn_min/max arrays use the same 1..L
 *           local coordinates. Caller can shift to envelope coordinates
 *           via p7pn_bands_to_cp9cm_bands.
 *
 * Args:     cm           - the covariance model (must have cp9map)
 *           tr           - parsetree from CYK alignment
 *           L            - length of subsequence aligned (sets the upper bound for pad clamp)
 *           pad          - band half-width on each side of each pinned position
 *           pn_min_m     - [0..M] OUTPUT: pre-allocated, will be filled
 *           pn_max_m     - [0..M] OUTPUT
 *           pn_min_i     - [0..M] OUTPUT
 *           pn_max_i     - [0..M] OUTPUT
 *           pn_min_d     - [0..M] OUTPUT
 *           pn_max_d     - [0..M] OUTPUT
 *
 * Returns:  eslOK on success.
 */
int
cm_BandsFromParsetree(CM_t *cm, Parsetree_t *tr, int L, int pad,
                      int *pn_min_m, int *pn_max_m,
                      int *pn_min_i, int *pn_max_i,
                      int *pn_min_d, int *pn_max_d)
{
  int  M = cm->cp9map->hmm_M;
  int  k, t, v, i, j;
  int  hn1, hs1, hn2, hs2;
  int  ka, kb;

  /* Initialize all bands to "unset" sentinel values */
  for(k = 0; k <= M; k++) {
    pn_min_m[k] = L + 2; pn_max_m[k] = -1;
    pn_min_i[k] = L + 2; pn_max_i[k] = -1;
    pn_min_d[k] = L + 2; pn_max_d[k] = -1;
  }

  /* Walk the parsetree, recording per-HMM-node emit positions */
  for(t = 0; t < tr->n; t++) {
    v = tr->state[t];
    i = tr->emitl[t];
    j = tr->emitr[t];

    /* Get HMM node(s) and state type(s) this CM state maps to */
    hn1 = cm->cp9map->cs2hn[v][0];
    hs1 = cm->cp9map->cs2hs[v][0];
    hn2 = cm->cp9map->cs2hn[v][1];
    hs2 = cm->cp9map->cs2hs[v][1];

    /* Update bands based on CM state type and which HMM node it maps to.
     * cp9map->cs2hs[v][n]: 0=MATCH, 1=INSERT, 2=DELETE
     * For MATP_MP: cs2hn[v][0] is the LEFT consensus position (HMM node), maps as MATCH
     *              cs2hn[v][1] is the RIGHT consensus position, also maps as MATCH
     *              The left position emits residue i, the right position emits residue j.
     * For MATP_ML: only the left HMM node, emits residue i (right is delete).
     * For MATP_MR: only the right HMM node, emits residue j (left is delete).
     * For MATL_ML/MATR_MR: only one HMM node, emits residue (i for ML, j for MR).
     * For inserts: emit one residue (i for IL, j for IR).
     * For deletes: emit nothing, but we still record the position around them.
     */
    if(hn1 >= 0 && hn1 <= M) {
      if(hs1 == 0) { /* MATCH */
        /* For MATP_MP, MATP_ML, MATL_ML: left HMM node emits residue i */
        /* For MATR_MR: only one HMM node, but it's stored in [v][0] and emits j */
        int pos = (cm->stid[v] == MATR_MR) ? j : i;
        if(pos < pn_min_m[hn1]) pn_min_m[hn1] = pos;
        if(pos > pn_max_m[hn1]) pn_max_m[hn1] = pos;
      }
      else if(hs1 == 1) { /* INSERT */
        /* IL emits i, IR emits j */
        int pos = (cm->sttype[v] == IR_st) ? j : i;
        if(pos < pn_min_i[hn1]) pn_min_i[hn1] = pos;
        if(pos > pn_max_i[hn1]) pn_max_i[hn1] = pos;
      }
      else if(hs1 == 2) { /* DELETE */
        /* No residue emitted; record the surrounding position.
         * Use i (the left bound) as a placeholder that gets refined by sweeps below. */
        if(i < pn_min_d[hn1]) pn_min_d[hn1] = i;
        if(i > pn_max_d[hn1]) pn_max_d[hn1] = i;
      }
    }
    /* Second HMM node (for MATP_MP, both left and right are matches; for MATP_D, both are deletes) */
    if(hn2 >= 0 && hn2 <= M) {
      if(hs2 == 0) {
        /* MATP_MP: hn2 is the right consensus position, emits residue j */
        int pos = j;
        if(pos < pn_min_m[hn2]) pn_min_m[hn2] = pos;
        if(pos > pn_max_m[hn2]) pn_max_m[hn2] = pos;
      }
      else if(hs2 == 1) {
        if(i < pn_min_i[hn2]) pn_min_i[hn2] = i;
        if(i > pn_max_i[hn2]) pn_max_i[hn2] = i;
      }
      else if(hs2 == 2) {
        if(j < pn_min_d[hn2]) pn_min_d[hn2] = j;
        if(j > pn_max_d[hn2]) pn_max_d[hn2] = j;
      }
    }
  }

  /* Apply pad and clamp to [1..L]. After this each filled band has half-width pad. */
  for(k = 0; k <= M; k++) {
    if(pn_min_m[k] != L + 2) {
      ka = pn_min_m[k] - pad; if(ka < 1) ka = 1;
      kb = pn_max_m[k] + pad; if(kb > L) kb = L;
      pn_min_m[k] = ka; pn_max_m[k] = kb;
    }
    if(pn_min_i[k] != L + 2) {
      ka = pn_min_i[k] - pad; if(ka < 1) ka = 1;
      kb = pn_max_i[k] + pad; if(kb > L) kb = L;
      pn_min_i[k] = ka; pn_max_i[k] = kb;
    }
    if(pn_min_d[k] != L + 2) {
      ka = pn_min_d[k] - pad; if(ka < 0) ka = 0;
      kb = pn_max_d[k] + pad; if(kb > L) kb = L;
      pn_min_d[k] = ka; pn_max_d[k] = kb;
    }
  }

  /* Sweeps: fill in unset bands by inheriting from neighbors.
   * Left-to-right: for each k with no band, set pn_*[k] to pn_*[k-1]'s range.
   * Right-to-left: similarly.
   * This handles HMM nodes that have no parsetree contribution
   * (e.g., delete-only nodes between matches). */
  /* Match bands left sweep */
  { int last_min = -1, last_max = -1;
    for(k = 1; k <= M; k++) {
      if(pn_min_m[k] != L + 2 && pn_max_m[k] != -1) {
        last_min = pn_min_m[k]; last_max = pn_max_m[k];
      } else if(last_min != -1) {
        pn_min_m[k] = last_min; pn_max_m[k] = last_max;
      }
    }
  }
  /* Match bands right sweep (catches early-k nodes that had no left neighbor) */
  { int last_min = -1, last_max = -1;
    for(k = M; k >= 1; k--) {
      if(pn_min_m[k] != L + 2 && pn_max_m[k] != -1
         && last_min == -1) {
        last_min = pn_min_m[k]; last_max = pn_max_m[k];
      } else if(pn_min_m[k] == L + 2 && last_min != -1) {
        pn_min_m[k] = last_min; pn_max_m[k] = last_max;
      }
    }
  }
  /* Insert bands: similarly */
  { int last_min = -1, last_max = -1;
    for(k = 0; k < M; k++) {
      if(pn_min_i[k] != L + 2 && pn_max_i[k] != -1) {
        last_min = pn_min_i[k]; last_max = pn_max_i[k];
      } else if(last_min != -1) {
        pn_min_i[k] = last_min; pn_max_i[k] = last_max;
      }
    }
  }
  { int last_min = -1, last_max = -1;
    for(k = M - 1; k >= 0; k--) {
      if(pn_min_i[k] != L + 2 && pn_max_i[k] != -1
         && last_min == -1) {
        last_min = pn_min_i[k]; last_max = pn_max_i[k];
      } else if(pn_min_i[k] == L + 2 && last_min != -1) {
        pn_min_i[k] = last_min; pn_max_i[k] = last_max;
      }
    }
  }
  /* Delete bands: similarly */
  { int last_min = -1, last_max = -1;
    for(k = 1; k <= M; k++) {
      if(pn_min_d[k] != L + 2 && pn_max_d[k] != -1) {
        last_min = pn_min_d[k]; last_max = pn_max_d[k];
      } else if(last_min != -1) {
        pn_min_d[k] = last_min; pn_max_d[k] = last_max;
      }
    }
  }
  { int last_min = -1, last_max = -1;
    for(k = M; k >= 1; k--) {
      if(pn_min_d[k] != L + 2 && pn_max_d[k] != -1
         && last_min == -1) {
        last_min = pn_min_d[k]; last_max = pn_max_d[k];
      } else if(pn_min_d[k] == L + 2 && last_min != -1) {
        pn_min_d[k] = last_min; pn_max_d[k] = last_max;
      }
    }
  }

  /* Node 0 special case: M_0 = begin state, must be at position 0 */
  pn_min_m[0] = 0; pn_max_m[0] = 0;

  /* CRITICAL: Ensure the bands include the full sequence range [1..L].
   * The downstream cp9_HMM2ijBands -> CM bands relies on these bounds:
   * - Some HMM node must be reachable at position 1 (start of seq)
   * - Some HMM node must be reachable at position L (end of seq)
   * - In particular, the LAST emitting node (typically near M) must have
   *   its match band include L, otherwise ROOT_S's j band won't include L
   *   and Inside DP will fail with "L is outside ROOT_S's j band".
   *
   * Walk left-to-right and ensure each node's band extends down to at most
   * pos 1 if it's near the start. Walk right-to-left and ensure each node's
   * band extends up to at least L if it's near the end. The simplest way:
   * widen the LAST set match band to include L, and the FIRST set match
   * band to include 1.
   */
  { int first_set = -1, last_set = -1;
    for(k = 1; k <= M; k++) {
      if(pn_min_m[k] != L + 2 && pn_max_m[k] != -1) {
        if(first_set == -1) first_set = k;
        last_set = k;
      }
    }
    if(first_set != -1) {
      /* Ensure first set node can be reached at position 1 */
      if(pn_min_m[first_set] > 1) pn_min_m[first_set] = 1;
    }
    if(last_set != -1) {
      /* Ensure last set node can be reached at position L */
      if(pn_max_m[last_set] < L) pn_max_m[last_set] = L;
    }
  }

  /* Convert any still-unset entries to -1 sentinel */
  for(k = 0; k <= M; k++) {
    if(pn_min_m[k] > pn_max_m[k] || pn_min_m[k] == L + 2) { pn_min_m[k] = -1; pn_max_m[k] = -1; }
    if(pn_min_i[k] > pn_max_i[k] || pn_min_i[k] == L + 2) { pn_min_i[k] = -1; pn_max_i[k] = -1; }
    if(pn_min_d[k] > pn_max_d[k] || pn_min_d[k] == L + 2) { pn_min_d[k] = -1; pn_max_d[k] = -1; }
  }
  pn_min_d[0] = pn_max_d[0] = -1; /* D_0 does not exist */

  return eslOK;
}


/* === ARCHIVED 2026-05-19: NOT HOOKED UP ===
 *
 * This function computes a per-state CYK pad array from CP9 HMM band
 * widths using a calibrated linear model:
 *   pad[v] = a[stt(v)] * mean_hmm_cells(v) + b[stt(v)]
 *
 * The cal_a[10] / cal_b[10] coefficients below were fit on
 * 26_0316 stage1.5d posterior-mass-mining compare.tsv (5 CMs, 6285
 * states with n_used >= 10 GT IO p99-pad samples). Recovering these
 * coefficients required real benchmarking work -- preserved here so
 * a future revival doesn't have to re-derive them.
 *
 * Empirical result (2026-05-12, benchmark-runs/cykperstate_walltime/):
 * perstate5 is 10.3x slower than uniform pad5 (sequence-weighted
 * aggregate). The long-tail per-state pads (LSU p99=90, max=160)
 * produce cell-count blowup that exceeds the cell savings from
 * tighter mean bands. Floor+cap variants (`f5_cap10`) still 1.39x
 * slower than pad5. p95 refit also failed (project memory
 * `project_cykperstate_p95_failed`).
 *
 * Decision (2026-05-19): kept in source as reference but not hooked
 * up. The `--cykbands-perstate`, `--cykperstate-maxpad`, `--cykpadfile`
 * CLI surface was removed along with the `cm->p7_cykbands_perstate`
 * and `cm->p7_cykperstate_maxpad` struct fields.
 *
 * To revive:
 *   1. Re-add a CLI flag (e.g., `--cykbands-perstate`) and a CM field
 *      to toggle this path.
 *   2. Call `cm_CYKPerstatePadCompute(cm, cp9b, additive_pad, maxpad)`
 *      in the cykbands caller in cm_alndata.c.
 *   3. Restore the per_state_pad parameter to the cykbands worker
 *      function (it was stripped on 2026-05-19 -- the conditional was
 *      a single line at step 2; see git history before the strip
 *      commit).
 *
 * === end ARCHIVED ===
 *
 * Function: cm_CYKPerstatePadCompute()
 * Date    : 2026-05-11
 *
 * Compute per-state CYK pad from current sequence's CP9 HMM band widths,
 * using a calibrated linear model pad = a[stt]*mean_hmm_cells + b[stt].
 *
 * Calibration: fit on 26_0316 stage1.5d posterior-mass-mining compare.tsv
 * (5 CMs, 6285 states with n_used>=10 ground-truth IO p99-pad samples).
 *
 * Returns malloc'd int[cm->M]; caller frees. Returns NULL on alloc failure.
 * additive_pad is added to each per-state pad (allows --cykpad N to act
 * as a floor when --cykbands-perstate is on). maxpad caps each per-state
 * pad at this value (0 = no cap).
 */
int *
cm_CYKPerstatePadCompute(CM_t *cm, CP9Bands_t *cp9b, int additive_pad, int maxpad)
{
  int  M = cm->M;
  int *pad_arr = malloc(sizeof(int) * M);
  if (pad_arr == NULL) return NULL;

  /* Calibration coefficients indexed by sttype constant:
   * D_st=0, MP_st=1, ML_st=2, MR_st=3, IL_st=4, IR_st=5, S_st=6, E_st=7, B_st=8, EL_st=9
   * E and EL have no calibration data (singletons / non-emitting); use conservative defaults.
   */
  static const double cal_a[10] = {
    0.0044,  /* D */   0.0140,  /* MP */  0.0225,  /* ML */  0.0088,  /* MR */
    0.0224,  /* IL */  0.0124,  /* IR */  0.0238,  /* S */
    0.0000,  /* E (no data) */
    0.0182,  /* B */
    0.0000   /* EL (no data) */
  };
  static const double cal_b[10] = {
    6.7158,  /* D */   0.8575,  /* MP */  0.9253,  /* ML */  2.3591,  /* MR */
    3.3514,  /* IL */  4.6203,  /* IR */ -0.1362,  /* S */
    5.0000,  /* E */
    0.9841,  /* B */
    5.0000   /* EL */
  };

  int n_capped = 0;
  int v;
  for (v = 0; v < M; v++) {
    int stt = cm->sttype[v];
    if (stt < 0 || stt >= 10) { pad_arr[v] = additive_pad; continue; }

    int cells = 1;
    if (cm->cp9map != NULL && cp9b != NULL) {
      int k = cm->cp9map->cs2hn[v][0];
      if (k >= 0 && k <= cp9b->hmm_M) {
        int iw = cp9b->imax[k] - cp9b->imin[k] + 1;
        int jw = cp9b->jmax[k] - cp9b->jmin[k] + 1;
        if (iw < 1) iw = 1;
        if (jw < 1) jw = 1;
        cells = iw * jw;
      }
    }

    double predicted = cal_a[stt] * (double)cells + cal_b[stt];
    int p = (int)ceil(predicted);
    if (p < 0) p = 0;
    p += additive_pad;
    if (maxpad > 0 && p > maxpad) { p = maxpad; n_capped++; }
    pad_arr[v] = p;
  }

  /* Distribution diagnostic so the sweep can tell "cap was active" from "no-op". */
  {
    int *sorted = malloc(sizeof(int) * M);
    if (sorted != NULL) {
      int i;
      for (i = 0; i < M; i++) sorted[i] = pad_arr[i];
      /* simple insertion sort -- M is small enough not to matter for diagnostic */
      for (i = 1; i < M; i++) {
        int key = sorted[i]; int j = i - 1;
        while (j >= 0 && sorted[j] > key) { sorted[j+1] = sorted[j]; j--; }
        sorted[j+1] = key;
      }
      int min_pad    = sorted[0];
      int max_pad    = sorted[M-1];
      int median_pad = sorted[M/2];
      int p90_pad    = sorted[(int)((M-1)*0.90)];
      int p99_pad    = sorted[(int)((M-1)*0.99)];
      fprintf(stderr, "#P7PB_PERSTATE M=%d min=%d median=%d p90=%d p99=%d max=%d capped=%d\n",
              M, min_pad, median_pad, p90_pad, p99_pad, max_pad, n_capped);
      free(sorted);
    }
  }

  return pad_arr;
}


/* Function: cm_BandsFromCYKParsetree()
 * Date    : 2026-04-08
 *
 * Purpose:  Derive per-CM-state bands (cp9b->imin/imax/jmin/jmax/hdmin/hdmax)
 *           directly from a CYK parsetree, without going through HMM-node
 *           pn bands. This avoids the linear HMM-node sweep contamination
 *           that breaks for permuted CMs (where the HMM-node order doesn't
 *           match the CM state order).
 *
 *           Algorithm:
 *           1. Walk parsetree, mark visited[v]=TRUE and set imin/imax/jmin/jmax
 *              from each occurrence of state v with absolute (i, j) coords.
 *           2. Apply pad to visited states, clamp to [i0..j0].
 *           3. For unvisited states, walk CM tree post-order (M-1 down to 0)
 *              and inherit from children: parent's bounds = child's bounds
 *              shifted by parent's StateLeftDelta/StateRightDelta.
 *              Bifurcations: i bound from left child, j bound from right child.
 *           4. Fill cp9b->imin/imax/jmin/jmax, then call cp9_GrowHDBands and
 *              ij2d_bands to compute hd bands.
 *
 *           For non-truncated mode, sets Jvalid[]=TRUE, others FALSE.
 *
 * Args:     cm        - the covariance model (must have cp9b allocated)
 *           errbuf    - error buffer
 *           tr        - parsetree from FastCYKScanHB_shmx; emitl/emitr in absolute dsq coords
 *           i0, j0    - envelope start/stop in absolute dsq coords (i0..j0)
 *           pad       - half-width pad to add on each side of visited bounds
 *           per_state_pad - if non-NULL, [0..M-1] per-state pad override; pad arg ignored for state v if per_state_pad[v] >= 0
 *           strict_unvisited - if TRUE, set Jvalid[v]=FALSE for unvisited states (--cykbands-strict).
 *                              This makes the F7 HB DP skip those states entirely, dramatically
 *                              reducing matrix cell count at the risk of parse failure or score loss.
 *           cp9b      - bands to fill (caller pre-allocated)
 *           pass_idx  - pipeline pass index (for truncation handling)
 *           debug     - if >0, print bands
 *
 * Returns:  eslOK on success.
 */
int
cm_BandsFromCYKParsetree(CM_t *cm, char *errbuf, Parsetree_t *tr,
                         int i0, int j0, int pad, const int *per_state_pad,
                         int strict_unvisited,
                         CP9Bands_t *cp9b, int pass_idx, int debug)
{
  int    status;
  int    M = cm->M;
  int    v, t, y, z, off;
  int    sdl, sdr;
  int    L = j0 - i0 + 1;
  int    do_trunc = cm_pli_PassAllowsTruncation(pass_idx);
  int   *visited = NULL;
  int   *skipped = NULL; /* Fix D (--cykskip-unvisited): parallel mark for unvisited states */
  int   *imin = cp9b->imin;
  int   *imax = cp9b->imax;
  int   *jmin = cp9b->jmin;
  int   *jmax = cp9b->jmax;

  ESL_ALLOC(visited, sizeof(int) * M);
  esl_vec_ISet(visited, M, FALSE);

  /* Initialize all bands to "unset" sentinel: imin/jmin = INT_MAX, imax/jmax = -1 */
  for(v = 0; v < M; v++) {
    imin[v] = INT_MAX; imax[v] = -1;
    jmin[v] = INT_MAX; jmax[v] = -1;
  }

  /* Step 1: Walk parsetree, set per-state bounds from each visit */
  for(t = 0; t < tr->n; t++) {
    v = tr->state[t];
    if(v < 0 || v >= M) continue;
    int i = tr->emitl[t];
    int j = tr->emitr[t];
    visited[v] = TRUE;
    if(i < imin[v]) imin[v] = i;
    if(i > imax[v]) imax[v] = i;
    if(j < jmin[v]) jmin[v] = j;
    if(j > jmax[v]) jmax[v] = j;
  }

  /* Step 2: Apply pad to visited states, clamp to [i0..j0].
   * If per_state_pad != NULL and per_state_pad[v] >= 0, use it instead of scalar pad. */
  for(v = 0; v < M; v++) {
    if(visited[v]) {
      int p = pad;
      if(per_state_pad != NULL && per_state_pad[v] >= 0) p = per_state_pad[v];
      imin[v] -= p; if(imin[v] < i0) imin[v] = i0;
      imax[v] += p; if(imax[v] > j0) imax[v] = j0;
      jmin[v] -= p; if(jmin[v] < i0) jmin[v] = i0;
      jmax[v] += p; if(jmax[v] > j0) jmax[v] = j0;
    }
  }

  /* Fix D (--cykskip-unvisited): mark unvisited states as SKIPPED.
   * Set bands to an empty range (jmin=1, jmax=0; imin=1, imax=0) so the
   * DP j-loop iterates over zero cells. Also mark visited[v]=TRUE so
   * Step 3 (the post-order inheritance walk below) leaves these states alone.
   * Final clamping is also gated to preserve the empty range.
   *
   * Aggressive optimization that commits to CYK MAP parsetree's subtree
   * choices: any truth alignment that visits states the CYK parse missed
   * cannot be recovered. Use with care; gated behind opt-in CLI flag. */
  if(cm->p7_cykskip_unvisited) {
    ESL_ALLOC(skipped, sizeof(int) * M);
    esl_vec_ISet(skipped, M, FALSE);
    for(v = 0; v < M; v++) {
      if(visited[v]) continue;
      imin[v] = 1; imax[v] = 0;
      jmin[v] = 1; jmax[v] = 0;
      skipped[v] = TRUE;
      visited[v] = TRUE; /* prevent Step 3 from inheriting */
    }
  }

  /* Step 3: For unvisited states, inherit from PARENT (top-down: 0 -> M-1).
   * This gives narrow bands near the MAP parse, unlike bottom-up which would
   * union over all children and degenerate to ~unbanded for sparse parsetrees.
   * Build parent[] map first. */
  int *parent = NULL;
  ESL_ALLOC(parent, sizeof(int) * M);
  for(v = 0; v < M; v++) parent[v] = -1;
  for(v = 0; v < M; v++) {
    if(cm->sttype[v] == B_st) {
      y = cm->cfirst[v]; if(y >= 0 && y < M) parent[y] = v;
      z = cm->cnum[v];   if(z >= 0 && z < M) parent[z] = v;
    } else if(cm->sttype[v] != E_st && cm->sttype[v] != EL_st) {
      for(off = 0; off < cm->cnum[v]; off++) {
        y = cm->cfirst[v] + off;
        if(y >= 0 && y < M && parent[y] == -1) parent[y] = v;
      }
    }
  }
  for(v = 0; v < M; v++) {
    if(visited[v]) continue;
    int p = parent[v];
    if(p < 0) {
      /* No parent (e.g., root) or unreachable; fall back to full envelope */
      imin[v] = i0; imax[v] = j0; jmin[v] = i0; jmax[v] = j0;
      continue;
    }
    /* Inherit from parent, shifted by parent's emit deltas:
     * If parent emits a left residue, child's i is parent's i + 1.
     * If parent emits a right residue, child's j is parent's j - 1.
     * For B parents: left child gets parent's i bounds, right child gets parent's j bounds. */
    if(cm->sttype[p] == B_st) {
      if(v == cm->cfirst[p]) {
        /* Left child of B: i bounds from parent's i, j is split point (we don't know exactly, use parent's range) */
        imin[v] = imin[p]; imax[v] = imax[p];
        jmin[v] = imin[p]; jmax[v] = jmax[p]; /* j can be anywhere parent's i to parent's j */
      } else {
        /* Right child of B: j bounds from parent's j, i is split point */
        imin[v] = imin[p]; imax[v] = jmax[p];
        jmin[v] = jmin[p]; jmax[v] = jmax[p];
      }
    } else {
      sdl = StateLeftDelta(cm->sttype[p]);
      sdr = StateRightDelta(cm->sttype[p]);
      imin[v] = imin[p] + sdl; if(imin[v] < i0) imin[v] = i0; if(imin[v] > j0) imin[v] = j0;
      imax[v] = imax[p] + sdl; if(imax[v] < i0) imax[v] = i0; if(imax[v] > j0) imax[v] = j0;
      jmin[v] = jmin[p] - sdr; if(jmin[v] < i0) jmin[v] = i0; if(jmin[v] > j0) jmin[v] = j0;
      jmax[v] = jmax[p] - sdr; if(jmax[v] < i0) jmax[v] = i0; if(jmax[v] > j0) jmax[v] = j0;
    }
  }
  free(parent);

  /* Final clamping and sanity: ensure imin<=imax, jmin<=jmax, all in [i0..j0].
   * Fix D: preserve empty bands for skipped[v] states. */
  for(v = 0; v < M; v++) {
    if(skipped != NULL && skipped[v]) continue;
    if(imin[v] == INT_MAX) { imin[v] = i0; imax[v] = j0; }
    if(jmin[v] == INT_MAX) { jmin[v] = i0; jmax[v] = j0; }
    if(imin[v] < i0) imin[v] = i0;
    if(imax[v] > j0) imax[v] = j0;
    if(jmin[v] < i0) jmin[v] = i0;
    if(jmax[v] > j0) jmax[v] = j0;
    if(imin[v] > imax[v]) imax[v] = imin[v];
    if(jmin[v] > jmax[v]) jmax[v] = jmin[v];
  }
  if(skipped != NULL) { free(skipped); skipped = NULL; }

  /* Set Jvalid/Lvalid/Rvalid/Tvalid for non-truncated mode.
   * In strict_unvisited mode, only visited states get Jvalid=TRUE; the F7 HB DP
   * gates on Jvalid[v] before computing alpha cells, so unvisited states are
   * skipped entirely. The EL state (M) is always kept valid (filler for end-locals).
   */
  if(!do_trunc) {
    if(strict_unvisited) {
      for(v = 0; v < M; v++) cp9b->Jvalid[v] = visited[v];
      cp9b->Jvalid[M] = TRUE; /* EL state */
    } else {
      esl_vec_ISet(cp9b->Jvalid, M + 1, TRUE);
    }
    esl_vec_ISet(cp9b->Lvalid, M + 1, FALSE);
    esl_vec_ISet(cp9b->Rvalid, M + 1, FALSE);
    esl_vec_ISet(cp9b->Tvalid, M + 1, FALSE);
  } else {
    /* Truncated bands not supported by this function path */
    esl_vec_ISet(cp9b->Jvalid, M + 1, TRUE);
    esl_vec_ISet(cp9b->Lvalid, M + 1, FALSE);
    esl_vec_ISet(cp9b->Rvalid, M + 1, FALSE);
    esl_vec_ISet(cp9b->Tvalid, M + 1, FALSE);
  }

  /* strict-alloc: shrink Jvalid=FALSE states' i,j bands to a single
   * sentinel cell so the matrix allocator doesn't reserve memory for
   * cells the DP will never compute. ij2d_bands below recomputes
   * hdmin/hdmax from the shrunk imin/imax, so the d-band collapses too.
   * The DP's existing Jvalid gating prevents the sentinel cells from
   * being read for actual scores.
   */
  if(!do_trunc && strict_unvisited) {
    for(v = 0; v < M; v++) {
      if(!cp9b->Jvalid[v]) {
        cp9b->imin[v] = cp9b->imax[v] = i0;
        cp9b->jmin[v] = cp9b->jmax[v] = i0;
      }
    }
  }

  cp9b->tau = cm->tau;

  /* Compute hdmin/hdmax */
  if((status = cp9_GrowHDBands(cp9b, errbuf)) != eslOK) goto ERROR;
  ij2d_bands(cm, cp9b, do_trunc, debug);

  free(visited);
  return eslOK;

 ERROR:
  if(visited) free(visited);
  if(skipped) free(skipped);
  return status;
}


/* Function: cp9_PredictStartAndEndPositionsP7B()
 * Date    : 2026-03-20
 *
 * Purpose:  Banded version of cp9_PredictStartAndEndPositions().
 *           Computes sp1/sp2 (first HMM nodes with significant occupancy)
 *           and ep1/ep2 (last such nodes) from banded CP9 posterior matrix,
 *           then derives Rmarg_imin/imax and Lmarg_jmin/jmax.
 *
 *           The only difference from the unbanded version is that the
 *           occupancy sum over rows for a given node k only visits rows
 *           where k is within the P7 band [kmin[i]..kmax[i]], and uses
 *           band-relative indexing to access the banded matrix.
 *
 * Args:     pmx   - banded CP9 posterior matrix (from cp9_PosteriorP7B)
 *           cp9b  - CP9 bands, with pn_min/pn_max already filled
 *           kmin  - [0..L] P7-derived min node for each row
 *           kmax  - [0..L] P7-derived max node for each row
 *           i0    - first position in original sequence coords
 *           j0    - last position in original sequence coords
 *
 * Returns:  void. cp9b->sp1/sp2/ep1/ep2 and Rmarg/Lmarg bounds are set.
 */
void
cp9_PredictStartAndEndPositionsP7B(CP9_MX *pmx, CP9Bands_t *cp9b, int *kmin, int *kmax, int i0, int j0)
{
  int i;
  int k;
  int L = j0-i0+1;
  int   iocc;
  float pocc;

  /* Part 1: Find sp1/sp2 — first nodes (left to right) with significant occupancy */
  k = 1;
  cp9b->sp1 = cp9b->sp2 = -1;
  while(k <= cp9b->hmm_M && (cp9b->sp1 == -1 || cp9b->sp2 == -1)) {
    if(cp9b->pn_min_m[k] == -1 && cp9b->pn_min_i[k] == -1 && cp9b->pn_min_d[k] == -1) {
      k++;
    }
    else {
      iocc = -INFTY;
      for(i = 0; i <= L; i++) {
	if(k >= kmin[i] && k <= kmax[i]) {
	  int kp = k - kmin[i];
	  iocc = ILogsum(iocc, ILogsum(pmx->mmx[i][kp], pmx->dmx[i][kp]));
	}
      }
      pocc = Score2Prob(iocc, 1.);
      if((cp9b->sp1 == -1) && (pocc > cp9b->thresh1)) cp9b->sp1 = k;
      if((cp9b->sp2 == -1) && (pocc > cp9b->thresh2)) cp9b->sp2 = k;
      k++;
    }
  }
  if(k == cp9b->hmm_M+1) {
    if(cp9b->sp1 == -1) { cp9b->sp1 = cp9b->hmm_M+1; }
    if(cp9b->sp2 == -1) { cp9b->sp2 = cp9b->hmm_M+1; }
  }

  /* Part 2: Find ep1/ep2 — last nodes (right to left) with significant occupancy */
  if((cp9b->sp1 == cp9b->hmm_M+1) &&
     (cp9b->sp2 == cp9b->hmm_M+1)) {
    cp9b->ep1 = 0;
    cp9b->ep2 = 0;
  }
  else {
    cp9b->ep1 = cp9b->ep2 = -1;
    k = cp9b->hmm_M;
    while(k >= 1 && (cp9b->ep1 == -1 || cp9b->ep2 == -1)) {
      if(cp9b->pn_min_m[k] == -1 && cp9b->pn_min_i[k] == -1 && cp9b->pn_min_d[k] == -1) {
	k--;
      }
      else {
	iocc = -INFTY;
	for(i = 0; i <= L; i++) {
	  if(k >= kmin[i] && k <= kmax[i]) {
	    int kp = k - kmin[i];
	    iocc = ILogsum(iocc, ILogsum(pmx->mmx[i][kp], pmx->dmx[i][kp]));
	  }
	}
	pocc = Score2Prob(iocc, 1.);
	if((cp9b->ep1 == -1) && (pocc > cp9b->thresh1)) cp9b->ep1 = k;
	if((cp9b->ep2 == -1) && (pocc > cp9b->thresh2)) cp9b->ep2 = k;
	k--;
      }
    }
    if(k == 0) {
      if(cp9b->ep1 == -1) { cp9b->ep1 = 0; }
      if(cp9b->ep2 == -1) { cp9b->ep2 = 0; }
    }
  }

  /* Parts 3-4: Derive Rmarg_imin/imax from sp1/sp2, Lmarg_jmin/jmax from ep1/ep2.
   * These only use pn_min/pn_max arrays (not the posterior matrix), so they are
   * identical to the unbanded version.
   */

  /* set cp9b->Rmarg_imin */
  if(cp9b->sp1 == cp9b->hmm_M+1) { cp9b->Rmarg_imin = i0; }
  else {
    cp9b->Rmarg_imin = INT_MAX;
    if(cp9b->sp1 != (cp9b->hmm_M+1) && cp9b->pn_min_m[cp9b->sp1] >= 0) cp9b->Rmarg_imin = ESL_MIN(cp9b->Rmarg_imin, cp9b->pn_min_m[cp9b->sp1]);
    if(cp9b->sp1 != (cp9b->hmm_M+1) && cp9b->pn_min_i[cp9b->sp1] >= 0) cp9b->Rmarg_imin = ESL_MIN(cp9b->Rmarg_imin, cp9b->pn_min_i[cp9b->sp1]);
    if(cp9b->sp1 != (cp9b->hmm_M+1) && cp9b->pn_min_d[cp9b->sp1] >= 0) cp9b->Rmarg_imin = ESL_MIN(cp9b->Rmarg_imin, cp9b->pn_min_d[cp9b->sp1]);
    if(cp9b->sp2 != (cp9b->hmm_M+1) && cp9b->pn_min_m[cp9b->sp2] >= 0) cp9b->Rmarg_imin = ESL_MIN(cp9b->Rmarg_imin, cp9b->pn_min_m[cp9b->sp2]);
    if(cp9b->sp2 != (cp9b->hmm_M+1) && cp9b->pn_min_i[cp9b->sp2] >= 0) cp9b->Rmarg_imin = ESL_MIN(cp9b->Rmarg_imin, cp9b->pn_min_i[cp9b->sp2]);
    if(cp9b->sp2 != (cp9b->hmm_M+1) && cp9b->pn_min_d[cp9b->sp2] >= 0) cp9b->Rmarg_imin = ESL_MIN(cp9b->Rmarg_imin, cp9b->pn_min_d[cp9b->sp2]);
    if(cp9b->Rmarg_imin == INT_MAX || cp9b->sp1 == (cp9b->hmm_M+1) || cp9b->sp2 == (cp9b->hmm_M+1)) cp9b->Rmarg_imin = i0;
    cp9b->Rmarg_imin = ESL_MAX(i0,   cp9b->Rmarg_imin);
    cp9b->Rmarg_imin = ESL_MIN(j0+1, cp9b->Rmarg_imin);
  }

  /* set cp9b->Rmarg_imax */
  if(cp9b->sp1 == cp9b->hmm_M+1) { cp9b->Rmarg_imax = j0; }
  else {
    cp9b->Rmarg_imax = INT_MIN;
    if(cp9b->sp1 != (cp9b->hmm_M+1) && cp9b->pn_max_m[cp9b->sp1] >= 0) cp9b->Rmarg_imax = ESL_MAX(cp9b->Rmarg_imax, cp9b->pn_max_m[cp9b->sp1]);
    if(cp9b->sp1 != (cp9b->hmm_M+1) && cp9b->pn_max_i[cp9b->sp1] >= 0) cp9b->Rmarg_imax = ESL_MAX(cp9b->Rmarg_imax, cp9b->pn_max_i[cp9b->sp1]);
    if(cp9b->sp1 != (cp9b->hmm_M+1) && cp9b->pn_max_d[cp9b->sp1] >= 0) cp9b->Rmarg_imax = ESL_MAX(cp9b->Rmarg_imax, cp9b->pn_max_d[cp9b->sp1]);
    if(cp9b->sp2 != (cp9b->hmm_M+1) && cp9b->pn_max_m[cp9b->sp2] >= 0) cp9b->Rmarg_imax = ESL_MAX(cp9b->Rmarg_imax, cp9b->pn_max_m[cp9b->sp2]);
    if(cp9b->sp2 != (cp9b->hmm_M+1) && cp9b->pn_max_i[cp9b->sp2] >= 0) cp9b->Rmarg_imax = ESL_MAX(cp9b->Rmarg_imax, cp9b->pn_max_i[cp9b->sp2]);
    if(cp9b->sp2 != (cp9b->hmm_M+1) && cp9b->pn_max_d[cp9b->sp2] >= 0) cp9b->Rmarg_imax = ESL_MAX(cp9b->Rmarg_imax, cp9b->pn_max_d[cp9b->sp2]);
    if(cp9b->Rmarg_imax == INT_MIN || cp9b->sp1 == (cp9b->hmm_M+1) || cp9b->sp2 == (cp9b->hmm_M+1)) cp9b->Rmarg_imax = j0+1;
    cp9b->Rmarg_imax = ESL_MAX(i0,   cp9b->Rmarg_imax);
    cp9b->Rmarg_imax = ESL_MIN(j0+1, cp9b->Rmarg_imax);
  }

  /* set cp9b->Lmarg_jmin */
  if(cp9b->ep1 == 0) { cp9b->Lmarg_jmin = i0-1; }
  else {
    cp9b->Lmarg_jmin = INT_MAX;
    if(cp9b->ep1 != 0 && cp9b->pn_min_m[cp9b->ep1] >= 0) cp9b->Lmarg_jmin = ESL_MIN(cp9b->Lmarg_jmin, cp9b->pn_min_m[cp9b->ep1]);
    if(cp9b->ep1 != 0 && cp9b->pn_min_i[cp9b->ep1] >= 0) cp9b->Lmarg_jmin = ESL_MIN(cp9b->Lmarg_jmin, cp9b->pn_min_i[cp9b->ep1]);
    if(cp9b->ep1 != 0 && cp9b->pn_min_d[cp9b->ep1] >= 0) cp9b->Lmarg_jmin = ESL_MIN(cp9b->Lmarg_jmin, cp9b->pn_min_d[cp9b->ep1]-1);
    if(cp9b->ep2 != 0 && cp9b->pn_min_m[cp9b->ep2] >= 0) cp9b->Lmarg_jmin = ESL_MIN(cp9b->Lmarg_jmin, cp9b->pn_min_m[cp9b->ep2]);
    if(cp9b->ep2 != 0 && cp9b->pn_min_i[cp9b->ep2] >= 0) cp9b->Lmarg_jmin = ESL_MIN(cp9b->Lmarg_jmin, cp9b->pn_min_i[cp9b->ep2]);
    if(cp9b->ep2 != 0 && cp9b->pn_min_d[cp9b->ep2] >= 0) cp9b->Lmarg_jmin = ESL_MIN(cp9b->Lmarg_jmin, cp9b->pn_min_d[cp9b->ep2]-1);
    if(cp9b->Lmarg_jmin == INT_MAX || cp9b->ep1 == 0 || cp9b->ep2 == 0) cp9b->Lmarg_jmin = i0-1;
    cp9b->Lmarg_jmin = ESL_MAX(i0-1, cp9b->Lmarg_jmin);
    cp9b->Lmarg_jmin = ESL_MIN(j0,   cp9b->Lmarg_jmin);
  }

  /* set cp9b->Lmarg_jmax */
  if(cp9b->ep1 == 0) { cp9b->Lmarg_jmax = j0; }
  else {
    cp9b->Lmarg_jmax = INT_MIN;
    if(cp9b->ep1 != 0 && cp9b->pn_max_m[cp9b->ep1] >= 0) cp9b->Lmarg_jmax = ESL_MAX(cp9b->Lmarg_jmax, cp9b->pn_max_m[cp9b->ep1]);
    if(cp9b->ep1 != 0 && cp9b->pn_max_i[cp9b->ep1] >= 0) cp9b->Lmarg_jmax = ESL_MAX(cp9b->Lmarg_jmax, cp9b->pn_max_i[cp9b->ep1]);
    if(cp9b->ep1 != 0 && cp9b->pn_max_d[cp9b->ep1] >= 0) cp9b->Lmarg_jmax = ESL_MAX(cp9b->Lmarg_jmax, cp9b->pn_max_d[cp9b->ep1]-1);
    if(cp9b->ep2 != 0 && cp9b->pn_max_m[cp9b->ep2] >= 0) cp9b->Lmarg_jmax = ESL_MAX(cp9b->Lmarg_jmax, cp9b->pn_max_m[cp9b->ep2]);
    if(cp9b->ep2 != 0 && cp9b->pn_max_i[cp9b->ep2] >= 0) cp9b->Lmarg_jmax = ESL_MAX(cp9b->Lmarg_jmax, cp9b->pn_max_i[cp9b->ep2]);
    if(cp9b->ep2 != 0 && cp9b->pn_max_d[cp9b->ep2] >= 0) cp9b->Lmarg_jmax = ESL_MAX(cp9b->Lmarg_jmax, cp9b->pn_max_d[cp9b->ep2]-1);
    if(cp9b->Lmarg_jmax == INT_MIN || cp9b->ep1 == 0 || cp9b->ep2 == 0) cp9b->Lmarg_jmax = j0;
    cp9b->Lmarg_jmax = ESL_MAX(i0-1, cp9b->Lmarg_jmax);
    cp9b->Lmarg_jmax = ESL_MIN(j0,   cp9b->Lmarg_jmax);
  }

  return;
}


/* Function: cp9_Seq2PosteriorsP7B
 * Date    : EPN, Tue Aug 19 13:12:12 2008
 *
 * Purpose:  Given a CM with precalc'ed CP9 HMM and CP9Map, and HMMER3 plan 7
 *           HMM bands for the CP9 HMM DP matrices, and a sequence,
 *           run HMM Forward and Backward algorithms, and return a CP9 posterior
 *           matrix.
 *           
 *           Note: this function was never updated to handle 
 *           truncated alignment (b/c it's currently not hooked up
 *           to any of the Infernal applications).
 *
 * Args:     cm           - the covariance model
 *           errbuf       - char buffer for error messages
 *           fmx          - CP9 dp matrix for Forward()
 *           bmx          - CP9 dp matrix for Backward()
 *           pmx          - CP9 dp matrix to fill with posteriors, can == bmx
 *           dsq          - sequence in digitized form
 *           L            - length of the sequence we're aligning (1..L)
 *           kmin         - P7 dervied band to enforce: [0.i..L] = k, min node k for residue i
 *           kmax         - P7 derived band to enforce: [0.i..L] = k, min node k for residue i
 *           debug_level  - verbosity level for debugging printf()s
 *           
 * Return:  eslOK on success
 */
int
cp9_Seq2PosteriorsP7B(CM_t *cm, char *errbuf, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *pmx, ESL_DSQ *dsq, int L, int *kmin, int *kmax, int debug_level)
{
  int status;
  float sc;
  CP9_t *cp9 = NULL;  /* ptr to cp9 HMM (this could be Lcp9, Rcp9, Tcp9 if we update this function to possibly handle truncated alignment) */

  /* Contract checks */
  if(dsq == NULL)        ESL_FAIL(eslEINCOMPAT, errbuf, "in cp9_Seq2Posteriors(), dsq is NULL.");
  if(cm->cp9map == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "in cp9_Seq2Posteriors, but cm->cp9map is NULL.\n");
  if((cm->align_opts & CM_ALIGN_HBANDED) && (cm->search_opts & CM_SEARCH_HBANDED)) 
    ESL_FAIL(eslEINCOMPAT, errbuf, "in cp9_Seq2Posteriors, CM_ALIGN_HBANDED and CM_SEARCH_HBANDED flags both up, exactly 1 must be up.\n");
  if((cm->search_opts & CM_SEARCH_HMMALNBANDS) && (! (cm->search_opts & CM_SEARCH_HBANDED))) 
    ESL_FAIL(eslEINCOMPAT, errbuf, "in cp9_Seq2Posteriors, CM_SEARCH_HMMALNBANDS flag raised, but not CM_SEARCH_HBANDED flag, this doesn't make sense\n");

  /* determine which cp9 HMM to use, if CM is has local begins use cp9 (its local too) else use cp9glb (its global) */
  cp9 = cm->cp9;
  if(cp9 == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cp9_Seq2Posteriors, relevant cp9 is NULL.\n");

  /* Step 1: Get HMM posteriors.*/
  if((status = cp9_ForwardP7B (cp9, errbuf, fmx, dsq, L, kmin, kmax, &sc)) != eslOK) return status;
  if(debug_level > 0) printf("CP9P7B Forward  score : %.4f\n", sc);

  if((status = cp9_BackwardP7B(cp9, errbuf, bmx, dsq, L, kmin, kmax, &sc)) != eslOK) return status;
  if(debug_level > 0) printf("CP9 Backward  score : %.4f\n", sc);

  if(cm->align_opts & CM_ALIGN_CHECKFB) {
    if((status = cp9_CheckFBP7B(fmx, bmx, cp9, errbuf, sc, 1, L, dsq, kmin, kmax)) != eslOK) return status;
    printf("Forward/Backward matrices checked.\n");
  }

  /* Get posteriors */
  if((status = cp9_PosteriorP7B(dsq, errbuf, L, cp9, fmx, bmx, pmx, kmin, kmax)) != eslOK) return status;

  return eslOK;
}


/* Function: cp9_PosteriorP7B()
 * based on Ian Holmes' hmmer/src/postprob.c::P7EmitterPosterior()
 *
 * Purpose:  Combines HMMER3 p7 banded Forward and Backward matrices into a 
 *           posterior probability matrix. For emitters (match and inserts) the 
 *           entries in row i of this matrix are the logs of the posterior 
 *           probabilities of each state emitting symbol i of the sequence. 
 *           For non-emitters the entries in row i of this matrix are the 
 *           logs of the posterior probabilities of each state being 'visited' 
 *           when the last emitted residue in the parse was symbol i of the
 *           sequence.
 *
 * Args:     dsq      - sequence in digitized form
 *           errbuf   - for error messages
 *           L        - length of target subsequence 1..L
 *           hmm      - the model
 *           forward  - pre-calculated forward matrix
 *           backward - pre-calculated backward matrix
 *           pm       - pre-allocated dynamic programming matrix for posteriors
 *           kmin     - P7 derived band to enforce: [0.i..L] = k, min node k for residue i
 *           kmax     - P7 derived band to enforce: [0.i..L] = k, min node k for residue i
 *           
 * Return:   eslOK on success;
 */
int
cp9_PosteriorP7B(ESL_DSQ *dsq, char *errbuf, int L, CP9_t *hmm, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *pmx, int *kmin, int *kmax)
{
  int i;
  int k;
  int sc;
  int M = hmm->M;
  int kp, kn, kx; 
  /*float temp_sc;*/

  if(bmx != pmx) GrowCP9Matrix(pmx, errbuf, L, M, kmin, kmax, NULL, NULL, NULL, NULL, NULL);

  /* parses must start/stop at (i = 1)/(j = L) */
  sc = bmx->mmx[0][0]; 

  /* note boundary conditions, i = 1 */
  assert(kmin[0] == 0);
  pmx->mmx[0][0] = fmx->mmx[0][0] + bmx->mmx[0][0] - sc; /* fmx->mmx[0][0] is 0, bmx->mmx[1][0] is overall score */
  pmx->imx[0][0] = -INFTY; /*need seq to get here*/
  pmx->dmx[0][0] = -INFTY; /*D_0 does not exist*/
  i  = 0;
  kn = ESL_MAX(kmin[i], 1);
  kx = kmax[i];
  kp = kn - kmin[i];
  for (k = kn; k <= kx; k++, kp++) {
    pmx->mmx[0][kp] = -INFTY; /*need seq to get here*/
    pmx->imx[0][kp] = -INFTY; /*need seq to get here*/
    pmx->dmx[0][kp] = fmx->dmx[0][kp] + bmx->dmx[0][kp] - sc;
  }

  for (i = 1; i <= L; i++) {
      k = 0;
      if(INBAND(i,0)) { 
	kp = k - kmin[i];
	assert(kp == 0);
	pmx->mmx[i][kp] = ESL_MAX(fmx->mmx[i][kp] + bmx->mmx[i][kp] - sc, -INFTY); /* M_0 doesn't emit */
	pmx->imx[i][kp] = ESL_MAX(fmx->imx[i][kp] + bmx->imx[i][kp] - hmm->isc[dsq[i]][0] - sc, -INFTY);
	/*hmm->isc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
	pmx->dmx[i][kp] = -INFTY; /* D_0 doesn't exist */
      }

      kn = ESL_MAX(kmin[i], 1);
      kx = kmax[i];
      kp = kn - kmin[i];
      for(k = kn; k <= kx; k++, kp++)
	{
	  pmx->mmx[i][kp] = ESL_MAX(fmx->mmx[i][kp] + bmx->mmx[i][kp] - hmm->msc[dsq[i]][k] - sc, -INFTY);
	  /*hmm->msc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
	  pmx->imx[i][kp] = ESL_MAX(fmx->imx[i][kp] + bmx->imx[i][kp] - hmm->isc[dsq[i]][k] - sc, -INFTY);
	  /*hmm->isc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
	  pmx->dmx[i][kp] = ESL_MAX(fmx->dmx[i][kp] + bmx->dmx[i][kp] - sc, -INFTY);
	}
    }	  

  /*
    float temp_sc;
    for(i = 0; i <= L; i++)
    {
    kp = 0;
    for(k = kmin[i]; k <= kmax[i]; k++, kp++)
    {
    temp_sc = Score2Prob(mx->mmx[i][kp], 1.);
    if(temp_sc > .0001)
    printf("mx->mmx[%3d][%3d]: %9d | %8f\n", i, k, mx->mmx[i][kp], temp_sc);
    temp_sc = Score2Prob(mx->imx[i][kp], 1.);
    if(temp_sc > .0001)
    printf("mx->imx[%3d][%3d]: %9d | %8f\n", i, k, mx->imx[i][kp], temp_sc);
    temp_sc = Score2Prob(mx->dmx[i][kp], 1.);
    if(temp_sc > .0001)
    printf("mx->dmx[%3d][%3d]: %9d | %8f\n", i, k, mx->dmx[i][kp], temp_sc);
    }
    }*/
  return eslOK;
}

/* Function: cp9_FB2HMMBandsP7B()
 * Date:     EPN, Fri Aug 15 14:00:59 2008
 *
 * Purpose: Determine the band on all HMM states given HMMER3 Plan 7 Banded 
 *          Forward and Backward matrices. Do this by calculating and summing 
 *          log posterior probabilities that each state emitted/was visited at each posn,
 *          starting at the band ends, and creeping in, until the half the
 *          maximum allowable probability excluded is reached on each side.
 *
 * Args:
 *
 * CP9_t hmm        the HMM
 * errbuf           char buffer for error messages
 * CP9_MX fmx:      forward DP matrix, already calc'ed
 * CP9_MX bmx:      backward DP matrix, already calc'ed
 * CP9_MX pmx:      DP matrix for posteriors, filled here, can == bmx
 * dsq              the digitized sequence
 * CP9Bands_t cp9b  CP9 bands data structure
 * int   L          length of target subsequence (1..L)
 * int   M          number of nodes in HMM (num columns of pmx matrix)
 * double p_thresh  the probability mass we're requiring is within each band
 * int do_old_hmm2ij TRUE if we'll use old cp9_HMM2ijBands_OLD() function downstream
 * int kmin         P7 dervied band to enforce: [0.i..L] = k, min node k for residue i
 * int kmax         P7 derived band to enforce: [0.i..L] = k, min node k for residue i
 * int debug_level  [0..3] tells the function what level of debugging print
 *                  statements to print.
 * 
 * Returns: eslOK on success;
 */
int
cp9_FB2HMMBandsP7B(CP9_t *hmm, char *errbuf, ESL_DSQ *dsq, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *pmx, CP9Bands_t *cp9b,
		   int L, int M, double p_thresh, int do_old_hmm2ij, int *kmin, int *kmax, int debug_level,
		   int do_pnmono, int do_pnmono_print)
{
  int status;
  int k;                                  /* counter over nodes of the model */
  int thresh = Prob2Score(((1. - p_thresh)/2.), 1.); /* allowable prob mass excluded on each side */

  /* *_m = match, *_i = insert, *_d = delete */
  int *kthresh_m, *kthresh_i, *kthresh_d; /* [0..k..hmm->M], individual thresholds for each state */
  int *nset_m, *nset_i, *nset_d;          /* [0..k..hmm->M], has minimum been set for this state? */
  int *xset_m, *xset_i, *xset_d;          /* [0..k..hmm->M], has maximum been set for this state? */
  int *mass_m, *mass_i, *mass_d;          /* [0..k..hmm->M], summed log prob of pmx->mx[i][k] from 0..k or k..L */
  int i;                                  /* actual position */
  int sc;                                 /* summed score of all parses (derived from backward matrix) 
					   * if(cm->search_opts & CM_SEARCH_HMMALNBANDS) Forward and Backward
					   * were run in 'scan mode' where each residue can be begin/end of a parse,
					   * so we have to sum up parses that end at each posn, 
					   * if ! (cm->search_opts & CM_SEARCH_HMMALNBANDS) we know we have 
					   * to start at residue 1 and end at residue L, so sc is simply bmx->mmx[0][0]
					   */
  int hmm_is_localized;                   /* TRUE if HMM has local begins, ends or ELs on */
  int kp, kn, kx;

  hmm_is_localized = ((hmm->flags & CPLAN9_LOCAL_BEGIN) || (hmm->flags & CPLAN9_LOCAL_END) || (hmm->flags & CPLAN9_EL)) ? TRUE : FALSE;

  if(bmx != pmx) GrowCP9Matrix(pmx, errbuf, L, M, kmin, kmax, NULL, NULL, NULL, NULL, NULL);

  /* allocations and initializations */
  ESL_ALLOC(nset_m, sizeof(int) * (M+1));
  ESL_ALLOC(nset_i, sizeof(int) * (M+1));
  ESL_ALLOC(nset_d, sizeof(int) * (M+1));
  ESL_ALLOC(xset_m, sizeof(int) * (M+1));
  ESL_ALLOC(xset_i, sizeof(int) * (M+1));
  ESL_ALLOC(xset_d, sizeof(int) * (M+1));
  ESL_ALLOC(mass_m, sizeof(int) * (M+1));
  ESL_ALLOC(mass_i, sizeof(int) * (M+1));
  ESL_ALLOC(mass_d, sizeof(int) * (M+1));  
  ESL_ALLOC(kthresh_m, sizeof(int) * (M+1));
  ESL_ALLOC(kthresh_i, sizeof(int) * (M+1));
  ESL_ALLOC(kthresh_d, sizeof(int) * (M+1));  

  esl_vec_ISet(mass_m, M+1, -INFTY);
  esl_vec_ISet(mass_i, M+1, -INFTY);
  esl_vec_ISet(mass_d, M+1, -INFTY);
  esl_vec_ISet(nset_m, M+1, FALSE);
  esl_vec_ISet(nset_i, M+1, FALSE);
  esl_vec_ISet(nset_d, M+1, FALSE);
  esl_vec_ISet(xset_m, M+1, FALSE);
  esl_vec_ISet(xset_i, M+1, FALSE);
  esl_vec_ISet(xset_d, M+1, FALSE);

  sc = bmx->mmx[0][0]; /* Forward/Backward run in 'align mode' parses must start at 1, end at L */
  /* sc is summed log prob of all possible parses of seq 1..L */

  /* note boundary conditions, i = 1 */
  assert(kmin[0] == 0);
  pmx->mmx[0][0] = fmx->mmx[0][0] + bmx->mmx[0][0] - sc; /* fmx->mmx[0][0] is 0, bmx->mmx[1][0] is overall score */
  pmx->imx[0][0] = -INFTY; /*need seq to get here*/
  pmx->dmx[0][0] = -INFTY; /*D_0 does not exist*/
  if((mass_m[0] = pmx->mmx[0][0]) > thresh) { 
    cp9b->pn_min_m[0] = 0; 
    nset_m[0] = TRUE; 
  }
  mass_i[0] = -INFTY; /* b/c pmx->imx[0][0] is -INFTY, set above */
  mass_d[0] = -INFTY; /* b/c pmx->dmx[0][0] is -INFTY, set above */

  i  = 0;
  kn = ESL_MAX(kmin[i], 1);
  kx = kmax[i];
  kp = kn - kmin[i];
  for (k = kn; k <= kx; k++, kp++) {
    pmx->mmx[0][kp] = -INFTY; /*need seq to get here*/
    pmx->imx[0][kp] = -INFTY; /*need seq to get here*/
    pmx->dmx[0][kp] = fmx->dmx[0][kp] + bmx->dmx[0][kp] - sc;
    /* mass_m[k] doesn't change b/c pmx->mmx[0][kp] is -INFTY */
    /* mass_i[k] doesn't change b/c pmx->imx[0][kp] is -INFTY */
    if((mass_d[k] = pmx->dmx[0][kp]) > thresh) { 
      cp9b->pn_min_d[k] = 0;
      nset_d[k] = TRUE; 
    }
  }

  /* Find minimum position in band for each state (M,I,D) of each node (0..M) */
  for (i = 1; i <= L; i++) {
      k = 0;
      if(INBAND(i,0)) { 
	kp = k - kmin[i];
	assert(kp == 0);
	pmx->mmx[i][kp] = ESL_MAX(fmx->mmx[i][kp] + bmx->mmx[i][kp] - sc, -INFTY); /* M_0 doesn't emit */
	if(! nset_m[0]) { 
	  if((mass_m[0] = ILogsum(mass_m[0], pmx->mmx[i][kp])) > thresh) { 
	    cp9b->pn_min_m[0] = i;
	    nset_m[0] = TRUE; 
	  }
	}
	
	pmx->imx[i][kp] = ESL_MAX(fmx->imx[i][kp] + bmx->imx[i][kp] - hmm->isc[dsq[i]][0] - sc, -INFTY);
	/*hmm->isc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
	if(! nset_i[0]) { 
	  if((mass_i[0] = ILogsum(mass_i[0], pmx->imx[i][kp])) > thresh) { 
	    cp9b->pn_min_i[0] = i;
	    nset_i[0] = TRUE; 
	  }
	}
	pmx->dmx[i][kp] = -INFTY; /* D_0 doesn't exist */
      }

      kn = ESL_MAX(kmin[i], 1);
      kx = kmax[i];
      kp = kn - kmin[i];
      for(k = kn; k <= kx; k++, kp++)
	{
	  pmx->mmx[i][kp] = ESL_MAX(fmx->mmx[i][kp] + bmx->mmx[i][kp] - hmm->msc[dsq[i]][k] - sc, -INFTY);
	  /*hmm->msc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
	  pmx->imx[i][kp] = ESL_MAX(fmx->imx[i][kp] + bmx->imx[i][kp] - hmm->isc[dsq[i]][k] - sc, -INFTY);
	  /*hmm->isc[dsq[i]][k] will have been counted in both fmx->mmx and bmx->mmx*/
	  pmx->dmx[i][kp] = ESL_MAX(fmx->dmx[i][kp] + bmx->dmx[i][kp] - sc, -INFTY);

	  if(! nset_m[k]) { 
	    if((mass_m[k] = ILogsum(mass_m[k], pmx->mmx[i][kp])) > thresh) { 
	      cp9b->pn_min_m[k] = i;
	      nset_m[k] = TRUE; 
	    }
	  }
	  if(! nset_i[k]) { 
	    if((mass_i[k] = ILogsum(mass_i[k], pmx->imx[i][kp])) > thresh) { 
	      cp9b->pn_min_i[k] = i;
	      nset_i[k] = TRUE; 
	    }
	  }
	  if(! nset_d[k]) { 
	    if((mass_d[k] = ILogsum(mass_d[k], pmx->dmx[i][kp])) > thresh) { 
	      cp9b->pn_min_d[k] = i;
	      nset_d[k] = TRUE; 
	    }
	  }
	}
    }	  

  esl_vec_ISet(mass_m, M+1, -INFTY);
  esl_vec_ISet(mass_i, M+1, -INFTY);
  esl_vec_ISet(mass_d, M+1, -INFTY);
  /* Find maximum position in band for each state (M,I,D) of each node (0..M)
   * by moving from L down to 1 */
  for (i = L; i >= 1; i--) /* i is the relative position in the seq */
    {
      kp = 0;
      for(k = kmin[i]; k <= kmax[i]; k++, kp++)
	{
	  if(! xset_m[k]) { 
	    if((mass_m[k] = ILogsum(mass_m[k], pmx->mmx[i][kp])) > thresh) { 
	      cp9b->pn_max_m[k] = i;
	      xset_m[k] = TRUE; 
	    }
	  }
	  if(! xset_i[k]) { 
	    if((mass_i[k] = ILogsum(mass_i[k], pmx->imx[i][kp])) > thresh) { 
	      cp9b->pn_max_i[k] = i;
	      xset_i[k] = TRUE; 
	    }
	  }
	  if(! xset_d[k]) { 
	    if((mass_d[k] = ILogsum(mass_d[k], pmx->dmx[i][kp])) > thresh) { 
	      cp9b->pn_max_d[k] = i;
	      xset_d[k] = TRUE; 
	    }
	  }
	}
    }	  
  /* note boundary conditions, i = 0 */
  if(INBAND(0,0)) { 
    assert(kmin[0] == 0);
    if(! xset_m[0]) { 
      if((mass_m[0] = ILogsum(mass_m[0], pmx->mmx[0][0])) > thresh) { 
	cp9b->pn_max_m[0] = 0; 
	xset_m[0] = TRUE; 
      }
    }
  }
  /* mass_i[0] is unchaged because b/c pmx->imx[0][0] is -INFTY, set above */
  /* mass_d[0] is unchaged because b/c pmx->dmx[0][0] is -INFTY, set above */
  kn = ESL_MAX(kmin[i], 1);
  kx = kmax[i];
  kp = kn - kmin[i];
  for(k = kn; k <= kx; k++, kp++) {
    /* mass_m[k] doesn't change b/c pmx->mmx[0][k] is -INFTY */
    /* mass_i[k] doesn't change b/c pmx->mmx[0][k] is -INFTY */
    if(!xset_d[k]) { 
      if((mass_d[k] = ILogsum(mass_d[k], pmx->dmx[0][kp])) > thresh) { 
	cp9b->pn_max_d[k] = 0;
	xset_d[k] = TRUE; 
      }
    }
  }	 

  /* new technique as of EPN, Sun Jan 27 08:48:34 2008 */
  /* Some states may not have had their min/max set. This occurs if the entire
   * state is outside the band (i.e. the summed probability the state is entered for ANY i
   * is less than our threshold. Current strategy in this situation is to set the
   * pn_min_* and pn_max_* values as special flags, (-2) so the function that
   * uses them to derive i and j bands knows this is the case and handles it
   * accordingly.
   */
  int mset;
  int dset;
  for(k = 0; k <= M; k++) { 
    mset = dset = TRUE;
    /* theoretically either nset_*[k] and xset_*[k] should be either both TRUE or both
     * FALSE, but I'm slightly worried about rare precision issues, so we check if one 
     * or the other is unset, and if so, we set both to argmax position */
    if(((! nset_m[k])) || (! xset_m[k]) || (cp9b->pn_max_m[k] < cp9b->pn_min_m[k])) { 
      cp9b->pn_min_m[k] = cp9b->pn_max_m[k] = -1;
      mset = FALSE;
    }
    if(((! nset_i[k])) || (! xset_i[k]) || (cp9b->pn_max_i[k] < cp9b->pn_min_i[k])) { 
      cp9b->pn_min_i[k] = cp9b->pn_max_i[k] = -1;
    }
    if(((! nset_d[k])) || (! xset_d[k]) || (cp9b->pn_max_d[k] < cp9b->pn_min_d[k])) { 
      cp9b->pn_min_d[k] = cp9b->pn_max_d[k] = -1;
      dset = FALSE;
    }
    if((!hmm_is_localized) && (mset == FALSE && dset == FALSE)) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "node: %d match nor delete HMM state bands were set in non-localized, non-scanning HMM, lower tau (should be << 0.5).\n", k);
  }
  
  cp9b->pn_min_d[0] = -1; /* D_0 doesn't exist */
  cp9b->pn_max_d[0] = -1; /* D_0 doesn't exist */

  if(do_pnmono) pn_match_bands_enforce_monotone(cp9b->pn_min_m, cp9b->pn_max_m, M, L, do_pnmono_print, "fb2hmm_p7b");

  /* Always print HMM bands for P7B path debugging */
  if(debug_level > 0) cp9_DebugPrintHMMBands(stdout, L, cp9b, (1.-p_thresh), 1);

  free(mass_m);
  free(mass_i);
  free(mass_d);
  free(nset_m);
  free(nset_i);
  free(nset_d);
  free(xset_m);
  free(xset_i);
  free(xset_d);
  free(kthresh_m);
  free(kthresh_i);
  free(kthresh_d);

  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
}


/* Function: cp9_FB2HMMBandsP7BF()
 *
 * Float-precision mirror of cp9_FB2HMMBandsP7B(). Reads CP9_FMX fmx/bmx,
 * computes per-cell posteriors into pmx (CP9_FMX), and derives per-state
 * pn_min/pn_max bands by sweeping for the position where accumulated
 * posterior mass crosses (1-p_thresh)/2 from each side.
 *
 * Logic identical to the int version. Differences are mechanical:
 * ILogsum -> p7_FLogsum, -INFTY -> -eslINFINITY, model scores via Scorify().
 */
int
cp9_FB2HMMBandsP7BF(CP9_t *hmm, char *errbuf, ESL_DSQ *dsq, CP9_FMX *fmx, CP9_FMX *bmx, CP9_FMX *pmx, CP9Bands_t *cp9b,
		    int L, int M, double p_thresh, int do_old_hmm2ij, int *kmin, int *kmax, int debug_level,
		    int do_pnmono, int do_pnmono_print)
{
  int status;
  int k;
  /* float threshold in log space (nats? bits? cp9 uses scoring as bits/INTSCALE; we keep log-prob) */
  float thresh = logf((1. - p_thresh) / 2.);

  float *nset_m, *nset_i, *nset_d;
  float *xset_m, *xset_i, *xset_d;
  float *mass_m, *mass_i, *mass_d;
  int   *_nset_m, *_nset_i, *_nset_d;
  int   *_xset_m, *_xset_i, *_xset_d;
  int    i;
  float  sc;
  int    hmm_is_localized;
  int    kp, kn, kx;

  hmm_is_localized = ((hmm->flags & CPLAN9_LOCAL_BEGIN) || (hmm->flags & CPLAN9_LOCAL_END) || (hmm->flags & CPLAN9_EL)) ? TRUE : FALSE;

  if(bmx != pmx) GrowCP9FMatrix(pmx, errbuf, L, M, kmin, kmax, NULL, NULL, NULL, NULL, NULL);

  /* mass_X are float log-prob accumulators; nset_X / xset_X are int boolean flags */
  ESL_ALLOC(_nset_m, sizeof(int) * (M+1));
  ESL_ALLOC(_nset_i, sizeof(int) * (M+1));
  ESL_ALLOC(_nset_d, sizeof(int) * (M+1));
  ESL_ALLOC(_xset_m, sizeof(int) * (M+1));
  ESL_ALLOC(_xset_i, sizeof(int) * (M+1));
  ESL_ALLOC(_xset_d, sizeof(int) * (M+1));
  ESL_ALLOC(mass_m, sizeof(float) * (M+1));
  ESL_ALLOC(mass_i, sizeof(float) * (M+1));
  ESL_ALLOC(mass_d, sizeof(float) * (M+1));

  /* keep the unused float aliases NULL to silence unused-var warnings; not all callers want them */
  nset_m = nset_i = nset_d = NULL;
  xset_m = xset_i = xset_d = NULL;
  (void)nset_m; (void)nset_i; (void)nset_d;
  (void)xset_m; (void)xset_i; (void)xset_d;

  esl_vec_FSet(mass_m, M+1, -eslINFINITY);
  esl_vec_FSet(mass_i, M+1, -eslINFINITY);
  esl_vec_FSet(mass_d, M+1, -eslINFINITY);
  esl_vec_ISet(_nset_m, M+1, FALSE);
  esl_vec_ISet(_nset_i, M+1, FALSE);
  esl_vec_ISet(_nset_d, M+1, FALSE);
  esl_vec_ISet(_xset_m, M+1, FALSE);
  esl_vec_ISet(_xset_i, M+1, FALSE);
  esl_vec_ISet(_xset_d, M+1, FALSE);

  sc = bmx->mmx[0][0];

  assert(kmin[0] == 0);
  pmx->mmx[0][0] = fmx->mmx[0][0] + bmx->mmx[0][0] - sc;
  pmx->imx[0][0] = -eslINFINITY;
  pmx->dmx[0][0] = -eslINFINITY;
  if((mass_m[0] = pmx->mmx[0][0]) > thresh) {
    cp9b->pn_min_m[0] = 0;
    _nset_m[0] = TRUE;
  }
  mass_i[0] = -eslINFINITY;
  mass_d[0] = -eslINFINITY;

  i  = 0;
  kn = ESL_MAX(kmin[i], 1);
  kx = kmax[i];
  kp = kn - kmin[i];
  for (k = kn; k <= kx; k++, kp++) {
    pmx->mmx[0][kp] = -eslINFINITY;
    pmx->imx[0][kp] = -eslINFINITY;
    pmx->dmx[0][kp] = fmx->dmx[0][kp] + bmx->dmx[0][kp] - sc;
    if((mass_d[k] = pmx->dmx[0][kp]) > thresh) {
      cp9b->pn_min_d[k] = 0;
      _nset_d[k] = TRUE;
    }
  }

  for (i = 1; i <= L; i++) {
      k = 0;
      if(INBAND(i,0)) {
	kp = k - kmin[i];
	assert(kp == 0);
	pmx->mmx[i][kp] = ESL_MAX(fmx->mmx[i][kp] + bmx->mmx[i][kp] - sc, -eslINFINITY);
	if(! _nset_m[0]) {
	  if((mass_m[0] = p7_FLogsum(mass_m[0], pmx->mmx[i][kp])) > thresh) {
	    cp9b->pn_min_m[0] = i;
	    _nset_m[0] = TRUE;
	  }
	}
	pmx->imx[i][kp] = ESL_MAX(fmx->imx[i][kp] + bmx->imx[i][kp] - Scorify(hmm->isc[dsq[i]][0]) - sc, -eslINFINITY);
	if(! _nset_i[0]) {
	  if((mass_i[0] = p7_FLogsum(mass_i[0], pmx->imx[i][kp])) > thresh) {
	    cp9b->pn_min_i[0] = i;
	    _nset_i[0] = TRUE;
	  }
	}
	pmx->dmx[i][kp] = -eslINFINITY;
      }

      kn = ESL_MAX(kmin[i], 1);
      kx = kmax[i];
      kp = kn - kmin[i];
      for(k = kn; k <= kx; k++, kp++) {
	pmx->mmx[i][kp] = ESL_MAX(fmx->mmx[i][kp] + bmx->mmx[i][kp] - Scorify(hmm->msc[dsq[i]][k]) - sc, -eslINFINITY);
	pmx->imx[i][kp] = ESL_MAX(fmx->imx[i][kp] + bmx->imx[i][kp] - Scorify(hmm->isc[dsq[i]][k]) - sc, -eslINFINITY);
	pmx->dmx[i][kp] = ESL_MAX(fmx->dmx[i][kp] + bmx->dmx[i][kp] - sc, -eslINFINITY);

	if(! _nset_m[k]) {
	  if((mass_m[k] = p7_FLogsum(mass_m[k], pmx->mmx[i][kp])) > thresh) {
	    cp9b->pn_min_m[k] = i;
	    _nset_m[k] = TRUE;
	  }
	}
	if(! _nset_i[k]) {
	  if((mass_i[k] = p7_FLogsum(mass_i[k], pmx->imx[i][kp])) > thresh) {
	    cp9b->pn_min_i[k] = i;
	    _nset_i[k] = TRUE;
	  }
	}
	if(! _nset_d[k]) {
	  if((mass_d[k] = p7_FLogsum(mass_d[k], pmx->dmx[i][kp])) > thresh) {
	    cp9b->pn_min_d[k] = i;
	    _nset_d[k] = TRUE;
	  }
	}
      }
  }

  esl_vec_FSet(mass_m, M+1, -eslINFINITY);
  esl_vec_FSet(mass_i, M+1, -eslINFINITY);
  esl_vec_FSet(mass_d, M+1, -eslINFINITY);
  for (i = L; i >= 1; i--) {
      kp = 0;
      for(k = kmin[i]; k <= kmax[i]; k++, kp++) {
	if(! _xset_m[k]) {
	  if((mass_m[k] = p7_FLogsum(mass_m[k], pmx->mmx[i][kp])) > thresh) {
	    cp9b->pn_max_m[k] = i;
	    _xset_m[k] = TRUE;
	  }
	}
	if(! _xset_i[k]) {
	  if((mass_i[k] = p7_FLogsum(mass_i[k], pmx->imx[i][kp])) > thresh) {
	    cp9b->pn_max_i[k] = i;
	    _xset_i[k] = TRUE;
	  }
	}
	if(! _xset_d[k]) {
	  if((mass_d[k] = p7_FLogsum(mass_d[k], pmx->dmx[i][kp])) > thresh) {
	    cp9b->pn_max_d[k] = i;
	    _xset_d[k] = TRUE;
	  }
	}
      }
  }
  if(INBAND(0,0)) {
    assert(kmin[0] == 0);
    if(! _xset_m[0]) {
      if((mass_m[0] = p7_FLogsum(mass_m[0], pmx->mmx[0][0])) > thresh) {
	cp9b->pn_max_m[0] = 0;
	_xset_m[0] = TRUE;
      }
    }
  }
  kn = ESL_MAX(kmin[i], 1);
  kx = kmax[i];
  kp = kn - kmin[i];
  for(k = kn; k <= kx; k++, kp++) {
    if(!_xset_d[k]) {
      if((mass_d[k] = p7_FLogsum(mass_d[k], pmx->dmx[0][kp])) > thresh) {
	cp9b->pn_max_d[k] = 0;
	_xset_d[k] = TRUE;
      }
    }
  }

  int mset, dset;
  for(k = 0; k <= M; k++) {
    mset = dset = TRUE;
    if(((! _nset_m[k])) || (! _xset_m[k]) || (cp9b->pn_max_m[k] < cp9b->pn_min_m[k])) {
      cp9b->pn_min_m[k] = cp9b->pn_max_m[k] = -1;
      mset = FALSE;
    }
    if(((! _nset_i[k])) || (! _xset_i[k]) || (cp9b->pn_max_i[k] < cp9b->pn_min_i[k])) {
      cp9b->pn_min_i[k] = cp9b->pn_max_i[k] = -1;
    }
    if(((! _nset_d[k])) || (! _xset_d[k]) || (cp9b->pn_max_d[k] < cp9b->pn_min_d[k])) {
      cp9b->pn_min_d[k] = cp9b->pn_max_d[k] = -1;
      dset = FALSE;
    }
    if((!hmm_is_localized) && (mset == FALSE && dset == FALSE)) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "node: %d match nor delete HMM state bands were set in non-localized, non-scanning HMM, lower tau (should be << 0.5).\n", k);
  }

  cp9b->pn_min_d[0] = -1;
  cp9b->pn_max_d[0] = -1;

  if(do_pnmono) pn_match_bands_enforce_monotone(cp9b->pn_min_m, cp9b->pn_max_m, M, L, do_pnmono_print, "fb2hmm_p7bf");

  if(debug_level > 0) cp9_DebugPrintHMMBands(stdout, L, cp9b, (1.-p_thresh), 1);

  free(mass_m); free(mass_i); free(mass_d);
  free(_nset_m); free(_nset_i); free(_nset_d);
  free(_xset_m); free(_xset_i); free(_xset_d);

  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "Memory allocation error.\n");
}


/* Function: cp9_PredictStartAndEndPositionsP7BF()
 *
 * Float-precision mirror of cp9_PredictStartAndEndPositionsP7B(). Reads
 * CP9_FMX posteriors. Per-node pocc[k] is summed in float space directly
 * (no Score2Prob detour through int log space). No median renormalization
 * — the float F/B is supposed to give pocc ~1.0 in glocal mode without a
 * heuristic scale factor, so the do_renorm hack from the int variant is
 * dropped here.
 */
void
cp9_PredictStartAndEndPositionsP7BF(CP9_FMX *pmx, CP9Bands_t *cp9b, int *kmin, int *kmax, int i0, int j0)
{
  int    i, k;
  int    L = j0-i0+1;
  float  pocc;
  float *pocc_arr;

  pocc_arr = malloc(sizeof(float) * (cp9b->hmm_M + 1));
  if(pocc_arr == NULL) cm_Fail("cp9_PredictStartAndEndPositionsP7BF(): malloc failed for pocc_arr");
  for(k = 0; k <= cp9b->hmm_M; k++) pocc_arr[k] = -1.0;
  for(k = 1; k <= cp9b->hmm_M; k++) {
    if(cp9b->pn_min_m[k] == -1 && cp9b->pn_min_i[k] == -1 && cp9b->pn_min_d[k] == -1) continue;
    pocc = 0.0;
    for(i = 0; i <= L; i++) {
      if(k >= kmin[i] && k <= kmax[i]) {
	int kp = k - kmin[i];
	pocc += expf(pmx->mmx[i][kp]);
	pocc += expf(pmx->dmx[i][kp]);
      }
    }
    pocc_arr[k] = pocc;
  }

  /* Part 1: sp1/sp2 — leftmost nodes with significant occupancy. */
  k = 1;
  cp9b->sp1 = cp9b->sp2 = -1;
  while(k <= cp9b->hmm_M && (cp9b->sp1 == -1 || cp9b->sp2 == -1)) {
    if(pocc_arr[k] < 0.0) { k++; }
    else {
      pocc = pocc_arr[k];
      if((cp9b->sp1 == -1) && (pocc > cp9b->thresh1)) cp9b->sp1 = k;
      if((cp9b->sp2 == -1) && (pocc > cp9b->thresh2)) cp9b->sp2 = k;
      k++;
    }
  }
  if(k == cp9b->hmm_M+1) {
    if(cp9b->sp1 == -1) { cp9b->sp1 = cp9b->hmm_M+1; }
    if(cp9b->sp2 == -1) { cp9b->sp2 = cp9b->hmm_M+1; }
  }

  /* Part 2: ep1/ep2 — rightmost nodes with significant occupancy. */
  if((cp9b->sp1 == cp9b->hmm_M+1) &&
     (cp9b->sp2 == cp9b->hmm_M+1)) {
    cp9b->ep1 = 0;
    cp9b->ep2 = 0;
  }
  else {
    cp9b->ep1 = cp9b->ep2 = -1;
    k = cp9b->hmm_M;
    while(k >= 1 && (cp9b->ep1 == -1 || cp9b->ep2 == -1)) {
      if(pocc_arr[k] < 0.0) { k--; }
      else {
	pocc = pocc_arr[k];
	if((cp9b->ep1 == -1) && (pocc > cp9b->thresh1)) cp9b->ep1 = k;
	if((cp9b->ep2 == -1) && (pocc > cp9b->thresh2)) cp9b->ep2 = k;
	k--;
      }
    }
    if(k == 0) {
      if(cp9b->ep1 == -1) { cp9b->ep1 = 0; }
      if(cp9b->ep2 == -1) { cp9b->ep2 = 0; }
    }
  }
  free(pocc_arr);

  /* Parts 3-4: Rmarg/Lmarg derivation — identical to int variant. */

  /* Rmarg_imin */
  if(cp9b->sp1 == cp9b->hmm_M+1) { cp9b->Rmarg_imin = i0; }
  else {
    cp9b->Rmarg_imin = INT_MAX;
    if(cp9b->sp1 != (cp9b->hmm_M+1) && cp9b->pn_min_m[cp9b->sp1] >= 0) cp9b->Rmarg_imin = ESL_MIN(cp9b->Rmarg_imin, cp9b->pn_min_m[cp9b->sp1]);
    if(cp9b->sp1 != (cp9b->hmm_M+1) && cp9b->pn_min_i[cp9b->sp1] >= 0) cp9b->Rmarg_imin = ESL_MIN(cp9b->Rmarg_imin, cp9b->pn_min_i[cp9b->sp1]);
    if(cp9b->sp1 != (cp9b->hmm_M+1) && cp9b->pn_min_d[cp9b->sp1] >= 0) cp9b->Rmarg_imin = ESL_MIN(cp9b->Rmarg_imin, cp9b->pn_min_d[cp9b->sp1]);
    if(cp9b->sp2 != (cp9b->hmm_M+1) && cp9b->pn_min_m[cp9b->sp2] >= 0) cp9b->Rmarg_imin = ESL_MIN(cp9b->Rmarg_imin, cp9b->pn_min_m[cp9b->sp2]);
    if(cp9b->sp2 != (cp9b->hmm_M+1) && cp9b->pn_min_i[cp9b->sp2] >= 0) cp9b->Rmarg_imin = ESL_MIN(cp9b->Rmarg_imin, cp9b->pn_min_i[cp9b->sp2]);
    if(cp9b->sp2 != (cp9b->hmm_M+1) && cp9b->pn_min_d[cp9b->sp2] >= 0) cp9b->Rmarg_imin = ESL_MIN(cp9b->Rmarg_imin, cp9b->pn_min_d[cp9b->sp2]);
    if(cp9b->Rmarg_imin == INT_MAX || cp9b->sp1 == (cp9b->hmm_M+1) || cp9b->sp2 == (cp9b->hmm_M+1)) cp9b->Rmarg_imin = i0;
    cp9b->Rmarg_imin = ESL_MAX(i0,   cp9b->Rmarg_imin);
    cp9b->Rmarg_imin = ESL_MIN(j0+1, cp9b->Rmarg_imin);
  }
  /* Rmarg_imax */
  if(cp9b->sp1 == cp9b->hmm_M+1) { cp9b->Rmarg_imax = j0; }
  else {
    cp9b->Rmarg_imax = INT_MIN;
    if(cp9b->sp1 != (cp9b->hmm_M+1) && cp9b->pn_max_m[cp9b->sp1] >= 0) cp9b->Rmarg_imax = ESL_MAX(cp9b->Rmarg_imax, cp9b->pn_max_m[cp9b->sp1]);
    if(cp9b->sp1 != (cp9b->hmm_M+1) && cp9b->pn_max_i[cp9b->sp1] >= 0) cp9b->Rmarg_imax = ESL_MAX(cp9b->Rmarg_imax, cp9b->pn_max_i[cp9b->sp1]);
    if(cp9b->sp1 != (cp9b->hmm_M+1) && cp9b->pn_max_d[cp9b->sp1] >= 0) cp9b->Rmarg_imax = ESL_MAX(cp9b->Rmarg_imax, cp9b->pn_max_d[cp9b->sp1]);
    if(cp9b->sp2 != (cp9b->hmm_M+1) && cp9b->pn_max_m[cp9b->sp2] >= 0) cp9b->Rmarg_imax = ESL_MAX(cp9b->Rmarg_imax, cp9b->pn_max_m[cp9b->sp2]);
    if(cp9b->sp2 != (cp9b->hmm_M+1) && cp9b->pn_max_i[cp9b->sp2] >= 0) cp9b->Rmarg_imax = ESL_MAX(cp9b->Rmarg_imax, cp9b->pn_max_i[cp9b->sp2]);
    if(cp9b->sp2 != (cp9b->hmm_M+1) && cp9b->pn_max_d[cp9b->sp2] >= 0) cp9b->Rmarg_imax = ESL_MAX(cp9b->Rmarg_imax, cp9b->pn_max_d[cp9b->sp2]);
    if(cp9b->Rmarg_imax == INT_MIN || cp9b->sp1 == (cp9b->hmm_M+1) || cp9b->sp2 == (cp9b->hmm_M+1)) cp9b->Rmarg_imax = j0+1;
    cp9b->Rmarg_imax = ESL_MAX(i0,   cp9b->Rmarg_imax);
    cp9b->Rmarg_imax = ESL_MIN(j0+1, cp9b->Rmarg_imax);
  }
  /* Lmarg_jmin */
  if(cp9b->ep1 == 0) { cp9b->Lmarg_jmin = i0-1; }
  else {
    cp9b->Lmarg_jmin = INT_MAX;
    if(cp9b->ep1 != 0 && cp9b->pn_min_m[cp9b->ep1] >= 0) cp9b->Lmarg_jmin = ESL_MIN(cp9b->Lmarg_jmin, cp9b->pn_min_m[cp9b->ep1]);
    if(cp9b->ep1 != 0 && cp9b->pn_min_i[cp9b->ep1] >= 0) cp9b->Lmarg_jmin = ESL_MIN(cp9b->Lmarg_jmin, cp9b->pn_min_i[cp9b->ep1]);
    if(cp9b->ep1 != 0 && cp9b->pn_min_d[cp9b->ep1] >= 0) cp9b->Lmarg_jmin = ESL_MIN(cp9b->Lmarg_jmin, cp9b->pn_min_d[cp9b->ep1]-1);
    if(cp9b->ep2 != 0 && cp9b->pn_min_m[cp9b->ep2] >= 0) cp9b->Lmarg_jmin = ESL_MIN(cp9b->Lmarg_jmin, cp9b->pn_min_m[cp9b->ep2]);
    if(cp9b->ep2 != 0 && cp9b->pn_min_i[cp9b->ep2] >= 0) cp9b->Lmarg_jmin = ESL_MIN(cp9b->Lmarg_jmin, cp9b->pn_min_i[cp9b->ep2]);
    if(cp9b->ep2 != 0 && cp9b->pn_min_d[cp9b->ep2] >= 0) cp9b->Lmarg_jmin = ESL_MIN(cp9b->Lmarg_jmin, cp9b->pn_min_d[cp9b->ep2]-1);
    if(cp9b->Lmarg_jmin == INT_MAX || cp9b->ep1 == 0 || cp9b->ep2 == 0) cp9b->Lmarg_jmin = i0-1;
    cp9b->Lmarg_jmin = ESL_MAX(i0-1, cp9b->Lmarg_jmin);
    cp9b->Lmarg_jmin = ESL_MIN(j0,   cp9b->Lmarg_jmin);
  }
  /* Lmarg_jmax */
  if(cp9b->ep1 == 0) { cp9b->Lmarg_jmax = j0; }
  else {
    cp9b->Lmarg_jmax = INT_MIN;
    if(cp9b->ep1 != 0 && cp9b->pn_max_m[cp9b->ep1] >= 0) cp9b->Lmarg_jmax = ESL_MAX(cp9b->Lmarg_jmax, cp9b->pn_max_m[cp9b->ep1]);
    if(cp9b->ep1 != 0 && cp9b->pn_max_i[cp9b->ep1] >= 0) cp9b->Lmarg_jmax = ESL_MAX(cp9b->Lmarg_jmax, cp9b->pn_max_i[cp9b->ep1]);
    if(cp9b->ep1 != 0 && cp9b->pn_max_d[cp9b->ep1] >= 0) cp9b->Lmarg_jmax = ESL_MAX(cp9b->Lmarg_jmax, cp9b->pn_max_d[cp9b->ep1]-1);
    if(cp9b->ep2 != 0 && cp9b->pn_max_m[cp9b->ep2] >= 0) cp9b->Lmarg_jmax = ESL_MAX(cp9b->Lmarg_jmax, cp9b->pn_max_m[cp9b->ep2]);
    if(cp9b->ep2 != 0 && cp9b->pn_max_i[cp9b->ep2] >= 0) cp9b->Lmarg_jmax = ESL_MAX(cp9b->Lmarg_jmax, cp9b->pn_max_i[cp9b->ep2]);
    if(cp9b->ep2 != 0 && cp9b->pn_max_d[cp9b->ep2] >= 0) cp9b->Lmarg_jmax = ESL_MAX(cp9b->Lmarg_jmax, cp9b->pn_max_d[cp9b->ep2]-1);
    if(cp9b->Lmarg_jmax == INT_MIN || cp9b->ep1 == 0 || cp9b->ep2 == 0) cp9b->Lmarg_jmax = j0;
    cp9b->Lmarg_jmax = ESL_MAX(i0-1, cp9b->Lmarg_jmax);
    cp9b->Lmarg_jmax = ESL_MIN(j0,   cp9b->Lmarg_jmax);
  }

  return;
}


/* Function: cp9_FBMatrices2BandsF()
 *
 * Float-precision mirror of cp9_FBMatrices2BandsP7B(). Phase 2 of the
 * P7-banded CP9 band derivation: given filled CP9_FMX fmx/bmx, compute
 * pn_min/pn_max bands via cp9_FB2HMMBandsP7BF, then convert to CM bands.
 * Used only by cp9_IterateSeq2BandsP7B when do_trunc=TRUE.
 */
int
cp9_FBMatrices2BandsF(CM_t *cm, char *errbuf, CP9_t *cp9, CP9_FMX *fmx, CP9_FMX *bmx, CP9_FMX *pmx,
		     ESL_DSQ *dsq, CP9Bands_t *cp9b, int *kmin, int *kmax,
		     int L, int i0, int j0, int pass_idx, int debug_level,
		     int do_pnmono, int do_pnmono_print)
{
  int status;
  int use_sums      = ((cm->align_opts & CM_ALIGN_SUMS) || (cm->search_opts & CM_SEARCH_SUMS)) ? TRUE : FALSE;
  int do_old_hmm2ij = ((cm->align_opts & CM_ALIGN_HMM2IJOLD) || (cm->search_opts & CM_SEARCH_HMM2IJOLD)) ? TRUE : FALSE;
  int do_trunc      = cm_pli_PassAllowsTruncation(pass_idx);

  if(use_sums) {
    printf("USE SUMS!\n");
    exit(1);
  }
  else {
    if((status = cp9_FB2HMMBandsP7BF(cp9, errbuf, dsq, fmx, bmx, pmx, cp9b, L, cp9b->hmm_M,
				     (1.-cm->tau), do_old_hmm2ij, kmin, kmax, debug_level,
				     do_pnmono, do_pnmono_print)) != eslOK) return status;
    cp9b->tau = cm->tau;
  }
  if(debug_level > 0) cp9_DebugPrintHMMBands(stdout, L, cp9b, cm->tau, 1);

  /* Shift HMM bands from 1..L to i0..j0 */
  if(i0 != 1) {
    int offset = i0 - 1;
    int k;
    for(k = 0; k <= cp9b->hmm_M; k++) {
      if(cp9b->pn_min_m[k] != -1) { cp9b->pn_min_m[k] += offset; cp9b->pn_max_m[k] += offset; }
      if(cp9b->pn_min_i[k] != -1) { cp9b->pn_min_i[k] += offset; cp9b->pn_max_i[k] += offset; }
      if(cp9b->pn_min_d[k] != -1) { cp9b->pn_min_d[k] += offset; cp9b->pn_max_d[k] += offset; }
    }
  }

  if(do_trunc) {
    /* float path: do_renorm always FALSE. The point of the float DP is that
     * pocc[k] arrives at ~1.0 in glocal mode without the median heuristic. */
    cp9_PredictStartAndEndPositionsP7BF(pmx, cp9b, kmin, kmax, i0, j0);
    /* brief 26_0430-149: in glocal alignment the full (J-mode) parse must always be
     * geometrically available. The thresh1 escalation in cp9_IterateSeq2BandsP7B
     * can retreat ep1 below clen (and, in principle, push sp1 above 1) on models
     * with a decaying posterior-occupancy tail -- e.g. pure-MATL VADR genome
     * models such as NC_001959, where every state's rpos == clen, so once
     * ep1 < clen every state gets Jvalid[v] = FALSE in
     * cp9_MarginalCandidatesFromStartEndPositions(). That excludes the full parse
     * and cm_TrInsideAlignHB() returns "no valid parsetree" in -g mode (local mode
     * survives via an EL escape). Floor sp1 <= 1 and ep1 >= clen so the entire
     * model stays J-valid. Scoped to glocal (CMH_LOCAL_BEGIN off) to leave local
     * mode byte-identical (it already has a valid root, and EL handles the tail).
     * This whole float-truncated band path is cmalign-only: cp9_IterateSeq2BandsP7B
     * is reached only from cm_alndata.c with doing_search == FALSE, so
     * (do_align && !doing_search) holds by construction. */
    if(! (cm->flags & CMH_LOCAL_BEGIN)) {
      if(cp9b->sp1 > 1)        cp9b->sp1 = 1;
      if(cp9b->ep1 < cm->clen) cp9b->ep1 = cm->clen;
    }
    if((status = cp9_MarginalCandidatesFromStartEndPositions(cm, cp9b, pass_idx, errbuf)) != eslOK) return status;
  }
  else {
    esl_vec_ISet(cp9b->Jvalid, cm->M+1, TRUE);
    esl_vec_ISet(cp9b->Lvalid, cm->M+1, FALSE);
    esl_vec_ISet(cp9b->Rvalid, cm->M+1, FALSE);
    esl_vec_ISet(cp9b->Tvalid, cm->M+1, FALSE);
  }

  if(do_old_hmm2ij) {
    /* brief 26_0430-148 (lands brief 26_0430-082): pass doing_search=FALSE so cp9_HMM2ijBands_OLD applies the
     * global-alignment ROOT_S span enforcement (hmmband.c). This float-truncated band path is
     * cmalign-only / do_trunc-only; without the pin the full J-mode parse is geometrically
     * excluded in -g mode -> cm_TrInsideAlignHB returns "no valid parsetree". */
    if((status = cp9_HMM2ijBands_OLD(cm, errbuf, cm->cp9b, cm->cp9map, i0, j0, FALSE, debug_level)) != eslOK) return status;
  }
  else {
    /* brief 26_0430-148 (lands brief 26_0430-082): doing_search=FALSE (see comment above). */
    if((status = cp9_HMM2ijBands(cm, errbuf, cp9, cm->cp9b, cm->cp9map, i0, j0, FALSE, do_trunc, debug_level)) != eslOK) return status;
  }
  if((status = cp9_GrowHDBands(cp9b, errbuf)) != eslOK) return status;
  ij2d_bands(cm, cp9b, do_trunc, debug_level);

  if(do_trunc && (! (cm->flags & CMH_LOCAL_BEGIN))) {
    /* brief 26_0430-185: brief 26_0430-149's sp1/ep1 floor (above) forces Jvalid[v] = TRUE for
     * essentially every state so the glocal J-parse stays geometrically available
     * for models with a decaying-but-real occupancy tail (its target case, e.g.
     * NC_001959). It does not widen the real per-state (j,d) bands, computed just
     * above by ij2d_bands() from the un-floored 1-tau threshold signal -- so on a
     * sequence with a genuinely, biologically missing region the floor also
     * flags states whose real band is empty as "J-valid": geometrically
     * unreachable "phantom valid" states that a truncated-alignment traceback
     * can walk into and die on (cm_TrInsideAlignHB() "no valid parsetree
     * found"). Veto Jvalid[v] back to FALSE for any state whose real band is
     * empty at every j in its jband, using the hd_min()/hd_max()
     * recompute-on-demand accessors (brief 26_0430-157) -- never reintroduce flat
     * hdmin[v][]/hdmax[v][] reads here, they're gone. */
    int v, jp, njp, found;
    for(v = 0; v < cp9b->cm_M; v++) {
      if(! cp9b->Jvalid[v]) continue;
      njp = cp9b->jmax[v] - cp9b->jmin[v] + 1;
      found = FALSE;
      for(jp = 0; jp < njp; jp++) {
        if(hd_min(cp9b, v, jp) <= hd_max(cp9b, v, jp)) { found = TRUE; break; }
      }
      if(! found) cp9b->Jvalid[v] = FALSE;
    }
  }

#if eslDEBUGLEVEL >= 1
  if((status = cp9_ValidateBands(cm, errbuf, cp9b, i0, j0, do_trunc)) != eslOK) return status;
#endif
  if(debug_level > 0) debug_print_ij_bands(cm);
  if(debug_level > 0) PrintDPCellsSaved_jd(cm, cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax, L);

  return eslOK;
}


/* Function: p7_Seq2Bands
 * Date    : EPN, Fri Aug 15 14:29:01 2008
 *
 * Purpose:  Given a CM with a valid Plan 7 HMM - run the MSV algorithm
 *           to determine bands to be used on the CP9 HMM parse.
 *           
 * Args:     cm          - the covariance model
 *           errbuf      - char buffer for reporting errors
 *           P7_PROFILE  - generic profile
 *           P7_GMX      - generic P7 dp matrix
 *           P7_BG       - P7 null model 
 *           P7_TR       - P7 trace
 *           dsq         - sequence in digitized form
 *           L           - length of sequence we're aligning (1..L)
 *           phi         - phi array, phi[k][v] is expected number of times (probability)
 *                         state v (0 = match, 1 insert, 2 = delete) in 
 *                         node k is *entered*. Node 0 is special, state 0 = B state, state 1 = N_state, state 2 = NULL
 *                         Calculated *without* taking insert->insert transitions into account.
 *           sc7         - minimum score to allow as a pin, 0. to allow any score
 *           len7        - min n-mer size, 1 to allow any size
 *           end7        - min distance from end to allow, prune away any others, 0 to not prune based on end proximity
 *           mprob7      - min match phi probability to allow in a pin 
 *           mcprob7     - min cumulative match phi probability to allow in a nmer pin
 *           iprob7      - max insert phi probability to allow in a pin 
 *           ilprob7     - max insert phi probability to allow in a state to the left of a pin 
 *           ret_i2k     - [0.i..L] = k, residue i emitted from node k's match state in MSV trace 
 *           ret_kmin    - [0.i..L] = k, min node k for residue i
 *           ret_kmax    - [0.i..L] = k, max node k for residue i
 *           ret_ncells  - number of cells within bands, to return
 *
 * Return:  eslOK on success;
 * 
 */
int
p7_Seq2Bands(CM_t *cm, char *errbuf, P7_PROFILE *gm, P7_GMX *gx, P7_BG *bg, P7_TRACE *p7_tr, ESL_DSQ *dsq, int L,
	     double **phi, float sc7, int len7, int end7, float mprob7, float mcprob7, float iprob7, float ilprob7, int pad7,
	     int **ret_i2k, int **ret_kmin, int **ret_kmax, int *ret_ncells)
{
  int   status;
  float usc, nullsc;
  int *k2i, *i2k;
  float *isc;
  int *iconflict;
  int *kmin, *kmax;
  int ncells;
  int M = gm->M;  /* model length; equals cm->mlp7->M and cm->clen for RNA profiles */
  ESL_STOPWATCH *s2b_watch = NULL;

  /* setup for all modes */
  p7_bg_SetLength(bg, L);
  p7_bg_NullOne(bg, dsq, L, &nullsc);

  /* generic mode setup */
  p7_gmx_GrowTo(gx, M, L);
  p7_ReconfigLength(gm, L);
  gx->M = M;
  gx->L = L;

  s2b_watch = esl_stopwatch_Create();

  /* Step 1: GMSV algorithm */
  esl_stopwatch_Start(s2b_watch);
  p7_GMSV(dsq, L, gm, gx, 2.0, &usc);
  esl_stopwatch_Stop(s2b_watch);

  /* Step 2: traceback MSV */
  esl_stopwatch_Start(s2b_watch);
  status = my_p7_GTraceMSV(dsq, L, gm, gx, p7_tr, &i2k, &k2i, &isc, &iconflict);
  esl_stopwatch_Stop(s2b_watch);

  /* Debug: print pins with isc scores before pruning */
  if(status == eslOK && getenv("MSVBAND_DBG")) {
    int dbg_i;
    for(dbg_i = 1; dbg_i <= L; dbg_i++) {
      if(i2k[dbg_i] != -1) {
        fprintf(stderr, "#MSVBAND_DBG: isc pin i=%5d k=%4d isc=%7.4f\n", dbg_i, i2k[dbg_i], isc[dbg_i]);
      }
    }
  }

  /* Step 3: prune pins */
  esl_stopwatch_Start(s2b_watch);
  if(status == eslOK) { /* trace is valid */
    prune_i2k(i2k, iconflict, isc, L, phi, sc7, len7, end7, mprob7, mcprob7, iprob7, ilprob7);
  }
  else if (status == eslEINCOMPAT) { /* trace was discontiguous, abort! remove all pins */
    esl_vec_ISet(k2i, (M+1), -1);
    esl_vec_ISet(i2k, (L+1), -1);
    esl_stopwatch_Destroy(s2b_watch);
    return status;
  }
  esl_stopwatch_Stop(s2b_watch);

  /* Step 4: pins -> bands */
  esl_stopwatch_Start(s2b_watch);
  if((status = p7_pins2bands(i2k, errbuf, L, M, pad7, &kmin, &kmax, &ncells)) != eslOK) { esl_stopwatch_Destroy(s2b_watch); return status; }
  esl_stopwatch_Stop(s2b_watch);
  esl_stopwatch_Destroy(s2b_watch);
  /*DumpP7Bands(stdout, i2k, kmin, kmax, L); */

  /* print gmx in heatmap format */ 
  /*
    ESL_DMATRIX *D;
    double min;
    double max;
    FILE *hfp;
    p7_gmx_Match2DMatrix(gx, TRUE, &D, &min, &max);
    hfp = fopen("cur.ps", "w");
    my_dmx_Visualize(hfp, D, 0.01, max, 0.01);
    fclose(hfp);
    esl_dmatrix_Destroy(D);
    */

  *ret_i2k  = i2k;
  *ret_kmin = kmin;
  *ret_kmax = kmax;
  *ret_ncells = ncells;

  free(iconflict);
  free(isc);
  free(k2i);

  return eslOK;
}


/* Function: p7_Seq2BandsVit()
 * Incept:   EPN*
 *
 * Synopsis: Derive p7 bands from glocal Viterbi alignment.
 *
 * Purpose:  Given a profile <gm> already configured in the desired mode
 *           (glocal, local, or truncated-glocal) and a digital sequence
 *           <dsq> of length <L>, run p7_GViterbi(), traceback the optimal
 *           alignment with p7_GTrace(), and derive per-residue kmin/kmax
 *           bands suitable for banded glocal Forward.
 *
 *           Unlike p7_Seq2Bands() (which uses GMSV + my_p7_GTraceMSV),
 *           the Viterbi trace is globally optimal and free of diagonal
 *           conflicts, so no pin-conflict pruning is performed.
 *
 *           Profile mode is the caller's responsibility: configure gm
 *           to p7_GLOCAL (or truncated equivalent) before calling, and
 *           restore it afterward if needed.
 *
 * Args:     errbuf    - for error messages
 *           gm        - profile, configured in desired mode by caller
 *           gx        - generic DP matrix (will be grown as needed)
 *           bg        - null model (unused; reserved for future)
 *           p7_tr     - pre-allocated trace (will be reused)
 *           dsq       - digital sequence, 1..L
 *           L         - length of dsq
 *           pad       - band half-width passed to p7_pins2bands()
 *           nodepad   - per-node pad array [0..M], or NULL to use pad
 *           hopback   - if >0, dilate band per pinned residue by min/max of
 *                       i2k across the 2*hopback+1-pin trace-order window.
 *                       Only applied when nodepad path is used (pin set is
 *                       monotone). 0 = off.
 *           vitend    - if >0, drop the first <vitend> and last <vitend>
 *                       Vit pins (M-state entries in i2k) before pins->bands.
 *                       Targets Mode-1 boundary truncation by widening the
 *                       band at the prefix/suffix to [0,M].
 *           ret_i2k   - RETURN: per-residue pin array (caller frees)
 *           ret_kmin  - RETURN: per-residue kmin array (caller frees)
 *           ret_kmax  - RETURN: per-residue kmax array (caller frees)
 *           ret_ncells- RETURN: total banded cells; 0 = fall back to unbanded
 *
 * Returns:  eslOK on success.
 *           eslFAIL if Viterbi finds no valid path; ret_* set to NULL,
 *           ret_ncells set to 0; caller should fall back to unbanded Forward.
 */
int
p7_Seq2BandsVit(char *errbuf, P7_PROFILE *gm, P7_GMX *gx, P7_BG *bg, P7_TRACE *p7_tr,
		ESL_DSQ *dsq, int L, int pad, int *nodepad, int hopback, int vitend,
		int **ret_i2k, int **ret_kmin, int **ret_kmax, int *ret_ncells)
{
  int    status;
  float  sc;
  int   *i2k  = NULL;
  int   *kmin = NULL;
  int   *kmax = NULL;
  int    ncells = 0;
  int    M = gm->M;
  int    tpos;

  /* Setup DP matrix */
  p7_gmx_GrowTo(gx, M, L);
  gx->M = M;
  gx->L = L;

  /* Step 1: Viterbi DP */
  if ((status = p7_GViterbi(dsq, L, gm, gx, &sc)) != eslOK)
    ESL_FAIL(status, errbuf, "p7_GViterbi() failed in p7_Seq2BandsVit()");

  /* Step 2: Traceback.
   * p7_GTrace() expects an empty trace (tr->N == 0); call Reuse() first.
   * If no valid path (eslFAIL), return eslOK with ncells=0 so the
   * caller falls back to unbanded Forward.
   */
  p7_trace_Reuse(p7_tr);
  status = p7_GTrace(dsq, L, gm, gx, p7_tr);
  if (status == eslFAIL) {
    /* Empty trace: no valid Viterbi path. Signal caller to fall back. */
    *ret_i2k    = NULL;
    *ret_kmin   = NULL;
    *ret_kmax   = NULL;
    *ret_ncells = 0;
    return eslOK;
  }
  if (status != eslOK) ESL_FAIL(status, errbuf, "p7_GTrace() failed in p7_Seq2BandsVit()");

  /* Step 3: Populate i2k from M states in trace.
   * Viterbi trace is non-conflicting by construction; no conflict handling needed.
   * For multi-hit (glocal multihit mode via J state), both hits populate i2k
   * at distinct sequence positions — no conflict.
   */
  ESL_ALLOC(i2k, sizeof(int) * (L + 1));
  esl_vec_ISet(i2k, (L + 1), -1);

  int nB = 0; /* count B states = number of hits in trace */
  int nM = 0; /* total M states */
  for (tpos = 0; tpos < p7_tr->N; tpos++) {
    if (p7_tr->st[tpos] == p7T_B) nB++;
    if (p7_tr->st[tpos] == p7T_M) {
      nM++;
      int i = p7_tr->i[tpos];
      int k = p7_tr->k[tpos];
      if (i >= 1 && i <= L && k >= 1 && k <= M)
	i2k[i] = k;
    }
  }
  if (getenv("VITBAND_MULTIHIT_DBG")) {
    fprintf(stderr, "#VITBAND_MULTIHIT L=%d M=%d nhit=%d nM=%d sc=%.2f\n", L, M, nB, nM, sc);
  }

  /* Step 3b: Optional vitend pruning. Drop the first <vitend> and last
   * <vitend> populated entries of i2k (walking forward, then backward,
   * over non--1 entries), setting them to -1. This widens the band at
   * the prefix/suffix to [0,M] (no pin -> no per-pin width contribution
   * and the bridge pass between consecutive pins doesn't reach those
   * edge rows). Targets Mode-1 boundary truncation: with end pins
   * removed, F5 backward can integrate posterior mass into the prefix
   * and suffix that would otherwise be locked to nodepad[k_first/last].
   * Skip if nM <= 2*vitend (would clear all pins; falls back to unbanded
   * via empty pin set, behaviour matches no-pin case).
   */
  if (vitend > 0 && nM > 0) {
    int dropped, i;
    /* Forward pass: drop first <vitend> non--1 entries. */
    dropped = 0;
    for (i = 1; i <= L && dropped < vitend; i++) {
      if (i2k[i] != -1) {
	i2k[i] = -1;
	dropped++;
      }
    }
    /* Backward pass: drop last <vitend> non--1 entries. */
    dropped = 0;
    for (i = L; i >= 1 && dropped < vitend; i--) {
      if (i2k[i] != -1) {
	i2k[i] = -1;
	dropped++;
      }
    }
  }

  /* Step 4: Pins -> bands */
  if (nodepad != NULL) {
    /* brief 26_0628-043: Viterbi-trace path doesn't expose a tunable ramp alpha --
     * pass brief 26_0628-042's validated default unconditionally. */
    if ((status = p7_pins2bands_nodepad(i2k, errbuf, L, M, nodepad, hopback, 0.75, &kmin, &kmax, &ncells)) != eslOK)
      goto ERROR;
  }
  else {
    if ((status = p7_pins2bands(i2k, errbuf, L, M, pad, &kmin, &kmax, &ncells)) != eslOK)
      goto ERROR;
  }

  *ret_i2k    = i2k;
  *ret_kmin   = kmin;
  *ret_kmax   = kmax;
  *ret_ncells = ncells;
  return eslOK;

 ERROR:
  if (i2k)  free(i2k);
  if (kmin) free(kmin);
  if (kmax) free(kmax);
  return status;
}


/**************************************************************
 * Function: CP9NodeForPosnP7B()
 * Incept:   EPN, Tue Aug 19 14:30:51 2008
 * 
 * Purpose:  Given a P7 banded CP9 posterior matrix,
 *           determine the node of the CP9 HMM that is most likely to 
 *           have emitted (from either its Match or Insert state)
 *           a given posn in the target sequence.
 *
 * Args:     hmm       - the CM plan 9 HMM
 *           errbuf    - for error messages
 *           x         - posn of target subsequence we're interested in
 *           L         - last position of target sequence 
 *           post      - the posterior matrix for the hmm
 *           kn        - min node k for residue x
 *           kx        - max node k for residue x
 *           ret_node  - RETURN: index of node with highest probability of emitting x
 *           ret_type  - RETURN: type of state in ret_node with highest probability 
 *           print_flag- TRUE to print out info on most likely node 
 *           
 *
 * Returns:  eslOK on success;
 *           eslEINVAL on contract violation.
 *           eslEINCOMPAT if kmin[x] >= 
 */
int 
CP9NodeForPosnP7B(CP9_t *hmm, char *errbuf, int x, CP9_MX *post, 
		  int kn, int kx, int *ret_node, int *ret_type, int print_flag)
{
  /* post->mmx[i][kp]: posterior probability that posn i was emitted from node k's 
     match state, k = kp + kmin[i] */  
  int  max_k;    /* node index with highest posterior probability of emitting posn x */
  int  max_type; /* type of state in max_k node with max probability '0' for match, 
		    '1' for insert */
  int  max_sc;   /* score (log probability) from post matrix for max_k node max_type state type */
  int  k;        /* counter over nodes */
  int  kp;       /* k': k offset in position x's band */
    
  if(kn > kx) ESL_FAIL(eslEINVAL, errbuf, "ERROR in CP9NodeForPosn(), kn (%d) > kx (%d)\n", kn, kx);
  
  kp = 0;
  k  = kn;
  if(post->mmx[x][0] > post->imx[x][0]) { 
    max_sc     = post->mmx[x][0];
    max_type   = 0; /* match */
  }
  else {
    max_sc     = post->imx[x][0];
    max_type   = 1; /* insert */
  }
  max_k = k; 

  /* move left to right through HMM nodes */
  for(k = kn+1, kp = 1; k <= kx; k++, kp++) {
    if(post->mmx[x][kp] > max_sc) {
      max_k  = k;
      max_sc = post->mmx[x][kp];
      max_type = 0; /* match */
    }
    if(post->imx[x][kp] > max_sc) {
      max_k  = k;
      max_sc = post->imx[x][kp];
      max_type = 1; /* insert */
    }
  }

  if(print_flag) { 
    if(max_type == 0) printf("MATCH  | mx->mmx[%3d][%3d]: %9d | %8f\n", x, max_k, post->mmx[x][max_k-kn],  Score2Prob(post->mmx[x][max_k-kn], 1.));
    else      	      printf("INSERT | mx->imx[%3d][%3d]: %9d | %8f\n", x, max_k, post->imx[x][max_k-kn], Score2Prob(post->imx[x][max_k-kn], 1.));
  }
  *ret_node = max_k;
  *ret_type = max_type;
  return eslOK;
}

/* Function: P7BandsAdjustForSubCM()
 * Incept:   EPN, Tue Aug 19 14:47:55 2008
 * 
 * Purpose:  Correct k bands kmin, kmax built from an original CM for it's sub CM model 
 *           that models spos..epos.
 *           
 * Args:     kmin     - [0.i..L] = k, min node k for residue i
 *           kmax     - [0.i..L] = k, max node k for residue i
 *           L        - length of current sequence
 *           spos     - min k valid in sub CM 
 *           epos     - max k valid in sub CM 
 *
 * Return:   <eslOK> on success.
 *
 */
int
P7BandsAdjustForSubCM(int *kmin, int *kmax, int L, int spos, int epos)
{
  int i;
  int M = epos - spos + 1;
  for(i = 0; i <= L; i++) { 
    kmin[i] = ESL_MAX(kmin[i] - (spos-1), 0);
    kmin[i] = ESL_MIN(kmin[i], M);

    kmax[i] = ESL_MAX(kmax[i] - (spos-1), 0);
    kmax[i] = ESL_MIN(kmax[i], M);
  }
  kmin[0] = 0; /* hard-coded, M_0 is begin state, it must emit full sequence */

  return eslOK;
}

/* Function: p7_kbands2gbands()
 * Date:     EPN*, Sun Mar 16 2026
 * 
 * Purpose:  Convert k-indexed bands (kmin[i], kmax[i]) to HMMER's 
 *           i-indexed P7_GBANDS structure. 
 *           
 *           Only includes positions where i2k[i] != -1 (MSV-aligned positions).
 *           Positions with i2k=-1 are excluded, creating gaps between segments.
 *           This allows unaligned regions to use N/J/C states freely.
 *           
 * Args:     i2k    - [0..i..L] = k or -1, which HMM node aligns to residue i
 *           kmin   - [0..i..L] = k, min node k for residue i
 *           kmax   - [0..i..L] = k, max node k for residue i
 *           L      - length of sequence
 *           M      - length of HMM model
 *           ret_bnd - RETURN: newly allocated P7_GBANDS structure
 *
 * Returns:  <eslOK> on success, <ret_bnd> points to new P7_GBANDS.
 *           <eslEMEM> on allocation failure.
 *           
 * Note:     Creates multiple segments for sparse alignments. Gaps between
 *           segments are handled by special state transitions (N/J/C).
 */
int
p7_kbands2gbands(int *i2k, int *kmin, int *kmax, int L, int M, P7_GBANDS **ret_bnd)
{
  P7_GBANDS *bnd = NULL;
  int        status;
  int        i;
  int        ka, kb;

  if ((bnd = p7_gbands_Create()) == NULL) { status = eslEMEM; goto ERROR; }
  
  /* Include all positions i=1..L, but widen bands for unaligned positions.
   * For positions with i2k=-1 (no MSV alignment), if the inherited band is
   * narrow, widen it to allow more flexibility for D states and special states.
   * Clamp bands to valid range [1..M] - HMM nodes are 1-indexed.
   */
  for (i = 1; i <= L; i++) {
    ka = ESL_MAX(1, kmin[i]);  /* ensure ka >= 1 */
    kb = ESL_MIN(M, kmax[i]);  /* ensure kb <= M */
    if (ka > kb) ka = kb;      /* ensure ka <= kb (can happen with non-uniform nodepad) */
    
    /* (widening for unaligned positions removed — was inflating bands) */
    
    if ((status = p7_gbands_Append(bnd, i, ka, kb)) != eslOK) goto ERROR;
  }
  
  /* Finalize the band structure */
  bnd->L = L;
  bnd->M = M;
  
  *ret_bnd = bnd;
  return eslOK;
  
 ERROR:
  if (bnd != NULL) p7_gbands_Destroy(bnd);
  *ret_bnd = NULL;
  return status;
}

/* Function: my_p7_GForwardBanded()
 * Date:     EPN*, Sun Mar 16 2026
 * 
 * Purpose:  Exact copy of p7_GForwardBanded() from HMMER.
 *           Copied from hmmer/src/generic_fwdback_banded.c:p7_GForwardBanded()
 *           GLOCAL mode only (entry at M_1, exit at M_M).
 *           
 * Args:     dsq    - digital sequence, 1..L
 *           L      - length of sequence
 *           gm     - profile
 *           gxb    - banded DP matrix (contains P7_GBANDS structure)
 *           opt_sc - optRETURN: Forward score in nats
 *           
 * Returns:  <eslOK> on success
 */
int
my_p7_GForwardBanded(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMXB *gxb, float *opt_sc)
{
  int         *bnd_ip = gxb->bnd->imem;          /* ptr to current ia, ib segment band in gxb->bnd */
  int         *bnd_kp = gxb->bnd->kmem;		 /* ptr to current ka, kb row band in gxb->bnd     */
  float       *dpc    = gxb->dp;	         /* ptr to current DP matrix cell */
  float       *xpc    = gxb->xmx;		 /* ptr to current special cell   */
  float const *tsc    = gm->tsc;		 /* sets up TSC() macro, access to profile's transitions */
  float const *rsc;				 /* will be set up for MSC(), ISC() macros for residue scores */
  float       *dpp;	                  	 /* ptr to previous DP matrix cell */
  float       *last_dpc;			 /* used to reinitialize dpp after each row        */
  int          ia, ib;				 /* current segment band is rows ia..ib            */
  int          last_ib;				 /* intersegment interval is last_ib+1..ia-1       */
  int          kac, kbc;			 /* current row band is kac..kbc                   */
  int          kap, kbp;			 /* previous row band is kap..kbp                  */
  int          kbc2;				 /* if kbc==M, kbc2=M-1, main loop goes kac..kbc2 and M is unrolled */
  float        xE, xN, xJ, xB, xC;               /* tmp scores on special states. only stored when in row bands */
  float        mvp, ivp, dvp;			 /* M,I,D cell values from previous row i-1     */
  float        dc;				 /* precalculated D(i,k+1) value on current row */
  float        sc;				 /* temporary score calculation M(i,k)          */
  int          g, i, k;				 /* indices running over segments, residues (rows) x_i, model positions (cols) k  */
  float        esc  = p7_profile_IsLocal(gm) ? 0 : -eslINFINITY;
  
  xN      = 0.0f;
  xJ      = -eslINFINITY;
  xC      = -eslINFINITY;
  last_ib = 0;


  for (g = 0; g < gxb->bnd->nseg; g++)
    {
      ia = *bnd_ip++;
      ib = *bnd_ip++;

      /* kap,kbp initialization for i=ia:
       *  left overhang dpp advance must always eval to 0, 
       *  {m,i,d}vp initialization must always eval to -eslINFINITY.
       */
      kap = kbp = gm->M+1;   
      dpp = dpc;		/* re-initialize dpp */
      
      /* re-initialization: specials for previous row ia-1 just outside banded segment.
       * Guard n*xsc[LOOP] against n=0 with xsc[LOOP]=-inf: IEEE 754 gives 0*(-inf)=NaN.
       * This occurs with truncated profiles (UNILOCAL/UNIGLOCAL) where LOOP = -inf.
       */
      xE  = -eslINFINITY;
      { int gap = ia - last_ib - 1;
        if (gap > 0) { xN = xN + gap * gm->xsc[p7P_N][p7P_LOOP];
	               xJ = xJ + gap * gm->xsc[p7P_J][p7P_LOOP];
	               xC = xC + gap * gm->xsc[p7P_C][p7P_LOOP]; } }
      xB  = p7_FLogsum( xN + gm->xsc[p7P_N][p7P_MOVE], xJ + gm->xsc[p7P_J][p7P_MOVE]);

      for (i = ia; i <= ib; i++)
	{
	  rsc      = gm->rsc[dsq[i]];   /* sets up MSC(k), ISC(k) residue scores for this row i */
	  dc       = -eslINFINITY;
	  xE       = -eslINFINITY;
	  last_dpc = dpc;

	  kac      = *bnd_kp++;         /* current row's band is cells k=kac..kbc  */
	  kbc      = *bnd_kp++;
	  kbc2     = (kbc == gm->M ? kbc-1 : kbc); /* a "do_M" flag works too, but this way we avoid an if statement */

	  /* DEBUG: check for dpc overflow before writing this row */

	  /* dpp must advance by any left overhang of previous row; but no more than the entire row */
	  dpp += (kac-1 > kap ? ESL_MIN(kac-kap-1, kbp-kap+1) * p7G_NSCELLS : 0);

	  if (kac > kap && kac-1 <= kbp) { mvp = *dpp++;       ivp = *dpp++;       dvp = *dpp++;       }
	  else                           { mvp = -eslINFINITY; ivp = -eslINFINITY; dvp = -eslINFINITY; }

	  for (k = kac; k <= kbc2; k++)
	    {
	      *dpc++ = sc = MSC(k) + p7_FLogsum( p7_FLogsum(mvp + TSC(p7P_MM, k-1), ivp + TSC(p7P_IM, k-1)),
						 p7_FLogsum(dvp + TSC(p7P_DM, k-1), xB  + TSC(p7P_BM, k-1)));


	      if (k >= kap && k <= kbp) {  mvp = *dpp++;       ivp = *dpp++;        dvp = *dpp++;       } 	      // an if seems unavoidable. alternatively, might unroll
	      else                      {  mvp = -eslINFINITY; ivp = -eslINFINITY;  dvp = -eslINFINITY; }	      // all possible (kap,kac)..(kbp,kbc) orderings, but this
                                                                                                                      // seems too complex

	      *dpc++ = ISC(k) + p7_FLogsum( mvp + TSC(p7P_MI, k), ivp + TSC(p7P_II, k));

	      xE     = p7_FLogsum( p7_FLogsum(sc + esc, dc + esc), xE);/* Mk->E accumulation      */

	      /* next D_k+1 */
	      *dpc++ = dc;
	      dc     = p7_FLogsum( sc + TSC(p7P_MD, k), dc + TSC(p7P_DD, k));	     
	    }

	  if (kbc2 < kbc) /* i.e., if kbc==M and we need to do the final M column: */
	    {
	      *dpc++ = sc = MSC(k) + p7_FLogsum( p7_FLogsum(mvp + TSC(p7P_MM, k-1), ivp + TSC(p7P_IM, k-1)),
						 p7_FLogsum(dvp + TSC(p7P_DM, k-1), xB  + TSC(p7P_BM, k-1)));
	      *dpc++ = -eslINFINITY; 
	      *dpc++ = dc;           
	      xE     = p7_FLogsum( p7_FLogsum(sc, dc), xE);
	    }

	  *xpc++ = xE;
	  *xpc++ = xN = xN + gm->xsc[p7P_N][p7P_LOOP];
	  *xpc++ = xJ = p7_FLogsum( xJ + gm->xsc[p7P_J][p7P_LOOP],  xE + gm->xsc[p7P_E][p7P_LOOP]);
	  *xpc++ = xB = p7_FLogsum( xJ + gm->xsc[p7P_J][p7P_MOVE],  xN + gm->xsc[p7P_N][p7P_MOVE]);
	  *xpc++ = xC = p7_FLogsum( xE + gm->xsc[p7P_E][p7P_MOVE],  xC + gm->xsc[p7P_C][p7P_LOOP]);

	  dpp = last_dpc;	/* this skips any right overhang on the previous row, so dpp advances (if necessary) to start of curr row */
	  kap = kac;
	  kbp = kbc;
	}
      last_ib = ib;
    }

  /* last_ib+1..L is outside any band segment, so it can only run through xC.
   * Guard against 0*(-inf)=NaN when tail=0 and xsc[C][LOOP]=-inf (truncated profiles).
   */
  if (opt_sc != NULL) { int tail = L - last_ib;
    *opt_sc = (tail > 0 ? xC + tail * gm->xsc[p7P_C][p7P_LOOP] : xC) + gm->xsc[p7P_C][p7P_MOVE]; }
  return eslOK;
}

/* Function: p7_GBackwardBanded()
 * Date:     EPN*, Sun Mar 16 2026
 * 
 * Purpose:  Banded Backward algorithm for P7 profile HMMs.
 *           Adapted from Sean Eddy's p7_GForwardBanded() in 
 *           hmmer/src/generic_fwdback_banded.c (2011) and p7_GBackward()
 *           in hmmer/src/generic_fwdback.c.
 *           
 *           Calculate the Backward matrix for sequence <dsq> of length <L>,
 *           using profile <gm>, constrained by bands in <gxb>.
 *           
 *           Works backwards through the sequence (i = L to 1) and through
 *           model nodes within each row's band, computing the probability
 *           of generating the rest of the sequence from each cell.
 *           
 * Args:     dsq    - digital sequence, 1..L
 *           L      - length of sequence
 *           gm     - profile
 *           gxb    - banded DP matrix (contains P7_GBANDS structure)
 *           opt_sc - optRETURN: Backward score in nats
 *           
 * Returns:  <eslOK> on success; <gxb> contains the Backward matrix,
 *           <opt_sc> (if non-NULL) contains Backward score.
 *           
 * Note:     Like p7_GBackward(), this calculates the probability of
 *           getting OUT of cell i,k, exclusive of emitting residue x_i.
 *           
 *           Emissions for row i are scored using dsq[i+1], because
 *           Backward looks at what comes AFTER the current cell.
 */
int
p7_GBackwardBanded(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMXB *gxb, float *opt_sc)
{
  int         *bnd_ip;                           /* ptr to segment band boundaries in gxb->bnd */
  int         *bnd_kp;                           /* ptr to row band boundaries in gxb->bnd */
  float       *dpc;                              /* ptr to current DP matrix cell */
  float       *xpc;                              /* ptr to current special states */
  float const *tsc    = gm->tsc;                 /* transition scores */
  float const *rsc;                              /* residue scores for current row */
  float       *dpn;                              /* ptr to next row DP matrix cell */
  int          ia, ib;                           /* segment boundaries */
  int          next_ia;                          /* next segment's ia (for handling gaps) */
  int          kac, kbc;                         /* current row band: kac..kbc */
  int          kan, kbn;                         /* next row band: kan..kbn */
  int          kbc2;                             /* if kbc==M, kbc2=M-1 */
  float        xE, xN, xJ, xB, xC;               /* special state scores */
  float        mnext, inext, dnext;              /* M,I,D scores from next row */
  float        dc;                               /* current D(i,k) being calculated backwards */
  float        sc;                               /* temporary score for M(i,k) */
  int          g, i, k;                          /* segment, row, node indices */
  float        esc  = p7_profile_IsLocal(gm) ? 0 : -eslINFINITY;

  /* We'll traverse segments backwards, and reconstruct band pointers backwards.
   * Position pointers at the END of the arrays, we'll decrement them.
   */
  bnd_ip = gxb->bnd->imem + (gxb->bnd->nseg * 2);
  bnd_kp = gxb->bnd->kmem + (gxb->bnd->nrow * 2);
  dpc    = gxb->dp  + (gxb->bnd->ncell * p7G_NSCELLS);
  xpc    = gxb->xmx + (gxb->bnd->nrow  * p7G_NXCELLS);
  
  /* Initialize: handle anything after the last segment */
  xC      = gm->xsc[p7P_C][p7P_MOVE];
  xE      = xC + gm->xsc[p7P_E][p7P_MOVE];
  xJ      = xB = xN = -eslINFINITY;
  next_ia = L+1;

  /* Main recursion: work backwards through segments */
  for (g = gxb->bnd->nseg-1; g >= 0; g--)
    {
      ib = *(--bnd_ip);
      ia = *(--bnd_ip);

      /* Initialize special states for rows beyond this segment (next_ia-1...ib+1) */
      if (next_ia > ib+1) {
        xC = xC + (next_ia - ib - 1) * gm->xsc[p7P_C][p7P_LOOP];
        xE = p7_FLogsum(xE, xC + gm->xsc[p7P_E][p7P_MOVE]);
        /* J,B,N remain -infinity for these positions */
      }

      /* Initialize next row bands for i=ib+1 (or end of sequence) */
      kan = kbn = gm->M+1;
      dpn = dpc;

      /* Work backwards through rows in this segment */
      for (i = ib; i >= ia; i--)
        {
          rsc      = (i < L) ? gm->rsc[dsq[i+1]] : NULL;

          /* Get current row's band */
          kbc      = *(--bnd_kp);
          kac      = *(--bnd_kp);
          kbc2     = (kbc == gm->M ? kbc-1 : kbc);

          /* ROW L SPECIAL CASE: Initialize row L separately
           * Following unbanded p7_GBackward() pattern where row L is
           * initialized completely before main recursion.
           */
          if (i == L) {
            /* Special states for row L */
            *(--xpc) = xC;                       /* C */
            *(--xpc) = xB = -eslINFINITY;        /* B */
            *(--xpc) = xJ = -eslINFINITY;        /* J */
            *(--xpc) = xN = -eslINFINITY;        /* N */
            *(--xpc) = xE;                       /* E */


            /* Initialize M_M, I_M, D_M */
            if (kbc == gm->M) {
              *(--dpc) = xE;              /* D_M <- E */
              *(--dpc) = -eslINFINITY;    /* I_M (doesn't exist) */
              *(--dpc) = xE;              /* M_M <- E */
            }

            /* Backwards sweep through k from M-1 (or kbc2) down to kac
             * Row L formula: MMX(L,k) = p7_FLogsum(xE + esc, DMX(L,k+1) + TSC(p7P_MD,k))
             *                DMX(L,k) = p7_FLogsum(xE + esc, DMX(L,k+1) + TSC(p7P_DD,k))
             *                IMX(L,k) = -eslINFINITY
             */
            dc = (kbc == gm->M) ? xE : -eslINFINITY;  /* dc starts as D(L,kbc2+1) */
            for (k = kbc2; k >= kac; k--) {
              /* Calculate D(L,k) using D(L,k+1) which is in dc */
              float dk = p7_FLogsum(xE + esc, dc + TSC(p7P_DD, k));
              float mk = p7_FLogsum(xE + esc, dc + TSC(p7P_MD, k));
              
              *(--dpc) = dk;                                          /* Store D(L,k) */
              *(--dpc) = -eslINFINITY;                                /* Store I(L,k) = -inf */
              *(--dpc) = mk;                                          /* Store M(L,k) using D(L,k+1) */
              
              
              dc = dk;  /* D(L,k) becomes D(L,k+1) for next iteration */
            }
          
          /* Debug output for complete row L */

            dpn = dpc;  /* dpn now points to start of row L for row L-1 to read from */
            kan = kac;
            kbn = kbc;
            continue;  /* Skip normal row processing */
          }

          /* NORMAL ROW (i < L): Calculate special states that depend on next row i+1 */
          
          xB = -eslINFINITY;
          for (k = ESL_MAX(1,kan); k <= ESL_MIN(gm->M, kbn); k++) {
            if (k >= kan && k <= kbn) {
              int offset = (k - kan) * p7G_NSCELLS;
              float m_next = dpn[offset];  /* M(i+1,k) */
              float contrib = m_next + TSC(p7P_BM, k-1) + (rsc ? MSC(k) : 0);
              xB = p7_FLogsum(xB, contrib);
            }
          }
          
          
          xJ = p7_FLogsum(xJ + gm->xsc[p7P_J][p7P_LOOP], xB + gm->xsc[p7P_J][p7P_MOVE]);
          xC = xC + gm->xsc[p7P_C][p7P_LOOP];
          xE = p7_FLogsum(xJ + gm->xsc[p7P_E][p7P_LOOP], xC + gm->xsc[p7P_E][p7P_MOVE]);
          xN = p7_FLogsum(xN + gm->xsc[p7P_N][p7P_LOOP], xB + gm->xsc[p7P_N][p7P_MOVE]);

          /* Store special states */
          *(--xpc) = xC;
          *(--xpc) = xB;
          *(--xpc) = xJ;
          *(--xpc) = xN;
          *(--xpc) = xE;

          /* Handle M_M state: M_M gets E state */
          if (kbc == gm->M) {
            *(--dpc) = xE;              /* D_M */
            *(--dpc) = -eslINFINITY;    /* I_M (doesn't exist) */
            *(--dpc) = xE;              /* M_M */
          }

          /* Main recursion: work backwards through k = kbc2...kac */
          /* Initialize dc to D(i,M): if kbc == M, then D(i,M) = xE; otherwise -inf */
          dc = (kbc == gm->M) ? xE : -eslINFINITY;
          for (k = kbc2; k >= kac; k--)
            {
              /* Get scores from next row i+1.
               * mnext (M(i+1,k+1)) and dnext (D(i+1,k+1)) require k+1 in the next
               * row's band; inext (I(i+1,k)) requires k in the next row's band.
               * These two conditions are INDEPENDENT: on a vertical insert run the
               * band can hold k but not k+1, so inext must be tested separately or
               * the insert chain is silently broken (loses all mass through it).
               */
              if (k+1 >= kan && k+1 <= kbn) {
                int offset = (k+1 - kan) * p7G_NSCELLS;
                mnext = dpn[offset]   + (rsc ? MSC(k+1) : 0);
                dnext = dpn[offset+2];
              } else {
                mnext = dnext = -eslINFINITY;
              }
              if (k >= kan && k <= kbn) {
                int curr_offset = (k - kan) * p7G_NSCELLS;
                inext = dpn[curr_offset+1] + (rsc ? ISC(k) : 0);
              } else {
                inext = -eslINFINITY;
              }


              /* D(i,k) - compute BEFORE storing 
               * D(i,k) = logsum(M(i+1,k+1) + TSC(DM,k) + MSC(k+1),
               *                 D(i,k+1) + TSC(DD,k),
               *                 xE + esc)
               * Note: D(i,k+1) is in dc from previous iteration!
               */
              float dk = p7_FLogsum(p7_FLogsum(mnext + TSC(p7P_DM, k), dc + TSC(p7P_DD, k)),
                                    xE + esc);
              

              *(--dpc) = dk;   /* Store D(i,k) */

              /* I(i,k) */
              *(--dpc) = p7_FLogsum(mnext + TSC(p7P_IM, k), inext + TSC(p7P_II, k));

              /* M(i,k) - uses dc which is D(i,k+1) from previous iteration */
              float mk = p7_FLogsum(p7_FLogsum(mnext + TSC(p7P_MM, k), inext + TSC(p7P_MI, k)),
                                        p7_FLogsum(xE + esc, dc + TSC(p7P_MD, k)));
              *(--dpc) = sc = mk;

              /* Update dc for next iteration */
              dc = dk;
              
            }
          
          /* Debug output for complete rows L and L-1 */

          dpn = dpc;  /* dpn now points to start of row i we just wrote, for next iteration's i-1 to read from */
          kan = kac;
          kbn = kbc;
        }
      next_ia = ia;
    }

  /* Handle row i=0: At i=0, only N,B are reachable 
   * Following unbanded p7_GBackward() pattern */
  if (opt_sc != NULL) {
    float xB_0, xN_0;
    float const *rsc = gm->rsc[dsq[1]];
    
    /* Account for N-loops from row next_ia-1 down to row 1 */
    if (next_ia > 1) {
      xN = xN + (next_ia - 1) * gm->xsc[p7P_N][p7P_LOOP];
    }
    
    /* Calculate B(0): sum over all M states in row 1 within bands
     * B(0) = sum_k M(1,k) + TSC(BM,k-1) + MSC(k)
     */
    xB_0 = -eslINFINITY;
    /* Find the bands for row i=1 (the first row after row 0) */
    int k1min, k1max;
    if (next_ia <= 1) {
      /* Row 1 is in a segment, find its bands */
      int *kp = gxb->bnd->kmem;  /* Start of band array */
      k1min = kp[0];
      k1max = kp[1];
      
      /* Sum over M(1,k) for k in band */
      for (int k = k1min; k <= k1max; k++) {
        if (k >= 1 && k <= gm->M) {
          /* Access M(1,k) from the banded matrix */
          int k_offset = (k - k1min) * p7G_NSCELLS;
          float m1k = gxb->dp[k_offset + p7G_M];  /* M(1,k) */
          xB_0 = p7_FLogsum(xB_0, m1k + TSC(p7P_BM, k-1) + MSC(k));
        }
      }
    } else {
      /* Row 1 is before any segment, so all M(1,k) = -inf, thus B(0) = -inf */
      xB_0 = -eslINFINITY;
    }
    
    /* Calculate N(0) = logsum(N(1) + N_LOOP, B(0) + N_MOVE) */
    xN_0 = p7_FLogsum(xN + gm->xsc[p7P_N][p7P_LOOP], xB_0 + gm->xsc[p7P_N][p7P_MOVE]);

    *opt_sc = xN_0;
  }

  return eslOK;
}


/* Function: p7_GDecodingBanded()
 * Date:     EPN*, Wed Mar 18 2026
 *
 * Purpose:  Banded posterior decoding. Given Forward and Backward banded
 *           matrices <fwd> and <bck> (both using the same P7_GBANDS),
 *           calculate posterior probabilities for each (i,k) cell.
 *
 *           <bck> may be overwritten with posteriors (pass bck == pp).
 *           <fwd> is read-only (needs fwd xmx at row i-1 for N/J/C).
 *
 * Args:     gm         - profile
 *           fwd        - banded Forward matrix
 *           bck        - banded Backward matrix
 *           pp         - RESULT: banded posterior matrix (may == bck)
 *           overall_sc - Forward score in nats (fwdsc)
 *
 * Returns:  eslOK on success.
 */
int
p7_GDecodingBanded(const P7_PROFILE *gm, const P7_GMXB *fwd, P7_GMXB *bck,
		   P7_GMXB *pp, float overall_sc)
{
  int         *bnd_ip = fwd->bnd->imem;
  int         *bnd_kp = fwd->bnd->kmem;
  float const *fwd_dp = fwd->dp;
  float const *bck_dp = bck->dp;
  float       *pp_dp  = pp->dp;
  float const *fwd_xp = fwd->xmx;
  float const *bck_xp = bck->xmx;
  float       *pp_xp  = pp->xmx;
  int          ia, ib;
  int          last_ib;
  int          kac, kbc;
  int          M = gm->M;
  float        fwd_xN_prev, fwd_xJ_prev, fwd_xC_prev;
  float        denom;
  int          g, i, k;
  int          nk;

  fwd_xN_prev = 0.0f;
  fwd_xJ_prev = -eslINFINITY;
  fwd_xC_prev = -eslINFINITY;
  last_ib     = 0;

  for (g = 0; g < fwd->bnd->nseg; g++)
    {
      ia = *bnd_ip++;
      ib = *bnd_ip++;

      /* Handle gap before this segment: advance prev specials through gap.
       * Guard against 0*(-inf)=NaN.
       */
      { int gap = ia - last_ib - 1;
        if (gap > 0) {
	  fwd_xN_prev = fwd_xN_prev + gap * gm->xsc[p7P_N][p7P_LOOP];
	  fwd_xJ_prev = fwd_xJ_prev + gap * gm->xsc[p7P_J][p7P_LOOP];
	  fwd_xC_prev = fwd_xC_prev + gap * gm->xsc[p7P_C][p7P_LOOP];
	}
      }

      for (i = ia; i <= ib; i++)
	{
	  kac = *bnd_kp++;
	  kbc = *bnd_kp++;
	  denom = 0.0f;

	  for (k = kac; k <= kbc; k++)
	    {
	      /* M posterior */
	      *pp_dp = expf(*fwd_dp + *bck_dp - overall_sc);
	      if (! isfinite(*pp_dp)) *pp_dp = 0.0f;  /* guard NaN from -inf + inf (cell unreachable in one direction) */
	      denom += *pp_dp;
	      pp_dp++; fwd_dp++; bck_dp++;

	      /* I posterior */
	      if (k < M) {
		*pp_dp = expf(*fwd_dp + *bck_dp - overall_sc);
		if (! isfinite(*pp_dp)) *pp_dp = 0.0f;  /* guard NaN from -inf + inf (cell unreachable in one direction) */
		denom += *pp_dp;
	      } else {
		*pp_dp = 0.0f;
	      }
	      pp_dp++; fwd_dp++; bck_dp++;

	      /* D posterior = 0 (D doesn't emit) */
	      *pp_dp = 0.0f;
	      pp_dp++; fwd_dp++; bck_dp++;
	    }

	  /* Special states: E=0, B=0; N/J/C use fwd(i-1) and bck(i) */
	  pp_xp[p7G_E] = 0.0f;

	  pp_xp[p7G_N] = expf(fwd_xN_prev + bck_xp[p7G_N] + gm->xsc[p7P_N][p7P_LOOP] - overall_sc);
	  if (! isfinite(pp_xp[p7G_N])) pp_xp[p7G_N] = 0.0f;  /* guard NaN from -inf + -inf */
	  denom += pp_xp[p7G_N];

	  pp_xp[p7G_J] = expf(fwd_xJ_prev + bck_xp[p7G_J] + gm->xsc[p7P_J][p7P_LOOP] - overall_sc);
	  if (! isfinite(pp_xp[p7G_J])) pp_xp[p7G_J] = 0.0f;
	  denom += pp_xp[p7G_J];

	  pp_xp[p7G_B] = 0.0f;

	  pp_xp[p7G_C] = expf(fwd_xC_prev + bck_xp[p7G_C] + gm->xsc[p7P_C][p7P_LOOP] - overall_sc);
	  if (! isfinite(pp_xp[p7G_C])) pp_xp[p7G_C] = 0.0f;
	  denom += pp_xp[p7G_C];

	  /* Row normalization */
	  if (denom > 0.0f) {
	    denom = 1.0f / denom;
	    nk = (kbc - kac + 1);
	    { float *p = pp_dp - nk * p7G_NSCELLS;
	      for (k = kac; k <= kbc; k++) {
		*p++ *= denom; /* M */
		*p++ *= denom; /* I */
		p++;           /* D stays 0 */
	      }
	    }
	    pp_xp[p7G_N] *= denom;
	    pp_xp[p7G_J] *= denom;
	    pp_xp[p7G_C] *= denom;
	  }

	  /* Save current row's fwd specials as prev for next row */
	  fwd_xN_prev = fwd_xp[p7G_N];
	  fwd_xJ_prev = fwd_xp[p7G_J];
	  fwd_xC_prev = fwd_xp[p7G_C];

	  fwd_xp += p7G_NXCELLS;
	  bck_xp += p7G_NXCELLS;
	  pp_xp  += p7G_NXCELLS;
	}
      last_ib = ib;
    }
  return eslOK;
}


/* Function: p7_GOptimalAccuracyBanded()
 * Date:     EPN*, Wed Mar 18 2026
 *
 * Purpose:  Banded optimal accuracy DP fill. Same structure as banded
 *           Forward, but uses ESL_MAX instead of p7_FLogsum, and adds
 *           posterior probability as an accuracy reward for emitting states.
 *
 * Args:     gm    - profile
 *           pp    - banded posterior matrix (from p7_GDecodingBanded)
 *           gx    - RESULT: banded OA DP matrix (may share banding with pp)
 *           ret_e - RETURN: OA score (expected # correct residues)
 *
 * Returns:  eslOK on success.
 */
int
p7_GOptimalAccuracyBanded(const P7_PROFILE *gm, const P7_GMXB *pp,
			  P7_GMXB *gx, float *ret_e)
{
  /* TSCDELTA: 1.0 if transition is possible, FLT_MIN if not */
#define TSCDELTA(s,k) ( (tsc[(k) * p7P_NTRANS + (s)] == -eslINFINITY) ? FLT_MIN : 1.0)

  int         *bnd_ip  = gx->bnd->imem;
  int         *bnd_kp  = gx->bnd->kmem;
  float       *dpc     = gx->dp;
  float       *xpc     = gx->xmx;
  float const *tsc     = gm->tsc;
  float const *pp_dp;
  float const *pp_xp   = pp->xmx;
  float       *dpp;
  float       *last_dpc;
  int          ia, ib;
  int          last_ib;
  int          kac, kbc;
  int          kap, kbp;
  int          kbc2;
  float        xE, xN, xJ, xB, xC;
  float        mvp, ivp, dvp;
  float        dc;
  float        sc;
  float        pp_m;
  int          g, i, k;
  int          M = gm->M;
  float        esc = p7_profile_IsLocal(gm) ? 1.0f : 0.0f;
  float        t1, t2;

  xN      = 0.0f;
  xJ      = -eslINFINITY;
  xC      = -eslINFINITY;
  last_ib = 0;
  pp_dp   = pp->dp;

  for (g = 0; g < gx->bnd->nseg; g++)
    {
      ia = *bnd_ip++;
      ib = *bnd_ip++;

      kap = kbp = M + 1;
      dpp = dpc;

      /* Re-init specials for gap before this segment.
       * In OA, N(i) = t1*(N(i-1) + pp_N(i)). For rows outside bands,
       * pp_N is not stored, but for UNILOCAL N_LOOP=-inf so t1=FLT_MIN
       * and the N contribution goes to ~0 anyway. Just mark as -inf for gap>0.
       */
      xE = -eslINFINITY;
      { int gap = ia - last_ib - 1;
	if (gap > 0) {
	  xN = -eslINFINITY;
	  xJ = -eslINFINITY;
	  xC = -eslINFINITY;
	}
      }
      xB = ESL_MAX( ((gm->xsc[p7P_N][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0f) * xN,
		     ((gm->xsc[p7P_J][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0f) * xJ);

      for (i = ia; i <= ib; i++)
	{
	  dc       = -eslINFINITY;
	  xE       = -eslINFINITY;
	  last_dpc = dpc;

	  kac = *bnd_kp++;
	  kbc = *bnd_kp++;
	  kbc2 = (kbc == M ? kbc - 1 : kbc);

	  /* dpp advance for left overhang */
	  dpp += (kac - 1 > kap ? ESL_MIN(kac - kap - 1, kbp - kap + 1) * p7G_NSCELLS : 0);

	  if (kac > kap && kac - 1 <= kbp) { mvp = *dpp++; ivp = *dpp++; dvp = *dpp++; }
	  else                              { mvp = -eslINFINITY; ivp = -eslINFINITY; dvp = -eslINFINITY; }

	  for (k = kac; k <= kbc2; k++)
	    {
	      pp_m = *pp_dp;  /* pp_M(i,k) */

	      /* M(i,k) = max(transitions) + pp_M(i,k) */
	      *dpc++ = sc = ESL_MAX(ESL_MAX(TSCDELTA(p7P_MM, k-1) * (mvp + pp_m),
					    TSCDELTA(p7P_IM, k-1) * (ivp + pp_m)),
				    ESL_MAX(TSCDELTA(p7P_DM, k-1) * (dvp + pp_m),
					    TSCDELTA(p7P_BM, k-1) * (xB  + pp_m)));

	      if (k >= kap && k <= kbp) { mvp = *dpp++; ivp = *dpp++; dvp = *dpp++; }
	      else                      { mvp = -eslINFINITY; ivp = -eslINFINITY; dvp = -eslINFINITY; }

	      /* I(i,k) = max(transitions) + pp_I(i,k) */
	      *dpc++ = ESL_MAX(TSCDELTA(p7P_MI, k) * (mvp + pp_dp[1]),
			       TSCDELTA(p7P_II, k) * (ivp + pp_dp[1]));

	      /* E update from M(i,k) */
	      xE = ESL_MAX(xE, esc * sc);

	      /* D(i,k) — no pp reward */
	      *dpc++ = dc;
	      dc = ESL_MAX(TSCDELTA(p7P_MD, k) * sc,
			   TSCDELTA(p7P_DD, k) * dc);

	      pp_dp += p7G_NSCELLS;
	    }

	  if (kbc2 < kbc) /* kbc == M: unrolled last node */
	    {
	      pp_m = *pp_dp;

	      *dpc++ = sc = ESL_MAX(ESL_MAX(TSCDELTA(p7P_MM, k-1) * (mvp + pp_m),
					    TSCDELTA(p7P_IM, k-1) * (ivp + pp_m)),
				    ESL_MAX(TSCDELTA(p7P_DM, k-1) * (dvp + pp_m),
					    TSCDELTA(p7P_BM, k-1) * (xB  + pp_m)));
	      *dpc++ = -eslINFINITY; /* no I_M */
	      *dpc++ = dc;           /* D_M */

	      /* E update: M_M and D_M always exit to E (no esc penalty at k=M) */
	      xE = ESL_MAX(xE, ESL_MAX(sc, dc));

	      pp_dp += p7G_NSCELLS;
	    }

	  /* Special states */
	  *xpc++ = xE;

	  t1 = ((gm->xsc[p7P_N][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0f);
	  *xpc++ = xN = t1 * (xN + pp_xp[p7G_N]);

	  t1 = ((gm->xsc[p7P_J][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0f);
	  t2 = ((gm->xsc[p7P_E][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0f);
	  *xpc++ = xJ = ESL_MAX(t1 * (xJ + pp_xp[p7G_J]), t2 * xE);

	  t1 = ((gm->xsc[p7P_N][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0f);
	  t2 = ((gm->xsc[p7P_J][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0f);
	  *xpc++ = xB = ESL_MAX(t1 * xN, t2 * xJ);

	  t1 = ((gm->xsc[p7P_C][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0f);
	  t2 = ((gm->xsc[p7P_E][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0f);
	  *xpc++ = xC = ESL_MAX(t1 * (xC + pp_xp[p7G_C]), t2 * xE);

	  pp_xp += p7G_NXCELLS;
	  dpp = last_dpc;
	  kap = kac;
	  kbp = kbc;
	}
      last_ib = ib;
    }

  if (ret_e != NULL) *ret_e = xC;
  return eslOK;

#undef TSCDELTA
}


/* Function: p7_GOATraceBanded()
 * Date:     EPN*, Wed Mar 18 2026
 *
 * Purpose:  Banded OA traceback. Traces back through banded OA matrix
 *           to find the optimal accuracy alignment path.
 *
 *           Builds a lookup table for random access into the banded
 *           matrices, then runs the same traceback logic as p7_GOATrace().
 *
 * Args:     gm  - profile
 *           pp  - banded posterior matrix
 *           gx  - banded OA DP matrix
 *           tr  - RESULT: OA trace (caller provides, possibly via Reuse)
 *
 * Returns:  eslOK on success.
 *           eslEMEM on allocation failure.
 */
int
p7_GOATraceBanded(const P7_PROFILE *gm, const P7_GMXB *pp, const P7_GMXB *gx,
		  P7_TRACE *tr)
{
#define TSCDELTA(s,k) ( (tsc[(k) * p7P_NTRANS + (s)] == -eslINFINITY) ? FLT_MIN : 1.0)

  float const *tsc = gm->tsc;
  P7_GBANDS   *bnd = gx->bnd;
  int          L    = bnd->L;
  int          M    = gm->M;
  int          nrow = bnd->nrow;

  /* Lookup tables for random access */
  int64_t *dp_off  = NULL;   /* dp offset for each banded row */
  int     *ka_arr  = NULL;   /* band ka for each banded row */
  int     *kb_arr  = NULL;   /* band kb for each banded row */
  int     *row_idx = NULL;   /* seq pos i -> banded row index r (-1 if outside) */

  int      i, k, r, g;
  int      ia, ib;
  int     *ip, *kp;
  int64_t  cum;
  float    postprob;
  int      sprv, scur;
  int      status;

  /* Allocate lookup tables */
  ESL_ALLOC(dp_off,  sizeof(int64_t) * nrow);
  ESL_ALLOC(ka_arr,  sizeof(int)     * nrow);
  ESL_ALLOC(kb_arr,  sizeof(int)     * nrow);
  ESL_ALLOC(row_idx, sizeof(int)     * (L + 1));

  /* Fill lookup tables */
  kp  = bnd->kmem;
  cum = 0;
  for (r = 0; r < nrow; r++) {
    ka_arr[r]  = *kp++;
    kb_arr[r]  = *kp++;
    dp_off[r]  = cum;
    cum       += (int64_t)(kb_arr[r] - ka_arr[r] + 1) * p7G_NSCELLS;
  }

  for (i = 0; i <= L; i++) row_idx[i] = -1;
  ip = bnd->imem;
  r  = 0;
  for (g = 0; g < bnd->nseg; g++) {
    ia = *ip++;
    ib = *ip++;
    for (i = ia; i <= ib; i++)
      row_idx[i] = r++;
  }

  /* --- Helper macros for banded cell access --- */
#define GXB_M(gxb, i, k) \
  ( ((i) >= 0 && (i) <= L && row_idx[(i)] >= 0 && (k) >= ka_arr[row_idx[(i)]] && (k) <= kb_arr[row_idx[(i)]]) \
    ? (gxb)->dp[dp_off[row_idx[(i)]] + ((k) - ka_arr[row_idx[(i)]]) * p7G_NSCELLS + p7G_M] \
    : -eslINFINITY )

#define GXB_I(gxb, i, k) \
  ( ((i) >= 0 && (i) <= L && row_idx[(i)] >= 0 && (k) >= ka_arr[row_idx[(i)]] && (k) <= kb_arr[row_idx[(i)]]) \
    ? (gxb)->dp[dp_off[row_idx[(i)]] + ((k) - ka_arr[row_idx[(i)]]) * p7G_NSCELLS + p7G_I] \
    : -eslINFINITY )

#define GXB_D(gxb, i, k) \
  ( ((i) >= 0 && (i) <= L && row_idx[(i)] >= 0 && (k) >= ka_arr[row_idx[(i)]] && (k) <= kb_arr[row_idx[(i)]]) \
    ? (gxb)->dp[dp_off[row_idx[(i)]] + ((k) - ka_arr[row_idx[(i)]]) * p7G_NSCELLS + p7G_D] \
    : -eslINFINITY )

  /* For xmx at i=0: OA fill initializes N(0)=0, B(0)=0, E/C/J(0)=-inf.
   * Row 0 is not in the banded matrix, so return these initial values directly.
   */
#define GXB_XMX(gxb, i, s) \
  ( ((i) == 0) \
    ? (((s) == p7G_N || (s) == p7G_B) ? 0.0f : -eslINFINITY) \
    : (((i) > 0 && (i) <= L && row_idx[(i)] >= 0) \
       ? (gxb)->xmx[row_idx[(i)] * p7G_NXCELLS + (s)] \
       : -eslINFINITY) )

  /* --- Traceback --- */
  i = L;
  k = 0;

  if ((status = p7_trace_AppendWithPP(tr, p7T_T, k, i, 0.0)) != eslOK) goto ERROR;
  if ((status = p7_trace_AppendWithPP(tr, p7T_C, k, i, 0.0)) != eslOK) goto ERROR;

  sprv = p7T_C;
  while (sprv != p7T_S)
    {
      switch (sprv) {

      case p7T_M: /* select_m: which state was predecessor of M(i,k)? */
	{
	  float path[4];
	  int   state[4] = { p7T_M, p7T_I, p7T_D, p7T_B };
	  path[0] = TSCDELTA(p7P_MM, k-1) * GXB_M(gx, i-1, k-1);
	  path[1] = TSCDELTA(p7P_IM, k-1) * GXB_I(gx, i-1, k-1);
	  path[2] = TSCDELTA(p7P_DM, k-1) * GXB_D(gx, i-1, k-1);
	  path[3] = TSCDELTA(p7P_BM, k-1) * GXB_XMX(gx, i-1, p7G_B);
	  scur = state[esl_vec_FArgMax(path, 4)];
	  k--; i--;
	}
	break;

      case p7T_D: /* select_d */
	{
	  float path[2];
	  path[0] = TSCDELTA(p7P_MD, k-1) * GXB_M(gx, i, k-1);
	  path[1] = TSCDELTA(p7P_DD, k-1) * GXB_D(gx, i, k-1);
	  scur = (path[0] >= path[1]) ? p7T_M : p7T_D;
	  k--;
	}
	break;

      case p7T_I: /* select_i */
	{
	  float path[2];
	  path[0] = TSCDELTA(p7P_MI, k) * GXB_M(gx, i-1, k);
	  path[1] = TSCDELTA(p7P_II, k) * GXB_I(gx, i-1, k);
	  scur = (path[0] >= path[1]) ? p7T_M : p7T_I;
	  i--;
	}
	break;

      case p7T_N: /* select_n */
	scur = (i == 0) ? p7T_S : p7T_N;
	break;

      case p7T_C: /* select_c */
	{
	  float t1c = ((gm->xsc[p7P_C][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0f);
	  float t2c = ((gm->xsc[p7P_E][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0f);
	  float path[2];
	  path[0] = t1c * (GXB_XMX(gx, i-1, p7G_C) + GXB_XMX(pp, i, p7G_C));
	  path[1] = t2c *  GXB_XMX(gx, i, p7G_E);
	  scur = (path[0] > path[1]) ? p7T_C : p7T_E;
	}
	break;

      case p7T_J: /* select_j */
	{
	  float t1j = ((gm->xsc[p7P_J][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0f);
	  float t2j = ((gm->xsc[p7P_E][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0f);
	  float path[2];
	  path[0] = t1j * (GXB_XMX(gx, i-1, p7G_J) + GXB_XMX(pp, i, p7G_J));
	  path[1] = t2j *  GXB_XMX(gx, i, p7G_E);
	  scur = (path[0] > path[1]) ? p7T_J : p7T_E;
	}
	break;

      case p7T_E: /* select_e: which k did E come from? */
	{
	  float max  = -eslINFINITY;
	  int   smax = -1;
	  int   kmax = -1;
	  int   ri   = row_idx[i];

	  if (! p7_profile_IsLocal(gm)) {
	    k = M;
	    scur = (GXB_M(gx, i, M) >= GXB_D(gx, i, M)) ? p7T_M : p7T_D;
	  } else {
	    if (ri >= 0) {
	      int ka = ka_arr[ri], kb = kb_arr[ri];
	      for (k = ka; k <= kb; k++) {
		float mv = GXB_M(gx, i, k);
		float dv = GXB_D(gx, i, k);
		if (mv >= max) { max = mv; smax = p7T_M; kmax = k; }
		if (dv >  max) { max = dv; smax = p7T_D; kmax = k; }
	      }
	    }
	    k    = kmax;
	    scur = smax;
	  }
	}
	break;

      case p7T_B: /* select_b */
	{
	  float t1b = ((gm->xsc[p7P_N][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0f);
	  float t2b = ((gm->xsc[p7P_J][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0f);
	  float path[2];
	  path[0] = t1b * GXB_XMX(gx, i, p7G_N);
	  path[1] = t2b * GXB_XMX(gx, i, p7G_J);
	  /* Brief 26_0430-131: under a unihit profile (E->J == -inf, set by
	   * p7_ProfileConfig for UNILOCAL/UNIGLOCAL), the J state is
	   * semantically unreachable, so B must enter from N. The OA fill
	   * floors forbidden transitions at FLT_MIN rather than hard -inf, so
	   * xJ can be 0 (not -inf); combined with the strict-'>' tie-break
	   * below, an xN==xJ==0 tie (which IBV's full-width boundary widening
	   * at rows 1/L-1/L induces) would otherwise default to a spurious
	   * p7T_J that map_new_msa() rejects ("J state unsupported"). Forcing N
	   * here is correct, not heuristic: with E->J blocked no residue can
	   * legitimately reach J. NOTE: do not "resync to stock" generic_optacc
	   * select_b and drop this guard -- stock's unbanded path never produces
	   * the boundary tie, but the IBV-banded path does. */
	  if (gm->xsc[p7P_E][p7P_LOOP] == -eslINFINITY) scur = p7T_N;
	  else scur = (path[0] > path[1]) ? p7T_N : p7T_J;
	}
	break;

      default:
	ESL_EXCEPTION(eslEINVAL, "bogus state in traceback");
      }

      if (scur == -1) ESL_EXCEPTION(eslEINVAL, "OA banded traceback choice failed");

      /* get_postprob */
      switch (scur) {
      case p7T_M: postprob = GXB_M(pp, i, k);                                      break;
      case p7T_I: postprob = GXB_I(pp, i, k);                                      break;
      case p7T_N: postprob = (scur == sprv) ? GXB_XMX(pp, i, p7G_N) : 0.0f;       break;
      case p7T_C: postprob = (scur == sprv) ? GXB_XMX(pp, i, p7G_C) : 0.0f;       break;
      case p7T_J: postprob = (scur == sprv) ? GXB_XMX(pp, i, p7G_J) : 0.0f;       break;
      default:    postprob = 0.0f;                                                  break;
      }

      if ((status = p7_trace_AppendWithPP(tr, scur, k, i, postprob)) != eslOK) goto ERROR;

      /* For NCJ self-loops, defer i decrement */
      if ((scur == p7T_N || scur == p7T_J || scur == p7T_C) && scur == sprv) i--;
      sprv = scur;
    }

  tr->M = M;
  tr->L = L;
  status = p7_trace_Reverse(tr);

  free(dp_off);
  free(ka_arr);
  free(kb_arr);
  free(row_idx);
  return status;

 ERROR:
  if (dp_off)  free(dp_off);
  if (ka_arr)  free(ka_arr);
  if (kb_arr)  free(kb_arr);
  if (row_idx) free(row_idx);
  return status;

#undef TSCDELTA
#undef GXB_M
#undef GXB_I
#undef GXB_D
#undef GXB_XMX
}


/*****************************************************************
 * Brief 26_0526-016: sqrt(M)-row checkpointed banded F / B / Decoding / OA
 *
 * The banded Forward/Backward matrices (P7_GMXB bxf/bxb) are the
 * dominant per-thread memory term at genome scale (brief 26_0526-014). This
 * block reduces that term by sqrt(nrow)-checkpointing the row axis:
 * instead of storing all <nrow> banded rows of the Forward (and OA)
 * matrix, store ~sqrt(nrow) "seed" rows and recompute the inter-seed
 * row blocks on demand. The Backward matrix is never materialized: it
 * streams (2 rows live) while its posterior is folded into the
 * resident posterior matrix <pp>.
 *
 * Construction mirrors 26_0610's linear-chain cm_CheckptAlignHB
 * (STEP A/B/OA/TB), adapted to the p7 banded DP whose checkpoint axis
 * is the sequence row i (not the CM state v). Byte-exact discipline:
 * the per-cell recurrences and p7_FLogsum / ESL_MAX call order are
 * copied VERBATIM from my_p7_GForwardBanded / p7_GBackwardBanded /
 * p7_GDecodingBanded / p7_GOptimalAccuracyBanded / p7_GOATraceBanded
 * above; only deck *storage* is windowed. p7_FLogsum is a symmetric
 * deterministic LUT, so identical inputs in identical order -> identical
 * bits. The existing non-checkpointed functions are NOT modified.
 *****************************************************************/

/* Brief 26_0526-017: the resident posterior <pp> stores only M and I per cell,
 * not the full p7G_NSCELLS=3 (M,I,D). pp_D is provably always 0 --
 * p7b_decode_row() set it to 0, the OA fill (p7b_oa_row) reads only
 * ppp[0]=M and ppp[1]=I, and the traceback reads pp only via PP_M/PP_I
 * (cells p7G_M=0, p7G_I=1) -- so no site ever reads pp_D. Storing 2
 * floats/cell instead of 3 shrinks the resident pp by ~1/3 (HSV genome
 * ~343 -> ~228 MB). This is a STORAGE-LAYOUT change only: the M and I
 * values written and read are bit-identical to before. The pp cell
 * stride is P7B_PP_NSCELLS; M stays at offset 0 (p7G_M), I at offset 1
 * (p7G_I). The Forward/Backward/OA scratch buffers keep p7G_NSCELLS. */
#define P7B_PP_NSCELLS 2   /* resident pp cell stride: M=0, I=1; D dropped */

/* p7b_pp_Create(): allocate a resident posterior matrix with the compact
 * 2-floats/cell dp (M,I).  Mirrors p7_gmxb_Create() (hmmer/, shared --
 * not edited) but sizes dp at P7B_PP_NSCELLS, not p7G_NSCELLS.  xmx is
 * unchanged (p7G_NXCELLS specials).  Freed with p7_gmxb_Destroy(). */
P7_GMXB *
p7b_pp_Create(P7_GBANDS *bnd)
{
  P7_GMXB *pp = NULL;
  int      status;

  ESL_ALLOC(pp, sizeof(P7_GMXB));
  pp->dp     = NULL;
  pp->xmx    = NULL;
  pp->bnd    = bnd;
  pp->dalloc = 0;
  pp->xalloc = 0;

  ESL_ALLOC(pp->dp,  sizeof(float) * bnd->ncell * P7B_PP_NSCELLS); /* M,I only (brief 26_0526-017) */
  ESL_ALLOC(pp->xmx, sizeof(float) * bnd->nrow  * p7G_NXCELLS);    /* ENJBC (0..4)        */
  pp->dalloc = bnd->ncell;
  pp->xalloc = bnd->nrow;
  return pp;

 ERROR:
  p7_gmxb_Destroy(pp);
  return NULL;
}

/* Function: p7_CheckptBandedOAMemNeeded()
 * Incept:   brief 26_0430-266
 *
 * Purpose:  Estimate the bytes actually allocated by the --hmm do_bandedoa
 *           (Mode 3, Viterbi-banded optimal accuracy) engine for banded
 *           region <bnd>, WITHOUT allocating anything. Three sub-cases,
 *           selected by <ckpt_mode> (brief 26_0628-081 added the third):
 *
 *             P7B_OAMEM_NOCKPT (0) (INFERNAL_HMM_CKPT_OFF set):
 *               the non-checkpointed engine -- 2 x p7_gmxb_Create(bnd),
 *               each a full banded dp (p7G_NSCELLS-wide) + xmx
 *               (p7G_NXCELLS-wide) allocation.  O(ncell).
 *
 *             P7B_OAMEM_CKPT (1) (the usual engine: this is what the
 *             caller picks whenever its estimate fits in --mxsize, and
 *             what INFERNAL_HMM_PPCKPT_OFF forces unconditionally):
 *               the singly-checkpointed engine -- p7b_pp_Create(bnd)'s
 *               resident 2-cell posterior + p7_GCheckptFBDecode_Banded()'s
 *               O(sqrt(nrow)) working set (p7b_forward_seeds()'s <=nblk
 *               seed[].dp row copies, p7b_backdecode()'s blkbuf (B rows)
 *               and bbuf0/bbuf1 (1 row each)).  Still O(ncell), because
 *               of the resident posterior.
 *
 *             P7B_OAMEM_CKPTPP (2) (the fallback the caller drops to when
 *             mode 1 would not fit in --mxsize; INFERNAL_HMM_PPCKPT_ON
 *             forces it):
 *               the DOUBLE-checkpointed engine, p7_GCheckptFBDecodeOA_Banded().
 *               No resident posterior at all: three ~nblk-entry seed
 *               arrays (Forward, Backward, OA) plus a handful of
 *               one-block (<=B+2 row) buffers -- the F, B and pp block
 *               buffers and the OA traceback window.  O(sqrt(nrow)*maxnc)
 *               in the DP term, plus O(nrow) for the P7B_GEO row tables
 *               and O(L) for the traceback's row_idx.  NOTE this term is
 *               independent of bnd->ncell: it is bounded by
 *               O(sqrt(nrow)*M) no matter how wide the band gets.
 *
 *           B and nblk mirror p7b_geo_Create()'s own derivation
 *           (B ~ round(sqrt(nrow)), nblk = ceil(nrow/B)); maxnc (widest
 *           row) is scanned from <bnd> directly since B/nblk/maxnc depend
 *           only on bnd, not on M/L/gm the way the full P7B_GEO struct does.
 *
 *           None of the three is O(M*L) -- that was the point of brief
 *           26_0430-266, fixing the --hmm preflight's over-estimate of
 *           do_bandedoa with the full-matrix formula meant for
 *           --hmmvit/--hmmnoband.  Modes 0 and 1 are still O(ncell)
 *           though, which is what mode 2 removes.
 *
 *           Note: the mode 0/1 estimates omit the O(nrow) P7B_GEO row
 *           tables (dominated by their O(ncell) term); they are left
 *           exactly as brief 26_0430-266 wrote them so this change
 *           cannot perturb an existing --mxsize decision.  Mode 2 counts
 *           them, since with the O(ncell) term gone they are no longer
 *           negligible.
 *
 * Args:     bnd       - the derived band (bnd->ncell, bnd->nrow, bnd->kmem set)
 *           ckpt_mode - P7B_OAMEM_NOCKPT | P7B_OAMEM_CKPT | P7B_OAMEM_CKPTPP.
 *                       (Historically an int flag do_ckpt; FALSE/TRUE still
 *                       mean modes 0/1, so old callers keep their meaning.)
 *           ret_bytes - RETURN: estimated peak bytes
 *
 * Returns:  eslOK on success.
 */
int
p7_CheckptBandedOAMemNeeded(const P7_GBANDS *bnd, int ckpt_mode, double *ret_bytes)
{
  int    *kp = bnd->kmem;
  int     nrow = bnd->nrow;
  int     maxnc = 0, r, ka, kb, nc;
  int     B, nblk;
  double  bytes;

  for (r = 0; r < nrow; r++) {
    ka = *kp++; kb = *kp++;
    nc = kb - ka + 1;
    if (nc > maxnc) maxnc = nc;
  }

  B    = (int) (sqrt((double) (nrow ? nrow : 1)) + 0.5);
  if (B < 1) B = 1;
  nblk = (nrow + B - 1) / B;

  if (ckpt_mode == P7B_OAMEM_CKPTPP) {
    /* brief 26_0628-081: no O(ncell) term at all.
     *   3 seed arrays (F, B, OA):   3 * nblk  * maxnc * p7G_NSCELLS floats
     *   F + B block buffers:        2 * (B+2) * maxnc * p7G_NSCELLS floats
     *   pp block buffer:                (B+2) * maxnc * P7B_PP_NSCELLS floats
     *   OA traceback window:            (B+2) * maxnc * p7G_NSCELLS floats
     *   pp + OA block specials:       2*(B+2) * p7G_NXCELLS floats
     *   P7B_GEO row tables:           nrow * (4 ints + 1 int64)
     *   traceback row_idx:            (L+1) ints
     */
    double nr = (double) (B + 2);
    bytes  = (double) sizeof(float) * (double) maxnc *
             ( p7G_NSCELLS   * (3.0 * (double) nblk + 3.0 * nr)
             + P7B_PP_NSCELLS * nr );
    bytes += (double) sizeof(float)   * p7G_NXCELLS * 2.0 * nr;
    bytes += (double) (4 * sizeof(int) + sizeof(int64_t)) * (double) nrow;
    bytes += (double) sizeof(int) * ((double) bnd->L + 1.0);
  }
  else if (ckpt_mode == P7B_OAMEM_CKPT) {
    bytes  = (double) sizeof(float) * ((double) bnd->ncell * P7B_PP_NSCELLS
					+ (double) bnd->nrow  * p7G_NXCELLS);
    bytes += (double) sizeof(float) * p7G_NSCELLS * (double) maxnc * (double) (nblk + B + 2);
  }
  else {
    bytes = 2.0 * (double) sizeof(float) * ((double) bnd->ncell * p7G_NSCELLS
					     + (double) bnd->nrow  * p7G_NXCELLS);
  }

  if (ret_bytes) *ret_bytes = bytes;
  return eslOK;
}

/* Per-call banded-row geometry, derived once from a P7_GBANDS. */
typedef struct {
  int      nrow;     /* number of banded rows                         */
  int      nseg;     /* number of band segments                       */
  int      M;        /* profile length                                */
  int      L;        /* sequence length                               */
  int     *irow;     /* [0..nrow-1] sequence position i of banded row */
  int     *ka;       /* [0..nrow-1] band low  k for the row           */
  int     *kb;       /* [0..nrow-1] band high k for the row           */
  int     *gap;      /* [0..nrow-1] # unbanded rows just below this    */
                     /*   row (i_r - i_{r-1} - 1; r=0: i_0 - 1).       */
                     /*   gap>0  <=>  row r starts a new segment.      */
  int64_t *dpoff;    /* [0..nrow-1] cell offset of row in a full dp    */
  int      B;        /* checkpoint block size ~ round(sqrt(nrow))     */
  int      nblk;     /* number of blocks = ceil(nrow/B)               */
  int      maxnc;    /* max (kb-ka+1) over rows (scratch sizing)       */
} P7B_GEO;

static void
p7b_geo_Destroy(P7B_GEO *g)
{
  if (g == NULL) return;
  if (g->irow)  free(g->irow);
  if (g->ka)    free(g->ka);
  if (g->kb)    free(g->kb);
  if (g->gap)   free(g->gap);
  if (g->dpoff) free(g->dpoff);
  free(g);
}

static P7B_GEO *
p7b_geo_Create(const P7_GBANDS *bnd, int M, int L)
{
  P7B_GEO *g = NULL;
  int      status;
  int     *ip = bnd->imem;
  int     *kp = bnd->kmem;
  int      nrow = bnd->nrow;
  int      seg, i, r = 0, prev_i = 0;
  int64_t  cum = 0;

  ESL_ALLOC(g, sizeof(P7B_GEO));
  g->irow = g->ka = g->kb = g->gap = NULL;
  g->dpoff = NULL;
  g->nrow = nrow; g->nseg = bnd->nseg; g->M = M; g->L = L;
  ESL_ALLOC(g->irow,  sizeof(int)     * (nrow ? nrow : 1));
  ESL_ALLOC(g->ka,    sizeof(int)     * (nrow ? nrow : 1));
  ESL_ALLOC(g->kb,    sizeof(int)     * (nrow ? nrow : 1));
  ESL_ALLOC(g->gap,   sizeof(int)     * (nrow ? nrow : 1));
  ESL_ALLOC(g->dpoff, sizeof(int64_t) * (nrow ? nrow : 1));

  g->maxnc = 0;
  for (seg = 0; seg < bnd->nseg; seg++) {
    int ia = *ip++;
    int ib = *ip++;
    for (i = ia; i <= ib; i++) {
      int ka = *kp++;
      int kb = *kp++;
      int nc = kb - ka + 1;
      g->irow[r]  = i;
      g->ka[r]    = ka;
      g->kb[r]    = kb;
      g->gap[r]   = i - prev_i - 1;     /* prev_i = 0 before first row */
      g->dpoff[r] = cum;
      cum        += (int64_t) nc * P7B_PP_NSCELLS;   /* brief 26_0526-017: dpoff indexes the compact 2-cell pp */
      if (nc > g->maxnc) g->maxnc = nc;
      prev_i = i;
      r++;
    }
  }
  g->B = (int) (sqrt((double) (nrow ? nrow : 1)) + 0.5);
  if (g->B < 1) g->B = 1;
  g->nblk = (nrow + g->B - 1) / g->B;
  return g;

 ERROR:
  p7b_geo_Destroy(g);
  return NULL;
}

/* p7b_fwd_row(): fill one banded Forward row.  Verbatim copy of the
 * per-row body of my_p7_GForwardBanded() (k-loop + special updates).
 *
 *   dpc  - OUT: row's dp cells, p7G_NSCELLS*(kb-ka+1) floats, k ascending
 *   dpp  - prev row's dp cells (k=kap..kbp), or NULL at a segment start
 *   xB_in- B-state value entering this row
 *   x{N,J,C}_in - N/J/C running specials entering this row
 *   o_x* - OUT: specials after this row (xE,xN,xJ,xB,xC)
 */
static void
p7b_fwd_row(const ESL_DSQ *dsq, const P7_PROFILE *gm, int M, float esc,
            int i, int kac, int kbc, const float *dpp, int kap, int kbp,
            float xB, float xN_in, float xJ_in, float xC_in,
            float *dpc,
            float *o_xE, float *o_xN, float *o_xJ, float *o_xB, float *o_xC)
{
  float const *tsc = gm->tsc;
  float const *rsc = gm->rsc[dsq[i]];
  float        mvp, ivp, dvp, dc, sc, xE;
  int          k;
  int          kbc2 = (kbc == M ? kbc - 1 : kbc);

  dc = -eslINFINITY;
  xE = -eslINFINITY;

  /* advance dpp by any left overhang of the previous row */
  if (dpp != NULL) dpp += (kac-1 > kap ? ESL_MIN(kac-kap-1, kbp-kap+1) * p7G_NSCELLS : 0);

  if (dpp != NULL && kac > kap && kac-1 <= kbp) { mvp = *dpp++; ivp = *dpp++; dvp = *dpp++; }
  else                                          { mvp = -eslINFINITY; ivp = -eslINFINITY; dvp = -eslINFINITY; }

  for (k = kac; k <= kbc2; k++)
    {
      *dpc++ = sc = MSC(k) + p7_FLogsum( p7_FLogsum(mvp + TSC(p7P_MM, k-1), ivp + TSC(p7P_IM, k-1)),
                                         p7_FLogsum(dvp + TSC(p7P_DM, k-1), xB  + TSC(p7P_BM, k-1)));

      if (dpp != NULL && k >= kap && k <= kbp) { mvp = *dpp++; ivp = *dpp++; dvp = *dpp++; }
      else                                     { mvp = -eslINFINITY; ivp = -eslINFINITY; dvp = -eslINFINITY; }

      *dpc++ = ISC(k) + p7_FLogsum( mvp + TSC(p7P_MI, k), ivp + TSC(p7P_II, k));

      xE     = p7_FLogsum( p7_FLogsum(sc + esc, dc + esc), xE);

      *dpc++ = dc;
      dc     = p7_FLogsum( sc + TSC(p7P_MD, k), dc + TSC(p7P_DD, k));
    }

  if (kbc2 < kbc) /* kbc==M: unrolled final column */
    {
      *dpc++ = sc = MSC(k) + p7_FLogsum( p7_FLogsum(mvp + TSC(p7P_MM, k-1), ivp + TSC(p7P_IM, k-1)),
                                         p7_FLogsum(dvp + TSC(p7P_DM, k-1), xB  + TSC(p7P_BM, k-1)));
      *dpc++ = -eslINFINITY;
      *dpc++ = dc;
      xE     = p7_FLogsum( p7_FLogsum(sc, dc), xE);
    }

  *o_xE = xE;
  *o_xN = xN_in + gm->xsc[p7P_N][p7P_LOOP];
  *o_xJ = p7_FLogsum( xJ_in + gm->xsc[p7P_J][p7P_LOOP], xE + gm->xsc[p7P_E][p7P_LOOP]);
  *o_xB = p7_FLogsum( *o_xJ + gm->xsc[p7P_J][p7P_MOVE],  *o_xN + gm->xsc[p7P_N][p7P_MOVE]);
  *o_xC = p7_FLogsum( xE + gm->xsc[p7P_E][p7P_MOVE],     xC_in + gm->xsc[p7P_C][p7P_LOOP]);
}

/* A Forward checkpoint seed: the state ENTERING a block (= the state
 * after the last row of the previous block).  block b (b>=1) starts
 * from seed[b]; block 0 starts from the initial state.
 */
typedef struct {
  float  *dp;        /* copy of the previous row's dp cells (NULL if none) */
  int     ka, kb;    /* that row's band (kap,kbp for the next row)         */
  float   xN, xJ, xC, xB;  /* running specials after that row             */
  int     valid;
} P7B_FSEED;

/* p7b_forward_seeds(): STEP A.  Stream banded Forward once, store the
 * entering state for each block, and return the Forward score.
 */
static int
p7b_forward_seeds(const ESL_DSQ *dsq, const P7_PROFILE *gm, const P7B_GEO *g,
                  P7B_FSEED *seed, float *ret_fwdsc)
{
  int    status;
  int    M = g->M, L = g->L, B = g->B;
  float  esc = p7_profile_IsLocal(gm) ? 0 : -eslINFINITY;
  float *buf0 = NULL, *buf1 = NULL;     /* prev / cur row dp scratch */
  float *prev = NULL, *cur;
  int    r, kap = M+1, kbp = M+1;
  float  xN = 0.0f, xJ = -eslINFINITY, xC = -eslINFINITY, xB;
  float  xE_o, xN_o, xJ_o, xB_o, xC_o;

  ESL_ALLOC(buf0, sizeof(float) * (g->maxnc * p7G_NSCELLS + 1));
  ESL_ALLOC(buf1, sizeof(float) * (g->maxnc * p7G_NSCELLS + 1));

  xB = -eslINFINITY;
  for (r = 0; r < g->nrow; r++)
    {
      int   i   = g->irow[r];
      int   kac = g->ka[r], kbc = g->kb[r];
      const float *dpp;
      int   seg_start = (g->gap[r] > 0) || (r == 0);

      cur = (prev == buf0) ? buf1 : buf0;

      if (seg_start) {
        if (g->gap[r] > 0) {
          xN += g->gap[r] * gm->xsc[p7P_N][p7P_LOOP];
          xJ += g->gap[r] * gm->xsc[p7P_J][p7P_LOOP];
          xC += g->gap[r] * gm->xsc[p7P_C][p7P_LOOP];
        }
        xB  = p7_FLogsum( xN + gm->xsc[p7P_N][p7P_MOVE], xJ + gm->xsc[p7P_J][p7P_MOVE]);
        dpp = NULL;
      } else {
        /* continue segment: xN/xJ/xC/xB already hold prev row's output */
        dpp = prev;
      }

      p7b_fwd_row(dsq, gm, M, esc, i, kac, kbc, dpp, kap, kbp,
                  xB, xN, xJ, xC, cur,
                  &xE_o, &xN_o, &xJ_o, &xB_o, &xC_o);

      xN = xN_o; xJ = xJ_o; xC = xC_o; xB = xB_o;
      kap = kac; kbp = kbc;
      prev = cur;

      /* snapshot a seed if the NEXT row begins a new block */
      if (r + 1 < g->nrow && ((r + 1) % B) == 0) {
        int b = (r + 1) / B;
        int nc = kbc - kac + 1;
        ESL_ALLOC(seed[b].dp, sizeof(float) * nc * p7G_NSCELLS);
        memcpy(seed[b].dp, cur, sizeof(float) * nc * p7G_NSCELLS);
        seed[b].ka = kac; seed[b].kb = kbc;
        seed[b].xN = xN;  seed[b].xJ = xJ; seed[b].xC = xC; seed[b].xB = xB;
        seed[b].valid = 1;
      }
    }

  /* Forward score: tail last_ib+1..L runs through C only. */
  { int last_ib = (g->nrow > 0) ? g->irow[g->nrow-1] : 0;
    int tail    = L - last_ib;
    *ret_fwdsc  = (tail > 0 ? xC + tail * gm->xsc[p7P_C][p7P_LOOP] : xC) + gm->xsc[p7P_C][p7P_MOVE];
  }

  free(buf0); free(buf1);
  return eslOK;

 ERROR:
  if (buf0) free(buf0);
  if (buf1) free(buf1);
  return status;
}

/* p7b_fwd_block(): recompute Forward rows [lo..hi] from seed[b], into
 * caller buffers: fdp[r-lo] points into a packed block buffer, and the
 * per-row running specials are returned in fxN/fxJ/fxC (specials AFTER
 * each row).  Also returns, in *s_xN/*s_xJ/*s_xC, the specials of the
 * row just below lo (i.e. the seed's specials) for the decoder's
 * fwd(i-1) at i=lo.  blkbuf must hold the block's cells contiguously.
 */
static void
p7b_fwd_block(const ESL_DSQ *dsq, const P7_PROFILE *gm, const P7B_GEO *g,
              const P7B_FSEED *seed, int b, int lo, int hi,
              float *blkbuf, int64_t *roff,
              float *fxN, float *fxJ, float *fxC,
              float *s_xN, float *s_xJ, float *s_xC)
{
  int    M = g->M;
  float  esc = p7_profile_IsLocal(gm) ? 0 : -eslINFINITY;
  int    r, kap, kbp;
  float  xN, xJ, xC, xB = -eslINFINITY;
  int64_t cum = 0;

  if (b == 0) { xN = 0.0f; xJ = -eslINFINITY; xC = -eslINFINITY; kap = kbp = M+1; }
  else        { xN = seed[b].xN; xJ = seed[b].xJ; xC = seed[b].xC; xB = seed[b].xB;
                kap = seed[b].ka; kbp = seed[b].kb; }

  *s_xN = xN; *s_xJ = xJ; *s_xC = xC;   /* seed specials = row (lo-1) */

  for (r = lo; r <= hi; r++)
    {
      int   i   = g->irow[r];
      int   kac = g->ka[r], kbc = g->kb[r];
      float xE_o, xN_o, xJ_o, xB_o, xC_o;
      float *dpc = blkbuf + cum;
      const float *dpp;
      int   no_prev   = (r == 0);                 /* global row 0 only */
      int   seg_start = (g->gap[r] > 0) || no_prev;

      if (seg_start) {
        if (g->gap[r] > 0) {
          xN += g->gap[r] * gm->xsc[p7P_N][p7P_LOOP];
          xJ += g->gap[r] * gm->xsc[p7P_J][p7P_LOOP];
          xC += g->gap[r] * gm->xsc[p7P_C][p7P_LOOP];
        }
        xB  = p7_FLogsum( xN + gm->xsc[p7P_N][p7P_MOVE], xJ + gm->xsc[p7P_J][p7P_MOVE]);
        dpp = NULL;
      } else {
        /* continue: xN/xJ/xC/xB hold prev row's output (seed for r==lo) */
        dpp = (r == lo) ? seed[b].dp : (blkbuf + roff[r-1-lo]);
      }

      p7b_fwd_row(dsq, gm, M, esc, i, kac, kbc, dpp, kap, kbp,
                  xB, xN, xJ, xC, dpc,
                  &xE_o, &xN_o, &xJ_o, &xB_o, &xC_o);

      roff[r-lo] = cum;
      fxN[r-lo] = xN_o; fxJ[r-lo] = xJ_o; fxC[r-lo] = xC_o;
      xN = xN_o; xJ = xJ_o; xC = xC_o; xB = xB_o;
      kap = kac; kbp = kbc;
      cum += (int64_t)(kbc - kac + 1) * p7G_NSCELLS;
    }
}

/* p7b_bwd_row(): fill one banded Backward row.  Verbatim copy of the
 * per-row body of p7_GBackwardBanded() (both the i==L special-case and
 * the normal i<L recurrence), writing dp cells in ascending-k [M,I,D]
 * layout.  Carried specials (xC,xJ,xN,xE,xB) are updated in place.
 *
 *   dpn       - next row (i+1) dp cells (k=kan..kbn), unused if segment top
 *   kan,kbn   - next row band (M+1,M+1 if next row outside band)
 *   dpc       - OUT: this row's dp cells (k ascending, [M,I,D])
 */
static void
p7b_bwd_row(const ESL_DSQ *dsq, const P7_PROFILE *gm, int M, float esc,
            int i, int L, int kac, int kbc,
            const float *dpn, int kan, int kbn,
            float *xC_io, float *xJ_io, float *xN_io, float *xE_io, float *xB_io,
            float *dpc)
{
  float const *tsc = gm->tsc;
  float        xC = *xC_io, xJ = *xJ_io, xN = *xN_io, xE = *xE_io, xB;
  float        dc;
  int          k;
  int          kbc2 = (kbc == M ? kbc - 1 : kbc);

  if (i == L)
    {
      /* row L special case (mirror lines ~6581-6622) */
      xB = -eslINFINITY; xJ = -eslINFINITY; xN = -eslINFINITY;  /* xC,xE carried */

      if (kbc == M) { int o=(M-kac)*p7G_NSCELLS; dpc[o]=xE; dpc[o+1]=-eslINFINITY; dpc[o+2]=xE; }
      dc = (kbc == M) ? xE : -eslINFINITY;
      for (k = kbc2; k >= kac; k--) {
        float dk = p7_FLogsum(xE + esc, dc + TSC(p7P_DD, k));
        float mk = p7_FLogsum(xE + esc, dc + TSC(p7P_MD, k));
        int   o  = (k-kac)*p7G_NSCELLS;
        dpc[o]   = mk; dpc[o+1] = -eslINFINITY; dpc[o+2] = dk;
        dc = dk;
      }
    }
  else
    {
      float const *rsc = gm->rsc[dsq[i+1]];

      xB = -eslINFINITY;
      for (k = ESL_MAX(1, kan); k <= ESL_MIN(M, kbn); k++) {
        int   off    = (k - kan) * p7G_NSCELLS;
        float m_next = dpn[off];
        xB = p7_FLogsum(xB, m_next + TSC(p7P_BM, k-1) + MSC(k));
      }

      xJ = p7_FLogsum(xJ + gm->xsc[p7P_J][p7P_LOOP], xB + gm->xsc[p7P_J][p7P_MOVE]);
      xC = xC + gm->xsc[p7P_C][p7P_LOOP];
      xE = p7_FLogsum(xJ + gm->xsc[p7P_E][p7P_LOOP], xC + gm->xsc[p7P_E][p7P_MOVE]);
      xN = p7_FLogsum(xN + gm->xsc[p7P_N][p7P_LOOP], xB + gm->xsc[p7P_N][p7P_MOVE]);

      if (kbc == M) { int o=(M-kac)*p7G_NSCELLS; dpc[o]=xE; dpc[o+1]=-eslINFINITY; dpc[o+2]=xE; }
      dc = (kbc == M) ? xE : -eslINFINITY;

      for (k = kbc2; k >= kac; k--)
        {
          float mnext, inext, dnext;
          if (k+1 >= kan && k+1 <= kbn) { int o=(k+1-kan)*p7G_NSCELLS; mnext = dpn[o] + MSC(k+1); dnext = dpn[o+2]; }
          else                          { mnext = dnext = -eslINFINITY; }
          if (k   >= kan && k   <= kbn) { int o=(k-kan)*p7G_NSCELLS;   inext = dpn[o+1] + ISC(k); }
          else                          { inext = -eslINFINITY; }

          float dk = p7_FLogsum(p7_FLogsum(mnext + TSC(p7P_DM, k), dc + TSC(p7P_DD, k)), xE + esc);
          float ik = p7_FLogsum(mnext + TSC(p7P_IM, k), inext + TSC(p7P_II, k));
          float mk = p7_FLogsum(p7_FLogsum(mnext + TSC(p7P_MM, k), inext + TSC(p7P_MI, k)),
                                p7_FLogsum(xE + esc, dc + TSC(p7P_MD, k)));
          int   o  = (k-kac)*p7G_NSCELLS;
          dpc[o] = mk; dpc[o+1] = ik; dpc[o+2] = dk;
          dc = dk;
          (void) dnext;
        }
    }

  *xC_io = xC; *xJ_io = xJ; *xN_io = xN; *xE_io = xE; *xB_io = xB;
}

/* p7b_decode_row(): fold one row's posterior into the caller's output
 * row buffers.  Verbatim copy of the per-row body of
 * p7_GDecodingBanded() (expf + per-row denom normalization + NaN
 * guards).
 *
 * brief 26_0628-081: the destination is now passed as explicit row
 * pointers <pdp>/<pxp> rather than (pp, r).  The resident-pp caller
 * passes pp->dp + g->dpoff[r] / pp->xmx + r*p7G_NXCELLS -- identical to
 * what this function used to compute itself; the double-checkpointed
 * caller passes a one-block posterior buffer.  No arithmetic changed. */
static void
p7b_decode_row(const P7B_GEO *g, const P7_PROFILE *gm, float fwdsc,
               int kac, int kbc,
               const float *fdp, const float *bdp,
               float bck_xN, float bck_xJ, float bck_xC,
               float fwd_xN_prev, float fwd_xJ_prev, float fwd_xC_prev,
               float *pdp, float *pxp)
{
  int    M = g->M;
  float  denom = 0.0f;
  int    k;

  for (k = kac; k <= kbc; k++) {
    int off  = (k - kac) * p7G_NSCELLS;     /* fdp/bdp: full 3-cell Forward/Backward rows */
    int poff = (k - kac) * P7B_PP_NSCELLS;  /* brief 26_0526-017: pp is compact 2-cell (M,I)      */
    pdp[poff] = expf(fdp[off] + bdp[off] - fwdsc);          /* M */
    if (! isfinite(pdp[poff])) pdp[poff] = 0.0f;  /* guard NaN from -inf + inf (cell unreachable in one direction) */
    denom += pdp[poff];
    if (k < M) {
      pdp[poff+1] = expf(fdp[off+1] + bdp[off+1] - fwdsc);  /* I */
      if (! isfinite(pdp[poff+1])) pdp[poff+1] = 0.0f;      /* guard NaN from -inf + inf */
      denom += pdp[poff+1];
    }
    else         pdp[poff+1] = 0.0f;
    /* pp_D dropped (brief 26_0526-017): was always 0 and never read downstream */
  }

  pxp[p7G_E] = 0.0f;

  pxp[p7G_N] = expf(fwd_xN_prev + bck_xN + gm->xsc[p7P_N][p7P_LOOP] - fwdsc);
  if (! isfinite(pxp[p7G_N])) pxp[p7G_N] = 0.0f;
  denom += pxp[p7G_N];

  pxp[p7G_J] = expf(fwd_xJ_prev + bck_xJ + gm->xsc[p7P_J][p7P_LOOP] - fwdsc);
  if (! isfinite(pxp[p7G_J])) pxp[p7G_J] = 0.0f;
  denom += pxp[p7G_J];

  pxp[p7G_B] = 0.0f;

  pxp[p7G_C] = expf(fwd_xC_prev + bck_xC + gm->xsc[p7P_C][p7P_LOOP] - fwdsc);
  if (! isfinite(pxp[p7G_C])) pxp[p7G_C] = 0.0f;
  denom += pxp[p7G_C];

  if (denom > 0.0f) {
    denom = 1.0f / denom;
    for (k = kac; k <= kbc; k++) { int poff = (k - kac) * P7B_PP_NSCELLS; pdp[poff] *= denom; pdp[poff+1] *= denom; }
    pxp[p7G_N] *= denom;
    pxp[p7G_J] *= denom;
    pxp[p7G_C] *= denom;
  }
}

/* p7b_backdecode(): STEP B.  Stream Backward (2 rows live) in descending
 * blocks, recompute each block's Forward from seeds, fold the posterior
 * into the resident <pp>.  The full Backward matrix is never stored.
 */
static int
p7b_backdecode(const ESL_DSQ *dsq, const P7_PROFILE *gm, const P7B_GEO *g,
               const P7B_FSEED *seed, float fwdsc, P7_GMXB *pp)
{
  int      status;
  int      M = g->M, L = g->L, B = g->B;
  float    esc = p7_profile_IsLocal(gm) ? 0 : -eslINFINITY;
  float   *blkbuf = NULL, *bbuf0 = NULL, *bbuf1 = NULL;
  int64_t *roff   = NULL;
  float   *fxN = NULL, *fxJ = NULL, *fxC = NULL;
  int      b, r;
  /* backward carry */
  float    xC, xJ, xN, xE, xB;
  const float *dpn = NULL;
  int      kan = M+1, kbn = M+1;
  float   *bnext = NULL, *bcur;

  ESL_ALLOC(blkbuf, sizeof(float) * ((int64_t) B * g->maxnc * p7G_NSCELLS + 1));
  ESL_ALLOC(roff,   sizeof(int64_t) * B);
  ESL_ALLOC(fxN,    sizeof(float) * B);
  ESL_ALLOC(fxJ,    sizeof(float) * B);
  ESL_ALLOC(fxC,    sizeof(float) * B);
  ESL_ALLOC(bbuf0,  sizeof(float) * (g->maxnc * p7G_NSCELLS + 1));
  ESL_ALLOC(bbuf1,  sizeof(float) * (g->maxnc * p7G_NSCELLS + 1));

  /* init backward specials (mirror p7_GBackwardBanded init) */
  xC = gm->xsc[p7P_C][p7P_MOVE];
  xE = xC + gm->xsc[p7P_E][p7P_MOVE];
  xJ = xB = xN = -eslINFINITY;

  for (b = g->nblk - 1; b >= 0; b--)
    {
      int   lo = b * B;
      int   hi = ESL_MIN((b+1)*B - 1, g->nrow - 1);
      float s_xN, s_xJ, s_xC;

      p7b_fwd_block(dsq, gm, g, seed, b, lo, hi, blkbuf, roff, fxN, fxJ, fxC, &s_xN, &s_xJ, &s_xC);

      for (r = hi; r >= lo; r--)
        {
          int   i   = g->irow[r];
          int   kac = g->ka[r], kbc = g->kb[r];
          int   top = (r == g->nrow - 1) || (g->gap[r+1] > 0);
          const float *use_dpn = dpn;
          int   use_kan = kan, use_kbn = kbn;
          float fwd_xN_prev, fwd_xJ_prev, fwd_xC_prev;

          if (top) {
            int gap_above = (r == g->nrow - 1) ? (L - i) : g->gap[r+1];
            if (gap_above > 0) { xC += gap_above * gm->xsc[p7P_C][p7P_LOOP];
                                 xE  = p7_FLogsum(xE, xC + gm->xsc[p7P_E][p7P_MOVE]); }
            use_dpn = NULL; use_kan = M+1; use_kbn = M+1;
          }

          bcur = (bnext == bbuf0) ? bbuf1 : bbuf0;
          p7b_bwd_row(dsq, gm, M, esc, i, L, kac, kbc, use_dpn, use_kan, use_kbn,
                      &xC, &xJ, &xN, &xE, &xB, bcur);

          /* fwd specials of row (r-1), gap-adjusted, for the decoder */
          if (r == lo) { fwd_xN_prev = s_xN; fwd_xJ_prev = s_xJ; fwd_xC_prev = s_xC; }
          else         { fwd_xN_prev = fxN[r-1-lo]; fwd_xJ_prev = fxJ[r-1-lo]; fwd_xC_prev = fxC[r-1-lo]; }
          if (g->gap[r] > 0) {
            fwd_xN_prev += g->gap[r] * gm->xsc[p7P_N][p7P_LOOP];
            fwd_xJ_prev += g->gap[r] * gm->xsc[p7P_J][p7P_LOOP];
            fwd_xC_prev += g->gap[r] * gm->xsc[p7P_C][p7P_LOOP];
          }

          p7b_decode_row(g, gm, fwdsc, kac, kbc, blkbuf + roff[r-lo], bcur,
                         xN, xJ, xC, fwd_xN_prev, fwd_xJ_prev, fwd_xC_prev,
                         pp->dp + g->dpoff[r], pp->xmx + (int64_t) r * p7G_NXCELLS);

          dpn = bcur; kan = kac; kbn = kbc; bnext = bcur;
        }
    }

  free(blkbuf); free(roff); free(fxN); free(fxJ); free(fxC); free(bbuf0); free(bbuf1);
  return eslOK;

 ERROR:
  if (blkbuf) free(blkbuf);
  if (roff)   free(roff);
  if (fxN)    free(fxN);
  if (fxJ)    free(fxJ);
  if (fxC)    free(fxC);
  if (bbuf0)  free(bbuf0);
  if (bbuf1)  free(bbuf1);
  return status;
}

/*****************************************************************
 * brief 26_0628-081: DOUBLE-CHECKPOINTING the Backward pass, so the
 * posterior <pp> never has to be resident.
 *
 * 26_0526-016/017 checkpointed the *compute* of this engine (F, B,
 * Decode and the OA fill all run out of ~sqrt(nrow) seed rows), but one
 * O(bnd->ncell) matrix survived: the resident posterior that
 * p7b_backdecode() fills and that the OA pass + traceback consume.  It
 * had to be resident because of a DIRECTION CONFLICT -- pp rows are
 * produced by Backward in DESCENDING row order, and consumed by OA in
 * ASCENDING row order.
 *
 * The fix mirrors what the Forward pass already does for itself:
 * checkpoint Backward too.  A second seed array (P7B_BSEED) stores the
 * Backward state ENTERING each block from above; a block's Backward
 * rows are then recomputed on demand into an O(B*maxnc) buffer, the pp
 * rows for that block are formed on the fly, and OA consumes them
 * ascending.  The traceback recomputes the same per-block pp window
 * alongside the OA window it already recomputed.
 *
 * Resident term drops from O(ncell) to O(sqrt(nrow)*maxnc); the cost is
 * two extra O(ncell) recompute passes (one Backward stream to lay the
 * seeds, plus per-block F+B recompute during OA/traceback).
 *
 * Byte-exact discipline, same as 26_0526-016: no per-cell arithmetic is
 * touched.  A Backward block resumed from a seed replays exactly the
 * same p7b_bwd_row() calls, on exactly the same float inputs, in
 * exactly the same order, as the uninterrupted stream would have -- so
 * the recomputed pp rows are bit-identical to the resident ones, and
 * everything downstream is bit-identical too.  The legacy resident-pp
 * path is left fully intact and reachable via INFERNAL_HMM_PPCKPT_OFF.
 *****************************************************************/

/* A Backward checkpoint seed: the Backward state ENTERING block <b>
 * from above, i.e. the state just after row (hi_b + 1) was processed.
 * For the last block there is no row above, so the seed is the Backward
 * initial state (dp == NULL).
 *
 * Note that <dp> doubles as the stored copy of row (hi_b+1)'s Backward
 * row, which the traceback's B+1-row window needs verbatim (see
 * p7b_bwd_block()).  xB is not carried: p7b_bwd_row() recomputes it
 * from scratch every row and never reads an incoming value.
 */
typedef struct {
  float  *dp;              /* copy of row (hi_b+1)'s Backward dp cells, or NULL */
  int     ka, kb;          /* that row's band (kan,kbn for the row below)       */
  float   xC, xJ, xN, xE;  /* carried specials AFTER that row                   */
  int     valid;
} P7B_BSEED;

/* p7b_backward_seeds(): STEP B1.  Stream banded Backward once (2 rows
 * live, no posterior formed, no Forward needed), storing the entering
 * state for each block.  Descending-row mirror of p7b_forward_seeds().
 */
static int
p7b_backward_seeds(const ESL_DSQ *dsq, const P7_PROFILE *gm, const P7B_GEO *g,
                   P7B_BSEED *bseed)
{
  int      status;
  int      M = g->M, L = g->L, B = g->B;
  float    esc = p7_profile_IsLocal(gm) ? 0 : -eslINFINITY;
  float   *bbuf0 = NULL, *bbuf1 = NULL, *bnext = NULL, *bcur;
  const float *dpn = NULL;
  int      kan = M+1, kbn = M+1;
  float    xC, xJ, xN, xE, xB;
  int      r;

  ESL_ALLOC(bbuf0, sizeof(float) * (g->maxnc * p7G_NSCELLS + 1));
  ESL_ALLOC(bbuf1, sizeof(float) * (g->maxnc * p7G_NSCELLS + 1));

  /* init backward specials (mirror p7b_backdecode / p7_GBackwardBanded init) */
  xC = gm->xsc[p7P_C][p7P_MOVE];
  xE = xC + gm->xsc[p7P_E][p7P_MOVE];
  xJ = xB = xN = -eslINFINITY;

  /* the topmost block enters from the Backward initial state.  (nblk==0
   * only when the band has no rows at all, in which case the loop below
   * is empty too and there is no block to seed.) */
  if (g->nblk > 0) {
    bseed[g->nblk-1].dp = NULL;
    bseed[g->nblk-1].ka = M+1;  bseed[g->nblk-1].kb = M+1;
    bseed[g->nblk-1].xC = xC;   bseed[g->nblk-1].xJ = xJ;
    bseed[g->nblk-1].xN = xN;   bseed[g->nblk-1].xE = xE;
    bseed[g->nblk-1].valid = 1;
  }

  for (r = g->nrow - 1; r >= 0; r--)
    {
      int   i   = g->irow[r];
      int   kac = g->ka[r], kbc = g->kb[r];
      int   top = (r == g->nrow - 1) || (g->gap[r+1] > 0);
      const float *use_dpn = dpn;
      int   use_kan = kan, use_kbn = kbn;

      if (top) {
        int gap_above = (r == g->nrow - 1) ? (L - i) : g->gap[r+1];
        if (gap_above > 0) { xC += gap_above * gm->xsc[p7P_C][p7P_LOOP];
                             xE  = p7_FLogsum(xE, xC + gm->xsc[p7P_E][p7P_MOVE]); }
        use_dpn = NULL; use_kan = M+1; use_kbn = M+1;
      }

      bcur = (bnext == bbuf0) ? bbuf1 : bbuf0;
      p7b_bwd_row(dsq, gm, M, esc, i, L, kac, kbc, use_dpn, use_kan, use_kbn,
                  &xC, &xJ, &xN, &xE, &xB, bcur);

      /* snapshot a seed if row r is the first row of a block: it is then
       * the row just above block (r/B - 1)'s last row. */
      if (r > 0 && (r % B) == 0) {
        int b = r / B - 1, nc = kbc - kac + 1;
        ESL_ALLOC(bseed[b].dp, sizeof(float) * nc * p7G_NSCELLS);
        memcpy(bseed[b].dp, bcur, sizeof(float) * nc * p7G_NSCELLS);
        bseed[b].ka = kac; bseed[b].kb = kbc;
        bseed[b].xC = xC;  bseed[b].xJ = xJ; bseed[b].xN = xN; bseed[b].xE = xE;
        bseed[b].valid = 1;
      }

      dpn = bcur; kan = kac; kbn = kbc; bnext = bcur;
    }

  free(bbuf0); free(bbuf1);
  return eslOK;

 ERROR:
  if (bbuf0) free(bbuf0);
  if (bbuf1) free(bbuf1);
  return status;
}

/* p7b_bwd_block(): recompute Backward rows [lo..hi] of block <b> from
 * bseed[b], into a packed block buffer (offsets in <roff>), plus the
 * per-row carried specials AFTER each row (bxN/bxJ/bxC) that
 * p7b_decode_row() wants.  Descending-row mirror of p7b_fwd_block().
 *
 * <lo> is always b*B.  <hi> may be hi_b (= the block's last row) or
 * hi_b+1 (the traceback window's +1 overlap row).  In the latter case
 * that extra row IS the seed row, so it is copied verbatim from
 * bseed[b].dp rather than recomputed -- which is also what makes it
 * bit-identical rather than merely equal.
 */
static void
p7b_bwd_block(const ESL_DSQ *dsq, const P7_PROFILE *gm, const P7B_GEO *g,
              const P7B_BSEED *bseed, int b, int lo, int hi,
              float *blkbuf, int64_t *roff, float *bxN, float *bxJ, float *bxC)
{
  int      M = g->M, L = g->L;
  float    esc = p7_profile_IsLocal(gm) ? 0 : -eslINFINITY;
  float    xC, xJ, xN, xE, xB = -eslINFINITY;
  const float *dpn;
  int      kan, kbn;
  int      r, rstart;
  int      hib = ESL_MIN((b+1) * g->B - 1, g->nrow - 1);
  int64_t  cum = 0;

  for (r = lo; r <= hi; r++) {
    roff[r-lo] = cum;
    cum += (int64_t)(g->kb[r] - g->ka[r] + 1) * p7G_NSCELLS;
  }

  xC  = bseed[b].xC; xJ = bseed[b].xJ; xN = bseed[b].xN; xE = bseed[b].xE;
  dpn = bseed[b].dp; kan = bseed[b].ka; kbn = bseed[b].kb;

  rstart = hi;
  if (hi > hib) {   /* the +1 overlap row is exactly the seed row */
    int nc = g->kb[hi] - g->ka[hi] + 1;
    memcpy(blkbuf + roff[hi-lo], bseed[b].dp, sizeof(float) * nc * p7G_NSCELLS);
    bxN[hi-lo] = xN; bxJ[hi-lo] = xJ; bxC[hi-lo] = xC;
    rstart = hib;
  }

  for (r = rstart; r >= lo; r--)
    {
      int   i   = g->irow[r];
      int   kac = g->ka[r], kbc = g->kb[r];
      int   top = (r == g->nrow - 1) || (g->gap[r+1] > 0);
      const float *use_dpn = dpn;
      int   use_kan = kan, use_kbn = kbn;
      float *bcur = blkbuf + roff[r-lo];

      if (top) {
        int gap_above = (r == g->nrow - 1) ? (L - i) : g->gap[r+1];
        if (gap_above > 0) { xC += gap_above * gm->xsc[p7P_C][p7P_LOOP];
                             xE  = p7_FLogsum(xE, xC + gm->xsc[p7P_E][p7P_MOVE]); }
        use_dpn = NULL; use_kan = M+1; use_kbn = M+1;
      }

      p7b_bwd_row(dsq, gm, M, esc, i, L, kac, kbc, use_dpn, use_kan, use_kbn,
                  &xC, &xJ, &xN, &xE, &xB, bcur);

      bxN[r-lo] = xN; bxJ[r-lo] = xJ; bxC[r-lo] = xC;
      dpn = bcur; kan = kac; kbn = kbc;
    }
}

/* A posterior-row source for the OA pass and traceback: either the
 * legacy resident matrix (<pp> non-NULL) or a recomputed one-block
 * window (<pp> NULL, rows [lo..hi] held in <ppblk>/<ppxmx>). */
typedef struct {
  const P7_GMXB *pp;       /* legacy resident pp, or NULL             */
  const float   *ppblk;    /* window mode: packed pp rows             */
  const int64_t *pproff;   /* window mode: per-row offsets into ppblk */
  const float   *ppxmx;    /* window mode: per-row specials           */
  int            lo, hi;   /* window mode: rows held                  */
} P7B_PPVIEW;

static const float *
p7b_ppview_dp(const P7B_PPVIEW *v, const P7B_GEO *g, int r)
{
  return (v->pp != NULL) ? v->pp->dp + g->dpoff[r]
                         : v->ppblk  + v->pproff[r - v->lo];
}

static const float *
p7b_ppview_xp(const P7B_PPVIEW *v, const P7B_GEO *g, int r)
{
  return (v->pp != NULL) ? v->pp->xmx + (int64_t) r * p7G_NXCELLS
                         : v->ppxmx   + (int64_t)(r - v->lo) * p7G_NXCELLS;
}

/* Scratch for recomputing one block's posterior rows on demand.  Every
 * buffer here is O(B*maxnc) or smaller -- this struct is what replaces
 * the O(ncell) resident pp. */
typedef struct {
  const ESL_DSQ   *dsq;
  const P7B_FSEED *fseed;   /* Forward  seeds (from p7b_forward_seeds)  */
  const P7B_BSEED *bseed;   /* Backward seeds (from p7b_backward_seeds) */
  float            fwdsc;
  int              nr;      /* rows a block window can hold (B+2)       */
  float   *fblk;  int64_t *froff;  float *fxN, *fxJ, *fxC;
  float   *bblk;  int64_t *broff;  float *bxN, *bxJ, *bxC;
  float   *ppblk; int64_t *pproff; float *ppxmx;
} P7B_PPWORK;

static void
p7b_ppwork_Destroy(P7B_PPWORK *w)
{
  if (w == NULL) return;
  if (w->fblk)   free(w->fblk);
  if (w->froff)  free(w->froff);
  if (w->fxN)    free(w->fxN);
  if (w->fxJ)    free(w->fxJ);
  if (w->fxC)    free(w->fxC);
  if (w->bblk)   free(w->bblk);
  if (w->broff)  free(w->broff);
  if (w->bxN)    free(w->bxN);
  if (w->bxJ)    free(w->bxJ);
  if (w->bxC)    free(w->bxC);
  if (w->ppblk)  free(w->ppblk);
  if (w->pproff) free(w->pproff);
  if (w->ppxmx)  free(w->ppxmx);
  free(w);
}

static P7B_PPWORK *
p7b_ppwork_Create(const ESL_DSQ *dsq, const P7B_GEO *g,
                  const P7B_FSEED *fseed, const P7B_BSEED *bseed, float fwdsc)
{
  P7B_PPWORK *w = NULL;
  int         status;
  int         nr = g->B + 2;   /* block (B rows) + traceback's +1 overlap row + slack */

  ESL_ALLOC(w, sizeof(P7B_PPWORK));
  w->fblk = w->fxN = w->fxJ = w->fxC = NULL;
  w->bblk = w->bxN = w->bxJ = w->bxC = NULL;
  w->ppblk = w->ppxmx = NULL;
  w->froff = w->broff = w->pproff = NULL;
  w->dsq = dsq; w->fseed = fseed; w->bseed = bseed; w->fwdsc = fwdsc; w->nr = nr;

  ESL_ALLOC(w->fblk,   sizeof(float)   * ((int64_t) nr * g->maxnc * p7G_NSCELLS + 1));
  ESL_ALLOC(w->froff,  sizeof(int64_t) * nr);
  ESL_ALLOC(w->fxN,    sizeof(float)   * nr);
  ESL_ALLOC(w->fxJ,    sizeof(float)   * nr);
  ESL_ALLOC(w->fxC,    sizeof(float)   * nr);
  ESL_ALLOC(w->bblk,   sizeof(float)   * ((int64_t) nr * g->maxnc * p7G_NSCELLS + 1));
  ESL_ALLOC(w->broff,  sizeof(int64_t) * nr);
  ESL_ALLOC(w->bxN,    sizeof(float)   * nr);
  ESL_ALLOC(w->bxJ,    sizeof(float)   * nr);
  ESL_ALLOC(w->bxC,    sizeof(float)   * nr);
  ESL_ALLOC(w->ppblk,  sizeof(float)   * ((int64_t) nr * g->maxnc * P7B_PP_NSCELLS + 1));
  ESL_ALLOC(w->pproff, sizeof(int64_t) * nr);
  ESL_ALLOC(w->ppxmx,  sizeof(float)   * ((int64_t) nr * p7G_NXCELLS));
  return w;

 ERROR:
  p7b_ppwork_Destroy(w);
  return NULL;
}

/* p7b_ppwork_fill(): recompute the posterior rows [lo..hi] of block <b>
 * -- Forward block from fseed[b], Backward block from bseed[b], then
 * p7b_decode_row() per row -- and point <view> at them.
 *
 * The per-row inputs handed to p7b_decode_row() are assembled exactly
 * as p7b_backdecode() assembles them (same seed-vs-fx seam at r==lo,
 * same gap[] adjustment of the Forward specials), so the rows produced
 * here are bit-identical to the resident ones.
 */
static void
p7b_ppwork_fill(P7B_PPWORK *w, const P7_PROFILE *gm, const P7B_GEO *g,
                int b, int lo, int hi, P7B_PPVIEW *view)
{
  float   s_xN, s_xJ, s_xC;
  int     r;
  int64_t cum = 0;

  p7b_fwd_block(w->dsq, gm, g, w->fseed, b, lo, hi, w->fblk, w->froff,
                w->fxN, w->fxJ, w->fxC, &s_xN, &s_xJ, &s_xC);
  p7b_bwd_block(w->dsq, gm, g, w->bseed, b, lo, hi, w->bblk, w->broff,
                w->bxN, w->bxJ, w->bxC);

  for (r = lo; r <= hi; r++)
    {
      int   kac = g->ka[r], kbc = g->kb[r];
      float fwd_xN_prev, fwd_xJ_prev, fwd_xC_prev;

      w->pproff[r-lo] = cum;

      if (r == lo) { fwd_xN_prev = s_xN; fwd_xJ_prev = s_xJ; fwd_xC_prev = s_xC; }
      else         { fwd_xN_prev = w->fxN[r-1-lo]; fwd_xJ_prev = w->fxJ[r-1-lo]; fwd_xC_prev = w->fxC[r-1-lo]; }
      if (g->gap[r] > 0) {
        fwd_xN_prev += g->gap[r] * gm->xsc[p7P_N][p7P_LOOP];
        fwd_xJ_prev += g->gap[r] * gm->xsc[p7P_J][p7P_LOOP];
        fwd_xC_prev += g->gap[r] * gm->xsc[p7P_C][p7P_LOOP];
      }

      p7b_decode_row(g, gm, w->fwdsc, kac, kbc,
                     w->fblk + w->froff[r-lo], w->bblk + w->broff[r-lo],
                     w->bxN[r-lo], w->bxJ[r-lo], w->bxC[r-lo],
                     fwd_xN_prev, fwd_xJ_prev, fwd_xC_prev,
                     w->ppblk + cum, w->ppxmx + (int64_t)(r-lo) * p7G_NXCELLS);

      cum += (int64_t)(kbc - kac + 1) * P7B_PP_NSCELLS;
    }

  view->pp     = NULL;
  view->ppblk  = w->ppblk;
  view->pproff = w->pproff;
  view->ppxmx  = w->ppxmx;
  view->lo     = lo;
  view->hi     = hi;
}


/* Function: p7_GCheckptFBDecode_Banded()
 * Incept:   Brief 26_0526-016
 *
 * Purpose:  sqrt(nrow)-row checkpointed banded Forward + Backward +
 *           posterior decoding. Drop-in replacement for the
 *               my_p7_GForwardBanded + p7_GBackwardBanded + p7_GDecodingBanded
 *           sequence, producing a byte-identical posterior matrix <pp>
 *           and Forward score, but never materializing either full F or
 *           B banded matrix (only the resident <pp> and O(sqrt(nrow))
 *           working set).
 *
 * Args:     dsq       - digital sequence 1..L
 *           L         - sequence length
 *           gm        - profile (glocal)
 *           pp        - RESULT: posterior matrix (created by caller from
 *                       the same P7_GBANDS); filled here.
 *           ret_fwdsc - RETURN: Forward score in nats
 *
 * Returns:  eslOK on success.
 */
int
p7_GCheckptFBDecode_Banded(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm,
                           P7_GMXB *pp, float *ret_fwdsc)
{
  int        status;
  P7B_GEO   *g = NULL;
  P7B_FSEED *seed = NULL;
  float      fwdsc = 0.;
  int        b;

  g = p7b_geo_Create(pp->bnd, gm->M, L);
  if (g == NULL) { status = eslEMEM; goto ERROR; }

  ESL_ALLOC(seed, sizeof(P7B_FSEED) * (g->nblk > 0 ? g->nblk : 1));
  for (b = 0; b < g->nblk; b++) { seed[b].dp = NULL; seed[b].valid = 0; }

  if ((status = p7b_forward_seeds(dsq, gm, g, seed, &fwdsc)) != eslOK) goto ERROR;
  if ((status = p7b_backdecode  (dsq, gm, g, seed, fwdsc, pp)) != eslOK) goto ERROR;

  if (getenv("INFERNAL_CKPT_VERBOSE"))
    fprintf(stderr, "# p7_GCheckptFBDecode_Banded: M=%d L=%d nrow=%d B=%d nblk=%d fwdsc=%.6f\n",
            gm->M, L, g->nrow, g->B, g->nblk, fwdsc);

  for (b = 0; b < g->nblk; b++) if (seed[b].dp) free(seed[b].dp);
  free(seed);
  p7b_geo_Destroy(g);
  if (ret_fwdsc) *ret_fwdsc = fwdsc;
  return eslOK;

 ERROR:
  if (seed) { for (b = 0; b < (g ? g->nblk : 0); b++) if (seed[b].dp) free(seed[b].dp); free(seed); }
  if (g) p7b_geo_Destroy(g);
  return status;
}


/*****************************************************************
 * STEP OA + TB: sqrt(nrow)-checkpointed banded optimal accuracy
 *****************************************************************/

#define P7B_TSCDELTA(s,k) ( (tsc[(k) * p7P_NTRANS + (s)] == -eslINFINITY) ? FLT_MIN : 1.0f )

/* p7b_oa_row(): fill one banded OA row.  Verbatim copy of the per-row
 * body of p7_GOptimalAccuracyBanded().  Reads the resident pp row
 * (pp_dp cells, pp_xp specials); writes OA dp (ascending k) + returns
 * specials. */
static void
p7b_oa_row(const P7_PROFILE *gm, int M, float esc,
           int kac, int kbc, const float *dpp, int kap, int kbp,
           const float *pp_dp, const float *pp_xp,
           float xB, float xN_in, float xJ_in, float xC_in,
           float *dpc,
           float *o_xE, float *o_xN, float *o_xJ, float *o_xB, float *o_xC)
{
  float const *tsc = gm->tsc;
  float  mvp, ivp, dvp, dc, sc, xE, pp_m, t1, t2;
  int    k;
  int    kbc2 = (kbc == M ? kbc - 1 : kbc);
  const float *ppp = pp_dp;

  dc = -eslINFINITY;
  xE = -eslINFINITY;

  if (dpp != NULL) dpp += (kac-1 > kap ? ESL_MIN(kac-kap-1, kbp-kap+1) * p7G_NSCELLS : 0);
  if (dpp != NULL && kac > kap && kac-1 <= kbp) { mvp = *dpp++; ivp = *dpp++; dvp = *dpp++; }
  else                                          { mvp = -eslINFINITY; ivp = -eslINFINITY; dvp = -eslINFINITY; }

  for (k = kac; k <= kbc2; k++)
    {
      pp_m = ppp[0];
      *dpc++ = sc = ESL_MAX(ESL_MAX(P7B_TSCDELTA(p7P_MM, k-1) * (mvp + pp_m),
                                    P7B_TSCDELTA(p7P_IM, k-1) * (ivp + pp_m)),
                            ESL_MAX(P7B_TSCDELTA(p7P_DM, k-1) * (dvp + pp_m),
                                    P7B_TSCDELTA(p7P_BM, k-1) * (xB  + pp_m)));

      if (dpp != NULL && k >= kap && k <= kbp) { mvp = *dpp++; ivp = *dpp++; dvp = *dpp++; }
      else                                     { mvp = -eslINFINITY; ivp = -eslINFINITY; dvp = -eslINFINITY; }

      *dpc++ = ESL_MAX(P7B_TSCDELTA(p7P_MI, k) * (mvp + ppp[1]),
                       P7B_TSCDELTA(p7P_II, k) * (ivp + ppp[1]));

      xE = ESL_MAX(xE, esc * sc);

      *dpc++ = dc;
      dc = ESL_MAX(P7B_TSCDELTA(p7P_MD, k) * sc, P7B_TSCDELTA(p7P_DD, k) * dc);

      ppp += P7B_PP_NSCELLS;   /* brief 26_0526-017: pp is compact 2-cell (M,I) */
    }

  if (kbc2 < kbc) /* kbc==M unrolled */
    {
      pp_m = ppp[0];
      *dpc++ = sc = ESL_MAX(ESL_MAX(P7B_TSCDELTA(p7P_MM, k-1) * (mvp + pp_m),
                                    P7B_TSCDELTA(p7P_IM, k-1) * (ivp + pp_m)),
                            ESL_MAX(P7B_TSCDELTA(p7P_DM, k-1) * (dvp + pp_m),
                                    P7B_TSCDELTA(p7P_BM, k-1) * (xB  + pp_m)));
      *dpc++ = -eslINFINITY;
      *dpc++ = dc;
      xE = ESL_MAX(xE, ESL_MAX(sc, dc));
      ppp += P7B_PP_NSCELLS;   /* brief 26_0526-017: pp is compact 2-cell (M,I) */
    }

  *o_xE = xE;
  t1 = (gm->xsc[p7P_N][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0f;
  *o_xN = t1 * (xN_in + pp_xp[p7G_N]);
  t1 = (gm->xsc[p7P_J][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0f;
  t2 = (gm->xsc[p7P_E][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0f;
  *o_xJ = ESL_MAX(t1 * (xJ_in + pp_xp[p7G_J]), t2 * xE);
  t1 = (gm->xsc[p7P_N][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0f;
  t2 = (gm->xsc[p7P_J][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0f;
  *o_xB = ESL_MAX(t1 * (*o_xN), t2 * (*o_xJ));
  t1 = (gm->xsc[p7P_C][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0f;
  t2 = (gm->xsc[p7P_E][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0f;
  *o_xC = ESL_MAX(t1 * (xC_in + pp_xp[p7G_C]), t2 * xE);
}

/* p7b_oa_enter(): compute the entering OA state for a row, factoring the
 * segment-start reset (OA resets N/J/C to -inf at a segment start; it
 * does NOT advance them through the gap the way Forward does). */
static void
p7b_oa_seg_start(const P7_PROFILE *gm, float *xN, float *xJ, float *xC, float *xB, int do_reset)
{
  float t1, t2;
  if (do_reset) { *xN = -eslINFINITY; *xJ = -eslINFINITY; *xC = -eslINFINITY; }
  t1 = (gm->xsc[p7P_N][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0f;
  t2 = (gm->xsc[p7P_J][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0f;
  *xB = ESL_MAX(t1 * (*xN), t2 * (*xJ));
}

/* p7b_oa_seeds(): STEP OA.  Stream banded OA forward once, store the
 * entering state per block, return the OA score.
 *
 * brief 26_0628-081: the posterior is read through <view>.  When <ppw>
 * is NULL, <view> is the legacy resident pp and nothing else changes.
 * When <ppw> is non-NULL the resident pp does not exist: at each block
 * boundary the block's pp rows are recomputed into <ppw>'s O(B*maxnc)
 * buffer and <view> is repointed at them.  The row loop, its order, and
 * the per-row arithmetic are untouched either way. */
static int
p7b_oa_seeds(const P7_PROFILE *gm, const P7B_GEO *g, P7B_PPVIEW *view,
             P7B_PPWORK *ppw, P7B_FSEED *seed, float *ret_oasc)
{
  int    status;
  int    M = g->M, B = g->B;
  float  esc = p7_profile_IsLocal(gm) ? 1.0f : 0.0f;
  float *buf0 = NULL, *buf1 = NULL, *prev = NULL, *cur;
  int    r, kap = M+1, kbp = M+1;
  float  xN = 0.0f, xJ = -eslINFINITY, xC = -eslINFINITY, xB = -eslINFINITY;
  float  xE_o, xN_o, xJ_o, xB_o, xC_o;

  ESL_ALLOC(buf0, sizeof(float) * (g->maxnc * p7G_NSCELLS + 1));
  ESL_ALLOC(buf1, sizeof(float) * (g->maxnc * p7G_NSCELLS + 1));

  for (r = 0; r < g->nrow; r++)
    {
      int   kac = g->ka[r], kbc = g->kb[r];
      const float *dpp;
      int   seg_start = (g->gap[r] > 0) || (r == 0);
      const float *pp_dp, *pp_xp;

      if (ppw != NULL && (r % B) == 0)   /* r starts block r/B: refill the pp window */
        p7b_ppwork_fill(ppw, gm, g, r / B, r, ESL_MIN(r + B - 1, g->nrow - 1), view);

      pp_dp = p7b_ppview_dp(view, g, r);
      pp_xp = p7b_ppview_xp(view, g, r);

      cur = (prev == buf0) ? buf1 : buf0;

      if (seg_start) { p7b_oa_seg_start(gm, &xN, &xJ, &xC, &xB, (g->gap[r] > 0)); dpp = NULL; }
      else           { dpp = prev; }

      p7b_oa_row(gm, M, esc, kac, kbc, dpp, kap, kbp, pp_dp, pp_xp,
                 xB, xN, xJ, xC, cur, &xE_o, &xN_o, &xJ_o, &xB_o, &xC_o);

      xN = xN_o; xJ = xJ_o; xC = xC_o; xB = xB_o;
      kap = kac; kbp = kbc;
      prev = cur;

      if (r + 1 < g->nrow && ((r + 1) % B) == 0) {
        int b = (r + 1) / B, nc = kbc - kac + 1;
        ESL_ALLOC(seed[b].dp, sizeof(float) * nc * p7G_NSCELLS);
        memcpy(seed[b].dp, cur, sizeof(float) * nc * p7G_NSCELLS);
        seed[b].ka = kac; seed[b].kb = kbc;
        seed[b].xN = xN; seed[b].xJ = xJ; seed[b].xC = xC; seed[b].xB = xB;
        seed[b].valid = 1;
      }
    }

  *ret_oasc = xC;   /* mirror p7_GOptimalAccuracyBanded: ret_e = xC after last row */
  free(buf0); free(buf1);
  return eslOK;

 ERROR:
  if (buf0) free(buf0);
  if (buf1) free(buf1);
  return status;
}

/* p7b_oa_fill_region(): recompute OA rows [lo..hi] from seed[b] into
 * caller buffers, producing BOTH dp (oadp, packed, offsets in oaoff)
 * and xmx (oaxmx, p7G_NXCELLS per row).  Used by the traceback window. */
static void
p7b_oa_fill_region(const P7_PROFILE *gm, const P7B_GEO *g, const P7B_PPVIEW *view,
                   const P7B_FSEED *seed, int b, int lo, int hi,
                   float *oadp, int64_t *oaoff, float *oaxmx)
{
  int    M = g->M;
  float  esc = p7_profile_IsLocal(gm) ? 1.0f : 0.0f;
  int    r, kap, kbp;
  float  xN, xJ, xC, xB = -eslINFINITY;
  int64_t cum = 0;

  if (b == 0) { xN = 0.0f; xJ = -eslINFINITY; xC = -eslINFINITY; kap = kbp = M+1; }
  else        { xN = seed[b].xN; xJ = seed[b].xJ; xC = seed[b].xC; xB = seed[b].xB;
                kap = seed[b].ka; kbp = seed[b].kb; }

  for (r = lo; r <= hi; r++)
    {
      int   kac = g->ka[r], kbc = g->kb[r];
      float xE_o, xN_o, xJ_o, xB_o, xC_o;
      float *dpc = oadp + cum;
      float *xpc = oaxmx + (int64_t)(r - lo) * p7G_NXCELLS;
      const float *dpp;
      int   seg_start = (g->gap[r] > 0) || (r == 0);
      const float *pp_dp = p7b_ppview_dp(view, g, r);
      const float *pp_xp = p7b_ppview_xp(view, g, r);

      if (seg_start) { p7b_oa_seg_start(gm, &xN, &xJ, &xC, &xB, (g->gap[r] > 0)); dpp = NULL; }
      else           { dpp = (r == lo) ? seed[b].dp : (oadp + oaoff[r-1-lo]); }

      p7b_oa_row(gm, M, esc, kac, kbc, dpp, kap, kbp, pp_dp, pp_xp,
                 xB, xN, xJ, xC, dpc, &xE_o, &xN_o, &xJ_o, &xB_o, &xC_o);

      xpc[p7G_E] = xE_o; xpc[p7G_N] = xN_o; xpc[p7G_J] = xJ_o; xpc[p7G_B] = xB_o; xpc[p7G_C] = xC_o;
      oaoff[r-lo] = cum;
      xN = xN_o; xJ = xJ_o; xC = xC_o; xB = xB_o;
      kap = kac; kbp = kbc;
      cum += (int64_t)(kbc - kac + 1) * p7G_NSCELLS;
    }
}

/* On-demand recomputed OA window for the checkpointed traceback.
 * brief 26_0628-081: when <ppw> is non-NULL the same window also holds
 * the block's recomputed posterior rows (in <view>), so the traceback's
 * PP_* accessors no longer need a resident pp either.  The window is
 * B+1 rows -- one full block plus the +1 overlap row -- which is what
 * lets the (i-1, i) row pairs in the C/J cases resolve without
 * thrashing (see the accessor-ordering comment in p7b_oa_trace). */
typedef struct {
  const P7B_GEO    *g;
  const int        *row_idx;
  const P7_PROFILE *gm;
  P7B_PPVIEW        view;
  P7B_PPWORK       *ppw;      /* NULL => view.pp is the resident pp */
  const P7B_FSEED  *seed;
  int     B, nrow;
  float  *oadp;  int64_t *oaoff;  float *oaxmx;
  int     wlo, whi;
} P7B_OAWIN;

static void
p7b_oawin_ensure(P7B_OAWIN *w, int r)
{
  if (r < w->wlo || r > w->whi) {
    int b = r / w->B;
    w->wlo = b * w->B;
    w->whi = ESL_MIN((b+1) * w->B, w->nrow - 1);
    if (w->ppw != NULL)   /* recompute this block's posterior rows first: OA reads them */
      p7b_ppwork_fill(w->ppw, w->gm, w->g, b, w->wlo, w->whi, &w->view);
    p7b_oa_fill_region(w->gm, w->g, &w->view, w->seed, b, w->wlo, w->whi, w->oadp, w->oaoff, w->oaxmx);
  }
}

static float
p7b_oawin_dp(P7B_OAWIN *w, int p, int k, int cell)
{
  int r;
  if (p < 0 || p > w->g->L) return -eslINFINITY;
  r = w->row_idx[p];
  if (r < 0) return -eslINFINITY;
  if (k < w->g->ka[r] || k > w->g->kb[r]) return -eslINFINITY;
  p7b_oawin_ensure(w, r);
  return w->oadp[ w->oaoff[r - w->wlo] + (int64_t)(k - w->g->ka[r]) * p7G_NSCELLS + cell ];
}

static float
p7b_oawin_x(P7B_OAWIN *w, int p, int s)
{
  int r;
  if (p == 0) return ((s == p7G_N || s == p7G_B) ? 0.0f : -eslINFINITY);
  if (p < 0 || p > w->g->L) return -eslINFINITY;
  r = w->row_idx[p];
  if (r < 0) return -eslINFINITY;
  p7b_oawin_ensure(w, r);
  return w->oaxmx[ (int64_t)(r - w->wlo) * p7G_NXCELLS + s ];
}

/* brief 26_0628-081: pp is read through the same window as OA.  With a
 * resident pp (w->ppw == NULL) no ensure() is needed and none is done,
 * so the legacy path's recompute count is unchanged. */
static float
p7b_pp_dp(P7B_OAWIN *w, int p, int k, int cell)
{
  const P7B_GEO *g = w->g;
  int r;
  if (p < 0 || p > g->L) return -eslINFINITY;
  r = w->row_idx[p];
  if (r < 0) return -eslINFINITY;
  if (k < g->ka[r] || k > g->kb[r]) return -eslINFINITY;
  if (w->ppw != NULL) p7b_oawin_ensure(w, r);
  /* brief 26_0526-017: 2-cell pp; cell in {p7G_M=0, p7G_I=1} */
  return p7b_ppview_dp(&w->view, g, r)[ (int64_t)(k - g->ka[r]) * P7B_PP_NSCELLS + cell ];
}

static float
p7b_pp_x(P7B_OAWIN *w, int p, int s)
{
  const P7B_GEO *g = w->g;
  int r;
  if (p == 0) return ((s == p7G_N || s == p7G_B) ? 0.0f : -eslINFINITY);
  if (p < 0 || p > g->L) return -eslINFINITY;
  r = w->row_idx[p];
  if (r < 0) return -eslINFINITY;
  if (w->ppw != NULL) p7b_oawin_ensure(w, r);
  return p7b_ppview_xp(&w->view, g, r)[s];
}

/* p7b_oa_trace(): STEP TB.  Byte-exact reproduction of p7_GOATraceBanded,
 * but the OA matrix is accessed through an on-demand recomputed window
 * (one block + 1 overlap row) instead of a full resident matrix.
 *
 * brief 26_0628-081: when <ppw> is non-NULL, the posterior is windowed
 * the same way (recomputed per block from the F/B seeds) instead of
 * being fully resident; when <ppw> is NULL, <view> carries the legacy
 * resident pp and behaviour is unchanged. */
static int
p7b_oa_trace(const P7_PROFILE *gm, const P7B_GEO *g, const P7B_PPVIEW *view,
             P7B_PPWORK *ppw, const P7B_FSEED *seed, P7_TRACE *tr)
{
  float const *tsc = gm->tsc;
  int     M = g->M, L = g->L, B = g->B, nrow = g->nrow;
  int     status;

  /* pp lookup: resident (legacy) or windowed (ckpt-pp) */
  int    *row_idx = NULL;   /* seq pos -> banded row r (-1 outside) */
  int     i, k, r;

  /* OA window state */
  float  *oadp  = NULL;     /* B+1 rows of OA dp, packed       */
  int64_t*oaoff = NULL;     /* [0..B] per-window-row dp offset */
  float  *oaxmx = NULL;     /* [0..B] * p7G_NXCELLS            */
  int     win_cap = B + 2;
  P7B_OAWIN w;

  int     sprv, scur;
  float   postprob;

  ESL_ALLOC(row_idx, sizeof(int) * (L + 1));
  for (i = 0; i <= L; i++) row_idx[i] = -1;
  for (r = 0; r < nrow; r++) row_idx[g->irow[r]] = r;

  ESL_ALLOC(oadp,  sizeof(float)   * ((int64_t) win_cap * g->maxnc * p7G_NSCELLS + 1));
  ESL_ALLOC(oaoff, sizeof(int64_t) * win_cap);
  ESL_ALLOC(oaxmx, sizeof(float)   * ((int64_t) win_cap * p7G_NXCELLS));

  w.g = g; w.row_idx = row_idx; w.gm = gm; w.view = *view; w.ppw = ppw; w.seed = seed;
  w.B = B; w.nrow = nrow; w.oadp = oadp; w.oaoff = oaoff; w.oaxmx = oaxmx;
  w.wlo = -1; w.whi = -2;

  /* OA accessors (windowed, recompute on demand); PP accessors (resident).
   * In the C/J cases that read both row i and row i-1, the LOWER row
   * (i-1) is evaluated first so the +1 overlap row keeps row i in the
   * window. */
#define OA_M(p,kk) p7b_oawin_dp(&w, (p), (kk), p7G_M)
#define OA_I(p,kk) p7b_oawin_dp(&w, (p), (kk), p7G_I)
#define OA_D(p,kk) p7b_oawin_dp(&w, (p), (kk), p7G_D)
#define OA_X(p,ss) p7b_oawin_x(&w, (p), (ss))
#define PP_M(p,kk) p7b_pp_dp(&w, (p), (kk), p7G_M)
#define PP_I(p,kk) p7b_pp_dp(&w, (p), (kk), p7G_I)
#define PP_X(p,ss) p7b_pp_x(&w, (p), (ss))

  i = L; k = 0;
  if ((status = p7_trace_AppendWithPP(tr, p7T_T, k, i, 0.0)) != eslOK) goto ERROR;
  if ((status = p7_trace_AppendWithPP(tr, p7T_C, k, i, 0.0)) != eslOK) goto ERROR;

  sprv = p7T_C;
  while (sprv != p7T_S)
    {
      switch (sprv) {
      case p7T_M:
        {
          if (k - 1 == 0) {
            /* Node 0 has no M/D state: the only valid entry into M1 is B->M1. */
            scur = p7T_B;
          } else {
            float path[4];
            int   state[4] = { p7T_M, p7T_I, p7T_D, p7T_B };
            path[0] = P7B_TSCDELTA(p7P_MM, k-1) * OA_M(i-1, k-1);
            path[1] = P7B_TSCDELTA(p7P_IM, k-1) * OA_I(i-1, k-1);
            path[2] = P7B_TSCDELTA(p7P_DM, k-1) * OA_D(i-1, k-1);
            path[3] = P7B_TSCDELTA(p7P_BM, k-1) * OA_X(i-1, p7G_B);
            scur = state[esl_vec_FArgMax(path, 4)];
          }
          k--; i--;
        }
        break;
      case p7T_D:
        {
          float path[2];
          path[0] = P7B_TSCDELTA(p7P_MD, k-1) * OA_M(i, k-1);
          path[1] = P7B_TSCDELTA(p7P_DD, k-1) * OA_D(i, k-1);
          scur = (path[0] >= path[1]) ? p7T_M : p7T_D;
          k--;
        }
        break;
      case p7T_I:
        {
          float path[2];
          path[0] = P7B_TSCDELTA(p7P_MI, k) * OA_M(i-1, k);
          path[1] = P7B_TSCDELTA(p7P_II, k) * OA_I(i-1, k);
          scur = (path[0] >= path[1]) ? p7T_M : p7T_I;
          i--;
        }
        break;
      case p7T_N:
        scur = (i == 0) ? p7T_S : p7T_N;
        break;
      case p7T_C:
        {
          float t1c = (gm->xsc[p7P_C][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0f;
          float t2c = (gm->xsc[p7P_E][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0f;
          float path[2];
          path[0] = t1c * (OA_X(i-1, p7G_C) + PP_X(i, p7G_C));
          path[1] = t2c *  OA_X(i, p7G_E);
          scur = (path[0] > path[1]) ? p7T_C : p7T_E;
        }
        break;
      case p7T_J:
        {
          float t1j = (gm->xsc[p7P_J][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0f;
          float t2j = (gm->xsc[p7P_E][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0f;
          float path[2];
          path[0] = t1j * (OA_X(i-1, p7G_J) + PP_X(i, p7G_J));
          path[1] = t2j *  OA_X(i, p7G_E);
          scur = (path[0] > path[1]) ? p7T_J : p7T_E;
        }
        break;
      case p7T_E:
        {
          float max  = -eslINFINITY;
          int   smax = -1, kmax = -1;
          int   ri   = row_idx[i];
          if (! p7_profile_IsLocal(gm)) {
            k = M;
            scur = (OA_M(i, M) >= OA_D(i, M)) ? p7T_M : p7T_D;
          } else {
            if (ri >= 0) {
              int ka = g->ka[ri], kb = g->kb[ri];
              for (k = ka; k <= kb; k++) {
                float mv = OA_M(i, k);
                float dv = OA_D(i, k);
                if (mv >= max) { max = mv; smax = p7T_M; kmax = k; }
                if (dv >  max) { max = dv; smax = p7T_D; kmax = k; }
              }
            }
            k = kmax; scur = smax;
          }
        }
        break;
      case p7T_B:
        {
          float t1b = (gm->xsc[p7P_N][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0f;
          float t2b = (gm->xsc[p7P_J][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0f;
          float path[2];
          path[0] = t1b * OA_X(i, p7G_N);
          path[1] = t2b * OA_X(i, p7G_J);
          if (gm->xsc[p7P_E][p7P_LOOP] == -eslINFINITY) scur = p7T_N;
          else scur = (path[0] > path[1]) ? p7T_N : p7T_J;
        }
        break;
      default:
        status = eslEINVAL; goto ERROR;
      }

      if (scur == -1) { status = eslEINVAL; goto ERROR; }

      switch (scur) {
      case p7T_M: postprob = PP_M(i, k);                                 break;
      case p7T_I: postprob = PP_I(i, k);                                 break;
      case p7T_N: postprob = (scur == sprv) ? PP_X(i, p7G_N) : 0.0f;     break;
      case p7T_C: postprob = (scur == sprv) ? PP_X(i, p7G_C) : 0.0f;     break;
      case p7T_J: postprob = (scur == sprv) ? PP_X(i, p7G_J) : 0.0f;     break;
      default:    postprob = 0.0f;                                       break;
      }

      if ((status = p7_trace_AppendWithPP(tr, scur, k, i, postprob)) != eslOK) goto ERROR;

      if ((scur == p7T_N || scur == p7T_J || scur == p7T_C) && scur == sprv) i--;
      sprv = scur;
    }

  tr->M = M; tr->L = L;
  status = p7_trace_Reverse(tr);

  free(row_idx); free(oadp); free(oaoff); free(oaxmx);
  return status;

#undef OA_M
#undef OA_I
#undef OA_D
#undef OA_X
#undef PP_M
#undef PP_I
#undef PP_X
 ERROR:
  if (row_idx) free(row_idx);
  if (oadp)    free(oadp);
  if (oaoff)   free(oaoff);
  if (oaxmx)   free(oaxmx);
  return status;
}

#undef P7B_TSCDELTA

/* Function: p7_GCheckptOA_Banded()
 * Incept:   Brief 26_0526-016
 *
 * Purpose:  sqrt(nrow)-row checkpointed banded optimal-accuracy alignment.
 *           Drop-in replacement for the
 *               p7_GOptimalAccuracyBanded + p7_GOATraceBanded
 *           pair, producing a byte-identical OA trace (alignment + PP)
 *           and OA score, but storing only ~sqrt(nrow) OA seed rows
 *           instead of a full OA banded matrix; the traceback recomputes
 *           OA row blocks from the seeds on demand.  Reads the resident
 *           posterior <pp> (from p7_GCheckptFBDecode_Banded).
 *
 * Args:     gm       - profile
 *           pp       - banded posterior matrix (resident)
 *           tr       - RESULT: OA trace (caller-provided, Reuse'd)
 *           ret_oasc - RETURN: OA score
 *
 * Returns:  eslOK on success.
 */
int
p7_GCheckptOA_Banded(const P7_PROFILE *gm, P7_GMXB *pp, P7_TRACE *tr, float *ret_oasc)
{
  int         status;
  P7B_GEO    *g = NULL;
  P7B_FSEED  *seed = NULL;
  P7B_PPVIEW  view;
  float       oasc = 0.;
  int         b;

  g = p7b_geo_Create(pp->bnd, gm->M, pp->bnd->L);
  if (g == NULL) { status = eslEMEM; goto ERROR; }

  ESL_ALLOC(seed, sizeof(P7B_FSEED) * (g->nblk > 0 ? g->nblk : 1));
  for (b = 0; b < g->nblk; b++) { seed[b].dp = NULL; seed[b].valid = 0; }

  /* resident-pp view (brief 26_0628-081): ppw==NULL selects the legacy path */
  view.pp = pp; view.ppblk = NULL; view.pproff = NULL; view.ppxmx = NULL;
  view.lo = 0;  view.hi = g->nrow - 1;

  if ((status = p7b_oa_seeds(gm, g, &view, NULL, seed, &oasc)) != eslOK) goto ERROR;
  if ((status = p7b_oa_trace(gm, g, &view, NULL, seed, tr))    != eslOK) goto ERROR;

  for (b = 0; b < g->nblk; b++) if (seed[b].dp) free(seed[b].dp);
  free(seed);
  p7b_geo_Destroy(g);
  if (ret_oasc) *ret_oasc = oasc;
  return eslOK;

 ERROR:
  if (seed) { for (b = 0; b < (g ? g->nblk : 0); b++) if (seed[b].dp) free(seed[b].dp); free(seed); }
  if (g) p7b_geo_Destroy(g);
  return status;
}

/* Function: p7_GCheckptFBDecodeOA_Banded()
 * Incept:   brief 26_0628-081
 *
 * Purpose:  Fully checkpointed banded Forward/Backward/Decode/optimal-
 *           accuracy/traceback -- the DOUBLE-checkpointed engine.  Drop-in
 *           replacement for the
 *               p7b_pp_Create + p7_GCheckptFBDecode_Banded + p7_GCheckptOA_Banded
 *           sequence, producing a bit-identical Forward score, OA score
 *           and OA trace (alignment + posterior annotation), but never
 *           allocating the O(bnd->ncell) resident posterior that that
 *           sequence needs.
 *
 *           Passes, in order:
 *             A  p7b_forward_seeds()  -- stream Forward,  lay ~sqrt(nrow) F seeds
 *             B1 p7b_backward_seeds() -- stream Backward, lay ~sqrt(nrow) B seeds
 *             OA p7b_oa_seeds()       -- per block: recompute F and B rows from
 *                                        their seeds, form that block's pp rows,
 *                                        run OA over them ascending, lay OA seeds
 *             TB p7b_oa_trace()       -- traceback over an on-demand recomputed
 *                                        (pp + OA) block window
 *
 *           Peak resident DP term: O(sqrt(nrow) * maxnc) -- three seed
 *           arrays plus a handful of one-block buffers.  Cost: two extra
 *           O(ncell) recompute passes relative to the resident-pp engine.
 *
 * Args:     dsq       - digital sequence 1..L
 *           L         - sequence length
 *           gm        - profile
 *           bnd       - the band (caller keeps ownership)
 *           tr        - RESULT: OA trace (caller-provided; Reuse'd by caller)
 *           ret_fwdsc - RETURN: Forward score in nats
 *           ret_oasc  - RETURN: OA score
 *
 * Returns:  eslOK on success.
 */
int
p7_GCheckptFBDecodeOA_Banded(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm,
                             P7_GBANDS *bnd, P7_TRACE *tr,
                             float *ret_fwdsc, float *ret_oasc)
{
  int         status;
  P7B_GEO    *g      = NULL;
  P7B_FSEED  *fseed  = NULL;
  P7B_BSEED  *bseed  = NULL;
  P7B_FSEED  *oaseed = NULL;
  P7B_PPWORK *ppw    = NULL;
  P7B_PPVIEW  view;
  float       fwdsc = 0., oasc = 0.;
  int         b;

  g = p7b_geo_Create(bnd, gm->M, L);
  if (g == NULL) { status = eslEMEM; goto ERROR; }

  ESL_ALLOC(fseed,  sizeof(P7B_FSEED) * (g->nblk > 0 ? g->nblk : 1));
  ESL_ALLOC(bseed,  sizeof(P7B_BSEED) * (g->nblk > 0 ? g->nblk : 1));
  ESL_ALLOC(oaseed, sizeof(P7B_FSEED) * (g->nblk > 0 ? g->nblk : 1));
  for (b = 0; b < g->nblk; b++) { fseed[b].dp  = NULL; fseed[b].valid  = 0; }
  for (b = 0; b < g->nblk; b++) { bseed[b].dp  = NULL; bseed[b].valid  = 0; }
  for (b = 0; b < g->nblk; b++) { oaseed[b].dp = NULL; oaseed[b].valid = 0; }

  if ((status = p7b_forward_seeds (dsq, gm, g, fseed, &fwdsc)) != eslOK) goto ERROR;
  if ((status = p7b_backward_seeds(dsq, gm, g, bseed))         != eslOK) goto ERROR;

  ppw = p7b_ppwork_Create(dsq, g, fseed, bseed, fwdsc);
  if (ppw == NULL) { status = eslEMEM; goto ERROR; }

  view.pp = NULL; view.ppblk = NULL; view.pproff = NULL; view.ppxmx = NULL;
  view.lo = 0; view.hi = -1;

  if ((status = p7b_oa_seeds(gm, g, &view, ppw, oaseed, &oasc))  != eslOK) goto ERROR;
  if ((status = p7b_oa_trace(gm, g, &view, ppw, oaseed, tr))     != eslOK) goto ERROR;

  if (getenv("INFERNAL_CKPT_VERBOSE"))
    fprintf(stderr, "# p7_GCheckptFBDecodeOA_Banded: M=%d L=%d nrow=%d B=%d nblk=%d maxnc=%d fwdsc=%.6f oasc=%.6f\n",
            gm->M, L, g->nrow, g->B, g->nblk, g->maxnc, fwdsc, oasc);

  for (b = 0; b < g->nblk; b++) { if (fseed[b].dp)  free(fseed[b].dp);
                                  if (bseed[b].dp)  free(bseed[b].dp);
                                  if (oaseed[b].dp) free(oaseed[b].dp); }
  free(fseed); free(bseed); free(oaseed);
  p7b_ppwork_Destroy(ppw);
  p7b_geo_Destroy(g);
  if (ret_fwdsc) *ret_fwdsc = fwdsc;
  if (ret_oasc)  *ret_oasc  = oasc;
  return eslOK;

 ERROR:
  for (b = 0; b < (g ? g->nblk : 0); b++) {
    if (fseed  && fseed[b].dp)  free(fseed[b].dp);
    if (bseed  && bseed[b].dp)  free(bseed[b].dp);
    if (oaseed && oaseed[b].dp) free(oaseed[b].dp);
  }
  if (fseed)  free(fseed);
  if (bseed)  free(bseed);
  if (oaseed) free(oaseed);
  if (ppw)    p7b_ppwork_Destroy(ppw);
  if (g)      p7b_geo_Destroy(g);
  return status;
}

/* Function: cm_nodepad_cmpint()
 * Helper for qsort in cm_ComputeP7CMNodePad().
 */
static int cm_nodepad_cmpint(const void *a, const void *b) {
  int x = *(const int *)a, y = *(const int *)b;
  return (x<y)?-1:(x>y);
}

/* Structure for work units passed through the work queue
 * in the threaded cm_ComputeP7CMNodePad() path.
 * A small recycling pool of these structs cycles between reader and workers.
 *
 * The master thread performs all RNG-consuming steps (EmitParsetree and
 * flanking-residue generation) in sample order, then hands the emitted
 * parsetree + sequence + flanked digitized sequence to a worker via this
 * struct. The worker does pure Viterbi DP and deficit recording, then
 * frees the cm_tr/esq/emb allocations (ownership transfers on pickup).
 * This makes results deterministic regardless of thread scheduling, and
 * bit-for-bit identical to the serial path.
 */
typedef struct {
  int           idx;    /* sample index to process; -1 = sentinel (stop signal) */
  Parsetree_t  *cm_tr;  /* emitted CM parsetree (filled by master; freed by worker) */
  ESL_SQ       *esq;    /* emitted sequence     (filled by master; freed by worker) */
  ESL_DSQ      *emb;    /* flanked digitized sequence (filled by master; freed by worker) */
  int           L;      /* emitted sequence length */
  int           L_emb;  /* total length of emb = 2*flank + L */
} CM_NODEPAD_WORK;

#ifdef HMMER_THREADS
/* Per-worker info for threaded cm_ComputeP7CMNodePad().
 * Each worker maintains its own per-node deficit arrays and p7 objects.
 * The main thread merges the deficit arrays after all workers finish.
 */
typedef struct {
  CM_t           *cm;         /* shared CM (read-only) */
  ESL_WORK_QUEUE *queue;      /* shared work queue */
  P7_BG          *bg;         /* per-worker p7 background model */
  P7_PROFILE     *gm;         /* per-worker p7 profile */
  P7_GMX         *gx;         /* per-worker p7 DP matrix */
  P7_TRACE       *p7tr;       /* per-worker p7 trace */
  int           **deficits;   /* deficits[k][0..def_n[k]-1]: per-node deficit values */
  int            *def_n;      /* def_n[k]: number of deficits recorded at node k */
  int            *def_alloc;  /* def_alloc[k]: allocated size of deficits[k] */
  int             M;          /* number of HMM nodes (copy for convenience) */
  int             status;     /* eslOK or error code set by worker */
  char            errbuf[eslERRBUFSIZE]; /* error message if status != eslOK */
} CM_NODEPAD_WINFO;

/* cm_nodepad_thread_worker()
 * Worker function for threaded cm_ComputeP7CMNodePad().
 * Each worker pulls work items from the work queue. Each item contains
 * an emitted CM parsetree, emitted sequence, and flanked digitized
 * sequence already prepared by the master thread (the only caller of
 * the RNG). The worker runs Viterbi banding, records deficits into its
 * own per-node arrays, then frees cm_tr/esq/emb (ownership transferred
 * on pickup). Stops on sentinel (idx == -1).
 *
 * This design guarantees bit-for-bit determinism across thread counts
 * and scheduling: all RNG consumption happens on the master thread in
 * strict sample order, identical to the serial path.
 */
static void
cm_nodepad_thread_worker(void *arg)
{
  ESL_THREADS      *obj = (ESL_THREADS *) arg;
  int               workeridx;
  CM_NODEPAD_WINFO *winfo;
  CM_NODEPAD_WORK  *work = NULL;
  void             *newwork;
  CM_t             *cm;
  int               M;
  const int         flank = 500;
  int               seen[10000];  /* per-worker local array; avoids static thread-safety issue */

  esl_threads_Started(obj, &workeridx);
  winfo = (CM_NODEPAD_WINFO *) esl_threads_GetData(obj, workeridx);

  cm = winfo->cm;
  M  = winfo->M;

  winfo->status = eslOK;

  esl_workqueue_WorkerUpdate(winfo->queue, NULL, &newwork);
  work = (CM_NODEPAD_WORK *) newwork;

  while (work->idx >= 0)   /* sentinel: idx == -1 means stop */
    {
      Parsetree_t *cm_tr = work->cm_tr;
      ESL_SQ      *esq   = work->esq;
      ESL_DSQ     *emb   = work->emb;
      int          L     = work->L;
      int          L_emb = work->L_emb;
      int         *i2k   = NULL;
      int         *kmin  = NULL;
      int         *kmax  = NULL;
      int          ncells = 0;
      int          p, t;
      int          distinct_k, emit_pins, kk;

      p7_ProfileConfig(cm->fp7, winfo->bg, winfo->gm, L_emb, p7_GLOCAL);
      p7_gmx_GrowTo(winfo->gx, M, L_emb);

      if (p7_Seq2BandsVit(winfo->errbuf, winfo->gm, winfo->gx, winfo->bg, winfo->p7tr,
                          emb, L_emb, /*pad=*/0, /*nodepad=*/NULL, /*hopback=*/0, /*vitend=*/0,
                          &i2k, &kmin, &kmax, &ncells) != eslOK) {
        if (i2k)  free(i2k);
        if (kmin) free(kmin);
        if (kmax) free(kmax);
        free(emb);
        FreeParsetree(cm_tr); esl_sq_Destroy(esq);
        work->cm_tr = NULL; work->esq = NULL; work->emb = NULL;
        goto next_item;
      }

      /* Reject Viterbi traces that don't find the embedded emit (degenerate alignments). */
      distinct_k = 0;
      for (kk = 0; kk <= M && kk < 10000; kk++) seen[kk] = 0;
      for (p = 1; p <= L_emb; p++) {
        if (i2k[p] >= 1 && i2k[p] <= M && i2k[p] < 10000 && !seen[i2k[p]]) { seen[i2k[p]] = 1; distinct_k++; }
      }
      if (distinct_k < M / 2) {
        free(i2k); free(kmin); free(kmax); free(emb);
        FreeParsetree(cm_tr); esl_sq_Destroy(esq);
        work->cm_tr = NULL; work->esq = NULL; work->emb = NULL;
        goto next_item;
      }
      emit_pins = 0;
      for (p = flank + 1; p <= flank + L; p++) {
        if (i2k[p] >= 1 && i2k[p] <= M) emit_pins++;
      }
      if (emit_pins < L / 4 || ncells == 0) {
        free(i2k); free(kmin); free(kmax); free(emb);
        FreeParsetree(cm_tr); esl_sq_Destroy(esq);
        work->cm_tr = NULL; work->esq = NULL; work->emb = NULL;
        goto next_item;
      }

      /* Walk parsetree; record per-node deficits into worker-local arrays. */
      for (t = 0; t < cm_tr->n; t++) {
        int v    = cm_tr->state[t];
        int ipos = cm_tr->emitl[t];
        int jpos = cm_tr->emitr[t];
        int hn1  = cm->cp9map->cs2hn[v][0];
        int hn2  = cm->cp9map->cs2hn[v][1];
        int hs1  = cm->cp9map->cs2hs[v][0];
        int hs2  = cm->cp9map->cs2hs[v][1];

        #define CMNP_W_RECORD(true_k, pos_local) do {                                    \
          int _pos = (pos_local) + flank;                                                \
          if ((true_k) >= 1 && (true_k) <= M && _pos >= 1 && _pos <= L_emb) {           \
            int _def;                                                                    \
            if (kmin[_pos] == -1 || kmax[_pos] == -1)          _def = M;                \
            else if ((true_k) >= kmin[_pos] && (true_k) <= kmax[_pos]) _def = 0;        \
            else if ((true_k) < kmin[_pos])                    _def = kmin[_pos] - (true_k); \
            else                                               _def = (true_k) - kmax[_pos]; \
            if (winfo->def_n[(true_k)] >= winfo->def_alloc[(true_k)]) {                 \
              int _newsz = winfo->def_alloc[(true_k)] ? winfo->def_alloc[(true_k)] * 2 : 16; \
              int *_tmp = realloc(winfo->deficits[(true_k)], sizeof(int) * _newsz);     \
              if (_tmp == NULL) {                                                        \
                free(i2k); free(kmin); free(kmax); free(emb);                            \
                FreeParsetree(cm_tr); esl_sq_Destroy(esq);                               \
                work->cm_tr = NULL; work->esq = NULL; work->emb = NULL;                 \
                winfo->status = eslEMEM; goto worker_done;                               \
              }                                                                          \
              winfo->deficits[(true_k)] = _tmp;                                         \
              winfo->def_alloc[(true_k)] = _newsz;                                      \
            }                                                                           \
            winfo->deficits[(true_k)][winfo->def_n[(true_k)]++] = _def;                \
          }                                                                             \
        } while(0)

        if (hn1 >= 0) {
          int pos1;
          if      (hs1 == 0) { pos1 = (cm->stid[v] == MATR_MR) ? jpos : ipos; CMNP_W_RECORD(hn1, pos1); }
          else if (hs1 == 1) { pos1 = (cm->sttype[v] == IR_st)  ? jpos : ipos; CMNP_W_RECORD(hn1, pos1); }
        }
        if (hn2 >= 0) {
          if      (hs2 == 0) CMNP_W_RECORD(hn2, jpos);
          else if (hs2 == 1) CMNP_W_RECORD(hn2, jpos);
        }
        #undef CMNP_W_RECORD
      }

      free(i2k); free(kmin); free(kmax); free(emb);
      FreeParsetree(cm_tr); esl_sq_Destroy(esq);
      work->cm_tr = NULL; work->esq = NULL; work->emb = NULL;

    next_item:
      esl_workqueue_WorkerUpdate(winfo->queue, work, &newwork);
      work = (CM_NODEPAD_WORK *) newwork;
      continue;
    }

 worker_done:
  esl_workqueue_WorkerUpdate(winfo->queue, work, NULL);
  esl_threads_Finished(obj, workeridx);
  return;
}
#endif /* HMMER_THREADS */

/* Function: cm_ComputeP7CMNodePad()
 * Synopsis: Compute per-HMM-node p7 band pads by Monte Carlo simulation.
 *
 * Purpose:  Empirically derive per-node band-pad widths for the p7 HMM in <cm>.
 *           Emits <nsamples> parsetrees from the CM, embeds each in random
 *           flanking residues (500 nt each side), runs Viterbi banding with
 *           pad=0, and for every emitting CM state records the band deficit
 *           at the true HMM node. The <quantile> of per-node deficit
 *           distributions becomes the stored pad. Algorithm matches
 *           p7bandsim.c (the standalone CLI equivalent).
 *
 *           On success, populates cm->p7_cm_nodepad[0..fp7->M], sets
 *           cm->p7_cm_nodepad_M = fp7->M, and raises CMH_P7NODEPAD.
 *
 *           Caller must ensure cm_Configure() has been called and that
 *           cm->fp7 and cm->cp9map are valid.
 *
 *           Uses GLOCAL p7 profile (matches p7bandsim default). Expensive:
 *           ~nsamples * O(M*L) Viterbi DPs.
 *
 *           When <ncpu> > 1 and HMMER_THREADS is defined, the Monte Carlo
 *           simulation is parallelized across <ncpu> worker threads using
 *           the esl_threads/esl_workqueue pattern. Each worker thread runs
 *           its own independent RNG (seeded from the master <r>) so results
 *           differ from the serial path but are statistically equivalent.
 *           When <ncpu> <= 1 the serial code path is used unchanged.
 *
 * Returns:  <eslOK> on success.
 *           <eslEINVAL> if cm->fp7 or cm->cp9map is NULL (with errbuf message).
 *           <eslEMEM> on allocation failure.
 */
int
cm_ComputeP7CMNodePad(CM_t *cm, ESL_RANDOMNESS *r, int nsamples, double quantile, int ncpu, char *errbuf)
{
  int         status;
  P7_HMM     *hmm   = NULL;
  P7_BG      *bg    = NULL;
  P7_PROFILE *gm    = NULL;
  P7_GMX     *gx    = NULL;
  P7_TRACE   *p7tr  = NULL;
  int       **deficits  = NULL;
  int        *def_n     = NULL;
  int        *def_alloc = NULL;
  int        *pad_out   = NULL;
  int         M;
  int         k;
  const int   flank = 500;  /* p7bandsim default; flanks emitted residue each side */

  if (cm->fp7    == NULL) ESL_XFAIL(eslEINVAL, errbuf, "cm_ComputeP7CMNodePad: CM has no fp7 filter HMM");
  if (cm->cp9map == NULL) ESL_XFAIL(eslEINVAL, errbuf, "cm_ComputeP7CMNodePad: CM has no cp9map (was cm_Configure called?)");

  hmm = cm->fp7;
  M   = hmm->M;

#ifdef HMMER_THREADS
  if (ncpu > 1)
    {
      /* Threaded path: distribute nsamples across ncpu worker threads.
       * Pattern mirrors cm_p7_Tau() in cm_p7_modelmaker.c.
       * Each worker runs its own RNG (seeded from master r), accumulates
       * per-node deficit arrays locally, and the main thread merges them.
       */
      ESL_THREADS       *threadObj  = NULL;
      ESL_WORK_QUEUE    *queue      = NULL;
      CM_NODEPAD_WINFO  *winfo      = NULL;
      CM_NODEPAD_WORK   *wpool      = NULL;   /* recycling pool of work items */
      P7_BG             *master_bg  = NULL;   /* master's bg for FChoose flank emission */
      int                npool;
      int                next_s;             /* next sample index to dispatch */
      int                sentinels_sent;
      int                j;
      void              *newptr;
      CM_NODEPAD_WORK   *work;

      /* 1. Create per-worker data. Workers consume RNG-free, so no per-worker RNG. */
      ESL_ALLOC(winfo, sizeof(CM_NODEPAD_WINFO) * ncpu);
      for (j = 0; j < ncpu; j++) {
        winfo[j].cm      = cm;
        winfo[j].queue   = NULL;   /* filled below */
        winfo[j].bg      = p7_bg_Create(hmm->abc);
        winfo[j].gm      = p7_profile_Create(hmm->M, hmm->abc);
        p7_ProfileConfig(hmm, winfo[j].bg, winfo[j].gm, 400, p7_GLOCAL);
        winfo[j].gx      = p7_gmx_Create(hmm->M, 400);
        winfo[j].p7tr    = p7_trace_Create();
        winfo[j].M       = M;
        winfo[j].status  = eslOK;
        winfo[j].errbuf[0] = '\0';
        ESL_ALLOC(winfo[j].deficits,  sizeof(int *) * (M + 1));
        ESL_ALLOC(winfo[j].def_n,     sizeof(int)   * (M + 1));
        ESL_ALLOC(winfo[j].def_alloc, sizeof(int)   * (M + 1));
        for (k = 0; k <= M; k++) { winfo[j].deficits[k] = NULL; winfo[j].def_n[k] = 0; winfo[j].def_alloc[k] = 0; }
        if (winfo[j].bg == NULL || winfo[j].gm == NULL ||
            winfo[j].gx == NULL || winfo[j].p7tr == NULL) { status = eslEMEM; goto THREADED_ERROR; }
      }

      /* Master's own bg for FChoose flank emission — keeps RNG use on one thread. */
      master_bg = p7_bg_Create(hmm->abc);
      if (master_bg == NULL) { status = eslEMEM; goto THREADED_ERROR; }

      /* 2. Create recycling pool of work items. */
      npool = ncpu * 2;
      ESL_ALLOC(wpool, sizeof(CM_NODEPAD_WORK) * npool);
      for (j = 0; j < npool; j++) {
        wpool[j].idx   = -1;
        wpool[j].cm_tr = NULL;
        wpool[j].esq   = NULL;
        wpool[j].emb   = NULL;
        wpool[j].L     = 0;
        wpool[j].L_emb = 0;
      }

      /* 3. Set up threads and work queue. */
      threadObj = esl_threads_Create(&cm_nodepad_thread_worker);
      queue     = esl_workqueue_Create(npool);
      for (j = 0; j < ncpu; j++) {
        winfo[j].queue = queue;
        esl_threads_AddThread(threadObj, &winfo[j]);
      }

      /* 4. Initialize queue with empty (recycled) work items. */
      for (j = 0; j < npool; j++)
        esl_workqueue_Init(queue, &wpool[j]);

      /* 5. Reader loop: master emits parsetree+flanks in sample order (the only
       *    RNG consumer), then hands each fully-built work unit to a worker.
       *    Sentinels (idx=-1) shut down workers when we run out of samples.
       */
      esl_workqueue_Reset(queue);
      esl_threads_WaitForStart(threadObj);

      next_s         = 0;
      sentinels_sent = 0;

      status = esl_workqueue_ReaderUpdate(queue, NULL, &newptr);
      if (status != eslOK) goto THREADED_ERROR;
      work = (CM_NODEPAD_WORK *) newptr;

      while (sentinels_sent < ncpu)
        {
          if (next_s < nsamples) {
            Parsetree_t *cm_tr_local = NULL;
            ESL_SQ      *esq_local   = NULL;
            ESL_DSQ     *emb_local   = NULL;
            int          L_local;
            int          L_emb_local;
            char         name[32];
            int          p;

            snprintf(name, sizeof(name), "sim%d", next_s);
            if ((status = EmitParsetree(cm, errbuf, r, name, TRUE, &cm_tr_local, &esq_local, &L_local)) != eslOK)
              goto THREADED_ERROR;
            L_emb_local = 2 * flank + L_local;
            emb_local = malloc(sizeof(ESL_DSQ) * (L_emb_local + 2));
            if (emb_local == NULL) {
              FreeParsetree(cm_tr_local); esl_sq_Destroy(esq_local);
              status = eslEMEM; goto THREADED_ERROR;
            }
            emb_local[0] = emb_local[L_emb_local + 1] = eslDSQ_SENTINEL;
            for (p = 1;                    p <= flank;         p++) emb_local[p] = esl_rnd_FChoose(r, master_bg->f, hmm->abc->K);
            for (p = 1;                    p <= L_local;       p++) emb_local[flank + p] = esq_local->dsq[p];
            for (p = flank + L_local + 1;  p <= L_emb_local;   p++) emb_local[p] = esl_rnd_FChoose(r, master_bg->f, hmm->abc->K);

            work->idx   = next_s++;
            work->cm_tr = cm_tr_local;
            work->esq   = esq_local;
            work->emb   = emb_local;
            work->L     = L_local;
            work->L_emb = L_emb_local;
          } else {
            work->idx   = -1;   /* sentinel */
            work->cm_tr = NULL;
            work->esq   = NULL;
            work->emb   = NULL;
            sentinels_sent++;
          }
          status = esl_workqueue_ReaderUpdate(queue, work, &newptr);
          if (status != eslOK) goto THREADED_ERROR;
          work = (CM_NODEPAD_WORK *) newptr;
        }

      esl_threads_WaitForFinish(threadObj);
      esl_workqueue_Complete(queue);

      /* 6. Check for errors from workers. */
      for (j = 0; j < ncpu; j++) {
        if (winfo[j].status != eslOK) {
          ESL_XFAIL(winfo[j].status, errbuf, "%s", winfo[j].errbuf);
        }
      }

      /* 7. Merge per-worker deficit arrays into combined arrays. */
      ESL_ALLOC(deficits,  sizeof(int *) * (M + 1));
      ESL_ALLOC(def_n,     sizeof(int)   * (M + 1));
      ESL_ALLOC(def_alloc, sizeof(int)   * (M + 1));
      for (k = 0; k <= M; k++) { deficits[k] = NULL; def_n[k] = 0; def_alloc[k] = 0; }

      for (j = 0; j < ncpu; j++) {
        for (k = 1; k <= M; k++) {
          int n = winfo[j].def_n[k];
          if (n == 0) continue;
          if (def_n[k] + n > def_alloc[k]) {
            int newsz = def_alloc[k] + n;
            int *tmp = realloc(deficits[k], sizeof(int) * newsz);
            if (tmp == NULL) { status = eslEMEM; goto THREADED_ERROR; }
            deficits[k]  = tmp;
            def_alloc[k] = newsz;
          }
          memcpy(deficits[k] + def_n[k], winfo[j].deficits[k], sizeof(int) * n);
          def_n[k] += n;
        }
      }

      /* 8. Clean up thread resources. */
      for (j = 0; j < ncpu; j++) {
        p7_bg_Destroy(winfo[j].bg);
        p7_profile_Destroy(winfo[j].gm);
        p7_gmx_Destroy(winfo[j].gx);
        p7_trace_Destroy(winfo[j].p7tr);
        if (winfo[j].deficits != NULL) {
          int kk;
          for (kk = 0; kk <= M; kk++) if (winfo[j].deficits[kk] != NULL) free(winfo[j].deficits[kk]);
          free(winfo[j].deficits);
        }
        if (winfo[j].def_n)     free(winfo[j].def_n);
        if (winfo[j].def_alloc) free(winfo[j].def_alloc);
      }
      free(winfo);   winfo = NULL;
      free(wpool);   wpool = NULL;
      if (master_bg) p7_bg_Destroy(master_bg);
      esl_workqueue_Destroy(queue);
      esl_threads_Destroy(threadObj);

      goto AGGREGATE;

    THREADED_ERROR:
      /* Clean up thread resources on error. Also free any emissions in
       * work items that workers haven't consumed yet (owned by master on error).
       */
      if (wpool != NULL) {
        for (j = 0; j < npool; j++) {
          if (wpool[j].cm_tr != NULL) FreeParsetree(wpool[j].cm_tr);
          if (wpool[j].esq   != NULL) esl_sq_Destroy(wpool[j].esq);
          if (wpool[j].emb   != NULL) free(wpool[j].emb);
        }
      }
      if (winfo != NULL) {
        for (j = 0; j < ncpu; j++) {
          if (winfo[j].bg      != NULL) p7_bg_Destroy(winfo[j].bg);
          if (winfo[j].gm      != NULL) p7_profile_Destroy(winfo[j].gm);
          if (winfo[j].gx      != NULL) p7_gmx_Destroy(winfo[j].gx);
          if (winfo[j].p7tr    != NULL) p7_trace_Destroy(winfo[j].p7tr);
          if (winfo[j].deficits != NULL) {
            int kk;
            for (kk = 0; kk <= M; kk++) if (winfo[j].deficits[kk] != NULL) free(winfo[j].deficits[kk]);
            free(winfo[j].deficits);
          }
          if (winfo[j].def_n)     free(winfo[j].def_n);
          if (winfo[j].def_alloc) free(winfo[j].def_alloc);
        }
        free(winfo);
      }
      if (wpool     != NULL) free(wpool);
      if (master_bg != NULL) p7_bg_Destroy(master_bg);
      if (queue     != NULL) esl_workqueue_Destroy(queue);
      if (threadObj != NULL) esl_threads_Destroy(threadObj);
      goto ERROR;
    }
  else
#endif /* HMMER_THREADS */
    {
      /* Serial path (ncpu <= 1): original behavior, unchanged. */
      int s, t;

      bg   = p7_bg_Create(hmm->abc);
      gm   = p7_profile_Create(hmm->M, hmm->abc);
      p7_ProfileConfig(hmm, bg, gm, 400, p7_GLOCAL);  /* length reconfigured per-sample below */
      gx   = p7_gmx_Create(hmm->M, 400);
      p7tr = p7_trace_Create();

      ESL_ALLOC(deficits,  sizeof(int *) * (M + 1));
      ESL_ALLOC(def_n,     sizeof(int)   * (M + 1));
      ESL_ALLOC(def_alloc, sizeof(int)   * (M + 1));
      for (k = 0; k <= M; k++) { deficits[k] = NULL; def_n[k] = 0; def_alloc[k] = 0; }

      for (s = 0; s < nsamples; s++) {
        Parsetree_t *cm_tr = NULL;
        ESL_SQ      *esq   = NULL;
        int          L;
        char         name[32];
        ESL_DSQ     *emb   = NULL;
        int         *i2k   = NULL;
        int         *kmin  = NULL;
        int         *kmax  = NULL;
        int          ncells = 0;
        int          L_emb;
        int          p;
        int          distinct_k, emit_pins, kk;
        int          seen[10000];  /* local array; safe for serial path */

        snprintf(name, sizeof(name), "sim%d", s);
        if ((status = EmitParsetree(cm, errbuf, r, name, TRUE, &cm_tr, &esq, &L)) != eslOK) goto ERROR;

        L_emb = 2 * flank + L;
        ESL_ALLOC(emb, sizeof(ESL_DSQ) * (L_emb + 2));
        emb[0] = emb[L_emb + 1] = eslDSQ_SENTINEL;
        for (p = 1;              p <= flank;       p++) emb[p] = esl_rnd_FChoose(r, bg->f, hmm->abc->K);
        for (p = 1;              p <= L;           p++) emb[flank + p] = esq->dsq[p];
        for (p = flank + L + 1;  p <= L_emb;       p++) emb[p] = esl_rnd_FChoose(r, bg->f, hmm->abc->K);

        p7_ProfileConfig(hmm, bg, gm, L_emb, p7_GLOCAL);
        p7_gmx_GrowTo(gx, M, L_emb);

        if (p7_Seq2BandsVit(errbuf, gm, gx, bg, p7tr, emb, L_emb, /*pad=*/0, /*nodepad=*/NULL, /*hopback=*/0, /*vitend=*/0,
                            &i2k, &kmin, &kmax, &ncells) != eslOK) {
          if (i2k)  free(i2k);
          if (kmin) free(kmin);
          if (kmax) free(kmax);
          free(emb);
          FreeParsetree(cm_tr); esl_sq_Destroy(esq);
          continue;
        }

        /* Reject Viterbi traces that don't find the embedded emit (degenerate alignments). */
        distinct_k = 0;
        for (kk = 0; kk <= M && kk < 10000; kk++) seen[kk] = 0;
        for (p = 1; p <= L_emb; p++) {
          if (i2k[p] >= 1 && i2k[p] <= M && i2k[p] < 10000 && !seen[i2k[p]]) { seen[i2k[p]] = 1; distinct_k++; }
        }
        if (distinct_k < M / 2) {
          free(i2k); free(kmin); free(kmax); free(emb);
          FreeParsetree(cm_tr); esl_sq_Destroy(esq);
          continue;
        }
        emit_pins = 0;
        for (p = flank + 1; p <= flank + L; p++) {
          if (i2k[p] >= 1 && i2k[p] <= M) emit_pins++;
        }
        if (emit_pins < L / 4 || ncells == 0) {
          free(i2k); free(kmin); free(kmax); free(emb);
          FreeParsetree(cm_tr); esl_sq_Destroy(esq);
          continue;
        }

        /* Walk parsetree; record per-node deficits. */
        for (t = 0; t < cm_tr->n; t++) {
          int v    = cm_tr->state[t];
          int ipos = cm_tr->emitl[t];
          int jpos = cm_tr->emitr[t];
          int hn1  = cm->cp9map->cs2hn[v][0];
          int hn2  = cm->cp9map->cs2hn[v][1];
          int hs1  = cm->cp9map->cs2hs[v][0];
          int hs2  = cm->cp9map->cs2hs[v][1];

          #define CMNP_RECORD(true_k, pos_local) do {                                     \
            int pos = (pos_local) + flank;                                                \
            if ((true_k) >= 1 && (true_k) <= M && pos >= 1 && pos <= L_emb) {             \
              int def;                                                                    \
              if (kmin[pos] == -1 || kmax[pos] == -1)         def = M;                    \
              else if ((true_k) >= kmin[pos] && (true_k) <= kmax[pos]) def = 0;           \
              else if ((true_k) < kmin[pos])                  def = kmin[pos] - (true_k); \
              else                                            def = (true_k) - kmax[pos]; \
              if (def_n[(true_k)] >= def_alloc[(true_k)]) {                               \
                int newsz = def_alloc[(true_k)] ? def_alloc[(true_k)] * 2 : 16;           \
                int *tmp = realloc(deficits[(true_k)], sizeof(int) * newsz);              \
                if (tmp == NULL) { status = eslEMEM; goto ERROR; }                        \
                deficits[(true_k)] = tmp;                                                 \
                def_alloc[(true_k)] = newsz;                                              \
              }                                                                           \
              deficits[(true_k)][def_n[(true_k)]++] = def;                                \
            }                                                                             \
          } while(0)

          if (hn1 >= 0) {
            int pos1;
            if      (hs1 == 0) { pos1 = (cm->stid[v] == MATR_MR) ? jpos : ipos; CMNP_RECORD(hn1, pos1); }
            else if (hs1 == 1) { pos1 = (cm->sttype[v] == IR_st) ? jpos : ipos; CMNP_RECORD(hn1, pos1); }
          }
          if (hn2 >= 0) {
            if      (hs2 == 0) CMNP_RECORD(hn2, jpos);
            else if (hs2 == 1) CMNP_RECORD(hn2, jpos);
          }
          #undef CMNP_RECORD
        }

        free(i2k); free(kmin); free(kmax); free(emb);
        FreeParsetree(cm_tr); esl_sq_Destroy(esq);
      }
    } /* end serial path */

 AGGREGATE:
  /* Aggregate: per-node quantile of deficits. */
  {
    int idx_q;
    ESL_ALLOC(pad_out, sizeof(int) * (M + 1));
    pad_out[0] = 0;
    for (k = 1; k <= M; k++) {
      if (def_n[k] == 0) { pad_out[k] = 0; continue; }
      qsort(deficits[k], def_n[k], sizeof(int), cm_nodepad_cmpint);
      idx_q = (int)(quantile * (def_n[k] - 1) + 0.5);
      if (idx_q >= def_n[k]) idx_q = def_n[k] - 1;
      pad_out[k] = deficits[k][idx_q];
    }
  }

  /* Publish to CM. Replace any existing pad array. */
  if (cm->p7_cm_nodepad != NULL) free(cm->p7_cm_nodepad);
  cm->p7_cm_nodepad   = pad_out;
  cm->p7_cm_nodepad_M = M;
  cm->flags       |= CMH_P7NODEPAD;
  pad_out = NULL;  /* ownership transferred */

  status = eslOK;

 ERROR:
  if (pad_out != NULL)  free(pad_out);
  if (deficits != NULL) { for (k = 0; k <= M; k++) if (deficits[k] != NULL) free(deficits[k]); free(deficits); }
  if (def_n    != NULL) free(def_n);
  if (def_alloc != NULL) free(def_alloc);
  if (p7tr != NULL) p7_trace_Destroy(p7tr);
  if (gx   != NULL) p7_gmx_Destroy(gx);
  if (gm   != NULL) p7_profile_Destroy(gm);
  if (bg   != NULL) p7_bg_Destroy(bg);
  return status;
}


/*****************************************************************
 * Pin+Bridge banding (v5: SW-on-diagonal scan + LSIS pin selection)
 *
 * Replaces full unbanded p7_GViterbi in the --p7band path with:
 *   1. SW-on-diagonal SSE scan -> raw pins
 *   2. Longest Score-weighted Increasing Subsequence (LSIS) pin selection
 *   3. Pin diagonal segments + bridge rectangles + pad widening -> band
 *   4. Banded p7 GViterbi inside that prefilter band -> trace
 *   5. Trace -> M-state pins -> existing p7_pins2bands_nodepad (unchanged)
 *
 * See scratch_pinbridge/REPORT_v5.md for prototype validation.
 *****************************************************************/

/* Pin: a (k,i) anchor with priority r (SW score in v5).
 * r widened to int16_t 2026-05-29 (brief-079 S7 pinbridge prototype)
 * to support 16-bit SSE SW scan without saturation at L<200.
 * r widened again to int32_t 2026-06-01 (brief 26_0430-094) for the 32-bit SSE SW
 * scan: self-alignment-scale path scores exceed the int16 range (~1e7 at
 * dengue scale), which is exactly the saturation the 32-bit kernel removes.
 * Existing consumers cast r to int already, so the widening is transparent. */
typedef struct {
  int     k;
  int     i;
  int32_t r;
  int32_t r_peak; /* peak SW value over the segment (diagnostic; for PB_DUMP_PINS) */
} PB_Pin;

static inline int pb_adaptive_T(int M)
{
  if (M < 200) return 20;
  if (M < 500) return 30;
  return 40;
}

/* Function: p7_GBandedViterbi()
 *
 * Purpose:  Banded p7 Viterbi DP. Structurally identical to
 *           my_p7_GForwardBanded() above; the only difference is
 *           operator: p7_FLogsum() (probability-space sum) is
 *           replaced with ESL_MAX() (max-of-paths).
 *
 *           See my_p7_GForwardBanded() for band structure and
 *           overhang/segment-boundary arithmetic. NaN guards for
 *           truncated profiles, dpp pointer-stride fix, E->J via
 *           p7P_LOOP, and local mode support are preserved.
 *
 * Args:     dsq    - digital sequence, 1..L
 *           L      - length of sequence
 *           gm     - profile (LOCAL or GLOCAL)
 *           gxb    - banded DP matrix (gxb->bnd carries the band)
 *           opt_sc - optRETURN: Viterbi score in nats
 *
 * Returns:  <eslOK> on success.
 */
int
p7_GBandedViterbi(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMXB *gxb, float *opt_sc)
{
  int         *bnd_ip = gxb->bnd->imem;
  int         *bnd_kp = gxb->bnd->kmem;
  float       *dpc    = gxb->dp;
  float       *xpc    = gxb->xmx;
  float const *tsc    = gm->tsc;
  float const *rsc;
  float       *dpp;
  float       *last_dpc;
  int          ia, ib;
  int          last_ib;
  int          kac, kbc;
  int          kap, kbp;
  int          kbc2;
  float        xE, xN, xJ, xB, xC;
  float        mvp, ivp, dvp;
  float        dc;
  float        sc;
  int          g, i, k;
  float        esc  = p7_profile_IsLocal(gm) ? 0 : -eslINFINITY;

  xN      = 0.0f;
  xJ      = -eslINFINITY;
  xC      = -eslINFINITY;
  last_ib = 0;

  for (g = 0; g < gxb->bnd->nseg; g++)
    {
      ia = *bnd_ip++;
      ib = *bnd_ip++;

      kap = kbp = gm->M+1;
      dpp = dpc;

      xE  = -eslINFINITY;
      { int gap = ia - last_ib - 1;
        if (gap > 0) { xN = xN + gap * gm->xsc[p7P_N][p7P_LOOP];
                       xJ = xJ + gap * gm->xsc[p7P_J][p7P_LOOP];
                       xC = xC + gap * gm->xsc[p7P_C][p7P_LOOP]; } }
      xB  = ESL_MAX( xN + gm->xsc[p7P_N][p7P_MOVE], xJ + gm->xsc[p7P_J][p7P_MOVE]);

      for (i = ia; i <= ib; i++)
        {
          rsc      = gm->rsc[dsq[i]];
          dc       = -eslINFINITY;
          xE       = -eslINFINITY;
          last_dpc = dpc;

          kac      = *bnd_kp++;
          kbc      = *bnd_kp++;
          kbc2     = (kbc == gm->M ? kbc-1 : kbc);

          dpp += (kac-1 > kap ? ESL_MIN(kac-kap-1, kbp-kap+1) * p7G_NSCELLS : 0);

          if (kac > kap && kac-1 <= kbp) { mvp = *dpp++;       ivp = *dpp++;       dvp = *dpp++;       }
          else                           { mvp = -eslINFINITY; ivp = -eslINFINITY; dvp = -eslINFINITY; }

          for (k = kac; k <= kbc2; k++)
            {
              *dpc++ = sc = MSC(k) + ESL_MAX( ESL_MAX(mvp + TSC(p7P_MM, k-1), ivp + TSC(p7P_IM, k-1)),
                                              ESL_MAX(dvp + TSC(p7P_DM, k-1), xB  + TSC(p7P_BM, k-1)));

              if (k >= kap && k <= kbp) {  mvp = *dpp++;       ivp = *dpp++;        dvp = *dpp++;       }
              else                      {  mvp = -eslINFINITY; ivp = -eslINFINITY;  dvp = -eslINFINITY; }

              *dpc++ = ISC(k) + ESL_MAX( mvp + TSC(p7P_MI, k), ivp + TSC(p7P_II, k));

              xE     = ESL_MAX( ESL_MAX(sc + esc, dc + esc), xE);

              *dpc++ = dc;
              dc     = ESL_MAX( sc + TSC(p7P_MD, k), dc + TSC(p7P_DD, k));
            }

          if (kbc2 < kbc) /* k==M unrolled */
            {
              *dpc++ = sc = MSC(k) + ESL_MAX( ESL_MAX(mvp + TSC(p7P_MM, k-1), ivp + TSC(p7P_IM, k-1)),
                                              ESL_MAX(dvp + TSC(p7P_DM, k-1), xB  + TSC(p7P_BM, k-1)));
              *dpc++ = -eslINFINITY;
              *dpc++ = dc;
              xE     = ESL_MAX( ESL_MAX(sc, dc), xE);
            }

          *xpc++ = xE;
          *xpc++ = xN = xN + gm->xsc[p7P_N][p7P_LOOP];
          *xpc++ = xJ = ESL_MAX( xJ + gm->xsc[p7P_J][p7P_LOOP],  xE + gm->xsc[p7P_E][p7P_LOOP]);
          *xpc++ = xB = ESL_MAX( xJ + gm->xsc[p7P_J][p7P_MOVE],  xN + gm->xsc[p7P_N][p7P_MOVE]);
          *xpc++ = xC = ESL_MAX( xE + gm->xsc[p7P_E][p7P_MOVE],  xC + gm->xsc[p7P_C][p7P_LOOP]);

          dpp = last_dpc;
          kap = kac;
          kbp = kbc;
        }
      last_ib = ib;
    }

  if (opt_sc != NULL) { int tail = L - last_ib;
    *opt_sc = (tail > 0 ? xC + tail * gm->xsc[p7P_C][p7P_LOOP] : xC) + gm->xsc[p7P_C][p7P_MOVE]; }
  return eslOK;
}


/* Per-row offset table for accessing a filled P7_GMXB.
 * Populated by pb_build_row_offsets(); consumed by p7_GBandedTrace().
 */
typedef struct {
  int *kac;        /* [0..L] band kac for row i (row i in band: nonzero entry) */
  int *kbc;        /* [0..L] band kbc for row i */
  long *dp_off;    /* [0..L] offset into gxb->dp  for row i (in floats) */
  long *x_off;     /* [0..L] offset into gxb->xmx for row i (in floats) */
  char *in_band;   /* [0..L] 1 if row i is in some segment */
  int   ia0;       /* first ia (for prefix accounting) */
  int   ib_last;   /* last ib */
  float x0[p7G_NXCELLS]; /* row-0 specials (initial state, before any residue) */
} PB_RowMap;

static void pb_rowmap_destroy(PB_RowMap *rm)
{
  if (!rm) return;
  if (rm->kac) free(rm->kac);
  if (rm->kbc) free(rm->kbc);
  if (rm->dp_off) free(rm->dp_off);
  if (rm->x_off)  free(rm->x_off);
  if (rm->in_band) free(rm->in_band);
  free(rm);
}

static int pb_build_row_offsets(P7_GMXB *gxb, int L, PB_RowMap **ret_rm)
{
  PB_RowMap *rm = NULL;
  int g, i;
  int *bnd_ip = gxb->bnd->imem;
  int *bnd_kp = gxb->bnd->kmem;
  long dp_off = 0, x_off = 0;
  int  ia, ib, kac, kbc;
  int  status;

  ESL_ALLOC(rm, sizeof(PB_RowMap));
  rm->kac = NULL; rm->kbc = NULL; rm->dp_off = NULL; rm->x_off = NULL; rm->in_band = NULL;
  ESL_ALLOC(rm->kac,     sizeof(int)  * (L+2));
  ESL_ALLOC(rm->kbc,     sizeof(int)  * (L+2));
  ESL_ALLOC(rm->dp_off,  sizeof(long) * (L+2));
  ESL_ALLOC(rm->x_off,   sizeof(long) * (L+2));
  ESL_ALLOC(rm->in_band, sizeof(char) * (L+2));
  for (i = 0; i <= L+1; i++) { rm->kac[i] = 0; rm->kbc[i] = 0; rm->dp_off[i] = -1; rm->x_off[i] = -1; rm->in_band[i] = 0; }
  for (i = 0; i < p7G_NXCELLS; i++) rm->x0[i] = -eslINFINITY;

  rm->ia0 = (gxb->bnd->nseg > 0 ? gxb->bnd->imem[0] : L+1);
  rm->ib_last = 0;

  for (g = 0; g < gxb->bnd->nseg; g++) {
    ia = *bnd_ip++;
    ib = *bnd_ip++;
    rm->ib_last = ib;
    for (i = ia; i <= ib; i++) {
      kac = *bnd_kp++;
      kbc = *bnd_kp++;
      rm->kac[i]    = kac;
      rm->kbc[i]    = kbc;
      rm->dp_off[i] = dp_off;
      rm->x_off[i]  = x_off;
      rm->in_band[i]= 1;
      dp_off += (long)(kbc - kac + 1) * p7G_NSCELLS;
      x_off  += p7G_NXCELLS;
    }
  }
  *ret_rm = rm;
  return eslOK;

 ERROR:
  pb_rowmap_destroy(rm);
  return eslEMEM;
}

/* Inline accessors. M(i,k), I(i,k), D(i,k) return -eslINFINITY when out of band. */
static inline float pb_M(P7_GMXB *gxb, PB_RowMap *rm, int i, int k) {
  if (!rm->in_band[i] || k < rm->kac[i] || k > rm->kbc[i]) return -eslINFINITY;
  return gxb->dp[rm->dp_off[i] + (k - rm->kac[i]) * p7G_NSCELLS + p7G_M];
}
static inline float pb_I(P7_GMXB *gxb, PB_RowMap *rm, int i, int k) {
  if (!rm->in_band[i] || k < rm->kac[i] || k > rm->kbc[i]) return -eslINFINITY;
  return gxb->dp[rm->dp_off[i] + (k - rm->kac[i]) * p7G_NSCELLS + p7G_I];
}
static inline float pb_D(P7_GMXB *gxb, PB_RowMap *rm, int i, int k) {
  if (!rm->in_band[i] || k < rm->kac[i] || k > rm->kbc[i]) return -eslINFINITY;
  return gxb->dp[rm->dp_off[i] + (k - rm->kac[i]) * p7G_NSCELLS + p7G_D];
}
/* Specials: stored only on rows in band (1..L). Row 0 specials (initial
 * state, before any residue) are stashed on PB_RowMap by the caller and
 * returned here when i==0. Inter-segment rows outside any band return -inf. */
static inline float pb_X(P7_GMXB *gxb, PB_RowMap *rm, int i, int which) {
  if (i == 0) return rm->x0[which];
  if (rm->in_band[i]) return gxb->xmx[rm->x_off[i] + which];
  return -eslINFINITY;
}


/* Function: p7_GBandedTrace()
 *
 * Purpose:  Reconstruction-style traceback over a banded Viterbi matrix
 *           filled by p7_GBandedViterbi(). Mirrors p7_GTrace() with band
 *           guards on the per-row k-search.
 *
 *           Profile mode: GLOCAL or LOCAL. esc handling matches
 *           my_p7_GForwardBanded()/p7_GBandedViterbi().
 *
 *           For our integration, the band is a single segment spanning
 *           rows 1..L (no holes) — same shape produced by pinbridge —
 *           so inter-segment N/J/C accumulation does not arise.
 */
int
p7_GBandedTrace(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMXB *gxb, P7_TRACE *tr)
{
  int          i   = L;
  int          k   = 0;
  int          M   = gm->M;
  int          status;
  float        tol = 1e-5;
  float const *tsc = gm->tsc;
  int          sprv, scur;
  PB_RowMap   *rm = NULL;
  float        esc = p7_profile_IsLocal(gm) ? 0 : -eslINFINITY;

  if ((status = pb_build_row_offsets(gxb, L, &rm)) != eslOK) goto ERROR;

  /* Row 0 specials: initial conditions before any residue. p7_GBandedViterbi
   * does not store row 0 in xmx; reconstruct from gm parameters using the
   * same recurrence the Viterbi uses at start-of-segment when last_ib=0. */
  rm->x0[p7G_E] = -eslINFINITY;
  rm->x0[p7G_N] = 0.0f;
  rm->x0[p7G_J] = -eslINFINITY;
  rm->x0[p7G_B] = ESL_MAX( rm->x0[p7G_N] + gm->xsc[p7P_N][p7P_MOVE],
                           rm->x0[p7G_J] + gm->xsc[p7P_J][p7P_MOVE]);
  rm->x0[p7G_C] = -eslINFINITY;

  if ((status = p7_trace_Append(tr, p7T_T, 0, 0)) != eslOK) goto ERROR;
  if ((status = p7_trace_Append(tr, p7T_C, 0, L)) != eslOK) goto ERROR;
  sprv = p7T_C;
  while (sprv != p7T_S) {
    float const *rsc = (i>0 ? gm->rsc[dsq[i]] : NULL);

    switch (sprv) {
    case p7T_C:
      {
        float xC_i  = pb_X(gxb, rm, i, p7G_C);
        float xC_p  = (i > 0 ? pb_X(gxb, rm, i-1, p7G_C) : -eslINFINITY);
        float xE_i  = pb_X(gxb, rm, i, p7G_E);
        if (xC_i == -eslINFINITY) do { status = eslFAIL; goto ERROR; } while(0);
        if      (esl_FCompare_old(xC_i, xC_p + gm->xsc[p7P_C][p7P_LOOP], tol) == eslOK) scur = p7T_C;
        else if (esl_FCompare_old(xC_i, xE_i + gm->xsc[p7P_E][p7P_MOVE], tol) == eslOK) scur = p7T_E;
        else do { status = eslFAIL; goto ERROR; } while(0);
      }
      break;

    case p7T_E:
      {
        float xE_i = pb_X(gxb, rm, i, p7G_E);
        if (xE_i == -eslINFINITY) do { status = eslFAIL; goto ERROR; } while(0);
        if (p7_profile_IsLocal(gm))
          {
            scur = p7T_M;
            int kac = rm->in_band[i] ? rm->kac[i] : 1;
            int kbc = rm->in_band[i] ? rm->kbc[i] : 0;
            for (k = kbc; k >= kac; k--) {
              float Mk = pb_M(gxb, rm, i, k);
              if (esl_FCompare_old(xE_i, Mk + esc, tol) == eslOK) break;
              float Dk = pb_D(gxb, rm, i, k);
              if (esl_FCompare_old(xE_i, Dk + esc, tol) == eslOK) { scur = p7T_D; break; }
            }
            if (k < kac) { status = eslFAIL; goto ERROR; }
          }
        else
          {
            float MM = pb_M(gxb, rm, i, M);
            float DM = pb_D(gxb, rm, i, M);
            if      (esl_FCompare_old(xE_i, MM, tol) == eslOK) { scur = p7T_M; k = M; }
            else if (esl_FCompare_old(xE_i, DM, tol) == eslOK) { scur = p7T_D; k = M; }
            else { status = eslFAIL; goto ERROR; }
          }
      }
      break;

    case p7T_M:
      {
        float Mik = pb_M(gxb, rm, i, k);
        if (Mik == -eslINFINITY) do { status = eslFAIL; goto ERROR; } while(0);
        float xB_p = (i > 0 ? pb_X(gxb, rm, i-1, p7G_B) : -eslINFINITY);
        float Mp   = (i > 0 ? pb_M(gxb, rm, i-1, k-1) : -eslINFINITY);
        float Ip   = (i > 0 ? pb_I(gxb, rm, i-1, k-1) : -eslINFINITY);
        float Dp   = (i > 0 ? pb_D(gxb, rm, i-1, k-1) : -eslINFINITY);
        float msck = (rsc != NULL ? MSC(k) : 0);  /* MSC macro indexes rsc[] */
        if      (esl_FCompare_old(Mik, xB_p + TSC(p7P_BM, k-1) + msck, tol) == eslOK) scur = p7T_B;
        else if (esl_FCompare_old(Mik, Mp   + TSC(p7P_MM, k-1) + msck, tol) == eslOK) scur = p7T_M;
        else if (esl_FCompare_old(Mik, Ip   + TSC(p7P_IM, k-1) + msck, tol) == eslOK) scur = p7T_I;
        else if (esl_FCompare_old(Mik, Dp   + TSC(p7P_DM, k-1) + msck, tol) == eslOK) scur = p7T_D;
        else do { status = eslFAIL; goto ERROR; } while(0);
        k--; i--;
      }
      break;

    case p7T_D:
      {
        float Dik = pb_D(gxb, rm, i, k);
        if (Dik == -eslINFINITY) do { status = eslFAIL; goto ERROR; } while(0);
        float Mp = pb_M(gxb, rm, i, k-1);
        float Dp = pb_D(gxb, rm, i, k-1);
        if      (esl_FCompare_old(Dik, Mp + TSC(p7P_MD, k-1), tol) == eslOK) scur = p7T_M;
        else if (esl_FCompare_old(Dik, Dp + TSC(p7P_DD, k-1), tol) == eslOK) scur = p7T_D;
        else do { status = eslFAIL; goto ERROR; } while(0);
        k--;
      }
      break;

    case p7T_I:
      {
        float Iik = pb_I(gxb, rm, i, k);
        if (Iik == -eslINFINITY) do { status = eslFAIL; goto ERROR; } while(0);
        float Mp = (i > 0 ? pb_M(gxb, rm, i-1, k) : -eslINFINITY);
        float Ip = (i > 0 ? pb_I(gxb, rm, i-1, k) : -eslINFINITY);
        float isck = (rsc != NULL ? ISC(k) : 0);
        if      (esl_FCompare_old(Iik, Mp + TSC(p7P_MI, k) + isck, tol) == eslOK) scur = p7T_M;
        else if (esl_FCompare_old(Iik, Ip + TSC(p7P_II, k) + isck, tol) == eslOK) scur = p7T_I;
        else do { status = eslFAIL; goto ERROR; } while(0);
        i--;
      }
      break;

    case p7T_N:
      {
        float xN_i = pb_X(gxb, rm, i, p7G_N);
        if (xN_i == -eslINFINITY && i > 0) do { status = eslFAIL; goto ERROR; } while(0);
        scur = ((i == 0) ? p7T_S : p7T_N);
      }
      break;

    case p7T_B:
      {
        float xB_i = pb_X(gxb, rm, i, p7G_B);
        if (xB_i == -eslINFINITY) do { status = eslFAIL; goto ERROR; } while(0);
        float xN_i = pb_X(gxb, rm, i, p7G_N);
        float xJ_i = pb_X(gxb, rm, i, p7G_J);
        if      (esl_FCompare_old(xB_i, xN_i + gm->xsc[p7P_N][p7P_MOVE], tol) == eslOK) scur = p7T_N;
        else if (esl_FCompare_old(xB_i, xJ_i + gm->xsc[p7P_J][p7P_MOVE], tol) == eslOK) scur = p7T_J;
        else do { status = eslFAIL; goto ERROR; } while(0);
      }
      break;

    case p7T_J:
      {
        float xJ_i = pb_X(gxb, rm, i, p7G_J);
        if (xJ_i == -eslINFINITY) do { status = eslFAIL; goto ERROR; } while(0);
        float xJ_p = (i > 0 ? pb_X(gxb, rm, i-1, p7G_J) : -eslINFINITY);
        float xE_i = pb_X(gxb, rm, i, p7G_E);
        if      (esl_FCompare_old(xJ_i, xJ_p + gm->xsc[p7P_J][p7P_LOOP], tol) == eslOK) scur = p7T_J;
        else if (esl_FCompare_old(xJ_i, xE_i + gm->xsc[p7P_E][p7P_LOOP], tol) == eslOK) scur = p7T_E;
        else do { status = eslFAIL; goto ERROR; } while(0);
      }
      break;

    default: do { status = eslFAIL; goto ERROR; } while(0);
    }

    if ((status = p7_trace_Append(tr, scur, k, i)) != eslOK) goto ERROR;
    if ( (scur == p7T_N || scur == p7T_J || scur == p7T_C) && scur == sprv) i--;
    sprv = scur;
  }

  tr->M = gm->M;
  tr->L = L;
  pb_rowmap_destroy(rm);
  return p7_trace_Reverse(tr);

 ERROR:
  pb_rowmap_destroy(rm);
  return status;
}


/*****************************************************************
 * Brief 26_0430-094 component knobs: env-gated defaults.
 *
 * Each of the three brief-094 mechanisms is gated behind an env var so it
 * can be A/B'd independently against the brief-091 baseline. The env var (if
 * set) overrides the compile-time default; with the env unset the default
 * applies. The defaults below are flipped ON in the final brief-094 commit;
 * setting PB_USE_32BIT_SW=0 / PB_LSIS_WINDOW=0 / PB_K_ADAPTIVE=0 reaches the
 * legacy brief-091 behavior at runtime.
 *****************************************************************/
#define PB_DEFAULT_USE_32BIT_SW 1    /* 1 = 32-bit SSE SW; 0 = legacy 16-bit  */
#define PB_DEFAULT_LSIS_WINDOW  200  /* >0 = limited-window LSIS (model-pos window); 0 = full O(N^2) */
#define PB_DEFAULT_K_ADAPTIVE   1    /* 1 = K = max(PB_DEFAULT_TOPK, L/(2*pad)); 0 = fixed PB_DEFAULT_TOPK */

/* Model-size gate for the brief-094 mechanisms (2026-06-01 validation).
 *
 * The 32-bit SW eliminates a band-size BLOWUP that only occurs when the
 * 16-bit SW accumulation saturates AND the model is large enough for the
 * resulting over-wide band to be expensive (dengue M=10735: 99M -> 0.48M
 * cells, 208x; HSV M=152222: enables pinbridge where 16-bit's saturated
 * band fails outright). At small M (rmark4e/4h, M <= ~2900) there is no
 * blowup to fix, and the 32-bit kernel instead (a) costs ~2x in the SW scan
 * (4 int32 lanes vs 8 int16) and (b) perturbs pin selection toward slightly
 * WIDER bands (rmark4e summed CM matrix 4.1 -> 5.0 GB, LSU 847 -> 1804 MB),
 * with a small accuracy cost (rmark4e 0.99242 -> 0.99124, rmark4h 0.97887 ->
 * 0.97688). So the mechanisms are gated ON only for M >= PB_M32_MIN, which
 * keeps the small-M path byte-identical to brief 26_0430-091 while capturing the
 * large-M win. An explicitly-set env var overrides the gate at any M (for
 * A/B and power users). The 2900..10735 crossover is uncalibrated; 8000 is a
 * conservative default that cleanly separates Rfam-scale from viral-genome
 * scale (the VADR target). */
#define PB_M32_MIN_DEFAULT 8000

/* Read an integer env var, falling back to <defval> when unset. */
static int
pb_env_int(const char *name, int defval)
{
  const char *s = getenv(name);
  return (s != NULL) ? atoi(s) : defval;
}

/* Resolve a brief-094 mechanism gate: an explicitly-set env var wins at any
 * M; otherwise the compile-time default applies only when m_ok (M-gate). For
 * the non-boolean window knob this returns <compiled_default> (e.g. 200) or 0. */
static int
pb_gate(const char *name, int compiled_default, int m_ok)
{
  const char *s = getenv(name);
  if (s != NULL) return atoi(s);
  return m_ok ? compiled_default : 0;
}

/*****************************************************************
 * CM_PB_OM32: Infernal-side 32-bit striped emission table (brief 26_0430-094).
 *
 * Mirrors the allocation/striping layout of HMMER's P7_OPROFILE->rwv
 * (impl_sse/p7_oprofile.c) but with int32 lanes (p7O_NQF = 4 lanes/vector)
 * instead of int16 (p7O_NQW = 8 lanes/vector). Built by READING the float
 * match-emission scores from P7_PROFILE->rsc (p7P_MSC) with the *same*
 * scale HMMER's wordify() uses (scale_w = 500/log2), so at small M (no DP
 * saturation) the 32-bit kernel produces pin scores byte-identical to the
 * 16-bit kernel. The only difference is at large M / self-alignment scale,
 * where the 16-bit accumulation saturates at +32767 and the 32-bit one does
 * not. HMMER is never modified; we only read P7_PROFILE.
 *****************************************************************/
struct cm_pb_om32_s {
  int       M;        /* model length                                        */
  int       Q;        /* p7O_NQF(M) = ceil(M/4), # of 4-lane int32 vectors    */
  int       Kp;       /* alphabet size (abc->Kp); riv has Kp rows            */
  float     scale_i;  /* nats -> 32-bit units (= HMMER scale_w = 500/log2)    */
  __m128i **riv;      /* [Kp][Q] striped 32-bit match-emission scores         */
  __m128i  *riv_mem;  /* backing allocation for riv                           */
};

/* Padding lanes (model position k > M) get a large-negative emission so the
 * SW recurrence max(0, mpv+sc) floors them to 0 and they never spawn a pin.
 * -1e9 dominates any reachable accumulation (< ~3e8 even at L=150K). */
#define PB_OM32_NEGINF (-1000000000)

static void
cm_pb_om32_Destroy(CM_PB_OM32 *om32)
{
  if (om32 == NULL) return;
  if (om32->riv_mem != NULL) free(om32->riv_mem);
  if (om32->riv     != NULL) free(om32->riv);
  free(om32);
}

static int
cm_pb_om32_Create(int M, const ESL_ALPHABET *abc, CM_PB_OM32 **ret_om32)
{
  CM_PB_OM32 *om32 = NULL;
  int         Q    = p7O_NQF(M);
  int         Kp   = abc->Kp;
  int         x;
  int         status;

  ESL_ALLOC(om32, sizeof(CM_PB_OM32));
  om32->M = M; om32->Q = Q; om32->Kp = Kp; om32->scale_i = 0.0f;
  om32->riv = NULL; om32->riv_mem = NULL;

  /* +15 slack so riv[0] can be bumped to a 16-byte boundary (mirrors HMMER). */
  ESL_ALLOC(om32->riv_mem, sizeof(__m128i) * Q * Kp + 15);
  ESL_ALLOC(om32->riv,     sizeof(__m128i *) * Kp);
  om32->riv[0] = (__m128i *) (((unsigned long int) om32->riv_mem + 15) & (~0xful));
  for (x = 1; x < Kp; x++) om32->riv[x] = om32->riv[0] + (x * Q);

  *ret_om32 = om32;
  return eslOK;

 ERROR:
  cm_pb_om32_Destroy(om32);
  *ret_om32 = NULL;
  return status;
}

/* Fill om32->riv from the LOCAL profile <gm>'s float match-emission scores.
 * Lane (q,z) holds model position k = (q+1) + z*Q (the same striping the SW
 * kernel and HMMER's rwv use). scale_i is set to HMMER's wordify scale so
 * the integer scores match the 16-bit table for k<=M positions. */
static int
cm_pb_om32_Build(const P7_PROFILE *gm, CM_PB_OM32 *om32)
{
  int M  = om32->M;
  int Q  = om32->Q;
  int Kp = om32->Kp;
  int x, q, z;
  union { __m128i v; int32_t i[4]; } tmp;

  om32->scale_i = 500.0 / eslCONST_LOG2;   /* == P7_OPROFILE scale_w */

  for (x = 0; x < Kp; x++) {
    for (q = 0; q < Q; q++) {
      for (z = 0; z < 4; z++) {
        int k = (q + 1) + z * Q;
        if (k <= M) {
          float sc = roundf(om32->scale_i * p7P_MSC(gm, k, x));
          /* defensive int32 clamp; emissions never approach this */
          if      (sc >=  2147483520.0f) tmp.i[z] =  2147483520;
          else if (sc <= -2147483520.0f) tmp.i[z] = -2147483520;
          else                           tmp.i[z] = (int32_t) sc;
        } else {
          tmp.i[z] = PB_OM32_NEGINF;
        }
      }
      om32->riv[x][q] = tmp.v;
    }
  }
  return eslOK;
}


/*****************************************************************
 * SW-on-diagonal pin scan + LSIS selection + band construction
 * (ports of scratch_pinbridge/pinbridge_proto.c v5)
 *****************************************************************/

__attribute__((target("sse4.1")))
static int
pb_sw_scan_collect_pins(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int T,
                        PB_Pin **ret_pins, int *ret_npins)
{
  int M = om->M;
  int Q = p7O_NQB(M);
  int i, q, z;

  int max_pins = (M + L) * 4;
  PB_Pin *pins = (PB_Pin *) malloc(max_pins * sizeof(PB_Pin));
  if (!pins) return eslEMEM;
  int npins = 0;

  __m128i *prev = NULL;
  __m128i *curr = NULL;
  if (posix_memalign((void **)&prev, 16, Q * sizeof(__m128i)) != 0) { free(pins); return eslEMEM; }
  if (posix_memalign((void **)&curr, 16, Q * sizeof(__m128i)) != 0) { free(prev); free(pins); return eslEMEM; }

  for (q = 0; q < Q; q++) { prev[q] = _mm_setzero_si128(); curr[q] = _mm_setzero_si128(); }

  __m128i biasv = _mm_set1_epi8((int8_t) om->bias_b);
  __m128i zero  = _mm_setzero_si128();

  for (i = 1; i <= L; i++) {
    __m128i const *rsc = om->rbv[dsq[i]];
    __m128i mpv = _mm_slli_si128(prev[Q - 1], 1);

    for (q = 0; q < Q; q++) {
      __m128i sc     = _mm_sub_epi8(rsc[q], biasv);
      __m128i cand   = _mm_adds_epi8(mpv, sc);
      __m128i newval = _mm_max_epi8(zero, cand);
      mpv      = prev[q];
      curr[q]  = newval;
    }

    /* Collect end-of-segment pins at row (i-1) */
    {
      union { __m128i v; uint8_t b[16]; } u_prev, u_curr;
      for (q = 0; q < Q; q++) {
        u_prev.v = prev[q];
        u_curr.v = curr[q];
        for (z = 0; z < 16; z++) {
          int k = (q + 1) + z * Q;
          if (k > M) break;
          if ((int8_t)u_prev.b[z] >= (int8_t)T && u_curr.b[z] == 0) {
            if (npins >= max_pins) {
              max_pins *= 2;
              PB_Pin *tmp = (PB_Pin *) realloc(pins, max_pins * sizeof(PB_Pin));
              if (!tmp) { free(pins); free(prev); free(curr); return eslEMEM; }
              pins = tmp;
            }
            pins[npins].k = k;
            pins[npins].i = i - 1;
            pins[npins].r = (int16_t)(int8_t)u_prev.b[z];
            npins++;
          }
        }
      }
    }

    __m128i *tmp = prev; prev = curr; curr = tmp;
    for (q = 0; q < Q; q++) curr[q] = _mm_setzero_si128();
  }

  /* Last-row pins */
  {
    union { __m128i v; uint8_t b[16]; } u_prev;
    for (q = 0; q < Q; q++) {
      u_prev.v = prev[q];
      for (z = 0; z < 16; z++) {
        int k = (q + 1) + z * Q;
        if (k > M) break;
        if ((int8_t)u_prev.b[z] >= (int8_t)T) {
          if (npins >= max_pins) {
            max_pins *= 2;
            PB_Pin *tmp2 = (PB_Pin *) realloc(pins, max_pins * sizeof(PB_Pin));
            if (!tmp2) { free(pins); free(prev); free(curr); return eslEMEM; }
            pins = tmp2;
          }
          pins[npins].k = k;
          pins[npins].i = L;
          pins[npins].r = (int16_t)(int8_t)u_prev.b[z];
          npins++;
        }
      }
    }
  }

  free(prev);
  free(curr);

  *ret_pins  = pins;
  *ret_npins = npins;
  return eslOK;
}

/*****************************************************************
 * Runtime pin-emission mode selector (brief 26_0430-101, 2026-06-02)
 *
 * PB_PIN_EMIT_MODE=trail (default): emit at trailing edge (i-1, k)
 *   when segment resets — production behavior from pinbridge-32bit.
 * PB_PIN_EMIT_MODE=peak: emit at peak position (i_peak, k) with
 *   score = peak_value, including an end-of-scan flush for sticky
 *   segments that never reset.
 *****************************************************************/
enum { PIN_MODE_TRAIL = 0, PIN_MODE_PEAK = 1 };

static int pb_get_emit_mode(void)
{
  static int cached = -1;
  if (cached < 0) {
    const char *s = getenv("PB_PIN_EMIT_MODE");
    if (s && strcmp(s, "peak") == 0) cached = PIN_MODE_PEAK;
    else                             cached = PIN_MODE_TRAIL;
  }
  return cached;
}

/* Brief 26_0430-102: gate end-of-scan flush separately to isolate from peak-on-reset. */
static int pb_get_flush_enabled(void)
{
  static int cached = -1;
  if (cached < 0) {
    const char *s = getenv("PB_PIN_FLUSH");
    cached = (s && strcmp(s, "1") == 0) ? 1 : 0;
  }
  return cached;
}

/*****************************************************************
 * 16-bit SSE SW scan (prototype 2026-05-29, brief 26_0430-079 S7 follow-up)
 *
 * The 8-bit version above saturates at +127 for L<200 on conserved
 * targets like RF00010, destroying LSIS's chain-discrimination
 * ability. This 16-bit variant uses the same striped SSE pattern
 * with epi16 ops + om->rwv (HMMER ViterbiFilter's emission table).
 * 8 lanes per register instead of 16 → ~2x wall, but >100x value
 * range, so saturation is not an issue at M <= ~10K typical scale.
 * Threshold T_w scaled to 16-bit emission units.
 *****************************************************************/
__attribute__((target("sse4.1")))
static int
pb_sw_scan_collect_pins_w(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int T_w,
                          PB_Pin **ret_pins, int *ret_npins)
{
  int M = om->M;
  int Q = p7O_NQW(M);    /* 8 lanes (16-bit words) per 128-bit register */
  int i, q, z;
  int pin_mode = pb_get_emit_mode();
  if (pin_mode == PIN_MODE_PEAK && L > 65000) pin_mode = PIN_MODE_TRAIL; /* 16-bit i position guard */

  int max_pins = (M + L) * 4;
  PB_Pin *pins = (PB_Pin *) malloc(max_pins * sizeof(PB_Pin));
  if (!pins) return eslEMEM;
  int npins = 0;

  __m128i *prev = NULL, *curr = NULL;
  __m128i *peak = NULL, *peak_i = NULL, *peak_prev = NULL, *peak_i_prev = NULL;

  if (posix_memalign((void **)&prev,       16, Q * sizeof(__m128i)) != 0) { free(pins); return eslEMEM; }
  if (posix_memalign((void **)&curr,       16, Q * sizeof(__m128i)) != 0) { free(prev); free(pins); return eslEMEM; }
  if (posix_memalign((void **)&peak,       16, Q * sizeof(__m128i)) != 0) { free(curr); free(prev); free(pins); return eslEMEM; }
  if (posix_memalign((void **)&peak_i,     16, Q * sizeof(__m128i)) != 0) { free(peak); free(curr); free(prev); free(pins); return eslEMEM; }
  if (posix_memalign((void **)&peak_prev,  16, Q * sizeof(__m128i)) != 0) { free(peak_i); free(peak); free(curr); free(prev); free(pins); return eslEMEM; }
  if (posix_memalign((void **)&peak_i_prev,16, Q * sizeof(__m128i)) != 0) { free(peak_prev); free(peak_i); free(peak); free(curr); free(prev); free(pins); return eslEMEM; }

  for (q = 0; q < Q; q++) {
    prev[q] = _mm_setzero_si128(); curr[q] = _mm_setzero_si128();
    peak[q] = _mm_setzero_si128(); peak_i[q] = _mm_setzero_si128();
    peak_prev[q] = _mm_setzero_si128(); peak_i_prev[q] = _mm_setzero_si128();
  }

  __m128i zero = _mm_setzero_si128();

  /* Brief 26_0430-104 diagnostic: dump the full SW accumulator matrix curr[i][k].
   * PB_DUMP_SW_MX=<path> writes one TSV line per (i,k) cell: i  k  curr.
   * Dumped after the DP update for row i, before pin emission.
   * Use --cpu 1 so concurrent threads don't interleave lines.
   * DIAGNOSTIC-ONLY: no effect on pin emission. */
  FILE *sw_dump_fp = NULL;
  {
    const char *sw_dump_path = getenv("PB_DUMP_SW_MX");
    if (sw_dump_path != NULL) {
      sw_dump_fp = fopen(sw_dump_path, "w");
      if (sw_dump_fp != NULL) {
        fprintf(sw_dump_fp, "# SW accumulator matrix dump (brief 104), kernel=w16\n");
        fprintf(sw_dump_fp, "# M=%d L=%d Q=%d T_w=%d\n", M, L, Q, T_w);
        fprintf(sw_dump_fp, "#i\tk\tcurr\n");
      }
    }
  }

  for (i = 1; i <= L; i++) {
    __m128i const *rsc = om->rwv[dsq[i]];                /* 16-bit emission scores */
    __m128i mpv = _mm_slli_si128(prev[Q - 1], 2);        /* shift left 2 bytes = 1 word */

    for (q = 0; q < Q; q++) {
      __m128i cand   = _mm_adds_epi16(mpv, rsc[q]);      /* SW extend (signed saturating) */
      __m128i newval = _mm_max_epi16(zero, cand);        /* SW floor at 0 */
      mpv     = prev[q];
      curr[q] = newval;
    }

    /* Brief 26_0430-104 diagnostic: dump row i of the SW accumulator (curr[]),
     * unstriped to per-k order. */
    if (sw_dump_fp != NULL) {
      union { __m128i v; int16_t b[8]; } u_dump;
      for (q = 0; q < Q; q++) {
        u_dump.v = curr[q];
        for (z = 0; z < 8; z++) {
          int k = (q + 1) + z * Q;
          if (k > M) break;
          fprintf(sw_dump_fp, "%d\t%d\t%d\n", i, k, (int)u_dump.b[z]);
        }
      }
    }

    /* Peak tracking: snapshot previous peak, update with curr, reset on segment end */
    for (q = 0; q < Q; q++) {
      peak_prev[q]   = peak[q];
      peak_i_prev[q] = peak_i[q];
      __m128i is_new_peak = _mm_cmpgt_epi16(curr[q], peak[q]);
      peak[q]   = _mm_max_epi16(peak[q], curr[q]);
      peak_i[q] = _mm_blendv_epi8(peak_i[q], _mm_set1_epi16((int16_t)i), is_new_peak);
      __m128i is_reset = _mm_cmpeq_epi16(curr[q], zero);
      peak[q]   = _mm_andnot_si128(is_reset, peak[q]);
      peak_i[q] = _mm_andnot_si128(is_reset, peak_i[q]);
    }

    /* Collect end-of-segment pins */
    if (pin_mode == PIN_MODE_TRAIL) {
      union { __m128i v; int16_t b[8]; } u_prev, u_curr, u_pk;
      for (q = 0; q < Q; q++) {
        u_prev.v = prev[q];
        u_curr.v = curr[q];
        u_pk.v   = peak_prev[q];
        for (z = 0; z < 8; z++) {
          int k = (q + 1) + z * Q;
          if (k > M) break;
          if (u_prev.b[z] >= (int16_t)T_w && u_curr.b[z] == 0) {
            if (npins >= max_pins) {
              max_pins *= 2;
              PB_Pin *tmp = (PB_Pin *) realloc(pins, max_pins * sizeof(PB_Pin));
              if (!tmp) { free(pins); free(prev); free(curr); free(peak); free(peak_i); free(peak_prev); free(peak_i_prev); return eslEMEM; }
              pins = tmp;
            }
            pins[npins].k      = k;
            pins[npins].i      = i - 1;
            pins[npins].r      = u_prev.b[z];
            pins[npins].r_peak = u_pk.b[z];
            npins++;
          }
        }
      }
    } else { /* PIN_MODE_PEAK */
      union { __m128i v; int16_t b[8]; } u_pk, u_pki, u_curr;
      for (q = 0; q < Q; q++) {
        u_pk.v  = peak_prev[q];
        u_pki.v = peak_i_prev[q];
        u_curr.v = curr[q];
        for (z = 0; z < 8; z++) {
          int k = (q + 1) + z * Q;
          if (k > M) break;
          if (u_pk.b[z] >= (int16_t)T_w && u_curr.b[z] == 0) {
            if (npins >= max_pins) {
              max_pins *= 2;
              PB_Pin *tmp = (PB_Pin *) realloc(pins, max_pins * sizeof(PB_Pin));
              if (!tmp) { free(pins); free(prev); free(curr); free(peak); free(peak_i); free(peak_prev); free(peak_i_prev); return eslEMEM; }
              pins = tmp;
            }
            pins[npins].k      = k;
            pins[npins].i      = (uint16_t)u_pki.b[z];
            pins[npins].r      = u_pk.b[z];
            pins[npins].r_peak = u_pk.b[z];
            npins++;
          }
        }
      }
    }

    __m128i *tmp = prev; prev = curr; curr = tmp;
    for (q = 0; q < Q; q++) curr[q] = _mm_setzero_si128();
  }

  /* Last-row pins (trail mode only; peak mode uses end-of-scan flush below) */
  if (pin_mode == PIN_MODE_TRAIL) {
    union { __m128i v; int16_t b[8]; } u_prev, u_pk;
    for (q = 0; q < Q; q++) {
      u_prev.v = prev[q];
      u_pk.v   = peak[q];
      for (z = 0; z < 8; z++) {
        int k = (q + 1) + z * Q;
        if (k > M) break;
        if (u_prev.b[z] >= (int16_t)T_w) {
          if (npins >= max_pins) {
            max_pins *= 2;
            PB_Pin *tmp2 = (PB_Pin *) realloc(pins, max_pins * sizeof(PB_Pin));
            if (!tmp2) { free(pins); free(prev); free(curr); free(peak); free(peak_i); free(peak_prev); free(peak_i_prev); return eslEMEM; }
            pins = tmp2;
          }
          pins[npins].k      = k;
          pins[npins].i      = L;
          pins[npins].r      = u_prev.b[z];
          pins[npins].r_peak = u_pk.b[z];
          npins++;
        }
      }
    }
  }

  /* End-of-scan flush (peak mode): emit for sticky segments still active at i=L */
  /* Brief 26_0430-102: flush is gated separately via PB_PIN_FLUSH to isolate from peak-on-reset. */
  if (pin_mode == PIN_MODE_PEAK && pb_get_flush_enabled()) {
    union { __m128i v; int16_t b[8]; } u_pk, u_pki, u_prev;
    for (q = 0; q < Q; q++) {
      u_pk.v   = peak[q];
      u_pki.v  = peak_i[q];
      u_prev.v = prev[q];    /* SW values at row L (after final swap) */
      for (z = 0; z < 8; z++) {
        int k = (q + 1) + z * Q;
        if (k > M) break;
        if (u_pk.b[z] >= (int16_t)T_w && u_prev.b[z] > 0) {
          if (npins >= max_pins) {
            max_pins *= 2;
            PB_Pin *tmp3 = (PB_Pin *) realloc(pins, max_pins * sizeof(PB_Pin));
            if (!tmp3) { free(pins); free(prev); free(curr); free(peak); free(peak_i); free(peak_prev); free(peak_i_prev); return eslEMEM; }
            pins = tmp3;
          }
          pins[npins].k      = k;
          pins[npins].i      = (uint16_t)u_pki.b[z];
          pins[npins].r      = u_pk.b[z];
          pins[npins].r_peak = u_pk.b[z];
          npins++;
        }
      }
    }
  }

  free(prev); free(curr); free(peak); free(peak_i); free(peak_prev); free(peak_i_prev);

  if (sw_dump_fp != NULL) fclose(sw_dump_fp);   /* brief 26_0430-104 SW matrix dump */

  *ret_pins  = pins;
  *ret_npins = npins;
  return eslOK;
}

/*****************************************************************
 * 32-bit SSE SW scan (brief 26_0430-094, 2026-06-01)
 *
 * The 16-bit version above saturates the SW *accumulation* at +32767 on
 * highly-conserved / self-alignment-scale targets (e.g. dengue self-align,
 * L=M=10735: all raw pins pegged at 32767). With every pin tied at the max,
 * score-based Top-K pruning degenerates to tie-break order and concentrates
 * pins in a tiny seq-position window, which the band coda then widens into a
 * near-full L*M matrix. This 32-bit variant uses int32 lanes (4 per register
 * via p7O_NQF) and a *non-saturating* _mm_add_epi32, so accumulated path
 * scores keep their true magnitude and Top-K ranking stays meaningful.
 *
 * Mirrors pb_sw_scan_collect_pins_w exactly except:
 *   - om->rwv (int16, 8 lanes)        -> om32->riv (int32, 4 lanes)
 *   - _mm_slli_si128(..., 2)          -> _mm_slli_si128(..., 4)   (1 int32 lane)
 *   - _mm_adds_epi16 (saturating)     -> _mm_add_epi32  (NON-saturating: the fix)
 *   - _mm_max_epi16                   -> _mm_max_epi32  (SSE4.1)
 *   - 8-lane union                    -> 4-lane union
 * At small M (no accumulation reaches 32767) it produces pins byte-identical
 * to the 16-bit kernel, since riv carries the same wordify-scaled scores.
 *****************************************************************/
__attribute__((target("sse4.1")))
static int
pb_sw_scan_collect_pins_i32(const ESL_DSQ *dsq, int L, const CM_PB_OM32 *om32, int T_i32,
                            PB_Pin **ret_pins, int *ret_npins)
{
  int M = om32->M;
  int Q = om32->Q;    /* 4 lanes (32-bit ints) per 128-bit register */
  int i, q, z;
  int pin_mode = pb_get_emit_mode();

  int max_pins = (M + L) * 4;
  PB_Pin *pins = (PB_Pin *) malloc(max_pins * sizeof(PB_Pin));
  if (!pins) return eslEMEM;
  int npins = 0;

  __m128i *prev = NULL, *curr = NULL;
  __m128i *peak = NULL, *peak_i = NULL, *peak_prev = NULL, *peak_i_prev = NULL;

  if (posix_memalign((void **)&prev,       16, Q * sizeof(__m128i)) != 0) { free(pins); return eslEMEM; }
  if (posix_memalign((void **)&curr,       16, Q * sizeof(__m128i)) != 0) { free(prev); free(pins); return eslEMEM; }
  if (posix_memalign((void **)&peak,       16, Q * sizeof(__m128i)) != 0) { free(curr); free(prev); free(pins); return eslEMEM; }
  if (posix_memalign((void **)&peak_i,     16, Q * sizeof(__m128i)) != 0) { free(peak); free(curr); free(prev); free(pins); return eslEMEM; }
  if (posix_memalign((void **)&peak_prev,  16, Q * sizeof(__m128i)) != 0) { free(peak_i); free(peak); free(curr); free(prev); free(pins); return eslEMEM; }
  if (posix_memalign((void **)&peak_i_prev,16, Q * sizeof(__m128i)) != 0) { free(peak_prev); free(peak_i); free(peak); free(curr); free(prev); free(pins); return eslEMEM; }

  for (q = 0; q < Q; q++) {
    prev[q] = _mm_setzero_si128(); curr[q] = _mm_setzero_si128();
    peak[q] = _mm_setzero_si128(); peak_i[q] = _mm_setzero_si128();
    peak_prev[q] = _mm_setzero_si128(); peak_i_prev[q] = _mm_setzero_si128();
  }

  __m128i zero = _mm_setzero_si128();

  for (i = 1; i <= L; i++) {
    __m128i const *rsc = om32->riv[dsq[i]];               /* 32-bit emission scores */
    __m128i mpv = _mm_slli_si128(prev[Q - 1], 4);         /* shift left 4 bytes = 1 int32 lane */

    for (q = 0; q < Q; q++) {
      __m128i cand   = _mm_add_epi32(mpv, rsc[q]);         /* SW extend (NON-saturating) */
      __m128i newval = _mm_max_epi32(zero, cand);          /* SW floor at 0 (SSE4.1) */
      mpv     = prev[q];
      curr[q] = newval;
    }

    /* Peak tracking: snapshot previous peak, update with curr, reset on segment end */
    for (q = 0; q < Q; q++) {
      peak_prev[q]   = peak[q];
      peak_i_prev[q] = peak_i[q];
      __m128i is_new_peak = _mm_cmpgt_epi32(curr[q], peak[q]);
      peak[q]   = _mm_max_epi32(peak[q], curr[q]);
      peak_i[q] = _mm_blendv_epi8(peak_i[q], _mm_set1_epi32((int32_t)i), is_new_peak);
      __m128i is_reset = _mm_cmpeq_epi32(curr[q], zero);
      peak[q]   = _mm_andnot_si128(is_reset, peak[q]);
      peak_i[q] = _mm_andnot_si128(is_reset, peak_i[q]);
    }

    /* Collect end-of-segment pins */
    if (pin_mode == PIN_MODE_TRAIL) {
      union { __m128i v; int32_t b[4]; } u_prev, u_curr, u_pk;
      for (q = 0; q < Q; q++) {
        u_prev.v = prev[q];
        u_curr.v = curr[q];
        u_pk.v   = peak_prev[q];
        for (z = 0; z < 4; z++) {
          int k = (q + 1) + z * Q;
          if (k > M) break;
          if (u_prev.b[z] >= T_i32 && u_curr.b[z] == 0) {
            if (npins >= max_pins) {
              max_pins *= 2;
              PB_Pin *tmp = (PB_Pin *) realloc(pins, max_pins * sizeof(PB_Pin));
              if (!tmp) { free(pins); free(prev); free(curr); free(peak); free(peak_i); free(peak_prev); free(peak_i_prev); return eslEMEM; }
              pins = tmp;
            }
            pins[npins].k      = k;
            pins[npins].i      = i - 1;
            pins[npins].r      = u_prev.b[z];
            pins[npins].r_peak = u_pk.b[z];
            npins++;
          }
        }
      }
    } else { /* PIN_MODE_PEAK */
      union { __m128i v; int32_t b[4]; } u_pk, u_pki, u_curr;
      for (q = 0; q < Q; q++) {
        u_pk.v   = peak_prev[q];
        u_pki.v  = peak_i_prev[q];
        u_curr.v = curr[q];
        for (z = 0; z < 4; z++) {
          int k = (q + 1) + z * Q;
          if (k > M) break;
          if (u_pk.b[z] >= T_i32 && u_curr.b[z] == 0) {
            if (npins >= max_pins) {
              max_pins *= 2;
              PB_Pin *tmp = (PB_Pin *) realloc(pins, max_pins * sizeof(PB_Pin));
              if (!tmp) { free(pins); free(prev); free(curr); free(peak); free(peak_i); free(peak_prev); free(peak_i_prev); return eslEMEM; }
              pins = tmp;
            }
            pins[npins].k      = k;
            pins[npins].i      = u_pki.b[z];
            pins[npins].r      = u_pk.b[z];
            pins[npins].r_peak = u_pk.b[z];
            npins++;
          }
        }
      }
    }

    __m128i *tmp = prev; prev = curr; curr = tmp;
    for (q = 0; q < Q; q++) curr[q] = _mm_setzero_si128();
  }

  /* Last-row pins (trail mode only; peak mode uses end-of-scan flush below) */
  if (pin_mode == PIN_MODE_TRAIL) {
    union { __m128i v; int32_t b[4]; } u_prev, u_pk;
    for (q = 0; q < Q; q++) {
      u_prev.v = prev[q];
      u_pk.v   = peak[q];
      for (z = 0; z < 4; z++) {
        int k = (q + 1) + z * Q;
        if (k > M) break;
        if (u_prev.b[z] >= T_i32) {
          if (npins >= max_pins) {
            max_pins *= 2;
            PB_Pin *tmp2 = (PB_Pin *) realloc(pins, max_pins * sizeof(PB_Pin));
            if (!tmp2) { free(pins); free(prev); free(curr); free(peak); free(peak_i); free(peak_prev); free(peak_i_prev); return eslEMEM; }
            pins = tmp2;
          }
          pins[npins].k      = k;
          pins[npins].i      = L;
          pins[npins].r      = u_prev.b[z];
          pins[npins].r_peak = u_pk.b[z];
          npins++;
        }
      }
    }
  }

  /* End-of-scan flush (peak mode): emit for sticky segments still active at i=L */
  /* Brief 26_0430-102: flush is gated separately via PB_PIN_FLUSH to isolate from peak-on-reset. */
  if (pin_mode == PIN_MODE_PEAK && pb_get_flush_enabled()) {
    union { __m128i v; int32_t b[4]; } u_pk, u_pki, u_prev;
    for (q = 0; q < Q; q++) {
      u_pk.v   = peak[q];
      u_pki.v  = peak_i[q];
      u_prev.v = prev[q];    /* SW values at row L (after final swap) */
      for (z = 0; z < 4; z++) {
        int k = (q + 1) + z * Q;
        if (k > M) break;
        if (u_pk.b[z] >= T_i32 && u_prev.b[z] > 0) {
          if (npins >= max_pins) {
            max_pins *= 2;
            PB_Pin *tmp3 = (PB_Pin *) realloc(pins, max_pins * sizeof(PB_Pin));
            if (!tmp3) { free(pins); free(prev); free(curr); free(peak); free(peak_i); free(peak_prev); free(peak_i_prev); return eslEMEM; }
            pins = tmp3;
          }
          pins[npins].k      = k;
          pins[npins].i      = u_pki.b[z];
          pins[npins].r      = u_pk.b[z];
          pins[npins].r_peak = u_pk.b[z];
          npins++;
        }
      }
    }
  }

  free(prev); free(curr); free(peak); free(peak_i); free(peak_prev); free(peak_i_prev);

  *ret_pins  = pins;
  *ret_npins = npins;
  return eslOK;
}

/*****************************************************************
 * Gap-aware LSIS gap costs (brief 26_0430-089, 2026-05-29)
 *
 * Plain score-sum LSIS chains pins to maximize the sum of SW pin
 * scores, ignoring the indel cost implied by the jump between two
 * consecutive pins. On fragmented RF00010 targets this lets LSIS
 * splice together pins from incompatible diagonals across a long
 * D-state run that Viterbi then scores as implausible (the remaining
 * 3 catastrophic failures after the 16-bit scan port).
 *
 * PB_GapData precomputes per-position transition penalties from the
 * generic profile gm->tsc (nats, <=0) converted to pinbridge 16-bit
 * units, so a closed-form gap cost between any two pins is O(1).
 *
 * Unit note: pin scores r come from om->rwv = wordify(sc) =
 * roundf(om->scale_w * sc) where sc is a nat-valued log-odds score
 * (p7P_MSC). The correct nats->pinbridge factor is therefore
 * om->scale_w itself (scale_w = 500/log2 already folds in the
 * nats->bits conversion). [Brief 26_0430-089 suggested scale_w/log2; that
 * double-counts the log2 and was corrected here, verified against
 * hmmer impl_sse/p7_oprofile.c wordify().]
 *
 * Per-position transition costs are clamped to PB_GAP_UNIT_CAP so an
 * impossible (-inf) transition cannot overflow the int prefix sums.
 *****************************************************************/
#define PB_GAP_UNIT_CAP 60000   /* max PB units for one transition penalty */

/* Brief 26_0430-091: cap the number of pins handed to the O(N^2) gap-aware LSIS.
 * Pins above the cap are pruned in score-descending order before LSIS.
 * 250 was selected from a K-sweep on LSU full (M=3400) + RF00010-frags +
 * SSU full / frags200: identical accuracy on RF00010-frags and on the 6
 * rmark4e families, with the LSU-full LSIS-stage cost dropping ~64x and
 * the total CYK bit-score actually IMPROVING by ~+70 on one LSU seq.
 * PB_TOPK env var overrides (set 0 to disable pruning entirely). */
#define PB_DEFAULT_TOPK 250

typedef struct {
  int    M;
  int   *S_MM;      /* S_MM[k] = sum_{j=1..k} cost(M_j->M_{j+1}), PB units, k=0..M */
  int   *S_DD;      /* S_DD[k] = sum_{j=1..k} cost(D_j->D_{j+1}), PB units, k=0..M */
  int   *md_at;     /* md_at[k] = cost(M_k->D_{k+1}), PB units                     */
  int   *dm_at;     /* dm_at[k] = cost(D_k->M_{k+1}), PB units                     */
  int   *ii_at;     /* ii_at[k]   = cost(I_k->I_k) loop, PB units                  */
  int   *ient_at;   /* ient_at[k] = cost(M_k->I_k + I_k->M_{k+1}) enter+exit, PB    */
  float  scale_natsToPB;  /* nats -> pinbridge units (= om->scale_w)               */
} PB_GapData;

static inline int
pb_nat_to_units(float natcost, float scale)
{
  float v;
  if (natcost < 0.0f) natcost = 0.0f;   /* numerical guard; cost is -(log prob) >= 0 */
  v = natcost * scale;
  if (v > (float) PB_GAP_UNIT_CAP) v = (float) PB_GAP_UNIT_CAP;
  return (int) roundf(v);
}

static void
pb_gap_data_free(PB_GapData *gd)
{
  if (gd == NULL) return;
  if (gd->S_MM)     free(gd->S_MM);
  if (gd->S_DD)     free(gd->S_DD);
  if (gd->md_at)    free(gd->md_at);
  if (gd->dm_at)    free(gd->dm_at);
  if (gd->ii_at)    free(gd->ii_at);
  if (gd->ient_at)  free(gd->ient_at);
  gd->S_MM = gd->S_DD = gd->md_at = gd->dm_at = gd->ii_at = gd->ient_at = NULL;
}

/* Precompute per-profile gap-cost tables. O(M). Call pb_gap_data_free()
 * to release. gm supplies transition log-probs (nats); om supplies the
 * scale to pinbridge 16-bit units. */
static int
pb_precompute_gap_data(const P7_PROFILE *gm, const P7_OPROFILE *om, PB_GapData *gd)
{
  int M = gm->M;
  int k;
  int status;

  gd->M    = M;
  gd->S_MM = gd->S_DD = gd->md_at = gd->dm_at = gd->ii_at = gd->ient_at = NULL;
  gd->scale_natsToPB = om->scale_w;

  ESL_ALLOC(gd->S_MM,    sizeof(int) * (M + 1));
  ESL_ALLOC(gd->S_DD,    sizeof(int) * (M + 1));
  ESL_ALLOC(gd->md_at,   sizeof(int) * (M + 1));
  ESL_ALLOC(gd->dm_at,   sizeof(int) * (M + 1));
  ESL_ALLOC(gd->ii_at,   sizeof(int) * (M + 1));
  ESL_ALLOC(gd->ient_at, sizeof(int) * (M + 1));

  gd->S_MM[0] = gd->S_DD[0] = 0;
  gd->md_at[0] = gd->dm_at[0] = gd->ii_at[0] = gd->ient_at[0] = 0;

  /* Transitions are defined for k = 1..M-1 (tsc hand-indexed [1..M-1]).
   * Position M has no outgoing core transition -> zero cost there. */
  for (k = 1; k <= M; k++) {
    int mm_pb = 0, dd_pb = 0;
    if (k <= M - 1) {
      mm_pb         = pb_nat_to_units(-p7P_TSC(gm, k, p7P_MM), gd->scale_natsToPB);
      dd_pb         = pb_nat_to_units(-p7P_TSC(gm, k, p7P_DD), gd->scale_natsToPB);
      gd->md_at[k]  = pb_nat_to_units(-p7P_TSC(gm, k, p7P_MD), gd->scale_natsToPB);
      gd->dm_at[k]  = pb_nat_to_units(-p7P_TSC(gm, k, p7P_DM), gd->scale_natsToPB);
      gd->ii_at[k]  = pb_nat_to_units(-p7P_TSC(gm, k, p7P_II), gd->scale_natsToPB);
      gd->ient_at[k]= pb_nat_to_units(-(p7P_TSC(gm, k, p7P_MI) +
                                        p7P_TSC(gm, k, p7P_IM)), gd->scale_natsToPB);
    } else {
      gd->md_at[k] = gd->dm_at[k] = gd->ii_at[k] = gd->ient_at[k] = 0;
    }
    gd->S_MM[k] = gd->S_MM[k-1] + mm_pb;
    gd->S_DD[k] = gd->S_DD[k-1] + dd_pb;
  }
  return eslOK;

 ERROR:
  pb_gap_data_free(gd);
  return status;
}

/* Closed-form gap cost (Option 2): O(1) approximate indel penalty for
 * the jump from pin (i_p,k_p) to pin (i_q,k_q), in pinbridge units.
 * Returns INT_MAX if the gap is not strictly forward in both i and k. */
static int
pb_gap_cost_closed_form(int i_p, int k_p, int i_q, int k_q, const PB_GapData *gd)
{
  int di = i_q - i_p;
  int dk = k_q - k_p;
  int mm_cost;

  if (di <= 0 || dk <= 0) return INT_MAX;   /* LSIS should not propose these */

  /* baseline: match-to-match transitions spanned, k_p .. k_q-1 (small) */
  mm_cost = gd->S_MM[k_q - 1] - gd->S_MM[k_p - 1];

  if (dk == di) {
    /* pure diagonal: di matches, MM transitions only */
    return mm_cost;
  }
  else if (dk > di) {
    /* net deletions d = dk - di: enter MD once + (d-1) DD loops + exit DM once */
    int  d    = dk - di;
    int  span = dk - 1; if (span < 1) span = 1;
    long dd_sum = (long) gd->S_DD[k_q - 1] - (long) gd->S_DD[k_p - 1];
    int  avg_dd = (int)(dd_sum / span);
    long cost = (long) mm_cost
              + (long)(d - 1) * (long) avg_dd
              + (long) gd->md_at[k_p]
              + (long) gd->dm_at[k_q - 1];
    if (cost > (long) INT_MAX) cost = INT_MAX;
    return (int) cost;
  }
  else {
    /* net insertions n = di - dk: enter MI once + (n-1) II loops + exit IM once,
     * at the cheapest insert state in the span. Charging the MI/IM boundary per
     * inserted residue (as a naive n*(MI+II+IM) would) over-penalizes long
     * insertions by ~(n-1)*(MI+IM) and wrongly truncates legitimate chains. */
    int n = di - dk;
    int best_ii = INT_MAX, best_ient = INT_MAX, kk;
    for (kk = k_p; kk <= k_q && kk <= gd->M; kk++) {
      if (gd->ii_at[kk]   < best_ii)   best_ii   = gd->ii_at[kk];
      if (gd->ient_at[kk] < best_ient) best_ient = gd->ient_at[kk];
    }
    if (best_ii   == INT_MAX) best_ii   = 0;
    if (best_ient == INT_MAX) best_ient = 0;
    long cost = (long) mm_cost + (long) best_ient + (long)(n - 1) * (long) best_ii;
    if (cost > (long) INT_MAX) cost = INT_MAX;
    return (int) cost;
  }
}

/* Mini-Viterbi gap cost (Option 3, brief 26_0430-089, flag-gated). Exact banded
 * p7 Viterbi over the rectangle [i_p,k_p]..[i_q,k_q]: the best-scoring
 * M/I/D path from pin p's match to pin q's match, scoring intermediate
 * transitions (gm->tsc) and emissions (gm->rsc, gm->isc) in nats, then
 * scaled to pinbridge units. O((di+1)*(dk+1)) cells per call. Returns
 * the positive cost, or INT_MAX on an invalid gap / no path / alloc fail.
 *
 * Relative indices: row a in [0..di] = absolute residue i_p+a; col b in
 * [0..dk] = absolute model position k_p+b. (0,0) is seeded at M_{k_p}
 * (pin p, residue i_p already emitted); the path must consume residues
 * i_p+1..i_q and terminate in a match at M_{k_q} (pin q). */
static int
pb_gap_cost_miniviterbi(int i_p, int k_p, int i_q, int k_q,
                        const P7_PROFILE *gm, const ESL_DSQ *dsq,
                        float scale_natsToPB)
{
  int di = i_q - i_p;
  int dk = k_q - k_p;
  int M  = gm->M;
  if (di <= 0 || dk <= 0) return INT_MAX;

  int nrow = di + 1;
  int ncol = dk + 1;
  float *Mx = (float *) malloc(sizeof(float) * nrow * ncol);
  float *Ix = (float *) malloc(sizeof(float) * nrow * ncol);
  float *Dx = (float *) malloc(sizeof(float) * nrow * ncol);
  if (!Mx || !Ix || !Dx) { free(Mx); free(Ix); free(Dx); return INT_MAX; }

  #define PBMX(a,b) Mx[(a)*ncol + (b)]
  #define PBIX(a,b) Ix[(a)*ncol + (b)]
  #define PBDX(a,b) Dx[(a)*ncol + (b)]
  const float NEG = -eslINFINITY;

  for (int a = 0; a < nrow; a++)
    for (int b = 0; b < ncol; b++) { PBMX(a,b) = PBIX(a,b) = PBDX(a,b) = NEG; }
  PBMX(0,0) = 0.0f;   /* sitting at pin p's match */

  for (int a = 0; a < nrow; a++) {
    int i_abs = i_p + a;
    for (int b = 0; b < ncol; b++) {
      int k_abs = k_p + b;
      if (a == 0 && b == 0) continue;

      /* M[a][b]: enter match k_abs emitting residue i_abs, from (a-1,b-1) */
      if (a >= 1 && b >= 1) {
        int   kprev = k_abs - 1;   /* in [k_p, k_q-1] subset of [1,M-1] */
        float best = NEG, v;
        v = PBMX(a-1,b-1) + p7P_TSC(gm, kprev, p7P_MM); if (v > best) best = v;
        v = PBIX(a-1,b-1) + p7P_TSC(gm, kprev, p7P_IM); if (v > best) best = v;
        v = PBDX(a-1,b-1) + p7P_TSC(gm, kprev, p7P_DM); if (v > best) best = v;
        if (best > NEG) PBMX(a,b) = best + p7P_MSC(gm, k_abs, dsq[i_abs]);
      }

      /* I[a][b]: insert at k_abs emitting residue i_abs, from (a-1,b).
       * Insert states exist only for k = 1..M-1. */
      if (a >= 1 && k_abs <= M - 1) {
        float best = NEG, v;
        v = PBMX(a-1,b) + p7P_TSC(gm, k_abs, p7P_MI); if (v > best) best = v;
        v = PBIX(a-1,b) + p7P_TSC(gm, k_abs, p7P_II); if (v > best) best = v;
        if (best > NEG) PBIX(a,b) = best + p7P_ISC(gm, k_abs, dsq[i_abs]);
      }

      /* D[a][b]: delete k_abs, no emission, from (a,b-1) */
      if (b >= 1) {
        int   kprev = k_abs - 1;
        float best = NEG, v;
        v = PBMX(a,b-1) + p7P_TSC(gm, kprev, p7P_MD); if (v > best) best = v;
        v = PBDX(a,b-1) + p7P_TSC(gm, kprev, p7P_DD); if (v > best) best = v;
        if (best > NEG) PBDX(a,b) = best;
      }
    }
  }

  float endsc = PBMX(nrow-1, ncol-1);   /* must end in match at pin q */
  int cost;
  if (endsc <= NEG) {
    cost = INT_MAX;
  } else {
    float c = -endsc * scale_natsToPB;  /* nat path score (<=0) -> positive PB cost */
    if (c < 0.0f) c = 0.0f;
    if (c > (float) INT_MAX) c = (float) INT_MAX;
    cost = (int) roundf(c);
  }

  #undef PBMX
  #undef PBIX
  #undef PBDX
  free(Mx); free(Ix); free(Dx);
  return cost;
}

/* LSIS over (i,k,r): longest score-weighted increasing subsequence.
 * Sort by (i asc, k desc) so same-i pins never chain.
 * Fenwick prefix-max over k stores (max_dp, achiever_idx). */
static int pb_pin_cmp_lsis(const void *a, const void *b)
{
  const PB_Pin *pa = (const PB_Pin *)a;
  const PB_Pin *pb = (const PB_Pin *)b;
  if (pa->i != pb->i) return pa->i - pb->i;
  return pb->k - pa->k;
}

/* (k asc, i asc) order for the limited-window LSIS (brief 26_0430-094). A valid chain
 * edge q->p requires k_q < k_p, so with k ascending every predecessor is
 * processed before p, and predecessors within a model-position window of p
 * form a contiguous block ending just before p -- enabling an early break. */
static int pb_pin_cmp_lsis_kasc(const void *a, const void *b)
{
  const PB_Pin *pa = (const PB_Pin *)a;
  const PB_Pin *pb = (const PB_Pin *)b;
  if (pa->k != pb->k) return pa->k - pb->k;
  return pa->i - pb->i;
}

/* Sort pins by score r descending. Used by Top-K pruning. */
static int pb_pin_cmp_r_desc(const void *a, const void *b)
{
  const PB_Pin *pa = (const PB_Pin *)a;
  const PB_Pin *pb = (const PB_Pin *)b;
  if (pa->r != pb->r) return (int)pb->r - (int)pa->r;
  /* Tie-break by (i, k) for deterministic ordering. */
  if (pa->i != pb->i) return pa->i - pb->i;
  return pa->k - pb->k;
}

/* Attack 1 from brief 26_0430-091: cap pin count at K by keeping top-K by score r.
 * In-place: reorders <pins> and returns the new count via <ret_n_out>.
 * If K <= 0 or npins <= K, no-op.
 * O(N log N) (qsort); negligible vs. the O(N^2) LSIS it feeds. */
static int
pb_prune_pins_topk(PB_Pin *pins, int npins, int K, int *ret_n_out)
{
  if (K <= 0 || npins <= K) { *ret_n_out = npins; return eslOK; }
  qsort(pins, npins, sizeof(PB_Pin), pb_pin_cmp_r_desc);
  *ret_n_out = K;
  return eslOK;
}

/* Legacy gap-blind LSIS (Fenwick prefix-max). Retained for comparison
 * and as a fallback; superseded by pb_lsis_select_gap_aware (brief 26_0430-089). */
static int
pb_lsis_select_legacy(PB_Pin *pins_in, int npins_in, int M,
               PB_Pin **ret_sel, int *ret_n_sel)
{
  *ret_sel = NULL;
  *ret_n_sel = 0;
  if (npins_in == 0) return eslOK;

  PB_Pin *pins = (PB_Pin *) malloc(npins_in * sizeof(PB_Pin));
  if (!pins) return eslEMEM;
  memcpy(pins, pins_in, npins_in * sizeof(PB_Pin));
  qsort(pins, npins_in, sizeof(PB_Pin), pb_pin_cmp_lsis);

  int *dp     = (int *) malloc(npins_in * sizeof(int));
  int *parent = (int *) malloc(npins_in * sizeof(int));
  if (!dp || !parent) { free(pins); free(dp); free(parent); return eslEMEM; }

  int  fen_size = M + 1;
  int *fen_val  = (int *) calloc(fen_size + 1, sizeof(int));
  int *fen_idx  = (int *) malloc((fen_size + 1) * sizeof(int));
  if (!fen_val || !fen_idx) { free(pins); free(dp); free(parent); free(fen_val); free(fen_idx); return eslEMEM; }
  for (int x = 0; x <= fen_size; x++) fen_idx[x] = -1;

  int best_dp = -1, best_idx = -1;
  for (int p = 0; p < npins_in; p++) {
    int kp = pins[p].k;
    int sp = (int) pins[p].r;
    int q_val = 0, q_idx = -1;
    for (int x = kp - 1; x > 0; x -= x & -x) {
      if (fen_val[x] > q_val) { q_val = fen_val[x]; q_idx = fen_idx[x]; }
    }
    dp[p]     = q_val + sp;
    parent[p] = q_idx;
    for (int x = kp; x <= fen_size; x += x & -x) {
      if (dp[p] > fen_val[x]) { fen_val[x] = dp[p]; fen_idx[x] = p; }
    }
    if (dp[p] > best_dp) { best_dp = dp[p]; best_idx = p; }
  }

  int *chain_idx = (int *) malloc(npins_in * sizeof(int));
  if (!chain_idx) { free(pins); free(dp); free(parent); free(fen_val); free(fen_idx); return eslEMEM; }
  int chain_len = 0;
  for (int cur = best_idx; cur >= 0; cur = parent[cur]) chain_idx[chain_len++] = cur;

  PB_Pin *sel = (PB_Pin *) malloc(chain_len * sizeof(PB_Pin));
  if (!sel) { free(pins); free(dp); free(parent); free(fen_val); free(fen_idx); free(chain_idx); return eslEMEM; }
  for (int j = 0; j < chain_len; j++) sel[j] = pins[chain_idx[chain_len - 1 - j]];

  free(pins); free(dp); free(parent); free(fen_val); free(fen_idx); free(chain_idx);

  *ret_sel = sel;
  *ret_n_sel = chain_len;
  return eslOK;
}

/* Gap-aware LSIS (brief 26_0430-089, Option 2). O(N^2) DP that prices the indel
 * cost of each pin->pin jump via pb_gap_cost_closed_form(), so the
 * objective is sum(pin r) - sum(gap cost) rather than sum(pin r) alone.
 * This stops LSIS from splicing pins across incompatible diagonals over
 * an implausible D-state run. Pins are sorted (i asc, k desc) as in the
 * legacy path; a chain edge q->p requires strictly increasing i and k.
 *
 * When use_vit_gaps != 0 (Option 3, --p7pinbridge-vitgaps), each edge's
 * gap cost is the exact mini-Viterbi over the gap rectangle instead of
 * the closed-form estimate; gm/dsq are then required (NULL otherwise).
 *
 * Brief 26_0430-094: when lsis_window > 0 the DP only considers predecessors within
 * <lsis_window> model positions of each pin, turning the O(N^2) inner loop
 * into O(N * window). This requires the (k asc, i asc) ordering so the
 * window is a contiguous, early-breakable block. When lsis_window <= 0 the
 * original (i asc, k desc) sort and full ascending scan are used verbatim,
 * so the result is byte-identical to the brief-091 path (gate-off A/B). */
static int
pb_lsis_select_gap_aware(PB_Pin *pins_in, int npins_in, int M,
                         const PB_GapData *gd,
                         int use_vit_gaps, const P7_PROFILE *gm, const ESL_DSQ *dsq,
                         int lsis_window,
                         PB_Pin **ret_sel, int *ret_n_sel)
{
  int  status;
  int *dp        = NULL;
  int *parent    = NULL;
  int *chain_idx = NULL;
  PB_Pin *pins   = NULL;
  PB_Pin *sel    = NULL;

  *ret_sel = NULL;
  *ret_n_sel = 0;
  if (npins_in == 0) return eslOK;

  ESL_ALLOC(pins,   npins_in * sizeof(PB_Pin));
  memcpy(pins, pins_in, npins_in * sizeof(PB_Pin));
  qsort(pins, npins_in, sizeof(PB_Pin),
        (lsis_window > 0) ? pb_pin_cmp_lsis_kasc : pb_pin_cmp_lsis);

  ESL_ALLOC(dp,     npins_in * sizeof(int));
  ESL_ALLOC(parent, npins_in * sizeof(int));

  long n_gapcost = 0;   /* # of gap costs evaluated (PB_DEBUG_LSIS) */
  int best_dp = INT_MIN, best_idx = -1;
  if (lsis_window <= 0) {
    /* Legacy full O(N^2) DP (brief 26_0430-091), ascending q -> byte-identical. */
    for (int p = 0; p < npins_in; p++) {
      dp[p]     = (int) pins[p].r;   /* chain starting fresh at this pin */
      parent[p] = -1;
      for (int q = 0; q < p; q++) {
        /* forward edge requires strictly increasing i and k */
        if (pins[q].i >= pins[p].i || pins[q].k >= pins[p].k) continue;
        n_gapcost++;
        int gap;
        if (use_vit_gaps)
          gap = pb_gap_cost_miniviterbi(pins[q].i, pins[q].k,
                                        pins[p].i, pins[p].k,
                                        gm, dsq, gd->scale_natsToPB);
        else
          gap = pb_gap_cost_closed_form(pins[q].i, pins[q].k,
                                        pins[p].i, pins[p].k, gd);
        if (gap == INT_MAX) continue;
        int cand = dp[q] + (int) pins[p].r - gap;
        if (cand > dp[p]) { dp[p] = cand; parent[p] = q; }
      }
      if (dp[p] > best_dp) { best_dp = dp[p]; best_idx = p; }
    }
  } else {
    /* Limited-window DP (brief 26_0430-094): (k asc, i asc) sort; scan predecessors
     * backward and break once the model-position gap exceeds the window. */
    for (int p = 0; p < npins_in; p++) {
      dp[p]     = (int) pins[p].r;
      parent[p] = -1;
      for (int q = p - 1; q >= 0; q--) {
        if ((pins[p].k - pins[q].k) > lsis_window) break;   /* k is monotone-ascending */
        /* forward edge requires strictly increasing i and k */
        if (pins[q].i >= pins[p].i || pins[q].k >= pins[p].k) continue;
        n_gapcost++;
        int gap;
        if (use_vit_gaps)
          gap = pb_gap_cost_miniviterbi(pins[q].i, pins[q].k,
                                        pins[p].i, pins[p].k,
                                        gm, dsq, gd->scale_natsToPB);
        else
          gap = pb_gap_cost_closed_form(pins[q].i, pins[q].k,
                                        pins[p].i, pins[p].k, gd);
        if (gap == INT_MAX) continue;
        int cand = dp[q] + (int) pins[p].r - gap;
        if (cand > dp[p]) { dp[p] = cand; parent[p] = q; }
      }
      if (dp[p] > best_dp) { best_dp = dp[p]; best_idx = p; }
    }
  }

  ESL_ALLOC(chain_idx, npins_in * sizeof(int));
  int chain_len = 0;
  for (int cur = best_idx; cur >= 0; cur = parent[cur]) chain_idx[chain_len++] = cur;

  ESL_ALLOC(sel, chain_len * sizeof(PB_Pin));
  for (int j = 0; j < chain_len; j++) sel[j] = pins[chain_idx[chain_len - 1 - j]];

  if (getenv("PB_DEBUG_LSIS") != NULL)
    fprintf(stderr, "#PB_LSIS_COST npins=%d n_gapcost=%ld vit_gaps=%d\n",
            npins_in, n_gapcost, use_vit_gaps);

  free(pins); free(dp); free(parent); free(chain_idx);
  *ret_sel = sel;
  *ret_n_sel = chain_len;
  return eslOK;

 ERROR:
  if (pins)      free(pins);
  if (dp)        free(dp);
  if (parent)    free(parent);
  if (chain_idx) free(chain_idx);
  if (sel)       free(sel);
  return status;
}


/* Build per-row kmin/kmax from selected pins: pin diagonal segments + bridge
 * rectangles + ±pad. Mirrors prototype's build_band_pinbridge() steps 4a-4e. */
static int
pb_build_band(PB_Pin *kpins, int nfinal, int L, int M, int pad,
              int **ret_kmin, int **ret_kmax)
{
  int *kmin = NULL, *kmax = NULL;
  int  i;
  int  status;

  ESL_ALLOC(kmin, sizeof(int) * (L + 2));
  ESL_ALLOC(kmax, sizeof(int) * (L + 2));
  for (i = 0; i <= L+1; i++) { kmin[i] = M + 1; kmax[i] = 0; }

  if (nfinal == 0) {
    for (i = 0; i <= L; i++) { kmin[i] = 1; kmax[i] = M; }
    *ret_kmin = kmin; *ret_kmax = kmax;
    return eslOK;
  }

  /* Preamble: i in [0..first_pin_diag_start - 1] -> [1, kpins[0].k] */
  {
    int first_diag_start = kpins[0].i - (int)kpins[0].r + 1;
    if (first_diag_start < 1) first_diag_start = 1;
    for (i = 0; i < first_diag_start; i++) {
      if (1            < kmin[i]) kmin[i] = 1;
      if (kpins[0].k   > kmax[i]) kmax[i] = kpins[0].k;
    }
  }

  /* Pin diagonal segments with ±pad */
  for (int p = 0; p < nfinal; p++) {
    int k_p = kpins[p].k;
    int i_p = kpins[p].i;
    int r   = (int)kpins[p].r;
    int j_start = i_p - r + 1; if (j_start < 1) j_start = 1;
    for (int j = j_start; j <= i_p; j++) {
      int k_diag = k_p - (i_p - j);
      int k_lo = k_diag - pad; if (k_lo < 1) k_lo = 1; if (k_lo > M) k_lo = M;
      int k_hi = k_diag + pad; if (k_hi > M) k_hi = M; if (k_hi < 1) k_hi = 1;
      if (k_lo < kmin[j]) kmin[j] = k_lo;
      if (k_hi > kmax[j]) kmax[j] = k_hi;
    }
  }

  /* Bridge rectangles between consecutive pins */
  for (int p = 0; p + 1 < nfinal; p++) {
    int k1 = kpins[p].k,     i1 = kpins[p].i;
    int k2 = kpins[p+1].k,   i2 = kpins[p+1].i;
    if (i2 <= i1 + 1) continue;
    int klo = (k1 < k2) ? k1 : k2;
    int khi = (k1 > k2) ? k1 : k2;
    for (int j = i1 + 1; j <= i2 - 1; j++) {
      if (klo < kmin[j]) kmin[j] = klo;
      if (khi > kmax[j]) kmax[j] = khi;
    }
  }

  /* Coda */
  {
    int last_k = kpins[nfinal-1].k;
    int last_i = kpins[nfinal-1].i;
    for (i = last_i + 1; i <= L; i++) {
      if (last_k < kmin[i]) kmin[i] = last_k;
      if (M      > kmax[i]) kmax[i] = M;
    }
  }

  /* Safety: any uncovered row -> [1, M] */
  for (i = 0; i <= L; i++) {
    if (kmin[i] > kmax[i]) { kmin[i] = 1; kmax[i] = M; }
  }

  *ret_kmin = kmin;
  *ret_kmax = kmax;
  return eslOK;

 ERROR:
  if (kmin) free(kmin);
  if (kmax) free(kmax);
  return eslEMEM;
}


/* p7_GBands_FromKminKmax():
 * Convert per-row kmin[1..L]/kmax[1..L] arrays into a single-segment
 * P7_GBANDS (no holes; rows 1..L are all in band).
 * The bnd struct is reused (caller may pass a freshly-Created or Reuse'd one).
 */
int
p7_GBands_FromKminKmax(int *kmin, int *kmax, int L, int M, P7_GBANDS *bnd)
{
  int i;
  int status;

  bnd->L = L;
  bnd->M = M;
  if ((status = p7_gbands_Reuse(bnd)) != eslOK) return status;

  /* Append rows 1..L as a single segment: p7_gbands_Append handles
   * segment opening/closing automatically when consecutive i's are
   * appended. */
  for (i = 1; i <= L; i++) {
    int ka = kmin[i], kb = kmax[i];
    if (ka < 1)   ka = 1;
    if (kb > M)   kb = M;
    if (ka > kb)  { ka = 1; kb = M; }   /* safety */
    if ((status = p7_gbands_Append(bnd, i, ka, kb)) != eslOK) return status;
  }
  return eslOK;
}


/* pb_load_pins_from_tsv(): Load a pin set from a TSV file (brief 26_0430-106 PB_LOAD_PINS hook).
 * Format: whitespace-separated "i k r" per line; lines starting with '#' are skipped;
 * lines that don't parse as three integers (e.g. headers) are skipped silently.
 * Validates 1 <= i <= L and 1 <= k <= M; aborts with eslFAIL on out-of-bound pins.
 * Caller is responsible for free()ing *ret_pins on eslOK return.
 */
static int
pb_load_pins_from_tsv(const char *path, int M, int L, PB_Pin **ret_pins, int *ret_npins)
{
  FILE   *fp     = NULL;
  PB_Pin *pins   = NULL;
  int     nalloc = 256;
  int     npins  = 0;
  char    line[512];

  *ret_pins  = NULL;
  *ret_npins = 0;

  if ((fp = fopen(path, "r")) == NULL) {
    fprintf(stderr, "PB_LOAD_PINS: cannot open '%s'\n", path);
    return eslFAIL;
  }
  if ((pins = (PB_Pin *) malloc(nalloc * sizeof(PB_Pin))) == NULL) {
    fclose(fp); return eslEMEM;
  }

  while (fgets(line, sizeof(line), fp) != NULL) {
    int i, k, r;
    if (line[0] == '#') continue;
    if (sscanf(line, "%d %d %d", &i, &k, &r) != 3) continue;
    if (i < 1 || i > L) {
      fprintf(stderr, "PB_LOAD_PINS: pin i=%d out of bounds [1,%d]\n", i, L);
      fclose(fp); free(pins); return eslFAIL;
    }
    if (k < 1 || k > M) {
      fprintf(stderr, "PB_LOAD_PINS: pin k=%d out of bounds [1,%d]\n", k, M);
      fclose(fp); free(pins); return eslFAIL;
    }
    if (npins == nalloc) {
      PB_Pin *tmp;
      nalloc *= 2;
      if ((tmp = (PB_Pin *) realloc(pins, nalloc * sizeof(PB_Pin))) == NULL) {
        fclose(fp); free(pins); return eslEMEM;
      }
      pins = tmp;
    }
    pins[npins].i      = i;
    pins[npins].k      = k;
    pins[npins].r      = (int32_t) r;
    pins[npins].r_peak = (int32_t) r;
    npins++;
  }
  fclose(fp);
  *ret_pins  = pins;
  *ret_npins = npins;
  return eslOK;
}


/* p7_Seq2BandsPinBridge():
 * Drop-in replacement for Steps 1-2 of p7_Seq2BandsVit():
 * SW-pinbridge prefilter + banded p7 Viterbi + banded trace.
 * Steps 3-4 (i2k extraction + p7_pins2bands_nodepad) are identical
 * and live in the caller (cm_alndata.c), unchanged.
 *
 * Args:
 *   gm        - profile (configured GLOCAL or LOCAL by caller)
 *   om        - optimized profile (for SSE rbv match scores, LOCAL config)
 *   om32      - 32-bit emission table for the 32-bit SW scan (brief 26_0430-094);
 *               may be NULL (then the 16-bit scan is always used)
 *   gxb       - banded Viterbi DP matrix scratch (will be Reinit'd to bnd)
 *   bnd       - GBANDS scratch (will be reused)
 *   tr        - p7 trace (will be Reuse'd; filled with banded-Vit traceback)
 *   dsq       - digital sequence, 1..L
 *   L         - sequence length
 *   pad       - diagonal pad for pinbridge bands (default 20 from caller)
 *   ret_sc    - optRETURN: banded Viterbi score
 *
 * Returns: eslOK on success; tr filled. eslFAIL if no valid trace
 *          (caller falls back to unbanded path).
 */
int
p7_Seq2BandsPinBridge(P7_PROFILE *gm, P7_OPROFILE *om, const CM_PB_OM32 *om32, P7_GMXB *gxb,
                      P7_GBANDS *bnd, P7_TRACE *tr,
                      const ESL_DSQ *dsq, int L, int pad, int use_vit_gaps, float *ret_sc,
                      double *ret_sw_ms, double *ret_lsis_ms,
                      double *ret_band_ms, double *ret_bvit_ms, double *ret_btrace_ms)
{
  int       status;
  int       M         = gm->M;
  PB_Pin   *raw_pins  = NULL;
  int       npins     = 0;
  PB_Pin   *sel_pins  = NULL;
  int       nsel      = 0;
  int      *kmin      = NULL;
  int      *kmax      = NULL;
  float     vit_sc    = -eslINFINITY;
  int       T         = pb_adaptive_T(M);
  /* 16-bit threshold: scale 8-bit T by (scale_w/scale_b) ~ 167.
   * Computed from om's actual scales for safety. */
  int       T_w       = (int)( (float)T * (om->scale_w / om->scale_b) );
  PB_GapData gd;
  int       gd_ok     = 0;
  struct timespec ta, tb;
  double    sw_ms = 0, lsis_ms = 0, band_ms = 0, bvit_ms = 0, btrace_ms = 0;

  /* Brief 26_0430-094 component gates: ON by default only for M >= PB_M32_MIN (the
   * regime where 16-bit saturation blows up the band); an explicitly-set env
   * var overrides the M-gate at any M. */
  int       large_m     = (M >= pb_env_int("PB_M32_MIN", PB_M32_MIN_DEFAULT));
  int       use_32bit   = pb_gate("PB_USE_32BIT_SW", PB_DEFAULT_USE_32BIT_SW, large_m) && (om32 != NULL);
  int       lsis_window = pb_gate("PB_LSIS_WINDOW",  PB_DEFAULT_LSIS_WINDOW,  large_m);
  int       k_adaptive  = pb_gate("PB_K_ADAPTIVE",   PB_DEFAULT_K_ADAPTIVE,   large_m);
  /* 32-bit threshold: same semantics as T_w, rescaled by scale_i/scale_w.
   * scale_i == scale_w by construction, so T_i32 == T_w. */
  int       T_i32       = (om32 != NULL) ? (int)( (float)T_w * (om32->scale_i / om->scale_w) ) : T_w;

  /* Step 1: SW scan -> raw pins. 32-bit (brief 26_0430-094) eliminates the 16-bit
   * accumulation saturation on conserved/self-align inputs; the legacy
   * 16-bit SSE scan (brief 26_0430-079 S7) is kept as the gate-off fallback. */
  clock_gettime(CLOCK_MONOTONIC, &ta);
  if (use_32bit) {
    if ((status = pb_sw_scan_collect_pins_i32(dsq, L, om32, T_i32, &raw_pins, &npins)) != eslOK) goto ERROR;
  } else {
    if ((status = pb_sw_scan_collect_pins_w(dsq, L, om, T_w, &raw_pins, &npins)) != eslOK) goto ERROR;
  }
  clock_gettime(CLOCK_MONOTONIC, &tb);
  sw_ms = (tb.tv_sec - ta.tv_sec)*1000.0 + (tb.tv_nsec - ta.tv_nsec)/1e6;

  /* PB_DUMP_PINS diagnostic (brief 26_0430-099/101): append one TSV line per raw pin */
  {
    const char *dump_path = getenv("PB_DUMP_PINS");
    if (dump_path != NULL) {
      FILE *fp = fopen(dump_path, "a");
      if (fp != NULL) {
        const char *tag    = getenv("PB_DUMP_TAG");
        const char *kernel = use_32bit ? "i32" : "w16";
        int _p;
        if (tag == NULL) tag = "-";
        for (_p = 0; _p < npins; _p++)
          fprintf(fp, "%s\t%s\tM=%d\tL=%d\tk=%d\ti=%d\tr=%d\tr_peak=%d\n",
                  tag, kernel, M, L,
                  raw_pins[_p].k, raw_pins[_p].i,
                  (int)raw_pins[_p].r, (int)raw_pins[_p].r_peak);
        fclose(fp);
      }
    }
  }

  /* PB_LOAD_PINS hook (brief 26_0430-106): replace SW kernel output with an offline-generated
   * pin set (TSV with "i k r" columns) so all downstream steps (LSIS, band-build,
   * CP9 F/B, CM DP) are identical across trail/peak/greedy variants. */
  {
    const char *load_path = getenv("PB_LOAD_PINS");
    if (load_path != NULL) {
      if (raw_pins) { free(raw_pins); raw_pins = NULL; }
      npins = 0;
      if ((status = pb_load_pins_from_tsv(load_path, M, L, &raw_pins, &npins)) != eslOK) goto ERROR;
    }
  }

  /* Step 1.5 (brief 26_0430-091, Attack 1): Top-K pruning to cap LSIS input size.
   * Default K = PB_DEFAULT_TOPK; PB_TOPK env var overrides (0 = off).
   * Brief 26_0430-094: when PB_K_ADAPTIVE is on, K floors at PB_DEFAULT_TOPK and
   * scales as L/(2*pad) (the band floor reached when K >= L/(2*pad)), so
   * large inputs get enough pins to span the full sequence rather than
   * concentrating in a tie-break window. PB_TOPK still overrides everything.
   * Cost folded into the LSIS timer below (it precedes the pure DP). */
  int       npins_lsis_in = npins;
  {
    const char *s = getenv("PB_TOPK");
    int K;
    if (s != NULL)              K = atoi(s);                                   /* explicit override */
    else if (k_adaptive && pad > 0) K = ESL_MAX(PB_DEFAULT_TOPK, L / (2 * pad));
    else                        K = PB_DEFAULT_TOPK;
    if (K > 0) (void) pb_prune_pins_topk(raw_pins, npins, K, &npins_lsis_in);
  }

  /* Step 2: gap-aware LSIS pin selection (brief 26_0430-089, Option 2).
   * Precompute model-aware gap-cost tables once, then run the DP (full
   * O(N^2), or O(N*window) when lsis_window > 0, brief 26_0430-094). */
  clock_gettime(CLOCK_MONOTONIC, &ta);
  if ((status = pb_precompute_gap_data(gm, om, &gd)) != eslOK) goto ERROR;
  gd_ok = 1;
  if ((status = pb_lsis_select_gap_aware(raw_pins, npins_lsis_in, M, &gd,
                                         use_vit_gaps, gm, dsq, lsis_window,
                                         &sel_pins, &nsel)) != eslOK) goto ERROR;
  clock_gettime(CLOCK_MONOTONIC, &tb);
  lsis_ms = (tb.tv_sec - ta.tv_sec)*1000.0 + (tb.tv_nsec - ta.tv_nsec)/1e6;

  /* Per-seq pin-count diagnostic (brief 26_0430-091).
   * npins_raw = pins from SW scan; npins_lsis_in = pins handed to LSIS (after
   * any pruning); nsel = chain length selected by LSIS. */
  fprintf(stderr, "#PB_NPINS M=%d L=%d npins_raw=%d npins_lsis_in=%d nsel=%d\n",
          M, L, npins, npins_lsis_in, nsel);
  fflush(stderr);
  /* Brief 26_0430-104 diagnostic: dump LSIS-selected pins to file.
   * PB_DUMP_LSIS_PINS=<path> writes TSV with header pin_order, i, k, r.
   * DIAGNOSTIC-ONLY: no effect on band building. */
  {
    const char *lsis_dump_path = getenv("PB_DUMP_LSIS_PINS");
    if (lsis_dump_path != NULL) {
      FILE *lsis_fp = fopen(lsis_dump_path, "w");
      if (lsis_fp != NULL) {
        int _lp;
        fprintf(lsis_fp, "pin_order\ti\tk\tr\n");
        for (_lp = 0; _lp < nsel; _lp++)
          fprintf(lsis_fp, "%d\t%d\t%d\t%d\n",
                  _lp, sel_pins[_lp].i, sel_pins[_lp].k, (int)sel_pins[_lp].r);
        fclose(lsis_fp);
      }
    }
  }

  /* Optional: dump LSIS-selected chain when PB_DEBUG_LSIS is set in env.
   * Used for gap-aware LSIS validation (brief 26_0430-089). */
  if (getenv("PB_DEBUG_LSIS") != NULL) {
    int _dbg_p;
    long long _dbg_sum = 0;
    for (_dbg_p = 0; _dbg_p < nsel; _dbg_p++) _dbg_sum += sel_pins[_dbg_p].r;
    fprintf(stderr, "#PB_LSIS_OUT nsel=%d raw=%d sum_r=%lld first_k=%d last_k=%d\n",
            nsel, npins, _dbg_sum,
            nsel > 0 ? sel_pins[0].k : -1,
            nsel > 0 ? sel_pins[nsel-1].k : -1);
    for (_dbg_p = 0; _dbg_p < nsel; _dbg_p++) {
      fprintf(stderr, "#PB_LSIS_PIN p=%d i=%d k=%d r=%d\n",
              _dbg_p, sel_pins[_dbg_p].i, sel_pins[_dbg_p].k, (int)sel_pins[_dbg_p].r);
    }
  }

  /* Step 3+4: build kmin/kmax + convert to GBANDS + (allocate or grow) gxb */
  clock_gettime(CLOCK_MONOTONIC, &ta);
  if ((status = pb_build_band(sel_pins, nsel, L, M, pad, &kmin, &kmax)) != eslOK) goto ERROR;
  /* Always widen boundary rows: a local profile can enter/exit at any k;
   * a glocal profile must enter at M_1 and exit at M_M, but the existing
   * preamble/coda already cover those, and widening costs only a handful
   * of cells on i=0,1,L. Trace can fail without this widening when the
   * pinbridge band misses the k values the optimal trace needs at i=1 / i=L. */
  /* Always widen boundary rows i=0,1,L to [1,M]. Even glocal profiles need
   * room at the boundaries: row L's M_M depends on row L-1's k=M-1, which
   * pinbridge's coda may not reach. Cost: a few extra cells per sequence;
   * the interior bandwidth (which dominates the total) is unaffected. */
  { int bi;
    for (bi = 0; bi <= 1 && bi <= L; bi++) { kmin[bi] = 1; kmax[bi] = M; }
    if (L >= 1) { kmin[L] = 1; kmax[L] = M; }
    if (L >= 2) { kmin[L-1] = 1; kmax[L-1] = M; }  /* k=M-1 reachable for M_M */
  }
  if ((status = p7_GBands_FromKminKmax(kmin, kmax, L, M, bnd)) != eslOK) goto ERROR;
  /* Caller may pass gxb with NULL dp (deferred allocation; p7_gmxb_Create
   * disallows zero-cell bnd, so the wrapper allocates the struct shell and
   * we fill dp/xmx here once bnd has cells). */
  if (gxb->dp == NULL) {
    if ((gxb->dp  = malloc(sizeof(float) * bnd->ncell * p7G_NSCELLS)) == NULL) { status = eslEMEM; goto ERROR; }
    if ((gxb->xmx = malloc(sizeof(float) * bnd->nrow  * p7G_NXCELLS)) == NULL) { status = eslEMEM; goto ERROR; }
    gxb->dalloc = bnd->ncell;
    gxb->xalloc = bnd->nrow;
    gxb->bnd    = bnd;
  } else {
    if ((status = p7_gmxb_Reinit(gxb, bnd)) != eslOK) goto ERROR;
  }
  clock_gettime(CLOCK_MONOTONIC, &tb);
  band_ms = (tb.tv_sec - ta.tv_sec)*1000.0 + (tb.tv_nsec - ta.tv_nsec)/1e6;

  /* Band-size diagnostic (brief 26_0430-094 validation; env-gated, uncommitted). */
  if (getenv("PB_BANDVIT_INPUT") != NULL)
    fprintf(stderr, "#PB_BANDVIT_INPUT L=%d M=%d bnd_nrow=%d bnd_ncell=%lld bnd_nseg=%d\n",
            L, M, bnd->nrow, (long long)bnd->ncell, bnd->nseg);

  /* Step 5: Banded p7 Viterbi inside the prefilter band */
  clock_gettime(CLOCK_MONOTONIC, &ta);
  if ((status = p7_GBandedViterbi(dsq, L, gm, gxb, &vit_sc)) != eslOK) goto ERROR;
  clock_gettime(CLOCK_MONOTONIC, &tb);
  bvit_ms = (tb.tv_sec - ta.tv_sec)*1000.0 + (tb.tv_nsec - ta.tv_nsec)/1e6;

  /* Step 6: Banded traceback */
  p7_trace_Reuse(tr);
  clock_gettime(CLOCK_MONOTONIC, &ta);
  status = p7_GBandedTrace(dsq, L, gm, gxb, tr);
  clock_gettime(CLOCK_MONOTONIC, &tb);
  btrace_ms = (tb.tv_sec - ta.tv_sec)*1000.0 + (tb.tv_nsec - ta.tv_nsec)/1e6;
  if (status != eslOK) {
    /* Trace failed inside band; signal caller to fall back. */
    if (raw_pins) free(raw_pins);
    if (sel_pins) free(sel_pins);
    if (kmin)     free(kmin);
    if (kmax)     free(kmax);
    if (gd_ok)    pb_gap_data_free(&gd);
    if (ret_sc) *ret_sc = vit_sc;
    if (ret_sw_ms)     *ret_sw_ms     = sw_ms;
    if (ret_lsis_ms)   *ret_lsis_ms   = lsis_ms;
    if (ret_band_ms)   *ret_band_ms   = band_ms;
    if (ret_bvit_ms)   *ret_bvit_ms   = bvit_ms;
    if (ret_btrace_ms) *ret_btrace_ms = btrace_ms;
    return eslFAIL;
  }

  if (ret_sc) *ret_sc = vit_sc;
  if (ret_sw_ms)     *ret_sw_ms     = sw_ms;
  if (ret_lsis_ms)   *ret_lsis_ms   = lsis_ms;
  if (ret_band_ms)   *ret_band_ms   = band_ms;
  if (ret_bvit_ms)   *ret_bvit_ms   = bvit_ms;
  if (ret_btrace_ms) *ret_btrace_ms = btrace_ms;

  if (raw_pins) free(raw_pins);
  if (sel_pins) free(sel_pins);
  if (kmin)     free(kmin);
  if (kmax)     free(kmax);
  if (gd_ok)    pb_gap_data_free(&gd);
  return eslOK;

 ERROR:
  if (raw_pins) free(raw_pins);
  if (sel_pins) free(sel_pins);
  if (kmin)     free(kmin);
  if (kmax)     free(kmax);
  if (gd_ok)    pb_gap_data_free(&gd);
  return status;
}


/* cm_p7_om_holder_Init() / cm_p7_om_holder_Reset():
 * Lifecycle helpers for the reusable LOCAL p7 profile + OPROFILE used by
 * the --p7pinbridge SW scan (brief 26_0430-090). Init zeroes the holder; the
 * profile/OPROFILE are built lazily on the first pinbridge sequence inside
 * p7_Seq2BandsPinBridgeWrap(). Reset frees the built objects (no-op if never
 * built). One holder per worker thread; never share across threads.
 */
void
cm_p7_om_holder_Init(CM_P7_OM_HOLDER *h)
{
  if (h == NULL) return;
  h->gm_local = NULL;
  h->om       = NULL;
  h->om32     = NULL;
  h->M        = 0;
  h->built    = FALSE;
}

void
cm_p7_om_holder_Reset(CM_P7_OM_HOLDER *h)
{
  if (h == NULL) return;
  if (h->om32     != NULL) cm_pb_om32_Destroy(h->om32);
  if (h->om       != NULL) p7_oprofile_Destroy(h->om);
  if (h->gm_local != NULL) p7_profile_Destroy(h->gm_local);
  h->gm_local = NULL;
  h->om       = NULL;
  h->om32     = NULL;
  h->M        = 0;
  h->built    = FALSE;
}

/* p7_Seq2BandsPinBridgeWrap():
 * Mirrors p7_Seq2BandsVit() signature exactly; only the band-derivation
 * step is swapped (SW-pinbridge prefilter + banded p7 Viterbi instead of
 * full unbanded p7_GViterbi). i2k extraction + p7_pins2bands_nodepad are
 * the same as p7_Seq2BandsVit.
 *
 * The SSE scan needs a LOCAL config of cm->fp7 (rbv requires LOCAL). That
 * LOCAL profile/OPROFILE depends only on the model, not the residues, so
 * when <om_holder> is non-NULL it is built once (lazily, on the first
 * sequence) and reused across the whole block/thread, with only a per-seq
 * p7_oprofile_ReconfigLength(). When <om_holder> is NULL the profile is
 * built and freed per call (single-sequence callers). The Viterbi runs
 * against gm in the caller's configured mode (GLOCAL or truncated LOCAL).
 */
int
p7_Seq2BandsPinBridgeWrap(CM_t *cm, char *errbuf, P7_PROFILE *gm,
                          P7_BG *bg, P7_TRACE *p7_tr,
                          ESL_DSQ *dsq, int L, int pad, int *nodepad,
                          int hopback, int vitend,
                          CM_P7_OM_HOLDER *om_holder,
                          int **ret_i2k, int **ret_kmin, int **ret_kmax, int *ret_ncells)
{
  int          status;
  float        sc;
  int         *i2k    = NULL;
  int         *kmin   = NULL;
  int         *kmax   = NULL;
  int          ncells = 0;
  int          M      = gm->M;
  int          tpos;
  P7_PROFILE  *gm_local = NULL;
  P7_OPROFILE *om       = NULL;
  CM_PB_OM32  *om32     = NULL;   /* 32-bit emission table for the 32-bit SW scan (brief 26_0430-094) */
  int          own_om   = FALSE; /* TRUE if we built gm_local/om/om32 locally (no holder) */
  P7_GMXB     *gxb      = NULL;
  P7_GBANDS   *bnd      = NULL;
  int          pb_pad   = (cm->p7_pinbridge_pad > 0) ? cm->p7_pinbridge_pad : 20;
  struct timespec ta, tb;
  double       om_ms = 0;
  double       sw_ms = 0, lsis_ms = 0, band_ms = 0, bvit_ms = 0, btrace_ms = 0;
  double       pins2bands_ms = 0;

  if (cm->fp7 == NULL) ESL_FAIL(eslEINVAL, errbuf, "p7_Seq2BandsPinBridgeWrap: cm->fp7 is NULL");

  /* PB_LOAD_BAND hook (brief 26_0430-117): replace the entire SW+LSIS+pb_build_band
   * pin-extraction pipeline with an offline-generated per-row band loaded
   * from a TSV file. Format: "i kmin kmax" per line, lines starting with
   * '#' are skipped. Rows i=1..L must be present; missing rows default to
   * [1, M]. i2k is returned as all -1 (downstream consumers use kmin/kmax
   * only; p7_kbands2gbands ignores i2k). */
  {
    const char *band_path = getenv("PB_LOAD_BAND");
    if (band_path != NULL) {
      FILE *fp = fopen(band_path, "r");
      char  line[512];
      int   ii;
      if (fp == NULL) ESL_FAIL(eslFAIL, errbuf, "PB_LOAD_BAND: cannot open '%s'", band_path);
      ESL_ALLOC(i2k,  sizeof(int) * (L + 1));
      ESL_ALLOC(kmin, sizeof(int) * (L + 1));
      ESL_ALLOC(kmax, sizeof(int) * (L + 1));
      for (ii = 0; ii <= L; ii++) { i2k[ii] = -1; kmin[ii] = 1; kmax[ii] = M; }
      /* B state at i=0: cp9_FB2HMMBandsP7BF asserts kmin[0]==0. */
      kmin[0] = 0; kmax[0] = 0;
      while (fgets(line, sizeof(line), fp) != NULL) {
        int li, lkmin, lkmax;
        if (line[0] == '#') continue;
        if (sscanf(line, "%d %d %d", &li, &lkmin, &lkmax) != 3) continue;
        if (li < 0 || li > L) continue;
        if (lkmin < 1) lkmin = 1;
        if (lkmax > M) lkmax = M;
        if (lkmin > lkmax) lkmin = lkmax;
        kmin[li] = lkmin;
        kmax[li] = lkmax;
      }
      fclose(fp);
      ncells = 0;
      for (ii = 1; ii <= L; ii++) ncells += kmax[ii] - kmin[ii] + 1;
      fprintf(stderr, "#P7BAND PB_LOAD_BAND M=%d L=%d ncells=%d avg_bw=%.1f\n",
              M, L, ncells, (float)ncells / (float)L);
      *ret_i2k    = i2k;
      *ret_kmin   = kmin;
      *ret_kmax   = kmax;
      *ret_ncells = ncells;
      return eslOK;
    }
  }

  /* Acquire the LOCAL p7 profile + OPROFILE for the SSE rbv scan.
   *
   * The LOCAL config of cm->fp7 depends only on the model, not the residues,
   * so when a holder is provided we build it once (lazily) and reuse it
   * across the block/thread, paying only a per-sequence ReconfigLength().
   * Without a holder (single-sequence callers) we build-and-free per call,
   * exactly as before. Either way the per-sequence om is byte-identical:
   * p7_oprofile_Convert produces L-independent core scores and
   * p7_oprofile_ReconfigLength(om, L) recomputes all L-dependent specials. */
  clock_gettime(CLOCK_MONOTONIC, &ta);
  if (om_holder != NULL) {
    if (! om_holder->built) {
      om_holder->gm_local = p7_profile_Create(M, cm->fp7->abc);
      if (om_holder->gm_local == NULL) ESL_FAIL(eslEMEM, errbuf, "p7_profile_Create failed");
      if ((status = p7_ProfileConfig(cm->fp7, bg, om_holder->gm_local, L, p7_LOCAL)) != eslOK)
        ESL_XFAIL(status, errbuf, "p7_ProfileConfig (LOCAL) failed");
      om_holder->om = p7_oprofile_Create(M, cm->fp7->abc);
      if (om_holder->om == NULL) ESL_XFAIL(eslEMEM, errbuf, "p7_oprofile_Create failed");
      if ((status = p7_oprofile_Convert(om_holder->gm_local, om_holder->om)) != eslOK)
        ESL_XFAIL(status, errbuf, "p7_oprofile_Convert failed");
      /* 32-bit emission table (brief 26_0430-094): L-independent (emission scores
       * only), so built once and reused across the block with no per-seq
       * reconfig. Built from the same LOCAL profile as om. */
      if ((status = cm_pb_om32_Create(M, cm->fp7->abc, &om_holder->om32)) != eslOK)
        ESL_XFAIL(status, errbuf, "cm_pb_om32_Create failed");
      if ((status = cm_pb_om32_Build(om_holder->gm_local, om_holder->om32)) != eslOK)
        ESL_XFAIL(status, errbuf, "cm_pb_om32_Build failed");
      om_holder->M     = M;
      om_holder->built = TRUE;
    }
    gm_local = om_holder->gm_local;
    om       = om_holder->om;
    om32     = om_holder->om32;
    own_om   = FALSE;
    /* Per-sequence length reconfig on the reused OPROFILE (the only L-dependent
     * step). gm_local is only the Convert source and is unused downstream, so
     * it needs no per-seq reconfig. */
    if ((status = p7_oprofile_ReconfigLength(om, L)) != eslOK)
      ESL_XFAIL(status, errbuf, "p7_oprofile_ReconfigLength failed");
  } else {
    /* No holder: build LOCAL p7 profile + OPROFILE for this call only. */
    gm_local = p7_profile_Create(M, cm->fp7->abc);
    if (gm_local == NULL) ESL_FAIL(eslEMEM, errbuf, "p7_profile_Create failed");
    if ((status = p7_ProfileConfig(cm->fp7, bg, gm_local, L, p7_LOCAL)) != eslOK)
      ESL_XFAIL(status, errbuf, "p7_ProfileConfig (LOCAL) failed");
    om = p7_oprofile_Create(M, cm->fp7->abc);
    if (om == NULL) ESL_XFAIL(eslEMEM, errbuf, "p7_oprofile_Create failed");
    if ((status = p7_oprofile_Convert(gm_local, om)) != eslOK)
      ESL_XFAIL(status, errbuf, "p7_oprofile_Convert failed");
    if ((status = p7_oprofile_ReconfigLength(om, L)) != eslOK)
      ESL_XFAIL(status, errbuf, "p7_oprofile_ReconfigLength failed");
    /* 32-bit emission table (brief 26_0430-094); freed with om/gm_local below. */
    if ((status = cm_pb_om32_Create(M, cm->fp7->abc, &om32)) != eslOK)
      ESL_XFAIL(status, errbuf, "cm_pb_om32_Create failed");
    if ((status = cm_pb_om32_Build(gm_local, om32)) != eslOK)
      ESL_XFAIL(status, errbuf, "cm_pb_om32_Build failed");
    own_om = TRUE;
  }

  /* Allocate banded matrix scratch (Reinit'd inside p7_Seq2BandsPinBridge once
   * the prefilter band is known). */
  bnd = p7_gbands_Create();
  if (bnd == NULL) ESL_XFAIL(eslEMEM, errbuf, "p7_gbands_Create failed");
  /* Allocate gxb shell only (NULL dp/xmx); p7_Seq2BandsPinBridge fills it
   * once bnd has cells. p7_gmxb_Create disallows zero-cell bnd. */
  if ((gxb = malloc(sizeof(P7_GMXB))) == NULL) ESL_XFAIL(eslEMEM, errbuf, "malloc P7_GMXB shell failed");
  gxb->dp = NULL; gxb->xmx = NULL; gxb->bnd = NULL; gxb->dalloc = 0; gxb->xalloc = 0;
  clock_gettime(CLOCK_MONOTONIC, &tb);
  om_ms = (tb.tv_sec - ta.tv_sec)*1000.0 + (tb.tv_nsec - ta.tv_nsec)/1e6;

  /* Step 1+2: SW-pinbridge prefilter + banded p7 Viterbi + banded trace.
   * cm->p7_pinbridge_vit_gaps selects mini-Viterbi gap costs (Option 3). */
  status = p7_Seq2BandsPinBridge(gm, om, om32, gxb, bnd, p7_tr, dsq, L, pb_pad,
                                 cm->p7_pinbridge_vit_gaps, &sc,
                                 &sw_ms, &lsis_ms, &band_ms, &bvit_ms, &btrace_ms);
  if (status != eslOK) {
    /* signal caller to fall back */
    *ret_i2k    = NULL;
    *ret_kmin   = NULL;
    *ret_kmax   = NULL;
    *ret_ncells = 0;
    p7_gmxb_Destroy(gxb);
    p7_gbands_Destroy(bnd);
    if (own_om) { cm_pb_om32_Destroy(om32); p7_oprofile_Destroy(om); p7_profile_Destroy(gm_local); }
    return eslOK;
  }

  /* PB_DUMP_P7VIT_TRACE diagnostic (brief 26_0430-103): dump the banded p7 Viterbi
   * trace (from p7_Seq2BandsPinBridge's banded traceback) to a TSV file.
   * Columns: tpos  state  i  k  residue_i (y/-) */
  { const char *_vit_path = getenv("PB_DUMP_P7VIT_TRACE");
    if (_vit_path != NULL) {
      FILE *_vfp = fopen(_vit_path, "w");
      if (_vfp != NULL) {
        int _t;
        fprintf(_vfp, "# p7 Viterbi trace dump (brief 103 / port of brief 100)\n");
        fprintf(_vfp, "# M=%d L=%d N=%d\n", M, L, p7_tr->N);
        fprintf(_vfp, "#tpos\tstate\ti\tk\tresidue_i\n");
        for (_t = 0; _t < p7_tr->N; _t++) {
          char _stc = '?';
          switch (p7_tr->st[_t]) {
            case p7T_M: _stc = 'M'; break; case p7T_D: _stc = 'D'; break;
            case p7T_I: _stc = 'I'; break; case p7T_S: _stc = 'S'; break;
            case p7T_N: _stc = 'N'; break; case p7T_B: _stc = 'B'; break;
            case p7T_E: _stc = 'E'; break; case p7T_C: _stc = 'C'; break;
            case p7T_T: _stc = 'T'; break; case p7T_J: _stc = 'J'; break;
            default: _stc = '?'; break;
          }
          fprintf(_vfp, "%d\t%c\t%d\t%d\t%s\n",
                  _t, _stc, p7_tr->i[_t], p7_tr->k[_t],
                  (p7_tr->i[_t] > 0) ? "y" : "-");
        }
        fclose(_vfp);
      }
    }
  }

  /* Step 3: i2k from M-state trace pins (identical to p7_Seq2BandsVit). */
  ESL_ALLOC(i2k, sizeof(int) * (L + 1));
  esl_vec_ISet(i2k, (L + 1), -1);

  for (tpos = 0; tpos < p7_tr->N; tpos++) {
    if (p7_tr->st[tpos] == p7T_M) {
      int it = p7_tr->i[tpos];
      int kt = p7_tr->k[tpos];
      if (it >= 1 && it <= L && kt >= 1 && kt <= M)
        i2k[it] = kt;
    }
  }

  /* Step 3b: vitend pruning (mirror p7_Seq2BandsVit) */
  if (vitend > 0) {
    int dropped, ii;
    dropped = 0;
    for (ii = 1; ii <= L && dropped < vitend; ii++) {
      if (i2k[ii] != -1) { i2k[ii] = -1; dropped++; }
    }
    dropped = 0;
    for (ii = L; ii >= 1 && dropped < vitend; ii--) {
      if (i2k[ii] != -1) { i2k[ii] = -1; dropped++; }
    }
  }

  /* Step 4: pins -> bands (identical to p7_Seq2BandsVit) */
  clock_gettime(CLOCK_MONOTONIC, &ta);
  if (nodepad != NULL) {
    /* brief 26_0628-043: pinbridge path doesn't expose a tunable ramp alpha -- pass
     * brief 26_0628-042's validated default unconditionally. */
    if ((status = p7_pins2bands_nodepad(i2k, errbuf, L, M, nodepad, hopback, 0.75, &kmin, &kmax, &ncells)) != eslOK)
      goto ERROR;
  } else {
    if ((status = p7_pins2bands(i2k, errbuf, L, M, pad, &kmin, &kmax, &ncells)) != eslOK)
      goto ERROR;
  }
  clock_gettime(CLOCK_MONOTONIC, &tb);
  pins2bands_ms = (tb.tv_sec - ta.tv_sec)*1000.0 + (tb.tv_nsec - ta.tv_nsec)/1e6;

  /* PB_DUMP_BAND diagnostic (brief 26_0430-103): dump final kmin/kmax band to file.
   * Env var PB_DUMP_BAND = path to output TSV. Columns: i, kmin[i], kmax[i].
   * This is the band returned to the caller (used for CP9 F/B posterior).
   * Appends to the file (so multiple sequences accumulate; use PB_DUMP_TAG
   * in the caller to separate runs, or point to per-seq files). */
  { const char *_pb_band_path = getenv("PB_DUMP_BAND");
    if (_pb_band_path != NULL) {
      FILE *_pb_band_fp = fopen(_pb_band_path, "w");
      if (_pb_band_fp != NULL) {
        int _pb_bi;
        fprintf(_pb_band_fp, "# M=%d L=%d\n", M, L);
        for (_pb_bi = 1; _pb_bi <= L; _pb_bi++)
          fprintf(_pb_band_fp, "%d\t%d\t%d\n", _pb_bi, kmin[_pb_bi], kmax[_pb_bi]);
        fclose(_pb_band_fp);
      }
    }
  }

  /* Emit per-stage timing for offline aggregation by CLEN bucket. */
  fprintf(stderr, "#P7PB_STAGE M=%d L=%d om_build=%.4f sw_scan=%.4f lsis=%.4f band_build=%.4f banded_vit=%.4f banded_trace=%.4f pins2bands=%.4f\n",
          M, L, om_ms/1000.0, sw_ms/1000.0, lsis_ms/1000.0, band_ms/1000.0, bvit_ms/1000.0, btrace_ms/1000.0, pins2bands_ms/1000.0);
  fflush(stderr);

  *ret_i2k    = i2k;
  *ret_kmin   = kmin;
  *ret_kmax   = kmax;
  *ret_ncells = ncells;

  p7_gmxb_Destroy(gxb);
  p7_gbands_Destroy(bnd);
  if (own_om) { cm_pb_om32_Destroy(om32); p7_oprofile_Destroy(om); p7_profile_Destroy(gm_local); }
  return eslOK;

 ERROR:
  if (i2k)  free(i2k);
  if (kmin) free(kmin);
  if (kmax) free(kmax);
  if (gxb)      p7_gmxb_Destroy(gxb);
  if (bnd)      p7_gbands_Destroy(bnd);
  /* Only free the LOCAL profile/OPROFILE/om32 if we own them; a holder's
   * objects are owned (and freed) by the caller via cm_p7_om_holder_Reset(),
   * which also cleans up a holder left partially built by a failure above. */
  if (own_om) {
    if (om32)     cm_pb_om32_Destroy(om32);
    if (om)       p7_oprofile_Destroy(om);
    if (gm_local) p7_profile_Destroy(gm_local);
  }
  return status;
}

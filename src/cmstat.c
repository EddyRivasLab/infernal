/************************************************************
 * @LICENSE@
 ************************************************************/

/* cmstat.c
 * EPN, Tue Aug 21 12:50:34 2007
 *
 * Display summary statistics for an CM or CM database 
 * (such as Rfam). 
 *
 * Based on SRE's hmmstat.c from HMMER3.
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "easel.h"
#include <esl_getopts.h>
#include <esl_histogram.h>
#include <esl_random.h>
#include <esl_stats.h>
#include <esl_stopwatch.h>
#include <esl_vectorops.h>
#include <esl_wuss.h>

#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* function declarations                */
#include "stats.h"
#include "cm_dispatch.h"	


static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs       incomp  help  docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",   1 },
  { "-1",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "use tabular output summary format, 1 line per CM", 1 },
  { "-t",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "perform and report on timing experiments", 1 },
  { "-L",        eslARG_INT,    "1000",NULL, "n>0",     NULL,      "-t",        NULL, "with -t, length of sequences for search    stats", 1 },
  { "-N",        eslARG_INT,    "25",  NULL, "n>0",     NULL,      "-t",        NULL, "with -t, number of sequences for alignment stats", 1 },
  { "--beta",    eslARG_REAL,   "1E-7",NULL, "0<x<1",   NULL,      NULL,        NULL, "set tail loss prob for QDB stats to <x>", 1 },
  { "--tau",     eslARG_REAL,   "1E-7",NULL, "0<x<1",   NULL,      NULL,       NULL,  "set tail loss prob for HMM banded stats to <x>", 1 },
  { "--exp",     eslARG_REAL,   NULL,  NULL, "0<x<=1.0",NULL,      NULL,        NULL, "exponentiate CM probabilities by <x> before calc'ing stats",  1 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <cmfile>";
static char banner[] = "display summary statistics for CMs";

static int    summarize_search(ESL_GETOPTS *go, CM_t *cm, ESL_RANDOMNESS *r, ESL_STOPWATCH *w); 
static int    summarize_alignment(ESL_GETOPTS *go, CM_t *cm);
static float  count_scan_dp_calcs(CM_t *cm, int L, int use_qdb);
static double cm_MeanMatchRelativeEntropy(const CM_t *cm);
static double cm_MeanMatchEntropy(const CM_t *cm);
static double cm_MeanMatchInfo(const CM_t *cm);

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go = NULL;   /* command line processing   */
  ESL_ALPHABET    *abc = NULL;  /* alphabet                  */
  ESL_RANDOMNESS  *r   = NULL;  /* source of randomness      */
  ESL_STOPWATCH   *w   = NULL;  /* for timings               */
  char            *cmfile;	/* name of input CM file     */ 
  CMFILE          *cmfp;	/* open input CM file stream */
  CM_t            *cm;          /* CM most recently read     */
  int              ncm;         /* CM index                  */

  /* Process command line options.
   */
  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK || 
      esl_opt_VerifyConfig(go)               != eslOK)
    {
      printf("Failed to parse command line: %s\n", go->errbuf);
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      cm_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\nwhere options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1=docgroup, 2 = indentation; 80=textwidth*/
      exit(0);
    }
  if (esl_opt_ArgNumber(go) != 1) 
    {
      puts("Incorrect number of command line arguments.");
      esl_usage(stdout, argv[0], usage);
      puts("\n  where options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      printf("\nTo see more help on other available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  cmfile     = esl_opt_GetArg(go, 1); 
  if(esl_opt_GetBoolean(go, "-t")) { 
    r = esl_randomness_CreateTimeseeded();
    w = esl_stopwatch_Create();
  }
  /* Initializations: open the CM file
   */
  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL)
    esl_fatal("Failed to open covariance model save file %s\n", cmfile);

  /* Main body: read CMs one at a time, print one line of stats.
   */
  printf("#\n");
  printf("# %-4s %-20s %8s %8s %6s %5s %5s %3s %6s\n", "idx",  "name",                 "nseq",     "eff_nseq", "clen",   "W",     "M",     "bif", "relent");
  printf("# %-4s %-20s %8s %8s %6s %5s %5s %3s %6s\n", "----", "--------------------", "--------", "--------", "------", "-----", "-----", "---", "------");

  ncm = 0;
  while (CMFileRead(cmfp, &abc, &cm))
  {
    ncm++;
    
    cm->beta = esl_opt_GetReal(go, "--beta");
    ConfigQDB(cm);
    
    printf("%6d %-20s %8d %8.2f %6d %5d %5d %3d %6.2f\n",
	   ncm,
	   cm->name,
	   cm->nseq,
	   cm->eff_nseq,
	   cm->clen,
	   cm->W,
	   cm->M,
	   CMCountStatetype(cm, B_st),
	   cm_MeanMatchRelativeEntropy(cm));
	   /*cm_MeanMatchInfo(cm));*/

    if(! esl_opt_GetBoolean(go, "-1")) {
      summarize_search(go, cm, r, w);
      summarize_alignment(go, cm);
    }
    FreeCM(cm);
  }
    
  esl_alphabet_Destroy(abc);
  CMFileClose(cmfp);
  esl_getopts_Destroy(go);
  exit(0);
}

/* Function:  summarize_search()
 * Incept:    EPN, Tue Aug 21 20:00:28 2007
 *
 * Purpose:   Summarize search statistics to varying extents
 *            based on command-line options.
 */
int
summarize_search(ESL_GETOPTS *go, CM_t *cm, ESL_RANDOMNESS *r, ESL_STOPWATCH *w) 
{
  int status;     
  int L = esl_opt_GetInteger(go, "-L"); /* length sequence to search */
  float dpc;       /* number of     mega-DP calcs for search of length L */
  float dpc_q;     /* number of QDB mega-DP calcs for search of length L */
  float th_acc;    /* theoretical QDB acceleration */

  /* optional, -t related variables */
  ESL_DSQ *dsq;    /* digitized sequence of length L for timings  */
  float t_c;       /* number of seconds (w->user) for        CYK search */
  float t_i;       /* number of seconds (w->user) for     Inside search */
  float t_cq;      /* number of seconds (w->user) for QDB    CYK search */
  float t_iq;      /* number of seconds (w->user) for QDB Inside search */

  if(L < cm->W) { L = cm->W; printf("\tL increased to minimum size of cm->W (%d)\n", L); }
  if(esl_opt_GetBoolean(go, "-t")) { 
    ESL_ALLOC(dsq, sizeof(ESL_DSQ) * L+2);
    esl_rnd_xfIID(r, cm->null, cm->abc->K, L, dsq);
  }

  /* estimate speedup due to QDB */
  dpc    = count_scan_dp_calcs(cm, L, FALSE) / 1000000.;
  dpc_q  = count_scan_dp_calcs(cm, L, TRUE)  / 1000000.;
  th_acc = dpc / dpc_q;

  /* print simple stats if -t not enabled */
  printf("\tL: %d DP megacalcs: non-banded: %6.2f QDB: %6.2f th. speedup: %6.2fX\n", L, dpc, dpc_q, th_acc);
 else /* -t enabled, do timings  */
    {
      /* cyk */
      esl_stopwatch_Start(w);
      CYKScan (cm, dsq, 1, L, cm->W, 0., NULL);
      esl_stopwatch_Stop(w);
      t_c = w->user;
      /* qdb cyk */
      esl_stopwatch_Start(w);
      CYKBandedScan (cm, dsq, cm->dmin, cm->dmax, 1, L, cm->W, 0., NULL);
      esl_stopwatch_Stop(w);
      t_cq = w->user;
      /* inside */
      esl_stopwatch_Start(w);
      iInsideScan (cm, dsq, 1, L, cm->W, 0., NULL);
      esl_stopwatch_Stop(w);
      t_i = w->user;
      /* inside */
      esl_stopwatch_Start(w);
      iInsideBandedScan (cm, dsq, cm->dmin, cm->dmax, 1, L, cm->W, 0., NULL);
      esl_stopwatch_Stop(w);
      t_iq = w->user;
      printf("\t   CYK:\t%6.2fX QDB spdup\t%6.2f megacalcs/s\t %6.0f res/s\n", (t_c/t_cq), (dpc_q/t_cq), (L/t_cq));
      printf("\tInside:\t%6.2fX QDB spdup\t%6.2f megacalcs/s\t %6.0f res/s\n", (t_i/t_iq), (dpc_q/t_iq), (L/t_iq));
  }
  return eslOK;

 ERROR:
  return status;
}

/* Function:  summarize_alignment()
 * Incept:    
 *
 * Purpose:   Summarize alignment statistics to varying extents
 *            based on command-line options.
 */
int
summarize_alignment(ESL_GETOPTS *go, CM_t *cm)
{
  /* HERE: do HMM banded alignment stats
   * sample N=100 seqs, and calculate posteriors, determine new
   * number of CYK DP calcs AND CP9 F/B calcs to get bands. */
  return eslOK;
}

/* Function: count_scan_dp_calcs()
 * Date:     EPN, Wed Aug 22 09:08:03 2007
 *
 * Purpose:  Count all DP calcs for a CM scan against a 
 *           sequence of length L. Similar to smallcyk.c's
 *           CYKDemands() but takes into account number of
 *           transitions from each state, and is concerned
 *           with a scanning dp matrix, not an alignment matrix.
 *
 * Args:     cm     - the model
 *           L      - length of sequence
 *           use_qdb- TRUE to enforce cm->dmin and cm->dmax for calculation
 *
 * Returns: (float) the total number of DP calculations, either using QDB or not.
 */
float
count_scan_dp_calcs(CM_t *cm, int L, int use_qdb)
{
  int v, j;
  float dpcalcs = 0.;
  float dpcalcs_bif = 0.;

  /* Contract check */
  if(cm->W > L) esl_fatal("ERROR in count_scan_dp_calcs(), cm->W: %d exceeds L: %d\n", cm->W, L);

  float  dpcells     = 0.;
  float  bif_dpcells = 0.;
  int d,w,y,kmin,kmax, bw;

  if(! use_qdb) 
    {
      dpcells = (cm->W * L) - (cm->W * (cm->W-1) * .5); /* fillable dp cells per state (deck) */
      for (j = 1; j < cm->W; j++)
	bif_dpcells += ((j    +2) * (j    +1) * .5) - 1;
      for (j = cm->W; j <= L; j++)
	bif_dpcells += ((cm->W+2) * (cm->W+1) * .5) - 1;
      dpcalcs_bif = CMCountStatetype(cm, B_st) * bif_dpcells;
      for(v = 0; v < cm->M; v++)
	{
	  if(cm->sttype[v] != B_st && cm->sttype[v] != E_st)
	    {
	      dpcalcs += dpcells * cm->cnum[v];

	      /* non-obvious subtractions that match implementation in scancyk.c::CYKScan() */
	      if(v == 0) dpcalcs  -= dpcells; 
	      if(cm->sttype[v] == MP_st) dpcalcs  -= L * cm->cnum[v];
	    }
	}
    }
  else /* use_qdb */
    {
      for(v = 0; v < cm->M; v++)
	{
	  if(cm->sttype[v] != B_st && cm->sttype[v] != E_st)
	    {
	      bw = cm->dmax[v] - cm->dmin[v] + 1; /* band width */
	      if(cm->dmin[v] == 0) bw--;
	      dpcalcs += ((bw * L) - (bw * (bw-1) * 0.5)) * cm->cnum[v];

	      /* non-obvious subtractions that match implementation in bandcyk.c::CYKBandedScan() */
	      if(v == 0) dpcalcs  -= ((bw * L) - (bw * (bw-1) * 0.5)); 
	      if(cm->sttype[v] == MP_st) dpcalcs  -= bw * cm->cnum[v];
	    }

	  else if(cm->sttype[v] == B_st)
	    {
	      w = cm->cfirst[v];
	      y = cm->cnum[v];
	      for (j = 1; j <= L; j++)
		{
		  d = (cm->dmin[v] > 0) ? cm->dmin[v] : 1;
		  for (; d <= cm->dmax[v] && d <= j; d++)
		    {
		      if(cm->dmin[y] > (d-cm->dmax[w])) kmin = cm->dmin[y];
		      else kmin = d-cm->dmax[w];
		      if(kmin < 0) kmin = 0;
		      if(cm->dmax[y] < (d-cm->dmin[w])) kmax = cm->dmax[y];
		      else kmax = d-cm->dmin[w];
		      if(kmin <= kmax)
			{
			  bw = (kmax - kmin + 1);
			  dpcalcs_bif += bw;
			}
		    }
		}
	    }
	}
    }
  /*printf("%d count_scan_dp_calcs dpc     %.0f\n", use_qdb, dpcalcs);
    printf("%d count_scan_dp_calcs dpc_bif %.0f\n", use_qdb, dpcalcs_bif);
    printf("%d count_scan_dp_calcs total   %.0f\n", use_qdb, dpcalcs + dpcalcs_bif);*/
  return dpcalcs + dpcalcs_bif;
}

/* Function:  cm_MeanMatchInfo()
 * Incept:    SRE, Fri May  4 11:43:56 2007 [Janelia]
 *
 * Purpose:   Calculate the mean information content per match state
 *            emission distribution, in bits:
 *            
 *            \[
 *              \frac{1}{M} \sum_{k=1}^{M}
 *                \left[ 
 *                   - \sum_x p_k(x) \log_2 p_k(x) 
 *                   + \sum_x f(x) \log_2 f(x)
 *                \right] 
 *            \]
 *            
 *            where $p_k(x)$ is emission probability for symbol $x$
 *            from match state $k$, and $f(x)$ is the null model's
 *            background emission probability for $x$.
 *            
 *            This statistic is used in "entropy weighting" to set the
 *            total sequence weight when model building.
 */
double
cm_MeanMatchInfo(const CM_t *cm)
{
  return esl_vec_FEntropy(cm->null, cm->abc->K) - cm_MeanMatchEntropy(cm);
}

/* Function:  cm_MeanMatchEntropy()
 * Incept:    SRE, Fri May  4 13:37:15 2007 [Janelia]
 *
 * Purpose:   Calculate the mean entropy per match state emission
 *            distribution, in bits:
 *            
 *            \[
 *              \frac{1}{clen} \sum_{v=0}^{M-1} -\sum_x p_v(x) \log_2 p_v(x)
 *            \]
 *       
 *            where $p_v(x)$ is emission probability for symbol $x$
 *            from MATL\_ML, MATR\_MR or MATP\_MP state $v$. For MATP\_MP
 *            states symbols $x$ are base pairs.
 */
double
cm_MeanMatchEntropy(const CM_t *cm)
{
  int    v;
  double H = 0.;

  for (v = 0; v < cm->M; v++)
    {
      if(cm->stid[v] == MATP_MP)
       H += esl_vec_FEntropy(cm->e[v], (cm->abc->K * cm->abc->K));
      else if(cm->stid[v] == MATL_ML || 
	      cm->stid[v] == MATR_MR)
	H += esl_vec_FEntropy(cm->e[v], cm->abc->K);
    }
  H /= (double) cm->clen;
  return H;
}


/* Function:  cm_MeanMatchRelativeEntropy()
 * Incept:    SRE, Fri May 11 09:25:01 2007 [Janelia]
 *
 * Purpose:   Calculate the mean relative entropy per match state emission
 *            distribution, in bits:
 *            
 *            \[
 *              \frac{1}{M} \sum_{v=0}^{M-1} \sum_x p_v(x) \log_2 \frac{p_v(x)}{f(x)}
 *            \]
 *       
 *            where $p_v(x)$ is emission probability for symbol $x$
 *            from MATL\_ML, MATR\_MR, or MATP\_MP state state $v$, 
 *            and $f(x)$ is the null model's background emission 
 *            probability for $x$. For MATP\_MP states, $x$ is a 
 *            base pair.
 */
double
cm_MeanMatchRelativeEntropy(const CM_t *cm)
{
  int    status;
  int    v;
  double KL = 0.;
  float *pair_null;
  int i,j;
  
  ESL_ALLOC(pair_null, (sizeof(float) * cm->abc->K * cm->abc->K));
  for(i = 0; i < cm->abc->K; i++)
    for(j = 0; j < cm->abc->K; j++)
      pair_null[(i * cm->abc->K) + j] = cm->null[i] * cm->null[j]; 
  
  for (v = 0; v < cm->M; v++)
    if(cm->stid[v] == MATP_MP)
      KL += esl_vec_FRelEntropy(cm->e[v], pair_null, (cm->abc->K * cm->abc->K));
    else if(cm->stid[v] == MATL_ML || 
	    cm->stid[v] == MATR_MR)
      KL += esl_vec_FRelEntropy(cm->e[v], cm->null, cm->abc->K);
  
  free(pair_null);

  KL /= (double) cm->clen;
  return KL;
  
 ERROR:
  esl_fatal("Memory allocation error.");
  return 0.; /* NOTREACHED */
}



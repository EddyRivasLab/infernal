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
#include <esl_sqio.h>
#include <esl_stats.h>
#include <esl_stopwatch.h>
#include <esl_vectorops.h>
#include <esl_wuss.h>

#include "funcs.h"		/* function declarations                */
#include "structs.h"		/* data structures, macros, #define's   */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs       incomp  help  docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",   1 },
  { "-1",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "use tabular output summary format, 1 line per CM", 1 },
  { "--beta",    eslARG_REAL,   "1E-7",NULL, "0<x<1",   NULL,      NULL,        NULL, "set tail loss prob for QDB stats to <x>", 2 },
  /* search stats options */
  { "-s",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "do search timing experiments", 2 },
  { "-L",        eslARG_INT,    "1000",NULL, "n>0",     NULL,      "-s",        NULL, "length of sequences for CM search stats", 2 },
  { "--cp9L",    eslARG_INT,    "100000",NULL,"n>0",    NULL,      "-s",        NULL, "length of sequences for CP9 HMM search stats", 2 },
  /* alignment stats options */
  { "-a",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "do alignment timing experiments", 2 },
  { "-N",        eslARG_INT,    "25",  NULL, "n>0",     NULL,      "-a",        NULL, "number of sequences for alignment stats", 1 },
  { "--tau",     eslARG_REAL,   "1E-7",NULL, "0<x<1",   NULL,      "-a",        NULL,  "set tail loss prob for HMM banded stats to <x>", 1 },
  { "--exp",     eslARG_REAL,   NULL,  NULL, "0<x",     NULL,      "-a",        NULL, "exponentiate CM probabilities by <x> before calc'ing stats",  1 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <cmfile>";
static char banner[] = "display summary statistics for CMs";

static int    summarize_search(ESL_GETOPTS *go, char *errbuf, CM_t *cm, ESL_RANDOMNESS *r, ESL_STOPWATCH *w); 
static int    summarize_alignment(ESL_GETOPTS *go, char *errbuf, CM_t *cm, ESL_RANDOMNESS *r, ESL_STOPWATCH *w); 
static float  count_align_dp_calcs(CM_t *cm, int L);
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
  char             errbuf[cmERRBUFSIZE];
  int              status;

  /* setup logsum lookups (could do this only if nec based on options, but this is safer) */
  init_ilogsum();
  FLogsumInit();

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
      puts("\nwhere basic options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1=docgroup, 2 = indentation; 80=textwidth*/
      puts("\noptions for statistics on search:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); /* 1=docgroup, 2 = indentation; 80=textwidth*/
      puts("\noptions for statistics on alignment:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); /* 1=docgroup, 2 = indentation; 80=textwidth*/
      exit(0);
    }
  if (esl_opt_ArgNumber(go) != 1) 
    {
      puts("Incorrect number of command line arguments.");
      esl_usage(stdout, argv[0], usage);
      puts("\n  where basic options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      printf("\nTo see more help on other available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  cmfile     = esl_opt_GetArg(go, 1); 
  r = esl_randomness_CreateTimeseeded();
  w = esl_stopwatch_Create();

  /* Initializations: open the CM file
   */
  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL)
    cm_Fail("Failed to open covariance model save file %s\n", cmfile);

  /* Main body: read CMs one at a time, print one line of stats.
   */
  ncm = 0;
  while (CMFileRead(cmfp, &abc, &cm))
  {
    if(ncm == 0 || (esl_opt_GetBoolean(go, "-s") || esl_opt_GetBoolean(go, "-a"))) { 
      printf("#\n");
      printf("# %-4s %-20s %8s %8s %6s %5s %5s %3s %6s\n", "idx",  "name",                 "nseq",     "eff_nseq", "clen",   "W",     "M",     "bif", "relent");
      printf("# %-4s %-20s %8s %8s %6s %5s %5s %3s %6s\n", "----", "--------------------", "--------", "--------", "------", "-----", "-----", "---", "------");
    }

    ncm++;
    
    cm->beta = esl_opt_GetReal(go, "--beta");
    cm->tau  = esl_opt_GetReal(go, "--tau");
    cm->config_opts |= CM_CONFIG_QDB;
    /*ConfigQDB(cm);*/
    ConfigCM(cm, NULL, NULL);
    
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

    if(esl_opt_GetBoolean(go, "-s")) { if((status = summarize_search(go, errbuf, cm, r, w))    != eslOK) cm_Fail(errbuf); }
    if(esl_opt_GetBoolean(go, "-a")) { if((status = summarize_alignment(go, errbuf, cm, r, w)) != eslOK) cm_Fail(errbuf); }
    FreeCM(cm);
  }    

  esl_alphabet_Destroy(abc);
  CMFileClose(cmfp);
  esl_getopts_Destroy(go);
  return 0;
}

/* Function:  summarize_search()
 * Incept:    EPN, Tue Aug 21 20:00:28 2007
 *
 * Purpose:   Summarize search statistics to varying extents
 *            based on command-line options.
 */
int
summarize_search(ESL_GETOPTS *go, char *errbuf, CM_t *cm, ESL_RANDOMNESS *r, ESL_STOPWATCH *w) 
{
  int status;     
  int L     = esl_opt_GetInteger(go, "-L");     /* length sequence to search */
  int L_cp9 = esl_opt_GetInteger(go, "--cp9L"); /* length sequence to search with CP9 */
  float dpc;       /* number of     mega-DP calcs for search of length L */
  float dpc_q;     /* number of QDB mega-DP calcs for search of length L */
  float th_acc;    /* theoretical QDB acceleration */
  float dpc_v;     /* number of CP9 mega-DP calcs for search of length L */
  int   minL = 0;  /* minimum length can safely scan with optimized Forward(), -1 ==> any length */
  int   be_safe;   /* should we be safe, and not use optimized Forward()? */

  /* optional, -t related variables */
  ESL_DSQ *dsq;    /* digitized sequence of length L for CM  timings  */
  ESL_DSQ *dsq_cp9;/* digitized sequence of length L for CP9 timings  */
  float t_c;       /* number of seconds (w->user) for        CYK search */
  float t_i;       /* number of seconds (w->user) for     Inside search */
  float t_cq;      /* number of seconds (w->user) for QDB    CYK search */
  float t_iq;      /* number of seconds (w->user) for QDB Inside search */
  float t_v;       /* number of seconds (w->user) for CP9 Viterbi search */
  float t_f;       /* number of seconds (w->user) for CP9 Forward search */

  if(L < cm->W) { L = cm->W; printf("\tL increased to minimum size of cm->W (%d)\n", L); }
  ESL_ALLOC(dsq,     sizeof(ESL_DSQ) * L    +2);
  ESL_ALLOC(dsq_cp9, sizeof(ESL_DSQ) * L_cp9+2);
  esl_rnd_xfIID(r, cm->null, cm->abc->K, L, dsq);
  esl_rnd_xfIID(r, cm->null, cm->abc->K, L_cp9, dsq_cp9);

  /* estimate speedup due to QDB */
  dpc    = CountScanDPCalcs(cm, L, FALSE) / 1000000.;
  dpc_q  = CountScanDPCalcs(cm, L, TRUE)  / 1000000.;
  th_acc = dpc / dpc_q;

  dpc_v  = (float) (cm->clen+1) * L_cp9 * 11; /* 11 transition's queried per HMM node: 9 main model, begin,end */ 
  dpc_v /= 1000000;

  /* First create scan info for non-QDB runs */
  int *tmp_dmin = cm->dmin;
  int *tmp_dmax = cm->dmax;
  cm->dmin = NULL;
  cm->dmax = NULL;
  cm->search_opts |= CM_SEARCH_NOQDB;
  cm_CreateScanMatrixForCM(cm, TRUE, TRUE);
  if(cm->smx == NULL) cm_Fail("summarize_search(), CreateScanMatrixForCM() call failed.");
  
  /* cyk */
  /*OLDFastCYKScan(cm, dsq, NULL, NULL, 1, L, cm->W, 0., NULL, NULL, NULL);*/
  esl_stopwatch_Start(w);
  if((status = FastCYKScan(cm, errbuf, dsq, 1, L, cm->W, 0., NULL, NULL, NULL)) != eslOK) goto ERROR;
  //CYKScan (cm, dsq, 1, L, cm->W, 0., NULL);
  esl_stopwatch_Stop(w);
  t_c = w->user;

  /* inside */
  cm->search_opts |= CM_SEARCH_INSIDE;
  esl_stopwatch_Start(w);
  if((status = FastIInsideScan(cm, errbuf, dsq, 1, L, cm->W, 0., NULL, NULL, NULL)) != eslOK) goto ERROR;
  /* iInsideScan (cm, dsq, 1, L, cm->W, 0., NULL); */
  esl_stopwatch_Stop(w);
  t_i = w->user;

  /* reset cm->dmin, cm->dmax, recalc scanmatrix */
  cm->dmin = tmp_dmin;
  cm->dmax = tmp_dmax;
  cm->search_opts &= ~CM_SEARCH_NOQDB;
  cm->search_opts &= ~CM_SEARCH_INSIDE;
  cm_FreeScanMatrixForCM(cm);
  cm_CreateScanMatrixForCM(cm, TRUE, TRUE);
  if(cm->smx == NULL) cm_Fail("summarize_search(), CreateScanMatrix() call failed.");

  /* qdb cyk */
  /*XFastCYKScan(cm, dsq, cm->dmin, cm->dmax, 1, L, cm->W, 0., NULL, NULL, NULL);*/
  esl_stopwatch_Start(w);
  /*CYKBandedScan (cm, dsq, cm->dmin, cm->dmax, 1, L, cm->W, 0., NULL); */
  if((status = FastCYKScan(cm, errbuf, dsq, 1, L, cm->W, 0., NULL, NULL, NULL)) != eslOK) goto ERROR;
  esl_stopwatch_Stop(w);
  t_cq = w->user;

  /* qdb inside */
  cm->search_opts |= CM_SEARCH_INSIDE;
  esl_stopwatch_Start(w);
  if((status = FastIInsideScan(cm, errbuf, dsq, 1, L, cm->W, 0., NULL, NULL, NULL)) != eslOK) goto ERROR;
  /*iInsideBandedScan (cm, dsq, cm->dmin, cm->dmax, 1, L, cm->W, 0., NULL);*/
  esl_stopwatch_Stop(w);
  t_iq = w->user;
  
  /* CP9 viterbi */
  esl_stopwatch_Start(w);
  if((status = cp9_FastViterbi(cm, errbuf, cm->cp9_mx, dsq_cp9, 1, L_cp9, cm->W, 0., NULL,
			       TRUE,   /* we're scanning */
			       FALSE,  /* we're not ultimately aligning */
			       TRUE,   /* be memory efficient */
			       NULL, NULL,
			       NULL,   /* don't want traces back */
			       NULL)) != eslOK) goto ERROR;

  esl_stopwatch_Stop(w);
  t_v = w->user;

  /* CP9 forward */
  esl_stopwatch_Start(w);
  /* determine the minimum length we can search safely with the optimized forward implementation. */
  if((status = cp9_WorstForward(cm, errbuf, cm->cp9_mx, -INFTY, TRUE, FALSE, &minL)) != eslOK) goto ERROR;
  /*CP9Forward(cm, dsq_cp9, 1, L_cp9, cm->W, 0., NULL, NULL, NULL,*/
  be_safe = FALSE;
  ESL_DPRINTF1(("minL: %d L: %d\n", minL, L));
  if(minL != -1 && minL <= L) be_safe = TRUE;
  esl_stopwatch_Start(w);
  if((status = Xcp9_FastForward(cm, errbuf, cm->cp9_mx, dsq_cp9, 1, L_cp9, cm->W, 0., NULL, 
				TRUE,   /* we are scanning */
				FALSE,  /* we are not ultimately aligning */
				TRUE,   /* be memory efficient */
				be_safe,
				NULL, NULL, NULL)) != eslOK) goto ERROR;
  esl_stopwatch_Stop(w);
  t_f = w->user;

  /* print results */
  float mc_s; /* million calcs / second */
  float kb_s; /* kilobases / second */
  float emp_acc; /* empirical acceleration from QDB */
  float L_kb     = (float) L / 1000.;
  float L_cp9_kb = (float) L_cp9 / 1000.;
  printf("#\n");
  printf("#\t\t\t Search statistics:\n");
  printf("#\t\t\t %7s %6s %6s %6s %8s %5s %5s\n",             "alg",     "qdb Mc", "L (kb)",  "Mc/s",     "kb/s",   "qdbXt", "qdbXa");
  printf("#\t\t\t %7s %6s %6s %6s %8s %5s %5s\n",             "-------", "------", "------","------", "--------", "------", "-----");
  mc_s = dpc_q / t_cq; 
  kb_s = ((float) L_kb) / t_cq; 
  emp_acc = t_c / t_cq; 
  printf(" \t\t\t %7s %6.1f %6.2f %6.1f %8.2f %5.1f %5.1f\n", "cyk",      dpc_q,   L_kb,    mc_s,     kb_s,       th_acc,   emp_acc);
  mc_s = dpc_q / t_iq; 
  kb_s = ((float) L_kb) / t_iq; 
  emp_acc = t_i / t_iq; 
  printf(" \t\t\t %7s %6s %6s %6.1f %8.2f %5s %5.1f\n",       "inside",   "\"",    "\"",    mc_s,     kb_s,        "\"",    emp_acc);
  mc_s = dpc_v / t_v; 
  kb_s = ((float) L_cp9_kb) / t_v; 
  printf(" \t\t\t %7s %6.1f %6.2f %6.1f %8.2f %5s %5s\n",     "viterbi",  dpc_v,   L_cp9_kb,mc_s,     kb_s,        "",      "");
  mc_s = dpc_v / t_f; 
  kb_s = ((float) L_cp9_kb) / t_f; 
  printf(" \t\t\t %7s %6s %6s %6.1f %8.2f %5s %5s\n",         "forward",  "\"",    "\"",    mc_s,     kb_s,        "",      "");
  
  free(dsq);
  free(dsq_cp9);
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
summarize_alignment(ESL_GETOPTS *go, char *errbuf, CM_t *cm, ESL_RANDOMNESS *r, ESL_STOPWATCH *w) 
{
  /* HERE: do HMM banded alignment stats
   * sample N=100 seqs, and calculate posteriors, determine new
   * number of CYK DP calcs AND CP9 F/B calcs to get bands. */
  int status;
  float dpc;  /* # DP calcs for non-banded alignment of consensus */
  CMConsensus_t *con = NULL;            /* consensus info for the CM */
  ESL_SQ *csq = NULL;
  float t_dc; /* user seconds time for D&C alignment */
  float t_hb; /* user seconds time for HMM banded alignment */
  float mc_s; /* million calcs/second */

  /* Create and align consensus sequence for D&C stats */
  CreateCMConsensus(cm, cm->abc, 3.0, 1.0, &con);
  if((csq = esl_sq_CreateFrom("consensus", con->cseq, NULL, NULL, NULL)) == NULL)
    { status = eslEMEM; goto ERROR; }
  esl_sq_Digitize(cm->abc, csq);
  dpc = count_align_dp_calcs(cm, csq->n) / 1000000.;

  /* cyk inside (score only) */
  esl_stopwatch_Start(w);
  /*CYKDivideAndConquer(cm, csq->dsq, csq->n, 0, 1, csq->n, NULL, NULL, NULL);*/
  CYKInsideScore(cm, csq->dsq, csq->n, 0, 1, csq->n, 
		 NULL, NULL); /* don't do QDB mode */
  esl_stopwatch_Stop(w);
  t_dc = w->user;
  mc_s = dpc / t_dc;

  /* HMM banded */
  /* Emit N seqs, and align them, to get total time up to reasonable level,
   * and to average out tightness of bands */
  int N = esl_opt_GetInteger(go, "-N");
  seqs_to_aln_t *seqs_to_aln = NULL;
  ESL_SQ **sq = NULL;
  ESL_ALLOC(sq, sizeof(ESL_SQ *) * N);
  int L; 
  float L_avg = 0.; 
  int i;
  for(i = 0; i < N; i++)
    {
      EmitParsetree(cm, r, "seq", TRUE, NULL, &(sq[i]), &L);
      /*esl_sqio_Write(stdout, sq[i], eslSQFILE_FASTA);*/
      L_avg += L;
    }
  L_avg /= (float) N;
  cm->align_opts |= CM_ALIGN_HBANDED;
  esl_stopwatch_Start(w);
  seqs_to_aln = CreateSeqsToAlnFromSq(sq, N, FALSE);
  if((status = ActuallyAlignTargets(cm, errbuf, seqs_to_aln, NULL, NULL, 0, 0, TRUE, NULL)) != eslOK) goto ERROR;
  esl_stopwatch_Stop(w);
  t_hb = w->user / (float) N;
  FreeSeqsToAln(seqs_to_aln);

  printf("#\n");
  printf("#\t\t\t Alignment statistics:\n");
  printf("#\t\t\t %7s %6s %6s %8s %8s %8s\n",             "alg",     "Mc",     "L",     "Mc/s",    "s/seq",        "accel");
  printf("#\t\t\t %7s %6s %6s %8s %8s %8s\n",             "-------", "------", "------","--------", "--------", "--------");
  /*mc_s = dpc / t_dc; */
  printf(" \t\t\t %7s %6.1f %6d %8.1f %8.3f %8s\n",       "cyk",      dpc,      csq->n,  mc_s,      t_dc,       "-");
  printf(" \t\t\t %7s %6s %6.0f %8s %8.3f %8.2f\n",       "hb cyk",   "?",      L_avg,   "?",       t_hb,       (t_dc/t_hb));

  esl_sq_Destroy(csq);
  FreeCMConsensus(con);
  return eslOK;
 ERROR:
  esl_fatal("ERROR code %d in summarize_stats().", status);
  return status; /* NOTREACHED */
}

/* Function: count_align_dp_calcs()
 * Date:     EPN, Wed Aug 22 09:08:03 2007
 *
 * Purpose:  Count all non-d&c inside DP calcs for a CM 
 *           alignment of a seq of length L. Similar to smallcyk.c's
 *           CYKDemands() but takes into account number of
 *           transitions from each state, and is concerned
 *           with a scanning dp matrix, not an alignment matrix.
 *
 * Args:     cm     - the model
 *           L      - length of sequence
 *
 * Returns: (float) the total number of DP calculations.
 */
float count_align_dp_calcs(CM_t *cm, int L)
{
  int v, j;
  float dpcalcs = 0.;
  float dpcalcs_bif = 0.;
  
  float  dpcells     = 0.;
  float  dpcells_bif = 0.;

  dpcells = (L+2) * (L+1) * 0.5; /* fillable dp cells per state (deck) */
  for (j = 0; j <= L; j++)
    dpcells_bif += (j+2) * (j+1) * .5;
  dpcalcs_bif = CMCountStatetype(cm, B_st) * dpcells_bif; /* no choice of transitions */
  for(v = 0; v < cm->M; v++)
    if(cm->sttype[v] != B_st && cm->sttype[v] != E_st)
      dpcalcs += dpcells * cm->cnum[v]; /* cnum choices of transitions */

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



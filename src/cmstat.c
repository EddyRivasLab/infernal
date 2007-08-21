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
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "verbose output", 1 },
  { "--beta",    eslARG_REAL,   "1E-7",NULL, "0<x<1",   NULL,      NULL,        NULL, "set tail loss prob for QDB stats to <x>", 1 },
  { "--tau",     eslARG_REAL,   "1E-7",NULL, "0<x<1",   NULL,      NULL,       NULL,  "set tail loss prob for HMM banded stats to <x>", 1 },
  { "--exp",     eslARG_REAL,   NULL,  NULL, "0<x<=1.0",NULL,      NULL,        NULL, "exponentiate CM probabilities by <x> before calc'ing stats",  1 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <cmfile>";
static char banner[] = "display summary statistics for CMs";

static double cm_MeanMatchRelativeEntropy(const CM_t *cm);
static double cm_MeanMatchEntropy(const CM_t *cm);
static double cm_MeanMatchInfo(const CM_t *cm);

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go = NULL;   /* command line processing   */
  ESL_ALPHABET    *abc = NULL;  /* alphabet                  */
  char            *cmfile;	/* name of input CM file     */ 
  CMFILE          *cmfp;	/* open input CM file stream */
  CM_t            *cm;          /* CM most recently read     */
  int              ncm;         /* CM index                  */
  float            s_dpc;       /*     DP cells for search   */
  float            s_dpc_qdb;   /* QDB DP cells for search   */
  float            s_spdup;     /* search speedup estimate   */

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
    if(esl_opt_GetBoolean(go, "-v"))
      {
	/* estimate speed up due to QDB */
	s_dpc     = CYKDemands(cm, cm->W, NULL,     NULL,     TRUE) / 1000000.; 
	s_dpc_qdb = CYKDemands(cm, cm->W, cm->dmin, cm->dmax, TRUE) / 1000000.; 
	s_spdup   = s_dpc / s_dpc_qdb;
	s_dpc    *= 1000000. / cm->W;
	s_dpc_qdb*= 1000000. / cm->W;
	printf("\tSearch stats:\n\t\t     %6.2f MC/MB\n\t\tQDB: %6.2f MC/MB\n\t\tACC: %6.2f\n",
	       s_dpc, s_dpc_qdb, s_spdup);
      }
    FreeCM(cm);
  }
    
  esl_alphabet_Destroy(abc);
  CMFileClose(cmfp);
  esl_getopts_Destroy(go);
  exit(0);
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


/*
 * EPN, Wed Nov 22 12:00:20 2006
 * Ported from RSEARCH-1.1
 *
 * stats.c
 *
 * Routines for calculating statistics in Infernal.  
 *
 * Robert J. Klein
 * Separate file started May 17, 2002
 */
#include "esl_config.h"
#include "p7_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <ctype.h>
#include <string.h>


#include "easel.h"
#include "esl_exponential.h"
#include "esl_histogram.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_ssi.h"
#include "esl_vectorops.h"

#include "cm_dispatch.h"
#include "funcs.h"		/* external functions                   */
#include "structs.h"		/* data structures, macros, #define's   */

/* Function: AllocCMStats()
 * Date:     EPN, Wed May  2 14:10:25 2007      
 *
 * Purpose:  Allocate a CMStats_t data structure given
 *           the number of partitions and relevant info needed
 *           for creating the HMMFilter_t objects.
 * 
 * Args:     np - number of partitions 
 */
CMStats_t *
AllocCMStats(int np)
{
  int status;
  CMStats_t  *cmstats;
  int i, p;

  ESL_ALLOC(cmstats, sizeof(struct cmstats_s));

  cmstats->np = np;
  ESL_ALLOC(cmstats->ps, sizeof(int) * cmstats->np);
  ESL_ALLOC(cmstats->pe, sizeof(int) * cmstats->np);
  ESL_ALLOC(cmstats->expAA, sizeof(ExpInfo_t **) * EXP_NMODES);
  ESL_ALLOC(cmstats->hfiA,  sizeof(HMMFilterInfo_t *) * FTHR_NMODES);
  for(i = 0; i < EXP_NMODES; i++) {
    ESL_ALLOC(cmstats->expAA[i], sizeof(ExpInfo_t *) * cmstats->np);
    for(p = 0; p < cmstats->np; p++) {
      cmstats->expAA[i][p] = CreateExpInfo();
      if(cmstats->expAA[i][p] == NULL) goto ERROR; /* memory error */
    }
  }
  for(i = 0; i < FTHR_NMODES; i++) {
    cmstats->hfiA[i] = CreateHMMFilterInfo();
    if(cmstats->hfiA[i] == NULL) goto ERROR; /* memory error */
  }
  return cmstats;

 ERROR:
  cm_Fail("AllocCMStats() memory allocation error.");
  return NULL; /* never reached */
}

/* Function: FreeCMStats()
 * Returns: (void) 
 */
void 
FreeCMStats(CMStats_t *cmstats)
{
  int i, p;
  for(i = 0; i < EXP_NMODES; i++)
    {
      for(p = 0; p < cmstats->np; p++)
	if(cmstats->expAA[i][p] != NULL) free(cmstats->expAA[i][p]);
      free(cmstats->expAA[i]);
    }
  free(cmstats->expAA);
  for(i = 0; i < FTHR_NMODES; i++)
    FreeHMMFilterInfo(cmstats->hfiA[i]);
  free(cmstats->hfiA);
  free(cmstats->ps);
  free(cmstats->pe);
  free(cmstats);
}  


/* Function: debug_print_cmstats
 */
int debug_print_cmstats(CM_t *cm, char *errbuf, CMStats_t *cmstats, int has_fthr)
{
  int status;
  int p;
  char *namedashes;
  int ni;
  int namewidth = strlen(cm->name); 
  ESL_ALLOC(namedashes, sizeof(char) * namewidth+1);
  namedashes[namewidth] = '\0';
  for(ni = 0; ni < namewidth; ni++) namedashes[ni] = '-';

  printf("Num partitions: %d\n", cmstats->np);
  for (p = 0; p < cmstats->np; p++)
    {
      printf("Partition %d: start: %d end: %d\n", p, cmstats->ps[p], cmstats->pe[p]);
      printf("cm_lc  exp tail:\t");
      debug_print_expinfo(cmstats->expAA[EXP_CM_LC][p]);
      printf("cm_gc  exp tail:\t");
      debug_print_expinfo(cmstats->expAA[EXP_CM_GC][p]);
      printf("cm_li  exp tail:\t");
      debug_print_expinfo(cmstats->expAA[EXP_CM_LI][p]);
      printf("cm_gi  exp tail:\t");
      debug_print_expinfo(cmstats->expAA[EXP_CM_GI][p]);
      printf("cp9_lv exp tail:\t");
      debug_print_expinfo(cmstats->expAA[EXP_CP9_LV][p]);
      printf("cp9_gv exp tail:\t");
      debug_print_expinfo(cmstats->expAA[EXP_CP9_GV][p]);
      printf("cp9_lf exp tail:\t");
      debug_print_expinfo(cmstats->expAA[EXP_CP9_LF][p]);
      printf("cp9_gf exp tail:\t");
      debug_print_expinfo(cmstats->expAA[EXP_CP9_GF][p]);
      printf("\n\n");
    }

  if(has_fthr)
    {
      printf("Filter CM_LC info:\n");
      if((status = DumpHMMFilterInfo(stdout, cmstats->hfiA[FTHR_CM_LC], errbuf, cm, EXP_CM_LC, EXP_CP9_LF, cmstats->hfiA[FTHR_CM_LC]->dbsize, 1, namewidth, namedashes)) != eslOK) return status;
      printf("Filter CM_LI info:\n");
      if((status = DumpHMMFilterInfo(stdout, cmstats->hfiA[FTHR_CM_LI], errbuf, cm, EXP_CM_LI, EXP_CP9_LF, cmstats->hfiA[FTHR_CM_LI]->dbsize, 1, namewidth, namedashes)) != eslOK) return status;
      printf("Filter CM_GC info:\n");
      if((status = DumpHMMFilterInfo(stdout, cmstats->hfiA[FTHR_CM_GC], errbuf, cm, EXP_CM_GC, EXP_CP9_GF, cmstats->hfiA[FTHR_CM_GC]->dbsize, 1, namewidth, namedashes)) != eslOK) return status;
      printf("Filter CM_GI info:\n");
      if((status = DumpHMMFilterInfo(stdout, cmstats->hfiA[FTHR_CM_GI], errbuf, cm, EXP_CM_GI, EXP_CP9_GF, cmstats->hfiA[FTHR_CM_GI]->dbsize, 1, namewidth, namedashes)) != eslOK) return status;
      printf("\n\n");
    }
  free(namedashes);
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "Memory allocation error in debug_print_cmstats().");
  return status; /* NEVERREACHED */
}

/* Function: debug_print_expinfo
 */
int debug_print_expinfo(ExpInfo_t *exp)
{
  if(exp->is_valid) printf("cur_eff_dbsize: %ld lambda: %f mu_extrap: %f mu_orig: %f dbsize: %ld nrandhits: %d tailp: %f (valid)\n", exp->cur_eff_dbsize, exp->lambda, exp->mu_extrap, exp->mu_orig, exp->dbsize, exp->nrandhits, exp->tailp);
  else              printf("invalid (not yet set)\n");
  fflush(stdout);
  return eslOK;
}

/* Function: get_gc_comp
 * Date:     EPN, Tue Aug  7 08:51:06 2007
 * Purpose:  Given a sequence and start and stop coordinates, returns 
 *           integerized GC composition of the region. This Easelfied
 *           version replaces RSEARCH's version.
 */
int get_gc_comp(const ESL_ALPHABET *abc, ESL_DSQ *dsq, int start, int stop) 
{
  int status;
  int i;
  float *ct = NULL;
  float  gc;

  /* contract check */
  if(abc == NULL) cm_Fail("get_gc_comp alphabet is NULL.");
  if(dsq == NULL) cm_Fail("get_gc_comp dsq is NULL.");
  if(abc->type != eslRNA && abc->type != eslDNA)  cm_Fail("get_gc_comp expects alphabet of RNA or DNA");

  if (start > stop) {
    i = start;
    start = stop;
    stop = i;
  }
  ESL_ALLOC(ct, sizeof(float) * abc->K);
  esl_vec_FSet(ct, abc->K, 0.);
  for (i = start; i <= stop; i++)
  {
    if(esl_abc_XIsGap(abc, dsq[i])) cm_Fail("in get_gc_comp, res %d is a gap!\n", i);
    esl_abc_FCount(abc, ct, dsq[i], 1.);
  }
  gc = ct[abc->inmap[(int) 'G']] + ct[abc->inmap[(int) 'C']];
  gc /= ((float) (stop-start+1));
  gc *= 100.;
  free(ct);
  return (int) (gc + 0.5);

 ERROR:
  cm_Fail("Memory allocation error.");
  return 0; /* never reached */
}


/* Function: get_alphabet_comp
 * Date:     EPN, Wed May  7 18:39:28 2008
 * Purpose:  Given a sequence and start and stop coordinates, allocate and fill 
 *           a abc->K sized vector with frequency of each residue.
 */
int get_alphabet_comp(const ESL_ALPHABET *abc, ESL_DSQ *dsq, int start, int stop, float **ret_freq) 
{
  int status;
  int i;
  float *freq;

  ESL_ALLOC(freq, sizeof(float) * abc->K);
  esl_vec_FSet(freq, abc->K, 0.0);

  /* contract check */
  if(abc == NULL) cm_Fail("get_alphabet_comp alphabet is NULL.");
  if(dsq == NULL) cm_Fail("get_alphabet_comp dsq is NULL.");
  if(abc->type != eslRNA && abc->type != eslDNA)  cm_Fail("get_alphabet_comp expects alphabet of RNA or DNA");
  if(ret_freq == NULL) cm_Fail("get_alphabet_comp ret_freq is NULL.");

  if (start > stop) {
    i = start;
    start = stop;
    stop = i;
  }
  for (i = start; i <= stop; i++)
  {
    if(esl_abc_XIsGap(abc, dsq[i])) cm_Fail("in get_alphabet_comp, res %d is a gap!\n", i);
    esl_abc_FCount(abc, freq, dsq[i], 1.);
  }
  esl_vec_FNorm(freq, abc->K);
  *ret_freq = freq;
  return eslOK;

 ERROR:
  cm_Fail("get_alphabet_comp() memory allocation error.");
  return 0; /* never reached */
}

/* Function: GetDBSize()
 *
 * Date:     EPN, Wed Apr  2 16:33:19 2008
 *
 * Purpose:  Given a sequence file ptr, determine the number, and summed 
 *           size of the seqs in the file. 
 *
 * Args:     sqfp          - open sequence file
 *           errbuf        - for writing error messages
 *           ret_N         - RETURN: total length (residues) or all seqs in sqfp
 *           ret_nseq      - RETURN: number of seqs in sqfp
 *           ret_namewidth - RETURN: max length of name in sqfp
 *
 * Returns:  eslOK on success, other status on failure, errbuf filled with error message.
 */
int GetDBSize (ESL_SQFILE *sqfp, char *errbuf, long *ret_N, int *ret_nseq, int *ret_namewidth)
{
  int     status;
  ESL_SQ *sq;
  long    N = 0;
  int     namewidth = 11; /* length of "target name" */
  int     nseq = 0;

  sq = esl_sq_Create(); 

  while ((status = esl_sqio_ReadInfo(sqfp, sq)) == eslOK) {
    N += sq->L;
    namewidth = ESL_MAX(namewidth, strlen(sq->name));
    esl_sq_Reuse(sq); 
    nseq++;
  } 
  if (status != eslEOF) { 
    ESL_FAIL(status, errbuf, "Parse failed, file %s:\n%s", 
	     sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
  }
  esl_sq_Destroy(sq); 
  esl_sqfile_Position(sqfp, (off_t) 0); /* rewind sequence file to beginning */

  if(ret_N != NULL)          *ret_N         = N;
  if(ret_nseq != NULL)       *ret_nseq      = nseq;
  if(ret_namewidth != NULL)  *ret_namewidth = namewidth;
  return eslOK;
}


/* Function: GetDBInfo()
 *
 * Date:     Easelification: EPN, Thu Dec  7 06:07:58 2006
 *           (initial - RSEARCH::get_dbinfo()) RJK, Thu Apr 11, 2002 [St. Louis]
 *
 * Purpose:  Given a sequence file name, determine the total size of the
 *           seqs in the file (DB) and GC content information.
 *
 * Args:     abc       - alphabet for sq file
 *           sqfp      - open sequence file
 *           errbuf    - for error messages
 *           ret_N     - RETURN: total length (residues) or all seqs in seqfile
 *           ret_nseq  - RETURN: number of seqs in seqfile
 *           ret_gc_ct - RETURN: gc_ct[x] observed 100-nt segments with GC% of x [0..100] 
 *           
 * Returns:  eslOK on success, other status on failure, errbuf filled with error message.
 */
int GetDBInfo (const ESL_ALPHABET *abc, ESL_SQFILE *sqfp, char *errbuf, long *ret_N, int *ret_nseq, double **ret_gc_ct)
{
  int               status;
  ESL_SQ           *sq;
  int               i, j, jp;  
  long              N = 0;
  double           *gc_ct;
  int               gc;
  int               nseq = 0;
  int               all_ambig_flag; /* used to check if curr DB chunk is all ambiguous characters 
				     * usually Ns, if it is, we don't count it towards the GC content info */
  if (abc       == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "GetDBInfo(), abc is NULL\n");
  if (ret_gc_ct == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "GetDBInfo(), ret_gc_ct is NULL\n");

  ESL_ALLOC(gc_ct, sizeof(double) * GC_SEGMENTS);
  for (i=0; i<GC_SEGMENTS; i++) gc_ct[i] = 0.;

  sq = esl_sq_CreateDigital(abc); 
  esl_sqfile_SetDigital(sqfp, abc);
  
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK) { 
    nseq++;
    N += sq->n;
    /*printf("new N: %d\n", N);*/
    for(i = 1; i <= sq->n; i += 100) {
      j = (i+99 <= sq->n) ? i+99 : sq->n;
      gc = get_gc_comp(abc, sq->dsq, i, j);
      all_ambig_flag = TRUE;
      for(jp = 0; jp < 100 && (jp+i) < sq->n; jp++) {
	if(sq->dsq[i+jp] < abc->K) {
	  all_ambig_flag = FALSE; 
	  break; 
	}
      }
      /*printf("N: %d i: %d gc: %d\n", N, i, gc);*/
      /* don't count GC content of chunks < 20 nt, very hacky;
       * don't count GC content of chunks that are all N, this
       * will be common in RepeatMasked genomes where poly-Ns could
       * skew the base composition stats of the genome */
      if(j > 20 && !all_ambig_flag)
	{
	  /*printf("j: %d i: %d adding 1 to gc_ct[%d]\n", j, i, ((int) gc));*/
	  gc_ct[(int) gc] += 1.;
	}
    }
    esl_sq_Reuse(sq); 
  } 
  if (status != eslEOF) 
    ESL_FAIL(status, errbuf, "Parse failed, file %s:\n%s",
	     sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
  esl_sq_Destroy(sq); 
  esl_sqfile_Position(sqfp, (off_t) 0);

  if(ret_N != NULL)      *ret_N     = N;
  if(ret_nseq != NULL)   *ret_nseq  = nseq;
  *ret_gc_ct = gc_ct;

  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "GetDBInfo(): memory allocation error.");
}

/* Function: E2MinScore()
 * Date:     EPN, Tue Dec 11 15:40:25 2007 
 *           (morphed from RSEARCH: RJK, Mon April 15 2002 [St. Louis])
 *
 * Purpose:  Given an E-value <E> and a CM with valid exp tail stats 
 *           determine the minimal bit score that will give an E-value 
 *           of <E> across all partitions. This will be a safe bit score
 *           cutoff to use when returning hits in DispatchSearch().
 *           Because no score less than this will have an E-value 
 *           less than E.
 *
 * Returns:  eslOK on success, <ret_sc> filled with bit score
 *           error code on failure, errbuf filled with message
 */
int E2MinScore (CM_t *cm, char *errbuf, int exp_mode, float E, float *ret_sc)
{
  int status;
  float low_sc, sc;
  int p;

  /* contract checks */
  if(!(cm->flags & CMH_EXPTAIL_STATS)) ESL_FAIL(eslEINCOMPAT, errbuf, "E2MinScore, CM's CMH_EXPTAIL_STATS flag is down.");
  if(ret_sc == NULL)                   ESL_FAIL(eslEINCOMPAT, errbuf, "E2MinScore, ret_sc is NULL");

  if((status = E2ScoreGivenExpInfo(cm->stats->expAA[exp_mode][0], errbuf, E, &low_sc)) != eslOK) return status;
  for(p = 1; p < cm->stats->np; p++) {
    if((status = E2ScoreGivenExpInfo(cm->stats->expAA[exp_mode][p], errbuf, E, &sc)) != eslOK) return status;
    low_sc = ESL_MIN(low_sc, sc);
  }
  *ret_sc = low_sc;
  return eslOK;
}

/* Function: E2ScoreGivenExpInfo()
 * Date:     EPN, Thu Apr  3 15:57:34 2008
 *
 * Purpose:  Given an E-value <E> and a exp tail stat structure
 *           determine the bit score that will give an E-value 
 *           of <E>.
 *
 * Returns:  eslOK on success, <ret_sc> filled with bit score
 *           error code on failure, errbuf filled with message
 */
int E2ScoreGivenExpInfo(ExpInfo_t *exp, char *errbuf, float E, float *ret_sc)
{
  float sc;
  /* contract checks */
  if(ret_sc == NULL)                   ESL_FAIL(eslEINCOMPAT, errbuf, "E2ScoreGivenExpInfo, ret_sc is NULL");
  sc  = exp->mu_extrap + ((log(E/exp->cur_eff_dbsize)) / (-1 * exp->lambda));
  *ret_sc = sc;
  return eslOK;
}


/* Function: Score2MaxE()
 * Date:     EPN, Wed Jan 16 14:25:44 2008
 *
 * Purpose:  Given a bit score <sc> and a CM with valid exp tail stats
 *           determine the maximal E-value that will be assigned to 
 *           bit score of <sc> across all partitions. 
 *           This will be a 'safe' E-value cutoff to use 
 *           to always return hits with bit score of <sc> or greater.
 *           because no E-value higher than this will be assigned to
 *           a bit score greater than <sc>.
 *
 * Returns:  eslOK on success, <ret_E> filled with E value
 *           error code on failure, errbuf filled with message
 */
int Score2MaxE (CM_t *cm, char *errbuf, int exp_mode, float sc, float *ret_E)
{
  float high_E, E;
  int p;

  /* contract checks */
  if(!(cm->flags & CMH_EXPTAIL_STATS)) ESL_FAIL(eslEINCOMPAT, errbuf, "Score2E(), CM's CMH_EXPTAIL_STATS flag is down.");
  if(ret_E == NULL)                   ESL_FAIL(eslEINCOMPAT, errbuf, "Score2E(), ret_E is NULL");

  high_E = Score2E(sc, cm->stats->expAA[exp_mode][0]->mu_extrap, cm->stats->expAA[exp_mode][0]->lambda, cm->stats->expAA[exp_mode][0]->cur_eff_dbsize);

  if(! cm->stats->expAA[exp_mode][0]->is_valid) ESL_FAIL(eslEINCOMPAT, errbuf, "Score2E(), CM's exp tail stats for mode: %d partition: %d are invalid.", exp_mode, p);
  for(p = 1; p < cm->stats->np; p++) {
    if(! cm->stats->expAA[exp_mode][p]->is_valid) ESL_FAIL(eslEINCOMPAT, errbuf, "Score2E(), CM's exp tail stats for mode: %d partition: %d are invalid.", exp_mode, p);
    E = Score2E(sc, cm->stats->expAA[exp_mode][p]->mu_extrap, cm->stats->expAA[exp_mode][p]->lambda, cm->stats->expAA[exp_mode][p]->cur_eff_dbsize);
    high_E = ESL_MAX(high_E, E);
  }
  *ret_E = high_E;
  return eslOK;
}


/* Function: Score2E()
 * Date:     EPN, Fri Feb 15 12:40:21 2008
 *
 * Purpose:  Given a bit score <x>, a mu and lambda that describe
 *           an exponential tail distribution, and an effective 
 *           database size, return the E-value of <sc>.
 *
 * Returns:  E value of <x>
 */
double Score2E (float x, double mu, double lambda, long eff_dbsize)
{
  return esl_exp_surv(x, mu, lambda) * (double) eff_dbsize;
}

/* Function: CM2ExpMode
 * Date:     EPN, Mon May  7 17:43:28 2007
 * Purpose:  Return the exp_mode for the CM and HMM
 *           given the flags and search options in the
 *           CM data structure.
 */
int CM2ExpMode(CM_t *cm, int search_opts, int *ret_cm_exp_mode, int *ret_cp9_exp_mode)
{
  int cm_exp_mode;
  int cp9_exp_mode;

  /* check contract */
  if(!(cm->flags & CMH_CP9) || cm->cp9 == NULL) cm_Fail("ERROR no CP9 in CM2ExpMode()\n");
  if(cm->flags & CMH_LOCAL_BEGIN) {
    if(search_opts & CM_SEARCH_INSIDE) cm_exp_mode = EXP_CM_LI;
    else               	               cm_exp_mode = EXP_CM_LC;
  }
  else {
    if(search_opts & CM_SEARCH_INSIDE) cm_exp_mode = EXP_CM_GI;
    else        	               cm_exp_mode = EXP_CM_GC;
  }

  if(cm->cp9->flags & CPLAN9_LOCAL_BEGIN) {
    if(search_opts & CM_SEARCH_HMMFORWARD) cp9_exp_mode = EXP_CP9_LF;
    else                                   cp9_exp_mode = EXP_CP9_LV;
  }
  else {
    if(search_opts & CM_SEARCH_HMMFORWARD) cp9_exp_mode = EXP_CP9_GF;
    else                                   cp9_exp_mode = EXP_CP9_GV;
  }

  if(ret_cm_exp_mode  != NULL) *ret_cm_exp_mode  = cm_exp_mode;
  if(ret_cp9_exp_mode != NULL) *ret_cp9_exp_mode = cp9_exp_mode;
  return eslOK;
}


/*
 * Function: CM2FthrMode
 * Date:     EPN, Tue Dec 11 13:16:35 2007
 * Purpose:  Return the filter threshold mode for the CM 
 *           given CM's flags and a passed in search options
 *           int.
 * 
 * Returns: eslOK on success, eslEINCOMPAT if search_opts indicate
 *          we're doing HMM search, errbuf is filled with error message.
 */
int CM2FthrMode(CM_t *cm, char *errbuf, int search_opts, int *ret_fthr_mode)
{
  int fthr_mode;

  /* check contract */
  if(search_opts & CM_SEARCH_HMMVITERBI) ESL_FAIL(eslEINCOMPAT, errbuf, "CM2FThrMode(), search_opts CM_SEARCH_HMMVITERBI flag raised.\n");
  if(search_opts & CM_SEARCH_HMMFORWARD) ESL_FAIL(eslEINCOMPAT, errbuf, "CM2FThrMode(), search_opts CM_SEARCH_HMMFORWARD flag raised.\n");
  if(ret_fthr_mode == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "CM2FThrMode(), ret_fthr_mode is NULL.");

  if(cm->flags & CMH_LOCAL_BEGIN) {
    if(search_opts & CM_SEARCH_INSIDE) fthr_mode = FTHR_CM_LI;
    else               	               fthr_mode = FTHR_CM_LC;
  }
  else {
    if(search_opts & CM_SEARCH_INSIDE) fthr_mode = FTHR_CM_GI;
    else        	               fthr_mode = FTHR_CM_GC;
  }
  if(ret_fthr_mode  != NULL) *ret_fthr_mode  = fthr_mode;
  return eslOK;
}

/* Function: ExpModeIsLocal()
 * Date:     EPN, Mon Dec 10 09:07:59 2007
 * Purpose:  Given a exp tail mode, return TRUE if it corresponds to 
 *           a local model configuration.
 *
 * Args:     exp_mode     - the mode 0..EXP_NMODES-1
 */
int
ExpModeIsLocal(int exp_mode)
{
  ESL_DASSERT1((exp_mode >= 0 && exp_mode < EXP_NMODES));

  switch (exp_mode) {
  case EXP_CM_LC: 
  case EXP_CM_LI: 
  case EXP_CP9_LV: 
  case EXP_CP9_LF: 
    return TRUE;
    break;
  case EXP_CM_GC:
  case EXP_CM_GI: 
  case EXP_CP9_GV: 
  case EXP_CP9_GF: 
    return FALSE;
    break;
  default: 
    cm_Fail("ExpModeIsLocal(): bogus exp_mode: %d\n", exp_mode);
  }
  return FALSE; /* never reached */
}


/* Function: ExpModeIsForCM()
 * Date:     EPN, Mon Dec 10 09:11:55 2007
 * Purpose:  Given a exp tail mode, return TRUE if it corresponds to 
 *           a CM (not a CP9 HMM).
 *
 * Args:     exp_mode     - the mode 0..EXP_NMODES-1
 */
int
ExpModeIsForCM(int exp_mode)
{
  ESL_DASSERT1((exp_mode >= 0 && exp_mode < EXP_NMODES));

  switch (exp_mode) {
  case EXP_CM_LC: 
  case EXP_CM_LI: 
  case EXP_CM_GC:
  case EXP_CM_GI: 
    return TRUE;
    break;
  case EXP_CP9_LV: 
  case EXP_CP9_LF: 
  case EXP_CP9_GV: 
  case EXP_CP9_GF: 
    return FALSE;
    break;
  default: 
    cm_Fail("ExpModeIsForCM(): bogus exp_mode: %d\n", exp_mode);
  }
  return FALSE; /* never reached */
}

/* Function: ExpModeToFthrMode()
 * Date:     EPN, Mon Dec 10 09:31:48 2007
 * Purpose:  Given a exp tail mode, return it's corresponding
 *           filter threshold mode, or -1 if there is none
 *
 * Args:     exp_mode     - the mode 0..EXP_NMODES-1
 */
int
ExpModeToFthrMode(int exp_mode)
{
  ESL_DASSERT1((exp_mode >= 0 && exp_mode < EXP_NMODES));

  switch (exp_mode) {
  case EXP_CM_GC: 
    return FTHR_CM_GC;
    break;
  case EXP_CM_GI: 
    return FTHR_CM_GI;
    break;
  case EXP_CM_LC: 
    return FTHR_CM_LC;
    break;
  case EXP_CM_LI: 
    return FTHR_CM_LI;
    break;
  case EXP_CP9_LV: 
  case EXP_CP9_LF: 
  case EXP_CP9_GV: 
  case EXP_CP9_GF: 
    return -1;
    break;
  default: 
    cm_Fail("ExpModeToFthrMode(): bogus exp_mode: %d\n", exp_mode);
  }
  return FALSE; /* never reached */
}

/* Function: CreateExpInfo()
 * Date:     EPN, Tue Dec 11 05:25:06 2007
 *
 * Purpose:  Allocate and minimally initialize a exp tail info object.
 *            
 * Returns:  Newly allocated ExpInfo_t object on success, NULL if some error occurs
 */
ExpInfo_t *
CreateExpInfo()
{
  int status;

  ExpInfo_t *exp = NULL;
  ESL_ALLOC(exp, sizeof(ExpInfo_t));

  exp->cur_eff_dbsize = 0;
  exp->lambda         = 0.;
  exp->mu_extrap      = 0.;
  exp->mu_orig        = 0.;
  exp->dbsize         = 0;
  exp->nrandhits      = 0;
  exp->tailp          = 0.;

  exp->is_valid = FALSE;
  return exp;

 ERROR:
  return NULL; /* reached if memory error */
}  


/* Function: SetExpInfo()
 * Date:     EPN, Wed Dec 12 13:43:36 2007
 *
 * Purpose:  Set parameters of a exp tail info object and raise it's is_valid 'flag'.
 *            
 * Returns:  void
 */
void 
SetExpInfo(ExpInfo_t *exp, double lambda, double mu_orig, long dbsize, int nrandhits, double tailp)
{
  exp->lambda    = lambda;
  exp->mu_orig   = mu_orig;
  exp->dbsize    = dbsize;
  exp->nrandhits = nrandhits;
  exp->tailp     = tailp;

  exp->mu_extrap = exp->mu_orig - log(1./exp->tailp) / exp->lambda;

  /* initialize exp->cur_eff_dbsize as effective database size as exp->nrandhits, this is
   * effective database size if we searched a database of size <exp->dbsize>, which we 
   * just calibrated for. we'll update this in cmsearch for the target database size. */
  exp->cur_eff_dbsize = (long) exp->nrandhits;

  exp->is_valid  = TRUE; /* we can now write Exp Info to a cm file */
  return;
}  

/* Function: DuplicateExpInfo()
 * Date:     EPN, Tue Dec 11 05:28:13 2007
 *
 * Purpose:  Duplicate a exp tail info object.
 *            
 * Returns:  Newly allocated ExpInfo_t object on success, NULL if some error occurs
 */
ExpInfo_t *
DuplicateExpInfo(ExpInfo_t *src)
{
  int status;

  ExpInfo_t *dest = NULL;
  ESL_ALLOC(dest, sizeof(ExpInfo_t));

  dest->cur_eff_dbsize = src->cur_eff_dbsize;
  dest->lambda         = src->lambda;
  dest->mu_orig        = src->mu_orig;
  dest->mu_extrap      = src->mu_extrap;
  dest->dbsize         = src->dbsize;
  dest->nrandhits      = src->nrandhits;
  dest->tailp          = src->tailp;
  dest->is_valid       = src->is_valid;
  return dest;

 ERROR:
  return NULL; /* reached if memory error */
}  


/* Function:  DescribeExpMode()
 * Incept:    EPN, Mon Jan  7 18:04:31 2008
 *
 * Purpose:   Returns the Exp Tail mode in text.
 *            For example, <DescribeExpMode(EXP_CM_GC)>
 *            returns "glocal CM  CYK".
 */
char *
DescribeExpMode(int exp_mode)
{
  switch (exp_mode) {
  case EXP_CP9_GV: return "hmm  glc  vit";
  case EXP_CP9_GF: return "hmm  glc  fwd";
  case EXP_CM_GC:  return " cm  glc  cyk";
  case EXP_CM_GI:  return " cm  glc  ins";
  case EXP_CP9_LV: return "hmm  loc  vit";
  case EXP_CP9_LF: return "hmm  loc  fwd";
  case EXP_CM_LC:  return " cm  loc  cyk";
  case EXP_CM_LI:  return " cm  loc  ins";
  default:     return "?";
  }
}


/* Function:  DescribeFthrMode()
 * Incept:    EPN, Mon Jan  7 18:41:41 2008
 *
 * Purpose:   Returns the Filter thresold mode in text.
 *            For example, <DescribeFThrMode(EXP_CM_GC)>
 *            returns "glocal CM CYK".
 */
char *
DescribeFthrMode(int fthr_mode)
{
  switch (fthr_mode) {
  case FTHR_CM_GC:  return "glc  cyk";
  case FTHR_CM_GI:  return "glc  ins";
  case FTHR_CM_LC:  return "loc  cyk";
  case FTHR_CM_LI:  return "loc  ins";
  default:     return "?";
  }
}

/* Function: UpdateExpsForDBSize()
 * Date:     EPN, Thu Jan 17 09:38:06 2008
 *
 * Purpose:  Update the <cur_eff_dbsize> parameter of the
 *           ExpInfo_t objects in a CM's cm->stats object 
 *           to reflect a database of size <dbsize>.
 *            
 * Returns:  eslOK on success, other Easel status code on contract 
 *           violation with informative error message in errbuf.
 */
int 
UpdateExpsForDBSize(CM_t *cm, char *errbuf, long dbsize)
{
  int i, p;
  /* contract checks */
  if(! (cm->flags & CMH_EXPTAIL_STATS)) ESL_FAIL(eslEINCOMPAT, errbuf, "UpdateExpsForDBSize(), cm does not have Exp stats.");

  for(i = 0; i < EXP_NMODES; i++) { 
    for(p = 0; p < cm->stats->np; p++) {
      cm->stats->expAA[i][p]->cur_eff_dbsize = (long) ((((double) dbsize / (double) cm->stats->expAA[i][p]->dbsize) * 
							((double) cm->stats->expAA[i][p]->nrandhits)) + 0.5);
    }
  }
  return eslOK;
}  

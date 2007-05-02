/*/************************************************************
 * @LICENSE@
 ************************************************************/
/* cmcalibrate.c
 * Score a CM and a CM Plan 9 HMM against random sequence 
 * data to set the statistical parameters for E-value determination,
 * and HMM filtering thresholds. 
 * 
 * EPN, Wed May  2 07:02:52 2007
 * based on HMMER-2.3.2's hmmcalibrate.c from SRE
 *  
 */
#include "config.h"	
#include "esl_config.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "structs.h"
#include "funcs.h"		/* external functions                   */
#include "stats.h"              /* EVD functions */
#include "cm_dispatch.h"	
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_histogram.h"
#include "esl_vectorops.h"
#include "esl_dmatrix.h"
#include "esl_ratematrix.h"
#include "esl_exponential.h"
#include "esl_gumbel.h"

#define ALGORITHMS "--cyk,--inside"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "show brief help on version and usage",   1 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "cmcalibrate [-options]";

int
main(int argc, char **argv)
{
  int status;			       /* status of a function call               */
  ESL_ALPHABET    *abc     = NULL;     /* sequence alphabet                       */
  ESL_GETOPTS     *go	   = NULL;     /* command line processing                 */
  ESL_RANDOMNESS  *r       = NULL;     /* source of randomness                    */
  CM_t            *cm      = NULL;     /* the covariance model                    */
  int              sc;		       /* a CYK or Inside score                   */
  int              nseq;	       /* counter over sequences                  */
  char             errbuf[eslERRBUFSIZE];
  char            *cmfile;             /* file to read CM(s) from                 */
  CMFILE          *cmfp;        /* open CM file for reading                 */
  enum { DO_CYK, DO_INSIDE } algorithm_choice;

  int  *partitions;             /* partition each GC % point seg goes to    */
  int   num_partitions = 1;     /* number of partitions                     */
  int   cm_num_samples;         /* # samples used to calculate  CM EVDs     */
  int   hmm_num_samples;        /* # samples used to calculate HMM EVDs     */
  int   sample_length;          /* sample len used for calc'ing stats (2*W) */
  int   do_qdb = TRUE;
  int   do_cp9_stats = TRUE;
  /*****************************************************************
   * Parse the command line
   *****************************************************************/
  go = esl_getopts_Create(options, usage);
  esl_opt_ProcessCmdline(go, argc, argv);
  esl_opt_VerifyConfig(go);
  if (esl_opt_IsSet(go, "-h")) {
    puts(usage);
    return eslOK;
  }

  if (esl_opt_ArgNumber(go) != 1) {
    puts("Incorrect number of command line arguments.");
    puts(usage);
    return eslFAIL;
  }
  cmfile = esl_opt_GetCmdlineArg(go, eslARG_STRING, NULL); /* NULL=no range checking */

  /*****************************************************************
   * Initializations, including opening and reading the HMM file 
   *****************************************************************/
  /*if ((r = esl_randomness_CreateTimeseeded()) == NULL)*/
  if ((r = esl_randomness_Create(33)) == NULL)
    esl_fatal("Failed to create random number generator: probably out of memory");

  /* currently set up for a single CM - temporary */
  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL)
    Die("Failed to open covariance model save file %s\n%s\n", cmfile, usage);
  if (! CMFileRead(cmfp, &cm))
    Die("Failed to read a CM from %s -- file corrupt?\n", cmfile);
  if (cm == NULL) 
    Die("%s empty?\n", cmfile);
  CMFileClose(cmfp);
  if(do_qdb)          cm->config_opts |= CM_CONFIG_QDB;
  if(do_cp9_stats)    cm->search_opts |= CM_SEARCH_CP9STATS;

  ConfigCM(cm, NULL, NULL);

  /****************************************************************
   * Get distribution of GC content from a long random sequence
   *****************************************************************/
  /* human_chr1.gc.code from ~/notebook/7_0502_inf_cmcalibrate/gc_distros/ */
  /* Replace with RFAMSEQ derived one */
  int *gc_ct = MallocOrDie(sizeof(int) * GC_SEGMENTS); /* this should really be doubles */
  gc_ct[0] = 79; 
  gc_ct[1] = 59; gc_ct[2] = 56; gc_ct[3] = 79; gc_ct[4] = 91; gc_ct[5] = 89;
  gc_ct[6] = 114; gc_ct[7] = 105; gc_ct[8] = 152; gc_ct[9] = 149; gc_ct[10] = 181;
  gc_ct[11] = 203; gc_ct[12] = 293; gc_ct[13] = 401; gc_ct[14] = 586; gc_ct[15] = 894;
  gc_ct[16] = 1256; gc_ct[17] = 1909; gc_ct[18] = 2901; gc_ct[19] = 4333; gc_ct[20] = 6421;
  gc_ct[21] = 9027; gc_ct[22] = 12494; gc_ct[23] = 16458; gc_ct[24] = 20906; gc_ct[25] = 26706;
  gc_ct[26] = 32624; gc_ct[27] = 38814; gc_ct[28] = 44862; gc_ct[29] = 51035; gc_ct[30] = 57425;
  gc_ct[31] = 62700; gc_ct[32] = 67732; gc_ct[33] = 72385; gc_ct[34] = 75656; gc_ct[35] = 78714;
  gc_ct[36] = 80693; gc_ct[37] = 81741; gc_ct[38] = 81731; gc_ct[39] = 81291; gc_ct[40] = 79600;
  gc_ct[41] = 76883; gc_ct[42] = 75118; gc_ct[43] = 71834; gc_ct[44] = 69564; gc_ct[45] = 68058;
  gc_ct[46] = 66882; gc_ct[47] = 64922; gc_ct[48] = 63846; gc_ct[49] = 61481; gc_ct[50] = 57994;
  gc_ct[51] = 54322; gc_ct[52] = 50144; gc_ct[53] = 46483; gc_ct[54] = 43104; gc_ct[55] = 40144;
  gc_ct[56] = 36987; gc_ct[57] = 33653; gc_ct[58] = 29844; gc_ct[59] = 25838; gc_ct[60] = 21861;
  gc_ct[61] = 17858; gc_ct[62] = 14685; gc_ct[63] = 12233; gc_ct[64] = 9981; gc_ct[65] = 8047;
  gc_ct[66] = 6495; gc_ct[67] = 5263; gc_ct[68] = 4210; gc_ct[69] = 3405; gc_ct[70] = 2653;
  gc_ct[71] = 2157; gc_ct[72] = 1786; gc_ct[73] = 1458; gc_ct[74] = 1230; gc_ct[75] = 1048;
  gc_ct[76] = 915; gc_ct[77] = 802; gc_ct[78] = 769; gc_ct[79] = 668; gc_ct[80] = 544;
  gc_ct[81] = 462; gc_ct[82] = 395; gc_ct[83] = 329; gc_ct[84] = 228; gc_ct[85] = 170;
  gc_ct[86] = 116; gc_ct[87] = 94; gc_ct[88] = 39; gc_ct[89] = 46; gc_ct[90] = 22;
  gc_ct[91] = 11; gc_ct[92] = 5; gc_ct[93] = 3; gc_ct[94] = 0; gc_ct[95] = 1;
  gc_ct[96] = 0; gc_ct[97] = 1; gc_ct[98] = 0; gc_ct[99] = 0; gc_ct[100] = 0;
  /* distro is skewed towards low GC, you can see it, thanks to fixed width font */

  /****************************************************************
   * Determine CM and CP9 EVDs 
   *****************************************************************/
  /* temporary: no partitions */
  partitions = MallocOrDie(sizeof(int) * GC_SEGMENTS+1);
  int i;
  int N; /* N is the database size (2 * length) */
  N = 2000000;
  for (i = 0; i < GC_SEGMENTS; i++) 
    partitions[i] = 0;

  sample_length = 2.0 * cm->W;
  cm_num_samples       = 1000;
  printf("%-40s ... ", "Determining CM  EVD"); fflush(stdout);
  //serial_make_histogram (gc_ct, partitions, num_partitions,
  //			 cm, cm_num_samples, sample_length, 
  //			 FALSE, /* we're not doing CP9 stats */
  //			 TRUE); /* use easel, discard eventually */
  for (i=0; i<GC_SEGMENTS; i++) 
    {
      cm->lambda[i] = 0.509238; //comment me
      cm->K[i] = .014367;    //comment me
      cm->mu[i] = log(cm->K[i]*N)/cm->lambda[i];
    }
  printf("done. [lambda:%f mu:%f K: %f]\n", cm->lambda[50], cm->mu[50], cm->K[50]);

  hmm_num_samples       = 5000;
  printf("%-40s ... ", "Determining CP9 EVD"); fflush(stdout);
  serial_make_histogram (gc_ct, partitions, num_partitions,
			 cm, cm_num_samples, sample_length, 
			 TRUE,  /* we are doing CP9 stats */
			 TRUE); /* use easel, discard eventually */
  for (i=0; i<GC_SEGMENTS; i++) 
    cm->cp9_mu[i] = log(cm->cp9_K[i]*N)/cm->cp9_lambda[i];
  printf("done. [lambda:%f mu:%f K: %f]\n", cm->cp9_lambda[50], cm->cp9_mu[50], cm->K[50]);

  /****************************************************************
   * Determine CP9 filtering thresholds
   *****************************************************************/
  float globE = 50.;
  double globP = 0.01;
  double globT;
  /* determine CM bit score that gives E value of globE */
  globT = esl_gumbel_invcdf(globP, cm->mu[50], cm->lambda[50]);
  printf("bit score that gives P value of %f: %f\n", globP, globT);
  printf("E value of score: %f is %f\n", globT, RJK_ExtremeValueE(globT, cm->mu[50], cm->lambda[50]));

  /* First, glocal CP9 to glocal CM */
  float fraction = 0.95;
  float thresh;
  thresh = 
    FindHMMFilterThreshold(cm, FALSE,  /* glocal CM  */
			   FALSE,      /* glocal CP9 */
			   fraction, globT, 1000);
  printf("Glocal glocal thresh: %.3f\n", thresh);
  thresh = 
    FindHMMFilterThreshold(cm, FALSE, /* glocal CM  */
			   TRUE,      /* local CP9 */
			   fraction, globT, 1000);
  printf("Glocal local thresh: %.3f\n", thresh);
  thresh = 
    FindHMMFilterThreshold(cm, TRUE,   /* local CM  */
			   FALSE,      /* glocal CP9 */
			   fraction, globT, 1000);
  printf("Local glocal thresh: %.3f\n", thresh);
  thresh = 
    FindHMMFilterThreshold(cm, TRUE,  /* local CM  */
			   TRUE,      /* local CP9 */
			   fraction, globT, 1000);
  printf("Local local thresh: %.3f\n", thresh);


  /* end of from cmemit.c */

  FreeCM(cm);
  return eslOK;

 ERROR:
  return status;
}


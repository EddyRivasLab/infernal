/************************************************************
 * @LICENSE@
 ************************************************************/

/* cmemit.c
 * EPN, 09.01.06 Janelia Farm
 * based on HMMER-2.3.2's hmmemit.c from SRE
 * Easelfied: EPN, Tue Aug 14 07:01:44 2007 
 *
 * Generate sequences from a CM.
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

#include "funcs.h"		/* function declarations                */
#include "structs.h"		/* data structures, macros, #define's   */

#define ALPHOPTS "--rna,--dna"                         /* Exclusive options for alphabet choice */
#define OUTOPTS  "-u,-c,-a,--hmmbuild"                 /* Exclusive options for output */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs       incomp  help  docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",   1 },
  { "-l",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "local; emit from a locally configured model",  1 },
  { "-s",        eslARG_INT,    NULL,  NULL, "n>0",     NULL,      NULL,        NULL, "set random number generator seed to <n>",  1 },
  { "-n",        eslARG_INT,    "10",  NULL, "n>0",     NULL,      NULL,        NULL, "generate <n> sequences",  1 },
  { "-q",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "quiet; suppress verbose banner",         1 },
  /* output options */
  { "-o",        eslARG_OUTFILE,NULL,  NULL, NULL,      NULL,      NULL,        NULL, "save sequences in file <f>", 2 },
  { "-u",        eslARG_NONE,"default",NULL, NULL,   OUTOPTS,     NULL,        NULL, "write generated sequences as unaligned FASTA",  2 },
  { "-a",        eslARG_NONE,   FALSE, NULL, NULL,   OUTOPTS,      NULL,        NULL, "write generated sequences as a STOCKHOLM alignment",  2 },
  { "-c",        eslARG_NONE,   FALSE, NULL, NULL,   OUTOPTS,      NULL,        NULL, "generate a single \"consensus\" sequence only",  2 },
  { "--hmmbuild",eslARG_NONE,   FALSE, NULL, NULL,   OUTOPTS,      NULL,   "--tfile", "build and output a ML CM Plan 9 HMM from generated alignment", 2 },
  { "--tfile",   eslARG_OUTFILE,NULL,  NULL, NULL,      NULL,      NULL,        NULL, "dump parsetrees to file <f>",  2 },
  { "--rna",     eslARG_NONE,"default",NULL, NULL,  ALPHOPTS,      NULL,        NULL, "output alignment as RNA sequence data", 2 },
  { "--dna",     eslARG_NONE,   FALSE, NULL, NULL,  ALPHOPTS,      NULL,        NULL, "output alignment as DNA (not RNA) sequence data", 2 },
  /* Miscellaneous expert options */
  { "--zeroinserts",eslARG_NONE,FALSE, NULL, NULL,      NULL,      NULL,        NULL, "zero insert emission scores (A,C,G,U equiprobable)", 3 },
  { "--exp",     eslARG_REAL,   NULL,  NULL, "x>0",     NULL,      NULL,        NULL, "exponentiate CM probabilities by <x> before emitting",  3 },
  { "--begin",   eslARG_INT,    NULL,  NULL, "n>=1",    NULL,"-a,--end",        NULL, "truncate alignment, begin at match column <n>", 3 },
  { "--end",     eslARG_INT,    NULL,  NULL, "n>=1",    NULL,"-a,--begin",      NULL, "truncate alignment,   end at match column <n>", 3 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

/* struct cfg_s : "Global" application configuration shared by all threads/processes
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 * This strategy is used despite the fact that a MPI version of cmemit does not
 * yet exist! 
 */
struct cfg_s {
  char         *cmfile;	        /* name of input CM file  */ 
  CMFILE       *cmfp;		/* open input CM file stream       */
  ESL_ALPHABET *abc;		/* digital alphabet for CM */
  ESL_ALPHABET *abc_out; 	/* digital alphabet for writing */
  FILE         *ofp;		/* output file (default is stdout) */
  FILE         *pfp;		/* optional output file for parsetrees */
  ESL_RANDOMNESS *r;            /* source of randomness */
  int           ncm;            /* number CM we're at in file */
};

static char usage[]  = "[-options] <cmfile>";
static char banner[] = "generate sequences from a covariance model";

static int  init_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);

static void master(const ESL_GETOPTS *go, struct cfg_s *cfg);

static int initialize_cm(const ESL_GETOPTS *go, const struct cfg_s *cfg, CM_t *cm, char *errbuf);
static int emit_unaligned(const ESL_GETOPTS *go, const struct cfg_s *cfg, CM_t *cm, char *errbuf);
static int emit_alignment(const ESL_GETOPTS *go, const struct cfg_s *cfg, CM_t *cm, char *errbuf);
static int emit_consensus(const ESL_GETOPTS *go, const struct cfg_s *cfg, CM_t *cm, char *errbuf);
static int build_cp9(const ESL_GETOPTS *go, const struct cfg_s *cfg, CM_t *cm, char *errbuf);
static int truncate_msa(const ESL_GETOPTS *go, const struct cfg_s *cfg, ESL_MSA *msa, char *errbuf);


int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go = NULL;   /* command line processing                     */
  /*ESL_STOPWATCH   *w  = esl_stopwatch_Create();*/
  struct cfg_s     cfg;

  /* setup logsum lookups (could do this only if nec based on options, but this is safer) */
  init_ilogsum();
  FLogsumInit();

  /*********************************************** 
   * Parse command line
   ***********************************************/

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
      puts("\nwhere general options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1=docgroup, 2 = indentation; 80=textwidth*/
      puts("\noutput options are:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\nmiscellaneous expert options:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); 
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
  /* Initialize what we can in the config structure (without knowing the alphabet yet).
   * We could assume RNA, but this HMMER3 based approach is more general.
   */
  cfg.cmfile     = esl_opt_GetArg(go, 1); 
  cfg.cmfp       = NULL;	           /* opened in init_cfg() */
  cfg.abc        = NULL;	           /* created in init_cfg() */
  cfg.abc_out    = NULL;	           /* created in init_cfg() */
  cfg.ofp        = NULL;	           /* opened in init_cfg() */
  cfg.pfp        = NULL;	           /* opened in init_cfg() */
  cfg.r          = NULL;	           /* created in init_cfg() */

  /* do work */
  /*esl_stopwatch_Start(w);*/
  master(go, &cfg);
  /*esl_stopwatch_Stop(w);
    esl_stopwatch_Display(cfg.ofp, w, "# CPU time: ");*/

  /* Clean up the cfg. 
   */
  if (! esl_opt_IsDefault(go, "-o")) { fclose(cfg.ofp); }
  if (cfg.pfp   != NULL) fclose(cfg.pfp);
  if (cfg.abc   != NULL) { esl_alphabet_Destroy(cfg.abc); cfg.abc = NULL; }
  if (cfg.abc_out != NULL) esl_alphabet_Destroy(cfg.abc_out);
  if (cfg.cmfp  != NULL) CMFileClose(cfg.cmfp);
  if (cfg.r     != NULL) esl_randomness_Destroy(cfg.r);

  esl_getopts_Destroy(go);
  /*esl_stopwatch_Destroy(w);*/
  return 0;
}

/* init_cfg()
 * Already set:
 *    cfg->cmfile  - command line arg 1
 * Sets: 
 *    cfg->cmfp    - open CM file
 *    cfg->abc_out - digital alphabet for output
 *    cfg->ofp     - open output alignment file                
 *    cfg->pfp     - optional output file for parsetrees
 *    cfg->r       - source of randomness
 */
static int
init_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  /* open CM file for reading */
  if ((cfg->cmfp = CMFileOpen(cfg->cmfile, NULL)) == NULL)
    ESL_FAIL(eslFAIL, errbuf, "Failed to open covariance model save file %s\n", cfg->cmfile);

  /* create output alphabet */
  if      (esl_opt_GetBoolean(go, "--rna"))     cfg->abc_out = esl_alphabet_Create(eslRNA);
  else if (esl_opt_GetBoolean(go, "--dna"))     cfg->abc_out = esl_alphabet_Create(eslDNA);

  /* open output file if necessary */
  if (esl_opt_GetString(go, "-o") != NULL) {
    if ((cfg->ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) 
      ESL_FAIL(eslFAIL, errbuf, "Failed to open -o output file %s\n", esl_opt_GetString(go, "-o"));
  } else cfg->ofp = stdout;

  /* open parsetree output file if necessary */
  if (esl_opt_GetString(go, "--tfile") != NULL) {
    if ((cfg->pfp = fopen(esl_opt_GetString(go, "--tfile"), "w")) == NULL)
      ESL_FAIL(eslFAIL, errbuf, "Failed to open --tfile output file %s\n", esl_opt_GetString(go, "--tfile"));
  }

  /* create RNG */
  if (! esl_opt_IsDefault(go, "-s")) 
    cfg->r = esl_randomness_Create((long) esl_opt_GetInteger(go, "-s"));
  else cfg->r = esl_randomness_CreateTimeseeded();

  if (cfg->abc_out == NULL) ESL_FAIL(eslEINVAL, errbuf, "Output alphabet creation failed.");
  if (cfg->r       == NULL) ESL_FAIL(eslEINVAL, errbuf, "Failed to create random number generator: probably out of memory");

  return eslOK;
}

/* master()
 * The serial version of cmemit. (There is no parallel version yet).
 * For each CM, emit sequences/alignment and create output.
 * 
 * We only return if successful. All errors are handled immediately and fatally with cm_Fail().
 */
static void
master(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int      status;
  char     errbuf[eslERRBUFSIZE];
  CM_t    *cm = NULL;

  if ((status = init_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);

  cfg->ncm = 0;

  while (CMFileRead(cfg->cmfp, &(cfg->abc), &cm))
  {
    if (cm == NULL) cm_Fail("Failed to read CM from %s -- file corrupt?\n", cfg->cmfile);
    cfg->ncm++;
    if((status = initialize_cm(go, cfg, cm, errbuf)) != eslOK) cm_Fail(errbuf);
    /* Pick 1 of 4 exclusive output options. Output is handled within each function. */
    if     (esl_opt_GetBoolean(go, "-u")) { 
      if((status = emit_unaligned(go, cfg, cm, errbuf)) != eslOK) cm_Fail(errbuf);
    }
    else if(esl_opt_GetBoolean(go, "-c")) {
      if((status = emit_consensus(go, cfg, cm, errbuf)) != eslOK) cm_Fail(errbuf);
    }
    else if(esl_opt_GetBoolean(go, "-a")) {
      if((status = emit_alignment(go, cfg, cm, errbuf)) != eslOK) cm_Fail(errbuf);
    }
    else if(esl_opt_GetBoolean(go, "--hmmbuild")) {
      if((status = build_cp9     (go, cfg, cm, errbuf)) != eslOK) cm_Fail(errbuf);
    }
    FreeCM(cm);
  }
  return;
}

/* initialize_cm()
 * Setup the CM based on the command-line options/defaults;
 * only set flags and a few parameters. ConfigCM() configures
 * the CM.
 */
static int
initialize_cm(const ESL_GETOPTS *go, const struct cfg_s *cfg, CM_t *cm, char *errbuf)
{
  /* Update cfg->cm->config_opts and cfg->cm->align_opts based on command line options */
  if(esl_opt_GetBoolean(go, "-l"))            cm->config_opts |= CM_CONFIG_LOCAL;
  if(esl_opt_GetBoolean(go, "--zeroinserts")) cm->config_opts |= CM_CONFIG_ZEROINSERTS;
  ConfigCM(cm, NULL, NULL);
  if(! esl_opt_IsDefault(go, "--exp"))        ExponentiateCM(cm, esl_opt_GetReal(go, "--exp"));
  
  return eslOK;
}

/* emit_unaligned
 * Given a configured CM, generate and output unaligned sequences.
 */
static int
emit_unaligned(const ESL_GETOPTS *go, const struct cfg_s *cfg, CM_t *cm, char *errbuf)
{
  /* Contract check, output alphabet must be identical to CM alphabet 
   * with sole exception that CM alphabet can be eslRNA with output alphabet eslDNA.
   */
  if(cm->abc->type != cfg->abc_out->type)
    if(! (cm->abc->type == eslRNA && cfg->abc_out->type != eslDNA))
      ESL_FAIL(eslFAIL, errbuf, "CM alphabet type must match output alphabet type (unless CM=RNA and output=DNA).");

  int status;
  Parsetree_t *tr = NULL;  /* generated parsetree */
  ESL_SQ *sq = NULL;
  char *name;
  int namelen;
  int i, L; 

  namelen = IntMaxDigits() + 1;  /* IntMaxDigits() returns number of digits in INT_MAX */
  if(cm->name != NULL) namelen += strlen(cm->name) + 1;
  ESL_ALLOC(name, sizeof(char) * namelen);

  for(i = 0; i < esl_opt_GetInteger(go, "-n"); i++)
    {
      if(cm->name != NULL) sprintf(name, "%s-%d", cm->name, i+1);
      else                 sprintf(name, "%d-%d", cfg->ncm, i+1);
      EmitParsetree(cm, cfg->r, name, TRUE, &tr, &sq, &L);
      sq->abc = cfg->abc_out;
      if((esl_sqio_Write(cfg->ofp, sq, eslSQFILE_FASTA)) != eslOK) 
	ESL_FAIL(eslFAIL, errbuf, "Error writing unaligned sequences.");
      if(cfg->pfp != NULL)
	{
	  fprintf(cfg->pfp, "> %s\n", sq->name);
	  ParsetreeDump(cfg->pfp, tr, cm, sq->dsq, NULL, NULL);
	  fprintf(cfg->pfp, "//\n");
	}
      FreeParsetree(tr);
      esl_sq_Destroy(sq); /* can't reuse b/c a new one is allocated in EmitParsetree() */
    }
  if(cfg->pfp != NULL) printf("Parsetrees saved to file %s.\n", esl_opt_GetString(go, "--tfile"));
  free(name);
  return eslOK;

 ERROR:
  return status;
}

/* emit_alignment
 * Given a configured CM, generate and output a MSA.
 */
static int
emit_alignment(const ESL_GETOPTS *go, const struct cfg_s *cfg, CM_t *cm, char *errbuf)
{
  /* Contract check, output alphabet must be identical to CM alphabet 
   * with sole exception that CM alphabet can be eslRNA with output alphabet eslDNA. */
  if(cm->abc->type != cfg->abc_out->type)
    if(! (cm->abc->type == eslRNA && cfg->abc_out->type != eslDNA))
      ESL_FAIL(eslFAIL, errbuf, "CM alphabet type must match output alphabet type (unless CM=RNA and output=DNA).");

  int status;
  Parsetree_t **tr = NULL;  /* generated parsetrees */
  ESL_SQ **sq = NULL;       /* generated sequences */
  char *name;
  int namelen;
  int i, L; 
  ESL_MSA *msa = NULL;      /* the MSA we're building */
  int nseq = esl_opt_GetInteger(go, "-n");
  int do_truncate;

  if((! esl_opt_IsDefault(go, "--begin")) && (! esl_opt_IsDefault(go, "--end"))) 
    do_truncate = TRUE;
  else do_truncate = FALSE;

  namelen = IntMaxDigits() + 1;
  if(cm->name != NULL) namelen += strlen(cm->name) + 1;
  ESL_ALLOC(name, sizeof(char) * namelen);

  ESL_ALLOC(sq, sizeof(ESL_SQ *) * nseq);
  ESL_ALLOC(tr, sizeof(Parsetree_t *) * nseq);

  for(i = 0; i < nseq; i++)
    {
      if(cm->name != NULL) sprintf(name, "%s-%d", cm->name, i+1);
      else                 sprintf(name, "%d-%d", cfg->ncm, i+1);
      EmitParsetree(cm, cfg->r, name, TRUE, &(tr[i]), &(sq[i]), &L);
      sq[i]->abc = cfg->abc_out;
      if(cfg->pfp != NULL)
	{
	  fprintf(cfg->pfp, "> %s\n", sq[i]->name);
	  ParsetreeDump(cfg->pfp, tr[i], cm, sq[i]->dsq, NULL, NULL);
	  fprintf(cfg->pfp, "//\n");
	}
    }
  if((status = Parsetrees2Alignment(cm, cfg->abc_out, sq, NULL, tr, nseq,
				    TRUE,  /* we want all match columns */
				    FALSE, /* we don't want ONLY match columns */
				    &msa) != eslOK))
    ESL_XFAIL(eslFAIL, errbuf, "Error generating alignment from parsetrees.");
  if(cm->name != NULL) 
    if((status = esl_strdup(cm->name, -1, &(msa->name))) != eslOK) goto ERROR;
  if((status = esl_strdup("Synthetic sequence alignment generated by cmemit", 
			  -1, &(msa->desc))) != eslOK)  goto ERROR;

  /* Truncate the alignment if nec */
  if(do_truncate)
    if((status = truncate_msa(go, cfg, msa, errbuf)) != eslOK) 
      ESL_XFAIL(status, errbuf, "Error truncating alignment.");

  /* Output the alignment */
  status = esl_msa_Write(cfg->ofp, msa, eslMSAFILE_STOCKHOLM);
  if      (status == eslEMEM) ESL_XFAIL(status, errbuf, "Memory error when outputting alignment\n");
  else if (status != eslOK)   ESL_XFAIL(status, errbuf, "Writing alignment file failed with error %d\n", status);

  if(cfg->pfp != NULL) printf("Parsetrees saved to file %s.\n", esl_opt_GetString(go, "--tfile"));
  for(i = 0; i < nseq; i++)
    {
      esl_sq_Destroy(sq[i]);
      FreeParsetree(tr[i]);
    }
  free(sq);
  free(tr);
  free(name);
  esl_msa_Destroy(msa);
  return eslOK;

 ERROR:
  
  if(sq != NULL) {
    for(i = 0; i < nseq; i++)
      if(sq[i] != NULL) esl_sq_Destroy(sq[i]);
    free(sq);
  }
  if(tr != NULL) {
    for(i = 0; i < nseq; i++)
      if(tr[i] != NULL) FreeParsetree(tr[i]);
    free(tr);
  }
  if(name != NULL) free(name);
  if(msa  != NULL) esl_msa_Destroy(msa);
  return status;
}

/* emit_consensus
 * Given a configured CM, generate and output the consensus sequence.
 */
static int
emit_consensus(const ESL_GETOPTS *go, const struct cfg_s *cfg, CM_t *cm, char *errbuf)
{
  int status;
  CMConsensus_t *con = NULL;            /* consensus info for the CM */
  ESL_SQ *csq = NULL;

  /* Determine consensus sequence */
  CreateCMConsensus(cm, cfg->abc_out, 3.0, 1.0, &con);
  if((csq = esl_sq_CreateFrom("CM generated consensus sequence [cmemit]", con->cseq, NULL, NULL, NULL)) == NULL)
    { status = eslEMEM; goto ERROR; }
  if((esl_sqio_Write(cfg->ofp, csq, eslSQFILE_FASTA)) != eslOK) 
    ESL_FAIL(status, errbuf, "Error writing unaligned sequences.");

  esl_sq_Destroy(csq);
  FreeCMConsensus(con);
  return eslOK;

 ERROR:
  if(csq != NULL) esl_sq_Destroy(csq);
  if(con != NULL) FreeCMConsensus(con);
  return status;
}

/* build_cp9
 * Given a configured CM, generete counts for a ML HMM
 * (no pseudocounts) by generating >= 1 MSA from the CM.
 * We use more than 1 MSA only to limit memory usage.
 */
static int
build_cp9(const ESL_GETOPTS *go, const struct cfg_s *cfg, CM_t *cm, char *errbuf)
{
  int status;
  Parsetree_t **tr = NULL;  /* generated parsetrees */
  ESL_SQ **sq = NULL;       /* generated sequences */
  char *name;
  int namelen;
  int i, L; 
  ESL_MSA *msa = NULL;      /* an MSA to pull counts from */
  int nseq = esl_opt_GetInteger(go, "-n");
  int *matassign = NULL;
  int nsampled = 0;                 /* number of sequences sampled thus far */
  int do_truncate;
  CP9_t  *shmm = NULL;
  CP9trace_t **cp9_tr;   /* fake tracebacks for each seq            */
  int bpos = 0;
  int epos = 0;
  int apos = 0;
  int msa_nseq = 1000;          /* number of seqs per MSA, current strategy is to 
				 * sample (nseq/nseq_per_msa) alignments from the CM, 
				 * and add counts from each to the shmm in counts form 
				 *(to limit memory) */
  
  /* Allocate and zero the new HMM we're going to build by sampling from the CM.
   */
  if((! esl_opt_IsDefault(go, "--begin")) && (! esl_opt_IsDefault(go, "--end"))) 
    {
      do_truncate = TRUE;
      bpos = esl_opt_GetInteger(go, "--begin");
      epos = esl_opt_GetInteger(go, "--end");
      shmm = AllocCPlan9((epos - bpos + 1), cm->abc);
    }
  else 
    {
      do_truncate = FALSE;
      shmm = AllocCPlan9(cm->clen, cm->abc);
    }
  ZeroCPlan9(shmm);

  namelen = IntMaxDigits() + 1;
  if(cm->name != NULL) namelen += strlen(cm->name) + 1;
  ESL_ALLOC(name, sizeof(char) * namelen);

  /* sample MSA(s) from the CM */
  ESL_ALLOC(sq,     sizeof(ESL_SQ *)      * msa_nseq);
  ESL_ALLOC(tr,     sizeof(Parsetree_t *) * msa_nseq);
  ESL_ALLOC(cp9_tr, sizeof(CP9trace_t *)  * msa_nseq);
  while(nsampled < nseq)
    {
      if(nsampled != 0) 
	{
	  /* clean up from previous MSA */
	  esl_msa_Destroy(msa);
	  free(matassign);
	  for (i = 0; i < msa_nseq; i++)
	    {
	      CP9FreeTrace(cp9_tr[i]);
	      FreeParsetree(tr[i]);
	      esl_sq_Reuse(sq[i]);
	    }
	}
      /* Emit msa_nseq parsetrees from the CM */
      if(nsampled + msa_nseq > nseq) msa_nseq = nseq - nsampled;
      for (i = 0; i < msa_nseq; i++)
	{
	  if(cm->name != NULL) sprintf(name, "%s-%d", cm->name, i+1);
	  else                 sprintf(name, "%d-%d", cfg->ncm, i+1);
	  EmitParsetree(cm, cfg->r, name, TRUE, &(tr[i]), &(sq[i]), &L);
	  sq[i]->abc = cfg->abc_out;
	}
      /* Build a new MSA from these parsetrees */
      if((status = Parsetrees2Alignment(cm, cfg->abc_out, sq, NULL, tr, nseq,
					TRUE,  /* we want all match columns */
					FALSE, /* we don't want ONLY match columns */
					&msa) != eslOK))
	ESL_XFAIL(eslFAIL, errbuf, "Error generating alignment from parsetrees during HMM construction.");

      /* Truncate the alignment if nec */
      if(do_truncate)
	if((status = truncate_msa(go, cfg, msa, errbuf)) != eslOK) 
	  ESL_XFAIL(status, errbuf, "Error truncating alignment during HMM construction.");

      /* Determine match assignment from RF annotation
       */
      ESL_ALLOC(matassign, sizeof(int) * (msa->alen+1));
      matassign[0] = 0;
      for (apos = 1; apos <= msa->alen; apos++)
	matassign[apos] = esl_abc_CIsGap(msa->abc, msa->rf[apos-1]) ? FALSE : TRUE;  

      /* Add the counts to the growing counts-based HMM */
      /* make fake tracebacks for each seq, first we need to digitize the MSA */
      esl_msa_Digitize(msa->abc, msa);
      CP9_fake_tracebacks(msa, matassign, &cp9_tr);
	  
      /* build model from tracebacks (code from HMMER's modelmakers.c::matassign2hmm() */
      for (i = 0; i < msa->nseq; i++) 
	CP9TraceCount(shmm, sq[i]->dsq, 1.0, cp9_tr[i]);
      nsampled += msa_nseq;
    }
      
  /* clean up from final MSA */
  esl_msa_Destroy(msa);
  free(matassign);
  for (i = 0; i < msa_nseq; i++)
    {
      CP9FreeTrace(cp9_tr[i]);
      FreeParsetree(tr[i]);
      esl_sq_Destroy(sq[i]);
    }
  free(cp9_tr);
  free(tr);
  free(sq);
  free(name);

  printf("Printing NON-normalized sampled HMM parameters:\n");
  debug_print_cp9_params(stdout, shmm, FALSE);
  printf("finished printing NON-normalized sampled HMM parameters.\n");

  CPlan9Renormalize(shmm);
  CP9Logoddsify(shmm);

  printf("Printing normalized sampled HMM parameters:\n");
  debug_print_cp9_params(stdout, shmm, TRUE);
  printf("finished printing normalized sampled HMM parameters.\n");

  FreeCPlan9(shmm);
  return eslOK;

 ERROR:
  if(shmm != NULL) FreeCPlan9(shmm);
  if(cp9_tr != NULL)
    {
      for(i = 0; i < msa_nseq; i++)
	CP9FreeTrace(cp9_tr[i]);
      free(cp9_tr);
    }
  if(tr != NULL)
    {
      for(i = 0; i < msa_nseq; i++)
	FreeParsetree(tr[i]);
      free(tr);
    }
  if(sq != NULL)
    {
      for(i = 0; i < msa_nseq; i++)
	esl_sq_Destroy(sq[i]);
      free(sq);
    }
  if(name != NULL) free(name);
  return status;
}

/* truncate_msa
 * Truncate a MSA outside begin..end consensus columns 
 * (non-gap RF chars) and return the alignment. Careful
 * to remove any consensus structure outside begin..end.
 */
static int
truncate_msa(const ESL_GETOPTS *go, const struct cfg_s *cfg, ESL_MSA *msa, char *errbuf)
{
  int status;
  int *useme = NULL;    /* 1..alen: keep this column? */
  int *ct    = NULL;    /* 1..alen base pair partners array */
  int apos = 0;
  int bpos = esl_opt_GetInteger(go, "--begin");
  int epos = esl_opt_GetInteger(go, "--end");
  int cc = 0;

  ESL_ALLOC(useme, sizeof(int) * (msa->alen+1));
  ESL_ALLOC(ct,    sizeof(int) * (msa->alen+1));

  int clen = 0;
  for (apos = 0, cc = 0; apos < msa->alen; apos++)
    if (!esl_abc_CIsGap(msa->abc, msa->rf[apos])) clen++;
  if(esl_opt_GetInteger(go, "--end") > clen)
    ESL_XFAIL(status, errbuf, "Error, with --end <n> option, <n> must be <= consensus length of CM (%d).\n", clen);

  /* remove pknots in place (actually unnec for CM ss_cons) */
  if((status = esl_wuss_nopseudo(msa->ss_cons, msa->ss_cons)) != eslOK) goto ERROR; 

  /* get a ct array from the structure */
  if((status = esl_wuss2ct(msa->ss_cons, msa->alen, ct)) != eslOK) goto ERROR;  

  /* Truncate the alignment prior to consensus column bpos and after 
   * consensus column epos.  */
  for (apos = 0, cc = 0; apos < msa->alen; apos++)
    {
      /* Careful here, placement of cc++ increment is impt, we want all 
       * inserts between cc=bpos-1 and cc=bpos, and b/t cc=epos and 
       * cc=epos+1. Also be careful: ct[] is index 1..alen, and 
       * msa->ss_cons is 0..alen-1. 
       */
      if(cc < (bpos-1) || cc > epos) {
	useme[apos] = 0;
	if(ct[(apos+1)] != 0) ct[ct[(apos+1)]] = 0;
	ct[(apos+1)] = 0;
      }
      else
	useme[apos] = 1;
      if (!esl_abc_CIsGap(msa->abc, msa->rf[apos])) { 
	cc++; 
	if(cc == (epos+1)){
	  useme[apos] = 0; 
	  /* we misassigned this guy, overwrite */ 
	  if(ct[(apos+1)] != 0) ct[ct[(apos+1)]] = 0;
	  ct[(apos+1)] = 0;
	}
      }
    }
  /* construct the new structure based on the cleaned ct array */
  if((status = esl_ct2wuss(ct, msa->alen, msa->ss_cons)) != eslOK) goto ERROR;
  
  /*printf("\n\nDEBUG PRINTING ORIG ALIGNMENT:\n");
    WriteStockholm(fp, msa);
    printf("\n\nDONE DEBUG PRINTING ORIG ALIGNMENT:\n");
    for(apos=0; apos < msa->alen; apos++)
    printf("useme[%d]: %d\n", apos, useme[apos]);
  */

  esl_msa_ColumnSubset(msa, useme);
  free(useme);
  free(ct);
  return eslOK;

 ERROR:
  if(useme != NULL) free(useme);
  if(ct    != NULL) free(ct);
  return status;
}






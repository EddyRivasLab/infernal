/* cmbuild.c
 * SRE, Thu Jul 27 13:19:43 2000 [StL]
 * SVN $Id$
 * 
 * Construct a CM from a given multiple sequence alignment.
 *  
 *****************************************************************
 * @LICENSE@
 *****************************************************************
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msaweight.h"
#include "esl_msacluster.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#include "structs.h"		/* data structures, macros, #define's   */
#include "prior.h"		/* mixture Dirichlet prior */
#include "funcs.h"		/* external functions                   */
#include "hmmband.h"
#include "stats.h"              /* for resolve_degenerate               */

#define WGTOPTS "--wgsc,--wblosum,--wpb,--wnone,--wgiven"      /* Exclusive options for relative weighting                    */
#define EFFOPTS "--eent,--enone"               /* Exclusive options for effective sequence number calculation */
#define ALPHOPTS "--rna"                       /* Exclusive options for alphabet choice */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs       incomp  help  docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",   1 },
  { "-n",        eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL,        NULL, "name the CM(s) <s>",                     1 },
  { "-o",        eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "direct summary output to file <f>, not stdout", 1 },
  { "-A",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "append this CM to <cmfile>",             1 },
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "force; allow overwriting of <cmfile>",   1 },
  { "-1",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "use tabular output summary format, 1 line per CM", 1 },
  { "--binary",  eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "save the model(s) in binary format",     2 },
  { "--rf",      eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,  "--rsearch", "use reference coordinate annotation to specify consensus", 2 },
  { "--informat",eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL,        NULL, "specify input alignment is in format <s>, not Stockholm",  2 },
  { "--gapthresh",eslARG_REAL,  "0.5", NULL, "0<=x<=1", NULL,      NULL,  "--rsearch", "fraction of gaps to allow in a consensus column [0..1]", 2 },
  { "--rsearch", eslARG_STRING, NULL,  NULL, NULL,      NULL, "--enone",        NULL,  "use RSEARCH parameterization with RIBOSUM matrix file <s>", 2 }, 
  { "--elself",  eslARG_REAL,  "0.94", NULL, "0<=x<=1", NULL,      NULL,        NULL, "set EL self transition prob to <x> [df: 0.94]", 2 },
  { "--nodetach",eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "do not 'detach' one of two inserts that model same column", 2 },
/* Alternate relative sequence weighting strategies */
  /* --wme not implemented in Infernal yet (b/c it's not in HMMER3) */
  { "--wgsc",    eslARG_NONE,"default",NULL, NULL,    WGTOPTS,    NULL,      NULL, "Gerstein/Sonnhammer/Chothia tree weights",         3},
  { "--wblosum", eslARG_NONE,  FALSE,  NULL, NULL,    WGTOPTS,    NULL,      NULL, "Henikoff simple filter weights",                   3},
  { "--wpb",     eslARG_NONE,  FALSE,  NULL, NULL,    WGTOPTS,    NULL,      NULL, "Henikoff position-based weights",                  3},
  { "--wnone",   eslARG_NONE,  FALSE,  NULL, NULL,    WGTOPTS,    NULL,      NULL, "don't do any relative weighting; set all to 1",    3},
  { "--wgiven",  eslARG_NONE,  FALSE,  NULL, NULL,    WGTOPTS,    NULL,      NULL, "use weights as given in MSA file",                 3},
  { "--pbswitch",eslARG_INT,  "1000",  NULL,"n>0",       NULL,    NULL,      NULL, "set failover to efficient PB wgts at > <n> seqs",  3},
  { "--wid",     eslARG_REAL, "0.62",  NULL,"0<=x<=1",   NULL,"--wblosum",   NULL, "for --wblosum: set identity cutoff",               3},
/* Alternate effective sequence weighting strategies */
  { "--eent",    eslARG_NONE,"default",NULL, NULL,    EFFOPTS,    NULL,      NULL, "adjust eff seq # to achieve relative entropy target", 4},
  { "--enone",   eslARG_NONE,  FALSE,  NULL, NULL,    EFFOPTS,    NULL,      NULL, "no effective seq # weighting: just use nseq",         4},
  { "--ere",     eslARG_REAL,  NULL,   NULL,"x>0",       NULL, "--eent",     NULL, "for --eent: set target relative entropy to <x>",      4},
  { "--eX",      eslARG_REAL,  "6.0",  NULL,"x>0",       NULL, "--eent",  "--ere", "for --eent: set minimum total rel ent param to <x>",  4}, 
/* Verbose output files */
  { "--cfile",   eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL,        NULL, "save count vectors to file <s>", 5 },
  { "--cmtbl",   eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL,        NULL, "save tabular description of CM topology to file <s>", 5 },
  { "--emap",    eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL,        NULL, "save consensus emit map to file <s>", 5 },
  { "--gtree",   eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL,        NULL, "save tree description of master tree to file <s>", 5 },
  { "--gtbl",    eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL,        NULL, "save tabular description of master tree to file <s>", 5 },
  { "--tfile",   eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL,        NULL, "dump individual sequence tracebacks to file <s>", 5 },
  { "--bfile",   eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL,        NULL, "save bands to file <f>, which can be read by cmsearch", 5 },
  { "--nobalance",eslARG_NONE,  FALSE, NULL, NULL,      NULL,      NULL,        NULL, "don't rebalance the CM; number in strict preorder", 6 },
  { "--regress",  eslARG_STRING, NULL, NULL, NULL,      NULL,      NULL,        NULL, "save regression test information to file <s>", 6 },  
  { "--treeforce",eslARG_NONE,  FALSE, NULL, NULL,      NULL,      NULL,        NULL, "score first seq in alignment and show parsetree", 6 },
  { "--ignorant", eslARG_NONE,  FALSE, NULL, NULL,      NULL,      NULL,        NULL, "strip the structural info from input alignment", 6 },
  { "--null",    eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL, "--rsearch", "read null (random sequence) model from file <s>", 7 },
  { "--prior",   eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL, "--rsearch", "read priors from file <s>", 7 },
  { "--ctarget", eslARG_INT,    "0",   NULL, "n>=0",    NULL,      NULL,    "--call", "build (at most) <n> CMs by partitioning MSA into <n> clusters", 8 },
  { "--cmindiff",eslARG_REAL,   "0.",  NULL,"0.<=x<=1.",NULL,      NULL,    "--call", "min difference b/t 2 clusters is <x>, each cluster -> CM", 8 }, 
  { "--call",    eslARG_NONE,  FALSE,  NULL, NULL,      NULL,      NULL,        NULL, "build a separate CM from every seq in MSA", 8 },
  { "--corig",   eslARG_NONE,  FALSE,  NULL, NULL,      NULL,      NULL,        NULL, "build an additional CM from the original, full MSA", 8 }, 
  { "--cdump",   eslARG_STRING, NULL,  NULL, NULL,      NULL,      NULL,        NULL, "dump an MSA for each cluster (CM) to file <s>", 8 },
/* Selecting the alphabet rather than autoguessing it */
  { "--rna",     eslARG_NONE,   FALSE, NULL, NULL,   ALPHOPTS,    NULL,     NULL, "input alignment is RNA sequence data",                  9},
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

/* struct cfg_s : "Global" application configuration shared by all threads/processes
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 * This strategy is used despite the fact that a MPI version of cmbuild does not
 * yet exist! 
 */
struct cfg_s {
  FILE         *ofp;		/* output file (default is stdout) */

  char         *alifile;	/* name of the alignment file we're building CMs from  */
  int           fmt;		/* format code for alifile */
  ESL_MSAFILE  *afp;            /* open alifile  */
  ESL_ALPHABET *abc;		/* digital alphabet */

  char         *cmfile;         /* file to write CM to                    */
  FILE         *cmfp;           /* CM output file handle                  */

  float        *null;		/* null model                              */
  Prior_t      *pri;		/* mixture Dirichlet prior for the HMM     */

  fullmat_t    *fullmat;        /* if --rsearch, the full RIBOSUM matrix */
  FILE         *cdfp;           /* if --cdump, output file handle for dumping clustered MSAs */

  int           be_verbose;	/* standard verbose output, as opposed to one-line-per-CM summary */
  int           nali;		/* which # alignment this is in file (only valid in serial mode)   */
};

static char usage[]  = "[-options] <cmfile output> <alignment file>";
static char banner[] = "build RNA covariance model(s) from alignment";

static int  init_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);

static void  master (const ESL_GETOPTS *go, struct cfg_s *cfg);

static int    process_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t **ret_cm, Parsetree_t ***opt_tr);
static int    output_result(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, int msaidx, ESL_MSA *msa, CM_t *cm);
static int    check_and_clean_msa(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa);
static int    set_relative_weights(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa);
static int    build_model(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t **ret_cm, Parsetree_t ***opt_tr);
static int    set_model_name(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t *cm);
static int    set_effective_seqnumber(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t *cm, const Prior_t *pri);
static int    parameterize(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, const Prior_t *prior);
static int    name_msa(const ESL_GETOPTS *go, ESL_MSA *msa, int nali);
static double default_target_relent(const ESL_ALPHABET *abc, int M, double eX);
static int    save_countvectors(char *cfile, CM_t *cm);
static void   model_trace_info_dump(FILE *ofp, CM_t *cm, Parsetree_t *tr, char *aseq);
static void   strip_wuss(char *ss);


int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go = NULL;   /* command line processing                     */
  ESL_STOPWATCH   *w  = esl_stopwatch_Create();
  struct cfg_s     cfg;

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
      puts("\nexpert model construction options:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\nsequence weighting options [default: GSC weighting]:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); 
      puts("\neffective sequence number related options:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80);
      puts("\nverbose output files, useful for detailed information about the CM:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80);
      puts("\ndebugging, experimentation:");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80);
      puts("\ncustomization of null model and priors:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80);
      puts("\noptions for building multiple CMs after clustering input MSA:");
      esl_opt_DisplayHelp(stdout, go, 8, 2, 80);
      puts("\n  options for selecting alphabet rather than guessing it:");
      esl_opt_DisplayHelp(stdout, go, 9, 2, 80);
      exit(0);
    }
  if (esl_opt_ArgNumber(go) != 2) 
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
  cfg.ofp        = NULL;	           /* opened in init_cfg() */
  cfg.alifile    = esl_opt_GetArg(go, 2);
  cfg.fmt        = eslMSAFILE_UNKNOWN;     /* autodetect alignment format by default. */ 
  cfg.afp        = NULL;	           /* created in init_cfg() */
  cfg.abc        = NULL;	           /* created in init_cfg() */
  cfg.cmfile     = esl_opt_GetArg(go, 1); 
  cfg.cmfp       = NULL;	           /* opened in init_cfg() */
  cfg.null       = NULL;	           /* created in init_cfg() */
  cfg.pri        = NULL;                   /* created in init_cfg() */
  cfg.fullmat    = NULL;                   /* read (possibly) in init_cfg() */
  cfg.cdfp       = NULL;	           /* opened (possibly) in init_cfg() */

  if (esl_opt_GetBoolean(go, "-1")) cfg.be_verbose = FALSE;        
  else                              cfg.be_verbose = TRUE;        
  cfg.nali       = 0;		           

  /* Start timing; do work; stop timing.*/
  esl_stopwatch_Start(w);
  master(go, &cfg);
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(cfg.ofp, w, "# CPU time: ");

  /* Clean up the cfg. 
   */
  if (! esl_opt_IsDefault(go, "-o")) { fclose(cfg.ofp); }
  if (cfg.ofp   != NULL) esl_msafile_Close(cfg.afp);
  if (cfg.abc   != NULL) esl_alphabet_Destroy(cfg.abc);
  if (cfg.cmfp  != NULL) fclose(cfg.cmfp);
  if (cfg.pri   != NULL) Prior_Destroy(cfg.pri);
  if (cfg.null  != NULL) free(cfg.null);
  if (cfg.cdfp  != NULL) fclose(cfg.cdfp);

  esl_getopts_Destroy(go);
  esl_stopwatch_Destroy(w);
  return 0;
}

/* init_cfg()
 * Already set:
 *    cfg->cmfile  - command line arg 1
 *    cfg->alifile - command line arg 2
 *    cfg->fmt     - format of alignment file
 * Sets: 
 *    cfg->afp     - open alignment file                
 *    cfg->abc     - digital alphabet
 *    cfg->cmfp    - open CM file
 *    cfg->null    - NULL model, used for all models
 *    cfg->pri     - prior, used for all models
 *    cfg->fullmat - RIBOSUM matrix used for all models (optional)
 *    cfg->cdfp    - open file to dump MSAs to (optional)
 */
static int
init_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  int status;

  if (esl_opt_GetString(go, "-o") != NULL) {
    if ((cfg->ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) 
      ESL_FAIL(eslFAIL, errbuf, "Failed to open -o output file %s\n", esl_opt_GetString(go, "-o"));
  } else cfg->ofp = stdout;

  status = esl_msafile_Open(cfg->alifile, cfg->fmt, NULL, &(cfg->afp));
  if (status == eslENOTFOUND)    ESL_FAIL(status, errbuf, "Alignment file %s doesn't exist or is not readable\n", cfg->alifile);
  else if (status == eslEFORMAT) ESL_FAIL(status, errbuf, "Couldn't determine format of alignment %s\n", cfg->alifile);
  else if (status != eslOK)      ESL_FAIL(status, errbuf, "Alignment file open failed with error %d\n", status);
  cfg->fmt = cfg->afp->format;

  /* Guess alphabet, then make sure it's RNA or DNA */
  int type;
  status = esl_msafile_GuessAlphabet(cfg->afp, &type);
  if (status == eslEAMBIGUOUS)    ESL_FAIL(status, errbuf, "Failed to guess the bio alphabet used in %s.\nUse --rna option to specify it.", cfg->alifile);
  else if (status == eslEFORMAT)  ESL_FAIL(status, errbuf, "Alignment file parse failed: %s\n", cfg->afp->errbuf);
  else if (status == eslENODATA)  ESL_FAIL(status, errbuf, "Alignment file %s is empty\n", cfg->alifile);
  else if (status != eslOK)       ESL_FAIL(status, errbuf, "Failed to read alignment file %s\n", cfg->alifile);
  /* We can read DNA/RNA but internally we treat it as RNA */
  if(! (type == eslRNA || type == eslDNA))
    ESL_FAIL(status, errbuf, "Alphabet is not DNA/RNA in %s\n", cfg->alifile);
  cfg->abc = esl_alphabet_Create(eslRNA);
  esl_msafile_SetDigital(cfg->afp, cfg->abc);

  /* open CM file for writing */
  if ((cfg->cmfp = fopen(cfg->cmfile, "w")) == NULL) ESL_FAIL(status, errbuf, "Failed to open CM file %s for writing", cfg->cmfile);

  /* Set up the prior */
  if (esl_opt_GetString(go, "--prior") != NULL)
    {
      FILE *pfp;
      if ((pfp = fopen(esl_opt_GetString(go, "--prior"), "r")) == NULL)
	esl_fatal("Failed to open prior file %s\n", esl_opt_GetString(go, "--prior"));
      if ((cfg->pri = Prior_Read(pfp)) == NULL)
	esl_fatal("Failed to parse prior file %s\n", esl_opt_GetString(go, "--prior"));
      fclose(pfp);
    }
  else 
    cfg->pri = Prior_Default();

  /* Set up the null/random seq model */
  if(esl_opt_GetString(go, "--null") != NULL) /* read freqs from a file and overwrite bg->f */
    {
      if((status = CMReadNullModel(cfg->abc, esl_opt_GetString(go, "--null"), &(cfg->null))) != eslOK)
	cm_Fail("Failure reading the null model, code: %d", status);
    }       
  else /* set up the default null model */
    {
      status = DefaultNullModel(cfg->abc, &(cfg->null)); /* default values, A,C,G,U = 0.25  */
      if(status != eslOK) cm_Fail("Failure creating the null model, code: %d", status);
    }

  /* if --rsearch was enabled, set up RIBOSUM matrix */
  if(esl_opt_GetString(go, "--rsearch") != NULL)
    {
      FILE *matfp;
      if ((matfp = MatFileOpen (esl_opt_GetString(go, "--rsearch"))) == NULL)
	cm_Fail("Failed to open matrix file %s\n", esl_opt_GetString(go, "--rsearch"));
      if (! (cfg->fullmat = ReadMatrix(cfg->abc, matfp)))
	cm_Fail("Failed to read matrix file %s\n", esl_opt_GetString(go, "--rsearch"));
      ribosum_calc_targets(cfg->fullmat); /* overwrite score matrix scores w/target probs */
      fclose(matfp);
    }

  /* if --cdump enabled, open output file for MSAs */
  if (esl_opt_GetString(go, "--cdump") != NULL)
    {
      if ((cfg->cdfp = fopen(esl_opt_GetString(go, "--cdump"), "w")) == NULL)
	cm_Fail("Failed to open output file %s for writing MSAs to", esl_opt_GetString(go, "--cdump"));
    }

  if (cfg->pri   == NULL) ESL_FAIL(eslEINVAL, errbuf, "alphabet initialization failed");
  if (cfg->null  == NULL) ESL_FAIL(eslEINVAL, errbuf, "null model initialization failed");

  return eslOK;
}

/* master()
 * The serial version of cmbuild. (There is no parallel version yet).
 * For each MSA, build at least one CM and save it.
 * 
 * We only return if successful. All errors are handled immediately and fatally with cm_Fail().
 */
static void
master(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int      status;
  char     errbuf[eslERRBUFSIZE];
  ESL_MSA *msa = NULL;
  CM_t    *cm = NULL;

  /* cluster option related variables */
  int          do_cluster; /* TRUE if --ctarget || --cmindiff || --call */
  int          ncm = 1;    /* number of CMs to be built for current MSA */
  int          c   = 0;    /* counter over CMs built for a single MSA */
  ESL_MSA    **cmsa;       /* pointer to cluster MSAs to build CMs from */
  
  if ((status = init_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);

  cfg->nali = 0;

  do_cluster = FALSE;
  if((esl_opt_GetInteger(go, "--ctarget"))  || (esl_opt_GetReal   (go, "--cmindiff")) || 
     (esl_opt_GetBoolean(go, "--call")))
    do_cluster = TRUE;

  while ((status = esl_msa_Read(cfg->afp, &msa)) != eslEOF)
    {
      if      (status == eslEFORMAT)  cm_Fail("Alignment file parse error, line %d of file %s:\n%s\nOffending line is:\n%s\n", cfg->afp->linenumber, cfg->afp->fname, cfg->afp->errbuf, cfg->afp->buf);
      else if (status != eslOK)       cm_Fail("Alignment file read unexpectedly failed with code %d\n", status);
      cfg->nali++;  

      /* if it's unnamed, name the MSA, we require a name (different from 
       * HMMER 3), because it will be used to name the CM. */
      if(name_msa(go, msa, cfg->nali) != eslOK) cm_Fail("Error (code: %d) naming MSA", status);
      if(msa->name == NULL)                     cm_Fail("Error naming MSA");
      ncm = 1;     /* default: only build 1 CM for each MSA in alignment file */

      if(do_cluster) /* divide input MSA into clusters, and build CM from each cluster */
	{
	  if((status = MSADivide(msa, esl_opt_GetBoolean(go, "--call"), esl_opt_GetInteger(go, "--ctarget"), 
				 esl_opt_GetReal(go, "--cmindiff"), esl_opt_GetBoolean(go, "--corig"), &ncm, &cmsa)) != eslOK)
	    cm_Fail("MSADivide error (code: %d)\n", status);
	  esl_msa_Destroy(msa); /* we've copied the master msa into cmsa[ncm], we can delete this copy */
	}
      for(c = 0; c < ncm; c++)
	{
	  if(do_cluster) {
	      msa = cmsa[c];
	      if(esl_opt_GetString(go, "--cdump") != NULL) esl_msa_Write(cfg->cdfp, msa, cfg->fmt); 
	  }
	  /* Print some stuff about what we're about to do.
	   */
	  if (cfg->be_verbose) {
	    fprintf(cfg->ofp, "Alignment:           %s\n",  msa->name);
	    fprintf(cfg->ofp, "Number of sequences: %d\n",  msa->nseq);
	    fprintf(cfg->ofp, "Number of columns:   %d\n",  msa->alen);
	    if(esl_opt_GetString(go, "--rsearch") != NULL)
	      printf ("RIBOSUM Matrix:      %s\n",  cfg->fullmat->name);
	    fputs("", cfg->ofp);
	    fflush(stdout);
	  }
	  puts("");
	  
	  if ((status = process_workunit(go, cfg, errbuf,            msa, &cm, NULL)) != eslOK) cm_Fail(errbuf);
	  if ((status = output_result(   go, cfg, errbuf, cfg->nali, msa,  cm))       != eslOK) cm_Fail(errbuf);
	  
	  if (cfg->be_verbose) {
	    puts("");
	    SummarizeCM(stdout, cm);  
	  }
	  
	  FreeCM(cm);
	  fflush(cfg->cmfp);
	  puts("//\n");
	  if(do_cluster && cmsa != NULL)
	    esl_msa_Destroy(cmsa[c]);
	}
    }
  if(do_cluster) free(cmsa);
  return;
}

/* A work unit consists of one multiple alignment, <msa>.
 * The job is to turn it into a new CM, returned in <*ret_cm>.
 */
static int
process_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t **ret_cm, Parsetree_t ***opt_tr)
{
  CM_t *cm = NULL;
  int status;

  if ((status =  check_and_clean_msa    (go, cfg, errbuf, msa))                         != eslOK) goto ERROR;
  if ((status =  set_relative_weights   (go, cfg, errbuf, msa))                         != eslOK) goto ERROR;
  if ((status =  build_model            (go, cfg, errbuf, msa, &cm, opt_tr))            != eslOK) goto ERROR;
  if ((status =  set_model_name         (go, cfg, errbuf, msa, cm))                     != eslOK) goto ERROR;
  if ((status =  set_effective_seqnumber(go, cfg, errbuf, msa, cm, cfg->pri))           != eslOK) goto ERROR;
  if ((status =  parameterize           (go, cfg, errbuf, cm, cfg->pri))                != eslOK) goto ERROR;

  *ret_cm = cm;
  return eslOK;

 ERROR:
  FreeCM(cm);
  *ret_cm = NULL;
  if (opt_tr != NULL) *opt_tr = NULL;
  return status;
}

static int
output_result(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, int msaidx, ESL_MSA *msa, CM_t *cm)
{
  int status;

  /* Special case: output the tabular results header. 
   * Arranged this way to keep the two fprintf()'s close together in the code,
   * so we can keep the data and labels properly sync'ed.
   */
  if (msa == NULL && ! cfg->be_verbose) 
    {
      fprintf(cfg->ofp, "# %3s %-20s %5s %5s %5s\n", "idx", "name",                 "nseq",  "alen",  "M");
      fprintf(cfg->ofp, "#%4s %-20s %5s %5s %5s\n", "----", "--------------------", "-----", "-----", "-----");
      return eslOK;
    }

  if ((status = cm_Validate(cm, 0.0001, errbuf))  != eslOK) return status;
  if ((status = CMFileWrite(cfg->cmfp, cm, esl_opt_GetBoolean(go, "--binary"))) != eslOK) return status;
  if (! cfg->be_verbose)	/* tabular output */
    {                    /* #   name nseq alen M */
      fprintf(cfg->ofp, "%-5d %-20s %5d %5d %5d\n",
	      msaidx,
	      (msa->name != NULL) ? msa->name : "",
	      msa->nseq,
	      msa->alen,
	      cm->clen);
    }
  return eslOK;
}

/* check_and_clean_msa
 * Ensure we can build a CM from the MSA.
 * This requires it has a name.
 */
static int
check_and_clean_msa(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa)
{
  int status;
  if (cfg->be_verbose) {
    fprintf(cfg->ofp, "%-40s ... ", "Checking MSA");  
    fflush(cfg->ofp); 
  }

  if (esl_opt_GetBoolean(go, "--rf") && msa->rf == NULL) 
    ESL_FAIL(eslFAIL, errbuf, "Alignment has no reference coord annotation.\n");
  if (msa->ss_cons == NULL) 
    ESL_FAIL(eslFAIL, errbuf, "Alignment did not contain consensus structure annotation.\n");
  if (! clean_cs(msa->ss_cons, msa->alen))
    ESL_FAIL(eslFAIL, errbuf, "Failed to parse consensus structure annotation\n");
  if (esl_opt_GetBoolean(go, "--ignorant")) strip_wuss(msa->ss_cons);
  
  /* MSA better have a name, we named it before */
  if(msa->name == NULL) goto ERROR;

  if (cfg->be_verbose) fprintf(cfg->ofp, "done.\n");
  return eslOK;

 ERROR:
  return status;
}

/* set_relative_weights():
 * Set msa->wgt vector, using user's choice of relative weighting algorithm.
 */
static int
set_relative_weights(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa)
{
  if (cfg->be_verbose) {
    fprintf(cfg->ofp, "%-40s ... ", "Relative sequence weighting");  
    fflush(cfg->ofp); 
  }

  if      (esl_opt_GetBoolean(go, "--wnone"))                  esl_vec_DSet(msa->wgt, msa->nseq, 1.);
  else if (esl_opt_GetBoolean(go, "--wgiven"))                 ;
  else if (msa->nseq >= esl_opt_GetInteger(go, "--pbswitch"))  esl_msaweight_PB(msa);
  else if (esl_opt_GetBoolean(go, "--wpb"))                    esl_msaweight_PB(msa);
  else if (esl_opt_GetBoolean(go, "--wgsc"))                   esl_msaweight_GSC(msa);
  else if (esl_opt_GetBoolean(go, "--wblosum"))                esl_msaweight_BLOSUM(msa, esl_opt_GetReal(go, "--wid"));

  if (cfg->be_verbose) fprintf(cfg->ofp, "done.\n");
  return eslOK;
}


/* build_model():
 * Given <msa>, collect counts;
 * upon return, <*ret_cm> is newly allocated and contains
 * relative-weighted observed counts.
 * Optionally, caller can request an array of inferred parsetrees for
 * the <msa> too.
 */
static int
build_model(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t **ret_cm, Parsetree_t ***opt_tr)
{
  int status;
  Parsetree_t     *tr;
  Parsetree_t     *mtr;
  int idx;
  CM_t *cm;
  char *aseq;                   

  if (cfg->be_verbose) {
    fprintf(cfg->ofp, "%-40s ... ", "Constructing model architecture"); 
    fflush(cfg->ofp);
  }

  HandModelmaker(msa, esl_opt_GetBoolean(go, "--rf"), esl_opt_GetReal(go, "--gapthresh"), &cm, &mtr);
  printf("done.\n");
  
  /* set the CM's null model, if rsearch mode, use the bg probs used to calc RIBOSUM */
  if(esl_opt_GetString(go, "--rsearch") != NULL) CMSetNullModel(cm, cfg->fullmat->g); 
  else CMSetNullModel(cm, cfg->null); 
  
  /* if we're using RSEARCH emissions (--rsearch) set the flag */
  if(esl_opt_GetString(go, "--rsearch") != NULL) cm->flags |= CM_RSEARCHEMIT;

  /* rebalance CM */
  if(!esl_opt_GetBoolean(go, "--nobalance"))
    {
      CM_t *new;
      new = CMRebalance(cm);
      FreeCM(cm);
      cm = new;
    }
  /* get counts */
  for (idx = (int) esl_opt_GetBoolean(go, "--treeforce"); idx < msa->nseq; idx++)
    {
      ESL_ALLOC(aseq, (msa->alen+1) * sizeof(char));
      esl_abc_Textize(msa->abc, msa->ax[idx], msa->alen, aseq);
      tr = Transmogrify(cm, mtr, msa->ax[idx], aseq, msa->alen);
      ParsetreeCount(cm, tr, msa->ax[idx], msa->wgt[idx]);
      FreeParsetree(tr);
      free(aseq);
    }
  cm->nseq     = msa->nseq;
  cm->eff_nseq = msa->nseq;
  if(esl_opt_GetBoolean(go, "--treeforce"))
    {
      cm->nseq--;
      cm->eff_nseq--;
    }      

  /* ensure the dual insert states we will detach were populated with 0 counts */
  if(!(esl_opt_GetBoolean(go, "--nodetach")))
    {
      printf("%-40s ... ", "Finding and checking dual inserts");
      cm_find_and_detach_dual_inserts(cm, 
				      TRUE,   /* Do check (END_E-1) insert states have 0 counts */
				      FALSE); /* Don't detach the states yet, wait til CM is priorified */
    }
  /* set the EL self transition probability */
  cm->el_selfsc = sreLOG2(esl_opt_GetReal(go, "--elself"));
  printf("done.\n");
  
  /* Before converting to probabilities, 
   * save a count vector file, if asked.
   * Used primarily for making data files for training priors.
   */
  if (esl_opt_GetString(go, "--cfile") != NULL) {
    printf("%-40s ... ", "Saving count vector file"); fflush(stdout);
    if (! save_countvectors(esl_opt_GetString(go, "--cfile"), cm)) printf("[FAILED]\n");
    else                                printf("done. [%s]\n", esl_opt_GetString(go, "--cfile"));
  }
  FreeParsetree(mtr);
  *ret_cm = cm;
  return eslOK;

 ERROR:
  if (cfg->be_verbose) fprintf(cfg->ofp, "FAILED.\n");
  return status;
}


/* set_model_name()
 * Give the model a name based on the MSA name.
 * 
 * if msa->name is unavailable, or -n was used,
 * a fatal error is thrown. 
 *
 * note: This is much simpler than how HMMER3 does
 *       this. The reason is that the --ctarg --cmindiff
 *       cluster options produce N > 1 CM per MSA,
 *       which are named <msa->name>.1 .. <msa->name>.N.
 * 
 */
static int
set_model_name(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t *cm)
{
  int status;

  if (cfg->be_verbose) {
    fprintf(cfg->ofp, "%-40s ... ", "Set model name");
    fflush(cfg->ofp);
  }

  if(cm_SetName(cm, msa->name) != eslOK) goto ERROR;
  if (cfg->be_verbose) fprintf(cfg->ofp, "done. [%s]\n", cm->name);
  return eslOK;

 ERROR:
  if (cfg->be_verbose) fprintf(cfg->ofp, "FAILED.\n");
  return status;
}

/* set_effective_seqnumber()
 * Incept:    EPN, Fri Jul 27 10:38:11 2007
 * <cm> comes in with weighted observed counts. It goes out with
 * those observed counts rescaled to sum to the "effective sequence
 * number". 
 *
 * <prior> is needed because we may need to parameterize test models
 * looking for the right relative entropy. (for --eent, the default)
 *
 * Based on HMMER3's hmmbuild func of same name, we don't allow
 * --eset or --eclust here though.
 */
static int
set_effective_seqnumber(const ESL_GETOPTS *go, const struct cfg_s *cfg,
			char *errbuf, ESL_MSA *msa, CM_t *cm, const Prior_t *pri)
{
  double neff;

  printf("%-40s ... ", "Set effective sequence number");
  fflush(stdout);

  if      (esl_opt_GetBoolean(go, "--enone") == TRUE) 
    {
      neff = msa->nseq;
      printf("done. [--enone: neff=nseq=%d]\n", msa->nseq);
    }
  else if (esl_opt_GetBoolean(go, "--eent") == TRUE)
    {
      double etarget; 
      int clen = 0;
      int nd;
      for(nd = 0; nd < cm->nodes; nd++)
	{
	  if(cm->ndtype[nd] == MATP_nd) clen += 2;
	  else if(cm->ndtype[nd] == MATL_nd) clen += 1;
	  else if(cm->ndtype[nd] == MATR_nd) clen += 1;
	}
      if (esl_opt_IsDefault(go, "--ere")) etarget = default_target_relent(cm->abc, clen, esl_opt_GetReal(go, "--eX"));
      else                                etarget = esl_opt_GetReal(go, "--ere");

      neff = CM_Eweight_RE(cm, pri, (float) msa->nseq, etarget, cm->null);
      cm->eff_nseq = neff;
      CMRescale(cm, neff / (float) msa->nseq);
      printf("done. [etarget %.2f bits; neff %.2f]\n", etarget, neff);
    }
  return eslOK;
}

/* parameterize()
 * Converts counts to probability parameters.
 */
static int
parameterize(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, const Prior_t *prior)
{

  if (cfg->be_verbose){
    fprintf(cfg->ofp, "%-40s ... ", "Converting counts to probabilities"); 
    fflush(cfg->ofp);
  }
  PriorifyCM(cm, prior);
  if(esl_opt_GetString(go, "--rsearch") != NULL)
    {
      /*rsearch_CMProbifyEmissions(cm, fullmat); *//* use those probs to set CM probs from cts */
      /*debug_print_cm_params(cm);*/
    }
  
  if(!esl_opt_GetBoolean(go, "--nodetach")) /* Detach dual inserts where appropriate, if
					     * we get here we've already checked these states */
    {
      cm_find_and_detach_dual_inserts(cm, 
				      FALSE, /* Don't check states have 0 counts (they won't due to priors) */
				      TRUE); /* Detach the states by setting trans probs into them as 0.0   */
    }
  CMRenormalize(cm);
  CMLogoddsify(cm);

  if (cfg->be_verbose) fprintf(cfg->ofp, "done.\n");
  return eslOK;
}

/* default_target_relent()
 * Incept:    EPN, Tue Jul 10 10:13:43 2007
 *            based on HMMER3's hmmbuild.c:default_target_relent()
 *            SRE, Fri May 25 15:14:16 2007 [Janelia]
 *
 * Purpose:   Implements a length-dependent calculation of the target relative entropy
 *            per position, attempting to ensure that the information content of
 *            the model is high enough to find local alignments; but don't set it
 *            below a hard alphabet-dependent limit (CM_ETARGET).
 *            notes.
 *            
 * Args:      clen - consensus length (2*MATP + MATL + MATR)
 *            eX - X parameter: minimum total rel entropy target
 *
 */
static double
default_target_relent(const ESL_ALPHABET *abc, int clen, double eX)
{
  double etarget;
  /* HMMER3 default eX = 6.0 as of Tue Jul 10 2007
   */
  etarget = 6.* (eX + log((double) ((clen * (clen+1)) / 2)) / log(2.))    / (double)(2*clen + 4);

  switch (abc->type) {
  case eslRNA:    if (etarget < DEFAULT_ETARGET)   etarget = DEFAULT_ETARGET;   break;
  default:        esl_fatal("ERROR in default_target_relent(), alphabet not RNA!\n");
  }
  return etarget;
}

/* strip_wuss() remove all base pair info from a SS string
 */
void
strip_wuss(char *ss)
{
  char *s;
  for (s = ss; *s != '\0'; s++)
    if ((*s != '~') && (*s != '.'))
      *s = ':';
  return;
}

/* name_msa() 
 *
 * Give a MSA a name if it doesn't have one,
 * Naming rule is the suffixless name of the file it came from,
 * plus a "-<X>" with <X> = number MSA in the file.
 *
 * For example the 3rd MSA in file "alignments.stk" would be
 * named "alignments-3".
 */
int
name_msa(const ESL_GETOPTS *go, ESL_MSA *msa, int nali)
{
  int status;
  char *name = NULL;
  void *tmp;
  int n;
  char buffer[50];
  if(msa != NULL && msa->name == NULL)  
    {
      esl_FileTail(esl_opt_GetArg(go, 2), TRUE, &name); /* TRUE=nosuffix */
      if (name == NULL) cm_Fail("Error getting file tail of the MSA.\n");
      else {
	n  = strlen(name);
	sprintf(buffer, "-%d", (nali));
	n += strlen(buffer);
	ESL_RALLOC(name, tmp, sizeof(char)*(n+1));
	esl_strcat(&name, -1, buffer, (n+1));
	ESL_ALLOC(msa->name, sizeof(char) * (strlen(name)+1));
	strcpy(msa->name, name);
	free(name);
	if ((status = esl_strchop(msa->name, n)) != eslOK) goto ERROR;
      }
    }
  return eslOK;
 ERROR:
  if(name != NULL) free(name);
  return status;
}

/* Function: save_countvectors()
 * Date:     SRE, Tue May  7 16:21:10 2002 [St. Louis]
 *
 * Purpose:  Save emission count vectors to a file.
 *           Used to gather data for training Dirichlet priors.
 *
 * Args:     cfile  - name of file to save vectors to.
 *           cm     - a model containing counts (before probify'ing)
 *
 */
static int
save_countvectors(char *cfile, CM_t *cm)
{
  FILE *fp;
  int   v,x;

  /* Print emission counts */
  if ((fp = fopen(cfile, "w")) == NULL) return 0;
  for (v = 0; v < cm->M; v++)
    {
      if (cm->sttype[v] == MP_st || 
	  cm->sttype[v] == ML_st || 
	  cm->sttype[v] == MR_st) 
	{
	  fprintf(fp, "E\t%-7s ", UniqueStatetype(cm->stid[v]));
	  if (cm->sttype[v] == MP_st) {
	    for (x = 0; x < cm->abc->K*cm->abc->K; x++)
	      fprintf(fp, "%8.3f ", cm->e[v][x]);
	  } else {
	    for (x = 0; x < cm->abc->K; x++)
	      fprintf(fp, "%8.3f ", cm->e[v][x]);
	  }
	  fprintf(fp, "\n");
	}
    }

  /* Print transition counts */
  for (v = 0; v < cm->M; v++)
    {
      if(cm->sttype[v] != B_st && cm->sttype[v] != E_st)
	{
	  fprintf(fp, "T\t%-7s : %-2d", UniqueStatetype(cm->stid[v]), cm->ndtype[(cm->ndidx[v] + 1)]);
	  for (x = 0; x < cm->cnum[v]; x++)
	    {
	      fprintf(fp, "%8.3f ", cm->t[v][x]);
	    }
	  fprintf(fp, "\n");
	}
    }

  fclose(fp);
  return 1;
}

/* EPN 08.18.05
 * model_trace_info_dump()
 * Function: model_trace_info_dump
 *
 * Purpose:  Given a trace from a sequence used to create the model, 
 *           print the subsequence length rooted at each start state.  
 *           The sequence positions in a Parsetree_t tr
 *           returned from Transmogrify refer to aligned positions.
 *           We want subsequence lengths that refer to unaligned lengths.
 * 
 * Args:    ofp      - filehandle to print to
 *          cm       - the CM
 *          tr       - the parsetree (trace)
 *          aseq     - the aligned sequence the trace corresponds to
 * Returns: (void) 
 */

static void
model_trace_info_dump(FILE *ofp, CM_t *cm, Parsetree_t *tr, char *aseq)
{
  int status;
  int a, i, j, tpos, d, l, r;
  int *map;

  ESL_ALLOC(map, sizeof(int) * strlen(aseq));
  
  a=0;
  for (i = 0; i < strlen(aseq); i++)
    if (! esl_abc_CIsGap(cm->abc, aseq[i])) map[i] = a++;
    else map[i] = -1;

  for (tpos = 0; tpos < tr->n; tpos++)
    if(cm->sttype[tr->state[tpos]] == S_st)
      {
	l = tr->emitl[tpos]-1;
	r = tr->emitr[tpos]-1;
	i = map[l];
	j = map[r];
	/* tr->emitl[tpos]-1 might map to a gap (root node emits the gaps
	 * also). So we look for first residue that exists in the unaligned
	 * seq.  Then we do the same for j, looking backwards.
	 */ 
	while (i == -1)
	  i = map[++l];
	while (j == -1)
	  j = map[--r];

	d = j-i+1;
	/* assume ofp is open (probably not good) */
	fprintf(ofp, "state:%d d:%d\n", tr->state[tpos], d);
	/*fprintf(ofp, "state:%d d:%d i:%d j:%d emitl:%d emitr:%d\n", tr->state[tpos], d, i, j, tr->emitl[tpos], tr->emitr[tpos]);*/
      }
  free(map);

 ERROR:
  esl_fatal("Memory allocation error.");
}

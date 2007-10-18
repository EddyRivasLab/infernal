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

#include "funcs.h"		/* external functions                   */
#include "structs.h"		/* data structures, macros, #define's   */

#define WGTOPTS "--wgsc,--wblosum,--wpb,--wnone,--wgiven"      /* Exclusive options for relative weighting                    */
#define EFFOPTS "--eent,--enone"               /* Exclusive options for effective sequence number calculation */
#define ALPHOPTS "--rna,--dna"                 /* Exclusive options for specifiying input MSA alphabet*/

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs       incomp  help  docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",   1 },
  { "-n",        eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL,        NULL, "name the CM(s) <s>",                     1 },
  { "-o",        eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "direct summary output to file <f>, not stdout", 1 },
  { "-A",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "append this CM to <cmfile>",             1 },
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "force; allow overwriting of <cmfile>",   1 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "verbose; print out extra information",   1 },
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
/* Debugging/experimentation */
  { "--nobalance",eslARG_NONE,  FALSE, NULL, NULL,      NULL,      NULL,        NULL, "don't rebalance the CM; number in strict preorder", 6 },
  { "--regress",  eslARG_STRING, NULL, NULL, NULL,      NULL,      NULL,        NULL, "save regression test information to file <s>", 6 },  
  { "--ignorant", eslARG_NONE,  FALSE, NULL, NULL,      NULL,      NULL,        NULL, "strip the structural info from input alignment", 6 },
/* Customizing null model or priors */
  { "--null",    eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL, "--rsearch", "read null (random sequence) model from file <s>", 7 },
  { "--prior",   eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL, "--rsearch", "read priors from file <s>", 7 },
/* Building multiple CMs after clustering input MSA */
  { "--ctarget", eslARG_INT,    "0",   NULL, "n>=0",    NULL,      NULL,    "--call", "build (at most) <n> CMs by partitioning MSA into <n> clusters", 8 },
  { "--cmindiff",eslARG_REAL,   "0.",  NULL,"0.<=x<=1.",NULL,      NULL,    "--call", "min difference b/t 2 clusters is <x>, each cluster -> CM", 8 }, 
  { "--call",    eslARG_NONE,  FALSE,  NULL, NULL,      NULL,      NULL,        NULL, "build a separate CM from every seq in MSA", 8 },
  { "--corig",   eslARG_NONE,  FALSE,  NULL, NULL,      NULL,      NULL,        NULL, "build an additional CM from the original, full MSA", 8 }, 
  { "--cdump",   eslARG_STRING, NULL,  NULL, NULL,      NULL,      NULL,        NULL, "dump an MSA for each cluster (CM) to file <s>", 8 },
/* Refining the seed alignment */
  { "--refine",  eslARG_OUTFILE, NULL,  NULL, NULL,       NULL,  NULL,       NULL, "refine input aln w/Expectation-Maximization, save to <s>", 9 },
  { "--gibbs",   eslARG_NONE,   FALSE,  NULL, NULL,       NULL,"--refine",   NULL, "w/--refine, use Gibbs sampling instead of EM", 9 },
  { "--seed",    eslARG_INT,     NULL,  NULL, "n>0",      NULL,"--gibbs",    NULL, "w/--gibbs, set random number generator seed to <n>",  9 },
  { "--hbanded", eslARG_NONE,   FALSE,  NULL, NULL,       NULL,"--refine",   NULL, "accelerate --refine using HMM banded CYK aln algorithm", 9 },
  { "--tau",     eslARG_REAL,   "1E-7", NULL, "0<x<1",    NULL,"--hbanded",  NULL, "set tail loss prob for --hbanded to <x>", 9 },
  { "--sub",     eslARG_NONE,   FALSE,  NULL, NULL,       NULL,"--refine",   NULL, "w/--refine, use sub CM for columns b/t HMM start/end points", 9 },
  { "--local",   eslARG_NONE,   FALSE,  NULL, NULL,       NULL,"--refine",   NULL, "w/--refine, align locally w.r.t the model", 9 },
/* Selecting the input MSA alphabet rather than autoguessing it */
  { "--rna",     eslARG_NONE,   FALSE, NULL, NULL,   ALPHOPTS,    NULL,     NULL, "input alignment is RNA sequence data", 10},
  { "--dna",     eslARG_NONE,   FALSE, NULL, NULL,   ALPHOPTS,    NULL,     NULL, "input alignment is DNA sequence data", 10},
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
  FILE         *trfp;           /* if --refine, output file handle for dumping refined MSAs */

  int           be_verbose;	/* standard verbose output, as opposed to one-line-per-CM summary */
  int           nali;		/* which # alignment this is in file (only valid in serial mode)   */
  ESL_RANDOMNESS *r;            /* source of randomness, only created if --gibbs enabled */
};

static char usage[]  = "[-options] <cmfile output> <alignment file>";
static char banner[] = "build RNA covariance model(s) from alignment";

static int  init_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);

static void  master (const ESL_GETOPTS *go, struct cfg_s *cfg);

static int    process_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t **ret_cm, Parsetree_t ***ret_msa_tr);
static int    output_result(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, int msaidx, ESL_MSA *msa, CM_t *cm);
static int    check_and_clean_msa(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa);
static int    set_relative_weights(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa);
static int    build_model(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t **ret_cm, Parsetree_t ***ret_msa_tr);
static int    set_model_name(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t *cm);
static int    set_effective_seqnumber(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t *cm, const Prior_t *pri);
static int    parameterize(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, const Prior_t *prior);
static int    name_msa(const ESL_GETOPTS *go, ESL_MSA *msa, int nali);
static double default_target_relent(const ESL_ALPHABET *abc, int M, double eX);
static int    save_countvectors(char *cfile, CM_t *cm);
static void   strip_wuss(char *ss);
static int    refine_msa(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *orig_cm, ESL_MSA *input_msa, Parsetree_t **input_msa_tr, CM_t **ret_cm, ESL_MSA **ret_msa, Parsetree_t ***ret_tr, int *ret_niter);
static int    get_unaln_seqs_from_msa(const ESL_MSA *msa, ESL_SQ ***ret_sq);
static int    convert_parsetrees_to_unaln_coords(Parsetree_t **tr, ESL_MSA *msa);
static int    initialize_cm(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm);

/*static void   model_trace_info_dump(FILE *ofp, CM_t *cm, Parsetree_t *tr, char *aseq);*/


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
      puts("\nexpert options for refining the input alignment:");
      esl_opt_DisplayHelp(stdout, go, 9, 2, 80);
      puts("\n options for selecting alphabet rather than guessing it:");
      esl_opt_DisplayHelp(stdout, go, 10, 2, 80);
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
  cfg.trfp       = NULL;	           /* opened (possibly) in init_cfg() */
  cfg.r          = NULL;	           /* created (possibly) in init_cfg() */

  if (esl_opt_GetBoolean(go, "-v")) cfg.be_verbose = TRUE;
  else                              cfg.be_verbose = FALSE;        
  cfg.nali       = 0;		           

  /* Start timing; do work; stop timing.*/
  esl_stopwatch_Start(w);
  master(go, &cfg);
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(cfg.ofp, w, "# CPU time: ");

  /* Clean up the cfg. 
   */
  if (! esl_opt_IsDefault(go, "-o")) { fclose(cfg.ofp); }
  if (cfg.afp   != NULL) esl_msafile_Close(cfg.afp);
  if (cfg.abc   != NULL) esl_alphabet_Destroy(cfg.abc);
  if (cfg.cmfp  != NULL) fclose(cfg.cmfp);
  if (cfg.pri   != NULL) Prior_Destroy(cfg.pri);
  if (cfg.null  != NULL) free(cfg.null);
  if (cfg.cdfp  != NULL) fclose(cfg.cdfp);
  if (cfg.trfp  != NULL) fclose(cfg.trfp);
  if (cfg.r     != NULL) esl_randomness_Destroy(cfg.r);

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
  if      (esl_opt_GetBoolean(go, "--rna")) type = eslRNA;
  else if (esl_opt_GetBoolean(go, "--dna")) type = eslDNA;
  else { 
    status = esl_msafile_GuessAlphabet(cfg->afp, &type);
    if (status == eslEAMBIGUOUS)    ESL_FAIL(status, errbuf, "Failed to guess the bio alphabet used in %s.\nUse --rna option to specify it.", cfg->alifile);
    else if (status == eslEFORMAT)  ESL_FAIL(status, errbuf, "Alignment file parse failed: %s\n", cfg->afp->errbuf);
    else if (status == eslENODATA)  ESL_FAIL(status, errbuf, "Alignment file %s is empty\n", cfg->alifile);
    else if (status != eslOK)       ESL_FAIL(status, errbuf, "Failed to read alignment file %s\n", cfg->alifile);
  }
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

  /* if --cdump enabled, open output file for cluster MSAs */
  if (esl_opt_GetString(go, "--cdump") != NULL)
    {
      if ((cfg->cdfp = fopen(esl_opt_GetString(go, "--cdump"), "w")) == NULL)
	cm_Fail("Failed to open output file %s for writing MSAs to", esl_opt_GetString(go, "--cdump"));
    }

  /* if --refine enabled, open output file for refined MSAs */
  if (esl_opt_GetString(go, "--refine") != NULL)
    {
      if ((cfg->trfp = fopen(esl_opt_GetString(go, "--refine"), "w")) == NULL)
	cm_Fail("Failed to open output file %s for writing MSAs to", esl_opt_GetString(go, "--refine"));
    }
  /* if --gibbs enabled, open output file for refined MSAs, and seed RNG */
  if(esl_opt_GetBoolean(go, "--gibbs"))
    {
      /* create RNG */
      if (! esl_opt_IsDefault(go, "--seed")) 
	cfg->r = esl_randomness_Create((long) esl_opt_GetInteger(go, "--seed"));
      else cfg->r = esl_randomness_CreateTimeseeded();
      if (cfg->r       == NULL) ESL_FAIL(eslEINVAL, errbuf, "Failed to create random number generator: probably out of memory");
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
  CM_t    *new_cm;
  Parsetree_t **tr;
  int          i = 0;
  int      niter = 0;
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
	  
	  /* msa -> cm */
	  if ((status = process_workunit(go, cfg, errbuf, msa, &cm, &tr)) != eslOK) cm_Fail(errbuf);
	  /* optionally, iterative over cm -> parsetrees -> msa -> cm ... until convergence, via EM (or Gibbs - not yet, but eventually) */
	  if (! esl_opt_IsDefault(go, "--refine")) {
	    if ((status = refine_msa(go, cfg, errbuf, cm, msa, tr, &new_cm, NULL, NULL, &niter)) != eslOK) cm_Fail(errbuf);
	    if (niter > 1) { /* if niter == 1, we didn't make a new CM (new_cm == cm), so we don't free it */
	      FreeCM(cm); 
	      cm = new_cm; 
	    } 
	  }	  
	  /* output cm */
	  if ((status = output_result(go, cfg, errbuf, cfg->nali, msa,  cm)) != eslOK) cm_Fail(errbuf);
	  
	  SummarizeCM(cfg->ofp, cm);  
	  if(cfg->ofp != stdout) SummarizeCM(stdout, cm);  
	  
	  FreeCM(cm);
	  fflush(cfg->cmfp);
	  puts("//\n");

	  if(tr != NULL) {
	    for(i = 0; i < msa->nseq; i++)
	      FreeParsetree(tr[i]);
	    free(tr);
	    tr = NULL;
	  }
	  esl_msa_Destroy(msa);
	}
    }
  if(do_cluster) free(cmsa);
  return;
}

/* A work unit consists of one multiple alignment, <msa>.
 * The job is to turn it into a new CM, returned in <*ret_cm>.
 * 
 */
static int
process_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t **ret_cm, Parsetree_t ***ret_msa_tr)
{
  CM_t *cm = NULL;
  int status;

  if ((status =  check_and_clean_msa    (go, cfg, errbuf, msa))                         != eslOK) goto ERROR;
  if ((status =  set_relative_weights   (go, cfg, errbuf, msa))                         != eslOK) goto ERROR;
  if ((status =  build_model            (go, cfg, errbuf, msa, &cm, ret_msa_tr))        != eslOK) goto ERROR;
  if ((status =  set_model_name         (go, cfg, errbuf, msa, cm))                     != eslOK) goto ERROR;
  if ((status =  set_effective_seqnumber(go, cfg, errbuf, msa, cm, cfg->pri))           != eslOK) goto ERROR;
  if ((status =  parameterize           (go, cfg, errbuf, cm, cfg->pri))                != eslOK) goto ERROR;
  
  *ret_cm = cm;
  return eslOK;

 ERROR:
  if(cm != NULL) FreeCM(cm);
  *ret_cm = NULL;
  if (ret_msa_tr != NULL) *ret_msa_tr = NULL;
  return status;
}

/* refine_msa: 
 * Refine the original (input) MSA using Expectation-Maximization or Gibbs sampling
 * by iterating over: build MSA of optimal parses, build CM from MSA,
 * until the summed scores of all parses converges.
 *
 * Note: input_msa_tr is modified, it's alignment coordinates are changed
 *       from aligned to unaligned.
 *
 */
static int
refine_msa(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *init_cm, ESL_MSA *input_msa, 
	   Parsetree_t **input_msa_tr, CM_t **ret_cm, ESL_MSA **ret_msa, Parsetree_t ***ret_tr, int *ret_niter)
{
  int              status;
  float            threshold   = 0.01;
  float            delta       = 1.;
  float            oldscore    = IMPOSSIBLE;
  float            totscore    = 0.;
  int              i           = 0;
  int              iter        = 0;
  ESL_SQ         **sq          = NULL;
  float           *sc          = NULL;
  int              nseq        = input_msa->nseq; 
  char            *msa_name    = NULL;
  CM_t            *cm          = NULL;
  seqs_to_aln_t   *seqs_to_aln = NULL; 
  ESL_MSA         *msa         = NULL;

  /* check contract */
  if(input_msa       == NULL) cm_Fail("in refine_msa, input_msa passed in as NULL");
  if(input_msa->name == NULL) cm_Fail("in refine_msa, input_msa must have a name");
  if(init_cm         == NULL) cm_Fail("in refine_msa, init_cm passed in as NULL");

  /* copy input MSA's name, we'll copy it to the MSA we create at each iteration */
  esl_strdup(input_msa->name, -1, &(msa_name));

  ESL_ALLOC(sc, sizeof(float) * nseq);
  esl_vec_FSet(sc, nseq, 0.);

  get_unaln_seqs_from_msa(input_msa, &sq); /* we need sqs for Parsetrees2Alignment */
  seqs_to_aln = CreateSeqsToAlnFromSq(sq, nseq, FALSE);

  /* determine scores of implicit parsetrees of input MSA seqs to initial CM */
  convert_parsetrees_to_unaln_coords(input_msa_tr, input_msa);
  if(cfg->be_verbose) fprintf(cfg->ofp, "iteration: %4d\n", iter);
  for(i = 0; i < nseq; i++) sc[i] = ParsetreeScore(init_cm, input_msa_tr[i], sq[i]->dsq, FALSE);
  oldscore = esl_vec_FSum(sc, nseq);
  fprintf(cfg->ofp, "iter: %4d (input alignment)     sc %10.4f (delta: N/A)\n", iter, oldscore);

  if(cfg->be_verbose) {
    fprintf(cfg->ofp, "INITIAL (input msa)\n");
    esl_msa_Write(stdout, input_msa, cfg->fmt); 
  }
     
  while(1)
    {
      iter++;
      if(iter == 1) { cm = init_cm; msa = input_msa; }
      
      /* 1. cm -> parsetrees */
      if(iter > 1) FreePartialSeqsToAln(seqs_to_aln, FALSE, TRUE, TRUE, TRUE, TRUE);
                                                  /* sq,    tr, cp9_tr, post, sc */ 
      /* initialize/configure CM, we may be doing HMM banded alignment for ex. */
      initialize_cm(go, cfg, errbuf, cm);
      actually_align_targets(cm, seqs_to_aln, NULL, NULL, 0, 0, (!cfg->be_verbose), cfg->r);
      
      /* sum parse scores and check for convergence */
      if(cfg->be_verbose) fprintf(cfg->ofp, "iteration: %4d\n", iter);
      totscore = esl_vec_FSum(seqs_to_aln->sc, nseq);
      delta    = (totscore - oldscore) / fabs(totscore);
      fprintf(cfg->ofp, "iter: %4d old sc %10.4f new sc %10.4f (delta: %10.4f)\n", iter, oldscore, totscore, delta);
      if(delta <= threshold && delta >= 0) break; /* only way out of while(1) loop */
      oldscore = totscore;

      /* 2. parsetrees -> msa */
      if( iter > 1) esl_msa_Destroy(msa);
      msa = NULL; /* even if iter == 1; we set msa to NULL, so we don't klobber input_msa */
      if((status = Parsetrees2Alignment(cm, cm->abc, seqs_to_aln->sq, NULL, seqs_to_aln->tr, nseq, FALSE, FALSE, &msa)) != eslOK) 
	cm_Fail("refine_msa() failed to make new MSA");
      esl_strdup(msa_name, -1, &(msa->name)); 
      esl_msa_Digitize(msa->abc, msa);

      if(cfg->be_verbose)
	esl_msa_Write(stdout, msa, cfg->fmt); 

      /* 3. msa -> cm */
      if(iter > 1) FreeCM(cm);
      cm = NULL; /* even if iter == 1; we set cm to NULL, so we don't klobber init_cm */
      if ((status = process_workunit(go, cfg, errbuf, msa, &cm, NULL))  != eslOK) cm_Fail(errbuf);

    }
  if(cfg->be_verbose) {
    fprintf(cfg->ofp, "FINAL (iter: %d)\n", iter);
    esl_msa_Write(stdout, msa, cfg->fmt); 
  }
  /* write out final alignment */
  esl_msa_Write(cfg->trfp, msa, cfg->fmt); 

  /* if CM was in local mode for aligning input MSA seqs, make it global so we can write it out */
  if((cm->flags & CM_LOCAL_BEGIN) || (cm->flags & CM_LOCAL_END))
    ConfigGlobal(cm);

  /* clean up */
  if(ret_cm == NULL) cm_Fail("ret_cm is NULL.");
  *ret_cm = cm;

  if(iter > 1) { /* if iter == 1, msa == init_msa, we don't want to free it, or overwrite it */
    if(ret_msa == NULL) esl_msa_Destroy(msa); 
    else *ret_msa = msa;
  }

  if(ret_tr == NULL) { FreeSeqsToAln(seqs_to_aln); }
  else { 
    FreePartialSeqsToAln(seqs_to_aln, TRUE, FALSE, TRUE,   TRUE, TRUE);
                                   /* sq,   tr,    cp9_tr, post, sc */ 
    free(seqs_to_aln->tr);
    free(seqs_to_aln);
  }

  *ret_niter = iter;

  free(sc);
  free(msa_name);

  return eslOK;

 ERROR:
  /* no cleanup, we die */
  cm_Fail("in refine_msa(), error, status: %d\n", status);
  return status; /* NEVERREACHED */
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
  if (! clean_cs(msa->ss_cons, msa->alen, FALSE))
    ESL_FAIL(eslFAIL, errbuf, "Failed to parse consensus structure annotation\n");
  if (esl_opt_GetBoolean(go, "--ignorant")) strip_wuss(msa->ss_cons);
  
  /* MSA better have a name, we named it before */
  if(msa->name == NULL) { sprintf(errbuf, "MSA is nameless"); goto ERROR; }

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
build_model(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t **ret_cm, Parsetree_t ***ret_msa_tr)
{
  int status;
  Parsetree_t     **tr;
  Parsetree_t     *mtr;
  int idx;
  CM_t *cm;
  char *aseq;                   

  if (cfg->be_verbose) {
    fprintf(cfg->ofp, "%-40s ... ", "Constructing model architecture"); 
    fflush(cfg->ofp);
  }

  HandModelmaker(msa, esl_opt_GetBoolean(go, "--rf"), esl_opt_GetReal(go, "--gapthresh"), &cm, &mtr);
  if(cfg->be_verbose) fprintf(cfg->ofp, "done.\n");
  
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
  ESL_ALLOC(tr, sizeof(Parsetree_t *) * (msa->nseq));
  for (idx = 0; idx < msa->nseq; idx++)
    {
      ESL_ALLOC(aseq, (msa->alen+1) * sizeof(char));
      esl_abc_Textize(msa->abc, msa->ax[idx], msa->alen, aseq);
      tr[idx] = Transmogrify(cm, mtr, msa->ax[idx], aseq, msa->alen);
      ParsetreeCount(cm, tr[idx], msa->ax[idx], msa->wgt[idx]);
      free(aseq);
    }
  if(ret_msa_tr == NULL) {
    for(idx = 0; idx < msa->nseq; idx++)
      FreeParsetree(tr[idx]);
    free(tr);
    tr = NULL;
  }
  else *ret_msa_tr = tr;

  cm->nseq     = msa->nseq;
  cm->eff_nseq = msa->nseq;

  /* ensure the dual insert states we will detach were populated with 0 counts */
  if(!(esl_opt_GetBoolean(go, "--nodetach")))
    {
      if(cfg->be_verbose) fprintf(cfg->ofp, "%-40s ... ", "Finding and checking dual inserts");
      cm_find_and_detach_dual_inserts(cm, 
				      TRUE,   /* Do check (END_E-1) insert states have 0 counts */
				      FALSE); /* Don't detach the states yet, wait til CM is priorified */
    }
  /* set the EL self transition probability */
  cm->el_selfsc = sreLOG2(esl_opt_GetReal(go, "--elself"));
  if(cfg->be_verbose) fprintf(cfg->ofp, "done.\n");
  
  /* Before converting to probabilities, 
   * save a count vector file, if asked.
   * Used primarily for making data files for training priors.
   */
  if (esl_opt_GetString(go, "--cfile") != NULL) {
    fprintf(cfg->ofp, "%-40s ... ", "Saving count vector file"); fflush(stdout);
    if (! save_countvectors(esl_opt_GetString(go, "--cfile"), cm)) fprintf(cfg->ofp, "[FAILED]\n");
    else                                fprintf(cfg->ofp, "done. [%s]\n", esl_opt_GetString(go, "--cfile"));
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

  if(cfg->be_verbose) fprintf(cfg->ofp, "%-40s ... ", "Set effective sequence number");
  fflush(stdout);

  if      (esl_opt_GetBoolean(go, "--enone") == TRUE) 
    {
      neff = msa->nseq;
      if(cfg->be_verbose) fprintf(cfg->ofp, "done. [--enone: neff=nseq=%d]\n", msa->nseq);
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
      if(cfg->be_verbose) fprintf(cfg->ofp, "done. [etarget %.2f bits; neff %.2f]\n", etarget, neff);
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

/* get_unaln_seqs_from_msa
 * Given a digitized MSA, allocate and create digitized versions
 * of the unaligned sequences within it.
 */
static int
get_unaln_seqs_from_msa(const ESL_MSA *msa, ESL_SQ ***ret_sq)
{
  int status;
  ESL_DSQ *uadsq = NULL;
  ESL_SQ **sq    = NULL;
  int nongap_len = 0;
  int i          = 0;
  int apos       = 1;
  int uapos      = 1;

  /* contract check */
  if(! (msa->flags & eslMSA_DIGITAL)) cm_Fail("get_unaln_seqs_from_msa() msa is not digitized.\n");

  ESL_ALLOC(sq, sizeof(ESL_SQ *) * msa->nseq);

  for (i = 0; i < msa->nseq; i++)
    {
      nongap_len = 0;
      for(apos = 1; apos <= msa->alen; apos++)
	nongap_len += (! esl_abc_XIsGap(msa->abc, msa->ax[i][apos]));
      ESL_ALLOC(uadsq, sizeof(ESL_DSQ) * (nongap_len + 2));
      uadsq[0] = uadsq[(nongap_len+1)] = eslDSQ_SENTINEL;

      uapos = 1;
      for(apos = 1; apos <= msa->alen; apos++)
	if(! esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) 
	  uadsq[uapos++] = msa->ax[i][apos];
      
      sq[i] = esl_sq_CreateDigitalFrom(msa->abc, msa->sqname[i], uadsq, nongap_len, NULL, NULL, NULL); 
      free(uadsq);
    }
  *ret_sq = sq;
  return eslOK;
  
 ERROR:
  cm_Fail("memory allocation error.");
  return status; /* NEVERREACHED */
}

/* convert_parsetrees_to_unaln_coords()
 *
 * Given a digitized MSA <msa> and parsetrees <tr> that correspond to 
 *  the ALIGNED coordinates in <msa>, modify tr[i]->emitl and tr[i]->emitr 
 * so they correspond with UNALIGNED coordinates. Written so we can call 
 * Parsetrees2Alignment() to make a  new msa, that will replace <msa> for 
 * training a CM.
 */
static int 
convert_parsetrees_to_unaln_coords(Parsetree_t **tr, ESL_MSA *msa)
{
  int status;
  int **map = NULL;
  int     i = 0;
  int     x = 0;
  int apos  = 1;
  int uapos = 1;
  /* contract check */
  if(! (msa->flags & eslMSA_DIGITAL)) cm_Fail("get_unaln_seqs_from_msa() msa is not digitized.\n");

  /* For each seq in the MSA, map the aligned sequences coords to 
   * the unaligned coords, we stay in digitized seq coords (1..alen)
   */
  ESL_ALLOC(map,   sizeof(int *)  * msa->nseq);
  for (i = 0; i < msa->nseq; i++) {
    ESL_ALLOC(map[i],   sizeof(int)  * (msa->alen+1));
    map[i][0] = -1; /* invalid */
    uapos = 1;
    for(apos = 1; apos <= msa->alen; apos++)
      if (!esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) map[i][apos] = uapos++;
      else                                             map[i][apos] = -1;
  }
  for (i = 0; i < msa->nseq; i++) {
    /* tr[i] is in alignment coords, convert it to unaligned coords, */
    for(x = 0; x < tr[i]->n; x++) {
      if(tr[i]->emitl[x] != -1)  tr[i]->emitl[x] = map[i][tr[i]->emitl[x]];
      if(tr[i]->emitr[x] != -1)  tr[i]->emitr[x] = map[i][tr[i]->emitr[x]];
    }
  }
  for (i = 0; i < msa->nseq; i++) free(map[i]);
  free(map);

  return eslOK;

 ERROR:
  cm_Fail("memory allocation error.");
  return status; /* NEVERREACHED */
}


/* initialize_cm()
 * Setup the CM based on the command-line options/defaults.
 * Configures the CM with a ConfigCM() call at end.
 */
static int
initialize_cm(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  /* set up params/flags/options of the CM */
  cm->tau    = esl_opt_GetReal(go, "--tau");  /* this will be DEFAULT_TAU unless changed at command line */

  /* update cm->align->opts */
  if(esl_opt_GetBoolean(go, "--gibbs"))       cm->align_opts  |= CM_ALIGN_SAMPLE;
  if(esl_opt_GetBoolean(go, "--hbanded")) {
    cm->align_opts  |= CM_ALIGN_HBANDED;
    cm->align_opts  |= CM_ALIGN_NOSMALL; 
  }
  if(esl_opt_GetBoolean(go, "--sub"))         cm->align_opts  |= CM_ALIGN_SUB;

  /* update cm->config_opts */
  if(esl_opt_GetBoolean(go, "--local"))
    {
      cm->config_opts |= CM_CONFIG_LOCAL;
      cm->config_opts |= CM_CONFIG_HMMLOCAL;
      cm->config_opts |= CM_CONFIG_HMMEL;
    }

  /* finally, configure the CM for alignment based on cm->config_opts and cm->align_opts.
   * this may make a cp9 HMM, for example.
   */
  ConfigCM(cm, NULL, NULL); 

  return eslOK;
}

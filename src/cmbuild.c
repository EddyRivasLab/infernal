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
#include "esl_distance.h"
#include "esl_dmatrix.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msaweight.h"
#include "esl_msacluster.h"
#include "esl_stack.h"
#include "esl_stopwatch.h"
#include "esl_tree.h"
#include "esl_vectorops.h"

#include "funcs.h"		/* external functions                   */
#include "structs.h"		/* data structures, macros, #define's   */

#define WGTOPTS "--wgsc,--wblosum,--wpb,--wnone,--wgiven"      /* Exclusive options for relative weighting                    */
#define EFFOPTS "--eent,--enone"               /* Exclusive options for effective sequence number calculation */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs       incomp  help  docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",   1 },
  { "-n",        eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL,        NULL, "name the CM(s) <s>, (only if single aln in file)", 1 },
  { "-A",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "append this CM to <cmfile>",             1 },
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "force; allow overwriting of <cmfile>",   1 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "be verbose with output", 1 },
  { "--iins",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "allow informative insert emissions, do not zero them", 1 },
  { "--Wbeta",   eslARG_REAL,   "1E-7",NULL, "x>0.0000000000000001",NULL,NULL,  NULL, "set tail loss prob for calc'ing W (max size of a hit) to <x>", 1 },
  { "--devhelp", eslARG_NONE,   NULL,  NULL, NULL,      NULL,      NULL,        NULL, "show list of undocumented developer options", 1 },
/* Expert model construction options */
  { "--rsearch", eslARG_INFILE, NULL,  NULL, NULL,      NULL,      NULL,        NULL,  "use RSEARCH parameterization with RIBOSUM matrix file <s>", 2 }, 
  { "--binary",  eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "save the model(s) in binary format",     2 },
  { "--rf",      eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,  "--rsearch", "use reference coordinate annotation to specify consensus", 2 },
  { "--gapthresh",eslARG_REAL,  "0.5", NULL, "0<=x<=1", NULL,      NULL,  "--rsearch", "fraction of gaps to allow in a consensus column [0..1]", 2 },
  { "--ignorant", eslARG_NONE,  FALSE, NULL, NULL,      NULL,      NULL,        NULL, "strip the structural info from input alignment", 2 },
/* Alternate relative sequence weighting strategies */
  /* --wme not implemented in Infernal yet (b/c it's not in HMMER3) */
  { "--wgsc",    eslARG_NONE,"default",NULL, NULL,    WGTOPTS,    NULL,      NULL, "Gerstein/Sonnhammer/Chothia tree weights",         3},
  { "--wblosum", eslARG_NONE,  FALSE,  NULL, NULL,    WGTOPTS,    NULL,      NULL, "Henikoff simple filter weights",                   3},
  { "--wpb",     eslARG_NONE,  FALSE,  NULL, NULL,    WGTOPTS,    NULL,      NULL, "Henikoff position-based weights",                  3},
  { "--wnone",   eslARG_NONE,  FALSE,  NULL, NULL,    WGTOPTS,    NULL,      NULL, "don't do any relative weighting; set all to 1",    3},
  { "--wgiven",  eslARG_NONE,  FALSE,  NULL, NULL,    WGTOPTS,    NULL,      NULL, "use weights as given in MSA file",                 3},
  { "--pbswitch",eslARG_INT,  "5000",  NULL,"n>0",       NULL,    NULL,      NULL, "set failover to efficient PB wgts at > <n> seqs",  3},
  { "--wid",     eslARG_REAL, "0.62",  NULL,"0<=x<=1",   NULL,"--wblosum",   NULL, "for --wblosum: set identity cutoff",               3},
/* Alternate effective sequence weighting strategies */
  { "--eent",    eslARG_NONE,"default",NULL, NULL,    EFFOPTS,    NULL,      NULL, "adjust eff seq # to achieve relative entropy target", 4},
  { "--enone",   eslARG_NONE,  FALSE,  NULL, NULL,    EFFOPTS,    NULL,      NULL, "no effective seq # weighting: just use nseq",         4},
  { "--ere",     eslARG_REAL,  NULL,   NULL,"x>0",       NULL, "--eent",     NULL, "for --eent: set CM target relative entropy to <x>",   4},
  { "--ehmmre",  eslARG_REAL,  NULL,   NULL,"x>0",       NULL, "--eent",     NULL, "for --eent: set minimum HMM relative entropy to <x>", 4}, 
/* Customizing null model or priors */
  { "--null",    eslARG_INFILE,  NULL, NULL, NULL,      NULL,      NULL, "--rsearch", "read null (random sequence) model from file <s>", 5 },
  { "--prior",   eslARG_INFILE,  NULL, NULL, NULL,      NULL,      NULL, "--rsearch", "read priors from file <s>", 5 },
/* Building multiple CMs after clustering input MSA */
  { "--ctarget", eslARG_INT,   NULL,   NULL, "n>0" ,    NULL,      NULL,    "--call", "build (at most) <n> CMs by partitioning MSA into <n> clusters", 6 },
  { "--cmaxid",  eslARG_REAL,  NULL,   NULL,"0.<x<1.",  NULL,      NULL,    "--call", "max fractional id b/t 2 clusters is <x>, each cluster -> CM", 6 }, 
  { "--call",    eslARG_NONE,  FALSE,  NULL, NULL,      NULL,      NULL,        NULL, "build a separate CM from every seq in MSA", 6 },
  { "--corig",   eslARG_NONE,  FALSE,  NULL, NULL,      NULL,      NULL,        NULL, "build an additional CM from the original, full MSA", 6 }, 
  { "--cdump",   eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "dump the MSA for each cluster (CM) to file <s>", 6 },
/* Refining the seed alignment */
  { "--refine",  eslARG_OUTFILE, NULL, NULL, NULL,       NULL,   NULL,          NULL, "refine input aln w/Expectation-Maximization, save to <s>", 7 },
  { "--gibbs",   eslARG_NONE,   FALSE, NULL, NULL,       NULL,"--refine",       NULL, "w/--refine, use Gibbs sampling instead of EM", 7 },
  { "-s",        eslARG_INT,     NULL, NULL, "n>0",      NULL,"--gibbs",        NULL, "w/--gibbs, set random number generator seed to <n>",  7 },
  { "-l",        eslARG_NONE,   FALSE, NULL, NULL,       NULL,"--refine",       NULL, "w/--refine, align locally w.r.t the model", 7 },
  { "-a",        eslARG_NONE,   FALSE, NULL, NULL,       NULL,"--refine",       NULL, "print individual sequence scores during MSA refinement", 7 },
  { "--optacc",  eslARG_NONE,"default",NULL,NULL,        NULL,      NULL,       NULL, "align with the Holmes/Durbin optimal accuracy algorithm", 201 },
  { "--cyk",     eslARG_NONE,   FALSE, NULL, NULL,       NULL,"--refine",       NULL, "w/--refine align w/the CYK algorithm, not optimal accuracy", 7 },
  { "--sub",     eslARG_NONE,   FALSE, NULL, NULL,       NULL,"--refine",       NULL, "w/--refine, use sub CM for columns b/t HMM start/end points", 7 },
  { "--nonbanded",eslARG_NONE,  FALSE, NULL, NULL,       NULL,"--refine",       NULL, "do not use bands to accelerate alignment with --refine", 7 },
  { "--tau",     eslARG_REAL,   "1E-7",NULL, "0<x<1",    NULL,"--refine","--nonbanded", "set tail loss prob for --hbanded to <x>", 7 },
  { "--fins",    eslARG_NONE,   FALSE, NULL, NULL,       NULL,"--refine",       NULL, "w/--refine, flush inserts left/right in alignments", 7 },
  { "--mxsize",  eslARG_REAL, "2048.0", NULL, "x>0.",    NULL,"--refine",       NULL, "set maximum allowable DP matrix size to <x> Mb", 7 },
  { "--rdump",   eslARG_OUTFILE, NULL,  NULL, NULL,      NULL,"--refine",       NULL, "w/--refine, print all intermediate alignments to <f>", 7 },

  /* All options below are developer options, only shown if --devhelp invoked */
  /* Developer debugging/experimentation */
  { "--nobalance",eslARG_NONE,  FALSE, NULL, NULL,      NULL,      NULL,        NULL, "don't rebalance the CM; number in strict preorder", 101 },
  { "--regress",  eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "save regression test information to file <s>", 101 },  
  { "--nodetach",eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "do not 'detach' one of two inserts that model same column", 101 },
  { "--elself",  eslARG_REAL,  "0.94", NULL, "0<=x<=1", NULL,      NULL,        NULL, "set EL self transition prob to <x>", 101 },
  { "--eX",      eslARG_REAL,  "6.0",  NULL,"x>0",      NULL,      "--eent", "--ere", "for --eent: set minimum total rel ent param to <x>",  101}, 
  { "--informat",eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL,        NULL, "specify input alignment is in format <s> (Stockholm or Pfam)",  101 },

  /* Developer verbose output options */
  { "--cfile",   eslARG_OUTFILE,  NULL, NULL, NULL,      NULL,      NULL,        NULL, "save count vectors to file <s>", 102 },
  { "--cmtbl",   eslARG_OUTFILE,  NULL, NULL, NULL,      NULL,      NULL,        NULL, "save tabular description of CM topology to file <s>", 102 },
  { "--emap",    eslARG_OUTFILE,  NULL, NULL, NULL,      NULL,      NULL,        NULL, "save consensus emit map to file <s>", 102 },
  { "--gtree",   eslARG_OUTFILE,  NULL, NULL, NULL,      NULL,      NULL,        NULL, "save tree description of master tree to file <s>", 102 },
  { "--gtbl",    eslARG_OUTFILE,  NULL, NULL, NULL,      NULL,      NULL,        NULL, "save tabular description of master tree to file <s>", 102 },
  { "--tfile",   eslARG_OUTFILE,  NULL, NULL, NULL,      NULL,      NULL,        NULL, "dump individual sequence tracebacks to file <s>", 102 },

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
  char         *alifile;	/* name of the alignment file we're building CMs from  */
  int           fmt;		/* format code for alifile */
  ESL_MSAFILE  *afp;            /* open alifile  */
  ESL_ALPHABET *abc;		/* digital alphabet */

  char         *cmfile;         /* file to write CM to                    */
  FILE         *cmfp;           /* CM output file handle                  */

  float        *null;		/* null model                              */
  Prior_t      *pri;		/* mixture Dirichlet prior for the HMM     */

  fullmat_t    *fullmat;        /* if --rsearch, the full RIBOSUM matrix */

  int           be_verbose;	/* standard verbose output, as opposed to one-line-per-CM summary */
  int           nali;		/* which # alignment this is in file */
  int           ncm_total;      /* which # CM this is that we're constructing (we may build > 1 per file) */
  int           namewidth;      /* max length of a CM name, nec for pretty tabular formatting */
  ESL_RANDOMNESS *r;            /* source of randomness, only created if --gibbs enabled */

  /* optional output files */
  FILE         *cfp;            /* for --cfile */
  FILE         *tblfp;          /* for --cmtbl */
  FILE         *efp;            /* for --emap */
  FILE         *gfp;            /* for --gtree */
  FILE         *gtblfp;         /* for --gtbl */
  FILE         *tracefp;        /* for --tfile */
  FILE         *cdfp;           /* if --cdump, output file handle for dumping clustered MSAs */
  FILE         *refinefp;       /* if --refine, output file handle for dumping refined MSAs */
  FILE         *rdfp;           /* if --rfile, output file handle for dumping intermediate MSAs during iterative refinement */

  ComLog_t      *comlog;       /* the comlog, same for all CMs, cfg.comlog serves as template 
				 * for all CMs, and is copied to each CM data structure */
  int           argc;          /* used to create the comlog */
  char        **argv;          /* used to create the comlog, be careful not to free this though, it's just a ptr */
};

static char usage[]  = "[-options] <cmfile output> <alignment file>";
static char banner[] = "build RNA covariance model(s) from alignment";

static int    init_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static void   master (const ESL_GETOPTS *go, struct cfg_s *cfg);
static int    process_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t **ret_cm, Parsetree_t **ret_mtr, Parsetree_t ***ret_msa_tr);
static int    output_result(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, int msaidx, int cmidx, ESL_MSA *msa, CM_t *cm, Parsetree_t *mtr, Parsetree_t **tr);
static int    check_and_clean_msa(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa);
static int    set_relative_weights(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa);
static int    build_model(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t **ret_cm, Parsetree_t **ret_mtr, Parsetree_t ***ret_msa_tr);
static int    set_model_name(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t *cm);
static int    set_model_cutoffs(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t *cm);
static int    set_effective_seqnumber(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t *cm, const Prior_t *pri);
static int    parameterize(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, const Prior_t *prior);
static int    name_msa(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa, int nali);
static double default_target_relent(const ESL_ALPHABET *abc, int M, double eX);
static void   strip_wuss(char *ss);
static int    refine_msa(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *orig_cm, ESL_MSA *input_msa, Parsetree_t **input_msa_tr, CM_t **ret_cm, ESL_MSA **ret_msa, Parsetree_t **ret_mtr, Parsetree_t ***ret_tr, int *ret_niter);
static int    get_unaln_seqs_from_msa(const ESL_MSA *msa, ESL_SQ ***ret_sq);
static int    convert_parsetrees_to_unaln_coords(Parsetree_t **tr, ESL_MSA *msa);
static int    initialize_cm(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm);
/* static void   model_trace_info_dump(FILE *ofp, CM_t *cm, Parsetree_t *tr, char *aseq); */
/* functions for dividing input MSA into clusters */
static int    select_node(ESL_TREE *T, double *diff, double mindiff, int **ret_clust, int *ret_nc, int *ret_best, char *errbuf);
static float  find_mindiff(ESL_TREE *T, double *diff, int target_nc, int **ret_clust, int *ret_nc, float *ret_mindiff, char *errbuf);
static int    MSADivide(ESL_MSA *mmsa, int do_all, int do_mindiff, int do_nc, float mindiff, int target_nc, int do_orig, int *ret_num_msa, ESL_MSA ***ret_cmsa, char *errbuf);
static int    write_cmbuild_info_to_comlog(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int    flatten_insert_emissions(CM_t *cm);
static int    print_run_info(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf);
static int    print_column_headings(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf);
static void   print_refine_column_headings(const ESL_GETOPTS *go, const struct cfg_s *cfg);
static int    print_countvectors(const struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int    get_namewidth(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go = NULL;   /* command line processing                     */
  ESL_STOPWATCH   *w  = esl_stopwatch_Create();
  if(w == NULL) cm_Fail("Memory allocation error, stopwatch could not be created.");
  esl_stopwatch_Start(w);
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
  if (esl_opt_GetBoolean(go, "--devhelp") == TRUE) 
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
      puts("\ncustomization of null model and priors:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80);
      puts("\noptions for building multiple CMs after clustering input MSA:");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80);
      puts("\nexpert options for refining the input alignment:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80);
      puts("\nundocumented developer options for debugging, experimentation:");
      esl_opt_DisplayHelp(stdout, go, 101, 2, 80);
      puts("\nundocumented developer options for verbose output/debugging:");
      esl_opt_DisplayHelp(stdout, go, 102, 2, 80);
      exit(0);
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
      puts("\ncustomization of null model and priors:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80);
      puts("\noptions for building multiple CMs after clustering input MSA:");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80);
      puts("\nexpert options for refining the input alignment:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80);
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
  cfg.cmfile     = esl_opt_GetArg(go, 1); 
  cfg.alifile    = esl_opt_GetArg(go, 2);
  cfg.afp        = NULL;	           /* created in init_cfg() */
  cfg.abc        = NULL;	           /* created in init_cfg() */
  cfg.cmfp       = NULL;	           /* opened in init_cfg() */
  cfg.null       = NULL;	           /* created in init_cfg() */
  cfg.pri        = NULL;                   /* created in init_cfg() */
  cfg.fullmat    = NULL;                   /* read (possibly) in init_cfg() */
  cfg.r          = NULL;	           /* created (possibly) in init_cfg() */
  cfg.comlog     = NULL;	           /* created in init_cfg() */
  cfg.namewidth  = 0;
  /* optional output files, opened in init_cfg(), if at all */
  cfg.cfp        = NULL;
  cfg.tblfp      = NULL;
  cfg.efp        = NULL;
  cfg.gfp        = NULL;
  cfg.gtblfp     = NULL;
  cfg.tracefp    = NULL;
  cfg.cdfp       = NULL;
  cfg.refinefp   = NULL;
  cfg.rdfp       = NULL;

  /* print the banner */
  cm_banner(stdout, argv[0], banner);

  if   (esl_opt_IsDefault(go, "--informat")) cfg.fmt = eslMSAFILE_UNKNOWN; /* autodetect sequence file format by default. */ 
  else { 
    cfg.fmt = esl_sqio_FormatCode(esl_opt_GetString(go, "--informat"));
    if(cfg.fmt == eslSQFILE_UNKNOWN)                                  cm_Fail("Can't recognize sequence file format: %s. valid options are: stockholm or pfam\n", esl_opt_GetString(go, "--informat"));
    if(cfg.fmt != eslMSAFILE_STOCKHOLM && cfg.fmt != eslMSAFILE_PFAM) cm_Fail("Sequence file format: %s is not accepted by cmbuild, valid options are: stockholm or pfam\n", esl_opt_GetString(go, "--informat"));
  }

  cfg.be_verbose = esl_opt_GetBoolean(go, "-v");
  cfg.nali       = 0;		           

  /* check if cmfile already exists, if it does and -F was not enabled then die */
  if (((! esl_opt_GetBoolean(go, "-F")) && (! esl_opt_GetBoolean(go, "-A"))) && esl_FileExists(cfg.cmfile))
    cm_Fail("CM file %s already exists. Either use -F to overwrite it, rename it, or delete it.", cfg.cmfile); 

  /* do work */
  master(go, &cfg);

  /* Clean up the cfg. 
   */
  /* close all output files */
  if (cfg.cfp != NULL) {
    printf("# Count vectors saved in file %s.\n", esl_opt_GetString(go, "--cfile"));
    fclose(cfg.cfp); 
  }
  if (cfg.tblfp != NULL) {
    printf("# CM topology description saved in file %s.\n", esl_opt_GetString(go, "--cmtbl"));
    fclose(cfg.tblfp); 
  }
  if (cfg.efp != NULL) {
    printf("# CM emit map saved in file %s.\n", esl_opt_GetString(go, "--emap"));
    fclose(cfg.efp); 
  }
  if (cfg.gfp != NULL) {
    printf("# Guide tree description saved in file %s.\n", esl_opt_GetString(go, "--gtree"));
    fclose(cfg.gfp); 
  }
  if (cfg.gtblfp != NULL) {
    printf("# Guide tree tabular description saved in file %s.\n", esl_opt_GetString(go, "--gtbl"));
    fclose(cfg.gtblfp); 
  }
  if (cfg.tracefp != NULL) {
    printf("# Implicit parsetrees of seqs from input alignment saved in file %s.\n", esl_opt_GetString(go, "--tfile"));
    fclose(cfg.tracefp); 
  }
  if (cfg.cdfp != NULL) {
    printf("# Alignments for each cluster saved in file %s.\n", esl_opt_GetString(go, "--cdump"));
    fclose(cfg.cdfp); 
  }
  if (cfg.refinefp != NULL) {
    printf("# Refined alignments used to build CMs saved in file %s.\n", esl_opt_GetString(go, "--refine"));
    fclose(cfg.refinefp); 
  }
  if (cfg.rdfp != NULL) {
    printf("# Intermediate alignments from MSA refinement saved in file %s.\n", esl_opt_GetString(go, "--rdump"));
    fclose(cfg.rdfp); 
  }
  if (cfg.afp   != NULL) esl_msafile_Close(cfg.afp);
  if (cfg.abc   != NULL) esl_alphabet_Destroy(cfg.abc);
  if (cfg.cmfp  != NULL) fclose(cfg.cmfp);
  if (cfg.pri   != NULL) Prior_Destroy(cfg.pri);
  if (cfg.null  != NULL) free(cfg.null);
  if (cfg.r     != NULL) esl_randomness_Destroy(cfg.r);
  if (cfg.comlog!= NULL) FreeComLog(cfg.comlog);

  esl_getopts_Destroy(go);
  esl_stopwatch_Stop(w);
  printf("#\n");
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
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
 *    cfg->comlog  - only allocated
 */
static int
init_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  int status;

  /* open input alignment file */
  status = esl_msafile_Open(cfg->alifile, cfg->fmt, NULL, &(cfg->afp));
  if      (status == eslENOTFOUND) ESL_FAIL(status, errbuf, "Alignment file %s doesn't exist or is not readable\n", cfg->alifile);
  else if (status == eslEFORMAT) ESL_FAIL(status, errbuf, "Couldn't determine format of alignment %s\n", cfg->alifile);
  else if (status != eslOK)      ESL_FAIL(status, errbuf, "Alignment file open failed with error %d\n", status);
  cfg->fmt = cfg->afp->format;

  /* Set the msafile alphabet as RNA, if it's DNA we're fine. 
   * If it's not RNA nor DNA, we can't deal with it anyway,
   * so we're hardcoded to RNA.
   */
  cfg->abc = esl_alphabet_Create(eslRNA);
  if(cfg->abc == NULL) ESL_FAIL(status, errbuf, "Failed to create alphabet for sequence file");
  esl_msafile_SetDigital(cfg->afp, cfg->abc);

  /* open CM file for writing */
  if (esl_opt_GetBoolean(go, "-A")) { /* we're appending to a CM file */
    if ((cfg->cmfp = fopen(cfg->cmfile, "a")) == NULL) ESL_FAIL(status, errbuf, "Failed to open CM file %s to append to", cfg->cmfile);
  }
  else { /* we're starting a new CM file */
    if ((cfg->cmfp = fopen(cfg->cmfile, "w")) == NULL) ESL_FAIL(status, errbuf, "Failed to open CM file %s for writing", cfg->cmfile);
  }
  /* Set up the prior */
  if (esl_opt_GetString(go, "--prior") != NULL)
    {
      FILE *pfp;
      if ((pfp = fopen(esl_opt_GetString(go, "--prior"), "r")) == NULL)
	cm_Fail("Failed to open prior file %s\n", esl_opt_GetString(go, "--prior"));
      if ((cfg->pri = Prior_Read(pfp)) == NULL)
	cm_Fail("Failed to parse prior file %s\n", esl_opt_GetString(go, "--prior"));
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

  /* if --corig enabled, make sure either --cmaxid, --ctarget, or --call also enabled */
  if (esl_opt_GetBoolean(go, "--corig"))
    if((esl_opt_IsDefault(go, "--ctarget")) && (esl_opt_IsDefault(go, "--cmaxid")) && (esl_opt_IsDefault(go, "--call")))
      cm_Fail("--corig only makes sense in combination with --ctarget, --cmaxid, OR --call");

  /* if --gibbs enabled, open output file for refined MSAs, and seed RNG */
  if(esl_opt_GetBoolean(go, "--gibbs"))
    {
      /* create RNG */
      if (! esl_opt_IsDefault(go, "-s")) 
	cfg->r = esl_randomness_Create((long) esl_opt_GetInteger(go, "-s"));
      else cfg->r = esl_randomness_CreateTimeseeded();
      if (cfg->r == NULL) ESL_FAIL(eslEINVAL, errbuf, "Failed to create random number generator: probably out of memory");
    }

  /* open output files */
  /* optionally, open count vector file */
  if (esl_opt_GetString(go, "--cfile") != NULL) {
    if ((cfg->cfp = fopen(esl_opt_GetString(go, "--cfile"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --cfile output file %s\n", esl_opt_GetString(go, "--cfile"));
    }
  /* optionally, open CM tabular file */
  if (esl_opt_GetString(go, "--cmtbl") != NULL) {
    if ((cfg->tblfp = fopen(esl_opt_GetString(go, "--cmtbl"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --cmtbl output file %s\n", esl_opt_GetString(go, "--cmtbl"));
    }
  /* optionally, open emit map file */
  if (esl_opt_GetString(go, "--emap") != NULL) {
    if ((cfg->efp = fopen(esl_opt_GetString(go, "--emap"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --emap output file %s\n", esl_opt_GetString(go, "--emap"));
    }
  /* optionally, open guide tree file */
  if (esl_opt_GetString(go, "--gtree") != NULL) {
    if ((cfg->gfp = fopen(esl_opt_GetString(go, "--gtree"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --gtree output file %s\n", esl_opt_GetString(go, "--gtree"));
    }
  /* optionally, open master tree file */
  if (esl_opt_GetString(go, "--gtbl") != NULL) {
    if ((cfg->gtblfp = fopen(esl_opt_GetString(go, "--gtbl"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --gtbl output file %s\n", esl_opt_GetString(go, "--gtbl"));
    }
  /* optionally, open trace file */
  if (esl_opt_GetString(go, "--tfile") != NULL) {
    if ((cfg->tracefp = fopen(esl_opt_GetString(go, "--tfile"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --tfile output file %s\n", esl_opt_GetString(go, "--tfile"));
    }
  /* if --refine enabled, open output file for refined MSAs */
  if (esl_opt_GetString(go, "--refine") != NULL)
    {
      if ((cfg->refinefp = fopen(esl_opt_GetString(go, "--refine"), "w")) == NULL)
	ESL_FAIL(eslFAIL, errbuf, "Failed to open output file %s for writing MSAs from --refine to", esl_opt_GetString(go, "--refine"));
    }
  /* optionally, open --rdump output alignment file */
  if (esl_opt_GetString(go, "--rdump") != NULL) {
    if ((cfg->rdfp = fopen(esl_opt_GetString(go, "--rdump"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --rdump output file %s\n", esl_opt_GetString(go, "--rdump"));
    }
  /* if --cdump enabled, open output file for cluster MSAs */
  if (esl_opt_GetString(go, "--cdump") != NULL)
    {
      /* check to make sure there's a reason for this option, --cmaxid, --ctarget or --call MUST also be enabled */
      if((esl_opt_IsDefault(go, "--ctarget")) && (esl_opt_IsDefault(go, "--cmaxid")) && (esl_opt_IsDefault(go, "--call")))
	cm_Fail("--cdump only makes sense in combination with --ctarget, --cmaxid, OR --call");
      if ((cfg->cdfp = fopen(esl_opt_GetString(go, "--cdump"), "w")) == NULL)
	cm_Fail("Failed to open output file %s for writing MSAs to", esl_opt_GetString(go, "--cdump"));
    }
  /* create the comlog */
  cfg->comlog = CreateComLog();
  if((status = write_cmbuild_info_to_comlog(go, cfg, errbuf)) != eslOK) return status;

  if (cfg->pri    == NULL) ESL_FAIL(eslEINVAL, errbuf, "alphabet initialization failed");
  if (cfg->null   == NULL) ESL_FAIL(eslEINVAL, errbuf, "null model initialization failed");
  if (cfg->comlog == NULL) ESL_FAIL(eslEINVAL, errbuf, "comlog initialization failed");

  cfg->nali = 0;
  cfg->ncm_total = 0;
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
  char     errbuf[cmERRBUFSIZE];
  ESL_MSA *msa = NULL;
  CM_t    *cm = NULL;
  Parsetree_t  *mtr;
  Parsetree_t **tr;
  int          i = 0;
  int      niter = 0;
  /* new_* data structures, created in refine_msa() if --refine enabled */
  CM_t         *new_cm;  
  Parsetree_t  *new_mtr;
  Parsetree_t **new_tr;
  ESL_MSA      *new_msa;
  /* cluster option related variables */
  int          do_cluster; /* TRUE if --ctarget || --cmaxid || --call */
  int          do_ctarget; /* TRUE if --ctarget */
  int          do_cmindiff; /* TRUE if --cmaxid  */
  int          do_call;    /* TRUE if --call */
  int          nc;         /* number of clusters, only != 0 if do_ctarget */
  float        mindiff;    /* minimum fractional diff b/t clusters, only != 0. if do_cmindiff */
  int          ncm = 1;    /* number of CMs to be built for current MSA */
  int          c   = 0;    /* counter over CMs built for a single MSA */
  ESL_MSA    **cmsa;       /* pointer to cluster MSAs to build CMs from */

  if ((status = init_cfg(go, cfg, errbuf))         != eslOK) cm_Fail(errbuf);
  if ((status = print_run_info (go, cfg, errbuf))  != eslOK) cm_Fail(errbuf);

  cfg->nali = 0;
  cfg->ncm_total = 0;

  do_ctarget  = (esl_opt_IsDefault(go, "--ctarget")) ? FALSE : TRUE;
  do_cmindiff = (esl_opt_IsDefault(go, "--cmaxid"))  ? FALSE : TRUE;
  do_call     = esl_opt_GetBoolean(go, "--call");
  do_cluster = (do_ctarget || do_cmindiff || do_call) ? TRUE : FALSE;
  if((do_ctarget + do_cmindiff + do_call) > TRUE) cm_Fail("More than one of --ctarget, --cmaxid, --call were enabled, shouldn't happen.");

  nc      = do_ctarget  ? esl_opt_GetInteger(go, "--ctarget")    : 0;
  mindiff = do_cmindiff ? (1. - esl_opt_GetReal(go, "--cmaxid")) : 0.;

  /* predict maximum length of CM name for pretty formatting */
  if((status = get_namewidth(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);

  while ((status = esl_msa_Read(cfg->afp, &msa)) != eslEOF)
    {
      if      (status == eslEFORMAT)  cm_Fail("Alignment file parse error, line %d of file %s:\n%s\nOffending line is:\n%s\n", cfg->afp->linenumber, cfg->afp->fname, cfg->afp->errbuf, cfg->afp->buf);
      else if (status != eslOK)       cm_Fail("Alignment file read unexpectedly failed with code %d\n", status);
      cfg->nali++;  

      /* if it's unnamed, name the MSA, we require a name (different from 
       * HMMER 3), because it will be used to name the CM. */
      if(name_msa(go, errbuf, msa, cfg->nali) != eslOK) cm_Fail(errbuf);
      if(msa->name == NULL)                             cm_Fail("Error naming MSA");
      ncm = 1;     /* default: only build 1 CM for each MSA in alignment file */

      if(do_cluster) /* divide input MSA into clusters, and build CM from each cluster */
	{
	  if((status = MSADivide(msa, do_call, do_cmindiff, do_ctarget, mindiff, nc,
				 esl_opt_GetBoolean(go, "--corig"), &ncm, &cmsa, errbuf)) != eslOK) cm_Fail(errbuf);
	  esl_msa_Destroy(msa); /* we've copied the master msa into cmsa[ncm], we can delete this copy */
	}
      for(c = 0; c < ncm; c++)
	{
	  cfg->ncm_total++;  
	  if(do_cluster) {
	      msa = cmsa[c];
	      if(esl_opt_GetString(go, "--cdump") != NULL) esl_msa_Write(cfg->cdfp, msa, cfg->fmt); 
	  }

	  /* if being verbose, print some stuff about what we're about to do.
	   */
	  if (cfg->be_verbose) {
	    fprintf(stdout, "Alignment:           %s\n",           msa->name);
	    fprintf(stdout, "Number of sequences: %d\n",           msa->nseq);
	    fprintf(stdout, "Number of columns:   %" PRId64 "\n",  msa->alen);
	    if(esl_opt_GetString(go, "--rsearch") != NULL)
	      printf ("RIBOSUM Matrix:      %s\n",  cfg->fullmat->name);
	    fputs("", stdout);
	    fflush(stdout);
	  }

	  /* msa -> cm */
	  if ((status = process_workunit(go, cfg, errbuf, msa, &cm, &mtr, &tr)) != eslOK) cm_Fail(errbuf);
	  /* optionally, iterate over cm -> parsetrees -> msa -> cm ... until convergence, via EM or Gibbs */
	  if (! esl_opt_IsDefault(go, "--refine")) {
	    fprintf(stdout, "#\n");
	    fprintf(stdout, "# Refining MSA for CM: %s (aln: %4d cm: %6d)\n", cm->name, cfg->nali, cfg->ncm_total);
	    if ((status = refine_msa(go, cfg, errbuf, cm, msa, tr, &new_cm, &new_msa, &new_mtr, &new_tr, &niter)) != eslOK) cm_Fail(errbuf);
	    if (niter > 1) { /* if niter == 1, we didn't make a new CM (new_cm == cm) mtr, or tr, so we don't free them */
	      FreeCM(cm); 
	      cm = new_cm; 
	      for(i = 0; i < msa->nseq; i++) FreeParsetree(tr[i]);
	      tr = new_tr;
	      FreeParsetree(mtr);
	      mtr = new_mtr;
	      esl_msa_Destroy(msa);
	      msa = new_msa;
	    } 
	  }	  
	  /* output cm */
	  if ((status = output_result(go, cfg, errbuf, cfg->nali, cfg->ncm_total, msa,  cm, mtr, tr)) != eslOK) cm_Fail(errbuf);
	  
	  if(cfg->be_verbose) { 
	    fprintf(stdout, "\n");
	    SummarizeCM(stdout, cm);  
	    fprintf(stdout, "//\n");
	  }

	  FreeCM(cm);
	  fflush(cfg->cmfp);

	  if(tr != NULL) {
	    for(i = 0; i < msa->nseq; i++) FreeParsetree(tr[i]);
	    free(tr);
	  }
	  if(mtr != NULL) FreeParsetree(mtr);

	  esl_msa_Destroy(msa);
	}
    }
  if(do_cluster) free(cmsa);
  if(cfg->fullmat != NULL) FreeMat(cfg->fullmat);
  return;
}

/* A work unit consists of one multiple alignment, <msa>.
 * The job is to turn it into a new CM, returned in <*ret_cm>.
 * 
 */
static int
process_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t **ret_cm, Parsetree_t **ret_mtr, Parsetree_t ***ret_msa_tr)
{
  CM_t *cm = NULL;
  int status;

  if ((status =  check_and_clean_msa    (go, cfg, errbuf, msa))                           != eslOK) goto ERROR;
  if ((status =  set_relative_weights   (go, cfg, errbuf, msa))                           != eslOK) goto ERROR;
  if ((status =  build_model            (go, cfg, errbuf, msa, &cm, ret_mtr, ret_msa_tr)) != eslOK) goto ERROR;
  if ((status =  set_model_name         (go, cfg, errbuf, msa, cm))                       != eslOK) goto ERROR;
  if ((status =  set_model_cutoffs      (go, cfg, errbuf, msa, cm))                       != eslOK) goto ERROR;
  if ((status =  set_effective_seqnumber(go, cfg, errbuf, msa, cm, cfg->pri))             != eslOK) goto ERROR;
  if ((status =  parameterize           (go, cfg, errbuf, cm, cfg->pri))                  != eslOK) goto ERROR;
  
  *ret_cm = cm;
  return eslOK;

 ERROR:
  if(cm != NULL) FreeCM(cm);
  *ret_cm = NULL;
  return status;
}

/* refine_msa() 
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
	   Parsetree_t **input_msa_tr, CM_t **ret_cm, ESL_MSA **ret_msa, Parsetree_t **ret_mtr, Parsetree_t ***ret_tr, 
	   int *ret_niter)
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
  Parsetree_t     *mtr         = NULL;
  Parsetree_t    **tr          = NULL;

  /* check contract */
  if(input_msa       == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "refine_msa(), input_msa passed in as NULL");
  if(input_msa->name == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "refine_msa(), input_msa must have a name");
  if(init_cm         == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "refine_msa(), init_cm passed in as NULL");
  if(ret_cm          == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "refine_msa(), ret_cm is NULL");
  if(ret_mtr         == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "refine_msa(), ret_mtr is NULL");
  if(ret_tr          == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "refine_msa(), ret_tr is NULL");

  /* copy input MSA's name, we'll copy it to the MSA we create at each iteration */
  if((status = esl_strdup(input_msa->name, -1, &(msa_name))) != eslOK) ESL_FAIL(eslEINCOMPAT, errbuf, "Memory allocation error.");

  ESL_ALLOC(sc, sizeof(float) * nseq);
  esl_vec_FSet(sc, nseq, 0.);

  get_unaln_seqs_from_msa(input_msa, &sq); /* we need sqs for Parsetrees2Alignment */
  seqs_to_aln = CreateSeqsToAlnFromSq(sq, nseq, FALSE);

  /* determine scores of implicit parsetrees of input MSA seqs to initial CM */
  convert_parsetrees_to_unaln_coords(input_msa_tr, input_msa);
  for(i = 0; i < nseq; i++) { 
    if((status = ParsetreeScore(init_cm, errbuf, input_msa_tr[i], sq[i]->dsq, FALSE, &(sc[i]), NULL)) != eslOK) return status;
  }
  oldscore = esl_vec_FSum(sc, nseq);

  /* print header for tabular output */
  print_refine_column_headings(go, cfg);
  fprintf(stdout, "  %5d %13.2f %10s\n", iter, oldscore, "-");

  /* print initial alignment to --rdump file, if --rdump was enabled */
  if(cfg->rdfp != NULL) 
    if((status = esl_msa_Write(cfg->rdfp, input_msa, cfg->fmt)) != eslOK) ESL_FAIL(status, errbuf, "refine_msa(), esl_msa_Write() call failed.");
  
  while(1)
    {
      iter++;
      if(iter == 1) { cm = init_cm; msa = input_msa; }
      
      /* 1. cm -> parsetrees */
      if(iter > 1) FreePartialSeqsToAln(seqs_to_aln, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE);
                                                  /* sq,    tr, cp9_tr, post, sc,   pp,   struct_sc */ 
      /* initialize/configure CM, we may be doing HMM banded alignment for ex. */
      initialize_cm(go, cfg, errbuf, cm);
      if((status = DispatchAlignments(cm, errbuf, seqs_to_aln, NULL, NULL, 0, 0, 0, (! esl_opt_GetBoolean(go, "-a")), TRUE, cfg->r, 
				      esl_opt_GetReal(go, "--mxsize"), stdout)) != eslOK) return status;
      
      /* sum parse scores and check for convergence */
      totscore = esl_vec_FSum(seqs_to_aln->sc, nseq);
      delta    = (totscore - oldscore) / fabs(totscore);
      if(esl_opt_GetBoolean(go, "-a")) print_refine_column_headings(go, cfg);
      fprintf(stdout, "  %5d %13.2f %10.3f\n", iter, totscore, delta);
      if(delta <= threshold && delta >= 0) break; /* only way out of while(1) loop */
      oldscore = totscore;

      /* 2. parsetrees -> msa */
      if( iter > 1) esl_msa_Destroy(msa);
      msa = NULL; /* even if iter == 1; we set msa to NULL, so we don't klobber input_msa */
      if((status = Parsetrees2Alignment(cm, errbuf, cm->abc, seqs_to_aln->sq, NULL, seqs_to_aln->tr, NULL, NULL, nseq, FALSE, FALSE, &msa)) != eslOK) 
	ESL_FAIL(status, errbuf, "refine_msa(), Parsetrees2Alignment() call failed.");
      if((status = esl_strdup(msa_name, -1, &(msa->name))) != eslOK) ESL_FAIL(status, errbuf, "refine_msa(), esl_strdup() call failed.");
      esl_msa_Digitize(msa->abc, msa);
      
      /* print intermediate alignment to --rdump file, if --rdump was enabled */
      if(cfg->rdfp != NULL) 
	if((status = esl_msa_Write(cfg->rdfp, msa, cfg->fmt)) != eslOK) ESL_FAIL(status, errbuf, "refine_msa(), esl_msa_Write() call failed.");

      /* 3. msa -> cm */
      if(iter > 1) { /* free previous iterations cm, mtr and tr */
	FreeCM(cm);
	FreeParsetree(mtr);
	for(i = 0; i < nseq; i++) FreeParsetree(tr[i]);
	free(tr);
      }
      cm = NULL; /* even if iter == 1; we set cm to NULL, so we don't klobber init_cm */
      mtr= NULL;
      tr = NULL;
      if ((status = process_workunit(go, cfg, errbuf, msa, &cm, &mtr, &tr))  != eslOK) cm_Fail(errbuf);
    }

  /* write out final alignment to --refine output file */
  if((status = esl_msa_Write(cfg->refinefp, msa, cfg->fmt)) != eslOK) ESL_FAIL(status, errbuf, "refine_msa(), esl_msa_Write() call failed.");

  /* if CM was in local mode for aligning input MSA seqs, make it global so we can write it out */
  if((cm->flags & CMH_LOCAL_BEGIN) || (cm->flags & CMH_LOCAL_END)) ConfigGlobal(cm);

  *ret_cm  = cm;
  *ret_msa = msa;
  *ret_tr  = tr;
  *ret_mtr = mtr;
  *ret_niter = iter;

  /* clean up */
  FreeSeqsToAln(seqs_to_aln);
  free(sc);
  free(msa_name);

  return eslOK;

 ERROR:
  /* no cleanup, we die */
  cm_Fail("in refine_msa(), error, status: %d\n", status);
  return status; /* NEVERREACHED */
}


/* Function: print_column_headings()
 * Date:     EPN, Fri Feb 29 10:08:25 2008
 *
 * Purpose:  Print column headings for tabular output to output file (stdout unless -o). 
 *
 * Returns:  eslOK on success
 */
static int
print_column_headings(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf)
{
  int status;
  char *namedashes;
  int ni;
  ESL_ALLOC(namedashes, sizeof(char) * cfg->namewidth+1);
  namedashes[cfg->namewidth] = '\0';
  for(ni = 0; ni < cfg->namewidth; ni++) namedashes[ni] = '-';

  fprintf(stdout, "# %-4s  %-6s  %-*s  %8s  %8s  %6s  %5s  %4s  %4s  %12s\n",    "",     "", cfg->namewidth, "",                     "",         "",         "",     "",      "", "", "rel entropy");
  fprintf(stdout, "# %-4s  %-6s  %-*s  %8s  %8s  %6s  %5s  %4s  %4s  %12s\n",    "",     "", cfg->namewidth, "",                     "",         "",         "",     "",      "", "", "------------");
  fprintf(stdout, "# %4s  %-6s  %-*s  %8s  %8s  %6s  %5s  %4s  %4s  %5s  %5s\n",  "aln",  "cm idx", cfg->namewidth, "name",                 "nseq",     "eff_nseq", "alen",   "clen",  "bps", "bifs",  "CM",     "HMM");
  fprintf(stdout, "# %-4s  %-6s  %-*s  %8s  %8s  %6s  %5s  %4s  %4s  %5s  %5s\n", "----", "------", cfg->namewidth, namedashes,             "--------", "--------", "------", "-----", "----", "----", "-----", "-----");

  free(namedashes);
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "Memory allocation error in print_column_headings()");
  return status; /* NEVERREACHED */
}

static int
output_result(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, int msaidx, int cmidx, ESL_MSA *msa, CM_t *cm, Parsetree_t *mtr, Parsetree_t **tr)
{
  int status;
  int i;
  float sc, struct_sc;

  if(msaidx == 1 && cmidx == 1) { 
    if((status = print_column_headings(go, cfg, errbuf)) != eslOK) return status;
  }

  /* copy the cmbuild command info to the CM */
  if ((status = CopyComLog(cfg->comlog, cm->comlog)) != eslOK) ESL_FAIL(eslFAIL, errbuf, "Problem copying com log info to CM. Probably out of memory.");
  if ((status = cm_Validate(cm, 0.0001, errbuf))     != eslOK) return status;
  if ((status = CMFileWrite(cfg->cmfp, cm, esl_opt_GetBoolean(go, "--binary"), errbuf)) != eslOK) return status;
  /* build the HMM, so we can print the CP9 relative entropy */
  if(!(build_cp9_hmm(cm, &(cm->cp9), &(cm->cp9map), FALSE, 0.0001, 0))) ESL_FAIL(eslFAIL, errbuf, "Couldn't build a CP9 HMM from the CM.");

  fprintf(stdout, "%6d  %6d  %-*s  %8d  %8.2f  %6" PRId64 "  %5d  %4d  %4d  %5.3f  %5.3f\n",
	  msaidx,
	  cmidx,
	  cfg->namewidth,
	  cm->name, 
	  msa->nseq,
	  cm->eff_nseq,
	  msa->alen,
	  cm->clen, 
	  CMCountStatetype(cm, MP_st), 
	  CMCountStatetype(cm, B_st), 
	  cm_MeanMatchRelativeEntropy(cm),
	  cp9_MeanMatchRelativeEntropy(cm));


  /* dump optional info to files: */
  if(cfg->tblfp != NULL) PrintCM(cfg->tblfp, cm); /* tabular description of CM topology */
  /* emit map */
  if(cfg->efp != NULL) {
    CMEmitMap_t *emap;
    emap = CreateEmitMap(cm);
    DumpEmitMap(cfg->efp, emap, cm);
    FreeEmitMap(emap);
  }
  /* save tabular description of guide tree topology, if nec */
  if(cfg->gtblfp != NULL) PrintParsetree(cfg->gtblfp, mtr);  
  /* save tree description of guide tree topology, if nec */
  if(cfg->gfp    != NULL) MasterTraceDisplay(cfg->gfp, mtr, cm);

  /* save parsetrees if nec */
  if(cfg->tracefp != NULL) { 
    for (i = 0; i < msa->nseq; i++) { 
      fprintf(cfg->tracefp, "> %s\n", msa->sqname[i]);

      if((status = ParsetreeScore(cm, errbuf, tr[i], msa->ax[i], FALSE, &sc, &struct_sc)) != eslOK) return status;
      fprintf(cfg->tracefp, "  %16s %.2f bits\n", "SCORE:", sc);
      fprintf(cfg->tracefp, "  %16s %.2f bits\n", "STRUCTURE SCORE:", struct_sc);
      ParsetreeDump(cfg->tracefp, tr[i], cm, msa->ax[i], NULL, NULL); /* NULLs are dmin, dmax */
      fprintf(cfg->tracefp, "//\n");
    }
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
  if (cfg->be_verbose) {
    fprintf(stdout, "%-40s ... ", "Checking MSA");  
    fflush(stdout); 
  }

  if (esl_opt_GetBoolean(go, "--rf") && msa->rf == NULL)        ESL_FAIL(eslFAIL, errbuf, "Alignment has no reference coord annotation.\n");
  if (msa->ss_cons == NULL)                                     ESL_FAIL(eslFAIL, errbuf, "Alignment did not contain consensus structure annotation.\n");
  if (! clean_cs(msa->ss_cons, msa->alen, (! cfg->be_verbose))) ESL_FAIL(eslFAIL, errbuf, "Failed to parse consensus structure annotation\n");
  if (esl_opt_GetBoolean(go, "--ignorant"))                     strip_wuss(msa->ss_cons); /* --ignorant, remove all bp info */

  if (! esl_opt_IsDefault(go, "--rsearch")) { 
    if(msa->nseq != 1) ESL_FAIL(eslEINCOMPAT, errbuf,"with --rsearch option, all of the input alignments must have exactly 1 sequence");
    /* We can't have ambiguous bases in the MSA, only A,C,G,U will do. The reason is that rsearch_CMProbifyEmissions() expects each
     * cm->e prob vector to have exactly 1.0 count for exactly 1 singlet or base pair. If we have ambiguous residues we'll have a 
     * fraction of a count for more than one residue/base pair for some v. 
     * ribosum_MSA_resolve_degeneracies() replaces ambiguous bases with most likely compatible base */
    ribosum_MSA_resolve_degeneracies(cfg->fullmat, msa); /* cm_Fails() if some error is encountered */
  }

  /* MSA better have a name, we named it before */
  if(msa->name == NULL) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "MSA is nameless, (we thought we named it...) shouldn't happen");

  if (cfg->be_verbose) fprintf(stdout, "done.\n");
  return eslOK;
}

/* set_relative_weights():
 * Set msa->wgt vector, using user's choice of relative weighting algorithm.
 */
static int
set_relative_weights(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa)
{
  if (cfg->be_verbose) {
    fprintf(stdout, "%-40s ... ", "Relative sequence weighting");  
    fflush(stdout); 
  }

  if      (esl_opt_GetBoolean(go, "--wnone"))                  esl_vec_DSet(msa->wgt, msa->nseq, 1.);
  else if (esl_opt_GetBoolean(go, "--wgiven"))                 ;
  else if (msa->nseq >= esl_opt_GetInteger(go, "--pbswitch"))  esl_msaweight_PB(msa);
  else if (esl_opt_GetBoolean(go, "--wpb"))                    esl_msaweight_PB(msa);
  else if (esl_opt_GetBoolean(go, "--wgsc"))                   esl_msaweight_GSC(msa);
  else if (esl_opt_GetBoolean(go, "--wblosum"))                esl_msaweight_BLOSUM(msa, esl_opt_GetReal(go, "--wid"));

  if (cfg->be_verbose) fprintf(stdout, "done.\n");
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
build_model(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t **ret_cm, Parsetree_t **ret_mtr, Parsetree_t ***ret_msa_tr)
{
  int status;
  Parsetree_t     **tr;
  Parsetree_t     *mtr;
  int idx;
  CM_t *cm;
  char *aseq;                   

  if (cfg->be_verbose) {
    fprintf(stdout, "%-40s ... ", "Constructing model architecture"); 
    fflush(stdout);
  }

  if((status = HandModelmaker(msa, errbuf, esl_opt_GetBoolean(go, "--rf"), esl_opt_GetReal(go, "--gapthresh"), &cm, &mtr)) != eslOK) return status;
  if(cfg->be_verbose) fprintf(stdout, "done.\n");
  
  /* set the CM's null model, if rsearch mode, use the bg probs used to calc RIBOSUM */
  if(! esl_opt_IsDefault(go, "--rsearch")) CMSetNullModel(cm, cfg->fullmat->g); 
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
  for (idx = 0; idx < msa->nseq; idx++) {
    ESL_ALLOC(aseq, (msa->alen+1) * sizeof(char));
    esl_abc_Textize(msa->abc, msa->ax[idx], msa->alen, aseq);
    tr[idx] = Transmogrify(cm, mtr, msa->ax[idx], aseq, msa->alen);
    ParsetreeCount(cm, tr[idx], msa->ax[idx], msa->wgt[idx]);
    free(aseq);
  }
  cm->nseq     = msa->nseq;
  cm->eff_nseq = msa->nseq;

  /* ensure the dual insert states we will detach were populated with 0 counts */
  if(!(esl_opt_GetBoolean(go, "--nodetach")))
    {
      if(cfg->be_verbose) fprintf(stdout, "%-40s ... ", "Finding and checking dual inserts");
      cm_find_and_detach_dual_inserts(cm, 
				      TRUE,   /* Do check (END_E-1) insert states have 0 counts */
				      FALSE); /* Don't detach the states yet, wait til CM is priorified */
      if (cfg->be_verbose) fprintf(stdout, "done.\n");
    }
  /* set the EL self transition probability */
  cm->el_selfsc = sreLOG2(esl_opt_GetReal(go, "--elself"));

  /* set the cm->beta_W parameter, which is not used in cmbuild, but is used by cmcalibrate and
   * cmsearch (and possibly others) to set cm->W */
  cm->beta_W = esl_opt_GetReal(go, "--Wbeta");

  /* Before converting to probabilities, save a count vector file, if asked.
   * Used primarily for making data files for training priors.
   */
  if (cfg->cfp != NULL) { 
    if ((status = print_countvectors(cfg, errbuf, cm)) != eslOK) goto ERROR;
  }

  *ret_cm  = cm;
  if(ret_mtr == NULL) FreeParsetree(mtr);
  else *ret_mtr = mtr;
  if(ret_msa_tr == NULL) {
    for(idx = 0; idx < msa->nseq; idx++)
      FreeParsetree(tr[idx]);
    free(tr);
    tr = NULL;
  }
  else *ret_msa_tr = tr;

  return eslOK;

 ERROR:
  return status;
}


/* set_model_name()
 * Give the model a name based on the MSA name.
 * 
 * We've ensured the msa has a name in name_msa() so if 
 * for some inconceivable reason it doesn't 
 * we die.
 *
 * note: This is much simpler than how HMMER3 does
 *       this. The reason is that the --ctarg --cmaxid
 *       cluster options produce N > 1 CM per MSA,
 *       which are named <msa->name>.1 .. <msa->name>.N.
 * 
 */
static int
set_model_name(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t *cm)
{
  int status;

  if (cfg->be_verbose) {
    fprintf(stdout, "%-40s ... ", "Set model name");
    fflush(stdout);
  }

  if(cm_SetName(cm, msa->name) != eslOK) goto ERROR;
  if (cfg->be_verbose) fprintf(stdout, "done. [%s]\n", cm->name);
  return eslOK;

 ERROR:
  if (cfg->be_verbose) fprintf(stdout, "FAILED.\n");
  return status;
}

/* set_model_cutoffs()
 * If the msa had them available, set the Rfam
 * cutoffs in the model.
 * 
 * Always returns eslOK;
 * 
 */
static int
set_model_cutoffs(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t *cm)
{
  int cutoff_was_set = FALSE;
  if(msa->cutset[eslMSA_TC1]) { 
    cm->tc = msa->cutoff[eslMSA_TC1];
    cm->flags |= CMH_TC;
    cutoff_was_set = TRUE;
  }
  if(msa->cutset[eslMSA_GA1]) { 
    cm->ga = msa->cutoff[eslMSA_GA1];
    cm->flags |= CMH_GA;
    cutoff_was_set = TRUE;
  }
  if(msa->cutset[eslMSA_NC1]) { 
    cm->nc = msa->cutoff[eslMSA_NC1];
    cm->flags |= CMH_NC;
    cutoff_was_set = TRUE;
  }
  if (cfg->be_verbose && cutoff_was_set ) {
    fprintf(stdout, "%-40s ... ", "Set model cutoffs");
    fprintf(stdout, "done.\n");
  }
  return eslOK;
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
  int status;
  double neff;
  int used_hmm_etarget = FALSE;
  if(cfg->be_verbose) fprintf(stdout, "%-40s ... ", "Set effective sequence number");
  fflush(stdout);

  if((esl_opt_GetBoolean(go, "--enone")) || (! esl_opt_IsDefault(go, "--rsearch")))
    {
      neff = msa->nseq;
      if(cfg->be_verbose) fprintf(stdout, "done. [--enone: neff=nseq=%d]\n", msa->nseq);
    }
  else if (esl_opt_GetBoolean(go, "--eent") == TRUE)
    {
      double etarget; 
      double hmm_etarget; 
      double hmm_re;
      int clen = 0;
      int nd;
      for(nd = 0; nd < cm->nodes; nd++) { 
	if(cm->ndtype[nd] == MATP_nd) clen += 2;
	else if(cm->ndtype[nd] == MATL_nd) clen += 1;
	else if(cm->ndtype[nd] == MATR_nd) clen += 1;
      }
      if (esl_opt_IsDefault(go, "--ere")) etarget = default_target_relent(cm->abc, clen, esl_opt_GetReal(go, "--eX"));
      else                                etarget = esl_opt_GetReal(go, "--ere");

      status = cm_EntropyWeight(cm, pri, etarget, FALSE, &hmm_re, &neff);
      /* if --ehmmre <x> enabled, ensure HMM relative entropy per match column is at least <x>, if not,
       * recalculate neff so HMM relative entropy of <x> is achieved.
       */
      if(! esl_opt_IsDefault(go, "--ehmmre")) { 
	hmm_etarget = esl_opt_GetReal(go, "--ehmmre"); 
	printf("cm hmm re: %f target: %f\n", hmm_re, hmm_etarget);
	if(hmm_re < hmm_etarget) { 
	  status = cm_EntropyWeight(cm, pri, hmm_etarget, TRUE, &hmm_re, &neff); /* TRUE says: pretend model is an HMM for entropy weighting */
	  if      (status == eslEMEM) ESL_FAIL(status, errbuf, "memory allocation failed");
	  else if (status != eslOK)   ESL_FAIL(status, errbuf, "internal failure in entropy weighting algorithm");
	  used_hmm_etarget = TRUE;
	}
      }
      if      (status == eslEMEM) ESL_FAIL(status, errbuf, "memory allocation failed");
      else if (status != eslOK)   ESL_FAIL(status, errbuf, "internal failure in entropy weighting algorithm");
      cm->eff_nseq = neff;
      cm_Rescale(cm, neff / (float) msa->nseq);

      if(cfg->be_verbose) { 
	if(used_hmm_etarget) fprintf(stdout, "done. [etarget (hmm) %.2f bits; neff %.2f]\n", hmm_etarget, neff);
	else                 fprintf(stdout, "done. [etarget (cm)  %.2f bits; neff %.2f]\n", etarget, neff);
      }
    }
  return eslOK;
}

/* parameterize()
 * Converts counts to probability parameters.
 */
static int
parameterize(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, const Prior_t *prior)
{
  int status; 

  if (cfg->be_verbose){
    fprintf(stdout, "%-40s ... ", "Converting counts to probabilities"); 
    fflush(stdout);
  }
  PriorifyCM(cm, prior);
  if(! (esl_opt_IsDefault(go, "--rsearch"))) {
    rsearch_CMProbifyEmissions(cm, cfg->fullmat); /* use those probs to set CM probs from cts */
    /*debug_print_cm_params(cm);*/
  }
  
  if(!esl_opt_GetBoolean(go, "--nodetach")) /* Detach dual inserts where appropriate, if
					     * we get here we've already checked these states */
    {
      cm_find_and_detach_dual_inserts(cm, 
				      FALSE, /* Don't check states have 0 counts (they won't due to priors) */
				      TRUE); /* Detach the states by setting trans probs into them as 0.0   */
    }
  if(! esl_opt_GetBoolean(go, "--iins")) { /* set all insert emission probabilities equal to the cm->null probabilities */ 
    if((status = flatten_insert_emissions(cm)) != eslOK) return status; 
    /* Note: flatten_insert_emissions() is purposefully a static function local to cmbuild.c b/c once CM files are calibrated
     * no other executable (i.e. cmsearch) should be able to modify the scores of the CM, as that would invalidate the Gumbels */
  }

  CMRenormalize(cm);
  CMLogoddsify(cm);

  if (cfg->be_verbose) fprintf(stdout, "done.\n");
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
  default:        cm_Fail("ERROR in default_target_relent(), alphabet not RNA!\n");
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
 * Give a MSA a name if it doesn't have one.
 * If -n <s> was enabled, name it <s>, else the 
 * naming rule is the suffixless name of the file it came from,
 * plus a "-<X>" with <X> = number MSA in the file.
 *
 * For example the 3rd MSA in file "alignments.stk" would be
 * named "alignments-3".
 */
int
name_msa(const ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa, int nali)
{
  int status;
  char *name = NULL;
  char *buffer = NULL;
  void *tmp;
  int n;
  int maxintlen;

  if(msa == NULL) ESL_FAIL(status, errbuf, "name_msa(), msa is NULL.");

  if(esl_opt_IsDefault(go, "-n") && msa->name != NULL) return eslOK; /* keep the msa's existing name */

  if(! (esl_opt_IsDefault(go, "-n"))) { /* give the msa the -n name */
    if(nali > 1) ESL_FAIL(eslEINCOMPAT, errbuf, "The -n option requires exactly 1 alignment, but the alignment file has > 1 alignments.");
    if((status = esl_strdup(esl_opt_GetString(go, "-n"), -1, &name)) != eslOK) ESL_FAIL(status, errbuf, "name_msa(), esl_strdup, memory allocation error.");
  }
  /* give the msa a name, the name of the file it comes from without the filetail, plus an appended "-X" where
   * X is the index of the alignment in the file 
   */
  else { 
    esl_FileTail(esl_opt_GetArg(go, 2), TRUE, &name); /* TRUE=nosuffix */
    if (name == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "Error getting file tail of the MSA.\n");
    n         = strlen(name);
    maxintlen = IntMaxDigits() + 1;  /* IntMaxDigits() returns number of digits in INT_MAX */
    ESL_ALLOC(buffer, sizeof(char) * maxintlen);
    sprintf(buffer, "-%d", (nali));
    n += strlen(buffer);
    ESL_RALLOC(name, tmp, sizeof(char)*(n+1));
    if((status = esl_strcat(&name, -1, buffer, -1)) != eslOK) goto ERROR;
  }

  if((status = esl_strdup(name, -1, &(msa->name))) != eslOK) ESL_FAIL(status, errbuf, "name_msa(), esl_strdup, memory allocation error.");

  if(name   != NULL) free(name);
  if(buffer != NULL) free(buffer);
  return eslOK;

 ERROR:
  if(name != NULL)   free(name);
  if(buffer != NULL) free(buffer);
  return status;
}

/* Function: print_countvectors()
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
print_countvectors(const struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int   v,x;

  if(cfg->cfp == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "save_countvectors(), but cfg->cfp is NULL, shouldn't happen.");

  /* Print emission counts */
  for (v = 0; v < cm->M; v++) {
    if (cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == MR_st) { 
      fprintf(cfg->cfp, "E\t%-7s ", UniqueStatetype(cm->stid[v]));
      if (cm->sttype[v] == MP_st) {
	for (x = 0; x < cm->abc->K*cm->abc->K; x++)
	  fprintf(cfg->cfp, "%8.3f ", cm->e[v][x]);
      } else {
	for (x = 0; x < cm->abc->K; x++)
	  fprintf(cfg->cfp, "%8.3f ", cm->e[v][x]);
      }
      fprintf(cfg->cfp, "\n");
    }
  }

  /* Print transition counts */
  for (v = 0; v < cm->M; v++) {
    if(cm->sttype[v] != B_st && cm->sttype[v] != E_st) {
      fprintf(cfg->cfp, "T\t%-7s : %-2d", UniqueStatetype(cm->stid[v]), cm->ndtype[(cm->ndidx[v] + 1)]);
      for (x = 0; x < cm->cnum[v]; x++) {
	fprintf(cfg->cfp, "%8.3f ", cm->t[v][x]);
      }
      fprintf(cfg->cfp, "\n");
    }
  }
  fprintf(cfg->cfp, "//\n");
  return eslOK;
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

void
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
  cm_Fail("Memory allocation error.");
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
      if(sq[i] == NULL) goto ERROR;
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
  else if(esl_opt_GetBoolean(go, "--optacc")) cm->align_opts  |= CM_ALIGN_OPTACC;

  if(esl_opt_GetBoolean(go, "--nonbanded"))   { 
    cm->align_opts  |= CM_ALIGN_SMALL; 
    cm->align_opts &= ~CM_ALIGN_OPTACC; /* turn optimal accuracy OFF */
  }
  else                                        cm->align_opts  |= CM_ALIGN_HBANDED;

  if(esl_opt_GetBoolean(go, "--sub"))         cm->align_opts  |= CM_ALIGN_SUB;
  if(esl_opt_GetBoolean(go, "--fins"))        cm->align_opts  |= CM_ALIGN_FLUSHINSERTS;

  /* update cm->config_opts */
  if(esl_opt_GetBoolean(go, "-l"))
    {
      cm->config_opts |= CM_CONFIG_LOCAL;
      cm->config_opts |= CM_CONFIG_HMMLOCAL;
      cm->config_opts |= CM_CONFIG_HMMEL;
    }

  /* finally, configure the CM for alignment based on cm->config_opts and cm->align_opts.
   * this may make a cp9 HMM, for example.
   */
  ConfigCM(cm, FALSE); /* FALSE says don't calc W */

  return eslOK;
}


/* Function: MSADivide()
 * EPN, Wed Mar 21 17:26:39 2007
 * 
 * Purpose:  Given an MSA, divide it into multiple MSAs, each with
 *           a different cluster of the original sequences. Each
 *           MSA will be used to construct a separate CM.
 *
 *           Different modes:
 *           
 *        1. if(do_all): each seq is its own cluster, so
 *           the number of new MSAs is number of seqs in input 
 *           master MSA. 
 *
 *        2. if(do_mindiff): define clusters
 *           such that we maximize the number of clusters while
 *           satisfying: minimum fractional difference b/t any 
 *           2 seqs in different clusters >= 'mindiff'. 
 *           The contract states that mindiff > 0. in this case.
 *           
 *        3. if(do_nc).: define clusters 
 *           such that we have exactly 'target_nc' clusters by
 *           searching for the 'mindiff' that gives exactly
 *           'target_nc' clusters. (We guarantee we can do this
 *           by rounding 'diff' fractional difference values b/t
 *           seqs to nearest 0.001). 
 *
 *        *. (NOT YET IMPLEMENTED)
 *           if(do_pickone): in mode 2 or 3, we select a single
 *           sequence from each cluster to represent that cluster. 
 *           The sequence is chosen that has the minimum average
 *           fractional difference with all other seqs in the cluster.
 *           (NOT YET IMPLEMENTED)
 *
 * Args:    
 * ESL_MSA *mmsa         - the master MSA, we cluster the seqs in this guy
 *                        and build a new MSA from each cluster
 * int     do_all       - TRUE (mode 1): each seq is its own cluster
 * int     do_mindiff   - TRUE (mode 2): satisfy clusters are at least mindiff different
 * int     do_nc        - TRUE (mode 3): set mindiff such that we get excatly target_nc clusters
 * float   mindiff      - the minimum fractional difference allowed
 *                        between 2 seqs of different clusters
 *                        (0. indicates mode 3) 
 * int     target_nc    - number of clusters to define (0 indicates mode 2)
 * int     do_orig      - TRUE to include the master MSA as one of the new MSAs
 * int    *ret_num_msa  - number of MSAs in ret_MSA
 * ESL_MSA  ***ret_cmsa - new MSAs, one for each cluster
 * char     *errbuf     - buffer for error messages
 *           
 * Return: ret_cmsa (alloc'ed here) and ret_num_msa
 */
int 
MSADivide(ESL_MSA *mmsa, int do_all, int do_mindiff, int do_nc, float mindiff, int target_nc,
	  int do_orig, int *ret_num_msa, ESL_MSA ***ret_cmsa, char *errbuf)
{
  int   status;        /* Easel status code */
  ESL_MSA **cmsa = NULL;/* the new MSAs we're creating from clusters of seqs in mmsa */
  int   i;             /* counter over sequences */
  int   m;             /* counter over new MSAs */
  int   n;             /* counter over tree nodes */
  ESL_TREE    *T = NULL;/* the tree, created by Single-Linkage Clustering */
  ESL_DMATRIX *D = NULL;/* the distance matrix */
  double *diff = NULL; /* [0..T->N-2], diff[n]= min distance between any leaf in right and
		        * left subtree of node n of tree T */
  double *minld = NULL;/* [0..T->N-2], min dist from node to any taxa in left  subtree */
  double *minrd = NULL;/* [0..T->N-2], min dist from node to any taxa in right subtree */
  int     nc;          /* number of clusters/MSAs  */
  int    *clust = NULL;/* [0..T->N-1], cluster number (0..nc-1) this seq is in */
  int    *csize = NULL;/* [0..nc-1], size of each cluster */
  int   **useme = NULL;/* [0.m.nc-1][0.i.N] TRUE to use seq i in new MSA m, FALSE not to */
  int     best;        /* 'best' node, returned by select_node() */
  void   *tmp;
  char   *buffer = NULL;
  int     ndigits;

  /* Contract check */
  if((do_all + do_nc + do_mindiff) != 1) ESL_FAIL(eslEINCOMPAT, errbuf, "MSADivide() exactly 1 of do_all, do_nc, do_mindiff must be TRUE.");
  if( do_nc && target_nc == 0)           ESL_FAIL(eslEINCOMPAT, errbuf, "MSADivide() target_nc is 0 but do_nc is TRUE!");
  if( do_mindiff && mindiff <= 0.)       ESL_FAIL(eslEINCOMPAT, errbuf, "MSADivide() mindiff is <= 0. but do_mindiff is TRUE!");
  if( do_mindiff && target_nc != 0)      ESL_FAIL(eslEINCOMPAT, errbuf, "MSADivide() do_mindiff is TRUE, but target_nc != 0");
  /* mmsa must be digital */
  if(!(mmsa->flags & eslMSA_DIGITAL))                 ESL_FAIL(eslEINCOMPAT, errbuf, "MSADivide() MSA is not digital.");

  if(do_nc) mindiff = 0.;

  /* Mode 1: Each seq becomes own MSA. Easy. */
  if(do_all) {
    ESL_ALLOC(clust, sizeof(int) * (mmsa->nseq));
    ESL_ALLOC(csize, sizeof(int) * (mmsa->nseq));
    nc = 0;
    /* each seq is its own cluster */
    for(i = 0; i < mmsa->nseq; i++) {
      clust[i] = nc++;
      csize[i] = 1;
    }
    printf("# Alignment split into %d clusters; each comprised of exactly 1 sequence\n", nc);
    printf("#\n");
   }
  else { /* Mode 2 or Mode 3 */ 
    /* Create distance matrix and infer tree by single linkage clustering */
    if((status = esl_dst_XDiffMx(mmsa->abc, mmsa->ax, mmsa->nseq, &D)) != eslOK) ESL_FAIL(status, errbuf, "esl_dst_XDiffMx() error, status: %d", status);
    if((status = esl_tree_SingleLinkage(D, &T)) != eslOK)                        ESL_FAIL(status, errbuf, "esl_tree_SingleLinkage() error, status: %d", status);
    if((status = esl_tree_SetTaxaParents(T)) != eslOK)                           ESL_FAIL(status, errbuf, "esl_tree_SetTaxaParentse() error, status: %d", status);
    /*esl_tree_WriteNewick(stdout, T);*/
    if((status = esl_tree_Validate(T, errbuf) != eslOK)) return status;
    
    /* determine the diff values: 
     * (use: n_child > n, unless n's children are taxa)
     * diff[n] is minimum distance between any taxa (leaf) in left subtree of 
     * n to any taxa in right subtree of n. 
     */
    ESL_ALLOC(diff,  (sizeof(double) * (T->N - 1)));  /* one for each node */
    ESL_ALLOC(minld, (sizeof(double) * (T->N - 1))); 
    ESL_ALLOC(minrd, (sizeof(double) * (T->N - 1))); 
    for (n = (T->N-2); n >= 0; n--) {
      minld[n] = T->ld[n] + ((T->left[n]  > 0) ? (minld[T->left[n]])  : 0);
      minrd[n] = T->rd[n] + ((T->right[n] > 0) ? (minrd[T->right[n]]) : 0);
      diff[n]  = minld[n] + minrd[n];
      diff[n] *= 1000.; 
      diff[n]  = (float) ((int) diff[n]);
      diff[n] /= 1000.; 
      /*printf("diff[n:%d]: %f\n", n, diff[n]);*/
    }
    free(minld); minld = NULL;
    free(minrd); minrd = NULL;
    /*for (n = 0; n < (T->N-1); n++)
      printf("diff[n:%d]: %f\n", n, diff[n]);
      for (n = 0; n < (T->N-1); n++)
      printf("left[n:%d]: %d right[n:%d]: %d\n", n, T->left[n], n, T->right[n]);*/
    
    if(do_mindiff) { /* Mode 2 */
      /* Define clusters that are at least mindiff different
       * from each other. */
      if((status = select_node(T, diff, mindiff, &clust, &nc, &best, errbuf)) != eslOK) return status;
      printf("# Alignment split into %d clusters; each will be used to train a CM.\n", nc);
      printf("# Maximum identity b/t any 2 seqs in different clusters: %.2f\n", (1.-mindiff));
      printf("#\n");
    }
    else { /* Mode 3, do_nc == TRUE, mindiff was set to 0.0 above */
      /* Find the minimum fractional difference (mindiff) that 
       * gives exactly target_nc clusters, also define clusters
       * based on that mindiff, this is all done with find_mindiff(),
       * which does a binary search for mindiff, we're guaranteed to 
       * find exactly target_nc clusters b/c diff values are rounded
       * to nearest 0.001. */
      if(target_nc > (T->N)) target_nc = T->N; /* max num clusters is num seqs */
      if((status = find_mindiff(T, diff, target_nc, &clust, &nc, &mindiff, errbuf)) != eslOK) return status;
      printf("# Alignment split into %d clusters; each will be used to train a CM.\n", nc);
      printf("# Maximum identity b/t any 2 seqs in different clusters: %.2f\n", (1.-mindiff));
      printf("#\n");
    }
    /* Determine the size of each cluster */
    ESL_ALLOC(csize, (sizeof(int) * (nc)));
    esl_vec_ISet(csize, nc, 0);
    for(i = 0; i < mmsa->nseq; i++)
      csize[clust[i]]++;
    
    /*printf("Distance matrix:\n");
      esl_dmatrix_Dump(stdout, D, NULL, NULL);*/
  }
  
  /* Create one new MSA for each cluster,
   * if(do_orig): keep the original MSA as cmsa[nc] */
  if(do_orig) {
    ESL_ALLOC(cmsa, (sizeof(ESL_MSA *) * (nc+1))); 
    for(m = 0; m < nc; m++) cmsa[m] = NULL;
  }
  else {
    ESL_ALLOC(cmsa, (sizeof(ESL_MSA *) * (nc)));
    for(m = 0; m < nc; m++) cmsa[m] = NULL;
  }
  ESL_ALLOC(useme, (sizeof(int *) * (nc+1)));
  for(m = 0; m <= nc; m++) {
    ESL_ALLOC(useme[m], (sizeof(int)) * mmsa->nseq);
    if(m < nc) esl_vec_ISet(useme[m], mmsa->nseq, FALSE);
    else       esl_vec_ISet(useme[m], mmsa->nseq, TRUE); /* keep all seqs in cmsa[nc]*/
  }
  
  for(i = 0; i < mmsa->nseq; i++)
    if(clust[i] != -1) 
      useme[clust[i]][i] = TRUE;
  ESL_ALLOC(buffer, sizeof(char) * (IntMaxDigits() + 1));  /* IntMaxDigits() returns number of digits in INT_MAX */
  for(m = 0; m < nc; m++) {
    if((status = esl_msa_SequenceSubset(mmsa, useme[m], &(cmsa[m]))) != eslOK) ESL_FAIL(status, errbuf, "MSADivide(), esl_msa_SequenceSubset error, status: %d.", status);
    /* rename the MSA it by adding ".<m+1>" */
    if(cmsa[m]->name == NULL) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "MSADivide(), an msa's name is NULL, shouldn't happen.");
    ndigits  = strlen(cmsa[m]->name);
    ndigits += sprintf(buffer, ".%d", (m+1));
    ESL_RALLOC(cmsa[m]->name, tmp, sizeof(char)*(ndigits+1));
    if ((status = esl_strcat(&cmsa[m]->name, -1, buffer, (ndigits+1))) != eslOK) goto ERROR;
    if ((status = esl_strchop(cmsa[m]->name, ndigits)) != eslOK) goto ERROR;
    free(useme[m]);
  }
  if(do_orig) {
    if((status = esl_msa_SequenceSubset(mmsa, useme[nc], &(cmsa[nc]))) != eslOK) ESL_FAIL(status, errbuf, "MSADivide(), esl_msa_SequenceSubset error, status: %d.", status);
  }
  free(useme[nc]);
  free(useme);
  
  if(do_orig) *ret_num_msa = nc+1;
  else        *ret_num_msa = nc;
  *ret_cmsa = cmsa;
  
  if(!do_all) { /* else we didn't allocate these structures */
    esl_tree_Destroy(T);
    esl_dmatrix_Destroy(D);
    free(diff);
    diff = NULL;
  }
  if(clust != NULL)  free(clust);
  if(csize != NULL)  free(csize);
  if(buffer!= NULL)  free(buffer);
  
  return eslOK;
  
 ERROR: 
  if(diff  != NULL) free(diff);
  if(minld != NULL) free(minld);
  if(minrd != NULL) free(minrd);
  if(clust != NULL) free(clust);
  if(csize != NULL) free(csize);
  if(buffer!= NULL) free(buffer);
  if(cmsa  != NULL) {
    for(m = 0; m < nc; m++)
      if(cmsa[m] != NULL) esl_msa_Destroy(cmsa[m]);
    free(cmsa);
  }
  return status;
}

/* Function: select_node()
 * EPN, Fri Mar 23 08:48:37 2007 
 * Adapted from SRE's select_node() in maketestset.c originally written
 * for the PROFMARK HMMER benchmark.
 * 
 * 
 * Purpose:  Define clusters of the taxa (seqs) in the tree such
 *           that minimum disparity b/t any 2 seqs in different 
 *           clusters is greater than <mindiff> and the number of
 *           clusters is maximized. <ret_best> is the index of the node
 *           of the tree under which the largest cluster belongs.
 *           <ret_nc> is the number of clusters after clustering, 
 *           <ret_clust> is an array [0..T->N-1] specifying which
 *           cluster each taxa belongs to.
 *           
 *           For high disparities, this cluster may contain all
 *           the sequences, and we'll return the root node (0).
 *
 * Args:    
 * ESL_TREE *T        - the tree
 * double   *diff     - [0..T->N-2]: for each node of the tree, the minimum
 *                      distance (sum of branch lengths) from any taxa (leaf)
 *                      in left subtree to any taxa in right subtree.
 * double    mindiff  - (see description above)
 * int     **ret_clust- [0..T->N-1] cluster number this seq is in, alloc'ed, filled here
 * int      *ret_nc   - number of clusters
 * int      *ret_best - RETURN: index of node of tree under which largest cluster belongs (see Purpose).
 * char     *errbuf   - buffer for error messages
 *
 * Returns: node index (as explained in Purpose)
 */
static int
select_node(ESL_TREE *T, double *diff, double mindiff, int **ret_clust, int *ret_nc, int *ret_best, char *errbuf)
{
  int status;     /* Easel status code */
  ESL_STACK *ns1; /* stack for traversing tree */
  ESL_STACK *ns2; /* another stack for traversing tree */
  int c;	  /* counter for clusters */
  int best;       /* index of current best node */
  int maxsize;    /* size of cluster for best node */
  int n, np;      /* counters over tree nodes */
  int *clust;     /* [1..T->N-1] cluster number this seq is in */

  /*printf("in selec_node mindiff: %f T->N: %d\n", mindiff, T->N);*/
  /* set tree cladesizes if not already set */
  if(T->cladesize == NULL) 
    if((status = esl_tree_SetCladesizes(T)) != eslOK) ESL_FAIL(status, errbuf, "select_node(), esl_tree_SetCladeSizes error, status: %d.", status);

  ESL_ALLOC(clust, (sizeof(int) * T->N));
  esl_vec_ISet(clust, T->N, 0);

  if((ns1 = esl_stack_ICreate()) == NULL) ESL_FAIL(status, errbuf, "select_node(), failed to create a stack, probably out of memory, status: %d.", status);
  if((ns2 = esl_stack_ICreate()) == NULL) ESL_FAIL(status, errbuf, "select_node(), failed to create a stack, probably out of memory, status: %d.", status);

  /* push root on stack to start */
  if((status = esl_stack_IPush(ns1, 0)) != eslOK) ESL_FAIL(status, errbuf, "select_node(), failed to push onto a stack, probably out of memory, status: %d.", status);	
  maxsize  = 0;
  best     = 0;
  c        = 0;
  while (esl_stack_IPop(ns1, &n) != eslEOD) {
    if ((n == 0 || diff[T->parent[n]] > mindiff) &&
	diff[n] <= mindiff) { /* we're at a cluster */
      if (T->cladesize[n] > maxsize) {
	maxsize = T->cladesize[n];
	best = n;
      }
      /* determine all taxa in the clade rooted at n*/
      esl_stack_IPush(ns2, n);	
      while (esl_stack_IPop(ns2, &np) != eslEOD) {
	/*printf("np: %d T->left[np]: %d\n", np, T->left[np]);*/
	if(T->left[np]  <= 0) clust[(-1*T->left[np])]  = c;
	else { if((status = esl_stack_IPush(ns2, T->left[np])) != eslOK) ESL_FAIL(status, errbuf, "select_node(), failed to push onto a stack, probably out of memory, status: %d.", status); }
	if(T->right[np] <= 0) clust[(-1*T->right[np])]  = c;
	else { if((status = esl_stack_IPush(ns2, T->right[np])) != eslOK) ESL_FAIL(status, errbuf, "select_node(), failed to push onto a stack, probably out of memory, status: %d.", status); }
      }
      c++;
    }
    else {		/* we're not a cluster, keep traversing */
      /*printf("n: %d T->left[n]: %d\n", n, T->left[n]);*/
      if(T->left[n]  <= 0) clust[(-1*T->left[n])]  = c++; /* single seq with its own cluster */
      else { if((status = esl_stack_IPush(ns1, T->left[n])) != eslOK) ESL_FAIL(status, errbuf, "select_node(), failed to push onto a stack, probably out of memory, status: %d.", status); }
      if(T->right[n] <= 0) clust[(-1*T->right[n])] = c++; /* single seq with its own cluster */
      else { if((status = esl_stack_IPush(ns1, T->right[n])) != eslOK) ESL_FAIL(status, errbuf, "select_node(), failed to push onto a stack, probably out of memory, status: %d.", status); }
    }
  }
  esl_stack_Destroy(ns1);
  esl_stack_Destroy(ns2);
  *ret_nc = c;
  *ret_clust = clust;
  *ret_best  = best;
  /*printf("nc: %d(%d) best: %d maxsize: %d nc: %d\n\n", *ret_nc, c, best, maxsize, c);
    for(n = 0; n < T->N; n++)
    printf("clust[%d]: %d\n", n, clust[n]);*/
  return eslOK;

 ERROR: 
  if(clust != NULL) free(clust);
  ESL_FAIL(status, errbuf, "select_node(), memory allocation error, status: %d.", status); 
}


/* Function: find_mindiff()
 * EPN, Fri Mar 23 18:59:42 2007
 * 
 * Purpose:  Given a tree resulting from single linkage clustering,
 *           find the min fractional difference (mindiff) that when used to
 *           define clusters (such that no seq in cluster A is less
 *           than mindiff different than any seq in cluster B), 
 *           gives >= target_nc.
 *
 * Args:    
 * ESL_TREE *T        - the tree
 * double   *diff     - [0..T->N-2]: for each node of the tree, the minimum
 *                      distance (sum of branch lengths) from any taxa (leaf)
 *                      in left subtree to any taxa in right subtree.
 * int      target_nc - number of clusters we want
 * int     **ret_clust- [0..T->N-1] cluster number this seq is in, alloc'ed, filled here
 * int      *ret_nc   - number of clusters
 * char     *errbuf   - buffer for error messages
 *
 * Returns: fractional difference (as explained in Purpose)
 */
static float
find_mindiff(ESL_TREE *T, double *diff, int target_nc, int **ret_clust, int *ret_nc, float *ret_mindiff, char *errbuf)
{
  int   status;
  float high       = 1.0;
  float low        = 0.0;
  int   high_nc    = 0;
  int   low_nc     = 0;
  float mindiff    = 0.5;
  int   curr_nc    = -1;
  int   curr_best  = -1;
  int   keep_going = TRUE;
  float thresh     = 0.001;
  int  *clust      = NULL;

  /* Contract check */
  if(target_nc > T->N) ESL_FAIL(eslEINCOMPAT, errbuf, "find_mindiff(), desired number of clusters is greater than number of seqs in the tree");

  while(keep_going) {
    if(clust != NULL) free(clust);
    if((status = select_node(T, diff, mindiff, &clust, &curr_nc, &curr_best, errbuf)) != eslOK) return status;
    if(curr_nc < target_nc) {
      high       = mindiff;
      high_nc    = curr_nc;
      mindiff   -= (mindiff - low) / 2.;
      if((fabs(high-0.) < thresh) && (fabs(low-0.) < thresh))  keep_going = FALSE; 
      /* stop, high and low have converged at 0. */
      /*printf("LOWER   nc: %d mindiff: %f low: %f high: %f\n", curr_nc, mindiff, low, high);*/
    }
    else {/* curr_nc >= target_nc */
      low        = mindiff;
      low_nc     = curr_nc;
      mindiff   += (high - mindiff) / 2.;
      if(fabs(high-low) < thresh)  keep_going = FALSE; /* stop, high and low have converged */
      /*printf("GREATER nc: %d mindiff: %f low: %f high: %f\n", curr_nc, mindiff, low, high);*/
    }
  }
  /* it's possible we can't reach our target, if so, set mindiff as minimum value that gives 
   * less than target_nc clusters. */
  if(curr_nc != target_nc) {
    /*printf("targ: %d curr: %d low: %d (%f) high: %d (%f)\n", target_nc, curr_nc, low_nc, low, high_nc, high);*/
    if(high_nc < target_nc) {
      mindiff = high;
      if((status = select_node(T, diff, mindiff, &clust, &curr_nc, &curr_best, errbuf)) != eslOK) return status;
    }
    else
      while(high_nc > target_nc) {
	high += thresh;
	if(high > 1.0)  ESL_FAIL(eslEINCONCEIVABLE, errbuf, "find_mindiff(), mindiff has risen above 1.0");
	mindiff = high;
	if((status = select_node(T, diff, mindiff, &clust, &curr_nc, &curr_best, errbuf)) != eslOK) return status;
	high_nc = curr_nc;
      }
  }
  /*printf("FINAL mindiff: %f\n", mindiff);  */
  *ret_nc    = curr_nc;
  *ret_clust = clust;
  *ret_mindiff = mindiff;

  return eslOK;
}

/* Function: write_cmbuild_info_to_comlog
 * Date:     EPN, Mon Dec 31 10:17:21 2007
 *
 * Purpose:  Set the cmbuild command info and creation date info 
 *           in a ComLog_t data structure.
 *           clog->bcom, clog->bdate should be NULL when we enter this function.
 *
 * Returns:  eslOK on success, eslEINCOMPAT on contract violation.
 */
static int
write_cmbuild_info_to_comlog(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  int status;
  int i;
  long seed;
  long temp;
  int  seedlen;
  char *seedstr = NULL;
  time_t date = time(NULL);

  if(cfg->comlog->bcom != NULL)  ESL_FAIL(eslEINCOMPAT, errbuf, "write_cmbuild_info_to_comlog(), cfg->comlog->bcom is non-NULL.");
  if(cfg->comlog->bdate != NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "write_cmbuild_info_to_comlog(), cfg->comlog->bcom is non-NULL.");

  /* Set the cmbuild command info, the cfg->comlog->bcom string */
  for (i = 0; i < go->optind; i++) { /* copy all command line options, but not the command line args yet, we may need to append '-s ' before the args */
    esl_strcat(&(cfg->comlog->bcom), -1, go->argv[i], -1);
    esl_strcat(&(cfg->comlog->bcom), -1, " ", 1);
  }
  /* if --gibbs enabled, and -s NOT enabled, we need to append the seed info also */
  if(esl_opt_GetBoolean(go, "--gibbs")) {
    if(cfg->r == NULL) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "write_cmbuild_info_to_comlog(), cfg->r is NULL but --gibbs enabled, shouldn't happen.");
    seed = esl_randomness_GetSeed(cfg->r);
    if(esl_opt_IsDefault(go, "-s")) {
      temp = seed; 
      seedlen = 1; 
      while(temp > 0) { temp/=10; seedlen++; } /* determine length of stringized version of seed */
      seedlen += 4; /* strlen(' -s ') */

      ESL_ALLOC(seedstr, sizeof(char) * (seedlen+1));
      sprintf(seedstr, " -s %ld ", seed);
      esl_strcat((&cfg->comlog->bcom), -1, seedstr, seedlen);
    }
    else { /* -s was enabled with --gibbs, we'll do a sanity check */
      if(seed != (long) esl_opt_GetInteger(go, "-s")) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "write_cmbuild_info_to_comlog(), cfg->r's seed is %ld, but -s was enabled with argument: %ld!, this shouldn't happen.", seed, (long) esl_opt_GetInteger(go, "-s"));
    }
  }

  for (i = go->optind; i < go->argc; i++) { /* copy all command line options, but not the command line args yet, we may need to append '-s ' before the args */
    esl_strcat(&(cfg->comlog->bcom), -1, go->argv[i], -1);
    if(i < (go->argc-1)) esl_strcat(&(cfg->comlog->bcom), -1, " ", 1);
  }

  /* Set the cmbuild creation date, the cfg->comlog->bdate string */
  if((status = esl_strdup(ctime(&date), -1, &(cfg->comlog->bdate))) != eslOK) goto ERROR;
  esl_strchop(cfg->comlog->bdate, -1); /* doesn't return anything but eslOK */

  if(seedstr != NULL) free(seedstr);
  return eslOK;

 ERROR:
  if(seedstr != NULL) free(seedstr);
  ESL_FAIL(status, errbuf, "write_cmbuild_info_to_comlog() error status: %d, probably out of memory.", status);
  return status; 
}

/* Function: flatten_insert_emissions()
 *
 * Purpose:  Set the insert emission *probabilities* of a CM to it's 
 *           null model probabilities. Subsequently in CMLogoddsify(),
 *           all insert emissions scores will become 0.0 bits.
 *           This option is called by default for all CMs unless --iins
 *           was enabled on the command line. It is impt that if we're
 *           going to zero insert scores, we do it within the CM file 
 *           (as opposed to previous versions of infernal which allowed us 
 *            to zero inserts with cmsearch, for example), because now we
 *           have E-values and if a Gumbel is fit to a CM with zeroed 
 *           inserts or with informative inserts, those insert scores should
 *           never change as long as that Gumbel is used.
 * 
 * Returns: eslOK on success.
 */
int
flatten_insert_emissions(CM_t *cm)
{
  int v;

  /* Contract check */
  if(cm->abc  == NULL) cm_Fail("flatten_insert_emissions(), cm->abc is NULL.\n");
  if(cm->null == NULL) cm_Fail("flatten_insert_emissions(), cm->null is NULL.\n");

  esl_vec_FNorm(cm->null, cm->abc->K);
  for (v = 0; v < cm->M; v++) {
    if(cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) {
      esl_vec_FSet(cm->e[v],    (cm->abc->K * cm->abc->K), 0.); /* zero them out */
      esl_vec_FCopy(cm->null, cm->abc->K, cm->e[v]); /* overwrite first cm->abc->K values (rest are irrelevant for non-MP states) with cm->null */
    }
  }
  return eslOK;
}

/* Function: print_run_info
 * Date:     EPN, Fri Feb 29 09:58:21 2008
 *
 * Purpose:  Print information on this run of cmbuild.
 *           Command used to run it, and execution date.
 *
 * Returns:  eslOK on success
 */
static int
print_run_info(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf)
{
  fprintf(stdout, "%-10s %s\n",  "# command:", cfg->comlog->bcom);
  fprintf(stdout, "%-10s %s\n",  "# date:",    cfg->comlog->bdate);
  fprintf(stdout, "#\n");
  return eslOK;
}

/* Function: print_refine_column_headings()
 * Date:     EPN, Fri Feb 29 10:36:03 2008
 *
 * Purpose:  Print column headings for tabular output of refine_msa() to 
 *           output file (stdout unless -o). 
 *
 * Returns:  eslOK on success
 */
static void
print_refine_column_headings(const ESL_GETOPTS *go, const struct cfg_s *cfg)
{
  fprintf(stdout, "#\n");
  fprintf(stdout, "# %-5s %-13s %10s\n", "iter",  "bit score sum", "fract diff");
  fprintf(stdout, "# %-5s %-13s %10s\n", "-----", "-------------", "----------");
  return;
}

/* Function: get_namewidth()
 * Date:     EPN, Fri May 23 05:45:32 2008
 *
 * Purpose:  Determine the maximum length of a CM name we'll create in the current
 *           cmbuild call from the MSA already opened cfg->afp.
 *           Sets cfg->namewidth as the max number of characters needed for
 *           all CM names.
 *
 * Returns:  eslOK on success
 *           eslEFORMAT on parse error of MSA in cfg->afp
 */
static int
get_namewidth(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  int status;
  int nali = 0;
  int do_cmaxid_or_call;
  int do_ctarget;
  int cur_namewidth;
  ESL_MSA *msa = NULL;

  cfg->namewidth = 6; /* length of "name", plus 2 spaces, just for looks */
  do_cmaxid_or_call = (! (esl_opt_IsDefault(go, "--cmaxid"))) || (esl_opt_GetBoolean(go, "--call")) ? TRUE : FALSE;
  do_ctarget        = (! esl_opt_IsDefault(go, "--ctarget")) ? TRUE : FALSE;

  /* if -n <s> enabled, set namewidth as length of <s> */
  if(!(esl_opt_IsDefault(go, "-n"))) { 
    cfg->namewidth = ESL_MAX(cfg->namewidth, strlen(esl_opt_GetString(go, "-n"))); 
    return eslOK;
  }

  /* else, get MSA names using either (1) stockholm GF ID markup name or (2) cmbuild's rules for naming the msa */
  while ((status = esl_msa_Read(cfg->afp, &msa)) != eslEOF) { 
      if      (status == eslEFORMAT)  ESL_FAIL(status, errbuf, "Alignment file parse error, line %d of file %s:\n%s\nOffending line is:\n%s\n", cfg->afp->linenumber, cfg->afp->fname, cfg->afp->errbuf, cfg->afp->buf);
      else if (status != eslOK)       ESL_FAIL(status, errbuf, "Alignment file read unexpectedly failed with code %d\n", status);
      nali++;

      /* name the msa, if it already has one from #=GF ID markup, name_msa() returns w/o modifying it */
      if((status = name_msa(go, errbuf, msa, nali)) != eslOK) return status;
      cur_namewidth = strlen(msa->name);
      if(do_cmaxid_or_call) cur_namewidth += 1 + IntDigits(msa->nseq); /* we could create as many as msa->nseq CMs from this msa */
      if(do_ctarget)        cur_namewidth += 1 + esl_opt_GetInteger(go, "--ctarget"); /* we'll make --ctarget <n> CMs from this msa */
      cfg->namewidth = ESL_MAX(cfg->namewidth, cur_namewidth);
      esl_msa_Destroy(msa);
  }
  /* close the MSA file and open it again, sloppy */
  esl_msafile_Close(cfg->afp);
  if   (esl_opt_IsDefault(go, "--informat")) cfg->fmt = eslMSAFILE_UNKNOWN; /* autodetect sequence file format by default. */ 
  else cfg->fmt = esl_sqio_FormatCode(esl_opt_GetString(go, "--informat"));
  status = esl_msafile_Open(cfg->alifile, cfg->fmt, NULL, &(cfg->afp));
  if      (status == eslENOTFOUND) ESL_FAIL(status, errbuf, "Alignment file %s doesn't exist or is not readable\n", cfg->alifile);
  else if (status == eslEFORMAT) ESL_FAIL(status, errbuf, "Couldn't determine format of alignment %s\n", cfg->alifile);
  else if (status != eslOK)      ESL_FAIL(status, errbuf, "Alignment file open failed with error %d\n", status);
  cfg->fmt = cfg->afp->format;
  esl_msafile_SetDigital(cfg->afp, cfg->abc);

  return eslOK;
}

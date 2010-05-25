/* cmalign.c
 * SRE, Thu Jul 25 11:28:03 2002 [St. Louis]
 * SVN $Id$
 * 
 * Align sequences to a CM.
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************
 */
#include "esl_config.h"
#include "p7_config.h"
#include "config.h"

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "easel.h"		/* general seq analysis library   */
#include "esl_alphabet.h"
#include "esl_getopts.h"		
#include "esl_mpi.h"
#include "esl_msa.h"
#include "esl_random.h"		
#include "esl_sq.h"		
#include "esl_sqio.h"
#include "esl_stack.h"
#include "esl_sse.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "funcs.h"		/* external functions                   */
#include "structs.h"		/* data structures, macros, #define's   */

#include "hmmer.h"
#if 0
#include "impl_sse.h"
#endif

#define ALGOPTS      "--cyk,--optacc,--viterbi,--sample"         /* Exclusive choice for algorithm */
#define BIGALGOPTS   "--cyk,--optacc,--viterbi,--sample,--small" /* Incompatible with --optacc,--sample (except their selves) */
#define ACCOPTS      "--nonbanded,--hbanded,--qdb"               /* Exclusive choice for acceleration strategies */
#define OUTALPHOPTS  "--rna,--dna"                               /* Exclusive choice for output alphabet */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs       incomp  help  docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",   1 },
  { "-o",        eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "output the alignment to file <f>, not stdout", 1 },
  { "-l",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "align locally w.r.t. the model",         1 },
  { "-q",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "quiet; suppress banner and scores, print only the alignment", 1 },
  { "-M",        eslARG_INFILE, NULL,  NULL, NULL,      NULL,      NULL,        NULL, "meta-cm mode: <cmfile> is a meta-cm built from aln in <f>", 1 },
  { "--ileaved", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,   "--chunk", "output alnment in interleaved Stockholm format (not 1 line/seq)",  1 },
  { "--no-prob", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "do not append posterior probabilities to alignment", 1 },
  { "--informat",eslARG_STRING, NULL,  NULL, NULL,      NULL,      NULL,        NULL, "specify the input file is in format <x>, not FASTA", 1 },
  { "--chunk",   eslARG_INT,  "1000",  NULL, NULL,      NULL,      NULL,        NULL, "num seqs for each temp alnment, for saving memory", 1 },
  { "--devhelp", eslARG_NONE,   NULL,  NULL, NULL,      NULL,      NULL,        NULL, "show list of undocumented developer options", 1 },
#ifdef HAVE_MPI
  { "--mpi",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,   "--sfile", "run as an MPI parallel program",                    1 },  
  { "--nseq-mpi",eslARG_INT,      "1", NULL, NULL,      NULL,   "--mpi",        NULL, "in MPI mode, num seqs in a workunit", 1 },  
#endif
  /* Algorithm options */
  { "--optacc",  eslARG_NONE,"default",NULL, NULL,     NULL,      NULL,  BIGALGOPTS, "align with the Holmes/Durbin optimal accuracy algorithm", 2 },
  { "--cyk",     eslARG_NONE,   FALSE, NULL, NULL,"--optacc",     NULL,     ALGOPTS, "align with the CYK algorithm", 2 },
  { "--sample",  eslARG_NONE,   FALSE, NULL, NULL,"--optacc",     NULL,  BIGALGOPTS, "sample alignment of each seq from posterior distribution", 2 },
  { "-s",        eslARG_INT,    "181", NULL, "n>=0",   NULL,"--sample",        NULL, "w/--sample, set RNG seed to <n> (if 0: one-time arbitrary seed)", 2 },
  { "--viterbi", eslARG_NONE,   FALSE, NULL, NULL,"--optacc","--no-prob",   ALGOPTS, "align to a CM Plan 9 HMM with the Viterbi algorithm",2 },
  { "--sub",     eslARG_NONE,   FALSE, NULL, NULL,     NULL,      NULL,        "-l", "build sub CM for columns b/t HMM predicted start/end points", 2 },
  { "--small",   eslARG_NONE,   FALSE, NULL, NULL,     NULL,"--cyk,--no-prob", "--hbanded", "use divide and conquer (d&c) alignment algorithm", 2 },
  { "--hbanded", eslARG_NONE,"default",NULL, NULL,     NULL,      NULL,     ACCOPTS, "accelerate using CM plan 9 HMM derived bands", 3 },
  { "--p7",      eslARG_NONE,   FALSE, NULL, NULL,     NULL,"--hbanded",    ACCOPTS, "really accelerate using plan 7 HMM derived bands", 3 },
  { "--nonbanded",eslARG_NONE,  FALSE, NULL, NULL,"--hbanded",    NULL,     ACCOPTS, "do not use bands to accelerate aln algorithm", 3 },
  { "--tau",     eslARG_REAL,   "1E-7",NULL, "0<x<1",  NULL,"--hbanded",       NULL, "set tail loss prob for --hbanded to <x>", 3 },
  { "--mxsize",  eslARG_REAL, "2048.0",NULL, "x>0.",   NULL,      NULL,"--small,--qdb", "set maximum allowable DP matrix size to <x> Mb", 3},
  /* Options that modify how the output alignment is created */
  { "--rna",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, OUTALPHOPTS, "output alignment as RNA sequence data", 4},
  { "--dna",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, OUTALPHOPTS, "output alignment as DNA (not RNA) sequence data", 4},
  { "--matchonly",eslARG_NONE,  FALSE, NULL, NULL,      NULL,      NULL,        NULL, "include only match columns in output alignment", 4 },
  /*  { "--resonly", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "include only match columns with >= 1 residues in output aln", 4 },*/
  /*  { "--fins",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "flush inserts left/right in output alignment", 4 }, */
  /* Including a preset alignment */
  { "--withali", eslARG_INFILE, NULL,  NULL, NULL,      NULL,      NULL,"--viterbi","incl. alignment in <f> (must be aln <cm file> was built from)", 5 },
  { "--withpknots",eslARG_NONE, NULL,  NULL, NULL,      NULL,"--withali",       NULL, "incl. structure (w/pknots) from <f> from --withali <f>", 5 },
  { "--rf",      eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--withali",       NULL, "--rf was originally used with cmbuild", 5 },
  { "--gapthresh",eslARG_REAL,  "0.5", NULL, "0<=x<=1", NULL,"--withali",       NULL, "--gapthresh <x> was originally used with cmbuild", 5 },
  /* Using only a single CM from a multi-CM file */
  { "--cm-idx",   eslARG_INT,   NULL, NULL,  "n>0",     NULL,      NULL, "--cm-name", "only align seqs with CM number <n>    in the CM file",  10 },
  { "--cm-name",  eslARG_STRING,NULL, NULL,  NULL,      NULL,      NULL,  "--cm-idx", "only align seqs with the CM named <s> in the CM file",  10 },
  /* Verbose output files */
  { "--tfile",   eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "dump individual sequence parsetrees to file <f>", 7 },
  { "--ifile",   eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "dump information on per-sequence inserts to file <f>", 7 },
  { "--elfile",  eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      "-l",        NULL, "dump information on per-sequence EL inserts to file <f>", 7 },
  { "--sfile",   eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "dump alignment score information to file <f>", 7 },
  /* options for experimental p7 HMM banding */
  { "--7pad",    eslARG_INT,       "0", NULL, "n>=0",   NULL,      NULL,        NULL, "w/p7 banding set pin pad to <n> residues", 99 },
  { "--7len",    eslARG_INT,       "4", NULL, "n>0",    NULL,      NULL,        NULL, "w/p7 banding set minimum length pin n-mer to <n>", 99 },
  { "--7sc",     eslARG_REAL,    "0.5", NULL, "x>-0.0001",NULL,    NULL,        NULL, "w/p7 banding set minimum pin score to <x>", 99 },
  { "--7end",    eslARG_INT,       "0", NULL, "n>=0",   NULL,      NULL,        NULL, "w/p7 banding remove pins within <n> residues of k-mer termini", 99 },
  { "--7mprob",  eslARG_REAL,    "0.0", NULL, "x>-0.0001",NULL,    NULL,        NULL, "w/p7 banding set min prob to enter match state pin to <x>", 99 },
  { "--7mcprob", eslARG_REAL,    "0.0", NULL, "x>-0.0001",NULL,    NULL,        NULL, "w/p7 banding set min cumulative prob to enter match state to <x>", 99 },
  { "--7iprob",  eslARG_REAL,    "1.0", NULL, "x<1.001",NULL,      NULL,        NULL, "w/p7 banding set max prob to enter insert state to <x>", 99 },
  { "--7ilprob", eslARG_REAL,    "1.0", NULL, "x<1.001",NULL,      NULL,        NULL, "w/p7 banding set max prob to enter left insert state to <x>", 99 },
  /* All options below are developer options, only shown if --devhelp invoked */
  /* Developer options related to alignment algorithm */
  { "--checkpost",eslARG_NONE,  FALSE, NULL, NULL,      NULL,      NULL, "--no-prob", "check that posteriors are correctly calc'ed", 101 },
  { "--no-null3",eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "turn OFF the NULL3 post hoc additional null model", 101 },
  /* developer options related to banded alignment */
  { "--checkfb", eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--hbanded",       "-l", "check that HMM posteriors for bands were correctly calc'ed", 102},
  { "--sums",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--hbanded",       NULL, "use posterior sums during HMM band calculation (widens bands)", 102 },
  { "--qdb",     eslARG_NONE,   FALSE, NULL, NULL,"--hbanded","--no-prob,--cyk",ACCOPTS, "use query dependent banded CYK alignment algorithm", 102 },
  { "--beta",    eslARG_REAL,   "1E-7",NULL, "0<x<1",   NULL,   "--qdb",        NULL, "set tail loss prob for --qdb to <x>", 102 },
  { "--hsafe",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--hbanded,--no-prob","--viterbi,--optacc", "realign (w/o bands) seqs with HMM banded CYK score < 0 bits", 102 },
  /* developer options related to output files and debugging */
  { "--regress", eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "save regression test data to file <f>", 103 },
  { "--bigmem",  eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "turn off memory saving mode, keep all seqs/parsetrees in memory", 103 },
  { "--banddump",eslARG_INT,    "0",   NULL, "0<=n<=3", NULL,      NULL,        NULL, "set verbosity of band info print statements to <n>", 103 },
  { "--dlev",    eslARG_INT,    "0",   NULL, "0<=n<=3", NULL,      NULL,        NULL, "set verbosity of debugging print statements to <n>", 103 },
  { "--stall",   eslARG_NONE,  FALSE, NULL, NULL,       NULL,      NULL,        NULL, "arrest after start: for debugging MPI under gdb", 103 },  
  /* Developer options related to experiment local begin/end modes */
  { "--pebegin", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      "-l",  "--pbegin", "set all local begins as equiprobable", 104 },
  { "--pfend",   eslARG_REAL,   NULL,  NULL, "0<x<1",   NULL,      "-l",    "--pend", "set all local end probs to <x>", 104 },
  { "--pbegin",  eslARG_REAL,  "0.05",NULL,  "0<x<1",   NULL,      "-l",        NULL, "set aggregate local begin prob to <x>", 104 },
  { "--pend",    eslARG_REAL,  "0.05",NULL,  "0<x<1",   NULL,      "-l",        NULL, "set aggregate local end prob to <x>", 104 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

/* struct cfg_s : "Global" application configuration shared by all threads/processes
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 */
struct cfg_s {
  char         *cmfile;	        /* name of input CM file  */ 
  char         *sqfile;	        /* name of sequence file  */ 
  ESL_SQFILE   *sqfp;           /* open sequence input file stream */
  int           fmt;		/* format code for seqfile */
  ESL_ALPHABET *abc;		/* digital alphabet for the CM */
  int           ncm;            /* number cm we're on */
  int           do_mpi;		/* TRUE if we're doing MPI parallelization */
  int           nproc;		/* how many MPI processes, total */
  int           my_rank;	/* who am I, in 0..nproc-1 */
  int           do_stall;	/* TRUE to stall the program until gdb attaches */
  ESL_RANDOMNESS *r;            /* source of randomness, only created if --sample enabled */
  int           nali;           /* number temporary alignment we're on */
  int           nseq;           /* number of sequences we've aligned thus far */

  /* Masters only (i/o streams) */
  CMFILE       *cmfp;		/* open input CM file stream       */
  FILE         *tmpfp;		/* the temporary output file where alignments are initially written */
  FILE         *ofp;		/* output file where alignments are ultimately written (default is stdout) */
  FILE         *tracefp;	/* optional output for parsetrees  */
  FILE         *insertfp;	/* optional output for insert info */
  FILE         *elfp;	        /* optional output for EL insert info */
  FILE         *scorefp;        /* optional output for alnment scores */
  FILE         *regressfp;	/* optional output for regression test  */
  ESL_MSAFILE  *withalifp;	/* optional input alignment to include */
  ESL_MSA      *withmsa;	/* MSA from withalifp to include */
  char         *withss_cons;	/* ss_cons string from withmsa (before knot stripping) */
  Parsetree_t  *withali_mtr;	/* guide tree for MSA from withalifp */
  ESL_ALPHABET *withali_abc;	/* digital alphabet for reading withali MSA */
  ESL_ALPHABET *abc_out;	/* digital alphabet for output */

  ESL_MSAFILE  *malifp;	        /* -M optional input alignment to include */
  ESL_MSA      **mali_msa;	/* MSAs from withalifp to include */
  Parsetree_t  **mali_mtr;	/* guide trees for MSAs from malifp */
  ESL_ALPHABET *mali_abc;	/* digital alphabet for reading mali MSA */
  int           mali_n;         /* number of alignments in mmsa */
};

static char usage[]  = "[-options] <cmfile> <sequence file>";
static char banner[] = "align sequences to an RNA CM";

static int  init_master_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
/* static int  init_shared_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf); */

static void  serial_master (const ESL_GETOPTS *go, struct cfg_s *cfg);
static void  serial_master_meta(const ESL_GETOPTS *go, struct cfg_s *cfg);

#ifdef HAVE_MPI
static int   mpi_master    (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int   mpi_worker    (const ESL_GETOPTS *go, struct cfg_s *cfg);
#endif

static int  process_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, seqs_to_aln_t *seqs_to_aln);
static int  output_result(const ESL_GETOPTS *go, struct cfg_s *cfg, int do_output_to_tmp, char *errbuf, CM_t *cm, seqs_to_aln_t *seqs_to_aln);

static int  initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int  check_withali(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, ESL_MSA **ret_msa, Parsetree_t **ret_mtr);
static int  include_withali(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, ESL_SQ ***ret_sq, Parsetree_t ***ret_tr, char ***ret_postcode, int *ret_nseq, char *errbuf);
static int  compare_cm_guide_trees(CM_t *cm1, CM_t *cm2);
static int  add_withali_pknots(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, ESL_MSA *newmsa);

static int  print_run_info(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf);
static void print_cm_info(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int  get_command(const ESL_GETOPTS *go, char *errbuf, char **ret_command);
static void print_info_file_header(FILE *fp, char *firstline, char *elstring);

/* Functions that enable memory efficiency by not storing all 
 * seqs/parsetrees from target file in memory at once.
 */
static int  create_and_output_final_msa(const ESL_GETOPTS *go, const struct cfg_s *cfg, FILE *ofp, char *errbuf, CM_t *cm, char *tmpfile);
static void update_maxins_and_maxel(ESL_MSA *msa, int clen, int64_t alen, int *maxins, int *maxel);
static int  determine_gap_columns_to_add(ESL_MSA *msa, int *maxins, int *maxel, int clen, int **ret_ngap_insA, int **ret_ngap_elA, int **ret_ngap_eitherA, char *errbuf);
static void inflate_gc_with_gaps_and_els(FILE *ofp, ESL_MSA *msa, int *ngap_insA, int *ngap_elA, char **ret_ss_cons2print, char **ret_rf2print);

/* meta-CM alignment functions, only used if -M enabled */
static int map_cpos_to_apos(ESL_MSA *msa, int **ret_c2a_map, int *ret_clen);
static int Parsetrees2Alignment_Minor2Major(CM_t *cm, const ESL_ALPHABET *abc, ESL_SQ **sq, float *wgt, 
					    Parsetree_t **tr, int nseq, int do_full, int do_matchonly, 
					    int *masteradd, ESL_MSA **ret_msa);
static int validate_meta(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t **cmlist, int ***ret_toadd2min, int **ret_maj_train_a2c_map);
static int major_alignment2minor_parsetrees(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t **cmlist, ESL_MSA *maj_target_msa, int *maj_train_a2c_map, int **ret_wcm);
static int majorfied_alignment2major_parsetrees(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *maj_cm, ESL_MSA *majed_min_msa, CMEmitMap_t *maj_emap, Parsetree_t ***ret_majed_tr);
/*
  static void print_stage_column_headings(const ESL_GETOPTS *go, const struct cfg_s *cfg);
  static int print_align_options(const struct cfg_s *cfg, CM_t *cm);
*/

#ifdef HAVE_MPI
static int determine_nseq_per_worker(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, int *ret_nseq_worker);
static int add_worker_seqs_to_master(seqs_to_aln_t *master_seqs, seqs_to_aln_t *worker_seqs, int offset);
#endif

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go = NULL;   /* command line processing                     */
  ESL_STOPWATCH   *w  = esl_stopwatch_Create();
  if(w == NULL) cm_Fail("Memory error, stopwatch not created.\n");
  esl_stopwatch_Start(w);
  struct cfg_s     cfg;
  int              m;           /* counter for cleaning up */
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
      puts("\nwhere basic options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1=docgroup, 2 = indentation; 80=textwidth*/
      puts("\nalignment algorithm related options:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\nbanded dynamic programming acceleration options:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80);
      puts("\noutput options:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80);
      puts("\noptions for including a fixed alignment within output alignment:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80);
      puts("\nusing a single CM from a multi-CM file:");
      esl_opt_DisplayHelp(stdout, go, 10, 2, 80);
      puts("\nverbose output files and debugging:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80);
      puts("\nexperimental options for plan7 banding using HMMER3 code:");
      esl_opt_DisplayHelp(stdout, go, 99, 2, 80);
      puts("\nundocumented developer algorithm options:");
      esl_opt_DisplayHelp(stdout, go, 101, 2, 80);
      puts("\nundocumented developer banded alignment options:");
      esl_opt_DisplayHelp(stdout, go, 102, 2, 80);
      puts("\nundocumented developer verbose output/debugging options:");
      esl_opt_DisplayHelp(stdout, go, 103, 2, 80);
      puts("\nundocumented developer options for experimental local begin/end modes:");
      esl_opt_DisplayHelp(stdout, go, 104, 2, 80);
      exit(0);
    }
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      cm_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\nwhere basic options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1=docgroup, 2 = indentation; 80=textwidth*/
      puts("\nalignment algorithm related options:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\nbanded dynamic programming acceleration options:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80);
      puts("\noutput options:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80);
      puts("\noptions for including a fixed alignment within output alignment:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80);
      puts("\nusing a single CM from a multi-CM file:");
      esl_opt_DisplayHelp(stdout, go, 10, 2, 80);
      puts("\nverbose output files and debugging:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80);
      exit(0);
    }
  if(esl_opt_ArgNumber(go) != 2) { 
      puts("Incorrect number of command line arguments.");
      esl_usage(stdout, argv[0], usage);
      puts("\nwhere basic options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      printf("\nTo see more help on other available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  /* Check for incompatible option combinations I don't know how to disallow with esl_getopts */
  /* --small requires EITHER --nonbanded or --qdb */
  if ((esl_opt_GetBoolean(go, "--small")) && (! ((esl_opt_GetBoolean(go, "--nonbanded")) || (esl_opt_GetBoolean(go, "--qdb"))))) { 
    esl_fatal("Error parsing options, --small is only allowed in combination with --nonbanded or --qdb.\n");
  }
  
  /* Initialize what we can in the config structure (without knowing the input alphabet yet).
   */
  cfg.cmfile     = esl_opt_GetArg(go, 1); 
  cfg.sqfile     = esl_opt_GetArg(go, 2); 
  cfg.sqfp       = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  if   (esl_opt_IsOn(go, "--informat")) { 
    cfg.fmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--informat"));
    if(cfg.fmt == eslSQFILE_UNKNOWN) cm_Fail("Can't recognize sequence file format: %s. valid options are: fasta, embl, genbank, ddbj, uniprot, stockholm, or pfam\n", esl_opt_GetString(go, "--informat"));
  } else
    cfg.fmt = eslSQFILE_UNKNOWN; /* autodetect sequence file format by default. */ 
  cfg.abc        = NULL;	           /* created in init_master_cfg() in masters, or in mpi_worker() in workers */
  if   (esl_opt_GetBoolean(go, "--dna")) cfg.abc_out = esl_alphabet_Create(eslDNA);
  else cfg.abc_out = esl_alphabet_Create(eslRNA); /* RNA. As it should be. */

  cfg.cmfp       = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.tmpfp      = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.ofp        = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.tracefp    = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.insertfp   = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.elfp       = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.scorefp    = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.regressfp  = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.withalifp  = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.withmsa    = NULL;	           /* filled in init_master_cfg() in masters, stays NULL for workers */
  cfg.withss_cons= NULL;	           /* filled in check_withali() in masters, stays NULL for workers */
  cfg.withali_mtr= NULL;	           /* filled in init_master_cfg() in masters, stays NULL for workers */
  cfg.withali_abc= NULL;	           /* created in init_master_cfg() in masters, stays NULL for workers */
  cfg.malifp     = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */

  cfg.mali_msa   = NULL;	           /* filled in init_master_cfg() in masters, stays NULL for workers */
  cfg.mali_mtr   = NULL;	           /* filled in init_master_cfg() in masters, stays NULL for workers */
  cfg.mali_abc   = NULL;	           /* created in init_master_cfg() in masters, stays NULL for workers */
  cfg.mali_n     = 0;	                   /* filled in serial_master_meta() if nec */

  cfg.ncm        = 0;
  cfg.nali       = 0;
  cfg.nseq       = 0;
  cfg.r          = NULL;	           /* created in init_master_cfg() for masters, mpi_worker() for workers*/

  cfg.do_mpi     = FALSE;	           /* this gets reset below, if we init MPI */
  cfg.nproc      = 0;		           /* this gets reset below, if we init MPI */
  cfg.my_rank    = 0;		           /* this gets reset below, if we init MPI */
  cfg.do_stall   = esl_opt_GetBoolean(go, "--stall");


  /* This is our stall point, if we need to wait until we get a
   * debugger attached to this process for debugging (especially
   * useful for MPI):
   */
  while (cfg.do_stall); 

  /* Figure out who we are, and send control there: 
   * we might be an MPI master, an MPI worker, or a serial program.
   */
#ifdef HAVE_MPI
  if (esl_opt_GetBoolean(go, "--mpi")) 
    {
      int              status;               /* easel status */
      char             errbuf[cmERRBUFSIZE]; /* for error messages in mpi_master() */
      cfg.do_mpi     = TRUE;

      MPI_Init(&argc, &argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &(cfg.my_rank));
      MPI_Comm_size(MPI_COMM_WORLD, &(cfg.nproc));

      if(cfg.nproc == 1) cm_Fail("MPI mode, but only 1 processor running...");

      if (cfg.my_rank > 0) { status = mpi_worker(go, &cfg); }
      else { 
	if(! esl_opt_GetBoolean(go, "-q")) cm_banner(stdout, argv[0], banner);
	status = mpi_master(go, &cfg, errbuf);
      }
      /* check status, if eslOK, we continue, else we exit. either way we call MPI_Finalize() */
      if(status == eslOK) { 
	esl_stopwatch_Stop(w);
	esl_stopwatch_MPIReduce(w, 0, MPI_COMM_WORLD);
	MPI_Finalize();
      }
      else { /* status != eslOK, master has error message in errbuf, worker does not */
	MPI_Finalize();
	if(cfg.my_rank == 0) cm_Fail(errbuf); /* master */
	else                 return 0;        /* worker */
      }
    }
  else
#endif /*HAVE_MPI*/
    {
      if(! esl_opt_GetBoolean(go, "-q")) cm_banner(stdout, argv[0], banner);
      if( esl_opt_IsOn(go, "-M"))        serial_master_meta(go, &cfg);
      serial_master(go, &cfg);
      esl_stopwatch_Stop(w);
    }
  /* Clean up the shared cfg. 
   */
  if (cfg.my_rank == 0) {
    if(cfg.tmpfp != NULL) { 
      fclose(cfg.tmpfp); 
    }
    if ( esl_opt_IsOn(go, "-o")) { 
      printf("# Alignment saved in file %s.\n", esl_opt_GetString(go, "-o"));
      fclose(cfg.ofp); 
    }
    if (cfg.tracefp   != NULL) { 
      printf("# Parsetrees saved in file %s.\n", esl_opt_GetString(go, "--tfile"));
      fclose(cfg.tracefp);
    }
    if (cfg.insertfp   != NULL) { 
      printf("# Per-sequence insert information saved in file %s.\n", esl_opt_GetString(go, "--ifile"));
      fclose(cfg.insertfp);
    }
    if (cfg.elfp   != NULL) { 
      printf("# Per-sequence EL insert information saved in file %s.\n", esl_opt_GetString(go, "--elfile"));
      fclose(cfg.elfp);
    }
    if (cfg.scorefp   != NULL) { 
      printf("# Alignment scores saved in file %s.\n", esl_opt_GetString(go, "--sfile"));
      fclose(cfg.scorefp);
    }
    if (cfg.regressfp   != NULL) {
      printf("# Regression data (alignment) saved in file %s.\n", esl_opt_GetString(go, "--regress"));
      fclose(cfg.regressfp);
    }
    if (cfg.cmfp      != NULL) CMFileClose(cfg.cmfp);
    if (cfg.sqfp      != NULL) esl_sqfile_Close(cfg.sqfp);
    if (cfg.withalifp != NULL) esl_msafile_Close(cfg.withalifp);
    if (cfg.withmsa   != NULL) esl_msa_Destroy(cfg.withmsa);
    if (cfg.withali_mtr != NULL) FreeParsetree(cfg.withali_mtr);
    if (cfg.withss_cons != NULL) free(cfg.withss_cons);
    if (cfg.malifp != NULL) esl_msafile_Close(cfg.malifp);
    if (cfg.mali_msa  != NULL) { for(m = 0; m < cfg.mali_n; m++) { esl_msa_Destroy(cfg.mali_msa[m]); } free(cfg.mali_msa); } 
    if (cfg.mali_mtr  != NULL) { for(m = 0; m < cfg.mali_n; m++) { FreeParsetree(cfg.mali_mtr[m]); }   free(cfg.mali_mtr); } 
  }
  if (cfg.r         != NULL) esl_randomness_Destroy(cfg.r);
  if (cfg.abc       != NULL) esl_alphabet_Destroy(cfg.abc);
  if (cfg.abc_out   != NULL) esl_alphabet_Destroy(cfg.abc_out);
  if (cfg.withali_abc != NULL) esl_alphabet_Destroy(cfg.withali_abc);
  if (cfg.my_rank == 0 && (! esl_opt_GetBoolean(go, "-q"))) { 
    printf("#\n");
    esl_stopwatch_Display(stdout, w, "# CPU time: ");
  }
  esl_getopts_Destroy(go);
  esl_stopwatch_Destroy(w);
  return 0;
}

/* init_master_cfg()
 * Called by masters, mpi or serial.
 * Already set:
 *    cfg->cmfile      - command line arg 1
 *    cfg->sqfile      - command line arg 2
 *    cfg->fmt         - format of output file
 * Sets: 
 *    cfg->sqfp        - open sequence file                
 *    cfg->tmpfp       - initial output file (temporary file)
 *    cfg->ofp         - ultimate output file (stdout by default)
 *    cfg->cmfp        - open CM file                
 *    cfg->abc         - digital input alphabet
 *    cfg->tracefp     - optional output file
 *    cfg->insertfp    - optional output file
 *    cfg->elfp        - optional output file
 *    cfg->scorefp     - optional output file
 *    cfg->regressfp   - optional output file
 *    cfg->withalifp   - optional input alignment file to include
 *    cfg->withmsa     - MSA from --withali file 
 *    cfg->withali_mtr - guide tree for MSA from --withali file 
 *    cfg->withali_abc - digital input alphabet for --withali file
 *    cfg->r           - source of randomness
 *                   
 * Errors in the MPI master here are considered to be "recoverable",
 * in the sense that we'll try to delay output of the error message
 * until we've cleanly shut down the worker processes. Therefore
 * errors return (code, errbuf) by the ESL_FAIL mech.
 */
static int
init_master_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  int  status;
  int  type;

  /* open input sequence file */
  status = esl_sqfile_Open(cfg->sqfile, cfg->fmt, NULL, &(cfg->sqfp));
  if (status == eslENOTFOUND)    ESL_FAIL(status, errbuf, "File %s doesn't exist or is not readable\n", cfg->sqfile);
  else if (status == eslEFORMAT) ESL_FAIL(status, errbuf, "Couldn't determine format of sequence file %s\n", cfg->sqfile);
  else if (status == eslEINVAL)  ESL_FAIL(status, errbuf, "Can’t autodetect stdin or .gz."); 
  else if (status != eslOK)      ESL_FAIL(status, errbuf, "Sequence file open failed with error %d\n", status);
  cfg->fmt = cfg->sqfp->format;

  /* Set the sqfile alphabet as RNA, if it's DNA we're fine. 
   * If it's not RNA nor DNA, we can't deal with it anyway,
   * so we're hardcoded to RNA.
   */
  cfg->abc = esl_alphabet_Create(eslRNA);
  if(cfg->abc == NULL) ESL_FAIL(status, errbuf, "Failed to create alphabet for sequence file");
  esl_sqfile_SetDigital(cfg->sqfp, cfg->abc);

  /* open CM file */
  if((cfg->cmfp = CMFileOpen(cfg->cmfile, NULL)) == NULL)
   ESL_FAIL(eslFAIL, errbuf, "Failed to open covariance model save file %s\n", cfg->cmfile);

  /* note, we don't open temporary output file, where we store intermediate alignments yet,
   * we do that in serial_master and mpi_master, because it needs to be rewritten for each
   * CM we align to (b/c we can align to more than one CM, though its probably uncommon) */

  /* open output file */
  if (esl_opt_GetString(go, "-o") != NULL) {
    if ((cfg->ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open -o output file %s\n", esl_opt_GetString(go, "-o"));
    } else cfg->ofp = stdout;

  /* seed master's RNG, this will only be used if --sample enabled, but we always initialize it for convenience (seeds always get sent to workers) */
  cfg->r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  /* optionally, open trace file */
  if (esl_opt_GetString(go, "--tfile") != NULL) {
    if ((cfg->tracefp = fopen(esl_opt_GetString(go, "--tfile"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --tfile output file %s\n", esl_opt_GetString(go, "--tfile"));
    }

  /* optionally, open insert info file */
  if (esl_opt_GetString(go, "--ifile") != NULL) {
    if ((cfg->insertfp = fopen(esl_opt_GetString(go, "--ifile"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --ifile output file %s\n", esl_opt_GetString(go, "--ifile"));
    print_info_file_header(cfg->insertfp, "Insert information file created by cmalign.", "");
  }

  /* optionally, open EL insert info file */
  if (esl_opt_GetString(go, "--elfile") != NULL) {
    if ((cfg->elfp = fopen(esl_opt_GetString(go, "--elfile"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --elfile output file %s\n", esl_opt_GetString(go, "--elfile"));
    print_info_file_header(cfg->elfp, "EL state (local end) insert information file created by cmalign.", "EL ");
  }

  /* optionally, open score info file */
  if (esl_opt_GetString(go, "--sfile") != NULL) {
    if ((cfg->scorefp = fopen(esl_opt_GetString(go, "--sfile"), "w")) == NULL) 
      ESL_FAIL(eslFAIL, errbuf, "Failed to open --sfile output file %s\n", esl_opt_GetString(go, "--file"));
  }

  /* optionally, open regression file */
  if (esl_opt_GetString(go, "--regress") != NULL) {
    if ((cfg->regressfp = fopen(esl_opt_GetString(go, "--regress"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --regress output file %s\n", esl_opt_GetString(go, "--regress"));
    }

  /* optionally, open withali file for reading */
  if(esl_opt_GetString(go, "--withali") != NULL)
    {
      status = esl_msafile_Open(esl_opt_GetString(go, "--withali"), eslMSAFILE_UNKNOWN, NULL, &(cfg->withalifp));
      if (status == eslENOTFOUND)    ESL_FAIL(status, errbuf, "--withali alignment file %s doesn't exist or is not readable\n", 
					      esl_opt_GetString(go, "--withali"));
      else if (status == eslEFORMAT) ESL_FAIL(status, errbuf, "Couldn't determine format of --withali alignment %s\n", 
					      esl_opt_GetString(go, "--withali"));
      else if (status != eslOK)      ESL_FAIL(status, errbuf, "Alignment file open failed with error %d\n", status);
      /* Guess the withali alphabet, if it's ambiguous, guess RNA,
       * we'll treat RNA and DNA both as RNA internally.
       * We can't handle any other alphabets, so this is hardcoded. */
      status = esl_msafile_GuessAlphabet(cfg->withalifp, &type);
      if (status == eslEAMBIGUOUS)    type = eslRNA; /* guess it's RNA, we'll fail downstream with an error message if it's not */
      else if (status == eslEFORMAT)  ESL_FAIL(status, errbuf, "Alignment file parse failed: %s\n", cfg->withalifp->errbuf);
      else if (status == eslENODATA)  ESL_FAIL(status, errbuf, "Alignment file %s is empty\n", esl_opt_GetString(go, "--withali"));
      else if (status != eslOK)       ESL_FAIL(status, errbuf, "Failed to read alignment file %s\n", esl_opt_GetString(go, "--withali"));
      /* we can read DNA/RNA but internally we treat it as RNA */
      if(! (type == eslRNA || type == eslDNA))
	ESL_FAIL(eslEFORMAT, errbuf, "Alphabet is not DNA/RNA in %s\n", esl_opt_GetString(go, "--withali"));
      cfg->withali_abc = esl_alphabet_Create(eslRNA);
      if(cfg->withali_abc == NULL) ESL_FAIL(status, errbuf, "Failed to create alphabet for --withali");
      esl_msafile_SetDigital(cfg->withalifp, cfg->withali_abc);
    }

  if(cfg->r == NULL) ESL_FAIL(eslEMEM, errbuf, "Failed to create master RNG.");
  return eslOK;
}

/* serial_master()
 * The serial version of cmalign.
 * Align each sequence to the CM.
 * 
 * A master can only return if it's successful. All errors are handled immediately and fatally with cm_Fail().
 */
static void
serial_master(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int      status;
  char     errbuf[cmERRBUFSIZE];
  CM_t     *cm;
  seqs_to_aln_t  *seqs_to_aln;       /* sequences to align, holds seqs, parsetrees, CP9 traces, postcodes */
  int       used_at_least_one_cm = FALSE; 
  int       chunksize;               /* number of seqs/parsetrees to hold in memory at once, size of temp alignments */
  int       do_small;                /* TRUE to try to save memory by outputting temp alignments, then merging */
  int       keep_reading = TRUE;
  int       do_output_to_tmp = TRUE; /* should we output to cfg->tmpfp (and then merge at end) or to cfg->ofp directly */
  int       created_tmpfile = FALSE; /* set to TRUE if we needed and opened the tmpfile */
  int       nseq_to_read;            /* number of sequences to read from sqfp for current chunk */
  int       do_read_all;             /* TRUE to read all seqs from current chunk, FALSE not to */

  chunksize = esl_opt_GetInteger(go, "--chunk");
  do_small  = ((esl_opt_GetBoolean(go, "--ileaved")) || (esl_opt_GetBoolean(go, "--bigmem"))) ? FALSE : TRUE;

  if ((status  = init_master_cfg(go, cfg, errbuf)) != eslOK)  cm_Fail(errbuf);
  if ((status  = print_run_info (go, cfg, errbuf))  != eslOK) cm_Fail(errbuf);
  
  while ((status = CMFileRead(cfg->cmfp, errbuf, &(cfg->abc), &cm)) == eslOK)
    {
      if (cm == NULL) cm_Fail("Failed to read CM from %s -- file corrupt?\n", cfg->cmfile);
      cfg->ncm++;
      cfg->nali = 0;

      if(! esl_opt_IsDefault(go, "--cm-idx")) { 
	if(cfg->ncm != esl_opt_GetInteger(go, "--cm-idx")) { FreeCM(cm); continue; }
      }	
      if(! esl_opt_IsDefault(go, "--cm-name")) { 
	if(strcmp(cm->name, esl_opt_GetString(go, "--cm-name")) != 0) { FreeCM(cm); continue; }
      }	
      used_at_least_one_cm = TRUE;

      char tmpfile[32] = "esltmpXXXXXX"; /* name of the tmpfile, needs to be here, so it's reset for each CM, else esl_tmpfile_named will fail on second CM */

      /* initialize the flags/options/params and configuration of the CM */
      if((status   = initialize_cm(go, cfg, errbuf, cm))                    != eslOK)    cm_Fail(errbuf);
      print_cm_info (go, cfg, errbuf, cm);
      
      /* read in chunksize seqs at a time, align each of them and output the result */
      keep_reading = TRUE;
      do_output_to_tmp = TRUE; /* until proven otherwise */
      while(keep_reading) { 
	seqs_to_aln = CreateSeqsToAln(chunksize, FALSE);

	/* determine the number of sequences to read from the input file: 
	 *    if   do_small                       : read <chunksize> seqs
	 *    else (--ileaved or --bigmem enabled): read all seqs
	 */
	nseq_to_read = (do_small) ? chunksize : 0; /* 0 is a flag for 'read all seqs' in ReadSeqsToAln() */

	status = ReadSeqsToAln(cfg->abc, cfg->sqfp, nseq_to_read, seqs_to_aln, 
			       FALSE); /* i'm not an mpi master */

	if(seqs_to_aln->nseq > 0) { 
	  /* determine if we're finished aligning all seqs we need to align (either all seqs or <n> from --eseq <n>) */ 
	  if(status == eslEOF) {
	    keep_reading = FALSE;
	    if(cfg->nali == 0) { 
	      /* If we get here, first alignment will be the full alignment b/c either n < chunksize 
	       * seqs in target file or --ileaved enabled, either way we don't write to tmpfile, 
	       * instead we output directly to final output file. 
	       */
	      do_output_to_tmp = FALSE; /* note, this only occurs if keep_reading was set to FALSE */
	    }
	  }
	  if (do_output_to_tmp && cfg->nali == 0) { 
	    /* first aln for temporary output file, open the file */	
	    if ((status = esl_tmpfile_named(tmpfile, &(cfg->tmpfp))) != eslOK) { 
	      cm_Fail("Failed to open temporary output file for cm %d (status %d)", cfg->ncm, status);
	    }
	    created_tmpfile = TRUE;
	  }
	  /* align current sequences */
	  if ((status = process_workunit(go, cfg, errbuf, cm, seqs_to_aln)) != eslOK) cm_Fail(errbuf);
	  if ((status = output_result   (go, cfg, do_output_to_tmp, errbuf, cm, seqs_to_aln)) != eslOK) cm_Fail(errbuf);
	  cfg->nali++;
	  cfg->nseq += seqs_to_aln->nseq;
	} /* end of if(seqs_to_aln->nseq > 0) */
	else { 
	  keep_reading = FALSE;
	  if(status != eslEOF) cm_Fail("Error reading sequence file."); 
	}
	FreeSeqsToAln(seqs_to_aln);
      } /* end of while(keep_reading) */

      if(do_output_to_tmp) { 
	fclose(cfg->tmpfp); /* we're done writing to tmpfp */
	cfg->tmpfp = NULL;
	/* merge all temporary alignments now in cfg->tmpfp, and output merged alignment */
	if((cfg->ofp == stdout) && (! esl_opt_GetBoolean(go, "-q"))) fprintf(cfg->ofp, "\n"); /* blank line to separate scores from aln if printing aln to stdout */
	if((status = create_and_output_final_msa(go, cfg, cfg->ofp, errbuf, cm, tmpfile)) != eslOK) cm_Fail(errbuf);
	/* if --regress, output alignment again, this time to ofp->regressfp */
	if(cfg->regressfp != NULL) { 
	  if((status = create_and_output_final_msa(go, cfg, cfg->regressfp, errbuf, cm, tmpfile)) != eslOK) cm_Fail(errbuf);
	}
      }
      /* else (! do_output_to_tmp): we outputted alignment directly to ofp (and possibly also to regressfp) in output_result()) */

      /* finish insert and el files */
      if(cfg->insertfp != NULL) { fprintf(cfg->insertfp, "//\n"); }
      if(cfg->elfp != NULL)     { fprintf(cfg->elfp,     "//\n"); }
	
      /* clean up */
      FreeCM(cm);
      if(created_tmpfile) remove(tmpfile); 
      esl_sqfile_Position(cfg->sqfp, (off_t) 0); /* we may be aligning the seqs in this file again with another CM */
    }
  if(status != eslEOF) cm_Fail(errbuf);
  if(! used_at_least_one_cm) { 
    if(! esl_opt_IsDefault(go, "--cm-idx"))  { cm_Fail("--cm-idx %d enabled, but only %d CMs in the cmfile.\n", esl_opt_GetInteger(go, "--cm-idx"), cfg->ncm); }
    if(! esl_opt_IsDefault(go, "--cm-name")) { cm_Fail("--cm-name %s enabled, but no CM named %s exists in the cmfile.\n", esl_opt_GetString(go, "--cm-name"), esl_opt_GetString(go, "--cm-name")); }
  }
  return;
}

#ifdef HAVE_MPI
/* mpi_master()
 * The MPI version of cmalign
 * Follows standard pattern for a master/worker load-balanced MPI program 
 * (SRE notes J1/78-79).
 * 
 * A master returns eslOK if it's successful. 
 * Errors in an MPI master come in two classes: recoverable and nonrecoverable.
 * If a recoverable error occurs, errbuf is filled with an error message
 * from the master or a worker, and it's sent back while returning a
 * non-eslOK error code.
 * 
 * Recoverable errors include (hopefully) all worker-side errors, and any
 * master-side error that do not affect MPI communication. Error
 * messages from recoverable messages are delayed until we've cleanly
 * shut down the workers.
 * 
 * Unrecoverable errors are master-side errors that may affect MPI
 * communication, meaning we cannot count on being able to reach the
 * workers and shut them down. Unrecoverable errors result in immediate
 * cm_Fail()'s, which will cause MPI to shut down the worker processes
 * uncleanly.
 */
static int
mpi_master(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  int      xstatus       = eslOK;	/* changes from OK on recoverable error */
  int      status;
  int      have_work     = TRUE;	/* TRUE while work remains  */
  int      nproc_working = 0;	        /* number of worker processes working, up to nproc-1 */
  int      wi;          	        /* rank of next worker to get an alignment to work on */
  char    *buf           = NULL;	/* input/output buffer, for packed MPI messages */
  int      bn            = 0;
  int      pos = 1;
  int      wi_error = 0;                /* worker index that sent back an error message, if an error occurs */
  CM_t    *cm;
  int used_at_least_one_cm = FALSE; 
  int nseq_per_worker  = 0;
  int nseq_this_worker = 0;
  int nseq_remaining   = 0;
  int nseq_sent        = 0;
  int do_small;                /* TRUE to try to save memory by outputting temp alignments, then merging */
  int chunksize;               /* number of seqs/parsetrees to hold in memory at once, size of temp alignments */
  int keep_reading = TRUE;
  int do_output_to_tmp = TRUE; /* should we output to cfg->tmpfp (and then merge at end) or to cfg->ofp directly */
  int  created_tmpfile = FALSE;      /* set to TRUE if we need and open the tmpfile */

  seqs_to_aln_t  *all_seqs_to_aln    = NULL;
  seqs_to_aln_t  *worker_seqs_to_aln = NULL;
  int            *seqidx         = NULL;
  uint32_t       *seedlist = NULL;       /* seeds for worker's RNGs, we send these to workers */
  
  MPI_Status mpistatus; 
  int      n;

  chunksize = esl_opt_GetInteger(go, "--chunk");
  do_small =  ((esl_opt_GetBoolean(go, "--ileaved")) || (esl_opt_GetBoolean(go, "--bigmem"))) ? FALSE : TRUE;

  /* Master initialization: including, figure out the alphabet type.
   * If any failure occurs, delay printing error message until we've shut down workers.
   */

  if (xstatus == eslOK) { if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) xstatus = status; }
  if (xstatus == eslOK) { bn = 4096; if ((buf = malloc(sizeof(char) * bn))             == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((seedlist       = malloc(sizeof(uint32_t) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((seqidx         = malloc(sizeof(int)  * cfg->nproc))     == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((status         = print_run_info(go, cfg, errbuf)) != eslOK) xstatus = status; }

  MPI_Bcast(&xstatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (xstatus != eslOK) return xstatus; /* errbuf was filled above */
  ESL_DPRINTF1(("MPI master is initialized\n"));

  for (wi = 0; wi < cfg->nproc; wi++) seqidx[wi] = 0;

  for (wi = 0; wi < cfg->nproc; wi++) {
    seedlist[wi] = esl_rnd_Roll(cfg->r, 1000000000); /* not sure what to use as max for seed */
    ESL_DPRINTF1(("wi %d seed: %ld\n", wi, seedlist[wi]));
  }

  /* Worker initialization:
   * Because we've already successfully initialized the master before we start
   * initializing the workers, we don't expect worker initialization to fail;
   * so we just receive a quick OK/error code reply from each worker to be sure,
   * and don't worry about an informative message. 
   */
  for (wi = 1; wi < cfg->nproc; wi++)
    MPI_Send(&(seedlist[wi]), 1, MPI_UNSIGNED_LONG, wi, 0, MPI_COMM_WORLD);
  MPI_Reduce(&xstatus, &status, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  if (status != eslOK) cm_Fail("One or more MPI worker processes failed to initialize.");
  ESL_DPRINTF1(("%d workers are initialized\n", cfg->nproc-1));
  free(seedlist);

  /* Main loop: combining load workers, send/receive, clear workers loops;
   * also, catch error states and die later, after clean shutdown of workers.
   * 
   * When a recoverable error occurs, have_work = FALSE, xstatus !=
   * eslOK, and errbuf is set to an informative message. No more
   * errbuf's can be received after the first one. We wait for all the
   * workers to clear their work units, then send them shutdown signals,
   * then finally print our errbuf and exit.
   * 
   * Unrecoverable errors just crash us out with cm_Fail().
   */

  while (xstatus == eslOK && ((status = CMFileRead(cfg->cmfp, errbuf, &(cfg->abc), &cm)) == eslOK))
    {
      if (cm == NULL) cm_Fail("Failed to read CM from %s -- file corrupt?\n", cfg->cmfile);
      cfg->ncm++;  
      cfg->nali = 0;
      ESL_DPRINTF1(("MPI master read CM number %d\n", cfg->ncm));

      if(! esl_opt_IsDefault(go, "--cm-idx")) { 
	if(cfg->ncm != esl_opt_GetInteger(go, "--cm-idx")) { FreeCM(cm); continue; }
      }	
      if(! esl_opt_IsDefault(go, "--cm-name")) { 
	if(strcmp(cm->name, esl_opt_GetString(go, "--cm-name")) != 0) { FreeCM(cm); continue; }
      }	
      used_at_least_one_cm = TRUE;

      char tmpfile[32] = "esltmpXXXXXX"; /* name of the tmpfile, needs to be here, so it's reset for each CM, else esl_tmpfile_named will fail on second CM */

      if((status = cm_master_MPIBcast(cm, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI broadcast CM failed.");
      
      /* initialize the flags/options/params of the CM */
      if((status   = initialize_cm(go, cfg, errbuf, cm))                    != eslOK) cm_Fail(errbuf);
      print_cm_info (go, cfg, errbuf, cm);
      determine_nseq_per_worker(go, cfg, cm, &nseq_per_worker); /* this func dies internally if there's some error */
      ESL_DPRINTF1(("nseq_per_worker: %d\n", nseq_per_worker));

      /* Read in chunksize seqs at a time from target file.
       * For each chunk, we send <nseq_per_worker> seqs to each worker,
       * once a chunk is completed by all workers, we create the alignment 
       * of the chunk and output it to a temp file.
       * A special case:
       * if (!do_small) (either --ileaved or --bigmem enabled), read all seqs (first chunk is only chunk, contains all seqs)
       */
      keep_reading     = TRUE;
      do_output_to_tmp = TRUE; /* until proven otherwise */
      while(keep_reading) { 
	have_work = TRUE;
	wi = 1;
	all_seqs_to_aln = CreateSeqsToAln(chunksize, TRUE);

	/* determine the number of sequences to read from the input file: 
	 *    if   do_small                       : read <chunksize> seqs
	 *    else (--ileaved or --bigmem enabled): read all seqs
	 */
	nseq_to_read = (do_small) ? chunksize : 0; /* 0 is a flag for 'read all seqs' in ReadSeqsToAln() */

	status = ReadSeqsToAln(cfg->abc, cfg->sqfp, nseq_to_read, all_seqs_to_aln,
			       TRUE); /* i am an mpi master */

	if(all_seqs_to_aln->nseq > 0) { 
	  nseq_remaining = all_seqs_to_aln->nseq;
	  nseq_sent = 0;

	  if(status == eslEOF) { 
	    keep_reading = FALSE;
	    if(cfg->nali == 0) { 
	      /* If we get here, first alignment will be the full alignment b/c either n < chunksize 
	       * seqs in target file or --ileaved enabled, either way we don't write to tmpfile, 
	       * instead we output directly to final output file. 
	       */
	      do_output_to_tmp = FALSE; /* note, this only occurs if keep_reading was set to FALSE */
	    }
	  }
	  if (do_output_to_tmp && cfg->nali == 0) { 
	    /* first aln for temporary output file, open the file */
	    if ((status = esl_tmpfile_named(tmpfile, &(cfg->tmpfp))) != eslOK) { 
	      cm_Fail("Failed to open temporary output file for cm %d (status %d)", cfg->ncm, status);
	    }
	  }
	  
	  /* Main send/receive loop: where we send <nseq_per_worker> sequences to workers and receive parsetrees back. 
	   * This loop is entered once for each chunk of seqs we read in the ReadSeqsToAln() call above. 
	   */
	  while (have_work || nproc_working)
	    {
	      if (have_work) 
		{
		  if(nseq_sent == all_seqs_to_aln->nseq) { 
		    have_work = FALSE;
		    ESL_DPRINTF1(("MPI master has sent all %d of its sequences (%d)\n", all_seqs_to_aln->nseq));
		  }
		  else { 
		    nseq_this_worker = ESL_MAX(nseq_per_worker, (all_seqs_to_aln->nseq - nseq_sent));
		  }
		}
	      
	      if ((have_work && nproc_working == cfg->nproc-1) || (!have_work && nproc_working > 0))
		{
		  if (MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &mpistatus) != 0) cm_Fail("mpi probe failed");
		  if (MPI_Get_count(&mpistatus, MPI_PACKED, &n)                != 0) cm_Fail("mpi get count failed");
		  wi = mpistatus.MPI_SOURCE;
		  ESL_DPRINTF1(("MPI master sees a result of %d bytes from worker %d\n", n, wi));
		  
		  if (n > bn) {
		    if ((buf = realloc(buf, sizeof(char) * n)) == NULL) cm_Fail("reallocation failed");
		    bn = n; 
		  }
		  if (MPI_Recv(buf, bn, MPI_PACKED, wi, 0, MPI_COMM_WORLD, &mpistatus) != 0) cm_Fail("mpi recv failed");
		  
		  /* If we're in a recoverable error state, we're only clearing worker results;
		   * just receive them, don't unpack them or print them.
		   * But if our xstatus is OK, go ahead and process the result buffer.
		   */
		  if (xstatus == eslOK) /* worker reported success. Get the result. */
		    {
		      pos = 0;
		      if (MPI_Unpack(buf, bn, &pos, &xstatus, 1, MPI_INT, MPI_COMM_WORLD)     != 0)     cm_Fail("mpi unpack failed");
		      if (xstatus == eslOK) /* worker reported success. Get the results. */
			{
			  ESL_DPRINTF1(("MPI master sees that the result buffer contains aligned sequences (seqidx: %d)\n", seqidx[wi]));
			  if ((status = cm_seqs_to_aln_MPIUnpack(cfg->abc, buf, bn, &pos, MPI_COMM_WORLD, &worker_seqs_to_aln)) != eslOK) cm_Fail("search results unpack failed");
			  ESL_DPRINTF1(("MPI master has unpacked search results\n"));
			  if ((status = add_worker_seqs_to_master(all_seqs_to_aln, worker_seqs_to_aln, seqidx[wi])) != eslOK) cm_Fail("adding worker results to master results failed");
			  /* careful not to free data from worker_seqs_to_aln we've
			   * just added to all_seqs_to_aln. we didn't copy it, we just
			   * had pointers in all_seqs_to_aln point to it. We can 
			   * free the worker's pointers to those pointers though */
			  if(worker_seqs_to_aln->sq       != NULL) free(worker_seqs_to_aln->sq);
			  if(worker_seqs_to_aln->tr       != NULL) free(worker_seqs_to_aln->tr);
			  if(worker_seqs_to_aln->cp9_tr   != NULL) free(worker_seqs_to_aln->cp9_tr);
			  if(worker_seqs_to_aln->postcode != NULL) free(worker_seqs_to_aln->postcode);
			  if(worker_seqs_to_aln->sc       != NULL) free(worker_seqs_to_aln->sc);
			  if(worker_seqs_to_aln->pp       != NULL) free(worker_seqs_to_aln->pp);
			  if(worker_seqs_to_aln->struct_sc!= NULL) free(worker_seqs_to_aln->struct_sc);
			  free(worker_seqs_to_aln);
			}
		      else	/* worker reported an error. Get the errbuf. */
			{
			  if (MPI_Unpack(buf, bn, &pos, errbuf, cmERRBUFSIZE, MPI_CHAR, MPI_COMM_WORLD) != 0) cm_Fail("mpi unpack of errbuf failed");
			  ESL_DPRINTF1(("MPI master sees that the result buffer contains an error message\n"));
			  have_work = FALSE;
			  wi_error  = wi;
			}
		    }
		  nproc_working--;
		}
	      
	      if (have_work)
		{   
		  /* send new alignment job */
		  ESL_DPRINTF1(("MPI master is sending sequence to search to worker %d\n", wi));
		  if ((status = cm_seqs_to_aln_MPISend(all_seqs_to_aln, nseq_sent, nseq_this_worker, wi, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI search job send failed");
		  seqidx[wi] = nseq_sent;
		  nseq_sent += nseq_this_worker;
		  wi++;
		  nproc_working++;
		}
	    }
	  /* if we've got valid results, output them */
	  if (xstatus == eslOK) { 
	    if ((status = output_result(go, cfg, do_output_to_tmp, errbuf, cm, all_seqs_to_aln)) != eslOK) cm_Fail(errbuf);
	  }
	  cfg->nali++;
	  cfg->nseq += all_seqs_to_aln->nseq;
	} /* end of if(all_seqs_to_aln->nseq > 0) */
	else { 
	  keep_reading = FALSE;
	  if(status != eslEOF) cm_Fail("Error reading sequence file."); 
	}
	FreeSeqsToAln(all_seqs_to_aln);
      } /* end of while(keep_reading) */


      ESL_DPRINTF1(("MPI master: done with this CM. Telling all workers\n"));
      /* send workers the message that we're done with this CM */
      for (wi = 1; wi < cfg->nproc; wi++) 
	if ((status = cm_seqs_to_aln_MPISend(NULL, 0, 0, wi, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("Shutting down a worker failed.");

      if(do_output_to_tmp) { 
	fclose(cfg->tmpfp); /* we're done writing to tmpfp */
	cfg->tmpfp = NULL;
	/* merge all temporary alignments now in cfg->tmpfp, and output merged alignment */
	if((cfg->ofp == stdout) && (! esl_opt_GetBoolean(go, "-q"))) fprintf(cfg->ofp, "\n"); /* blank line to separate scores from aln if printing aln to stdout */
	if((status = create_and_output_final_msa(go, cfg, cfg->ofp, errbuf, cm, tmpfile)) != eslOK) cm_Fail(errbuf);
	/* if --regress, output alignment again, this time to ofp->regressfp */
	if(cfg->regressfp != NULL) { 
	  if((status = create_and_output_final_msa(go, cfg, cfg->regressfp, errbuf, cm, tmpfile)) != eslOK) cm_Fail(errbuf);
	}
      }
      /* else (! do_output_to_tmp): we outputted alignment directly to ofp (and possibly also to regressfp) in output_result()) */

      /* finish insert and el files */
      if(cfg->insertfp != NULL) { fprintf(cfg->insertfp, "//\n"); }
      if(cfg->elfp != NULL)     { fprintf(cfg->elfp,     "//\n"); }

      FreeCM(cm);
      if(created_tmpfile) remove(tmpfile);
      esl_sqfile_Position(cfg->sqfp, (off_t) 0); /* we may be aligning this file again with another CM */
    }

  /* On success or recoverable errors:
   * Shut down workers cleanly. 
   */
  ESL_DPRINTF1(("MPI master is done. Shutting down all the workers cleanly\n"));
  if((cm_master_MPIBcast(NULL, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI broadcast CM failed.");
  free(buf);
  
  if     (xstatus != eslOK) { fprintf(stderr, "Worker: %d had a problem.\n", wi_error); return xstatus; }
  else if((! used_at_least_one_cm) && (! esl_opt_IsDefault(go, "--cm-idx")))  { ESL_FAIL(eslEINVAL, errbuf, "--cm-idx %d enabled, but only %d CMs in the cmfile.\n", esl_opt_GetInteger(go, "--cm-idx"), cfg->ncm); }
  else if((! used_at_least_one_cm) && (! esl_opt_IsDefault(go, "--cm-name"))) { ESL_FAIL(eslEINVAL, errbuf, "--cm-name %s enabled, but no CM named %s exists in the cmfile.\n", esl_opt_GetString(go, "--cm-name"), esl_opt_GetString(go, "--cm-name")); }
  else if(status != eslEOF) return status;
  else                      return eslOK; 
}

static int
mpi_worker(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int           xstatus = eslOK;
  int           status;
  CM_t         *cm  = NULL;
  char         *wbuf = NULL;	/* packed send/recv buffer  */
  int           wn   = 0;	/* allocation size for wbuf */
  int           sz, n;		/* size of a packed message */
  int           pos;
  char          errbuf[cmERRBUFSIZE];
  seqs_to_aln_t *seqs_to_aln = NULL;
  uint32_t       seed;                  /* seed for RNG, rec'd from master */
  int           i;
  MPI_Status mpistatus;           /* MPI status... */

  /* After master initialization: master broadcasts its status.
   */
  MPI_Bcast(&xstatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (xstatus != eslOK) return xstatus; /* master saw an error code; workers do an immediate normal shutdown. */
  ESL_DPRINTF1(("worker %d: sees that master has initialized\n", cfg->my_rank));
  
  /* Master now sends worker initialization information (RNG seed) 
   * Workers returns their status post-initialization.
   * Initial allocation of wbuf must be large enough to guarantee that
   * we can pack an error result into it, because after initialization,
   * errors will be returned as packed (code, errbuf) messages.
   */
  if (MPI_Recv(&seed, 1, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD, &mpistatus) != 0) ESL_XEXCEPTION(eslESYS, "mpi recv failed");
  if (xstatus == eslOK) { if((cfg->r = esl_randomness_Create(seed)) == NULL)          xstatus = eslEMEM; }
  if (xstatus == eslOK) { wn = 4096;  if ((wbuf = malloc(wn * sizeof(char))) == NULL) xstatus = eslEMEM; }
  MPI_Reduce(&xstatus, &status, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD); /* everyone sends xstatus back to master */
  if (xstatus != eslOK) {
    if (wbuf != NULL) free(wbuf);
    return xstatus; /* shutdown; we passed the error back for the master to deal with. */
  }
  ESL_DPRINTF1(("worker %d: initialized seed: %ld\n", cfg->my_rank, seed));

  /* source = 0 (master); tag = 0 */
  while ((status = cm_worker_MPIBcast(0, MPI_COMM_WORLD, &wbuf, &wn, &(cfg->abc), &cm)) == eslOK)
    {
      ESL_DPRINTF1(("Worker %d succesfully received CM, num states: %d num nodes: %d\n", cfg->my_rank, cm->M, cm->nodes));
      
      /* initialize the flags/options/params of the CM */
      if((status   = initialize_cm(go, cfg, errbuf, cm))                    != eslOK) goto ERROR;
      
      while((status = cm_seqs_to_aln_MPIRecv(cfg->abc, 0, 0, MPI_COMM_WORLD, &wbuf, &wn, &seqs_to_aln)) == eslOK)
	{
	  ESL_DPRINTF1(("worker %d: has received alignment job, nseq: %d\n", cfg->my_rank, seqs_to_aln->nseq));
	  /* align all sequences */
	  if ((status = process_workunit(go, cfg, errbuf, cm, seqs_to_aln)) != eslOK) goto ERROR;
	  ESL_DPRINTF1(("worker %d: has gathered alignment results\n", cfg->my_rank));

	  /* free the sequences, master already has a copy so we don't send them back */
	  for(i = 0; i < seqs_to_aln->nseq; i++) esl_sq_Destroy(seqs_to_aln->sq[i]);
	  free(seqs_to_aln->sq);
	  seqs_to_aln->sq = NULL;
	  n = 0;
	  if (MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &sz) != 0) /* room for the status code */
	    ESL_XFAIL(eslESYS, errbuf, "mpi pack size failed"); 
	  n += sz;
	  if (cm_seqs_to_aln_MPIPackSize(seqs_to_aln, 0, seqs_to_aln->nseq, MPI_COMM_WORLD, &sz) != eslOK) 
	    ESL_XFAIL(eslFAIL, errbuf, "cm_seqs_to_aln_MPIPackSize() call failed"); 
	  n += sz;  
	  if (n > wn) {
	    void *tmp;
	    ESL_RALLOC(wbuf, tmp, sizeof(char) * n);
	    wn = n;
	  }
	  ESL_DPRINTF1(("worker %d: has calculated the alignment results will pack into %d bytes\n", cfg->my_rank, n));
	  status = eslOK;

	  pos = 0;
	  if (MPI_Pack(&status, 1, MPI_INT, wbuf, wn, &pos, MPI_COMM_WORLD) != 0) 
	    ESL_XFAIL(eslESYS, errbuf, "mpi pack failed.");
	  if (cm_seqs_to_aln_MPIPack(seqs_to_aln, 0, seqs_to_aln->nseq, wbuf, wn, &pos, MPI_COMM_WORLD) != eslOK) 
	    ESL_XFAIL(eslFAIL, errbuf, "cm_seqs_to_aln_MPIPack() call failed"); 
	  MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);
	  ESL_DPRINTF1(("worker %d: has sent results to master in message of %d bytes\n", cfg->my_rank, pos));
	  
	  FreeSeqsToAln(seqs_to_aln);
	}
      if(status == eslEOD) ESL_DPRINTF1(("worker %d: has seen message to stop with this CM.\n", cfg->my_rank));
      else ESL_XFAIL(eslFAIL, errbuf, "within CM loop, unexpected status code: %d received from cm_seqs_to_aln_MPIRecv()\n", status);

      FreeCM(cm);
      cm = NULL;
    }
  if (status == eslEOD) ESL_DPRINTF1(("worker %d told CMs are done.\n", cfg->my_rank));
  else ESL_FAIL(eslFAIL, errbuf, "outside CM loop, unexpected status code: %d received from cm_seqs_to_aln_MPIRecv()\n", status);

  if (wbuf != NULL) free(wbuf);
  return eslOK;

 ERROR:

  ESL_DPRINTF1(("worker %d: fails, is sending an error message, as follows:\n%s\n", cfg->my_rank, errbuf));
  pos = 0;
  MPI_Pack(&status, 1,               MPI_INT,  wbuf, wn, &pos, MPI_COMM_WORLD);
  MPI_Pack(errbuf,  cmERRBUFSIZE,    MPI_CHAR, wbuf, wn, &pos, MPI_COMM_WORLD);
  MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);

  /* if we get here this worker failed and sent an error message, now the master knows a worker
   * failed but it has to send the message to all other workers (besides this one) to abort so they 
   * can be shut down cleanly. As currently implemented, this means we have to wait here for that 
   * signal which comes in the form of a special 'empty' work packet that tells us we're done with
   * the current CM, and then a 'empty' CM broadcast that tells us we're done with all CMs in the file.
   */
  status = cm_seqs_to_aln_MPIRecv(cfg->abc, 0, 0, MPI_COMM_WORLD, &wbuf, &wn, &seqs_to_aln);
  status = cm_worker_MPIBcast(0, MPI_COMM_WORLD, &wbuf, &wn, &(cfg->abc), &cm);
  /* status after each of the above calls should be eslEOD, but if it isn't we can't really do anything 
   * about it b/c we've already sent our error message, so in that scenario the MPI will break uncleanly 
   */
  return eslFAIL; /* recoverable error, master has error message and will print it */
}
#endif /*HAVE_MPI*/

static int
output_result(const ESL_GETOPTS *go, struct cfg_s *cfg, int do_output_to_tmp, char *errbuf, CM_t *cm, seqs_to_aln_t *seqs_to_aln)
{
  int status;
  ESL_MSA *msa = NULL;
  int i, imax;
  float sc, struct_sc;

  /* contract check */
  if((do_output_to_tmp) && (esl_opt_GetBoolean(go, "--ileaved"))) ESL_FAIL(eslEINVAL, errbuf, "--ileaved enabled, but trying to output to temporary alignment file. This shouldn't happen.");
  if((! do_output_to_tmp) && (cfg->nali != 0)) ESL_FAIL(eslEINVAL, errbuf, "in output_result(), not outputting to tmp file but also not on first alignment. This shouldn't happen.");

  /* print per-CM info to insertfp and elfp, if nec */
  if((cfg->nali == 0) && (cfg->insertfp != NULL)) { fprintf(cfg->insertfp, "%s %d\n", cm->name, cm->clen); } 
  if((cfg->nali == 0) && (cfg->elfp != NULL))     { fprintf(cfg->elfp,     "%s %d\n", cm->name, cm->clen); } 

#ifdef HAVE_MPI
  /* if --mpi and ! -q, output the scores */
  if(esl_opt_GetBoolean(go, "--mpi") && (!esl_opt_GetBoolean(go, "-q"))) { 
    char *namedashes;
    int ni;
    int namewidth = 8; /* length of 'seq name' */
    /* determine the longest name in seqs_to_aln */
    for(ni = 0; ni < seqs_to_aln->nseq; ni++) namewidth = ESL_MAX(namewidth, strlen(seqs_to_aln->sq[ni]->name));
    ESL_ALLOC(namedashes, sizeof(char) * namewidth+1);
    namedashes[namewidth] = '\0';
    for(ni = 0; ni < namewidth; ni++) namedashes[ni] = '-';
    if(seqs_to_aln->struct_sc == NULL || (! NOT_IMPOSSIBLE(seqs_to_aln->struct_sc[0]))) { 
      if(cm->align_opts & CM_ALIGN_OPTACC) { 
	fprintf(stdout, "#\n");
	fprintf(stdout, "# %9s  %-*s  %5s  %8s  %8s\n", "seq idx",  namewidth, "seq name",   "len",  "bit sc",   "avg prob");
	fprintf(stdout, "# %9s  %-*s  %5s  %8s  %8s\n",  "---------", namewidth, namedashes, "-----", "--------", "--------");
      }
      else { 
	fprintf(stdout, "#\n");
	fprintf(stdout, "# %9s  %-*s  %5s  %8s\n",  "seq idx", namewidth, "seq name",   "len",  "bit sc");
	fprintf(stdout, "# %9s  %-*s  %5s  %8s\n",  "---------", namewidth, namedashes, "-----", "--------");
      }
    }
    else { /* we have struct scores */
      if(cm->align_opts & CM_ALIGN_OPTACC) { 
	fprintf(stdout, "#\n");
	fprintf(stdout, "# %9s  %-*s  %5s  %18s  %8s\n", "",         namewidth, "",                      "",       "    bit scores    ",   "");
	fprintf(stdout, "# %9s  %-*s  %5s  %18s  %8s\n", "",         namewidth, "",                      "",       "------------------",   "");
	fprintf(stdout, "# %9s  %-*s  %5s  %8s  %8s  %8s\n", "seq idx",  namewidth, "seq name",   "len", "total",    "struct",   "avg prob");
	fprintf(stdout, "# %9s  %-*s  %5s  %8s  %8s  %8s\n",  "---------", namewidth, namedashes, "-----", "--------", "--------", "--------");
      }
      else { 
	fprintf(stdout, "#\n");
	fprintf(stdout, "# %9s  %-*s  %5s  %18s\n", "", namewidth,       "",                  "",       "    bit scores    ");
	fprintf(stdout, "# %9s  %-*s  %5s  %18s\n", "", namewidth,       "",                  "",       "------------------");
	fprintf(stdout, "# %9s  %-*s  %5s  %8s  %8s\n",  "seq idx", namewidth,  "seq name",  "len",  "total",   "struct");
	    fprintf(stdout, "# %9s  %-*s  %5s  %8s  %8s\n",  "---------", namewidth, namedashes, "-----", "--------", "--------");
      }
    }
    free(namedashes);
    
    int imin;
    float null3_correction = 0.;
    for (i = 0; i < seqs_to_aln->nseq; i++) {
      if(!(esl_opt_GetBoolean(go, "--no-null3"))) { ScoreCorrectionNull3CompUnknown(cm->abc, cm->null, seqs_to_aln->sq[i]->dsq, 1, seqs_to_aln->sq[i]->n, &null3_correction); }
      if(! NOT_IMPOSSIBLE(seqs_to_aln->struct_sc[i])) { 
	if(cm->align_opts & CM_ALIGN_OPTACC) fprintf(stdout, "  %9d  %-*s  %5" PRId64 "  %8.2f  %8.3f\n", (cfg->nseq+i+1), namewidth, seqs_to_aln->sq[i]->name, seqs_to_aln->sq[i]->n, seqs_to_aln->sc[i] - null3_correction, seqs_to_aln->pp[i]);
	else                                 fprintf(stdout, "  %9d  %-*s  %5" PRId64 "  %8.2f\n",        (cfg->nseq+i+1), namewidth, seqs_to_aln->sq[i]->name, seqs_to_aln->sq[i]->n, seqs_to_aln->sc[i] - null3_correction);
      }
      else { /* we have struct scores */
	/* EPN, Thu Dec 17 14:08:20 2009 
	 * Following line would correct struct_sc for null3 penalty in an inexact way.
	 * I'm not sure if this is a good idea. For now its not done, but this line is 
	 * left here, commented out, for reference.
	 *
	 * adjust struct_sc for NULL3 correction, this is inexact 
	 * if(!(esl_opt_GetBoolean(go, "--no-null3"))) seqs_to_aln->struct_sc[i] -= ((float) ParsetreeCountMPEmissions(cm, seqs_to_aln->tr[i]) / (float) seqs_to_aln->sq[i]->n) * null3_correction; 
	 */ 
	if(cm->align_opts & CM_ALIGN_OPTACC) fprintf(stdout, "  %9d  %-*s  %5" PRId64 "  %8.2f  %8.2f  %8.3f\n", (cfg->nseq+i+1), namewidth, seqs_to_aln->sq[i]->name, seqs_to_aln->sq[i]->n, seqs_to_aln->sc[i] - null3_correction, seqs_to_aln->struct_sc[i], seqs_to_aln->pp[i]);
	else                                 fprintf(stdout, "  %9d  %-*s  %5" PRId64 "  %8.2f  %8.2f\n",        (cfg->nseq+i+1), namewidth, seqs_to_aln->sq[i]->name, seqs_to_aln->sq[i]->n, seqs_to_aln->sc[i] - null3_correction, seqs_to_aln->struct_sc[i]);
      }
    }
    fflush(stdout);
  }
#endif

  /* if nec, output the traces */
  if(cfg->tracefp != NULL) { 
    for (i = 0; i < seqs_to_aln->nseq; i++) { 
      fprintf(cfg->tracefp, "> %s\n", seqs_to_aln->sq[i]->name);
      if(esl_opt_GetBoolean(go,"--viterbi")) { 
	fprintf(cfg->tracefp, "  SCORE : %.2f bits\n", CP9TraceScore(cm->cp9, seqs_to_aln->sq[i]->dsq, seqs_to_aln->cp9_tr[i]));
	CP9PrintTrace(cfg->tracefp, seqs_to_aln->cp9_tr[i], cm->cp9, seqs_to_aln->sq[i]->dsq);
      }
      else { 
	if((status = ParsetreeScore(cm, NULL, errbuf, seqs_to_aln->tr[i], seqs_to_aln->sq[i]->dsq, FALSE, &sc, &struct_sc, NULL, NULL, NULL)) != eslOK) return status;
	fprintf(cfg->tracefp, "  %16s %.2f bits\n", "SCORE:", sc);
	fprintf(cfg->tracefp, "  %16s %.2f bits\n", "STRUCTURE SCORE:", struct_sc);
	ParsetreeDump(cfg->tracefp, seqs_to_aln->tr[i], cm, seqs_to_aln->sq[i]->dsq, NULL, NULL); /* NULLs are dmin, dmax */
      }
      fprintf(cfg->tracefp, "//\n");
    }
  }
  /* Optionally include a fixed alignment provided with --withali,
   * this has already been checked to see it matches the CM structure.
   * Note, we only do this if cfg->nali is 0, which means this is the 
   * first temporary alignment we're creating. (All temporary alignments
   * are merged at the end of cmalign's run). 
   */
  if((cfg->nali == 0) && (esl_opt_GetString(go, "--withali") != NULL))
    {
      /* grow the seqs_to_aln object */
      imax = seqs_to_aln->nseq;
      if((seqs_to_aln->nseq + cfg->withmsa->nseq) > seqs_to_aln->nalloc) 
	GrowSeqsToAln(seqs_to_aln, seqs_to_aln->nseq + cfg->withmsa->nseq - seqs_to_aln->nalloc, FALSE);
      if((status = include_withali(go, cfg, cm, &(seqs_to_aln->sq), &(seqs_to_aln->tr), &(seqs_to_aln->postcode), &(seqs_to_aln->nseq), errbuf)) != eslOK)
	ESL_FAIL(status, errbuf, "--withali alignment file %s doesn't have SS_cons annotation compatible with the CM\n", esl_opt_GetString(go, "--withali"));
    }
  
  /* create the msa, we free a parsetrees (or cp9 trace) as soon as we create each aligned sequence, to save memory */
  if(esl_opt_GetBoolean(go, "--viterbi"))
    {
      ESL_DASSERT1((seqs_to_aln->cp9_tr != NULL));
      if((status = CP9Traces2Alignment(cm, cfg->abc_out, seqs_to_aln->sq, NULL, seqs_to_aln->nseq, seqs_to_aln->cp9_tr, 
				       TRUE, esl_opt_GetBoolean(go, "--matchonly"), &msa)) != eslOK)
	goto ERROR;
    }
  else
    {
      assert(seqs_to_aln->tr != NULL);
      /*for(i = 0; i < seqs_to_aln->nseq; i++) { tr_Mb += SizeofParsetree(seqs_to_aln->tr[i]); }
	printf("All %d parsetrees are total size: %.6f Mb\n", seqs_to_aln->nseq, tr_Mb);*/
      if((status = Parsetrees2Alignment(cm, errbuf, cfg->abc_out, seqs_to_aln->sq, NULL, seqs_to_aln->tr, seqs_to_aln->postcode, 
					seqs_to_aln->nseq, cfg->insertfp, cfg->elfp, TRUE, esl_opt_GetBoolean(go, "--matchonly"), TRUE, &msa)) != eslOK)
	esl_fatal("Error creating the alignment: %s", errbuf);
      seqs_to_aln->tr = NULL; /* these were freed by Parsetrees2Alignment() */
    }
  
  if((! do_output_to_tmp) && (cfg->ofp == stdout) && (! esl_opt_GetBoolean(go, "-q"))) fprintf(cfg->ofp, "\n"); /* blank line to separate scores from aln if printing aln to stdout */
  
  /* if nec, replace msa->ss_cons with ss_cons from withmsa alignment */
  if((cfg->nali == 0) && esl_opt_GetBoolean(go, "--withpknots")) {
    if((status = add_withali_pknots(go, cfg, errbuf, cm, msa)) != eslOK) return status;
  }

  status = esl_msa_Write(do_output_to_tmp ? cfg->tmpfp : cfg->ofp, msa, (esl_opt_GetBoolean(go, "--ileaved") ? eslMSAFILE_STOCKHOLM : eslMSAFILE_PFAM));
  /* note that the contract asserted that if --ileaved then do_output_to_tmp must be FALSE */
  if      (status == eslEMEM) ESL_FAIL(status, errbuf, "Memory error when outputting alignment\n");
  else if (status != eslOK)   ESL_FAIL(status, errbuf, "Writing alignment file failed with error %d\n", status);

  /* if nec and (!do_output_to_tmp), then output the alignment to the regression file
   * (if do_output_to_tmp), we'll handle this in serial_master */
  if ((! do_output_to_tmp) && cfg->regressfp != NULL) {
    /* Must delete author info from msa, because it contains version
     * and won't diff clean in regression tests. */
    if(msa->au != NULL) free(msa->au); msa->au = NULL;
    status = esl_msa_Write(cfg->regressfp, msa, (esl_opt_GetBoolean(go, "--ileaved") ? eslMSAFILE_STOCKHOLM : eslMSAFILE_PFAM));
    if (status == eslEMEM)    ESL_FAIL(status, errbuf, "Memory error when outputting regression file\n");
    else if (status != eslOK) ESL_FAIL(status, errbuf, "Writing regression file failed with error %d\n", status);
  }
  
  if(msa != NULL) esl_msa_Destroy(msa);
  return eslOK;

 ERROR:
  if(msa != NULL) esl_msa_Destroy(msa);
  return status;
}

/* An alignment work unit consists a seqs_to_aln_t object which contains sequences to align, 
 * and space for their parsetrees, or CP9 traces, and possibly posterior code strings.
 * The job is to align the sequences and create parsetrees or cp9 traces and maybe posterior code strings.
 */
static int
process_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, 
		 seqs_to_aln_t *seqs_to_aln)
{
  int status;
  int be_quiet = esl_opt_GetBoolean(go, "-q");

#ifdef HAVE_MPI
  if(esl_opt_GetBoolean(go, "--mpi")) be_quiet = TRUE;
#endif

  if((status = DispatchAlignments(cm, errbuf, seqs_to_aln,
				  NULL, NULL, 0,  /* we're not aligning search hits */
				  esl_opt_GetInteger(go, "--banddump"),
				  esl_opt_GetInteger(go, "--dlev"), be_quiet, 
				  (! esl_opt_GetBoolean(go, "--no-null3")), cfg->r,
				  esl_opt_GetReal(go, "--mxsize"), stdout, cfg->scorefp, cfg->nseq+1,
				  esl_opt_GetInteger(go, "--7pad"), 
				  esl_opt_GetInteger(go, "--7len"), 
				  esl_opt_GetReal(go, "--7sc"), 
				  esl_opt_GetInteger(go, "--7end"), 
				  esl_opt_GetReal(go, "--7mprob"), 
				  esl_opt_GetReal(go, "--7mcprob"), 
				  esl_opt_GetReal(go, "--7iprob"), 
				  esl_opt_GetReal(go, "--7ilprob"))) != eslOK) goto ERROR;
  return eslOK;
  
 ERROR:
  ESL_DPRINTF1(("worker %d: has caught an error in process_search_workunit\n", cfg->my_rank));
  FreeCM(cm);
  return status;
}

/* initialize_cm()
 * Setup the CM based on the command-line options/defaults.
 * Configures the CM with a ConfigCM() call at end.
 */
static int
initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int status;
  int nstarts, nexits, nd;

  /* set up params/flags/options of the CM */
  cm->beta_qdb = esl_opt_GetReal(go, "--beta");
  cm->tau      = esl_opt_GetReal(go, "--tau");  /* this will be DEFAULT_TAU unless changed at command line */

  /* update cm->config_opts */
  if(esl_opt_GetBoolean(go, "-l"))
    {
      cm->config_opts |= CM_CONFIG_LOCAL;
      cm->config_opts |= CM_CONFIG_HMMLOCAL;
      cm->config_opts |= CM_CONFIG_HMMEL;
    }

  /* update cm->align_opts */
  /* optimal accuracy alignment is default */
  if(esl_opt_GetBoolean(go, "--optacc"))      cm->align_opts  |= CM_ALIGN_OPTACC;
  if(esl_opt_GetBoolean(go, "--sample"))      cm->align_opts  |= CM_ALIGN_SAMPLE;
  if(esl_opt_GetBoolean(go, "--hbanded"))     cm->align_opts  |= CM_ALIGN_HBANDED;
  if(esl_opt_GetBoolean(go, "--p7"))          cm->align_opts  |= CM_ALIGN_P7BANDED;
  if(esl_opt_GetBoolean(go, "--nonbanded"))   cm->align_opts  &= ~CM_ALIGN_HBANDED;
  if(esl_opt_GetBoolean(go, "--sub"))         cm->align_opts  |= CM_ALIGN_SUB;
  if(esl_opt_GetBoolean(go, "--viterbi"))     cm->align_opts  |= CM_ALIGN_HMMVITERBI;
  if(esl_opt_GetBoolean(go, "--small"))       cm->align_opts  |= CM_ALIGN_SMALL;
  if(esl_opt_GetBoolean(go, "--hsafe"))       cm->align_opts  |= CM_ALIGN_HMMSAFE;
  if(esl_opt_GetBoolean(go, "--checkpost"))   cm->align_opts  |= CM_ALIGN_CHECKINOUT;
  if(esl_opt_GetBoolean(go, "--checkfb"))     cm->align_opts  |= CM_ALIGN_CHECKFB;
  if(esl_opt_GetBoolean(go, "--sums"))        cm->align_opts  |= CM_ALIGN_SUMS;
  /* EPN, Mon Dec 14 18:22:18 2009 Deprecated --fins, it doesn't work with new memory efficiency 
   * if(esl_opt_GetBoolean(go, "--fins"))        cm->align_opts  |= CM_ALIGN_FLUSHINSERTS; */

  /* We can only do posteriors if we're not doing D&C (--small), --viterbi, 
   * nor --hsafe (which falls over to D&C if score is too low). Currently, these 
   * are all already enforced by ESL_GETOPTS, by having --small, --viterbi, 
   * and --hsafe all requiring --no-prob. This is a second line of defense in case 
   * the ESL_GETOPTS ever changes to remove one of those requirements. 
  */
  if((! esl_opt_GetBoolean(go, "--no-prob")) &&
     (! esl_opt_GetBoolean(go, "--small")) &&
     (! esl_opt_GetBoolean(go, "--viterbi")) &&
     (! esl_opt_GetBoolean(go, "--hsafe"))) { 
    cm->align_opts  |= CM_ALIGN_POST;
  }
  /* config QDB? */
  if(esl_opt_GetBoolean(go, "--qdb"))          
    { 
      cm->align_opts  |= CM_ALIGN_QDB;
      cm->config_opts |= CM_CONFIG_QDB;
    }

  /* BEGIN (POTENTIALLY) TEMPORARY BLOCK */
  /* set aggregate local begin/end probs, set with --pbegin, --pend, defaults are DEFAULT_PBEGIN, DEFAULT_PEND */
  cm->pbegin = esl_opt_GetReal(go, "--pbegin");
  cm->pend   = esl_opt_GetReal(go, "--pend");
  /* possibly overwrite local begin probs such that all begin points are equiprobable (--pebegin) */
  if(esl_opt_GetBoolean(go, "--pebegin")) {
    nstarts = 0;
    for (nd = 2; nd < cm->nodes; nd++) 
      if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd || cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BIF_nd) 
	nstarts++;
    /* printf("nstarts: %d\n", nstarts); */
    cm->pbegin = 1.- (1./(1+nstarts));
    /* printf("pbegin: %.5f\n", cm->pbegin); */
  }
  /* possibly overwrite cm->pend so that local end prob from all legal states is fixed,
   * this is strange in that cm->pend may be placed as a number greater than 1., this number
   * is then divided by nexits in ConfigLocalEnds() to get the prob for each v --> EL transition,
   * this is guaranteed by the way we calculate it to be < 1.,  it's the argument from --pfend */
  if( esl_opt_IsOn(go, "--pfend")) {
    nexits = 0;
    for (nd = 1; nd < cm->nodes; nd++) {
      if ((cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
	   cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BEGL_nd ||
	   cm->ndtype[nd] == BEGR_nd) && 
	  cm->ndtype[nd+1] != END_nd)
	nexits++;
    }
    cm->pend = nexits * esl_opt_GetReal(go, "--pfend");
  }
  /* END (POTENTIALLY) TEMPORARY BLOCK */


  /* finally, configure the CM for alignment based on cm->config_opts and cm->align_opts.
   * set local mode, make cp9 HMM, calculate QD bands etc. 
   */
  if((status = ConfigCM(cm, errbuf, FALSE, NULL, NULL)) != eslOK) return status; /* FALSE says do not calculate W unless nec b/c we're using QDBs */

  /* if(cfg->my_rank == 0) printf("CM %d: %s\n", (cfg->ncm), cm->name); 
   * debug_print_cm_params(stdout, cm);
   */

  /* if we're master and we're trying to include an alignment, make sure it is consistent with CM structure */
  if(cfg->my_rank == 0 && (esl_opt_IsOn(go, "--withali"))) { 
    if((status = check_withali(go, cfg, cm, &(cfg->withmsa), &(cfg->withali_mtr))) != eslOK)
      ESL_FAIL(status, errbuf, "--withali alignment file %s doesn't have a SS_cons compatible with the CM\n", esl_opt_GetString(go, "--withali"));
    }
  return eslOK;
}

/* Function: compare_cm_guide_trees()
 * EPN, Tue Mar  6 08:32:12 2007
 *
 * Purpose:  Given two CMs, cm1 and cm2, compare them, returning TRUE 
 *           iff they have the same guide tree (same node architecture).
 *
 * Args:     cm1          - covariance model number 1
 *           cm2          - covariance model number 2
 * 
 * Returns:  TRUE if CMs have same guide tree, FALSE otherwise
 */
static int compare_cm_guide_trees(CM_t *cm1, CM_t *cm2)
{
  int          nd; 
  if(cm1->nodes != cm2->nodes) return FALSE;
  for(nd = 0; nd < cm1->nodes; nd++)
    if(cm1->ndtype[nd] != cm2->ndtype[nd]) return FALSE;
  return TRUE;
}

/* Function: check_withali()
 * EPN, Tue Mar  6 06:25:02 2007
 *
 * Purpose:  Ensure that the alignment to include has a secondary
 *           structure that matches our CM. Pass the alignment back
 *           as *ret_msa.
 *
 * Returns:  <eslOK> on success.
 *           <eslEINCOMPAT> if alignment doesn't match the CM 
 */
static int check_withali(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, ESL_MSA **ret_msa, Parsetree_t **ret_mtr)
{
  int           status;
  ESL_MSA      *msa      = NULL; /* alignment we're including  */
  CM_t         *new_cm   = NULL; /* CM built from MSA, we check it has same guide tree as 'cm' */
  Parsetree_t  *mtr      = NULL; /* master structure tree from the alignment*/
  char          errbuf[cmERRBUFSIZE];

  /* cfg->withalifp is open */
  status = esl_msa_Read(cfg->withalifp, &msa);
  if      (status == eslEFORMAT) cm_Fail("--withali alignment file parse error:\n%s\n", cfg->withalifp->errbuf);
  else if (status == eslEINVAL)  cm_Fail("--withali alignment file parse error:\n%s\n", cfg->withalifp->errbuf);
  else if (status == eslEOF)     cm_Fail("--withali alignment file %s empty?\n",        cfg->withalifp->fname);
  else if (status != eslOK)      cm_Fail("--withali alignment file read failed with error code %d\n", status);

  /* Some input data cleaning. */
  if (esl_opt_GetBoolean(go, "--rf") && msa->rf == NULL) 
    ESL_FAIL(eslFAIL, errbuf, "--rf invoked but --withali alignment has no reference coord annotation.\n");
  if (msa->ss_cons == NULL) 
    ESL_FAIL(eslFAIL, errbuf, "--withali alignment did not contain consensus structure annotation.\n");
  if (esl_opt_GetBoolean(go, "--withpknots")) /* copy the original secondary structure */
    esl_strdup(msa->ss_cons, -1, &(cfg->withss_cons));
 if (! clean_cs(msa->ss_cons, msa->alen, TRUE))
    ESL_FAIL(eslFAIL, errbuf, "Failed to parse consensus structure annotation for --withali alignment\n");

  /* Build a CM from a master guide tree built from the msa, 
   * then check to make sure this CM has same emit map as the CM
   * we've had passed in. This is fragile and hopefully temporary. 
   * Another solution would be to use a checksum, but CM files don't 
   * have checksums yet.
   */
  if((status = HandModelmaker(msa, errbuf, esl_opt_GetBoolean(go, "--rf"), esl_opt_GetReal(go, "--gapthresh"), &new_cm, &mtr)) != eslOK) return status;
  if(!(compare_cm_guide_trees(cm, new_cm))) {
    /* no need to try rebalancing, that doesn't change the guidetree (seriously this is from cm.c::CMRebalance():
     * for (x = 0; x < new->nodes; x++) new->ndtype[x]  = cm->ndtype[x];
     */
    status = eslEINCOMPAT;
    goto ERROR;
  }

  /* if we get here, the CM guide trees match */
  if(new_cm   != NULL) FreeCM(new_cm);
  *ret_mtr = mtr;
  *ret_msa = msa;
  return eslOK;

 ERROR:
  if(msa != NULL)      esl_msa_Destroy(msa);
  if(new_cm   != NULL) FreeCM(new_cm);
  if(mtr != NULL)      FreeParsetree(mtr);
  return eslEINCOMPAT;
}

/* Function: include_withali()
 * EPN, Tue Mar  6 06:25:02 2007
 *
 * Purpose:  Determine the implicit parses of sequences in an
 *           MSA to a CM and append them to passed in data structures.
 *           We've already checked to make sure the MSA's consensus
 *           structure matches the CM.
 *
 * Args:     go           - command line options
 *           cfg          - cmalign configuration, includes msa to add
 *           cm           - CM we're aligning to 
 *           ret_sq       - pre-existing sequences, to append msa seqs to
 *           ret_tr       - pre-existing parsetrees, to append to
 *           ret_postcode - pre-existing posterior codes, to append to
 *           ret_nseq     - number of exisiting seqs, updated at end of function
 *           errbuf       - easel error message
 * 
 * Returns:  eslOK on success, eslEMEM on memory error
 *           Also new, realloc'ed arrays for sq, tr in ret_seq, ret_tr; 
 *           *ret_nseq is increased by number of seqs in cfg->withmsa.
*/
static int include_withali(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, ESL_SQ ***ret_sq, Parsetree_t ***ret_tr, char ***ret_postcode, int *ret_nseq, char *errbuf)
{
  int           status;
  void         *tmp;      /* for ESL_RALLOC() */
  int           i;	  /* counter over aseqs       */
  int           ip;	  /* offset counter over aseqs */
  char        **uaseq;    /* unaligned seqs, dealigned from the MSA */
  char         **aseq;    /*   aligned text seqs */
  int           apos;     /*   aligned position index */
  int           uapos;    /* unaligned position index */
  int           x;        /* counter of parsetree nodes */
  int         **map;      /* [0..msa->nseq-1][0..msa->alen] map from aligned
			   * positions to unaligned (non-gap) positions */
  int           do_post;  /* TRUE if we need to worry about post codes */
  /* for swapping pts at end of func so seqs from withali appear at top of aln */
  Parsetree_t **tmp_tr = NULL;  
  ESL_SQ      **tmp_sq = NULL;
  char        **tmp_postcode = NULL;

  /* Contract check */
  if(cfg->withmsa == NULL)                     ESL_XFAIL(eslEINVAL, errbuf, "Error in include_withali() withmsa is NULL.\n");
  if(! (cfg->withmsa->flags & eslMSA_DIGITAL)) ESL_XFAIL(eslEINVAL, errbuf, "Error in include_withali() withmsa is not digitized.\n");
  do_post = (*ret_postcode != NULL) ? TRUE : FALSE;

  /* For each seq in the MSA, map the aligned sequences coords to 
   * the unaligned coords, we stay in digitized seq coords (1..alen),
   * we need this for converting parsetrees from Transmogrify (which
   * have emitl and emitr in aligned coords) to unaligned coords, so 
   * we can call Parsetrees2Alignment() with them. */
  ESL_ALLOC(map,   sizeof(int *)  * cfg->withmsa->nseq);
  ESL_ALLOC(uaseq, sizeof(char *) * cfg->withmsa->nseq);
  ESL_ALLOC(aseq,  sizeof(char *) * cfg->withmsa->nseq);
  for (i = 0; i < cfg->withmsa->nseq; i++)
    {
      ESL_ALLOC(map[i],   sizeof(int)  * (cfg->withmsa->alen+1));
      ESL_ALLOC(aseq[i],  sizeof(char) * (cfg->withmsa->alen+1));
      map[i][0] = -1; /* invalid */
      uapos = 1;
      for(apos = 1; apos <= cfg->withmsa->alen; apos++)
	{
	  if (!esl_abc_XIsGap(cfg->withmsa->abc, cfg->withmsa->ax[i][apos]))
	    map[i][apos] = uapos++;
	  else
	    map[i][apos] = -1;
	}
      /* we need digitized AND text seqs for Transmogrify */
      esl_abc_Textize(cfg->withmsa->abc, cfg->withmsa->ax[i], cfg->withmsa->alen, aseq[i]);
      esl_strdup(aseq[i], -1, &(uaseq[i]));
      esl_strdealign(uaseq[i], uaseq[i], "-_.", NULL);
    }
  ESL_RALLOC((*ret_tr),  tmp, (sizeof(Parsetree_t *)  * (*ret_nseq + cfg->withmsa->nseq)));
  ESL_RALLOC((*ret_sq),  tmp, (sizeof(ESL_SQ *)       * (*ret_nseq + cfg->withmsa->nseq)));
  /* if do_post, check to see if withmsa has posterior annotation, if so, store it */
  if(do_post) { 
    ESL_RALLOC((*ret_postcode),  tmp, (sizeof(char *) * (*ret_nseq + cfg->withmsa->nseq)));
    for (i = *ret_nseq; i < (*ret_nseq + cfg->withmsa->nseq); i++) { 
      (*ret_postcode)[i] = NULL;
    }
    if(cfg->withmsa->pp != NULL) { 
      ip = i - *ret_nseq;
      if(cfg->withmsa->pp[ip] != NULL) { 
	esl_strdup(cfg->withmsa->pp[ip], -1, &((*ret_postcode)[i]));
	esl_strdealign((*ret_postcode)[i], aseq[ip], "-_.~", NULL);
      }
    }
  }

  /* Transmogrify each aligned seq to get a parsetree */
  for (i = *ret_nseq; i < (*ret_nseq + cfg->withmsa->nseq); i++)
    {
      ip = i - *ret_nseq;
      (*ret_tr)[i] = Transmogrify(cm, cfg->withali_mtr, cfg->withmsa->ax[ip], aseq[ip], cfg->withmsa->alen);
      /* ret_tr[i] is in alignment coords, convert it to unaligned coords, */
      for(x = 0; x < (*ret_tr)[i]->n; x++)
	{
	  if((*ret_tr)[i]->emitl[x] != -1)
	    (*ret_tr)[i]->emitl[x] = map[ip][(*ret_tr)[i]->emitl[x]];
	  if((*ret_tr)[i]->emitr[x] != -1)
	    (*ret_tr)[i]->emitr[x] = map[ip][(*ret_tr)[i]->emitr[x]];
	}
      (*ret_sq)[i]      = esl_sq_CreateFrom(cfg->withmsa->sqname[ip], uaseq[ip], NULL, NULL, NULL);
      esl_sq_Digitize(cm->abc, (*ret_sq)[i]);
    }

  /* Swap some pointers so the included alignment appears at the top of the output 
   * alignment instead of the bottom. */
  ESL_ALLOC(tmp_tr, sizeof(Parsetree_t *) * (*ret_nseq + cfg->withmsa->nseq));
  ESL_ALLOC(tmp_sq, sizeof(ESL_SQ *)      * (*ret_nseq + cfg->withmsa->nseq));
  if(do_post) { 
    ESL_ALLOC(tmp_postcode, sizeof(char *) * (*ret_nseq + cfg->withmsa->nseq));
  }
  for(i = 0; i < (*ret_nseq + cfg->withmsa->nseq); i++)
    {
      tmp_tr[i] = (*ret_tr)[i];
      tmp_sq[i] = (*ret_sq)[i];
      if(do_post) { 
	tmp_postcode[i] = (*ret_postcode)[i];
      }
    }
  for(i = 0; i < *ret_nseq; i++)
    {
      ip = i + cfg->withmsa->nseq;
      (*ret_tr)[ip] = tmp_tr[i];
      (*ret_sq)[ip] = tmp_sq[i];
      if(do_post) { 
	(*ret_postcode)[ip] = tmp_postcode[i];
      }
    }
  for(i = *ret_nseq; i < (*ret_nseq + cfg->withmsa->nseq); i++)
    {
      ip = i - *ret_nseq;
      (*ret_tr)[ip] = tmp_tr[i];
      (*ret_sq)[ip] = tmp_sq[i];
      if(do_post) { 
	(*ret_postcode)[ip] = tmp_postcode[i];
      }
    }
  /* update *ret_nseq */
  *ret_nseq    += cfg->withmsa->nseq;

  /* Clean up and exit. */
  esl_Free2D((void **) map,   cfg->withmsa->nseq);
  esl_Free2D((void **) uaseq, cfg->withmsa->nseq);
  esl_Free2D((void **) aseq,  cfg->withmsa->nseq);
  free(tmp_tr);
  free(tmp_sq);
  if(do_post) { 
    free(tmp_postcode);
  }
  return eslOK;

 ERROR:
  esl_Free2D((void **) map,   cfg->withmsa->nseq);
  esl_Free2D((void **) uaseq, cfg->withmsa->nseq);
  if(tmp_sq != NULL) free(tmp_sq);
  if(tmp_tr != NULL) free(tmp_tr);
  if(tmp_tr != NULL) free(tmp_tr);
  if(tmp_postcode != NULL) free(tmp_postcode);

  return status;
}

/* Function: add_withali_pknots()
 * EPN, Wed Oct 17 18:24:33 2007
 *
 * Purpose:  Determine the pseudoknots that were in consensus columns of
 *           the --withali alignment and add them to newmsa->ss_cons.
 *
 * Args:     go           - command line options
 *           cfg          - cmalign configuration, includes msa to add
 *           errbuf       - easel error message
 *           cm           - the CM, only used to get cm->clen 
 *           newmsa       - MSA from Parsetrees2Alignment(), we want to add to it's ss_cons
 * 
 * Returns:  eslOK on success, eslEMEM on memory error
 *           eslEINVAL on unpredicted situation
*/
static int add_withali_pknots(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, ESL_MSA *newmsa)
{
  int           status;
  int           apos;     /*   aligned position index */
  int           cpos;     /* consensus column index */
  int           ngaps;
  float         gapthresh;
  int           withmsa_clen = 0;
  int           idx;     /* sequence index */
  int           i, j, i_cpos, j_cpos; /* residue position indices */
  /* Contract check */
  if(cfg->withmsa == NULL) cm_Fail("ERROR in add_withali_pknots() cfg->withmsa is NULL.\n");
  if(cfg->withss_cons  == NULL) cm_Fail("ERROR in add_withali_pknots() cfg->withss_cons is NULL.\n");
  if(! (cfg->withmsa->flags & eslMSA_DIGITAL)) cm_Fail("ERROR in add_withali_pknots() cfg->withmsa is not digitized.\n");

  /* 10 easy, convoluted steps. One reason for so many steps is 
   * we can't build ss_cons strings from pseudoknotted ct arrays,
   * so we have to work around that. 
   */

  /* 1. determine consensus columns of withmsa.  
   * If we've gotten this far, there should be same number as
   * cm->clen. code block stolen from cm_modelmaker.c matassign is
   * 1..alen. Values are 1 if a match (consensus) column, 0 if insert
   * column.
   */
  gapthresh = esl_opt_GetReal(go, "--gapthresh");
  int *matassign;
  ESL_ALLOC(matassign, sizeof(int) * (cfg->withmsa->alen+1));
  /* Watch for off-by-one. rf is [0..alen-1]; matassign is [1..alen] */
  if (esl_opt_GetBoolean(go, "--rf")) {
    for (apos = 1; apos <= cfg->withmsa->alen; apos++) 
      matassign[apos] = (esl_abc_CIsGap(cfg->withmsa->abc, cfg->withmsa->rf[apos-1])? FALSE : TRUE);
  }
  else { /* --rf not enabled */
    for (apos = 1; apos <= cfg->withmsa->alen; apos++) {
      for (ngaps = 0, idx = 0; idx < cfg->withmsa->nseq; idx++)
	if (esl_abc_XIsGap(cfg->withmsa->abc, cfg->withmsa->ax[idx][apos])) ngaps++;
      matassign[apos] = ((double) ngaps / (double) cfg->withmsa->nseq > gapthresh) ? 0 : 1;
    }
  }
  for (apos = 1; apos <= cfg->withmsa->alen; apos++) withmsa_clen += matassign[apos];
  if(withmsa_clen != cm->clen) ESL_FAIL(eslFAIL, errbuf, "withmsa consensus length != cm consensus length. A previous check for this passed, this is a coding error.");

  /* 2. get ct array for consensus structure of withmsa BEFORE we stripped away it's pknots,
   * this was saved in check_withali().
   */
  int *ct;
  ESL_ALLOC(ct, (cfg->withmsa->alen+1) * sizeof(int));
  if (esl_wuss2ct(cfg->withss_cons, cfg->withmsa->alen, ct) != eslOK)  
    ESL_FAIL(eslFAIL, errbuf, "withmsa original ss_cons inconsistent. A previous check for this passed, this is a coding error.");

  /* 3. also get a ct with no pknots, we'll need this to figure out where the pknots go */
  int *ct_noknots;
  ESL_ALLOC(ct_noknots, (cfg->withmsa->alen+1) * sizeof(int));
  if (esl_wuss2ct(cfg->withmsa->ss_cons, cfg->withmsa->alen, ct_noknots) != eslOK)  
    ESL_FAIL(eslFAIL, errbuf, "withmsa original ss_cons inconsistent. A previous check for this passed, this is a coding error.");

  /* 4. remove any base pairs (i,j) from ct and ct_noknots for which i or j are non-consensus columns */
  for (apos = 1; apos <= cfg->withmsa->alen; apos++) {
    if(! matassign[apos]) { /* apos is not a consensus column */
      if(ct[apos] != 0) ct[ct[apos]] = 0;
      ct[apos] = 0;
      if(ct_noknots[apos] != 0) ct_noknots[ct_noknots[apos]] = 0;
      ct_noknots[apos] = 0;
    }
  }

  /* 5. get a map from alignment coords to consensus coords */
  int *a2c_map;
  ESL_ALLOC(a2c_map, sizeof(int) * (cfg->withmsa->alen + 1));
  a2c_map[0] = -1;
  cpos = 1;
  for(apos = 1; apos <= cfg->withmsa->alen; apos++) {
    if(matassign[apos]) a2c_map[apos] = cpos++;
    else a2c_map[apos] = 0;
  }

  /* 6. use a2c_map to create c_ct_noknots, which is just ct_noknots
   * changed from aligned coordinates to consensus column
   * coordinates. remember no non-consensus column cpos should have
   * ct_noknots[cpos] != 0, because we stripped those bps.
   */
  int *c_ct_noknots;
  ESL_ALLOC(c_ct_noknots, sizeof(int) * (cm->clen+1));
  esl_vec_ISet(c_ct_noknots, (cm->clen+1), 0);
  for(apos = 1; apos <= cfg->withmsa->alen; apos++) {
    i = apos; j = ct_noknots[i];
    if(j != 0 && i < j) { /* if i > j, we've already covered it */
      if(a2c_map[i] == 0) ESL_FAIL(eslFAIL, errbuf, "withmsa apos: %d has structure, but is not consensus, this should never happen.", i);
      if(a2c_map[j] == 0) ESL_FAIL(eslFAIL, errbuf, "withmsa apos: %d has structure, but is not consensus, this should never happen.", j);
      c_ct_noknots[a2c_map[i]] = a2c_map[j];
      c_ct_noknots[a2c_map[j]] = a2c_map[i];
    }      
  }
  for(cpos = 1; cpos <= cm->clen; cpos++)
    printf("ct[%d]: %d\n", cpos, c_ct_noknots[cpos]);
  
  /* 7. build new consensus structure, with no knots and only consensus columns */
  char *c_sscons;
  ESL_ALLOC(c_sscons, sizeof(char) * (cm->clen + 1));
  if((status = esl_ct2wuss(c_ct_noknots, cm->clen, c_sscons)) != eslOK) cm_Fail("ct2wuss failed with (supposedly) no knots");
  
  /* 8. add back in consensus knots, using a2c_map */
  for(apos = 1; apos <= cfg->withmsa->alen; apos++) {
    if(matassign[apos]) {
      i = apos; j = ct[i];
      if(ct[i] != 0 && i < j) { /* if i > j, we've already updated it */
	i_cpos = a2c_map[i];
	j_cpos = a2c_map[j];
	/* printf("\t\tct bp i: %3d i_cpos: %3d j_cpos: %d j_cpos: %3d\n", i, i_cpos, j, j_cpos);
	   printf("\t\tc_sscons[i_cpos-1]: %c\n", c_sscons[(i_cpos-1)]);
	   printf("\t\tc_sscons[j_cpos-1]: %c\n", c_sscons[(j_cpos-1)]);
	   printf("\t\tcfg->withss_cons[(i-1)]: %c\n", cfg->withss_cons[(i-1)]);
	   printf("\t\tcfg->withss_cons[(j-1)]: %c\n", cfg->withss_cons[(j-1)]);
	*/
	c_sscons[(i_cpos-1)] = cfg->withss_cons[(i-1)]; /* add pknot annotation for left bp */
	c_sscons[(j_cpos-1)] = cfg->withss_cons[(j-1)]; /* add pknot annotation for right bp */
      }
    }
  }

  /* 9. c_sscons is the pknotted structure, but limited to the consensus columns,
   * final step is to overwrite newmsa->ss_cons characters for consensus columns 
   * only, by simply replacing them with characters from c_sscons.
   * we need a new map from consensus columns to newmsa align coords first,
   * we use the fact that newmsa->rf columns that are non-gapped are consensus
   * columns.
   */
  int *new_c2a_map;
  ESL_ALLOC(new_c2a_map, sizeof(int) * (cm->clen + 1));
  new_c2a_map[0] = -1;
  cpos = 0;
  for(apos = 1; apos <= newmsa->alen; apos++) 
    if(! esl_abc_CIsGap(newmsa->abc, newmsa->rf[(apos-1)])) new_c2a_map[++cpos] = apos;

  /* 10. overwrite newmsa->ss_cons */
  for(cpos = 1; cpos <= cm->clen; cpos++) 
    newmsa->ss_cons[(new_c2a_map[cpos]-1)] = c_sscons[(cpos-1)];
  
  /* free memory and return */
  free(new_c2a_map);
  free(c_sscons);
  free(a2c_map);
  free(ct_noknots);
  free(c_ct_noknots);
  free(ct);
  free(matassign);
  return eslOK;

 ERROR:
  return status;
}

/* Function: print_run_info
 * Date:     EPN, Thu Feb 28 14:26:42 2008
 *
 * Purpose:  Print information on this run of cmalign.
 *           Command used to run it, and execution date.
 *
 * Returns:  eslOK on success
 */
static int
print_run_info(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf)
{
  int status;
  char *command;
  char *date;

  if(esl_opt_GetBoolean(go, "-q")) return eslOK;

  if((status = get_command(go, errbuf, &command)) != eslOK) return status;
  if((status = GetDate    (errbuf, &date))        != eslOK) return status;

  fprintf(stdout, "%-10s %s\n",  "# command:", command);
  fprintf(stdout, "%-10s %s\n",  "# date:",    date);
  if(cfg->nproc > 1) fprintf(stdout, "# %-8s %d\n", "nproc:", cfg->nproc);
  if(esl_opt_GetBoolean(go, "--sample")) fprintf(stdout, "%-10s %" PRIu32 "\n", "# seed:", esl_randomness_GetSeed(cfg->r));
  fprintf(stdout, "#\n");

  if(cfg->scorefp != NULL) { 
    fprintf(cfg->scorefp, "%-10s %s\n",  "# command:", command);
    fprintf(cfg->scorefp, "%-10s %s\n",  "# date:",    date);
    if(cfg->nproc > 1) fprintf(cfg->scorefp, "# %-8s %d\n", "nproc:", cfg->nproc);
    if(esl_opt_GetBoolean(go, "--sample")) fprintf(cfg->scorefp, "%-10s %" PRIu32 "\n", "# seed:", esl_randomness_GetSeed(cfg->r));
    fprintf(cfg->scorefp, "#\n");
  }

  free(command);
  free(date);
  return eslOK;
}

/* Function: print_cm_info
 * Date:     EPN, Thu Feb 28 14:44:17 2008
 *
 * Purpose:  Print per-CM info to output file (stdout unless -o). 
 */
static void
print_cm_info(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm)
{

  if(esl_opt_GetBoolean(go, "-q")) return;

  int do_hbanded = (cm->align_opts & CM_ALIGN_HBANDED) ? TRUE : FALSE;
  int do_qdb     = (cm->align_opts & CM_ALIGN_QDB)     ? TRUE : FALSE;
  if(cm->align_opts & CM_ALIGN_HMMVITERBI) { do_hbanded = do_qdb = FALSE; }

  fprintf(stdout, "# %-25s  %9s  %6s  %3s  %5s  %6s\n", "cm name",                   "algorithm", "config", "sub", "bands", (do_hbanded) ? "tau" : ((do_qdb) ? "beta" : "")); 
  fprintf(stdout, "# %-25s  %9s  %6s  %3s  %5s  %6s\n", "-------------------------", "---------", "------", "---", "-----", (do_hbanded || do_qdb) ? "------" : ""); 
  fprintf(stdout, "# %-25.25s  %9s  %6s  %3s", 
	  cm->name,
	  ((esl_opt_GetBoolean(go, "--cyk")) ? "cyk" : ((esl_opt_GetBoolean(go, "--sample")) ? "sample" : (esl_opt_GetBoolean(go, "--viterbi") ? "hmm vit" : "opt acc"))), 
	  (esl_opt_GetBoolean(go, "-l")) ? "local" : "global",
	  (esl_opt_GetBoolean(go, "--sub")) ? "yes" : "no");
  /* bands and beta/tau */
  if     (do_hbanded)    fprintf(stdout, "  %5s  %6.0e\n", "hmm", cm->tau);
  else if(do_qdb)        fprintf(stdout, "  %5s  %6.0e\n", "qdb", cm->beta_qdb);
  else                   fprintf(stdout, "  %5s  %6s\n", "none", "");

  if(cfg->scorefp != NULL) { 
    fprintf(cfg->scorefp, "# %-25s  %9s  %6s  %3s  %5s  %6s\n", "cm name",                   "algorithm", "config", "sub", "bands", (do_hbanded) ? "tau" : ((do_qdb) ? "beta" : "")); 
    fprintf(cfg->scorefp, "# %-25s  %9s  %6s  %3s  %5s  %6s\n", "-------------------------", "---------", "------", "---", "-----", (do_hbanded || do_qdb) ? "------" : ""); 
    fprintf(cfg->scorefp, "# %-25.25s  %9s  %6s  %3s", 
	    cm->name,
	    ((esl_opt_GetBoolean(go, "--cyk")) ? "cyk" : ((esl_opt_GetBoolean(go, "--sample")) ? "sample" : (esl_opt_GetBoolean(go, "--viterbi") ? "hmm vit" : "opt acc"))), 
	    (esl_opt_GetBoolean(go, "-l")) ? "local" : "global",
	    (esl_opt_GetBoolean(go, "--sub")) ? "yes" : "no");
    /* bands and beta/tau */
    if     (do_hbanded)    fprintf(cfg->scorefp, "  %5s  %6.0e\n", "hmm", cm->tau);
    else if(do_qdb)        fprintf(cfg->scorefp, "  %5s  %6.0e\n", "qdb", cm->beta_qdb);
    else                   fprintf(cfg->scorefp, "  %5s  %6s\n", "none", "");
  }
  return;
}

/* Function: get_command
 * Date:     EPN, Fri Jan 25 13:56:10 2008
 *
 * Purpose:  Return the command used to call cmscore
 *           in <ret_command>.
 *
 * Returns:  eslOK on success; eslEMEM on allocation failure.
 */
int 
get_command(const ESL_GETOPTS *go, char *errbuf, char **ret_command)
{
  int status;
  int i;
  char *command = NULL;

  for (i = 0; i < go->argc; i++) { /* copy all command line options and args */
    if((status = esl_strcat(&(command),  -1, go->argv[i], -1)) != eslOK) goto ERROR;
    if(i < (go->argc-1)) if((status = esl_strcat(&(command), -1, " ", 1)) != eslOK) goto ERROR;
  }
  *ret_command = command;

  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "get_command(): memory allocation error.");
  return status;
}

/* Function: print_info_file_header
 * Date:     EPN, Fri Dec  4 08:15:31 2009
 *
 * Purpose:  Print the header section of an insert or EL insert
 *           (--ifile, --elfile) information file.
 *
 * Returns:  void
 */
void
print_info_file_header(FILE *fp, char *firstline, char *elstring)
{
  fprintf(fp, "# %s\n", firstline);
  fprintf(fp, "# This file includes 2+<nseq> non-'#' pre-fixed lines per model used for alignment,\n");
  fprintf(fp, "# where <nseq> is the number of sequences in the target file.\n");
  fprintf(fp, "# The first non-'#' prefixed line per model includes 2 tokens, separated by a single space (' '):\n");
  fprintf(fp, "# The first token is the model name and the second is the consensus length of the model (<clen>).\n");
  fprintf(fp, "# The following <nseq> lines include (4+3*<n>) whitespace delimited tokens per line.\n");
  fprintf(fp, "# The format for these <nseq> lines is:\n");
  fprintf(fp, "#   <seqname> <seqlen> <spos> <epos> <c_1> <u_1> <i_1> <c_2> <u_2> <i_2> .... <c_x> <u_x> <i_x> .... <c_n> <u_n> <i_n>\n");
  fprintf(fp, "#   indicating <seqname> has >= 1 %sinserted residues after <n> different consensus positions,\n", elstring);
  fprintf(fp, "#   <seqname> is the name of the sequence\n");
  fprintf(fp, "#   <seqlen>  is the unaligned length of the sequence\n");
  fprintf(fp, "#   <spos>    is the first (5'-most) consensus position filled by a nongap for this sequence (-1 if 0 nongap consensus posns)\n");
  fprintf(fp, "#   <epos>    is the final (3'-most) consensus position filled by a nongap for this sequence (-1 if 0 nongap consensus posns)\n");
  fprintf(fp, "#   <seqlen>  is the unaligned length of the sequence\n");
  fprintf(fp, "#   <c_x> is a consensus position (between 0 and <clen>; if 0: inserts before 1st consensus posn)\n");
  fprintf(fp, "#   <u_x> is the *unaligned* position (b/t 1 and <seqlen>) in <seqname> of the first %sinserted residue after <c_x>.\n", elstring);
  fprintf(fp, "#   <i_x> is the number of %sinserted residues after position <c_x> for <seqname>.\n", elstring);
  fprintf(fp, "# Lines for sequences with 0 %sinserted residues will include only <seqname> <seqlen> <spos> <epos>.\n", elstring);
  fprintf(fp, "# The final non-'#' prefixed line per model includes only '//', indicating the end of info for a model.\n");
  fprintf(fp, "#\n");

  return;
}


/* Function: create_and_output_final_msa
 * Incept:   EPN, Mon Dec 14 05:35:51 2009
 *
 * Purpose:  Read the >=1 MSAs that were written to a temporary file,
 *           merge them and output the merged MSA to a file without
 *           storing any of the full MSAs (incl. the final one) in
 *           memory.  To accomplish this a first pass of reading is
 *           done to determine how many gap columns must be added to
 *           each MSA to create the merged MSA during which only non
 *           per-sequence information is stored. After this pass, with
 *           the size of the merged alignment known, a second pass
 *           occurs during which only GS annotation is regurgitated
 *           (if any exists in at least 1 aln). Then a final pass
 *           occurs during which all other per-sequence data (PPs,
 *           aligned seqs) are regurgitated, taking care to add gap
 *           columns as necessary to make each input alignment the
 *           correct width of the merged alignment.
 *
 * Args:     go  - options
 *           cfg - cmalign config
 *           ofp - file to write to, this is usually cfg->ofp, but
 *                 can be cfg->regressfp (if --regress)
 *           errbuf - for error messages
 *           cm  - CM used for alignment, useful for cm->clen
 *           tmpfile - name of temporary file with alignments to merge
 * 
 * Returns:   <eslOK> on success. 
 *            Returns <eslEOF> if there are no more alignments in <afp>.
 *            <eslEFORMAT> if parse fails because of a file format problem,
 *            in which case afp->errbuf is set to contain a formatted message 
 *            that indicates the cause of the problem. <eslEMEM> on allocation
 *            error.
 *
 * Xref:      /groups/eddy/home/nawrockie/notebook/9_1211_inf_cmalign_memeff/
 */
int 
create_and_output_final_msa(const ESL_GETOPTS *go, const struct cfg_s *cfg, FILE *ofp, char *errbuf, CM_t *cm, char *tmpfile) 
{
  int           status;
  int           ai;                            /* counters over alignments */
  int           nseq_tot;                      /* number of sequences in all alignments */
  int           nseq_cur;                      /* number of sequences in current alignment */
  int64_t       alen_cur;                      /* length of current alignment */
  int64_t      *alenA = NULL;                  /* [0..nali_tot-1] alignment length of input msas (even after 
						* potentially removingeinserts (--rfonly)) */
  ESL_MSA     **msaA = NULL;                   /* [0..nali_tot-1] all msas read from all files */
  int          *maxins = NULL;                 /* [0..cpos..cm->clen+1] max number of inserts 
						* before each consensus position in all alignments */
  int          *maxel = NULL;                  /* [0..cpos..cm->clen+1] max number of EL inserts ('~' missing data symbols) 
						* before each consensus position in all alignments */
  int           cur_clen;                      /* consensus length (non-gap #=GC RF length) of current alignment */
  int           apos;                          /* alignment position */
  ESL_MSA      *fmsa = NULL;                   /* the merged alignment created by merging all alignments in msaA */
  int           alen_fmsa;                     /* number of columns in merged MSA */
  int          *ngap_insA = NULL;               /* [0..alen] number of insert gap columns to add after each alignment column when merging */
  int          *ngap_elA = NULL;                /* [0..alen] number of missing data ('~') gap columns to add after each alignment column when merging */
  int          *ngap_eitherA = NULL;            /* [0..apos..alen] = ngap_insA[apos] + ngap_elA[apos] */
  char         *rf2print = NULL;                /* #=GC RF annotation for final alignment */
  char         *ss_cons2print = NULL;           /* #=GC SS_cons annotation for final alignment */

  /* variables only used in small mode (--savemem) */
  int           ngs_cur;                       /* number of GS lines in current alignment (only used if do_small) */
  int           gs_exists = FALSE;             /* set to TRUE if do_small and any input aln has >= 1 GS line */
  int           maxname, maxgf, maxgc, maxgr;  /* max length of seqname, GF tag, GC tag, GR tag in all input alignments */
  int           maxname_cur, maxgf_cur, maxgc_cur, maxgr_cur; /* max length of seqname, GF tag, GC tag, GR tag in current input alignment */
  int           margin = 0;                    /* total margin length for output msa */
  int           regurg_header = FALSE;         /* set to TRUE if we're printing out header */
  int           regurg_gf     = FALSE;         /* set to TRUE if we're printing out GF (we won't if ofp == cfg->regressfp (i.e. we're printing to regress file */
  ESL_MSAFILE  *afp;

  /* Allocate and initialize */
  ESL_ALLOC(msaA,   sizeof(ESL_MSA *) * cfg->nali);
  ESL_ALLOC(alenA,  sizeof(int64_t) * cfg->nali);

  /****************************************************************************
   * Read alignments one at a time, storing all non-sequence info, separately *
   ****************************************************************************/
  if((status = esl_msafile_Open(tmpfile, eslMSAFILE_PFAM, NULL, &afp)) != eslOK) cm_Fail("unable to open temp file %s for reading", tmpfile);

  ai = 0;
  nseq_tot = 0;
  maxname = maxgf = maxgc = maxgr = 0;

  /* allocate maxins */
  ESL_ALLOC(maxins, sizeof(int) * (cm->clen+1)); 
  esl_vec_ISet(maxins, (cm->clen+1), 0);
  /* allocate maxel */
  ESL_ALLOC(maxel, sizeof(int) * (cm->clen+1)); 
  esl_vec_ISet(maxel, (cm->clen+1), 0); /* these will all stay 0 unless we see '~' in the alignments */

  /* read all alignments, there should be cfg->nali of them */
  for(ai = 0; ai < cfg->nali; ai++) { 
    status = esl_msa_ReadInfoPfam(afp, NULL, cfg->abc, -1, NULL, NULL, &(msaA[ai]), &nseq_cur, &alen_cur, &ngs_cur, &maxname_cur, &maxgf_cur, &maxgc_cur, &maxgr_cur, NULL, NULL, NULL, NULL, NULL);
    if      (status == eslEFORMAT) cm_Fail("Rereading alignment %d for merging, parse error:\n%s\n", ai+1, afp->errbuf);
    else if (status == eslEINVAL)  cm_Fail("Rereading alignment %d for merging, parse error:\n%s\n", ai+1, afp->errbuf);
    else if (status != eslOK)      cm_Fail("Rereading alignment %d for merging, parse error:\n%s\n", ai+1, afp->errbuf);

    msaA[ai]->abc = cfg->abc; 
    if(msaA[ai]->rf == NULL) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "When rereading alignment %d for merging, no RF annotation found.", ai+1);
    cur_clen = 0;
    for(apos = 0; apos < (int) alen_cur; apos++) { 
      if((! esl_abc_CIsGap(msaA[ai]->abc, msaA[ai]->rf[apos])) && (! esl_abc_CIsMissing(msaA[ai]->abc, msaA[ai]->rf[apos]))) cur_clen++;
    }
    if(cur_clen != cm->clen) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "When rereading alignment %d for merging, consensus length wrong (%d, when %d was expected)", ai, cur_clen, cm->clen);
    maxname = ESL_MAX(maxname, maxname_cur); 
    maxgf   = ESL_MAX(maxgf, maxgf_cur); 
    maxgc   = ESL_MAX(maxgc, maxgc_cur); 
    maxgr   = ESL_MAX(maxgr, maxgr_cur); 
    msaA[ai]->alen = alen_cur;
    alenA[ai]      = alen_cur; /* to remember total width of aln to expect in second pass */
    nseq_tot += nseq_cur;
    if(ngs_cur > 0) gs_exists = TRUE; 
      
    /* determine max number inserts and ELs between each position */
    update_maxins_and_maxel(msaA[ai], cm->clen, msaA[ai]->alen, maxins, maxel);
  }
  /* final check, make sure we've read all msas from the file, we should have, we only printed cfg->nali */
  status = esl_msa_ReadInfoPfam(afp, NULL, cfg->abc, -1, NULL, NULL, NULL, 
				      NULL, NULL, NULL, NULL, NULL, NULL, NULL, 
				      NULL, NULL, NULL, NULL, NULL);
  if(status != eslEOF) ESL_FAIL(status, errbuf, "More alignments in temp file than expected.");
  esl_msafile_Close(afp);
  
  /*********************************************
   *  Merge all alignments into the merged MSA *
   *********************************************/

  /* We allocate space for all sequences, but leave sequences as NULL
   * (nseq = -1).  We didn't store the sequences on the first pass
   * through the alignment files, and we'll never allocate space for
   * the sequences in fmsa, we'll just output them as we reread them
   * on another pass through the individual alignments. If we read >=
   * 1 GS line in any of the temporary alignments, we need to do an
   * additional pass through them, outputting only GS data. Then, in a
   * final (3rd) pass we'll output aligned data.
   */     
  fmsa = esl_msa_Create(nseq_tot, -1); 
  alen_fmsa = cm->clen + esl_vec_ISum(maxins, (cm->clen+1)); 

  /* if there was any GS annotation in any of the individual alignments,
   * do second pass through alignment files, outputting GS annotation as we go. */
  if(gs_exists) { 
    if((status = esl_msafile_Open(tmpfile, eslMSAFILE_PFAM, NULL, &afp)) != eslOK) cm_Fail("unable to open temp file %s for reading on second pass", tmpfile);
    for(ai = 0; ai < cfg->nali; ai++) { 
      regurg_header = (ai == 0) ? TRUE : FALSE;
      regurg_gf     = ((ofp != cfg->regressfp) && (ai == 0)) ? TRUE : FALSE;
      status = esl_msa_RegurgitatePfam(afp, ofp, 
				       maxname, maxgf, maxgc, maxgr, /* max width of a seq name, gf tag, gc tag, gr tag */
				       regurg_header, /* regurgitate stockholm header ? */
				       FALSE,         /* regurgitate // trailer ? */
				       regurg_header, /* regurgitate blank lines */
				       regurg_header, /* regurgitate comments */
				       regurg_gf,     /* regurgitate GF ? */
				       TRUE,          /* regurgitate GS ? */
				       FALSE,         /* regurgitate GC ? */
				       FALSE,         /* regurgitate GR ? */
				       FALSE,         /* regurgitate aseq ? */
				       NULL,          /* output all seqs, not just those stored in a keyhash */
				       NULL,          /* output all seqs, don't skip those listed in a keyhash */
				       NULL,          /* useme,  irrelevant, we're only outputting GS */
				       NULL,          /* add2me, irrelevant, we're only outputting GS */
				       alenA[ai], /* alignment length, as we read it in first pass (inserts may have been removed since then) */
				       '.', NULL, NULL);
      if(status == eslEOF) cm_Fail("Second pass, error out of temp alignments too soon, when trying to read aln %d", ai);
      if(status != eslOK)  cm_Fail("Second pass, error reading temp alignment %d %s", ai, afp->errbuf); 
      fflush(ofp);
    }
    esl_msafile_Close(afp);
    fprintf(ofp, "\n"); /* a single blank line to separate GS annotation from aligned data */
  }
  /* do another (either second or third) pass through alignment files, outputting aligned sequence data (and GR) as we go */

  if((status = esl_msafile_Open(tmpfile, eslMSAFILE_PFAM, NULL, &afp)) != eslOK) cm_Fail("unable to open temp file %s for reading on second (or third) pass", tmpfile);

  for(ai = 0; ai < cfg->nali; ai++) { 
    /* determine how many all gap columns to insert after each alignment position
     * of the temporary msa when copying it to the merged msa */
    if((status = determine_gap_columns_to_add(msaA[ai], maxins, maxel, cm->clen, &(ngap_insA), &(ngap_elA), &(ngap_eitherA), errbuf)) != eslOK) 
      cm_Fail("error determining number of all gap columns to add to temp alignment %d", ai);
    regurg_header = ((! gs_exists) && (ai == 0)) ? TRUE : FALSE;
    regurg_gf     = ((ofp != cfg->regressfp) && (! gs_exists) && (ai == 0)) ? TRUE : FALSE;

    status = esl_msa_RegurgitatePfam(afp, ofp,
				     maxname, maxgf, maxgc, maxgr, /* max width of a seq name, gf tag, gc tag, gr tag */
				     regurg_header,  /* regurgitate stockholm header ? */
				     FALSE,          /* regurgitate // trailer ? */
				     regurg_header,  /* regurgitate blank lines */
				     regurg_header,  /* regurgitate comments */
				     regurg_gf,      /* regurgitate GF ? */
				     FALSE,          /* regurgitate GS ? */
				     FALSE,          /* regurgitate GC ? */
				     TRUE,           /* regurgitate GR ? */
				     TRUE,           /* regurgitate aseq ? */
				     NULL,           /* output all seqs, not just those stored in a keyhash */
				     NULL,           /* output all seqs, don't skip those stored in a keyhash */
				     NULL,           /* useme, not nec b/c we want to keep all columns */
				     ngap_eitherA,   /* number of all gap columns to add after each apos */
				     alenA[ai],      /* alignment length, as we read it in first pass, not strictly necessary */
				     '.', NULL, NULL);
    if(status == eslEOF) cm_Fail("Second pass, error out of alignments too soon, when trying to read temp aln %d", ai);
    if(status != eslOK)  cm_Fail("Second pass, error reading temp alignment %d: %s", ai, afp->errbuf); 
    if(ai == 0) { 
      /* create the GC SS_cons and GC RF to print from the first alignment,
       * we use the first alignment b/c this is the one potentially with rewritten pknots
       * from --withpknots.
       */
      inflate_gc_with_gaps_and_els(ofp, msaA[ai], ngap_insA, ngap_elA, &ss_cons2print, &rf2print);
    }
    free(ngap_insA);
    free(ngap_elA);
    free(ngap_eitherA);
    
    esl_msa_Destroy(msaA[ai]);
    msaA[ai] = NULL;
    fflush(ofp);
  }
  /* output SS_cons and RF */
  margin = maxname+1;
  if (maxgc > 0 && maxgc+6         > margin) margin = maxgc+6;
  if (maxgr > 0 && maxname+maxgr+7 > margin) margin = maxname+maxgr+7; 
  fprintf(ofp, "#=GC %-*s %s\n", margin-6, "SS_cons", ss_cons2print);
  fprintf(ofp, "#=GC %-*s %s\n", margin-6, "RF", rf2print);
  fprintf(ofp, "//\n");

  esl_msafile_Close(afp);

  if(ss_cons2print != NULL) free(ss_cons2print);
  if(rf2print != NULL) free(rf2print);
  if(alenA != NULL)  free(alenA);
  if(msaA != NULL)   free(msaA);
  if(maxins != NULL) free(maxins);
  if(maxel != NULL)  free(maxel);
  if(fmsa != NULL)   esl_msa_Destroy(fmsa);
  return eslOK;

 ERROR: 
  esl_fatal("Out of memory. Reformat to Pfam with esl-reformat and try esl-alimerge --savemem.");
  return eslEMEM; /*NEVERREACHED*/
}

/* Function: update_maxins_and_maxel
 * Date:     EPN, Sun Nov 22 09:40:48 2009
 * 
 * Update maxins[] and maxel[], arrays that keeps track of the
 * max number of inserted ('.' gap #=GC RF) columns and inserted EL
 * emissions ('~' gap #=GC RF) columns before each cpos (consensus
 * (non-gap #=GC RF) column)
 *
 * Consensus columns are index [0..cpos..clen].
 * 
 * max{ins,el}[0]      is number of {IL/IR inserts, EL inserts} before 1st cpos.
 * max{ins,el}[clen-1] is number of {IL/IR inserts, EL inserts} before final cpos.
 * max{ins,el}[clen]   is number of {IL/IR inserts, EL inserts} after  final cpos.
 * 
 * Caller has already checked that msa->rf != NULL
 * and its non-gap length is clen. If we find either
 * of these is not true, we die (but this shouldn't happen).
 * 
 * Returns: void.
 */
void
update_maxins_and_maxel(ESL_MSA *msa, int clen, int64_t alen, int *maxins, int *maxel) 
{
  int apos;
  int cpos = 0;
  int nins = 0;
  int nel = 0;

  for(apos = 0; apos < alen; apos++) { 
    if(esl_abc_CIsGap(msa->abc, msa->rf[apos])) { 
      nins++;
    }
    else if (esl_abc_CIsMissing(msa->abc, msa->rf[apos])) { 
      nel++;
    }
    else {
      maxins[cpos] = ESL_MAX(maxins[cpos], nins);
      maxel[cpos]  = ESL_MAX(maxel[cpos], nel);
      cpos++;
      nins = 0;
      nel = 0;
    }
  }
      
  /* update final value, max{ins,el}[clen+1], the number of inserts
   * after the final consensus position */
  maxins[cpos] = ESL_MAX(maxins[cpos], nins);
  maxel[cpos]  = ESL_MAX(maxel[cpos], nel);
  if(cpos != clen) cm_Fail("Unexpected error in update_maxins_and_maxel(), expected clen (%d) not equal to actual clen (%d).\n", clen, cpos);

  return;
}

/* determine_gap_columns_to_add
 *                   
 * Given <maxins> and <maxel>, two arrays of the number of gap RF
 * (inserts) positions and '~' RF (EL inserts) after each non-gap RF 
 * (consensus) position in the eventual final merged alignment, 
 * calculate how many inserts and missing data inserts
 * we need to add at each position of <msa> to expand it out to the 
 * appropriate size of the eventual merged alignment.
 * 
 * max{ins,el}[0]      is number of inserts,ELs before 1st cpos in merged aln
 * max{ins,el}[cpos]   is number of inserts,ELs after  final cpos in merged aln
 *                             for cpos = 1..clen 
 * clen is the number of non-gap RF positions in msa (and in eventual merged msa).             
 * 
 * We allocate fill and return ret_ngap_insA[0..msa->alen], ret_ngap_elA[0..msa->alen], 
 * and ret_ngap_eitherA[0..msa->alen] here.
 *
 * ret_n{ins,el}gapA[0]      is number of inserts,ELs to add before 1st position of msa 
 * ret_n{ins,el}gapA[apos]   is number of inserts,ELs to add after alignment position apos
 *                             for apos = 1..msa->alen
 * 
 * ret_ngap_eitherA[apos] = ngap_insA[apos] + ngap_elA[apos]
 * 
 * This is similar to the esl_msa.c helper function of the same name,
 * but that function does not bother with missing data '~'.
 * 
 * Returns eslOK on success.
 *         eslEMEM on memory alloaction error 
 *         eslERANGE if a value exceeds what we expected (based on earlier
 *                   checks before this function was entered).
 *         if !eslOK, errbuf if filled.
 */
int
determine_gap_columns_to_add(ESL_MSA *msa, int *maxins, int *maxel, int clen, int **ret_ngap_insA, int **ret_ngap_elA, int **ret_ngap_eitherA, char *errbuf)
{
  int status;
  int apos;
  int prv_cpos = 0;  /* alignment position corresponding to consensus position cpos-1 */
  int cpos = 0;
  int nins = 0;
  int nel = 0;
  int *ngap_insA = NULL;
  int *ngap_elA = NULL;
  int *ngap_eitherA = NULL;
  int insert_before_el_flag = FALSE; /* set to TRUE for a cpos if we observe an insert column 5' of a missing data column between cpos-1 and cpos */
  int el_before_insert_flag = FALSE; /* set to TRUE for a cpos if we observe a missing data column 5' of an insert column between cpos-1 and cpos */

  /* contract check */
  if(maxel[0] != 0) ESL_FAIL(eslEINVAL, errbuf, "missing characters exist prior to first cpos, this shouldn't happen.\n");
  if(msa->ss_cons == NULL) ESL_FAIL(eslEINVAL, errbuf, "MSA's SS_cons is null in determine_gap_columns_to_add.\n");

  ESL_ALLOC(ngap_insA, sizeof(int) * (msa->alen+1));
  ESL_ALLOC(ngap_elA, sizeof(int) * (msa->alen+1));
  ESL_ALLOC(ngap_eitherA, sizeof(int) * (msa->alen+1));
  esl_vec_ISet(ngap_insA, (msa->alen+1), 0);
  esl_vec_ISet(ngap_elA, (msa->alen+1), 0);
  esl_vec_ISet(ngap_eitherA, (msa->alen+1), 0);
  
  for(apos = 0; apos < msa->alen; apos++) { 
    if(esl_abc_CIsMissing(msa->abc, msa->rf[apos])) { 
      nel++;
      if(nins > 0) insert_before_el_flag = TRUE;
    }
    else if(esl_abc_CIsGap(msa->abc, msa->rf[apos])) { 
      nins++;
      if(nel > 0) el_before_insert_flag = TRUE;
    }
    else { /* a consensus position */
      /* a few sanity checks */
      if(nins  > maxins[cpos])  ESL_FAIL(eslEINCONCEIVABLE, errbuf, "%d inserts before cpos %d greater than max expected (%d).\n", nins, cpos, maxins[cpos]); 
      if(nel > maxel[cpos]) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "%d EL inserts before cpos %d greater than max expected (%d).\n", nel, cpos, maxel[cpos]); 
      if(insert_before_el_flag && el_before_insert_flag) ESL_XFAIL(eslEINVAL, errbuf, "found inserts then missing chars '~' then inserts in RF before cpos %d.\n", cpos);

      if (cpos == 0) { 
	if(nel != 0) ESL_FAIL(eslEINVAL, errbuf, "found missing chars prior to first cpos, shouldn't happen\n");
	ngap_insA[prv_cpos]  = maxins[cpos] - nins; /* inserts before first position: flush right (so add all-gap columns after leftmost column) */
	/* we already checked that maxel[0] is 0 during contract check above */
      }
      else {
	/* Determine where to place inserts and/or missing data.
	 * Handle each of 4 possibilities separately, note that if 
	 * maxins[cpos] == 0 then nins == 0, and if maxel[cpos] == 0 then nel == 0 (see sanity check above). 
	 */
	if(maxins[cpos] >  0 && maxel[cpos] == 0) { /* most common case */
	  ngap_insA[prv_cpos + 1 + (nins/2)] = maxins[cpos] - nins; /* internal cpos: split */
	}
	else if(maxins[cpos] == 0 && maxel[cpos] > 0) { 
	  ngap_elA[prv_cpos + 1 + (nel/2)] = maxel[cpos] - nel; /* internal cpos: split */
	}
	else if(maxins[cpos] >  0 && maxel[cpos] > 0) { 
	  /* We have to determine 5'->3' order: it could be inserts 5'
	   * (before) missing data or missing data 5' (before)
	   * inserts.  Rule for this is convoluted but we can
	   * determine which order based on the consensus structure:
	   * ff cpos-1 and cpos form a base pair (SS_cons[cpos-1] ==
	   * '<' && SS_cons[cpos] == '>'), then we do inserts 5' of
	   * missing data, else we do missing data 5' of inserts.  
	   * 
	   * The reason is b/c all ELs (missing data) occur at the end
	   * of a stem (prior to END_nd's left consensus position).
	   * Inserts emitted by states immediately prior to an END_nd
	   * are involved in the lone case of ambiguity in the CM
	   * grammar because it indicates an insert position to which
	   * two insert states (one IL and one IR) emit. This
	   * ambiguity is dealt with by 'detaching' (see
	   * cm_modelmaker.c:cm_find_and_detach_dual_inserts()) insert
	   * states with index equal to the state index of an END_E
	   * state *minus 1*. These detached states are always IL
	   * states that would emit inserts 5' of any ELs (missing
	   * data) (e.g. a MATL_IL just before an END_E) *except* when
	   * a MATP_nd is immediately prior to an END_E node, in which
	   * case the detached state is an IR state (always MATP_IR
	   * just before END_E). Because these states are detached we
	   * know that the *other* state that inserts at that position
	   * must be used.  Either an IR (ROOT_IR, MATR_IR or MATP_IR)
	   * in the former case or a MATP_IL in the latter. The IR
	   * will insert 3' of an EL, and the MATP_IL will insert 5'
	   * of an EL. We can determine the MATP_IL cases by asking if
	   * cpos-1 and cpos form a basepair in the SS_cons, b/c
	   * that's the only case in which a MATP node exists
	   * immediately before an END_nd. Note that this basepair
	   * will always be a '<' '>' basepair (as opposed to a '('
	   * ')' or '{' '}' for example) b/c it is the lowest nesting
	   * order (enclosing basepair of stem. (Told you it was
	   * convoluted.)
	   */
	  if((msa->ss_cons[prv_cpos] == '<') && (msa->ss_cons[apos] == '>')) { /* special case */
	    /* inserts should be 5' of missing data (ELs) */
	    ngap_insA[prv_cpos + 1        + (nins/2)] = maxins[cpos] - nins; /* internal cpos: split */
	    ngap_elA[prv_cpos  + 1 + nins + (nel/2)] = maxel[cpos] - nel; /* internal cpos: split */
	  }
	  else { 
	    /* missing data (ELs) should be 5' of inserts */
	    ngap_elA[prv_cpos + 1 + (nel/2)] = maxel[cpos] - nel; /* internal cpos: split */
	    ngap_insA[prv_cpos + 1 + nel + (nins/2)] = maxins[cpos] - nins; /* internal cpos: split */
	  }
	}
	/* final case is if (maxins[cpos] == 0 && maxel[cpos] == 0) 
	 * in this case we do nothing. 
	 */
      }

      cpos++;
      prv_cpos = apos;
      nins = 0;
      nel = 0;
      el_before_insert_flag = FALSE;
      insert_before_el_flag = FALSE;
    }
  }
  cpos = clen;
  apos = msa->alen;
  /* deal with inserts after final consensus position, we could have both inserts and missing data, but only in case where 
   * CM has 0 bps (model is ROOT, MATL, MATL, ..., MATL, END) so we know that MATL_IL has been detached, ROOT_IR emits inserts
   * after final cpos, and thus missing data (ELs) should be 5' of inserts 
   */
  if(maxins[cpos] > 0 && maxel[cpos] == 0) { /* most common case */
    ngap_insA[apos] = maxins[cpos] - nins; /* flush left inserts (no ELs) */
  }
  else if(maxins[cpos] == 0 && maxel[cpos] > 0) { 
    ngap_elA[apos] = maxel[cpos] - nel; /* flush left ELs (no inserts) */
  }
  else if(maxins[cpos] > 0 && maxel[cpos] > 0) { 
    if((msa->ss_cons[prv_cpos] == '<') && (msa->ss_cons[apos] == '>')) ESL_FAIL(eslEINVAL, errbuf, "final position needs ELs and inserts after it, and is right half of a base pair, this shouldn't happen");
    /* missing data (ELs) should be 5' of inserts */
    ngap_elA[apos]      = maxel[cpos] - nel; /* flush left */
    ngap_insA[apos+nel] = maxins[cpos] - nins; /* flush left, after ELs */
  }

  /* determine ngap_eitherA[], the number of gaps due to either inserts or missing data after each apos */
  for(apos = 0; apos <= msa->alen; apos++) { 
    ngap_eitherA[apos] = ngap_insA[apos] + ngap_elA[apos];
  }

  /* validate that clen is what it should be */
  if(cpos != clen) { 
    if(ngap_insA != NULL) free(ngap_insA);
    if(ngap_elA != NULL) free(ngap_elA);
    if(ngap_eitherA != NULL) free(ngap_eitherA);
    ESL_FAIL(eslEINCONCEIVABLE, errbuf, "consensus length (%d) is not the expected length (%d).", cpos, clen);
  }

  *ret_ngap_insA  = ngap_insA;
  *ret_ngap_elA = ngap_elA;
  *ret_ngap_eitherA = ngap_eitherA;

  return eslOK;

 ERROR: 
  if(ngap_insA  != NULL) free(ngap_insA);
  if(ngap_elA != NULL) free(ngap_elA);
  if(ngap_eitherA != NULL) free(ngap_eitherA);
  ESL_FAIL(status, errbuf, "Memory allocation error.");
  return status; /*NEVERREACHED*/
}

/* inflate_gc_with_gaps_and_els
 *                   
 * Given an MSA and two arrays specifying the number of inserts '.'
 * and EL inserts '~' to add after each position, create the 
 * SS_cons and RF strings to output for the merged alignment
 * after adding the gaps and missing data symbols ('~') and
 * return them in ret_ss_cons2print and ret_rf2print. Caller
 * is responsible for freeing them.
 *
 * Returns void. If something unexpected occurs, including an 
 * allocation error, we die here and print error message.
 */
void
inflate_gc_with_gaps_and_els(FILE *ofp, ESL_MSA *msa, int *ngap_insA, int *ngap_elA, char **ret_ss_cons2print, char **ret_rf2print) 
{
  int status;
  int apos  = 0;
  int apos2print  = 0;
  int i, j;
  int el_before_ins = FALSE;
  int prv_cpos = -1;
  int alen2print = 0;
  char *rf2print;
  char *ss_cons2print;

  alen2print = msa->alen + esl_vec_ISum(ngap_insA, msa->alen+1) + esl_vec_ISum(ngap_elA, msa->alen+1);
  ESL_ALLOC(rf2print,      sizeof(char) * (alen2print+1));
  ESL_ALLOC(ss_cons2print, sizeof(char) * (alen2print+1));
  rf2print[alen2print] = '\0';
  ss_cons2print[alen2print] = '\0';

  if(msa->ss_cons == NULL) cm_Fail("Error trying to add inserts to SS_cons, but unexpectedly it doesn't exist.");
  if(msa->rf      == NULL) cm_Fail("Error trying to add inserts to SS_cons, but unexpectedly it doesn't exist.");
  for(apos = 0; apos <= msa->alen; apos++) { 
    el_before_ins = TRUE; /* until proven otherwise */
    if(ngap_insA[apos] > 0 && ngap_elA[apos] > 0) { 
      /* Rare case: we need to add at least one '.' and '~'. We need
       * to determine which to add first. The only possible way we do
       * '.'s and then '~'s is if we're inserting between the two
       * match positions of a MATP node and the following node is an
       * END (search for 'convoluted' in comments in
       * determine_gap_columns_to_add() for explanation of why).
       * This will only occur if previously seen consensus SS_cons 
       * char was a '<' and next consensus SS_cons char is a '>'.
       * First, check if prev consensus SS_cons char is '<',
       * if it is check if following one is a '>' by searching for
       * it. This is inefficient but should be rare so it's okay 
       * (we'll only do this at most once per END node in model).
       */
      if((prv_cpos == -1) || (msa->ss_cons[prv_cpos] != '<')) el_before_ins = TRUE; /* normal case */
      else { /* prv_cpos is '<', check if next one is a '>' */
	j = i+1; 
	while(j < msa->alen && msa->ss_cons[j] == '.') j++;
	if((j < msa->alen) && (msa->ss_cons[j] == '>')) el_before_ins = FALSE;
      }
    }
    if(el_before_ins) { 
      for(i = 0; i < ngap_elA[apos]; i++) { 
	rf2print[apos2print] = '~';
	ss_cons2print[apos2print++] = '~';
      }
      for(i = 0; i < ngap_insA[apos]; i++) { 
	rf2print[apos2print] = '.';
	ss_cons2print[apos2print++] = '.';
      }
    } 
    else { /* insert before els, rare, only occurs if el_before_ins set to FALSE above */
      for(i = 0; i < ngap_insA[apos]; i++) { 
	rf2print[apos2print] = '.';
	ss_cons2print[apos2print++] = '.';
      }
      for(i = 0; i < ngap_elA[apos]; i++) { 
	rf2print[apos2print] = '~';
	ss_cons2print[apos2print++] = '~';
      }
    }
    if(apos < msa->alen) { 
      rf2print[apos2print]        = msa->rf[apos];
      ss_cons2print[apos2print++] = msa->ss_cons[apos];
      if(! esl_abc_CIsGap(msa->abc, msa->rf[apos])) prv_cpos = apos;
    }	
  }    

  *ret_ss_cons2print = ss_cons2print;
  *ret_rf2print = rf2print;

  return;

 ERROR:
  cm_Fail("Allocation error when creating final alignment RF and SS_cons.");
  return; /* NEVERREACHED */
}


#ifdef HAVE_MPI
/* determine_nseq_per_worker()
 * Given a CM, return the number of sequences we think we should send
 * to each worker (we don't know the number of sequences in the file).
 * The calculation is based on trying to get a worker to spend 
 * a specific amount of time MPI_WORKER_ALIGN_TARGET_SEC, a constant
 * from structs.h. 
 */
static int
determine_nseq_per_worker(const ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm, int *ret_nseq_worker)
{
  *ret_nseq_worker = esl_opt_GetInteger(go, "--nseq-mpi");
  return eslOK;
}

/* add_worker_seqs_to_master
 * Add results (parstrees or CP9 traces, and possibly postcodes) from a
 * worker's seqs_to_aln object to a master seqs_to_aln object.
 */
static int
add_worker_seqs_to_master(seqs_to_aln_t *master_seqs, seqs_to_aln_t *worker_seqs, int offset)
{
  int x;

  if(worker_seqs->sq != NULL) cm_Fail("add_worker_seqs_to_master(), worker_seqs->sq non-NULL.");
  if(master_seqs->nseq < (offset + worker_seqs->nseq)) cm_Fail("add_worker_seqs_to_master(), master->nseq: %d, offset %d, worker->nseq: %d\n", master_seqs->nseq, offset, worker_seqs->nseq);

  if(worker_seqs->tr != NULL) {
    if(master_seqs->tr == NULL) cm_Fail("add_worker_seqs_to_master(), worker returned parsetrees, master->tr is NULL.");
    for(x = offset; x < (offset + worker_seqs->nseq); x++) {
      assert(master_seqs->tr[x] == NULL); 
      master_seqs->tr[x] = worker_seqs->tr[(x-offset)];
    }
  }

  if(worker_seqs->cp9_tr != NULL) {
    if(master_seqs->cp9_tr == NULL) cm_Fail("add_worker_seqs_to_master(), worker returned cp9 traces, master->cp9_tr is NULL.");
    for(x = offset; x < (offset + worker_seqs->nseq); x++) {
      assert(master_seqs->cp9_tr[x] == NULL); 
      master_seqs->cp9_tr[x] = worker_seqs->cp9_tr[(x-offset)];
    }
  }

  if(worker_seqs->postcode != NULL) {
    if(master_seqs->postcode == NULL) cm_Fail("add_worker_seqs_to_master(), worker returned postcodes, master->postcode is NULL.");
    for(x = offset; x < (offset + worker_seqs->nseq); x++) {
      assert(master_seqs->postcode[x] == NULL); 
      master_seqs->postcode[x] = worker_seqs->postcode[(x-offset)];
    }
  }

  if(worker_seqs->sc != NULL) {
    if(master_seqs->sc == NULL) cm_Fail("add_worker_seqs_to_master(), worker returned scores, master->sc is NULL.");
    for(x = offset; x < (offset + worker_seqs->nseq); x++) {
      assert(!(NOT_IMPOSSIBLE(master_seqs->sc[x])));
      master_seqs->sc[x] = worker_seqs->sc[(x-offset)];
    }
  }

  if(worker_seqs->pp != NULL) {
    if(master_seqs->pp == NULL) cm_Fail("add_worker_seqs_to_master(), worker returned post probs, master->pp is NULL.");
    for(x = offset; x < (offset + worker_seqs->nseq); x++) {
      assert(!(NOT_IMPOSSIBLE(master_seqs->pp[x])));
      master_seqs->pp[x] = worker_seqs->pp[(x-offset)];
    }
  }

  if(worker_seqs->struct_sc != NULL) {
    if(master_seqs->struct_sc == NULL) cm_Fail("add_worker_seqs_to_master(), worker returned post probs, master->struct_sc is NULL.");
    for(x = offset; x < (offset + worker_seqs->nseq); x++) {
      assert(!(NOT_IMPOSSIBLE(master_seqs->struct_sc[x])));
      master_seqs->struct_sc[x] = worker_seqs->struct_sc[(x-offset)];
    }
  }

  return eslOK;
}
#endif /* of #ifdef HAVE_MPI */

/* serial_master_meta()
 * The serial version of cmalign in meta mode (-M enabled)
 * 1. Read all CMs in CM file
 * 2. Read all the alignments from the meta-cm training alignment file
 * 3. Validate the alignments and CMs make up a valid meta-CM
 * 4. Align all seqs to the major CM 
 * 5. Impose master alignment of each seq onto all other CMs, and determine implicit score
 * 6. Select winning CM, and realign to that CM.
 * 
 * A master can only return if it's successful. All errors are handled immediately and fatally with cm_Fail().
 */
static void
serial_master_meta(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int      status;
  char     errbuf[cmERRBUFSIZE];
  seqs_to_aln_t  *seqs_to_aln;  /* sequences to align, holds seqs, parsetrees, CP9 traces, postcodes */
  int    m,i;         /* counters */
  int    nalloc = 1;  /* number of CMs we've allocated */
  CM_t **cmlist;      /* [0..m..cfg->ncm-1] CM read from cmfile */
  void *tmp;          /* for ESL_RALLOC */
  int **toadd2min = NULL;/* [0..m..cfg->ncm][0..c..cmlist[m]->clen] = x, when morphing minor parsetrees to a major alignment, add x major consensus columns after minor consensus column c of CM m, 
		          *                                              these x major consensus columns DO NOT MAP to a consensus column in minor CM m */
  int *maj_train_a2c_map = NULL; /* [0..a..cfg->mali_msa[0]->alen] = x, alignment column a from the first training alignment (cfg->mali_msa[0]) maps to major consensus column x */
  ESL_MSA *majed_min_msa = NULL; /* temporary majorfied (inferred major) alignment from minor CM parsetrees */
  int *wcm;           /* [0..i..seqs_to_aln->nseq-1] = m, cmlist[m] is 'winning' (highest scoring) CM for seq i */
  int *maj_train_emitl; /* saved major alignment guidetree's emitl vector */
  int *maj_train_emitr; /* saved major alignment guidetree's emitr vector */
  seqs_to_aln_t *min_seqs_to_aln; /* temporary seqs_to_aln for each minor CM */
  int mct;                        /* number of winning seqs to align to each minor CM */
  int min_i, maj_i;               /* seq indices in min_seqs_to_aln and seqs_to_aln respectively */
  CMEmitMap_t *maj_emap;          /* emitmap for the major CM */
  Parsetree_t **majed_tr = NULL;  /* majorfied parsetrees */

  if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);
  if ((status = print_run_info (go, cfg, errbuf))  != eslOK) cm_Fail(errbuf);
  
  /* 1. Read all CMs in CM file */
  cfg->ncm = 0;
  ESL_ALLOC(cmlist, sizeof(CM_t *) * nalloc);
  status = eslOK;
  while (status == eslOK) { 
    status = CMFileRead(cfg->cmfp, errbuf, &(cfg->abc), &(cmlist[cfg->ncm]));
    if(status == eslOK) { 
      cfg->ncm++;
      if(cfg->ncm == nalloc) { 
	nalloc++; /* could be += n==5 or so */
	ESL_RALLOC(cmlist, tmp, sizeof(CM_t *) * nalloc);
      }	
      /*printf("cm %4d clen: %d\n", cfg->ncm, cmlist[cfg->ncm-1]->clen);*/
      if((status   = initialize_cm(go, cfg, errbuf, cmlist[cfg->ncm-1])) != eslOK) cm_Fail(errbuf);
      if(cfg->ncm == 1) print_cm_info (go, cfg, errbuf, cmlist[cfg->ncm-1]);
    }
  }
  if(status != eslEOF) cm_Fail(errbuf);

  /* 2. Read all the alignments from the -M <f> meta-cm training alignment file */
  status = eslOK;
  ESL_ALLOC(cfg->mali_msa, sizeof(ESL_MSA *) * cfg->ncm);
  ESL_MSA *tmp_msa; /* so we can ensure we have correct number of alignments in file, should be equal to number of CMs we just read */
  while (status == eslOK) { 
    status = esl_msa_Read(cfg->malifp, &tmp_msa);
    if(status == eslOK) { 
      if(cfg->mali_n >= cfg->ncm) cm_Fail("with -M, read %d CMs, but %s has > %d alignments.", cfg->ncm, esl_opt_GetString(go, "-M"), cfg->ncm);
      cfg->mali_msa[cfg->mali_n] = tmp_msa;
      cfg->mali_n++;
      /*printf("ali %4d alen: %" PRId64 "\n", cfg->mali_n, cfg->mali_msa[(cfg->mali_n-1)]->alen);*/
    }
    else { esl_msa_Destroy(tmp_msa); } 
  }
  if(cfg->mali_n != cfg->ncm) cm_Fail("with -M, read %d CMs, but %s has %d alignments.", cfg->ncm, esl_opt_GetString(go, "-M"), cfg->mali_n);

  /* 3. Validate the alignments and CMs make up a valid meta-CM */
  if((status = validate_meta(go, cfg, errbuf, cmlist, &toadd2min, &maj_train_a2c_map)) != eslOK) cm_Fail(errbuf);

  /* initialize the flags/options/params and configuration of the CMs */
  for(m = 0; m < cfg->ncm; m++) if((status   = initialize_cm(go, cfg, errbuf, cmlist[m])) != eslOK)    cm_Fail(errbuf);

  /* 4. Align all seqs to the major CM */
  /* read in all sequences */
  seqs_to_aln = CreateSeqsToAln(100, FALSE);
  if((status = ReadSeqsToAln(cfg->abc, cfg->sqfp, 0, seqs_to_aln, FALSE)) != eslEOF) cm_Fail("Error reading sqfile: %s\n", cfg->sqfile);
  /* align sequences to CM to get parsetrees  */
  if ((status = process_workunit(go, cfg, errbuf, cmlist[0], seqs_to_aln)) != eslOK) cm_Fail(errbuf);
  /* convert parsetrees to alignment */
  ESL_MSA *maj_target_msa;
  if((status = Parsetrees2Alignment(cmlist[0], errbuf, cfg->abc_out, seqs_to_aln->sq, NULL, seqs_to_aln->tr, NULL, seqs_to_aln->nseq, NULL, NULL, TRUE, esl_opt_GetBoolean(go, "--matchonly"), FALSE, &maj_target_msa)) != eslOK)
    cm_Fail("serial_master_meta(), error generating major alignment from parsetrees.");
  /* printf("\n"); status = esl_msa_Write(stdout, maj_target_msa, eslMSAFILE_STOCKHOLM); printf("\n"); */

  /* 5. Determine the implicit minor alignments defined by the major alignment, and determine implicit scores */
  if((status = major_alignment2minor_parsetrees(go, cfg, errbuf, cmlist, maj_target_msa, maj_train_a2c_map, &wcm)) != eslOK) cm_Fail(errbuf);

  /* 6. Realign each seq to it's winning minor CM */
  /* copy the training alignment master CM's guidetree emitl and emitr vectors, we'll overwrite them as we convert the minor
   * alignments to major coords */
  ESL_ALLOC(maj_train_emitl, sizeof(int) * cfg->mali_mtr[0]->n);
  ESL_ALLOC(maj_train_emitr, sizeof(int) * cfg->mali_mtr[0]->n);
  esl_vec_ICopy(cfg->mali_mtr[0]->emitl, cfg->mali_mtr[0]->n, maj_train_emitl);
  esl_vec_ICopy(cfg->mali_mtr[0]->emitr, cfg->mali_mtr[0]->n, maj_train_emitr);
  maj_emap = CreateEmitMap(cmlist[0]); /* used to convert major CM guidetree coords from training alignment coords to majorfied minor CM alignment coords */

  /* free major traces to all seqs that didn't have the major CM as the winner (we don't have to realign those) */
  for(i = 0; i < maj_target_msa->nseq; i++) if(wcm[i] != 0) { FreeParsetree(seqs_to_aln->tr[i]); seqs_to_aln->tr[i] = NULL; }
  for(m = 1; m < cfg->ncm; m++) { 
    mct = 0;
    for(i = 0; i < maj_target_msa->nseq; i++) if(wcm[i] == m) mct++;
    if(mct == 0) continue;
    min_seqs_to_aln = CreateSeqsToAln(mct, FALSE); 
    for(i = 0; i < maj_target_msa->nseq; i++) { 
      if(wcm[i] == m) min_seqs_to_aln->sq[min_seqs_to_aln->nseq++] = seqs_to_aln->sq[i]; 
    }      
    /* align the sequences to the minor CM to get minor parsetrees */
    if ((status = process_workunit(go, cfg, errbuf, cmlist[m], min_seqs_to_aln)) != eslOK) cm_Fail(errbuf);
    /* minor parsetrees -> implicit major alignment (majorfied alignment) */
    if((status = Parsetrees2Alignment_Minor2Major(cmlist[m], cfg->abc_out, min_seqs_to_aln->sq, NULL, min_seqs_to_aln->tr, min_seqs_to_aln->nseq, TRUE, esl_opt_GetBoolean(go, "--matchonly"), toadd2min[m], &majed_min_msa)) != eslOK)
      cm_Fail("serial_master_meta(), error generating alignment from parsetrees to major CM.");
    /* status = esl_msa_Write(stdout, majed_min_msa, eslMSAFILE_STOCKHOLM); 
       DumpEmitMap(stdout, maj_emap, cmlist[0]);
      ParsetreeDump(stdout, cfg->mali_mtr[0], cmlist[0], maj_target_msa->ax[0], NULL, NULL); */

    /* majorfied alignment -> major parsetrees */
    if((status = majorfied_alignment2major_parsetrees(go, cfg, errbuf, cmlist[0], majed_min_msa, maj_emap, &majed_tr)) != eslOK) cm_Fail(errbuf);

    /* swap corresponding pointers of orig seqstoaln->tr to majed_tr */
    min_i = 0;
    for(maj_i = 0; maj_i < seqs_to_aln->nseq; maj_i++) { 
      if(wcm[maj_i] == m) seqs_to_aln->tr[maj_i] = majed_tr[min_i++];
    }
    free(majed_tr); /* only free the top level ptr, seqs_to_aln->tr now points to the actual parsetrees */
    esl_msa_Destroy(majed_min_msa);
    /* careful, don't free the sq, seqs_to_aln->sq still points to them, and we need them */
    FreePartialSeqsToAln(min_seqs_to_aln, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE);
    /*                                       sq,   tr,cp9_tr,post,   sc,   pp, struct_sc */
    free(min_seqs_to_aln->sq);
    free(min_seqs_to_aln);
  }
  /* output the alignment to the master CM */
  if ((status = output_result(go, cfg, FALSE, errbuf, cmlist[0], seqs_to_aln)) != eslOK) cm_Fail(errbuf);

  /* clean up */
  FreeSeqsToAln(seqs_to_aln);
  for(m = 0; m < cfg->ncm; m++) { 
    FreeCM(cmlist[m]);
    free(toadd2min[m]);
  }
  esl_msa_Destroy(maj_target_msa);
  free(cmlist);
  free(toadd2min);
  free(maj_train_a2c_map);
  free(wcm);
  free(maj_train_emitl);
  free(maj_train_emitr);
  FreeEmitMap(maj_emap);

  esl_sqfile_Position(cfg->sqfp, (off_t) 0); /* we may be aligning these seqs another CM */
  return;

 ERROR:
  cm_Fail("Memory allocation error in serial_master_meta().");
  return;
}

/* map_cpos_to_apos
 *                   
 * Given an MSA, determine the alignment position each
 * consensus position refers to. 
 */
static int map_cpos_to_apos(ESL_MSA *msa, int **ret_c2a_map, int *ret_clen)
{
  int status;
  int clen = 0;
  int *c2a_map = NULL;
  int cpos = 0;
  int apos = 0;
  /* contract check */
  if(msa->rf == NULL) { status = eslEINVAL; goto ERROR; }

  /* count consensus columns */
  for(apos = 1; apos <= msa->alen; apos++)
    if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) clen++;

  /* build map */
  ESL_ALLOC(c2a_map, sizeof(int) * (clen+1));
  c2a_map[0] = -1;
  for(apos = 1; apos <= msa->alen; apos++) 
    if(! esl_abc_CIsGap(msa->abc, msa->rf[(apos-1)])) c2a_map[++cpos] = apos;

  *ret_c2a_map = c2a_map;
  *ret_clen    = clen;
  return eslOK;

 ERROR:
  if(c2a_map != NULL) free(c2a_map);
  return status;
}

/* Function : Parsetrees2Alignment_Minor2Major()
 *
 * Purpose:   Creates a MSA from a set of parsetrees and a CM.
 * 
 * Args:     cm         - the CM the CP9 was built from, needed to get emitmap,
 *                        so we know where to put EL transitions
 *           abc        - alphabet to use to create the return MSA
 *           sq         - sequences, must be digitized (we check for it)
 *           wgt        - weights for seqs, NULL for none
 *           nseq       - number of sequences
 *           tr         - array of tracebacks
 *           do_full    - TRUE to always include all match columns in alignment
 *           do_matchonly - TRUE to ONLY include match columns
 *           masteradd  - [0..c..cm->clen] number of master consensus columns to add after (minor) consensus column c
 *           ret_msa    - MSA, alloc'ed/created here
 *
 * Return:   eslOK on succes, eslEMEM on memory error.
 *           MSA structure in ret_msa, caller responsible for freeing.
 *
 * Returns:   eslOK on success, eslEMEM on memory error, 
 *            Also ret_msa is filled with a new MSA.
 *
 */
static int
Parsetrees2Alignment_Minor2Major(CM_t *cm, const ESL_ALPHABET *abc, ESL_SQ **sq, float *wgt, 
				 Parsetree_t **tr, int nseq, int do_full, int do_matchonly, 
				 int *masteradd, ESL_MSA **ret_msa)
{
  char errbuf[eslERRBUFSIZE];

  /* Contract check. We allow the caller to specify the alphabet they want the 
   * resulting MSA in, but it has to make sense (see next few lines). */
  if(cm->abc->type == eslRNA)
    { 
      if(abc->type != eslRNA && abc->type != eslDNA)
	cm_Fail("ERROR in Parsetrees2Alignment_Minor2Major(), cm alphabet is RNA, but requested output alphabet is neither DNA nor RNA.");
    }
  else if(cm->abc->K != abc->K)
    cm_Fail("ERROR in Parsetrees2Alignment_Minor2Major(), cm alphabet size is %d, but requested output alphabet size is %d.", cm->abc->K, abc->K);

  int          status = eslOK;/* easel status flag */
  ESL_MSA     *msa   = NULL; /* multiple sequence alignment */
  CMEmitMap_t *emap  = NULL; /* consensus emit map for the CM */
  int          i;            /* counter over traces */
  int          v, nd;        /* state, node indices */
  int          cpos;         /* counter over consensus positions (0)1..clen */
  int         *matuse= NULL; /* TRUE if we need a cpos in mult alignment */
  int         *iluse = NULL; /* # of IL insertions after a cpos for 1 trace */
  int         *eluse = NULL; /* # of EL insertions after a cpos for 1 trace */
  int         *iruse = NULL; /* # of IR insertions after a cpos for 1 trace */
  int         *maxil = NULL; /* max # of IL insertions after a cpos */
  int         *maxel = NULL; /* max # of EL insertions after a cpos */
  int         *maxir = NULL; /* max # of IR insertions after a cpos */
  int	      *matmap= NULL; /* apos corresponding to a cpos */
  int         *ilmap = NULL; /* first apos for an IL following a cpos */
  int         *elmap = NULL; /* first apos for an EL following a cpos */
  int         *irmap = NULL; /* first apos for an IR following a cpos */
  int         *mastermap = NULL; /* first apos for a master cpos not in the minor CM emap */
  int          alen;	     /* length of msa in columns */
  int          apos;	     /* position in an aligned sequence in MSA */
  int          rpos;	     /* position in an unaligned sequence in dsq */
  int          tpos;         /* position in a parsetree */
  int          el_len;	     /* length of an EL insertion in residues */
  CMConsensus_t *con = NULL; /* consensus information for the CM */
  int          prvnd = -1;   /* keeps track of previous node for EL */
  int          nins;         /* insert counter used for splitting inserts */

  emap = CreateEmitMap(cm);

  ESL_ALLOC(matuse, sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(iluse,  sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(eluse,  sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(iruse,  sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(maxil,  sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(maxel,  sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(maxir,  sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(matmap, sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(ilmap,  sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(elmap,  sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(irmap,  sizeof(int)*(emap->clen+1));   
  ESL_ALLOC(mastermap, sizeof(int)*(emap->clen+1));   
  
  for (cpos = 0; cpos <= emap->clen; cpos++) 
    {
      if(!do_full || cpos == 0)
	matuse[cpos] = 0;
      else
	matuse[cpos] = 1;
      maxil[cpos] = maxel[cpos] = maxir[cpos] = 0;
      ilmap[cpos] = elmap[cpos] = irmap[cpos] = mastermap[cpos] = 0;
    }

  /* Look at all the traces; find maximum length of
   * insert needed at each of the clen+1 possible gap
   * points. (There are three types of insert, IL/EL/IR.)
   * Also find whether we don't need some of the match
   * (consensus) columns.
   */
  for (i = 0; i < nseq; i++) 
    {
      for (cpos = 0; cpos <= emap->clen; cpos++) 
	iluse[cpos] = eluse[cpos] = iruse[cpos] = 0;

      for (tpos = 0; tpos < tr[i]->n; tpos++)
	{
	  v  = tr[i]->state[tpos];
	  if (cm->sttype[v] == EL_st) nd = prvnd;
	  else                        nd = cm->ndidx[v];
	  
	  switch (cm->sttype[v]) {
	  case MP_st: 
	    matuse[emap->lpos[nd]] = 1;
	    matuse[emap->rpos[nd]] = 1;
	    break;
	  case ML_st:
	    matuse[emap->lpos[nd]] = 1;
	    break;
	  case MR_st:
	    matuse[emap->rpos[nd]] = 1;
	    break;
	  case IL_st:
	    iluse[emap->lpos[nd]]++;
	    break;
	  case IR_st:		
            /* remember, convention on rpos is that IR precedes this
             * cpos. Make it after the previous cpos, hence the -1. 
	     */
	    iruse[emap->rpos[nd]-1]++;
	    break;
	  case EL_st:
	    el_len = tr[i]->emitr[tpos] - tr[i]->emitl[tpos] + 1;
	    eluse[emap->epos[nd]] = el_len;
            /* not possible to have >1 EL in same place; could assert this */
	    break;
	  }

	  prvnd = nd;
	} /* end looking at trace i */

      for (cpos = 0; cpos <= emap->clen; cpos++) 
	{
	  if (iluse[cpos] > maxil[cpos]) maxil[cpos] = iluse[cpos];
	  if (eluse[cpos] > maxel[cpos]) maxel[cpos] = eluse[cpos];
	  if (iruse[cpos] > maxir[cpos]) maxir[cpos] = iruse[cpos];
	}
    } /* end calculating lengths used by all traces */
  
  /* Now we can calculate the total length of the multiple alignment, alen;
   * and the maps ilmap, elmap, and irmap that turn a cpos into an apos
   * in the multiple alignment: e.g. for an IL that follows consensus position
   * cpos, put it at or after apos = ilmap[cpos] in aseq[][].
   * IR's are filled in backwards (3'->5') and rightflushed.
   */
  alen = 0;
  for (cpos = 0; cpos <= emap->clen; cpos++)
    {
      if (matuse[cpos]) {
	matmap[cpos] = alen; 
	alen++;
      } else 
	matmap[cpos] = -1;

      ilmap[cpos]     = alen; alen += maxil[cpos];
      elmap[cpos]     = alen; alen += maxel[cpos];
      mastermap[cpos] = alen; alen += masteradd[cpos];
      alen += maxir[cpos]; irmap[cpos] = alen-1; 
    }

  /* We're getting closer.
   * Now we can allocate for the MSA.
   */
  msa = esl_msa_Create(nseq, alen);
  if(msa == NULL) goto ERROR;
  msa->nseq = nseq;
  msa->alen = alen;
  msa->abc  = (ESL_ALPHABET *) abc;

  for (i = 0; i < nseq; i++)
    {
      /* Contract check */
      if(sq[i]->dsq == NULL) cm_Fail("ERROR in Parsetrees2Alignment_Minor2Major(), sq %d is not digitized.\n", i);

      /* Initialize the aseq with all pads '.' (in insert cols) 
       * and deletes '-' (in match cols).
       */
      for (apos = 0; apos < alen; apos++)
	msa->aseq[i][apos] = '.';
      for (cpos = 0; cpos <= emap->clen; cpos++) { 
	if (matmap[cpos] != -1) msa->aseq[i][matmap[cpos]] = '-';
	for(apos = mastermap[cpos]; apos < mastermap[cpos] + masteradd[cpos]; apos++) msa->aseq[i][apos] = '-'; /* add master consensus column deletes, not in this cm's consensus */
      }
      msa->aseq[i][alen] = '\0';

      /* Traverse this guy's trace, and place all his
       * emitted residues.
       */
      for (cpos = 0; cpos <= emap->clen; cpos++)
	iluse[cpos] = iruse[cpos] = 0;

      for (tpos = 0; tpos < tr[i]->n; tpos++) 
	{
	  v  = tr[i]->state[tpos];
	  if (cm->sttype[v] == EL_st) nd = prvnd;
	  else                        nd = cm->ndidx[v];

	  switch (cm->sttype[v]) {
	  case MP_st:
	    cpos = emap->lpos[nd];
	    apos = matmap[cpos];
	    rpos = tr[i]->emitl[tpos];
	    msa->aseq[i][apos] = abc->sym[sq[i]->dsq[rpos]];

	    cpos = emap->rpos[nd];
	    apos = matmap[cpos];
	    rpos = tr[i]->emitr[tpos];
	    msa->aseq[i][apos] = abc->sym[sq[i]->dsq[rpos]];
	    break;
	    
	  case ML_st:
	    cpos = emap->lpos[nd];
	    apos = matmap[cpos];
	    rpos = tr[i]->emitl[tpos];
	    msa->aseq[i][apos] = abc->sym[sq[i]->dsq[rpos]];
	    break;

	  case MR_st:
	    cpos = emap->rpos[nd];
	    apos = matmap[cpos];
	    rpos = tr[i]->emitr[tpos];
	    msa->aseq[i][apos] = abc->sym[sq[i]->dsq[rpos]];
	    break;

	  case IL_st:
	    cpos = emap->lpos[nd];
	    apos = ilmap[cpos] + iluse[cpos];
	    rpos = tr[i]->emitl[tpos];
	    msa->aseq[i][apos] = tolower((int) abc->sym[sq[i]->dsq[rpos]]);
	    iluse[cpos]++;
	    break;

	  case EL_st: 
            /* we can assert eluse[cpos] always == 0 when we enter,
	     * because we can only have one EL insertion event per 
             * cpos. If we ever decide to regularize (split) insertions,
             * though, we'll want to calculate eluse in the rpos loop.
             */
	    cpos = emap->epos[nd]; 
	    apos = elmap[cpos]; 
	    for (rpos = tr[i]->emitl[tpos]; rpos <= tr[i]->emitr[tpos]; rpos++)
	      {
		msa->aseq[i][apos] = tolower((int) abc->sym[sq[i]->dsq[rpos]]);
		apos++;
	      }
	    break;

	  case IR_st: 
	    cpos = emap->rpos[nd]-1;  /* -1 converts to "following this one" */
	    apos = irmap[cpos] - iruse[cpos];  /* writing backwards, 3'->5' */
	    rpos = tr[i]->emitr[tpos];
	    msa->aseq[i][apos] = tolower((int) abc->sym[sq[i]->dsq[rpos]]);
	    iruse[cpos]++;
	    break;

	  case D_st:
	    if (cm->stid[v] == MATP_D || cm->stid[v] == MATL_D) 
	      {
		cpos = emap->lpos[nd];
		if (matuse[cpos]) msa->aseq[i][matmap[cpos]] = '-';
	      }
	    if (cm->stid[v] == MATP_D || cm->stid[v] == MATR_D) 
	      {
		cpos = emap->rpos[nd];
		if (matuse[cpos]) msa->aseq[i][matmap[cpos]] = '-';
	      }
	    break;

	  } /* end of the switch statement */


	  prvnd = nd;
	} /* end traversal over trace i. */

      /* IL/EL Insertions are currently flush-left and IR insertions are currently flush-right.
       * This is pre-1.0 Infernal behavior. If(cm->align_opts & CM_ALIGN_FLUSHINSERTS) we leave them all alone,
       * otherwise we regularize (split) the internal inserts, we flush the 5' inserts right and the 3'
       * inserts left (note: pre 1.0 behavior does the opposite, flushes 5' left (assuming they're ROOT_ILs)
       * and flushes 3' right (assuming they're ROOT_IRs).
       *
       * We have to be careful about EL's. We don't want to group IL/IR's and EL's together and then split them
       * because we need to annotate IL/IR's as '.'s in the consensus structure and EL's as '~'. So we split
       * each group separately. There should only be either IL or IR's at any position (b/c assuming we've
       * detached the CM grammar ambiguity (which is default in cmbuild)). But we don't count on it here.
       */
      if(! (cm->align_opts & CM_ALIGN_FLUSHINSERTS)) /* default behavior, split insert in half */
	{
	  /* Deal with inserts before first consensus position, ILs, then ELs, then IRs
	   * IL's are flush left, we want flush right */
	  rightjustify(msa->abc, msa->aseq[i], maxil[0]);
	  /* EL's are flush left, we want flush right I think these are impossible, but just in case... */
	  rightjustify(msa->abc, msa->aseq[i]+maxil[0], maxel[0]);
	  /* IR's are flush right, we want flush right, do nothing */

	  /* split all internal insertions */
	  for (cpos = 1; cpos < emap->clen; cpos++) 
	    {
	      if(maxil[cpos] > 1) /* we're flush LEFT, want to split */
		{
		  apos = matmap[cpos]+1;
		  for (nins = 0; islower((int) (msa->aseq[i][apos])); apos++)
		    nins++;
		  nins /= 2;		/* split the insertion in half */
		  rightjustify(msa->abc, msa->aseq[i]+matmap[cpos]+1+nins, maxil[cpos]-nins);
		}
	      if(maxel[cpos] > 1) /* we're flush LEFT, want to split */
		{
		  apos = matmap[cpos]+1 + maxil[cpos];
		  for (nins = 0; islower((int) (msa->aseq[i][apos])); apos++)
		    nins++;
		  nins /= 2;		/* split the insertion in half */
		  rightjustify(msa->abc, msa->aseq[i]+matmap[cpos]+1+maxil[cpos]+nins, maxel[cpos]-nins);
		}
	      if(maxir[cpos] > 1) /* we're flush RIGHT, want to split */
		{
		  apos = matmap[cpos+1]-1;
		  for (nins = 0; islower((int) (msa->aseq[i][apos])); apos--)
		    nins++;
		  nins ++; nins /= 2;		/* split the insertion in half (++ makes it same behavior as IL/EL */
		  leftjustify(msa->abc, msa->aseq[i]+matmap[cpos]+1 + maxil[cpos] + maxel[cpos], maxir[cpos]-nins);
		}
	    }
	  /* Deal with inserts after final consensus position, IL's then EL's, then IR's
	   * IL's are flush left, we want flush left, do nothing 
	   * EL's are flush left, we want flush left, do nothing 
	   * IR's are flush right, we want flush left */
	  leftjustify(msa->abc, msa->aseq[i]+matmap[emap->clen]+1 + maxil[emap->clen] + maxel[emap->clen], maxir[emap->clen]);
	}
    } /* end loop over all parsetrees */


  /* Gee, wasn't that easy?
   * Add the rest of the ("optional") information to the MSA.
   */
  CreateCMConsensus(cm, abc, 3.0, 1.0, &con);

  /* "author" info */
  ESL_ALLOC(msa->au, sizeof(char) * (strlen(PACKAGE_VERSION)+10));
  sprintf(msa->au, "Infernal %s", PACKAGE_VERSION);

  for (i = 0; i < nseq; i++)
    {
      if((status = esl_strdup(sq[i]->name, -1, &(msa->sqname[i]))) != eslOK) goto ERROR;
      /* TODO: individual SS annotations
       */
      if (wgt == NULL) msa->wgt[i] = 1.0;
      else             msa->wgt[i] = wgt[i];
    }

  /* Construct the primary sequence consensus/reference coordinate system line,
   * msa->rf. Actually it's not strictly necessary b/c we'll only use this MSA
   * to infer the implicity master parsetrees, but we can do it easily, so we do.
   * 
   * Note, we don't construct the secondary structure consensus line, msa->ss_cons.
   * We don't need to because we're only going to use this alignment to pass to
   * Alignment2Parsetrees() to get the implicit master CM parsetrees. If we wanted
   * the SS_cons line, we'd need to pass the master CM into this function, and have
   * a map from master to minor consensus columns (though we could probably infer this
   * from masteradd)
   *
   */
  ESL_ALLOC(msa->rf,      (sizeof(char) * (alen+1)));
  for (cpos = 0; cpos <= emap->clen; cpos++) 
    {
      if (matuse[cpos]) 
	{ /* CMConsensus is off-by-one right now, 0..clen-1 relative to cpos's 1..clen */

	  /* bug i1, xref STL7 p.12. Before annotating something as a base pair,
	   * make sure the paired column is also present.
	   */
	  if (con->ct[cpos-1] != -1 && matuse[con->ct[cpos-1]+1] == 0) {
	    msa->rf[matmap[cpos]]      = con->cseq[cpos-1];
	  } else {
	    msa->rf[matmap[cpos]]      = con->cseq[cpos-1];
	  }
	}
      if (maxil[cpos] > 0) 
	for (apos = ilmap[cpos]; apos < ilmap[cpos] + maxil[cpos]; apos++)
	  {
	    msa->rf[apos] = '.';
	  }
      if (maxel[cpos] > 0)
	for (apos = elmap[cpos]; apos < elmap[cpos] + maxel[cpos]; apos++)
	  {
	    msa->rf[apos] = '~';
	  }
      if (masteradd[cpos] > 0) { /* add master consensus columns not in this cm's consensus as '+' */
	for (apos = mastermap[cpos]; apos < mastermap[cpos] + masteradd[cpos]; apos++) 
	  { 
	    msa->rf[apos] = '+';
	  }
      }	
      if (maxir[cpos] > 0)	/* remember to write backwards */
	for (apos = irmap[cpos]; apos > irmap[cpos] - maxir[cpos]; apos--)
	  {
	    msa->rf[apos] = '.';
	  }
    }
  msa->rf[alen] = '\0';

  /* If we only want the match columns, shorten the alignment
   * by getting rid of the inserts. (Alternatively we could probably
   * simplify the building of the alignment, but all that pretty code
   * above already existed, so we do this post-msa-building shortening).
   */
  if(do_matchonly)
    {
      int *useme;
      ESL_ALLOC(useme, sizeof(int) * (msa->alen));
      esl_vec_ISet(useme, msa->alen, FALSE);
      for(cpos = 0; cpos <= emap->clen; cpos++)
	if(matmap[cpos] != -1) useme[matmap[cpos]] = TRUE;
      if((status = esl_msa_ColumnSubset(msa, errbuf, useme)) != eslOK) cm_Fail(errbuf);
      free(useme);
    }

  FreeCMConsensus(con);
  FreeEmitMap(emap);
  free(matuse);
  free(iluse);
  free(eluse);
  free(iruse);
  free(maxil);
  free(maxel);
  free(maxir);
  free(matmap);
  free(ilmap);
  free(elmap);
  free(irmap);
  free(mastermap);
  *ret_msa = msa;
  return eslOK;

 ERROR:
  if(con   != NULL)  FreeCMConsensus(con);
  if(emap  != NULL)  FreeEmitMap(emap);
  if(matuse!= NULL)  free(matuse);
  if(iluse != NULL)  free(iluse);
  if(eluse != NULL)  free(eluse);
  if(iruse != NULL)  free(iruse);
  if(maxil != NULL)  free(maxil);
  if(maxel != NULL)  free(maxel);
  if(maxir != NULL)  free(maxir);
  if(matmap!= NULL)  free(matmap);
  if(ilmap != NULL)  free(ilmap);
  if(elmap != NULL)  free(elmap);
  if(irmap != NULL)  free(irmap);
  if(msa   != NULL)  esl_msa_Destroy(msa);
  return status;
}


/* Function: validate_meta()
 * 
 * Incept:   EPN, Wed Jul 30 15:52:46 2008
 * 
 * Purpose:  Validate that the CMs read from the CM file in <cmlist> were
 *           indeed built in order from the alignments in cfg->mali_msa.
 * 
 *           This is done rather inelegantly by actually building new CMs
 *           from each alignment and verifying that the guidetree is identical
 *           to the guidetree for each CM. This is not 100% robust either,
 *           eventually we should use some type of checksum, or just take
 *           an multiple multiple alignment (the cmbuild input used to build the
 *           CMs) as input to cmalign, and just build them within cmalign.
 * 
 * Args:     go      - getopts 
 *           cfg     - the model
 *           errbuf  - for error messages 
 *           cmlist  - [0..m..cfg->ncm-1] the CMs 
 *           ret_toadd2min   - RETURN: (see comments in code)
 *           ret_maj_train_a2c_map - RETURN: (see comments in code)
 */
int
validate_meta(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t **cmlist, int ***ret_toadd2min, int **ret_maj_train_a2c_map)
{
  int          status;
  CM_t *tmp_cm = NULL;
  CM_t *tmp2_cm = NULL;
  int m, j, i;
  int **maj2min_cmap; /* [0..m..cfg->ncm][0..c..cmlist[0]->clen] = x, consensus column c of cmlist[0] (the major CM) maps to consensus column x of CM cmlist[m], (minor CM) */
  int **min2maj_cmap; /* [0..m..cfg->ncm][0..c..cmlist[m]->clen] = x, consensus column c of cmlist[m] (a minor CM)   maps to consensus column x of CM cmlist[0], (major CM) */
  int **toadd2min;    /* [0..m..cfg->ncm][0..c..cmlist[m]->clen] = x, when morphing minor parsetrees to a major alignment, add x major consensus columns after minor consensus column c of CM m, 
		       *                                              these x major consensus columns DO NOT MAP to a consensus column in minor CM m */
  int *maj_train_a2c_map;   /* [0..a..cfg->mali_msa[0]->alen] = x, alignment column a from the first training alignment (cfg->mali_msa[0]) maps to major consensus column x */
  int cpos, apos;
  int major_is_consensus;
  int *cposA; /* temporary only, [0..m..cfg->ncm-1] = x, x is current consensus column for CM m */

  /* Validation 1: Build CMs from each alignment, the guidetree should match the corresponding CM guidetree we read from the file */
  ESL_ALLOC(cfg->mali_mtr, sizeof(Parsetree_t *) * cfg->ncm);
  for(m = 0; m < cfg->ncm; m++) { 
    if(cfg->mali_msa[m]->rf == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "with -M, all alignments from %s must have #=GC RF notation, but alignment %d does not.", esl_opt_GetString(go, "-M"), m+1);
    if((status = HandModelmaker(cfg->mali_msa[m], errbuf, TRUE, 0.5, &tmp_cm, &(cfg->mali_mtr[m]))) != eslOK) return status;
    /*                              !use RF! */
    if(!(CompareCMGuideTrees(cmlist[m], tmp_cm))) { 
      tmp2_cm = CMRebalance(tmp_cm);
      FreeCM(tmp_cm);
      tmp_cm = NULL;
      if(!(CompareCMGuideTrees(cmlist[m], tmp2_cm))) cm_Fail("with -M, CM %d could not have been built by aligment %d. Did you remember to use --rf to cmbuild?", m+1, m+1);
      FreeCM(tmp2_cm);
      tmp2_cm = NULL;
    }
    if(tmp_cm != NULL)  FreeCM(tmp_cm);
  }
  
  /* Validation 2: the major (first) alignment must contain all the sequences each of the other alignments, 
   *               no additional sequences can exist in the 2nd->Nth (minor) alignments that are not in the major aln,
   *               and all sequences in the major alignment must exist in exactly 1 (not 0) minor alignments. 
   *               They must be in order also. 
   *               Further, all the alignments must be the same length, and the aligned sequences in each of the 
   *               minor alignments must exactly match (be the exact same aligned sequence) as the corresponding 
   *               sequence in the major alignment. 
   */
  m = 1;
  j = 0;
  for(i = 0; i < cfg->mali_msa[0]->nseq; i++) { 
    if(j == cfg->mali_msa[m]->nseq) { /* minor alignment m validated, move to next alignment */
      m++; 
      if(m >= cfg->mali_n) ESL_FAIL(eslEINCOMPAT, errbuf, "with -M, ran out of minor alignments, didn't account for all the sequences in the major alignment.");
      j = 0; 
      if(cfg->mali_msa[m]->alen != cfg->mali_msa[0]->alen) ESL_FAIL(eslEINCOMPAT, errbuf, "with -M, all alignments must have the same number of columns, but alignment 1 has %" PRId64 " columns, and alignment %d has %" PRId64 "  columns\n", cfg->mali_msa[0]->alen, m, cfg->mali_msa[m]->alen);
    }
    /* compare sequence name */
    if(strcmp(cfg->mali_msa[0]->sqname[i], cfg->mali_msa[m]->sqname[j]) != 0)
      ESL_FAIL(eslEINCOMPAT, errbuf, "with -M, sequence names for alignment 1 seq %d (%s) and alignment %d seq %d (%s) didn't match.", i+1, cfg->mali_msa[0]->sqname[i], m+1, j+1, cfg->mali_msa[m]->sqname[j]);
    /* compare actual sequence, digitized ax */
    if (memcmp(cfg->mali_msa[0]->ax[i], cfg->mali_msa[m]->ax[j], sizeof(ESL_DSQ) * (cfg->mali_msa[0]->alen)) != 0) 
      ESL_FAIL(eslEINCOMPAT, errbuf, "with -M, aligned sequence mismatch, alignment 1 seq %d (%s) and alignment %d seq %d (%s) are not identical.", i+1, cfg->mali_msa[0]->sqname[i], m+1, j+1, cfg->mali_msa[m]->sqname[j]);
    j++; /* move to next seq in alignment m */
  }

  /* Validation 3: The consensus columns of the major CM must be a superset of the consensus columns of all the minor CMs.
   *               As we check this we create maps from the major to each of the minor guide trees and vice versa, 
   *               We already checked all the training alignments are the same length so that we can easily map the coordinates.
   *               In fact these maps are currently unnecessary, but I left the code in b/c they may prove necessary in 
   *               the future. All we need currently is (1) toadd2min array which specifies the number of major consensus
   *               columns that do not exist in the minor consensus that should appear after each minor consensus column.
   *               and (2) maj_train_a2c_map which maps the columns of the major CM training alignment onto major consensus columns.
   */

  ESL_ALLOC(cposA, sizeof(int) * cfg->ncm);
  esl_vec_ISet(cposA, cfg->ncm, 0);

  ESL_ALLOC(maj2min_cmap, sizeof(int *) * cfg->ncm);
  ESL_ALLOC(min2maj_cmap, sizeof(int *) * cfg->ncm);
  ESL_ALLOC(toadd2min,    sizeof(int *) * cfg->ncm);
  for(m = 0; m < cfg->ncm; m++) { 
    ESL_ALLOC(maj2min_cmap[m], sizeof(int) * (cmlist[0]->clen+1));
    ESL_ALLOC(min2maj_cmap[m], sizeof(int) * (cmlist[m]->clen+1));
    ESL_ALLOC(toadd2min[m],    sizeof(int) * (cmlist[m]->clen+1));
    esl_vec_ISet(maj2min_cmap[m], (cmlist[0]->clen+1), -1);
    esl_vec_ISet(min2maj_cmap[m], (cmlist[m]->clen+1), -1);
    esl_vec_ISet(toadd2min[m],    (cmlist[m]->clen+1), 0);
  }    
  /* maj2min_cmap[0] and min2maj_cmap[0] is a map of the major CM to itself, this is trival, useless, and could just as easily be set to NULL */
  for(cpos = 1; cpos <= cmlist[0]->clen; cpos++) maj2min_cmap[0][cpos] = min2maj_cmap[0][cpos];

  ESL_ALLOC(maj_train_a2c_map, sizeof(int) * (cfg->mali_msa[0]->alen+1));
  esl_vec_ISet(maj_train_a2c_map, (cfg->mali_msa[0]->alen+1), -1);
  
  /* map major 2 minor consensus columns */
  for(apos = 1; apos <= cfg->mali_msa[0]->alen; apos++) { /* remember all the alignments are the same length, we checked */
    major_is_consensus = FALSE;
    for(m = 0; m < cfg->ncm; m++) { 
      if(! esl_abc_CIsGap(cfg->mali_msa[m]->abc, cfg->mali_msa[m]->rf[(apos-1)])) { 
	cposA[m]++;
	if(m == 0) { 
	  major_is_consensus = TRUE;
	  maj_train_a2c_map[apos] = cposA[0];
	}
	if(!major_is_consensus) ESL_FAIL(eslEINCOMPAT, errbuf, "with -M, alignment %d, position %d is a consensus column %d, but position %d is not consensus in the major (first) alignment from %s.", m, apos, cposA[m], apos, esl_opt_GetString(go, "-M"));
	maj2min_cmap[m][cposA[0]] = cposA[m];
	min2maj_cmap[m][cposA[m]] = cposA[0];
      }
      else if(major_is_consensus) { /* apos is major consensus column cposA[0], apos is not a minor consensus column in CM m */
	/* could assert maj2min_cmap[0][cposA[0]] == -1 */
	toadd2min[m][cposA[m]]++;
      }
    }
  }

  *ret_toadd2min   = toadd2min;
  *ret_maj_train_a2c_map = maj_train_a2c_map;

  /* free the maj2min and min2maj arrays, we didn't use them for anything,
   * they're left in for possible future use */
  for(m = 0; m < cfg->ncm; m++) { 
    free(maj2min_cmap[m]);
    free(min2maj_cmap[m]);
  }
  free(maj2min_cmap);
  free(min2maj_cmap);
  free(cposA);

  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "validate_meta(), memory allocation error.");
  return status; /* never reached */
}


/* Function: major_alignment2minor_parsetrees()
 * 
 * Incept:   EPN, Wed Jul 30 17:28:54 2008
 * 
 * Purpose:  Convert the major alignment into it's implicit minor parsetrees,
 *           and determine the highest scoring minor CM for each sequence.
 * 
 * Args:     go      - getopts 
 *           cfg     - the model
 *           errbuf  - for error messages 
 *           cmlist  - [0..m..cfg->ncm-1] the CMs 0 is major, 1..ncm-1 are minor
 *           maj_msa - major alignment
 *           maj_train_a2c_map - [0..a..cfg->mali_msa[0]->alen] = x, alignment column a from the first training alignment (cfg->mali_msa[0]) maps to major consensus column x 
 *           ret_wcm - RETURN: [0..i..maj_msa->nseq-1] highest scoring minor CM for seq i
 *
 * Returns:  eslOK on success; ret_wcm filled; 
 *
 * Throws:   eslEINCOMPAT on contract violation; eslEMEM on memory error.
 */
int
major_alignment2minor_parsetrees(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t **cmlist, ESL_MSA *maj_target_msa, int *maj_train_a2c_map, int **ret_wcm)
{
  int          status;
  int           m,i,x;     /* counters */
  int *maj_target_c2a_map; /* [0..c..cmlist[0]->clen-1] = x, major consensus column c maps to target (major alignment of target seqs) alignment column x */
  int  maj_train_cpos;     /* consensus position in major training alignment */
  int target_apos;         /* alignment position in (implied) minor target alignment */
  int tmp_clen;            /* temporary consensus length */
  ESL_SQ **tmp_sq;         /* temporary sequences */
  Parsetree_t **tmp_tr;    /* temporary parsetrees */
  float **totscAA;         /* [0..m..cfg->ncm-1][0..i..nseq-1] = x, CM m, sequence i total parsetree bit score = x */
  float **pscAA;           /* [0..m..cfg->ncm-1][0..i..nseq-1] = x, CM m, sequence i primary sequence emission score = x */
  float **sscAA;           /* [0..m..cfg->ncm-1][0..i..nseq-1] = x, CM m, sequence i structural score = x */
  float **tscAA;           /* [0..m..cfg->ncm-1][0..i..nseq-1] = x, CM m, sequence i summed transition bit score = x */
  float  *n3scA;           /* [0..i..nseq-1] = x, sequece i null3 penalty = x (same for all CMs) */
  float  *wsc;             /* [0..i..nseq-1] winning score, max score (if --sub: max psc, else: max totsc) for seq i over all minor parsetrees */
  int    *wcm;             /* [0..i..nseq-1] minor CM with highest score for seq i */
  int     do_sub;          /* TRUE if --sub enabled, FALSE if not */
  float   sc;              /* a temporary score */
  int     namewidth;       /* max strlen of all names in the alignment */
  char   *namedashes;      /* string of exactly namewidth '-'s */
  int     ni;              /* ctr */
  do_sub = (esl_opt_GetBoolean(go, "--sub")) ? TRUE : FALSE;

  /* map the major consensus positions to alignment positions in maj_target_msa (this is different from maj_train_a2c_map which maps alignment positions of training
   * alignment to consensus columns of major CM */
  if((status = map_cpos_to_apos(maj_target_msa, &maj_target_c2a_map, &tmp_clen))  != eslOK) ESL_FAIL(status, errbuf, "major_alignment2minor_parsetrees(), problem mapping major consensus positions to alignment positions.");
  if(tmp_clen != cmlist[0]->clen) ESL_FAIL(status, errbuf, "major_alignment2minor_parsetrees(), major target alignment clen %d != major cm clen: %d\n", tmp_clen, cmlist[0]->clen);

  /* convert the train guidetree coords for each CM to major CM coords for the alignment in maj_target_msa */
  for(m = 0; m < cfg->ncm; m++) { 
    for(x = 0; x < cfg->mali_mtr[m]->n; x++) { 
      /*printf("m: %3d x: %3d emitl: %3d emitr: %3d", m, x, cfg->mali_mtr[m]->emitl[x], cfg->mali_mtr[m]->emitr[x]);*/
      maj_train_cpos             = maj_train_a2c_map[cfg->mali_mtr[m]->emitl[x]];
      target_apos                = maj_target_c2a_map[maj_train_cpos];
      cfg->mali_mtr[m]->emitl[x] = target_apos;

      maj_train_cpos             = maj_train_a2c_map[cfg->mali_mtr[m]->emitr[x]];
      target_apos                = maj_target_c2a_map[maj_train_cpos];
      cfg->mali_mtr[m]->emitr[x] = target_apos;
      /*printf("  newl: %3d newr: %3d\n", cfg->mali_mtr[m]->emitl[x], cfg->mali_mtr[m]->emitr[x]);*/
    }
  }

  /* for each minor CM, transmogrify the major CM alignment into their implicit minor parsetrees */
  ESL_ALLOC(wsc, sizeof(float) * maj_target_msa->nseq);
  ESL_ALLOC(wcm, sizeof(int) * maj_target_msa->nseq);
  esl_vec_FSet(wsc, maj_target_msa->nseq, IMPOSSIBLE);
  esl_vec_ISet(wcm, maj_target_msa->nseq, -1);

  ESL_ALLOC(totscAA, sizeof(float *) * cfg->ncm);
  ESL_ALLOC(pscAA,   sizeof(float *) * cfg->ncm);
  ESL_ALLOC(sscAA,   sizeof(float *) * cfg->ncm);
  ESL_ALLOC(tscAA,   sizeof(float *) * cfg->ncm);
  for(m = 0; m < cfg->ncm; m++) { 
    ESL_ALLOC(totscAA[m], sizeof(float) * maj_target_msa->nseq);
    ESL_ALLOC(pscAA[m],   sizeof(float) * maj_target_msa->nseq);
    ESL_ALLOC(sscAA[m],   sizeof(float) * maj_target_msa->nseq);
    ESL_ALLOC(tscAA[m],   sizeof(float) * maj_target_msa->nseq);
  }
  ESL_ALLOC(n3scA,  sizeof(float) * maj_target_msa->nseq);

  if((status = esl_msa_Digitize(cmlist[0]->abc, maj_target_msa, NULL)) != eslOK) ESL_FAIL(eslEINCOMPAT, errbuf, "Failure digitizing the major target CM alignment.");
  for(m = 0; m < cfg->ncm; m++) { 
    if((status = Alignment2Parsetrees(maj_target_msa, cmlist[m], cfg->mali_mtr[m], errbuf, &tmp_sq, &tmp_tr)) != eslOK) return status;
    for(i = 0; i < maj_target_msa->nseq; i++) { 
      if((status = ParsetreeScore(cmlist[m], NULL, errbuf, tmp_tr[i], tmp_sq[i]->dsq, FALSE, &(totscAA[m][i]), &(sscAA[m][i]), &(pscAA[m][i]), NULL, NULL)) != eslOK) return status;
      //if(m == 0) if((status = ParsetreeScoreCorrectionNull3(cmlist[m], errbuf, tmp_tr[i], tmp_sq[i]->dsq, 0, &(n3scA[i]))) != eslOK) return status;
      tscAA[m][i] = totscAA[m][i] + n3scA[i] - pscAA[m][i] - sscAA[m][i];
      sc = do_sub ? pscAA[m][i] : totscAA[m][i];
      if(sc > wsc[i]) { wsc[i] = sc; wcm[i] = m; } 
      esl_sq_Destroy(tmp_sq[i]);
      FreeParsetree(tmp_tr[i]);
    }
    free(tmp_sq);
    free(tmp_tr);
  }

  /* output scores in tabular format */
  namewidth = 8; /* length of 'seq name' */
  /* determine the longest name in msa */
  for(ni = 0; ni < maj_target_msa->nseq; ni++) namewidth = ESL_MAX(namewidth, strlen(maj_target_msa->sqname[ni]));
  ESL_ALLOC(namedashes, sizeof(char) * namewidth+1);
  namedashes[namewidth] = '\0';
  for(ni = 0; ni < namewidth; ni++) namedashes[ni] = '-';

  printf("\n\n\n");
  printf("# %7s  %-*s", "seq idx",  namewidth, "seq name");
  for(m = 0; m < cfg->ncm; m++) { printf("  %5s %-3d", "CM", m); } 
  printf("\n");
  printf("# %7s  %-*s", "-------",  namewidth, namedashes);
  for(m = 0; m < cfg->ncm; m++) { printf("  %9s", "---------"); } 
  printf("\n");

  for(i = 0; i < maj_target_msa->nseq; i++) { 
    printf("  %7d  %-*s  ", i+1, namewidth, maj_target_msa->sqname[i]);
    for(m = 0; m < (cfg->ncm-1); m++) { 
      printf("%1s", (wcm[i] == m) ? "*" : " ");
      if(do_sub) printf("%8.3f  ", pscAA[m][i]);
      else       printf("%8.3f  ", totscAA[m][i]);
    }
    printf("%1s", (wcm[i] == m) ? "*" : " ");
    if(do_sub) printf("%8.3f\n", pscAA[m][i]);
    else       printf("%8.3f\n", totscAA[m][i]);
  }
  for(m = 0; m < cfg->ncm; m++) {
    free(totscAA[m]);
    free(pscAA[m]);
    free(sscAA[m]);
    free(tscAA[m]);
  }
  free(totscAA);
  free(pscAA);
  free(sscAA);
  free(tscAA);
  free(n3scA);
  free(wsc);
  free(namedashes);
  free(maj_target_c2a_map);

  *ret_wcm = wcm;
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "major_alignment2minor_parsetrees(), memory allocation error.");
  return status; /* never reached */
}  

/* Function: majorfied_alignment2major_parsetrees()
 * 
 * Incept:   EPN, Thu Jul 31 09:23:30 2008
 * 
 * Purpose:  Convert an inferred major alignment (a majorfied alignment determined
 *           from minor parsetrees) into it's implied MSA.
 * 
 * Args:     go      - getopts 
 *           cfg     - the model
 *           errbuf  - for error messages 
 *           maj_cm  - the major CM
 *           majed_min_msa - majorfied minor alignment
 *           maj_emap - emit map for maj_cm
 *           ret_majed_tr - RETURN: [0..i..majed_min_msa->nseq-1] major parsetree for each sequence
 *
 * Returns:  eslOK on success; ret_majed_tr filled; 
 *
 * Throws:   eslEINCOMPAT on contract violation; 
 */
int
majorfied_alignment2major_parsetrees(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *maj_cm, ESL_MSA *majed_min_msa, CMEmitMap_t *maj_emap, Parsetree_t ***ret_majed_tr)
{
  int          status;
  int           m = 0;     /* counter */
  int           x = 0;     /* counter */
  int *majed_min_c2a_map;  /* [0..c..maj_clen] = a, consensus position c of majed_min_msa is alignment position a */
  int  majed_min_clen;     /* number of consensus columns parsed in majed_min_msa, this better = maj_cm->clen */
  Parsetree_t **majed_tr;  /* major parsetrees */
    
  if((status = esl_msa_Digitize(maj_cm->abc, majed_min_msa, NULL)) != eslOK) ESL_FAIL(eslEINCOMPAT, errbuf, "Failure digitizing the minor CM alignment for CM %d of winning seqs.", m);
  if((status = map_cpos_to_apos(majed_min_msa, &majed_min_c2a_map, &majed_min_clen))  != eslOK) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "Problem mapping consensus positions to alignment positions for majorfied minor alignment %d.", m);
  if(majed_min_clen != maj_cm->clen) ESL_FAIL(eslEINCOMPAT, errbuf, "Majorfied minor alignment has clen != major clen: %d\n", maj_cm->clen);

  /* convert the major CM guidetree coords from the training alignment coords, to the majorfied minor CM alignment coords of current seqs */
  /* handle the ROOT (x == 0) node special, it must have emitl == 1, emitr == msa->alen */
  cfg->mali_mtr[0]->emitl[0] = 1;
  cfg->mali_mtr[0]->emitr[0] = majed_min_msa->alen;
  for(x = 1; x < cfg->mali_mtr[0]->n; x++) { 
    //printf("x: %4d ol: %4d or: %4d ", x, cfg->mali_mtr[0]->emitl[x], cfg->mali_mtr[0]->emitr[x]);
    /* careful: we have to correct for an off-by-one b/t how non-MATL non-MATP nodes emitmap's lpos (in CreateEmitMap()) and guidetree emitl's (in HandModelMaker()) are calculated */
    if(maj_cm->ndtype[x] == MATL_nd || maj_cm->ndtype[x] == MATP_nd) 
      cfg->mali_mtr[0]->emitl[x] = majed_min_c2a_map[maj_emap->lpos[x]];
    else if (maj_cm->ndtype[x] == BEGR_nd) 
      cfg->mali_mtr[0]->emitl[x] = majed_min_c2a_map[maj_emap->lpos[x]]+1;
    else 
      cfg->mali_mtr[0]->emitl[x] = majed_min_c2a_map[maj_emap->lpos[x]+1]; 

    /* careful: we have to correct for an off-by-one b/t how non-MATR non-MATP nodes emitmap's rpos (in CreateEmitMap()) and guidetree emitr's (in HandModelMaker()) are calculated */
    if(maj_cm->ndtype[x] == MATR_nd || maj_cm->ndtype[x] == MATP_nd) 
      cfg->mali_mtr[0]->emitr[x] = majed_min_c2a_map[maj_emap->rpos[x]];
    else 
      cfg->mali_mtr[0]->emitr[x] = majed_min_c2a_map[maj_emap->rpos[x]-1];
    //printf(" nl: %4d nr: %4d\n", cfg->mali_mtr[0]->emitl[x], cfg->mali_mtr[0]->emitr[x]);
  }
  free(majed_min_c2a_map);

  /* esl_msa_Write(stdout, majed_min_msa, eslMSAFILE_STOCKHOLM); */
  if((status = Alignment2Parsetrees(majed_min_msa, maj_cm, cfg->mali_mtr[0], errbuf, NULL, &majed_tr)) != eslOK) return status;
  
  *ret_majed_tr = majed_tr;
  return eslOK;
}


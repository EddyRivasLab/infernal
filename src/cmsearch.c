/* cmsearch.c
 * SRE, Fri May  3 13:58:18 2002
 * SVN $Id$
 * 
 * Search sequences with a CM.
 * 
 *****************************************************************
 * @LICENSE@
 ***************************************************************** 
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <time.h>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "easel.h"              
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_mpi.h"
#include "esl_msa.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_sqio.h"
#include "esl_sqio_ascii.h"
#include "esl_ssi.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#include "funcs.h"		/* external functions                   */
#include "structs.h"		/* data structures, macros, #define's   */

#define STRATOPTS1  "--cyk,--inside,--viterbi,--forward"               /* incompatible with --cyk, --inside (besides themselves) */
#define STRATOPTS2  "--cyk,--inside,--viterbi,--forward,--fil-hmm"     /* incompatible with --viterbi,--forward  (besides themselves) */
#define ALPHOPTS    "--rna,--dna"                                      /* exclusive choice for output alphabet */

#define CUTOPTS1    "-E,-T,--ga,--tc,--nc"                             /* incompatible with -E, -T (besides themselves) */
#define CUTOPTS2    "-E,-T,--ga,--tc,--nc,--viterbi,--forward"         /* incompatible with --ga, --tc, --nc (besides themselves) */
#define HMMONLYOPTS "--viterbi,--forward"                              /* with these options, there are no filters, use only HMM */
								     
static ESL_OPTIONS options[] = {
  /* name             type            default     env  range               toggles       reqs       incomp  help                                      docgroup*/
  /* basic options */
  { "-h",             eslARG_NONE,    FALSE,     NULL, NULL,                  NULL,      NULL,        NULL, "show brief help on version and usage",   1 },
  { "-o",             eslARG_OUTFILE, NULL,      NULL, NULL,                  NULL,      NULL,        NULL, "direct output to file <f>, not stdout", 1 },
  { "-g",             eslARG_NONE,    FALSE,     NULL, NULL,                  NULL,      NULL,        NULL, "configure CM/HMM for glocal alignment [default: local]", 1 },
  { "-p",             eslARG_NONE,    FALSE,     NULL, NULL,                  NULL,"--aln-hbanded", "--noalign", "append posterior probabilities to hit alignments", 1 },
  { "-x",             eslARG_NONE,    FALSE,     NULL, NULL,                  NULL,      NULL,"--noalign,-v","annotate non-compensatory bps in output alignments with 'x'", 1 },
  { "-v",             eslARG_NONE,    FALSE,     NULL, NULL,                  NULL,      NULL,"--noalign,-x","annotate negative scoring non-canonical bps with 'v'", 1 },
  { "-Z",             eslARG_REAL,    FALSE,     NULL, NULL,                  NULL,      NULL,        NULL, "set Z (database size in *Mb*) to <x> for E-value calculations", 1},
  { "--toponly",      eslARG_NONE,    FALSE,     NULL, NULL,                  NULL,      NULL,        NULL, "only search the top strand", 1 },
  { "--bottomonly",   eslARG_NONE,    FALSE,     NULL, NULL,                  NULL,      NULL,        NULL, "only search the bottom strand", 1 },
  { "--forecast",     eslARG_INT,     NULL,      NULL, NULL,                  NULL,      NULL,        NULL, "don't do search, forecast running time with <n> processors", 1 },
  { "--informat",     eslARG_STRING,  NULL,      NULL, NULL,                  NULL,      NULL,        NULL, "specify the input file is in format <x>, not FASTA", 1 },
  { "--mxsize",       eslARG_REAL,    "2048.0",  NULL, "x>0.",                NULL,      NULL,        NULL, "set maximum allowable HMM banded DP matrix size to <x> Mb", 1 },
  { "--devhelp",      eslARG_NONE,    NULL,      NULL, NULL,                  NULL,      NULL,        NULL, "show list of undocumented developer options", 1 },
#ifdef HAVE_MPI
  { "--mpi",          eslARG_NONE,    FALSE,     NULL, NULL,                  NULL,      NULL,        NULL, "run as an MPI parallel program", 1 },  
#endif
  /* options for algorithm for final round of search */
  { "--inside",       eslARG_NONE,    "default", NULL, NULL,                  NULL,      NULL,    STRATOPTS1, "use scanning CM Inside algorithm", 2 },
  { "--cyk",          eslARG_NONE,    FALSE,     NULL, NULL,            "--inside",      NULL,    STRATOPTS1, "use scanning CM CYK algorithm", 2 },
  { "--forward",      eslARG_NONE,    FALSE,     NULL, NULL, "--fil-hmm,--fil-qdb",      NULL,    STRATOPTS2, "use scanning HMM Forward algorithm", 2 },
  { "--viterbi",      eslARG_NONE,    FALSE,     NULL, NULL, "--fil-hmm,--fil-qdb",      NULL,    STRATOPTS2, "use scanning HMM Viterbi algorithm", 2 },
  /* CM cutoff options */
  /* IMPORTANT: Default values for -E and -T must remain non-NULL, the option processing logic below depends on it */
  { "-E",             eslARG_REAL,    "1.0",     NULL, "x>0.",                NULL,      NULL,    CUTOPTS1, "use cutoff E-value of <x> for final round of search", 3 },
  { "-T",             eslARG_REAL,    "0.0",     NULL, NULL,                  NULL,      NULL,    CUTOPTS1, "use cutoff bit score of <x> for final round of search", 3 },
  { "--nc",           eslARG_NONE,    FALSE,     NULL, NULL,                  NULL,      NULL,    CUTOPTS2, "use CM Rfam NC noise cutoff as cutoff bit score", 3 },
  { "--ga",           eslARG_NONE,    FALSE,     NULL, NULL,                  NULL,      NULL,    CUTOPTS2, "use CM Rfam GA gathering threshold as cutoff bit score", 3 },
  { "--tc",           eslARG_NONE,    FALSE,     NULL, NULL,                  NULL,      NULL,    CUTOPTS2, "use CM Rfam TC trusted cutoff as cutoff bit score", 3 },
  /* banded options (for final round of searching) */
  { "--no-qdb",       eslARG_NONE,    FALSE,     NULL, NULL,                  NULL,      NULL,  HMMONLYOPTS, "do not use QDBs in final round of searching", 4 },
  { "--beta",         eslARG_REAL,    "1e-15",   NULL, "0<x<1",               NULL,      NULL,  HMMONLYOPTS, "set tail loss prob for QDB calculation to <x>", 4 },
  { "--hbanded",      eslARG_NONE,    FALSE,     NULL, NULL,                  NULL,      NULL,  HMMONLYOPTS, "calculate and use HMM bands in final round of CM search", 4 },
  { "--tau",          eslARG_REAL,    "1e-7",    NULL, "0<x<1",               NULL,"--hbanded", HMMONLYOPTS, "set tail loss prob for --hbanded to <x>", 4 },
  /* filtering options, by default do HMM, then CYK filter */
  { "--fil-hmm",      eslARG_NONE,    "default", NULL, NULL,                  NULL,      NULL,"--fil-no-hmm", "filter with HMM Forward algorithm", 201 },
  { "--fil-no-hmm",   eslARG_NONE,    FALSE,     NULL, NULL,           "--fil-hmm",      NULL,          NULL, "do not filter with HMM Forward algorithm", 5 },
  { "--fil-qdb",      eslARG_NONE,    "default", NULL, NULL,                  NULL,      NULL,"--fil-no-qdb", "filter with CM QDB (banded) CYK algorithm", 201 },
  { "--fil-no-qdb",   eslARG_NONE,    FALSE,     NULL, NULL,           "--fil-qdb",      NULL,          NULL, "do not filter with CM banded CYK", 5 },
  { "--fil-beta",     eslARG_REAL,    "1e-10",   NULL, "x>0",                 NULL,      NULL,"--fil-no-qdb", "set tail loss prob for QDB filter to <x>", 5 },
  /* filter cutoff options */
  { "--fil-T-qdb",    eslARG_REAL,    "0.0",     NULL, NULL,                  NULL,      NULL, "--fil-E-qdb,--fil-no-qdb",           "set QDB CM filter cutoff bit score as <x>", 6 },
  { "--fil-T-hmm",    eslARG_REAL,    "3.0",     NULL, NULL,                  NULL,      NULL, "--fil-E-hmm,--fil-S-hmm,--fil-no-hmm","set HMM filter cutoff bit score as <x>", 6 },
  { "--fil-E-qdb",    eslARG_REAL,    NULL,      NULL, "x>0.999",             NULL,      NULL, "--fil-T-qdb,--fil-no-qdb",            "set QDB CM filter cutoff E-value as <x>", 6 },
  { "--fil-E-hmm",    eslARG_REAL,    NULL,      NULL, "x>0.999",             NULL,      NULL, "--fil-T-hmm,--fil-S-hmm,--fil-no-hmm","set HMM filter cutoff E-value as <x>", 6 }, 
  { "--fil-S-hmm",    eslARG_REAL,    NULL,      NULL, "0<x<1.001",           NULL,      NULL, "--fil-E-hmm,--fil-T-hmm,--fil-no-hmm","set HMM filter predicted surv fract as <x>", 6 }, 
  { "--fil-Xmin-hmm", eslARG_REAL,    NULL,      NULL, "x>1.0999",            NULL,      NULL, "--fil-no-hmm",                        "set min HMM surv fract such that total time is <x> * an HMM", 106 },
  { "--fil-Smax-hmm", eslARG_REAL,    "0.5",     NULL, "0<x<1.001",           NULL,      NULL, "--fil-T-hmm,--fil-E-hmm,--fil-S-hmm,--fil-no-hmm", "set maximum HMM survival fraction as <x>", 6 },
  { "--fil-Smin-hmm", eslARG_REAL,    "0.02",    NULL, "0<x<1.001",           NULL,      NULL, "--fil-T-hmm,--fil-E-hmm,--fil-S-hmm,--fil-no-hmm", "set minimum HMM survival fraction as <x>", 6 },
  { "--fil-A-hmm",    eslARG_NONE,    FALSE,     NULL, NULL,                  NULL,      NULL, "--fil-no-hmm",                        "always filter w/HMM w/surv fract <= <x> from --fil-Smax-hmm", 6 },
  { "--fil-finE-hmm", eslARG_REAL,    NULL,      NULL, "x>0.",                NULL,      NULL, "--fil-finT-hmm",                      "pretend final E cutoff=<x> for HMM filter cutoff calc", 106 }, 
  { "--fil-finT-hmm", eslARG_REAL,    NULL,      NULL, NULL,                  NULL,      NULL, "--fil-finE-hmm",                      "pretend final bit sc cutoff=<x> for HMM filter cutoff calc", 106 }, 
  { "--fil-finE-qdb", eslARG_REAL,    NULL,      NULL, "x>0.",                NULL,      NULL," --fil-finT-qdb",                      "pretend final E cutoff=<x> for QDB filter cutoff calc", 106 }, 
  { "--fil-finT-qdb", eslARG_REAL,    NULL,      NULL, NULL,                  NULL,      NULL, "--fil-finE-qdb",                      "pretend final bit sc cutoff=<x> for QDB filter cutoff calc", 106 }, 
  /* W definition options (require --viterbi or --forward) */
  { "--hmm-W",        eslARG_INT,     NULL,      NULL, "n>1",                 NULL,      NULL,  "--hmm-cW", "set HMM window size as <n>", 7 },
  { "--hmm-cW",       eslARG_REAL,    NULL,      NULL, "x>0.01",              NULL,      NULL,   "--hmm-W", "set HMM window size as <x> * consensus length", 7 },
  /* alignment options */
  { "--noalign",      eslARG_NONE,    FALSE,     NULL, NULL,                  NULL,      NULL,       NULL,        "find start/stop/score only; don't do alignments", 8 },
  { "--aln-hbanded",  eslARG_NONE,    FALSE,     NULL, NULL,                  NULL,      NULL,"--noalign",        "use HMM bands to align hits", 8 },
  { "--aln-optacc",   eslARG_NONE,    FALSE,     NULL, NULL,                  NULL, "--aln-hbanded", "--noalign", "align hits with the optimal accuracy algorithm, not CYK", 8 },
  /* Using only a single CM from a multi-CM file */
  { "--cm-idx",       eslARG_INT,     NULL,      NULL, "n>0",                 NULL,      NULL, "--cm-name", "only search with CM number <n>    in the CM file",  11 },
  { "--cm-name",      eslARG_STRING,  NULL,      NULL, NULL,                 NULL,      NULL,  "--cm-idx", "only search with the CM named <s> in the CM file",  11 },
  /* Specifying a range of sequences to search, instead of the entire sequence file */
  { "--sseq",         eslARG_INT,     NULL,      NULL,"n>0",                 NULL,      NULL,       NULL, "first seq to search in <seqfile> is seq number <n>",  12 },
  { "--eseq",         eslARG_INT,     NULL,      NULL,"n>0",                 NULL,      NULL,       NULL, "final seq to search in <seqfile> is seq number <n>",  12 },
  /* verbose output files */
  { "--tabfile",      eslARG_OUTFILE, NULL,      NULL, NULL,                 NULL,      NULL,"--forecast", "save hits in tabular format to file <f>", 9 },
  { "--gcfile",       eslARG_OUTFILE, NULL,      NULL, NULL,                 NULL,      NULL,"--sseq,--eseq", "save GC content stats of target sequence file to <f>", 9 },
  /* Setting output alphabet */
  { "--rna",          eslARG_NONE,"default",     NULL, NULL,             ALPHOPTS,      NULL,        NULL, "output hit alignments as RNA sequence data", 10 },
  { "--dna",          eslARG_NONE,   FALSE,      NULL, NULL,             ALPHOPTS,      NULL,        NULL, "output hit alignments as DNA (not RNA) sequence data", 10 },
  /* All options below are developer options, only shown if --devhelp invoked */
  { "--lambda",       eslARG_REAL,   NULL,       NULL, NULL,                 NULL,      NULL,          NULL,         "overwrite lambdas in <cmfile> to <x> for E-value calculations", 101}, 
  { "--aln2bands",    eslARG_NONE,   FALSE,      NULL, NULL,                 NULL, "--hbanded",   HMMONLYOPTS,       "w/-hbanded, derive HMM bands w/o scanning Forward/Backward", 101 },
  { "--rtrans",       eslARG_NONE,   FALSE,      NULL, NULL,                 NULL,      NULL, "--viterbi,--forward", "replace CM transition scores from <cmfile> with RSEARCH scores", 101 },
  { "--sums",         eslARG_NONE,   FALSE,      NULL, NULL,                 NULL,"--hbanded",       NULL,           "use posterior sums during HMM band calculation (widens bands)", 101 },
  { "--null2",        eslARG_NONE,   FALSE,      NULL, NULL,                 NULL,"--no-null3", "--noalign",         "turn on the post hoc second null model", 101 },
  { "--no-null3",     eslARG_NONE,   FALSE,      NULL, NULL,                 NULL,      NULL,        NULL,           "turn OFF the NULL3 post hoc additional null model", 101 },
  { "-s",             eslARG_INT,    "181",      NULL, "n>=0",               NULL,      NULL,          NULL,         "set RNG seed to <n> (if 0: one-time arbitrary seed)", 101 },
  { "--stall",        eslARG_NONE,   FALSE,      NULL, NULL,                 NULL,      NULL,        NULL,           "arrest after start: for debugging MPI under gdb", 101 },  
  /* Developer options related to experimental local begin/end modes */
  { "--pebegin",      eslARG_NONE,   FALSE,      NULL, NULL,                 NULL,      NULL, "-g,--pbegin","set all local begins as equiprobable", 102 },
  { "--pfend",        eslARG_REAL,   NULL,       NULL, "0<x<1",              NULL,      NULL, "-g,--pend",  "set all local end probs to <x>", 102 },
  { "--pbegin",       eslARG_REAL,   "0.05",     NULL, "0<x<1",              NULL,      NULL,        "-g",  "set aggregate local begin prob to <x>", 102 },
  { "--pend",         eslARG_REAL,   "0.05",     NULL, "0<x<1",              NULL,      NULL,        "-g",  "set aggregate local end prob to <x>", 102 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

/* struct cfg_s : "Global" application configuration shared by all threads/processes
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 * NOTE: MPI not yet implemented.
 */
struct cfg_s {
  char         *cmfile;	        /* name of input CM file  */ 
  char         *sqfile;	        /* name of sequence file  */ 
  ESL_SQFILE   *sqfp;           /* open sequence input file stream */
  FILE         *ofp;            /* output file (default is stdout) */
  int           fmt;		/* format code for seqfile */
  ESL_ALPHABET *abc;		/* digital alphabet for input */
  long          dbsize;         /* database size in nucleotides for E values (doubled if doing rev comp) */
  int           namewidth;      /* maximum name length in the database, imperative for obsessively pretty formatting in tab file output */
  int           nseq;           /* number of seqs in the target db */
  int           ncm;            /* number CM we're at in file */
  int           do_rc;          /* should we search reverse complement? (for convenience */
  int           init_rci;       /* initial strand to search 0 for top, 1 for bottom (only 1 if --bottomonly enabled) */
  float         avg_hit_len;    /* average CM hit length, calc'ed using QDB calculation algorithm */
  FILE         *gcfp;           /* optional output file for --gcfile */
  FILE         *tfp;            /* optional output file for --tab */

  int           do_mpi;		/* TRUE if we're doing MPI parallelization */
  int           nproc;		/* how many MPI processes, total */
  int           my_rank;	/* who am I, in 0..nproc-1 */
  int           do_stall;	/* TRUE to stall the program until gdb attaches */

  /* Masters only (mainly i/o streams) */
  CMFILE       *cmfp;		/* open input CM file stream       */
  ESL_ALPHABET *abc_out; 	/* digital alphabet for writing */
};

static char usage[]  = "[-options] <cmfile> <sequence file>";
static char banner[] = "search a sequence database with an RNA CM";

static int  init_master_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
/* static int  init_shared_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf); */

static void  serial_master (const ESL_GETOPTS *go, struct cfg_s *cfg);
#ifdef HAVE_MPI
static int   mpi_master    (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int   mpi_worker    (const ESL_GETOPTS *go, struct cfg_s *cfg);
#endif
static int initialize_cm                        (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int read_next_search_seq                 (const ESL_ALPHABET *abc, ESL_SQFILE *seqfp, int do_revcomp, dbseq_t **ret_dbseq);
static int print_run_info                       (const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf);
extern int get_command                          (const ESL_GETOPTS *go, char *errbuf, char **ret_command);
static int set_searchinfo_for_calibrated_cm     (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int set_searchinfo_for_uncalibrated_cm   (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int print_searchinfo_for_calibrated_cm   (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float *cm_surv_fractA, int *cm_nhitsA, double in_asec, double in_total_psec, double *ret_total_psec);
static int print_searchinfo_for_uncalibrated_cm (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float *cm_surv_fractA, int *cm_nhitsA, double in_asec);
static int estimate_search_time_for_round       (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, int stype, int search_opts, ScanMatrix_t *smx, ESL_RANDOMNESS *r, double *ret_sec_per_res);
static int dump_gc_info                         (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int overwrite_lambdas                    (const ESL_GETOPTS *go, const struct cfg_s *cfg, CM_t *cm, char *errbuf);

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go = NULL;   /* command line processing                     */
  ESL_STOPWATCH   *w  = esl_stopwatch_Create();
  if(w == NULL) cm_Fail("Memory error, stopwatch not created.\n");
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
      puts("\nalgorithm for final round of search (after >= 0 filters): [default: --inside]");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\ncutoff options for final round of search (after >= 0 filters):");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80);
      puts("\noptions for banded DP in final round of search:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 
      puts("\nfiltering options:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80);
      puts("\nfilter cutoff options (survival fractions are predicted, not guaranteed):");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80);
      puts("\ndefining window size (W) for HMM only searches (require --forward or --viterbi):");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80);
      puts("\noptions for returning alignments of search hits:");
      esl_opt_DisplayHelp(stdout, go, 8, 2, 80);
      puts("\nusing a single CM from a multi-CM file:");
      esl_opt_DisplayHelp(stdout, go, 11, 2, 80);
      puts("\nspecifying a range of sequences to search, instead of the full target file:");
      esl_opt_DisplayHelp(stdout, go, 12, 2, 80);
      puts("\nverbose output files:");
      esl_opt_DisplayHelp(stdout, go, 9, 2, 80);
      puts("\noptions for selecting output alphabet:");
      esl_opt_DisplayHelp(stdout, go,10, 2, 80);
      puts("\nundocumented developer miscellaneous options:");
      esl_opt_DisplayHelp(stdout, go, 101, 2, 80);
      puts("\nundocumented developer options related to experimental local begin/end modes:");
      esl_opt_DisplayHelp(stdout, go, 102, 2, 80);
      puts("\nundocumented developer options related to filtering:");
      esl_opt_DisplayHelp(stdout, go, 106, 2, 80);
      exit(0);
    }
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      cm_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\nwhere general options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1=docgroup, 2 = indentation; 80=textwidth*/
      puts("\nalgorithm for final round of search (after >= 0 filters): [default: --inside]");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\ncutoff options for final round of search (after >= 0 filters):");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80);
      puts("\noptions for banded DP in final round of search:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 
      puts("\nfiltering options:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80);
      puts("\nfilter cutoff options (survival fractions are predicted, not guaranteed):");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80);
      puts("\ndefining window size (W) for HMM only searches (require --forward or --viterbi):");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80);
      puts("\noptions for returning alignments of search hits:");
      esl_opt_DisplayHelp(stdout, go, 8, 2, 80);
      puts("\nusing a single CM from a multi-CM file:");
      esl_opt_DisplayHelp(stdout, go, 11, 2, 80);
      puts("\nspecifying a range of sequences to search, instead of the full target file:");
      esl_opt_DisplayHelp(stdout, go, 12, 2, 80);
      puts("\nverbose output files:");
      esl_opt_DisplayHelp(stdout, go, 9, 2, 80);
      puts("\noptions for selecting output alphabet:");
      esl_opt_DisplayHelp(stdout, go,10, 2, 80);
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
  /* Check for incompatible option combinations I don't know how to disallow with esl_getopts */
  if ( esl_opt_IsOn(go, "--hmm-W") && (! ((esl_opt_GetBoolean(go, "--viterbi")) || (esl_opt_GetBoolean(go, "--forward"))))) { 
    printf("Error parsing options, --hmm-W <n> only makes sense in combination with --forward or --viterbi.\n");
    exit(1);
  }
  if ( esl_opt_IsOn(go, "--hmm-cW") && (! ((esl_opt_GetBoolean(go, "--viterbi")) || (esl_opt_GetBoolean(go, "--forward"))))) { 
    printf("Error parsing options, --hmm-cW <x> only makes sense in combination with --forward or --viterbi.\n");
    exit(1);
  }
  if ((esl_opt_IsUsed(go, "--fil-Smin-hmm")) && (esl_opt_IsUsed(go, "--fil-Smax-hmm"))) { 
    if (((esl_opt_GetReal(go, "--fil-Smin-hmm")) - (esl_opt_GetReal(go, "--fil-Smax-hmm"))) > eslSMALLX1) { 
      printf("Error parsing options, --fil-Smin-hmm <x> (%f) must be less than --fil-Smax-hmm <x> (%f).\n", esl_opt_GetReal(go, "--fil-Smin-hmm"), esl_opt_GetReal(go, "--fil-Smax-hmm"));
      exit(1);
    }
  }
  if ((esl_opt_IsOn(go, "--sseq")) && (esl_opt_IsOn(go, "--eseq"))) { 
    if((esl_opt_GetInteger(go, "--sseq")) > (esl_opt_GetInteger(go, "--eseq"))) { 
      printf("Error parsing options, both --sseq <x> and --eseq <y> are used, but <x> > <y>\n");
      exit(1);
    }
  }      
  
  /* Initialize what we can in the config structure (without knowing the input alphabet yet).
   */
  cfg.cmfile     = esl_opt_GetArg(go, 1); 
  cfg.sqfile     = esl_opt_GetArg(go, 2); 
  cfg.ofp        = NULL;
  cfg.gcfp       = NULL;
  cfg.tfp        = NULL;
  cfg.sqfp       = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */

  if   (esl_opt_IsOn(go, "--informat")) { 
    cfg.fmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--informat"));
    if(cfg.fmt == eslSQFILE_UNKNOWN) cm_Fail("Can't recognize sequence file format: %s. valid options are: fasta, embl, genbank, ddbj, uniprot, stockholm, or pfam\n", esl_opt_GetString(go, "--informat"));
  } else
    cfg.fmt = eslSQFILE_UNKNOWN; /* autodetect sequence file format by default. */ 

  cfg.abc        = NULL;	           /* created in init_master_cfg() in masters, or in mpi_worker() in workers */
  if      (esl_opt_GetBoolean(go, "--rna")) cfg.abc_out = esl_alphabet_Create(eslRNA);
  else if (esl_opt_GetBoolean(go, "--dna")) cfg.abc_out = esl_alphabet_Create(eslDNA);
  else    cm_Fail("Can't determine output alphabet");
  cfg.dbsize     = 0;                      /* db size used for E-values */
  cfg.namewidth  = 0;                      /* max name length in database */
  cfg.nseq       = 0;                      /* number of seqs in database */
  cfg.ncm        = 0;                      /* what number CM we're on, updated in masters, stays 0 (irrelevant) for workers */
  cfg.cmfp       = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.do_rc      = (! esl_opt_GetBoolean(go, "--toponly")); 
  cfg.init_rci   = esl_opt_GetBoolean(go, "--bottomonly") ? 1 : 0; 
  cfg.avg_hit_len= 0.;

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
      int              status;
      char             errbuf[cmERRBUFSIZE];
      if( esl_opt_IsOn(go, "--forecast")) cm_Fail("--forecast is incompatible with --mpi.");
      cfg.do_mpi     = TRUE;
      MPI_Init(&argc, &argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &(cfg.my_rank));
      MPI_Comm_size(MPI_COMM_WORLD, &(cfg.nproc));

      if(cfg.nproc == 1) cm_Fail("MPI mode, but only 1 processor running... (did you execute mpirun?)");

      if (cfg.my_rank > 0)  { status = mpi_worker(go, &cfg); }
      else {
	cm_banner(stdout, argv[0], banner);
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
      cm_banner(stdout, argv[0], banner);
      serial_master(go, &cfg);
      esl_stopwatch_Stop(w);
    }

  /* Clean up the shared cfg. 
   */
  if (cfg.my_rank == 0) {
    if ( esl_opt_IsOn(go, "-o")) { 
      printf("# Search results saved in file %s.\n", esl_opt_GetString(go, "-o"));
      fclose(cfg.ofp); 
    }
    if (cfg.cmfp      != NULL) CMFileClose(cfg.cmfp);
    if (cfg.sqfp      != NULL) esl_sqfile_Close(cfg.sqfp);
    if (cfg.tfp       != NULL) { 
      fclose(cfg.tfp);
      printf("# Tabular version of hit list saved in file %s.\n", esl_opt_GetString(go, "--tabfile"));
    }
    if (cfg.gcfp      != NULL) printf("# GC content stats of %s saved in file %s.\n", cfg.sqfile, esl_opt_GetString(go, "--gcfile"));
  }
  if (cfg.abc       != NULL) esl_alphabet_Destroy(cfg.abc);
  if (cfg.abc_out   != NULL) esl_alphabet_Destroy(cfg.abc_out);
  esl_getopts_Destroy(go);
  if (cfg.my_rank == 0) { 
    printf("#\n");
    esl_stopwatch_Display(stdout, w, "# CPU time: ");
  }
  esl_stopwatch_Destroy(w);
  return 0;
}

/* init_master_cfg()
 * Called by masters, mpi or serial.
 * Already set:
 *    cfg->cmfile      - command line arg 1
 *    cfg->sqfile      - command line arg 2
 *    cfg->fmt         - format of output file
 * Allocates/Sets: 
 *    cfg->sqfp        - open sequence file                
 *    cfg->cmfp        - open CM file                
 *
 * Errors in the MPI master here are considered to be "recoverable",
 * in the sense that we'll try to delay output of the error message
 * until we've cleanly shut down the worker processes. Therefore
 * errors return (code, errbuf) by the ESL_FAIL mech.
 */
static int
init_master_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  int status;

  /* open output file, or set to stdout if none */
  if (esl_opt_GetString(go, "-o") != NULL) {
    if ((cfg->ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) 
      ESL_FAIL(eslFAIL, errbuf, "Failed to open -o output file %s\n", esl_opt_GetString(go, "-o"));
  } else cfg->ofp = stdout;

  /* open input sequence file */
  status = esl_sqfile_Open(cfg->sqfile, cfg->fmt, NULL, &(cfg->sqfp));
  if (status == eslENOTFOUND)    ESL_FAIL(status, errbuf, "File %s doesn't exist or is not readable\n", cfg->sqfile);
  else if (status == eslEFORMAT) ESL_FAIL(status, errbuf, "Couldn't determine format of sequence file %s\n", cfg->sqfile);
  else if (status == eslEINVAL)  ESL_FAIL(status, errbuf, "Can’t autodetect stdin or .gz."); 
  else if (status != eslOK)      ESL_FAIL(status, errbuf, "Sequence file open failed with error %d\n", status);
  cfg->fmt = cfg->sqfp->format;

  /* if we need to, open the SSI index, if it exists */
  if(esl_opt_IsOn(go, "--sseq") || esl_opt_IsOn(go, "--eseq")) { 
    status = esl_sqfile_OpenSSI(cfg->sqfp, NULL);
    if(status != eslOK) { 
      if      (status == eslEFORMAT)   ESL_FAIL(status, errbuf, "SSI index is in incorrect format\n");
      else if (status == eslERANGE)    ESL_FAIL(status, errbuf, "SSI index is in 64-bit format and we can't read it\n");
      else if (status != eslENOTFOUND) ESL_FAIL(status, errbuf, "Failed to open existing SSI index\n");
    }
    /* we don't die if <seqfile>.ssi does not exist (in which case eslENOTFOUND is returned) */
  }

  /* Set the sqfile alphabet as RNA, if it's DNA we're fine. 
   * If it's not RNA nor DNA, we can't deal with it anyway,
   * so we're hardcoded to RNA.
   */
  cfg->abc = esl_alphabet_Create(eslRNA);
  if(cfg->abc == NULL) ESL_FAIL(status, errbuf, "Failed to create alphabet for sequence file");
  esl_sqfile_SetDigital(cfg->sqfp, cfg->abc);

  /* GetDBSize() reads all sequences, rewinds seq file and returns db size */
  if((status = GetDBSize(cfg->sqfp, errbuf, 
			 (esl_opt_IsOn(go, "--sseq")) ? esl_opt_GetInteger(go, "--sseq")-1 : -1,
			 (esl_opt_IsOn(go, "--eseq")) ? esl_opt_GetInteger(go, "--eseq")-1 : -1,
			 &(cfg->dbsize), &(cfg->nseq), &(cfg->namewidth))) != eslOK) return status;  
  if((! esl_opt_GetBoolean(go, "--toponly")) && (! esl_opt_GetBoolean(go, "--bottomonly"))) cfg->dbsize *= 2;

  /* overwrite dbsize if -Z enabled */
  if( esl_opt_IsOn(go, "-Z")) cfg->dbsize = (long) (esl_opt_GetReal(go, "-Z") * 1000000.); /* convert to Mb then to a long */

  /* if nec, open output file for --gcfile, and print to it */
  if ( esl_opt_IsOn(go, "--gcfile")) { 
    if ((cfg->gcfp = fopen(esl_opt_GetString(go, "--gcfile"), "w")) == NULL) 
      ESL_FAIL(eslFAIL, errbuf, "Failed to open --gcfile output file %s\n", esl_opt_GetString(go, "--gcfile"));
    if((status = dump_gc_info(go, cfg, errbuf)) != eslOK) return status;
    fclose(cfg->gcfp);
  }

  /* if nec, open output file for --tabfile */
  if ( esl_opt_IsOn(go, "--tabfile")) { 
    if ((cfg->tfp = fopen(esl_opt_GetString(go, "--tabfile"), "w")) == NULL) 
      ESL_FAIL(eslFAIL, errbuf, "Failed to open --tabfile output file %s\n", esl_opt_GetString(go, "--tabfile"));
  }

  /* open CM file */
  if ((cfg->cmfp = CMFileOpen(cfg->cmfile, NULL)) == NULL)
    ESL_FAIL(eslFAIL, errbuf, "Failed to open covariance model save file %s\n", cfg->cmfile);

  return eslOK;
}

/* serial_master()
 * The serial version of cmsearch.
 * 
 * 
 * A master can only return if it's successful. All errors are handled immediately and fatally with cm_Fail().
 */
static void
serial_master(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int            status;
  int            status2;
  char           errbuf[cmERRBUFSIZE];
  CM_t          *cm = NULL;
  CMConsensus_t *cons = NULL;     /* precalculated consensus info for display purposes */
  int            using_e_cutoff;
  int            using_sc_cutoff;
  int            rci;
  dbseq_t       *dbseq = NULL;
  int            do_top;
  float         *cm_surv_fractA  = NULL; /* 0..n..cm->si->nrounds fraction of db surviving round n for current CM */
  float         *seq_surv_fractA = NULL; /* 0..n..cm->si->nrounds fraction of db surviving round n for current seq */
  int           *cm_nhitsA = NULL;      /* 0..n..cm->si->nrounds number of hits reported for round n for current CM */
  int           *seq_nhitsA = NULL;     /* 0..n..cm->si->nrounds number of hits reported for round n for current seq */
  int            n, h;
  double         cm_psec;               /* predicted number of seconds for current CM versus full DB */
  double         total_psec = 0.;       /* predicted number of seconds for all CMs versus full DB */
  char           time_buf[128];	        /* for printing predicted time if --forecast only */
  ESL_STOPWATCH *w  = esl_stopwatch_Create();
  int            cm_namewidth;          /* length for printing model name field to tab file */
  char          *namedashes = NULL;     /* string of dashes for underlining 'target name' column in tab output */
  char          *cm_namedashes = NULL;  /* string of dashes for underlining 'model name' column in tab output */
  int            ni;                    /* index for filling dashes strings */
  int            used_at_least_one_cm = FALSE; /* only used if --cm-idx and --cm-name options are enabled */
  int            seqidx = 0;            /* index of current sequence in the sequence file (not nec number of seqs we've searched (if --sseq)) */
  
  if(w == NULL) cm_Fail("serial_master(): memory error, stopwatch not created.\n");

  if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);
  /*if ((status = init_shared_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);*/
  if ((status = print_run_info (go, cfg, errbuf))  != eslOK) cm_Fail(errbuf);
  do_top = (cfg->init_rci == 0) ? TRUE : FALSE; 

  /* create namedashes string, only used if --tabfile */
  ESL_ALLOC(namedashes, sizeof(char) * (cfg->namewidth+1));
  namedashes[cfg->namewidth] = '\0';
  for(ni = 0; ni < cfg->namewidth; ni++) namedashes[ni] = '-';

  while ((status = CMFileRead(cfg->cmfp, errbuf, &(cfg->abc), &cm)) == eslOK)
    {
      if (cm == NULL) cm_Fail("Failed to read CM from %s -- file corrupt?\n", cfg->cmfile);
      if((! (cm->flags & CMH_EXPTAIL_STATS)) && esl_opt_IsOn(go, "--forecast")) cm_Fail("--forecast only works with calibrated CM files. Run cmcalibrate (please)."); 
      /* potentially overwrite lambdas in cm->stats */
      if ( esl_opt_IsOn(go, "--lambda")) if((status = overwrite_lambdas(go, cfg, cm, errbuf)) != eslOK) cm_Fail(errbuf);
      cfg->ncm++;

      if(! esl_opt_IsDefault(go, "--cm-idx")) { 
	if(cfg->ncm != esl_opt_GetInteger(go, "--cm-idx")) { FreeCM(cm); continue; }
      }	
      if(! esl_opt_IsDefault(go, "--cm-name")) { 
	if(strcmp(cm->name, esl_opt_GetString(go, "--cm-name")) != 0) { FreeCM(cm); continue; }
      }	
      used_at_least_one_cm = TRUE;

      /* create cm_namedashes string, only used if --tabfile */
      cm_namewidth = ESL_MAX(strlen(cm->name), strlen("model name"));
      if(cm_namedashes != NULL) free(cm_namedashes); 
      ESL_ALLOC(cm_namedashes, sizeof(char) * (cm_namewidth+1));
      cm_namedashes[cm_namewidth] = '\0';
      for(ni = 0; ni < cm_namewidth; ni++) cm_namedashes[ni] = '-';

      /* initialize the flags/options/params and configuration of the CM */
      if((  status = initialize_cm(go, cfg, errbuf, cm))                    != eslOK) cm_Fail(errbuf);
      if((  status = CreateCMConsensus(cm, cfg->abc_out, 3.0, 1.0, &cons))  != eslOK) cm_Fail(errbuf);
      if(cm->flags & CMH_EXPTAIL_STATS) { 
	if((status = cm_GetAvgHitLen(cm, errbuf, &(cfg->avg_hit_len)))      != eslOK) cm_Fail(errbuf);
	if((status = UpdateExpsForDBSize(cm, errbuf, cfg->dbsize))          != eslOK) cm_Fail(errbuf);
	if((status = set_searchinfo_for_calibrated_cm(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);
      }
      else { if((status = set_searchinfo_for_uncalibrated_cm(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf); }

      ESL_ALLOC(cm_surv_fractA, sizeof(float) * (cm->si->nrounds+1));
      ESL_ALLOC(cm_nhitsA,      sizeof(int) * (cm->si->nrounds+1));
      esl_vec_FSet(cm_surv_fractA, (cm->si->nrounds+1), 0.);
      esl_vec_ISet(cm_nhitsA,     (cm->si->nrounds+1), 0);
      if(cm->flags & CMH_EXPTAIL_STATS) { if((status = print_searchinfo_for_calibrated_cm  (go, cfg, errbuf, cm, NULL, NULL, 0., 0., &cm_psec)) != eslOK) cm_Fail(errbuf); }
      else                              { if((status = print_searchinfo_for_uncalibrated_cm(go, cfg, errbuf, cm, NULL, NULL, 0.)) != eslOK) cm_Fail(errbuf); }

      if( esl_opt_IsOn(go, "--forecast")) { /* special mode, we don't do the search, just print the predicted timings */
	total_psec += cm_psec;
	free(cm_surv_fractA);
	free(cm_nhitsA);
	continue;
      }

      fprintf(cfg->ofp, "\nCM: %s\n", cm->name);
      if(cfg->tfp != NULL) { 
	fprintf(cfg->tfp, "#\n");
        fprintf(cfg->tfp, "# CM: %s\n", cm->name);
	/*fprintf(cfg->tfp, "# Predicted average hit length: %.2f\n", cfg->avg_hit_len);
	  fprintf(cfg->tfp, "# CM->W: %d (subtract (W-1) from stop and add (W-1) to start, and merge overlapping hits to simulate filter)\n", cm->W);*/
	fprintf(cfg->tfp, "# %-*s  %-*s  %22s  %12s  %8s  %8s  %3s\n", cm_namewidth, "", cfg->namewidth, "", "target coord", "query coord", "", "", "");
	fprintf(cfg->tfp, "# %-*s  %-*s  %22s  %12s  %8s  %8s  %3s\n", cm_namewidth, "", cfg->namewidth, "", "----------------------", "------------", "", "", "");
	fprintf(cfg->tfp, "# %-*s  %-*s  %10s  %10s  %5s  %5s  %8s  %8s  %3s\n", cm_namewidth, "model name", cfg->namewidth, "target name", "start", "stop", "start", "stop", "bit sc", "E-value", "GC%");
	fprintf(cfg->tfp, "# %-*s  %-*s  %10s  %10s  %5s  %5s  %8s  %8s  %3s\n", cm_namewidth, cm_namedashes, cfg->namewidth, namedashes, "----------", "----------", "-----", "-----", "--------", "--------", "---");
      }
      using_e_cutoff  = (cm->si->cutoff_type[cm->si->nrounds] == E_CUTOFF)     ? TRUE : FALSE;
      using_sc_cutoff = (cm->si->cutoff_type[cm->si->nrounds] == SCORE_CUTOFF) ? TRUE : FALSE;

      seqidx = 0;
      /* If nec, position the file to the first seq specified by --sseq */
      if(esl_opt_IsOn(go, "--sseq")) { 
	if((status = PositionSqFileByNumber(cfg->sqfp, esl_opt_GetInteger(go, "--sseq")-1, errbuf)) != eslOK) cm_Fail(errbuf);
	seqidx = esl_opt_GetInteger(go, "--sseq")-1;
      }
	 
      esl_stopwatch_Start(w);
      while ((status2 = read_next_search_seq(cfg->abc, cfg->sqfp, cfg->do_rc, &dbseq)) == eslOK)
	{
	  for(rci = cfg->init_rci; rci <= cfg->do_rc; rci++) {
	    /*printf("SEARCHING >%s %d\n", dbseq->sq[rci]->name, rci);*/
	    if ((status = ProcessSearchWorkunit(cm, errbuf, dbseq->sq[rci]->dsq, dbseq->sq[rci]->n, &dbseq->results[rci], esl_opt_GetReal(go, "--mxsize"), cfg->my_rank, &seq_surv_fractA, &seq_nhitsA)) != eslOK) cm_Fail(errbuf);
	    for(n = 0; n < cm->si->nrounds; n++) { 
	      cm_surv_fractA[n] += (dbseq->sq[rci]->n * seq_surv_fractA[n]);
	      cm_nhitsA[n]      += seq_nhitsA[n];
	    }
	    free(seq_surv_fractA);
	    free(seq_nhitsA);
	    RemoveOverlappingHits(dbseq->results[rci], 1, dbseq->sq[rci]->n);

	    /* write the final round nhits and surv_fract with values reflecting what we'll print, after overlaps have been removed */
	    cm_nhitsA[cm->si->nrounds] += dbseq->results[rci]->num_results;
	    for(h = 0; h < dbseq->results[rci]->num_results; h++) cm_surv_fractA[cm->si->nrounds] += fabs( (float) (dbseq->results[rci]->data[h].stop - dbseq->results[rci]->data[h].start + 1));
	  }
	  PrintResults (cm, cfg->ofp, cfg->tfp, cm->si, cfg->abc_out, cons, dbseq, do_top, cfg->do_rc, esl_opt_GetBoolean(go, "-x"), esl_opt_GetBoolean(go, "-v"), cfg->namewidth);
	  for(rci = 0; rci <= cfg->do_rc; rci++) { /* we can free results for top strand even if cfg->init_rci is 1 due to --bottomonly */
	    FreeResults(dbseq->results[rci]);
	    esl_sq_Destroy(dbseq->sq[rci]);
	  }
	  free(dbseq);
	  seqidx++;
	  if((esl_opt_IsOn(go, "--eseq")) && (seqidx == esl_opt_GetInteger(go, "--eseq"))) break; /* other way out of the loop */
 	}
      esl_stopwatch_Stop(w);
      if ((esl_opt_IsOn(go, "--eseq")) && (status2 != eslOK)) { /* we ran out of seqs too early, we never got to --eseq */
	cm_Fail("Ran out of seqs before getting to final seq %d (from --eseq)", esl_opt_GetInteger(go, "--eseq"));
      }
      else if ((! esl_opt_IsOn(go, "--eseq")) && status2 != eslEOF) cm_Fail("Parse failed!, file %s:\n%s", 
									    cfg->sqfp->filename, esl_sqfile_GetErrorBuf(cfg->sqfp));

      /* convert cm_surv_fractA[] values from residue counts into fractions */
      for(n = 0; n <= cm->si->nrounds; n++) cm_surv_fractA[n] /= (double) (cfg->dbsize);
      if(cm->flags & CMH_EXPTAIL_STATS) { if((status = print_searchinfo_for_calibrated_cm  (go, cfg, errbuf, cm, cm_surv_fractA, cm_nhitsA, w->elapsed, cm_psec, NULL)) != eslOK) cm_Fail(errbuf); }
      else                              { if((status = print_searchinfo_for_uncalibrated_cm(go, cfg, errbuf, cm, cm_surv_fractA, cm_nhitsA, w->elapsed)) != eslOK) cm_Fail(errbuf); }
      fprintf(cfg->ofp, "//\n");
      FreeCM(cm);
      FreeCMConsensus(cons);
      free(cm_surv_fractA);
      free(cm_nhitsA);
      esl_sqfile_Position(cfg->sqfp, (off_t) 0); /* we may be searching this file again with another CM */
    }
  if(status != eslEOF) cm_Fail(errbuf);
  if(! used_at_least_one_cm) { 
    if(! esl_opt_IsDefault(go, "--cm-idx"))  { cm_Fail("--cm-idx %d enabled, but only %d CMs in the cmfile.\n", esl_opt_GetInteger(go, "--cm-idx"), cfg->ncm); }
    if(! esl_opt_IsDefault(go, "--cm-name")) { cm_Fail("--cm-name %s enabled, but no CM named %s exists in the cmfile.\n", esl_opt_GetString(go, "--cm-name"), esl_opt_GetString(go, "--cm-name")); }
  }

  if(cfg->ncm > 1 &&  esl_opt_IsOn(go, "--forecast")) { 
    fprintf(stdout, "#\n");
    fprintf(stdout, "# %20s\n", "predicted total time");
    fprintf(stdout, "# %20s\n", "--------------------");
    FormatTimeString(time_buf, total_psec, FALSE);
    fprintf(stdout, "  %20s\n", time_buf);
  }
      
  if(namedashes != NULL)    free(namedashes);
  if(cm_namedashes != NULL) free(cm_namedashes);
  esl_stopwatch_Destroy(w);
  return;

 ERROR:
  cm_Fail("serial_master: memory allocation error.");
  return; /* NEVERREACHED */
}

#ifdef HAVE_MPI
/* mpi_master()
 * The MPI version of cmsearch
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
  int      using_e_cutoff; 
  int      using_sc_cutoff; 
  int      wi_error = 0;                /* worker index that sent back an error message, if an error occurs */

  CM_t *cm;
  CMConsensus_t *cons = NULL;     /* precalculated consensus info for display purposes */

  int si      = 0;        /* sequence index */
  int si_recv = 1;        /* sequence index of the sequence we've just received results for from a worker */

  /* properties of the workers, indexed 1..wi..nproc-1 */
  int *silist = NULL;     /* [0..wi..nproc-1], the sequence index worker wi is working on */
  int in_rc = FALSE;      /* are we currently on the reverse complement? */
  int *rclist = NULL;     /* [0..wi..nproc-1] 0 if worker wi is searching top strand, 1 if wi is searching bottom strand */
  int rci;                /* index that ranges from 0 to 1 */
  int seqpos = 1;         /* sequence position in the current sequence */
  int *seqposlist = NULL; /* [0..wi..nproc-1] the first position of the sequence that worker wi is searching */
  int len;                /* length of chunk */
  int *lenlist = NULL;    /* [0..wi..nproc-1] length of chunk worker wi is searching */

  /* properties of the sequences currently being worked on, we can have at most 1 per worker, so these are of size 
   * cfg->nproc, but indexed by si = 0..nproc-2, cfg->nproc-1 is never used. */
  int *sentlist = NULL;   /* [0..si..nproc-1] TRUE if all chunks for sequence index si have been sent, FALSE otherwise */
  int ndbseq = 0;         /* ndbseq is the number of currently active sequences, we can read a new seq IFF ndbseq < (cfg->nproc-1) */
  dbseq_t **dbseqlist= NULL; /* pointers to the dbseq_t objects that hold the actual sequence data, and the results data */
  dbseq_t  *dbseq = NULL;    /* a database sequence */
  double    cm_psec;                /* predicted number of seconds for current CM versus full DB (ON MASTER PROC BUT WE ASSUME
				     * OTHER PROCS ARE THE SAME SPEED!) */
  float    *cm_surv_fractA = NULL;  /* 0..n..cm->si->nrounds fraction of db that survived round n for current CM */
  int      *cm_nhitsA = NULL;       /* 0..n..cm->si->nrounds number of hits reported for round n for current CM */
  int       cm_namewidth;           /* length for printing model name field to tab file */
  char     *namedashes = NULL;      /* string of dashes for underlining 'target name' column in tab output */
  char     *cm_namedashes = NULL;   /* string of dashes for underlining 'model name' column in tab output */
  int       ni;                     /* index for filling dashes strings */
  int       used_at_least_one_cm = FALSE; /* only used if --cm-idx and --cm-name options are enabled */
  int       seqidx = 0;             /* index of current sequence in the sequence file (not nec number of seqs we've searched (if --sseq)) */

  ESL_STOPWATCH *w  = esl_stopwatch_Create();
  if(w == NULL) cm_Fail("mpi_master(): memory error, stopwatch not created.\n");  

  MPI_Status mpistatus; 
  int      n, h;

  int need_seq = TRUE;
  int chunksize;
  search_results_t *worker_results;
  float            *worker_surv_fractA = NULL; /* 0..n..cm->si->nrounds fraction of db surviving round n for current seq */
  int              *worker_nhitsA      = NULL; /* 0..n..cm->si->nrounds num hits surviving round n for current seq */

  /* Master initialization: including, figure out the alphabet type.
   * If any failure occurs, delay printing error message until we've shut down workers.
   */
  if (xstatus == eslOK) { if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) xstatus = status; }
  if (xstatus == eslOK) { bn = 4096; if ((buf = malloc(sizeof(char) * bn)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((silist     = malloc(sizeof(int) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((rclist     = malloc(sizeof(int) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((seqposlist = malloc(sizeof(int) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((lenlist    = malloc(sizeof(int) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((dbseqlist  = malloc(sizeof(dbseq_t *) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((sentlist   = malloc(sizeof(int) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((status = print_run_info(go, cfg, errbuf))  != eslOK) xstatus = status; }

  MPI_Bcast(&xstatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (xstatus != eslOK) return xstatus; /* errbuf was filled above */
  ESL_DPRINTF1(("MPI master is initialized\n"));

  /* create namedashes string, only used if --tabfile */
  ESL_ALLOC(namedashes, sizeof(char) * (cfg->namewidth+1));
  namedashes[cfg->namewidth] = '\0';
  for(ni = 0; ni < cfg->namewidth; ni++) namedashes[ni] = '-';

  for (wi = 0; wi < cfg->nproc; wi++) 
  { 
    silist[wi] = rclist[wi] = seqposlist[wi] = lenlist[wi] = -1;
    dbseqlist[wi] = NULL;
    sentlist[wi] = FALSE;
  }
  /* Worker initialization:
   * Because we've already successfully initialized the master before we start
   * initializing the workers, we don't expect worker initialization to fail;
   * so we just receive a quick OK/error code reply from each worker to be sure,
   * and don't worry about an informative message. 
   */
  MPI_Bcast(&(cfg->dbsize), 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Reduce(&xstatus, &status, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  if (status != eslOK) cm_Fail("One or more MPI worker processes failed to initialize.");
  ESL_DPRINTF1(("%d workers are initialized\n", cfg->nproc-1));

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
      ESL_DPRINTF1(("MPI master read CM number %d\n", cfg->ncm));
      if(! esl_opt_IsDefault(go, "--cm-idx")) { 
	if(cfg->ncm != esl_opt_GetInteger(go, "--cm-idx")) { FreeCM(cm); continue; }
      }	
      if(! esl_opt_IsDefault(go, "--cm-name")) { 
	if(strcmp(cm->name, esl_opt_GetString(go, "--cm-name")) != 0) { FreeCM(cm); continue; }
      }	
      used_at_least_one_cm = TRUE;

      /* potentially overwrite lambdas in cm->stats */
      if (esl_opt_IsOn(go, "--lambda")) if((status = overwrite_lambdas(go, cfg, cm, errbuf)) != eslOK) cm_Fail(errbuf);

      if((status = cm_master_MPIBcast(cm, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI broadcast CM failed.");

      /* create cm_namedashes string, only used if --tabfile */
      cm_namewidth = ESL_MAX(strlen(cm->name), strlen("model name"));
      if(cm_namedashes != NULL) free(cm_namedashes); 
      ESL_ALLOC(cm_namedashes, sizeof(char) * (cm_namewidth+1));
      cm_namedashes[cm_namewidth] = '\0';
      for(ni = 0; ni < cm_namewidth; ni++) cm_namedashes[ni] = '-';
      
      /* initialize the flags/options/params of the CM */
      if((status   = initialize_cm(go, cfg, errbuf, cm))                    != eslOK) cm_Fail(errbuf);
      if((status   = cm_GetAvgHitLen(cm, errbuf, &(cfg->avg_hit_len)))      != eslOK) cm_Fail(errbuf);
      if((status   = CreateCMConsensus(cm, cfg->abc_out, 3.0, 1.0, &cons))  != eslOK) cm_Fail(errbuf);
      if(cm->flags & CMH_EXPTAIL_STATS) {
	if((status = UpdateExpsForDBSize(cm, errbuf, cfg->dbsize))          != eslOK) cm_Fail(errbuf);
	if((status = set_searchinfo_for_calibrated_cm(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);
      }
      else { if((status = set_searchinfo_for_uncalibrated_cm(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf); }

      using_e_cutoff  = (cm->si->cutoff_type[cm->si->nrounds] == E_CUTOFF)     ? TRUE : FALSE;
      using_sc_cutoff = (cm->si->cutoff_type[cm->si->nrounds] == SCORE_CUTOFF) ? TRUE : FALSE;

      ESL_ALLOC(cm_surv_fractA, sizeof(float) * (cm->si->nrounds+1));
      ESL_ALLOC(cm_nhitsA,      sizeof(int) * (cm->si->nrounds+1));
      esl_vec_FSet(cm_surv_fractA, (cm->si->nrounds+1), 0.);
      esl_vec_ISet(cm_nhitsA,      (cm->si->nrounds+1), 0);
      if(cm->flags & CMH_EXPTAIL_STATS) { if((status = print_searchinfo_for_calibrated_cm  (go, cfg, errbuf, cm, NULL, NULL, 0., 0., &cm_psec)) != eslOK) cm_Fail(errbuf); }
      else                              { if((status = print_searchinfo_for_uncalibrated_cm(go, cfg, errbuf, cm, NULL, NULL, 0.)) != eslOK) cm_Fail(errbuf); }

      fprintf(cfg->ofp, "CM: %s\n", cm->name);
      if(cfg->tfp != NULL) { 
	fprintf(cfg->tfp, "#\n");
        fprintf(cfg->tfp, "# CM: %s\n", cm->name);
	fprintf(cfg->tfp, "# %-*s  %-*s  %22s  %12s  %8s  %8s  %3s\n", cm_namewidth, "", cfg->namewidth, "", "target coord", "query coord", "", "", "");
	fprintf(cfg->tfp, "# %-*s  %-*s  %22s  %12s  %8s  %8s  %3s\n", cm_namewidth, "", cfg->namewidth, "", "----------------------", "------------", "", "", "");
	fprintf(cfg->tfp, "# %-*s  %-*s  %10s  %10s  %5s  %5s  %8s  %8s  %3s\n", cm_namewidth, "model name", cfg->namewidth, "target name", "start", "stop", "start", "stop", "bit sc", "E-value", "GC%");
	fprintf(cfg->tfp, "# %-*s  %-*s  %10s  %10s  %5s  %5s  %8s  %8s  %3s\n", cm_namewidth, cm_namedashes, cfg->namewidth, namedashes, "----------", "----------", "-----", "-----", "--------", "--------", "---");
      }

      /* reset vars for searching with current CM */
      wi = 1;
      ndbseq = 0;
      need_seq = TRUE;
      have_work = TRUE;	/* TRUE while work remains  */
      seqpos = 1;
      in_rc = FALSE;
      seqidx = 0;
      
      /* If nec, position the file to the first seq specified by --sseq */
      if(esl_opt_IsOn(go, "--sseq")) { 
	if((status = PositionSqFileByNumber(cfg->sqfp, esl_opt_GetInteger(go, "--sseq")-1, errbuf)) != eslOK) cm_Fail(errbuf);
	seqidx = esl_opt_GetInteger(go, "--sseq")-1;
      }

      esl_stopwatch_Start(w);
      while (have_work || nproc_working)
	{
	  if (need_seq) 
	    {
	      need_seq = FALSE;
	      /* read a new seq */
	      if     ((esl_opt_IsOn(go, "--eseq")) && (seqidx == esl_opt_GetInteger(go, "--eseq"))) 
		{ 
		  /* we've already read and processed the final seq we'll search */
		  have_work = FALSE;
		}
	      else if((status = read_next_search_seq(cfg->abc, cfg->sqfp, cfg->do_rc, &dbseq)) == eslOK) 
		{
		  ndbseq++;
		  ESL_DASSERT1((ndbseq < cfg->nproc));

		  dbseq->chunks_sent = 0;
		  dbseq->alignments_sent = -1;     /* None sent yet */
		  for(rci = 0; rci <= cfg->do_rc; rci++) {
		    dbseq->results[rci] = CreateResults(INIT_RESULTS);
		  }
		  in_rc = (cfg->init_rci == 0) ? FALSE : TRUE; /* if --bottomonly --> cfg->init_rci = 1, and we only search bottom strand */
		  seqpos = 1;
		  
		  si = 0;
		  while(dbseqlist[si] != NULL) si++;
		  ESL_DASSERT1((si < cfg->nproc));
		  dbseqlist[si] = dbseq;
		  sentlist[si]  = FALSE;
		  have_work = TRUE;
		  chunksize = DetermineSeqChunksize(cfg->nproc, dbseq->sq[0]->n, cm->W);
		  seqidx++;
		  ESL_DPRINTF1(("L: %ld chunksize: %d\n", dbseq->sq[0]->n, chunksize));
		}
	      else if(status == eslEOF) have_work = FALSE;
	      else goto ERROR;
	    }
	
	  if ((have_work && nproc_working == cfg->nproc-1) || (!have_work && nproc_working > 0))
	    {
	      /* we're waiting to receive */
	      if (MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &mpistatus) != 0) cm_Fail("mpi probe failed");
	      if (MPI_Get_count(&mpistatus, MPI_PACKED, &n)                != 0) cm_Fail("mpi get count failed");
	      wi = mpistatus.MPI_SOURCE;
	      ESL_DPRINTF1(("MPI master sees a result of %d bytes from worker %d\n", n, wi));
	      
	      if (n > bn) {
		if ((buf = realloc(buf, sizeof(char) * n)) == NULL) cm_Fail("reallocation failed");
		bn = n; 
	      }
	      if (MPI_Recv(buf, bn, MPI_PACKED, wi, 0, MPI_COMM_WORLD, &mpistatus) != 0) cm_Fail("mpi recv failed");
	      ESL_DPRINTF1(("MPI master has received the buffer\n"));
	      
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
		      si_recv = silist[wi];
		      ESL_DPRINTF1(("MPI master sees that the result buffer contains search results (si_recv:%d)\n", si_recv));

		      /* first unpack the surv_fractA and nhitsA, that give survival fractions and num hits surviving for each round */
		      ESL_ALLOC(worker_surv_fractA, sizeof(float) * (cm->si->nrounds+1));
		      ESL_ALLOC(worker_nhitsA,      sizeof(int)   * (cm->si->nrounds+1));
		      if (MPI_Unpack(buf, bn, &pos, worker_surv_fractA, (cm->si->nrounds+1), MPI_FLOAT, MPI_COMM_WORLD) != 0) cm_Fail("mpi unpack failed");
		      if (MPI_Unpack(buf, bn, &pos, worker_nhitsA,      (cm->si->nrounds+1), MPI_INT,   MPI_COMM_WORLD) != 0) cm_Fail("mpi unpack failed");

		      if ((status = cm_search_results_MPIUnpack(buf, bn, &pos, MPI_COMM_WORLD, &worker_results)) != eslOK)    cm_Fail("search results unpack failed");
		      ESL_DPRINTF1(("MPI master has unpacked search results\n"));
		      
		      /* update cm_surv_fractA[] and cm_nhitsA[] which holds number of residues and hits that survived each round of searching/filtering, don't update for final round, we do that later */
		      for(n = 0; n < cm->si->nrounds; n++) { 
			cm_surv_fractA[n] += (lenlist[wi] * worker_surv_fractA[n]);
			cm_nhitsA[n]      += worker_nhitsA[n];
		      }
		      free(worker_surv_fractA);

		      /* worker_results will be NULL if 0 results (hits) sent back */
		      int x;
		      if(worker_results != NULL) { 
			/* add results to dbseqlist[si_recv]->results[rclist[wi]] */
			if(! esl_opt_GetBoolean(go, "--noalign")) { 
			  for(x = 0; x < worker_results->num_results; x++) {
			    assert(worker_results->data[x].tr != NULL);
			    assert(worker_results->data[x].tr->n > 0);
			  }
			}
			AppendResults(worker_results, dbseqlist[si_recv]->results[rclist[wi]], seqposlist[wi]);
			/* careful, dbseqlist[si_recv]->results[rclist[wi]] now points to the traces and posterior code strings in worker_results->data,
			 * don't free those (don't use FreeResults(worker_results)) */
			free(worker_results->data);
			free(worker_results);
			worker_results = NULL;
		      }
		      dbseqlist[si_recv]->chunks_sent--;
		      if(sentlist[si_recv] && dbseqlist[si_recv]->chunks_sent == 0)
			{
			  for(rci = 0; rci <= cfg->do_rc; rci++) {
			    RemoveOverlappingHits(dbseqlist[si_recv]->results[rci], 1, dbseqlist[si_recv]->sq[rci]->n);

			    /* write the final round nhits and surv_fract with values reflecting what we'll print, after overlaps have been removed, this wasn't done in the earlier filter rounds
			     * b/c of recursive nature of implementation (each worker searched it's own survivors after filter was applied). THUS DIFFERENCES WILL EXIST BETWEEN MPI AND SERIAL 
			     * IN SURVIVAL FRACTIONS AND NHITS OF FILTER ROUNDS! */
			    cm_nhitsA[cm->si->nrounds] += dbseqlist[si_recv]->results[rci]->num_results;
			    for(h = 0; h < dbseqlist[si_recv]->results[rci]->num_results; h++) cm_surv_fractA[cm->si->nrounds] += fabs( (float) (dbseqlist[si_recv]->results[rci]->data[h].stop - dbseqlist[si_recv]->results[rci]->data[h].start + 1.));
			  }
			  PrintResults(cm, cfg->ofp, cfg->tfp, cm->si, cfg->abc_out, cons, dbseqlist[si_recv], TRUE, cfg->do_rc, esl_opt_GetBoolean(go, "-x"), esl_opt_GetBoolean(go, "-v"), cfg->namewidth);
			  for(rci = 0; rci <= cfg->do_rc; rci++) {
			    esl_sq_Destroy(dbseqlist[si_recv]->sq[rci]);
			    FreeResults(dbseqlist[si_recv]->results[rci]);
			  }
			  free(dbseqlist[si_recv]);
			  dbseqlist[si_recv] = NULL;
			  ndbseq--;
			}
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
	      /* send new search job */
	      len = (chunksize < (dbseqlist[si]->sq[0]->n - seqpos + 1)) ? chunksize : (dbseqlist[si]->sq[0]->n - seqpos + 1);
	      ESL_DPRINTF1(("MPI master is sending sequence i0..j0 %d..%d to search to worker %d\n", seqpos, seqpos+len-1, wi));
	      assert(seqpos > 0);
	      if ((status = cm_dsq_MPISend(dbseqlist[si]->sq[in_rc]->dsq+seqpos-1, len, wi, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI search job send failed");
	      
	      silist[wi]      = si;
	      seqposlist[wi]  = seqpos;
	      lenlist[wi]     = len;
	      rclist[wi]      = in_rc;
	      dbseqlist[si]->chunks_sent++;
	      
	      wi++;
	      nproc_working++;
	      
	      if(len == chunksize) seqpos += len - cm->W + 1;
	      else if(cfg->do_rc && !in_rc) {
		in_rc = TRUE;
		seqpos = 1; 
	      }
	      else {
		need_seq     = TRUE;
		sentlist[si] = TRUE; /* we've sent all chunks from this seq */
	      }
	    }
	}
      esl_stopwatch_Stop(w);
      ESL_DPRINTF1(("MPI master: done with this CM. Telling all workers\n"));
      for (wi = 1; wi < cfg->nproc; wi++) 
	if ((status = cm_dsq_MPISend(NULL, 0, wi, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("Shutting down a worker failed.");
      /* convert cm_surv_fractA[] values from residue counts into fractions */
      for(n = 0; n <= cm->si->nrounds; n++) cm_surv_fractA[n] /= (double) (cfg->dbsize);
      if(cm->flags & CMH_EXPTAIL_STATS) { if((status = print_searchinfo_for_calibrated_cm  (go, cfg, errbuf, cm, cm_surv_fractA, cm_nhitsA, w->elapsed, cm_psec, NULL)) != eslOK) cm_Fail(errbuf); }
      else                              { if((status = print_searchinfo_for_uncalibrated_cm(go, cfg, errbuf, cm, cm_surv_fractA, cm_nhitsA, w->elapsed)) != eslOK) cm_Fail(errbuf); }
      fprintf(cfg->ofp, "//\n");
      free(cm_surv_fractA);
      free(cm_nhitsA);
      FreeCM(cm);
      FreeCMConsensus(cons);
      esl_sqfile_Position(cfg->sqfp, (off_t) 0); /* we may be searching this file again with another CM */
    }

  /* On success or recoverable errors:
   * Shut down workers cleanly. 
   */
  ESL_DPRINTF1(("MPI master is done. Shutting down all the workers cleanly\n"));
  if((cm_master_MPIBcast(NULL, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI broadcast CM failed.");
  free(buf);
  
  if(namedashes != NULL)    free(namedashes);
  if(cm_namedashes != NULL) free(cm_namedashes);
  esl_stopwatch_Destroy(w);

  if     (xstatus != eslOK) { fprintf(stderr, "Worker: %d had a problem.\n", wi_error); return xstatus; }
  else if((! used_at_least_one_cm) && (! esl_opt_IsDefault(go, "--cm-idx")))  { ESL_FAIL(eslEINVAL, errbuf, "--cm-idx %d enabled, but only %d CMs in the cmfile.\n", esl_opt_GetInteger(go, "--cm-idx"), cfg->ncm); }
  else if((! used_at_least_one_cm) && (! esl_opt_IsDefault(go, "--cm-name"))) { ESL_FAIL(eslEINVAL, errbuf, "--cm-name %s enabled, but no CM named %s exists in the cmfile.\n", esl_opt_GetString(go, "--cm-name"), esl_opt_GetString(go, "--cm-name")); }
  else if(status != eslEOF) return status;  /* problem reading CM file */
  else                      return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "mpi_master() memory allocation error.");
  return eslOK; /* NOTREACHED */
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
  ESL_DSQ      *dsq = NULL;
  int           L;
  search_results_t *results = NULL;
  float         *surv_fractA = NULL; /* 0..n..cm->si->nrounds fraction of db surviving round n for current seq */
  int           *nhitsA      = NULL; /* 0..n..cm->si->nrounds number of hits surviving round n for current seq */
  char           errbuf[cmERRBUFSIZE];
  /* After master initialization: master broadcasts its status.
   */
  MPI_Bcast(&xstatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (xstatus != eslOK) return xstatus; /* master saw an error code; workers do an immediate normal shutdown. */
  ESL_DPRINTF1(("worker %d: sees that master has initialized\n", cfg->my_rank));
  
  /* Master now broadcasts worker initialization information (db size N) 
   * Workers returns their status post-initialization.
   * Initial allocation of wbuf must be large enough to guarantee that
   * we can pack an error result into it, because after initialization,
   * errors will be returned as packed (code, errbuf) messages.
   */
  MPI_Bcast(&(cfg->dbsize), 1, MPI_LONG, 0, MPI_COMM_WORLD);
  if (xstatus == eslOK) { wn = 4096;  if ((wbuf = malloc(wn * sizeof(char))) == NULL) xstatus = eslEMEM; }
  MPI_Reduce(&xstatus, &status, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD); /* everyone sends xstatus back to master */
  if (xstatus != eslOK) {
    if (wbuf != NULL) free(wbuf);
    return xstatus; /* shutdown; we passed the error back for the master to deal with. */
  }
  ESL_DPRINTF1(("worker %d: initialized N: %ld\n", cfg->my_rank, cfg->dbsize));
  
  /* source = 0 (master); tag = 0 */
  while ((status = cm_worker_MPIBcast(0, MPI_COMM_WORLD, &wbuf, &wn, &(cfg->abc), &cm)) == eslOK)
    {
      ESL_DPRINTF1(("Worker %d succesfully received CM, num states: %d num nodes: %d\n", cfg->my_rank, cm->M, cm->nodes));
      
      /* initialize the flags/options/params of the CM */
      if((status   = initialize_cm(go, cfg, errbuf, cm))                    != eslOK) goto ERROR;
      if((status   = cm_GetAvgHitLen(cm, errbuf, &(cfg->avg_hit_len)))      != eslOK) goto ERROR;
      if(cm->flags & CMH_EXPTAIL_STATS) {
	if((status = UpdateExpsForDBSize(cm, errbuf, cfg->dbsize))          != eslOK) goto ERROR;
	if((status = set_searchinfo_for_calibrated_cm(go, cfg, errbuf, cm)) != eslOK) goto ERROR;
      }
      else { if((status = set_searchinfo_for_uncalibrated_cm(go, cfg, errbuf, cm)) != eslOK) goto ERROR; }
      
      while((status = cm_dsq_MPIRecv(0, 0, MPI_COMM_WORLD, &wbuf, &wn, &dsq, &L)) == eslOK)
	{
	  ESL_DPRINTF1(("worker %d: has received search job, length: %d\n", cfg->my_rank, L));
	  if ((status = ProcessSearchWorkunit(cm, errbuf, dsq, L, &results, esl_opt_GetReal(go, "--mxsize"), cfg->my_rank, &surv_fractA, &nhitsA)) != eslOK) goto ERROR;
	  ESL_DPRINTF1(("worker %d: has gathered search results\n", cfg->my_rank));
	  
	  n = 0;
	  if (MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &sz) != 0) /* room for the status code */
	    ESL_XFAIL(eslESYS, errbuf, "mpi pack size failed"); 
	  n += sz;
	  if (MPI_Pack_size((cm->si->nrounds+1), MPI_FLOAT, MPI_COMM_WORLD, &sz) != 0) /* room for the surv_fractA array */
	    ESL_XFAIL(eslESYS, errbuf, "mpi pack size failed"); 
	  n += sz;
	  if (MPI_Pack_size((cm->si->nrounds+1), MPI_INT, MPI_COMM_WORLD, &sz) != 0) /* room for the nhitsA array */
	    ESL_XFAIL(eslESYS, errbuf, "mpi pack size failed"); 
	  n += sz;
	  if (cm_search_results_MPIPackSize(results, MPI_COMM_WORLD, &sz) != eslOK)
	    ESL_XFAIL(eslFAIL, errbuf, "cm_serch_results_MPIPackSize() call failed"); 
	  n += sz;  

	  if (n > wn) {
	    void *tmp;
	    ESL_RALLOC(wbuf, tmp, sizeof(char) * n);
	    wn = n;
	  }
	  ESL_DPRINTF1(("worker %d: has calculated the search results will pack into %d bytes\n", cfg->my_rank, n));
	  status = eslOK;

	  pos = 0;
	  if (MPI_Pack(&status, 1, MPI_INT, wbuf, wn, &pos, MPI_COMM_WORLD) != 0) 
	    ESL_XFAIL(eslESYS, errbuf, "mpi pack failed.");
	  if (MPI_Pack(surv_fractA, (cm->si->nrounds+1), MPI_FLOAT, wbuf, wn, &pos, MPI_COMM_WORLD) != 0) 
	    ESL_XFAIL(eslESYS, errbuf, "mpi pack failed.");
	  if (MPI_Pack(nhitsA, (cm->si->nrounds+1), MPI_INT,   wbuf, wn, &pos, MPI_COMM_WORLD) != 0) 
	    ESL_XFAIL(eslESYS, errbuf, "mpi pack failed.");
	  if (cm_search_results_MPIPack(results, wbuf, wn, &pos, MPI_COMM_WORLD) != eslOK)
	    ESL_XFAIL(eslFAIL, errbuf, "cm_search_results_MPIPack() call failed"); 
	  MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);
	  ESL_DPRINTF1(("worker %d: has sent results to master in message of %d bytes\n", cfg->my_rank, pos));

	  free(surv_fractA);
	  free(nhitsA);
	  FreeResults(results);
	  free(dsq);
	}
      if(status == eslEOD)ESL_DPRINTF1(("worker %d: has seen message to stop with this CM.\n", cfg->my_rank));
      else ESL_XFAIL(eslFAIL, errbuf, "within CM loop, unexpected status code: %d received from cm_dsq_MPIRecv()\n", status);

      FreeCM(cm);
      cm = NULL;
    }
  if (status == eslEOD) ESL_DPRINTF1(("Worker %d told CMs are done.\n", cfg->my_rank));
  else ESL_XFAIL(eslFAIL, errbuf, "outside CM loop, unexpected status code: %d received from cm_seqs_to_aln_MPIRecv()\n", status);
  
  if (wbuf != NULL) free(wbuf);
  return eslOK;

 ERROR:
  ESL_DPRINTF1(("worker %d: fails, is sending an error message, as follows:\n%s\n", cfg->my_rank, errbuf));
  pos = 0;
  MPI_Pack(&status, 1,                MPI_INT,  wbuf, wn, &pos, MPI_COMM_WORLD);
  MPI_Pack(errbuf,  cmERRBUFSIZE,    MPI_CHAR, wbuf, wn, &pos, MPI_COMM_WORLD);
  MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);

  /* if we get here this worker failed and sent an error message, now the master knows a worker
   * failed but it has to send the message to all other workers (besides this one) to abort so they 
   * can be shut down cleanly. As currently implemented, this means we have to wait here for that 
   * signal which comes in the form of a special 'empty' work packet that tells us we're done with
   * the current CM, and then a 'empty' CM broadcast that tells us we're done with all CMs in the file.
   */
  status = cm_dsq_MPIRecv(0, 0, MPI_COMM_WORLD, &wbuf, &wn, &dsq, &L);
  status = cm_worker_MPIBcast(0, MPI_COMM_WORLD, &wbuf, &wn, &(cfg->abc), &cm);
  /* status after each of the above calls should be eslEOD, but if it isn't we can't really do anything 
   * about it b/c we've already sent our error message, so in that scenario the MPI will break uncleanly 
   */
  return eslFAIL; /* recoverable error, master has error message and will print it */
}
#endif /*HAVE_MPI*/

/* initialize_cm()
 * Setup the CM based on the command-line options/defaults;
 * only set flags and a few parameters. ConfigCM() configures
 * the CM. 
 */
static int
initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int status;
  int use_hmmonly;
  int config_qdb;
  int hmm_W_set;
  int nstarts, nexits, nd;
  
  /* set up CM parameters that are option-changeable */
  cm->beta_qdb = esl_opt_GetReal(go, "--beta");
  cm->tau      = esl_opt_GetReal(go, "--tau");  /* this will be DEFAULT_TAU unless changed at command line */

  use_hmmonly = (esl_opt_GetBoolean(go, "--viterbi") || esl_opt_GetBoolean(go, "--forward")) ? TRUE : FALSE;
  config_qdb  = (use_hmmonly || esl_opt_GetBoolean(go, "--no-qdb")) ? FALSE : TRUE; 

  /* Update cm->config_opts and cm->align_opts based on command line options */

  /* config_opts */
  if(! esl_opt_GetBoolean(go, "-g")) { 
    cm->config_opts |= CM_CONFIG_LOCAL;
    cm->config_opts |= CM_CONFIG_HMMLOCAL;
    cm->config_opts |= CM_CONFIG_HMMEL;
  }
  /* config QDB for final round of search? yes, unless --no-qdb, --viterbi or --forward */
  if(config_qdb) cm->config_opts |= CM_CONFIG_QDB;

  /* search_opts */
  if(! use_hmmonly) 
    if(  esl_opt_GetBoolean(go, "--inside"))    cm->search_opts |= CM_SEARCH_INSIDE;
  if(  esl_opt_GetBoolean(go, "--noalign"))     cm->search_opts |= CM_SEARCH_NOALIGN;
  if(  esl_opt_GetBoolean(go, "--no-qdb"))      cm->search_opts |= CM_SEARCH_NOQDB;
  if(  esl_opt_GetBoolean(go, "--hbanded"))     cm->search_opts |= CM_SEARCH_HBANDED;
  if(  esl_opt_GetBoolean(go, "--aln2bands"))   cm->search_opts |= CM_SEARCH_HMMALNBANDS;
  if(  esl_opt_GetBoolean(go, "--null2"))       cm->search_opts |= CM_SEARCH_NULL2;
  if(! esl_opt_GetBoolean(go, "--no-null3"))    cm->search_opts |= CM_SEARCH_NULL3;
  if(  esl_opt_GetBoolean(go, "--viterbi"))  { 
    cm->search_opts |= CM_SEARCH_HMMVITERBI;
    cm->search_opts |= CM_SEARCH_NOQDB;
  }
  if(  esl_opt_GetBoolean(go, "--forward"))  { 
    cm->search_opts |= CM_SEARCH_HMMFORWARD;
    cm->search_opts |= CM_SEARCH_NOQDB;
  }

  /* align_opts, by default, DO NOT align with HMM bands */
  if(esl_opt_GetBoolean(go, "--aln-hbanded"))  { 
    cm->align_opts |= CM_ALIGN_HBANDED;
    /* these guys are only available if --aln-hbanded also enabled */
    if(esl_opt_GetBoolean(go, "--aln-optacc")) cm->align_opts |= CM_ALIGN_OPTACC;
    if(esl_opt_GetBoolean(go, "-p"))           cm->align_opts |= CM_ALIGN_POST;
  }
  else cm->align_opts |= CM_ALIGN_SMALL;

  /*******************************
   * Begin developer options block 
   *******************************/
  /* handle special developer's options, not recommend for normal users */
  if(  esl_opt_GetBoolean(go, "--rtrans"))      cm->flags       |= CM_RSEARCHTRANS;
  if(  esl_opt_GetBoolean(go, "--sums"))        cm->search_opts |= CM_SEARCH_SUMS;

  /* set aggregate local begin/end probs, set with --pbegin, --pend, defaults are DEFAULT_PBEGIN, DEFAULT_PEND */
  cm->pbegin = esl_opt_GetReal(go, "--pbegin");
  cm->pend   = esl_opt_GetReal(go, "--pend");
  /* possibly overwrite local begin probs such that all begin points are equiprobable (--pebegin) */
  if(esl_opt_GetBoolean(go, "--pebegin")) {
    nstarts = 0;
    for (nd = 2; nd < cm->nodes; nd++) 
      if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd || cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BIF_nd) 
	nstarts++;
    cm->pbegin = 1.- (1./(1+nstarts));
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
  /*******************************
   * End developer options block 
   *******************************/

  /* finally, configure the CM for search based on cm->config_opts and cm->align_opts.
   * set local mode, make cp9 HMM, calculate QD bands etc. 
   */
  /* if --hmm-W <n> or --hmm-cW was set on command line in combination with --viterbi or --forward, 
   * DO NOT use the QDB band definition alg to calculate W, use <n> from --hmm-W
   */
  hmm_W_set = ( esl_opt_IsOn(go, "--hmm-W") || esl_opt_IsOn(go, "--hmm-cW")) ? TRUE : FALSE;
  if(use_hmmonly && hmm_W_set) { 
    if((status = ConfigCM(cm, errbuf, FALSE, NULL, NULL)) != eslOK) return status;  /* FALSE says: DON'T calculate W */
    if(esl_opt_IsOn(go, "--hmm-W")) cm->W = esl_opt_GetInteger(go, "--hmm-W"); 
    else                            cm->W = (int) (esl_opt_GetReal(go, "--hmm-cW") * cm->clen); 
  }
  else { 
    if((status = ConfigCM(cm, errbuf, TRUE, NULL, NULL)) != eslOK) return status;  /* TRUE says: calculate W */
  }

  /* Setup ScanMatrix for CYK/Inside scanning functions, we can't 
   * do it in initialize_cm(), b/c it's W dependent; W was just set.
   * We don't need it if we're only using an HMM though.
   */
  if(use_hmmonly) cm->smx = NULL;
  else { 
    int do_float = TRUE;
    int do_int   = FALSE;
    if(cm->search_opts & CM_SEARCH_INSIDE) { do_float = FALSE; do_int = TRUE; }
    cm_CreateScanMatrixForCM(cm, do_float, do_int);
    if(cm->smx == NULL) cm_Fail("initialize_cm(), use_hmmonly is FALSE, CreateScanMatrixForCM() call failed, mx is NULL.");
  }

  ESL_DPRINTF1(("cm->pbegin: %.3f\n", cm->pbegin));
  ESL_DPRINTF1(("cm->pend: %.3f\n", cm->pend));

  return eslOK;
}

/* Function: set_searchinfo_for_calibrated_cm()
 * Date:     EPN, Mon Jan 21 08:56:04 2008 (updated)
 * 
 * Purpose:  For a CM with exponential tail and filter thresholds statistics 
 *           from cmcalibrate: determine how many rounds of searching we will 
 *           do (all rounds but last round are filters), and set the relevant 
 *           info in the SearchInfo_t <cm->si> object, including cutoffs.
 *
 *           We only enter this function if the CM file has exponential tail
 *           and filter stats. If we don't we enter 
 *           set_searchinfo_for_uncalibrated_cm().
 * 
 * ------------------------------------------------------------------------------------
 * How HMM filters thresholds are set for calibrated CMs based on command line options:
 * ------------------------------------------------------------------------------------
 * Note that --fil-no-hmm, --fil-T-hmm, --fil-S-hmm, fil-E-hmm are mutually 
 * exclusive (enforced by esl_getopts).
 * 
 * if --fil-no-hmm:     don't filter with an HMM.
 * if --fil-T-hmm <x>:  set HMM filter bit score cutoff as <x>.
 * if --fil-S-hmm <x>:  set HMM filter bit score/E-value cutoff as that which gives
 *                      predicted survival fraction of <x>              
 * if --fil-E-hmm <x>:  set HMM filter E-value cutoff as <x>
 * 
 * If none of these options are selected, default thresholding occurs as follows:
 * The appropriate HMM filter threshold cutoff (that will achieve a predicted
 * sensitivity of F (0.993 by default) as determined by an empirical simulation
 * in cmcalibrate) for the given final CM threshold is chosen from the CM file.
 * 
 * The following options supplement those above and get printed with -h:
 * --fil-Smax-hmm <x>:  if predicted survival fraction (for HMM threshold as set 
 *                      following rules above) is > <x> turn HMM filter off, unless 
 *                      --fil-A-hmm (see below).
 * 
 * --fil-Smin-hmm <x>:  if predicted survival fraction (for HMM threshold as set 
 *                      following rules above) is < <x> reset HMM threshold as that
 *                      which gives predicted survival fraction of <x>.
 * 
 * --fil-A-hmm:         always use an HMM filter. If the predicted survival fraction
 *                      (for HMM threshold as set following rules above) is > <x>
 *                      from --fil-Smax-hmm <x>, reset HMM threshold as that which
 *                      gives predicted survival fraction of <x>.
 *
 * Expert options (only printed with a --devhelp)
 * --fil-Xmin-hmm <x>:  set minimum HMM filter threshold as that which will yield a total
 *                      run time of <x> times the predicted time an HMM only scan would
 *                      take (calc'ed based on number of DP calcs). If the HMM threshold is 
 *                      less than this, then reset it to the threshold that achieves <x>.
 * 
 * --fil-finE-hmm <x>:  for the HMM filter threshold calculation pretend the final E-value 
 *                      calculation is <x>. For ex: with 'cmsearch -E 100 --fil-finE-hmm 1'
 *                      the filter is set as if the final E-value threhsold is 1, but the
 *                      program will report any final hit with E-value < 100.
 *
 * --fil-finT-hmm <x>:  same as --fil-finE-hmm, but <x> is a bit score.
 *                     
 *
 * ---------------------------------------------------------------------------------------
 * How QDB CYK filter thresholds are set for calibrated CMs based on command line options:
 * ---------------------------------------------------------------------------------------
 * Note that --fil-no-qdb, --fil-T-qdb, and --fil-E-qdb are mutually 
 * exclusive (enforced by esl_getopts).
 * 
 * if --fil-no-qdb:     don't filter with QDB CYK.
 * if --fil-T-qdb <x>:  set QDB CYK filter bit score cutoff as <x>.
 * if --fil-E-qdb <x>:  set QDB CYK filter E-value cutoff as <x>
 * 
 * If none of these options are selected, default threhsolding occurs as follows:
 * The QDB threshold E-value is set as that which is 100 times the final algorithm E-value
 * threshold. If this threshold results in the CYK filter requiring less than a 
 * predicted 3\% of the total DP calculations in the full search, then the CYK
 * E-value threshold is raised (made less strict) such that the predicted number
 * of CYK filter DP calcs is exactly 3\% the total number of DP calcs for the full
 * search.
 * 
 * Expert QDB CYK filter options (only printed with a --devhelp)
 * --fil-finE-qdb <x>:  for the QDB filter threshold calculation pretend the final E-value 
 *                      calculation is <x>. For ex: with 'cmsearch -E 100 --fil-finE-qdb 1'
 *                      the filter is set as if the final E-value threhsold is 1, but the
 *                      program will report any final hit with E-value < 100.
 *
 * --fil-finT-qdb <x>:  same as --fil-finE-qdb, but <x> is a bit score.
 * 
 * beta, the tail loss probability for the QDB calculation is set as <x> from --fil-beta <x>.
 *
 *
 * ------------------------------------------------------------------------------
 * How final thresholds are set for calibrated CMs based on command line options:
 * ------------------------------------------------------------------------------
 *
 * How final algorithm is set:
 * --inside:     search with CM inside (TRUE by default)
 * --cyk:        search with CM CYK 
 * --forward:    search with HMM Forward
 * --viterbi:    search with HMM Viterbi
 * 
 * How final algorithm threshold is set:
 * -T <x>:       set final algorithm threshold as bit score <x>
 * -E <x>:       set final algorithm threshold as E-value <x>
 * --ga:         use Rfam gathering threshold (bit sc) from CM file 
 * --tc:         use Rfam trusted cutoff      (bit sc) from CM file
 * --nc:         use Rfam noise cutoff        (bit sc) from CM file
 *
 * NOTE: --ga, --nc, --tc are incompatible with --forward and --viterbi
 *
 */
static int
set_searchinfo_for_calibrated_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int status;                 /* easel status code */
  int n;                      /* counter over rounds */
  int stype;                  /* type of filter */
  int search_opts;            /* search_opts for filter */
  int use_hmmonly;            /* TRUE if --viterbi or --forward */
  int fthr_mode;              /* filter threshold mode */
  float final_E  = -1.;       /* final round E-value cutoff, stays -1 if !CMH_EXPTAIL_STATS */
  float fqdb_E   = -1.;       /* QDB filter round E-value cutoff, stays -1 if !CMH_EXPTAIL_STATS */
  float fhmm_E   = -1.;       /* HMM filter round E-value cutoff, stays -1 if !CMH_EXPTAIL_STATS */
  float final_sc = -1.;       /* final round bit score cutoff */
  float fqdb_sc  = -1.;       /* QDB filter round bit score cutoff */
  float fhmm_sc  = -1.;       /* HMM filter round bit score cutoff */
  int   final_ctype = E_CUTOFF;/* final round cutoff type, SCORE_CUTOFF or E_CUTOFF, init'ed only to silence compiler warnings */
  int   fqdb_ctype = E_CUTOFF;/* QDB filter round cutoff type, SCORE_CUTOFF or E_CUTOFF, init'ed only to silence compiler warnings */
  int   fhmm_ctype = E_CUTOFF;/* HMM filter round cutoff type, SCORE_CUTOFF or E_CUTOFF, init'ed only to silence compiler warnings */
  float final_S = -1;         /* predicted survival fraction from final round */
  float fqdb_S = -1;          /* predicted survival fraction from qdb filter round */
  float fhmm_S = -1;          /* predicted survival fraction from HMM filter round */
  float fhmm_ncalcs_per_res = 0.;  /* number of millions of filter HMM DP calcs predicted per residue */
  float fqdb_ncalcs_per_res = 0.;  /* number of millions of filter QDB DP calcs predicted per residue */
  float final_ncalcs_per_res = 0.; /* number of millions of final stage DP calcs predicted per residue */
  float all_filters_ncalcs_per_res = 0.;/* number of millions of DP calcs predicted for all filter rounds */
  float fqdb_Smin = 1.;       /* minimally useful survival fraction for qdb filter round */
  float fqdb_Emin = 0.;       /* minimally useful E-value cutoff for qdb filter round, calc'ed from fqdb_Smin */
  float fqdb2final_Efactor;   /* fqdb_E is set as max(fqdb_Smin, final_E * fqdb2final_Efactor) */
  float xfil = 0.03;          /* used to set *_Smin values, minimal fraction of filter dp calcs to do in the final round */
  int   do_qdb_filter = TRUE; /* TRUE to add QDB filter, FALSE not to */
  int   do_hmm_filter = TRUE; /* TRUE to add HMM filter, FALSE not to */
  double fqdb_beta_qdb;       /* beta for QDBs in QDB filter round */
  double fqdb_beta_W;         /* beta for W in QDB filter round */
  int    fqdb_W;              /* W for QDB filter round */
  int   *fqdb_dmin, *fqdb_dmax;/* d bands (QDBs) for QDB filter round */
  int safe_windowlen;         /* used to get QDBs */
  ScanMatrix_t *fqdb_smx = NULL;/* the scan matrix for the QDB filter round */
  int cm_mode, hmm_mode;      /* CM exp tail mode and CP9 HMM exp tail mode for E-value statistics */
  int qdb_mode;               /* CM exp tail mode during QDB filter round */
  int cut_point;              /* HMM forward E-value cut point from filter threshold stats */
  float xhmm;                 /* for filtered search, multiplier of QDB ncalcs times predicted HMM dp calcs */
  float final_sc_acc2hmmfil = -1.; /* final round bit score cutoff according to HMM filter, used to calc filter thresholds */
  float final_E_acc2hmmfil  = -1.; /* final round bit score cutoff according to HMM filter, used to calc filter thresholds */
  float final_sc_acc2qdbfil = -1.; /* final round bit score cutoff according to QDB filter, used to calc filter thresholds */
  float final_E_acc2qdbfil  = -1.; /* final round bit score cutoff according to QDB filter, used to calc filter thresholds */

  if(cm->si != NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "set_searchinfo_for_calibrated_cm(), cm->si is not NULL, shouldn't happen.\n");
  if(! (cm->flags & CMH_EXPTAIL_STATS)) ESL_FAIL(eslEINCOMPAT, errbuf, "set_searchinfo_for_calibrated_cm(): but cm: %s has no exp tail stats.", cm->name);
  if(! (cm->flags & CMH_FILTER_STATS))  ESL_FAIL(eslEINCOMPAT, errbuf, "set_searchinfo_for_calibrated_cm(): but cm: %s has no filter stats.", cm->name);

  /* Create SearchInfo, specifying no filtering, we change the threshold below */
  CreateSearchInfo(cm, SCORE_CUTOFF, 0., -1.);
  if(cm->si == NULL) cm_Fail("set_searchinfo_for_calibrated_cm(), CreateSearchInfo() call failed.");
  SearchInfo_t *si = cm->si; 
  if(si->nrounds > 0) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_searchinfo_for_calibrated_cm(), si->nrounds (%d) > 0\n", si->nrounds);
  
  /* First, set up cutoff for final round, this will be round si->nrounds == 0, b/c no filters have been added yet */
  n           = si->nrounds;
  stype       = si->stype[n];
  search_opts = si->search_opts[n];
  use_hmmonly = ((search_opts & CM_SEARCH_HMMVITERBI) || (search_opts & CM_SEARCH_HMMFORWARD));
  if(use_hmmonly) do_hmm_filter = do_qdb_filter = FALSE; /* don't filter if we're searching only with the HMM */
  if((! use_hmmonly) && (stype != SEARCH_WITH_CM)) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_searchinfo_for_calibrated_cm(), final round search_opts HMMVITERBI or HMMFORWARD flags down, but type != SEARCH_WITH_CM.");
  /* determine configuration of CM and CP9 HMM based on cm->flags & cm->search_opts */
  CM2ExpMode(cm, search_opts, &cm_mode, &hmm_mode); 

  final_S = final_E = final_sc = -1.;
  if(!use_hmmonly) { if((status = cm_CountSearchDPCalcs(cm, errbuf, 10*cm->smx->W, cm->smx->dmin, cm->smx->dmax, cm->smx->W, TRUE,  NULL, &(final_ncalcs_per_res))) != eslOK) return status; }
  /* set up final round cutoff, either 0 or 1 of 5 options is enabled. 
   * (note, -E and -T use IsUsed() because they are by default 'on') */
  if( (! esl_opt_IsUsed(go, "-E")) && 
      (! esl_opt_IsUsed(go, "-T"))   && 
      (! esl_opt_IsOn(go, "--ga")) && 
      (! esl_opt_IsOn(go, "--tc")) && 
      (! esl_opt_IsOn(go, "--nc"))) { 
    /* No relevant options enabled, cutoff is default E value cutoff */
    final_ctype = E_CUTOFF;
    final_E     = esl_opt_GetReal(go, "-E");
    if((status  = E2MinScore(cm, errbuf, (use_hmmonly ? hmm_mode : cm_mode), final_E, &final_sc)) != eslOK) return status;
  }
  else if( esl_opt_IsUsed(go, "-E")) { /* -E enabled, use that */
    final_ctype = E_CUTOFF;
    final_E     = esl_opt_GetReal(go, "-E");
    if((status = E2MinScore(cm, errbuf, (use_hmmonly ? hmm_mode : cm_mode), final_E, &final_sc)) != eslOK) return status;
  }
  else if ( esl_opt_IsUsed(go, "-T")) { /* -T enabled, use that */
    final_ctype = SCORE_CUTOFF;
    final_sc    = esl_opt_GetReal(go, "-T");
  }
  else if ( esl_opt_IsOn(go, "--ga")) { /* --ga enabled, use that, if available, else die */
    if(use_hmmonly) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "--ga is incompatible with --viterbi and --forward.");
    if(! (cm->flags & CMH_GA)) ESL_FAIL(eslEINVAL, errbuf, "No GA gathering threshold in CM file, can't use --ga.");
    final_ctype = SCORE_CUTOFF;
    final_sc    = cm->ga;
  }
  else if ( esl_opt_IsOn(go, "--tc")) { /* --tc enabled, use that, if available, else die */
    if(use_hmmonly) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "--tc is incompatible with --viterbi and --forward.");
    if(! (cm->flags & CMH_TC)) ESL_FAIL(eslEINVAL, errbuf, "No TC trusted cutoff in CM file, can't use --tc.");
    final_ctype = SCORE_CUTOFF;
    final_sc    = cm->tc;
  }
  else if ( esl_opt_IsOn(go, "--nc")) { /* --nc enabled, use that, if available, else die */
    if(use_hmmonly) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "--nc is incompatible with --viterbi and --forward.");
    if(! (cm->flags & CMH_NC)) ESL_FAIL(eslEINVAL, errbuf, "No NC noise cutoff in CM file, can't use --nc.");
    final_ctype = SCORE_CUTOFF;
    final_sc    = cm->nc;
  }
  else ESL_FAIL(eslEINCONCEIVABLE, errbuf, "No final round cutoff selected. This shouldn't happen.");

  /* Determine the E-value for bit sc cutoff, or bit sc cutoff for E-value regardless of user options,
   * we'll print these to stdout eventually. Note, we may be repeating calculations here... */
  if(final_ctype == SCORE_CUTOFF) { /* determine max E-value that corresponds to final_sc bit sc cutoff  across all partitions */
    if((status = Score2MaxE(cm, errbuf, (use_hmmonly ? hmm_mode : cm_mode), final_sc, &final_E)) != eslOK) return status;
  }
  else if(final_ctype == E_CUTOFF) { /* determine min bit sc that corresponds to final_E E-val cutoff across all partitions */
    if((status  = E2MinScore(cm, errbuf, (use_hmmonly ? hmm_mode : cm_mode), final_E, &final_sc)) != eslOK) return status;
  }
  final_S = E2SurvFract(final_E, cm->W, cfg->avg_hit_len, cfg->dbsize, FALSE); /* FALSE says don't add a W pad, we're not filtering in final round */
  /* update the search info, which holds the thresholds for final round */
  UpdateSearchInfoCutoff(cm, cm->si->nrounds, final_ctype, final_sc, final_E);   
  ValidateSearchInfo(cm, cm->si);

  /* Handle case where --fil-finE-hmm <x> or --fil-finT-hmm <x> was set on command line, if this is the case, 
   * we pretend <x> is final E/final T cutoff for purposes of setting HMM filter threshold.
   */
  if(esl_opt_IsOn(go, "--fil-finE-hmm")) { 
    final_E_acc2hmmfil = esl_opt_GetReal(go, "--fil-finE-hmm");
    if((status  = E2MinScore(cm, errbuf, (use_hmmonly ? hmm_mode : cm_mode), final_E_acc2hmmfil, &final_sc_acc2hmmfil)) != eslOK) return status;
  }
  else if (esl_opt_IsOn(go, "--fil-finT-hmm")) { 
    final_sc_acc2hmmfil = esl_opt_GetReal(go, "--fil-finT-hmm");
    if((status = Score2MaxE(cm, errbuf, (use_hmmonly ? hmm_mode : cm_mode), final_sc_acc2hmmfil, &final_E_acc2hmmfil)) != eslOK) return status;
  }
  else { /* neither --fil-finE-hmm nor --fil-finT-hmm were set on command line */
    final_sc_acc2hmmfil = final_sc;
    final_E_acc2hmmfil  = final_E;
  }
  /* Handle case where --fil-finE-qdb <x> or --fil-finT-qdb <x> was set on command line, if this is the case, 
   * we pretend <x> is final E/final T cutoff for purposes of setting QDB filter threshold.
   */
  if(esl_opt_IsOn(go, "--fil-finE-qdb")) { 
    final_E_acc2qdbfil = esl_opt_GetReal(go, "--fil-finE-qdb");
    if((status  = E2MinScore(cm, errbuf, (use_hmmonly ? hmm_mode : cm_mode), final_E_acc2qdbfil, &final_sc_acc2qdbfil)) != eslOK) return status;
  }
  else if (esl_opt_IsOn(go, "--fil-finT-qdb")) { 
    final_sc_acc2qdbfil = esl_opt_GetReal(go, "--fil-finT-qdb");
    if((status = Score2MaxE(cm, errbuf, (use_hmmonly ? hmm_mode : cm_mode), final_sc_acc2qdbfil, &final_E_acc2qdbfil)) != eslOK) return status;
  }
  else { /* neither --fil-finE-qdb nor --fil-finT-qdb were set on command line */
    final_sc_acc2qdbfil = final_sc;
    final_E_acc2qdbfil  = final_E;
  }

  /* DumpSearchInfo(cm->si); */
  /* done with threshold for final round */

  /* Set up the filters and their thresholds 
   * A. determine thresholds/stats for HMM filter to add
   * B. determine thresholds/stats for CM QDB filter to add
   * C. add QDB filter, if necessary (before HMM filter, filters added like a stack, HMM filter is added last but used first)
   * D. add HMM filter, if necessary (after QDB filter)
   */

  /* 1. determine thresholds/stats for HMM filter */
  do_hmm_filter = esl_opt_GetBoolean(go, "--fil-hmm");
  fhmm_S = fhmm_E = fhmm_sc = -1.;
  if((status = cp9_GetNCalcsPerResidue(cm->cp9, errbuf, &fhmm_ncalcs_per_res)) != eslOK) return status;
  /* reset HMM mode to either EXP_CP9_LF or EXP_CP9_GF for local or glocal HMM forward filtering */
  hmm_mode = ExpModeIsLocal(cm_mode) ? EXP_CP9_LF : EXP_CP9_GF;

  if(do_hmm_filter) { /* determine thresholds for HMM forward filter */
    if(use_hmmonly) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_searchinfo_for_calibrated_cm(), --fil-hmm enabled, along with --viterbi or --forward, shouldn't happen.");
    if( (! esl_opt_IsUsed(go, "--fil-E-hmm")) && (! esl_opt_IsUsed(go, "--fil-T-hmm")) && (! esl_opt_IsUsed(go, "--fil-S-hmm"))) { /* default: use HMM filter threshold stats, if they exist in cmfile, else use default bit score cutoff */
      /* No relevant options selected. Set HMM filter cutoff as appropriate HMM E value cutoff from cmfile */
      /* determine filter threshold mode, the mode of final stage of searching, either FTHR_CM_LC,
       * FTHR_CM_LI, FTHR_CM_GC, FTHR_CM_GI (can't be an HMM mode b/c --viterbi and --forward toggle --fil-hmm off)
       */

      /* determine fthr mode */
      if((status = CM2FthrMode(cm, errbuf, cm->search_opts, &fthr_mode)) != eslOK) return status;
      HMMFilterInfo_t *hfi_ptr = cm->stats->hfiA[fthr_mode]; /* for convenience */
      if(hfi_ptr->is_valid == FALSE) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_searchinfo_for_calibrated_cm(), cm's CMH_FILTER_STATS is raised but best filter info for fthr_mode %d is invalid.", fthr_mode);
      fhmm_ctype = E_CUTOFF;
      /* determine the appropriate filter cut point <cut_point> to use, for details
       * see comments in function: searchinfo.c:GetHMMFilterFwdECutGivenCME() for details */
      if(final_ctype == SCORE_CUTOFF) { /* final round has bit score cutoff */
	if((status = GetHMMFilterFwdECutGivenCMBitScore(hfi_ptr, errbuf, final_sc_acc2hmmfil, cfg->dbsize, &cut_point, cm, cm_mode)) != eslOK) return status; 
      }
      else { /* final round has E-value cutoff */
	if((status = GetHMMFilterFwdECutGivenCME(hfi_ptr, errbuf, final_E_acc2hmmfil, cfg->dbsize, &cut_point)) != eslOK) return status; 
      }
      if(cut_point != -1 ) { /* it's worth it to filter according to the CM file */
	fhmm_E = hfi_ptr->fwd_E_cut[cut_point] * ((double) cfg->dbsize / (double) hfi_ptr->dbsize); 
	fhmm_S = E2SurvFract(fhmm_E, cm->W, cfg->avg_hit_len, cfg->dbsize, TRUE);
	/* check if --fil-Xmin-hmm applies */
	if(esl_opt_IsOn(go, "--fil-Xmin-hmm")) { 
	  xhmm = GetHMMFilterXHMM(hfi_ptr, cut_point, cm->W, cfg->avg_hit_len, final_ncalcs_per_res, fhmm_ncalcs_per_res);
	  while((xhmm < esl_opt_GetReal(go, "--fil-Xmin-hmm")) && (cut_point >= 0)) { 
	    cut_point--;
	    xhmm = GetHMMFilterXHMM(hfi_ptr, cut_point, cm->W, cfg->avg_hit_len, final_ncalcs_per_res, fhmm_ncalcs_per_res);
	    fhmm_E = hfi_ptr->fwd_E_cut[cut_point] * ((double) cfg->dbsize / (double) hfi_ptr->dbsize); 
	    fhmm_S = E2SurvFract(fhmm_E, cm->W, cfg->avg_hit_len, cfg->dbsize, TRUE);
	  }
	}
	/* check if --fil-Smax-hmm applies */
	if(fhmm_S > esl_opt_GetReal(go, "--fil-Smax-hmm")) { /* predicted survival fraction exceeds maximum allowed, turn filter off, unless --fil-A-hmm enabled */
	  if(esl_opt_GetBoolean(go, "--fil-A-hmm")) { /* we always filter, set threshold as that which gives predicted surv fract of <x> from --fil-Smax-hmm <x> */
	    fhmm_E = SurvFract2E(esl_opt_GetReal(go, "--fil-Smax-hmm"), cm->W, cfg->avg_hit_len, cfg->dbsize);
	  }
	  else { /* --fil-A-hmm not enabled */
	    /* it's not worth it to filter, our HMM filter cutoff would be so low, 
	     * letting so much of the db survive, the filter is a waste of time */
	    do_hmm_filter = FALSE;
	    ESL_DPRINTF1(("cut_point not -1, fhmm_S exceeds max allowed by --fil_Smax-hmm\n"));
	  }
	}
	/* check if --fil-Smin-hmm applies */
	if((esl_opt_IsOn(go, "--fil-Smin-hmm")) && (fhmm_S < esl_opt_GetReal(go, "--fil-Smin-hmm"))) { /* predicted survival fraction is below minimum allowed, set E cutoff as E value that gives min allowed survival fraction */
	  fhmm_E = SurvFract2E(esl_opt_GetReal(go, "--fil-Smin-hmm"), cm->W, cfg->avg_hit_len, cfg->dbsize);
	}
      }
      else { /* cut_point is -1, CM file told us it's not worth it to filter, we only filter if --fil-A-hmm is enabled */
	if(esl_opt_GetBoolean(go, "--fil-A-hmm")) { 
	  fhmm_E = SurvFract2E(esl_opt_GetReal(go, "--fil-Smax-hmm"), cm->W, cfg->avg_hit_len, cfg->dbsize);
	}
	else { /* --fil-A-hmm not enabled and CM file told us it's not worth it to filter */
	  do_hmm_filter = FALSE;
	  /* it's not worth it to filter, our HMM filter cutoff would be so low, 
	   * letting so much of the db survive, the filter is a waste of time */
	  ESL_DPRINTF1(("cut_point -1, always_better FALSE\n"));
	}
      }
    } /* end of if( !esl_opt_IsUsed(go, "--fil-E-hmm") && !esl_opt_IsUsed(go, "--fil-T-hmm") && !esl_opt_IsUsed(go, "--fil-S-hmm")) */
    else if( esl_opt_IsUsed(go, "--fil-E-hmm")) {
      fhmm_ctype = E_CUTOFF;
      fhmm_E     = esl_opt_GetReal(go, "--fil-E-hmm");
      /*fhmm_E     = ESL_MIN(fhmm_E, 1.0);*/ /* we don't allow filter E cutoffs below 1. */
    }
    else if( esl_opt_IsUsed(go, "--fil-T-hmm")) {
      fhmm_ctype = SCORE_CUTOFF;
      fhmm_sc    = esl_opt_GetReal(go, "--fil-T-hmm");
    }
    else if(esl_opt_IsUsed(go, "--fil-S-hmm")) {
      fhmm_ctype = E_CUTOFF;
      fhmm_E = SurvFract2E(esl_opt_GetReal(go, "--fil-S-hmm"), cm->W, cfg->avg_hit_len, cfg->dbsize);
    }
    else ESL_FAIL(eslEINCONCEIVABLE, errbuf, "No HMM filter cutoff selected. This shouldn't happen.");

    if(do_hmm_filter) { 
      if(fhmm_ctype == SCORE_CUTOFF) { /* determine max E-value that corresponds to fhmm_sc bit sc cutoff across all partitions */
	if((status = Score2MaxE(cm, errbuf, hmm_mode, fhmm_sc, &fhmm_E)) != eslOK) return status;
      }
      else if(fhmm_ctype == E_CUTOFF) { /* determine min bit sc that corresponds to fhmm_E E-val cutoff across all partitions */
	if((status  = E2MinScore(cm, errbuf, hmm_mode, fhmm_E, &fhmm_sc)) != eslOK) return status;
      }
      fhmm_S = E2SurvFract(fhmm_E, cm->W, cfg->avg_hit_len, cfg->dbsize, TRUE);
    }
  }

  /* B. determine thresholds/stats for CM QDB filter to add (do this after HMM stats b/c we use fhmm_S when setting fqdb_E */
  qdb_mode = (ExpModeIsLocal(cm_mode)) ? EXP_CM_LC : EXP_CM_GC; /* always do CYK with QDB filter, only question is local or glocal? */
  fqdb_S = fqdb_E = fqdb_sc = -1.;

  if(do_qdb_filter && esl_opt_GetBoolean(go, "--fil-qdb")) { /* determine thresholds, beta for qdb cyk filter */
    /* build the ScanMatrix_t for the QDB filter round, requires calcing dmin, dmax */

    fqdb_beta_qdb = esl_opt_GetReal(go, "--fil-beta");

    safe_windowlen = cm->W * 3;
    while(!(BandCalculationEngine(cm, safe_windowlen, fqdb_beta_qdb, FALSE, &fqdb_dmin, &fqdb_dmax, NULL, NULL))) {
      free(fqdb_dmin);
      free(fqdb_dmax);
      fqdb_dmin = NULL;
      fqdb_dmax = NULL;
      safe_windowlen *= 2;
      if(safe_windowlen > (cm->clen * 1000)) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_searchinfo_for_calibrated_cm(), band calculation safe_windowlen big: %d\n", safe_windowlen);
    }
    /* tricky step here, W is calculated for filter using maximum of fqdb_beta_qdb and cm->beta_W, this is b/c 
     * model was (possibly) calibrated with W as calc'ed with cm->beta_W and we don't want a bigger hit
     * to be possible then was during calibration (could overestimate significance of scores), 
     * W can be less than cm->beta_W though because that would only lead to possible underestimation 
     * of significance of scores (the good direction).
     */
    if(cm->beta_W > fqdb_beta_qdb) { 
      fqdb_beta_W = cm->beta_W;
      fqdb_W      = cm->W;
    }
    else { 
      fqdb_beta_W = fqdb_beta_qdb;
      fqdb_W      = fqdb_dmax[0];
    }
    fqdb_smx = cm_CreateScanMatrix(cm, fqdb_W, fqdb_dmin, fqdb_dmax, fqdb_beta_W, fqdb_beta_qdb, TRUE, TRUE, FALSE);
    if((status = cm_CountSearchDPCalcs(cm, errbuf, 10*fqdb_smx->W, fqdb_smx->dmin, fqdb_smx->dmax, fqdb_smx->W, TRUE,  NULL, &(fqdb_ncalcs_per_res))) != eslOK) return status;

    if( (! esl_opt_IsOn(go, "--fil-E-qdb")) && (! esl_opt_IsUsed(go, "--fil-T-qdb"))) {
      /* No relevant options selected, use default method for setting CM QDB filter */
      fqdb_ctype = E_CUTOFF;

      /* Determine the 'minimally useful' QDB CYK filter round survival fraction. 
      * This is the survival fraction <fqdb_Smin> that we predict would require the 
      * *final* round to perform <xfil> * x millions of DP calcs, where x is the number of
      * predicted DP calcs predicted for all filter rounds. This calculation assumes the 
      * final round will have to search fqdb_Smin fraction of the database.
      */
      if(do_hmm_filter) all_filters_ncalcs_per_res = fhmm_ncalcs_per_res + (fqdb_ncalcs_per_res * fhmm_S);
      else              all_filters_ncalcs_per_res = fqdb_ncalcs_per_res;
      fqdb_Smin = xfil * (all_filters_ncalcs_per_res / final_ncalcs_per_res); 
      fqdb_Emin = SurvFract2E(fqdb_Smin, fqdb_smx->W, cfg->avg_hit_len, cfg->dbsize);
      
      fqdb2final_Efactor = 100.;
      fqdb_E = ESL_MAX(final_E_acc2qdbfil * fqdb2final_Efactor, fqdb_Emin);
    }
    else if (esl_opt_IsUsed(go, "--fil-E-qdb")) { /* survival fraction for QDB filter set on command line, use that */
      fqdb_ctype = E_CUTOFF;
      fqdb_E     = esl_opt_GetReal(go, "--fil-E-qdb");
    }
    else if (esl_opt_IsUsed(go, "--fil-T-qdb")) {
      fqdb_ctype = SCORE_CUTOFF;
      fqdb_sc    = esl_opt_GetReal(go, "--fil-T-qdb");
    }
    else ESL_FAIL(eslEINCONCEIVABLE, errbuf, "No CM filter cutoff selected. This shouldn't happen.");

    if(fqdb_ctype == SCORE_CUTOFF) { /* determine max E-value that corresponds to fqdb_sc bit sc cutoff  across all partitions */
      if((status = Score2MaxE(cm, errbuf, qdb_mode, fqdb_sc, &fqdb_E)) != eslOK) return status;
    }
    else if(fqdb_ctype == E_CUTOFF) { /* determine min bit sc that corresponds to fqdb_E E-val cutoff across all partitions */
      if((status  = E2MinScore(cm, errbuf, qdb_mode, fqdb_E, &fqdb_sc)) != eslOK) return status;
    }

    fqdb_S = E2SurvFract(fqdb_E, cm->W, cfg->avg_hit_len, cfg->dbsize, TRUE);
    /* if HMM is predicted to be a better filter, turn QDB filter OFF*/
    if(do_hmm_filter && (fqdb_S > fhmm_S)) do_qdb_filter = FALSE;
    /* if fqdb_S > 0.9, turn QDB filter off, it's not worth it */
    if(fqdb_S > 0.9) do_qdb_filter = FALSE;

  }
  else do_qdb_filter = FALSE;

  /* C. add QDB filter, if necessary (before HMM filter, filters added like a stack, HMM filter is added last but used first) */
  if(do_qdb_filter) { 
    AddFilterToSearchInfo(cm, TRUE, FALSE, FALSE, FALSE, FALSE, fqdb_smx, NULL, fqdb_ctype, fqdb_sc, fqdb_E, (! esl_opt_GetBoolean(go, "--no-null3")));
    /* DumpSearchInfo(cm->si); */
  }
  else if (fqdb_smx != NULL) cm_FreeScanMatrix(cm, fqdb_smx); 
  /* D. add HMM filter, if necessary (after QDB filter, filters added like a stack, HMM filter is added last but used first) */
  if (do_hmm_filter) { 
    AddFilterToSearchInfo(cm, FALSE, FALSE, FALSE, TRUE, FALSE, NULL, NULL, fhmm_ctype, fhmm_sc, fhmm_E, (! esl_opt_GetBoolean(go, "--no-null3")));
    /*DumpSearchInfo(cm->si); */
  }
  ValidateSearchInfo(cm, cm->si);
  /* DumpSearchInfo(cm->si); */
  return eslOK;
}

/* Function: set_searchinfo_for_uncalibrated_cm()
 * Date:     EPN, Thu Mar  6 05:29:16 2008
 * 
 * Purpose:  For a CM WITHOUT exponential tail and filter thresholds statistics 
 *           from cmcalibrate: determine how many rounds of searching we will 
 *           do (all rounds but last round are filters), and set the relevant 
 *           info in the SearchInfo_t <cm->si> object, including cutoffs.
 *
 * 
 * --------------------------------------------------------------------------------------
 * How HMM filters thresholds are set for UNcalibrated CMs based on command line options:
 * --------------------------------------------------------------------------------------
 * Note that --fil-no-hmm, --fil-T-hmm, --fil-S-hmm, fil-E-hmm are mutually 
 * exclusive (enforced by esl_getopts).
 * 
 * if --fil-no-hmm:     don't filter with an HMM.
 * if --fil-T-hmm <x>:  set HMM filter bit score cutoff as <x>.
 * if --fil-S-hmm <x>:  ERROR, requires calibration; esl_getopts detects and exits
 * if --fil-E-hmm <x>:  ERROR, requires calibration; esl_getopts detects and exits
 * 
 * If none of these options are selected, default thresholding occurs as follows:
 * HMM filter threshold is set as 3.0 bits.
 * 
 * -----------------------------------------------------------------------------------------
 * How QDB CYK filter thresholds are set for UNcalibrated CMs based on command line options:
 * -----------------------------------------------------------------------------------------
 * Note that --fil-no-qdb, --fil-T-qdb, and --fil-E-qdb are mutually 
 * exclusive (enforced by esl_getopts).
 * 
 * if --fil-no-qdb:     don't filter with QDB CYK.
 * if --fil-T-qdb <x>:  set QDB CYK filter bit score cutoff as <x>.
 * if --fil-E-qdb <x>:  ERROR, requires calibration; esl_getopts detects and exits
 * 
 * If none of these options are selected, default threhsolding occurs as follows:
 * QDB filter threshold is set as 0.0 bits.
 * 
 * Expert QDB CYK filter options (only printed with a --devhelp)
 * --fil-finE-qdb <x>:  ERROR, requires calibration; esl_getopts detects and exits
 * --fil-finT-qdb <x>:  ERROR, requires calibration; esl_getopts detects and exits
 * 
 * beta, the tail loss probability for the QDB calculation is set as <x> from --fil-beta <x>.
 *
 * --------------------------------------------------------------------------------
 * How final thresholds are set for UNcalibrated CMs based on command line options:
 * --------------------------------------------------------------------------------
 *
 * How final algorithm is set:
 * --inside:     search with CM inside (TRUE by default)
 * --cyk:        search with CM CYK 
 * --forward:    search with HMM Forward
 * --viterbi:    search with HMM Viterbi
 * 
 * How final algorithm threshold is set:
 * -T <x>:       set final algorithm threshold as bit score <x>
 * -E <x>:       ERROR, esl_getopts detects this and exits
 * --ga:         use Rfam gathering threshold (bit sc) from CM file 
 * --tc:         use Rfam trusted cutoff      (bit sc) from CM file
 * --nc:         use Rfam noise cutoff        (bit sc) from CM file
 *
 * NOTE: --ga, --nc, --tc are incompatible with --forward and --viterbi
 *
 */
static int
set_searchinfo_for_uncalibrated_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int n;                      /* counter over rounds */
  int stype;                  /* type of filter */
  int search_opts;            /* search_opts for filter */
  int use_hmmonly;            /* TRUE if --viterbi or --forward */
  float final_sc = -1.;       /* final round bit score cutoff */
  float fqdb_sc  = -1.;       /* QDB filter round bit score cutoff */
  float fhmm_sc  = -1.;       /* HMM filter round bit score cutoff */
  int   do_qdb_filter = TRUE; /* TRUE to add QDB filter, FALSE not to */
  int   do_hmm_filter = TRUE; /* TRUE to add HMM filter, FALSE not to */
  double fqdb_beta_qdb;       /* beta for QDBs in QDB filter round */
  double fqdb_beta_W;         /* beta for W in QDB filter round */
  int    fqdb_W;              /* W for QDB filter round */
  int   *fqdb_dmin, *fqdb_dmax; /* d bands (QDBs) for QDB filter round */
  int safe_windowlen;         /* used to get QDBs */
  ScanMatrix_t *fqdb_smx = NULL;/* the scan matrix for the QDB filter round */
  int cm_mode, hmm_mode;      /* CM exp tail mode and CP9 HMM exp tail mode for E-value statistics */
  int qdb_mode;               /* CM exp tail mode during QDB filter round */

  if(cm->si != NULL)                ESL_FAIL(eslEINCOMPAT, errbuf, "set_searchinfo_for_uncalibrated_cm(), cm->si is not NULL, shouldn't happen.\n");
  if(cm->flags & CMH_EXPTAIL_STATS) ESL_FAIL(eslEINCOMPAT, errbuf, "set_searchinfo_for_uncalibrated_cm(): but cm: %s has exp tail stats.", cm->name);
  if(cm->flags & CMH_FILTER_STATS)  ESL_FAIL(eslEINCOMPAT, errbuf, "set_searchinfo_for_uncalibrated_cm(): but cm: %s has filter stats.", cm->name);

  if(esl_opt_IsUsed(go, "--fil-finE-hmm")) { /* can't deal with this b/c we don't have exp tail stats */
    ESL_FAIL(eslEINVAL, errbuf, "--fil-finE-hmm requires exp tail statistics in <cm file>. Use cmcalibrate to get exp tail stats.");
  }
  if(esl_opt_IsUsed(go, "--fil-finT-hmm")) { /* can't deal with this b/c we don't have exp tail stats */
    ESL_FAIL(eslEINVAL, errbuf, "--fil-finT-hmm requires exp tail statistics in <cm file>. Use cmcalibrate to get exp tail stats.");
  }
  if(esl_opt_IsUsed(go, "--fil-finE-qdb")) { /* can't deal with this b/c we don't have exp tail stats */
    ESL_FAIL(eslEINVAL, errbuf, "--fil-finE-qdb requires exp tail statistics in <cm file>. Use cmcalibrate to get exp tail stats.");
  }
  if(esl_opt_IsUsed(go, "--fil-finT-qdb")) { /* can't deal with this b/c we don't have exp tail stats */
    ESL_FAIL(eslEINVAL, errbuf, "--fil-finT-qdb requires exp tail statistics in <cm file>. Use cmcalibrate to get exp tail stats.");
  }

  /* Create SearchInfo, specifying no filtering, we change the threshold below */
  CreateSearchInfo(cm, SCORE_CUTOFF, 0., -1.);
  if(cm->si == NULL) cm_Fail("set_searchinfo_for_uncalibrated_cm(), CreateSearchInfo() call failed.");
  SearchInfo_t *si = cm->si; 
  if(si->nrounds > 0) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_searchinfo_for_uncalibrated_cm(), si->nrounds (%d) > 0\n", si->nrounds);
  
  /* First, set up cutoff for final round, this will be round si->nrounds == 0, b/c no filters have been added yet */
  n           = si->nrounds;
  stype       = si->stype[n];
  search_opts = si->search_opts[n];
  use_hmmonly = ((search_opts & CM_SEARCH_HMMVITERBI) || (search_opts & CM_SEARCH_HMMFORWARD));
  if(use_hmmonly) do_hmm_filter = do_qdb_filter = FALSE; /* don't filter if we're searching only with the HMM */
  if((! use_hmmonly) && (stype != SEARCH_WITH_CM)) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_searchinfo_for_uncalibrated_cm(), final round search_opts HMMVITERBI or HMMFORWARD flags down, but type != SEARCH_WITH_CM.");

  /* determine configuration of CM and CP9 HMM based on cm->flags & cm->search_opts */
  CM2ExpMode(cm, search_opts, &cm_mode, &hmm_mode); 

  /* set up final round cutoff, either 0 or 1 of 5 options is enabled. 
   * (note, -E and -T use IsUsed() because they are by default 'on') */
  if( (! esl_opt_IsUsed(go, "-E"))   && 
      (! esl_opt_IsUsed(go, "-T"))   && 
      (! esl_opt_IsOn(go, "--ga")) && 
      (! esl_opt_IsOn(go, "--tc")) && 
      (! esl_opt_IsOn(go, "--nc"))) { 
    /* No relevant options enabled, cutoff is default bit score cutoff */
    final_sc    = esl_opt_GetReal(go, "-T");
  }
  else if( esl_opt_IsUsed(go, "-E")) { /* -E enabled, error b/c we don't have exp tail stats */
    ESL_FAIL(eslEINVAL, errbuf, "-E requires exp tail statistics in <cm file>. Use cmcalibrate to get exp tail stats.");
  }
  else if( esl_opt_IsUsed(go, "-T")) { /* -T enabled, use that */
    final_sc    = esl_opt_GetReal(go, "-T");
  }
  else if( esl_opt_IsOn(go, "--ga")) { /* --ga enabled, use that, if available, else die */
    if(use_hmmonly) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "--ga is incompatible with --viterbi and --forward.");
    if(! (cm->flags & CMH_GA)) ESL_FAIL(eslEINVAL, errbuf, "No GA gathering threshold in CM file, can't use --ga.");
    final_sc    = cm->ga;
  }
  else if( esl_opt_IsOn(go, "--tc")) { /* --tc enabled, use that, if available, else die */
    if(use_hmmonly) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "--tc is incompatible with --viterbi and --forward.");
    if(! (cm->flags & CMH_TC)) ESL_FAIL(eslEINVAL, errbuf, "No TC trusted cutoff in CM file, can't use --tc.");
    final_sc    = cm->tc;
  }
  else if( esl_opt_IsOn(go, "--nc")) { /* --nc enabled, use that, if available, else die */
    if(use_hmmonly) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "--nc is incompatible with --viterbi and --forward.");
    if(! (cm->flags & CMH_NC)) ESL_FAIL(eslEINVAL, errbuf, "No NC noise cutoff in CM file, can't use --nc.");
    final_sc    = cm->nc;
  }
  else ESL_FAIL(eslEINCONCEIVABLE, errbuf, "No final round cutoff selected. This shouldn't happen.");
  /* update the search info, which holds the thresholds for final round */
  UpdateSearchInfoCutoff(cm, cm->si->nrounds, SCORE_CUTOFF, final_sc, -1.);   
  ValidateSearchInfo(cm, cm->si);
  /* DumpSearchInfo(cm->si); */

  /* done with threshold for final round */

  /* Set up the filters and their thresholds 
   * A. determine thresholds/stats for HMM filter to add
   * B. determine thresholds/stats for CM QDB filter to add
   * C. add QDB filter, if necessary (before HMM filter, filters added like a stack, HMM filter is added last but used first)
   * D. add HMM filter, if necessary (after QDB filter)
   */

  /* 1. determine thresholds/stats for HMM filter */
  do_hmm_filter = esl_opt_GetBoolean(go, "--fil-hmm");
  fhmm_sc = -1.;
  if(do_hmm_filter) { /* determine thresholds for HMM forward filter */
    if(use_hmmonly) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_searchinfo_for_uncalibrated_cm(), --fil-hmm enabled, along with --viterbi or --forward, shouldn't happen.");
    if( (! esl_opt_IsUsed(go, "--fil-E-hmm")) && (! esl_opt_IsUsed(go, "--fil-T-hmm"))) { 
      /* No relevant options selected, set cutoff as default HMM filter bit score. */
      fhmm_sc    = esl_opt_GetReal(go, "--fil-T-hmm");
    } /* end of if( !esl_opt_IsUsed(go, "--fil-E-hmm") && !esl_opt_IsUsed(go, "--fil-T-hmm")) */
    else if(esl_opt_IsUsed(go, "--fil-E-hmm")) { /* can't deal with this b/c we don't have exp tail stats */
      ESL_FAIL(eslEINVAL, errbuf, "--fil-E-hmm requires exp tail statistics in <cm file>. Use cmcalibrate to get exp tail stats.");
    }
    else if(esl_opt_IsUsed(go, "--fil-S-hmm")) { /* can't deal with this b/c we don't have exp tail stats */
      ESL_FAIL(eslEINVAL, errbuf, "--fil-S-hmm requires exp tail statistics in <cm file>. Use cmcalibrate to get exp tail stats.");
    }
    else if(esl_opt_IsUsed(go, "--fil-T-hmm")) { /* HMM filter bit score threshold was set on command line, use that */
      fhmm_sc = esl_opt_GetReal(go, "--fil-T-hmm");
    }
    else ESL_FAIL(eslEINCONCEIVABLE, errbuf, "No HMM filter cutoff selected. This shouldn't happen.");
  }
  if(esl_opt_IsUsed(go, "--fil-Smin-hmm")) { /* can't deal with this b/c we don't have exp tail stats */
      ESL_FAIL(eslEINVAL, errbuf, "--fil-Smin-hmm requires exp tail statistics in <cm file>. Use cmcalibrate to get exp tail stats.");
  }
  if(esl_opt_IsUsed(go, "--fil-Smax-hmm")) { /* can't deal with this b/c we don't have exp tail stats */
      ESL_FAIL(eslEINVAL, errbuf, "--fil-Smax-hmm requires exp tail statistics in <cm file>. Use cmcalibrate to get exp tail stats.");
  }

  /* B. determine thresholds/stats for CM QDB filter to add */
  qdb_mode = (ExpModeIsLocal(cm_mode)) ? EXP_CM_LC : EXP_CM_GC; /* always do CYK with QDB filter, only question is local or glocal? */
  fqdb_sc = -1.;
  if(do_qdb_filter && esl_opt_GetBoolean(go, "--fil-qdb")) { /* determine thresholds, beta for qdb cyk filter */
    /* build the ScanMatrix_t for the QDB filter round, requires calcing dmin, dmax */
    fqdb_beta_qdb = esl_opt_GetReal(go, "--fil-beta");

    safe_windowlen = cm->W * 3;
    while(!(BandCalculationEngine(cm, safe_windowlen, fqdb_beta_qdb, FALSE, &fqdb_dmin, &fqdb_dmax, NULL, NULL))) {
      free(fqdb_dmin);
      free(fqdb_dmax);
      fqdb_dmin = NULL;
      fqdb_dmax = NULL;
      safe_windowlen *= 2;
      if(safe_windowlen > (cm->clen * 1000)) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_searchinfo_for_uncalibrated_cm(), band calculation safe_windowlen big: %d\n", safe_windowlen);
    }
    /* tricky step here, W is calculated for filter using maximum of fqdb_beta_qdb and cm->beta_W, this is b/c 
     * model was (possibly) calibrated with W as calc'ed with cm->beta_W and we don't want a bigger hit
     * to be possible then was during calibration (could overestimate significance of scores), 
     * W can be less than cm->beta_W though because that would only lead to possible underestimation 
     * of significance of scores (the good direction).
     */
    if(cm->beta_W > fqdb_beta_qdb) { 
      fqdb_beta_W = cm->beta_W;
      fqdb_W      = cm->W;
    }
    else { 
      fqdb_beta_W = fqdb_beta_qdb;
      fqdb_W      = fqdb_dmax[0];
    }
    fqdb_smx = cm_CreateScanMatrix(cm, fqdb_W, fqdb_dmin, fqdb_dmax, fqdb_beta_W, fqdb_beta_qdb, TRUE, TRUE, FALSE);
    
    if((!esl_opt_IsUsed(go, "--fil-E-qdb")) && (!esl_opt_IsUsed(go, "--fil-T-qdb"))) {
      /* No relevant options selected, set CM QDB filter as default bit score. */
      fqdb_sc    = esl_opt_GetReal(go, "--fil-T-qdb");
    }
    else if( esl_opt_IsUsed(go, "--fil-E-qdb")) { /* can't deal with this b/c we don't have exp tail stats */
      ESL_FAIL(eslEINVAL, errbuf, "--fil-E-qdb requires exp tail statistics in <cm file>. Use cmcalibrate to get exp tail stats.");
    }
    else if( esl_opt_IsUsed(go, "--fil-T-qdb")) { /* CM CYK bit score threshold was set on command line, use that */
      fqdb_sc    = esl_opt_GetReal(go, "--fil-T-qdb");
    }
    else ESL_FAIL(eslEINCONCEIVABLE, errbuf, "No CM filter cutoff selected. This shouldn't happen.");
  }
  else do_qdb_filter = FALSE;

  /* C. add QDB filter, if necessary (before HMM filter, filters added like a stack, HMM filter is added last but used first) */
  if(do_qdb_filter) { 
    AddFilterToSearchInfo(cm, TRUE, FALSE, FALSE, FALSE, FALSE, fqdb_smx, NULL, SCORE_CUTOFF, fqdb_sc, -1., (! esl_opt_GetBoolean(go, "--no-null3")));
    /* DumpSearchInfo(cm->si); */
  }
  else if (fqdb_smx != NULL) cm_FreeScanMatrix(cm, fqdb_smx); 
  /* D. add HMM filter, if necessary (after QDB filter, filters added like a stack, HMM filter is added last but used first) */
  if (do_hmm_filter) { 
    AddFilterToSearchInfo(cm, FALSE, FALSE, FALSE, TRUE, FALSE, NULL, NULL, SCORE_CUTOFF, fhmm_sc, -1., (! esl_opt_GetBoolean(go, "--no-null3")));
    /* DumpSearchInfo(cm->si); */
  }
  ValidateSearchInfo(cm, cm->si);
  return eslOK;
}

/*
 * Function: read_next_search_seq
 *
 * Date:     RJK, Wed May 29, 2002 [St. Louis]
 *           easeled: EPN, Fri Dec  8 11:40:20 2006
 *
 * Purpose:  Given a dbfp and whether or not to take the reverse complement,
 *           reads in the next sequence and prepares reverse complement if nec.
 *
 * Returns:  eslOK on success; eslEOF if end of file, 
 *           some other status code from esl_sqio_Read() if an error occurs.
 */
int read_next_search_seq (const ESL_ALPHABET *abc, ESL_SQFILE *dbfp, int do_revcomp, dbseq_t **ret_dbseq) 
{
  int status;
  dbseq_t *dbseq = NULL;

  ESL_ALLOC(dbseq, sizeof(dbseq_t));
  dbseq->sq[0] = NULL;
  dbseq->sq[1] = NULL;

  dbseq->sq[0] = esl_sq_CreateDigital(abc);

  while((status = esl_sqio_Read(dbfp, dbseq->sq[0])) == eslOK && (dbseq->sq[0]->n == 0)) /* skip zero length seqs */
    esl_sq_Reuse(dbseq->sq[0]);

  if(status != eslOK) goto ERROR;
  if (do_revcomp) {
    /* make a new ESL_SQ object, to store the reverse complement */
    if((dbseq->sq[1] = esl_sq_CreateDigitalFrom(abc, dbseq->sq[0]->name, dbseq->sq[0]->dsq, 
						dbseq->sq[0]->n, dbseq->sq[0]->desc, 
						dbseq->sq[0]->acc, dbseq->sq[0]->ss)) == NULL) goto ERROR;
    /* reverse complement it in place */
    revcomp(dbseq->sq[1]->abc, dbseq->sq[1], dbseq->sq[1]);
  }
  dbseq->results[0] = NULL;
  dbseq->results[1] = NULL;

  *ret_dbseq = dbseq;
  return eslOK;

 ERROR:
  if(dbseq->sq[0] != NULL) esl_sq_Destroy(dbseq->sq[0]);
  if(dbseq->sq[1] != NULL) esl_sq_Destroy(dbseq->sq[1]);
  if(dbseq != NULL) free(dbseq);
  return status;
}



/* Function: print_run_info
 * Date:     EPN, Mon Mar  3 11:09:18 2008
 *
 * Purpose:  Print information on this run of cmsearch.
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

  if((status = get_command(go, errbuf, &command)) != eslOK) return status;
  if((status = GetDate    (errbuf, &date))    != eslOK) return status;

  fprintf(stdout, "%-13s %s\n",  "# command:", command);
  fprintf(stdout, "%-13s %s\n",  "# date:",    date);
  if(cfg->nproc > 1) fprintf(stdout, "%-13s %d\n", "# nproc:", cfg->nproc);
  fprintf(stdout, "%-13s %d\n",  "# num seqs:",   cfg->nseq);
  fprintf(stdout, "%-13s %.6f\n",  "# dbsize(Mb):", (double) cfg->dbsize / 1000000.);

  if(cfg->tfp != NULL) { 
    fprintf(cfg->tfp, "%-13s %s\n",  "# command:", command);
    fprintf(cfg->tfp, "%-13s %s\n",  "# date:",    date);
    if(cfg->nproc > 1) fprintf(cfg->tfp, "%-13s %d\n", "# nproc:", cfg->nproc);
    fprintf(cfg->tfp, "%-13s %d\n",  "# num seqs:",   cfg->nseq);
    fprintf(cfg->tfp, "%-13s %.6f\n",  "# dbsize(Mb):", (double) cfg->dbsize / 1000000.);
  }

  free(command);
  free(date);
  return eslOK;
}

/* Function: get_command
 * Date:     EPN, Mon Mar  3 11:10:36 2008
 *
 * Purpose:  Return the command used to call cmsearch
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


/* Function: print_searchinfo_for_calibrated_cm
 * Date:     EPN, Thu May 17 14:47:36 2007
 * Purpose:  Print info about search (cutoffs, algorithm, etc.). 
 *           with a CM that has exp tail stats. Can be called in 
 *           2 different modes, mode 1 is 'pre-search', called prior
 *           before a search, mode 2 is 'post-search', called after
 *           a search is done. We're in mode 1 iff cm_surv_fractA == NULL
 *           and cm_nhitsA == NULL.
 *
 *            
 */
int print_searchinfo_for_calibrated_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float *cm_surv_fractA, int *cm_nhitsA, double in_asec, double in_total_psec, double *ret_total_psec)
{
  int status;
  int n;
  float surv_fract = 1.;
  float prv_surv_fract = 1.;
  int cutoff_type;
  float sc_cutoff;
  float e_cutoff;
  int using_filters;
  int cm_mode;
  int cp9_mode;
  int exp_mode;             /* index to use for exp tail stats in cm->stats->gumAA[], exp_mode = (stype == SEARCH_WITH_CM) cm_mode : cp9_mode; */
  int stype;
  int search_opts;
  ScanMatrix_t *smx;
  HybridScanInfo_t *hsi;
  int   use_qdb;            /* are we using qdb for current round? */
  char  time_buf[128];	    /* for printing predicted times */
  int   do_pad;             /* TRUE to add W pad onto survival fraction (FALSE only in final round of searching) */
  int   pre_search_mode;    /* TRUE if this function was called before the search was run */
  double sec_per_res;        /* seconds required to search 1 residue with current model, current round of searching */
  double psec;              /* predicted seconds for current cm, current round */
  double total_psec = 0.;   /* predicted number of seconds for full round */
  ESL_RANDOMNESS *r = NULL; 

  pre_search_mode = (cm_surv_fractA == NULL && cm_nhitsA == NULL) ? TRUE : FALSE;
  if(pre_search_mode && ret_total_psec == NULL) ESL_FAIL(eslEINVAL, errbuf, "print_searchinfo_for_calibrated_cm, pre-search mode, but ret_total_psec is NULL.");

  r = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  if(r == NULL) ESL_FAIL(eslEMEM, errbuf, "print_searchinfo_for_calibrated_cm, memory error, couldn't create randomness object.");

  /* Could use ESL_GETOPTS here, but using the CM flags assures we're reporting
   * on how the CM is actually config'ed, not how we want it to be
   */

  if(! (cm->flags & CMH_EXPTAIL_STATS)) ESL_FAIL(eslEINCOMPAT, errbuf, "print_searchinfo_for_calibrated_cm(): cm: %s does not have exp tail stats.", cm->name);
  using_filters = (cm->si->nrounds > 0) ? TRUE : FALSE;

  if(pre_search_mode) { 
    fprintf(stdout, "#\n# Pre-search info for CM %d: %s\n", cfg->ncm, cm->name);
    if(cfg->tfp != NULL) fprintf(cfg->tfp, "#\n# Pre-search info for CM %d: %s\n", cfg->ncm, cm->name);
  }
  else { 
    fprintf(stdout, "#\n# Post-search info for CM %d: %s\n", cfg->ncm, cm->name);
    if(cfg->tfp != NULL) fprintf(cfg->tfp, "#\n# Post-search info for CM %d: %s\n", cfg->ncm, cm->name);
  }
  fprintf(stdout, "#\n");
  if(cfg->tfp != NULL) fprintf(cfg->tfp, "#\n");

  for(n = 0; n <= cm->si->nrounds; n++) {
    stype       = cm->si->stype[n];
    search_opts = cm->si->search_opts[n];
    cutoff_type = cm->si->cutoff_type[n];
    sc_cutoff   = cm->si->sc_cutoff[n];
    e_cutoff    = cm->si->e_cutoff[n];
    smx         = cm->si->smx[n];
    hsi         = cm->si->hsi[n];

    /* Determine configuration of CM and CP9 based on cm->flags & cm->search_opts */
    CM2ExpMode(cm, search_opts, &cm_mode, &cp9_mode); 
    exp_mode = (stype == SEARCH_WITH_CM) ? cm_mode : cp9_mode; 

    prv_surv_fract = surv_fract;
    do_pad = (n == cm->si->nrounds) ? FALSE : TRUE;
    surv_fract = E2SurvFract(e_cutoff, (stype == SEARCH_WITH_CM) ? smx->W : cm->W, cfg->avg_hit_len, cfg->dbsize, do_pad);

    if(pre_search_mode) { 
      if((status = estimate_search_time_for_round(go, cfg, errbuf, cm, stype, search_opts, smx, r, &sec_per_res)) != eslOK) return status;
      psec = prv_surv_fract * (double) cfg->dbsize * sec_per_res;
      /*printf("psec: %f\nprv_surv_fract: %f\n", psec, prv_surv_fract);*/
#ifdef HAVE_MPI
    /* if we're paralellized take that into account */
    if(esl_opt_GetBoolean(go, "--mpi") && cfg->nproc > 1) psec /= (cfg->nproc-1);
#endif
    /* if we're only forecasting the time, divide by <n> from --forecast <n>, which is theoretical number of processors */
    if( esl_opt_IsOn(go, "--forecast")) { 
      if(esl_opt_GetInteger(go, "--forecast") > 1) { 
	psec /= (esl_opt_GetInteger(go, "--forecast") - 1);
      }
    }
    total_psec += psec;

    FormatTimeString(time_buf, psec, TRUE);
      if(n == 0) { 
	fprintf(stdout, "# %3s  %3s  %3s  %3s  %5s  %19s  %26s\n",               ""    , "",    "",    "",    "",  "      cutoffs      ",   "       predictions        ");
	fprintf(stdout, "# %3s  %3s  %3s  %3s  %5s  %19s  %26s\n",               "",     "",    "",    "",    "",  "-------------------",   "--------------------------");
	fprintf(stdout, "# %3s  %3s  %3s  %3s  %5s  %10s  %7s  %7s  %17s\n", "rnd",  "mod", "alg", "cfg", "beta",  "E value",    "bit sc",  "surv",    "time (hr:min:sec)");
	fprintf(stdout, "# %3s  %3s  %3s  %3s  %5s  %10s  %7s  %7s  %17s\n", "---",  "---", "---", "---", "-----", "----------", "-------", "-------", "-----------------");
	if(cfg->tfp != NULL) { 
	  fprintf(cfg->tfp, "# %3s  %3s  %3s  %3s  %5s  %19s  %26s\n",               ""    , "",    "",    "",    "",  "      cutoffs      ",   "       predictions        ");
	  fprintf(cfg->tfp, "# %3s  %3s  %3s  %3s  %5s  %19s  %26s\n",               "",     "",    "",    "",    "",  "-------------------",   "--------------------------");
	  fprintf(cfg->tfp, "# %3s  %3s  %3s  %3s  %5s  %10s  %7s  %7s  %17s\n", "rnd",  "mod", "alg", "cfg", "beta",  "E value",    "bit sc",  "surv",    "time (hr:min:sec)");
	  fprintf(cfg->tfp, "# %3s  %3s  %3s  %3s  %5s  %10s  %7s  %7s  %17s\n", "---",  "---", "---", "---", "-----", "----------", "-------", "-------", "-----------------");
	}
      }
      
      fprintf(stdout, "  %3d", (n+1));
      if(cfg->tfp != NULL) fprintf(cfg->tfp, "# %3d", (n+1)); /* note: only line that's printed differently in tfp and stdout, add a prefix \# */
      if(stype == SEARCH_WITH_CM) { 
	fprintf(stdout, "  %3s  %3s  %3s  ", "cm", ((search_opts & CM_SEARCH_INSIDE) ? "ins" : "cyk"), ((cm->flags & CMH_LOCAL_BEGIN) ? "loc" : "glc"));
	if(cfg->tfp != NULL) fprintf(cfg->tfp, "  %3s  %3s  %3s  ", "cm", ((search_opts & CM_SEARCH_INSIDE) ? "ins" : "cyk"), ((cm->flags & CMH_LOCAL_BEGIN) ? "loc" : "glc"));
	use_qdb  = (smx->dmin == NULL && smx->dmax == NULL) ? FALSE : TRUE;
	if(use_qdb) { 
	  fprintf(stdout, "%5g", smx->beta_qdb);
	  if(cfg->tfp != NULL) fprintf(cfg->tfp, "%5g", smx->beta_qdb);
	}
	else {
	  fprintf(stdout, "%5s", "-");
	  if(cfg->tfp != NULL) fprintf(cfg->tfp, "%5s", "-");
	}
      }
      else { 
	fprintf(stdout, "  %3s  %3s  %3s  %5s", "hmm", ((search_opts & CM_SEARCH_HMMFORWARD) ? "fwd" : "vit"), ((cm->cp9->flags & CPLAN9_LOCAL_BEGIN) ? "loc" : "glc"), "-");
	if(cfg->tfp != NULL) fprintf(cfg->tfp, "  %3s  %3s  %3s  %5s", "hmm", ((search_opts & CM_SEARCH_HMMFORWARD) ? "fwd" : "vit"), ((cm->cp9->flags & CPLAN9_LOCAL_BEGIN) ? "loc" : "glc"), "-");
      }
      if(e_cutoff < -0.1)  if((status = Score2MaxE(cm, errbuf, exp_mode, sc_cutoff, &e_cutoff)) != eslOK) return status;
      if(e_cutoff < 0.01)  { 
	fprintf(stdout, "  %10.1e", e_cutoff);
	if(cfg->tfp != NULL) fprintf(cfg->tfp, "  %10.1e", e_cutoff);
      }
      else {
	fprintf(stdout, "  %10.3f", e_cutoff);
  	if(cfg->tfp != NULL) fprintf(cfg->tfp, "  %10.3f", e_cutoff);
      }
      
      fprintf(stdout, "  %7.2f", sc_cutoff);
      if(cfg->tfp != NULL) fprintf(cfg->tfp, "  %7.2f", sc_cutoff);
      if(surv_fract < 0.0001) { 
	fprintf(stdout, "  %7.1e", surv_fract);
	if(cfg->tfp != NULL) fprintf(cfg->tfp, "  %7.1e", surv_fract); 
      }
      else { 
	fprintf(stdout, "  %7.4f", surv_fract);
	if(cfg->tfp != NULL) fprintf(cfg->tfp, "  %7.4f", surv_fract);
      }
      fprintf(stdout, "  %17s\n", time_buf);
      if(cfg->tfp != NULL) fprintf(cfg->tfp, "  %17s\n", time_buf);

      if(n != 0 && n == cm->si->nrounds) { /* print total expected run time */
	FormatTimeString(time_buf, total_psec, TRUE);
	fprintf(stdout, "# %3s  %3s  %3s  %3s  %5s  %10s  %7s  %7s  %17s\n", "---",  "---", "---", "---", "-----", "----------", "-------", "-------", "-----------------");
	fprintf(stdout, "  %3s  %3s  %3s  %3s  %5s  %10s  %7s  %7s  %17s\n", "all",  "-",   "-",   "-",   "-",     "-",          "-",       "-",       time_buf);
	fprintf(stdout, "#\n");
	if(cfg->tfp != NULL) { 
	  fprintf(cfg->tfp, "# %3s  %3s  %3s  %3s  %5s  %10s  %7s  %7s  %17s\n", "---",  "---", "---", "---", "-----", "----------", "-------", "-------", "-----------------");
	  fprintf(cfg->tfp, "# %3s  %3s  %3s  %3s  %5s  %10s  %7s  %7s  %17s\n", "all",  "-",   "-",   "-",   "-",     "-",          "-",       "-",       time_buf);
	  fprintf(cfg->tfp, "#\n");
	}
      }
    }
    else { /* search is done */
      if(n == 0) { 
	fprintf(stdout, "# %3s  %3s  %3s  %3s  %5s  %19s  %17s\n",               ""    , "",    "",    "",    "",  "  number of hits   ",   "  surv fraction  ");
	fprintf(stdout, "# %3s  %3s  %3s  %3s  %5s  %19s  %17s\n",               "",     "",    "",    "",    "",  "-------------------",   "-----------------");
	fprintf(stdout, "# %3s  %3s  %3s  %3s  %5s  %10s  %7s  %8s  %7s\n", "rnd",  "mod", "alg", "cfg", "beta",  "expected",    "actual", "expected", "actual");
	fprintf(stdout, "# %3s  %3s  %3s  %3s  %5s  %10s  %7s  %8s  %7s\n", "---",  "---", "---", "---", "-----", "----------", "-------", "--------", "-------");
	if(cfg->tfp != NULL) { 
	  fprintf(cfg->tfp, "# %3s  %3s  %3s  %3s  %5s  %19s  %17s\n",               ""    , "",    "",    "",    "",  "  number of hits   ",   "  surv fraction  ");
	  fprintf(cfg->tfp, "# %3s  %3s  %3s  %3s  %5s  %19s  %17s\n",               "",     "",    "",    "",    "",  "-------------------",   "-----------------");
	  fprintf(cfg->tfp, "# %3s  %3s  %3s  %3s  %5s  %10s  %7s  %8s  %7s\n", "rnd",  "mod", "alg", "cfg", "beta",  "expected",    "actual", "expected", "actual");
	  fprintf(cfg->tfp, "# %3s  %3s  %3s  %3s  %5s  %10s  %7s  %8s  %7s\n", "---",  "---", "---", "---", "-----", "----------", "-------", "--------", "-------");
	}
      }
      
      fprintf(stdout, "  %3d", (n+1));
      if(cfg->tfp != NULL) fprintf(cfg->tfp, "# %3d", (n+1)); /* note: only line that's printed differently in tfp and stdout, add a prefix \# */
      if(stype == SEARCH_WITH_CM) { 
	fprintf(stdout, "  %3s  %3s  %3s  ", "cm", ((search_opts & CM_SEARCH_INSIDE) ? "ins" : "cyk"), ((cm->flags & CMH_LOCAL_BEGIN) ? "loc" : "glc"));
	if(cfg->tfp != NULL) fprintf(cfg->tfp, "  %3s  %3s  %3s  ", "cm", ((search_opts & CM_SEARCH_INSIDE) ? "ins" : "cyk"), ((cm->flags & CMH_LOCAL_BEGIN) ? "loc" : "glc"));
	use_qdb  = (smx->dmin == NULL && smx->dmax == NULL) ? FALSE : TRUE;
	if(use_qdb) { 
	  fprintf(stdout, "%5g", smx->beta_qdb);
	  if(cfg->tfp != NULL) fprintf(cfg->tfp, "%5g", smx->beta_qdb);
	}
	else {
	  fprintf(stdout, "%5s", "-");
	  if(cfg->tfp != NULL) fprintf(cfg->tfp, "%5s", "-");
	}
      }
      else { 
	fprintf(stdout, "  %3s  %3s  %3s  %5s", "hmm", ((search_opts & CM_SEARCH_HMMFORWARD) ? "fwd" : "vit"), ((cm->cp9->flags & CPLAN9_LOCAL_BEGIN) ? "loc" : "glc"), "-");
	if(cfg->tfp != NULL) fprintf(cfg->tfp, "  %3s  %3s  %3s  %5s", "hmm", ((search_opts & CM_SEARCH_HMMFORWARD) ? "fwd" : "vit"), ((cm->cp9->flags & CPLAN9_LOCAL_BEGIN) ? "loc" : "glc"), "-");
      }
      if(e_cutoff < -0.1)  if((status = Score2MaxE(cm, errbuf, exp_mode, sc_cutoff, &e_cutoff)) != eslOK) return status;
      if(e_cutoff < 0.01)  { 
	fprintf(stdout, "  %10.1e", e_cutoff);
	if(cfg->tfp != NULL) fprintf(cfg->tfp, "  %10.1e", e_cutoff);
      }
      else {
	fprintf(stdout, "  %10.3f", e_cutoff);
	if(cfg->tfp != NULL) fprintf(cfg->tfp, "  %10.3f", e_cutoff);
      }
      
      fprintf(stdout, "  %7d", cm_nhitsA[n]);
      if(cfg->tfp != NULL) fprintf(cfg->tfp, "  %7d", cm_nhitsA[n]);
      if(surv_fract < 0.0001) { 
	fprintf(stdout, "  %8.1e", surv_fract);
	if(cfg->tfp != NULL) fprintf(cfg->tfp, "  %8.1e", surv_fract);
      }
      else {
	fprintf(stdout, "  %8.4f", surv_fract);
	if(cfg->tfp != NULL) fprintf(cfg->tfp, "  %8.4f", surv_fract);
      }
      if(cm_surv_fractA[n] < 0.0001) { 
	fprintf(stdout, "  %7.1e\n", cm_surv_fractA[n]);
	if(cfg->tfp != NULL) fprintf(cfg->tfp, "  %7.1e\n", cm_surv_fractA[n]);
      }
      else {
	fprintf(stdout, "  %7.4f\n", cm_surv_fractA[n]);
	if(cfg->tfp != NULL) fprintf(cfg->tfp, "  %7.4f\n", cm_surv_fractA[n]);
      }

      if(n == cm->si->nrounds) { /* print total expected run time */
	fprintf(stdout, "#\n");
	if(cfg->tfp != NULL) fprintf(cfg->tfp, "#\n");
	/*fprintf(stdout, "# %24s\n",       "   total search time    ");*/
	/*fprintf(stdout, "# %24s\n",       "------------------------");*/
	fprintf(stdout, "# %13s  %13s  (hr:min:sec)\n", "expected time",  "actual time");
	fprintf(stdout, "# %13s  %13s\n", "-------------",  "-------------");
	if(cfg->tfp != NULL) { 
	  fprintf(cfg->tfp, "# %13s  %13s  (hr:min:sec)\n", "expected time",  "actual time");
	  fprintf(cfg->tfp, "# %13s  %13s\n", "-------------",  "-------------");
	}
	FormatTimeString(time_buf, in_total_psec, TRUE);
	fprintf(stdout, "  %13s", time_buf);
	if(cfg->tfp != NULL) fprintf(cfg->tfp, "# %13s", time_buf); /* note: 1 of only 2 line that's printed differently in tfp and stdout, add a prefix \# */
	FormatTimeString(time_buf, in_asec, TRUE);
	fprintf(stdout, "  %13s\n", time_buf);
	if(cfg->tfp != NULL) fprintf(cfg->tfp, "  %13s\n", time_buf);
      }
    }
  }
  if(cfg->tfp != NULL) fflush(cfg->tfp);
  esl_randomness_Destroy(r);
  
  if(ret_total_psec != NULL) *ret_total_psec = total_psec;
  return eslOK;
}  

/* Function: print_searchinfo_for_uncalibrated_cm
 * Date:     EPN, Thu May 17 14:47:36 2007
 * Purpose:  Print info about search (cutoffs, algorithm, etc.). 
 *           with a CM that does not have exp tail stats. Can be called in 
 *           2 different modes, mode 1 is 'pre-search', called prior
 *           before a search, mode 2 is 'post-search', called after
 *           a search is done. We're in mode 1 iff cm_surv_fractA == NULL
 *           and cm_nhitsA == NULL.
 */
int print_searchinfo_for_uncalibrated_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float *cm_surv_fractA, int *cm_nhitsA, double in_asec)
{
  int n;
  float sc_cutoff;
  int stype;
  int search_opts;
  ScanMatrix_t *smx;
  int   use_qdb;            /* are we using qdb for current round? */
  int   pre_search_mode;    /* TRUE if this function was called before the search was run, FALSE if called after */
  char  time_buf[128];	    /* for printing run time */

  /* Could use ESL_GETOPTS here, but using the CM flags assures we're reporting
   * on how the CM is actually config'ed, not how we want it to be
   */
  pre_search_mode = (cm_surv_fractA == NULL && cm_nhitsA == NULL) ? TRUE : FALSE;

  if(cm->flags & CMH_EXPTAIL_STATS) ESL_FAIL(eslEINCOMPAT, errbuf, "print_searchinfo_for_uncalibrated_cm(): but cm: %s has exp tail stats.", cm->name);

  if(pre_search_mode) { 
    fprintf(stdout, "#\n# Pre-search info for CM %d: %s\n", cfg->ncm, cm->name);
    if(cfg->tfp != NULL) fprintf(cfg->tfp, "#\n# Pre-search info for CM %d: %s\n", cfg->ncm, cm->name);
  }
  else {
    fprintf(stdout, "#\n# Post-search info for CM %d: %s\n", cfg->ncm, cm->name);
    if(cfg->tfp != NULL) fprintf(cfg->tfp, "#\n# Post-search info for CM %d: %s\n", cfg->ncm, cm->name);
  }
  fprintf(stdout, "#\n");
  if(cfg->tfp != NULL) fprintf(cfg->tfp, "#\n");

  for(n = 0; n <= cm->si->nrounds; n++) {
    stype       = cm->si->stype[n];
    search_opts = cm->si->search_opts[n];
    sc_cutoff   = cm->si->sc_cutoff[n];
    smx         = cm->si->smx[n];

    use_qdb     = (smx == NULL || (smx->dmin == NULL && smx->dmax == NULL)) ? FALSE : TRUE;

    if(n == 0) { 
      if(pre_search_mode) { 
	fprintf(stdout, "# %3s  %3s  %3s  %3s  %5s  %10s\n", "rnd",  "mod", "alg", "cfg", "beta",  "bit sc cut");
	fprintf(stdout, "# %3s  %3s  %3s  %3s  %5s  %10s\n", "---",  "---", "---", "---", "-----", "----------");
	if(cfg->tfp != NULL) { 
	  fprintf(cfg->tfp, "# %3s  %3s  %3s  %3s  %5s  %10s\n", "rnd",  "mod", "alg", "cfg", "beta",  "bit sc cut");
	  fprintf(cfg->tfp, "# %3s  %3s  %3s  %3s  %5s  %10s\n", "---",  "---", "---", "---", "-----", "----------");
	}
      }
      else { /* post search */
	fprintf(stdout, "# %3s  %3s  %3s  %3s  %5s  %10s  %8s  %10s\n", "rnd",  "mod", "alg", "cfg", "beta",  "bit sc cut", "num hits", "surv fract");
	fprintf(stdout, "# %3s  %3s  %3s  %3s  %5s  %10s  %8s  %10s\n", "---",  "---", "---", "---", "-----", "----------", "--------", "----------");
	if(cfg->tfp != NULL) { 
	  fprintf(cfg->tfp, "# %3s  %3s  %3s  %3s  %5s  %10s  %8s  %10s\n", "rnd",  "mod", "alg", "cfg", "beta",  "bit sc cut", "num hits", "surv fract");
	  fprintf(cfg->tfp, "# %3s  %3s  %3s  %3s  %5s  %10s  %8s  %10s\n", "---",  "---", "---", "---", "-----", "----------", "--------", "----------");
	}
      }
    }
    fprintf(stdout, "  %3d", (n+1));
    if(cfg->tfp != NULL) fprintf(cfg->tfp, "# %3d", (n+1)); /* note: 1 of only 2 lines that's printed differently in tfp and stdout, add a prefix \# */
    if(stype == SEARCH_WITH_CM) { 
      fprintf(stdout, "  %3s  %3s  %3s  ", "cm", ((search_opts & CM_SEARCH_INSIDE) ? "ins" : "cyk"), ((cm->flags & CMH_LOCAL_BEGIN) ? "loc" : "glc"));
      if(cfg->tfp != NULL) fprintf(cfg->tfp, "  %3s  %3s  %3s  ", "cm", ((search_opts & CM_SEARCH_INSIDE) ? "ins" : "cyk"), ((cm->flags & CMH_LOCAL_BEGIN) ? "loc" : "glc"));
      if(use_qdb) { 
	fprintf(stdout, "%5g", smx->beta_qdb);
	if(cfg->tfp != NULL) fprintf(cfg->tfp, "%5g", smx->beta_qdb);
      }
      else {
	fprintf(stdout, "%5s", "-");
	if(cfg->tfp != NULL) fprintf(cfg->tfp, "%5s", "-");
      }
    }
    else { 
      fprintf(stdout, "  %3s  %3s  %3s  %5s", "hmm", ((search_opts & CM_SEARCH_HMMFORWARD) ? "fwd" : "vit"), ((cm->cp9->flags & CPLAN9_LOCAL_BEGIN) ? "loc" : "glc"), "-");
      if(cfg->tfp != NULL) fprintf(cfg->tfp, "  %3s  %3s  %3s  %5s", "hmm", ((search_opts & CM_SEARCH_HMMFORWARD) ? "fwd" : "vit"), ((cm->cp9->flags & CPLAN9_LOCAL_BEGIN) ? "loc" : "glc"), "-");
    }
    fprintf(stdout, "  %10.2f", sc_cutoff);
    if(cfg->tfp != NULL) fprintf(cfg->tfp, "  %10.2f", sc_cutoff);
    if(pre_search_mode) { 
      fprintf(stdout, "\n");
      if(cfg->tfp != NULL) fprintf(cfg->tfp, "\n");
    }
    else { /* post search */
      if(cm_surv_fractA[n] < 0.0001) { 
	fprintf(stdout, "  %8d  %10.1e\n", cm_nhitsA[n], cm_surv_fractA[n]);
	if(cfg->tfp != NULL) fprintf(cfg->tfp, "  %8d  %10.1e\n", cm_nhitsA[n], cm_surv_fractA[n]);
      }
      else {
	fprintf(stdout, "  %8d  %10.4f\n", cm_nhitsA[n], cm_surv_fractA[n]);
	if(cfg->tfp != NULL) fprintf(cfg->tfp, "  %8d  %10.4f\n", cm_nhitsA[n], cm_surv_fractA[n]);
      }
      if(n == cm->si->nrounds) { 
	FormatTimeString(time_buf, in_asec, FALSE);
	fprintf(stdout, "#\n");
	fprintf(stdout, "# %11s  (hr:min:sec)\n", "run time");
	fprintf(stdout, "# %11s\n", "-----------");
	fprintf(stdout, "  %11s\n", time_buf);
	if(cfg->tfp != NULL) { 
	  fprintf(cfg->tfp, "#\n");
	  fprintf(cfg->tfp, "# %11s  (hr:min:sec)\n", "run time");
	  fprintf(cfg->tfp, "# %11s\n", "-----------");
	  fprintf(cfg->tfp, "# %11s\n", time_buf); /* note: 1 of only 2 line that's printed differently in tfp and stdout, add a prefix \# */
	}
      }
    }
  }
  if(cfg->tfp != NULL) fflush(cfg->tfp);
  return eslOK;
}

/* Function: estimate_search_time_for_round
 * Date:     EPN, Tue Mar  4 13:18:52 2008
 * Purpose:  Estimate search time for each round of searching.
 *           This is done by actually searching a sequence with the 
 *           appropriate algorithm. The length of the sequence to
 *           search is set such that it should take about 0.1 seconds.
 */
int estimate_search_time_for_round(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, int stype, int search_opts, ScanMatrix_t *smx, ESL_RANDOMNESS *r, double *ret_sec_per_res)
{
  int    status;
  int    L;                /* length of sequence we'll generate and search to get time estimate */
  double psec_per_Mc;      /* rough prediction at seconds per Mc based on empirical run times I've witnessed, 
                            * doesn't need to be very accurate as we just use it to set length of seq to search to get real prediction */
  float  Mc;               /* millions of DP calculations we're going to do */
  float  Mc_per_res;       /* millions of dp calcs per residue, if searching with CM, corrects for first W residues requiring less dp calcs */
  int    irrelevant_W;     /* temporary W */
  int    orig_search_opts; /* cm->search_opts when function was entered */
  float  sec_per_res;      /* seconds per residue */
  float  targ_sec = 0.1;   /* target number of seconds our timing expt will take */
  int    Lmin = 400;       /* minimum number of residues to search to get timing */

  ESL_DSQ *dsq;
  ESL_STOPWATCH *w  = esl_stopwatch_Create();

  if(w == NULL)               ESL_FAIL(eslEMEM,   errbuf, "estimate_search_time_for_round(): memory error, stopwatch not created.\n");
  if(ret_sec_per_res == NULL) ESL_FAIL(eslEINVAL, errbuf, "estimate_search_time_for_round(): ret_sec_per_res is NULL");

  if(stype == SEARCH_WITH_CM) {
    if(smx == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "estimate_search_time_for_round(), stype is SEARCH_WITH_CM, but smx is NULL");
    int use_qdb     = (smx->dmin == NULL && smx->dmax == NULL) ? FALSE : TRUE;
    if(use_qdb) { if((status = cm_GetNCalcsPerResidueForGivenBeta(cm, errbuf, FALSE, smx->beta_qdb, &Mc_per_res, &irrelevant_W))  != eslOK) return status; }
    else        { if((status = cm_GetNCalcsPerResidueForGivenBeta(cm, errbuf, TRUE,  cm->beta_W,    &Mc_per_res, &irrelevant_W))  != eslOK) return status; }
    psec_per_Mc = (search_opts & CM_SEARCH_INSIDE) ? (1. /  75.) : (1. / 275.);  /*  75 Mc/S inside;  275 Mc/S CYK */
    /* determine L that will take about targ_sec seconds */
    L = targ_sec / (psec_per_Mc * Mc_per_res);
    L = ESL_MAX(L, Lmin); /* we have to search at least <Lmin> residues */
    /* now determine exactly how many dp calculations we'd do if we search L residues, 
     * this won't be the same as Mc_per_res * L b/c Mc_per_res was passed in from caller
     * and was calculated after correcting for the fact that the first W residues have fewer
     * DP calcs than all other residues, b/c d < W for first W residues.
     */  
    if((status = cm_CountSearchDPCalcs(cm, errbuf, L, smx->dmin, smx->dmax, smx->W, FALSE,  NULL, &Mc)) != eslOK) return status;
    /* FALSE says don't correct for fewer dp calcs for first W residues, we want to know how many total DP calcs
     * there will be in L residues */
    Mc *= L; /* Mc was for 1 residue, multiply by L to get Mc for L residues */
  }

  else { 
    if((status = cp9_GetNCalcsPerResidue(cm->cp9, errbuf, &Mc_per_res)) != eslOK) return status;
    psec_per_Mc = (search_opts & CM_SEARCH_HMMFORWARD) ? (1. / 175.) : (1. / 380.);  /* 175 Mc/S forward; 380 Mc/S viterbi */
    /* determine L that will take about 0.1 seconds */
    L  = 0.1 / (psec_per_Mc * Mc_per_res);
    /* how many millions of DP cells will it be? */
    Mc = Mc_per_res * L;
  }

  orig_search_opts = cm->search_opts;
  cm->search_opts = search_opts; /* we'll restore cm->search_opts to orig_search_opts at end of the function */

  ESL_ALLOC(dsq,  sizeof(ESL_DSQ) * (L +2));
  esl_rsq_xfIID(r, cm->null, cm->abc->K, L, dsq);

  esl_stopwatch_Start(w);
  if(stype == SEARCH_WITH_CM) { 
    if(search_opts & CM_SEARCH_INSIDE) { 
      if((status = FastIInsideScan(cm, errbuf, smx, dsq, 1, L, 0., NULL, (! esl_opt_GetBoolean(go, "--no-null3")), NULL, NULL)) != eslOK) return status;
    }
    else 
      if((status = FastCYKScan    (cm, errbuf, smx, dsq, 1, L, 0., NULL, (! esl_opt_GetBoolean(go, "--no-null3")), NULL, NULL)) != eslOK) return status;
  }
  else { /* search with HMM */
    if(search_opts & CM_SEARCH_HMMFORWARD) { /* forward */
      if((status = cp9_Forward(cm, errbuf, cm->cp9_mx, dsq, 1, L, cm->W, 0., NULL,
			       TRUE,   /* we're scanning */
			       FALSE,  /* we're not ultimately aligning */
			       TRUE,   /* be memory efficient */
			       (! esl_opt_GetBoolean(go, "--no-null3")),
			       NULL, NULL, NULL)) != eslOK) return status;
    }
    else { /* viterbi */
      if((status = cp9_Viterbi(cm, errbuf, cm->cp9_mx, dsq, 1, L, cm->W, 0., NULL,
			       TRUE,   /* we're scanning */
			       FALSE,  /* we're not ultimately aligning */
			       TRUE,   /* be memory efficient */
			       (! esl_opt_GetBoolean(go, "--no-null3")),
			       NULL, NULL,
			       NULL,   /* don't want traces back */
			       NULL)) != eslOK) return status;
    }
  }
  esl_stopwatch_Stop(w);
  free(dsq);

  cm->search_opts = orig_search_opts;
  sec_per_res = w->user * (Mc_per_res / Mc);
  *ret_sec_per_res = sec_per_res;

  ESL_DPRINTF1(("L: %d\n", L));
  ESL_DPRINTF1(("w->user: %f\n", w->user));
  ESL_DPRINTF1(("sec_per_res: %f\n", sec_per_res));
  ESL_DPRINTF1(("Mc_per_res: %f\n", Mc_per_res));
  ESL_DPRINTF1(("Mc: %f\n", Mc));
  esl_stopwatch_Destroy(w);

  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "estimate_search_time_for_round(): memory error.\n");
  return status; /* NEVERREACHED */
}

/* Function: dump_gc_info
 * Date:     EPN, Sun Mar 23 09:57:55 2008
 * Purpose:  Dump GC stats for the target database to a file.
 */
int dump_gc_info(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  int status;
  double *gc_freq;
  ESL_ALPHABET *tmp_abc = NULL;
  long dbsize;
  int i, j;
  int nbars;
  int nseq;

  tmp_abc = esl_alphabet_Create(eslRNA);
  if((status = GetDBInfo(tmp_abc, cfg->sqfp, errbuf, &dbsize, &nseq, &gc_freq)) != eslOK) return status; 
  esl_vec_DNorm(gc_freq, GC_SEGMENTS);
  esl_alphabet_Destroy(tmp_abc);

  fprintf(cfg->gcfp, "# %-25s %s\n",  "seqfile:", cfg->sqfile);
  fprintf(cfg->gcfp, "# %-25s %d\n",  "number of sequences:", nseq);
  fprintf(cfg->gcfp, "# %-25s %ld\n", "dbsize (nt, one strand):", dbsize);
  fprintf(cfg->gcfp, "#\n");
  fprintf(cfg->gcfp, "# %10s  %22s\n", "GC percent", "freq of 100 nt windows");
  fprintf(cfg->gcfp, "# %10s  %22s\n", "----------", "--------------");
  for (i=0; i<GC_SEGMENTS; i++) { 
    nbars = (int) (gc_freq[i] * 400);
    fprintf(cfg->gcfp, "  %10d  %22.12f  ", i, gc_freq[i]);
    for(j = 0; j < nbars && j < 40; j++) fprintf(cfg->gcfp, "=");
    fprintf(cfg->gcfp, "\n");
  }
  free(gc_freq);
  return eslOK;
}  


/* Function: overwrite_lambdas
 * Date:     EPN, Fri May  9 09:06:15 2008
 *
 * Purpose:  If --lambda <x> was enabled we overwrite the lambdas 
 *           we read from the CM file as <x>.
 *
 * Returns:  eslOK on success
 */
static int
overwrite_lambdas(const ESL_GETOPTS *go, const struct cfg_s *cfg, CM_t *cm, char *errbuf)
{
  double lambda;
  int i, p;

  if(!esl_opt_IsOn(go, "--lambda"))     ESL_FAIL(eslEINCOMPAT, errbuf, "overwrite_lambdas(), but --lambda was not enabled, shouldn't happen.\n");
  if(! (cm->flags & CMH_EXPTAIL_STATS)) ESL_FAIL(eslEINCOMPAT, errbuf, "--lambda only works with calibrated CM files. Run cmcalibrate (please).");
  if(cm->stats == NULL)                 ESL_FAIL(eslEINCOMPAT, errbuf, "overwrite_lambdas(), cm->stats is NULL, shouldn't happen.\n");

  lambda = esl_opt_GetReal(go, "--lambda");
  for(i = 0; i < EXP_NMODES; i++) { 
    for(p = 0; p < cm->stats->np; p++) { 
      cm->stats->expAA[i][p]->lambda = lambda;
      /* mu_extrap, the extrapolated mu value is lambda dependent, so we have to update it */
      cm->stats->expAA[i][p]->mu_extrap = cm->stats->expAA[i][p]->mu_orig - log(1./cm->stats->expAA[i][p]->tailp) / cm->stats->expAA[i][p]->lambda;
    }
  }
  return eslOK;
}

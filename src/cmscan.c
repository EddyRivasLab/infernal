/* cmscan: search sequence(s) against a covariance model database
 * 
 * EPN, Tue Jun 28 04:33:27 2011
 * SRE, Mon Oct 20 08:28:05 2008 [Janelia] (hmmscan.c)
 * SVN $Id$
 */
#include "esl_config.h"
#include "p7_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#ifdef HAVE_MPI
#include "mpi.h"
#include "esl_mpi.h"
#endif /*HAVE_MPI*/

#ifdef SUB_WITH_HMMER_THREADS
#include <unistd.h>
#include "esl_threads.h"
#include "esl_workqueue.h"
#endif /*SUB_WITH_HMMER_THREADS*/

#include "hmmer.h"
#include "funcs.h"		/* external functions                   */
#include "structs.h"		/* data structures, macros, #define's   */

typedef struct {
#ifdef SUB_WITH_HMMER_THREADS
  ESL_WORK_QUEUE   *queue;
#endif /*SUB_WITH_HMMER_THREADS*/
  ESL_SQ           *qsq;
  CM_PIPELINE      *pli;         /* work pipeline                           */
  CM_TOPHITS       *th;          /* top hit results                         */
} WORKER_INFO;

#define REPOPTS         "-E,-T,--cut_ga,--cut_nc,--cut_tc"
#define INCOPTS         "--incE,--incT,--cut_ga,--cut_nc,--cut_tc"
#define THRESHOPTS      "-E,-T,--incE,--incT,--cut_ga,--cut_nc,--cut_tc"
#define XFASTMIDMAXOPTS "--max,--F1,--F2,--F3,--F4,--F5,--F6,--noF1,--noF2,--noF3,--nogfwd,--nocyk,--nohmm,--hmm"
#define TIMINGOPTS      "--time-F1,--time-F2,--time-F3,--time-dF3,--time-bfil,--time-F4"

#if defined (SUB_WITH_HMMER_THREADS) && defined (HAVE_MPI)
#define CPUOPTS     "--mpi"
#define MPIOPTS     "--cpu"
#else
#define CPUOPTS     NULL
#define MPIOPTS     NULL
#endif

#ifdef HAVE_MPI
#define DAEMONOPTS  "-o,--tblout,--domtblout,--mpi,--stall"
#else
#define DAEMONOPTS  "-o,--tblout,--domtblout"
#endif

static ESL_OPTIONS options[] = {
  /* name           type          default  env  range toggles  reqs   incomp                         help                                           docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "show brief help on version and usage",                          1 },
  /* Control of output */
  { "-o",           eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "direct output to file <f>, not stdout",                         2 },
  { "--tblout",     eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save parseable table of per-sequence hits to file <s>",         2 },
  { "--acc",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "prefer accessions over names in output",                        2 },
  { "--noali",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "don't output alignments, so output is smaller",                 2 },
  { "--notextw",    eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL, "--textw",        "unlimit ASCII text output line width",                          2 },
  { "--textw",      eslARG_INT,    "120", NULL, "n>=120",NULL,  NULL, "--notextw",      "set max width of ASCII text output lines",                      2 },
  /* Control of reporting thresholds */
  { "-E",           eslARG_REAL,  "10.0", NULL, "x>0",   NULL,  NULL,  REPOPTS,         "report models <= this E-value threshold in output",             4 },
  { "-T",           eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  REPOPTS,         "report models >= this score threshold in output",               4 },
  /* Control of inclusion (significance) thresholds: */
  { "--incE",       eslARG_REAL,  "0.01", NULL, "x>0",   NULL,  NULL,  INCOPTS,         "consider models <= this E-value threshold as significant",      5 },
  { "--incT",       eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  INCOPTS,         "consider models >= this score threshold as significant",        5 },
  /* Model-specific thresholding for both reporting and inclusion */
  { "--cut_ga",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's GA gathering cutoffs to set all thresholding",    6 },
  { "--cut_nc",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's NC noise cutoffs to set all thresholding",        6 },
  { "--cut_tc",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's TC trusted cutoffs to set all thresholding",      6 },
  /* Control of acceleration pipeline */
  /* Control of acceleration pipeline */
  { "--fZ",         eslARG_REAL,    NULL, NULL, NULL,    NULL,  NULL,  "--rfam",        "set heuristic filters to defaulst used for a db of size <x> Mb", 7 },
  { "--rfam",       eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  "--fZ",          "set heuristic filters at Rfam-level (more speed, less power)", 7 },
  { "--max",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  XFASTMIDMAXOPTS, "turn all heuristic filters off  (less speed, more power)",     7 },
  { "--noenvdef",   eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  "--hmm",         "do not define domains inside windows prior to CM stages",      7 },
  { "--msvmerge",   eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  "--noF1,--nohmm","merge MSV windows prior to Viterbi filter",                   7 },
  { "--pad",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  "--hmm,--noenvdef","pad domains i..j to j-W+1..i+W-1",                             7 },
  { "--nohmm",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--max,--hmm",    "skip all HMM filter stages (MSV/Vit/Fwd/gFwd/envdef)",         7 },
  { "--noF1",       eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--max",          "skip the MSV filter stage",                                    7 },
  { "--shortmsv",   eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--max,--noF1",  "run MSV on short 2*W subseqs, not longer subseqs",             7 },
  { "--noF2",       eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--max",          "skip the Viterbi filter stage",                                7 },
  { "--noF3",       eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--max",          "skip the Forward filter stage",                                7 },
  { "--noF4",       eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--max",          "skip the glocal Forward filter stage",                         7 },
  { "--noF6",       eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--hmm",          "skip the CYK filter stage",                                    7 },
  { "--doF1b",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--max,--noF1",  "turn on the MSV composition bias filter",                      7 },
  { "--noF2b",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--max,--noF2",  "turn on the Vit composition bias filter",                      7 },
  { "--noF3b",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--max,--noF3",  "turn on the Fwd composition bias filter",                      7 },
  { "--noF4b",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--max,--noF3",  "turn on the glocal Fwd composition bias filter",               7 },
  { "--noF5b",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--max,--noenvdef", "turn on the per-domain composition bias filter",               7 },
  { "--F1",         eslARG_REAL,  "0.35", NULL, NULL,    NULL,  NULL, "--max",          "Stage 1 (MSV) threshold: promote hits w/ P <= F1",             7 },
  { "--F1b",        eslARG_REAL,  "0.35", NULL, NULL,    NULL,"--domsvbias", "--max",   "Stage 1 (MSV) bias threshold: promote hits w/ P <= F1b",       7 },
  { "--F2",         eslARG_REAL,  "0.20", NULL, NULL,    NULL,  NULL, "--max",          "Stage 2 (Vit) threshold: promote hits w/ P <= F2",             7 },
  { "--F2b",        eslARG_REAL,  "0.20", NULL, NULL,    NULL,  NULL, "--noF2b,--max","Stage 2 (Vit) bias threshold: promote hits w/ P <= F2b",   7 },
  { "--F3",         eslARG_REAL,  "0.003", NULL, NULL,    NULL,  NULL, "--max",          "Stage 3 (Fwd) threshold: promote hits w/ P <= F3",             7 },
  { "--F3b",        eslARG_REAL,  "0.003", NULL, NULL,    NULL,  NULL, "--noF3b,--max","Stage 3 (Fwd) bias threshold: promote hits w/ P <= F3b",       7 },
  { "--F4",         eslARG_REAL,  "0.003", NULL, NULL,    NULL,  NULL, "--max",          "Stage 4 (gFwd) glocal threshold: promote hits w/ P <= F4", 7 },
  { "--F4b",        eslARG_REAL,  "0.003", NULL, NULL,    NULL,  NULL, "--nogfwdbias,--max","Stage 4 (gFwd) glocal bias thr: promote hits w/ P <= F4b", 7 },
  { "--F5",         eslARG_REAL,  "0.003", NULL, NULL,    NULL,  NULL, "--max",          "Stage 5 (env defn) threshold: promote hits w/ P <= F5", 7 },
  { "--F5b",        eslARG_REAL,  "0.003", NULL, NULL,    NULL,  NULL, "--noedefbias,--max",  "Stage 5 (env defn) bias thr: promote hits w/ P <= F5b", 7 },
  { "--F6",         eslARG_REAL,  "1e-4", NULL, NULL,    NULL,  NULL, "--max,--nocyk,--hmm","Stage 6 (CYK) threshold: promote hits w/ P <= F6",         7 },
  { "--cykenvx",    eslARG_INT,     "10", NULL, "n>=1",  NULL,  NULL, "--max,--nocyk,--hmm","CYK envelope redefinition threshold multiplier, <n> * F4", 7 },
  { "--nocykenv",   eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--max,--nocyk,--hmm","Do not redefine envelopes after stage 4 based on CYK hits",7 },
  { "--time-F1",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, TIMINGOPTS,        "abort after Stage 1 MSV; for timings",                        7 },
  { "--time-F2",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, TIMINGOPTS,        "abort after Stage 2 Vit; for timings",                        7 },
  { "--time-F3",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, TIMINGOPTS,        "abort after Stage 3 Fwd; for timings",                        7 },
  { "--time-F4",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, TIMINGOPTS,        "abort after Stage 4 glocal Fwd; for timings",                 7 },
  { "--time-F5",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, TIMINGOPTS,        "abort after Stage 5 envelope def; for timings",               7 },
  { "--time-F6",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, TIMINGOPTS,        "abort after Stage 6 CYK; for timings",                        7 },
  { "--rt1",        eslARG_REAL,  "0.25", NULL, NULL,    NULL,  NULL, "--nohmm,--noenvdef","set domain definition rt1 parameter as <x>",                  7 },
  { "--rt2",        eslARG_REAL,  "0.10", NULL, NULL,    NULL,  NULL, "--nohmm,--noenvdef","set domain definition rt2 parameter as <x>",                  7 },
  { "--rt3",        eslARG_REAL,  "0.20", NULL, NULL,    NULL,  NULL, "--nohmm,--noenvdef","set domain definition rt3 parameter as <x>",                  7 },
  { "--ns",         eslARG_INT,   "200",  NULL, NULL,    NULL,  NULL, "--nohmm,--noenvdef","set number of domain tracebacks to <n>",                      7 },
  { "--localenv",   eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--nohmm,--noenvdef","define domains with HMM in local (not glocal) mode",          7 },
  { "--wnosplit",   eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--noF1",         "do not split windows after MSV stage", 7 },
  { "--wmult",      eslARG_REAL,   "3.0", NULL, NULL,    NULL,  NULL, "--wnosplit,--noF1",     "scalar multiplier for flagging window to split (if --wsplit)", 7 },
  { "--wcorr",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, NULL,              "use window size correction for Vit/Fwd filters", 7 },
  { "--nocorr",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, NULL,              "use no  correction for domain definition", 7 },
  { "--oldcorr",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--nocorr",        "use old correction for domain definition", 7 },
  { "--envhitbias", eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--noedefbias",     "calc domain bias for only the domain, not entire window", 7 },
  { "--nogreedy",   eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,             "do not resolve hits with greedy algorithm, use optimal one", 7 },
  { "--doml",       eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,             "filter with a ML p7 HMM in addition to additional p7 HMMs in the file", 7 },
  { "--noadd",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,             "do not filter with any additional p7 HMMs in the file", 7 },
  { "--filcmW",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,             "use CM's window length for all HMM filters", 7 },
  { "--glen",       eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, NULL,              "use length dependent glocal p7 filter P-value thresholds", 7},
  { "--glN",        eslARG_INT,   "201",  NULL, NULL,    NULL,"--glen",NULL,             "minimum value to start len-dependent glocal threshold", 7},
  { "--glX",        eslARG_INT,   "500",  NULL, NULL,    NULL,"--glen",NULL,             "maximum value for len-dependent glocal threshold", 7},
  { "--glstep",     eslARG_INT,   "100",  NULL, NULL,    NULL,"--glen",NULL,             "for len-dependent glocal thr, step size for halving thr", 7},
  /* Other options */
  { "--null2",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "turn on biased composition score corrections",               12 },
  { "-Z",           eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL,  NULL,            "set # of comparisons done, for E-value calculation",           12 },
  { "--seed",       eslARG_INT,    "42",  NULL, "n>=0",  NULL,  NULL,  NULL,            "set RNG seed to <n> (if 0: one-time arbitrary seed)",          12 },
  { "--qformat",    eslARG_STRING,  NULL, NULL, NULL,    NULL,  NULL,  NULL,            "assert input <seqfile> is in format <s>: no autodetection",    12 },
  { "--daemon",     eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL,  DAEMONOPTS,      "run program as a daemon",                                      12 },
  /* options affecting the alignment of hits */
  { "--aln-cyk",      eslARG_NONE, FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "align hits with CYK", 8 },
  { "--aln-nonbanded",eslARG_NONE, FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "do not use HMM bands when aligning hits", 8 },
  { "--aln-sizelimit",eslARG_REAL,"128.", NULL, "x>0",   NULL,  NULL,  NULL,            "set maximum allowed size of DP matrices for hit alignment to <x> Mb", 8 },
  /* Options taken from infernal 1.0.2 cmsearch */
  /* options for algorithm for final round of search */
  { "-g",             eslARG_NONE,    FALSE,     NULL, NULL,    NULL,        NULL,            NULL, "configure CM for glocal alignment [default: local]", 1 },
  { "--cyk",          eslARG_NONE,    FALSE,     NULL, NULL,    NULL,        NULL, "--nocyk,--hmm", "use scanning CM CYK algorithm", 20 },
  { "--hmm",          eslARG_NONE,    FALSE,     NULL, NULL,    NULL,        NULL,         "--cyk",  "do not use the CM, use only the HMM", 20},
  /* banded options for CYK filter round of searching */
  { "--ftau",         eslARG_REAL,    "1e-4",    NULL, "0<x<1", NULL,        NULL,"--fqdb,--hmm,--nocyk",  "set tail loss prob for --fhbanded to <x>", 20 },
  { "--fsums",        eslARG_NONE,    FALSE,     NULL, NULL,    NULL,        NULL,"--fqdb,--hmm,--nocyk",  "w/--fhbanded use posterior sums (widens bands)", 20 },
  { "--fqdb",         eslARG_NONE,    FALSE,     NULL, NULL,    NULL,        NULL,"--hmm,--nocyk",  "use QDBs in CYK filter round, not HMM bands", 20 },
  { "--fbeta",        eslARG_REAL,    "1e-9",    NULL, "0<x<1", NULL,    "--fqdb","--hmm,--nocyk",  "set tail loss prob for CYK filter QDB calculation to <x>", 20 },
  { "--fnonbanded",   eslARG_NONE,    FALSE,     NULL, NULL,    NULL,    "--fqdb","--hmm,--nocyk",  "do not use any bands for CYK filter round", 20},
  /* banded options for final round of searching */
  { "--tau",          eslARG_REAL,   "1e-7",     NULL, "0<x<1", NULL,        NULL,       "--qdb,--nonbanded,--hmm", "set tail loss prob for --hbanded to <x>", 20 },
  { "--sums",         eslARG_NONE,    FALSE,     NULL, NULL,    NULL,        NULL,       "--qdb,--nonbanded,--hmm", "w/--hbanded use posterior sums (widens bands)", 20 },
  { "--qdb",          eslARG_NONE,    FALSE,     NULL, NULL,    NULL,        NULL,"--nonbanded,--hmm,--tau", "use QDBs (instead of HMM bands) in final Inside round", 20 },
  { "--beta",         eslARG_REAL,   "1e-15",    NULL, "0<x<1",  NULL,     "--qdb",        "--hmm", "set tail loss prob for final Inside QDB calculation to <x>", 20 },
  { "--nonbanded",    eslARG_NONE,    FALSE,     NULL, NULL,    NULL,        NULL,"--hmm,--tau,--sums,--beta", "do not use QDBs or HMM bands in final Inside round of CM search", 20 },
  /* other infernal 1.0.2 options */
  { "--toponly",      eslARG_NONE,    FALSE,     NULL, NULL,    NULL,        NULL,            NULL, "only search the top strand", 20 },
  { "--bottomonly",   eslARG_NONE,    FALSE,     NULL, NULL,    NULL,        NULL,            NULL, "only search the bottom strand", 20 },
  { "--nonull3",      eslARG_NONE,    FALSE,     NULL, NULL,    NULL,        NULL,            NULL, "turn OFF the NULL3 post hoc additional null model", 20 },
  /* experimental options */
  { "--cp9noel",      eslARG_NONE,    FALSE,     NULL, NULL,    NULL,        NULL,            "-g", "turn OFF local ends in cp9 HMMs", 20 },
  { "--cp9gloc",      eslARG_NONE,    FALSE,     NULL, NULL,    NULL,        NULL,  "-g,--cp9noel", "configure CP9 HMM in glocal mode", 20 },
#ifdef SUB_WITH_HMMER_THREADS
  { "--cpu",        eslARG_INT, NULL,"HMMER_NCPU","n>=0",NULL,  NULL,  CPUOPTS,         "number of parallel CPU workers to use for multithreads",       12 },
#endif
#ifdef HAVE_MPI
  { "--stall",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,"--mpi", NULL,            "arrest after start: for debugging MPI under gdb",              12 },  
  { "--mpi",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  MPIOPTS,         "run as an MPI parallel program",                               12 },
#endif
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};


/* struct cfg_s : "Global" application configuration shared by all threads/processes
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 */
struct cfg_s {
  char            *seqfile;           /* query sequence file                             */
  char            *cmfile;            /* database CM file                                */

  int              do_mpi;            /* TRUE if we're doing MPI parallelization         */
  int              nproc;             /* how many MPI processes, total                   */
  int              my_rank;           /* who am I, in 0..nproc-1                         */

  /* TEMPORARY?: Decide how handle Z */
  int64_t          Z;                /* database size, in number of models               */
  int              Z_setby;          /* how Z was set: CM_ZSETBY_SSIINFO, CM_ZSETBY_OPTION, CM_ZSETBY_FILEINFO */

};

static char usage[]  = "[-options] <cmdb> <seqfile>";
static char banner[] = "search sequence(s) against a profile database";

static int  serial_master(ESL_GETOPTS *go, struct cfg_s *cfg);
static int  serial_loop  (ESL_GETOPTS *go, WORKER_INFO *info, CM_FILE *cmfp);
#ifdef SUB_WITH_HMMER_THREADS
#define BLOCK_SIZE 1000

static int  thread_loop(ESL_THREADS *obj, ESL_WORK_QUEUE *queue, CM_FILE *cmfp);
static void pipeline_thread(void *arg);
#endif /*SUB_WITH_HMMER_THREADS*/

#ifdef HAVE_MPI
static int  mpi_master   (ESL_GETOPTS *go, struct cfg_s *cfg);
static int  mpi_worker   (ESL_GETOPTS *go, struct cfg_s *cfg);
#endif /*HAVE_MPI*/

/* process_commandline()
 * 
 * Processes the commandline, filling in fields in <cfg> and creating and returning
 * an <ESL_GETOPTS> options structure. The help page (hmmsearch -h) is formatted
 * here.
 */
static void
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_cmfile, char **ret_seqfile)
{
  ESL_GETOPTS *go = NULL;

  if ((go = esl_getopts_Create(options))     == NULL)     cm_Fail("Internal failure creating options object");
  if (esl_opt_ProcessEnvironment(go)         != eslOK)  { printf("Failed to process environment: %s\n", go->errbuf); goto ERROR; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { printf("Failed to parse command line: %s\n",  go->errbuf); goto ERROR; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { printf("Failed to parse command line: %s\n",  go->errbuf); goto ERROR; }
 
  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      cm_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\nBasic options:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/

      puts("\nOptions controlling output:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 

      puts("\nOptions controlling reporting thresholds:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 

      puts("\nOptions controlling inclusion (significance) thresholds:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80); 

      puts("\nOptions for model-specific thresholding:");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80); 

      puts("\nOptions controlling acceleration heuristics:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80); 

      puts("\nOther expert options:");
      esl_opt_DisplayHelp(stdout, go, 12, 2, 80); 
      exit(0);
    }

  if (esl_opt_ArgNumber(go)                 != 2)      { puts("Incorrect number of command line arguments.");      goto ERROR; }
  if ((*ret_cmfile  = esl_opt_GetArg(go, 1)) == NULL)  { puts("Failed to get <cmdb> argument on command line");   goto ERROR; }
  if ((*ret_seqfile = esl_opt_GetArg(go, 2)) == NULL)  { puts("Failed to get <seqfile> argument on command line"); goto ERROR; }

  /* Validate any attempted use of stdin streams */
  if (strcmp(*ret_cmfile, "-") == 0) {
    puts("cmscan cannot read <cm database> from stdin stream, because it must have cmpress'ed auxfiles");
    goto ERROR;
  }

  *ret_go = go;
  return;
  
 ERROR:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  puts("\nwhere most common options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
  printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
  exit(1);  
}


static int
output_header(FILE *ofp, ESL_GETOPTS *go, char *cmfile, char *seqfile)
{
  cm_banner(ofp, go->argv[0], banner);
  
  fprintf(ofp, "# query sequence file:             %s\n", seqfile);
  fprintf(ofp, "# target CM database:              %s\n", cmfile);
  if (esl_opt_IsUsed(go, "-o"))          fprintf(ofp, "# output directed to file:         %s\n",      esl_opt_GetString(go, "-o"));
  if (esl_opt_IsUsed(go, "--tblout"))    fprintf(ofp, "# per-seq hits tabular output:     %s\n",      esl_opt_GetString(go, "--tblout"));
  if (esl_opt_IsUsed(go, "--acc"))       fprintf(ofp, "# prefer accessions over names:    yes\n");
  if (esl_opt_IsUsed(go, "--noali"))     fprintf(ofp, "# show alignments in output:       no\n");
  if (esl_opt_IsUsed(go, "--notextw"))   fprintf(ofp, "# max ASCII text line length:      unlimited\n");
  if (esl_opt_IsUsed(go, "--textw"))     fprintf(ofp, "# max ASCII text line length:      %d\n",            esl_opt_GetInteger(go, "--textw"));  
  if (esl_opt_IsUsed(go, "-E"))          fprintf(ofp, "# profile reporting threshold:     E-value <= %g\n", esl_opt_GetReal(go, "-E"));
  if (esl_opt_IsUsed(go, "-T"))          fprintf(ofp, "# profile reporting threshold:     score >= %g\n",   esl_opt_GetReal(go, "-T"));
  if (esl_opt_IsUsed(go, "--incE"))      fprintf(ofp, "# profile inclusion threshold:     E-value <= %g\n", esl_opt_GetReal(go, "--incE"));
  if (esl_opt_IsUsed(go, "--incT"))      fprintf(ofp, "# profile inclusion threshold:     score >= %g\n",   esl_opt_GetReal(go, "--incT"));
  if (esl_opt_IsUsed(go, "--cut_ga"))    fprintf(ofp, "# model-specific thresholding:     GA cutoffs\n"); 
  if (esl_opt_IsUsed(go, "--cut_nc"))    fprintf(ofp, "# model-specific thresholding:     NC cutoffs\n"); 
  if (esl_opt_IsUsed(go, "--cut_tc"))    fprintf(ofp, "# model-specific thresholding:     TC cutoffs\n"); 

  if (esl_opt_IsUsed(go, "--cyk"))        fprintf(ofp, "# use CYK for final search stage         on\n");
  if (esl_opt_IsUsed(go, "--fqdb"))       fprintf(ofp, "# QDBs (CYK filter stage)                on\n");
  if (esl_opt_IsUsed(go, "--fsums"))      fprintf(ofp, "# HMM bands from sums (filter)           on\n");
  if (esl_opt_IsUsed(go, "--qdb"))        fprintf(ofp, "# QDBs (final stage)                     on\n");
  if (esl_opt_IsUsed(go, "--nonbanded"))  fprintf(ofp, "# No bands (final stage)                 on\n");
  if (esl_opt_IsUsed(go, "--sums"))       fprintf(ofp, "# HMM bands from sums (final)            on\n");
  if (esl_opt_IsUsed(go, "--max"))        fprintf(ofp, "# Max sensitivity mode:                  on [all heuristic filters off]\n");
  if (esl_opt_IsUsed(go, "--fZ"))         fprintf(ofp, "# Filters set as if DB size in Mb is:    %f\n", esl_opt_GetReal(go, "--fZ"));
  if (esl_opt_IsUsed(go, "--rfam"))       fprintf(ofp, "# Rfam pipeline mode:                    on\n");
  if (esl_opt_IsUsed(go, "--noenvdef"))   fprintf(ofp, "# Envelope definition prior to CM search:off\n");
  if (esl_opt_IsUsed(go, "--pad"))        fprintf(ofp, "# hit padding strategy:                  on\n");
  if (esl_opt_IsUsed(go, "--noF1"))       fprintf(ofp, "# HMM MSV filter:                        off\n");
  if (esl_opt_IsUsed(go, "--noF2"))       fprintf(ofp, "# HMM Vit filter:                        off\n");
  if (esl_opt_IsUsed(go, "--noF3"))       fprintf(ofp, "# HMM Fwd filter:                        off\n");
  if (esl_opt_IsUsed(go, "--noF4"))       fprintf(ofp, "# HMM glocal Fwd filter:                 off\n");
  if (esl_opt_IsUsed(go, "--nohmm"))      fprintf(ofp, "# HMM filters (MSV/bias/Vit/Fwd):        off\n");
  if (esl_opt_IsUsed(go, "--noF6"))       fprintf(ofp, "# CM CYK filter:                         off\n");
  if (esl_opt_IsUsed(go, "--doF1b"))      fprintf(ofp, "# HMM MSV biased comp filter:            on\n");
  if (esl_opt_IsUsed(go, "--noF2b"))      fprintf(ofp, "# HMM Vit biased comp filter:            off\n");
  if (esl_opt_IsUsed(go, "--noF3b"))      fprintf(ofp, "# HMM Fwd biased comp filter:            off\n");
  if (esl_opt_IsUsed(go, "--noF4b"))      fprintf(ofp, "# HMM gFwd biased comp filter:           off\n");
  if (esl_opt_IsUsed(go, "--noF5b"))      fprintf(ofp, "# HMM per-envelope biased comp filter:   off\n");
  if (esl_opt_IsUsed(go, "--F1"))         fprintf(ofp, "# HMM MSV filter P threshold:            <= %g\n", esl_opt_GetReal(go, "--F1"));
  if (esl_opt_IsUsed(go, "--F1b"))        fprintf(ofp, "# HMM MSV bias P threshold:              <= %g\n", esl_opt_GetReal(go, "--F1b"));
  if (esl_opt_IsUsed(go, "--F2"))         fprintf(ofp, "# HMM Vit filter P threshold:            <= %g\n", esl_opt_GetReal(go, "--F2"));
  if (esl_opt_IsUsed(go, "--F2b"))        fprintf(ofp, "# HMM Vit bias P threshold:              <= %g\n", esl_opt_GetReal(go, "--F2b"));
  if (esl_opt_IsUsed(go, "--F3"))         fprintf(ofp, "# HMM Fwd filter P threshold:            <= %g\n", esl_opt_GetReal(go, "--F3"));
  if (esl_opt_IsUsed(go, "--F3b"))        fprintf(ofp, "# HMM Fwd bias P threshold:              <= %g\n", esl_opt_GetReal(go, "--F3b"));
  if (esl_opt_IsUsed(go, "--F4"))         fprintf(ofp, "# HMM glocal Fwd filter P threshold:     <= %g\n", esl_opt_GetReal(go, "--F4"));
  if (esl_opt_IsUsed(go, "--F4b"))        fprintf(ofp, "# HMM glocal Fwd bias P threshold:       <= %g\n", esl_opt_GetReal(go, "--F4b"));
  if (esl_opt_IsUsed(go, "--F5"))         fprintf(ofp, "# HMM env defn filter P threshold:       <= %g\n", esl_opt_GetReal(go, "--F5"));
  if (esl_opt_IsUsed(go, "--F5b"))        fprintf(ofp, "# HMM env defn bias   P threshold:       <= %g\n", esl_opt_GetReal(go, "--F5b"));
  if (esl_opt_IsUsed(go, "--F6"))         fprintf(ofp, "# CM CYK filter P threshold:             <= %g\n", esl_opt_GetReal(go, "--F6"));
  if (esl_opt_IsUsed(go, "--rt1"))        fprintf(ofp, "# Domain definition rt1 parameter        %g\n", esl_opt_GetReal(go, "--rt1"));
  if (esl_opt_IsUsed(go, "--rt2"))        fprintf(ofp, "# Domain definition rt2 parameter        %g\n", esl_opt_GetReal(go, "--rt2"));
  if (esl_opt_IsUsed(go, "--rt3"))        fprintf(ofp, "# Domain definition rt3 parameter        %g\n", esl_opt_GetReal(go, "--rt3"));
  if (esl_opt_IsUsed(go, "--ns"))         fprintf(ofp, "# Number of envelope tracebacks sampled  %d\n", esl_opt_GetInteger(go, "--ns"));
  if (esl_opt_IsUsed(go, "--localenv"))   fprintf(ofp, "# Define envelopes in local mode         on\n");
  if (esl_opt_IsUsed(go, "--null2"))      fprintf(ofp, "# null2 bias corrections:                on\n");
  if (esl_opt_IsUsed(go, "--nonull3"))    fprintf(ofp, "# null3 bias corrections:                off\n");
  if (esl_opt_IsUsed(go, "--toponly"))    fprintf(ofp, "# search top-strand only:                on\n");
  if (esl_opt_IsUsed(go, "--bottomonly")) fprintf(ofp, "# search bottom-strand only:             on\n");

  if (esl_opt_IsUsed(go, "--null2"))      fprintf(ofp, "# null2 bias corrections:                on\n");
  if (esl_opt_IsUsed(go, "-Z"))          fprintf(ofp, "# sequence search space set to:    %.0f\n",    esl_opt_GetReal(go, "-Z"));
  if (esl_opt_IsUsed(go, "--seed"))  {
    if (esl_opt_GetInteger(go, "--seed")==0)fprintf(ofp, "# random number seed:              one-time arbitrary\n");
    else                                    fprintf(ofp, "# random number seed set to:       %d\n", esl_opt_GetInteger(go, "--seed"));
  }
  if (esl_opt_IsUsed(go, "--qformat"))   fprintf(ofp, "# input seqfile format asserted:   %s\n", esl_opt_GetString(go, "--qformat"));
  if (esl_opt_IsUsed(go, "--daemon"))    fprintf(ofp, "run as a daemon process\n");
#ifdef SUB_WITH_HMMER_THREADS
  if (esl_opt_IsUsed(go, "--cpu"))       fprintf(ofp, "# number of worker threads:        %d\n", esl_opt_GetInteger(go, "--cpu"));  
#endif
#ifdef HAVE_MPI
  if (esl_opt_IsUsed(go, "--mpi"))       fprintf(ofp, "# MPI:                             on\n");
#endif
  fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
  return eslOK;
}


int
main(int argc, char **argv)
{
  int              status   = eslOK;

  ESL_GETOPTS     *go  = NULL;	/* command line processing                 */
  struct cfg_s     cfg;         /* configuration data                      */

  /* Set processor specific flags */
  impl_Init();

  /* Initialize what we can in the config structure (without knowing the alphabet yet) 
   */
  cfg.cmfile     = NULL;
  cfg.seqfile    = NULL;

  cfg.do_mpi     = FALSE;	           /* this gets reset below, if we init MPI */
  cfg.nproc      = 0;		           /* this gets reset below, if we init MPI */
  cfg.my_rank    = 0;		           /* this gets reset below, if we init MPI */

  /* Initializations */
  init_ilogsum();
  FLogsumInit();
  p7_FLogsumInit();		/* we're going to use table-driven Logsum() approximations at times */
  process_commandline(argc, argv, &go, &cfg.cmfile, &cfg.seqfile);    

  /* Figure out who we are, and send control there: 
   * we might be an MPI master, an MPI worker, or a serial program.
   */
#ifdef HAVE_MPI
  /* TEMP */
  pid_t pid;
  /* get the process id */
  pid = getpid();
  printf("The process id is %d\n", pid);
  fflush(stdout);
  /* TEMP */

  /* pause the execution of the programs execution until the user has a
   * chance to attach with a debugger and send a signal to resume execution
   * i.e. (gdb) signal SIGCONT
   */
  if (esl_opt_GetBoolean(go, "--stall")) pause();

  if (esl_opt_GetBoolean(go, "--mpi")) 
    {
      cfg.do_mpi     = TRUE;
      MPI_Init(&argc, &argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &(cfg.my_rank));
      MPI_Comm_size(MPI_COMM_WORLD, &(cfg.nproc));

      if (cfg.my_rank > 0)  status = mpi_worker(go, &cfg);
      else 		    status = mpi_master(go, &cfg);

      MPI_Finalize();
    }
  else
#endif /*HAVE_MPI*/
    {
      status = serial_master(go, &cfg);
    }

  esl_getopts_Destroy(go);

  return status;
}


/* serial_master()
 * The serial version of cmscan.
 * For each query sequence in <seqfile> search the database of CMs for hits.
 * 
 * A master can only return if it's successful. All errors are handled immediately and fatally with p7_Fail().
 */
static int
serial_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  FILE            *ofp      = stdout;	         /* output file for results (default stdout)        */
  FILE            *tblfp    = NULL;		 /* output stream for tabular per-seq (--tblout)    */
  int              seqfmt   = eslSQFILE_UNKNOWN; /* format of seqfile                               */
  ESL_SQFILE      *sqfp     = NULL;              /* open seqfile                                    */
  CM_FILE         *cmfp     = NULL;		 /* open CM database file                           */
  ESL_ALPHABET    *abc      = NULL;              /* sequence alphabet                               */
  ESL_STOPWATCH   *w        = NULL;              /* timing                                          */
  ESL_SQ          *qsq      = NULL;		 /* query sequence                                  */
  int              nquery   = 0;
  int              textw;
  int              status   = eslOK;
  int              hstatus  = eslOK;
  int              sstatus  = eslOK;
  int              i;

  int              ncpus    = 0;

  int              infocnt  = 0;
  WORKER_INFO     *info     = NULL;
#ifdef SUB_WITH_HMMER_THREADS
  P7_OM_BLOCK     *block    = NULL;
  ESL_THREADS     *threadObj= NULL;
  ESL_WORK_QUEUE  *queue    = NULL;
#endif
  char             errbuf[eslERRBUFSIZE];

  w = esl_stopwatch_Create();

  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");

  /* If caller declared an input format, decode it */
  if (esl_opt_IsOn(go, "--qformat")) {
    seqfmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--qformat"));
    if (seqfmt == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized input sequence file format\n", esl_opt_GetString(go, "--qformat"));
  }

  /* validate options if running as a daemon */
  if (esl_opt_IsOn(go, "--daemon")) {

    /* running as a daemon, the input format must be type daemon */
    if (seqfmt != eslSQFILE_UNKNOWN && seqfmt != eslSQFILE_DAEMON) 
      esl_fatal("Input format %s not supported.  Must be daemon\n", esl_opt_GetString(go, "--qformat"));
    seqfmt = eslSQFILE_DAEMON;

    if (strcmp(cfg->seqfile, "-") != 0) esl_fatal("Query sequence file must be '-'\n");
  }
  /* Open the target CM database and read 1 CM, but only to get the sequence alphabet */
  status = cm_file_Open(cfg->cmfile, CMDBENV, FALSE, &cmfp, errbuf);
  if      (status == eslENOTFOUND) cm_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", cfg->cmfile, errbuf);
  else if (status == eslEFORMAT)   cm_Fail("File format problem in trying to open HMM file %s.\n%s\n",                cfg->cmfile, errbuf);
  else if (status != eslOK)        cm_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",               status, cfg->cmfile, errbuf);  
  if      (cmfp->do_gzip)          cm_Fail("Reading gzipped CM files is not supported");
  if      (cmfp->do_stdin)         cm_Fail("Reading CM files from stdin is not supported");
  // TODO:cmpress? if (! cmfp->is_pressed)           cm_Fail("Failed to open binary auxfiles for %s: use cmpress first\n",             hfp->fname);

  hstatus = cm_file_Read(cmfp, &abc, NULL);
  if(hstatus == eslEFORMAT)  cm_Fail("bad file format in CM file %s\n%s",           cfg->cmfile, cmfp->errbuf);
  else if (hstatus != eslOK) cm_Fail("Unexpected error in reading CMs from %s\n%s", cfg->cmfile, cmfp->errbuf); 

  /* Determine database size: #CMs in the CM database, currently requires <cmfile>.ssi SSI file (RETHINK IF YOU IMPLEMENT cmpress) */
  if(esl_opt_IsUsed(go, "-Z")) { /* enabling -Z is the only way to bypass the SSI index file requirement */
    cfg->Z       = (int64_t) esl_opt_GetReal(go, "-Z");
    cfg->Z_setby = CM_ZSETBY_OPTION; 
  }
  else { 
    if(cmfp->ssi == NULL) cm_Fail("Failed to open SSI index for CM file: %s\n", cmfp->fname);
    cfg->Z = (int64_t) cmfp->ssi->nprimary;
    cfg->Z_setby = CM_ZSETBY_SSIINFO; 
  }
  cm_file_Close(cmfp);

  /* Open the query sequence database */
  status = esl_sqfile_OpenDigital(abc, cfg->seqfile, seqfmt, NULL, &sqfp);
  if      (status == eslENOTFOUND) cm_Fail("Failed to open sequence file %s for reading\n",      cfg->seqfile);
  else if (status == eslEFORMAT)   cm_Fail("Sequence file %s is empty or misformatted\n",        cfg->seqfile);
  else if (status == eslEINVAL)    cm_Fail("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        cm_Fail("Unexpected error %d opening sequence file %s\n", status, cfg->seqfile);
  qsq = esl_sq_CreateDigital(abc);

  /* Open the results output files */
  if (esl_opt_IsOn(go, "-o"))          { if ((ofp      = fopen(esl_opt_GetString(go, "-o"),          "w")) == NULL)  esl_fatal("Failed to open output file %s for writing\n",                 esl_opt_GetString(go, "-o")); }
  if (esl_opt_IsOn(go, "--tblout"))    { if ((tblfp    = fopen(esl_opt_GetString(go, "--tblout"),    "w")) == NULL)  esl_fatal("Failed to open tabular per-seq output file %s for writing\n", esl_opt_GetString(go, "--tblfp")); }
 
  output_header(ofp, go, cfg->cmfile, cfg->seqfile);

#ifdef SUB_WITH_HMMER_THREADS
  /* initialize thread data */
  if (esl_opt_IsOn(go, "--cpu")) ncpus = esl_opt_GetInteger(go, "--cpu");
  else                           esl_threads_CPUCount(&ncpus);

  if (ncpus > 0)
    {
      threadObj = esl_threads_Create(&pipeline_thread);
      queue = esl_workqueue_Create(ncpus * 2);
    }
#endif

  infocnt = (ncpus == 0) ? 1 : ncpus;
  ESL_ALLOC(info, sizeof(*info) * infocnt);

  for (i = 0; i < infocnt; ++i)
    {
      //info[i].bg    = p7_bg_Create(abc);
#ifdef SUB_WITH_HMMER_THREADS
      info[i].queue = queue;
#endif
    }

#ifdef SUB_WITH_HMMER_THREADS
  for (i = 0; i < ncpus * 2; ++i)
    {
      block = p7_oprofile_CreateBlock(BLOCK_SIZE);
      if (block == NULL)    esl_fatal("Failed to allocate sequence block");

      status = esl_workqueue_Init(queue, block);
      if (status != eslOK)  esl_fatal("Failed to add block to work queue");
    }
#endif

  /* Outside loop: over each query sequence in <seqfile>. */
  while ((sstatus = esl_sqio_Read(sqfp, qsq)) == eslOK)
    {
      nquery++;
      esl_stopwatch_Start(w);	                          

      /* Open the target profile database */
      status = cm_file_Open(cfg->cmfile, CMDBENV, FALSE, &cmfp, NULL);
      if (status != eslOK)        cm_Fail("Unexpected error %d in opening cm file %s.\n%s",           status, cfg->cmfile, cmfp->errbuf);  
  
#ifdef SUB_WITH_HMMER_THREADS
      /* if we are threaded, create a lock to prevent multiple readers */
      if (ncpus > 0)
	{
	  status = cm_file_CreateLock(cmfp);
	  if (status != eslOK) cm_Fail("Unexpected error %d creating lock\n", status);
	}
#endif

      fprintf(ofp, "Query:       %s  [L=%ld]\n", qsq->name, (long) qsq->n);
      if (qsq->acc[0]  != 0) fprintf(ofp, "Accession:   %s\n", qsq->acc);
      if (qsq->desc[0] != 0) fprintf(ofp, "Description: %s\n", qsq->desc);

      for (i = 0; i < infocnt; ++i)
	{
	  /* Create processing pipeline and hit list */
	  info[i].th  = cm_tophits_Create(); 
	  info[i].pli = cm_pipeline_Create(go, 100, 100, cfg->Z, cfg->Z_setby, CM_SCAN_MODELS); /* M_hint = 100, L_hint = 100 are just dummies for now */
	  info[i].pli->cmfp = cmfp;  /* for two-stage input, pipeline needs <hfp> */

	  cm_pli_NewSeq(info[i].pli, qsq, 0); /* 0 is sequence index, only used to remove overlaps */
	  info[i].qsq = qsq;

#ifdef SUB_WITH_HMMER_THREADS
	  if (ncpus > 0) esl_threads_AddThread(threadObj, &info[i]);
#endif
	}

#ifdef SUB_WITH_HMMER_THREADS
      if (ncpus > 0)  hstatus = thread_loop(threadObj, queue, hfp);
      else	      hstatus = serial_loop(go, info, hfp);
#else
      hstatus = serial_loop(go, info, cmfp);
#endif
      switch(hstatus)
	{
	case eslEFORMAT:
	  cm_Fail("bad file format in CM file %s",             cfg->cmfile);
	  break;
	case eslEINCOMPAT:
	  cm_Fail("CM file %s contains different alphabets",   cfg->cmfile);
	  break;
	case eslEMEM:
	  cm_Fail("Out of memory");
	  break;
	case eslEOF:
	  /* do nothing */
	  break;
	default:
	  cm_Fail("Unexpected error in reading CMs from %s",   cfg->cmfile); 
	}

      /* merge the results of the search results */
      for (i = 1; i < infocnt; ++i)
	{
	  cm_tophits_Merge(info[0].th, info[i].th);
	  cm_pipeline_Merge(info[0].pli, info[i].pli);

	  cm_pipeline_Destroy(info[i].pli, NULL);
	  cm_tophits_Destroy(info[i].th);
	}

      /* Print results */
      cm_tophits_SortByScore(info->th);
      cm_tophits_Threshold(info->th, info->pli);

      /* tally up total number of hits and target coverage */
      info->pli->n_output = info->pli->pos_output = 0;
      for (i = 0; i < info->th->N; i++) {
	if ((info->th->hit[i]->flags & CM_HIT_IS_REPORTED) || (info->th->hit[i]->flags & CM_HIT_IS_INCLUDED)) { 
	  info->pli->n_output++;
	  info->pli->pos_output += abs(info->th->hit[i]->stop - info->th->hit[i]->start) + 1;
	}
      }
      cm_tophits_Targets(ofp, info->th, info->pli, textw); fprintf(ofp, "\n\n");
      fprintf(ofp, "\n\n");
      if(info->pli->do_alignments) {
	if((status = cm_tophits_HitAlignments(ofp, info->th, info->pli, textw)) != eslOK) esl_fatal("Out of memory");
	fprintf(ofp, "\n\n");
      }

      if (tblfp)    cm_tophits_TabularTargets(tblfp,    qsq->name, qsq->acc, info->th, info->pli, (nquery == 1));

      esl_stopwatch_Stop(w);
      cm_pli_Statistics(ofp, info->pli, w);
      fprintf(ofp, "//\n"); fflush(ofp);

      cm_file_Close(cmfp);
      cm_pipeline_Destroy(info->pli, NULL);
      cm_tophits_Destroy(info->th);
      esl_sq_Reuse(qsq);
    }
  if      (sstatus == eslEFORMAT) esl_fatal("Parse failed (sequence file %s):\n%s\n",
					    sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
  else if (sstatus != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s",
					    sstatus, sqfp->filename);

  /* Terminate outputs - any last words?
   */
  if (tblfp)    cm_tophits_TabularTail(tblfp,    "cmscan", CM_SCAN_MODELS, cfg->seqfile, cfg->cmfile, go);
  if (ofp)      fprintf(ofp, "[ok]\n");

  /* Cleanup - prepare for successful exit
   */
  for (i = 0; i < infocnt; ++i)
    {
      //p7_bg_Destroy(info[i].bg);
    }

#ifdef SUB_WITH_HMMER_THREADS
  if (ncpus > 0)
    {
      esl_workqueue_Reset(queue);
      while (esl_workqueue_Remove(queue, (void **) &block) == eslOK)
	{
	  p7_oprofile_DestroyBlock(block);
	}
      esl_workqueue_Destroy(queue);
      esl_threads_Destroy(threadObj);
    }
#endif

  free(info);

  esl_sq_Destroy(qsq);
  esl_stopwatch_Destroy(w);
  esl_alphabet_Destroy(abc);
  esl_sqfile_Close(sqfp);

  if (ofp != stdout) fclose(ofp);
  if (tblfp)         fclose(tblfp);

  return eslOK;

 ERROR:
  return eslFAIL;
}

#ifdef HAVE_MPI

/* Define common tags used by the MPI master/slave processes */
#define HMMER_ERROR_TAG          1
#define HMMER_HMM_TAG            2
#define HMMER_SEQUENCE_TAG       3
#define HMMER_BLOCK_TAG          4
#define HMMER_PIPELINE_TAG       5
#define HMMER_TOPHITS_TAG        6
#define HMMER_HIT_TAG            7
#define HMMER_TERMINATING_TAG    8
#define HMMER_READY_TAG          9

/* mpi_failure()
 * Generate an error message.  If the clients rank is not 0, a
 * message is created with the error message and sent to the
 * master process for handling.
 */
static void
mpi_failure(char *format, ...)
{
  va_list  argp;
  int      status = eslFAIL;
  int      len;
  int      rank;
  char     str[512];

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* format the error mesg */
  va_start(argp, format);
  len = vsnprintf(str, sizeof(str), format, argp);
  va_end(argp);

  /* make sure the error string is terminated */
  str[sizeof(str)-1] = '\0';

  /* if the caller is the master, print the results and abort */
  if (rank == 0)
    {
      fprintf(stderr, "\nError: ");
      fprintf(stderr, "%s", str);
      fprintf(stderr, "\n");
      fflush(stderr);

      MPI_Abort(MPI_COMM_WORLD, status);
      exit(1);
    }
  else
    {
      MPI_Send(str, len, MPI_CHAR, 0, HMMER_ERROR_TAG, MPI_COMM_WORLD);
      pause();
    }
}

#define MAX_BLOCK_SIZE (512*1024)

typedef struct {
  uint64_t  offset;
  uint64_t  length;
  uint64_t  count;
} MSV_BLOCK;

typedef struct {
  int        complete;
  int        size;
  int        current;
  int        last;
  MSV_BLOCK *blocks;
} BLOCK_LIST;

/* this routine parses the database keeping track of the blocks
 * offset within the file, number of sequences and the length
 * of the block.  These blocks are passed as work units to the
 * MPI workers.  If multiple hmm's are in the query file, the
 * blocks are reused without parsing the database a second time.
 */
int next_block(CM_FILE *cmfp, BLOCK_LIST *list, MSV_BLOCK *block)
{
  P7_OPROFILE   *om       = NULL;
  ESL_ALPHABET  *abc      = NULL;
  int            status   = eslOK;

  /* if the list has been calculated, use it instead of parsing the database */
  if (list->complete)
    {
      if (list->current == list->last)
	{
	  block->offset = 0;
	  block->length = 0;
	  block->count  = 0;

	  status = eslEOF;
	}
      else
	{
	  int inx = list->current++;

	  block->offset = list->blocks[inx].offset;
	  block->length = list->blocks[inx].length;
	  block->count  = list->blocks[inx].count;

	  status = eslOK;
	}

      return status;
    }

  block->offset = 0;
  block->length = 0;
  block->count = 0;

  while (block->length < MAX_BLOCK_SIZE && (status = p7_oprofile_ReadInfoMSV(hfp, &abc, &om)) == eslOK)
    {
      if (block->count == 0) block->offset = om->roff;
      block->length = om->eoff - block->offset + 1;
      block->count++;
      p7_oprofile_Destroy(om);
    }

  if (status == eslEOF && block->count > 0) status = eslOK;
  if (status == eslEOF) list->complete = 1;

  /* add the block to the list of known blocks */
  if (status == eslOK)
    {
      int inx;

      if (list->last >= list->size)
	{
	  void *tmp;
	  list->size += 500;
	  ESL_RALLOC(list->blocks, tmp, sizeof(MSV_BLOCK) * list->size);
	}

      inx = list->last++;
      list->blocks[inx].offset = block->offset;
      list->blocks[inx].length = block->length;
      list->blocks[inx].count  = block->count;
    }

  return status;

 ERROR:
  return eslEMEM;
}

/* mpi_master()
 * The MPI version of hmmscan
 * Follows standard pattern for a master/worker load-balanced MPI program (J1/78-79).
 * 
 * A master can only return if it's successful. 
 * Errors in an MPI master come in two classes: recoverable and nonrecoverable.
 * 
 * Recoverable errors include all worker-side errors, and any
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
mpi_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  FILE            *ofp      = stdout;	         /* output file for results (default stdout)        */
  FILE            *tblfp    = NULL;		 /* output stream for tabular per-seq (--tblout)    */
  int              seqfmt   = eslSQFILE_UNKNOWN; /* format of seqfile                               */
  P7_BG           *bg       = NULL;	         /* null model                                      */
  ESL_SQFILE      *sqfp     = NULL;              /* open seqfile                                    */
  P7_HMMFILE      *hfp      = NULL;		 /* open HMM database file                          */
  ESL_ALPHABET    *abc      = NULL;              /* sequence alphabet                               */
  P7_OPROFILE     *om       = NULL;		 /* target profile                                  */
  ESL_STOPWATCH   *w        = NULL;              /* timing                                          */
  ESL_SQ          *qsq      = NULL;		 /* query sequence                                  */
  int              nquery   = 0;
  int              textw;
  int              status   = eslOK;
  int              hstatus  = eslOK;
  int              sstatus  = eslOK;
  int              dest;

  char            *mpi_buf  = NULL;              /* buffer used to pack/unpack structures */
  int              mpi_size = 0;                 /* size of the allocated buffer */
  BLOCK_LIST      *list     = NULL;
  MSV_BLOCK        block;

  int              i;
  int              size;
  MPI_Status       mpistatus;
  char             errbuf[eslERRBUFSIZE];

  w = esl_stopwatch_Create();
  
  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");

  /* If caller declared an input format, decode it */
  if (esl_opt_IsOn(go, "--qformat")) {
    seqfmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--qformat"));
    if (seqfmt == eslSQFILE_UNKNOWN) mpi_failure("%s is not a recognized input sequence file format\n", esl_opt_GetString(go, "--qformat"));
  }

  /* Open the target profile database to get the sequence alphabet */
  status = p7_hmmfile_OpenE(cfg->hmmfile, p7_HMMDBENV, &hfp, errbuf);
  if      (status == eslENOTFOUND) mpi_failure("File existence/permissions problem in trying to open HMM file %s.\n%s\n", cfg->hmmfile, errbuf);
  else if (status == eslEFORMAT)   mpi_failure("File format problem in trying to open HMM file %s.\n%s\n",                cfg->hmmfile, errbuf);
  else if (status != eslOK)        mpi_failure("Unexpected error %d in opening HMM file %s.\n%s\n",               status, cfg->hmmfile, errbuf);  
  if (! hfp->is_pressed)           mpi_failure("Failed to open binary dbs for HMM file %s: use hmmpress first\n",         hfp->fname);
  
  hstatus = p7_oprofile_ReadMSV(hfp, &abc, &om);
  if      (hstatus == eslEFORMAT)   mpi_failure("bad file format in HMM file %s",             cfg->hmmfile);
  else if (hstatus == eslEINCOMPAT) mpi_failure("HMM file %s contains different alphabets",   cfg->hmmfile);
  else if (hstatus != eslOK)        mpi_failure("Unexpected error in reading HMMs from %s",   cfg->hmmfile); 

  p7_oprofile_Destroy(om);
  p7_hmmfile_Close(hfp);

  /* Open the query sequence database */
  status = esl_sqfile_OpenDigital(abc, cfg->seqfile, seqfmt, NULL, &sqfp);
  if      (status == eslENOTFOUND) mpi_failure("Failed to open sequence file %s for reading\n",      cfg->seqfile);
  else if (status == eslEFORMAT)   mpi_failure("Sequence file %s is empty or misformatted\n",        cfg->seqfile);
  else if (status == eslEINVAL)    mpi_failure("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        mpi_failure("Unexpected error %d opening sequence file %s\n", status, cfg->seqfile);

  /* Open the results output files */
  if (esl_opt_IsOn(go, "-o")          && (ofp      = fopen(esl_opt_GetString(go, "-o"),          "w")) == NULL)
    mpi_failure("Failed to open output file %s for writing\n",                 esl_opt_GetString(go, "-o"));
  if (esl_opt_IsOn(go, "--tblout")    && (tblfp    = fopen(esl_opt_GetString(go, "--tblout"),    "w")) == NULL)
    mpi_failure("Failed to open tabular per-seq output file %s for writing\n", esl_opt_GetString(go, "--tblfp"));
 
  ESL_ALLOC(list, sizeof(MSV_BLOCK));
  list->complete = 0;
  list->size     = 0;
  list->current  = 0;
  list->last     = 0;
  list->blocks   = NULL;

  output_header(ofp, go, cfg->hmmfile, cfg->seqfile);
  qsq = esl_sq_CreateDigital(abc);
  bg = p7_bg_Create(abc);

  /* Outside loop: over each query sequence in <seqfile>. */
  while ((sstatus = esl_sqio_Read(sqfp, qsq)) == eslOK)
    {
      P7_PIPELINE     *pli     = NULL;		/* processing pipeline                      */
      P7_TOPHITS      *th      = NULL;        	/* top-scoring sequence hits                */

      nquery++;

      esl_stopwatch_Start(w);	                          

      /* seqfile may need to be rewound (multiquery mode) */
      if (nquery > 1) list->current = 0;

      /* Open the target profile database */
      status = p7_hmmfile_OpenE(cfg->hmmfile, p7_HMMDBENV, &hfp, NULL);
      if (status != eslOK) mpi_failure("Unexpected error %d in opening hmm file %s.\n", status, cfg->hmmfile);  
  
      fprintf(ofp, "Query:       %s  [L=%ld]\n", qsq->name, (long) qsq->n);
      if (qsq->acc[0]  != 0) fprintf(ofp, "Accession:   %s\n", qsq->acc);
      if (qsq->desc[0] != 0) fprintf(ofp, "Description: %s\n", qsq->desc);

      /* Create processing pipeline and hit list */
      th  = p7_tophits_Create(); 
      pli = p7_pipeline_Create(go, 100, 100, FALSE, p7_SCAN_MODELS); /* M_hint = 100, L_hint = 100 are just dummies for now */
      pli->hfp = hfp;  /* for two-stage input, pipeline needs <hfp> */

      p7_pli_NewSeq(pli, qsq);

      /* Main loop: */
      while ((hstatus = next_block(hfp, list, &block)) == eslOK)
	{
	  if (MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus) != 0) 
	    mpi_failure("MPI error %d receiving message from %d\n", mpistatus.MPI_SOURCE);

	  MPI_Get_count(&mpistatus, MPI_PACKED, &size);
	  if (mpi_buf == NULL || size > mpi_size) {
	    void *tmp;
	    ESL_RALLOC(mpi_buf, tmp, sizeof(char) * size);
	    mpi_size = size; 
	  }

	  dest = mpistatus.MPI_SOURCE;
	  MPI_Recv(mpi_buf, size, MPI_PACKED, dest, mpistatus.MPI_TAG, MPI_COMM_WORLD, &mpistatus);

	  if (mpistatus.MPI_TAG == HMMER_ERROR_TAG)
	    mpi_failure("MPI client %d raised error:\n%s\n", dest, mpi_buf);
	  if (mpistatus.MPI_TAG != HMMER_READY_TAG)
	    mpi_failure("Unexpected tag %d from %d\n", mpistatus.MPI_TAG, dest);
      
	  MPI_Send(&block, 3, MPI_LONG_LONG_INT, dest, HMMER_BLOCK_TAG, MPI_COMM_WORLD);
	}
      switch(hstatus)
	{
	case eslEFORMAT:
	  mpi_failure("bad file format in HMM file %s",              cfg->hmmfile);
	  break;
	case eslEINCOMPAT:
	  mpi_failure("HMM file %s contains different alphabets",    cfg->hmmfile);
	  break;
	case eslEOF:
	  /* do nothing */
	  break;
	default:
	  mpi_failure("Unexpected error %d in reading HMMs from %s", hstatus, cfg->hmmfile); 
	}

      block.offset = 0;
      block.length = 0;
      block.count  = 0;

      /* wait for all workers to finish up their work blocks */
      for (i = 1; i < cfg->nproc; ++i)
	{
	  if (MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus) != 0) 
	    mpi_failure("MPI error %d receiving message from %d\n", mpistatus.MPI_SOURCE);

	  MPI_Get_count(&mpistatus, MPI_PACKED, &size);
	  if (mpi_buf == NULL || size > mpi_size) {
	    void *tmp;
	    ESL_RALLOC(mpi_buf, tmp, sizeof(char) * size);
	    mpi_size = size; 
	  }

	  dest = mpistatus.MPI_SOURCE;
	  MPI_Recv(mpi_buf, size, MPI_PACKED, dest, mpistatus.MPI_TAG, MPI_COMM_WORLD, &mpistatus);

	  if (mpistatus.MPI_TAG == HMMER_ERROR_TAG)
	    mpi_failure("MPI client %d raised error:\n%s\n", dest, mpi_buf);
	  if (mpistatus.MPI_TAG != HMMER_READY_TAG)
	    mpi_failure("Unexpected tag %d from %d\n", mpistatus.MPI_TAG, dest);
	}

      /* merge the results of the search results */
      for (dest = 1; dest < cfg->nproc; ++dest)
	{
	  P7_PIPELINE     *mpi_pli   = NULL;
	  P7_TOPHITS      *mpi_th    = NULL;

	  /* send an empty block to signal the worker they are done */
	  MPI_Send(&block, 3, MPI_LONG_LONG_INT, dest, HMMER_BLOCK_TAG, MPI_COMM_WORLD);

	  /* wait for the results */
	  if ((status = p7_tophits_MPIRecv(dest, HMMER_TOPHITS_TAG, MPI_COMM_WORLD, &mpi_buf, &mpi_size, &mpi_th)) != eslOK)
	    mpi_failure("Unexpected error %d receiving tophits from %d", status, dest);

	  if ((status = p7_pipeline_MPIRecv(dest, HMMER_PIPELINE_TAG, MPI_COMM_WORLD, &mpi_buf, &mpi_size, go, &mpi_pli)) != eslOK)
	    mpi_failure("Unexpected error %d receiving pipeline from %d", status, dest);

	  p7_tophits_Merge(th, mpi_th);
	  p7_pipeline_Merge(pli, mpi_pli);

	  p7_pipeline_Destroy(mpi_pli, NULL);
	  p7_tophits_Destroy(mpi_th);
	}

      /* Print the results.  */
      p7_tophits_SortBySortkey(th);
      p7_tophits_Threshold(th, pli);
      p7_tophits_Targets(ofp, th, pli, textw); fprintf(ofp, "\n\n");
      p7_tophits_Domains(ofp, th, pli, textw); fprintf(ofp, "\n\n");

      if (tblfp)    p7_tophits_TabularTargets(tblfp,    qsq->name, qsq->acc, th, pli, (nquery == 1));
      if (domtblfp) p7_tophits_TabularDomains(domtblfp, qsq->name, qsq->acc, th, pli, (nquery == 1));

      esl_stopwatch_Stop(w);
      p7_pli_Statistics(ofp, pli, w);
      fprintf(ofp, "//\n");

      p7_hmmfile_Close(hfp);
      p7_pipeline_Destroy(pli, NULL);
      p7_tophits_Destroy(th);
      esl_sq_Reuse(qsq);
    }
  if (sstatus == eslEFORMAT) 
    mpi_failure("Parse failed (sequence file %s):\n%s\n", sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
  else if (sstatus != eslEOF)     
    mpi_failure("Unexpected error %d reading sequence file %s", sstatus, sqfp->filename);

  /* monitor all the workers to make sure they have ended */
  for (i = 1; i < cfg->nproc; ++i)
    {
      if (MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus) != 0) 
	mpi_failure("MPI error %d receiving message from %d\n", mpistatus.MPI_SOURCE);

      MPI_Get_count(&mpistatus, MPI_PACKED, &size);
      if (mpi_buf == NULL || size > mpi_size) {
	void *tmp;
	ESL_RALLOC(mpi_buf, tmp, sizeof(char) * size);
	mpi_size = size; 
      }

      dest = mpistatus.MPI_SOURCE;
      MPI_Recv(mpi_buf, size, MPI_PACKED, dest, mpistatus.MPI_TAG, MPI_COMM_WORLD, &mpistatus);

      if (mpistatus.MPI_TAG == HMMER_ERROR_TAG)
	mpi_failure("MPI client %d raised error:\n%s\n", dest, mpi_buf);
      if (mpistatus.MPI_TAG != HMMER_TERMINATING_TAG)
	mpi_failure("Unexpected tag %d from %d\n", mpistatus.MPI_TAG, dest);
    }

 /* Terminate outputs - any last words?
   */
  if (tblfp)    p7_tophits_TabularTail(tblfp,    "hmmscan", p7_SCAN_MODELS, cfg->seqfile, cfg->hmmfile, go);
  if (domtblfp) p7_tophits_TabularTail(domtblfp, "hmmscan", p7_SCAN_MODELS, cfg->seqfile, cfg->hmmfile, go);
  if (ofp)      fprintf(ofp, "[ok]\n");

  /* Cleanup - prepare for successful exit
   */
  free(list);
  if (mpi_buf != NULL) free(mpi_buf);

  p7_bg_Destroy(bg);

  esl_sq_Destroy(qsq);
  esl_stopwatch_Destroy(w);
  esl_alphabet_Destroy(abc);
  esl_sqfile_Close(sqfp);

  if (ofp != stdout) fclose(ofp);
  if (tblfp)         fclose(tblfp);
  if (domtblfp)      fclose(domtblfp);

  return eslOK;

 ERROR:
  return eslFAIL;
}


static int
mpi_worker(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int              seqfmt   = eslSQFILE_UNKNOWN; /* format of seqfile                               */
  P7_BG           *bg       = NULL;	         /* null model                                      */
  ESL_SQFILE      *sqfp     = NULL;              /* open seqfile                                    */
  P7_HMMFILE      *hfp      = NULL;		 /* open HMM database file                          */
  ESL_ALPHABET    *abc      = NULL;              /* sequence alphabet                               */
  P7_OPROFILE     *om       = NULL;		 /* target profile                                  */
  ESL_STOPWATCH   *w        = NULL;              /* timing                                          */
  ESL_SQ          *qsq      = NULL;		 /* query sequence                                  */
  int              status   = eslOK;
  int              hstatus  = eslOK;
  int              sstatus  = eslOK;

  char            *mpi_buf  = NULL;              /* buffer used to pack/unpack structures */
  int              mpi_size = 0;                 /* size of the allocated buffer */

  MPI_Status       mpistatus;
  char             errbuf[eslERRBUFSIZE];

  w = esl_stopwatch_Create();

  /* Open the target profile database to get the sequence alphabet */
  status = p7_hmmfile_OpenE(cfg->hmmfile, p7_HMMDBENV, &hfp, errbuf);
  if      (status == eslENOTFOUND) mpi_failure("File existence/permissions problem in trying to open HMM file %s.\n%s\n", cfg->hmmfile, errbuf);
  else if (status == eslEFORMAT)   mpi_failure("File format problem in trying to open HMM file %s.\n%s\n",                cfg->hmmfile, errbuf);
  else if (status != eslOK)        mpi_failure("Unexpected error %d in opening HMM file %s.\n%s\n",               status, cfg->hmmfile, errbuf);  
  if (! hfp->is_pressed)           mpi_failure("Failed to open binary dbs for HMM file %s: use hmmpress first\n",         hfp->fname);

  hstatus = p7_oprofile_ReadMSV(hfp, &abc, &om);
  if      (hstatus == eslEFORMAT)   mpi_failure("bad file format in HMM file %s",             cfg->hmmfile);
  else if (hstatus == eslEINCOMPAT) mpi_failure("HMM file %s contains different alphabets",   cfg->hmmfile);
  else if (hstatus != eslOK)        mpi_failure("Unexpected error in reading HMMs from %s",   cfg->hmmfile); 

  p7_oprofile_Destroy(om);
  p7_hmmfile_Close(hfp);

  /* Open the query sequence database */
  status = esl_sqfile_OpenDigital(abc, cfg->seqfile, seqfmt, NULL, &sqfp);
  if      (status == eslENOTFOUND) mpi_failure("Failed to open sequence file %s for reading\n",      cfg->seqfile);
  else if (status == eslEFORMAT)   mpi_failure("Sequence file %s is empty or misformatted\n",        cfg->seqfile);
  else if (status == eslEINVAL)    mpi_failure("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        mpi_failure("Unexpected error %d opening sequence file %s\n", status, cfg->seqfile);

  qsq = esl_sq_CreateDigital(abc);
  bg = p7_bg_Create(abc);

  /* Outside loop: over each query sequence in <seqfile>. */
  while ((sstatus = esl_sqio_Read(sqfp, qsq)) == eslOK)
    {
      P7_PIPELINE     *pli     = NULL;		/* processing pipeline                      */
      P7_TOPHITS      *th      = NULL;        	/* top-scoring sequence hits                */

      MSV_BLOCK        block;

      esl_stopwatch_Start(w);

      status = 0;
      MPI_Send(&status, 1, MPI_INT, 0, HMMER_READY_TAG, MPI_COMM_WORLD);

      /* Open the target profile database */
      status = p7_hmmfile_OpenE(cfg->hmmfile, p7_HMMDBENV, &hfp, NULL);
      if (status != eslOK) mpi_failure("Unexpected error %d in opening hmm file %s.\n", status, cfg->hmmfile);  
  
      /* Create processing pipeline and hit list */
      th  = p7_tophits_Create(); 
      pli = p7_pipeline_Create(go, 100, 100, FALSE, p7_SCAN_MODELS); /* M_hint = 100, L_hint = 100 are just dummies for now */
      pli->hfp = hfp;  /* for two-stage input, pipeline needs <hfp> */

      p7_pli_NewSeq(pli, qsq);

      /* receive a sequence block from the master */
      MPI_Recv(&block, 3, MPI_LONG_LONG_INT, 0, HMMER_BLOCK_TAG, MPI_COMM_WORLD, &mpistatus);
      while (block.count > 0)
	{
	  uint64_t length = 0;
	  uint64_t count  = block.count;

	  hstatus = p7_oprofile_Position(hfp, block.offset);
	  if (hstatus != eslOK) mpi_failure("Cannot position optimized model to %ld\n", block.offset);

	  while (count > 0 && (hstatus = p7_oprofile_ReadMSV(hfp, &abc, &om)) == eslOK)
	    {
	      length = om->eoff - block.offset + 1;

	      p7_pli_NewModel(pli, om, bg);
	      p7_bg_SetLength(bg, qsq->n);
	      p7_oprofile_ReconfigLength(om, qsq->n);
	      
	      p7_Pipeline(pli, om, bg, qsq, th);
	      
	      p7_oprofile_Destroy(om);
	      p7_pipeline_Reuse(pli);

	      --count;
	    }

	  /* check the status of reading the hmm */

	  /* lets do a little bit of sanity checking here to make sure the blocks are the same */
	  if (count > 0)              
	    {
	      switch(hstatus)
		{
		case eslEFORMAT:
		  mpi_failure("bad file format in HMM file %s",              cfg->hmmfile);
		  break;
		case eslEINCOMPAT:
		  mpi_failure("HMM file %s contains different alphabets",    cfg->hmmfile);
		  break;
		case eslOK:
		case eslEOF:
		  mpi_failure("Block count mismatch - expected %ld found %ld at offset %ld\n", block.count, block.count-count, block.offset);
		  break;
		default:
		  mpi_failure("Unexpected error %d in reading HMMs from %s", hstatus, cfg->hmmfile); 
		}
	    }
	  if (block.length != length) 
	    mpi_failure("Block length mismatch - expected %ld found %ld at offset %ld\n", block.length, length, block.offset);

	  /* inform the master we need another block of sequences */
	  status = 0;
	  MPI_Send(&status, 1, MPI_INT, 0, HMMER_READY_TAG, MPI_COMM_WORLD);

	  /* wait for the next block of sequences */
	  MPI_Recv(&block, 3, MPI_LONG_LONG_INT, 0, HMMER_BLOCK_TAG, MPI_COMM_WORLD, &mpistatus);
	}

      esl_stopwatch_Stop(w);

      /* Send the top hits back to the master. */
      p7_tophits_MPISend(th, 0, HMMER_TOPHITS_TAG, MPI_COMM_WORLD,  &mpi_buf, &mpi_size);
      p7_pipeline_MPISend(pli, 0, HMMER_PIPELINE_TAG, MPI_COMM_WORLD,  &mpi_buf, &mpi_size);

      p7_hmmfile_Close(hfp);
      p7_pipeline_Destroy(pli, NULL);
      p7_tophits_Destroy(th);
      esl_sq_Reuse(qsq);
    } /* end outer loop over query HMMs */
  if (sstatus == eslEFORMAT) 
    mpi_failure("Parse failed (sequence file %s):\n%s\n", sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
  else if (sstatus != eslEOF)
    mpi_failure("Unexpected error %d reading sequence file %s", sstatus, sqfp->filename);

  status = 0;
  MPI_Send(&status, 1, MPI_INT, 0, HMMER_TERMINATING_TAG, MPI_COMM_WORLD);

  if (mpi_buf != NULL) free(mpi_buf);

  p7_bg_Destroy(bg);

  esl_sq_Destroy(qsq);
  esl_stopwatch_Destroy(w);
  esl_alphabet_Destroy(abc);
  esl_sqfile_Close(sqfp);

  return eslOK;
}
#endif /*HAVE_MPI*/

static int
serial_loop(ESL_GETOPTS *go, WORKER_INFO *info, CM_FILE *cmfp)
{
  int            status;

  CM_t *cm;
  ESL_ALPHABET  *abc = NULL;
  P7_BG           **bgA = NULL;         /* null models                              */
  P7_OPROFILE     **omA = NULL;         /* optimized query profile HMMs            */
  P7_PROFILE      **gmA = NULL;         /* generic   query profile HMMs            */
  P7_HMM **hmmA = NULL;
  float           **p7_evparamAA = NULL;/* [0..nhmm-1][0..CM_p7_NEVPARAM] E-value parameters */
  int      m, z, z_offset;
  int      use_mlp7_hmm = FALSE;
  int             nhmm = 0;
  CMConsensus_t  *cmcons = NULL;
  ESL_STOPWATCH *w = NULL;

  w = esl_stopwatch_Create();

  /* Main loop: */
  while ((status = cm_file_Read(cmfp, &abc, &cm)) == eslOK)
    {
      esl_stopwatch_Start(w);

      /***************************************************/
      /* TEMP: eventually encapsulate this in a function */

      /* Convert HMMs to optimized models.
       * First, determine how many and which HMMs we'll be using *
       * by default, we filter only with the <nhmm> (== cm->nap7) HMMs written in the CM file */
      /* but, if --doml selected or there's no HMMs in the file, use a ML P7 HMM built from the CM */
      if (esl_opt_GetBoolean(go, "--doml") || (cm->nap7 == 0)) { 
	use_mlp7_hmm = TRUE;
	nhmm = 1 + cm->nap7;
      }
      else if(esl_opt_GetBoolean(go, "--noadd")) { 
	/* or if --noadd is selected, only use a ML p7 HMM */
	use_mlp7_hmm = TRUE;
	nhmm = 1;
      }
      else { 
	nhmm = cm->nap7;
      }

      /* Now we know how many HMMs, allocate and fill the necessary data structures for each */
      ESL_ALLOC(hmmA, sizeof(P7_HMM *) * nhmm);
      /* use the ML p7 HMM if nec */
      if (esl_opt_GetBoolean(go, "--doml") || (cm->nap7 == 0)) { 
	hmmA[0] = cm->mlp7;
      }
      /* copy the HMMs from the file into the CM data structure */
      if(! esl_opt_GetBoolean(go, "--noadd")) { 
	z_offset = use_mlp7_hmm ? 1 : 0;
	for(z = 0; z < cm->nap7; z++) { 
	  hmmA[z+z_offset] = cm->ap7A[z];
	}
      }
  
      /* allocate space for profiles, bg and EV params for each HMM */
      ESL_ALLOC(gmA, sizeof(P7_PROFILE *)     * nhmm);
      ESL_ALLOC(omA, sizeof(P7_OPROFILE *)    * nhmm);
      ESL_ALLOC(bgA, sizeof(P7_BG *)          * nhmm);
      ESL_ALLOC(p7_evparamAA, sizeof(float *) * nhmm);
      for(m = 0; m < nhmm; m++) { 
	gmA[m] = p7_profile_Create (hmmA[m]->M, abc);
	omA[m] = p7_oprofile_Create(hmmA[m]->M, abc);
	bgA[m] = p7_bg_Create(abc);
	p7_ProfileConfig(hmmA[m], bgA[m], gmA[m], 100, p7_LOCAL);  /* 100 is a dummy length for now; and MSVFilter requires local mode */
	p7_oprofile_Convert(gmA[m], omA[m]);                       /* <om> is now p7_LOCAL, multihit */
	/* after omA[m]'s been created, convert gmA[m] to glocal, to define envelopes in cm_pipeline() */
	p7_ProfileConfig(hmmA[m], bgA[m], gmA[m], 100, p7_GLOCAL);
	/* copy evalue parameters */
	ESL_ALLOC(p7_evparamAA[m], sizeof(float) * CM_p7_NEVPARAM);
	if(m == 0 && use_mlp7_hmm) { esl_vec_FCopy(cm->mlp7_evparam,       CM_p7_NEVPARAM, p7_evparamAA[m]); }
	else if(use_mlp7_hmm)      { esl_vec_FCopy(cm->ap7_evparamAA[m-1], CM_p7_NEVPARAM, p7_evparamAA[m]); }
	else                       { esl_vec_FCopy(cm->ap7_evparamAA[m],   CM_p7_NEVPARAM, p7_evparamAA[m]); }
      }
      /* END OF TEMP: eventually of block that goes into a function */

      /* BEG OF TEMP: eventually encapsulate this in a function */
      /* Configure the CM */
      if(! esl_opt_GetBoolean(go, "-g")) { 
	if(! esl_opt_GetBoolean(go, "--cp9gloc")) { 
	  cm->config_opts |= CM_CONFIG_LOCAL;
	  cm->config_opts |= CM_CONFIG_HMMLOCAL;
	  if(! esl_opt_GetBoolean(go, "--cp9noel")) cm->config_opts |= CM_CONFIG_HMMEL; 
	}
      }
#if 0
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, "TIMING HMM filter setup");
      esl_stopwatch_Start(w);
      if((status = ConfigCM(cm, errbuf, FALSE, NULL, NULL)) != eslOK) cm_Fail(errbuf); 
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, "TIMING ConfigCM()");
      esl_stopwatch_Start(w);
      if((status = CreateCMConsensus(cm, abc, 3.0, 1.0, &(cmcons)))!= eslOK) cm_Fail("Failed to create CMConsensus data structure.\n");
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, "TIMING CreateCMConsensus()");
      esl_stopwatch_Start(w);
      /* END OF TEMP: end of block */
#endif      
      
      cmcons = NULL;
      if((status = cm_pli_NewModel(info->pli, cm, 
				   FALSE, FALSE,           /* need_fsmx, need_smx */
				   NULL, NULL, NULL, NULL, /* fcyk_dmin, fcyk_dmax, final_dmin, final_dmax */
				   omA, bgA, nhmm)) != eslOK) cm_Fail(info->pli->errbuf);
      
      if((status = cm_Pipeline(info->pli, cm, omA, gmA, bgA, p7_evparamAA, nhmm, info->qsq, info->th, &(cmcons))) != eslOK)
	cm_Fail("cm_pipeline() failed unexpected with status code %d\n%s", status, info->pli->errbuf);
      
      esl_stopwatch_Stop(w);
      esl_stopwatch_Display(stdout, w, "TIMING cm_Pipeline()");
      cm_pipeline_Reuse(info->pli);
      
      if(gmA != NULL)          { for(m = 0; m < nhmm; m++) { if(gmA[m] != NULL)          p7_profile_Destroy(gmA[m]);  } free(gmA); }
      if(omA != NULL)          { for(m = 0; m < nhmm; m++) { if(omA[m] != NULL)          p7_oprofile_Destroy(omA[m]); } free(omA); }
      if(bgA != NULL)          { for(m = 0; m < nhmm; m++) { if(bgA[m] != NULL)          p7_bg_Destroy(bgA[m]);       } free(bgA); }
      if(p7_evparamAA != NULL) { for(m = 0; m < nhmm; m++) { if(p7_evparamAA[m] != NULL) free(p7_evparamAA[m]);       } free(p7_evparamAA); }
      if(cm != NULL)     FreeCM(cm);
      if(cmcons != NULL) FreeCMConsensus(cmcons);
    } /* end of while(cm_file_Read() == eslOK) */

  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);

  return status;

 ERROR:
  if(gmA != NULL)          { for(m = 0; m < nhmm; m++) { if(gmA[m] != NULL)          p7_profile_Destroy(gmA[m]);  } free(gmA); }
  if(omA != NULL)          { for(m = 0; m < nhmm; m++) { if(omA[m] != NULL)          p7_oprofile_Destroy(omA[m]); } free(omA); }
  if(bgA != NULL)          { for(m = 0; m < nhmm; m++) { if(bgA[m] != NULL)          p7_bg_Destroy(bgA[m]);       } free(bgA); }
  if(p7_evparamAA != NULL) { for(m = 0; m < nhmm; m++) { if(p7_evparamAA[m] != NULL) free(p7_evparamAA[m]);       } free(p7_evparamAA); }
  if(cm != NULL)     FreeCM(cm);
  if(cmcons != NULL) FreeCMConsensus(cmcons);
  return status;
}

#ifdef SUB_WITH_HMMER_THREADS
static int
thread_loop(ESL_THREADS *obj, ESL_WORK_QUEUE *queue, P7_HMMFILE *hfp)
{
  int  status   = eslOK;
  int  sstatus  = eslOK;
  int  eofCount = 0;
  P7_OM_BLOCK   *block;
  ESL_ALPHABET  *abc = NULL;
  void          *newBlock;

  esl_workqueue_Reset(queue);
  esl_threads_WaitForStart(obj);

  status = esl_workqueue_ReaderUpdate(queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue reader failed");
      
  /* Main loop: */
  while (sstatus == eslOK)
    {
      block = (P7_OM_BLOCK *) newBlock;
      sstatus = p7_oprofile_ReadBlockMSV(hfp, &abc, block);
      if (sstatus == eslEOF)
	{
	  if (eofCount < esl_threads_GetWorkerCount(obj)) sstatus = eslOK;
	  ++eofCount;
	}
	  
      if (sstatus == eslOK)
	{
	  status = esl_workqueue_ReaderUpdate(queue, block, &newBlock);
	  if (status != eslOK) esl_fatal("Work queue reader failed");
	}
    }

  status = esl_workqueue_ReaderUpdate(queue, block, NULL);
  if (status != eslOK) esl_fatal("Work queue reader failed");

  if (sstatus == eslEOF)
    {
      /* wait for all the threads to complete */
      esl_threads_WaitForFinish(obj);
      esl_workqueue_Complete(queue);  
    }
  
  esl_alphabet_Destroy(abc);
  return sstatus;
}

static void 
pipeline_thread(void *arg)
{
  int i;
  int status;
  int workeridx;
  WORKER_INFO   *info;
  ESL_THREADS   *obj;
  P7_OM_BLOCK   *block;
  void          *newBlock;
  
#ifdef HAVE_FLUSH_ZERO_MODE
  /* In order to avoid the performance penalty dealing with sub-normal
   * values in the floating point calculations, set the processor flag
   * so sub-normals are "flushed" immediately to zero.
   * On OS X, need to reset this flag for each thread
   * (see TW notes 05/08/10 for details)
   */
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif

  obj = (ESL_THREADS *) arg;
  esl_threads_Started(obj, &workeridx);

  info = (WORKER_INFO *) esl_threads_GetData(obj, workeridx);

  status = esl_workqueue_WorkerUpdate(info->queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue worker failed");

  /* loop until all blocks have been processed */
  block = (P7_OM_BLOCK *) newBlock;
  while (block->count > 0)
    {
      /* Main loop: */
      for (i = 0; i < block->count; ++i)
	{
	  P7_OPROFILE *om = block->list[i];

	  p7_pli_NewModel(info->pli, om, info->bg);
	  p7_bg_SetLength(info->bg, info->qsq->n);
	  p7_oprofile_ReconfigLength(om, info->qsq->n);

	  p7_Pipeline(info->pli, om, info->bg, info->qsq, info->th);

	  p7_oprofile_Destroy(om);
	  p7_pipeline_Reuse(info->pli);

	  block->list[i] = NULL;
	}

      status = esl_workqueue_WorkerUpdate(info->queue, block, &newBlock);
      if (status != eslOK) esl_fatal("Work queue worker failed");

      block = (P7_OM_BLOCK *) newBlock;
    }

  status = esl_workqueue_WorkerUpdate(info->queue, block, NULL);
  if (status != eslOK) esl_fatal("Work queue worker failed");

  esl_threads_Finished(obj, workeridx);
  return;
}
#endif   /* SUB_WITH_HMMER_THREADS */


/*****************************************************************
 * @LICENSE@
 *****************************************************************/


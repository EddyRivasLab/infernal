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

#ifdef HMMER_THREADS
#include <unistd.h>
#include "esl_threads.h"
#include "esl_workqueue.h"
#endif /*HMMER_THREADS*/

#include "hmmer.h"
#include "funcs.h"		/* external functions                   */
#include "structs.h"		/* data structures, macros, #define's   */

typedef struct {
#ifdef HMMER_THREADS
  ESL_WORK_QUEUE   *queue;
#endif /*HMMER_THREADS*/
  ESL_SQ           *qsq;
  P7_BG            *bg;	             /* null model                            */
  CM_PIPELINE      *pli;             /* work pipeline                         */
  CM_TOPHITS       *th;              /* top hit results                       */
  float            *p7_evparam;      /* E-value parameters for p7 filter      */
  int               cm_config_opts;  /* set based on command-line options     */
  int               in_rc;           /* TRUE if qsq is currently revcomp'ed   */
} WORKER_INFO;

#define REPOPTS         "-E,-T,--cut_ga,--cut_nc,--cut_tc"
#define INCOPTS         "--incE,--incT,--cut_ga,--cut_nc,--cut_tc"
#define THRESHOPTS      "-E,-T,--incE,--incT,--cut_ga,--cut_nc,--cut_tc"
#define XFASTMIDMAXOPTS "--max,--F1,--F2,--F3,--F4,--F5,--F6,--noF1,--noF2,--noF3,--nogfwd,--nocyk,--nohmm,--hmm"
#define TIMINGOPTS      "--time-F1,--time-F2,--time-F3,--time-dF3,--time-bfil,--time-F4"

#if defined (HMMER_THREADS) && defined (HAVE_MPI)
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
  { "--filcmW",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,             "use CM's window length for all HMM filters", 7 },
  { "--glen",       eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, NULL,              "use length dependent glocal p7 filter P-value thresholds", 7},
  { "--glN",        eslARG_INT,   "201",  NULL, NULL,    NULL,"--glen",NULL,             "minimum value to start len-dependent glocal threshold", 7},
  { "--glX",        eslARG_INT,   "500",  NULL, NULL,    NULL,"--glen",NULL,             "maximum value for len-dependent glocal threshold", 7},
  { "--glstep",     eslARG_INT,   "100",  NULL, NULL,    NULL,"--glen",NULL,             "for len-dependent glocal thr, step size for halving thr", 7},
  { "--noends",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,             "don't search for local envelopes in first/final cm->W residues", 7},
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
  { "--aln-scanbands",eslARG_NONE, FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "use HMM bands from final search stage for alignment of hits, don't recalc", 8},
  { "--aln-newbands", eslARG_NONE, FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "recalculate HMM bands for alignment of hits, don't use scan bands", 8},
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
  { "--tau",          eslARG_REAL,   "5e-6",     NULL, "0<x<1", NULL,        NULL,       "--qdb,--nonbanded,--hmm", "set tail loss prob for --hbanded to <x>", 20 },
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
#ifdef HMMER_THREADS
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
static int  serial_loop  (WORKER_INFO *info, CM_FILE *cmfp);
#if 0 
static int  setup_cm     (CM_t *cm, const ESL_ALPHABET *abc, WORKER_INFO *info, char *errbuf, 
			  int *opt_nhmm, P7_HMM ***opt_hmmA, P7_BG ***opt_bgA, P7_OPROFILE ***opt_omA, 
			  P7_PROFILE ***opt_gmA, float ***opt_p7evpAA);
#endif
static int  determine_config_opts(const ESL_GETOPTS *go);
static void copy_subseq(const ESL_SQ *src_sq, ESL_SQ *dest_sq, int64_t i, int64_t L, int in_rc);

#ifdef HMMER_THREADS
#define BLOCK_SIZE 25

static int  thread_loop(ESL_THREADS *obj, ESL_WORK_QUEUE *queue, CM_FILE *cmfp);
static void pipeline_thread(void *arg);
#endif /*HMMER_THREADS*/

#ifdef HAVE_MPI
static int  mpi_master   (ESL_GETOPTS *go, struct cfg_s *cfg);
static int  mpi_worker   (ESL_GETOPTS *go, struct cfg_s *cfg);
#endif /*HAVE_MPI*/

/* process_commandline()
 * 
 * Processes the commandline, filling in fields in <cfg> and creating and returning
 * an <ESL_GETOPTS> options structure. The help page (cmscan -h) is formatted
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
#ifdef HMMER_THREADS
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
  ESL_STOPWATCH   *mw       = NULL;              /* timing                                          */
  ESL_SQ          *qsq      = NULL;		 /* query sequence                                  */
  int              nquery   = 0;
  int              textw;
  int              status   = eslOK;
  int              hstatus  = eslOK;
  int              sstatus  = eslOK;
  int              i;
  int              in_rc; 

  int              ncpus    = 0;

  int              infocnt  = 0;
  WORKER_INFO     *info     = NULL;
#ifdef HMMER_THREADS
  CM_P7_OM_BLOCK  *block    = NULL;
  ESL_THREADS     *threadObj= NULL;
  ESL_WORK_QUEUE  *queue    = NULL;
#endif
  char             errbuf[eslERRBUFSIZE];

  w = esl_stopwatch_Create();
  mw = esl_stopwatch_Create();

  esl_stopwatch_Start(mw);

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
  if      (status == eslENOTFOUND) cm_Fail("File existence/permissions problem in trying to open CM file %s.\n%s\n", cfg->cmfile, errbuf);
  else if (status == eslEFORMAT)   cm_Fail("File format problem in trying to open CM file %s.\n%s\n",                cfg->cmfile, errbuf);
  else if (status != eslOK)        cm_Fail("Unexpected error %d in opening CM file %s.\n%s\n",               status, cfg->cmfile, errbuf);  
  if      (cmfp->do_gzip)          cm_Fail("Reading gzipped CM files is not supported");
  if      (cmfp->do_stdin)         cm_Fail("Reading CM files from stdin is not supported");
  if (! cmfp->is_pressed)          cm_Fail("Failed to open binary auxfiles for %s: use cmpress first\n",             cmfp->fname);

  hstatus = cm_file_Read(cmfp, FALSE, &abc, NULL);
  if(hstatus == eslEFORMAT)  cm_Fail("bad file format in CM file %s\n%s",           cfg->cmfile, cmfp->errbuf);
  else if (hstatus != eslOK) cm_Fail("Unexpected error in reading CMs from %s\n%s", cfg->cmfile, cmfp->errbuf); 

  /* Determine database size: default is to updated as we read target CMs */
  if(esl_opt_IsUsed(go, "-Z")) { 
    cfg->Z       = (int64_t) esl_opt_GetReal(go, "-Z");
    cfg->Z_setby = CM_ZSETBY_OPTION; 
  }
  else { 
    if(cmfp->ssi == NULL) cm_Fail("Failed to open SSI index for CM file: %s\n", cmfp->fname);
    cfg->Z = (int64_t) cmfp->ssi->nprimary;
    cfg->Z_setby = CM_ZSETBY_SSI_AND_QLENGTH; /* we will multiply Z by each query sequence length */
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

#ifdef HMMER_THREADS
  /* initialize thread data */
  if (esl_opt_IsOn(go, "--cpu")) ncpus = esl_opt_GetInteger(go, "--cpu");
  else                           esl_threads_CPUCount(&ncpus);
  printf("NCPUS: %d\n", ncpus);
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
      info[i].bg    = p7_bg_Create(abc);
      info[i].in_rc = FALSE;
      ESL_ALLOC(info[i].p7_evparam, sizeof(float) * CM_p7_NEVPARAM);
      info[i].cm_config_opts = determine_config_opts(go);
#ifdef HMMER_THREADS
      info[i].queue = queue;
#endif
    }

#ifdef HMMER_THREADS
  for (i = 0; i < ncpus * 2; ++i)
    {
      block = cm_p7_oprofile_CreateBlock(BLOCK_SIZE);
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

      fprintf(ofp, "Query:       %s  [L=%ld]\n", qsq->name, (long) qsq->n);
      if (qsq->acc[0]  != 0) fprintf(ofp, "Accession:   %s\n", qsq->acc);
      if (qsq->desc[0] != 0) fprintf(ofp, "Description: %s\n", qsq->desc);

      for (i = 0; i < infocnt; ++i)
	{
	  /* Create processing pipeline and hit list */
	  info[i].th  = cm_tophits_Create(); 
	  info[i].pli = cm_pipeline_Create(go, abc, 100, 100, cfg->Z * qsq->n, cfg->Z_setby, CM_SCAN_MODELS); /* M_hint = 100, L_hint = 100 are just dummies for now */
	  info[i].pli->nseqs++;
	  info[i].qsq = qsq;
	}

      /* scan all target CMs twice, once with the top strand of the query and once with the bottom strand */
      for(in_rc = 0; in_rc <= 1; in_rc++) { 
	if(in_rc == 0 && (! info->pli->do_top)) continue; /* skip top strand */
	if(in_rc == 1 && (! info->pli->do_bot)) continue; /* skip bottom strand */
	if(in_rc == 1) { 
	  if(qsq->abc->complement == NULL) { 
	    if(! info->pli->do_top) cm_Fail("Trying to only search bottom strand but sequence cannot be complemented"); 
	    else continue; /* skip bottom strand b/c we can't complement the sequence */
	  }
	  esl_sq_ReverseComplement(qsq);
	}

	/* open the target profile database */
	if ((status = cm_file_Open(cfg->cmfile, CMDBENV, FALSE, &cmfp, NULL)) != eslOK) cm_Fail("Unexpected error %d in opening cm file %s.\n%s", status, cfg->cmfile, cmfp->errbuf);  
#ifdef HMMER_THREADS
	if (ncpus > 0) { /* if we are threaded, create locks to prevent multiple readers */
	  if ((status = cm_file_CreateLock(cmfp)) != eslOK) cm_Fail("Unexpected error %d creating lock\n", status);
	}
#endif
	for (i = 0; i < infocnt; ++i) {
	  info[i].pli->cmfp = cmfp;           /* for four-stage input, pipeline needs <cmfp> */
	  cm_pli_NewSeq(info[i].pli, qsq, 0); /* 0 is sequence index, it's irrelevant in this context (only used to remove overlaps) */
	  info[i].in_rc = in_rc;
#ifdef HMMER_THREADS
	if (ncpus > 0) esl_threads_AddThread(threadObj, &info[i]);
#endif
	}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
#ifdef HMMER_THREADS
	if (ncpus > 0)  hstatus = thread_loop(threadObj, queue, cmfp);
	else	        hstatus = serial_loop(info, cmfp);
#else
	hstatus = serial_loop(info, cmfp);
#endif
	switch(hstatus) 
	  { 
	  case eslEFORMAT:   cm_Fail("bad file format in CM file %s",             cfg->cmfile);  break;
	  case eslEINCOMPAT: cm_Fail("CM file %s contains different alphabets",   cfg->cmfile);  break;
	  case eslEMEM:      cm_Fail("Out of memory");	                                         break;
	  case eslEOF:       /* do nothing */                                                    break;
	  default:           cm_Fail("Unexpected error in reading CMs from %s",   cfg->cmfile);  break; 
	  }

	cm_file_Close(cmfp);
      } /* end of 'for(in_rc == 0; in_rc <= 1; in_rc++)' */

      /***********************************************************************/
      /* merge the results of the search results */
      for (i = 1; i < infocnt; ++i)
	{
	  cm_tophits_Merge(info[0].th, info[i].th);
	  cm_pipeline_Merge(info[0].pli, info[i].pli);

	  cm_pipeline_Destroy(info[i].pli, NULL);
	  cm_tophits_Destroy(info[i].th);
	}

      if(info[0].pli->do_top && info[0].pli->do_bot) { 
	/* A brutish hack: we've searched all models versus each
	 * sequence then reverse complement it and search all models
	 * versus it again, so we think we've double counted all
	 * models
	 */
	info[0].pli->nmodels /= 2;
	info[0].pli->nnodes  /= 2;
      }

      if(info[0].pli->research_ends) { 
	/* We may have overlaps so sort by sequence index/position and remove duplicates */
	cm_tophits_SortBySeqIdx(info[0].th);
	cm_tophits_RemoveDuplicates(info[0].th);
      }

      /* Sort by score and enforce threshold. */
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

  esl_stopwatch_Stop(mw);
  esl_stopwatch_Display(stdout, mw, "Total runtime:");

  /* Cleanup - prepare for successful exit
   */
  for (i = 0; i < infocnt; ++i)
    {
      free(info[i].p7_evparam);
      p7_bg_Destroy(info[i].bg);
    }

#ifdef HMMER_THREADS
  if (ncpus > 0)
    {
      esl_workqueue_Reset(queue);
      while (esl_workqueue_Remove(queue, (void **) &block) == eslOK)
	{
	  cm_p7_oprofile_DestroyBlock(block);
	}
      esl_workqueue_Destroy(queue);
      esl_threads_Destroy(threadObj);
    }
#endif

  free(info);

  esl_sq_Destroy(qsq);
  esl_stopwatch_Destroy(w);
  esl_stopwatch_Destroy(mw);
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
#define INFERNAL_ERROR_TAG          1
#define INFERNAL_HMM_TAG            2
#define INFERNAL_SEQUENCE_TAG       3
#define INFERNAL_BLOCK_TAG          4
#define INFERNAL_PIPELINE_TAG       5
#define INFERNAL_TOPHITS_TAG        6
#define INFERNAL_HIT_TAG            7
#define INFERNAL_TERMINATING_TAG    8
#define INFERNAL_READY_TAG          9

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
      MPI_Send(str, len, MPI_CHAR, 0, INFERNAL_ERROR_TAG, MPI_COMM_WORLD);
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

/* This routine parses the database keeping track of the blocks
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
  off_t          prv_offset = 0;

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
  if((prv_offset = ftello(cmfp->ffp)) < 0) return eslESYS;

  while (block->length < MAX_BLOCK_SIZE && 
	 (status = cm_p7_oprofile_ReadMSV(cmfp, FALSE, &abc, NULL, NULL, NULL, NULL, NULL, &om)) == eslOK)
    {
      if (block->count == 0) block->offset = prv_offset;
      if((prv_offset = ftello(cmfp->ffp)) < 0) return eslESYS;
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
 * The MPI version of cmscan
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
  CM_FILE         *cmfp     = NULL;		 /* open CMM database file                          */
  ESL_ALPHABET    *abc      = NULL;              /* sequence alphabet                               */
  ESL_STOPWATCH   *w        = NULL;              /* timing                                          */
  ESL_STOPWATCH   *mw       = NULL;              /* timing                                          */
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
  mw = esl_stopwatch_Create();

  esl_stopwatch_Start(mw);
  
  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");

  /* If caller declared an input format, decode it */
  if (esl_opt_IsOn(go, "--qformat")) {
    seqfmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--qformat"));
    if (seqfmt == eslSQFILE_UNKNOWN) mpi_failure("%s is not a recognized input sequence file format\n", esl_opt_GetString(go, "--qformat"));
  }

  /* Open the target profile database to get the sequence alphabet */
  status = cm_file_Open(cfg->cmfile, CMDBENV, FALSE, &cmfp, errbuf);
  if      (status == eslENOTFOUND) mpi_failure("File existence/permissions problem in trying to open CM file %s.\n%s\n", cfg->cmfile, errbuf);
  else if (status == eslEFORMAT)   mpi_failure("File format problem in trying to open CM file %s.\n%s\n",                cfg->cmfile, errbuf);
  else if (status != eslOK)        mpi_failure("Unexpected error %d in opening CM file %s.\n%s\n",               status, cfg->cmfile, errbuf);  
  if      (cmfp->do_gzip)          mpi_failure("Reading gzipped CM files is not supported");
  if      (cmfp->do_stdin)         mpi_failure("Reading CM files from stdin is not supported");
  if (! cmfp->is_pressed)          mpi_failure("Failed to open binary auxfiles for %s: use cmpress first\n",             cmfp->fname);

  hstatus = cm_file_Read(cmfp, FALSE, &abc, NULL);
  if(hstatus == eslEFORMAT)  mpi_failure("bad file format in CM file %s\n%s",           cfg->cmfile, cmfp->errbuf);
  else if (hstatus != eslOK) mpi_failure("Unexpected error in reading CMs from %s\n%s", cfg->cmfile, cmfp->errbuf); 

  /* Determine database size: default is to updated as we read target CMs */
  if(esl_opt_IsUsed(go, "-Z")) { 
    cfg->Z       = (int64_t) esl_opt_GetReal(go, "-Z");
    cfg->Z_setby = CM_ZSETBY_OPTION; 
  }
  else { 
    if(cmfp->ssi == NULL) mpi_failure("Failed to open SSI index for CM file: %s\n", cmfp->fname);
    cfg->Z = (int64_t) cmfp->ssi->nprimary;
    cfg->Z_setby = CM_ZSETBY_SSI_AND_QLENGTH; /* we will multiply Z by each query sequence length */
  }
  cm_file_Close(cmfp);

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
    mpi_failure("Failed to open tabular per-seq output file %s for writing\n", esl_opt_GetString(go, "--tblout"));
 
  ESL_ALLOC(list, sizeof(MSV_BLOCK));
  list->complete = 0;
  list->size     = 0;
  list->current  = 0;
  list->last     = 0;
  list->blocks   = NULL;

  output_header(ofp, go, cfg->cmfile, cfg->seqfile);
  qsq = esl_sq_CreateDigital(abc);
  bg = p7_bg_Create(abc);

  /* Outside loop: over each query sequence in <seqfile>. */
  while ((sstatus = esl_sqio_Read(sqfp, qsq)) == eslOK)
    {
      CM_PIPELINE     *pli     = NULL;		/* processing pipeline                      */
      CM_TOPHITS      *th      = NULL;        	/* top-scoring sequence hits                */

      nquery++;

      esl_stopwatch_Start(w);	                          

      /* seqfile may need to be rewound (multiquery mode) */
      if (nquery > 1) list->current = 0;

      /* Open the target profile database */
      status = cm_file_Open(cfg->cmfile, CMDBENV, FALSE, &cmfp, errbuf);
      if (status != eslOK) mpi_failure("Unexpected error %d in opening cm file %s.\n%s", status, cfg->cmfile, errbuf);  
  
      fprintf(ofp, "Query:       %s  [L=%ld]\n", qsq->name, (long) qsq->n);
      if (qsq->acc[0]  != 0) fprintf(ofp, "Accession:   %s\n", qsq->acc);
      if (qsq->desc[0] != 0) fprintf(ofp, "Description: %s\n", qsq->desc);

      /* Create processing pipeline and hit list */
      th  = cm_tophits_Create(); 
      pli = cm_pipeline_Create(go, abc, 100, 100, cfg->Z * qsq->n, cfg->Z_setby, CM_SCAN_MODELS); /* M_hint = 100, L_hint = 100 are just dummies for now */
      pli->cmfp = cmfp;  /* for four-stage input, pipeline needs <cmfp> */

      cm_pli_NewSeq(pli, qsq, nquery-1);

      /* Main loop: */
      while ((hstatus = next_block(cmfp, list, &block)) == eslOK)
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

	  if (mpistatus.MPI_TAG == INFERNAL_ERROR_TAG)
	    mpi_failure("MPI client %d raised error:\n%s\n", dest, mpi_buf);
	  if (mpistatus.MPI_TAG != INFERNAL_READY_TAG)
	    mpi_failure("Unexpected tag %d from %d\n", mpistatus.MPI_TAG, dest);
      
	  MPI_Send(&block, 3, MPI_LONG_LONG_INT, dest, INFERNAL_BLOCK_TAG, MPI_COMM_WORLD);
	}
      switch(hstatus)
	{
	case eslEFORMAT:   mpi_failure("bad file format in CM file %s",           cfg->cmfile); break;
	case eslEINCOMPAT: mpi_failure("CM file %s contains different alphabets", cfg->cmfile); break;
	case eslEOF: 	   /* do nothing */	                                                break;
	default:	   mpi_failure("Unexpected error %d in reading CMs from %s", hstatus, cfg->cmfile); break;
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

	  if (mpistatus.MPI_TAG == INFERNAL_ERROR_TAG)
	    mpi_failure("MPI client %d raised error:\n%s\n", dest, mpi_buf);
	  if (mpistatus.MPI_TAG != INFERNAL_READY_TAG)
	    mpi_failure("Unexpected tag %d from %d\n", mpistatus.MPI_TAG, dest);
	}

      /* merge the results of the search results */
      for (dest = 1; dest < cfg->nproc; ++dest)
	{
	  CM_PIPELINE     *mpi_pli   = NULL;
	  CM_TOPHITS      *mpi_th    = NULL;

	  /* send an empty block to signal the worker they are done */
	  MPI_Send(&block, 3, MPI_LONG_LONG_INT, dest, INFERNAL_BLOCK_TAG, MPI_COMM_WORLD);

	  /* wait for the results */
	  if ((status = cm_tophits_MPIRecv(dest, INFERNAL_TOPHITS_TAG, MPI_COMM_WORLD, &mpi_buf, &mpi_size, &mpi_th)) != eslOK)
	    mpi_failure("Unexpected error %d receiving tophits from %d", status, dest);

	  if ((status = cm_pipeline_MPIRecv(dest, INFERNAL_PIPELINE_TAG, MPI_COMM_WORLD, &mpi_buf, &mpi_size, go, &mpi_pli)) != eslOK)
	    mpi_failure("Unexpected error %d receiving pipeline from %d", status, dest);

	  cm_tophits_Merge(th, mpi_th);
	  cm_pipeline_Merge(pli, mpi_pli);

	  cm_pipeline_Destroy(mpi_pli, NULL);
	  cm_tophits_Destroy(mpi_th);
	}

      /* Print the results.  */
      cm_tophits_SortByScore(th);
      cm_tophits_Threshold(th, pli);

      /* tally up total number of hits and target coverage */
      pli->n_output = pli->pos_output = 0;
      for (i = 0; i < th->N; i++) {
	if ((th->hit[i]->flags & CM_HIT_IS_REPORTED) || (th->hit[i]->flags & CM_HIT_IS_INCLUDED)) { 
	  pli->n_output++;
	  pli->pos_output += abs(th->hit[i]->stop - th->hit[i]->start) + 1;
	}
      }

      cm_tophits_Targets(ofp, th, pli, textw); fprintf(ofp, "\n\n");
      if(pli->do_alignments) {
	if((status = cm_tophits_HitAlignments(ofp, th, pli, textw)) != eslOK) esl_fatal("Out of memory");
	fprintf(ofp, "\n\n");
      }

      if (tblfp)    cm_tophits_TabularTargets(tblfp,    qsq->name, qsq->acc, th, pli, (nquery == 1));

      esl_stopwatch_Stop(w);
      cm_pli_Statistics(ofp, pli, w);
      fprintf(ofp, "//\n");

      cm_file_Close(cmfp);
      cm_pipeline_Destroy(pli, NULL);
      cm_tophits_Destroy(th);
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

      if (mpistatus.MPI_TAG == INFERNAL_ERROR_TAG)
	mpi_failure("MPI client %d raised error:\n%s\n", dest, mpi_buf);
      if (mpistatus.MPI_TAG != INFERNAL_TERMINATING_TAG)
	mpi_failure("Unexpected tag %d from %d\n", mpistatus.MPI_TAG, dest);
    }

 /* Terminate outputs - any last words?
   */
  if (tblfp)    cm_tophits_TabularTail(tblfp,    "cmscan", CM_SCAN_MODELS, cfg->seqfile, cfg->cmfile, go);
  if (ofp)      fprintf(ofp, "[ok]\n");

  esl_stopwatch_Stop(mw);
  esl_stopwatch_Display(stdout, mw, "Total runtime:");

  /* Cleanup - prepare for successful exit
   */
  free(list);
  if (mpi_buf != NULL) free(mpi_buf);

  p7_bg_Destroy(bg);

  esl_sq_Destroy(qsq);
  esl_stopwatch_Destroy(w);
  esl_stopwatch_Destroy(mw);
  esl_alphabet_Destroy(abc);
  esl_sqfile_Close(sqfp);

  if (ofp != stdout) fclose(ofp);
  if (tblfp)         fclose(tblfp);

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
  CM_FILE         *cmfp     = NULL;		 /* open CM database file                           */
  CM_t            *cm       = NULL;              /* the CM                                          */
  CMConsensus_t   *cmcons   = NULL;              /* CM consensus information                        */
  ESL_ALPHABET    *abc      = NULL;              /* sequence alphabet                               */
  P7_OPROFILE     *om       = NULL;		 /* target profile                                  */
  P7_PROFILE      *gm       = NULL;              /* generic query profile HMM                       */
  ESL_STOPWATCH   *w        = NULL;              /* timing                                          */
  ESL_SQ          *qsq      = NULL;		 /* query sequence                                  */
  int              status   = eslOK;
  int              hstatus  = eslOK;
  int              sstatus  = eslOK;
  int              cm_clen, cm_W;          /* consensus, window length for current CM */        
  float            gfmu, gflambda;         /* glocal fwd mu, lambda for current hmm filter */
  off_t            cm_offset;              /* file offset for current CM */
  int              cm_config_opts;         /* sent to the pipeline so we can appropriately configure CMs */
  float           *p7_evparam;             /* E-value parameters for the p7 filter */
  int              prv_ntophits;

  char            *mpi_buf  = NULL;              /* buffer used to pack/unpack structures */
  int              mpi_size = 0;                 /* size of the allocated buffer */
  int              nquery   = 0;

  MPI_Status       mpistatus;
  char             errbuf[eslERRBUFSIZE];

  w = esl_stopwatch_Create();

  /* Open the target profile database to get the sequence alphabet */
  /* Open the target profile database to get the sequence alphabet */
  status = cm_file_Open(cfg->cmfile, CMDBENV, FALSE, &cmfp, errbuf);
  if      (status == eslENOTFOUND) mpi_failure("File existence/permissions problem in trying to open CM file %s.\n%s\n", cfg->cmfile, errbuf);
  else if (status == eslEFORMAT)   mpi_failure("File format problem in trying to open CM file %s.\n%s\n",                cfg->cmfile, errbuf);
  else if (status != eslOK)        mpi_failure("Unexpected error %d in opening CM file %s.\n%s\n",               status, cfg->cmfile, errbuf);  
  if      (cmfp->do_gzip)          mpi_failure("Reading gzipped CM files is not supported");
  if      (cmfp->do_stdin)         mpi_failure("Reading CM files from stdin is not supported");
  if (! cmfp->is_pressed)          mpi_failure("Failed to open binary auxfiles for %s: use cmpress first\n",             cmfp->fname);

  hstatus = cm_file_Read(cmfp, FALSE, &abc, NULL);
  if(hstatus == eslEFORMAT)  mpi_failure("bad file format in CM file %s\n%s",           cfg->cmfile, cmfp->errbuf);
  else if (hstatus != eslOK) mpi_failure("Unexpected error in reading CMs from %s\n%s", cfg->cmfile, cmfp->errbuf); 

  /* Determine database size: default is to updated as we read target CMs */
  if(esl_opt_IsUsed(go, "-Z")) { 
    cfg->Z       = (int64_t) esl_opt_GetReal(go, "-Z");
    cfg->Z_setby = CM_ZSETBY_OPTION; 
  }
  else { 
    if(cmfp->ssi == NULL) mpi_failure("Failed to open SSI index for CM file: %s\n", cmfp->fname);
    cfg->Z = (int64_t) cmfp->ssi->nprimary;
    cfg->Z_setby = CM_ZSETBY_SSI_AND_QLENGTH; /* we will multiply Z by each query sequence length */
  }
  cm_file_Close(cmfp);

  /* Open the query sequence database */
  status = esl_sqfile_OpenDigital(abc, cfg->seqfile, seqfmt, NULL, &sqfp);
  if      (status == eslENOTFOUND) mpi_failure("Failed to open sequence file %s for reading\n",      cfg->seqfile);
  else if (status == eslEFORMAT)   mpi_failure("Sequence file %s is empty or misformatted\n",        cfg->seqfile);
  else if (status == eslEINVAL)    mpi_failure("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        mpi_failure("Unexpected error %d opening sequence file %s\n", status, cfg->seqfile);

  qsq            = esl_sq_CreateDigital(abc);
  bg             = p7_bg_Create(abc);
  cm_config_opts = determine_config_opts(go);

  ESL_ALLOC(p7_evparam, sizeof(float) * CM_p7_NEVPARAM);

  /* Outside loop: over each query sequence in <seqfile>. */
  while ((sstatus = esl_sqio_Read(sqfp, qsq)) == eslOK)
    {
      CM_PIPELINE     *pli     = NULL;		/* processing pipeline                      */
      CM_TOPHITS      *th      = NULL;        	/* top-scoring sequence hits                */

      MSV_BLOCK        block;

      nquery++;

      esl_stopwatch_Start(w);

      status = 0;
      MPI_Send(&status, 1, MPI_INT, 0, INFERNAL_READY_TAG, MPI_COMM_WORLD);

      /* Open the target profile database */
      status = cm_file_Open(cfg->cmfile, CMDBENV, FALSE, &cmfp, errbuf);
      if (status != eslOK) mpi_failure("Unexpected error %d in opening cm file %s.\n%s", status, cfg->cmfile, errbuf);  
  
      /* Create processing pipeline and hit list */
      th  = cm_tophits_Create(); 
      pli = cm_pipeline_Create(go, abc, 100, 100, cfg->Z * qsq->n, cfg->Z_setby, CM_SCAN_MODELS); /* M_hint = 100, L_hint = 100 are just dummies for now */
      pli->cmfp = cmfp;  /* for four-stage input, pipeline needs <cmfp> */

      cm_pli_NewSeq(pli, qsq, nquery-1);

      /* receive a sequence block from the master */
      MPI_Recv(&block, 3, MPI_LONG_LONG_INT, 0, INFERNAL_BLOCK_TAG, MPI_COMM_WORLD, &mpistatus);
      while (block.count > 0)
	{
	  uint64_t length = 0;
	  uint64_t count  = block.count;

	  hstatus = cm_p7_oprofile_Position(cmfp, block.offset);
	  if (hstatus != eslOK) mpi_failure("Cannot position optimized model to %ld\n", block.offset);

	  while (count > 0 && 
		 (hstatus = cm_p7_oprofile_ReadMSV(cmfp, TRUE, &abc, &cm_offset, &cm_clen, &cm_W, &gfmu, &gflambda, &om)) == eslOK)
	    {
	      length = om->eoff - block.offset + 1;

	      esl_vec_FCopy(om->evparam, p7_NEVPARAM, p7_evparam);
	      p7_evparam[CM_p7_GFMU]     = gfmu;
	      p7_evparam[CM_p7_GFLAMBDA] = gflambda;
	      gm     = NULL; /* this will get filled in cm_Pipeline() only if necessary */
	      cm     = NULL; /* this will get filled in cm_Pipeline() only if necessary */
	      cmcons = NULL; /* this will get filled in cm_Pipeline() only if necessary */
	      if((status = cm_pli_NewModel(pli, CM_NEWMODEL_MSV, 
					   cm,                                   /* this is NULL b/c we don't have one yet */
					   cm_clen, cm_W,                        /* we read these in cm_p7_oprofile_ReadMSV() */
					   FALSE, FALSE, NULL, NULL, NULL, NULL, /* all these are irrelevant in CM_NEWMODEL_MSV mode */
					   om, bg)) != eslOK) mpi_failure(pli->errbuf);

	      prv_ntophits = th->N;
	      if((status = cm_Pipeline(pli, cm_offset, cm_config_opts, om, bg, p7_evparam, qsq, FALSE, th, &gm, &cm, &cmcons)) != eslOK)
		mpi_failure("cm_pipeline() failed unexpected with status code %d\n%s", status, pli->errbuf);

	      if(th->N != prv_ntophits) { 
		cm_tophits_ComputeEvalues(th, (double) cm->expA[pli->final_cm_exp_mode]->cur_eff_dbsize, prv_ntophits);
	      }
	      cm_pipeline_Reuse(pli);

	      if(cmcons != NULL) FreeCMConsensus(cmcons);
	      if(pli->fsmx != NULL && cm != NULL) { cm_FreeScanMatrix(cm, pli->fsmx); pli->fsmx = NULL; }
	      if(pli->smx  != NULL && cm != NULL) { cm_FreeScanMatrix(cm, pli->smx);  pli->smx  = NULL; }
	      if(cm  != NULL) FreeCM(cm);
	      if(om  != NULL) p7_oprofile_Destroy(om);
	      if(gm  != NULL) p7_profile_Destroy(gm);

	      --count;
	    }
	  /* check the status of reading the msv filter */
	  if (count > 0)              
	    {
	      switch(hstatus)
		{
		case eslEFORMAT:    mpi_failure("bad file format in HMM file %s\n%s",          cfg->cmfile, cmfp->errbuf); break;
		case eslEINCOMPAT:  mpi_failure("HMM file %s contains different alphabets",    cfg->cmfile); break;
		case eslOK:         
		case eslEOF:        mpi_failure("Block count mismatch - expected %ld found %ld at offset %ld\n", block.count, block.count-count, block.offset); break;
		default:  	    mpi_failure("Unexpected error %d in reading HMMs from %s\n%s", hstatus, cfg->cmfile, cmfp->errbuf); break;
		}
	    }

	  /* lets do a little bit of sanity checking here to make sure the blocks are the same */
	  if (block.length != length) mpi_failure("Block length mismatch - expected %ld found %ld at offset %ld\n", block.length, length, block.offset);

	  /* inform the master we need another block of sequences */
	  status = 0;
	  MPI_Send(&status, 1, MPI_INT, 0, INFERNAL_READY_TAG, MPI_COMM_WORLD);

	  /* wait for the next block of sequences */
	  MPI_Recv(&block, 3, MPI_LONG_LONG_INT, 0, INFERNAL_BLOCK_TAG, MPI_COMM_WORLD, &mpistatus);
	}

      esl_stopwatch_Stop(w);

      /* Send the top hits back to the master. */
      cm_tophits_MPISend(th, 0, INFERNAL_TOPHITS_TAG, MPI_COMM_WORLD,  &mpi_buf, &mpi_size);
      cm_pipeline_MPISend(pli, 0, INFERNAL_PIPELINE_TAG, MPI_COMM_WORLD,  &mpi_buf, &mpi_size);

      cm_file_Close(cmfp);
      cm_pipeline_Destroy(pli, NULL);
      cm_tophits_Destroy(th);
      esl_sq_Reuse(qsq);
    } /* end outer loop over query HMMs */
  if (sstatus == eslEFORMAT) 
    mpi_failure("Parse failed (sequence file %s):\n%s\n", sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
  else if (sstatus != eslEOF)
    mpi_failure("Unexpected error %d reading sequence file %s", sstatus, sqfp->filename);

  status = 0;
  MPI_Send(&status, 1, MPI_INT, 0, INFERNAL_TERMINATING_TAG, MPI_COMM_WORLD);

  if (mpi_buf != NULL) free(mpi_buf);

  p7_bg_Destroy(bg);

  esl_sq_Destroy(qsq);
  esl_stopwatch_Destroy(w);
  esl_alphabet_Destroy(abc);
  esl_sqfile_Close(sqfp);

  return eslOK;

 ERROR: 
  mpi_failure("out of memory");
  return status; /* NEVER REACHED */
}
#endif /*HAVE_MPI*/

static int
serial_loop(WORKER_INFO *info, CM_FILE *cmfp)
{
  int            status;

  CM_t             *cm         = NULL; 
  ESL_ALPHABET     *abc        = NULL;
  P7_OPROFILE      *om         = NULL;      /* optimized query profile HMM            */
  P7_PROFILE       *gm         = NULL;      /* generic   query profile HMM            */
  CMConsensus_t    *cmcons     = NULL;
  int               init_ntophits;          /* number of top hits before first cm_Pipeline() call */
  int               init_npastfwd;          /* number of hits past local Fwd filter before first cm_Pipeline() call */
  int               prv_ntophits;           /* number of top hits before each cm_Pipeline() call */
  int               cm_clen, cm_W;          /* consensus, window length for current CM */        
  float             gfmu, gflambda;         /* glocal fwd mu, lambda for current hmm filter */
  off_t             cm_offset;              /* file offset for current CM */
  ESL_SQ           *termsq  = esl_sq_CreateDigital(info->qsq->abc); /* terminal sequence, first or final ESL_MIN(pli->maxW, dbsq->n) residues of dbsq */

  /* Main loop: */
  while ((status = cm_p7_oprofile_ReadMSV(cmfp, TRUE, &abc, &cm_offset, &cm_clen, &cm_W, &gfmu, &gflambda, &om)) == eslOK)
    {
      esl_vec_FCopy(om->evparam, p7_NEVPARAM, info->p7_evparam);
      info->p7_evparam[CM_p7_GFMU]     = gfmu;
      info->p7_evparam[CM_p7_GFLAMBDA] = gflambda;
      gm     = NULL; /* this will get filled in cm_Pipeline() only if necessary */
      cm     = NULL; /* this will get filled in cm_Pipeline() only if necessary */
      cmcons = NULL; /* this will get filled in cm_Pipeline() only if necessary */
      if((status = cm_pli_NewModel(info->pli, CM_NEWMODEL_MSV, 
				   cm,                                   /* this is NULL b/c we don't have one yet */
				   cm_clen, cm_W,                        /* we read these in cm_p7_oprofile_ReadMSV() */
				   FALSE, FALSE, NULL, NULL, NULL, NULL, /* all these are irrelevant in CM_NEWMODEL_MSV mode */
				   om, info->bg)) != eslOK) cm_Fail(info->pli->errbuf);

      init_ntophits = info->th->N;
      init_npastfwd = info->pli->n_past_fwd;
      if((status = cm_Pipeline(info->pli, cm_offset, info->cm_config_opts, om, info->bg, info->p7_evparam, info->qsq, /*i_am_terminal=*/FALSE, info->th, &gm, &cm, &cmcons)) != eslOK)
	cm_Fail("cm_pipeline() failed unexpected with status code %d\n%s", status, info->pli->errbuf);
      cm_pipeline_Reuse(info->pli); /* prepare for next search */
      if(info->in_rc && info->th->N != init_ntophits) cm_tophits_UpdateHitPositions(info->th, init_ntophits, info->qsq->start, info->in_rc);

      if(info->pli->research_ends && (init_npastfwd != info->pli->n_past_fwd)) { 
	/* research first info->pli->maxW residues with local envelope defn */
	copy_subseq(info->qsq, termsq, 1, ESL_MIN(info->pli->maxW, info->qsq->n), info->in_rc);
	prv_ntophits = info->th->N;
	if((status = cm_Pipeline(info->pli, cm_offset, info->cm_config_opts, om, info->bg, info->p7_evparam, termsq, /*i_am_terminal=*/TRUE, info->th, &gm, &cm, &cmcons)) != eslOK) 
	  cm_Fail("cm_pipeline() failed unexpected with status code %d\n%s\n", status, info->pli->errbuf);
	cm_pipeline_Reuse(info->pli); /* prepare for next search */
	if(info->in_rc && info->th->N != prv_ntophits) cm_tophits_UpdateHitPositions(info->th, prv_ntophits, termsq->start, info->in_rc);

	/* research final info->pli->maxW residues with local envelope defn */
	if(info->qsq->n > info->pli->maxW) { 
	  copy_subseq(info->qsq, termsq, info->qsq->n - info->pli->maxW + 1, info->pli->maxW, info->in_rc);
	  prv_ntophits = info->th->N;
	  if((status = cm_Pipeline(info->pli, cm_offset, info->cm_config_opts, om, info->bg, info->p7_evparam, termsq, /*i_am_terminal=*/TRUE, info->th, &gm, &cm, &cmcons)) != eslOK) 
	    cm_Fail("cm_pipeline() failed unexpected with status code %d\n%s\n", status, info->pli->errbuf);
	  cm_pipeline_Reuse(info->pli); /* prepare for next search */
	  if(info->th->N != prv_ntophits) cm_tophits_UpdateHitPositions(info->th, prv_ntophits, termsq->start, info->in_rc);
	}
      }	

      if(info->th->N != init_ntophits) { 
	cm_tophits_ComputeEvalues(info->th, (double) cm->expA[info->pli->final_cm_exp_mode]->cur_eff_dbsize, init_ntophits);
      }

      if(cmcons != NULL) FreeCMConsensus(cmcons);
      if(info->pli->fsmx != NULL && cm != NULL) { cm_FreeScanMatrix(cm, info->pli->fsmx); info->pli->fsmx = NULL; }
      if(info->pli->smx  != NULL && cm != NULL) { cm_FreeScanMatrix(cm, info->pli->smx);  info->pli->smx  = NULL; }
      if(cm  != NULL) FreeCM(cm);
      if(om  != NULL) p7_oprofile_Destroy(om);
      if(gm  != NULL) p7_profile_Destroy(gm);
    } /* end of while(cm_p7_oprofile_ReadMSV() == eslOK) */

  esl_alphabet_Destroy(abc);
  esl_sq_Destroy(termsq);

  return status;
}

#ifdef HMMER_THREADS
static int
thread_loop(ESL_THREADS *obj, ESL_WORK_QUEUE *queue, CM_FILE *cmfp)
{
  int  status   = eslOK;
  int  sstatus  = eslOK;
  int  eofCount = 0;
  CM_P7_OM_BLOCK  *block;
  ESL_ALPHABET    *abc = NULL;
  void            *newBlock;

  esl_workqueue_Reset(queue);
  esl_threads_WaitForStart(obj);

  status = esl_workqueue_ReaderUpdate(queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue reader failed");
      
  /* Main loop: */
  while (sstatus == eslOK)
    {
      block = (CM_P7_OM_BLOCK *) newBlock;
      sstatus = cm_p7_oprofile_ReadBlockMSV(cmfp, &abc, block);
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
  WORKER_INFO     *info;
  ESL_THREADS     *obj;
  CM_P7_OM_BLOCK  *block;
  void            *newBlock;

  CM_t             *cm         = NULL; 
  ESL_ALPHABET     *abc        = NULL;
  P7_PROFILE       *gm         = NULL;      /* generic   query profile HMM            */
  CMConsensus_t    *cmcons     = NULL;
  int               init_ntophits;          /* number of top hits before first cm_Pipeline() call */
  int               init_npastfwd;          /* number of hits past local Fwd filter before first cm_Pipeline() call */
  int               prv_ntophits;           /* number of top hits before each cm_Pipeline() call */
  int               cm_clen, cm_W;          /* consensus, window length for current CM */        
  float             gfmu, gflambda;         /* glocal fwd mu, lambda for current hmm filter */
  off_t             cm_offset;              /* file offset for current CM */
  ESL_SQ           *termsq = NULL;          /* terminal sequence, first or final ESL_MIN(pli->maxW, dbsq->n) residues of dbsq */
  
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
  termsq  = esl_sq_CreateDigital(info->qsq->abc);

  /* loop until all blocks have been processed */
  block = (CM_P7_OM_BLOCK *) newBlock;
  while (block->count > 0)
    {
      /* Main loop: */
      for (i = 0; i < block->count; ++i)
	{
	  P7_OPROFILE *om = block->list[i];
	  cm_offset       = block->cm_offsetA[i];
	  cm_clen         = block->cm_clenA[i];
	  cm_W            = block->cm_WA[i];
	  gfmu            = block->gfmuA[i];
	  gflambda        = block->gflambdaA[i];

	  esl_vec_FCopy(om->evparam, p7_NEVPARAM, info->p7_evparam);
	  info->p7_evparam[CM_p7_GFMU]     = gfmu;
	  info->p7_evparam[CM_p7_GFLAMBDA] = gflambda;
	  gm     = NULL; /* this will get filled in cm_Pipeline() only if necessary */
	  cm     = NULL; /* this will get filled in cm_Pipeline() only if necessary */
	  cmcons = NULL; /* this will get filled in cm_Pipeline() only if necessary */
	  if((status = cm_pli_NewModel(info->pli, CM_NEWMODEL_MSV, 
				       cm,                                   /* this is NULL b/c we don't have one yet */
				       cm_clen, cm_W,                        /* we read these in cm_p7_oprofile_ReadMSV() */
				       FALSE, FALSE, NULL, NULL, NULL, NULL, /* all these are irrelevant in CM_NEWMODEL_MSV mode */
				       om, info->bg)) != eslOK) cm_Fail(info->pli->errbuf);

	  init_ntophits = info->th->N;
	  init_npastfwd = info->pli->n_past_fwd;
	  if((status = cm_Pipeline(info->pli, cm_offset, info->cm_config_opts, om, info->bg, info->p7_evparam, info->qsq, /*i_am_terminal=*/FALSE, info->th, &gm, &cm, &cmcons)) != eslOK)
	    cm_Fail("cm_pipeline() failed unexpected with status code %d\n%s", status, info->pli->errbuf);
	  cm_pipeline_Reuse(info->pli);
	  if(info->in_rc && info->th->N != init_ntophits) cm_tophits_UpdateHitPositions(info->th, init_ntophits, info->qsq->start, info->in_rc);

	  if(info->pli->research_ends && (init_npastfwd != info->pli->n_past_fwd)) { 
	    /* research first info->pli->maxW residues with local envelope defn */
	    copy_subseq(info->qsq, termsq, 1, ESL_MIN(info->pli->maxW, info->qsq->n), info->in_rc);
	    prv_ntophits = info->th->N;
	    if((status = cm_Pipeline(info->pli, cm_offset, info->cm_config_opts, om, info->bg, info->p7_evparam, termsq, /*i_am_terminal=*/TRUE, info->th, &gm, &cm, &cmcons)) != eslOK) 
	      cm_Fail("cm_pipeline() failed unexpected with status code %d\n%s\n", status, info->pli->errbuf);
	    cm_pipeline_Reuse(info->pli); /* prepare for next search */
	    if(info->in_rc && info->th->N != prv_ntophits) cm_tophits_UpdateHitPositions(info->th, prv_ntophits, termsq->start, info->in_rc);
	    
	    /* research final info->pli->maxW residues with local envelope defn */
	    if(info->qsq->n > info->pli->maxW) { 
	      copy_subseq(info->qsq, termsq, info->qsq->n - info->pli->maxW + 1, info->pli->maxW, info->in_rc);
	      prv_ntophits = info->th->N;
	      if((status = cm_Pipeline(info->pli, cm_offset, info->cm_config_opts, om, info->bg, info->p7_evparam, termsq, /*i_am_terminal=*/TRUE, info->th, &gm, &cm, &cmcons)) != eslOK) 
		cm_Fail("cm_pipeline() failed unexpected with status code %d\n%s\n", status, info->pli->errbuf);
	      cm_pipeline_Reuse(info->pli); /* prepare for next search */
	      if(info->th->N != prv_ntophits) cm_tophits_UpdateHitPositions(info->th, prv_ntophits, termsq->start, info->in_rc);
	    }
	  }

	  if(info->th->N != init_ntophits) { 
	    cm_tophits_ComputeEvalues(info->th, (double) cm->expA[info->pli->final_cm_exp_mode]->cur_eff_dbsize, init_ntophits);
	  }
		
	  if(cmcons != NULL) FreeCMConsensus(cmcons);
	  if(info->pli->fsmx != NULL && cm != NULL) { cm_FreeScanMatrix(cm, info->pli->fsmx); info->pli->fsmx = NULL; }
	  if(info->pli->smx  != NULL && cm != NULL) { cm_FreeScanMatrix(cm, info->pli->smx);  info->pli->smx  = NULL; }
	  if(cm != NULL) FreeCM(cm);
	  if(om != NULL) p7_oprofile_Destroy(om);
	  if(gm != NULL) p7_profile_Destroy(gm);
	  
	  block->list[i] = NULL;
	  block->cm_offsetA[i] = 0;
	  block->cm_clenA[i]   = 0;
	  block->cm_WA[i]      = 0;
	  block->gfmuA[i]      = 0.;
	  block->gflambdaA[i]  = 0.;
	}
      
      status = esl_workqueue_WorkerUpdate(info->queue, block, &newBlock);
      if (status != eslOK) esl_fatal("Work queue worker failed");
      
      block = (CM_P7_OM_BLOCK *) newBlock;
    }
  
  status = esl_workqueue_WorkerUpdate(info->queue, block, NULL);
  if (status != eslOK) esl_fatal("Work queue worker failed");

  esl_threads_Finished(obj, workeridx);

  esl_alphabet_Destroy(abc);
  esl_sq_Destroy(termsq);

  return;
}
#endif   /* HMMER_THREADS */

 
/* Function:  copy_subseq()
 * Incept:    EPN, Tue Aug  9 11:13:02 2011
 *
 * Purpose: Copy a subsequence of an existing sequence <src_sq>
 *           starting at position <i>, of length <L> to another
 *           sequence object <dest_sq>. Copy only residues from
 *           <i>..<i>+<L>-1. <dest_sq> must be pre-allocated.
 *
 * Returns: eslOK on success. 
 */
void
copy_subseq(const ESL_SQ *src_sq, ESL_SQ *dest_sq, int64_t i, int64_t L, int in_rc)
{ 
  ESL_DASSERT1((((src_sq->start <= src_sq->end) && (! in_rc)) || 
		((src_sq->start >= src_sq->end) && (  in_rc))));
  assert       (((src_sq->start <= src_sq->end) && (! in_rc)) || 
		((src_sq->start >= src_sq->end) && (  in_rc)));

  esl_sq_Reuse(dest_sq);
  esl_sq_GrowTo(dest_sq, L);
  memcpy((void*)(dest_sq->dsq+1), src_sq->dsq+i, L * sizeof(ESL_DSQ));
  dest_sq->dsq[0] = dest_sq->dsq[L+1] = eslDSQ_SENTINEL;
  dest_sq->n     = L;

  if(in_rc) { 
    dest_sq->start = src_sq->start  - i + 1;
    dest_sq->end   = dest_sq->start - L + 1;
  }
  else { 
    dest_sq->start = src_sq->start  + i - 1;
    dest_sq->end   = dest_sq->start + L - 1;
  }

  esl_sq_SetName     (dest_sq, src_sq->name);
  esl_sq_SetAccession(dest_sq, src_sq->acc);
  esl_sq_SetDesc     (dest_sq, src_sq->desc);

  return;
}

/* Function:  determine_config_opts()
 * Incept:    EPN, Wed Jul  6 13:45:44 2011
 *
 * Purpose:  Determine what a CM's <config_opts> flags 
 *           should be based on command-line options.
 *
 * Returns: void.
 */
int
determine_config_opts(const ESL_GETOPTS *go)
{ 
  int config_opts = 0;

  if(! esl_opt_GetBoolean(go, "-g")) { 
    if(! esl_opt_GetBoolean(go, "--cp9gloc")) { 
      config_opts |= CM_CONFIG_LOCAL;
      config_opts |= CM_CONFIG_HMMLOCAL;
      if(! esl_opt_GetBoolean(go, "--cp9noel")) config_opts |= CM_CONFIG_HMMEL; 
    }
  }

  return config_opts;
}

/*****************************************************************
 * @LICENSE@
 *****************************************************************/


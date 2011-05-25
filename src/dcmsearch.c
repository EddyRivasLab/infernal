/* dcmsearch: search CM(s) against a nucleotide sequence database,
 *            using profile HMM(s) to prefilter the database.
 * 
 * Based on HMMER 3's nhmmer.c, which was based on hmmsearch.c.
 * EPN, Fri Sep 24 10:58:08 2010
 */
#include "esl_config.h"
#include "p7_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_getopts.h"
#include "esl_scorematrix.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_ssi.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#ifdef HMMER_THREADS
#include <unistd.h>
#include "esl_threads.h"
#include "esl_workqueue.h"
#endif /*HMMER_THREADS*/

#include "hmmer.h"
#include "funcs.h"		/* external functions                   */
#include "structs.h"		/* data structures, macros, #define's   */

/* set the max residue count to 1 meg when reading a block */
#define CMSEARCH_MAX_RESIDUE_COUNT 100000 /* differs from HMMER's default which is MAX_RESIDUE_COUNT from esl_sqio_(ascii|ncbi).c */

typedef struct {
#ifdef HMMER_THREADS
  ESL_WORK_QUEUE   *queue;
#endif /*HMMER_THREADS*/
  CM_PIPELINE      *pli;         /* work pipeline                           */
  CM_TOPHITS       *th;          /* top hit results                         */
  CM_t             *cm;          /* a covariance model                      */
  P7_BG           **bgA;         /* null models                              */
  P7_OPROFILE     **omA;         /* optimized query profile HMMs            */
  P7_PROFILE      **gmA;         /* generic   query profile HMMs            */
  float           **p7_evparamAA;/* [0..nhmm-1][0..CM_p7_NEVPARAM] E-value parameters */
  int               nhmm;        /* number of HMM filters, size of omA, gmA, bgA */
} WORKER_INFO;

#define REPOPTS     "-E,-T,--cut_ga,--cut_nc,--cut_tc"
#define DOMREPOPTS  "--domE,--domT,--cut_ga,--cut_nc,--cut_tc"
#define INCOPTS     "--incE,--incT,--cut_ga,--cut_nc,--cut_tc"
#define INCDOMOPTS  "--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
#define THRESHOPTS  "-E,-T,--domE,--domT,--incE,--incT,--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
#define XFASTMIDMAXOPTS "--max,--F1,--F2,--F3,--F4,--F5,--F6,--noF1,--noF2,--noF3,--nogfwd,--nocyk,--nohmm,--hmm"
#define TIMINGOPTS  "--time-F1,--time-F2,--time-F3,--time-dF3,--time-bfil,--time-F4"

#define CPUOPTS     NULL
#define MPIOPTS     NULL

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles   reqs   incomp              help                                                      docgroup*/
  { "-h",           eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "show brief help on version and usage",                         1 },
  /* Control of output */
  { "-o",           eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "direct output to file <f>, not stdout",                        2 },
  { "-A",           eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save multiple alignment of all hits to file <s>",              2 },
  { "--tblout",     eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save parseable table of per-sequence hits to file <s>",        2 },
//  { "--domtblout",  eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save parseable table of per-domain hits to file <s>",          2 },
  { "--acc",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "prefer accessions over names in output",                       2 },
  { "--noali",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "don't output alignments, so output is smaller",                2 },
  { "--notextw",    eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL, "--textw",        "unlimit ASCII text output line width",                         2 },
  { "--textw",      eslARG_INT,    "120", NULL, "n>=120",NULL,  NULL, "--notextw",      "set max width of ASCII text output lines",                     2 },
  /* Control of scoring system */
  { "--popen",      eslARG_REAL,  "0.02", NULL, "0<=x<0.5",NULL,  NULL,  NULL,              "gap open probability",                                         3 },
  { "--pextend",    eslARG_REAL,   "0.4", NULL, "0<=x<1",  NULL,  NULL,  NULL,              "gap extend probability",                                       3 },
  { "--mxfile",     eslARG_INFILE,  NULL, NULL, NULL,      NULL,  NULL,  NULL,              "substitution score matrix [default: BLOSUM62]",                3 },
  /* Control of reporting thresholds */
  { "-E",           eslARG_REAL,  "10.0", NULL, "x>0",   NULL,  NULL,  REPOPTS,         "report sequences <= this E-value threshold in output",         4 },
  { "-T",           eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  REPOPTS,         "report sequences >= this score threshold in output",           4 },
  /* Control of inclusion (significance) thresholds */
  { "--incE",       eslARG_REAL,  "0.01", NULL, "x>0",   NULL,  NULL,  INCOPTS,         "consider sequences <= this E-value threshold as significant",  5 },
  { "--incT",       eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  INCOPTS,         "consider sequences >= this score threshold as significant",    5 },
  /* Model-specific thresholding for both reporting and inclusion */
  { "--cut_ga",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's GA gathering cutoffs to set all thresholding",   6 },
  { "--cut_nc",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's NC noise cutoffs to set all thresholding",       6 },
  { "--cut_tc",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use profile's TC trusted cutoffs to set all thresholding",     6 },
  /* Control of acceleration pipeline */
  { "--fZ",         eslARG_REAL,    NULL, NULL, NULL,    NULL,  NULL,  "--rfam",        "set heuristic filters to defaults used for a database size of <x> Mb", 7 },
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
  { "-Z",           eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL,  NULL,            "set database size in *Mb* to <x> for E-value calculations",   12 },
  { "--seed",       eslARG_INT,   "181",  NULL, "n>=0",  NULL,  NULL,  NULL,            "set RNG seed to <n> (if 0: one-time arbitrary seed)",         12 },
  { "--w_beta",     eslARG_REAL,   NULL,  NULL, NULL,    NULL,  NULL,  NULL,            "tail mass at which window length is determined",               12 },
  { "--w_length",   eslARG_INT,    NULL,  NULL, NULL,    NULL,  NULL,  NULL,            "window length ",                                              12 },
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

/* Not used, but retained because esl option-handling code errors if it isn't kept here.  Placed in group 99 so it doesn't print to help*/
  { "--domZ",       eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL,  NULL,            "Not used",   99 },
  { "--domE",       eslARG_REAL,  "10.0", NULL, "x>0",   NULL,  NULL,  DOMREPOPTS,      "Not used",   99 },
  { "--domT",       eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  DOMREPOPTS,      "Not used",   99 },
  { "--incdomE",    eslARG_REAL,  "0.01", NULL, "x>0",   NULL,  NULL,  INCDOMOPTS,      "Not used",    99 },
  { "--incdomT",    eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  INCDOMOPTS,      "Not used",      99 },
/* will eventually bring these back, but store in group 99 for now, so they don't print to help*/
  { "--qformat",    eslARG_STRING,  NULL, NULL, NULL,    NULL,  NULL,  NULL,            "assert query <seqfile> is in format <s>: no autodetection",   99 },
  { "--tformat",    eslARG_STRING,  NULL, NULL, NULL,    NULL,  NULL,  NULL,            "assert target <seqfile> is in format <s>>: no autodetection", 99 },

#ifdef HMMER_THREADS
  { "--cpu",        eslARG_INT, NULL,"HMMER_NCPU","n>=0",NULL,  NULL,  CPUOPTS,         "number of parallel CPU workers to use for multithreads",      12 },
#endif
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};


/* struct cfg_s : "Global" application configuration shared by all threads/processes
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 */
struct cfg_s {
  char            *dbfile;           /* target sequence database file                   */
  char            *cmfile;           /* query HMM file                                  */
  int64_t          Z;                /* database size, in nucleotides                   */

  int              do_mpi;            /* TRUE if we're doing MPI parallelization         */
  int              nproc;             /* how many MPI processes, total                   */
  int              my_rank;           /* who am I, in 0..nproc-1                         */
};

static char usage[]  = "[options] <query cmfile> <target seqfile>";
static char banner[] = "search a sequence database with an RNA CM";

static int  serial_master(ESL_GETOPTS *go, struct cfg_s *cfg);
static int  serial_loop  (WORKER_INFO *info, ESL_SQFILE *dbfp);

static int  determine_qdbs(CM_t *cm, double beta, int **ret_dmin, int **ret_dmax);
static int  update_hit_positions(CM_TOPHITS *th, int istart, int64_t offset);

#ifdef HMMER_THREADS
#define BLOCK_SIZE 1000

static int  thread_loop(WORKER_INFO *info, ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQFILE *dbfp);
static void pipeline_thread(void *arg);
#endif /*HMMER_THREADS*/


static void
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_cmfile, char **ret_seqfile)
{
  ESL_GETOPTS *go = NULL;

  if ((go = esl_getopts_Create(options))     == NULL)     p7_Die("Internal failure creating options object");
  if (esl_opt_ProcessEnvironment(go)         != eslOK)  { printf("Failed to process environment: %s\n", go->errbuf); goto ERROR; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { printf("Failed to parse command line: %s\n", go->errbuf); goto ERROR; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { printf("Failed to parse command line: %s\n", go->errbuf); goto ERROR; }
 
  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) {
    cm_banner(stdout, argv[0], banner);
    esl_usage(stdout, argv[0], usage);
    puts("\nBasic options:");
    esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 120=textwidth*/
    
    puts("\nOptions directing output:");
    esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
    
    //      puts("\nOptions controlling scoring system:");
    //     esl_opt_DisplayHelp(stdout, go, 3, 2, 80);
    
    puts("\nOptions controlling reporting thresholds:");
    esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 
    
    puts("\nOptions controlling inclusion (significance) thresholds:");
    esl_opt_DisplayHelp(stdout, go, 5, 2, 80); 
    
    puts("\nOptions controlling acceleration heuristics:");
    esl_opt_DisplayHelp(stdout, go, 7, 2, 80);

    puts("\nOptions controlling p7 HMM statistics calibrations:");
    esl_opt_DisplayHelp(stdout, go, 8, 2, 80);
    
    puts("\nOptions for reading a p7 HMM from a file:");
    esl_opt_DisplayHelp(stdout, go, 13, 2, 80); 

    puts("\nOptions from Infernal 1.0.2 cmsearch:");
    esl_opt_DisplayHelp(stdout, go, 20, 2, 80); 

    puts("\nOther expert options:");
    esl_opt_DisplayHelp(stdout, go, 12, 2, 80); 
    exit(0);

  }

  if (esl_opt_ArgNumber(go)                  != 2)     { puts("Incorrect number of command line arguments.");      goto ERROR; }
  if ((*ret_cmfile = esl_opt_GetArg(go, 1))  == NULL)  { puts("Failed to get <cmfile> argument on command line"); goto ERROR; }
  if ((*ret_seqfile = esl_opt_GetArg(go, 2)) == NULL)  { puts("Failed to get <seqfile> argument on command line"); goto ERROR; }
  *ret_go = go;
  return;
  
 ERROR:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  puts("\nwhere basic options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
  printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
  exit(1);  
}

static int
output_header(FILE *ofp, const ESL_GETOPTS *go, char *cmfile, char *seqfile)
{
  cm_banner(ofp, go->argv[0], banner);
  
                                         fprintf(ofp, "# query CM file:                         %s\n", cmfile);
                                         fprintf(ofp, "# target sequence database:              %s\n", seqfile);
  if (esl_opt_IsUsed(go, "-o"))          fprintf(ofp, "# output directed to file:               %s\n",      esl_opt_GetString(go, "-o"));
  if (esl_opt_IsUsed(go, "-A"))          fprintf(ofp, "# MSA of all hits saved to file:         %s\n",      esl_opt_GetString(go, "-A"));
  if (esl_opt_IsUsed(go, "--tblout"))    fprintf(ofp, "# per-seq hits tabular output:           %s\n",      esl_opt_GetString(go, "--tblout"));
  if (esl_opt_IsUsed(go, "--acc"))       fprintf(ofp, "# prefer accessions over names:          yes\n");
  if (esl_opt_IsUsed(go, "--noali"))     fprintf(ofp, "# show alignments in output:             no\n");
  if (esl_opt_IsUsed(go, "--notextw"))   fprintf(ofp, "# max ASCII text line length:            unlimited\n");
  if (esl_opt_IsUsed(go, "--textw"))     fprintf(ofp, "# max ASCII text line length:            %d\n",             esl_opt_GetInteger(go, "--textw"));
  if (esl_opt_IsUsed(go, "--popen"))     fprintf(ofp, "# gap open probability:                  %f\n",             esl_opt_GetReal  (go, "--popen"));
  if (esl_opt_IsUsed(go, "--pextend"))   fprintf(ofp, "# gap extend probability:                %f\n",             esl_opt_GetReal  (go, "--pextend"));
  if (esl_opt_IsUsed(go, "--mxfile"))    fprintf(ofp, "# subst score matrix:                    %s\n",             esl_opt_GetString(go, "--mxfile"));
  if (esl_opt_IsUsed(go, "-E"))          fprintf(ofp, "# sequence reporting threshold:          E-value <= %g\n",  esl_opt_GetReal(go, "-E"));
  if (esl_opt_IsUsed(go, "-T"))          fprintf(ofp, "# sequence reporting threshold:          score >= %g\n",    esl_opt_GetReal(go, "-T"));
  if (esl_opt_IsUsed(go, "--incE"))      fprintf(ofp, "# sequence inclusion threshold:          E-value <= %g\n",  esl_opt_GetReal(go, "--incE"));
  if (esl_opt_IsUsed(go, "--incT"))      fprintf(ofp, "# sequence inclusion threshold:          score >= %g\n",    esl_opt_GetReal(go, "--incT"));
  if (esl_opt_IsUsed(go, "--cut_ga"))    fprintf(ofp, "# model-specific thresholding:           GA cutoffs\n");
  if (esl_opt_IsUsed(go, "--cut_nc"))    fprintf(ofp, "# model-specific thresholding:           NC cutoffs\n");
  if (esl_opt_IsUsed(go, "--cut_tc"))    fprintf(ofp, "# model-specific thresholding:           TC cutoffs\n");
  if (esl_opt_IsUsed(go, "--cyk"))       fprintf(ofp, "# use CYK for final search stage         on\n");
  if (esl_opt_IsUsed(go, "--fqdb"))      fprintf(ofp, "# QDBs (CYK filter stage)                on\n");
  if (esl_opt_IsUsed(go, "--fsums"))     fprintf(ofp, "# HMM bands from sums (filter)           on\n");
  if (esl_opt_IsUsed(go, "--qdb"))       fprintf(ofp, "# QDBs (final stage)                     on\n");
  if (esl_opt_IsUsed(go, "--nonbanded")) fprintf(ofp, "# No bands (final stage)                 on\n");
  if (esl_opt_IsUsed(go, "--sums"))      fprintf(ofp, "# HMM bands from sums (final)            on\n");
  if (esl_opt_IsUsed(go, "--max"))       fprintf(ofp, "# Max sensitivity mode:                  on [all heuristic filters off]\n");
  if (esl_opt_IsUsed(go, "--fZ"))        fprintf(ofp, "# Filters set as if DB size in Mb is:    %f\n", esl_opt_GetReal(go, "--fZ"));
  if (esl_opt_IsUsed(go, "--rfam"))      fprintf(ofp, "# Rfam pipeline mode:                    on\n");
  if (esl_opt_IsUsed(go, "--noenvdef"))  fprintf(ofp, "# Envelope definition prior to CM search:off\n");
  if (esl_opt_IsUsed(go, "--pad"))       fprintf(ofp, "# hit padding strategy:                  on\n");
  if (esl_opt_IsUsed(go, "--noF1"))      fprintf(ofp, "# HMM MSV filter:                        off\n");
  if (esl_opt_IsUsed(go, "--noF2"))      fprintf(ofp, "# HMM Vit filter:                        off\n");
  if (esl_opt_IsUsed(go, "--noF3"))      fprintf(ofp, "# HMM Fwd filter:                        off\n");
  if (esl_opt_IsUsed(go, "--noF4"))      fprintf(ofp, "# HMM glocal Fwd filter:                 off\n");
  if (esl_opt_IsUsed(go, "--nohmm"))     fprintf(ofp, "# HMM filters (MSV/bias/Vit/Fwd):        off\n");
  if (esl_opt_IsUsed(go, "--noF6"))      fprintf(ofp, "# CM CYK filter:                         off\n");
  if (esl_opt_IsUsed(go, "--doF1b"))     fprintf(ofp, "# HMM MSV biased comp filter:            on\n");
  if (esl_opt_IsUsed(go, "--noF2b"))     fprintf(ofp, "# HMM Vit biased comp filter:            off\n");
  if (esl_opt_IsUsed(go, "--noF3b"))     fprintf(ofp, "# HMM Fwd biased comp filter:            off\n");
  if (esl_opt_IsUsed(go, "--noF4b"))     fprintf(ofp, "# HMM gFwd biased comp filter:           off\n");
  if (esl_opt_IsUsed(go, "--noF5b"))     fprintf(ofp, "# HMM per-envelope biased comp filter:   off\n");
  if (esl_opt_IsUsed(go, "--F1"))        fprintf(ofp, "# HMM MSV filter P threshold:            <= %g\n", esl_opt_GetReal(go, "--F1"));
  if (esl_opt_IsUsed(go, "--F1b"))       fprintf(ofp, "# HMM MSV bias P threshold:              <= %g\n", esl_opt_GetReal(go, "--F1b"));
  if (esl_opt_IsUsed(go, "--F2"))        fprintf(ofp, "# HMM Vit filter P threshold:            <= %g\n", esl_opt_GetReal(go, "--F2"));
  if (esl_opt_IsUsed(go, "--F2b"))       fprintf(ofp, "# HMM Vit bias P threshold:              <= %g\n", esl_opt_GetReal(go, "--F2b"));
  if (esl_opt_IsUsed(go, "--F3"))        fprintf(ofp, "# HMM Fwd filter P threshold:            <= %g\n", esl_opt_GetReal(go, "--F3"));
  if (esl_opt_IsUsed(go, "--F3b"))       fprintf(ofp, "# HMM Fwd bias P threshold:              <= %g\n", esl_opt_GetReal(go, "--F3b"));
  if (esl_opt_IsUsed(go, "--F4"))        fprintf(ofp, "# HMM glocal Fwd filter P threshold:     <= %g\n", esl_opt_GetReal(go, "--F4"));
  if (esl_opt_IsUsed(go, "--F4b"))       fprintf(ofp, "# HMM glocal Fwd bias P threshold:       <= %g\n", esl_opt_GetReal(go, "--F4b"));
  if (esl_opt_IsUsed(go, "--F5"))        fprintf(ofp, "# HMM env defn filter P threshold:       <= %g\n", esl_opt_GetReal(go, "--F5"));
  if (esl_opt_IsUsed(go, "--F5b"))       fprintf(ofp, "# HMM env defn bias   P threshold:       <= %g\n", esl_opt_GetReal(go, "--F5b"));
  if (esl_opt_IsUsed(go, "--F6"))        fprintf(ofp, "# CM CYK filter P threshold:             <= %g\n", esl_opt_GetReal(go, "--F6"));
  if (esl_opt_IsUsed(go, "--rt1"))       fprintf(ofp, "# Domain definition rt1 parameter        %g\n", esl_opt_GetReal(go, "--rt1"));
  if (esl_opt_IsUsed(go, "--rt2"))       fprintf(ofp, "# Domain definition rt2 parameter        %g\n", esl_opt_GetReal(go, "--rt2"));
  if (esl_opt_IsUsed(go, "--rt3"))       fprintf(ofp, "# Domain definition rt3 parameter        %g\n", esl_opt_GetReal(go, "--rt3"));
  if (esl_opt_IsUsed(go, "--ns"))        fprintf(ofp, "# Number of envelope tracebacks sampled  %d\n", esl_opt_GetInteger(go, "--ns"));
  if (esl_opt_IsUsed(go, "--localenv"))  fprintf(ofp, "# Define envelopes in local mode           on\n");
  if (esl_opt_IsUsed(go, "--null2"))     fprintf(ofp, "# null2 bias corrections:                on\n");
  if (esl_opt_IsUsed(go, "--nonull3"))   fprintf(ofp, "# null3 bias corrections:                off\n");
  if (esl_opt_IsUsed(go, "--toponly"))   fprintf(ofp, "# search top-strand only:                on\n");
  if (esl_opt_IsUsed(go, "--bottomonly"))fprintf(ofp, "# search bottom-strand only:             on\n");

  if (esl_opt_IsUsed(go, "--time-F1"))   fprintf(ofp, "# abort after Stage 1 MSV (for timing)   on\n");
  if (esl_opt_IsUsed(go, "--time-F2"))   fprintf(ofp, "# abort after Stage 2 Vit (for timing)   on\n");
  if (esl_opt_IsUsed(go, "--time-F3"))   fprintf(ofp, "# abort after Stage 3 Fwd (for timing)   on\n");
  if (esl_opt_IsUsed(go, "--time-F4"))   fprintf(ofp, "# abort after Stage 4 gFwd (for timing)  on\n");
  if (esl_opt_IsUsed(go, "--time-F5"))   fprintf(ofp, "# abort after Stage 5 env defn (for timing) on\n");
  if (esl_opt_IsUsed(go, "--time-F6"))   fprintf(ofp, "# abort after Stage 6 CYK (for timing)   on\n");

  if (esl_opt_IsUsed(go, "-Z"))          fprintf(ofp, "# database size is set to:               %.1f Mb\n",    esl_opt_GetReal(go, "-Z"));
  if (esl_opt_IsUsed(go, "--seed"))  {
    if (esl_opt_GetInteger(go, "--seed") == 0) fprintf(ofp, "# random number seed:                    one-time arbitrary\n");
    else                                       fprintf(ofp, "# random number seed set to:             %d\n", esl_opt_GetInteger(go, "--seed"));
  }
//  if (esl_opt_IsUsed(go, "--qformat")) fprintf(ofp, "# query <seqfile> format asserted:       %s\n",     esl_opt_GetString(go, "--qformat"));
  if (esl_opt_IsUsed(go, "--tformat"))   fprintf(ofp, "# targ <seqfile> format asserted:        %s\n", esl_opt_GetString(go, "--tformat"));
  if (esl_opt_IsUsed(go, "--w_beta"))
                                         fprintf(ofp, "# window length beta value:              %g\n", esl_opt_GetReal(go, "--w_beta"));
  if (esl_opt_IsUsed(go, "--w_length") )
                                         fprintf(ofp, "# window length :                        %d\n", esl_opt_GetInteger(go, "--w_length"));
#ifdef HMMER_THREADS
  if (esl_opt_IsUsed(go, "--cpu"))       fprintf(ofp, "# number of worker threads:              %d\n", esl_opt_GetInteger(go, "--cpu"));  
#endif
  fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
  return eslOK;
}

int
main(int argc, char **argv)
{

  int              status   = eslOK;

  ESL_GETOPTS     *go  = NULL;    /* command line processing                 */
  struct cfg_s     cfg;         /* configuration data                      */

  /* Set processor specific flags */
  impl_Init();

  /* Initialize what we can in the config structure (without knowing the alphabet yet)
   */
  cfg.cmfile     = NULL;
  cfg.dbfile     = NULL;

  cfg.do_mpi     = FALSE;               /* this gets reset below, if we init MPI */
  cfg.nproc      = 0;                   /* this gets reset below, if we init MPI */
  cfg.my_rank    = 0;                   /* this gets reset below, if we init MPI */

  /* Initializations */
  init_ilogsum();
  FLogsumInit();
  p7_FLogsumInit();        /* we're going to use table-driven Logsum() approximations at times */
  process_commandline(argc, argv, &go, &cfg.cmfile, &cfg.dbfile);    

  status = serial_master(go, &cfg);

  esl_getopts_Destroy(go);

  return status;
}

/* serial_master()
 * The serial version of cmsearch.
 * For each query CM in <cmfile> search the database for hits.
 * 
 * A master can only return if it's successful. All errors are handled immediately and fatally with p7_Fail().
 */
static int
serial_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  FILE            *ofp      = stdout;            /* results output file (-o)                        */
  FILE            *afp      = NULL;              /* alignment output file (-A)                      */
  FILE            *tblfp    = NULL;              /* output stream for tabular per-seq (--tblout)    */

  CMFILE          *cmfp;		         /* open input CM file stream                       */
  CM_t            *cm      = NULL;               /* one CM query                                    */
  P7_HMM         **hmmA    = NULL;               /* HMMs, used for filtering                        */
  int              nhmm    = 0;                  /* number of HMMs used for filtering */

  int              dbformat = eslSQFILE_UNKNOWN; /* format of dbfile                                */
  ESL_SQFILE      *dbfp     = NULL;              /* open input sequence file                        */

  ESL_ALPHABET    *abc      = NULL;              /* digital alphabet                                */
  ESL_ALPHABET   **hmm_abcA = NULL;              /* digital alphabets for hmms                      */
  ESL_STOPWATCH   *w;
  int              textw    = 0;
  int              nquery   = 0;
  int              status   = eslOK;
  int              qhstatus = eslOK;
  int              sstatus  = eslOK;
  int              i, m, z, z_offset;
  int64_t          L;    
  
  int              ncpus    = 0;

  int              infocnt  = 0;
  WORKER_INFO     *info     = NULL;

#ifdef HMMER_THREADS
  ESL_SQ_BLOCK    *block    = NULL;
  ESL_THREADS     *threadObj= NULL;
  ESL_WORK_QUEUE  *queue    = NULL;
#endif
  char             errbuf[cmERRBUFSIZE];

  double window_beta = -1.0 ;
  int window_length  = -1;

  int             *fcyk_dmin  = NULL;
  int             *fcyk_dmax  = NULL;
  int             *final_dmin = NULL;
  int             *final_dmax = NULL;

  if (esl_opt_IsUsed(go, "--w_beta")) { if (  ( window_beta   = esl_opt_GetReal(go, "--w_beta") )  < 0 || window_beta > 1  ) esl_fatal("Invalid window-length beta value\n"); }
  if (esl_opt_IsUsed(go, "--w_length")) { if (( window_length = esl_opt_GetInteger(go, "--w_length")) < 4  ) esl_fatal("Invalid window length value\n"); }

  w = esl_stopwatch_Create();

  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");

  if (esl_opt_IsOn(go, "--tformat")) {
    dbformat = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--tformat"));
    if (dbformat == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized sequence database file format\n", esl_opt_GetString(go, "--tformat"));
  }

#if 0 
  /* Open the target sequence database */
  status = esl_sqfile_Open(cfg->dbfile, dbformat, p7_SEQDBENV, &dbfp);
  if      (status == eslENOTFOUND) esl_fatal("Failed to open target sequence database %s for reading\n",      cfg->dbfile);
  else if (status == eslEFORMAT)   esl_fatal("Target sequence database file %s is empty or misformatted\n",   cfg->dbfile);
  else if (status == eslEINVAL)    esl_fatal("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        esl_fatal("Unexpected error %d opening target sequence database file %s\n", status, cfg->dbfile);
#endif

  /* Open the target sequence database */
  status = esl_sqfile_Open(cfg->dbfile, dbformat, p7_SEQDBENV, &dbfp);

  /* Open the SSI index for retrieval */
  if (dbfp->data.ascii.do_gzip)  esl_fatal("Reading gzipped sequence files is not supported.");
  if (dbfp->data.ascii.do_stdin) esl_fatal("Reading sequence files from stdin is not supported.");
  status = esl_sqfile_OpenSSI(dbfp, NULL);
  if      (status == eslEFORMAT)   esl_fatal("SSI index for database file is in incorrect format\n");
  else if (status == eslERANGE)    esl_fatal("SSI index for database file is in 64-bit format and we can't read it\n");
  else if (status != eslOK)        esl_fatal("Failed to open SSI index for database file\n");

  /* Determine database size */
  if(esl_opt_IsUsed(go, "-Z")) { 
    cfg->Z = (int64_t) (esl_opt_GetReal(go, "-Z") * 1000000.);
  }
  else { /* step through sequence SSI file to get database size */
    cfg->Z = 0;
    for(i = 0; i < dbfp->data.ascii.ssi->nprimary; i++) { 
      esl_ssi_FindNumber(dbfp->data.ascii.ssi, i, NULL, NULL, NULL, &L, NULL);
      cfg->Z += L;
    }
  }
  printf("cfg->Z: %" PRId64 " residues\n", cfg->Z);

  /* Open the query CM file */
  if ((cmfp = CMFileOpen(cfg->cmfile, NULL)) == NULL)
    esl_fatal("Failed to open covariance model save file %s\n", cfg->cmfile);

  /* Open the results output files */
  if (esl_opt_IsOn(go, "-o"))          { if ((ofp      = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) p7_Fail("Failed to open output file %s for writing\n",    esl_opt_GetString(go, "-o")); }
  if (esl_opt_IsOn(go, "-A"))          { if ((afp      = fopen(esl_opt_GetString(go, "-A"), "w")) == NULL) p7_Fail("Failed to open alignment file %s for writing\n", esl_opt_GetString(go, "-A")); }
  if (esl_opt_IsOn(go, "--tblout"))    { if ((tblfp    = fopen(esl_opt_GetString(go, "--tblout"),    "w")) == NULL)  esl_fatal("Failed to open tabular per-seq output file %s for writing\n", esl_opt_GetString(go, "--tblout")); }

#ifdef HMMER_THREADS
  /* initialize thread data */
  if (esl_opt_IsOn(go, "--cpu")) ncpus = esl_opt_GetInteger(go, "--cpu");
  else                                   esl_threads_CPUCount(&ncpus);
  /* EPN TEMP: until I write cm_clone to copy CMs to different infos */
  /* ncpus = 1; */
  printf("NCPUS: %d\n", ncpus);
  if (ncpus > 0) {
      threadObj = esl_threads_Create(&pipeline_thread);
      queue = esl_workqueue_Create(ncpus * 2);
  }
#endif

  infocnt = (ncpus == 0) ? 1 : ncpus;
  ESL_ALLOC(info, sizeof(*info) * infocnt);

  /* <abc> is not known 'til first CM is read. Could be DNA or RNA*/
  qhstatus = CMFileRead(cmfp, errbuf, &abc, &cm);

  if (qhstatus == eslOK) {
    /* One-time initializations after alphabet <abc> becomes known */
    output_header(ofp, go, cfg->cmfile, cfg->dbfile);
    dbfp->abc = abc;
    
    for (i = 0; i < infocnt; ++i)    {
      info[i].pli   = NULL;
      info[i].th    = NULL;
      info[i].cm    = NULL;
      info[i].omA   = NULL;
      info[i].p7_evparamAA = NULL;
      info[i].bgA   = NULL;
      info[i].nhmm  = 0;
#ifdef HMMER_THREADS
      info[i].queue = queue;
#endif
    }

#ifdef HMMER_THREADS    
    for (i = 0; i < ncpus * 2; ++i) {
      block = esl_sq_CreateDigitalBlock(BLOCK_SIZE, abc);
      if (block == NULL)           esl_fatal("Failed to allocate sequence block");
      
      status = esl_workqueue_Init(queue, block);
      if (status != eslOK)          esl_fatal("Failed to add block to work queue");
    }
#endif
  }
  
  /* Outer loop: over each query CM in <cmfile>. */
  while (qhstatus == eslOK) {
    nquery++;
    esl_stopwatch_Start(w);

    /* Make sure we have E-value stats for both the CM and the p7, if not we can't run the pipeline */
    if(! (cm->flags & CMH_EXPTAIL_STATS)) cm_Fail("no E-value parameters were read for CM: %s\n", cm->name);
    if(! (cm->flags & CMH_MLP7_STATS))    cm_Fail("no plan7 HMM E-value parameters were read for CM: %s\n", cm->name);

    /* Determine QDBs, if necessary. By default, they aren't. */
    fcyk_dmin  = fcyk_dmax  = NULL; 
    final_dmin = final_dmax = NULL; 
    if(esl_opt_GetBoolean(go, "--fqdb")) { /* determine CYK filter QDBs */
      determine_qdbs(cm, esl_opt_GetReal(go, "--fbeta"), &(fcyk_dmin), &(fcyk_dmax));
    }
    if(esl_opt_GetBoolean(go, "--qdb")) { /* determine final stage QDBs */
      determine_qdbs(cm, esl_opt_GetReal(go, "--beta"), &(final_dmin), &(final_dmax));
    }

    /* seqfile may need to be rewound (multiquery mode) */
    if (nquery > 1) {
      if (! esl_sqfile_IsRewindable(dbfp))
	esl_fatal("Target sequence file %s isn't rewindable; can't search it with multiple queries", cfg->dbfile);
      esl_sqfile_Position(dbfp, 0);
    }

    fprintf(ofp, "Query:       %s  [CLEN=%d]\n", cm->name, cm->clen);
    if (cm->acc)  fprintf(ofp, "Accession:   %s\n", cm->acc);
    if (cm->desc) fprintf(ofp, "Description: %s\n", cm->desc);
    
    /* Convert HMMs to optimized models */
    /* TODO: simplify this code */
    /* by default, we filter only with the cm->nap7 HMMs written in the CM file */
    nhmm = cm->nap7;
    z_offset = 0;
    /* but, if --doml selected or there's no HMMs in the file, use a  ML P7 HMM built from the CM */
    if (esl_opt_GetBoolean(go, "--doml") || (cm->nap7 == 0)) { 
      nhmm = 1 + cm->nap7;
      z_offset = 1;
    }
    else if(esl_opt_GetBoolean(go, "--noadd")) { 
      /* or if --noadd is selected, only use a ML p7 HMM */
      nhmm = 1;
      z_offset = 1;
    }      

    ESL_ALLOC(hmmA,     sizeof(P7_HMM *) * nhmm);
    ESL_ALLOC(hmm_abcA, sizeof(ESL_ALPHABET *) * nhmm);

    /* use the ML p7 HMM if nec */
    if (esl_opt_GetBoolean(go, "--doml") || (cm->nap7 == 0)) { 
      hmmA[0] = cm->mlp7;
      hmm_abcA[0] = esl_alphabet_Create(cm->mlp7->abc->type);
    }
    /* copy the HMMs from the file into the CM data structure */
    if(! esl_opt_GetBoolean(go, "--noadd")) { 
      for(z = 0; z < cm->nap7; z++) { 
	hmmA[z+z_offset]     = cm->ap7A[z];
	hmm_abcA[z+z_offset] = esl_alphabet_Create(cm->ap7A[z]->abc->type);
      }
    }

    P7_PROFILE      **gmA      = NULL;
    P7_OPROFILE     **omA      = NULL;       /* optimized query profile                  */
    P7_BG           **bgA      = NULL;

    ESL_ALLOC(gmA, sizeof(P7_PROFILE *)  * nhmm);
    ESL_ALLOC(omA, sizeof(P7_OPROFILE *) * nhmm);
    ESL_ALLOC(bgA, sizeof(P7_BG *)  * nhmm);
    for(m = 0; m < nhmm; m++) { 
      gmA[m] = p7_profile_Create (hmmA[m]->M, abc);
      omA[m] = p7_oprofile_Create(hmmA[m]->M, abc);
      bgA[m] = p7_bg_Create(hmm_abcA[m]);
      p7_ProfileConfig(hmmA[m], bgA[m], gmA[m], 100, p7_LOCAL);  /* 100 is a dummy length for now; and MSVFilter requires local mode */
      p7_oprofile_Convert(gmA[m], omA[m]);                       /* <om> is now p7_LOCAL, multihit */
      /* now convert gm to glocal, to define envelopes in cm_pipeline() */
      p7_ProfileConfig(hmmA[m], bgA[m], gmA[m], 100, p7_GLOCAL);
    }

    /* Create processing pipeline and hit list */
    for (i = 0; i < infocnt; ++i) {
      if((status = CloneCMJustReadFromFile(cm, errbuf, &(info[i].cm))) != eslOK) esl_fatal(errbuf);
      info[i].th  = cm_tophits_Create();
      info[i].nhmm = nhmm;
      ESL_ALLOC(info[i].gmA, sizeof(P7_PROFILE *)  * nhmm);
      ESL_ALLOC(info[i].omA, sizeof(P7_OPROFILE *) * nhmm);
      ESL_ALLOC(info[i].bgA, sizeof(P7_BG *)  * nhmm);
      ESL_ALLOC(info[i].p7_evparamAA, sizeof(float *)  * nhmm);

      /* now configure each CM */
      if(! esl_opt_GetBoolean(go, "-g")) { 
	if(! esl_opt_GetBoolean(go, "--cp9gloc")) { 
	  info[i].cm->config_opts |= CM_CONFIG_LOCAL;
	  info[i].cm->config_opts |= CM_CONFIG_HMMLOCAL;
	  if(! esl_opt_GetBoolean(go, "--cp9noel")) info[i].cm->config_opts |= CM_CONFIG_HMMEL; 
	}
      }
      if((status = ConfigCM(info[i].cm, errbuf, 
			    (i == 0) ? TRUE : FALSE, /* calculate W ? */
			    NULL, NULL)) != eslOK) cm_Fail("Error configuring CM: %s\n", errbuf);
      if(i != 0) info[i].cm->W = info[0].cm->W;

      for(m = 0; m < nhmm; m++) { 
	info[i].gmA[m]  = p7_profile_Clone(gmA[m]);
	info[i].omA[m]  = p7_oprofile_Clone(omA[m]);
	info[i].bgA[m]  = p7_bg_Create(hmm_abcA[m]);
	ESL_ALLOC(info[i].p7_evparamAA[m], sizeof(float) * CM_p7_NEVPARAM);

	/* HERE HERE HERE TODO: SIMPLFIY THIS CODE: */

	if(m == 0 && ((esl_opt_GetBoolean(go, "--doml")) || (cm->nap7 == 0))) { 
	  esl_vec_FCopy(cm->mlp7_evparam, CM_p7_NEVPARAM, info[i].p7_evparamAA[m]); 
	}
	else { 
	  esl_vec_FCopy(cm->ap7_evparamAA[m], CM_p7_NEVPARAM, info[i].p7_evparamAA[m]); 
	}
      }
      info[i].pli = cm_pipeline_Create(go, omA[0]->M, 100, cfg->Z, p7_SEARCH_SEQS); /* L_hint = 100 is just a dummy for now */
      cm_pli_NewModel(info[i].pli, info[i].cm, fcyk_dmin, fcyk_dmax, final_dmin, final_dmax, info[i].omA, info[i].bgA, info[i].nhmm);
      info[i].pli->do_top = (esl_opt_GetBoolean(go, "--bottomonly")) ? FALSE : TRUE;
      info[i].pli->do_bot = (esl_opt_GetBoolean(go, "--toponly"))    ? FALSE : TRUE;

#ifdef HMMER_THREADS
      if (ncpus > 0) esl_threads_AddThread(threadObj, &info[i]);
#endif
    }

#ifdef HMMER_THREADS
    if (ncpus > 0)  sstatus = thread_loop(info, threadObj, queue, dbfp);
    else            sstatus = serial_loop(info, dbfp);
#else
      sstatus = serial_loop(info, dbfp);
#endif
      switch(sstatus) {
        case eslEFORMAT:
          esl_fatal("Parse failed (sequence file %s):\n%s\n",
                    dbfp->filename, esl_sqfile_GetErrorBuf(dbfp));
          break;
        case eslEOF:
          /* do nothing */
          break;
        default:
          esl_fatal("Unexpected error %d reading sequence file %s", sstatus, dbfp->filename);
      }

      /* need to re-compute e-values before merging (when list will be sorted) */
      double eff_dbsize = 0.;
      eff_dbsize = (double) ((((double) cfg->Z / (double) cm->stats->expAA[info[0].pli->final_cm_exp_mode][0]->dbsize) * 
			      ((double) cm->stats->expAA[info[0].pli->final_cm_exp_mode][0]->nrandhits)) + 0.5);
	
      for (i = 0; i < infocnt; ++i) { 
	/* TO DO: compute E-values differently if --hmm (p7_tophits_ComputeNhmmerEvalues?) */
	cm_tophits_ComputeEvalues(info[i].th, eff_dbsize);
      }

      /* merge the results of the search results */
      for (i = 1; i < infocnt; ++i) {
          cm_tophits_Merge(info[0].th, info[i].th);
          cm_pipeline_Merge(info[0].pli, info[i].pli);

          cm_pipeline_Destroy(info[i].pli, cm);
          cm_tophits_Destroy(info[i].th);
	  for(m = 0; m < info[i].nhmm; m++) { 
	    p7_oprofile_Destroy(info[i].omA[m]);
	    p7_profile_Destroy(info[i].gmA[m]);
	    p7_bg_Destroy(info[i].bgA[m]);
	    free(info[i].p7_evparamAA[m]);
	  }
      }

      /* Print the results.  */
      cm_tophits_Sort(info->th);
      cm_tophits_Threshold(info->th, info->pli);

      /* tally up total number of hits and target coverage */
      info->pli->n_output = info->pli->pos_output = 0;
      for (i = 0; i < info->th->N; i++) {
	if ((info->th->hit[i]->flags & CM_HIT_IS_REPORTED) || (info->th->hit[i]->flags & CM_HIT_IS_INCLUDED)) { 
	  info->pli->n_output++;
	  info->pli->pos_output += abs(info->th->hit[i]->stop - info->th->hit[i]->start) + 1;
	}
      }

      /* Output */
      cm_tophits_Targets(ofp, info->th, info->pli, textw);
      fprintf(ofp, "\n\n");

      if(info->pli->show_alignments) {
	cm_tophits_HitAlignments(ofp, info->th, info->pli, textw);
	fprintf(ofp, "\n\n");
      }

      if (tblfp != NULL) { 
	cm_tophits_TabularTargets(tblfp, cm->name, cm->acc, info->th, info->pli, (nquery == 1)); 
      }
  
      esl_stopwatch_Stop(w);
      cm_pli_Statistics(ofp, info->pli, w);
      fprintf(ofp, "//\n");

      /* Output the results in an MSA (-A option) */
      if (afp) {
	esl_fatal("-A option not yet implemented.");

#if 0
	ESL_MSA *msa = NULL;

          if (cm_tophits_Alignment(info->th, abc, NULL, NULL, 0, p7_DEFAULT, &msa) == eslOK) {
            if (textw > 0) esl_msa_Write(afp, msa, eslMSAFILE_STOCKHOLM);
            else           esl_msa_Write(afp, msa, eslMSAFILE_PFAM);
	    
            fprintf(ofp, "# Alignment of %d hits satisfying inclusion thresholds saved to: %s\n", msa->nseq, esl_opt_GetString(go, "-A"));
          }  else {
	    fprintf(ofp, "# No hits satisfy inclusion thresholds; no alignment saved\n");
          }

          esl_msa_Destroy(msa);
#endif
      }

      cm_pipeline_Destroy(info->pli, cm);
      cm_tophits_Destroy(info->th);

      for(m = 0; m < info->nhmm; m++) { 
	p7_oprofile_Destroy(info->omA[m]);
	p7_profile_Destroy(info->gmA[m]);
	p7_bg_Destroy(info->bgA[m]);
	free(info->p7_evparamAA[m]);

	p7_oprofile_Destroy(omA[m]);
	p7_profile_Destroy(gmA[m]);
	p7_bg_Destroy(bgA[m]);
      }

      FreeCM(cm);

      qhstatus = CMFileRead(cmfp, errbuf, &abc, &cm);
  } /* end outer loop over query CMs */

  switch(qhstatus) {
    case eslEFORMAT:
      cm_Fail("bad file format in CM file %s", cfg->cmfile);
      break;
    case eslEOF:
      /* do nothing */
      break;
    default:
      cm_Fail("Unexpected error (%d: %s) in reading CMs from %s", qhstatus, errbuf, cfg->cmfile);
  }
  
#ifdef HMMER_THREADS
  if (ncpus > 0) {
      esl_workqueue_Reset(queue);
      while (esl_workqueue_Remove(queue, (void **) &block) == eslOK) {
          esl_sq_DestroyBlock(block);
      }
      esl_workqueue_Destroy(queue);
      esl_threads_Destroy(threadObj);
  }
#endif

  free(info);

  CMFileClose(cmfp);
  esl_sqfile_Close(dbfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);

  if (ofp != stdout) fclose(ofp);
  if (afp)           fclose(afp);
  if (tblfp)         fclose(tblfp);
  //if (domtblfp)      fclose(domtblfp);

  return eslOK;

 ERROR:
  return eslFAIL;
}

/* TODO: MPI code needs to be added here */

static int
serial_loop(WORKER_INFO *info, ESL_SQFILE *dbfp)
{
  int       status;
  int       wstatus;
  int       prev_hit_cnt;
  ESL_SQ   *dbsq   =  esl_sq_CreateDigital(info->cm->abc);

  wstatus = esl_sqio_ReadWindow(dbfp, 0, CMSEARCH_MAX_RESIDUE_COUNT, dbsq);

  while (wstatus == eslOK ) {
    
    cm_pli_NewSeq(info->pli, info->cm, dbsq);
    info->pli->nres -= dbsq->C; /* to account for overlapping region of windows */
    
    if (info->pli->do_top) { 
      prev_hit_cnt = info->th->N;
      if((status = cm_Pipeline(info->pli, info->cm, info->omA, info->gmA, info->bgA, info->p7_evparamAA, info->nhmm, dbsq, info->th)) != eslOK) cm_Fail("cm_pipeline() failed unexpected with status code %d\n", status);
      cm_pipeline_Reuse(info->pli); /* prepare for next search */
      
      /* modify hit positions to account for the position of the window in the full sequence */
      update_hit_positions(info->th, prev_hit_cnt, dbsq->start-1);
    }

    /* reverse complement */
    if (info->pli->do_bot && dbsq->abc->complement != NULL) { 
      prev_hit_cnt = info->th->N;
      esl_sq_ReverseComplement(dbsq);
      if((status = cm_Pipeline(info->pli, info->cm, info->omA, info->gmA, info->bgA, info->p7_evparamAA, info->nhmm, dbsq, info->th)) != eslOK) cm_Fail("cm_pipeline() failed unexpected with status code %d\n", status);
      cm_pipeline_Reuse(info->pli); /* prepare for next search */
      
      /* modify hit positions to account for the position of the window in the full sequence */
      update_hit_positions(info->th, prev_hit_cnt, dbsq->start-1);

      info->pli->nres += dbsq->W; /* not sure why this is here (originally from nhmmer.c) - to account for W overlap? */
    }
    
    wstatus = esl_sqio_ReadWindow(dbfp, info->omA[0]->max_length, CMSEARCH_MAX_RESIDUE_COUNT, dbsq);
    if (wstatus == eslEOD) { /* no more left of this sequence ... move along to the next sequence. */
      info->pli->nseqs++;
      esl_sq_Reuse(dbsq);
      wstatus = esl_sqio_ReadWindow(dbfp, 0, CMSEARCH_MAX_RESIDUE_COUNT, dbsq);
    }
  }
  esl_sq_Destroy(dbsq);
 
  return wstatus;
}

#ifdef HMMER_THREADS
static int
thread_loop(WORKER_INFO *info, ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQFILE *dbfp)
{

  /*int      wstatus, wstatus_next;*/
  int  status  = eslOK;
  int  sstatus = eslOK;
  int  eofCount = 0;
  ESL_SQ_BLOCK *block;
  void         *newBlock;
  int i;

  esl_workqueue_Reset(queue);
  esl_threads_WaitForStart(obj);

  status = esl_workqueue_ReaderUpdate(queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue reader failed");
  ((ESL_SQ_BLOCK *)newBlock)->complete = TRUE;

  /* Main loop: */
  while (sstatus == eslOK ) {
      block = (ESL_SQ_BLOCK *) newBlock;

      /* reset block as an empty vessel, possibly keeping the first sq intact for reading in the next window */
      if (block->count > 0 && block->complete) { 
	esl_sq_Reuse(block->list);
      }
      for (i=1; i<block->count; i++) { 
	esl_sq_Reuse(block->list + i);
      }

      sstatus = esl_sqio_ReadBlock(dbfp, block, CMSEARCH_MAX_RESIDUE_COUNT, TRUE);

      info->pli->nseqs += block->count - (block->complete ? 0 : 1); /* if there's an incomplete sequence read into the block wait to count it until it's complete. */

      if (sstatus == eslEOF) {
          if (eofCount < esl_threads_GetWorkerCount(obj)) sstatus = eslOK;
          ++eofCount;
      }

      if (sstatus == eslOK) {
          status = esl_workqueue_ReaderUpdate(queue, block, &newBlock);
          if (status != eslOK) esl_fatal("Work queue reader failed");
      }

      /* newBlock needs all this information so the next ReadBlock call will know what to do */
      ((ESL_SQ_BLOCK *)newBlock)->complete = block->complete;
      if (!block->complete) {
	/* the final sequence on the block was a probably-incomplete window of the active sequence, 
	 * so prep the next block to read in the next window */
	esl_sq_Copy(block->list + (block->count - 1) , ((ESL_SQ_BLOCK *)newBlock)->list);
	((ESL_SQ_BLOCK *)newBlock)->list->C = info->omA[0]->max_length;
      }
  }

  status = esl_workqueue_ReaderUpdate(queue, block, NULL);
  if (status != eslOK) esl_fatal("Work queue reader failed");

  if (sstatus == eslEOF) {
    /* wait for all the threads to complete */
    esl_threads_WaitForFinish(obj);
    esl_workqueue_Complete(queue);  
  }

  return sstatus;
}

static void 
pipeline_thread(void *arg)
{
  int prev_hit_cnt;
  int i;
  int status;
  int workeridx;
  WORKER_INFO   *info;
  ESL_THREADS   *obj;

  ESL_SQ_BLOCK  *block = NULL;
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
  block = (ESL_SQ_BLOCK *) newBlock;
  while (block->count > 0) { 
    /* Main loop: */
    for (i = 0; i < block->count; ++i) { 
      ESL_SQ *dbsq = block->list + i;
	
      cm_pli_NewSeq(info->pli, info->cm, dbsq);
      
      info->pli->nres -= dbsq->C; /* to account for overlapping region of windows */
      
      if (info->pli->do_top) { 
	prev_hit_cnt = info->th->N;
	if((status = cm_Pipeline(info->pli, info->cm, info->omA, info->gmA, info->bgA, info->p7_evparamAA, info->nhmm, dbsq, info->th)) != eslOK) cm_Fail("cm_pipeline() failed unexpected with status code %d\n", status);
	cm_pipeline_Reuse(info->pli); /* prepare for next search */
	
	/* modify hit positions to account for the position of the window in the full sequence */
	update_hit_positions(info->th, prev_hit_cnt, dbsq->start-1);
      }
	
      /* reverse complement */
      if (info->pli->do_bot && dbsq->abc->complement != NULL) {
	prev_hit_cnt = info->th->N;
	esl_sq_ReverseComplement(dbsq);
	if((status = cm_Pipeline(info->pli, info->cm, info->omA, info->gmA, info->bgA, info->p7_evparamAA, info->nhmm, dbsq, info->th)) != eslOK) cm_Fail("cm_pipeline() failed unexpected with status code %d\n", status);
	cm_pipeline_Reuse(info->pli); /* prepare for next search */

	/* modify hit positions to account for the position of the window in the full sequence */
	update_hit_positions(info->th, prev_hit_cnt, dbsq->start-1);

	info->pli->nres += dbsq->W; /* not sure why this is here (originally from nhmmer.c) - to account for W overlap? */
      }
    }
  
    status = esl_workqueue_WorkerUpdate(info->queue, block, &newBlock);
    if (status != eslOK) esl_fatal("Work queue worker failed");
    
    block = (ESL_SQ_BLOCK *) newBlock;
  }
  
  status = esl_workqueue_WorkerUpdate(info->queue, block, NULL);
  if (status != eslOK) esl_fatal("Work queue worker failed");

  esl_threads_Finished(obj, workeridx);
  return;
}
#endif   /* HMMER_THREADS */


/* Function:  determine_qdbs
 * Synopsis:  Calculate QDBs for a CM
 * Incept:    EPN, Tue May 24 10:25:28 2011
 *
 * Purpose:   Calculate QDBs (dmin/dmax) for a CM given a beta value.
 *
 * Returns:   <eslOK> on success.
 */
int
determine_qdbs(CM_t *cm, double beta, int **ret_dmin, int **ret_dmax)
{
  int safe_W;
  int *dmin = NULL; 
  int *dmax = NULL;

  safe_W = cm->clen * 3;
  while(!(BandCalculationEngine(cm, safe_W, beta, FALSE, &dmin, &dmax, NULL, NULL))) { 
    free(dmin);
    free(dmax);
    safe_W *= 2;
    if(safe_W > (cm->clen * 1000)) cm_Fail("Unable to calculate QDBs");
  }

  *ret_dmin = dmin;
  *ret_dmax = dmax;

  return eslOK;
}

/* Function:  update_hit_positions
 * Synopsis:  Update sequence positions in hit lits.
 * Incept:    EPN, Wed May 25 09:12:52 2011
 *
 * Purpose:   For hits <istart..th->N-1> in a hit list, update their
 *            positions by adding <offset> to start, stop and
 *            ad->sqfrom, ad->sqto.  This is necessary when we've
 *            searched a chunk of subsequence that originated
 *            somewhere within a larger sequence.
 *
 * Returns:   <eslOK> on success.
 */
int
update_hit_positions(CM_TOPHITS *th, int istart, int64_t offset)
{ 
  if(offset == 0) return eslOK;

  int i;
  for (i=istart; i < th->N ; i++) {
    th->unsrt[i].start += offset;
    th->unsrt[i].stop  += offset;
    if(th->unsrt[i].ad != NULL) { 
      th->unsrt[i].ad->sqfrom += offset;
      th->unsrt[i].ad->sqto   += offset;
    }
  }
  return eslOK;
}
/*****************************************************************
 * @LICENSE@
 *****************************************************************/

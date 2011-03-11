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
#ifdef P7_IMPL_DUMMY_INCLUDED
#define NHMMER_MAX_RESIDUE_COUNT (1024 * 100)
#else
#define NHMMER_MAX_RESIDUE_COUNT MAX_RESIDUE_COUNT   /*from esl_sqio_(ascii|ncbi).c*/
#endif

typedef struct {
#if 0
  /*#ifdef HMMER_THREADS*/
  ESL_WORK_QUEUE   *queue;
#endif /*HMMER_THREADS*/
  CM_PIPELINE      *pli;         /* work pipeline                           */
  P7_TOPHITS       *th;          /* top hit results                         */
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
#define XFASTMIDMAXOPTS "--fast,--mid,--max,--F1,--F2,--F3,--dF3,--F4,--nomsv,--novit,--nofwd,--nocyk,--nohmm,--hmm"
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
  { "--fast",       eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  XFASTMIDMAXOPTS, "set heuristic filters at strict-level (more speed, less power)", 7 },
  { "--mid",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  XFASTMIDMAXOPTS, "set heuristic filters at mid-level (mid-speed, mid-power)",    7 },
  { "--max",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  XFASTMIDMAXOPTS, "turn all heuristic filters off  (less speed, more power)",     7 },
  { "--noddef",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  "--hmm",         "do not define domains inside windows prior to CM stages",      7 },
  { "--msvmerge",   eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  "--nomsv,--nohmm","merge MSV windows prior to Viterbi filter",                   7 },
  { "--pad",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  "--hmm,--noddef","pad domains i..j to j-W+1..i+W-1",                             7 },
  { "--nomsv",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--max",          "skip the MSV filter stage",                                    7 },
  { "--shortmsv",   eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--max,--nomsv",  "run MSV on short 2*W subseqs, not longer subseqs",             7 },
  { "--novit",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--max",          "skip the Viterbi filter stage",                                7 },
  { "--nofwd",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--max",          "skip the Forward filter stage",                                7 },
  { "--nogfwd",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--max",          "skip the glocal Forward filter stage",                                7 },
  { "--nohmm",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--max,--hmm",    "skip all HMM filter stages (MSV/bias/Vit/Fwd)",                7 },
  { "--nocyk",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--hmm",          "skip the CYK filter stage",                                    7 },
  { "--domsvbias",  eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--max,--nomsv",  "turn on the MSV composition bias filter",                      7 },
  { "--novitbias",  eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--max,--novit",  "turn on the Vit composition bias filter",                     7 },
  { "--nofwdbias",  eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--max,--nofwd",  "turn on the Fwd composition bias filter",                     7 },
  { "--nogfwdbias", eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--max,--nofwd",  "turn on the glocal Fwd composition bias filter",              7 },
  { "--nodombias",  eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--max,--noddef", "turn on the per-domain composition bias filter",              7 },
  { "--donF3",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--max,--nofwd",  "filter based on glocal Forward score prior to doing ddef",     7 },
  { "--F1",         eslARG_REAL,  "0.40", NULL, NULL,    NULL,  NULL, "--max",          "Stage 1 (MSV) threshold: promote hits w/ P <= F1",             7 },
  { "--F1b",        eslARG_REAL,  "0.40", NULL, NULL,    NULL, "--domsvbias", "--max",  "Stage 1 (MSV) bias threshold: promote hits w/ P <= F1b",       7 },
  { "--F2",         eslARG_REAL,  "0.15", NULL, NULL,    NULL,  NULL, "--max",          "Stage 2 (Vit) threshold: promote hits w/ P <= F2",             7 },
  { "--F2b",        eslARG_REAL,  "0.15", NULL, NULL,    NULL, "--dovitbias",  "--max",  "Stage 2 (Vit) bias threshold: promote hits w/ P <= F2b",       7 },
  { "--F3",         eslARG_REAL,  "0.02", NULL, NULL,    NULL,  NULL, "--max",          "Stage 3 (Fwd) threshold: promote hits w/ P <= F3",             7 },
  { "--F3b",        eslARG_REAL,  "0.02", NULL, NULL,    NULL,"--dofwdbias",  "--max",  "Stage 3 (Fwd) bias threshold: promote hits w/ P <= F3b",       7 },
  { "--gF3",        eslARG_REAL,  "0.02", NULL, NULL,    NULL,  NULL, "--max",          "Stage 3 (gFwd) glocal threshold: promote hits w/ P <= gF3", 7 },
  { "--gF3b",       eslARG_REAL,  "0.02", NULL, NULL,    NULL,"--dogfwdbias",  "--max",  "Stage 3 (gFwd) glocal bias thr: promote hits w/ P <= gF3b", 7 },
  { "--dF3",        eslARG_REAL,  "0.02", NULL, NULL,    NULL,  NULL, "--max",          "Stage 3 (Fwd) per-domain threshold: promote hits w/ P <= dF3", 7 },
  { "--dF3b",       eslARG_REAL,  "0.02", NULL, NULL,    NULL,  NULL, "--nodombias,--max",  "Stage 3 (Fwd) per domain bias thr: promote hits w/ P <= dF3b", 7 },
  { "--F4",         eslARG_REAL,  "5e-4", NULL, NULL,    NULL,  NULL, "--max,--nocyk,--hmm","Stage 4 (CYK) threshold: promote hits w/ P <= F4",         7 },
  { "--E4",         eslARG_REAL,   NULL,  NULL, NULL,    NULL,  NULL, "--max,--nocyk,--hmm,--F4","Stage 4 (CYK) threshold: promote hits w/ E <= F4",    7 },
  { "--time-F1",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, TIMINGOPTS,        "abort after Stage 1 MSV; for timings",                        7 },
  { "--envF4x",     eslARG_INT,     "10", NULL, "n>=1",  NULL,  NULL, "--max,--nocyk,--hmm","CYK envelope redefinition threshold multiplier, <n> * F4", 7 },
  { "--noenvF4",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--max,--nocyk,--hmm","Do not redefine envelopes after stage 4 based on CYK hits",7 },
  { "--time-F2",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, TIMINGOPTS,        "abort after Stage 2 Vit; for timings",                        7 },
  { "--time-F3",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, TIMINGOPTS,        "abort after Stage 3 Fwd; for timings",                        7 },
  { "--time-gF3",   eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, TIMINGOPTS,        "abort after Stage 3 glocal Fwd; for timings",                 7 },
  { "--time-dF3",   eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, TIMINGOPTS,        "abort after domain def; for timings",                         7 },
  { "--time-bfil",  eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, TIMINGOPTS,        "abort after bfil; for timings",                               7 },
  { "--time-F4",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, TIMINGOPTS,        "abort after Stage 4 CYK; for timings",                        7 },
  { "--rt1",        eslARG_REAL,  "0.25", NULL, NULL,    NULL,  NULL, "--nohmm,--noddef","set domain definition rt1 parameter as <x>",                  7 },
  { "--rt2",        eslARG_REAL,  "0.10", NULL, NULL,    NULL,  NULL, "--nohmm,--noddef","set domain definition rt2 parameter as <x>",                  7 },
  { "--rt3",        eslARG_REAL,  "0.20", NULL, NULL,    NULL,  NULL, "--nohmm,--noddef","set domain definition rt3 parameter as <x>",                  7 },
  { "--ns",         eslARG_INT,   "200",  NULL, NULL,    NULL,  NULL, "--nohmm,--noddef","set number of domain tracebacks to <n>",                      7 },
  { "--localdom",   eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--nohmm,--noddef","define domains with HMM in glocal (not local) mode",          7 },
  { "--wnosplit",   eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--nomsv",         "do not split windows after MSV stage", 7 },
  { "--wmult",      eslARG_REAL,   "3.0", NULL, NULL,    NULL,  NULL, "--wnosplit,--nomsv",     "scalar multiplier for flagging window to split (if --wsplit)", 7 },
  { "--wcorr",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, NULL,              "use window size correction for Vit/Fwd filters", 7 },
  { "--nocorr",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, NULL,              "use no  correction for domain definition", 7 },
  { "--oldcorr",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--nocorr",        "use old correction for domain definition", 7 },
  { "--domhitbias", eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--nodombias",     "calc domain bias for only the domain, not entire window", 7 },
  { "--dosfwdbias", eslARG_NONE,   FALSE, NULL, NULL,    NULL,"--dofwdbias", NULL,       "use traceback samplings for fwd bias calculation", 7 },
  { "--fbns",       eslARG_INT,    "50",  NULL, NULL,    NULL,"--dosfwdbias",NULL,       "sample <n> tracebacks for fwd bias calculation", 7 },
  { "--gmsv",       eslARG_NONE,   FALSE, NULL, NULL,    NULL,"--shortmsv,--localdom", NULL,        "use generic MSV", 7 },
  { "--nogreedy",   eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,             "do not resolve hits with greedy algorithm, use optimal one", 7 },
  { "--noml",       eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,             "do not filter with a ML p7 HMM", 7 },
  { "--noadd",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,"--noml",           "do not filter with any additional p7 HMMs in the file", 7 },
  { "--filcmW",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,             "use CM's window length for all HMM filters", 7 },
  { "--glen",       eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, NULL,              "use length dependent glocal p7 filter P-value thresholds", 7},
  { "--glN",        eslARG_INT,   "201",  NULL, NULL,    NULL,"--glen",NULL,             "minimum value to start len-dependent glocal threshold", 7},
  { "--glX",        eslARG_INT,   "500",  NULL, NULL,    NULL,"--glen",NULL,             "maximum value for len-dependent glocal threshold", 7},
  { "--glstep",     eslARG_INT,   "100",  NULL, NULL,    NULL,"--glen",NULL,             "for len-dependent glocal thr, step size for halving thr", 7},
  { "--envF3",      eslARG_NONE,   FALSE, NULL, NULL,    NULL, NULL,   NULL,             "redefine envelope after local forward stage", 7 },
  { "--eF3X",       eslARG_REAL,  "1.1",  NULL, NULL,    NULL, NULL,   NULL,             "max avg num of passes through model for F3 env redefn", 7 },
  { "--eF3ns",      eslARG_INT,     "50",  NULL, NULL,    NULL, NULL,   NULL,             "number of samples for F3 env redefn", 7 },
  { "--envgF3",     eslARG_NONE,   FALSE, NULL, NULL,    NULL, NULL,   NULL,             "redefine envelope after glocal forward stage", 7 },
  { "--egF3X",      eslARG_REAL,  "1.1",  NULL, NULL,    NULL, NULL,   NULL,             "max avg num of passes through model for gF3 env redefn", 7 },
  { "--egF3ns",     eslARG_INT,    "200",  NULL, NULL,    NULL, NULL,   NULL,             "number of samples for gF3 env redefn", 7 },
  { "--egF3S",      eslARG_NONE,   FALSE,  NULL, NULL,    NULL, NULL,   NULL,             "be strict in gF3 env redef", 7 },
  { "--tmp",   eslARG_NONE,   FALSE, NULL,NULL, NULL,    NULL, "--nohmm,--noddef","use generic local", 7 },
/* Other options */
  { "--null2",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "turn on biased composition score corrections",               12 },
  { "-Z",           eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL,  NULL,            "set database size in *Mb* to <x> for E-value calculations",   12 },
  { "--seed",       eslARG_INT,    "42",  NULL, "n>=0",  NULL,  NULL,  NULL,            "set RNG seed to <n> (if 0: one-time arbitrary seed)",         12 },
  { "--w_beta",     eslARG_REAL,   NULL,  NULL, NULL,    NULL,  NULL,  NULL,            "tail mass at which window length is determined",               12 },
  { "--w_length",   eslARG_INT,    NULL,  NULL, NULL,    NULL,  NULL,  NULL,            "window length ",                                              12 },
  /* Options taken from infernal 1.0.2 cmsearch */
  /* options for algorithm for final round of search */
  { "-g",             eslARG_NONE,    FALSE,     NULL, NULL,    NULL,        NULL,            NULL, "configure CM for glocal alignment [default: local]", 1 },
  { "--cyk",          eslARG_NONE,    FALSE,     NULL, NULL,    NULL,        NULL, "--nocyk,--hmm", "use scanning CM CYK algorithm", 20 },
  { "--hmm",          eslARG_NONE,    FALSE,     NULL, NULL,    NULL,        NULL,         "--cyk",  "do not use the CM, use only the HMM", 20},
  /* banded options for CYK filter round of searching */
  { "--fnoqdb",       eslARG_NONE,    FALSE,     NULL, NULL,    NULL,        NULL,"--hmm,--nocyk",  "do not use QDBs in CYK filter round", 20 },
  { "--fbeta",        eslARG_REAL,    "1e-9",    NULL, "0<x<1", NULL,        NULL,"--hmm,--nocyk",  "set tail loss prob for CYK filter QDB calculation to <x>", 20 },
  { "--fhbanded",     eslARG_NONE,    FALSE,     NULL, NULL,    NULL,        NULL,"--hmm,--nocyk",  "calculate and use HMM bands in CYK filter round of CM search", 20 },
  { "--ftau",         eslARG_REAL,    "5e-6",    NULL, "0<x<1", NULL,"--fhbanded","--hmm,--nocyk",  "set tail loss prob for --fhbanded to <x>", 20 },
  { "--fsums",        eslARG_NONE,    FALSE,     NULL, NULL,    NULL,"--fhbanded","--hmm,--nocyk",  "w/--fhbanded use posterior sums (widens bands)", 20 },
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

/* Not used, but retained because esl option-handling code errors if it isn't kept here.  Placed in group 99 so it doesn't print to help*/
  { "--domZ",       eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL,  NULL,            "Not used",   99 },
  { "--domE",       eslARG_REAL,  "10.0", NULL, "x>0",   NULL,  NULL,  DOMREPOPTS,      "Not used",   99 },
  { "--domT",       eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  DOMREPOPTS,      "Not used",   99 },
  { "--incdomE",    eslARG_REAL,  "0.01", NULL, "x>0",   NULL,  NULL,  INCDOMOPTS,      "Not used",    99 },
  { "--incdomT",    eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  INCDOMOPTS,      "Not used",      99 },
/* will eventually bring these back, but store in group 99 for now, so they don't print to help*/
  { "--qformat",    eslARG_STRING,  NULL, NULL, NULL,    NULL,  NULL,  NULL,            "assert query <seqfile> is in format <s>: no autodetection",   99 },
  { "--tformat",    eslARG_STRING,  NULL, NULL, NULL,    NULL,  NULL,  NULL,            "assert target <seqfile> is in format <s>>: no autodetection", 99 },


#if 0
  /*#ifdef HMMER_THREADS*/
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
  char            *dbfile;            /* target sequence database file                   */
  char            *cmfile;           /* query HMM file                                  */

  int              do_mpi;            /* TRUE if we're doing MPI parallelization         */
  int              nproc;             /* how many MPI processes, total                   */
  int              my_rank;           /* who am I, in 0..nproc-1                         */
};

static char usage[]  = "[options] <query cmfile> <target seqfile>";
static char banner[] = "search a sequence database with an RNA CM";

static int  serial_master(ESL_GETOPTS *go, struct cfg_s *cfg);
static int  serial_loop  (WORKER_INFO *info, ESL_SQFILE *dbfp);

/* TEMPORARY, *_cm_pipeline() copies of p7_tophits.c functions, 
 * modified only to accept CM_PIPELINE object instead of P7_PIPELINE object 
 */
static int p7_tophits_Threshold_cm_pipeline(P7_TOPHITS *th, CM_PIPELINE *pli);
static int p7_tophits_Targets_cm_pipeline(FILE *ofp, P7_TOPHITS *th, CM_PIPELINE *pli, int textw);
static int p7_tophits_Domains_cm_pipeline(FILE *ofp, P7_TOPHITS *th, CM_PIPELINE *pli, int textw);
static int p7_tophits_TabularTargets_cm_pipeline(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, CM_PIPELINE *pli, int show_header);
static int p7_tophits_ComputeCMEvalues(P7_TOPHITS *th, double eff_dbsize);


#if 0
/*#ifdef HMMER_THREADS*/
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
  if (esl_opt_IsUsed(go, "--fnoqdb"))    fprintf(ofp, "# QDBs (CYK filter stage)                off\n");
  if (esl_opt_IsUsed(go, "--fhbanded"))  fprintf(ofp, "# HMM bands (filter stage)               on\n");
  if (esl_opt_IsUsed(go, "--fsums"))     fprintf(ofp, "# HMM bands from sums (filter)           on\n");
  if (esl_opt_IsUsed(go, "--qdb"))       fprintf(ofp, "# QDBs (final stage)                     on\n");
  if (esl_opt_IsUsed(go, "--nonbanded")) fprintf(ofp, "# No bands (final stage)                 on\n");
  if (esl_opt_IsUsed(go, "--sums"))      fprintf(ofp, "# HMM bands from sums (final)            on\n");
  if (esl_opt_IsUsed(go, "--max"))       fprintf(ofp, "# Max sensitivity mode:                  on [all heuristic filters off]\n");
  if (esl_opt_IsUsed(go, "--mid"))       fprintf(ofp, "# Mid-level filtering mode:              on\n");
  if (esl_opt_IsUsed(go, "--fast"))      fprintf(ofp, "# Strict-level filtering mode:           on\n");
  if (esl_opt_IsUsed(go, "--noddef"))    fprintf(ofp, "# Domain definition prior to CM search:  off\n");
  if (esl_opt_IsUsed(go, "--pad"))       fprintf(ofp, "# Domain padding strategy:               on\n");
  if (esl_opt_IsUsed(go, "--nomsv"))     fprintf(ofp, "# MSV filter:                            off\n");
  if (esl_opt_IsUsed(go, "--novit"))     fprintf(ofp, "# Vit filter:                            off\n");
  if (esl_opt_IsUsed(go, "--nofwd"))     fprintf(ofp, "# Fwd filter:                            off\n");
  if (esl_opt_IsUsed(go, "--nohmm"))     fprintf(ofp, "# HMM filters (MSV/bias/Vit/Fwd):        off\n");
  if (esl_opt_IsUsed(go, "--nocyk"))     fprintf(ofp, "# CYK filter:                            off\n");
  if (esl_opt_IsUsed(go, "--domsvbias")) fprintf(ofp, "# MSV biased comp HMM filter:            on\n");
  if (esl_opt_IsUsed(go, "--novitbias")) fprintf(ofp, "# Vit biased comp HMM filter:            off\n");
  if (esl_opt_IsUsed(go, "--nofwdbias")) fprintf(ofp, "# Fwd biased comp HMM filter:            off\n");
  if (esl_opt_IsUsed(go, "--nogfwdbias"))fprintf(ofp, "# gFwd biased comp HMM filter:           off\n");
  if (esl_opt_IsUsed(go, "--nodombias")) fprintf(ofp, "# Per-domain biased comp HMM filter:     off\n");
  if (esl_opt_IsUsed(go, "--F1"))        fprintf(ofp, "# MSV filter P threshold:                <= %g\n", esl_opt_GetReal(go, "--F1"));
  if (esl_opt_IsUsed(go, "--F2"))        fprintf(ofp, "# Vit filter P threshold:                <= %g\n", esl_opt_GetReal(go, "--F2"));
  if (esl_opt_IsUsed(go, "--F3"))        fprintf(ofp, "# Fwd filter P threshold:                <= %g\n", esl_opt_GetReal(go, "--F3"));
  if (esl_opt_IsUsed(go, "--F4"))        fprintf(ofp, "# CYK filter P threshold:                <= %g\n", esl_opt_GetReal(go, "--F4"));
  if (esl_opt_IsUsed(go, "--E4"))        fprintf(ofp, "# CYK filter E threshold:                <= %g\n", esl_opt_GetReal(go, "--E4"));
  if (esl_opt_IsUsed(go, "--F1b"))       fprintf(ofp, "# MSV bias P threshold:                  <= %g\n", esl_opt_GetReal(go, "--F1b"));
  if (esl_opt_IsUsed(go, "--F2b"))       fprintf(ofp, "# Vit bias P threshold:                  <= %g\n", esl_opt_GetReal(go, "--F2b"));
  if (esl_opt_IsUsed(go, "--F3b"))       fprintf(ofp, "# Fwd bias P threshold:                  <= %g\n", esl_opt_GetReal(go, "--F3b"));
  if (esl_opt_IsUsed(go, "--dF3b"))      fprintf(ofp, "# Domain defn P threshold:               <= %g\n", esl_opt_GetReal(go, "--dF3b"));
  if (esl_opt_IsUsed(go, "--rt1"))       fprintf(ofp, "# Domain definition rt1 parameter        %g\n", esl_opt_GetReal(go, "--rt1"));
  if (esl_opt_IsUsed(go, "--rt2"))       fprintf(ofp, "# Domain definition rt2 parameter        %g\n", esl_opt_GetReal(go, "--rt2"));
  if (esl_opt_IsUsed(go, "--rt3"))       fprintf(ofp, "# Domain definition rt3 parameter        %g\n", esl_opt_GetReal(go, "--rt3"));
  if (esl_opt_IsUsed(go, "--ns"))        fprintf(ofp, "# Number of domain tracebacks sampled    %d\n", esl_opt_GetInteger(go, "--ns"));
  if (esl_opt_IsUsed(go, "--localdom"))  fprintf(ofp, "# Define domains in local mode           on\n");
  if (esl_opt_IsUsed(go, "--null2"))     fprintf(ofp, "# null2 bias corrections:                on\n");
  if (esl_opt_IsUsed(go, "--nonull3"))   fprintf(ofp, "# null3 bias corrections:                off\n");
  if (esl_opt_IsUsed(go, "--toponly"))   fprintf(ofp, "# search top-strand only:                on\n");
  if (esl_opt_IsUsed(go, "--bottomonly"))fprintf(ofp, "# search bottom-strand only:             on\n");

  if (esl_opt_IsUsed(go, "--time-F1"))   fprintf(ofp, "# abort after Stage 1 MSV (for timing)   on\n");
  if (esl_opt_IsUsed(go, "--time-F2"))   fprintf(ofp, "# abort after Stage 2 Vit (for timing)   on\n");
  if (esl_opt_IsUsed(go, "--time-F3"))   fprintf(ofp, "# abort after Stage 3 Fwd (for timing)   on\n");
  if (esl_opt_IsUsed(go, "--time-dF3"))  fprintf(ofp, "# abort after domain defn (for timing)   on\n");
  if (esl_opt_IsUsed(go, "--time-bfil")) fprintf(ofp, "# abort after bfil calc   (for timing)   on\n");
  if (esl_opt_IsUsed(go, "--time-F4"))   fprintf(ofp, "# abort after Stage 4 CYK (for timing)   on\n");

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
#if 0
  /*#ifdef HMMER_THREADS*/
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
  
  int              ncpus    = 0;

  int              infocnt  = 0;
  WORKER_INFO     *info     = NULL;
#if 0
  /*#ifdef HMMER_THREADS*/
  ESL_SQ_BLOCK    *block    = NULL;
  ESL_THREADS     *threadObj= NULL;
  ESL_WORK_QUEUE  *queue    = NULL;
#endif
  char             errbuf[cmERRBUFSIZE];

  double window_beta = -1.0 ;
  int window_length  = -1;
  if (esl_opt_IsUsed(go, "--w_beta")) { if (  ( window_beta   = esl_opt_GetReal(go, "--w_beta") )  < 0 || window_beta > 1  ) esl_fatal("Invalid window-length beta value\n"); }
  if (esl_opt_IsUsed(go, "--w_length")) { if (( window_length = esl_opt_GetInteger(go, "--w_length")) < 4  ) esl_fatal("Invalid window length value\n"); }

  w = esl_stopwatch_Create();

  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");

  if (esl_opt_IsOn(go, "--tformat")) {
    dbformat = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--tformat"));
    if (dbformat == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized sequence database file format\n", esl_opt_GetString(go, "--tformat"));
  }

  /* Open the target sequence database */
  status = esl_sqfile_Open(cfg->dbfile, dbformat, p7_SEQDBENV, &dbfp);
  if      (status == eslENOTFOUND) esl_fatal("Failed to open target sequence database %s for reading\n",      cfg->dbfile);
  else if (status == eslEFORMAT)   esl_fatal("Target sequence database file %s is empty or misformatted\n",   cfg->dbfile);
  else if (status == eslEINVAL)    esl_fatal("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        esl_fatal("Unexpected error %d opening target sequence database file %s\n", status, cfg->dbfile);

  /* Open the query CM file */
  if ((cmfp = CMFileOpen(cfg->cmfile, NULL)) == NULL)
    esl_fatal("Failed to open covariance model save file %s\n", cfg->cmfile);

  /* Open the results output files */
  if (esl_opt_IsOn(go, "-o"))          { if ((ofp      = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) p7_Fail("Failed to open output file %s for writing\n",    esl_opt_GetString(go, "-o")); }
  if (esl_opt_IsOn(go, "-A"))          { if ((afp      = fopen(esl_opt_GetString(go, "-A"), "w")) == NULL) p7_Fail("Failed to open alignment file %s for writing\n", esl_opt_GetString(go, "-A")); }
  if (esl_opt_IsOn(go, "--tblout"))    { if ((tblfp    = fopen(esl_opt_GetString(go, "--tblout"),    "w")) == NULL)  esl_fatal("Failed to open tabular per-seq output file %s for writing\n", esl_opt_GetString(go, "--tblout")); }

#if 0
  /*#ifdef HMMER_THREADS*/
  /* initialize thread data */
  if (esl_opt_IsOn(go, "--cpu")) ncpus = esl_opt_GetInteger(go, "--cpu");
  else                                   esl_threads_CPUCount(&ncpus);
  /* EPN TEMP: until I write cm_clone to copy CMs to different infos */
  ncpus = 1;

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
#if 0
	  /*#ifdef HMMER_THREADS*/
          info[i].queue = queue;
#endif
      }

#if 0
      /*#ifdef HMMER_THREADS*/
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

    /* Make sure we have E-value stats for both the CM and the p7, if not we can't run the pipeline */
    if(! (cm->flags & CMH_EXPTAIL_STATS)) cm_Fail("no E-value parameters were read for CM: %s\n", cm->name);
    if(! (cm->flags & CMH_MLP7_STATS))    cm_Fail("no plan7 HMM E-value parameters were read for CM: %s\n", cm->name);

    P7_PROFILE      **gmA      = NULL;
    P7_OPROFILE     **omA      = NULL;       /* optimized query profile                  */
    P7_BG           **bgA      = NULL;
    int              safe_W;
    int             *fcyk_dmin  = NULL;
    int             *fcyk_dmax  = NULL;
    int             *final_dmin = NULL;
    int             *final_dmax = NULL;

    if(! esl_opt_GetBoolean(go, "-g")) { 
      if(! esl_opt_GetBoolean(go, "--cp9gloc")) { 
	cm->config_opts |= CM_CONFIG_LOCAL;
	cm->config_opts |= CM_CONFIG_HMMLOCAL;
	if(! esl_opt_GetBoolean(go, "--cp9noel")) cm->config_opts |= CM_CONFIG_HMMEL; 
      }
    }
    if((status = ConfigCM(cm, errbuf, 
			  TRUE, /* do calculate W */
			  NULL, NULL)) != eslOK) cm_Fail("Error configuring CM: %s\n", errbuf);

    /*********************************************/
    /* EPN TEMP: need to find a better spot for this code block, problem is I want dmin/dmax but I don't want them
     * stored in the CM data structure, just need the arrays so I can pass it to cm_pli_NewModel() and it will
     * create scan matrices for each set (fcyk_d{min,max}, final_d{min,max}) */

    /* calculate QDBs for filter and/or final rounds */
    if(esl_opt_GetBoolean(go, "--fnoqdb") || esl_opt_GetBoolean(go, "--nocyk")) { 
      fcyk_dmin = fcyk_dmax = NULL; 
    }
    else { 
      safe_W = cm->clen * 3;
      while(!(BandCalculationEngine(cm, safe_W, esl_opt_GetReal(go, "--fbeta"), FALSE, &(fcyk_dmin), &(fcyk_dmax), NULL, NULL))) { 
	free(fcyk_dmin);
	free(fcyk_dmax);
	safe_W *= 2;
	if(safe_W > (cm->clen * 1000)) cm_Fail("Unable to calculate QDBs");
      }
    }

    if(! esl_opt_GetBoolean(go, "--qdb")) { 
      final_dmin = final_dmax = NULL; 
    }
    else { 
      safe_W = cm->clen * 3;
      while(!(BandCalculationEngine(cm, safe_W, esl_opt_GetReal(go, "--beta"), FALSE, &(final_dmin), &(final_dmax), NULL, NULL))) { 
	free(final_dmin);
	free(final_dmax);
	safe_W *= 2;
	if(safe_W > (cm->clen * 1000)) cm_Fail("Unable to calculate QDBs");
      }
    }
    /*********************************************/

    nquery++;
    esl_stopwatch_Start(w);

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
    if(esl_opt_GetBoolean(go, "--noml")) { 
      if(cm->nap7 == 0) cm_Fail("--noml only makes sense if CM file has additional p7 HMMs, cm: %s does not\n", cm->name);
      nhmm = cm->nap7;
      z_offset = 0;
    }
    else if (! esl_opt_GetBoolean(go, "--noadd")) { 
      nhmm = 1 + cm->nap7;
      z_offset = 1;
    }
    else { 
      nhmm = 1;
    }
    ESL_ALLOC(hmmA,     sizeof(P7_HMM *) * nhmm);
    ESL_ALLOC(hmm_abcA, sizeof(ESL_ALPHABET *) * nhmm);
    printf("cm->nap7: %d\n", cm->nap7);
    if(! esl_opt_GetBoolean(go, "--noml")) { 
      hmmA[0] = cm->mlp7;
      hmm_abcA[0] = esl_alphabet_Create(cm->mlp7->abc->type);
    }
    if(! esl_opt_GetBoolean(go, "--noadd")) { 
      for(z = 0; z < cm->nap7; z++) { 
	hmmA[z+z_offset]     = cm->ap7A[z];
	hmm_abcA[z+z_offset] = esl_alphabet_Create(cm->ap7A[z]->abc->type);
      }
    }

    ESL_ALLOC(gmA, sizeof(P7_PROFILE *)  * nhmm);
    ESL_ALLOC(omA, sizeof(P7_OPROFILE *) * nhmm);
    ESL_ALLOC(bgA, sizeof(P7_BG *)  * nhmm);
    for(m = 0; m < nhmm; m++) { 
      gmA[m] = p7_profile_Create (hmmA[m]->M, abc);
      omA[m] = p7_oprofile_Create(hmmA[m]->M, abc);
      bgA[m] = p7_bg_Create(hmm_abcA[m]);
      p7_ProfileConfig(hmmA[m], bgA[m], gmA[m], 100, p7_LOCAL); /* 100 is a dummy length for now; and MSVFilter requires local mode */
      p7_oprofile_Convert(gmA[m], omA[m]);                        /* <om> is now p7_LOCAL, multihit */

      if((! esl_opt_GetBoolean(go, "--tmp")) && (! esl_opt_GetBoolean(go, "--gmsv"))) { 
	p7_ProfileConfig(hmmA[m], bgA[m], gmA[m], 100, p7_GLOCAL); /* this will be used to define domains in cm_pipeline() 
								      (we'll alternate b/t multi/uni later when processing domains) */
      }
    }

    for (i = 0; i < infocnt; ++i) {
      /* Create processing pipeline and hit list */
      info[i].th  = p7_tophits_Create();
      info[i].cm  = cm;
      info[i].nhmm = nhmm;
      ESL_ALLOC(info[i].gmA, sizeof(P7_PROFILE *)  * nhmm);
      ESL_ALLOC(info[i].omA, sizeof(P7_OPROFILE *) * nhmm);
      ESL_ALLOC(info[i].bgA, sizeof(P7_BG *)  * nhmm);
      ESL_ALLOC(info[i].p7_evparamAA, sizeof(float *)  * nhmm);
      for(m = 0; m < nhmm; m++) { 
	info[i].gmA[m]  = p7_profile_Clone(gmA[m]);
	info[i].omA[m]  = p7_oprofile_Clone(omA[m]);
	info[i].bgA[m]  = p7_bg_Create(hmm_abcA[m]);
	ESL_ALLOC(info[i].p7_evparamAA[m], sizeof(float) * CM_p7_NEVPARAM);

	if(m == 0 && (! esl_opt_GetBoolean(go, "--noml"))) { 
	  esl_vec_FCopy(cm->mlp7_evparam, CM_p7_NEVPARAM, info[i].p7_evparamAA[m]); 
	}
	else { 
	  if(esl_opt_GetBoolean(go, "--noml")) { 
	    esl_vec_FCopy(cm->ap7_evparamAA[m], CM_p7_NEVPARAM, info[i].p7_evparamAA[m]); 
	  }
	  else {
	    esl_vec_FCopy(cm->ap7_evparamAA[m-1], CM_p7_NEVPARAM, info[i].p7_evparamAA[m]); 
	  }
	}
      }
      info[i].pli = cm_pipeline_Create(go, omA[0]->M, 100, p7_SEARCH_SEQS); /* L_hint = 100 is just a dummy for now */
      cm_pli_NewModel(info[i].pli, cm, fcyk_dmin, fcyk_dmax, final_dmin, final_dmax, info[i].omA, info[i].bgA, info[i].nhmm);
      info[i].pli->do_top = (esl_opt_GetBoolean(go, "--bottomonly")) ? FALSE : TRUE;
      info[i].pli->do_bot = (esl_opt_GetBoolean(go, "--toponly"))    ? FALSE : TRUE;

#if 0
      /*#ifdef HMMER_THREADS*/
          if (ncpus > 0) esl_threads_AddThread(threadObj, &info[i]);
#endif
      }

#if 0
    /*#ifdef HMMER_THREADS*/
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

      //need to re-compute e-values before merging (when list will be sorted)
      double dbsize     = 0.;
      double eff_dbsize = 0.;
      if (esl_opt_IsUsed(go, "-Z")) { 
	dbsize     = 1000000*esl_opt_GetReal(go, "-Z");
      }
      else { 
	for (i = 0; i < infocnt; ++i) dbsize += info[i].pli->nres;
      }
      eff_dbsize = (double) (((dbsize / (double) cm->stats->expAA[info[0].pli->final_cm_exp_mode][0]->dbsize) * 
			      ((double) cm->stats->expAA[info[0].pli->final_cm_exp_mode][0]->nrandhits)) + 0.5);
	
      for (i = 0; i < infocnt; ++i) { 
	if(esl_opt_GetBoolean(go, "--hmm")) { p7_tophits_ComputeNhmmerEvalues(info[i].th, dbsize);     }
	else                                { p7_tophits_ComputeCMEvalues    (info[i].th, eff_dbsize); }
      }

      /* merge the results of the search results */
      for (i = 1; i < infocnt; ++i) {
          p7_tophits_Merge(info[0].th, info[i].th);
          cm_pipeline_Merge(info[0].pli, info[i].pli);

          cm_pipeline_Destroy(info[i].pli, cm);
          p7_tophits_Destroy(info[i].th);
	  for(m = 0; m < info[i].nhmm; m++) { 
	    p7_oprofile_Destroy(info[i].omA[m]);
	    p7_profile_Destroy(info[i].gmA[m]);
	    p7_bg_Destroy(info[i].bgA[m]);
	    free(info[i].p7_evparamAA[m]);
	  }
      }

      /* Print the results.  */
      p7_tophits_Sort(info->th);
      p7_tophits_RemoveDuplicates(info->th);

      p7_tophits_Threshold_cm_pipeline(info->th, info->pli);

      /* tally up total number of hits and target coverage */
      info->pli->n_output = info->pli->pos_output = 0;
      for (i = 0; i < info->th->N; i++) {
          if (info->th->hit[i]->dcl[0].is_reported || info->th->hit[i]->dcl[0].is_included) {
              info->pli->n_output++;
              info->pli->pos_output += abs(info->th->hit[i]->dcl[0].jali - info->th->hit[i]->dcl[0].iali) + 1;
          }
      }

      p7_tophits_Targets_cm_pipeline(ofp, info->th, info->pli, textw); fprintf(ofp, "\n\n");
      p7_tophits_Domains_cm_pipeline(ofp, info->th, info->pli, textw); fprintf(ofp, "\n\n");

      if (tblfp != NULL) { 
	p7_tophits_TabularTargets_cm_pipeline(tblfp, cm->name, cm->acc, info->th, info->pli, (nquery == 1)); 
      }
  
      esl_stopwatch_Stop(w);
      cm_pli_Statistics(ofp, info->pli, w);
      fprintf(ofp, "//\n");

      /* Output the results in an MSA (-A option) */
      if (afp) {
          ESL_MSA *msa = NULL;

          if (p7_tophits_Alignment(info->th, abc, NULL, NULL, 0, p7_DEFAULT, &msa) == eslOK) {
            if (textw > 0) esl_msa_Write(afp, msa, eslMSAFILE_STOCKHOLM);
            else           esl_msa_Write(afp, msa, eslMSAFILE_PFAM);
	    
            fprintf(ofp, "# Alignment of %d hits satisfying inclusion thresholds saved to: %s\n", msa->nseq, esl_opt_GetString(go, "-A"));
          }  else {
	    fprintf(ofp, "# No hits satisfy inclusion thresholds; no alignment saved\n");
          }

          esl_msa_Destroy(msa);
      }

      cm_pipeline_Destroy(info->pli, cm);
      p7_tophits_Destroy(info->th);

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
  
#if 0
  /*#ifdef HMMER_THREADS*/
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

//TODO: MPI code needs to be added here

static int
serial_loop(WORKER_INFO *info, ESL_SQFILE *dbfp)
{
  int       status;
  int      wstatus;
  int i;
  int prev_hit_cnt;
  ESL_SQ   *dbsq   =  esl_sq_CreateDigital(info->cm->abc);

  wstatus = esl_sqio_ReadWindow(dbfp, 0, NHMMER_MAX_RESIDUE_COUNT, dbsq);

  while (wstatus == eslOK ) {
    
    P7_DOMAIN *dcl;
    cm_pli_NewSeq(info->pli, info->cm, dbsq);
    info->pli->nres -= dbsq->C; // to account for overlapping region of windows
    
    if (info->pli->do_top) { 
      prev_hit_cnt = info->th->N;
      if((status = cm_Pipeline(info->pli, info->cm, info->omA, info->gmA, info->bgA, info->p7_evparamAA, info->nhmm, dbsq, info->th)) != eslOK) cm_Fail("cm_pipeline() failed unexpected with status code %d\n", status);
      cm_pipeline_Reuse(info->pli); // prepare for next search
      
      // modify hit positions to account for the position of the window in the full sequence
      for (i=prev_hit_cnt; i < info->th->N ; i++) {
	dcl = info->th->unsrt[i].dcl;
	dcl->ienv += dbsq->start - 1;
	dcl->jenv += dbsq->start - 1;
	dcl->iali += dbsq->start - 1;
	dcl->jali += dbsq->start - 1;
	dcl->ad->sqfrom += dbsq->start - 1;
	dcl->ad->sqto += dbsq->start - 1;
      }
    }
#ifdef eslAUGMENT_ALPHABET
    //reverse complement
    if (info->pli->do_bot && dbsq->abc->complement != NULL )
      {
	prev_hit_cnt = info->th->N;
	esl_sq_ReverseComplement(dbsq);
	if((status = cm_Pipeline(info->pli, info->cm, info->omA, info->gmA, info->bgA, info->p7_evparamAA, info->nhmm, dbsq, info->th)) != eslOK) cm_Fail("cm_pipeline() failed unexpected with status code %d\n", status);
	cm_pipeline_Reuse(info->pli); // prepare for next search
	  
	for (i=prev_hit_cnt; i < info->th->N ; i++) {
	  dcl = info->th->unsrt[i].dcl;
	  // modify hit positions to account for the position of the window in the full sequence
	  dcl->ienv = dbsq->start - dcl->ienv + 1;
	  dcl->jenv = dbsq->start - dcl->jenv + 1;
	  dcl->iali = dbsq->start - dcl->iali + 1;
	  dcl->jali = dbsq->start - dcl->jali + 1;
	  dcl->ad->sqfrom = dbsq->start - dcl->ad->sqfrom + 1;
	  dcl->ad->sqto = dbsq->start - dcl->ad->sqto + 1;
	}
	  
	info->pli->nres += dbsq->W;
      }
  }
#endif /*eslAUGMENT_ALPHABET*/

  wstatus = esl_sqio_ReadWindow(dbfp, info->omA[0]->max_length, NHMMER_MAX_RESIDUE_COUNT, dbsq);
  if (wstatus == eslEOD) { // no more left of this sequence ... move along to the next sequence.
    info->pli->nseqs++;
    esl_sq_Reuse(dbsq);
    wstatus = esl_sqio_ReadWindow(dbfp, 0, NHMMER_MAX_RESIDUE_COUNT, dbsq);
  }
  esl_sq_Destroy(dbsq);

  return wstatus;

}

#if 0
/*#ifdef HMMER_THREADS*/
static int
thread_loop(WORKER_INFO *info, ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQFILE *dbfp)
{

  //int      wstatus, wstatus_next;
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

      //reset block as an empty vessel, possibly keeping the first sq intact for reading in the next window
      if (block->count > 0 && block->complete)
          esl_sq_Reuse(block->list);
      for (i=1; i<block->count; i++)
          esl_sq_Reuse(block->list + i);

      sstatus = esl_sqio_ReadBlock(dbfp, block, NHMMER_MAX_RESIDUE_COUNT, TRUE);

      info->pli->nseqs += block->count - (block->complete ? 0 : 1);// if there's an incomplete sequence read into the block wait to count it until it's complete.

      if (sstatus == eslEOF) {
          if (eofCount < esl_threads_GetWorkerCount(obj)) sstatus = eslOK;
          ++eofCount;
      }

      if (sstatus == eslOK) {
          status = esl_workqueue_ReaderUpdate(queue, block, &newBlock);
          if (status != eslOK) esl_fatal("Work queue reader failed");
      }

      //newBlock needs all this information so the next ReadBlock call will know what to do
      ((ESL_SQ_BLOCK *)newBlock)->complete = block->complete;
      if (!block->complete) {
          // the final sequence on the block was a probably-incomplete window of the active sequence,
          // so prep the next block to read in the next window
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
  int i, j;
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
      P7_DOMAIN *dcl;
      
      info->pli->nres -= dbsq->C; // to account for overlapping region of windows
      
      if (info->pli->do_top) { 
	prev_hit_cnt = info->th->N;
	if((status = cm_Pipeline(info->pli, info->cm, info->omA, info->gmA, info->bgA, info->p7_evparamAA, info->nhmm, dbsq, info->th)) != eslOK) cm_Fail("cm_pipeline() failed unexpected with status code %d\n", status);
	cm_pipeline_Reuse(info->pli); // prepare for next search
	
	// modify hit positions to account for the position of the window in the full sequence
	for (j=prev_hit_cnt; j < info->th->N ; ++j) {
          dcl = info->th->unsrt[j].dcl;
          dcl->ienv += dbsq->start - 1;
          dcl->jenv += dbsq->start - 1;
          dcl->iali += dbsq->start - 1;
          dcl->jali += dbsq->start - 1;
          if(dcl->ad != NULL) { 
	    dcl->ad->sqfrom += dbsq->start - 1;
	    dcl->ad->sqto += dbsq->start - 1;
	  }
	}
      }
	
#ifdef eslAUGMENT_ALPHABET
      //reverse complement
      if (info->pli->do_bot && dbsq->abc->complement != NULL) {
	prev_hit_cnt = info->th->N;
	esl_sq_ReverseComplement(dbsq);
	if((status = cm_Pipeline(info->pli, info->cm, info->omA, info->gmA, info->bgA, info->p7_evparamAA, info->nhmm, dbsq, info->th)) != eslOK) cm_Fail("cm_pipeline() failed unexpected with status code %d\n", status);
	cm_pipeline_Reuse(info->pli); // prepare for next search
	
	for (j=prev_hit_cnt; j < info->th->N ; ++j) {
	  dcl = info->th->unsrt[j].dcl;
	  // modify hit positions to account for the position of the window in the full sequence
	  dcl->ienv = dbsq->start - dcl->ienv + 1;
	  dcl->jenv = dbsq->start - dcl->jenv + 1;
	  dcl->iali = dbsq->start - dcl->iali + 1;
	  dcl->jali = dbsq->start - dcl->jali + 1;
	  if(dcl->ad != NULL) { 
	    dcl->ad->sqfrom = dbsq->start - dcl->ad->sqfrom + 1;
	    dcl->ad->sqto = dbsq->start - dcl->ad->sqto + 1;
	  }
	}
	info->pli->nres += dbsq->W;
      }
#endif /*eslAUGMENT_ALPHABET*/
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



/*****************************************************************
 * @LICENSE@
 *****************************************************************/

/* Function:  p7_tophits_Threshold_cm_pipeline()
 * Synopsis:  Apply score and E-value thresholds to a hitlist before output.
 * Incept:    SRE, Tue Dec  9 09:04:55 2008 [Janelia]
 *
 * Purpose:   Identical to p7_pipeline.c::p7_tophits_Threshold() but takes 
 *            a CM_PIPELINE instead of a P7_PIPELINE argument, and 
 *            pli->long_targets is implicitly TRUE. 
 *         
 *            After a pipeline has completed, go through it and mark all
 *            the targets and domains that are "significant" (satisfying
 *            the reporting thresholds set for the pipeline). 
 *            
 *            Also sets the final total number of reported and
 *            included targets, the number of reported and included
 *            targets in each target, and the size of the search space
 *            for per-domain conditional E-value calculations,
 *            <pli->domZ>. By default, <pli->domZ> is the number of
 *            significant targets reported.
 *
 *            If model-specific thresholds were used in the pipeline,
 *            we cannot apply those thresholds now. They were already
 *            applied in the pipeline. In this case all we're
 *            responsible for here is counting them (setting
 *            nreported, nincluded counters).
 *            
 * Returns:   <eslOK> on success.
 */
int
p7_tophits_Threshold_cm_pipeline(P7_TOPHITS *th, CM_PIPELINE *pli)
{
  int h, d;    /* counters over sequence hits, domains in sequences */
  
  /* Flag reported, included targets (if we're using general thresholds) */
  if (! pli->use_bit_cutoffs) { 
    for (h = 0; h < th->N; h++) { 
      if (cm_pli_TargetReportable(pli, th->hit[h]->score, th->hit[h]->pvalue)) { 
	th->hit[h]->flags |= p7_IS_REPORTED;
	if (cm_pli_TargetIncludable(pli, th->hit[h]->score, th->hit[h]->pvalue))
	  th->hit[h]->flags |= p7_IS_INCLUDED;
      }
    }
  }

  /* Second pass is over domains, flagging reportable/includable ones. 
   * Note how this enforces a hierarchical logic of 
   * (sequence|domain) must be reported to be included, and
   * domain can only be (reported|included) if whole sequence is too.
   */
  if (! pli->use_bit_cutoffs) {
    for (h = 0; h < th->N; h++) {
      if (th->hit[h]->flags & p7_IS_REPORTED) {
	for (d = 0; d < th->hit[h]->ndom; d++) {
          if (cm_pli_TargetReportable(pli, th->hit[h]->dcl[d].bitscore, th->hit[h]->dcl[d].pvalue))
            th->hit[h]->dcl[d].is_reported = TRUE;
          if ((th->hit[h]->flags & p7_IS_INCLUDED) &&
              cm_pli_TargetIncludable(pli, th->hit[h]->dcl[d].bitscore, th->hit[h]->dcl[d].pvalue))
            th->hit[h]->dcl[d].is_included = TRUE;
	  /*printf("TH h: %3d  d: %3d  report? %s  include? %s\n", h, d, 
	    ((th->hit[h]->dcl[d].is_reported == TRUE) ? "TRUE" : "FALSE"),
	    ((th->hit[h]->dcl[d].is_included == TRUE) ? "TRUE" : "FALSE")); */
        }
      }
    }
  }

  /* Count reported, included targets */
  th->nreported = 0;
  th->nincluded = 0;
  for (h = 0; h < th->N; h++)
  {
      if (th->hit[h]->flags & p7_IS_REPORTED)  th->nreported++;
      if (th->hit[h]->flags & p7_IS_INCLUDED)  th->nincluded++;
  }

  /* Count the reported, included domains */
  for (h = 0; h < th->N; h++)  
    for (d = 0; d < th->hit[h]->ndom; d++)
      {
        if (th->hit[h]->dcl[d].is_reported) th->hit[h]->nreported++;
        if (th->hit[h]->dcl[d].is_included) th->hit[h]->nincluded++;
      }

  /* EPN: We can't call this b/c it's local to p7_pipeline.c */
  /*workaround_bug_h74(th);*/  /* blech. This function is defined above; see commentary and crossreferences there. */

  return eslOK;
}


/* Function:  p7_tophits_Targets_cm_pipeline()
 * Synopsis:  Standard output format for a top target hits list.
 * Incept:    SRE, Tue Dec  9 09:10:43 2008 [Janelia]
 *
 * Purpose:   Identical to p7_pipeline.c::p7_tophits_Targets() but takes 
 *            a CM_PIPELINE instead of a P7_PIPELINE argument, and 
 *            pli->long_targets is implicitly TRUE. 
 *
 *            Output a list of the reportable top target hits in <th> 
 *            in human-readable ASCII text format to stream <ofp>, using
 *            final pipeline accounting stored in <pli>. 
 * 
 *            The tophits list <th> should already be sorted (see
 *            <p7_tophits_Sort()> and thresholded (see
 *            <p7_tophits_Threshold>).
 *
 * Returns:   <eslOK> on success.
 */
int
p7_tophits_Targets_cm_pipeline(FILE *ofp, P7_TOPHITS *th, CM_PIPELINE *pli, int textw)
{
  char   newness;
  int    h;
  int    d;
  int    namew;
  int    posw;
  int    descw;
  char  *showname;
  int    have_printed_incthresh = FALSE;

  /* when --acc is on, we'll show accession if available, and fall back to name */
  if (pli->show_accessions) namew = ESL_MAX(8, p7_tophits_GetMaxShownLength(th));
  else                      namew = ESL_MAX(8, p7_tophits_GetMaxNameLength(th));

  if (textw >  0)           descw = ESL_MAX(32, textw - namew - 61); /* 61 chars excluding desc is from the format: 2 + 22+2 +22+2 +8+2 +<name>+1 */
  else                      descw = 0;                               /* unlimited desc length is handled separately */

  posw = ESL_MAX(6, p7_tophits_GetMaxPositionLength(th));

  fprintf(ofp, "Scores for complete hit%s:\n",     pli->mode == p7_SEARCH_SEQS ? "s" : "");
  /* The minimum width of the target table is 111 char: 47 from fields, 8 from min name, 32 from min desc, 13 spaces */
  fprintf(ofp, "  %9s %6s %5s  %-*s %*s %*s %s\n", "E-value", " score", " bias", namew, (pli->mode == p7_SEARCH_SEQS ? "Sequence":"Model"), posw, "start", posw, "end", "Description");
  fprintf(ofp, "  %9s %6s %5s  %-*s %*s %*s %s\n", "-------", "------", "-----", namew, "--------", posw, "-----", posw, "-----", "-----------");
  
  for (h = 0; h < th->N; h++)
    if (th->hit[h]->flags & p7_IS_REPORTED)
    {
        d    = th->hit[h]->best_domain;

        if (! (th->hit[h]->flags & p7_IS_INCLUDED) && ! have_printed_incthresh) {
          fprintf(ofp, "  ------ inclusion threshold ------\n");
          have_printed_incthresh = TRUE;
        }

        if (pli->show_accessions)
          {   /* the --acc option: report accessions rather than names if possible */
            if (th->hit[h]->acc != NULL && th->hit[h]->acc[0] != '\0') showname = th->hit[h]->acc;
            else                                                       showname = th->hit[h]->name;
          }
        else
          showname = th->hit[h]->name;

        if      (th->hit[h]->flags & p7_IS_NEW)     newness = '+';
        else if (th->hit[h]->flags & p7_IS_DROPPED) newness = '-';
        else                                        newness = ' ';
	fprintf(ofp, "%c %9.2g %6.1f %5.1f  %-*s %*d %*d ",
                newness,
                th->hit[h]->pvalue, // * pli->Z,
                th->hit[h]->score,
                eslCONST_LOG2R * th->hit[h]->dcl[d].dombias, // domain bias - seems like the right one to use, no?
                //th->hit[h]->pre_score - th->hit[h]->score, /* bias correction */
                namew, showname,
                posw, th->hit[h]->dcl[d].iali,
                posw, th->hit[h]->dcl[d].jali);

        if (textw > 0) fprintf(ofp, "%-.*s\n", descw, th->hit[h]->desc == NULL ? "" : th->hit[h]->desc);
        else           fprintf(ofp, "%s\n",           th->hit[h]->desc == NULL ? "" : th->hit[h]->desc);
        /* do NOT use *s with unlimited (INT_MAX) line length. Some systems
	 * have an fprintf() bug here (we found one on an Opteron/SUSE Linux
         * system (#h66)
         */
    }
  if (th->nreported == 0) fprintf(ofp, "\n   [No hits detected that satisfy reporting thresholds]\n");
  return eslOK;
}


/* Function:  p7_tophits_Domains_cm_pipeline()
 * Synopsis:  Standard output format for top domain hits and alignments.
 * Incept:    SRE, Tue Dec  9 09:32:32 2008 [Janelia]
 *
 * Purpose:   Identical to p7_pipeline.c::p7_tophits_Domains() but takes 
 *            a CM_PIPELINE instead of a P7_PIPELINE argument, and 
 *            pli->long_targets is implicitly TRUE. 
 *
 *            For each reportable target sequence, output a tabular summary
 *            of reportable domains found in it, followed by alignments of
 *            each domain.
 * 
 *            Similar to <p7_tophits_Targets()>; see additional notes there.
 */
int
p7_tophits_Domains_cm_pipeline(FILE *ofp, P7_TOPHITS *th, CM_PIPELINE *pli, int textw)
{
  int h, d;
  int nd;
  int namew, descw;
  char *showname;

  fprintf(ofp, "Annotation for each hit %s:\n",
	  pli->show_alignments ? " (and alignments)" : "");

  for (h = 0; h < th->N; h++)
    if (th->hit[h]->flags & p7_IS_REPORTED)
    {
        if (pli->show_accessions && th->hit[h]->acc != NULL && th->hit[h]->acc[0] != '\0')
        {
            showname = th->hit[h]->acc;
            namew    = strlen(th->hit[h]->acc);
        }
        else
        {
            showname = th->hit[h]->name;
            namew = strlen(th->hit[h]->name);
        }

        if (textw > 0)
        {
          descw = ESL_MAX(32, textw - namew - 5);
          fprintf(ofp, ">> %s  %-.*s\n", showname, descw, (th->hit[h]->desc == NULL ? "" : th->hit[h]->desc));
        }
        else
        {
          fprintf(ofp, ">> %s  %s\n",    showname,        (th->hit[h]->desc == NULL ? "" : th->hit[h]->desc));
        }

        if (th->hit[h]->nreported == 0)
        {
          fprintf(ofp,"   [No individual domains that satisfy reporting thresholds (although complete target did)]\n\n");
          continue;
        }

        /* The domain table is 101 char wide:
              #     score  bias    Evalue hmmfrom   hmmto    alifrom  ali to    envfrom  env to     acc
             ---   ------ ----- --------- ------- -------    ------- -------    ------- -------    ----
               1 ?  123.4  23.1    6.8e-9       3    1230 ..       1     492 []       2     490 .] 0.90
             123 ! 1234.5 123.4 123456789 1234567 1234567 .. 1234567 1234567 [] 1234567 1234568 .] 0.12
        */

	fprintf(ofp, "   %6s %5s %9s %7s %7s %2s %7s %7s %2s %7s %7s %2s %4s\n",  "score",  "bias",  "  Evalue", "hmmfrom",  "hmm to", "  ", "alifrom",  "ali to", "  ", "envfrom",  "env to", "  ",  "acc");
	fprintf(ofp, "   %6s %5s %9s %7s %7s %2s %7s %7s %2s %7s %7s %2s %4s\n",  "------", "-----", "---------", "-------", "-------", "  ", "-------", "-------", "  ", "-------", "-------", "  ", "----");

        nd = 0;
        for (d = 0; d < th->hit[h]->ndom; d++)
          if (th->hit[h]->dcl[d].is_reported)
            {
              nd++;
	      if(th->hit[h]->dcl[d].ad == NULL) { 
		fprintf(ofp, " %c %6.1f %5.1f %9.2g %7s %7s %c%c %7s %7s %c%c %7d %7d %c%c %4.2f\n",
			th->hit[h]->dcl[d].is_included ? '!' : '?',
			th->hit[h]->dcl[d].bitscore,
			th->hit[h]->dcl[d].dombias * eslCONST_LOG2R, /* convert NATS to BITS at last moment */
			th->hit[h]->dcl[d].pvalue,
			"-", "-", '-', '-', "-", "-", '-', '-', 
			th->hit[h]->dcl[d].ienv,
			th->hit[h]->dcl[d].jenv,
			'-', '-', 
			(th->hit[h]->dcl[d].oasc / (1.0 + fabs((float) (th->hit[h]->dcl[d].jenv - th->hit[h]->dcl[d].ienv)))));
	      }
	      else { 
		fprintf(ofp, " %c %6.1f %5.1f %9.2g %7d %7d %c%c %7ld %7ld %c%c %7d %7d %c%c %4.2f\n",
			th->hit[h]->dcl[d].is_included ? '!' : '?',
			th->hit[h]->dcl[d].bitscore,
			th->hit[h]->dcl[d].dombias * eslCONST_LOG2R, /* convert NATS to BITS at last moment */
			th->hit[h]->dcl[d].pvalue,
			th->hit[h]->dcl[d].ad->hmmfrom,
			th->hit[h]->dcl[d].ad->hmmto,
			(th->hit[h]->dcl[d].ad->hmmfrom == 1) ? '[' : '.',
			(th->hit[h]->dcl[d].ad->hmmto   == th->hit[h]->dcl[d].ad->M) ? ']' : '.',
			th->hit[h]->dcl[d].ad->sqfrom,
			th->hit[h]->dcl[d].ad->sqto,
			(th->hit[h]->dcl[d].ad->sqfrom == 1) ? '[' : '.',
			(th->hit[h]->dcl[d].ad->sqto   == th->hit[h]->dcl[d].ad->L) ? ']' : '.',
			th->hit[h]->dcl[d].ienv,
			th->hit[h]->dcl[d].jenv,
			(th->hit[h]->dcl[d].ienv == 1) ? '[' : '.',
			(th->hit[h]->dcl[d].jenv == th->hit[h]->dcl[d].ad->L) ? ']' : '.',
			(th->hit[h]->dcl[d].oasc / (1.0 + fabs((float) (th->hit[h]->dcl[d].jenv - th->hit[h]->dcl[d].ienv)))));
	      }
	    }

        if (pli->show_alignments) {
	  fprintf(ofp, "\n  Alignment:\n");
	  for (d = 0; d < th->hit[h]->ndom; d++) {
	    if (th->hit[h]->dcl[d].is_reported) { 
	      nd++;
	      fprintf(ofp, "  score: %.1f bits\n", th->hit[h]->dcl[d].bitscore);
              p7_alidisplay_Print(ofp, th->hit[h]->dcl[d].ad, 40, textw, pli->show_accessions);
              fprintf(ofp, "\n");
            }
	  }
	}
	else {
          fprintf(ofp, "\n");
	}
    }
  if (th->nreported == 0) { fprintf(ofp, "\n   [No targets detected that satisfy reporting thresholds]\n"); return eslOK; }
  return eslOK;
}


/* Function:  p7_tophits_TabularTargets_cm_pipeline()
 * Synopsis:  Output parsable table of per-sequence hits.
 * Incept:    SRE, Wed Mar 18 15:26:17 2009 [Janelia]
 *
 * Purpose:   Identical to p7_pipeline.c::p7_tophits_TabularTargets() but
 *            takes a CM_PIPELINE instead of a P7_PIPELINE argument, and
 *            pli->long_targets is implicitly TRUE.
 *
 *            Output a parseable table of reportable per-sequence hits
 *            in sorted tophits list <th> in an easily parsed ASCII
 *            tabular form to stream <ofp>, using final pipeline
 *            accounting stored in <pli>.
 *            
 *            Designed to be concatenated for multiple queries and
 *            multiple top hits list.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_tophits_TabularTargets_cm_pipeline(FILE *ofp, char *qname, char *qacc, P7_TOPHITS *th, CM_PIPELINE *pli, int show_header)
{
  int qnamew = ESL_MAX(20, strlen(qname));
  int tnamew = ESL_MAX(20, p7_tophits_GetMaxNameLength(th));
  int qaccw  = ((qacc != NULL) ? ESL_MAX(10, strlen(qacc)) : 10);
  int taccw  = ESL_MAX(10, p7_tophits_GetMaxAccessionLength(th));
  int posw;
  posw = ESL_MAX(7, p7_tophits_GetMaxPositionLength(th));

  int h,d;

  if (show_header)
  {
    fprintf(ofp, "#%-*s %-*s %-*s %-*s %s %s %*s %*s %9s %6s %5s %s\n",
	    tnamew-1, " target name",        taccw, "accession",  qnamew, "query name",           qaccw, "accession", "hmmfrom", "hmm to", posw, "alifrom", posw, "ali to",  "  E-value", " score", " bias", "description of target");
    fprintf(ofp, "#%*s %*s %*s %*s %s %s %*s %*s %9s %6s %5s %s\n",
	    tnamew-1, "-------------------", taccw, "----------", qnamew, "--------------------", qaccw, "----------", "-------", "-------", posw, "-------", posw, "-------", "---------", "------", "-----", "---------------------");
  }

  for (h = 0; h < th->N; h++)
    if (th->hit[h]->flags & p7_IS_REPORTED)    {
        d    = th->hit[h]->best_domain;
	if(th->hit[h]->dcl[d].ad != NULL) { 
	  fprintf(ofp, "%-*s %-*s %-*s %-*s %7d %7d %*d %*d %9.2g %6.1f %5.1f %s\n",
		  tnamew, th->hit[h]->name,
		  taccw,  th->hit[h]->acc ? th->hit[h]->acc : "-",
		  qnamew, qname,
		  qaccw,  ( (qacc != NULL && qacc[0] != '\0') ? qacc : "-"),
		  th->hit[h]->dcl[d].ad->hmmfrom,
		  th->hit[h]->dcl[d].ad->hmmto,
		  posw, th->hit[h]->dcl[d].iali,
		  posw, th->hit[h]->dcl[d].jali,
		  th->hit[h]->pvalue,
		  th->hit[h]->score,
		  th->hit[h]->dcl[d].dombias * eslCONST_LOG2R, /* convert NATS to BITS at last moment */
		  th->hit[h]->desc ?  th->hit[h]->desc : "-");
	}
	else { /* th->hit[h]->dcl[d].ad == NULL) */
	  fprintf(ofp, "%-*s %-*s %-*s %-*s %7s %7s %*d %*d %9.2g %6.1f %5.1f %s\n",
		  tnamew, th->hit[h]->name,
		  taccw,  th->hit[h]->acc ? th->hit[h]->acc : "-",
		  qnamew, qname,
		  qaccw,  ( (qacc != NULL && qacc[0] != '\0') ? qacc : "-"),
		  "-", "-", 
		  posw, th->hit[h]->dcl[d].iali,
		  posw, th->hit[h]->dcl[d].jali,
		  th->hit[h]->pvalue,
		  th->hit[h]->score,
		  th->hit[h]->dcl[d].dombias * eslCONST_LOG2R, /* convert NATS to BITS at last moment */
		  th->hit[h]->desc ?  th->hit[h]->desc : "-");
	}
    }
  
  return eslOK;
}

/* Function:  p7_tophits_ComputeCMEvalues()
 * Synopsis:  Compute CM E-values
 * Incept:    EPN, Tue Sep 28 05:26:20 2010
 *
 * Purpose:   After cmsearch pipeline has completed, the TopHits object contains
 *            objects with p-values that haven't yet been converted to e-values.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_tophits_ComputeCMEvalues(P7_TOPHITS *th, double eff_dbsize)
{
  int i; 

  for (i = 0; i < th->N ; i++) { 
    /*printf("CM P: %g  ",  th->unsrt[i].pvalue);*/
    th->unsrt[i].pvalue *= eff_dbsize;
    /*printf("E: %g  (eff_dbsize: %12.6f\n",  th->unsrt[i].pvalue, eff_dbsize);*/
    th->unsrt[i].sortkey = -th->unsrt[i].pvalue;
  }
  return eslOK;
}

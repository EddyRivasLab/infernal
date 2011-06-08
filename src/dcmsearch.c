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

/* set the max residue count to 1 meg when reading a block */
#define CMSEARCH_MAX_RESIDUE_COUNT 100000 /* differs from HMMER's default which is MAX_RESIDUE_COUNT from esl_sqio_(ascii|ncbi).c */

typedef struct {
#ifdef HMMER_THREADS
  ESL_WORK_QUEUE   *queue;
#endif /*HMMER_THREADS*/
  CM_PIPELINE      *pli;         /* work pipeline                           */
  CM_TOPHITS       *th;          /* top hit results                         */
  CM_t             *cm;          /* a covariance model                      */
  CMConsensus_t    *cmcons;      /* CM consensus info, for display purposes */
  P7_BG           **bgA;         /* null models                              */
  P7_OPROFILE     **omA;         /* optimized query profile HMMs            */
  P7_PROFILE      **gmA;         /* generic   query profile HMMs            */
  float           **p7_evparamAA;/* [0..nhmm-1][0..CM_p7_NEVPARAM] E-value parameters */
  int               nhmm;        /* number of HMM filters, size of omA, gmA, bgA */
  int               use_mlp7_hmm;/* TRUE if the CM's mlp7 will be used to filter */
  /* QDBs for CM, only necessary if --qdb or --fqdb are enabled, else they are NULL */
  int              *fcyk_dmin;
  int              *fcyk_dmax;
  int              *final_dmin;
  int              *final_dmax;

} WORKER_INFO;

#define REPOPTS     "-E,-T,--cut_ga,--cut_nc,--cut_tc"
#define DOMREPOPTS  "--domE,--domT,--cut_ga,--cut_nc,--cut_tc"
#define INCOPTS     "--incE,--incT,--cut_ga,--cut_nc,--cut_tc"
#define INCDOMOPTS  "--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
#define THRESHOPTS  "-E,-T,--domE,--domT,--incE,--incT,--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
#define XFASTMIDMAXOPTS "--max,--F1,--F2,--F3,--F4,--F5,--F6,--noF1,--noF2,--noF3,--nogfwd,--nocyk,--nohmm,--hmm"
#define TIMINGOPTS  "--time-F1,--time-F2,--time-F3,--time-dF3,--time-bfil,--time-F4"

#if defined (HMMER_THREADS) && defined (HAVE_MPI)
#define CPUOPTS     "--mpi"
#define MPIOPTS     "--cpu"
#else
#define CPUOPTS     NULL
#define MPIOPTS     NULL
#endif

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
#ifdef HAVE_MPI
  { "--stall",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,"--mpi", NULL,            "arrest after start: for debugging MPI under gdb",             12 },  
  { "--mpi",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  MPIOPTS,         "run as an MPI parallel program",                              12 },
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

/* Functions to avoid code duplication for common tasks */
static int   open_dbfile(ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, ESL_SQFILE **ret_dbfp);
static int   open_dbfile_ssi(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_SQFILE *dbfp, char *errbuf);
static int   determine_dbsize_using_ssi(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_SQFILE *dbfp, char *errbuf);
static WORKER_INFO *create_info();
static int   clone_info(ESL_GETOPTS *go, WORKER_INFO *src_info, WORKER_INFO *dest_infoA, int dest_infocnt, char *errbuf);
static void  free_info(WORKER_INFO *info);
static int   setup_cm(ESL_GETOPTS *go, WORKER_INFO *info, int do_calculate_W, char *errbuf);
static int   setup_qdbs(ESL_GETOPTS *go, WORKER_INFO *info, char *errbuf);
static int   setup_hmm_filters(ESL_GETOPTS *go, WORKER_INFO *info, const ESL_ALPHABET *abc, char *errbuf, P7_HMM ***ret_hmmA);
static int   update_hit_positions(CM_TOPHITS *th, int hit_start, int64_t seq_start, int in_revcomp);

#ifdef HMMER_THREADS
#define BLOCK_SIZE 1000

static int  thread_loop(WORKER_INFO *info, ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQFILE *dbfp);
static void pipeline_thread(void *arg);
#endif /*HMMER_THREADS*/

#ifdef HAVE_MPI
static int  mpi_master   (ESL_GETOPTS *go, struct cfg_s *cfg);
static int  mpi_worker   (ESL_GETOPTS *go, struct cfg_s *cfg);
#endif /*HAVE_MPI*/


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

    puts("\nOptions controlling alignment of hits:");
    esl_opt_DisplayHelp(stdout, go, 8, 2, 80); 

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

  /* Figure out who we are, and send control there: 
   * we might be an MPI master, an MPI worker, or a serial program.
   */
#ifdef HAVE_MPI

  /* TEMP */
  pid_t pid;
  /* get the process id */
  pid = getpid();
  //printf("The process id is %d\n", pid);
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
  P7_HMM         **hmmA    = NULL;               /* HMMs, used for filtering                        */
  ESL_SQFILE      *dbfp     = NULL;              /* open input sequence file                        */
  ESL_ALPHABET    *abc      = NULL;              /* digital alphabet                                */
  ESL_STOPWATCH   *w;
  int              textw    = 0;
  int              nquery   = 0;
  int              status   = eslOK;
  int              qhstatus = eslOK;
  int              sstatus  = eslOK;
  int              i;
  double           eff_dbsize;                   /* the effective database size, max # hits expected in db of this size */
  
  int              ncpus    = 0;

  WORKER_INFO     *tinfo;                        /* the template info, cloned to make the worked info */
  WORKER_INFO     *info          = NULL;         /* the worker info */
  int              infocnt       = 0;

#ifdef HMMER_THREADS
  ESL_SQ_BLOCK    *block    = NULL;
  ESL_THREADS     *threadObj= NULL;
  ESL_WORK_QUEUE  *queue    = NULL;
#endif
  char             errbuf[cmERRBUFSIZE];

  w = esl_stopwatch_Create();

  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");

  /* Open the database file */
  if((status = open_dbfile(go, cfg, errbuf, &dbfp)) != eslOK) esl_fatal(errbuf); 
  /* Determine database size */
  if(esl_opt_IsUsed(go, "-Z")) { /* enabling -Z is the only way to bypass the SSI index file requirement */
    cfg->Z = (int64_t) (esl_opt_GetReal(go, "-Z") * 1000000.); 
  }
  else { 
    if((status = open_dbfile_ssi           (go, cfg, dbfp, errbuf)) != eslOK) esl_fatal(errbuf); 
    if((status = determine_dbsize_using_ssi(go, cfg, dbfp, errbuf)) != eslOK) esl_fatal(errbuf); 
  }

  /* Open the query CM file */
  if ((cmfp = CMFileOpen(cfg->cmfile, NULL)) == NULL) cm_Fail("Failed to open covariance model save file %s\n", cfg->cmfile);

  /* Open the results output files */
  if (esl_opt_IsOn(go, "-o"))          { if ((ofp      = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) p7_Fail("Failed to open output file %s for writing\n",    esl_opt_GetString(go, "-o")); }
  if (esl_opt_IsOn(go, "-A"))          { if ((afp      = fopen(esl_opt_GetString(go, "-A"), "w")) == NULL) p7_Fail("Failed to open alignment file %s for writing\n", esl_opt_GetString(go, "-A")); }
  if (esl_opt_IsOn(go, "--tblout"))    { if ((tblfp    = fopen(esl_opt_GetString(go, "--tblout"),    "w")) == NULL)  esl_fatal("Failed to open tabular per-seq output file %s for writing\n", esl_opt_GetString(go, "--tblout")); }

#ifdef HMMER_THREADS
  /* initialize thread data */
  if (esl_opt_IsOn(go, "--cpu")) ncpus = esl_opt_GetInteger(go, "--cpu");
  else                                   esl_threads_CPUCount(&ncpus);
  printf("NCPUS: %d\n", ncpus);
  if (ncpus > 0) {
      threadObj = esl_threads_Create(&pipeline_thread);
      queue = esl_workqueue_Create(ncpus * 2);
  }
#endif

  infocnt = (ncpus == 0) ? 1 : ncpus;
  ESL_ALLOC(info, sizeof(WORKER_INFO) * infocnt);
  ESL_ALLOC(tinfo, sizeof(WORKER_INFO));

  /* <abc> is not known 'til first CM is read. Could be DNA or RNA*/
  qhstatus = CMFileRead(cmfp, errbuf, &abc, &(tinfo->cm));

  if (qhstatus == eslOK) {
    /* One-time initializations after alphabet <abc> becomes known */
    output_header(ofp, go, cfg->cmfile, cfg->dbfile);
    dbfp->abc = abc;
    
    for (i = 0; i < infocnt; ++i)    {
      info[i].pli          = NULL;
      info[i].th           = NULL;
      info[i].cm           = NULL;
      info[i].omA          = NULL;
      info[i].bgA          = NULL;
      info[i].p7_evparamAA = NULL;
      info[i].nhmm         = 0;
      info[i].use_mlp7_hmm = FALSE;
      info[i].fcyk_dmin    = NULL;
      info[i].fcyk_dmax    = NULL;
      info[i].final_dmin   = NULL;
      info[i].final_dmax   = NULL;
#ifdef HMMER_THREADS
      info[i].queue        = queue;
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

    tinfo->pli = NULL; 
    tinfo->th  = NULL;
    tinfo->cmcons = NULL;
    tinfo->fcyk_dmin = tinfo->fcyk_dmax = tinfo->final_dmin = tinfo->final_dmax = NULL;

    /* Make sure we have E-value stats for both the CM and the p7, if not we can't run the pipeline */
    if(! (tinfo->cm->flags & CMH_EXPTAIL_STATS)) cm_Fail("no E-value parameters were read for CM: %s\n", tinfo->cm->name);
    if(! (tinfo->cm->flags & CMH_MLP7_STATS))    cm_Fail("no plan7 HMM E-value parameters were read for CM: %s\n,", tinfo->cm->name);

    /* seqfile may need to be rewound (multiquery mode) */
    if (nquery > 1) {
      if (! esl_sqfile_IsRewindable(dbfp))
	esl_fatal("Target sequence file %s isn't rewindable; can't search it with multiple queries", cfg->dbfile);
      esl_sqfile_Position(dbfp, 0);
    }

    fprintf(ofp, "Query:       %s  [CLEN=%d]\n", tinfo->cm->name, tinfo->cm->clen);
    if (tinfo->cm->acc)  fprintf(ofp, "Accession:   %s\n", tinfo->cm->acc);
    if (tinfo->cm->desc) fprintf(ofp, "Description: %s\n", tinfo->cm->desc);
    
    /* Determine QDBs, if necessary. By default they are not needed. */
    if((status = setup_qdbs(go, tinfo, errbuf)) != eslOK) cm_Fail(errbuf);
    /* Setup HMM filters */
    if((status = setup_hmm_filters(go, tinfo, abc, errbuf, &hmmA)) != eslOK) cm_Fail(errbuf);
    /* Clone all data structures in tinfo into the WORKER_INFO's in info */
    if((status = clone_info(go, tinfo, info, infocnt, errbuf)) != eslOK) cm_Fail(errbuf);
    
    /* Create processing pipeline and hit list */
    tinfo[i].pli = NULL;
    tinfo[i].th  = NULL;
    for (i = 0; i < infocnt; ++i) {
      info[i].pli = cm_pipeline_Create(go, info[i].omA[0]->M, 100, cfg->Z, CM_SEARCH_SEQS); /* L_hint = 100 is just a dummy for now */
      cm_pli_NewModel(info[i].pli, info[i].cm, info[i].fcyk_dmin, info[i].fcyk_dmax, info[i].final_dmin, info[i].final_dmax, info[i].omA, info[i].bgA, info[i].nhmm);

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

      /* we need to re-compute e-values before merging (when list will be sorted) */
      eff_dbsize = (double) ((((double) cfg->Z / (double) info[0].cm->stats->expAA[info[0].pli->final_cm_exp_mode][0]->dbsize) * 
			      ((double) info[0].cm->stats->expAA[info[0].pli->final_cm_exp_mode][0]->nrandhits)) + 0.5);
      for (i = 0; i < infocnt; ++i) { 
	/* TO DO: compute E-values differently if --hmm (p7_tophits_ComputeNhmmerEvalues?) */
	cm_tophits_ComputeEvalues(info[i].th, eff_dbsize);
      }

      /* merge the results of the search results */
      for (i = 1; i < infocnt; ++i) {
          cm_tophits_Merge(info[0].th, info[i].th);
          cm_pipeline_Merge(info[0].pli, info[i].pli);
	  free_info(&info[i]);
      }

      /* Sort and enforce threshold */
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
	if((status = cm_tophits_HitAlignments(ofp, info->th, info->pli, textw)) != eslOK) esl_fatal("Out of memory");
	fprintf(ofp, "\n\n");
      }

      if (tblfp != NULL) { 
	cm_tophits_TabularTargets(tblfp, tinfo->cm->name, tinfo->cm->acc, info->th, info->pli, (nquery == 1)); 
      }
  
      esl_stopwatch_Stop(w);
      cm_pli_Statistics(ofp, info->pli, w);
      fprintf(ofp, "//\n");

      /* Output the results in an MSA (-A option) */
      if (afp) {
	esl_fatal("-A option not yet implemented.");

#if 0
	/* WHEN -A IS IMPLEMENTED, UNCOMMENT THIS BLOCK */
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
      free_info(tinfo);
      free_info(info);
      free(hmmA);
      
      qhstatus = CMFileRead(cmfp, errbuf, &abc, &(tinfo->cm));
  } /* end outer loop over query CMs */

  switch(qhstatus) {
  case eslEFORMAT:   cm_Fail("bad file format in CM file %s",             cfg->cmfile);      break;
  case eslEOF:       /* do nothing. EOF is what we want. */                                    break;
  default:           cm_Fail("Unexpected error (%d) in reading CMs from %s", qhstatus, cfg->cmfile);
  }

  /* Terminate outputs... any last words?
   */
  if (tblfp)    cm_tophits_TabularTail(tblfp,    "cmsearch", CM_SEARCH_SEQS, cfg->cmfile, cfg->dbfile, go);
  if (ofp)      fprintf(ofp, "[ok]\n");

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
  free(tinfo);

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

#ifdef HAVE_MPI

#define MAX_BLOCK_SIZE (512*1024)
typedef struct {
  /* A SEQINFO_BLOCK includes sufficient information for a worker to
   * determine exactly where to start and stop reading its sequence
   * block B from the file. We take advantage of SSI index information,
   * which is required in MPI mode. 
   *
   * A worker will start reading sequences at position <first_from> of
   * the sequence with primary key number <first_idx> and finish
   * reading at position <final_to> of the sequence with primary key
   * number <final_idx>, reading complete sequences with primary keys
   * <first_idx>+1..<final_idx>-1. 
   *
   * Note that <first_idx> may equal <final_idx>
   * and <first_from> and <final_to> may be interior sequence
   * coordinates (!= 1, and != L) indicating the full block is a
   * single subsequence of a single database sequence.
   *
   * The workers don't have to worry about reading overlapping chunks
   * of the database, the master has already assured that each block
   * is of an appropriate size and overlaps appropriately with other
   * blocks sent to (possibly) other workers.
   *
   */
  int64_t   first_idx;   /* first SSI primary key number in block */
  int64_t   final_idx;   /* final SSI primary key number in block */
  int64_t   first_from;  /* first sequence position in first sequence in block */
  int64_t   final_to;    /* final sequence position in final sequence in block */
  int64_t   blockL;      /* the number of total residues that exist in this block */
  int       complete;    /* TRUE if this SEQINFO_BLOCK is complete, all data should be valid */
} SEQINFO_BLOCK;

typedef struct {
  SEQINFO_BLOCK  **blocks;    /* array of SEQINFO_BLOCKS */
  int              N;         /* number of SEQINFO_BLOCK elements in blocks */
  int              nalloc;    /* number of blocks elements allocated blocks */
} SEQINFO_BLOCK_LIST;

static SEQINFO_BLOCK      *create_seqinfo_block();
static void                dump_seqinfo_block(FILE *fp, SEQINFO_BLOCK *block);
static SEQINFO_BLOCK_LIST *create_seqinfo_block_list();
static void                dump_seqinfo_block_list(FILE *fp, SEQINFO_BLOCK_LIST *list);
static void                free_seqinfo_block_list(SEQINFO_BLOCK_LIST *list);

static int  seqinfo_block_mpi_send(SEQINFO_BLOCK *block, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
static int  seqinfo_block_mpi_recv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, SEQINFO_BLOCK **ret_block);
static int  seqinfo_block_mpi_pack_size(SEQINFO_BLOCK *block, MPI_Comm comm, int *ret_n);
static int  seqinfo_block_mpi_pack(SEQINFO_BLOCK *block, char *buf, int n, int *pos, MPI_Comm comm);
static int  seqinfo_block_mpi_unpack(char *buf, int n, int *pos, MPI_Comm comm, SEQINFO_BLOCK **ret_block);

static int add_blocks(ESL_SQFILE *dbfp, ESL_SQ *sq, int64_t ncontext, char *errbuf, SEQINFO_BLOCK_LIST *block_list, 
		      int64_t *ret_new_idx, int64_t *ret_noverlap, int64_t *ret_nseq);
static int inspect_next_sequence_using_ssi(ESL_SQFILE *dbfp, ESL_SQ *sq, int64_t ncontext, int64_t pkey_idx, char *errbuf, 
					   SEQINFO_BLOCK_LIST *block_list, int64_t *ret_noverlap);

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

/* mpi_master()
 * The MPI version of cmsearch.
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
 * p7_Fail()'s, which will cause MPI to shut down the worker processes
 * uncleanly.
 */
static int
mpi_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  FILE            *ofp      = stdout;            /* results output file (-o)                        */
  FILE            *afp      = NULL;              /* alignment output file (-A)                      */
  FILE            *tblfp    = NULL;              /* output stream for tabular per-seq (--tblout)    */

  CMFILE          *cmfp;		         /* open input CM file stream                       */
  P7_HMM         **hmmA    = NULL;               /* HMMs, used for filtering                        */
  ESL_SQFILE      *dbfp     = NULL;              /* open input sequence file                        */
  ESL_SQ          *dbsq     = NULL;              /* one target sequence (digital)                   */
  ESL_ALPHABET    *abc      = NULL;              /* digital alphabet                                */
  ESL_STOPWATCH   *w;
  int              textw    = 0;
  int              nquery   = 0;
  int              status   = eslOK;
  int              hstatus  = eslOK;
  int              sstatus  = eslOK;
  int              dest;

  WORKER_INFO     *info          = NULL;         /* contains CM, HMMs, p7 profiles, cmcons, etc. */

  int              i;

  char            *mpi_buf  = NULL;              /* buffer used to pack/unpack structures */
  int              mpi_size = 0;                 /* size of the allocated buffer */
  SEQINFO_BLOCK_LIST  *block_list = NULL;

  int              size;
  MPI_Status       mpistatus;
  char             errbuf[eslERRBUFSIZE];
  int64_t          pkey_idx;
  SEQINFO_BLOCK   *cur_block = NULL;
  int64_t          nseq = 0;                 /* number of sequences in one block */
  int64_t          tot_nseq = 0;             /* number of sequences in all blocks for one CM */
  int64_t          noverlap = 0;             /* number of overlapping residues in one block */
  int64_t          tot_noverlap = 0;         /* number of overlapping residues in all blocks for one CM */

  w = esl_stopwatch_Create();

  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");

  /* Open the database file */
  if((status = open_dbfile    (go, cfg, errbuf, &dbfp)) != eslOK) mpi_failure(errbuf); 
  /* MPI requires SSI indexing */
  if((status = open_dbfile_ssi(go, cfg, dbfp, errbuf)) != eslOK)  mpi_failure(errbuf); 
  /* Determine database size */
  if(esl_opt_IsUsed(go, "-Z")) { /* enabling -Z is the only way to bypass the SSI index file requirement */
    cfg->Z = (int64_t) (esl_opt_GetReal(go, "-Z") * 1000000.); 
  }
  else { 
    if((status = determine_dbsize_using_ssi(go, cfg, dbfp, errbuf)) != eslOK) mpi_failure(errbuf); 
  }
  printf("MPIM cfg->Z: %" PRId64 " residues\n", cfg->Z);

  /* Open the query CM file */
  if ((cmfp = CMFileOpen(cfg->cmfile, NULL)) == NULL)
    esl_fatal("Failed to open covariance model save file %s\n", cfg->cmfile);

  /* Open the results output files */
  if (esl_opt_IsOn(go, "-o") && (ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL)
    mpi_failure("Failed to open output file %s for writing\n",    esl_opt_GetString(go, "-o"));

  if (esl_opt_IsOn(go, "-A") && (afp = fopen(esl_opt_GetString(go, "-A"), "w")) == NULL) 
    mpi_failure("Failed to open alignment file %s for writing\n", esl_opt_GetString(go, "-A"));

  if (esl_opt_IsOn(go, "--tblout") && (tblfp = fopen(esl_opt_GetString(go, "--tblout"), "w")) == NULL)
    mpi_failure("Failed to open tabular per-seq output file %s for writing\n", esl_opt_GetString(go, "--tblfp"));

  /* allocate and initialize <info> which will hold the CMs, HMMs, etc. */
  if((info = create_info()) == NULL) mpi_failure("Out of memory");

  /* <abc> is not known 'til first CM is read. */
  hstatus = CMFileRead(cmfp, errbuf, &abc, &(info->cm));
  if (hstatus == eslOK)
    {
      /* One-time initializations after alphabet <abc> becomes known */
      output_header(ofp, go, cfg->cmfile, cfg->dbfile);
      dbsq = esl_sq_CreateDigital(abc);
    }

  /* Outer loop: over each query CM in <cmfile>. */
  while (hstatus == eslOK) 
    {
      nquery++;
      esl_stopwatch_Start(w);

      /* Make sure we have E-value stats for both the CM and the p7, if not we can't run the pipeline */
      if(! (info->cm->flags & CMH_EXPTAIL_STATS)) mpi_failure("no E-value parameters were read for CM: %s\n", info->cm->name);
      if(! (info->cm->flags & CMH_MLP7_STATS))    mpi_failure("no plan7 HMM E-value parameters were read for CM: %s\n", info->cm->name);

      fprintf(ofp, "Query:       %s  [CLEN=%d]\n", info->cm->name, info->cm->clen);
      if (info->cm->acc)  fprintf(ofp, "Accession:   %s\n", info->cm->acc);
      if (info->cm->desc) fprintf(ofp, "Description: %s\n", info->cm->desc);

      /* Setup the CM and calculate QDBs, if nec */
      if((status = setup_cm(go, info, TRUE, errbuf)) != eslOK) mpi_failure(errbuf);
      if((status = setup_qdbs(go, info, errbuf))     != eslOK) mpi_failure(errbuf);
      /* Setup HMM filters */
      if((status = setup_hmm_filters(go, info, abc, errbuf, &hmmA)) != eslOK) mpi_failure(errbuf);

      /* Create processing pipeline and hit list */
      info->th  = cm_tophits_Create(); 
      info->pli = cm_pipeline_Create(go, info->omA[0]->M, 100, cfg->Z, CM_SEARCH_SEQS); /* L_hint = 100 is just a dummy for now */
      cm_pli_NewModel(info->pli, info->cm, info->fcyk_dmin, info->fcyk_dmax, info->final_dmin, info->final_dmax, info->omA, info->bgA, info->nhmm);

      /* Create block_list and add an empty block to it */
      block_list = create_seqinfo_block_list();
      block_list->blocks[0] = create_seqinfo_block();
      block_list->N = 1;

      /* Main loop: */
      pkey_idx = 0;
      sstatus = eslOK;
      tot_noverlap = noverlap = 0;
      tot_nseq = nseq = 0;
      while((pkey_idx < dbfp->data.ascii.ssi->nprimary) && 
	    ((sstatus = add_blocks(dbfp, dbsq, info->cm->W, errbuf, block_list, &pkey_idx, &noverlap, &nseq)) == eslOK))
	{
	  tot_noverlap += noverlap;
	  tot_nseq     += nseq;

	  cur_block = block_list->blocks[0];
	  while(cur_block != NULL && cur_block->complete) { 
	    /* cur_block points to a valid, complete block, send it to a worker */
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
	    
	    seqinfo_block_mpi_send(cur_block, dest, INFERNAL_BLOCK_TAG, MPI_COMM_WORLD, &mpi_buf, &mpi_size);
	    
	    /* Free the block we just sent, and remove it from block_list 
	     * by swapping pointers and decrementing N. Then update cur_block. */
	    free(cur_block);
	    for(i = 0; i < (block_list->N-1); i++) { 
	      block_list->blocks[i] = block_list->blocks[i+1];
	    }
	    block_list->blocks[(block_list->N-1)] = NULL;
	    block_list->N--;
	    cur_block = block_list->blocks[0];
	  }
	}
      if(sstatus != eslOK) mpi_failure(errbuf);

      /* create an empty block to send to workers */
      cur_block = create_seqinfo_block();
      cur_block->complete = TRUE;

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
	  seqinfo_block_mpi_send(cur_block, dest, INFERNAL_BLOCK_TAG, MPI_COMM_WORLD, &mpi_buf, &mpi_size);
	  //MPI_Send(&block, 3, MPI_LONG_LONG_INT, dest, INFERNAL_BLOCK_TAG, MPI_COMM_WORLD);

	  /* wait for the results */
	  if ((status = cm_tophits_MPIRecv(dest, INFERNAL_TOPHITS_TAG, MPI_COMM_WORLD, &mpi_buf, &mpi_size, &mpi_th)) != eslOK)
	    mpi_failure("Unexpected error %d receiving tophits from %d", status, dest);

	  if ((status = cm_pipeline_MPIRecv(dest, INFERNAL_PIPELINE_TAG, MPI_COMM_WORLD, &mpi_buf, &mpi_size, go, &mpi_pli)) != eslOK)
	    mpi_failure("Unexpected error %d receiving pipeline from %d", status, dest);

	  cm_tophits_Merge(info->th, mpi_th);
	  cm_pipeline_Merge(info->pli, mpi_pli);

	  cm_pipeline_Destroy(mpi_pli, NULL);
	  cm_tophits_Destroy(mpi_th);
	}
      /* Set number of seqs, and subtract number of overlapping residues searched from total */
      info->pli->nseqs = tot_nseq;
      if(info->pli->do_top && info->pli->do_bot) tot_noverlap *= 2; /* count overlaps twice on each strand, if nec */
      printf("nres: %ld\nnoverlap: %ld\n", info->pli->nres, noverlap);
      info->pli->nres -= tot_noverlap;
      printf("nres: %ld\n", info->pli->nres);

      /* Sort and enforce threshold */
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

      cm_tophits_Targets(ofp, info->th, info->pli, textw); 
      fprintf(ofp, "\n\n");
      if(info->pli->show_alignments) {
	if((status = cm_tophits_HitAlignments(ofp, info->th, info->pli, textw)) != eslOK) esl_fatal("Out of memory");
	fprintf(ofp, "\n\n");
      }
      if (tblfp != NULL) { 
	cm_tophits_TabularTargets(tblfp, info->cm->name, info->cm->acc, info->th, info->pli, (nquery == 1)); 
      }
      esl_stopwatch_Stop(w);
      cm_pli_Statistics(ofp, info->pli, w);
      fprintf(ofp, "//\n");

#if 0
	/* WHEN -A IS IMPLEMENTED, UNCOMMENT THIS BLOCK */
      /* Output the results in an MSA (-A option) */
      if (afp) {
	ESL_MSA *msa = NULL;

	if (cm_tophits_Alignment(th, abc, NULL, NULL, 0, p7_ALL_CONSENSUS_COLS, &msa) == eslOK)
	  {
	    if (textw > 0) esl_msa_Write(afp, msa, eslMSAFILE_STOCKHOLM);
	    else           esl_msa_Write(afp, msa, eslMSAFILE_PFAM);
	  
	    fprintf(ofp, "# Alignment of %d hits satisfying inclusion thresholds saved to: %s\n", msa->nseq, esl_opt_GetString(go, "-A"));
	  } 
	else fprintf(ofp, "# No hits satisfy inclusion thresholds; no alignment saved\n");
	  
	esl_msa_Destroy(msa);
      }
#endif
      free_info(info);
      free(info);
      free_seqinfo_block_list(block_list);
      block_list = NULL;

      hstatus = CMFileRead(cmfp, errbuf, &abc, &(info->cm));
      if(hstatus == eslOK) { 
	if((info = create_info()) == NULL) mpi_failure("Out of memory"); /* for the next model */
      }
    } /* end outer loop over query HMMs */
  
  switch(hstatus) {
  case eslEFORMAT:   mpi_failure("bad file format in CM file %s",             cfg->cmfile);      break;
  case eslEOF:       /* EOF is good, that's what we expect here */                                break;
  default:           mpi_failure("Unexpected error (%d) in reading HMMs from %s", hstatus, cfg->cmfile);
  }
  
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

  /* Terminate outputs... any last words?
   */
  if (tblfp)    cm_tophits_TabularTail(tblfp,    "cmsearch", CM_SEARCH_SEQS, cfg->cmfile, cfg->dbfile, go);
  if (ofp)      fprintf(ofp, "[ok]\n");

  /* Cleanup - prepare for exit
   */
  free(block_list);
  if (mpi_buf != NULL) free(mpi_buf);

  CMFileClose(cmfp);
  esl_sqfile_Close(dbfp);
  esl_alphabet_Destroy(abc);
  esl_sq_Destroy(dbsq);
  esl_stopwatch_Destroy(w);

  if (ofp != stdout) fclose(ofp);
  if (afp)           fclose(afp);
  if (tblfp)         fclose(tblfp);

  return eslOK;

 ERROR:
  return eslEMEM;
}


static int
mpi_worker(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  ESL_SQ          *dbsq     = NULL;              /* one target sequence (digital)                   */
  ESL_ALPHABET    *abc      = NULL;              /* digital alphabet                                */
  CMFILE          *cmfp;		         /* open input CM file stream                       */
  P7_HMM         **hmmA    = NULL;               /* HMMs, used for filtering                        */
  ESL_SQFILE      *dbfp     = NULL;              /* open input sequence file                        */
  ESL_STOPWATCH   *w;
  int              status   = eslOK;
  int              hstatus  = eslOK;
  int              prev_hit_cnt;
  double           eff_dbsize;                   /* the effective database size, max # hits expected in db of this size */

  WORKER_INFO     *info          = NULL;         /* contains CM, HMMs, p7 profiles, cmcons, etc. */

  char            *mpi_buf  = NULL;              /* buffer used to pack/unpack structures           */
  int              mpi_size = 0;                 /* size of the allocated buffer                    */

  char             errbuf[eslERRBUFSIZE];
  int64_t          pkey_idx;
  int64_t          L;
  int64_t          seq_from, seq_to;
  int64_t          readL;
  char            *pkey;

  w = esl_stopwatch_Create();

  /* Open the database file */
  if((status = open_dbfile    (go, cfg, errbuf, &dbfp)) != eslOK) mpi_failure(errbuf); 
  /* MPI requires SSI indexing */
  if((status = open_dbfile_ssi(go, cfg, dbfp, errbuf))  != eslOK) mpi_failure(errbuf); 
  /* Determine database size */
  if(esl_opt_IsUsed(go, "-Z")) { /* enabling -Z is the only way to bypass the SSI index file requirement */
    cfg->Z = (int64_t) (esl_opt_GetReal(go, "-Z") * 1000000.); 
  }
  else { 
    if((status = determine_dbsize_using_ssi(go, cfg, dbfp, errbuf)) != eslOK) mpi_failure(errbuf); 
  }
  printf("MPIW cfg->Z: %" PRId64 " residues\n", cfg->Z);

  /* Open the query CM file */
  if ((cmfp = CMFileOpen(cfg->cmfile, NULL)) == NULL)
    esl_fatal("Failed to open covariance model save file %s\n", cfg->cmfile);

  /* allocate and initialize <info> which will hold the CMs, HMMs, etc. */
  if((info = create_info()) == NULL) mpi_failure("Out of memory");

  /* <abc> is not known 'til first CM is read. */
  hstatus = CMFileRead(cmfp, errbuf, &abc, &(info->cm));
  if (hstatus == eslOK)
    {
      /* One-time initializations after alphabet <abc> becomes known */
      dbsq = esl_sq_CreateDigital(abc);
    }

  /* Outer loop: over each query CM in <cmfile>. */
  while (hstatus == eslOK) 
    {
      SEQINFO_BLOCK *block;

      esl_stopwatch_Start(w);

      status = 0;
      MPI_Send(&status, 1, MPI_INT, 0, INFERNAL_READY_TAG, MPI_COMM_WORLD);

      /* Setup the CM and calculate QDBs, if nec */
      if((status = setup_cm  (go, info, TRUE, errbuf)) != eslOK) mpi_failure(errbuf);
      if((status = setup_qdbs(go, info, errbuf))       != eslOK) mpi_failure(errbuf);
      /* Setup HMM filters */
      if((status = setup_hmm_filters(go, info, abc, errbuf, &hmmA)) != eslOK) mpi_failure(errbuf);

      /* Create processing pipeline and hit list */
      info->th  = cm_tophits_Create(); 
      info->pli = cm_pipeline_Create(go, info->omA[0]->M, 100, cfg->Z, CM_SEARCH_SEQS); /* L_hint = 100 is just a dummy for now */
      cm_pli_NewModel(info->pli, info->cm, info->fcyk_dmin, info->fcyk_dmax, info->final_dmin, info->final_dmax, info->omA, info->bgA, info->nhmm);

      /* receive a sequence info block from the master */
      if((status = seqinfo_block_mpi_recv(0, INFERNAL_BLOCK_TAG, MPI_COMM_WORLD, &mpi_buf, &mpi_size, &block)) != eslOK) mpi_failure("Failed to receive sequence block, error status code: %d\n", status); 

      while(block->first_idx != -1) { /* receipt of a block with first_idx == -1 signals us that we're done with the database */
	readL = 0;
	pkey_idx = block->first_idx;
	while(pkey_idx <= block->final_idx) 
	  { 
	    /* determine the primary key and length of the sequence */
	    status = esl_ssi_FindNumber(dbfp->data.ascii.ssi, pkey_idx, NULL, NULL, NULL, &L, &pkey);
	    
	    /* determine if we're reading the full sequence or a subsequence */
	    if((pkey_idx != block->first_idx) && (pkey_idx != block->final_idx)) { 
	      /* read full sequence */
	      if((status = esl_sqio_Fetch(dbfp, pkey, dbsq)) != eslOK) mpi_failure("unable to fetch sequence %s\n", pkey);
	      seq_from = 1;
	      seq_to   = L;
	    }
	    else {
	      /* read a subsequence */
	      seq_from = (pkey_idx == block->first_idx) ? block->first_from : 1;
	      seq_to   = (pkey_idx == block->final_idx) ? block->final_to   : L;
	      if((status = esl_sqio_FetchSubseq(dbfp, pkey, seq_from, seq_to, dbsq)) != eslOK) mpi_failure("unable to fetch subsequence %ld..%ld of %s\n", seq_from, seq_to, pkey);
	      /* FetchSubseq renames the sequence, set it back */
	      esl_sq_SetName(dbsq, dbsq->source);
	    }
	    /* tell pipeline we've got a new sequence (this updates info->pli->nres) */
	    cm_pli_NewSeq(info->pli, info->cm, dbsq);
	    readL += dbsq->n;

	    /* search top strand */
	    if (info->pli->do_top) { 
	      prev_hit_cnt = info->th->N;
	      if((status = cm_Pipeline(info->pli, info->cm, info->cmcons, info->omA, info->gmA, info->bgA, info->p7_evparamAA, info->nhmm, dbsq, info->th)) != eslOK) 
		mpi_failure("cm_pipeline() failed unexpected with status code %d\n", status);
	      cm_pipeline_Reuse(info->pli); /* prepare for next search */

	      /* If we're a subsequence, update hit positions so they're relative 
	       * to the full-length sequence. */
	      if(seq_from != 1) update_hit_positions(info->th, prev_hit_cnt, seq_from, FALSE);
	    }
	    
	    /* search reverse complement */
	    if(info->pli->do_bot && dbsq->abc->complement != NULL) { 
	      esl_sq_ReverseComplement(dbsq);
	      prev_hit_cnt = info->th->N;
	      if((status = cm_Pipeline(info->pli, info->cm, info->cmcons, info->omA, info->gmA, info->bgA, info->p7_evparamAA, info->nhmm, dbsq, info->th)) != eslOK) 
		mpi_failure("cm_pipeline() failed unexpected with status code %d\n", status);
	      cm_pipeline_Reuse(info->pli); /* prepare for next search */
	      info->pli->nres += dbsq->n; /* add dbsq->n residues, the reverse complement we just searched */

	      /* Hit positions will be relative to the reverse-complemented sequence
	       * (i.e. start < end), which may be a subsequence. Update hit positions 
	       * so they're relative to the full-length original sequence (start > end). */
	      update_hit_positions(info->th, prev_hit_cnt, seq_from, TRUE);
	    }
	    esl_sq_Reuse(dbsq);
	    pkey_idx++;
	  }
	/* sanity check to make sure the blocks are the same */
	if (block->blockL != readL) mpi_failure("Block length mismatch - expected %ld found %ld for block idx:%ld from:%ld\n", block->blockL, readL, block->first_idx, block->first_from);

	/* inform the master we need another block of sequences */
	status = 0;
	MPI_Send(&status, 1, MPI_INT, 0, INFERNAL_READY_TAG, MPI_COMM_WORLD);
	
	/* wait for the next block of sequences */
	if((status = seqinfo_block_mpi_recv(0, INFERNAL_BLOCK_TAG, MPI_COMM_WORLD, &mpi_buf, &mpi_size, &block)) != eslOK) mpi_failure("Failed to receive sequence block, error status code: %d\n", status); 
      }   
 
      esl_stopwatch_Stop(w);
      
      /* compute e-values before sending back to master */
      eff_dbsize = (double) ((((double) cfg->Z / (double) info->cm->stats->expAA[info->pli->final_cm_exp_mode][0]->dbsize) * 
			      ((double) info->cm->stats->expAA[info->pli->final_cm_exp_mode][0]->nrandhits)) + 0.5);
      cm_tophits_ComputeEvalues(info->th, eff_dbsize);
      
      cm_tophits_MPISend(info->th,   0, INFERNAL_TOPHITS_TAG,  MPI_COMM_WORLD,  &mpi_buf, &mpi_size);
      cm_pipeline_MPISend(info->pli, 0, INFERNAL_PIPELINE_TAG, MPI_COMM_WORLD,  &mpi_buf, &mpi_size);
      
      free_info(info);
      free(info);
      
      hstatus = CMFileRead(cmfp, errbuf, &abc, &(info->cm));
      if(hstatus == eslOK) { 
	if((info = create_info()) == NULL) mpi_failure("Out of memory"); /* info is for the next model */
      }
    } /* end outer loop over query CMs */
  
  switch(hstatus) {
  case eslEFORMAT:   mpi_failure("bad file format in CM file %s",             cfg->cmfile);      break;
  case eslEOF:       /* EOF is good, that's what we expect here */                                break;
  default:           mpi_failure("Unexpected error (%d) in reading HMMs from %s", hstatus, cfg->cmfile);
  }

  status = 0;
  MPI_Send(&status, 1, MPI_INT, 0, INFERNAL_TERMINATING_TAG, MPI_COMM_WORLD);

  if (mpi_buf != NULL) free(mpi_buf);

  CMFileClose(cmfp);
  esl_sqfile_Close(dbfp);

  esl_sq_Destroy(dbsq);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);

  return eslOK;
}
#endif /*HAVE_MPI*/

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
      if((status = cm_Pipeline(info->pli, info->cm, info->cmcons, info->omA, info->gmA, info->bgA, info->p7_evparamAA, info->nhmm, dbsq, info->th)) != eslOK) cm_Fail("cm_pipeline() failed unexpected with status code %d\n", status);
      cm_pipeline_Reuse(info->pli); /* prepare for next search */
      
      /* modify hit positions to account for the position of the window in the full sequence */
      update_hit_positions(info->th, prev_hit_cnt, dbsq->start, FALSE);
    }

    /* reverse complement */
    if (info->pli->do_bot && dbsq->abc->complement != NULL) { 
      prev_hit_cnt = info->th->N;
      esl_sq_ReverseComplement(dbsq);
      if((status = cm_Pipeline(info->pli, info->cm, info->cmcons, info->omA, info->gmA, info->bgA, info->p7_evparamAA, info->nhmm, dbsq, info->th)) != eslOK) cm_Fail("cm_pipeline() failed unexpected with status code %d\n", status);
      cm_pipeline_Reuse(info->pli); /* prepare for next search */
      
      /* modify hit positions to account for the position of the window in the full sequence */
      update_hit_positions(info->th, prev_hit_cnt, dbsq->start, TRUE);

      info->pli->nres += dbsq->W; /* add dbsq->W residues, the reverse complement we just searched */
    }
    
    wstatus = esl_sqio_ReadWindow(dbfp, info->pli->maxW, CMSEARCH_MAX_RESIDUE_COUNT, dbsq);
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
	((ESL_SQ_BLOCK *)newBlock)->list->C = info->pli->maxW;
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
	if((status = cm_Pipeline(info->pli, info->cm, info->cmcons, info->omA, info->gmA, info->bgA, info->p7_evparamAA, info->nhmm, dbsq, info->th)) != eslOK) cm_Fail("cm_pipeline() failed unexpected with status code %d\n", status);
	cm_pipeline_Reuse(info->pli); /* prepare for next search */
	
	/* modify hit positions to account for the position of the window in the full sequence */
	update_hit_positions(info->th, prev_hit_cnt, dbsq->start, FALSE);
      }
	
      /* reverse complement */
      if (info->pli->do_bot && dbsq->abc->complement != NULL) {
	prev_hit_cnt = info->th->N;
	esl_sq_ReverseComplement(dbsq);
	if((status = cm_Pipeline(info->pli, info->cm, info->cmcons, info->omA, info->gmA, info->bgA, info->p7_evparamAA, info->nhmm, dbsq, info->th)) != eslOK) cm_Fail("cm_pipeline() failed unexpected with status code %d\n", status);
	cm_pipeline_Reuse(info->pli); /* prepare for next search */

	/* modify hit positions to account for the position of the window in the full sequence */
	update_hit_positions(info->th, prev_hit_cnt, dbsq->start, TRUE);

	info->pli->nres += dbsq->W; /* add dbsq->W residues, the reverse complement we just searched */
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

/* Function:  open_dbfile
 * Synopsis:  Open the database file.
 * Incept:    EPN, Mon Jun  6 09:13:08 2011
 *
 * Returns: eslOK on success. ESL_SQFILE *dbfp in *ret_dbfp.
 *          Upon error, fills errbuf with error message and
 *          returns error status code.
 */
static int
open_dbfile(ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, ESL_SQFILE **ret_dbfp)
{
  int status;
  int dbfmt    = eslSQFILE_UNKNOWN; /* format of dbfile                                */

  if (esl_opt_IsOn(go, "--tformat")) {
    dbfmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--tformat"));
    if (dbfmt == eslSQFILE_UNKNOWN) ESL_FAIL(status, errbuf, "%s is not a recognized sequence database file format\n", esl_opt_GetString(go, "--tformat"));
  }
  status = esl_sqfile_Open(cfg->dbfile, dbfmt, p7_SEQDBENV, ret_dbfp);
  if      (status == eslENOTFOUND) ESL_FAIL(status, errbuf, "Failed to open sequence file %s for reading\n",          cfg->dbfile);
  else if (status == eslEFORMAT)   ESL_FAIL(status, errbuf, "Sequence file %s is empty or misformatted\n",            cfg->dbfile);
  else if (status == eslEINVAL)    ESL_FAIL(status, errbuf, "Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        ESL_FAIL(status, errbuf, "Unexpected error %d opening sequence file %s\n", status, cfg->dbfile);  

  return eslOK;
}

/* Function:  open_dbfile_ssi
 * Synopsis:  Open database file's SSI index.
 * Incept:    EPN, Wed Jun  8 07:19:49 2011
 *
 * Returns: eslOK on succes.
 *          Upon error, error status is returned, and errbuf is filled.
 */
static int
open_dbfile_ssi(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_SQFILE *dbfp, char *errbuf)
{
  int status;

  /* Open the database file's SSI index for retrieval.  Some error
   * messages differ depending on if we're in MPI mode or not, because
   * SSI indexing is required in MPI mode, but not required in serial
   * mode if -Z is also used.
   */
#ifdef HAVE_MPI
  if (cfg->do_mpi) 
    { 
      if (dbfp->data.ascii.do_gzip)   ESL_FAIL(status, errbuf, "Reading gzipped sequence files is not supported with --mpi.");
      if (dbfp->data.ascii.do_stdin)  ESL_FAIL(status, errbuf, "Reading sequence files from stdin is not supported with --mpi.");
    }
  else
#endif
    {
      if (dbfp->data.ascii.do_gzip)   ESL_FAIL(status, errbuf, "Reading gzipeed sequence files is not supported unless -Z is used.");
      if (dbfp->data.ascii.do_stdin)  ESL_FAIL(status, errbuf, "Reading sequence files from stdin is not supported unless -Z is used.");
    }
  status = esl_sqfile_OpenSSI(dbfp, NULL);
  if      (status == eslEFORMAT) ESL_FAIL(status, errbuf, "SSI index for database file is in incorrect format\n");
  else if (status == eslERANGE)  ESL_FAIL(status, errbuf, "SSI index for database file is in 64-bit format and we can't read it\n");
  else if (status != eslOK) { 
#ifdef HAVE_MPI
    if (cfg->do_mpi) 
      { 
	ESL_FAIL(status, errbuf, "Failed to open SSI index, use esl-sfetch to index the database file.\n");
      }
    else
#endif
      {
	ESL_FAIL(status, errbuf, "Failed to open SSI index, use -Z or esl-sfetch to index the database file.\n");
      }
  }
  return status;
}

/* Function:  determine_dbsize_using_ssi
 * Synopsis:  Determine size of the database using SSI index.
 * Incept:    EPN, Mon Jun  6 09:17:32 2011
 *
 * Returns:   eslOK on success. cfg->Z is set.
 *            Upon error, error status is returned and errbuf is filled.
 *          
 */
static int
determine_dbsize_using_ssi(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_SQFILE *dbfp, char *errbuf)
{
  int status;
  int64_t L;    
  int i;
  
  if(dbfp->data.ascii.ssi == NULL) ESL_FAIL(status, errbuf, "SSI index failed to open");

  /* step through sequence SSI file to get database size */
  cfg->Z = 0;
  for(i = 0; i < dbfp->data.ascii.ssi->nprimary; i++) { 
    status = esl_ssi_FindNumber(dbfp->data.ascii.ssi, i, NULL, NULL, NULL, &L, NULL);
    if(status == eslENOTFOUND) ESL_FAIL(status, errbuf, "unable to find sequence %d in SSI index file, try re-indexing with esl-sfetch.", i);
    if(status == eslEFORMAT)   ESL_FAIL(status, errbuf, "SSI index for database file is in incorrect format.");
    if(status != eslOK)        ESL_FAIL(status, errbuf, "proble with SSI index for database file.");
    cfg->Z += L;
  }

  return eslOK;
}

/* Function:  create_info
 * Incept:    EPN, Mon Jun  6 14:52:29 2011
 *
 * Purpose:  Create (allocate and initalize) a WORKER_INFO object.
 *
 * Returns:  The new WORKER_INFO is returned. NULL is returned upon an error.
 */
WORKER_INFO *
create_info()
{ 
  int status;
  WORKER_INFO *info = NULL;

  ESL_ALLOC(info, sizeof(WORKER_INFO));
  info->pli          = NULL;
  info->th           = NULL;
  info->cm           = NULL;
  info->omA          = NULL;
  info->bgA          = NULL;
  info->p7_evparamAA = NULL;
  info->nhmm         = 0;
  info->use_mlp7_hmm = FALSE;
  info->fcyk_dmin    = NULL;
  info->fcyk_dmax    = NULL;
  info->final_dmin   = NULL;
  info->final_dmax   = NULL;

  return info;

 ERROR: 
  if(info != NULL) free(info);
  return NULL;
}

/* Function:  clone_info
 * Incept:    EPN, Mon Jun  6 10:24:21 2011
 *
 * Purpose: Given a template WORKER_INFO <tinfo>, clone it into the
 *          <dest_infocnt> WORKER_INFOs in <dest_infoA>. After cloning the CM,
 *          configure it.
 *
 * Returns: <eslOK> on success. Dies immediately upon an error.
 */
int
clone_info(ESL_GETOPTS *go, WORKER_INFO *src_info, WORKER_INFO *dest_infoA, int dest_infocnt, char *errbuf)
{ 
  int status;
  int i, m;

  for (i = 0; i < dest_infocnt; ++i) {
    /* Clone the CM in src_info, it should have just been read from the file (it should not have been
     * configured yet). It's difficult to clone a configured CM due to extra data structures that 
     * get added during configuration. 
     */
    if((status = CloneCMJustReadFromFile(src_info->cm, errbuf, &(dest_infoA[i].cm))) != eslOK) return status;
    dest_infoA[i].th           = cm_tophits_Create();
    dest_infoA[i].nhmm         = src_info->nhmm;
    dest_infoA[i].use_mlp7_hmm = src_info->use_mlp7_hmm;
    ESL_ALLOC(dest_infoA[i].gmA, sizeof(P7_PROFILE *)  * src_info->nhmm);
    ESL_ALLOC(dest_infoA[i].omA, sizeof(P7_OPROFILE *) * src_info->nhmm);
    ESL_ALLOC(dest_infoA[i].bgA, sizeof(P7_BG *)       * src_info->nhmm);
    ESL_ALLOC(dest_infoA[i].p7_evparamAA, sizeof(float *)  * src_info->nhmm);

    /* configure the CM */
    if((status = setup_cm(go, (&dest_infoA[i]), 
			  (i == 0) ? TRUE : FALSE, /* calculate W? */
			  errbuf)) != eslOK) return status;
    if(i > 0) dest_infoA[i].cm->W = dest_infoA[0].cm->W;

    for(m = 0; m < src_info->nhmm; m++) { 
      dest_infoA[i].gmA[m]  = p7_profile_Clone(src_info->gmA[m]);
      dest_infoA[i].omA[m]  = p7_oprofile_Clone(src_info->omA[m]);
      dest_infoA[i].bgA[m]  = p7_bg_Create(src_info->bgA[m]->abc);
      ESL_ALLOC(dest_infoA[i].p7_evparamAA[m], sizeof(float) * CM_p7_NEVPARAM);
      if(m == 0 && dest_infoA[i].use_mlp7_hmm) { esl_vec_FCopy(src_info->cm->mlp7_evparam,     CM_p7_NEVPARAM, dest_infoA[i].p7_evparamAA[m]); }
      else                                     { esl_vec_FCopy(src_info->cm->ap7_evparamAA[m], CM_p7_NEVPARAM, dest_infoA[i].p7_evparamAA[m]); }
    }

    if(src_info->fcyk_dmin != NULL) { 
      ESL_ALLOC(dest_infoA[i].fcyk_dmin, sizeof(int) * src_info->cm->M);
      esl_vec_ICopy(src_info->fcyk_dmin, src_info->cm->M, dest_infoA[i].fcyk_dmin);
    }
    if(src_info->fcyk_dmax != NULL) { 
      ESL_ALLOC(dest_infoA[i].fcyk_dmax, sizeof(int) * src_info->cm->M);
      esl_vec_ICopy(src_info->fcyk_dmax, src_info->cm->M, dest_infoA[i].fcyk_dmax);
    }
    if(src_info->final_dmin != NULL) { 
      ESL_ALLOC(dest_infoA[i].final_dmin, sizeof(int) * src_info->cm->M);
      esl_vec_ICopy(src_info->final_dmin, src_info->cm->M, dest_infoA[i].final_dmin);
    }
    if(src_info->final_dmax != NULL) { 
      ESL_ALLOC(dest_infoA[i].final_dmax, sizeof(int) * src_info->cm->M);
      esl_vec_ICopy(src_info->final_dmax, src_info->cm->M, dest_infoA[i].final_dmax);
    }
  }
  return eslOK;
  
  ERROR:
  ESL_FAIL(status, errbuf, "Out of memory");
}

/* Function:  free_info
 * Incept:    EPN, Mon Jun  6 10:46:27 2011
 *
 * Purpose:  Free a WORKER_INFO object.
 *
 * Returns: void. Dies immediately upon an error.
 */
void
free_info(WORKER_INFO *info)
{ 
  int m;
  
  if(info->pli != NULL) cm_pipeline_Destroy(info->pli, info->cm);
  if(info->th  != NULL) cm_tophits_Destroy(info->th);
  
  if(info->cm     != NULL) FreeCM(info->cm);
  if(info->cmcons != NULL) FreeCMConsensus(info->cmcons);
  
  for(m = 0; m < info->nhmm; m++) { 
    p7_oprofile_Destroy(info->omA[m]);
    p7_profile_Destroy(info->gmA[m]);
    p7_bg_Destroy(info->bgA[m]);
    free(info->p7_evparamAA[m]);
  }
  free(info->omA);
  free(info->gmA);
  free(info->bgA);
  free(info->p7_evparamAA);

  if(info->fcyk_dmin  != NULL) free(info->fcyk_dmin);
  if(info->fcyk_dmax  != NULL) free(info->fcyk_dmax);
  if(info->final_dmin != NULL) free(info->final_dmin);
  if(info->final_dmin != NULL) free(info->final_dmax);
  
  return;
}

/* Function:  setup_cm()
 * Incept:    EPN, Mon Jun  6 12:06:38 2011
 *
 * Purpose:  Given a WORKER_INFO <info> with a valid CM just read
 *           from a file, configure the CM and create the CMConsensus data. 
 *
 * Returns: eslOK on success. Upon an error, fills errbuf with
 *          message and returns appropriate error status code.
 */
int
setup_cm(ESL_GETOPTS *go, WORKER_INFO *info, int do_calculate_W, char *errbuf)
{ 
  int status;

  /* Configure the CM */
  if(! esl_opt_GetBoolean(go, "-g")) { 
    if(! esl_opt_GetBoolean(go, "--cp9gloc")) { 
      info->cm->config_opts |= CM_CONFIG_LOCAL;
      info->cm->config_opts |= CM_CONFIG_HMMLOCAL;
      if(! esl_opt_GetBoolean(go, "--cp9noel")) info->cm->config_opts |= CM_CONFIG_HMMEL; 
    }
  }
  if((status = ConfigCM(info->cm, errbuf, do_calculate_W, NULL, NULL)) != eslOK) { 
    return status;
  }
  if((status = CreateCMConsensus(info->cm, info->cm->abc, 3.0, 1.0, &(info->cmcons)))!= eslOK) {
    ESL_FAIL(status, errbuf, "Failed to create CMConsensus data structure.\n");
  }      

  return eslOK;
}

/* Function:  setup_qdbs()
 * Incept:    EPN, Mon Jun  6 14:04:36 2011
 *
 * Purpose:  Given a WORKER_INFO <info> with a valid CM 
 *           determine QDBs, if necessary.
 *
 * Returns: eslOK on success. Upon an error, fills errbuf with
 *          message and returns appropriate error status code.
 */
int
setup_qdbs(ESL_GETOPTS *go, WORKER_INFO *info, char *errbuf)
{
  int safe_W;

  /* Determine QDBs, if necessary. By default they are not needed. */
  info->fcyk_dmin  = info->fcyk_dmax  = NULL; 
  info->final_dmin = info->final_dmax = NULL; 

  if(esl_opt_GetBoolean(go, "--fqdb")) { /* determine CYK filter QDBs */
    safe_W = info->cm->clen * 3;
    while(!(BandCalculationEngine(info->cm, safe_W, esl_opt_GetReal(go, "--fbeta"), FALSE, &(info->fcyk_dmin), &(info->fcyk_dmax), NULL, NULL))) { 
      free(info->fcyk_dmin);
      free(info->fcyk_dmax);
      safe_W *= 2;
      if(safe_W > (info->cm->clen * 1000)) ESL_FAIL(eslEINVAL, errbuf, "Unable to calculate QDBs");
    }
  }

  if(esl_opt_GetBoolean(go, "--qdb")) { /* determine CYK filter QDBs */
    safe_W = info->cm->clen * 3;
    while(!(BandCalculationEngine(info->cm, safe_W, esl_opt_GetReal(go, "--beta"), FALSE, &(info->final_dmin), &(info->final_dmax), NULL, NULL))) { 
      free(info->final_dmin);
      free(info->final_dmax);
      safe_W *= 2;
      if(safe_W > (info->cm->clen * 1000)) ESL_FAIL(eslEINVAL, errbuf, "Unable to calculate QDBs");
    }
  }

  return eslOK;
}


/* Function:  setup_hmm_filters()
 * Incept:    EPN, Mon Jun  6 11:31:42 2011
 *
 * Purpose:  Given a WORKER_INFO <info> with a valid CM just read
 *           from a file, set up the HMM filter related data in
 *           <info>.
 *
 * Returns: <eslOK> on success. Dies immediately upon an error.
 */
int
setup_hmm_filters(ESL_GETOPTS *go, WORKER_INFO *info, const ESL_ALPHABET *abc, char *errbuf, P7_HMM ***ret_hmmA)
{
  int      status;
  P7_HMM **hmmA = NULL;
  int      m, z, z_offset;

  /* Convert HMMs to optimized models.
   * First, determine how many and which HMMs we'll be using *
   * by default, we filter only with the <nhmm> (== cm->nap7) HMMs written in the CM file */
  info->use_mlp7_hmm = FALSE;
  info->nhmm = info->cm->nap7;
  /* but, if --doml selected or there's no HMMs in the file, use a ML P7 HMM built from the CM */
  if (esl_opt_GetBoolean(go, "--doml") || (info->cm->nap7 == 0)) { 
    info->use_mlp7_hmm = TRUE;
    info->nhmm = 1 + info->cm->nap7;
  }
  else if(esl_opt_GetBoolean(go, "--noadd")) { 
    /* or if --noadd is selected, only use a ML p7 HMM */
    info->use_mlp7_hmm = TRUE;
    info->nhmm = 1;
  }      

  /* Now we know how many HMMs, allocate and fill the necessary data structures for each */
  ESL_ALLOC(hmmA, sizeof(P7_HMM *) * info->nhmm);
  /* use the ML p7 HMM if nec */
  if (esl_opt_GetBoolean(go, "--doml") || (info->cm->nap7 == 0)) { 
    hmmA[0] = info->cm->mlp7;
  }
  /* copy the HMMs from the file into the CM data structure */
  if(! esl_opt_GetBoolean(go, "--noadd")) { 
    z_offset = info->use_mlp7_hmm ? 1 : 0;
    for(z = 0; z < info->cm->nap7; z++) { 
      hmmA[z+z_offset] = info->cm->ap7A[z];
    }
  }
  
  /* allocate space for profiles, bg and EV params for each HMM in info */
  ESL_ALLOC(info->gmA, sizeof(P7_PROFILE *)     * info->nhmm);
  ESL_ALLOC(info->omA, sizeof(P7_OPROFILE *)    * info->nhmm);
  ESL_ALLOC(info->bgA, sizeof(P7_BG *)          * info->nhmm);
  ESL_ALLOC(info->p7_evparamAA, sizeof(float *) * info->nhmm);
  for(m = 0; m < info->nhmm; m++) { 
    info->gmA[m] = p7_profile_Create (hmmA[m]->M, abc);
    info->omA[m] = p7_oprofile_Create(hmmA[m]->M, abc);
    info->bgA[m] = p7_bg_Create(abc);
    p7_ProfileConfig(hmmA[m], info->bgA[m], info->gmA[m], 100, p7_LOCAL);  /* 100 is a dummy length for now; and MSVFilter requires local mode */
    p7_oprofile_Convert(info->gmA[m], info->omA[m]);                       /* <om> is now p7_LOCAL, multihit */
    /* after omA[m]'s been created, convert gmA[m] to glocal, to define envelopes in cm_pipeline() */
    p7_ProfileConfig(hmmA[m], info->bgA[m], info->gmA[m], 100, p7_GLOCAL);
    /* copy evalue parameters */
    ESL_ALLOC(info->p7_evparamAA[m], sizeof(float) * CM_p7_NEVPARAM);
    if(m == 0 && info->use_mlp7_hmm) { esl_vec_FCopy(info->cm->mlp7_evparam,     CM_p7_NEVPARAM, info->p7_evparamAA[m]); }
    else                             { esl_vec_FCopy(info->cm->ap7_evparamAA[m], CM_p7_NEVPARAM, info->p7_evparamAA[m]); }
  }
  *ret_hmmA = hmmA;
  
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "Out of memory.");
}


/* Function:  update_hit_positions
 * Synopsis:  Update sequence positions in hit lits.
 * Incept:    EPN, Wed May 25 09:12:52 2011
 *
 * Purpose: For hits <hit_start..th->N-1> in a hit list, update the
 *          positions of start, stop and ad->sqfrom, ad->sqto.  This
 *          is necessary when we've searched a chunk of subsequence
 *          that originated somewhere within a larger sequence.
 *
 * Returns: <eslOK> on success.
 */
 int
update_hit_positions(CM_TOPHITS *th, int hit_start, int64_t seq_start, int in_revcomp)
{ 
  if((seq_start == 1) && (in_revcomp == FALSE)) return eslOK;

  int i;
  if(in_revcomp) { 
    for (i = hit_start; i < th->N ; i++) {
      th->unsrt[i].start = seq_start - th->unsrt[i].start + 1;
      th->unsrt[i].stop  = seq_start - th->unsrt[i].stop  + 1;
      if(th->unsrt[i].ad != NULL) { 
	th->unsrt[i].ad->sqfrom = seq_start - th->unsrt[i].ad->sqfrom + 1;
	th->unsrt[i].ad->sqto   = seq_start - th->unsrt[i].ad->sqto + 1;
	th->unsrt[i].ad->L      = seq_start - th->unsrt[i].ad->L + 1;
      }
    }

  }
  else { 
    for (i = hit_start; i < th->N ; i++) {
      th->unsrt[i].start += seq_start-1;
      th->unsrt[i].stop  += seq_start-1;
      if(th->unsrt[i].ad != NULL) { 
	th->unsrt[i].ad->sqfrom += seq_start-1;
	th->unsrt[i].ad->sqto   += seq_start-1;
	th->unsrt[i].ad->L      += seq_start-1;
      }
    }
  }
  return eslOK;
}

/* Function:  determine_nres_to_read()
 * Synopsis:  Determine how many residues to read for esl_sqio_ReadWindow()
 * Incept:    EPN, Tue Jun  7 12:52:28 2011
 *
 * Returns: Number of residues to read in <ret_nres_to_read>
 */
void
determine_nres_to_read(uint64_t readL, uint64_t ncontext, int64_t *ret_nres_to_read)
{ 
  int64_t nres_to_read;
  nres_to_read = ESL_MAX((5. * ncontext), (CMSEARCH_MAX_RESIDUE_COUNT - readL));

  //printf("in determine_nres_to_read:\n");
  //printf("\tncontext: %ld\n", ncontext);
  //printf("\treadL:    %ld\n", readL);
  //printf("\tMAX:      %d\n",  CMSEARCH_MAX_RESIDUE_COUNT);
  //printf("\treturn:   %ld\n", nres_to_read);

  *ret_nres_to_read = nres_to_read;
  return;
}
 
#ifdef HAVE_MPI

/* Function:  seqinfo_block_mpi_send()
 * Synopsis:  Send a SEQINFO_BLOCK in an MPI buffer.
 * Incept:    EPN, Mon Jun  6 19:49:44 2011
 */
int
seqinfo_block_mpi_send(SEQINFO_BLOCK *block, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc)
{
  int   status;
  int   pos;
  int   n;

  /* calculate the buffer size needed to hold the SEQINFO_BLOCK */
  if ((status = seqinfo_block_mpi_pack_size(block, comm, &n)) != eslOK) goto ERROR;

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n);
    *nalloc = n; 
  }

  /* Pack the SEQINFO_BLOCK into the buffer */
  pos  = 0;
  if ((status = seqinfo_block_mpi_pack(block, *buf, n, &pos, comm)) != eslOK) goto ERROR;

  /* Send the packed SEQINFO_BLOCK to the destination. */
  if (MPI_Send(*buf, n, MPI_PACKED, dest, tag, comm) != 0)  ESL_EXCEPTION(eslESYS, "mpi send failed");

  return eslOK;

 ERROR:
  return status;

}

/* Function:  seqinfo_block_mpi_recv()
 * Synopsis:  Receive a SEQINFO_BLOCK in an MPI buffer.
 * Incept:    EPN, Mon Jun  6 19:47:20 2011
 */
int
seqinfo_block_mpi_recv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, SEQINFO_BLOCK **ret_block)
{
  int         n;
  int         status;
  int         pos;
  MPI_Status  mpistatus;

  /* Probe first, because we need to know if our buffer is big enough.
   */
  MPI_Probe(source, tag, comm, &mpistatus);
  MPI_Get_count(&mpistatus, MPI_PACKED, &n);

  /* make sure we are getting the tag we expect and from whom we expect if from */
  if (tag    != MPI_ANY_TAG    && mpistatus.MPI_TAG    != tag) {
    status = eslFAIL;
    goto ERROR;
  }
  if (source != MPI_ANY_SOURCE && mpistatus.MPI_SOURCE != source) {
    status = eslFAIL;
    goto ERROR;
  }

  /* set the source and tag */
  tag = mpistatus.MPI_TAG;
  source = mpistatus.MPI_SOURCE;

  /* Make sure the buffer is allocated appropriately */
  if (*buf == NULL || n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n); 
    *nalloc = n; 
  }

  /* Receive the packed top hits */
  MPI_Recv(*buf, n, MPI_PACKED, source, tag, comm, &mpistatus);

  /* Unpack it - watching out for the EOD signal of M = -1. */
  pos = 0;
  if ((status = seqinfo_block_mpi_unpack(*buf, n, &pos, comm, ret_block)) != eslOK) goto ERROR;
  
  return eslOK;

 ERROR:
  return status;
}

/* Function:  seqinfo_block_mpi_pack_size()
 * Synopsis:  Calculates size needed to pack a sequence block.
 * Incept:    EPN, Mon Jun  6 19:18:15 2011
 *
 * Returns:   <eslOK> on success, and <*ret_n> contains the answer.
 *
 * Throws:    <eslESYS> if an MPI call fails, and <*ret_n> is 0.
 */
int
seqinfo_block_mpi_pack_size(SEQINFO_BLOCK *block, MPI_Comm comm, int *ret_n)
{
  int   status;
  int   n = 0;
  int   sz;

  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed"); n += sz;  /* first_idx  */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed"); n += sz;  /* final_idx  */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed"); n += sz;  /* first_from */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed"); n += sz;  /* final_to   */
  if (MPI_Pack_size(1, MPI_LONG_LONG_INT, comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed"); n += sz;  /* blockL     */
  if (MPI_Pack_size(1, MPI_INT,           comm, &sz) != 0) ESL_XEXCEPTION(eslESYS, "pack size failed"); n += sz;  /* complete  */

  *ret_n = n;
  return eslOK;

 ERROR:
  *ret_n = 0;
  return status;
}

/* Function:  seqinfo_block_mpi_pack()
 * Synopsis:  Packs a SEQINFO_BLOCK into MPI buffer.
 * Incept:    EPN, Mon Jun  6 19:42:01 2011
 *
 * Returns:   <eslOK> on success; <buf> now contains the
 *            packed <hit>, and <*position> is set to the byte
 *            immediately following the last byte of the HIT
 *            in <buf>. 
 *
 * Throws:    <eslESYS> if an MPI call fails; or <eslEMEM> if the
 *            buffer's length <n> was overflowed in trying to pack
 *            <msa> into <buf>. In either case, the state of
 *            <buf> and <*position> is undefined, and both should
 *            be considered to be corrupted.
 */
int
seqinfo_block_mpi_pack(SEQINFO_BLOCK *block, char *buf, int n, int *pos, MPI_Comm comm)
{
  int             status;

  if (MPI_Pack(&block->first_idx,  1, MPI_LONG_LONG_INT, buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&block->final_idx,  1, MPI_LONG_LONG_INT, buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&block->first_from, 1, MPI_LONG_LONG_INT, buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&block->final_to,   1, MPI_LONG_LONG_INT, buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&block->blockL,     1, MPI_LONG_LONG_INT, buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 
  if (MPI_Pack(&block->complete,   1, MPI_INT,           buf, n, pos, comm) != 0) ESL_XEXCEPTION(eslESYS, "pack failed"); 

  if (*pos > n) ESL_EXCEPTION(eslEMEM, "buffer overflow");
  return eslOK;
  
 ERROR:
  return status;
}

/* Function:  seqinfo_block_mpi_unpack()
 * Synopsis:  Unpacks a SEQINFO_BLOCK from an MPI buffer.
 * Incept:    EPN, Mon Jun  6 19:43:49 2011
 *
 * Returns:   <eslOK> on success. <*pos> is updated to the position of
 *            the next element in <buf> to unpack (if any). <*ret_block>
 *            contains a newly allocated SEQINFO_BLOCK, which the caller is
 *            responsible for free'ing.
 *            
 * Throws:    <eslESYS> on an MPI call failure. <eslEMEM> on allocation failure.
 *            In either case, <*ret_hit> is <NULL>, and the state of <buf>
 *            and <*pos> is undefined and should be considered to be corrupted.
 */
int
seqinfo_block_mpi_unpack(char *buf, int n, int *pos, MPI_Comm comm, SEQINFO_BLOCK **ret_block)
{
  int  status;
  SEQINFO_BLOCK *block;

  ESL_ALLOC(block, sizeof(SEQINFO_BLOCK));

  if (MPI_Unpack(buf, n, pos, &block->first_idx,  1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if (MPI_Unpack(buf, n, pos, &block->final_idx,  1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if (MPI_Unpack(buf, n, pos, &block->first_from, 1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if (MPI_Unpack(buf, n, pos, &block->final_to,   1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if (MPI_Unpack(buf, n, pos, &block->blockL,     1, MPI_LONG_LONG_INT, comm) != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 
  if (MPI_Unpack(buf, n, pos, &block->complete,   1, MPI_INT, comm)           != 0) ESL_XEXCEPTION(eslESYS, "mpi unpack failed"); 

  *ret_block = block;

  return eslOK;

 ERROR:
  return status;
}

SEQINFO_BLOCK *
create_seqinfo_block()
{ 
  int status;

  SEQINFO_BLOCK *block = NULL;
  ESL_ALLOC(block, sizeof(SEQINFO_BLOCK));

  block->first_idx  = -1;
  block->final_idx  = -1;
  block->first_from = 0;
  block->final_to   = 0;
  block->blockL     = 0;
  block->complete   = FALSE;

  return block;

 ERROR: 
  if(block != NULL) free(block);
  return NULL;
}

void
dump_seqinfo_block(FILE *fp, SEQINFO_BLOCK *block)
{ 

  fprintf(fp, "SEQINFO_BLOCK:\n");
  fprintf(fp, "%-15s: %ld\n", "first_idx",  block->first_idx);
  fprintf(fp, "%-15s: %ld\n", "final_idx",  block->final_idx);
  fprintf(fp, "%-15s: %ld\n", "first_from", block->first_from);
  fprintf(fp, "%-15s: %ld\n", "final_to",   block->final_to);
  fprintf(fp, "%-15s: %d\n",  "complete",   block->complete);

  return;
}

void
dump_seqinfo_block_list(FILE *fp, SEQINFO_BLOCK_LIST *list)
{ 
  int i;

  fprintf(fp, "SEQINFO_BLOCK_LIST:\n");
  fprintf(fp, "%-15s: %d\n", "N",         list->N);
  fprintf(fp, "%-15s: %d\n", "nalloc",    list->nalloc);
  for(i = 0; i < list->N; i++) { 
    dump_seqinfo_block(fp, list->blocks[i]);
  }

  return;
}
 
SEQINFO_BLOCK_LIST *
create_seqinfo_block_list()
{ 
  int status;

  SEQINFO_BLOCK_LIST *list = NULL;
  ESL_ALLOC(list, sizeof(SEQINFO_BLOCK_LIST));
  
  list->blocks = NULL;
  ESL_ALLOC(list->blocks, sizeof(SEQINFO_BLOCK *) * (100));
  list->N      = 0;
  list->nalloc = 100;

  return list;
  
 ERROR: 
  if(list != NULL) free(list);
  return NULL;
}

void
free_seqinfo_block_list(SEQINFO_BLOCK_LIST *list)
{ 
  int i;
  if(list->blocks != NULL) { 
    for(i = 0; i < list->N; i++) { 
      free(list->blocks[i]);
    }
    free(list->blocks);
  }
  free(list);

  return;
}

/* Function:  add_blocks()
 * Synopsis:  Adds >= 1 new blocks to a SEQINFO_BLOCK_LIST.
 * Incept:    EPN, Tue Jun  7 09:27:34 2011
 * 
 * Purpose:   Determines at least 1 new SEQINFO_BLOCK to send to a MPI
 *            worker by inspecting (via SSI info) least 1 more
 *            sequence in the input file. If the first sequence
 *            inspected contains more than one block, then all blocks
 *            necessary to include that full sequence will be added to
 *            <block_list>. If the first sequence inspected does not
 *            contain enough residues to fill a block, more sequences
 *            will be looked at until at least one block is filled.
 *
 *            It is probable that the final block added to
 *            <block_list> will be incomplete, in which case the next
 *            call to this function will complete it. Its also
 *            possible that <block_list> 
 *     
 *            Likewise, upon entering <block_list> should not be 
 *            empty, but rather contain exactly 1 incomplete block,
 *            that has either been started to be calculated by a 
 *            previous call to next_block(), or that has just been
 *            initialized (which happens when this function is called
 *            for the first time, for example).
 *
 * Returns:   <eslOK> on success and <block_list->blocks[0]> is 
 *            a complete block.
 *            (i.e <block_list->blocks[0]->complete == TRUE>)
 * 
 *            <eslEOF> if we're at the end of a file, in this 
 *            case <block_list->blocks[0]> may still be valid,
 *            if it is non-empty. But it's also possible that 
 *            there was no sequence left, in which case 
 *            <block_list> will be freed and set as NULL.
 * 
 *            If we return any other status code besides eslOK or 
 *            eslEOF, errbuf will be filled and caller should
 *            exit with mpi_failure().
 */
int 
add_blocks(ESL_SQFILE *dbfp, ESL_SQ *sq, int64_t ncontext, char *errbuf, SEQINFO_BLOCK_LIST *block_list, int64_t *ret_pkey_idx, int64_t *ret_noverlap, int64_t *ret_nseq)
{
  int status = eslOK; /* error code status */
  int64_t pkey_idx;
  int64_t noverlap = 0;
  int64_t tot_noverlap = 0;
  int64_t nseq = 0;

  /* Contract check, SSI should be valid */
  if(dbfp->data.ascii.ssi == NULL) ESL_FAIL(status, errbuf, "No SSI index available (it should've been opened earlier)");
  /* we should have exactly 1 non-complete block */
  if(block_list->N != 1 || block_list->blocks == NULL || block_list->blocks[0]->complete) { 
    ESL_FAIL(eslFAIL, errbuf, "contract violation in next_block, block_list does not contain exactly 1 incomplete block");
  }

  status = eslOK;
  pkey_idx = *ret_pkey_idx;
  /* read sequences until we have at least 1 complete block to send to a worker */
  while((pkey_idx < dbfp->data.ascii.ssi->nprimary) && (! block_list->blocks[0]->complete)) { 
    status = inspect_next_sequence_using_ssi(dbfp, sq, ncontext, pkey_idx, errbuf, block_list, &noverlap);
    nseq++;
    if(block_list == NULL) ESL_FAIL(eslFAIL, errbuf, "error parsing database in add_blocks()");
    pkey_idx++;
    tot_noverlap += noverlap;
  }
  if(pkey_idx == dbfp->data.ascii.ssi->nprimary && (! block_list->blocks[block_list->N-1]->complete)) { 
    /* We reached the end of the file and didn't finish the final
     * block.  If its non-empty then it's valid, else it's not. */
    if(block_list->blocks[block_list->N-1]->blockL > 0) { /* valid */
      block_list->blocks[block_list->N-1]->complete = TRUE; 
    }
    else {
      /* final block was initialized but doesn't correspond to any sequence yet, destroy it */
      free(block_list->blocks[block_list->N-1]);
      block_list->blocks[block_list->N-1] = NULL;
      block_list->N--;
    }
  }

  *ret_pkey_idx = pkey_idx;
  *ret_noverlap = tot_noverlap;
  *ret_nseq     = nseq;

  return status;

}

/* Function:  inspect_next_sequence_using_ssi()
 * Synopsis:  Inspect a database sequence and update SEQINFO_BLOCK_LIST.
 * Incept:    EPN, Tue Jun  7 07:41:02 2011
 *
 * Purpose:   Uses SSI indexing to gather information on the next
 *            sequence in <dbfp> and add that information to a 
 *            <block_list>. This will update the information in
 *            at least one block in the list, potentially finishing
 *            one block and creating and potentially finishing 
 *            more blocks, depending on the length of the sequence.
 *
 * Returns:   <eslOK> on success.
 *            <eslEOF> if we're at the end of a file.
 *            In both cases block_list will be updated.
 * 
 *            If we return any other status code besides eslOK or 
 *            eslEOF, errbuf will be filled and caller should
 *            exit with mpi_failure().

 */
int 
inspect_next_sequence_using_ssi(ESL_SQFILE *dbfp, ESL_SQ *sq, int64_t ncontext, int64_t pkey_idx, char *errbuf, SEQINFO_BLOCK_LIST *block_list, int64_t *ret_noverlap)
{
  int     status;           /* error code status */
  int64_t L;                /* length of sequence */
  int64_t ntoadd;           /* max length of sequence we could still add to a incomplete block */
  int64_t seq_from;         /* start index of subsequence */
  int64_t nremaining;       /* number of residues for the sequence not yet assigned to any block */
  int64_t noverlap = 0;     /* number of residues that exist in more than one block due to overlaps */
  SEQINFO_BLOCK *cur_block = NULL; /* pointer to current block we are working on */        

  /* Contract check */
  if(dbfp->data.ascii.ssi == NULL) ESL_FAIL(status, errbuf, "No SSI index available (it should've been opened earlier)");

  /* Get length of 'next' sequence (next sequence in list of 
   * primary keys in SSI index, probably not the next sequence 
   * in the file) */
  status = esl_ssi_FindNumber(dbfp->data.ascii.ssi, pkey_idx, NULL, NULL, NULL, &L, NULL);
  if(status == eslENOTFOUND) ESL_FAIL(status, errbuf, "unable to find sequence %ld in SSI index file, try re-indexing with esl-sfetch.", pkey_idx);
  if(status == eslEFORMAT)   ESL_FAIL(status, errbuf, "SSI index for database file is in incorrect format.");
  if(status != eslOK)        ESL_FAIL(status, errbuf, "proble with SSI index for database file.");

  /* Create a new block list if necessary */
  if(block_list == NULL) { 
    /* This should only happen on first call to this function, or if
     * by chance the previous sequence (inspected with the previous
     * call to next_sequence_using_ssi() filled its final block
     * exactly full, with no remaining residues with which to start a
     * new block. */
    block_list = create_seqinfo_block_list();
    /* add first block */
    block_list->blocks[0] = create_seqinfo_block();
    block_list->N++;
  }
    
  /* Now fill in >= 1 blocks with this sequence.
   * Note: we fill blocks to be exactly CMSEARCH_MAX_RESIDUE_COUNT
   * whenever possible (i.e. only the final block including the end of
   * the final sequence can not equal this size). We could be smarter
   * about not allowing stupid overlaps (i.e. including 1..L-1 in one
   * block and L-ncontext+1..L in the next block), but that's not
   * how serial mode operates (it relies of esl_sqio_ReadWindow()
   * and is ignorant of how much sequence is left after the current
   * window), so we do it the same way here for consistency. Plus
   * this way is simpler to implement and understand. 
   */
  cur_block  = block_list->blocks[0];
  nremaining = L;

  while(nremaining > 0) { 
    ntoadd   = CMSEARCH_MAX_RESIDUE_COUNT - cur_block->blockL; /* max number residues we can add to cur_block */
    ntoadd   = ESL_MIN(ntoadd, nremaining);
    if(cur_block->blockL == 0) { /* this block was just created */
      cur_block->first_idx = pkey_idx; 
    }
    cur_block->final_idx = pkey_idx; /* this may be changed by a subsequent call to this function */

    if(nremaining <= ntoadd) { 
      /* Case 1: We can fit all remaining sequence in current block */
      if(cur_block->blockL == 0) { 
	cur_block->first_from = (L - nremaining) + 1;
      }
      cur_block->final_to  = L;
      cur_block->blockL   += nremaining;
      if(cur_block->blockL == CMSEARCH_MAX_RESIDUE_COUNT) { 
	/* rare (nremaining == ntoadd), but certainly possible */
	cur_block->complete = TRUE;
	/* we'll make a new block below (under the else() below 
	 * this block will remain in its initial state until
	 * we enter this function the next time. */
      }
      /* update nremaining, this will break us out of the while()  */
      nremaining = 0;
    }
    else { 
      /* Case 2: We can't fit all of the remaining sequence in current
       *         block, add as much of the remaining sequence as we
       *         can to the current block, and then make a new block
       *         for the remaining sequence (new block is made below,
       *         outside of this else).
       */
      seq_from = L - nremaining + 1;
      if(cur_block->blockL == 0) { 
	cur_block->first_from = seq_from;
      }
      cur_block->final_to  = seq_from + ntoadd - 1;
      cur_block->blockL   += ntoadd;
      cur_block->complete  = TRUE;

      /* update nremaining */
      nremaining -= (ntoadd - ncontext); /* <noverlap> residues will overlap between this block and the previous one */
      noverlap   += ncontext;
    }
    
    /* create a new block if necessary */
    if(cur_block->complete) { 
      if(block_list->N == block_list->nalloc) { 
	ESL_REALLOC(block_list->blocks, sizeof(SEQINFO_BLOCK *) * (block_list->nalloc+100));
	block_list->nalloc += 100;
      }
      cur_block = create_seqinfo_block();
      block_list->blocks[block_list->N++] = cur_block;
    }
  }
  *ret_noverlap = noverlap;

  return status;

 ERROR:
  ESL_FAIL(status, errbuf, "out of memory");
}

#endif


/*****************************************************************
 * @LICENSE@
 *****************************************************************/

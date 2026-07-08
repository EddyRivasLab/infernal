/* cmalign: align sequences to a CM.
 * 
 * EPN, Fri Dec 30 10:13:34 2011 [Updated for v1.1]
 * SRE, Thu Jul 25 11:28:03 2002 [St. Louis]
 */

#include <esl_config.h>
#include <p7_config.h>
#include "config.h"	

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <inttypes.h>

#include "easel.h"		/* general seq analysis library   */
#include "esl_alphabet.h"
#include "esl_getopts.h"		
#include "esl_mpi.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_msafile2.h"
#include "esl_msaweight.h"
#include "esl_random.h"		
#include "esl_sq.h"		
#include "esl_sqio.h"
#include "esl_stack.h"
#include "esl_sse.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

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

#include "infernal.h"

/* Max number of sequences per tmp alignment, if seq file exceeds
 * either of these, final output alignment will be in 1 line/seq Pfam
 * format. Max number of residues is CM_MAX_RESIDUE_COUNT from infernal.h
 * where it's currently defined as (1024 * 1024).
 */
#define CMALIGN_MAX_NSEQ  10000  /* 10k sequences, average parsetree is 25 bytes/position this means ~250Mb for all parsetrees */

#define DEBUGSERIAL 0
#define DEBUGMPI    0

/* Brief 26_0526-017 (Part A): scale-aware default for the IBV D&C base-case slab
 * on the --hmm --p7ibv path. When the user has NOT set --p7ibv-base-slab,
 * the deriver otherwise picks the adaptive 256 MB-capped slab, which is
 * far above the memory knee at common scale (brief 26_0526-015: base_slab 1024 for
 * norovirus, 372 for sars). Brief 26_0526-015's sweep showed base_slab ~= 64 is the
 * memory knee at common scale (norovirus 237->70 MB, sars 476->265 MB for
 * ~+0.4-1.6 s wall) AND is at/below the adaptive value at genome scale
 * (HSV adaptive ~72, and 64 is below the matrix peak so peak RSS is
 * unchanged there). A flat 64 default is therefore robust across scales.
 * The deriver's OUTPUT is byte-invariant to base_slab (brief 26_0526-017 gate A1),
 * so this is a silent memory-only default; an explicit --p7ibv-base-slab
 * overrides it. */
#define HMM_P7IBV_KNEE_BASE_SLAB 64

typedef struct {
#ifdef HMMER_THREADS
  ESL_WORK_QUEUE   *queue;
#endif /*HMMER_THREADS*/
  CM_t             *cm;          /* a covariance model */
  CM_ALNDATA      **dataA;       /* array of CM_ALNDATA objects with ptrs to sqs, parsetrees, scores */
  int               n;           /* size of outdataA   */
  float             mxsize;      /* max size (Mb) of allowable DP mx */
  int               pass_idx;    /* pipeline pass index, controls truncation bit sc penalty */
  ESL_STOPWATCH    *w;           /* stopwatch for timing stages (band calc, alignment) */
  ESL_STOPWATCH    *w_tot;       /* stopwatch for timing total time for processing 1 seq */
  int               do_failover; /* TRUE if we're trying to do HMM banded truncated alignment,
				  * and bands obscure all possible alignments (very rare) to
				  * failover into HMM banded standard alignment.
				  */
  /* HMM-only alignment fields (--hmm mode, used by hmm_pipeline_thread) */
  P7_PROFILE       *gm;           /* thread-local p7 profile (NULL if not --hmm) */
  P7_BG            *bg;           /* brief 26_0430-182: thread-local null model, needed to re-run p7_ProfileConfig() per-seq under Tgm (NULL if not --hmm) */
  P7_HMM           *hmm;          /* ptr to p7 HMM, shared read-only (NULL if not --hmm) */
  P7_GMX           *gx;           /* Viterbi DP matrix (NULL if not --hmm or unbanded-only) */
  P7_GMX           *gxf;          /* Forward matrix (NULL unless --hmm --hmmnoband) */
  P7_GMX           *gxb;          /* Backward matrix (NULL unless --hmm --hmmnoband) */
  P7_TRACE        **hmm_tr;       /* shared trace array; worker writes tr[seqidx] (NULL if not --hmm) */
  int               do_hmmvit;    /* TRUE for --hmm --hmmvit */
  int               do_hmmnoband; /* TRUE for --hmm --hmmnoband */
  int               do_p7ibv;     /* TRUE for --hmm --p7ibv (banded OA via IBV deriver) */
  int               ibv_delta;    /* IBV Delta milli-bits (--p7ibv-delta) */
  int               ibv_base_slab;/* IBV D&C base-case slab; 0=auto (--p7ibv-base-slab) */
  int               do_trunc;     /* brief 26_0430-182: TRUE if CM_ALIGN_TRUNC set (drives Tgm gm config + IBV deriver do_trunc arg) */
} WORKER_INFO;

#define ACCOPTS      "--hbanded,--nonbanded,--p7band"         /* Exclusive choice for acceleration or not */
#define ALGOPTS      "--cyk,--optacc,--sample"               /* Exclusive choice for algorithm */
#if defined (HMMER_THREADS) && defined (HAVE_MPI)
#define CPUOPTS     "--mpi"
#define MPIOPTS     "--cpu"
#else
#define CPUOPTS     NULL
#define MPIOPTS     NULL
#endif

static ESL_OPTIONS options[] = {
  /* name                   type       default env          range    toggles         reqs         incomp  help  docgroup*/
  { "-h",            eslARG_NONE,       FALSE, NULL,        NULL,       NULL,        NULL,          NULL, "show brief help on version and usage",               1 },
  { "-o",         eslARG_OUTFILE,        NULL, NULL,        NULL,       NULL,        NULL,          NULL, "output the alignment to file <f>, not stdout",       1 },
  { "-g",            eslARG_NONE,       FALSE, NULL,        NULL,       NULL,        NULL,          NULL, "configure CM for global alignment [default: local]", 1 },
  /* options controlling the alignment algorithm */
  { "--optacc",      eslARG_NONE,   "default", NULL,        NULL,    ALGOPTS,        NULL,     "--small", "use the Holmes/Durbin optimal accuracy algorithm  [default]",     2 },
  { "--cyk",         eslARG_NONE,       FALSE, NULL,        NULL,    ALGOPTS,        NULL,          NULL, "use the CYK algorithm",                                           2 },
  { "--sample",      eslARG_NONE,       FALSE, NULL,        NULL,    ALGOPTS,        NULL,     "--small", "sample alignment of each seq from posterior distribution",        2 },
  { "--seed",         eslARG_INT,       "181", NULL,      "n>=0",       NULL,  "--sample",          NULL, "w/--sample, set RNG seed to <n> (if 0: one-time arbitrary seed)", 2 },
  { "--notrunc",     eslARG_NONE,       FALSE, NULL,        NULL,       NULL,        NULL,          NULL, "do not use truncated alignment algorithm",                        2 },
  { "--sub",         eslARG_NONE,       FALSE, NULL,        NULL,       NULL,"--notrunc,-g",        NULL, "build sub CM for columns b/t HMM predicted start/end points",     2 },
  { "--hmm",         eslARG_NONE,       FALSE, NULL,        NULL,       NULL,        NULL, "--sub,--small", "use the p7 HMM only to align (no CM alignment)",                  2 },
  { "--hmmvit",      eslARG_NONE,       FALSE, NULL,        NULL,       NULL,   "--hmm", "--hmmnoband", "w/--hmm, use Viterbi traces (faster, less accurate)",               2 },
  { "--hmmnoband",   eslARG_NONE,       FALSE, NULL,        NULL,       NULL,   "--hmm",   "--hmmvit", "w/--hmm, do not use Viterbi bands for OA alignment",              2 },
  /* options affecting speed and memory */
  { "--hbanded",     eslARG_NONE,   "default", NULL,        NULL,    ACCOPTS,        NULL,                     NULL, "accelerate using CM plan 9 HMM derived bands",               3 },
  { "--tau",         eslARG_REAL,      "1e-7", NULL, "1e-18<x<1",       NULL,        NULL,            "--nonbanded", "set tail loss prob for HMM bands to <x>",                    3 },
  { "--mxsize",      eslARG_REAL,    "1024.0", NULL,      "x>0.",       NULL,        NULL,                     NULL, "set maximum allowable DP matrix size to <x> Mb",             3 },
  { "--fixedtau",    eslARG_NONE,       FALSE, NULL,        NULL,       NULL,        NULL,            "--nonbanded", "do not adjust tau (tighten bands) until mx < limit", 3 },
  { "--maxtau",      eslARG_REAL,      "0.05", NULL,   "0<x<0.5",       NULL,        NULL, "--fixedtau,--nonbanded", "set max tau <x> when tightening HMM bands",                  3 },
  { "--nonbanded",   eslARG_NONE,       FALSE, NULL,        NULL,    ACCOPTS,        NULL,                     NULL, "do not use HMM bands for faster alignment",                  3 },
  { "--p7band",      eslARG_NONE,       FALSE, NULL,        NULL,    ACCOPTS,        NULL,                     NULL, "use p7 Viterbi-derived bands for faster alignment",          3 },
  { "--p7padplus",    eslARG_INT,         "7", NULL,      "n>=0",       NULL,   "--p7band",                    NULL, "add <n> to every per-node p7 band pad [default 7]",          3 },
  { "--p7pinbridge", eslARG_NONE,       FALSE, NULL,        NULL,       NULL,   "--p7band",                    NULL, "SW-pinbridge prefilter + banded p7 Viterbi (--p7band)", 3 },
  { "--p7pbpad",      eslARG_INT,        "20", NULL,      "n>=0",       NULL, "--p7pinbridge",                 NULL, "diagonal pad for SW-pinbridge prefilter [default 20]",  3 },
  { "--p7pinbridge-vitgaps", eslARG_NONE, FALSE, NULL,     NULL,       NULL, "--p7pinbridge",                 NULL, "exact mini-Viterbi gap costs in gap-aware LSIS (Opt 3)", 3 },
  { "--p7ibv",       eslARG_NONE,       FALSE, NULL,        NULL,       NULL,        NULL,  "--p7pinbridge", "use F+B direct-band derivation (w/--p7band or --hmm)",        3 },
  { "--p7ibv-delta",  eslARG_INT,     "20000", NULL,      "n>=0",       NULL,     "--p7ibv",              NULL, "IBV Delta milli-bits",                                       3 },
  { "--p7ibv-mode",  eslARG_STRING, "delta", NULL,        NULL,       NULL,     "--p7ibv",              NULL, "IBV band mode: delta|fixed|hybrid (brief 140)",             3 },
  { "--p7ibv-width", eslARG_INT,       "20", NULL,      "n>=0",       NULL,     "--p7ibv",              NULL, "fixed-width pad W around argmax-k pin (fixed/hybrid)",       3 },
  { "--p7ibv-mem",   eslARG_NONE,       FALSE, NULL,        NULL,       NULL,     "--p7ibv",              NULL, "use D&C O(M*logL) band deriver (brief 124)",                 3 },
  { "--p7ibv-base-slab", eslARG_INT,      "0", NULL,      "n>=0",       NULL, "--p7ibv-mem",              NULL, "D&C base-case slab size; 0=auto (mem-capped)",               3 },
  { "--p7ibv-ckpt",  eslARG_NONE,       FALSE, NULL,        NULL,       NULL, "--p7ibv-mem",              NULL, "checkpoint Pass-2 banded CP9 F/B (low mem; brief 146)",      3 },
  { "--p7ibv-wv",    eslARG_NONE,       FALSE, NULL,        NULL,       NULL,     "--p7ibv",              NULL, "windowed-Viterbi band: i2k +/- F+B-halfwidth pad (brief 169)",3 },
  { "--p7wv-nsamp",  eslARG_INT,        "40", NULL,       "n>0",       NULL, "--p7ibv-wv",              NULL, "WV pad calibration: # CM-emitted samples",                  3 },
  { "--p7wv-q",      eslARG_REAL,     "0.99", NULL,    "0<x<=1",       NULL, "--p7ibv-wv",              NULL, "WV pad calibration: half-width quantile",                   3 },
  { "--p7wv-floor",  eslARG_INT,         "2", NULL,      "n>=0",       NULL, "--p7ibv-wv",              NULL, "WV pad calibration: floor pad",                             3 },
  { "--p7wv-seed",   eslARG_INT,       "181", NULL,      "n>=0",       NULL, "--p7ibv-wv",              NULL, "WV pad calibration: RNG seed",                              3 },
  { "--p7wvpad-dump",eslARG_OUTFILE,   NULL,  NULL,        NULL,       NULL, "--p7wv-calib","--p7wvpad-file", "brief172: dump calibrated WV pad to <f> (amortize calib)",   3 },
  { "--p7wvpad-file",eslARG_INFILE,    NULL,  NULL,        NULL,       NULL, "--p7ibv-wv",   "--p7wv-nsamp,--p7wv-pad", "brief172: load WV pad from <f> (skip per-run calib)",  3 },
  { "--p7wv-pad",    eslARG_INT,        "30", NULL,      "n>=0",       NULL, "--p7ibv-wv",   "--p7wv-calib", "brief173: constant WV band half-width (no calibration)",     3 },
  { "--p7wv-calib",  eslARG_NONE,       FALSE, NULL,        NULL,       NULL, "--p7ibv-wv",   "--p7wvpad-file", "brief173: opt back in to per-node WV pad calibration",       3 },
  { "--p7kmeranchor",eslARG_NONE,       FALSE, NULL,        NULL,       NULL,        NULL, "--p7ibv,--p7pinbridge,--p7kmerchain", "k-mer best-window anchor bands (--p7band/--hmm)",   3 },
  { "--p7kmerchain", eslARG_NONE,       FALSE, NULL,        NULL,       NULL,        NULL, "--p7ibv,--p7pinbridge,--p7kmeranchor", "genome-scale k-mer seed+chain bands (--p7band/--hmm)", 3 },
  { "--p7kmerchain-alpha", eslARG_REAL, "0.75", NULL,      "x>=0",       NULL, "--p7kmerchain",                   NULL, "brief 043: kmerchain ramp-slack alpha [default 0.75]",       3 },
  { "--p7kmerchain-mink", eslARG_INT,      "0", NULL,      "n>=0",       NULL,        NULL,                     NULL, "brief 046: gate kmeranchor/kmerchain if k>=<n> tier finds 0 hits [default 0=off]", 3 },
  { "--p7kmerchain-mgate", eslARG_INT,     "0", NULL,      "n>=0",       NULL,        NULL,                     NULL, "brief 26_0628-047: gate kmeranchor/kmerchain if M < <n> [default 0=off]",          3 },
  { "--p7kmerchain-fbvit", eslARG_NONE, FALSE, NULL,  NULL,       NULL,        NULL,                     NULL, "brief 26_0628-047: gate fallback uses old Vit-trace band, not --p7ibv",           3 },
  { "--cykbands",    eslARG_NONE,       FALSE, NULL,        NULL,       NULL,   "--p7band",                    NULL, "run CYK pre-pass and tighten bands before Inside/Outside",   3 },
  { "--cykpad",       eslARG_INT,         "2", NULL,      "n>=0",       NULL,  "--cykbands",                   NULL, "pad <n> for parsetree band tightening [default 2]",  3 },
  { "--cykskip-unvisited", eslARG_NONE, FALSE, NULL,        NULL,       NULL,  "--cykbands",                   NULL, "skip CM states not visited by CYK parsetree (aggressive)",    3 },
  { "--dump-bands",    eslARG_OUTFILE,     NULL, NULL,        NULL,       NULL,   "--p7band",                    NULL, "dump per-(state,j) band TSV to <f> before cm_AlignHB",      3 },
  { "--small",       eslARG_NONE,       FALSE, NULL,        NULL,       NULL,        NULL,                "--mxsize", "use small memory divide and conquer (d&c) algorithm",       3 },  /* for --small, required opts are enforced below */
  { "--ckpt",        eslARG_NONE,       FALSE, NULL,        NULL,       NULL,        "-g","--cyk,--sample,--nonbanded,--small,--sub", "use checkpointed sqrt(M)-memory HMM-banded optacc engines", 3 },
  /* options controlling optional output */
  { "--sfile",    eslARG_OUTFILE,        NULL, NULL,        NULL,       NULL,        NULL,          NULL, "dump alignment score information to file <f>",            4 },
  { "--tfile",    eslARG_OUTFILE,        NULL, NULL,        NULL,       NULL,        NULL,          NULL, "dump individual sequence parsetrees to file <f>",         4 },
  { "--ifile",    eslARG_OUTFILE,        NULL, NULL,        NULL,       NULL,        NULL,          NULL, "dump information on per-sequence inserts to file <f>",    4 },
  { "--elfile",   eslARG_OUTFILE,        NULL, NULL,        NULL,       NULL,        NULL,          "-g", "dump information on per-sequence EL inserts to file <f>", 4 },
  /* other expert options */
  { "--mapali",    eslARG_INFILE,        NULL, NULL,        NULL,       NULL,        NULL,          NULL, "include alignment in file <f> (same ali that CM came from)", 5 },
  { "--mapstr",      eslARG_NONE,        NULL, NULL,        NULL,       NULL,  "--mapali",          NULL, "include structure (w/pknots) from <f> from --mapali <f>",    5 },
  { "--noss",        eslARG_NONE,        NULL, NULL,        NULL,       NULL,  "--mapali",    "--mapstr", "cmbuild --noss option was used w/aln from --mapali <f>",     5 },
  { "--informat",  eslARG_STRING,        NULL, NULL,        NULL,       NULL,        NULL,          NULL, "assert <seqfile> is in format <s>: no autodetection",        5 },
  { "--outformat", eslARG_STRING, "Stockholm", NULL,        NULL,       NULL,        NULL,          NULL, "output alignment in format <s>",                             5 },
  { "--dnaout",      eslARG_NONE,       FALSE, NULL,        NULL,       NULL,        NULL,          NULL, "output alignment as DNA (not RNA) sequence data",            5 },
  { "--noprob",      eslARG_NONE,       FALSE, NULL,        NULL,       NULL,        NULL,          NULL, "do not include posterior probabilities in the alignment",    5 },
  { "--matchonly",   eslARG_NONE,       FALSE, NULL,        NULL,       NULL,        NULL,          NULL, "include only match columns in output alignment",             5 },
  { "--miss",        eslARG_NONE,       FALSE, NULL,        NULL,       NULL,        NULL,          NULL, "mark seqs w/terminal gaps as fragments w/missing (~) chars", 5 },
  { "--ileaved",     eslARG_NONE,       FALSE, NULL,        NULL,       NULL,        NULL, "--outformat","force output in interleaved Stockholm format",                5 },
  { "--flanktoins",  eslARG_REAL,        NULL, NULL,   "0<x<0.4",       NULL,"--flankselfins",      NULL, "change transition probs into ROOT_IL/IR to <x> (e.g. 0.1)",  5 }, 
  { "--flankselfins",eslARG_REAL,        NULL, NULL,   "0<x<0.9",       NULL,"--flanktoins",        NULL, "change self transit probs for ROOT_IL/IR to <x> (e.g. 0.8)", 5 }, 
  { "--regress",  eslARG_OUTFILE,        NULL, NULL,        NULL,       NULL, "--ileaved",    "--mapali", "save regression test data to file <f>",                      5 }, 
  { "--verbose",     eslARG_NONE,       FALSE, NULL,        NULL,       NULL,        NULL,          NULL, "report extra information; mainly useful for debugging",      5 },
  /*{ "--noannot",   eslARG_NONE,       FALSE, NULL,        NULL,       NULL,        NULL,          NULL, "do not add cmalign execution annotation to the alignment",   5 },*/
#ifdef HMMER_THREADS 
  { "--cpu",          eslARG_INT,      CMNCPU, "INFERNAL_NCPU","n>=0",  NULL,        NULL,       CPUOPTS, "number of parallel CPU workers to use for multithreads",     5 },
#endif
#ifdef HAVE_MPI
  { "--mpi",         eslARG_NONE,       FALSE, NULL,        NULL,       NULL,        NULL,       MPIOPTS, "run as an MPI parallel program",                             5 },  
  { "--stall",       eslARG_NONE,       FALSE, NULL,        NULL,       NULL,        NULL,          NULL, "arrest after start: for debugging MPI under gdb",            5 },  
#endif
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

struct cfg_s {
  char            *cmfile;      /* name of input CM file  */ 
  char            *sqfile;	/* name of sequence file  */ 
  CM_FILE         *cmfp;	/* open input CM file stream       */
  ESL_SQFILE      *sqfp;        /* open sequence input file stream */
  ESL_ALPHABET    *abc;         /* alphabet for input */
  ESL_ALPHABET    *abc_out;     /* alphabet for output */
  int              infmt;       /* input alignment format */
  int              outfmt;      /* output alignment format */
  int              be_verbose;  /* TRUE if --verbose used */
  int              do_oneblock; /* TRUE to force output of full alignment 
				 * in one block, if input file is really big,
				 * we'll fail and tell the user to pick a
				 * different format.
				 */
  /* mpi */
  int              do_mpi;      
  int              my_rank;
  int              nproc;
  int              do_stall;    /* TRUE to stall the program until gdb attaches */

  /* Masters only */
  FILE            *tmpfp;	/* the temporary output file where alignments are initially written if !do_oneblock */
  FILE            *ofp;	        /* output file where alignments are ultimately written (default is stdout) */
  FILE            *tfp;         /* optional output for parsetrees  */
  FILE            *ifp;	        /* optional output for insert info */
  FILE            *efp;	        /* optional output for EL insert info */
  FILE            *sfp;         /* optional output for alignment scores */
  FILE            *rfp;         /* optional output for --regress alignment */
};

/* brief 26_0628-035: temporary peak-RSS attribution instrumentation, gated by
 * BRIEF035_MEMPOINT. Reads /proc/self/status VmRSS. Revert before finishing
 * if not worth keeping (duplicated from cm_p7_band.c's static copy). */
static long
brief035_rss_kb(void)
{
  FILE *fp = fopen("/proc/self/status", "r");
  char line[256];
  long rss = -1;
  if (fp == NULL) return -1;
  while (fgets(line, sizeof(line), fp) != NULL) {
    if (strncmp(line, "VmRSS:", 6) == 0) { sscanf(line+6, "%ld", &rss); break; }
  }
  fclose(fp);
  return rss;
}

static char usage[]  = "[-options] <cmfile> <seqfile>";
static char banner[] = "align sequences to a CM";

/* brief 26_0628-047: shared kmer-gate (M-gate/N-gate/no-anchor) fallback deriver,
 * used by all 3 kmeranchor/kmerchain call sites below (serial, threaded
 * worker, MPI worker) in place of the old hardcoded p7_Seq2BandsVit-shaped
 * fallback. Defaults to --p7ibv's D&C deriver (p7_Seq2BandsIBV_dnc), the
 * same entry point --hmm --p7ibv itself calls, with --p7ibv-mode/-width's
 * own CLI defaults (cm->p7_ibv_mode/width are guaranteed to still hold
 * their cm.c defaults here, since --p7ibv and --p7kmeranchor/--p7kmerchain
 * are mutually exclusive CLI options -- see the "reqs"/"incompat" fields on
 * the --p7kmeranchor/--p7kmerchain option lines).
 * ibv_delta/ibv_base_slab are passed in already resolved: ibv_base_slab to
 * the brief-017 memory knee, exactly as the --hmm --p7ibv call sites resolve
 * it (the caller may not have `go`, e.g. the threaded worker). ibv_delta is
 * cm->p7_ibv_delta (struct default 3000), deliberately NOT --p7ibv-delta's
 * own CLI default (20000, what the direct --hmm --p7ibv call sites use) --
 * brief 26_0628-050 found 20000 provably worse than <=10000 on mir-2807 and
 * SNORA16 (a higher-forward-score but structurally wrong HMM registration
 * only becomes reachable at wide deltas); cm->p7_ibv_delta is never
 * CLI-overridden here since --p7ibv-delta requires --p7ibv, which is
 * mutually exclusive with --p7kmeranchor/--p7kmerchain. Returns ncells=0
 * (not a hard failure) if
 * cm->fp7 is unusable, mirroring the derivers' own ncells==0 "fall back
 * further" convention -- caller should still fall through to the old
 * unbanded-OA Forward/Backward safety net in that vanishingly rare case. */
static int
kmer_gate_p7ibv_fallback(CM_t *cm, char *errbuf, const ESL_DSQ *dsq, int L, int do_trunc,
                          int ibv_delta, int ibv_base_slab,
                          int **ret_i2k, int **ret_kmin, int **ret_kmax, int *ret_ncells)
{
  *ret_i2k = NULL; *ret_kmin = NULL; *ret_kmax = NULL; *ret_ncells = 0;
  if (cm->fp7 == NULL) return eslOK;   /* caller falls back further */
  return p7_Seq2BandsIBV_dnc(cm, errbuf, dsq, L, ibv_delta, ibv_base_slab,
                              FALSE, /* do_boundary_widen: P135B_FORCE_WIDEN override, same default as --hmm --p7ibv */
                              FALSE, /* do_kband: unbanded D&C, same as --hmm --p7ibv */
                              do_trunc, cm->p7_ibv_mode, cm->p7_ibv_width,
                              ret_i2k, ret_kmin, ret_kmax, ret_ncells);
}

static void serial_master(ESL_GETOPTS *go, struct cfg_s *cfg);
static int  serial_loop  (WORKER_INFO *info, char *errbuf, ESL_SQ_BLOCK *sq_block, ESL_RANDOMNESS *r);
static void hmm_alignment(ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm);
static void output_hmm_insert_info(FILE *ifp, CM_t *cm, P7_HMM *hmm, ESL_SQ **sqarr, P7_TRACE **tr, int nseq);

#ifdef HMMER_THREADS
static int  thread_loop(WORKER_INFO *info, char *errbuf, ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQ_BLOCK *sq_block);
static void pipeline_thread(void *arg);
static int  hmm_thread_loop(WORKER_INFO *info, ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQ **sqarr, int nseq);
static void hmm_pipeline_thread(void *arg);
#endif /*HMMER_THREADS*/

#if HAVE_MPI 
static int  mpi_master   (ESL_GETOPTS *go, struct cfg_s *cfg);
static int  mpi_worker   (ESL_GETOPTS *go, struct cfg_s *cfg);
static void mpi_failure  (char *format, ...);
#define INFERNAL_ERROR_TAG          1
#define INFERNAL_DSQ_TAG            2
#define INFERNAL_INITIALREADY_TAG   3
#define INFERNAL_ALNDATA_TAG        4
#endif

/* Functions to avoid code duplication for common tasks */
static void process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_cmfile, char **ret_sqfile, int *ret_infmt, int *ret_outfmt);
static int  output_header(FILE *ofp, const ESL_GETOPTS *go, char *cmfile, char *sqfile, CM_t *cm, int ncpus);
static int  init_master_cfg (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int  init_shared_cfg (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int  initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int  map_alignment(const char *msafile, CM_t *cm, int noss_used, char *errbuf, CM_ALNDATA ***ret_dataA, int *ret_ndata, char **ret_ss);
static int  output_alignment(ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, FILE *ofp, CM_ALNDATA **dataA, int ndata, char *map_sscons);
static void output_info_file_header(FILE *fp, char *firstline, char *elstring);
static int  output_scores(FILE *ofp, CM_t *cm, char *errbuf, CM_ALNDATA **dataA, int ndata, int first_idx, int be_verbose);
/*static int  add_annotation_to_msa(ESL_GETOPTS *go, char *errbuf, ESL_MSA *msa);*/

/* Functions that enable memory efficiency by storing only a fraction
 * of the seqs/parsetrees from target file in memory at once.
 */
static int  create_and_output_final_msa(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nali, char *tmpfile);
static void update_maxins_and_maxel(ESL_MSA *msa, int clen, int64_t alen, int *maxins, int *maxel);
static int  determine_gap_columns_to_add(ESL_MSA *msa, int *maxins, int *maxel, int clen, int **ret_ngap_insA, int **ret_ngap_elA, int **ret_ngap_eitherA, char *errbuf);
static void inflate_gc_with_gaps_and_els(FILE *ofp, ESL_MSA *msa, int *ngap_insA, int *ngap_elA, char **ret_ss_cons2print, char **ret_rf2print);
static void configure_root_inserts(CM_t *cm, float to_insert_prob, float self_insert_prob);

int
main(int argc, char **argv)
{
  int              status   = eslOK;

  ESL_GETOPTS     *go  = NULL;    /* command line processing                 */
  struct cfg_s     cfg;           /* configuration data                      */

  /* start stopwatch */
  ESL_STOPWATCH   *w   = NULL;    /* for overall timing                      */
  if((w = esl_stopwatch_Create()) == NULL) cm_Fail("out of memory, trying to create stopwatch");
  esl_stopwatch_Start(w);

  /* Set processor specific flags */
  impl_Init();

  /* Initialize what we can in the config structure (without knowing the alphabet yet)
   */
  cfg.cmfile      = NULL;
  cfg.sqfile      = NULL;
  cfg.cmfp        = NULL; 
  cfg.sqfp        = NULL; 
  cfg.do_mpi      = FALSE;               /* this gets reset below, if we init MPI */
  cfg.nproc       = 0;                   /* this gets reset below, if we init MPI */
  cfg.my_rank     = 0;                   /* this gets reset below, if we init MPI */
  cfg.abc         = NULL; 

  cfg.tmpfp       = NULL;	         /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.ofp         = NULL;	         /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.tfp         = NULL;	         /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.ifp         = NULL;	         /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.efp         = NULL;	         /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.sfp         = NULL;	         /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.rfp         = NULL;	         /* opened in init_master_cfg() in masters, stays NULL for workers */
 

  cfg.infmt       = eslSQFILE_UNKNOWN;    /* reset below in process_commandline() */
  cfg.outfmt      = eslMSAFILE_STOCKHOLM; /* reset below in process_commandline() */
  cfg.do_oneblock = FALSE;                /* reset below after process_commandline() call */

  /* Initializations */
  init_ilogsum();
  FLogsumInit();
  process_commandline(argc, argv, &go, &(cfg.cmfile), &(cfg.sqfile), &(cfg.infmt), &(cfg.outfmt));

  /* Determine if we need to output the alignment all at once in a
   * single block. If not, and we can tell the alignment is going to
   * be big, we'll output temporary alignments of one block of
   * sequences at a time in Pfam format then go back and merge them
   * all at the end. This saves memory by only requiring we keep 1
   * block in memory at a time.
   */
  if((cfg.outfmt != eslMSAFILE_STOCKHOLM && cfg.outfmt != eslMSAFILE_PFAM) || 
     (esl_opt_GetBoolean(go, "--ileaved"))) { 
    /* format is not Stockholm, nor Pfam OR --ileaved enabled for interleaved alignment */
    cfg.do_oneblock = TRUE;
  }
  else { 
    cfg.do_oneblock = FALSE;
  }

  /* update cfg now that we have go */
  cfg.abc_out    = esl_opt_GetBoolean(go, "--dnaout") ? esl_alphabet_Create(eslDNA) : esl_alphabet_Create(eslRNA);

  /* Figure out who we are, and send control there: 
   * we might be an MPI master, an MPI worker, or a serial program.
   */
#ifdef HAVE_MPI

#if eslDEBUGLEVEL >= 1
  pid_t pid;
  /* get the process id */
  pid = getpid();
  printf("#DEBUG: The process id is %d\n", pid);
  fflush(stdout);
#endif

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

      if(cfg.nproc == 1) cm_Fail("MPI mode, but only 1 processor running... (did you execute mpirun?)");

      if (cfg.my_rank > 0)  status = mpi_worker(go, &cfg);
      else 		    status = mpi_master(go, &cfg);

      MPI_Finalize();
    }
  else
#endif /*HAVE_MPI*/
    {
      serial_master(go, &cfg);
    }
  esl_stopwatch_Stop(w);

  /* Close output files */
  if(cfg.ofp   != NULL && esl_opt_IsUsed(go, "-o")) { 
    fclose(cfg.ofp); 
    /* print timing to stdout, only if aln was not output to stdout */
    printf("#\n");
    esl_stopwatch_Display(stdout, w, "# CPU time: ");
  }
  if(cfg.tmpfp != NULL) fclose(cfg.tmpfp); 
  if(cfg.tfp   != NULL) fclose(cfg.tfp); 
  if(cfg.ifp   != NULL) fclose(cfg.ifp); 
  if(cfg.efp   != NULL) fclose(cfg.efp); 
  if(cfg.sfp   != NULL) fclose(cfg.sfp); 
  if(cfg.cmfp  != NULL) cm_file_Close(cfg.cmfp);
  if(cfg.sqfp  != NULL) esl_sqfile_Close(cfg.sqfp);

  if(cfg.abc     != NULL) esl_alphabet_Destroy(cfg.abc);
  if(cfg.abc_out != NULL) esl_alphabet_Destroy(cfg.abc_out);
  esl_stopwatch_Destroy(w);
  esl_getopts_Destroy(go);

  return status;
}
/* serial_master()
 * The serial version of cmalign.
 * 
 * A master can only return if it's successful. All errors are handled immediately and fatally with cm_Fail().
 */
static void
serial_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int             status;                /* Easel status */
  char            errbuf[eslERRBUFSIZE]; /* for printing error messages */
  CM_t           *cm = NULL;             /* a CM */
  int             i, k;                  /* counter over parsetrees and workers */
  int             nali;                  /* index of the (possibly temporary) alignment we are working on */
  int             nseq_cur;              /* number of sequences in current alignment */
  int             nseq_aligned;          /* number of sequences so far aligned */
  int             do_sample;             /* TRUE if we're sampling alignments (--sample) */
  ESL_RANDOMNESS *r = NULL;              /* RNG, used only if --sample */

  /* variables related to output, we may use a tmpfile if seqfile is large */
  int      use_tmpfile;              /* print out current alignment to tmpfile? */
  int      created_tmpfile = FALSE;  /* TRUE if we've created a tmp file for current CM */
  char tmpfile[32] = "esltmpXXXXXX"; /* name of the tmpfile */
  CM_ALNDATA **merged_dataA = NULL;  /* array of all CM_ALNDATA pointers for current alignment */
  int          merged_data_idx = 0;  /* index in merged_dataA */
  int          nmerged;              /* size of merged_dataA */

  /* variables related to reading sequence blocks */
  int            sstatus = eslOK;  /* status from esl_sq_ReadBlock() */
  ESL_SQ_BLOCK  *sq_block;         /* a sequence block */
  ESL_SQ_BLOCK  *nxt_sq_block;     /* sequence block for next loop iteration */
  int            reached_eof;      /* TRUE if we've reached EOF in target sequence file */

  /* variables related to --mapali */
  char         *map_file   = NULL; /* name of alignment file from --mapali */
  CM_ALNDATA  **map_dataA  = NULL; /* array of CM_ALNDATA pointers for mapali alignment */
  int           nmap_data  = 0;    /* number of CM_ALNDATA ptrs in map_dataA */
  int           nmap_cur   = 0;    /* number of CM_ALNDATA ptrs to include in current iteration, 0 unless nali==0 */
  char         *map_sscons = NULL; /* SS_cons from mapali, only used if --mapstr */

  /* variables related to threaded implementation */
  int              ncpus     = 0;    /* number of CPUs working */
  ESL_SQ         **init_sqA  = NULL; /* for initializing workers */
  WORKER_INFO     *info      = NULL; /* the worker info */
  int              infocnt   = 0;    /* number of worker infos */
#ifdef HMMER_THREADS
  ESL_THREADS     *threadObj = NULL;
  ESL_WORK_QUEUE  *queue     = NULL;
#endif
  
  /* General notes on {serial,mpi}_master()'s strategy: 
   * 
   * Ideally, we'd read in all sequences, align them all, and output
   * the alignment. But we're worried about running out of memory, so
   * we read in sequence blocks (sq_block) of at most CMALIGN_MAX_NRES
   * (10 Mb) at a time, and process each in turn. (A parsetree is
   * about 25 bytes per residue, so that should be about 250 Mb). If
   * there's more than one such block, we output each to a tmp file
   * and free the parsetrees afterwards so we don't require too much
   * memory. Once finished, we go through the tmp file and merge all
   * the alignments within it into a single one (without ever storing
   * all of them simultaneously) and output it to the standard output
   * file. In this case we need to use Pfam (1 line/seq) format so we
   * don't have to store the full alignment/set of parsetrees at
   * once. If there's only one block we just output it to the standard
   * output file in interleaved format (which is what previous
   * versions of cmalign did), no tmp file is needed.
   * 
   * Each sequence in the current block is processed independently,
   * i.e. a workunit for a threaded/MPI worker is a single
   * sequence. We could make a workunit a smaller sequence block so
   * workers wouldn't have to update as much, but empirically it seems
   * there's not too much overhead to the updates and sequence
   * alignment times vary significantly so it's advantageous to have a
   * worker process a single sequence at a time, lest they have
   * multiple difficult sequences in a single block. I originally
   * implemented the strategy of splitting big blocks into smaller
   * ones for the threaded implementation, each of which was an
   * independent unit (r3808) but it was significantly more complex
   * with little to no advantage in speed over this implementation.
   */

  if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);
  do_sample  = esl_opt_GetBoolean(go, "--sample") ? TRUE : FALSE;
  if(do_sample) { 
    if((r = esl_randomness_Create(esl_opt_GetInteger(go, "--seed"))) == NULL) cm_Fail("out of memory, trying to create RNG");
  }

#ifdef HMMER_THREADS
  /* initialize thread data */
  ncpus = ESL_MIN(esl_opt_GetInteger(go, "--cpu"), esl_threads_GetCPUCount());
  if (ncpus > 0) {
      threadObj = esl_threads_Create(&pipeline_thread);
      queue = esl_workqueue_Create(ncpus * 2);
  }
#endif

  infocnt = (ncpus == 0) ? 1 : ncpus;
  ESL_ALLOC(info, sizeof(WORKER_INFO) * infocnt);

  /* Read one CM, and make sure there's only one. This fills cfg->abc. */
  status = cm_file_Read(cfg->cmfp, TRUE, &(cfg->abc), &cm);
  if(status != eslOK) cm_Fail(cfg->cmfp->errbuf);
  status = cm_file_Read(cfg->cmfp, TRUE, &(cfg->abc), NULL);
  if(status != eslEOF) cm_Fail("CM file %s does not contain just one CM\n", cfg->cmfp->fname);

  if(cfg->ofp != stdout) output_header(stdout, go, cfg->cmfile, cfg->sqfile, cm, ncpus);

  for (k = 0; k < infocnt; ++k)    {
    info[k].cm          = NULL;
    info[k].dataA       = NULL;
    info[k].n           = 0;
    info[k].mxsize      = esl_opt_GetReal(go, "--mxsize");
    info[k].pass_idx    = esl_opt_GetBoolean(go, "--notrunc") ? PLI_PASS_STD_ANY : PLI_PASS_5P_AND_3P_FORCE;
    info[k].w           = esl_stopwatch_Create();
    info[k].w_tot       = esl_stopwatch_Create();
    info[k].do_failover = (esl_opt_GetBoolean(go, "--hbanded")  && (! esl_opt_GetBoolean(go, "--notrunc"))) ? TRUE : FALSE;
#ifdef HMMER_THREADS
    info[k].queue  = queue;
#endif
  }
  
#ifdef HMMER_THREADS    
  ESL_ALLOC(init_sqA, sizeof(ESL_SQ *) * ESL_MAX(1, ncpus * 2)); // avoid malloc of 0
  for (k = 0; k < ncpus * 2; k++) {
    init_sqA[k] = NULL;
    if((init_sqA[k] = esl_sq_CreateDigital(cfg->abc)) == NULL)          cm_Fail("Failed to allocate a sequence");
    if((status      = esl_workqueue_Init(queue, init_sqA[k])) != eslOK) cm_Fail("Failed to add sequence to work queue");
  }
#endif

  /* initialization */
  nali = nseq_cur = nseq_aligned = 0;
  if((status = initialize_cm(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);

  /* --hmm mode: HMM-only alignment, bypass normal CM alignment pipeline */
  if(esl_opt_GetBoolean(go, "--hmm")) {
    hmm_alignment(go, cfg, cm);
    /* clean up and return; hmm_alignment handles all output */
    for(k = 0; k < infocnt; ++k) { 
      if(info[k].w     != NULL) esl_stopwatch_Destroy(info[k].w);
      if(info[k].w_tot != NULL) esl_stopwatch_Destroy(info[k].w_tot);
    }
    free(info);
#ifdef HMMER_THREADS
    if (ncpus > 0) {
      esl_workqueue_Reset(queue); 
      if(init_sqA != NULL) { 
        for (k = 0; k < ncpus * 2; k++) { 
          if(init_sqA[k] != NULL) esl_sq_Destroy(init_sqA[k]);
        }
        free(init_sqA);
        init_sqA = NULL;
      }
      esl_workqueue_Destroy(queue);
      esl_threads_Destroy(threadObj);
    }
    if(init_sqA != NULL) free(init_sqA);
#endif
    FreeCM(cm);
    return;
  }

  for (k = 0; k < infocnt; ++k) {
    if((status = cm_Clone(cm, errbuf, &(info[k].cm))) != eslOK) cm_Fail(errbuf);
  }
  reached_eof = FALSE;

  /* include the mapali, if nec */
  if((map_file = esl_opt_GetString(go, "--mapali")) != NULL) { 
    if((status = map_alignment(map_file, cm, esl_opt_GetBoolean(go, "--noss"), errbuf, &map_dataA, &nmap_data, &map_sscons)) != eslOK) cm_Fail(errbuf);
    if(esl_opt_GetBoolean(go, "--mapstr") && map_sscons == NULL) cm_Fail("Failed to read SS_cons for --mapstr from %s", map_file);
  }

  /* Our main loop will loop over reading a single large block
   * (<sq_block>) of sequences, up to CM_MAX_RESIDUE_COUNT
   * (100000) residues, and up to CMALIGN_MAX_NSEQ
   * sequences (10,000), but potentially less if we reach the end of
   * the sequence file first.
   */

  /* Read the first block */
  sq_block = esl_sq_CreateDigitalBlock(CMALIGN_MAX_NSEQ, cfg->abc);
  sstatus = esl_sqio_ReadBlock(cfg->sqfp, sq_block, -1, -1, /*max_init_window=*/FALSE, FALSE); /* FALSE says: read complete sequences */
  nxt_sq_block = sq_block; /* special case of first block read */

  while(sstatus == eslOK) { 
#ifdef HMMER_THREADS
    if (ncpus > 0) { 
      for (k = 0; k < infocnt; ++k) esl_threads_AddThread(threadObj, &info[k]);
    }
#endif
    sq_block = nxt_sq_block; /* our current sq_block becomes the one we read on the previous iteration */
    sq_block->first_seqidx = nseq_aligned;
    nseq_cur = sq_block->count;

    /* Before we do any aligning, read the next sequence block, so we
     * can determine if we've reached the end of the seqfile. We need
     * to know this for two reasons:
     *
     * (1) if the first block read above included all sequences (which
     * we won't know until we try to read another block), we don't
     * need to go into memory-saving mode and output to a tmpfile, we
     * can output (in interleaved mode) to the final output file.
     *
     * (2) if do_oneblock (we're trying to output the full alignment
     * as a single block) we need to fail if we still have sequences
     * left, because the sequence file exceeded the size limits. And
     * we want to fail *before* we align all the sequences, so the
     * user isn't cross when the job fails after seemingly going along
     * fine for a while.
     */
    nxt_sq_block = esl_sq_CreateDigitalBlock(CMALIGN_MAX_NSEQ, cfg->abc);
    sstatus = esl_sqio_ReadBlock(cfg->sqfp, nxt_sq_block, -1, -1, /*max_init_window=*/FALSE, FALSE); /* FALSE says: read complete sequences */
    if(sstatus == eslEOF) { 
      reached_eof = TRUE; /* nxt_sq_block will not have been filled */
      esl_sq_DestroyBlock(nxt_sq_block); 
    }
    if((! reached_eof) && cfg->do_oneblock) esl_fatal("Error: the sequence file is too big (has > %d seqs or %d residues) for --ileaved or output\nformat other than Pfam. Use esl-reformat to reformat alignment later.", CMALIGN_MAX_NSEQ, CM_MAX_RESIDUE_COUNT);

    /* align the sequences in the block */
#ifdef HMMER_THREADS
    if (ncpus > 0)  status = thread_loop(info, errbuf, threadObj, queue, sq_block);
    else            status = serial_loop(info, errbuf, sq_block, r);
#else
    status = serial_loop(info, errbuf, sq_block, r);
#endif
    if(status != eslOK) cm_Fail(errbuf);

    /* create a single array of all CM_ALNDATA objects, in original (input) order */
    nmap_cur = (nali == 0) ? nmap_data : 0;
    nmerged  = nseq_cur + nmap_cur;
    ESL_ALLOC(merged_dataA, sizeof(CM_ALNDATA *) * ESL_MAX(1, nmerged)); // avoid malloc of 0
    /* prepend mapali data if nec */
    if(nmap_cur > 0) {
      for(i = 0; i < nmap_cur; i++) merged_dataA[i] = map_dataA[i];
      free(map_dataA); /* don't free the CM_ALNDATA objects, merged_dataA is pointing at them */
      map_dataA = NULL;
    }
    for(k = 0; k < infocnt; ++k) { 
      for(i = 0; i < info[k].n; i++) { 
	merged_data_idx = (info[k].dataA[i]->idx - nseq_aligned) + nmap_cur;
	merged_dataA[merged_data_idx] = info[k].dataA[i];
      }
      /* free dataA pointer from info, but not actual CM_ALNDATA objects, merged_dataA is pointing at them  */
      if(info[k].dataA != NULL) { 
	free(info[k].dataA); 
	info[k].dataA = NULL;
      }
      info[k].n = 0;
    }

    /* output alignment (if do_oneblock we died above if we didn't reach EOF yet) */
    use_tmpfile = (reached_eof && (! created_tmpfile)) ? FALSE : TRUE; /* output to tmpfile only if this is the first alignment and we've aligned all seqs */
    if(use_tmpfile && (! created_tmpfile)) { 
      /* first aln for temporary output file, open the file */	
      if ((status = esl_tmpfile_named(tmpfile, &(cfg->tmpfp))) != eslOK) cm_Fail("Failed to open temporary output file (status %d)", status);
      created_tmpfile = TRUE;
    }
    if((status   = output_alignment(go, cfg, errbuf, cm, (use_tmpfile ? cfg->tmpfp : cfg->ofp), merged_dataA, nseq_cur + nmap_cur, map_sscons)) != eslOK) cm_Fail(errbuf);
    /* optionally output same alignment to regress file */
    if(cfg->rfp != NULL) { 
      if((status = output_alignment(go, cfg, errbuf, cm, cfg->rfp,                              merged_dataA, nseq_cur + nmap_cur, map_sscons)) != eslOK) cm_Fail(errbuf);
    }    
    nali++;
    nseq_aligned += nseq_cur;

    /* output scores to stdout, if -o used */
    if(cfg->ofp != stdout) { 
      if((status =  output_scores(stdout,   cm, errbuf, merged_dataA, nseq_cur + nmap_cur, nmap_cur, cfg->be_verbose)) != eslOK) cm_Fail(errbuf);
    }
    /* output scores to scores file, if --sfp used */
    if(cfg->sfp != NULL) { 
      if(nali == 1) output_header(stdout, go, cfg->cmfile, cfg->sqfile, cm, ncpus);
      if((status =  output_scores(cfg->sfp, cm, errbuf, merged_dataA, nseq_cur + nmap_cur, nmap_cur, cfg->be_verbose)) != eslOK) cm_Fail(errbuf);
    }

    /* free block and worker data */
    esl_sq_DestroyBlock(sq_block);
    sq_block = NULL;
    for(i = 0; i < nmap_cur; i++) { /* free the mapali seqs if nec */
      if(merged_dataA[i]->sq != NULL) esl_sq_Destroy(merged_dataA[i]->sq); 
    }
    for(i = 0; i < nmerged; i++) { 
      cm_alndata_Destroy(merged_dataA[i], FALSE); /* FALSE: don't free sq's, we just free'd them by destroying the block */
    }
    free(merged_dataA);
  } /* end of outer while loop 'while(sstatus == eslOK)' */
  if     (sstatus == eslEFORMAT) cm_Fail("Parse failed (sequence file %s):\n%s\n", cfg->sqfp->filename, esl_sqfile_GetErrorBuf(cfg->sqfp));
  else if(sstatus == eslEMEM)    cm_Fail("Out of memory");
  else if(sstatus != eslEOF)     cm_Fail("Unexpected error while reading sequence file");

  /* if nec, close tmpfile then merge all alignments in it */
  if(created_tmpfile) { 
    fclose(cfg->tmpfp); /* we're done writing to tmpfp */
    cfg->tmpfp = NULL;
    /* merge all temporary alignments now in cfg->tmpfp, and output merged alignment */
    if((status = create_and_output_final_msa(go, cfg, errbuf, cm, nali, tmpfile)) != eslOK) cm_Fail(errbuf);
    remove(tmpfile); 
  }
    
  /* finish insert and el files */
  if(cfg->ifp != NULL) { fprintf(cfg->ifp, "//\n"); }
  if(cfg->efp != NULL) { fprintf(cfg->efp, "//\n"); }

  /* clean up */
#ifdef HMMER_THREADS
  if (ncpus > 0) {
    esl_workqueue_Reset(queue); 
    if(init_sqA != NULL) { 
      for (k = 0; k < ncpus * 2; k++) { 
	if(init_sqA[k] != NULL) esl_sq_Destroy(init_sqA[k]);
      }
      free(init_sqA);
      init_sqA = NULL;
    }
    esl_workqueue_Destroy(queue);
    esl_threads_Destroy(threadObj);
  }
  if(init_sqA != NULL) free(init_sqA); /* init_sqA will be NULL if and only if ncpus == 0 */
#endif
  if(r != NULL) esl_randomness_Destroy(r);
  for(k = 0; k < infocnt; ++k) { 
    if(info[k].cm    != NULL) FreeCM(info[k].cm);
    if(info[k].dataA != NULL) free(info[k].dataA);
    if(info[k].w     != NULL) esl_stopwatch_Destroy(info[k].w);
    if(info[k].w_tot != NULL) esl_stopwatch_Destroy(info[k].w_tot);
  }
  free(info);

  if(map_sscons != NULL) free(map_sscons);
  FreeCM(cm);

  return;
  
  ERROR:
  cm_Fail("Memory allocation error.");
  return;
}

/* output_hmm_insert_info()
 *
 * Emit the per-sequence insert-information file (--ifile) for the
 * --hmm alignment path, derived from the p7 traces <tr[]>.
 *
 * This is the trace-based analog of the CM path's ifile writer
 * (Parsetrees2Alignment() -> insertfp emission in cm_parsetree.c).
 * It is ADDITIVE: the CM path is untouched. The output format and
 * semantics are mirrored from the CM path so VADR's
 * vdr_CmalignParseInsertFile() parses both identically:
 *
 *   model line:  "<cm->name> <cm->clen>\n"   (mirrors output_alignment())
 *   per-seq line: "<name> <L> <spos> <epos>  [<mdlpos> <uapos> <inslen>]...\n"
 *   closing line: "//\n"
 *
 * spos/epos = first/last consensus (match) model position the sequence
 * occupies (match-residue based, like the CM path; deletes don't count).
 * Each insert triplet: <mdlpos> = model position after which the insert
 * occurs (0..clen; 0 = before first consensus, clen = after last),
 * <uapos> = unaligned position (1..L) of the first inserted residue,
 * <inslen> = number of inserted residues.
 *
 * The model-position convention exactly matches p7_tracealign_Seqs()'s
 * map_new_msa() (tracealign.c), which produces the .stk this run wrote:
 *   - a run of emitting p7T_N states  -> mdlpos 0   (5' flanking)
 *   - a run of p7T_I states at node k -> mdlpos k   (insert after match k)
 *   - a run of emitting p7T_C states  -> mdlpos clen (3' flanking)
 * so the ifile is a faithful encoding of the alignment in the .stk.
 *
 * Glocal vs local entry/exit is handled implicitly: spos/epos track the
 * first/last p7T_M regardless of how the model was entered (B->M_k or
 * via leading deletes), and the N/C flanking residues become the
 * mdlpos=0 / mdlpos=clen inserts. VADR runs this under -g (UNIGLOCAL).
 *
 * Only the serial/threaded path is covered (traces land in the shared
 * <tr[]> regardless of --cpu). MPI is out of scope (HAVE_MPI undefined).
 */
static void
output_hmm_insert_info(FILE *ifp, CM_t *cm, P7_HMM *hmm, ESL_SQ **sqarr, P7_TRACE **tr, int nseq)
{
  int idx, z;
  int M = hmm->M;   /* consensus length; == cm->clen for the ML p7 HMM */

  /* model line: byte-identical to the CM path's output_alignment() emission */
  fprintf(ifp, "%s %d\n", cm->name, cm->clen);

  for (idx = 0; idx < nseq; idx++) {
    P7_TRACE *t    = tr[idx];
    int       spos = -1;
    int       epos = -1;

    /* spos/epos: first/last match-state model position (residue-bearing) */
    for (z = 0; z < t->N; z++) {
      if (t->st[z] == p7T_M) {
        if (spos == -1) spos = t->k[z];
        epos = t->k[z];
      }
    }

    fprintf(ifp, "%s %" PRId64 " %d %d", sqarr[idx]->name, sqarr[idx]->n, spos, epos);

    /* Walk the trace once, emitting insert triplets in increasing-mdlpos
     * order (N-term=0, then I_k for k=1..M-1, then C-term=M), which is
     * exactly trace order. Mute (non-emitting) N/C states have i==0 and
     * are skipped; the first emitting state of each flanking run is the
     * one whose predecessor shares the same state type (matching
     * map_new_msa()'s "if (st[z-1]==p7T_N/C)" counting).
     */
    z = 0;
    while (z < t->N) {
      int st = t->st[z];

      if (st == p7T_N && z > 0 && t->st[z-1] == p7T_N) {
        /* 5' flanking run -> mdlpos 0 */
        int first = t->i[z];
        int len   = 0;
        while (z < t->N && t->st[z] == p7T_N) { if (t->i[z] > 0) len++; z++; }
        if (len > 0) fprintf(ifp, "  %d %d %d", 0, first, len);
      }
      else if (st == p7T_I) {
        /* insert after match position k -> mdlpos k */
        int k     = t->k[z];
        int first = t->i[z];
        int len   = 0;
        while (z < t->N && t->st[z] == p7T_I && t->k[z] == k) { len++; z++; }
        fprintf(ifp, "  %d %d %d", k, first, len);
      }
      else if (st == p7T_C && z > 0 && t->st[z-1] == p7T_C) {
        /* 3' flanking run -> mdlpos M (== clen) */
        int first = t->i[z];
        int len   = 0;
        while (z < t->N && t->st[z] == p7T_C) { if (t->i[z] > 0) len++; z++; }
        if (len > 0) fprintf(ifp, "  %d %d %d", M, first, len);
      }
      else z++;
    }

    fprintf(ifp, "\n");
  }

  /* closing line, mirrors the CM path (cmalign.c serial_master end) */
  fprintf(ifp, "//\n");
}

/* hmm_alignment()
 * 
 * HMM-only alignment mode (--hmm). Bypasses the CM alignment pipeline
 * entirely. Uses the CM's embedded p7 HMM to align sequences, producing
 * Stockholm output via p7_tracealign_Seqs().
 *
 * Three sub-modes:
 *   --hmm --hmmvit:    Viterbi traces (fastest, least accurate)
 *   --hmm --hmmnoband: Full (unbanded) optimal accuracy alignment
 *   --hmm (default):   Viterbi-banded optimal accuracy alignment
 */
static void
hmm_alignment(ESL_GETOPTS *go, struct cfg_s *cfg, CM_t *cm)
{
  int           status;
  P7_HMM       *hmm     = NULL;   /* the p7 HMM from the CM */
  P7_BG        *bg      = NULL;   /* null model */
  P7_PROFILE   *gm      = NULL;   /* generic profile */
  P7_GMX       *gx      = NULL;   /* generic DP matrix (Viterbi) */
  P7_GMX       *gxf     = NULL;   /* Forward matrix (unbanded OA) */
  P7_GMX       *gxb     = NULL;   /* Backward matrix (unbanded OA) */
  ESL_SQ      **sqarr   = NULL;   /* array of sequences */
  P7_TRACE    **tr      = NULL;   /* array of traces */
  ESL_MSA      *msa     = NULL;   /* output alignment */
  int           nseq    = 0;      /* number of sequences */
  int           nalloc  = 0;      /* allocated size of sqarr/tr */
  int           idx;
  int           p7mode;           /* p7 profile mode */
  float         sc, fwdsc, oasc;
  char          errbuf[eslERRBUFSIZE];

  int do_hmmvit    = (cm->align_opts & CM_ALIGN_P7HMMVIT)    ? TRUE : FALSE;
  int do_hmmnoband = (cm->align_opts & CM_ALIGN_P7HMMNOBAND) ? TRUE : FALSE;
  int do_bandedoa  = (! do_hmmvit && ! do_hmmnoband)         ? TRUE : FALSE;
  /* --p7ibv (w/--hmm): derive banded-OA bands via the D&C IBV deriver instead
   * of a full p7_GViterbi + trace, skipping the O(M*L) P7_GMX allocation.
   * CLI validation guarantees --p7ibv only reaches here in banded-OA mode.
   */
  int do_p7ibv     = (cm->p7_use_ibv && do_bandedoa)         ? TRUE : FALSE;
  /* brief 26_0430-182: mirror p7_ibv.c:1793's expression exactly -- the same test
   * used to calibrate p7bpad/node-pad at align-time. Drives both the Tgm
   * profile config (Part A) and the IBV deriver's do_trunc arg (Part B).
   */
  int do_trunc     = (cm->align_opts & CM_ALIGN_TRUNC)       ? TRUE : FALSE;

  /* banded functions declared in cm_p7_band.c */
  extern int p7_kbands2gbands(int *i2k, int *kmin, int *kmax, int L, int M, P7_GBANDS **ret_bnd);
  extern int my_p7_GForwardBanded(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMXB *gxb, float *opt_sc);
  extern int p7_GBackwardBanded(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMXB *gxb, float *opt_sc);
  extern int p7_GDecodingBanded(const P7_PROFILE *gm, const P7_GMXB *fwd, P7_GMXB *bck, P7_GMXB *pp, float overall_sc);
  extern int p7_GOptimalAccuracyBanded(const P7_PROFILE *gm, const P7_GMXB *pp, P7_GMXB *gx, float *ret_e);
  extern int p7_GOATraceBanded(const P7_PROFILE *gm, const P7_GMXB *pp, const P7_GMXB *gx, P7_TRACE *tr);
  extern int p7_GCheckptFBDecode_Banded(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMXB *pp, float *ret_fwdsc); /* brief 26_0526-016 */
  extern int p7_GCheckptOA_Banded(const P7_PROFILE *gm, P7_GMXB *pp, P7_TRACE *tr, float *ret_oasc);                    /* brief 26_0526-016 */
  extern P7_GMXB *p7b_pp_Create(P7_GBANDS *bnd);                                                                       /* brief 26_0526-017: compact 2-cell resident pp */

  /* Verify the CM has a valid p7 HMM */
  if (! (cm->flags & CMH_MLP7)) cm_Fail("--hmm requires a CM file with an embedded p7 HMM (use cmconvert)");
  hmm = cm->mlp7;

  /* Set up null model and profile mode (shared across serial/threaded paths) */
  bg = p7_bg_Create(hmm->abc);
  p7mode = esl_opt_GetBoolean(go, "-g") ? p7_UNIGLOCAL : p7_UNILOCAL;

  /* Read all sequences into an array */
  nalloc = 256;
  ESL_ALLOC(sqarr, sizeof(ESL_SQ *) * nalloc);
  nseq = 0;
  while (1) {
    ESL_SQ *sq = esl_sq_CreateDigital(cfg->abc);
    status = esl_sqio_Read(cfg->sqfp, sq);
    if (status == eslEOF) { esl_sq_Destroy(sq); break; }
    if (status != eslOK)  cm_Fail("Error reading sequence file %s: %s", cfg->sqfile, esl_sqfile_GetErrorBuf(cfg->sqfp));
    if (nseq >= nalloc) {
      nalloc *= 2;
      ESL_REALLOC(sqarr, sizeof(ESL_SQ *) * nalloc);
    }
    sqarr[nseq++] = sq;
  }
  if (nseq == 0) cm_Fail("No sequences found in %s", cfg->sqfile);

  /* Allocate trace array */
  ESL_ALLOC(tr, sizeof(P7_TRACE *) * nseq);
  for (idx = 0; idx < nseq; idx++)
    tr[idx] = do_hmmvit ? p7_trace_Create() : p7_trace_CreateWithPP();

  /* ---- Compute traces for each sequence (threaded or serial) ---- */
  {
    int ncpus = 0;
#ifdef HMMER_THREADS
    ncpus = ESL_MIN(esl_opt_GetInteger(go, "--cpu"), esl_threads_GetCPUCount());
#endif

    if (ncpus > 0) {
#ifdef HMMER_THREADS
      /* Threaded path: one sequence per worker, following cmalign's thread_loop pattern */
      ESL_THREADS    *threadObj = NULL;
      ESL_WORK_QUEUE *queue     = NULL;
      WORKER_INFO    *winfo     = NULL;
      ESL_SQ        **init_sqA  = NULL;
      int             k;

      threadObj = esl_threads_Create(&hmm_pipeline_thread);
      queue     = esl_workqueue_Create(ncpus * 2);

      ESL_ALLOC(winfo,    sizeof(WORKER_INFO) * ncpus);
      ESL_ALLOC(init_sqA, sizeof(ESL_SQ *)    * ncpus * 2);

      /* Initialize work queue with token sequences */
      for (k = 0; k < ncpus * 2; k++) {
	init_sqA[k] = esl_sq_CreateDigital(cfg->abc);
	esl_workqueue_Init(queue, init_sqA[k]);
      }

      /* Initialize per-thread WORKER_INFO with HMM-specific fields */
      for (k = 0; k < ncpus; k++) {
	winfo[k].queue       = queue;
	winfo[k].bg          = p7_bg_Create(hmm->abc);
	winfo[k].gm          = p7_profile_Create(hmm->M, hmm->abc);
	/* brief 26_0430-182 Part A: Tgm (5'+3' truncation-aware local profile) when
	 * do_trunc, mirroring cm_alndata.c:459-461's proven --p7band pattern.
	 * Per-sequence length is set later by p7_ReconfigLength() (non-trunc)
	 * or the Tgm-aware re-setup in hmm_pipeline_thread() (trunc). */
	if (do_trunc) {
	  p7_ProfileConfig(hmm, bg, winfo[k].gm, 400, p7_LOCAL);
	  p7_ProfileConfig5PrimeAnd3PrimeTrunc(winfo[k].gm, 400);
	} else {
	  p7_ProfileConfig(hmm, bg, winfo[k].gm, 400, p7mode);
	}
	winfo[k].hmm         = hmm;
	/* Under --p7ibv the banded-OA path needs no full Viterbi P7_GMX. */
	winfo[k].gx          = (do_hmmvit || (do_bandedoa && ! do_p7ibv)) ? p7_gmx_Create(hmm->M, 400) : NULL;
	winfo[k].gxf         = do_hmmnoband ? p7_gmx_Create(hmm->M, 400) : NULL;
	winfo[k].gxb         = do_hmmnoband ? p7_gmx_Create(hmm->M, 400) : NULL;
	winfo[k].hmm_tr      = tr;  /* shared trace array; worker writes to tr[seqidx] */
	winfo[k].do_hmmvit   = do_hmmvit;
	winfo[k].do_hmmnoband = do_hmmnoband;
	winfo[k].do_p7ibv    = do_p7ibv;
	winfo[k].do_trunc    = do_trunc;
	winfo[k].ibv_delta   = esl_opt_GetInteger(go, "--p7ibv-delta");
	/* brief 26_0526-017 Part A: default base_slab to the knee (memory-only; byte-invariant
	 * per gate A1); honor an explicit --p7ibv-base-slab unchanged. */
	winfo[k].ibv_base_slab = esl_opt_IsDefault(go, "--p7ibv-base-slab")
	                         ? HMM_P7IBV_KNEE_BASE_SLAB
	                         : esl_opt_GetInteger(go, "--p7ibv-base-slab");
	/* CM only needed by the IBV/kmeranchor/kmerchain derivers (for cm->fp7);
	 * else unused in --hmm mode. brief 26_0628-032: kmeranchor/kmerchain also need it. */
	winfo[k].cm          = (do_p7ibv || cm->p7_use_kmeranchor || cm->p7_use_kmerchain) ? cm : NULL;
	winfo[k].dataA       = NULL;
	winfo[k].n           = 0;
	winfo[k].mxsize      = esl_opt_GetReal(go, "--mxsize");
	winfo[k].pass_idx    = 0;
	winfo[k].w           = NULL;
	winfo[k].w_tot       = NULL;
	winfo[k].do_failover = FALSE;

	esl_threads_AddThread(threadObj, &winfo[k]);
      }

      /* Distribute sequences to workers (one seq per work unit) */
      hmm_thread_loop(&winfo[0], threadObj, queue, sqarr, nseq);

      /* Clean up threaded resources */
      esl_workqueue_Reset(queue);
      for (k = 0; k < ncpus * 2; k++) esl_sq_Destroy(init_sqA[k]);
      free(init_sqA);
      esl_workqueue_Destroy(queue);
      esl_threads_Destroy(threadObj);
      for (k = 0; k < ncpus; k++) {
	p7_profile_Destroy(winfo[k].gm);
	p7_bg_Destroy(winfo[k].bg);
	if (winfo[k].gx)  p7_gmx_Destroy(winfo[k].gx);
	if (winfo[k].gxf) p7_gmx_Destroy(winfo[k].gxf);
	if (winfo[k].gxb) p7_gmx_Destroy(winfo[k].gxb);
      }
      free(winfo);
#endif /* HMMER_THREADS */
    }
    else {
      /* Serial path: single profile and matrices */
      gm = p7_profile_Create(hmm->M, hmm->abc);
      /* brief 26_0430-182 Part A: Tgm when do_trunc (see winfo[k].gm comment above). */
      if (do_trunc) {
	p7_ProfileConfig(hmm, bg, gm, 400, p7_LOCAL);
	p7_ProfileConfig5PrimeAnd3PrimeTrunc(gm, 400);
      } else {
	p7_ProfileConfig(hmm, bg, gm, 400, p7mode);
      }

      if (do_hmmvit || (do_bandedoa && ! do_p7ibv)) gx  = p7_gmx_Create(hmm->M, 400);
      if (do_hmmnoband)           { gxf = p7_gmx_Create(hmm->M, 400); gxb = p7_gmx_Create(hmm->M, 400); }

      for (idx = 0; idx < nseq; idx++) {
	ESL_SQ *sq = sqarr[idx];

	/* brief 26_0430-182 Part A: p7_ReconfigLength() unconditionally overwrites
	 * xsc[N/C/J][MOVE|LOOP], which would clobber the Tgm -eslINFINITY
	 * N->N/C->C loop-disable set by p7_ProfileConfig5PrimeAnd3PrimeTrunc().
	 * Re-run the Tgm setup per-sequence (sq->n) instead when do_trunc. */
	if (do_trunc) {
	  p7_ProfileConfig(hmm, bg, gm, sq->n, p7_LOCAL);
	  p7_ProfileConfig5PrimeAnd3PrimeTrunc(gm, sq->n);
	} else {
	  p7_ReconfigLength(gm, sq->n);
	}

	/* preflight: check HMM matrix size vs --mxsize before GrowTo.
	 * Skipped under --p7ibv: the full P7_GMX is never allocated.
	 * brief 26_0628-032: also skipped under kmeranchor/kmerchain -- like --p7ibv,
	 * the full P7_GMX is only touched on the rare ncells==0 fallback, not
	 * on the genome-scale success path this brief's memory story depends on. */
	if (! do_p7ibv && ! cm->p7_use_kmeranchor && ! cm->p7_use_kmerchain) {
	  double single_bytes = (double) sizeof(float) * (double)(hmm->M + 1) * (double)(sq->n + 1) * (double) p7G_NSCELLS;
	  int    nmat         = do_hmmnoband ? 2 : 1;
	  double needed_mb    = (single_bytes * (double) nmat) / (1024.0 * 1024.0);
	  double mxsize_limit = esl_opt_GetReal(go, "--mxsize");
	  if (needed_mb > mxsize_limit) {
	    int recommended_mxsize = (int)(ceil(needed_mb / 1024.0) * 1024.0);
	    cm_Fail("HMM-only alignment mx needs %.2f Mb > %.2f Mb limit. Use --mxsize %d.",
		    needed_mb, mxsize_limit, recommended_mxsize);
	  }
	}

	if (do_hmmvit) {
	  /* --- Mode 1: Viterbi trace --- */
	  p7_gmx_GrowTo(gx, hmm->M, sq->n);
	  p7_GViterbi(sq->dsq, sq->n, gm, gx, &sc);
	  p7_trace_Reuse(tr[idx]);
	  if ((status = p7_GTrace(sq->dsq, sq->n, gm, gx, tr[idx])) != eslOK)
	    cm_Fail("p7_GTrace() failed for sequence %s", sq->name);
	}
	else if (do_hmmnoband) {
	  /* --- Mode 2: Unbanded optimal accuracy --- */
	  p7_gmx_GrowTo(gxf, hmm->M, sq->n);
	  p7_gmx_GrowTo(gxb, hmm->M, sq->n);

	  p7_GForward (sq->dsq, sq->n, gm, gxf, &fwdsc);
	  p7_GBackward(sq->dsq, sq->n, gm, gxb, NULL);
	  p7_GDecoding(gm, gxf, gxb, gxb);
	  p7_GOptimalAccuracy(gm, gxb, gxf, &oasc);
	  p7_trace_Reuse(tr[idx]);
	  p7_GOATrace(gm, gxb, gxf, tr[idx]);
	}
	else {
	  /* --- Mode 3 (default): Viterbi-banded optimal accuracy --- */
	  int     *i2k   = NULL;
	  int     *kmin  = NULL;
	  int     *kmax  = NULL;
	  int      ncells = 0;
	  int      pad   = 30;
	  P7_GBANDS *bnd = NULL;
	  P7_GMXB *bxf   = NULL;
	  P7_GMXB *bxb   = NULL;
	  P7_TRACE *vtr  = NULL;
	  float    bwdsc = 0.;                                                    /* brief 26_0430-135b: capture backward total */
	  int      p7ibv_delta = esl_opt_GetInteger(go, "--p7ibv-delta");        /* brief 26_0430-135b */
	  int      do_widen = (getenv("P135B_FORCE_WIDEN") != NULL) ? TRUE : FALSE; /* brief 26_0430-135b widen override */

	  if (do_p7ibv) {
	    /* IBV D&C deriver: bands straight from cm->fp7, no full P7_GMX. */
	    if (cm->fp7 == NULL || cm->fp7->M != hmm->M)
	      cm_Fail("--hmm --p7ibv requires cm->fp7 with M matching the ML p7 HMM");
	    /* brief 26_0526-017 Part A: default base_slab to the knee (memory-only; byte-invariant
	     * per gate A1); honor an explicit --p7ibv-base-slab unchanged. */
	    if ((status = p7_Seq2BandsIBV_dnc(cm, errbuf, sq->dsq, sq->n,
					      p7ibv_delta,
					      (esl_opt_IsDefault(go, "--p7ibv-base-slab")
					       ? HMM_P7IBV_KNEE_BASE_SLAB
					       : esl_opt_GetInteger(go, "--p7ibv-base-slab")),
					      do_widen, /* brief 26_0430-135b: P135B_FORCE_WIDEN override; default FALSE (non-truncated --hmm) */
					      FALSE,    /* brief 26_0430-172: do_kband (unbanded D&C on --hmm path) */
					      do_trunc, /* brief 26_0430-182 Part B: was hardcoded FALSE; --hmm defaults to truncated (CM_ALIGN_TRUNC set unless --notrunc), so this must track it like p7_ibv.c:1793 */
					      cm->p7_ibv_mode, cm->p7_ibv_width, /* brief 26_0430-140 */
					      &i2k, &kmin, &kmax, &ncells)) != eslOK)
	      cm_Fail("p7_Seq2BandsIBV_dnc() failed for sequence %s: %s", sq->name, errbuf);
	  }
	  else {
	    /* brief 26_0628-032: k-mer anchor/chain deriver, opt-in via --p7kmeranchor/
	     * --p7kmerchain (mirrors cm_alndata.c:536-566's --p7band dispatch).
	     * brief 26_0628-038: do_trunc now threaded through, mirroring cm_alndata.c's
	     * CM-mode dispatch (cm->align_opts & CM_ALIGN_TRUNC). */
	    int did_kmer = FALSE;
	    int *local_nodepad = NULL;
	    if ((cm->p7_use_kmeranchor || cm->p7_use_kmerchain) && (cm->flags & CMH_P7NODEPAD)) {
	      int k;
	      ESL_ALLOC(local_nodepad, sizeof(int) * (hmm->M + 1));
	      for (k = 0; k <= hmm->M; k++) local_nodepad[k] = cm->p7_cm_nodepad[k] + cm->p7bpad;
	    }
	    if (cm->p7_use_kmeranchor) {
	      did_kmer = TRUE;
	      if ((status = p7_Seq2BandsKmerAnchor(cm, errbuf, sq->dsq, sq->n, local_nodepad,
						   do_trunc, /* brief 26_0628-038: track CM_ALIGN_TRUNC like cm_alndata.c:541 */
						   &i2k, &kmin, &kmax, &ncells)) != eslOK)
		cm_Fail("p7_Seq2BandsKmerAnchor() failed for sequence %s: %s", sq->name, errbuf);
	    }
	    else if (cm->p7_use_kmerchain) {
	      did_kmer = TRUE;
	      if ((status = p7_Seq2BandsKmerChain(cm, errbuf, sq->dsq, sq->n, local_nodepad,
						  do_trunc, /* brief 26_0628-038: track CM_ALIGN_TRUNC like cm_alndata.c:558 */
						  &i2k, &kmin, &kmax, &ncells)) != eslOK)
		cm_Fail("p7_Seq2BandsKmerChain() failed for sequence %s: %s", sq->name, errbuf);
	    }
	    if (local_nodepad) free(local_nodepad);

	    if (did_kmer && ncells == 0 && ! cm->p7_kmerchain_fallback_vit) {
	      /* brief 26_0628-047: M-gate/N-gate fired, or no anchor found -- default
	       * fallback target is now --p7ibv's D&C deriver instead of a
	       * Vit-trace band (mir-2807: the old Vit-trace fallback itself
	       * landed on the wrong alignment even when the gate correctly
	       * fired; --p7ibv is known more accurate). brief 26_0628-050: use
	       * cm->p7_ibv_delta (struct default 3000, same value cm_alndata.c's
	       * --p7band fallback already uses) rather than --p7ibv-delta's own
	       * CLI default (20000) -- 20000 is provably worse than <=10000 on
	       * mir-2807 and slightly worse on SNORA16 (a different, higher-
	       * scoring-but-wrong HMM registration only becomes reachable at wide
	       * deltas). --p7ibv-delta can never be set here anyway (mutually
	       * exclusive with --p7kmeranchor/--p7kmerchain), so this only
	       * changes this fallback's own behavior. */
	      int p7ibv_base_slab = (esl_opt_IsDefault(go, "--p7ibv-base-slab")
				      ? HMM_P7IBV_KNEE_BASE_SLAB
				      : esl_opt_GetInteger(go, "--p7ibv-base-slab"));
	      if ((status = kmer_gate_p7ibv_fallback(cm, errbuf, sq->dsq, sq->n, do_trunc,
						      cm->p7_ibv_delta, p7ibv_base_slab,
						      &i2k, &kmin, &kmax, &ncells)) != eslOK)
		cm_Fail("kmer_gate_p7ibv_fallback() failed for sequence %s: %s", sq->name, errbuf);
	    }
	    if (! did_kmer || ncells == 0) {
	      /* Default Mode-3 path (no kmer flag set); the kmeranchor/kmerchain
	       * ncells==0 fallback when --p7kmerchain-fbvit reverts to
	       * the old behavior; and the safety net when the --p7ibv fallback
	       * above itself also found nothing usable (mirrors cm_alndata.c's
	       * kmeranchor->vitband / kmerchain->vitband fallback shape). */
	      vtr = p7_trace_Create();
	      p7_gmx_GrowTo(gx, hmm->M, sq->n);
	      p7_GViterbi(sq->dsq, sq->n, gm, gx, &sc);
	      p7_trace_Reuse(vtr);
	      status = p7_GTrace(sq->dsq, sq->n, gm, gx, vtr);

	      if (status != eslOK || vtr->N == 0) {
		P7_GMX *fallback_gxf = p7_gmx_Create(hmm->M, sq->n);
		P7_GMX *fallback_gxb = p7_gmx_Create(hmm->M, sq->n);
		p7_GForward (sq->dsq, sq->n, gm, fallback_gxf, &fwdsc);
		p7_GBackward(sq->dsq, sq->n, gm, fallback_gxb, NULL);
		p7_GDecoding(gm, fallback_gxf, fallback_gxb, fallback_gxb);
		p7_GOptimalAccuracy(gm, fallback_gxb, fallback_gxf, &oasc);
		p7_trace_Reuse(tr[idx]);
		p7_GOATrace(gm, fallback_gxb, fallback_gxf, tr[idx]);
		p7_gmx_Destroy(fallback_gxf);
		p7_gmx_Destroy(fallback_gxb);
		p7_trace_Destroy(vtr);
		continue;
	      }

	      {
		int tpos;
		ESL_ALLOC(i2k, sizeof(int) * (sq->n + 1));
		esl_vec_ISet(i2k, (sq->n + 1), -1);
		for (tpos = 0; tpos < vtr->N; tpos++) {
		  if (vtr->st[tpos] == p7T_M) {
		    int i = vtr->i[tpos];
		    int k = vtr->k[tpos];
		    if (i >= 1 && i <= sq->n && k >= 1 && k <= hmm->M)
		      i2k[i] = k;
		  }
		}
	      }

	      if ((status = p7_pins2bands(i2k, errbuf, sq->n, hmm->M, pad, &kmin, &kmax, &ncells)) != eslOK)
		cm_Fail("p7_pins2bands() failed for sequence %s: %s", sq->name, errbuf);
	    }
	  }
	  if ((status = p7_kbands2gbands(i2k, kmin, kmax, sq->n, hmm->M, &bnd)) != eslOK)
	    cm_Fail("p7_kbands2gbands() failed for sequence %s", sq->name);
	  if (getenv("BRIEF035_MEMPOINT") != NULL)
	    fprintf(stderr, "#MEMPOINT after_gbands seq=%s L=%d rss_kb=%ld\n", sq->name, (int) sq->n, brief035_rss_kb());

	  if (getenv("INFERNAL_HMM_CKPT_OFF") == NULL) {
	    /* brief 26_0526-016: sqrt(nrow)-checkpointed F/B/Decode/OA/traceback.
	     * Default-on (byte-exact vs the full path at norovirus/dengue/
	     * sars/HSV); set INFERNAL_HMM_CKPT_OFF to force the full path.
	     * bxb holds the resident posterior; no full F, B, or OA matrix
	     * is ever materialized (bxf is not allocated). brief 26_0526-017: bxb uses
	     * the compact 2-cell (M,I) pp allocator, ~1/3 smaller than 3-cell. */
	    bxb = p7b_pp_Create(bnd);
	    if (getenv("BRIEF035_MEMPOINT") != NULL)
	      fprintf(stderr, "#MEMPOINT after_cp9alloc_ckpt seq=%s L=%d rss_kb=%ld\n", sq->name, (int) sq->n, brief035_rss_kb());
	    if ((status = p7_GCheckptFBDecode_Banded(sq->dsq, sq->n, gm, bxb, &fwdsc)) != eslOK)
	      cm_Fail("p7_GCheckptFBDecode_Banded() failed for sequence %s", sq->name);
	    p7_trace_Reuse(tr[idx]);
	    if ((status = p7_GCheckptOA_Banded(gm, bxb, tr[idx], &oasc)) != eslOK)
	      cm_Fail("p7_GCheckptOA_Banded() failed for sequence %s", sq->name);
	  }
	  else {
	  bxf = p7_gmxb_Create(bnd);
	  bxb = p7_gmxb_Create(bnd);
	  if (getenv("BRIEF035_MEMPOINT") != NULL)
	    fprintf(stderr, "#MEMPOINT after_cp9alloc_nockpt seq=%s L=%d rss_kb=%ld\n", sq->name, (int) sq->n, brief035_rss_kb());

	  if ((status = my_p7_GForwardBanded(sq->dsq, sq->n, gm, bxf, &fwdsc)) != eslOK)
	    cm_Fail("my_p7_GForwardBanded() failed for sequence %s", sq->name);
	  if ((status = p7_GBackwardBanded(sq->dsq, sq->n, gm, bxb, &bwdsc)) != eslOK)
	    cm_Fail("p7_GBackwardBanded() failed for sequence %s", sq->name);
	  if (getenv("P135B_FB_INSTRUMENT") != NULL)
	    fprintf(stderr, "#P135B_FBTOTAL seq=%s M=%d L=%d delta=%d widen=%d ncells=%d fwd=%.6f bwd=%.6f gap=%.6f\n",
		    sq->name, hmm->M, (int) sq->n, p7ibv_delta, do_widen, ncells, fwdsc, bwdsc, fwdsc - bwdsc);
	  if ((status = p7_GDecodingBanded(gm, bxf, bxb, bxb, fwdsc)) != eslOK)
	    cm_Fail("p7_GDecodingBanded() failed for sequence %s", sq->name);
	  if ((status = p7_GOptimalAccuracyBanded(gm, bxb, bxf, &oasc)) != eslOK)
	    cm_Fail("p7_GOptimalAccuracyBanded() failed for sequence %s", sq->name);

	  p7_trace_Reuse(tr[idx]);
	  if ((status = p7_GOATraceBanded(gm, bxb, bxf, tr[idx])) != eslOK)
	    cm_Fail("p7_GOATraceBanded() failed for sequence %s", sq->name);
	  }

	  if (getenv("BRIEF035_MEMPOINT") != NULL)
	    fprintf(stderr, "#MEMPOINT alignment_peak seq=%s L=%d rss_kb=%ld\n", sq->name, (int) sq->n, brief035_rss_kb());
	  free(i2k);
	  free(kmin);
	  free(kmax);
	  if (vtr) p7_trace_Destroy(vtr);
	  p7_gbands_Destroy(bnd);
	  if (bxf) p7_gmxb_Destroy(bxf);
	  p7_gmxb_Destroy(bxb);
	}
      } /* end serial for loop */

      if (gx  != NULL) p7_gmx_Destroy(gx);  gx  = NULL;
      if (gxf != NULL) p7_gmx_Destroy(gxf); gxf = NULL;
      if (gxb != NULL) p7_gmx_Destroy(gxb); gxb = NULL;
      p7_profile_Destroy(gm); gm = NULL;
    } /* end serial/threaded branch */
  }

  /* ---- Convert traces to MSA ---- */
  if ((status = p7_tracealign_Seqs(sqarr, tr, nseq, hmm->M, p7_ALL_CONSENSUS_COLS, hmm, &msa)) != eslOK)
    cm_Fail("p7_tracealign_Seqs() failed");

  /* Add SS_cons from CM to the MSA if available.
   * The p7 RF annotation marks consensus columns (M states), which
   * correspond 1:1 to CM consensus positions. We map CM SS_cons
   * onto the MSA's consensus columns.
   */
  if (cm->cmcons != NULL && cm->cmcons->cstr != NULL && msa->rf != NULL) {
    int cpos, apos;
    ESL_ALLOC(msa->ss_cons, sizeof(char) * (msa->alen + 1));
    cpos = 0;
    for (apos = 0; apos < msa->alen; apos++) {
      if (msa->rf[apos] != '.' && msa->rf[apos] != '~') {
        /* consensus column */
        msa->ss_cons[apos] = (cpos < cm->clen) ? cm->cmcons->cstr[cpos] : '.';
        cpos++;
      } else {
        msa->ss_cons[apos] = '.';
      }
    }
    msa->ss_cons[msa->alen] = '\0';
  }

  /* Convert to DNA if --dnaout */
  if (esl_opt_GetBoolean(go, "--dnaout") && cfg->abc_out->type == eslDNA) {
    int i2, apos;
    for (i2 = 0; i2 < msa->nseq; i2++) {
      for (apos = 0; apos < msa->alen; apos++) {
        if (msa->aseq[i2][apos] == 'U') msa->aseq[i2][apos] = 'T';
        if (msa->aseq[i2][apos] == 'u') msa->aseq[i2][apos] = 't';
      }
    }
  }

  /* Write the MSA */
  status = esl_msafile_Write(cfg->ofp, msa, cfg->outfmt);
  if (status != eslOK) cm_Fail("Failed to write alignment");

  /* Emit per-sequence insert info (--ifile), derived from the p7 traces.
   * Gated on cfg->ifp (open iff --ifile was given), mirroring the CM path.
   * The CM-path ifile writer is left untouched; this is a parallel
   * trace-based emitter (see output_hmm_insert_info() above). */
  if (cfg->ifp != NULL) output_hmm_insert_info(cfg->ifp, cm, hmm, sqarr, tr, nseq);

  /* Clean up */
  esl_msa_Destroy(msa);
  for (idx = 0; idx < nseq; idx++) {
    p7_trace_Destroy(tr[idx]);
    esl_sq_Destroy(sqarr[idx]);
  }
  free(tr);
  free(sqarr);
  if (gx  != NULL) p7_gmx_Destroy(gx);
  if (gxf != NULL) p7_gmx_Destroy(gxf);
  if (gxb != NULL) p7_gmx_Destroy(gxb);
  if (gm  != NULL) p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);

  return;

 ERROR:
  cm_Fail("Memory allocation error in hmm_alignment()");
  return;
}

/* serial_loop(): 
 * 
 * Align all sequences in a sequence block and store parsetrees.
 * 
 * serial_loop() unlike thread_loop() gets a ESL_RANDOMNESS <r> passed
 * in, it is required for sampling alignments with --sample.
 * serial_loop() will always be called if --sample is used, because we
 * enforce that if HMMER_THREADS is defined --cpu 0 must accompany
 * --sample. The reason for this is otherwise the sampled alignments
 * would be affected by the number of threads, since each thread
 * requires its own (separately-seeded) RNG.
 */
static int
serial_loop(WORKER_INFO *info, char *errbuf, ESL_SQ_BLOCK *sq_block, ESL_RANDOMNESS *r)
{
  int status;
  int i;  /* counter over sequences */
  ESL_SQ  *sqp = NULL; /* ptr to a ESL_SQ, only used if there's an error */
  CM_P7_OM_HOLDER om_holder; /* reusable --p7pinbridge LOCAL profile/OPROFILE (brief 26_0430-090) */

  /* allocate dataA */
  info->n = sq_block->count;
  ESL_ALLOC(info->dataA, sizeof(CM_ALNDATA *) * info->n);
  for(i = 0; i < info->n; i++) info->dataA[i] = NULL;

  cm_p7_om_holder_Init(&om_holder);
  for(i = 0; i < info->n; i++) {
    status = DispatchSqAlignment(info->cm, errbuf, sq_block->list + i, sq_block->first_seqidx + i, info->mxsize,
				 TRMODE_UNKNOWN, info->pass_idx, FALSE, /* FALSE: info->cm->cp9b not valid */
				 info->w, info->w_tot, r, &om_holder, &(info->dataA[i]));
    /* If alignment failed: potentially retry alignment in HMM banded
     * std (non-truncated) mode. We will only possibly do this if our
     * initial try was HMM banded truncated alignment (if not,
     * info->do_failover will be FALSE).
     */
    if(status == eslEAMBIGUOUS && info->do_failover == TRUE) { 
      assert(info->cm->align_opts & CM_ALIGN_TRUNC);
      info->cm->align_opts &= ~CM_ALIGN_TRUNC; /* lower truncated alignment flag, just for this sequence */
      status = DispatchSqAlignment(info->cm, errbuf, sq_block->list + i, sq_block->first_seqidx + i, info->mxsize,
				   TRMODE_UNKNOWN, PLI_PASS_STD_ANY, FALSE, /* USE PLI_PASS_STD_ANY; FALSE: info->cm->cp9b not valid */
				   info->w, info->w_tot, r, &om_holder, &(info->dataA[i]));
      info->cm->align_opts |= CM_ALIGN_TRUNC; /* reraise truncated alignment flag */
    }
    if(status != eslOK) {
      sqp = (sq_block->list + i);
      fprintf(stderr, "Problem during alignment of sequence %s\n", sqp->name);
      cm_Fail(errbuf);
    }
  }
  cm_p7_om_holder_Reset(&om_holder);
  return eslOK;
  
 ERROR: 
  ESL_FAIL(status, errbuf, "out of memory");
  return status; /* NEVERREACHED */
}
 
#ifdef HMMER_THREADS
static int
thread_loop(WORKER_INFO *info, char *errbuf, ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQ_BLOCK *sq_block)
{
  int      status = eslOK;
  int      i, k;           /* counter over sequences, workers */
  ESL_SQ  *sq;
  void    *new_sq;
  ESL_SQ  *empty_sq;
  int      nworkers = esl_threads_GetWorkerCount(obj);

  esl_workqueue_Reset(queue);
#if DEBUGSERIAL
  printf("master threads reset\n");
#endif

  esl_threads_WaitForStart(obj);

#if DEBUGSERIAL
  printf("master threads started\n");
#endif

  status = esl_workqueue_ReaderUpdate(queue, NULL, &new_sq);
  if (status != eslOK) cm_Fail("Work queue reader failed");

#if DEBUGSERIAL
  printf("master initial update\n");
#endif 

  /* main loop: */
  for(i = 0; i < sq_block->count; i++) { 
    sq    = (ESL_SQ *) new_sq;
    sq    = sq_block->list + i;
    sq->W = sq_block->first_seqidx + i; 
    /* we overload sq->W w/seqidx (the original value is irrelevant in this context) */
    status = esl_workqueue_ReaderUpdate(queue, sq, &new_sq);
    if (status != eslOK) cm_Fail("Work queue reader failed");

#if DEBUGSERIAL
    printf("master internal update\n");
#endif
  }

  /* now send a empty sq to all workers signaling them to stop */
  empty_sq = esl_sq_Create();
  for(k = 0; k < nworkers; k++) { 
    status = esl_workqueue_ReaderUpdate(queue, empty_sq, &new_sq);
    if (status != eslOK) cm_Fail("Work queue reader failed");
#if DEBUGSERIAL
    printf("master termination update\n");
#endif
  }

  status = esl_workqueue_ReaderUpdate(queue, sq, NULL);
#if DEBUGSERIAL
  printf("master final update\n");
#endif

  /* wait for all the threads to complete */
  esl_threads_WaitForFinish(obj);
#if DEBUGSERIAL
  printf("master got finish\n");
#endif

  esl_workqueue_Complete(queue);  
#if DEBUGSERIAL
  printf("master completed\n");
#endif

  esl_sq_Destroy(empty_sq);
  return status;
}

/* pipeline_thread()
 * 
 * Receive a block of sequences from the master, 
 * align them and store their parsetrees.
 */
static void 
pipeline_thread(void *arg)
{
  int           status;
  int           i, j;
  int           workeridx;
  WORKER_INFO  *info;
  ESL_THREADS  *obj;
  ESL_SQ       *sq = NULL;
  void         *new_sq = NULL;
  char          errbuf[eslERRBUFSIZE];
  int           nalloc    = 0;
  int           allocsize = 1000;
  CM_P7_OM_HOLDER om_holder; /* reusable --p7pinbridge LOCAL profile/OPROFILE, per worker thread (brief 26_0430-090) */
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

#if DEBUGSERIAL
  printf("started thread %d\n", workeridx);
#endif

  info = (WORKER_INFO *) esl_threads_GetData(obj, workeridx);

#if DEBUGSERIAL
  printf("got data %d\n", workeridx);
#endif

  status = esl_workqueue_WorkerUpdate(info->queue, NULL, &new_sq);
  if (status != eslOK) cm_Fail("Work queue worker failed\n");

#if DEBUGSERIAL
  printf("initial update %d\n", workeridx);
#endif

  /* loop until all sequences have been processed */
  sq = (ESL_SQ *) new_sq;
  i = 0;
  cm_p7_om_holder_Init(&om_holder);
  while (sq->L != -1) {
    /* reallocate info->dataA if necessary */
    if(info->n == nalloc) { 
      ESL_REALLOC(info->dataA, sizeof(CM_ALNDATA *) * (nalloc + allocsize));
      for(j = nalloc; j < info->n + allocsize; j++) info->dataA[j] = NULL;
      nalloc += allocsize;
    }
    status = DispatchSqAlignment(info->cm, errbuf, sq, sq->W, info->mxsize,
				 TRMODE_UNKNOWN, info->pass_idx, FALSE, /* FALSE: info->cm->cp9b not valid */
				 info->w, info->w_tot, NULL, &om_holder, &(info->dataA[i]));
    /* sq->W has been overloaded (its original value is irrelevant in this context).
     * It is now the sequence index, defined in thread_loop() 
     */

    /* If alignment failed: potentially retry alignment in HMM banded
     * std (non-truncated) mode. We will only possibly do this if our
     * initial try was HMM banded truncated alignment (if not,
     * info->do_failover will be FALSE).
     */
    if(status == eslEAMBIGUOUS && info->do_failover == TRUE) { 
      assert(info->cm->align_opts & CM_ALIGN_TRUNC);
      info->cm->align_opts &= ~CM_ALIGN_TRUNC; /* lower truncated alignment flag, just for this sequence */
      status = DispatchSqAlignment(info->cm, errbuf, sq, sq->W, info->mxsize,
				   TRMODE_UNKNOWN, PLI_PASS_STD_ANY, FALSE, /* USE PLI_PASS_STD_ANY; FALSE: info->cm->cp9b not valid */
				   info->w, info->w_tot, NULL, &om_holder, &(info->dataA[i]));
      info->cm->align_opts |= CM_ALIGN_TRUNC; /* reraise truncated alignment flag */
    }
    if(status != eslOK) { 
      fprintf(stderr, "Problem during alignment of sequence %s\n", sq->name);
      cm_Fail(errbuf);
    }

    i++;
    info->n++;

    status = esl_workqueue_WorkerUpdate(info->queue, sq, &new_sq);
    if (status != eslOK) cm_Fail("Work queue worker failed");
    sq = (ESL_SQ *) new_sq;

#if DEBUGSERIAL
    printf("internal update %d sq->L: %" PRId64 "\n", workeridx, sq->L);
#endif
  }
  cm_p7_om_holder_Reset(&om_holder);

  status = esl_workqueue_WorkerUpdate(info->queue, sq, NULL);
  if (status != eslOK) cm_Fail("Work queue worker failed");
  
#if DEBUGSERIAL
  printf("final update %d\n", workeridx);
#endif

  esl_threads_Finished(obj, workeridx);
  return;

 ERROR: 
  cm_Fail("out of memory");
  return;  /* NEVERREACHED */
}

/* hmm_thread_loop()
 * 
 * Distribute sequences from sqarr[] to worker threads via work queue,
 * one sequence per work unit. Follows thread_loop() pattern but works
 * with a pre-read ESL_SQ** array instead of ESL_SQ_BLOCK.
 * Sequence index is passed via the sq->W overload.
 */
static int
hmm_thread_loop(WORKER_INFO *info, ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQ **sqarr, int nseq)
{
  int      status = eslOK;
  int      i, k;
  ESL_SQ  *sq;
  void    *new_sq;
  ESL_SQ  *empty_sq;
  int      nworkers = esl_threads_GetWorkerCount(obj);

  esl_workqueue_Reset(queue);
  esl_threads_WaitForStart(obj);

  status = esl_workqueue_ReaderUpdate(queue, NULL, &new_sq);
  if (status != eslOK) cm_Fail("Work queue reader failed");

  /* main loop: send each sequence to a worker */
  for (i = 0; i < nseq; i++) {
    sq    = sqarr[i];
    sq->W = i;  /* overload W with sequence index */
    status = esl_workqueue_ReaderUpdate(queue, sq, &new_sq);
    if (status != eslOK) cm_Fail("Work queue reader failed");
  }

  /* send empty sq to all workers signaling them to stop */
  empty_sq = esl_sq_Create();
  for (k = 0; k < nworkers; k++) {
    status = esl_workqueue_ReaderUpdate(queue, empty_sq, &new_sq);
    if (status != eslOK) cm_Fail("Work queue reader failed");
  }

  status = esl_workqueue_ReaderUpdate(queue, empty_sq, NULL);

  /* wait for all threads to complete */
  esl_threads_WaitForFinish(obj);
  esl_workqueue_Complete(queue);
  esl_sq_Destroy(empty_sq);

  return status;
}

/* hmm_pipeline_thread()
 *
 * Worker thread for --hmm alignment. Receives sequences from work queue,
 * computes p7 traces (Viterbi, unbanded OA, or Viterbi-banded OA),
 * and writes them directly into the shared tr[] array at tr[seqidx].
 * Follows pipeline_thread() pattern.
 */
static void 
hmm_pipeline_thread(void *arg)
{
  int           status;
  int           workeridx;
  WORKER_INFO  *info;
  ESL_THREADS  *obj;
  ESL_SQ       *sq = NULL;
  void         *new_sq = NULL;
  float         sc, fwdsc, oasc;
  char          errbuf[eslERRBUFSIZE];

  /* banded functions declared in cm_p7_band.c */
  extern int p7_kbands2gbands(int *i2k, int *kmin, int *kmax, int L, int M, P7_GBANDS **ret_bnd);
  extern int my_p7_GForwardBanded(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMXB *gxb, float *opt_sc);
  extern int p7_GBackwardBanded(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMXB *gxb, float *opt_sc);
  extern int p7_GDecodingBanded(const P7_PROFILE *gm, const P7_GMXB *fwd, P7_GMXB *bck, P7_GMXB *pp, float overall_sc);
  extern int p7_GOptimalAccuracyBanded(const P7_PROFILE *gm, const P7_GMXB *pp, P7_GMXB *gx, float *ret_e);
  extern int p7_GOATraceBanded(const P7_PROFILE *gm, const P7_GMXB *pp, const P7_GMXB *gx, P7_TRACE *tr);
  extern int p7_GCheckptFBDecode_Banded(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMXB *pp, float *ret_fwdsc); /* brief 26_0526-016 */
  extern int p7_GCheckptOA_Banded(const P7_PROFILE *gm, P7_GMXB *pp, P7_TRACE *tr, float *ret_oasc);                    /* brief 26_0526-016 */
  extern P7_GMXB *p7b_pp_Create(P7_GBANDS *bnd);                                                                       /* brief 26_0526-017: compact 2-cell resident pp */

#ifdef HAVE_FLUSH_ZERO_MODE
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif
  obj = (ESL_THREADS *) arg;
  esl_threads_Started(obj, &workeridx);
  info = (WORKER_INFO *) esl_threads_GetData(obj, workeridx);

  status = esl_workqueue_WorkerUpdate(info->queue, NULL, &new_sq);
  if (status != eslOK) cm_Fail("Work queue worker failed");

  /* loop until all sequences have been processed */
  sq = (ESL_SQ *) new_sq;
  while (sq->L != -1) {
    int idx = sq->W;  /* sequence index, overloaded by hmm_thread_loop */

    /* Reconfigure profile for this sequence length.
     * brief 26_0430-182 Part A: p7_ReconfigLength() unconditionally overwrites
     * xsc[N/C/J][MOVE|LOOP], clobbering the Tgm N->N/C->C loop-disable;
     * re-run the Tgm setup per-sequence instead when do_trunc. */
    if (info->do_trunc) {
      p7_ProfileConfig(info->hmm, info->bg, info->gm, sq->n, p7_LOCAL);
      p7_ProfileConfig5PrimeAnd3PrimeTrunc(info->gm, sq->n);
    } else {
      p7_ReconfigLength(info->gm, sq->n);
    }

    /* preflight: check HMM matrix size vs --mxsize before GrowTo.
     * Skipped under --p7ibv: the full P7_GMX is never allocated.
     * brief 26_0628-032: also skipped under kmeranchor/kmerchain (see serial-path
     * comment in hmm_alignment() for reasoning). */
    if (! info->do_p7ibv && ! (info->cm != NULL && (info->cm->p7_use_kmeranchor || info->cm->p7_use_kmerchain))) {
      double single_bytes = (double) sizeof(float) * (double)(info->hmm->M + 1) * (double)(sq->n + 1) * (double) p7G_NSCELLS;
      int    nmat         = info->do_hmmnoband ? 2 : 1;
      double needed_mb    = (single_bytes * (double) nmat) / (1024.0 * 1024.0);
      if (needed_mb > (double) info->mxsize) {
	int recommended_mxsize = (int)(ceil(needed_mb / 1024.0) * 1024.0);
	cm_Fail("HMM-only alignment mx needs %.2f Mb > %.2f Mb limit. Use --mxsize %d.",
		needed_mb, (double) info->mxsize, recommended_mxsize);
      }
    }

    if (info->do_hmmvit) {
      /* --- Viterbi trace --- */
      p7_gmx_GrowTo(info->gx, info->hmm->M, sq->n);
      p7_GViterbi(sq->dsq, sq->n, info->gm, info->gx, &sc);
      p7_trace_Reuse(info->hmm_tr[idx]);
      if ((status = p7_GTrace(sq->dsq, sq->n, info->gm, info->gx, info->hmm_tr[idx])) != eslOK)
	cm_Fail("p7_GTrace() failed for sequence %s", sq->name);
    }
    else if (info->do_hmmnoband) {
      /* --- Unbanded optimal accuracy --- */
      p7_gmx_GrowTo(info->gxf, info->hmm->M, sq->n);
      p7_gmx_GrowTo(info->gxb, info->hmm->M, sq->n);

      p7_GForward (sq->dsq, sq->n, info->gm, info->gxf, &fwdsc);
      p7_GBackward(sq->dsq, sq->n, info->gm, info->gxb, NULL);
      p7_GDecoding(info->gm, info->gxf, info->gxb, info->gxb);
      p7_GOptimalAccuracy(info->gm, info->gxb, info->gxf, &oasc);
      p7_trace_Reuse(info->hmm_tr[idx]);
      p7_GOATrace(info->gm, info->gxb, info->gxf, info->hmm_tr[idx]);
    }
    else {
      /* --- Viterbi-banded optimal accuracy --- */
      int     *i2k   = NULL;
      int     *kmin  = NULL;
      int     *kmax  = NULL;
      int      ncells = 0;
      int      pad   = 30;
      P7_GBANDS *bnd = NULL;
      P7_GMXB *bxf   = NULL;
      P7_GMXB *bxb   = NULL;
      P7_TRACE *vtr  = NULL;
      float    bwdsc = 0.;                                                    /* brief 26_0430-135b: capture backward total */
      int      p7ibv_delta = info->ibv_delta;                                 /* brief 26_0430-135b */
      int      do_widen = (getenv("P135B_FORCE_WIDEN") != NULL) ? TRUE : FALSE; /* brief 26_0430-135b widen override */

      if (info->do_p7ibv) {
	/* IBV D&C deriver: bands straight from cm->fp7, no full P7_GMX. */
	if (info->cm == NULL || info->cm->fp7 == NULL || info->cm->fp7->M != info->hmm->M)
	  cm_Fail("--hmm --p7ibv requires cm->fp7 with M matching the ML p7 HMM");
	if ((status = p7_Seq2BandsIBV_dnc(info->cm, errbuf, sq->dsq, sq->n,
					  p7ibv_delta, info->ibv_base_slab,
					  do_widen, /* brief 26_0430-135b: P135B_FORCE_WIDEN override; default FALSE (non-truncated --hmm) */
					  FALSE,    /* brief 26_0430-172: do_kband (unbanded D&C on --hmm path) */
					  info->do_trunc, /* brief 26_0430-182 Part B: was hardcoded FALSE; track CM_ALIGN_TRUNC like p7_ibv.c:1793 */
					  info->cm->p7_ibv_mode, info->cm->p7_ibv_width, /* brief 26_0430-140 */
					  &i2k, &kmin, &kmax, &ncells)) != eslOK)
	  cm_Fail("p7_Seq2BandsIBV_dnc() failed for sequence %s: %s", sq->name, errbuf);
      }
      else {
	/* brief 26_0628-032: k-mer anchor/chain deriver, opt-in via --p7kmeranchor/
	 * --p7kmerchain (mirrors cm_alndata.c:536-566's --p7band dispatch).
	 * brief 26_0628-038: do_trunc now threaded through, mirroring cm_alndata.c's
	 * CM-mode dispatch (cm->align_opts & CM_ALIGN_TRUNC). */
	int did_kmer = FALSE;
	int *local_nodepad = NULL;
	if (info->cm != NULL && (info->cm->p7_use_kmeranchor || info->cm->p7_use_kmerchain)
	    && (info->cm->flags & CMH_P7NODEPAD)) {
	  int k;
	  ESL_ALLOC(local_nodepad, sizeof(int) * (info->hmm->M + 1));
	  for (k = 0; k <= info->hmm->M; k++) local_nodepad[k] = info->cm->p7_cm_nodepad[k] + info->cm->p7bpad;
	}
	if (info->cm != NULL && info->cm->p7_use_kmeranchor) {
	  did_kmer = TRUE;
	  if ((status = p7_Seq2BandsKmerAnchor(info->cm, errbuf, sq->dsq, sq->n, local_nodepad,
					       info->do_trunc, /* brief 26_0628-038: track CM_ALIGN_TRUNC like cm_alndata.c:541 */
					       &i2k, &kmin, &kmax, &ncells)) != eslOK)
	    cm_Fail("p7_Seq2BandsKmerAnchor() failed for sequence %s: %s", sq->name, errbuf);
	}
	else if (info->cm != NULL && info->cm->p7_use_kmerchain) {
	  did_kmer = TRUE;
	  if ((status = p7_Seq2BandsKmerChain(info->cm, errbuf, sq->dsq, sq->n, local_nodepad,
					      info->do_trunc, /* brief 26_0628-038: track CM_ALIGN_TRUNC like cm_alndata.c:558 */
					      &i2k, &kmin, &kmax, &ncells)) != eslOK)
	    cm_Fail("p7_Seq2BandsKmerChain() failed for sequence %s: %s", sq->name, errbuf);
	}
	if (local_nodepad) free(local_nodepad);

	if (did_kmer && ncells == 0 && ! info->cm->p7_kmerchain_fallback_vit) {
	  /* brief 26_0628-047: M-gate/N-gate fired, or no anchor found -- default
	   * fallback target is --p7ibv's D&C deriver instead of a Vit-trace
	   * band (see serial hmm_alignment()'s matching comment above).
	   * brief 26_0628-050: use info->cm->p7_ibv_delta (struct default 3000),
	   * not info->ibv_delta (--p7ibv-delta's CLI default 20000) -- see the
	   * serial hmm_alignment() comment above for why. */
	  if ((status = kmer_gate_p7ibv_fallback(info->cm, errbuf, sq->dsq, sq->n, info->do_trunc,
						  info->cm->p7_ibv_delta, info->ibv_base_slab,
						  &i2k, &kmin, &kmax, &ncells)) != eslOK)
	    cm_Fail("kmer_gate_p7ibv_fallback() failed for sequence %s: %s", sq->name, errbuf);
	}
	if (! did_kmer || ncells == 0) {
	  /* Default banded-OA path (no kmer flag set); the kmeranchor/kmerchain
	   * ncells==0 fallback when --p7kmerchain-fbvit reverts to the
	   * old behavior; and the safety net when the --p7ibv fallback above
	   * itself also found nothing usable (mirrors cm_alndata.c's
	   * kmeranchor->vitband / kmerchain->vitband fallback shape). */
	  vtr = p7_trace_Create();
	  p7_gmx_GrowTo(info->gx, info->hmm->M, sq->n);
	  p7_GViterbi(sq->dsq, sq->n, info->gm, info->gx, &sc);
	  p7_trace_Reuse(vtr);
	  status = p7_GTrace(sq->dsq, sq->n, info->gm, info->gx, vtr);

	  if (status != eslOK || vtr->N == 0) {
	    /* Viterbi failed; fall back to unbanded OA */
	    P7_GMX *fallback_gxf = p7_gmx_Create(info->hmm->M, sq->n);
	    P7_GMX *fallback_gxb = p7_gmx_Create(info->hmm->M, sq->n);
	    p7_GForward (sq->dsq, sq->n, info->gm, fallback_gxf, &fwdsc);
	    p7_GBackward(sq->dsq, sq->n, info->gm, fallback_gxb, NULL);
	    p7_GDecoding(info->gm, fallback_gxf, fallback_gxb, fallback_gxb);
	    p7_GOptimalAccuracy(info->gm, fallback_gxb, fallback_gxf, &oasc);
	    p7_trace_Reuse(info->hmm_tr[idx]);
	    p7_GOATrace(info->gm, fallback_gxb, fallback_gxf, info->hmm_tr[idx]);
	    p7_gmx_Destroy(fallback_gxf);
	    p7_gmx_Destroy(fallback_gxb);
	    p7_trace_Destroy(vtr);
	    goto HMM_NEXT_SQ;
	  }

	  {
	    int tpos;
	    ESL_ALLOC(i2k, sizeof(int) * (sq->n + 1));
	    esl_vec_ISet(i2k, (sq->n + 1), -1);
	    for (tpos = 0; tpos < vtr->N; tpos++) {
	      if (vtr->st[tpos] == p7T_M) {
		int i = vtr->i[tpos];
		int k = vtr->k[tpos];
		if (i >= 1 && i <= sq->n && k >= 1 && k <= info->hmm->M)
		  i2k[i] = k;
	      }
	    }
	  }

	  if ((status = p7_pins2bands(i2k, errbuf, sq->n, info->hmm->M, pad, &kmin, &kmax, &ncells)) != eslOK)
	    cm_Fail("p7_pins2bands() failed for sequence %s: %s", sq->name, errbuf);
	}
      }
      if ((status = p7_kbands2gbands(i2k, kmin, kmax, sq->n, info->hmm->M, &bnd)) != eslOK)
	cm_Fail("p7_kbands2gbands() failed for sequence %s", sq->name);
      if (getenv("BRIEF035_MEMPOINT") != NULL)
	fprintf(stderr, "#MEMPOINT after_gbands_threaded seq=%s L=%d rss_kb=%ld\n", sq->name, (int) sq->n, brief035_rss_kb());

      if (getenv("INFERNAL_HMM_CKPT_OFF") == NULL) {
	/* brief 26_0628-036: port of brief 26_0526-016's sqrt(nrow)-checkpointed F/B/Decode/OA/
	 * traceback into the threaded worker (mirrors hmm_alignment()'s serial-path
	 * gate above; see cm_p7_band.c:8641/9146). bxb holds the resident posterior;
	 * bxf is never allocated. */
	bxb = p7b_pp_Create(bnd);
	if (getenv("BRIEF035_MEMPOINT") != NULL)
	  fprintf(stderr, "#MEMPOINT after_cp9alloc_ckpt_threaded seq=%s L=%d rss_kb=%ld\n", sq->name, (int) sq->n, brief035_rss_kb());
	if ((status = p7_GCheckptFBDecode_Banded(sq->dsq, sq->n, info->gm, bxb, &fwdsc)) != eslOK)
	  cm_Fail("p7_GCheckptFBDecode_Banded() failed for sequence %s", sq->name);
	p7_trace_Reuse(info->hmm_tr[idx]);
	if ((status = p7_GCheckptOA_Banded(info->gm, bxb, info->hmm_tr[idx], &oasc)) != eslOK)
	  cm_Fail("p7_GCheckptOA_Banded() failed for sequence %s", sq->name);
      }
      else {
      bxf = p7_gmxb_Create(bnd);
      bxb = p7_gmxb_Create(bnd);
      if (getenv("BRIEF035_MEMPOINT") != NULL)
	fprintf(stderr, "#MEMPOINT after_cp9alloc_nockpt_threaded seq=%s L=%d rss_kb=%ld\n", sq->name, (int) sq->n, brief035_rss_kb());

      if ((status = my_p7_GForwardBanded(sq->dsq, sq->n, info->gm, bxf, &fwdsc)) != eslOK)
	cm_Fail("my_p7_GForwardBanded() failed for sequence %s", sq->name);
      if ((status = p7_GBackwardBanded(sq->dsq, sq->n, info->gm, bxb, &bwdsc)) != eslOK)
	cm_Fail("p7_GBackwardBanded() failed for sequence %s", sq->name);
      if (getenv("P135B_FB_INSTRUMENT") != NULL)
	fprintf(stderr, "#P135B_FBTOTAL seq=%s M=%d L=%d delta=%d widen=%d ncells=%d fwd=%.6f bwd=%.6f gap=%.6f\n",
		sq->name, info->hmm->M, (int) sq->n, p7ibv_delta, do_widen, ncells, fwdsc, bwdsc, fwdsc - bwdsc);
      if ((status = p7_GDecodingBanded(info->gm, bxf, bxb, bxb, fwdsc)) != eslOK)
	cm_Fail("p7_GDecodingBanded() failed for sequence %s", sq->name);
      if ((status = p7_GOptimalAccuracyBanded(info->gm, bxb, bxf, &oasc)) != eslOK)
	cm_Fail("p7_GOptimalAccuracyBanded() failed for sequence %s", sq->name);

      p7_trace_Reuse(info->hmm_tr[idx]);
      if ((status = p7_GOATraceBanded(info->gm, bxb, bxf, info->hmm_tr[idx])) != eslOK)
	cm_Fail("p7_GOATraceBanded() failed for sequence %s", sq->name);
      }
      if (getenv("BRIEF035_MEMPOINT") != NULL)
	fprintf(stderr, "#MEMPOINT alignment_peak_threaded seq=%s L=%d rss_kb=%ld\n", sq->name, (int) sq->n, brief035_rss_kb());

      free(i2k);
      free(kmin);
      free(kmax);
      if (vtr) p7_trace_Destroy(vtr);
      p7_gbands_Destroy(bnd);
      if (bxf) p7_gmxb_Destroy(bxf);
      p7_gmxb_Destroy(bxb);
    }

  HMM_NEXT_SQ:
    status = esl_workqueue_WorkerUpdate(info->queue, sq, &new_sq);
    if (status != eslOK) cm_Fail("Work queue worker failed");
    sq = (ESL_SQ *) new_sq;
  }

  status = esl_workqueue_WorkerUpdate(info->queue, sq, NULL);
  if (status != eslOK) cm_Fail("Work queue worker failed");

  esl_threads_Finished(obj, workeridx);
  return;

 ERROR:
  cm_Fail("out of memory");
  return;  /* NEVERREACHED */
}
#endif   /* HMMER_THREADS */

#if HAVE_MPI

/* P7_TRACE MPI serialization for --hmm mode.
 * Pack/unpack a P7_TRACE + sequence index into an MPI buffer.
 * Only the fields needed by p7_tracealign_Seqs() are transmitted:
 * N, st[], k[], i[], pp[] (if present), L, M.
 */
static int
hmm_trace_MPIPackSize(P7_TRACE *tr, MPI_Comm comm, int *ret_n)
{
  int n = 0, sz;
  MPI_Pack_size(1,     MPI_LONG_LONG_INT, comm, &sz); n += sz; /* idx     */
  MPI_Pack_size(1,     MPI_INT,           comm, &sz); n += sz; /* N       */
  MPI_Pack_size(1,     MPI_INT,           comm, &sz); n += sz; /* has_pp  */
  MPI_Pack_size(1,     MPI_INT,           comm, &sz); n += sz; /* L       */
  MPI_Pack_size(1,     MPI_INT,           comm, &sz); n += sz; /* M       */
  MPI_Pack_size(tr->N, MPI_CHAR,          comm, &sz); n += sz; /* st[]    */
  MPI_Pack_size(tr->N, MPI_INT,           comm, &sz); n += sz; /* k[]     */
  MPI_Pack_size(tr->N, MPI_INT,           comm, &sz); n += sz; /* i[]     */
  if (tr->pp != NULL) {
    MPI_Pack_size(tr->N, MPI_FLOAT,       comm, &sz); n += sz; /* pp[]    */
  }
  *ret_n = n;
  return eslOK;
}

static int
hmm_trace_MPIPack(P7_TRACE *tr, int64_t idx, char *buf, int n, int *pos, MPI_Comm comm)
{
  int has_pp = (tr->pp != NULL) ? 1 : 0;
  MPI_Pack(&idx,     1,     MPI_LONG_LONG_INT, buf, n, pos, comm);
  MPI_Pack(&(tr->N), 1,     MPI_INT,           buf, n, pos, comm);
  MPI_Pack(&has_pp,  1,     MPI_INT,           buf, n, pos, comm);
  MPI_Pack(&(tr->L), 1,     MPI_INT,           buf, n, pos, comm);
  MPI_Pack(&(tr->M), 1,     MPI_INT,           buf, n, pos, comm);
  MPI_Pack(tr->st,   tr->N, MPI_CHAR,          buf, n, pos, comm);
  MPI_Pack(tr->k,    tr->N, MPI_INT,           buf, n, pos, comm);
  MPI_Pack(tr->i,    tr->N, MPI_INT,           buf, n, pos, comm);
  if (has_pp) MPI_Pack(tr->pp, tr->N, MPI_FLOAT, buf, n, pos, comm);
  return eslOK;
}

static int
hmm_trace_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, P7_TRACE **ret_tr, int64_t *ret_idx)
{
  int        status;
  int64_t    idx;
  int        N, has_pp, L, M;
  P7_TRACE  *tr = NULL;

  MPI_Unpack(buf, n, pos, &idx,    1, MPI_LONG_LONG_INT, comm);
  MPI_Unpack(buf, n, pos, &N,      1, MPI_INT,           comm);
  MPI_Unpack(buf, n, pos, &has_pp, 1, MPI_INT,           comm);
  MPI_Unpack(buf, n, pos, &L,      1, MPI_INT,           comm);
  MPI_Unpack(buf, n, pos, &M,      1, MPI_INT,           comm);

  tr = has_pp ? p7_trace_CreateWithPP() : p7_trace_Create();
  if (tr == NULL) return eslEMEM;
  if ((status = p7_trace_GrowTo(tr, N)) != eslOK) { p7_trace_Destroy(tr); return status; }
  tr->N = N;
  tr->L = L;
  tr->M = M;

  MPI_Unpack(buf, n, pos, tr->st, N, MPI_CHAR,  comm);
  MPI_Unpack(buf, n, pos, tr->k,  N, MPI_INT,   comm);
  MPI_Unpack(buf, n, pos, tr->i,  N, MPI_INT,   comm);
  if (has_pp) MPI_Unpack(buf, n, pos, tr->pp, N, MPI_FLOAT, comm);

  *ret_tr  = tr;
  *ret_idx = idx;
  return eslOK;
}

static int
hmm_trace_MPISend(P7_TRACE *tr, int64_t idx, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc)
{
  int   n = 0;
  int   pos;

  hmm_trace_MPIPackSize(tr, comm, &n);
  if (n > *nalloc) {
    void *tmp;
    ESL_RALLOC(*buf, tmp, sizeof(char) * n);
    *nalloc = n;
  }
  pos = 0;
  hmm_trace_MPIPack(tr, idx, *buf, n, &pos, comm);
  MPI_Send(*buf, n, MPI_PACKED, dest, tag, comm);
  return eslOK;

 ERROR:
  return eslEMEM;
}

/* mpi_master()
 * The MPI version of cmalign.
 * Follows standard pattern for a master/worker load-balanced MPI program 
 * (SRE notes J1/78-79).
 * 
 * A master returns eslOK if it's successful.  Errors in an MPI master
 * come in two classes: recoverable and nonrecoverable.  
 * 
 * Recoverable errors include most worker-side errors, and any
 * master-side error that do not affect MPI communication. Error
 * messages from recoverable messages are delayed until we've cleanly
 * shut down the workers. The 
 * 
 * Some worker side errors (such as ESL_ALLOCs) are likely to be 
 * unrecoverable and will almost certainly cause MPI to crash
 * uncleanly, they're only here because I couldn't find a way around
 * them without massive reimplementation. Hopefully they rarely occur.
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
  int             status;                /* Easel status */
  char            errbuf[eslERRBUFSIZE]; /* for printing error messages */
  CM_t           *cm = NULL;             /* a CM */
  int             i;                     /* counter over parsetrees */
  int             si;                    /* sequence index */
  int             nali;                  /* index of the (possibly temporary) alignment we are working on */
  int             nseq_cur;              /* number of sequences in current alignment */
  int             nseq_aligned;          /* number of sequences so far aligned */
  /* MPI is incompatible with --sample, b/c it would not be reproducible */

  /* variables related to output, we may use a tmpfile if seqfile is large */
  int      use_tmpfile;              /* print out current alignment to tmpfile? */
  int      created_tmpfile = FALSE;  /* TRUE if we've created a tmp file for current CM */
  char tmpfile[32] = "esltmpXXXXXX"; /* name of the tmpfile */
  CM_ALNDATA **merged_dataA = NULL;  /* array of all CM_ALNDATA pointers for current alignment */
  int          merged_data_idx = 0;  /* index in merged_dataA */
  int          nmerged;              /* size of merged_dataA */

  /* variables related to reading sequence blocks */
  int            sstatus = eslOK;  /* status from esl_sq_ReadBlock() */
  ESL_SQ_BLOCK  *sq_block;         /* a sequence block */
  ESL_SQ_BLOCK  *nxt_sq_block;     /* sequence block for next loop iteration */
  int            reached_eof;      /* TRUE if we've reached EOF in target sequence file */

  /* variables related to --mapali */
  char         *map_file   = NULL; /* name of alignment file from --mapali */
  CM_ALNDATA  **map_dataA  = NULL; /* array of CM_ALNDATA pointers for mapali alignment */
  int           nmap_data  = 0;    /* number of CM_ALNDATA ptrs in map_dataA */
  int           nmap_cur   = 0;    /* number of CM_ALNDATA ptrs to include in current iteration, 0 unless nali==0 */
  char         *map_sscons = NULL; /* SS_cons from mapali, only used if --mapstr */

  /* variables related to MPI implementation */
  MPI_Status       mpistatus;       /* the mpi status */
  char            *mpibuf  = NULL;  /* buffer used to pack/unpack structures */
  int              mpibuf_size = 0; /* current size of mpibuf_size */
  int              buf_size;        /* size of received buffer */
  int              pos;             /* for packing/unpacking an MPI buffer */
  int              wi;              /* worker index that we're about to send to or receive from */
  int              nworkers;        /* number of workers */
  int              nworking;        /* number of workers currently doing work */
  int              have_work;       /* TRUE while work remains (sqs remain in sq_block) */
  CM_ALNDATA      *wkr_data = NULL; /* data recieved from a worker */
  ESL_SQ          *sq = NULL;       /* sequence to send to a worker */

  /* See 'General notes on {serial,mpi}_master()'s strategy' for details
   * on the code organization here. 
   */

  if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) mpi_failure(errbuf);
  if(esl_opt_GetBoolean(go, "--sample")) mpi_failure("--sample does not work with in MPI mode (b/c results would not be exactly reproducible)");

  /* Read one CM, and make sure there's only one. This fills cfg->abc. */
  status = cm_file_Read(cfg->cmfp, TRUE, &(cfg->abc), &cm);
  if(status != eslOK) mpi_failure(cfg->cmfp->errbuf);
  status = cm_file_Read(cfg->cmfp, TRUE, &(cfg->abc), NULL);
  if(status != eslEOF) mpi_failure("CM file %s does not contain just one CM\n", cfg->cmfp->fname);

  nworkers  = cfg->nproc - 1;
  if(cfg->ofp != stdout) output_header(stdout, go, cfg->cmfile, cfg->sqfile, cm, nworkers+1);

  /* initialization */
  nali = nseq_cur = nseq_aligned = 0;
  if((status = initialize_cm(go, cfg, errbuf, cm)) != eslOK) mpi_failure(errbuf);

  /* --hmm mode: HMM-only alignment via MPI.
   * Master reads all sequences, distributes dsqs to workers,
   * workers compute p7 traces and send them back, master
   * collects all traces and calls p7_tracealign_Seqs().
   */
  if(esl_opt_GetBoolean(go, "--hmm")) {
    ESL_SQ      **sqarr  = NULL;
    P7_TRACE    **tr     = NULL;
    int           nseq   = 0;
    int           nalloc = 256;
    int           idx2;
    int           do_hmmvit = (cm->align_opts & CM_ALIGN_P7HMMVIT) ? TRUE : FALSE;

    /* Read all sequences */
    ESL_ALLOC(sqarr, sizeof(ESL_SQ *) * nalloc);
    while (1) {
      ESL_SQ *sq_tmp = esl_sq_CreateDigital(cfg->abc);
      status = esl_sqio_Read(cfg->sqfp, sq_tmp);
      if (status == eslEOF) { esl_sq_Destroy(sq_tmp); break; }
      if (status != eslOK)  mpi_failure("Error reading sequence file");
      if (nseq >= nalloc) { nalloc *= 2; ESL_REALLOC(sqarr, sizeof(ESL_SQ *) * nalloc); }
      sqarr[nseq++] = sq_tmp;
    }
    if (nseq == 0) mpi_failure("No sequences found");

    /* Allocate trace array */
    ESL_ALLOC(tr, sizeof(P7_TRACE *) * nseq);
    for (idx2 = 0; idx2 < nseq; idx2++) tr[idx2] = NULL;

    /* Distribute sequences to workers and collect traces */
    have_work = TRUE;
    nworking  = 0;
    si        = 0;
    while(have_work || nworking > 0) {
      if (MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus) != 0)
	mpi_failure("MPI error receiving message");
      if (MPI_Get_count(&mpistatus, MPI_PACKED, &buf_size) != 0)
	mpi_failure("MPI get count failed");
      if (mpibuf == NULL || buf_size > mpibuf_size) {
	ESL_REALLOC(mpibuf, sizeof(char) * buf_size);
	mpibuf_size = buf_size;
      }
      wi = mpistatus.MPI_SOURCE;
      MPI_Recv(mpibuf, buf_size, MPI_PACKED, wi, mpistatus.MPI_TAG, MPI_COMM_WORLD, &mpistatus);

      if (mpistatus.MPI_TAG == INFERNAL_ALNDATA_TAG) {
	/* Receive trace from worker */
	int64_t recv_idx;
	P7_TRACE *recv_tr = NULL;
	pos = 0;
	status = hmm_trace_MPIUnpack(mpibuf, buf_size, &pos, MPI_COMM_WORLD, &recv_tr, &recv_idx);
	if (status != eslOK) mpi_failure("problem unpacking trace from worker %d", wi);
	tr[recv_idx] = recv_tr;
	nworking--;
      }
      else if (mpistatus.MPI_TAG == INFERNAL_ERROR_TAG) {
	mpi_failure("MPI client %d raised error:\n%s\n", wi, mpibuf);
      }
      else if (mpistatus.MPI_TAG != INFERNAL_INITIALREADY_TAG) {
	mpi_failure("Unexpected tag %d from %d\n", mpistatus.MPI_TAG, wi);
      }

      if (have_work) {
	sq = sqarr[si];
	status = cm_dsq_MPISend(sq->dsq, sq->L, (int64_t)si, wi, INFERNAL_DSQ_TAG, MPI_COMM_WORLD, &mpibuf, &mpibuf_size);
	if (status != eslOK) mpi_failure("problem sending dsq to worker %d", wi);
	nworking++;
	si++;
	if (si == nseq) have_work = FALSE;
      }
    }

    /* Tell workers we're done: send NULL dsq for end-of-block, then again for end-of-file */
    for (wi = 1; wi < cfg->nproc; wi++)
      cm_dsq_MPISend(NULL, -1, -1, wi, INFERNAL_DSQ_TAG, MPI_COMM_WORLD, &mpibuf, &mpibuf_size);
    for (wi = 1; wi < cfg->nproc; wi++)
      cm_dsq_MPISend(NULL, -1, -1, wi, INFERNAL_DSQ_TAG, MPI_COMM_WORLD, &mpibuf, &mpibuf_size);

    /* Build MSA from traces (reuses hmm_alignment's post-processing logic) */
    {
      P7_HMM  *hmm_mpi = cm->mlp7;
      ESL_MSA *msa = NULL;

      if ((status = p7_tracealign_Seqs(sqarr, tr, nseq, hmm_mpi->M, p7_ALL_CONSENSUS_COLS, hmm_mpi, &msa)) != eslOK)
	mpi_failure("p7_tracealign_Seqs() failed");

      /* Add SS_cons */
      if (cm->cmcons != NULL && cm->cmcons->cstr != NULL && msa->rf != NULL) {
	int cpos, apos;
	ESL_ALLOC(msa->ss_cons, sizeof(char) * (msa->alen + 1));
	cpos = 0;
	for (apos = 0; apos < msa->alen; apos++) {
	  if (msa->rf[apos] != '.' && msa->rf[apos] != '~') {
	    msa->ss_cons[apos] = (cpos < cm->clen) ? cm->cmcons->cstr[cpos] : '.';
	    cpos++;
	  } else {
	    msa->ss_cons[apos] = '.';
	  }
	}
	msa->ss_cons[msa->alen] = '\0';
      }

      /* --dnaout */
      if (esl_opt_GetBoolean(go, "--dnaout") && cfg->abc_out->type == eslDNA) {
	int i3, apos;
	for (i3 = 0; i3 < msa->nseq; i3++)
	  for (apos = 0; apos < msa->alen; apos++) {
	    if (msa->aseq[i3][apos] == 'U') msa->aseq[i3][apos] = 'T';
	    if (msa->aseq[i3][apos] == 'u') msa->aseq[i3][apos] = 't';
	  }
      }

      esl_msafile_Write(cfg->ofp, msa, cfg->outfmt);
      esl_msa_Destroy(msa);
    }

    /* Clean up */
    for (idx2 = 0; idx2 < nseq; idx2++) {
      if (tr[idx2] != NULL) p7_trace_Destroy(tr[idx2]);
      esl_sq_Destroy(sqarr[idx2]);
    }
    free(tr);
    free(sqarr);
    if (mpibuf != NULL) free(mpibuf);
    FreeCM(cm);
    return eslOK;
  }

  reached_eof = FALSE;

  /* include the mapali, if nec */
  if((map_file = esl_opt_GetString(go, "--mapali")) != NULL) { 
    if((status = map_alignment(map_file, cm, esl_opt_GetBoolean(go, "--noss"), errbuf, &map_dataA, &nmap_data, &map_sscons)) != eslOK) mpi_failure(errbuf);
    if(esl_opt_GetBoolean(go, "--mapstr") && map_sscons == NULL) mpi_failure("Failed to read SS_cons for --mapstr from %s", map_file);
  }

  /* Our main loop will loop over reading a single large block
   * (<sq_block>) of sequences, up to MAX_RESIDUE_COUNT
   * (1024*1024=1048576) residues, and up to CMALIGN_MAX_NSEQ
   * sequences (10,000), but potentially less if we reach the end of
   * the sequence file first.
   */

  /* Read the first block */
  sq_block = esl_sq_CreateDigitalBlock(CMALIGN_MAX_NSEQ, cfg->abc);
  sstatus = esl_sqio_ReadBlock(cfg->sqfp, sq_block, -1, -1, /*max_init_window=*/FALSE, FALSE); /* FALSE says: read complete sequences */
  nxt_sq_block = sq_block; /* special case of first block read */

#if DEBUGMPI
  printf("master read the first block\n");
#endif 

  while(sstatus == eslOK) { 
    sq_block = nxt_sq_block; /* our current sq_block becomes the one we read on the previous iteration */
    sq_block->first_seqidx = nseq_aligned;
    nseq_cur = sq_block->count;

    /* Before we do any aligning, read the next sequence block, so we
     * can determine if we've reached the end of the seqfile. We need
     * to know this for two reasons:
     *
     * (1) if the first block read above included all sequences (which
     * we won't know until we try to read another block), we don't
     * need to go into memory-saving mode and output to a tmpfile, we
     * can output (in interleaved mode) to the final output file.
     *
     * (2) if do_oneblock (we're trying to output the full alignment
     * as a single block) we need to fail if we still have sequences
     * left, because the sequence file exceeded the size limits. And
     * we want to fail *before* we align all the sequences, so the
     * user isn't cross when the job fails after seemingly going along
     * fine for a while.
     */
    nxt_sq_block = esl_sq_CreateDigitalBlock(CMALIGN_MAX_NSEQ, cfg->abc);
    sstatus = esl_sqio_ReadBlock(cfg->sqfp, nxt_sq_block, -1, -1, /*max_init_window=*/FALSE, FALSE); /* FALSE says: read complete sequences */
    if(sstatus == eslEOF) { 
      reached_eof = TRUE; /* nxt_sq_block will not have been filled */
      esl_sq_DestroyBlock(nxt_sq_block); 
    }
    if((! reached_eof) && cfg->do_oneblock) esl_fatal("Error: the sequence file is too big (has > %d seqs or %d residues) for --ileaved or\noutput format other than Pfam. Use esl-reformat to reformat alignment later.", CMALIGN_MAX_NSEQ, MAX_RESIDUE_COUNT);

    /* allocate an array for all CM_ALNDATA objects we'll receive from workers */
    nmap_cur = (nali == 0) ? nmap_data : 0;
    nmerged = nseq_cur + nmap_cur;
    ESL_ALLOC(merged_dataA, sizeof(CM_ALNDATA *) * ESL_MAX(1, nseq_cur + nmap_cur)); // avoid malloc of 0
    /* prepend mapali data if nec */
    if(nmap_cur > 0) {
      for(i = 0; i < nmap_cur; i++) merged_dataA[i] = map_dataA[i];
      free(map_dataA); /* don't free the CM_ALNDATA objects, merged_dataA is pointing at them */
      map_dataA = NULL;
    }
#if DEBUGMPI
    printf("master about to enter main loop\n");
#endif 

    /* main send/recv loop: send sequences to workers and receive their results */
    have_work = TRUE;
    nworking  = 0;
    si        = 0; /* sequence index */
    while(have_work || nworking > 0) { 
#if DEBUGMPI
      printf("master waiting for a message from any worker\n");
#endif	
      /* wait for message (results, ready tag or error) from any worker */
      if (MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus) != 0) mpi_failure("MPI error %d receiving message from %d\n", mpistatus.MPI_SOURCE);
      if (MPI_Get_count(&mpistatus, MPI_PACKED, &buf_size)                   != 0) mpi_failure("MPI get count failed");;
      if (mpibuf == NULL || buf_size > mpibuf_size) {
	ESL_REALLOC(mpibuf, sizeof(char) * buf_size);
	mpibuf_size = buf_size; 
      }
      wi = mpistatus.MPI_SOURCE;
      MPI_Recv(mpibuf, buf_size, MPI_PACKED, wi, mpistatus.MPI_TAG, MPI_COMM_WORLD, &mpistatus);
      
#if DEBUGMPI
      printf("master received message from worker %d, tag %d\n", wi, mpistatus.MPI_TAG);
#endif	
      /* tag should be either:
       * INFERNAL_INITIALREADY_TAG: worker just initialized, and is ready for work 
       * INFERNAL_ALNDATA_TAG:      worker finished work, and sent us results 
       * INFERNAL_ERROR_TAG:        worker sent us an error
       */
      if (mpistatus.MPI_TAG == INFERNAL_ALNDATA_TAG) { 
	/* receive CM_ALNDATA result from worker, and point merged_dataA at it */
	pos = 0;
	status = cm_alndata_MPIUnpack(mpibuf, buf_size, &pos, MPI_COMM_WORLD, cfg->abc, &wkr_data);
	if(status != eslOK) mpi_failure("problem with alignment results received from worker %d", wi);
	merged_data_idx = wkr_data->idx - nseq_aligned + nmap_cur;
	merged_dataA[merged_data_idx]     = wkr_data; /* wkr_data we received had everything we need except the sq */
	merged_dataA[merged_data_idx]->sq = sq_block->list + (wkr_data->idx - sq_block->first_seqidx); /* point sq at appropriate sequence */
	nworking--; /* one less worker is working now */
#if DEBUGMPI
	printf("received results from worker %d, %d workers now working\n", wi, nworking);
#endif	
      }
      else if(mpistatus.MPI_TAG == INFERNAL_ERROR_TAG) { 
	mpi_failure("MPI client %d raised error:\n%s\n", wi, mpibuf);
      }
      else if (mpistatus.MPI_TAG != INFERNAL_INITIALREADY_TAG) { 
	mpi_failure("Unexpected tag %d from %d\n", mpistatus.MPI_TAG, wi);
      }
      
      if(have_work) { /* send new sequence: si's dsq, L, and seqidx to the worker */
	sq = sq_block->list + si;
	status = cm_dsq_MPISend(sq->dsq, sq->L, sq_block->first_seqidx + si, wi, INFERNAL_DSQ_TAG, MPI_COMM_WORLD, &mpibuf, &mpibuf_size);
	if(status != eslOK) mpi_failure("problem sending dsq to worker %d", wi);
	nworking++; /* one more worker is working now */
	si++;       /* move onto next sequence */
#if DEBUGMPI
	printf("master sent dsq si: %d/%d to worker %d, %d workers now working\n", si, sq_block->count, wi, nworking);
#endif
	if(si == sq_block->count) { 
	  have_work = FALSE;
	  ESL_DPRINTF1(("#DEBUG: MPI master has sent all %d of its sequences\n", sq_block->count));
#if DEBUGMPI
	  printf("master is out of work\n");
#endif 
	}
      }
    }
    if(nworking != 0) mpi_failure("%d workers still working when all should be idle", nworking);

    /* we're done with sq_block, tell the workers by sending a NULL dsq */
#if DEBUGMPI
    printf("master sending NULL dsq to all workers signalling we're done with the file\n");
#endif
    for (wi = 1; wi < cfg->nproc; wi++) {
      if((status = cm_dsq_MPISend(NULL, -1, -1, wi, INFERNAL_DSQ_TAG, MPI_COMM_WORLD, &mpibuf, &mpibuf_size)) != eslOK) mpi_failure(errbuf);
    }

    /* output alignment (if do_oneblock we died above if we didn't reach EOF yet) */
    use_tmpfile = (reached_eof && (! created_tmpfile)) ? FALSE : TRUE; /* output to tmpfile only if this is the first alignment and we've aligned all seqs */
    if(use_tmpfile && (! created_tmpfile)) { 
      /* first aln for temporary output file, open the file */	
      if ((status = esl_tmpfile_named(tmpfile, &(cfg->tmpfp))) != eslOK) mpi_failure("Failed to open temporary output file (status %d)", status);
      created_tmpfile = TRUE;
    }
    if((status   = output_alignment(go, cfg, errbuf, cm, (use_tmpfile ? cfg->tmpfp : cfg->ofp), merged_dataA, nseq_cur + nmap_cur, map_sscons)) != eslOK) cm_Fail(errbuf);
    /* optionally output same alignment to regress file */
    if(cfg->rfp != NULL) { 
      if((status = output_alignment(go, cfg, errbuf, cm, cfg->rfp,                              merged_dataA, nseq_cur + nmap_cur, map_sscons)) != eslOK) mpi_failure(errbuf);
    }    
    nali++;
    nseq_aligned += nseq_cur;

    /* output scores to stdout, if -o used */
    if(cfg->ofp != stdout) { 
      if((status =  output_scores(stdout,   cm, errbuf, merged_dataA, nseq_cur + nmap_cur, nmap_cur, cfg->be_verbose)) != eslOK) mpi_failure(errbuf);
    }
    /* output scores to scores file, if --sfp used */
    if(cfg->sfp != NULL) { 
      if(nali == 1) output_header(stdout, go, cfg->cmfile, cfg->sqfile, cm, nworkers+1);
      if((status =  output_scores(cfg->sfp, cm, errbuf, merged_dataA, nseq_cur + nmap_cur, nmap_cur, cfg->be_verbose)) != eslOK) mpi_failure(errbuf);
    }

    /* free block and worker data */
    esl_sq_DestroyBlock(sq_block);
    sq_block = NULL;
    for(i = 0; i < nmerged; i++) { 
      cm_alndata_Destroy(merged_dataA[i], FALSE); /* FALSE: don't free sq's, we just free'd them by destroying the block */
    }
    free(merged_dataA);
  } /* end of outer while loop 'while(sstatus == eslOK)' */

  /* done with all sequences/blocks in the sequence file, tell the workers */
#if DEBUGMPI
  printf("master sending NULL dsq to all workers signalling we're done with the file\n");
#endif
  for (wi = 1; wi < cfg->nproc; wi++) {
    if((status = cm_dsq_MPISend(NULL, -1, -1, wi, INFERNAL_DSQ_TAG, MPI_COMM_WORLD, &mpibuf, &mpibuf_size)) != eslOK) mpi_failure(errbuf);
  }
  
  if     (sstatus == eslEFORMAT) mpi_failure("Parse failed (sequence file %s):\n%s\n", cfg->sqfp->filename, esl_sqfile_GetErrorBuf(cfg->sqfp));
  else if(sstatus == eslEMEM)    mpi_failure("Out of memory");
  else if(sstatus != eslEOF)     mpi_failure("Unexpected error while reading sequence file");
    
  /* if nec, close tmpfile then merge all alignments in it */
  if(created_tmpfile) { 
    fclose(cfg->tmpfp); /* we're done writing to tmpfp */
    cfg->tmpfp = NULL;
    /* merge all temporary alignments now in cfg->tmpfp, and output merged alignment */
    if((status = create_and_output_final_msa(go, cfg, errbuf, cm, nali, tmpfile)) != eslOK) mpi_failure(errbuf);
    remove(tmpfile); 
  }
    
  /* finish insert and el files */
  if(cfg->ifp != NULL) { fprintf(cfg->ifp, "//\n"); }
  if(cfg->efp != NULL) { fprintf(cfg->efp, "//\n"); }

  /* clean up */
  if(map_sscons != NULL) free(map_sscons);
  FreeCM(cm);
  if(mpibuf != NULL) free(mpibuf);

  return eslOK;
  
  ERROR:
  mpi_failure("Memory allocation error.");
  return status;
}

/* mpi_worker()
 * 
 * Receive sequences from the master, align them and 
 * send CM_ALNDATA results back to master.
 */
static int
mpi_worker(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int             status;                 /* Easel status */
  char            errbuf[eslERRBUFSIZE];  /* for printing error messages */
  CM_t           *cm          = NULL;     /* the CM */
  WORKER_INFO     info;                   /* the worker info */
  ESL_SQ         *sq          = NULL;     /* sequence we're aligning */
  ESL_DSQ        *dsq         = NULL;     /* digitial sequence, rec'd from master */
  CM_ALNDATA     *data        = NULL;     /* data we fill for each seq and send back to master */
  int64_t         L, idx;                 /* sequence length, index, rec'd from master */
  char           *mpibuf      = NULL;     /* buffer used to pack/unpack structures */
  int             mpibuf_size = 0;        /* size of the mpibuf                    */
  int             blocks_remain_in_file;  /* set to FALSE to break outer loop over blocks */
  int             seqs_remain_in_block;   /* set to FALSE to break inner loop over seqs  */
  CM_P7_OM_HOLDER om_holder;              /* reusable --p7pinbridge LOCAL profile/OPROFILE, per MPI worker (brief 26_0430-090) */

  if ((status = init_shared_cfg(go, cfg, errbuf)) != eslOK) mpi_failure(errbuf);
  if(esl_opt_GetBoolean(go, "--sample")) mpi_failure("--sample does not work with in MPI mode (b/c results would not be exactly reproducible)");

  /* Read one CM, and make sure there's only one. This fills cfg->abc. */
  status = cm_file_Read(cfg->cmfp, TRUE, &(cfg->abc), &cm);
  if(status != eslOK) mpi_failure(cfg->cmfp->errbuf);
  status = cm_file_Read(cfg->cmfp, TRUE, &(cfg->abc), NULL);
  if(status != eslEOF) mpi_failure("CM file %s does not contain just one CM\n", cfg->cmfp->fname);

  if((status = initialize_cm(go, cfg, errbuf, cm)) != eslOK) mpi_failure(errbuf);

  /* --hmm mode: HMM-only alignment worker.
   * Receives dsqs from master, computes p7 traces, sends them back.
   */
  if(esl_opt_GetBoolean(go, "--hmm")) {
    P7_HMM     *hmm_w   = cm->mlp7;
    P7_BG      *bg_w    = p7_bg_Create(hmm_w->abc);
    P7_PROFILE *gm_w    = p7_profile_Create(hmm_w->M, hmm_w->abc);
    P7_GMX     *gx_w    = NULL;
    P7_GMX     *gxf_w   = NULL;
    P7_GMX     *gxb_w   = NULL;
    int         p7mode_w = esl_opt_GetBoolean(go, "-g") ? p7_UNIGLOCAL : p7_UNILOCAL;
    int         do_hmmvit_w    = (cm->align_opts & CM_ALIGN_P7HMMVIT)    ? TRUE : FALSE;
    int         do_hmmnoband_w = (cm->align_opts & CM_ALIGN_P7HMMNOBAND) ? TRUE : FALSE;
    int         do_bandedoa_w  = (! do_hmmvit_w && ! do_hmmnoband_w);
    int         do_trunc_w     = (cm->align_opts & CM_ALIGN_TRUNC)       ? TRUE : FALSE; /* brief 26_0430-182 Part A */
    float       sc_w, fwdsc_w, oasc_w;

    /* banded functions declared in cm_p7_band.c */
    extern int p7_kbands2gbands(int *i2k, int *kmin, int *kmax, int L, int M, P7_GBANDS **ret_bnd);
    extern int my_p7_GForwardBanded(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMXB *gxb, float *opt_sc);
    extern int p7_GBackwardBanded(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMXB *gxb, float *opt_sc);
    extern int p7_GDecodingBanded(const P7_PROFILE *gm, const P7_GMXB *fwd, P7_GMXB *bck, P7_GMXB *pp, float overall_sc);
    extern int p7_GOptimalAccuracyBanded(const P7_PROFILE *gm, const P7_GMXB *pp, P7_GMXB *gx, float *ret_e);
    extern int p7_GOATraceBanded(const P7_PROFILE *gm, const P7_GMXB *pp, const P7_GMXB *gx, P7_TRACE *tr);

    /* brief 26_0430-182 Part A: Tgm when do_trunc_w (mirrors hmm_alignment()'s serial-path setup). */
    if (do_trunc_w) {
      p7_ProfileConfig(hmm_w, bg_w, gm_w, 400, p7_LOCAL);
      p7_ProfileConfig5PrimeAnd3PrimeTrunc(gm_w, 400);
    } else {
      p7_ProfileConfig(hmm_w, bg_w, gm_w, 400, p7mode_w);
    }
    if (do_hmmvit_w || do_bandedoa_w) gx_w  = p7_gmx_Create(hmm_w->M, 400);
    if (do_hmmnoband_w)             { gxf_w = p7_gmx_Create(hmm_w->M, 400); gxb_w = p7_gmx_Create(hmm_w->M, 400); }

    /* Signal ready, then enter recv/compute/send loop */
    status = eslOK;
    MPI_Send(&status, 1, MPI_INT, 0, INFERNAL_INITIALREADY_TAG, MPI_COMM_WORLD);

    status = cm_dsq_MPIRecv(0, INFERNAL_DSQ_TAG, MPI_COMM_WORLD, &mpibuf, &mpibuf_size, &dsq, &L, &idx);
    if (status == eslEOD) dsq = NULL; /* termination signal */

    while (dsq != NULL) {
      P7_TRACE *wtr = do_hmmvit_w ? p7_trace_Create() : p7_trace_CreateWithPP();

      /* brief 26_0430-182 Part A: re-run Tgm setup per-seq (p7_ReconfigLength() would
       * clobber the Tgm N->N/C->C loop-disable; see hmm_pipeline_thread comment). */
      if (do_trunc_w) {
	p7_ProfileConfig(hmm_w, bg_w, gm_w, L, p7_LOCAL);
	p7_ProfileConfig5PrimeAnd3PrimeTrunc(gm_w, L);
      } else {
	p7_ReconfigLength(gm_w, L);
      }

      /* preflight: check HMM matrix size vs --mxsize before GrowTo.
       * brief 26_0628-032: skipped under kmeranchor/kmerchain -- the full P7_GMX is
       * only touched on the rare ncells==0 fallback, not on the genome-scale
       * success path (see serial-path comment in hmm_alignment()). */
      if (! cm->p7_use_kmeranchor && ! cm->p7_use_kmerchain) {
	double single_bytes = (double) sizeof(float) * (double)(hmm_w->M + 1) * (double)(L + 1) * (double) p7G_NSCELLS;
	int    nmat         = do_hmmnoband_w ? 2 : 1;
	double needed_mb    = (single_bytes * (double) nmat) / (1024.0 * 1024.0);
	double mxsize_limit = esl_opt_GetReal(go, "--mxsize");
	if (needed_mb > mxsize_limit) {
	  int recommended_mxsize = (int)(ceil(needed_mb / 1024.0) * 1024.0);
	  mpi_failure("HMM-only alignment mx needs %.2f Mb > %.2f Mb limit. Use --mxsize %d.",
		      needed_mb, mxsize_limit, recommended_mxsize);
	}
      }

      if (do_hmmvit_w) {
	p7_gmx_GrowTo(gx_w, hmm_w->M, L);
	p7_GViterbi(dsq, L, gm_w, gx_w, &sc_w);
	p7_GTrace(dsq, L, gm_w, gx_w, wtr);
      }
      else if (do_hmmnoband_w) {
	p7_gmx_GrowTo(gxf_w, hmm_w->M, L);
	p7_gmx_GrowTo(gxb_w, hmm_w->M, L);
	p7_GForward (dsq, L, gm_w, gxf_w, &fwdsc_w);
	p7_GBackward(dsq, L, gm_w, gxb_w, NULL);
	p7_GDecoding(gm_w, gxf_w, gxb_w, gxb_w);
	p7_GOptimalAccuracy(gm_w, gxb_w, gxf_w, &oasc_w);
	p7_GOATrace(gm_w, gxb_w, gxf_w, wtr);
      }
      else {
	/* Viterbi-banded OA */
	int     *i2k_w  = NULL;
	int     *kmin_w = NULL;
	int     *kmax_w = NULL;
	int      ncells_w = 0;
	int      pad_w  = 30;
	P7_GBANDS *bnd_w = NULL;
	P7_GMXB *bxf_w  = NULL;
	P7_GMXB *bxb_w  = NULL;
	P7_TRACE *vtr_w = NULL;
	int       tpos_w;

	/* brief 26_0628-032: k-mer anchor/chain deriver, opt-in via --p7kmeranchor/
	 * --p7kmerchain (mirrors cm_alndata.c:536-566's --p7band dispatch).
	 * brief 26_0628-038: do_trunc now threaded through, mirroring cm_alndata.c's
	 * CM-mode dispatch (cm->align_opts & CM_ALIGN_TRUNC); do_trunc_w already
	 * computed above (brief 26_0430-182 Part A) for the Tgm setup in this same scope. */
	int did_kmer_w = FALSE;
	int *local_nodepad_w = NULL;
	if ((cm->p7_use_kmeranchor || cm->p7_use_kmerchain) && (cm->flags & CMH_P7NODEPAD)) {
	  int k;
	  ESL_ALLOC(local_nodepad_w, sizeof(int) * (hmm_w->M + 1));
	  for (k = 0; k <= hmm_w->M; k++) local_nodepad_w[k] = cm->p7_cm_nodepad[k] + cm->p7bpad;
	}
	if (cm->p7_use_kmeranchor) {
	  did_kmer_w = TRUE;
	  if (p7_Seq2BandsKmerAnchor(cm, errbuf, dsq, L, local_nodepad_w,
				     do_trunc_w, /* brief 26_0628-038: track CM_ALIGN_TRUNC like cm_alndata.c:541 */
				     &i2k_w, &kmin_w, &kmax_w, &ncells_w) != eslOK)
	    mpi_failure("p7_Seq2BandsKmerAnchor() failed: %s", errbuf);
	}
	else if (cm->p7_use_kmerchain) {
	  did_kmer_w = TRUE;
	  if (p7_Seq2BandsKmerChain(cm, errbuf, dsq, L, local_nodepad_w,
				    do_trunc_w, /* brief 26_0628-038: track CM_ALIGN_TRUNC like cm_alndata.c:558 */
				    &i2k_w, &kmin_w, &kmax_w, &ncells_w) != eslOK)
	    mpi_failure("p7_Seq2BandsKmerChain() failed: %s", errbuf);
	}
	if (local_nodepad_w) free(local_nodepad_w);

	if (did_kmer_w && ncells_w == 0 && ! cm->p7_kmerchain_fallback_vit) {
	  /* brief 26_0628-047: M-gate/N-gate fired, or no anchor found -- default
	   * fallback target is --p7ibv's D&C deriver instead of a Vit-trace
	   * band (see serial hmm_alignment()'s matching comment).
	   * brief 26_0628-050: use cm->p7_ibv_delta (struct default 3000), not
	   * --p7ibv-delta's CLI default (20000) -- see serial hmm_alignment(). */
	  int p7ibv_base_slab_w = (esl_opt_IsDefault(go, "--p7ibv-base-slab")
				    ? HMM_P7IBV_KNEE_BASE_SLAB
				    : esl_opt_GetInteger(go, "--p7ibv-base-slab"));
	  if (kmer_gate_p7ibv_fallback(cm, errbuf, dsq, L, do_trunc_w,
					cm->p7_ibv_delta, p7ibv_base_slab_w,
					&i2k_w, &kmin_w, &kmax_w, &ncells_w) != eslOK)
	    mpi_failure("kmer_gate_p7ibv_fallback() failed: %s", errbuf);
	}
	if (! did_kmer_w || ncells_w == 0) {
	  /* Default Viterbi-banded OA path (no kmer flag set); the kmeranchor/
	   * kmerchain ncells==0 fallback when --p7kmerchain-fbvit
	   * reverts to the old behavior; and the safety net when the --p7ibv
	   * fallback above itself also found nothing usable (mirrors
	   * cm_alndata.c's kmeranchor->vitband / kmerchain->vitband fallback
	   * shape). */
	  vtr_w = p7_trace_Create();
	  p7_gmx_GrowTo(gx_w, hmm_w->M, L);
	  p7_GViterbi(dsq, L, gm_w, gx_w, &sc_w);
	  p7_GTrace(dsq, L, gm_w, gx_w, vtr_w);

	  if (vtr_w->N == 0) {
	    /* fallback to unbanded */
	    P7_GMX *fb_gxf = p7_gmx_Create(hmm_w->M, L);
	    P7_GMX *fb_gxb = p7_gmx_Create(hmm_w->M, L);
	    p7_GForward (dsq, L, gm_w, fb_gxf, &fwdsc_w);
	    p7_GBackward(dsq, L, gm_w, fb_gxb, NULL);
	    p7_GDecoding(gm_w, fb_gxf, fb_gxb, fb_gxb);
	    p7_GOptimalAccuracy(gm_w, fb_gxb, fb_gxf, &oasc_w);
	    p7_GOATrace(gm_w, fb_gxb, fb_gxf, wtr);
	    p7_gmx_Destroy(fb_gxf);
	    p7_gmx_Destroy(fb_gxb);
	    p7_trace_Destroy(vtr_w);
	    goto HMM_MPI_SEND;
	  }

	  ESL_ALLOC(i2k_w, sizeof(int) * (L + 1));
	  esl_vec_ISet(i2k_w, (L + 1), -1);
	  for (tpos_w = 0; tpos_w < vtr_w->N; tpos_w++) {
	    if (vtr_w->st[tpos_w] == p7T_M) {
	      int ii = vtr_w->i[tpos_w];
	      int kk = vtr_w->k[tpos_w];
	      if (ii >= 1 && ii <= L && kk >= 1 && kk <= hmm_w->M)
		i2k_w[ii] = kk;
	    }
	  }

	  p7_pins2bands(i2k_w, errbuf, L, hmm_w->M, pad_w, &kmin_w, &kmax_w, &ncells_w);
	}
	p7_kbands2gbands(i2k_w, kmin_w, kmax_w, L, hmm_w->M, &bnd_w);
	bxf_w = p7_gmxb_Create(bnd_w);
	bxb_w = p7_gmxb_Create(bnd_w);

	my_p7_GForwardBanded(dsq, L, gm_w, bxf_w, &fwdsc_w);
	p7_GBackwardBanded(dsq, L, gm_w, bxb_w, NULL);
	p7_GDecodingBanded(gm_w, bxf_w, bxb_w, bxb_w, fwdsc_w);
	p7_GOptimalAccuracyBanded(gm_w, bxb_w, bxf_w, &oasc_w);
	p7_GOATraceBanded(gm_w, bxb_w, bxf_w, wtr);

	free(i2k_w);
	free(kmin_w);
	free(kmax_w);
	p7_trace_Destroy(vtr_w);
	p7_gbands_Destroy(bnd_w);
	p7_gmxb_Destroy(bxf_w);
	p7_gmxb_Destroy(bxb_w);
      }

    HMM_MPI_SEND:
      /* Send trace back to master */
      hmm_trace_MPISend(wtr, idx, 0, INFERNAL_ALNDATA_TAG, MPI_COMM_WORLD, &mpibuf, &mpibuf_size);
      p7_trace_Destroy(wtr);
      free(dsq);

      /* Receive next dsq */
      status = cm_dsq_MPIRecv(0, INFERNAL_DSQ_TAG, MPI_COMM_WORLD, &mpibuf, &mpibuf_size, &dsq, &L, &idx);
      if (status == eslEOD) dsq = NULL;
    }

    /* Receive end-of-file signal (second NULL dsq) */
    status = cm_dsq_MPIRecv(0, INFERNAL_DSQ_TAG, MPI_COMM_WORLD, &mpibuf, &mpibuf_size, &dsq, &L, &idx);

    /* Clean up */
    if (gx_w  != NULL) p7_gmx_Destroy(gx_w);
    if (gxf_w != NULL) p7_gmx_Destroy(gxf_w);
    if (gxb_w != NULL) p7_gmx_Destroy(gxb_w);
    p7_profile_Destroy(gm_w);
    p7_bg_Destroy(bg_w);
    FreeCM(cm);
    if (mpibuf != NULL) free(mpibuf);
    return eslOK;
  }

  /* initialize our worker info */
  info.cm          = cm;
  info.dataA       = NULL;
  info.n           = 0;
  info.mxsize      = esl_opt_GetReal(go, "--mxsize");
  info.pass_idx    = esl_opt_GetBoolean(go, "--notrunc") ? PLI_PASS_STD_ANY : PLI_PASS_5P_AND_3P_FORCE;
  info.w           = esl_stopwatch_Create();
  info.w_tot       = esl_stopwatch_Create();
  info.do_failover = (esl_opt_GetBoolean(go, "--hbanded")  && (! esl_opt_GetBoolean(go, "--notrunc"))) ? TRUE : FALSE;

  /* Main loop: actually two nested while loops, over sequence blocks
   * (while(blocks_remain_in_file)) and over sequences within blocks
   * (while(seqs_remain_in_block)). We exit each loop when the master
   * tells us to, by sending a NULL dsq. If we receive a NULL dsq,
   * we'll exit the inner loop over sequences, and if we immediately
   * receive another NULL dsq we'll exit the outer loop over blocks.
   */
  cm_p7_om_holder_Init(&om_holder);
  blocks_remain_in_file = TRUE;
  while(blocks_remain_in_file) {
    /* inform the master that we're ready for our first seq of the block */
    status = eslOK;
    MPI_Send(&status, 1, MPI_INT, 0, INFERNAL_INITIALREADY_TAG, MPI_COMM_WORLD);

#if DEBUGMPI
    printf("worker %d sent initial ready tag to master, waiting for 1st dsq\n", cfg->my_rank);
#endif 

    /* receive first dsq in block, if it's NULL, we know we're out of blocks */
    status = cm_dsq_MPIRecv(0, INFERNAL_DSQ_TAG, MPI_COMM_WORLD, &mpibuf, &mpibuf_size, &dsq, &L, &idx);
    if     (status == eslOK  && dsq    == NULL)   mpi_failure("problem receiving 1st dsq");
    else if(status == eslEOD && dsq    != NULL)   mpi_failure("problem receiving termination signal");
    else if(status != eslOK  && status != eslEOD) mpi_failure("problem receiving 1st dsq");

    if(dsq == NULL) { /* master is telling us that we're finished with the sequence file */
      blocks_remain_in_file = FALSE;
      seqs_remain_in_block  = FALSE;
#if DEBUGMPI
      printf("worker %d received NULL 1st dsq from master, shutting down\n", cfg->my_rank);
#endif 
    }
    else { /* dsq is valid */
      blocks_remain_in_file = TRUE;
      seqs_remain_in_block  = TRUE;
    }

    while(seqs_remain_in_block) { 
#if DEBUGMPI
      printf("worker %d dsq %" PRId64 " of length %" PRId64 " received from master\n", cfg->my_rank, idx, L);
#endif 
      /* create a sequence object from dsq */
      if ((sq = esl_sq_CreateDigitalFrom(cfg->abc, "irrelevant", dsq, L, NULL, NULL, NULL)) == NULL) mpi_failure("out of memory");
      free(dsq); /* esl_sq_CreateDigitalFrom() makes a copy of dsq */
      
      /* align the sequence */
      status = DispatchSqAlignment(info.cm, errbuf, sq, idx, info.mxsize, TRMODE_UNKNOWN, info.pass_idx, FALSE, /* FALSE: cm->cp9b not valid */
				   info.w, info.w_tot, NULL, &om_holder, &data);
      
      /* If alignment failed: potentially retry alignment in HMM banded
       * std (non-truncated) mode. We will only possibly do this if our
       * initial try was HMM banded truncated alignment (if not,
       * info.do_failover will be FALSE).
       */
      if(status == eslEAMBIGUOUS && info.do_failover == TRUE) { 
	assert(info.cm->align_opts & CM_ALIGN_TRUNC);
	info.cm->align_opts &= ~CM_ALIGN_TRUNC; /* lower truncated alignment flag, just for this sequence */
	status = DispatchSqAlignment(info.cm, errbuf, sq, idx, info.mxsize,
				     TRMODE_UNKNOWN, PLI_PASS_STD_ANY, FALSE, /* USE PLI_PASS_STD_ANY; FALSE: info->cm->cp9b not valid */
				     info.w, info.w_tot, NULL, &om_holder, &data);
	info.cm->align_opts |= CM_ALIGN_TRUNC; /* reraise truncated alignment flag */
      }
      if(status != eslOK) { 
        fprintf(stderr, "Problem during alignment of sequence %s\n", sq->name);
        mpi_failure(errbuf);
      }

      /* pack up the data and send it back to the master (FALSE: don't send data->sq) */
      status = cm_alndata_MPISend(data, FALSE, errbuf, 0, INFERNAL_ALNDATA_TAG, MPI_COMM_WORLD, &mpibuf, &mpibuf_size);
      if(status != eslOK) mpi_failure(errbuf);

      /* clean up old sequence */
      esl_sq_Destroy(sq);
      cm_alndata_Destroy(data, FALSE); /* don't free data->sq, it was pointing at the sq we just free'd */
      
      /* receive next sequence from the master, if it's null that's our signal to stop with this block */
      status = cm_dsq_MPIRecv(0, INFERNAL_DSQ_TAG, MPI_COMM_WORLD, &mpibuf, &mpibuf_size, &dsq, &L, &idx);
      if     (status == eslOK  && dsq    == NULL)   mpi_failure("problem receiving dsq");
      else if(status == eslEOD && dsq    != NULL)   mpi_failure("problem receiving termination signal");
      else if(status != eslOK  && status != eslEOD) mpi_failure("problem receiving dsq");
      
      if(dsq == NULL) seqs_remain_in_block = FALSE; /* we're done with this block */
#if DEBUGMPI
      if(dsq == NULL) printf("worker received NULL dsq from master, block is done");
#endif
      
    } /* end of 'while(seqs_remain_in_block)' */
  } /* end of 'while(blocks_remain_in_file)' */
  cm_p7_om_holder_Reset(&om_holder);

  if(info.cm    != NULL) FreeCM(info.cm);
  if(info.dataA != NULL) free(info.dataA);
  if(info.w     != NULL) esl_stopwatch_Destroy(info.w);
  if(info.w_tot != NULL) esl_stopwatch_Destroy(info.w_tot);
  if(mpibuf     != NULL) free(mpibuf);

  return eslOK;
}

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
#endif /*HAVE_MPI*/

static void
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_cmfile, char **ret_sqfile, int *ret_infmt, int *ret_outfmt)
{
  ESL_GETOPTS *go      = NULL;
  int          infmt   = eslSQFILE_UNKNOWN;
  int          outfmt  = eslMSAFILE_STOCKHOLM;

  if ((go = esl_getopts_Create(options))     == NULL)     cm_Fail("Internal failure creating options object");
  if (esl_opt_ProcessEnvironment(go)         != eslOK)  { printf("Failed to process environment: %s\n", go->errbuf); goto ERROR; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { printf("Failed to parse command line: %s\n", go->errbuf); goto ERROR; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { printf("Failed to parse command line: %s\n", go->errbuf); goto ERROR; }
 
  /* help format: */
  if (esl_opt_GetBoolean(go, "-h")) { 
    cm_banner(stdout, argv[0], banner);
    esl_usage(stdout, argv[0], usage);
    puts("\nBasic options:");
    esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
    puts("\nOptions controlling alignment algorithm:");
    esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
    puts("\nOptions controlling speed and memory requirements:");
    esl_opt_DisplayHelp(stdout, go, 3, 2, 80); 
    puts("\nOptional output files:");
    esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 
    puts("\nOther options:");
    esl_opt_DisplayHelp(stdout, go, 5, 2, 80); 
    puts("\nSequence input formats:   FASTA, GenBank");
    puts("Alignment output formats: Stockholm, Pfam, AFA (aligned FASTA), A2M, Clustal, PHYLIP\n");
    exit(0);
  } 

  if (esl_opt_ArgNumber(go)                 != 2)     { puts("Incorrect number of command line arguments.");      goto ERROR; }
  if ((*ret_cmfile = esl_opt_GetArg(go, 1)) == NULL)  { puts("Failed to get <cmfile> argument on command line");  goto ERROR; }
  if ((*ret_sqfile = esl_opt_GetArg(go, 2)) == NULL)  { puts("Failed to get <seqfile> argument on command line"); goto ERROR; }

  if (strcmp(*ret_cmfile, "-") == 0 && strcmp(*ret_sqfile, "-") == 0) { 
    puts("\nERROR: Either <cmfile> or <seqfile> may be '-' (to read from stdin), but not both.\n");
    goto ERROR;
  }

  /* If caller declared an input format, decode it */
  if (esl_opt_IsOn(go, "--informat")) {
    infmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--informat"));
    if (infmt == eslSQFILE_UNKNOWN) { 
      printf("\nERROR: %s is not a recognized input sequence file format\n\n", esl_opt_GetString(go, "--informat"));
      goto ERROR;
    }
  }

  /* Determine output alignment file format */
  outfmt = esl_msafile_EncodeFormat(esl_opt_GetString(go, "--outformat"));
  if (outfmt == eslMSAFILE_UNKNOWN) {
    printf("\nERROR: %s is not a recognized output MSA file format\n\n", esl_opt_GetString(go, "--outformat"));
    goto ERROR;
  }

  /* Check for incompatible option combinations too complex for esl_getopts to enforce during declaration */

#ifdef HMMER_THREADS
  /* if --sample, enforce that --cpu 0 is used if HMMER_THREADS, otherwise number of threads would
   * affect the sampled alignments (each thread requires its own RNG) 
   */
  if (esl_opt_GetBoolean(go, "--sample")) { 
    if((! esl_opt_IsUsed(go, "--cpu")) || 
       (  esl_opt_IsUsed(go, "--cpu") && (esl_opt_GetInteger(go, "--cpu") != 0))) { 
      puts("\nERROR: --sample requires --cpu 0\n");
      goto ERROR;
    }
  }
#endif /* HMMER_THREADS */
#ifdef HAVE_MPI
  /* --sample is incompatible with --mpi b/c sampled parsetrees would be dependent
   * on number of workers (each of which needs its own (separately seeded) RNG)
   */
  if (esl_opt_GetBoolean(go, "--sample") && esl_opt_IsUsed(go, "--mpi")) {
    puts("\nERROR: --sample is incompatible with --mpi\n");
    goto ERROR;
  }	
#endif /* HAVE_MPI */  

  /* --verbose only makes sense in combination with -o or --sfile, 
   * because if neither is used, scores are not output.
   */
  if (esl_opt_GetBoolean(go, "--verbose") && (! esl_opt_IsUsed(go, "-o")) && (! esl_opt_IsUsed(go, "--sfile"))) {
    puts("\nERROR: --verbose only makes sense in combination with -o or --sfile\n");
    goto ERROR;
  }	

  /* Finally, check for incompatible option combinations that
   * esl_getopts can handle, but that would require an error message
   * like: "Option 'x' is incompatible with options
   * y1,y2,y3,y4....yn", where there's so many y's that the message is
   * truncated because errbuf runs out of space. As a workaround we
   * laboriously check for all incompatible options of that type here.
   */
  if(esl_opt_IsUsed(go, "--small")) {
    if((! esl_opt_IsUsed(go, "--cyk")) || (! esl_opt_IsUsed(go, "--noprob")) || (! esl_opt_IsUsed(go, "--nonbanded")) || (! esl_opt_IsUsed(go, "--notrunc"))) {
      puts("Failed to parse command line: Option --small requires --cyk, --noprob, --nonbanded, --notrunc");
      goto ERROR;
    }
  }

  /* --p7ibv only derives bands; it needs an anchor mode to use them:
   * --p7band (CM-side banded alignment) or --hmm (HMM-only banded OA).
   */
  if(esl_opt_GetBoolean(go, "--p7ibv") && (! esl_opt_GetBoolean(go, "--p7band")) && (! esl_opt_GetBoolean(go, "--hmm"))) {
    puts("\nERROR: --p7ibv requires --p7band or --hmm\n");
    goto ERROR;
  }
  /* --hmm --p7ibv is the banded-OA HMM sub-mode; reject the other --hmm
   * sub-modes (Viterbi-trace and unbanded full OA) in combination with it.
   */
  if(esl_opt_GetBoolean(go, "--hmm") && esl_opt_GetBoolean(go, "--p7ibv")) {
    if(esl_opt_GetBoolean(go, "--hmmvit")) {
      puts("\nERROR: --hmmvit incompatible with --p7ibv (pins-to-trace not implemented yet)\n");
      goto ERROR;
    }
    if(esl_opt_GetBoolean(go, "--hmmnoband")) {
      puts("\nERROR: --hmmnoband incompatible with --p7ibv (--hmmnoband means no bands)\n");
      goto ERROR;
    }
  }

  /* brief 26_0628-032: --p7kmeranchor/--p7kmerchain only derive bands; they need an
   * anchor mode to use them: --p7band (CM-side) or --hmm (HMM-only banded OA),
   * same requirement as --p7ibv above.
   */
  if((esl_opt_GetBoolean(go, "--p7kmeranchor") || esl_opt_GetBoolean(go, "--p7kmerchain")) &&
     (! esl_opt_GetBoolean(go, "--p7band")) && (! esl_opt_GetBoolean(go, "--hmm"))) {
    puts("\nERROR: --p7kmeranchor/--p7kmerchain require --p7band or --hmm\n");
    goto ERROR;
  }
  /* brief 26_0628-038: --p7kmeranchor/--p7kmerchain no longer require --notrunc --
   * do_trunc is now threaded through both the --p7band (brief 26_0628-033) and
   * --hmm (brief 26_0628-038) call sites, mirroring cm_alndata.c's CM-mode dispatch. */
  /* --hmm --p7kmeranchor/--p7kmerchain is the k-mer-banded-OA HMM sub-mode;
   * reject the other --hmm sub-modes (Viterbi-trace and unbanded full OA) in
   * combination with it, mirroring the --p7ibv incompatibility above.
   */
  if(esl_opt_GetBoolean(go, "--hmm") && (esl_opt_GetBoolean(go, "--p7kmeranchor") || esl_opt_GetBoolean(go, "--p7kmerchain"))) {
    if(esl_opt_GetBoolean(go, "--hmmvit")) {
      puts("\nERROR: --hmmvit incompatible with --p7kmeranchor/--p7kmerchain (pins-to-trace not implemented yet)\n");
      goto ERROR;
    }
    if(esl_opt_GetBoolean(go, "--hmmnoband")) {
      puts("\nERROR: --hmmnoband incompatible with --p7kmeranchor/--p7kmerchain (--hmmnoband means no bands)\n");
      goto ERROR;
    }
  }
  /* brief 26_0628-046: --p7kmerchain-mink only means something if kmeranchor/kmerchain
   * is actually in use; not expressible as an esl_getopts "reqs" OR (either
   * flag suffices), so check manually, mirroring the --p7kmeranchor/--p7kmerchain
   * "requires --p7band or --hmm" check above. */
  if(esl_opt_IsOn(go, "--p7kmerchain-mink") && esl_opt_GetInteger(go, "--p7kmerchain-mink") > 0 &&
     (! esl_opt_GetBoolean(go, "--p7kmeranchor")) && (! esl_opt_GetBoolean(go, "--p7kmerchain"))) {
    puts("\nERROR: --p7kmerchain-mink requires --p7kmeranchor or --p7kmerchain\n");
    goto ERROR;
  }
  /* brief 26_0628-047: --p7kmerchain-mgate/-fallback-vit only mean something if
   * kmeranchor/kmerchain is actually in use; same manual-check shape as
   * --p7kmerchain-mink above (not expressible as esl_getopts "reqs" OR). */
  if(esl_opt_IsOn(go, "--p7kmerchain-mgate") && esl_opt_GetInteger(go, "--p7kmerchain-mgate") > 0 &&
     (! esl_opt_GetBoolean(go, "--p7kmeranchor")) && (! esl_opt_GetBoolean(go, "--p7kmerchain"))) {
    puts("\nERROR: --p7kmerchain-mgate requires --p7kmeranchor or --p7kmerchain\n");
    goto ERROR;
  }
  if(esl_opt_GetBoolean(go, "--p7kmerchain-fbvit") &&
     (! esl_opt_GetBoolean(go, "--p7kmeranchor")) && (! esl_opt_GetBoolean(go, "--p7kmerchain"))) {
    puts("\nERROR: --p7kmerchain-fbvit requires --p7kmeranchor or --p7kmerchain\n");
    goto ERROR;
  }

  *ret_go     = go;
  *ret_infmt  = infmt;
  *ret_outfmt = outfmt;

  return;
  
 ERROR:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  puts("\nwhere basic options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
  printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
  exit(1);  
}

/* output_header(): 
 *
 * In contrast with other Infernal applications, which output header
 * to stdout, we output the header to stdout only if the user has
 * specified a non-stdout output file for the alignment. Otherwise,
 * the alignment will be printed to stdout without a header because
 * we want it to be a valid Stockholm format (or other) alignment.
 */

static int
output_header(FILE *ofp, const ESL_GETOPTS *go, char *cmfile, char *sqfile, CM_t *cm, int ncpus)
{
  cm_banner(ofp, go->argv[0], banner);
                                            fprintf(ofp, "# CM file:                                     %s\n", cmfile);
			                    fprintf(ofp, "# sequence file:                               %s\n", sqfile);
                                            fprintf(ofp, "# CM name:                                     %s\n", cm->name);
  if (esl_opt_IsUsed(go, "-o"))          {  fprintf(ofp, "# saving alignment to file:                    %s\n", esl_opt_GetString(go, "-o")); }
  if (esl_opt_IsUsed(go, "-g"))          {  fprintf(ofp, "# model configuration:                         global\n"); }

  if (esl_opt_IsUsed(go, "--optacc"))    {  fprintf(ofp, "# alignment algorithm:                         optimal accuracy\n"); }
  if (esl_opt_IsUsed(go, "--cyk"))       {  fprintf(ofp, "# alignment algorithm:                         CYK\n"); }
  if (esl_opt_IsUsed(go, "--sample"))    {  fprintf(ofp, "# sampling aln from posterior distribution:    yes\n"); }
  if (esl_opt_IsUsed(go, "--seed"))      {
    if (esl_opt_GetInteger(go, "--seed") == 0) fprintf(ofp, "# random number seed:                          one-time arbitrary\n");
    else                                       fprintf(ofp, "# random number seed set to:                   %d\n", esl_opt_GetInteger(go, "--seed"));
  }
  if (esl_opt_IsUsed(go, "--notrunc"))   {  fprintf(ofp, "# truncated sequence alignment mode:           off\n"); }
  if (esl_opt_IsUsed(go, "--sub"))       {  fprintf(ofp, "# alternative truncated seq alignment mode:    on\n"); }
  if (esl_opt_IsUsed(go, "--hmm"))       {  fprintf(ofp, "# alignment method:                            p7 HMM only (no CM)\n"); }
  if (esl_opt_IsUsed(go, "--hmmvit"))    {  fprintf(ofp, "# HMM alignment algorithm:                     Viterbi\n"); }
  if (esl_opt_IsUsed(go, "--hmmnoband")) {  fprintf(ofp, "# HMM alignment banding:                       off (full OA)\n"); }

  if (esl_opt_IsUsed(go, "--mxsize"))    {  fprintf(ofp, "# maximum total DP matrix size set to:         %.2f Mb\n", esl_opt_GetReal(go, "--mxsize")); }
  if (esl_opt_IsUsed(go, "--hbanded"))   {  fprintf(ofp, "# using HMM bands for acceleration:            yes\n"); }
  if (esl_opt_IsUsed(go, "--tau"))       {  fprintf(ofp, "# tail loss probability for HMM bands set to:  %g\n", esl_opt_GetReal(go, "--tau")); }
  if (esl_opt_IsUsed(go, "--fixedtau"))  {  fprintf(ofp, "# tighten HMM bands when necessary:            no\n"); }
  if (esl_opt_IsUsed(go, "--maxtau"))    {  fprintf(ofp, "# maximum tau allowed during band tightening:  %g\n", esl_opt_GetReal(go, "--maxtau")); }
  if (esl_opt_IsUsed(go, "--nonbanded")) {  fprintf(ofp, "# using HMM bands for acceleration:            no\n"); }
  if (esl_opt_IsUsed(go, "--small"))     {  fprintf(ofp, "# small memory D&C alignment algorithm:        on\n"); }

  if (esl_opt_IsUsed(go, "--sfile"))     {  fprintf(ofp, "# saving alignment score info to file:         %s\n", esl_opt_GetString(go, "--sfile")); }
  if (esl_opt_IsUsed(go, "--tfile"))     {  fprintf(ofp, "# saving parsetrees to file:                   %s\n", esl_opt_GetString(go, "--tfile")); }
  if (esl_opt_IsUsed(go, "--ifile"))     {  fprintf(ofp, "# saving insert information to file:           %s\n", esl_opt_GetString(go, "--ifile")); }
  if (esl_opt_IsUsed(go, "--elfile"))    {  fprintf(ofp, "# saving local end information to file:        %s\n", esl_opt_GetString(go, "--elfile")); }

  if (esl_opt_IsUsed(go, "--mapali"))    {  fprintf(ofp, "# including alignment from file:               %s\n", esl_opt_GetString(go, "--mapali")); }
  if (esl_opt_IsUsed(go, "--mapstr"))    {  fprintf(ofp, "# including structure from alnment from file:  %s\n", esl_opt_GetString(go, "--mapali")); }
  if (esl_opt_IsUsed(go, "--informat"))  {  fprintf(ofp, "# input sequence file format specified as:     %s\n", esl_opt_GetString(go, "--informat")); }
  if (esl_opt_IsUsed(go, "--outformat")) {  fprintf(ofp, "# output alignment format specified as:        %s\n", esl_opt_GetString(go, "--outformat")); }
  if (esl_opt_IsUsed(go, "--dnaout"))    {  fprintf(ofp, "# output alignment alphabet:                   DNA\n"); }
  if (esl_opt_IsUsed(go, "--noprob"))    {  fprintf(ofp, "# posterior probability annotation:            off\n"); }
  if (esl_opt_IsUsed(go, "--matchonly")) {  fprintf(ofp, "# include alignment insert columns:            no\n"); }
  if (esl_opt_IsUsed(go, "--ileaved"))   {  fprintf(ofp, "# forcing interleaved Stockholm output aln:    yes\n"); }
  if (esl_opt_IsUsed(go, "--regress"))   {  fprintf(ofp, "# saving alignment without author info to:     %s\n", esl_opt_GetString(go, "--regress")); }

  /* output number of processors being used, always (this differs from H3 which only does this if --cpu) */
  int output_ncpu = FALSE;
#ifdef HAVE_MPI
  if (esl_opt_IsUsed(go, "--mpi"))       {  fprintf(ofp, "# MPI:                                         on [%d processors]\n", ncpus); output_ncpu = TRUE; }
#endif 
#ifdef HMMER_THREADS
  if (! output_ncpu)                     {  fprintf(ofp, "# number of worker threads:                    %d%s\n", ncpus, (esl_opt_IsUsed(go, "--cpu") ? " [--cpu]" : "")); output_ncpu = TRUE; }
#endif 
  if (! output_ncpu)                     {  fprintf(ofp, "# number of worker threads:                    0 [serial mode; threading unavailable]\n"); }
  fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");

  return eslOK;
}

/* init_master_cfg()
 * Called by masters, mpi or serial.
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
  
  /* initialize cfg variables used by masters and workers */
  if((status = init_shared_cfg(go, cfg, errbuf)) != eslOK) return status;

  /* open output files */
  cfg->ofp = stdout;
  if (esl_opt_IsUsed(go, "-o")) { 
    if ((cfg->ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL)      ESL_FAIL(eslFAIL, errbuf, "Failed to open -o output file %s\n", esl_opt_GetString(go, "-o"));
  } 
  if (esl_opt_IsUsed(go, "--tfile")) { 
    if ((cfg->tfp = fopen(esl_opt_GetString(go, "--tfile"), "w")) == NULL) ESL_FAIL(eslFAIL, errbuf, "Failed to open --tfile output file %s\n", esl_opt_GetString(go, "--tfile"));
  }
  if (esl_opt_IsUsed(go, "--ifile")) { 
    if ((cfg->ifp = fopen(esl_opt_GetString(go, "--ifile"), "w")) == NULL) ESL_FAIL(eslFAIL, errbuf, "Failed to open --ifile output file %s\n", esl_opt_GetString(go, "--ifile"));
    output_info_file_header(cfg->ifp, "Insert information file created by cmalign.", "");
  }
  if (esl_opt_IsUsed(go, "--elfile")) { 
    if ((cfg->efp = fopen(esl_opt_GetString(go, "--elfile"), "w")) == NULL) ESL_FAIL(eslFAIL, errbuf, "Failed to open --elfile output file %s\n", esl_opt_GetString(go, "--elfile"));
    output_info_file_header(cfg->efp, "EL state (local end) insert information file created by cmalign.", "EL ");
  }
  if (esl_opt_IsUsed(go, "--sfile")) { 
    if ((cfg->sfp = fopen(esl_opt_GetString(go, "--sfile"), "w")) == NULL) ESL_FAIL(eslFAIL, errbuf, "Failed to open --sfile output file %s\n", esl_opt_GetString(go, "--sfile"));
  }
  if (esl_opt_IsUsed(go, "--regress")) { 
    if ((cfg->rfp = fopen(esl_opt_GetString(go, "--regress"), "w")) == NULL) ESL_FAIL(eslFAIL, errbuf, "Failed to open --regress output file %s\n", esl_opt_GetString(go, "--regress"));
  }
  return eslOK;
}

/* init_shared_cfg()
 * Called by serial masters and mpi workers and masters.
 *
 * Errors in the MPI master here are considered to be "recoverable",
 * in the sense that we'll try to delay output of the error message
 * until we've cleanly shut down the worker processes. Therefore
 * errors return (code, errbuf) by the ESL_FAIL mech.
 */
static int
init_shared_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  int status;

  /* open CM file */
  status = cm_file_Open(cfg->cmfile, NULL, FALSE, &(cfg->cmfp), errbuf);
  if      (status == eslENOTFOUND) return status;
  else if (status == eslEFORMAT)   return status;
  else if (status != eslOK)        return status;

  /* open sequence file */
  status = esl_sqfile_Open(cfg->sqfile, cfg->infmt, p7_SEQDBENV, &(cfg->sqfp));
  if      (status == eslENOTFOUND) ESL_FAIL(status, errbuf, "Failed to open sequence file %s for reading\n",          cfg->sqfile);
  else if (status == eslEFORMAT)   ESL_FAIL(status, errbuf, "Sequence file %s is empty or misformatted\n",            cfg->sqfile);
  else if (status == eslEINVAL && cfg->infmt == eslSQFILE_UNKNOWN) ESL_FAIL(status, errbuf, "Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        ESL_FAIL(status, errbuf, "Unexpected error %d opening sequence file %s\n", status, cfg->sqfile);  

  cfg->be_verbose = esl_opt_GetBoolean(go, "--verbose");

  return eslOK;
}

/* initialize_cm()
 * Setup the CM based on the command-line options/defaults;
 * set flags and a few parameters. cm_Configure configures
 * the CM.
 */
static int
initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int status;

  /* set up alignment options in cm->align_opts */
  if     (  esl_opt_GetBoolean(go, "--cyk"))    cm->align_opts |= CM_ALIGN_CYK;
  else if(  esl_opt_GetBoolean(go, "--sample")) cm->align_opts |= CM_ALIGN_SAMPLE;
  else                                          cm->align_opts |= CM_ALIGN_OPTACC;
  if(  esl_opt_GetBoolean(go, "--hbanded"))     cm->align_opts |= CM_ALIGN_HBANDED;
  if(  esl_opt_GetBoolean(go, "--nonbanded"))   cm->align_opts |= CM_ALIGN_NONBANDED;
  if(  esl_opt_GetBoolean(go, "--p7band"))    { cm->align_opts |= CM_ALIGN_HBANDED; cm->align_opts |= CM_ALIGN_P7BANDED; }
  if(! esl_opt_GetBoolean(go, "--noprob"))      cm->align_opts |= CM_ALIGN_POST;
  if(! esl_opt_GetBoolean(go, "--notrunc"))     cm->align_opts |= CM_ALIGN_TRUNC;
  if(  esl_opt_GetBoolean(go, "--sub"))         cm->align_opts |= CM_ALIGN_SUB;   /* --sub requires --notrunc */
  if(  esl_opt_GetBoolean(go, "--small"))       cm->align_opts |= CM_ALIGN_SMALL; /* --small requires --noprob --nonbanded --cyk */
  if(  esl_opt_GetBoolean(go, "--hmm"))         cm->align_opts |= CM_ALIGN_P7HMM;
  if(  esl_opt_GetBoolean(go, "--hmmvit"))       cm->align_opts |= CM_ALIGN_P7HMMVIT;
  if(  esl_opt_GetBoolean(go, "--hmmnoband"))    cm->align_opts |= CM_ALIGN_P7HMMNOBAND;
  if(  esl_opt_GetBoolean(go, "--ckpt"))        cm->align_opts |= CM_ALIGN_CHECKPT; /* --ckpt requires -g; sqrt(M)-mem optacc in both truncated (default) and --notrunc modes */
  if((! esl_opt_GetBoolean(go, "--fixedtau")) &&
     (  esl_opt_GetBoolean(go, "--hbanded"))) { 
    cm->align_opts |= CM_ALIGN_XTAU;
  }

  /* set up configuration options in cm->config_opts */
  if(  esl_opt_GetBoolean(go, "--nonbanded"))   cm->config_opts |= CM_CONFIG_NONBANDEDMX;
  if(! esl_opt_GetBoolean(go, "--notrunc"))     cm->config_opts |= CM_CONFIG_TRUNC;
  if(  esl_opt_GetBoolean(go, "--sub"))         cm->config_opts |= CM_CONFIG_SUB;   /* --sub requires --notrunc */
  if(! esl_opt_GetBoolean(go, "-g")) { 
    cm->config_opts |= CM_CONFIG_LOCAL;
    cm->config_opts |= CM_CONFIG_HMMLOCAL;
    cm->config_opts |= CM_CONFIG_HMMEL;
  }
  
  cm->tau    = esl_opt_GetReal(go, "--tau");
  cm->maxtau = esl_opt_GetReal(go, "--maxtau");
  if(esl_opt_GetBoolean(go, "--p7band")) cm->p7bpad = esl_opt_GetInteger(go, "--p7padplus");
  if(esl_opt_GetBoolean(go, "--p7pinbridge")) {
    cm->p7_use_pinbridge = TRUE;
    cm->p7_pinbridge_pad = esl_opt_GetInteger(go, "--p7pbpad");
    if(esl_opt_GetBoolean(go, "--p7pinbridge-vitgaps")) cm->p7_pinbridge_vit_gaps = TRUE;
  }
  if(esl_opt_GetBoolean(go, "--p7ibv")) {
    cm->p7_use_ibv   = TRUE;
    cm->p7_ibv_delta = esl_opt_GetInteger(go, "--p7ibv-delta");
    cm->p7_ibv_width = esl_opt_GetInteger(go, "--p7ibv-width");  /* brief 26_0430-140 */
    {                                                            /* brief 26_0430-140: parse --p7ibv-mode */
      const char *ibvmode = esl_opt_GetString(go, "--p7ibv-mode");
      if      (strcmp(ibvmode, "delta")  == 0) cm->p7_ibv_mode = P7IBV_MODE_DELTA;
      else if (strcmp(ibvmode, "fixed")  == 0) cm->p7_ibv_mode = P7IBV_MODE_FIXED;
      else if (strcmp(ibvmode, "hybrid") == 0) cm->p7_ibv_mode = P7IBV_MODE_HYBRID;
      else cm_Fail("--p7ibv-mode must be one of: delta, fixed, hybrid (got '%s')", ibvmode);
    }
    if(esl_opt_GetBoolean(go, "--p7ibv-mem")) {
      cm->p7_ibv_mem       = TRUE;
      cm->p7_ibv_base_slab = esl_opt_GetInteger(go, "--p7ibv-base-slab");
      if(esl_opt_GetBoolean(go, "--p7ibv-ckpt")) cm->p7_ibv_ckpt = TRUE;
    }
    if(esl_opt_GetBoolean(go, "--p7ibv-wv")) cm->p7_ibv_wv = TRUE;  /* brief 26_0430-169 */
  }
  if(esl_opt_GetBoolean(go, "--p7kmeranchor")) cm->p7_use_kmeranchor = TRUE;  /* brief 26_0628-026 */
  if(esl_opt_GetBoolean(go, "--p7kmerchain"))  cm->p7_use_kmerchain  = TRUE;  /* brief 26_0628-027 */
  cm->p7_kmerchain_ramp_alpha = esl_opt_GetReal(go, "--p7kmerchain-alpha");  /* brief 26_0628-043; req="--p7kmerchain" so only meaningful there */
  cm->p7_kmerchain_mink = esl_opt_GetInteger(go, "--p7kmerchain-mink");     /* brief 26_0628-046; 0 = disabled (default) */
  cm->p7_kmerchain_mgate = esl_opt_GetInteger(go, "--p7kmerchain-mgate");   /* brief 26_0628-047; 0 = disabled (default) */
  cm->p7_kmerchain_fallback_vit = esl_opt_GetBoolean(go, "--p7kmerchain-fbvit"); /* brief 26_0628-047; default FALSE (--p7ibv fallback) */
  if(esl_opt_GetBoolean(go, "--cykbands")) {
    cm->p7_use_cykbands = TRUE;
    cm->p7_cykbands_pad = esl_opt_GetInteger(go, "--cykpad");
    cm->p7_cykskip_unvisited = esl_opt_GetBoolean(go, "--cykskip-unvisited");
  }
  if(esl_opt_IsUsed(go, "--dump-bands")) {
    cm->p7_dump_bands_file = (char *) esl_opt_GetString(go, "--dump-bands");
  }

  if((esl_opt_IsUsed(go, "--flanktoins")) && (esl_opt_IsUsed(go, "--flankselfins"))) { 
    configure_root_inserts(cm, esl_opt_GetReal(go, "--flanktoins"), esl_opt_GetReal(go, "--flankselfins"));
  }
  
  /* configure */
  if((status = cm_Configure(cm, errbuf, -1)) != eslOK) return status;

  /* Brief 26_0430-169: with --p7ibv-wv, calibrate the windowed-Viterbi per-node pad
   * (F+B-halfwidth quantile) ONCE per CM here -- single-threaded, after
   * cm_Configure populated cm->fp7 and before any worker threads spawn -- and
   * cache it on the CM (workers read it read-only).  This is the align-time
   * calibration: works on existing CMs (only needs cm->fp7), no rebuild. */
  if(cm->p7_ibv_wv) {
    int wk;
    if(cm->fp7 == NULL) ESL_FAIL(eslEINVAL, errbuf, "--p7ibv-wv requires cm->fp7 (ML p7 filter)");
    /* Brief 26_0430-172 Phase B (pad amortization): the genome WV pad calibration runs
     * nsamp full-length deriver passes (prohibitive at genome). --p7wvpad-file
     * loads a once-computed pad (skip per-run calib); --p7wvpad-dump writes the
     * calibrated pad for reuse.  This is the pad-storage mechanism the brief
     * requires; serializing it onto the CM file (tag P7WVPAD, cmbuild --p7wv-q)
     * is the production form and is a mechanical follow-up (see summary).  The
     * pad is per-consensus-column [0..fp7->M] (fp7->M == clen). */
    if(esl_opt_IsOn(go, "--p7wvpad-file")) {
      FILE *pf = fopen(esl_opt_GetString(go, "--p7wvpad-file"), "r");
      int   padM = 0, kk, vv;
      char  line[256];
      if(pf == NULL) ESL_FAIL(eslFAIL, errbuf, "failed to open --p7wvpad-file %s", esl_opt_GetString(go, "--p7wvpad-file"));
      cm->p7_wv_nodepad = malloc(sizeof(int) * (cm->fp7->M + 1));
      if(cm->p7_wv_nodepad == NULL) { fclose(pf); ESL_FAIL(eslEMEM, errbuf, "malloc failed for --p7wvpad-file"); }
      for(wk = 0; wk <= cm->fp7->M; wk++) cm->p7_wv_nodepad[wk] = 0;
      while(fgets(line, sizeof(line), pf) != NULL) {
        if(line[0] == '#') continue;
        if(sscanf(line, "%d %d", &kk, &vv) == 2 && kk >= 0 && kk <= cm->fp7->M) { cm->p7_wv_nodepad[kk] = vv; if(kk > padM) padM = kk; }
      }
      fclose(pf);
      if(padM != cm->fp7->M) ESL_FAIL(eslEINCOMPAT, errbuf, "--p7wvpad-file max index %d != fp7->M %d", padM, cm->fp7->M);
      cm->p7_wv_nodepad_M = cm->fp7->M;
    } else if(! esl_opt_GetBoolean(go, "--p7wv-calib")) {
      /* Brief 26_0430-173 Part A (DEFAULT): a constant band half-width of 30 ties the
       * per-node calibrated p95 pad in aggregate (brief 26_0430-174), so the default WV
       * path skips Monte-Carlo calibration entirely -- the post-172 genome
       * dominator (~29-50 min cm_ComputeP7WVNodePad) vanishes.  Fill every node
       * with --p7wv-pad's value (index 0 = 0, matching the calibrator).  The
       * per-node calibration machinery is preserved (opt-in via --p7wv-calib /
       * --p7wvpad-file) for the later tighter-band optimization phase. */
      int padval = esl_opt_GetInteger(go, "--p7wv-pad");
      cm->p7_wv_nodepad = malloc(sizeof(int) * (cm->fp7->M + 1));
      if(cm->p7_wv_nodepad == NULL) ESL_FAIL(eslEMEM, errbuf, "malloc failed for --p7wv-pad constant pad");
      cm->p7_wv_nodepad[0] = 0;
      for(wk = 1; wk <= cm->fp7->M; wk++) cm->p7_wv_nodepad[wk] = padval;
      cm->p7_wv_nodepad_M = cm->fp7->M;
    } else {
      ESL_RANDOMNESS *wv_r = esl_randomness_Create((uint32_t) esl_opt_GetInteger(go, "--p7wv-seed"));
      if(wv_r == NULL) ESL_FAIL(eslEMEM, errbuf, "failed to allocate RNG for --p7ibv-wv pad calibration");
      status = cm_ComputeP7WVNodePad(cm, errbuf, wv_r,
                                     esl_opt_GetInteger(go, "--p7wv-nsamp"),
                                     esl_opt_GetReal(go,    "--p7wv-q"),
                                     esl_opt_GetInteger(go, "--p7ibv-delta"),
                                     esl_opt_GetInteger(go, "--p7wv-floor"),
                                     &(cm->p7_wv_nodepad));
      esl_randomness_Destroy(wv_r);
      if(status != eslOK) return status;
      cm->p7_wv_nodepad_M = cm->fp7->M;
    }
    if(esl_opt_IsOn(go, "--p7wvpad-dump")) {
      FILE *df = fopen(esl_opt_GetString(go, "--p7wvpad-dump"), "w");
      if(df == NULL) ESL_FAIL(eslFAIL, errbuf, "failed to open --p7wvpad-dump %s", esl_opt_GetString(go, "--p7wvpad-dump"));
      fprintf(df, "# brief172 WV per-node pad  M=%d  (k pad)\n", cm->fp7->M);
      for(wk = 0; wk <= cm->fp7->M; wk++) fprintf(df, "%d %d\n", wk, cm->p7_wv_nodepad[wk]);
      fclose(df);
    }
  }

  return eslOK;
}


/* map_alignment()
 *                   
 * Called if the --mapali <f> option is used. Open and read a 
 * MSA from <f>, confirm it was the same alignment used to 
 * build the CM, and convert its aligned sequences to 
 * parsetrees. Return data (sequences and parsetrees) gets 
 * populated into <ret_dataA>.
 * 
 * Also, return the dealigned (cm-> clen length) SS_cons
 * from the alignment in <ret_ss>. This will be used to 
 * overwrite the output alignment's SS_cons if --mapstr 
 * was used. 
 */
static int
map_alignment(const char *msafile, CM_t *cm, int noss_used, char *errbuf, CM_ALNDATA ***ret_dataA, int *ret_ndata, char **ret_ss)
{
  int            status;
  ESL_MSAFILE   *afp       = NULL;
  ESL_MSA       *msa       = NULL;
  ESL_ALPHABET  *abc       = (ESL_ALPHABET *) cm->abc; /* removing const'ness to make compiler happy. Safe. */
  uint32_t       chksum    = 0;
  int            i, x;              /* counters */
  int            apos, uapos, cpos; /* counter over aligned, unaligned, consensus positions */
  int           *a2u_map   = NULL;  /* map from aligned to unaligned positions */
  Parsetree_t   *mtr       = NULL;  /* the guide tree for mapali */
  CM_ALNDATA   **dataA     = NULL;  /* includes ptrs to sq and parsetrees */
  char          *aseq      = NULL;  /* aligned sequence, req'd by Transmogrify() */
  char          *ss        = NULL;  /* msa's SS_cons, if there is one, dealigned to length cm->clen */
  int           *used_el   = NULL;  /* [1..msa->alen] used_el[apos] = TRUE if apos is modeled by EL state, else FALSE */
  /* variables used for copying the structure and possibly removing broken basepairs (for which exactly 1 of the 2 paired positions is a consensus column */
  int            opos;              /* position that apos pairs with */
  int           *i_am_rf   = NULL;  /* [1..msa->alen] i_am_rf[apos] = 1 if alignment position apos is a consensus (RF) position, else 0 */
  int           *msa_ct    = NULL;  /* [1..msa->alen] msa_ct[apos] = x; x==0 if apos is unpaired, x==opos if apos is paired to opos in msa->ss_cons */
  char          *msa_ss_cons_copy = NULL; /* copy of the msa's SS_cons we remove broken basepair halves from */
  
  status = esl_msafile_Open(&abc, msafile, NULL, eslMSAFILE_UNKNOWN, NULL, &afp);
  if (status != eslOK) esl_msafile_OpenFailure(afp, status);

  status = esl_msafile_Read(afp, &msa);
  if (status != eslOK) esl_msafile_ReadFailure(afp, status);

  if (! (cm->flags & CMH_CHKSUM))  cm_Fail("CM has no checksum. --mapali unreliable without it.");
  if (! (cm->flags & CMH_MAP))     cm_Fail("CM has no map. --mapali can't work without it.");
  esl_msa_Checksum(msa, &chksum);
  if (cm->checksum != chksum)      cm_Fail("--mapali MSA %s isn't same as the one CM came from (checksum mismatch)", msafile);

  /* allocate and initialize dataA */
  ESL_ALLOC(dataA, sizeof(CM_ALNDATA *) * msa->nseq);
  for(i = 0; i < msa->nseq; i++) dataA[i] = cm_alndata_Create();

  /* allocated msa_ct */
  ESL_ALLOC(msa_ct, sizeof(int) * (msa->alen+1));

  /* if --noss used, potentially remove ss_cons and replace with no basepairs */
  if(noss_used) { 
    if(msa->ss_cons != NULL) { free(msa->ss_cons); msa->ss_cons = NULL; }
    ESL_ALLOC(msa->ss_cons, sizeof(char) * (msa->alen+1)); msa->ss_cons[msa->alen] = '\0'; 
    memset(msa->ss_cons,  '.', msa->alen);
  }  

  /* get SS_cons from the msa possibly for --mapstr, important to do it here, before it is potentially deknotted in HandModelmaker() */
  if(msa->ss_cons != NULL) { 
    /* post 1.1.1 release modification [EPN, Tue Jul 28 15:34:45 2015] 
     * Be careful to deal with basepairs where exactly one of the two paired positions is a consensus 
     * position (other is an insert). The way we deal is to remove the structure annotation for the 
     * one that is a consensus position. Replace it with a '.'.
     */
    /* set i_am_rf array, i_am_rf[apos] = 1 if apos is a consensus position, else it's 0 */
    ESL_ALLOC(i_am_rf, sizeof(int) * (msa->alen+1));
    esl_vec_ISet(i_am_rf, (msa->alen+1), 0);
    for(cpos = 1; cpos <= cm->clen; cpos++) i_am_rf[cm->map[cpos]] = 1;

    /* get CT array that describes all basepairs in the ss_cons (ct array is 1..alen, not 0..alen-1 */
    if((status = esl_strdup(msa->ss_cons, msa->alen, &msa_ss_cons_copy)) != eslOK) cm_Fail("Out of memory");
    if((status = esl_wuss2ct(msa_ss_cons_copy, msa->alen, msa_ct)) != eslOK) cm_Fail("Problem including structure from --mapali, maybe out of memory");
    for(apos = 1; apos <= msa->alen; apos++) { 
      if(i_am_rf[apos] && msa_ct[apos] != 0) { 
        /* apos is a consensus position that is part of a pair, make
         * sure it's mate is also consensus, if not, remove the
         * annotation of both of them from the consensus structure
         */
        opos = msa_ct[apos];
        if(! i_am_rf[opos]) { 
          msa_ss_cons_copy[apos-1] = '.';
          msa_ss_cons_copy[opos-1] = '.';
        }
      }
    }

    ESL_ALLOC(ss, sizeof(char) * (cm->clen+1));
    ss[cm->clen] = '\0';
    for(cpos = 1; cpos <= cm->clen; cpos++) ss[cpos-1] = msa_ss_cons_copy[cm->map[cpos]-1];
  }
  else { 
    cm_Fail("--mapali MSA in %s does not have any SS_cons annotation, use --noss if you used --noss with cmbuild", msafile);
  }

  /* add RF annotation to the msa, so we can use it in HandModelMaker() */
  if(msa->rf != NULL) free(msa->rf);
  ESL_ALLOC(msa->rf, sizeof(char) * (msa->alen+1));
  /* init to all inserts, then set match states based on cm->map */
  for (apos = 0; apos <  msa->alen; apos++) msa->rf[apos] = '.';
  for (cpos = 1; cpos <= cm->clen;  cpos++) msa->rf[cm->map[cpos]-1] = 'x'; /* note off by one */

  /* create a guide tree, which we'll need to convert aligned sequences to parsetrees */
  status = HandModelmaker(msa, errbuf, 
			  TRUE,  /* use_rf */
			  FALSE, /* use_el, no */
			  FALSE, /* use_wts, irrelevant */
			  0.5,   /* symfrac, irrelevant */
			  NULL,  /* returned CM, irrelevant */
			  &mtr); /* guide tree */
  if(status != eslOK) return status;

  /* create a parsetree from each aligned sequence */
  ESL_ALLOC(used_el, sizeof(int)  * (msa->alen+1));
  used_el[0] = FALSE; /* invalid */
  for(apos = 0; apos < msa->alen; apos++) { 
    used_el[apos+1] = (msa->rf[apos] == '~') ? TRUE : FALSE;
  }
  ESL_ALLOC(a2u_map, sizeof(int)  * (msa->alen+1));
  a2u_map[0] = -1; /* invalid */
  for (i = 0; i < msa->nseq; i++) { 
    if((status = Transmogrify(cm, errbuf, mtr, msa->ax[i], used_el, msa->alen, &(dataA[i]->tr))) != eslOK) return status;
    /* dataA[i]->tr is in alignment coords, convert it to unaligned coords.
     * First we construct a map of aligned to unaligned coords, then
     * we use it to convert. 
     */
    uapos = 1;
    for(apos = 1; apos <= msa->alen; apos++) { 
      a2u_map[apos] = (esl_abc_XIsGap(msa->abc, msa->ax[i][apos])) ? -1 : uapos++; 
    }
    for(x = 0; x < dataA[i]->tr->n; x++) { 
      if(dataA[i]->tr->emitl[x] != -1) dataA[i]->tr->emitl[x] = a2u_map[dataA[i]->tr->emitl[x]];
      if(dataA[i]->tr->emitr[x] != -1) dataA[i]->tr->emitr[x] = a2u_map[dataA[i]->tr->emitr[x]];
    }
  }

  /* get sequences */
  for (i = 0; i < msa->nseq; i++) esl_sq_FetchFromMSA(msa, i, &(dataA[i]->sq));

  *ret_dataA = dataA;
  *ret_ndata = msa->nseq;
  *ret_ss    = ss;

  esl_msafile_Close(afp);
  esl_msa_Destroy(msa);
  FreeParsetree(mtr);
  free(a2u_map);
  free(used_el);
  if(i_am_rf          != NULL) free(i_am_rf);
  if(msa_ct           != NULL) free(msa_ct);
  if(msa_ss_cons_copy != NULL) free(msa_ss_cons_copy);

  return eslOK;

 ERROR:
  *ret_ndata = 0;
  *ret_dataA = NULL;
  if (dataA     != NULL) { 
    for(i = 0; i < msa->nseq; i++) cm_alndata_Destroy(dataA[i], TRUE); 
    dataA = NULL;
  }
  if (afp             != NULL) esl_msafile_Close(afp);
  if (msa             != NULL) esl_msa_Destroy(msa);
  if (a2u_map         != NULL) free(a2u_map);
  if (aseq            != NULL) free(aseq);  
  if(i_am_rf          != NULL) free(i_am_rf);
  if(msa_ct           != NULL) free(msa_ct);
  if(msa_ss_cons_copy != NULL) free(msa_ss_cons_copy);

  ESL_FAIL(status, errbuf, "out of memory");
}

static int
output_alignment(ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, FILE *ofp, CM_ALNDATA **dataA, int ndata, char *map_sscons)
{
  int           status;
  ESL_MSA      *msa = NULL;
  int           j;
  ESL_SQ      **sqpA = NULL;   /*  array of sequence pointers,  only nec for Parsetrees2Alignment() */
  Parsetree_t **trA  = NULL;   /*  array of Parsetree pointers, only nec for Parsetrees2Alignment() */
  char        **ppstrA = NULL; /*  array of PP string pointers, only nec for Parsetrees2Alignment() */
  float         sc;
  float         struct_sc;
  int           first_ali = (dataA[0]->idx == 0) ? TRUE : FALSE;
  int           cpos, apos;   /* counters over consensus positions, alignment positions */

  /* contract check */
  if(ofp == cfg->tmpfp && esl_opt_GetBoolean(go, "--ileaved")) ESL_FAIL(eslEINVAL, errbuf, "--ileaved enabled, but trying to output to temporary alignment file. This shouldn't happen.");
  if(ofp == cfg->tmpfp && cfg->outfmt != eslMSAFILE_PFAM && cfg->outfmt != eslMSAFILE_STOCKHOLM) ESL_FAIL(eslEINVAL, errbuf, "output format not Stockholm, nor Pfam, but trying to output to temporary alignment file. This shouldn't happen.");
  if(esl_opt_GetBoolean(go, "--mapstr") && map_sscons == NULL) ESL_FAIL(eslEINVAL, errbuf, "--mapstr enabled, but SS_cons not read from the --mapali alignment.");

  /* output the parsetrees, if nec */
  if(cfg->tfp != NULL) { 
    for (j = 0; j < ndata; j++) { 
      if((status = ParsetreeScore(cm, NULL, errbuf, dataA[j]->tr, dataA[j]->sq->dsq, FALSE, &sc, &struct_sc, NULL, NULL, NULL)) != eslOK) return status;
      fprintf(cfg->tfp, ">%s\n", dataA[j]->sq->name);
      fprintf(cfg->tfp, "  %16s %.2f bits\n", "SCORE:", sc);
      fprintf(cfg->tfp, "  %16s %.2f bits\n", "STRUCTURE SCORE:", struct_sc);
      ParsetreeDump(cfg->tfp, dataA[j]->tr, cm, dataA[j]->sq->dsq);
      fprintf(cfg->tfp, "//\n");
    }
  }

  /* print per-CM info to insertfp and elfp, if nec */
  if(first_ali && cfg->ifp != NULL) { fprintf(cfg->ifp, "%s %d\n", cm->name, cm->clen); } 
  if(first_ali && cfg->efp != NULL) { fprintf(cfg->efp, "%s %d\n", cm->name, cm->clen); } 

  /* create the alignment */
  ESL_ALLOC(sqpA,   sizeof(ESL_SQ *)      * ndata); for(j = 0; j < ndata; j++) sqpA[j]   = dataA[j]->sq;
  ESL_ALLOC(trA,    sizeof(Parsetree_t *) * ndata); for(j = 0; j < ndata; j++) trA[j]    = dataA[j]->tr;
  ESL_ALLOC(ppstrA, sizeof(char *)        * ndata); for(j = 0; j < ndata; j++) ppstrA[j] = dataA[j]->ppstr;
  if((status = Parsetrees2Alignment(cm, errbuf, cfg->abc_out, sqpA, NULL, trA, ppstrA, ndata, cfg->ifp, cfg->efp, 
                                    /*do_full=*/TRUE, 
                                    /*do_matchonly=*/esl_opt_GetBoolean(go, "--matchonly"), 
                                    /*allow_trunc=*/esl_opt_GetBoolean(go, "--miss"), 
                                    &msa)) != eslOK) return status;

  if(ofp == cfg->rfp) { /* --regress file, remove GF author annotation */
    free(msa->au);
    msa->au = NULL;
  }

  /* rewrite SS_cons if --mapstr used */
  if(esl_opt_GetBoolean(go, "--mapstr")) { 
    /* step along the existing SS_cons, overwriting consensus positions in place */
    cpos = 0; /* span 0..clen-1 */
    for(apos = 0; apos < msa->alen; apos++) { /* span 0..alen-1 */
      if((! esl_abc_CIsGap    (cm->abc, msa->rf[apos])) && 
	 (! esl_abc_CIsMissing(cm->abc, msa->rf[apos]))) { 
	msa->ss_cons[apos] = map_sscons[cpos++];    
      }
    }
  }

  /* Determine format: if we're printing to a tmpfile we must use
   * Pfam format, so we can go back later and merge all alignments
   * in the tmpfile. If we're not printing to a tmpfile, then we
   * are about to output the full alignment in one block, and 
   * we do that in the output format cfg->outfmt. We've checked
   * that this all makes sense earlier in the program, and the
   * contract of this function asserted so (see above).
   */
  status = esl_msafile_Write(ofp, msa, (ofp == cfg->tmpfp ? eslMSAFILE_PFAM : cfg->outfmt));
  if      (status == eslEMEM) ESL_FAIL(status, errbuf, "Memory error when outputting alignment\n");
  else if (status != eslOK)   ESL_FAIL(status, errbuf, "Writing alignment file failed with error %d\n", status);

  if(msa    != NULL) esl_msa_Destroy(msa);
  if(sqpA   != NULL) free(sqpA);
  if(trA    != NULL) free(trA);
  if(ppstrA != NULL) free(ppstrA);

  return eslOK;

 ERROR:
  if(msa    != NULL) esl_msa_Destroy(msa);
  if(sqpA   != NULL) free(sqpA);
  if(trA    != NULL) free(trA);
  if(ppstrA != NULL) free(ppstrA);
  return status;
}


/* Function: output_info_file_header
 * Date:     EPN, Fri Dec  4 08:15:31 2009
 *
 * Purpose:  Print the header section of an insert or EL insert
 *           (--ifile, --elfile) information file.
 *
 * Returns:  void
 */
void
output_info_file_header(FILE *fp, char *firstline, char *elstring)
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
  fprintf(fp, "#   <c_x> is a consensus position (between 0 and <clen>; if 0: inserts before 1st consensus posn)\n");
  fprintf(fp, "#   <u_x> is the *unaligned* position (b/t 1 and <seqlen>) in <seqname> of the first %sinserted residue after <c_x>.\n", elstring);
  fprintf(fp, "#   <i_x> is the number of %sinserted residues after position <c_x> for <seqname>.\n", elstring);
  fprintf(fp, "# Lines for sequences with 0 %sinserted residues will include only <seqname> <seqlen> <spos> <epos>.\n", elstring);
  fprintf(fp, "# The final non-'#' prefixed line per model includes only '//', indicating the end of info for a model.\n");
  fprintf(fp, "#\n");

  return;
}

/* Function: output_scores()
 * Date:     EPN, Tue Jan  3 14:49:56 2012
 *
 * Purpose:  Print scores and other information to a scores file.
 *
 * Returns:  eslOK on success.
 *           eslEMEM if out of memory.
 */
int
output_scores(FILE *ofp, CM_t *cm, char *errbuf, CM_ALNDATA **dataA, int ndata, int first_idx, int be_verbose)
{
  int   status;               /* easel status */
  int   i;                    /* counter */
  int   namewidth = 8;        /* length of 'seq name' */
  char *namedashes = NULL;    /* namewidth-long string of dashes */
  int   idxwidth;    ;        /* length of max index */
  char *idxdashes = NULL;     /* idxwidth-long string of dashes */
  int64_t maxidx;             /* maximum index */
  
  /* alignment options */
  int do_nonbanded = (cm->align_opts & CM_ALIGN_NONBANDED) ? TRUE : FALSE;
  int do_post      = (cm->align_opts & CM_ALIGN_POST)      ? TRUE : FALSE;
  int do_sub       = (cm->align_opts & CM_ALIGN_SUB)       ? TRUE : FALSE;
  int do_trunc     = (cm->align_opts & CM_ALIGN_TRUNC)     ? TRUE  : FALSE;

  for(i = first_idx; i < ndata; i++) namewidth = ESL_MAX(namewidth, strlen(dataA[i]->sq->name));

  maxidx = dataA[ndata-1]->idx+1;
  idxwidth = 0; do { idxwidth++; maxidx/=10; } while (maxidx); /* poor man's (int)log_10(maxidx)+1 */
  idxwidth = ESL_MAX(idxwidth, 3);

  ESL_ALLOC(namedashes, sizeof(char) * (namewidth+1));
  namedashes[namewidth] = '\0';
  for(i = 0; i < namewidth; i++) namedashes[i] = '-';

  ESL_ALLOC(idxdashes, sizeof(char) * (idxwidth+1));
  idxdashes[idxwidth] = '\0';
  for(i = 0; i < idxwidth; i++) idxdashes[i] = '-';

  fprintf(ofp, "# %*s  %-*s  %6s  %7s  %7s  %5s  %8s  %6s  %-30s  %8s",    idxwidth, "",          namewidth,         "",      " ",        "",        "",      "",         "",       "", "       running time (s)",         "");
  if(be_verbose) fprintf(ofp, "  %7s  %7s  %7s  %8s", "", "", "", "");
  fprintf(ofp, "\n");

  fprintf(ofp, "# %*s  %-*s  %6s  %7s  %7s  %5s  %8s  %6s  %30s  %8s",     idxwidth, "",          namewidth,         "",      " ",        "",        "",      "",         "",       "", "-------------------------------", "");
  if(be_verbose) fprintf(ofp, "  %7s  %7s  %7s  %8s", "", "", "", "");
  fprintf(ofp, "\n");

  fprintf(ofp, "# %*s  %-*s  %6s  %7s  %7s  %5s  %8s  %6s  %9s  %9s  %9s  %8s", idxwidth, "idx",   namewidth, "seq name", "length", "cm from",   "cm to", "trunc",   "bit sc", "avg pp", "band calc", "alignment", "total", "mem (Mb)");
  if(be_verbose) fprintf(ofp, "  %7s  %7s  %7s  %8s", "tau", "thresh1", "thresh2", "failover");
  fprintf(ofp, "\n");

  fprintf(ofp, "# %*s  %-*s  %6s  %7s  %7s  %5s  %8s  %6s  %9s  %9s  %9s  %8s", idxwidth, idxdashes, namewidth, namedashes, "------", "-------", "-------", "-----", "--------", "------", "---------", "---------", "---------", "--------");
  if(be_verbose) fprintf(ofp, "  %7s  %7s  %7s  %8s", "-------", "-------", "-------", "--------");
  fprintf(ofp, "\n");

  for(i = first_idx; i < ndata; i++) { 
    fprintf(ofp, "  %*" PRId64 "  %-*s  %6" PRId64 "  %7d  %7d", idxwidth, dataA[i]->idx+1, namewidth, dataA[i]->sq->name, dataA[i]->sq->n, dataA[i]->spos, dataA[i]->epos);
    if(do_sub) { 
      if     (dataA[i]->spos != 1 && dataA[i]->epos != cm->clen) fprintf(ofp, "  %5s", "5'&3'");
      else if(dataA[i]->spos == 1 && dataA[i]->epos != cm->clen) fprintf(ofp, "  %5s", "3'");
      else if(dataA[i]->spos != 1 && dataA[i]->epos == cm->clen) fprintf(ofp, "  %5s", "5'");
      else if(dataA[i]->spos == 1 && dataA[i]->epos == cm->clen) fprintf(ofp, "  %5s", "no");
    }
    else { 
      if     (dataA[i]->tr->mode[0] == TRMODE_T) fprintf(ofp, "  %5s", "5'&3'");
      else if(dataA[i]->tr->mode[0] == TRMODE_L) fprintf(ofp, "  %5s", "3'");
      else if(dataA[i]->tr->mode[0] == TRMODE_R) fprintf(ofp, "  %5s", "5'");
      else                                       fprintf(ofp, "  %5s", "no");
    }
    fprintf(ofp, "  %8.2f", dataA[i]->sc);
    if(do_post)        fprintf(ofp, "  %6.3f", dataA[i]->pp);
    else               fprintf(ofp, "  %6s",   "-");
    if(! do_nonbanded) fprintf(ofp, "  %9.2f", dataA[i]->secs_bands);
    else               fprintf(ofp, "  %9s",   "-");
    fprintf(ofp, "  %9.2f  %9.2f", dataA[i]->secs_aln, dataA[i]->secs_tot);
    fprintf(ofp, "  %8.2f", dataA[i]->mb_tot);
    if(be_verbose) { 
      if(dataA[i]->tau > -0.5) fprintf(ofp, "  %7.2g", dataA[i]->tau); /* tau is -1. if aln did not use HMM bands */
      else                     fprintf(ofp, "  %7s", "-");
      if(do_trunc)             fprintf(ofp, "  %7.2f  %7.2f", dataA[i]->thresh1, dataA[i]->thresh2); 
      else                     fprintf(ofp, "  %7s  %7s", "-", "-");
      if(do_trunc && (! do_nonbanded)) { 
	fprintf(ofp, "  %8s", (dataA[i]->tr->is_std) ? "yes" : "no");
      }
      else {
	fprintf(ofp, "  %8s", "-");
      }
    }
    fprintf(ofp, "\n");
  }

  if(namedashes != NULL) free(namedashes);
  if(idxdashes  != NULL) free(idxdashes);
  return eslOK;

 ERROR: 
  ESL_FAIL(status, errbuf, "out of memory");
  return status; /* NEVER REACHED */
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
 * Args:     go      - options
 *           cfg     - cmalign config
 *           errbuf  - for error messages
 *           cm      - CM used for alignment, useful for cm->clen
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
create_and_output_final_msa(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nali, char *tmpfile) 
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
  /*int           alen_fmsa;*/                  /* number of columns in merged MSA */
  int          *ngap_insA = NULL;               /* [0..alen] number of insert gap columns to add after each alignment column when merging */
  int          *ngap_elA = NULL;                /* [0..alen] number of missing data ('~') gap columns to add after each alignment column when merging */
  int          *ngap_eitherA = NULL;            /* [0..apos..alen] = ngap_insA[apos] + ngap_elA[apos] */
  char         *rf2print = NULL;                /* #=GC RF annotation for final alignment */
  char         *ss_cons2print = NULL;           /* #=GC SS_cons annotation for final alignment */

  /* variables only used in small mode */
  int           ngs_cur;                       /* number of GS lines in current alignment (only used if do_small) */
  int           gs_exists = FALSE;             /* set to TRUE if do_small and any input aln has >= 1 GS line */
  int           maxname, maxgf, maxgc, maxgr;  /* max length of seqname, GF tag, GC tag, GR tag in all input alignments */
  int           maxname_cur, maxgf_cur, maxgc_cur, maxgr_cur; /* max length of seqname, GF tag, GC tag, GR tag in current input alignment */
  int           margin = 0;                    /* total margin length for output msa */
  int           regurg_header = FALSE;         /* set to TRUE if we're printing out header */
  int           regurg_gf     = FALSE;         /* set to TRUE if we're printing out GF */
  ESL_MSAFILE2 *afp;

  /* Allocate and initialize */
  ESL_ALLOC(msaA,   sizeof(ESL_MSA *) * nali);
  ESL_ALLOC(alenA,  sizeof(int64_t) * nali);

  /****************************************************************************
   * Read alignments one at a time, storing all non-sequence info, separately *
   ****************************************************************************/
  if((status = esl_msafile2_Open(tmpfile, NULL, &afp)) != eslOK) cm_Fail("unable to open temp file %s for reading", tmpfile);

  ai = 0;
  nseq_tot = 0;
  maxname = maxgf = maxgc = maxgr = 0;

  /* allocate maxins */
  ESL_ALLOC(maxins, sizeof(int) * (cm->clen+1)); 
  esl_vec_ISet(maxins, (cm->clen+1), 0);
  /* allocate maxel */
  ESL_ALLOC(maxel, sizeof(int) * (cm->clen+1)); 
  esl_vec_ISet(maxel, (cm->clen+1), 0); /* these will all stay 0 unless we see '~' in the alignments */

  /* read all alignments, there should be nali of them */
  for(ai = 0; ai < nali; ai++) { 
    status = esl_msafile2_ReadInfoPfam(afp, NULL, cfg->abc, -1, NULL, NULL, &(msaA[ai]), &nseq_cur, &alen_cur, &ngs_cur, &maxname_cur, &maxgf_cur, &maxgc_cur, &maxgr_cur, NULL, NULL, NULL, NULL, NULL);
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
  /* final check, make sure we've read all msas from the file, we should have, we only printed nali */
  status = esl_msafile2_ReadInfoPfam(afp, NULL, cfg->abc, -1, NULL, NULL, NULL, 
				     NULL, NULL, NULL, NULL, NULL, NULL, NULL, 
				     NULL, NULL, NULL, NULL, NULL);
  if(status != eslEOF) ESL_FAIL(status, errbuf, "More alignments in temp file than expected.");
  esl_msafile2_Close(afp);
  
  /********************************************
   * Merge all alignments into the merged MSA *
   ********************************************/

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
  /*alen_fmsa = cm->clen + esl_vec_ISum(maxins, (cm->clen+1));*/

  /* if there was any GS annotation in any of the individual alignments,
   * do second pass through alignment files, outputting GS annotation as we go. */
  if(gs_exists) { 
    if((status = esl_msafile2_Open(tmpfile, NULL, &afp)) != eslOK) cm_Fail("unable to open temp file %s for reading on second pass", tmpfile);
    for(ai = 0; ai < nali; ai++) { 
      regurg_header = (ai == 0) ? TRUE : FALSE;
      regurg_gf     = (ai == 0) ? TRUE : FALSE;
      status = esl_msafile2_RegurgitatePfam(afp, cfg->ofp, 
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
      fflush(cfg->ofp);
    }
    esl_msafile2_Close(afp);
    fprintf(cfg->ofp, "\n"); /* a single blank line to separate GS annotation from aligned data */
  }
  /* do another (either second or third) pass through alignment files, outputting aligned sequence data (and GR) as we go */

  if((status = esl_msafile2_Open(tmpfile, NULL, &afp)) != eslOK) cm_Fail("unable to open temp file %s for reading on second (or third) pass", tmpfile);

  for(ai = 0; ai < nali; ai++) { 
    /* determine how many all gap columns to insert after each alignment position
     * of the temporary msa when copying it to the merged msa */
    if((status = determine_gap_columns_to_add(msaA[ai], maxins, maxel, cm->clen, &(ngap_insA), &(ngap_elA), &(ngap_eitherA), errbuf)) != eslOK) 
      cm_Fail("error determining number of all gap columns to add to temp alignment %d\n%s", ai, errbuf);
    regurg_header = ((! gs_exists) && (ai == 0)) ? TRUE : FALSE;
    regurg_gf     = ((! gs_exists) && (ai == 0)) ? TRUE : FALSE;

    status = esl_msafile2_RegurgitatePfam(afp, cfg->ofp,
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
      inflate_gc_with_gaps_and_els(cfg->ofp, msaA[ai], ngap_insA, ngap_elA, &ss_cons2print, &rf2print);
    }
    free(ngap_insA);
    free(ngap_elA);
    free(ngap_eitherA);
    
    esl_msa_Destroy(msaA[ai]);
    msaA[ai] = NULL;
    fflush(cfg->ofp);
  }
  /* output SS_cons and RF */
  margin = maxname+1;
  if (maxgc > 0 && maxgc+6         > margin) margin = maxgc+6;
  if (maxgr > 0 && maxname+maxgr+7 > margin) margin = maxname+maxgr+7; 
  fprintf(cfg->ofp, "#=GC %-*s %s\n", margin-6, "SS_cons", ss_cons2print);
  fprintf(cfg->ofp, "#=GC %-*s %s\n", margin-6, "RF", rf2print);
  fprintf(cfg->ofp, "//\n");

  esl_msafile2_Close(afp);

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

  /* contract check */
  if(maxel[0]     != 0)    ESL_FAIL(eslEINVAL, errbuf, "missing characters exist prior to first cpos, this shouldn't happen.\n");
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
      if(nins > 0) ESL_FAIL(eslEINVAL, errbuf, "after nongap RF pos %d, %d gap columns precede a missing data column (none should)", cpos, nins);
    }
    else if(esl_abc_CIsGap(msa->abc, msa->rf[apos])) { 
      nins++;
    }
    else { /* a consensus position */
      /* a few sanity checks */
      if(nins > maxins[cpos])  ESL_FAIL(eslEINCONCEIVABLE, errbuf, "%d inserts before cpos %d greater than max expected (%d).\n", nins, cpos, maxins[cpos]); 
      if(nel  > maxel[cpos]) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "%d EL inserts before cpos %d greater than max expected (%d).\n", nel, cpos, maxel[cpos]); 

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
	  /* Rule is (as of SVN r4271, and version 1.1rc2) that 
	   * ELs always come before (5' of) insertions, even for the
	   * rare case of a MATP node immediately prior to an END node.
	   */
	  ngap_elA[prv_cpos  + 1 +       (nel/2)]  = maxel[cpos]  - nel;  /* internal cpos: split */
	  ngap_insA[prv_cpos + 1 + nel + (nins/2)] = maxins[cpos] - nins; /* internal cpos: split */
	}
	/* final case is if (maxins[cpos] == 0 && maxel[cpos] == 0) 
	 * in this case we do nothing. 
	 */
      }

      cpos++;
      prv_cpos = apos;
      nins = 0;
      nel = 0;
    }
  }
  /* first, validate that clen is what it should be */
  if(cpos != clen) { 
    if(ngap_insA != NULL) free(ngap_insA);
    if(ngap_elA != NULL) free(ngap_elA);
    if(ngap_eitherA != NULL) free(ngap_eitherA);
    ESL_FAIL(eslEINCONCEIVABLE, errbuf, "consensus length (%d) is not the expected length (%d).", cpos, clen);
  }
  
  if(maxins[cpos] > 0 && maxel[cpos] == 0) { /* most common case */
    ngap_insA[prv_cpos + 1 + nins] = maxins[cpos] - nins; /* flush left inserts (no missing) */
  }
  else if(maxins[cpos] == 0 && maxel[cpos] > 0) { 
    ngap_elA[prv_cpos + 1 + nel] = maxel[cpos] - nel; /* flush left ELs (no gaps) */
  }
  else if(maxins[cpos] > 0 && maxel[cpos] > 0) { 
    /* missing data (ELs) is always 5' of gaps */
    ngap_elA[prv_cpos + 1 + nel]         = maxel[cpos] - nel; /* flush left */
    ngap_insA[prv_cpos + 1 + nel + nins] = maxins[cpos] - nins; /* flush left, after ELs */
  }

  /* determine ngap_eitherA[], the number of gaps due to either inserts or missing data after each apos */
  for(apos = 0; apos <= msa->alen; apos++) { 
    ngap_eitherA[apos] = ngap_insA[apos] + ngap_elA[apos];
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
  int i;
  int alen2print = 0;
  char *rf2print;
  char *ss_cons2print;

  alen2print = msa->alen + esl_vec_ISum(ngap_insA, msa->alen+1) + esl_vec_ISum(ngap_elA, msa->alen+1);
  ESL_ALLOC(rf2print,      sizeof(char) * (alen2print+1));
  ESL_ALLOC(ss_cons2print, sizeof(char) * (alen2print+1));
  rf2print[alen2print] = '\0';
  ss_cons2print[alen2print] = '\0';

  if(msa->ss_cons == NULL) cm_Fail("Error: trying to add inserts to SS_cons, but unexpectedly it doesn't exist.");
  if(msa->rf      == NULL) cm_Fail("Error: trying to add inserts to SS_cons, but unexpectedly it doesn't exist.");
  for(apos = 0; apos <= msa->alen; apos++) { 
    /* ELs always come before (5' of) inserts */
    for(i = 0; i < ngap_elA[apos]; i++) { 
      rf2print[apos2print] = '~';
      ss_cons2print[apos2print++] = '~';
    }
    for(i = 0; i < ngap_insA[apos]; i++) { 
      rf2print[apos2print] = '.';
      ss_cons2print[apos2print++] = '.';
    }
    if(apos < msa->alen) { 
      rf2print[apos2print]        = msa->rf[apos];
      ss_cons2print[apos2print++] = msa->ss_cons[apos];
    }	
  }    
  
  *ret_ss_cons2print = ss_cons2print;
  *ret_rf2print = rf2print;
  
  return;
  
 ERROR:
  cm_Fail("Allocation error when creating final alignment RF and SS_cons.");
  return; /* NEVERREACHED */
}

/* configure_root_inserts
 *                   
 * Modify the transition probabilities into and out of the 
 * ROOT_IL and ROOT_IR states. 
 * The motivation is to allow cmalign to more accurately
 * align sequences that have extra nonhomologous sequence
 * on the ends. Defaultly-paramaterized models (especially
 * those with zero basepairs) tend to mess up the alignment
 * at the ends if there are extra nucleotides. 
 *
 * The value of 'prob' was limited to be 0. < prob < 0.4 by
 * getopts but we also do a sanity check here.
 * 
 * Returns void.
 */
void
configure_root_inserts(CM_t *cm, float to_insert_prob, float self_insert_prob)
{
  float state0_sum = 0.; /* will store cumulative prob of transitinos out of state 1 */
  float state1_sum = 0.; /* will store cumulative prob of transitinos out of state 1 */
  float state2_sum = 0.; /* will store cumulative prob of transitinos out of state 1 */

  if((to_insert_prob <= 0.) || to_insert_prob > 0.4) { 
    cm_Fail("ERROR with --flanktoins <x>, <x> should be > 0. and < 0.4");
  }
  if((to_insert_prob <= 0.) || to_insert_prob > 0.4) { 
    cm_Fail("ERROR with --flankselfins <x>, <x> should be > 0. and < 0.9");
  }
  if((to_insert_prob + self_insert_prob) > 0.95) { 
    cm_Fail("ERROR with --flanktoins <x1> and --flankselfins <x2>, <x1> + <x2> must be less than 0.95");
  }

  /* Deal with transitions out of ROOT_S first 
   * first  transition out of ROOT_S is always to ROOT_IL 
   * second transition out of ROOT_S is always to ROOT_IR 
   * third  transition out of ROOT_S is always to 'match' state in split set of next node, e.g. MATL_ML, MATR_MR, MATP_MP or BIF_B 
   * remaining number of transitions depend on next node type, but are all 'delete' states (unless BIF) 
   */
  cm->t[0][0] = to_insert_prob;   /* ROOT_S  -> ROOT_IL */
  cm->t[0][1] = to_insert_prob;   /* ROOT_S  -> ROOT_IR */
  state0_sum = to_insert_prob + to_insert_prob;

  cm->t[1][0] = self_insert_prob; /* ROOT_IL -> ROOT_IL */
  cm->t[1][1] = to_insert_prob;   /* ROOT_IL -> ROOT_IR */
  state1_sum = self_insert_prob + to_insert_prob;

  cm->t[2][0] = self_insert_prob; /* ROOT_IR -> ROOT_IR */
  state2_sum = self_insert_prob;
  
  /* 3/4 of the remaining to_insert_probability goes to the match state, unless BIF_B */
  if(cm->ndtype[1] == BIF_nd) { 
    cm->t[0][2] = 1. - state0_sum; /* ROOT_S  -> BIF_B */

    cm->t[1][2] = 1. - state1_sum; /* ROOT_IL -> BIF_B */

    cm->t[2][1] = 1. - state2_sum; /* ROOT_IR -> BIF_B */
  }
  else if(cm->ndtype[1] == MATP_nd) { 
    cm->t[0][2] =  (1. - state0_sum) * 0.75;       /* ROOT_S -> MATP_MP */
    cm->t[0][3] = ((1. - state0_sum) * 0.25) / 3.; /* ROOT_S -> MATP_ML */
    cm->t[0][4] = ((1. - state0_sum) * 0.25) / 3.; /* ROOT_S -> MATP_MR */
    cm->t[0][5] = ((1. - state0_sum) * 0.25) / 3.; /* ROOT_S -> MATP_D */

    cm->t[1][2] =  (1. - state1_sum) * 0.75;       /* ROOT_IL -> MATP_MP */
    cm->t[1][3] = ((1. - state1_sum) * 0.25) / 3.; /* ROOT_IL -> MATP_ML */
    cm->t[1][4] = ((1. - state1_sum) * 0.25) / 3.; /* ROOT_IL -> MATP_MR */
    cm->t[1][5] = ((1. - state1_sum) * 0.25) / 3.; /* ROOT_IL -> MATP_D */

    cm->t[2][1] =  (1. - state2_sum) * 0.75;       /* ROOT_IR -> MATP_MP */
    cm->t[2][2] = ((1. - state2_sum) * 0.25) / 3.; /* ROOT_IR -> MATP_ML */
    cm->t[2][3] = ((1. - state2_sum) * 0.25) / 3.; /* ROOT_IR -> MATP_MR */
    cm->t[2][4] = ((1. - state2_sum) * 0.25) / 3.; /* ROOT_IR -> MATP_MR */
  }  
  else if((cm->ndtype[1] == MATL_nd) || (cm->ndtype[1] == MATR_nd)) { 
    cm->t[0][2] = (1. - state0_sum) * 0.75; /* ROOT_S -> MAT{L,R}_M{L,R} */
    cm->t[0][3] = (1. - state0_sum) * 0.25; /* ROOT_S -> MAT{L,R}_D */

    cm->t[1][2] = (1. - state1_sum) * 0.75; /* ROOT_IL -> MAT{L,R}_M{L,R} */
    cm->t[1][3] = (1. - state1_sum) * 0.25; /* ROOT_IL -> MAT{L,R}_D */

    cm->t[2][1] = (1. - state2_sum) * 0.75; /* ROOT_IR -> MAT{L,R}_M{L,R} */
    cm->t[2][2] = (1. - state2_sum) * 0.25; /* ROOT_IR -> MAT{L,R}_D */
  }  
  else { 
    cm_Fail("ERROR, with --flankins, unexpected second node type, not one of BIF, MATP, MATL or MATR");
  }

  /* should already be normalized, but to be safe: */
  esl_vec_FNorm(cm->t[0], cm->cnum[0]);
  esl_vec_FNorm(cm->t[1], cm->cnum[1]);
  esl_vec_FNorm(cm->t[2], cm->cnum[2]);

  return;
}


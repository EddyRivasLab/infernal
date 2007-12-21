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

#include "easel.h"              /* better general sequence analysis library */
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_mpi.h"
#include "esl_msa.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"

#include "funcs.h"		/* external functions                   */
#include "structs.h"		/* data structures, macros, #define's   */

#define STRATOPTS1 "--cyk,--inside,--hmmviterbi,--hmmforward" /* incompatible with --cyk, --inside */
#define STRATOPTS2 "--cyk,--inside,--hmmviterbi,--hmmforward,--fgiven" /* incompatible with --hmmviterbi,--hmmforward */
#define ALPHOPTS   "--rna,--dna"                              /* exclusive choice for output alphabet */

#define I_CMCUTOPTS   "-E,-T,--ga,--tc,--nc,--hmmviterbi,--hmmforward" /* exclusive choice for CM cutoff */
#define I_HMMCUTOPTS1  "--hmmcalcthr,--hmmE,--hmmT,--hmmviterbi,--hmmforward" /* exclusive choice for HMM cutoff set 1 */
#define I_HMMCUTOPTS2  "--hmmcalcthr,--hmmE,--hmmT"            /* exclusive choice for HMM cutoff set 2 */
#define FOPTS0        "--fhmmviterbi,--fhmmforward,--hmmviterbi,--hmmforward"  /* incompatible with --fgiven */
#define FOPTS1        "--fcyk,--finside,--hmmviterbi,--hmmforward"           /* incompatible with --fcyk and --finside */
#define FOPTS2        "--fhmmviterbi,--fhmmforward,--fgiven,--hmmviterbi,--hmmforward" /* incompatible with --fhmmviterbi and --fhmmforward */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs       incomp  help  docgroup*/
  /* basic options */
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",   1 },
  { "-g",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--hmmviterbi,--hmmforward", "configure CM for glocal alignment [default: local]", 1 },
  { "--informat",eslARG_STRING, NULL,  NULL, NULL,      NULL,      NULL,        NULL, "specify the input file is in format <x>, not FASTA", 1 },
  { "--toponly", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "only search the top strand", 1 },
  { "--bottomonly", eslARG_NONE,FALSE, NULL, NULL,      NULL,      NULL,        NULL, "only search the bottom strand", 1 },
  { "--window",  eslARG_INT,    NULL,  NULL, "n>0",     NULL,      NULL,        NULL, "set scanning window size to <n> [default: calculated]", 1 },
  { "--null2",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--hmmviterbi,--hmmforward", "turn on the post hoc second null model", 1 },
  { "--iins",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "allow informative insert emissions, do not zero them", 1 },
  { "--rtrans",  eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--hmmviterbi,--hmmforward", "replace CM transition scores from <cmfile> with RSEARCH scores", 1 },
  { "--greedy",  eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--hmmviterbi,--hmmforward", "resolve overlapping hits with a greedy algorithm a la RSEARCH", 1 },
  /* 4 --p* options below are hopefully temporary b/c if we have E-values for the CM using a certain cm->pbegin, cm->pend,
   * changing those values in cmsearch invalidates the E-values, so we should pick hard-coded values for cm->pbegin cm->pend */
  { "--pebegin", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "-g,--pbegin","set all local begins as equiprobable", 1 },
  { "--pfend",   eslARG_REAL,   NULL,  NULL, "0<x<1",   NULL,      NULL, "-g,--pend",  "set all local end probs to <x>", 1 },
  { "--pbegin",  eslARG_REAL,  "0.05",NULL,  "0<x<1",   NULL,      NULL,        "-g", "set aggregate local begin prob to <x>", 1 },
  { "--pend",    eslARG_REAL,  "0.05",NULL,  "0<x<1",   NULL,      NULL,        "-g", "set aggregate local end prob to <x>", 1 },
  /* options for algorithm for final round of search */
  { "--cyk",       eslARG_NONE,"default",NULL,NULL,     NULL,      NULL,    STRATOPTS1, "use scanning CM CYK algorithm", 2 },
  { "--inside",    eslARG_NONE, FALSE, NULL, NULL,      NULL,      NULL,    STRATOPTS1, "use scanning CM Inside algorithm", 2 },
  { "--hmmviterbi",eslARG_NONE, FALSE, NULL, NULL,      NULL,      NULL,    STRATOPTS2, "use scanning HMM Viterbi algorithm", 2 },
  { "--hmmforward",eslARG_NONE, FALSE, NULL, NULL,      NULL,      NULL,    STRATOPTS2, "use scanning HMM Forward algorithm", 2 },
  /* options for filtering with a CM */
  { "--fgiven",     eslARG_NONE, FALSE, NULL, NULL,     NULL,      NULL,       FOPTS0, "use filtering info from CM file", 14 },
  { "--fcyk",       eslARG_NONE, FALSE, NULL, NULL,     NULL,      NULL,       FOPTS1, "filter with CM CYK algorithm", 14 },
  { "--finside",    eslARG_NONE, FALSE, NULL, NULL,     NULL,      NULL,       FOPTS1, "filter with CM Inside algorithm", 14 },
  { "--fhmmviterbi",eslARG_NONE, FALSE, NULL, NULL,     NULL,      NULL,       FOPTS2, "filter with HMM Viterbi algorithm", 14 },
  { "--fhmmforward",eslARG_NONE, FALSE, NULL, NULL,     NULL,      NULL,       FOPTS2, "filter with HMM Forward algorithm", 14 },
  /* CM cutoff options */
  { "-E",        eslARG_REAL,   "0.1", NULL, "x>0.",    NULL,      NULL,  I_CMCUTOPTS, "use cutoff E-value of <x> for final round of CM search", 3 },
  { "-T",        eslARG_REAL,   "0.0", NULL, NULL,      NULL,      NULL,  I_CMCUTOPTS, "use cutoff bit score of <x> for final round of CM search", 3 },
  { "--ga",      eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,  I_CMCUTOPTS, "use CM Rfam GA gathering threshold as cutoff bit score", 3 },
  { "--tc",      eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,  I_CMCUTOPTS, "use CM Rfam TC trusted cutoff as cutoff bit score", 3 },
  { "--nc",      eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,  I_CMCUTOPTS, "use CM Rfam NC noise cutoff as cutoff bit score", 3 },
  { "--fE",      eslARG_REAL,   "100.",NULL, "x>0.",    NULL,      NULL,       "--fT", "use cutoff E-value of <x> for CM-filter round", 3 },
  { "--fT",      eslARG_REAL,   "0.0", NULL, NULL,      NULL,      NULL,       "--fE", "use cutoff bit score of <x> for CM-filter round", 3 },
  { "--fgreedy", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,         NULL, "resolve overlapping hits for CM-filter greedily", 1 },
  /* HMM cutoff options */
  { "--hmmcalcthr",eslARG_NONE, FALSE, NULL, NULL,     NULL,       NULL, I_HMMCUTOPTS1,"calculate HMM filter threshold by sampling from CM", 4 },
  { "--hmmE",   eslARG_REAL,   "50.", NULL, "x>0.",    NULL,       NULL, I_HMMCUTOPTS2,"use cutoff E-value of <x> for CP9 HMM filter/search", 4 },
  { "--hmmT",   eslARG_REAL,   "0.0", NULL, NULL,      NULL,       NULL, I_HMMCUTOPTS2,"use cutoff bit score of <x> for CP9 HMM filter/search", 4 },
  /* QDB related options */
  { "--beta",    eslARG_REAL,   "1E-7",NULL, "x>0",     NULL,      NULL, "--hmmviterbi,--hmmforward", "set tail loss prob for QDB and window size calculation to <x>", 5 },
  { "--fbeta",   eslARG_REAL,   "1E-4",NULL, "x>0",     NULL,      NULL,        NULL, "set QDB tail loss prob for --fcyk, --finside CM filters to <x>", 5 },
  { "--noqdb",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--hmmviterbi,--hmmforward", "DO NOT use query dependent banding (QDB) for acceleration", 5 },
  { "--qdbfile", eslARG_STRING, NULL,  NULL, NULL,      NULL,      NULL,"--hmmviterbi,--hmmforward,--noqdb","read QDBs from file <s> (outputted from cmbuild)", 5 },
  /* HMM filtering options */
  { "--hbanded", eslARG_NONE,   FALSE, NULL, NULL,      NULL,       NULL,       NULL, "calculate and use HMM bands in CM search", 6 },
  { "--tau",     eslARG_REAL,   "1e-7",NULL, "0<x<1",   NULL,"--hbanded",       NULL, "set tail loss prob for --hbanded to <x>", 6 },
  { "--scan2bands",eslARG_NONE, FALSE, NULL, NULL,      NULL,"--hbanded",       NULL, "derive HMM bands from scanning Forward/Backward", 6 },
  { "--sums",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--hbanded",       NULL, "use posterior sums during HMM band calculation (widens bands)", 6 },
  /* HMM configuration options */
  { "--hmmglocal",eslARG_NONE,  FALSE, NULL, NULL,      NULL,      NULL,        NULL, "configure HMM for glocal alignment [default: local]", 7 },
  { "--hmmnoel",  eslARG_NONE,  FALSE, NULL, NULL,      NULL,      NULL,"-g,--hmmglocal", "DO NOT enable HMM EL local ends that mirror CM", 7 },
  { "--hmmgreedy",eslARG_NONE,  FALSE, NULL, NULL,      NULL,      NULL,        NULL, "resolve HMM overlapping hits with a greedy algorithm a la RSEARCH", 7 },
  /* filter threshold calculation options */
  { "--seed",    eslARG_INT,    NULL,  NULL, "n>0",     NULL,"--hmmcalcthr",    NULL, "set random number generator seed to <n>", 8 },
  { "--N",       eslARG_INT,   "1000", NULL, "n>0",     NULL,"--hmmcalcthr",    NULL, "number of emitted sequences for HMM filter threshold calc", 8 },
  { "--F",       eslARG_REAL,  "0.95", NULL, "0<x<=1",  NULL,"--hmmcalcthr",    NULL, "required fraction of seqs that survive HMM filter", 8 },
  { "--fstep",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--hmmcalcthr",    NULL, "step from F to 1.0 while S < Starg", 8 },
  { "--starg",   eslARG_REAL,  "0.01", NULL, "0<x<=1",  NULL,"--hmmcalcthr",    NULL, "target filter survival fraction", 8 },
  { "--spad",    eslARG_REAL,  "1.0",  NULL, "0<=x<=1", NULL,"--hmmcalcthr",    NULL, "fraction of (sc(S) - sc(Starg)) to add to sc(S)", 8 },
  { "--fastfil", eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--hmmcalcthr",    NULL, "calculate filter thr quickly, assume parsetree sc is optimal", 8 },
  { "--gemit",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--hmmcalcthr",    NULL, "when calc'ing filter thresholds, always emit globally from CM", 8 },
  /* alignment options */
  { "--noalign", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,       NULL, "find start/stop/score only; don't do alignments", 9 },
  { "--optacc",  eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,"--noalign", "align hits with the Holmes/Durbin optimal accuracy algorithm", 9 },
  { "--post",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,"--noalign", "append posterior probabilities to hit alignments", 9 },
  /* Enforcing a subsequence */
  { "--enfstart",eslARG_INT,    FALSE, NULL, "n>0",     NULL,"--enfseq",        NULL, "enforce MATL stretch starting at consensus position <n>", 10 },
  { "--enfseq",  eslARG_STRING, NULL,  NULL, NULL,      NULL,"--enfstart",      NULL, "enforce MATL stretch starting at --enfstart <n> emits seq <s>", 10 },
  { "--enfnohmm",eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--enfstart",      NULL, "DO NOT filter first w/an HMM that only enforces --enfseq <s>", 10 },
  /* verbose output files */
  { "--tfile",   eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "dump parsetrees for each hit to file <f>", 11 },
  { "--gcfile",  eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "save GC content stats of target sequence file to <f>", 11 },
  { "--bfile",   eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "save bands for each state to file <f>", 11 },
  { "--filhfile",eslARG_OUTFILE, NULL, NULL, NULL,      NULL,"--hmmcalcthr",    NULL, "save CP9 filter threshold histogram(s) to file <s>", 11 },
  { "--filrfile",eslARG_OUTFILE, NULL, NULL, NULL,      NULL,"--hmmcalcthr",    NULL, "save CP9 filter threshold information file <s>", 11 },
/* Setting output alphabet */
  { "--rna",     eslARG_NONE,"default",NULL, NULL,  ALPHOPTS,      NULL,        NULL, "output alignment as RNA sequence data", 12 },
  { "--dna",     eslARG_NONE,   FALSE, NULL, NULL,  ALPHOPTS,      NULL,        NULL, "output alignment as DNA (not RNA) sequence data", 12 },
/* Other options */
  { "--stall",   eslARG_NONE,  FALSE, NULL, NULL,       NULL,      NULL,        NULL, "arrest after start: for debugging MPI under gdb",   13 },  
  { "--hmmmaxE", eslARG_REAL,   NULL, NULL, "x>0.",     NULL,"--fgiven",        NULL, "with --fgiven, set maximum HMM filter E-value to <x>", 13 },
#ifdef HAVE_MPI
  { "--mpi",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,  "--qdbfile","run as an MPI parallel program", 13 },  
#endif
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
  int           fmt;		/* format code for seqfile */
  ESL_ALPHABET *abc;		/* digital alphabet for input */
  long          N;              /* database size in nucleotides (doubled if doing rev comp) */
  int           ncm;            /* number CM we're at in file */
  int           do_rc;          /* should we search reverse complement? (for convenience */
  int           init_rci;       /* initial strand to search 0 for top, 1 for bottom (only 1 if --bottomonly enabled) */
  float        *avglen;         /* [0..v..M-1] average hit len for subtree rooted at each state v for current CM */

  int           do_mpi;		/* TRUE if we're doing MPI parallelization */
  int           nproc;		/* how many MPI processes, total */
  int           my_rank;	/* who am I, in 0..nproc-1 */
  int           do_stall;	/* TRUE to stall the program until gdb attaches */

  /* Masters only (mainly i/o streams) */
  CMFILE       *cmfp;		/* open input CM file stream       */
  FILE         *tfp;	        /* optional output for parsetrees  */
  FILE         *bfp;	        /* optional output for qdbs */
  FILE         *filhfp;	        /* optional output for filter thr calc histgram */
  FILE         *filrfp;	        /* optional output for filter thr calc R info file */
  ESL_ALPHABET *abc_out; 	/* digital alphabet for writing */
  int          *preset_dmin;    /* remains NULL unless --qdbfile, which is incompatible with --mpi */
  int          *preset_dmax;    /* remains NULL unless --qdbfile, which is incompatible with --mpi */
};

static char usage[]  = "[-options] <cmfile> <sequence file>";
static char banner[] = "align sequences to an RNA CM";

static int  init_master_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
/* static int  init_shared_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf); */

static void  serial_master (const ESL_GETOPTS *go, struct cfg_s *cfg);
#ifdef HAVE_MPI
static void  mpi_master    (const ESL_GETOPTS *go, struct cfg_s *cfg);
static void  mpi_worker    (const ESL_GETOPTS *go, struct cfg_s *cfg);
#endif
static int process_search_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, ESL_DSQ *dsq, int L, search_results_t **ret_results);
static int initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int update_avg_hit_len(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int set_gumbels(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int set_searchinfo(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int print_searchinfo(const ESL_GETOPTS *go, struct cfg_s *cfg, FILE *fp, CM_t *cm, long N, char *errbuf);
static int set_window(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int calc_filter_threshold(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float *ret_Smin);

static int read_qdb_file(FILE *fp, CM_t *cm, int *dmin, int *dmax);
static int is_integer(char *s);
static int read_next_search_seq(const ESL_ALPHABET *abc, ESL_SQFILE *seqfp, int do_revcomp, dbseq_t **ret_dbseq);
#if HAVE_MPI
static int determine_seq_chunksize(struct cfg_s *cfg, int L, int W);
#endif 

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go = NULL;   /* command line processing                     */
  ESL_STOPWATCH   *w  = esl_stopwatch_Create();
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
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      cm_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\nwhere general options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1=docgroup, 2 = indentation; 80=textwidth*/
      puts("\nstrategy choice: (exclusive) [default: CM only]");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\nCM cutoff options (exclusive) [default: E value of 0.1]");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); 
      puts("\nHMM cutoff options (exclusive) [default: bit score of 0.0]");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80);
      puts("\nquery dependent banding (QDB) related options:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80);
      puts("\nHMM filtering options: (require --fhmm)");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80);
      puts("\nHMM configuration options:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80);
      puts("\nfilter threshold calculation options: (require --hmmcalcthr)");
      esl_opt_DisplayHelp(stdout, go, 8, 2, 80);
      puts("\noptions for returning alignments of search hits:");
      esl_opt_DisplayHelp(stdout, go, 9, 2, 80);
      puts("\noptions for enforcing a single-stranded subsequence:");
      esl_opt_DisplayHelp(stdout, go, 10, 2, 80);
      puts("\nverbose output files:");
      esl_opt_DisplayHelp(stdout, go, 11, 2, 80);
      puts("\noptions for selecting output alphabet:");
      esl_opt_DisplayHelp(stdout, go, 12, 2, 80);
      puts("\nother options:");
      esl_opt_DisplayHelp(stdout, go, 13, 2, 80);
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
  /* Initialize what we can in the config structure (without knowing the input alphabet yet).
   */
  cfg.cmfile     = esl_opt_GetArg(go, 1); 
  cfg.sqfile     = esl_opt_GetArg(go, 2); 
  cfg.sqfp       = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  if   (esl_opt_IsDefault(go, "--informat")) cfg.fmt = eslSQFILE_UNKNOWN; /* autodetect sequence file format by default. */ 
  else { 
    cfg.fmt = esl_sqio_FormatCode(esl_opt_GetString(go, "--informat"));
    if(cfg.fmt == eslSQFILE_UNKNOWN) cm_Fail("Can't recognize sequence file format: %s. valid options are: fasta, embl, genbank, ddbj, uniprot, stockholm, or pfam\n", esl_opt_GetString(go, "--informat"));
  }
  cfg.abc        = NULL;	           /* created in init_master_cfg() in masters, or in mpi_worker() in workers */
  if      (esl_opt_GetBoolean(go, "--rna")) cfg.abc_out = esl_alphabet_Create(eslRNA);
  else if (esl_opt_GetBoolean(go, "--dna")) cfg.abc_out = esl_alphabet_Create(eslDNA);
  else    cm_Fail("Can't determine output alphabet");
  cfg.N          = 0;                      /* db size  */
  cfg.ncm        = 0;                      /* what number CM we're on, updated in masters, stays 0 (irrelevant) for workers */
  cfg.cmfp       = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.tfp        = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.bfp        = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.filhfp     = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.filrfp     = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.preset_dmin= NULL;                   /* filled in initialize_cm() only if --qdbfile, which conflicts with --mpi */
  cfg.preset_dmax= NULL;                   /* filled in initialize_cm() only if --qdbfile, which conflicts with --mpi */
  cfg.do_rc      = (! esl_opt_GetBoolean(go, "--toponly")); 
  cfg.init_rci   = esl_opt_GetBoolean(go, "--bottomonly") ? 1 : 0; 
  cfg.avglen     = NULL;                   /* filled in init_master_cfg() in masters, stays NULL for workers */

  cfg.do_mpi     = FALSE;	           /* this gets reset below, if we init MPI */
  cfg.nproc      = 0;		           /* this gets reset below, if we init MPI */
  cfg.my_rank    = 0;		           /* this gets reset below, if we init MPI */
  cfg.do_stall   = esl_opt_GetBoolean(go, "--stall");
  
  /* This is our stall point, if we need to wait until we get a
   * debugger attached to this process for debugging (especially
   * useful for MPI):
   */
  while (cfg.do_stall); 

  /* Start timing. */
  esl_stopwatch_Start(w);

  /* Figure out who we are, and send control there: 
   * we might be an MPI master, an MPI worker, or a serial program.
   */
#ifdef HAVE_MPI
  if (esl_opt_GetBoolean(go, "--mpi")) 
    {
      cfg.do_mpi     = TRUE;
      MPI_Init(&argc, &argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &(cfg.my_rank));
      MPI_Comm_size(MPI_COMM_WORLD, &(cfg.nproc));

      if(cfg.nproc == 1) cm_Fail("MPI mode, but only 1 processor running... (did you execute mpirun?)");

      if (cfg.my_rank > 0)  mpi_worker(go, &cfg);
      else 		    mpi_master(go, &cfg);

      esl_stopwatch_Stop(w);
      esl_stopwatch_MPIReduce(w, 0, MPI_COMM_WORLD);
      MPI_Finalize();
    }
  else
#endif /*HAVE_MPI*/
    {
      serial_master(go, &cfg);
      esl_stopwatch_Stop(w);
    }
  if (cfg.my_rank == 0) esl_stopwatch_Display(stdout, w, "# CPU time: ");

  /* Clean up the shared cfg. 
   */
  if (cfg.my_rank == 0) {
    if (cfg.cmfp      != NULL) CMFileClose(cfg.cmfp);
    if (cfg.sqfp      != NULL) esl_sqfile_Close(cfg.sqfp);
    if (cfg.tfp       != NULL) fclose(cfg.tfp);
    if (cfg.bfp       != NULL) fclose(cfg.bfp);
    if (cfg.filhfp    != NULL) fclose(cfg.filhfp);
    if (cfg.filrfp    != NULL) fclose(cfg.filrfp);
  }
  if (cfg.abc       != NULL) esl_alphabet_Destroy(cfg.abc);
  if (cfg.abc_out   != NULL) esl_alphabet_Destroy(cfg.abc_out);
  if (cfg.avglen    != NULL) free(cfg.avglen);
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
 * Allocates/Sets: 
 *    cfg->sqfp        - open sequence file                
 *    cfg->cmfp        - open CM file                
 *    cfg->tfp         - optional output file
 *    cfg->bfp         - optional output file
 *    cfg->filhfp      - optional output file
 *    cfg->filrfp      - optional output file
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

  /* open input sequence file */
  status = esl_sqfile_Open(cfg->sqfile, cfg->fmt, NULL, &(cfg->sqfp));
  if (status == eslENOTFOUND)    ESL_FAIL(status, errbuf, "File %s doesn't exist or is not readable\n", cfg->sqfile);
  else if (status == eslEFORMAT) ESL_FAIL(status, errbuf, "Couldn't determine format of sequence file %s\n", cfg->sqfile);
  else if (status == eslEINVAL)  ESL_FAIL(status, errbuf, "Can’t autodetect stdin or .gz."); 
  else if (status != eslOK)      ESL_FAIL(status, errbuf, "Sequence file open failed with error %d\n", status);
  cfg->fmt = cfg->sqfp->format;

  /* GetDBInfo() reads all sequences, rewinds seq file and returns db size */
  GetDBInfo(NULL, cfg->sqfp, &(cfg->N), NULL);  
  if ((! esl_opt_GetBoolean(go, "--toponly")) && (! esl_opt_GetBoolean(go, "--bottomonly"))) cfg->N *= 2;

  /* open CM file */
  if ((cfg->cmfp = CMFileOpen(cfg->cmfile, NULL)) == NULL)
    ESL_FAIL(eslFAIL, errbuf, "Failed to open covariance model save file %s\n", cfg->cmfile);

  /* optionally, open trace file */
  if (esl_opt_GetString(go, "--tfile") != NULL) {
    if ((cfg->tfp = fopen(esl_opt_GetString(go, "--tfile"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --tfile output file %s\n", esl_opt_GetString(go, "--tfile"));
    }

  /* optionally, open bands file */
  if (esl_opt_GetString(go, "--bfile") != NULL) {
    if ((cfg->bfp = fopen(esl_opt_GetString(go, "--bfile"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --bfile output file %s\n", esl_opt_GetString(go, "--bfile"));
    }

  /* optionally, open filter threshold calc histogram file */
  if (esl_opt_GetString(go, "--filhfile") != NULL) {
    if ((cfg->filhfp = fopen(esl_opt_GetString(go, "--filhfile"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --filhfile output file %s\n", esl_opt_GetString(go, "--filhfile"));
    }

  /* optionally, open filter threshold calc info file */
  if (esl_opt_GetString(go, "--filrfile") != NULL) {
    if ((cfg->filrfp = fopen(esl_opt_GetString(go, "--filrfile"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --filrfile output file %s\n", esl_opt_GetString(go, "--filrfile"));
    }

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
  char           errbuf[cmERRBUFSIZE];
  CM_t          *cm = NULL;
  CMConsensus_t *cons = NULL;     /* precalculated consensus info for display purposes */
  float          Smin;
  int            using_e_cutoff;
  int            rci;
  dbseq_t       *dbseq = NULL;
  int            do_top;


  if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);
  /*if ((status = init_shared_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);*/
  do_top = (cfg->init_rci == 0) ? TRUE : FALSE; 

  while (CMFileRead(cfg->cmfp, &(cfg->abc), &cm))
    {
      if (cm == NULL) cm_Fail("Failed to read CM from %s -- file corrupt?\n", cfg->cmfile);
      cfg->ncm++;

      /* initialize the flags/options/params and configuration of the CM */
      if((  status = initialize_cm(go, cfg, errbuf, cm))                    != eslOK) cm_Fail(errbuf);
      if((  status = update_avg_hit_len(go, cfg, errbuf, cm))               != eslOK) cm_Fail(errbuf);
      if((  status = CreateCMConsensus(cm, cfg->abc_out, 3.0, 1.0, &cons))  != eslOK) cm_Fail(errbuf);
      if(cm->flags & CMH_GUMBEL_STATS) 
	if((status = set_gumbels(go, cfg, errbuf, cm))                      != eslOK) cm_Fail(errbuf);
      if((status = set_window (go, cfg, errbuf, cm))                        != eslOK) cm_Fail(errbuf);
      if((status = set_searchinfo(go, cfg, errbuf, cm))                     != eslOK) cm_Fail(errbuf);
      if(esl_opt_GetBoolean(go, "--hmmcalcthr"))
	if((status = calc_filter_threshold(go, cfg, errbuf, cm, &Smin))     != eslOK) cm_Fail(errbuf);
      print_searchinfo(go, cfg, stdout, cm, cfg->N, errbuf);
      using_e_cutoff = (cm->si->cutoff_type[cm->si->nrounds] == E_CUTOFF) ? TRUE : FALSE;
	 
      while ((status = read_next_search_seq(cfg->abc, cfg->sqfp, cfg->do_rc, &dbseq)) == eslOK)
	{
	  for(rci = cfg->init_rci; rci <= cfg->do_rc; rci++) {
	    /*printf("SEARCHING >%s %d\n", dbseq->sq[reversed]->name, reversed);*/
	    if ((status = process_search_workunit(go, cfg, errbuf, cm, dbseq->sq[rci]->dsq, dbseq->sq[rci]->n, &dbseq->results[rci])) != eslOK) cm_Fail(errbuf);
	    remove_overlapping_hits(dbseq->results[rci], 1, dbseq->sq[rci]->n);
	    if(using_e_cutoff) remove_hits_over_e_cutoff(cm, cm->si, dbseq->results[rci], dbseq->sq[rci]); 
	  }
	  print_results (cm, cm->si, cfg->abc_out, cons, dbseq, do_top, cfg->do_rc);
	  for(rci = 0; rci <= cfg->do_rc; rci++) { /* we can free results for top strand even if cfg->init_rci is 1, due to --bottomonly */
	    FreeResults(dbseq->results[rci]);
	    esl_sq_Destroy(dbseq->sq[rci]);
	  }
	  free(dbseq);
	}
      if (status != eslEOF) cm_Fail("Parse failed, line %d, file %s:\n%s", 
				    cfg->sqfp->linenumber, cfg->sqfp->filename, cfg->sqfp->errbuf);
      FreeCM(cm);
      FreeCMConsensus(cons);
      esl_sqio_Rewind(cfg->sqfp); /* we may be searching this file again with another CM */
    }
}

#ifdef HAVE_MPI
/* mpi_master()
 * The MPI version of cmsearch
 * Follows standard pattern for a master/worker load-balanced MPI program 
 * (SRE notes J1/78-79).
 * 
 * A master can only return if it's successful. 
 * Errors in an MPI master come in two classes: recoverable and nonrecoverable.
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
static void
mpi_master(const ESL_GETOPTS *go, struct cfg_s *cfg)
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
  int      wi_error = 0;                /* worker index that sent back an error message, if an error occurs */

  CM_t *cm;
  CMConsensus_t *cons = NULL;     /* precalculated consensus info for display purposes */
  float          Smin;

  int si      = 0;
  int si_recv = 1;
  int *silist = NULL;

  int in_rc = FALSE;
  int *rclist = NULL;
  int rci;

  int seqpos = 1;
  int *seqposlist = NULL;

  int len;
  int *lenlist = NULL;

  int *sentlist = NULL;

  int ndbseq = 0;
  dbseq_t **dbseqlist    = NULL;
  dbseq_t *dbseq = NULL;
  
  char     errbuf[cmERRBUFSIZE];
  MPI_Status mpistatus; 
  int      n;

  int need_seq = TRUE;
  int chunksize;
  search_results_t *worker_results;


  /* TEMPORARY */
  if(esl_opt_GetBoolean(go, "--hmmcalcthr"))
    cm_Fail("ERROR, --hmmcalcthr and --mpi not yet implemented.");

  /* Master initialization: including, figure out the alphabet type.
   * If any failure occurs, delay printing error message until we've shut down workers.
   */
  if (xstatus == eslOK) { if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) xstatus = status; }
  /*if (xstatus == eslOK) { if ((status = init_shared_cfg(go, cfg, errbuf)) != eslOK) xstatus = status; }*/
  if (xstatus == eslOK) { bn = 4096; if ((buf = malloc(sizeof(char) * bn)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((silist     = malloc(sizeof(int) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((rclist     = malloc(sizeof(int) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((seqposlist = malloc(sizeof(int) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((lenlist    = malloc(sizeof(int) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((dbseqlist  = malloc(sizeof(dbseq_t *) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((sentlist   = malloc(sizeof(int) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }

  MPI_Bcast(&xstatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (xstatus != eslOK) cm_Fail(errbuf);
  ESL_DPRINTF1(("MPI master is initialized\n"));

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
  MPI_Bcast(&(cfg->N), 1, MPI_LONG, 0, MPI_COMM_WORLD);
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

  while (xstatus == eslOK && CMFileRead(cfg->cmfp, &(cfg->abc), &cm))
    {
      cfg->ncm++;  
      ESL_DPRINTF1(("MPI master read CM number %d\n", cfg->ncm));
      if((status = cm_master_MPIBcast(cm, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI broadcast CM failed.");
      
      /* initialize the flags/options/params of the CM */
      if((status   = initialize_cm(go, cfg, errbuf, cm))                    != eslOK) cm_Fail(errbuf);
      if((status   = update_avg_hit_len(go, cfg, errbuf, cm))               != eslOK) cm_Fail(errbuf);
      if((status   = CreateCMConsensus(cm, cfg->abc_out, 3.0, 1.0, &cons))  != eslOK) cm_Fail(errbuf);
      if(cm->flags & CMH_GUMBEL_STATS) 
	if((status = set_gumbels(go, cfg, errbuf, cm))                      != eslOK) cm_Fail(errbuf);
      if((status = set_window (go, cfg, errbuf, cm))                        != eslOK) cm_Fail(errbuf);
      if((status = set_searchinfo(go, cfg, errbuf, cm))                     != eslOK) cm_Fail(errbuf);
      if(esl_opt_GetBoolean(go, "--hmmcalcthr"))
	if((status = calc_filter_threshold(go, cfg, errbuf, cm, &Smin))     != eslOK) cm_Fail(errbuf);

      print_searchinfo(go, cfg, stdout, cm, cfg->N, errbuf);
      using_e_cutoff = (cm->si->cutoff_type[cm->si->nrounds] == E_CUTOFF) ? TRUE : FALSE;

      wi = 1;
      ndbseq = 0;
      while (have_work || nproc_working)
	{
	  if (need_seq) /* see mpifuncs.c:search_enqueue*/
	    {
	      need_seq = FALSE;
	      /* read a new seq */
	      if((status = read_next_search_seq(cfg->abc, cfg->sqfp, cfg->do_rc, &dbseq)) == eslOK) 
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
		  chunksize = determine_seq_chunksize(cfg, dbseq->sq[0]->n, cm->W);
		  ESL_DPRINTF1(("L: %d chunksize: %d\n", dbseq->sq[0]->n, chunksize));
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
		      if ((status = cm_search_results_MPIUnpack(buf, bn, &pos, MPI_COMM_WORLD, &worker_results)) != eslOK) cm_Fail("search results unpack failed");
		      ESL_DPRINTF1(("MPI master has unpacked search results\n"));
		      
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
			/* careful, dbseqlist[si_recv]->results[rclist[wi]] now points to the nodes in worker_results->data,
			 * don't free those (don't use FreeResults(worker_results)) */
			free(worker_results);
			worker_results = NULL;
		      }
		      dbseqlist[si_recv]->chunks_sent--;
		      if(sentlist[si_recv] && dbseqlist[si_recv]->chunks_sent == 0)
			{
			  for(rci = 0; rci <= cfg->do_rc; rci++) {
			    remove_overlapping_hits(dbseqlist[si_recv]->results[rci], 1, dbseqlist[si_recv]->sq[rci]->n);
			    if(using_e_cutoff) remove_hits_over_e_cutoff(cm, cm->si, dbseqlist[si_recv]->results[rci], dbseqlist[si_recv]->sq[rci]);
			    
			    
			  }					      
			  print_results(cm, cm->si, cfg->abc_out, cons, dbseqlist[si_recv], TRUE, cfg->do_rc);
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
      ESL_DPRINTF1(("MPI master: done with this CM. Telling all workers\n"));
      for (wi = 1; wi < cfg->nproc; wi++) 
	if ((status = cm_dsq_MPISend(NULL, 0, wi, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("Shutting down a worker failed.");
      FreeCM(cm);
      FreeCMConsensus(cons);
      esl_sqio_Rewind(cfg->sqfp); /* we may be searching this file again with another CM */
    }
  
  /* On success or recoverable errors:
   * Shut down workers cleanly. 
   */
  ESL_DPRINTF1(("MPI master is done. Shutting down all the workers cleanly\n"));
  if((status = cm_master_MPIBcast(NULL, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI broadcast CM failed.");
  free(buf);
  
  if (xstatus != eslOK) { fprintf(stderr, "Worker: %d had a problem.\n", wi_error); cm_Fail(errbuf); }
  else                  return;

 ERROR: 
  cm_Fail("memory allocation error.");
  return; /* NOTREACHED */
}


static void
mpi_worker(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int           xstatus = eslOK;
  int           status;
  /*int           type;*/
  CM_t         *cm  = NULL;
  char         *wbuf = NULL;	/* packed send/recv buffer  */
  int           wn   = 0;	/* allocation size for wbuf */
  int           sz, n;		/* size of a packed message */
  int           pos;
  char          errbuf[cmERRBUFSIZE];
  /*float         Smin;*/
  /*MPI_Status  mpistatus;*/
  ESL_DSQ      *dsq = NULL;
  int           L;
  search_results_t *results = NULL;

  /* After master initialization: master broadcasts its status.
   */
  MPI_Bcast(&xstatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (xstatus != eslOK) return; /* master saw an error code; workers do an immediate normal shutdown. */
  ESL_DPRINTF1(("worker %d: sees that master has initialized\n", cfg->my_rank));
  
  /* Master now broadcasts worker initialization information (db size N) 
   * Workers returns their status post-initialization.
   * Initial allocation of wbuf must be large enough to guarantee that
   * we can pack an error result into it, because after initialization,
   * errors will be returned as packed (code, errbuf) messages.
   */
  MPI_Bcast(&(cfg->N), 1, MPI_LONG, 0, MPI_COMM_WORLD);
  if (xstatus == eslOK) { wn = 4096;  if ((wbuf = malloc(wn * sizeof(char))) == NULL) xstatus = eslEMEM; }
  MPI_Reduce(&xstatus, &status, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD); /* everyone sends xstatus back to master */
  if (xstatus != eslOK) {
    if (wbuf != NULL) free(wbuf);
    return; /* shutdown; we passed the error back for the master to deal with. */
  }
  ESL_DPRINTF1(("worker %d: initialized N: %ld\n", cfg->my_rank, cfg->N));
  
  /* source = 0 (master); tag = 0 */
  while ((status = cm_worker_MPIBcast(0, MPI_COMM_WORLD, &wbuf, &wn, &(cfg->abc), &cm)) == eslOK)
    {
      ESL_DPRINTF1(("Worker %d succesfully received CM, num states: %d num nodes: %d\n", cfg->my_rank, cm->M, cm->nodes));
      
      /* initialize the flags/options/params of the CM */
      if((status   = initialize_cm(go, cfg, errbuf, cm))                    != eslOK) goto ERROR;
      if((status   = update_avg_hit_len(go, cfg, errbuf, cm))               != eslOK) goto ERROR;
      if(cm->flags & CMH_GUMBEL_STATS) 
	if((status = set_gumbels(go, cfg, errbuf, cm))                      != eslOK) goto ERROR;
      if((status = set_window(go, cfg, errbuf, cm))                         != eslOK) goto ERROR;
      if((status = set_searchinfo(go, cfg, errbuf, cm))                     != eslOK) goto ERROR;
      
      /* print_searchinfo(go, cfg, stdout, cm, cm_mode, cp9_mode, cfg->N, errbuf); */
      
      while((status = cm_dsq_MPIRecv(0, 0, MPI_COMM_WORLD, &wbuf, &wn, &dsq, &L)) == eslOK)
	{
	  ESL_DPRINTF1(("worker %d: has received search job, length: %d\n", cfg->my_rank, L));
	  if ((status = process_search_workunit(go, cfg, errbuf, cm, dsq, L, &results)) != eslOK) goto ERROR;
	  ESL_DPRINTF1(("worker %d: has gathered search results\n", cfg->my_rank));
	  
	  n = 0;
	  if (MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &sz) != 0) /* room for the status code */
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
	  if (cm_search_results_MPIPack(results, wbuf, wn, &pos, MPI_COMM_WORLD) != eslOK)
	    ESL_XFAIL(eslFAIL, errbuf, "cm_search_results_MPIPack() call failed"); 
	  MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);
	  ESL_DPRINTF1(("worker %d: has sent results to master in message of %d bytes\n", cfg->my_rank, pos));

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
  return;

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

  return;
}
#endif /*HAVE_MPI*/

/* A search work unit consists of a CM, digitized sequence dsq, and indices i and j.
 * The job is to search dsq from i..j and return search results in <*ret_results>.
 */
static int
process_search_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, ESL_DSQ *dsq, int L, search_results_t **ret_results)
{
  int status;
  search_results_t **results;
  int n;

  if(cm->si == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cm->si is NULL in process_search_workunit()\n");

  ESL_ALLOC(results, sizeof(search_results_t *) * (cm->si->nrounds+1));
  for(n = 0; n <= cm->si->nrounds; n++) results[n] = CreateResults(INIT_RESULTS);

  if((status = DispatchSearch(cm, errbuf, 0, dsq, 1, L, results, NULL, NULL)) != eslOK) goto ERROR;

  /* we only care about the final results, that survived all the rounds (all the filtering rounds plus the final round) */
  *ret_results = results[cm->si->nrounds];
  /* free the results describing what survived each round of filtering (if any) */
  for(n = 0; n < cm->si->nrounds; n++) FreeResults(results[n]);
  free(results);

  return eslOK;
  
 ERROR:
  ESL_DPRINTF1(("worker %d: has caught an error in process_search_workunit\n", cfg->my_rank));
  FreeCM(cm);
  return status;
}

/* A CP9 filter work unit consists of a CM and an int (nseq).
 * The job is to emit nseq sequences with a score better than cutoff (rejecting
 * those that are worse), and then search those seqs with a CP9, returning the scores of the
 * best CP9 hit within each sequence.
 */
static int
process_cp9filter_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int nseq)
{
  /*int status;*/
  cm_Fail("WRITE process_cp9filter_workunit()");
  return eslOK;
  
  /* ERROR:
  ESL_DPRINTF1(("worker %d: has caught an error in process_cp9filter_workunit\n", cfg->my_rank));
  return status;*/
}

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
  int nstarts, nexits, nd;
  
  /* set up CM parameters that are option-changeable */
  cm->beta   = esl_opt_GetReal(go, "--beta"); /* this will be DEFAULT_BETA unless changed at command line */
  cm->tau    = esl_opt_GetReal(go, "--tau");  /* this will be DEFAULT_TAU unless changed at command line */

  use_hmmonly = (esl_opt_GetBoolean(go, "--hmmviterbi") || esl_opt_GetBoolean(go, "--hmmforward")) ? TRUE : FALSE;

  /* Update cm->config_opts and cm->align_opts based on command line options */

  /* config_opts */
  if(! esl_opt_GetBoolean(go, "-g"))            cm->config_opts |= CM_CONFIG_LOCAL;
  if(! esl_opt_GetBoolean(go, "--hmmglocal"))   cm->config_opts |= CM_CONFIG_HMMLOCAL;
  if(! esl_opt_GetBoolean(go, "--hmmnoel"))     cm->config_opts |= CM_CONFIG_HMMEL;
  if(! esl_opt_GetBoolean(go, "--iins"))        cm->config_opts |= CM_CONFIG_ZEROINSERTS;
  /* config QDB? Yes, unless --noqdb || --hmmviterbi || --hmmforward enabled */
  if(! (esl_opt_GetBoolean(go, "--noqdb") || (use_hmmonly)))
    cm->config_opts |= CM_CONFIG_QDB;
  /* are we enforcing a subseq? */
  if(! esl_opt_IsDefault (go, "--enfseq")) {
    cm->config_opts |= CM_CONFIG_ENFORCE;
    cm_Fail("EPN, Mon Dec  3 13:14:43 2007, if you want to keep --enfseq, you'll have to add a HMM filter here somehow.");
    if((! use_hmmonly) && (! esl_opt_GetBoolean(go, "--enfnohmm")))
      /* We want to filter with special enforced CP9 HMM for the enforced subseq */
      ;//cm->search_opts |= CM_SEARCH_HMMFILTER;
    cm->enf_start = EnforceFindEnfStart(cm, esl_opt_GetInteger(go, "--enfstart"));
    /* --enfstart MUST have been enabled, --enfseq requires it */
    cm->enf_seq = esl_opt_GetString(go, "--enfseq");
  }

  /* search_opts */
  if(  esl_opt_GetBoolean(go, "--inside"))      cm->search_opts |= CM_SEARCH_INSIDE;
  if(  esl_opt_GetBoolean(go, "--toponly"))     cm->search_opts |= CM_SEARCH_TOPONLY;
  if(  esl_opt_GetBoolean(go, "--noalign"))     cm->search_opts |= CM_SEARCH_NOALIGN;
  if(  esl_opt_GetBoolean(go, "--null2"))       cm->search_opts |= CM_SEARCH_NULL2;
  if(  esl_opt_GetBoolean(go, "--greedy"))      cm->search_opts |= CM_SEARCH_CMGREEDY;
  if(  esl_opt_GetBoolean(go, "--hmmgreedy"))   cm->search_opts |= CM_SEARCH_HMMGREEDY;
  if(  esl_opt_GetBoolean(go, "--noqdb"))       cm->search_opts |= CM_SEARCH_NOQDB;
  if(  esl_opt_GetBoolean(go, "--hbanded"))     cm->search_opts |= CM_SEARCH_HBANDED;
  if(  esl_opt_GetBoolean(go, "--scan2bands"))  cm->search_opts |= CM_SEARCH_HMMSCANBANDS;
  if(  esl_opt_GetBoolean(go, "--sums"))        cm->search_opts |= CM_SEARCH_SUMS;
  if(  esl_opt_GetBoolean(go, "--hmmviterbi"))  { 
    cm->search_opts |= CM_SEARCH_HMMVITERBI;
    cm->search_opts |= CM_SEARCH_NOQDB;
  }
  if(  esl_opt_GetBoolean(go, "--hmmforward"))  { 
    cm->search_opts |= CM_SEARCH_HMMFORWARD;
    cm->search_opts |= CM_SEARCH_NOQDB;
  }

  /* align_opts */
  cm->align_opts |= CM_ALIGN_HBANDED;
  if(esl_opt_GetBoolean(go, "--optacc"))        cm->align_opts |= CM_ALIGN_OPTACC;
  if(esl_opt_GetBoolean(go, "--post"))          cm->align_opts |= CM_ALIGN_POST;

  /* flags */
  if(  esl_opt_GetBoolean(go, "--rtrans"))      cm->flags       |= CM_RSEARCHTRANS;

  /* read in QDBs if nec */
  if(! esl_opt_IsDefault(go, "--qdbfile"))
    {
      /* can't be in MPI mode, --mpi is incompatible with --qdbfile */
      FILE *qdb_fp;
      if(cfg->preset_dmin != NULL || cfg->preset_dmax != NULL)
	cm_Fail("ERROR: trying to read QDBs for more than one CM. With --qdbfile, <cm file> must have exactly 1 CM in it.");
      ESL_ALLOC(cfg->preset_dmin, sizeof(int) * cm->M);
      ESL_ALLOC(cfg->preset_dmax, sizeof(int) * cm->M);
      if ((qdb_fp = fopen(esl_opt_GetString(go, "--qdbfile"), "r")) == NULL)
	cm_Fail("failed to open QDB file %s", esl_opt_GetString(go, "--qdbfile"));
      if(!(read_qdb_file(qdb_fp, cm, cfg->preset_dmin, cfg->preset_dmax)))
	cm_Fail("ERROR reading QDB file: %s.\nDoes it correspond (same number of states) to this model?\n", esl_opt_GetString(go, "--qdbfile"));
      fclose(qdb_fp);
    }

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
  if(! esl_opt_IsDefault(go, "--pfend")) {
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
  /* finally, configure the CM for alignment based on cm->config_opts and cm->align_opts.
   * set local mode, make cp9 HMM, calculate QD bands etc. 
   */
  ConfigCM(cm, cfg->preset_dmin, cfg->preset_dmax); /* preset_d* usually NULL, unless --qdbfile */
  if(cm->config_opts & CM_CONFIG_ENFORCE) ConfigCMEnforce(cm);

  ESL_DPRINTF1(("cm->pbegin: %.3f\n", cm->pbegin));
  ESL_DPRINTF1(("cm->pend: %.3f\n", cm->pend));
  /* print qdbs to file if nec */
  if(! esl_opt_IsDefault(go, "--bfile")) {
    fprintf(cfg->bfp, "beta:%f\n", cm->beta);
    debug_print_bands(cfg->bfp, cm, cm->dmin, cm->dmax);
    fprintf(cfg->bfp, "beta:%f\n", cm->beta);
  }

  if(cfg->my_rank == 0) printf("CM %d: %s\n", (cfg->ncm), cm->name);
  return eslOK;
  
 ERROR:
  return status;
}


/* Function: update_avg_hit_len()
 * Date:     EPN, Sun Dec  9 15:50:39 2007
 * 
 * Purpose:  Calculate the average subseq length rooted at each state
 *           using the QDB calculation.
 *
 * Returns:  eslOK on success;
 */
int
update_avg_hit_len(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int safe_windowlen;
  float *avglen = NULL;

  safe_windowlen = cm->W * 2;
  while(!(BandCalculationEngine(cm, safe_windowlen, 1E-15, TRUE, NULL, NULL, NULL, &avglen))) {
    safe_windowlen *= 2;
    if(safe_windowlen > (cm->clen * 1000)) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "update_avg_hit_len(), band calculation safe_windowlen big: %d\n", safe_windowlen);
  }
  if(cfg->avglen != NULL) free(cfg->avglen);
  cfg->avglen = avglen;
  return eslOK;
}

/* set_gumbels()
 * Setup Gumbel distribution parameters for CM and HMM based on DB size.
 */
static int
set_gumbels(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  double tmp_K;                 /* for converting mu from cmfile to mu for N*/
  int i, p;

  /* Determine K from mu, lambda, L, then set CM mu for N */
  for(i = 0; i < GUM_NMODES; i++)
    for(p = 0; p < cm->stats->np; p++)
      {
	tmp_K = exp(cm->stats->gumAA[i][p]->mu * cm->stats->gumAA[i][p]->lambda) / 
	  cm->stats->gumAA[i][p]->L;
	cm->stats->gumAA[i][p]->mu = log(tmp_K * ((double) cfg->N)) /
	  cm->stats->gumAA[i][p]->lambda;
	cm->stats->gumAA[i][p]->L = cfg->N; /* update L, the seq size the stats correspond to */
      }
  ESL_DPRINTF1(("CM/CP9 statistics read from CM file\n"));
  if (cm->stats->np == 1) 
    ESL_DPRINTF1(("No partition points\n"));
  else {
    ESL_DPRINTF1(("Partition points are: "));
    for (p=0; p < cm->stats->np; p++)
      ESL_DPRINTF1(("%d %d..%d", p, cm->stats->ps[p], cm->stats->pe[p]));
  }

  return eslOK;
}


/* set_searchinfo()
 * Determine how many rounds of searching we will do (all rounds but last
 * round are filters), and set the relevant info in the SearchInfo_t <cm->si>
 * object, including cutoffs.
 */
static int
set_searchinfo(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int status;
  int n;
  int stype;
  int add_cyk_filter     = FALSE;
  int add_inside_filter  = FALSE;
  int add_viterbi_filter = FALSE;
  int add_forward_filter = FALSE;
  int search_opts;
  int use_hmmonly;
  int fthr_mode;
  int cutoff_type;
  float sc_cutoff = -1.;
  float e_cutoff = -1.;
  int  *dmin, *dmax; /* these become QDBs if we add a CM_FILTER */
  ScanMatrix_t *fsmx; 
  int safe_windowlen;
  float surv_fract;
  int cm_mode, cp9_mode;

  if(cm->si != NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "set_searchinfo(), cm->si is not NULL, shouldn't happen.\n");

  /* Create SearchInfo, specifying no filtering, we change the threshold below */
  CreateSearchInfo(cm, SCORE_CUTOFF, 0., -1.);
  if(cm->si == NULL) cm_Fail("set_searchinfo(), CreateSearchInfo() call failed.");
  SearchInfo_t *si = cm->si; 
  if(si->nrounds > 0) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_search_info(), si->nrounds (%d) > 0\n", si->nrounds);

  /*************************************************************************************
   * Filter related options:
   *
   * User can specify 0 to 2 rounds of filtering and cutoffs on the command line, 
   * 0 or 1 rounds can be CM  filters, with --fcm,  --fcmT,  and --fcmE (req's Gumbels)
   *                                        --fcminside specifies use inside, not CYK
   * 0 or 1 rounds can by HMM filters, with --fhmm, --fhmmT, and --fhmmE (req's Gumbels)
   *                                        --fhmmforward specifies use forward, not viterbi
   * Or user can specify that an HMM or hybrid filter as described in the the CM file 
   * be used with option --fgiven. --fgiven is incompatible with --fhmmviterbi and --fhmmforward
   * but not with --fcyk and --finside
   *
   *************************************************************************************
   * Final round related options (after all filtering is complete):
   *
   * --cyk:        search with CM CYK (TRUE by default)
   * --inside:     search with CM inside 
   * -T:           CM bit score threshold
   * -E:           CM E-value threshold (requires Gumbel info in CM file)
   * --ga:         use Rfam gathering threshold (bit sc) from CM file
   * --tc:         use Rfam trusted cutoff      (bit sc) from CM file
   * --nc:         use Rfam noise cutoff        (bit sc) from CM file
   *
   * --hmmviterbi: search with HMM viterbi
   * --hmmforward: search with HMM forward
   * --hmmT:       bit score threshold for --hmmviterbi or --hmmforward
   * --hmmE:       E-value threshold (requires Gumbel info in CM file)
   *
   *************************************************************************************
   */
  
  /* First, set up cutoff for final round, this will be round 0, unless filter info was read from the CM file */
  n           = si->nrounds;
  stype       = si->stype[n];
  search_opts = si->search_opts[n];

  /* determine configuration of CM and CP9 based on cm->flags & cm->search_opts */
  CM2Gumbel_mode(cm, search_opts, &cm_mode, &cp9_mode); 

  use_hmmonly = ((search_opts & CM_SEARCH_HMMVITERBI) || (search_opts & CM_SEARCH_HMMFORWARD));
  if(! use_hmmonly) {
    if(stype != SEARCH_WITH_CM) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_searchinfo(), search_opts for final round of search does not have HMMVITERBI or HMMFORWARD flags raised, but is not of type SEARCH_WITH_CM.");
    /* set up CM cutoff, either 0 or 1 of 6 options is enabled. 
     * esl_opt_IsDefault() returns FALSE even if option is enabled with default value 
     * We will NOT use this if --hmmviterbi
     */
    if(esl_opt_IsDefault(go, "-E") && 
       esl_opt_IsDefault(go, "-T") && 
       esl_opt_IsDefault(go, "--ga") && 
       esl_opt_IsDefault(go, "--tc") && 
       esl_opt_IsDefault(go, "--nc")) { 
      /* Choose from, in order of priority:
       * 1. default CM E value if CM file has Gumbel stats
       * 3. default CM bit score
       */
      if(cm->flags & CMH_GUMBEL_STATS) { /* use default CM E-value cutoff */
	cutoff_type = E_CUTOFF;
	e_cutoff    = esl_opt_GetReal(go, "-E");
	if((status = E2Score(cm, errbuf, cm_mode, e_cutoff, &sc_cutoff)) != eslOK) return status;
      }
      else { /* no Gumbel stats in CM file, use default bit score cutoff */
	cutoff_type = SCORE_CUTOFF;
	sc_cutoff   = esl_opt_GetReal(go, "-T");
	e_cutoff    = -1.; /* invalid, we'll never use it */  
      }
    }
    else if(! esl_opt_IsDefault(go, "-E")) {
      if(! (cm->flags & CMH_GUMBEL_STATS))
	ESL_FAIL(eslEINVAL, errbuf, "-E requires Gumbel statistics in <cm file>. Use cmcalibrate to get Gumbel stats.");
      cutoff_type = E_CUTOFF;
      e_cutoff    = esl_opt_GetReal(go, "-E");
      if((status = E2Score(cm, errbuf, cm_mode, e_cutoff, &sc_cutoff)) != eslOK) return status;
    }
    else if(! esl_opt_IsDefault(go, "-T")) {
      cutoff_type = SCORE_CUTOFF;
      sc_cutoff   = esl_opt_GetReal(go, "-T");
      e_cutoff    = -1.; /* invalid, we'll never use it */  
      if((sc_cutoff < 0.) && (! esl_opt_GetBoolean(go, "--greedy"))) ESL_FAIL(eslEINVAL, errbuf, "with -T <x> option, <x> can only be less than 0. if --greedy also enabled.");
    }
    else if(! esl_opt_IsDefault(go, "--ga")) {
      if(! (cm->flags & CMH_GA))
	ESL_FAIL(eslEINVAL, errbuf, "No GA gathering threshold in CM file, can't use --ga.");
      cutoff_type = SCORE_CUTOFF;
      sc_cutoff   = esl_opt_GetReal(go, "--ga");
      e_cutoff    = -1.; /* we'll never use it */
    }
    else if(! esl_opt_IsDefault(go, "--tc")) {
      if(! (cm->flags & CMH_TC))
	ESL_FAIL(eslEINVAL, errbuf, "No TC trusted cutoff in CM file, can't use --tc.");
      cutoff_type = SCORE_CUTOFF;
      sc_cutoff   = esl_opt_GetReal(go, "--tc");
      e_cutoff    = -1.; /* we'll never use it */
    }
    else if(! esl_opt_IsDefault(go, "--nc")) {
      if(! (cm->flags & CMH_NC))
	ESL_FAIL(eslEINVAL, errbuf, "No NC noise cutoff in CM file, can't use --nc.");
      cutoff_type = SCORE_CUTOFF;
      sc_cutoff   = esl_opt_GetReal(go, "--nc");
      e_cutoff    = -1.; /* we'll never use it */
    }
    else ESL_FAIL(eslEINCONCEIVABLE, errbuf, "No CM cutoff selected. This shouldn't happen.");
  } /* end of if(! use_hmmonly) */
  else { 
    if(stype != SEARCH_WITH_HMM) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "search_opts for final round of search has HMMVITERBI or HMMFORWARD flags raised, but is not of type SEARCH_WITH_HMM.");
    /* Set up CP9 HMM cutoff, either 0 or 1 of 2 options is enabled 
     * esl_opt_IsDefault() returns FALSE even if option is enabled with default value 
     */
    if(esl_opt_IsDefault(go, "--hmmE") && 
       esl_opt_IsDefault(go, "--hmmT")) {
      /* Choose from, in order of priority:
       * 1. default CP9 E value if CM file has Gumbel stats
       * 2. default CP9 bit score
       */
      if(cm->flags & CMH_GUMBEL_STATS) { /* use default CP9 E-value cutoff */
	cutoff_type = E_CUTOFF;
	e_cutoff    = esl_opt_GetReal(go, "--hmmE"); 
	if((status = E2Score(cm, errbuf, cp9_mode, e_cutoff, &sc_cutoff)) != eslOK) return status;
      }
      else { /* no Gumbel stats in CM file, use default bit score cutoff */
	cutoff_type = SCORE_CUTOFF;
	sc_cutoff   = esl_opt_GetReal(go, "--hmmT");
	e_cutoff    = -1; /* we'll never use it */
      }
    }
    else if(! esl_opt_IsDefault(go, "--hmmE")) {
      if(! (cm->flags & CMH_GUMBEL_STATS))
	ESL_FAIL(eslEINVAL, errbuf, "--hmmE requires Gumbel statistics in <cm file>. Use cmcalibrate to get Gumbel stats.");
      cutoff_type = E_CUTOFF;
      e_cutoff    = esl_opt_GetReal(go, "--hmmE");
      if((status  = E2Score(cm, errbuf, cp9_mode, e_cutoff, &sc_cutoff)) != eslOK) return status; 
    }
    else if(! esl_opt_IsDefault(go, "--hmmT")) {
      cutoff_type = SCORE_CUTOFF;
      sc_cutoff   = esl_opt_GetReal(go, "--hmmT");
      e_cutoff    = -1.; /* we'll never use this */
      if((sc_cutoff < 0.) && (! esl_opt_GetBoolean(go, "--hmmgreedy"))) ESL_FAIL(eslEINVAL, errbuf, "with --hmmT <x> option, <x> can only be less than 0. if --hmmgreedy also enabled.");
    }
  }
  /* update the search info, which holds the thresholds */
  UpdateSearchInfoCutoff(cm, cm->si->nrounds, cutoff_type, sc_cutoff, e_cutoff);   
  ValidateSearchInfo(cm, cm->si);
  /* DumpSearchInfo(cm->si); */
  /* done with threshold for final round */

  /* Set up the filters and their thresholds 
   * 1. add a CM  filter, if necessary
   * 2. add a HMM filter, if necessary
   */

  /* CM filter */
  add_cyk_filter    = esl_opt_GetBoolean(go, "--fcyk");
  add_inside_filter = esl_opt_GetBoolean(go, "--finside");
  ESL_DASSERT1((!(add_cyk_filter && add_inside_filter))); /* should be enforced by getopts */
  if(add_cyk_filter || add_inside_filter) { /* determine thresholds for filters */
    /* set up CM filter cutoff, either 0 or 1 of 2 options is enabled. 
     * esl_opt_IsDefault() returns FALSE even if option is enabled with default value 
     */
    if(esl_opt_IsDefault(go, "--fE") && 
       esl_opt_IsDefault(go, "--fT")) {
      /* Choose from, in order of priority:
       * 1. default CM filter E value if CM file has Gumbel stats
       * 3. default CM filter bit score
       */
      if(cm->flags & CMH_GUMBEL_STATS) { /* use default CM E-value cutoff */
	cutoff_type = E_CUTOFF;
	e_cutoff    = esl_opt_GetReal(go, "--fE");
	if((status  = E2Score(cm, errbuf, cm_mode, e_cutoff, &sc_cutoff)) != eslOK) return status;
      }
      else { /* no Gumbel stats in CM file, use default bit score cutoff */
	cutoff_type = SCORE_CUTOFF;
	sc_cutoff   = esl_opt_GetReal(go, "--fT");
	e_cutoff    = -1.; /* invalid, we'll never use it */  
      }
    }
    else if(! esl_opt_IsDefault(go, "--fE")) {
      if(! (cm->flags & CMH_GUMBEL_STATS))
	ESL_FAIL(eslEINVAL, errbuf, "--fE requires Gumbel statistics in <cm file>. Use cmcalibrate to get Gumbel stats.");
      cutoff_type = E_CUTOFF;
      e_cutoff    = esl_opt_GetReal(go, "--fE");
      if((status  = E2Score(cm, errbuf, cm_mode, e_cutoff, &sc_cutoff)) != eslOK) return status;
    }
    else if(! esl_opt_IsDefault(go, "--fT")) {
      cutoff_type = SCORE_CUTOFF;
      sc_cutoff   = esl_opt_GetReal(go, "--fT");
      if((sc_cutoff < 0.) && (! esl_opt_GetBoolean(go, "--fgreedy"))) ESL_FAIL(eslEINVAL, errbuf, "with --fT <x> option, <x> can only be less than 0. if --fgreedy also enabled.");
      e_cutoff    = -1.; /* we'll never use it */
    }
    else ESL_FAIL(eslEINCONCEIVABLE, errbuf, "No CM filter cutoff selected. This shouldn't happen.");
    
    /* build the ScanMatrix_t for this round, requires calcing dmin, dmax */
    safe_windowlen = cm->W * 2;
    while(!(BandCalculationEngine(cm, safe_windowlen, esl_opt_GetReal(go, "--fbeta"), FALSE, &dmin, &dmax, NULL, NULL))) {
      free(dmin);
      free(dmax);
      dmin = NULL;
      dmax = NULL;
      safe_windowlen *= 2;
      if(safe_windowlen > (cm->clen * 1000)) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_searchinfo(), band calculation safe_windowlen big: %d\n", safe_windowlen);
    }
    fsmx = cm_CreateScanMatrix(cm, dmax[0], dmin, dmax, esl_opt_GetReal(go, "--fbeta"), TRUE, add_cyk_filter, add_inside_filter);
    /* add the filter */
    AddFilterToSearchInfo(cm, add_cyk_filter, add_inside_filter, FALSE, FALSE, FALSE, fsmx, NULL, cutoff_type, sc_cutoff, e_cutoff);
    ValidateSearchInfo(cm, cm->si);
    /* DumpSearchInfo(cm->si); */
  }
  else if (! esl_opt_IsDefault(go, "--fbeta")) ESL_FAIL(eslEINCOMPAT, errbuf, "--fbeta has an effect with --fcyk or --finside");

  /* HMM filter */
  /* if --fgiven was enabled, --fhmmviterbi, --fhmmforward could NOT have been selected, 
   * so we won't enter any of the loops below.
   */
  add_viterbi_filter = esl_opt_GetBoolean(go, "--fhmmviterbi");
  add_forward_filter = esl_opt_GetBoolean(go, "--fhmmforward");
  ESL_DASSERT1((!(add_viterbi_filter && add_forward_filter))); /* should be enforced by getopts */

  if(esl_opt_GetBoolean(go, "--fgiven")) {
    /* determine filter threshold mode, the mode of final stage of searching, either FTHR_CM_LC,
     * FTHR_CM_LI, FTHR_CM_GC, FTHR_CM_GI, (can't be an HMM mode b/c getopts enforces --fgiven incompatible with
     * --hmmviterbi and --hmmforward). 
     */
    if((status = CM2FthrMode(cm, errbuf, cm->search_opts, &fthr_mode)) != eslOK) return status;
    if(!(cm->flags & CMH_FILTER_STATS))              ESL_FAIL(eslEINCOMPAT, errbuf,      "set_searchinfo(), --fgiven enabled, but cm's CMH_FILTER_STATS flag is down.");
    if(cm->stats->bfA[fthr_mode]->is_valid == FALSE) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_searchinfo(), --fgiven enabled, cm's CMH_FILTER_STATS is raised, but best filter info for fthr_mode %d is invalid.", fthr_mode);
    ESL_DASSERT1((!add_viterbi_filter)); /* --fhmmviterbi is incompatible with --fgiven, enforced by getopts */
    ESL_DASSERT1((!add_forward_filter)); /* --fhmmforward is incompatible with --fgiven, enforced by getopts */
    cutoff_type = E_CUTOFF;
    e_cutoff    = cm->stats->bfA[fthr_mode]->e_cutoff * ((double) cfg->N / (double) cm->stats->bfA[fthr_mode]->db_size); 
    /* check if --hmmmaxE applies */
    if((! esl_opt_IsDefault(go, "--hmmmaxE")) && (cm->stats->bfA[fthr_mode]->ftype == FILTER_WITH_HMM_VITERBI || cm->stats->bfA[fthr_mode]->ftype == FILTER_WITH_HMM_FORWARD)) {
      e_cutoff = ESL_MIN(e_cutoff, esl_opt_GetReal(go, "--hmmmaxE"));
    }
    if((status  = E2Score(cm, errbuf, cm_mode, e_cutoff, &sc_cutoff)) != eslOK) return status; /* note: use cm_mode, not fthr_mode */

    /* Predict survival fraction from filter based on E-value, assume average hit length is cfg->avglen[0] (from QD band calc) */
    surv_fract = (e_cutoff * ((2. * cm->W) - cfg->avglen[0])) / ((double) cfg->N); 
    if(surv_fract < 0.9) { /* else filtering is not worth it */
      if(cm->stats->bfA[fthr_mode]->ftype == FILTER_WITH_HMM_VITERBI) add_viterbi_filter = TRUE;
      if(cm->stats->bfA[fthr_mode]->ftype == FILTER_WITH_HMM_FORWARD) add_forward_filter = TRUE;
      if(cm->stats->bfA[fthr_mode]->ftype == FILTER_WITH_HYBRID)      /*add_hybrid_filter = TRUE;*/
	ESL_FAIL(eslEINCOMPAT, errbuf, "set_searchinfo(), --fgiven enabled and hybrid filter set in CM file, but we can't deal with this yet.\n");
    }
  }
  else if(add_viterbi_filter || add_forward_filter) { /* determine thresholds for filters */
    /* Set up HMM cutoff, either 0 or 1 of 3 options is enabled */
    ESL_DASSERT1((! use_hmmonly)); /* should be enforced by getopts */
    if(esl_opt_IsDefault(go, "--hmmcalcthr") && 
       esl_opt_IsDefault(go, "--hmmE") && 
       esl_opt_IsDefault(go, "--hmmT")) {
      /* Choose from, in order of priority:
       * 1. default CP9 E value if CM file has Gumbel stats
       * 2. default CP9 bit score
       */
      if(cm->flags & CMH_GUMBEL_STATS) { /* use default CM E-value cutoff */
	cutoff_type = E_CUTOFF;
	e_cutoff    = esl_opt_GetReal(go, "--hmmE");
	if((status  = E2Score(cm, errbuf, cp9_mode, e_cutoff, &sc_cutoff)) != eslOK) return status; 
      }
      else { /* no Gumbel stats in CM file, use default bit score cutoff */
	cutoff_type = SCORE_CUTOFF;
	sc_cutoff   = esl_opt_GetReal(go, "--hmmT");
	e_cutoff    = -1.; /* we'll never use this */
      }
    }
    else if(! esl_opt_IsDefault(go, "--hmmcalcthr")) {
      if(! (cm->flags & CMH_GUMBEL_STATS))
	ESL_FAIL(eslEINVAL, errbuf, "--hmmcalcthr requires Gumbel statistics in <cm file>. Use cmcalibrate to get Gumbel stats.");
      cutoff_type = E_CUTOFF;
      /* this gets overwritten later after threshold is calculated */
      e_cutoff = esl_opt_GetReal(go, "--hmmE");
      if((status  = E2Score(cm, errbuf, cp9_mode, e_cutoff, &sc_cutoff)) != eslOK) return status; 
    }
    else if(! esl_opt_IsDefault(go, "--hmmE")) {
      if(! (cm->flags & CMH_GUMBEL_STATS))
	ESL_FAIL(eslEINVAL, errbuf, "--hmmE requires Gumbel statistics in <cm file>. Use cmcalibrate to get Gumbel stats.");
      cutoff_type = E_CUTOFF;
      e_cutoff    = esl_opt_GetReal(go, "--hmmE");
      if((status  = E2Score(cm, errbuf, cp9_mode, e_cutoff, &sc_cutoff)) != eslOK) return status; 
    }
    else if(! esl_opt_IsDefault(go, "--hmmT")) {
      cutoff_type = SCORE_CUTOFF;
      sc_cutoff   = esl_opt_GetReal(go, "--hmmT");
      e_cutoff    = -1.; /* we'll never use this */
      if((sc_cutoff < 0.) && (! esl_opt_GetBoolean(go, "--hmmgreedy"))) ESL_FAIL(eslEINVAL, errbuf, "with --hmmT <x> option, <x> can only be less than 0. if --hmmgreedy also enabled.");
    }
    else ESL_FAIL(eslEINCONCEIVABLE, errbuf, "No HMM filter cutoff selected. This shouldn't happen.");
  }
  if(add_viterbi_filter || add_forward_filter) {
    /* add the filter */
    AddFilterToSearchInfo(cm, FALSE, FALSE, add_viterbi_filter, add_forward_filter, FALSE, NULL, NULL, cutoff_type, sc_cutoff, e_cutoff);
    ValidateSearchInfo(cm, cm->si);
    /*DumpSearchInfo(cm->si); */
  }

  return eslOK;
}

/* set_window()
 * Set cm->W, the window size for scanning.
 *
 * 1. cm->W is overwritten here if --window enabled.
 * 2. cm->W is set to dmax[0] if --noqdb, --hmmviterbi, or --hmmforward enabled after calc'ing QDBs for 
 *    sole purpose of determining cm->W. 
 * 3. else cm->W was to cm->dmax[0] in ConfigCM()'s call to ConfigQDB(), 
 *    which is what it should be.
 */
static int
set_window(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int use_hmmonly;
  int do_float;
  int do_int;

  use_hmmonly = (esl_opt_GetBoolean(go, "--hmmviterbi") || esl_opt_GetBoolean(go, "--hmmforward")) ? TRUE : FALSE;

  if(! esl_opt_IsDefault(go, "--window")) {
    if((! esl_opt_GetBoolean(go, "--noqdb")) && (! use_hmmonly))
      ESL_FAIL(eslEINCOMPAT, errbuf, "--window only makes sense with --noqdb, --hmmviterbi, or --hmmforward enabled. Use smaller --beta values to decrease window size.\n");
    cm->W = esl_opt_GetInteger(go, "--window");
  }
  else if(esl_opt_GetBoolean(go, "--noqdb") || use_hmmonly) {
    if(cm->dmin != NULL || cm->dmax != NULL) 
      ESL_FAIL(eslEINCONCEIVABLE, errbuf, "--hmmviterbi, --hmmforward or --noqdb enabled, but cm->dmin and cm->dmax non-null. This shouldn't happen.");
    int *dmin;
    int *dmax;
    int safe_windowlen = cm->clen * 2;
    while(!(BandCalculationEngine(cm, safe_windowlen, cm->beta, 0, &(dmin), &(dmax), NULL, NULL)))
      {
	free(dmin);
	free(dmax);
	safe_windowlen *= 2;
	if(safe_windowlen > (cm->clen * 1000))
	  ESL_FAIL(eslEINVAL, errbuf, "ERROR in set_window, safe_windowlen big: %d\n", safe_windowlen);
      }
    cm->W = dmax[0];
    free(dmin);
    free(dmax);
    CMLogoddsify(cm); /* QDB calculation invalidates log odds scores */
  }

  /* Setup ScanMatrix for CYK/Inside scanning functions, we can't 
   * do it in initialize_cm(), b/c it's W dependent; W was just set.
   * We don't need it if we're only using an HMM though.
   */
  if(use_hmmonly) cm->smx = NULL;
  else { 
    do_float = TRUE;
    do_int   = FALSE;
    if(cm->search_opts & CM_SEARCH_INSIDE) { do_float = FALSE; do_int = TRUE; }
    cm_CreateScanMatrixForCM(cm, do_float, do_int);
    if(cm->smx == NULL) cm_Fail("set_window(), use_hmmonly is FALSE, CreateScanMatrixForCM() call failed, mx is NULL.");
  }

  return eslOK;
}

/* calc_filter_threshold()
 * Calculate the filter threshold for the CP9 HMM.
 */
static int
calc_filter_threshold(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float *ret_Smin)
{
  *ret_Smin = 0.;
  return eslOK;
}

/* read_qdb_file()
 * Read QDBs from a file outputted from cmbuild. Only useful for testing/debugging,
 */
static int  
read_qdb_file(FILE *fp, CM_t *cm, int *dmin, int *dmax)
{
  int     status;
  char   *buf;
  int     n;			/* length of buf */
  char   *s;
  int     M;			/* number of states in model */
  int     v;		        /* counter for states */
  char   *tok;
  int     toklen;
  int     read_v;

  /* format of QDB file: 
   * line  1        :<cm->M>
   * lines 2 -> M+1 :<v> <dmin> <dmax> */

  buf = NULL;
  n   = 0;
  if (feof(fp) || (status = esl_fgets(&buf, &n, fp)) != eslOK) goto ERROR;

  s   = buf;
  if ((status = esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto ERROR;
  if (! is_integer(tok))                                    goto ERROR;
  M = atoi(tok);
  if(M != cm->M) goto ERROR;

  v = 0;
  while ((status = esl_fgets(&buf, &n, fp)) == eslOK) 
    {
      if (strncmp(buf, "//", 2) == 0) 
	break;
      s   = buf;
      if ((status = esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto ERROR;      
      if (! is_integer(tok)) { status = eslEINVAL;                    goto ERROR; }
      read_v = atoi(tok);
      if(v != read_v) goto ERROR;

      if ((status = esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto ERROR;      
      if (! is_integer(tok)) { status = eslEINVAL;                    goto ERROR; }
      dmin[v] = atoi(tok);

      if ((status = esl_strtok(&s, " \t\n", &tok, &toklen)) != eslOK) goto ERROR;      
      if (! is_integer(tok)) {                                        goto ERROR; }
      dmax[v] = atoi(tok);

      v++;
    }
  if(v != M) { status = eslEINVAL; goto ERROR; }
  if(status != eslOK) goto ERROR;

  if (buf != NULL) free(buf);
  return eslOK;

 ERROR:
  if (cm != NULL)  FreeCM(cm);
  if (buf != NULL) free(buf);
  return status;
}

/* EPN, Thu Aug 23 15:43:13 2007
 * is_integer() savagely ripped verbatim out
 * of Easel's esl_getopts.c, where it was private.
 */

/* Function: is_integer()
 * 
 * Returns TRUE if <s> points to something that atoi() will parse
 * completely and convert to an integer.
 */
static int
is_integer(char *s)
{
  int hex = 0;

  if (s == NULL) return 0;
  while (isspace((int) (*s))) s++;      /* skip whitespace */
  if (*s == '-' || *s == '+') s++;      /* skip leading sign */
				        /* skip leading conversion signals */
  if ((strncmp(s, "0x", 2) == 0 && (int) strlen(s) > 2) ||
      (strncmp(s, "0X", 2) == 0 && (int) strlen(s) > 2))
    {
      s += 2;
      hex = 1;
    }
  else if (*s == '0' && (int) strlen(s) > 1)
    s++;
				/* examine remainder for garbage chars */
  if (!hex)
    while (*s != '\0')
      {
	if (!isdigit((int) (*s))) return 0;
	s++;
      }
  else
    while (*s != '\0')
      {
	if (!isxdigit((int) (*s))) return 0;
	s++;
      }
  return 1;
}

/*
 * Function: print_searchinfo
 * Date:     EPN, Thu May 17 14:47:36 2007
 * Purpose:  Print info about search (cutoffs, algorithm, etc.) to file or stdout 
 */
int print_searchinfo(const ESL_GETOPTS *go, struct cfg_s *cfg, FILE *fp, CM_t *cm, long N, char *errbuf)
{
  int p, n;
  float surv_fract;
  int cutoff_type;
  float sc_cutoff;
  float e_cutoff;
  int using_filters;
  int cm_mode;
  int cp9_mode;
  int stype;
  int search_opts;
  ScanMatrix_t *smx;
  HybridScanInfo_t *hsi;

  /* contract check */
  if(cm->si == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "set_searchinfo(), cm->si is NULL, shouldn't happen.\n");
  SearchInfo_t *si = cm->si;

  /* Could use ESL_GETOPTS here, but using the CM flags assures we're reporting
   * on how the CM is actually config'ed, not how we want it to be
   */
    
  using_filters = (si->nrounds > 0) ? TRUE : FALSE;

  fprintf(fp,"------------------------------\n");
  for(n = 0; n <= si->nrounds; n++) {
    if(!using_filters)                        fprintf(fp, "No filtering.\n");
    else if(n < si->nrounds && using_filters) fprintf(fp, "Filter round %d of %d:\n", (n+1), si->nrounds); 
    else                                      fprintf(fp, "Final round:\n");
    
    stype       = si->stype[n];
    search_opts = si->search_opts[n];
    cutoff_type = si->cutoff_type[n];
    sc_cutoff   = si->sc_cutoff[n];
    e_cutoff    = si->e_cutoff[n];
    smx         = si->smx[n];
    hsi         = si->hsi[n];
    
    /* Determine configuration of CM and CP9 based on cm->flags & cm->search_opts */
    CM2Gumbel_mode(cm, search_opts, &cm_mode, &cp9_mode); 
    
    if(stype == SEARCH_WITH_CM) { /* using the full CM this round */
      if(cutoff_type == E_CUTOFF) { /* we use the stats in cm->stats->gumAA[cm_mode], 
				     * alternatively could store stats in search_info_t si, but I think that's 
				     * more trouble than it's worse */
	if(!(cm->flags & CMH_GUMBEL_STATS)) ESL_FAIL(eslEINCOMPAT, errbuf, "trying to use E-value for CM cutoff, but CM has no Gumbel stats.");
	fprintf(fp, "CM  cutoff (E value): %.2f\n", e_cutoff);
	for(p = 0; p < cm->stats->np; p++)
	  fprintf(fp, "   GC %2d-%3d bit sc:  %.2f mu: %.5f lambda: %.5f\n", cm->stats->ps[p], cm->stats->pe[p], 
		  (cm->stats->gumAA[cm_mode][p]->mu - (log(e_cutoff) / cm->stats->gumAA[cm_mode][p]->lambda)), 
		  cm->stats->gumAA[cm_mode][p]->mu, cm->stats->gumAA[cm_mode][p]->lambda);
      }
      else if (cutoff_type == SCORE_CUTOFF) fprintf(fp, "CM cutoff (bit sc):   %.2f\n", sc_cutoff);
      fprintf (fp, "CM search algorithm:  ");
      if(search_opts & CM_SEARCH_INSIDE) fprintf(fp, "Inside\n");
      else                               fprintf(fp, "CYK\n");
      fprintf (fp, "CM configuration:     ");
      if(cm->flags & CMH_LOCAL_BEGIN) fprintf(fp, "Local\n");
      else                            fprintf(fp, "Glocal\n");
      if(smx->dmin == NULL && smx->dmax == NULL) fprintf(fp, "QDB OFF\n");
      else                                       fprintf(fp, "QDB ON, beta:        %6g\n", smx->beta);
    }
    else if (stype == SEARCH_WITH_HMM) { /* using the HMM this round */
      if(cutoff_type == E_CUTOFF) {
	if(!(cm->flags & CMH_GUMBEL_STATS)) ESL_FAIL(eslEINCOMPAT, errbuf, "trying to use E-value for HMM cutoff, but CM has no Gumbel stats.");
	fprintf(fp, "CP9 cutoff (E value): %.2f\n", e_cutoff);
	if(n < si->nrounds) { /* we're filtering */
	  /* HMM filtering sends j-W..i+W to be re-searched with CM for HMM hits i..j */
	  /* Predict survival fraction from filter based on E-value, assume average hit length is cfg->avglen[0] (from QD band calc) */
	  surv_fract = (e_cutoff * ((2. * cm->W) - cfg->avglen[0])) / ((double) N); 
	  fprintf(fp, "   Predicted survival fraction: %.5f (1/%.3f)\n", surv_fract, (1./surv_fract));
	}
	for(p = 0; p < cm->stats->np; p++)
	fprintf(fp, "   GC %2d-%3d bit sc:  %.2f mu: %.5f lambda: %.5f\n", cm->stats->ps[p], cm->stats->pe[p], 
		(cm->stats->gumAA[cp9_mode][p]->mu - 
		 (log(e_cutoff) / cm->stats->gumAA[cp9_mode][p]->lambda)), 
		cm->stats->gumAA[cp9_mode][p]->mu, cm->stats->gumAA[cp9_mode][p]->lambda);
      }
      else if (cutoff_type == SCORE_CUTOFF)        fprintf (fp, "CP9 cutoff (bit sc):  %.2f\n", sc_cutoff);
      if      (search_opts & CM_SEARCH_HMMVITERBI) fprintf (fp, "CP9 search algorithm: Viterbi\n");
      else if (search_opts & CM_SEARCH_HMMFORWARD) fprintf (fp, "CP9 search algorithm: Forward\n");
      printf ("CP9 configuration:    ");
      if(cm->cp9->flags & CPLAN9_LOCAL_BEGIN) { 
	if(cm->cp9->flags & CPLAN9_EL) fprintf(fp, "Local (EL on)\n");
	else                           fprintf(fp, "Local (EL off)\n");
      }
      else                                    fprintf(fp, "Glocal\n");
    }
    else if(stype == SEARCH_WITH_HYBRID) {
      if(cutoff_type != SCORE_CUTOFF) ESL_FAIL(eslEINCOMPAT, errbuf, "search round %d is hybrid strategy, but with an E-value cutoff, not allowed.\n", n);
      printf("Hybrid scanner cutoff (bit sc):    %.2f\n", sc_cutoff);
      printf("TO DO: print more useful info here, like which sub CM roots are used and predicted speedup\n");
    }
    fprintf(fp,"------------------------------\n");
    if(n == si->nrounds) { 
      fprintf(fp, "DB size, nt (N):      %ld\n", N);
      fprintf(fp,"------------------------------\n");
    }
  }
  fprintf(fp, "\n");
  fflush(fp);
  return eslOK;
}

/*
 * Function: read_next_search_seq
 *
 * Date:     RJK, Wed May 29, 2002 [St. Louis]
 *           easeled: EPN, Fri Dec  8 11:40:20 2006
 *
 * Purpose:  Given a dbfp and whether or not to take the reverse complement,
 *           reads in the next sequence and prepares reverse complement.
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

#if HAVE_MPI
/* determine_seq_chunksize()
 * From RSEARCH, with one change, ideal situation is considered
 * when we put 1 chunk for each STRAND of each seq on each proc.
 *
 * Set the chunk size as follows:
 * 1.  Ideally take smallest multiple of cm->W that gives result greater than
 *     (seqlen + (cm->W * (num_procs-2))) / (num_procs-1)
 *     This should put one chunk for EACH STRAND on each processor.
 * 2.  If this is less than MPI_MIN_CHUNK_W_MULTIPLIER * cm->W, use that value.
 * 3.  If this is greater than MPI_MAX_CHUNK_SIZE, use that.
 */
static int
determine_seq_chunksize(struct cfg_s *cfg, int L, int W)
{
  int chunksize;
  chunksize = ((L + (W * (cfg->nproc-2))) / (cfg->nproc)) + 1;
  chunksize = ((chunksize / W) + 1) * W;
  chunksize = ESL_MAX(chunksize, W * MPI_MIN_CHUNK_W_MULTIPLIER); 
  chunksize = ESL_MIN(chunksize, MPI_MAX_CHUNK_SIZE);
  return chunksize;
}
#endif

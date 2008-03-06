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
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_mpi.h"
#include "esl_msa.h"
#include "esl_sqio.h"
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
  /* name           type      default  env  range     toggles      reqs       incomp  help  docgroup*/
  /* basic options */
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",   1 },
  { "-o",        eslARG_OUTFILE,NULL,  NULL, NULL,      NULL,      NULL,        NULL, "direct output to file <f>, not stdout", 1 },
  { "-g",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "configure CM/HMM for glocal alignment [default: local]", 1 },
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "force; allow search with non-calibrated cmfile",   1 },
  { "--informat",eslARG_STRING, NULL,  NULL, NULL,      NULL,      NULL,        NULL, "specify the input file is in format <x>, not FASTA", 1 },
  { "--toponly", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "only search the top strand", 1 },
  { "--bottomonly", eslARG_NONE,FALSE, NULL, NULL,      NULL,      NULL,        NULL, "only search the bottom strand", 1 },
  { "--null2",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, HMMONLYOPTS, "turn on the post hoc second null model", 1 },
  { "--forecast",eslARG_INT,    NULL,  NULL, NULL,      NULL,      NULL,        NULL, "don't do search, forecast running time with <n> processors", 1 },
  /* 4 --p* options below are hopefully temporary b/c if we have E-values for the CM using a certain cm->pbegin, cm->pend,
   * changing those values in cmsearch invalidates the E-values, so we should pick hard-coded values for cm->pbegin cm->pend */
  { "--pebegin", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "-g,--pbegin","set all local begins as equiprobable", 1 },
  { "--pfend",   eslARG_REAL,   NULL,  NULL, "0<x<1",   NULL,      NULL, "-g,--pend",  "set all local end probs to <x>", 1 },
  { "--pbegin",  eslARG_REAL,  "0.05",NULL,  "0<x<1",   NULL,      NULL,        "-g", "set aggregate local begin prob to <x>", 1 },
  { "--pend",    eslARG_REAL,  "0.05",NULL,  "0<x<1",   NULL,      NULL,        "-g", "set aggregate local end prob to <x>", 1 },
  /* options for algorithm for final round of search */
  { "--inside",  eslARG_NONE,  "default",NULL,NULL,    NULL,       NULL,    STRATOPTS1, "use scanning CM Inside algorithm", 2 },
  { "--cyk",     eslARG_NONE,  FALSE, NULL, NULL, "--inside",      NULL,    STRATOPTS1, "use scanning CM CYK algorithm", 2 },
  { "--viterbi", eslARG_NONE,  FALSE, NULL, NULL, "--fil-hmm,--fil-qdb",     NULL,    STRATOPTS2, "use scanning HMM Viterbi algorithm", 2 },
  { "--forward", eslARG_NONE,  FALSE, NULL, NULL, "--fil-hmm,--fil-qdb",     NULL,    STRATOPTS2, "use scanning HMM Forward algorithm", 2 },
  /* CM cutoff options */
  { "-E",        eslARG_REAL,   "0.1", NULL, "x>0.",    NULL,      NULL,    CUTOPTS1, "use cutoff E-value of <x> for final round of search", 3 },
  { "-T",        eslARG_REAL,   "0.0", NULL, NULL,      NULL,      NULL,    CUTOPTS1, "use cutoff bit score of <x> for final round of search", 3 },
  { "--ga",      eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    CUTOPTS2, "use CM Rfam GA gathering threshold as cutoff bit score", 3 },
  { "--tc",      eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    CUTOPTS2, "use CM Rfam TC trusted cutoff as cutoff bit score", 3 },
  { "--nc",      eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    CUTOPTS2, "use CM Rfam NC noise cutoff as cutoff bit score", 3 },
  /* banded options (for final round of searching) */
  { "--no-qdb",  eslARG_NONE,  FALSE, NULL, NULL,       NULL,      NULL,  HMMONLYOPTS, "do not use QDBs in final round of searching", 4 },
  { "--beta",    eslARG_REAL,  "1e-15",NULL, "0<x<1",   NULL,      NULL,  HMMONLYOPTS, "set tail loss prob for QDB calculation to <x>", 4 },
  { "--hbanded", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,  HMMONLYOPTS, "calculate and use HMM bands in final round of CM search", 4 },
  { "--tau",     eslARG_REAL,   "1e-7",NULL, "0<x<1",   NULL,"--hbanded", HMMONLYOPTS, "set tail loss prob for --hbanded to <x>", 4 },
  { "--aln2bands",eslARG_NONE, FALSE, NULL, NULL,      NULL, "--hbanded", HMMONLYOPTS, "w/--hbanded derive HMM bands w/o scanning Forward/Backward", 4 },
  /* filtering options, by default do HMM, then CYK filter */
  { "--fil-qdb",   eslARG_NONE, "default", NULL, NULL,  NULL,      NULL,"--fil-no-qdb", "filter with CM QDB (banded) CYK algorithm", 5 },
  { "--fil-beta",  eslARG_REAL, NULL,      NULL, "x>0", NULL,      NULL,"--fil-no-qdb", "set tail loss prob for filter QDB and W calculation to <x>", 4 },
  { "--fil-no-qdb",eslARG_NONE, FALSE,     NULL, NULL,  "--fil-qdb",NULL,         NULL, "do not filter with CM banded CYK", 5 },
  { "--fil-hmm",   eslARG_NONE, "default", NULL, NULL,  NULL,      NULL,"--fil-no-hmm", "filter with HMM forward algorithm", 5 },
  { "--fil-no-hmm",eslARG_NONE, FALSE,     NULL, NULL,  "--fil-hmm",NULL,         NULL, "do not filter with HMM forward algorithm", 5 },
  /* filter cutoff options */
  { "--fil-S-qdb",eslARG_REAL,  "0.02",NULL, "0<x<1.",  NULL,      NULL, "--fil-T-qdb", "set QDB CM filter cutoff to achieve survival fraction <x>", 6 },
  { "--fil-S-hmm",eslARG_REAL,  "0.02",NULL, "0<x<1",   NULL,      NULL, "--fil-T-hmm", "set HMM filter cutoff to achieve survival fraction <x>", 6 },
  { "--fil-T-qdb",eslARG_REAL,  "0.0", NULL, NULL,      NULL,      NULL, "--fil-S-qdb", "use cutoff bit score of <x> for QDB CM filter", 6 },
  { "--fil-T-hmm",eslARG_REAL,  "3.0", NULL, NULL,      NULL,      NULL, "--fil-S-hmm", "use cutoff bit score of <x> for HMM filter", 6 },
  ///  { "--fil-E-qdb",eslARG_REAL,  "0.02",NULL, "x>0.",    NULL,      NULL, "--fil-T-qdb", "use E-value of <x> QDB CM filter", 6 },
  ///  { "--fil-E-hmm",eslARG_REAL,  "0.02",NULL, "x>0.",    NULL,      NULL, "--fil-T-hmm", "use E-value cut of <x> for HMM filter", 6 }, 
  { "--fil-Smax-hmm",eslARG_REAL,NULL, NULL, "0<x<1",    NULL,      NULL,"--fil-T-hmm,--fil-S-hmm", "set maximum HMM survival fraction as <x>", 6 },
  /* alignment options */
  { "-p",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,"--noalign", "append posterior probabilities to hit alignments", 7 },
  { "--noalign", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,       NULL, "find start/stop/score only; don't do alignments", 7 },
  { "--optacc",  eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,"--noalign", "align hits with the Holmes/Durbin optimal accuracy algorithm", 7 },
  { "--addx",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,"--noalign", "add line to output alnments marking non-compensatory bps with 'x'", 7 },
  /* verbose output files */
  { "--glbf",    eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "save hits in simple 'glbf' format to file <f>", 8 },
  { "--gcfile",  eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "save GC content stats of target sequence file to <f>", 8 },
  /* Setting output alphabet */
  { "--rna",     eslARG_NONE,"default",NULL, NULL,  ALPHOPTS,      NULL,        NULL, "output alignment as RNA sequence data", 9 },
  { "--dna",     eslARG_NONE,   FALSE, NULL, NULL,  ALPHOPTS,      NULL,        NULL, "output alignment as DNA (not RNA) sequence data", 9 },
  /* Other options */
  { "--stall",   eslARG_NONE,  FALSE, NULL, NULL,       NULL,      NULL,        NULL, "arrest after start: for debugging MPI under gdb", 10 },  
  { "--mxsize",  eslARG_REAL, "256.0", NULL, "x>0.",    NULL,      NULL,        NULL, "set maximum allowable HMM banded DP matrix size to <x> Mb", 10 },
#ifdef HAVE_MPI
  { "--mpi",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "run as an MPI parallel program", 10 },  
#endif
  /* Development options */
#ifdef HAVE_DEVOPTS
  { "--rtrans",  eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--viterbi,--forward", "replace CM transition scores from <cmfile> with RSEARCH scores", 11 },
  { "--sums",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--hbanded",       NULL, "use posterior sums during HMM band calculation (widens bands)", 11 },
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
  FILE         *ofp;            /* output file (default is stdout) */
  int           fmt;		/* format code for seqfile */
  ESL_ALPHABET *abc;		/* digital alphabet for input */
  long          dbsize;         /* database size in nucleotides (doubled if doing rev comp) */
  int           ncm;            /* number CM we're at in file */
  int           do_rc;          /* should we search reverse complement? (for convenience */
  int           init_rci;       /* initial strand to search 0 for top, 1 for bottom (only 1 if --bottomonly enabled) */
  float         avg_hit_len;    /* average CM hit length, calc'ed using QDB calculation algorithm */

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
static void  mpi_master    (const ESL_GETOPTS *go, struct cfg_s *cfg);
static void  mpi_worker    (const ESL_GETOPTS *go, struct cfg_s *cfg);
#endif
static int initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int read_next_search_seq(const ESL_ALPHABET *abc, ESL_SQFILE *seqfp, int do_revcomp, dbseq_t **ret_dbseq);
static int print_run_info(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf);
extern int get_command(const ESL_GETOPTS *go, char *errbuf, char **ret_command);
static int set_searchinfo_for_calibrated_cm     (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int set_searchinfo_for_uncalibrated_cm   (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int print_searchinfo_for_calibrated_cm   (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, FILE *fp, CM_t *cm, float *cm_surv_fractA, int *cm_nhitsA, double in_asec, double in_total_psec, double *ret_total_psec);
static int print_searchinfo_for_uncalibrated_cm (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, FILE *fp, CM_t *cm, float *cm_surv_fractA, int *cm_nhitsA, double in_asec);
static int estimate_search_time_for_round(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, int stype, int search_opts, ScanMatrix_t *smx, ESL_RANDOMNESS *r, double *ret_sec_per_res);


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
      puts("\noptions for returning alignments of search hits:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80);
      puts("\nverbose output files:");
      esl_opt_DisplayHelp(stdout, go, 8, 2, 80);
      puts("\noptions for selecting output alphabet:");
      esl_opt_DisplayHelp(stdout, go, 9, 2, 80);
      puts("\nother options:");
      esl_opt_DisplayHelp(stdout, go, 10, 2, 80);
#ifdef HAVE_DEVOPTS
      puts("\nnon-standard options used only in development:");
      esl_opt_DisplayHelp(stdout, go, 11, 2, 80);
#endif 
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
  cfg.ofp        = NULL;
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
  cfg.dbsize     = 0;                      /* db size  */
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
      if(! esl_opt_IsDefault(go, "--forecast")) cm_Fail("--forecast is incompatible with --mpi.");
      cfg.do_mpi     = TRUE;
      MPI_Init(&argc, &argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &(cfg.my_rank));
      MPI_Comm_size(MPI_COMM_WORLD, &(cfg.nproc));

      if(cfg.nproc == 1) cm_Fail("MPI mode, but only 1 processor running... (did you execute mpirun?)");

      if (cfg.my_rank > 0)  mpi_worker(go, &cfg);
      else {
	cm_banner(stdout, argv[0], banner);
	mpi_master(go, &cfg);
      }

      esl_stopwatch_Stop(w);
      esl_stopwatch_MPIReduce(w, 0, MPI_COMM_WORLD);
      MPI_Finalize();
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
    if (! esl_opt_IsDefault(go, "-o")) { 
      printf("# Search results saved in file %s.\n", esl_opt_GetString(go, "-o"));
      fclose(cfg.ofp); 
    }
    if (cfg.cmfp      != NULL) CMFileClose(cfg.cmfp);
    if (cfg.sqfp      != NULL) esl_sqfile_Close(cfg.sqfp);
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

  /* GetDBInfo() reads all sequences, rewinds seq file and returns db size */
  GetDBInfo(NULL, cfg->sqfp, &(cfg->dbsize), NULL);  
  if ((! esl_opt_GetBoolean(go, "--toponly")) && (! esl_opt_GetBoolean(go, "--bottomonly"))) cfg->dbsize *= 2;

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
  char           errbuf[cmERRBUFSIZE];
  CM_t          *cm = NULL;
  CMConsensus_t *cons = NULL;     /* precalculated consensus info for display purposes */
  int            using_e_cutoff;
  int            rci;
  dbseq_t       *dbseq = NULL;
  int            do_top;
  float         *cm_surv_fractA  = NULL; /* 0..n..cm->si->nrounds fraction of db surviving round n for current CM */
  float         *seq_surv_fractA = NULL; /* 0..n..cm->si->nrounds fraction of db surviving round n for current seq */
  int           *cm_nhitsA = NULL;      /* 0..n..cm->si->nrounds number of hits reported for round n for current CM */
  int           *seq_nhitsA = NULL;     /* 0..n..cm->si->nrounds number of hits reported for round n for current seq */
  int            n;
  double         cm_psec;               /* predicted number of seconds for current CM versus full DB */
  double         total_psec = 0.;       /* predicted number of seconds for all CMs versus full DB */
  char           time_buf[128];	        /* for printing predicted time if --forecast only */
  ESL_STOPWATCH *w  = esl_stopwatch_Create();
  if(w == NULL) cm_Fail("serial_master(): memory error, stopwatch not created.\n");

  if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);
  /*if ((status = init_shared_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);*/
  if ((status = print_run_info (go, cfg, errbuf))  != eslOK) cm_Fail(errbuf);
  do_top = (cfg->init_rci == 0) ? TRUE : FALSE; 

  while (CMFileRead(cfg->cmfp, &(cfg->abc), &cm))
    {
      if (cm == NULL) cm_Fail("Failed to read CM from %s -- file corrupt?\n", cfg->cmfile);
      if((! (cm->flags & CMH_EXPTAIL_STATS)) && (! esl_opt_IsDefault(go, "--forecast"))) cm_Fail("--forecast only works with calibrated CM files. Run cmcalibrate (please)."); 
      cfg->ncm++;
      if(cfg->ncm == 1) { /* check if we have exp stats */
	if((! (cm->flags & CMH_EXPTAIL_STATS)) && (! esl_opt_GetBoolean(go, "-F"))) {
	  cm_Fail("%s has not been calibrated with cmcalibrate.\n       Once calibrated, E-values will be available and HMM filter thresholds\n       will be appropriately set to accelerate the search.\n       To override this error without calibrating, use -F.", cfg->cmfile);
	}
      }

      /* initialize the flags/options/params and configuration of the CM */
      if((  status = initialize_cm(go, cfg, errbuf, cm))                    != eslOK) cm_Fail(errbuf);

      if((  status = cm_GetAvgHitLen(cm, errbuf, &(cfg->avg_hit_len)))      != eslOK) cm_Fail(errbuf);
      if((  status = CreateCMConsensus(cm, cfg->abc_out, 3.0, 1.0, &cons))  != eslOK) cm_Fail(errbuf);
      if(cm->flags & CMH_EXPTAIL_STATS) { 
	if((status = UpdateExpsForDBSize(cm, errbuf, cfg->dbsize))          != eslOK) cm_Fail(errbuf);
	if((status = set_searchinfo_for_calibrated_cm(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);
      }
      else { if((status = set_searchinfo_for_uncalibrated_cm(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf); }

      ESL_ALLOC(cm_surv_fractA, sizeof(float) * (cm->si->nrounds+1));
      ESL_ALLOC(cm_nhitsA,      sizeof(int) * (cm->si->nrounds+1));
      esl_vec_FSet(cm_surv_fractA, (cm->si->nrounds+1), 0.);
      esl_vec_ISet(cm_nhitsA,     (cm->si->nrounds+1), 0);
      if(cm->flags & CMH_EXPTAIL_STATS) { if((status = print_searchinfo_for_calibrated_cm   (go, cfg, errbuf, stdout, cm, NULL, NULL, 0., 0., &cm_psec)) != eslOK) cm_Fail(errbuf); }
      else                              { if((status = print_searchinfo_for_uncalibrated_cm(go, cfg, errbuf, stdout, cm, NULL, NULL, 0.)) != eslOK) cm_Fail(errbuf); }

      if(! esl_opt_IsDefault(go, "--forecast")) { /* special mode, we don't do the search, just print the predicting timings */
	total_psec += cm_psec;
	free(cm_surv_fractA);
	free(cm_nhitsA);
	continue;
      }

      fprintf(cfg->ofp, "Model: %s\n", cm->name);

      using_e_cutoff = (cm->si->cutoff_type[cm->si->nrounds] == E_CUTOFF) ? TRUE : FALSE;
	 
      esl_stopwatch_Start(w);
      while ((status = read_next_search_seq(cfg->abc, cfg->sqfp, cfg->do_rc, &dbseq)) == eslOK)
	{
	  for(rci = cfg->init_rci; rci <= cfg->do_rc; rci++) {
	    /*printf("SEARCHING >%s %d\n", dbseq->sq[reversed]->name, reversed);*/
	    if ((status = ProcessSearchWorkunit(cm, errbuf, dbseq->sq[rci]->dsq, dbseq->sq[rci]->n, &dbseq->results[rci], esl_opt_GetReal(go, "--mxsize"), cfg->my_rank, &seq_surv_fractA, &seq_nhitsA)) != eslOK) cm_Fail(errbuf);
	    for(n = 0; n <= cm->si->nrounds; n++) { 
	      cm_surv_fractA[n] += (dbseq->sq[rci]->n * seq_surv_fractA[n]);
	      cm_nhitsA[n]      += seq_nhitsA[n];
	    }
	    free(seq_surv_fractA);
	    free(seq_nhitsA);
	    RemoveOverlappingHits(dbseq->results[rci], 1, dbseq->sq[rci]->n);
	    if(using_e_cutoff) RemoveHitsOverECutoff(cm, cm->si, dbseq->results[rci], dbseq->sq[rci]); 
	  }
	  PrintResults (cm, cfg->ofp, cm->si, cfg->abc_out, cons, dbseq, do_top, cfg->do_rc, esl_opt_GetBoolean(go, "--addx"));
	  for(rci = 0; rci <= cfg->do_rc; rci++) { /* we can free results for top strand even if cfg->init_rci is 1, due to --bottomonly */
	    FreeResults(dbseq->results[rci]);
	    esl_sq_Destroy(dbseq->sq[rci]);
	  }
	  free(dbseq);
	}
      esl_stopwatch_Stop(w);
      if (status != eslEOF) cm_Fail("Parse failed, line %d, file %s:\n%s", 
				    cfg->sqfp->linenumber, cfg->sqfp->filename, cfg->sqfp->errbuf);
      /* convert cm_surv_fractA[] values from residue counts into fractions */
      for(n = 0; n <= cm->si->nrounds; n++) cm_surv_fractA[n] /= (double) (cfg->dbsize);
      if(cm->flags & CMH_EXPTAIL_STATS) { if((status = print_searchinfo_for_calibrated_cm   (go, cfg, errbuf, stdout, cm, cm_surv_fractA, cm_nhitsA, w->elapsed, cm_psec, NULL)) != eslOK) cm_Fail(errbuf); }
      else                              { if((status = print_searchinfo_for_uncalibrated_cm(go, cfg, errbuf, stdout, cm, cm_surv_fractA, cm_nhitsA, w->elapsed)) != eslOK) cm_Fail(errbuf); }
      fprintf(cfg->ofp, "//\n");
      FreeCM(cm);
      FreeCMConsensus(cons);
      free(cm_surv_fractA);
      free(cm_nhitsA);
      esl_sqio_Rewind(cfg->sqfp); /* we may be searching this file again with another CM */
    }

  if(cfg->ncm > 1 && (! esl_opt_IsDefault(go, "--forecast"))) { 
    fprintf(stdout, "#\n");
    fprintf(stdout, "# %20s\n", "predicted total time");
    fprintf(stdout, "# %20s\n", "--------------------");
    FormatTimeString(time_buf, total_psec, FALSE);
    fprintf(stdout, "  %20s\n", time_buf);
  }
      
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
  dbseq_t  *dbseq = NULL;  /* a database sequence */
  double    cm_psec;                /* predicted number of seconds for current CM versus full DB (ON MASTER PROC BUT WE ASSUME
				     * OTHER PROCS ARE THE SAME SPEED!) */
  float    *cm_surv_fractA = NULL;  /* 0..n..cm->si->nrounds fraction of db that survived round n for current CM */
  int      *cm_nhitsA = NULL;       /* 0..n..cm->si->nrounds number of hits reported for round n for current CM */
  ESL_STOPWATCH *w  = esl_stopwatch_Create();
  if(w == NULL) cm_Fail("mpi_master(): memory error, stopwatch not created.\n");  

  char     errbuf[cmERRBUFSIZE];
  MPI_Status mpistatus; 
  int      n;

  int need_seq = TRUE;
  int chunksize;
  search_results_t *worker_results;
  float            *worker_surv_fractA = NULL; /* 0..n..cm->si->nrounds fraction of db surviving round n for current seq */
  int              *worker_nhitsA      = NULL; /* 0..n..cm->si->nrounds num hits surviving round n for current seq */

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
  if (xstatus == eslOK) { if ((status = print_run_info(go, cfg, errbuf))  != eslOK) xstatus = status; }

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

  while (xstatus == eslOK && CMFileRead(cfg->cmfp, &(cfg->abc), &cm))
    {
      cfg->ncm++;  
      ESL_DPRINTF1(("MPI master read CM number %d\n", cfg->ncm));

      if(cfg->ncm == 1) { /* check if we have exp stats */
	if((! (cm->flags & CMH_EXPTAIL_STATS)) && (! esl_opt_GetBoolean(go, "-F"))) { 
	  cm_Fail("%s has not been calibrated with cmcalibrate.\n       Once calibrated, E-values will be available and HMM filter thresholds\n       will be appropriately set to accelerate the search.\n       To override this warning without calibrating, use -F.", cfg->cmfile);
	}
      }

      if((status = cm_master_MPIBcast(cm, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI broadcast CM failed.");
      
      /* initialize the flags/options/params of the CM */
      if((status   = initialize_cm(go, cfg, errbuf, cm))                    != eslOK) cm_Fail(errbuf);
      if((status   = cm_GetAvgHitLen(cm, errbuf, &(cfg->avg_hit_len)))      != eslOK) cm_Fail(errbuf);
      if((status   = CreateCMConsensus(cm, cfg->abc_out, 3.0, 1.0, &cons))  != eslOK) cm_Fail(errbuf);
      if(cm->flags & CMH_EXPTAIL_STATS) {
	if((status = UpdateExpsForDBSize(cm, errbuf, cfg->dbsize))          != eslOK) cm_Fail(errbuf);
	if((status = set_searchinfo_for_calibrated_cm(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf);
      }
      else { if((status = set_searchinfo_for_uncalibrated_cm(go, cfg, errbuf, cm)) != eslOK) cm_Fail(errbuf); }

      using_e_cutoff = (cm->si->cutoff_type[cm->si->nrounds] == E_CUTOFF) ? TRUE : FALSE;

      ESL_ALLOC(cm_surv_fractA, sizeof(float) * (cm->si->nrounds+1));
      ESL_ALLOC(cm_nhitsA,      sizeof(int) * (cm->si->nrounds+1));
      esl_vec_FSet(cm_surv_fractA, (cm->si->nrounds+1), 0.);
      esl_vec_ISet(cm_nhitsA,      (cm->si->nrounds+1), 0);
      if(cm->flags & CMH_EXPTAIL_STATS) { if((status = print_searchinfo_for_calibrated_cm  (go, cfg, errbuf, stdout, cm, NULL, NULL, 0., 0., &cm_psec)) != eslOK) cm_Fail(errbuf); }
      else                              { if((status = print_searchinfo_for_uncalibrated_cm(go, cfg, errbuf, stdout, cm, NULL, NULL, 0.)) != eslOK) cm_Fail(errbuf); }

      fprintf(cfg->ofp, "Model: %s\n", cm->name);

      /* reset vars for searching with current CM */
      wi = 1;
      ndbseq = 0;
      need_seq = TRUE;
      have_work = TRUE;	/* TRUE while work remains  */
      seqpos = 1;
      in_rc = FALSE;
      esl_stopwatch_Start(w);
      while (have_work || nproc_working)
	{
	  if (need_seq) 
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
		  chunksize = DetermineSeqChunksize(cfg->nproc, dbseq->sq[0]->n, cm->W);
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

		      /* first unpack the surv_fractA and nhitsA, that give survival fractions and num hits surviving for each round */
		      ESL_ALLOC(worker_surv_fractA, sizeof(float) * (cm->si->nrounds+1));
		      ESL_ALLOC(worker_nhitsA,      sizeof(int)   * (cm->si->nrounds+1));
		      if (MPI_Unpack(buf, bn, &pos, worker_surv_fractA, (cm->si->nrounds+1), MPI_FLOAT, MPI_COMM_WORLD) != 0) cm_Fail("mpi unpack failed");
		      if (MPI_Unpack(buf, bn, &pos, worker_nhitsA,      (cm->si->nrounds+1), MPI_INT,   MPI_COMM_WORLD) != 0) cm_Fail("mpi unpack failed");

		      if ((status = cm_search_results_MPIUnpack(buf, bn, &pos, MPI_COMM_WORLD, &worker_results)) != eslOK)    cm_Fail("search results unpack failed");
		      ESL_DPRINTF1(("MPI master has unpacked search results\n"));
		      
		      /* update cm_surv_fractA[] and cm_nhitsA[] which holds number of residues and hits that survived each round of searching/filtering */
		      for(n = 0; n <= cm->si->nrounds; n++) { 
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
			/* careful, dbseqlist[si_recv]->results[rclist[wi]] now points to the nodes in worker_results->data,
			 * don't free those (don't use FreeResults(worker_results)) */
			free(worker_results);
			worker_results = NULL;
		      }
		      dbseqlist[si_recv]->chunks_sent--;
		      if(sentlist[si_recv] && dbseqlist[si_recv]->chunks_sent == 0)
			{
			  for(rci = 0; rci <= cfg->do_rc; rci++) {
			    RemoveOverlappingHits(dbseqlist[si_recv]->results[rci], 1, dbseqlist[si_recv]->sq[rci]->n);
			    if(using_e_cutoff) RemoveHitsOverECutoff(cm, cm->si, dbseqlist[si_recv]->results[rci], dbseqlist[si_recv]->sq[rci]);
			  }					      
			  PrintResults(cm, cfg->ofp, cm->si, cfg->abc_out, cons, dbseqlist[si_recv], TRUE, cfg->do_rc, esl_opt_GetBoolean(go, "--addx"));
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
      if(cm->flags & CMH_EXPTAIL_STATS) { if((status = print_searchinfo_for_calibrated_cm   (go, cfg, errbuf, stdout, cm, cm_surv_fractA, cm_nhitsA, w->elapsed, cm_psec, NULL)) != eslOK) cm_Fail(errbuf); }
      else                              { if((status = print_searchinfo_for_uncalibrated_cm(go, cfg, errbuf, stdout, cm, cm_surv_fractA, cm_nhitsA, w->elapsed)) != eslOK) cm_Fail(errbuf); }
      fprintf(cfg->ofp, "//\n");
      free(cm_surv_fractA);
      free(cm_nhitsA);
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
  
  esl_stopwatch_Destroy(w);

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
  float         *surv_fractA = NULL; /* 0..n..cm->si->nrounds fraction of db surviving round n for current seq */
  int           *nhitsA      = NULL; /* 0..n..cm->si->nrounds number of hits surviving round n for current seq */

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
  MPI_Bcast(&(cfg->dbsize), 1, MPI_LONG, 0, MPI_COMM_WORLD);
  if (xstatus == eslOK) { wn = 4096;  if ((wbuf = malloc(wn * sizeof(char))) == NULL) xstatus = eslEMEM; }
  MPI_Reduce(&xstatus, &status, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD); /* everyone sends xstatus back to master */
  if (xstatus != eslOK) {
    if (wbuf != NULL) free(wbuf);
    return; /* shutdown; we passed the error back for the master to deal with. */
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

/* initialize_cm()
 * Setup the CM based on the command-line options/defaults;
 * only set flags and a few parameters. ConfigCM() configures
 * the CM. 
 */
static int
initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  int use_hmmonly;
  int nstarts, nexits, nd;
  
  /* set up CM parameters that are option-changeable */
  cm->beta_qdb = esl_opt_GetReal(go, "--beta");
  cm->tau      = esl_opt_GetReal(go, "--tau");  /* this will be DEFAULT_TAU unless changed at command line */

  use_hmmonly = (esl_opt_GetBoolean(go, "--viterbi") || esl_opt_GetBoolean(go, "--forward")) ? TRUE : FALSE;

  /* Update cm->config_opts and cm->align_opts based on command line options */

  /* config_opts */
  if(! esl_opt_GetBoolean(go, "-g")) { 
    cm->config_opts |= CM_CONFIG_LOCAL;
    cm->config_opts |= CM_CONFIG_HMMLOCAL;
    cm->config_opts |= CM_CONFIG_HMMEL;
  }
  /* config QDB for final round of search? yes, unless --no-qdb, otherwise we only use QDB to filter */
  if(! esl_opt_GetBoolean(go, "--no-qdb")) cm->config_opts |= CM_CONFIG_QDB;

  /* search_opts */
  if(! use_hmmonly) 
    if(  esl_opt_GetBoolean(go, "--inside"))    cm->search_opts |= CM_SEARCH_INSIDE;
  if(  esl_opt_GetBoolean(go, "--noalign"))     cm->search_opts |= CM_SEARCH_NOALIGN;
  if(  esl_opt_GetBoolean(go, "--null2"))       cm->search_opts |= CM_SEARCH_NULL2;
  if(  esl_opt_GetBoolean(go, "--no-qdb"))      cm->search_opts |= CM_SEARCH_NOQDB;
  if(  esl_opt_GetBoolean(go, "--hbanded"))     cm->search_opts |= CM_SEARCH_HBANDED;
  if(  esl_opt_GetBoolean(go, "--aln2bands"))   cm->search_opts |= CM_SEARCH_HMMALNBANDS;
  if(  esl_opt_GetBoolean(go, "--viterbi"))  { 
    cm->search_opts |= CM_SEARCH_HMMVITERBI;
    cm->search_opts |= CM_SEARCH_NOQDB;
  }
  if(  esl_opt_GetBoolean(go, "--forward"))  { 
    cm->search_opts |= CM_SEARCH_HMMFORWARD;
    cm->search_opts |= CM_SEARCH_NOQDB;
  }

  /* align_opts, by default, align with HMM bands */
  cm->align_opts |= CM_ALIGN_HBANDED;
  if(esl_opt_GetBoolean(go, "--optacc"))       cm->align_opts |= CM_ALIGN_OPTACC;
  if(esl_opt_GetBoolean(go, "-p"))             cm->align_opts |= CM_ALIGN_POST;

  /* handle special developer's options, not recommend for normal users */
#ifdef HAVE_DEVOPTS
  if(  esl_opt_GetBoolean(go, "--rtrans"))      cm->flags       |= CM_RSEARCHTRANS;
  if(  esl_opt_GetBoolean(go, "--sums"))        cm->search_opts |= CM_SEARCH_SUMS;
#endif

  /* TEMPORARY BLOCK */
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
  /* END TEMPORARY BLOCK */

  /* finally, configure the CM for search based on cm->config_opts and cm->align_opts.
   * set local mode, make cp9 HMM, calculate QD bands etc. 
   */
  ConfigCM(cm, TRUE);  /* TRUE says: calculate W */

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

/* Function: set_searchinfo_for_calibrated_cm)
 * Date:     EPN, Mon Jan 21 08:56:04 2008 (updated)
 * 
 * Purpose:  For a CM WITH exponential tail and filter thresholds statistics 
 *           from cmcalibrate: determine how many rounds of searching we will 
 *           do (all rounds but last round are filters), and set the relevant 
 *           info in the SearchInfo_t <cm->si> object, including cutoffs.
 *
 *************************************************************************
 * FIX EVERYTHING BELOW THIS LINE !!!!
 *************************************************************************
 * Filters:
 *
 * User can specify 0 or 1 round of HMM filtering with Forward, followed by 0 or 
 * 1 round of CM filtering with QDB CYK. This is determined as described below. 
 * 
 * HMM filtering:
 * ('none' below means none of --fil-no-hmm, --fil-T-hmm, --fil-S-hmm, --fil-Smax-hmm) 
 * 
 *                  exp tail &
 *                  filter thr 
 *   options        in cmfile?    filter/don't filter and cutoff determination
 *   -------       ----------    ------------------------------------------------
 * 1. none                 yes    filter (usually). automatically determine appropriate HMM 
 *                                filter cutoff for final round cutoff using cmfile's filter 
 *                                threshold info from cmcalibrate. Sometimes this info may
 *                                tell us it's inappropriate to use an HMM filter, in which
 *                                case, don't filter.
 * 
 * 2. none                  no    filter. set filter bit score cutoff to default value <x> for 
 *                                --fil-T-hmm (3.0 bits)
 *                               
 * 3. --fil-no-hmm      yes/no    don't filter.
 * 
 * 4. --fil-T-hmm <x>   yes/no    filter. set filter bit score cutoff as <x>
 * 
 * 5. --fil-S-hmm <x>      yes    filter. set E-value filter cutoff to E-value that gives
 *                                predicted survival fraction of <x>. NOTE: larger models,
 *                                with larger Ws, will have smaller E-value cutoffs for
 *                                same <x> (a somewhat unsettling behavior).
 * 
 * 6. --fil-S-hmm <x>       no    die. we can't deal with this combo, tell the user.
 * 
 * 7. --fil-Smax-hmm <x>   yes    filter, same as in case 1 above, but if the 'appropriate'
 *                                CM E-value cutoff gives a predicted survival fraction > <x>,
 *                                reset it as the E-value cutoff that gives <x> (even if 
 *                                cmfile told us to turn filtering off).
 * 
 * 8. --fil-Smax-hmm <x>    no    die. we can't deal with this combo, tell the user.
 *
 * 
 * CM filtering: 
 * ('none' below means none of --fil-no-qdb, --fil-T-qdb, --fil-S-qdb) 
 *                               
 *                                final 
 *                  exp tail &    round 
 *                  filter thr    cutoff
 *   options        in cmfile?    type    filter/don't filter and cutoff determination
 *   -------       -----------    ------  -----------------------------------------
 * 1. none                 yes    E val   filter. set E-value filter cutoff to E-value that gives
 *                                        predicted survival fraction of <x>. NOTE: larger models,
 *                                        with larger Ws, will have smaller E-value cutoffs for
 *                                        same <x> (a somewhat unsettling behavior).
 * 
 * 2. none                  no    filter. set filter bit score cutoff to default value <x> for 
 *                                --fil-T-qdb (0.0 bits)
 *                               
 * 3. --fil-no-qdb      yes/no    don't filter. 
 * 
 * 4. --fil-T-qdb <x>   yes/no    filter. set filter bit score cutoff as <x>
 * 
 * 5. --fil-S-qdb <x>      yes    filter. set E-value filter cutoff to E-value that gives
 *                                predicted survival fraction of <x>. NOTE: larger models,
 *                                with larger Ws, will have smaller E-value cutoffs for
 *                                same <x> (a somewhat unsettling behavior).
 * 
 * 6. --fil-S-qdb <x>       no    die. we can't deal with this combo, tell the user.
 *
 * ('none' below means none of --fil-no-qdb, --fil-T-qdb, --fil-S-qdb) 
 * 
 * beta, the tail loss probability for the QDB calculation is determined as follows:
 * If --fil-beta <x> is not enabled, cm->beta_qdb is set as cm->beta_W as read in the cmfile
 *               and initially set by cmbuild.
 * If --fil-beta <x> is enabled, <x> is used as cm->beta_qdb
 *
 * Final round related options (after all filtering is complete):
 *
 * --inside:     search with CM inside (TRUE by default)
 * --cyk:        search with CM CYK 
 * --forward:    search with HMM Forward
 * --viterbi:    search with HMM Viterbi
 * -T:           CM/HMM bit score threshold
 * -E:           CM/HMM E-value threshold (requires exp tail info in CM file)
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
  int   final_ctype;          /* final round cutoff type, SCORE_CUTOFF or E_CUTOFF */
  int   fqdb_ctype;           /* QDB filter round cutoff type, SCORE_CUTOFF or E_CUTOFF */
  int   fhmm_ctype;           /* HMM filter round cutoff type, SCORE_CUTOFF or E_CUTOFF */
  float final_S = -1;         /* predicted survival fraction from final round */
  float fqdb_S = -1;          /* predicted survival fraction from qdb filter round */
  float fhmm_S = -1;          /* predicted survival fraction from HMM filter round */
  float fhmm_ncalcs_per_res;  /* number of millions of filter HMM DP calcs predicted per residue */
  float fqdb_ncalcs_per_res;  /* number of millions of filter QDB DP calcs predicted per residue */
  float final_ncalcs_per_res; /* number of millions of final stage DP calcs predicted per residue */
  float all_filters_ncalcs_per_res;/* number of millions of DP calcs predicted for all filter rounds */
  float fhmm_Smin = 0.;       /* minimally useful survival fraction for hmm filter round */
  float fqdb_Smin = 0.;       /* minimally useful survival fraction for qdb filter round */
  float xfil = 0.01;          /* used to set *_Smin values, minimal fraction of filter dp calcs to do in the final round */
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
  int cut_point;              /* HMM forward E-value cut point from filter threshold stats */

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
  if((! use_hmmonly) && (stype != SEARCH_WITH_CM)) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_searchinfo_for_calibrated_cm(), search_opts for final round of search does not have HMMVITERBI or HMMFORWARD flags raised, but is not of type SEARCH_WITH_CM.");
  /* determine configuration of CM and CP9 HMM based on cm->flags & cm->search_opts */
  CM2ExpMode(cm, search_opts, &cm_mode, &hmm_mode); 

  final_S = final_E = final_sc = -1.;
  if((status = cm_CountSearchDPCalcs(cm, errbuf, 10*cm->smx->W, cm->smx->dmin, cm->smx->dmax, cm->smx->W, TRUE,  NULL, &(final_ncalcs_per_res))) != eslOK) return status;
  /* set up final round cutoff, either 0 or 1 of 5 options is enabled. 
   * esl_opt_IsDefault() returns FALSE even if option is enabled with default value.
   */
  if(esl_opt_IsDefault(go, "-E") && 
     esl_opt_IsDefault(go, "-T") && 
     esl_opt_IsDefault(go, "--ga") && 
     esl_opt_IsDefault(go, "--tc") && 
     esl_opt_IsDefault(go, "--nc")) { 
    /* No relevant options enabled, cutoff is default E value cutoff */
    final_ctype = E_CUTOFF;
    final_E     = esl_opt_GetReal(go, "-E");
    if((status  = E2MinScore(cm, errbuf, (use_hmmonly ? hmm_mode : cm_mode), final_E, &final_sc)) != eslOK) return status;
  }
  else if(! esl_opt_IsDefault(go, "-E")) { /* -E enabled, use that */
    final_ctype = E_CUTOFF;
    final_E     = esl_opt_GetReal(go, "-E");
    if((status = E2MinScore(cm, errbuf, (use_hmmonly ? hmm_mode : cm_mode), final_E, &final_sc)) != eslOK) return status;
  }
  else if(! esl_opt_IsDefault(go, "-T")) { /* -T enabled, use that */
    final_ctype = SCORE_CUTOFF;
    final_sc    = esl_opt_GetReal(go, "-T");
    final_E     = -1.; /* invalid, we'll never use it */  
  }
  else if(! esl_opt_IsDefault(go, "--ga")) { /* --ga enabled, use that, if available, else die */
    if(use_hmmonly) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "--ga is incompatible with --viterbi and --forward.");
    if(! (cm->flags & CMH_GA)) ESL_FAIL(eslEINVAL, errbuf, "No GA gathering threshold in CM file, can't use --ga.");
    final_ctype = SCORE_CUTOFF;
    final_sc    = cm->ga;
    final_E     = -1.; /* we'll never use it */
  }
  else if(! esl_opt_IsDefault(go, "--tc")) { /* --tc enabled, use that, if available, else die */
    if(use_hmmonly) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "--tc is incompatible with --viterbi and --forward.");
    if(! (cm->flags & CMH_TC)) ESL_FAIL(eslEINVAL, errbuf, "No TC trusted cutoff in CM file, can't use --tc.");
    final_ctype = SCORE_CUTOFF;
    final_sc    = cm->tc;
    final_E     = -1.; /* we'll never use it */
  }
  else if(! esl_opt_IsDefault(go, "--nc")) { /* --nc enabled, use that, if available, else die */
    if(use_hmmonly) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "--nc is incompatible with --viterbi and --forward.");
    if(! (cm->flags & CMH_NC)) ESL_FAIL(eslEINVAL, errbuf, "No NC noise cutoff in CM file, can't use --nc.");
    final_ctype = SCORE_CUTOFF;
    final_sc    = cm->nc;
    final_E     = -1.; /* we'll never use it */
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
  /* DumpSearchInfo(cm->si); */
  if(final_S >= 0.09) do_qdb_filter = FALSE; /* we want QDB filter to let through 10X survival fraction of final round, so if final round S > 0.09, qdb filter would only filter out 10% of database, so we turn it off */
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

  if(do_hmm_filter) { /* determine thresholds for HMM forward filter */
    if(use_hmmonly) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_searchinfo_for_calibrated_cm(), --fil-hmm enabled, along with --viterbi or --forward, shouldn't happen.");
    if(esl_opt_IsDefault(go, "--fil-S-hmm") && esl_opt_IsDefault(go, "--fil-T-hmm")) { /* default: use HMM filter threshold stats, if they exist in cmfile, else use default bit score cutoff */
      /* No relevant options selected. Set HMM filter cutoff as appropriate HMM E value cutoff from cmfile */
      /* determine filter threshold mode, the mode of final stage of searching, either FTHR_CM_LC,
       * FTHR_CM_LI, FTHR_CM_GC, FTHR_CM_GI (can't be an HMM mode b/c --viterbi and --forward toggle --fil-hmm off)
       */

      /* Determine the 'minimally useful' HMM filter round survival fraction. 
       * This is the survival fraction <fhmm_Smin> that we predict would require the 
       * final round (if no more filtering occured) to perform x millions of 
       * DP calcs, where x is 1 + <xfil> * the number of predicted DP calcs predicted 
       * for the HMM filter round.
       */
      fhmm_Smin = xfil * (fhmm_ncalcs_per_res / final_ncalcs_per_res); 
      /* if HMM filter survival fraction == fhmm_Smin, total number of DP calcs for the 
       * next round == <xfil> * fhmm_calcs_per_res, so we predict the final round will 
       * require at most <xfil> fraction of the DP calcs that the HMM filter round requires. */

      if((status = CM2FthrMode(cm, errbuf, cm->search_opts, &fthr_mode)) != eslOK) return status;
      HMMFilterInfo_t *hfi_ptr = cm->stats->hfiA[fthr_mode]; /* for convenience */
      if(hfi_ptr->is_valid == FALSE) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_searchinfo_for_calibrated_cm(), cm's CMH_FILTER_STATS is raised but best filter info for fthr_mode %d is invalid.", fthr_mode);
      fhmm_ctype = E_CUTOFF;
      /* determine the appropriate filter cut point <cut_point> to use, for details
       * see comments in function: searchinfo.c:GetHMMFilterFwdECutGivenCME() for details */
      if(final_ctype == SCORE_CUTOFF) { /* final round has bit score cutoff */
	if((status = GetHMMFilterFwdECutGivenCMBitScore(hfi_ptr, errbuf, final_sc, cfg->dbsize, &cut_point, cm, cm_mode)) != eslOK) return status; 
      }
      else { /* final round has E-value cutoff */
	if((status = GetHMMFilterFwdECutGivenCME(hfi_ptr, errbuf, final_E, cfg->dbsize, &cut_point)) != eslOK) return status; 
      }
      if(cut_point != -1) { /* it's worth it to filter */
	fhmm_E = hfi_ptr->fwd_E_cut[cut_point] * ((double) cfg->dbsize / (double) hfi_ptr->dbsize); 
	fhmm_S = E2SurvFract(fhmm_E, cm->W, cfg->avg_hit_len, cfg->dbsize, TRUE);
	/* check if --fil-Smax-hmm applies */
	if(! esl_opt_IsDefault(go, "--fil-Smax-hmm")) { 
	  if(fhmm_S > esl_opt_GetReal(go, "--fil-Smax-hmm")) { /* predicted survival fraction exceeds maximum allowed, set E cutoff as E value that gives max allowed survival fraction */
	    fhmm_E = SurvFract2E(esl_opt_GetReal(go, "--fil-Smax-hmm"), cm->W, cfg->avg_hit_len, cfg->dbsize);
	  }
	}
	if((status  = E2MinScore(cm, errbuf, cm_mode, fhmm_E, &fhmm_sc)) != eslOK) return status; /* note: use cm_mode, not fthr_mode */
      }
      else { 
	do_hmm_filter = FALSE;
	/* it's not worth it to filter, our HMM filter cutoff would be so low, 
	 * letting so much of the db survive, the filter is a waste of time */
	ESL_DPRINTF1(("cut_point -1, always_better FALSE\n"));
      }
    } /* end of if(esl_opt_IsDefault(go, "--fil-S-hmm") && esl_opt_IsDefault(go, "--fil-T-hmm")) */
    else if(! esl_opt_IsDefault(go, "--fil-S-hmm")) {
      fhmm_ctype = E_CUTOFF;
      fhmm_E     = SurvFract2E(esl_opt_GetReal(go, "--fil-S-hmm"), cm->W, cfg->avg_hit_len, cfg->dbsize);
      fhmm_E     = ESL_MIN(fhmm_E, 1.0); /* we don't allow filter E cutoffs below 1. */
    }
    else if(! esl_opt_IsDefault(go, "--fil-T-hmm")) {
      fhmm_ctype = SCORE_CUTOFF;
      fhmm_sc    = esl_opt_GetReal(go, "--fil-T-hmm");
    }
    else ESL_FAIL(eslEINCONCEIVABLE, errbuf, "No HMM filter cutoff selected. This shouldn't happen.");

    if(do_hmm_filter) { 
      if(fhmm_ctype == SCORE_CUTOFF) { /* determine max E-value that corresponds to fhmm_sc bit sc cutoff  across all partitions */
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
    if(! esl_opt_IsDefault(go, "--fil-beta")) fqdb_beta_qdb = esl_opt_GetReal(go, "--fil-beta");
    else fqdb_beta_qdb = cm->beta_W; /* use beta used to calc W in CM file by default */
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

    /* Determine the 'minimally useful' QDB CYK filter round survival fraction. 
     * This is the survival fraction <fqdb_Smin> that we predict would require the 
     * *final* round to perform <xfil> * x millions of DP calcs, where x is the number of
     * predicted DP calcs predicted for all filter rounds.
     */
    if(do_hmm_filter) all_filters_ncalcs_per_res = fhmm_ncalcs_per_res + (fqdb_ncalcs_per_res * fhmm_S);
    else              all_filters_ncalcs_per_res = fqdb_ncalcs_per_res;

    fqdb_Smin = xfil * (all_filters_ncalcs_per_res / final_ncalcs_per_res); 
    /* if QDB filter survival fraction == fqdb_Smin, total number of DP calcs for the 
     * final round == <xfil> * all_filters_ncalcs_per_res, so we predict the final round will 
     * require <xfil> fraction of the DP calcs that the filter rounds require. */

    if(esl_opt_IsDefault(go, "--fil-S-qdb") && esl_opt_IsDefault(go, "--fil-T-qdb")) {
      /* No relevant options selected, set CM QDB filter as E value that gives default survival fraction. */
      /* HERE HERE HERE, do the following for the CM QDB cutoff is: 
       * check. fqdb_S = ESL_MIN(fqdb_S, 1000. * final_S);
       * check. fqdb_S gives E value cutoff of at least 1
       * check. fqdb_S <= fhmm_S 
       * time for fqdb is at least 1.1X time for predicted HMM scan
       */
      fqdb_ctype = E_CUTOFF;
      fqdb_S     = esl_opt_GetReal(go, "--fil-S-qdb");
      fqdb_S     = ESL_MAX(fqdb_S, final_S * 10.);   /* final_S must be < 0.09 (enforced above), or else do_qdb_filter was set to FALSE */
      fqdb_S     = ESL_MIN(fqdb_S, final_S * 1000.); /* final_S must be < 0.09 (enforced above), or else do_qdb_filter was set to FALSE */
      if(do_hmm_filter) fqdb_S = ESL_MIN(fqdb_S, fhmm_S); /* predicted survival fraction from QDB filter can't be more than predicted survival fraction from the HMM filter */
      fqdb_S     = ESL_MAX(fqdb_S, fqdb_Smin); /* we never want to exceed our minimally useful survival fraction, that will make the final round require 10% DP calcs of all filter rounds */
#if 0
      if(do_hmm_filter && fqdb_S > fhmm_S) { 
	  /* predicted survival fraction from QDB filter is higher than predicted survival fraction from the HMM filter,
	   * this means HMM should be a better filter, turn OFF QDB filter */
	  do_qdb_filter = FALSE; 
	}
	else { 
#endif
	  fqdb_E     = SurvFract2E(fqdb_S, fqdb_smx->W, cfg->avg_hit_len, cfg->dbsize);
	  fqdb_E     = ESL_MAX(fqdb_E, 1.0); /* we don't allow filter E cutoffs below 1. */
#if 0
	}
#endif
    }
    else if(! esl_opt_IsDefault(go, "--fil-S-qdb")) { /* survival fraction for QDB filter set on command line, use that */
      fqdb_ctype = E_CUTOFF;
      fqdb_S     = esl_opt_GetReal(go, "--fil-S-qdb");
      //fqdb_S     = ESL_MAX(fqdb_S, final_S * 10.); /* final_S must be < 0.09 (enforced above), or else do_qdb_filter was set to FALSE */
      if(do_hmm_filter) fqdb_S = ESL_MIN(fqdb_S, fhmm_S); /* predicted survival fraction from QDB filter can't be more than predicted survival fraction from the HMM filter */
      fqdb_E     = SurvFract2E(fqdb_S, fqdb_smx->W, cfg->avg_hit_len, cfg->dbsize);
      fqdb_E     = ESL_MAX(fqdb_E, 1.0); /* we don't allow filter E cutoffs below 1. */
    }
    else if(! esl_opt_IsDefault(go, "--fil-T-qdb")) {
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
    if(fqdb_S > 0.9) do_qdb_filter = FALSE; /* turn off QDB filter if predicted survival fraction is crappy */
  }
  else do_qdb_filter = FALSE;

  /* C. add QDB filter, if necessary (before HMM filter, filters added like a stack, HMM filter is added last but used first) */
  if(do_qdb_filter) { 
    AddFilterToSearchInfo(cm, TRUE, FALSE, FALSE, FALSE, FALSE, fqdb_smx, NULL, fqdb_ctype, fqdb_sc, fqdb_E);
    /* DumpSearchInfo(cm->si); */
  }
  else if (fqdb_smx != NULL) cm_FreeScanMatrix(cm, fqdb_smx); 
  /* D. add HMM filter, if necessary (after QDB filter, filters added like a stack, HMM filter is added last but used first) */
  if (do_hmm_filter) { 
    AddFilterToSearchInfo(cm, FALSE, FALSE, FALSE, TRUE, FALSE, NULL, NULL, fhmm_ctype, fhmm_sc, fhmm_E);
    /* DumpSearchInfo(cm->si); */
  }
  ValidateSearchInfo(cm, cm->si);
  return eslOK;
}

/* Function: set_searchinfo_for_uncalibrated_cm)
 * Date:     EPN, Thu Mar  6 05:29:16 2008
 * 
 * Purpose:  For a CM WITHOUT exponential tail and filter thresholds statistics 
 *           from cmcalibrate: determine how many rounds of searching we will 
 *           do (all rounds but last round are filters), and set the relevant 
 *           info in the SearchInfo_t <cm->si> object, including cutoffs.
 *
 *************************************************************************
 * FIX EVERYTHING BELOW THIS LINE !!!!
 *************************************************************************
 * Filters:
 *
 * User can specify 0 or 1 round of HMM filtering with Forward, followed by 0 or 
 * 1 round of CM filtering with QDB CYK. This is determined as described below. 
 * The default behavior depends on whether or not the CM file has exp tail and 
 * filter threshold information from cmcalibrate.
 * 
 * HMM filtering:
 * ('none' below means none of --fil-no-hmm, --fil-T-hmm, --fil-S-hmm, --fil-Smax-hmm) 
 * 
 *                  exp tail &
 *                  filter thr 
 *   options        in cmfile?    filter/don't filter and cutoff determination
 *   -------       ----------    ------------------------------------------------
 * 1. none                 yes    filter (usually). automatically determine appropriate HMM 
 *                                filter cutoff for final round cutoff using cmfile's filter 
 *                                threshold info from cmcalibrate. Sometimes this info may
 *                                tell us it's inappropriate to use an HMM filter, in which
 *                                case, don't filter.
 * 
 * 2. none                  no    filter. set filter bit score cutoff to default value <x> for 
 *                                --fil-T-hmm (3.0 bits)
 *                               
 * 3. --fil-no-hmm      yes/no    don't filter.
 * 
 * 4. --fil-T-hmm <x>   yes/no    filter. set filter bit score cutoff as <x>
 * 
 * 5. --fil-S-hmm <x>      yes    filter. set E-value filter cutoff to E-value that gives
 *                                predicted survival fraction of <x>. NOTE: larger models,
 *                                with larger Ws, will have smaller E-value cutoffs for
 *                                same <x> (a somewhat unsettling behavior).
 * 
 * 6. --fil-S-hmm <x>       no    die. we can't deal with this combo, tell the user.
 * 
 * 7. --fil-Smax-hmm <x>   yes    filter, same as in case 1 above, but if the 'appropriate'
 *                                CM E-value cutoff gives a predicted survival fraction > <x>,
 *                                reset it as the E-value cutoff that gives <x> (even if 
 *                                cmfile told us to turn filtering off).
 * 
 * 8. --fil-Smax-hmm <x>    no    die. we can't deal with this combo, tell the user.
 *
 * 
 * CM filtering: 
 * ('none' below means none of --fil-no-qdb, --fil-T-qdb, --fil-S-qdb) 
 *                               
 *                                final 
 *                  exp tail &    round 
 *                  filter thr    cutoff
 *   options        in cmfile?    type    filter/don't filter and cutoff determination
 *   -------       -----------    ------  -----------------------------------------
 * 1. none                 yes    E val   filter. set E-value filter cutoff to E-value that gives
 *                                        predicted survival fraction of <x>. NOTE: larger models,
 *                                        with larger Ws, will have smaller E-value cutoffs for
 *                                        same <x> (a somewhat unsettling behavior).
 * 
 * 2. none                  no    filter. set filter bit score cutoff to default value <x> for 
 *                                --fil-T-qdb (0.0 bits)
 *                               
 * 3. --fil-no-qdb      yes/no    don't filter. 
 * 
 * 4. --fil-T-qdb <x>   yes/no    filter. set filter bit score cutoff as <x>
 * 
 * 5. --fil-S-qdb <x>      yes    filter. set E-value filter cutoff to E-value that gives
 *                                predicted survival fraction of <x>. NOTE: larger models,
 *                                with larger Ws, will have smaller E-value cutoffs for
 *                                same <x> (a somewhat unsettling behavior).
 * 
 * 6. --fil-S-qdb <x>       no    die. we can't deal with this combo, tell the user.
 *
 * ('none' below means none of --fil-no-qdb, --fil-T-qdb, --fil-S-qdb) 
 * 
 * beta, the tail loss probability for the QDB calculation is determined as follows:
 * If --fil-beta <x> is not enabled, cm->beta_qdb is set as cm->beta_W as read in the cmfile
 *               and initially set by cmbuild.
 * If --fil-beta <x> is enabled, <x> is used as cm->beta_qdb
 *
 * Final round related options (after all filtering is complete):
 *
 * --inside:     search with CM inside (TRUE by default)
 * --cyk:        search with CM CYK 
 * --forward:    search with HMM Forward
 * --viterbi:    search with HMM Viterbi
 * -T:           CM/HMM bit score threshold
 * -E:           CM/HMM E-value threshold (requires exp tail info in CM file)
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
  if((! use_hmmonly) && (stype != SEARCH_WITH_CM)) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_searchinfo_for_uncalibrated_cm(), search_opts for final round of search does not have HMMVITERBI or HMMFORWARD flags raised, but is not of type SEARCH_WITH_CM.");
  /* determine configuration of CM and CP9 HMM based on cm->flags & cm->search_opts */
  CM2ExpMode(cm, search_opts, &cm_mode, &hmm_mode); 

  /* set up final round cutoff, either 0 or 1 of 5 options is enabled. 
   * esl_opt_IsDefault() returns FALSE even if option is enabled with default value.
   */
  if(esl_opt_IsDefault(go, "-E") && 
     esl_opt_IsDefault(go, "-T") && 
     esl_opt_IsDefault(go, "--ga") && 
     esl_opt_IsDefault(go, "--tc") && 
     esl_opt_IsDefault(go, "--nc")) { 
    /* No relevant options enabled, cutoff is default bit score cutoff */
    final_sc    = esl_opt_GetReal(go, "-T");
  }
  else if(! esl_opt_IsDefault(go, "-E")) { /* -E enabled, error b/c we don't have exp tail stats */
    ESL_FAIL(eslEINVAL, errbuf, "-E requires exp tail statistics in <cm file>. Use cmcalibrate to get exp tail stats.");
  }
  else if(! esl_opt_IsDefault(go, "-T")) { /* -T enabled, use that */
    final_sc    = esl_opt_GetReal(go, "-T");
  }
  else if(! esl_opt_IsDefault(go, "--ga")) { /* --ga enabled, use that, if available, else die */
    if(use_hmmonly) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "--ga is incompatible with --viterbi and --forward.");
    if(! (cm->flags & CMH_GA)) ESL_FAIL(eslEINVAL, errbuf, "No GA gathering threshold in CM file, can't use --ga.");
    final_sc    = cm->ga;
  }
  else if(! esl_opt_IsDefault(go, "--tc")) { /* --tc enabled, use that, if available, else die */
    if(use_hmmonly) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "--tc is incompatible with --viterbi and --forward.");
    if(! (cm->flags & CMH_TC)) ESL_FAIL(eslEINVAL, errbuf, "No TC trusted cutoff in CM file, can't use --tc.");
    final_sc    = cm->tc;
  }
  else if(! esl_opt_IsDefault(go, "--nc")) { /* --nc enabled, use that, if available, else die */
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
    if(esl_opt_IsDefault(go, "--fil-S-hmm") && esl_opt_IsDefault(go, "--fil-T-hmm")) { 
      /* No relevant options selected, set cutoff as default HMM filter bit score. */
      fhmm_sc    = esl_opt_GetReal(go, "--fil-T-hmm");
    } /* end of if(esl_opt_IsDefault(go, "--fil-S-hmm") && esl_opt_IsDefault(go, "--fil-T-hmm")) */
    else if(! esl_opt_IsDefault(go, "--fil-S-hmm")) { /* can't deal with this b/c we don't have exp tail stats */
      ESL_FAIL(eslEINVAL, errbuf, "--fil-S-hmm requires exp tail statistics in <cm file>. Use cmcalibrate to get exp tail stats.");
    }
    else if(! esl_opt_IsDefault(go, "--fil-T-hmm")) { /* HMM filter bit score threshold was set on command line, use that */
      fhmm_sc    = esl_opt_GetReal(go, "--fil-T-hmm");
    }
    else ESL_FAIL(eslEINCONCEIVABLE, errbuf, "No HMM filter cutoff selected. This shouldn't happen.");
  }

  /* B. determine thresholds/stats for CM QDB filter to add */
  qdb_mode = (ExpModeIsLocal(cm_mode)) ? EXP_CM_LC : EXP_CM_GC; /* always do CYK with QDB filter, only question is local or glocal? */
  fqdb_sc = -1.;
  if(do_qdb_filter && esl_opt_GetBoolean(go, "--fil-qdb")) { /* determine thresholds, beta for qdb cyk filter */
    /* build the ScanMatrix_t for the QDB filter round, requires calcing dmin, dmax */
    if(! esl_opt_IsDefault(go, "--fil-beta")) fqdb_beta_qdb = esl_opt_GetReal(go, "--fil-beta");
    else fqdb_beta_qdb = cm->beta_W; /* use beta used to calc W in CM file by default */
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
    
    if(esl_opt_IsDefault(go, "--fil-S-qdb") && esl_opt_IsDefault(go, "--fil-T-qdb")) {
      /* No relevant options selected, set CM QDB filter as default bit score. */
      fqdb_sc    = esl_opt_GetReal(go, "--fil-T-qdb");
    }
    else if(! esl_opt_IsDefault(go, "--fil-S-qdb")) { /* can't deal with this b/c we don't have exp tail stats */
      ESL_FAIL(eslEINVAL, errbuf, "--fil-S-qdb requires exp tail statistics in <cm file>. Use cmcalibrate to get exp tail stats.");
    }
    else if(! esl_opt_IsDefault(go, "--fil-T-qdb")) { /* CM CYK bit score threshold was set on command line, use that */
      fqdb_sc    = esl_opt_GetReal(go, "--fil-T-qdb");
    }
    else ESL_FAIL(eslEINCONCEIVABLE, errbuf, "No CM filter cutoff selected. This shouldn't happen.");
  }
  else do_qdb_filter = FALSE;

  /* C. add QDB filter, if necessary (before HMM filter, filters added like a stack, HMM filter is added last but used first) */
  if(do_qdb_filter) { 
    AddFilterToSearchInfo(cm, TRUE, FALSE, FALSE, FALSE, FALSE, fqdb_smx, NULL, SCORE_CUTOFF, fqdb_sc, -1.);
    /* DumpSearchInfo(cm->si); */
  }
  else if (fqdb_smx != NULL) cm_FreeScanMatrix(cm, fqdb_smx); 
  /* D. add HMM filter, if necessary (after QDB filter, filters added like a stack, HMM filter is added last but used first) */
  if (do_hmm_filter) { 
    AddFilterToSearchInfo(cm, FALSE, FALSE, FALSE, TRUE, FALSE, NULL, NULL, SCORE_CUTOFF, fhmm_sc, -1.);
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
  fprintf(stdout, "%-13s %ld\n",  "# dbsize(nt):", cfg->dbsize);

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
int print_searchinfo_for_calibrated_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, FILE *fp, CM_t *cm, float *cm_surv_fractA, int *cm_nhitsA, double in_asec, double in_total_psec, double *ret_total_psec)
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
  float avg_hit_len;        /* average length of a hit, from qdb calculation */
  int   use_qdb;            /* are we using qdb for current round? */
  char  time_buf[128];	    /* for printing predicted times */
  int   do_pad;             /* TRUE to add W pad onto survival fraction (FALSE only in final round of searching) */
  int   pre_search_mode;    /* TRUE if this function was called before the search was run */
  double sec_per_res;        /* seconds required to search 1 residue with current model, current round of searching */
  double psec;              /* predicted seconds for current cm, current round */
  double total_psec = 0.;   /* predicted number of seconds for full round */


  pre_search_mode = (cm_surv_fractA == NULL && cm_nhitsA == NULL) ? TRUE : FALSE;
  if(pre_search_mode && ret_total_psec == NULL) ESL_FAIL(status, errbuf, "print_searchinfo_for_calibrated_cm, pre-search mode, but ret_total_psec is NULL.");

  ESL_RANDOMNESS *r = esl_randomness_Create((long) 33);
  if(r == NULL) ESL_FAIL(status, errbuf, "print_searchinfo_for_calibrated_cm, memory error, couldn't create randomness object.");

  /* Could use ESL_GETOPTS here, but using the CM flags assures we're reporting
   * on how the CM is actually config'ed, not how we want it to be
   */

  if(! (cm->flags & CMH_EXPTAIL_STATS)) ESL_FAIL(eslEINCOMPAT, errbuf, "print_searchinfo_for_calibrated_cm(): cm: %s does not have exp tail stats.", cm->name);
  using_filters = (cm->si->nrounds > 0) ? TRUE : FALSE;

  if(pre_search_mode) fprintf(stdout, "#\n# Pre-search info for CM %d: %s\n", cfg->ncm, cm->name);
  else                fprintf(stdout, "#\n# Post-search info for CM %d: %s\n", cfg->ncm, cm->name);

  fprintf(stdout, "#\n");

  if((status = cm_GetAvgHitLen        (cm,      errbuf, &avg_hit_len))        != eslOK) return status;

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
    surv_fract = E2SurvFract(e_cutoff, (stype == SEARCH_WITH_CM) ? smx->W : cm->W, avg_hit_len, cfg->dbsize, do_pad);

    if(pre_search_mode) { 
      if((status = estimate_search_time_for_round(go, cfg, errbuf, cm, stype, search_opts, smx, r, &sec_per_res)) != eslOK) return status;
      psec = prv_surv_fract * (double) cfg->dbsize * sec_per_res;
      /*printf("psec: %f\nprv_surv_fract: %f\n", psec, prv_surv_fract);*/
#ifdef HAVE_MPI
    /* if we're paralellized take that into account */
    if(esl_opt_GetBoolean(go, "--mpi") && cfg->nproc > 1) psec /= (cfg->nproc-1);
#endif
    /* if we're only forecasting the time, divide by <n> from --forecast <n>, which is theoretical number of processors */
    if(! esl_opt_IsDefault(go, "--forecast")) { 
      if(esl_opt_GetInteger(go, "--forecast") > 1) { 
	psec /= (esl_opt_GetInteger(go, "--forecast") - 1);
      }
    }
    total_psec += psec;

    FormatTimeString(time_buf, psec, TRUE);
      if(n == 0) { 
	fprintf(stdout, "# %3s  %3s  %3s  %3s  %5s  %19s  %20s\n",               ""    , "",    "",    "",    "",  "      cutoffs      ",   "    predictions     ");
	fprintf(stdout, "# %3s  %3s  %3s  %3s  %5s  %19s  %20s\n",               "",     "",    "",    "",    "",  "-------------------",   "--------------------");
	fprintf(stdout, "# %3s  %3s  %3s  %3s  %5s  %10s  %7s  %7s  %11s\n", "rnd",  "mod", "alg", "cfg", "beta",  "E value",    "bit sc",  "surv",    "run time");
	fprintf(stdout, "# %3s  %3s  %3s  %3s  %5s  %10s  %7s  %7s  %11s\n", "---",  "---", "---", "---", "-----", "----------", "-------", "-------", "-----------");
      }
      
      fprintf(stdout, "  %3d", (n+1));
      if(stype == SEARCH_WITH_CM) { 
	fprintf(stdout, "  %3s  %3s  %3s  ", "cm", ((search_opts & CM_SEARCH_INSIDE) ? "ins" : "cyk"), ((cm->flags & CMH_LOCAL_BEGIN) ? "loc" : "glc"));
	use_qdb  = (smx->dmin == NULL && smx->dmax == NULL) ? FALSE : TRUE;
	if(use_qdb) fprintf(stdout, "%5g", smx->beta_qdb);
	else        fprintf(stdout, "%5s", "-");
      }
      else { 
	fprintf(stdout, "  %3s  %3s  %3s  %5s", "hmm", ((search_opts & CM_SEARCH_HMMFORWARD) ? "fwd" : "vit"), ((cm->cp9->flags & CPLAN9_LOCAL_BEGIN) ? "loc" : "glc"), "-");
      }
      if(e_cutoff < -0.1)  if((status = Score2MaxE(cm, errbuf, exp_mode, sc_cutoff, &e_cutoff)) != eslOK) return status;
      if(e_cutoff < 0.01)  fprintf(stdout, "  %10.1e", e_cutoff);
      else                 fprintf(stdout, "  %10.3f", e_cutoff);
      
      fprintf(stdout, "  %7.2f", sc_cutoff);
      if(surv_fract < 0.0001) fprintf(stdout, "  %7.1e", surv_fract);
      else                    fprintf(stdout, "  %7.4f", surv_fract);
      fprintf(stdout, "  %11s\n", time_buf);

      if(n != 0 && n == cm->si->nrounds) { /* print total expected run time */
	fprintf(stdout, "# %3s  %3s  %3s  %3s  %5s  %10s  %7s  %7s  %11s\n", "---",  "---", "---", "---", "-----", "----------", "-------", "-------", "-----------");
	FormatTimeString(time_buf, total_psec, TRUE);
	fprintf(stdout, "  %3s  %3s  %3s  %3s  %5s  %10s  %7s  %7s  %11s\n", "all",  "-",   "-",   "-",   "-",     "-",          "-",       "-",       time_buf);
      }
    }
    else { /* search is done */
      if(n == 0) { 
	fprintf(stdout, "# %3s  %3s  %3s  %3s  %5s  %19s  %17s\n",               ""    , "",    "",    "",    "",  "  number of hits   ",   "  surv fraction  ");
	fprintf(stdout, "# %3s  %3s  %3s  %3s  %5s  %19s  %17s\n",               "",     "",    "",    "",    "",  "-------------------",   "-----------------");
	fprintf(stdout, "# %3s  %3s  %3s  %3s  %5s  %10s  %7s  %8s  %7s\n", "rnd",  "mod", "alg", "cfg", "beta",  "expected",    "actual", "expected", "actual");
	fprintf(stdout, "# %3s  %3s  %3s  %3s  %5s  %10s  %7s  %8s  %7s\n", "---",  "---", "---", "---", "-----", "----------", "-------", "--------", "-------");
      }
      
      fprintf(stdout, "  %3d", (n+1));
      if(stype == SEARCH_WITH_CM) { 
	fprintf(stdout, "  %3s  %3s  %3s  ", "cm", ((search_opts & CM_SEARCH_INSIDE) ? "ins" : "cyk"), ((cm->flags & CMH_LOCAL_BEGIN) ? "loc" : "glc"));
	use_qdb  = (smx->dmin == NULL && smx->dmax == NULL) ? FALSE : TRUE;
	if(use_qdb) fprintf(stdout, "%5g", smx->beta_qdb);
	else        fprintf(stdout, "%5s", "-");
      }
      else { 
	fprintf(stdout, "  %3s  %3s  %3s  %5s", "hmm", ((search_opts & CM_SEARCH_HMMFORWARD) ? "fwd" : "vit"), ((cm->cp9->flags & CPLAN9_LOCAL_BEGIN) ? "loc" : "glc"), "-");
      }
      if(e_cutoff < -0.1)  if((status = Score2MaxE(cm, errbuf, exp_mode, sc_cutoff, &e_cutoff)) != eslOK) return status;
      if(e_cutoff < 0.01)  fprintf(stdout, "  %10.1e", e_cutoff);
      else                 fprintf(stdout, "  %10.3f", e_cutoff);
      
      fprintf(stdout, "  %7d", cm_nhitsA[n]);
      if(surv_fract < 0.0001) fprintf(stdout, "  %8.1e", surv_fract);
      else                    fprintf(stdout, "  %8.4f", surv_fract);
      if(cm_surv_fractA[n] < 0.0001) fprintf(stdout, "  %7.1e\n", cm_surv_fractA[n]);
      else                           fprintf(stdout, "  %7.4f\n", cm_surv_fractA[n]);

      if(n == cm->si->nrounds) { /* print total expected run time */
	fprintf(stdout, "#\n");
	/*fprintf(stdout, "# %24s\n",       "   total search time    ");*/
	/*fprintf(stdout, "# %24s\n",       "------------------------");*/
	fprintf(stdout, "# %13s  %13s\n", "expected time",  "actual time");
	fprintf(stdout, "# %13s  %13s\n", "-------------",  "-------------");
	FormatTimeString(time_buf, in_total_psec, TRUE);
	fprintf(stdout, "  %13s", time_buf);
	FormatTimeString(time_buf, in_asec, TRUE);
	fprintf(stdout, "  %13s\n", time_buf);
      }
    }
  }
  fflush(fp);
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
int print_searchinfo_for_uncalibrated_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, FILE *fp, CM_t *cm, float *cm_surv_fractA, int *cm_nhitsA, double in_asec)
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

  if(pre_search_mode) fprintf(stdout, "#\n# Pre-search info for CM %d: %s\n", cfg->ncm, cm->name);
  else                fprintf(stdout, "#\n# Post-search info for CM %d: %s\n", cfg->ncm, cm->name);
  fprintf(stdout, "#\n");

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
      }
      else { /* post search */
	fprintf(stdout, "# %3s  %3s  %3s  %3s  %5s  %10s  %8s  %10s\n", "rnd",  "mod", "alg", "cfg", "beta",  "bit sc cut", "num hits", "surv fract");
	fprintf(stdout, "# %3s  %3s  %3s  %3s  %5s  %10s  %8s  %10s\n", "---",  "---", "---", "---", "-----", "----------", "--------", "----------");
      }
    }
    fprintf(stdout, "  %3d", (n+1));
    if(stype == SEARCH_WITH_CM) { 
      fprintf(stdout, "  %3s  %3s  %3s  ", "cm", ((search_opts & CM_SEARCH_INSIDE) ? "ins" : "cyk"), ((cm->flags & CMH_LOCAL_BEGIN) ? "loc" : "glc"));
      if(use_qdb) fprintf(stdout, "%5g", smx->beta_qdb);
      else        fprintf(stdout, "%5s", "-");
    }
    else { 
      fprintf(stdout, "  %3s  %3s  %3s  %5s", "hmm", ((search_opts & CM_SEARCH_HMMFORWARD) ? "fwd" : "vit"), ((cm->cp9->flags & CPLAN9_LOCAL_BEGIN) ? "loc" : "glc"), "-");
    }
    fprintf(stdout, "  %10.2f", sc_cutoff);
    if(pre_search_mode) fprintf(stdout, "\n");
    else { /* post search */
      if(cm_surv_fractA[n] < 0.0001) fprintf(stdout, "  %8d  %10.1e\n", cm_nhitsA[n], cm_surv_fractA[n]);
      else                           fprintf(stdout, "  %8d  %10.4f\n", cm_nhitsA[n], cm_surv_fractA[n]);
      if(n == cm->si->nrounds) { 
	fprintf(stdout, "#\n");
	fprintf(stdout, "# %11s\n", "run time");
	fprintf(stdout, "# %11s\n", "-----------");
	FormatTimeString(time_buf, in_asec, FALSE);
	fprintf(stdout, "  %11s\n", time_buf);
      }
    }
  }
  fflush(fp);
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

  ESL_DSQ *dsq;
  ESL_STOPWATCH *w  = esl_stopwatch_Create();

  if(w == NULL)               ESL_FAIL(status, errbuf, "estimate_search_time_for_round(): memory error, stopwatch not created.\n");
  if(ret_sec_per_res == NULL) ESL_FAIL(status, errbuf, "estimate_search_time_for_round(): ret_sec_per_res is NULL");

  if(stype == SEARCH_WITH_CM) {
    if(smx == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "estimate_search_time_for_round(), stype is SEARCH_WITH_CM, but smx is NULL");
    int use_qdb     = (smx->dmin == NULL && smx->dmax == NULL) ? FALSE : TRUE;
    if(use_qdb) { if((status = cm_GetNCalcsPerResidueForGivenBeta(cm, errbuf, FALSE, smx->beta_qdb, &Mc_per_res, &irrelevant_W))  != eslOK) return status; }
    else        { if((status = cm_GetNCalcsPerResidueForGivenBeta(cm, errbuf, TRUE,  cm->beta_W,    &Mc_per_res, &irrelevant_W))  != eslOK) return status; }
    psec_per_Mc = (search_opts & CM_SEARCH_INSIDE) ? (1. /  75.) : (1. / 275.);  /*  75 Mc/S inside;  275 Mc/S CYK */
    /* determine L that will take about targ_sec seconds */
    L = targ_sec / (psec_per_Mc * Mc_per_res);
    L = ESL_MAX(L, (int) (((float) cm->W)/5.)); /* we have to search at least (cm->W/5.) residues */
    /* now determine exactly how many dp calculations we'd do if we search L residues, 
     * this won't be the same as Mc_per_res * L b/c Mc_per_res was passed in from caller
     * and was calculated after correcting for the fact that the first W residues have fewer
     * DP calcs than all other residues, b/c d < W for first W residues.
     */  
    if((status = cm_CountSearchDPCalcs(cm, errbuf, L, NULL, NULL, smx->W, FALSE,  NULL, &Mc)) != eslOK) return status;
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
  esl_rnd_xfIID(r, cm->null, cm->abc->K, L, dsq);

  esl_stopwatch_Start(w);
  if(stype == SEARCH_WITH_CM) { 
    if(search_opts & CM_SEARCH_INSIDE) { 
      if((status = FastIInsideScan(cm, errbuf, smx, dsq, 1, L, 0., NULL, NULL, NULL)) != eslOK) return status;
    }
    else 
      if((status = FastCYKScan(cm, errbuf, smx, dsq, 1, L, 0., NULL, NULL, NULL)) != eslOK) return status;
  }
  else { /* search with HMM */
    if(search_opts & CM_SEARCH_HMMFORWARD) { /* forward */
      if((status = cp9_Forward(cm, errbuf, cm->cp9_mx, dsq, 1, L, cm->W, 0., NULL,
			       TRUE,   /* we're scanning */
			       FALSE,  /* we're not ultimately aligning */
			       TRUE,   /* be memory efficient */
			       NULL, NULL, NULL)) != eslOK) return status;
    }
    else { /* viterbi */
      if((status = cp9_Viterbi(cm, errbuf, cm->cp9_mx, dsq, 1, L, cm->W, 0., NULL,
			     TRUE,   /* we're scanning */
			     FALSE,  /* we're not ultimately aligning */
			     TRUE,   /* be memory efficient */
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

  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "estimate_search_time_for_round(): memory error.\n");
  return status; /* NEVERREACHED */
}

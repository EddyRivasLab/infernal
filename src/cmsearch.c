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

#define STRATOPTS  "--cmonly,--hmmfilter,--hmmviterbi,--hmmforward" /* exclusive choice for search strategy */
#define ALPHOPTS   "--rna,--dna"                                    /* exclusive choice for output alphabet */

#define I_CMCUTOPTS   "-E,-T,--ga,--tc,--nc,--hmmviterbi,--hmmforward" /* exclusive choice for CM cutoff */
#define I_HMMCUTOPTS1 "--hmmthr,--hmmcalcthr,--hmmE,--hmmT"            /* exclusive choice for HMM cutoff set 1 */
#define I_HMMCUTOPTS2 "--hmmthr,--hmmcalcthr,--hmmE,--hmmT,--hmmviterbi,--hmmforward"  /* exclusive choice for HMM cutoff set 2 */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs       incomp  help  docgroup*/
  /* basic options */
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",   1 },
  { "-g",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--hmmviterbi,--hmmforward", "configure CM for glocal alignment [default: local]", 1 },
  { "-i",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--hmmviterbi,--hmmforward", "use scanning Inside algorithm instead of CYK", 1 },
  { "--informat",eslARG_STRING, NULL,  NULL, NULL,      NULL,      NULL,        NULL, "specify the input file is in format <x>, not FASTA", 1 },
  { "--toponly", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "only search the top strand", 1 },
  { "--window",  eslARG_INT,    NULL,  NULL, "n>0",     NULL,      NULL,        NULL, "set scanning window size to <n> [default: calculated]", 1 },
  { "--null2",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--hmmviterbi,--hmmforward", "turn on the post hoc second null model", 1 },
  { "--iins",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "allow informative insert emissions, do not zero them", 1 },
  { "--rtrans",  eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--hmmviterbi,--hmmforward", "replace CM transition scores from <cmfile> with RSEARCH scores", 1 },
  { "--greedy",  eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--hmmviterbi,--hmmforward", "resolve overlapping hits with a greedy algorithm a la RSEARCH", 1 },
  /* strategy choice */
  { "--cmonly",   eslARG_NONE,"default",NULL,NULL,      STRATOPTS, NULL,        NULL, "search only with CM, do not filter [default]", 2 },
  { "--hmmfilter",eslARG_NONE,  FALSE, NULL, NULL,      STRATOPTS, NULL,        NULL, "subseqs j-W+1..i+W-1 survive (j=end from Fwd, i=start from Bwd)", 2 },
  { "--hmmviterbi",eslARG_NONE, FALSE, NULL, NULL,      STRATOPTS, NULL,        NULL, "do not use CM at all, just scan with HMM (Viterbi alg)", 2 },
  { "--hmmforward",eslARG_NONE, FALSE, NULL, NULL,      STRATOPTS, NULL,        NULL, "do not use CM at all, just scan with HMM (Forward alg)", 2 },
  /* CM cutoff options */
  { "-E",        eslARG_REAL,   "0.1", NULL, "x>0.",    NULL,      NULL,  I_CMCUTOPTS, "use cutoff E-value of <x> for CM search", 3 },
  { "-T",        eslARG_REAL,   "0.0", NULL, NULL,      NULL,      NULL,  I_CMCUTOPTS, "use cutoff bit score of <x> for CM search", 3 },
  { "--ga",      eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,  I_CMCUTOPTS, "use CM Rfam GA gathering threshold as cutoff bit score", 3 },
  { "--tc",      eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,  I_CMCUTOPTS, "use CM Rfam TC trusted cutoff as cutoff bit score", 3 },
  { "--nc",      eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,  I_CMCUTOPTS, "use CM Rfam NC noise cutoff as cutoff bit score", 3 },
  /* HMM cutoff options */
  { "--hmmthr",  eslARG_NONE,   FALSE, NULL, NULL,    NULL,"--hmmfilter",I_HMMCUTOPTS1,"use HMM filter from cmcalibrate (in <cm file>)", 4 },
  { "--hmmcalcthr",eslARG_NONE, FALSE, NULL, NULL,    NULL,"--hmmfilter",I_HMMCUTOPTS1,"calculate HMM filter threshold by sampling from CM", 4 },
  { "--hmmE",    eslARG_REAL,   "50.", NULL, "x>0.",    NULL,      NULL, I_HMMCUTOPTS2,"use cutoff E-value of <x> for CP9 HMM filter/search", 4 },
  { "--hmmT",    eslARG_REAL,   "0.0", NULL, NULL,      NULL,      NULL, I_HMMCUTOPTS2,"use cutoff bit score of <x> for CP9 HMM filter/search", 4 },
   /* QDB related options */
  { "--beta",    eslARG_REAL,   "1e-7",NULL, "x>0",     NULL,      NULL, "--hmmviterbi,--hmmforward", "set tail loss prob for QDB and window size calculation to <x>", 5 },
  { "--noqdb",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--hmmviterbi,--hmmforward", "DO NOT use query dependent banding (QDB) for acceleration", 5 },
  { "--qdbfile", eslARG_STRING, NULL,  NULL, NULL,      NULL,      NULL,"--hmmviterbi,--hmmforward,--noqdb","read QDBs from file <s> (outputted from cmbuild)", 5 },
  /* HMM filtering options */
  { "--hbanded", eslARG_NONE,   FALSE, NULL, NULL,      NULL,       NULL,       NULL, "calculate and use HMM bands in CM search", 6 },
  { "--tau",     eslARG_REAL,   "1e-7",NULL, "0<x<1",   NULL,"--hbanded",       NULL, "set tail loss prob for --hbanded to <x>", 6 },
  { "--scan2bands",eslARG_NONE, FALSE, NULL, NULL,      NULL,"--hbanded",       NULL, "derive HMM bands from scanning Forward/Backward", 6 },
  { "--sums",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--hbanded",       NULL, "use posterior sums during HMM band calculation (widens bands)", 6 },
  /* HMM configuration options */
  { "--hmmglocal",eslARG_NONE,  FALSE, NULL, NULL,      NULL,      NULL,  "--cmonly", "configure HMM for glocal alignment [default: local]", 7 },
  { "--hmmnoel",  eslARG_NONE,  FALSE, NULL, NULL,      NULL,      NULL,"-g,--cmonly,--hmmglocal", "DO NOT enable HMM EL local ends that mirror CM", 7 },
  { "--hmmgreedy",eslARG_NONE,  FALSE, NULL, NULL,      NULL,      NULL,  "--cmonly", "resolve HMM overlapping hits with a greedy algorithm a la RSEARCH", 7 },
  { "--hmmrescan",eslARG_NONE,  FALSE, NULL, NULL,      NULL,      NULL,  "--cmonly", "rescan subseq hits w/Forward (auto ON if --enfseq)", 7 },
  /* filter threshold calculation options */
  { "--seed",    eslARG_INT,    FALSE, NULL, "n>0",     NULL,"--hmmcalcthr",    NULL, "set random number generator seed to <n>", 8 },
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
  { "--stall",   eslARG_NONE,  FALSE, NULL, NULL,       NULL,      NULL,    NULL, "arrest after start: for debugging MPI under gdb",   13 },  
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
static int process_search_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, ESL_DSQ *dsq, 
				   int L, float min_cm_cutoff, float min_cp9_cutoff, search_results_t **ret_results);
static int initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int set_gumbels(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int set_cutoffs(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, int *ret_cm_mode, int *ret_cp9_mode,
		       float *ret_min_cm_cutoff, float *ret_min_cp9_cutoff, int *ret_using_e_cutoff);
static int set_window(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int calc_filter_threshold(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, float *ret_Smin);
static int print_search_info(FILE *fp, CM_t *cm, int cm_mode, int cp9_mode, long N, char *errbuf);

static int read_qdb_file(FILE *fp, CM_t *cm, int *dmin, int *dmax);
static int is_integer(char *s);
#if HAVE_MPI
static int determine_cm_min_max_chunksize(struct cfg_s *cfg, CM_t *cm, int *ret_min_chunksize, int *ret_max_chunksize);
static int determine_seq_chunksize(struct cfg_s *cfg, int L, int min_chunksize, int chunksize);
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
      puts("\nHMM filtering options: (require --hmmfilter)");
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
  cfg.fmt        = eslSQFILE_UNKNOWN;      /* autodetect sequence file format by default. */ 
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

      if(cfg.nproc == 1) cm_Fail("ERROR, MPI mode, but only 1 processor running...");

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
  if (! esl_opt_GetBoolean(go, "--toponly")) cfg->N *= 2;

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
  int            cm_mode  = -1;   /* CM algorithm mode                        */
  int            cp9_mode = -1;   /* CP9 algorithm mode                       */
  float          Smin;
  float          min_cm_cutoff;
  float          min_cp9_cutoff;
  int            using_e_cutoff;
  int            rci;
  dbseq_t       *dbseq = NULL;

  if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);
  /*if ((status = init_shared_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);*/

  while (CMFileRead(cfg->cmfp, &(cfg->abc), &cm))
    {
      if (cm == NULL) cm_Fail("Failed to read CM from %s -- file corrupt?\n", cfg->cmfile);
      cfg->ncm++;

      /* initialize the flags/options/params and configuration of the CM */
      if((  status = initialize_cm(go, cfg, errbuf, cm))                    != eslOK) cm_Fail(errbuf);
      if((  status = CreateCMConsensus(cm, cfg->abc_out, 3.0, 1.0, &cons))  != eslOK) cm_Fail(errbuf);
      if(cm->flags & CMH_GUMBEL_STATS) 
	if((status = set_gumbels(go, cfg, errbuf, cm))                      != eslOK) cm_Fail(errbuf);
      if((  status = set_cutoffs(go, cfg, errbuf, cm, &cm_mode, &cp9_mode, &min_cm_cutoff,
				 &min_cp9_cutoff, &using_e_cutoff))         != eslOK) cm_Fail(errbuf);
      if((  status = set_window(go, cfg, errbuf, cm))                       != eslOK) cm_Fail(errbuf);
      if(esl_opt_GetBoolean(go, "--hmmcalcthr"))
	if((status = calc_filter_threshold(go, cfg, errbuf, cm, &Smin))     != eslOK) cm_Fail(errbuf);
      print_search_info(stdout, cm, cm_mode, cp9_mode, cfg->N, errbuf);

      while ((status = read_next_search_seq(cfg->abc, cfg->sqfp, cfg->do_rc, &dbseq)) == eslOK)
	{
	  for(rci = 0; rci <= cfg->do_rc; rci++) {
	    /*printf("SEARCHING >%s %d\n", dbseq->sq[reversed]->name, reversed);*/
	    if ((status = process_search_workunit(go, cfg, errbuf, cm, dbseq->sq[rci]->dsq, dbseq->sq[rci]->n, 
						  min_cm_cutoff, min_cp9_cutoff, &dbseq->results[rci])) != eslOK) cm_Fail(errbuf);
	    remove_overlapping_hits(dbseq->results[rci], 1, dbseq->sq[rci]->n);
	    if(using_e_cutoff) remove_hits_over_e_cutoff(cm, dbseq->results[rci], dbseq->sq[rci], 
							 ((cm->search_opts & CM_SEARCH_HMMVITERBI) || (cm->search_opts & CM_SEARCH_HMMFORWARD))); /* HMM hits? */
	  }
	  print_results (cm, cfg->abc_out, cons, dbseq, cfg->do_rc,
			 ((cm->search_opts & CM_SEARCH_HMMVITERBI) || (cm->search_opts & CM_SEARCH_HMMFORWARD))); /* use HMM stats? */
	  for(rci = 0; rci <= cfg->do_rc; rci++) {
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
 * EPN: GOAL OF IMPLEMENTATION FOLLOWS IN LOWERCASE.
 * IT IS NOT YET ACHIEVED.
 * TO ACHIEVE WE'LL NEED ALL FUNCS CALLED BY MPI TO
 * RETURN CLEANLY ALWAYS - BIG TASK TO REWRITE THOSE.
 * CURRENTLY NEARLY ALL ERRORS ARE UNRECOVERABLE, BUT THESE
 * ARE NOT LIMITED TO MPI COMMUNICATION ERRORS.
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
  float    min_cm_cutoff;
  float    min_cp9_cutoff;
  int      using_e_cutoff; 

  CM_t *cm;
  CMConsensus_t *cons = NULL;     /* precalculated consensus info for display purposes */
  int            cm_mode  = -1;   /* CM algorithm mode                        */
  int            cp9_mode = -1;   /* CP9 algorithm mode                       */
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
  int min_chunksize;
  int max_chunksize;
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
  if (xstatus == eslOK) { if ((silist = malloc(sizeof(int) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((rclist = malloc(sizeof(int) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((seqposlist = malloc(sizeof(int) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((lenlist = malloc(sizeof(int) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((dbseqlist = malloc(sizeof(dbseq_t *) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  if (xstatus == eslOK) { if ((sentlist = malloc(sizeof(int) * cfg->nproc)) == NULL) { sprintf(errbuf, "allocation failed"); xstatus = eslEMEM; } }
  for (wi = 0; wi < cfg->nproc; wi++) 
  { 
    silist[wi] = rclist[wi] = seqposlist[wi] = lenlist[wi] = -1;
    dbseqlist[wi] = NULL;
    sentlist[wi] = FALSE;
  }
  MPI_Bcast(&xstatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (xstatus != eslOK) cm_Fail(errbuf);
  ESL_DPRINTF1(("MPI master is initialized\n"));

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

  while (CMFileRead(cfg->cmfp, &(cfg->abc), &cm))
    {
      cfg->ncm++;  
      ESL_DPRINTF1(("MPI master read CM number %d\n", cfg->ncm));
      if((status = cm_master_MPIBcast(cm, 0, MPI_COMM_WORLD, &buf, &bn)) != eslOK) cm_Fail("MPI broadcast CM failed.");
      
      /* initialize the flags/options/params of the CM */
      if((status   = initialize_cm(go, cfg, errbuf, cm))                    != eslOK) cm_Fail(errbuf);
      if((status   = CreateCMConsensus(cm, cfg->abc_out, 3.0, 1.0, &cons))  != eslOK) cm_Fail(errbuf);
      if(cm->flags & CMH_GUMBEL_STATS) 
	if((status = set_gumbels(go, cfg, errbuf, cm))                     != eslOK)  cm_Fail(errbuf);
      if((  status = set_cutoffs(go, cfg, errbuf, cm, &cm_mode, &cp9_mode, &min_cm_cutoff,
				 &min_cp9_cutoff, &using_e_cutoff))         != eslOK) cm_Fail(errbuf);
      if((  status = set_window(go, cfg, errbuf, cm))                       != eslOK) cm_Fail(errbuf);
      if(esl_opt_GetBoolean(go, "--hmmcalcthr"))
	if((status = calc_filter_threshold(go, cfg, errbuf, cm, &Smin))     != eslOK) cm_Fail(errbuf);
      print_search_info(stdout, cm, cm_mode, cp9_mode, cfg->N, errbuf);
      determine_cm_min_max_chunksize(cfg, cm, &min_chunksize, &max_chunksize);
      printf("min_chunksize: %d\n", min_chunksize);
      printf("max_chunksize: %d\n", max_chunksize);

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
		  in_rc = FALSE;
		  seqpos = 1;
		  
		  si = 0;
		  while(dbseqlist[si] != NULL) si++;
		  ESL_DASSERT1((si < cfg->nproc));
		  dbseqlist[si] = dbseq;
		  sentlist[si]  = FALSE;
		  have_work = TRUE;
		  chunksize = determine_seq_chunksize(cfg, dbseq->sq[0]->n, min_chunksize, max_chunksize);
		  printf("L: %d chunksize: %d\n", dbseq->sq[0]->n, chunksize);
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
		  si_recv = silist[wi];
		  ESL_DPRINTF1(("MPI master sees that the result buffer contains search results (si_recv:%d)\n", si_recv));
		  if ((status = cm_search_results_MPIUnpack(buf, bn, &pos, MPI_COMM_WORLD, &worker_results)) != eslOK) cm_Fail("search results unpack failed");
		  ESL_DPRINTF1(("MPI master has unpacked search results\n"));
		      
		  /* worker_results will be NULL if 0 results (hits) sent back */
		  int x;
		  if(worker_results != NULL) { 
		    /* add results to dbseqlist[si_recv]->results[rclist[wi]] */
		    for(x = 0; x < worker_results->num_results; x++)
		      { 
			assert(worker_results->data[x].tr != NULL);
			assert(worker_results->data[x].tr->n > 0);
		      }
		    AppendResults(worker_results, dbseqlist[si_recv]->results[rclist[wi]], seqposlist[wi]);
		    /* careful, dbseqlist[si_recv]->results[rclist[wi]] now points to the nodes in worker_results->data,
		     * don't free those (don't use FreeResults(worker_results)) */
		    free(worker_results);
		    worker_results = NULL;
		  }
		  for(x = 0; x < dbseqlist[si_recv]->results[rclist[wi]]->num_results; x++)
		    { 
		      assert(dbseqlist[si_recv]->results[rclist[wi]]->data[x].tr != NULL);
		      assert(dbseqlist[si_recv]->results[rclist[wi]]->data[x].tr->n > 0);
		    }
		  dbseqlist[si_recv]->chunks_sent--;
		  if(sentlist[si_recv] && dbseqlist[si_recv]->chunks_sent == 0)
		    {
		      for(rci = 0; rci <= cfg->do_rc; rci++) {
			remove_overlapping_hits(dbseqlist[si_recv]->results[rci], 1, dbseqlist[si_recv]->sq[rci]->n);
			if(using_e_cutoff) remove_hits_over_e_cutoff(cm, dbseqlist[si_recv]->results[rci], dbseqlist[si_recv]->sq[rci], 
								     ((cm->search_opts & CM_SEARCH_HMMVITERBI) || (cm->search_opts & CM_SEARCH_HMMFORWARD)); /* HMM hits? */

		      }					      
		      print_results(cm, cfg->abc_out, cons, dbseqlist[si_recv], cfg->do_rc, 
				    ((cm->search_opts & CM_SEARCH_HMMVITERBI) || (cm->search_opts & CM_SEARCH_HMMFORWARD)); /* use HMM stats (used HMM only) */
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
		}
	      nproc_working--;
	    }
	  
	  if (have_work)
	    {   
	      /* send new search job */
	      ESL_DPRINTF1(("MPI master is sending sequence to search to worker %d\n", wi));
	      len = (chunksize < (dbseqlist[si]->sq[0]->n - seqpos + 1)) ? chunksize : (dbseqlist[si]->sq[0]->n - seqpos + 1);
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
  
  if (xstatus != eslOK) cm_Fail(errbuf);
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
  int           cm_mode  = -1;   /* CM algorithm mode                        */
  int           cp9_mode = -1;   /* CP9 algorithm mode                       */
  /*float         Smin;*/
  /*MPI_Status  mpistatus;*/
  ESL_DSQ      *dsq = NULL;
  int           L;
  search_results_t *results = NULL;
  float         min_cm_cutoff;
  float         min_cp9_cutoff;
  int           using_e_cutoff;

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
      if(cm->flags & CMH_GUMBEL_STATS) 
	if((status = set_gumbels(go, cfg, errbuf, cm))                      != eslOK) goto ERROR;
      if((  status = set_cutoffs(go, cfg, errbuf, cm, &cm_mode, &cp9_mode, &min_cm_cutoff, 
				 &min_cp9_cutoff, &using_e_cutoff))         != eslOK) goto ERROR;
      if((  status = set_window(go, cfg, errbuf, cm))                       != eslOK) goto ERROR;
      
      /* print_search_info(stdout, cm, cm_mode, cp9_mode, cfg->N, errbuf); */
      
      while((status = cm_dsq_MPIRecv(0, 0, MPI_COMM_WORLD, &wbuf, &wn, &dsq, &L)) == eslOK)
	{
	  ESL_DPRINTF1(("worker %d: has received search job, length: %d\n", cfg->my_rank, L));
	  if ((status = process_search_workunit(go, cfg, errbuf, cm, dsq, L, min_cm_cutoff, min_cp9_cutoff,
						&results)) != eslOK) goto ERROR;
	  ESL_DPRINTF1(("worker %d: has gathered search results\n", cfg->my_rank));
	  
	  n = 0;
	  if (cm_search_results_MPIPackSize(results, MPI_COMM_WORLD, &sz) != eslOK) goto ERROR; n += sz;  
	  if (n > wn) {
	    void *tmp;
	    ESL_RALLOC(wbuf, tmp, sizeof(char) * n);
	    wn = n;
	  }
	  ESL_DPRINTF1(("worker %d: has calculated the search results will pack into %d bytes\n", cfg->my_rank, n));
	  status = eslOK;

	  pos = 0;
	  if (cm_search_results_MPIPack(results, wbuf, wn, &pos, MPI_COMM_WORLD) != eslOK) goto ERROR;
	  MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);
	  ESL_DPRINTF1(("worker %d: has sent results to master in message of %d bytes\n", cfg->my_rank, pos));

	  FreeResults(results);
	  free(dsq);
	}
      if(status == eslEOD)ESL_DPRINTF1(("worker %d: has seen message to stop with this CM.\n", cfg->my_rank));
      else goto ERROR;
      FreeCM(cm);
      cm = NULL;
    }
  if (status == eslEOD) ESL_DPRINTF1(("Worker %d told CMs are done.\n", cfg->my_rank));
  else goto ERROR;
  
  if (wbuf != NULL) free(wbuf);
  return;

 ERROR:
  ESL_DPRINTF1(("worker %d: fails, is sending an error message, as follows:\n%s\n", cfg->my_rank, errbuf));
  pos = 0;
  MPI_Pack(&status, 1,                MPI_INT,  wbuf, wn, &pos, MPI_COMM_WORLD);
  MPI_Pack(errbuf,  cmERRBUFSIZE,    MPI_CHAR, wbuf, wn, &pos, MPI_COMM_WORLD);
  MPI_Send(wbuf, pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD);
  return;
}
#endif /*HAVE_MPI*/

/* A search work unit consists of a CM, digitized sequence dsq, and indices i and j.
 * The job is to search dsq from i..j and return search results in <*ret_results>.
 */
static int
process_search_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, 
			ESL_DSQ *dsq, int L, float min_cm_cutoff, float min_cp9_cutoff, search_results_t **ret_results)
{
  int status;
  search_results_t **results;
  int n;

  if(cm->fi == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "cm->fi is NULL in process_search_workunit()\n");

  ESL_ALLOC(results, sizeof(search_results_t *) * (cm->fi->nrounds+1));
  for(n = 0; n <= cm->fi->nrounds; n++) results[n] = CreateResults(INIT_RESULTS);

  if((status = ActuallySearchTarget(cm, errbuf, 0, dsq, 1, L, results, NULL, NULL)) != eslOK) goto ERROR;

  /* we only care about the final results, that survived all the rounds (all the filtering rounds plus the final round) */
  *ret_results = results[cm->fi->nrounds];
  /* free the results describing what survived each round of filtering (if any) */
  for(n = 0; n < cm->fi->nrounds; n++) FreeResults(results[n]);
  free(results);

  return eslOK;
  
 ERROR:
  ESL_DPRINTF1(("worker %d: has caught an error in process_search_workunit\n", cfg->my_rank));
  FreeCM(cm);
  return status;
}

/* A CP9 filter work unit consists of a CM and an int (nseq).
 * The job is to emit nseq sequences with a score better than cm->cutoff (rejecting
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
    if((! use_hmmonly) && (! esl_opt_GetBoolean(go, "--enfnohmm")))
      /* We want to filter with special enforced CP9 HMM for the enforced subseq */
      cm->search_opts |= CM_SEARCH_HMMFILTER;
    cm->enf_start = EnforceFindEnfStart(cm, esl_opt_GetInteger(go, "--enfstart"));
    /* --enfstart MUST have been enabled, --enfseq requires it */
    cm->enf_seq = esl_opt_GetString(go, "--enfseq");
  }

  /* search_opts */
  if(  esl_opt_GetBoolean(go, "-i"))            cm->search_opts |= CM_SEARCH_INSIDE;
  if(  esl_opt_GetBoolean(go, "--toponly"))     cm->search_opts |= CM_SEARCH_TOPONLY;
  if(  esl_opt_GetBoolean(go, "--noalign"))     cm->search_opts |= CM_SEARCH_NOALIGN;
  if(  esl_opt_GetBoolean(go, "--null2"))       cm->search_opts |= CM_SEARCH_NULL2;
  if(  esl_opt_GetBoolean(go, "--greedy"))      cm->search_opts |= CM_SEARCH_CMGREEDY;
  if(  esl_opt_GetBoolean(go, "--hmmgreedy"))   cm->search_opts |= CM_SEARCH_HMMGREEDY;
  if(  esl_opt_GetBoolean(go, "--noqdb"))       cm->search_opts |= CM_SEARCH_NOQDB;
  if(  esl_opt_GetBoolean(go, "--hmmfilter"))   cm->search_opts |= CM_SEARCH_HMMFILTER;
  if(  esl_opt_GetBoolean(go, "--hmmviterbi"))  cm->search_opts |= CM_SEARCH_HMMVITERBI;
  if(  esl_opt_GetBoolean(go, "--hmmforward"))  cm->search_opts |= CM_SEARCH_HMMFORWARD;
  if(  esl_opt_GetBoolean(go, "--hbanded"))     cm->search_opts |= CM_SEARCH_HBANDED;
  if(  esl_opt_GetBoolean(go, "--scan2bands"))  cm->search_opts |= CM_SEARCH_HMMSCANBANDS;
  if(  esl_opt_GetBoolean(go, "--sums"))        cm->search_opts |= CM_SEARCH_SUMS;

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

  /* finally, configure the CM for alignment based on cm->config_opts and cm->align_opts.
   * set local mode, make cp9 HMM, calculate QD bands etc. 
   */
  ConfigCM(cm, cfg->preset_dmin, cfg->preset_dmax); /* preset_d* usually NULL, unless --qdbfile */
  if(cm->config_opts & CM_CONFIG_ENFORCE) ConfigCMEnforce(cm);

  /* print qdbs to file if nec */
  if(! esl_opt_IsDefault(go, "--bfile")) {
    fprintf(cfg->bfp, "beta:%f\n", cm->beta);
    debug_print_bands(cfg->bfp, cm, cm->dmin, cm->dmax);
    fprintf(cfg->bfp, "beta:%f\n", cm->beta);
  }

  /* create FilterInfo, specifying no filtering, we change this later if necessary */
  cm_CreateFilterInfo(cm, 0.);
  if(cm->fi == NULL) cm_Fail("initialize_cm(), CreateFilterInfo() call failed.");
  cm_ValidateFilterInfo(cm->fi);
  cm_DumpFilterInfo(cm->fi);

  if(cfg->my_rank == 0) printf("CM %d: %s\n", (cfg->ncm), cm->name);
  return eslOK;
  
 ERROR:
  return status;
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
  for(i = 0; i < NGUMBELMODES; i++)
    for(p = 0; p < cm->stats->np; p++)
      {
	tmp_K = exp(cm->stats->gumAA[i][p]->mu * cm->stats->gumAA[i][p]->lambda) / 
	  cm->stats->gumAA[i][p]->L;
	cm->stats->gumAA[i][p]->mu = log(tmp_K * ((double) cfg->N)) /
	  cm->stats->gumAA[i][p]->lambda;
	cm->stats->gumAA[i][p]->L = cfg->N; /* update L, the seq size the stats correspond to */
      }
  printf ("CM/CP9 statistics read from CM file\n");
  if (cm->stats->np == 1) 
    printf ("No partition points\n");
  else 
    {
      printf ("Partition points are: ");
      for (p=0; p < cm->stats->np; p++)
	printf ("%d %d..%d", p, cm->stats->ps[p], cm->stats->pe[p]);
    }

  return eslOK;
}

/* set_cutoffs()
 * Determine cutoffs for the CM and HMM.
 */
static int
set_cutoffs(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, int *ret_cm_mode, int *ret_cp9_mode,
	    float *ret_min_cm_cutoff, float *ret_min_cp9_cutoff, int *ret_using_e_cutoff)
{
  int           cm_mode  = -1;   /* CM algorithm mode                        */
  int           cp9_mode = -1;   /* CP9 algorithm mode                       */
  int           use_hmmonly;     /* True if --hmmviterbi or --hmmforward enabled */

  /* Determine configuration of CM and CP9 based on cm->flags & cm->search_opts */
  CM2Gumbel_mode(cm, &cm_mode, &cp9_mode); 

  /* Set up CM cutoff, either 0 or 1 of 6 options is enabled. 
   * esl_opt_IsDefault() returns FALSE even if option is enabled with default value 
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
      cm->cutoff_type = E_CUTOFF;
      cm->cutoff      = esl_opt_GetReal(go, "-E");
    }
    else { /* no Gumbel stats in CM file, use default bit score cutoff */
      cm->cutoff_type = SCORE_CUTOFF;
      cm->cutoff      = esl_opt_GetReal(go, "-T");
    }
  }
  else if(! esl_opt_IsDefault(go, "-E")) {
    if(! (cm->flags & CMH_GUMBEL_STATS))
      ESL_FAIL(eslEINVAL, errbuf, "-E requires Gumbel statistics in <cm file>. Use cmcalibrate to get Gumbel stats.");
    cm->cutoff_type = E_CUTOFF;
    cm->cutoff      = esl_opt_GetReal(go, "-E");
  }
  else if(! esl_opt_IsDefault(go, "-T")) {
    cm->cutoff_type = SCORE_CUTOFF;
    cm->cutoff      = esl_opt_GetReal(go, "-T");
    if((cm->cutoff < 0.) && (! esl_opt_GetBoolean(go, "--greedy"))) ESL_FAIL(eslEINVAL, errbuf, "with -T <x> option, <x> can only be less than 0. if --greedy also enabled.");
  }
  else if(! esl_opt_IsDefault(go, "--ga")) {
    if(! (cm->flags & CMH_GA))
      ESL_FAIL(eslEINVAL, errbuf, "No GA gathering threshold in CM file, can't use --ga.");
    cm->cutoff_type = SCORE_CUTOFF;
    cm->cutoff      = esl_opt_GetReal(go, "--ga");
  }
  else if(! esl_opt_IsDefault(go, "--tc")) {
    if(! (cm->flags & CMH_TC))
      ESL_FAIL(eslEINVAL, errbuf, "No TC trusted cutoff in CM file, can't use --tc.");
    cm->cutoff_type = SCORE_CUTOFF;
    cm->cutoff      = esl_opt_GetReal(go, "--tc");
  }
  else if(! esl_opt_IsDefault(go, "--nc")) {
    if(! (cm->flags & CMH_NC))
      ESL_FAIL(eslEINVAL, errbuf, "No NC noise cutoff in CM file, can't use --nc.");
    cm->cutoff_type = SCORE_CUTOFF;
    cm->cutoff      = esl_opt_GetReal(go, "--nc");
  }
  else ESL_FAIL(eslEINCONCEIVABLE, errbuf, "No CM cutoff selected. This shouldn't happen.");

  /* Set up CP9 HMM cutoff, either 0 or 1 of 5 options is enabled 
   * esl_opt_IsDefault() returns FALSE even if option is enabled with default value 
   */
  if(esl_opt_IsDefault(go, "--hmmthr") && 
     esl_opt_IsDefault(go, "--hmmcalcthr") && 
     esl_opt_IsDefault(go, "--hmmE") && 
     esl_opt_IsDefault(go, "--hmmT"))
    {
      /* Choose from, in order of priority:
       * 1. filter threshold in CM file (if ! --hmmviterbi && ! --hmmforward) 
       * 2. default CP9 E value if CM file has Gumbel stats
       * 3. default CP9 bit score
       */
      use_hmmonly = (esl_opt_GetBoolean(go, "--hmmviterbi") || esl_opt_GetBoolean(go, "--hmmforward")) ? TRUE : FALSE;
      if((! use_hmmonly) && cm->flags & CMH_FTHR_STATS)  /* if !hmm_only use CP9 filter threshold from CM file */
	{
	  cm->cp9_cutoff_type = E_CUTOFF;
	  if     (cp9_mode == CP9_L) cm->cp9_cutoff = cm->stats->fthrA[cm_mode]->l_eval;
	  else if(cp9_mode == CP9_G) cm->cp9_cutoff = cm->stats->fthrA[cm_mode]->g_eval;
	  /* correct for new db size */
	  cm->cp9_cutoff *= (double) cfg->N / (double) cm->stats->fthrA[cm_mode]->db_size; 
	}
      else if(cm->flags & CMH_GUMBEL_STATS) /* use default CP9 E-value cutoff */
	{
	  cm->cp9_cutoff_type = E_CUTOFF;
	  cm->cp9_cutoff      = esl_opt_GetReal(go, "--hmmE"); 
	}
      else /* no filter stats nor Gumbel stats in CM file, use default bit score cutoff */
	{
	  cm->cp9_cutoff_type = SCORE_CUTOFF;
	  cm->cp9_cutoff      = esl_opt_GetReal(go, "--hmmT");
	}
    }
  else if(! esl_opt_IsDefault(go, "--hmmthr")) {
    if(! (cm->flags & CMH_FTHR_STATS))
      ESL_FAIL(eslEINVAL, errbuf, "--hmmthr requires filter threshold statistics in <cm file>. Use cmcalibrate to get CP9 filter threshold stats.");
    cm->cp9_cutoff_type = E_CUTOFF;
    if     (cp9_mode == CP9_L) cm->cp9_cutoff = cm->stats->fthrA[cm_mode]->l_eval;
    else if(cp9_mode == CP9_G) cm->cp9_cutoff = cm->stats->fthrA[cm_mode]->g_eval;
    /* correct for new db size */
    cm->cp9_cutoff *= (double) cfg->N / (double) cm->stats->fthrA[cm_mode]->db_size; 
  }
  else if(! esl_opt_IsDefault(go, "--hmmcalcthr")) {
    if(! (cm->flags & CMH_GUMBEL_STATS))
      ESL_FAIL(eslEINVAL, errbuf, "--hmmcalcthr requires Gumbel statistics in <cm file>. Use cmcalibrate to get Gumbel stats.");
    cm->cp9_cutoff_type = E_CUTOFF;
    /* this gets overwritten later after threshold is calculated */
    cm->cp9_cutoff = esl_opt_GetReal(go, "--hmmE");
  }
  else if(! esl_opt_IsDefault(go, "--hmmE")) {
    if(! (cm->flags & CMH_GUMBEL_STATS))
      ESL_FAIL(eslEINVAL, errbuf, "--hmmE requires Gumbel statistics in <cm file>. Use cmcalibrate to get Gumbel stats.");
    cm->cp9_cutoff_type = E_CUTOFF;
    cm->cp9_cutoff      = esl_opt_GetReal(go, "--hmmE");
  }
  else if(! esl_opt_IsDefault(go, "--hmmT")) {
    cm->cp9_cutoff_type = SCORE_CUTOFF;
    cm->cp9_cutoff      = esl_opt_GetReal(go, "--hmmT");
    if((cm->cp9_cutoff < 0.) && (! esl_opt_GetBoolean(go, "--hmmgreedy"))) ESL_FAIL(eslEINVAL, errbuf, "with -hmmT <x> option, <x> can only be less than 0. if --hmmgreedy also enabled.");
  }
  else ESL_FAIL(eslEINCONCEIVABLE, errbuf, "No CP9 cutoff selected. This shouldn't happen.");
     
  *ret_cm_mode  = cm_mode;
  *ret_cp9_mode = cp9_mode;

  *ret_min_cm_cutoff  = MinCMScCutoff(cm);
  *ret_min_cp9_cutoff = MinCP9ScCutoff(cm);

  if (((!use_hmmonly) && (cm->cutoff_type     == E_CUTOFF))  || 
      (( use_hmmonly) && (cm->cp9_cutoff_type == E_CUTOFF)))
    *ret_using_e_cutoff = TRUE;
  else *ret_using_e_cutoff = FALSE;

  /* finally, update fault FilterInfo, specifying no filtering, we'll change this later if nec */
  if(cm->fi == NULL) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "set_cutoffs(), cm->fi is NULL. Shouldn't happen");
  if(use_hmmonly) cm_UpdateFilterInfoCutoff(cm, cm->fi->nrounds, *ret_min_cp9_cutoff); /* update CP9 bit score cutoff */
  else            cm_UpdateFilterInfoCutoff(cm, cm->fi->nrounds, *ret_min_cm_cutoff);  /* use CM  bit score cutoff */
  cm_ValidateFilterInfo(cm->fi);
  cm_DumpFilterInfo(cm->fi);

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

  /* Setup ScanMatrix for CYK/Inside scanning functions, this is dependent on W, so we can't 
   * do it in initialize_cm(), but we have to wait til W is set, which just happened. 
   */
  int do_float = TRUE;
  int do_int   = FALSE;
  if(cm->search_opts & CM_SEARCH_INSIDE) { do_float = FALSE; do_int = TRUE; }
  cm_CreateScanMatrixForCM(cm, do_float, do_int);
  if(cm->smx == NULL) cm_Fail("initialize_cm(), CreateScanMatrixForCM() call failed.");


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
 * Function: print_search_info
 * Date:     EPN, Thu May 17 14:47:36 2007
 * Purpose:  Print info about search (cutoffs, algorithm, etc.) to file or stdout 
 */
int print_search_info(FILE *fp, CM_t *cm, int cm_mode, int cp9_mode, long N, char *errbuf)
{
  int p;
  float surv_fract;
  float avg_hit_len;
  int use_hmmonly;

  /* Could use ESL_GETOPTS here, but using the CM flags assures we're reporting
   * on how the CM is actually config'ed, not how we want it to be, also
   * not using ESL_GETOPTS makes this function portable (which isn't
   * really impt).
   */

  use_hmmonly = ((cm->search_opts & CM_SEARCH_HMMVITERBI) || (cm->search_opts & CM_SEARCH_HMMFORWARD));
  if(! use_hmmonly) {
    if(cm->cutoff_type == E_CUTOFF) {
      fprintf(fp, "CM cutoff (E value):  %.2f\n", cm->cutoff);
      for(p = 0; p < cm->stats->np; p++)
	fprintf(fp, "   GC %2d-%3d bit sc:  %.2f mu: %.5f lambda: %.5f\n", cm->stats->ps[p], cm->stats->pe[p], 
		(cm->stats->gumAA[cm_mode][p]->mu - (log(cm->cutoff) / cm->stats->gumAA[cm_mode][p]->lambda)), 
		cm->stats->gumAA[cm_mode][p]->mu, cm->stats->gumAA[cm_mode][p]->lambda);
    }		       
    else if (cm->cutoff_type == SCORE_CUTOFF) fprintf(fp, "CM cutoff (bit sc):   %.2f\n", cm->cutoff);
      printf ("CM search algorithm:  ");
      if(cm->search_opts & CM_SEARCH_INSIDE) fprintf(fp, "Inside\n");
      else                                   fprintf(fp, "CYK\n");
      printf ("CM configuration:     ");
      if(cm->flags & CMH_LOCAL_BEGIN) fprintf(fp, "Local\n");
      else                            fprintf(fp, "Glocal\n");
    }
  else 
    fprintf(fp, "Scanning with CP9 HMM only\n");
  if (cm->search_opts & CM_SEARCH_HMMFILTER)
    fprintf(fp, "Filtering with a CP9 HMM\n");
  
  if(use_hmmonly || cm->search_opts & CM_SEARCH_HMMFILTER) {
    if(cm->cp9_cutoff_type == E_CUTOFF) {
      if(!(cm->flags & CMH_GUMBEL_STATS))
	ESL_FAIL(eslEINCONCEIVABLE, errbuf, "trying to use E-value for CM cutoff, but CM has no Gumbel stats.");
      
      /* Predict survival fraction from filter based on E-value, consensus length, W and N */
      if(cp9_mode == CP9_G) avg_hit_len = cm->clen;       /* should be weighted sum of gamma[0] from QDB calc */
      if(cp9_mode == CP9_L) avg_hit_len = cm->clen * 0.5; /* should be weighted sum of gamma[0] from QDB calc */
      surv_fract = (cm->cp9_cutoff * ((2. * cm->W) - avg_hit_len)) / ((double) N); 
      /* HMM filtering sends j-W..i+W to be re-searched with CM for HMM hits i..j */
      fprintf(fp, "CP9 cutoff (E value): %.2f\n", cm->cp9_cutoff);
      fprintf(fp, "   Predicted survival fraction: %.5f (1/%.3f)\n", surv_fract, (1./surv_fract));
      for(p = 0; p < cm->stats->np; p++)
	fprintf(fp, "   GC %2d-%3d bit sc:  %.2f mu: %.5f lambda: %.5f\n", cm->stats->ps[p], cm->stats->pe[p], 
		(cm->stats->gumAA[cp9_mode][p]->mu - 
		 (log(cm->cp9_cutoff) / cm->stats->gumAA[cp9_mode][p]->lambda)), 
		cm->stats->gumAA[cp9_mode][p]->mu, cm->stats->gumAA[cp9_mode][p]->lambda);
    }
    else if (cm->cp9_cutoff_type == SCORE_CUTOFF) 
      fprintf(fp, "CP9 cutoff (bit sc):  %.2f\n", cm->cp9_cutoff);
    if     (cm->search_opts & CM_SEARCH_HMMVITERBI) printf ("CP9 search algorithm: Viterbi\n");
    else if(cm->search_opts & CM_SEARCH_HMMFORWARD) printf ("CP9 search algorithm: Forward\n");
    printf ("CP9 configuration:    ");
      if(cm->cp9->flags & CPLAN9_LOCAL_BEGIN) { 
	if(cm->cp9->flags & CPLAN9_EL) fprintf(fp, "Local (EL on)\n");
	else                           fprintf(fp, "Local (EL off)\n");
      }
      else                                    fprintf(fp, "Glocal\n");
    }
  printf     ("DB size, nt (N):      %ld\n\n", N);
  fflush(stdout);
  return eslOK;
}

#if HAVE_MPI
/* determine_cm_min_max_chunksize()
 * Given a CM, return the minimum and maximum subseq length to 
 * send to each process (min_chunksize & max_chunksize) based on 
 * minimum and maximum time we want a worker to spend working on a single
 * job, these are MPI_WORKER_MIN_SEC and MPI_WORKER_MAX_SEC constants from 
 * structs.h. 
 */
static int
determine_cm_min_max_chunksize(struct cfg_s *cfg, CM_t *cm, int *ret_min_chunksize, int *ret_max_chunksize)
{
  int   L = 1000000;
  float dpc; /* number of million dp calcs for L = 1 MB with or without QDB (depending on cm->search_opts) */
  /* TO DO: update CountScanDPCalcs() for HMM filtering, not just QDB */
  dpc = CountScanDPCalcs(cm, L, (! (cm->search_opts & CM_SEARCH_NOQDB))) / 1000000;
  printf("dpc: %f\n", dpc);
  *ret_min_chunksize = (int) (((float) MPI_WORKER_MIN_SEC * (float) L) / (dpc / MDPC_SEC));
  *ret_max_chunksize = (int) (((float) MPI_WORKER_MAX_SEC * (float) L) / (dpc / MDPC_SEC));
  return eslOK;
}

/* determine_seq_chunksize()
 *
 * Given a sequence length, return the appropriate subseq length
 * to send to each process (chunksize). Ideally we send exactly 1 chunk
 * to each processor, but only if the time required is within a range
 * given by MPI_WORKER_MIN_SEC and MPI_WORKER_MAX_SEC (in structs.h). 
 * The chunksizes corresponding to those min/max times have been
 * precomputed for the current CM, and are sent in as min_chunksize
 * and max_chunksize.
 */
static int
determine_seq_chunksize(struct cfg_s *cfg, int L, int min_chunksize, int max_chunksize)
{
  int chunksize;
  if(cfg->do_rc) L *= 2;
  chunksize = (int) ((float) L / (float) (cfg->nproc - 1));
  if(chunksize < min_chunksize) return min_chunksize;
  if(chunksize > max_chunksize) return max_chunksize;
  return chunksize;
}
#endif

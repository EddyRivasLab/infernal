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
#include <string.h>
#include <ctype.h>
#include <time.h>

#include "easel.h"              /* better general sequence analysis library */
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_msa.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"

#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* external functions                   */
#include "hmmband.h"
#include "stats.h"
#include "mpifuncs.h"
#include "cm_dispatch.h"

#define STRATOPTS  "--cmonly,--hmmfilter,--hmmonly"                 /* Exclusive choice for search strategy */
#define CMCUTOPTS  "-E,-T,--ga,--tc,--nc,--negsc"                   /* Exclusive choice for CM cutoff */
#define HMMCUTOPTS "--hmmthr,--hmmcalcthr,--hmmE,--hmmT,--hmmnegsc" /* Exclusive choice for HMM cutoff */
#define ALPHOPTS   "--rna,--dna"                                    /* Exclusive choice for output alphabet */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs       incomp  help  docgroup*/
  /* basic options */
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",   1 },
  { "-g",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--hmmonly", "configure CM for glocal alignment [default: local]", 1 },
  { "-i",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--hmmonly", "use scanning Inside algorithm instead of CYK", 1 },
  { "--informat",eslARG_STRING, NULL,  NULL, NULL,      NULL,      NULL,        NULL, "specify the input file is in format <x>, not FASTA", 1 },
  { "--toponly", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "only search the top strand", 1 },
  { "--noalign", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "find start/stop/score only; don't do alignments", 1 },
  { "--window",  eslARG_INT,    NULL,  NULL, "n>0",     NULL,      NULL,        NULL, "set scanning window size to <n> [default: calculated]", 1 },
  { "--null2",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--hmmonly", "turn on the post hoc second null model", 1 },
  { "--iins",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "allow informative insert emissions, do not zero them", 1 },
  { "--elsilent",eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "disallow CM local end (EL) emissions", 1 },
  { "--rtrans",  eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--hmmonly", "replace CM transition scores from <cmfile> with RSEARCH scores", 1 },
  { "--greedy",  eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--hmmonly", "resolve overlapping hits with a greedy algorithm a la RSEARCH", 1 },
  /* strategy choice */
  { "--cmonly",   eslARG_NONE,"default",NULL,NULL,      STRATOPTS, NULL,        NULL, "search only with CM, do not filter [default]", 2 },
  { "--hmmfilter",eslARG_NONE,  FALSE, NULL, NULL,      STRATOPTS, NULL,        NULL, "subseqs j-W+1..i+W-1 survive (j=end from Fwd, i=start from Bwd)", 2 },
  { "--hmmonly",  eslARG_NONE,  FALSE, NULL, NULL,      STRATOPTS, NULL,        NULL, "do not use CM at all, just scan with HMM (Forward +  Backward)", 2 },
  /* CM cutoff options */
  { "-E",        eslARG_REAL,   "0.1", NULL, "x>0.",    CMCUTOPTS, NULL, "--hmmonly", "use cutoff E-value of <x> for CM search", 3 },
  { "-T",        eslARG_REAL,   "0.0", NULL, "x>=0.",   CMCUTOPTS, NULL, "--hmmonly", "use cutoff bit score of <x> for CM search", 3 },
  { "--ga",      eslARG_NONE,   FALSE, NULL, NULL,      CMCUTOPTS, NULL, "--hmmonly", "use CM \"GA\" gathering threshold as cutoff bit score", 3 },
  { "--tc",      eslARG_NONE,   FALSE, NULL, NULL,      CMCUTOPTS, NULL, "--hmmonly", "use CM \"TC\" trusted cutoff as cutoff bit score", 3 },
  { "--nc",      eslARG_NONE,   FALSE, NULL, NULL,      CMCUTOPTS, NULL, "--hmmonly", "use CM \"NC\" noise cutoff as cutoff bit score", 3 },
  { "--negsc",   eslARG_REAL,   FALSE, NULL, NULL,      CMCUTOPTS, NULL, "--hmmonly", "set minimum CM bit score to report as <x> < 0 (experimental!)", 3 },
  /* HMM cutoff options */
  { "--hmmthr",  eslARG_NONE,   FALSE, NULL, NULL,     HMMCUTOPTS,"--hmmfilter",NULL, "use HMM filter from cmcalibrate (in <cm file>)", 4 },
  { "--hmmcalcthr",eslARG_NONE, FALSE, NULL, NULL,     HMMCUTOPTS,"--hmmfilter",NULL, "calculate HMM filter threshold by sampling from CM", 4 },
  { "--hmmE",    eslARG_REAL,   "1.0", NULL, "x>0.",   HMMCUTOPTS, NULL,  "--cmonly", "use cutoff E-value of <x> for CP9 HMM filter/search", 4 },
  { "--hmmT",    eslARG_REAL,   "0.0", NULL, "x>=0.",  HMMCUTOPTS, NULL,  "--cmonly", "use cutoff bit score of <x> for CP9 HMM filter/search", 4 },
  { "--hmmnegsc",eslARG_REAL,   FALSE, NULL, NULL,     HMMCUTOPTS, NULL,  "--cmonly", "set minimum HMM bit score to report as <x> < 0 (experimental!)", 4 },
   /* QDB related options */
  { "--beta",    eslARG_REAL,   "1e-7",NULL, "x>0",     NULL,      NULL, "--hmmonly", "set tail loss prob for QDB and window size calculation to <x>", 5 },
  { "--noqdb",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL, "--hmmonly", "DO NOT use query dependent banding (QDB) for acceleration", 5 },
  { "--qdbfile", eslARG_STRING, NULL,  NULL, NULL,      NULL,      NULL,"--hmmonly,--noqdb","read QDBs from file <s> (outputted from cmbuild)", 5 },
  /* HMM filtering options */
  { "--hmmpad",  eslARG_INT,    NULL,  NULL, NULL,      NULL,"--hmmfilter",     NULL, "subseqs \'i-<n>..j+<n>\' survive", 6 },
  { "--hbanded", eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--hmmfilter",     NULL, "calculate and use HMM bands in CM search", 6 },
  { "--tau",     eslARG_REAL,   "1E-7",NULL, "0<x<1",   NULL,"--hbanded",       NULL, "set tail loss prob for --hbanded to <x>", 6 },
  { "--scan2bands",eslARG_NONE, FALSE, NULL, NULL,      NULL,"--hbanded",       NULL, "derive HMM bands from scanning Forward/Backward", 6 },
  { "--sums",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--hbanded",       NULL, "use posterior sums during HMM band calculation (widens bands)", 6 },
  /* HMM configuration options */
  { "--hmmlocal",eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,  "--cmonly", "configure HMM for local alignment [default: glocal]", 7 },
  { "--hmmnoel", eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,  "--cmonly", "DO NOT enable HMM EL local ends that mirror CM", 7 },
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
  /* Enforcing a subsequence */
  { "--enfstart",eslARG_INT,    FALSE, NULL, "n>0",     NULL,"--enfseq",        NULL, "enforce MATL stretch starting at consensus position <n>", 9 },
  { "--enfseq",  eslARG_STRING, NULL,  NULL, NULL,      NULL,"--enfstart",      NULL, "enforce MATL stretch starting at --enfstart <n> emits seq <s>", 9 },
  { "--enfnohmm",eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--enfstart",      NULL, "DO NOT filter first w/an HMM that only enforces --enfseq <s>", 9 },
  /* verbose output files */
  { "--tfile",   eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "dump parsetrees for each hit to file <f>", 10 },
  { "--gcfile",  eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "save GC content stats of target sequence file to <f>", 10 },
  { "--bfile",   eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "save bands for each state to file <f>", 10 },
  { "--filhfile",eslARG_OUTFILE, NULL, NULL, NULL,      NULL,"--hmmcalcthr",    NULL, "save CP9 filter threshold histogram(s) to file <s>", 10 },
  { "--filrfile",eslARG_OUTFILE, NULL, NULL, NULL,      NULL,"--hmmcalcthr",    NULL, "save CP9 filter threshold information file <s>", 10 },
/* Setting output alphabet */
  { "--rna",     eslARG_NONE,"default",NULL, NULL,  ALPHOPTS,      NULL,        NULL, "output alignment as RNA sequence data", 11 },
  { "--dna",     eslARG_NONE,   FALSE, NULL, NULL,  ALPHOPTS,      NULL,        NULL, "output alignment as DNA (not RNA) sequence data", 11 },
/* Other options */
  { "--stall",   eslARG_NONE,  FALSE, NULL, NULL,       NULL,      NULL,    NULL, "arrest after start: for debugging MPI under gdb",   12 },  
#ifdef HAVE_MPI
  { "--mpi",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,  "--qdbfile","run as an MPI parallel program", 12 },  
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
  int           ncm;            /* number CM we're at in file */

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
static int read_qdb_file(FILE *fp, CM_t *cm, int *dmin, int *dmax);
static int initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int is_integer(char *s);

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go = NULL;   /* command line processing                     */
  ESL_STOPWATCH   *w  = esl_stopwatch_Create();
  struct cfg_s     cfg;

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
      puts("\noptions for enforcing a single-stranded subsequence:");
      esl_opt_DisplayHelp(stdout, go, 9, 2, 80);
      puts("\nverbose output files:");
      esl_opt_DisplayHelp(stdout, go, 10, 2, 80);
      puts("\noptions for selecting output alphabet:");
      esl_opt_DisplayHelp(stdout, go, 11, 2, 80);
      puts("\nother options:");
      esl_opt_DisplayHelp(stdout, go, 12, 2, 80);
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
  else    esl_fatal("Can't determine output alphabet");
  cfg.ncm        = 0;
  cfg.cmfp       = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.tfp        = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.bfp        = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.filhfp     = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.filrfp     = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.preset_dmin= NULL;                   /* filled in initialize_cm() only if --qdbfile, which conflicts with --mpi */
  cfg.preset_dmax= NULL;                   /* filled in initialize_cm() only if --qdbfile, which conflicts with --mpi */

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
      cfg.be_verbose = FALSE;
      MPI_Init(&argc, &argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &(cfg.my_rank));
      MPI_Comm_size(MPI_COMM_WORLD, &(cfg.nproc));

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
 *    cfg->abc         - digital input alphabet
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

  /* open CM file */
  if ((cfg->cmfp = CMFileOpen(cfg->cmfile, NULL)) == NULL)
    ESL_FAIL(eslFAIL, NULL, "Failed to open covariance model save file %s\n", cfg->cmfile);

  /* optionally, open trace file */
  if (esl_opt_GetString(go, "--tfile") != NULL) {
    if ((cfg->tfp = fopen(esl_opt_GetString(go, "--tfile"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --tfile output file %s\n", esl_opt_GetString(go, "--tfile"));
    }

  /* optionally, open bands file */
  if (esl_opt_GetString(go, "--bfile") != NULL) {
    if ((cfg->tfp = fopen(esl_opt_GetString(go, "--bfile"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --bfile output file %s\n", esl_opt_GetString(go, "--bfile"));
    }

  /* optionally, open filter threshold calc histogram file */
  if (esl_opt_GetString(go, "--filhfile") != NULL) {
    if ((cfg->tfp = fopen(esl_opt_GetString(go, "--filhfile"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --filhfile output file %s\n", esl_opt_GetString(go, "--filhfile"));
    }

  /* optionally, open filter threshold calc info file */
  if (esl_opt_GetString(go, "--filrfile") != NULL) {
    if ((cfg->tfp = fopen(esl_opt_GetString(go, "--filrfile"), "w")) == NULL) 
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
  int      status;
  char     errbuf[eslERRBUFSIZE];
  CM_t    *cm = NULL;
  
  if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);
  /* init_shared_cfg UNNEC? */
  /*if ((status = init_shared_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);*/

  while (CMFileRead(cfg->cmfp, &(cfg->abc), &cm))
    {
      if (cm == NULL) cm_Fail("Failed to read CM from %s -- file corrupt?\n", cfg->cmfile);
      cfg->ncm++;

      /* initialize the flags/options/params of the CM*/
      if((status = initialize_cm(go, cfg, errbuf, cm) != eslOK)) esl_fatal("ERROR status: %d in initialize_cm()", status);
	 /* Configure the CM for alignment based on cm->config_opts and cm->align_opts.
	  * set local mode, make cp9 HMM, calculate QD bands etc. */
      ConfigCM(cm, cfg->preset_dmin, cfg->preset_dmax); /* preset_d* usually NULL, unless --qdbfile */
      if(cm->config_opts & CM_CONFIG_ENFORCE) ConfigCMEnforce(cm);
      FreeCM(cm);
    }	 
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

  /* set up CM parameters that are option-changeable */
  cm->beta   = esl_opt_GetReal(go, "--beta"); /* this will be DEFAULT_BETA unless changed at command line */
  cm->tau    = esl_opt_GetReal(go, "--tau");  /* this will be DEFAULT_TAU unless changed at command line */
  if(! esl_opt_IsDefault(go, "--negsc"))    cm->sc_boost     = -1. * esl_opt_GetReal(go, "--negsc");
  if(! esl_opt_IsDefault(go, "--hmmnegsc")) cm->cp9_sc_boost = -1. * esl_opt_GetReal(go, "--hmmnegsc");
  if(! esl_opt_IsDefault(go, "--hmmpad"))   cm->hmmpad       =       esl_opt_GetInteger(go, "--hmmpad");

  /* Update cm->config_opts and cm->align_opts based on command line options */

  /* config_opts */
  if(! esl_opt_GetBoolean(go, "-g"))            cm->config_opts |= CM_CONFIG_LOCAL;
  if(  esl_opt_GetBoolean(go, "--hmmlocal"))    cm->config_opts |= CM_CONFIG_HMMLOCAL;
  if(! esl_opt_GetBoolean(go, "--hmmnoel"))     cm->config_opts |= CM_CONFIG_HMMEL;
  if(! esl_opt_GetBoolean(go, "--iins"))        cm->config_opts |= CM_CONFIG_ZEROINSERTS;
  if(  esl_opt_GetBoolean(go, "--elsilent"))    cm->config_opts |= CM_CONFIG_ELSILENT;
  /* set up QDB? Yes, unless --window and --noqdb both enabled */
  if(! ((! esl_opt_IsDefault(go, "--window")) && (esl_opt_GetBoolean(go, "--noqdb"))))
    cm->config_opts |= CM_CONFIG_QDB;
  /* are we enforcing a subseq? */
  if(! esl_opt_IsDefault (go, "--enfseq")) {
    cm->config_opts |= CM_CONFIG_ENFORCE;
    if((! esl_opt_GetBoolean(go, "--hmmonly")) && (! esl_opt_GetBoolean(go, "--enfnohmm")))
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
  if(  esl_opt_GetBoolean(go, "--hmmonly"))     cm->search_opts |= CM_SEARCH_HMMONLY;
  if(! esl_opt_IsDefault (go, "--hmmpad"))      cm->search_opts |= CM_SEARCH_HMMPAD;
  if(  esl_opt_GetBoolean(go, "--hbanded"))     cm->search_opts |= CM_SEARCH_HBANDED;
  if(  esl_opt_GetBoolean(go, "--sums"))        cm->search_opts |= CM_SEARCH_SUMS;

  /* If do_enforce set do_hmm_rescan to TRUE if we're filtering or scanning with an HMM,
   * this way only subseqs that include the enf_subseq should pass the filter */
  if((esl_opt_GetBoolean(go, "--hmmrescan")) ||    
      (((! esl_opt_IsDefault(go, "--enfseq")) && (! esl_opt_GetBoolean(go, "--enfnohmm"))) ||
       ((! esl_opt_IsDefault(go, "--enfseq")) && (esl_opt_GetBoolean(go, "--hmmfilter")))))
    cm->search_opts |= CM_SEARCH_HMMRESCAN; 
  
  /* flags */
  if(  esl_opt_GetBoolean(go, "--rtrans"))      cm->flags       |= CM_RSEARCHTRANS;

  /* finally, read in QDBs if nec */
  if(! esl_opt_IsDefault(go, "--qdbfile"))
    {
      /* can't be in MPI mode, --mpi is incompatible with --qdbfile */
      FILE *qdb_fp;
      if(cfg->preset_dmin != NULL || cfg->preset_dmax != NULL)
	esl_fatal("ERROR: trying to read QDBs for more than one CM. With --qdbfile, <cm file> must have exactly 1 CM in it.");
      ESL_ALLOC(cfg->preset_dmin, sizeof(int) * cm->M);
      ESL_ALLOC(cfg->preset_dmax, sizeof(int) * cm->M);
      if ((qdb_fp = fopen(esl_opt_GetString(go, "--qdbfile"), "r")) == NULL)
	esl_fatal("failed to open QDB file %s", esl_opt_GetString(go, "--qdbfile"));
      if(!(read_qdb_file(qdb_fp, cm, cfg->preset_dmin, cfg->preset_dmax)))
	esl_fatal("ERROR reading QDB file: %s.\nDoes it correspond (same number of states) to this model?\n", esl_opt_GetString(go, "--qdbfile"));
      fclose(qdb_fp);
    }

  printf("CM %d: %s\n", (cfg->ncm), cm->name);

  /* HEREHEREHEREHERE*/

  return eslOK;

 ERROR: 
  return status;
}

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


/* cmalign.c
 * SRE, Thu Jul 25 11:28:03 2002 [St. Louis]
 * SVN $Id$
 * 
 * Align sequences to a CM.
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
#include <float.h>

#include "easel.h"		/* general seq analysis library   */
#include "esl_alphabet.h"
#include "esl_getopts.h"		
#include "esl_msa.h"
#include "esl_sqio.h"		
#include "esl_stopwatch.h"

#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* external functions                   */
#include "hmmband.h"         
#include "cplan9.h"
#include "cm_postprob.h"
#include "mpifuncs.h"
#include "cm_dispatch.h"

#define ALGOPTS  "--cyk,--inside,--outside,--hmmonly"        /* Exclusive choice for scoring algorithms */
#define MEMOPTS  "--small,--nosmall"                         /* Exclusive choice for memory choice */
#define ACCOPTS  "--nonbanded,--hbanded,--qdb"               /* Exclusive choice for acceleration strategies */
#define ALPHOPTS "--rna,--dna"                               /* Exclusive choice for output alphabet */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs       incomp  help  docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",   1 },
  { "-l",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "align locally w.r.t. the model",         1 },
  { "-o",        eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "output the alignment to file <f>, not stdout", 1 },
  { "-q",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "quiet; suppress verbose banner",         1 },
  { "-1",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "use tabular output summary format, 1 line per sequence", 1 },
#ifdef HAVE_MPI
  { "--mpi",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "run as an MPI parallel program",                    1 },  
#endif
  /* Miscellaneous expert options */
  { "--informat",eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL,        NULL, "specify input alignment is in format <s>, not Stockholm",  2 },
  { "--time",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "print timings for alignment, band calculation, etc.", 2 },
  { "--regress", eslARG_OUTFILE,FALSE, NULL, NULL,      NULL,      NULL,        NULL, "save regresion test data to file <f>", 2 },
  { "--full",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "include all match columns in output alignment", 2 },
  { "--zeroinserts",eslARG_NONE,FALSE, NULL, NULL,      NULL,      NULL,        NULL, "set insert emission scores to 0", 2 },
  { "--elsilent",eslARG_NONE,   FALSE, NULL, NULL,      NULL,      "-l",        NULL, "disallow local end (EL) emissions", 2 },
  /* Algorithm options */
  { "--cyk",     eslARG_NONE,"default",NULL, NULL,   ALGOPTS,      NULL,        NULL, "don't align; return scores from the Inside algorithm", 3 },
  { "--hmmonly", eslARG_NONE,   FALSE, NULL, NULL,   ALGOPTS,      NULL,"--hbanded,--qdb", "align using Viterbi to CM plan 9 HMM only",3 },
  { "--inside",  eslARG_NONE,   FALSE, NULL, NULL,   ALGOPTS,   "--cyk",        NULL, "don't align; return scores from the Inside algorithm", 3 },
  { "--outside", eslARG_NONE,   FALSE, NULL, NULL,   ALGOPTS,      NULL,        NULL, "don't align; return scores from the Outside algorithm", 3 },
  { "--post",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "align with CYK and append posterior probabilities", 3 },
  { "--checkpost",eslARG_NONE,  FALSE, NULL, NULL,      NULL,  "--post",        NULL, "check that posteriors are correctly calc'ed", 3 },
  { "--sub",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "build sub CM for columns b/t HMM predicted start/end points", 3 },
  /* Memory options */
  { "--small",   eslARG_NONE,"default",  NULL, NULL,      NULL,      NULL, "--nosmall", "use divide and conquer (d&c) alignment algorithm", 4 },
  { "--nosmall", eslARG_NONE,   FALSE,  NULL, NULL,      NULL,      NULL,   "--small", "use normal alignment algorithm, not d&c", 4 },
  /* Banded alignment */
  { "--nonbanded",eslARG_NONE,"default",NULL, NULL,   ACCOPTS,      NULL,      NULL,   "accelerate using CM plan 9 HMM banded CYK aln algorithm", 5 },
  { "--hbanded", eslARG_NONE,   FALSE, NULL, NULL,"--nonbanded,--qdb,--small",NULL,NULL, "accelerate using CM plan 9 HMM banded CYK aln algorithm", 5 },
  { "--tau",     eslARG_REAL,   "1E-7",NULL, "0<x<1",   NULL,"--hbanded",       NULL, "set tail loss prob for --hbanded to <x>", 5 },
  { "--hsafe",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--hbanded",       NULL, "realign (non-banded) seqs with HMM banded CYK score < 0 bits", 5 },
  { "--sums",    eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--hbanded",       NULL, "use posterior sums during HMM band calculation (widens bands)", 5 },
  { "--qdb",     eslARG_NONE,   FALSE, NULL, NULL,   ACCOPTS,      NULL,        NULL, "use query dependent banded CYK alignment algorithm", 5 },
  { "--beta",    eslARG_REAL,   "1E-7",NULL, "0<x<1",   NULL,   "--qdb",        NULL, "set tail loss prob for --qdb to <x>", 5 },
  /* Including a preset alignment */
  { "--withali", eslARG_STRING, NULL,  NULL, NULL,      NULL,    "--cyk",        NULL, "enforce MATL stretch starting at --enfstart <n> emits seq <s>", 6 },
  { "--rf",      eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--withali",       NULL, "--rf was originally used with cmbuild", 6 },
  { "--gapthresh",eslARG_REAL,  "0.5", NULL, "0<=x<=1", NULL,"--withali",       NULL, "--gapthresh <x> was originally used with cmbuild", 6 },
  /* Verbose output files/debugging */
  { "--tfile",   eslARG_OUTFILE,FALSE, NULL, NULL,      NULL,      NULL,        NULL, "dump individual sequence tracebacks to file <f>", 7 },
  { "--banddump",eslARG_NONE,   FALSE, NULL, "1<=n<=3", NULL,      NULL,        NULL, "set verbosity of band info print statements to <n> [1..3]", 7 },
  { "--dlev",    eslARG_NONE,   FALSE, NULL, "1<=n<=3", NULL,      NULL,        NULL, "set verbosity of debugging print statements to <n> [1..3]", 7 },
  /* Enforcing a subsequence */
  { "--enfstart",eslARG_INT,    FALSE, NULL, "n>0",     NULL,"--enfseq",        NULL, "enforce MATL stretch starting at consensus position <n>", 8 },
  { "--enfseq",  eslARG_STRING, NULL,  NULL, NULL,      NULL,"--enfstart",      NULL, "enforce MATL stretch starting at --enfstart <n> emits seq <s>", 8 },
/* Setting output alphabet */
  { "--rna",     eslARG_NONE,"default",NULL, NULL,  ALPHOPTS,      NULL,        NULL, "output alignment as RNA sequence data [default]", 9},
  { "--dna",     eslARG_NONE,   FALSE, NULL, NULL,  ALPHOPTS,      NULL,        NULL, "output alignment as DNA (not RNA) sequence data", 9},
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
  int           be_verbose;	/* one-line-per-seq summary */
  int           nseq;           /* which number sequence this is in file (only valid in serial mode) */
  CM_t         *cm;             /* the CM to align to */

  int           do_mpi;		/* TRUE if we're doing MPI parallelization */
  int           nproc;		/* how many MPI processes, total */
  int           my_rank;	/* who am I, in 0..nproc-1 */
  int           do_stall;	/* TRUE to stall the program until gdb attaches */

  /* Masters only (i/o streams) */
  CMFILE       *cmfp;		/* open input CM file stream       */
  FILE         *ofp;		/* output file (default is stdout) */
  FILE         *tracefp;	/* optional output for parsetrees  */
  FILE         *regressfp;	/* optional output for regression test  */
  ESL_MSAFILE  *withalifp;	/* optional input alignment to include */
  ESL_MSA      *withmsa;	/* MSA from withalifp to include */
  ESL_ALPHABET *abc_out; 	/* digital alphabet for writing */
};

static char usage[]  = "[-options] <cmfile> <sequence file>";
static char banner[] = "align sequences to an RNA CM";

static int  init_master_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int  init_shared_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);

static void  serial_master (const ESL_GETOPTS *go, struct cfg_s *cfg);
#ifdef HAVE_MPI
static void  mpi_master    (const ESL_GETOPTS *go, struct cfg_s *cfg);
static void  mpi_worker    (const ESL_GETOPTS *go, struct cfg_s *cfg);
#endif

static int    process_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t **ret_cm, Parsetree_t ***opt_tr);
static int output_result(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, ESL_SQ **sq, Parsetree_t **tr, CP9trace_t **cp9_tr, char **postcode);

static int initialize_cm(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf);
static int check_withali(const ESL_GETOPTS *go, const struct cfg_s *cfg, ESL_MSA **ret_msa);
static int include_withali(const ESL_GETOPTS *go, struct cfg_s *cfg, ESL_SQ ***ret_sq, Parsetree_t ***ret_tr, char *errbuf);
static int compare_cms(CM_t *cm1, CM_t *cm2);
static int make_aligned_string(char *aseq, char *gapstring, int alen, char *ss, char **ret_s);

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
      puts("\nexpert miscellaneous options:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\noptions specifying alignment algorithm:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); 
      puts("\nmemory related options:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80);
      puts("\nbanded dynamic programming acceleration options:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80);
      puts("\noptions for including a fixed alignment within output alignment:");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80);
      puts("\nverbose output files and debugging:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80);
      puts("\noptions for enforcing alignment of a single-stranded subsequence:");
      esl_opt_DisplayHelp(stdout, go, 8, 2, 80);
      puts("\noptions for selecting output alphabet::");
      esl_opt_DisplayHelp(stdout, go, 9, 2, 80);
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
  if (esl_opt_GetBoolean(go, "-1")) cfg.be_verbose = FALSE;        
  else                              cfg.be_verbose = TRUE;        
  if      (esl_opt_GetBoolean(go, "--rna")) cfg.abc_out = esl_alphabet_Create(eslRNA);
  else if (esl_opt_GetBoolean(go, "--dna")) cfg.abc_out = esl_alphabet_Create(eslDNA);
  else    esl_fail("Can't determine output alphabet");
  cfg.cmfp       = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.ofp        = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.tracefp    = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.regressfp  = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.withalifp  = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.withmsa    = NULL;	           /* filled in init_master_cfg() in masters, stays NULL for workers */
  cfg.nseq       = 0;		           /* this counter is incremented in masters */

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
  if (cfg.my_rank == 0) esl_stopwatch_Display(cfg.ofp, w, "# CPU time: ");

  /* Clean up the shared cfg. 
   */
  if (cfg.my_rank == 0) {
    if (! esl_opt_IsDefault(go, "-o")) { fclose(cfg.ofp); }
    if (cfg.sqfp      != NULL) esl_sqfile_Close(cfg.sqfp);
    if (cfg.abc       != NULL) esl_alphabet_Destroy(cfg.abc);
    if (cfg.tracefp   != NULL) fclose(cfg.tracefp);
    if (cfg.regressfp != NULL) fclose(cfg.regressfp);
    if (cfg.withalifp != NULL) esl_msafile_Close(cfg.withalifp);
    if (cfg.withmsa   != NULL) esl_msa_Destroy(cfg.withmsa);
  }
  esl_getopts_Destroy(go);
  esl_stopwatch_Destroy(w);
  return 0;
}

/* init_master_cfg()
 * Called by masters, mpi or serial.
 * Already set:
 *    cfg->cmfile  - command line arg 1
 *    cfg->sqfile  - command line arg 2
 *    cfg->fmt     - format of output file
 * Sets: 
 *    cfg->sqfp      - open sequence file                
 *    cfg->ofp       - output file (stdout by default)
 *    cfg->cmfp      - open CM file                
 *    cfg->abc       - digital input alphabet
 *    cfg->tracefp   - optional output file
 *    cfg->regressfp - optional output file
 *    cfg->withalifp - optional input alignment file to include
 *    cfg->withmsa   - MSA from --withali file 
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
  CMFILE *cmfp; 
  char *cmfile;

  cmfile = esl_opt_GetArg(go, 1);

  /* open input sequence file */
  status = esl_sqfile_Open(cfg->sqfile, cfg->fmt, NULL, &(cfg->sqfp));
  if (status == eslENOTFOUND)    ESL_FAIL(status, errbuf, "File %s doesn't exist or is not readable\n", cfg->sqfile);
  else if (status == eslEFORMAT) ESL_FAIL(status, errbuf, "Couldn't determine format of sequence file %s\n", cfg->sqfile);
  else if (status == eslEINVAL)  ESL_FAIL(status, errbuf, "Can’t autodetect stdin or .gz."); 
  else if (status != eslOK)      ESL_FAIL(status, errbuf, "Sequence file open failed with error %d\n", status);
  cfg->fmt = cfg->sqfp->format;

  /* open output file */
  if (esl_opt_GetString(go, "-o") != NULL) {
    if ((cfg->ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open -o output file %s\n", esl_opt_GetString(go, "-o"));
    } else cfg->ofp = stdout;

  /* open CM file and read (only) the first CM in it. */
  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL)
    ESL_FAIL(eslFAIL, NULL, "Failed to open covariance model save file %s\n", cmfile);
  if(!CMFileRead(cmfp, NULL, &(cfg->cm))) cm_Fail(errbuf);
  ESL_FAIL(eslFAIL, NULL, "Failed to read a CM from %s -- file corrupt?\n", cmfile);
  if (cfg->cm == NULL) ESL_FAIL(eslFAIL, NULL, "Failed to read a CM from %s -- file empty?\n", cmfile);
  CMFileClose(cmfp);

  /* optionally, open trace file */
  if (esl_opt_GetString(go, "--tfile") != NULL) {
    if ((cfg->tracefp = fopen(esl_opt_GetString(go, "--tfile"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --tfile output file %s\n", esl_opt_GetString(go, "--tfile"));
    }

  /* optionally, open regression file */
  if (esl_opt_GetString(go, "--regress") != NULL) {
    if ((cfg->regressfp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) 
	ESL_FAIL(eslFAIL, errbuf, "Failed to open --regress output file %s\n", esl_opt_GetString(go, "--regress"));
    }

  /* optionally, open withali file for reading */
  if(esl_opt_GetString(go, "--withali") != NULL)
    {
      status = esl_msafile_Open(esl_opt_GetString(go, "--withali"), eslMSAFILE_UNKNOWN, NULL, &(cfg->withalifp));
      if (status == eslENOTFOUND)    ESL_FAIL(status, errbuf, "--withali alignment file %s doesn't exist or is not readable\n", 
					      esl_opt_GetString(go, "--withali"));
      else if (status == eslEFORMAT) ESL_FAIL(status, errbuf, "Couldn't determine format of --withali alignment %s\n", 
					      esl_opt_GetString(go, "--withali"));
      else if (status != eslOK)      ESL_FAIL(status, errbuf, "Alignment file open failed with error %d\n", status);

      if((status = check_withali(go, cfg, &(cfg->withmsa))) != eslOK)
	ESL_FAIL(status, errbuf, "--withali alignment file %s doesn't have a SS_cons compatible with the CM\n", status);
    }
  return eslOK;
}

/* serial_master()
 * The serial version of cmbuild.
 * For each MSA, build a CM and save it.
 * 
 * A master can only return if it's successful. All errors are handled immediately and fatally with cm_Fail().
 */
static void
serial_master(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int      status;
  char     errbuf[eslERRBUFSIZE];
  ESL_SQ **sq;
  CM_t    *cm;
  int       i;
  Parsetree_t    **tr;          /* parse trees for the sequences            */
  CP9trace_t     **cp9_tr;      /* CP9 traces for the sequences             */
  char           **postcode;    /* posterior decode array of strings        */

  
  if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);
  /* init_shared_cfg UNNEC? */
  /*if ((status = init_shared_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);*/

  /* initialize the flags/options/params of the CM*/
  initialize_cm(go, cfg, errbuf);
  /* Configure the CM for alignment based on cm->config_opts and cm->align_opts.
   * set local mode, make cp9 HMM, calculate QD bands etc. */
  ConfigCM(cfg->cm, NULL, NULL);
  if(cm->config_opts & CM_CONFIG_ENFORCE) ConfigCMEnforce(cm);

  /* Align the sequences, there's no 'process_workunit() function' but rather different functions for serial/mpi */ 
  if(esl_opt_GetBoolean(go, "--hmmonly"))
    {
      serial_align_targets(cfg->sqfp, cm, &sq, NULL,  &postcode, &cp9_tr, &cfg->nseq, 
			   esl_opt_GetInteger(go, "--banddump"), esl_opt_GetInteger(go, "--dlev"), esl_opt_GetBoolean(go, "-q"));
      if ((status = output_result(go, cfg, errbuf, sq, tr,     NULL, postcode))        != eslOK) cm_Fail(errbuf);
    }
  else
    {
      serial_align_targets(cfg->sqfp, cm, &sq, &tr,   &postcode, NULL,    &cfg->nseq,
			   esl_opt_GetInteger(go, "--banddump"), esl_opt_GetInteger(go, "--dlev"), esl_opt_GetBoolean(go, "-q"));
      if ((status = output_result(go, cfg, errbuf, sq, NULL, cp9_tr, postcode))        != eslOK) cm_Fail(errbuf);
    }

  /* clean up */
  if(esl_opt_GetBoolean(go, "--cyk")) {
    for (i = 0; i < cfg->nseq; i++) FreeParsetree(tr[i]);
    free(tr);
  }
  if(esl_opt_GetBoolean(go, "--hmmonly")) {
    for (i = 0; i < cfg->nseq; i++) CP9FreeTrace(cp9_tr[i]);
    free(cp9_tr);
  }
  for (i = 0; i < cfg->nseq; i++) esl_sq_Destroy(sq[i]);
  free(sq);
  return;

 ERROR:
  cm_Fail("Reallocation error.\n");
}

static int
output_result(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, ESL_SQ **sq, Parsetree_t **tr, CP9trace_t **cp9_tr, char **postcode)
{
  int status;
  ESL_MSA *msa = NULL;
  int i, ip;

  /* Open the output file set up ofp
  /* Output the tabular results header. 
   */
  if (esl_opt_GetBoolean(go, "-1"))
    {
      printf("# %3s %-20s %5s %5s %5s\n", "idx", "name",                 "bitsc",  "len",  "time?");
      printf("#%4s %-20s %5s %5s %5s\n", "----", "--------------------", "-----", "-----", "-----");
    }

  /* create a new MSA, if we didn't do either --inside or --outside */
  if(esl_opt_GetBoolean(go, "--cyk") || esl_opt_GetBoolean(go, "--hmmonly"))
  {
    /* optionally include a fixed alignment provided with --withali,
     * this has already been checked to see it matches the CM structure */
    if(esl_opt_GetString(go, "--withali") != NULL)
      {
	if((status = include_withali(go, cfg, &sq, &tr, errbuf)) != eslOK)
	  ESL_FAIL(status, errbuf, "--withali alignment file %s doesn't have a SS_cons compatible with the CM\n", status);
      }

    if(esl_opt_GetBoolean(go, "--hmmonly"))
      if((status = CP9Traces2Alignment(cfg->cm, sq, NULL, cfg->nseq, cp9_tr, esl_opt_GetBoolean(go, "--full"), &msa)) != eslOK)
	goto ERROR;
    else
      if((status = Parsetrees2Alignment(cfg->cm, sq, NULL, tr, cfg->nseq, esl_opt_GetBoolean(go, "--full"), &msa)) != eslOK)
	goto ERROR;
    
    if(esl_opt_GetBoolean(go, "--post")) 
      {                                                                              
	char *apostcode;   /* aligned posterior decode array */
	if(postcode == NULL) 
	  cm_Fail("ERROR --post enabled, but {serial,parallel}_align_targets() did not return post codes.\n");
	for (i = cfg->withmsa->nseq; i < cfg->nseq; i++)                                                   
	  {                                                                          
	    ip = i - cfg->withmsa->nseq;
	    if((status =make_aligned_string(msa->aseq[i], "-_.", msa->alen, postcode[ip], &apostcode)) != eslOK)
	      ESL_FAIL(status, errbuf, "error creating posterior string\n", status);
	    esl_msa_AppendGR(msa, "POST", i, apostcode);                                  
	    free(apostcode);                                                         
	    free(postcode[ip]);                                                       
	  }                                                                          
	free(postcode);                                                              
      }                                                                              
    status = esl_msa_Write(cfg->ofp, msa, eslMSAFILE_STOCKHOLM);
    if      (status == eslEMEM) ESL_FAIL(status, errbuf, "Memory error when outputting alignment\n");
    else if (status != eslOK)   ESL_FAIL(status, errbuf, "Writing alignment file failed with error %d\n", status);

    /* Detailed traces for debugging training set. */
    if(cfg->tracefp != NULL)
      {
	if(esl_opt_GetBoolean(go,"--hmmonly")) { printf("%-40s ... ", "Saving CP9 HMM traces"); fflush(stdout); }
	else                                   { printf("%-40s ... ", "Saving CM parsetrees");  fflush(stdout); }
	for (i = 0; i < msa->nseq; i++) 
	  {
	    fprintf(cfg->tracefp, "> %s\n", sq[i]->name);
	    if(esl_opt_GetBoolean(go,"--hmmonly")) 
	      {
		fprintf(cfg->tracefp, "  SCORE : %.2f bits\n", CP9TraceScore(cfg->cm->cp9, sq[i]->dsq, cp9_tr[i]));
		CP9PrintTrace(cfg->tracefp, cp9_tr[i], cfg->cm->cp9, sq[i]->dsq);
	      }
	    else
	      {
		fprintf(cfg->tracefp, "  SCORE : %.2f bits\n", ParsetreeScore(cfg->cm, tr[i], sq[i]->dsq, FALSE));
		ParsetreeDump(cfg->tracefp, tr[i], cfg->cm, sq[i]->dsq);
	      }
	    fprintf(cfg->tracefp, "//\n");
	  }
	printf("done. [%s]\n", esl_opt_GetString(go, "--tfile"));
      }
    if (cfg->regressfp != NULL)
      {
	/* Must delete author info from msa, because it contains version
	 * and won't diff clean in regression tests. */
	if(msa->au != NULL) free(msa->au); msa->au = NULL;
	status = esl_msa_Write(cfg->ofp, msa, eslMSAFILE_STOCKHOLM);
	if (status == eslEMEM)    ESL_FAIL(status, errbuf, "Memory error when outputting regression file\n");
	else if (status != eslOK) ESL_FAIL(status, errbuf, "Writing regression file failed with error %d\n", status);
      }
    }
  return eslOK;

 ERROR:
  if(msa != NULL) esl_msa_Destroy(msa);
  return status;
}

/* initialize_cm()
 * Setup the CM based on the command-line options/defaults;
 * only set flags and a few parameters. ConfigCM() configures
 * the CM.
 */
static int
initialize_cm(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf)
{
  if (cfg->be_verbose) {
    fprintf(cfg->ofp, "%-40s ... ", "Initializing CM ");  
    fflush(cfg->ofp); 
  }

  /* set up params/flags/options of the CM */
  cfg->cm->beta   = esl_opt_GetReal(go, "--beta"); /* this will be DEFAULT_BETA unless changed at command line */
  cfg->cm->tau    = esl_opt_GetReal(go, "--tau");  /* this will be DEFAULT_TAU unless changed at command line */

  /* Update cfg->cm->config_opts and cfg->cm->align_opts based on command line options */
  if(esl_opt_GetBoolean(go, "-l"))
    {
      cfg->cm->config_opts |= CM_CONFIG_LOCAL;
      cfg->cm->config_opts |= CM_CONFIG_HMMLOCAL;
      cfg->cm->config_opts |= CM_CONFIG_HMMEL;
    }
  if(esl_opt_GetBoolean(go, "--elsilent"))    cfg->cm->config_opts |= CM_CONFIG_ELSILENT;
  if(esl_opt_GetBoolean(go, "--zeroinserts")) cfg->cm->config_opts |= CM_CONFIG_ZEROINSERTS;
  if(esl_opt_GetBoolean(go, "--hbanded"))     cfg->cm->align_opts  |= CM_ALIGN_HBANDED;
  if(esl_opt_GetBoolean(go, "--sums"))        cfg->cm->align_opts  |= CM_ALIGN_SUMS;
  if(esl_opt_GetBoolean(go, "--sub"))         cfg->cm->align_opts  |= CM_ALIGN_SUB;
  if(esl_opt_GetBoolean(go, "--hmmonly"))     cfg->cm->align_opts  |= CM_ALIGN_HMMONLY;
  if(esl_opt_GetBoolean(go, "--inside"))      cfg->cm->align_opts  |= CM_ALIGN_INSIDE;
  if(esl_opt_GetBoolean(go, "--outside"))     cfg->cm->align_opts  |= CM_ALIGN_OUTSIDE;
  if(esl_opt_GetBoolean(go, "--nosmall"))     cfg->cm->align_opts  |= CM_ALIGN_NOSMALL;
  if(esl_opt_GetBoolean(go, "--post"))        cfg->cm->align_opts  |= CM_ALIGN_POST;
  if(esl_opt_GetBoolean(go, "--time"))        cfg->cm->align_opts  |= CM_ALIGN_TIME;
  if(esl_opt_GetBoolean(go, "--checkpost"))   cfg->cm->align_opts  |= CM_ALIGN_CHECKINOUT;
  if(esl_opt_GetBoolean(go, "--hsafe"))       cfg->cm->align_opts  |= CM_ALIGN_HMMSAFE;
  if(esl_opt_GetBoolean(go, "--enfseq"))
    {
      cfg->cm->config_opts |= CM_CONFIG_ENFORCE;
      cfg->cm->enf_start    = EnforceFindEnfStart(cfg->cm, esl_opt_GetInteger(go, "--enfstart"));
      cfg->cm->enf_seq      = esl_opt_GetString(go, "--enfseq");
    }
  if(esl_opt_GetBoolean(go, "--qdb"))          
    { 
      cfg->cm->align_opts  |= CM_ALIGN_QDB;
      cfg->cm->config_opts |= CM_CONFIG_QDB;
    }
  if (cfg->be_verbose) fprintf(cfg->ofp, "done.\n");
  return eslOK;
}

/* Function: compare_cms()
 * EPN, Tue Mar  6 08:32:12 2007
 *
 * Purpose:  Given two CMs, cm1 and cm2, compare them, returning TRUE 
 *           iff they have the same guide tree (same node architecture).
 *
 * Args:     cm1          - covariance model number 1
 *           cm2          - covariance model number 2
 * 
 * Returns:  TRUE if CMs have same guide tree, FALSE otherwise
 */
static int compare_cms(CM_t *cm1, CM_t *cm2)
{
  int          nd; 
  if(cm1->nodes != cm2->nodes) return FALSE;
  for(nd = 0; nd < cm1->nodes; nd++)
    if(cm1->ndtype[nd] != cm2->ndtype[nd]) return FALSE;
  return TRUE;
}

/* Function: check_withali()
 * EPN, Tue Mar  6 06:25:02 2007
 *
 * Purpose:  Ensure that the alignment to include has a secondary
 *           structure that matches our CM. Pass the alignment back
 *           as *ret_msa.
 *
 * Returns:  <eslOK> on success.
 *           <eslEINCOMPAT> if alignment doesn't match the CM 
 */
static int check_withali(const ESL_GETOPTS *go, const struct cfg_s *cfg, ESL_MSA **ret_msa)
{
  int           status;
  ESL_MSA      *msa      = NULL; /* alignment we're including  */
  CM_t         *new_cm   = NULL; /* CM built from MSA, we check it has same guide tree as 'cm' */
  CM_t         *newer_cm = NULL; /* used briefly if we rebalance new_cm */
  Parsetree_t  *mtr      = NULL; /* master structure tree from the alignment*/
  char          errbuf[eslERRBUFSIZE];
  CM_BG *bg              = NULL;

  /* cfg->withalifp is open */
  status = esl_msa_Read(cfg->withalifp, &msa);
  if (status == eslEFORMAT)  cm_Fail("--withali alignment file parse error, line %d of file %s:\n%s\nOffending line is:\n%s\n", cfg->withalifp->linenumber, cfg->withalifp->fname, cfg->withalifp->errbuf, cfg->withalifp->buf);
  else if (status != eslOK)       cm_Fail("--withali alignment file read unexpectedly failed with code %d\n", status);

  /* Some input data cleaning. */
  if (esl_opt_GetBoolean(go, "--rf") && msa->rf == NULL) 
    ESL_FAIL(eslFAIL, errbuf, "--rf invoked but --withali alignment has no reference coord annotation.\n");
  if (msa->ss_cons == NULL) 
    ESL_FAIL(eslFAIL, errbuf, "--withali alignment did not contain consensus structure annotation.\n");
  if (! clean_cs(msa->ss_cons, msa->alen))
    ESL_FAIL(eslFAIL, errbuf, "Failed to parse consensus structure annotation for --withali alignment\n");

  /* Build a CM from a master guide tree built from the msa, 
   * then check to make sure this CM has same emit map as the CM
   * we've had passed in. This is fragile and hopefully temporary. 
   * Another solution would be to use a checksum, but CM files don't 
   * have checksums yet.
   */
  bg = cm_bg_Create(msa->abc); /* default values, A,C,G,U = 0.25  */
  HandModelmaker(msa, bg, esl_opt_GetBoolean(go, "--rf"), esl_opt_GetReal(go, "--gapthresh"), &new_cm, &mtr);
  if(!(compare_cm_guide_trees(cfg->cm, new_cm)))
    {
      CM_t *newer_cm;
      newer_cm = CMRebalance(new_cm);
      FreeCM(new_cm);
      new_cm = NULL;
      if(!(compare_cm_guide_trees(cfg->cm, newer_cm)))
	{
	  status = eslEINCOMPAT;
	  goto ERROR;
	}
    }

  /* if we get here, the CM guide trees match */
  if(new_cm   != NULL) FreeCM(new_cm);
  if(newer_cm != NULL) FreeCM(newer_cm);
  if(bg       != NULL) cm_bg_Destroy(bg);
  FreeParsetree(mtr);
  *ret_msa = msa;
  return eslOK;

 ERROR:
  if(msa != NULL)      esl_msa_Destroy(msa);
  if(new_cm   != NULL) FreeCM(new_cm);
  if(newer_cm != NULL) FreeCM(newer_cm);
  if(mtr != NULL)      FreeParsetree(mtr);
  if(bg  != NULL)      cm_bg_Destroy(bg);
  return eslEINCOMPAT;
}

/* Function: include_withali()
 * EPN, Tue Mar  6 06:25:02 2007
 *
 * Purpose:  Infer the implicit parses of sequences in an
 *           MSA to a CM and append to passed in data structures
 *           We've already checked to make sure the MSA's consensus
 *           structure matches the CM.
 *
 * Args:     go           - command line options
 *           cfg          - cmalign configuration, includes msa to add
 *           ret_sq       - pre-existing sequences, to append msa seqs to
 *           ret_tr       - pre-existing parsetrees, to append to
 *           errbuf       - easel error message
 * 
 * Returns:  eslOK on success, eslEMEM on memory error
 *           Also new, realloc'ed arrays for sq, tr in ret_seq, ret_tr; 
 *           cfg->nseq is increased by number of seqs in cfg->withmsa.
*/
static int include_withali(const ESL_GETOPTS *go, struct cfg_s *cfg, ESL_SQ ***ret_sq, Parsetree_t ***ret_tr, char *errbuf)
{
  int           status;
  void         *tmp;      /* for ESL_RALLOC() */
  int           i;	  /* counter over aseqs       */
  int           ip;	  /* offset counter over aseqs */
  char        **uaseq;    /* unaligned seqs, dealigned from the MSA */
  int           apos;     /*   aligned position index */
  int           uapos;    /* unaligned position index */
  int           x;        /* counter of parsetree nodes */
  int         **map;      /* [0..msa->nseq-1][0..msa->alen] map from aligned
			   * positions to unaligned (non-gap) positions */
  char *aseq;                   
  Parsetree_t *mtr;

  /* Contract check */
  if(cfg->withmsa == NULL) esl_fatal("ERROR in include_withali() withmsa is NULL.\n");

  /* For each seq in the MSA, map the aligned sequences coords to 
   * the unaligned coords, we stay in digitized seq coords (1..alen),
   * we need this for converting parsetrees from Transmogrify, which
   * have emitl and emitr in aligned coords to unaligned coords, so 
   * we can call Parsetrees2Alignment() with them. */
  ESL_ALLOC(map,   sizeof(int *)  * cfg->withmsa->nseq);
  ESL_ALLOC(uaseq, sizeof(char *) * cfg->withmsa->nseq);
  for (i = 0; i < cfg->withmsa->nseq; i++)
    {
      ESL_ALLOC(map[i],   sizeof(int)  * (cfg->withmsa->alen+1));
      ESL_ALLOC(uaseq[i], sizeof(char) * (cfg->withmsa->alen+1));
      map[i][0] = -1; /* invalid */
      uapos = 1;
      for(apos = 0; apos < cfg->withmsa->alen; apos++)
	{
	  if (!esl_abc_CIsGap(cfg->withmsa->abc, cfg->withmsa->aseq[i][(apos+1)]))
	    map[i][(apos+1)] = uapos++;
	  else
	    map[i][(apos+1)] = -1;
	}
      esl_strdup(cfg->withmsa->aseq[i], -1, &(uaseq[i]));
      esl_sq_Dealign(uaseq[i], uaseq[i], "-_.", cfg->withmsa->alen);
    }
  ESL_RALLOC((*ret_tr), tmp, (sizeof(Parsetree_t *) * (cfg->nseq + cfg->withmsa->nseq)));
  ESL_RALLOC((*ret_sq), tmp, (sizeof(ESL_SQ *)      * (cfg->nseq + cfg->withmsa->nseq)));

  /* Swap some pointers so the included alignment appears at the top of the output 
   * alignment instead of the bottom. */
  for(i = 0; i < cfg->nseq; i++)
    {
      ip = i + cfg->withmsa->nseq;
      (*ret_tr)[ip] = (*ret_tr)[i];
      (*ret_sq)[ip] = (*ret_sq)[i];
    }

  /* Transmogrify each aligned seq to get a parsetree */
  for (i = 0; i < cfg->withmsa->nseq; i++)
    {
      esl_abc_Textize(cfg->withmsa->abc, cfg->withmsa->ax[i], cfg->withmsa->alen, aseq);
      (*ret_tr)[i] = Transmogrify(cfg->cm, mtr, cfg->withmsa->ax[i], cfg->withmsa->aseq[i], cfg->withmsa->alen);
      free(aseq);
      /* ret_tr[i] is in alignment coords, convert it to unaligned coords, */
      for(x = 0; x < (*ret_tr)[i]->n; x++)
	{
	  if((*ret_tr)[i]->emitl[x] != -1)
	    (*ret_tr)[i]->emitl[x] = map[i][(*ret_tr)[i]->emitl[x]];
	  if((*ret_tr)[i]->emitr[x] != -1)
	    (*ret_tr)[i]->emitr[x] = map[i][(*ret_tr)[i]->emitr[x]];
	}
      (*ret_sq)[i]      = esl_sq_CreateFrom(cfg->withmsa->sqname[i], uaseq[i], NULL, NULL, NULL);
      esl_sq_Digitize(cfg->cm->abc, (*ret_sq)[i]);
    }
  cfg->nseq += cfg->withmsa->nseq; /* we added these seqs to sq, tr */
  /* Clean up and exit. */
  esl_Free2D((void **) map,   cfg->withmsa->nseq);
  esl_Free2D((void **) uaseq, cfg->withmsa->nseq);
  return eslOK;

 ERROR:
  esl_Free2D((void **) map,   cfg->withmsa->nseq);
  esl_Free2D((void **) uaseq, cfg->withmsa->nseq);
  return status;
}


/* Function: make_aligned_string() 
 * Incept:   EPN, Thu Aug  2 18:24:49 2007
 *           stolen from Squid:alignio.c:MakeAlignedString()
 * 
 * Purpose:  Given a raw string of some type (secondary structure, say),
 *           align it to a given aseq by putting gaps wherever the
 *           aseq has gaps. 
 *           
 * Args:     aseq:  template for alignment
 *           gapstring: defines all gap chars ex: "-_."
 *           alen:  length of aseq
 *           ss:    raw string to align to aseq
 *           ret_s: RETURN: aligned ss
 *           
 * Return:   eslOK on success, eslEMEM on memory allocation error,
 *           eslEINCONCEIVABLE on strange error,
 *           ret_ss is alloc'ed here and must be free'd by caller.
 */
int
make_aligned_string(char *aseq, char *gapstring, int alen, char *ss, char **ret_s)
{
  int status;
  char *new; 
  int   apos, rpos;

  ESL_ALLOC(new, (sizeof(char) * (alen+1)));
  for (apos = rpos = 0; apos < alen; apos++)
    {
      if (strchr(gapstring, aseq[apos]) != NULL)
	new[apos] = aseq[apos];
      else
	new[apos] = ss[rpos++];
    }
  new[apos] = '\0';

  if (rpos != strlen(ss))
    {
      if(new != NULL) free(new);
      return eslEINCONCEIVABLE;
    }
  *ret_s = new;
  return eslOK;

 ERROR:
  if(new != NULL) free(new);
  return status;
}

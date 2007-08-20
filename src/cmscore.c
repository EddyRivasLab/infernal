/* SRE, Thu Aug  3 17:08:45 2000 [StL]
 * SVN $Id$
 * 
 * Score a CM against unaligned sequence examples.
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
#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"		
#include "esl_msa.h"
#include "esl_sqio.h"		
#include "esl_stopwatch.h"

#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* external functions                   */
#include "cm_dispatch.h"	/* alignment functions                   */

#define ALGOPTS  "--std,--qdb,--qdbsmall,--qdbboth,--hbanded,--hmmonly"  /* Exclusive choice for scoring algorithms */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs       incomp  help  docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",   1 },
  { "-l",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,     "--sub", "align locally w.r.t. the model",         1 },
  { "-i",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "print individual timings & scores, not just summary", 1 },
  { "-1",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "use tabular output summary format, 1 line per sequence", 1 },
#ifdef HAVE_MPI
  { "--mpi",     eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "run as an MPI parallel program",                    1 },  
#endif
  /* Miscellaneous expert options */
  { "--stringent",eslARG_NONE,  FALSE, NULL, NULL,      NULL,      NULL,        NULL, "require the two parsetrees to be identical", 2 },
  { "--sub",      eslARG_NONE,  FALSE, NULL, NULL,      NULL,      NULL,        "-l", "build sub CM for columns b/t HMM predicted start/end points", 2 },
  { "--zeroinserts",eslARG_NONE,FALSE, NULL, NULL,      NULL,      NULL,        NULL, "set insert emission scores to 0", 2 },
  { "--elsilent",eslARG_NONE,   FALSE, NULL, NULL,      NULL,      "-l",        NULL, "disallow local end (EL) emissions", 2 },
  /* Output options */
  { "--regress", eslARG_OUTFILE, NULL, NULL, NULL,      NULL,      NULL,        NULL, "save regression test data to file <f>", 3 },
  { "--tfile",   eslARG_OUTFILE,NULL,  NULL, NULL,      NULL,      NULL,        NULL, "dump parsetrees to file <f>",  3 },
  /* Stage 2 algorithm options */
  { "--std",     eslARG_NONE,"default",NULL, NULL,   ALGOPTS,      NULL,        NULL, "compare divide and conquer versus standard CYK", 4 },
  { "--qdb",     eslARG_NONE,   FALSE, NULL, NULL,   ALGOPTS,      NULL,        NULL, "compare non-banded d&c versus QDB standard CYK", 4 },
  { "--qdbsmall",eslARG_NONE,   FALSE, NULL, NULL,   ALGOPTS,      NULL,        NULL, "compare non-banded d&c versus QDB d&c", 4 },
  { "--qdbboth", eslARG_NONE,   FALSE, NULL, NULL,   ALGOPTS,      NULL,        NULL, "compare        QDB d&c versus QDB standard CYK", 4 },
  { "--beta",    eslARG_REAL,   "1E-7",NULL, "0<x<1",   NULL,      NULL,        NULL, "set tail loss prob for QDB to <x>", 4 },
  { "--hbanded", eslARG_NONE,   FALSE,  NULL, NULL,  ALGOPTS,      NULL,        NULL, "accelerate using CM plan 9 HMM banded CYK aln algorithm", 4 },
  { "--tau",     eslARG_REAL,   "1E-7",NULL, "0<x<1",   NULL,"--hbanded",       NULL, "set tail loss prob for --hbanded to <x>", 4 },
  { "--hsafe",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--hbanded",       NULL, "realign (non-banded) seqs with HMM banded CYK score < 0 bits", 4 },
  { "--hmmonly", eslARG_NONE,   FALSE, NULL, NULL,   ALGOPTS,      NULL,        NULL, "align to a CM Plan 9 HMM with the Viterbi algorithm", 4 },
  { "--scoreonly",eslARG_NONE,  FALSE, NULL, NULL,   ALGOPTS,      NULL,        NULL, "for standard CYK stage, do only score, save memory", 4 },
  /* Options for testing multiple rounds of banded alignment, stage 2->N alignment */
  { "--betas",   eslARG_INT,  NULL,    NULL, "0<n<50",   NULL, "--betae",       NULL, "set initial (stage 2) tail loss prob to 10E-<x> for QDB", 5 },
  { "--betae",   eslARG_INT,  NULL,    NULL, "0<n<50",   NULL, "--betas",       NULL, "set final   (stage N) tail loss prob to 10E-<x> for QDB", 5 },
  { "--taus",    eslARG_INT,  NULL,    NULL, "0<n<50",   NULL,"--hbanded,--taue",NULL,"set initial (stage 2) tail loss prob to 10E-<x> for HMM banding", 5 },
  { "--taue",    eslARG_INT,  NULL,    NULL, "0<n<50",   NULL,"--hbanded,--taus",NULL,"set final   (stage N) tail loss prob to 10E-<x> for HMM banding", 5 },
/* Other options */
  { "--stall",   eslARG_NONE,  FALSE, NULL, NULL,       NULL,      NULL,    NULL, "arrest after start: for debugging MPI under gdb",   10 },  
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
  int           fmt;		/* format code for seqfile */
  ESL_ALPHABET *abc;		/* digital alphabet for input */
  int           be_verbose;	/* one-line-per-seq summary */
  int           nseq;           /* which number sequence this is in file (only valid in serial mode) */
  CM_t         *cm;             /* the CM to align to */
  int           nstages;        /* number of stages of alignment we'll do */
  int           s;              /* which stage we're on [0..nstages-1], 0 = stage 1, 1 = stage 2 etc. */
  double       *beta;           /* [0..nstages-1] beta for each stage, NULL if not doing QDB */
  double       *tau;            /* [0..nstages-1] tau  for each stage, NULL if not doing HMM banding */
  float        *s1_sc;          /* [0..cfg->nseq] scores for stage 1 parses, filled in 1st output_result() call*/
  ESL_STOPWATCH *s1_w;          /* stopwatch for timing stage 1 */
  ESL_STOPWATCH *w;             /* stopwatch for timing stages 2+ */

  int           do_mpi;		/* TRUE if we're doing MPI parallelization */
  int           nproc;		/* how many MPI processes, total */
  int           my_rank;	/* who am I, in 0..nproc-1 */
  int           do_stall;	/* TRUE to stall the program until gdb attaches */

  /* Masters only (i/o streams) */
  CMFILE       *cmfp;		/* open input CM file stream       */
  FILE         *tracefp;	/* optional output for parsetrees  */
  FILE         *regressfp;	/* optional output for regression test  */
};

static char usage[]  = "[-options] <cmfile> <sequence file>";
static char banner[] = "score RNA covariance model against sequences";

static int  init_master_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
/* static int  init_shared_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf); */

static void  serial_master (const ESL_GETOPTS *go, struct cfg_s *cfg);
#ifdef HAVE_MPI
static void  mpi_master    (const ESL_GETOPTS *go, struct cfg_s *cfg);
static void  mpi_worker    (const ESL_GETOPTS *go, struct cfg_s *cfg);
#endif

/* static int process_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t **ret_cm, Parsetree_t ***opt_tr); */
static int output_result(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, ESL_SQ **sq, Parsetree_t **tr, CP9trace_t **cp9_tr,
			 float *scoreonly_sc);
static int initialize_cm(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf);
static int summarize_align_options(CM_t *cm);

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go = NULL;   /* command line processing                     */
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
      puts("\noutput options are:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); 
      puts("\nstage 2 alignment options, to compare to stage 1 (D&C non-banded), are:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 
      puts("\noptions for testing multiple tau/beta values for --hbanded/--qdb:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80); 
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
  /* Check for incompatible option combinations I don't know how to disallow with esl_getopts */
  /* Can't have --betas and --betae without a --qdb* option */
  if ((! esl_opt_IsDefault(go, "--betas")) && (! esl_opt_IsDefault(go, "--betae")))
    if(! ((esl_opt_GetBoolean(go, "--qdb")) || (esl_opt_GetBoolean(go, "--qdbsmall"))))
    {
	printf("Error parsing options, --betas and --betae combination requires either --qdb or --qdbsmall.\n");
	exit(1);
      }
  if ((! esl_opt_IsDefault(go, "--betas")) && (! esl_opt_IsDefault(go, "--betae")))
    if((esl_opt_GetInteger(go, "--betas")) > (esl_opt_GetInteger(go, "--betae")))
      {
	printf("Error parsing options, --betas <n> argument must be less than --betae <n> argument.\n");
	exit(1);
      }
  if ((! esl_opt_IsDefault(go, "--taus")) && (! esl_opt_IsDefault(go, "--taue")))
    if((esl_opt_GetInteger(go, "--taus")) > (esl_opt_GetInteger(go, "--taue")))
      {
	printf("Error parsing options, --taus <n> argument must be less than --taue <n> argument.\n");
	exit(1);
      }


  /* Initialize what we can in the config structure (without knowing the input alphabet yet).
   */
  cfg.cmfile     = esl_opt_GetArg(go, 1); 
  cfg.sqfile     = esl_opt_GetArg(go, 2); 
  cfg.fmt        = eslSQFILE_UNKNOWN;      /* autodetect sequence file format by default. */ 
  cfg.abc        = NULL;	           /* created in init_master_cfg() in masters, or in mpi_worker() in workers */
  if (esl_opt_GetBoolean(go, "-1")) cfg.be_verbose = FALSE;        
  else                              cfg.be_verbose = TRUE;        
  cfg.cmfp       = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.tracefp    = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.regressfp  = NULL;	           /* opened in init_master_cfg() in masters, stays NULL for workers */
  cfg.nseq       = 0;		           /* this counter is incremented in masters */
  cfg.cm         = NULL;                   /* created in init_master_cfg() in masters, or in mpi_worker() in workers */
  cfg.nstages    = 0;                      /* set in init_master_cfg() in masters, stays 0 for workers */
  cfg.s          = 0;                      /* initialized to 0 in init_master_cfg() in masters, stays 0 for workers */
  cfg.beta       = NULL;                   /* alloc'ed and filled in init_master_cfg() in masters, stays NULL in workers */
  cfg.tau        = NULL;                   /* alloc'ed and filled in init_master_cfg() in masters, stays NULL in workers */
  cfg.s1_sc      = NULL;                   /* alloc'ed and filled in first call to output_result() in masters, stays NULL in workers */
  cfg.s1_w       = NULL;                   /* created in init_master_cfg in masters, stays NULL in workers */
  cfg.w          = NULL;                   /* created in init_master_cfg in masters, stays NULL in workers */


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
    }

  /* Clean up the shared cfg. 
   */
  if (cfg.my_rank == 0) {
    if (cfg.tracefp   != NULL) {
      printf("Parsetrees saved in file %s.\n", esl_opt_GetString(go, "--tfile"));
      fclose(cfg.tracefp);
    }
    if (cfg.regressfp   != NULL) {
      printf("Regression data (parsetrees) saved in file %s.\n", esl_opt_GetString(go, "--regress"));
      fclose(cfg.regressfp);
    }
    if (cfg.s1_sc     != NULL) free(cfg.s1_sc);
  }
  if (cfg.cm        != NULL) FreeCM(cfg.cm);
  if (cfg.abc       != NULL) esl_alphabet_Destroy(cfg.abc);
  if (cfg.beta      != NULL) free(cfg.beta);
  if (cfg.tau       != NULL) free(cfg.tau);
  if (cfg.s1_w      != NULL) esl_stopwatch_Destroy(cfg.s1_w);
  if (cfg.w         != NULL) esl_stopwatch_Destroy(cfg.w);
  esl_getopts_Destroy(go);
  return 0;
}

/* init_master_cfg()
 * Called by masters, mpi or serial.
 * Already set:
 *    cfg->cmfile      - command line arg 1
 *    cfg->sqfile      - command line arg 2
 *    cfg->fmt         - format of output file
 * Sets: 
 *    cfg->cmfp        - open CM file                
 *    cfg->abc         - digital input alphabet
 *    cfg->tracefp     - optional output file
 *    cfg->regressfp   - optional output file
 *    cfg->nstages     - number of alignment stages
 *    cfg->beta        - beta values for each stage
 *    cfg->tau         - tau values for each stage
 *    cfg->s1_w        - stopwatch for timing stage 1
 *    cfg->w           - stopwatch for timing stage 2+
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

  /* open CM file and read (only) the first CM in it. */
  if ((cmfp = CMFileOpen(cfg->cmfile, NULL)) == NULL)
    cm_Fail("Failed to open covariance model save file %s\n", cfg->cmfile);
  if(!CMFileRead(cmfp, &(cfg->abc), &(cfg->cm)))
    cm_Fail("Failed to read a CM from %s -- file corrupt?\n", cfg->cmfile);
  if (cfg->cm == NULL) cm_Fail("Failed to read a CM from %s -- file empty?\n", cfg->cmfile);
  CMFileClose(cmfp);

  /* optionally, open trace file */
  if (esl_opt_GetString(go, "--tfile") != NULL) {
    if ((cfg->tracefp = fopen(esl_opt_GetString(go, "--tfile"), "w")) == NULL) 
	cm_Fail("Failed to open --tfile output file %s\n", esl_opt_GetString(go, "--tfile"));
    }

  /* optionally, open regression file */
  if (esl_opt_GetString(go, "--regress") != NULL) {
    if ((cfg->regressfp = fopen(esl_opt_GetString(go, "--regress"), "w")) == NULL) 
	cm_Fail("Failed to open --regress output file %s\n", esl_opt_GetString(go, "--regress"));
    }

  /* determine the number of stages and beta and tau values for each stage */
  cfg->beta = NULL;
  cfg->tau  = NULL;
  int s;
  if((! (esl_opt_IsDefault(go, "--betas"))) && (! (esl_opt_IsDefault(go, "--betae"))))
    {
      cfg->nstages = esl_opt_GetInteger(go, "--betae") - esl_opt_GetInteger(go, "--betas") + 2;
      ESL_ALLOC(cfg->beta, sizeof(double) * cfg->nstages);
      cfg->beta[0] = 0.; /* this won't matter b/c stage 1 is non-QDB */
      cfg->beta[1] = pow(10., (-1. * esl_opt_GetInteger(go, "--betas")));
      for(s = 2; s < cfg->nstages; s++)
	cfg->beta[s] = cfg->beta[(s-1)] / 10.;
    }
  else if((! (esl_opt_IsDefault(go, "--taus"))) && (! (esl_opt_IsDefault(go, "--taue"))))
    {
      cfg->nstages = esl_opt_GetInteger(go, "--taue") - esl_opt_GetInteger(go, "--taus") + 2;
      ESL_ALLOC(cfg->tau, sizeof(double) * cfg->nstages);
      cfg->tau[1] = pow(10., (-1. * esl_opt_GetInteger(go, "--taus")));
      for(s = 2; s < cfg->nstages; s++)
	cfg->tau[s] = cfg->tau[(s-1)] / 10.;
    }
  else
    {
      cfg->nstages = 2;
      if(esl_opt_GetBoolean(go, "--qdb") || esl_opt_GetBoolean(go, "--qdbsmall"))
	{
	  ESL_ALLOC(cfg->beta, sizeof(double) * cfg->nstages);
	  cfg->beta[0] = 0.; /* this won't matter b/c stage 1 is non-QDB */
	  cfg->beta[1] = esl_opt_GetReal(go, "--beta");
	}
      else if(esl_opt_GetBoolean(go, "--hbanded"))
	{
	  ESL_ALLOC(cfg->tau, sizeof(double) * cfg->nstages);
	  cfg->tau[0] = 0.; /* this won't matter b/c stage 1 is non-banded */
	  cfg->tau[1] = esl_opt_GetReal(go, "--tau");
	}
    }  
  cfg->s1_sc = NULL; /* alloc'ed and filled in first call of output_result */
  cfg->s1_w  = esl_stopwatch_Create();
  cfg->w     = esl_stopwatch_Create();
  return eslOK;

 ERROR:
  return status;
}
/* serial_master()
 * The serial version of cmscore.
 * Score each sequence against the CM with specified
 * scoring algorithms.
 * 
 * A master can only return if it's successful. All errors are handled immediately and fatally with cm_Fail().
 */
static void
serial_master(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int      status;
  char     errbuf[eslERRBUFSIZE];
  ESL_SQ **sq;
  int       i;
  Parsetree_t    **tr;          /* parse trees for the sequences for stage 2 -> nstages */
  CP9trace_t     **cp9_tr;      /* CP9 traces for the sequences for stage 2 -> nstages */
  ESL_SQFILE      *sqfp;        /* open sequence file stream, we wastefully read the sequences
				 * at each stage, opening and closing the seq file  */

  
  if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);
  /* init_shared_cfg UNNEC? */
  /*if ((status = init_shared_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);*/

  /* Align sequences cfg->nstages times */
  for(cfg->s = 0; cfg->s < cfg->nstages; cfg->s++) 
    {
      /* open input sequence file */
      status = esl_sqfile_Open(cfg->sqfile, cfg->fmt, NULL, &(sqfp));
      if (status == eslENOTFOUND)    cm_Fail("File %s doesn't exist or is not readable\n", cfg->sqfile);
      else if (status == eslEFORMAT) cm_Fail("Couldn't determine format of sequence file %s\n", cfg->sqfile);
      else if (status == eslEINVAL)  cm_Fail("Can’t autodetect stdin or .gz."); 
      else if (status != eslOK)      cm_Fail("Sequence file open failed with error %d\n", status);
      cfg->fmt = sqfp->format;
      
      /* Start timing. */
      if(cfg->s == 0) esl_stopwatch_Start(cfg->s1_w);
      else            esl_stopwatch_Start(cfg->w);

      /* initialize the flags/options/params of the CM for the current stage */
      initialize_cm(go, cfg, errbuf);
      /* Configure the CM for alignment based on cm->config_opts and cm->align_opts.
       * set local mode, make cp9 HMM, calculate QD bands etc. */
      ConfigCM(cfg->cm, NULL, NULL);
      /* Align the sequences for this stage, there's no 'process_workunit() function' 
       * but rather different functions for serial/mpi */ 
      if(cfg->s > 0 && esl_opt_GetBoolean(go, "--scoreonly")) /* only score the seqs, don't get parsetrees */
	{
	  float *scoreonly_sc;
	  serial_align_targets(sqfp, cfg->cm, &sq, NULL, NULL, NULL, &cfg->nseq, 
			       &scoreonly_sc, 0, 0, (! esl_opt_GetBoolean(go, "-i")));
	  
	  /* stop timing */
	  if(cfg->s == 0) esl_stopwatch_Stop(cfg->s1_w);
	  else            esl_stopwatch_Stop(cfg->w);
	  if ((status = output_result(go, cfg, errbuf, sq, tr, NULL, scoreonly_sc))        != eslOK) cm_Fail(errbuf);
	  free(scoreonly_sc);
	}
      else if(cfg->s > 0 && ((esl_opt_GetBoolean(go, "--hmmonly")))) /* hmm alignment */
	{
	  serial_align_targets(sqfp, cfg->cm, &sq, NULL, NULL, &cp9_tr, &cfg->nseq, 
			       NULL, 0, 0, (! esl_opt_GetBoolean(go, "-i")));
	  /* stop timing */
	  if(cfg->s == 0) esl_stopwatch_Stop(cfg->s1_w);
	  else            esl_stopwatch_Stop(cfg->w);
	  if ((status = output_result(go, cfg, errbuf, sq, NULL, cp9_tr, NULL))        != eslOK) cm_Fail(errbuf);
	  for (i = 0; i < cfg->nseq; i++) CP9FreeTrace(cp9_tr[i]);
	  free(cp9_tr);
	}
      else /* CYK alignment of some sort */
	{
	  serial_align_targets(sqfp, cfg->cm, &sq, &tr, NULL, NULL, &cfg->nseq, 
			       NULL, 0, 0, (! esl_opt_GetBoolean(go, "-i")));
	  /* stop timing */
	  if(cfg->s == 0) esl_stopwatch_Stop(cfg->s1_w);
	  else            esl_stopwatch_Stop(cfg->w);
	  if ((status = output_result(go, cfg, errbuf, sq, tr, NULL, NULL))        != eslOK) cm_Fail(errbuf);
	  for (i = 0; i < cfg->nseq; i++) FreeParsetree(tr[i]);
	  free(tr);
	}

      /* free sequences, we wastefully read them in separately (within serial_align_targets()) for each stage */
      for (i = 0; i < cfg->nseq; i++) esl_sq_Destroy(sq[i]);
      free(sq);
      esl_sqfile_Close(sqfp);
    }
}


static int
output_result(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, ESL_SQ **sq, Parsetree_t **tr, CP9trace_t **cp9_tr, 
	      float *scoreonly_sc)
{
  int status;
  int i;

  /* print the parsetrees to regression file or parse file */
  for(i = 0; i < cfg->nseq; i++)
    {
      for(i = 0; i < cfg->nseq; i++)
	{
	  if (cfg->regressfp != NULL) 
	    {
	      fprintf(cfg->regressfp, "> %s\n", sq[i]->name);
	      if(esl_opt_GetBoolean(go,"--hmmonly")) 
		{
		  fprintf(cfg->regressfp, "  SCORE : %.2f bits\n", CP9TraceScore(cfg->cm->cp9, sq[i], cp9_tr[i]));
		  CP9PrintTrace(cfg->regressfp, cp9_tr[i], cfg->cm->cp9, sq[i]);
		}
	      else
		{
		  fprintf(cfg->regressfp, "  SCORE : %.2f bits\n", ParsetreeScore(cfg->cm, tr[i], sq[i], FALSE));
		  ParsetreeDump(cfg->regressfp, tr[i], cfg->cm, sq[i], NULL, NULL); /* NULLs are dmin, dmax */
		}
	      fprintf(cfg->regressfp, "//\n");
	    }
	  if (cfg->tracefp   != NULL) 
	    {
	      fprintf(cfg->tracefp, "> %s\n", sq[i]->name);
	      if(esl_opt_GetBoolean(go,"--hmmonly")) 
		{
		  fprintf(cfg->tracefp, "  SCORE : %.2f bits\n", CP9TraceScore(cfg->cm->cp9, sq[i], cp9_tr[i]));
		  CP9PrintTrace(cfg->tracefp, cp9_tr[i], cfg->cm->cp9, sq[i]);
		}
	      else
		{
		  fprintf(cfg->tracefp, "  SCORE : %.2f bits\n", ParsetreeScore(cfg->cm, tr[i], sq[i], FALSE));
		  ParsetreeDump(cfg->tracefp, tr[i], cfg->cm, sq[i], NULL, NULL); /* NULLs are dmin, dmax */
		}
	      fprintf(cfg->tracefp, "//\n");
	    }
	}
    }

  /* print info about scores of parsetrees */
  if(cfg->s == 0) /* store the scores, only */
    {
      ESL_ALLOC(cfg->s1_sc, sizeof(float) * cfg->nseq);
      for(i = 0; i < cfg->nseq; i++)
	cfg->s1_sc[i] = ParsetreeScore(cfg->cm, tr[i], sq[i], FALSE);
    }
  else /* if(cfg->s > 0) we don't do the comparison test for stage 0 */
    {
      /* Compare parsetrees from stage 1 and stage s (current stage) and collect stats */
      float *sc;           /* parse scores for current stage of alignment */
      double diff_sc = 0.; /* difference in summed parse scores for this stage versus stage 1 */
      int    diff_ct = 0.; /* number of parses different between this stage and stage 1 */
      double spdup;        /* speed-up versus stage 1 */

      ESL_ALLOC(sc, (sizeof (float) * cfg->nseq));
      for(i = 0; i < cfg->nseq; i++)
	{
	  /* TO DO: write function that in actually_align_targets(), takes
	   * a CP9 parse, and converts it to a CM parsetree */
	  if(esl_opt_GetBoolean(go, "--hmmonly"))
	    sc[i] = CP9TraceScore(cfg->cm->cp9, sq[i], cp9_tr[i]);
	  else if (esl_opt_GetBoolean(go, "--scoreonly"))
	    sc[i] = scoreonly_sc[i];
	  else
	    sc[i] = ParsetreeScore(cfg->cm, tr[i], sq[i], FALSE);
	  if(esl_opt_GetBoolean(go, "-i"))
	    printf("%-12s S1: %.3f S%d: %.3f diff: %.3f\n", sq[i]->name, cfg->s1_sc[i], (cfg->s+1), sc[i], (fabs(cfg->s1_sc[i] - sc[i])));
	  if(fabs(cfg->s1_sc[i] -  sc[i]) > 0.0001)
	    {
	      diff_ct++;
	      diff_sc += fabs(cfg->s1_sc[i] - sc[i]);
	    }
	}
      /* Print summary for this stage versus stage 1 */ 
      printf("Results summary for stage %d:\n", (cfg->s+1));
      printf("---------------------------------\n");
      printf("Number seqs aligned:     %d\n", cfg->nseq);
      esl_stopwatch_Display(stdout, cfg->s1_w, "Stage  1 time:           ");
      printf("Stage %2d time:           ", (cfg->s+1));
      esl_stopwatch_Display(stdout, cfg->w, "");
      spdup = cfg->s1_w->user / cfg->w->user;
      printf("Speedup (user):          %.2f\n", spdup);

      if(! esl_opt_GetBoolean(go, "--hmmonly"))
	{
	  printf("Avg bit score diff:      %.2f\n", (diff_sc / ((float) cfg->nseq)));
	  if(diff_ct == 0)
	    printf("Avg sc diff(>1e-4):      %.2f\n", 0.);
	  else
	    printf("Avg sc diff(>1e-4):      %.2f\n", (diff_sc / ((float) diff_ct)));
	  printf("Num   diff (>1e-4):      %d\n", (diff_ct));
	  printf("Fract diff (>1e-4):      %.5f\n", (((float) diff_ct) / ((float) cfg->nseq)));
	  printf("\n\n");
	}
      free(sc);
    }
  return eslOK;

 ERROR:
  return status;
}

/* initialize_cm()
 * Setup the CM based on the command-line options/defaults
 * for the specified stage alignment. We only set flags and 
 * a few parameters. ConfigCM() configures the CM.
 */
static int
initialize_cm(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf)
{
  /* set up params/flags/options of the CM */
  if(cfg->beta != NULL) cfg->cm->beta = cfg->beta[cfg->s];
  if(cfg->tau  != NULL) cfg->cm->tau  = cfg->tau[cfg->s];

  cfg->cm->align_opts = 0;  /* clear alignment options from previous stage */
  cfg->cm->config_opts = 0; /* clear configure options from previous stage */
  /* Clear QDBs if they exist */
  if(cfg->cm->flags & CM_QDB)
    {
      free(cfg->cm->dmin);
      free(cfg->cm->dmax);
      cfg->cm->dmin = NULL;
      cfg->cm->dmax = NULL;
      cfg->cm->flags &= ~CM_QDB;
    }

  /* enable option to check parsetree score against the alignment score */
  cfg->cm->align_opts  |= CM_ALIGN_CHECKPARSESC;

  /* Update cfg->cm->config_opts and cfg->cm->align_opts based on command line options */
  if(esl_opt_GetBoolean(go, "-l"))
    {
      cfg->cm->config_opts |= CM_CONFIG_LOCAL;
      cfg->cm->config_opts |= CM_CONFIG_HMMLOCAL;
      cfg->cm->config_opts |= CM_CONFIG_HMMEL;
    }
  if(esl_opt_GetBoolean(go, "--sub"))         
    { 
      cfg->cm->align_opts  |=  CM_ALIGN_SUB;
      cfg->cm->align_opts  &= ~CM_ALIGN_CHECKPARSESC; /* parsetree score won't match aln score */
    }
  if(esl_opt_GetBoolean(go, "-i"))            cfg->cm->align_opts  |= CM_ALIGN_TIME;
  if(esl_opt_GetBoolean(go, "--zeroinserts")) cfg->cm->config_opts |= CM_CONFIG_ZEROINSERTS;
  if(esl_opt_GetBoolean(go, "--elsilent"))    cfg->cm->config_opts |= CM_CONFIG_ELSILENT;

  if(cfg->s == 0) /* set up stage 1 alignment we'll compare all other stages to */
    {
      /* only one option allows cmscore NOT to do standard CYK as first stage aln */
      if(esl_opt_GetBoolean(go, "--qdbboth")) { 
	cfg->cm->align_opts  |= CM_ALIGN_QDB;
	cfg->cm->config_opts |= CM_CONFIG_QDB;
      }
    }      
  else { /* cfg->s > 0, we're at least on stage 2 */
    if(esl_opt_GetBoolean(go, "--hbanded"))     cfg->cm->align_opts  |= CM_ALIGN_HBANDED;
    if(esl_opt_GetBoolean(go, "--hmmonly"))     cfg->cm->align_opts  |= CM_ALIGN_HMMONLY;
    if(esl_opt_GetBoolean(go, "--hsafe"))       cfg->cm->align_opts  |= CM_ALIGN_HMMSAFE;
    if(esl_opt_GetBoolean(go, "--scoreonly"))   cfg->cm->align_opts  |= CM_ALIGN_SCOREONLY;
    if(esl_opt_GetBoolean(go, "--qdb") || esl_opt_GetBoolean(go, "--qdbsmall"))                    
      { 
	cfg->cm->align_opts  |= CM_ALIGN_QDB;
	cfg->cm->config_opts |= CM_CONFIG_QDB;
      }
    /* only 1 way stage 2+ alignment will be D&C, if --qdbsmall was enabled */
    if(! esl_opt_GetBoolean(go, "--qdbsmall"))
      cfg->cm->align_opts  |= CM_ALIGN_NOSMALL;
  }
  printf("Stage %2d alignment:\n", (cfg->s+1));
  summarize_align_options(cfg->cm);
  return eslOK;
}

/* Function: summarize_align_options
 * Date:     EPN, Wed Jan 17 09:08:18 2007
 * Purpose:  Print out alignment options in pretty format. 
 */
int summarize_align_options(CM_t *cm)
{
  printf("---------------------------------\n");
  /* Algorithm */
  if(cm->align_opts & CM_ALIGN_INSIDE)
    printf("Algorithm:               Inside\n");
  else if(cm->align_opts & CM_ALIGN_HMMONLY) 
    printf("Algorithm:               CP9 HMM Viterbi\n");
  else if(cm->align_opts & CM_ALIGN_SCOREONLY)
    printf("Algorithm:               CYK Standard (score only)\n");
  else if(cm->align_opts & CM_ALIGN_NOSMALL)
    printf("Algorithm:               CYK Standard\n");
  else 
    printf("Algorithm:               CYK D&C\n");

  /* Bands */
  if(cm->align_opts & CM_ALIGN_HBANDED)
    {
      if(cm->align_opts & CM_ALIGN_SUMS)
	printf("Bands:                   CP9 HMM (sums)\n");
      else
	printf("Bands:                   CP9 HMM\n"); 
      printf("Tail loss:               %g\n", cm->tau);
    }
  else if(cm->align_opts & CM_ALIGN_QDB)
    {
      printf("Bands:                   QDB\n"); 
      printf("Tail loss:               %g\n", cm->beta);
    }
  else
    {
      printf("Bands:                   None\n"); 
      printf("Tail loss:               0.0\n");
    }

  /* Locality */
  if(cm->config_opts & CM_CONFIG_LOCAL)
    printf("Local:                   Yes\n");
  else
    printf("Local:                   No\n");

  /* Sub mode? */
  if(cm->align_opts & CM_ALIGN_SUB)
    printf("Sub mode.\n");
  if(cm->align_opts & CM_ALIGN_FSUB)
    printf("Full Sub mode.\n");

  printf("---------------------------------\n");
  return eslOK;
}




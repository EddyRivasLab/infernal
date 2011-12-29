/* cmcalibrate: score a CM against random sequence data to set 
 *              the statistical parameters for E-value determination.
 * 
 * EPN, Thu Dec 15 19:40:16 2011 [Updated for v1.1 (removed CP9 and fthr)]
 * EPN, Wed May  2 07:02:52 2007
 * based on HMMER-2.3.2's hmmcalibrate.c from SRE
 */

#include "esl_config.h"
#include "p7_config.h"
#include "config.h"	

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <time.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_histogram.h"
#include "esl_mpi.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_ratematrix.h"
#include "esl_stack.h"
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

#define REALLYSMALLX        1e-20
#define EXPTAIL_CHUNKLEN    10000 /* sequence chunk length for random sequence searches */
///#define EXPTAIL_CHUNKLEN    1000 /* sequence chunk length for random sequence searches */
#define DEBUGMPI 0

typedef struct {
#ifdef HMMER_THREADS
  ESL_WORK_QUEUE   *queue;
#endif /*HMMER_THREADS*/
  CM_t             *cm;        /* a covariance model        */
  float            *scA;       /* vector of scores of hits  */
  int64_t           nhits;     /* number of hits in scA     */
  float             cutoff;    /* minimum hit score to keep */
} WORKER_INFO;

#if defined (HMMER_THREADS) && defined (HAVE_MPI)
#define CPUOPTS     "--mpi"
#define MPIOPTS     "--cpu"
#else
#define CPUOPTS     NULL
#define MPIOPTS     NULL
#endif

static ESL_OPTIONS options[] = {
  /* name                  type   default   env        range      toggles      reqs       incomp  help  docgroup*/
  { "-h",           eslARG_NONE,    FALSE,  NULL,      NULL,      NULL,        NULL,        NULL, "show brief help on version and usage",   1 },
  { "-L",           eslARG_REAL,    "1.6",  NULL,"0.01<=x<=100.",  NULL,        NULL,         NULL, "set random seq length to search in Mb to <x>", 2 },
  { "--forecast",   eslARG_INT,     NULL,   NULL,      NULL,      NULL,        NULL,        NULL, "don't do calibration, forecast running time with <n> processors", 1 },
  { "--seed",       eslARG_INT,     "181",  NULL,      "n>=0",    NULL,        NULL,        NULL, "set RNG seed to <n> (if 0: one-time arbitrary seed)", 1 },
  { "--beta",       eslARG_REAL,    "1E-15",NULL,      "x>0",     NULL,        NULL,"--nonbanded","set tail loss prob for QDB to <x>", 1 },
  { "--nonbanded",  eslARG_NONE,    FALSE,  NULL,      NULL,      NULL,        NULL,        NULL, "do not use QDBs for calibrating CM search modes", 1 },
  { "--gtailn",     eslARG_INT,     "250",  NULL,   "n>=100",     NULL,        NULL,   "--tailp", "fit the top <n> hits/Mb in histogram for  CM glocal modes", 1 },
  { "--ltailn",     eslARG_INT,     "750",  NULL,   "n>=100",     NULL,        NULL,   "--tailp", "fit the top <n> hits/Mb in histogram for  CM  local modes", 1 },
  { "--tailp",      eslARG_REAL,    NULL,   NULL,"0.0<x<0.6",     NULL,        NULL,        NULL, "set fraction of histogram tail to fit to exp tail to <x>", 1 },
  { "--hfile",      eslARG_OUTFILE, NULL,   NULL,      NULL,      NULL,        NULL,        NULL, "save fitted score histogram(s) to file <f>", 1 },
  { "--sfile",      eslARG_OUTFILE, NULL,   NULL,      NULL,      NULL,        NULL,        NULL, "save survival plot to file <f>", 1 },
  { "--qqfile",     eslARG_OUTFILE, NULL,   NULL,      NULL,      NULL,        NULL,        NULL, "save Q-Q plot for score histograms to file <f>", 1 },
  { "--ffile",      eslARG_OUTFILE, NULL,   NULL,      NULL,      NULL,        NULL,        NULL, "save lambdas for different tail fit probs to file <f>", 1 },
  { "--devhelp",    eslARG_NONE,    NULL,   NULL,      NULL,      NULL,        NULL,        NULL, "show list of undocumented developer options", 1 },
#ifdef HMMER_THREADS 
  { "--cpu",        eslARG_INT,     NULL,"HMMER_NCPU", "n>=0",    NULL,        NULL,     CPUOPTS, "number of parallel CPU workers to use for multithreads", 1 },
#endif
#ifdef HAVE_MPI
  { "--mpi",        eslARG_NONE,    FALSE,  NULL,      NULL,      NULL,        NULL,     MPIOPTS, "run as an MPI parallel program", 1 },  
  { "--stall",      eslARG_NONE,    FALSE,  NULL,      NULL,      NULL,        NULL,        NULL, "arrest after start: for debugging MPI under gdb", 1 },  
#endif
  /* Developer options the average user doesn't need to know about (only shown if --devhelp) */
  { "-x",           eslARG_NONE,   FALSE,   NULL, NULL,     NULL,        NULL, "--forecast", "print arguably interesting info",  101},
  { "-v",           eslARG_NONE,   FALSE,   NULL, NULL,     NULL,        NULL, "--forecast", "print arguably interesting info",  101},
  { "-T",           eslARG_REAL,    NULL,   NULL, NULL,     NULL,        NULL,        NULL, "set bit sc cutoff for exp tail fitting to <x> [df: -INFTY]", 101 },
  { "--random",     eslARG_NONE,    NULL,   NULL, NULL,     NULL,        NULL,        NULL, "use GC content of random null background model of CM",  101},
  { "--gc",         eslARG_INFILE,  NULL,   NULL, NULL,     NULL,        NULL,        NULL, "use GC content distribution from file <f>",  101},
  { "--nonull3",    eslARG_NONE,   FALSE,   NULL, NULL,      NULL,        NULL,        NULL, "turn OFF the NULL3 post hoc additional null model", 101 },
#ifdef HAVE_MPI
  /* Developer option, for debugging */
#endif
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

struct cfg_s {
  char              *cmfile;	         /* name of input CM file  */ 
  ESL_RANDOMNESS    *r;                  /* random number generator, for sequences to search */
  ESL_RANDOMNESS    *r_est;              /* random number generator, for sequences used in time estimates */
  ESL_ALPHABET      *abc;                /* alphabet */
  ESL_STOPWATCH     *w_stage;            /* stopwatch for each stage */
  double            *gc_freq;            /* gc frequence [0..100], only used if --gc */
  ExpInfo_t       ***expAA;              /* the exponential tail info, 1st dim: 1 for each CM, 2nd dim: EXP_NMODES */
  int                ncm;                /* what number CM we're on */
  int                be_verbose;	 /* print extra info? */
  int                cmalloc;            /* number of expAA we have allocated (1st dim) */
  char              *tmpfile;            /* tmp file we're writing to */
  char              *mode;               /* write mode, "w" or "wb"                     */
  double            *dnull;              /* double version of cm->null, for generating random seqs */
  float              sc_cutoff;          /* minimum score of a hit we'll consider (-eslINFINITY by default) */

  /* the HMM that generates sequences for exponential tail fitting */
  int                ghmm_nstates;       /* number of states in the HMM */
  double            *ghmm_sA;            /* start probabilities [0..ghmm_nstates-1] */
  double           **ghmm_tAA;           /* transition probabilities [0..nstates-1][0..nstates-1] */
  double           **ghmm_eAA;           /* emission probabilities   [0..nstates-1][0..abc->K-1] */

  /* number of sequences and the length of each seq for exp tail
   * fitting, set such that: exp_cmN is the number of 10 Kb
   * seqs we'll search for CM local/glocal exponential tail fitting:
   *
   * N = (esl_opt_GetBoolean(go, "-L")  * 10^6) / EXPTAIL_CHUNKLEN(10000); 
   *
   * We don't search just 1 long sequence (i.e. 1.5 Mb) b/c using
   * sequence lengths above 10 Kb for exp tail calibration can yield
   * millions of hits (for CM searches) before overlaps are removed,
   * which requires a lot of memory.
   */
  int              N;        /* number of 10 Kb seqs for  local CM exp tail fitting */
  int              L;        /* the size of seq chunks to search, set as 10,000 (10 Kb) */

  /* info for the comlog we'll add to the cmfiles */
  char            *ccom;               /* command line used in this execution of cmcalibrate */
  char            *cdate;              /* date of this execution of cmcalibrate */

  /* mpi */
  int              do_mpi;
  int              my_rank;
  int              nproc;
  int              do_stall;          /* TRUE to stall the program until gdb attaches */

  /* Masters only (i/o streams) */
  CM_FILE         *cmfp;	      /* open input CM file stream       */
  FILE            *exphfp;            /* optional output for exp tail histograms */
  FILE            *expsfp;            /* optional output for exp tail survival plot */
  FILE            *expqfp;            /* optional output for exp tail QQ file */
  FILE            *exptfitfp;         /* optional output for exp tail fit file */
};

static char usage[]  = "[-options] <cmfile>";
static char banner[] = "fit exponential tails for CM E-values";

static void serial_master(ESL_GETOPTS *go, struct cfg_s *cfg);
static int  serial_loop  (WORKER_INFO *info, char *errbuf, ESL_SQ_BLOCK *sq_block);

/* Functions to avoid code duplication for common tasks */
static int  init_master_cfg (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int  init_shared_cfg (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int  initialize_cm   (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, int do_local);
static int  initialize_stats(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int  fit_histogram(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, float *scores, int nscores, int exp_mode, double *ret_mu, double *ret_lambda, int *ret_nrandhits, float *ret_tailp);
static int  get_random_dsq(const struct cfg_s *cfg, char *errbuf, CM_t *cm, int L, ESL_RANDOMNESS *r, ESL_DSQ **ret_dsq);
static int  get_cmcalibrate_comlog_info(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);
static int  update_comlog(const ESL_GETOPTS *go, char *errbuf, char *ccom, char *cdate, CM_t *cm);
static int  set_dnull(struct cfg_s *cfg, CM_t *cm, char *errbuf);
static int  print_run_info(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf);
static int  get_command(const ESL_GETOPTS *go, char *errbuf, char **ret_command);
static int  print_per_cm_column_headings(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm);
static int  print_per_cm_summary        (const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, double psec, double asec);
static int  print_exp_line(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, int exp_mode, int N, int L, double psec);
static int  print_post_calibration_info (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, FILE *fp, CM_t *cm, double *psecA, double *asecA);
static int  estimate_time_for_exp_round (const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, int exp_mode, double *ret_sec_per_res);
static int  process_search_workunit(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float cutoff, CM_TOPHITS **ret_th);
static int  generate_sequences(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, ESL_SQ_BLOCK **ret_sq_block);
static int  output_header(FILE *ofp, const ESL_GETOPTS *go, char *cmfile);
static void process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_cmfile);

#ifdef HMMER_THREADS
#define BLOCK_SIZE 150

static int  thread_loop(WORKER_INFO *info, char *errbuf, ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQ_BLOCK *sq_block);
static void pipeline_thread(void *arg);
#endif /*HMMER_THREADS*/

#ifdef HAVE_MPI
static int  mpi_master   (ESL_GETOPTS *go, struct cfg_s *cfg);
static int  mpi_worker   (ESL_GETOPTS *go, struct cfg_s *cfg);
#endif /*HAVE_MPI*/


#ifdef HAVE_MPI

/* Define common tags used by the MPI master/slave processes */
#define INFERNAL_ERROR_TAG          1
#define INFERNAL_SEQIDX_TAG         2
#define INFERNAL_SCORES_TAG         3
#define INFERNAL_TERMINATING_TAG    4
#define INFERNAL_READY_TAG          5

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
 * The MPI version of cmcalibrate.
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
  int      status;                /* Easel status */
  int      qhstatus = eslOK;      /* status from cm_file_Read() */
  char     errbuf[cmERRBUFSIZE];  /* for printing error messages */
  CM_t    *cm    = NULL;          /* the CM */
  CM_t    *nc_cm = NULL;          /* a non-configured copy of the CM */
  int      i, h;                  /* counters */
  int      cmi;                   /* CM index, which model we're working on */
  char     time_buf[128];	  /* string for printing elapsed time (safely holds up to 10^14 years) */
  int      exp_mode;              /* ctr over exp tail modes */
  double   cm_psec;               /* predicted number of seconds for calibrating current CM */
  double   cm_asec;               /* predicted number of seconds for calibrating current CM */
  double   total_psec = 0.;       /* predicted number of seconds for calibrating all CMs */
  double   total_asec = 0.;       /* predicted number of seconds for calibrating all CMs */
  double   asecA[EXP_NMODES];     /* stores actual timings for each exp tail fit stage, for each CM */
  double   psecA[EXP_NMODES];     /* stores predicted timings for each exp tail fit stage, for each CM */
  double   psec;                  /* predicted seconds */
  int64_t  merged_nhits = 0;      /* number of hits reported thus far, for all seqs */
  float   *merged_scA = NULL;     /* [0..merged_nhits-1] hit scores for all seqs */
  double   tmp_mu, tmp_lambda;    /* temporary mu and lambda used for setting exp tails */
  int      tmp_nrandhits;         /* temporary number of rand hits found */
  float    tmp_tailp;             /* temporary tail mass probability fit to an exponential */

  MPI_Status       mpistatus;       /* the mpi status */
  char            *mpibuf  = NULL;  /* buffer used to pack/unpack structures */
  int              mpibuf_size = 0; /* current size of mpibuf_size */
  int              buf_size;        /* size of received buffer */
  int              dest;            /* destination for MPISend() */
  int              pos;             /* for packing/unpacking an MPI buffer */
  int64_t          nhits;           /* nhits received from a worker */
  
  if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);
  if ((status = print_run_info (go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);

  /* <cfg->abc> is not known 'til first CM is read. Could be DNA or RNA*/
  qhstatus = cm_file_Read(cfg->cmfp, TRUE, &(cfg->abc), &cm);

  if (qhstatus == eslOK) {
    /* One-time initializations after alphabet <abc> becomes known */
    output_header(stdout, go, cfg->cmfile);

    /* master needs to be able to generate sequences only for timing predictions */
    if((status = CreateGenomicHMM(cfg->abc, errbuf, &(cfg->ghmm_sA), &(cfg->ghmm_tAA), &(cfg->ghmm_eAA), &cfg->ghmm_nstates)) != eslOK) cm_Fail("unable to create generative HMM\n%s", errbuf);
    if((status = set_dnull(cfg, cm, errbuf)) != eslOK) cm_Fail("unable to create set_dnull\n%s\n", errbuf);
  }

  /* Outer loop: over each query CM in <cmfile>. */
  while (qhstatus == eslOK) {
    /* clone the non-configured CM we just read, we'll come back to it when we switch from global to local */
    if((status = cm_Clone(cm, errbuf, &nc_cm)) != eslOK) cm_Fail("unable to clone CM");

    cfg->ncm++;
    if(cfg->ncm == cfg->cmalloc) { /* expand our memory */
      cfg->cmalloc  += 128;
      ESL_REALLOC(cfg->expAA, sizeof(ExpInfo_t **) * cfg->cmalloc);
    }
    cmi = cfg->ncm-1;
    cm_asec = cm_psec = 0.;

    if((status = initialize_cm(go, cfg, errbuf, cm, FALSE))          != eslOK) cm_Fail(errbuf);
    if((status = initialize_stats(go, cfg, errbuf))                  != eslOK) cm_Fail(errbuf);
    if((status = print_per_cm_column_headings(go, cfg, errbuf, cm))  != eslOK) cm_Fail(errbuf);

    for(exp_mode = 0; exp_mode < EXP_NMODES; exp_mode++) {
      /* do we need to switch from global configuration to local? */
      if(exp_mode > 0 && (! ExpModeIsLocal(exp_mode-1)) && ExpModeIsLocal(exp_mode)) {
	/* switch from global to local by copying the current exptail stats from <cm>
	 * into <nc_cm> and then configure <nc_cm> for local mode. We do it this
	 * way because as a rule we don't allow reconfiguration of CMs (to limit
	 * execution paths through configuration functions)
	 */
	FreeCM(cm);
	cm = nc_cm;
	if((status = initialize_cm(go, cfg, errbuf, cm, TRUE)) != eslOK) cm_Fail(errbuf);
      }
      /* update search info for round 0 (final round) for exp tail mode */
      if(ExpModeIsInside(exp_mode)) cm->search_opts |= CM_SEARCH_INSIDE;
      else                          cm->search_opts &= ~CM_SEARCH_INSIDE;

      /* estimate time for this round */
      if((status = estimate_time_for_exp_round(go, cfg, errbuf, cm, exp_mode, &psec)) != eslOK) cm_Fail(errbuf); 
      psec *= cfg->N * cfg->L; /* psec was per residue */
      if(cfg->nproc > 1) psec /= (cfg->nproc-1); /* parallelization will speed us up */
      psecA[exp_mode] = psec;
      cm_psec    += psec;
      total_psec += psec;
      print_exp_line(go, cfg, errbuf, exp_mode, cfg->N, cfg->L, psec);

      esl_stopwatch_Start(cfg->w_stage);
      fflush(stdout);

      /* for each sequence i from 0..N-1, tell the appropriate worker (widx = (i+1) % nworkers) to search it */
      for(i = 0; i < cfg->N; i++) { 
	if (MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus) != 0) mpi_failure("MPI error %d receiving message from %d\n", mpistatus.MPI_SOURCE);
	    
	MPI_Get_count(&mpistatus, MPI_PACKED, &buf_size);
	if (mpibuf == NULL || buf_size > mpibuf_size) {
	  ESL_REALLOC(mpibuf, sizeof(char) * buf_size);
	  mpibuf_size = buf_size; 
	}
	dest = mpistatus.MPI_SOURCE;
	MPI_Recv(mpibuf, buf_size, MPI_PACKED, dest, mpistatus.MPI_TAG, MPI_COMM_WORLD, &mpistatus);
	
	if (mpistatus.MPI_TAG == INFERNAL_ERROR_TAG)  mpi_failure("MPI client %d raised error:\n%s\n", dest, mpibuf);
	if (mpistatus.MPI_TAG != INFERNAL_READY_TAG)  mpi_failure("Unexpected tag %d from %d\n", mpistatus.MPI_TAG, dest);

#if DEBUGMPI
	printf("master rec'd INFERNAL_READY_TAG from worker %d\n", dest);
#endif
	
	/* send the int 'i' to the worker, telling it to search his copy of the i'th sequence for this CM/mode pair */
	if (MPI_Send(&i, 1, MPI_INT, dest, INFERNAL_SEQIDX_TAG, MPI_COMM_WORLD) != 0) mpi_failure("mpi send failed");
#if DEBUGMPI
	printf("master sent i: %d to worker %d\n", i, dest);
#endif
      }

      /* wait for all workers to finish up their work blocks */
      for (i = 1; i < cfg->nproc; ++i) { 
	if (MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus) != 0) mpi_failure("MPI error %d receiving message from %d\n", mpistatus.MPI_SOURCE);

	MPI_Get_count(&mpistatus, MPI_PACKED, &buf_size);
	if (mpibuf == NULL || buf_size > mpibuf_size) {
	  ESL_REALLOC(mpibuf, sizeof(char) * buf_size);
	  mpibuf_size = buf_size; 
	}
	dest = mpistatus.MPI_SOURCE;
	MPI_Recv(mpibuf, buf_size, MPI_PACKED, dest, mpistatus.MPI_TAG, MPI_COMM_WORLD, &mpistatus);

	if (mpistatus.MPI_TAG == INFERNAL_ERROR_TAG)  mpi_failure("MPI client %d raised error:\n%s\n", dest, mpibuf);
	if (mpistatus.MPI_TAG != INFERNAL_READY_TAG)  mpi_failure("Unexpected tag %d from %d\n", mpistatus.MPI_TAG, dest);
#if DEBUGMPI
	printf("master done sending indices, recived INFERNAL_READY_TAG from worker %d\n", dest);
#endif 
      }

      ESL_STOPWATCH *merge_w = esl_stopwatch_Create();
      esl_stopwatch_Start(merge_w);
      /* merge the arrays of scores we get back from the workers */
      for (dest = 1; dest < cfg->nproc; ++dest) { 
	/* send -1 to the worker, telling it we're done with this CM/mode pair */
	i = -1;
	if (MPI_Send(&i, 1, MPI_INT, dest, INFERNAL_SEQIDX_TAG, MPI_COMM_WORLD) != 0) ESL_EXCEPTION(eslESYS, "mpi send failed");
#if DEBUGMPI
	printf("master sent -1 to worker %d\n", dest);
#endif 

	/* wait for the results */
	if (MPI_Probe(dest, MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus) != 0) mpi_failure("MPI error %d receiving message from %d\n", mpistatus.MPI_SOURCE);
	MPI_Get_count(&mpistatus, MPI_PACKED, &buf_size);
	if (mpibuf == NULL || buf_size > mpibuf_size) {
	  ESL_REALLOC(mpibuf, sizeof(char) * buf_size);
	  mpibuf_size = buf_size; 
	}

	/* receive and unpack the array of hit scores */
	MPI_Recv(mpibuf, buf_size, MPI_PACKED, dest, mpistatus.MPI_TAG, MPI_COMM_WORLD, &mpistatus);
#if DEBUGMPI
	printf("master received results from worker %d\n", dest);
#endif 

	if (mpistatus.MPI_TAG == INFERNAL_ERROR_TAG)  mpi_failure("MPI client %d raised error:\n%s\n", dest, mpibuf);
	if (mpistatus.MPI_TAG != INFERNAL_SCORES_TAG) mpi_failure("Unexpected tag %d from %d\n", mpistatus.MPI_TAG, dest);
	pos = 0;
	if (MPI_Unpack(mpibuf, buf_size, &pos, &nhits, 1, MPI_LONG_LONG_INT, MPI_COMM_WORLD) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed");
	if(nhits > 0) { 
	  if(merged_nhits == 0) ESL_ALLOC  (merged_scA, sizeof(float) * (merged_nhits + nhits));
	  else                  ESL_REALLOC(merged_scA, sizeof(float) * (merged_nhits + nhits));
	  for(h = 0; h < nhits; h++) { 
	    if (MPI_Unpack(mpibuf, buf_size, &pos, &(merged_scA[merged_nhits+h]), 1, MPI_FLOAT, MPI_COMM_WORLD) != 0) ESL_XEXCEPTION(eslESYS, "unpack failed");
	  }
	  merged_nhits += nhits;
#if DEBUGMPI
	  printf("master added %" PRId64 " hits to merged_nhits, now %" PRId64 " total\n", nhits, merged_nhits);
#endif 
	}
      }
      esl_stopwatch_Stop(merge_w);
      esl_stopwatch_Display(stdout, merge_w, "# Merge CPU time: ");

      if(cfg->exptfitfp != NULL) { 
	fprintf(cfg->exptfitfp, "# CM: %s\n", cm->name);
	fprintf(cfg->exptfitfp, "# mode: %12s\n", DescribeExpMode(exp_mode));
      }
#if DEBUGMPI
      printf("about to call fit_histogram, %" PRId64 " hits\n", merged_nhits);
#endif 
      if((status = fit_histogram(go, cfg, errbuf, merged_scA, merged_nhits, exp_mode, &tmp_mu, &tmp_lambda, &tmp_nrandhits, &tmp_tailp)) != eslOK) mpi_failure(errbuf);
      SetExpInfo(cfg->expAA[cmi][exp_mode], tmp_lambda, tmp_mu, (long) (cfg->L * cfg->N), tmp_nrandhits, tmp_tailp);
	
      esl_stopwatch_Stop(cfg->w_stage);
      asecA[exp_mode] = cfg->w_stage->elapsed;
      cm_asec    += cfg->w_stage->elapsed;
      total_asec += cfg->w_stage->elapsed;
      FormatTimeString(time_buf, cfg->w_stage->elapsed, FALSE);
      printf("  %10s\n", time_buf);
      fflush(stdout);

      /* clear hits */
      if(merged_scA != NULL) { free(merged_scA); merged_scA = NULL; }
      merged_nhits = 0;
    } /* end of for(exp_mode = 0; exp_mode < EXP_NMODES; exp_mode++) */

    if(cfg->be_verbose) if((status = debug_print_expinfo_array(cm, errbuf, cfg->expAA[cmi])) != eslOK) mpi_failure(errbuf);
    print_per_cm_summary(go, cfg, errbuf, cm, cm_psec, cm_asec);
    if((status = print_post_calibration_info(go, cfg, errbuf, stdout, cm, psecA, asecA)) != eslOK) mpi_failure(errbuf); 
    FreeCM(cm);

    printf("//\n");
    fflush(stdout);

    /* read the next CM */
    qhstatus = cm_file_Read(cfg->cmfp, TRUE, &(cfg->abc), &cm);
  } /* end of while(qhstatus == eslOK) */
  if(qhstatus != eslEOF) mpi_failure(cfg->cmfp->errbuf);
  
  if(cfg->ncm > 1 && (esl_opt_IsOn(go, "--forecast"))) { 
    fprintf(stdout, "#\n");
    FormatTimeString(time_buf, total_psec, FALSE);
    fprintf(stdout, "# total predicted time for all %d CMs: %s\n", cfg->ncm, time_buf);
    fprintf(stdout, "#\n");
  }
  return eslOK;
      
 ERROR:
  mpi_failure("Memory allocation error.");
  return eslEMEM;
}

/* mpi_worker()
 * 
 * Receive information from the master on which sequences to search
 * and send back hit scores for, and do it.
 */


static int
mpi_worker(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int      status;                /* Easel status */
  int      qhstatus = eslOK;      /* status from cm_file_Read() */
  char     errbuf[cmERRBUFSIZE];  /* for printing error messages */
  CM_t    *cm    = NULL;          /* the CM */
  CM_t    *nc_cm = NULL;          /* a non-configured copy of the CM */
  int      i, h;                  /* counters */
  int      cmi;                   /* CM index, which model we're working on */
  int      exp_mode;              /* ctr over exp tail modes */

  char    *mpibuf  = NULL;        /* buffer used to pack/unpack structures */
  int      mpibuf_size = 0;       /* size of the mpibuf                    */
  ESL_SQ_BLOCK    *sq_block  = NULL;
  WORKER_INFO     info;           /* the worker info */

  MPI_Status       mpistatus;        /* the mpi status */
  int              buf_size, n, pos; /* for packing/unpacking an MPI buffer */
  CM_TOPHITS *th = NULL;

  if ((status = init_shared_cfg(go, cfg, errbuf)) != eslOK) mpi_failure(errbuf);

  /* <cfg->abc> is not known 'til first CM is read. Could be DNA or RNA*/
  qhstatus = cm_file_Read(cfg->cmfp, TRUE, &(cfg->abc), &cm);

  if (qhstatus == eslOK) {
    /* One-time initializations after alphabet <abc> becomes known */
    if((status = CreateGenomicHMM(cfg->abc, errbuf, &(cfg->ghmm_sA), &(cfg->ghmm_tAA), &(cfg->ghmm_eAA), &cfg->ghmm_nstates)) != eslOK) cm_Fail("unable to create generative HMM\n%s", errbuf);
    if((status = set_dnull(cfg, cm, errbuf)) != eslOK) cm_Fail("unable to create set_dnull\n%s\n", errbuf);

    info.cm     = NULL;
    info.scA    = NULL;
    info.nhits  = 0;
    info.cutoff = cfg->sc_cutoff;
  }

  /* Outer loop: over each query CM in <cmfile>. */
  while (qhstatus == eslOK) {
    /* clone the non-configured CM we just read, we'll come back to it when we switch from global to local */
    if((status = cm_Clone(cm, errbuf, &nc_cm)) != eslOK) cm_Fail("unable to clone CM");

    cfg->ncm++;
    if(cfg->ncm == cfg->cmalloc) { /* expand our memory */
      cfg->cmalloc  += 128;
      ESL_REALLOC(cfg->expAA, sizeof(ExpInfo_t **) * cfg->cmalloc);
    }
    cmi = cfg->ncm-1;

    if((status = initialize_cm(go, cfg, errbuf, cm, FALSE)) != eslOK) mpi_failure(errbuf);

    /* get the sequences to search for all modes for this CM (we'll only search some of them, the ones the master tells us to) */
    if((status = generate_sequences(go, cfg, errbuf, cm, &sq_block)) != eslOK) mpi_failure(errbuf);

    for(exp_mode = 0; exp_mode < EXP_NMODES; exp_mode++) {
      /* do we need to switch from global configuration to local? */
      if(exp_mode > 0 && (! ExpModeIsLocal(exp_mode-1)) && ExpModeIsLocal(exp_mode)) {
	/* switch from global to local by copying the current exptail stats from <cm>
	 * into <nc_cm> and then configure <nc_cm> for local mode. We do it this
	 * way because as a rule we don't allow reconfiguration of CMs (to limit
	 * execution paths through configuration functions)
	 */
	FreeCM(cm);
	cm = nc_cm;
	if((status = initialize_cm(go, cfg, errbuf, cm, TRUE)) != eslOK) cm_Fail(errbuf);
      }
      /* update search info for round 0 (final round) for exp tail mode */
      if(ExpModeIsInside(exp_mode)) cm->search_opts |= CM_SEARCH_INSIDE;
      else                          cm->search_opts &= ~CM_SEARCH_INSIDE;
      
      if((status = cm_Clone(cm, errbuf, &(info.cm))) != eslOK) cm_Fail(errbuf);

      /* inform the master that we're ready for our first index */
      status = 0;
      MPI_Send(&status, 1, MPI_INT, 0, INFERNAL_READY_TAG, MPI_COMM_WORLD);

#if DEBUGMPI
      printf("worker %d waiting for 1st index from master\n", cfg->my_rank);
#endif 
      /* receive a sequence index from the master */
      if((status = MPI_Recv(&i, 1, MPI_INT, 0, INFERNAL_SEQIDX_TAG, MPI_COMM_WORLD, &mpistatus)) != 0) mpi_failure("failed to receive sequence index to search");

      while(i != -1) { 
#if DEBUGMPI
	printf("worker %d received i: %d from master\n", cfg->my_rank, i);
#endif 

	if((status = process_search_workunit(info.cm, errbuf, sq_block->list[i].dsq, sq_block->list[i].L, info.cutoff, &th)) != eslOK) mpi_failure("(worker) %s", errbuf);
#if DEBUGMPI
	printf("worker %d finished with i: %d\n", cfg->my_rank, i);
#endif 
	/* append copy of hit scores to scA */
	if(th->N > 0) { 
	  ESL_REALLOC(info.scA, sizeof(float) * (info.nhits + th->N)); /* this works even if info.scA is NULL */
	  for(h = 0; h < th->N; h++) info.scA[(info.nhits+h)] = th->unsrt[h].score;
	  info.nhits += th->N;
	}
#if DEBUGMPI
	printf("worker %d just added %" PRId64" hits to scA, now %" PRId64 " hits\n", cfg->my_rank, th->N, info.nhits);
#endif 
	cm_tophits_Destroy(th);

	/* inform the master we need another sequence index */
	status = 0;
	MPI_Send(&status, 1, MPI_INT, 0, INFERNAL_READY_TAG, MPI_COMM_WORLD);
	
#if DEBUGMPI
	printf("worker %d sent ready tag back to master\n", cfg->my_rank);
#endif 

	/* wait for the next sequence index */
	if((status = MPI_Recv(&i, 1, MPI_INT, 0, INFERNAL_SEQIDX_TAG, MPI_COMM_WORLD, &mpistatus)) != 0) mpi_failure("failed to receive sequence index to search");
      }

#if DEBUGMPI
      printf("worker %d received -1 seq idx from master, about to pack up hits\n", cfg->my_rank);
#endif 

      /* done with this CM/mode pair, pack up results */
      buf_size = 0;
      status = MPI_Pack_size(1, MPI_LONG_LONG_INT, MPI_COMM_WORLD, &n); buf_size += n; if (status != 0) mpi_failure("pack size failed");
      for(h = 0; h < info.nhits; h++) { 
	status = MPI_Pack_size(1, MPI_FLOAT, MPI_COMM_WORLD, &n); buf_size += n; if (status != 0) mpi_failure("pack size failed");
      }
      if (mpibuf == NULL || buf_size > mpibuf_size) {
	ESL_REALLOC(mpibuf, sizeof(char) * buf_size);
	mpibuf_size = buf_size;
      }
      pos = 0;
      if (MPI_Pack(&(info.nhits),    1, MPI_LONG_LONG_INT, mpibuf, buf_size, &pos, MPI_COMM_WORLD) != 0) mpi_failure("mpi pack failed"); 
      for(h = 0; h < info.nhits; h++) { 
	if (MPI_Pack(&(info.scA[h]), 1, MPI_FLOAT,         mpibuf, buf_size, &pos, MPI_COMM_WORLD) != 0) mpi_failure("mpi pack failed"); 
      }
#if DEBUGMPI
      printf("worker %d about to send %d byte-sized buf back to master\n", cfg->my_rank, n);
#endif 

      /* send results back to master */
      if (MPI_Send(mpibuf, buf_size, MPI_PACKED, 0, INFERNAL_SCORES_TAG, MPI_COMM_WORLD)    != 0) mpi_failure("mpi send failed");

      if(info.scA != NULL) { free(info.scA); info.scA = NULL; }
      info.nhits = 0;
    } /* end of for loop over exp_mode */
    FreeCM(cm);
    if(sq_block != NULL) esl_sq_DestroyBlock(sq_block);
    /* read the next CM */
    qhstatus = cm_file_Read(cfg->cmfp, TRUE, &(cfg->abc), &cm);
  }
  if(qhstatus != eslEOF) mpi_failure(cfg->cmfp->errbuf);

  return eslOK;

 ERROR:
  mpi_failure("Memory allocation error.");
  return eslEMEM;
}
#endif /*HAVE_MPI*/

static void
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_cmfile)
{
  ESL_GETOPTS *go = NULL;

  if ((go = esl_getopts_Create(options))     == NULL)     cm_Fail("Internal failure creating options object");
  if (esl_opt_ProcessEnvironment(go)         != eslOK)  { printf("Failed to process environment: %s\n", go->errbuf); goto ERROR; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { printf("Failed to parse command line: %s\n", go->errbuf); goto ERROR; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { printf("Failed to parse command line: %s\n", go->errbuf); goto ERROR; }
 
  /* help format: */
  if (esl_opt_GetBoolean(go, "--devhelp") == TRUE) {
    cm_banner(stdout, argv[0], banner);
    esl_usage(stdout, argv[0], usage);
    puts("\nStandard options:");
    esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 120=textwidth*/
    puts("\nUndocumented developer options:");
    esl_opt_DisplayHelp(stdout, go, 101, 2, 80); /* 1= group; 2 = indentation; 120=textwidth*/
  }    
  if (esl_opt_GetBoolean(go, "-h") == TRUE) {
    cm_banner(stdout, argv[0], banner);
    esl_usage(stdout, argv[0], usage);
    puts("\nOptions:");
    esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 120=textwidth*/
  }    

  if (esl_opt_ArgNumber(go)                  != 1)     { puts("Incorrect number of command line arguments.");      goto ERROR; }
  if ((*ret_cmfile = esl_opt_GetArg(go, 1))  == NULL)  { puts("Failed to get <cmfile> argument on command line"); goto ERROR; }
  
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
output_header(FILE *ofp, const ESL_GETOPTS *go, char *cmfile)
{
  cm_banner(ofp, go->argv[0], banner);
  
                                               fprintf(ofp, "# CM file:                                     %s\n", cmfile);
  if (esl_opt_IsUsed(go, "--forecast"))  {     fprintf(ofp, "# forecast-mode:                               on (%d CPUs)\n",  esl_opt_GetInteger(go, "--forecast")); }
  if (esl_opt_IsUsed(go, "--seed"))      {
    if (esl_opt_GetInteger(go, "--seed") == 0) fprintf(ofp, "# random number seed:                          one-time arbitrary\n");
    else                                       fprintf(ofp, "# random number seed set to:                   %d\n", esl_opt_GetInteger(go, "--seed"));
  }
  if (esl_opt_IsUsed(go, "-L"))          {     fprintf(ofp, "# total sequence length to search per mode:    %g Mb\n", esl_opt_GetReal(go, "-L")); }
  if (esl_opt_IsUsed(go, "--gtailn"))    {     fprintf(ofp, "# number of hits/Mb to fit (glocal):           %d\n", esl_opt_GetInteger(go, "--gtailn")); }
  if (esl_opt_IsUsed(go, "--ltailn"))    {     fprintf(ofp, "# number of hits/Mb to fit (local):            %d\n", esl_opt_GetInteger(go, "--gtailn")); }
  if (esl_opt_IsUsed(go, "--tailp"))     {     fprintf(ofp, "# fraction of histogram tail to fit:           %g\n", esl_opt_GetReal(go, "--tailp")); }
  if (esl_opt_IsUsed(go, "--beta"))      {     fprintf(ofp, "# tail loss probability for QDBs:              %g\n", esl_opt_GetReal(go, "--beta")); }
  if (esl_opt_IsUsed(go, "--nonbanded")) {     fprintf(ofp, "# query dependent bands (QDBs):                off\n"); }
  if (esl_opt_IsUsed(go, "--hfile"))     {     fprintf(ofp, "# saving fitted score histograms to file:      %s\n", esl_opt_GetString(go, "--hfile")); }
  if (esl_opt_IsUsed(go, "--sfile"))     {     fprintf(ofp, "# saving survival plot to file:                %s\n", esl_opt_GetString(go, "--sfile")); }
  if (esl_opt_IsUsed(go, "--qqfile"))    {     fprintf(ofp, "# saving Q-Q plot for histograms to file:      %s\n", esl_opt_GetString(go, "--qqfile")); }
  if (esl_opt_IsUsed(go, "--ffile"))     {     fprintf(ofp, "# saving lambdas for tail fit probs to file:   %s\n", esl_opt_GetString(go, "--ffile")); }
  if (esl_opt_IsUsed(go, "-v"))          {     fprintf(ofp, "# verbose output mode:                         on\n"); }
  if (esl_opt_IsUsed(go, "-T"))          {     fprintf(ofp, "# bit score threshold for histograms:          %g\n", esl_opt_GetReal(go, "-T")); }
  if (esl_opt_IsUsed(go, "--random"))    {     fprintf(ofp, "# generating seqs with cm->null (usually iid): yes\n"); }
  if (esl_opt_IsUsed(go, "--gc"))        {     fprintf(ofp, "# nucleotide distribution from file:           %s\n", esl_opt_GetString(go, "--gc")); }
  if (esl_opt_IsUsed(go, "--nonull3"))   {     fprintf(ofp, "# null3 bias corrections:                      off\n"); }
#ifdef HMMER_THREADS
  if (esl_opt_IsUsed(go, "--cpu"))       {     fprintf(ofp, "# number of worker threads:                    %d\n", esl_opt_GetInteger(go, "--cpu")); }
#endif
#ifdef HAVE_MPI
  if (esl_opt_IsUsed(go, "--mpi"))       {     fprintf(ofp, "# MPI:                                         on\n"); }
#endif
  fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
  return eslOK;
}

int
main(int argc, char **argv)
{
  int              status   = eslOK;

  ESL_GETOPTS     *go  = NULL;    /* command line processing                 */
  struct cfg_s     cfg;           /* configuration data                      */
  char             errbuf[cmERRBUFSIZE];
  int              i;
  ESL_STOPWATCH   *w  = esl_stopwatch_Create();
  if(w == NULL) cm_Fail("Memory allocation error, stopwatch could not be created.");
  esl_stopwatch_Start(w);
  

  /* Set processor specific flags */
  impl_Init();

  /* Initialize what we can in the config structure (without knowing the alphabet yet)
   */
  cfg.cmfile       = NULL;
  cfg.do_mpi       = FALSE;               /* this gets reset below, if we init MPI */
  cfg.nproc        = 0;                   /* this gets reset below, if we init MPI */
  cfg.my_rank      = 0;                   /* this gets reset below, if we init MPI */
  cfg.r            = NULL; 
  cfg.abc          = NULL; 
  cfg.w_stage      = NULL; 
  cfg.be_verbose   = FALSE;
  cfg.cmalloc      = 128;
  cfg.tmpfile      = NULL;
  cfg.mode         = NULL;
  cfg.dnull        = NULL;
  cfg.ghmm_nstates = 0;
  cfg.ghmm_sA      = NULL;
  cfg.ghmm_tAA     = NULL;
  cfg.ghmm_eAA     = NULL;
  cfg.ncm          = 0;
  cfg.ccom         = NULL; /* created in get_cmcalibrate_comlog_info() for masters, stays NULL in workers */
  cfg.cdate        = NULL; /* created in get_cmcalibrate_comlog_info() for masters, stays NULL in workers */
  cfg.L            = EXPTAIL_CHUNKLEN; /* 10 Kb chunks are searched */
  cfg.N            = 0;    /* gets set below, after go is setup */
  cfg.cmfp         = NULL; /* remains NULL for mpi workers */
  cfg.exphfp       = NULL; /* remains NULL for mpi workers */
  cfg.expsfp       = NULL; /* remains NULL for mpi workers */
  cfg.expqfp       = NULL; /* remains NULL for mpi workers */
  cfg.exptfitfp    = NULL; /* remains NULL for mpi workers */

  cfg.cmalloc      = 128;
  ESL_ALLOC(cfg.expAA, sizeof(ExpInfo_t **) * cfg.cmalloc); /* this will grow if needed */

  ESL_DASSERT1((EXP_CM_GC  == 0));
  ESL_DASSERT1((EXP_CM_GI  == 1));
  ESL_DASSERT1((EXP_CM_LC  == 2));
  ESL_DASSERT1((EXP_CM_LI  == 3));
  ESL_DASSERT1((EXP_NMODES == 4));

  /* Initializations */
  init_ilogsum();
  FLogsumInit();
  process_commandline(argc, argv, &go, &(cfg.cmfile));
  cfg.N = (int) (((esl_opt_GetReal(go, "-L") * 1000000.) / (float) cfg.L) + 0.5);
  ///cfg.N = 1;

  /* Figure out who we are, and send control there: 
   * we might be an MPI master, an MPI worker, or a serial program.
   */
#ifdef HAVE_MPI

  /* TEMP */
  pid_t pid;
  /* get the process id */
  pid = getpid();
  printf("The process id is %d\n", pid);
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

  if(! esl_opt_IsOn(go, "--forecast") && cfg.my_rank == 0) { /* master, serial or mpi */
    /* Rewind the CM file for a second pass.
     * Write a temporary CM file with new stats information in it
     */
    int   status;
    int   cmi;
    CM_t *cm;
    FILE *outfp;
    sigset_t blocksigs;  /* list of signals to protect from */

    cm_file_Position(cfg.cmfp, (off_t) 0);

    if (esl_FileExists(cfg.tmpfile))                    cm_Fail("Ouch. Temporary file %s appeared during the run.", cfg.tmpfile);
    if ((outfp = fopen(cfg.tmpfile, cfg.mode)) == NULL) cm_Fail("Ouch. Temporary file %s couldn't be opened for writing.", cfg.tmpfile); 
    
    for (cmi = 0; cmi < cfg.ncm; cmi++) {
      if ((status = cm_file_Read(cfg.cmfp, TRUE, &(cfg.abc), &cm)) != eslOK) cm_Fail("Ran out of CMs too early in pass 2");
      if (cm == NULL)                                                        cm_Fail("CM file %s was corrupted? Parse failed in pass 2", cfg.cmfile);

      /* update the cm->comlog info */
      if((status = update_comlog(go, errbuf, cfg.ccom, cfg.cdate, cm)) != eslOK) cm_Fail(errbuf);
	
      if(cm->expA != NULL) { 
	for(i = 0; i < EXP_NMODES; i++)  free(cm->expA[i]); free(cm->expA);
      }
      ESL_ALLOC(cm->expA, sizeof(ExpInfo_t *) * EXP_NMODES);

      cm->expA   = cfg.expAA[cmi];
      cm->flags |= CMH_EXPTAIL_STATS; 
      if(cfg.cmfp->is_binary) { 
	if ((status = cm_file_WriteBinary(outfp, -1, cm, NULL)) != eslOK) ESL_FAIL(status, errbuf, "binary CM save failed");
      }
      else { 
	if ((status = cm_file_WriteASCII(outfp, -1, cm)) != eslOK) ESL_FAIL(status, errbuf, "CM save failed");
      }
      FreeCM(cm);
    } /* end of from idx = 0 to ncm */
    
    /* Now, carefully remove original file and replace it
     * with the tmpfile. Note the protection from signals;
     * we wouldn't want a user to ctrl-C just as we've deleted
     * their CM file but before the new one is moved.
     */
    cm_file_Close(cfg.cmfp);
    if (fclose(outfp)   != 0)                            cm_Fail("system error during rewrite of CM file");
    if (sigemptyset(&blocksigs) != 0)                    cm_Fail("system error during rewrite of CM file.");
    if (sigaddset(&blocksigs, SIGINT) != 0)              cm_Fail("system error during rewrite of CM file.");
    if (sigprocmask(SIG_BLOCK, &blocksigs, NULL) != 0)   cm_Fail("system error during rewrite of CM file.");
    if (remove(cfg.cmfile) != 0)                         cm_Fail("system error during rewrite of CM file.");
    if (rename(cfg.tmpfile, cfg.cmfile) != 0)            cm_Fail("system error during rewrite of CM file.");
    if (sigprocmask(SIG_UNBLOCK, &blocksigs, NULL) != 0) cm_Fail("system error during rewrite of CM file.");
    free(cfg.tmpfile);
    
    /* master specific cleaning */
    if (cfg.exphfp   != NULL) { 
      fclose(cfg.exphfp);
      printf("# Histogram of high scoring hits in random seqs saved to file %s.\n", esl_opt_GetString(go, "--hfile"));
    }
    if (cfg.expsfp   != NULL) { 
      fclose(cfg.expsfp);
      printf("# Survival plot for exponential tails saved to file %s.\n", esl_opt_GetString(go, "--sfile"));
    }
    if (cfg.expqfp   != NULL) { 
      fclose(cfg.expqfp);
      printf("# Exponential tail QQ plots saved to file %s.\n", esl_opt_GetString(go, "--qqfile"));
    }
    if (cfg.exptfitfp   != NULL) { 
      fclose(cfg.exptfitfp);
      printf("# Exponential tail fit points saved to file %s.\n", esl_opt_GetString(go, "--ffile"));
    }
    if (cfg.ccom  != NULL) free(cfg.ccom);
    if (cfg.cdate != NULL) free(cfg.cdate);
    if (cfg.expAA != NULL) free(cfg.expAA);
  }
  esl_stopwatch_Stop(w);

  /* clean up */
  if (cfg.abc       != NULL) esl_alphabet_Destroy(cfg.abc);
  if (cfg.w_stage   != NULL) esl_stopwatch_Destroy(cfg.w_stage);
  if (cfg.r         != NULL) esl_randomness_Destroy(cfg.r);
  if (cfg.r_est     != NULL) esl_randomness_Destroy(cfg.r_est);
  if (cfg.my_rank == 0) { 
    printf("#\n");
    esl_stopwatch_Display(stdout, w, "# CPU time: ");
  }
  if(cfg.ghmm_eAA != NULL) { 
    for(i = 0; i < cfg.ghmm_nstates; i++) free(cfg.ghmm_eAA[i]); 
    free(cfg.ghmm_eAA);
  }
  if(cfg.ghmm_tAA != NULL) { 
    for(i = 0; i < cfg.ghmm_nstates; i++) free(cfg.ghmm_tAA[i]); 
    free(cfg.ghmm_tAA);
  }
  if(cfg.ghmm_sA != NULL) free(cfg.ghmm_sA);

  esl_stopwatch_Destroy(w);
  esl_getopts_Destroy(go);
  return eslOK;

 ERROR: 
  cm_Fail("Memory allocation error.");
  return status; /* NEVERREACHED */
}

/* init_master_cfg()
 * Called by masters, mpi or serial.
 * Allocates/sets: 
 *    cfg->cmfp        - open CM file                
 *    cfg->exphfp      - optional output file
 *    cfg->expsfp      - optional output file
 *    cfg->expqfp      - optional output file
 *    cfg->gc_freq     - observed GC freqs (if --gc invoked)
 *    cfg->expAA       - the exp tail stats, allocated only
 *    cfg->r           - source of randomness
 *    cfg->tmpfile     - temp file for rewriting cm file
 *    cfg->be_verbose  - print extra info? 
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

  /* Masters only initializations: */
  /* open optional output files */
  if (esl_opt_GetString(go, "--hfile") != NULL) {
    if ((cfg->exphfp = fopen(esl_opt_GetString(go, "--hfile"), "w")) == NULL)
      ESL_FAIL(eslFAIL, errbuf, "Failed to open exp tail histogram save file %s for writing\n", esl_opt_GetString(go, "--hfile"));
  }
  if (esl_opt_GetString(go, "--sfile") != NULL) { 
    if ((cfg->expsfp = fopen(esl_opt_GetString(go, "--sfile"), "w")) == NULL)
      ESL_FAIL(eslFAIL, errbuf, "Failed to open survival plot save file %s for writing\n", esl_opt_GetString(go, "--sfile"));
  }
  if (esl_opt_GetString(go, "--qqfile") != NULL) {
    if ((cfg->expqfp = fopen(esl_opt_GetString(go, "--qqfile"), "w")) == NULL)
      ESL_FAIL(eslFAIL, errbuf, "Failed to open exp tail QQ plot save file %s for writing\n", esl_opt_GetString(go, "--qqfile"));
  }
  if (esl_opt_GetString(go, "--ffile") != NULL) {
    if ((cfg->exptfitfp = fopen(esl_opt_GetString(go, "--ffile"), "w")) == NULL)
      ESL_FAIL(eslFAIL, errbuf, "Failed to open exp tail save file %s for writing\n", esl_opt_GetString(go, "--ffile"));
  }

  cfg->be_verbose = FALSE;
  if (esl_opt_GetBoolean(go, "-v")) cfg->be_verbose = TRUE;        
  
  /* create the stopwatch */
  cfg->w_stage = esl_stopwatch_Create();
  if(cfg->w_stage == NULL) ESL_FAIL(eslEMEM, errbuf, "Failed to create stopwatch.");

  /* From HMMER 2.4X hmmcalibrate.c:
   * Generate calibrated CM(s) in a tmp file in the current
   * directory. When we're finished, we delete the original
   * CM file and rename() this one. That way, the worst
   * effect of a catastrophic failure should be that we
   * leave a tmp file lying around, but the original CM
   * file remains uncorrupted. tmpnam() doesn't work portably here,
   * because it'll put the file in /tmp and we won't
   * necessarily be able to rename() it from there.
   */
  ESL_ALLOC(cfg->tmpfile, (sizeof(char) * (strlen(cfg->cmfile) + 5)));
  strcpy(cfg->tmpfile, cfg->cmfile);
  strcat(cfg->tmpfile, ".xxx");	/* could be more inventive here... */
  if (esl_FileExists(cfg->tmpfile)) ESL_FAIL(eslFAIL, errbuf, "temporary file %s already exists; please delete it first", cfg->tmpfile);
  if (cfg->cmfp->is_binary) cfg->mode = "wb";
  else                      cfg->mode = "w"; 

  /* fill cfg->ccom, and cfg->cdate */
  if((status = get_cmcalibrate_comlog_info(go, cfg, errbuf)) != eslOK) return status;

  return eslOK;

 ERROR:
  ESL_FAIL(eslEMEM, errbuf, "init_master_cfg(), memory allocation error."); 
  return eslEMEM; /* NEVER REACHED */
}


/* init_shared_cfg()
 * Called by serial masters and mpi workers and masters.
 * Allocates/sets: 
 *    cfg->cmfp        - open CM file                
 *    cfg->gc_freq     - observed GC freqs (if --gc invoked)
 *    cfg->r           - source of randomness
 *    cfg->sc_cutoff   - cutoff score for collecting hits 
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

  /* optionally, get distribution of GC content from an input database (default is use cm->null for GC distro) */
  if(esl_opt_GetString(go, "--gc") != NULL) {
    ESL_ALPHABET *tmp_abc = NULL;
    tmp_abc = esl_alphabet_Create(eslRNA);
    ESL_SQFILE   *dbfp;             
    status = esl_sqfile_Open(esl_opt_GetString(go, "--gc"), eslSQFILE_UNKNOWN, NULL, &dbfp);
    if (status == eslENOTFOUND)    ESL_FAIL(status, errbuf, "No such file: %s.", esl_opt_GetString(go, "--gc")); 
    else if (status == eslEFORMAT) ESL_FAIL(status, errbuf, "file: %s format unrecognized.", esl_opt_GetString(go, "--gc")); 
    else if (status != eslOK)      ESL_FAIL(status, errbuf, "Failed to open sequence database file %s, code %d.", esl_opt_GetString(go, "--gc"), status); 
    if((status = GetDBInfo(tmp_abc, dbfp, errbuf, NULL, NULL, &(cfg->gc_freq))) != eslOK) return status; 
    esl_vec_DNorm(cfg->gc_freq, GC_SEGMENTS);
    esl_alphabet_Destroy(tmp_abc);
    esl_sqfile_Close(dbfp); 
  }
  
  /* seed RNG */
  cfg->r     = esl_randomness_Create(esl_opt_GetInteger(go, "--seed"));
  cfg->r_est = esl_randomness_Create(esl_opt_GetInteger(go, "--seed"));

  if(cfg->r     == NULL) ESL_FAIL(eslEMEM, errbuf, "Failed to create RNG.");
  if(cfg->r_est == NULL) ESL_FAIL(eslEMEM, errbuf, "Failed to create RNG.");

  if(esl_opt_IsOn(go, "-T")) cfg->sc_cutoff = esl_opt_GetReal(go, "-T");
  else                       cfg->sc_cutoff = -eslINFINITY;

  return eslOK;
}


/* serial_master()
 * The serial version of cmcalibrate.
 * 
 * A master can only return if it's successful. All errors are handled immediately and fatally with cm_Fail().
 */
static void
serial_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int      status;                /* Easel status */
  int      qhstatus = eslOK;      /* status from cm_file_Read() */
  char     errbuf[cmERRBUFSIZE];  /* for printing error messages */
  CM_t    *cm    = NULL;          /* the CM */
  CM_t    *nc_cm = NULL;          /* a non-configured copy of the CM */
  int      i, h;                  /* counters */
  int      cmi;                   /* CM index, which model we're working on */
  char     time_buf[128];	  /* string for printing elapsed time (safely holds up to 10^14 years) */
  int      exp_mode;              /* ctr over exp tail modes */
  double   cm_psec;               /* predicted number of seconds for calibrating current CM */
  double   cm_asec;               /* predicted number of seconds for calibrating current CM */
  double   total_psec = 0.;       /* predicted number of seconds for calibrating all CMs */
  double   total_asec = 0.;       /* predicted number of seconds for calibrating all CMs */
  double   asecA[EXP_NMODES];     /* stores actual timings for each exp tail fit stage, for each CM */
  double   psecA[EXP_NMODES];     /* stores predicted timings for each exp tail fit stage, for each CM */
  double   psec;                  /* predicted seconds */
  int      ncpus        = 0;      /* number of CPUs working */
  int64_t  merged_nhits = 0;      /* number of hits reported thus far, for all seqs */
  float   *merged_scA = NULL;     /* [0..merged_nhits-1] hit scores for all seqs */
  double   tmp_mu, tmp_lambda;    /* temporary mu and lambda used for setting exp tails */
  int      tmp_nrandhits;         /* temporary number of rand hits found */
  float    tmp_tailp;             /* temporary tail mass probability fit to an exponential */

  /* variables needed for threaded implementation */
  WORKER_INFO     *info      = NULL; /* the worker info */
  int              infocnt   = 0;    /* number of worker infos */
#ifdef HMMER_THREADS
  ESL_SQ_BLOCK    *sq_block  = NULL;
  ESL_THREADS     *threadObj = NULL;
  ESL_WORK_QUEUE  *queue     = NULL;
#endif
  
  if ((status = init_master_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);
  if ((status = print_run_info (go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);

#ifdef HMMER_THREADS
  /* initialize thread data */
  if      (esl_opt_IsUsed(go, "--forecast")) ncpus = 0;
  else if (esl_opt_IsOn  (go, "--cpu"))      ncpus = esl_opt_GetInteger(go, "--cpu");
  else                                       esl_threads_CPUCount(&ncpus);
  printf("NCPUS: %d\n", ncpus);
  if (ncpus > 0) {
      threadObj = esl_threads_Create(&pipeline_thread);
      queue = esl_workqueue_Create(ncpus * 2);
  }
#endif

  infocnt = (ncpus == 0) ? 1 : ncpus;
  ESL_ALLOC(info, sizeof(WORKER_INFO) * infocnt);

  /* <cfg->abc> is not known 'til first CM is read. Could be DNA or RNA*/
  qhstatus = cm_file_Read(cfg->cmfp, TRUE, &(cfg->abc), &cm);

  if (qhstatus == eslOK) {
    /* One-time initializations after alphabet <abc> becomes known */
    output_header(stdout, go, cfg->cmfile);

    if((status = CreateGenomicHMM(cfg->abc, errbuf, &(cfg->ghmm_sA), &(cfg->ghmm_tAA), &(cfg->ghmm_eAA), &cfg->ghmm_nstates)) != eslOK) cm_Fail("unable to create generative HMM\n%s", errbuf);
    if((status = set_dnull(cfg, cm, errbuf)) != eslOK) cm_Fail("unable to create set_dnull\n%s\n", errbuf);
    
    for (i = 0; i < infocnt; ++i)    {
      info[i].cm     = NULL;
      info[i].scA    = NULL;
      info[i].nhits  = 0;
      info[i].cutoff = cfg->sc_cutoff;
#ifdef HMMER_THREADS
      info[i].queue      = queue;
#endif
    }

#ifdef HMMER_THREADS    
    for (i = 0; i < ncpus * 2; ++i) {
      if((sq_block = esl_sq_CreateDigitalBlock(BLOCK_SIZE, cfg->abc)) == NULL)  cm_Fail("Failed to allocate sequence block");
      if((status   = esl_workqueue_Init(queue, sq_block))             != eslOK) cm_Fail("Failed to add block to work queue");
    }
#endif
  }
  
  /* Outer loop: over each query CM in <cmfile>. */
  while (qhstatus == eslOK) {
    /* clone the non-configured CM we just read, we'll come back to it when we switch from global to local */
    if((status = cm_Clone(cm, errbuf, &nc_cm)) != eslOK) cm_Fail("unable to clone CM");

    cfg->ncm++;
    if(cfg->ncm == cfg->cmalloc) { /* expand our memory */
      cfg->cmalloc  += 128;
      ESL_REALLOC(cfg->expAA, sizeof(ExpInfo_t **) * cfg->cmalloc);
    }
    cmi = cfg->ncm-1;
    cm_asec = cm_psec = 0.;

    if((status = initialize_cm(go, cfg, errbuf, cm, FALSE))          != eslOK) cm_Fail(errbuf);
    if((status = initialize_stats(go, cfg, errbuf))                  != eslOK) cm_Fail(errbuf);
    if((status = print_per_cm_column_headings(go, cfg, errbuf, cm))  != eslOK) cm_Fail(errbuf);

    /* get the sequences to search */
    if((status = generate_sequences(go, cfg, errbuf, cm, &sq_block)) != eslOK) cm_Fail(errbuf);
      
    for(exp_mode = 0; exp_mode < EXP_NMODES; exp_mode++) {
      /* do we need to switch from global configuration to local? */
      if(exp_mode > 0 && (! ExpModeIsLocal(exp_mode-1)) && ExpModeIsLocal(exp_mode)) {
	/* switch from global to local by copying the current exptail stats from <cm>
	 * into <nc_cm> and then configure <nc_cm> for local mode. We do it this
	 * way because as a rule we don't allow reconfiguration of CMs (to limit
	 * execution paths through configuration functions)
	 */
	FreeCM(cm);
	cm = nc_cm;
	if((status = initialize_cm(go, cfg, errbuf, cm, TRUE)) != eslOK) cm_Fail(errbuf);
      }
      /* update search info for round 0 (final round) for exp tail mode */
      if(ExpModeIsInside(exp_mode)) cm->search_opts |= CM_SEARCH_INSIDE;
      else                          cm->search_opts &= ~CM_SEARCH_INSIDE;

      /* estimate time for this round */
      if((status = estimate_time_for_exp_round(go, cfg, errbuf, cm, exp_mode, &psec)) != eslOK) cm_Fail(errbuf); 
      psec *= cfg->N * cfg->L; /* psec was per residue */
      /* update psec due to parallelization */
      if(esl_opt_IsUsed(go, "--forecast")) { 
	psec /= ESL_MAX(1, (esl_opt_GetInteger(go, "--forecast") - 1)); 
	/* Note, this is correct for MPI but will be slightly off for
	 * threaded - because all threads do searches, whereas in MPI
	 * master doesn't do searches.
	 */
      }
      else if(ncpus > 0) { 
	psec /= ncpus;
      }
      
      psecA[exp_mode] = psec;
      cm_psec    += psec;
      total_psec += psec;
      print_exp_line(go, cfg, errbuf, exp_mode, cfg->N, cfg->L, psec);
      if(esl_opt_IsUsed(go, "--forecast")) continue; /* special mode, we don't do the calibration, just print the predicting timings */

      esl_stopwatch_Start(cfg->w_stage);
      fflush(stdout);
	
      /* clone CM for each thread */
      for (i = 0; i < infocnt; ++i) {
	if((status = cm_Clone(cm, errbuf, &(info[i].cm))) != eslOK) cm_Fail(errbuf);
	if(esl_opt_GetBoolean(go, "-x")) info[i].cm = cm;
	info[i].scA    = NULL;
	info[i].nhits  = 0;
	info[i].cutoff = cfg->sc_cutoff;

#ifdef HMMER_THREADS
	if (ncpus > 0) esl_threads_AddThread(threadObj, &info[i]);
#endif
      }

#ifdef HMMER_THREADS
      if (ncpus > 0)  status = thread_loop(info, errbuf, threadObj, queue, sq_block);
      else            status = serial_loop(info, errbuf, sq_block);
#else
      status = serial_loop(info, errbuf, sq_block);
#endif
      if(status != eslOK) cm_Fail(errbuf);

      merged_nhits = 0;
      if(merged_scA != NULL) { free(merged_scA); merged_scA = NULL; }
      for (i = 0; i < infocnt; ++i) {
	if(info[i].nhits > 0) { 
	  ESL_REALLOC(merged_scA, sizeof(float) * (merged_nhits + info[i].nhits)); /* this works even if merged_scA == NULL */
	  for(h = 0; h < info[i].nhits; h++) merged_scA[(merged_nhits+h)] = info[i].scA[h];
	  merged_nhits += info[i].nhits;
	}
      }

      if(cfg->exptfitfp != NULL) { 
	fprintf(cfg->exptfitfp, "# CM: %s\n", cm->name);
	fprintf(cfg->exptfitfp, "# mode: %12s\n", DescribeExpMode(exp_mode));
      }
      printf("about to call fit_histogram, %" PRId64 " hits\n", merged_nhits);
      if((status = fit_histogram(go, cfg, errbuf, merged_scA, merged_nhits, exp_mode, &tmp_mu, &tmp_lambda, &tmp_nrandhits, &tmp_tailp)) != eslOK) cm_Fail(errbuf);
      SetExpInfo(cfg->expAA[cmi][exp_mode], tmp_lambda, tmp_mu, (long) (cfg->L * cfg->N), tmp_nrandhits, tmp_tailp);
	
      esl_stopwatch_Stop(cfg->w_stage);
      asecA[exp_mode] = cfg->w_stage->elapsed;
      cm_asec    += cfg->w_stage->elapsed;
      total_asec += cfg->w_stage->elapsed;
      FormatTimeString(time_buf, cfg->w_stage->elapsed, FALSE);
      printf("  %10s\n", time_buf);
      fflush(stdout);
    } /* end of for(exp_mode = 0; exp_mode < EXP_NMODES; exp_mode++) */

    if(cfg->be_verbose) if((status = debug_print_expinfo_array(cm, errbuf, cfg->expAA[cmi])) != eslOK) cm_Fail(errbuf);
    print_per_cm_summary(go, cfg, errbuf, cm, cm_psec, cm_asec);
    if(!esl_opt_IsOn(go, "--forecast")) { if((status = print_post_calibration_info(go, cfg, errbuf, stdout, cm, psecA, asecA)) != eslOK) cm_Fail(errbuf); }
    FreeCM(cm);
    if(sq_block != NULL) esl_sq_DestroyBlock(sq_block);

    printf("//\n");
    fflush(stdout);

    /* read the next CM */
    qhstatus = cm_file_Read(cfg->cmfp, TRUE, &(cfg->abc), &cm);
  } /* end of while(qhstatus == eslOK) */
  if(qhstatus != eslEOF) cm_Fail(cfg->cmfp->errbuf);
  
  if(cfg->ncm > 1 && (esl_opt_IsOn(go, "--forecast"))) { 
    fprintf(stdout, "#\n");
    FormatTimeString(time_buf, total_psec, FALSE);
    fprintf(stdout, "# total predicted time for all %d CMs: %s\n", cfg->ncm, time_buf);
    fprintf(stdout, "#\n");
  }
  return;
      
 ERROR:
  cm_Fail("Memory allocation error.");
  return;
}

/* serial_loop(): 
 * 
 * Search each sequence and collect scores of hits. 
 */
static int
serial_loop(WORKER_INFO *info, char *errbuf, ESL_SQ_BLOCK *sq_block)
{
  int status;
  CM_TOPHITS *th = NULL;
  int  i, h;             /* counters */

  /* reinitialize array of hit scores */
  info->nhits = 0;
  if(info->scA != NULL) free(info->scA);

  for(i = 0; i < sq_block->count; i++) { 
    if((status = process_search_workunit(info->cm, errbuf, sq_block->list[i].dsq, sq_block->list[i].L, info->cutoff, &th)) != eslOK) cm_Fail(errbuf);
    ///int zz; for(zz = 0; zz < th->N; zz++) {
    ///printf("hit %3d  %4" PRId64 "..%4" PRId64 "  (%4" PRId64 ")  %6.2f bits\n", zz, th->unsrt[zz].start, th->unsrt[zz].stop, th->unsrt[zz].stop - th->unsrt[zz].start + 1, th->unsrt[zz].score);
    ///}
    ///if(info->cm->search_opts & CM_SEARCH_INSIDE) cm_Fail("done");
    /* append copy of hit scores to scA */
    if(th->N > 0) { 
      ESL_REALLOC(info->scA, sizeof(float) * (info->nhits + th->N)); /* this works even if info->scA == NULL */
      for(h = 0; h < th->N; h++) info->scA[(info->nhits+h)] = th->unsrt[h].score;
      info->nhits += th->N;
    }
    cm_tophits_Destroy(th);
  }
  return eslOK;
  
 ERROR: 
  ESL_FAIL(status, errbuf, "out of memory");
  return status; /* NEVERREACHED */
}
 
#ifdef HMMER_THREADS
static int
thread_loop(WORKER_INFO *info, char *errbuf, ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQ_BLOCK *sq_block)
{
  int     status   = eslOK;
  ESL_SQ *sq       = NULL;
  void   *new_sq   = NULL;
  ESL_SQ *empty_sq = NULL;
  int     nworkers = esl_threads_GetWorkerCount(obj);
  int     i, j;       

  esl_workqueue_Reset(queue);
  esl_threads_WaitForStart(obj);

  status = esl_workqueue_ReaderUpdate(queue, NULL, &new_sq);
  if (status != eslOK) cm_Fail("Work queue reader failed");

  /* Main loop: */
  for(i = 0; i < sq_block->count; i++) { 
    sq = (ESL_SQ *) new_sq;
    sq = sq_block->list + i;
    status = esl_workqueue_ReaderUpdate(queue, sq, &new_sq);
    if (status != eslOK) cm_Fail("Work queue reader failed");
  }

  /* now send a empty sq to all workers signaling them to stop */
  empty_sq = esl_sq_Create();
  for(j = 0; j < nworkers; j++) { 
    status = esl_workqueue_ReaderUpdate(queue, empty_sq, &new_sq);
    if (status != eslOK) cm_Fail("Work queue reader failed");
  }

  status = esl_workqueue_ReaderUpdate(queue, sq, NULL);

  /* wait for all the threads to complete */
  esl_threads_WaitForFinish(obj);
  esl_workqueue_Complete(queue);  

  esl_sq_Destroy(empty_sq);
  return status;
}

/* pipeline_thread()
 * 
 * Receive a sequence(s) from the master, search it
 * with the CM, and add scores of hits to info->scA.
 */

static void 
pipeline_thread(void *arg)
{
  int          status;
  int          workeridx;
  WORKER_INFO *info;
  ESL_THREADS *obj;
  ESL_SQ      *sq = NULL;
  void        *new_sq;
  int          h;
  CM_TOPHITS  *th = NULL;
  char         errbuf[cmERRBUFSIZE];
  
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

  status = esl_workqueue_WorkerUpdate(info->queue, NULL, &new_sq);
  if (status != eslOK) cm_Fail("Work queue worker failed");

  /* loop until all sequences have been processed */
  sq = (ESL_SQ *) new_sq;
  while (sq->L != -1) { 
    if((status = process_search_workunit(info->cm, errbuf, sq->dsq, sq->L, info->cutoff, &th)) != eslOK) cm_Fail(errbuf);
    /* append copy of hit scores to scA */
    if(th->N > 0) { 
      ESL_REALLOC(info->scA, sizeof(float) * (info->nhits + th->N)); /* this works even if info->scA == NULL */
      for(h = 0; h < th->N; h++) info->scA[(info->nhits+h)] = th->unsrt[h].score;
      info->nhits += th->N;
    }
    cm_tophits_Destroy(th);

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


/* initialize_cm()
 * Setup the CM based on the command-line options/defaults;
 * only set flags and a few parameters. ConfigCM() configures
 * the CM.
 */
static int
initialize_cm(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, int do_local)
{
  int status;

  /* config QDB? */
  if(esl_opt_GetBoolean(go, "--nonbanded")) { 
    cm->search_opts |= CM_SEARCH_NONBANDED; /* don't use QDB to search */
    /* no need to recalculate QDBs, don't raise CM_CONFIG_QDB */
  }
  else {
    cm->search_opts |= CM_SEARCH_QDB; /* use QDB to search */
    /* check if cm->qdbinfo->beta2 == <x> from --beta, if so we don't need to recalculate QDBs */
    if(CheckCMQDBInfo(cm->qdbinfo, 0., FALSE, esl_opt_GetReal(go, "--beta"), TRUE) != eslOK) {
      /* we'll use beta2 for calibration, setting them both as equal makes it slightly more efficient */
      cm->config_opts |= CM_CONFIG_QDB;   /* configure QDB */
      cm->qdbinfo->beta1 = esl_opt_GetReal(go, "--beta"); 
      cm->qdbinfo->beta2 = esl_opt_GetReal(go, "--beta"); 
    }
  }

  cm->search_opts |= CM_SEARCH_NOALIGN;

  if(! esl_opt_GetBoolean(go, "--nonull3")) cm->search_opts |= CM_SEARCH_NULL3;

  /* ALWAYS use the greedy overlap resolution algorithm to return hits for exp calculation
   * it's irrelevant for filter threshold stats, we return best score per seq for that */
  /* don't turn on CM_SEARCH_CMNOTGREEDY */

  if(do_local) { 
    cm->config_opts |= CM_CONFIG_LOCAL;
    cm->config_opts |= CM_CONFIG_HMMLOCAL;
    cm->config_opts |= CM_CONFIG_HMMEL;
  }
  /* we'll need a scan matrix too */
  cm->config_opts |= CM_CONFIG_SCANMX;

  /* configure */
  if((status = cm_Configure(cm, errbuf, -1)) != eslOK) return status; 

  if(cm->smx == NULL) ESL_FAIL(eslEINVAL, errbuf, "unable to create scan matrix for CM");
  return eslOK;
}

/* initialize_stats()
 * Allocate and initialize cfg->expAA */
static int
initialize_stats(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  int status;
  int i;
  int cmi = cfg->ncm-1;

  ESL_DPRINTF1(("initializing cmstats\n"));

  ESL_ALLOC(cfg->expAA[cmi], sizeof(ExpInfo_t *) * EXP_NMODES);
  for(i = 0; i < EXP_NMODES; i++) cfg->expAA[cmi][i] = CreateExpInfo();
  
  return eslOK;

 ERROR: 
  return status;
}


/* fit_histogram()
 * Create, fill and fit the tail of a histogram to an exponential tail. Data to fill the histogram
 * is given as <scores>.
 */
static int
fit_histogram(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, float *scores, int nscores, int exp_mode, double *ret_mu, double *ret_lambda, int *ret_nrandhits, float *ret_tailp)
{
  int status;
  double mu;
  double lambda;
  int i;
  double *xv;         /* raw data from histogram */
  int     n,z;  
  float tailp;
  double  params[2];
  int     nrandhits; 
  float   a;
  float   nhits_to_fit;

  ESL_HISTOGRAM *h = NULL;       /* histogram of scores */

  /* Initialize histogram; these numbers are guesses */
  if((h = esl_histogram_CreateFull(-100., 100., .1)) == NULL) return eslEMEM;    

  /* fill histogram */
  for(i = 0; i < nscores; i++) {
    if((status = esl_histogram_Add(h, scores[i])) != eslOK) ESL_FAIL(status, errbuf, "fit_histogram(), esl_histogram_Add() call returned non-OK status: %d\n", status);
    /* printf("%4d %.3f\n", i, scores[i]); */
  }

  /* fit scores to an exponential tail */
  if(cfg->exptfitfp != NULL) { 
    /* fit to 41 different tailp values and print lambda, mu to a save file*/
    fprintf(cfg->exptfitfp, "# %11s  %10s  %10s  %12s\n", "tail pmass",  "lambda",     "mu",         "nhits");
    fprintf(cfg->exptfitfp, "# %11s  %10s  %10s  %12s\n", "-----------", "----------", "----------", "------------");
    for(a = 0.; a >= -4.; a -= 0.1) { 
      tailp = pow(10., a);
      esl_histogram_GetTailByMass(h, tailp, &xv, &n, &z); 
      if(n > 1) { 
	esl_exp_FitComplete(xv, n, &(params[0]), &(params[1]));
	esl_histogram_SetExpectedTail(h, params[0], tailp, &esl_exp_generic_cdf, &params);
	fprintf(cfg->exptfitfp, "  %.9f  %10.6f  %10.4f  %12d\n", tailp, params[1], params[0], n);
      }
      else { 
	fprintf(cfg->exptfitfp, "  %.9f  %10s  %10s  %12d\n", tailp, "N/A", "N/A", n);
      }
    }
    fprintf(cfg->exptfitfp, "//\n");
  }
  /* end of if cfg->exptfitfp != NULL) */

  /* determine the fraction of the tail to fit, if --tail-p, it's easy */
  if(esl_opt_IsOn(go, "--tailp")) { 
    tailp = esl_opt_GetReal(go, "--tailp");
    tailp = ESL_MIN(tailp, ((float) esl_opt_GetInteger(go, "--tailxn") / (float) h->n)); /* ensure we don't exceed our max nhits in tail */
  }
  else { /* number of hits is per Mb and specific to local or glocal fits */
    if(ExpModeIsLocal(exp_mode)) { /* local mode */
      nhits_to_fit = (float) esl_opt_GetInteger(go, "--ltailn") * ((cfg->N * cfg->L) / 1000000.);
      tailp = nhits_to_fit / (float) h->n;
      if(tailp > 1.) ESL_FAIL(eslERANGE, errbuf, "--ltailn <n>=%d cannot be used, there's only %.3f hits per Mb in the histogram! Lower <n> or use --tailp.", esl_opt_GetInteger(go, "--ltailn"), (h->n / ((float) cfg->N * ((float) cfg->L) / 1000000.)));
    }
    else { /* glocal mode */
      nhits_to_fit = (float) esl_opt_GetInteger(go, "--gtailn") * ((cfg->N * cfg->L) / 1000000.);
      tailp = nhits_to_fit / (float) h->n;
      if(tailp > 1.) ESL_FAIL(eslERANGE, errbuf, "--gtailn <n>=%d cannot be used, there's only %.3f hits per Mb in the histogram! Lower <n> or use --tailp.", esl_opt_GetInteger(go, "--gtailn"), (h->n / ((float) cfg->N * ((float) cfg->L) / 1000000.)));
    }
  }

  esl_histogram_GetTailByMass(h, tailp, &xv, &n, &z); /* fit to right 'tailfit' fraction, 0.01 by default */
  if(n <= 1) ESL_FAIL(eslERANGE, errbuf, "fit_histogram(), too few points in right tailfit: %f fraction of histogram. Increase <x> with -L <x>.", tailp);
  esl_exp_FitComplete(xv, n, &(params[0]), &(params[1]));
  esl_histogram_SetExpectedTail(h, params[0], tailp, &esl_exp_generic_cdf, &params);

  /* printf("# Exponential fit to %.7f%% tail: lambda = %f\n", tailp*100.0, params[1]); */
  mu = params[0];
  lambda = params[1];
  if(isnan(lambda)) ESL_FAIL(eslERANGE, errbuf, "fit_histogram(), exp tail fit lambda is NaN, too few hits in histogram. Increase <x> with -L <x>.");
  if(isinf(lambda)) ESL_FAIL(eslERANGE, errbuf, "fit_histogram(), exp tail fit lambda is inf, too few hits in histogram. Increase <x> with -L <x>.");
  nrandhits = h->n; /* total number of hits in the histogram */

  /* print to output files if nec */
  if(cfg->exphfp != NULL) esl_histogram_Plot(cfg->exphfp, h);
  if(cfg->expqfp != NULL) esl_histogram_PlotQQ(cfg->expqfp, h, &esl_exp_generic_invcdf, params);
  if (cfg->expsfp != NULL) {
    esl_histogram_PlotSurvival(cfg->expsfp, h);
    esl_exp_Plot(cfg->expsfp, (params[0] - log(1./tailp) / params[1]), 0.693147, esl_exp_surv, h->xmin - 5., h->xmax + 5., 0.1); /* extrapolate mu */
  }
  esl_histogram_Destroy(h);

  *ret_mu     = mu;
  *ret_lambda = lambda;
  *ret_nrandhits = nrandhits;
  *ret_tailp = tailp;
  return eslOK;
}

/* Function: get_random_dsq()
 * Date:     EPN, Tue Sep 11 08:31:47 2007
 * 
 * Purpose:  Generate a random digitized seq and return it.
 *           Two possible modes:
 *           1. if(cfg->gc_freq == NULL && dnull != NULL) 
 *              use dnull disto (a double version of cm->null) to generate
 *           2. if(cfg->gc_freq != NULL && dnull == NULL) 
 *              use choose a GC frequency from cfg->gc_freq
 *              and generate with that
 *
 * Returns:  eslOK on success, ret_dsq filled with newly alloc'ed ESL_DSQ *,
 *           some other status code on failure.
 */
int
get_random_dsq(const struct cfg_s *cfg, char *errbuf, CM_t *cm, int L, ESL_RANDOMNESS *r, ESL_DSQ **ret_dsq)
{
  int      status;
  double   gc_comp;
  double  *distro = NULL;
  int      do_free_distro = FALSE;
  ESL_DSQ *dsq = NULL;

  /* contract check, make sure we're in a valid mode */
  if(cfg->gc_freq == NULL && cfg->dnull == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "get_random_dsq(), cfg->gc_freq == NULL and cfg->dnull == NULL");
  if(cfg->gc_freq != NULL && cfg->dnull != NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "get_random_dsq(), cfg->gc_freq != NULL and cfg->dnull != NULL");

  /* generate sequence */
  if      (cfg->gc_freq == NULL && cfg->dnull != NULL) distro = cfg->dnull;
  else if (cfg->gc_freq != NULL && cfg->dnull == NULL) {
    assert(cm->abc->K == 4);
    ESL_ALLOC(distro, sizeof(double) * cm->abc->K);
    do_free_distro = TRUE;
    gc_comp = 0.01 * esl_rnd_DChoose(r, cfg->gc_freq, GC_SEGMENTS);
    distro[1] = distro[2] = 0.5 * gc_comp;
    distro[0] = distro[3] = 0.5 * (1. - gc_comp);
  }
  /* generate sequence */
  ESL_ALLOC(dsq, sizeof(ESL_DSQ) * (L+2));
  if ((status = esl_rsq_xIID(r, distro, cm->abc->K, L, dsq) != eslOK)) return status;

  if (do_free_distro) free(distro);
  *ret_dsq = dsq;
  return eslOK;

 ERROR:
  return status;
}

/* Function: get_cmcalibrate_comlog_info
 * Date:     EPN, Mon Dec 31 14:59:52 2007
 *
 * Purpose:  Create the cmcalibrate command info and creation date info 
 *           to eventually be set in the CM's ComLog_t data structure.
 *
 * Returns:  eslOK on success, eslEINCOMPAT on contract violation.
 */
static int
get_cmcalibrate_comlog_info(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  int status;
  int i;
  uint32_t seed;
  int  temp;
  int  seedlen;
  char *seedstr;
  time_t date = time(NULL);

  if(cfg->ccom  != NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "get_cmcalibrate_comlog_info(), cfg->ccom  is non-NULL.");
  if(cfg->cdate != NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "get_cmcalibrate_comlog_info(), cfg->cdate is non-NULL.");
  
  /* Set the cmcalibrate command info, the cfg->ccom string */
  for (i = 0; i < go->optind; i++) { /* copy all command line options, but not the command line args yet, we may need to append '-s ' before the args */
    esl_strcat(&(cfg->ccom),  -1, go->argv[i], -1);
    esl_strcat(&(cfg->ccom),  -1, " ", 1);
  }
  /* if -s NOT enabled, we need to append the seed info also */
  seed = esl_randomness_GetSeed(cfg->r);
  if(esl_opt_IsUsed(go, "--seed")) { /* -s was used on command line, we'll do a sanity check */
    if(seed != esl_opt_GetInteger(go, "--seed")) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "get_cmcalibrate_comlog_info(), cfg->r's seed is %" PRIu32 ", but -s was enabled with argument: %" PRIu32 "!, this shouldn't happen.", seed, esl_opt_GetInteger(go, "--seed"));
  } else {
    temp = seed; 
    seedlen = 1; 
    while(temp > 0) { temp/=10; seedlen++; } /* determine length of stringized version of seed */
    seedlen += 4; /* strlen(' -s ') */
    ESL_ALLOC(seedstr, sizeof(char) * (seedlen+1));
    sprintf(seedstr, " -s %" PRIu32 " ", seed);
    esl_strcat((&cfg->ccom), -1, seedstr, seedlen);
    free(seedstr);
  }

  for (i = go->optind; i < go->argc; i++) { /* copy command line args yet */
    esl_strcat(&(cfg->ccom), -1, go->argv[i], -1);
    if(i < (go->argc-1)) esl_strcat(&(cfg->ccom), -1, " ", 1);
  }
  
  /* Set the cmcalibrate call date, the cfg->cdate string */
  if((status = esl_strdup(ctime(&date), -1, &(cfg->cdate))) != eslOK) goto ERROR;
  esl_strchop(cfg->cdate, -1); /* doesn't return anything but eslOK */

  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "get_cmcalibrate_comlog_info() error status: %d, probably out of memory.", status);
  return status; 
}


/* Function: update_comlog
 * Date:     EPN, Mon Dec 31 15:14:26 2007
 *
 * Purpose:  Update the CM's comlog info to reflect the current
 *           cmcalibrate call.
 *
 * Returns:  eslOK on success, eslEINCOMPAT on contract violation.
 */
static int
update_comlog(const ESL_GETOPTS *go, char *errbuf, char *ccom, char *cdate, CM_t *cm)
{
  int status;
  if(ccom  == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "update_comlog(), ccom  is non-NULL.");
  if(cdate == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "update_comlog(), cdate is non-NULL.");

  /* free all cmcalibrate comlog info, we're about to overwrite any information that any previous cmcalibrate
   * call could have written to the cm file.
   */
  if(cm->comlog->ccom  != NULL)  { free(cm->comlog->ccom);  cm->comlog->ccom = NULL;  }
  if(cm->comlog->cdate != NULL)  { free(cm->comlog->cdate); cm->comlog->cdate = NULL; }
  
  if((status = esl_strdup(ccom, -1, &(cm->comlog->ccom)))  != eslOK) goto ERROR; 
  if((status = esl_strdup(cdate,-1, &(cm->comlog->cdate))) != eslOK) goto ERROR; 
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "update_comlog() error status: %d, probably out of memory.", status);
  return status; 
}


/* Function: set_dnull
 * Date:     EPN, Thu Jan 24 09:48:54 2008
 *
 * Purpose:  Allocate, fill and return dnull, a double version of cm->null used
 *           for generating random seqs.
 *
 * Returns:  eslOK on success
 */
static int
set_dnull(struct cfg_s *cfg, CM_t *cm, char *errbuf)
{
  int status;
  int i;

  if(cfg->dnull != NULL) { free(cfg->dnull); }
  ESL_ALLOC(cfg->dnull, sizeof(double) * cm->abc->K);
  for(i = 0; i < cm->abc->K; i++) cfg->dnull[i] = (double) cm->null[i];
  esl_vec_DNorm(cfg->dnull, cm->abc->K);    

  return eslOK;

 ERROR:
  ESL_FAIL(eslEINCOMPAT, errbuf, "set_dnull(), memory allocation error.");
}


/* Function: print_run_info
 * Date:     EPN, Sun Mar  2 16:57:25 2008
 *
 * Purpose:  Print information on this run of cmcalibrate.
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
  if((status = GetDate    (errbuf, &date))        != eslOK) return status;

  fprintf(stdout, "%-10s %s\n",  "# command:", command);
  fprintf(stdout, "%-10s %s\n",  "# date:",    date);
  fprintf(stdout, "%-10s %" PRIu32 "\n", "# seed:", esl_randomness_GetSeed(cfg->r));
  if(cfg->nproc > 1) fprintf(stdout, "# %-8s %d\n", "nproc:", cfg->nproc);

  free(command);
  free(date);
  return eslOK;
}

/* Function: get_command
 * Date:     EPN, Fri Jan 25 13:56:10 2008
 *
 * Purpose:  Return the command used to call cmscore
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

/* Function: print_post_calibration_info
 * Date:     EPN, Wed Mar  5 05:25:02 2008
 * Purpose:  Print info about calibration for a CM we just calibrated including
 *           timings. 
 */
int print_post_calibration_info(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, FILE *fp, CM_t *cm, double *exp_psecA, double *exp_asecA)
{
  char  time_buf[128];	      /* for printing run time */
  int   exp_mode;             /* counter over exp tail modes */
  double total_psec = 0.;     /* predicted number of seconds for cm, all stages */
  double total_asec = 0.;     /* actual number of seconds for cm, all stages */
  float      L_Mb;            /* total seq length we'll calibrate exp tails on in Mb */
  ExpInfo_t *exp;             /* pointer to current exp tail info, for convenience */

  if(exp_psecA == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "print_post_calibration_info, exp_psecA is NULL");
  if(exp_asecA == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "print_post_calibration_info, exp_asecA is NULL");

  fprintf(fp, "#\n");
  fprintf(fp, "# Post-calibration info for CM %d: %s\n", cfg->ncm, cm->name);
  fprintf(fp, "#\n");
  fprintf(fp, "# Exponential tail fitting:\n");
  fprintf(fp, "#\n");
  fprintf(fp, "# %-3s  %3s  %3s %7s %6s %6s %7s %21s\n",       "",    "",    "",  "",   "",       "",        "",            "    running time     ");
  fprintf(fp, "# %-3s  %3s  %3s %7s %6s %6s %7s %21s\n",       "",    "",    "",  "",   "",       "",        "",            "---------------------");
  fprintf(fp, "# %-3s  %3s  %3s %7s %6s %6s %7s %10s %10s\n","mod", "cfg", "alg", "L (Mb)",  "mu",     "lambda",   "nhits", "predicted",     "actual");
  fprintf(fp, "# %3s  %3s  %3s %7s %6s %6s %7s %10s %10s\n", "---", "---", "---", "-------", "------", "------", "-------", "----------", "----------");

  for(exp_mode = 0; exp_mode < EXP_NMODES; exp_mode++) {
    L_Mb = ((float) cfg->N * (float) cfg->L) / 1000000.;
    total_psec += exp_psecA[exp_mode];
    total_asec += exp_asecA[exp_mode];
    FormatTimeString(time_buf, exp_psecA[exp_mode], FALSE);
    exp = cfg->expAA[cfg->ncm-1][exp_mode];
    fprintf(fp, "  %-12s %7.2f %6.2f %6.3f %7d %10s", DescribeExpMode(exp_mode), L_Mb, exp->mu_orig, exp->lambda, exp->nrandhits, time_buf);
    FormatTimeString(time_buf, exp_asecA[exp_mode], FALSE);
    fprintf(fp, " %10s\n", time_buf);
  }

  fprintf(fp, "#\n");
  FormatTimeString(time_buf, total_psec, FALSE);
  fprintf(fp, "# total predicted time for CM: %s\n", time_buf);
  FormatTimeString(time_buf, total_asec, FALSE);
  fprintf(fp, "# total actual    time for CM: %s\n", time_buf);
  
  fflush(fp);
  return eslOK;
}
  
/* Function: estimate_time_for_exp_round
 * Date:     EPN, Wed Mar  5 05:46:45 2008
 * Purpose:  Estimate search time for round of exp tail fitting
 *           of exp mode <exp_mode>. This is done by actually 
 *           searching a sequence with the appropriate algorithm. 
 *           The length of the sequence to search is set such 
 *           that it should take about <targ_sec> seconds.
 */
int estimate_time_for_exp_round(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, CM_t *cm, int exp_mode, double *ret_sec_per_res)
{
  int    status;
  int    L;                /* length of sequence we'll generate and search to get time estimate */
  double psec_per_Mc;      /* rough prediction at seconds per Mc based on empirical run times I've witnessed, 
                            * doesn't need to be very accurate as we just use it to set length of seq to search to get real prediction */
  float  Mc;               /* millions of DP calculations we're going to do */
  float  Mc_per_res_corr;  /* millions of dp calcs per residue, if searching with CM, corrects         for first W residues requiring less dp calcs */
  float  Mc_per_res;       /* millions of dp calcs per residue, if searching with CM, does not correct for first W residues requiring less dp calcs */
  int    orig_search_opts; /* cm->search_opts when function was entered */
  float  sec_per_res;      /* seconds per residue */
  float  targ_sec = 0.1;   /* target number of seconds our timing expt will take */
  int    Lmin = 100;       /* minimum number of residues to search to get timing */
  int    Lmax = 10000;     /* maximum number of residues to search to get timing */
  int    use_qdb;          /* TRUE if we're using QDB, FALSE if not */
  ESL_DSQ *dsq;            /* the random seq we'll create and search to get predicted time */
  ESL_STOPWATCH *w  = esl_stopwatch_Create(); /* for timings */

  if(w == NULL) ESL_FAIL(status, errbuf, "estimate_time_for_exp_round(): memory error, stopwatch not created.\n");

  orig_search_opts = cm->search_opts; /* we'll restore this at end of func, just to be safe */
  if(ExpModeIsInside(exp_mode)) cm->search_opts |= CM_SEARCH_INSIDE;
  else                          cm->search_opts &= ~CM_SEARCH_INSIDE;

  if(cm->smx == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "estimate_time_for_exp_round(), cm->smx is NULL");
  use_qdb  = (cm->search_opts & CM_SEARCH_QDB) ? TRUE : FALSE;

  if(use_qdb) { if((status = cm_CountSearchDPCalcs(cm, errbuf, 10*cm->W, NULL,               NULL,               cm->W, TRUE,  NULL, &Mc_per_res_corr)) != eslOK) return status; }
  else        { if((status = cm_CountSearchDPCalcs(cm, errbuf, 10*cm->W, cm->qdbinfo->dmin2, cm->qdbinfo->dmax2, cm->W, TRUE,  NULL, &Mc_per_res_corr)) != eslOK) return status; }

  psec_per_Mc = (cm->search_opts & CM_SEARCH_INSIDE) ? (1. /  75.) : (1. / 275.);  /*  75 Mc/S inside;  275 Mc/S CYK */
  /* determine L that will take about <targ_sec> seconds */
  L = targ_sec / (psec_per_Mc * Mc_per_res_corr);
  L = ESL_MAX(L, Lmin); /* we have to search at least <Lmin> residues */
  L = ESL_MIN(L, Lmax); /* we want to search at most  <Lmax> residues */
  /* now determine exactly how many dp calculations we'd do if we search L residues, 
   * this won't be the same as Mc_per_res * L b/c Mc_per_res from cm_GetNCalcsPerResidueGivenBeta
   * b/c that was calculated after correcting for the fact that the first W residues have fewer
   * DP calcs than all other residues, b/c d < W for first W residues.
   */  
  if(use_qdb) { if((status = cm_CountSearchDPCalcs(cm, errbuf, L, NULL,               NULL,               cm->W, FALSE,  NULL, &Mc_per_res)) != eslOK) return status; }
  else        { if((status = cm_CountSearchDPCalcs(cm, errbuf, L, cm->qdbinfo->dmin2, cm->qdbinfo->dmax2, cm->W, FALSE,  NULL, &Mc_per_res)) != eslOK) return status; }

  /* FALSE says don't correct for fewer dp calcs for first W residues, we want to know how many total DP calcs
   * there will be in L residues */
  Mc = Mc_per_res * L;

  esl_stopwatch_Start(w);
  /* simulate a workunit, generate a sequence, search it, and remove overlaps */
  /*printf("exptL: %d\n", L);*/
  if(esl_opt_GetBoolean(go, "--random")) { 
    if((status = get_random_dsq(cfg, errbuf, cm, L, cfg->r_est, &dsq)) != eslOK) return status;
  }
  else { 
    if((status = SampleGenomicSequenceFromHMM(cfg->r_est, cm->abc, errbuf, cfg->ghmm_sA, cfg->ghmm_tAA, cfg->ghmm_eAA, cfg->ghmm_nstates, L, &dsq)) != eslOK) cm_Fail(errbuf);
  }

  if((status = process_search_workunit(cm, errbuf, dsq, L, cfg->sc_cutoff, NULL)) != eslOK) cm_Fail(errbuf);
  esl_stopwatch_Stop(w);

  cm->search_opts = orig_search_opts;
  sec_per_res = w->user * (Mc_per_res_corr / Mc);
  ESL_DPRINTF1(("L: %d\n", L));
  ESL_DPRINTF1(("w->user: %f\n", w->user));
  ESL_DPRINTF1(("sec_per_res: %f\n", sec_per_res));
  ESL_DPRINTF1(("Mc_per_res: %f\n", Mc_per_res));
  ESL_DPRINTF1(("Mc: %f\n", Mc));

  /*printf("L: %d\n", L);
    printf("w->user: %f\n", w->user);
    printf("sec_per_res: %f\n", sec_per_res);
    printf("Mc_per_res: %f\n", Mc_per_res);
    printf("Mc: %f\n", Mc);*/
  
  if(w != NULL) esl_stopwatch_Destroy(w);

  free(dsq);
  
  if(ret_sec_per_res != NULL) *ret_sec_per_res = sec_per_res;
  return eslOK;
}


/* Function: print_per_cm_column_headings
 * Date:     EPN, Tue Jan  8 05:51:47 2008
 *
 * Purpose:  Print per-CM info to stdout. 
 *
 * Returns:  eslOK on success
 */
static int
print_per_cm_column_headings(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
  printf("#\n");
  if(! esl_opt_IsOn(go, "--forecast")) { 
    printf("# Calibrating CM %d: %s\n", cfg->ncm, cm->name);
    printf("#\n");
    printf("# %8s  %3s  %3s  %3s  %9s %22s\n",           "",      "",    "",    "",    "", "     running time    ");
    printf("# %8s  %3s  %3s  %3s  %9s %22s\n",           "",      "",    "",    "",    "", "----------------------");
    printf("# %-8s  %3s  %3s  %3s  %9s %10s  %10s\n", "stage",    "mod", "cfg", "alg", "L (Mb)", "predicted", "actual");
    printf("# %8s  %3s  %3s  %3s  %9s %10s  %10s\n",  "--------", "---", "---", "---", "---------", "----------", "----------");
  }
  else { 
    printf("# Forecasting time for %d processor(s) to calibrate CM %d: %s\n", esl_opt_GetInteger(go, "--forecast"), cfg->ncm, cm->name);
    printf("#\n");
    printf("# %-8s  %3s  %3s  %3s  %9s %14s\n", "stage",    "mod", "cfg", "alg", "L (Mb)", "predicted time");
    printf("# %8s  %3s  %3s  %3s  %9s %14s\n",  "--------", "---", "---", "---", "---------", "--------------");
  }
  return eslOK;
}


/* Function: print_per_cm_summary
 * Date:     EPN, Thu Mar  6 10:35:49 2008
 *
 * Purpose:  Print per-CM summary to stdout. 
 *
 * Returns:  eslOK on success
 */
static int
print_per_cm_summary(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, double psec, double asec)
{
  char  time_buf[128];	      /* for printing run time */
  if(! esl_opt_IsOn(go, "--forecast")) { 
    FormatTimeString(time_buf, psec, FALSE);
    printf("# %8s  %3s  %3s  %3s  %9s %10s  %10s\n", "--------", "---", "---", "---", "---------", "----------", "----------");
    printf("# %-8s  %3s  %3s  %3s  %9s %10s",       "all",        "-",   "-",   "-",         "-", time_buf);
    FormatTimeString(time_buf, asec, FALSE);
    printf("  %10s  (hr:min:sec)\n", time_buf);
  }
  else { 
    FormatTimeString(time_buf, psec, FALSE);
    printf("# %8s  %3s  %3s  %3s  %9s %14s\n", "--------", "---", "---", "---",  "---------", "--------------");
    printf("# %-8s  %3s  %3s  %3s  %9s %14s  (hr:min:sec)\n", "all",        "-",   "-",   "-",         "-",     time_buf);
  }
  return eslOK;
}

/* Function: print_exp_line
 * Date:     EPN, Thu Mar  6 13:38:01 2008
 *
 * Purpose:  Print a line describing exp tail fitting for a given mode.
 */
static int
print_exp_line(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, int exp_mode, int N, int L, double psec)
{
  char  time_buf[128];	      /* for printing run time */
  float L_Mb;              /* exp seq length in Mb */

  FormatTimeString(time_buf, psec, FALSE);

  L_Mb =  (float) N * (float) L; 
  L_Mb /= 1000000.;

  if(esl_opt_IsOn(go, "--forecast")) { 
    printf("  %-8s  %-12s  %9.2f %14s\n",              "exp tail", DescribeExpMode(exp_mode), L_Mb, time_buf);
  }
  else { 
    printf("  %-8s  %-12s  %9.2f %10s",              "exp tail", DescribeExpMode(exp_mode), L_Mb, time_buf);
  }
  fflush(stdout);
  return eslOK;
}

/* Function: generate_sequences()
 * Date:     EPN, Fri Dec 16 09:41:32 2011
 *
 * Purpose:  Generate all sequences to be used for fitting a single 
 *           model in a single mode, and return them as a single
 *           ESL_SQ_BLOCK in <*ret_sq_block>.
 *
 * Returns:  eslOK on success, filled block is in *ret_sq_block.
 *           eslEMEM if out of memory.
 */
int
generate_sequences(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, ESL_SQ_BLOCK **ret_sq_block)
{
  int           status;
  ESL_SQ_BLOCK *sq_block = NULL;
  ESL_SQ       *sq       = NULL;
  ESL_DSQ      *dsq      = NULL;
  int           i;
  int           namelen = strlen("irrelevant");
  sq_block = esl_sq_CreateDigitalBlock(cfg->N, cfg->abc);

  for(i = 0; i < cfg->N; i++) { 
    /* generate the dsq */
    if(esl_opt_GetBoolean(go, "--random")) {
      if((status = get_random_dsq(cfg, errbuf, cm, cfg->L, cfg->r, &dsq)) != eslOK) goto ERROR;
    }
    else { 
      if((status = SampleGenomicSequenceFromHMM(cfg->r, cfg->abc, errbuf, cfg->ghmm_sA, cfg->ghmm_tAA, cfg->ghmm_eAA, cfg->ghmm_nstates, cfg->L, &dsq)) != eslOK) goto ERROR;
    }

    /* TEMP */
    ///int x; for(x = 1; x <= cfg->L; x++) dsq[x] = x % 4;

    /* Copy dsq we just created into sq_block->list+i 
     * We can't use esl_sq_CreateDigitalFrom() bc sq_block->list already contains a contiguous set of ESL_SQ objects 
     */
    sq = sq_block->list + i;
    if ((status = esl_abc_dsqdup(dsq, cfg->L, &(sq->dsq))) != eslOK) goto ERROR;
    sq->L = cfg->L; 
    /* set the name (we don't actually use it, but all valid ESL_SQ objects are supposed to have it) */
    ESL_ALLOC(sq->name, sizeof(char) * (namelen+1));
    strcpy(sq->name, "irrelevant");
    sq->nalloc = namelen;
    /*printf("sq->L: %" PRId64 " sq_block->list[i].L: %" PRId64 " first 10 letters: %d%d%d%d%d%d%d%d%d%d\n", sq->L, sq_block->list[i].L, 
      sq->dsq[1], sq->dsq[2], sq->dsq[3], sq->dsq[4], sq->dsq[5], sq->dsq[6], sq->dsq[7], sq->dsq[8], sq->dsq[9], sq->dsq[10]);*/

    sq_block->count++;
  }
  sq_block->first_seqidx = 0;
  sq_block->complete     = TRUE;

  *ret_sq_block = sq_block;
  return eslOK;

 ERROR: 
  if(sq_block != NULL) {
    esl_sq_DestroyBlock(sq_block); 
  }
  *ret_sq_block = NULL;
  return status; 
}  

/* Function: process_search_workunit()
 * Date:     EPN, Thu Dec  8 13:48:02 2011
 *
 * Purpose:  Perform search workunit, which consists of a CM, digitized sequence
 *           and indices i and j. The job is to search dsq from i..j and return 
 *           search results in <*ret_results>. Called by cmsearch and cmcalibrate,
 *           which is why it's here and not local in cmsearch.c.
 *
 * Args:     cm              - the covariance model, must have valid searchinfo (si).
 *           errbuf          - char buffer for reporting errors
 *           dsq             - the digitized sequence
 *           L               - length of target sequence 
 *           cutoff          - minimum bit score cutoff to report
 *           ret_th          - search_results_t to create and fill
 *
 * Returns:  eslOK on succes;
 *           <ret_th> is filled with a newly alloc'ed and filled CM_TOPHITS structure, must be freed by caller
 */
int
process_search_workunit(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float cutoff, CM_TOPHITS **ret_th)
{
  int status;
  CM_TOPHITS *th = NULL;
  int use_qdbs = (cm->search_opts & CM_SEARCH_QDB) ? TRUE : FALSE;

  th = cm_tophits_Create();
  if(th == NULL) ESL_FAIL(eslEMEM, errbuf, "out of memory");

  if(cm->search_opts & CM_SEARCH_INSIDE) { 
    if((status = FastIInsideScan(cm, errbuf, cm->smx,                   
				 use_qdbs ? SMX_QDB2_LOOSE : SMX_NOQDB, /* qdbidx, indicates which QDBs to use */
				 dsq, 1, L,                             /* sequence, bounds */
				 cutoff,                                /* minimum score to report */
				 th,                                    /* hitlist to add to */
				 cm->search_opts & CM_SEARCH_NULL3,     /* do the NULL3 correction? */
				 0., NULL, NULL,                        /* vars for redefining envelopes, which we won't do */
				 NULL, NULL))                           /* ret_vsc, ret_sc, irrelevant here */
       != eslOK) return status;
  }
  else { 
   if((status = FastCYKScan(cm, errbuf, cm->smx, 
			    use_qdbs ? SMX_QDB2_LOOSE : SMX_NOQDB, /* qdbidx, indicates which QDBs to use */
			    dsq, 1, L,                             /* sequence, bounds */
			    cutoff,                                /* minimum score to report */
			    th,                                    /* hitlist to add to */
			    cm->search_opts & CM_SEARCH_NULL3,     /* do the NULL3 correction? */
			    0., NULL, NULL,                        /* vars for redefining envelopes, which we won't do */
			    NULL, NULL))                           /* ret_vsc, ret_sc, irrelevant here */
       != eslOK) return status;
  }
  /* we don't have to remove overlaps, that's already been done in
   * FastCYKScan() or FastIInsideScan() 
   */
  
  if(ret_th != NULL) *ret_th = th;
  else                cm_tophits_Destroy(th);
  return eslOK;
}


/*****************************************************************
 * @LICENSE@
 *****************************************************************/

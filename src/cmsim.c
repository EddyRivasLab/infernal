/* cmsim: scoring CMs against simulated sequences.
 * 
 * Main testbed for exploring the statistical behavior of Infernal
 * scores on random sequences, and importance sampling.
 * 
 * EPN, Fri Apr 29 14:00:58 2011
 */
#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "easel.h"
#include <esl_getopts.h>
#include <esl_histogram.h>
#include <esl_random.h>
#include <esl_stats.h>
#include <esl_stopwatch.h>
#include <esl_vectorops.h>
#include <esl_wuss.h>

#include "hmmer.h"

#include "funcs.h"		/* function declarations                */
#include "structs.h"		/* data structures, macros, #define's   */

#define ALPHOPTS "--rna,--dna"                         /* Exclusive options for alphabet choice */
#define OUTOPTS  "-u,-c,-a,--ahmm,--shmm"              /* Exclusive options for output */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles  reqs  incomp  help  docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL, NULL, "show brief help on version and usage",   1 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL, NULL, "verbose: print scores",                             1 },
  { "-i",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL, NULL, "do importance sampling",                            1 },
  { "-g",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL, NULL, "put CM in glocal mode",                             1 },
  { "-L",        eslARG_INT,  "10000", NULL, "n>0",     NULL,  NULL, NULL, "length of random target seqs",                      1 },
  { "-N",        eslARG_INT,     "10", NULL, "n>0",     NULL,  NULL, NULL, "number of random target seqs",                      1 },
  { "--beta",    eslARG_REAL,  "1e-7", NULL, "0<x<1",   NULL,  NULL, NULL,     "set tail loss prob for QDB calculation to <x>", 1 },
  { "--noqdb",   eslARG_NONE,   FALSE, NULL, NULL,      NULL,  NULL, "--beta", "do not use QDBs", 1 },
  { "--exp",     eslARG_REAL,   NULL,  NULL, "x>0",     NULL,  NULL, NULL, "exponentiate CM probabilities by <x> before sampling",  1 },

  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

/* struct cfg_s : "Global" application configuration shared by all threads/processes
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 * This strategy is used despite the fact that a MPI version of cmemit does not
 * yet exist! 
 */
struct cfg_s {
  char         *cmfile;	        /* name of input CM file  */ 
  CMFILE       *cmfp;		/* open input CM file stream       */
  ESL_ALPHABET *abc;		/* digital alphabet for CM */
  ESL_RANDOMNESS *r;            /* source of randomness */
  int           ncm;            /* number CM we're at in file */
};

static char usage[]  = "[-options] <cmfile>";
static char banner[] = "score random sequences with a covariance model";

static int  init_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf);

static void master(const ESL_GETOPTS *go, struct cfg_s *cfg);

static int initialize_cm(const ESL_GETOPTS *go, const struct cfg_s *cfg, CM_t *cm, char *errbuf);
static int print_run_info(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf);
static int get_command(const ESL_GETOPTS *go, char *errbuf, char **ret_command);

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go = NULL;   /* command line processing                     */
  ESL_STOPWATCH   *w  = esl_stopwatch_Create();
  if(w == NULL) cm_Fail("Memory allocation error, stopwatch could not be created.");
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
      puts("\nmiscellaneous output options are:");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\nexpert options:");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); 
      exit(0);
    }
  if (esl_opt_ArgNumber(go) != 1) 
    {
      puts("Incorrect number of command line arguments.");
      esl_usage(stdout, argv[0], usage);
      puts("\n  where basic options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      printf("\nTo see more help on other available options, do %s -h\n\n", argv[0]);
      exit(1);
    }
  /* Initialize what we can in the config structure (without knowing the alphabet yet).
   * We could assume RNA, but this HMMER3 based approach is more general.
   */
  cfg.cmfile     = esl_opt_GetArg(go, 1); 
  cfg.cmfp       = NULL;	           /* opened in init_cfg() */
  cfg.abc        = NULL;	           /* created in init_cfg() */
  cfg.r          = NULL;	           /* created in init_cfg() */

  cm_banner(stdout, argv[0], banner);

  /* do work */
  master(go, &cfg);

  /* Clean up the cfg. 
   */
  if (cfg.abc   != NULL) { esl_alphabet_Destroy(cfg.abc); cfg.abc = NULL; }
  if (cfg.cmfp  != NULL) CMFileClose(cfg.cmfp);
  if (cfg.r     != NULL) esl_randomness_Destroy(cfg.r);

  esl_getopts_Destroy(go);
  esl_stopwatch_Stop(w);
  printf("#\n");
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  esl_stopwatch_Destroy(w);
  return 0;
}

/* init_cfg()
 * Already set:
 *    cfg->cmfile  - command line arg 1
 * Sets: 
 *    cfg->cmfp    - open CM file
 *    cfg->r       - source of randomness
 */
static int
init_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
  /* open CM file for reading */
  if ((cfg->cmfp = CMFileOpen(cfg->cmfile, NULL)) == NULL)
    ESL_FAIL(eslFAIL, errbuf, "Failed to open covariance model save file %s\n", cfg->cmfile);

  /* create RNG */
  cfg->r = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));

  if (cfg->r       == NULL) ESL_FAIL(eslEINVAL, errbuf, "Failed to create random number generator: probably out of memory");

  return eslOK;
}

/* master()
 * The serial version of cmsim. (There is no parallel version yet).
 * For each CM, generate random sequences and search them.
 * 
 * We only return if successful. All errors are handled immediately and fatally with cm_Fail().
 */
static void
master(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int      status;
  char     errbuf[cmERRBUFSIZE];
  CM_t    *cm = NULL;

  if ((status = init_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);
  if ((status = print_run_info (go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);

  cfg->ncm = 0;

  while ((status = CMFileRead(cfg->cmfp, errbuf, &(cfg->abc), &cm)) == eslOK)
  {
    if (cm == NULL) cm_Fail("Failed to read CM from %s -- file corrupt?\n", cfg->cmfile);
    cfg->ncm++;
    if((status = initialize_cm(go, cfg, cm, errbuf)) != eslOK) cm_Fail(errbuf);

    printf("CM 1: %s\n", cm->name);
    FreeCM(cm);
  }
  if(status != eslEOF) cm_Fail(errbuf);
  return;
}

/* initialize_cm()
 * Setup the CM based on the command-line options/defaults;
 * only set flags and a few parameters. ConfigCM() configures
 * the CM.
 */
static int
initialize_cm(const ESL_GETOPTS *go, const struct cfg_s *cfg, CM_t *cm, char *errbuf)
{
  int status;

  /* Update cfg->cm->config_opts and cfg->cm->align_opts based on command line options */
  cm->beta_qdb = esl_opt_GetReal(go, "--beta");
  if(! esl_opt_GetBoolean(go, "-g")) {
    cm->config_opts |= CM_CONFIG_LOCAL;
  }
  if(! esl_opt_GetBoolean(go, "--noqdb")) { 
    cm->config_opts |= CM_CONFIG_QDB;
  }

  if((status = ConfigCM(cm, errbuf, FALSE, NULL, NULL)) != eslOK) return status; /* FALSE says do not calculate W unless nec b/c we're using QDBs */

  if(esl_opt_IsOn(go, "--exp"))        ExponentiateCM(cm, esl_opt_GetReal(go, "--exp"));
  
  return eslOK;
}

/* Function: print_run_info
 * Date:     EPN, Mon Mar  3 06:01:13 2008
 *
 * Purpose:  Print information on this run of cmsim.
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

  fprintf(stdout, "#\n");
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





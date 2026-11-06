/* cmconvert: converting covariance model files to Infernal-1.1 CM format.
 *
 * EPN, Fri Jul  1 05:11:15 2011
 * SRE, Thu Oct 16 08:57:43 2008 [janelia] (hmmconvert.c)
 */
#include <esl_config.h>
#include <p7_config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#ifdef HMMER_THREADS
#include "esl_threads.h"
#endif

#include "hmmer.h"

#include "infernal.h"

#define OUTOPTS "-a,-b,-1,--mlhmm,--fhmm,--cp9ascm"

#if defined (HMMER_THREADS)
#define CPUOPTS     NULL
#else
#define CPUOPTS     "--cpu"
#endif

static ESL_OPTIONS options[] = {
  /* name               type  default   env  range   toggles        reqs      incomp  help                                                         docgroup */
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,       NULL,       NULL, "show brief help and exit",                                         1 },
  { "--version", eslARG_NONE,   FALSE, NULL, NULL,      NULL,       NULL,       NULL, "show version info and exit",                                       1 },
  { "--devhelp", eslARG_NONE,   FALSE, NULL, NULL,      NULL,       NULL,       NULL, "show list of otherwise hidden developer/expert options",           1 },
  { "-a",        eslARG_NONE,"default",NULL, NULL,   OUTOPTS,       NULL,       NULL, "ascii:  output models in Infernal ASCII format",                   1 },
  { "-b",        eslARG_NONE,   FALSE, NULL, NULL,   OUTOPTS,       NULL,       NULL, "binary: output models in Infernal binary format",                  1 },
  { "-1",        eslARG_NONE,   FALSE, NULL, NULL,   OUTOPTS,       NULL,       NULL, "output backward-compatible v0.7-v1.0.2 ASCII format",   1 },
  { "-o",        eslARG_OUTFILE,FALSE, NULL, NULL,      NULL,       NULL,       NULL, "save CM file to file <f>, not stdout",                             1 },
  { "--mlhmm",   eslARG_NONE,   FALSE, NULL, NULL,   OUTOPTS,       NULL,       NULL, "output maximum likelihood HMM for CM in HMMER3 format",            1 },
  { "--fhmm",    eslARG_NONE,   FALSE, NULL, NULL,   OUTOPTS,       NULL,       NULL, "output filter HMM for CM in HMMER3 format (1/e: the ML HMM)",     1 },
  { "--cp9ascm", eslARG_NONE,   FALSE, NULL, NULL,   OUTOPTS,       NULL,       NULL, "output CM's ML CP9 HMM as an unstructured (MATL-chain) CM", 1 },
  { "--outfmt",  eslARG_STRING, NULL,  NULL, NULL,      NULL,       NULL,"-1,--mlhmm,--fhmm,--cp9ascm", "output format 1/a or 1/e (default: 1/e, or 1/a for a 1/a CM)", 1 },
  { "--dropknots",eslARG_NONE,   FALSE, NULL, NULL,      NULL,       NULL,       NULL, "allow dropping pseudoknots the output format can't store", 1 },
  { "--upgrade", eslARG_NONE,   FALSE, NULL, NULL,      NULL,       NULL,"-1,--mlhmm,--fhmm,--cp9ascm,--outfmt", "convert older-format CMs to the current format", 1 },
  { "--prior",   eslARG_INFILE, NULL,  NULL, NULL,      NULL,"--upgrade", "--noprior", "with --upgrade: CMs were built with prior file <f>",      1 },
  { "--noprior", eslARG_NONE,   FALSE, NULL, NULL,      NULL,"--upgrade",  "--prior", "with --upgrade: convert --prior-built CMs without their prior", 1 },
  { "--cpu",     eslARG_INT,  CMNCPU,"INFERNAL_NCPU","n>=0",NULL,       NULL,  CPUOPTS, "CPUs for calibration (--upgrade, v1.0 input)",  1 },
  /* options for controlling p7 per-node band pad computation: shown only by --devhelp */
  { "--p7pad",      eslARG_NONE,    FALSE, NULL, NULL,    NULL,  NULL,    "--no-p7pad", "compute p7 per-node band pads (O(M^2) memory)",   2 },
  { "--no-p7pad",   eslARG_NONE,    FALSE, NULL, NULL,    NULL,  NULL,       "--p7pad", "don't compute p7 per-node band pads (default)",   2 },
  { "--p7pad-N",    eslARG_INT,    "1000", NULL, "n>0",   NULL,  "--p7pad",       NULL, "number of samples for p7 pad simulation",        2 },
  { "--p7pad-q",    eslARG_REAL,   "0.99", NULL, "0<x<=1",NULL,  "--p7pad",       NULL, "quantile for p7 pad deficit distribution",       2 },
  { "--p7pad-seed", eslARG_INT,      "42", NULL, "n>=0",  NULL,  "--p7pad",       NULL, "RNG seed for p7 pad simulation (0=arbitrary)",   2 },
  { "--p7pad-cpu",  eslARG_INT,    CMNCPU, NULL, "n>=0",  NULL,  "--p7pad",       NULL, "number of CPUs for p7 pad simulation (0=serial)", 2 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <cmfile>";
static char banner[] = "convert CM file to a different Infernal format";

static int  configure_model(CM_t *cm, char *errbuf, int set_filter);
static int  upgrade_model(const ESL_GETOPTS *go, CM_t *cm, int format, const Prior_t *user_pri, char *errbuf);
static int  comlog_cmbuild_option(const char *comlog, const char *opt, char **opt_arg);
static int  comlog_cmbuild_number(const char *comlog, const char *opt, double *ret_x);

int
main(int argc, char **argv)
{
  ESL_GETOPTS   *go      = NULL;
  ESL_ALPHABET  *abc     = NULL;
  char          *cmfile  = NULL;
  CM_FILE       *cmfp    = NULL;
  CM_t          *cm      = NULL;
  FILE          *ofp     = NULL;
  int            fmtcode = -1;	/* -1 = cm_file_DefaultFormat(): INFERNAL1/e, or INFERNAL1/a for a CM with a stored filter HMM; set by --outfmt */
  int            status;
  char           errbuf[eslERRBUFSIZE];
  Prior_t       *user_pri = NULL;   /* --prior <f>: the prior the CMs were built with, for --upgrade */

  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK ||
      esl_opt_VerifyConfig(go)               != eslOK)
    {
      printf("Failed to parse command line: %s\n", go->errbuf);
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }
  if (esl_opt_GetBoolean(go, "-h") || esl_opt_GetBoolean(go, "--devhelp"))
    {
      cm_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      puts("\nOptions:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      if (esl_opt_GetBoolean(go, "--devhelp")) {
        puts("\nOptions for p7 per-node band pad computation:");
        esl_opt_DisplayHelp(stdout, go, 2, 2, 80);
      }
      exit(0);
    }
  if (esl_opt_GetBoolean(go, "--version"))
    {
      if (argc != 2) esl_fatal("Incorrect usage: to get version info, use --version alone");
      esl_printf("%s %s\n", "cmconvert", INFERNAL_VERSION);
      exit(0);
    }
  if (esl_opt_ArgNumber(go) != 1)
    {
      puts("Incorrect number of command line arguments.");
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }
  cmfile = esl_opt_GetArg(go, 1);

  /* --outfmt <s>: choose the 1.x output format for -a/-b. Default (NULL) writes
   * the default format (fmtcode -1, cm_file_DefaultFormat(): INFERNAL1/e, or for
   * a CM read with a stored filter HMM, INFERNAL1/a, unless --upgrade). 1/a
   * predates the consensus pseudoknot annotation, so writing to it drops
   * pknots (handled below). The development formats 1/b, 1/c and 1/d were
   * never released, and can't be selected (brief 26_0821-185).
   */
  {
    char *outfmt = esl_opt_GetString(go, "--outfmt");
    if (outfmt != NULL) {
      if      (strcmp(outfmt, "1/a") == 0) fmtcode = CM_FILE_1a;
      else if (strcmp(outfmt, "1/e") == 0) fmtcode = CM_FILE_1e;
      else    cm_Fail("No such output format code %s (try 1/a or 1/e).\n", outfmt);
    }
    if (fmtcode == CM_FILE_1a && esl_opt_GetBoolean(go, "--p7pad"))
      cm_Fail("--p7pad can't be used with --outfmt 1/a: the INFERNAL1/a format can't store p7 band pads.\n");
  }

  status = cm_file_Open(cmfile, NULL, TRUE, &cmfp, errbuf); /* TRUE says: allow CM file to be in v1.0 --> v1.0.2 format */
  if      (status == eslENOTFOUND) cm_Fail("File existence/permissions problem in trying to open CM file %s.\n%s\n", cmfile, errbuf);
  else if (status == eslEFORMAT)   cm_Fail("File format problem in trying to open CM file %s.\n%s\n",                cmfile, errbuf);
  else if (status != eslOK)        cm_Fail("Unexpected error %d in opening CM file %s.\n%s\n",                       status, cmfile, errbuf);

  /* --prior <f> (with --upgrade) */
  if (esl_opt_IsOn(go, "--prior")) {
    FILE *pfp;
    if ((pfp = fopen(esl_opt_GetString(go, "--prior"), "r")) == NULL) cm_Fail("Failed to open prior file %s\n", esl_opt_GetString(go, "--prior"));
    if ((user_pri = Prior_Read(pfp)) == NULL)                          cm_Fail("Failed to parse prior file %s\n", esl_opt_GetString(go, "--prior"));
    fclose(pfp);
  }

  /* open output file for writing, if nec */
  if ( esl_opt_IsOn(go, "-o") ) {
    if ((ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) ESL_FAIL(eslFAIL, errbuf, "Failed to open output file %s", esl_opt_GetString(go, "-o"));
  }
  else ofp = stdout;

  while ((status = cm_file_Read(cmfp, TRUE, &abc, &cm)) == eslOK)
    {
      CM_t *cp9cm = NULL;
      int   drops_pknot = FALSE; /* TRUE if the output format can't store this CM's pknot annotation */

      /* The consensus pseudoknot annotation only exists in INFERNAL1/e.
       * Converting to an older format that predates it (1/a via --outfmt,
       * or the legacy v0.7-v1.0.2 -1 format) would lose it, so that's an
       * error unless --dropknots allows it (brief 26_0821-135). Checked
       * before any other work on the CM. */
      if (cm->flags & CMH_PKNOT) {
        if ((esl_opt_GetBoolean(go, "-a") || esl_opt_GetBoolean(go, "-b")) &&
            fmtcode == CM_FILE_1a) {
          drops_pknot = TRUE;
          if (! esl_opt_GetBoolean(go, "--dropknots"))
            cm_Fail("CM %s has consensus pseudoknot annotation, which the INFERNAL1/a format cannot store.\n"
                    "Use --dropknots to convert it anyway without that annotation, or omit --outfmt.",
                    cm->name);
        }
        else if (esl_opt_GetBoolean(go, "-1")) {
          drops_pknot = TRUE;
          if (! esl_opt_GetBoolean(go, "--dropknots"))
            cm_Fail("CM %s has consensus pseudoknot annotation, which the v0.7-v1.0.2 (-1) format cannot store.\n"
                    "Use --dropknots to convert it anyway without that annotation.", cm->name);
        }
      }

      /* --upgrade, or a v1.0-format CM written in the default format:
       * give it derived, calibrated filter HMMs (current format) instead
       * of storing its calibrated ML p7 HMM as a filter HMM block (brief
       * 26_0821-171). 1/e models are already current. */
      int   is_v1p0_default = (cmfp->format == CM_FILE_1 && fmtcode == -1 &&
                               (esl_opt_GetBoolean(go, "-a") || esl_opt_GetBoolean(go, "-b"))) ? TRUE : FALSE;
      int   do_upgrade  = ((esl_opt_GetBoolean(go, "--upgrade") || is_v1p0_default) && ! (cm->flags & CMH_FP7_DERIVED)) ? TRUE : FALSE;

      /* A CM with a filter HMM stored in an older format is written as
       * INFERNAL1/a unless it is upgraded. --outfmt 1/e asks for the
       * current format, which needs --upgrade. --p7pad needs the current
       * format, to store the pads. Both are checked before any work on
       * the CM (brief 26_0821-185). */
      if ((esl_opt_GetBoolean(go, "-a") || esl_opt_GetBoolean(go, "-b")) && ! do_upgrade) {
        if (fmtcode == CM_FILE_1e && ! (cm->flags & CMH_FP7_DERIVED))
          cm_Fail("CM %s is in an older format. Use --upgrade, not --outfmt 1/e, to convert it to the current format.", cm->name);
        if (fmtcode == -1 &&
            cm_file_Check1aStorable(cm, esl_opt_GetBoolean(go, "--p7pad") ? CMH_P7NODEPAD : 0, errbuf) != eslOK)
          cm_Fail(errbuf);
      }

      /* Upgrading a CM built with cmbuild --prior <f> needs that prior
       * to derive its HMM-only search HMM; without it, that HMM can only
       * be the filter HMM (HP7 ml). That's an error unless --noprior
       * allows it, as --dropknots does for pknots (brief 26_0821-179).
       * Checked before any other work on the CM. --prior and --noprior
       * require --upgrade, so a v1.0 CM upgraded by default is told to
       * add it. */
      if (do_upgrade && user_pri == NULL && ! esl_opt_GetBoolean(go, "--noprior") &&
          comlog_cmbuild_option(cm->comlog, "--prior", NULL))
        cm_Fail("CM %s was built with cmbuild --prior, and converting it needs that prior file.\n"
                "Use %s--prior <f> to give it, or %s--noprior to convert it anyway with its filter HMM as its HMM-only search HMM.",
                cm->name,
                esl_opt_GetBoolean(go, "--upgrade") ? "" : "--upgrade ",
                esl_opt_GetBoolean(go, "--upgrade") ? "" : "--upgrade ");

      if(cmfp->format == CM_FILE_1 || esl_opt_GetBoolean(go, "--mlhmm") || esl_opt_GetBoolean(go, "--cp9ascm")) {
	/* if format == CM_FILE_1, we need to calculate QDBs
	 * (cm->dmin, cm->dmax), cm->W, cm->consensus. These are
	 * calculated during model configuration. If --mlhmm, we
	 * need E-value params for the ML p7 HMM, we calc those
	 * in configure_model(). If --cp9ascm, we need cm->cp9 to be
	 * built (globally, i.e. without CM_CONFIG_LOCAL), which
	 * configure_model() also does as a side effect of cm_Configure()
	 * (CP9_2_CM() gives the new CM its own filter HMMs, so the input
	 * CM's filter isn't needed). With --upgrade, the filter HMMs are
	 * derived and calibrated below instead.
	 */
	if ((status = configure_model(cm, errbuf, ! do_upgrade && ! esl_opt_GetBoolean(go, "--cp9ascm"))) != eslOK) cm_Fail(errbuf);
      }

      /* --upgrade: replace the filter HMM stored in an older-format CM
       * with the derived filter and HMM-only HMMs, calibrated as cmbuild
       * would, so the CM is written in the current format (brief
       * 26_0821-165). Before the pads, so that --p7pad computes them
       * for the derived filter HMM. */
      if (do_upgrade) {
        if ((status = upgrade_model(go, cm, cmfp->format, user_pri, errbuf)) != eslOK) cm_Fail(errbuf);
      }

      /* Compute p7 per-node band pads and embed in CM, only with --p7pad and
       * only if the output format supports them (not -1, --mlhmm, --fhmm, --cp9ascm).
       * Pads are opt-in: they cost O(M^2) memory, and the only default reader
       * is cmalign's kmerchain banding above CLEN 4000, which uses a uniform
       * --p7padplus pad without them (briefs 26_0821-155/26_0821-156). Without --p7pad, pads
       * already in the input CM are kept as they are. --no-p7pad is still
       * accepted, as a no-op, for existing scripts.
       */
      if (  esl_opt_GetBoolean(go, "--p7pad")     &&
          ! esl_opt_GetBoolean(go, "-1")           &&
          ! esl_opt_GetBoolean(go, "--mlhmm")      &&
          ! esl_opt_GetBoolean(go, "--fhmm")       &&
          ! esl_opt_GetBoolean(go, "--cp9ascm"))
        {
          if (cm->fp7 == NULL) cm_Fail("CM %s has no filter HMM; cannot compute p7 node pads\n", cm->name);
          /* cm_ComputeP7CMNodePad needs cm->cp9map. For inputs that skipped
           * configure_model() (1/a, 1/e), cp9map is NULL — build just
           * the map cheaply rather than running full cm_Configure.
           */
          if (cm->cp9map == NULL) {
            cm->cp9map = AllocCP9Map(cm);
            CP9_map_cm2hmm(cm, cm->cp9map, 0);
          }
          {
#ifdef HMMER_THREADS
            int pad_ncpu = ESL_MIN(esl_opt_GetInteger(go, "--p7pad-cpu"), esl_threads_GetCPUCount());
#else
            int pad_ncpu = 0;
#endif
            if (pad_ncpu > esl_opt_GetInteger(go, "--p7pad-N")) pad_ncpu = esl_opt_GetInteger(go, "--p7pad-N");
            ESL_RANDOMNESS *pad_r = esl_randomness_Create((uint32_t) esl_opt_GetInteger(go, "--p7pad-seed"));
            if ((status = cm_ComputeP7CMNodePad(cm, pad_r,
                                              esl_opt_GetInteger(go, "--p7pad-N"),
                                              esl_opt_GetReal   (go, "--p7pad-q"),
                                              pad_ncpu,
                                              errbuf)) != eslOK) cm_Fail(errbuf);
            esl_randomness_Destroy(pad_r);
          }
        }

      if (esl_opt_GetBoolean(go, "--cp9ascm")) {
	if ((status = CP9_2_CM(cm, errbuf, &cp9cm)) != eslOK) cm_Fail(errbuf);
      }

      /* append command line info to the appropriate comlog */
      if (esl_opt_GetBoolean(go, "--mlhmm")) {
	if((status = p7_hmm_AppendComlog (cm->mlp7, go->argc, go->argv)) != eslOK) cm_Fail("Failed to record command log");
      }
      else if (esl_opt_GetBoolean(go, "--fhmm")) {
	if((status = p7_hmm_AppendComlog (cm->fp7,  go->argc, go->argv)) != eslOK) cm_Fail("Failed to record command log");
      }
      else if (esl_opt_GetBoolean(go, "--cp9ascm")) {
	if((status = cm_AppendComlog (cp9cm, go->argc, go->argv, FALSE , 0)) != eslOK) cm_Fail("Failed to record command log");
      }
      else {
	if((status = cm_AppendComlog (cm, go->argc, go->argv, FALSE , 0)) != eslOK) cm_Fail("Failed to record command log");
      }

      /* --dropknots, and this CM's pknot annotation can't be written (checked
       * above): drop it silently. For -a/-b targeting 1/a we must
       * actually clear CMH_PKNOT and free cm->pknot, because the binary pknot
       * block is flag-gated (not format-gated) and would otherwise be emitted;
       * the -1 legacy writer cm_file_Write1p0ASCII() ignores cm->pknot. */
      if (drops_pknot && ! esl_opt_GetBoolean(go, "-1")) {
        cm->flags &= ~CMH_PKNOT;
        if (cm->pknot != NULL) { free(cm->pknot); cm->pknot = NULL; }
      }

      if      (esl_opt_GetBoolean(go, "-a")       == TRUE) cm_file_WriteASCII (ofp, fmtcode, cm);
      else if (esl_opt_GetBoolean(go, "-b")       == TRUE) cm_file_WriteBinary(ofp, fmtcode, cm, NULL, NULL);
      else if (esl_opt_GetBoolean(go, "-1")       == TRUE) cm_file_Write1p0ASCII(ofp, cm);
      else if (esl_opt_GetBoolean(go, "--mlhmm")  == TRUE) p7_hmmfile_WriteASCII(ofp, -1, cm->mlp7); /* -1 = write the current default format */
      else if (esl_opt_GetBoolean(go, "--fhmm")   == TRUE) p7_hmmfile_WriteASCII(ofp, -1, cm->fp7);  /* -1 = write the current default format */
      else if (esl_opt_GetBoolean(go, "--cp9ascm")== TRUE) cm_file_WriteASCII(ofp, -1, cp9cm);        /* -1 = write the current default format */

      if (cp9cm != NULL) FreeCM(cp9cm);
      FreeCM(cm);
    }
  if      (status == eslEFORMAT)   cm_Fail("bad file format in CM file %s\n%s",             cmfile, cmfp->errbuf);
  else if (status == eslEINCOMPAT) cm_Fail("CM file %s contains different alphabets\n%s",   cmfile, cmfp->errbuf);
  else if (status != eslEOF)       cm_Fail("Unexpected error in reading CMs from %s\n%s",   cmfile, cmfp->errbuf);

  cm_file_Close(cmfp);

  if(esl_opt_IsOn(go, "-o")) fclose(ofp);
  if(user_pri != NULL) Prior_Destroy(user_pri);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}


/* configure_model()
 * Configure the model. This determines QDBs and W, which
 * the new file format includes, but v1.0-->v1.0.2 did not,
 * thus we have to calculate them for models read from
 * v1.0-->v1.0.2 cm files.
 */
static int
configure_model(CM_t *cm, char *errbuf, int set_filter)
{
  int status;
  int lmsvL, lvitL, lfwdL, gfwdL;
  int lmsvN, lvitN, lfwdN, gfwdN;
  float lftailp, gftailp;
  double fil_gfmu, fil_gflambda;

  /* Configure the model, we must calculate QDBs so we can write them to the CM file */
  cm->config_opts |= CM_CONFIG_QDB;
  if ((status = cm_Configure(cm, errbuf, -1)) != eslOK) return status;
  if ((status = cm_SetConsensus(cm, cm->cmcons, NULL)) != eslOK) ESL_FAIL(status, errbuf, "Failed to calculate consensus sequence");
  if (! set_filter) return eslOK;

  /* We'll define the filter HMM as the ML p7 HMM because that's the
   * only option available (by default, in cmbuild, a filter HMM gets
   * built separately that's different from the ml p7, but that
   * requires the cmbuild input alignment).  The cm->mlp7 HMM was
   * created in cm_Configure(), and we calibrate it here. There are
   * more options than this in cmbuild, but here we enforce
   * defaults. See cmbuild.c::build_and_calibrate_p7_filter(). */
  lmsvL = lvitL = 200;
  lfwdL = 100;
  gfwdL = ESL_MAX(100, 2.*cm->clen);
  lmsvN = lvitN = lfwdN = gfwdN = 200;
  lftailp = 0.055;
  gftailp = 0.065;

  /* Calibrate the ML p7 hmm */
  if((status = cm_p7_Calibrate(cm->mlp7, errbuf,
			       lmsvL, lvitL, lfwdL, gfwdL, /* length of sequences to search for local (lL) and glocal (gL) modes */
			       lmsvN, lvitN, lfwdN, gfwdN, /* number of seqs to search for each alg */
			       lftailp,                    /* fraction of tail mass to fit for local Fwd */
			       gftailp,                    /* fraction of tail mass to fit for glocal Fwd */
			       42,                         /* seed: default */
			       0,                          /* ncpus: 0 = serial */
			       cm_p7_UseP7mlPredictor(cm, TRUE), /* filter is the ML HMM (brief 26_0824-085) */
			       FALSE,                      /* do_sim: predict, don't simulate, glocal Fwd stats (brief 26_0824-087) */
			       &fil_gfmu, &fil_gflambda))
     != eslOK) return status;
  if((status = cm_SetFilterHMM(cm, cm->mlp7, fil_gfmu, fil_gflambda)) != eslOK) ESL_FAIL(status, errbuf, "Unable to set the HMM filter for the CM");

  return eslOK;
}


/* upgrade_model()
 * Incept:    EPN, Fri Oct 30 2026 (w/Claude)
 *
 * Purpose:   cmconvert --upgrade (brief 26_0821-165): give <cm>, read
 *            from a CM file in older format <format>, the filter HMMs
 *            that cmbuild gives a CM now, so that it is written in the
 *            current format: discard its stored filter HMM, derive the
 *            filter HMM (the CM's ML p7 HMM) and the HMM-only search
 *            HMM from the CM, and calibrate both, with
 *            cm_p7_CalibrateDerivedFilterHMMs() and cmbuild's default
 *            calibration settings.
 *
 *            The HMM-only HMM depends on how the CM was built, which a
 *            CM file doesn't record except in its COM lines. An old-
 *            format file was written by a cmbuild older than this
 *            one, so these rules describe how those cmbuilds built a
 *            CM, not how cmbuild builds one now; they are kept here,
 *            not shared with cmbuild:
 *            - The CM's transition prior: <user_pri> (cmconvert
 *              --prior) if given; else the v0.56-v1.0.2 prior for a
 *              v1.0 format file, or if the cmbuild COM line has --p56
 *              or --v1p0; else the default prior. If the COM line has
 *              --prior <f> and <user_pri> isn't given, the prior is
 *              unknown and the HMM-only HMM is the filter HMM (HP7 ml);
 *              the caller allows that only with cmconvert --noprior.
 *              (A recorded path isn't trusted to be readable, or to be
 *              the same file.)
 *            - A CM with no basepairs built with the default zero-
 *              basepair (HMMER-like) prior, i.e. not from a v1.0 file
 *              and without --p56/--v1p0/--noh3pri/--prior on its
 *              cmbuild COM line, nor cmconvert --prior: HP7 ml.
 *              cmbuild before 1.2 (and 1.2 development versions
 *              before brief 26_0821-182) built such a CM as if it
 *              were an HMM, with an HMMER3-like prior, and searched it
 *              with its ML HMM. cmbuild now builds it like any other
 *              CM, and writes INFERNAL1/e, with an HP7 line, so its
 *              files never reach this function. (A
 *              1.2 0-bp CM written in an old format with cmconvert
 *              --outfmt and then upgraded gets HP7 ml, because the
 *              old format can't say which cmbuild built it.)
 *            - The filter HMMs' glocal Forward E-value parameters are
 *              simulated, not predicted, where cmbuild would simulate
 *              them: --Esim, a --prior build, --eset 0, or an EFFN at
 *              EFORCE_MAX_NEFF (a missed --eforce target).
 *            - The calibration settings are cmbuild's defaults,
 *              changed by any of cmbuild's calibration options on its
 *              COM line (--ElL, --EgL, --Elcmult, --Egcmult, --EmN,
 *              --EvN, --ElfN, --EgfN, --Elftp, --Egftp, --Efitlam,
 *              --Eseed), as cmbuild applies them. (Model libraries of
 *              long CMs are often built with small --EmN etc.)
 *            The CM's own E-value parameters (ECM lines, from
 *            cmcalibrate) belong to the CM and are kept.
 *
 * Returns:   eslOK on success, errbuf filled on failure.
 */
static int
upgrade_model(const ESL_GETOPTS *go, CM_t *cm, int format, const Prior_t *user_pri, char *errbuf)
{
  int            status;
  Prior_t       *pri       = NULL;
  const Prior_t *pri2use   = NULL;
  char          *eset_arg  = NULL;
  int            com_prior = comlog_cmbuild_option(cm->comlog, "--prior", NULL);
  int            is_v1p0   = (format == CM_FILE_1 || comlog_cmbuild_option(cm->comlog, "--p56", NULL) || comlog_cmbuild_option(cm->comlog, "--v1p0", NULL)) ? TRUE : FALSE;
  int            built_as_hmm;      /* TRUE if a pre-1.2 cmbuild built this 0-bp CM as an HMM */
  int            do_sim;
  int            ncpus     = 0;
  double         x;
  CM_P7_CALOPTS  co;

  built_as_hmm = (CMCountNodetype(cm, MATP_nd) == 0 && user_pri == NULL && ! com_prior && ! is_v1p0 &&
                  ! comlog_cmbuild_option(cm->comlog, "--noh3pri", NULL)) ? TRUE : FALSE;

  if      (user_pri != NULL) pri2use = user_pri;
  else if (com_prior)        pri2use = NULL;  /* unknown prior file */
  else if (is_v1p0)          pri2use = pri = Prior_Default_v0p56_through_v1p02();
  else                       pri2use = pri = Prior_Default();

  if      (built_as_hmm)                                           cm->hp7_xform = CM_HP7_ML;
  else if (pri2use != NULL && Prior_HMMOnlyAlpha(pri2use, cm->hp7_alpha)) cm->hp7_xform = CM_HP7_H3PRI;
  else                                                             cm->hp7_xform = CM_HP7_ML;
  if (cm->hp7_xform == CM_HP7_ML) cm_p7_DefaultHMMOnlyAlpha(cm->hp7_alpha); /* unused; keep the in-memory default */

  do_sim = (user_pri != NULL || com_prior || cm->eff_nseq >= EFORCE_MAX_NEFF || comlog_cmbuild_option(cm->comlog, "--Esim", NULL)) ? TRUE : FALSE;
  if (comlog_cmbuild_option(cm->comlog, "--eset", &eset_arg) && eset_arg != NULL && atof(eset_arg) == 0.) do_sim = TRUE;
  if (eset_arg != NULL) free(eset_arg);

  /* calibration settings: cmbuild's, as in its build_and_calibrate_p7_filter() */
  cm_p7_CalOptsSetDefaults(&co, cm->clen, do_sim);
  if      (comlog_cmbuild_number(cm->comlog, "--ElL",     &x)) co.lmsvL = co.lvitL = co.lfwdL = (int) x;
  else if (comlog_cmbuild_number(cm->comlog, "--Elcmult", &x)) co.lmsvL = co.lvitL = co.lfwdL = x * cm->clen;
  if      (comlog_cmbuild_number(cm->comlog, "--EgL",     &x)) co.gfwdL = (int) x;
  else if (comlog_cmbuild_number(cm->comlog, "--Egcmult", &x)) co.gfwdL = x * cm->clen;
  if (comlog_cmbuild_option(cm->comlog, "--Efitlam", NULL)) {
    co.lmsvN = co.lvitN = co.lfwdN = 10000;
    co.lftailp = 0.01;
    if (! do_sim) { co.gfwdN = 10000; co.gftailp = 0.01; }
  }
  if (comlog_cmbuild_number(cm->comlog, "--EmN",   &x)) co.lmsvN   = (int) x;
  if (comlog_cmbuild_number(cm->comlog, "--EvN",   &x)) co.lvitN   = (int) x;
  if (comlog_cmbuild_number(cm->comlog, "--ElfN",  &x)) co.lfwdN   = (int) x;
  if (comlog_cmbuild_number(cm->comlog, "--EgfN",  &x)) co.gfwdN   = (int) x;
  if (comlog_cmbuild_number(cm->comlog, "--Elftp", &x)) co.lftailp = x;
  if (comlog_cmbuild_number(cm->comlog, "--Egftp", &x)) co.gftailp = x;
  if (comlog_cmbuild_number(cm->comlog, "--Eseed", &x)) co.seed    = (int) x;
#ifdef HMMER_THREADS
  ncpus = ESL_MIN(esl_opt_GetInteger(go, "--cpu"), esl_threads_GetCPUCount());
#endif
  if ((status = cm_p7_CalibrateDerivedFilterHMMs(cm, &co, ncpus, errbuf)) != eslOK) goto ERROR;

  if (pri != NULL) Prior_Destroy(pri);
  return eslOK;

 ERROR:
  if (pri != NULL) Prior_Destroy(pri);
  return status;
}

/* comlog_cmbuild_option()
 * Incept:    EPN, Fri Oct 30 2026 (w/Claude)
 *
 * Purpose:   Return TRUE if the first cmbuild command in a CM's
 *            command log <comlog> (its COM lines) has option <opt>,
 *            and if <opt_arg> isn't NULL, return the token after it
 *            in <*opt_arg> (newly allocated; NULL if none). Return
 *            FALSE if <comlog> is NULL or has no cmbuild command.
 *            Used by cmconvert --upgrade (brief 26_0821-165).
 */
static int
comlog_cmbuild_option(const char *comlog, const char *opt, char **opt_arg)
{
  char *lines = NULL;
  char *line, *lp;
  char *tok,  *tp;
  char *prog;
  int   found = FALSE;
  int   done  = FALSE;

  if (opt_arg != NULL) *opt_arg = NULL;
  if (comlog == NULL) return FALSE;
  if (esl_strdup(comlog, -1, &lines) != eslOK) return FALSE;

  for (line = strtok_r(lines, "\n", &lp); line != NULL && ! done; line = strtok_r(NULL, "\n", &lp)) {
    /* the program is the first token, possibly with a path */
    if ((tok = strtok_r(line, " \t", &tp)) == NULL) continue;
    prog = strrchr(tok, '/');
    prog = (prog == NULL) ? tok : prog+1;
    if (strcmp(prog, "cmbuild") != 0) continue;
    done = TRUE; /* only the first cmbuild command */
    while ((tok = strtok_r(NULL, " \t", &tp)) != NULL) {
      if (strcmp(tok, opt) == 0) {
        found = TRUE;
        if (opt_arg != NULL && (tok = strtok_r(NULL, " \t", &tp)) != NULL) esl_strdup(tok, -1, opt_arg);
        break;
      }
      if (strncmp(tok, opt, strlen(opt)) == 0 && tok[strlen(opt)] == '=') { /* --opt=arg */
        found = TRUE;
        if (opt_arg != NULL) esl_strdup(tok + strlen(opt) + 1, -1, opt_arg);
        break;
      }
    }
  }
  free(lines);
  return found;
}

/* comlog_cmbuild_number()
 * Incept:    EPN, Fri Oct 30 2026 (w/Claude)
 *
 * Purpose:   As comlog_cmbuild_option(), for an option with a numeric
 *            argument: return TRUE and the argument in <*ret_x> if
 *            the CM's first cmbuild command has option <opt> with a
 *            number after it.
 */
static int
comlog_cmbuild_number(const char *comlog, const char *opt, double *ret_x)
{
  char *arg   = NULL;
  int   found = FALSE;

  if (comlog_cmbuild_option(comlog, opt, &arg) && arg != NULL && esl_str_IsReal(arg)) {
    *ret_x = atof(arg);
    found  = TRUE;
  }
  if (arg != NULL) free(arg);
  return found;
}

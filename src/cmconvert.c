/* cmconvert: converting covariance model files to Infernal-1.1 CM format.
 * 
 * EPN, Fri Jul  1 05:11:15 2011
 * SRE, Thu Oct 16 08:57:43 2008 [janelia] (hmmconvert.c)
 * SVN $Id$
 */
#include "esl_config.h"
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"

#include "hmmer.h"
#include "funcs.h"		/* external functions                   */
#include "structs.h"		/* data structures, macros, #define's   */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,       NULL,    NULL, "show brief help on version and usage",                             0 },
  { "-a",        eslARG_NONE,"default",NULL, NULL, "-a,-b,-0,-1",   NULL,    NULL, "ascii:  output models in INFERNAL 1.1 ASCII format",               0 },
  { "-b",        eslARG_NONE,   FALSE, NULL, NULL, "-a,-b,-0,-1",   NULL,    NULL, "binary: output models in INFERNAL 1.1 binary format",              0 },
  { "-1",        eslARG_NONE,   FALSE, NULL, NULL, "-a,-b,-0,-1",   NULL,    NULL, "output backward compatible Infernal v1.0-->v1.0.2 ASCII format",   0 },
  { "-0",        eslARG_NONE,   FALSE, NULL, NULL, "-a,-b,-0,-1",   NULL,    NULL, "output backward compatible Infernal v0.1-->v0.81  ASCII format",   0 },
  { "--outfmt",  eslARG_STRING, NULL,  NULL, NULL,      NULL,       NULL, "-0,-1", "choose output legacy 1.x file formats by name, such as '1/a'",     0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <cmfile>";
static char banner[] = "convert CM file to a different Infernal format";

static int  configure_model(CM_t *cm, char *errbuf);

int 
main(int argc, char **argv)
{
  ESL_GETOPTS   *go      = cm_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  ESL_ALPHABET  *abc     = NULL;
  char          *cmfile  = esl_opt_GetArg(go, 1);
  CM_FILE       *cmfp    = NULL;
  CM_t          *cm      = NULL;
  FILE          *ofp     = stdout;
  char          *outfmt  = esl_opt_GetString(go, "--outfmt");
  int            fmtcode = -1;	/* -1 = write the current default format */
  int            status;
  char           errbuf[eslERRBUFSIZE];

  if (outfmt != NULL) {
    if      (strcmp(outfmt, "1/a") == 0) fmtcode = CM_FILE_1a;
    else    cm_Fail("No such 1.x output format code %s.\n", outfmt);
  }

  status = cm_file_Open(cmfile, NULL, TRUE, &cmfp, errbuf); /* TRUE says: allow CM file to be in v1.0 --> v1.0.2 format */
  if      (status == eslENOTFOUND) cm_Fail("File existence/permissions problem in trying to open CM file %s.\n%s\n", cmfile, errbuf);
  else if (status == eslEFORMAT)   cm_Fail("File format problem in trying to open CM file %s.\n%s\n",                cmfile, errbuf);
  else if (status != eslOK)        cm_Fail("Unexpected error %d in opening CM file %s.\n%s\n",                       status, cmfile, errbuf);  

  while ((status = cm_file_Read(cmfp, TRUE, &abc, &cm)) == eslOK)
    {
      if(cmfp->format == CM_FILE_1) { 
	/* we need to calculate QDBs (cm->dmin, cm->dmax), cm->W, cm->consensus, 
	 * as well as E-value parameters for the ML p7 HMM  */
	if ((status = configure_model(cm, errbuf))       != eslOK) cm_Fail(errbuf);
      }	
      if      (esl_opt_GetBoolean(go, "-a") == TRUE) cm_file_WriteASCII (ofp, fmtcode, cm);
      else if (esl_opt_GetBoolean(go, "-b") == TRUE) cm_file_WriteBinary(ofp, fmtcode, cm, NULL);

      FreeCM(cm);
    }
  if      (status == eslEFORMAT)   cm_Fail("bad file format in CM file %s\n%s",             cmfile, cmfp->errbuf);
  else if (status == eslEINCOMPAT) cm_Fail("CM file %s contains different alphabets\n%s",   cmfile, cmfp->errbuf);
  else if (status != eslEOF)       cm_Fail("Unexpected error in reading CMs from %s\n%s",   cmfile, cmfp->errbuf);

  cm_file_Close(cmfp);
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
configure_model(CM_t *cm, char *errbuf)
{
  int status; 
  CMConsensus_t *cons    = NULL;
  int lmsvL, lvitL, lfwdL, gfwdL;
  int lmsvN, lvitN, lfwdN, gfwdN;
  float lftailp, gftailp;
  double fil_gfmu, fil_gflambda;

  /* Configure the model, we must calculate QDBs so we can write them to the CM file */
  cm->config_opts |= CM_CONFIG_QDB;   
  if((status = cm_Configure(cm, errbuf, -1)) != eslOK) return status;

  CreateCMConsensus(cm, cm->abc, 3.0, 1.0, &(cons));
  if ((status = cm_SetConsensus  (cm, cons, NULL)) != eslOK) ESL_FAIL(status, errbuf, "Failed to calculate consensus sequence");
  FreeCMConsensus(cons);

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
			       &fil_gfmu, &fil_gflambda))  
     != eslOK) return status;
  if((status = cm_SetFilterHMM(cm, cm->mlp7, fil_gfmu, fil_gflambda)) != eslOK) ESL_FAIL(status, errbuf, "Unable to set the HMM filter for the CM");

  return eslOK;
}
/*****************************************************************
 * @LICENSE@
 *****************************************************************/

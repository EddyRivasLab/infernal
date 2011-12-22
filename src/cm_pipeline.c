/* Infernal's accelerated seq/profile comparison pipeline
 *  
 * Contents:
 *   1. CM_PIPELINE: allocation, initialization, destruction
 *   2. Pipeline API
 *   TODO  3. Example 1: search mode (in a sequence db)
 *   TODO: 4. Example 2: scan mode (in an HMM db)
 *   5. Copyright and license information
 * 
 * SRE, Fri Dec  5 10:09:39 2008 [Janelia] [BSG3, Bear McCreary]
 * SVN $Id: p7_pipeline.c 3352 2010-08-24 20:56:01Z eddys $
 */
#include "esl_config.h"
#include "p7_config.h"
#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h> 

#include "easel.h"
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "funcs.h"
#include "structs.h"

#define DOPRINT  1
#define DOPRINT2 0
#define DOPRINT3 0

static int  merge_windows_from_two_lists(int64_t *ws1, int64_t *we1, double *wp1, int *wl1, int nwin1, int64_t *ws2, int64_t *we2, double *wp2, int *wl2, int nwin2, int64_t **ret_mws, int64_t **ret_mwe, double **ret_mwp, int **ret_mwl, int *ret_nmwin);
static void copy_subseq(const ESL_SQ *src_sq, ESL_SQ *dest_sq, int64_t i, int64_t L);

/*****************************************************************
 * 1. The CM_PIPELINE object: allocation, initialization, destruction.
 *****************************************************************/

/* Function:  cm_pipeline_Create()
 * Synopsis:  Create a new accelerated comparison pipeline.
 * Incept:    EPN, Fri Sep 24 16:14:39 2010
 *            SRE, Fri Dec  5 10:11:31 2008 [Janelia] (p7_pipeline_Create())
 *
 * Purpose:   Given an application configuration structure <go>
 *            containing certain standardized options (described
 *            below), some initial guesses at the model size <M_hint>
 *            and sequence length <L_hint> that will be processed,
 *            and a <mode> that can be either <cm_SCAN_MODELS> or
 *            <cm_SEARCH_SEQS> depending on whether we're searching one sequence
 *            against a model database (cmscan mode) or one model
 *            against a sequence database (cmsearch mode); create new
 *            pipeline object.
 *
 *            In search mode, we would generally know the length of
 *            our query profile exactly, and would pass <cm->clen> as <M_hint>;
 *            in scan mode, we generally know the length of our query
 *            sequence exactly, and would pass <sq->n> as <L_hint>.
 *            Targets will come in various sizes as we read them,
 *            and the pipeline will resize any necessary objects as
 *            needed, so the other (unknown) length is only an
 *            initial allocation.
 *
 *            <Z> is passed as the database size, in residues, if
 *            known. If unknown, 0 should be passed as <Z>.
 *            
 *            The configuration <go> must include settings for the 
 *            following options:
 *            
 *            || option      ||            description                     || usually  ||
 *            | --noali      |  don't output alignments (smaller output)    |   FALSE   |
 *            | -E           |  report hits <= this E-value threshold       |    10.0   |
 *            | -T           |  report hits >= this bit score threshold     |    NULL   |
 *            | -Z           |  set initial hit search space size           |    NULL   |
 *            | --incE       |  include hits <= this E-value threshold      |    0.01   |
 *            | --incT       |  include hits >= this bit score threshold    |    NULL   |
 *            | --cut_ga     |  model-specific thresholding using GA        |   FALSE   |
 *            | --cut_nc     |  model-specific thresholding using NC        |   FALSE   |
 *            | --cut_tc     |  model-specific thresholding using TC        |   FALSE   |
 *            | --max        |  turn all heuristic filters off              |   FALSE   |
 *            | --fZ <x>     |  set filter thr as if dbsize were <x> Mb     |     OFF   |
 *            | --F1         |  Stage 1 (MSV) thresh: promote hits P <= F1  |    0.35   |
 *            | --F2         |  Stage 2 (Vit) thresh: promote hits P <= F2  |    0.10   |
 *            | --F3         |  Stage 3 (Fwd) thresh: promote hits P <= F3  |    0.02   |
 *            | --null2      |  turn ON biased comp score correction        |   FALSE   |
 *            | --seed       |  RNG seed (0=use arbitrary seed)             |     181   |
 *            | --acc        |  prefer accessions over names in output      |   FALSE   |
 * *** opts below this line are unique to cm_pipeline.c ***
 *            | --nonull3    |  turn off NULL3 correction                   |   FALSE   |
 *            | -g           |  configure the CM for glocal alignment       |   FALSE   |
 *            | --domsvbias  |  turn ON composition bias filter HMM for MSV |   FALSE   |
 *            | --novitbias  |  turn OFF composition bias filter HMM for Vit |   FALSE   |
 *            | --nofwdbias  |  turn OFF composition bias filter HMM for Fwd |   FALSE   |
 *            | --nogfwdbias |  turn OFF composition bias filter HMM for gFwd |   FALSE   |
 *            | --noedefbias  |  turn OFF composition bias filter HMM for ddef| FALSE   |
 *            | --F1         |  Stage 1 (MSV) filter thresh                 |   0.40    |
 *            | --F1b        |  Stage 1 (MSV) bias filter thresh            |   OFF     |
 *            | --F2         |  Stage 2 (Vit) filter thresh                 |   0.15    |
 *            | --F2b        |  Stage 2 (Vit) bias filter thresh            |   0.15    |
 *            | --F3         |  Stage 3 (Fwd) filter thresh                 |   0.02    |
 *            | --F3b        |  Stage 3 (Fwd) bias filter thresh            |   0.02    |
 *            | --F4         |  Stage 4 (gFwd) filter thresh                |   0.02    |
 *            | --F4b        |  Stage 4 (gFwd) bias filter thresh           |   0.02    |
 *            | --F5         |  Stage 5 (env def) filter thresh             |   0.02    |
 *            | --F5b        |  Stage 5 (env def) bias filter thresh        |   0.02    |
 *            | --F6         |  Stage 5 (CYK) filter thresh                 |   0.0005  |
 *            | --fast       |  set filters at strict-level                 |   FALSE   |
 *            | --mid        |  set filters at mid-level                    |   FALSE   |
 *            | --cyk        |  set final search stage as CYK, not inside   |   FALSE   |
 *            | --hmm        |  don't use CM at all, HMM only               |   FALSE   |
 *            | --noenvdef   |  don't define envelopes before CM stages     |   FALSE   |
 *            | --nomsv      |  skip MSV filter stage                       |   FALSE   |
 *            | --shortmsv   |  use standard MSV, not LongTarget variant    |   FALSE   |
 *            | --novit      |  skip Viterbi filter stage                   |   FALSE   |
 *            | --nofwd      |  skip Forward filter stage                   |   FALSE   |
 *            | --nohmm      |  skip all HMM filter stages                  |   FALSE   |
 *            | --nocyk      |  don't filter with CYK, skip to final stage  |   FALSE   |
 *            | --fnoqdb     |  don't use QDBs in CYK filter                |   FALSE   |
 *            | --fbeta      |  set beta for QDBs in CYK filter to this     |   1e-9    |
 *            | --fhbanded   |  use HMM bands for CYK filter                |   FALSE   |
 *            | --ftau       |  set tau for HMM bands in CYK filter to this |   5e-6    |
 *            | --fsums      |  use sums to get CYK filter HMM bands        |   FALSE   |
 *            | --tau        |  set tau for final HMM bands to this         |   FALSE   |
 *            | --sums       |  use sums to get final HMM bands             |   FALSE   |
 *            | --qdb        |  use QDBs, not HMM bands in final stage      |   FALSE   |
 *            | --beta       |  set beta for QDBs in final stage            |   1e-15   |
 *            | --nonbanded  |  don't use QDBs or HMM bands in final stage  |   FALSE   |
 *            | --rt1        |  set envelope def rt1 parameter as <x>         |   0.25    |
 *            | --rt2        |  set envelope def rt2 parameter as <x>         |   0.10    |
 *            | --rt3        |  set envelope def rt3 parameter as <x>         |   0.20    |
 *            | --localweak  |  rescore low-scoring envelopes in local mode   |   FALSE   |
 *            | --localenv   |  define envelopes in local mode                |   FALSE   |
 *            | --ns         |  set number of samples for envelope traceback  |    200    |
 *            | --wsplit     |  split windows after MSV                     |   FALSE   |
 *            | --wmult      |  multiplier for splitting windows if --wsplit|   3.0     |
 * Returns:   ptr to new <cm_PIPELINE> object on success. Caller frees this
 *            with <cm_pipeline_Destroy()>.
 *
 * Throws:    <NULL> on allocation failure.
 */
CM_PIPELINE *
cm_pipeline_Create(ESL_GETOPTS *go, ESL_ALPHABET *abc, int clen_hint, int L_hint, int64_t Z, enum cm_zsetby_e Z_setby, enum cm_pipemodes_e mode)
{
  CM_PIPELINE *pli  = NULL;
  int          seed = esl_opt_GetInteger(go, "--seed");
  int          status;
  double       Z_Mb; /* database size in Mb */
  int          pass_idx; /* counter over passes */

  ESL_ALLOC(pli, sizeof(CM_PIPELINE));

  /* allocate matrices */
  if ((pli->fwd  = p7_omx_Create(clen_hint, L_hint, L_hint)) == NULL) goto ERROR;
  if ((pli->bck  = p7_omx_Create(clen_hint, L_hint, L_hint)) == NULL) goto ERROR;
  if ((pli->oxf  = p7_omx_Create(clen_hint, 0,      L_hint)) == NULL) goto ERROR;
  if ((pli->oxb  = p7_omx_Create(clen_hint, 0,      L_hint)) == NULL) goto ERROR;     
  if ((pli->gfwd = p7_gmx_Create(clen_hint, L_hint))         == NULL) goto ERROR;
  if ((pli->gbck = p7_gmx_Create(clen_hint, L_hint))         == NULL) goto ERROR;
  if ((pli->gxf  = p7_gmx_Create(clen_hint, L_hint))         == NULL) goto ERROR;
  if ((pli->gxb  = p7_gmx_Create(clen_hint, L_hint))         == NULL) goto ERROR;     

  /* intialize model-dependent parameters, these are invalid until we call cm_pli_NewModel() */
  pli->maxW = 0;
  pli->cmW  = 0;
  pli->clen = 0;

  /* Normally, we reinitialize the RNG to the original seed every time we're
   * about to collect a stochastic trace ensemble. This eliminates run-to-run
   * variability. As a special case, if seed==0, we choose an arbitrary one-time 
   * seed: time() sets the seed, and we turn off the reinitialization.
   */
  pli->r                  =  esl_randomness_CreateFast(seed);
  pli->do_reseeding       = (seed == 0) ? FALSE : TRUE;
  pli->ddef               = p7_domaindef_Create(pli->r);
  pli->ddef->do_reseeding = pli->do_reseeding;

  /* Configure reporting thresholds */
  pli->by_E            = TRUE;
  pli->E               = esl_opt_GetReal(go, "-E");
  pli->T               = 0.0;
  pli->use_bit_cutoffs = FALSE;
  if (esl_opt_IsOn(go, "-T")) 
    {
      pli->T    = esl_opt_GetReal(go, "-T"); 
      pli->by_E = FALSE;
    } 

  /* Configure inclusion thresholds */
  pli->inc_by_E           = TRUE;
  pli->incE               = esl_opt_GetReal(go, "--incE");
  pli->incT               = 0.0;
  if (esl_opt_IsOn(go, "--incT")) 
    {
      pli->incT     = esl_opt_GetReal(go, "--incT"); 
      pli->inc_by_E = FALSE;
    } 

  /* Configure for one of the model-specific thresholding options */
  if (esl_opt_GetBoolean(go, "--cut_ga"))
    {
      pli->T        = 0.0;
      pli->by_E     = FALSE;
      pli->incT     = 0.0;
      pli->inc_by_E = FALSE;
      pli->use_bit_cutoffs = CMH_GA;
    }
  if (esl_opt_GetBoolean(go, "--cut_nc"))
    {
      pli->T        = 0.0;
      pli->by_E     = FALSE;
      pli->incT     = 0.0;
      pli->inc_by_E = FALSE;
      pli->use_bit_cutoffs = CMH_NC;
    }
  if (esl_opt_GetBoolean(go, "--cut_tc"))
    {
      pli->T        = 0.0;
      pli->by_E     = FALSE;
      pli->incT     = 0.0;
      pli->inc_by_E = FALSE;
      pli->use_bit_cutoffs = CMH_TC;
    }

  /* Configure envelope definition parameters */
  pli->rt1 = esl_opt_GetReal(go, "--rt1");
  pli->rt2 = esl_opt_GetReal(go, "--rt2");
  pli->rt3 = esl_opt_GetReal(go, "--rt3");
  pli->ns  = esl_opt_GetInteger(go, "--ns");
  pli->ddef->rt1 = pli->rt1;
  pli->ddef->rt2 = pli->rt2;
  pli->ddef->rt3 = pli->rt3;
  pli->ddef->nsamples = pli->ns;
  pli->do_localenv      = esl_opt_GetBoolean(go, "--localenv");
  pli->do_wsplit        = (esl_opt_GetBoolean(go, "--wnosplit")) ? FALSE : TRUE;
  pli->wmult            = esl_opt_GetReal(go, "--wmult");
  pli->do_wcorr         = esl_opt_GetBoolean(go, "--wcorr");
  pli->do_nocorr        = esl_opt_GetBoolean(go, "--nocorr");
  pli->do_oldcorr       = esl_opt_GetBoolean(go, "--oldcorr");
  pli->do_envwinbias    = (! esl_opt_GetBoolean(go, "--envhitbias"));
  pli->do_filcmW        = esl_opt_GetBoolean(go, "--filcmW");
  pli->do_glen          = esl_opt_GetBoolean(go, "--glen");
  pli->glen_min         = esl_opt_GetInteger(go, "--glN");
  pli->glen_max         = esl_opt_GetInteger(go, "--glX");
  pli->glen_step        = esl_opt_GetInteger(go, "--glstep");
  pli->research_ends    = (esl_opt_GetBoolean(go, "--noends"))  ? FALSE : TRUE;
  pli->do_trunc_ends    = (esl_opt_GetBoolean(go, "--notrunc")) ? FALSE : TRUE;
  pli->do_force_ends    = (esl_opt_GetBoolean(go, "--noforce")) ? FALSE : TRUE;
  pli->do_local_ends    = (esl_opt_GetBoolean(go, "--locends")) ? TRUE  : FALSE;
  pli->xtau             = esl_opt_GetReal(go, "--xtau");

  /* Initialize per-sequence information */
  pli->rs5term = FALSE;
  pli->rs3term = FALSE;

  pli->tro = CreateTruncOpts();

  /* Configure acceleration pipeline thresholds */
  pli->do_cm         = TRUE;
  pli->do_hmm        = TRUE;
  pli->do_max        = FALSE;
  pli->do_rfam       = FALSE;
  pli->do_envelopes  = TRUE;
  pli->do_pad        = esl_opt_GetBoolean(go, "--pad");
  pli->do_msv        = TRUE;
  pli->do_shortmsv   = FALSE;
  pli->do_msvbias    = FALSE;
  pli->do_vit        = TRUE;
  pli->do_vitbias    = TRUE;
  pli->do_fwd        = TRUE;
  pli->do_fwdbias    = TRUE;
  pli->do_gfwd       = TRUE;
  pli->do_gfwdbias   = TRUE;
  pli->do_edefbias   = TRUE;
  pli->do_cyk        = TRUE;
  pli->do_null2      = FALSE;
  pli->do_null3      = TRUE;

  /* Configure search space sizes for E value calculations and filter thresholds.
   * This is database size dependent. Which was passed in. If -Z <x> enabled, we overwrite
   * the passed in value with <x>.
   */
  if (esl_opt_IsUsed(go, "-Z")) { 
      pli->Z_setby = CM_ZSETBY_OPTION;
      pli->Z       = (int64_t) (esl_opt_GetReal(go, "-Z") * 1000000.); 
  }
  else { 
    pli->Z       = Z;       /* Z is an input variable to this function */
    pli->Z_setby = Z_setby; /* so is Z_setby */
  }
  /* set filter thresholds dependent on Z */
  if(esl_opt_IsUsed(go, "--fZ")) { 
    Z_Mb = esl_opt_GetReal(go, "--fZ"); 
  }
  else { 
    Z_Mb = pli->Z / 1000000.;
  }
  if(Z_Mb >= (100000. - eslSMALLX1)) { /* Z >= 100 Gb */
    pli->F1 = pli->F1b = 0.05;
    pli->F2 = pli->F2b = 0.04;
    pli->F3 = pli->F3b = 0.0004;
    pli->F4 = pli->F4b = 0.0004;
    pli->F5 = pli->F5b = 0.0004;
    pli->F6 = 0.0001;
  }
  else if(Z_Mb >= (10000. - eslSMALLX1)) { /* 100 Gb > Z >= 10 Gb */
    pli->F1 = pli->F1b = 0.06;
    pli->F2 = pli->F2b = 0.05;
    pli->F3 = pli->F3b = 0.0005;
    pli->F4 = pli->F4b = 0.0005;
    pli->F5 = pli->F5b = 0.0005;
    pli->F6 = 0.0001;
  }
  else if(Z_Mb >= (1000. - eslSMALLX1)) { /* 10 Gb > Z >= 1 Gb */
    pli->F1 = pli->F1b = 0.06;
    pli->F2 = pli->F2b = 0.15;
    pli->F3 = pli->F3b = 0.0005;
    pli->F4 = pli->F4b = 0.0005;
    pli->F5 = pli->F5b = 0.0005;
    pli->F6 = 0.0001;
  }
  else if(Z_Mb >= (100. - eslSMALLX1)) { /* 1 Gb  > Z >= 100 Mb */
    pli->F1 = pli->F1b = 0.30;
    pli->F2 = pli->F2b = 0.15;
    pli->F3 = pli->F3b = 0.002;
    pli->F4 = pli->F4b = 0.002;
    pli->F5 = pli->F5b = 0.002;
    pli->F6 = 0.0001;
  }
  else if(Z_Mb >= (10. - eslSMALLX1)) { /* 100 Mb  > Z >= 10 Mb */
    pli->F1 = pli->F1b = 0.35;
    pli->F2 = pli->F2b = 0.20;
    pli->F3 = pli->F3b = 0.003;
    pli->F4 = pli->F4b = 0.003;
    pli->F5 = pli->F5b = 0.003;
    pli->F6 = 0.0001;
  }
  else if(Z_Mb >= (1. - eslSMALLX1)) { /* 10 Mb  > Z >= 1 Mb */
    pli->F1 = pli->F1b = 0.35;
    pli->F2 = pli->F2b = 0.20;
    pli->F3 = pli->F3b = 0.015;
    pli->F4 = pli->F4b = 0.015;
    pli->F5 = pli->F5b = 0.015;
    pli->F6 = 0.0001;
  }
  else { /* 1 Mb  > Z */
    pli->do_msv = FALSE;
    pli->F1 = pli->F1b = 1.00; /* this is irrelevant */
    pli->F2 = pli->F2b = 0.25;
    pli->F3 = pli->F3b = 0.02;
    pli->F4 = pli->F4b = 0.02;
    pli->F5 = pli->F5b = 0.02;
    pli->F6 = 0.0001;
  }

  pli->F6env     = ESL_MIN(1.0, pli->F6 * (float) esl_opt_GetInteger(go, "--cykenvx"));
  pli->do_cykenv = (esl_opt_GetBoolean(go, "--nocykenv")) ? FALSE : TRUE;
    
  /* obey manually-defined filter thresholds */
  if(esl_opt_IsUsed(go, "--F1"))  { pli->do_msv      = TRUE; pli->F1  = esl_opt_GetReal(go, "--F1");  }
  if(esl_opt_IsUsed(go, "--F1b")) { pli->do_msvbias  = TRUE; pli->F1b = esl_opt_GetReal(go, "--F1b"); }
  if(esl_opt_IsUsed(go, "--F2"))  { pli->do_vit      = TRUE; pli->F2  = esl_opt_GetReal(go, "--F2");  }
  if(esl_opt_IsUsed(go, "--F2b")) { pli->do_vitbias  = TRUE; pli->F2b = esl_opt_GetReal(go, "--F2b"); }
  if(esl_opt_IsUsed(go, "--F3"))  { pli->do_fwd      = TRUE; pli->F3  = esl_opt_GetReal(go, "--F3");  }
  if(esl_opt_IsUsed(go, "--F3b")) { pli->do_fwdbias  = TRUE; pli->F3b = esl_opt_GetReal(go, "--F3b"); }
  if(esl_opt_IsUsed(go, "--F4"))  { pli->do_gfwd     = TRUE; pli->F4  = esl_opt_GetReal(go, "--F4");  }
  if(esl_opt_IsUsed(go, "--F4b")) { pli->do_gfwdbias = TRUE; pli->F4b = esl_opt_GetReal(go, "--F4b"); }
  if(esl_opt_IsUsed(go, "--F5"))  {                          pli->F5  = esl_opt_GetReal(go, "--F5");  }
  if(esl_opt_IsUsed(go, "--F5b")) { pli->do_edefbias = TRUE; pli->F5b = esl_opt_GetReal(go, "--F5b"); }
  if(esl_opt_IsUsed(go, "--F6"))  { pli->do_cyk      = TRUE; pli->F6  = esl_opt_GetReal(go, "--F6");  }
  
  pli->orig_F4  = pli->F4;
  pli->orig_F4b = pli->F4b;
  pli->orig_F5  = pli->F5;
  pli->orig_F5b = pli->F5b;

  if(esl_opt_GetBoolean(go, "--noF1"))     pli->do_msv        = FALSE; 
  if(esl_opt_GetBoolean(go, "--shortmsv")) pli->do_shortmsv   = TRUE;
  if(esl_opt_GetBoolean(go, "--noF2"))     pli->do_vit        = FALSE; 
  if(esl_opt_GetBoolean(go, "--noF3"))     pli->do_fwd        = FALSE; 
  if(esl_opt_GetBoolean(go, "--noF4"))     pli->do_gfwd       = FALSE; 
  if(esl_opt_GetBoolean(go, "--noF6"))     pli->do_cyk        = FALSE; 
  if(esl_opt_GetBoolean(go, "--doF1b"))    pli->do_msvbias    = TRUE;
  if(esl_opt_GetBoolean(go, "--noF2b"))    pli->do_vitbias    = FALSE;
  if(esl_opt_GetBoolean(go, "--noF3b"))    pli->do_fwdbias    = FALSE;
  if(esl_opt_GetBoolean(go, "--noF4b"))    pli->do_gfwdbias   = FALSE;
  if(esl_opt_GetBoolean(go, "--noF5b"))    pli->do_edefbias   = FALSE;
  if(esl_opt_GetBoolean(go, "--hmm")) { 
    pli->do_cm  = FALSE;
    pli->do_cyk = FALSE; 
  }
  if(esl_opt_GetBoolean(go, "--nohmm")) { 
    pli->do_hmm        = FALSE; 
    pli->do_msv        = FALSE; 
    pli->do_vit        = FALSE;
    pli->do_fwd        = FALSE; 
    pli->do_gfwd       = FALSE; 
    pli->do_envelopes  = FALSE; 
    pli->do_msvbias    = FALSE;
    pli->do_vitbias    = FALSE;
    pli->do_fwdbias    = FALSE;
    pli->do_gfwdbias   = FALSE;
    pli->do_edefbias   = FALSE;
  }
  if(esl_opt_GetBoolean(go, "--max")) { /* turn off all filters */
    pli->do_max = TRUE;
    pli->do_hmm = pli->do_msv = pli->do_vit = pli->do_fwd = pli->do_cyk = FALSE; 
    pli->do_msvbias = pli->do_vitbias = pli->do_fwdbias = pli->do_gfwdbias = pli->do_edefbias = FALSE;
    pli->do_cykenv = FALSE;
    pli->F1 = pli->F2 = pli->F3 = pli->F4 = pli->F5 = 1.0;
  }
  if(esl_opt_GetBoolean(go, "--rfam")) { /* set up strict-level filtering, these are the same as the defaults
					  * for a 100 Gb database or larger */
    pli->do_rfam = TRUE;
    pli->F1 = pli->F1b = 0.05;
    pli->F2 = pli->F2b = 0.04;
    pli->F3 = pli->F3b = 0.0004;
    pli->F4 = pli->F4b = 0.0004;
    pli->F5 = pli->F5b = 0.0004;
    pli->F6 = 0.0001;
  }

  pli->do_time_F1   = esl_opt_GetBoolean(go, "--time-F1");
  pli->do_time_F2   = esl_opt_GetBoolean(go, "--time-F2");
  pli->do_time_F3   = esl_opt_GetBoolean(go, "--time-F3");
  pli->do_time_F4   = esl_opt_GetBoolean(go, "--time-F4");
  pli->do_time_F5   = esl_opt_GetBoolean(go, "--time-F5");
  pli->do_time_F6   = esl_opt_GetBoolean(go, "--time-F6");

  if (esl_opt_GetBoolean(go, "--null2"))   pli->do_null2      = TRUE;
  if (esl_opt_GetBoolean(go, "--nonull3")) pli->do_null3      = FALSE;

  /* Configure options for the CM stages */
  pli->fcyk_cm_search_opts  = 0;
  pli->final_cm_search_opts = 0;
  pli->fcyk_beta  = esl_opt_GetReal(go, "--fbeta");
  pli->fcyk_tau   = esl_opt_GetReal(go, "--ftau");
  pli->final_beta = esl_opt_GetReal(go, "--beta");
  pli->final_tau  = esl_opt_GetReal(go, "--tau");

  if(! esl_opt_GetBoolean(go, "--hmm")) { 
    /* set up CYK filter and final round parameters */
    /* three options for banding in CYK filter and final round: 
     * --{f}nonbanded: non-banded
     * --{f}qdb:       QDB
     * neither:        HMM bands
     */
    /* CYK filter settings */
    if     (esl_opt_GetBoolean(go, "--fnonbanded"))  pli->fcyk_cm_search_opts  |= CM_SEARCH_NONBANDED;
    else if(esl_opt_GetBoolean(go, "--fqdb"))        pli->fcyk_cm_search_opts  |= CM_SEARCH_QDB;
    else                                             pli->fcyk_cm_search_opts  |= CM_SEARCH_HBANDED;

    if(  esl_opt_GetBoolean(go, "--fsums"))          pli->fcyk_cm_search_opts  |= CM_SEARCH_SUMS;
    if(! esl_opt_GetBoolean(go, "--nonull3"))        pli->fcyk_cm_search_opts  |= CM_SEARCH_NULL3;

    /* set up final round parameters */
    if(! esl_opt_GetBoolean(go, "--cyk"))            pli->final_cm_search_opts |= CM_SEARCH_INSIDE;

    if     (esl_opt_GetBoolean(go, "--nonbanded"))   pli->final_cm_search_opts |= CM_SEARCH_NONBANDED;
    else if(esl_opt_GetBoolean(go, "--qdb"))         pli->final_cm_search_opts |= CM_SEARCH_QDB;
    else                                             pli->final_cm_search_opts |= CM_SEARCH_HBANDED;

    if(  esl_opt_GetBoolean(go, "--sums"))           pli->final_cm_search_opts |= CM_SEARCH_SUMS;
    if(! esl_opt_GetBoolean(go, "--nonull3"))        pli->final_cm_search_opts |= CM_SEARCH_NULL3;
    if(esl_opt_GetBoolean(go, "--nogreedy"))         pli->final_cm_search_opts |= CM_SEARCH_CMNOTGREEDY;
  }

  /* Determine cm->config_opts we'll use to configure CMs after reading within a SCAN pipeline */
  pli->cm_config_opts = 0;
  /* should we configure CM/CP9 into local mode? */
  if(! esl_opt_GetBoolean(go, "-g")) { 
    pli->cm_config_opts |= CM_CONFIG_LOCAL;
    if(! esl_opt_GetBoolean(go, "--cp9gloc")) { 
      pli->cm_config_opts |= CM_CONFIG_HMMLOCAL;
      if(! esl_opt_GetBoolean(go, "--cp9noel")) pli->cm_config_opts |= CM_CONFIG_HMMEL; 
    }
  }
  /* should we setup for truncated alignments? */
  if(! esl_opt_GetBoolean(go, "--notrunc")) pli->cm_config_opts |= CM_CONFIG_TRUNC; 
  if(  esl_opt_GetBoolean(go, "--noforce")) pli->cm_config_opts |= CM_CONFIG_TRUNC_NOFORCE;

  /* will we be requiring a CM_SCAN_MX? a CM_TR_SCAN_MX? */
  if(esl_opt_GetBoolean(go, "--fqdb") || esl_opt_GetBoolean(go, "--qdb")) { 
    pli->cm_config_opts |= CM_CONFIG_SCANMX;
    if(! esl_opt_GetBoolean(go, "--notrunc")) { 
      pli->cm_config_opts |= CM_CONFIG_TRSCANMX; 
    }
  }

  /* Determine statistics modes for CM stages */
  pli->do_glocal_cm_stages = (esl_opt_GetBoolean(go, "-g")) ? TRUE : FALSE;
  pli->fcyk_cm_exp_mode  = pli->do_glocal_cm_stages ? EXP_CM_GC : EXP_CM_LC;
  if(pli->final_cm_search_opts & CM_SEARCH_INSIDE) { 
    pli->final_cm_exp_mode = pli->do_glocal_cm_stages ? EXP_CM_GI : EXP_CM_LI;
  }
  else {
    pli->final_cm_exp_mode = pli->do_glocal_cm_stages ? EXP_CM_GC : EXP_CM_LC;
  }

  /* Accounting as we collect results */
  pli->nmodels = 0;
  pli->nseqs   = 0;
  pli->nnodes  = 0;
  for(pass_idx = 0; pass_idx < NPLI_PASSES; pass_idx++) { 
    cm_pli_ZeroAccounting(&(pli->acct[pass_idx]));
  }

  pli->mode             = mode;
  pli->do_top           = (esl_opt_GetBoolean(go, "--bottomonly"))    ? FALSE : TRUE;
  pli->do_bot           = (esl_opt_GetBoolean(go, "--toponly"))       ? FALSE : TRUE;
  pli->show_accessions  = (esl_opt_GetBoolean(go, "--acc")            ? TRUE  : FALSE);
  pli->do_alignments    = (esl_opt_GetBoolean(go, "--noali")          ? FALSE : TRUE);
  pli->align_cyk        = (esl_opt_GetBoolean(go, "--aln-cyk")        ? TRUE  : FALSE);
  pli->align_hbanded    = (esl_opt_GetBoolean(go, "--aln-nonbanded")  ? FALSE : TRUE);
  pli->hb_size_limit    = esl_opt_GetReal    (go, "--aln-sizelimit");
  pli->do_hb_recalc     = esl_opt_GetBoolean(go, "--aln-newbands")    ? TRUE  : FALSE;

  pli->cur_cm_idx   = -1;
  pli->cur_seq_idx  = -1;
  pli->cur_pass_idx = -1;

  pli->abc              = abc;
  pli->cmfp             = NULL;
  pli->errbuf[0]        = '\0';

  return pli;

 ERROR:
  cm_pipeline_Destroy(pli, NULL);
  return NULL;
}

/* Function:  cm_pipeline_Reuse()
 * Synopsis:  Reuse a pipeline for next target.
 * Incept:    EPN, Fri Sep 24 16:28:55 2010
 *            SRE, Fri Dec  5 10:31:36 2008 [Janelia] (p7_pipeline_Reuse())
 *
 * Purpose:   Reuse <pli> for next target sequence (search mode)
 *            or model (scan mode). 
 *            
 *            May eventually need to distinguish from reusing pipeline
 *            for next query, but we're not really focused on multiquery
 *            use of hmmscan/hmmsearch/phmmer for the moment.
 */
int
cm_pipeline_Reuse(CM_PIPELINE *pli)
{
  p7_omx_Reuse(pli->oxf);
  p7_omx_Reuse(pli->oxb);
  p7_omx_Reuse(pli->fwd);
  p7_omx_Reuse(pli->bck);
  p7_domaindef_Reuse(pli->ddef);
  /* TO DO, write scanmatrixreuse */
  return eslOK;
}


/* Function:  cm_pipeline_Destroy()
 * Synopsis:  Free a <CM_PIPELINE> object.
 * Incept:    EPN, Fri Sep 24 16:29:34 2010
 *            SRE, Fri Dec  5 10:30:23 2008 [Janelia] (p7_pipeline_Destroy())
 *
 * Purpose:   Free a <CM_PIPELINE> object. Requires a cm (sigh) to free the scan matrix.
 */
void
cm_pipeline_Destroy(CM_PIPELINE *pli, CM_t *cm)
{
  if (pli == NULL) return;
  
  p7_omx_Destroy(pli->oxf);
  p7_omx_Destroy(pli->oxb);
  p7_omx_Destroy(pli->fwd);
  p7_omx_Destroy(pli->bck);
  p7_gmx_Destroy(pli->gfwd);
  p7_gmx_Destroy(pli->gbck);
  p7_gmx_Destroy(pli->gxf);
  p7_gmx_Destroy(pli->gxb);
  esl_randomness_Destroy(pli->r);
  p7_domaindef_Destroy(pli->ddef);
  if(pli->tro != NULL) free(pli->tro);
  free(pli);
}
/*---------------- end, CM_PIPELINE object ----------------------*/


/*****************************************************************
 * 2. The pipeline API.
 *****************************************************************/

/* Function:  cm_pli_TargetReportable
 * Synopsis:  Returns TRUE if target score meets reporting threshold.
 * Incept:    EPN, Fri Sep 24 16:34:37 2010
 *            SRE, Tue Dec  9 08:57:26 2008 [Janelia] (p7_pli_TargetReportable())
 *
 * Purpose:   Returns <TRUE> if the bit score <score> and/or 
 *            E-value <Eval> meets per-target reporting thresholds 
 *            for the processing pipeline.
 */
int
cm_pli_TargetReportable(CM_PIPELINE *pli, float score, double Eval)
{
  if      (  pli->by_E && Eval  <= pli->E) return TRUE;
  else if (! pli->by_E && score >= pli->T) return TRUE;
  return FALSE;
}

/* Function:  cm_pli_TargetIncludable()
 * Synopsis:  Returns TRUE if target score meets inclusion threshold.
 * Incept:    EPN, Fri Sep 24 16:32:43 2010
 *            SRE, Fri Jan 16 11:18:08 2009 [Janelia] (p7_pli_TargetIncludable())
 */
int
cm_pli_TargetIncludable(CM_PIPELINE *pli, float score, double Eval)
{
  if      (  pli->by_E && Eval  <= pli->incE) return TRUE;
  else if (! pli->by_E && score >= pli->incT) return TRUE;

  return FALSE;
}

/* Function:  cm_pli_NewModel()
 * Synopsis:  Prepare pipeline for a new CM/HMM
 * Incept:    EPN, Fri Sep 24 16:35:35 2010
 *            SRE, Fri Dec  5 10:35:37 2008 [Janelia] (p7_pli_NewModel())
 *
 * Purpose:   Caller has a new model.
 *            Prepare the pipeline <pli> to receive this model as either 
 *            a query or a target.
 *
 *            The information of the model we may receive varies, as
 *            indicated by <modmode> and <pli->mode>. This is enforced 
 *            by a contract check upon entrance, failure causes immediate
 *            return of eslEINCOMPAT.
 *
 *      case  <pli->mode>     <modmode>         <cm>      <om> and <bg>
 *      ----  --------------  ---------------   --------  -------------
 *         1  CM_SEARCH_SEQS  CM_NEWMODEL_CM    non-null  non-null     
 *         2  CM_SCAN_SEQS    CM_NEWMODEL_MSV   NULL      non-null     
 *         3  CM_SCAN_SEQS    CM_NEWMODEL_CM    non-null  NULL         
 *
 *            <cm_clen> and <cm_W> are always valid, but are only 
 *            necessary for case 2.
 *
 *            Note that if we're in SEARCH mode, <modmode> will 
 *            always be CM_NEWMODEL_CM. In SCAN mode, we may only
 *            call this function once with <modmode> == CM_NEWMODEL_MSV
 *            (case 2 above, this happens if no hit from the query
 *            sequence survives to the CM stage). Also, if we
 *            are in SCAN mode with modmode==CM_NEWMODEL_CM we must
 *            have entered this function previously for the same
 *            'model' with modmode==CM_NEWMODEL_MSV.
 *
 *            The pipeline may alter the null model in <bg> in a
 *            model-specific way (if we're using composition bias
 *            filter HMMs in the pipeline).
 *
 * Returns:   <eslOK> on success.
 *
 *            <eslEINCOMPAT> in contract is not met.
 * 
 *            <eslEINVAL> if pipeline expects to be able to use a
 *            model's bit score thresholds, but this model does not
 *            have the appropriate ones set.
 */
int
cm_pli_NewModel(CM_PIPELINE *pli, int modmode, CM_t *cm, int cm_clen, int cm_W, P7_OPROFILE *om, P7_BG *bg, int64_t cur_cm_idx)
{
  int status = eslOK;
  int i, nsteps;
  float T;

  /* check contract */
  if(pli->mode == CM_SEARCH_SEQS) { /* case 1 */
    if(modmode != CM_NEWMODEL_CM) 
    if(cm == NULL)                 ESL_FAIL(eslEINCOMPAT, pli->errbuf, "cm_pli_NewModel(), contract violated, SEARCH mode and CM is NULL"); 
    if(cm->clen != cm_clen)        ESL_FAIL(eslEINCOMPAT, pli->errbuf, "cm_pli_NewModel(), contract violated, cm->clen != cm_clen"); 
    if(cm->W    != cm_W)           ESL_FAIL(eslEINCOMPAT, pli->errbuf, "cm_pli_NewModel(), contract violated, cm->W != cm_W"); 
    if(om == NULL)                 ESL_FAIL(eslEINCOMPAT, pli->errbuf, "cm_pli_NewModel(), contract violated, SEARCH mode and om is NULL"); 
    if(bg == NULL)                 ESL_FAIL(eslEINCOMPAT, pli->errbuf, "cm_pli_NewModel(), contract violated, SEARCH mode and bg is NULL"); 
  }
  else if(pli->mode == CM_SCAN_MODELS) { 
    if(modmode == CM_NEWMODEL_MSV) { /* case 2 */
      if(cm != NULL)              ESL_FAIL(eslEINCOMPAT, pli->errbuf, "cm_pli_NewModel(), contract violated, SCAN/MSV mode, and CM is non-NULL"); 
      if(om == NULL)              ESL_FAIL(eslEINCOMPAT, pli->errbuf, "cm_pli_NewModel(), contract violated, SCAN/MSV mode, and om is NULL"); 
      if(bg == NULL)              ESL_FAIL(eslEINCOMPAT, pli->errbuf, "cm_pli_NewModel(), contract violated, SCAN/MSV mode, and bg is NULL"); 
    }
    else if(modmode == CM_NEWMODEL_CM) { /* case 3 */
      if(cm == NULL)              ESL_FAIL(eslEINCOMPAT, pli->errbuf, "cm_pli_NewModel(), contract violated, SCAN/CM mode, and CM is NULL"); 
      if(cm->clen != cm_clen)     ESL_FAIL(eslEINCOMPAT, pli->errbuf, "cm_pli_NewModel(), contract violated, cm->clen != cm_clen");
      if(cm->W    != cm_W)        ESL_FAIL(eslEINCOMPAT, pli->errbuf, "cm_pli_NewModel(), contract violated, cm->W != cm_W");
      if(om != NULL)              ESL_FAIL(eslEINCOMPAT, pli->errbuf, "cm_pli_NewModel(), contract violated, SCAN/CM mode, and om is non-NULL"); 
      if(bg != NULL)              ESL_FAIL(eslEINCOMPAT, pli->errbuf, "cm_pli_NewModel(), contract violated, SCAN/CM mode, and bg is non-NULL"); 
    }
  }

  pli->cur_cm_idx = cur_cm_idx;

  /* Two sets (A and B) of value updates: 
   * case 1: we do both sets 
   * case 2: we do set A only
   * case 3: we do set B only
   */
  if(pli->mode == CM_SEARCH_SEQS || modmode == CM_NEWMODEL_MSV) { 
    /* set A updates: case 1 and 2 do these */
    pli->nmodels++;
    pli->nnodes += cm_clen;

    if (pli->do_msvbias || pli->do_vitbias || pli->do_fwdbias || pli->do_gfwdbias || pli->do_edefbias) { 
      p7_bg_SetFilter(bg, om->M, om->compo);
    }
    /* copy some values from the model */
    pli->cmW  = cm_W;
    pli->clen = cm_clen;

    /* determine pli->maxW, this will be the number of residues
     * that must overlap between adjacent windows on a
     * single sequence, this is MAX of cm->W and om->max_length
     */
    pli->maxW = ESL_MAX(cm_W, om->max_length);

    /* reset consensus length specific thresholds  */
    pli->F4  = pli->orig_F4;
    pli->F4b = pli->orig_F4b;
    pli->F5  = pli->orig_F5;
    pli->F5b = pli->orig_F5b;
    
    /* update consensus length specific thresholds, if nec */
    if(pli->do_glen) { 
      if(pli->clen >= pli->glen_min) { 
	nsteps = (int) (1 + ((ESL_MAX(pli->clen, pli->glen_max) - (pli->glen_min-1)) / pli->glen_step));
	for(i = 0; i < nsteps; i++) { 
	  pli->F4  /= 2.;
	  pli->F4b /= 2.;
	  pli->F5  /= 2.;
	  pli->F5b /= 2.;
	}
      }
    }
  }
  if(pli->mode == CM_SEARCH_SEQS || modmode == CM_NEWMODEL_CM) { 
    /* set B updates: case 1 and 3 do these (they require a valid CM) */

    /* Update the current effective database size so it pertains to the
     * new model. Also, If we're using an E-value threshold determine
     * the bit score for this model that pertains to that E-value.
     */
    if((status = UpdateExpsForDBSize(cm, pli->errbuf, (long) pli->Z)) != eslOK) return status;
    if(pli->by_E) { 
      if((status = E2ScoreGivenExpInfo(cm->expA[pli->final_cm_exp_mode], pli->errbuf, pli->E, &T)) != eslOK) ESL_FAIL(status, pli->errbuf, "problem determining min score for E-value %6g for model %s\n", pli->E, cm->name);
      pli->T = (double) T;
    }

    /* if we're using Rfam GA, NC, or TC cutoffs, update them for this model */
    if (pli->use_bit_cutoffs) { 
      if((status = cm_pli_NewModelThresholds(pli, cm)) != eslOK) return status;
    }
  }
  return eslOK;
}

/* Function:  cm_pli_NewModelThresholds()
 * Synopsis:  Set reporting and inclusion bit score thresholds on a new model.
 * Incept:    EPN, Wed Jun 15 14:40:46 2011
 *            SRE, Sat Oct 17 12:07:43 2009 [Janelia] (p7_pli_NewModelThresholds()
 *
 * Purpose:   Set the bit score thresholds on a new model, if we're 
 *            either using Rfam GA, TC, or NC cutoffs for reporting or
 *            inclusion, and/or if we already know the total 
 *            database size.
 *            
 *            In a "search" pipeline, this only needs to be done once
 *            per query model, so <cm_pli_NewModelThresholds()> gets 
 *            called by <cm_pli_NewModel()>.
 *            
 *            In a "scan" pipeline, this needs to be called for each
 *            target model.
 *
 * Returns:   <eslOK> on success. 
 *            
 *            <eslEINVAL> if pipeline expects to be able to use a
 *            model's bit score thresholds, but this model does not
 *            have the appropriate ones set.
 *
 * Xref:      Written to fix bug #h60 (p7_pli_NewModelThreshold())
 */
int
cm_pli_NewModelThresholds(CM_PIPELINE *pli, CM_t *cm)
{

  if (pli->use_bit_cutoffs) { 
    if(pli->use_bit_cutoffs == CMH_GA) {
      if (! (cm->flags & CMH_GA)) ESL_FAIL(eslEINVAL, pli->errbuf, "GA bit threshold unavailable for model %s\n", cm->name);
      pli->T = pli->incT = cm->ga;
    }
    else if(pli->use_bit_cutoffs == CMH_TC) {
      if (! (cm->flags & CMH_TC)) ESL_FAIL(eslEINVAL, pli->errbuf, "TC bit threshold unavailable for model %s\n", cm->name);
      pli->T = pli->incT = cm->tc;
    }
    else if(pli->use_bit_cutoffs == CMH_NC) {
      if (! (cm->flags & CMH_NC)) ESL_FAIL(eslEINVAL, pli->errbuf, "NC bit threshold unavailable for model %s\n", cm->name);
      pli->T = pli->incT = cm->nc;
    }
  }
  return eslOK;
}

/* Function:  cm_pli_NewSeq()
 * Synopsis:  Prepare pipeline for a new sequence (target or query)
 * Incept:    EPN, Fri Sep 24 16:39:52 2010
 *            SRE, Fri Dec  5 10:57:15 2008 [Janelia] (p7_pli_NewSeq())
 *
 * Purpose:   Caller has a new sequence <sq>. Prepare the pipeline <pli>
 *            to receive this model as either a query or a target.
 *
 * Returns:   <eslOK> on success.
 */
 int
 cm_pli_NewSeq(CM_PIPELINE *pli, const ESL_SQ *sq, int64_t cur_seq_idx)
{
  /* Update number of residues read/searched in standard pipeline passes */
  pli->acct[PLI_PASS_STD].nres += sq->n;

  /* Update cur_seq_idx, which is a unique identifier for the sequence, so 
   * we can reliably remove overlaps. This index is copied to all the 
   * hit objects for all hits found by the pipeline when searching this sequence. */
  pli->cur_seq_idx = cur_seq_idx;

  /* Note that we don't update pli->Z, that must be set at beginning
   * of a search, this differs from hmmsearch and nhmmer, which by
   * default update Z as sequences are read.
   */
  return eslOK;
}

/* Function:  cm_pipeline_Merge()
 * Synopsis:  Merge the pipeline statistics
 * Incept:    EPN, Fri Sep 24 16:41:17 2010   
 *
 * Purpose:   Caller has a new model <om>. Prepare the pipeline <pli>
 *            to receive this model as either a query or a target.
 *
 *            The pipeline may alter the null model <bg> in a model-specific
 *            way (if we're using a composition bias filter HMM in the
 *            pipeline).
 *
 * Returns:   <eslOK> on success.
 * 
 *            <eslEINVAL> if pipeline expects to be able to use a
 *            model's bit score thresholds, but this model does not
 *            have the appropriate ones set.
 */
int
cm_pipeline_Merge(CM_PIPELINE *p1, CM_PIPELINE *p2)
{
  /* if we are searching a sequence database, we need to keep track of the
   * number of sequences and residues processed.
   */
  int p; /* counter over pipeline passes */

  if (p1->mode == CM_SEARCH_SEQS)
    {
      p1->nseqs   += p2->nseqs;
      for(p = 0; p < NPLI_PASSES; p++) p1->acct[p].nres += p2->acct[p].nres;
    }
  else
    {
      p1->nmodels += p2->nmodels;
      p1->nnodes  += p2->nnodes;
    }

  for(p = 0; p < NPLI_PASSES; p++) { 
    p1->acct[p].n_past_msv  += p2->acct[p].n_past_msv;
    p1->acct[p].n_past_vit  += p2->acct[p].n_past_vit;
    p1->acct[p].n_past_fwd  += p2->acct[p].n_past_fwd;
    p1->acct[p].n_past_gfwd += p2->acct[p].n_past_gfwd;
    p1->acct[p].n_past_edef += p2->acct[p].n_past_edef;
    p1->acct[p].n_past_cyk  += p2->acct[p].n_past_cyk;
    p1->acct[p].n_past_ins  += p2->acct[p].n_past_ins;
    p1->acct[p].n_output    += p2->acct[p].n_output;

    p1->acct[p].n_past_msvbias += p2->acct[p].n_past_msvbias;
    p1->acct[p].n_past_vitbias += p2->acct[p].n_past_vitbias;
    p1->acct[p].n_past_fwdbias += p2->acct[p].n_past_fwdbias;
    p1->acct[p].n_past_gfwdbias+= p2->acct[p].n_past_gfwdbias;
    p1->acct[p].n_past_edefbias += p2->acct[p].n_past_edefbias;

    p1->acct[p].pos_past_msv  += p2->acct[p].pos_past_msv;
    p1->acct[p].pos_past_vit  += p2->acct[p].pos_past_vit;
    p1->acct[p].pos_past_fwd  += p2->acct[p].pos_past_fwd;
    p1->acct[p].pos_past_gfwd += p2->acct[p].pos_past_gfwd;
    p1->acct[p].pos_past_edef += p2->acct[p].pos_past_edef;
    p1->acct[p].pos_past_cyk  += p2->acct[p].pos_past_cyk;
    p1->acct[p].pos_past_ins  += p2->acct[p].pos_past_ins;
    p1->acct[p].pos_output    += p2->acct[p].pos_output;

    p1->acct[p].pos_past_msvbias += p2->acct[p].pos_past_msvbias;
    p1->acct[p].pos_past_vitbias += p2->acct[p].pos_past_vitbias;
    p1->acct[p].pos_past_fwdbias += p2->acct[p].pos_past_fwdbias;
    p1->acct[p].pos_past_gfwdbias+= p2->acct[p].pos_past_gfwdbias;
    p1->acct[p].pos_past_edefbias += p2->acct[p].pos_past_edefbias;

    p1->acct[p].n_overflow_fcyk  += p2->acct[p].n_overflow_fcyk;
    p1->acct[p].n_overflow_final += p2->acct[p].n_overflow_final;
    p1->acct[p].n_aln_hb         += p2->acct[p].n_aln_hb;
    p1->acct[p].n_aln_dccyk      += p2->acct[p].n_aln_dccyk;
  }
  if (p1->Z_setby == p7_ZSETBY_NTARGETS) { 
    p1->Z += (p1->mode == p7_SCAN_MODELS) ? p2->nmodels : p2->nseqs;
  }

  return eslOK;
}

/* Function:  cm_pli_p7Filter()
 * Synopsis:  The accelerated p7 comparison pipeline: MSV through Forward filter.
 * Incept:    EPN, Wed Nov 24 13:07:02 2010
 *            TJW, Fri Feb 26 10:17:53 2018 [Janelia] (p7_Pipeline_Longtargets())
 *
 * Purpose:   Run the accelerated pipeline to compare profile <om>
 *            against sequence <sq>. Some combination of the MSV,
 *            Viterbi and Forward algorithms are used, based on 
*            option flags set in <pli>. 
 *
 *            In a normal pipeline run, this function call should be
 *            followed by a call to cm_pli_p7BoundaryDef().
 *
 * Returns:   <eslOK> on success. For the <ret_nwin> windows that
 *            survive all filters, the start and end positions of the
 *            windows are stored and * returned in <ret_ws> and
 *            <ret_we> respectively.

 *            <eslEINVAL> if (in a scan pipeline) we're supposed to
 *            set GA/TC/NC bit score thresholds but the model doesn't
 *            have any.
 *
 *            <eslERANGE> on numerical overflow errors in the
 *            optimized vector implementations; particularly in
 *            posterior decoding. I don't believe this is possible for
 *            multihit local models, but I'm set up to catch it
 *            anyway. We may emit a warning to the user, but cleanly
 *            skip the problematic sequence and continue.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 * Xref:      J4/25.
 */
int
cm_pli_p7Filter(CM_PIPELINE *pli, P7_OPROFILE *om, P7_BG *bg, float *p7_evparam, FM_HMMDATA *fm_hmmdata, const ESL_SQ *sq, int64_t **ret_ws, int64_t **ret_we, int *ret_nwin)
{
  int              status;
  float            mfsc, vfsc, fwdsc;/* filter scores                           */
  float            filtersc;         /* HMM null filter score                   */
  int              have_filtersc;    /* TRUE if filtersc has been calc'ed for current window */
  float            nullsc;           /* null model score */
  float            wsc;              /* the corrected bit score for a window */
  double           P;                /* P-value of a hit */
  int              i, i2;            /* counters */
  int              wlen;             /* length of current window */
  void            *p;                /* for ESL_RALLOC */
  int             *useme = NULL;     /* used when merging overlapping windows */
  int              overlap;          /* number of overlapping positions b/t 2 adjacent windows */
  float            wcorr;            /* bit sc correction for a window */
  int            **survAA = NULL;    /* [0..s..Np7_SURV-1][0..i..nwin-1] TRUE if window i survived stage s */
  int              nalloc;           /* currently allocated size for ws, we */
  int64_t         *ws = NULL;        /* [0..nwin-1] window start positions */
  int64_t         *we = NULL;        /* [0..nwin-1] window end   positions */
  double          *wp = NULL;        /* [0..nwin-1] window P-values, P-value of furthest-reached filter algorithm */
  int              nwin;             /* number of windows */
  int64_t         *new_ws = NULL;    /* used when copying/modifying ws */
  int64_t         *new_we = NULL;    /* used when copying/modifying we */
  int              nsurv_fwd;        /* number of windows that survive fwd filter */
  int              new_nsurv_fwd;    /* used when merging fwd survivors */
  ESL_DSQ         *subdsq;           /* a ptr to the first position of a window */
  int              have_rest;        /* do we have the full <om> read in? */

  if (sq->n == 0) return eslOK;    /* silently skip length 0 seqs; they'd cause us all sorts of weird problems */
  p7_omx_GrowTo(pli->oxf, om->M, 0, sq->n);    /* expand the one-row omx if needed */
  have_rest = (om->mode == p7_NO_MODE) ? FALSE : TRUE; /* we use om->mode as a flag to tell us whether we already have read the full <om> from disk or not */

  /* Set false target length. This is a conservative estimate of the length of window that'll
   * soon be passed on to later phases of the pipeline;  used to recover some bits of the score
   * that we would miss if we left length parameters set to the full target length */

  if(pli->do_filcmW) { 
    p7_oprofile_ReconfigMSVLength(om, pli->cmW);
    om->max_length = pli->cmW;
  }
  else { 
    p7_oprofile_ReconfigMSVLength(om, om->max_length); /* nhmmer's way */
  }

  #if DOPRINT
  printf("\nPIPELINE p7Filter() %s  %" PRId64 " residues\n", sq->name, sq->n);
  #endif

  /* initializations */
  nsurv_fwd = 0;
  nwin = 0;

  /***************************************************/
  /* Filter 1: MSV, long target-variant, with p7 HMM */
  if(pli->do_msv && (! pli->do_shortmsv)) { 
    /* ws_int, we_int: workaround for p7_MSVFilter_longtarget() taking integers, instead of int64_t */
    int *ws_int;
    int *we_int;
    p7_MSVFilter_longtarget(sq->dsq, sq->n, om, pli->oxf, bg, pli->F1, &ws_int, &we_int, &nwin);
    ESL_ALLOC(ws, sizeof(int64_t) * nwin);
    ESL_ALLOC(we, sizeof(int64_t) * nwin);
    for(i = 0; i < nwin; i++) { ws[i] = ws_int[i]; we[i] = we_int[i]; }
    //free(ws_int); free(we_int);
#if 0
    /* UNCOMMENT THIS BLOCK WHEN UPGRADING TO NEW PROTOTYPE FOR MSVFILTER */
    FM_WINDOWLIST    windowlist;       /* list of windows, structure taken by p7_MSVFilter_longtarget() */
    fm_initWindows(&windowlist);
    status = p7_MSVFilter_longtarget(sq->dsq, sq->n, om, pli->oxf, fm_hmmdata, bg, pli->F1, &windowlist, FALSE);
    ESL_ALLOC(ws, sizeof(int64_t) * windowlist.count);
    ESL_ALLOC(we, sizeof(int64_t) * windowlist.count);
    nwin = windowlist.count;
    for(i = 0; i < nwin; i++) { 
      ws[i] =         windowlist.windows[i].n;
      we[i] = ws[i] + windowlist.windows[i].length - 1;
    }
    free(windowlist.windows);
#endif

    if(pli->do_wsplit) { 
      /* split up windows > (pli->wmult * W) into length 2W-1, with W-1 overlapping residues */
      nalloc = nwin + 100;
      ESL_ALLOC(new_ws, sizeof(int64_t) * nalloc);
      ESL_ALLOC(new_we, sizeof(int64_t) * nalloc);
      for (i = 0, i2 = 0; i < nwin; i++, i2++) {
	wlen = we[i] - ws[i] + 1;
	if((i2+1) == nalloc) { 
	  nalloc += 100;
	  ESL_RALLOC(new_ws, p, sizeof(int64_t) * nalloc);
	  ESL_RALLOC(new_we, p, sizeof(int64_t) * nalloc);
	}
	if(wlen > (pli->wmult * pli->cmW)) { 
	  /* split this window */
	  new_ws[i2] = ws[i]; 
	  new_we[i2] = ESL_MIN((new_ws[i2] + (2 * pli->cmW) - 1), we[i]);
	  while(new_we[i2] < we[i]) { 
	    i2++;
	    if((i2+1) == nalloc) { 
	      nalloc += 100;
	      ESL_RALLOC(new_ws, p, sizeof(int64_t) * nalloc);
	      ESL_RALLOC(new_we, p, sizeof(int64_t) * nalloc);
	    }
	    new_ws[i2] = ESL_MIN(new_ws[i2-1] + pli->cmW, we[i]);
	    new_we[i2] = ESL_MIN(new_we[i2-1] + pli->cmW, we[i]);
	  }	    
	}
	else { /* do not split this window */
	  new_ws[i2] = ws[i]; 
	  new_we[i2] = we[i];
	}
      }
      free(ws);
      free(we);
      ws = new_ws;
      we = new_we;
      nwin = i2;
    }
  }
  else { /* do_msv is FALSE or do_shortmsv is TRUE, divide up into windows of 2*W, overlapping by W-1 residues */
    nwin = 1; /* first window */
    if(sq->n > (2 * pli->cmW)) { 
      nwin += (int) (sq->n - (2 * pli->cmW)) / ((2 * pli->cmW) - (pli->cmW - 1));
      /*            (L     -  first window)/(number of unique residues per window) */
      if(((sq->n - (2 * pli->cmW)) % ((2 * pli->cmW) - (pli->cmW - 1))) > 0) { 
	nwin++; /* if the (int) cast in previous line removed any fraction of a window, we add it back here */
      }
    }
    ESL_ALLOC(ws, sizeof(int64_t) * nwin);
    ESL_ALLOC(we, sizeof(int64_t) * nwin);
    for(i = 0; i < nwin; i++) { 
      ws[i] = 1 + (i * (pli->cmW + 1));
      we[i] = ESL_MIN((ws[i] + (2*pli->cmW) - 1), sq->n);
      /*printf("window %5d/%5d  %10d..%10d (L=%10" PRId64 ")\n", i+1, nwin, ws[i], we[i], sq->n);*/
    }
  }      
  if(! pli->do_shortmsv) pli->acct[pli->cur_pass_idx].n_past_msv += nwin;

  ESL_ALLOC(wp, sizeof(double) * nwin);
  for(i = 0; i < nwin; i++) { wp[i] = pli->F1; } /* TEMP (?) p7_MSVFilter_longtarget() does not return P-values */

  /*********************************************/
  /* allocate and initialize survAA, which will keep track of number of windows surviving each stage */
  ESL_ALLOC(survAA, sizeof(int *) * Np7_SURV);
  for (i = 0; i < Np7_SURV; i++) { 
    ESL_ALLOC(survAA[i], sizeof(int) * nwin);
    esl_vec_ISet(survAA[i], nwin, FALSE);
  }
    
  for (i = 0; i < nwin; i++) {
    subdsq = sq->dsq + ws[i] - 1;
    have_filtersc = FALSE;
    wlen = we[i] - ws[i] + 1;

    if(pli->do_wcorr) { 
      wcorr = (2 * log(2. / (pli->cmW+2))) - (2. * log(2. / (wlen+2)));
      p7_bg_SetLength(bg, pli->cmW);
      p7_bg_NullOne  (bg, subdsq, pli->cmW, &nullsc);
    }
    else { 
      wcorr = 0.;
      p7_bg_SetLength(bg, wlen);
      p7_bg_NullOne  (bg, subdsq, wlen, &nullsc);
    }

    if(pli->do_shortmsv) { 
      p7_oprofile_ReconfigMSVLength(om, wlen);
      p7_MSVFilter(subdsq, wlen, om, pli->oxf, &mfsc);
      wsc = (mfsc + wcorr - nullsc) / eslCONST_LOG2;
      P = esl_gumbel_surv(wsc,  p7_evparam[CM_p7_LMMU],  p7_evparam[CM_p7_LMLAMBDA]);
      wp[i] = P;
      if (P > pli->F1) continue;
      pli->acct[pli->cur_pass_idx].n_past_msv++;
    }

#if DOPRINT
    printf("SURVIVOR window %5d [%10" PRId64 "..%10" PRId64 "] survived MSV       %6.2f bits  P %g\n", i, ws[i], we[i], (pli->do_shortmsv) ? wsc : 0., wp[i]);
#endif
    survAA[p7_SURV_F1][i] = TRUE;
    
    if (pli->do_msv && pli->do_msvbias) {
      /******************************************************************************/
      /* Filter 1B: Bias filter with p7 HMM */
      /* Have to run msv again, to get the full score for the window.
	 (using the standard "per-sequence" msv filter this time). */
      p7_oprofile_ReconfigMSVLength(om, wlen);
      p7_MSVFilter(subdsq, wlen, om, pli->oxf, &mfsc);
      if(pli->do_wcorr) p7_bg_FilterScore(bg, subdsq, pli->cmW,     &filtersc);
      else              p7_bg_FilterScore(bg, subdsq, wlen, &filtersc);
      have_filtersc = TRUE;
      
      wsc = (mfsc + wcorr - filtersc) / eslCONST_LOG2;
      P   = esl_gumbel_surv(wsc,  p7_evparam[CM_p7_LMMU],  p7_evparam[CM_p7_LMLAMBDA]);
      wp[i] = P;

      if (P > pli->F1b) continue;

      /******************************************************************************/
    }
    pli->acct[pli->cur_pass_idx].n_past_msvbias++;
    survAA[p7_SURV_F1b][i] = TRUE;

#if DOPRINT
    printf("SURVIVOR window %5d [%10" PRId64 "..%10" PRId64 "] survived MSV-Bias  %6.2f bits  P %g\n", i, ws[i], we[i], (pli->do_shortmsv) ? wsc : 0., wp[i]);
#endif      
    if(pli->do_time_F1) return eslOK;
    
    /* In scan mode, if it passes the MSV filter (or if MSV filter is off) read the rest of the profile */
    if (pli->mode == CM_SCAN_MODELS && (! have_rest)) {
      if (pli->cmfp) p7_oprofile_ReadRest(pli->cmfp->hfp, om);
      /* Note: we don't call cm_pli_NewModelThresholds() yet (as p7_pipeline() 
       * does at this point), because we don't yet have the CM */
      have_rest = TRUE;
    }
    if(pli->do_msv && pli->do_msvbias) { /* we already called p7_oprofile_ReconfigMSVLength() above */
      p7_oprofile_ReconfigRestLength(om, wlen);
    }
    else { /* we did not call p7_oprofile_ReconfigMSVLength() above */
      p7_oprofile_ReconfigLength(om, wlen);
    }

    if (pli->do_vit) { 
      /******************************************************************************/
      /* Filter 2: Viterbi with p7 HMM */
      /* Second level filter: ViterbiFilter(), multihit with <om> */
      p7_ViterbiFilter(subdsq, wlen, om, pli->oxf, &vfsc);
      wsc = (vfsc + wcorr - nullsc) / eslCONST_LOG2; 
      P   = esl_gumbel_surv(wsc,  p7_evparam[CM_p7_LVMU],  p7_evparam[CM_p7_LVLAMBDA]);
      wp[i] = P;
      if (P > pli->F2) continue;
    }
    pli->acct[pli->cur_pass_idx].n_past_vit++;
    survAA[p7_SURV_F2][i] = TRUE;

#if DOPRINT
    printf("SURVIVOR window %5d [%10" PRId64 "..%10" PRId64 "] survived Vit       %6.2f bits  P %g\n", i, ws[i], we[i], wsc, wp[i]);
#endif

    /********************************************/
    if (pli->do_vit && pli->do_vitbias) { 
      if(! have_filtersc) { 
	if(pli->do_wcorr) p7_bg_FilterScore(bg, subdsq, pli->cmW,     &filtersc);
	else              p7_bg_FilterScore(bg, subdsq, wlen, &filtersc);
      }
      have_filtersc = TRUE;
      wsc = (vfsc + wcorr - filtersc) / eslCONST_LOG2;
      P = esl_gumbel_surv(wsc,  p7_evparam[CM_p7_LVMU],  p7_evparam[CM_p7_LVLAMBDA]);
      wp[i] = P;
      if (P > pli->F2b) continue;
      /******************************************************************************/
    }
    pli->acct[pli->cur_pass_idx].n_past_vitbias++;
    survAA[p7_SURV_F2b][i] = TRUE;

#if DOPRINT
    printf("SURVIVOR window %5d [%10" PRId64 "..%10" PRId64 "] survived Vit-Bias  %6.2f bits  P %g\n", i, ws[i], we[i], wsc, wp[i]);
#endif
    if(pli->do_time_F2) continue; 
    /********************************************/

    if(pli->do_fwd) { 
      /******************************************************************************/
      /* Filter 3: Forward with p7 HMM */
      /* Parse it with Forward and obtain its real Forward score. */
      p7_ForwardParser(subdsq, wlen, om, pli->oxf, &fwdsc);
      wsc = (fwdsc + wcorr - nullsc) / eslCONST_LOG2; 
      P = esl_exp_surv(wsc,  p7_evparam[CM_p7_LFTAU],  p7_evparam[CM_p7_LFLAMBDA]);
      wp[i] = P;
      if (P > pli->F3) continue;
    }
    /******************************************************************************/
    pli->acct[pli->cur_pass_idx].n_past_fwd++;
    survAA[p7_SURV_F3][i] = TRUE;

#if DOPRINT
    printf("SURVIVOR window %5d [%10" PRId64 "..%10" PRId64 "] survived Fwd       %6.2f bits  P %g\n", i, ws[i], we[i], wsc, wp[i]);
#endif

    if (pli->do_fwd && pli->do_fwdbias) { 
      if (! have_filtersc) { 
	if(pli->do_wcorr) p7_bg_FilterScore(bg, subdsq, pli->cmW, &filtersc);
	else              p7_bg_FilterScore(bg, subdsq, wlen,     &filtersc);
      }
      have_filtersc = TRUE;
      wsc = (fwdsc + wcorr - filtersc) / eslCONST_LOG2;
      P = esl_exp_surv(wsc,  p7_evparam[CM_p7_LFTAU],  p7_evparam[CM_p7_LFLAMBDA]);
      wp[i] = P;
      if (P > pli->F3b) continue;
      /******************************************************************************/
    }
    pli->acct[pli->cur_pass_idx].n_past_fwdbias++;
    nsurv_fwd++;
    survAA[p7_SURV_F3b][i] = TRUE;

#if DOPRINT
    printf("SURVIVOR window %5d [%10" PRId64 "..%10" PRId64 "] survived Fwd-Bias  %6.2f bits  P %g\n", i, ws[i], we[i], wsc, wp[i]);
#endif
    if(pli->do_time_F3) continue;
  }

  /* Go back through all windows, and tally up total number of
   * residues that survived each stage, without double-counting
   * overlapping residues.  Note, based on the way the windows were
   * split, we know that any overlapping residues must occur in
   * adjacent windows and we exploit that here.
   */
  for(i = 0; i < nwin; i++) {
    wlen = we[i] - ws[i] + 1;
    
    if(survAA[p7_SURV_F1][i])  pli->acct[pli->cur_pass_idx].pos_past_msv     += wlen; 
    if(survAA[p7_SURV_F1b][i]) pli->acct[pli->cur_pass_idx].pos_past_msvbias += wlen; 
    if(survAA[p7_SURV_F2][i])  pli->acct[pli->cur_pass_idx].pos_past_vit     += wlen; 
    if(survAA[p7_SURV_F2b][i]) pli->acct[pli->cur_pass_idx].pos_past_vitbias += wlen; 
    if(survAA[p7_SURV_F3][i])  pli->acct[pli->cur_pass_idx].pos_past_fwd     += wlen; 
    if(survAA[p7_SURV_F3b][i]) pli->acct[pli->cur_pass_idx].pos_past_fwdbias += wlen; 

    /* now subtract residues we've double counted */
    if(i > 0) { 
      overlap = we[i-1] - ws[i] + 1;
      if(overlap > 0) { 
	if(survAA[p7_SURV_F1][i]  && survAA[p7_SURV_F1][i-1])  pli->acct[pli->cur_pass_idx].pos_past_msv     -= overlap;
	if(survAA[p7_SURV_F1b][i] && survAA[p7_SURV_F1b][i-1]) pli->acct[pli->cur_pass_idx].pos_past_msvbias -= overlap;
	if(survAA[p7_SURV_F2][i]  && survAA[p7_SURV_F2][i-1])  pli->acct[pli->cur_pass_idx].pos_past_vit     -= overlap;
	if(survAA[p7_SURV_F2b][i] && survAA[p7_SURV_F2b][i-1]) pli->acct[pli->cur_pass_idx].pos_past_vitbias -= overlap;
	if(survAA[p7_SURV_F3][i]  && survAA[p7_SURV_F3][i-1])  pli->acct[pli->cur_pass_idx].pos_past_fwd     -= overlap;
	if(survAA[p7_SURV_F3b][i] && survAA[p7_SURV_F3b][i-1]) pli->acct[pli->cur_pass_idx].pos_past_fwdbias -= overlap;
      }
    }
  }
 
  /* Finally, create list of just those that survived fwd, and merge any overlapping windows together */
  if(nsurv_fwd > 0) { 
    ESL_ALLOC(new_ws, sizeof(int64_t) * nsurv_fwd);
    ESL_ALLOC(new_we, sizeof(int64_t) * nsurv_fwd);
    for (i = 0, i2 = 0; i < nwin; i++) { 
      if(survAA[p7_SURV_F3b][i]) { 
	new_ws[i2] = ws[i];
	new_we[i2] = we[i];
	/*printf("window %5d  %10d..%10d\n", i2+1, new_ws[i2], new_we[i2]);*/
	i2++;
      }
    }
    if(pli->do_wsplit || (! pli->do_msv)) { 
      /* we could have overlapping windows, merge those that do overlap */
      new_nsurv_fwd = 0;
      ESL_ALLOC(useme, sizeof(int) * nsurv_fwd);
      esl_vec_ISet(useme, nsurv_fwd, FALSE);
      i2 = 0;
      for(i = 0, i2 = 0; i < nsurv_fwd; i++) { 
	useme[i] = TRUE;
	i2 = i+1;
	while((i2 < nsurv_fwd) && ((new_we[i]+1) >= (new_ws[i2]))) { 
	  useme[i2] = FALSE;
	  new_we[i] = new_we[i2]; /* merged i with i2, rewrite end for i */
	  i2++;
	}
	i = i2-1;
      }
      i2 = 0;
      for(i = 0; i < nsurv_fwd; i++) { 
	if(useme[i]) { 
	  new_ws[i2] = new_ws[i];
	  new_we[i2] = new_we[i];
	  i2++;
	}
      }
      nsurv_fwd = i2;
      free(useme);
    }
    free(ws);
    free(we);
    free(wp);
    ws = new_ws;
    we = new_we;
  }
  else { 
    if(ws != NULL) free(ws); ws = NULL;
    if(we != NULL) free(we); we = NULL;
    if(wp != NULL) free(wp); wp = NULL;
  }

  if(survAA != NULL) { 
    for (i = 0; i < Np7_SURV; i++) free(survAA[i]);
    free(survAA);
  }

  *ret_ws   = ws;
  *ret_we   = we;
  *ret_nwin = nsurv_fwd;

  return eslOK;

 ERROR:
  ESL_EXCEPTION(eslEMEM, "Error allocating memory for hit list in pipeline\n");

}

/* Function:  cm_pli_p7EnvelopeDef()
 * Synopsis:  Envelope definition of hits surviving Forward, prior to passing to CYK.
 * Incept:    EPN, Wed Nov 24 13:18:54 2010
 *
 * Purpose:   For each window x, from <ws[x]>..<we[x]>, determine
 *            the envelope boundaries for any hits within it using 
 *            a p7 profile. 
 *
 *            In the SCAN pipeline we may enter this function with
 *            *opt_gm == NULL because we haven't yet read it from the
 *            HMM file. In that case, read it and return it in
 *            <*opt_gm>.  Otherwise, <*opt_gm> is valid upon entering.
 *
 *            If the P-value of any detected envelopes is sufficiently
 *            high (above pli->F5), we skip them (i.e. envelope defn
 *            acts as a filter too). Further, in glocal mode
 *            (<do_glocal>==TRUE) we skip any window for which the
 *            glocal Forward P-value is too high (above pli->F4).
 *
 *            In a normal pipeline run, this function call should be
 *            just after a call to cm_pli_p7Filter() and just
 *            before a call to cm_pli_CMStage().
 *
 * Returns:   <eslOK on success. For all <ret_nenv> envelopes not
 *            filtered out, return their envelope boundaries in
 *            <ret_es> and <ret_ee>.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
cm_pli_p7EnvelopeDef(CM_PIPELINE *pli, P7_OPROFILE *om, P7_BG *bg, float *p7_evparam, const ESL_SQ *sq, int64_t *ws, int64_t *we, int nwin, P7_PROFILE **opt_gm, P7_PROFILE **opt_Rgm, P7_PROFILE **opt_Lgm, P7_PROFILE **opt_Tgm, int64_t **ret_es, int64_t **ret_ee, int *ret_nenv)
{
  int              status;                     
  double           P;                          /* P-value of a hit */
  int              d, i;                       /* counters */
  void            *p;                          /* for ESL_RALLOC */
  int              ali_len, env_len, env_wlen; /* lengths of alignment, envelope, envelope window */
  float            env_sc, env_nullsc;         /* envelope bit score, and null1 score */
  float            sc_for_pvalue;              /* score for window, used for calc'ing P value */
  float            env_edefbias;                /* null2 correction for envelope */
  float            env_sc_for_pvalue;          /* corrected score for envelope, used for calc'ing P value */
  int64_t          wlen;             /* window length of current window */
  int64_t         *es = NULL;        /* [0..nenv-1] envelope start positions */
  int64_t         *ee = NULL;        /* [0..nenv-1] envelope end   positions */
  int              nenv;             /* number of surviving envelopes */
  int              nenv_alloc;       /* current size of es, ee */
  ESL_DSQ         *subdsq;           /* a ptr to the first position of a window */
  ESL_SQ          *seq = NULL;       /* a copy of a window */
  int64_t          estart, eend;     /* envelope start/end positions */
  float            nullsc, filtersc, fwdsc, bcksc;
  P7_PROFILE      *gm  = NULL;        /* a ptr to *opt_gm, for convenience */
  int              do_local_envdef;  /* TRUE if we define envelopes with p7 in local mode, FALSE for glocal */

  /* variables related to forcing first and/or final residue within truncated hits */
  P7_PROFILE      *Rgm = NULL;        /* a ptr to *Ropt_gm, for convenience */
  P7_PROFILE      *Lgm = NULL;        /* a ptr to *Lopt_gm, for convenience */
  P7_PROFILE      *Tgm = NULL;        /* a ptr to *Topt_gm, for convenience */
  float            safe_lfwdsc;       /* a score <= local Forward score, determined via a correction to a Forward score with Rgm or Lgm */
  int              use_Rgm, use_Lgm, use_Tgm; /* only one of these can be TRUE, should we use the specially configured T, R, or L profiles
					       * because we're defining envelopes for hits possibly truncated 5', 3' or both 5' and 3'? */
  float            Rgm_correction;    /* nat score correction for windows and envelopes defined with Rgm */
  float            Lgm_correction;    /* nat score correction for windows and envelopes defined with Lgm */

  if (sq->n == 0) return eslOK;    /* silently skip length 0 seqs; they'd cause us all sorts of weird problems */
  if (nwin == 0)  { 
    *ret_es = es;
    *ret_ee = ee;
    *ret_nenv = 0;
    return eslOK;    /* if there's no envelopes to search in, return */
  }
  do_local_envdef = (pli->do_localenv || (pli->do_local_ends && (pli->rs5term || pli->rs3term))) ? TRUE : FALSE;

  if(pli->do_force_ends && (pli->rs5term || pli->rs3term) && nwin > 1) ESL_FAIL(eslFAIL, pli->errbuf, "researching terminus but more than 1 window exist");

  nenv_alloc = nwin;
  ESL_ALLOC(es, sizeof(int64_t) * nenv_alloc);
  ESL_ALLOC(ee, sizeof(int64_t) * nenv_alloc);
  nenv = 0;
  seq = esl_sq_CreateDigital(sq->abc);

#if DOPRINT
  printf("\nPIPELINE p7EnvelopeDef() %s  %" PRId64 " residues\n", sq->name, sq->n);
#endif

  /* if we're in SCAN mode, we don't yet have the full generic model, read the HMM and create it */
  if (pli->mode == CM_SCAN_MODELS && (do_local_envdef) && (*opt_gm) == NULL) { 
    P7_HMM       *hmm = NULL;
    if (pli->cmfp      == NULL) ESL_FAIL(status, pli->errbuf, "No file available to read HMM from in cm_pli_p7EnvelopeDef()");
    if (pli->cmfp->hfp == NULL) ESL_FAIL(status, pli->errbuf, "No file available to read HMM from in cm_pli_p7EnvelopeDef()");
    if((status = cm_p7_hmmfile_Read(pli->cmfp, pli->abc, om->offs[p7_MOFFSET], &hmm)) != eslOK) ESL_FAIL(status, pli->errbuf, pli->cmfp->errbuf);
    *opt_gm = p7_profile_Create(hmm->M, pli->abc);
    p7_ProfileConfig(hmm, bg, *opt_gm, 100, p7_GLOCAL);

    if(pli->do_force_ends) { 
      if(pli->rs5term) { 
	*opt_Rgm = p7_profile_Clone(*opt_gm);
	p7_ProfileConfig5PrimeTrunc(*opt_Rgm, wlen);
      }
      if(pli->rs3term) { 
	*opt_Lgm = p7_profile_Clone(*opt_gm);
	p7_ProfileConfig3PrimeTrunc(hmm, *opt_Lgm, wlen);
      }
      if(pli->rs5term && pli->rs3term) { 
	*opt_Tgm = p7_profile_Clone(*opt_gm);
	p7_ProfileConfig(hmm, bg, *opt_Tgm, wlen, p7_LOCAL);
	p7_ProfileConfig5PrimeAnd3PrimeTrunc(*opt_Tgm, wlen);
      }  
    }
    p7_hmm_Destroy(hmm);
  }
  gm  = *opt_gm;
  Rgm = *opt_Rgm;
  Lgm = *opt_Lgm;
  Tgm = *opt_Tgm;

  use_Tgm = use_Rgm = use_Lgm = FALSE;
  if((! do_local_envdef) && pli->do_force_ends) { 
    if     (pli->rs5term && pli->rs3term) use_Tgm = TRUE;
    else if(pli->rs5term)                 use_Rgm = TRUE;
    else if(pli->rs3term)                 use_Lgm = TRUE;
  }

  for (i = 0; i < nwin; i++) {
#if DOPRINT    
    printf("p7 envdef win: %4d of %4d [%6" PRId64 "..%6" PRId64 "] rs5term: %d rs3term: %d\n", i, nwin, ws[i], we[i], pli->rs5term, pli->rs3term);
#endif
    if(pli->do_force_ends && pli->rs5term && ws[i] != 1)     ESL_FAIL(eslFAIL, pli->errbuf, "researching 5' terminus but residue 1 outside window");
    if(pli->do_force_ends && pli->rs3term && we[i] != sq->n) ESL_FAIL(eslFAIL, pli->errbuf, "researching 3' terminus but residue L outside window");

    wlen   = we[i]   - ws[i] + 1;
    subdsq = sq->dsq + ws[i] - 1;
    
    /* set up seq object for domaindef function */
    esl_sq_GrowTo(seq, wlen);
    memcpy((void*)(seq->dsq), subdsq, (wlen+1) * sizeof(uint8_t)); 
    seq->dsq[0] = seq->dsq[wlen+1] = eslDSQ_SENTINEL;
    seq->n = wlen;

    p7_bg_SetLength(bg, wlen);
    p7_bg_NullOne(bg, seq->dsq, wlen, &nullsc);

    if(do_local_envdef) { /* we can use optimized matrices and, consequently, p7_domaindef_ByPosteriorHeuristics */
      p7_oprofile_ReconfigLength(om, wlen);
      p7_ForwardParser(seq->dsq, wlen, om, pli->oxf, NULL);
      p7_omx_GrowTo(pli->oxb, om->M, 0, wlen);
      p7_BackwardParser(seq->dsq, wlen, om, pli->oxf, pli->oxb, NULL);
      status = p7_domaindef_ByPosteriorHeuristics (seq, om, pli->oxf, pli->oxb, pli->fwd,  pli->bck, pli->ddef, NULL); 
    }
    else { /* we're defining envelopes in glocal mode, so we need to fill 
	    * generic fwd/bck matrices and pass them to p7_domaindef_GlocalByPosteriorHeuristics() */
      if(use_Tgm) { 
	/* no length reconfiguration necessary */
	p7_gmx_GrowTo(pli->gxf, Tgm->M, wlen);
	p7_GForward (seq->dsq, wlen, Tgm, pli->gxf, &fwdsc);
	/*printf("Tfwdsc: %.4f\n", fwdsc);*/
	/* we use local Fwd statistics to determine statistical significance of this score,
	 * it has already had basically a 1/log(M*(M+1)) penalty for equiprobable local begins and ends */
	sc_for_pvalue = (fwdsc - nullsc) / eslCONST_LOG2;
	P = esl_exp_surv (sc_for_pvalue,  p7_evparam[CM_p7_LFTAU],  p7_evparam[CM_p7_LFLAMBDA]);
      }
      else if(use_Rgm) { 
	p7_ReconfigLength5PrimeTrunc(Rgm, wlen);
	p7_gmx_GrowTo(pli->gxf, Rgm->M, wlen);
	p7_GForward (seq->dsq, wlen, Rgm, pli->gxf, &fwdsc);
	/*printf("Rfwdsc: %.4f\n", fwdsc);*/
	/* we use local Fwd statistics to determine significance of the score, but we need to correct for lack of equiprobable ends in fwdsc */
	Rgm_correction = (log(2./(Rgm->M * (Rgm->M+1))) - log(1./Rgm->M)); /* GForward penalized 0. for ends and log(1/Rgm->M) for begins into any state */
	safe_lfwdsc = fwdsc + Rgm_correction;
	/* safe_lfwdsc is now less than or equal to the local foward score we would get with a truly local model */
	sc_for_pvalue = (safe_lfwdsc - nullsc) / eslCONST_LOG2;
	P = esl_exp_surv (sc_for_pvalue,  p7_evparam[CM_p7_LFTAU],  p7_evparam[CM_p7_LFLAMBDA]);
      }
      else if(use_Lgm) { 
	p7_ReconfigLength3PrimeTrunc(Lgm, wlen);
	p7_gmx_GrowTo(pli->gxf, Lgm->M, wlen);
	p7_GForward (seq->dsq, wlen, Lgm, pli->gxf, &fwdsc);
	/*printf("Lfwdsc: %.4f\n", fwdsc);*/
	/* we use local Fwd statistics to determine significance of the score, but we need to correct for lack of equiprobable begins and ends in fwdsc */
	Lgm_correction = log(2./(Lgm->M * (Lgm->M+1))); /* GForward penalized 0. for ends and 0. for begins into M1 */
	safe_lfwdsc = fwdsc + Lgm_correction; 
	/* safe_lfwdsc is now less than or equal to the local foward score we would get with a truly local model */
	sc_for_pvalue = (safe_lfwdsc - nullsc) / eslCONST_LOG2;
	P = esl_exp_surv (sc_for_pvalue,  p7_evparam[CM_p7_LFTAU],  p7_evparam[CM_p7_LFLAMBDA]);
      }
      else { /* normal case, not looking for truncated hits */
	p7_ReconfigLength(gm, wlen);
	p7_gmx_GrowTo(pli->gxf, gm->M, wlen);
	p7_GForward (seq->dsq, wlen, gm, pli->gxf, &fwdsc);
	/*printf(" fwdsc: %.4f\n", fwdsc);*/
	sc_for_pvalue = (fwdsc - nullsc) / eslCONST_LOG2;
	/* Does this score exceed our glocal forward filter threshold? If not, move on to next seq */
	P = esl_exp_surv (sc_for_pvalue,  p7_evparam[CM_p7_GFMU],  p7_evparam[CM_p7_GFLAMBDA]);
      }


#if DOPRINT	
      if(P > pli->F4) { 
	printf("KILLED   window %5d [%10" PRId64 "..%10" PRId64 "]          gFwd      %6.2f bits  P %g\n", i, ws[i], we[i], sc_for_pvalue, P);
      }
#endif      
      if(P > pli->F4) continue;

      pli->acct[pli->cur_pass_idx].n_past_gfwd++;
      pli->acct[pli->cur_pass_idx].pos_past_gfwd += wlen;

#if DOPRINT	
      printf("SURVIVOR window %5d [%10" PRId64 "..%10" PRId64 "] survived gFwd      %6.2f bits  P %g\n", i, ws[i], we[i], sc_for_pvalue, P);
#endif      

      if(pli->do_gfwdbias) {
	/* calculate bias filter score for entire window */
	p7_bg_FilterScore(bg, seq->dsq, wlen, &filtersc);
	sc_for_pvalue = (fwdsc - filtersc) / eslCONST_LOG2;
	P = esl_exp_surv (sc_for_pvalue,  p7_evparam[CM_p7_GFMU],  p7_evparam[CM_p7_GFLAMBDA]);
	if(P > pli->F4b) continue;
#if DOPRINT	
	printf("SURVIVOR window %5d [%10" PRId64 "..%10" PRId64 "] survived gFwdBias  %6.2f bits  P %g\n", i, ws[i], we[i], sc_for_pvalue, P);
#endif 
	pli->acct[pli->cur_pass_idx].n_past_gfwdbias++;
	pli->acct[pli->cur_pass_idx].pos_past_gfwdbias += wlen;
      }
      if(pli->do_time_F4) continue; 

      /* this block needs to match up with if..else if...else if...else block calling p7_GForward above */
      if(use_Tgm) { 
	/* no length reconfiguration necessary */
	p7_gmx_GrowTo(pli->gxb, Tgm->M, wlen);
	p7_GBackward(seq->dsq, wlen, Tgm, pli->gxb, &bcksc);
	if((status = p7_domaindef_GlocalByPosteriorHeuristics(seq, Tgm, pli->gxf, pli->gxb, pli->gfwd, pli->gbck, pli->ddef, pli->do_null2)) != eslOK) ESL_FAIL(status, pli->errbuf, "unexpected failure during glocal envelope defn"); 
	/*printf("Tbcksc: %.4f\n", bcksc);*/
      }
      else if(use_Rgm) { 
	p7_gmx_GrowTo(pli->gxb, Rgm->M, wlen);
	p7_GBackward(seq->dsq, wlen, Rgm, pli->gxb, &bcksc);
	if((status = p7_domaindef_GlocalByPosteriorHeuristics(seq, Rgm, pli->gxf, pli->gxb, pli->gfwd, pli->gbck, pli->ddef, pli->do_null2)) != eslOK) ESL_FAIL(status, pli->errbuf, "unexpected failure during glocal envelope defn");; 
	/*printf("Rbcksc: %.4f\n", bcksc);*/
      }
      else if(use_Lgm) { 
	p7_gmx_GrowTo(pli->gxb, Lgm->M, wlen);
	p7_GBackward(seq->dsq, wlen, Lgm, pli->gxb, &bcksc);
	if((status = p7_domaindef_GlocalByPosteriorHeuristics(seq, Lgm, pli->gxf, pli->gxb, pli->gfwd, pli->gbck, pli->ddef, pli->do_null2)) != eslOK) ESL_FAIL(status, pli->errbuf, "unexpected failure during glocal envelope defn");
	/*printf("Lbcksc: %.4f\n", bcksc);*/
      }
      else { /* normal case, not looking for truncated hits */
	p7_gmx_GrowTo(pli->gxb, gm->M, wlen);
	p7_GBackward(seq->dsq, wlen, gm, pli->gxb, &bcksc);
	if((status = p7_domaindef_GlocalByPosteriorHeuristics(seq, gm, pli->gxf, pli->gxb, pli->gfwd, pli->gbck, pli->ddef, pli->do_null2)) != eslOK) ESL_FAIL(status, pli->errbuf, "unexpected failure during glocal envelope defn");
	/*printf(" bcksc: %.4f\n", bcksc);*/
      }
    }
    
    if (status != eslOK) ESL_FAIL(status, pli->errbuf, "envelope definition workflow failure"); /* eslERANGE can happen */
    if (pli->ddef->nregions   == 0)  continue; /* score passed threshold but there's no discrete domains here       */
    if (pli->ddef->nenvelopes == 0)  continue; /* rarer: region was found, stochastic clustered, no envelopes found */
    
    /* For each domain found in the p7_domaindef_*() function, determine if it passes our criteria */
    for(d = 0; d < pli->ddef->ndom; d++) { 
      
      if(do_local_envdef) { /* we called p7_domaindef_ByPosteriorHeuristics() above, which fills pli->ddef->dcl[d].ad, but we don't need it */
	p7_alidisplay_Destroy(pli->ddef->dcl[d].ad);
      }

      env_len = pli->ddef->dcl[d].jenv - pli->ddef->dcl[d].ienv +1;
      env_sc  = pli->ddef->dcl[d].envsc;
      
      /* Make a correction to the score based on the fact that our envsc 
       * is just a segment of the full window, there are three options here: */
      if(pli->do_nocorr) { 
	/* no correction */
	env_sc_for_pvalue = (env_sc - nullsc) / eslCONST_LOG2;
      }
      else if(pli->do_oldcorr) { 
	/* correction from nhmmer's p7_pipeline.c::p7_Pipeline_LongTarget() */
	env_len = pli->ddef->dcl[d].jenv - pli->ddef->dcl[d].ienv +1;
	ali_len = pli->ddef->dcl[d].jali - pli->ddef->dcl[d].iali +1;
	env_wlen = ESL_MIN(pli->cmW, wlen); //see notes, ~/notebook/20100716_hmmer_score_v_eval_bug/, end of Thu Jul 22 13:36:49 EDT 2010
	env_sc  = pli->ddef->dcl[d].envsc;
	/* For these modifications, see notes, ~/notebook/20100716_hmmer_score_v_eval_bug/, end of Thu Jul 22 13:36:49 EDT 2010 */
	env_sc -= 2 * log(2. / (wlen+2)) +   (env_len-ali_len) * log((float) wlen / (wlen+2));
	env_sc += 2 * log(2. / (env_wlen+2)) ;
	/* handle extremely rare case that the env_len is actually larger than om->max_length */
	/* (I don't think its very rare for dcmsearch b/c MSV P value threshold is so high) */
	env_sc +=  (ESL_MAX(env_wlen, env_len) - ali_len) * log((float) env_wlen / (float) (env_wlen+2));
	env_nullsc = (float) env_wlen * log((float)env_wlen/(env_wlen+1)) + log(1./(env_wlen+1));
	env_sc_for_pvalue = (env_sc - env_nullsc) / eslCONST_LOG2;
      }
      else { 
	/* default correction,  from hmmsearch's p7_pipeline() */
	/* here is the p7_pipeline code, verbatim: 
	 *  Ld = hit->dcl[d].jenv - hit->dcl[d].ienv + 1;
	 *  hit->dcl[d].bitscore = hit->dcl[d].envsc + (sq->n-Ld) * log((float) sq->n / (float) (sq->n+3)); 
	 *  hit->dcl[d].dombias  = (pli->do_null2 ? p7_FLogsum(0.0, log(bg->omega) + hit->dcl[d].domcorrection) : 0.0); 
	 *  hit->dcl[d].bitscore = (hit->dcl[d].bitscore - (nullsc + hit->dcl[d].dombias)) / eslCONST_LOG2; 
	 *  hit->dcl[d].pvalue   = esl_exp_surv (hit->dcl[d].bitscore,  om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);
	 *  
	 *  and here is the same code, simplified with our different names for the variables, etc 
	 * (we don't use hit->dcl the way p7_pipeline does after this):  */
	env_sc            = env_sc + (wlen - env_len) * log((float) wlen / (float) (wlen+3)); /* NATS, for the moment... */
	env_edefbias       = (pli->do_null2 ? p7_FLogsum(0.0, log(bg->omega) + pli->ddef->dcl[d].domcorrection) : 0.0); /* NATS, and will stay so */
	env_sc_for_pvalue = (env_sc - (nullsc + env_edefbias)) / eslCONST_LOG2; /* now BITS, as it should be */
      }
      if(use_Rgm) env_sc_for_pvalue += Rgm_correction / eslCONST_LOG2; /* glocal env def penalized 0. for ends and log(1/Lgm->M) for begins into any state */
      if(use_Lgm) env_sc_for_pvalue += Lgm_correction / eslCONST_LOG2; /* glocal env def penalized 0. for ends and 0. for begins into M1 */

      if(do_local_envdef || use_Tgm || use_Rgm || use_Lgm) P = esl_exp_surv (env_sc_for_pvalue,  p7_evparam[CM_p7_LFTAU], p7_evparam[CM_p7_LFLAMBDA]);
      else                                                 P = esl_exp_surv (env_sc_for_pvalue,  p7_evparam[CM_p7_GFMU],  p7_evparam[CM_p7_GFLAMBDA]);
      /***************************************************/
      
      /* check if we can skip this envelope based on its P-value or bit-score */
      if(P > pli->F5) continue;
    
      /* Define envelope to search with CM:
       *
       * if (! pli->do_pad): just pass existing envelope
       * from the p7DomainDef() function.
       *
       * if (pli->do_pad):
       * 
       *                envelope/domain
       *   ----------xxxxxxxxxxxxxxxx---------
       *   |<-----------------------|
       *                W
       *             |---------------------->|
       *                         W
       *
       * On rare occasions, the envelope size may exceed W, in that case we
       * define dstarts-W-1..dend+W-1 as the envelope.
       *
       */
      if(! pli->do_pad) { 
	estart = pli->ddef->dcl[d].ienv;
	eend   = pli->ddef->dcl[d].jenv;
      }
      else { /* pli->do_pad is TRUE */
	if(env_len <= pli->cmW) { /* envelope size is less than or equal to W */
	  estart = ESL_MAX(1,    pli->ddef->dcl[d].jenv - (pli->cmW-1));
	  eend   = ESL_MIN(wlen, pli->ddef->dcl[d].ienv + (pli->cmW-1));
	}
	else { /* envelope size > W, pad W-1 residues to each side */
	  estart = ESL_MAX(1,    pli->ddef->dcl[d].ienv - (int) ((pli->cmW - 1)));
	  eend   = ESL_MIN(wlen, pli->ddef->dcl[d].jenv + (int) ((pli->cmW - 1)));
	}
      }
#if DOPRINT
      printf("SURVIVOR envelope     [%10" PRId64 "..%10" PRId64 "] survived F5       %6.2f bits  P %g\n", pli->ddef->dcl[d].ienv + ws[i] - 1, pli->ddef->dcl[d].jenv + ws[i] - 1, env_sc_for_pvalue, P);
#endif
      pli->acct[pli->cur_pass_idx].n_past_edef++;
      pli->acct[pli->cur_pass_idx].pos_past_edef += env_len;

      /* if we're doing a bias filter on envelopes - check if we skip envelope due to that */
      if(pli->do_edefbias) {
	if(pli->do_envwinbias) { 
	  /* calculate bias filter score for entire window */
	  p7_bg_FilterScore(bg, seq->dsq, wlen, &filtersc);
	}
	else {
	  /* calculate bias filter score for only the residues in the envelope */
	  p7_bg_FilterScore(bg, seq->dsq + pli->ddef->dcl[d].ienv - 1, env_len, &filtersc);
	  /* add in null1-only score for residues in the window but outside the envelope */
	  filtersc += (float) (wlen-env_len) * logf(bg->p1);
	}
	env_sc_for_pvalue = (env_sc - filtersc) / eslCONST_LOG2;
	if(do_local_envdef) P = esl_exp_surv (env_sc_for_pvalue,  p7_evparam[CM_p7_LFTAU], p7_evparam[CM_p7_LFLAMBDA]);
	else                P = esl_exp_surv (env_sc_for_pvalue,  p7_evparam[CM_p7_GFMU],  p7_evparam[CM_p7_GFLAMBDA]);
	if(P > pli->F5b) continue;
      }
#if DOPRINT
      printf("SURVIVOR envelope     [%10" PRId64 "..%10" PRId64 "] survived F5-bias  %6.2f bits  P %g\n", pli->ddef->dcl[d].ienv + ws[i] - 1, pli->ddef->dcl[d].jenv + ws[i] - 1, env_sc_for_pvalue, P);
#endif
      pli->acct[pli->cur_pass_idx].n_past_edefbias++;
      pli->acct[pli->cur_pass_idx].pos_past_edefbias += env_len;
	
      if(pli->do_time_F5) { continue; }

      /* if we get here, the envelope has survived, add it to the growing list */
      if((nenv+1) == nenv_alloc) { 
	nenv_alloc *= 2;
	ESL_RALLOC(es, p, sizeof(int64_t) * nenv_alloc);
	ESL_RALLOC(ee, p, sizeof(int64_t) * nenv_alloc);
      }
      es[nenv] = estart + ws[i] - 1;
      ee[nenv] = eend   + ws[i] - 1;
      nenv++;
    }
    pli->ddef->ndom = 0; /* reset for next use */
  }

  /* clean up, set return variables, and return */
  if(seq != NULL) esl_sq_Destroy(seq);

  *ret_es = es;
  *ret_ee = ee;
  *ret_nenv = nenv;

  return eslOK;

 ERROR:
  ESL_EXCEPTION(eslEMEM, "Error: out of memory");
}


/* Function:  cm_pli_CMStage()
 * Synopsis:  Final stage of pipeline, CYK filter, then final scoring with Inside.
 * Incept:    EPN, Sun Nov 28 13:47:31 2010
 *
 * Purpose:   For each envelope x, from <es[x]>..<ee[x]>, run 
 *            CYK to see if any hits above threshold exist, and
 *            if so pass to Inside for final hit definition.
 *
 *            In a normal pipeline run, this function call should be
 *            just after a call to cm_pli_p7EnvelopeDef().
 *
 *            If pli->mode is CM_SCAN_MODELS, it's possible that we
 *            haven't yet read our CM from the file. This is true when
 *            (*opt_cm == NULL). If so, we read the CM from the file
 *            after positioning it to position <cm_offset> and
 *            configure the CM after setting cm->config_opts to
 *            <pli->cm_config_opts>. In this case, <*opt_cmcons>
 *            should also be NULL and we create a CMConsensus_t object
 *            and return it in <*opt_cmcons>. Otherwise, if <*opt_cm>
 *            is valid (non-NULL), <*opt_cmcons> should be as well.
 *
 * Returns:   <eslOK> on success. If a significant hit is obtained,
 *            its information is added to the growing <hitlist>.
 *
 *            <eslEINVAL> if (in a scan pipeline) we're supposed to
 *            set GA/TC/NC bit score thresholds but the model doesn't
 *            have any.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
cm_pli_CMStage(CM_PIPELINE *pli, off_t cm_offset, const ESL_SQ *sq, int64_t *es, int64_t *ee, int nenv, CM_TOPHITS *hitlist, CM_t **opt_cm, CMConsensus_t **opt_cmcons)
{
  int              status;
  char             errbuf[cmERRBUFSIZE];   /* for error messages */
  CM_HIT          *hit     = NULL;         /* ptr to the current hit output data   */
  float            cyksc, inssc, finalsc;  /* bit scores                           */
  char             cykmode, insmode;       /* marginal alignment modes */
  double           P;                      /* P-value of a hit */
  int              i, h;                   /* counters */
  int              nhit;                   /* number of hits reported */
  double           save_tau;               /* CM's tau upon entering function */
  float            save_cp9b_thresh1;      /* cm->cp9b's thresh1 upon entering function */
  float            save_cp9b_thresh2;      /* cm->cp9b's thresh2 upon entering function */
  int64_t          cyk_envi, cyk_envj;     /* cyk_envi..cyk_envj is new envelope as defined by CYK hits */
  float            cyk_env_cutoff;         /* bit score cutoff for envelope redefinition */
  int              do_hbanded_filter_scan; /* use HMM bands for filter stage? */
  int              do_hbanded_final_scan;  /* use HMM bands for final stage? */
  int              do_qdb_or_nonbanded_filter_scan; /* use QDBs or no bands for filter stage (! do_hbanded_filter_scan) */
  int              do_qdb_or_nonbanded_final_scan;  /* use QDBs or no bands for final  stage (! do_hbanded_final_scan) */
  int              do_final_greedy;        /* TRUE to use greedy hit resolution in final stage, FALSE not to */
  CM_t            *cm = NULL;              /* ptr to *opt_cm, for convenience only */
  CP9Bands_t      *scan_cp9b = NULL;       /* a copy of the HMM bands derived in the final CM search stage, if its HMM banded */
  float            hbmx_Mb;                /* approximate size in Mb for HMM banded matrix for current hit */
  int              do_trunc;               /* TRUE if we are allowing L and/or R marginal alignments, FALSE if not */
  int              check_fcyk_beta;        /* TRUE if we just read a CM and we need to check if its beta1 == pli->fcyk_beta */
  int              check_final_beta;       /* TRUE if we just read a CM and we need to check if its beta2 == pli->final_beta */

  if (sq->n == 0) return eslOK;    /* silently skip length 0 seqs; they'd cause us all sorts of weird problems */
  if (nenv == 0)  return eslOK;    /* if there's no envelopes to search in, return */

  do_trunc = (pli->do_trunc_ends && (pli->tro->allowR || pli->tro->allowL)) ? TRUE : FALSE;

  /* if we're in SCAN mode, and we don't yet have a CM, read it and configure it */
  if (pli->mode == CM_SCAN_MODELS) { 
    if(*opt_cm == NULL) { 
#ifdef HMMER_THREADS
      /* lock the mutex to prevent other threads from reading the file at the same time */
      if (pli->cmfp->syncRead) { 
	if (pthread_mutex_lock (&pli->cmfp->readMutex)      != 0) ESL_FAIL(eslESYS, pli->errbuf, "mutex lock failed");
      }
#endif
      cm_file_Position(pli->cmfp, cm_offset);
      if((status = cm_file_Read(pli->cmfp, FALSE, &(pli->abc), opt_cm)) != eslOK) ESL_FAIL(status, pli->errbuf, pli->cmfp->errbuf);
#ifdef HMMER_THREADS
      if (pli->cmfp->syncRead) { 
	if (pthread_mutex_unlock (&pli->cmfp->readMutex)      != 0) ESL_EXCEPTION(eslESYS, "mutex unlock failed");
      }
#endif    
      /*printf("CONFIGURING CM: %s (%d)\n", (*opt_cm)->name, pli->cm_config_opts);*/
      (*opt_cm)->config_opts = pli->cm_config_opts;
      /* check if we need to recalculate QDBs prior to building the scan matrix in cm_Configure() 
       * (we couldn't do this until we read the CM file to find out what cm->qdbinfo->beta1/beta2 were 
       */
      check_fcyk_beta  = (pli->fcyk_cm_search_opts  & CM_SEARCH_QDB) ? TRUE : FALSE;
      check_final_beta = (pli->final_cm_search_opts & CM_SEARCH_QDB) ? TRUE : FALSE;
      if((status = CheckCMQDBInfo(cm->qdbinfo, pli->fcyk_beta, check_fcyk_beta, pli->final_beta, check_final_beta)) != eslOK) { 
	(*opt_cm)->config_opts   |= CM_CONFIG_QDB;
	(*opt_cm)->qdbinfo->beta1 = pli->fcyk_beta;
	(*opt_cm)->qdbinfo->beta2 = pli->final_beta;
      }
      /* else we don't have to change (*opt_cm)->qdbinfo->beta1/beta2 */

      if((status = cm_Configure(*opt_cm, pli->errbuf)) != eslOK) return status;
      /* update the pipeline about the model */
      if((status = cm_pli_NewModel(pli, CM_NEWMODEL_CM, *opt_cm, (*opt_cm)->clen, (*opt_cm)->W, 
				   NULL, NULL, pli->cur_cm_idx)) /* om, bg */
	 != eslOK); 
    }
    if(*opt_cmcons == NULL) { 
      if((status = CreateCMConsensus(*opt_cm, (*opt_cm)->abc, 3.0, 1.0, opt_cmcons))!= eslOK) ESL_FAIL(status, pli->errbuf, "In cm_pli_CMStage() failed to create CMConsensus data structure.\n");
    }
  }
  else { /* not SCAN mode, *opt_cm and *opt_cmconsensus should be valid */
    if(opt_cmcons == NULL || *opt_cmcons == NULL) ESL_FAIL(eslEINCOMPAT, pli->errbuf, "Entered cm_pli_CMStage() with invalid CMConsensus structure"); 
    if(opt_cm     == NULL || *opt_cm     == NULL) ESL_FAIL(eslEINCOMPAT, pli->errbuf, "Entered cm_pli_CMStage() with invalid CM"); 
  }

  cm = *opt_cm;
  save_tau          = cm->tau;
  save_cp9b_thresh1 = cm->cp9b->thresh1;
  save_cp9b_thresh2 = cm->cp9b->thresh2;
  do_final_greedy = (pli->final_cm_search_opts & CM_SEARCH_CMNOTGREEDY) ? FALSE : TRUE;

  nhit = hitlist->N;
  /* Determine bit score cutoff for CYK envelope redefinition, 
   * any residue that exists in a CYK hit that reaches this threshold will be included
   * in the redefined envelope, any that doesn't will not be.
   */
  cyk_env_cutoff = cm->expA[pli->fcyk_cm_exp_mode]->mu_extrap + (log(pli->F6env) / (-1 * cm->expA[pli->fcyk_cm_exp_mode]->lambda));

#if DOPRINT
  printf("\nPIPELINE CMStage() %s  %" PRId64 " residues\n", sq->name, sq->n);
#endif
  

  for (i = 0; i < nenv; i++) {
    cm->tau           = save_tau;
    cm->cp9b->thresh1 = save_cp9b_thresh1;
    cm->cp9b->thresh2 = save_cp9b_thresh2;
      
#if DOPRINT
    printf("\nSURVIVOR Envelope %5d [%10ld..%10ld] being passed to CYK   allowR: %d allowL: %d force_i0_RT: %d force_j0_LT: %d\n", i, es[i], ee[i],
	   pli->tro->allowR, pli->tro->allowL, pli->tro->force_i0_RT, pli->tro->force_j0_LT); 
#endif
    do_hbanded_filter_scan          = (pli->fcyk_cm_search_opts  & CM_SEARCH_HBANDED) ? TRUE  : FALSE;
    do_qdb_or_nonbanded_filter_scan = (pli->fcyk_cm_search_opts  & CM_SEARCH_HBANDED) ? FALSE : TRUE;
    do_hbanded_final_scan           = (pli->final_cm_search_opts & CM_SEARCH_HBANDED) ? TRUE  : FALSE;
    do_qdb_or_nonbanded_final_scan  = (pli->final_cm_search_opts & CM_SEARCH_HBANDED) ? FALSE : TRUE;

    if(pli->do_cyk) { 
      /******************************************************************************/
      /* Filter 4: CYK with CM */
      cm->search_opts  = pli->fcyk_cm_search_opts;
      cm->tau          = pli->fcyk_tau;

      /* Two options, dependent on pli->fcyk_cm_search_opts.
       * 1. Do HMM banded CYK.
       * 2. Do QDB or non-banded CYK.
       *
       * If we try option 1 and it will require too much memory (> 1 Gb) 
       * we fail over to option 2.
       */
      if(do_hbanded_filter_scan) { /* get HMM bands */
	/* Calculate HMM bands. We'll tighten tau and recalculate bands until 
	 * the resulting HMM banded matrix is under our size limit.
	 */
	while(1) { 
	  if((status = cp9_Seq2Bands(cm, pli->errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, sq->dsq, es[i], ee[i], cm->cp9b, TRUE, pli->tro, 0)) != eslOK) return status;
	  if(do_trunc) {
	    if((status = cm_tr_hb_mx_SizeNeeded(cm, pli->errbuf, cm->cp9b, ee[i]-es[i]+1, NULL, NULL, NULL, NULL, &hbmx_Mb)) != eslOK) return status; 
	  }
	  else { 
	    if((status = cm_hb_mx_SizeNeeded   (cm, pli->errbuf, cm->cp9b, ee[i]-es[i]+1, NULL, &hbmx_Mb)) != eslOK) return status; 
	  }
	  if(hbmx_Mb < pli->hb_size_limit) break; /* our matrix will be small enough, break out of while(1) */
	  if(cm->tau > 0.01)               break; /* our bands have gotten too tight, break out of while(1) */
	  /* printf("CYK 0 tau: %10g  hbmx_Mb: %10.2f\n", cm->tau, hbmx_Mb);*/
	  cm->tau *= pli->xtau;
	  /* cm->cp9b->thresh1 *= 2.; */
	  /* cm->cp9b->thresh2 -= (1.0-cp9b->thresh2); */
	  /* cm->cp9b->thresh1 = ESL_MIN(0.25, cm->cp9b->thresh1); */
	  /* cm->cp9b->thresh2 = ESL_MAX(0.25, cm->cp9b->thresh2); */
	}	  
	/* printf("CYK 1 tau: %10g  hbmx_Mb: %10.2f\n", cm->tau, hbmx_Mb); */

	/* reset search opts */
	cm->search_opts  = pli->fcyk_cm_search_opts;

	if(do_trunc) { 
	  status = TrCYKScanHB(cm, errbuf, cm->trhbmx, pli->hb_size_limit, pli->tro, sq->dsq, es[i], ee[i], 
			       0.,                                  /* minimum score to report, irrelevant */
			       NULL,                                /* hitlist to add to, irrelevant here */
			       pli->do_null3,                       /* do the NULL3 correction? */
			       cyk_env_cutoff,                      /* bit score == F6env P value, cutoff for envelope redefinition */
			       (pli->do_cykenv) ? &cyk_envi : NULL, /* envelope start, derived from CYK hits */
			       (pli->do_cykenv) ? &cyk_envj : NULL, /* envelope stop,  derived from CYK hits */
			       &cykmode,                            /* best alignment mode */
			       &cyksc);                             /* best score, irrelevant here */
	  /*printf("\n\n!!! Returned from TrCYKScanHB(): %.4f bits\n\n\n", cyksc);*/
	}
	else { 
	  status = FastCYKScanHB(cm, errbuf, cm->hbmx, pli->hb_size_limit, sq->dsq, es[i], ee[i], 
				 0.,                                  /* minimum score to report, irrelevant */
				 NULL,                                /* hitlist to add to, irrelevant here */
				 pli->do_null3,                       /* do the NULL3 correction? */
				 cyk_env_cutoff,                      /* bit score == F6env P value, cutoff for envelope redefinition */
				 (pli->do_cykenv) ? &cyk_envi : NULL, /* envelope start, derived from CYK hits */
				 (pli->do_cykenv) ? &cyk_envj : NULL, /* envelope stop,  derived from CYK hits */
				 &cyksc);                             /* best score, irrelevant here */
	}
	if     (status == eslERANGE) { do_qdb_or_nonbanded_filter_scan = TRUE; }
	else if(status != eslOK)     { printf("ERROR: %s\n", errbuf); return status; }
      }
      if(do_qdb_or_nonbanded_filter_scan) { /* careful, different from just an 'else', b/c we may have just set this as true if status == eslERANGE */
	/* make sure we have a CM_SCAN_MX for this CM, if not build one,
	 * this should be okay because we should only very rarely need to do this */
	if(cm->smx == NULL) { printf("FIX ME! ADD TRUNC VERSION! cm->smx is NULL, probably overflow sized hb mx\n"); continue; }
	pli->acct[pli->cur_pass_idx].n_overflow_fcyk++;
	if((status = FastCYKScan(cm, errbuf, cm->smx, SMX_QDB1_TIGHT, sq->dsq, es[i], ee[i],
				 0.,                                 /* minimum score to report, irrelevant */
				 NULL,                               /* hitlist to add to, irrelevant here */
				 pli->do_null3,                      /* do the NULL3 correction? */
				 cyk_env_cutoff,                     /* bit score == envF6 P value, cutoff for envelope redefinition */
				 (pli->do_cykenv) ? &cyk_envi : NULL, /* envelope start, derived from CYK hits */
				 (pli->do_cykenv) ? &cyk_envj : NULL, /* envelope stop,  derived from CYK hits */
				 NULL,                               /* ret_vsc, irrelevant here */
				 &cyksc)) != eslOK) { 
	  printf("ERROR: %s\n", errbuf); return status;  }
      }
      /* update envelope boundaries, if nec */
      if(pli->do_cykenv && (cyk_envi != -1 && cyk_envj != -1)) { 
	es[i] = cyk_envi;
	ee[i] = cyk_envj;
      }
      P = esl_exp_surv(cyksc, cm->expA[pli->fcyk_cm_exp_mode]->mu_extrap, cm->expA[pli->fcyk_cm_exp_mode]->lambda);
      if (P > pli->F6) continue;
      /******************************************************************************/
    }	
    pli->acct[pli->cur_pass_idx].n_past_cyk++;
    pli->acct[pli->cur_pass_idx].pos_past_cyk += ee[i] - es[i] + 1;

#if DOPRINT
    printf("SURVIVOR envelope     [%10" PRId64 "..%10" PRId64 "] survived CYK       %6.2f bits  P %g\n", es[i], ee[i], cyksc, P);
#endif
    if(pli->do_time_F6) continue; 

    /******************************************************************************/
    /* Final stage: Inside/CYK with CM, report hits to our hitlist */
    cm->search_opts  = pli->final_cm_search_opts;
    cm->tau          = pli->final_tau;
        
    /*******************************************************************
     * Determine if we're doing a HMM banded scan, if so, we may already have HMM bands 
     * if our CYK filter also used them. 
     *******************************************************************/
    if(do_hbanded_final_scan) { /* use HMM bands */
	/* Calculate HMM bands. We'll tighten tau and recalculate bands until 
	 * the resulting HMM banded matrix is under our size limit.
	 */
	while(1) { 
	  if((status = cp9_Seq2Bands(cm, pli->errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, sq->dsq, es[i], ee[i], cm->cp9b, TRUE, pli->tro, 0)) != eslOK) return status;
	  if(do_trunc) { 
	    if((status = cm_tr_hb_mx_SizeNeeded(cm, pli->errbuf, cm->cp9b, ee[i]-es[i]+1, NULL, NULL, NULL, NULL, &hbmx_Mb)) != eslOK) return status; 
	  }
	  else { 
	    if((status = cm_hb_mx_SizeNeeded   (cm, pli->errbuf, cm->cp9b, ee[i]-es[i]+1, NULL, &hbmx_Mb)) != eslOK) return status; 
	  }
	  if(hbmx_Mb < pli->hb_size_limit) break; /* our matrix will be small enough, break out of while(1) */
	  if(cm->tau > 0.01)               break; /* our bands have gotten too tight, break out of while(1) */
	  /* printf("INS 0 tau: %10g  hbmx_Mb: %10.2f\n", cm->tau, hbmx_Mb); */
	  cm->tau *= pli->xtau;
	}	  
	/* printf("INS 1 tau: %10g  hbmx_Mb: %10.2f\n", cm->tau, hbmx_Mb); */
	if(cm->search_opts & CM_SEARCH_INSIDE) { /* final algorithm is HMM banded Inside */
	  if(do_trunc) { 
	    status = FTrInsideScanHB(cm, errbuf, cm->trhbmx, pli->hb_size_limit, pli->tro, sq->dsq, es[i], ee[i], 
				     pli->T,            /* minimum score to report */
				     hitlist,           /* our hitlist */
				     pli->do_null3,     /* do the NULL3 correction? */
				     0., NULL, NULL,    /* redefined envelope cutoff,start,end, irrelevant here */
				     &insmode,          /* best mode, irrelevant here (it's also stored in any reported hits) */
				     &inssc);           /* best score, irrelevant here */
	    /*printf("\n\n!!! Returned from FTrInsideScanHB(): %.4f bits\n\n\n", inssc);*/
	  }
	  else { 
	    status = FastFInsideScanHB(cm, errbuf, cm->hbmx, pli->hb_size_limit, sq->dsq, es[i], ee[i], 
				       pli->T,            /* minimum score to report */
				       hitlist,           /* our hitlist */
				       pli->do_null3,     /* do the NULL3 correction? */
				       0., NULL, NULL,    /* redefined envelope cutoff,start,end, irrelevant here */
				       &inssc);           /* best score, irrelevant here */
	  }
	  /* if status == eslERANGE: HMM banded scan was skipped b/c mx needed to be too large, 
	   * we'll repeat the scan with QDBs or without bands below */
	  if     (status == eslERANGE) { do_qdb_or_nonbanded_final_scan = TRUE; }
	  else if(status != eslOK)     { printf("ERROR: %s\n", errbuf); return status; }
	}
	else { /* final algorithm is HMM banded CYK */
	  /*printf("calling HMM banded CYK scan\n");*/
	  if(do_trunc) { 
	    status = TrCYKScanHB(cm, errbuf, cm->trhbmx, pli->hb_size_limit, pli->tro, sq->dsq, es[i], ee[i], 
				 pli->T,                             /* minimum score to report */
				 hitlist,                            /* our hitlist */
				 pli->do_null3,                      /* do the NULL3 correction? */
				 0., NULL, NULL,                     /* envelope redefinition parameters, irrelevant here */
				 &cykmode,                           /* best alignment mode, irrevelant here (it's also stored in any reported hits) */
				 &cyksc);                            /* best score, irrelevant here */
	    /*printf("\n\n!!! Returned from TrCYKScanHB(): %.4f bits\n\n\n", cyksc);*/
	  }
	  else { 
	    status = FastCYKScanHB(cm, errbuf, cm->hbmx, pli->hb_size_limit, sq->dsq, es[i], ee[i], 
				   pli->T,                             /* minimum score to report */
				   hitlist,                            /* our hitlist */
				   pli->do_null3,                      /* do the NULL3 correction? */
				   0., NULL, NULL,                     /* envelope redefinition parameters, irrelevant here */
				   &cyksc);                            /* best score, irrelevant here */
	  }
	  /* if status == eslERANGE: HMM banded scan was skipped b/c mx needed to be too large, 
	   * we'll repeat the scan with QDBs or without bands below */
	  if     (status == eslERANGE) { do_qdb_or_nonbanded_final_scan = TRUE; }
	  else if(status != eslOK)     { printf("ERROR: %s\n", errbuf); return status; }
	}
    }
    if(do_qdb_or_nonbanded_final_scan) { /* careful, different from just an 'else', b/c we may have just set this as true if status == eslERANGE */
      /*******************************************************************
       * Run non-HMM banded (probably qdb) version of CYK or Inside *
       *******************************************************************/
      if(cm->smx == NULL) { printf("FIX ME! ADD TRUNC VERSION! cm->smx is NULL, probably overflow sized hb mx\n"); continue; }
      pli->acct[pli->cur_pass_idx].n_overflow_final++;
      if(cm->search_opts & CM_SEARCH_INSIDE) { /* final algorithm is Inside */
	if((status = FastIInsideScan(cm, errbuf, cm->smx, SMX_QDB2_LOOSE, sq->dsq, es[i], ee[i],
				     pli->T,            /* minimum score to report */
				     hitlist,           /* our hitlist */
				     pli->do_null3,     /* apply the null3 correction? */
				     0., NULL, NULL,    /* envelope redefinition parameters, irrelevant here */
				     NULL,              /* ret_vsc, irrelevant here */
				     &finalsc)) != eslOK) { /* best score, irrelevant here */
	  printf("ERROR: %s\n", errbuf); return status; }
      }
      else { /* final algorithm is CYK */
	if((status = FastCYKScan(cm, errbuf, cm->smx, SMX_QDB2_LOOSE, sq->dsq, es[i], ee[i],
				 pli->T,            /* minimum score to report */
				 hitlist,           /* our hitlist */
				 pli->do_null3,     /* apply the null3 correction? */
				 0., NULL, NULL,    /* envelope redefinition parameters, irrelevant here */
				 NULL,              /* ret_vsc, irrelevant here */
				 &finalsc)) != eslOK) { /* best score, irrelevant here */
	  printf("ERROR: %s\n", errbuf); return status; }
      }
    }

    /* save a copy of the bands we calculated for the final search stage */
    if(pli->do_alignments && (! do_qdb_or_nonbanded_final_scan)) { 
      scan_cp9b = cp9_CloneBands(cm->cp9b, pli->errbuf);
      if(scan_cp9b == NULL) return eslEMEM;
#if eslDEBUGLEVEL >= 1
      if((status = cp9_ValidateBands(cm, pli->errbuf, cm->cp9b, es[i], ee[i])) != eslOK) return status;
      ESL_DPRINTF1(("original bands validated.\n"));
      if((status = cp9_ValidateBands(cm, pli->errbuf, scan_cp9b, es[i], ee[i])) != eslOK) return status;
      ESL_DPRINTF1(("cloned bands validated.\n"));
#endif
    }
    else { 
      scan_cp9b = NULL;
    }
  
    /* add info to each hit DP scanning functions didn't have access to, and align the hits if nec */
    for (h = nhit; h < hitlist->N; h++) { 
      hit = &(hitlist->unsrt[h]);
      hit->cm_idx   = pli->cur_cm_idx;
      hit->seq_idx  = pli->cur_seq_idx;
      hit->pass_idx = pli->cur_pass_idx;
      hit->maxW     = pli->maxW;
      hit->pvalue   = esl_exp_surv(hit->score, cm->expA[pli->final_cm_exp_mode]->mu_extrap, cm->expA[pli->final_cm_exp_mode]->lambda);
      hit->srcL     = sq->L; /* this may be -1, in which case it will be updated by caller (cmsearch or cmscan) when full length is known */

      /* initialize remaining values we don't know yet */
      hit->evalue  = 0.;
      hit->ad      = NULL;
	  
      if (pli->mode == CM_SEARCH_SEQS) { 
	if (                       (status  = esl_strdup(sq->name, -1, &(hit->name)))  != eslOK) ESL_FAIL(eslEMEM, pli->errbuf, "allocation failure");
	if (sq->acc[0]  != '\0' && (status  = esl_strdup(sq->acc,  -1, &(hit->acc)))   != eslOK) ESL_FAIL(eslEMEM, pli->errbuf, "allocation failure");
        if (sq->desc[0] != '\0' && (status  = esl_strdup(sq->desc, -1, &(hit->desc)))  != eslOK) ESL_FAIL(eslEMEM, pli->errbuf, "allocation failure");
      } 
      else {
	if ((status  = esl_strdup(cm->name, -1, &(hit->name)))  != eslOK) ESL_FAIL(eslEMEM, pli->errbuf, "allocation failure");
	if ((status  = esl_strdup(cm->acc,  -1, &(hit->acc)))   != eslOK) ESL_FAIL(eslEMEM, pli->errbuf, "allocation failure");
        if ((status  = esl_strdup(cm->desc, -1, &(hit->desc)))  != eslOK) ESL_FAIL(eslEMEM, pli->errbuf, "allocation failure");
      }
#if DOPRINT
      printf("SURVIVOR envelope     [%10ld..%10ld] survived Inside    %6.2f bits  P %g\n", hit->start, hit->stop, hit->score, hit->pvalue);
#endif
      /* Get an alignment of the hit, if nec */
      if(pli->do_alignments) { 
	if((status = cm_pli_AlignHit(pli, cm, *opt_cmcons, sq, do_trunc, hit, 
				     (h == nhit) ? TRUE : FALSE,    /* TRUE if this is the first hit we're aligning (h == nhit) */
				     scan_cp9b))                    /* a copy of the HMM bands determined in the last search stage, NULL if HMM bands not used */
	   != eslOK) return status;
      }

      /* Finally, if we're using model-specific bit score thresholds,
       * determine if the significance of the hit (is it reported
       * and/or included?)  Adapted from Sean's comments at an
       * analogous point in p7_pipeline():
       *
       * If we're using model-specific bit score thresholds (GA | TC | NC)
       * and we're in a cmscan pipeline (mode = CM_SCAN_MODELS), then we
       * *must* apply those reporting or inclusion thresholds now, because
       * this model is about to go away; we won't have its thresholds
       * after all targets have been processed.
       * 
       * If we're using E-value thresholds and we don't know the
       * search space size (Z_setby == CM_ZSETBY_NTARGETS), we 
       * *cannot* apply those thresholds now, and we *must* wait 
       * until all targets have been processed (see cm_tophits_Threshold()).
       * 
       * For any other thresholding, it doesn't matter whether we do
       * it here (model-specifically) or at the end (in
       * cm_tophits_Threshold()). 
       * 
       * What we actually do, then, is to set the flags if we're using
       * model-specific score thresholds (regardless of whether we're
       * in a scan or a search pipeline); otherwise we leave it to 
       * cm_tophits_Threshold(). cm_tophits_Threshold() is always
       * responsible for *counting* the reported, included sequences.
       * 
       * [xref J5/92]
       */

      if (pli->use_bit_cutoffs) { 
	if (cm_pli_TargetReportable(pli, hit->score, hit->evalue)) { /* evalue is invalid, but irrelevant if pli->use_bit_cutoffs */
	  hit->flags |= CM_HIT_IS_REPORTED;
	  if (cm_pli_TargetIncludable(pli, hit->score, hit->evalue)) /* ditto */
	    hit->flags |= CM_HIT_IS_INCLUDED;
	}
      }
    } /* end of 'for(h = nhit'... */
    if(scan_cp9b != NULL) { 
      FreeCP9Bands(scan_cp9b);
      scan_cp9b = NULL;
    }

    nhit = hitlist->N;
  } /* end of for each envelope loop */
  cm->tau = save_tau;
  /* free the scan matrices if we just allocated them */
  if(scan_cp9b != NULL) FreeCP9Bands(scan_cp9b);

  return eslOK;
}

/* Function:  cm_pli_AlignHit()
 * Synopsis:  Align a hit that survives all stages of the pipeline to a CM.
 * Incept:    EPN, Mon Aug  8 10:46:21 2011
 *
 * Purpose:   For a given hit <hit> in sequence <sq> spanning
 *            <hit->start> to <hit->stop>, align it to a CM and create a
 *            CM_ALIDISPLAY object and store it in <hit->ad>. The
 *            alignment algorithm used is HMM banded optimal accuracy
 *            unless that requires too much memory. See comments
 *            within code for more detail.
 * 
 *            If <do_trunc> is TRUE, we use truncated alignment
 *            algorithms, else we use standard ones.
 *
 * Returns: <eslOK> on success.
 */
int
cm_pli_AlignHit(CM_PIPELINE *pli, CM_t *cm, CMConsensus_t *cmcons, const ESL_SQ *sq, int do_trunc, CM_HIT *hit, int first_hit, CP9Bands_t *scan_cp9b)
{
  int                 status;               /* Easel status code */
  Parsetree_t        *tr = NULL;            /* pointer to the pointer to the parsetree we're currently creating */
  char               *ppstr = NULL;         /* posterior decode array of strings */
  ESL_DSQ            *subdsq;               /* ptr to start of a hit */
  int                 hitlen;               /* hit length */
  float               hbmx_Mb;              /* approximate size in Mb for HMM banded matrix for current hit */
  float               shmx_Mb;              /* approximate size in Mb for HMM banded shadow matrix for current hit */
  float               total_Mb;             /* approximate size in Mb needed for alignment of a hit */
  int                 do_optacc;            /* TRUE to do optimal accuracy alignment */
  int                 do_ppstr;             /* TRUE to derive posteriors for hit alignments */
  int                 do_hbanded;           /* TRUE to use HMM bands for alignment */
  int                 do_nonbanded;         /* TRUE to align without bands (instead of using HMM bands) */
  ESL_STOPWATCH      *watch = NULL;         /* stopwatch for timing alignment step */
  float               optacc_sc, ins_sc, cyk_sc;  /* optimal accuracy score, inside score, CYK score */
  float               null3_correction;     /* null 3 bit score penalty, for CYK score */

  watch = esl_stopwatch_Create();
  if(! watch) ESL_FAIL(eslEMEM, pli->errbuf, "out of memory");

  esl_stopwatch_Start(watch);  
	
  /* set defaults, these may be changed below */
  optacc_sc    = 0.;
  do_optacc    = FALSE;
  do_ppstr     = FALSE;
  do_hbanded   = FALSE;
  do_nonbanded = FALSE;

  hitlen = hit->stop - hit->start + 1;
  subdsq = sq->dsq + hit->start - 1;

  if(pli->align_hbanded) { 
    /* Align with HMM bands, if we can do it in the allowed amount of memory */

    /* Either recalculate the HMM bands or adjust those that we already have */
    if(scan_cp9b == NULL || pli->do_hb_recalc) { 
      /* Calculate HMM bands. We'll tighten tau and recalculate bands until 
       * the resulting HMM banded matrix is under our size limit.
       */
      cm->tau = pli->final_tau;
      while(1) { 
	if((status = cp9_Seq2Bands(cm, pli->errbuf, cm->cp9_mx, cm->cp9_bmx, cm->cp9_bmx, subdsq, 1, hitlen, cm->cp9b, FALSE, pli->tro, 0)) != eslOK) return status; 
	if(do_trunc) { 
	  if((status = cm_tr_hb_mx_SizeNeeded(cm, pli->errbuf, cm->cp9b, hitlen, NULL, NULL, NULL, NULL, &hbmx_Mb)) != eslOK) return status; 
	}
	else { 
	  if((status = cm_hb_mx_SizeNeeded   (cm, pli->errbuf, cm->cp9b, hitlen, NULL, &hbmx_Mb)) != eslOK) return status; 
	}
	if(hbmx_Mb < pli->hb_size_limit) break; /* our matrix will be small enough, break out of while(1) */
	if(cm->tau > 0.01)               break; /* our bands have gotten too tight, break out of while(1) */
	/* printf("ARC 0 tau: %10g  hbmx_Mb: %10.2f\n", cm->tau, hbmx_Mb); */
	cm->tau *= pli->xtau;
      }	  
      /* printf("ARC 1 tau: %10g  hbmx_Mb: %10.2f\n", cm->tau, hbmx_Mb); */
    }
    else { 
      /* shift existing CP9 HMM bands by a fixed offset, this guarantees our alignment will be the same hit our search found */
      if(! first_hit) { 
	/* This is not the first hit in this envelope, we need to clone the scan_cp9b bands that were passed in, 
	 * before shifting them (because they were shifted for the previous hit in this envelope and are no longer usable) 
	 */
	if(cm->cp9b != NULL) FreeCP9Bands(cm->cp9b);
	cm->cp9b = cp9_CloneBands(scan_cp9b, pli->errbuf);
	if(cm->cp9b == NULL) return status;
      }
      cp9_ShiftCMBands(cm, hit->start, hit->stop, do_trunc);
      /* NOTE: cm->cp9b bands would currently fail a cp9_ValidateBands() check..., but they'll work for our purposes here */
    }
    
    /* Determine the number of cells needed in the CM_HB_MX's required
     * for alignment. It it's above our limit, we won't use HMM bands.
     *
     *           matrix sizes                 pli->align_cyk   algorithm to use  HMM banded? posteriors?
     *          ----------------------------  -------------  ----------------  ----------- -----------
     * if       hbmx_Mb < pli->hb_size_limit  FALSE          optimal accuracy  yes         yes
     * else if  hbmx_Mb < pli->hb_size_limit  TRUE           CYK               yes         yes
     * else if  hbmx_Mb > pli->hb_size_limit  TRUE/FALSE     D&C CYK           no          no
     */     
    if(do_trunc) { 
      if((status = cm_tr_hb_mx_SizeNeeded       (cm, pli->errbuf, cm->cp9b, hitlen, NULL, NULL, NULL, NULL, &hbmx_Mb)) != eslOK) return status; 
      if((status = cm_tr_hb_shadow_mx_SizeNeeded(cm, pli->errbuf, cm->cp9b, NULL,   NULL, NULL, NULL, NULL, NULL, NULL, &shmx_Mb)) != eslOK) return status; 
    }
    else { 
      if((status = cm_hb_mx_SizeNeeded          (cm, pli->errbuf, cm->cp9b, hitlen, NULL, &hbmx_Mb)) != eslOK) return status; 
      if((status = cm_hb_shadow_mx_SizeNeeded   (cm, pli->errbuf, cm->cp9b, NULL,   NULL, &shmx_Mb)) != eslOK) return status; 
    }
    /* printf("ALN 1 tau: %10g  hbmx_Mb: %10.2f\n", cm->tau, hbmx_Mb); */
    
    if(hbmx_Mb < pli->hb_size_limit) { 
      pli->acct[pli->cur_pass_idx].n_aln_hb++;
      do_optacc    = (pli->align_cyk) ? FALSE : TRUE;
      total_Mb     = (pli->align_cyk) ? (hbmx_Mb + shmx_Mb) : (2*hbmx_Mb + shmx_Mb);
      do_ppstr     = TRUE;
      do_hbanded   = TRUE;
      do_nonbanded = FALSE;
    }
    else { /* not enough memory for HMM banded alignment, fall back to D&C CYK */
      pli->acct[pli->cur_pass_idx].n_aln_dccyk++;
      do_optacc    = FALSE;
      do_ppstr     = FALSE;
      ppstr        = NULL;
      do_hbanded   = FALSE;
      do_nonbanded = TRUE;
    }
  } /* end of if(pli->do_hbanded) */
  else { 
    do_nonbanded = TRUE;
  }
  
  if(do_hbanded) { 
    if(do_trunc) { 
      status = cm_TrAlignHB(cm, pli->errbuf, subdsq, hitlen, 
			    pli->hb_size_limit,                            /* limit for a single CM_HB_MX, so this is safe */
			    hit->mode,                                     /* alignment mode, determined in truncated scanner function */
			    do_optacc,                                     /* use optimal accuracy alg? */
			    FALSE,                                         /* don't sample aln from Inside matrix */
			    cm->trhbmx, cm->trshhbmx,                      /* inside, shadow matrices */
			    (do_ppstr || do_optacc) ? cm->trohbmx : NULL,  /* outside DP matrix */
			    (do_ppstr || do_optacc) ? cm->trehbmx : NULL,  /* emit matrix */
			    NULL,                                          /* ESL_RANDOMNESS, unneeded b/c we're not sampling */
			    (do_ppstr  ? &ppstr  : NULL),                  /* posterior string */
			    (do_optacc ? &ins_sc : NULL),                  /* inside score, NULL if we're not doing opt acc */
			    &tr,                                           /* parsetree */
			    NULL,                                          /* mode of optimal alignment, irrelevant */
			    (do_optacc ? &optacc_sc : &cyk_sc));           /* optimal accuracy or CYK score */
      /* TEMP BLOCK */
      ParsetreeDump(stdout, tr, cm, subdsq, NULL, NULL);
      float parsetree_sc;
      ParsetreeScore(cm, NULL, NULL, tr, subdsq, FALSE, &parsetree_sc, NULL, NULL, NULL, NULL);
      printf("Parsetree score      : %.4f           (TROPTACC OR TRCYK)\n", parsetree_sc);
      /* END TEMP BLOCK */
    }
    else { 
      status = cm_AlignHB(cm, pli->errbuf, subdsq, hitlen, 
			  pli->hb_size_limit,                          /* limit for a single CM_HB_MX, so this is safe */
			  do_optacc,                                   /* use optimal accuracy alg? */
			  FALSE,                                       /* don't sample aln from Inside matrix */
			  cm->hbmx, cm->shhbmx,                        /* inside, shadow matrices */
			  (do_ppstr || do_optacc) ? cm->ohbmx : NULL,  /* outside DP matrix */
			  (do_ppstr || do_optacc) ? cm->ehbmx : NULL,  /* emit matrix */
			  NULL,                                        /* ESL_RANDOMNESS, unneeded b/c we're not sampling */
			  (do_ppstr  ? &ppstr  : NULL),                /* posterior string */
			  (do_optacc ? &ins_sc : NULL),                /* inside score, NULL if we're not doing opt acc */
			  &tr,                                         /* parsetree */
			  (do_optacc ? &optacc_sc : &cyk_sc));         /* optimal accuracy or CYK score */
    }
    if(status == eslERANGE) { 
      /* matrix was too big, despite our pre-check (should be rare), use D&C CYK */
      do_optacc    = FALSE;
      do_ppstr     = FALSE;
      do_hbanded   = FALSE; 
      do_nonbanded = TRUE;
    }
    else if(status != eslOK) return status;
    esl_stopwatch_Stop(watch); /* we started it above before we calc'ed the CP9 bands */ 
  }
  if((! do_hbanded) && do_nonbanded) { /* use non-banded D&C CYK (slow!) */
    esl_stopwatch_Start(watch);  
    cyk_sc = CYKDivideAndConquer(cm, subdsq, hitlen, 0, 1, hitlen, &tr, NULL, NULL); 
    esl_stopwatch_Stop(watch);  
    total_Mb = CYKNonQDBSmallMbNeeded(cm, hitlen);
    ppstr = NULL;
  }
    
  /* add null3 correction to cyk_sc if necessary (optacc_sc doesn't get null3-corrected, it's an avg PP) */
  if((! do_optacc) && pli->do_null3) { 
    ScoreCorrectionNull3CompUnknown(cm->abc, cm->null, subdsq, 1, hitlen, cm->null3_omega, &null3_correction);
    cyk_sc -= null3_correction;
  }

  hit->ad = cm_alidisplay_Create(cm->abc, tr, cm, cmcons, sq, hit->start, ppstr, 
				 (do_optacc) ? optacc_sc : cyk_sc, 
				 do_optacc, do_hbanded, total_Mb, watch->elapsed);
  /*cm_alidisplay_Dump(stdout, hit->ad);*/

  /* clean up and return */
  FreeParsetree(tr);
  if(do_ppstr && ppstr != NULL) free(ppstr);
  if(watch  != NULL) { esl_stopwatch_Destroy(watch); }

  return eslOK;
}

/* Function:  merge_windows_from_two_lists()
 * Synopsis:  Merge two sorted lists of windows into one, collapsing overlapping windows together. 
 * Incept:    EPN, Mon Nov 29 14:47:40 2010
 *
 * Purpose:   Given two lists of windows, merge them into one list
 *            by collapsing any overlapping windows into a single window. 
 *
 *            Each window in list 1 is defined as:
 *            <ws1[i]>..we1[i] for i=0..<nwin1-1>
 *
 *            Each window in list 2 is defined as:
 *            <ws2[i]>..we2[i] for i=0..<nwin2-1>
 *
 *            The windows within each list must be sorted in
 *            increasing order and be non-overlapping, i.e. the
 *            following must hold:
 *
 *            ws1[i] <  ws1[j] for all i < j.
 *            ws1[i] <= we1[i] for all i.
 *
 *            ws2[i] <  ws2[j] for all i < j.
 *            ws2[i] <= we2[i] for all i.
 *
 *            Note, we check that these conditions hold at the
 *            beginning of the function, and return eslEINVAL if not.
 *
 *            A new list is created and returned, as <mws>,
 *            <mwe>, and <nmwin>, start, end positions and 
 *            number of windows respectively. 
 *
 *            The lowest P-value of any of the merged windows is also
 *            kept, and returned in <ret_mwp>, <ret_mwp[i]> is the
 *            lowest P-value of any of the P-values from <wp1>/<wp2>
 *            for windows that were merged together to create
 *            <mws>..<mwe>. If the P-value was from a window in list 1
 *            (<wp1>), <ret_mwl[i]> is set as 1, if it is from a
 *            window from list 2 (<wp2>), <ret_mwl[i]> is set as 2.
 * 
 * Returns:   New merged list in <ret_mws>, <ret_mwe>, <ret_mwp>,
 *            <ret_mwl>, <ret_nmwin>, storing start position, end
 *            position, P-value and list origin of each of the
 *            <ret_nmwin> windows in the merged list.
 */
int
merge_windows_from_two_lists(int64_t *ws1, int64_t *we1, double *wp1, int *wl1, int nwin1, int64_t *ws2, int64_t *we2, double *wp2, int *wl2, int nwin2, int64_t **ret_mws, int64_t **ret_mwe, double **ret_mwp, int **ret_mwl, int *ret_nmwin)
{
  int status;
  int64_t *mws = NULL;
  int64_t *mwe = NULL;
  double  *mwp = NULL;
  int     *mwl = NULL;
  int nmwin;
  int mi, i1, i2;
  int nalloc;

  /* check that our required conditions hold */
  for(i1 = 0; i1 < nwin1-1; i1++) if(ws1[i1] >= ws1[i1+1]) return eslEINVAL;
  for(i2 = 0; i2 < nwin2-1; i2++) if(ws2[i2] >= ws2[i2+1]) return eslEINVAL;
  for(i1 = 0; i1 < nwin1; i1++) if(ws1[i1] > we1[i1]) return eslEINVAL;
  for(i2 = 0; i2 < nwin2; i2++) if(ws2[i2] > we2[i2]) return eslEINVAL;

  nalloc = nwin1 + nwin2; /* we'll never exceed this */
  ESL_ALLOC(mws, sizeof(int64_t) * nalloc);
  ESL_ALLOC(mwe, sizeof(int64_t) * nalloc);
  ESL_ALLOC(mwp, sizeof(double)  * nalloc);
  ESL_ALLOC(mwl, sizeof(int)     * nalloc);

  i1 = 0;
  i2 = 0;
  mi = 0;
  while(i1 < nwin1 || i2 < nwin2) { 
    if((i1 < nwin1) && ((i2 == nwin2) || (ws1[i1] < ws2[i2]))) {
      /**************************************************
       * case 1: our next hit, in order is from list 1 *
       **************************************************/
      /* initialize merged hit as copy of hit from list 1 */
      mws[mi] = ws1[i1]; /* our merged window begins at ws1[i1], this won't change */
      mwp[mi] = wp1[i1]; /* for now, our best P-value is wp1[i1] */
      mwl[mi] = wl1[i1]; /* for now, our best P-value comes from list 1 */
      mwe[mi] = we1[i1]; /* for now, our merged window ends at we1[i1] */
      /* merge with all overlapping windows in list 2 */
      while((i2 < nwin2) && (mwe[mi] >= ws2[i2])) { /* check all windows in list 2 that overlap with our current list 1 window */
	mwe[mi] = ESL_MAX(mwe[mi], we2[i2]); /* our merged window now ends at larger of mwe[me] and we2[i2] */
	if(wp2[i2] < mwp[mi]) { 
	  mwp[mi] = wp2[i2];  /* our best P-value now is wp2[i2] */
	  mwl[mi] = wl2[i2];  /* our best P-value now comes from list 2 */
	}
	i2++;
      }
      i1++;
      /* finally, if we merged any windows from list 2, our window is possibly larger, 
       * so we merge with all now-overlapping windows in list 1 */
      while((i1 < nwin1) && (mwe[mi] >= ws1[i1])) { 
	mwe[mi] = ESL_MAX(mwe[mi], we1[i1]); /* our merged window now ends at larger of mwe[me] and we2[i2] */
	if(wp1[i1] < mwp[mi]) { 
	  mwp[mi] = wp1[i1];  /* our best P-value now is wp1[i1] */
	  mwl[mi] = wl1[i1];  /* our best P-value now comes from list 1 */
	}
	i1++;
      }
      mi++;
    }
    else { 
      /****************************************************************************************************************
       * case 2: our next hit, in order is from list 2, same code as case 1, with list 1 and list 2 variables inverted *
       ****************************************************************************************************************/
      /* initialize merged hit as copy of hit from list 2 */
      mws[mi] = ws2[i2]; /* our merged window begins at ws1[i2], this won't change */
      mwp[mi] = wp2[i2]; /* for now, our best P-value is wp2[i2] */
      mwl[mi] = wl2[i2]; /* for now, our best P-value comes from list 2 */
      mwe[mi] = we2[i2]; /* for now, our merged window ends at we2[i2] */
      /* merge with all overlapping windows in list 1 */
      while((i1 < nwin1) && (mwe[mi] >= ws1[i1])) { /* check all windows in list 1 that overlap with our current list 2 window */
	mwe[mi] = ESL_MAX(mwe[mi], we1[i1]); /* our merged window now ends at larger of mwe[me] and we1[i1] */
	if(wp1[i1] < mwp[mi]) { 
	  mwp[mi] = wp1[i1];  /* our best P-value now is wp1[i1] */
	  mwl[mi] = wl1[i1];  /* our best P-value now comes from list 1 */
	}
	i1++;
      }
      i2++;
      /* finally, if we merged any windows from list 1, our window is possibly larger, 
       * so we merge with all now-overlapping windows in list 2 */
      while((i2 < nwin2) && (mwe[mi] >= ws2[i2])) { /* check all windows in list 2 that overlap with our current list 1 window */
	mwe[mi] = ESL_MAX(mwe[mi], we2[i2]); /* our merged window now ends at larger of mwe[me] and we2[i2] */
	if(wp2[i2] < mwp[mi]) { 
	  mwp[mi] = wp2[i2];  /* our best P-value now is wp2[i2] */
	  mwl[mi] = wl2[i2];  /* our best P-value now comes from list 2 */
	}
	i2++;
      }
      mi++;
    }
  }
  nmwin = mi;

  /*
  for(i1 = 0; i1 < nwin1; i1++) printf("list1 win %5d  %10" PRId64 "..%10" PRId64 " P: %10g\n", i1, ws1[i1], we1[i1], wp1[i1]);
  printf("\n");
  for(i2 = 0; i2 < nwin2; i2++) printf("list2 win %5d  %10" PRId64 "..%10" PRId64 " P: %10g\n", i2, ws2[i2], we2[i2], wp2[i2]);
  printf("\n");
  for(mi = 0; mi < nmwin; mi++) printf("mlist win %5d  %10" PRId64 "..%10" PRId64 " P: %10g  %d\n", mi, mws[mi], mwe[mi], mwp[mi], mwl[mi]);
  printf("\n");
  */

  *ret_mws = mws;
  *ret_mwe = mwe;
  *ret_mwp = mwp;
  *ret_mwl = mwl;
  *ret_nmwin = nmwin;

  return eslOK;

 ERROR:
  return status;
}

/* Function:  cm_Pipeline()
 * Synopsis:  The accelerated seq/profile comparison pipeline using HMMER3 scanning.
 * Incept:    EPN, Fri Sep 24 16:42:21 2010
 *            TJW, Fri Feb 26 10:17:53 2018 [Janelia] (p7_Pipeline_Longtargets())
 *
 * Purpose:   Run the accelerated pipeline to compare profile <om>
 *            against sequence <sq>. This function calls three other
 *            functions: cm_pli_p7Filter(), cm_pli_p7EnvelopeDef, and
 *            cm_pli_CMStage().  This organization allows us to
 *            use multiple p7 HMMs as filters.
 *
 *            <research_5term> and <research_3term> indicate whether we
 *            are re-searching a terminal chunk of the sequence either
 *            5' (in which case we have the first residue) and/or 3'
 *            in which case we have the final residue. 
 *
 * Returns:   <eslOK> on success. If a significant hit is obtained,
 *            its information is added to the growing <hitlist>.
 *
 *            <eslEINVAL> if (in a scan pipeline) we're supposed to
 *            set GA/TC/NC bit score thresholds but the model doesn't
 *            have any.
 *
 *            <eslERANGE> on numerical overflow errors in the
 *            optimized vector implementations; particularly in
 *            posterior decoding. I don't believe this is possible for
 *            multihit local models, but I'm set up to catch it
 *            anyway. We may emit a warning to the user, but cleanly
 *            skip the problematic sequence and continue.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 * Xref:      J4/25.
 */
int
cm_Pipeline(CM_PIPELINE *pli, off_t cm_offset, P7_OPROFILE *om, P7_BG *bg, float *p7_evparam, FM_HMMDATA *fm_hmmdata, ESL_SQ *sq, CM_TOPHITS *hitlist, P7_PROFILE **opt_gm, P7_PROFILE **opt_Rgm, P7_PROFILE **opt_Lgm, P7_PROFILE **opt_Tgm, CM_t **opt_cm, CMConsensus_t **opt_cmcons)
{
  int       status;
  int       nwin = 0;   /* number of windows surviving MSV & Vit & lFwd, filled by cm_pli_p7Filter() */
  int64_t  *ws = NULL;  /* [0..i..nwin-1] window start positions, filled by cm_pli_p7Filter() */
  int64_t  *we = NULL;  /* [0..i..nwin-1] window end   positions, filled by cm_pli_p7Filter() */
  int       nenv = 0;   /* number of envelopes surviving MSV & Vit & lFwd & gFwd & EnvDef, filled by cm_pli_p7EnvelopeDef */
  int64_t  *es  = NULL; /* [0..i..nenv-1] window start positions, filled by cm_pli_p7EnvelopeDef() */
  int64_t  *ee  = NULL; /* [0..i..nenv-1] window end   positions, filled by cm_pli_p7EnvelopeDef() */

  /* variables necessary for re-searching sequence ends */
  ESL_SQ   *sq2search = NULL;  /* a pointer to the sequence to search on current pass */
  ESL_SQ   *term5sq   = NULL;  /* a copy of the 5'-most pli->maxW residues from sq */
  ESL_SQ   *term3sq   = NULL;  /* a copy of the 3'-most pli->maxW residues from sq */
  int       p;                 /* counter over passes */
  int       nwin_pass1 = 0;    /* number of windows that survived pass 1 */
  int       do_pass2;          /* should we do pass 2? re-search the 5'-most pli->maxW residues */
  int       do_pass3;          /* should we do pass 3? re-search the 3'-most pli->maxW residues */
  int       do_pass4;          /* should we do pass 4? re-search the full sequence allowing for T truncated hits */
  int       not5term;          /* TRUE if sq definitely does not contain the 5'-most pli->maxW residues of its source sequence, we can skip passes 2 and 4 */
  int       not3term;          /* TRUE if sq definitely does not contain the 3'-most pli->maxW residues of its source sequence, we can skip passes 3 and 4 */
  int       h;                 /* counter over hits */
  int       prv_ntophits;      /* number of hits */
  int64_t   start_offset;      /* offset to add to start/stop coordinates of hits found in pass 3, in which we re-search the 3' terminus */

  if (sq->n == 0) return eslOK;    /* silently skip length 0 seqs; they'd cause us all sorts of weird problems */

#if 0
  printf("F1: %f\n", pli->F1);
  printf("F2: %f\n", pli->F2);
  printf("F3: %f\n", pli->F3);
  printf("F4: %f\n", pli->F4);
  printf("F5: %f\n", pli->F5);
  printf("F6: %f\n", pli->F6);

  printf("F1b: %f\n", pli->F1b);
  printf("F2b: %f\n", pli->F2b);
  printf("F3b: %f\n", pli->F3b);
  printf("F4b: %f\n", pli->F4b);
  printf("F5b: %f\n", pli->F5b);
#endif

  not5term = not3term = FALSE; /* by default, we may have 5' terminus and we may have 3' terminus */
  /* Determine, if we can, that we definitely do not have the 5' or 3' terminii.
   * The two if statements involving sq->L are commented out. This has the effect
   * that sq->L is always treated as unknown and so we always re-search the 3' terminus.
   * We do this because we want mpi versus serial/threaded cmsearch to give identical
   * pipeline statistics, and while in mpi we know sq->L (b/c we've read it from an 
   * SSI index), in serial/threaded we don't. So to get identical results we ignore
   * the knowledge of sq->L in cmsearch mpi. In cmscan, we always read complete sequences
   * so we always know L, but we will also always research the 3' terminus, thus 
   * commenting out these two lines does not affect cmscan.
   */
  if(sq->start <= sq->end) { /* not in reverse complement (or 1 residue sequence, in revcomp) */
    if(sq->start != 1)                  not5term = TRUE;
    /*if(sq->L != -1 && sq->end != sq->L) not3term = TRUE;*/
  }
  else { /* in reverse complement */
    if(sq->end != 1)                       not3term = TRUE;
    /*if(sq->L != -1 && sq->start != sq->L)  not5term = TRUE;*/
  }

#if DOPRINT
  /*printf("\nPIPELINE ENTRANCE %s  %s  %" PRId64 " residues (pli->maxW: %d om->max_length: %d cm->W: %d)\n", sq->name, sq->desc, sq->n, pli->maxW, om->max_length, (*opt_cm)->W);*/
  printf("\nPIPELINE ENTRANCE %-15s  (n: %6" PRId64 " start: %6" PRId64 " end: %6" PRId64 " C: %6" PRId64 " W: %6" PRId64 " L: %6" PRId64 " not5term: %d not3term: %d have5term: %5s have3term: %5s)\n",
	 sq->name, sq->n, sq->start, sq->end, sq->C, sq->W, sq->L, not5term, not3term, not5term ? "No" : "Maybe", not3term ? "No" : "Maybe");
#endif

  do_pass2 = ((! not5term) && pli->research_ends) ? TRUE : FALSE;
  do_pass3 = ((! not3term) && pli->research_ends) ? TRUE : FALSE;
  do_pass4 = ((! not5term) && (! not3term) && pli->research_ends && sq->n <= pli->maxW) ? TRUE : FALSE;

  for(p = 1; p <= 4; p++) { 
    if(p == 2 && (! do_pass2))    continue;
    if(p == 3 && (! do_pass3))    continue;
    if(p == 4 && (! do_pass4))    continue;

    /* Update nres for non-standard pass (this is done by caller via cm_pli_NewSeq() for pass 1)
     * It's important to do this precisely here, in between the 3 'continue's above, but 
     * prior the one below that 'continue's only because we know no windows pass local Fwd (F3) 
     */
    if(p > 1) pli->acct[p].nres += ESL_MIN(pli->maxW, sq->n);

    /* if we know there's no windows that pass local Fwd (F3), our terminal seqs won't have any either, continue */
    if(p > 1 && nwin_pass1 == 0) continue; 

    /* set sq2search and remember start_offset */
    start_offset = 0;
    if(p == 1 || p == 4 || sq->n <= pli->maxW) { 
      sq2search = sq;
    }
    else if(p == 2) { /* research first (5') pli->maxW residues */
      term5sq = esl_sq_CreateDigital(bg->abc);
      copy_subseq(sq, term5sq, 1, pli->maxW);
      sq2search = term5sq;
    }
    else { /* p must be 3 */
      term3sq = esl_sq_CreateDigital(bg->abc);
      copy_subseq(sq, term3sq, sq->n - pli->maxW + 1, pli->maxW);
      sq2search = term3sq;
      start_offset = sq->n - pli->maxW;
    }

    /* Set rs5term and rs3term and tro, this informs the pipeline stages about truncated hit/alignment options */
    pli->rs5term          = (p == 2 || p == 4) ? TRUE : FALSE;
    pli->rs3term          = (p == 3 || p == 4) ? TRUE : FALSE;
    pli->tro->allowR      = (pli->do_trunc_ends && pli->rs5term)                       ? TRUE : FALSE;
    pli->tro->allowL      = (pli->do_trunc_ends && pli->rs3term)                       ? TRUE : FALSE;
    pli->tro->force_i0_RT = (pli->do_trunc_ends && pli->rs5term && pli->do_force_ends) ? TRUE : FALSE;
    pli->tro->force_j0_LT = (pli->do_trunc_ends && pli->rs3term && pli->do_force_ends) ? TRUE : FALSE;

    pli->cur_pass_idx = p;

    /* execute the pipeline */
#if DOPRINT
    printf("\nPIPELINE calling p7Filter() %s  %" PRId64 " residues (pass: %d rs5term: %d rs3term: %d allowR: %d allowL: %d force_i0_RT: %d force_j0_LT: %d)\n", sq2search->name, sq2search->n, p, pli->rs5term, pli->rs3term, pli->tro->allowR, pli->tro->allowL, pli->tro->force_i0_RT, pli->tro->force_j0_LT);
#endif
    if((status = cm_pli_p7Filter(pli, om, bg, p7_evparam, fm_hmmdata, sq2search, &ws, &we, &nwin)) != eslOK) return status;
    if(p == 1) nwin_pass1 = nwin;
    if(pli->do_time_F1 || pli->do_time_F2 || pli->do_time_F3) return status;

#if DOPRINT
    printf("\nPIPELINE calling p7EnvelopeDef() %s  %" PRId64 " residues (pass: %d rs5term: %d rs3term: %d allowR: %d allowL: %d force_i0_RT: %d force_j0_LT: %d)\n", sq2search->name, sq2search->n, p, pli->rs5term, pli->rs3term, pli->tro->allowR, pli->tro->allowL, pli->tro->force_i0_RT, pli->tro->force_j0_LT);
#endif
    if((status = cm_pli_p7EnvelopeDef(pli, om, bg, p7_evparam, sq2search, ws, we, nwin, opt_gm, opt_Rgm, opt_Lgm, opt_Tgm, &es, &ee, &nenv)) != eslOK) return status;

#if DOPRINT
    printf("\nPIPELINE calling CMStage() %s  %" PRId64 " residues (pass: %d rs5term: %d rs3term: %d allowR: %d allowL: %d force_i0_RT: %d force_j0_LT: %d)\n", sq2search->name, sq2search->n, p, pli->rs5term, pli->rs3term, pli->tro->allowR, pli->tro->allowL, pli->tro->force_i0_RT, pli->tro->force_j0_LT);
#endif
    prv_ntophits = hitlist->N;
    if((status = cm_pli_CMStage(pli, cm_offset, sq2search, es, ee, nenv, hitlist, opt_cm, opt_cmcons)) != eslOK) return status;
    if(pli->do_time_F6) return status;

    /* Two ways we have to update any hits we've just found: */
    if(hitlist->N > prv_ntophits) { 
      /* 1. if we're in a terminus researching pass (i.e. not pass 1), raise CM_HIT_FROM_TERMINUS_RESEARCH flag */
      if(p != 1) { 
	/*printf("setting hits %4d..%" PRId64 " flags to CM_HIT_FROM_TERMINUS_RESEARCH\n", prv_ntophits, hitlist->N);*/
	for(h = prv_ntophits; h < hitlist->N; h++) hitlist->unsrt[h].flags |= CM_HIT_FROM_TERMINUS_RESEARCH;
      }
      /* 2. if we're researching a 3' terminus, adjust the start/stop positions so they are relative to the actual 5' start */
      if(start_offset != 0) { /* this will only be non-zero if we're in pass3 */
	for(h = prv_ntophits; h < hitlist->N; h++) { 
	  hitlist->unsrt[h].start += start_offset;
	  hitlist->unsrt[h].stop  += start_offset;
	  if(hitlist->unsrt[h].ad != NULL) { 
	    hitlist->unsrt[h].ad->sqfrom += start_offset;
	    hitlist->unsrt[h].ad->sqto   += start_offset;
	  }
	}
      }
    }
  /*****************************************************************/

    if(ws  != NULL) free(ws);
    if(we  != NULL) free(we);
    if(es  != NULL) free(es);
    if(ee  != NULL) free(ee);
  }
  
  if(term5sq != NULL) esl_sq_Destroy(term5sq);
  if(term3sq != NULL) esl_sq_Destroy(term3sq);

  return eslOK;
}
  
/* Function:  cm_pli_Statistics()
 * Synopsis:  Final statistics output from a processing pipeline.
 * Incept:    EPN, Fri Sep 24 16:48:06 2010 
 *            SRE, Tue Dec  9 10:19:45 2008 [Janelia] (p7_pli_Statistics())
 *
 * Purpose:   Print a standardized report of the internal statistics of
 *            a finished processing pipeline <pli> to stream <ofp>.
 *            
 *            If stopped, non-<NULL> stopwatch <w> is provided for a
 *            stopwatch that was timing the pipeline, then the report
 *            includes timing information.
 *
 * Returns:   <eslOK> on success.
 */
int
cm_pli_Statistics(FILE *ofp, CM_PIPELINE *pli, int pass_idx, ESL_STOPWATCH *w)
{
  double  ntargets; 
  int64_t nwin_fcyk  = 0; /* number of windows CYK filter evaluated */
  int64_t nwin_final = 0; /* number of windows final stage evaluated */

  CM_PLI_ACCT *pli_acct = &(pli->acct[pass_idx]);

  fprintf(ofp, "Internal pipeline statistics summary: (%s)\n", cm_pli_DescribePass(pli, pass_idx));
  fprintf(ofp, "-------------------------------------\n");
  if (pli->mode == CM_SEARCH_SEQS) {
    fprintf(ofp,   "Query model(s):                                    %15" PRId64 "  (%" PRId64 " consensus positions)\n",     pli->nmodels, pli->nnodes);
    fprintf(ofp,   "Target sequences:                                  %15" PRId64 "  (%" PRId64 " residues)\n",  pli->nseqs,   pli_acct->nres);
    ntargets = pli->nseqs;
  } else {
    fprintf(ofp,   "Query sequence(s):                                 %15" PRId64 "  (%" PRId64 " residues)\n",  pli->nseqs,   pli_acct->nres);
    fprintf(ofp,   "Target model(s):                                   %15" PRId64 "  (%" PRId64 " consensus positions)\n",     pli->nmodels, pli->nnodes);
    ntargets = pli->nmodels;
  }

  if(pli->do_msv) { 
    fprintf(ofp, "Windows   passing  local HMM MSV           filter: %15" PRId64 "  (%.4g); expected (%.4g)\n",
	    pli_acct->n_past_msv,
	    (double) pli_acct->pos_past_msv / pli_acct->nres,
	    pli->F1 * pli->nmodels);
    nwin_fcyk = nwin_final = pli_acct->n_past_msv;
  }
  else { 
    fprintf(ofp, "Windows   passing  local HMM MSV           filter: %15s  (off)\n", "");
  }

  if(pli->do_msvbias) { 
    fprintf(ofp, "Windows   passing  local HMM MSV      bias filter: %15" PRId64 "  (%.4g); expected (%.4g)\n",
	    pli_acct->n_past_msvbias,
	    (double) pli_acct->pos_past_msvbias / pli_acct->nres,
	    pli->F1b * pli->nmodels);
    nwin_fcyk = nwin_final = pli_acct->n_past_msvbias;
  }

  if(pli->do_vit) { 
    fprintf(ofp, "Windows   passing  local HMM Viterbi       filter: %15" PRId64 "  (%.4g); expected (%.4g)\n",
	    pli_acct->n_past_vit,
	    (double) pli_acct->pos_past_vit / pli_acct->nres,
	    pli->F2 * pli->nmodels);
    nwin_fcyk = nwin_final = pli_acct->n_past_vit;
  }
  else { 
    fprintf(ofp, "Windows   passing  local HMM Viterbi       filter: %15s  (off)\n", "");
  }

  if(pli->do_vitbias) { 
    fprintf(ofp, "Windows   passing  local HMM Viterbi  bias filter: %15" PRId64 "  (%.4g); expected (%.4g)\n",
	    pli_acct->n_past_vitbias,
	    (double) pli_acct->pos_past_vitbias / pli_acct->nres,
	    pli->F2b * pli->nmodels);
    nwin_fcyk = nwin_final = pli_acct->n_past_vitbias;
  }
  else { 
    fprintf(ofp, "Windows   passing  local HMM Viterbi  bias filter: %15s  (off)\n", "");
  }

  if(pli->do_fwd) { 
    fprintf(ofp, "Windows   passing  local HMM Forward       filter: %15" PRId64 "  (%.4g); expected (%.4g)\n",
	    pli_acct->n_past_fwd,
	    (double) pli_acct->pos_past_fwd / pli_acct->nres,
	    pli->F3 * pli->nmodels);
    nwin_fcyk = nwin_final = pli_acct->n_past_fwd;
  }
  else { 
    fprintf(ofp, "Windows   passing  local HMM Forward       filter: %15s  (off)\n", "");
  }

  if(pli->do_fwdbias) { 
    fprintf(ofp, "Windows   passing  local HMM Forward  bias filter: %15" PRId64 "  (%.4g); expected (%.4g)\n",
	    pli_acct->n_past_fwdbias,
	    (double) pli_acct->pos_past_fwdbias / pli_acct->nres,
	    pli->F3b * pli->nmodels);
    nwin_fcyk = nwin_final = pli_acct->n_past_fwdbias;
  }
  else { 
    fprintf(ofp, "Windows   passing  local HMM Forward  bias filter: %15s  (off)\n", "");
  }

  if(pli->do_gfwd) { 
    fprintf(ofp, "Windows   passing glocal HMM Forward       filter: %15" PRId64 "  (%.4g); expected (%.4g)\n",
	    pli_acct->n_past_gfwd,
	    (double) pli_acct->pos_past_gfwd / pli_acct->nres,
	    pli->F4 * pli->nmodels);
    nwin_fcyk = nwin_final = pli_acct->n_past_gfwd;
  }
  else { 
    fprintf(ofp, "Windows   passing glocal HMM Forward       filter: %15s  (off)\n", "");
  }
  if(pli->do_gfwdbias) { 
    fprintf(ofp, "Windows   passing glocal HMM Forward  bias filter: %15" PRId64 "  (%.4g); expected (%.4g)\n",
	    pli_acct->n_past_gfwdbias,
	    (double) pli_acct->pos_past_gfwdbias / pli_acct->nres,
	    pli->F4b * pli->nmodels);
    nwin_fcyk = nwin_final = pli_acct->n_past_gfwdbias;
  }
  else { 
    fprintf(ofp, "Windows   passing glocal HMM Forward  bias filter: %15s  (off)\n", "");
  }

  if(pli->do_envelopes) { 
    fprintf(ofp, "Envelopes passing glocal HMM envelope defn filter: %15" PRId64 "  (%.4g); expected (%.4g)\n",
	    pli_acct->n_past_edef,
	    (double) pli_acct->pos_past_edef / pli_acct->nres,
	    pli->F5 * pli->nmodels);
    nwin_fcyk = nwin_final = pli_acct->n_past_edef;
  }
  else { 
    fprintf(ofp, "Envelopes passing glocal HMM envelope defn filter: %15s  (off)\n", "");
  }

  if(pli->do_edefbias) { 
    fprintf(ofp, "Envelopes passing glocal HMM envelope bias filter: %15" PRId64 "  (%.4g); expected (%.4g)\n",
	    pli_acct->n_past_edefbias,
	    (double) pli_acct->pos_past_edefbias / pli_acct->nres,
	    pli->F5b * pli->nmodels);
    nwin_fcyk = nwin_final = pli_acct->n_past_edefbias;
  }
  else { 
    fprintf(ofp, "Windows passing glocal HMM envelope bias filter: %15s  (off)\n", "");
  }

  if(pli->do_cyk) { 
    fprintf(ofp, "Envelopes passing %6s CM  CYK           filter: %15" PRId64 "  (%.4g); expected (%.4g)\n",
	    (pli->do_glocal_cm_stages) ? "glocal" : "local",
	    pli_acct->n_past_cyk,
	    (double) pli_acct->pos_past_cyk / pli_acct->nres,
	    pli->F6 * pli->nmodels);
    nwin_final = pli_acct->n_past_cyk;
  }
  else { 
    fprintf(ofp, "Envelopes passing %6s CM  CYK           filter: %15s  (off)\n", 
	    (pli->do_glocal_cm_stages) ? "glocal" : "local",
	    "");
  }

  fprintf(ofp, "Total hits reported:                               %15d  (%.4g)\n",
          (int)pli_acct->n_output,
          (double) pli_acct->pos_output / pli_acct->nres );

  fprintf(ofp, "\n");
  if(nwin_fcyk > 0) { 
    fprintf(ofp, "%-6s filter stage scan matrix overflows:         %15" PRId64 "  (%.4g)\n", 
	    "CYK", 
	    pli_acct->n_overflow_fcyk,
	    (double) pli_acct->n_overflow_fcyk / (double) nwin_fcyk);
  }
  else { 
    fprintf(ofp, "%-6s filter stage scan matrix overflows:         %15d  (%.4g)\n", 
	    "CYK", 
	    0, 0.);
  }
  if(nwin_final > 0) { 
    fprintf(ofp, "%-6s final  stage scan matrix overflows:         %15" PRId64 "  (%.4g)\n", 
	    (pli->final_cm_search_opts & CM_SEARCH_INSIDE) ? "Inside" : "CYK",
	    pli_acct->n_overflow_final,
	    (double) pli_acct->n_overflow_final / (double) nwin_final);
  }
  else { 
    fprintf(ofp, "%-6s final  stage scan matrix overflows:         %15d  (%.4g)\n", 
	    (pli->final_cm_search_opts & CM_SEARCH_INSIDE) ? "Inside" : "CYK",
	    0, 0.);
  }
  if (w != NULL) {
    esl_stopwatch_Display(ofp, w, "# CPU time: ");
    fprintf(ofp, "# Mc/sec: %.2f\n", 
	    (double) pli_acct->nres * (double) pli->nnodes / (w->elapsed * 1.0e6));
  }

  return eslOK;
}

  
/* Function:  cm_pli_ZeroAccounting()
 * Synopsis:  Zero a set of pipeline accounting statistics.
 * Incept:    EPN, Mon Nov 28 18:40:00 2011
 *
 * Returns:   <eslOK> on success.
 */
int
cm_pli_ZeroAccounting(CM_PLI_ACCT *pli_acct)
{
  pli_acct->nres              = 0;
  pli_acct->n_past_msv        = 0;
  pli_acct->n_past_vit        = 0;
  pli_acct->n_past_fwd        = 0;
  pli_acct->n_past_gfwd       = 0;
  pli_acct->n_past_edef       = 0;
  pli_acct->n_past_cyk        = 0;
  pli_acct->n_past_ins        = 0;
  pli_acct->n_output          = 0;
  pli_acct->n_past_msvbias    = 0;
  pli_acct->n_past_vitbias    = 0;
  pli_acct->n_past_fwdbias    = 0;
  pli_acct->n_past_gfwdbias   = 0;
  pli_acct->n_past_edefbias   = 0;
  pli_acct->pos_past_msv      = 0;
  pli_acct->pos_past_vit      = 0;
  pli_acct->pos_past_fwd      = 0;
  pli_acct->pos_past_gfwd     = 0;
  pli_acct->pos_past_edef     = 0;
  pli_acct->pos_past_cyk      = 0;
  pli_acct->pos_past_ins      = 0;      
  pli_acct->pos_output        = 0;
  pli_acct->pos_past_msvbias  = 0;
  pli_acct->pos_past_vitbias  = 0;
  pli_acct->pos_past_fwdbias  = 0;
  pli_acct->pos_past_gfwdbias = 0;
  pli_acct->pos_past_edefbias = 0;
  
  pli_acct->n_overflow_fcyk    = 0;
  pli_acct->n_overflow_final   = 0;
  pli_acct->n_aln_hb           = 0;
  pli_acct->n_aln_dccyk        = 0;

  return eslOK;
}


/* Function:  cm_pli_DescribePass()
 * Date:      EPN, Tue Nov 29 04:39:38 2011
 *
 * Purpose:   Translate internal flags for pipeline pass index
 *            into human-readable strings, for clearer output.
 * 
 * Args:      pass_idx - a pipeline pass index
 *                       PLI_PASS_SUMMED, PLI_PASS_STD, PLI_PASS_5P, PLI_PASS_3P PLI_PASS_5P_AND_3P
 *
 * Returns:   the appropriate string
 */
char *
cm_pli_DescribePass(CM_PIPELINE *pli, int pass_idx) 
{
  switch (pass_idx) {
    /*case PLI_PASS_SUMMED:    return "full sequences and terminii allowing for 5', 3' or 5' and 3' truncations";
      case PLI_PASS_STD:       return "full sequences, not allowing for truncations at terminii";
      case PLI_PASS_5P:        return "5' terminal sequence regions, specifically looking for 5' truncations";
      case PLI_PASS_3P:        return "3' terminal sequence regions, specifically looking for 3' truncations";
      case PLI_PASS_5P_AND_3P: return "short full sequences, specifically looking for hits that are both 5' and 3' truncated";
    */
  case PLI_PASS_SUMMED:    return "full sequences and terminii";
  case PLI_PASS_STD:       return "full sequences";
  case PLI_PASS_5P:        return "5' terminal sequence regions";
  case PLI_PASS_3P:        return "3' terminal sequence regions";
  case PLI_PASS_5P_AND_3P: return "full sequences short enough to contain a 5' and 3' truncated hit";
  default: cm_Fail("bogus pipeline pass index %d\n", pass_idx);
  }
  /* Note: we could return different strings depending on pli->research_ends, pli->do_trunc_ends, pli->do_force_ends,
   * but we don't currently.
   */
  return "";
}

/*------------------- end, pipeline API -------------------------*/

 
/* Function:  copy_subseq()
 * Incept:    EPN, Tue Aug  9 11:13:02 2011
 *
 * Purpose: Copy a subsequence of an existing sequence <src_sq>
 *           starting at position <i>, of length <L> to another
 *           sequence object <dest_sq>. Copy only residues from
 *           <i>..<i>+<L>-1. <dest_sq> must be pre-allocated.
 *
 * Returns: eslOK on success. 
 */
void
copy_subseq(const ESL_SQ *src_sq, ESL_SQ *dest_sq, int64_t i, int64_t L)
{ 
  /*printf("entering copy_subseq i: %" PRId64 " j: %" PRId64 " L: %" PRId64 " start-end: %" PRId64 "... %" PRId64 "\n",
    i, i+L-1, L, src_sq->start, src_sq->end);
    fflush(stdout);*/

  esl_sq_Reuse(dest_sq);
  esl_sq_GrowTo(dest_sq, L);
  memcpy((void*)(dest_sq->dsq+1), src_sq->dsq+i, L * sizeof(ESL_DSQ));
  dest_sq->dsq[0] = dest_sq->dsq[L+1] = eslDSQ_SENTINEL;
  dest_sq->n      = L;
  dest_sq->L      = src_sq->L;

  if(src_sq->start <= src_sq->end) { 
    ESL_DASSERT1((L <= (src_sq->end - src_sq->start + 1)));
    assert(L <= (src_sq->end - src_sq->start + 1));
    dest_sq->start = src_sq->start  + i - 1;
    dest_sq->end   = dest_sq->start + L - 1;
  }
  else { 
    ESL_DASSERT1((L <= (src_sq->start - src_sq->end + 1)));
    assert(L <= (src_sq->start - src_sq->end + 1));
    dest_sq->start = src_sq->end    + L - 1;
    dest_sq->end   = dest_sq->start - L + 1;
  }
  /*printf("leaving copy_subseq dest_sq->start..end: %" PRId64 " %" PRId64 " n: %" PRId64 "\n",
    dest_sq->start, dest_sq->end, dest_sq->n);*/

  esl_sq_SetName     (dest_sq, src_sq->name);
  esl_sq_SetAccession(dest_sq, src_sq->acc);
  esl_sq_SetDesc     (dest_sq, src_sq->desc);

  return;
}

/*****************************************************************
 * 3. Example 1: "search mode" in a sequence db
 *****************************************************************/
/*****************************************************************
 * 4. Example 2: "scan mode" in an HMM db
 *****************************************************************/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/

/************************************************************
 * @LICENSE@
 ************************************************************/

/* cm_wrappers.c
 * EPN, Wed Dec  6 06:11:46 2006
 * 
 * Wrapper functions for aligning and searching seqs
 * with a CM, with helper functions.
 * 
 */

#include "config.h"
#include "squidconf.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <math.h>

#include "easel.h"
#include "squid.h"		/* general sequence analysis library    */
#include "msa.h"                /* squid's multiple alignment i/o       */
#include "stopwatch.h"          /* squid's process timing module        */

#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* external functions                   */
#include "sre_stack.h"
#include "hmmband.h"         
#include "cm_postprob.h"
#include "stats.h"
#include "esl_gumbel.h"
#include "mpifuncs.h"

static db_seq_t *read_next_seq (ESL_SQFILE *dbfp, int do_revcomp);
static void print_results (CM_t *cm, CMConsensus_t *cons, db_seq_t *dbseq,
			   int do_complement, int do_stats, double *mu, 
			   double *lambda) ;
static int get_gc_comp(char *seq, int start, int stop);
static void remove_hits_over_e_cutoff (scan_results_t *results, char *seq,
				       float cutoff, double *lambda, double *mu);

/* coordinate -- macro that checks if it's reverse complement and if so 
   returns coordinate in original strand
   a = true if revcomp, false if not
   b = the position in current seq
   c = length of the seq
*/
#define coordinate(a,b,c) ( a ? -1*b+c+1 : b)

/* EPN, Tue Dec  5 14:25:02 2006
 * 
 * Function: AlignSeqsWrapper()
 * 
 * Purpose:  Given a CM, digitized sequences, and a slew of options, 
 *           do preliminaries, call the correct CYK function and return
 *           parsetrees and optionally postal codes (if do_post).
 * 
 * Args:     CM           - the covariance model
 *           dsq          - digitized sequences to align
 *           sqinfo       - info on the seq's we're aligning
 *           nseq         - number of seqs we're aligning
 *           ret_tr       - RETURN: parsetrees (pass NULL if trace isn't wanted)
 *           do_local     - TRUE to do local alignment, FALSE not to.
 *           do_small     - TRUE to use D&C CYK, FALSE not to
 *           do_qdb       - TRUE to use query dependet bands
 *           qdb_beta     - tail loss prob for QDB calculation
 *           do_hbanded   - TRUE use CP9 hmm derived target dependent bands
 *           use_sums     - TRUE to fill and use the posterior sums for CP9 band calculation 
 *           cp9bandp     - tail loss probability for CP9 hmm bands 
 *           do_sub       - TRUE to build and use a sub CM for alignment
 *           do_fullsub   - TRUE to build and use a full sub CM for alignment
 *           fsub_pmass   - probability mass to require in fullsub mode 
 *           do_hmmonly   - TRUE to align to the CP9 HMM, with viterbi
 *           do_inside    - TRUE to do Inside, and not return a parsetree
 *           do_outside   - TRUE to do Outside, and not return a parsetree
 *           do_check     - TRUE to check Inside and Outside probabilities
 *           do_post      - TRUE to do a posterior decode instead of CYK 
 *           ret_postcode - RETURN: postal code string, (NULL if do_post = FALSE)
 *           do_timings   - TRUE to report timings for alignment 
 *           bdump_level  - verbosity level for band related print statements
 *           debug_level  - verbosity level for debugging print statements
 *           silent_mode  - TRUE to not print anything, FALSE to print scores 
 *           do_enforce   - TRUE to read .enforce file and enforce MATL stretch 
 *           enf_start    - if (do_enforce), first MATL node to enforce each parse enter
 *           enf_end      - if (do_enforce), last  MATL node to enforce each parse enter
 *           do_elsilent  - disallow EL emissions
 * 
 *   Last 6 args are specific to partial-test.c (temporary?) these are usually NULL
 *           actual_spos  - [0..nseq-1] start consensus posn for truncated (partial) seq
 *           actual_epos  - [0..nseq-1] end   consensus posn for truncated (partial) seq
 *           ret_post_spos- [0..nseq-1] posterior probability from HMM of spos being start
 *           ret_post_epos- [0..nseq-1] posterior probability from HMM of epos being end
 *           ret_dist_spos- [0..nseq-1] distance (+/-) of max post start from spos
 *           ret_dist_epos- [0..nseq-1] distance (+/-) of max post end   from epos
 */
void
AlignSeqsWrapper(CM_t *cm, char **dsq, SQINFO *sqinfo, int nseq, Parsetree_t ***ret_tr, int do_local, 
		 int do_small, int do_qdb, double qdb_beta,
		 int do_hbanded, int use_sums, double hbandp, int do_sub, int do_fullsub, float fsub_pmass,
		 int do_hmmonly, int do_inside, int do_outside, int do_check, int do_post, 
		 char ***ret_postcode, int do_timings, int bdump_level, int debug_level, int silent_mode, 
		 int do_enforce, int enf_start, int enf_end, int do_elsilent,
		 int *actual_spos, int *actual_epos, float **ret_post_spos, float **ret_post_epos,
		 int **ret_dist_spos, int **ret_dist_epos)
{
  Stopwatch_t  *watch1, *watch2;      /* for timings */
  int i;                              /* counter over sequences */
  int v;                              /* state counter */
  char           **postcode;    /* posterior decode array of strings        */
  Parsetree_t    **tr;          /* parse trees for the sequences */
  float            sc;		/* score for one sequence alignment */
  float            maxsc;	/* max score in all seqs */
  float            minsc;	/* min score in all seqs */
  float            avgsc;	/* avg score over all seqs */
  int              nd;          /* counter over nodes */
  /* variables related to CM Plan 9 HMMs */
  struct cplan9_s       *hmm;           /* constructed CP9 HMM */
  CP9Bands_t *cp9b;                     /* data structure for hmm bands (bands on the hmm states) 
				         * and arrays for CM state bands, derived from HMM bands*/
  CP9Map_t              *cp9map;        /* maps the hmm to the cm and vice versa */
  struct cp9_dpmatrix_s *cp9_mx;        /* growable DP matrix for viterbi                       */
  struct cp9_dpmatrix_s *cp9_fwd;       /* growable DP matrix for forward                       */
  struct cp9_dpmatrix_s *cp9_bck;       /* growable DP matrix for backward                      */
  struct cp9_dpmatrix_s *cp9_posterior; /* growable DP matrix for posterior decode              */
  float                  swentry;	/* S/W aggregate entry probability       */
  float                  swexit;        /* S/W aggregate exit probability        */
  float forward_sc; 
  float backward_sc; 

  /* variables related to the do_sub option */
  CM_t              *sub_cm;       /* sub covariance model                      */
  CMSubMap_t        *submap;
  CP9Bands_t        *sub_cp9b;     /* data structure for hmm bands (bands on the hmm states) 
				    * and arrays for CM state bands, derived from HMM bands */
  CM_t              *orig_cm;      /* the original, template covariance model the sub CM was built from */
  int                spos;         /* HMM node most likely to have emitted posn 1 of target seq */
  int                spos_state;   /* HMM state type for curr spos 0=match or 1=insert */
  int                epos;         /* HMM node most likely to have emitted posn L of target seq */
  int                epos_state;   /* HMM state type for curr epos 0=match or  1=insert */
  Parsetree_t     *orig_tr;        /* parsetree for the orig_cm; created from the sub_cm parsetree */

  struct cplan9_s *sub_hmm;        /* constructed CP9 HMM; written to hmmfile              */
  CP9Map_t        *sub_cp9map;     /* maps the sub_hmm to the sub_cm and vice versa */
  struct cplan9_s *orig_hmm;       /* original CP9 HMM built from orig_cm */
  CP9Map_t        *orig_cp9map;    

  /* variables related to query dependent banding (qdb) */
  int    expand_flag;           /* TRUE if the dmin and dmax vectors have just been 
				 * expanded (in which case we want to recalculate them 
				 * before we align a new sequence), and FALSE if not*/
  double **gamma;               /* P(subseq length = n) for each state v    */
  int     *dmin;                /* minimum d bound for state v, [0..v..M-1] */
  int     *dmax;                /* maximum d bound for state v, [0..v..M-1] */
  int *orig_dmin;               /* original dmin values passed in */
  int *orig_dmax;               /* original dmax values passed in */
  int safe_windowlen; 

  /* variables related to inside/outside */
  float           ***alpha;     /* alpha DP matrix for Inside() */
  float           ***beta;      /* beta DP matrix for Inside() */
  float           ***post;      /* post DP matrix for Inside() */

  /* partial-test variables */
  int do_ptest;                 /* TRUE to fill partial-test variables */
  float *post_spos;
  float *post_epos;
  int   *dist_spos;
  int   *dist_epos;

  if(do_fullsub)
    do_sub = TRUE;

  /*printf("in AlignSeqsWrapper() do_local: %d do_sub: %d do_fullsub: %d\n", do_local, do_sub, do_fullsub);*/

  do_ptest = FALSE;
  if(ret_post_spos != NULL)
    do_ptest = TRUE;
  if(do_ptest && (ret_post_epos == NULL || ret_dist_spos == NULL || ret_dist_epos == NULL))
    Die("ERROR partial-test arrays must either all be NULL or non-NULL\n");

  /* Allocate partial-test arrays */
  if(ret_post_spos != NULL)
    post_spos = MallocOrDie(sizeof(float) * nseq);
  if(ret_post_epos != NULL)
    post_epos = MallocOrDie(sizeof(float) * nseq);
  if(ret_dist_spos != NULL)
    dist_spos = MallocOrDie(sizeof(int  ) * nseq);
  if(ret_dist_epos != NULL)
    dist_epos = MallocOrDie(sizeof(int  ) * nseq);

  tr    = MallocOrDie(sizeof(Parsetree_t) * nseq);
  minsc = FLT_MAX;
  maxsc = -FLT_MAX;
  avgsc = 0;

  watch1 = StopwatchCreate(); /* watch1 is used to time each step individually */
  watch2 = StopwatchCreate(); /* watch2 times the full alignment (including band calc)
				 for each seq */

  if(do_hbanded || do_sub || do_ptest) /* We need a CP9 HMM to build sub_cms */
    {
      /* Ensure local begins and ends in the CM are off */
      /* TO DO: write a function that puts CM back in glocal mode */

      if(!build_cp9_hmm(cm, &hmm, &cp9map, FALSE, 0.0001, debug_level))
	Die("Couldn't build a CP9 HMM from the CM\n");
      cp9_mx  = CreateCPlan9Matrix(1, hmm->M, 25, 0);

      /* Keep this data for the original CM safe; we'll be doing
       * pointer swapping to ease the sub_cm alignment implementation. */
      orig_hmm = hmm;
      orig_cp9map = cp9map;
      if(do_hbanded)
	cp9b = AllocCP9Bands(cm, hmm);

      StopwatchZero(watch2);
      StopwatchStart(watch2);
    }

  /* Relocated ConfigLocal() call to here, below the CM Plan 9 construction.
   * Otherwise its impossible to make a CM Plan 9 HMM from the local CM
   * that passes the current tests to ensure the HMM is "close enough" to
   * the CM. This is something to look into later.
   */
  if (do_local)
    { 
      if(do_enforce)
	ConfigLocalEnforce(cm, 0.5, 0.5, enf_start, enf_end);
      else
	ConfigLocal(cm, 0.5, 0.5);
      CMLogoddsify(cm);
      /*CMHackInsertScores(cm);*/	/* "TEMPORARY" fix for bad priors */
    }

  if(do_elsilent) 
    ConfigLocal_DisallowELEmissions(cm);

  /* the --enforce option, added specifically for enforcing the template region of
   * telomerase RNA */

    {
      printf("Enforcing MATL stretch from %d to %d.\n", enf_start, enf_end);
      /* Configure local alignment so the MATL stretch is unavoidable */
      ConfigLocalEnforce(cm, 0.5, 0.5, enf_start, enf_end);
      CMLogoddsify(cm);
      printf("Done enforcing.\n");
    }

  if((do_local && do_hbanded) && !do_sub)
      {
	/*printf("configuring the CM plan 9 HMM for local alignment.\n");*/
	CPlan9SWConfig(hmm, 0.5, 0.5);
	CP9Logoddsify(hmm);

      }
  if(do_sub || do_ptest) /* to get spos and epos for the sub_cm, 
			  * we config the HMM to local mode with equiprobable start/end points.*/
      {
	/*printf("configuring the CM plan 9 HMM for local alignment.\n");*/
	swentry= ((hmm->M)-1.)/hmm->M; /* all start pts equiprobable, including 1 */
	swexit = ((hmm->M)-1.)/hmm->M; /* all end   pts equiprobable, including M */
	CPlan9SWConfig(hmm, swentry, swexit);
	CP9Logoddsify(hmm);
	orig_tr    = MallocOrDie(sizeof(Parsetree_t));
      }
  /* set up the query dependent bands, this has to be done after the ConfigLocal() call */
  if(do_qdb || bdump_level > 0)
    {
      safe_windowlen = cm->W * 2;
      while(!(BandCalculationEngine(cm, safe_windowlen, qdb_beta, 0, &dmin, &dmax, &gamma, do_local)))
	{
	  FreeBandDensities(cm, gamma);
	  free(dmin);
	  free(dmax);
	  safe_windowlen *= 2;
	  /*printf("ERROR BandCalculationEngine returned false, windowlen adjusted to %d\n", safe_windowlen);*/
	}
      /* If we're enforcing a subsequence, we need to reenforce it b/c BandCalculationEngine() 
       * changes the local end probabilities */
      if(do_enforce && do_local)
	{
	  ConfigLocalEnforce(cm, 0.5, 0.5, enf_start, enf_end);
	  CMLogoddsify(cm);
	}
      if(bdump_level > 1) 
	  /*printf("qdb_beta:%f\n", qdb_beta);*/
	  debug_print_bands(cm, dmin, dmax);
      expand_flag = FALSE;
      /* Copy dmin and dmax, so we can replace them after expansion */
      orig_dmin = MallocOrDie(sizeof(int) * cm->M);
      orig_dmax = MallocOrDie(sizeof(int) * cm->M);
      for(v = 0; v < cm->M; v++)
	{
	  orig_dmin[v] = dmin[v];
	  orig_dmax[v] = dmax[v];
	}
    }	  
  if(do_post)
    postcode = malloc(sizeof(char *) * nseq);

  orig_cm = cm;

  /*****************************************************************
   *  Collect parse trees for each sequence
   *****************************************************************/

  for (i = 0; i < nseq; i++)
    {
      StopwatchZero(watch1);
      StopwatchStart(watch1);
      StopwatchZero(watch2);
      StopwatchStart(watch2);
      
      if (sqinfo[i].len == 0) Die("ERROR: sequence named %s has length 0.\n", sqinfo[i].name);

      /* Potentially, do HMM calculations. */
      if(do_hbanded || do_sub || do_ptest)
	{
	  /* We want HMM posteriors for this sequence to the full length (non-sub) HMM */
	  StopwatchZero(watch1);
	  StopwatchStart(watch1);

	  /* Step 1: Get HMM posteriors.*/
	  /*sc = CP9Viterbi(dsq[i], 1, sqinfo[i].len, hmm, cp9_mx);*/
	  forward_sc = CP9Forward(dsq[i], 1, sqinfo[i].len, orig_hmm, &cp9_fwd);
	  if(debug_level > 0) printf("CP9 i: %d | forward_sc : %.2f\n", i, forward_sc);
	  backward_sc = CP9Backward(dsq[i], 1, sqinfo[i].len, orig_hmm, &cp9_bck);
	  if(debug_level > 0) printf("CP9 i: %d | backward_sc: %.2f\n", i, backward_sc);
	  
	  /*debug_check_CP9_FB(cp9_fwd, cp9_bck, hmm, forward_sc, 1, sqinfo[i].len, dsq[i]);*/
	  cp9_posterior = cp9_bck;
	  CP9FullPosterior(dsq[i], 1, sqinfo[i].len, orig_hmm, cp9_fwd, cp9_bck, cp9_posterior);
	}

      if(do_ptest) /* determine the posterior probability from HMM of the correct start posn
		    * and end posn, as well as distance from max posterior */
	{
	  /* Determine HMM post probability of actual start and end */
	  /*	  CP9NodeForPosn(orig_hmm, 1, sqinfo[i].len, actual_spos[i], cp9_posterior, 
			 &spos, &spos_state, &(post_spos[i]), debug_level);
	  CP9NodeForPosn(orig_hmm, 1, sqinfo[i].len, actual_epos[i], cp9_posterior, 
	  &epos, &spos_state, &(post_epos[i]), debug_level);*/
	  printf("(%4d) s: %3d post: %.2f\n", i, actual_spos[i], 
		 Score2Prob(cp9_posterior->mmx[1][actual_spos[i]], 1.));
	  printf("(%4d) e: %3d post: %.2f\n", i, actual_epos[i], 
		 Score2Prob(cp9_posterior->mmx[sqinfo[i].len][actual_epos[i]], 1.));
	  /* Determine HMM post probability of most likely start and end */
	  /*CP9NodeForPosn(orig_hmm, 1, sqinfo[i].len, 1, cp9_posterior, 
			 &spos, &spos_state, NULL, debug_level);
	  CP9NodeForPosn(orig_hmm, 1, sqinfo[i].len, sqinfo[i].len, cp9_posterior, 
	  &epos, &epos_state, NULL, debug_level);*/
	}
      
      /* If we're in sub mode:
       * (1) Get HMM posteriors. (we already did this above)
       * (2) Infer the start (spos) and end (epos) HMM states by 
       *     looking at the posterior matrix.
       * (3) Build the sub_cm from the original CM.
       *
       * If we're also doing HMM banded alignment:
       * (4) Build a new CP9 HMM from the sub CM.
       * (5) Do Forward/Backward again, and get a new posterior matrix.
       */
      if(do_sub)
	{
	  /* (2) infer the start and end HMM states by looking at the posterior matrix.
	   * Remember: we're necessarily in local mode, the --sub option turns local mode on. 
	   */
	  CP9NodeForPosn(orig_hmm, 1, sqinfo[i].len, 1,             cp9_posterior, &spos, &spos_state, 
			 do_fullsub, fsub_pmass, TRUE, debug_level);
	  CP9NodeForPosn(orig_hmm, 1, sqinfo[i].len, sqinfo[i].len, cp9_posterior, &epos, &epos_state, 
			 do_fullsub, fsub_pmass, FALSE, debug_level);
	  /* If the most likely state to have emitted the first or last residue
	   * is the insert state in node 0, it only makes sense to start modelling
	   * at consensus column 1. */
	  if(spos == 0 && spos_state == 1) 
	      spos = 1;
	  if(epos == 0 && epos_state == 1) 
	      epos = 1;
	  if(epos < spos) /* This is a possible but hopefully rarely encountered situation. */
	    epos = spos;
	  
	  /* (3) Build the sub_cm from the original CM. */
	  if(!(build_sub_cm(orig_cm, &sub_cm, 
			    spos, epos,         /* first and last col of structure kept in the sub_cm  */
			    &submap,            /* maps from the sub_cm to cm and vice versa           */
			    do_fullsub,         /* build or not build a sub CM that models all columns */
			    debug_level)))      /* print or don't print debugging info                 */
	    Die("Couldn't build a sub CM from the CM\n");
	  cm    = sub_cm; /* orig_cm still points to the original CM */
	  /* If the sub_cm models the full consensus length of the orig_cm, with only
	   * structure removed, we configure it for local alignment to allow it to 
	   * skip the single stranded regions at the beginning and end. But only 
	   * if we don't need to build a CP9 HMM from the sub_cm to do banded alignment.*/
	  if(do_fullsub && !do_hbanded)
	    {
	      printf("calling ConfigLocal_fullsub_post()\n");
	      /* FIX THIS WHOLE THING */
	      ConfigLocal_fullsub_post(sub_cm, orig_cm, orig_cp9map, submap, cp9_posterior, sqinfo[i].len);
	      /*ConfigLocal_fullsub(cm, 0.5, 0.5, orig_cp9map->pos2nd[submap->sstruct],
		orig_cp9map->pos2nd[submap->estruct]);*/
	      /*ConfigLocal(sub_cm, 0.5, 0.5);*/
	      /*printf("DEBUG PRINTING CM PARAMS AFTER CONFIGLOCAL_FULLSUB_POST CALL\n");
		debug_print_cm_params(cm);
		printf("DONE DEBUG PRINTING CM PARAMS AFTER CONFIGLOCAL_FULLSUB_POST CALL\n");*/
	      CMLogoddsify(cm);
	      do_local = TRUE; /* we wait til we get here to set do_local, if we 
				* configure for local alignment earlier it would've 
				* screwed up CP9 construction. */
	    }	  
	  if(do_hbanded) /* we're doing HMM banded alignment to the sub_cm */
	    {
	      /* (4) Build a new CP9 HMM from the sub CM. */
	      /* Eventually, I think we can do this by just adjusting the parameters of the original HMM 
		 CP9_2sub_cp9(hmm, &sub_hmm2, spos, epos, orig_phi);
	      */
	      if(!build_cp9_hmm(sub_cm, &sub_hmm, &sub_cp9map, FALSE, 0.0001, debug_level))
		Die("Couldn't build a sub CP9 HMM from the sub CM\n");

	      /* Allocate HMM banding data structures for use with the sub CM and sub HMM */
	      sub_cp9b = AllocCP9Bands(sub_cm, sub_hmm);
	      
	      if(do_fullsub)
		{
		  /* FIX THIS WHOLE THING! */
		  CPlan9SWConfig(sub_hmm, 0.5, 0.5);
		  CP9Logoddsify(sub_hmm);
		  ConfigLocal_fullsub(sub_cm, 0.5, 0.5, sub_cp9map->pos2nd[submap->sstruct],
				      sub_cp9map->pos2nd[submap->estruct]);
		  /*ConfigLocal(sub_cm, 0.5, 0.5);*/
		  /*printf("debug printing sub cm params after config local full sub:\n");
		  debug_print_cm_params(sub_cm);
		  printf("done debug printing sub cm params after config local full sub:\n");*/
		  
		  CMLogoddsify(cm);
		  do_local = TRUE;
		}
	      /* (5) Do Forward/Backward again, and get a new posterior matrix. 
	       * We have to free cp9_fwd and cp9_posterior because we used them 
	       * to find spos and epos. */

	      FreeCPlan9Matrix(cp9_fwd);
	      FreeCPlan9Matrix(cp9_posterior);
	      forward_sc = CP9Forward(dsq[i], 1, sqinfo[i].len, sub_hmm, &cp9_fwd);
	      if(debug_level) printf("CP9 i: %d | forward_sc : %.2f\n", i, forward_sc);
	      backward_sc = CP9Backward(dsq[i], 1, sqinfo[i].len, sub_hmm, &cp9_bck);
	      if(debug_level) printf("CP9 i: %d | backward_sc: %.2f\n", i, backward_sc);
	      /*debug_check_CP9_FB(cp9_fwd, cp9_bck, hmm, forward_sc, 1, sqinfo[i].len, dsq[i]);*/
	      cp9_posterior = cp9_bck;
	      CP9FullPosterior(dsq[i], 1, sqinfo[i].len, sub_hmm, cp9_fwd, cp9_bck, cp9_posterior);
	      /* cp9_posterior has the posteriors for the sub_hmm */

	      /* Change some pointers so that the functions that create bands use the
	       * sub_* data structures. The orig_* data structures will still point
	       * to the original CM versions. */
	      hmm           = sub_hmm;    
	      cp9map        = sub_cp9map;
	      cp9b          = sub_cp9b;
	    }
	}
      if(do_hbanded)
	{
	  StopwatchStop(watch1);
	  if(do_timings) StopwatchDisplay(stdout, "CP9 Forward/Backward CPU time: ", watch1);
	  StopwatchZero(watch1);
	  StopwatchStart(watch1);
      
	  /* Align the current seq to the cp9 HMM, we don't care
	   * about the trace, just the posteriors.
	   * Step 1: Get HMM posteriors. (if do_sub, we already did this above,
	   *                              the posteriors are for the sub_hmm)
	   * Step 2: posteriors -> HMM bands.
	   * Step 3: HMM bands  ->  CM bands.
	   */
	  
	  /* Step 2: posteriors -> HMM bands.*/
	  if(use_sums)
	    CP9_ifill_post_sums(cp9_posterior, 1, sqinfo[i].len, cp9b->hmm_M,
				cp9b->isum_pn_m, cp9b->isum_pn_i, cp9b->isum_pn_d);
	  /* match states */
	  CP9_hmm_band_bounds(cp9_posterior->mmx, 1, sqinfo[i].len, cp9b->hmm_M, 
			      cp9b->isum_pn_m, cp9b->pn_min_m, cp9b->pn_max_m,
			      (1.-hbandp), HMMMATCH, use_sums, debug_level);
	  /* insert states */
	  CP9_hmm_band_bounds(cp9_posterior->imx, 1, sqinfo[i].len, cp9b->hmm_M,
			      cp9b->isum_pn_i, cp9b->pn_min_i, cp9b->pn_max_i,
			      (1.-hbandp), HMMINSERT, use_sums, debug_level);
	  /* delete states (note: delete_flag set to TRUE) */
	  CP9_hmm_band_bounds(cp9_posterior->dmx, 1, sqinfo[i].len, cp9b->hmm_M,
			      cp9b->isum_pn_d, cp9b->pn_min_d, cp9b->pn_max_d,
			      (1.-hbandp), HMMDELETE, use_sums, debug_level);

	  if(debug_level != 0)
	    {
	      printf("printing hmm bands\n");
	      print_hmm_bands(stdout, sqinfo[i].len, cp9b->hmm_M, cp9b->pn_min_m, 
			      cp9b->pn_max_m, cp9b->pn_min_i, cp9b->pn_max_i, 
			      cp9b->pn_min_d, cp9b->pn_max_d, hbandp, debug_level);
	    }
	  
	  /* Step 3: HMM bands  ->  CM bands. */
	  hmm2ij_bands(cm, cp9map, 1, sqinfo[i].len, cp9b->pn_min_m, cp9b->pn_max_m, 
		       cp9b->pn_min_i, cp9b->pn_max_i, cp9b->pn_min_d, cp9b->pn_max_d, 
		       cp9b->imin, cp9b->imax, cp9b->jmin, cp9b->jmax, debug_level);
	  
	  StopwatchStop(watch1);
	  if(do_timings) StopwatchDisplay(stdout, "CP9 Band calculation CPU time: ", watch1);
	  /* Use the CM bands on i and j to get bands on d, specific to j. */
	  for(v = 0; v < cm->M; v++)
	    {
	      cp9b->hdmin[v] = malloc(sizeof(int) * (cp9b->jmax[v] - cp9b->jmin[v] + 1));
	      cp9b->hdmax[v] = malloc(sizeof(int) * (cp9b->jmax[v] - cp9b->jmin[v] + 1));
	    }
	  ij2d_bands(cm, sqinfo[i].len, cp9b->imin, cp9b->imax, cp9b->jmin, cp9b->jmax,
		     cp9b->hdmin, cp9b->hdmax, -1);
	  
	  if(debug_level != 0)
	    PrintDPCellsSaved_jd(cm, cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax, 
				 sqinfo[i].len);
	  
	  FreeCPlan9Matrix(cp9_fwd);
	  FreeCPlan9Matrix(cp9_posterior);
	  /* Done with the HMM. On to the CM. */
	}
      
      /* Determine which CYK alignment algorithm to use, based
       * on command-line options AND memory requirements.
       */
      if(do_hbanded)
	{
	  /* write a function to determine size of jd banded memory
	   * req'd, and set do_small to true if its > thresh.
	   if(do_small) * We're only going to band on d in memory, but 
	   * we need to calculate safe_hd bands on the d dimension. 
	   {
	  */
	}
      
      if(do_qdb)
	{
	  /*Check if we need to reset the query dependent bands b/c they're currently expanded. */
	  if(expand_flag)
	    {
	      for(v = 0; v < cm->M; v++)
		{
		  dmin[v] = orig_dmin[v];
		  dmax[v] = orig_dmax[v];
		}
	      expand_flag = FALSE;
	    }
	  if((sqinfo[i].len < dmin[0]) || (sqinfo[i].len > dmax[0]))
	    {
	      /* the seq we're aligning is outside the root band, so we expand.*/
	      ExpandBands(cm, sqinfo[i].len, dmin, dmax);
	      if(debug_level > 0) printf("Expanded bands for seq : %s\n", sqinfo[i].name);
	      if(bdump_level > 2) 
		{
		  printf("printing expanded bands :\n");
		  debug_print_bands(cm, dmin, dmax);
		}
	      expand_flag = TRUE;
	    }
	}

      if(!silent_mode) 
	{
	  if(do_sub) 
	    printf("Aligning (to a sub CM) %-20s", sqinfo[i].name);
	  else
	    printf("Aligning %-30s", sqinfo[i].name);
	}
      if (do_inside)
	{
	  if(do_hbanded)
	    sc = FInside_b_jd_me(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
				 BE_PARANOID,	/* non memory-saving mode */
				 NULL, NULL,	/* manage your own matrix, I don't want it */
				 NULL, NULL,	/* manage your own deckpool, I don't want it */
				 do_local,        /* TRUE to allow local begins */
				 cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax); /* j and d bands */
	  else
	    sc = FInside(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
			 BE_EFFICIENT,	/* memory-saving mode */
			 NULL, NULL,	/* manage your own matrix, I don't want it */
			 NULL, NULL,	/* manage your own deckpool, I don't want it */
			 do_local);       /* TRUE to allow local begins */
	}
      else if(do_outside)
	{	
	  if(do_hbanded)
	    {
	      sc = FInside_b_jd_me(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
				   BE_PARANOID,	/* save full alpha so we can run outside */
				   NULL, &alpha,	/* fill alpha, and return it, needed for FOutside() */
				   NULL, NULL,	/* manage your own deckpool, I don't want it */
				   do_local,        /* TRUE to allow local begins */
				   cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax); /* j and d bands */
	      /*do_check = TRUE;*/
	      sc = FOutside_b_jd_me(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
				    BE_PARANOID,	/* save full beta */
				    NULL, NULL,	/* manage your own matrix, I don't want it */
				    NULL, NULL,	/* manage your own deckpool, I don't want it */
				    do_local,       /* TRUE to allow local begins */
				    alpha,          /* alpha matrix from FInside_b_jd_me() */
				    NULL,           /* don't save alpha */
				    do_check,       /* TRUE to check Outside probs agree with Inside */
				    cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax); /* j and d bands */
	    }
	  else
	    {
	      sc = FInside(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
			   BE_PARANOID,	/* save full alpha so we can run outside */
			   NULL, &alpha,	/* fill alpha, and return it, needed for FOutside() */
			   NULL, NULL,	/* manage your own deckpool, I don't want it */
			   do_local);       /* TRUE to allow local begins */
	      sc = FOutside(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
			    BE_PARANOID,	/* save full beta */
			    NULL, NULL,	/* manage your own matrix, I don't want it */
			    NULL, NULL,	/* manage your own deckpool, I don't want it */
			    do_local,       /* TRUE to allow local begins */
			    alpha,          /* alpha matrix from FInside() */
			    NULL,           /* don't save alpha */
			    do_check);      /* TRUE to check Outside probs agree with Inside */
	    }
	}
      else if (do_small) 
	{
	  if(do_qdb)
	    {
	      sc = CYKDivideAndConquer(cm, dsq[i], sqinfo[i].len, 0, 1, sqinfo[i].len, 
				       &(tr[i]), dmin, dmax);
	      if(bdump_level > 0)
 		qdb_trace_info_dump(cm, tr[i], dmin, dmax, bdump_level);
	    }
	  else if(do_hbanded) /*j and d bands not tight enough to allow HMM banded full CYK*/
	    {
	      /* Calc the safe d bands */
	      hd2safe_hd_bands(cm->M, cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax, 
			       cp9b->safe_hdmin, cp9b->safe_hdmax);
	      if(debug_level > 3)
		{
		  printf("\nprinting hd bands\n\n");
		  debug_print_hd_bands(cm, cp9b->hdmin, cp9b->hdmax, cp9b->jmin, cp9b->jmax);
		  printf("\ndone printing hd bands\n\n");
		}
	      /* Note the following CYK call will not enforce j bands, even
	       * though user specified --hbanded. */
	      sc = CYKDivideAndConquer(cm, dsq[i], sqinfo[i].len, 0, 1, sqinfo[i].len, 
				       &(tr[i]), cp9b->safe_hdmin, cp9b->safe_hdmax);
	      if(bdump_level > 0)
		qdb_trace_info_dump(cm, tr[i], dmin, dmax, bdump_level);
	    }
	  else
	    {
	      /*printf("DEBUG PRINTING CM PARAMS BEFORE D&C CALL\n");
		debug_print_cm_params(cm);
		printf("DONE DEBUG PRINTING CM PARAMS BEFORE D&C CALL\n");*/

	      sc = CYKDivideAndConquer(cm, dsq[i], sqinfo[i].len, 0, 1, sqinfo[i].len, &(tr[i]),
				       NULL, NULL); /* we're not in QDB mode */
	      if(bdump_level > 0)
		{
		  /* We want band info but --banded wasn't used.  Useful if you're curious
		   * why a banded parse is crappy relative to non-banded parse, e.g. allows you 
		   * to see where the non-banded parse went outside the bands.
		   */
		  qdb_trace_info_dump(cm, tr[i], dmin, dmax, bdump_level);
		}
	    }
	}
      else if(do_qdb)
	{
	  sc = CYKInside(cm, dsq[i], sqinfo[i].len, 0, 1, sqinfo[i].len, &(tr[i]), dmin, dmax);
	  if(bdump_level > 0)
	    qdb_trace_info_dump(cm, tr[i], dmin, dmax, bdump_level);
	}
      else if(do_hbanded)
	{
	  sc = CYKInside_b_jd(cm, dsq[i], sqinfo[i].len, 0, 1, sqinfo[i].len, &(tr[i]), cp9b->jmin, 
			      cp9b->jmax, cp9b->hdmin, cp9b->hdmax, cp9b->safe_hdmin, cp9b->safe_hdmax);
	  if(bdump_level > 0)
	    qdb_trace_info_dump(cm, tr[i], cp9b->safe_hdmin, cp9b->safe_hdmax, bdump_level);
	}
      else
	{
	  sc = CYKInside(cm, dsq[i], sqinfo[i].len, 0, 1, sqinfo[i].len, &(tr[i]), NULL, NULL);
	  if(bdump_level > 0)
	    {
	      /* We want band info but --hbanded wasn't used.  Useful if you're curious
	       * why a banded parse is crappy relative to non-banded parse, e.g. allows you 
	       * to see where the non-banded parse went outside the bands.
	       */
	      qdb_trace_info_dump(cm, tr[i], dmin, dmax, bdump_level);
	    }
	}
      if(do_post) /* Do Inside() and Outside() runs and use alpha and beta to get posteriors */
	{	
	  /*alpha = MallocOrDie(sizeof(float **) * (cm->M));
	  beta  = MallocOrDie(sizeof(float **) * (cm->M+1));
	  */
	  post  = MallocOrDie(sizeof(float **) * (cm->M+1));
	  /*
	  for (v = 0; v < cm->M; v++) alpha[v] = NULL;
	  for (v = 0; v < cm->M+1; v++) beta[v] = NULL;
	  */
	  if(do_hbanded)
	    {
	      for (v = 0; v < cm->M; v++)
		{
		  post[v] = NULL;
		  post[v] = alloc_jdbanded_vjd_deck(sqinfo[i].len, 1, sqinfo[i].len, cp9b->jmin[v], 
						    cp9b->jmax[v], cp9b->hdmin[v], cp9b->hdmax[v]);
		}
	      post[cm->M] = NULL;
	      post[cm->M] = alloc_vjd_deck(sqinfo[i].len, 1, sqinfo[i].len);
	      sc = FInside_b_jd_me(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
				   BE_PARANOID,	/* save full alpha so we can run outside */
				   NULL, &alpha,	/* fill alpha, and return it, needed for FOutside() */
				   NULL, NULL,	/* manage your own deckpool, I don't want it */
				   do_local,       /* TRUE to allow local begins */
				   cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax); /* j and d bands */
	      sc = FOutside_b_jd_me(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
				    BE_PARANOID,	/* save full beta */
				    NULL, &beta,	/* fill beta, and return it, needed for CMPosterior() */
				    NULL, NULL,	/* manage your own deckpool, I don't want it */
				    do_local,       /* TRUE to allow local begins */
				    alpha, &alpha,  /* alpha matrix from FInside(), and save it for CMPosterior*/
				    do_check,      /* TRUE to check Outside probs agree with Inside */
				    cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax); /* j and d bands */
	      CMPosterior_b_jd_me(sqinfo[i].len, cm, alpha, NULL, beta, NULL, post, &post,
				  cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax);
	      postcode[i] = CMPostalCode_b_jd_me(cm, sqinfo[i].len, post, tr[i],
						 cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax);
	    }
	  else
	    {
	      for (v = 0; v < cm->M+1; v++)
		{
		  post[v] = NULL;
		  post[v] = alloc_vjd_deck(sqinfo[i].len, 1, sqinfo[i].len);
		}
	      sc = FInside(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
			   BE_PARANOID,	/* save full alpha so we can run outside */
			   NULL, &alpha,	/* fill alpha, and return it, needed for FOutside() */
			   NULL, NULL,	/* manage your own deckpool, I don't want it */
			   do_local);       /* TRUE to allow local begins */
	      sc = FOutside(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
			    BE_PARANOID,	/* save full beta */
			    NULL, &beta,	/* fill beta, and return it, needed for CMPosterior() */
			    NULL, NULL,	/* manage your own deckpool, I don't want it */
			    do_local,       /* TRUE to allow local begins */
			    alpha, &alpha,  /* alpha matrix from FInside(), and save it for CMPosterior*/
			    do_check);      /* TRUE to check Outside probs agree with Inside */
	      CMPosterior(sqinfo[i].len, cm, alpha, NULL, beta, NULL, post, &post);
	      if(do_check)
		{
		  CMCheckPosterior(sqinfo[i].len, cm, post);
		  printf("\nPosteriors checked.\n\n");
		}
	      postcode[i] = CMPostalCode(cm, sqinfo[i].len, post, tr[i]);
	    }

	  /* free post */
	  if(post != NULL)
	    {
	      for (v = 0; v <= (cm->M); v++)
		if (post[v] != NULL) { free_vjd_deck(post[v], 1, sqinfo[i].len); post[v] = NULL;}
	      free(post);
	    }
	}
      avgsc += sc;
      if (sc > maxsc) maxsc = sc;
      if (sc < minsc) minsc = sc;
      
      if(!silent_mode) printf("    score: %10.2f bits\n", sc);
      
      /* If debug level high enough, print out the parse tree */
      if(debug_level > 2)
	{
	  fprintf(stdout, "  SCORE : %.2f bits\n", ParsetreeScore(cm, tr[i], dsq[i], FALSE));;
	  ParsetreeDump(stdout, tr[i], cm, dsq[i]);
	  fprintf(stdout, "//\n");
	}
      /* Dump the trace with info on i, j and d bands
       * if bdump_level is high enough */
      if(bdump_level > 0 && do_hbanded)
	ijd_banded_trace_info_dump(cm, tr[i], cp9b->imin, cp9b->imax, cp9b->jmin, cp9b->jmax, 
				   cp9b->hdmin, cp9b->hdmax, 1);
      
      /* Clean up the structures we use calculating HMM bands, that are allocated
       * differently for each sequence. 
       */
      if(do_hbanded)
	{
	  for(v = 0; v < cm->M; v++)
	    { 
	      free(cp9b->hdmin[v]); 
	      free(cp9b->hdmax[v]);
	    }
	  StopwatchStop(watch2);
	  if(do_timings) 
	    { 
	      StopwatchDisplay(stdout, "band calc and jd CYK CPU time: ", watch2);
	      printf("\n");
	    }
	}
      if(do_sub && !(do_inside || do_outside))
	{
	  /* Convert the sub_cm parsetree to a full CM parsetree */
	  if(debug_level > 0)
	    ParsetreeDump(stdout, tr[i], cm, dsq[i]);
	  if(!(sub_cm2cm_parsetree(orig_cm, sub_cm, &orig_tr, tr[i], submap, do_fullsub, debug_level)))
	    {
	      printf("\n\nIncorrectly converted original trace:\n");
	      ParsetreeDump(stdout, orig_tr, orig_cm, dsq[i]);
	      exit(1);
	    }
	  if(debug_level > 0)
	    {
	      printf("\n\nConverted original trace:\n");
	      ParsetreeDump(stdout, orig_tr, orig_cm, dsq[i]);
	    }
	  /* Replace the sub_cm trace with the converted orig_cm trace. */
	  FreeParsetree(tr[i]);
	  tr[i] = orig_tr;
	  
	  FreeSubMap(submap);
	  FreeCM(sub_cm); /* cm and sub_cm now point to NULL */
	  if(do_hbanded)
	    {
	      FreeCP9Map(sub_cp9map);
	      FreeCPlan9(sub_hmm);
	      FreeCP9Bands(sub_cp9b);
	    }
	}
    }
  /* Clean up. */
  if((do_sub && !do_hbanded) || (do_hbanded && !do_sub)) /* ha! */
    FreeCP9Map(cp9map);
  if(do_hbanded && !do_sub)
    FreeCP9Bands(cp9b);

  if(do_hbanded || do_sub || do_ptest)
    {
      FreeCPlan9Matrix(cp9_mx);
      FreeCPlan9(orig_hmm);
    }
  if (do_qdb)
    {
      FreeBandDensities(cm, gamma);
      free(dmin);
      free(dmax);
      free(orig_dmin);
      free(orig_dmax);
    }
  StopwatchFree(watch1);
  StopwatchFree(watch2);
  
  *ret_tr = tr; 
  if (ret_postcode != NULL) *ret_postcode = postcode; 
  
  if(ret_post_spos != NULL)
    *ret_post_spos = post_spos;
  if(ret_post_epos != NULL)
    *ret_post_epos = post_epos;
  if(ret_dist_spos != NULL)
    *ret_dist_spos = dist_spos;
  if(ret_dist_epos != NULL)
    *ret_dist_epos = dist_epos;
}

/*
 * Function: serial_search_database
 * Date:     EPN, Thu Dec  7 18:07:24 2006
 *           original code by RJK, Tue May 28, 2002 [St. Louis]
 * Purpose:  Given an open database file, a model, and various parameters, does
 *           the search using CYKScan and then determines and prints out the 
 *           alignments.
 *
 * Parameters:        dbfp         the database
 *                    cm           the model
 *                    cons         precalc'ed consensus info for display
 *                    W            maximum size of hit
 *                    cutoff_type  either E_CUTOFF or SCORE_CUTOFF
 *                    cutoff       min. score to report
 *                    do_revcomp   search complementary strand
 *                    do_align     calculate and do alignment
 *                    do_stats     calculate statistics
 *                    mu           for stats
 *                    lambda       for stats
 *                    dmin         minimum bound on d for state v; 0..M (NULL if non-banded)
 *                    dmax         maximum bound on d for state v; 0..M (NULL if non-banded)         
 *                    do_inside    TRUE to scan with inside instead of CYK
 */
void serial_search_database (ESL_SQFILE *dbfp, CM_t *cm, CMConsensus_t *cons,
			     int W, int cutoff_type, float cutoff, 
			     int do_revcomp, int do_align, int do_stats,
			     double *mu, double *lambda, int *dmin, int *dmax,
			     int do_inside) 
{

  int reversed;                /* Am I currently doing reverse complement? */
  int i,a;
  db_seq_t *dbseq;
  float min_cutoff;
  /*printf("in serial_search database do_align: %d do_revcomp: %d\n", do_align, do_revcomp);*/
  
  if (cutoff_type == SCORE_CUTOFF) 
    min_cutoff = cutoff;
  else 
    min_cutoff = e_to_score (cutoff, mu, lambda);
  
  while ((dbseq = read_next_seq(dbfp, do_revcomp)))
    {
      for (reversed = 0; reversed <= do_revcomp; reversed++) 
	{
	  /* Scan */
	  dbseq->results[reversed] = CreateResults(INIT_RESULTS);
	  if(dmin == NULL && dmax == NULL)
	    if(do_inside)
	      iInsideScan(cm, dbseq->sq[reversed]->dsq, 1, dbseq->sq[reversed]->n, W,
			  min_cutoff, 0, dbseq->results[reversed]);
	    else /* !do_inside */
	      CYKScan (cm, dbseq->sq[reversed]->dsq, 1, dbseq->sq[reversed]->n, W,
			min_cutoff, 0, dbseq->results[reversed]);
	  else
	    if(do_inside)
	      iInsideBandedScan(cm, dbseq->sq[reversed]->dsq, dmin, dmax, 
			       1, dbseq->sq[reversed]->n, W,
			       min_cutoff, 0, dbseq->results[reversed]);
	    else /* !do_inside */
	    CYKBandedScan (cm, dbseq->sq[reversed]->dsq, dmin, dmax, 1, 
			   dbseq->sq[reversed]->n, W, min_cutoff, 0, 
			   dbseq->results[reversed]);
	  

	  /* Align results */
	  if (do_align) {
	    for (i=0; i<dbseq->results[reversed]->num_results; i++) {
	      CYKDivideAndConquer
		(cm, dbseq->sq[reversed]->dsq, dbseq->sq[reversed]->n,
		 dbseq->results[reversed]->data[i].bestr,
		 dbseq->results[reversed]->data[i].start, 
		 dbseq->results[reversed]->data[i].stop, 
		 &(dbseq->results[reversed]->data[i].tr),
		 dmin, dmax); /* dmin and dmax will be NULL if non-banded,
			       * alternatively, could always pass NULL to 
			       * always do non-banded alignment. */
	      
	      /* Now, subtract out the starting point of the result so 
		 that it can be added in later.  This makes the print_alignment
		 routine compatible with the parallel version, while not needing
		 to send the entire database seq over for each alignment. */
	      for (a=0; a<dbseq->results[reversed]->data[i].tr->n; a++) {
		dbseq->results[reversed]->data[i].tr->emitl[a] -= 
		  (dbseq->results[reversed]->data[i].start - 1);
		dbseq->results[reversed]->data[i].tr->emitr[a] -= 
		  (dbseq->results[reversed]->data[i].start - 1);
	      }
	    }
	  }
	}
      /* Print results */
      print_results (cm, cons, dbseq, do_revcomp, do_stats,
		     mu, lambda);
      fflush (stdout);
      
      FreeResults(dbseq->results[0]);
      esl_sq_Destroy(dbseq->sq[0]);
      if (do_revcomp) {
	esl_sq_Destroy(dbseq->sq[1]);
	FreeResults(dbseq->results[1]);
      }
    }
}

#ifdef USE_MPI
/*
 * Function: parallel_search_database 
 * Date:     original (RSEARCH): RJK, Tue May 28, 2002 [St. Louis] 
 *           infernalized: EPN, Sat Dec  9 12:57:42 2006
 * Purpose:  Given the same parameters as serial_search_database, does
 *           the database search with alignments and printing results, but
 *           in a parallel fashion.
 *     
 *           It the master node performs the following tasks:
 *           while (!eof(dbfp) && slaves_working) {
 *             For each empty process, send next job
 *             Wait for a result
 *             For each empty process, send next job
 *             Process result
 *           }
 *            
 *
 *           The slave processes do the following:
 *           while (1) {
 *             receive_job
 *             if (terminate code) return;
 *             do the job
 *             send the results
 *           }
 */
void parallel_search_database (ESL_SQFILE *dbfp, CM_t *cm, CMConsensus_t *cons,
			       int W, int cutoff_type, float cutoff, 
			       int do_revcomp, int do_align, int do_stats,
			       double *mu, double *lambda, int *dmin, int *dmax,
			       int do_inside, 
			       int mpi_my_rank, int mpi_master_rank, 
			       int mpi_num_procs) 
{
  char job_type;
  int seqlen;
  char *seq;
  scan_results_t *results;
  db_seq_t **active_seqs;
  job_t **process_status;
  int eof = FALSE;
  job_t *job_queue = NULL;
  int proc_index, active_seq_index;
  int bestr;       /* Best root state -- for alignments */
  Parsetree_t *tr;
  float min_cutoff;
  /*do_align = FALSE;
    dmin = NULL;
    dmax = NULL;*/
  /*printf("B PSD rank: %4d mast: %4d\n", mpi_my_rank, mpi_master_rank);*/

  if (cutoff_type == SCORE_CUTOFF) min_cutoff = cutoff;
  else min_cutoff = e_to_score (cutoff, mu, lambda);

  if (mpi_my_rank == mpi_master_rank) 
    {
      /* Set up arrays to hold pointers to active seqs and jobs on
	 processes */
      active_seqs = MallocOrDie(sizeof(db_seq_t *) * mpi_num_procs);
      process_status = MallocOrDie(sizeof(job_t *) * mpi_num_procs);
      for (active_seq_index=0; active_seq_index<mpi_num_procs; 
	   active_seq_index++) 
	active_seqs[active_seq_index] = NULL;
      for (proc_index = 0; proc_index < mpi_num_procs; proc_index++)
	process_status[proc_index] = NULL;
      
      do
	{
	  /* Check for idle processes.  Send jobs */
	  for (proc_index=0; proc_index<mpi_num_procs; proc_index++) {
	    if (proc_index == mpi_master_rank) continue;  /* Skip master process */
	    if (process_status[proc_index] == NULL) {         
	      /* I'm idle -- need a job */
	      if (job_queue == NULL) {           /* Queue is empty */
		/* Find next non-master open process */
		for (active_seq_index=0; active_seq_index<mpi_num_procs; 
		     active_seq_index++) {
		  if (active_seqs[active_seq_index] == NULL) break;
		}
		if (active_seq_index == mpi_num_procs) {
		  Die ("Tried to read more than %d seqs at once\n", mpi_num_procs);
		}
		active_seqs[active_seq_index] = read_next_seq(dbfp, do_revcomp);
		if (active_seqs[active_seq_index] == NULL) {
		  eof = TRUE;
		  break;            /* Queue is empty and no more seqs */
		}
		else
		  job_queue = enqueue (active_seqs[active_seq_index], 
				       active_seq_index, W, do_revcomp, 
				       STD_SCAN_WORK);
	      }
	      if (job_queue != NULL)
		send_next_job (&job_queue, process_status + proc_index, 
			       proc_index);
	    } 
	  } 
	  /* Wait for next reply */
	  if (procs_working(process_status, mpi_num_procs, mpi_master_rank)) 
	    {
	      active_seq_index = check_results (active_seqs, process_status, W);
	      if (active_seqs[active_seq_index]->chunks_sent == 0) 
		{
		  if (cutoff_type == E_CUTOFF)
		    remove_hits_over_e_cutoff 
		      (active_seqs[active_seq_index]->results[0],
		       active_seqs[active_seq_index]->sq[0]->seq,
		       cutoff, lambda, mu);
		  if (do_revcomp) 
		    {
		      if (cutoff_type == E_CUTOFF)
			remove_hits_over_e_cutoff 
			  (active_seqs[active_seq_index]->results[1],
			   active_seqs[active_seq_index]->sq[1]->seq,
			   cutoff, lambda, mu);
		    }
	      /* Check here if doing alignments and queue them or check if 
		 all done */
		  if (do_align && 
		      active_seqs[active_seq_index]->alignments_sent == -1) {
		    enqueue_alignments (&job_queue, active_seqs[active_seq_index],
					active_seq_index, do_revcomp, ALIGN_WORK);
		  }
		  if (!do_align || 
		      active_seqs[active_seq_index]->alignments_sent == 0) {
		    print_results (cm, cons, active_seqs[active_seq_index], 
				   do_revcomp, do_stats, mu, lambda);
		    if (do_revcomp) {
		      FreeResults(active_seqs[active_seq_index]->results[1]);
		      esl_sq_Destroy(active_seqs[active_seq_index]->sq[1]);
		    }
		    FreeResults(active_seqs[active_seq_index]->results[0]);
		    esl_sq_Destroy(active_seqs[active_seq_index]->sq[0]);
		    active_seqs[active_seq_index] = NULL;
		  }
		}
	    }
	} while (!eof || job_queue != NULL || 
		 (procs_working(process_status, mpi_num_procs, mpi_master_rank)));
      for (proc_index=0; proc_index<mpi_num_procs; proc_index++) {
	if (proc_index != mpi_master_rank) {
	  send_terminate (proc_index);
	}
      }
      free(active_seqs);
      free(process_status);
    } 
  else  /* not the master node */
    {
      seq = NULL;
      do 
	{
	  job_type = receive_job (&seqlen, &seq, &bestr, mpi_master_rank);
	  if (job_type == STD_SCAN_WORK) 
	    {
	      /* Do the scan */
	      results = CreateResults(INIT_RESULTS);
	      if(dmin == NULL && dmax == NULL)
		if(do_inside)
		  iInsideScan(cm, seq, 1, seqlen, W,
			     min_cutoff, 0, results);
		else /* !do_inside */
		  CYKScan (cm, seq, 1, seqlen, W,
			   min_cutoff, 0, results);
	      else
		if(do_inside)
		  iInsideBandedScan(cm, seq, dmin, dmax, 1, seqlen, W,
				   min_cutoff, 0, results);
		else /* !do_inside */
		  CYKBandedScan (cm, seq, dmin, dmax, 1, seqlen, W,
				 min_cutoff, 0, results);

	      send_scan_results (results, mpi_master_rank);
	      FreeResults(results);
	    } 
	  else if (job_type == ALIGN_WORK && do_align) 
	    {
	      CYKDivideAndConquer(cm, seq, seqlen, bestr, 1, seqlen, &tr,
				  dmin, dmax);
	      send_align_results (tr, mpi_master_rank);
	      FreeParsetree(tr);
	    }
	  if (seq != NULL)
	    free(seq);
	  seq = NULL;
	} while (job_type != TERMINATE_WORK);
    }
  MPI_Barrier(MPI_COMM_WORLD);
  /*printf("E PSD rank: %4d mast: %4d\n", mpi_my_rank, mpi_master_rank);*/
}
#endif


/*
 * Function: read_next_seq
 * Date:     RJK, Wed May 29, 2002 [St. Louis]
 *           easeled: EPN, Fri Dec  8 11:40:20 2006
 * Purpose:  Given a dbfp and whether or not to take the reverse complement,
 *           reads in the sequence and prepares reverse complement.
 */
db_seq_t *read_next_seq (ESL_SQFILE *dbfp, int do_revcomp) 
{
  db_seq_t *ret_dbseq;
  int status;
  char *tmp_seq;

  ret_dbseq = MallocOrDie(sizeof(db_seq_t));

  ret_dbseq->sq[0] = esl_sq_Create();
  status = (esl_sqio_Read(dbfp, ret_dbseq->sq[0]) == eslOK);

  while(status && ret_dbseq->sq[0]->n == 0) /* skip zero length seqs */
    {
      esl_sq_Reuse(ret_dbseq->sq[0]);
      status = (esl_sqio_Read(dbfp, ret_dbseq->sq[0]) == eslOK);
    }
  if(!status)
    return NULL;

  s2upper(ret_dbseq->sq[0]->seq);                      /* Makes uppercase */
  /* Following line will be unnecessary once ESL_SQ objects have dsq's implemented
   * (i.e. allocated and filled within a esl_sqio_Read() call */
  ret_dbseq->sq[0]->dsq = DigitizeSequence (ret_dbseq->sq[0]->seq, ret_dbseq->sq[0]->n);

  if (do_revcomp)
    {
      /* make a new SQ object, to store the reverse complement */
      tmp_seq = MallocOrDie(sizeof(char) * (ret_dbseq->sq[0]->n+1));
      revcomp(tmp_seq, ret_dbseq->sq[0]->seq);
      ret_dbseq->sq[1] = esl_sq_CreateFrom(ret_dbseq->sq[0]->name, tmp_seq, 
					ret_dbseq->sq[0]->desc, ret_dbseq->sq[0]->acc, 
					ret_dbseq->sq[0]->ss);
      /* Following line will be unnecessary once ESL_SQ objects have dsq's implemented
       * (i.e. allocated and filled within a esl_sq_CreateFrom() call */
      ret_dbseq->sq[1]->dsq = DigitizeSequence (ret_dbseq->sq[1]->seq, ret_dbseq->sq[1]->n);
    }
  ret_dbseq->results[0] = NULL;
  ret_dbseq->results[1] = NULL;

  return(ret_dbseq);
}

/*
 * Function: print_results
 * Date:     RJK, Wed May 29, 2002 [St. Louis]
 *           easelfied: EPN, Fri Dec  8 08:29:05 2006 
 * Purpose:  Given the needed information, prints the results.
 *
 *           cm                  the model
 *           cons                consensus seq for model (query seq)
 *           dbseq               the database seq
 *           name                sequence name
 *           len                 length of the sequence
 *           in_revcomp       are we doing the minus strand
 *           do_stats            should we calculate stats?
 *           mu, lambda          for statistics
 */
void print_results (CM_t *cm, CMConsensus_t *cons, db_seq_t *dbseq,
		    int do_complement, int do_stats, double *mu, 
		    double *lambda) 
{
  int i;
  char *name;
  int len;
  scan_results_t *results;
  Fancyali_t *ali;
  int in_revcomp;
  int header_printed = 0;
  int gc_comp;

  name = dbseq->sq[0]->name;
  len = dbseq->sq[0]->n;

  for (in_revcomp = 0; in_revcomp <= do_complement; in_revcomp++) {
    results = dbseq->results[in_revcomp];
    if (results == NULL || results->num_results == 0) continue;

    if (!header_printed) {
      header_printed = 1;
      printf (">%s\n\n", name);
    }
    printf ("  %s strand results:\n\n", in_revcomp ? "Minus" : "Plus");

    for (i=0; i<results->num_results; i++) {
      gc_comp = get_gc_comp (dbseq->sq[in_revcomp]->seq, 
			     results->data[i].start, results->data[i].stop);
      printf (" Query = %d - %d, Target = %d - %d\n", 
	      cons->lpos[cm->ndidx[results->data[i].bestr]]+1,
	      cons->rpos[cm->ndidx[results->data[i].bestr]]+1,
	      coordinate(in_revcomp, results->data[i].start, len), 
	      coordinate(in_revcomp, results->data[i].stop, len));
      if (do_stats) {
	printf (" Score = %.2f, E = %.4g, P = %.4g, GC = %3d\n", results->data[i].score,
		RJK_ExtremeValueE(results->data[i].score, mu[gc_comp], 
				  lambda[gc_comp]),
		esl_gumbel_surv((double) results->data[i].score, mu[gc_comp], 
				lambda[gc_comp]), gc_comp);
	/*printf("  Mu[gc=%d]: %f, Lambda[gc=%d]: %f\n", gc_comp, mu[gc_comp], gc_comp,
	  lambda[gc_comp]);*/
	/*ExtremeValueP(results->data[i].score, mu[gc_comp], 
	  lambda[gc_comp]));*/
      } else {
	printf (" Score = %.2f, GC = %3d\n", results->data[i].score, gc_comp);
      }
      printf ("\n");
      if (results->data[i].tr != NULL) {
	ali = CreateFancyAli (results->data[i].tr, cm, cons, 
			      dbseq->sq[in_revcomp]->dsq +
			      (results->data[i].start-1));
	PrintFancyAli(stdout, ali);
	printf ("\n");
      }
    }
  }
  fflush(stdout);
}

/* Function: get_gc_comp
 * Date:     RJK, Mon Oct 7, 2002 [St. Louis]
 * Purpose:  Given a sequence and start and stop coordinates, returns 
 *           integer GC composition of the region 
 */
int get_gc_comp(char *seq, int start, int stop) {
  int i;
  int gc_ct;
  char c;

  if (start > stop) {
    i = start;
    start = stop;
    stop = i;
  }
  gc_ct = 0;
  for (i=start; i<=stop; i++) {
    c = resolve_degenerate(seq[i]);
    if (c=='G' || c == 'C')
      gc_ct++;
  }
  return ((int)(100.*gc_ct/(stop-start+1)));
}
 
/*
 * Function: remove_hits_over_e_cutoff
 * Date:     RJK, Tue Oct 8, 2002 [St. Louis]
 * Purpose:  Given an E-value cutoff, lambdas, mus, a sequence, and
 *           a list of results, calculates GC content for each hit, 
 *           calculates E-value, and decides wheter to keep hit or not.
 */
void remove_hits_over_e_cutoff (scan_results_t *results, char *seq,
				float cutoff, double *lambda, double *mu) 
{
  int gc_comp;
  int i, x;
  scan_result_node_t swap;

  if (results == NULL)
    return;

  for (i=0; i<results->num_results; i++) {
    gc_comp = get_gc_comp (seq, results->data[i].start, results->data[i].stop);
    if (RJK_ExtremeValueE(results->data[i].score, 
			  mu[gc_comp], lambda[gc_comp])> cutoff) {
      results->data[i].start = -1;
    }
  }

  for (x=0; x<results->num_results; x++) {
    while (results->num_results > 0 && 
	   results->data[results->num_results-1].start == -1)
      results->num_results--;
    if (x<results->num_results && results->data[x].start == -1) {
      swap = results->data[x];
      results->data[x] = results->data[results->num_results-1];
      results->data[results->num_results-1] = swap;
      results->num_results--;
    }
  }
  while (results->num_results > 0 &&
	 results->data[results->num_results-1].start == -1)
    results->num_results--;
  sort_results(results);
}  

int  
EnforceSubsequence(CM_t *cm, int enf_start, char *enf_seq)
{
  int nd;
  float small_chance = 1e-15; /* any parse not including the enforced path includes
			       * an emission or transition with a -45 bit score */
  char *enf_dsq;
  int   enf_end;
  int v;
  int a;

  enf_end = enf_start + strlen(enf_seq) - 1;
  /*printf("in EnforceSubsequence, start posn: %d enf_seq: %s\n", enf_start, enf_seq);*/
  for(nd = (enf_start-1); nd <= enf_end; nd++)
    {
      if(cm->ndtype[nd] != MATL_nd)
	Die("ERROR, trying to enforce a non-MATL stretch (node: %d not MATL).\n", nd);
    }

  /* Go through each node and enforce the template by changing the emission and
   * transition probabilities as appropriate. */
  
  /* First deal with node before enf_start, we want to ensure that enf_start is
   * entered. We know enf_start - 1 and enf_start are both MATL nodes */
  nd = enf_start - 1;
  v  = cm->nodemap[nd];       /* MATL_ML*/
  cm->t[v][2] = small_chance; /* ML->D  */
  v++;                        /* MATL_D */
  cm->t[v][2] = small_chance; /*  D->D  */
  v++;                        /* MATL_IL*/
  cm->t[v][2] = small_chance; /*  IL->D */

  /* Now move on to the MATL nodes we're enforcing emits the enf_seq */
  enf_dsq = DigitizeSequence(enf_seq, (strlen(enf_seq)));
  for(a = 1; a <= strlen(enf_seq); a++)
    if(enf_dsq[a] > 3) 
      Die("ERROR enforced sequence must be contain only A,C,G,U.\n");

  for(nd = enf_start; nd <= enf_end; nd++) 
    {
      /* Enforce the transitions, unless we're the last node of the stretch */
      v  = cm->nodemap[nd];       /* MATL_ML*/
      if(nd < enf_end)
	{
	  cm->t[v][0] = small_chance; /* ML->IL */
	  cm->t[v][2] = small_chance; /* ML->D  */
	}
      /* Enforce the emission. */
      for(a = 0; a < MAXABET; a++)
	{
	  if(a != enf_dsq[(nd-enf_start+1)])
	    cm->e[v][a] = small_chance;
	  else
	    cm->e[v][a] = 1. - (3 * small_chance);
	}
    }
  CMRenormalize(cm);
  CMLogoddsify(cm);
  return 1;
}

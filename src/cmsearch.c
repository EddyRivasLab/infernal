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

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "squid.h"		/* general sequence analysis library    */
#include "msa.h"                /* squid's multiple alignment i/o       */
#include "stopwatch.h"          /* squid's process timing module        */

#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* external functions                   */
#include "hmmer_funcs.h"
#include "hmmer_structs.h"
#include "hmmband.h"

static void debug_print_bands(CM_t *cm, int *dmin, int *dmax);

static char banner[] = "cmsearch - search a sequence database with an RNA covariance model";

static char usage[]  = "\
Usage: cmsearch [-options] <cmfile> <sequence file>\n\
The sequence file is expected to be in FASTA format.\n\
  Most commonly used options are:\n\
   -h     : help; print brief help on version and usage\n\
   -W <n> : set scanning window size to <n> (default: precalc'd in cmbuild)\n\
";

static char experts[] = "\
  Expert, in development, or infrequently used options are:\n\
   --informat <s>: specify that input alignment is in format <s>, not FASTA\n\
   --toponly     : only search the top strand\n\
   --local       : do local alignment\n\
   --noalign     : find start/stop only; don't do alignments\n\
   --dumptrees   : dump verbose parse tree information for each hit\n\
   --banded      : use experimental banded CYK scanning algorithm\n\
   --bandp <f>   : tail loss prob for --banded [df: 0.0001]\n\
   --banddump    : print bands for each state\n\
   --thresh <f>  : reporting bit score threshold (try 0 before < 0) [df: 0]\n\
   --X           : project X!\n\
   --inside      : scan with Inside, not CYK (caution ~5X slower(!))\n\
\n\
  * HMM banded search/alignment related options:\n\
   --hmmonly     : scan with a CM plan 9 HMM, not the CM\n\
   --hbanded     : use banded CYK with bands from CM plan 9 HMM scan\n\
   --hbandp <f>  : tail loss prob for --hbanded (default:0.0001)\n\
   --sums        : use posterior sums during HMM band calculation (widens bands)\n\
";

static struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE }, 
  { "-W", TRUE, sqdARG_INT }, 
  { "--dumptrees",  FALSE, sqdARG_NONE },
  { "--informat",   FALSE, sqdARG_STRING },
  { "--local",      FALSE, sqdARG_NONE },
  { "--noalign",    FALSE, sqdARG_NONE },
  { "--toponly",    FALSE, sqdARG_NONE },
  { "--banded",     FALSE, sqdARG_NONE },
  { "--bandp",      FALSE, sqdARG_FLOAT},
  { "--banddump",   FALSE, sqdARG_NONE},
  { "--thresh",     FALSE, sqdARG_FLOAT},
  { "--X",          FALSE, sqdARG_NONE },
  { "--hmmonly",    FALSE, sqdARG_NONE },
  { "--inside",     FALSE, sqdARG_NONE },
  { "--hbanded",    FALSE, sqdARG_NONE },
  { "--hbandp",     FALSE, sqdARG_FLOAT},
  { "--sums",       FALSE, sqdARG_NONE},
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv)
{
  char            *cmfile;      /* file to read CM from */	
  char            *seqfile;     /* file to read sequences from */
  int              format;      /* format of sequence file */
  CMFILE          *cmfp;        /* open CM file for reading */
  SQFILE	  *sqfp;        /* open seqfile for reading */
  CM_t            *cm;          /* a covariance model       */
  char            *seq;         /* RNA sequence */
  SQINFO           sqinfo;      /* optional info attached to seq */
  char            *dsq;         /* digitized RNA sequence */
  Stopwatch_t     *watch;
  int              i;
  int              reversed;    /* TRUE when we're doing the reverse complement strand */
  int              maxlen;
  Parsetree_t     *tr;		/* parse of an individual hit */
  CMConsensus_t   *cons;	/* precalculated consensus info for display */
  Fancyali_t      *ali;         /* alignment, formatted for display */

  double  **gamma;              /* P(subseq length = n) for each state v    */
  int     *dmin;		/* minimum d bound for state v, [0..v..M-1] */
  int     *dmax; 		/* maximum d bound for state v, [0..v..M-1] */
  double   bandp;		/* tail loss probability for banding */

  int    nhits;			/* number of hits in a seq */
  int   *hitr;			/* initial states for hits */
  int   *hiti;                  /* start positions of hits */
  int   *hitj;                  /* end positions of hits */
  float *hitsc;			/* scores of hits */

  int    windowlen;		/* maximum len of hit; scanning window size */
  int    do_revcomp;		/* true to do reverse complement too */
  int    do_local;		/* TRUE to do local alignment */
  int    do_align;              /* TRUE to calculate and show alignments */
  int    do_dumptrees;		/* TRUE to dump parse trees */
  int    do_banded;		/* TRUE to do banded CYK */
  int    do_projectx;           /* TRUE to activate special in-progress testing code */
  int    do_bdump;              /* TRUE to print out bands */
  /*EPN 08.18.05*/
  int    set_window;            /* TRUE to set window length due to -W option*/

  char *optname;                /* name of option found by Getopt()        */
  char *optarg;                 /* argument found by Getopt()              */
  int   optind;                 /* index in argv[]                         */


  /*EPN 11.11.05 */
  int      safe_windowlen;	/* initial windowlen (W) used for calculating bands
				 * in BandCalculationEngine().
				 * this needs to be significantly bigger than what
				 * we expect dmax[0] to be, for truncation error
				 * handling.
				 */

  float thresh;                 /* bit score threshold, report scores > thresh */
  int   use_hmm;                /* TRUE to scan with a CM Plan 9 HMM */
  int   do_hmmonly;             /* TRUE to scan with a CM Plan 9 HMM ONLY!*/

  /* CM Plan 9 HMM data structures */
  struct cplan9_s       *cp9_hmm;       /* constructed CP9 HMM; written to hmmfile              */
  int                    cp9_M;         /* number of nodes in CP9 HMM (MATL+MATR+2*MATP)        */
  struct cp9_dpmatrix_s *cp9_mx;        /* growable DP matrix for viterbi                       */
  struct cp9_dpmatrix_s *cp9_fwd;       /* growable DP matrix for forward                       */
  struct cp9_dpmatrix_s *cp9_bck;       /* growable DP matrix for backward                      */
  struct cp9_dpmatrix_s *cp9_posterior; /* growable DP matrix for posterior decode              */
  float                  fwd_sc;        /* score for Forward() run */
  float                  fb_sc;         /* score from Forward or Backward */
  float                  sc;            /* score from CYK */

  /* data structures for hmm bands (bands on the hmm states) */
  int     *pn_min_m;          /* HMM band: minimum position node k match state will emit  */
  int     *pn_max_m;          /* HMM band: maximum position node k match state will emit  */
  int     *pn_min_i;          /* HMM band: minimum position node k insert state will emit */
  int     *pn_max_i;          /* HMM band: maximum position node k insert state will emit */
  int     *pn_min_d;          /* HMM band: minimum position node k delete state will emit */
  int     *pn_max_d;          /* HMM band: maximum position node k delete state will emit */
  double   hbandp;            /* tail loss probability for hmm bands */
  int      use_sums;          /* TRUE to fill and use the posterior sums, false not to. */
  int     *isum_pn_m;         /* [1..k..M] sum over i of log post probs from post->mmx[i][k]*/
  int     *isum_pn_i;         /* [1..k..M] sum over i of log post probs from post->imx[i][k]*/
  int     *isum_pn_d;         /* [1..k..M] sum over i of log post probs from post->dmx[i][k]*/

  /* data structures for mapping HMM to CM */
  int *node_cc_left;    /* consensus column each CM node's left emission maps to
			 * [0..(cm->nodes-1)], -1 if maps to no consensus column*/
  int *node_cc_right;   /* consensus column each CM node's right emission maps to
			 * [0..(cm->nodes-1)], -1 if maps to no consensus column*/
  int *cc_node_map;     /* node that each consensus column maps to (is modelled by)
			 * [1..hmm_nmc] */
  int **cs2hn_map;      /* 2D CM state to HMM node map, 1st D - CM state index
		         * 2nd D - 0 or 1 (up to 2 matching HMM states), value: HMM node
		         * that contains state that maps to CM state, -1 if none.*/
  int **cs2hs_map;      /* 2D CM state to HMM node map, 1st D - CM state index
		         * 2nd D - 2 elements for up to 2 matching HMM states, 
		         * value: HMM STATE (0(M), 1(I), 2(D) that maps to CM state,
		         * -1 if none.
		         * For example: HMM node cs2hn_map[v][0], state cs2hs_map[v][0]
                         * maps to CM state v.*/
  int ***hns2cs_map;    /* 3D HMM node-state to CM state map, 1st D - HMM node index, 2nd D - 
		         * HMM state (0(M), 1(I), 2(D)), 3rd D - 2 elements for up to 
		         * 2 matching CM states, value: CM states that map, -1 if none.
		         * For example: CM states hsn2cs_map[k][0][0] and hsn2cs_map[k][0][1]
		         * map to HMM node k's match state.*/

  /* arrays for CM state bands, derived from HMM bands */
  int *imin;            /* [1..M] imin[v] = first position in band on i for state v*/
  int *imax;            /* [1..M] imax[v] = last position in band on i for state v*/
  int *jmin;            /* [1..M] jmin[v] = first position in band on j for state v*/
  int *jmax;            /* [1..M] jmax[v] = last position in band on j for state v*/
  int **hdmin;          /* [v=1..M][0..(jmax[v]-jmin[v])] 
			 * hdmin[v][j0] = first position in band on d for state v, and position
			 * j = jmin[v] + j0.*/
  int **hdmax;          /* [v=1..M][0..(jmax[v]-jmin[v])] 
			 * hdmin[v][j0] = last position in band on d for state v, and position
			 * j = jmin[v] + j0.*/
  int *safe_hdmin;      /* [1..M] safe_hdmin[v] = min_d (hdmin[v][j0]) (over all valid j0) */
  int *safe_hdmax;      /* [1..M] safe_hdmax[v] = max_d (hdmax[v][j0]) (over all valid j0) */

  unsigned char   *p7dsq;     /* digitized RNA sequences (plan 7 version)*/
  int    ks;            /* Counter over HMM state types (0 (match), 1(ins) or 2 (del))*/
  int    v;             /* counter over states of the CM */
  int    k;             /* counter over HMM nodes */
  int    debug_level;   /* verbosity level for debugging printf() statements,
			 * passed to many functions. */
  float hmm_thresh;     /* bit score threshold for reporting hits to HMM */

  int do_inside;        /* TRUE to use scanning Inside algorithm instead of CYK */
  int do_hbanded;       /* TRUE to first scan with a CP9 HMM to derive bands for a CYK scan */
  /*********************************************** 
   * Parse command line
   ***********************************************/

  format            = SQFILE_UNKNOWN;
  windowlen         = 200;
  set_window        = FALSE;
  do_revcomp        = TRUE;
  do_local          = FALSE;
  do_align          = TRUE;
  do_dumptrees      = FALSE;
  do_banded         = FALSE;
  bandp             = 0.0001;
  do_projectx       = FALSE;
  do_bdump          = FALSE;
  thresh            = 0.;
  do_hmmonly        = FALSE;
  use_hmm           = FALSE;
  do_inside         = FALSE;
  hmm_thresh        = 0;
  do_hbanded        = TRUE;
  hbandp            = 0.0001;
  use_sums          = FALSE;
  debug_level = 0;
  
  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if       (strcmp(optname, "-W")          == 0) 
      { windowlen    = atoi(optarg); set_window = TRUE; }
    else if  (strcmp(optname, "--dumptrees") == 0) do_dumptrees = TRUE;
    else if  (strcmp(optname, "--local")     == 0) do_local     = TRUE;
    else if  (strcmp(optname, "--noalign")   == 0) do_align     = FALSE;
    else if  (strcmp(optname, "--toponly")   == 0) do_revcomp   = FALSE;
    else if  (strcmp(optname, "--banded")    == 0) do_banded    = TRUE;
    else if  (strcmp(optname, "--bandp")     == 0) bandp        = atof(optarg);
    else if  (strcmp(optname, "--X")         == 0) do_projectx  = TRUE;
    else if  (strcmp(optname, "--banddump")  == 0) do_bdump     = TRUE;
    else if  (strcmp(optname, "--thresh")    == 0) thresh       = atof(optarg);
    else if  (strcmp(optname, "--hmmonly")   == 0) { use_hmm      = TRUE; do_hmmonly = TRUE; } 
    else if  (strcmp(optname, "--inside")    == 0) do_inside    = TRUE;
    else if  (strcmp(optname, "--hbanded")   == 0) {use_hmm = TRUE; do_hbanded  = TRUE;}
    else if  (strcmp(optname, "--hbandp")    == 0) hbandp       = atof(optarg);
    else if  (strcmp(optname, "--sums")      == 0) use_sums     = TRUE;
    else if  (strcmp(optname, "--informat")  == 0) {
      format = String2SeqfileFormat(optarg);
      if (format == SQFILE_UNKNOWN) 
	Die("unrecognized sequence file format \"%s\"", optarg);
    }
    else if (strcmp(optname, "-h") == 0) {
      MainBanner(stdout, banner);
      puts(usage);
      puts(experts);
      exit(EXIT_SUCCESS);
    }
  }

  if (argc - optind != 2) Die("Incorrect number of arguments.\n%s\n", usage);
  cmfile = argv[optind++];
  seqfile = argv[optind++]; 

  /*********************************************** 
   * Preliminaries: open our files for i/o; get a CM
   ***********************************************/

  watch = StopwatchCreate();

  if ((sqfp = SeqfileOpen(seqfile, format, NULL)) == NULL)
    Die("Failed to open sequence database file %s\n%s\n", seqfile, usage);
  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL)
    Die("Failed to open covariance model save file %s\n%s\n", cmfile, usage);

  if (! CMFileRead(cmfp, &cm))
    Die("Failed to read a CM from %s -- file corrupt?\n", cmfile);
  if (cm == NULL) 
    Die("%s empty?\n", cmfile);

  /* EPN 08.18.05 */
  if (! (set_window)) windowlen = cm->W;
  /*printf("***cm->W : %d***\n", cm->W);*/

  if (do_local) ConfigLocal(cm, 0.5, 0.5);
  CMLogoddsify(cm);
  /*CMHackInsertScores(cm);*/	/* make insert emissions score zero. "TEMPORARY" FIX. */
  cons = CreateCMConsensus(cm, 3.0, 1.0); 

  if (do_banded || do_projectx || do_bdump)
    {
      safe_windowlen = windowlen * 2;
      while(!(BandCalculationEngine(cm, safe_windowlen, bandp, 0, &dmin, &dmax, &gamma, do_local)))
	{
	  FreeBandDensities(cm, gamma);
	  free(dmin);
	  free(dmax);
	  safe_windowlen *= 2;
	}

      /* EPN 11.11.05 
       * An important design decision.
       * We're changing the windowlen value here. By default,
       * windowlen is read from the cm file (set to cm->W). 
       * Here we're doing a banded scan though. Its pointless to allow
       * a windowlen that's greater than the largest possible banded hit 
       * (which is dmax[0]). So we reset windowlen to dmax[0].
       * Its also possible that BandCalculationEngine() returns a dmax[0] that 
       * is > cm->W. This should only happen if the bandp we're using now is < 1E-7 
       * (1E-7 is the bandp value used to determine cm->W in cmbuild). If this 
       * happens, the current implementation reassigns windowlen to this larger value.
       * NOTE: if W was set at the command line, the command line value is 
       *       always used.
       */
      if(!(set_window))
	{
	  windowlen = dmax[0];
	}
      if(do_bdump) 
	{
	  printf("bandp:%f\n", bandp);
	  debug_print_bands(cm, dmin, dmax);
	}
    }
  
  if (use_hmm)
    {
      /* build a CM Plan 9 HMM, and use it to scan. */
      Alphabet_type = hmmNOTSETYET;
      SetAlphabet(hmmNUCLEIC); /* Set up the hmmer_alphabet global variable */
      cp9_M = 0;
      for(v = 0; v <= cm->M; v++)
	{
	  if(cm->stid[v] ==  MATP_MP)
	    cp9_M += 2;
	  else if(cm->stid[v] == MATL_ML || cm->stid[v] == MATR_MR)
	    cp9_M ++;
	}
      /* build the HMM data structure */
      cp9_hmm = AllocCPlan9(cp9_M);
      ZeroCPlan9(cp9_hmm);
      
      /* Get information mapping the HMM to the CM and vice versa, used
       * for mapping bands. 
       * This is relevant if we're building (read_cp9_flag == FALSE)
       * or reading from a file (read_cp9_flag == TRUE)
       */
      node_cc_left  = malloc(sizeof(int) * cm->nodes);
      node_cc_right = malloc(sizeof(int) * cm->nodes);
      cc_node_map   = malloc(sizeof(int) * (cp9_hmm->M + 1));
      map_consensus_columns(cm, cp9_hmm->M, node_cc_left, node_cc_right,
			    cc_node_map, debug_level);
      
      cs2hn_map     = malloc(sizeof(int *) * (cm->M+1));
      for(v = 0; v <= cm->M; v++)
	cs2hn_map[v]     = malloc(sizeof(int) * 2);
      
      cs2hs_map     = malloc(sizeof(int *) * (cm->M+1));
      for(v = 0; v <= cm->M; v++)
	cs2hs_map[v]     = malloc(sizeof(int) * 2);
      
      hns2cs_map    = malloc(sizeof(int **) * (cp9_hmm->M+1));
      for(k = 0; k <= cp9_hmm->M; k++)
	{
	  hns2cs_map[k]    = malloc(sizeof(int *) * 3);
	  for(ks = 0; ks < 3; ks++)
	    hns2cs_map[k][ks]= malloc(sizeof(int) * 2);
	}
      
      CP9_map_cm2hmm_and_hmm2cm(cm, cp9_hmm, node_cc_left, node_cc_right, 
				cc_node_map, cs2hn_map, cs2hs_map, 
				hns2cs_map, debug_level);
      
      /* fill in parameters of HMM using the CM and some ideas/formulas/tricks
       * from Zasha Weinberg's thesis (~p.123) (see CP9_cm2wrhmm.c code) */
      if(!(CP9_cm2wrhmm(cm, cp9_hmm, node_cc_left, node_cc_right, cc_node_map, cs2hn_map,
			cs2hs_map, hns2cs_map, debug_level)))
	Die("Couldn't build a CM Plan 9 HMM from the CM.\n");

      /*debug_print_cp9_params(cp9_hmm);
	exit(1);*/

    }
  /* Allocate data structures for use with HMM banding strategy (CP9) */
  if(use_hmm)
    {
      pn_min_m    = malloc(sizeof(int) * (cp9_hmm->M+1));
      pn_max_m    = malloc(sizeof(int) * (cp9_hmm->M+1));
      pn_min_i    = malloc(sizeof(int) * (cp9_hmm->M+1));
      pn_max_i    = malloc(sizeof(int) * (cp9_hmm->M+1));
      pn_min_d    = malloc(sizeof(int) * (cp9_hmm->M+1));
      pn_max_d    = malloc(sizeof(int) * (cp9_hmm->M+1));
      imin        = malloc(sizeof(int) * cm->M);
      imax        = malloc(sizeof(int) * cm->M);
      jmin        = malloc(sizeof(int) * cm->M);
      jmax        = malloc(sizeof(int) * cm->M);
      hdmin       = malloc(sizeof(int *) * cm->M);
      hdmax       = malloc(sizeof(int *) * cm->M);
      safe_hdmin  = malloc(sizeof(int) * cm->M);
      safe_hdmax  = malloc(sizeof(int) * cm->M);
      isum_pn_m   = malloc(sizeof(int) * cp9_hmm->M+1);
      isum_pn_i   = malloc(sizeof(int) * cp9_hmm->M+1);
      isum_pn_d   = malloc(sizeof(int) * cp9_hmm->M+1);
    }      
  
  StopwatchZero(watch);
  StopwatchStart(watch);
  
  /* EPN 11.18.05 Now that know what windowlen is, we need to ensure that
   * cm->el_selfsc * W >= IMPOSSIBLE (cm->el_selfsc is the score for an EL self transition)
   * This is done because we are potentially multiply cm->el_selfsc * W, and adding
   * that to IMPOSSIBLE. To avoid underflow issues this value must be less than
   * 3 * IMPOSSIBLE. Here we guarantee its less than 2 * IMPOSSIBLE (to be safe).
   */
  if((cm->el_selfsc * windowlen) < IMPOSSIBLE)
    cm->el_selfsc = (IMPOSSIBLE / (windowlen+1));

  maxlen   = 0;
  reversed = FALSE;
  while (reversed || ReadSeq(sqfp, sqfp->format, &seq, &sqinfo))
    {
      if (sqinfo.len == 0) continue; 	/* silently skip len 0 seqs */
      if (sqinfo.len > maxlen) maxlen = sqinfo.len;
      dsq = DigitizeSequence(seq, sqinfo.len);

      hmm_thresh = thresh;
      if (use_hmm) /* either scan only with CP9 HMM, or use it to infer bands for CYK */
	{
	  p7dsq = hmmer_DigitizeSequence(seq, sqinfo.len);

	  if(do_hmmonly) /* scan only with CP9 HMM */
	    {
	      fwd_sc = CP9ForwardScan(p7dsq, sqinfo.len, windowlen, cp9_hmm, &cp9_fwd, 
				      &nhits, &hitr, &hiti, &hitj, &hitsc, hmm_thresh);
	      /*printf("forward  sc: %f\n", fwd_sc);*/
	      /*printf("hmmer_Alphabet: %s\n", hmmer_Alphabet);
		sc = CP9Forward(p7dsq, sqinfo.len, cp9_hmm, &cp9_fwd);*/
	    }
	  else /* use CP9 HMM to infer bands for CYK */
	    {
	      /* Scan the sequence with "scanning" Forward and Backward algorithms,
	       * Forward gives likely endpoints (j's), Backward gives likely start points
	       * (i's). Then combine them to get likely i and j pairs (see 
	       * code in CP9_scan.c. 
	       */
	      fb_sc = CP9ForwardBackwardScan(p7dsq, sqinfo.len, windowlen, cp9_hmm, &cp9_fwd, &cp9_bck,
					     &nhits, &hitr, &hiti, &hitj, &hitsc, hmm_thresh);
	      /*printf("forward/backward sc: %f\n", fb_sc);*/
	      FreeCPlan9Matrix(cp9_fwd);
	      FreeCPlan9Matrix(cp9_bck);

	      /* For each hit defined by an i-j pair, use non-scanning, traditional
	       * Forward and Backward algorithms to get hmm bands, we'll use to 
	       * get j and d bands (a la HMM banded cmalign).
	       * Step 1: Get HMM posteriors.
	       * Step 2: posteriors -> HMM bands.
	       * Step 3: HMM bands  ->  CM bands.
	       */
	      for (i = 0; i < nhits; i++)
		{
		  /* Step 1: Get HMM posteriors.*/
		  fb_sc = CP9Forward(p7dsq, hiti[i], hitj[i], cp9_hmm, &cp9_fwd);
		  /*printf("hit: %d i: %d j: %d forward_sc : %f\n", i, hiti[i], hitj[i], fb_sc);*/
		  fb_sc = CP9Backward(p7dsq, hiti[i], hitj[i], cp9_hmm, &cp9_bck);
		  /*printf("CP9 i: %d | backward_sc: %f\n", i, fb_sc);*/

		  /*debug_check_CP9_FB(cp9_fwd, cp9_bck, cp9_hmm, fb_sc, hiti[i], hitj[i], p7dsq);*/
		  cp9_posterior = cp9_bck;
		  CP9FullPosterior(p7dsq, hiti[i], hitj[i], cp9_hmm, cp9_fwd, cp9_bck, cp9_posterior);

		  /* Step 2: posteriors -> HMM bands.
		   * NOTE: HMM bands have offset sequence indices, from 1 to W, 
		   * with W = hitj[i] - hiti[i] + 1.
		   */
		  if(!(use_sums))
		    {
		      /* match states */
		      CP9_hmm_band_bounds(cp9_posterior->mmx, hiti[i], hitj[i], cp9_hmm->M,
					  NULL, pn_min_m, pn_max_m, (1.-hbandp), HMMMATCH, debug_level);
		      /* insert states */
		      CP9_hmm_band_bounds(cp9_posterior->imx, hiti[i], hitj[i], cp9_hmm->M,
					  NULL, pn_min_i, pn_max_i, (1.-hbandp), HMMINSERT, debug_level);
		      /* delete states (note: delete_flag set to TRUE) */
		      CP9_hmm_band_bounds(cp9_posterior->dmx, hiti[i], hitj[i], cp9_hmm->M,
					  NULL, pn_min_d, pn_max_d, (1.-hbandp), HMMDELETE, debug_level);
		    }
		  else
		    {
		      CP9_ifill_post_sums(cp9_posterior, hiti[i], hitj[i], cp9_hmm->M,
					  isum_pn_m, isum_pn_i, 
					  isum_pn_d);
		      /* match states */
		      CP9_hmm_band_bounds(cp9_posterior->mmx, hiti[i], hitj[i], cp9_hmm->M,
					  isum_pn_m, pn_min_m, pn_max_m, (1.-hbandp), HMMMATCH, debug_level);
		      /* insert states */
		      CP9_hmm_band_bounds(cp9_posterior->imx, hiti[i], hitj[i], cp9_hmm->M,
					  isum_pn_i, pn_min_i, pn_max_i, (1.-hbandp), HMMINSERT, debug_level);
		      /* delete states */
		      CP9_hmm_band_bounds(cp9_posterior->dmx, hiti[i], hitj[i], cp9_hmm->M,
					  isum_pn_d, pn_min_d, pn_max_d, (1.-hbandp), HMMDELETE, debug_level);
		    }
		  if(debug_level != 0)
		    {
		      printf("printing hmm bands\n");
		      print_hmm_bands(stdout, sqinfo.len, cp9_hmm->M, pn_min_m, pn_max_m, pn_min_i,
				      pn_max_i, pn_min_d, pn_max_d, hbandp, debug_level);
		    }
		  
		  /* Step 3: HMM bands  ->  CM bands. */
		  hmm2ij_bands(cm, cp9_hmm->M, node_cc_left, node_cc_right, cc_node_map, 
			       hiti[i], hitj[i], pn_min_m, pn_max_m, pn_min_i, pn_max_i, 
			       pn_min_d, pn_max_d, imin, imax, jmin, jmax, cs2hn_map,
			       debug_level);
	  
		  /* Use the CM bands on i and j to get bands on d, specific to j. */
		  for(v = 0; v < cm->M; v++)
		    {
		      hdmin[v] = malloc(sizeof(int) * (jmax[v] - jmin[v] + 1));
		      hdmax[v] = malloc(sizeof(int) * (jmax[v] - jmin[v] + 1));
		    }
		  ij2d_bands(cm, (hitj[i] - hiti[i] + 1), imin, imax, jmin, jmax,
			     hdmin, hdmax, -1);

		  /*if(debug_level != 0)*/
		  /*PrintDPCellsSaved_jd(cm, jmin, jmax, hdmin, hdmax, (hitj[i] - hiti[i] + 1));*/

		  /*debug_print_hd_bands(cm, hdmin, hdmax, jmin, jmax);*/

		  FreeCPlan9Matrix(cp9_fwd);
		  FreeCPlan9Matrix(cp9_bck);

		  if (do_align) 
		    {
		      printf("Aligning subseq from i: %d to j: %d using HMM bands\n", hiti[i], hitj[i]);
		      sc = CYKInside_b_jd(cm, dsq, sqinfo.len, 0, hiti[i], hitj[i], &tr, 
					  jmin, jmax, hdmin, hdmax, safe_hdmin, safe_hdmax);
		      /*if(bdump_level > 0)
			banded_trace_info_dump(cm, tr[i], safe_hdmin, safe_hdmax, bdump_level);
		      */
		      if(sc >= thresh) /* only print alignments with CYK scores > reporting thresh */
			{
			  printf("\tCYK i: %5d j: %5d sc: %10.2f bits\n", hiti[i], hitj[i], sc);
			  ali = CreateFancyAli(tr, cm, cons, dsq);
			  PrintFancyAli(stdout, ali);
			  FreeFancyAli(ali);
			  FreeParsetree(tr);
			}
		    }
		}
	    }
	}
      else if (do_banded)
	if (do_inside)
	  InsideBandedScan(cm, dsq, dmin, dmax, sqinfo.len, windowlen, 
			   &nhits, &hitr, &hiti, &hitj, &hitsc, thresh);
	else
	  CYKBandedScan(cm, dsq, dmin, dmax, sqinfo.len, windowlen, 
		      &nhits, &hitr, &hiti, &hitj, &hitsc, thresh);
      else if (do_inside)
	InsideScan(cm, dsq, sqinfo.len, windowlen, 
		   &nhits, &hitr, &hiti, &hitj, &hitsc, thresh);
      else
	CYKScan(cm, dsq, sqinfo.len, windowlen, 
		&nhits, &hitr, &hiti, &hitj, &hitsc, thresh);
      if (! reversed) printf("sequence: %s\n", sqinfo.name);

      if(!(use_hmm && !do_hmmonly))
	{
	  for (i = 0; i < nhits; i++)
	    {
	      printf("hit %-4d: %6d %6d %8.2f bits\n", i, 
		     reversed ? sqinfo.len - hiti[i] + 1 : hiti[i], 
		     reversed ? sqinfo.len - hitj[i] + 1 : hitj[i],
		     hitsc[i]);
	      if (do_align) 
		{
		  CYKDivideAndConquer(cm, dsq, sqinfo.len, 
				      hitr[i], hiti[i], hitj[i], &tr);
		  ali = CreateFancyAli(tr, cm, cons, dsq);
		  PrintFancyAli(stdout, ali);
		  
		  if (do_dumptrees) {
		    ParsetreeDump(stdout, tr, cm, dsq);
		    printf("\tscore = %.2f\n", ParsetreeScore(cm,tr,dsq));
		  }
		  if (do_projectx) {
		    BandedParsetreeDump(stdout, tr, cm, dsq, gamma, windowlen, dmin, dmax);
		  }
		  
		  FreeFancyAli(ali);
		  FreeParsetree(tr);
		}
	    }
	}
      free(hitr);
      free(hiti);
      free(hitj);
      free(hitsc);
      
      free(dsq);
      if(use_hmm) free(p7dsq);
      if (! reversed && do_revcomp) {
	revcomp(seq,seq);
	reversed = TRUE;
      } else {
	reversed = FALSE;
	FreeSequence(seq, &sqinfo);
      }
    }
  StopwatchStop(watch);
  StopwatchDisplay(stdout, "\nCPU time: ", watch);
  if(use_hmm) printf("CP9 Forward memory:   %8.2f MB\n", CP9ForwardScanRequires(cp9_hmm, maxlen, windowlen));
  printf("CYK memory        :   %8.2f MB\n\n", CYKScanRequires(cm, maxlen, windowlen));

  if (do_banded)
    {
      FreeBandDensities(cm, gamma);
      free(dmin);
      free(dmax);
    }
  FreeCMConsensus(cons);
  FreeCM(cm);
  CMFileClose(cmfp);
  SeqfileClose(sqfp);
  StopwatchFree(watch);
  SqdClean();
  return EXIT_SUCCESS;
}

/* EPN 05.09.05
  debug_print_bands()
 * Function: debug_print_bands
 *
 * Purpose:  Print bands for each state.
 */

static void
debug_print_bands(CM_t *cm, int *dmin, int *dmax)
{
  int v;

  printf("\nPrinting bands :\n");
  printf("****************\n");
  for(v = 0; v < cm->M; v++)
   {
     printf("band state:%d type:%d min:%d max:%d\n", v, cm->sttype[v], dmin[v], dmax[v]);
   }
  printf("****************\n\n");
}


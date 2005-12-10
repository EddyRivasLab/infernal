/* cmalign.c
 * SRE, Thu Jul 25 11:28:03 2002 [St. Louis]
 * CVS $Id$
 * 
 * Align sequences to a CM.
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************
 */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>

#include "squid.h"		/* general sequence analysis library    */
#include "msa.h"                /* squid's multiple alignment i/o       */
#include "stopwatch.h"          /* squid's process timing module        */

#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* external functions                   */

MSA *Parsetrees2Alignment(CM_t *cm, char **dsq, SQINFO *sqinfo, float *wgt, 
			  Parsetree_t **tr, int nseq);

MSA *Parsetrees2Alignment_full(CM_t *cm, char **dsq, SQINFO *sqinfo, float *wgt, 
			       Parsetree_t **tr, int nseq);

static void debug_print_bands(CM_t *cm, int *dmin, int *dmax);
static void ExpandBands(CM_t *cm, int qlen, int *dmin, int *dmax);
static void banded_trace_info_dump(CM_t *cm, Parsetree_t *tr, int *dmin, int *dmax, int bdump_level);

static char banner[] = "cmalign - align sequences to an RNA CM";

static char usage[]  = "\
Usage: cmalign [-options] <cmfile> <sequence file>\n\
  Most commonly used options are:\n\
   -h     : help; print brief help on version and usage\n\
   -l     : local; align locally w.r.t. the model\n\
   -o <f> : output the alignment file to file <f>\n\
   -q     : quiet; suppress verbose banner\n\
";

static char experts[] = "\
  Expert, in development, or infrequently used options are:\n\
   --informat <s>: specify that input alignment is in format <s>\n\
   --nosmall     : use normal alignment algorithm, not d&c\n\
   --regress <f> : save regression test data to file <f>\n\
   --banded      : use experimental banded CYK alignment algorithm\n\
   --bandp <f>   : tail loss prob for --banded (default:0.0001)\n\
   --bandexpand  : naively expand bands if target sequence is outside root band\n\
   --banddump <n>: turn band info print statements to verbosity level <n> [1-3]\n\
    -W <n>       : window size for calculating bands (default: precalc'd in cmbuild)\n\
   --full        : include all match columns in output alignment\n\
";

static struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE }, 
  { "-l", TRUE, sqdARG_NONE },
  { "-o", TRUE, sqdARG_STRING },
  { "-q", TRUE, sqdARG_NONE },
  { "--informat",   FALSE, sqdARG_STRING },
  { "--nosmall",    FALSE, sqdARG_NONE },
  { "--regress",    FALSE, sqdARG_STRING },
  { "--banded",     FALSE, sqdARG_NONE },
  { "--bandp",      FALSE, sqdARG_FLOAT},
  { "--bandexpand", FALSE, sqdARG_NONE},
  { "--banddump"  , FALSE, sqdARG_INT},
  { "-W", TRUE, sqdARG_INT },
  { "--full", FALSE, sqdARG_NONE },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv)
{
  char            *cmfile;      /* file to read CM from */	
  CMFILE          *cmfp;	/* open CM file */
  char            *seqfile;     /* file to read sequences from */
  int              format;      /* format of sequence file */
  CM_t            *cm;          /* a covariance model       */
  char           **rseq;        /* RNA sequences to align */
  int              nseq;	/* number of rseq */
  SQINFO          *sqinfo;      /* optional info attached to sequences */
  char           **dsq;         /* digitized RNA sequences */
  Parsetree_t    **tr;          /* parse trees for the sequences */
  MSA             *msa;         /* alignment that's created */
  char            *outfile;	/* optional output file name */
  FILE            *ofp;         /* an open output file */
  int              i;		/* counter over sequences/parsetrees */

  float            sc;		/* score for one sequence alignment */
  float            maxsc;	/* max score in all seqs */
  float            minsc;	/* min score in all seqs */
  float            avgsc;	/* avg score over all seqs */
  
  char  *regressfile;           /* regression test data file */
  int    be_quiet;		/* TRUE to suppress verbose output & banner */
  int    do_local;		/* TRUE to config the model in local mode   */
  int    do_small;		/* TRUE to do divide and conquer alignments */
  int    do_banded;             /* TRUE to do banded CYKDivideAndConquer */
  int    windowlen;             /* window length for calculating bands */
  int    do_full;               /* TRUE to output all match columns in output alignment */
  int    do_bandexpand;         /* TRUE to naively expand bands when necessary */
  int    bdump_level;           /* verbosity level for --banddump option, 0 is OFF */
  int    expand_flag;           /* TRUE if the dmin and dmax vectors have 
                                   just been expanded (in which case we want
				   to recalculate them before we align a new
				   sequence), and FALSE if they're as calculated
				   in BandBounds given gamma from BandDistribution */
  /*EPN 08.18.05*/
  int    set_window;            /* TRUE to set window length due to -W option*/


  double  **gamma;              /* P(subseq length = n) for each state v    */
  int     *dmin;		/* minimum d bound for state v, [0..v..M-1] */
  int     *dmax; 		/* maximum d bound for state v, [0..v..M-1] */
  double   bandp;		/* tail loss probability for banding */
  
  char *optname;                /* name of option found by Getopt()        */
  char *optarg;                 /* argument found by Getopt()              */
  int   optind;                 /* index in argv[]                         */

  /*EPN 11.13.05 */
  int      safe_windowlen;	/* initial windowlen (W) used for calculating bands
				 * in BandCalculationEngine().
				 * this needs to be significantly bigger than what
				 * we expect dmax[0] to be, for truncation error
				 * handling.
				 */
  Stopwatch_t     *watch;

  
  /*********************************************** 
   * Parse command line
   ***********************************************/

  format      = SQFILE_UNKNOWN;
  be_quiet    = FALSE;
  do_local    = FALSE;
  do_small    = TRUE;
  outfile     = NULL;
  regressfile = NULL;
  windowlen   = 200;
  set_window  = FALSE;
  do_banded   = FALSE;
  bandp       = 0.0001;
  do_full     = FALSE;
  do_bandexpand = FALSE;
  bdump_level = 0;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-l")          == 0) do_local    = TRUE;
    else if (strcmp(optname, "-o")          == 0) outfile     = optarg;
    else if (strcmp(optname, "-q")          == 0) be_quiet    = TRUE;
    else if (strcmp(optname, "--nosmall")   == 0) do_small    = FALSE;
    else if (strcmp(optname, "--regress")   == 0) regressfile = optarg;
    else if (strcmp(optname, "--banded")    == 0) do_banded   = TRUE;
    else if (strcmp(optname, "--bandp")    == 0) bandp        = atof(optarg);
    else if (strcmp(optname, "-W")          == 0) 
      { windowlen    = atoi(optarg); set_window = TRUE; }
    else if (strcmp(optname, "--full")      == 0) do_full      = TRUE;
    else if (strcmp(optname, "--bandexpand")== 0) do_bandexpand= TRUE;
    else if (strcmp(optname, "--banddump")== 0)  bdump_level = atoi(optarg);
    else if (strcmp(optname, "--informat")  == 0) {
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

  if (bdump_level > 3) Die("Highest available --banddump verbosity level is 3\n%s", usage);
  if (do_bandexpand && (!(do_banded))) Die("Doesn't make sense to use --bandexpand option with --banded option\n", usage);
  if (argc - optind != 2) Die("Incorrect number of arguments.\n%s\n", usage);
  cmfile = argv[optind++];
  seqfile = argv[optind++]; 


  /* Try to work around inability to autodetect from a pipe or .gz:
   * assume FASTA format
   */
  if (format == SQFILE_UNKNOWN &&
      (Strparse("^.*\\.gz$", seqfile, 0) || strcmp(seqfile, "-") == 0))
    format = SQFILE_FASTA;
  
  /*****************************************************************
   * Input and configure the CM
   *****************************************************************/

  watch = StopwatchCreate();

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL)
    Die("Failed to open covariance model save file %s\n%s\n", cmfile, usage);
  if (! CMFileRead(cmfp, &cm))
    Die("Failed to read a CM from %s -- file corrupt?\n", cmfile);
  if (cm == NULL) 
    Die("%s empty?\n", cmfile);
  CMFileClose(cmfp);

  /* EPN 08.18.05 */
  if (! (set_window)) windowlen = cm->W;
  /*printf("***cm->W : %d***\n", cm->W);*/
  printf("***cm->el_selfsc: %f\n", cm->el_selfsc);
  /* EPN 11.18.05 Now that know what windowlen is, we need to ensure that
   * cm->el_selfsc * W >= IMPOSSIBLE (cm->el_selfsc is the score for an EL self transition)
   * This is done because we are potentially multiply cm->el_selfsc * W, and adding
   * that to IMPOSSIBLE. To avoid underflow issues this value must be less than
   * 3 * IMPOSSIBLE. Here we guarantee its less than 2 * IMPOSSIBLE (to be safe).
   */
  if((cm->el_selfsc * windowlen) < IMPOSSIBLE)
    cm->el_selfsc = (IMPOSSIBLE / (windowlen+1));

  if (do_local) ConfigLocal(cm, 0.5, 0.5);
  CMLogoddsify(cm);
  CMHackInsertScores(cm);	/* "TEMPORARY" fix for bad priors */

  /*****************************************************************
   * Input and digitize the unaligned sequences
   *****************************************************************/

  if (! ReadMultipleRseqs(seqfile, format, &rseq, &sqinfo, &nseq))
    Die("Failed to read any sequences from file %s", seqfile);
  dsq = MallocOrDie(sizeof(char *) * nseq);
  for (i = 0; i < nseq; i++) 
    dsq[i] = DigitizeSequence(rseq[i], sqinfo[i].len);

  /*********************************************** 
   * Show the banner
   ***********************************************/

  if (! be_quiet) 
    {
      MainBanner(stdout, banner);
      printf(   "CM file:              %s\n", cmfile);
      printf(   "Sequence file:        %s\n", seqfile);
      printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
    }

  /*****************************************************************
   * Do the work: 
   *        collect parse trees for each sequence
   *        run them thru Parsetrees2Alignment().
   *****************************************************************/
  
  tr    = MallocOrDie(sizeof(Parsetree_t) * nseq);
  minsc = FLT_MAX;
  maxsc = -FLT_MAX;
  avgsc = 0;

  if(do_banded || bdump_level > 0)
    {
      safe_windowlen = windowlen * 2;
      while(!(BandCalculationEngine(cm, safe_windowlen, bandp, 0, &dmin, &dmax, &gamma, do_local)))
	{
	  FreeBandDensities(cm, gamma);
	  free(dmin);
	  free(dmax);
	  safe_windowlen *= 2;
	}

      /* EPN 11.13.05 
       * An important design decision.
       * We're changing the windowlen value here. By default,
       * windowlen is read from the cm file (set to cm->W). 
       * Here we're doing a banded alignment though. Its pointless to allow
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
      if(bdump_level > 1) 
	{
	  printf("bandp:%f\n", bandp);
	  debug_print_bands(cm, dmin, dmax);
	}
    }
  expand_flag = FALSE;

  StopwatchZero(watch);
  StopwatchStart(watch);

  for (i = 0; i < nseq; i++) 
    {
      if(do_bandexpand)
	{
	  if(expand_flag)
	    {
	      /* The bands were expanded from aligning the previous
	       * sequence, we want to return them to the default,
	       * the current, inefficient approach is to recalculate
	       * gamma and the bands.
	       */
	      FreeBandDensities(cm, gamma);
	      free(dmin);
	      free(dmax);
	      while(!(BandCalculationEngine(cm, safe_windowlen, bandp, 0, &dmin, &dmax, &gamma, do_local)))
		{
		  FreeBandDensities(cm, gamma);
		  free(dmin);
		  free(dmax);
		  safe_windowlen *= 2;
		}
	      expand_flag = FALSE;
	    }
	  if((sqinfo[i].len < dmin[0]) || (sqinfo[i].len > dmax[0]))
	    {
	      /* the query sequence we're aligning is longer than
		 the upper limit band on the root node, or is 
		 shorter than the lower limit on the band on
		 the root node, so we expand */
	      ExpandBands(cm, sqinfo[i].len, dmin, dmax);
	      printf("Expanded bands for seq : %s\n", sqinfo[i].name);
	      if(bdump_level > 2) 
		{
		  printf("printing expanded bands :\n");
		  debug_print_bands(cm, dmin, dmax);
		}
	      expand_flag = TRUE;
	    }
	}
      else 
	{
	  if (do_banded && (sqinfo[i].len < dmin[0] || sqinfo[i].len > dmax[0]))
	    {
	      /* the query sequence we're aligning is longer than
		 the upper limit band on the root node, or is 
		 shorter than the lower limit on the band on
		 the root node, so we quit and tell user
	         they can try naive band expansion if they want.*/
	      Die("Length of sequence to align (%d nt) lies outside the root band.\ndmin[0]: %d and dmax[0]: %d\nImpossible to align with banded CYK unless you try --bandexpand.\n%s", sqinfo[i].len, dmin[0], dmax[0], usage);
	    }
	}
      printf("Aligning %s\n", sqinfo[i].name);
      /*debug_print_bands(cm, dmin, dmax);*/

      if (do_small) 
	{
	  if(do_banded)
	    {
	      sc = CYKDivideAndConquer_b(cm, dsq[i], sqinfo[i].len, 0, 1, sqinfo[i].len, 
					      &(tr[i]), dmin, dmax);
	      if(bdump_level > 0)
		banded_trace_info_dump(cm, tr[i], dmin, dmax, bdump_level);
	    }
	  else
	    {
	      sc = CYKDivideAndConquer(cm, dsq[i], sqinfo[i].len, 0, 1, sqinfo[i].len, &(tr[i]));
	      if(bdump_level > 0)
		{
		  /* We want band info but --banded wasn't used.  Useful if you're curious
		   * why a banded parse is crappy relative to non-banded parse, e.g. allows you 
		   * to see where the non-banded parse went outside the bands.
		   */
		  banded_trace_info_dump(cm, tr[i], dmin, dmax, bdump_level);
		}
	    }
	}
      else if(do_banded)
	{
	  sc = CYKInside_b(cm, dsq[i], sqinfo[i].len, 0, 1, sqinfo[i].len, &(tr[i]), dmin, dmax);
	  if(bdump_level > 0)
	    banded_trace_info_dump(cm, tr[i], dmin, dmax, bdump_level);
	}
      else
	{
	  sc = CYKInside(cm, dsq[i], sqinfo[i].len, 0, 1, sqinfo[i].len, &(tr[i]));
	  if(bdump_level > 0)
	    {
	      /* We want band info but --banded wasn't used.  Useful if you're curious
	       * why a banded parse is crappy relative to non-banded parse, e.g. allows you 
	       * to see where the non-banded parse went outside the bands.
	       */
	      banded_trace_info_dump(cm, tr[i], dmin, dmax, bdump_level);
	    }
	}
      avgsc += sc;
      if (sc > maxsc) maxsc = sc;
      if (sc < minsc) minsc = sc;

    }
  avgsc /= nseq;

  if(do_full)
    {
      msa = Parsetrees2Alignment_full(cm, dsq, sqinfo, NULL, tr, nseq);
    }
  else
    {
      msa = Parsetrees2Alignment(cm, dsq, sqinfo, NULL, tr, nseq);
    }
  /*****************************************************************
   * Output the alignment.
   *****************************************************************/

  if (outfile != NULL && (ofp = fopen(outfile, "w")) != NULL) 
    {
      WriteStockholm(ofp, msa);
      printf("Alignment saved in file %s\n", outfile);
      fclose(ofp);
    }
  else
    WriteStockholm(stdout, msa);
  
  if (regressfile != NULL && (ofp = fopen(regressfile, "w")) != NULL) 
    {
      /* Must delete author info from msa, because it contains version
       * and won't diff clean in regression tests.
       */
      free(msa->au); msa->au = NULL;
      WriteStockholm(ofp, msa);
      fclose(ofp);
    }

  /*****************************************************************
   * Clean up and exit.
   *****************************************************************/

  StopwatchStop(watch);
  StopwatchDisplay(stdout, "\nCPU time: ", watch);

  for (i = 0; i < nseq; i++) 
    {
      FreeParsetree(tr[i]);
      free(dsq[i]);
      FreeSequence(rseq[i], &(sqinfo[i]));
    }

  if (do_banded)
    {
      FreeBandDensities(cm, gamma);
      free(dmin);
      free(dmax);
    }

  MSAFree(msa);
  FreeCM(cm);
  free(rseq);
  free(sqinfo);
  free(dsq);
  free(tr);
  
  StopwatchFree(watch);
  SqdClean();
  return 0;
}




MSA *
Parsetrees2Alignment(CM_t *cm, char **dsq, SQINFO *sqinfo, float *wgt, 
		     Parsetree_t **tr, int nseq)
{
  MSA         *msa;          /* multiple sequence alignment */
  CMEmitMap_t *emap;         /* consensus emit map for the CM */
  int          i;            /* counter over traces */
  int          v, nd;        /* state, node indices */
  int          cpos;         /* counter over consensus positions (0)1..clen */
  int         *matuse;       /* TRUE if we need a cpos in mult alignment */
  int         *iluse;        /* # of IL insertions after a cpos for 1 trace */
  int         *eluse;        /* # of EL insertions after a cpos for 1 trace */
  int         *iruse;        /* # of IR insertions after a cpos for 1 trace */
  int         *maxil;        /* max # of IL insertions after a cpos */
  int         *maxel;        /* max # of EL insertions after a cpos */
  int         *maxir;        /* max # of IR insertions after a cpos */
  int	      *matmap;       /* apos corresponding to a cpos */
  int         *ilmap;        /* first apos for an IL following a cpos */
  int         *elmap;        /* first apos for an EL following a cpos */
  int         *irmap;        /* first apos for an IR following a cpos */
  int          alen;	     /* length of msa in columns */
  int          apos;	     /* position in an aligned sequence in MSA */
  int          rpos;	     /* position in an unaligned sequence in dsq */
  int          tpos;         /* position in a parsetree */
  int          el_len;	     /* length of an EL insertion in residues */
  CMConsensus_t *con;        /* consensus information for the CM */
  int          prvnd;	     /* keeps track of previous node for EL */

  emap = CreateEmitMap(cm);

  matuse = MallocOrDie(sizeof(int)*(emap->clen+1));   
  iluse  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  eluse  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  iruse  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  maxil  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  maxel  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  maxir  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  matmap = MallocOrDie(sizeof(int)*(emap->clen+1));   
  ilmap  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  elmap  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  irmap  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  
  for (cpos = 0; cpos <= emap->clen; cpos++) 
    {
      matuse[cpos] = 0;
      maxil[cpos] = maxel[cpos] = maxir[cpos] = 0;
      ilmap[cpos] = elmap[cpos] = irmap[cpos] = 0;
    }

  /* Look at all the traces; find maximum length of
   * insert needed at each of the clen+1 possible gap
   * points. (There are three types of insert, IL/EL/IR.)
   * Also find whether we don't need some of the match
   * (consensus) columns.
   */
  for (i = 0; i < nseq; i++) 
    {
      for (cpos = 0; cpos <= emap->clen; cpos++) 
	iluse[cpos] = eluse[cpos] = iruse[cpos] = 0;

      for (tpos = 0; tpos < tr[i]->n; tpos++)
	{
	  v  = tr[i]->state[tpos];
	  if (cm->sttype[v] == EL_st) nd = prvnd;
	  else                        nd = cm->ndidx[v];
	  
	  switch (cm->sttype[v]) {
	  case MP_st: 
	    matuse[emap->lpos[nd]] = 1;
	    matuse[emap->rpos[nd]] = 1;
	    break;
	  case ML_st:
	    matuse[emap->lpos[nd]] = 1;
	    break;
	  case MR_st:
	    matuse[emap->rpos[nd]] = 1;
	    break;
	  case IL_st:
	    iluse[emap->lpos[nd]]++;
	    break;
	  case IR_st:		
            /* remember, convention on rpos is that IR precedes this
             * cpos. Make it after the previous cpos, hence the -1. 
	     */
	    iruse[emap->rpos[nd]-1]++;
	    break;
	  case EL_st:
	    el_len = tr[i]->emitr[tpos] - tr[i]->emitl[tpos] + 1;
	    eluse[emap->epos[nd]] = el_len;
            /* not possible to have >1 EL in same place; could assert this */
	    break;
	  }

	  prvnd = nd;
	} /* end looking at trace i */

      for (cpos = 0; cpos <= emap->clen; cpos++) 
	{
	  if (iluse[cpos] > maxil[cpos]) maxil[cpos] = iluse[cpos];
	  if (eluse[cpos] > maxel[cpos]) maxel[cpos] = eluse[cpos];
	  if (iruse[cpos] > maxir[cpos]) maxir[cpos] = iruse[cpos];
	}
    } /* end calculating lengths used by all traces */
  

  /* Now we can calculate the total length of the multiple alignment, alen;
   * and the maps ilmap, elmap, and irmap that turn a cpos into an apos
   * in the multiple alignment: e.g. for an IL that follows consensus position
   * cpos, put it at or after apos = ilmap[cpos] in aseq[][].
   * IR's are filled in backwards (3'->5') and rightflushed.
   */
  alen = 0;
  for (cpos = 0; cpos <= emap->clen; cpos++)
    {
      if (matuse[cpos]) {
	matmap[cpos] = alen; 
	alen++;
      } else 
	matmap[cpos] = -1;

      ilmap[cpos] = alen; alen += maxil[cpos];
      elmap[cpos] = alen; alen += maxel[cpos];
      alen += maxir[cpos]; irmap[cpos] = alen-1; 
    }

  /* We're getting closer.
   * Now we can allocate for the MSA.
   */
  msa = MSAAlloc(nseq, alen);
  msa->nseq = nseq;
  msa->alen = alen;

  for (i = 0; i < nseq; i++)
    {
      /* Initialize the aseq with all pads '.' (in insert cols) 
       * and deletes '-' (in match cols).
       */
      for (apos = 0; apos < alen; apos++)
	msa->aseq[i][apos] = '.';
      for (cpos = 0; cpos <= emap->clen; cpos++)
	if (matmap[cpos] != -1) msa->aseq[i][matmap[cpos]] = '-';
      msa->aseq[i][alen] = '\0';

      /* Traverse this guy's trace, and place all his
       * emitted residues.
       */
      for (cpos = 0; cpos <= emap->clen; cpos++)
	iluse[cpos] = iruse[cpos] = 0;

      for (tpos = 0; tpos < tr[i]->n; tpos++) 
	{
	  v  = tr[i]->state[tpos];
	  if (cm->sttype[v] == EL_st) nd = prvnd;
	  else                        nd = cm->ndidx[v];

	  switch (cm->sttype[v]) {
	  case MP_st:
	    cpos = emap->lpos[nd];
	    apos = matmap[cpos];
	    rpos = tr[i]->emitl[tpos];
	    msa->aseq[i][apos] = Alphabet[(int) dsq[i][rpos]];

	    cpos = emap->rpos[nd];
	    apos = matmap[cpos];
	    rpos = tr[i]->emitr[tpos];
	    msa->aseq[i][apos] = Alphabet[(int) dsq[i][rpos]];
	    break;
	    
	  case ML_st:
	    cpos = emap->lpos[nd];
	    apos = matmap[cpos];
	    rpos = tr[i]->emitl[tpos];
	    msa->aseq[i][apos] = Alphabet[(int) dsq[i][rpos]];
	    break;

	  case MR_st:
	    cpos = emap->rpos[nd];
	    apos = matmap[cpos];
	    rpos = tr[i]->emitr[tpos];
	    msa->aseq[i][apos] = Alphabet[(int) dsq[i][rpos]];
	    break;

	  case IL_st:
	    cpos = emap->lpos[nd];
	    apos = ilmap[cpos] + iluse[cpos];
	    rpos = tr[i]->emitl[tpos];
	    msa->aseq[i][apos] = tolower((int) Alphabet[(int) dsq[i][rpos]]);
	    iluse[cpos]++;
	    break;

	  case EL_st: 
            /* we can assert eluse[cpos] always == 0 when we enter,
	     * because we can only have one EL insertion event per 
             * cpos. If we ever decide to regularize (split) insertions,
             * though, we'll want to calculate eluse in the rpos loop.
             */
	    cpos = emap->epos[nd]; 
	    apos = elmap[cpos]; 
	    for (rpos = tr[i]->emitl[tpos]; rpos <= tr[i]->emitr[tpos]; rpos++)
	      {
		msa->aseq[i][apos] = tolower((int) Alphabet[(int) dsq[i][rpos]]);
		apos++;
	      }
	    break;

	  case IR_st: 
	    cpos = emap->rpos[nd]-1;  /* -1 converts to "following this one" */
	    apos = irmap[cpos] - iruse[cpos];  /* writing backwards, 3'->5' */
	    rpos = tr[i]->emitr[tpos];
	    msa->aseq[i][apos] = tolower((int) Alphabet[(int) dsq[i][rpos]]);
	    iruse[cpos]++;
	    break;

	  case D_st:
	    if (cm->stid[v] == MATP_D || cm->stid[v] == MATL_D) 
	      {
		cpos = emap->lpos[nd];
		if (matuse[cpos]) msa->aseq[i][matmap[cpos]] = '-';
	      }
	    if (cm->stid[v] == MATP_D || cm->stid[v] == MATR_D) 
	      {
		cpos = emap->rpos[nd];
		if (matuse[cpos]) msa->aseq[i][matmap[cpos]] = '-';
	      }
	    break;

	  } /* end of the switch statement */


	  prvnd = nd;
	} /* end traversal over trace i. */

      /* Here is where we could put some insert-regularization code
       * a la HMMER: reach into each insert, find a random split point,
       * and shove part of it flush-right. But, for now, don't bother.
       */

    } /* end loop over all parsetrees */


  /* Gee, wasn't that easy?
   * Add the rest of the ("optional") information to the MSA.
   */
  con = CreateCMConsensus(cm, 3.0, 1.0);

  /* "author" info */
  msa->au   = MallocOrDie(sizeof(char) * (strlen(PACKAGE_VERSION)+10));
  sprintf(msa->au, "Infernal %s", PACKAGE_VERSION);

  for (i = 0; i < nseq; i++)
    {
      msa->sqname[i] = sre_strdup(sqinfo[i].name, -1);
      msa->sqlen[i]  = sqinfo[i].len;
      if (sqinfo[i].flags & SQINFO_ACC)
        MSASetSeqAccession(msa, i, sqinfo[i].acc);
      if (sqinfo[i].flags & SQINFO_DESC)
        MSASetSeqDescription(msa, i, sqinfo[i].desc);

      /* TODO: individual SS annotations
       */
      
      if (wgt == NULL) msa->wgt[i] = 1.0;
      else             msa->wgt[i] = wgt[i];
    }

  /* Construct the secondary structure consensus line, msa->ss_cons:
   *       IL, IR are annotated as .
   *       EL is annotated as ~
   *       and match columns use the structure code.
   * Also the primary sequence consensus/reference coordinate system line,
   * msa->rf.
   */
  msa->ss_cons = MallocOrDie(sizeof(char) * (alen+1));
  msa->rf = MallocOrDie(sizeof(char) * (alen+1));
  for (cpos = 0; cpos <= emap->clen; cpos++) 
    {
      if (matuse[cpos]) 
	{ /* CMConsensus is off-by-one right now, 0..clen-1 relative to cpos's 1..clen */

	  /* bug i1, xref STL7 p.12. Before annotating something as a base pair,
	   * make sure the paired column is also present.
	   */
	  if (con->ct[cpos-1] != -1 && matuse[con->ct[cpos-1]+1] == 0) {
	    msa->ss_cons[matmap[cpos]] = '.';
	    msa->rf[matmap[cpos]]      = con->cseq[cpos-1];
	  } else {
	    msa->ss_cons[matmap[cpos]] = con->cstr[cpos-1];	
	    msa->rf[matmap[cpos]]      = con->cseq[cpos-1];
	  }
	}
      if (maxil[cpos] > 0) 
	for (apos = ilmap[cpos]; apos < ilmap[cpos] + maxil[cpos]; apos++)
	  {
	    msa->ss_cons[apos] = '.';
	    msa->rf[apos] = '.';
	  }
      if (maxel[cpos] > 0)
	for (apos = elmap[cpos]; apos < elmap[cpos] + maxel[cpos]; apos++)
	  {
	    msa->ss_cons[apos] = '~';
	    msa->rf[apos] = '~';
	  }
      if (maxir[cpos] > 0)	/* remember to write backwards */
	for (apos = irmap[cpos]; apos > irmap[cpos] - maxir[cpos]; apos--)
	  {
	    msa->ss_cons[apos] = '.';
	    msa->rf[apos] = '.';
	  }
    }
  msa->ss_cons[alen] = '\0';
  msa->rf[alen] = '\0';

  FreeCMConsensus(con);
  FreeEmitMap(emap);
  free(matuse);
  free(iluse);
  free(eluse);
  free(iruse);
  free(maxil);
  free(maxel);
  free(maxir);
  free(matmap);
  free(ilmap);
  free(elmap);
  free(irmap);
  return msa;
}

/***********************************************************************************
 * EPN 06.08.05
 * Function : Parsetrees2Alignment_full()
 * Purpose  : Create a MSA with ALL consensus match columns included, even if they
 *            are all gaps.  Called if the --full option is used.  Function was
 *            created to make it possible to 'merge' alignment files outputted from
 *            different cmalign runs that used the same CM (useful if we're splitting
 *            up big alignments onto multiple nodes of the cluster).  
 * 
 *            For notes see ~nawrocki/lab/rRNA/infernal_0426/src/output_full_RF_060805.
 * 
 *            This function is derived from Parsetrees2Alignment() (above) and
 *            only has minimal changes to achieve the desired effect.  The capability
 *            to output a 'full' alignment could easily be added to the original
 *            Parsetrees2Alignment() function, but to be safe, I left that function
 *            unscathed (cause it works, and I didn't write it), and just made this
 *            new function so as not to screw up cmalign when the --full option is
 *            not used.  This is a temporary fix to adding --full functionality.
 ***********************************************************************************/

MSA *
Parsetrees2Alignment_full(CM_t *cm, char **dsq, SQINFO *sqinfo, float *wgt, 
			  Parsetree_t **tr, int nseq)
{
  MSA         *msa;          /* multiple sequence alignment */
  CMEmitMap_t *emap;         /* consensus emit map for the CM */
  int          i;            /* counter over traces */
  int          v, nd;        /* state, node indices */
  int          cpos;         /* counter over consensus positions (0)1..clen */
  int         *matuse;       /* TRUE if we need a cpos in mult alignment */
  int         *iluse;        /* # of IL insertions after a cpos for 1 trace */
  int         *eluse;        /* # of EL insertions after a cpos for 1 trace */
  int         *iruse;        /* # of IR insertions after a cpos for 1 trace */
  int         *maxil;        /* max # of IL insertions after a cpos */
  int         *maxel;        /* max # of EL insertions after a cpos */
  int         *maxir;        /* max # of IR insertions after a cpos */
  int	      *matmap;       /* apos corresponding to a cpos */
  int         *ilmap;        /* first apos for an IL following a cpos */
  int         *elmap;        /* first apos for an EL following a cpos */
  int         *irmap;        /* first apos for an IR following a cpos */
  int          alen;	     /* length of msa in columns */
  int          apos;	     /* position in an aligned sequence in MSA */
  int          rpos;	     /* position in an unaligned sequence in dsq */
  int          tpos;         /* position in a parsetree */
  int          el_len;	     /* length of an EL insertion in residues */
  CMConsensus_t *con;        /* consensus information for the CM */
  int          prvnd;	     /* keeps track of previous node for EL */

  emap = CreateEmitMap(cm);

  matuse = MallocOrDie(sizeof(int)*(emap->clen+1));   
  iluse  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  eluse  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  iruse  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  maxil  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  maxel  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  maxir  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  matmap = MallocOrDie(sizeof(int)*(emap->clen+1));   
  ilmap  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  elmap  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  irmap  = MallocOrDie(sizeof(int)*(emap->clen+1));   
  
  for (cpos = 0; cpos <= emap->clen; cpos++) 
    {
      /* EPN 06.08.05 main difference between this function 
	 and Parsetrees2Alignment() is the following line */
      /* original line : matuse[cpos] = 0;*/
      if(cpos == 0)
	/*consensus positions are 1..clen*/
	matuse[cpos] = 0;
      else
	matuse[cpos] = 1;
      maxil[cpos] = maxel[cpos] = maxir[cpos] = 0;
      ilmap[cpos] = elmap[cpos] = irmap[cpos] = 0;
    }

  /* Look at all the traces; find maximum length of
   * insert needed at each of the clen+1 possible gap
   * points. (There are three types of insert, IL/EL/IR.)
   * Also find whether we don't need some of the match
   * (consensus) columns.
   */
  for (i = 0; i < nseq; i++) 
    {
      for (cpos = 0; cpos <= emap->clen; cpos++) 
	iluse[cpos] = eluse[cpos] = iruse[cpos] = 0;

      for (tpos = 0; tpos < tr[i]->n; tpos++)
	{
	  v  = tr[i]->state[tpos];
	  if (cm->sttype[v] == EL_st) nd = prvnd;
	  else                        nd = cm->ndidx[v];
	  
	  switch (cm->sttype[v]) {
	  case MP_st: 
	    matuse[emap->lpos[nd]] = 1;
	    matuse[emap->rpos[nd]] = 1;
	    break;
	  case ML_st:
	    matuse[emap->lpos[nd]] = 1;
	    break;
	  case MR_st:
	    matuse[emap->rpos[nd]] = 1;
	    break;
	  case IL_st:
	    iluse[emap->lpos[nd]]++;
	    break;
	  case IR_st:		
            /* remember, convention on rpos is that IR precedes this
             * cpos. Make it after the previous cpos, hence the -1. 
	     */
	    iruse[emap->rpos[nd]-1]++;
	    break;
	  case EL_st:
	    el_len = tr[i]->emitr[tpos] - tr[i]->emitl[tpos] + 1;
	    eluse[emap->epos[nd]] = el_len;
            /* not possible to have >1 EL in same place; could assert this */
	    break;
	  }

	  prvnd = nd;
	} /* end looking at trace i */

      for (cpos = 0; cpos <= emap->clen; cpos++) 
	{
	  if (iluse[cpos] > maxil[cpos]) maxil[cpos] = iluse[cpos];
	  if (eluse[cpos] > maxel[cpos]) maxel[cpos] = eluse[cpos];
	  if (iruse[cpos] > maxir[cpos]) maxir[cpos] = iruse[cpos];
	}
    } /* end calculating lengths used by all traces */
  

  /* Now we can calculate the total length of the multiple alignment, alen;
   * and the maps ilmap, elmap, and irmap that turn a cpos into an apos
   * in the multiple alignment: e.g. for an IL that follows consensus position
   * cpos, put it at or after apos = ilmap[cpos] in aseq[][].
   * IR's are filled in backwards (3'->5') and rightflushed.
   */
  alen = 0;
  for (cpos = 0; cpos <= emap->clen; cpos++)
    {
      if (matuse[cpos]) {
	matmap[cpos] = alen; 
	alen++;
      } else 
	matmap[cpos] = -1;

      ilmap[cpos] = alen; alen += maxil[cpos];
      elmap[cpos] = alen; alen += maxel[cpos];
      alen += maxir[cpos]; irmap[cpos] = alen-1; 
    }

  /* We're getting closer.
   * Now we can allocate for the MSA.
   */
  msa = MSAAlloc(nseq, alen);
  msa->nseq = nseq;
  msa->alen = alen;

  for (i = 0; i < nseq; i++)
    {
      /* Initialize the aseq with all pads '.' (in insert cols) 
       * and deletes '-' (in match cols).
       */
      for (apos = 0; apos < alen; apos++)
	msa->aseq[i][apos] = '.';
      for (cpos = 0; cpos <= emap->clen; cpos++)
	if (matmap[cpos] != -1) msa->aseq[i][matmap[cpos]] = '-';
      msa->aseq[i][alen] = '\0';

      /* Traverse this guy's trace, and place all his
       * emitted residues.
       */
      for (cpos = 0; cpos <= emap->clen; cpos++)
	iluse[cpos] = iruse[cpos] = 0;

      for (tpos = 0; tpos < tr[i]->n; tpos++) 
	{
	  v  = tr[i]->state[tpos];
	  if (cm->sttype[v] == EL_st) nd = prvnd;
	  else                        nd = cm->ndidx[v];

	  switch (cm->sttype[v]) {
	  case MP_st:
	    cpos = emap->lpos[nd];
	    apos = matmap[cpos];
	    rpos = tr[i]->emitl[tpos];
	    msa->aseq[i][apos] = Alphabet[(int) dsq[i][rpos]];

	    cpos = emap->rpos[nd];
	    apos = matmap[cpos];
	    rpos = tr[i]->emitr[tpos];
	    msa->aseq[i][apos] = Alphabet[(int) dsq[i][rpos]];
	    break;
	    
	  case ML_st:
	    cpos = emap->lpos[nd];
	    apos = matmap[cpos];
	    rpos = tr[i]->emitl[tpos];
	    msa->aseq[i][apos] = Alphabet[(int) dsq[i][rpos]];
	    break;

	  case MR_st:
	    cpos = emap->rpos[nd];
	    apos = matmap[cpos];
	    rpos = tr[i]->emitr[tpos];
	    msa->aseq[i][apos] = Alphabet[(int) dsq[i][rpos]];
	    break;

	  case IL_st:
	    cpos = emap->lpos[nd];
	    apos = ilmap[cpos] + iluse[cpos];
	    rpos = tr[i]->emitl[tpos];
	    msa->aseq[i][apos] = tolower((int) Alphabet[(int) dsq[i][rpos]]);
	    iluse[cpos]++;
	    break;

	  case EL_st: 
            /* we can assert eluse[cpos] always == 0 when we enter,
	     * because we can only have one EL insertion event per 
             * cpos. If we ever decide to regularize (split) insertions,
             * though, we'll want to calculate eluse in the rpos loop.
             */
	    cpos = emap->epos[nd]; 
	    apos = elmap[cpos]; 
	    for (rpos = tr[i]->emitl[tpos]; rpos <= tr[i]->emitr[tpos]; rpos++)
	      {
		msa->aseq[i][apos] = tolower((int) Alphabet[(int) dsq[i][rpos]]);
		apos++;
	      }
	    break;

	  case IR_st: 
	    cpos = emap->rpos[nd]-1;  /* -1 converts to "following this one" */
	    apos = irmap[cpos] - iruse[cpos];  /* writing backwards, 3'->5' */
	    rpos = tr[i]->emitr[tpos];
	    msa->aseq[i][apos] = tolower((int) Alphabet[(int) dsq[i][rpos]]);
	    iruse[cpos]++;
	    break;

	  case D_st:
	    if (cm->stid[v] == MATP_D || cm->stid[v] == MATL_D) 
	      {
		cpos = emap->lpos[nd];
		if (matuse[cpos]) msa->aseq[i][matmap[cpos]] = '-';
	      }
	    if (cm->stid[v] == MATP_D || cm->stid[v] == MATR_D) 
	      {
		cpos = emap->rpos[nd];
		if (matuse[cpos]) msa->aseq[i][matmap[cpos]] = '-';
	      }
	    break;

	  } /* end of the switch statement */


	  prvnd = nd;
	} /* end traversal over trace i. */

      /* Here is where we could put some insert-regularization code
       * a la HMMER: reach into each insert, find a random split point,
       * and shove part of it flush-right. But, for now, don't bother.
       */

    } /* end loop over all parsetrees */


  /* Gee, wasn't that easy?
   * Add the rest of the ("optional") information to the MSA.
   */
  con = CreateCMConsensus(cm, 3.0, 1.0);

  /* "author" info */
  msa->au   = MallocOrDie(sizeof(char) * (strlen(PACKAGE_VERSION)+10));
  sprintf(msa->au, "Infernal %s", PACKAGE_VERSION);

  for (i = 0; i < nseq; i++)
    {
      msa->sqname[i] = sre_strdup(sqinfo[i].name, -1);
      msa->sqlen[i]  = sqinfo[i].len;
      if (sqinfo[i].flags & SQINFO_ACC)
        MSASetSeqAccession(msa, i, sqinfo[i].acc);
      if (sqinfo[i].flags & SQINFO_DESC)
        MSASetSeqDescription(msa, i, sqinfo[i].desc);

      /* TODO: individual SS annotations
       */
      
      if (wgt == NULL) msa->wgt[i] = 1.0;
      else             msa->wgt[i] = wgt[i];
    }

  /* Construct the secondary structure consensus line, msa->ss_cons:
   *       IL, IR are annotated as .
   *       EL is annotated as ~
   *       and match columns use the structure code.
   * Also the primary sequence consensus/reference coordinate system line,
   * msa->rf.
   */
  msa->ss_cons = MallocOrDie(sizeof(char) * (alen+1));
  msa->rf = MallocOrDie(sizeof(char) * (alen+1));
  for (cpos = 0; cpos <= emap->clen; cpos++) 
    {
      if (matuse[cpos]) 
	{ /* CMConsensus is off-by-one right now, 0..clen-1 relative to cpos's 1..clen */

	  /* bug i1, xref STL7 p.12. Before annotating something as a base pair,
	   * make sure the paired column is also present.
	   */
	  if (con->ct[cpos-1] != -1 && matuse[con->ct[cpos-1]+1] == 0) {
	    msa->ss_cons[matmap[cpos]] = '.';
	    msa->rf[matmap[cpos]]      = con->cseq[cpos-1];
	  } else {
	    msa->ss_cons[matmap[cpos]] = con->cstr[cpos-1];	
	    msa->rf[matmap[cpos]]      = con->cseq[cpos-1];
	  }
	}
      if (maxil[cpos] > 0) 
	for (apos = ilmap[cpos]; apos < ilmap[cpos] + maxil[cpos]; apos++)
	  {
	    msa->ss_cons[apos] = '.';
	    msa->rf[apos] = '.';
	  }
      if (maxel[cpos] > 0)
	for (apos = elmap[cpos]; apos < elmap[cpos] + maxel[cpos]; apos++)
	  {
	    msa->ss_cons[apos] = '~';
	    msa->rf[apos] = '~';
	  }
      if (maxir[cpos] > 0)	/* remember to write backwards */
	for (apos = irmap[cpos]; apos > irmap[cpos] - maxir[cpos]; apos--)
	  {
	    msa->ss_cons[apos] = '.';
	    msa->rf[apos] = '.';
	  }
    }
  msa->ss_cons[alen] = '\0';
  msa->rf[alen] = '\0';

  FreeCMConsensus(con);
  FreeEmitMap(emap);
  free(matuse);
  free(iluse);
  free(eluse);
  free(iruse);
  free(maxil);
  free(maxel);
  free(maxir);
  free(matmap);
  free(ilmap);
  free(elmap);
  free(irmap);
  return msa;
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

/* EPN 07.22.05
 * ExpandBands()
 * Function: ExpandBands
 *
 * Purpose:  Called when the sequence we are about to align 
 *           using bands is either shorter in length than
 *           the dmin on the root state, or longer in length
 *           than the dmax on the root state.
 *            
 *           This function expands the bands on ALL states
 *           v=1..cm->M-1 in the following manner :
 *           
 *           case 1 : target len < dmin[0]
 *                    subtract (dmin[0]-target len) from
 *                    dmin of all states, and ensure
 *                    dmin[v]>=0 for all v.
 *                    Further :
 *                    if cm->sttype[v] == MP_st ensure dmin[v]>=2;
 *                    if cm->sttype[v] == IL_st || ML_st ensure dmin[v]>=1;
 *                    if cm->sttype[v] == IR_st || MR_st ensure dmin[v]>=1;
 *                        
 *           case 2 : target len > dmax[0]
 *                    add (target len-dmax[0] to dmax
 *                    of all states.
 *
 *           Prior to handling such situtations with this
 *           hack, the program would choke and die.  This
 *           hacky approach is used as a simple, inefficient
 *           not well thought out, but effective way to 
 *           solve this problem.
 * 
 * Args:    cm       - the CM
 *          tlen     - length of target sequence about to be aligned
 *          dmin     - minimum d bound for each state v; [0..v..M-1]
 *                     may be modified in this function
 *          dmax     - maximum d bound for each state v; [0..v..M-1]
 *                     may be modified in this function
 *
 * Returns: (void) 
 */

static void
ExpandBands(CM_t *cm, int tlen, int *dmin, int *dmax)
{
  int v;
  int diff;
  int root_min;
  int root_max;
  int M = cm->M;
  root_min = dmin[0];
  root_max = dmax[0];

  if(tlen < root_min)
    {
      diff = root_min - tlen;
      for(v=0; v<M; v++)
	{
	  dmin[v] -= diff;
	  if((cm->sttype[v] == MP_st) && (dmin[v] < 2)) 
	    dmin[v] = 2;
	  else if(((cm->sttype[v] == IL_st) || (cm->sttype[v] == ML_st)) 
		  && (dmin[v] < 1)) 
	    dmin[v] = 1;
	  else if(((cm->sttype[v] == IR_st) || (cm->sttype[v] == MR_st)) 
		  && (dmin[v] < 1)) 
	    dmin[v] = 1;
	  else if(dmin[v] < 0) 
	    dmin[v] = 0;
	}
    }
  else if(tlen > root_max)
    {
      diff = tlen - root_min;
      for(v=0; v<M; v++)
	{
	  dmax[v] += diff;
	}
    }
  printf("Expanded bands : \n");
}

/* EPN 08.15.05
 * banded_trace_info_dump()
 * Function: banded_trace_info_dump
 *
 * Purpose:  Called when the user has enabled the --banddump
 *           options.  This function determines how close the
 *           trace was to the bands at each state in the trace,
 *           and prints out that information in differing levels
 *           of verbosity depending on an input parameter 
 *           (bdump_level).
 * 
 * Args:    tr       - the parsetree (trace)
 *          dmin     - minimum d bound for each state v; [0..v..M-1]
 *                     may be modified in this function
 *          dmax     - maximum d bound for each state v; [0..v..M-1]
 *                     may be modified in this function
 *          bdump_level - level of verbosity
 * Returns: (void) 
 */

static void
banded_trace_info_dump(CM_t *cm, Parsetree_t *tr, int *dmin, int *dmax, int bdump_level)
{
  char **sttypes;
  char **nodetypes;
  int v, i, j, d, tpos;
  int mindiff;            /* d - dmin[v] */
  int maxdiff;            /* dmax[v] - d */

  sttypes = malloc(sizeof(char *) * 10);
  sttypes[0] = "D";
  sttypes[1] = "MP";
  sttypes[2] = "ML";
  sttypes[3] = "MR";
  sttypes[4] = "IL";
  sttypes[5] = "IR";
  sttypes[6] = "S";
  sttypes[7] = "E";
  sttypes[8] = "B";
  sttypes[9] = "EL";

  nodetypes = malloc(sizeof(char *) * 8);
  nodetypes[0] = "BIF";
  nodetypes[1] = "MATP";
  nodetypes[2] = "MATL";
  nodetypes[3] = "MATR";
  nodetypes[4] = "BEGL";
  nodetypes[5] = "BEGR";
  nodetypes[6] = "ROOT";
  nodetypes[7] = "END";

  for (tpos = 0; tpos < tr->n; tpos++)
    {
      v  = tr->state[tpos];
      i = tr->emitl[tpos];
      j = tr->emitr[tpos];
      d = j-i+1;

      if(cm->sttype[v] != EL_st)
	{
	  mindiff = d-dmin[v];
	  maxdiff = dmax[v]-d;
	  if(bdump_level > 1 || ((mindiff < 0) || (maxdiff < 0)))
	    printf("%-4s %-3s v: %4d | d: %4d | dmin: %4d | dmax: %4d | %3d | %3d |\n", nodetypes[(int) cm->ndtype[cm->ndidx[v]]], sttypes[(int) cm->sttype[v]], v, d, dmin[v], dmax[v], mindiff, maxdiff);
	}
      else
	{
	  if(bdump_level > 1)
	    printf("%-8s v: %4d | d: %4d |\n", sttypes[(int) cm->sttype[v]], v, d);
	}
    }
  free(sttypes);
  free(nodetypes);
}

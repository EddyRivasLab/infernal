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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>

#include "squid.h"		/* general sequence analysis library    */
#include "msa.h"                /* squid's multiple alignment i/o       */

#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* external functions                   */
#include "version.h"            /* versioning info for Infernal         */

MSA *Parsetrees2Alignment(CM_t *cm, char **dsq, SQINFO *sqinfo, float *wgt, 
			  Parsetree_t **tr, int nseq);

static char banner[] = "cmalign - align sequences to an RNA CM";

static char usage[]  = "\
Usage: cmalign [-options] <cmfile> <sequence file>\n\
  Available options are:\n\
   -h     : help; print brief help on version and usage\n\
   -l     : local; align locally w.r.t. the model\n\
   -o <f> : output the alignment file to file <f>\n\
   -q     : quiet; suppress verbose banner\n\
";

static char experts[] = "\
   --informat <s>: specify that input alignment is in format <s>\n\
   --nosmall     : use normal alignment algorithm, not d&c\n\
   --regress <f> : save regression test data to file <f>\n\
";

static struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE }, 
  { "-l", TRUE, sqdARG_NONE },
  { "-o", TRUE, sqdARG_STRING },
  { "-q", TRUE, sqdARG_NONE },
  { "--informat",   FALSE, sqdARG_STRING },
  { "--nosmall",    FALSE, sqdARG_NONE },
  { "--regress",    FALSE, sqdARG_STRING },
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

  char *optname;                /* name of option found by Getopt()        */
  char *optarg;                 /* argument found by Getopt()              */
  int   optind;                 /* index in argv[]                         */

  /*********************************************** 
   * Parse command line
   ***********************************************/

  format      = SQFILE_UNKNOWN;
  be_quiet    = FALSE;
  do_local    = FALSE;
  do_small    = TRUE;
  outfile     = NULL;
  regressfile = NULL;
  
  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-l")          == 0) do_local    = TRUE;
    else if (strcmp(optname, "-o")          == 0) outfile     = optarg;
    else if (strcmp(optname, "-q")          == 0) be_quiet    = TRUE;
    else if (strcmp(optname, "--nosmall")   == 0) do_small    = FALSE;
    else if (strcmp(optname, "--regress")   == 0) regressfile = optarg;
    else if (strcmp(optname, "--informat")  == 0) {
      format = String2SeqfileFormat(optarg);
      if (format == SQFILE_UNKNOWN) 
	Die("unrecognized sequence file format \"%s\"", optarg);
    }
    else if (strcmp(optname, "-h") == 0) {
      Banner(stdout, banner);
      puts(usage);
      puts(experts);
      exit(EXIT_SUCCESS);
    }
  }

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

  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL)
    Die("Failed to open covariance model save file %s\n%s\n", cmfile, usage);
  if (! CMFileRead(cmfp, &cm))
    Die("Failed to read a CM from %s -- file corrupt?\n", cmfile);
  if (cm == NULL) 
    Die("%s empty?\n", cmfile);
  CMFileClose(cmfp);

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
      Banner(stdout, banner);
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
  for (i = 0; i < nseq; i++) 
    {
      if (do_small) 
	sc = CYKDivideAndConquer(cm, dsq[i], sqinfo[i].len, 0, 1, sqinfo[i].len, &(tr[i]));
      else
	sc = CYKInside(cm, dsq[i], sqinfo[i].len, 0, 1, sqinfo[i].len, &(tr[i]));
      
      avgsc += sc;
      if (sc > maxsc) maxsc = sc;
      if (sc < minsc) minsc = sc;
    }
  avgsc /= nseq;

  msa = Parsetrees2Alignment(cm, dsq, sqinfo, NULL, tr, nseq);

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

  for (i = 0; i < nseq; i++) 
    {
      FreeParsetree(tr[i]);
      free(dsq[i]);
      FreeSequence(rseq[i], &(sqinfo[i]));
    }
  MSAFree(msa);
  FreeCM(cm);
  free(rseq);
  free(sqinfo);
  free(dsq);
  free(tr);
  
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
  msa->au   = MallocOrDie(sizeof(char) * (strlen(RELEASE)+10));
  sprintf(msa->au, "Infernal %s", RELEASE);

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

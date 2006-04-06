/************************************************************
 * @LICENSE@
 ************************************************************/

/* hmmalign.c
 * SRE, Thu Dec 18 16:05:29 1997 [St. Louis]
 * 
 * main() for aligning a set of sequences to an HMM.
 * CVS $Id: hmmalign.c 912 2003-10-03 19:12:59Z eddy $
 */ 

#include "hmmer_config.h"		/* compile-time configuration constants */
#include "squidconf.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hmmer_structs.h"		/* data structures, macros, #define's   */
#include "hmmer_funcs.h"		/* function declarations                */
#include "globals.h"		/* alphabet global variables            */
#include "squid.h"		/* general sequence analysis library    */
#include "msa.h"		/* squid's multiple alignment i/o       */

static char banner[] = "hmmalign - align sequences to an HMM profile";

static char usage[]  = "\
Usage: hmmalign [-options] <hmm file> <sequence file>\n\
Available options are:\n\
   -h     : help; print brief help on version and usage\n\
   -m     : only print symbols aligned to match states\n\
   -o <f> : save alignment in file <f>\n\
   -q     : quiet - suppress verbose banner\n\
   -p     : append posterior probabilities to output alignment\n\
";

static char experts[] = "\
   --informat <s>  : sequence file is in format <s>\n\
   --mapali <f>    : include alignment in file <f> using map in HMM\n\
   --oneline       : output Stockholm fmt with 1 line/seq, not interleaved\n\
   --outformat <s> : output alignment in format <s> [default: Stockholm]\n\
                       formats include: MSF, Clustal, Phylip, SELEX\n\
   --withali <f>   : include alignment to (fixed) alignment in file <f>\n\
\n";

static struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE   }, 
  { "-m", TRUE, sqdARG_NONE   } ,
  { "-o", TRUE, sqdARG_STRING },
  { "-q", TRUE, sqdARG_NONE   },
  { "-p", TRUE, sqdARG_NONE   }, 
  { "--informat",  FALSE, sqdARG_STRING },
  { "--mapali",    FALSE, sqdARG_STRING },
  { "--oneline",   FALSE, sqdARG_NONE },
  { "--outformat", FALSE, sqdARG_STRING },
  { "--withali",   FALSE, sqdARG_STRING },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

static void include_alignment(char *seqfile, struct plan7_s *hmm, int do_mapped,
			      char ***rseq, unsigned char ***dsq, SQINFO **sqinfo, 
			      struct p7trace_s ***tr, int *nseq);

static void actually_write_stockholm_post(FILE *fp, MSA *msa, int cpl, int do_post);

int
main(int argc, char **argv) 
{
  char            *hmmfile;	/* file to read HMMs from                  */
  HMMFILE         *hmmfp;       /* opened hmmfile for reading              */
  struct plan7_s  *hmm;         /* HMM to align to                         */ 
  char            *seqfile;     /* file to read target sequence from       */ 
  int              format;	/* format of seqfile                       */
  char           **rseq;        /* raw, unaligned sequences                */ 
  SQINFO          *sqinfo;      /* info associated with sequences          */
  unsigned char  **dsq;         /* digitized raw sequences                 */
  int              nseq;        /* number of sequences                     */  
  float           *wgt;		/* weights to assign to alignment          */
  MSA             *msa;         /* alignment that's created                */    
  int              i;
  struct dpmatrix_s *mx;        /* growable DP matrix                      */
  struct dpmatrix_s *fwd;       /* growable DP matrix for forwards         */
  struct dpmatrix_s *bck;       /* growable DP matrix for backwards        */
  struct dpmatrix_s *posterior; /* growable DP matrix for backwards        */
  struct p7trace_s **tr;        /* traces for aligned sequences            */

  char *optname;                /* name of option found by Getopt()         */
  char *optarg;                 /* argument found by Getopt()               */
  int   optind;                 /* index in argv[]                          */
  int   be_quiet;		/* TRUE to suppress verbose banner          */
  int   matchonly;		/* TRUE to show only match state syms       */
  char *outfile;                /* optional alignment output file           */
  int   outfmt;			/* code for output alignment format         */
  int   do_oneline;             /* TRUE to do one line/seq, no interleaving */
  FILE *ofp;                    /* handle on alignment output file          */
  char *withali;                /* name of additional alignment file to align */
  char *mapali;                 /* name of additional alignment file to map   */
  int  do_decode;               /* TRUE to do posterior decode              */

  char **postcode;              /* posterior decode array of strings        */
  char *apostcode;              /* aligned posterior decode array           */
  /*********************************************** 
   * Parse command line
   ***********************************************/
  
  format     = SQFILE_UNKNOWN;	  /* default: autodetect input format     */
  outfmt     = MSAFILE_STOCKHOLM; /* default: output in Stockholm format  */
  do_oneline = FALSE;		  /* default: interleaved format          */
  matchonly  = FALSE;
  outfile    = NULL;		  /* default: output alignment to stdout  */
  be_quiet   = FALSE;
  withali    = NULL;
  mapali     = NULL;
  do_decode  = FALSE;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-m")        == 0) matchonly  = TRUE;
    else if (strcmp(optname, "-o")        == 0) outfile    = optarg;
    else if (strcmp(optname, "-q")        == 0) be_quiet   = TRUE; 
    else if (strcmp(optname, "-p")        == 0) do_decode  = TRUE; 
    else if (strcmp(optname, "--mapali")  == 0) mapali     = optarg;
    else if (strcmp(optname, "--oneline") == 0) do_oneline = TRUE;
    else if (strcmp(optname, "--withali") == 0) withali    = optarg;
    else if (strcmp(optname, "--informat") == 0) 
      {
	format = String2SeqfileFormat(optarg);
	if (format == SQFILE_UNKNOWN) 
	  Die("unrecognized sequence file format \"%s\"", optarg);
      }
    else if (strcmp(optname, "--outformat") == 0) 
      {
	outfmt = String2SeqfileFormat(optarg);
	if (outfmt == MSAFILE_UNKNOWN) 
	  Die("unrecognized output alignment file format \"%s\"", optarg);
	if (! IsAlignmentFormat(outfmt))
	  Die("\"%s\" is not a multiple alignment format", optarg);
	if ((do_decode) && (outfmt != MSAFILE_STOCKHOLM))
	  Die("Uh... can't output posterior probs on format other than Stockholm\nNo good reason really ... just laziness");
      }
    else if (strcmp(optname, "-h") == 0) 
      {
	HMMERBanner(stdout, banner);
	puts(usage);
	puts(experts);
	exit(EXIT_SUCCESS);
      }
  }
  if (argc - optind != 2)
    Die("Incorrect number of arguments.\n%s\n", usage);

  hmmfile = argv[optind++];
  seqfile = argv[optind++]; 

  /* Try to work around inability to autodetect from a pipe or .gz:
   * assume FASTA format
   */
  if (format == SQFILE_UNKNOWN &&
      (Strparse("^.*\\.gz$", seqfile, 0) || strcmp(seqfile, "-") == 0))
    format = SQFILE_FASTA;

 /*********************************************** 
  * Open HMM file (might be in HMMERDB or current directory).
  * Read a single HMM from it.
  * 
  * Currently hmmalign disallows the J state and
  * only allows one domain per sequence. To preserve
  * the S/W entry information, the J state is explicitly
  * disallowed, rather than calling a Plan7*Config() function.
  * this is a workaround in 2.1 for the 2.0.x "yo!" bug.
  ***********************************************/


  if ((hmmfp = HMMFileOpen(hmmfile, "HMMERDB")) == NULL)
    Die("Failed to open HMM file %s\n%s", hmmfile, usage);
  printf("hmmer_Alphabet is %s\n", hmmer_Alphabet);
  if (!HMMFileRead(hmmfp, &hmm)) 
    Die("Failed to read any HMMs from %s\n", hmmfile);
  HMMFileClose(hmmfp);
  printf("hmmer_Alphabet is %s\n", hmmer_Alphabet);

  if (hmm == NULL) 
    Die("HMM file %s corrupt or in incorrect format? Parse failed", hmmfile);
  hmm->xt[XTE][MOVE] = 1.;	      /* only 1 domain/sequence ("global" alignment) */
  hmm->xt[XTE][LOOP] = 0.;
  P7Logoddsify(hmm, TRUE);
				/* do we have the map we might need? */


  if (mapali != NULL && ! (hmm->flags & PLAN7_MAP))
    Die("HMMER: HMM file %s has no map; you can't use --mapali.", hmmfile);

  /*********************************************** 
   * Open sequence file in current directory.
   * Read all seqs from it.
   ***********************************************/

  if (! ReadMultipleRseqs(seqfile, format, &rseq, &sqinfo, &nseq))
    Die("Failed to read any sequences from file %s", seqfile);

  /*********************************************** 
   * Show the banner
   ***********************************************/

  if (! be_quiet) 
    {
      HMMERBanner(stdout, banner);
      printf(   "HMM file:             %s\n", hmmfile);
      printf(   "Sequence file:        %s\n", seqfile);
      printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
    }

  /*********************************************** 
   * Do the work
   ***********************************************/

  /* Allocations and initializations.
   */
  dsq = MallocOrDie(sizeof(unsigned char *)    * nseq);
  tr  = MallocOrDie(sizeof(struct p7trace_s *) * nseq);
  mx  = CreatePlan7Matrix(1, hmm->M, 25, 0);
  if(do_decode)
    {
      /*P7Forward() and P7Backward() allocate a new matrix
	in each call.  I *should* change them to be more
	efficient, but I don't want to mess with core_algorithms.c */
      /*      fwd  = CreatePlan7Matrix(1, hmm->M, 25, 0);
	      bck  = CreatePlan7Matrix(1, hmm->M, 25, 0); */

      postcode = malloc(sizeof(char *) * nseq);
    }

  /* Align each sequence to the model, collect traces
   */
  for (i = 0; i < nseq; i++)
    {
      printf("rseq[%d] is %s\n", i, rseq[i]);
      dsq[i] = hmmer_DigitizeSequence(rseq[i], sqinfo[i].len);

      if (P7ViterbiSpaceOK(sqinfo[i].len, hmm->M, mx))
	(void) P7Viterbi(dsq[i], sqinfo[i].len, hmm, mx, &(tr[i]));
      else
	(void) P7SmallViterbi(dsq[i], sqinfo[i].len, hmm, mx, &(tr[i]));

      if(do_decode)
	{
	  printf("Alphabet_type: %d\n", Alphabet_type);
	  P7Forward(dsq[i], sqinfo[i].len, hmm, &fwd);
	  P7Backward(dsq[i], sqinfo[i].len, hmm, &bck);
	  posterior = bck;
	  P7EmitterPosterior(sqinfo[i].len, hmm, fwd, bck, posterior);
	  postcode[i] = PostalCode(sqinfo[i].len, posterior, tr[i]);
	  FreePlan7Matrix(fwd);
	  FreePlan7Matrix(bck);
	}
      
    }
  FreePlan7Matrix(mx);

  /* Include an aligned alignment, if desired.
   */
  if (mapali != NULL)
    include_alignment(mapali, hmm, TRUE, &rseq, &dsq, &sqinfo, &tr, &nseq);
  if (withali != NULL) 
    include_alignment(withali, hmm, FALSE, &rseq, &dsq, &sqinfo, &tr, &nseq);

  /* Turn traces into a multiple alignment
   */ 
  wgt = MallocOrDie(sizeof(float) * nseq);
  FSet(wgt, nseq, 1.0);
  
  msa = P7Traces2Alignment(dsq, sqinfo, wgt, nseq, hmm->M, tr, matchonly, NULL);

  /* Append the POST markup lines to the msa */
  if(do_decode)
    {
      for (i = 0; i < nseq; i++)
	{
	  MakeAlignedString(msa->aseq[i], msa->alen, postcode[i], &apostcode); 
	  MSAAppendGR(msa, "POST", i, apostcode);
	  free(apostcode);
	  free(postcode[i]);
	}
      free(postcode);
    }

  /*********************************************** 
   * Output the alignment
   ***********************************************/
  
  if (do_decode)
    {
      if (outfile != NULL && (ofp = fopen(outfile, "w")) != NULL)
	{
	  /* use our revoltingly hacked actually_write_stockholm_post() 
	     function */
	  if(!do_oneline)
	    actually_write_stockholm_post(ofp, msa, 50, do_decode);
	  else
	    actually_write_stockholm_post(ofp, msa, msa->alen, do_decode);
	  printf("Alignment saved in file %s\n", outfile);
	  fclose(ofp);
	}
      else
	{
	  /* same thing but to stdout */
	  if(!do_oneline)
	    actually_write_stockholm_post(stdout, msa, 50, do_decode);
	  else
	    actually_write_stockholm_post(stdout, msa, msa->alen, do_decode);
	}
    }
  else if (outfile != NULL && (ofp = fopen(outfile, "w")) != NULL)
    {
      MSAFileWrite(ofp, msa, outfmt, do_oneline);
      printf("Alignment saved in file %s\n", outfile);
      fclose(ofp);
    }
  else
    MSAFileWrite(stdout, msa, outfmt, do_oneline);


  /*********************************************** 
   * Cleanup and exit
   ***********************************************/
  
  for (i = 0; i < nseq; i++) 
    {
      P7FreeTrace(tr[i]);
      FreeSequence(rseq[i], &(sqinfo[i]));
      free(dsq[i]);
    }
  MSAFree(msa);
  FreePlan7(hmm);
  free(sqinfo);
  free(rseq);
  free(dsq);
  free(wgt);
  free(tr);

  SqdClean();
  return 0;
}


/* Function: include_alignment()
 * Date:     SRE, Sun Jul  5 15:25:13 1998 [St. Louis]
 *
 * Purpose:  Given the name of a multiple alignment file,
 *           align that alignment to the HMM, and add traces
 *           to an existing array of traces. If do_mapped
 *           is TRUE, we use the HMM's map file. If not,
 *           we use P7ViterbiAlignAlignment().
 *
 * Args:     seqfile  - name of alignment file
 *           hmm      - model to align to
 *           do_mapped- TRUE if we're to use the HMM's alignment map
 *           rsq      - RETURN: array of rseqs to add to
 *           dsq      - RETURN: array of dsq to add to
 *           sqinfo   - RETURN: array of SQINFO to add to
 *           tr       - RETURN: array of traces to add to
 *           nseq     - RETURN: number of seqs           
 *
 * Returns:  new, realloc'ed arrays for rsq, dsq, sqinfo, tr; nseq is
 *           increased to nseq+ainfo.nseq.
 */
static void
include_alignment(char *seqfile, struct plan7_s *hmm, int do_mapped,
		  char ***rsq, unsigned char ***dsq, SQINFO **sqinfo, 
		  struct p7trace_s ***tr, int *nseq)
{
  int format;			/* format of alignment file */
  MSA   *msa;			/* alignment to align to    */
  MSAFILE *afp;
  SQINFO  *newinfo;		/* sqinfo array from msa */
  unsigned char **newdsq;
  char **newrseq;
  int   idx;			/* counter over aseqs       */
  struct p7trace_s *master;     /* master trace             */
  struct p7trace_s **addtr;     /* individual traces for aseq */

  format = MSAFILE_UNKNOWN;	/* invoke Babelfish */
  if ((afp = MSAFileOpen(seqfile, format, NULL)) == NULL)
    Die("Alignment file %s could not be opened for reading", seqfile);
  if ((msa = MSAFileRead(afp)) == NULL)
    Die("Failed to read an alignment from %s\n", seqfile);
  MSAFileClose(afp);
  for (idx = 0; idx < msa->nseq; idx++)
    s2upper(msa->aseq[idx]);
  newinfo = MSAToSqinfo(msa);

				/* Verify checksums before mapping */
  if (do_mapped && GCGMultchecksum(msa->aseq, msa->nseq) != hmm->checksum)
    Die("The checksums for alignment file %s and the HMM alignment map don't match.", 
	seqfile);
				/* Get a master trace */
  if (do_mapped) master = MasterTraceFromMap(hmm->map, hmm->M, msa->alen);
  else           master = P7ViterbiAlignAlignment(msa, hmm);

				/* convert to individual traces */
  ImposeMasterTrace(msa->aseq, msa->nseq, master, &addtr);
				/* add those traces to existing ones */
  *tr = MergeTraceArrays(*tr, *nseq, addtr, msa->nseq);
  
				/* additional bookkeeping: add to dsq, sqinfo */
  *rsq = ReallocOrDie((*rsq), sizeof(char *) * (*nseq + msa->nseq));
  DealignAseqs(msa->aseq, msa->nseq, &newrseq);
  for (idx = *nseq; idx < *nseq + msa->nseq; idx++)
    (*rsq)[idx] = newrseq[idx - (*nseq)];
  free(newrseq);

  *dsq = ReallocOrDie((*dsq), sizeof(unsigned char *) * (*nseq + msa->nseq));
  hmmer_DigitizeAlignment(msa, &newdsq);
  for (idx = *nseq; idx < *nseq + msa->nseq; idx++)
    (*dsq)[idx] = newdsq[idx - (*nseq)];
  free(newdsq);
			/* unnecessarily complex, but I can't be bothered... */
  *sqinfo = ReallocOrDie((*sqinfo), sizeof(SQINFO) * (*nseq + msa->nseq));
  for (idx = *nseq; idx < *nseq + msa->nseq; idx++)
    SeqinfoCopy(&((*sqinfo)[idx]), &(newinfo[idx - (*nseq)]));
  free(newinfo);
  
  *nseq = *nseq + msa->nseq;

				/* Cleanup */
  P7FreeTrace(master);
  MSAFree(msa);
				/* Return */
  return;
}


/* Function: actually_write_stockholm_post()
 * Date : 06.02.05
 *   
 * Purpose : A vicious hack.  The actually_write_stockholm() function
 *           was torn out of squid's stockholm.c.  This version
 *           can deal with the POST markup, whereas the original
 *           could not.  Changing the original would've been messy
 *           add I'd have to make sure that squid behaved and such,
 *           so I took the easier but uglier route and added this
 *           function here which is called only if the -p option
 *           is used.

 *           NOTES FROM actually_write_stockholm() :
 * Date:     SRE, Fri May 21 17:39:22 1999 [St. Louis]
 *
 * Purpose:  Write an alignment in Stockholm format to 
 *           an open file. This is the function that actually
 *           does the work. The API's WriteStockholm()
 *           and WriteStockholmOneBlock() are wrappers.
 *
 * Args:     fp    - file that's open for writing
 *           msa   - alignment to write        
 *           cpl   - characters to write per line in alignment block
 *           do_post - 1 if there are post lines, 0 otherwise
 *                   0 should cause reversion back to original function
 *                   behavior.
 * Returns:  (void)
 */
static void
actually_write_stockholm_post(FILE *fp, MSA *msa, int cpl, int do_post)
{
  int  i, j;
  int  len = 0;
  int  namewidth;
  int  typewidth = 0;		/* markup tags are up to 5 chars long */
  int  markupwidth = 0;		/* #=GR, #=GC are four char wide + 1 space */
  char *buf;
  int  currpos;
  char *s, *tok;
  
  /* Figure out how much space we need for name + markup
   * to keep the alignment in register. Required by Stockholm
   * spec, even though our Stockholm parser doesn't care (Erik's does).
   */
  namewidth = 0;
  for (i = 0; i < msa->nseq; i++)
    if ((len = strlen(msa->sqname[i])) > namewidth) 
      namewidth = len;

  /* Figure out how much space we need for markup tags
   *   markupwidth = always 4 if we're doing markup:  strlen("#=GR")
   *   typewidth   = longest markup tag
   */

  if (msa->ss      != NULL) { markupwidth = 4; typewidth = 2; }
  if (msa->sa      != NULL) { markupwidth = 4; typewidth = 2; }
  for (i = 0; i < msa->ngr; i++)
    if ((len = strlen(msa->gr_tag[i])) > typewidth) typewidth = len;

  if (msa->rf      != NULL) { markupwidth = 4; if (typewidth < 2) typewidth = 2; }
  if (msa->ss_cons != NULL) { markupwidth = 4; if (typewidth < 7) typewidth = 7; }
  if (msa->sa_cons != NULL) { markupwidth = 4; if (typewidth < 7) typewidth = 7; }
  
  /*EPN 06.02.05 - here's the one change */
  if (do_post) { markupwidth = 4; if (typewidth < 4) typewidth = 4; }

  for (i = 0; i < msa->ngc; i++)
    if ((len = strlen(msa->gc_tag[i])) > typewidth) typewidth = len;
  
  buf = MallocOrDie(sizeof(char) * (cpl+namewidth+typewidth+markupwidth+61)); 

  /* Magic Stockholm header
   */
  fprintf(fp, "# STOCKHOLM 1.0\n");

  /* Free text comments
   */
  for (i = 0;  i < msa->ncomment; i++)
    fprintf(fp, "# %s\n", msa->comment[i]);
  if (msa->ncomment > 0) fprintf(fp, "\n");

  /* GF section: per-file annotation
   */
  if (msa->name  != NULL)       fprintf(fp, "#=GF ID    %s\n", msa->name);
  if (msa->acc   != NULL)       fprintf(fp, "#=GF AC    %s\n", msa->acc);
  if (msa->desc  != NULL)       fprintf(fp, "#=GF DE    %s\n", msa->desc);
  if (msa->au    != NULL)       fprintf(fp, "#=GF AU    %s\n", msa->au);
  
  /* Thresholds are hacky. Pfam has two. Rfam has one.
   */
  if      (msa->cutoff_is_set[MSA_CUTOFF_GA1] && msa->cutoff_is_set[MSA_CUTOFF_GA2])
    fprintf(fp, "#=GF GA    %.1f %.1f\n", msa->cutoff[MSA_CUTOFF_GA1], msa->cutoff[MSA_CUTOFF_GA2]);
  else if (msa->cutoff_is_set[MSA_CUTOFF_GA1])
    fprintf(fp, "#=GF GA    %.1f\n", msa->cutoff[MSA_CUTOFF_GA1]);
  if      (msa->cutoff_is_set[MSA_CUTOFF_NC1] && msa->cutoff_is_set[MSA_CUTOFF_NC2])
    fprintf(fp, "#=GF NC    %.1f %.1f\n", msa->cutoff[MSA_CUTOFF_NC1], msa->cutoff[MSA_CUTOFF_NC2]);
  else if (msa->cutoff_is_set[MSA_CUTOFF_NC1])
    fprintf(fp, "#=GF NC    %.1f\n", msa->cutoff[MSA_CUTOFF_NC1]);
  if      (msa->cutoff_is_set[MSA_CUTOFF_TC1] && msa->cutoff_is_set[MSA_CUTOFF_TC2])
    fprintf(fp, "#=GF TC    %.1f %.1f\n", msa->cutoff[MSA_CUTOFF_TC1], msa->cutoff[MSA_CUTOFF_TC2]);
  else if (msa->cutoff_is_set[MSA_CUTOFF_TC1])
    fprintf(fp, "#=GF TC    %.1f\n", msa->cutoff[MSA_CUTOFF_TC1]);

  for (i = 0; i < msa->ngf; i++)
    fprintf(fp, "#=GF %-5s %s\n", msa->gf_tag[i], msa->gf[i]); 
  fprintf(fp, "\n");


  /* GS section: per-sequence annotation
   */
  if (msa->flags & MSA_SET_WGT) 
    {
      for (i = 0; i < msa->nseq; i++) 
	fprintf(fp, "#=GS %-*.*s WT    %.2f\n", namewidth, namewidth, msa->sqname[i], msa->wgt[i]);
      fprintf(fp, "\n");
    }
  if (msa->sqacc != NULL) 
    {
      for (i = 0; i < msa->nseq; i++) 
	if (msa->sqacc[i] != NULL)
	  fprintf(fp, "#=GS %-*.*s AC    %s\n", namewidth, namewidth, msa->sqname[i], msa->sqacc[i]);
      fprintf(fp, "\n");
    }
  if (msa->sqdesc != NULL) 
    {
      for (i = 0; i < msa->nseq; i++) 
	if (msa->sqdesc[i] != NULL)
	  fprintf(fp, "#=GS %*.*s DE    %s\n", namewidth, namewidth, msa->sqname[i], msa->sqdesc[i]);
      fprintf(fp, "\n");
    }
  for (i = 0; i < msa->ngs; i++)
    {
      /* Multiannotated GS tags are possible; for example, 
       *     #=GS foo DR PDB; 1xxx;
       *     #=GS foo DR PDB; 2yyy;
       * These are stored, for example, as:
       *     msa->gs[0][0] = "PDB; 1xxx;\nPDB; 2yyy;"
       * and must be decomposed.
       */
      for (j = 0; j < msa->nseq; j++)
	if (msa->gs[i][j] != NULL)
	  {
	    s = msa->gs[i][j];
	    while ((tok = sre_strtok(&s, "\n", NULL)) != NULL)
	      fprintf(fp, "#=GS %*.*s %5s %s\n", namewidth, namewidth,
		      msa->sqname[j], msa->gs_tag[i], tok);
	  }
      fprintf(fp, "\n");
    }

  /* Alignment section:
   * contains aligned sequence, #=GR annotation, and #=GC annotation
   */
  for (currpos = 0; currpos < msa->alen; currpos += cpl)
    {
      if (currpos > 0) fprintf(fp, "\n");
      for (i = 0; i < msa->nseq; i++)
	{
	  strncpy(buf, msa->aseq[i] + currpos, cpl);
	  buf[cpl] = '\0';	      
	  fprintf(fp, "%-*.*s  %s\n", namewidth+typewidth+markupwidth+1, namewidth+typewidth+markupwidth+1, 
		  msa->sqname[i], buf);

	  if (msa->ss != NULL && msa->ss[i] != NULL) {
	    strncpy(buf, msa->ss[i] + currpos, cpl);
	    buf[cpl] = '\0';	 
	    fprintf(fp, "#=GR %-*.*s SS     %s\n", namewidth, namewidth, msa->sqname[i], buf);
	  }
	  if (msa->sa != NULL && msa->sa[i] != NULL) {
	    strncpy(buf, msa->sa[i] + currpos, cpl);
	    buf[cpl] = '\0';
	    fprintf(fp, "#=GR %-*.*s SA     %s\n", namewidth, namewidth, msa->sqname[i], buf);
	  }
	  for (j = 0; j < msa->ngr; j++)
	    if (msa->gr[j][i] != NULL) {
	      strncpy(buf, msa->gr[j][i] + currpos, cpl);
	      buf[cpl] = '\0';
	      fprintf(fp, "#=GR %-*.*s %-*.*s %s\n", 
		      namewidth, namewidth, msa->sqname[i], markupwidth, markupwidth, msa->gr_tag[j], buf);
	    }
	}
      if (msa->ss_cons != NULL) {
	strncpy(buf, msa->ss_cons + currpos, cpl);
	buf[cpl] = '\0';
	fprintf(fp, "#=GC %-*.*s %s\n", namewidth+typewidth+1, namewidth+typewidth+1, "SS_cons", buf);
      }

      if (msa->sa_cons != NULL) {
	strncpy(buf, msa->sa_cons + currpos, cpl);
	buf[cpl] = '\0';
	fprintf(fp, "#=GC %-*.*s %s\n", namewidth+typewidth, namewidth+typewidth, "SA_cons", buf);
      }

      if (msa->rf != NULL) {
	strncpy(buf, msa->rf + currpos, cpl);
	buf[cpl] = '\0';
	fprintf(fp, "#=GC %-*.*s %s\n", namewidth+typewidth+1, namewidth+typewidth+1, "RF", buf);
      }
      for (j = 0; j < msa->ngc; j++) {
	strncpy(buf, msa->gc[j] + currpos, cpl);
	buf[cpl] = '\0';
	fprintf(fp, "#=GC %-*.*s %s\n", namewidth+typewidth+1, namewidth+typewidth+1, 
		msa->gc_tag[j], buf);
      }
    }
  fprintf(fp, "//\n");
  free(buf);
}


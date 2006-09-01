/************************************************************
 * @LICENSE@
 ************************************************************/

/* cmemit.c
 * EPN, 09.01.06 Janelia Farm
 * based on HMMER-2.3.2's hmmemit.c from SRE
 *  
 * main() for generating sequences from an CM
 */

#include "config.h"		/* compile-time configuration constants */
#include "squidconf.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* function declarations                */
#include "globals.h"		/* alphabet global variables            */
#include "squid.h"		/* general sequence analysis library    */
#include "msa.h"		/* squid's multiple sequence i/o        */

static char banner[] = "cmemit - generate sequences from a covariance model";

static char usage[]  = "\
Usage: cmemit [-options] <hmm file>\n\
Available options are:\n\
   -a     : write generated sequences as an alignment, not FASTA\n\
   -c     : generate a single \"consensus\" sequence\n\
   -h     : help; print brief help on version and usage\n\
   -n <n> : emit <n> sequences (default 10)\n\
   -o <f> : save sequences in file <f>\n\
   -q     : quiet - suppress verbose banner\n\
";

static char experts[] = "\
   --seed <n>     : set random number seed to <n>\n\
   --full         : include all match columns in output alignment\n\
";

static struct opt_s OPTIONS[] = {
  { "-a",        TRUE,  sqdARG_NONE },  
  { "-c",        TRUE,  sqdARG_NONE },  
  { "-h",        TRUE,  sqdARG_NONE }, 
  { "-n",        TRUE,  sqdARG_INT},  
  { "-o",        TRUE,  sqdARG_STRING},
  { "-q",        TRUE,  sqdARG_NONE},  
  { "--seed",    FALSE, sqdARG_INT},
  { "--full",    FALSE, sqdARG_NONE},
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv) 
{
  char            *cmfile;      /* file to read CM from */	
  CMFILE          *cmfp;	/* open CM file */
  char            *seqfile;     /* file to read sequences from */
  CM_t            *cm;          /* CM to generate from */

  FILE            *fp;          /* output file handle                      */
  char            *ofile;       /* output sequence file                    */
  int              L;		/* length of a sequence                    */
  int              i;		/* counter over sequences                  */

  int              nseq;	/* number of seqs to sample                */
  int              seed;	/* random number generator seed            */
  int              be_quiet;	/* TRUE to silence header/footer           */
  int              do_alignment;/* TRUE to output in aligned format        */ 
  int              do_consensus;/* TRUE to do a single consensus seq       */
  int              do_full;     /* TRUE to output all match columns in aln */

  char *optname;                /* name of option found by Getopt()         */
  char *optarg;                 /* argument found by Getopt()               */
  int   optind;                 /* index in argv[]                          */


  /*********************************************** 
   * Parse command line
   ***********************************************/

  nseq         = 10;
  seed         = time ((time_t *) NULL);
  be_quiet     = FALSE;
  do_alignment = FALSE;  
  do_consensus = FALSE;
  do_full      = FALSE;
  ofile        = NULL;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-a")     == 0) do_alignment = TRUE;
    else if (strcmp(optname, "-c")     == 0) do_consensus = TRUE;
    else if (strcmp(optname, "-n")     == 0) nseq         = atoi(optarg); 
    else if (strcmp(optname, "-o")     == 0) ofile        = optarg;
    else if (strcmp(optname, "-q")     == 0) be_quiet     = TRUE;
    else if (strcmp(optname, "--seed") == 0) seed         = atoi(optarg);
    else if (strcmp(optname, "--full") == 0) do_full      = TRUE;
    else if (strcmp(optname, "-h") == 0) 
      {
	MainBanner(stdout, banner);
	puts(usage);
	puts(experts);
	exit(EXIT_SUCCESS);
      }
  }
  if (argc - optind != 1)
    Die("Incorrect number of arguments.\n%s\n", usage);

  cmfile = argv[optind++];

  sre_srandom(seed);

  if (do_alignment && do_consensus)
    Die("Sorry, -a and -c are incompatible.\nUsage:\n%s", usage); 
  if (nseq != 10 && do_consensus)
    Warn("-c (consensus) overrides -n (# of sampled seqs)");
  if (do_full && !do_alignment)
    Die("--full only makes sense with -a\n");

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
  CMLogoddsify(cm);

  /*********************************************** 
   * Open output file, if needed.
   ***********************************************/

   if (ofile == NULL) fp = stdout;
   else {
     if ((fp = fopen(ofile, "w")) == NULL)
       Die("Failed to open output file %s for writing", ofile);
   }

  /*********************************************** 
   * Show the options banner
   ***********************************************/

  if (! be_quiet) 
    {
      MainBanner(stdout, banner);
      printf("CM file:             %s\n", cmfile);
      if (! do_consensus) {
	printf("Number of seqs:       %d\n", nseq);
	printf("Random seed:          %d\n", seed);
      } else {
	printf("Generating consensus sequence.\n");
      }
      printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
    }

    /*********************************************** 
     * Do the work.
     * If we're generating an alignment, we have to collect
     * all our traces, then output. If we're generating unaligned
     * sequences, we can emit one at a time.
     ***********************************************/

  if (do_consensus) 
    {
      CMConsensus_t *con;            /* consensus info for the CM */
      char    *seq;
      SQINFO   sqinfo;      /* info about sequence (name/desc)        */
	
      /* Determine consensus sequence */
      con = CreateCMConsensus(cm, 3.0, 1.0);
      
      seq = MallocOrDie(sizeof(char) * (con->clen+1));
      strcpy(seq, con->cseq);
      L = con->clen;
      strcpy(sqinfo.name, cm->name);
      strcpy(sqinfo.desc, "CM generated consensus sequence [cmemit]");
      
      sqinfo.len = L;
      sqinfo.flags = SQINFO_NAME | SQINFO_DESC | SQINFO_LEN;
      
      WriteSeq(fp, SQFILE_FASTA, seq, &sqinfo);
      free(seq);
      FreeCMConsensus(con);
    }
  else if(do_alignment)
    {
      Parsetree_t **tr;             /* Parsetrees of emitted aligned sequences */
      char    **dsq;                /* digitized sequences                     */
      SQINFO            *sqinfo;    /* info about sequences (name/desc)        */
      MSA               *msa;       /* alignment */
      float             *wgt;
      
      dsq    = MallocOrDie(sizeof(char *)    * nseq);
      tr     = MallocOrDie(sizeof(Parsetree_t) * nseq);
      sqinfo = MallocOrDie(sizeof(SQINFO)             * nseq);
      wgt    = MallocOrDie(sizeof(float)              * nseq);
      FSet(wgt, nseq, 1.0);
      
      for (i = 0; i < nseq; i++)
	{
	    EmitParsetree(cm, &(tr[i]), NULL, &(dsq[i]), &L);
	    sprintf(sqinfo[i].name, "seq%d", i+1);
	    sqinfo[i].len   = L;
	    sqinfo[i].flags = SQINFO_NAME | SQINFO_LEN;
	  }
	
	msa = Parsetrees2Alignment(cm, dsq, sqinfo, NULL, tr, nseq, do_full);
	msa->name = sre_strdup(cm->name, -1);
	msa->desc = sre_strdup("Synthetic sequence alignment generated by cmemit", -1);
	
	/* Output the alignment */
	WriteStockholm(fp, msa);
	
	/* Free memory */
	for (i = 0; i < nseq; i++) 
	  {
	    FreeParsetree(tr[i]);
	    free(dsq[i]);
	  }
	free(sqinfo);
	free(dsq);
	free(wgt);
	MSAFree(msa);
	free(tr);
      }
    else				/* unaligned sequence output */
      {
	Parsetree_t      *tr;         /* generated trace                        */
	char             *dsq;        /* digitized sequence                     */
	char             *seq;        /* alphabetic sequence                    */
	SQINFO            sqinfo;     /* info about sequence (name/len)         */

	for (i = 0; i < nseq; i++)
	  {
	    EmitParsetree(cm, &tr, &seq, &dsq, &L);
	    sprintf(sqinfo.name, "%s-%d", cm->name, i+1);
	    sqinfo.len   = L;
	    sqinfo.flags = SQINFO_NAME | SQINFO_LEN;

	    WriteSeq(fp, SQFILE_FASTA, seq, &sqinfo);
	  
	    FreeParsetree(tr);
	    free(dsq);
	    free(seq);
	  }
      }
    FreeCM(cm);
    
    /* We're done; clean up and exit.
     */
    if (ofile != NULL) {
      fclose(fp);
      if (!be_quiet) printf("Output saved in file %s\n", ofile);
    }
    SqdClean();
    return 0;
}


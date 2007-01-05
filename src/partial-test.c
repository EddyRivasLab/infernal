/* partial-test.c
 * EPN, Thu Nov 30 14:16:04 2006
 * 
 * Test the alignment of partial sequences.
 * Emit sequences from a CM, truncate them and align 
 * them back to it.
 *
 *****************************************************************
 * @LICENSE@
 *****************************************************************  
 */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "squid.h"
#include "sqfuncs.h"
#include "dirichlet.h"
#include "sre_stack.h"

#include "structs.h"
#include "funcs.h"

static void partial_test_AlignSeqsWrapper(cm, pdsq, psqinfo, nseq, &ptr, do_local, do_small, do_qdb,
					  beta, do_hbanded, use_sums, hbandp,
					  do_sub, do_fullsub, fsub_pmass, FALSE, FALSE, FALSE, 
					  FALSE, FALSE, NULL, FALSE, 0, 0, TRUE,
					  spos, epos, &post_spos, &post_epos, &dist_spos,
					  &dist_epos);

static char banner[] = "partial-test - test alignment of partial sequences";

static char usage[] = "\
Usage: partial-test [-options] <cmfile>\n\
  where options are:\n\
  -h     : help; print brief help on version and usage\n\
  -n <n> : number of seqs to emit, truncate and align [default: 100]\n\
  -s <n> : set random number seed to <n> \n\
";

static char experts[] = "\
  --read <f> : read seqs to truncate from file <f>\n\ 
  --sub        : aln w/sub CM for columns b/t HMM predicted start/end points\n\
  --fsub       : aln w/sub CM for structure b/t HMM predicted start/end points\n\
  --cp9        : aln w/CM plan 9 HMM\n\
  --global     : run alignment in global mode [default: local]\n\
  --post <f>   : minimum posterior prob from CP9 F/B to include in sub CM\n\
  --fixlen <n> : fix length of partial seqs to <n>\n\
  --minlen <n> : set minimum length of partial seqs as <n>\n\
  --debug  <n> : set verbosity of debugging print statements to <n> [1..3]\n\
  --histo  <n> : build histogram of HMM posterior probability of start/end\n\
\n\
  * HMM banded alignment related options:\n\
   --hbanded     : use experimental CM plan 9 HMM banded CYK aln algorithm\n\
   --hbandp <f>  : tail loss prob for --hbanded [default: 0.0001]\n\
   --sums        : use posterior sums during HMM band calculation (widens bands)\n\
";

static struct opt_s OPTIONS[] = { 
  { "-h", TRUE, sqdARG_NONE }, 
  { "-n", TRUE, sqdARG_INT },
  { "-s", TRUE, sqdARG_INT },
  { "--read",    FALSE, sqdARG_STRING},
  { "--sub",       FALSE, sqdARG_NONE},
  { "--fsub",      FALSE, sqdARG_FLOAT},
  { "--cp9",       FALSE, sqdARG_NONE },
  { "--global",    FALSE, sqdARG_NONE },
  { "--post",      FALSE, sqdARG_FLOAT },
  { "--fixlen",    FALSE, sqdARG_INT },
  { "--minlen",    FALSE, sqdARG_INT },
  { "--hbanded",   FALSE, sqdARG_NONE },
  { "--hbandp",    FALSE, sqdARG_FLOAT},
  { "--debug",     FALSE, sqdARG_INT},
  { "--sums",      FALSE, sqdARG_NONE},
  { "--histo",     FALSE, sqdARG_NONE}

};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv)
{
  char    *cmfile;		/* file to read CM from */	
  CMFILE  *cmfp;		/* open CM file for reading */
  CM_t    *cm;			/* a covariance model       */
  CM_t    *sub_cm;              /* sub covariance model     */
  CMSubMap_t *submap;
  CMSubInfo_t *subinfo;
  int      do_local;
  int      nseq;                /* number of seqs to aln */
  int      v;			/* counter over states */
  int      sstruct;             /* start position for sub CM */
  int      estruct;             /* end position for sub CM */
  int      ncols;               /* number of consensus (match) columns in CM */
  int      i;                   /* counter over sub CMs */
  int      j;                   /* counter */
  int      temp;

  double   pthresh;		
  int      seed;		/* random number seed for MC */

  char *optname;                /* name of option found by Getopt()        */
  char *optarg;                 /* argument found by Getopt()              */
  int   optind;                 /* index in argv[]                         */
  int   do_sub;                 /* TRUE to build a sub CM that models only between 
				 * predicted start/end points. */
  int   do_fullsub;             /* TRUE to build sub CM(s) that model same number of columns
				 * as the template CM, with structure outside sstruct..estruct
				 * removed.                          */
  float fsub_pmass;             /* probability mass from HMM posteriors req'd of start before sstruct
				 * and end after estruct */
  int debug_level;              /* verbosity of debugging print statements */
  int L;
  int do_fixlen;                /* TRUE to set fixed length of partial sequences */
  int fixlen;                   /* if not 0, fixed length of partial sequences */
  int do_minlen;                /* TRUE to set minimum length of partial sequences */
  int minlen;                   /* if not 0, minimum length of partial sequences */
  int passed;

  int Lp;                       /* Length of the partial sequence */
  int pred_spos;                /* predicted consensus column res 1  of partial seq aligns to */
  int pred_epos;                /* predicted consensus column res Lp of partial seq aligns to */
  int s_ct;                     /* num times pred_spos = spos */
  int e_ct;                     /* num times pred_epos = epos */

  /* posterior histogram related variables */
  int    do_histo;              /* TRUE to build histograms */
  float *post_spos;             /* [0..nseq-1] posterior probability from HMM of spos being start */
  float *post_epos;             /* [0..nseq-1] posterior probability from HMM of epos being end */
  int   *dist_spos;             /* [0..nseq-1] distance (+/-) of max post start from actual spos */
  int   *dist_epos;             /* [0..nseq-1] distance (+/-) of max post end   from actual epos */

  /* CYK modes */
  int do_small;                 /* TRUE to use D&C; FALSE to use full CYKInside() */
  int do_qdb;                   /* TRUE to use query dependent bands (QDB)        */
  double qdb_beta;              /* tail loss prob for QDB                         */
  int do_hbanded;               /* TRUE to use CP9 HMM bands for alignment        */
  double hbandp;                /* tail loss prob for HMM bands                   */
  int    use_sums;              /* TRUE to fill and use the posterior sums, false not to. */
  double beta;                  /* tail loss prob for QDB                         */

  /* --read related variables */
  int              do_read;     /* TRUE to read seqs from a file, instead of emitting */
  char            *seqfile;     /* file to read sequences from */
  int              format;      /* format of sequence file */

  
  Parsetree_t **tr;             /* Parsetrees of emitted sequence         */
  Parsetree_t **ptr;            /* Parsetrees of partial emitted sequence */
  char        **seq;            /* actual sequence                        */
  char        **dsq;            /* digitized sequences                    */
  SQINFO       *sqinfo;         /* info about sequences (name/desc)       */
  char        **temp_dsq;       /* test digitized sequence                */
  SQINFO       *temp_sqinfo;    /* info about test sequence (name/desc)   */
  Parsetree_t **temp_tr;        /* Parsetrees of test emitted sequence    */
  
  SQINFO       *psqinfo;        /* info about test sequence (name/desc)   */
  
  int attempts = 0;
  int eq = 0;
  MSA               *msa;       /* alignment */
  int *useme;
  int apos;
  int cc;
  int do_full = FALSE;
  int *spos;
  int *epos;
  char **partial;       /* dealigned seq after truncation            */
  char **pseq;		/* dealigned seqs after truncation           */
  char **pdsq;		/* partial digitized sequences               */
  int x;
  
  CMEmitMap_t *emap; 
  
  /*********************************************** 
   * Parse command line
   ***********************************************/

  seed           = (int) time ((time_t *) NULL);
  pthresh        = 0.1;
  do_sub         = FALSE;
  do_fullsub     = FALSE;
  nseq           =  100;
  debug_level    = 0;
  do_local       = FALSE;
  do_fixlen      = FALSE;
  fixlen         = 0;
  do_minlen      = FALSE;
  minlen         = 0;
  s_ct           = 0;
  e_ct           = 0;
  do_small       = TRUE;
  do_qdb         = FALSE;
  do_hbanded     = FALSE;
  hbandp         = 0.0001;
  beta           = 0.0000001;
  use_sums       = FALSE;
  do_histo       = FALSE;
  fsub_pmass     = 0.;
  do_read        = FALSE;
  seqfile        = NULL;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
		&optind, &optname, &optarg))  {
    if      (strcmp(optname, "-n") == 0) nseq         = atoi(optarg);
    else if (strcmp(optname, "-s") == 0) seed           = atoi(optarg);
    else if (strcmp(optname, "--read")      == 0) { do_read = TRUE; seqfile = optarg; }
    else if (strcmp(optname, "--sub")       == 0) do_sub = TRUE;
    else if (strcmp(optname, "--fsub")      == 0) 
      { do_sub = TRUE; do_fullsub = TRUE; fsub_pmass = atof(optarg); }
    else if (strcmp(optname, "--global")    == 0) do_local = FALSE;
    else if (strcmp(optname, "--post")      == 0) pthresh  = atof(optarg); 
    else if (strcmp(optname, "--fixlen")    == 0) { do_fixlen = TRUE; fixlen   = atoi(optarg); }
    else if (strcmp(optname, "--minlen")    == 0) { do_minlen = TRUE; minlen   = atoi(optarg); }
    else if (strcmp(optname, "--hbanded")   == 0) { do_hbanded = TRUE; do_small = FALSE; }
    else if (strcmp(optname, "--hbandp")    == 0) hbandp       = atof(optarg);
    else if (strcmp(optname, "--sums")      == 0) use_sums     = TRUE;
    else if (strcmp(optname, "--debug")     == 0) debug_level = atoi(optarg);
    else if (strcmp(optname, "--histo")     == 0) do_histo = TRUE;
    else if (strcmp(optname, "-h") == 0) {
      MainBanner(stdout, banner);
      puts(usage);
      puts(experts);
      exit(EXIT_SUCCESS);
    }
  }
  
  if (argc - optind != 1) Die("Incorrect number of arguments.\n%s\n", usage);
  cmfile = argv[optind++];
  
  if(do_histo)
    {
      post_spos = MallocOrDie(sizeof(float) * nseq);
      post_epos = MallocOrDie(sizeof(float) * nseq);
      dist_spos = MallocOrDie(sizeof(int  ) * nseq);
      dist_epos = MallocOrDie(sizeof(int  ) * nseq);
    }

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

  if(do_local) ConfigLocal(cm, 0.5, 0.5);
  CMLogoddsify(cm);

  /* Determine number of consensus columns modelled by CM */
  ncols = 0;
  for(v = 0; v <= cm->M; v++)
    {
      if(cm->stid[v] ==  MATP_MP) ncols += 2;
      else if(cm->stid[v] == MATL_ML || cm->stid[v] == MATR_MR)	ncols++;
    }

  if(do_minlen && do_fixlen)
    Die("Pick either --fixlen or --minlen, not both.\n");
  if(do_minlen && (minlen <= 0 || minlen > ncols))
    Die("--minlen enabled, but non-sensical, must be >= 0 and <= %d (# match cols)\n", ncols);
  if(do_fixlen && (fixlen <= 0 || fixlen > ncols))
    Die("--fixlen enabled, but non-sensical, must be >= 0 and <= %d (# match cols)\n", ncols);

  /*********************************************** 
   * Open output file, if needed.
   ***********************************************/

  /*   if (ofile == NULL) fp = stdout;
   else {
     if ((fp = fopen(ofile, "w")) == NULL)
       Die("Failed to open output file %s for writing", ofile);
       }*/

  /*********************************************** 
   * Show the options banner
   ***********************************************/

  MainBanner(stdout, banner);
  printf("CM file:             %s\n", cmfile);
  printf("Number of seqs:       %d\n", nseq);
  printf("Random seed:          %d\n", seed);
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
  
  /****************************************************** 
   * Do the work. Read or emit seqs and align them.
   ******************************************************/
  if(do_read)
    {
      /* read the sequences from the input file */
      format = SQFILE_FASTA;
      printf("opening file: %s\n", seqfile);
      if (! ReadMultipleRseqs(seqfile, format, &seq, &sqinfo, &nseq))
	Die("Failed to read any sequences from file %s, expecting FASTA.", seqfile);
      dsq = MallocOrDie(sizeof(char *) * nseq);
      for (i = 0; i < nseq; i++) 
	dsq[i] = DigitizeSequence(seq[i], sqinfo[i].len);
    }      

  /* Allocate and initialize */
  if(!do_read)
    {
      dsq    = MallocOrDie(sizeof(char *)      * nseq);
      seq    = MallocOrDie(sizeof(char *)      * nseq);
      sqinfo = MallocOrDie(sizeof(SQINFO)      * nseq);
    }
  tr     = MallocOrDie(sizeof(Parsetree_t) * nseq);
  i = 0;
  
  pdsq    = MallocOrDie(sizeof(char *)      * nseq);
  psqinfo = MallocOrDie(sizeof(SQINFO)      * nseq);
  ptr     = MallocOrDie(sizeof(Parsetree_t) * nseq);
  
  temp_sqinfo = MallocOrDie(sizeof(SQINFO)      * 1);
  temp_dsq    = MallocOrDie(sizeof(char *)      * 1);
  temp_tr     = MallocOrDie(sizeof(Parsetree_t) * 1);

  spos = MallocOrDie(sizeof(int) * nseq);
  epos = MallocOrDie(sizeof(int) * nseq);
  
  emap = CreateEmitMap(cm);

  /* Either use read seqs (if do_read) or emit nseq seqs from the CM,
   * align them to the CM using the desired CYK algorithm (potentially
   * banded), and truncate them to partial sequences. Then realign
   * with the same CYK algorithm and determine how close the partial
   * seq alignment is with the corresponding subalignment of the full
   * sequence.
   */

  for(i = 0; i < nseq; i++)
    {
      if(debug_level >= 0) printf("i: %d\n", i);
      if(!(do_read))
	{
	  attempts++;
	  /* Emit a sequence from the CM, we don't care about the parsetree,
	   * we'll determine this via alignment */
	  EmitParsetree(cm, NULL, &(seq[i]), &(dsq[i]), &L);
	  while(L == 0)
	    {
	      free(seq[i]);
	      free(dsq[i]);
	      EmitParsetree(cm, NULL, &(seq[i]), &(dsq[i]), &L);
	    }
	  strcpy(sqinfo[i].name, cm->name);
	  sqinfo[i].len   = L;
	  sqinfo[i].flags = SQINFO_NAME | SQINFO_LEN;
	  strcpy(temp_sqinfo[0].name, cm->name);
	  temp_sqinfo[0].len   = L;
	} /* end of if(!(do_read)) */
      else
	{
	  strcpy(temp_sqinfo[0].name, sqinfo[i].name);
	  temp_sqinfo[0].len   = sqinfo[i].len;
	}
      temp_tr[0]  = tr[i];
      temp_dsq[0] = dsq[i];
      temp_sqinfo[0].flags = SQINFO_NAME | SQINFO_LEN;
      
      /* Align the sequence, setting do_sub to FALSE */
      AlignSeqsWrapper(cm, temp_dsq, temp_sqinfo, 1, &temp_tr, do_local, 
		       do_small, do_qdb, beta, do_hbanded, use_sums, 
		       hbandp, 
		       FALSE, FALSE, fsub_pmass, /* these are do_sub and do_fullsub */
		       FALSE, FALSE, FALSE, 
		       FALSE, FALSE, NULL, FALSE, 0, 0, TRUE,
		       NULL, NULL, NULL, NULL, NULL, NULL);
      
      
      strcpy(temp_sqinfo[0].name, cm->name);
      msa = Parsetrees2Alignment(cm, temp_dsq, temp_sqinfo, NULL, temp_tr, 1, TRUE);
      /* important, the final variable do_full must be set to TRUE, we want
       * all match columns in the alignment.
       */
      /*msa->name = sre_strdup(cm->name, -1);
	msa->desc = sre_strdup("Synthetic sequence alignment", -1);*/
      
      /* Step 1: pick random start and end consensus position
       *         for truncation
       */
      passed = FALSE;
      while(!passed)
	{
	  spos[i] = ((int) (sre_random() * ncols)) + 1;
	  epos[i] = ((int) (sre_random() * ncols)) + 1;
	  if(spos[i] > epos[i])
	    {
	      temp = spos[i];
	      spos[i] = epos[i];
	      epos[i] = temp;
	      passed = TRUE;
	    }
	  if(!do_minlen && !do_fixlen)
	    passed = TRUE;
	  if(do_minlen && ((epos[i]-spos[i]+1) < minlen))
	    passed = FALSE;
	  if(do_fixlen)
	    if((spos[i] + fixlen - 1) > ncols)
	      passed = FALSE;
	    else
	      epos[i] = spos[i] + fixlen - 1;
	}	      
      /* Truncate the alignment prior to consensus column spos[i] and after 
	 consensus column epos[i] */
      useme = (int *) MallocOrDie (sizeof(int) * (msa->alen+1));
      for (apos = 0, cc = 0; apos < msa->alen; apos++)
	{
	  /* Careful here, placement of cc++ increment is impt, 
	   * we want all inserts between cc=spos[i]-1 and cc=spos[i],
	   * and between cc=epos[i] and cc=epos[i]+1.
	   */
	  if(cc < (spos[i]-1) || cc > epos[i])
	    useme[apos] = 0;
	  else
	    useme[apos] = 1;
	  if (!isgap(msa->rf[apos])) 
	    { 
	      cc++; 
	      if(cc == (epos[i]+1))
		useme[apos] = 0; 
	    }
	}
      /*printf("\n\nDEBUG PRINTING ORIG ALIGNMENT:\n");
	WriteStockholm(fp, msa);
	printf("\n\nDONE DEBUG PRINTING ORIG ALIGNMENT:\n");
	for(apos=0; apos < msa->alen; apos++)
	printf("useme[%d]: %d\n", apos, useme[apos]);
      */
      MSAShorterAlignment(msa, useme);
      free(useme);
      
      /* Get the dealigned sequence, and save it as pdsq[i] */
      DealignAseqs(msa->aseq, msa->nseq, &pseq);
      if(debug_level > 0)
	{
	  printf("i: %d seq: %s\n", i, seq[i]);
	  printf("i: %d pseq: %s\n", i, pseq[0]);
	}
      Lp = strlen(pseq[0]);
      if(Lp == 0) /* don't want 0 len seqs */
	{
	  /* we handle this oddly, decrement i and wait
	   * for next iteration of the loop to increment i
	   * and refill this seq.
	   */
	  free(pseq);
	  i--;
	}
      else
	{
	  pdsq[i] = DigitizeSequence(pseq[0], Lp);
	  psqinfo[i].len   = Lp;
	  strcpy(psqinfo[i].name, cm->name);
	  psqinfo[i].flags = SQINFO_NAME | SQINFO_LEN;
	}
    }

  /* Align all the partial sequences to the CM */
  if(do_histo)
    partial_test_AlignSeqsWrapper(cm, pdsq, psqinfo, nseq, &ptr, do_local, do_small, do_qdb,
				  beta, do_hbanded, use_sums, hbandp,
				  do_sub, do_fullsub, fsub_pmass, FALSE, FALSE, FALSE, 
				  FALSE, FALSE, NULL, FALSE, 0, 0, TRUE,
				  spos, epos, &post_spos, &post_epos, &dist_spos,
				  &dist_epos);
  else
    partial_test_AlignSeqsWrapper(cm, pdsq, psqinfo, nseq, &ptr, do_local, do_small, do_qdb,
				  beta, do_hbanded, use_sums, hbandp,
				  do_sub, do_fullsub, fsub_pmass, FALSE, FALSE, FALSE, 
				  FALSE, FALSE, NULL, FALSE, 0, 0, TRUE,
				  NULL, NULL, NULL, NULL, NULL, NULL);
  s_ct = 0;
  e_ct = 0;
  /* For each sequence, compare the partial alignment with
   * the full alignment */
  for(i = 0; i < nseq; i++)
    {
      Lp = psqinfo[i].len;
      /*printf("dumping parsetree for partial seq of length %d\n", Lp);
	ParsetreeDump(stdout, ptr[i], cm, pdsq[i]);*/

      /* Determine the consensus column the first and last residue 
       * of the partial sequence were aligned to */
      pred_spos = pred_epos = -1;
      for(x = 0; x < ptr[i]->n; x++)
	{
	  /* find consensus column that residue 1 aligns to */
	  if(ptr[i]->emitl[x] == 1 && (cm->sttype[ptr[i]->state[x]] == ML_st ||
				    cm->sttype[ptr[i]->state[x]] == MP_st))
	    pred_spos = emap->lpos[cm->ndidx[ptr[i]->state[x]]];
	  if(ptr[i]->emitl[x] == 1 && cm->sttype[ptr[i]->state[x]] == IL_st)
	    pred_spos = emap->lpos[cm->ndidx[ptr[i]->state[x]]] + 1;
	  
	  if(ptr[i]->emitr[x] == 1 && (cm->sttype[ptr[i]->state[x]] == MR_st ||
				    cm->sttype[ptr[i]->state[x]] == MP_st))
	    pred_spos = emap->rpos[cm->ndidx[ptr[i]->state[x]]];
	  if(ptr[i]->emitr[x] == 1 && cm->sttype[ptr[i]->state[x]] == IR_st)
	    pred_spos = emap->rpos[cm->ndidx[ptr[i]->state[x]]];

	  /* find consensus column that residue Lp aligns to */
	  if(ptr[i]->emitl[x] == Lp && (cm->sttype[ptr[i]->state[x]] == ML_st ||
					cm->sttype[ptr[i]->state[x]] == MP_st))
	    pred_epos = emap->lpos[cm->ndidx[ptr[i]->state[x]]];
	  if(ptr[i]->emitl[x] == Lp && cm->sttype[ptr[i]->state[x]] == IL_st)
	    pred_epos = emap->lpos[cm->ndidx[ptr[i]->state[x]]];
	  
	  if(ptr[i]->emitr[x] == Lp && (cm->sttype[ptr[i]->state[x]] == MR_st ||
					cm->sttype[ptr[i]->state[x]] == MP_st))
	    pred_epos = emap->rpos[cm->ndidx[ptr[i]->state[x]]];
	  if(ptr[i]->emitr[x] == Lp && cm->sttype[ptr[i]->state[x]] == IR_st)
	    {
	      if(cm->ndidx[ptr[i]->state[x]] == 0)
		pred_epos = 0; /* ROOT is special */
	      else
		pred_epos = emap->rpos[cm->ndidx[ptr[i]->state[x]]] - 1;
	    }
	}  
      if(pred_spos == -1)
	Die("pred_spos is still -1!\n");
      if(pred_epos == -1)
	Die("pred_epos is still -1!\n");
      printf("(%4d) S: %4d %4d %4d E: %4d %4d %4d \n", i, spos[i], pred_spos, spos[i]-pred_spos, epos[i], pred_epos, epos[i]-pred_epos);
      
      if(pred_spos == spos[i])
	s_ct++;
      if(pred_epos == epos[i])
	e_ct++;
      /*WriteStockholm(stdout, msa);*/
      /*printf("attempts: %3d passed (i=%3d)\n", attempts, i);*/
    }

  /*if (tr[0] != NULL) FreeParsetree(tr[0]);  
    if (ptr[i] != NULL) FreeParsetree(ptr[i]); 
    free(dsq[0]);*/
  printf("N %d S %d E %d\n", nseq, s_ct, e_ct);
  /*free(sqinfo);*/
  return EXIT_SUCCESS;
}
  
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
partial_test_AlignSeqsWrapper(CM_t *cm, char **dsq, SQINFO *sqinfo, int nseq, Parsetree_t ***ret_tr, int do_local, 
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
  /*float           ***alpha;*/     /* alpha DP matrix for Inside() */
  /*float           ***beta; */     /* beta DP matrix for Inside() */
  /*float           ***post; */     /* post DP matrix for Inside() */
  int             ***alpha;    /* alpha DP matrix for Inside() */
  int             ***beta;     /* beta DP matrix for Inside() */
  int             ***post;     /* post DP matrix for Inside() */


  /* partial-test variables */
  int do_ptest;                 /* TRUE to fill partial-test variables */
  float *post_spos;
  float *post_epos;
  int   *dist_spos;
  int   *dist_epos;

  if(do_fullsub)
    do_sub = TRUE;

  printf("in AlignSeqsWrapper() do_local: %d do_sub: %d do_fullsub: %d do_enforce: %d\n", do_local, do_sub, do_fullsub, do_enforce);

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
  if(do_enforce)
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
	    {
	      sc = IInside_b_jd_me(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
				   BE_PARANOID,	/* non memory-saving mode */
				   NULL, NULL,	/* manage your own matrix, I don't want it */
				   NULL, NULL,	/* manage your own deckpool, I don't want it */
				   do_local,        /* TRUE to allow local begins */
				   cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax); /* j and d bands */
	    }
	  else
	    {
	      sc = IInside(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
			   BE_EFFICIENT,	/* memory-saving mode */
			   NULL, NULL,	/* manage your own matrix, I don't want it */
			   NULL, NULL,	/* manage your own deckpool, I don't want it */
			   do_local);       /* TRUE to allow local begins */
	    }

	}
      else if(do_outside)
	{	
	  if(do_hbanded)
	    {
	      
	      sc = IInside_b_jd_me(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
				   BE_PARANOID,	/* save full alpha so we can run outside */
				   NULL, &alpha,	/* fill alpha, and return it, needed for FOutside() */
				   NULL, NULL,	/* manage your own deckpool, I don't want it */
				   do_local,        /* TRUE to allow local begins */
				   cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax); /* j and d bands */
	      /*do_check = TRUE;*/
	      sc = IOutside_b_jd_me(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
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
	      sc = IInside(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
			   BE_PARANOID,	/* save full alpha so we can run outside */
			   NULL, &alpha,	/* fill alpha, and return it, needed for FOutside() */
			   NULL, NULL,	/* manage your own deckpool, I don't want it */
			   do_local);       /* TRUE to allow local begins */
	      sc = IOutside(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
			    BE_PARANOID,	/* save full beta */
			    NULL, NULL,	/* manage your own matrix, I don't want it */
			    NULL, NULL,	/* manage your own deckpool, I don't want it */
			    do_local,       /* TRUE to allow local begins */
			    alpha,         /* alpha matrix from IInside() */
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
	  post  = MallocOrDie(sizeof(int **) * (cm->M+1));
	  /*
	  for (v = 0; v < cm->M; v++) alpha[v] = NULL;
	  for (v = 0; v < cm->M+1; v++) beta[v] = NULL;
	  */
	  if(do_hbanded)
	    {
	      for (v = 0; v < cm->M; v++)
		{
		  post[v] = NULL;
		  post[v] = Ialloc_jdbanded_vjd_deck(sqinfo[i].len, 1, sqinfo[i].len, cp9b->jmin[v], 
						      cp9b->jmax[v], cp9b->hdmin[v], cp9b->hdmax[v]);
		}
	      post[cm->M] = NULL;
	      post[cm->M] = alloc_vjd_deck(sqinfo[i].len, 1, sqinfo[i].len);
	      sc = IInside_b_jd_me(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
				   BE_PARANOID,	/* save full alpha so we can run outside */
				   NULL, &alpha,	/* fill alpha, and return it, needed for IOutside() */
				   NULL, NULL,	/* manage your own deckpool, I don't want it */
				   do_local,       /* TRUE to allow local begins */
				   cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax); /* j and d bands */
	      sc = IOutside_b_jd_me(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
				    BE_PARANOID,	/* save full beta */
				    NULL, &beta,	/* fill beta, and return it, needed for ICMPosterior() */
				    NULL, NULL,	/* manage your own deckpool, I don't want it */
				    do_local,       /* TRUE to allow local begins */
				    alpha, &alpha,  /* alpha matrix from IInside(), and save it for CMPosterior*/
				    do_check,      /* TRUE to check Outside probs agree with Inside */
				    cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax); /* j and d bands */
	      ICMPosterior_b_jd_me(sqinfo[i].len, cm, alpha, NULL, beta, NULL, post, &post,
				   cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax);
	      postcode[i] = ICMPostalCode_b_jd_me(cm, sqinfo[i].len, post, tr[i],
						  cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax);
	      /*postcode[i] = CMPostalCode_b_jd_me(cm, sqinfo[i].len, post, tr[i],
		cp9b->jmin, cp9b->jmax, cp9b->hdmin, cp9b->hdmax);*/
	    }
	  else
	    {
	      for (v = 0; v < cm->M+1; v++)
		{
		  post[v] = NULL;
		  post[v] = alloc_vjd_deck(sqinfo[i].len, 1, sqinfo[i].len);
		  post[v] = NULL;
		  post[v] = Ialloc_vjd_deck(sqinfo[i].len, 1, sqinfo[i].len);
		}
	      sc = IInside(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
			   BE_PARANOID,	/* save full alpha so we can run outside */
			   NULL, &alpha,	/* fill alpha, and return it, needed for IOutside() */
			   NULL, NULL,	/* manage your own deckpool, I don't want it */
			   do_local);       /* TRUE to allow local begins */
	      sc = IOutside(cm, dsq[i], sqinfo[i].len, 1, sqinfo[i].len,
			    BE_PARANOID,	/* save full beta */
			    NULL, &beta,	/* fill beta, and return it, needed for CMPosterior() */
			    NULL, NULL,	/* manage your own deckpool, I don't want it */
			    do_local,       /* TRUE to allow local begins */
			    alpha, &alpha,  /* alpha matrix from IInside(), and save it for CMPosterior*/
			    do_check);      /* TRUE to check Outside probs agree with Inside */
	      ICMPosterior(sqinfo[i].len, cm, alpha, NULL, beta, NULL, post, &post);
	      if(do_check || TRUE)
		{
		  ICMCheckPosterior(sqinfo[i].len, cm, post);
		  printf("\nPosteriors checked (I).\n\n");
		}
	      postcode[i] = ICMPostalCode(cm, sqinfo[i].len, post, tr[i]);
	      /*postcode[i] = CMPostalCode(cm, sqinfo[i].len, post, tr[i]);*/
	    }

	  /* free post */
	  if(post != NULL)
	    {
	      for (v = 0; v <= (cm->M); v++)
		if (post[v] != NULL) { Ifree_vjd_deck(post[v], 1, sqinfo[i].len); post[v] = NULL;}
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

  printf("leaving AlignSeqsWrapper()\n");
}

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

static char banner[] = "partial-test - test alignment of partial sequences";

static char usage[] = "\
Usage: partial-test [-options] <cmfile>\n\
  where options are:\n\
  -h     : help; print brief help on version and usage\n\
  -n <n> : number of seqs to emit, truncate and align [default: 100]\n\
  -s <n> : set random number seed to <n> \n\
";

static char experts[] = "\
  --sub        : aln w/sub CM for columns b/t HMM predicted start/end points\n\
  --fsub       : aln w/sub CM for structure b/t HMM predicted start/end points\n\
  --cp9        : aln w/CM plan 9 HMM\n\
  --global     : run alignment in global mode [default: local]\n\
  --post <f>   : minimum posterior prob from CP9 F/B to include in sub CM\n\
  --fixlen <n> : fix length of partial seqs to <n>\n\
  --minlen <n> : set minimum length of partial seqs as <n>\n\
  --debug  <n> : set verbosity of debugging print statements to <n> [1..3]\n\
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
  { "--sub",       FALSE, sqdARG_NONE},
  { "--fsub",      FALSE, sqdARG_NONE},
  { "--cp9",       FALSE, sqdARG_NONE },
  { "--global",    FALSE, sqdARG_NONE },
  { "--post",      FALSE, sqdARG_FLOAT },
  { "--fixlen",    FALSE, sqdARG_INT },
  { "--minlen",    FALSE, sqdARG_INT },
  { "--hbanded",   FALSE, sqdARG_NONE },
  { "--hbandp",    FALSE, sqdARG_FLOAT},
  { "--debug",     FALSE, sqdARG_INT},
  { "--sums",      FALSE, sqdARG_NONE}

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

  /* CYK modes */
  int do_small;                 /* TRUE to use D&C; FALSE to use full CYKInside() */
  int do_qdb;                   /* TRUE to use query dependent bands (QDB)        */
  double qdb_beta;              /* tail loss prob for QDB                         */
  int do_hbanded;               /* TRUE to use CP9 HMM bands for alignment        */
  double hbandp;                /* tail loss prob for HMM bands                   */
  int    use_sums;              /* TRUE to fill and use the posterior sums, false not to. */

  
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
  use_sums       = FALSE;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
		&optind, &optname, &optarg))  {
    if      (strcmp(optname, "-n") == 0) nseq         = atoi(optarg);
    else if (strcmp(optname, "-s") == 0) seed           = atoi(optarg);
    else if (strcmp(optname, "--sub")       == 0) do_sub = TRUE;
    else if (strcmp(optname, "--fsub")      == 0) do_fullsub = TRUE;
    else if (strcmp(optname, "--global")    == 0) do_local = FALSE;
    else if (strcmp(optname, "--post")      == 0) pthresh  = atof(optarg); 
    else if (strcmp(optname, "--fixlen")    == 0) { do_fixlen = TRUE; fixlen   = atoi(optarg); }
    else if (strcmp(optname, "--minlen")    == 0) { do_minlen = TRUE; minlen   = atoi(optarg); }
    else if (strcmp(optname, "--hbanded")   == 0) { do_hbanded = TRUE; do_small = FALSE; }
    else if (strcmp(optname, "--hbandp")    == 0) hbandp       = atof(optarg);
    else if (strcmp(optname, "--sums")      == 0) use_sums     = TRUE;
    else if (strcmp(optname, "--debug")     == 0) debug_level = atoi(optarg);
    else if (strcmp(optname, "-h") == 0) {
      MainBanner(stdout, banner);
      puts(usage);
      puts(experts);
      exit(EXIT_SUCCESS);
    }
  }
  
  if (argc - optind != 1) Die("Incorrect number of arguments.\n%s\n", usage);
  cmfile = argv[optind++];
  
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
  
  /*********************************************** 
   * Do the work.
   * Emit one sequence at a time, align it to the CM, IFF the optimal parsetree
   * is identical to the emitted parsetree, truncate it and realign it.
   ***********************************************/
  
  Parsetree_t **tr;             /* Parsetrees of emitted sequence         */
  Parsetree_t **ptr;            /* Parsetrees of partial emitted sequence */
  char        **seq;            /* actual sequence                        */
  char        **dsq;            /* digitized sequences                    */
  SQINFO       *sqinfo;         /* info about sequences (name/desc)       */
  char        **temp_dsq;       /* test digitized sequence                */
  SQINFO       *temp_sqinfo;    /* info about test sequence (name/desc)   */
  Parsetree_t **temp_tr;        /* Parsetrees of test emitted sequence    */

  SQINFO       *psqinfo;    /* info about test sequence (name/desc)   */

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

  
  dsq    = MallocOrDie(sizeof(char *)      * nseq);
  seq    = MallocOrDie(sizeof(char *)      * nseq);
  tr     = MallocOrDie(sizeof(Parsetree_t) * nseq);
  i = 0;
  sqinfo = MallocOrDie(sizeof(SQINFO)      * nseq);

  pdsq    = MallocOrDie(sizeof(char *)      * nseq);
  psqinfo = MallocOrDie(sizeof(SQINFO)      * nseq);
  ptr     = MallocOrDie(sizeof(Parsetree_t) * nseq);

  temp_sqinfo = MallocOrDie(sizeof(SQINFO)      * 1);
  temp_dsq    = MallocOrDie(sizeof(char *)      * 1);
  temp_tr     = MallocOrDie(sizeof(Parsetree_t) * 1);
  
  emap = CreateEmitMap(cm);

  spos = MallocOrDie(sizeof(int) * nseq);
  epos = MallocOrDie(sizeof(int) * nseq);

  /* Emit nseq seqs from the CM, align them to the CM using
   * the desired CYK algorithm (potentially banded), and truncate them
   * to partial sequences. Then realign with the same CYK algorithm
   * and determine how close the partial seq alignment is with 
   * the corresponding subalignment of the full sequence.
   */
  for(i = 0; i < nseq; i++)
    {
      if(debug_level > 0) printf("i: %d\n", i);
      attempts++;
      /* Emit a sequence from the CM, we don't care about the parsetree,
       * we'll determine this via alignment */
      EmitParsetree(cm, NULL, &(seq[i]), &(dsq[i]), &L);
      sqinfo[i].len   = L;
      sqinfo[i].flags = SQINFO_NAME | SQINFO_LEN;
      
      temp_tr[0] = tr[i];
      temp_dsq[0] = dsq[i];
      temp_sqinfo[0].len   = L;
      temp_sqinfo[0].flags = SQINFO_NAME | SQINFO_LEN;

      /* Align the sequence, setting do_sub to FALSE */
      AlignSeqsWrapper(cm, temp_dsq, temp_sqinfo, 1, &temp_tr, do_local, 
		       do_small, do_qdb, 0.0000001, do_hbanded, use_sums, 
		       hbandp, 
		       FALSE, FALSE,  /* these are do_sub and do_fullsub */
		       FALSE, FALSE, FALSE, 
		       FALSE, FALSE, NULL, FALSE, 0, 0, TRUE);


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
      pdsq[i] = DigitizeSequence(pseq[0], Lp);
      psqinfo[i].len   = Lp;
      strcpy(psqinfo[i].name, cm->name);
      psqinfo[i].flags = SQINFO_NAME | SQINFO_LEN;
    }

  /* Align all the partial sequences to the CM */
  AlignSeqsWrapper(cm, pdsq, psqinfo, nseq, &ptr, do_local, do_small, do_qdb,
		   0.0000001, do_hbanded, use_sums, hbandp,
		   do_sub, do_fullsub, FALSE, FALSE, FALSE, 
		   FALSE, FALSE, NULL, FALSE, 0, 0, TRUE);

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
      if(pred_spos != spos[i] || pred_epos != epos[i])
	printf("(i:%4d) S: %4d (A) %4d (P) E: %4d (A) %4d (P)\n", i, spos[i], pred_spos, epos[i], pred_epos);
      
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
  printf("n: %d | s_ct: %d | e_ct: %d\n", nseq, s_ct, e_ct);
  /*free(sqinfo);*/
}
  

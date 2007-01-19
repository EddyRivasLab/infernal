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
#include <ctype.h>
#include <float.h>

#include "squid.h"
#include "easel.h"
#include "esl_vectorops.h"
#include "esl_random.h"
#include "sqfuncs.h"
#include "dirichlet.h"
#include "sre_stack.h"

#include "easel.h"
#include "structs.h"
#include "funcs.h"
#include "stopwatch.h"          /* squid's process timing module        */
#include "hmmband.h"
#include "cm_postprob.h"

static void
pt_AlignSeqsWrapper(CM_t *cm, char **dsq, SQINFO *sqinfo, int nseq, Parsetree_t ***ret_tr, 
		    float fsub_pmass, int bdump_level, int debug_level, int silent_mode, 
		    int *actual_spos, int *actual_epos, float **ret_post_spos, 
		    float **ret_post_epos, int **ret_dist_spos, int **ret_dist_epos);

static char banner[] = "partial-test - test alignment of partial sequences";

static char usage[] = "\
Usage: partial-test [-options] <cmfile>\n\
  where options are:\n\
  -h     : help; print brief help on version and usage\n\
  -n <n> : number of seqs to emit, truncate and align [default: 100]\n\
  -s <n> : set random number seed to <n> \n\
";

static char experts[] = "\
  --read <f>   : read seqs to truncate from file <f>\n\
  --distro <f> : truncate based on inferred start/ends of seqs in file <f>\n\
  --sub        : aln w/sub CM for columns b/t HMM predicted start/end points\n\
  --fsub       : aln w/sub CM for structure b/t HMM predicted start/end points\n\
  --cp9        : aln w/CM plan 9 HMM\n\
  --global     : run alignment in global mode [default: local]\n\
  --post <f>   : minimum posterior prob from CP9 F/B to include in sub CM\n\
  --fixlen <n> : fix length of partial seqs to <n>\n\
  --minlen <n> : set minimum length of partial seqs as <n>\n\
  --debug  <n> : set verbosity of debugging print statements to <n> [1..3]\n\
  --histo  <n> : build histogram of HMM posterior probability of start/end\n\
  --repeat <n> : randomly truncate each seq <n> times (default 1)\n\
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
  { "--read",      FALSE, sqdARG_STRING},
  { "--distro",    FALSE, sqdARG_STRING},
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
  { "--histo",     FALSE, sqdARG_NONE},
  { "--repeat",    FALSE, sqdARG_INT}

};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv)
{
  char    *cmfile;		/* file to read CM from */	
  CMFILE  *cmfp;		/* open CM file for reading */
  CM_t    *cm;			/* a covariance model       */
  int      do_local;
  int      nseq;                /* number of seqs to aln */
  int      nrepeats;            /* number of times to truncate each seq */
  int      v;			/* counter over states */
  int      ncols;               /* number of consensus (match) columns in CM */
  int      i;                   /* counter over sub CMs */
  int      temp;

  double   pthresh;		

  ESL_RANDOMNESS *r;            /* Easel's random object */
  long     seed;		/* random number seed */

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
  int do_hbanded;               /* TRUE to use CP9 HMM bands for alignment        */
  double hbandp;                /* tail loss prob for HMM bands                   */
  int    use_sums;              /* TRUE to fill and use the posterior sums, false not to. */
  double beta;                  /* tail loss prob for QDB                         */

  /* --read related variables */
  int              do_read;     /* TRUE to read seqs from a file, instead of emitting */
  char            *read_seqfile;/* file to read sequences from */
  int              format;      /* format of sequence file */

  /* --distro related variables */
  int              do_distro;   /* TRUE to read seqs from a file, and build probability
				 * distro of start/ends based on their inferred start/ends */
  char            *distro_seqfile; /* file to read sequences from */
  CP9_dpmatrix_t  *cp9_post;    /* growable DP matrix for CP9 posteriors     */
  int              distro_nseq; /* number of seqs to aln */
  char           **distro_seq;  /* actual sequence                        */
  char           **distro_dsq;  /* digitized sequences                    */
  SQINFO          *distro_sqinfo; /* info about sequences (name/desc)       */
  float           *sdistro;     /* distribution of inferred starts from seqs in distro_seqfile */
  float           *edistro;     /* distribution of inferred ends   from seqs in distro_seqfile */
  int              distro_spos; /* current inferred start pos for seq from distro_seqfile */
  int              distro_epos; /* current inferred end pos for seq from distro_seqfile */
  int              distro_spos_state; /* current inferred state type (0=match, 1=insert) */
  int              distro_epos_state; /* current inferred state type (0=match, 1=insert) */
  float            swentry, swexit;  /* for configuring the CP9 HMM for local aln */

  Parsetree_t **tr;             /* Parsetrees of emitted sequence         */
  Parsetree_t **ptr;            /* Parsetrees of partial emitted sequence */
  char        **seq;            /* actual sequence                        */
  char        **dsq;            /* digitized sequences                    */
  SQINFO       *sqinfo;         /* info about sequences (name/desc)       */
  char        **temp_dsq;       /* test digitized sequence                */
  SQINFO       *temp_sqinfo;    /* info about test sequence (name/desc)   */
  Parsetree_t **temp_tr;        /* Parsetrees of test emitted sequence    */
  
  SQINFO       *psqinfo;        /* info about test sequence (name/desc)   */
  
  MSA               *msa;       /* alignment */
  int *useme;
  int apos;
  int cc;
  int *spos;
  int *epos;
  char **pseq;		/* dealigned seqs after truncation           */
  char **pdsq;		/* partial digitized sequences               */
  
  CMEmitMap_t *emap; 
  
  float *emit_sdistro;
  float *emit_edistro;
  int x;

  /*********************************************** 
   * Parse command line
   ***********************************************/

  seed           = (int) time ((time_t *) NULL);
  pthresh        = 0.1;
  do_sub         = FALSE;
  do_fullsub     = FALSE;
  nseq           = 100;
  nrepeats       = 1;
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
  read_seqfile   = NULL;
  do_distro      = FALSE;
  distro_seqfile = NULL;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
		&optind, &optname, &optarg))  {
    if      (strcmp(optname, "-n") == 0) nseq         = atoi(optarg);
    else if (strcmp(optname, "-s") == 0) seed           = (long) atoi(optarg);
    else if (strcmp(optname, "--read")      == 0) { do_read = TRUE;   read_seqfile = optarg;   }
    else if (strcmp(optname, "--distro")    == 0) { do_distro = TRUE; distro_seqfile = optarg; }
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
    else if (strcmp(optname, "--repeat")    == 0) nrepeats = atoi(optarg);
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

  /* Update cm->config_opts and cm->align_opts based on command line options */
  if(do_local)        cm->config_opts |= CM_CONFIG_LOCAL;
  if(do_hbanded)      cm->align_opts  |= CM_ALIGN_HBANDED;
  if(use_sums)        cm->align_opts  |= CM_ALIGN_SUMS;
  if(do_sub)          cm->align_opts  |= CM_ALIGN_SUB;
  if(do_fullsub)      cm->align_opts  |= CM_ALIGN_FSUB;
  if(!do_small)       cm->align_opts  |= CM_ALIGN_NOSMALL;
  if(do_qdb)          
    { 
      cm->align_opts  |= CM_ALIGN_QDB;
      cm->config_opts |= CM_CONFIG_QDB;
    }

  /* Configure the CM for alignment based on cm->config_opts and cm->align_opts.
   * set local mode, make cp9 HMM, calculate QD bands etc. */
  ConfigCM(cm, NULL, NULL);

  /* Determine number of consensus columns modelled by CM */
  ncols = 0;
  for(v = 0; v <= cm->M; v++)
    {
      if(cm->stid[v] ==  MATP_MP) ncols += 2;
      else if(cm->stid[v] == MATL_ML || cm->stid[v] == MATR_MR)	ncols++;
    }

  /* check for incompatible/misused options */
  if(do_minlen && do_fixlen)
    Die("Pick either --fixlen or --minlen, not both.\n");
  if(do_minlen && (minlen <= 0 || minlen > ncols))
    Die("--minlen enabled, but non-sensical, must be >= 0 and <= %d (# match cols)\n", ncols);
  if(do_fixlen && (fixlen <= 0 || fixlen > ncols))
    Die("--fixlen enabled, but non-sensical, must be >= 0 and <= %d (# match cols)\n", ncols);

  /****************************************************** 
   * If do_distro, open and read sequences from distro file,
   * We'll build probability distributions of the predicted
   * start/ends of these sequences, from which to pick
   * start/ends for truncation events later.
   ******************************************************/
  if(do_distro)
    {
      /* read the sequences from the input file */
      format = SQFILE_FASTA;
      printf("opening file: %s\n", distro_seqfile);
      if (! ReadMultipleRseqs(distro_seqfile, format, &distro_seq, &distro_sqinfo, &distro_nseq))
	Die("Failed to read any sequences from file %s, expecting FASTA.", distro_seqfile);
      distro_dsq = MallocOrDie(sizeof(char *) * distro_nseq);

      /* allocate and initialize the distributions */
      sdistro    = MallocOrDie(sizeof(float) * (ncols));
      edistro    = MallocOrDie(sizeof(float) * (ncols));
      emit_sdistro = MallocOrDie(sizeof(int) * ncols);
      emit_edistro = MallocOrDie(sizeof(int) * ncols);
      for(i = 0; i < ncols; i++)
	{
	  sdistro[i] = 0.;
	  edistro[i] = 0.;
	  emit_sdistro[i] = 0.;
	  emit_edistro[i] = 0.;
	}
      /* build a CP9 for the CM if nec. (we may have already if --hbanded enabled) */
      if(!(cm->flags & CM_CP9))
	{
	  if(!build_cp9_hmm(cm, &(cm->cp9), &(cm->cp9map), FALSE, 0.0001, 0))
	    Die("Couldn't build a CP9 HMM from the CM\n");
	  cm->flags |= CM_CP9; /* raise the CP9 flag */
	  /* configure the CP9 HMM for local alignment, all positions equiprobable */
	  swentry= ((cm->cp9->M)-1.)/cm->cp9->M; /* all start pts equiprobable, including 1 */
	  swexit = ((cm->cp9->M)-1.)/cm->cp9->M; /* all end   pts equiprobable, including M */
	  CPlan9SWConfig(cm->cp9, swentry, swexit);
	  CP9Logoddsify(cm->cp9);
	}
      for (i = 0; i < distro_nseq; i++) 
	{
	  distro_dsq[i] = DigitizeSequence(distro_seq[i], distro_sqinfo[i].len);
	  CP9_seq2posteriors(cm, distro_dsq[i], 1, distro_sqinfo[i].len, &cp9_post, 0); 
	  /* infer the start and end HMM nodes (consensus cols) from posterior matrix.
	   * Remember: we're necessarily in local mode, the --sub option turns local mode on. 
	   */
	  CP9NodeForPosn(cm->cp9, 1, distro_sqinfo[i].len, 1, 
			 cp9_post, &distro_spos, &distro_spos_state, do_fullsub, 0., TRUE, 0);
	  CP9NodeForPosn(cm->cp9, 1, distro_sqinfo[i].len, distro_sqinfo[i].len, 
			 cp9_post, &distro_epos, &distro_epos_state, do_fullsub, 0., TRUE, 0);
	  printf("i: %d ncols: %d S: %d E: %d\n", i, ncols, distro_spos, distro_epos);
	  if(distro_spos == 0 && distro_spos_state == 1) distro_spos = 1;
	  if(distro_epos == 0 && distro_epos_state == 1) distro_epos = 1;
	  sdistro[distro_spos-1] += 1.;
	  edistro[distro_epos-1] += 1.;
	}
      /* We've got all we need for the distribution, free distro_* data structures */
      for(i = 0; i < distro_nseq; i++)
	{
	  free(distro_seq[i]);
	  free(distro_dsq[i]);
	}
      free(distro_sqinfo);
      FreeCPlan9Matrix(cp9_post);

      /* Normalize the distros */
      esl_vec_FNorm(sdistro, ncols);
      esl_vec_FNorm(edistro, ncols);
      printf("\n\n");
    }      

  /*********************************************** 
   * Open output file, if needed.
   ***********************************************/

  /*   if (ofile == NULL) fp = stdout;
   else {
     if ((fp = fopen(ofile, "w")) == NULL)
       Die("Failed to open output file %s for writing", ofile);
       }*/
  
  /****************************************************** 
   * Do the work. Read or emit seqs and align them.
   ******************************************************/
  if(do_read)
    {
      /* read the sequences from the input file */
      format = SQFILE_FASTA;
      printf("opening file: %s\n", read_seqfile);
      if (! ReadMultipleRseqs(read_seqfile, format, &seq, &sqinfo, &nseq))
	Die("Failed to read any sequences from file %s, expecting FASTA.", read_seqfile);
      dsq = MallocOrDie(sizeof(char *) * nseq);
      for (i = 0; i < nseq; i++) 
	dsq[i] = DigitizeSequence(seq[i], sqinfo[i].len);
    }      

  /*********************************************** 
   * Show the options banner
   ***********************************************/

  MainBanner(stdout, banner);
  printf("CM file:                     %s\n", cmfile);
  printf("Number of seqs:              %d\n", nseq);
  printf("Number of truncs per seq:    %d\n", nrepeats);
  printf("Random seed:                 %ld\n", seed);
  printf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");

  /* Allocate and initialize */
  r = esl_randomness_Create(seed); 
  if(!do_read)
    {
      dsq    = MallocOrDie(sizeof(char *)      * (nseq));
      seq    = MallocOrDie(sizeof(char *)      * (nseq));
      sqinfo = MallocOrDie(sizeof(SQINFO)      * (nseq));
      tr     = MallocOrDie(sizeof(Parsetree_t) * (nseq));
    }
  i = 0;
  
  pdsq    = MallocOrDie(sizeof(char *)      * (nrepeats * nseq));
  psqinfo = MallocOrDie(sizeof(SQINFO)      * (nrepeats * nseq));
  ptr     = MallocOrDie(sizeof(Parsetree_t) * (nrepeats * nseq));
  
  temp_sqinfo = MallocOrDie(sizeof(SQINFO)      * 1);
  temp_dsq    = MallocOrDie(sizeof(char *)      * 1);
  temp_tr     = MallocOrDie(sizeof(Parsetree_t) * 1);

  spos = MallocOrDie(sizeof(int) * (nrepeats * nseq));
  epos = MallocOrDie(sizeof(int) * (nrepeats * nseq));
  
  emap = CreateEmitMap(cm);

  /* Either use read seqs (if do_read) or emit nseq seqs from the CM,
   * align them to the CM using the desired CYK algorithm (potentially
   * banded), and truncate them to partial sequences. Then realign
   * with the same CYK algorithm and determine how close the partial
   * seq alignment is with the corresponding subalignment of the full
   * sequence.
   */
  if(!(do_read))
    {
      for(i = 0; i < nseq; i++)
	{
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
	}
    }
  /* Align the sequences */
  /* turn off the do_sub option for initial alignment */
  printf("%-40s ... ", "Aligning full length sequences");
  if(do_sub) cm->align_opts &= ~CM_ALIGN_SUB;
  pt_AlignSeqsWrapper(cm, dsq, sqinfo, nseq, &tr, 
		      fsub_pmass, 0, 0, TRUE,
		      NULL, NULL, NULL, NULL, NULL, NULL);
  printf("done.\n");
  
  for(i = 0; i < (nrepeats * nseq); i++)
    {
      strcpy(temp_sqinfo[0].name, sqinfo[(i%nseq)].name);
      temp_sqinfo[0].len   = sqinfo[(i%nseq)].len;
      temp_tr[0]  = tr[(i%nseq)];
      temp_dsq[0] = dsq[(i%nseq)];
      temp_sqinfo[0].flags = SQINFO_NAME | SQINFO_LEN;
	  
      strcpy(temp_sqinfo[0].name, cm->name);
      msa = Parsetrees2Alignment(cm, temp_dsq, temp_sqinfo, NULL, temp_tr, 1, TRUE);
      /* important, the final variable do_full must be set to TRUE, we want
       * all match columns in the alignment.  */
      
      /* Step 1: pick random start and end consensus position
       *         for truncation, or if(do_distro) pick from
       *         the empirical distribution. */
      passed = FALSE;
      while(!passed)
	{
	  if(do_distro)
	    {
	      spos[i] = esl_rnd_FChoose(r, sdistro, ncols) + 1;
	      epos[i] = esl_rnd_FChoose(r, edistro, ncols) + 1;
	      emit_sdistro[(spos[i]-1)] += 1.;
	      emit_edistro[(epos[i]-1)] += 1.;
	    }
	  else
	    {
	      spos[i] = ((int) (esl_random(r) * ncols)) + 1;
	      epos[i] = ((int) (esl_random(r) * ncols)) + 1;
	    }
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
	    {
	      if((spos[i] + fixlen - 1) > ncols)
		passed = FALSE;
	      else
		epos[i] = spos[i] + fixlen - 1;
	    }	  
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
  /*****************************************************
   * Align the partial seqs to the CM and collect stats 
   * on how often we correctly get spos and epos.
   *****************************************************/

  /* Turn do_sub option back on if nec. */
  if(do_sub) cm->align_opts |= CM_ALIGN_SUB;
  printf("do_sub: %d\n", (cm->align_opts & CM_ALIGN_SUB));

  /* Align all the partial sequences to the CM */
  if(do_histo)
    pt_AlignSeqsWrapper(cm, pdsq, psqinfo, (nseq*nrepeats), &ptr, fsub_pmass, 0, 0, TRUE, 
			spos, epos, &post_spos, &post_epos, &dist_spos, &dist_epos);
  else
    pt_AlignSeqsWrapper(cm, pdsq, psqinfo, (nseq*nrepeats), &ptr, fsub_pmass, 0, 0, TRUE,
			NULL, NULL, NULL, NULL, NULL, NULL);
  s_ct = 0;
  e_ct = 0;
  /* For each sequence, compare the partial alignment with
   * the full alignment */
  for(i = 0; i < (nrepeats * nseq); i++)
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
		{
		  if(ptr[i]->state[x] == 2) /* ROOT_IR */
		    pred_epos = emap->rpos[cm->ndidx[ptr[i]->state[x]]] - 1;
		  else
		    pred_epos = 0; /* ROOT_IL */
		}
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
  if(do_distro)
    {
      /* Normalize the distros */
      esl_vec_FNorm(emit_sdistro, ncols);
      esl_vec_FNorm(emit_edistro, ncols);
      for(i = 0; i < ncols; i++)
	printf("Start[%4d]: %.4f %.4f\n", i, sdistro[i], emit_sdistro[i]);
      printf("\n\n");
      for(i = 0; i < ncols; i++)
	printf("End  [%4d]: %.4f %.4f\n", i, edistro[i], emit_edistro[i]);
    }
  /* Clean up. */
  if(do_distro)
    {
      free(sdistro);
      free(edistro);
      free(emit_sdistro);
      free(emit_edistro);
    }

  FreeCM(cm);
  /*if (tr[0] != NULL) FreeParsetree(tr[0]);  
    if (ptr[i] != NULL) FreeParsetree(ptr[i]); 
    free(dsq[0]);*/
  printf("N %d S %d E %d\n", nseq, s_ct, e_ct);
  /*free(sqinfo);*/
  for (i = 0; i < nseq; i++) 
    {
      free(dsq[i]);
      free(seq[i]);
      FreeParsetree(tr[i]);
    }
  for (i = 0; i < (nrepeats*nseq); i++) 
    {
      free(pdsq[i]);
      FreeParsetree(ptr[i]);
    }
  free(sqinfo);
  free(psqinfo);
  free(spos);
  free(epos);
  free(tr);
  free(ptr);
  free(pdsq);
  free(dsq);
  free(seq);
  FreeEmitMap(emap);
  return EXIT_SUCCESS;
}
  
/* EPN, Tue Dec  5 14:25:02 2006
 * 
 * Function: AlignSeqsWrapper()
 * 
 * Purpose:  Given a CM, digitized sequences, and a slew of options, 
 *           do preliminaries, call the correct CYK function and return
 *           parsetrees. Potentially collect stats for partial-test.
 * 
 * Args:     CM           - the covariance model
 *           dsq          - digitized sequences to align
 *           sqinfo       - info on the seq's we're aligning
 *           nseq         - number of seqs we're aligning
 *           ret_tr       - RETURN: parsetrees (pass NULL if trace isn't wanted)
 *           fsub_pmass   - probability mass to require in fullsub mode 
 *           bdump_level  - verbosity level for band related print statements
 *           debug_level  - verbosity level for debugging print statements
 *           silent_mode  - TRUE to not print anything, FALSE to print scores 
 *
 *           actual_spos  - [0..nseq-1] start consensus posn for truncated (partial) seq
 *           actual_epos  - [0..nseq-1] end   consensus posn for truncated (partial) seq
 *           ret_post_spos- [0..nseq-1] posterior probability from HMM of spos being start
 *           ret_post_epos- [0..nseq-1] posterior probability from HMM of epos being end
 *           ret_dist_spos- [0..nseq-1] distance (+/-) of max post start from spos
 *           ret_dist_epos- [0..nseq-1] distance (+/-) of max post end   from epos
 */
void
pt_AlignSeqsWrapper(CM_t *cm, char **dsq, SQINFO *sqinfo, int nseq, Parsetree_t ***ret_tr, 
		    float fsub_pmass, int bdump_level, int debug_level, int silent_mode, 
		    int *actual_spos, int *actual_epos, float **ret_post_spos, 
		    float **ret_post_epos, int **ret_dist_spos, int **ret_dist_epos)
{
  Stopwatch_t  *watch1, *watch2;      /* for timings */
  int i;                              /* counter over sequences */
  int v;                              /* state counter */
  Parsetree_t    **tr;          /* parse trees for the sequences */
  float            sc;		/* score for one sequence alignment */
  float            maxsc;	/* max score in all seqs */
  float            minsc;	/* min score in all seqs */
  float            avgsc;	/* avg score over all seqs */

  /* variables related to CM Plan 9 HMMs */
  struct cplan9_s       *hmm;           /* constructed CP9 HMM */
  CP9Bands_t            *cp9b;          /* data structure for hmm bands (bands on the hmm states) 
				         * and arrays for CM state bands, derived from HMM bands*/
  CP9Map_t              *cp9map;        /* maps the hmm to the cm and vice versa */
  struct cp9_dpmatrix_s *cp9_post;      /* growable DP matrix for posterior decode              */

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
  CP9Bands_t      *orig_cp9b; 

  /* variables related to query dependent banding (qdb) */
  int    expand_flag;           /* TRUE if the dmin and dmax vectors have just been 
				 * expanded (in which case we want to recalculate them 
				 * before we align a new sequence), and FALSE if not*/
  int     *dmin;                /* minimum d bound for state v, [0..v..M-1] */
  int     *dmax;                /* maximum d bound for state v, [0..v..M-1] */
  int *orig_dmin;               /* original dmin values passed in */
  int *orig_dmax;               /* original dmax values passed in */

  /* variables related to inside/outside */
  /*float           ***alpha;*/     /* alpha DP matrix for Inside() */
  /*float           ***beta; */     /* beta DP matrix for Inside() */
  /*float           ***post; */     /* post DP matrix for Inside() */
  int             ***alpha;    /* alpha DP matrix for Inside() */


  /* partial-test variables */
  int do_ptest;                 /* TRUE to fill partial-test variables */
  float *post_spos;
  float *post_epos;
  int   *dist_spos;
  int   *dist_epos;

  int do_hmmonly = FALSE;
  int do_elsilent = FALSE;
  int do_timings = FALSE;
  int do_check = FALSE;
  int do_post = FALSE;
  
  int do_local    = cm->config_opts  & CM_CONFIG_LOCAL;
  int do_hbanded  = cm->align_opts   & CM_ALIGN_HBANDED;
  int do_sub      = cm->align_opts   & CM_ALIGN_SUB;
  int do_fullsub  = cm->align_opts   & CM_ALIGN_FSUB;
  int do_qdb      = cm->align_opts   & CM_ALIGN_QDB;
  int do_inside   = cm->align_opts   & CM_ALIGN_INSIDE;
  int do_outside  = cm->align_opts   & CM_ALIGN_OUTSIDE;
  int do_small    = !(cm->align_opts & CM_ALIGN_NOSMALL);
  
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
      /* Keep this data for the original CM safe; we'll be doing
       * pointer swapping to ease the sub_cm alignment implementation. */
      hmm         = cm->cp9;
      cp9map      = cm->cp9map;

      orig_hmm = hmm;
      orig_cp9map = cp9map;
      if(do_hbanded)
	cp9b = AllocCP9Bands(cm, hmm);

      StopwatchZero(watch2);
      StopwatchStart(watch2);
    }

  /* Copy the QD bands in case we expand them. */
  if(do_qdb)
    {
      if(bdump_level > 1) 
	  /*printf("cm->beta:%f\n", cm->beta);*/
	  debug_print_bands(cm, cm->dmin, cm->dmax);
      expand_flag = FALSE;
      /* Copy dmin and dmax, so we can replace them after expansion */
      orig_dmin = MallocOrDie(sizeof(int) * cm->M);
      orig_dmax = MallocOrDie(sizeof(int) * cm->M);
      for(v = 0; v < cm->M; v++)
	{
	  orig_dmin[v] = cm->dmin[v];
	  orig_dmax[v] = cm->dmax[v];
	}
    }	  

  if(do_elsilent) 
    ConfigLocal_DisallowELEmissions(cm);

  if(do_hbanded)
    {
      cp9b = AllocCP9Bands(cm, cm->cp9);
      orig_cp9b = cp9b; 
    }
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
      if(do_hbanded)
	{
	  if(do_sub || do_ptest)
	    CP9_seq2bands(orig_cm, dsq[i], 1, sqinfo[i].len, orig_cp9b, 
			  &cp9_post, /* we DO want the posterior matrix back */
			  debug_level);
	  else
	    CP9_seq2bands(orig_cm, dsq[i], 1, sqinfo[i].len, orig_cp9b, 
			  NULL, /* we don't want the posterior matrix back */
			  debug_level);
	}
      if(do_sub && !do_hbanded)
	{
	  /* (1) Get HMM posteriors (if do_hbanded, we already have them) */
	  CP9_seq2posteriors(orig_cm, dsq[i], 1, sqinfo[i].len, &cp9_post, debug_level); 
	}
      if(do_ptest) /* determine the posterior probability from HMM of the correct start posn
		    * and end posn, as well as distance from max posterior */
	{
	  /* Determine HMM post probability of actual start and end */
	  /*	  CP9NodeForPosn(orig_hmm, 1, sqinfo[i].len, actual_spos[i], cp9_post, 
			 &spos, &spos_state, &(post_spos[i]), debug_level);
	  CP9NodeForPosn(orig_hmm, 1, sqinfo[i].len, actual_epos[i], cp9_post, 
	  &epos, &spos_state, &(post_epos[i]), debug_level);*/
	  printf("(%4d) s: %3d post: %.2f\n", i, actual_spos[i], 
		 Score2Prob(cp9_post->mmx[1][actual_spos[i]], 1.));
	  printf("(%4d) e: %3d post: %.2f\n", i, actual_epos[i], 
		 Score2Prob(cp9_post->mmx[sqinfo[i].len][actual_epos[i]], 1.));
	  /* Determine HMM post probability of most likely start and end */
	  /*CP9NodeForPosn(orig_hmm, 1, sqinfo[i].len, 1, cp9_post, 
			 &spos, &spos_state, NULL, debug_level);
	  CP9NodeForPosn(orig_hmm, 1, sqinfo[i].len, sqinfo[i].len, cp9_post, 
	  &epos, &epos_state, NULL, debug_level);*/
	}
      
      /* If we're in sub mode:
       * (1) Get HMM posteriors (we've already done this 
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
	  CP9NodeForPosn(orig_hmm, 1, sqinfo[i].len, 1,             cp9_post, &spos, &spos_state, 
			 do_fullsub, fsub_pmass, TRUE, debug_level);
	  CP9NodeForPosn(orig_hmm, 1, sqinfo[i].len, sqinfo[i].len, cp9_post, &epos, &epos_state, 
			 do_fullsub, fsub_pmass, FALSE, debug_level);
	  /* If the most likely state to have emitted the first or last residue
	   * is the insert state in node 0, it only makes sense to start modelling
	   * at consensus column 1. */
	  if(spos == 0 && spos_state == 1) 
	      spos = 1;
	  if(epos == 0 && epos_state == 1) 
	      epos = 1;
	  if(epos < spos) /* This is a possible but hopefully rarely encountered situation, we
			   * build a sub_cm identical to the CM, be setting spos = 1, epos = L */
	    {
	      spos = 1;
	      epos = sqinfo[i].len;
	    }
	  /* (3) Build the sub_cm from the original CM. */
	  if(!(build_sub_cm(orig_cm, &sub_cm, 
			    spos, epos,         /* first and last col of structure kept in the sub_cm  */
			    &submap,            /* maps from the sub_cm to cm and vice versa           */
			    do_fullsub,         /* build or not build a sub CM that models all columns */
			    debug_level)))      /* print or don't print debugging info                 */
	    Die("Couldn't build a sub CM from the CM\n");
	  /* Configure the sub_cm, the same as the cm, this will build a CP9 HMM if (do_hbanded) */
	  ConfigCM(sub_cm, NULL, NULL);
	  cm    = sub_cm; /* orig_cm still points to the original CM */

	  if(do_hbanded) /* we're doing HMM banded alignment to the sub_cm */
	    {
	      /* Get the HMM bands for the sub_cm */
	      sub_hmm = sub_cm->cp9;
	      sub_cp9map = sub_cm->cp9map;
	      sub_cp9b   = AllocCP9Bands(sub_cm, sub_cm->cp9);
	      CP9_seq2bands(sub_cm, dsq[i], 1, sqinfo[i].len, sub_cp9b, 
			    NULL, /* we don't want the posterior matrix back */
			    debug_level);
	      hmm           = sub_hmm;    
	      cp9map        = sub_cp9map;
	      cp9b          = sub_cp9b;
	    }
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
	  if(!(sub_cm2cm_parsetree(orig_cm, sub_cm, &orig_tr, tr[i], submap, (cm->align_opts & CM_ALIGN_FSUB), debug_level)))
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
	  FreeCPlan9Matrix(cp9_post);
	  FreeCM(sub_cm); /* cm and sub_cm now point to NULL */
	}
    }
  if(do_hbanded)
    FreeCP9Bands(cp9b);
  if (do_qdb)
    {
      free(orig_dmin);
      free(orig_dmax);
    }
  StopwatchFree(watch1);
  StopwatchFree(watch2);
  
  *ret_tr = tr; 
  
  if(ret_post_spos != NULL)
    *ret_post_spos = post_spos;
  if(ret_post_epos != NULL)
    *ret_post_epos = post_epos;
  if(ret_dist_spos != NULL)
    *ret_dist_spos = dist_spos;
  if(ret_dist_epos != NULL)
    *ret_dist_epos = dist_epos;

  /*printf("leaving AlignSeqsWrapper()\n");*/
}

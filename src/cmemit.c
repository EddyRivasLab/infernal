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
#include "squid.h"		/* general sequence analysis library    */
#include "msa.h"		/* squid's multiple sequence i/o        */
#include "stats.h"
#include "cm_dispatch.h"	
#include "sre_random.h"
#include <esl_vectorops.h>
#include <esl_histogram.h>
#include <esl_stats.h>
#include <esl_random.h>

static char banner[] = "cmemit - generate sequences from a covariance model";

static char usage[]  = "\
Usage: cmemit [-options] <cm file>\n\
Available options are:\n\
   -l     : local; emit from a locally configured CM\n\
   -s <n> : set random number seed to <n>\n\
   -a     : write generated sequences as an alignment, not FASTA\n\
   -c     : generate a single \"consensus\" sequence\n\
   -h     : help; print brief help on version and usage\n\
   -n <n> : emit <n> sequences (default 10)\n\
   -o <f> : save sequences in file <f>\n\
   -q     : quiet - suppress verbose banner\n\
";

static char experts[] = "\
   --begin <n>    : truncate alignment, begin at match column <n>\n\
   --end   <n>    : truncate alignment, end   at match column <n>\n\
   --zeroinserts  : set insert emission scores to 0\n\
   --hmmbuild     : build a ML CM Plan 9 HMM from the samples\n\
   --hmmscore     : score samples with a CM Plan 9 HMM\n\
   --hmmlocal     : w/hmmscore, search with CP9 HMM in local mode\n\
   --hmmfast      : w/hmmscore, assume parse tree scores are optimal\n\
   --cmeval       : w/hmmscore, min CM E-value to consider\n\
   --cmN <n>      : w/hmmscore & cmeval, db size E-val corresponds with\n\
   --inside       : w/hmmscore, score CM seqs with inside\n\
   --ptest        : do parse test, are generated parses optimal?\n\
   --exp <x>      : exponentiate CM probs by <x> prior to emitting\n\
";

static struct opt_s OPTIONS[] = {
  { "-h",        TRUE,  sqdARG_NONE }, 
  { "-s",        TRUE,  sqdARG_INT},
  { "-l",        TRUE,  sqdARG_NONE },
  { "-a",        TRUE,  sqdARG_NONE },  
  { "-c",        TRUE,  sqdARG_NONE },  
  { "-n",        TRUE,  sqdARG_INT},  
  { "-o",        TRUE,  sqdARG_STRING},
  { "-q",        TRUE,  sqdARG_NONE},  
  { "--begin",   FALSE, sqdARG_INT },
  { "--end",     FALSE, sqdARG_INT },
  { "--zeroinserts",FALSE, sqdARG_NONE},
  { "--hmmbuild",FALSE, sqdARG_NONE },
  { "--hmmscore",FALSE, sqdARG_NONE },
  { "--hmmlocal",FALSE, sqdARG_NONE },
  { "--hmmfast", FALSE, sqdARG_NONE },
  { "--cmeval",  FALSE, sqdARG_FLOAT },
  { "--cmN",     FALSE, sqdARG_FLOAT },
  { "--ptest",   FALSE, sqdARG_NONE },
  { "--exp",     FALSE, sqdARG_FLOAT },
  { "--inside",  FALSE, sqdARG_NONE }
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv) 
{
  char            *cmfile;      /* file to read CM from */	
  CMFILE          *cmfp;	/* open CM file */
  ESL_RANDOMNESS  *r = NULL;    /* source of randomness                    */
  CM_t            *cm;          /* CM to generate from */

  FILE            *fp;          /* output file handle                      */
  char            *ofile;       /* output sequence file                    */
  int              L;		/* length of a sequence                    */
  int              i;		/* counter over sequences                  */

  int              nseq;	/* number of seqs to sample                */
  long             seed;	/* random number generator seed            */
  int              seed_set;	/* TRUE if -s set on command line      */
  int              do_qdb;	/* TRUE to config QDBs (on by default)     */
  int              do_local;	/* TRUE to config the model in local mode  */
  int              be_quiet;	/* TRUE to silence header/footer           */
  int              do_alignment;/* TRUE to output in aligned format        */ 
  int              do_consensus;/* TRUE to do a single consensus seq       */

  char *optname;                /* name of option found by Getopt()        */
  char *optarg;                 /* argument found by Getopt()              */
  int   optind;                 /* index in argv[]                         */
  int   begin_set;              /* TRUE if --begin entered at command line */
  int   end_set;                /* TRUE if --end entered at command line   */
  int   spos;                   /* start match column for MSA              */
  int   epos;                   /* end match column for MSA                */
  int   ncols;                  /* number of match columns modelled by CM  */
  int   v;                      /* counter over states                     */
  int   do_zero_inserts;         /* TRUE to zero insert emission scores    */
  int   build_cp9;              /* TRUE to build a ML CP9 HMM from the 
				 * sampled parses of the CM                */
  int   do_score_cp9;           /* TRUE to score each seq against CP9 HMM  */
  int   do_hmmlocal;           /* if score_cp9, put CP9 in local mode     */
  int   do_fastfil;             /* TRUE to assume parse tree score is opt  */
  float cm_sc_cutoff = IMPOSSIBLE;  
  float cm_e_cutoff  = 10.;
  int   cm_cutoff_type = SCORE_CUTOFF; 
  long  cmN = 0;                /* if nec, db size associated with cm_e_cutoff */
  int   cmN_set;                /* TRUE if --cmN enabled                   */
  int   do_ptest;               /* TRUE:check if emitted parses are optimal*/
  int   do_exp;                 /* TRUE to exponentiate CM before emitting */
  double exp_factor;            /* factor to exponentiate CM params by     */
  int   do_inside;              /* TRUE to search emitted seqs w/inside    */
  /*********************************************** 
   * Parse command line
   ***********************************************/

  nseq         = 10;
  seed         = time ((time_t *) NULL);
  be_quiet     = FALSE;
  seed_set     = FALSE;
  do_qdb       = TRUE;  /* on by default */
  do_local     = FALSE;
  do_alignment = FALSE;  
  do_consensus = FALSE;
  begin_set    = FALSE;
  end_set      = FALSE;
  do_zero_inserts=FALSE;
  build_cp9    = FALSE;
  do_score_cp9 = FALSE;
  do_fastfil   = FALSE;
  do_hmmlocal = FALSE;
  do_ptest     = FALSE;
  do_exp       = FALSE;
  do_inside    = FALSE;
  cmN_set      = FALSE;
  ofile        = NULL;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-l")      == 0) do_local     = TRUE;
    else if (strcmp(optname, "-s")      == 0) { seed_set = TRUE; seed = (long) atoi(optarg); } 
    else if (strcmp(optname, "-a")      == 0) do_alignment = TRUE;
    else if (strcmp(optname, "-c")      == 0) do_consensus = TRUE;
    else if (strcmp(optname, "-n")      == 0) nseq         = atoi(optarg); 
    else if (strcmp(optname, "-o")      == 0) ofile        = optarg;
    else if (strcmp(optname, "-q")      == 0) be_quiet     = TRUE;
    else if (strcmp(optname, "--begin") == 0) { begin_set  = TRUE; spos = atoi(optarg); }
    else if (strcmp(optname, "--end")   == 0) { end_set    = TRUE; epos = atoi(optarg); }
    else if (strcmp(optname, "--zeroinserts")== 0) do_zero_inserts = TRUE;
    else if (strcmp(optname, "--hmmbuild") == 0) build_cp9 = TRUE;
    else if (strcmp(optname, "--hmmscore") == 0) do_score_cp9 = TRUE;
    else if (strcmp(optname, "--hmmlocal") == 0) do_hmmlocal = TRUE;
    else if (strcmp(optname, "--hmmfast")  == 0) do_fastfil   = TRUE;
    else if (strcmp(optname, "--cmeval")   == 0) { cm_cutoff_type = E_CUTOFF; cm_e_cutoff = atof(optarg); } 
    else if (strcmp(optname, "--cmN")      == 0) { cmN_set = TRUE; cmN = (long) atoi(optarg); } 
    else if (strcmp(optname, "--ptest")    == 0) do_ptest = TRUE; 
    else if (strcmp(optname, "--exp")      == 0) { do_exp = TRUE; exp_factor = atof(optarg); }
    else if (strcmp(optname, "--inside")   == 0) do_inside = TRUE; 
    else if (strcmp(optname, "-h")      == 0) 
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

  /* Check for incompatible option combos */
  if (do_alignment && do_consensus)
    Die("Sorry, -a and -c are incompatible.\nUsage:\n%s", usage); 
  if (nseq != 10 && do_consensus)
    Warn("-c (consensus) overrides -n (# of sampled seqs)");
  if((begin_set && !end_set) || (!begin_set && end_set))
    Die("Must use both --begin and --end or neither.\n");
  if(begin_set && (!do_alignment && !build_cp9))
    Die("--begin and --end only work with -a or --cp9\n");
  if(build_cp9 && do_alignment)
    Die("Sorry, --cp9 and -a are incompatible.\nUsage:\n%s", usage);
  if(do_hmmlocal && !do_score_cp9)
    Die("ERROR --hmmlocal only makes sense in combination with --hmmscore.\n");
  if(do_score_cp9 && (do_alignment || do_consensus))
    Die("ERROR --hmmscore does not work in combination with -a or -c.\n");
  if(do_fastfil && do_ptest)
    Die("ERROR --hmmfast does not work in combination with --ptest.\n");
  if((!cmN_set) && (cm_cutoff_type == E_CUTOFF))
    Die("ERROR --cmN and --cmeval must both be used if one is used.\n");
  
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

  if(do_qdb)          cm->config_opts |= CM_CONFIG_QDB;
  if(do_local)        cm->config_opts |= CM_CONFIG_LOCAL;
  if(do_hmmlocal)     cm->config_opts |= CM_CONFIG_HMMLOCAL;
  if(do_inside)       cm->search_opts |= CM_SEARCH_INSIDE;
  if(do_zero_inserts) cm->config_opts |= CM_CONFIG_ZEROINSERTS;
  ConfigCM(cm, NULL, NULL);

  /* Determine number of consensus columns modelled by CM */
  ncols = 0;
  for(v = 0; v <= cm->M; v++)
    {
      if(cm->stid[v] ==  MATP_MP)
	ncols += 2;
      else if(cm->stid[v] == MATL_ML || cm->stid[v] == MATR_MR)
	ncols++;
    }

  if(begin_set && end_set)
    {
      if(spos < 1) Die("ERROR, when using --begin <n>, <n> must be >= 1\n");
      if(epos > ncols) Die("ERROR, --end %d selected; there's only %d match columns.\n", epos, ncols);
    }


  /**********************
   * Creae and seed RNG *
   **********************/
  if (!(seed_set)) 
    seed = time ((time_t *) NULL);
  if ((r = esl_randomness_Create(seed)) == NULL) /* we want to know what seed is, this is why
						  * we don't use esl_randomness_CreateTimeseeded(),
						  * b/c we lose the seed in that function. */
    esl_fatal("Failed to create random number generator: probably out of memory");

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
	printf("Random seed:          %ld\n", seed);
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

  if (do_exp)
    ExponentiateCM(cm, exp_factor);

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
      int *useme;
      int apos;
      int cc;
      int *ct;		/* 0..alen-1 base pair partners array         */

      dsq    = MallocOrDie(sizeof(char *)      * nseq);
      tr     = MallocOrDie(sizeof(Parsetree_t) * nseq);
      sqinfo = MallocOrDie(sizeof(SQINFO)      * nseq);
      wgt    = MallocOrDie(sizeof(float)       * nseq);
      FSet(wgt, nseq, 1.0);
      
      for (i = 0; i < nseq; i++)
	{
	  EmitParsetree(cm, r, &(tr[i]), NULL, &(dsq[i]), &L);
	  sprintf(sqinfo[i].name, "seq%d", i+1);
	  sqinfo[i].len   = L;
	  sqinfo[i].flags = SQINFO_NAME | SQINFO_LEN;
	}
      
      msa = Parsetrees2Alignment(cm, dsq, sqinfo, NULL, tr, nseq, 
				 TRUE); /* we want all match columns in alignment */
      msa->name = sre_strdup(cm->name, -1);
      msa->desc = sre_strdup("Synthetic sequence alignment generated by cmemit", -1);
      
      if(begin_set && end_set)
	{
	  WUSS2ct(msa->ss_cons, msa->alen, FALSE, &ct);  
	  /* Truncate the alignment prior to consensus column spos and after 
	     consensus column epos */
	  useme = (int *) MallocOrDie (sizeof(int) * (msa->alen+1));
	  for (apos = 0, cc = 0; apos < msa->alen; apos++)
	    {
	      /* Careful here, placement of cc++ increment is impt, 
	       * we want all inserts between cc=spos-1 and cc=spos,
	       * and between cc=epos and cc=epos+1.
	       * Also be careful: ct[] is index 1..alen,
	       * and msa->ss_cons is 0..alen-1.
	       */
	      if(cc < (spos-1) || cc > epos)
		{
		  useme[apos] = 0;
		  if(ct[(apos+1)] != 0) ct[ct[(apos+1)]] = 0;
		  ct[(apos+1)] = 0;
		}
	      else
		{
		  useme[apos] = 1;
		}
	      if (!isgap(msa->rf[apos])) 
		{ 
		  cc++; 
		  if(cc == (epos+1))
		    {
		      useme[apos] = 0; 
		      /* we misassigned this guy, overwrite */ 
		      if(ct[(apos+1)] != 0) ct[ct[(apos+1)]] = 0;
		      ct[(apos+1)] = 0;
		    }
		}
	    }
	  /* construct the new structure based on the ct array,
	   * we don't do full WUSS notation (just laziness) */
	  for (apos = 0; apos < msa->alen; apos++)
	    {
	      if(ct[(apos+1)] == 0)
		{
		  if(msa->ss_cons[apos] != '.')
		    msa->ss_cons[apos] = ':';
		}
	      else if (ct[apos+1]  > (apos+1)) msa->ss_cons[apos] = '<';
	      else if (ct[apos+1]  < (apos+1)) msa->ss_cons[apos] = '>';
	    }	    
	  free(ct);
	  
	  /*rintf("\n\nDEBUG PRINTING ORIG ALIGNMENT:\n");
	    WriteStockholm(fp, msa);
	    printf("\n\nDONE DEBUG PRINTING ORIG ALIGNMENT:\n");
	    for(apos=0; apos < msa->alen; apos++)
	    printf("useme[%d]: %d\n", apos, useme[apos]);
	  */
	  
	  MSAShorterAlignment(msa, useme);
	  free(useme);
	}
      
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

  else if(build_cp9)
    {
      Parsetree_t **tr;             /* Parsetrees of emitted aligned sequences */
      char    **dsq;                /* digitized sequences                     */
      char    **seq;                /* actual sequences (real letters)         */
      SQINFO            *sqinfo;    /* info about sequences (name/desc)        */
      MSA               *msa;       /* alignment */
      float             *wgt;
      int apos;
      int cc;

      int *matassign;
      int *useme;
      int msa_nseq;                 /* this is the number of sequences per MSA,
				     * current strategy is to sample (nseq/nseq_per_msa)
				     * alignments from the CM, and add counts from
				     * each to the shmm in counts form (to limit memory)
				     */
      int nsampled;                 /* number of sequences sampled thus far */

      CP9_t  *shmm;
      CP9trace_t **cp9_tr;   /* fake tracebacks for each seq            */
      int idx;

      msa_nseq = 1000;
      /* Allocate and zero the new HMM we're going to build by sampling from
       * the CM.
       */
      if(begin_set && end_set)
	shmm = AllocCPlan9(epos - spos + 1);
      else
	shmm = AllocCPlan9(ncols);

      ZeroCPlan9(shmm);

      /* sample MSA(s) from the CM */
      nsampled = 0;
      dsq    = MallocOrDie(sizeof(char *)             * msa_nseq);
      seq    = MallocOrDie(sizeof(char *)             * msa_nseq);
      tr     = MallocOrDie(sizeof(Parsetree_t)        * msa_nseq);
      sqinfo = MallocOrDie(sizeof(SQINFO)             * msa_nseq);
      wgt    = MallocOrDie(sizeof(float)              * msa_nseq);
      FSet(wgt, msa_nseq, 1.0);
      
      while(nsampled < nseq)
	{
	  /*printf("nsampled: %d\n", nsampled);*/
	  if(nsampled != 0)
	    {
	      /* clean up from previous MSA */
	      MSAFree(msa);
	      free(matassign);
	      if(begin_set && end_set)
		{
		  free(useme);
		  for (i = 0; i < msa_nseq; i++)
		    free(seq[i]);
		}
	      for (i = 0; i < msa_nseq; i++)
		{
		  CP9FreeTrace(cp9_tr[i]);
		  FreeParsetree(tr[i]);
		  free(dsq[i]);
		}
	      free(cp9_tr);
	    }
	  /* Emit msa_nseq parsetrees from the CM */
	  if(nsampled + msa_nseq > nseq)
	    msa_nseq = nseq - nsampled;
	  for (i = 0; i < msa_nseq; i++)
	    {
	      EmitParsetree(cm, r, &(tr[i]), NULL, &(dsq[i]), &L);
	      sprintf(sqinfo[i].name, "seq%d", i+1);
	      sqinfo[i].len   = L;
	      sqinfo[i].flags = SQINFO_NAME | SQINFO_LEN;
	    }
	  /* Build a new MSA from these parsetrees */
	  msa = Parsetrees2Alignment(cm, dsq, sqinfo, NULL, tr, msa_nseq, TRUE);
	  
	  /* Add the counts to the growing counts-based HMM */

	  /* If necessary, truncate the alignment prior to consensus column spos 
	   * and after consensus column epos */
	  if(begin_set && end_set)
	    {
	      useme = (int *) MallocOrDie (sizeof(int) * (msa->alen+1));
	      for (apos = 0, cc = 0; apos < msa->alen; apos++)
		{
		  /* Careful here, placement of cc++ increment is impt, 
		   * we want all inserts between cc=spos-1 and cc=spos,
		   * and between cc=epos and cc=epos+1.
		   */
		  if(cc < (spos-1) || cc > epos)
		    useme[apos] = 0;
		  else
		    useme[apos] = 1;
		  if (!isgap(msa->rf[apos])) 
		    { 
		      cc++; 
		      if(cc == (epos+1))
			useme[apos] = 0; 
		      /* we misassigned this guy, overwrite */ 
		    }
		}
	      MSAShorterAlignment(msa, useme);
	      
	      /* Shorten the dsq's to match the shortened MSA */
	      for (i = 0; i < msa_nseq; i++)
		{
		  MakeDealignedString(msa->aseq[i], msa->alen, msa->aseq[i], &(seq[i])); 
		  free(dsq[i]);
		  dsq[i] = DigitizeSequence(seq[i], strlen(seq[i]));
		}
	    }
	  /* Determine match assignment from RF annotation
	   */
	  matassign = (int *) MallocOrDie (sizeof(int) * (msa->alen+1));
	  matassign[0] = 0;
	  for (apos = 0; apos < msa->alen; apos++)
	    {
	      matassign[apos+1] = 0;
	      if (!isgap(msa->rf[apos])) 
		matassign[apos+1] = 1;
	    }
	  /* make fake tracebacks for each seq */
	  CP9_fake_tracebacks(msa->aseq, msa->nseq, msa->alen, matassign, &cp9_tr);
	  
	  /* build model from tracebacks (code from HMMER's modelmakers.c::matassign2hmm() */
	  for (idx = 0; idx < msa->nseq; idx++) {
	    CP9TraceCount(shmm, dsq[idx], msa->wgt[idx], cp9_tr[idx]);
	  }
	  nsampled += msa_nseq;
	}
      
      /* clean up from final MSA */
      MSAFree(msa);
      free(matassign);
      if(begin_set && end_set)
	{
	  free(useme);
	  for (i = 0; i < msa_nseq; i++)
	    free(seq[i]);
	}
      for (i = 0; i < msa_nseq; i++)
	{
	  CP9FreeTrace(cp9_tr[i]);
	  FreeParsetree(tr[i]);
	  free(dsq[i]);
	}
      free(cp9_tr);

      printf("PRINTING NON-NORM SAMPLED HMM PARAMS:\n");
      debug_print_cp9_params(shmm);
      printf("DONE PRINTING NON-NORM SAMPLED HMM PARAMS:\n");

      CPlan9Renormalize(shmm);
      CP9Logoddsify(shmm);

      printf("PRINTING NORM SAMPLED HMM PARAMS:\n");
      debug_print_cp9_params(shmm);
      printf("DONE PRINTING NORM SAMPLED HMM PARAMS:\n");
      
      FreeCPlan9(shmm);
      free(seq);
      free(dsq);
    }      
  else if(do_score_cp9)				/* emit unaligned seqs and score 
						 * each with a CP9 HMM */
    {
      Parsetree_t      *tr;         /* generated trace                        */
      char             *dsq;        /* digitized sequence                     */
      char             *seq;        /* alphabetic sequence                    */
      float             cm_sc;
      float            *hmm_sc;
      int               max_attempts = 500 * nseq;
      int               nattempts = 0;
      float             f;
      float             hmm_cutoff = 0.;
      double           *xv;
      int               n;
      double            mean, var;
      float             diff = 0.;  /* total diff b/t optimal/emitted parse scs */
      int               ndiff = 0;
      float             esc;        /* emitted parse tree score               */
      SQINFO            sqinfo;     /* info about sequence (name/len)         */
      int               cm_mode;    /* EVD mode to use for CM */
      int               cp9_mode;   /* EVD mode to use for CP9 */
      float             tmp_K;

      SetCMCutoff(cm, cm_cutoff_type, cm_sc_cutoff, cm_e_cutoff);
      SetCP9Cutoff(cm, SCORE_CUTOFF, 0., 0., cm->cutoff);
      CM2EVD_mode(cm, &cm_mode, &cp9_mode); 

      if(cm->cutoff_type == E_CUTOFF)
	{
	  /* Working with first partition */
	  if(cm->stats->np > 1)
	    Die("ERROR CM EVD stats for > 1 partition, not yet implemented.\n");
	  /* First update mu based on --cmN argument */
	  tmp_K = exp(cm->stats->evdAA[cm_mode][0]->mu * cm->stats->evdAA[cm_mode][0]->lambda) / 
	    cm->stats->evdAA[cm_mode][0]->L;
	  cm->stats->evdAA[cm_mode][0]->mu = log(tmp_K * ((double) cmN)) /
	    cm->stats->evdAA[cm_mode][0]->lambda;
	  cm->stats->evdAA[cm_mode][0]->L = cmN; /* update L, the seq size stats correspond to */
	  cm_sc_cutoff = (cm->stats->evdAA[cm_mode][0]->mu - 
			  (log(cm->cutoff) / cm->stats->evdAA[cm_mode][0]->lambda));

	  tmp_K = exp(cm->stats->evdAA[cp9_mode][0]->mu * cm->stats->evdAA[cp9_mode][0]->lambda) / 
	    cm->stats->evdAA[cp9_mode][0]->L;
	  cm->stats->evdAA[cp9_mode][0]->mu = log(tmp_K * ((double) cmN)) /
	    cm->stats->evdAA[cp9_mode][0]->lambda;
	  cm->stats->evdAA[cp9_mode][0]->L = cmN; /* update L, the seq size stats correspond to */

	  PrintSearchInfo(stdout, cm, cm_mode, cp9_mode, cmN);

	}
      else
	{
	  cm_sc_cutoff = IMPOSSIBLE;
	  PrintSearchInfo(stdout, cm, cm_mode, cp9_mode, 0);
	}
      
      ESL_HISTOGRAM *h = esl_histogram_CreateFull(100, 1000, 0.2);
      hmm_sc = MallocOrDie(sizeof(float) * nseq);
      sre_srandom(33);

      printf("do fastfil: %d\n", do_fastfil);
      for (i = 0; i < nseq; i++)
	{
	  if(nattempts++ > max_attempts) 
	    Die("ERROR number of attempts exceeded 500 times number of seqs.\n");
	  EmitParsetree(cm, r, &tr, &seq, &dsq, &L);
	  /*sprintf(sqinfo.name, "seq%d", i+1);
	    sqinfo.len   = L;
	    sqinfo.flags = SQINFO_NAME | SQINFO_LEN;
	    printf("nattempts: %d\n", nattempts);
	    WriteSeq(stdout, SQFILE_FASTA, seq, &sqinfo);*/

	  cm->search_opts &= ~CM_SEARCH_HMMONLY;
	  if(do_fastfil) 
	    cm_sc = ParsetreeScore(cm, tr, dsq, FALSE);
	  else
	    cm_sc = actually_search_target(cm, dsq, 1, L,
					   0.,    /* cutoff is 0 bits (actually we'll find highest
						   * negative score if it's < 0.0) */
					   0.,    /* CP9 cutoff is 0 bits */
					   NULL,  /* don't keep results */
					   FALSE, /* don't filter with a CP9 HMM */
					   FALSE, /* we're not calcing CM  stats */
					   FALSE, /* we're not calcing CP9 stats */
					   NULL); /* filter fraction N/A */
	  while(cm_sc < cm_sc_cutoff)
	    {
	      FreeParsetree(tr);
	      free(dsq);
	      free(seq);
	      EmitParsetree(cm, r, &tr, &seq, &dsq, &L);

	      /*sprintf(sqinfo.name, "seq%d", i+1);
		sqinfo.len   = L;
		sqinfo.flags = SQINFO_NAME | SQINFO_LEN;
		printf("nattempts: %d\n", nattempts);
		WriteSeq(stdout, SQFILE_FASTA, seq, &sqinfo);*/

	      if(do_fastfil) 
		cm_sc = ParsetreeScore(cm, tr, dsq, FALSE);
	      else
		cm_sc = actually_search_target(cm, dsq, 1, L,
					       0.,    /* cutoff is 0 bits (actually we'll find highest
						       * negative score if it's < 0.0) */
					       0.,    /* CP9 cutoff is 0 bits */
					       NULL,  /* don't keep results */
					       FALSE, /* don't filter with a CP9 HMM */
					       FALSE, /* we're not calcing CM  stats */
					       FALSE, /* we're not calcing CP9 stats */
					       NULL); /* filter fraction N/A */
	      if(nattempts++ > max_attempts)
		Die("ERROR number of attempts exceeded 50 times number of seqs.\n");
	    }
	  esl_histogram_Add(h, cm_sc);
	  if(do_ptest) 
	    {
	      esc   = ParsetreeScore(cm, tr, dsq, FALSE);
	      diff += fabs(cm_sc - esc);
	    }
	  cm->search_opts |= CM_SEARCH_HMMONLY;
	  hmm_sc[i] = CP9Forward(cm, dsq, 1, L, cm->W, 0., 
				 NULL,   /* don't return scores of hits */
				 NULL,   /* don't return posns of hits */
				 NULL,   /* don't keep track of hits */
				 TRUE,   /* we're scanning */
				 FALSE,  /* we're not ultimately aligning */
				 FALSE,  /* we're not rescanning */
				 TRUE,   /* be memory efficient */
				 NULL);  /* don't want the DP matrix back */
	  FreeParsetree(tr);
	  free(dsq);
	  free(seq);


	  printf("i: %4d cm sc: %10.4f hmm sc: %10.4f ", i, cm_sc, hmm_sc[i]);
	  if(do_ptest) printf("D: %10.3f (E: %10.3f)\n", (cm_sc - esc), esc);
	  if(fabs(cm_sc-esc) > 0.0001) ndiff++;
	  else printf("\n");
	}
      if(do_ptest) printf("Summary: Num diff: %d Bit diff: %.3f Avg per seq: %.3f\n", ndiff, diff, (diff/nseq));
      /* Sort the HMM scores with quicksort */
      esl_vec_FSortIncreasing(hmm_sc, nseq);

      /*esl_histogram_Print(stdout, h);*/

      esl_histogram_GetData(h, &xv, &n);
      esl_stats_Mean(xv, n, &mean, &var);
      printf("N: %d\nMean:   %f\nVar:    %f\nSt dev: %f\n", n, mean, var, sqrt(var));

      esl_histogram_Destroy(h);
      
      f = 0.95;
      hmm_cutoff = hmm_sc[(int) ((1. - f) * (float) nseq)];
      printf("HMM glocal filter bit score threshold for finding %.2f CM hits > %.2f bits: %.4f\n", f, cm_sc_cutoff, hmm_cutoff);
      printf("\n\nnattempts: %d\n", nattempts);

      float eval = RJK_ExtremeValueE(hmm_cutoff, 
				     cm->stats->evdAA[cp9_mode][0]->mu,
				     cm->stats->evdAA[cp9_mode][0]->lambda);
      printf("05.21.07 %d %d %f %f\n", cm_mode, cp9_mode, hmm_cutoff, eval);
      free(hmm_sc);
    }
  else				/* unaligned sequence output */
    {
      Parsetree_t      *tr;         /* generated trace                        */
      char             *dsq;        /* digitized sequence                     */
      char             *seq;        /* alphabetic sequence                    */
      SQINFO            sqinfo;     /* info about sequence (name/len)         */
      float             diff = 0.;  /* total diff b/t optimal/emitted parse scs */
      float             esc;        /* emitted parse tree score               */
      float             osc;        /* optimal parse tree score               */
      int               ndiff = 0;
      for (i = 0; i < nseq; i++)
	{
	  EmitParsetree(cm, r, &tr, &seq, &dsq, &L);

	  sprintf(sqinfo.name, "%s-%d", cm->name, i+1);
	  sqinfo.len   = L;
	  sqinfo.flags = SQINFO_NAME | SQINFO_LEN;
	  WriteSeq(fp, SQFILE_FASTA, seq, &sqinfo);
	  
	  if(do_ptest) /* Check score diff, if any w/optimal parse */
	    {
	      
	      esc = ParsetreeScore(cm, tr, dsq, FALSE);
	      osc = actually_search_target(cm, dsq, 1, L,
					   0.,    /* cutoff is 0 bits (actually we'll find highest
						   * negative score if it's < 0.0) */
					   0.,    /* CP9 cutoff is 0 bits */
					   NULL,  /* don't keep results */
					   FALSE, /* don't filter with a CP9 HMM */
					   FALSE, /* we're not calcing CM  stats */
					   FALSE, /* we're not calcing CP9 stats */
					   NULL); /* filter fraction N/A */
	      diff += fabs(osc-esc);
	      printf("i: %5d %10.3f (E: %10.3f O: %10.3f)\n", i, (osc-esc), esc, osc);
	      if(fabs(osc-esc) > 0.0001) ndiff++;
	    }
	  FreeParsetree(tr);
	  free(dsq);
	  free(seq);
	}
      if(do_ptest) printf("Summary: Num diff: %d Bit diff: %.3f Avg per seq: %.3f\n", ndiff, diff, (diff/nseq));
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


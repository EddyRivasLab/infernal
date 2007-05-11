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
#include "cm_dispatch.h"	
#include <esl_vectorops.h>
#include <esl_histogram.h>
#include <esl_stats.h>

static char banner[] = "cmemit - generate sequences from a covariance model";

static char usage[]  = "\
Usage: cmemit [-options] <cm file>\n\
Available options are:\n\
   -l     : local; emit from a locally configured CM\n\
   -a     : write generated sequences as an alignment, not FASTA\n\
   -c     : generate a single \"consensus\" sequence\n\
   -h     : help; print brief help on version and usage\n\
   -n <n> : emit <n> sequences (default 10)\n\
   -o <f> : save sequences in file <f>\n\
   -q     : quiet - suppress verbose banner\n\
";

static char experts[] = "\
   --seed <n>     : set random number seed to <n>\n\
   --begin <n>    : truncate alignment, begin at match column <n>\n\
   --end   <n>    : truncate alignment, end   at match column <n>\n\
   --hmmbuild     : build a ML CM Plan 9 HMM from the samples\n\
   --hmmscore     : score samples with a CM Plan 9 HMM\n\
   --hmmlocal     : w/hmmscore, search with CP9 HMM in local mode\n\
   --cmbitsc      : w/hmmscore, min CM bit score to consider\n\
";

static struct opt_s OPTIONS[] = {
  { "-h",        TRUE,  sqdARG_NONE }, 
  { "-l",        TRUE, sqdARG_NONE },
  { "-a",        TRUE,  sqdARG_NONE },  
  { "-c",        TRUE,  sqdARG_NONE },  
  { "-n",        TRUE,  sqdARG_INT},  
  { "-o",        TRUE,  sqdARG_STRING},
  { "-q",        TRUE,  sqdARG_NONE},  
  { "--seed",    FALSE, sqdARG_INT},
  { "--begin",   FALSE, sqdARG_INT },
  { "--end",     FALSE, sqdARG_INT },
  { "--hmmbuild",FALSE, sqdARG_NONE },
  { "--hmmscore",FALSE, sqdARG_NONE },
  { "--hmmlocal",FALSE, sqdARG_NONE },
  { "--cmbitsc", FALSE, sqdARG_FLOAT }
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv) 
{
  char            *cmfile;      /* file to read CM from */	
  CMFILE          *cmfp;	/* open CM file */
  CM_t            *cm;          /* CM to generate from */

  FILE            *fp;          /* output file handle                      */
  char            *ofile;       /* output sequence file                    */
  int              L;		/* length of a sequence                    */
  int              i;		/* counter over sequences                  */

  int              nseq;	/* number of seqs to sample                */
  int              seed;	/* random number generator seed            */
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
  int   build_cp9;              /* TRUE to build a ML CP9 HMM from the 
				 * sampled parses of the CM                */
  int   do_score_cp9;           /* TRUE to score each seq against CP9 HMM  */
  int   do_local_cp9;           /* if score_cp9, put CP9 in local mode     */
  float cm_minsc = IMPOSSIBLE;  

  /*********************************************** 
   * Parse command line
   ***********************************************/

  nseq         = 10;
  seed         = time ((time_t *) NULL);
  be_quiet     = FALSE;
  do_local     = FALSE;
  do_alignment = FALSE;  
  do_consensus = FALSE;
  begin_set    = FALSE;
  end_set      = FALSE;
  build_cp9    = FALSE;
  do_score_cp9 = FALSE;
  do_local_cp9 = FALSE;
  ofile        = NULL;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-l")      == 0) do_local     = TRUE;
    else if (strcmp(optname, "-a")      == 0) do_alignment = TRUE;
    else if (strcmp(optname, "-c")      == 0) do_consensus = TRUE;
    else if (strcmp(optname, "-n")      == 0) nseq         = atoi(optarg); 
    else if (strcmp(optname, "-o")      == 0) ofile        = optarg;
    else if (strcmp(optname, "-q")      == 0) be_quiet     = TRUE;
    else if (strcmp(optname, "--seed")  == 0) seed         = atoi(optarg);
    else if (strcmp(optname, "--begin") == 0) { begin_set  = TRUE; spos = atoi(optarg); }
    else if (strcmp(optname, "--end")   == 0) { end_set    = TRUE; epos = atoi(optarg); }
    else if (strcmp(optname, "--hmmbuild") == 0) build_cp9 = TRUE;
    else if (strcmp(optname, "--hmmscore") == 0) do_score_cp9 = TRUE;
    else if (strcmp(optname, "--hmmlocal") == 0) do_local_cp9 = TRUE;
    else if (strcmp(optname, "--cmbitsc")  == 0) cm_minsc  = atof(optarg);
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

  sre_srandom(seed);

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
  if(do_local_cp9 && !do_score_cp9)
    Die("ERROR --hmmlocal only makes sense in combination with --hmmscore.\n");
  if(do_score_cp9 && (do_alignment || do_consensus))
    Die("ERROR --hmmscore does not work in combination with -a or -c.\n");
  
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

  if(do_local)        cm->config_opts |= CM_CONFIG_LOCAL;
  if(do_score_cp9)    cm->search_opts |= CM_SEARCH_HMMONLY;
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
	  EmitParsetree(cm, &(tr[i]), NULL, &(dsq[i]), &L);
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
	      EmitParsetree(cm, &(tr[i]), NULL, &(dsq[i]), &L);
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
      float             *hmm_sc;
      int               max_attempts = 500 * nseq;
      int               nattempts = 0;
      float             f;
      float             hmm_cutoff = 0.;
      double           *xv;
      int               n;
      double            mean, var;
      ESL_HISTOGRAM *h = esl_histogram_CreateFull(100, 1000, 0.2);

      hmm_sc = MallocOrDie(sizeof(float) * nseq);

      /* Put the CP9 in global alignment mode, unless --hmmlocal was specified */
      if(!do_local_cp9)
	{
	  CPlan9GlobalConfig(cm->cp9);
	  CP9Logoddsify(cm->cp9);
	}
      for (i = 0; i < nseq; i++)
	{
	  if(nattempts++ > max_attempts) 
	    Die("ERROR number of attempts exceeded 500 times number of seqs.\n");
	  EmitParsetree(cm, &tr, &seq, &dsq, &L);
	  cm_sc = ParsetreeScore(cm, tr, dsq, FALSE);
	  while(cm_sc < cm_minsc)
	    {
	      FreeParsetree(tr);
	      free(dsq);
	      free(seq);
	      EmitParsetree(cm, &tr, &seq, &dsq, &L);
	      cm_sc = ParsetreeScore(cm, tr, dsq, FALSE);
	      if(nattempts++ > max_attempts)
		Die("ERROR number of attempts exceeded 50 times number of seqs.\n");
	    }
	  esl_histogram_Add(h, cm_sc);
	  hmm_sc[i] = actually_search_target(cm, dsq, 1, L,
					     cm->cutoff, cm->cp9_cutoff,
					     NULL,  /* don't report hits to a results structure */
					     FALSE, /* we're not filtering with a CP9 HMM */
					     FALSE, /* we're not building a histogram for CM stats  */
					     FALSE, /* we're not building a histogram for CP9 stats */
					     NULL); /* filter fraction, irrelevant here */
	  FreeParsetree(tr);
	  free(dsq);
	  free(seq);
	  printf("i: %4d cm sc: %10.4f hmm sc: %10.4f\n", i, cm_sc, hmm_sc[i]);
	}
      /* Sort the HMM scores with quicksort */
      esl_vec_FSortIncreasing(hmm_sc, nseq);

      esl_histogram_Print(stdout, h);

      esl_histogram_GetData(h, &xv, &n);
      esl_stats_Mean(xv, n, &mean, &var);
      printf("N: %d\nMean:   %f\nVar:    %f\nSt dev: %f\n", n, mean, var, sqrt(var));

      esl_histogram_Destroy(h);
      
      f = 0.9;
      while(f < 1.001)
	{
	  hmm_cutoff = hmm_sc[(int) ((1. - f) * (float) nseq)];
	  printf("HMM glocal filter bit score threshold for finding %.2f CM hits > %.2f bits: %.4f\n", f, cm_minsc, hmm_cutoff);
	  f += 0.01;
	} 
      printf("\n\nnattempts: %d\n", nattempts);
      free(hmm_sc);
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


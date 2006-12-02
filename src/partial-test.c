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
Usage: partial-testt [-options] <cmfile>\n\
  where options are:\n\
  -h     : help; print brief help on version and usage\n\
  -n <n> : number of seqs to emit, truncate and align [default: 100]\n\
  -s <n> : set random number seed to <n> \n\
";

static char experts[] = "\
  --sub      : aln w/sub CM for columns b/t HMM predicted start/end points\n\
  --fsub     : aln w/sub CM for structure b/t HMM predicted start/end points\n\
  --cp9      : aln w/CM plan 9 HMM\n\
  --global   : run alignment in global mode [default: local]\n\
  --post <f> : minimum posterior prob from CP9 F/B to include in sub CM\n\
";

static struct opt_s OPTIONS[] = { 
  { "-h", TRUE, sqdARG_NONE }, 
  { "-n", TRUE, sqdARG_INT },
  { "-s", TRUE, sqdARG_INT },
  { "--sub",       FALSE, sqdARG_NONE},
  { "--fsub",      FALSE, sqdARG_NONE},
  { "--cp9",       FALSE, sqdARG_NONE },
  { "--global",    FALSE, sqdARG_NONE },
  { "--post",      FALSE, sqdARG_FLOAT }
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
  int print_flag;               /* TRUE to print debug statements */
  int L;

  /*********************************************** 
   * Parse command line
   ***********************************************/

  seed           = (int) time ((time_t *) NULL);
  pthresh        = 0.1;
  do_sub         = FALSE;
  do_fullsub     = FALSE;
  nseq           =  100;
  print_flag     = FALSE;
  do_local       = FALSE;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
		&optind, &optname, &optarg))  {
    if      (strcmp(optname, "-n") == 0) nseq         = atoi(optarg);
    else if (strcmp(optname, "-s") == 0) seed           = atoi(optarg);
    else if (strcmp(optname, "--sub")       == 0) do_sub = TRUE;
    else if (strcmp(optname, "--fsub")      == 0) do_fullsub = TRUE;
    else if (strcmp(optname, "--global")     == 0) do_local = FALSE;
    else if (strcmp(optname, "--debug")     == 0) print_flag = TRUE;
    else if (strcmp(optname, "--post")      == 0) pthresh  = atof(optarg); 
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
  
  Parsetree_t *tr1;             /* Parsetree of emitted sequence */
  Parsetree_t *tr2;             /* Parsetree of aligned sequence */
  char    *dsq;                /* digitized sequence                     */
  MSA               *msa;       /* alignment */
  int *useme;
  int apos;
  int cc;
  int *ct;		/* 0..alen-1 base pair partners array         */
  int attempts = 0;
  int eq = 0;
  while(i < nseq)
    {
      attempts++;
      EmitParsetree(cm, &tr1, NULL, &dsq, &L);
      CYKDivideAndConquer(cm, dsq, L, 0, 1, L, &tr2, NULL, NULL);
      if(ParsetreeCompare(tr1, tr2))
	{
	  i++;
	  printf("attempts: %3d passed (i=%3d)\n", attempts, i);
	}
      if (tr1 != NULL) FreeParsetree(tr1);  
      if (tr2 != NULL) FreeParsetree(tr2); 
      free(dsq);
    }	  
}

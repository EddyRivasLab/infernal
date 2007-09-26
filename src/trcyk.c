#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_alphabet.c"
#include "esl_msa.h"
#include "esl_sqio.h"

#include "structs.h"
#include "funcs.h"

int
main(int argc, char **argv)
{
   char          *cmfile;
   ESL_ALPHABET  *abc;
   char          *seqfile;
   ESL_SQFILE    *sqfp;
   int            format;
   CMFILE        *cmfp;
   CM_t          *cm;
   ESL_SQ        *seq;
   float          sc;
   Parsetree_t   *tr;
   Fancyali_t    *fali;
   CMConsensus_t *cons;

   int do_local;

   /* int status;    */
   /* char *optname; */
   /* char *optarg; */
   int   optind;

   cmfile = seqfile = NULL;
   abc = NULL;
   sqfp = NULL;
   cmfp = NULL;
   cm = NULL;
   seq = NULL;
   tr = NULL;
   fali = NULL;
   cons = NULL;
   format = eslSQFILE_UNKNOWN;
   do_local = TRUE;

   /* Should process options, but for now assume none and set optind */
   optind = 1;

   if ( argc - optind != 2 ) cm_Die("Incorrect number of arguments\n");
   cmfile = argv[optind++];
   seqfile = argv[optind++];

   if ( (cmfp = CMFileOpen(cmfile, NULL)) == NULL )
      cm_Die("Failed to open covariance model save file\n");
   if (! CMFileRead(cmfp, &abc, &cm))
      cm_Die("Failed to read a CM from cm file\n");
   if (cm == NULL)
      cm_Die("CM file empty?\n");
   CMFileClose(cmfp);

   if ( esl_sqfile_Open(seqfile, format, NULL, &sqfp) != eslOK )
      cm_Die("Failed to open sequence database file\n");

   if (do_local) cm->config_opts |= CM_CONFIG_LOCAL;

   ConfigCM(cm, NULL, NULL);
   CreateCMConsensus(cm, cm->abc, 3.0, 1.0, &cons);

   seq = esl_sq_Create();
   while ( esl_sqio_Read(sqfp, seq) == eslOK )
   {
      if (seq->n == 0) continue;

      printf("sequence: %s\n", seq->name);

      int i0 = 1;
      int j0 = seq->n;
      
      if (seq->dsq == NULL) 
         esl_sq_Digitize(abc, seq);
      /* Do alignment */
      sc = TrCYK_Inside(cm, seq->dsq, seq->n, 0, i0, j0, &tr);
      /* Print alignment */
      printf("score:    %.2f\n",sc);
      fali = CreateFancyAli(tr, cm, cons, seq->dsq, cm->abc);
      PrintFancyAli(stdout, fali, 
		    0,      /* offset in seq index */
		    FALSE); /* not on reverse complement strand */
      FreeFancyAli(fali);

      printf("sequence: %s (reversed)\n", seq->name);
      revcomp(abc, seq, seq);
      sc = TrCYK_Inside(cm,seq->dsq, seq->n, 0, i0, j0, &tr);
      printf("score:    %.2f\n",sc);
      fali = CreateFancyAli(tr, cm, cons,seq->dsq, cm->abc);
      PrintFancyAli(stdout, fali, 0, FALSE);

      FreeFancyAli(fali);
   }
   esl_sq_Destroy(seq);

   FreeCMConsensus(cons);
   FreeCM(cm);
   esl_sqfile_Close(sqfp);

   return EXIT_SUCCESS;
}

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
   float          sc, rev_sc;
   Parsetree_t   *tr;
   Fancyali_t    *fali;
   Fancyali_t    *rev_fali;
   CMConsensus_t *cons;

   int do_local;

   /* int status;    */
   /* char *optname; */
   /* char *optarg; */
   int   optind;

   int status;
   char errbuf[cmERRBUFSIZE];

   cmfile = seqfile = NULL;
   abc = NULL;
   sqfp = NULL;
   cmfp = NULL;
   cm = NULL;
   seq = NULL;
   tr = NULL;
   fali = NULL;
   rev_fali = NULL;
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
   if ((status = CMFileRead(cmfp, errbuf, &abc, &cm)) != eslOK)
      cm_Die("Failed to read a CM from cm file\n");
   if (cm == NULL)
      cm_Die("CM file empty?\n");
   CMFileClose(cmfp);

   if ( esl_sqfile_Open(seqfile, format, NULL, &sqfp) != eslOK )
      cm_Die("Failed to open sequence database file\n");

   if (do_local) cm->config_opts |= CM_CONFIG_LOCAL;

   ConfigCM(cm, TRUE); /* TRUE says: calculate W */
   CreateCMConsensus(cm, cm->abc, 3.0, 1.0, &cons);
   SetMarginalScores(cm);

   seq = esl_sq_Create();
   while ( esl_sqio_Read(sqfp, seq) == eslOK )
   {
      if (seq->n == 0) continue;

      int i0 = 1;
      int j0 = seq->n;
      
      if (seq->dsq == NULL) 
         esl_sq_Digitize(abc, seq);
      sc = TrCYK_DnC(cm, seq->dsq, seq->n, 0, i0, j0, &tr);
      fali = CreateFancyAli(cm->abc, tr, cm, cons, seq->dsq, NULL, NULL);
      FreeParsetree(tr);

      revcomp(abc, seq, seq);
      rev_sc = TrCYK_DnC(cm,seq->dsq, seq->n, 0, i0, j0, &tr);
      rev_fali = CreateFancyAli(cm->abc, tr, cm, cons,seq->dsq, NULL, NULL);
      FreeParsetree(tr);

      if (sc > rev_sc)
      {
         printf("sequence: %s\n", seq->name);
         printf("score:    %.2f\n",sc);
         PrintFancyAli(stdout, fali, 0, FALSE, FALSE);
      }
      else
      {
         printf("sequence: %s (reversed)\n", seq->name);
         printf("score:    %.2f\n",rev_sc);
         PrintFancyAli(stdout, fali, seq->n, TRUE, FALSE);
      }

      FreeFancyAli(fali);
      FreeFancyAli(rev_fali);

   esl_sq_Destroy(seq);
   seq = esl_sq_Create();

   }
   esl_sq_Destroy(seq);

   FreeCMConsensus(cons);
   FreeCM(cm);
   esl_sqfile_Close(sqfp);

   return EXIT_SUCCESS;
}

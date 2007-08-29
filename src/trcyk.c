#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "squid.h"

#include "structs.h"
#include "funcs.h"

int
main(int argc, char **argv)
{
   char          *cmfile;
   char          *seqfile;
   SQFILE        *sqfp;
   int            format;
   CMFILE        *cmfp;
   CM_t          *cm;
   char          *seq;
   SQINFO         sqinfo;
   char          *dsq;
   float          sc;
   Parsetree_t   *tr;
   Fancyali_t    *fali;
   CMConsensus_t *cons;

   int do_local;

   /* int status;    */
   /* char *optname; */
   /* char *optarg; */
   int   optind;

   format = SQFILE_UNKNOWN;
   do_local = TRUE;

   /* Should process options, but for now assume none and set optind */
   optind = 1;

   if ( argc - optind != 2 ) Die("Incorrect number of arguments\n");
   cmfile = argv[optind++];
   seqfile = argv[optind++];

   if ( (cmfp = CMFileOpen(cmfile, NULL)) == NULL )
      Die("Failed to open covariance model save file\n");
   if (! CMFileRead(cmfp, &cm))
      Die("Failed to read a CM from cm file\n");
   if (cm == NULL)
      Die("CM file empty?\n");
   CMFileClose(cmfp);

   if ( (sqfp = SeqfileOpen(seqfile, format, NULL)) == NULL)
      Die("Failed to open sequence database file\n");

   if (do_local) cm->config_opts |= CM_CONFIG_LOCAL;

   ConfigCM(cm, NULL, NULL);
   CreateCMConsensus(cm, 3.0, 1.0, &con);

   while ( ReadSeq(sqfp, sqfp->format, &seq, &sqinfo) )
   {
      if (sqinfo.len == 0) continue;
      dsq = DigitizeSequence(seq, sqinfo.len);

      printf("sequence: %s\n", sqinfo.name);

int i0 = 1;
int j0 = sqinfo.len;
      /* Do alignment */
      sc = TrCYKInside(cm, dsq, sqinfo.len, 0, i0, j0, &tr, NULL, NULL);
      /* Print alignment */
/* ParsetreeDump(stdout,tr,cm,dsq); printf("\n\n"); */
      printf("score:    %.2f\n",sc);
      fali = CreateFancyAli(tr, cm, cons, dsq);
      PrintFancyAli(stdout, fali, 
		    0,      /* offset in seq index */
		    FALSE); /* not on reverse complement strand */
      FreeFancyAli(fali);
      free(dsq);
      FreeSequence(seq, &sqinfo);
   }

   FreeCMConsensus(cons);
   FreeCM(cm);
   SeqfileClose(sqfp);

   return EXIT_SUCCESS;
}

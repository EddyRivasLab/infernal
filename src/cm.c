/* cm.c
 * SRE, Sat Jul 29 09:01:20 2000 [St. Louis]
 * CVS $Id$
 * 
 * Routines for dealing with the CM data structure.
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************  
 */

#include <stdio.h>
#include <stdlib.h>

#include "squid.h"

#include "structs.h"
#include "funcs.h"


/* Function: CreateCM()
 * Date:     SRE, Sat Jul 29 09:02:16 2000 [St. Louis]
 *
 * Purpose:  Create a covariance model, given the number of states 
 *           that should be in it.
 *
 * Args:     nstates = number of states in the model
 *
 * Returns:  ptr to allocated cm. 
 *           Caller is responsible for free'ing the cm.
 */
CM_t *
CreateCM(int nstates)
{
  CM_t *cm;

  cm = MallocOrDie(sizeof(CM_t));

				/* general information: added later */
  cm->name = NULL;
  cm->acc  = NULL;
  cm->desc = NULL;
  cm->M    = nstates;
				/* structural information */
  cm->sttype = MallocOrDie(nstates * sizeof(char));
  cm->ndidx  = MallocOrDie(nstates * sizeof(int));
  cm->ndtype = MallocOrDie(nstates * sizeof(char));
  cm->stid   = MallocOrDie(nstates * sizeof(char));
  cm->cfirst = MallocOrDie(nstates * sizeof(int));
  cm->cnum   = MallocOrDie(nstates * sizeof(int));
  cm->plast  = MallocOrDie(nstates * sizeof(int));
  cm->pnum   = MallocOrDie(nstates * sizeof(int));
				/* parameter information */
  cm->t      = FMX2Alloc(nstates, MAXCONNECT);
  cm->e      = FMX2Alloc(nstates, Alphabet_size*Alphabet_size);
  cm->tsc    = FMX2Alloc(nstates, MAXCONNECT);
  cm->esc    = FMX2Alloc(nstates, Alphabet_size*Alphabet_size);

  return cm;
}

/* Function: FreeCM()
 * Date:     SRE, Sat Jul 29 11:22:32 2000 [St. Louis]
 *
 * Purpose:  Free a CM data structure.
 *
 * Args:     cm - the model to free. (duh).
 *
 * Returns:  (void)
 */
void
FreeCM(CM_t *cm)
{
  if (cm->name != NULL) free(cm->name);
  if (cm->acc  != NULL) free(cm->acc);
  if (cm->desc != NULL) free(cm->desc);

  free(cm->sttype);
  free(cm->ndidx);
  free(cm->ndtype);
  free(cm->stid);
  free(cm->cfirst);
  free(cm->cnum);
  free(cm->plast);
  free(cm->pnum);
  FMX2Free(cm->t);
  FMX2Free(cm->e);
  FMX2Free(cm->tsc);
  FMX2Free(cm->esc);
  free(cm);
}



/* Function: PrintCM()
 * Date:     SRE, Sat Jul 29 10:55:16 2000 [St. Louis]
 *
 * Purpose:  Debugging: show a tabular representation of a CM structure.
 *
 * Args:     fp - output stream (e.g. stdout)
 *           cm - the CM to show
 *
 * Returns:  (void)
 */
void
PrintCM(FILE *fp, CM_t *cm)
{
  int x;

  fprintf(fp, "%5s %6s %5s %6s %7s %6s %5s %5s %5s\n",
	  " idx ","sttype", "ndidx", "ndtype", "  stid ", "cfirst", " cnum", "plast", " pnum");
  fprintf(fp, "%5s %6s %5s %6s %7s %5s %5s %5s %5s\n",
	  "-----", "------", "-----", "------","-------","------","-----", "-----", "-----");
  
  for (x = 0; x < cm->M; x++)
    {
      fprintf(fp, "%5d %-6s %5d %6s %-7s %6d %5d %5d %5d\n",
	      x, Statetype(cm->sttype[x]), cm->ndidx[x], 
	      Nodetype(cm->ndtype[x]), UniqueStatetype(cm->stid[x]),
	      cm->cfirst[x], cm->cnum[x],
	      cm->plast[x], cm->pnum[x]);
    }
}

/* Function: SummarizeCM()
 * Date:     SRE, Sat Jul 29 12:19:31 2000 [St. Louis]
 *
 * Purpose:  Print some summary information about a new CM;
 *           called by cmbuild after each new model construction.
 *
 * Args:     fp - output stream (e.g. stdout)
 *           cm - cm to summarize
 *
 * Returns:  (void)
 */
void
SummarizeCM(FILE *fp, CM_t *cm)
{
  int x;
  int count[UNIQUESTATES];

  for (x = 0; x < UNIQUESTATES; x++) count[x] = 0;

  for (x = 0; x < cm->M; x++)
    count[cm->stid[x]]++;
  
  fprintf(fp, "Summary report for a covariance model:\n");
  fprintf(fp, "------------------------------------\n");
  fprintf(fp, "Total states:       %d\n", cm->M);
  fprintf(fp, "Total nodes:        %d\n", cm->nodes);
  fprintf(fp, "Bifurcations:       %d\n", count[BIF_B]);
  fprintf(fp, "Consensus pairs:    %d\n", count[MATP_MP]);
  fprintf(fp, "Consensus singlets: %d\n", count[MATL_ML]+count[MATR_MR]);
}

/* Functions: Statetype(), Nodetype(), UniqueStatetype()
 * Date:      SRE, Sat Jul 29 11:07:47 2000 [St. Louis]
 *
 * Purpose:   Translate internal flags into human-readable strings, 
 *            for clearer debugging output.
 * 
 * Args:      type - a state type, node type, or unique statetype
 *
 * Returns:   an appropriate string
 */
char *
Statetype(int type) 
{
  switch (type) {
  case D_st:  return "D";
  case MP_st: return "MP";
  case ML_st: return "ML";
  case MR_st: return "MR";
  case IL_st: return "IL";
  case IR_st: return "IR";
  case S_st:  return "S";
  case E_st:  return "E";
  case B_st:  return "B";
  default: Die("bogus state type %d\n", type);
  }
  return "";
}
char *
Nodetype(int type) 
{
  switch (type) {
  case DUMMY_nd: return "-";
  case BIF_nd:   return "BIF";
  case MATP_nd:  return "MATP";
  case MATL_nd:  return "MATL";
  case MATR_nd:  return "MATR";
  case BEGL_nd:  return "BEGL";
  case BEGR_nd:  return "BEGR";
  case ROOT_nd:  return "ROOT";
  case END_nd:   return "END";
  default: Die("bogus node type %d\n", type);
  }
  return "";
}
char *
UniqueStatetype(int type)
{
  switch (type) {
  case DUMMY:   return "DUMMY";   
  case ROOT_S:  return "ROOT_S";
  case ROOT_IL: return "ROOT_IL";
  case ROOT_IR: return "ROOT_IR";
  case BEGL_S : return "BEGL_S";
  case BEGR_S : return "BEGR_S";
  case BEGR_IL: return "BEGR_IL";
  case MATP_D : return "MATP_D";
  case MATP_MP: return "MATP_MP";
  case MATP_ML: return "MATP_ML";
  case MATP_MR: return "MATP_MR";
  case MATP_IL: return "MATP_IL";
  case MATP_IR: return "MATP_IR";
  case MATL_D : return "MATL_D";
  case MATL_ML: return "MATL_ML";
  case MATL_IL: return "MATL_IL";
  case MATR_D : return "MATR_D";
  case MATR_MR: return "MATR_MR";
  case MATR_IR: return "MATR_IR";
  case END_E  : return "END_E";
  case BIF_B  : return "BIF_B";
  default: Die("bogus unique state type %d\n", type);
  }
  return "";
}


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
#include <math.h>

#include "squid.h"

#include "structs.h"
#include "funcs.h"


/* Function: CreateCM()
 * Date:     SRE, Sat Jul 29 09:02:16 2000 [St. Louis]
 *
 * Purpose:  Create a covariance model, given the number of states 
 *           that should be in it.
 *
 * Args:     nnodes  =  number of nodes in the model
 *           nstates = number of states in the model
 *
 * Returns:  ptr to allocated cm. 
 *           Caller is responsible for free'ing the cm.
 */
CM_t *
CreateCM(int nnodes, int nstates)
{
  CM_t *cm;

  cm = MallocOrDie(sizeof(CM_t));

				/* general information: added later */
  cm->name = NULL;
  cm->acc  = NULL;
  cm->desc = NULL;
  cm->M    = nstates;
				/* null model information */
  cm->null   = MallocOrDie(Alphabet_size * sizeof(float));
				/* structural information */
  cm->sttype = MallocOrDie(nstates * sizeof(char));
  cm->ndidx  = MallocOrDie(nstates * sizeof(int));
  cm->stid   = MallocOrDie(nstates * sizeof(char));
  cm->cfirst = MallocOrDie(nstates * sizeof(int));
  cm->cnum   = MallocOrDie(nstates * sizeof(int));
  cm->plast  = MallocOrDie(nstates * sizeof(int));
  cm->pnum   = MallocOrDie(nstates * sizeof(int));
				/* node->state map information */
  cm->nodes  = nnodes;
  cm->nodemap= MallocOrDie(nnodes  * sizeof(int));
  cm->ndtype = MallocOrDie(nnodes  * sizeof(char));
				/* parameter information */
  cm->t      = FMX2Alloc(nstates, MAXCONNECT);
  cm->e      = FMX2Alloc(nstates, Alphabet_size*Alphabet_size);
  cm->tsc    = FMX2Alloc(nstates, MAXCONNECT);
  cm->esc    = FMX2Alloc(nstates, Alphabet_size*Alphabet_size);

  return cm;
}

/* Function: CMZero()
 * Date:     SRE, Mon Jul 31 19:14:31 2000 [St. Louis]
 *
 * Purpose:  Initialize the probability parameters of a CM to zero.
 *
 * Returns:  (void)
 */
void 
CMZero(CM_t *cm)
{
  int v;			/* counter over states                 */
  int x;			/* counter over symbols or transitions */

  for (v = 0; v < cm->M; v++) {
    for (x = 0; x < Alphabet_size * Alphabet_size; x++) cm->e[v][x] = 0.0;
    for (x = 0; x < MAXCONNECT; x++)                    cm->t[v][x] = 0.0;
  }
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

  free(cm->null);
  free(cm->sttype);
  free(cm->ndidx);
  free(cm->stid);
  free(cm->cfirst);
  free(cm->cnum);
  free(cm->plast);
  free(cm->pnum);
  free(cm->nodemap);
  free(cm->ndtype);
  FMX2Free(cm->t);
  FMX2Free(cm->e);
  FMX2Free(cm->tsc);
  FMX2Free(cm->esc);
  free(cm);
}


/* Function: CMSetDefaultNullModel()
 * Date:     SRE, Tue Aug  1 15:31:52 2000 [St. Louis]
 *
 * Purpose:  Initialize the null model to equiprobable (e.g. 0.25)
 */
void
CMSetDefaultNullModel(CM_t *cm)
{
  int x;
  for (x = 0; x < Alphabet_size; x++)
    cm->null[x] = 1./(float)Alphabet_size;
}


/* Function: CMSimpleProbify()
 * Date:     SRE, Tue Aug  1 11:07:17 2000 [St. Louis]
 *
 * Purpose:  Convert a counts-based CM to probability form, using
 *           a plus-one Laplace prior.
 */
void
CMSimpleProbify(CM_t *cm)
{
  int v,x;

  for (v = 0; v < cm->M; v++) 
    {
      /* Transitions. B, E have no transition probabilities.
       */
      if (cm->sttype[v] != B_st && cm->sttype[v] != E_st) 
	{
	  for (x = 0; x < cm->cnum[v]; x++) cm->t[v][x] += 1.0; /* Laplace prior */
	  FNorm(cm->t[v], cm->cnum[v]);	                        /* normalize to a probability */
	}

      /* Emissions.
       */
      if (cm->sttype[v] == MP_st) 
	{
	  for (x = 0; x < Alphabet_size*Alphabet_size; x++) cm->e[v][x] += 1.0;
	  FNorm(cm->e[v], Alphabet_size*Alphabet_size);
	}
      else if (cm->sttype[v] == ML_st || cm->sttype[v] == MR_st || 
	       cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) 
	{
	  for (x = 0; x < Alphabet_size; x++) cm->e[v][x] += 1.0;
	  FNorm(cm->e[v], Alphabet_size);
	}
    }
}

/* Function: CMLogoddsify()
 * Date:     SRE, Tue Aug  1 15:18:26 2000 [St. Louis]
 *
 * Purpose:  Convert the probabilities in a CM to log-odds form.
 */
void
CMLogoddsify(CM_t *cm)
{
  int v, x, y;

  for (v = 0; v < cm->M; v++)
    {
      if (cm->sttype[v] != B_st && cm->sttype[v] != E_st)
	for (x = 0; x < cm->cnum[v]; x++)
	  cm->tsc[v][x] = log(cm->t[v][x]);
      
      if (cm->sttype[v] == MP_st)
	for (x = 0; x < Alphabet_size; x++)
	  for (y = 0; y < Alphabet_size; y++)
	    cm->esc[v][x*Alphabet_size+y] = log(cm->e[v][x*Alphabet_size+y] / (cm->null[x]*cm->null[y]));

      if (cm->sttype[v] == ML_st || cm->sttype[v] == MR_st ||
	  cm->sttype[v] == IL_st || cm->sttype[v] == IR_st)
	for (x = 0; x < Alphabet_size; x++)
	  cm->esc[v][x] = log(cm->e[v][x] / cm->null[x]);
    }
}

/* Function: CMCountStatetype()
 * Date:     SRE, Wed Aug  2 09:15:00 2000 [St. Louis]
 *
 * Purpose:  A convenience for counting the # of occurrences
 *           of a particular state type in a CM. Useful for
 *           "how many bifurcations does this model have", etc.
 *
 * Args:     cm   - the model
 *           type - a state type (e.g. E_st or MP_st)    
 *
 * Returns:  how many states of that type are in the model
 */
int
CMCountStatetype(CM_t *cm, char type)
{
  int v;
  int count = 0;
  for (v = 0; v < cm->M; v++)
    if (cm->sttype[v] == type) count++;
  return count;
}


/* Function: CalculateStateIndex()
 * Date:     SRE, Mon Jul 31 15:37:55 2000 [St. Louis]
 *
 * Purpose:  Given a node index and a unique state type, use the CM's
 *           nodemap to calculate and return a state index in the CM.
 *
 *           Doesn't check that the node type matches what's implied
 *           by the utype! (e.g., if you pass utype==MATP_MP, the node
 *           had better be a MATP.)
 *
 * Args:     cm     - the covariance model
 *           node   - node index, 0..cm->nodes-1
 *           utype  - unique statetype, e.g. MATP_MP
 *
 * Returns:  a state index, 0..cm->M-1
 *
 * Used in:  modelmaker.c:transmogrify() 
 */
int
CalculateStateIndex(CM_t *cm, int node, char utype)
{
  int base;

  base = cm->nodemap[node];
  switch (utype) {
  case ROOT_S:  return base;
  case ROOT_IL: return base+1;
  case ROOT_IR: return base+2;
  case BEGL_S:  return base;
  case BEGR_S:  return base;
  case BEGR_IL: return base+1;
  case MATP_MP: return base;
  case MATP_ML: return base+1;
  case MATP_MR: return base+2;
  case MATP_D:  return base+3;  
  case MATP_IL: return base+4;
  case MATP_IR: return base+5; 
  case MATL_ML: return base;
  case MATL_D:  return base+1;
  case MATL_IL: return base+2;
  case MATR_MR: return base;
  case MATR_D:  return base+1;
  case MATR_IR: return base+2;
  case END_E:   return base;
  case BIF_B:   return base;
  default: Die("bogus utype %d in CalculateStateIndex()", utype);
  }
  return base;			/* not used */
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
	      Nodetype(cm->ndtype[cm->ndidx[x]]), UniqueStatetype(cm->stid[x]),
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
    count[(int) cm->stid[x]]++;
  
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
  case MATP_MP: return "MATP_MP";
  case MATP_ML: return "MATP_ML";
  case MATP_MR: return "MATP_MR";
  case MATP_D : return "MATP_D";
  case MATP_IL: return "MATP_IL";
  case MATP_IR: return "MATP_IR";
  case MATL_ML: return "MATL_ML";
  case MATL_D : return "MATL_D";
  case MATL_IL: return "MATL_IL";
  case MATR_MR: return "MATR_MR";
  case MATR_D : return "MATR_D";
  case MATR_IR: return "MATR_IR";
  case END_E  : return "END_E";
  case BIF_B  : return "BIF_B";
  default: Die("bogus unique state type %d\n", type);
  }
  return "";
}


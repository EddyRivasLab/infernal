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
#include <string.h>
#include <math.h>

#include "squid.h"
#include "vectorops.h"

#include "nstack.h"
#include "structs.h"
#include "funcs.h"


/* Function: CreateCM(); CreateCMShell(); CreateCMBody()
 * Date:     SRE, Sat Jul 29 09:02:16 2000 [St. Louis]
 *
 * Purpose:  Create a covariance model, given the number of states 
 *           and nodes that should be in it.
 *
 *           Allocation is usually one step: CreateCM(N, M).
 *           
 *           In cmio.c, allocation is two-step: CreateCMShell()
 *           then CreateCMBody(cm, N, M). This way we can create
 *           a model, start storing header info in it (including M
 *           and N themselves), and only after we read M and N
 *           do we allocate the bulk of the model.
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

  cm = CreateCMShell();
  if (cm != NULL) CreateCMBody(cm, nnodes, nstates);
  return cm;
}
CM_t *
CreateCMShell(void)
{
  CM_t *cm;

  cm = MallocOrDie(sizeof(CM_t));

				/* general information: added later */
  cm->name   = NULL;
  cm->acc    = NULL;
  cm->desc   = NULL;
  cm->annote = NULL;

				/* null model information */
  cm->null   = MallocOrDie(Alphabet_size * sizeof(float));
				/* structural information */
  cm->M      = 0;
  cm->sttype = NULL;
  cm->ndidx  = NULL;
  cm->stid   = NULL;
  cm->cfirst = NULL;
  cm->cnum   = NULL;
  cm->plast  = NULL;
  cm->pnum   = NULL;
				/* node->state map information */
  cm->nodes  = 0;
  cm->nodemap= NULL;
  cm->ndtype = NULL;
				/* parameter information */
  cm->t      = NULL;
  cm->e      = NULL;
  cm->begin  = NULL;
  cm->end    = NULL;
  cm->tsc    = NULL;
  cm->esc    = NULL;
  cm->beginsc= NULL;
  cm->endsc  = NULL;

  cm->flags  = 0;
  return cm;
}
void
CreateCMBody(CM_t *cm, int nnodes, int nstates)
{
				/* structural information */
  cm->M      = nstates;
  cm->sttype = MallocOrDie((nstates+1) * sizeof(char));
  cm->ndidx  = MallocOrDie(nstates * sizeof(int));
  cm->stid   = MallocOrDie((nstates+1) * sizeof(char));
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
  cm->begin  = MallocOrDie(nstates * sizeof(float));
  cm->end    = MallocOrDie(nstates * sizeof(float));
  cm->tsc    = FMX2Alloc(nstates, MAXCONNECT);
  cm->esc    = FMX2Alloc(nstates, Alphabet_size*Alphabet_size);
  cm->beginsc= MallocOrDie(nstates * sizeof(float));
  cm->endsc  = MallocOrDie(nstates * sizeof(float));

  /* the EL state at M is special: we only need state
   * type info recorded, so functions looking at parsetrees  
   * can interpret what an "M" index means.
   */
  cm->sttype[cm->M] = EL_st;
  cm->stid[cm->M]   = END_EL;

  cm->flags  = 0;
}


/* Function: CMZero()
 * Date:     SRE, Mon Jul 31 19:14:31 2000 [St. Louis]
 *
 * Purpose:  Initialize the probability parameters and scores of a CM to zero.
 *
 * Returns:  (void)
 */
void 
CMZero(CM_t *cm)
{
  int v;			/* counter over states                 */
  int x;			/* counter over symbols or transitions */

  for (v = 0; v < cm->M; v++) {
    for (x = 0; x < Alphabet_size * Alphabet_size; x++) cm->e[v][x]   = 0.0;
    for (x = 0; x < MAXCONNECT; x++)                    cm->t[v][x]   = 0.0;
    for (x = 0; x < Alphabet_size * Alphabet_size; x++) cm->esc[v][x] = 0.0;
    for (x = 0; x < MAXCONNECT; x++)                    cm->tsc[v][x] = 0.0;
    cm->begin[v] = cm->end[v] = 0.;
    cm->beginsc[v] = cm->endsc[v] = 0.;
  }
}

/* Function:  CMRenormalize()
 * Incept:    SRE, Wed Aug 14 14:16:55 2002 [St. Louis]
 *
 * Purpose:   Renormalize all the probability distributions in a CM.
 *            Used by cmio.c's flatfile parser, for example.
 *
 * Xref:      STL6 p.108
 */
void
CMRenormalize(CM_t *cm)
{
  int v;

  FNorm(cm->null, Alphabet_size);
  for (v = 0; v < cm->M; v++)
    {
      if (cm->cnum[v] > 0 && cm->sttype[v] != B_st)
	FNorm(cm->t[v], cm->cnum[v]);
      
      if (cm->sttype[v] == ML_st || cm->sttype[v] == MR_st || cm->sttype[v] == IL_st || cm->sttype[v] == IR_st)
	FNorm(cm->e[v], Alphabet_size);
      if (cm->sttype[v] == MP_st)
	FNorm(cm->e[v], Alphabet_size * Alphabet_size);
    }
  if (cm->flags & CM_LOCAL_BEGIN) FNorm(cm->begin, cm->M);
  if (cm->flags & CM_LOCAL_END)   Die("Renormalization of models in local end mode not supported yet");
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
  if (cm->name   != NULL) free(cm->name);
  if (cm->acc    != NULL) free(cm->acc);
  if (cm->desc   != NULL) free(cm->desc);
  if (cm->annote != NULL) free(cm->annote);

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
  free(cm->begin);
  free(cm->end);
  FMX2Free(cm->tsc);
  FMX2Free(cm->esc);
  free(cm->beginsc);
  free(cm->endsc);
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
	  cm->tsc[v][x] = sreLOG2(cm->t[v][x]);
      
      if (cm->sttype[v] == MP_st)
	for (x = 0; x < Alphabet_size; x++)
	  for (y = 0; y < Alphabet_size; y++)
	    cm->esc[v][x*Alphabet_size+y] = sreLOG2(cm->e[v][x*Alphabet_size+y] / (cm->null[x]*cm->null[y]));

      if (cm->sttype[v] == ML_st || cm->sttype[v] == MR_st ||
	  cm->sttype[v] == IL_st || cm->sttype[v] == IR_st)
	for (x = 0; x < Alphabet_size; x++)
	  cm->esc[v][x] = sreLOG2(cm->e[v][x] / cm->null[x]);

      /* These work even if begin/end distributions are inactive 0's,
       * sreLOG2 will set beginsc, endsc to -infinity.
       */
      cm->beginsc[v] = sreLOG2(cm->begin[v]);
      cm->endsc[v]   = sreLOG2(cm->end[v]);
    }
}

/* Function:  CMHackInsertScores()
 * Incept:    SRE, Wed Jul 24 09:48:22 2002 [St. Louis]
 *
 * Purpose:   Temporary (I hope): make all insert scores 0.
 *            If you let inserts train on the data, you can get
 *            positive insert emission scores. Local alignments,
 *            in particular, can then consist of just a couple of
 *            consensus states and a long string of insert 
 *            states, hitting base-composition-biased sequence
 *            with very high score. This is a Bad Thing.
 *            
 *            The long term solution for this problem will
 *            go in with mixture Dirichlet priors, but for now
 *            (with only Laplace coded), this'll appease the
 *            pitchfork and torches mob at Cambridge.
 *
 * Args:      cm - the model 
 *
 * Returns:   (void)
 *
 * Xref:      STL6 p.93.
 */
void
CMHackInsertScores(CM_t *cm)
{
  int v, x;
  for (v = 0; v < cm->M; v++)
    {
      if (cm->sttype[v] == IL_st || cm->sttype[v] == IR_st)
	for (x = 0; x < Alphabet_size; x++)
	  cm->esc[v][x] = 0.;
    }
}


/* Function: CMCountStatetype(), CMSubtreeCountStatetype(), CMSegmentCountStatetype
 * Date:     SRE, Wed Aug  2 09:15:00 2000 [St. Louis]
 *
 * Purpose:  Conveniences for counting the # of occurrences
 *           of a particular state type in a CM. Useful for
 *           "how many bifurcations does this model have", etc.
 *          
 *           CMSubtreeCountStatetype() only counts underneath     
 *           a particular subtree rooted at state v
 *
 * Args:     cm   - the model
 *           r    - the root of the subtree to start from (inclusive)
 *           z    - end of the subtree to stop at (inclusive) 
 *           type - a state type (e.g. E_st or MP_st)    
 *
 * Returns:  how many states of that type are in the model
 */
int
CMSegmentCountStatetype(CM_t *cm, int r, int z, char type)
{
  int count = 0;
  int v;
  for (v = r; v <= z; v++) 
    if (cm->sttype[v] == type) count++;
  return count;
}
int
CMSubtreeCountStatetype(CM_t *cm, int v, char type)
{
  int unsatisfied_starts = 1;
  int count = 0;

  while (unsatisfied_starts) {
    if (cm->sttype[v] == B_st) unsatisfied_starts++;
    if (cm->sttype[v] == E_st) unsatisfied_starts--; 
    if (cm->sttype[v] == type) count++;
    v++;
  }
  return count;
}
int
CMCountStatetype(CM_t *cm, char type)
{
  return CMSubtreeCountStatetype(cm, 0, type);
}
int 
CMSubtreeFindEnd(CM_t *cm, int r)
{
  int unsatisfied_starts = 1;

  while (unsatisfied_starts) {
    if (cm->sttype[r] == B_st) unsatisfied_starts++;
    if (cm->sttype[r] == E_st) unsatisfied_starts--; 
    r++;
  }
  return (r-1);
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

/* Function:  TotalStatesInNode(), SplitStatesInNode(), InsertStatesInNode()
 * Incept:    SRE, Thu Aug  8 09:57:59 2002 [St. Louis]
 *
 * Purpose:   Returns the number of states in a node type.
 *
 * Args:      ndtype  - type of node (cm->ndtype[])
 */
int
TotalStatesInNode(int ndtype)
{
  switch (ndtype) {
  case BIF_nd:   return 1;
  case MATP_nd:  return 6;
  case MATL_nd:  return 3;
  case MATR_nd:  return 3;
  case BEGL_nd:  return 1;
  case BEGR_nd:  return 2;
  case ROOT_nd:  return 3;
  case END_nd:   return 1;
  default:       Die("Bogus node type %d", ndtype);
  }
  return 0;/*NOTREACHED*/
}
int
SplitStatesInNode(int ndtype)
{
  switch (ndtype) {
  case BIF_nd:   return 1;
  case MATP_nd:  return 4;
  case MATL_nd:  return 2;
  case MATR_nd:  return 2;
  case BEGL_nd:  return 1;
  case BEGR_nd:  return 1;
  case ROOT_nd:  return 1;
  case END_nd:   return 1;
  default:       Die("Bogus node type %d", ndtype);
  }
  return 0;/*NOTREACHED*/
}
int
InsertStatesInNode(int ndtype)
{
  switch (ndtype) {
  case BIF_nd:   return 0;
  case MATP_nd:  return 2;
  case MATL_nd:  return 1;
  case MATR_nd:  return 1;
  case BEGL_nd:  return 0;
  case BEGR_nd:  return 1;
  case ROOT_nd:  return 2;
  case END_nd:   return 0;
  default:       Die("Bogus node type %d", ndtype);
  }
  return 0;/*NOTREACHED*/
}




/* Function:  StateDelta(), StateLeftDelta(), StateRightDelta()
 * Incept:    SRE, Thu Oct  9 11:23:13 2003 [St. Louis]
 *
 * Purpose:   Convenience functions, mirroring some notation in Durbin et al.
 *            and elsewhere. \Delta notation simplifies some expositions
 *            of dynamic programming code.
 *            
 *            \Delta^R_v = 1 if the state emits right; else 0
 *            \Delta^L_v = 1 if the state emits left;  else 0
 *            \Delta_v   = 2 for pairwise, 1 for singlet, 0 for mute states.
 *            
 *            B_st, EL_st are special cases - Delta is returned as zero,
 *            but can't be used the same way.                                
 *
 * Args:      sttype   - state type code, e.g. MP_st
 *
 * Returns:   (see above)
 */
int
StateDelta(int sttype)
{
  switch (sttype) {
  case D_st:  return 0;
  case MP_st: return 2;
  case ML_st: return 1;
  case MR_st: return 1;
  case IL_st: return 1;
  case IR_st: return 1;
  case S_st:  return 0;
  case E_st:  return 0;
  case B_st:  return 0;
  case EL_st: return 0;
  default: Die("bogus state type %d\n", sttype);
  }
}
int
StateLeftDelta(int sttype)
{
  switch (sttype) {
  case D_st:  return 0;
  case MP_st: return 1;
  case ML_st: return 1;
  case MR_st: return 0;
  case IL_st: return 1;
  case IR_st: return 0;
  case S_st:  return 0;
  case E_st:  return 0;
  case B_st:  return 0;
  case EL_st: return 0;
  default: Die("bogus state type %d\n", sttype);
  }
}
int
StateRightDelta(int sttype)
{
  switch (sttype) {
  case D_st:  return 0;
  case MP_st: return 1;
  case ML_st: return 0;
  case MR_st: return 1;
  case IL_st: return 0;
  case IR_st: return 1;
  case S_st:  return 0;
  case E_st:  return 0;
  case B_st:  return 0;
  case EL_st: return 0;
  default: Die("bogus state type %d\n", sttype);
  }
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
  
  fprintf(fp, "Summary report for CM structure:\n");
  fprintf(fp, "--------------------------------------\n");
  fprintf(fp, "Total states:       %d\n", cm->M);
  fprintf(fp, "Total nodes:        %d\n", cm->nodes);
  fprintf(fp, "Bifurcations:       %d\n", count[BIF_B]);
  fprintf(fp, "MATP nodes:         %d\n", count[MATP_MP]);
  fprintf(fp, "MATL nodes:         %d\n", count[MATL_ML]);
  fprintf(fp, "MATR nodes:         %d\n", count[MATR_MR]);
  fprintf(fp, "Consensus columns:  %d    (2*MATP+MATL+MATR)\n",
	  count[MATP_MP]*2+count[MATL_ML]+count[MATR_MR]);
  fprintf(fp, "Base pairs:         %d    (MATP)\n", count[MATP_MP]);
  fprintf(fp, "Single stranded:    %d    (MATL+MATR)\n", count[MATL_ML]+count[MATR_MR]);

}

/* Functions: Statetype(), Nodetype(), UniqueStatetype();
 *            StateCode(), NodeCode(), UniqueStateCode()
 * Date:      SRE, Sat Jul 29 11:07:47 2000 [St. Louis]
 *
 * Purpose:   Translate internal flags into human-readable strings, 
 *            for clearer debugging output (*type functions);
 *            or vice versa (*Code functions)
 * 
 * Args:      type - a state type, node type, or unique statetype code
 *            s    - string representing a code
 *
 * Returns:   the appropriate string; or the appropriate code.
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
  case EL_st: return "EL";
  default: Die("bogus state type %d\n", type);
  }
  return ""; /*NOTREACHED*/
}
int
StateCode(char *s)
{
  if      (strcmp(s, "D")  == 0) return D_st;
  else if (strcmp(s, "MP") == 0) return MP_st;
  else if (strcmp(s, "ML") == 0) return ML_st;
  else if (strcmp(s, "MR") == 0) return MR_st;
  else if (strcmp(s, "IL") == 0) return IL_st;
  else if (strcmp(s, "IR") == 0) return IR_st;
  else if (strcmp(s, "S")  == 0) return S_st;
  else if (strcmp(s, "E")  == 0) return E_st;
  else if (strcmp(s, "B")  == 0) return B_st;
  else if (strcmp(s, "EL") == 0) return EL_st;
  return -1;
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
int
NodeCode(char *s)
{
  if      (strcmp(s, "BIF")  == 0) return BIF_nd;
  else if (strcmp(s, "MATP") == 0) return MATP_nd;
  else if (strcmp(s, "MATL") == 0) return MATL_nd;
  else if (strcmp(s, "MATR") == 0) return MATR_nd;
  else if (strcmp(s, "BEGL") == 0) return BEGL_nd;
  else if (strcmp(s, "BEGR") == 0) return BEGR_nd;
  else if (strcmp(s, "ROOT") == 0) return ROOT_nd;
  else if (strcmp(s, "END")  == 0) return END_nd;
  return -1;
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
  case END_EL : return "END_EL";
  default: Die("bogus unique state type %d\n", type);
  }
  return "";
}
int
UniqueStateCode(char *s)
{
  if      (strcmp(s, "ROOT_S")  == 0) return ROOT_S;
  else if (strcmp(s, "ROOT_IL") == 0) return ROOT_IL;  
  else if (strcmp(s, "ROOT_IR") == 0) return ROOT_IR;  
  else if (strcmp(s, "BEGL_S")  == 0) return BEGL_S;  
  else if (strcmp(s, "BEGR_S")  == 0) return BEGR_S;  
  else if (strcmp(s, "BEGR_IL") == 0) return BEGR_IL;  
  else if (strcmp(s, "MATP_MP") == 0) return MATP_MP;  
  else if (strcmp(s, "MATP_ML") == 0) return MATP_ML;  
  else if (strcmp(s, "MATP_MR") == 0) return MATP_MR;  
  else if (strcmp(s, "MATP_D")  == 0) return MATP_D;  
  else if (strcmp(s, "MATP_IL") == 0) return MATP_IL;  
  else if (strcmp(s, "MATP_IR") == 0) return MATP_IR;  
  else if (strcmp(s, "MATL_ML") == 0) return MATL_ML;  
  else if (strcmp(s, "MATL_D")  == 0) return MATL_D;  
  else if (strcmp(s, "MATL_IL") == 0) return MATL_IL;  
  else if (strcmp(s, "MATR_MR") == 0) return MATR_MR;  
  else if (strcmp(s, "MATR_D")  == 0) return MATR_D;  
  else if (strcmp(s, "MATR_IR") == 0) return MATR_IR;  
  else if (strcmp(s, "BIF_B")   == 0) return BIF_B;
  else if (strcmp(s, "END_E")   == 0) return END_E;
  else if (strcmp(s, "END_EL")  == 0) return END_EL;
  else Die("bogus unique statetype %s\n", s);
  return 0; /*NOTREACHED*/
}
int
DeriveUniqueStateCode(int ndtype, int sttype)
{
  switch (ndtype) {
  case BIF_nd:   
    switch (sttype) {
    case B_st:  return BIF_B;
    default:    return -1;
    }
  case MATP_nd:  
    switch (sttype) {
    case D_st:  return MATP_D;
    case MP_st: return MATP_MP;
    case ML_st: return MATP_ML;
    case MR_st: return MATP_MR;
    case IL_st: return MATP_IL;
    case IR_st: return MATP_IR;
    default:    return -1;
    }
  case MATL_nd:  
    switch (sttype) {
    case D_st:  return MATL_D;
    case ML_st: return MATL_ML;
    case IL_st: return MATL_IL;
    default:    return -1;
    }
  case MATR_nd:  
    switch (sttype) {
    case D_st:  return MATR_D;
    case MR_st: return MATR_MR;
    case IR_st: return MATR_IR;
    default:    return -1;
    }
  case BEGL_nd:  
    switch (sttype) {
    case S_st:  return BEGL_S;
    default:    return -1;
    }
  case BEGR_nd:  
    switch (sttype) {
    case S_st:  return BEGR_S;
    case IL_st: return BEGR_IL;
    default:    return -1;
    }
  case ROOT_nd:  
    switch (sttype) {
    case S_st:  return ROOT_S;
    case IL_st: return ROOT_IL;
    case IR_st: return ROOT_IR;
    default:    return -1;
    }
  case END_nd:   
    switch (sttype) {
    case E_st:  return END_E;
    default:    return -1;
    }
  default: 
    return -1;
  }
}

/* Function: CMRebalance()
 * Date:     SRE, Mon Apr  8 11:40:46 2002 [St. Louis]
 *
 * Purpose:  Rebalance a CM tree to guarantee O(N^2 log N) memory in
 *           smallcyk.c's divide and conquer algorithm.
 * 
 *           Input: a CM that's numbered in preorder traversal: 
 *           visit root, visit left, visit right. (e.g., left
 *           child S always visited before right child S, 
 *           cfirst[w] < cnum[y], as produced by modelmaker.c).
 *           
 *           Output: a renumbered CM, in a modified preorder traversal:
 *           visit root, visit min weight child, visit max weight child,
 *           where weight is the # of extra CYK decks that'll need to
 *           be held in memory to calculate this subgraph.
 *           
 * Args:     cm - the old CM
 *
 * Returns:  A new CM. 
 *           Caller is responsible for free'ing this with FreeCM().
 */
CM_t *
CMRebalance(CM_t *cm)
{
  Nstack_t *pda;          /* stack used for traversing old CM */
  CM_t     *new;          /* new CM we're creating */
  int      *wgt;          /* # of extra CYK decks required to calc subgraphs */
  int      *newidx;       /* newidx[v] = old CM state v's new index in new CM */
  int       v, w, y,z;	  /* state indices in old CM */
  int       nv;		  /* state index in new CM */
  int       x;		  /* counter over transitions, residues, nodes */

  /* Create the new model. Copy information that's unchanged by
   * renumbering the CM.
   */
  new = CreateCM(cm->nodes, cm->M);
  new->name = sre_strdup(cm->name, -1);
  new->acc  = sre_strdup(cm->acc,  -1);
  new->desc = sre_strdup(cm->desc, -1);
  for (x = 0; x < Alphabet_size; x++) new->null[x] = cm->null[x];

  /* Calculate "weights" (# of required extra decks) on every B and S state.
   * Recursive rule here is: 1 + min(wgt[left], wgt[right]).
   */
  wgt = MallocOrDie(sizeof(int) * cm->M);
  for (v = cm->M-1; v >= 0; v--) 
    {
      if      (cm->sttype[v] == E_st) /* initialize unbifurcated segments with 1 */
	wgt[v] = 1; 
      else if (cm->sttype[v] == B_st) /* "cfirst"=left S child. "cnum"=right S child. */
	wgt[v] = 1 + MIN(wgt[cm->cfirst[v]], wgt[cm->cnum[v]]);
      else 
	wgt[v] = wgt[v+1];            /* all other states propagate up to S */
    }

  /* Now, preorder traverse the new CM. At each bifurcation, we want
   * to visit the S with minimum weight first. v is an index on the
   * old CM, and we hop it around using this traversal order and a
   * pushdown stack. nv is an index on the new CM, which just moves
   * in preorder traversal 0..cm->M-1.
   * 
   */
  v = 0;
  z = cm->M-1;
  pda = CreateNstack();
  newidx = MallocOrDie(sizeof(int) * cm->M);
  for (nv = 0; nv < cm->M; nv++)
    {    
      /* Keep a map of where the old states are going in new CM 
       * old state v becomes newidx[v] in the new model.
       * This is guaranteed to be a one to one map.
       */
      newidx[v] = nv;		

      /* Copy old v to new nv. 
       * First, the easy stuff, that's unaffected by renumbering.
       */
      new->sttype[nv] = cm->sttype[v];
      new->ndidx[nv]  = cm->ndidx[v];
      new->stid[nv]   = cm->stid[v];
      new->pnum[nv]   = cm->pnum[v];
      for (x = 0; x < MAXCONNECT; x++) {
	new->t[nv][x]   = cm->t[v][x];
	new->tsc[nv][x] = cm->t[v][x];
      }
      for (x = 0; x < Alphabet_size*Alphabet_size; x++) {
	new->e[nv][x] = cm->e[v][x];
	new->esc[nv][x] = cm->esc[v][x];
      }

      /* Slightly harder - the plast connection for nv, to the last
       * of 1-6 parent states. We use the newidx map to get it from plast[v].
       */
      if (nv != 0) new->plast[nv] = newidx[cm->plast[v]];
      else         new->plast[nv] = -1;	/* ROOT. */

      /* Now, figure out next v, and make cfirst, cnum connections.
       * 
       * If we're a B, then traverse to the lighter child S state first.
       * Remember the overload in CM struct: cfirst = idx of left child; 
       * cnum = idx of right child. So if we visit left w first, cfirst=nv+1; 
       * if we visit right y first, cnum=nv+1. Getting the second child
       * index is a little tricky: we rely on knowing that 
       * the # of states in the first subgraph we visit is y-w,
       * so we know the second child index is nv+y-w+1.
       * 
       * If we're an E, pop the next v off the stack. cfirst=-1,cnum=0, because
       * it has no children.
       * 
       * Else, the next v is just v++. cfirst for new nv can be calculated by using the
       * offset in the old model: e.g. nv + (cfirst[v] - v). cnum is unchanged.
       * 
       */
      if (cm->sttype[v] == B_st) 
	{
	  w = cm->cfirst[v];	/* left child of v*/
	  y = cm->cnum[v];	/* right child of v*/

	  if (wgt[w] < wgt[y])	/* left (w) lighter? visit w first, defer y */
	    { 
	      PushNstack(pda, y); 
	      PushNstack(pda, z);
	      v = w; 
	      z = y-1;
	      new->cfirst[nv] = nv+1;     /* left child is nv+1 */
	      new->cnum[nv]   = nv+y-w+1; 
	    }  
	  else			/* right (y) lighter? visit y first, defer w */
	    { 
	      PushNstack(pda, w); 
	      PushNstack(pda, y-1);
	      v = y;		/* z unchanged. */
	      new->cfirst[nv] = nv+z-y+2; 
	      new->cnum[nv]   = nv+1;     /* right child is nv+1 */
	    }
	}
      else if (cm->sttype[v] == E_st) 
	{
	  new->cfirst[nv] = -1;
	  new->cnum[nv]   = 0;
	  PopNstack(pda, &z);
	  PopNstack(pda, &v);
	}
      else	
	{
	  new->cfirst[nv] = nv + (cm->cfirst[v]-v); /* use offset in old model */
	  new->cnum[nv]   = cm->cnum[v];            /* cnum unchanged. */
	  v++;
	}
    }

  /* Deal with the renumbered begin and end transition distributions,
   * using the newidx[v] map.
   */
  for (v = 0; v < cm->M; v++)
    {
      new->begin[newidx[v]] = cm->begin[v];
      new->end[newidx[v]]   = cm->end[v];
    }

  /* Guide tree numbering is unchanged - still in preorder.
   * Associate nodes with new state numbering.
   */
  for (x = 0; x < new->nodes; x++) 
    {
      new->nodemap[x] = newidx[cm->nodemap[x]];
      new->ndtype[x]  = cm->ndtype[x];
    }

  free(wgt);
  free(newidx);
  FreeNstack(pda);
  return new;
}

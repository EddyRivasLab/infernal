/* The core CM data structure.
 * 
 * Based on SRE HMMER3:p7_hmm.c
 *
 * Contents:
 *   1. The CM_CORE object: allocation, initialization, destruction.
 *   2. Convenience routines for setting fields in a CM.
 *   3. Convenience routines for getting information about a CM.
 *   4. Renormalization and rescaling counts in core CMs.
 *   5. Debugging and development code.
 *   6. Other routines in the API.
 *   7. Unit tests.
 *   8. Test driver. 
 *   9. Copyright and license.
 * 
 * EPN, Wed Aug  1 08:24:13 2007 [DC]
 * based on HMMER3's p7_hmm.c (SRE, Mon Jan  1 16:20:29 2007)
 * SVN $Id$
 */

#include "cm_config.h"		/* must be included first */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_vectorops.h"
#include "esl_random.h"
#include "esl_dirichlet.h"

#include "infernal.h"


/*****************************************************************
 * 1. The CM_CORE object: allocation, initialization, destruction.
 *****************************************************************/
/* Function: cm_core_Create(); 
 * Incept:   SRE, Sat Jul 29 09:02:16 2000 [St. Louis]
 *           
 * Purpose:  Allocate a <CM_CORE>, given the number of states 
 *           and nodes that should be in it, for the alphabet
 *           <abc>, and return a pointer to it.
 *
 *           The CM only keeps a copy of the <abc> alphabet
 *           pointer. The caller is responsible for providing the
 *           alphabet, keeping it around while the CM is in use,
 *           and (eventually) free'ing the alphabet when it's
 *           not needed any more. (Basically, just a step removed
 *           from keeping the alphabet as a global.)
 *
 * Args:     nnodes  =  number of nodes in the model
 *           nstates = number of states in the model
 *           clen    = consensus length of the model
 *           abc     = pointer to alphabet to use
 *
 * Returns:  ptr to allocated cm. 
 *           Caller is responsible for free'ing the cm.
 *
 * Throws:   <NULL> on allocation failure.
 */
CM_CORE *
cm_core_Create(int nnodes, int nstates, int clen, const ESL_ALPHABET *abc) 
{
  CM_CORE *cm = NULL;

  if ((cm = cm_core_CreateShell()) == NULL) return NULL;
  cm_core_CreateBody(cm, nnodes, nstates, clen, abc);
  return cm;
}  

/* Function:  cm_core_CreateShell()
 * Incept:    SRE, Sat Jul 29 09:02:16 2000 [St. Louis]
 *
 * Purpose:   Allocate the shell of a <CM_CORE>: everything that
 *            doesn't depend on knowing the number of states/nodes.
 *            
 *            CM input (<cmio.c>) uses two-step shell/body
 *            allocation because it has to read for a ways from the
 *            CM file before it reads the model size or the
 *            alphabet type.
 *
 * Returns:   a pointer to the new <CM_CORE> on success.
 *
 * Throws:    <NULL> on allocation failure.
 */
CM_CORE *
cm_core_CreateShell(void) 
{
  CM_CORE *cm = NULL;
  int     status;

  ESL_ALLOC(cm, sizeof(CM_CORE));

				/* general information: added later */
  cm->name     = NULL;
  cm->cs       = NULL;
  cm->acc      = NULL;
  cm->desc     = NULL;
  cm->rf       = NULL;
  cm->comlog   = NULL; 
  cm->nseq     = 0;
  cm->eff_nseq = 0.;
  cm->ctime    = NULL;
  cm->checksum = 0;

				/* structural information */
  cm->M      = 0;
  cm->clen   = 0;
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
  cm->el     = 0.;

  cm->ga = 0.;
  cm->tc = 0.;
  cm->nc = 0.;

  cm->offset   = 0;
  cm->flags    = 0;
  cm->abc      = NULL;
  cm->bg       = NULL;
  return cm;

 ERROR:
  return NULL;
}  

/* Function:  cm_core_CreateBody()
 * Incept:    SRE, Sat Jul 29 09:02:16 2000 [St. Louis]
 *
 * Purpose:   Given an allocated shell <cm>, and a now-known number
 *            of nodes and states, and alphabet <abc>, allocate
 *            the remainder of it for that many states/nodes.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; in this case, the CM
 *            is likely corrupted, and the caller should destroy it.
 */
int
cm_core_CreateBody(CM_CORE *cm, int nnodes, int nstates, int clen, 
		 const ESL_ALPHABET *abc) 
{
  int k;
  int status;

  cm->abc = abc;
  cm->clen= clen;
  cm->M   = nstates;

				/* structural information */
  ESL_ALLOC(cm->sttype, (nstates+1) * sizeof(char));
  ESL_ALLOC(cm->nididx,  nstates    * sizeof(int));
  ESL_ALLOC(cm->stid,   (nstates+1) * sizeof(char));
  ESL_ALLOC(cm->cfirst,  nstates    * sizeof(int));
  ESL_ALLOC(cm->cnum,    nstates    * sizeof(int));
  ESL_ALLOC(cm->plast,   nstates    * sizeof(int));
  ESL_ALLOC(cm->pnum,    nstates    * sizeof(int));

				/* node->state map information */
  cm->nodes  = nnodes;
  ESL_ALLOC(cm->nodemap, nnodes  * sizeof(int));
  ESL_ALLOC(cm->ndtype,  nnodes  * sizeof(int));

				/* parameter information */
  /* level 1 */
  ESL_ALLOC(cm->t,    (nstates) * sizeof(float *));
  ESL_ALLOC(cm->e,    (nstates) * sizeof(float *));
  cm->t[0]   = NULL;
  cm->e[0]   = NULL;

  /* level 2 */
  ESL_ALLOC(cm->t[0],   (cmH_MAXTRANSITIONS*(nstates))      * sizeof(float));
  ESL_ALLOC(cm->e[0],   (cm->abc->K * cm->abc->K*(nstates)) * sizeof(float));
  for (v = 0; v < nstates; v++) {
    cm->e[v] = cm->e[0] + v * (cm->abc->K * cm->abc->K);
    cm->t[v] = cm->t[0] + v * cmH_MAXTRANSITIONS;
  }

  /* the EL state at M is special: we only need state
   * type info recorded, so functions looking at parsetrees  
   * can interpret what an "M" index means.
   */
  cm->sttype[cm->M] = EL_st;
  cm->stid[cm->M]   = END_EL;

  cm->flags         = 0;

  if ((status = cm_core_Zero(cm)) != eslOK) goto ERROR;

  ESL_ALLOC(cm->cs,  (cm->clen+2) * sizeof(char));
  /* Optional allocation, status flag dependent */
  if (cm->flags & cmH_RF)  ESL_ALLOC(cm->rf,  (cm->clen+2) * sizeof(char));
  
  return eslOK;

 ERROR:
  return status;
}  


/* Function:  cm_core_Destroy()
 * Incept:    SRE, Sat Jul 29 11:22:32 2000 [St. Louis]
 *
 * Purpose:   Frees both the shell and body of an <cm>.
 *            Works even if the <cm> is damaged (incompletely allocated)
 *            or even <NULL>.
 *
 * Note:      Remember, leave reference pointers like abc, gm, and
 *            bg alone. These are under the application's control not ours.
 *
 * Returns:   (void).
 */
void
cm_core_Destroy(CM_CORE *cm)
{
  if (cm == NULL) return;

  if (cm->e     != NULL) {
    if (cm->e[0] != NULL) free(cm->e[0]);
    free(cm->e);
  }
  if (cm->t != NULL) {
    if (cm->t[0] != NULL) free(cm->t[0]);
    free(cm->t);
  }

  if (cm->name    != NULL) free(cm->name);
  if (cm->acc     != NULL) free(cm->acc);
  if (cm->desc    != NULL) free(cm->desc);
  if (cm->rf      != NULL) free(cm->rf);
  if (cm->cs      != NULL) free(cm->cs);
  if (cm->comlog  != NULL) free(cm->comlog);
  if (cm->ctime   != NULL) free(cm->ctime);

  if(cm->sttype != NULL)   free(cm->sttype);
  if(cm->ndidx  != NULL)   free(cm->ndidx);
  if(cm->stid   != NULL)   free(cm->stid);
  if(cm->cfirst != NULL)   free(cm->cfirst);
  if(cm->cnum   != NULL)   free(cm->cnum);
  if(cm->plast  != NULL)   free(cm->plast);
  if(cm->pnum   != NULL)   free(cm->pnum);
  if(cm->nodemap!= NULL)   free(cm->nodemap);
  if(cm->ndtype != NULL)   free(cm->ndtype);

  free(cm);
  return;
}

/* Function:  cm_core_CopyParameters()
 * Incept:    SRE, Fri May  4 14:10:17 2007 [Janelia]
 *            Infernalized: EPN, Wed Aug  1 09:25:09 2007
 *
 * Purpose:   Copy parameters of <src> to <dest>. The CM <dest> must
 *            be allocated by the caller for the same 
 *            alphabet and M as <src>. 
 *            
 *            No annotation is copied.  This is because several
 *            annotation fields are variable-length strings that
 *            require individual allocations.  The
 *            <cm_core_CopyParameters()> function is for cases where we
 *            have to repeatedly reset the parameters of a model - for
 *            example, in entropy weighting.
 *
 * Returns:   <eslOK> on success.
 */
int
cm_core_CopyParameters(const CM_CORE *src, CM_CORE *dest)
{
  int k;
  for (k = 0; k < src->M; k++) {
    esl_vec_FCopy(src->t[k],   cmH_NTRANSITIONS,           dest->t[k]);
    esl_vec_FCopy(src->e[k],   (src->abc->K * src->abc->K),dest->e[k]);
  }
  return eslOK;
}

/* Function:  cm_core_Duplicate()
 * Incept:    SRE, Fri Jan 26 15:34:42 2007 [Janelia]
 *
 * Purpose:   Duplicates a cm.
 * 
 *            Note: does not duplicate the objects the CM refers to,
 *            if any (profile, null model, or alphabet); only copies
 *            the reference pointers.
 * 
 * Returns:   a pointer to the duplicate.
 *
 * Throws:    <NULL> on allocation failure.
 */
CM_CORE *
cm_core_Duplicate(const CM_CORE *cm)
{
  int     status;
  CM_CORE *new = NULL;

  if ((new = cm_core_Create(cm->M, cm->nodes, cm->clen, cm->abc)) == NULL) goto ERROR;
  cm_core_CopyParameters(cm, new);

  esl_vec_ICopy(cm->sttype, (cm->M+1), new->sttype);
  esl_vec_ICopy(cm->nididx,  cm->M,    new->ndidx);
  esl_vec_ICopy(cm->stid,   (cm->M+1), new->stid);
  esl_vec_ICopy(cm->cfirst,  cm->M,    new->cfirst);
  esl_vec_ICopy(cm->cnum,    cm->M,    new->cnum);
  esl_vec_ICopy(cm->plast,   cm->M,    new->plast);
  esl_vec_ICopy(cm->pnum,    cm->M,    new->pnum);

  esl_vec_ICopy(cm->nodemap,   cm->nodes,  new->nodemap);
  esl_vec_ICopy(cm->ndtype,    cm->nodes,  new->ndtype);
  
  if (cm->name != NULL    && (status = esl_strdup(cm->name,   -1, &(new->name)))   != eslOK) goto ERROR;
  if (cm->acc  != NULL    && (status = esl_strdup(cm->acc,    -1, &(new->acc)))    != eslOK) goto ERROR;
  if (cm->cs   != NULL    && (status = esl_strdup(cm->cs,     -1, &(new->cs)))     != eslOK) goto ERROR;
  if (cm->desc != NULL    && (status = esl_strdup(cm->desc,   -1, &(new->desc)))   != eslOK) goto ERROR;
  if (cm->flags & cmH_RF  && (status = esl_strdup(cm->rf,     -1, &(new->rf)))     != eslOK) goto ERROR;
  if (cm->comlog != NULL  && (status = esl_strdup(cm->comlog, -1, &(new->comlog))) != eslOK) goto ERROR;
  if (cm->ctime  != NULL  && (status = esl_strdup(cm->ctime,  -1, &(new->ctime)))  != eslOK) goto ERROR;

  new->nseq     = cm->nseq;
  new->eff_nseq = cm->eff_nseq;
  new->checksum = cm->checksum;
  new->ga       = cm->ga;
  new->tc       = cm->tc;
  new->nc       = cm->nc;
  new->offset   = cm->offset;
  new->flags    = cm->flags;
  new->abc      = cm->abc;
  return new;

 ERROR:
  if (new != NULL) cm_core_Destroy(new);
  return NULL;
}

/* Function:  cm_core_Scale()
 * Incept:    SRE, Fri May  4 14:19:33 2007 [Janelia]
 *            Infernalized: EPN, Wed Aug  1 09:35:26 2007
 *
 * Purpose:   Given a counts-based model <cm>, scale core
 *            by a multiplicative factor of <scale>. Used in
 *            absolute sequence weighting or effective
 *            sequence weighting.
 *
 * Returns:   <eslOK> on success.
 */
int
cm_core_Scale(CM_CORE *cm, double scale)
{
  int k;

  for (v = 0; v < cm->M; v++) {
    esl_vec_FScale(cm->t[v], cm->cnum[v],                             scale);  
    esl_vec_FScale(cm->e[v], cm_core_NEmitAlph(cm->abc, cm->sttype[v]), scale);  
  }
  return eslOK;
}


/* Function:  cm_core_Zero()
 * Incept:    SRE, Mon Jan  1 16:32:59 2007 [Casa de Gatos]
 *            Infernalized: EPN, Wed Aug  1 09:35:31 2007
 *
 * Purpose:   Zeroes the counts/probabilities fields in core model.
 *
 * Returns:   <eslOK> on success.
 */
int
cm_core_Zero(CM_CORE *cm)
{
  int k;

  for (k = 0; k <= cm->M; k++) {
    esl_vec_FSet(cm->t[k],   cmH_NTRANSITIONS,       0.);  
    esl_vec_FSet(cm->e[k],  (cm->abc->K*cm->abc->K), 0.);  
  }
  return eslOK;
}



/* Function:  cm_core_DescribeStatetype()
 * Incept:    SRE, Mon Jan  1 18:47:34 2007 [Casa de Gatos]
 *            Infernalized: EPN, Wed Aug  1 09:46:29 2007
 *
 * Purpose:   Returns the state type in text, as a string of length 2 
 *            (3 if you count NUL). States with only 1 letter have an
 *            'S' concatenated to the end. For example, <cm_Statetype(D_st)>
 *            returns "DS", and <cm_Statetype(IL_st)> returns "IL".
 */
char *
cm_core_DescribeStatetype(int st)
{
  switch (st) {
  case cmT_M: return "DS";
  case cmT_D: return "MP";
  case cmT_I: return "ML";
  case cmT_S: return "MR";
  case cmT_N: return "IL";
  case cmT_B: return "IR";
  case cmT_E: return "SS";
  case cmT_C: return "ES";
  case cmT_T: return "BS";
  case cmT_J: return "EL";
  default:     return "?S";
  }
}




/*****************************************************************
 * 2. Convenience routines for setting fields in an CM.
 *****************************************************************/ 

/* Function: cm_core_SetName()
 * Incept:   SRE, Mon Jan  1 16:53:23 2007 [Casa de Gatos]
 *           Infernalized: EPN, Wed Aug  1 09:43:55 2007
 * 
 * Purpose:  Set or change the name of a CM to <name>.
 *           Any trailing whitespace (including newline) is chopped off.     
 *      
 * Returns:  <eslOK> on success.
 *
 * Throws:   <eslEMEM> on allocation error, and original name (if any) 
 *           remains.
 */
int
cm_core_SetName(CM_CORE *cm, char *name)
{
  int   status;
  void *tmp;
  int   n;

  if (name == NULL) {
    if (cm->name != NULL) free(cm->name); 
    cm->name = NULL;
  } else {
    n = strlen(name);
    ESL_RALLOC(cm->name, tmp, sizeof(char)*(n+1));
    strcpy(cm->name, name);
    if ((status = esl_strchop(cm->name, n)) != eslOK) goto ERROR;
  }
  return eslOK;

 ERROR:
  return status;
}

/* Function: cm_core_SetAccession()
 * Incept:   SRE, Mon Jan  1 16:53:53 2007 [Casa de Gatos]
 *           Infernalized: EPN, Wed Aug  1 09:44:04 2007
 * 
 * Purpose:  Set or change the accession number of a CM to <acc>,
 *           and raise the <CM_ACC> flag. Trailing whitespace (including newline) 
 *           is chopped.  
 *           
 *           If <acc> is <NULL>, unset the CM's accession (if any) and drop 
 *           the <CM_ACC> flag.
 *
 * Returns:  <eslOK> on success.
 *
 * Throws:   <eslEMEM> on allocation error, and original name (if any) 
 *           remains.
 */
int
cm_core_SetAccession(CM_CORE *cm, char *acc)
{
  int   status;
  void *tmp;
  int   n;

  if (acc == NULL) {
    if (cm->acc != NULL) free(cm->acc); 
    cm->acc = NULL;
    cm->flags &= ~cmH_ACC;
  } else {
    n = strlen(acc);
    ESL_RALLOC(cm->acc, tmp, sizeof(char)*(n+1));
    strcpy(cm->acc, acc);
    if ((status = esl_strchop(cm->acc, n)) != eslOK) goto ERROR;
    cm->flags |= cmH_ACC;
  }
  return eslOK;

 ERROR:
  return status;
}

/* Function: cm_core_SetDescription()
 * Incept:   SRE, Mon Jan  1 16:59:28 2007 [Casa de Gatos]
 *           Infernalized: EPN, Wed Aug  1 09:44:11 2007
 * 
 * Purpose:  Set or change the description line of a CM. 
 *           Trailing whitespace (including newline) is chopped.
 */
int
cm_core_SetDescription(CM_CORE *cm, char *desc)
{
  int   status;
  void *tmp;
  int   n;

  if (desc == NULL) 
    {
      if (cm->desc != NULL) free(cm->desc); 
      cm->desc   = NULL;
      cm->flags &= ~cmH_DESC;
    }
  else
    {
      n = strlen(desc);
      ESL_RALLOC(cm->desc, tmp, sizeof(char)*(n+1));
      strcpy(cm->desc, desc);
      if ((status = esl_strchop(cm->desc, n)) != eslOK) goto ERROR;
      cm->flags |= cmH_DESC;
    }
  return eslOK;

 ERROR:
  return status;
}

/* Function: cm_core_AppendComlog()
 * Incept:   SRE, Mon Jan  1 18:23:42 2007 [Casa de Gatos] 
 *           Infernalized: EPN, Wed Aug  1 09:44:18 2007
 * 
 * Purpose:  Concatenate command line options and append as a new line in the
 *           command line log. Command line log is multiline, with each line
 *           ending in newline char, except for last line.
 *           
 * Returns:  <eslOK> on success.
 * 
 * Throws:   <eslEMEM> on allocation failure.          
 */
int
cm_core_AppendComlog(CM_CORE *cm, int argc, char **argv)
{
  int   status;
  void *tmp;
  int   n;
  int   i;

  /* figure out length of added command line, and (re)allocate comlog */
  n = argc-1;	/* account for 1 space per arg, except last one */
  for (i = 0; i < argc; i++)
    n += strlen(argv[i]);

  if (cm->comlog != NULL) {
    n += strlen(cm->comlog) + 1; /* +1 for the \n we're going to add to the old comlog */
    ESL_RALLOC(cm->comlog, tmp, sizeof(char)* (n+1));
    strcat(cm->comlog, "\n");
  } else {
    ESL_ALLOC(cm->comlog, sizeof(char)* (n+1));
    *(cm->comlog) = '\0'; /* need this to make strcat work */
  }

  for (i = 0; i < argc-1; i++)
    {
      strcat(cm->comlog, argv[i]);
      strcat(cm->comlog, " ");
    }
  strcat(cm->comlog, argv[argc-1]);
  return eslOK;

 ERROR:
  return status;
}

/* Function: cm_core_SetCtime()
 * Incept:   SRE, Wed Oct 29 11:53:19 1997 [TWA 721 over the Atlantic]
 *           Infernalized: EPN, Wed Aug  1 09:44:26 2007
 * 
 * Purpose:  Set the <ctime> field in a new CM to the current time.
 *
 *           This function is not reentrant and not threadsafe, because
 *           it calls the nonreentrant ANSI C ctime() function.
 * 
 * Returns:  <eslOK> on success.
 * 
 * Throws:   <eslEMEM> on allocation failure. <eslESYS> if the <time()>
 *           system call fails to obtain the calendar time.
 */
int
cm_core_SetCtime(CM_CORE *cm)
{
  int    status;
  char  *s = NULL;
  time_t date;

  if ((date   = time(NULL))                       == -1) { status = eslESYS; goto ERROR; }
  if ((status = esl_strdup(ctime(&date), -1, &s)) != eslOK) goto ERROR;
  if ((status = esl_strchop(s, -1))               != eslOK) goto ERROR;
  
  if (cm->ctime != NULL) free(cm->ctime);
  cm->ctime = s;
  return eslOK;

 ERROR:
  if (s != NULL) free(s);
  return status;
}
/*---------------- end, internal-setting routines ---------------*/

/*****************************************************************
 * 3. Convenience routines for getting information about a CM.
 *****************************************************************/ 

/* Function: cm_core_CountStatetype(), cm_core_SubtreeCountStatetype(), 
 *           cm_core_SegmentCountStatetype
 * Date:     SRE, Wed Aug  2 09:15:00 2000 [St. Louis]
 *
 * Purpose:  Conveniences for counting the # of occurrences
 *           of a particular state type in a CM. Useful for
 *           "how many bifurcations does this model have", etc.
 *          
 *           cm_core_SubtreeCountStatetype() only counts underneath     
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
cm_core_SegmentCountStatetype(CM_CORE *cm, int r, int z, char type)
{
  int count = 0;
  int v;
  for (v = r; v <= z; v++) 
    if (cm->sttype[v] == type) count++;
  return count;
}
int
cm_core_SubtreeCountStatetype(CM_CORE *cm, int v, char type)
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
cm_core_CountStatetype(CM_CORE *cm, char type)
{
  return CMSubtreeCountStatetype(cm, 0, type);
}
int 
cm_core_SubtreeFindEnd(CM_CORE *cm, int r)
{
  int unsatisfied_starts = 1;

  while (unsatisfied_starts) {
    if (cm->sttype[r] == B_st) unsatisfied_starts++;
    if (cm->sttype[r] == E_st) unsatisfied_starts--; 
    r++;
  }
  return (r-1);
}

/* Function: cm_core_CalculateStateIndex()
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
cm_core_CalculateStateIndex(CM_CORE *cm, int node, char utype)
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
  default: esl_fatal("bogus utype %d in CalculateStateIndex()", utype);
  }
  return base;			/* not used */
}

/* Function:  cm_core_TotalStatesInNode(), cm_core_SplitStatesInNode(), 
 *            cm_core_InsertStatesInNode()
 * Incept:    SRE, Thu Aug  8 09:57:59 2002 [St. Louis]
 *
 * Purpose:   Returns the number of states in a node type.
 *
 * Args:      ndtype  - type of node (cm->ndtype[])
 */
int
cm_core_TotalStatesInNode(int ndtype)
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
  default:       esl_fatal("Bogus node type %d", ndtype);
  }
  return 0;/*NOTREACHED*/
}
int
cm_core_SplitStatesInNode(int ndtype)
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
  default:       esl_fatal("Bogus node type %d", ndtype);
  }
  return 0;/*NOTREACHED*/
}
int
cm_core_InsertStatesInNode(int ndtype)
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
  default:       esl_fatal("Bogus node type %d", ndtype);
  }
  return 0;/*NOTREACHED*/
}

/* Function:  cm_core_StateDelta(), cm_core_StateLeftDelta(), cm_core_StateRightDelta()
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
cm_core_StateDelta(int sttype)
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
  default: esl_fatal("bogus state type %d\n", sttype);
  }
  /*NOTREACHED*/
  return 0;
}
int
cm_core_StateLeftDelta(int sttype)
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
  default: esl_fatal("bogus state type %d\n", sttype);
  }
  /*NOTREACHED*/
  return 0;
}
int
cm_core_StateRightDelta(int sttype)
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
  default: esl_fatal("bogus state type %d\n", sttype);
  }
  /*NOTREACHED*/
  return 0;
}
/* Function:  cm_core_NEmitAlph()
 * Incept:    EPN, Wed Aug  1 10:26:09 2007
 *
 * Purpose:   Return the number of emissions for a given stateype
 */
int
cm_core_NEmitAlph(const ESL_ALPHABET *abc, int sttype)
{
  switch (sttype) {
  case D_st:  return 0;
  case MP_st: return abc->K * abc->K;
  case ML_st: return abc->K;
  case MR_st: return abc->K;
  case IL_st: return abc->K;
  case IR_st: return abc->K;
  case S_st:  return 0;
  case E_st:  return 0;
  case B_st:  return 0;
  case EL_st: return abc->K;
  default: esl_fatal("bogus state type %d\n", sttype);
  }
  /*NOTREACHED*/
  return 0;
}

/*****************************************************************
 * 4. Renormalization of core CMs.
 *****************************************************************/ 

/* Function: cm_core_Renormalize()
 * Incept:   SRE, Mon Jan  1 18:39:42 2007 [Casa de Gatos]
 * 
 * Purpose:  Take a core CM in counts form, and renormalize
 *           all probability vectors in the core probability model.
 *
 *           Leaves other flags (stats and profile) alone, so caller
 *           needs to be wary. Renormalizing a probability model that
 *           has stats and profile scores wouldn't usually invalidate
 *           those data; and if we're renormalizing a counts model, we
 *           shouldn't have stats or profile scores yet anyway.
 *           
 * Args:     cm - the model to renormalize.
 *                 
 * Return:   <eslOK> on success.
 */                          
int
cm_core_Renormalize(CM_CORE *cm)
{
  int   v;			/* counter for states */

  for (v = 0; v < cm->M; v++) {
    if((nemit = cm_core_NEmitAlph(cm->abc, cm->sttype[v])) > 0)
      esl_vec_FNorm(cm->e[k], nemit);
    if(cm->sttype[v] != B_st && cm->sttype[v] != E_st)
      esl_vec_FNorm(cm->t[v], cm->cnum[v]);
  }
  return eslOK;
}
  
/*****************************************************************
 * 4. Debugging and development code
 *****************************************************************/

/* Function:  cm_core_Dump()
 * Incept:    SRE, Mon Jan  1 18:44:15 2007 [Casa de Gatos]
 *
 * Purpose:   Debugging: dump the probabilities (or counts) from a core CM.
 * 
 * Returns:   <eslOK> on success.
 */
int
cm_core_Dump(FILE *fp, CM_CORE *cm)
{
  int v;			/* counter for states */
  int x;			/* counter for symbols */
  int ts;			/* counter for state transitions */
  
  for (k = 0; k <= cm->M; k++)
    {				/* Line 1: v, match emissions */
      fprintf(fp, " %5d ", v);
      for (x = 0; x < NEmitAlph(cm->abc, cm->sttype[v]); x++) 
        fprintf(fp, "%9.4f ", cm->e[v][x]);
      fputs("\n", fp);
				/* Line 2: transition probs */
      fprintf(fp, "       ");
      for (ts = 0; ts < cm->cnum[v]; ts++)
	fprintf(fp, "%9.4f ", cm->t[v][ts]); 
      fputs("\n", fp);
    }
  fputs("//\n", fp);
  return eslOK;
}

/* Function:  cm_core_Compare()
 * Incept:    SRE, Sat Jan  6 14:14:58 2007 [Casa de Gatos]
 *            Infernalized: EPN, Wed Aug  1 11:07:21 2007
 *
 * Purpose:   Compare two CMs <cm1> and <cm2> to each other;
 *            return <eslOK> if they're identical, and <eslFAIL>
 *            if they differ. Floating-point probabilities are 
 *            compared for equality within a fractional tolerance
 *            <tol>. 
 */
int
cm_core_Compare(CM_CORE *h1, CM_CORE *cm2, float tol)
{
  int k;
  
  if (cm1->abc->type != cm2->abc->type) return eslFAIL;
  if (cm1->M         != cm2->M)         return eslFAIL;
  if (cm1->nodes     != cm2->nodes)     return eslFAIL;
  if (cm1->clen      != cm2->clen)      return eslFAIL;
  if (cm1->flags     != cm2->flags)     return eslFAIL;

  if(esl_vec_ICompare(cm1->sttype, cm2->sttype, (cm1->M+1)) != eslOK) return eslFAIL;
  if(esl_vec_ICompare(cm1->ndidx,  cm2->ndidx,   cm1->M)    != eslOK) return eslFAIL;
  if(esl_vec_ICompare(cm1->stid,   cm2->stid,   (cm1->M+1)) != eslOK) return eslFAIL;
  /* next 4 lines are probably unnec, same stids should guarantee same connectivity */
  if(esl_vec_ICompare(cm1->cfirst, cm2->cfirst,  cm1->M)    != eslOK) return eslFAIL;
  if(esl_vec_ICompare(cm1->cnum,   cm2->cnum,    cm1->M)    != eslOK) return eslFAIL;
  if(esl_vec_ICompare(cm1->plast,  cm2->plast,   cm1->M)    != eslOK) return eslFAIL;
  if(esl_vec_ICompare(cm1->pnum,   cm2->pnum,    cm1->M)    != eslOK) return eslFAIL;

  if(esl_vec_ICompare(cm1->nodemap, cm2->nodemap, cm1->nodes) != eslOK) return eslFAIL;
  if(esl_vec_ICompare(cm1->nodemap, cm2->nodemap, cm1->nodes) != eslOK) return eslFAIL;
  if(esl_vec_ICompare(cm1->ndtype,  cm2->ndtype,  cm1->nodes) != eslOK) return eslFAIL;

  for (v = 0; v < cm1->M; v++)	/* (it's safe to include 0 here.) */
    {
      if (esl_vec_FCompare(cm1->e[k], cm2->e[k], cm_core_NEmitAlph(cm1->abc, cm1->sttype[v]), tol) != eslOK) return eslFAIL;
      if (esl_vec_FCompare(cm1->t[k], cm2->e[k], cm1->cnum[v], tol) != eslOK) return eslFAIL;
    }

  if (strcmp(cm1->name,   cm2->name)   != 0) return eslFAIL;
  if (strcmp(cm1->comlog, cm2->comlog) != 0) return eslFAIL;
  if (strcmp(cm1->ctime,  cm2->ctime)  != 0) return eslFAIL;
  if (cm1->nseq     != cm2->nseq)            return eslFAIL;
  if (cm1->eff_nseq != cm2->eff_nseq)        return eslFAIL;
  if (cm1->checksum != cm2->checksum)        return eslFAIL;

  if ((cm1->flags & cmH_CS)   && strcmp(cm1->cs,   cm2->cs)   != 0) return eslFAIL;
  if ((cm1->flags & cmH_ACC)  && strcmp(cm1->acc,  cm2->acc)  != 0) return eslFAIL;
  if ((cm1->flags & cmH_DESC) && strcmp(cm1->desc, cm2->desc) != 0) return eslFAIL;
  if ((cm1->flags & cmH_RF)   && strcmp(cm1->rf,   cm2->rf)   != 0) return eslFAIL;

  if (cm1->flags & cmH_GA) 
    if (esl_FCompare(cm1->ga, cm2->ga, tol) != eslOK) return eslFAIL;

  if (cm1->flags & cmH_TC) 
    if (esl_FCompare(cm1->tc, cm2->tc, tol) != eslOK) return eslFAIL;

  if (cm1->flags & cmH_NC) 
    if (esl_FCompare(cm1->nc, cm2->nc, tol) != eslOK) return eslFAIL;

  return eslOK;
}

/* Function:  cm_core_Validate()
 * Incept:    SRE, Sat Jan  6 14:43:00 2007 [Casa de Gatos]
 *           
 * Purpose:   Validates the internals of the CM structure <cm>.
 * 
 *            Probability vectors are validated to sum up to
 *            within a fractional tolerance <tol> of 1.0.
 *
 *            Probably only useful for debugging and development,
 *            not production code.
 *
 * Returns:   <eslOK> if <cm> internals look fine.
 *            Returns <eslFAIL> if something is wrong.
 */
int
cm_core_Validate(CM_CORE *cm, float tol, char *errbuf)
{
  int status;
  int k;

  if (cm            == NULL)       ESL_XFAIL(eslFAIL, errbuf, "CM is a null pointer");
  if (cm->M         <  1)          ESL_XFAIL(eslFAIL, errbuf, "CM has M < 1");
  if (cm->nodes     <  1)          ESL_XFAIL(eslFAIL, errbuf, "CM has nodes < 1");
  if (cm->abc       == NULL)       ESL_XFAIL(eslFAIL, errbuf, "CM has no alphabet reference");
  if (cm->abc->type == eslUNKNOWN) ESL_XFAIL(eslFAIL, errbuf, "CM's alphabet is set to unknown");
  
  for (v = 0; v < cm->M; k++)
    {
      if (esl_vec_FValidate(cm->e[v], cm_core_NEmitAlph(cm->abc, cm->sttype[v]), tol, NULL) != eslOK)
	    ESL_XFAIL(eslFAIL, errbuf, "e[%d] fails pvector validation", v); 
      if (cm->sttype[v] != B_st && cm->sttype[v] != E_st)
	if(esl_vec_FValidate(cm->t[v], cm->cnum[v], tol, NULL) != eslOK) 
	  ESL_XFAIL(eslFAIL, errbuf, "t[%d] fails pvector validation", v);
    }

  /* Don't be strict about mandatory name, comlog, ctime, cs for now in development */
  /*  if (cm->name     == NULL) return eslFAIL; */
  /*  if (cm->comlog   == NULL) return eslFAIL; */
  /*  if (cm->ctime    == NULL) return eslFAIL;  */
  /*  if (cm->cs       == NULL) return eslFAIL;  */
  if (cm->nseq     <  0 )   ESL_XFAIL(eslFAIL, errbuf, "invalid nseq");
  if (cm->eff_nseq <  0 )   ESL_XFAIL(eslFAIL, errbuf, "invalid eff_nseq");
  if (cm->checksum <  0 )   ESL_XFAIL(eslFAIL, errbuf, "invalid checksum");

  if (  (cm->flags & cmH_ACC)  && cm->acc  == NULL) ESL_XFAIL(eslFAIL, errbuf, "accession null but cmH_ACC flag is up");
  if (! (cm->flags & cmH_ACC)  && cm->acc  != NULL) ESL_XFAIL(eslFAIL, errbuf, "accession present but cmH_ACC flag is down");
  if (  (cm->flags & cmH_DESC) && cm->desc == NULL) ESL_XFAIL(eslFAIL, errbuf, "description null but cmH_DESC flag is up");
  if (! (cm->flags & cmH_DESC) && cm->desc != NULL) ESL_XFAIL(eslFAIL, errbuf, "description present but cmH_DESC flag is down");
  if (cm->flags & cmH_RF) {
    if (cm->rf == NULL || strlen(cm->rf) != cm->M+1) ESL_XFAIL(eslFAIL, errbuf, "cmH_RF flag up, but rf string is invalid");
  } else 
    if (cm->rf != NULL)                                ESL_XFAIL(eslFAIL, errbuf, "cmH_RF flag down, but rf string is present");

  return eslOK;

 ERROR:
  return status;
}
/*------------- end of debugging/development code ----------------*/

/*****************************************************************
 * 5. Test driver.
 *****************************************************************/

#ifdef cmCM_TESTDRIVE

#include <cm_config.h>
#include <infernal.h>

int
main(int argc, char **argv)
{
  utest_foo();
  exit(0); /* success */
}

#endif /*cmCM_TESTDRIVE*/
/*-------------------- end of test driver ---------------------*/


/************************************************************
 * @LICENSE@
 ************************************************************/


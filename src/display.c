/* display.c
 * SRE, Thu May 23 08:18:05 2002 [St. Louis]
 * CVS $Id$
 * 
 * Routines for formatting and displaying parse trees
 * for output.
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************  
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "squid.h"
#include "nstack.h"
#include "structs.h"
#include "funcs.h"

static int *createMultifurcationOrderChart(CM_t *cm);
static void createFaceCharts(CM_t *cm, int **ret_inface, int **ret_outface);

/* Flags used to define PDA moves
 */
#define PDA_RESIDUE 0
#define PDA_STATE   1
#define PDA_MARKER  2

/* Function:  CreateFancyAli()
 * Incept:    SRE, Thu May 23 13:46:09 2002 [St. Louis]
 *
 * Purpose:   Given a trace (and the model and sequence it corresponds
 *            to), create an alignment for display; return in a Fancyali_t
 *            structure.
 *
 * Args:      tr    - parsetree for cm aligned to dsq
 *            cm    - model
 *            cons  - consensus information for cm; see CreateCMConsensus()
 *            dsq   - digitized sequence
 *
 * Returns:   fancy alignment structure.
 *            Caller frees, with FreeFancyAli(ali).
 *
 * Xref:      STL6 p.58
 */
Fancyali_t *
CreateFancyAli(Parsetree_t *tr, CM_t *cm, CMConsensus_t *cons, char *dsq)
{
  Fancyali_t *ali;              /* alignment structure we're building, and returning */
  Nstack_t   *pda;              /* pushdown automaton used to traverse the trace. */
  int         type;		/* type of move on automaton: PDA_RESIDUE or PDA_STATE */
  int         v;		/* state index       */
  int         nd;		/* node index        */
  int         ti;		/* position in trace */
  int         ninset;		/* # nt in an EL     */
  int         pos;		/* position in growing ali */
  int         lc, rc;		/* indices for left, right pos in consensus */
  int         symi, symj;
  int         lannote, rannote; /* characters to put in annotation line; left, right */
  int         lstr, rstr;	/* chars to put in structure line; left, right */
  int         lcons, rcons;	/* chars to put in model consensus line; left, right */
  int         lmid, rmid;	/* chars to put in middle alignment quality line; left, right */
  int         lseq, rseq;	/* chars to put in aligned target seq line; left, right */
  int         do_left, do_right;/* flags to generate left, right */

  ali = MallocOrDie(sizeof(Fancyali_t));
  
  /* Calculate length of the alignment display.
   *   MATP node        : +2
   *   MATL, MATR node  : +1
   *   IL, IR state     : +1
   *   EL:              : 4 + width of length display : "*[nn]*"
   *   anything else    : 0.
   */
  ali->len = 0;
  for (ti = 0; ti < tr->n; ti++)
    {
      v  = tr->state[ti];
      if (v == cm->M) {  /* special case: local exit into EL */
	ninset     = tr->emitr[ti] - tr->emitl[ti] + 1;
	ali->len += 4;
	do { ali->len++; ninset /= 10; } while (ninset);
	continue;
      } else {
	nd = cm->ndidx[v];
	if      (cm->sttype[v]  == IL_st   || cm->sttype[v]  == IR_st)   ali->len += 1;
	else if (cm->ndtype[nd] == MATL_nd || cm->ndtype[nd] == MATR_nd) ali->len += 1;
	else if (cm->ndtype[nd] == MATP_nd)                              ali->len += 2;
      }	
    }

  /* Allocate and initialize.
   * Blank the annotation lines (memset calls) - only needed
   * because of the way we deal w/ EL.                                            
   */
  if (cm->annote != NULL ) ali->annote = MallocOrDie(sizeof(char) * (ali->len+1));
  else                     ali->annote = NULL;
  ali->cstr     = MallocOrDie(sizeof(char) * (ali->len+1));
  ali->cseq     = MallocOrDie(sizeof(char) * (ali->len+1));
  ali->mid      = MallocOrDie(sizeof(char) * (ali->len+1));
  ali->aseq     = MallocOrDie(sizeof(char) * (ali->len+1));

  if (cm->annote != NULL) memset(ali->annote, ' ', ali->len);
  memset(ali->cstr, ' ', ali->len);
  memset(ali->cseq, ' ', ali->len);
  memset(ali->mid,  ' ', ali->len);
  memset(ali->aseq, ' ', ali->len);

  /* Fill in the lines: traverse the traceback.
   */
  pos = 0;
  pda = CreateNstack();
  PushNstack(pda, 0);
  PushNstack(pda, PDA_STATE);
  while (PopNstack(pda, &type))
    {
      if (type == PDA_RESIDUE) {
	if (cm->annote != NULL) { PopNstack(pda, &rannote); ali->annote[pos] = rannote;}
	PopNstack(pda, &rstr); 	ali->cstr[pos] = rstr;
	PopNstack(pda, &rcons);	ali->cseq[pos] = rcons;
	PopNstack(pda, &rmid);	ali->mid[pos]  = rmid;
	PopNstack(pda, &rseq);	ali->aseq[pos] = rseq;
	pos++;
	continue;
      }
	
      /* Else, we're PDA_STATE - e.g. dealing with a trace node.
       */
      PopNstack(pda, &ti);
      v = tr->state[ti];

      /* Deal with EL (local ends, state M) as a special case.
       * We get away with only writing into aseq because we've
       * memset() the display strings to blank.
       */
      if (v == cm->M) { 
	ninset = tr->emitr[ti] - tr->emitl[ti] + 1;
	sprintf(ali->aseq+pos, "*[%d]*", ninset);
	pos += 4;
	do { pos++; ninset /= 10; } while (ninset);
	continue;
      }

      /* Fetch some info into tmp variables, for "clarity"
       */
      nd = cm->ndidx[v];	  /* what CM node we're in */
      lc   = cons->lpos[nd];	  /* where this CM node aligns to in the consensus */
      rc   = cons->rpos[nd];
      symi = dsq[tr->emitl[ti]];  /* residue indices this trace node is aligned to */
      symj = dsq[tr->emitr[ti]];

      /* Calculate four of the five lines: annote, str, cons, and seq.
       */
      do_left = do_right = FALSE;
      if (cm->sttype[v] == IL_st) {
	do_left = TRUE;
	if (cm->annote != NULL) lannote = ' ';
	lstr    = ' ';
	lcons   = '.';
	lseq    = tolower(Alphabet[symi]);
      } else if (cm->sttype[v] == IR_st) {
	do_right = TRUE;
	if (cm->annote != NULL) lannote = cm->annote[rc];
	rstr    = ' ';
	rcons   = '.';
	rseq    = tolower(Alphabet[symj]);
      } else {
	if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd) {
	  do_left = TRUE;
	  if (cm->annote != NULL) lannote = cm->annote[lc];
	  lstr  = cons->cstr[lc];
	  lcons = cons->cseq[lc];
	  if (cm->sttype[v] == MP_st || cm->sttype[v] == ML_st) lseq  = Alphabet[symi];
	  else                                                  lseq  = '-';
	}
	if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATR_nd) {
	  do_right = TRUE;
	  if (cm->annote != NULL) rannote = cm->annote[rc];
	  rstr  = cons->cstr[rc];
	  rcons = cons->cseq[rc];
	  if (cm->sttype[v] == MP_st || cm->sttype[v] == MR_st) rseq  = Alphabet[symj];
	  else                                                  rseq  = '-';
	}
      }

      /* Use emission p and score to set lmid, rmid line for emitting states.
       */
      lmid = rmid = ' ';
      if (cm->sttype[v] == MP_st) {
	if (lseq == toupper(lcons) && rseq == toupper(rcons))
	  {
	    lmid = lseq;
	    rmid = rseq;
	  }
	else if (IsCompensatory(cm->e[v], symi, symj)) 
	  lmid = rmid = ':';
	else if (DegeneratePairScore(cm->esc[v], symi, symj) >= 0)
	  lmid = rmid = '+';
      } else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) {
	if (lseq == toupper(lcons)) 
	  lmid = lseq;
	else if (DegenerateSingletScore(cm->esc[v], symi) > 0)
	  lmid = '+';
      } else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) {
	if (rseq == toupper(rcons)) 
	  rmid = rseq;
	else if (DegenerateSingletScore(cm->esc[v], symj) > 0)
	  rmid = '+';
      }

      /* If we're storing a residue leftwise - just do it.
       * If rightwise - push it onto stack.
       */
      if (do_left) {
	if (cm->annote != NULL) ali->annote[pos] = lannote;
	ali->cstr[pos] = lstr;
	ali->cseq[pos] = lcons;
	ali->mid[pos]  = lmid;
	ali->aseq[pos] = lseq;
	pos++;
      }
      if (do_right) {
	PushNstack(pda, (int) rseq);
	PushNstack(pda, (int) rmid);
	PushNstack(pda, (int) rcons);
	PushNstack(pda, (int) rstr);
	if (cm->annote != NULL) PushNstack(pda, (int) rannote);
	PushNstack(pda, PDA_RESIDUE);
      }

      /* Push the child trace nodes onto the PDA;
       * right first, so it pops last.
       */
      if (tr->nxtr[ti] != -1) {
	PushNstack(pda, tr->nxtr[ti]);
	PushNstack(pda, PDA_STATE);
      }
      if (tr->nxtl[ti] != -1) {
	PushNstack(pda, tr->nxtl[ti]);
	PushNstack(pda, PDA_STATE);
      }
    } /* end loop over the PDA; PDA now empty */
	 
  if (cm->annote != NULL) ali->annote[ali->len] = '\0';
  ali->cstr[ali->len] = '\0';
  ali->cseq[ali->len] = '\0';
  ali->mid[ali->len]  = '\0';
  ali->aseq[ali->len] = '\0';

  FreeNstack(pda);
  return ali;
}

void
PrintFancyAli(FILE *fp, Fancyali_t *ali)
{
  char *buf;
  int   pos;
  int   linelength;

  linelength = 60;
  buf = MallocOrDie(sizeof(char) * (linelength + 1));
  buf[linelength] = '\0';

  for (pos = 0; pos < ali->len; pos += linelength)
    {
      if (ali->annote != NULL) {
	strncpy(buf, ali->annote+pos, linelength);
	fprintf(fp, "%s\n", buf);
      }
      strncpy(buf, ali->cstr+pos, linelength);  fprintf(fp, "%s\n", buf);
      strncpy(buf, ali->cseq+pos, linelength);  fprintf(fp, "%s\n", buf);
      strncpy(buf, ali->mid+pos,  linelength);  fprintf(fp, "%s\n", buf);
      strncpy(buf, ali->aseq+pos, linelength);  fprintf(fp, "%s\n", buf);
      fprintf(fp, "\n");
    }
  free(buf);
}


/* Function:  FreeFancyAli()
 * Incept:    SRE, Fri May 24 15:37:30 2002 [St. Louis]
 */
void
FreeFancyAli(Fancyali_t *ali)
{
  if (ali->annote != NULL) free(ali->annote);
  if (ali->cstr   != NULL) free(ali->cstr);
  if (ali->cseq   != NULL) free(ali->cseq);
  if (ali->mid    != NULL) free(ali->mid);
  if (ali->aseq   != NULL) free(ali->aseq);
  free(ali);
}


/* Function:  CreateCMConsensus()
 * Incept:    SRE, Thu May 23 10:39:39 2002 [St. Louis]
 *
 * Purpose:   Create displayable strings for consensus sequence
 *            and consensus structure; also create map of 
 *            nodes -> right and left consensus positions.
 *            
 *            Consensus sequence shows maximum scoring residue(s)
 *            for each emitting node. If score < pthresh (for pairs)
 *            or < sthresh (for singlets), the residue is shown
 *            in lower case. (That is, "strong" consensus residues
 *            are in upper case.)
 *            
 *            Consensus structure annotates
 *            base pairs according to "multifurcation order" (how
 *            many multifurcation loops are beneath this pair).
 *               terminal stems:  <>
 *               order 1:         ()
 *               order 2:         []
 *               order >2:        {}
 *            Singlets are annotated . if external, _ if hairpin,
 *            - if bulge or interior loop, and , for multifurcation loop.
 *               
 *            Example:
 *                ..(((,,<<<__>>>,<<<__>>->,,)))..
 *                AAGGGAACCCTTGGGTGGGTTCCACAACCCAA   
 *
 * Args:      cm         - the model
 *            pthresh    - bit score threshold for base pairs to be lowercased
 *            sthresh    - bit score threshold for singlets to be lowercased
 *            
 * Returns:   CMConsensus_t structure.
 *            Caller frees w/ FreeCMConsensus().
 *
 * Xref:      STL6 p.58.
 */
CMConsensus_t *
CreateCMConsensus(CM_t *cm, float pthresh, float sthresh)
{
  CMConsensus_t *con;           /* growing consensus info */
  char     *cseq;               /* growing consensus sequence display string   */
  char     *cstr;               /* growing consensus structure display string  */
  int      *lpos, *rpos;        /* maps node->consensus position, [0..nodes-1] */
  int       cpos;		/* current position in cseq, cstr              */
  int       nalloc;		/* current allocated length of cseq, cstr      */
  Nstack_t *pda;                /* pushdown automaton used to traverse model   */
  int      *multiorder;         /* "height" of each node (multifurcation order), [0..nodes-1] */
  int      *inface;             /* face count for each node, inside */
  int      *outface;            /* face count for each node, outside */
  int       v, nd;
  int       type;
  char      lchar, rchar;
  char      lstruc, rstruc;
  int       x;

  lpos    = MallocOrDie(sizeof(int) * cm->nodes);
  rpos    = MallocOrDie(sizeof(int) * cm->nodes);
  cseq    = MallocOrDie(sizeof(char) * 100);
  cstr    = MallocOrDie(sizeof(char) * 100);
  nalloc  = 100;
  cpos    = 0;

  for (nd = 0; nd < cm->nodes; nd++) 
    lpos[nd] = rpos[nd] = -1;

  multiorder = createMultifurcationOrderChart(cm);
  createFaceCharts(cm, &inface, &outface);

  pda = CreateNstack();
  PushNstack(pda, 0);
  PushNstack(pda, PDA_STATE);
  while (PopNstack(pda, &type)) {
    if (type == PDA_RESIDUE) 
      {
	PopNstack(pda, &x); rchar  = (char) x;
	PopNstack(pda, &x); rstruc = (char) x;
	PopNstack(pda, &nd);
	rpos[nd]   = cpos;
	cseq[cpos] = rchar;
	cstr[cpos] = rstruc;
	cpos++;
      }
    else if (type == PDA_MARKER) 
      {
	PopNstack(pda, &nd);
	rpos[nd]   = cpos;
      }
    else if (type == PDA_STATE) 
      {
	PopNstack(pda, &v);
	nd    = cm->ndidx[v];
	lchar = rchar = lstruc = rstruc = 0;

	/* Determine what we emit: 
	 * MATP, MATL, MATR consensus states only.
	 */
	if (cm->stid[v] == MATP_MP) 
	  {
	    x = FMax(cm->esc[v], Alphabet_size*Alphabet_size);
	    lchar = Alphabet[x / Alphabet_size];
	    rchar = Alphabet[x % Alphabet_size];
	    if (cm->esc[v][x] < pthresh) {
	      lchar = tolower(lchar);
	      rchar = tolower(rchar);
	    }
	    switch (multiorder[nd]) {
	    case 0:  lstruc = '<'; rstruc = '>'; break;
	    case 1:  lstruc = '('; rstruc = ')'; break;
	    case 2:  lstruc = '['; rstruc = ']'; break;
	    default: lstruc = '{'; rstruc = '}'; break;
	    }
	} else if (cm->stid[v] == MATL_ML) {
	  x = FMax(cm->esc[v], Alphabet_size);
	  lchar = Alphabet[x];
	  if (cm->esc[v][x] < sthresh) lchar = tolower(lchar);
	  if      (outface[nd] == 0)                    lstruc = '.'; /* external ss */
	  else if (inface[nd] == 0 && outface[nd] == 1) lstruc = '_'; /* hairpin loop */
	  else if (inface[nd] == 1 && outface[nd] == 1) lstruc = '-'; /* bulge/interior */
	  else                                          lstruc = ','; /* multiloop */
	  rstruc = ' ';
	} else if (cm->stid[v] == MATR_MR) {
	  x = FMax(cm->esc[v], Alphabet_size);
	  rchar = Alphabet[x];
	  if (cm->esc[v][x] < sthresh) rchar = tolower(rchar);
	  if      (outface[nd] == 0)                    rstruc = '.'; /* external ss */
	  else if (inface[nd] == 0 && outface[nd] == 1) rstruc = '?'; /* doesn't happen */
	  else if (inface[nd] == 1 && outface[nd] == 1) rstruc = '-'; /* bulge/interior */
	  else                                          rstruc = ','; /* multiloop */
	  lstruc = ' ';
	}

	/* Emit. A left base, we can do now; 
	 * a right base, we defer onto PDA.
	 */
	lpos[nd]   = cpos;	/* we always set lpos, rpos even for nonemitters */
	if (lchar) {
	  cseq[cpos] = lchar;
	  cstr[cpos] = lstruc;
	  cpos++;
	}
	if (rchar) {
	  PushNstack(pda, nd);
	  PushNstack(pda, rstruc);
	  PushNstack(pda, rchar);
	  PushNstack(pda, PDA_RESIDUE);
	} else {
	  PushNstack(pda, nd);
	  PushNstack(pda, PDA_MARKER);
	}

	/* Transit - to consensus states only.
	 * The obfuscated idiom finds the index of the next consensus
	 * state without making assumptions about numbering or connectivity.
	 */
	if (cm->sttype[v] == B_st) {
	  PushNstack(pda, cm->cnum[v]);     /* right S  */
	  PushNstack(pda, PDA_STATE);
	  PushNstack(pda, cm->cfirst[v]);   /* left S */
	  PushNstack(pda, PDA_STATE);
	} else if (cm->sttype[v] != E_st) {
	  v = cm->nodemap[cm->ndidx[cm->cfirst[v] + cm->cnum[v] - 1]];
	  PushNstack(pda, v);
	  PushNstack(pda, PDA_STATE);
	}
      } /*end PDA_STATE block*/

    if (cpos == nalloc) {
      nalloc += 100;
      cseq = ReallocOrDie(cseq, sizeof(char) * nalloc);
      cstr = ReallocOrDie(cstr, sizeof(char) * nalloc);
    }
  }/* PDA now empty... done generating cseq, cstr, and node->consensus residue map */
  cseq[cpos] = '\0';
  cstr[cpos] = '\0';

  FreeNstack(pda);
  free(multiorder);
  free(inface);
  free(outface);

  con = MallocOrDie(sizeof(CMConsensus_t));
  con->cseq = cseq;
  con->cstr = cstr;
  con->lpos = lpos;
  con->rpos = rpos;
  con->clen = cpos;
  return con;
}

void
FreeCMConsensus(CMConsensus_t *con)
{
  if (con->cseq != NULL) free(con->cseq);
  if (con->cstr != NULL) free(con->cstr);
  if (con->lpos != NULL) free(con->lpos);
  if (con->rpos != NULL) free(con->rpos);
  free(con);
}

/* Function:  createMultifurcationOrderChart()
 * Incept:    SRE, Thu May 23 09:48:33 2002 [St. Louis] 
 *
 * Purpose:   Calculates the degree of multifurcation beneath
 *            the master subtree rooted at every node n.
 *            Returns [0..nodes-1] array of these values.
 *
 *            Terminal stems have value 0. All nodes n starting with
 *            the BEG node for a terminal stem have height[n] = 0.
 *            
 *            A stem "above" a multifurcation into all terminal stems
 *            has value 1; all nodes n starting with BEG and ending
 *            with BIF have height[n] = 1.
 * 
 *            And so on, for "higher order" (deeper) multifurcations.
 * 
 *            Used for figuring out what characters we'll display a
 *            consensus pair with.
 *            
 *            THIS FUNCTION IS BUGGY (Sat Jun  1 12:24:23 2002)
 *            
 * Args:      cm - the model.
 *
 * Returns:   [0..cm->nodes-1] array of multifurcation orders, for each node.
 *            This array is allocated here; caller free's w/ free().
 *
 * xref:     STL6 p.58.
 */
static int *
createMultifurcationOrderChart(CM_t *cm)
{
  int  v, nd, left, right;
  int *height;
  int *seg_has_pairs;

  height        = MallocOrDie(sizeof(int) * cm->nodes);
  seg_has_pairs = MallocOrDie(sizeof(int) * cm->nodes);
  for (nd = cm->nodes-1; nd >= 0; nd--)
    {
      v = cm->nodemap[nd];

      if       (cm->stid[v] == MATP_MP) seg_has_pairs[nd] = TRUE;
      else if  (cm->stid[v] == END_E)   seg_has_pairs[nd] = FALSE;
      else if  (cm->stid[v] == BIF_B)   seg_has_pairs[nd] = FALSE;
      else                              seg_has_pairs[nd] = seg_has_pairs[nd+1];

      if (cm->stid[v] == END_E) 
	height[nd]        = 0;
      else if (cm->stid[v] == BIF_B) 
	{
	  left  = cm->ndidx[cm->cfirst[v]]; 
	  right = cm->ndidx[cm->cnum[v]];
	  height[nd] = MAX(height[left] + seg_has_pairs[left],
			   height[right] + seg_has_pairs[right]);
	}
      else
	height[nd] = height[nd+1]; 
    }
  free(seg_has_pairs);
  return height;
}	
	
     
/* Function:  createFaceCharts()
 * Incept:    SRE, Thu May 23 12:40:04 2002 [St. Louis]
 *
 * Purpose:   Calculate "inface" and "outface" for each node
 *            in the consensus (master) structure of the CM.
 *            These can be used to label nodes:
 *                                inface       outface    
 *                             ------------   ----------
 *             external ss         any           0                   
 *             hairpin loop         0            1
 *             bulge/interior       1            1
 *             multifurc           >1            1  
 *             multifurc            1           >1   
 *             doesn't happen       0           >1
 *             
 *             hairpin closing bp   0            1
 *             extern closing bp    1            0
 *             stem bp              1            1
 *             multifurc close bp  >1            1
 *             multifurc close bp   1           >1
 *             doesn't happen       0           >1
 *
 * Args:       cm          - the model
 *             ret_inface  - RETURN: inface[0..nodes-1]
 *             ret_outface - RETURN: outface[0..nodes-1]         
 *
 * Returns:    inface, outface; 
 *             they're alloc'ed here. Caller free's with free()
 *
 * Xref:       STL6 p.58
 */
static void
createFaceCharts(CM_t *cm, int **ret_inface, int **ret_outface)
{
  int *inface;
  int *outface;
  int  nd, left, right, parent;
  int  v,w,y;

  inface  = MallocOrDie(sizeof(int) * cm->nodes);
  outface = MallocOrDie(sizeof(int) * cm->nodes);

  /* inface - the number of faces below us in descendant
   *          subtrees. if 0, we're either external, or
   *          a closing basepair, or we're in a hairpin loop. 
   *          inface is exclusive of current pair - so we
   *          can easily detect closing base pairs.
   */
  for (nd = cm->nodes-1; nd >= 0; nd--)
    {
      v = cm->nodemap[nd];
      if      (cm->ndtype[nd] == END_nd) inface[nd] = 0;
      else if (cm->ndtype[nd] == BIF_nd) {
	left  = cm->ndidx[cm->cfirst[v]];
	right = cm->ndidx[cm->cnum[v]];
	inface[nd] = inface[left] + inface[right];
      } else {
	if (cm->ndtype[nd+1] == MATP_nd) inface[nd] = 1;
	else                             inface[nd] = inface[nd+1];
      }
    }

  /* outface - the number of faces above us in the tree
   *           excluding our subtree. if 0, we're external.
   *           Like inface, outface is exclusive of current
   *           pair.
   */
  for (nd = 0; nd < cm->nodes; nd++)
    {
      v = cm->nodemap[nd];
      if      (cm->ndtype[nd] == ROOT_nd) outface[nd] = 0;
      else if (cm->ndtype[nd] == BEGL_nd) 
	{
	  parent = cm->ndidx[cm->plast[v]];
	  y      = cm->nodemap[parent];
	  right  = cm->ndidx[cm->cnum[y]];
	  outface[nd] = outface[parent] + inface[right];
	}
      else if (cm->ndtype[nd] == BEGR_nd)
	{
	  parent = cm->ndidx[cm->plast[v]];
	  w      = cm->nodemap[parent];
	  left   = cm->ndidx[cm->cfirst[y]];
	  outface[nd] = outface[parent] + inface[left];
	}
      else 
	{
	  parent = nd-1;
	  if (cm->ndtype[parent] == MATP_nd) outface[nd] = 1;
	  else                               outface[nd] = outface[parent];
	}
    }
  
  *ret_inface  = inface;
  *ret_outface = outface;
}
	
	    

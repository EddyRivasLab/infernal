/* Formatting, transmitting, and printing single alignments to a
 * profile.
 * 
 * Contents:
 *   1. The CM_ALIDISPLAY object.
 *   2. The CM_ALIDISPLAY API.
 *   3. Debugging/dev code.
 *   4. Benchmark driver.
 *   5. Unit tests.
 *   6. Test driver.
 *   7. Example.
 *   8. Copyright and license.
 */
#include "esl_config.h"
#include "p7_config.h"
#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "easel.h"
#include "esl_stack.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "funcs.h"
#include "structs.h"

static int  bp_is_canonical(char lseq, char rseq);

/*****************************************************************
 * 1. The CM_ALIDISPLAY object
 *****************************************************************/


/* Function:  cm_alidisplay_Create()
 * Incept:    EPN, Wed May 25 05:38:12 2011
 *            SRE, Thu May 23 13:46:09 2002 [St. Louis] (display.c:CreateFancyAli())
 *
 * Synopsis:  Create an alignment display, from parsetree and model.
 *
 * Purpose:   Creates and returns an alignment display <which> in
 *            traceback <tr>, where the traceback corresponds to an
 *            alignment of optimized profile <om> to digital sequence
 *            <dsq>, and the unique name of that target sequence <dsq>
 *            is <sqname>. The <which> index starts at 0.
 *
 * Args:      tr     - parsetree
 *            cm     - covariance model (query)
 *            sq     - digital sequence (target)
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <NULL> on allocation failure, or if something's internally corrupt 
 *            in the data.
 * Purpose:   Given a trace (and the model and sequence it corresponds
 *            to), create a pairwise alignment for display; return in a Fancyali_t
 *            structure.
 *
 * Args:      abc   - alphabet to create alignment with (often cm->abc)
 *            tr    - parsetree for cm aligned to dsq
 *            cm    - model
 *            cons  - consensus information for cm; see CreateCMConsensus()
 *            sq    - the sequence, parsetree corresponds to subsequence beginning at spos
 *            spos  - position in sq which corresponds to first position in tr
 *            ppstr - posterior probability string 
 *            aln_sc       - if(used_optacc) avg post prob of all aligned residues, else CYK score
 *            used_optacc  - TRUE if aln algorithm used was optimal accuracy 
 *            used_hbands  - TRUE if HMM bands were used for alignment
 *            matrix_Mb    - size of DP matrix in Mb used for alignment
 *            elapsed_secs - time (seconds) required for alignment
 *
 * Returns:   CM_ALIDISPLAY structure.
 *            Caller frees, with cm_alidisplay_Destroy().
 *
 * Xref:      STL6 p.58
 */
CM_ALIDISPLAY *
cm_alidisplay_Create(const ESL_ALPHABET *abc, Parsetree_t *tr, CM_t *cm, CMConsensus_t *cons, const ESL_SQ *sq, int64_t seqoffset, char *ppstr, float aln_sc, int used_optacc, int used_hbands, float matrix_Mb, double elapsed_secs)
{
  int            status;
  CM_ALIDISPLAY *ad = NULL;      /* alidisplay structure we're building       */
  ESL_STACK  *pda = NULL;        /* pushdown automaton used to traverse trace */
  int         type;		 /* type of pda move: PDA_RESIDUE, PDA_STATE  */
  int         v;		 /* state index       */
  int         nd;		 /* node index        */
  int         ti;		 /* position in trace */
  int         qinset, tinset;	 /* # consensus nt skipped by an EL, or in an EL */
  int         ninset;		 /* max # nt in an EL     */
  int         pos;		 /* position in growing ali */
  int         lc, rc;		 /* indices for left, right pos in consensus */
  int         symi, symj;
  int         d;
  char        mode;
  int         lrf, rrf;  /* chars in annotation line; left, right     */
  int         lstr, rstr;	 /* chars in structure line; left, right      */
  int         lcons, rcons;	 /* chars in consensus line; left, right      */
  int         lmid, rmid;	 /* chars in ali quality line; left, right    */
  int         lnnc, rnnc;	 /* chars in negative scoring noncanonical line; left right */
  int         lseq, rseq;	 /* chars in aligned target line; left, right */
  int         lpost, rpost;	 /* chars in aligned posteriors, left, right  */
  int         do_left, do_right; /* flags to generate left, right             */
  int         cpos_l, cpos_r;    /* positions in consensus (1..clen)          */
  int         spos_l, spos_r;	 /* positions in dsq (1..L)                   */
  int        *scoord = NULL;     /* [0..N-1] coordinates for aligned residues */
  int        *ccoord = NULL;     /* [0..N-1] coordinates for model positions  */
  int         cm_namelen, cm_acclen, cm_desclen;
  int         sq_namelen, sq_acclen, sq_desclen;
  int         len, n;

  /* Contract check. We allow the caller to specify the alphabet they want the 
   * resulting MSA in, but it has to make sense (see next few lines). */
  if(cm->abc->type == eslRNA) { 
    if(abc->type != eslRNA && abc->type != eslDNA)
      cm_Fail("ERROR in cm_alidisplay_Create(), cm alphabet is RNA, but requested output alphabet is neither DNA nor RNA.");
  }
  else if(cm->abc->K != abc->K) {
    cm_Fail("ERROR in cm_alidisplay_Create(), cm alphabet size is %d, but requested output alphabet size is %d.", cm->abc->K, abc->K);
  }
  
  /* Calculate length of the alignment display (len).
   *   MATP node        : +2
   *   MATL, MATR node  : +1
   *   IL, IR state     : +1
   *   EL:              : 4 + width of length display : "*[nn]*"
   *   anything else    : 0.
   */
  len = 0;
  for (ti = 0; ti < tr->n; ti++) { 
    v  = tr->state[ti];
    if (v == cm->M) {  /* special case: local exit into EL */
      nd = cm->ndidx[tr->state[ti-1]]; /* calculate node that EL replaced */
      qinset     = cons->rpos[nd] - cons->lpos[nd] + 1;
      tinset     = tr->emitr[ti]  - tr->emitl[ti]  + 1;
      ninset     = ESL_MAX(qinset,tinset);
      len += 4;
      do { len++; ninset/=10; } while (ninset); /* poor man's (int)log_10(ninset)+1 */
      continue;
    } else {
      nd = cm->ndidx[v];
      if      (cm->sttype[v]  == IL_st   || cm->sttype[v]  == IR_st)   len += 1;
      else if (cm->ndtype[nd] == MATL_nd || cm->ndtype[nd] == MATR_nd) len += 1;
      else if (cm->ndtype[nd] == MATP_nd)                              len += 2;
    }	
    /* Catch marginal-type local ends and treat them like EL for output */
    if ((tr->nxtl[ti] == -1) && (cm->sttype[v] != E_st)) {
      nd = cm->ndidx[tr->state[ti]];
      qinset     = cons->rpos[nd] - cons->lpos[nd] + 1;
      tinset     = tr->emitr[ti]  - tr->emitl[ti]  + 1;
      if (tinset > 0) tinset--;
      ninset     = ESL_MAX(qinset,tinset);
      len += 4;
      do { len++; ninset/=10; } while (ninset); /* poor man's (int)log_10(ninset)+1 */
    }
  }

  /* Now we know the length of all arrays (len), determine total amount of memory required, and allocate it */
  /* Allocate the char arrays */

  n = (len+1) * 5; /* model, csline, mline, aseq, nline mandatory */
  if(ppstr  != NULL) n += len+1;
  if(cm->rf != NULL) n += len+1;
  cm_namelen = strlen(cm->name);                           n += cm_namelen + 1;
  cm_acclen  = (cm->acc  != NULL ? strlen(cm->acc)  : 0);  n += cm_acclen  + 1; 
  cm_desclen = (cm->desc != NULL ? strlen(cm->desc) : 0);  n += cm_desclen + 1; 
  sq_namelen = strlen(sq->name);                           n += sq_namelen + 1;
  sq_acclen  = strlen(sq->acc);                            n += sq_acclen  + 1; /* sq->acc is "\0" when unset */
  sq_desclen = strlen(sq->desc);                           n += sq_desclen + 1; /* sq->desc is "\0" when unset */

  ESL_ALLOC(ad, sizeof(CM_ALIDISPLAY));
  ad->mem          = NULL;
  ad->memsize      = sizeof(char) * n;
  ad->aln_sc       = aln_sc;
  ad->used_optacc  = used_optacc;
  ad->used_hbands  = used_hbands;
  ad->matrix_Mb    = matrix_Mb;
  ad->elapsed_secs = elapsed_secs;
  ad->N            = len;

  pos = 0;
  ESL_ALLOC(ad->mem, ad->memsize);
  if (cm->rf != NULL) { ad->rfline = ad->mem + pos;  pos += len+1;} else { ad->rfline = NULL; }
  ad->nline   = ad->mem + pos;  pos += len+1;
  ad->csline  = ad->mem + pos;  pos += len+1;
  ad->model   = ad->mem + pos;  pos += len+1;
  ad->mline   = ad->mem + pos;  pos += len+1;
  ad->aseq    = ad->mem + pos;  pos += len+1;
  if (ppstr  != NULL) { ad->ppline = ad->mem + pos;  pos += len+1;} else { ad->ppline = NULL; }
  ad->cmname  = ad->mem + pos;  pos += cm_namelen +1;
  ad->cmacc   = ad->mem + pos;  pos += cm_acclen +1;
  ad->cmdesc  = ad->mem + pos;  pos += cm_desclen +1;
  ad->sqname  = ad->mem + pos;  pos += sq_namelen +1;
  ad->sqacc   = ad->mem + pos;  pos += sq_acclen +1;  
  ad->sqdesc  = ad->mem + pos;  pos += sq_desclen +1; 

  /* Set name, acc, desc char arrays */
  strcpy(ad->cmname, cm->name);
  if (cm->acc  != NULL) strcpy(ad->cmacc,  cm->acc);  else ad->cmacc[0]  = 0;
  if (cm->desc != NULL) strcpy(ad->cmdesc, cm->desc); else ad->cmdesc[0] = 0;
  strcpy(ad->sqname,  sq->name);
  strcpy(ad->sqacc,   sq->acc);
  strcpy(ad->sqdesc,  sq->desc);

  ad->clen = cm->clen;
  ad->L    = sq->n;

  /* Allocate and initialize.
   * Blank the annotation lines (memset calls) - only needed
   * because of the way we deal w/ EL. 
   */
  if (cm->rf != NULL) memset(ad->rfline, ' ', ad->N);
  memset(ad->nline,   ' ', ad->N);
  memset(ad->csline,  ' ', ad->N);
  memset(ad->model,   ' ', ad->N);
  memset(ad->mline,   ' ', ad->N);
  memset(ad->aseq,    ' ', ad->N);
  if(ppstr != NULL)   memset(ad->ppline, ' ', ad->N);

  ESL_ALLOC(ccoord, sizeof(int) * len);
  ESL_ALLOC(scoord, sizeof(int) * len);
  esl_vec_ISet(ccoord, ad->N, 0);
  esl_vec_ISet(scoord, ad->N, 0);

  /* Fill in the lines: traverse the traceback.
   */
  pos = 0;
  if((pda = esl_stack_ICreate()) == NULL) goto ERROR;
  if((status = esl_stack_IPush(pda, 0)) != eslOK) goto ERROR;
  if((status = esl_stack_IPush(pda, PDA_STATE)) != eslOK) goto ERROR;

  while (esl_stack_IPop(pda, &type) != eslEOD)
    {
      if (type == PDA_RESIDUE) {
	if (cm->rf != NULL) { 
	  esl_stack_IPop(pda, &rrf); 
	  ad->rfline[pos] = rrf;
	}
	esl_stack_IPop(pda, &rnnc);	  ad->nline[pos]  = rnnc;
	esl_stack_IPop(pda, &rstr); 	  ad->csline[pos] = rstr;
	esl_stack_IPop(pda, &rcons);	  ad->model[pos]  = rcons;
	esl_stack_IPop(pda, &rmid);	  ad->mline[pos]  = rmid;
	esl_stack_IPop(pda, &rseq);       ad->aseq[pos]   = rseq;
	if(ppstr != NULL) {
	  esl_stack_IPop(pda, &rpost);    ad->ppline[pos] = rpost;
	}
	esl_stack_IPop(pda, &cpos_r);     ccoord[pos] = cpos_r;
	esl_stack_IPop(pda, &spos_r);     scoord[pos] = spos_r;
	pos++;
	continue;
      }
	
      /* Else, we're PDA_STATE - e.g. dealing with a trace node.
       */
      esl_stack_IPop(pda, &ti);
      v = tr->state[ti];

      /* Deal with EL (local ends, state M) as a special case.
       * We get away with only writing into aseq because we've
       * memset() the display strings to blank.
       */
      if (v == cm->M) { 
	int numwidth;		/* number of chars to leave for displaying width numbers */

	nd = 1 + cm->ndidx[tr->state[ti-1]]; /* calculate node that EL replaced */
	qinset     = cons->rpos[nd] - cons->lpos[nd] + 1;
	tinset     = tr->emitr[ti]  - tr->emitl[ti]  + 1;
	ninset     = ESL_MAX(qinset,tinset);
	numwidth = 0; do { numwidth++; ninset/=10; } while (ninset); /* poor man's (int)log_10(ninset)+1 */
	memset(ad->csline+pos,  '~', numwidth+4);
	sprintf(ad->model+pos, "*[%*d]*", numwidth, qinset);
	sprintf(ad->aseq+pos, "*[%*d]*", numwidth, tinset);
	/* do nothing for posteriors here, they'll stay as they were init'ed, as ' ' */
	pos += 4 + numwidth;
	continue;
      }

      /* Fetch some info into tmp variables, for "clarity"
       */
      nd   = cm->ndidx[v];	  /* what CM node we're in */
      lc   = cons->lpos[nd];	  /* where CM node aligns to in consensus (left) */
      rc   = cons->rpos[nd];      /* where CM node aligns to in consensus (right) */
      symi = sq->dsq[tr->emitl[ti] + (seqoffset-1)];  /* residue indices that node is aligned to (left) */
      symj = sq->dsq[tr->emitr[ti] + (seqoffset-1)];  /* residue indices that node is aligned to (right) */
      if(ppstr != NULL) { /* posterior codes are indexed 0..alen-1, off-by-one w.r.t dsq */
	lpost = '.'; /* init to gap, if it corresponds to a residue, we'll reset it below */
	rpost = '.'; /* init to gap, if it corresponds to a residue, we'll reset it below */
      }
      d = tr->emitr[ti] - tr->emitl[ti] + 1;
      mode = tr->mode[ti];

      /* Calculate four of the six lines: rfline, csline, model, and aseq.
       */
      do_left = do_right = FALSE;
      if (cm->sttype[v] == IL_st) {
	do_left = TRUE;
	if (cm->rf != NULL) lrf = '.';
	lstr    = '.';
	lcons   = '.';
	lseq = (mode == TRMODE_J || mode == TRMODE_L) ? tolower((int) abc->sym[symi]) : '~';
	cpos_l  = 0;
	spos_l  = tr->emitl[ti] + seqoffset-1;
	if(ppstr != NULL) { 
	  lpost = (mode == TRMODE_J || mode == TRMODE_L) ? ppstr[tr->emitl[ti]-1] : '~'; /* watch off-by-one b/t ppstr and dsq */
	}
      } else if (cm->sttype[v] == IR_st) {
	do_right = TRUE;
	if (cm->rf != NULL) rrf = '.';
	rstr    = '.';
	rcons   = '.';
	rseq = (mode == TRMODE_J || mode == TRMODE_R) ? tolower((int) abc->sym[symj]) : '~';
	cpos_r  = 0;
	spos_r  = tr->emitr[ti] + seqoffset-1;
	if(ppstr != NULL) { 
	  rpost = (mode == TRMODE_J || mode == TRMODE_R) ? ppstr[tr->emitr[ti]-1] : '~'; /* watch off-by-one b/t ppstr and dsq */
	}
      } else {
	if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd) {
	  do_left = TRUE;
	  if (cm->rf != NULL) lrf = cm->rf[lc+1];
	  lstr   = cons->cstr[lc];
	  lcons  = (cm->flags & CMH_CONS) ? cm->consensus[(lc+1)] : cons->cseq[lc];
	  cpos_l = lc+1;
	  if (cm->sttype[v] == MP_st || cm->sttype[v] == ML_st) {
	    lseq = (mode == TRMODE_J || mode == TRMODE_L) ? abc->sym[symi] : '~';
	    spos_l = tr->emitl[ti] + seqoffset-1;
	    if(ppstr != NULL) { 
	      lpost = (mode == TRMODE_J || mode == TRMODE_L) ? ppstr[tr->emitl[ti]-1] : '~'; /* watch off-by-one b/t ppstr and dsq */
	    }
	  } else {
	    lseq   = (mode == TRMODE_J || mode == TRMODE_L) ? '-' : '~';
	    spos_l = 0;
	    if(ppstr != NULL) { 
	      lpost  = (mode == TRMODE_J || mode == TRMODE_L) ? '.' : '~';
	    }
	  }
	}
	if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATR_nd) {
	  do_right = TRUE;
	  if (cm->rf != NULL) rrf = cm->rf[rc+1];
	  rstr   = cons->cstr[rc];
	  rcons  = (cm->flags & CMH_CONS) ? cm->consensus[(rc+1)] : cons->cseq[rc];
	  cpos_r = rc+1;
	  if (cm->sttype[v] == MP_st || cm->sttype[v] == MR_st) {
	    rseq = (mode == TRMODE_J || mode == TRMODE_R) ? abc->sym[symj] : '~';
	    spos_r = tr->emitr[ti] + seqoffset-1;
	    if(ppstr != NULL) { 
	      rpost = (mode == TRMODE_J || mode == TRMODE_R) ? ppstr[tr->emitr[ti]-1] : '~'; /* watch off-by-one b/t ppstr and dsq */
	    }
	  } else {
	    rseq   = (mode == TRMODE_J || mode == TRMODE_R) ? '-' : '~';
	    spos_r = 0;
	    if(ppstr != NULL) { 
	      rpost  = (mode == TRMODE_J || mode == TRMODE_R) ? '.' : '~';
	    }
	  }
	}
      }
      
      /* Use emission p and score to set lmid, rmid line for emitting states.
       */
      lmid = rmid = ' ';
      lnnc = rnnc = ' ';
      if (cm->sttype[v] == MP_st) {
	if (lseq == toupper(lcons) && rseq == toupper(rcons))
	  {
	    lmid = lseq;
	    rmid = rseq;
	  }
        else if (mode != 3)
          {
            if (mode == TRMODE_L && lseq == toupper(lcons)) lmid = lseq;
            if (mode == TRMODE_R && rseq == toupper(rcons)) rmid = rseq;
          }
	else if (DegeneratePairScore(cm->abc, cm->esc[v], symi, symj) >= 0) 
	  lmid = rmid = ':';
	
	/* determine lnnc, rnnc for optional negative scoring non-canonical annotation, they are 'v' 
	 * if lseq and rseq are a negative scoring non-canonical (not a AU,UA,GC,CG,GU,UG) pair */
	if ((mode == TRMODE_J) && (DegeneratePairScore(cm->abc, cm->esc[v], symi, symj) < 0) && (! bp_is_canonical(lseq, rseq))) {
	  lnnc = rnnc = 'v';
	}
      } else if (cm->sttype[v] == ML_st || cm->sttype[v] == IL_st) {
	if (lseq == toupper(lcons)) 
	  lmid = lseq;
        else if ( (mode != 3) && (mode != 2) )
          ;
	else if(esl_abc_FAvgScore(cm->abc, symi, cm->esc[v]) > 0)
	  lmid = '+';
      } else if (cm->sttype[v] == MR_st || cm->sttype[v] == IR_st) {
	if (rseq == toupper(rcons)) 
	  rmid = rseq;
        else if ( (mode != 3) && (mode != 1) )
          ;
	else if(esl_abc_FAvgScore(cm->abc, symj, cm->esc[v]) > 0)
	  rmid = '+';
      }
      if(cm->stid[v] == MATP_ML || cm->stid[v] == MATP_MR) { 
	if(mode == TRMODE_J) { 
	  lnnc = rnnc = 'v'; /* mark non-truncated half base-pairs (MATP_ML or MATP_MR) with 'v' */
	}
      }
      /* If we're storing a residue leftwise - just do it.
       * If rightwise - push it onto stack.
       */
      if (do_left) {
	if (cm->rf != NULL) ad->rfline[pos] = lrf;
	ad->nline[pos]   = lnnc;
	ad->csline[pos]  = lstr;
	ad->model[pos]   = lcons;
	ad->mline[pos]   = lmid;
	ad->aseq[pos]    = lseq;
	if(ppstr != NULL)  ad->ppline[pos] = lpost;
	ccoord[pos] = cpos_l;
	scoord[pos] = spos_l;
	pos++;
      }
      if (do_right) {
	if ((status = esl_stack_IPush(pda, spos_r)) != eslOK) goto ERROR;
	if ((status = esl_stack_IPush(pda, cpos_r)) != eslOK) goto ERROR;
	if(ppstr != NULL) {
	  if ((status = esl_stack_IPush(pda, (int) rpost)) != eslOK) goto ERROR;
	}
	if ((status = esl_stack_IPush(pda, (int) rseq))  != eslOK) goto ERROR;
	if ((status = esl_stack_IPush(pda, (int) rmid))  != eslOK) goto ERROR;
	if ((status = esl_stack_IPush(pda, (int) rcons)) != eslOK) goto ERROR;
	if ((status = esl_stack_IPush(pda, (int) rstr))  != eslOK) goto ERROR;
	if ((status = esl_stack_IPush(pda, (int) rnnc))  != eslOK) goto ERROR; 
	if (cm->rf != NULL) if ((status = esl_stack_IPush(pda, (int) rrf)) != eslOK) goto ERROR;
	if ((status = esl_stack_IPush(pda, PDA_RESIDUE)) != eslOK) goto ERROR;
      }
      
      /* Push the child trace nodes onto the PDA;
       * right first, so it pops last.
       */
      if (tr->nxtr[ti] != -1) {
	if ((status = esl_stack_IPush(pda, tr->nxtr[ti])) != eslOK) goto ERROR;
	if ((status = esl_stack_IPush(pda, PDA_STATE)) != eslOK) goto ERROR;
      }
      if (tr->nxtl[ti] != -1) {
	if ((status = esl_stack_IPush(pda, tr->nxtl[ti])) != eslOK) goto ERROR;
	if ((status = esl_stack_IPush(pda, PDA_STATE)) != eslOK) goto ERROR;
      }
      else if (cm->sttype[v] != E_st) {
	; /* marginal type local end, do nothing */
#if 0
        /* Catch marginal-type local ends, treat like EL for output */
	int numwidth;		/* number of chars to leave for displaying width numbers */

	nd = 1 + cm->ndidx[tr->state[ti]]; /* calculate node that EL replaced */
	qinset     = cons->rpos[nd] - cons->lpos[nd] + 1;
	tinset     = tr->emitr[ti]  - tr->emitl[ti]  + 1;
        if (tinset > 0) tinset--;
	ninset     = ESL_MAX(qinset,tinset);
	numwidth = 0; do { numwidth++; ninset/=10; } while (ninset); /* poor man's (int)log_10(ninset)+1 */
	memset(ad->csline+pos,  '~', numwidth+4);
	sprintf(ad->model+pos, "*[%*d]*", numwidth, qinset);
	sprintf(ad->aseq+pos, "*[%*d]*", numwidth, tinset);
	/* do nothing for posteriors here, they'll stay as they were init'ed, as ' ' */
	pos += 4 + numwidth;
#endif
      }
    } /* end loop over the PDA; PDA now empty */
	 
  if(cm->rf != NULL) ad->rfline[ad->N] = '\0';
  ad->nline[ad->N]  = '\0';
  ad->csline[ad->N] = '\0';
  ad->model[ad->N]  = '\0';
  ad->mline[ad->N]  = '\0';
  ad->aseq[ad->N]   = '\0';
  if(ppstr != NULL)  ad->ppline[ad->N] = '\0'; 

  /* Laboriously determine the maximum bounds. */
  ad->sqfrom = 0;
  for (pos = 0; pos < ad->N; pos++)
    if (scoord[pos] != 0) {
      ad->sqfrom = scoord[pos];
      break;
    }
  ad->sqto = 0;
  for (pos = 0; pos < ad->N; pos++)
    if (scoord[pos] != 0) ad->sqto = scoord[pos];
  ad->cfrom = 0; 
  for (pos = 0; pos < ad->N; pos++)
    if (ccoord[pos] != 0) {
      ad->cfrom = ccoord[pos];
      break;
    }
  ad->cto = 0;
  for (pos = 0; pos < ad->N; pos++)
    if (ccoord[pos] != 0) ad->cto = ccoord[pos];

  if(scoord != NULL) free(scoord);
  if(ccoord != NULL) free(ccoord);
  esl_stack_Destroy(pda);

  return ad;

 ERROR:
  if(pda != NULL)    esl_stack_Destroy(pda);
  if(ad != NULL)     cm_alidisplay_Destroy(ad);
  if(scoord != NULL) free(scoord);
  if(ccoord != NULL) free(ccoord);
  return NULL; 
}

/* Function:  cm_alidisplay_Clone()
 * Synopsis:  Make a duplicate of an ALIDISPLAY.
 *
 * Purpose:   Create a duplicate of alignment display <ad>.
 *            Return a pointer to the duplicate. Caller
 *            is responsible for freeing the new object.
 *
 * Returns:   pointer to new <CM_ALIDISPLAY>
 *
 * Throws:    <NULL> on allocation failure.
 */
CM_ALIDISPLAY *
cm_alidisplay_Clone(const CM_ALIDISPLAY *ad)
{
  CM_ALIDISPLAY *ad2 = NULL;
  int status;

  ESL_ALLOC(ad2, sizeof(CM_ALIDISPLAY));
  ad2->model   = ad2->rfline = ad2->nline   = ad2->csline  = ad2->model = ad2->mline = ad2->aseq = ad2->ppline = NULL;
  ad2->cmname  = ad2->cmacc  = ad2->cmdesc  = NULL;
  ad2->sqname  = ad2->sqacc  = ad2->sqdesc  = NULL;
  ad2->mem     = NULL;
  ad2->memsize = 0;

  if (ad->memsize) 		/* serialized */
    {
      ESL_ALLOC(ad2->mem, sizeof(char) * ad->memsize);
      ad2->memsize = ad->memsize;
      memcpy(ad2->mem, ad->mem, ad->memsize);

      ad2->rfline = (ad->rfline ? ad2->mem + (ad->rfline - ad->mem) : NULL );
      ad2->nline  = (ad->nline  ? ad2->mem + (ad->nline - ad->mem)  : NULL );
      ad2->csline = ad2->mem + (ad->csline - ad->mem);
      ad2->model  = ad2->mem + (ad->model  - ad->mem);
      ad2->mline  = ad2->mem + (ad->mline  - ad->mem);
      ad2->aseq   = ad2->mem + (ad->aseq   - ad->mem);
      ad2->ppline = (ad->ppline ? ad2->mem + (ad->ppline - ad->mem) : NULL );
      ad2->N      = ad->N;

      ad2->cmname = ad2->mem + (ad->cmname - ad->mem);
      ad2->cmacc  = ad2->mem + (ad->cmacc  - ad->mem);
      ad2->cmdesc = ad2->mem + (ad->cmdesc - ad->mem);
      ad2->cfrom  = ad->cfrom;
      ad2->cto    = ad->cto;
      ad2->clen   = ad->clen;

      ad2->sqname  = ad2->mem + (ad->sqname - ad->mem);
      ad2->sqacc   = ad2->mem + (ad->sqacc  - ad->mem);
      ad2->sqdesc  = ad2->mem + (ad->sqdesc - ad->mem);
      ad2->sqfrom  = ad->sqfrom;
      ad2->sqto    = ad->sqto;
      ad2->L       = ad->L;
    }
  else				/* deserialized */
    {
      if ( esl_strdup(ad->rfline, -1, &(ad2->rfline)) != eslOK) goto ERROR;
      if ( esl_strdup(ad->nline,  -1, &(ad2->nline))  != eslOK) goto ERROR;
      if ( esl_strdup(ad->csline, -1, &(ad2->csline)) != eslOK) goto ERROR;
      if ( esl_strdup(ad->model,  -1, &(ad2->model))  != eslOK) goto ERROR;
      if ( esl_strdup(ad->mline,  -1, &(ad2->mline))  != eslOK) goto ERROR;
      if ( esl_strdup(ad->aseq,   -1, &(ad2->aseq))   != eslOK) goto ERROR;
      if ( esl_strdup(ad->ppline, -1, &(ad2->ppline)) != eslOK) goto ERROR;
      ad2->N = ad->N;

      if ( esl_strdup(ad->cmname, -1, &(ad2->cmname)) != eslOK) goto ERROR;
      if ( esl_strdup(ad->cmacc,  -1, &(ad2->cmacc))  != eslOK) goto ERROR;
      if ( esl_strdup(ad->cmdesc, -1, &(ad2->cmdesc)) != eslOK) goto ERROR;
      ad2->cfrom  = ad->cfrom;
      ad2->cto    = ad->cto;
      ad2->clen   = ad->clen;

      if ( esl_strdup(ad->sqname,  -1, &(ad2->sqname)) != eslOK) goto ERROR;
      if ( esl_strdup(ad->sqacc,   -1, &(ad2->sqacc))  != eslOK) goto ERROR;
      if ( esl_strdup(ad->sqdesc,  -1, &(ad2->sqdesc)) != eslOK) goto ERROR;
      ad2->sqfrom  = ad->sqfrom;
      ad2->sqto    = ad->sqto;
      ad2->L       = ad->L;      
    }

  return ad2;

 ERROR:
  if (ad2) cm_alidisplay_Destroy(ad2);
  return NULL;
}


/* Function:  cm_alidisplay_Sizeof()
 * Synopsis:  Returns the total size of a CM_ALIDISPLAY, in bytes.
 *
 * Purpose:   Return the total size of <CM_ALIDISPLAY> <ad>, in bytes.
 *
 *            Note that <ad->memsize = cm_alidisplay_Sizeof(ad) - sizeof(CM_ALIDISPLAY)>,
 *            for a serialized object, because <ad->memsize> only refers to the sum
 *            of the variable-length allocated fields.
 *
 * Args:      ad - CM_ALIDISPLAY to get the size of
 *
 * Returns:   size of <ad> in bytes
 */
size_t
cm_alidisplay_Sizeof(const CM_ALIDISPLAY *ad)
{
  size_t n = sizeof(CM_ALIDISPLAY);

  if (ad->rfline) n += ad->N+1; /* +1 for \0 */
  n += 5 * (ad->N+1);           /* nline, csline, model, mline, aseq */
  if (ad->ppline) n += ad->N+1; 
  n += 1 + strlen(ad->cmname);	  
  n += 1 + strlen(ad->cmacc);	/* optional acc, desc fields: when not present, just "" ("\0") */
  n += 1 + strlen(ad->cmdesc);
  n += 1 + strlen(ad->sqname);
  n += 1 + strlen(ad->sqacc);  
  n += 1 + strlen(ad->sqdesc); 
 
  return n;
}

/* Function:  cm_alidisplay_Serialize()
 * Synopsis:  Serialize a CM_ALIDISPLAY, using internal memory.
 *
 * Purpose:   Serialize the <CM_ALIDISPLAY> <ad>, internally converting
 *            all its variable-length allocations to a single
 *            contiguous memory allocation. Serialization aids
 *            interprocess communication.
 *            
 *            If <ad> is already serialized, do nothing.
 *
 * Args:      ad  - alidisplay to serialize
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure, and <ad> is restored to
 *            its original (deserialized) state.
 */
int
cm_alidisplay_Serialize(CM_ALIDISPLAY *ad)
{
  int pos;
  int n;
  int status;

  if (ad->mem) return eslOK;	/* already serialized, so no-op */
  ad->memsize = cm_alidisplay_Sizeof(ad) - sizeof(CM_ALIDISPLAY);
  ESL_ALLOC(ad->mem, ad->memsize);

  /* allow no exceptions past this point, because API guarantees restore of original state upon error */

  pos = 0;
  if (ad->rfline) { memcpy(ad->mem+pos, ad->rfline, ad->N+1); free(ad->rfline); ad->rfline = ad->mem+pos;  pos += ad->N+1; }
  memcpy(ad->mem+pos, ad->nline,  ad->N+1); free(ad->nline);  ad->nline  = ad->mem+pos; pos += ad->N+1; 
  memcpy(ad->mem+pos, ad->csline, ad->N+1); free(ad->csline); ad->csline = ad->mem+pos; pos += ad->N+1; 
  memcpy(ad->mem+pos, ad->model,  ad->N+1); free(ad->model);  ad->model  = ad->mem+pos; pos += ad->N+1; 
  memcpy(ad->mem+pos, ad->mline,  ad->N+1); free(ad->mline);  ad->mline  = ad->mem+pos; pos += ad->N+1; 
  memcpy(ad->mem+pos, ad->aseq,   ad->N+1); free(ad->aseq);   ad->aseq   = ad->mem+pos; pos += ad->N+1; 
  if (ad->ppline) { memcpy(ad->mem+pos, ad->ppline, ad->N+1); free(ad->ppline); ad->ppline = ad->mem+pos;  pos += ad->N+1; }
  n = 1 + strlen(ad->cmname);   memcpy(ad->mem + pos, ad->cmname, n);  free(ad->cmname);  ad->cmname  = ad->mem+pos; pos += n;
  n = 1 + strlen(ad->cmacc);    memcpy(ad->mem + pos, ad->cmacc,  n);  free(ad->cmacc);   ad->cmacc   = ad->mem+pos; pos += n;
  n = 1 + strlen(ad->cmdesc);   memcpy(ad->mem + pos, ad->cmdesc, n);  free(ad->cmdesc);  ad->cmdesc  = ad->mem+pos; pos += n;
  n = 1 + strlen(ad->sqname);   memcpy(ad->mem + pos, ad->sqname,  n); free(ad->sqname);  ad->sqname  = ad->mem+pos; pos += n;
  n = 1 + strlen(ad->sqacc);    memcpy(ad->mem + pos, ad->sqacc,   n); free(ad->sqacc);   ad->sqacc   = ad->mem+pos; pos += n;
  n = 1 + strlen(ad->sqdesc);   memcpy(ad->mem + pos, ad->sqdesc,  n); free(ad->sqdesc);  ad->sqdesc  = ad->mem+pos; pos += n;
  
  return eslOK;

 ERROR:
  if (ad->mem) free(ad->mem); ad->mem = NULL;
  return status;
}

/* Function:  cm_alidisplay_Deserialize()
 * Synopsis:  Deserialize a CM_ALIDISPLAY, using internal memory.
 *
 * Purpose:   Deserialize the <CM_ALIDISPLAY> <ad>, converting its internal
 *            allocations from a single contiguous memory chunk to individual
 *            variable-length allocations. Deserialization facilitates 
 *            reallocation/editing of individual elements of the display.
 *            
 *            If <ad> is already deserialized, do nothing.
 *
 * Args:      ad - alidisplay to serialize
 *
 * Returns:   <eslOK> on success
 *
 * Throws:    <eslEMEM> on allocation failure, and <ad> is restored to
 *            its original (serialized) state.
 */
int
cm_alidisplay_Deserialize(CM_ALIDISPLAY *ad)
{
  int pos;
  int n;
  int status;

  if (ad->mem == NULL) return eslOK; /* already deserialized, so no-op */

  pos = 0;
  if (ad->rfline) { ESL_ALLOC(ad->rfline, sizeof(char) * ad->N+1); memcpy(ad->rfline, ad->mem+pos, ad->N+1); pos += ad->N+1; }
  ESL_ALLOC(ad->nline, sizeof(char) * ad->N+1); memcpy(ad->nline,  ad->mem+pos, ad->N+1); pos += ad->N+1; 
  ESL_ALLOC(ad->csline,sizeof(char) * ad->N+1); memcpy(ad->csline, ad->mem+pos, ad->N+1); pos += ad->N+1; 
  ESL_ALLOC(ad->model, sizeof(char) * ad->N+1); memcpy(ad->model,  ad->mem+pos, ad->N+1); pos += ad->N+1; 
  ESL_ALLOC(ad->mline, sizeof(char) * ad->N+1); memcpy(ad->mline,  ad->mem+pos, ad->N+1); pos += ad->N+1; 
  ESL_ALLOC(ad->aseq,  sizeof(char) * ad->N+1); memcpy(ad->aseq,   ad->mem+pos, ad->N+1); pos += ad->N+1; 
  if (ad->ppline) { ESL_ALLOC(ad->ppline, sizeof(char) * ad->N+1); memcpy(ad->ppline, ad->mem+pos, ad->N+1); pos += ad->N+1; }
  n = 1 + strlen(ad->mem+pos);  ESL_ALLOC(ad->cmname,   sizeof(char) * n); memcpy(ad->cmname,  ad->mem+pos, n); pos += n;
  n = 1 + strlen(ad->mem+pos);  ESL_ALLOC(ad->cmacc,    sizeof(char) * n); memcpy(ad->cmacc,   ad->mem+pos, n); pos += n;
  n = 1 + strlen(ad->mem+pos);  ESL_ALLOC(ad->cmdesc,   sizeof(char) * n); memcpy(ad->cmdesc,  ad->mem+pos, n); pos += n;
  n = 1 + strlen(ad->mem+pos);  ESL_ALLOC(ad->sqname,   sizeof(char) * n); memcpy(ad->sqname,  ad->mem+pos, n); pos += n;
  n = 1 + strlen(ad->mem+pos);  ESL_ALLOC(ad->sqacc,    sizeof(char) * n); memcpy(ad->sqacc,   ad->mem+pos, n); pos += n;
  n = 1 + strlen(ad->mem+pos);  ESL_ALLOC(ad->sqdesc,   sizeof(char) * n); memcpy(ad->sqdesc,  ad->mem+pos, n); pos += n;

  free(ad->mem);
  ad->mem     = NULL;
  ad->memsize = 0;
  return eslOK;
  
 ERROR:
  /* restore serialized state, if an alloc fails. tedious, if not nontrivial. */
  /* the pointers are non-NULL whether we just allocated them or if they're pointing into mem, so we have to check against mem+pos */
  pos = 0;
  if (ad->rfline) { if (ad->rfline != ad->mem+pos) { free(ad->rfline); ad->rfline = ad->mem+pos; }  pos += ad->N+1; }
  if (ad->nline  != ad->mem+pos) { free(ad->nline);  ad->nline  = ad->mem+pos; }  pos += ad->N+1; 
  if (ad->csline != ad->mem+pos) { free(ad->csline); ad->csline = ad->mem+pos; }  pos += ad->N+1; 
  if (ad->model  != ad->mem+pos) { free(ad->model);  ad->model  = ad->mem+pos; }  pos += ad->N+1; 
  if (ad->mline  != ad->mem+pos) { free(ad->mline);  ad->mline  = ad->mem+pos; }  pos += ad->N+1; 
  if (ad->aseq   != ad->mem+pos) { free(ad->aseq);   ad->aseq   = ad->mem+pos; }  pos += ad->N+1; 
  if (ad->ppline) { if (ad->ppline != ad->mem+pos) { free(ad->ppline); ad->ppline = ad->mem+pos; }  pos += ad->N+1; } 

  n = 1 + strlen(ad->cmname);  if (ad->cmname != ad->mem+pos) { free(ad->cmname); ad->cmname = ad->mem+pos;  }  pos += n;
  n = 1 + strlen(ad->cmacc);   if (ad->cmacc  != ad->mem+pos) { free(ad->cmacc);  ad->cmacc  = ad->mem+pos;  }  pos += n;
  n = 1 + strlen(ad->cmname);  if (ad->cmdesc != ad->mem+pos) { free(ad->cmdesc); ad->cmdesc = ad->mem+pos;  }  pos += n;
  n = 1 + strlen(ad->sqname);  if (ad->sqname != ad->mem+pos) { free(ad->sqname); ad->sqname = ad->mem+pos;  }  pos += n;
  n = 1 + strlen(ad->sqacc);   if (ad->sqacc  != ad->mem+pos) { free(ad->sqacc);  ad->sqacc  = ad->mem+pos;  }  pos += n;
  n = 1 + strlen(ad->sqname);  if (ad->sqdesc != ad->mem+pos) { free(ad->sqdesc); ad->sqdesc = ad->mem+pos;  }  pos += n;
  return status;
}


/* Function:  cm_alidisplay_Destroy()
 * Synopsis:  Frees a <CM_ALIDISPLAY>
 */
void
cm_alidisplay_Destroy(CM_ALIDISPLAY *ad)
{
  if (ad == NULL) return;
  if (ad->mem)
    {	/* serialized form */
      free(ad->mem);
    }
  else
    {	/* deserialized form */
      if (ad->rfline)  free(ad->rfline);
      if (ad->nline)   free(ad->nline);
      if (ad->csline)  free(ad->csline);
      if (ad->model)   free(ad->model);
      if (ad->mline)   free(ad->mline);
      if (ad->aseq)    free(ad->aseq);
      if (ad->ppline)  free(ad->ppline);
      if (ad->cmname)  free(ad->cmname);
      if (ad->cmacc)   free(ad->cmacc);
      if (ad->cmdesc)  free(ad->cmdesc);
      if (ad->sqname)  free(ad->sqname);
      if (ad->sqacc)   free(ad->sqacc);
      if (ad->sqdesc)  free(ad->sqdesc);
    }
  free(ad);
}


/* Function: bp_is_canonical
 * Date:     EPN, Wed Oct 14 06:17:27 2009
 *
 * Purpose:  Determine if two residues form a canonical base pair or not.
 *           Works for RNA or DNA (because for some reason cmsearch allows
 *           the user to format output as DNA (with --dna)).
 *
 * Returns:  TRUE if:
 *           lseq  rseq
 *           ----  ----
 *            A     U
 *            U     A
 *            C     G
 *            G     C
 *            G     U
 *            U     G
 *            A     T
 *            T     A
 *            G     T
 *            T     G
 *            Else, return FALSE.
 */
int 
bp_is_canonical(char lseq, char rseq)
{
  switch (toupper(lseq)) { 
  case 'A':
    switch (toupper(rseq)) { 
    case 'U': return TRUE; break;
    case 'T': return TRUE; break;
    default: break;
    }
    break;
  case 'C':
    switch (toupper(rseq)) { 
    case 'G': return TRUE; break;
    default: break;
    }
    break;
  case 'G':
    switch (toupper(rseq)) { 
    case 'C': return TRUE; break;
    case 'U': return TRUE; break;
    case 'T': return TRUE; break;
    default: break;
    }
    break;
  case 'U':
    switch (toupper(rseq)) { 
    case 'A': return TRUE; break;
    case 'G': return TRUE; break;
    default: break;
    }
    break;
  case 'T':
    switch (toupper(rseq)) { 
    case 'A': return TRUE; break;
    case 'G': return TRUE; break;
    default: break;
    }
    break;
  default: break;
  }

  return FALSE;
}

/*---------------- end, alidisplay object -----------------------*/



/*****************************************************************
 * 2. The CM_ALIDISPLAY API
 *****************************************************************/

static int
integer_textwidth(long n)
{
  int w = (n < 0)? 1 : 0;
  while (n != 0) { n /= 10; w++; }
  return w;
}

/* Function:  cm_alidisplay_EncodePostProb()
 * Synopsis:  Convert a posterior probability to a char code.
 *
 * Purpose:   Convert the posterior probability <p> to
 *            a character code suitable for Stockholm format
 *            <#=GC PP_cons> and <#=GR seqname PP> annotation
 *            lines. HMMER uses the same codes in alignment
 *            output.
 *            
 *            Characters <0-9*> are used; $0.0 \leq p < 0.05$
 *            is coded as 0, $0.05 \leq p < 0.15$ is coded as
 *            1, ... and so on ..., $0.85 \leq p < 0.95$ is
 *            coded as 9, and $0.95 \leq p \leq 1.0$ is coded
 *            as '*'.
 *
 * Returns:   the encoded character.
 */
char
cm_alidisplay_EncodePostProb(float p)
{
  return (p + 0.05 >= 1.0) ? '*' :  (char) ((p + 0.05) * 10.0) + '0';
}


/* Function:  cm_alidisplay_DecodePostProb()
 * Synopsis:  Convert a char code post prob to an approx float.
 *
 * Purpose:   Convert posterior probability code <pc>, which
 *            is [0-9*], to an approximate floating point probability.
 *            
 *            The result is crude, because <pc> has already discretized
 *            with loss of precision. We require that 
 *            <cm_alidisplay_EncodePostProb(cm_alidisplay_DecodePostProb(pc)) == pc>,
 *            and that <pc=='0'> decodes to a nonzero probability just to
 *            avoid any possible absorbing-zero artifacts.
 *
 * Returns:   the decoded real-valued approximate probability.
 */
float
cm_alidisplay_DecodePostProb(char pc)
{
  if      (pc == '0') return 0.01;
  else if (pc == '*') return 1.0;
  else if (pc == '.') return 0.0;
  else                return ((float) (pc - '0') / 10.);
}



/* Function:  cm_alidisplay_Print()
 * Synopsis:  Human readable output of <CM_ALIDISPLAY>
 *
 * Purpose:   Prints alignment <ad> to stream <fp>.
 *            
 *            Put at least <min_aliwidth> alignment characters per
 *            line; try to make lines no longer than <linewidth>
 *            characters, including name, coords, and spacing.  The
 *            width of lines may exceed <linewidth>, if that's what it
 *            takes to put a name, coords, and <min_aliwidth>
 *            characters of alignment on a line.
 *            
 *            As a special case, if <linewidth> is negative or 0, then
 *            alignments are formatted in a single block of unlimited
 *            line length.
 *
 * Returns:   eslOK on success
 *            eslEINVAL if ad->aseq or ad->model are invalid,
 *            specifically if local end formatting is invalid.
 *            All local ends should begin with '*[' and end with
 *            ']*' with the intervening characters being 0 or
 *            more whitespace characters followed by an integer.
 */
int
cm_alidisplay_Print(FILE *fp, CM_ALIDISPLAY *ad, int min_aliwidth, int linewidth, int show_accessions, int do_noncanonicals)
{
  char *buf          = NULL;
  char *show_cmname  = NULL;
  char *show_seqname = NULL;
  int   namewidth, coordwidth, aliwidth, cur_aliwidth;
  int   pos;
  int   status;
  int   ni, nk;
  int   z, zp;
  long  i1,i2;
  int   k1,k2;
  int   ni_toadd, nk_toadd; 

  /* implement the --acc option for preferring accessions over names in output  */
  show_cmname  = (show_accessions && ad->cmacc[0] != '\0') ? ad->cmacc : ad->cmname;
  show_seqname = (show_accessions && ad->sqacc[0] != '\0') ? ad->sqacc : ad->sqname;
      
  /* dynamically size the output lines */
  namewidth  = ESL_MAX(strlen(show_cmname), strlen(show_seqname));
  coordwidth = ESL_MAX(ESL_MAX(integer_textwidth(ad->cfrom),
			       integer_textwidth(ad->cto)),
		       ESL_MAX(integer_textwidth(ad->sqfrom),
			       integer_textwidth(ad->sqto)));
  aliwidth   = (linewidth > 0) ? linewidth - namewidth - 2*coordwidth - 5 : ad->N;
  if (aliwidth < ad->N && aliwidth < min_aliwidth) aliwidth = min_aliwidth; /* at least, regardless of some silly linewidth setting */
  ESL_ALLOC(buf, sizeof(char) * (aliwidth+1));
  buf[aliwidth] = '\0';

  /* Break the alignment into multiple blocks of width aliwidth for printing */
  i1 = ad->sqfrom;
  k1 = ad->cfrom ;
  cur_aliwidth = aliwidth; 

  for (pos = 0; pos < ad->N; pos += cur_aliwidth)
    {
      if (pos > 0) fprintf(fp, "\n"); /* blank line betweeen blocks */

      ni = nk = 0; 
      cur_aliwidth = aliwidth; /* this will change if a aliwidth-wide block will end in the middle of a local end display */

      for (z = pos; z < pos + cur_aliwidth && z < ad->N; z++) {
	if (ad->aseq[z]  == '*' && ad->model[z] == '*') { 
	  /* we're at the beginning of a local end display (example: "*[ 7]*"), process it */
	  nk_toadd = ni_toadd = 0;
	  if(ad->aseq[z+1] != '[' || ad->model[z+1] != '[') { status = eslEINVAL; goto ERROR; }

	  zp = z+2;
	  while(ad->model[zp] == ' ') { zp++; } /* chew up any whitespace */
	  while(ad->model[zp] != ']') { nk_toadd *= 10; nk_toadd += ad->model[zp] - '0'; zp++; } /* determine size of local end in model */

	  zp = z+2;
	  while(ad->aseq[zp] == ' ')  { zp++; } /* chew up any whitespace */
	  while(ad->aseq[zp] != ']')  { ni_toadd *= 10; ni_toadd += ad->aseq[zp] - '0';  zp++; } /* determine size of local end in aseq */

	  if((zp+1) >= (pos + aliwidth)) { /* the local end display will not fit completely on this block, save it for next block */
	    cur_aliwidth = z - pos; 
	  }
	  else { 
	    nk += nk_toadd;
	    ni += ni_toadd;
	    z = zp+1; /* position z at end of local end display (one char past the ']', on the '*') */
	  }
	}
	else { /* normal case, we're not at the beginning of a local end */
	  if (ad->model[z] != '.') nk++; /* k advances except on insert states */
	  if (ad->aseq[z]  != '-') ni++; /* i advances except on delete states */
	}
      }

      if(aliwidth != cur_aliwidth) buf[cur_aliwidth] = '\0';

      k2 = k1+nk-1;
      if (ad->sqfrom < ad->sqto) { i2 = i1+ni-1; }
      else                       { i2 = i1-ni+1; }

      if (ad->rfline != NULL) { strncpy(buf, ad->rfline+pos, cur_aliwidth); fprintf(fp, "  %*s %s %*sRF\n", namewidth+coordwidth+1, "", buf, aliwidth-cur_aliwidth, ""); }
      if (do_noncanonicals)   { strncpy(buf, ad->nline+pos,  cur_aliwidth); fprintf(fp, "  %*s %s %*sNC\n", namewidth+coordwidth+1, "", buf, aliwidth-cur_aliwidth, ""); }
      strncpy(buf, ad->csline+pos, cur_aliwidth); fprintf(fp, "  %*s %s %*sCS\n", namewidth+coordwidth+1, "", buf, aliwidth-cur_aliwidth, "");
      strncpy(buf, ad->model+pos,  cur_aliwidth); fprintf(fp, "  %*s %*d %s %*s%-*d\n", namewidth,  show_cmname, coordwidth, k1, buf, aliwidth-cur_aliwidth, "", coordwidth, k2);
      strncpy(buf, ad->mline+pos,  cur_aliwidth); fprintf(fp, "  %*s %s\n", namewidth+coordwidth+1, " ", buf);
      if (ni > 0) { strncpy(buf, ad->aseq+pos, cur_aliwidth); fprintf(fp, "  %*s %*ld %s %*s%-*ld\n", namewidth, show_seqname, coordwidth, i1,  buf, aliwidth-cur_aliwidth, "", coordwidth, i2);  }
      else        { strncpy(buf, ad->aseq+pos, cur_aliwidth); fprintf(fp, "  %*s %*s %s %*s%*s\n",    namewidth, show_seqname, coordwidth, "-", buf, aliwidth-cur_aliwidth, "", coordwidth, "-"); }
      if (ad->ppline != NULL) { strncpy(buf, ad->ppline+pos, cur_aliwidth); fprintf(fp, "  %*s %s %*sPP\n", namewidth+coordwidth+1, "", buf, aliwidth-cur_aliwidth, ""); }

      k1 += nk;
      if (ad->sqfrom < ad->sqto) { i1 += ni; }
      else                       { i1 -= ni; } /* revcomp hit for DNA */
    }
  fflush(fp);
  free(buf);
  return eslOK;

 ERROR:
  if (buf != NULL) free(buf);
  return status;
}  
/*------------------- end, alidisplay API -----------------------*/


/*****************************************************************
 * 3. Debugging/dev code
 *****************************************************************/

/* Function:  cm_alidisplay_Dump()
 * Synopsis:  Print contents of CM_ALIDISPLAY for inspection.
 *
 * Purpose:   Print contents of the <CM_ALIDISPLAY> <ad> to
 *            stream <fp> for inspection. Includes all elements
 *            of the structure, whether the object is allocated
 *            in serialized or deserialized form, and the total
 *            size of the object in bytes.
 *
 * Returns:   <eslOK>
 */
int
cm_alidisplay_Dump(FILE *fp, const CM_ALIDISPLAY *ad)
{
  fprintf(fp, "CM_ALIDISPLAY dump\n");
  fprintf(fp, "------------------\n");

  fprintf(fp, "rfline  = %s\n", ad->rfline ? ad->rfline : "[none]");
  fprintf(fp, "nline   = %s\n", ad->nline  ? ad->nline : "[none]");
  fprintf(fp, "csline  = %s\n", ad->csline ? ad->csline : "[none]");
  fprintf(fp, "model   = %s\n", ad->model);
  fprintf(fp, "mline   = %s\n", ad->mline);
  fprintf(fp, "ppline  = %s\n", ad->ppline ? ad->ppline : "[none]");
  fprintf(fp, "aseq    = %s\n", ad->aseq);
  fprintf(fp, "N       = %d\n", ad->N);
  fprintf(fp, "\n");

  fprintf(fp, "cmname = %s\n", ad->cmname);
  fprintf(fp, "cmacc  = %s\n", ad->cmacc[0]  == '\0' ? "[none]" : ad->cmacc);
  fprintf(fp, "cmdesc = %s\n", ad->cmdesc[0] == '\0' ? "[none]" : ad->cmdesc);
  fprintf(fp, "cfrom  = %d\n", ad->cfrom);
  fprintf(fp, "cto    = %d\n", ad->cto);
  fprintf(fp, "clen   = %d\n", ad->clen);
  fprintf(fp, "\n");

  fprintf(fp, "sqname  = %s\n",  ad->sqname);
  fprintf(fp, "sqacc   = %s\n",  ad->sqacc[0]  == '\0' ? "[none]" : ad->sqacc);
  fprintf(fp, "sqdesc  = %s\n",  ad->sqdesc[0] == '\0' ? "[none]" : ad->sqdesc);
  fprintf(fp, "sqfrom  = %ld\n", ad->sqfrom);
  fprintf(fp, "sqto    = %ld\n", ad->sqto);
  fprintf(fp, "L       = %ld\n", ad->L);
  fprintf(fp, "\n");

  fprintf(fp, "aln_sc  = %.2f\n",ad->aln_sc);
  fprintf(fp, "optacc  = %s\n",  ad->used_optacc ? "TRUE" : "FALSE");
  fprintf(fp, "CYK     = %s\n",  ad->used_optacc ? "FALSE" : "TRUE");
  fprintf(fp, "hbanded = %s\n",  ad->used_hbands ? "TRUE" : "FALSE");
  fprintf(fp, "mx Mb   = %.6f\n",ad->matrix_Mb);
  fprintf(fp, "seconds = %.6f\n",ad->elapsed_secs);
  fprintf(fp, "\n");

  fprintf(fp, "size    = %d bytes\n",  (int) cm_alidisplay_Sizeof(ad));
  fprintf(fp, "%s\n", ad->mem ? "serialized" : "not serialized");
  return eslOK;
}

/* Function:  cm_alidisplay_Compare()
 * Synopsis:  Compare two <CM_ALIDISPLAY> objects for equality
 *
 * Purpose:   Compare alignment displays <ad1> and <ad2> for 
 *            equality. Return <eslOK> if they have identical 
 *            contents; <eslFAIL> if not.
 *            
 *            Only contents matter, not serialization status;
 *            a serialized and deserialized version of the same
 *            alidisplay will compare identical.
 */
int
cm_alidisplay_Compare(const CM_ALIDISPLAY *ad1, const CM_ALIDISPLAY *ad2)
{
  if (ad1->mem && ad2->mem)	/* both objects serialized */
    {
      if (ad1->memsize != ad2->memsize)                  return eslFAIL;
      if (memcmp(ad1->mem, ad2->mem, ad1->memsize) != 0) return eslFAIL;
    }
  
  if (esl_strcmp(ad1->rfline,  ad2->rfline)  != eslOK) return eslFAIL;
  if (esl_strcmp(ad1->csline,  ad2->csline)  != eslOK) return eslFAIL;
  if (esl_strcmp(ad1->model,   ad2->model)   != eslOK) return eslFAIL;
  if (esl_strcmp(ad1->mline,   ad2->mline)   != eslOK) return eslFAIL;
  if (esl_strcmp(ad1->aseq,    ad2->aseq)    != eslOK) return eslFAIL;
  if (esl_strcmp(ad1->ppline,  ad2->ppline)  != eslOK) return eslFAIL;
  if (ad1->N != ad2->N)                                return eslFAIL;

  if (esl_strcmp(ad1->cmname, ad2->cmname) != eslOK) return eslFAIL;
  if (esl_strcmp(ad1->cmacc,  ad2->cmacc)  != eslOK) return eslFAIL;
  if (esl_strcmp(ad1->cmdesc, ad2->cmdesc) != eslOK) return eslFAIL;
  if (ad1->cfrom  != ad2->cfrom)                     return eslFAIL;
  if (ad1->cto    != ad2->cto)                       return eslFAIL;
  if (ad1->clen   != ad2->clen)                      return eslFAIL;

  if (esl_strcmp(ad1->sqname,  ad2->sqname)  != eslOK) return eslFAIL;
  if (esl_strcmp(ad1->sqacc,   ad2->sqacc)   != eslOK) return eslFAIL;
  if (esl_strcmp(ad1->sqdesc,  ad2->sqdesc)  != eslOK) return eslFAIL;
  if (ad1->sqfrom != ad2->sqfrom)                      return eslFAIL;
  if (ad1->sqto   != ad2->sqto)                        return eslFAIL;
  if (ad1->clen   != ad2->clen)                        return eslFAIL;
  
  return eslOK;
}


/*-------------- end, debugging/dev code ------------------------*/



/*****************************************************************
 * 4. Benchmark driver.
 *****************************************************************/
/****************************************************************
 * 5. Unit tests.
 ****************************************************************/
/*****************************************************************
 * 6. Test driver.
 *****************************************************************/
/*****************************************************************
 * 7. Example.
 *****************************************************************/
/*****************************************************************
 * @LICENSE@
 *****************************************************************/

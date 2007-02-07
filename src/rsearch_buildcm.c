/* buildcm.c
 * RJK, March 8, 2002
 * Modified from cmbuild,c, 
 *   SRE, Thu Jul 27 13:19:43 2000 [StL]
 * CVS $Id: buildcm.c,v 1.18 2003/04/08 16:07:29 rjklein Exp $
 * 
 * Given a single sequence and the substitution matrix, constructs a
 * covariance model.
 *  
 *****************************************************************
 * @LICENSE@
 ***************************************************************** 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* external functions                   */
#include "squid.h"		/* general sequence analysis library    */
#include "msa.h"                /* squid's multiple alignment i/o       */
#include "rnamat.h"

/* The six classes */
#define M_cl 0
#define IL_cl 1
#define DL_cl 2
#define IR_cl 3
#define DR_cl 4
#define DB_cl 5

/*
 * calculate_gap_penalty
 *
 * Given the from state, the to state, and the gap parameters, returns
 * the gap penalty.
 */
float calculate_gap_penalty (char from_state, char to_state, 
			     int from_node, int to_node, 
			     float input_alpha, float input_beta, 
			     float input_alphap, float input_betap) {
  int from_class, to_class;
  double alpha, beta;          /* Alpha or beta values to use */

  /* There are potentially 400 different combinations of state pairs here.
     To make it manageable, break down into the 6 classes on p. 8 of lab
     book 7, numbered as follows
     0    M    ROOT_S, BEGL_S, BEGR_S, MATP_MP, MATL_ML, MATR_MR, END_E, BIF_B
     1    IL   ROOT_IL, BEGR_IL, MATP_IL, MATL_IL
     2    DL   MATP_MR, MATL_D
     3    IR   ROOT_IR, MATP_IR, MATR_IR
     4    DR   MATP_ML, MATR_D
     5    DB   MATP_D
  */
  switch (from_state) {
  case MATP_D:
    from_class = DB_cl;
    break;
  case MATP_ML:
  case MATR_D:
    from_class = DR_cl;
    break;
  case ROOT_IR:
  case MATP_IR:
  case MATR_IR:
    from_class = IR_cl;
    break;
  case MATP_MR:
  case MATL_D:
    from_class = DL_cl;
    break;
  case ROOT_IL:
  case BEGR_IL:
  case MATP_IL:
  case MATL_IL:
    from_class = IL_cl;
    break;
  default:
    from_class = M_cl;
  }

  switch (to_state) {
  case MATP_D:
    to_class = DB_cl;
    break;
  case MATP_ML:
  case MATR_D:
    to_class = DR_cl;
    break;
  case ROOT_IR:
  case MATP_IR:
  case MATR_IR:
    to_class = IR_cl;
    break;
  case MATP_MR:
  case MATL_D:
    to_class = DL_cl;
    break;
  case ROOT_IL:
  case BEGR_IL:
  case MATP_IL:
  case MATL_IL:
    to_class = IL_cl;
    break;
  default:
    to_class = M_cl;
  }

  /* Now set alpha and beta according to state classes and nodes */
  /* Alpha is alpha' for MATP->MATP, alpha otherwise */
  if (from_node == MATP_nd && to_node == MATP_nd)
    alpha = input_alphap;
  else
    alpha = input_alpha;
  /* Beta is beta' iff from_cl is DB and MATP->MATP */
  if (from_class == DB_cl && from_node == MATP_nd && to_node == MATP_nd)
    beta = input_betap;
  else 
    beta = input_beta;

  /* Now that we have the proper class, return the appropriate gap penalty */
  if (from_class == M_cl) {
    if (to_class == M_cl) {
      return (0.);
    } else if (to_class == DB_cl) {
      return (alpha);
    } else {
      return (0.5*alpha);
    }
  } else if (from_class == IL_cl) {
    if (to_class == M_cl) {
      return (beta + 0.5*alpha);
    } else if (to_class == IL_cl) {
      return (beta);
    } else if (to_class == DB_cl) {
      return (beta+1.5*alpha);
    } else {
      return (beta + alpha);
    }
  } else if (from_class == DL_cl) {
    if (to_class == M_cl) {
      return (beta + 0.5*alpha);
    } else if (to_class == DL_cl) {
      return (beta);
    } else if (to_class == DB_cl) {
      return (beta + 0.5*alpha);
    } else {
      return (beta + alpha);
    }
  } else if (from_class == IR_cl) {
    if (to_class == M_cl) {
      return (beta + 0.5*alpha);
    } else if (to_class == IR_cl) {
      return (beta);
    } else if (to_class == DB_cl) {
      return (beta+1.5*alpha);
    } else {
      return (beta + alpha);
    }
  } else if (from_class == DR_cl) {
    if (to_class == M_cl) {
      return (beta + 0.5*alpha);
    } else if (to_class == DR_cl) {
      return (beta);
    } else if (to_class == DB_cl) {
      return (beta + 0.5*alpha);
    } else {
      return (beta + alpha);
    }
  } else {                /* DB_cl */
    if (to_class == IL_cl || to_class == IR_cl) {
      return (2*beta + 1.5*alpha);
    } else if (to_class == M_cl) {
      return (2*beta + alpha);
    } else if (to_class == DB_cl) {
      return (2*beta);
    } else {
      return (2*beta + 0.5*alpha);
    }
  }
  return (0);
}

/*
 * SingleSequenceLogoddsify
 *
 * Given a cm and the full matrix, calculates the log-odds scores to use.
 *
 * Assumes that only one emission for each state is greater than 0.
 */
void SingleSequenceLogoddsify (CM_t *cm, fullmat_t *fullmat, float alpha, 
			       float beta, float alphap, float betap) {
  int v, x, y;
  int cur_emission;

  for (v=0; v<cm->M; v++) {
    if (cm->sttype[v] != B_st && cm->sttype[v] != E_st) {
      for (x=0; x<cm->cnum[v]; x++) {
	cm->tsc[v][x] = -1. * calculate_gap_penalty 
	  (cm->stid[v], cm->stid[cm->cfirst[v]+x], 
	   cm->ndtype[cm->ndidx[v]], cm->ndtype[cm->ndidx[cm->cfirst[v]+x]],
	   alpha, beta, alphap, betap);
	/* alphas and betas were positive -- gap score is a penalty, so
	   multiply by -1 */
      }
    }

    if (cm->stid[v] == MATP_MP) {
      /* First, figure out which letter was in the query */
      for (x=0; x<Alphabet_size; x++) {
	for (y=0; y<Alphabet_size; y++) {
	  if (cm->e[v][x*Alphabet_size+y] > 0) {
	    cur_emission = numbered_basepair(Alphabet[x], Alphabet[y]);
	    x = Alphabet_size;
	    y = Alphabet_size;
	  }
	}
      }
      /* Now, calculate the scores from that */
      for (x=0; x<Alphabet_size*Alphabet_size; x++) {
	  cm->esc[v][x] = \
	    fullmat->paired->matrix[matrix_index(cur_emission, x)];
      }
    } else if (cm->stid[v] == MATL_ML || cm->stid[v] == MATR_MR) {
      /* Which letter was in the query */
      for (x=0; x<Alphabet_size; x++) {
	if (cm->e[v][x] > 0) {
	  cur_emission = numbered_nucleotide(Alphabet[x]);
	  x = Alphabet_size;
	}
      }
      for (x=0; x<Alphabet_size; x++) {
	cm->esc[v][x] = \
	  fullmat->unpaired->matrix[matrix_index(cur_emission, x)];
      }
    } else if (cm->stid[v] == MATP_ML || cm->stid[v] == MATP_MR) {
      /* Which letter was in the query */
      for (x=0; x<Alphabet_size; x++) {
	for (y=0; y<Alphabet_size; y++) {
	  /* Get emisison count for MATP state of this state's node */
          /* Depends on states in a given node being ordered always in 
             same order, which is same as how they're indexed in the 
             defines for UNIQUESTATES */
	  if (cm->e[v-(cm->stid[v]-MATP_MP)][x*Alphabet_size + y] > 0) {
	    if (cm->stid[v] == MATP_ML) {
	      cur_emission = numbered_nucleotide(Alphabet[x]);
	    } else {
	      cur_emission = numbered_nucleotide(Alphabet[y]);
	    }
	    x = Alphabet_size;
	    y = Alphabet_size;
	  }
	}
      }
      for (x=0; x<Alphabet_size; x++) {
	cm->esc[v][x] = \
	  fullmat->unpaired->matrix[matrix_index(cur_emission, x)];
      }
    } else if (cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) {
      /* Don't give any score for emissions matching to an Insert state */
      if (cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) {
	for (x = 0; x < Alphabet_size; x++) {
	  cm->esc[v][x] = 0;
	}
      }
    }
  }
  fflush(stdout);
}

void PrintFullCM (CM_t *cm) {
  int i, j;
  printf ("Model nstates: %d\n", cm->M);
  printf ("Model nnodes: %d\n", cm->nodes);
  printf ("Null model:\n");
  for (i = 0; i<Alphabet_size; i++)
    printf ("\t%d: %.2f\n", i, cm->null[i]);
  printf ("--------------\n");
  printf ("state\tsttype\tndidx\tstid\tcfirst\tcnum\tplast\tpnum\n");
  for (i=0; i<cm->M; i++)
    printf ("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", i, cm->sttype[i], cm->ndidx[i], cm->stid[i], cm->cfirst[i], cm->cnum[i], cm->plast[i], cm->pnum[i]);
  printf ("----------------\n");
  printf ("node\tnodemap\tndtype\n");
  for (i=0; i<cm->nodes; i++)
    printf ("%d\t%d\t%d\n", i, cm->nodemap[i], cm->ndtype[i]);
  printf ("-----------------\n");
  for (i=0; i<cm->M; i++) {
    printf ("State %d:\n", i);
    printf ("\tTRANSITION\n");
    printf ("\tnum\tt\ttsc\n");
    for (j=0; j<cm->cnum[i]; j++)
      printf ("\t%d\t%.2f\t%.2f\n", j, cm->t[i][j], cm->tsc[i][j]);
    printf ("\tEMISSION\n");
    printf ("\tnum\te\tesc\n");
    for (j=0; j<Alphabet_size*Alphabet_size; j++)
      printf ("\t%d\t%.2f\t%.2f\n", j, cm->e[i][j], cm->esc[i][j]);
  }
  fflush(stdout);
}

/*
 * Function: SingleSequenceLocalScores
 * Date:     RJK, Tue Jun 25, 2002 [St. Louis]
 * Purpose:  Given a covariance model and a beginsc and endsc value,
 *           fills in beginsc and endsc for each state where begin to or
 *           end from is allowed with appropriate value.
 *           Begin to is only allowed to first state of MATP, MATl, MATR, and
 *              BIF nodes
 *           End from is only allowed from first state of MATP, MATL, MATR,
 *             BEGL, or BEGR node that is not at the end.
 */
void SingleSequenceLocalScores (CM_t *cm, float beginsc, float endsc) {
  int nd,v;

  /* First, initialize all to IMPOSSIBLE */
  for (v=cm->M - 1; v>=0; v--) {
    cm->beginsc[v] = IMPOSSIBLE;
    cm->endsc[v] = IMPOSSIBLE;
  }

  /* beginsc states */
  for (nd = 2; nd < cm->nodes; nd++) {
    if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
        cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BIF_nd)
      cm->beginsc[cm->nodemap[nd]] = beginsc;
  }

  /* endsc states */
  for (nd = 1; nd < cm->nodes; nd++) {
    if ((cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
         cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BEGL_nd ||
         cm->ndtype[nd] == BEGR_nd) &&
        cm->ndtype[nd+1] != END_nd)
      cm->endsc[cm->nodemap[nd]] = endsc;
  }

  cm->flags |= CM_LOCAL_BEGIN;
  cm->flags |= CM_LOCAL_END;
}

/*
 * Function: build_cm()
 * Date:     RJK, Wed Apr 10, 2002
 * Purpose:  Build a CM from the first seq in a multiple alignment file
 *
 * This is taken from the code in cmbuild.c, and therefore is passed 
 * a multiple-sequence alignment structure.  However, at least for now,
 * this structure should only contain 1 sequence.  Specifically, the downstream
 * code from HandmodelMaker expects that for any given state, all emission 
 * probabilities are 0 except for one which is 1.  This requirment will be 
 * met by HandModelmaker for MSAs containing one sequence with weight 1.  
 * This is explicitly checked for here.
 */
CM_t *build_cm (MSA *msa, fullmat_t *fullmat, int *querylen,
		float alpha, float beta, float alphap, float betap,
		float beginsc, float endsc) {
  char           **dsq;		/* digitized aligned sequences             */

  Parsetree_t     *mtr;         /* master structure tree from the alignment*/
  Parsetree_t     *tr;          /* parse tree for the single seq */
  CM_t            *cm;          /* a covariance model                      */

  /* Make sure there's a consensus secondary structure */
  if (msa->ss_cons == NULL) {
    if (msa->ss == NULL || msa->ss[0] == NULL) {
      Die ("No secondary structure given\n");
    } else {
      msa->ss_cons =\
	MallocOrDie (sizeof(char)*(strlen(msa->ss[0]) + 1));
      strncpy (msa->ss_cons, msa->ss[0], strlen(msa->ss[0]));
    }
  }
  
  /* Check the MSA */
  if (msa->nseq != 1 || msa->wgt == NULL)
    Die ("MSA failed to meet single-sequence criteria, %d seqs\n", msa->nseq);
  msa->wgt[0] = 1.;

  /* Set the query length */
  *querylen =      msa->alen;

  /* Digitize the alignment: this takes care of
   * case sensivitity (A vs. a), speeds all future
   * array indexing, and deals with the poor fools
   * who would give us horrid DNA (T) instead of
   * lovely RNA (U). It does cause one wee problem:
   * you need to keep in mind that a digitized seq
   * is indexed 1..alen, but msa (and its annotation!!)
   * is indexed 0..alen-1.
   */
  dsq = DigitizeAlignment(msa->aseq, msa->nseq, msa->alen);
;
  /* Construct a model -- use parameter of 0 for use_rf and 1 for
     gapthresh to force all sequences to "aligned" state.
   */
  HandModelmaker(msa, dsq, 0, 1.0, &cm, &mtr);
  tr = Transmogrify (cm, mtr, dsq[0], msa->aseq[0], msa->alen);
  ParsetreeCount (cm, tr, dsq[0], msa->wgt[0]);
  FreeParsetree(tr);

  /* Now, create log-odds scores that are actually used via matrix 
     and gap penalties*/
  SingleSequenceLogoddsify (cm, fullmat, alpha, beta, alphap, betap);

  /* Infernal specific function: we're writing a CM file for the
   * RSEARCH parameterized CM, so we need valid probabilities (cm->e and
   * cm->t) which we call prob2ascii() on when writing the CM file. To
   * get these we go backwards from the log odds to the probabilities.
   */
  /*SingleSequenceProbify(cm);*/

  /* Put in beginsc and endsc parameters */
  SingleSequenceLocalScores (cm, beginsc, endsc);

  FreeParsetree(mtr);
  Free2DArray((void**)dsq, msa->nseq);

  /*PrintFullCM(cm);*/
  
  return (cm);
}

/* 
 * Function: read_cm
 * Date: RJK, Thu Feb 27, 2003 [St. Louis]
 * Purpose: Given a file name for an INFERNAL (version 0.54) cm, reads it
 *          into the CM structure (in place of the RSEARCH method)
 *          Method lifted from cmsearch.c
 */
CM_t *read_cm (char *queryfile) {
  CM_t *cm;
  CMFILE *cmfp;

  if ((cmfp = CMFileOpen(queryfile, NULL)) == NULL)
    Die("Failed to open covariance model save file %s\n%s\n", queryfile);
  if (! CMFileRead(cmfp, &cm))
    Die("Failed to read a CM from %s -- file corrupt?\n", queryfile);
  if (cm == NULL) 
    Die("%s empty?\n", queryfile);

  ConfigLocal(cm, 0.5, 0.5);
  CMLogoddsify(cm);
  CMHackInsertScores(cm);	/* make insert emissions score zero. "TEMPORARY" FIX. */
  return(cm);
}


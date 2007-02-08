/*
 * rnamat.c
 *
 * Routines for dealing with RNA subsitution/transition matrix.
 *
 * Robert J. Klein
 * February 25, 2002
 */

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "squid.h"
#include "msa.h"
#include "structs.h"
#include "rnamat.h"
#include "sre_stack.h"

/*
 * Maps c as follows:
 * A->0
 * C->1
 * G->2
 * T->3
 * U->3
 * else->-1
 */
int numbered_nucleotide (char c) {
  switch (c) {
  case 'A':
  case 'a':
    return (0);
  case 'C':
  case 'c':
    return (1);
  case 'G':
  case 'g':
    return (2);
  case 'T':
  case 't':
  case 'U':
  case 'u':
    return (3);
  }
  return (-1);
}

/*
 * Maps base pair c,d as follows:
 *
 * AA -> 0
 * AC -> 1
 * ....
 * TG -> 15
 * TT -> 16 (T==U)
 * Anything else maps to -1
 */
int numbered_basepair (char c, char d) {
  int c_num, d_num;
  c_num = numbered_nucleotide (c);
  d_num = numbered_nucleotide (d);
  if (c_num < 0 || d_num < 0) {
    return (-1);
  } else {
    return ((c_num << 2) | d_num);
  }
}

/* Function: rjk_KHS2ct()
 * Incept:   SRE 29 Feb 2000 [Seattle]; from COVE 1.0 code
 * Modified: RJK 27 Feb 2002 [St. Louis]; from Infernal code (rna_ops.c)
 * Purpose:  Convert a secondary structure string (0..len-1) to an array of 
 *           integers representing what position each position is base-paired
 *           to (0..len-1) or -1 if none.  This is a change from what Sean
 *           did in the Infernal code back towards the original way it was
 *           done in the Squid code (compstruct_main.c).  In this case, the
 *           numbering scheme does not match Zuker's .ct files, but does
 *           match the way the MSA is stored using the SQUID library
 *           functions.
 *           
 *           This version does not allow pseudoknots.  Thus ">" and "<" are
 *           used for base pairs, and all other characters, including white
 *           space, are taken to mean unpaired nucleotides.
 *
 * Return:   ret_ct is allocated here and must be free'd by caller.
 *           Returns pointer to ret_ct, or NULL if ss is somehow inconsistent.
 */
int *rjk_KHS2ct(char *ss, int len) {
  Nstack_t *pda;                 
  int      *ct;
  int       pos, pair;

 /* Initialization: always initialize the main pda (0),
   */
  pda = CreateNstack();

  ct = MallocOrDie (len * sizeof(int));
  for (pos = 0; pos < len; pos++)
    ct[pos] = -1;

  for (pos = 0; pos < len; pos++) {
      if (!isprint(ss[pos])) {   /* armor against garbage strings */
	free (ct);
	FreeNstack(pda);
	return (NULL);
      } else if (ss[pos] == '>') {  /* left side of a pair: push onto stack */
        PushNstack(pda, pos);
      } else if (ss[pos] == '<') { /* right side of a pair; resolve pair */
	if (! PopNstack(pda, &pair)) {
	  free (ct);
	  FreeNstack(pda);
	  return (NULL);
	} else {
	  ct[pos]  = pair;
	  ct[pair] = pos;
	}
      }
  }
                                /* nothing should be left on stacks */
  if (! NstackIsEmpty(pda)) {
    free (ct);
    FreeNstack(pda);
    return (NULL);
  }
  FreeNstack(pda);
  return (ct);
}

/*
 * Setup the matrix by allocating matrix in two dimensions as triangle.
 * Initialize to 0.0
 */
matrix_t *setup_matrix (int size) {
  int c;
  matrix_t *mat;

  mat = MallocOrDie(sizeof(matrix_t));
  mat->edge_size = size;
  mat->full_size = matrix_index((size-1),(size-1)) + 1;
  mat->E = 0.0;
  mat->H = 0.0;
  mat->matrix = MallocOrDie (sizeof(double) * mat->full_size);
  for (c=0; c<mat->full_size; c++) {
    mat->matrix[c] = 0.0;
  }
  return(mat);
}

/*
 * middle_of_stem
 *
 * Boolean function, returns TRUE if the gap is in the middle of a stem,
 * false if otherwise.
 *
 * Inputs:
 * msa -- the msa
 * i, j -- indices of the two seqs
 * pos -- the position of the gap in question
 * ct -- array of ct arrays
 */
int middle_of_stem (MSA *msa, int i, int j, int pos, int **ct) {
  int gap_start_pos, gap_stop_pos;

  if (is_rna_gap (msa, i, pos) && is_rna_gap (msa, j, pos)) 
    /* Double gap -- exit */
    return (0);

  if (is_rna_gap (msa, j, pos)) {
    /* Swap so gap is in i */
    gap_start_pos = i;
    i = j;
    j = gap_start_pos;
  }

  /* Now, find start and end positions of the gap */
  for (gap_start_pos = pos; gap_start_pos >= 0; gap_start_pos--) {
    if (!is_rna_gap (msa, i, gap_start_pos)) {
      break;
    }
  }
  for (gap_stop_pos = pos; gap_stop_pos < msa->alen; gap_stop_pos++) {
    if (!is_rna_gap (msa, i, gap_stop_pos)) {
      break;
    }
  }
  if (gap_start_pos < 0 || gap_start_pos >= msa->alen)
    /* Gap at end of alignemtn; can't be internal to stem */
    return (0);

  if (ct[i][gap_start_pos] == 0 || ct[j][gap_start_pos] == 0 || 
      ct[i][gap_stop_pos] == 0 || ct[j][gap_stop_pos] == 0)
    /* Either of the ends not paired */
    return (0);

  if (ct[i][gap_start_pos] != ct[j][gap_start_pos] || 
      ct[i][gap_stop_pos] != ct[j][gap_stop_pos])
    /* The ends not paired homologously */
    return (0);

  return (1);
}

/* TAKEN FROM SQUID's weight.c's simple_distance, but rewritten to
 *  be simple_identity
 * Function: simple_identity()
 * 
 * Purpose:  For two identical-length null-terminated strings, return
 *           the fractional identity between them. (0..1)
 *           (Gaps don't count toward anything.)
 */
static float
simple_identity(char *s1, char *s2)
{
  int diff  = 0;
  int valid = 0;

  for (; *s1 != '\0'; s1++, s2++)
    {
      if (isgap(*s1) || isgap(*s2)) continue;
      if (*s1 == *s2) diff++;
      valid++;
    }
  return (valid > 0 ? ((float) diff / (float) valid) : 0.0);
}
    
  
/*
 * count_matrix
 *
 * Given an MSA and two matrices, counts pairs of column(s) from two sequences
 * at a time into the matrices using the BLOSUM algorithm.
 *
 * Fills in paired matrix (for basepairs), unpaired, background nt count in
 * aligned regions
 *
 * If product weights is false, weight of a pair is average of their weights.
 * If it is true, weight of a pair is product of their weights
 *
 * Each nucleotide at each position can be:
 *    one of gap, unknown character (not ACGTU), known character
 *    one of left bp, right bp, not base paired
 * If both characters are gaps:
 *   Skip
 * If both characters are known:
 *   If both are left bps and match to same right bps
 *     continue
 *   If both are right bps and match to same left bps
 *     add to pairedmat
 *   Otherwise
 *     Add to unpairedmat
 *
 */
void count_matrix (MSA *msa, fullmat_t *fullmat, double *background_nt,
		   int cutoff_perc, int product_weights) {
  int i, j;            /* Seqs we're checking */
  int pos;             /* Column we're counting */
  int prev_pos;        /* Last column we counted */
  float cur_wgt;
  int **ct;            /* ct matrix for all the seqs */
#ifdef REPORT_COUNTED
  int totpair = 0;
  int totsing = 0;
#endif

  /*****************************************
   * 1.  Allocate and fill in ct array
   *****************************************/
  ct = MallocOrDie (sizeof(int *)*msa->nseq);
  if (msa->ss_cons != NULL) {
    ct[0] = rjk_KHS2ct (msa->ss_cons, msa->alen);
  } else {
    ct[0] = rjk_KHS2ct (msa->ss[0], msa->alen);
  }
  for (i=1; i<msa->nseq; i++) {
    if (msa->ss_cons != NULL) {
      ct[i] = ct[0];
    } else {
      ct[i] = rjk_KHS2ct (msa->ss[i], msa->alen);
    }
  }
  for (i=0; i<msa->nseq; i++) {
    if (ct[i] == NULL) {
      Die ("CT string %d is NULL\n", i);
    }
  }

  /**********************************
   * 2.  Count
   **********************************/
  for (i=0; i<msa->nseq; i++) {
    for (j=0; j<i; j++) {
      /* First, make sure it's above the cutoff */
      if (simple_identity(msa->aseq[i], msa->aseq[j]) < 
	  0.01*(float)cutoff_perc)
	continue;       /* Not above cutoff */

      if (product_weights == 1) {
	cur_wgt = msa->wgt[i] * msa->wgt[j];
      } else {
	cur_wgt = (msa->wgt[i] + msa->wgt[j])/2.;
      }

      /* First, skip starting double gaps */
      for (prev_pos = 0; prev_pos <msa->alen &&
	     is_rna_gap (msa, i, prev_pos) && is_rna_gap (msa, j, prev_pos);
	   prev_pos++);

      for (pos=prev_pos; pos<msa->alen; pos++) {
	if (is_defined_rna_nucleotide(msa, i, pos) &&
		   is_defined_rna_nucleotide (msa, j, pos)) {
	  /* If both positions are defined nucleotides */
	  /* If both are left bps and match to same right bps, continue 
             If both are right bps and match to same left bps, add to \
	     pairedmat.  Otherwise, add to unpairedmat */
	  if (ct[i][pos] >= 0 && ct[j][pos] >= 0) {        /* Base pairs */
	    if (ct[i][pos] == ct[j][pos]) {              /* Both equal */
	      if (is_defined_rna_nucleotide(msa, i, ct[i][pos]) && \
		  is_defined_rna_nucleotide (msa, j, ct[j][pos])) {   
		/* Both are RNA nucleotides */
		if (pos < ct[i][pos] && pos <ct[j][pos]) { /* Both left bps */
		  continue;
		} else {                                   /* Both right bps */
		  fullmat->paired->matrix\
		    [matrix_index(numbered_basepair(msa->aseq[i][ct[i][pos]], 
						    msa->aseq[i][pos]), 
				  numbered_basepair(msa->aseq[j][ct[j][pos]], 
						    msa->aseq[j][pos]))] += 
		    cur_wgt;
		  background_nt[numbered_nucleotide
			       (msa->aseq[i][ct[i][pos]])] += cur_wgt;
		  background_nt[numbered_nucleotide
			       (msa->aseq[i][pos])] += cur_wgt;
		  background_nt[numbered_nucleotide
			       (msa->aseq[j][ct[j][pos]])] += cur_wgt;
		  background_nt[numbered_nucleotide
			       (msa->aseq[j][pos])] += cur_wgt;
#ifdef REPORT_COUNTED
		  totpair++;
#endif
		  continue;
		}
	      }
	    }
	  }

	  fullmat->unpaired->matrix
	    [matrix_index(numbered_nucleotide(msa->aseq[i][pos]), 
			  numbered_nucleotide(msa->aseq[j][pos]))] += cur_wgt;
	  background_nt[numbered_nucleotide(msa->aseq[i][pos])] += cur_wgt;
	  background_nt[numbered_nucleotide(msa->aseq[j][pos])] += cur_wgt;
#ifdef REPORT_COUNTED
	  totsing++;
#endif
	}
      }
#ifdef REPORT_COUNTED
      printf ("%d pairs counted\n%d singles counted\n", totpair, totsing);
      totpair = 0;
      totsing = 0;
#endif
    }
  }

  /* Free ct arrays */
  if (ct[0] == ct[1]) {
    free (ct[0]);
  } else {
    for (i=0; i<msa->nseq; i++) {
      free (ct[i]);
    }
  }
  free (ct);
}


/*
 * print_matrix
 *
 * Dumps the paired and unpaired matrices and gap penalties
 */
void print_matrix (FILE *fp, fullmat_t *fullmat) { 
      
  int i, j;

  fprintf (fp, "%s\n\n", fullmat->name);

  fprintf (fp, "    ");
  for (i=0; i<sizeof(RNA_ALPHABET)-1; i++) {
    fprintf (fp, "%c         ", RNA_ALPHABET[i]);
  }
  fprintf (fp, "\n");
  for (i=0; i<sizeof(RNA_ALPHABET)-1; i++) {
    fprintf (fp, "%c   ", RNA_ALPHABET[i]);
    for (j=0; j<=i; j++) {
	fprintf (fp, "%-9.2f ", fullmat->unpaired->matrix[matrix_index(numbered_nucleotide(RNA_ALPHABET[i]), numbered_nucleotide(RNA_ALPHABET[j]))]);
     }
    fprintf (fp, "\n");
  }
  if (strstr (fullmat->name, "RIBOPROB") == NULL)    /* Not probability mat */
    fprintf (fp, "H: %.4f\nE: %.4f\n", fullmat->unpaired->H, fullmat->unpaired->E);

  fprintf (fp, "\n    ");
  for (i=0; i<sizeof(RNAPAIR_ALPHABET)-1; i++) {
    fprintf (fp, "%c%c        ", RNAPAIR_ALPHABET[i], RNAPAIR_ALPHABET2[i]);
  }
  fprintf (fp, "\n");
  for (i=0; i<sizeof(RNAPAIR_ALPHABET)-1; i++) {
    fprintf (fp, "%c%c  ", RNAPAIR_ALPHABET[i], RNAPAIR_ALPHABET2[i]);
    for (j=0; j<=i; j++) {
      fprintf (fp, "%-9.2f ", fullmat->paired->matrix[matrix_index(numbered_basepair(RNAPAIR_ALPHABET[i], RNAPAIR_ALPHABET2[i]), numbered_basepair (RNAPAIR_ALPHABET[j], RNAPAIR_ALPHABET2[j]))]);
    }
    fprintf (fp, "\n");
  }

  if (strstr (fullmat->name, "RIBOPROB") == NULL)    /* Not probability mat */ 
    fprintf (fp, "H: %.4f\nE: %.4f\n", fullmat->paired->H, fullmat->paired->E);
  fprintf (fp, "\n");
}


/*
 * Read the matrix from a file
 */
fullmat_t *ReadMatrix(FILE *matfp) {
  char linebuf[256];
  char fullbuf[16384];
  int fullbuf_used = 0;
  fullmat_t *fullmat;
  int i;
  char *cp, *end_mat_pos;

  fullmat = MallocOrDie (sizeof(fullmat_t));
  fullmat->unpaired = setup_matrix (RNA_ALPHABET_SIZE);
  fullmat->paired = setup_matrix (RNA_ALPHABET_SIZE*RNA_ALPHABET_SIZE);

  while (fgets (linebuf, 255, matfp)) {
    strncpy (fullbuf+fullbuf_used, linebuf, 16384-fullbuf_used-1);
    fullbuf_used += strlen(linebuf);
    if (fullbuf_used >= 16384) {
      Die ("ERROR: Matrix file bigger than 16kb\n");
    }
  }

  /* First, find RIBO, and copy matrix name to fullmat->name */
  cp = strstr (fullbuf, "RIBO");
  for (i = 0; cp[i] && !isspace(cp[i]); i++);   /* Find space after RIBO */
  fullmat->name = MallocOrDie(sizeof(char)*(i+1));
  strncpy (fullmat->name, cp, i);
  fullmat->name[i] = '\0';
  cp = cp + i;

  /* Now, find the first A */
  cp = strchr (cp, 'A');
  fullmat->unpaired->edge_size = 0;
  /* And count how edge size of the matrix */
  while (*cp != '\n' && cp-fullbuf < fullbuf_used) {
    if (!isspace (cp[0]) && isspace (cp[1])) {
      fullmat->unpaired->edge_size++;
    }
    cp++;
  }

  /* Find next A */
  while (*cp != 'A' && (cp-fullbuf) < fullbuf_used) cp++;
  
  /* Take numbers until we hit the H: */
  end_mat_pos = strstr (cp, "H:");
  for (i=0; cp - fullbuf < end_mat_pos-fullbuf; i++) {
    while (!isdigit(*cp) && *cp != '-' && *cp != '.' && \
	   cp-fullbuf < fullbuf_used && cp != end_mat_pos) { 
	cp++;
    }
    if (cp == end_mat_pos)
      break;
    if (cp-fullbuf < fullbuf_used) {
      fullmat->unpaired->matrix[i] = atof(cp);
      while ((isdigit (*cp) || *cp == '-' || *cp == '.') &&\
	     (cp-fullbuf <fullbuf_used)) {
	cp++;
      }
    }
  }
  fullmat->unpaired->full_size = i;

  /* Skip the H: */
  cp += 2;
  fullmat->unpaired->H = atof(cp);

  /* Now, go past the E: */
  cp = strstr (cp, "E:") + 2;
  fullmat->unpaired->E = atof(cp);

  /********* PAIRED MATRIX ************/
  /* Now, find the first A */
  cp = strchr (cp, 'A');
  fullmat->paired->edge_size = 0;
  /* And count how edge size of the matrix */
  while (*cp != '\n') {
    if (!isspace (cp[0]) && isspace (cp[1])) {
      fullmat->paired->edge_size++;
    }
    cp++;
  }

  /* Find next A */
  while (*cp != 'A' && (cp-fullbuf) < fullbuf_used) cp++;

  /* Take numbers until we hit the H: */
  end_mat_pos = strstr (cp, "H:");
  for (i=0; cp - fullbuf < end_mat_pos-fullbuf; i++) {
    while (!isdigit(*cp) && *cp != '-' && *cp != '.' && \
	   cp-fullbuf < fullbuf_used && cp != end_mat_pos) { 
	cp++;
    }
    if (cp == end_mat_pos)
      break;
    if (cp-fullbuf < fullbuf_used) {
      fullmat->paired->matrix[i] = atof(cp);
      while ((isdigit (*cp) || *cp == '-' || *cp == '.') &&\
	     (cp-fullbuf <fullbuf_used)) {
	cp++;
      }
    }
  }
  fullmat->paired->full_size = i;

  /* Skip the H: */
  cp += 2;
  fullmat->paired->H = atof(cp);

  /* Now, go past the E: */
  cp = strstr (cp, "E:") + 2;
  fullmat->paired->E = atof(cp);

  return (fullmat);
}

/*
 * MatFileOpen
 *
 * Given three strings, tries combinations to open the matrix file
 * as follows:
 *
 * MATRIX_DIR = default matrix directory provided at compile-time through
 *   ./configure and setting of data directory ($prefix/share/rsearch/matrices)
 * deflt = default matrix name
 * matdir = matrix directory from RNAMAT environment variable
 * matfile = filenane/matrix name override from -m parameter
 *
 * Order to test:
 * 1.  matfile
 * 2.  matfile.mat
 * 3.  matdir/matfile
 * 4.  matdir/matfile.mat
 * 5.  matdir/default
 * 6.  matdir/default.mat
 * 7.  MATRIX_DIR/matfile
 * 8.  MATRIX_DIR/matfile.mat
 * 9.  MATRIX_DIR/default
 * 10. MATRIX_DIR/default.mat
 */
FILE *MatFileOpen (char *deflt, char *matdir, char *matfile) {
     char buf[1024];
     FILE *fp;
     
     /* Only do 1-4 if matfile defined */
     if (matfile[0] != '\0') {
       /* 1.  matfile */
       fp = fopen (matfile, "r");
       if (fp != NULL) return (fp);

       /* 2.  matfile.mat */
       snprintf (buf, 1023, "%s.mat", matfile);
       buf[1023] = '\0';
       fp = fopen (buf, "r");
       if (fp != NULL) return (fp);

       /* 3. matdir/matfile */
       snprintf (buf, 1023, "%s/%s", matdir, matfile);
       buf[1023] = '\0';
       fp = fopen (buf, "r");
       if (fp != NULL) return (fp);

       /* 4.  matdir/matfile.mat */
       snprintf (buf, 1023, "%s/%s.mat", matdir, matfile);
       buf[1023] = '\0';
       fp = fopen (buf, "r");
       if (fp != NULL) return (fp);

       /* EPN added, we couldn't open the matfile specified
	* at the command line, we should die. 
	*/
       return (NULL);
     }

     /* 5.  matdir/default */
     snprintf (buf, 1023, "%s/%s", matdir, deflt);
     buf[1023] = '\0';
     fp = fopen (buf, "r");
     if (fp != NULL) return (fp);

     /* 6.  matdir/default.mat */
     snprintf (buf, 1023, "%s/%s.mat", matdir, deflt);
     buf[1023] = '\0';
     fp = fopen (buf, "r");
     if (fp != NULL) return (fp);

     /* Only do 7-8 if matfile defined */
     if (matfile[0] != '\0') {
      /* 7. MATRIX_DIR/matfile */
       snprintf (buf, 1023, "%s/%s", MATRIX_DIR, matfile);
       buf[1023] = '\0';
       fp = fopen (buf, "r");
       if (fp != NULL) return (fp);

       /* 8.  MATRIX_DIR/matfile.mat */
       snprintf (buf, 1023, "%s/%s.mat", MATRIX_DIR, matfile);
       buf[1023] = '\0';
       fp = fopen (buf, "r");
       if (fp != NULL) return (fp);
     }

     /* 9.  MATRIX_DIR/default */
     snprintf (buf, 1023, "%s/%s", MATRIX_DIR, deflt);
     buf[1023] = '\0';
     fp = fopen (buf, "r");
     if (fp != NULL) return (fp);

     /* 10.  MATRIX_DIR/default.mat */
     snprintf (buf, 1023, "%s/%s.mat", MATRIX_DIR, deflt);
     buf[1023] = '\0';
     fp = fopen (buf, "r");
     if (fp != NULL) return (fp);

     return (NULL);
}
     
/*
 * Function: get_min_alpha_beta_sum()
 * Date:     RJK, Mon Apr 29, 2002 [St. Louis]
 * Purpose:  Given a full matrix, reports minimum sum allowed 
 *           for alpha and beta (or alpha' and beta')
 *           The maximum for alpha+beta is found by
 *           min Sp(i,j)(k,l) - max{Su(i,k), Su(j,l)}
 *            for all i,j,k.l
 */
float get_min_alpha_beta_sum (fullmat_t *fullmat) {
  float max_sum = 9999999.9;              /* max allowed value of alpha+beta */
  float cur_dif;
  int i,j,k,l;
  int pair_ij, pair_kl;
  
  for (i=0; i<Alphabet_size; i++)
    for (j=0; j<Alphabet_size; j++)
      for (k=0; k<Alphabet_size; k++)
	for (l=0; l<Alphabet_size; l++) {
	  pair_ij = numbered_basepair(Alphabet[i],Alphabet[j]);
	  pair_kl = numbered_basepair(Alphabet[k],Alphabet[l]);
	  /* First check for i,k paired */
	  cur_dif = fullmat->paired->matrix[matrix_index(pair_ij, pair_kl)] -
	    fullmat->unpaired->matrix[matrix_index(i,k)];
	  if (cur_dif < max_sum)
	    max_sum = cur_dif;
	  /* And repeat for j,l */
	  cur_dif = fullmat->paired->matrix[matrix_index(pair_ij, pair_kl)] -
	    fullmat->unpaired->matrix[matrix_index(j,l)];
	  if (cur_dif < max_sum)
	    max_sum = cur_dif;
	}
  return (-1. * max_sum);
}

/* EPN, Tue Feb  6 15:34:00 2007 */
void FreeMat(fullmat_t *fullmat) 
{
  if(fullmat->unpaired != NULL)
    {
      free(fullmat->unpaired->matrix);
      free(fullmat->unpaired);
    }
  if(fullmat->paired != NULL)
    {
      free(fullmat->paired->matrix);
      free(fullmat->paired);
    }
  if(fullmat->name != NULL)
    free(fullmat->name);
  free(fullmat);
}

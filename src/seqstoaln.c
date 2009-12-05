/************************************************************
 * @LICENSE@
 ************************************************************/

/* seqstoaln.c
 * 
 * Support for the seqs_to_aln_t data structure.
 * Originally developed to ease MPI implementation of cmalign.
 * 
 * Note: these functions all originated in dispatch.c.
 * dispatch.c: EPN, Wed Dec  6 06:11:46 2006
 *
 * EPN, Wed Dec  5 13:43:21 2007
 */

#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <math.h>

#include "easel.h"
#include "esl_gumbel.h"
#include "esl_msa.h"         
#include "esl_random.h"         
#include "esl_randomseq.h"         
#include "esl_stack.h"
#include "esl_stopwatch.h"   

#include "funcs.h"		/* external functions                   */
#include "structs.h"		/* data structures, macros, #define's   */

/* Function: CreateSeqsToAln()
 * Date:     EPN, Sat Sep  1 10:51:28 2007
 *
 * Purpose:  Allocate and return a seqs_to_aln_t data structure.
 *
 *           If(i_am_mpi_master), we allocate the seqs_to_aln data
 *           structure differently, b/c we need to keep valid pointers 
 *           for the results (parsetrees, cp9 traces, postcodes) that may
 *           come back from the workers in any order.
 *
 * Returns:  An initialized and allocated (for nalloc seqs) 
 *           seqs_to_aln_t object.
 *           Dies immediately on a memory error.
 */
seqs_to_aln_t *CreateSeqsToAln(int size, int i_am_mpi_master)
{
  int status;
  int i;
  seqs_to_aln_t *seqs_to_aln;

  ESL_ALLOC(seqs_to_aln, sizeof(seqs_to_aln_t));
  ESL_ALLOC(seqs_to_aln->sq,      sizeof(ESL_SQ *)      * size);
  seqs_to_aln->tr       = NULL;
  seqs_to_aln->cp9_tr   = NULL;
  seqs_to_aln->postcode = NULL;
  seqs_to_aln->sc       = NULL;
  seqs_to_aln->pp       = NULL;
  seqs_to_aln->struct_sc= NULL;
  seqs_to_aln->nalloc   = size;
  seqs_to_aln->nseq     = 0;

  if(i_am_mpi_master) {
    ESL_ALLOC(seqs_to_aln->tr,       sizeof(Parsetree_t *) * size);
    ESL_ALLOC(seqs_to_aln->cp9_tr,   sizeof(CP9trace_t)    * size);
    ESL_ALLOC(seqs_to_aln->postcode, sizeof(char **)       * size);
    ESL_ALLOC(seqs_to_aln->sc,       sizeof(float)         * size);
    ESL_ALLOC(seqs_to_aln->pp,       sizeof(float)         * size);
    ESL_ALLOC(seqs_to_aln->struct_sc,sizeof(float)         * size);
    for(i = 0; i < size; i++) {
      seqs_to_aln->tr[i]       = NULL;
      seqs_to_aln->cp9_tr[i]   = NULL;
      seqs_to_aln->postcode[i] = NULL;
      seqs_to_aln->sc[i]       = IMPOSSIBLE;
      seqs_to_aln->pp[i]       = IMPOSSIBLE;
      seqs_to_aln->struct_sc[i]= IMPOSSIBLE;
    }
  }
  return seqs_to_aln;

 ERROR:
  cm_Fail("Memory error.");
  return NULL; /* NEVERREACHED */
}

/* Function: CreateSeqsToAlnFromSq()
 * Date:     EPN, Wed Sep  5 18:13:11 2007
 *
 * Purpose:  Allocate and return a seqs_to_aln_t data structure, setting
 *           seqs_to_aln->sq ptr to an input ptr to ESL_SQs.
 *
 *           If(i_am_mpi_master), we allocate the seqs_to_aln data
 *           structure differently, b/c we need to keep valid pointers 
 *           for the results (parsetrees, cp9 traces, postcodes) that may
 *           come back from the workers in any order.
 *
 * Returns:  An initialized and allocated (for nalloc seqs) 
 *           seqs_to_aln_t object.
 *           Dies immediately on a memory error.
 */
seqs_to_aln_t *CreateSeqsToAlnFromSq(ESL_SQ **sq, int size, int i_am_mpi_master)
{
  int status;
  int i;
  seqs_to_aln_t *seqs_to_aln;

  ESL_ALLOC(seqs_to_aln, sizeof(seqs_to_aln_t));
  seqs_to_aln->sq       = sq;
  seqs_to_aln->tr       = NULL;
  seqs_to_aln->cp9_tr   = NULL;
  seqs_to_aln->postcode = NULL;
  seqs_to_aln->sc       = NULL;
  seqs_to_aln->pp       = NULL;
  seqs_to_aln->struct_sc= NULL;
  seqs_to_aln->nalloc   = size;
  seqs_to_aln->nseq     = size;

  if(i_am_mpi_master) {
    ESL_ALLOC(seqs_to_aln->tr,       sizeof(Parsetree_t *) * size);
    ESL_ALLOC(seqs_to_aln->cp9_tr,   sizeof(CP9trace_t)    * size);
    ESL_ALLOC(seqs_to_aln->postcode, sizeof(char **)       * size);
    ESL_ALLOC(seqs_to_aln->sc,       sizeof(float)         * size);
    ESL_ALLOC(seqs_to_aln->pp,       sizeof(float)         * size);
    ESL_ALLOC(seqs_to_aln->struct_sc,sizeof(float)         * size);
    for(i = 0; i < size; i++) {
      seqs_to_aln->tr[i]       = NULL;
      seqs_to_aln->cp9_tr[i]   = NULL;
      seqs_to_aln->postcode[i] = NULL;
      seqs_to_aln->sc[i]       = IMPOSSIBLE;
      seqs_to_aln->pp[i]       = IMPOSSIBLE;
      seqs_to_aln->struct_sc[i]= IMPOSSIBLE;
    }
  }
  return seqs_to_aln;

 ERROR:
  cm_Fail("Memory error.");
  return NULL; /* NEVERREACHED */
}

/* Function: GrowSeqsToAln()
 * Date:     EPN, Sat Sep  1 11:10:22 2007
 *
 * Purpose:  Grow a seqs_to_aln_t object by <new_alloc>.
 *
 *           If(i_am_mpi_master), we allocate the seqs_to_aln data
 *           structure differently, b/c we need to keep valid pointers 
 *           for the results (parsetrees, cp9 traces, postcodes) that may
 *           come back from the workers in any order.
 *
 * Returns:  eslOK;
 */
int GrowSeqsToAln(seqs_to_aln_t *seqs_to_aln, int new_alloc, int i_am_mpi_master)
{
  int status;
  void *tmp;
  int i;

  ESL_RALLOC(seqs_to_aln->sq, tmp, sizeof(ESL_SQ *) * (seqs_to_aln->nalloc + new_alloc)); 

  if(i_am_mpi_master) {
    ESL_RALLOC(seqs_to_aln->tr,       tmp, sizeof(Parsetree_t *) * (seqs_to_aln->nalloc + new_alloc));
    ESL_RALLOC(seqs_to_aln->cp9_tr,   tmp, sizeof(CP9trace_t)    * (seqs_to_aln->nalloc + new_alloc));
    ESL_RALLOC(seqs_to_aln->postcode, tmp, sizeof(char **)       * (seqs_to_aln->nalloc + new_alloc));
    ESL_RALLOC(seqs_to_aln->sc,       tmp, sizeof(float)         * (seqs_to_aln->nalloc + new_alloc));
    ESL_RALLOC(seqs_to_aln->pp,       tmp, sizeof(float)         * (seqs_to_aln->nalloc + new_alloc));
    ESL_RALLOC(seqs_to_aln->struct_sc,tmp, sizeof(float)         * (seqs_to_aln->nalloc + new_alloc));
    for(i = seqs_to_aln->nalloc; i < (seqs_to_aln->nalloc + new_alloc); i++) {
      seqs_to_aln->tr[i]       = NULL;
      seqs_to_aln->cp9_tr[i]   = NULL;
      seqs_to_aln->postcode[i] = NULL;
      seqs_to_aln->sc[i]       = IMPOSSIBLE;
      seqs_to_aln->pp[i]       = IMPOSSIBLE;
      seqs_to_aln->struct_sc[i]= IMPOSSIBLE;
    }
  }
  
  seqs_to_aln->nalloc += new_alloc;
  return eslOK;

 ERROR:
  cm_Fail("Memory reallocation error.");
  return status; /* NEVERREACHED */
}


/* Function: FreeSeqsToAln()
 *
 * Date:     EPN, Sat Sep  1 11:18:39 2007
 *
 * Purpose:  Free a seqs_to_aln_t object.
 *
 * Returns:  void
 *
 */
void FreeSeqsToAln(seqs_to_aln_t *s) 
{
  int i;
  
  if(s->sq != NULL) /* with MPI workers, we sometimes free the sequences outside this function */
    {
      for (i=0; i < s->nseq; i++) 
	if(s->sq[i] != NULL) esl_sq_Destroy(s->sq[i]);
      free(s->sq);
    }

  if(s->tr != NULL) {
    for (i=0; i < s->nseq; i++)
      if(s->tr[i] != NULL) FreeParsetree(s->tr[i]);
    free(s->tr);
  }

  if(s->cp9_tr != NULL) {
    for (i=0; i < s->nseq; i++)
      if(s->cp9_tr[i] != NULL) CP9FreeTrace(s->cp9_tr[i]);
    free(s->cp9_tr);
  }
 
  if(s->postcode != NULL) {
    for (i=0; i < s->nseq; i++)
      if(s->postcode[i] != NULL) free(s->postcode[i]);
    free(s->postcode);
  }

  if(s->sc != NULL) free(s->sc);
  if(s->pp != NULL) free(s->pp);
  if(s->struct_sc != NULL) free(s->struct_sc);
  
  free(s);
}

/*
 * Function: FreePartialSeqsToAln()
 *
 * Date:     EPN, Wed Sep  5 06:58:39 2007
 *
 * Purpose:  Free specified parts of a seqs_to_aln_t object. 
 *
 * Returns:  void
 *
 */
void FreePartialSeqsToAln(seqs_to_aln_t *s, int do_free_sq, int do_free_tr, int do_free_cp9_tr, int do_free_post, int do_free_sc, int do_free_pp, int do_free_struct_sc) 
{
  int i;
  
  if(do_free_sq && s->sq != NULL) {
    for (i=0; i < s->nseq; i++) 
      if(s->sq[i] != NULL) esl_sq_Destroy(s->sq[i]);
    free(s->sq);
    s->sq = NULL;
  }

  if(do_free_tr && s->tr != NULL) {
    for (i=0; i < s->nseq; i++)
      if(s->tr[i] != NULL) FreeParsetree(s->tr[i]);
    free(s->tr);
    s->tr = NULL;
  }

  if(do_free_cp9_tr && s->cp9_tr != NULL) {
    for (i=0; i < s->nseq; i++)
      if(s->cp9_tr[i] != NULL) CP9FreeTrace(s->cp9_tr[i]);
    free(s->cp9_tr);
    s->cp9_tr = NULL;
  }
 
  if(do_free_post && s->postcode != NULL) {
    for (i=0; i < s->nseq; i++)
      if(s->postcode[i] != NULL) free(s->postcode[i]);
    free(s->postcode);
    s->postcode = NULL;
  }

  if(do_free_sc && s->sc != NULL) {
    free(s->sc);
    s->sc = NULL;
  }

  if(do_free_pp && s->pp != NULL) {
    free(s->pp);
    s->pp = NULL;
  }

  if(do_free_struct_sc && s->struct_sc != NULL) {
    free(s->struct_sc);
    s->struct_sc = NULL;
  }
}

/* Function: ReadSeqsToAln()
 * Date:     EPN, Fri Aug 31 15:20:37 2007
 *
 * Purpose:  Given a pointer to a seq file we're reading seqs to align
 *           from, read in nseq seqs from the seq file, or 
 *           if nseq == 0 && do_real_all == TRUE, read all the seqs.
 *           Add the sequences to a growing seqs_to_aln_t object,
 *           a pointer to which is passed in.
 *
 *           If(i_am_mpi_master), we allocate the seqs_to_aln data
 *           structure differently, b/c we need to keep valid pointers 
 *           for the results (parsetrees, cp9 traces, postcodes) that may
 *           come back from the workers in any order.
 *
 * Returns:  <eslOK> on success with <*ret_seqs_to_aln_t> filled with 
 *           seqs to align, *ret_seqs_to_aln_t->nseq gives number of seqs.
 *           Dies immediately on failure with informative error message.
 */
int ReadSeqsToAln(const ESL_ALPHABET *abc, ESL_SQFILE *seqfp, int nseq, int do_read_all, seqs_to_aln_t *seqs_to_aln, int i_am_mpi_master) 
{
  int status;
  int keep_reading = TRUE;
  int i;
  int nseq_orig;

  /* contract check */
  if(  do_read_all && nseq != 0) cm_Fail("if do_read_all is TRUE,  nseq must be zero.");
  if(! do_read_all && nseq <= 0) cm_Fail("if do_read_all is FALSE, nseq must be at least 1.");

  nseq_orig = seqs_to_aln->nseq;
  i         = seqs_to_aln->nseq;
  if(i == seqs_to_aln->nalloc) GrowSeqsToAln(seqs_to_aln, 100, i_am_mpi_master);

  seqs_to_aln->sq[i] = esl_sq_CreateDigital(abc);
  if(seqs_to_aln->sq[i] == NULL) cm_Fail("Memory allocation error.");
  while (keep_reading && (status = esl_sqio_Read(seqfp, seqs_to_aln->sq[i])) == eslOK) {
    if(seqs_to_aln->sq[i]->n == 0) { esl_sq_Reuse(seqs_to_aln->sq[i]); continue; }
    i++;
    if(i == seqs_to_aln->nalloc) GrowSeqsToAln(seqs_to_aln, 100, i_am_mpi_master);
    if(! do_read_all && (i - nseq_orig) == nseq)   keep_reading = FALSE; 
    seqs_to_aln->sq[i] = esl_sq_CreateDigital(abc);
    if(seqs_to_aln->sq[i] == NULL) cm_Fail("Memory allocation error.");
  }
  /* destroy the last sequence that was alloc'ed but not filled */
  esl_sq_Destroy(seqs_to_aln->sq[i]);
  if ((  do_read_all && status  != eslEOF) || 
      (! do_read_all && (status != eslEOF && status != eslOK)))
    cm_Fail("Parse failed, line %d, file %s:\n%s", 
	    seqfp->linenumber, seqfp->filename, seqfp->errbuf);

  seqs_to_aln->nseq = i;
  return status;

}

/* Function: CMEmitSeqsToAln()
 * Date:     EPN, Tue Sep  4 13:22:11 2007   
 *
 * Purpose:  Create a seqs_to_aln object by generating sequences
 *           from a CM.
 *
 * Note:     Sequences are allocated slightly different if the MPI master
 *           calls this function, to allow us to store them after receiving
 *           them back from workers in any order.
 *
 * Args:     r               - source of randomness
 *           cm              - CM to emit from
 *           ncm             - number for CM (only for naming seqs if cm->name == NULL)
 *           nseq            - number of seqs to emit
 *           padW            - pad W-L residues on either side of the sequence to simulate
 *                             what happens when a hit survives a filter
 *           pdist           - probability distribution to use for to padW, can be NULL if padW == FALSE
 *           i_am_mpi_master - TRUE if called from MPI master (see Note)
 *
 * Returns:  Ptr to a newly allocated seqs_to_aln object with nseq sequences to align.
 *           Dies immediately on failure with informative error message.
 */
seqs_to_aln_t *CMEmitSeqsToAln(ESL_RANDOMNESS *r, CM_t *cm, int ncm, int nseq, int padW, double *pdist, int i_am_mpi_master) 
{
  int status;
  seqs_to_aln_t *seqs_to_aln = NULL;
  char *name = NULL;
  int namelen;
  int L;
  int i;
  char errbuf[cmERRBUFSIZE];
  int padL;
  int half_padL;
  int n, np;
  ESL_DSQ *randdsq = NULL;
  ESL_DSQ *newdsq = NULL;

  if(padW == TRUE && pdist == NULL) cm_Fail("CMEmitSeqsToAln(), padW is TRUE, but pdist is NULL, shouldn't happen\n");

  seqs_to_aln = CreateSeqsToAln(nseq, i_am_mpi_master);

  namelen = IntMaxDigits() + 1;  /* IntMaxDigits() returns number of digits in INT_MAX */
  if(cm->name != NULL) namelen += strlen(cm->name) + 1;
  ESL_ALLOC(name, sizeof(char) * namelen);

  for(i = 0; i < nseq; i++)
    {
      if(cm->name != NULL) sprintf(name, "%s-%d", cm->name, i+1);
      else                 sprintf(name, "%d-%d", ncm, i+1);
      if((status = EmitParsetree(cm, errbuf, r, name, TRUE, NULL, &(seqs_to_aln->sq[i]), &L)) != eslOK) cm_Fail(errbuf);
      while(L == 0 || L > cm->W) { /* if L > cm->W we skip the seq and sample again, this is a hack, but avoids downstream problems with requiring huge HMM banded matrices b/c we went EL emission crazy in one seq */
	esl_sq_Destroy(seqs_to_aln->sq[i]); 
	if((status = EmitParsetree(cm, errbuf, r, name, TRUE, NULL, &(seqs_to_aln->sq[i]), &L)) != eslOK) cm_Fail(errbuf); 
      }
      if(padW) { /* pad the sequence equally on both sides, so that the full length of the entire sequence is (2*cm->W)-L */
	padL      = (cm->W - L) * 2;
	half_padL = padL / 2;
	if(padL > 0) { 
	  ESL_ALLOC(randdsq, sizeof(ESL_DSQ)* (padL+2));
	  ESL_ALLOC(newdsq,  sizeof(ESL_DSQ)* (L+padL+2));
	  if (esl_rsq_xIID(r, pdist, cm->abc->K, padL, randdsq)   != eslOK) cm_Fail("CMEmitSeqsToAln(): failure creating random sequence.");
	  for(n = 0; n <= half_padL; n++) { 
	    np = n;
	    newdsq[np] = randdsq[n];
	  }
	  for(n = 1; n <= L; n++) { 
	    np = half_padL + n;
	    newdsq[np] = seqs_to_aln->sq[i]->dsq[n];
	  }
	  for(n = half_padL+1; n <= (padL+1); n++) { 
	    np = L + n;
	    newdsq[np] = randdsq[n];
	  }	  
	  free(seqs_to_aln->sq[i]->dsq);
	  seqs_to_aln->sq[i]->dsq = newdsq;
	  seqs_to_aln->sq[i]->n = L+padL;
	  free(randdsq);
	}
      }
    }
  seqs_to_aln->nseq = nseq;

  free(name);
  return seqs_to_aln;


 ERROR:
  cm_Fail("memory allocation error");
  return NULL;
}

/* Function: RandomEmitSeqsToAln()
 * Date:     EPN, Tue Sep  4 13:42:16 2007
 *
 * Purpose:  Create a seqs_to_aln object by generating sequences
 *           randomly from a background distro.
 *
 * Note:     Sequences are allocated slightly different if the MPI master
 *           calls this function, to allow us to store them after receiving
 *           them back from workers in any order.
 *
 * Args:     r               - source of randomness
 *           abc             - alphabet 
 *           pdist           - probability distribution to use for emitting
 *           extranum        - use this as first part of sequence name (could be ncm)
 *           nseq            - number of seqs to emit
 *           L_distro        - length distribution (0..Lmax) to draw L (lengths of random seqs) from
 *           Lmax            - maximum length in L_distro
 *           i_am_mpi_master - TRUE if called from MPI master (see Note)
 *
 * Returns:  Ptr to a newly allocated seqs_to_aln object with nseq sequences to align.
 *           Dies immediately on failure with informative error message.
 */
seqs_to_aln_t *RandomEmitSeqsToAln(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, double *pdist, int extranum, int nseq, double *L_distro, int Lmax, int i_am_mpi_master) 
{
  int status;
  seqs_to_aln_t *seqs_to_aln = NULL;
  char *name = NULL;
  int namelen;
  int i;
  int L;

  ESL_DSQ *randdsq = NULL;

  seqs_to_aln = CreateSeqsToAln(nseq, i_am_mpi_master);

  namelen = IntMaxDigits() + 1;  /* IntMaxDigits() returns number of digits in INT_MAX */
  namelen *= 2; /* we'll use two ints in the name below */ 
  namelen += 2; /* for the two '-'s */
  namelen += 7; /* for the two 'randseq's */
  ESL_ALLOC(name, sizeof(char) * namelen);

  for(i = 0; i < nseq; i++)
    {
      sprintf(name, "randseq-%d-%d", extranum, i+1);
      L = esl_rnd_DChoose(r, L_distro, Lmax+1);
      ESL_ALLOC(randdsq, sizeof(ESL_DSQ)* (L+2));
      if (esl_rsq_xIID(r, pdist, abc->K, L, randdsq)  != eslOK) cm_Fail("RandomEmitSeqsToAln(): failure creating random sequence.");
      if((seqs_to_aln->sq[i] = esl_sq_CreateDigitalFrom(abc, name, randdsq, L, NULL, NULL, NULL)) == NULL) 
	 cm_Fail("RandomEmitSeqsToAln() error.");
      free(randdsq);
    }
  seqs_to_aln->nseq = nseq;

  free(name);
  return seqs_to_aln;

 ERROR:
  cm_Fail("memory allocation error");
  return NULL;
}


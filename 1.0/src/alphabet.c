/* alphabet.c
 * SRE, Tue Aug  1 10:31:11 2000 [St. Louis]
 * SVN $Id$
 * 
 * Stuff having to do with manipulating symbols in the (RNA) alphabet.
 * 
 *****************************************************************
 *    This copyrighted source code is freely distributed 
 *    under the terms of the GNU General Public License. See
 *    the files COPYRIGHT and LICENSE for details.
 *****************************************************************  
 */

#include "esl_config.h"
#include "config.h"

#include <string.h>
#include <ctype.h>
#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_vectorops.h"

#include "funcs.h"
#include "structs.h"

/* Function: PairCount()
 * Date:     SRE, Tue Aug  1 10:34:20 2000 [St. Louis]
 *
 * Purpose:  Given a possibly degenerate symbol code for left
 *           and right symbols in a pair, increment a symbol
 *           counter array appropriately.
 *           
 * Args:     abc      - pointer to the internal alphabet
 *           counters - vector to count into [0..abc->K^2-1]
 *           syml     - index of left symbol  [0..abc->sym_iupac-1]
 *           symr     - index of right symbol [0..abc->sym_iupac-1]
 *           wt       - weight to use for the count (often 1.0).          
 *
 * Returns:  void
 */
void
PairCount(const ESL_ALPHABET *abc, float *counters, ESL_DSQ syml, ESL_DSQ symr, float wt)
{
  if (syml < abc->K && symr < abc->K) 
    counters[(int) (syml * abc->K + symr)] += wt;
  else {
    /*    float left = NULL;
	  float right = NULL;
	  ESL_ALLOC(left,  sizeof(float) * esl->abc->K);
	  ESL_ALLOC(right, sizeof(float) * esl->abc->K);*/
    float left[MAXABET];
    float right[MAXABET];

    int   l,r;
    
    esl_vec_FSet(left, MAXABET, 0.);
    esl_vec_FSet(right, MAXABET, 0.);
    esl_abc_FCount(abc, left,  syml, wt);
    esl_abc_FCount(abc, right, symr, wt);

    for (l = 0; l < abc->K; l++)
      for (r = 0; r < abc->K; r++)
	counters[l*abc->K +r] += left[l] * right[r];
  }
  return;
}
float
DegeneratePairScore(const ESL_ALPHABET *abc, float *esc, ESL_DSQ syml, ESL_DSQ symr)
{
  float left[MAXABET], right[MAXABET];
  int l,r;
  float sc;

  if (syml < abc->K && symr < abc->K) 
    return esc[(int) (syml*abc->K+symr)];

  esl_vec_FSet(left, MAXABET, 0.);
  esl_vec_FSet(right, MAXABET, 0.);
  esl_abc_FCount(abc, left,  syml, 1.);
  esl_abc_FCount(abc, right, symr, 1.);
  
  sc = 0.;
  for (l = 0; l < abc->K; l++)
    for (r = 0; r < abc->K; r++)
      sc += esc[l*abc->K+r] * left[l] * right[r];
  return sc;
}
int
iDegeneratePairScore(const ESL_ALPHABET *abc, int *iesc, ESL_DSQ syml, ESL_DSQ symr)
{
  float left[MAXABET], right[MAXABET];
  int l,r;
  float sc;

  if (syml < abc->K && symr < abc->K) 
    return iesc[(int) (syml*abc->K+symr)];

  esl_vec_FSet(left, MAXABET, 0.);
  esl_vec_FSet(right, MAXABET, 0.);
  esl_abc_FCount(abc, left,  syml, 1.);
  esl_abc_FCount(abc, right, symr, 1.);

  sc = 0.;
  for (l = 0; l < MAXABET; l++)
    for (r = 0; r < MAXABET; r++)
      sc += iesc[l*abc->K+r] * left[l] * right[r];
  return (int) sc;
}

/* EPN, Wed Aug 20 13:44:16 2008
 * FastPairScore*() functions: 
 * Written to calculate base pairs scores involving 1 or 2 canonical 
 * residues as efficiently as I know how. 
 * Calculating 'optimized' emission scores, which calculates all possible
 * base pair scores (including non-canonicals) in cm_mx.c:{F,I}CalcOptimizedEmitScores()
 * was a time bottleneck in alignment with P7 and CP9 bands during sub CM construction.
 * sub CMs are created for each target seq, and optimized scores must be 
 * calc'ed for each sub CM.
 */

/* Function:  FastPairScoreBothDegenerate() (float version)
 *           iFastPairScoreBothDegenerate() (int   version)
 * Incept:    EPN, Wed Aug 20 13:18:15 2008
 *
 * Purpose:   Score a base pair given left and right vectors denoting
 *            fractional weight of each canonical for left residue 
 *            and right residue. This function should only be called if 
 *            at least 2 values in both left and right are non-zero,
 *            that is, if both left and right correspond to non-canonical
 *            residues.
 *
 *            Example: if we're calculating for left char 'S' (C or G), 
 *                     and right char 'N' (A,C,G,U), then we'd have:
 *                     left  = [0.00, 0.50, 0.50, 0.00]
 *                     right = [0.25, 0.25, 0.25, 0.25]
 *
 *            K     - alphabet size, abc->K
 *            esc   - emission vector, canonical pair scores are used to calc <sc>
 *            left  - [0..l..K-1] fraction of each canonical residue
 *            right - [0..r..K-1] fraction of each canonical residue
 * 
 * Returns:   <sc> the score of the pair
 */

float
FastPairScoreBothDegenerate(int K, float *esc, float *left, float *right)
{
  int l,r;
  float sc;

  sc = 0.;
  for (l = 0; l < K; l++)
    for (r = 0; r < K; r++)
      sc += esc[l*K+r] * left[l] * right[r];
  return sc;
}
int
iFastPairScoreBothDegenerate(int K, int *iesc, float *left, float *right)
{
  int l,r;
  float sc;

  sc = 0.;
  for (l = 0; l < K; l++)
    for (r = 0; r < K; r++)
      sc += iesc[l*K+r] * left[l] * right[r];
  return (int) sc;
}
/* Function:  FastPairScoreLeftOnlyDegenerate() (float version)
 *           iFastPairScoreLeftOnlyDegenerate() (int   version)
 * Incept:    EPN, Wed Aug 20 13:18:15 2008
 *
 * Purpose:   Score a base pair with left residue non-canonical and
 *            right residue canonical, given a <left> vector denoting
 *            fractional weight of each canonical for left residue. 
 *            This function should only be called if at least 2 values 
 *            left are non-zero, that is, if left corresponds to a 
 *            non-canonical residue.
 *
 *            Example: if we're calculating for left char 'S' (C or G), 
 *                     left  = [0.00, 0.50, 0.50, 0.00]
 *
 *            K     - alphabet size, abc->K
 *            esc   - emission vector, canonical pair scores are used to calc <sc>
 *            left  - [0..l..K-1] fraction of each canonical residue
 *            symr  - canonical right residue must be in range 0..K-1
 * 
 * Returns:   <sc> the score of the pair
 */
float
FastPairScoreLeftOnlyDegenerate(int K, float *esc, float *left, ESL_DSQ symr)
{
  int l;
  float sc;

  sc = 0.;
  for (l = 0; l < K; l++) sc += esc[l*K+symr] * left[l];
  return sc;
}
int
iFastPairScoreLeftOnlyDegenerate(int K, int *iesc, float *left, ESL_DSQ symr)
{
  int l;
  float sc;

  sc = 0.;
  for (l = 0; l < K; l++) sc += iesc[l*K+symr] * left[l];
  return (int) sc;
}
/* Function:  FastPairScoreRightOnlyDegenerate() (float version)
 *           iFastPairScoreRightOnlyDegenerate() (int   version)
 * Incept:    EPN, Wed Aug 20 13:18:15 2008
 *
 * Purpose:   Score a base pair with right residue non-canonical and
 *            right residue canonical, given a <right> vector denoting
 *            fractional weight of each canonical for right residue. 
 *            This function should only be called if at least 2 values 
 *            left are non-zero, that is, if right corresponds to a 
 *            non-canonical residue.
 *
 *            Example: if we're calculating for right char 'M' (A or C),
 *                     right  = [0.50, 0.50, 0.00, 0.00]
 *
 *            K     - alphabet size, abc->K
 *            esc   - emission vector, canonical pair scores are used to calc <sc>
 *            right - [0..l..K-1] fraction of each canonical residue
 *            syml  - canonical left residue, must be in range 0..K-1
 * 
 * Returns:   <sc> the score of the pair
 */
float
FastPairScoreRightOnlyDegenerate(int K, float *esc, float *right, ESL_DSQ syml)
{
  int r;
  float sc;

  sc = 0.;
  for (r = 0; r < K; r++) sc += esc[syml*K+r] * right[r];
  return sc;
}
float
iFastPairScoreRightOnlyDegenerate(int K, int *iesc, float *right, ESL_DSQ syml)
{
  int r;
  float sc;

  sc = 0.;
  for (r = 0; r < K; r++) sc += iesc[syml*K+r] * right[r];
  return (int) sc;
}

/* Function: LeftMarginalScore()
 * Author:   DLK
 *
 * Purpose:  Calculate marginal probability for left half
 *           of an emission pair.  Implicitly assumes
 *           a uniform background distribution
 */
float
LeftMarginalScore(const ESL_ALPHABET *abc, float *esc, ESL_DSQ dres)
{
   float left[MAXABET];
   int i;
   float sc;

   if (dres < abc->K) 
   {
      sc = esl_vec_FLogSum(&(esc[dres*abc->K]),abc->K);
      sc -= sreLOG2(abc->K);
   }
   else /* degenerate */
   {
      esl_vec_FSet(left, MAXABET, 0.);
      esl_abc_FCount(abc, left, dres, 1.);

      sc = 0.;
      for (i = 0; i < MAXABET; i++)
      {
         sc += esl_vec_FLogSum(&(esc[i*abc->K]),abc->K)*left[i];
         sc -= sreLOG2(abc->K)*left[i];
      }
   }

   return sc;
}

/* Function: RightMarginalScore()
 * Author:   DLK
 *
 * Purpose:  Calculate marginal probability for right half
 *           of an emission pair.  Implicitly assumes
 *           a uniform background distribution
 */
float
RightMarginalScore(const ESL_ALPHABET *abc, float *esc, ESL_DSQ dres)
{
   float right[MAXABET];
   int i,j;
   float sc;
   float row[abc->K];

   if (dres < abc->K)
   {
      for (i=0; i<abc->K; i++)
         row[i] = esc[i*abc->K+dres];
      sc = esl_vec_FLogSum(row,abc->K);
      sc -= sreLOG2(abc->K);
   }
   else /* degenerate */
   {
      esl_vec_FSet(right, MAXABET, 0.);
      esl_abc_FCount(abc, right, dres, 1.);

      sc = 0.;
      for (i=0; i < MAXABET; i++)
      {
         for (j=0; j<abc->K; j++)
            row[j] = esc[j*abc->K+dres];
         sc += esl_vec_FLogSum(row,abc->K)*right[i];
         sc -= sreLOG2(abc->K)*right[i];
      }
   }

   return sc;
}

/* Following funcs from RSEARCH, only used by RSEARCH code */

static char random_from_string (ESL_RANDOMNESS *r, char *s);
/* Function: resolve_degenerate
 * Date:     September, 1998 (from hmmgcc)
 * This function resolves "degnerate" nucleotides by selecting a random 
 * A, C, G, or T as appropriate by the code present there.  Returns
 * the character passed in if that character does not represent a
 * non-degnerate nucleotide (either A, C, G, or T or not representative
 * at all of a nucleotide.
 *
 * The degenerate code used here is:
 * (taken from http://www.neb.com/neb/products/REs/RE_code.html
 *
 *                         R = G or A
 *                         K = G or T
 *                         B = not A (C or G or T)
 *                         V = not T (A or C or G)
 *                         Y = C or T
 *                         S = G or C
 *                         D = not C (A or G or T)
 *                         N = A or C or G or T
 *                         M = A or C
 *                         W = A or T
 *                         H = not G (A or C or T)
 *
 * This function assumes all letters are already uppercased via toupper
 * before calling.  In other words, it will return a "n" if passed an "n"
 * because it will assume that the symbol for all nucleotides will be passed
 * in as "N".
 */
char resolve_degenerate (ESL_RANDOMNESS *r, char c) {
  c = toupper(c);
  switch (c) {
    case 'A' : return(c);
    case 'C' : return(c);
    case 'G' : return(c);
    case 'T' : return(c);
    case 'R' : return(random_from_string(r, "GA"));
    case 'K' : return(random_from_string(r, "GT"));
    case 'B' : return(random_from_string(r, "CGT"));
    case 'V' : return(random_from_string(r, "ACG"));
    case 'Y' : return(random_from_string(r, "CT"));
    case 'S' : return(random_from_string(r, "GC"));
    case 'D' : return(random_from_string(r, "AGT"));
    case 'N' : return(random_from_string(r, "ACGT"));
    case 'M' : return(random_from_string(r, "AC"));
    case 'W' : return(random_from_string(r, "AT"));
    case 'H' : return(random_from_string(r, "ACT"));
  }
  return(c);
}

 
/*
 * Function: random_from_string
 * Date:     September, 1998 (approx.) -- from hmmgcc
 * This function returns a character randomly chosen from the string.
 * Used in conjunction with the function below that resolves degenerate code
 * nucleotides.
 */
char random_from_string (ESL_RANDOMNESS *r, char *s) {
  int i;
  do 
    {
      /*i = (int) ((float)(strlen(s)-1)*esl_random(r)/(RAND_MAX+1.0));*/
      i = (int) ((float)(strlen(s)) * esl_random(r));
    } while (i<0 || i>=strlen(s));
  return(s[i]);
}


/* Function: revcomp()
 * Incept:   EPN, Tue Aug  7 10:05:14 2007
 *           based on Squid's revcomp()
 *
 * Purpose:  Reverse complement ESL_SQ seq; store in comp.
 *           Can revcomp "in place" (revcomp(seq, seq)).
 *           sq can be in digital or text form.
 *
 * Args:     comp  - destination for reverse complement of sq
 *           seq   - sequence to reverse complement
 *
 * Returns:  eslOK on success;
 *           Dies immediately if any error occurs.
 */
int
revcomp(const ESL_ALPHABET *abc, ESL_SQ *comp, ESL_SQ *sq)
{
  int status;
  int do_digital = FALSE;
  int i;

  /* contract checks */
  if (comp == NULL)
    cm_Fail("ERROR in revcomp, comp is NULL.");
  if(sq == NULL)
    cm_Fail("ERROR in revcomp, sq is NULL.");
  if(sq->dsq != NULL && comp->dsq == NULL)
    cm_Fail("ERROR in revcomp, sq is digital, comp is not.");
  if(sq->dsq == NULL && comp->dsq != NULL)
    cm_Fail("ERROR in revcomp, comp is digital, sq is not.");
  if(abc->type != eslRNA && abc->type != eslDNA)
    cm_Fail("ERROR in revcomp, alphabet type must be RNA or DNA.");
  if(comp->n < sq->n)
    cm_Fail("ERROR in revcomp, comp->n is smaller than sq->n.");

  if(sq->dsq != NULL) do_digital = TRUE;

  if(do_digital) {
    if((status = esl_rsq_XReverse(sq->dsq, sq->n, comp->dsq)) != eslOK) 
      goto ERROR; 
  }
  else {
    if((status = esl_rsq_CReverse(sq->seq, comp->seq) != eslOK))
      goto ERROR; 
  } 

  if(do_digital)
    {
      for(i = 1; i <= sq->n; i++)
	{
	  if(sq->dsq[i] >= abc->Kp) { status = eslEINVAL; goto ERROR; }
	  switch (abc->sym[sq->dsq[i]]) {
	  case 'A': 
	    if(abc->type == eslRNA) 
	      comp->dsq[i] = abc->inmap[(int) 'U']; 
	    else
	      comp->dsq[i] = abc->inmap[(int) 'T']; 
	    break;
	  case 'C': comp->dsq[i] = abc->inmap[(int) 'G']; break;
	  case 'G': comp->dsq[i] = abc->inmap[(int) 'C']; break;
	  case 'T': comp->dsq[i] = abc->inmap[(int) 'A']; break;
	  case 'U': comp->dsq[i] = abc->inmap[(int) 'A']; break;
	  case 'R': comp->dsq[i] = abc->inmap[(int) 'Y']; break;
	  case 'Y': comp->dsq[i] = abc->inmap[(int) 'R']; break;
	  case 'M': comp->dsq[i] = abc->inmap[(int) 'K']; break;
	  case 'K': comp->dsq[i] = abc->inmap[(int) 'M']; break;
	  case 'S': comp->dsq[i] = abc->inmap[(int) 'S']; break;
	  case 'W': comp->dsq[i] = abc->inmap[(int) 'W']; break;
	  case 'H': comp->dsq[i] = abc->inmap[(int) 'D']; break;
	  case 'D': comp->dsq[i] = abc->inmap[(int) 'H']; break;
	  case 'B': comp->dsq[i] = abc->inmap[(int) 'V']; break;
	  case 'V': comp->dsq[i] = abc->inmap[(int) 'B']; break;
	  default:  break;		/* anything else? leave it; it's prob a gap or an X */
	  }
	}
    }
  else
    {
      for(i = 0; i < sq->n; i++)
	{
	  if(islower(sq->seq[i])) { status = eslEINVAL; goto ERROR; }
	     switch (sq->seq[i]) {
	     case 'A': 
	       if(abc->type == eslRNA) 
		 comp->seq[i] = 'U';
	       else
		 comp->seq[i] = 'T'; 
	       break;
	     case 'C': comp->seq[i] = 'G'; break;
	     case 'G': comp->seq[i] = 'C'; break;
	     case 'T': comp->seq[i] = 'A'; break;
	     case 'U': comp->seq[i] = 'A'; break;
	     case 'R': comp->seq[i] = 'Y'; break;
	     case 'Y': comp->seq[i] = 'R'; break;
	     case 'M': comp->seq[i] = 'K'; break;
	     case 'K': comp->seq[i] = 'M'; break;
	     case 'S': comp->seq[i] = 'S'; break;
	     case 'W': comp->seq[i] = 'W'; break;
	     case 'H': comp->seq[i] = 'D'; break;
	     case 'D': comp->seq[i] = 'H'; break;
	     case 'B': comp->seq[i] = 'V'; break;
	     case 'V': comp->seq[i] = 'B'; break;
	     default:  break;		/* anything else? leave it; it's prob a gap or an X */
	     }
	}
    }
  return eslOK;

 ERROR: 
  cm_Fail("Unexpected error code: %d in revcomp().", status);
  return status; /* NOTREACHED */
}
    

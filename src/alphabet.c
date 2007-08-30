/* alphabet.c
 * SRE, Tue Aug  1 10:31:11 2000 [St. Louis]
 * SVN $Id$
 * 
 * Stuff having to do with manipulating symbols in the (RNA) alphabet.
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************  
 */

#include "esl_config.h"
#include "config.h"

#include <string.h>
#include <ctype.h>
#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"
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
PairCount(const ESL_ALPHABET *abc, float *counters, char syml, char symr, float wt)
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
DegeneratePairScore(const ESL_ALPHABET *abc, float *esc, char syml, char symr)
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
iDegeneratePairScore(const ESL_ALPHABET *abc, int *iesc, char syml, char symr)
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

/* Function: LeftMarginalScore()
 * Author:   DLK
 *
 * Purpose:  Calculate marginal probability for left half
 *           of an emission pair.  Implicitly assumes
 *           a uniform background distribution
 */
float
LeftMarginalScore(const ESL_ALPHABET *abc, float *esc, int dres)
{
   float sc;
   sc = esl_vec_FLogSum(&(esc[dres*abc->K]),abc->K);
   sc -= sreLOG2(abc->K);
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
RightMarginalScore(const ESL_ALPHABET *abc, float *esc, int dres)
{
   int i;
   float sc;
   float row[abc->K];
   for (i=0; i<abc->K; i++)
      row[i] = esc[i*abc->K+dres];
   sc = esl_vec_FLogSum(row,abc->K);
   sc -= sreLOG2(abc->K);
   return sc;
}


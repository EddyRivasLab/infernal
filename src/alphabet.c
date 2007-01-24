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

#include <string.h>
#include <ctype.h>

#include "squid.h"
#include "vectorops.h"

#include "easel.h"
#include "esl_vectorops.h"

#include "structs.h"

/* Function: SymbolIndex()
 * 
 * Purpose:  Convert a symbol to its index in Alphabet[].
 *           Bogus characters are silently converted to 'N'.
 *           More robust than the SYMIDX() macro but
 *           presumably slower.
 */ 
char
SymbolIndex(char sym)
{
  char *s;
  return ((s = strchr(Alphabet, (char) toupper((int) sym))) == NULL) ?
	  (char) (Alphabet_iupac-1) : (char) (s - Alphabet);
} 


/* Function: SingletCount()
 * Date:     SRE, Tue Aug  1 10:33:25 2000 [St. Louis]
 *
 * Purpose:  Given a possibly degenerate symbol code, increment
 *           a symbol counter array (generally an emission
 *           probability vector in counts form) appropriately.
 *           
 * Args:     counters:  vector to count into. [0..Alphabet_size-1]
 *           symidx:    symbol index to count: [0..Alphabet_iupac-1]
 *           wt:        weight to use for the count; often 1.0
 *           
 * Return:   (void)    
 */
void
SingletCount(float *counters, char symidx, float wt)
{
  int x;

  if (symidx < Alphabet_size) 
    counters[(int) symidx] += wt;
  else
    for (x = 0; x < Alphabet_size; x++) {
      if (Degenerate[(int) symidx][x])
	counters[x] += wt / (float) DegenCount[(int) symidx];
    }
}

/* Function: PairCount()
 * Date:     SRE, Tue Aug  1 10:34:20 2000 [St. Louis]
 *
 * Purpose:  Given a possibly degenerate symbol code for left
 *           and right symbols in a pair, increment a symbol
 *           counter array appropriately.
 *           
 * Args:     counters - vector to count into [0..Alphabet_size^2-1]
 *           syml     - index of left symbol  [0..Alphabet_iupac-1]
 *           symr     - index of right symbol [0..Alphabet_iupac-1]
 *           wt       - weight to use for the count (often 1.0).          
 *
 * Returns:  void
 */
void
PairCount(float *counters, char syml, char symr, float wt)
{
  if (syml < Alphabet_size && symr < Alphabet_size) 
    counters[(int) (syml * Alphabet_size + symr)] += wt;
  else {
    float left[MAXABET],right[MAXABET];
    int   l,r;
    
    FSet(left, MAXABET, 0.);
    FSet(right, MAXABET, 0.);
    SingletCount(left, syml, wt);
    SingletCount(right, symr, wt);

    for (l = 0; l < Alphabet_size; l++)
      for (r = 0; r < Alphabet_size; r++)
	counters[l*Alphabet_size +r] += left[l] * right[r];
  }
}
  

float
DegeneratePairScore(float *esc, char syml, char symr)
{
  float left[MAXABET], right[MAXABET];
  int l,r;
  float sc;

  if (syml < Alphabet_size && symr < Alphabet_size) 
    return esc[(int) (syml*Alphabet_size+symr)];

  FSet(left, MAXABET, 0.);
  FSet(right, MAXABET, 0.);
  SingletCount(left, syml, 1.);
  SingletCount(right, symr, 1.);

  sc = 0.;
  for (l = 0; l < Alphabet_size; l++)
    for (r = 0; r < Alphabet_size; r++)
      sc += esc[l*Alphabet_size+r] * left[l] * right[r];
  return sc;
}
int
iDegeneratePairScore(int *iesc, char syml, char symr)
{
  float left[MAXABET], right[MAXABET];
  int l,r;
  float sc;

  if (syml < Alphabet_size && symr < Alphabet_size) 
    return iesc[(int) (syml*Alphabet_size+symr)];

  FSet(left, MAXABET, 0.);
  FSet(right, MAXABET, 0.);
  SingletCount(left, syml, 1.);
  SingletCount(right, symr, 1.);

  sc = 0.;
  for (l = 0; l < Alphabet_size; l++)
    for (r = 0; r < Alphabet_size; r++)
      sc += iesc[l*Alphabet_size+r] * left[l] * right[r];
  return (int) sc;
}
float 
DegenerateSingletScore(float *esc, char sym)
{
  float nt[MAXABET];		
  float sc;
  int   x;

  if (sym < Alphabet_size) return esc[(int) sym];

  FSet(nt, MAXABET, 0.);
  SingletCount(nt, sym, 1.);
  sc = 0.;
  for (x = 0; x < Alphabet_size; x++)
    sc += esc[x] * nt[x];
  return sc;
}
int
iDegenerateSingletScore(int *iesc, char sym)
{
  float nt[MAXABET];		
  float sc;
  int   x;

  if (sym < Alphabet_size) return iesc[(int) sym];

  FSet(nt, MAXABET, 0.);
  SingletCount(nt, sym, 1.);
  sc = 0.;
  for (x = 0; x < Alphabet_size; x++)
    sc += iesc[x] * nt[x];
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
LeftMarginalScore(float *esc, int dres)
{
   float sc;
   sc = esl_vec_FLogSum(&(esc[dres*Alphabet_size]),Alphabet_size);
   sc -= log(Alphabet_size);
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
RightMarginalScore(float *esc, int dres)
{
   int i;
   float sc;
   float row[Alphabet_size];
   for (i=0; i<Alphabet_size; i++)
      row[i] = esc[i*Alphabet_size+dres];
   sc = esl_vec_FLogSum(row,Alphabet_size);
   sc -= log(Alphabet_size);
   return sc;
}

/* Function: DigitizeSequence()
 * Date:     SRE, Wed Aug  2 13:05:49 2000 [St. Louis]
 *
 * Purpose:  Digitize a sequence in preparation for a DP algorithm.
 *           a dsq is 1..L, with 0 and L+1 filled with flag bytes.
 *             values in dsq:  
 *             0..Alphabet_iupac-1:      Symbol index.
 *             DIGITAL_GAP      (126):   gap symbol.
 *             DIGITAL_SENTINEL (127):   end bytes 0,L+1
 */
char *
DigitizeSequence(char *seq, int L)
{
  char *dsq;
  int   i;
  char  c;

  dsq = MallocOrDie(sizeof(char) * (L+2));
  dsq[0] = dsq[L+1] = DIGITAL_SENTINEL;
  for (i = 0; i < L; i++) {
    if (isgap(seq[i])) dsq[i+1] = DIGITAL_GAP; 
    else {
      c = toupper((int) seq[i]);
      if (c == 'T') c = 'U';	/* it's RNA, dammit. */
      dsq[i+1] = SymbolIndex(c);
    }
  }
  return dsq;
}

/* Function: DigitizeAlignment()
 * Date:     SRE, Sat Aug  5 18:16:21 2000 [St. Louis]
 *
 * Purpose:  Convert aseqs to digitized sequences.
 *           As with unaligned seqs, while the raw sequences
 *           are 0..L-1, the digitized seqs are 1..L.
 *           Gaps are digitized as a 126 (127 is used for
 *           the sentinel byte at 0 and L+1).
 *
 * Args:     aseq -- [0..nseq-1][0..alen-1] array of seqs
 *           nseq -- number of sequences
 *           alen -- length of alignment
 *
 * Returns:  **dsq -- digitized aligned sequences
 */
char **
DigitizeAlignment(char **aseq, int nseq, int alen)
{
  char **dsq;
  int    idx;
  
  dsq = MallocOrDie(sizeof(char *) * nseq);
  for (idx = 0; idx < nseq; idx++)
    dsq[idx] = DigitizeSequence(aseq[idx], alen);
  return dsq;
}

/* alphabet.c
 * SRE, Tue Aug  1 10:31:11 2000 [St. Louis]
 * CVS $Id$
 * 
 * Stuff having to do with manipulating symbols in the (RNA) alphabet.
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************  
 */

#include <string.h>
#include <ctype.h>

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
    
    SingletCount(left, syml, wt);
    SingletCount(right, symr, wt);

    for (l = 0; l < Alphabet_size; l++)
      for (r = 0; r < Alphabet_size; r++)
	counters[l*Alphabet_size +r] += left[l] * right[r];
  }
}
  

/* globals.c
 * SRE 28 Feb 2000
 * RCS $Id$
 * 
 * Settings of the global variables used in INFERNAL,
 * such as the sequence alphabet. 
 *
 * These globals are supposed
 * to enable run-time configuration *once* -- they are not
 * to be touched after initial program startup and configuration.
 * for example, they are all assumed to be threadsafe.
 ************************************************************
 * @LICENSE@
 ************************************************************
 */

#include "squid.h"
#include "structs.h"

/* The symbol alphabet.
 * The package is designed to be configurable for analysis of
 * a different alphabet just by changing these globals. I doubt that 
 * it would be *useful* to apply it to something other than RNA,
 * but hey, HMMER got used to do signal processing on auto engine
 * electronics and musical compositions, so what the hell do I know.
 *
 * Must deal with IUPAC degeneracies. Nondegenerate symbols
 * come first in Alphabet[], followed by degenerate symbols.
 * We also have to deal with some non-IUPAC stuff:
 * T (bah!) and X (often misused for N).
 *
 * Parts of the code assume that the last symbol is a
 * symbol for an unknown residue, i.e. 'N'.
 */
int   Alphabet_type  = kRNA;
int   Alphabet_size  = 4;
int   Alphabet_iupac = 17;
char *Alphabet       = "ACGUTXRYMKSWHBVDN";

char  Degenerate[17][4] = { 
  /* A  C  G  T */
  {  1, 0, 0, 0 },  /* A */
  {  0, 1, 0, 0 },  /* C */
  {  0, 0, 1, 0 },  /* G */
  {  0, 0, 0, 1 },  /* U */
  {  0, 0, 0, 1 },  /* T */
  {  1, 1, 1, 1 },  /* X */  
  {  1, 0, 1, 0 },  /* R */  
  {  0, 1, 0, 1 },  /* Y */
  {  1, 1, 0, 0 },  /* M */
  {  0, 0, 1, 1 },  /* K */
  {  0, 1, 1, 0 },  /* S */
  {  1, 0, 0, 1 },  /* W */
  {  1, 1, 0, 1 },  /* H */
  {  0, 1, 1, 1 },  /* B */
  {  1, 1, 1, 0 },  /* V */
  {  1, 0, 1, 1 },  /* D */
  {  1, 1, 1, 1 },  /* N */
};
                     /*A  C  G  T  U  X  R  Y  M  K  S  W  H  B  V  D  N*/
int DegenCount[17] = { 1, 1, 1, 1, 1, 4, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4 };


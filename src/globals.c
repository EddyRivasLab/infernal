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

/* CMTransitionIndex[][]
 * Defines a lookup table that lets us turn two unique state types into
 * a transition vector index in a cm --
 * e.g. the value of t_v(x) is in cm->t[v][CMTransitionIndex[v][x]]
 */
int CMTransitionIndex[UNIQUESTATES][UNIQUESTATES] = {
  /* ROOT:S IL IR BEGL:S BEGR: S IL MATP:MP ML MR D IL IR MATL:ML D IL MATR: MR D IR END:E BIF:B
  {      -1, 0, 1,    -1,     -1,-1,      2, 3, 4,5,-1,-1,      2,3,-1,       2,3,-1,    2,    2 }, /* ROOT_S  */
  {      -1, 0, 1,    -1,     -1,-1,      2, 3, 4,5,-1,-1,      2,3,-1,       2,3,-1,    2,    2 }, /* ROOT_IL */
  {      -1,-1, 0,    -1,     -1,-1,      1, 2, 3,4,-1,-1,      1,2,-1,       1,2,-1,    1,    1 }, /* ROOT_IR */
  {      -1,-1,-1,    -1,     -1,-1,      0, 1, 2,3,-1,-1,      0,1,-1,       0,1,-1,    0,    0 }, /* BEGL_S  */
  {      -1,-1,-1,    -1,     -1, 0,      1, 2, 3,4,-1,-1,      1,2,-1,       1,2,-1,    1,    1 }, /* BEGR_S  */ 
  {      -1,-1,-1,    -1,     -1, 0,      1, 2, 3,4,-1,-1,      1,2,-1,       1,2,-1,    1,    1 }, /* BEGR_IL */ 
  {      -1,-1,-1,    -1,     -1,-1,      2, 3, 4,5, 0, 1,      2,3,-1,       2,3,-1,    2,    2 }, /* MATP_MP */
  {      -1,-1,-1,    -1,     -1,-1,      2, 3, 4,5, 0, 1,      2,3,-1,       2,3,-1,    2,    2 }, /* MATP_ML */
  {      -1,-1,-1,    -1,     -1,-1,      2, 3, 4,5, 0, 1,      2,3,-1,       2,3,-1,    2,    2 }, /* MATP_MR */
  {      -1,-1,-1,    -1,     -1,-1,      2, 3, 4,5, 0, 1,      2,3,-1,       2,3,-1,    2,    2 }, /* MATP_D  */
  {      -1,-1,-1,    -1,     -1,-1,      2, 3, 4,5, 0, 1,      2,3,-1,       2,3,-1,    2,    2 }, /* MATP_IL */
  {      -1,-1,-1,    -1,     -1,-1,      1, 2, 3,4,-1, 0,      1,2,-1,       1,2,-1,    1,    1 }, /* MATP_IR */
  {      -1,-1,-1,    -1,     -1,-1,      1, 2, 3, 4,-1,-1,     1,2,0,        1,2,-1,    1,    1 }, /* MATL_ML */
  {      -1,-1,-1,    -1,     -1,-1,      1, 2, 3, 4,-1,-1,     1,2,0,        1,2,-1,    1,    1 }, /* MATL_D  */
  {      -1,-1,-1,    -1,     -1,-1,      1, 2, 3, 4,-1,-1,     1,2,0,        1,2,-1,    1,    1 }, /* MATL_IL */
  {      -1,-1,-1,    -1,     -1,-1,      1, 2, 3, 4,-1,-1,     1,2,-1,       1,2,0,     1,    1 }, /* MATR_MR */
  {      -1,-1,-1,    -1,     -1,-1,      1, 2, 3, 4,-1,-1,     1,2,-1,       1,2,0,     1,    1 }, /* MATR_D  */
  {      -1,-1,-1,    -1,     -1,-1,      1, 2, 3, 4,-1,-1,     1,2,-1,       1,2,0,     1,    1 }, /* MATR_IR */
  {      -1,-1,-1,    -1,     -1,-1,     -1,-1,-1,-1,-1,-1,    -1,-1,-1,     -1,-1,-1,  -1,   -1 }, /* END_E   */
  {      -1,-1,-1,     0,      0,-1,     -1,-1,-1,-1,-1,-1,    -1,-1,-1,     -1,-1,-1,  -1,   -1 }, /* BIF_B   */
};

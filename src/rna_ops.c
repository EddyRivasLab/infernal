/* rna_ops.c
 * SRE, Fri Jul 28 09:37:18 2000 [St. Louis]
 * CVS $Id$
 * 
 * Common functions for manipulating RNA information.
 * 
 ***************************************************************** 
 * @LICENSE@
 ***************************************************************** 
 */

#include <ctype.h>

#include "squid.h"
#include "nstack.h"

/* Function: KHS2ct()
 * Incept:   SRE 29 Feb 2000 [Seattle]; from COVE 1.0 code
 * 
 * Purpose:  Convert a secondary structure string to an array of integers
 *           representing what position each position is base-paired 
 *           to (0..len-1), or -1 if none. This is off-by-one from a
 *           Zuker .ct file representation.
 *           
 *           The structure string contains "><" for base pairs and gap symbols
 *           (usually '.') for unpaired bases.
 * 
 *           The .ct representation can accomodate pseudoknots, but the 
 *           secondary structure string cannot in general. 
 *           However, the structure string can contain "Aa", "Bb", etc. pairs as
 *           a representation of up to a 26th order pseudoknot. 
 *           (The highest order biological pseudoknot I know of is the 
 *           S1 pseudoknot, which is 3rd order: order in the Rivas and
 *           Eddy, JMB 1998 sense, referring to the depth of nesting.)

 *           Other symbols are ignored. If allow_pseudoknots is FALSE,
 *           the pseudoknot symbols will be ignored and these positions
 *           will be treated as single stranded.
 *           
 * Return:   ret_ct is allocated here and must be free'd by caller.
 *           Returns 1 on success, 0 if ss is somehow inconsistent.
 */
int 
KHS2ct(char *ss, int len, int allow_pseudoknots, int **ret_ct)
{
  Nstack_t *pda[27];                 /* 1 secondary structure + up to 26 levels of pk's */
  int      *ct;
  int       i;
  int       pos, pair;
  int       status = 1;              /* success or failure return status */

 /* Initialization: always initialize the main pda (0),
   *                 and also the 26 pk pda's if we're allowing pseudoknots.
   */
  pda[0] = CreateNstack();
  if (allow_pseudoknots) for (i = 1; i < 27; i++) pda[i] = CreateNstack();

  ct = MallocOrDie (len * sizeof(int));
  for (pos = 0; pos < len; pos++)
    ct[pos] = -1;

  for (pos = 0; pos < len; pos++)
    {
      if (!isprint(ss[pos])) status = 0; /* armor against garbage strings */

      else if (ss[pos] == '>')  /* left side of a pair: push onto stack 0 */
        PushNstack(pda[0], pos);
      else if (ss[pos] == '<')  /* right side of a pair; resolve pair */
        {
          if (! PopNstack(pda[0], &pair))
            { status = 0; }
          else
            {
              ct[pos]  = pair;
              ct[pair] = pos;
            }
        }
                                /* same stuff for pseudoknots */
      else if (isupper((int) ss[pos])) {
	if (allow_pseudoknots) PushNstack(pda[ss[pos] - 'A' + 1], pos);
      }
      else if (islower((int) ss[pos])) {
	if (allow_pseudoknots) {
          if (! PopNstack(pda[ss[pos] - 'a' + 1], &pair))
            { status = 0; }
          else
            {
              ct[pos]  = pair;
              ct[pair] = pos;
            }
	}
      }
      else if (!isgap(ss[pos])) status = 0; /* bogus character */
    }
                                /* nothing should be left on stacks */
  if (! NstackIsEmpty(pda[0])) status = 0;
  FreeNstack(pda[0]);
  if (allow_pseudoknots) {
    for (i = 1; i < 27; i++) 
      {
	if (! NstackIsEmpty(pda[i])) status = 0;
	FreeNstack(pda[i]);
      }
  }

  *ret_ct = ct;
  return status;
}


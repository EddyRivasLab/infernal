/* emit.c
 * SRE, Mon Oct 13 12:36:35 2003 [St. Louis]
 * 
 * Sampling sequences or parse trees from a CM. 
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************  
 * CVS $Id$
 */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "squid.h"
#include "sre_stack.h"

#include "structs.h"
#include "funcs.h"

/* Function:  EmitParsetree()
 * Incept:    SRE, Mon Oct 13 12:42:43 2003 [St. Louis]
 *
 * Purpose:   
 * 
 *            Shares similar structure with display.c:CreateCMConsensus().
 *
 * Args:      
 *
 * Returns:   
 *
 * Xref:      
 */


Parsetree_t *
EmitParsetree(CM_t *cm)
{
  int v;			/* counter for current state */
  int y;			/* counter for next state    */
  int type;			/* PDA_RESIDUE, PDA_STATE, or PDA_MARKER */
  Nstack_t *pda;                
  Parsetree_t *tr;

  tr  = CreateParsetree();
  pda = CreateNstack();

  PushNstack(pda, 0);
  PushNstack(pda, PDA_STATE);
  tpos = InsertTraceNode(tr, -1, TRACE_LEFT_CHILD, -1, -1, 
  while (PopNstack(pda, &type)) 
    {



      if (type == PDA_STATE) {
	PopNstack(pda, &v);


	/* OK, we've popped state v off the PDA. 
	 * Pick a production (emission from v if possible; then transition to new state y).
	 * Push productions as appropriate.
	 */
	switch (cm->sttype[v]) {
	case D_st:
	case S_st:
	  y = FChoose(cm->t[v], cm->cnum[v]);

	  if (! PushNstack(pda, y))           Die("shouldn't fail here");
	  if (! PushNstack(pda, PDA_STATE))   Die("shouldn't fail here");
	  break;

	case MP_st:
	  y = FChoose(cm->t[v], cm->cnum[v]);

	  x     = FChoose(cm->e[v], Alphabet_size*Alphabet_size);
	  lchar = Alphabet[x / Alphabet_size];
	  rchar = Alphabet[x % Alphabet_size];
	  
	  if (! PushNstack(pda, rchar))       Die("shouldn't fail here");
	  if (! PushNstack(pda, PDA_RESIDUE)) Die("shouldn't fail here");
	  if (! PushNstack(pda, y))           Die("shouldn't fail here");
	  if (! PushNstack(pda, PDA_STATE))   Die("shouldn't fail here");
	  if (! PushCstack(gsq, lchar))       Die("shouldn't fail here");	
	  break;

	case ML_st:
	case IL_st:
	  y = FChoose(cm->t[v], cm->cnum[v]);

	  x     = FChoose(cm->e[v], Alphabet_size);
	  lchar = Alphabet[x];

	  if (! PushNstack(pda, y))           Die("shouldn't fail here");
	  if (! PushNstack(pda, PDA_STATE))   Die("shouldn't fail here");
	  if (! PushCstack(gsq, lchar))       Die("shouldn't fail here");	
	  break;

	case MR_st:
	case IR_st:
	  y = FChoose(cm->t[v], cm->cnum[v]);

	  x     = FChoose(cm->e[v], Alphabet_size);
	  rchar = Alphabet[x];

	  if (! PushNstack(pda, rchar))       Die("shouldn't fail here");
	  if (! PushNstack(pda, PDA_RESIDUE)) Die("shouldn't fail here");
	  if (! PushNstack(pda, y))           Die("shouldn't fail here");
	  if (! PushNstack(pda, PDA_STATE))   Die("shouldn't fail here");

	  break;
	  
	case B_st:
	  y = cm->cfirst[v];	/* left child  */
	  z = cm->cnum[v];	/* right child */
	  
	  if (! PushNstack(pda, z))           Die("shouldn't fail here");
	  if (! PushNstack(pda, PDA_STATE))   Die("shouldn't fail here");
	  if (! PushNstack(pda, y))           Die("shouldn't fail here");
	  if (! PushNstack(pda, PDA_STATE))   Die("shouldn't fail here");
	  break;
	
	case E_st:
	  break;

	case EL_st:
	  Die("Not yet built to handle the EL state, sorry.");
	  /* NOTREACHED */
	  exit(1);

	default:
	  Die("invalid state");
	  /* NOTREACHED */
	  exit(1);
	}
  
    }
  FreeNstack(pda);

  

/* emitmap.c
 * SRE, Thu Aug  8 12:57:59 2002
 * 
 * Implementation of CMEmitMap_t structure:
 * map of a CM's nodes onto consensus sequence positions.
 * Structure is declared in structs.h.
 * 
 * Used for constructing multiple alignments.
 * 
 *   clen              : consensus length.
 *       clen is 2* n(MATP) + n(MATL) + n(MATR).
 *       The consensus sequence is indexed 1..clen.
 *       0 and clen+1 are also used, as boundary conditions.
 *                       
 *   lpos[0..nodes-1]  : has values 0 to clen+1.
 *       Any left match emission from node nd is placed in lpos[nd].
 *       Any left insert emission from node nd *follows* lpos[nd].
 *       
 *   rpos[0..nodes-1]  : has values 0..clen+1    
 *       Any right match emission from node nd is placed in rpos[nd].
 *       Any right insert emission from node nd *precedes* rpos[nd]
 *       
 *   epos[0..nodes-1]  : has values 0..clen+1.
 *       Any EL insertion from a nd->EL transition *follows* epos[nd].
 *       
 * There are no dummy values; all lpos, rpos, epos are valid coords
 * 0..clen+1, as described above, even for END_nd's.
 *
 * For nonemitting nodes, rpos and lpos give a noninclusive bound:
 * for example, lpos[0] = 0 and rpos[0] = clen+1 by definition.
 * 
 * Insertions occur between consensus positions. An inter-consensus-position
 * space may contain more than one type of insertion: an IL insertion and 
 * an EL insertion, an IR insertion and an EL insertion; or (in
 * a single absurd case of a model with a consensus length of 0) 
 * all three insertion types. Insertions are placed in order IL/EL/IR.
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************  
 * CVS $Id$
 */

#include <stdio.h>
#include <stdlib.h>

#include "squid.h"
#include "sre_stack.h"

#include "structs.h"
#include "funcs.h"

CMEmitMap_t *
CreateEmitMap(CM_t *cm)
{
  CMEmitMap_t *map;
  Nstack_t    *pda;
  int          cpos;
  int          nd;
  int          on_right;
  
  map       = MallocOrDie(sizeof(CMEmitMap_t));
  map->lpos = MallocOrDie(sizeof(int) * cm->nodes);
  map->rpos = MallocOrDie(sizeof(int) * cm->nodes);
  map->epos = MallocOrDie(sizeof(int) * cm->nodes);

  for (nd = 0; nd < cm->nodes; nd++)
    map->lpos[nd] = map->rpos[nd] = map->epos[nd] = -1;

  cpos = 0;
  nd   = 0;
  pda  = CreateNstack();
  PushNstack(pda, 0);		/* 0 = left side. 1 would = right side. */
  PushNstack(pda, nd);
  while (PopNstack(pda, &nd))
    {
      PopNstack(pda, &on_right);

      if (on_right) 
	{
	  map->rpos[nd] = cpos+1;
	  if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATR_nd) cpos++;
	}
      else
	{
	  if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd) cpos++;
	  map->lpos[nd] = cpos;

	  if (cm->ndtype[nd] == BIF_nd) 
	    {
				/* push the BIF back on for its right side  */
	      PushNstack(pda, 1);
	      PushNstack(pda, nd);
                            /* push node index for right child */
	      PushNstack(pda, 0);
	      PushNstack(pda, cm->ndidx[cm->cnum[cm->nodemap[nd]]]);   
                            /* push node index for left child */
	      PushNstack(pda, 0);
	      PushNstack(pda, cm->ndidx[cm->cfirst[cm->nodemap[nd]]]); 
	    }
	  else
	    {
				/* push the node back on for right side */
	      PushNstack(pda, 1);
	      PushNstack(pda, nd);
				/* push child node on */
	      if (cm->ndtype[nd] != END_nd) {
		PushNstack(pda, 0);
		PushNstack(pda, nd+1);
	      }
	    }
	}
    }

  /* Construct the epos map: if we do a v->EL transition,
   * the EL follows what consensus position (and its IL insertions,
   * if any)
   */
  for (nd = cm->nodes-1; nd >= 0; nd--) {
    if (cm->ndtype[nd] == END_nd) 
      cpos = map->lpos[nd];
    else if (cm->ndtype[nd] == BIF_nd) /* propagate epos for *right* branch. */
      cpos = map->epos[cm->ndidx[cm->cnum[cm->nodemap[nd]]]];

    map->epos[nd] = cpos;
  }

  map->clen = map->rpos[0]-1;
  FreeNstack(pda);
  return map;
}
  
void
DumpEmitMap(FILE *fp, CMEmitMap_t *map, CM_t *cm)
{
  int nd;

  fprintf(fp, "CM to consensus emit map; consensus length = %d \n",
	  map->clen);
  fprintf(fp, "%4s %9s %4s %4s %4s\n", 
	  "Node", "Node type", "lpos", "rpos", "epos");
  fprintf(fp, "%4s %9s %4s %4s %4s\n", 
	  "----", "---------", "----", "----", "----");
  for (nd = 0; nd < cm->nodes; nd++)
    fprintf(fp, "%4d %9s %4d %4d %4d\n", 
	    nd, Nodetype(cm->ndtype[nd]), 
	    map->lpos[nd], map->rpos[nd], map->epos[nd]);
}

void
FreeEmitMap(CMEmitMap_t *map)
{
  free(map->lpos);
  free(map->rpos);
  free(map->epos);
  free(map);
}

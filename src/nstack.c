/* nstack.c
 * SRE 1 March 2000 [Seattle]
 * An efficient pushdown stack for integers. 
 * CVS $Id$
 *
 * Basic API:
 *   say I want to push the numbers 42, 7, and 3 onto a stack,
 *   then pop them off and print them: 
 *   
 *   Nstack_t *ns;
 *   int       x;
 *
 *   ns = CreateNstack();
 *   PushNstack(ns, 42);
 *   PushNstack(ns, 7);
 *   PushNstack(ns, 3);
 *   while (PopNstack(ns, &x)) 
 *      printf("%d\n", x);
 *   FreeNstack(ns);   
 * 
 * Diagnostics:
 *   CreateNstack() returns NULL on an allocation failure.
 *   PushNstack() returns 0 on an allocation failure, else 1.
 *   PopNstack() returns 0 when the stack is empty, else 1.
 *****************************************************************
 *
 * Other functions:
 *   NstackIsEmpty(ns)      :  returns TRUE if stack is empty, else FALSE.
 *   NstackSetBlocksize(ns) :  change the chunk size for reallocation to
 *                             something other than the default 100.
 *****************************************************************
 * Implementation notes:
 *   The "stack" is kept as a growable array, ns->data. We add
 *   elements to this array in chunks, where the number of
 *   elements per chunk is set in ns->memblock. memblock starts
 *   arbitrarily at 100.
 *****************************************************************
 * @LICENSE@
 *****************************************************************
 */ 

#include <stdlib.h>
#include "nstack.h"

Nstack_t *
CreateNstack(void)
{
  Nstack_t *ns;
  
  ns           = malloc(sizeof(Nstack_t));
  if (ns == NULL) return NULL;
  ns->memblock = 100;		/* optimize if you want; hardcoded for now */
  ns->nalloc   = ns->memblock;
  ns->data     = malloc(sizeof(int) * ns->nalloc);
  if (ns->data == NULL) { free(ns); return NULL; }
  ns->n        = 0;
  return ns;
}
int
PushNstack(Nstack_t *ns, int x)
{
  int *ptr;

  if (ns->n == ns->nalloc) {
    ns->nalloc += ns->memblock;
    ptr = realloc(ns->data, sizeof(int) * ns->nalloc);
    if (ptr == NULL) return 0; else ns->data = ptr;
  }
  ns->data[ns->n] = x;
  ns->n++;
  return 1;
}
int
PopNstack(Nstack_t *ns, int *x)
{
  if (ns->n == 0) {*x = 0; return 0;}
  ns->n--;
  *x = ns->data[ns->n];
  return 1;
}
void
FreeNstack(Nstack_t *ns)
{
  free(ns->data);
  free(ns);
}
int 
NstackIsEmpty(Nstack_t *ns)
{
  if (ns->n == 0) return 1;
  else            return 0;
}
void
NstackSetBlocksize(Nstack_t *ns, int newsize)
{
  if (newsize > 0) ns->memblock = newsize;
}


Mstack_t *
CreateMstack(void)
{
  Mstack_t *ms;
  
  ms           = malloc(sizeof(Mstack_t));
  if (ms == NULL) return NULL;
  ms->memblock = 100;		/* optimize if you want; hardcoded for now */
  ms->nalloc   = ms->memblock;
  ms->data     = malloc(sizeof(void *) * ms->nalloc);
  if (ms->data == NULL) { free(ms); return NULL; }
  ms->n        = 0;
  return ms;
}
int
PushMstack(Mstack_t *ms, void *object)
{
  int *ptr;

  if (ms->n == ms->nalloc) {
    ms->nalloc += ms->memblock;
    ptr = realloc(ms->data, sizeof(void *) * ms->nalloc);
    if (ptr == NULL) return 0; else ms->data = ptr;
  }
  ms->data[ms->n] = object;
  ms->n++;
  return 1;
}
void *
PopMstack(Mstack_t *ns)
{
  if (ms->n == 0) return NULL;
  ms->n--;
  return ms->data[ms->n];
}
void
FreeMstack(Mstack_t *ms)
{
  free(ms->data);
  free(ms);
}
int 
MstackIsEmpty(Mstack_t *ms)
{
  if (ms->n == 0) return 1;
  else            return 0;
}
void
MstackSetBlocksize(Mstack_t *ns, int newsize)
{
  if (newsize > 0) ms->memblock = newsize;
}

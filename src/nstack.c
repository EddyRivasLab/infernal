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
 *****************************************************************
 *
 * Other functions:
 *   NstackIsEmpty(ns):   returns TRUE if stack is empty, else FALSE.
 *
 *****************************************************************
 * Implementation notes:
 *   The "stack" is kept as a growable array, ns->data. We add
 *   elements to this array in chunks, where the number of
 *   elements per chunk is set in ns->memblock. memblock is
 *   currently hardcoded to 100 but could be made settable at
 *   some point, and can even be dynamically adjusted in an
 *   active stack.
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

/* nstack.h
 * An efficient pushdown stack for integers.
 * SRE 1 March 2000. [Seattle]
 * CVS $Id$
 *****************************************************************
 * @LICENSE@
 *****************************************************************
 */
#ifndef NSTACKH_INCLUDED
#define NSTACKH_INCLUDED

typedef struct nstack_s {
  int *data;			/* the data stack                           */
  int  n;			/* current (topmost) elem in data           */
  int  nalloc;			/* # of elems allocated right now           */
  int  memblock;		/* memory allocation block size, # of elems */
} Nstack_t;

Nstack_t *CreateNstack(void);
int       PushNstack(Nstack_t *ns, int x);
int       PopNstack(Nstack_t *ns,  int *x);
void      FreeNstack(Nstack_t *ns);
int       NstackIsEmpty(Nstack_t *ns);

#endif NSTACKH_INCLUDED


/* nstack.h
 * Pushdown stack implementations, for integers and for objects.
 * nstack - SRE 1 March 2000. [Seattle]
 * mstack - SRE, Fri Oct 10 10:18:16 2003 [St. Louis]
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

typedef struct mstack_s {
  void **data;			/* the data stack                           */
  int    n;			/* current (topmost) elem in data           */
  int    nalloc;		/* # of elems allocated right now           */
  int    memblock;		/* memory allocation block size, # of elems */
} Mstack_t;

extern Nstack_t *CreateNstack(void);
extern int       PushNstack(Nstack_t *ns, int x);
extern int       PopNstack(Nstack_t *ns,  int *x);
extern void      FreeNstack(Nstack_t *ns);
extern int       NstackIsEmpty(Nstack_t *ns);
extern void      NstackSetBlocksize(Nstack_t *ns, int newsize);

extern Mstack_t *CreateMstack(void);
extern int       PushMstack(Mstack_t *ms, void *obj);
extern void *    PopMstack(Mstack_t *ms);
extern void      FreeMstack(Mstack_t *ms);
extern int       MstackIsEmpty(Mstack_t *ms);
extern void      MstackSetBlocksize(Mstack_t *ms, int newsize);

#endif /*NSTACKH_INCLUDED*/


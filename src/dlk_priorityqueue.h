/* dlk_priorityqueue.h
 *
 * Heap priority queue implementation for objects
 *
 *****************************************************************
 * @LICENSE@
 *****************************************************************
 */
#ifndef DLK_PQH_INCLUDED
#define DLK_PQH_INCLUDED

/* Heap operations */
extern void HeapUp(  void **data, float *priority, int top, int bottom);
extern void HeapDown(void **data, float *priority, int top, int bottom);

/* Priority queue structure and functions */
typedef struct pq_s {
  void **data;			/* the data heap                            */
  float *priority;		/* Associated priorities for the data       */
  /* NOTE: ALL operations must be done in parallel to maintain association  */
  int    n;			/* # of elems in heap                       */
  int    nalloc;		/* # of elems allocated right now           */
  int    memblock;		/* memory allocation block size, # of elems */
} PriorityQueue_t;

extern PriorityQueue_t *CreatePQ(void);
extern int              EnqueuePQ(PriorityQueue_t *pq, void *obj, float priority);
extern void            *DequeuePQ(PriorityQueue_t *pq);
extern void             FreePQ(PriorityQueue_t *pq);
extern int              IsEmptyPQ(PriorityQueue_t *pq);
extern void             SetBlocksizePQ(PriorityQueue_t *pq, int newsize);

#endif /*DLK_PQH_INCLUDED*/


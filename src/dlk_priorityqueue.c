/* dlk_priorityqueue.c
 *
 * Heap implementations of priority queue for objects
 *
 * Heap is kept as a growable array. The heap's memory is
 * grown when necessary by adding some block size. The 
 * initial allocation and the block size are set to 100
 * by default. The block size can be changed by the caller.
 * 
 *****************************************************************
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dlk_priorityqueue.h"

PriorityQueue_t *
CreatePQ(void)
{
  PriorityQueue_t *pq;
  
  pq           = malloc(sizeof(PriorityQueue_t));
  if (pq == NULL) return NULL;
  pq->memblock = 100;		/* optimize if you want; hardcoded for now */
  pq->nalloc   = pq->memblock;
  pq->data     = malloc(sizeof(void *) * pq->nalloc);
  if (pq->data == NULL) { free(pq); return NULL; }
  pq->priority = malloc(sizeof(float ) * pq->nalloc);
  if (pq->priority == NULL) { free(pq->data); free(pq); return NULL; }
  pq->n        = 0;
  return pq;
}
int
EnqueuePQ(PriorityQueue_t *pq, void *object, float priority)
{
  void **ptr1;
  float *ptr2;

  if (pq->n == pq->nalloc) {
    pq->nalloc += pq->memblock;
    ptr1 = realloc(pq->data, sizeof(void *) * pq->nalloc);
    if (ptr1 == NULL) return 0; else pq->data = ptr1;
    ptr2 = realloc(pq->priority, sizeof(float) * pq->nalloc);
    if (ptr2 == NULL) return 0; else pq->priority = ptr2;
  }

  pq->data[pq->n] = object;
  pq->priority[pq->n] = priority;
  pq->n++;
  HeapUp(pq->data,pq->priority,0,pq->n-1);
  return 1;
}
void *
DequeuePQ(PriorityQueue_t *pq)
{
  void *rtn_ptr;

  if (pq->n == 0) return NULL;
  rtn_ptr = pq->data[0];
  pq->n--;
  pq->data[0] = pq->data[pq->n];
  pq->priority[0] = pq->priority[pq->n];
  HeapDown(pq->data,pq->priority,0,pq->n-1);
  return rtn_ptr;
}
void
FreePQ(PriorityQueue_t *pq)
{
  free(pq->priority);
  free(pq->data);
  free(pq);
}
int 
IsEmptyPQ(PriorityQueue_t *pq)
{
  if (pq->n == 0) return 1;
  else            return 0;
}
void
SetBlocksizePQ(PriorityQueue_t *pq, int newsize)
{
  if (newsize > 0) pq->memblock = newsize;
}

/* Basic heap operations */
void
HeapUp(void **data, float *priority, int top, int bottom)
{
  int parent;
  void *tmp_ptr;
  float tmp;

  if (bottom > top)
  {
    parent = (bottom-1)/2;
    if (priority[parent] < priority[bottom])
    {
      tmp_ptr = data[parent]; data[parent] = data[bottom]; data[bottom] = tmp_ptr;
      tmp = priority[parent]; priority[parent] = priority[bottom]; priority[bottom] = tmp;
      HeapUp(data,priority,top,parent);
    }
  }

  return;
}

void
HeapDown(void **data, float *priority, int top, int bottom)
{
  int max_child;
  void *tmp_ptr;
  float tmp;

  if ((top*2+1) <= bottom)
  {
    max_child = top*2+1;
    if ((max_child+1) <= bottom && priority[max_child+1] > priority[max_child])
      max_child++;
    if (priority[top] < priority[max_child])
    {
      tmp_ptr = data[top]; data[top] = data[max_child]; data[max_child] = tmp_ptr;
      tmp = priority[top]; priority[top] = priority[max_child]; priority[max_child] = tmp;
      HeapDown(data,priority,max_child,bottom);
    }
  }

  return;
}

/* Should replace the following with a test of PQ */

#ifdef DLK_PQ_TESTDRIVE

/* Test driver for priority queue (heap implementation).
 * To compile:
 *    gcc -g -Wall -DDLK_PQ_TESTDRIVE -o test dlk_priorityqueue.c
 * To run:
 *    ./test
 * Returns 0 (success) w/ no output, or returns 1 and says why.
 */   
int 
main(void)
{
  PriorityQueue_t *pq;
  int      *obj;
  int       n1, n2;
  int       i;

  /* Exercises of the PQ functions/API.
   * 
   * Set the blocksize immediately to 50, so it'll allocate once
   * for 100, then next for 50 more. Putting 200 on the stack
   * therefore forces two reallocations after the inital 100 alloc.
   * 
   * Put 200 "objects" on the stack and pop them off;
   * do this twice, once with a "while pop" loop, and once
   * with a "while stack not empty" loop.
   */
  if ((pq = CreatePQ()) == NULL) {
    fprintf(stderr, "memory allocation failed\n");
    return EXIT_FAILURE;
  }
  SetBlocksizePQ(pq, 50);
  n1 = 200;
  for (i = 0; i < n1; i++)
    {
      if ((obj = malloc(sizeof(int) * 64)) == NULL) {
	fprintf(stderr, "memory allocation failed\n");
	return EXIT_FAILURE;
      }
      if (! EnqueuePQ(pq, obj, i)) {
	fprintf(stderr, "memory allocation failed in EnqueuePQ()\n");
	return EXIT_FAILURE;
      }
    }
  n2 = 0;
  while ((obj = DequeuePQ(pq)) != NULL) {
    free(obj); 
    n2++; 
  }
  if (n1 != n2){ 
    fprintf(stderr, "Put %d objects on; got %d off\n", n1, n2); 
    return EXIT_FAILURE;
  }
  for (i = 0; i < n1; i++)
    {
      if ((obj = malloc(sizeof(int) * 64)) == NULL)
	return EXIT_FAILURE;
      if (! EnqueuePQ(pq, obj, i)) {
	fprintf(stderr, "memory allocation failed in PushPriorityQueue()\n");
	return EXIT_FAILURE;
      }
    }
  n2 = 0;
  while (! IsEmptyPQ(pq)) {
    if ((obj = DequeuePQ(pq)) == NULL) {
      fprintf(stderr, "dequeue was NULL\n");
      return EXIT_FAILURE;
    }
    free(obj); 
    n2++; 
  }
  if (n1 != n2){ 
    fprintf(stderr, "Put %d objects on; got %d off\n", n1, n2); 
    return EXIT_FAILURE;
  }
  FreePQ(pq);

  return EXIT_SUCCESS;
}
#endif /*DLK_PQ_TESTDRIVE*/

/*****************************************************************  
 * @LICENSE@
 *****************************************************************/

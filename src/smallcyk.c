/* smallcyk.c
 * SRE, Wed Aug  2 08:42:49 2000 [St. Louis]
 * CVS $Id$
 * 
 * Divide and conquer CYK alignment.
 * 
 *****************************************************************
 * @LICENSE@
 *****************************************************************  
 */

#include <stdio.h>

#include "squid.h"

#include "structs.h"
#include "nstack.h"
#include "funcs.h"

static float **alloc_deck(int L);


/* Functions: deckpool_*()
 * Date:      SRE, Wed Aug  2 10:43:17 2000 [St. Louis]
 *
 * Purpose:   Implementation of a pushdown stack for storing decks
 *            of the inside dynamic programming matrix, with the
 *            usual _create, _push, _pop, and _free API.
 */
struct deckpool_s {
  float ***pool;
  int      n;
  int      nalloc;
  int      block;
};
static struct deckpool_s *
deckpool_create(void)
{
  struct deckpool_s *dpool;

  dpool = MallocOrDie(sizeof(struct deckpool_s));
  dpool->block  = 10;		/* configurable if you want */
  dpool->pool   = MallocOrDie(sizeof(float **) * dpool->block);
  dpool->nalloc = dpool->block;;
  dpool->n      = 0;
  return dpool;
}
static void 
deckpool_push(struct deckpool_s *dpool, float **deck)
{
  if (dpool->n == dpool->nalloc) {
    dpool->nalloc += dpool->block;
    dpool->pool = ReallocOrDie(dpool->pool, sizeof(float **) * dpool->nalloc);
  }
  dpool->pool[dpool->n] = deck;
  dpool->n++;
}
static float **
deckpool_pop(struct deckpool_s *d)
{
  if (d->n == 0) return NULL;
  d->n--;
  return d->pool[d->n];
}
static void
deckpool_free(struct deckpool_s *d)
{
  free(d->pool);
  free(d);
}

int
CYKDeckCount(CM_t *cm)
{
  Nstack_t *pda;
  int       v;
  int       max;
  int       ndecks;

  pda = CreateNstack();

  ndecks = 1;
  max    = 1;
  PushNstack(pda, 0);		/* push root state onto stack       */
  while (PopNstack(pda, &v)) 
    {
      ndecks--;
      for (; cm->sttype[v] != B_st && cm->sttype[v] != E_st; v++)
	{
	  if (ndecks + cm->cnum[v] + 1 > max) max = ndecks + cm->cnum[v] + 1;
	}
      if (cm->sttype[v] == B_st) 
	{
	  PushNstack(pda, cm->cnum[v]);	     /* right child S */
	  PushNstack(pda, cm->cfirst[v]);    /* left child S  */
	  ndecks+=2;
	  if (ndecks > max) max = ndecks;
	}
    }
  FreeNstack(pda);
  printf("Maximum # of decks = %d\n", max);
  return max;
}


static float **
alloc_deck(int L)
{
  float **a;
  int     j;
  
  a = MallocOrDie(sizeof(float *) * (L+1));
  for (j = 0; j <= L; j++)
    a[j] = MallocOrDie(sizeof(float) * (j+1));
  return a;
}


void
CYKInside(CM_t *cm, char *dsq, int L)
{
  float  **end;                 /* we re-use the end deck. */
  int      nends;               /* counter that tracks when we can release end deck to the pool */
  float ***alpha;               /* the scoring cube [v=0..M-1][j=0..L][d=0..j] */
  struct deckpool_s *dpool;     /* a pool of decks we can reuse */
  int     *touch;               /* keeps track of how many higher decks still need this deck */
  int      v,y,z;		/* indices for states  */
  int      j,d,i,k;		/* indices in sequence dimensions */
  float    sc;			/* a temporary variable holding a score */
  int      yoffset;		/* y=base+offset -- counter in child states that v can transit to */
  
  /* Allocations and initializations
   */
  end   = alloc_deck(L);
  nends = CMCountStatetype(cm, E_st);
  for (j = 0; j <= L; j++) {
    end[j][0] = 0.;
    for (d = 1; d <= j; d++) end[j][d] = IMPOSSIBLE;
  }
  
  alpha = MallocOrDie(sizeof(float **) * cm->M);
  for (v = 0; v < cm->M; v++) alpha[v] = NULL;

  dpool = deckpool_create();

  touch = MallocOrDie(sizeof(int) * cm->M);
  for (v = 0; v < cm->M; v++)
    touch[v] = cm->pnum[v];

  /* Main recursion
   */
  for (v = cm->M-1; v >= 0; v--) 
    {
      /* First we need a deck to fill in.
       * 1. if we're an E, reuse the end deck (and it's already calculated)
       * 2. else, see if we can take something from the pool
       * 3. else, allocate a new deck.
       */
      if (cm->sttype[v] == E_st) { 
	alpha[v] = end; continue; 
      }
      if ((alpha[v] = deckpool_pop(dpool)) == NULL)
	alpha[v] = alloc_deck(L);

      if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	{
	  for (j = 0; j <= L; j++) 
	    for (d = 0; d <= j; d++)
	      {
		y = cm->cfirst[v];
		alpha[v][j][d] = alpha[y][j][d] + cm->tsc[v][0];
		for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++) 
		  if ((sc = alpha[y+yoffset][j][d] + cm->tsc[v][yoffset]) >  alpha[v][j][d])
		    alpha[v][j][d] = sc;
		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
	      }
	}
      else if (cm->sttype[v] == B_st)
	{
	  for (j = 0; j <= L; j++) 
	    for (d = 0; d <= j; d++)
	      {
		y = cm->cfirst[v];
		z = cm->cnum[v];
		  
		alpha[v][j][d] = alpha[y][j][d] + alpha[z][j][0];
		for (k = 1; k <= d; k++)
		  if ((sc = alpha[y][j-k][d-k] + alpha[z][j][k]) > alpha[v][j][d])
		    alpha[v][j][d] = sc;
		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
	      }
	}
      else if (cm->sttype[v] == MP_st)
	{
	  for (j = 0; j <= L; j++) {
	    alpha[v][j][0] = IMPOSSIBLE;
	    if (j > 0) alpha[v][j][1] = IMPOSSIBLE;
	    for (d = 2; d <= j; d++) 
	      {
		y = cm->cfirst[v];
		alpha[v][j][d] = alpha[y][j-1][d-2] + cm->tsc[v][0];
		for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++) 
		  if ((sc = alpha[y+yoffset][j-1][d-2] + cm->tsc[v][yoffset]) >  alpha[v][j][d])
		    alpha[v][j][d] = sc;
		
		i = j-d+1;
		if (dsq[i] < Alphabet_size && dsq[j] < Alphabet_size)
		  alpha[v][j][d] += cm->esc[v][(int) (dsq[i]*Alphabet_size+dsq[j])];
		else
		  alpha[v][j][d] += DegeneratePairScore(cm->esc[v], dsq[i], dsq[j]);

		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
	      }
	  }
	}
      else if (cm->sttype[v] == IL_st || cm->sttype[v] == ML_st)
	{
	  for (j = 0; j <= L; j++) {
	    alpha[v][j][0] = IMPOSSIBLE;
	    for (d = 1; d <= j; d++)
	      {
		y = cm->cfirst[v];
		alpha[v][j][d] = alpha[y][j][d-1] + cm->tsc[v][0];
		for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++) 
		  if ((sc = alpha[y+yoffset][j][d-1] + cm->tsc[v][yoffset]) >  alpha[v][j][d])
		    alpha[v][j][d] = sc;
		
		i = j-d+1;
		if (dsq[i] < Alphabet_size)
		  alpha[v][j][d] += cm->esc[v][(int) dsq[i]];
		else
		  alpha[v][j][d] += DegenerateSingletScore(cm->esc[v], dsq[i]);
		
		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
	      }
	  }
	}
      else if (cm->sttype[v] == IR_st || cm->sttype[v] == MR_st)
	{
	  for (j = 0; j <= L; j++) {
	    alpha[v][j][0] = IMPOSSIBLE;
	    for (d = 1; d <= j; d++)
	      {
		y = cm->cfirst[v];
		alpha[v][j][d] = alpha[y][j-1][d-1] + cm->tsc[v][0];
		for (yoffset = 1; yoffset < cm->cnum[v]; yoffset++) 
		  if ((sc = alpha[y+yoffset][j-1][d-1] + cm->tsc[v][yoffset]) > alpha[v][j][d])
		    alpha[v][j][d] = sc;
		
		if (dsq[j] < Alphabet_size)
		  alpha[v][j][d] += cm->esc[v][(int) dsq[j]];
		else
		  alpha[v][j][d] += DegenerateSingletScore(cm->esc[v], dsq[j]);
		
		if (alpha[v][j][d] < IMPOSSIBLE) alpha[v][j][d] = IMPOSSIBLE;
	      }
	  }
	}				/* finished calculating deck v. */
      
      
      /* Now we worry about reuse of memory.
       * Look at our children; if they're fully released, take their deck
       * into the pool for reuse.
       */
      if (cm->sttype[v] == B_st) 
	{ /* we can definitely release the S children of a bifurc. */
	  y = cm->cfirst[v];
	  z = cm->cnum[v];
	  deckpool_push(dpool, alpha[y]);
	  deckpool_push(dpool, alpha[z]);
	  alpha[y] = NULL;
	  alpha[z] = NULL;
	}
      else
	{
	  for (y = cm->cfirst[v]; y < cm->cfirst[v]+cm->cnum[v]; y++)
	    {
	      touch[y]--;
	      if (touch[y] == 0) {
		if (cm->sttype[y] == E_st) { 
		  nends--; 
		  if (nends == 0) { deckpool_push(dpool, alpha[y]); alpha[y] = NULL; }
		} else {
		  deckpool_push(dpool, alpha[y]);
		  alpha[y] = NULL;
		}
	      }
	    }
	}
    } /* end loop over v */
  
  printf("Sir, I am proud to tell you my final score is: %.2f bits\n", 
	 alpha[0][L][L]/0.693);

    
  deckpool_free(dpool);
}

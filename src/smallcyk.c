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
static int
deckpool_pop(struct deckpool_s *d, float ***ret_deck)
{
  if (d->n == 0) { *ret_deck = NULL; return 0;}
  d->n--;
  *ret_deck = d->pool[d->n];
  return 1;
}
static void
deckpool_free(struct deckpool_s *d)
{
  free(d->pool);
  free(d);
}

void
CYKDemands(CM_t *cm, int L)
{
  float Mb_per_deck;		/* megabytes per deck */
  int   bif_decks;		/* bifurcation decks  */
  int   nends;			/* end decks (only need 1, even for multiple E's */
  int   maxdecks;		/* maximum # of decks needed by CYKInside() */
  float smallmemory;		/* how much memory small version of CYKInside() needs */
  float bigmemory;		/* how much memory a full CYKInside() would take */
  float dpcells;		/* # of dp cells */
  float bifcalcs;		/* # of inner loops executed for bifurcation calculations */
  float dpcalcs;		/* # of inner loops executed for non-bif calculations */
  int   j;

  Mb_per_deck = (float)(L+2)*(float)(L+1)*2. / 1000000.;
  bif_decks   = CMCountStatetype(cm, B_st);
  nends       = CMCountStatetype(cm, E_st);
  maxdecks    = CYKDeckCount(cm);
  smallmemory = (float) maxdecks * Mb_per_deck;
  bigmemory   = (float) (cm->M - nends +1) * Mb_per_deck;
  dpcells     = (float) (L+2)*(float)(L+1)*0.5*(float) (cm->M - nends +1);
  
  bifcalcs = 0.;
  for (j = 0; j <= L; j++)
    bifcalcs += (float)(j+1)*(float)(j+2)/2.;
  bifcalcs *= (float) bif_decks;
  
  dpcalcs = (float) (L+2)*(float)(L+1)*0.5*(float) (cm->M - bif_decks - nends +1);

  printf("CYK cpu/memory demand estimates:\n");
  printf("Mb per cyk deck:                %.2f\n", Mb_per_deck);
  printf("# of decks (M):                 %d\n",   cm->M);
  printf("# of decks needed in small CYK: %d\n",   maxdecks);
  printf("RAM needed for full CYK, Mb:    %.2f\n", bigmemory);
  printf("RAM needed for small CYK, Mb:   %.2f\n", smallmemory);
  printf("# of dp cells, total:           %.0f\n", dpcells);
  printf("# of non-bifurc dp cells:       %.0f\n", dpcalcs);
  printf("# of bifurcations:              %d\n",   bif_decks);
  printf("# of bifurc dp inner loop calcs:%.0f\n", bifcalcs);
  printf("# of dp inner loops:            %.0f\n", dpcalcs+bifcalcs);
}

int
CYKDeckCount(CM_t *cm)
{
  Nstack_t *pda;		/* pushdown stack simulating the deck pool */
  int       v,y,z;		/* state indices */
  int       nends;
  int       ndecks;
  int      *touch;		/* keeps track of how many higher decks still need this deck */

  /* Initializations, mirroring key parts of CYKInside()
   */
  ndecks = 1;			/* the end deck, which we always need. */
  nends  = CMCountStatetype(cm, E_st);
  pda    = CreateNstack();
  touch = MallocOrDie(sizeof(int) * cm->M);
  for (v = 0; v < cm->M; v++)
    touch[v] = cm->pnum[v];

  for (v = cm->M-1; v >= 0; v--)
    {
      if (cm->sttype[v] != E_st) {
	if (! PopNstack(pda, &y)) ndecks++; /* simulated allocation of a new deck */
      }
      
      if (cm->sttype[v] == B_st) { /* release both S children of a bifurc */
	y = cm->cfirst[v];
	z = cm->cnum[v];
	PushNstack(pda, y);
	PushNstack(pda, z);
      } else {
	for (y = cm->cfirst[v]; y < cm->cfirst[v]+cm->cnum[v]; y++)
	  {
	    touch[y]--;
	    if (touch[y] == 0) 
	      {
		if (cm->sttype[y] == E_st) { 
		  nends--; 
		  if (nends == 0) { PushNstack(pda, cm->M-1); }
		} else 
		  PushNstack(pda, y);
	      }
	  }
      }
    }
  FreeNstack(pda);
  free(touch);
  printf("Maximum # of decks = %d\n", ndecks);
  return ndecks;
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
static void
free_deck(float **d, int L)
{
  int j;
  for (j = 0; j <= L; j++)
    free(d[j]);
  free(d);
}
static float **
alloc_subseq_deck(int L, int i, int j)
{
  float **a;
  int     k;
  
  a = MallocOrDie(sizeof(float *) * (L+1)); /* always alloc 0..L rows, same as alloc_deck */
  for (k = 0;   k < i-1;    k++) a[k] = NULL;
  for (k = j+1; k <= L;     k++) a[k] = NULL;
  for (k = 0;   k <= j-i+1; k++) a[k] = MallocOrDie(sizeof(float) * (k+1));
  return a;
}
static void
free_subseq_deck(float **a, int i, int j)
{
  int k;
  for (k = 0; k <= j-i+1; k++) free(a[k]);
  free(a);
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
  int      ndecks;		/* total decks that had to be used in memory */
  
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
      if (! deckpool_pop(dpool, &(alpha[v]))) 
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
	      if (touch[y] == 0) 
		{
		  if (cm->sttype[y] == E_st) { 
		    nends--; 
		    if (nends == 0) { deckpool_push(dpool, end); end = NULL;}
		  } else 
		    deckpool_push(dpool, alpha[y]);
		  alpha[y] = NULL;
		}
	    }
	}
    } /* end loop over v */
  
  printf("Sir, I am proud to tell you my final score is: %.2f bits\n", 
	 alpha[0][L][L]/0.693);

  /* Now we free our memory. 
   * In principle, in alpha[], we have one root deck, and everything else is NULL;
   *               in the dpool, we have some number of decks;
   *               and end should be NULL.
   */
  free_deck(alpha[0], L);
  for (v = 1; v <= cm->M-1; v++) if (alpha[v] != NULL) Die("oi, stray deck at level %d", v);
  if (end != NULL) Die("oi, the end deck wasn't nullified");

  ndecks = 0;
  while (deckpool_pop(dpool, &end)) {
    free_deck(end, L);
    ndecks++;
  }
  deckpool_free(dpool);
}


/* Function: CYKOutside()
 * Date:     SRE, Mon Aug  7 07:45:37 2000 [St. Louis]
 */
void
CYKOutside(CM_t *cm, char *dsq, int L, float ***alpha)
{
  float ***beta;		/* the scoring cube [v=0..M-1][j=0..L][d=0..j]*/
  int      v,y,z;		/* indices for states */
  int      j,d,i,k;		/* indices in sequence dimensions */
  float    sc;			/* a temporary variable holding a score */
  struct deckpool_s *dpool;     /* a pool of decks for beta that we can reuse */
  int     *touch;               /* keeps track of how many lower decks still need this deck */
  float    escore;		/* an emission score, tmp variable */

  /* Allocations and initializations
   */
  beta = MallocOrDie(sizeof(float **) * cm->M);
  for (v = 0; v < cm->M; v++) beta[v] = NULL;

  dpool = deckpool_create();

  touch = MallocOrDie(sizeof(int) * cm->M);
  for (v = 0; v < cm->M; v++)
    if (cm->sttype[v] == B_st) touch[v] = 2;
    else                       touch[v] = cm->cnum[v];
				
  for (j = 0; j <= L; j++)
    for (d = 0; d <= j; j++)
      beta[0][j][d] = IMPOSSIBLE; /* can prob speed this initialization up */
  beta[0][L][L] = 0;		
  
  /* Main loop down through the decks
   */
  for (v = 2; v < cm->M; v++)
    {
      /* First we need to fetch a deck of memory to fill in;
       * we try to reuse a deck but if one's not available we allocate
       * a fresh one.
       */
      if (! deckpool_pop(dpool, &(beta[v])))
	beta[v] = alloc_deck(L);

      /* main recursion:
       */
      for (j = L; j >= 0; j--)
	for (d = j; d >= 0; d--) 
	  {
	    if (cm->stid[v] == BEGL_S) 
	      {
		y = cm->plast[v];	/* the parent bifurcation    */
		z = cm->cnum[y];	/* the other (right) S state */

		beta[v][j][d] = beta[y][j][d] + alpha[z][j][0]; /* init on k=0 */
		for (k = 1; k <= L-j; k++)
		  if ((sc = beta[y][j+k][d+k] + alpha[z][j+k][k]) > beta[v][j][d])
		    beta[v][j][d] = sc;
		if (beta[v][j][d] < IMPOSSIBLE) beta[v][j][d] == IMPOSSIBLE;
	      }
	    else if (cm->stid[v] == BEGR_S) 
	      {
		y = cm->plast[v];	        /* the parent bifurcation    */
		z = cm->cfirst[y];	/* the other (left) S state */

		beta[v][j][d] = beta[y][j][d] + alpha[z][j-d][0];	/* init on k=0 */
		for (k = 1; k <= j-d; k++) 
		  if ((sc = beta[y][j][d+k] + alpha[z][j-d][k]) > beta[v][j][d])
		    beta[v][j][d] = sc;
		if (beta[v][j][d] < IMPOSSIBLE) beta[v][j][d] == IMPOSSIBLE;
	      }
	    else
	      {
		alpha[v][j][d] = IMPOSSIBLE;
		i = j-d+1;
		for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--) {
		  switch(cm->sttype[j]) {
		  case MP_st: 
		    if (d == j || d == j-1) continue; /* boundary condition */

		    if (dsq[i-1] < Alphabet_size && dsq[j+1] > Alphabet_size)
		      escore = cm->esc[y][(int) (dsq[i-1]*Alphabet_size+dsq[j+1])];
		    else
		      escore = DegeneratePairScore(cm->esc[y], dsq[i-1], dsq[j+1]);

		    if ((sc = beta[y][j+1][d+2] + cm->tsc[y][v] + escore) > beta[v][j][d])
		      beta[v][j][d] = sc;
		    break;

		  case ML_st:
		  case IL_st: 
		    if (d == j) continue;	/* boundary condition (note when j=0, d=0*/

		    if (dsq[i-1] < Alphabet_size) 
		      escore = cm->esc[y][(int) dsq[i-1]];
		    else
		      escore = DegenerateSingletScore(cm->esc[y], dsq[i-1]);
		  
		    if ((sc = beta[y][j][d+1] + cm->tsc[y][v] + escore) > beta[v][j][d])
		      beta[v][j][d] = sc;
		    break;
		  
		  case MR_st:
		  case IR_st:
		    if (d == j || j == L) continue;
		  
		    if (dsq[j+1] < Alphabet_size) 
		      escore = cm->esc[y][(int) dsq[j+1]];
		    else
		      escore = DegenerateSingletScore(cm->esc[y], dsq[j+1]);

		    if ((sc = beta[y][j+1][d+1] + cm->tsc[y][v] + escore) > beta[v][j][d])
		      beta[v][j][d] = sc;
		    break;
		  
		  case B_st:
		  case E_st:
		  case D_st:
		    if ((sc = beta[y][j][d] + cm->tsc[y][v]) > beta[v][j][d])
		      beta[v][j][d] = sc;
		    break;

		  default: Die("bogus parent state %d\n", cm->sttype[y]);
		  }/* end switch over states*/
		}
	      }/*ends our handling of beta[v][j][d] */
	  }

      /* Finished deck v.
       * now worry about reuse of memory in beta:
       */
      for (y = cm->plast[v]; y > cm->plast[v]-cm->pnum[v]; y--)
	{
	  touch[y]--;
	  if (touch[y] == 0) {
	    deckpool_push(dpool, beta[y]);
	    beta[y] = NULL;
	  }
	}
    } /* end loop over decks v. */

  free(touch);
  /*dpool*/
  /*beta*/
}


/* Function: cyk_inside_engine()
 * Date:     SRE, Mon Aug  7 13:15:37 2000 [St. Louis]
 *
 * Purpose:  Run the inside phase of a CYK alignment algorithm.
 *
 * Args:     
 *
 * Returns:  
 *
 * Example:  
 */
float 
cyk_inside_engine(CM_t *cm, char *dsq, int L, int vroot, int vend, int i0, int j0, 
		  float ***ret_alpha, struct deckpool_s **ret_dpool)
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
  int      ndecks;		/* total decks that had to be used in memory */
  int      W;			/* subsequence length */
  
  /* Allocations and initializations
   */
  W     = j0-i0+1;
  end   = alloc_subseq_deck(L, i0, j0);
  nends = CMSubtreeCountStatetype(cm, vroot, E_st);
  for (j = i0-1; j <= j0; j++) {
    end[j][0] = 0.;
    for (d = 1; d <= W; d++) end[j][d] = IMPOSSIBLE;
  }
  
  alpha = MallocOrDie(sizeof(float **) * cm->M);
  for (v = 0; v < cm->M; v++) alpha[v] = NULL;

  dpool = deckpool_create();

  touch = MallocOrDie(sizeof(int) * cm->M);
  for (v = 0;     v < vroot; v++) touch[v] = 0;
  for (v = vroot; v <= vend; v++) touch[v] = cm->pnum[v];
  for (v = vend;  v < cm->M; v++) touch[v] = 0;

  /* Main recursion
   */
  for (v = vend; v >= vroot; v--) 
    {
      /* First we need a deck to fill in.
       * 1. if we're an E, reuse the end deck (and it's already calculated)
       * 2. else, see if we can take something from the pool
       * 3. else, allocate a new deck.
       */
      if (cm->sttype[v] == E_st) { 
	alpha[v] = end; continue; 
      }
      if (! deckpool_pop(dpool, &(alpha[v]))) 
	alpha[v] = alloc_subseq_deck(L, i0, j0);

      if (cm->sttype[v] == D_st || cm->sttype[v] == S_st) 
	{
	  for (j = i0-1; j <= j0; j++) 
	    for (d = 0; d <= W; d++)
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
	  for (j = i0-1; j <= j0; j++) 
	    for (d = 0; d <= W; d++)
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
	  for (j = i0-1; j <= j0; j++) {
	    alpha[v][j][0] = IMPOSSIBLE;
	    if (j > 0) alpha[v][j][1] = IMPOSSIBLE;
	    for (d = 2; d <= W; d++) 
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
	      if (touch[y] == 0) 
		{
		  if (cm->sttype[y] == E_st) { 
		    nends--; 
		    if (nends == 0) { deckpool_push(dpool, end); end = NULL;}
		  } else 
		    deckpool_push(dpool, alpha[y]);
		  alpha[y] = NULL;
		}
	    }
	}
    } /* end loop over v */
  
  printf("Sir, I am proud to tell you my final score is: %.2f bits\n", 
	 alpha[0][L][L]/0.693);

  /* Now we free our memory. 
   * In principle, in alpha[], we have one root deck, and everything else is NULL;
   *               in the dpool, we have some number of decks;
   *               and end should be NULL.
   */
  free_deck(alpha[0], L);
  for (v = 1; v <= cm->M-1; v++) if (alpha[v] != NULL) Die("oi, stray deck at level %d", v);
  if (end != NULL) Die("oi, the end deck wasn't nullified");

  ndecks = 0;
  while (deckpool_pop(dpool, &end)) {
    free_deck(end, L);
    ndecks++;
  }
  deckpool_free(dpool);
}


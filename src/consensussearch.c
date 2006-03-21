/* wordminisearch.c
 * 
 ***************************************************************** 
 * @LICENSE@
 ***************************************************************** 
 */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "squid.h"
#include "msa.h"

#include "structs.h"
#include "funcs.h"

#include "r-what.h"

static char banner[] = "consensussearch - search algorithm; no indels relative to CM consensus";

static char usage[]  = "\
Usage: consensussearch [-options] <cmfile> <sequence file>\n\
The sequence file is expected to be in FASTA format.\n\
Most commonly used options are:\n\
  -h	 : help\n\
  -W <n> : set scanning window size to <n> (default: 200)\n\
  -S <n> : set score threshold for reporting hits to <n> (default: 10)\n\
  -T <n> : set word score threshold to <n> (default: -INF, no threshold)\n\
  -X <n> : set dropoff score for extension to <n> (default -10)\n\
  --wordlen <n> : set word length (in number of states) to <n> (default 4)\n\
";

static struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE },
  { "-W", TRUE, sqdARG_INT },
  { "-X", TRUE, sqdARG_FLOAT },
  { "-S", TRUE, sqdARG_FLOAT },
  { "-T", TRUE, sqdARG_FLOAT },
  { "--wordlen", FALSE, sqdARG_INT },
  { "--toponly", FALSE, sqdARG_NONE },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

/* function declarations for this file */
float BCMiniCYKOut(CM_t *cm, int *consensus_d, char *dsq, int L, float dropoff_sc, BPA_t **master_align);
float BCMiniCYKIn(CM_t *cm, int *consensus_d, char *dsq, int L, float dropoff_sc, BPA_t *master_align);


/* Function: BCMiniCYKOut()
 * Author:   DLK
 *
 * Purpose:  Given a seed, extend a CM hit in the Outside direction
 *           Extension is gapless (relative to consensus model)
 *
 * Args:     cm		- the covariance model
 *           consensus_d- array of consensus subsequence lengths (subscript v states)
 *           dsq	- digitized sequence to search; 1..L
 *           L		- length of the sequence
 *           dropoff_sc - equiv. to BLAST parameter X: continue extension
 *                        until score drops from max by X or more
 *           master_align - branched partial alignment structure.  Pointer
 *                        to structure must be modifiable, as we may add
 *                        parent chunks as we extend outward.  On call,
 *                        should be allocated and contain initial alignment
 *                        information.
 *
 * Returns:  master_align
 *           return value: score over the extension (seed score not included)
 */
float
BCMiniCYKOut(CM_t *cm, int *consensus_d, char *dsq, int L, float dropoff_sc, BPA_t **master_align)
{
  int i,j,d,v;
  int y;
  float tsc, esc;
  float max_sc = 0.0;
  float delta_sc = 0.0;
  BPA_t *temp_align = NULL;

  /* Get initial info from master_align */
  v = (*master_align)->chunk->init_v;
  j = (*master_align)->chunk->init_j;
  d = (*master_align)->chunk->init_d;

  while (delta_sc > dropoff_sc) { /* Continue extension */
    /* Get new v and calculate tsc */
    y = cm->plast[v] - cm->pnum[v] + 1;
    tsc = cm->tsc[y][v-cm->cfirst[y]];
    v = y;
    
    /* Get new j,d and Calculate esc */
    if ((cm->stid[v] == MATP_MP) && (j<(L-1)) && ((j-d+1)>0)) 
    {
      j++; d+=2;
      i = j-d+1;
      if (dsq[i] < Alphabet_size && dsq[j] < Alphabet_size)
	esc = cm->esc[v][(int) (dsq[i]*Alphabet_size+dsq[j])];
      else
	esc = DegeneratePairScore(cm->esc[v],dsq[i],dsq[j]);
    }
    else if ((cm->stid[v] == MATL_ML) && ((j-d+1)>0))
    {
      d++;
      i = j-d+1;
      if (dsq[i] < Alphabet_size)
	esc = cm->esc[v][(int) dsq[i]];
      else
	esc = DegenerateSingletScore(cm->esc[v],dsq[i]);
    }
    else if ((cm->stid[v] == MATR_MR) && (j<(L-1)))
    {
      j++; d++;
      if (dsq[j] < Alphabet_size)
	esc = cm->esc[v][(int) dsq[j]];
      else
	esc = DegenerateSingletScore(cm->esc[v],dsq[j]);
    }
    else if (cm->stid[v] == BEGL_S)
    {
      if (cm->pnum[v] != 1) fprintf(stderr,"WARNING: BEGL_S has more than one parent?\n");
      y = cm->plast[v];
      if (cm->stid[y] != BIF_B) fprintf(stderr,"WARNING: BEGL_S parent not BIF_B?\n");
      temp_align = malloc(sizeof(BPA_t));
      temp_align->chunk = malloc(sizeof(PA_t));
      temp_align->chunk->init_v = temp_align->chunk->cur_v = y;
      temp_align->chunk->init_d = temp_align->chunk->cur_d = consensus_d[y];
      temp_align->chunk->init_j = temp_align->chunk->cur_j = j + consensus_d[cm->cnum[y]];
      temp_align->left_child = *master_align;
      temp_align->right_child = malloc(sizeof(BPA_t));
      temp_align->right_child->chunk = malloc(sizeof(PA_t));
      temp_align->right_child->chunk->init_v = temp_align->right_child->chunk->cur_v = cm->cnum[y];
      temp_align->right_child->chunk->init_d = temp_align->right_child->chunk->cur_d = consensus_d[cm->cnum[y]];
      temp_align->right_child->chunk->init_j = temp_align->right_child->chunk->cur_j = j + consensus_d[cm->cnum[y]];

      esc = BCMiniCYKIn(cm,consensus_d,dsq,L,dropoff_sc,temp_align->right_child);
      v = temp_align->chunk->cur_v;
      j = temp_align->chunk->cur_j;
      d = temp_align->chunk->cur_d;
    }
    else if (cm->stid[v] == BEGR_S)
    {
      if (cm->pnum[v] != 1) fprintf(stderr,"WARNING: BEGR_S has more than one parent?\n");
      y = cm->plast[v];
      if (cm->stid[y] != BIF_B) fprintf(stderr,"WARNING: BEGR_S parent not BIF_B?\n");
      temp_align = malloc(sizeof(BPA_t));
      temp_align->chunk = malloc(sizeof(PA_t));
      temp_align->chunk->init_v = temp_align->chunk->cur_v = y;
      temp_align->chunk->init_d = temp_align->chunk->cur_d = consensus_d[y];
      temp_align->chunk->init_j = temp_align->chunk->cur_j = j + consensus_d[cm->cfirst[y]];
      temp_align->right_child = *master_align;
      temp_align->left_child = malloc(sizeof(BPA_t));
      temp_align->left_child->chunk = malloc(sizeof(PA_t));
      temp_align->left_child->chunk->init_v = temp_align->left_child->chunk->cur_v = cm->cfirst[y];
      temp_align->left_child->chunk->init_d = temp_align->left_child->chunk->cur_d = consensus_d[cm->cfirst[y]];
      temp_align->left_child->chunk->init_j = temp_align->left_child->chunk->cur_j = j + consensus_d[cm->cfirst[y]];

      esc = BCMiniCYKIn(cm,consensus_d,dsq,L,dropoff_sc,temp_align->right_child);
      v = temp_align->chunk->cur_v;
      j = temp_align->chunk->cur_j;
      d = temp_align->chunk->cur_d;
    }
    else 	/* Nonpermitted indel state *
		 * Extension terminates     */
    {
      break;
    }

    if (d != consensus_d[v])
      fprintf(stderr,"WARNING: departed from consensus d in BCMiniCYKOut !\n");

    /* Check for improvement over max_sc */
    delta_sc = delta_sc + tsc + esc;
    if (delta_sc > 0)
    {
      max_sc += delta_sc;
      delta_sc = 0.0;
      if (temp_align == NULL)
      {
        (*master_align)->chunk->init_v = v;
        (*master_align)->chunk->init_j = j;
        (*master_align)->chunk->init_d = d;
      }
      else
      {
        *master_align = temp_align;
        (*master_align)->chunk->init_v = v;
        (*master_align)->chunk->init_j = j;
        (*master_align)->chunk->init_d = d;
        temp_align = NULL;
/* CHECK! do the last v,j,d before we hit the bifurcation ever get written to the appropriate child ? */
      }
    }
  }

  if (temp_align != NULL) BPA_Free(temp_align);

  return max_sc;
}

/* Function: BCMiniCYKIn()
 * Author:   DLK
 *
 * Purpose:  Given a seed, extend a CM hit in the Inside direction
 *           Extension is gapless (relative to consensus model)
 *
 * Args:     cm		- the covariance model
 *           consensus_d- array of consensus subsequence lengths (subscript v states)
 *           dsq	- digitized sequence to search; 1..L
 *           L		- length of the sequence
 *           dropoff_sc - equiv. to BLAST parameter X: continue extension
 *                        until score drops from max by X or more
 *           master_align - branched partial alignment structure.  For Inside
 *                        extension, only structure (not reference to it) needs
 *                        to be modified
 *
 * Returns:  master_align
 *           return value: score over the extension (seed score not included)
 */
float
BCMiniCYKIn(CM_t *cm, int *consensus_d, char *dsq, int L, float dropoff_sc, BPA_t *master_align)
{
  int i,j,d,v;
  int y, yoffset, z;
  float tsc, esc;
  BPA_t *left_child;
  BPA_t *right_child;
  float max_sc = 0.0;
  float delta_sc = 0.0;

  int bif_flag = 0;

  /* Get initial info from master_align */
  v = master_align->chunk->cur_v;
  j = master_align->chunk->cur_j;
  d = master_align->chunk->cur_d;

  if (d != consensus_d[v]) fprintf(stderr,"WARNING: d doesn't match consensus_d in BCMiniCYKIn\n");

  while ((delta_sc > dropoff_sc) && (!bif_flag)){ /* Continue extension */
    /* Get new v and calculate tsc */
    y = cm->cfirst[v];
    yoffset = 0;
    z = cm->stid[y+yoffset];
    while (yoffset+1 < cm->cnum[v] && z != MATP_MP && z != MATL_ML && z != MATR_MR)
    {
      yoffset++;
      z = cm->stid[y+yoffset];
    }
    tsc = cm->tsc[v][yoffset];
    v = y+yoffset;
    
    /* Get new j,d and Calculate esc */
    if (cm->stid[v] == MATP_MP) 
    {
      j--;
      d-=2;
      i = j-d+1;
      if (i >= j) { break; }
      if (dsq[i] < Alphabet_size && dsq[j] < Alphabet_size)
	esc = cm->esc[v][(int) (dsq[i]*Alphabet_size+dsq[j])];
      else
	esc = DegeneratePairScore(cm->esc[v],dsq[i],dsq[j]);
    }
    else if (cm->stid[v] == MATL_ML)
    {
      d--;
      i = j-d+1;
      if (i >= j) { break; }
      if (dsq[i] < Alphabet_size)
	esc = cm->esc[v][(int) dsq[i]];
      else
	esc = DegenerateSingletScore(cm->esc[v],dsq[i]);
    }
    else if (cm->stid[v] == MATR_MR)
    {
      j--; d--;
      i = j-d+1;
      if (i >= j) { break; }
      if (dsq[j] < Alphabet_size)
	esc = cm->esc[v][(int) dsq[j]];
      else
	esc = DegenerateSingletScore(cm->esc[v],dsq[j]);
    }
    else if (cm->stid[v] == BIF_B)
    {
    bif_flag = 1;
    }
    else 	/* Nonpermitted indel state
		 * Extension terminates */
    {
      break;
    }

    if (d != consensus_d[v]) fprintf(stderr,"WARNING: d departed from consensus_d in BCMiniCYKIn\n");

    /* Check for improvement over max_sc */
    delta_sc = delta_sc + tsc + esc;
    if (delta_sc > 0)
    {
      max_sc += delta_sc;
      delta_sc = 0.0;
      master_align->chunk->cur_v = v;
      master_align->chunk->cur_j = j;
      master_align->chunk->cur_d = d;
    }
  }

  if (bif_flag)
  {
    y = cm->cfirst[v];
    z = cm->cnum[v];

    left_child = malloc(sizeof(BPA_t));
    left_child->chunk = malloc(sizeof(PA_t));
    left_child->chunk->init_v = y;
    left_child->chunk->init_j = j - consensus_d[z];
    left_child->chunk->init_d = consensus_d[y];
    left_child->left_child = NULL;
    left_child->right_child = NULL;

    right_child = malloc(sizeof(BPA_t));
    right_child->chunk = malloc(sizeof(PA_t));
    right_child->chunk->init_v = y;
    right_child->chunk->init_j = j - consensus_d[z];
    right_child->chunk->init_d = consensus_d[y];
    right_child->left_child = NULL;
    right_child->right_child = NULL;

    delta_sc += BCMiniCYKIn(cm,consensus_d,dsq,L,dropoff_sc,left_child);
    delta_sc += BCMiniCYKIn(cm,consensus_d,dsq,L,dropoff_sc,right_child);

    if (delta_sc > 0)
    {
      max_sc += delta_sc;
      delta_sc = 0.0;
      master_align->chunk->cur_v = v;
      master_align->chunk->cur_j = j;
      master_align->chunk->cur_d = d;
      master_align->left_child = left_child;
      master_align->right_child = right_child;
    }
    else
    {
      BPA_Free(left_child);
      BPA_Free(right_child);
    }
  }

  return max_sc;
}

/* Function: WordInitiate()
 * Author:   DLK
 *
 * Purpose:  Given a v,j,d position, create an initial word hit
 *           based at that location.  Derived from MiniCYKIn.
 *
 * Args:     cm		- the covariance model
 *           dsq	- the digitized sequence
 *           L		- length of the sequence
 *           v		- starting state
 *           j		- starting j position
 *           d		- d = j-i+1 - relative starting i position
 *           wordlen	- target length of word hit
 *           iscored	- TRUE/FALSE: has residue i been emitted
 *           jscored	- TRUE/FALSE: has residue j been emitted
 *           ret_v	- RETURN: ending state
 *           ret_j	- RETURN: ending j
 *           ret_d	- RETURN: ending d
 *
 * Returns:  iscored, jscored, ret_v, ret_j, ret_d
 *           return value: score of word hit
 */
float
WordInitiate(CM_t *cm, char *dsq, int L, int v, int j, int d, int wordlen,
             int *iscored, int *jscored, int *ret_v, int *ret_j, int *ret_d)
{
  int i;
  int x;
  int y, yoffset, z;
  float tsc, esc, seed_sc;

  /* Keep track of whether initial i,j have been emitted yet */
  *iscored = FALSE; *jscored = FALSE;

  *ret_v = v;
  *ret_j = j;
  *ret_d = d;
  seed_sc = 0;
  
  /* esc for starting state */
  if (cm->stid[v] == MATP_MP)
  {
    i = j - d + 1;
    if (dsq[i] < Alphabet_size && dsq[j] < Alphabet_size)
      seed_sc += cm->esc[v][(int) (dsq[i]*Alphabet_size+dsq[j])];
    else
      seed_sc += DegeneratePairScore(cm->esc[v],dsq[i],dsq[j]);
    *iscored = TRUE;
    *jscored = TRUE;
  }
  else if (cm->stid[v] == MATL_ML)
  {
    i = j - d + 1;
    if (dsq[i] < Alphabet_size)
      seed_sc += cm->esc[v][(int) dsq[i]];
    else
      seed_sc += DegenerateSingletScore(cm->esc[v],dsq[i]);
    *iscored = TRUE;
  }
  else if (cm->stid[v] == MATR_MR)
  {
    i = j - d + 1;
    if (dsq[j] < Alphabet_size)
      seed_sc += cm->esc[v][(int) dsq[j]];
    else
      seed_sc += DegenerateSingletScore(cm->esc[v],dsq[j]);
    *jscored = TRUE;
  }
  else /* Should never reach here */
  { fprintf(stderr,"WARNING: Illegal state on word initiation v=%d stid=%d\n",v,cm->stid[v]); }

  /* Add tsc+esc for additional states */
  for (x=1; x<wordlen; x++)	/* Here, x is counter for states currently in word */
  {
    y = cm->cfirst[v];
    yoffset = 0;
    z = cm->stid[y+yoffset];
    while (yoffset+1 < cm->cnum[v] && z != MATP_MP && z != MATL_ML && z != MATR_MR)
    {
      yoffset++;
      z = cm->stid[y+yoffset];
    }
    tsc = cm->tsc[v][yoffset];
    v = y+yoffset;

    if (cm->stid[v] == MATP_MP)
    {
      if (*iscored && *jscored) { j--; d-=2; }
      else if (*iscored) { d--; *jscored = TRUE; }
      else if (*jscored) { j--; d--; *iscored = TRUE; }
      i = j - d + 1;
      if (i >= j) { break; }
      if (dsq[i] < Alphabet_size && dsq[j] < Alphabet_size)
        esc = cm->esc[v][(int) (dsq[i]*Alphabet_size+dsq[j])];
      else
	esc = DegeneratePairScore(cm->esc[v],dsq[i],dsq[j]);
    }
    else if (cm->stid[v] == MATL_ML)
    {
      if (*iscored) { d--; }
      else         { *iscored = TRUE; }
      i = j - d + 1;
      if (i >= j) { break; }
      if (dsq[i] > Alphabet_size)
	esc = cm->esc[v][(int) dsq[i]];
      else
	esc = DegenerateSingletScore(cm->esc[v],dsq[i]);
    }
    else if (cm->stid[v] == MATR_MR)
    {
      if (*jscored) { j--; d--; }
      else         { *jscored = TRUE; }
      i = j - d + 1;
      if (i >= j) { break; }
      if (dsq[j] > Alphabet_size)
	esc = cm->esc[v][(int) dsq[j]];
      else
	esc = DegenerateSingletScore(cm->esc[v],dsq[j]);
    }
    else	/* Nonpermitted state - terminate */
      break;

    seed_sc += tsc + esc;
    *ret_v = v;
    *ret_j = j;
    *ret_d = d;
  }

  /* Force wordlength - given impossible score if not */
  if ( x<wordlen ) { seed_sc = IMPOSSIBLE; }
    
  return seed_sc;
}

/* Function: ConsensusScan()
 * Author:   DLK
 *
 * Purpose:  Prototype BLAST-like algorithm (brute force)
 *           For all (i,j) and all v that are MP, attempt
 *           ungapped extension; report hits
 *
 */
void
ConsensusScan(CM_t *cm, char *dsq, int L, int W, int wordlen, float word_sc,
         float dropoff_sc, float report_sc)
{
  int v, inner_v, outer_v;	/* state indices, 0..M-1 */
  int i, inner_i, outer_i;	/* left sequence indices */
  int j, inner_j, outer_j;	/* right sequence indices */
  int d;			/* subsequence length, d = j-1+1 */
  float seed_sc,in_sc, out_sc, total_sc;
  				/* score components */
  int    nhits = 0;		/* # of hits */
  int   *hit_outer_v;		/* Outer state index of hits */
  int   *hit_inner_v;		/* Inner state index of hits */
  int   *hit_outer_i;		/* Outer left seq. positions */
  int   *hit_inner_i;		/* Inner left seq. positions */
  int   *hit_outer_j;		/* Outer right seq. positions */
  int   *hit_inner_j;		/* Inner right seq. positions */
  float *hit_sc;		/* Scores of hits */
  int   *hit_index;		/* index of sorted hits */
  int    alloc_nhits = 1000;	/* Amount of space allocated */
  int    x;			/* counter variable */
  int    duplicate;		/* TRUE/FALSE: hit is already in list */

  float  max_sc;
  int    max_i;
  int    tmp,y;

  int   *consensus_d;

  int iscored, jscored;		/* TRUE/FALSE: has residue i/j been emitted yet */
  int temp_v, temp_i, temp_j, temp_d;

  int    duplicate_low_index;	/* Start looking for duplicates here; below this is impossible */

  hit_outer_v = MallocOrDie(sizeof(int)   * alloc_nhits);
  hit_inner_v = MallocOrDie(sizeof(int)   * alloc_nhits);
  hit_outer_i = MallocOrDie(sizeof(int)   * alloc_nhits);
  hit_inner_i = MallocOrDie(sizeof(int)   * alloc_nhits);
  hit_outer_j = MallocOrDie(sizeof(int)   * alloc_nhits);
  hit_inner_j = MallocOrDie(sizeof(int)   * alloc_nhits);
  hit_sc      = MallocOrDie(sizeof(float) * alloc_nhits);

  consensus_d = MallocOrDie(sizeof(int)   * cm->M);

  duplicate_low_index = 0;
  for (j=1; j<=L; j++)
  {
    for (v=cm->M-1; v>0; v--)
    {
      if (cm->stid[v] == MATP_MP || cm->stid[v] == MATL_ML || cm->stid[v] == MATR_MR)
      {
        d = consensus_d[v];
	seed_sc = WordInitiate(cm,dsq,L,v,j,d,wordlen,&iscored,&jscored,&temp_v,&temp_j,&temp_d);

	if (seed_sc > word_sc)
	{
	  /* If the seed score meets the threshold, extend both out and in */
	  /* NOTE: unscored i/j residues are always assigned to the inner subsequence */
	  /* MiniCYKIn() deals with them if it can, while MiniCYKOut() ignores them   */
	  in_sc  = MiniCYKIn( cm,dsq,L,temp_v,temp_j,temp_d,dropoff_sc,&iscored,&jscored,&inner_v,&inner_i,&inner_j);
          out_sc = MiniCYKOut(cm,dsq,L,     v,     j,     d,dropoff_sc,                  &outer_v,&outer_i,&outer_j);
          total_sc = seed_sc + in_sc + out_sc;

          if (!iscored && outer_i==j-d+1)	/* Entire hit is single-stranded, never emits left residue */
          {
            inner_i = 0;
            outer_i = 0;
          }
          else if (!iscored)			/* outer_i to i-1 were emitted by MiniCYKOut, i not emitted */
          {
            inner_i--;
          }
          else if (!jscored && outer_j==j)	/* Entire hit is single-stranded, never emits right residue */
          {
            inner_j = 0;
            outer_j = 0;
          }
          else if (!jscored)			/* j+1 to outer_j were emitted by MiniCYKOut, j not emitted */
          {
            inner_j++;
          }

          if (total_sc > report_sc)	/* Add to hitlist */
          {
            /* If extended score meets reporting threshold, check for duplicates */
            duplicate = FALSE;
            x = duplicate_low_index;
            while (x<nhits && !duplicate)
            {
              if (outer_i == hit_outer_i[x])
                if (inner_i == hit_inner_i[x])
                  if (outer_j == hit_outer_j[x])
                    if (inner_j == hit_inner_j[x])
                      if (outer_v == hit_outer_v[x])
                        if (inner_v == hit_inner_v[x])
                          duplicate = TRUE;
              if (j > hit_outer_j[x] + W)
                duplicate_low_index = x+1;
              x++;
            }

            if (! duplicate)
            {
              /* If not a duplicate, record the hit */
              hit_outer_v[nhits] = outer_v;
              hit_inner_v[nhits] = inner_v;
              hit_outer_i[nhits] = outer_i;
              hit_inner_i[nhits] = inner_i;
              hit_outer_j[nhits] = outer_j;
              hit_inner_j[nhits] = inner_j;
              hit_sc[nhits]      = total_sc;
              nhits++;

              if (nhits == alloc_nhits)
              {
                hit_outer_v = ReallocOrDie(hit_outer_v, sizeof(int)   * (alloc_nhits + 1000));
                hit_inner_v = ReallocOrDie(hit_inner_v, sizeof(int)   * (alloc_nhits + 1000));
                hit_outer_i = ReallocOrDie(hit_outer_i, sizeof(int)   * (alloc_nhits + 1000));
                hit_inner_i = ReallocOrDie(hit_inner_i, sizeof(int)   * (alloc_nhits + 1000));
                hit_outer_j = ReallocOrDie(hit_outer_j, sizeof(int)   * (alloc_nhits + 1000));
                hit_inner_j = ReallocOrDie(hit_inner_j, sizeof(int)   * (alloc_nhits + 1000));
                hit_sc      = ReallocOrDie(hit_sc,      sizeof(float) * (alloc_nhits + 1000));
                alloc_nhits += 1000;
              }
            }
          }
        }
      }
    }
  }
  printf("outer v\tinner v\touter i\tinner i\t");
  printf("inner j\touter j\ttotal sc\n");

  hit_index = MallocOrDie(sizeof(int) * nhits);
  for (x=0; x<nhits; x++) { hit_index[x] = x; }
  i=0;
  while (i<nhits)
  {
    /* Select top scoring hit and print it */
    max_i  = i;
    max_sc = hit_sc[hit_index[i]];
    for (j=i+1; j<nhits; j++)
    {
      if (hit_sc[hit_index[j]] > max_sc)
      {
	max_sc = hit_sc[hit_index[j]];
	max_i  = j;
      }
    }
    tmp = hit_index[i];
    hit_index[i] = hit_index[max_i];
    hit_index[max_i] = tmp;

    x = hit_index[i];
    printf("%d\t%d\t",      hit_outer_v[x],hit_inner_v[x]);
    printf("%d\t%d\t",      hit_outer_i[x],hit_inner_i[x]);
    printf("%d\t%d\t%.3f\n",hit_inner_j[x],hit_outer_j[x],hit_sc[x]);

    /* Eliminate lower-scoring hits that overlap */
    for (j=i+1; j<nhits; j++)
    {
      y = hit_index[j];
      /* Lazy overlap checking! 
       * Should also test for a segment of the lower hit completely surrounding
       * a segment of the higher hit
       */
      if ( (hit_outer_i[y] >= hit_outer_i[x] && hit_outer_i[y] <= hit_inner_i[x] ) ||
           (hit_inner_i[y] >= hit_outer_i[x] && hit_inner_i[y] <= hit_inner_i[x] ) ||
           (hit_inner_j[y] >= hit_outer_i[x] && hit_inner_j[y] <= hit_inner_i[x] ) ||
           (hit_outer_j[y] >= hit_outer_i[x] && hit_outer_j[y] <= hit_inner_i[x] ) ||
           (hit_outer_i[y] >= hit_inner_j[x] && hit_outer_i[y] <= hit_outer_j[x] ) ||
           (hit_inner_i[y] >= hit_inner_j[x] && hit_inner_i[y] <= hit_outer_j[x] ) ||
           (hit_inner_j[y] >= hit_inner_j[x] && hit_inner_j[y] <= hit_outer_j[x] ) ||
           (hit_outer_j[y] >= hit_inner_j[x] && hit_outer_j[y] <= hit_outer_j[x] )  )
      {
	tmp = hit_index[j];
	hit_index[j] = hit_index[nhits-1];
	hit_index[nhits-1] = tmp;
	nhits--;
	j--;	/*Back up and check the one we just swapped in */
      }
    }
    i++;
  }
  
  free(hit_outer_v);
  free(hit_inner_v);
  free(hit_outer_i);
  free(hit_inner_i);
  free(hit_outer_j);
  free(hit_inner_j);
  free(hit_sc);
  free(hit_index);

  free(consensus_d);

  return;
}

int
main(int argc, char **argv)
{
  char			*cmfile;	/* file to read CM from */
  char			*seqfile;	/* file to read sequences from */
  int			 format;	/* format of sequence file */
  CMFILE		*cmfp;		/* CM file for reading */
  SQFILE		*sqfp;		/* seqfile for reading */
  CM_t			*cm;		/* covariance model */
  char			*seq;		/* RNA sequence */
  SQINFO		 sqinfo;	/* optional info attached to seq */
  char			*dsq;		/* digitized RNA sequence */
  int			 reversed;	/* TRUE when we're doing the reverse complement strand */

  int	windowlen;			/* scanning window size */
  float dropoff_sc;			/* Dropoff for termination of hit extension */
  float report_sc;			/* Threshold for reporting hits */
  float word_sc;			/* Threshold for words to initiate extension */
  int   wordlen;			/* Length of word in # of states.  # of nucs may vary */
  int	do_revcomp;			/* TRUE to do reverse complement too */

  char	*optname;			/* name of option found by Getopt() */
  char	*optarg;			/* argument found by Getopt() */
  int	 optind;			/* index in argv[] */

  /* float **msc; */

  /* Parse command line */
  format = SQFILE_UNKNOWN;
  windowlen = 200;
  dropoff_sc = -10;
  report_sc = 0;
  word_sc = IMPOSSIBLE;
  wordlen = 4;
  do_revcomp = TRUE;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
	        &optind, &optname, &optarg)) {
    if      (strcmp(optname, "-W")        == 0) windowlen       = atoi(optarg);
    else if (strcmp(optname, "-X")        == 0) dropoff_sc      = atof(optarg);
    else if (strcmp(optname, "-S")        == 0) report_sc       = atof(optarg);
    else if (strcmp(optname, "-T")        == 0) word_sc         = atof(optarg);
    else if (strcmp(optname, "--wordlen") == 0) wordlen         = atof(optarg);
    else if (strcmp(optname, "--toponly") == 0) do_revcomp      = FALSE;
    else if (strcmp(optname, "-h")        == 0) {
      MainBanner(stdout, banner);
      puts(usage);
      exit(EXIT_SUCCESS);
    }
  }


  if (argc - optind != 2) Die("Incorrect number of arguments.\n%s\n", usage);
  cmfile = argv[optind++];
  seqfile = argv[optind++];

  /* File input, get a CM */

  if ((sqfp = SeqfileOpen(seqfile, format, NULL)) == NULL)
    Die("Failed to open sequence database file %s\n%s\n", seqfile, usage);
  if ((cmfp = CMFileOpen(cmfile, NULL)) == NULL)
    Die("Failed to open covariance model save file %s\n%s\n", cmfile, usage);

  if (! CMFileRead(cmfp, &cm))
    Die("Failed to read a CM from %s -- file corrupt?\n", cmfile);
  if (cm == NULL)
    Die("Failed to read a CM from %s -- file empty?\n", cmfile);


  CMLogoddsify(cm);
  CMHackInsertScores(cm);	/* make insert emissions score zero -- "TEMPORARY" FIX */

  reversed = FALSE;
  while (reversed || ReadSeq(sqfp, sqfp->format, &seq, &sqinfo))
  {
    if (sqinfo.len == 0) continue;	/* silently skip len 0 seqs */
    dsq = DigitizeSequence(seq, sqinfo.len);

    if (! reversed) printf("sequence: %s\n", sqinfo.name);
    else            printf("reversed sequence: %s %d\n", sqinfo.name, sqinfo.len);
    ConsensusScan(cm, dsq, sqinfo.len, windowlen, wordlen, word_sc, dropoff_sc, report_sc);

    free(dsq);
    if (! reversed && do_revcomp)
    {
      revcomp(seq,seq);
      reversed = TRUE;
    }
    else
    {
      reversed = FALSE;
      FreeSequence(seq, &sqinfo);
    }
  }

  FreeCM(cm);
  CMFileClose(cmfp);
  SeqfileClose(sqfp);
  SqdClean();
  return EXIT_SUCCESS;
}

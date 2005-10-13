/* minicyk.c
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

static char banner[] = "minisearch - prototype BLAST-like search algorithm";

static char usage[]  = "\
Usage: minisearch [-options] <cmfile> <sequence file>\n\
The sequence file is expected to be in FASTA format.\n\
Most commonly used options are:\n\
  -h	 : help\n\
  -W <n> : set scanning window size to <n> (default: 200)\n\
  -S <n> : set score threshold for reporting hits to <n> (default: 10)\n\
  -X <n> : set dropoff score for extension to <n> (default -10)\n\
";

static struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE },
  { "-W", TRUE, sqdARG_INT },
  { "-X", TRUE, sqdARG_FLOAT },
  { "-S", TRUE, sqdARG_FLOAT },
  { "--toponly", FALSE, sqdARG_NONE },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))


/* Function: MiniCYKOut()
 * Author:   DLK
 *
 * Purpose:  Given a seed, extend a CM hit in the Outside direction
 *           Extension is gapless (relative to consensus model)
 *           and unbifurcated.
 *
 * Args:     cm		- the covariance model
 *           dsq	- digitized sequence to search; 1..L
 *           L		- length of the sequence
 *           v		- starting state for extension
 *           j		- starting j position for extension
 *           d		- d = j-i+1 - relative starting i position
 *           dropoff_sc - equiv. to BLAST parameter X: continue extension
 *                        until score drops from max by X or more
 *           ret_out_v  - RETURN: ending state
 *           ret_out_i  - RETURN: ending i
 *           ret_out_j  - RETURN: ending j
 *
 * Returns:  ret_out_v, ret_out_i, ret_out_j allocated and freed by caller
 *           return value: score over the extension (seed score not included)
 */
float
MiniCYKOut(CM_t *cm, char *dsq, int L, int v, int j, int d, float dropoff_sc,
           int *ret_out_v, int *ret_out_i, int *ret_out_j)
{
  int i;
  int y;
  float tsc, esc;
  float max_sc = 0.0;
  float delta_sc = 0.0;

  *ret_out_v = v;
  *ret_out_i = j-d+1;
  *ret_out_j = j;

  /* Given: memory allocations, v,j,d */
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
    else 	/* Nonpermitted indel or bifucation state
		 * Extension terminates */
    {
      break;
    }

    /* Check for improvement over max_sc */
    delta_sc = delta_sc + tsc + esc;
    if (delta_sc > 0)
    {
      max_sc += delta_sc;
      delta_sc = 0.0;
      *ret_out_v = v;
      *ret_out_i = j-d+1;
      *ret_out_j = j;
    }
  }

  return max_sc;
}

/* Function: MiniCYKIn()
 * Author:   DLK
 *
 * Purpose:  Given a seed, extend a CM hit in the Inside direction
 *           Extension is gapless (relative to consensus model)
 *           and unbifurcated.
 *
 * Args:     cm		- the covariance model
 *           dsq	- digitized sequence to search; 1..L
 *           L		- length of the sequence
 *           v		- starting state for extension
 *           j		- starting j position for extension
 *           d		- d = j-i+1 - relative starting i position
 *           dropoff_sc - equiv. to BLAST parameter X: continue extension
 *                        until score drops from max by X or more
 *           ret_in_v   - RETURN: ending state
 *           ret_in_i   - RETURN: ending i
 *           ret_in_j   - RETURN: ending j
 *
 * Returns:  ret_in_v, ret_in_i, ret_in_j allocated and freed by caller
 *           return value: score over the extension (seed score not included)
 */
float
MiniCYKIn(CM_t *cm, char *dsq, int L, int v, int j, int d, float dropoff_sc,
          int *ret_in_v, int *ret_in_i, int *ret_in_j)
{
  int i;
  int y, yoffset, z;
  float tsc, esc;
  float max_sc = 0.0;
  float delta_sc = 0.0;

  *ret_in_v = v;
  *ret_in_i = j-d+1;
  *ret_in_j = j;

  /* Given: memory allocations, v,j,d */
  while (delta_sc > dropoff_sc) { /* Continue extension */
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
      j--; d-=2;
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
    else 	/* Nonpermitted indel or bifucation state
		 * Extension terminates */
    {
      break;
    }

    /* Check for improvement over max_sc */
    delta_sc = delta_sc + tsc + esc;
    if (delta_sc > 0)
    {
      max_sc += delta_sc;
      delta_sc = 0.0;
      *ret_in_v = v;
      *ret_in_i = j-d+1;
      *ret_in_j = j;
    }
  }

  return max_sc;
}

/* Function: MiniScan()
 * Author:   DLK
 *
 * Purpose:  Prototype BLAST-like algorithm (brute force)
 *           For all (i,j) and all v that are MP, attempt
 *           ungapped extension; report hits
 *
 */
void
MiniScan(CM_t *cm, char *dsq, int L, int W, float dropoff_sc, float report_sc
        )
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
  float *hit_sc;			/* Scores of hits */
  int    alloc_nhits = 10;	/* Amount of space allocated */
  int    x;			/* counter variable */
  int    duplicate;		/* TRUE/FALSE: hit is already in list */

  hit_outer_v = MallocOrDie(sizeof(int)   * alloc_nhits);
  hit_inner_v = MallocOrDie(sizeof(int)   * alloc_nhits);
  hit_outer_i = MallocOrDie(sizeof(int)   * alloc_nhits);
  hit_inner_i = MallocOrDie(sizeof(int)   * alloc_nhits);
  hit_outer_j = MallocOrDie(sizeof(int)   * alloc_nhits);
  hit_inner_j = MallocOrDie(sizeof(int)   * alloc_nhits);
  hit_sc      = MallocOrDie(sizeof(float) * alloc_nhits);

  for (j=1; j<=L; j++)
  {
    for (v=cm->M-1; v>0; v--)
    {
      if (cm->sttype[v] == MP_st)
      {
        for (d=1; d<=W && d<=j; d++)
        {
	  i=j-d+1;
	  if (dsq[i] < Alphabet_size && dsq[j] < Alphabet_size)
	    seed_sc = cm->esc[v][(int) (dsq[i]*Alphabet_size+dsq[j])];
	  else
	    seed_sc = DegeneratePairScore(cm->esc[v],dsq[i],dsq[j]);
	  in_sc  = MiniCYKIn( cm,dsq,L,v,j,d,dropoff_sc,&inner_v,&inner_i,&inner_j);
          out_sc = MiniCYKOut(cm,dsq,L,v,j,d,dropoff_sc,&outer_v,&outer_i,&outer_j);
	  total_sc = seed_sc + in_sc + out_sc;

	  if (total_sc > report_sc)	/* Add to hitlist */
	  {
	    duplicate = FALSE;
	    x = 0;
	    while (x<nhits && !duplicate)
	    {
	      if (outer_v == hit_outer_v[x])
		if (inner_v == hit_inner_v[x])
		  if (outer_i == hit_outer_i[x])
		    if (inner_i == hit_inner_i[x])
		      if (outer_j == hit_outer_j[x])
			if (inner_j == hit_inner_j[x])
			  duplicate = TRUE;
	      x++;
	    }

            if (! duplicate)
	    {
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
                hit_outer_v = ReallocOrDie(hit_outer_v, sizeof(int)   * (alloc_nhits + 10));
                hit_inner_v = ReallocOrDie(hit_inner_v, sizeof(int)   * (alloc_nhits + 10));
                hit_outer_i = ReallocOrDie(hit_outer_i, sizeof(int)   * (alloc_nhits + 10));
                hit_inner_i = ReallocOrDie(hit_inner_i, sizeof(int)   * (alloc_nhits + 10));
                hit_outer_j = ReallocOrDie(hit_outer_j, sizeof(int)   * (alloc_nhits + 10));
                hit_inner_j = ReallocOrDie(hit_inner_j, sizeof(int)   * (alloc_nhits + 10));
                hit_sc      = ReallocOrDie(hit_sc,      sizeof(float) * (alloc_nhits + 10));
		alloc_nhits += 10;
	      }
	    }
	  }
	}
      }
    }
  }
  printf("outer v\touter i\touter j\t");
  printf("inner v\tinner i\tinner j\ttotal sc\n");
  for (x = 0; x<nhits; x++)
  {
    printf("%d\t%d\t%d\t",hit_outer_v[x],hit_outer_i[x],hit_outer_j[x]);
    printf("%d\t%d\t%d\t%.3f\n",hit_inner_v[x],hit_inner_i[x],hit_inner_j[x],hit_sc[x]);
  }

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
  int	do_revcomp;			/* TRUE to do reverse complement too */

  char	*optname;			/* name of option found by Getopt() */
  char	*optarg;			/* argument found by Getopt() */
  int	 optind;			/* index in argv[] */

  /* Parse command line */
  format = SQFILE_UNKNOWN;
  windowlen = 200;
  dropoff_sc = -10;
  report_sc = 10;
  do_revcomp = TRUE;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
	        &optind, &optname, &optarg)) {
    if      (strcmp(optname, "-W")        == 0) windowlen       = atoi(optarg);
    else if (strcmp(optname, "-X")        == 0) dropoff_sc      = atof(optarg);
    else if (strcmp(optname, "-S")        == 0) report_sc       = atof(optarg);
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
    else            printf("reversed sequence: %s\n", sqinfo.name);
    MiniScan(cm, dsq, sqinfo.len, windowlen, dropoff_sc, report_sc);

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

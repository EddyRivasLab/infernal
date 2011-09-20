/* The "brute" integration test. Based on hmmer's test of the
 * same name.
 * 
 * Create an entirely hand-specified covariance model from given
 * parameters; enumerate all paths and calculate either the sum or max
 * in either standard or truncated mode (Inside, TrInside, CYK, or
 * TrCYK) and compare it to results from our implementations of those
 * algorithms (cm_InsideAlign(), cm_TrInsideAlign(), cm_CYKAlign(), 
 * cm_TrCYKAlign()).
 * 
 * The test currently only exercises a simple four node model composed
 * of one ROOT, two MATL and one End node.
 * 
 * Besides the hand-specified model, the integration test samples many
 * more "brute" CMs randomly, peppering them with zero probability
 * transitions where possible.
 * 
 * CYK scores (hand enumerated vs. cm_CYKAlign() or cm_TrCYKAlign())
 * should match exactly (within machine precision). Inside scores
 * should match "closely", with some error introduced by the
 * discretization in FLogsum()'s table lookup.
 *
 * This is an important test of correctness for the non-banded CYK and
 * Inside implementations. HMM banded versions can then be tested 
 * against those, with the knowledge that the bands can affect the
 * optimal score.
 * 
 * EPN, Tue Sep 20 04:37:17 2011
 */

/*  gcc   -o itest_brute     -std=gnu99 -g -Wall -I. -L. -I../hmmer/src -L../hmmer/src -I../easel -L../easel itest_brute.c -linfernal -lhmmer -leasel -lm
 */
#include "esl_config.h"
#include "p7_config.h"
#include "config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_vectorops.h"

#include "funcs.h"		/* function declarations                */
#include "structs.h"		/* data structures, macros, #define's   */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-o",        eslARG_OUTFILE, NULL, NULL, NULL,  NULL,  NULL, NULL, "save each tested CM to file <f>",                0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose",                                     0 },
  { "-N",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "number of randomly sampled CMs",                 0 },
  { "--vv",      eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be very verbose",                                0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <cmfile>";
static char banner[] = "the brute force CM integration test";

struct cm_brute_matl_param_s {
  /* Node 1: ROOT node 
   *                    State 0: ROOT_S
   *                    State 1: ROOT_IL
   *                    State 2: ROOT_IR
   * Node 2: MATL node  
   *                    State 3: MATL_ML
   *                    State 4: MATL_D
   *                    State 5: MATL_IL
   * Node 3: MATL node  
   *                    State 6: MATL_ML
   *                    State 7: MATL_D
   *                    State 8: MATL_IL (detached insert, unreachable to avoid ambiguity) 
   * Node 4: END node  
   *                    State 9: END_E
   */
  double t0t1;      	/* cm->t[0][0] ROOT_S  (0) -> ROOT_IL (1) */
  double t0t2;      	/* cm->t[0][1] ROOT_S  (0) -> ROOT_IR (2) */
  double t0t3;      	/* cm->t[0][2] ROOT_S  (0) -> MATL_ML (3) */
  double t0t4;      	/* cm->t[0][3] ROOT_S  (0) -> MATL_D  (4) */

  double t1t1;      	/* cm->t[1][0] ROOT_IL (1) -> ROOT_IL (1) */
  double t1t2;      	/* cm->t[1][1] ROOT_IL (1) -> ROOT_IR (2) */
  double t1t3;      	/* cm->t[1][2] ROOT_IL (1) -> MATL_ML (3) */
  double t1t4;      	/* cm->t[1][3] ROOT_IL (1) -> MATL_D  (4) */

  double t2t2;      	/* cm->t[2][0] ROOT_IR (2) -> ROOT_IR (2) */
  double t2t3;      	/* cm->t[2][1] ROOT_IR (2) -> MATL_ML (3) */
  double t2t4;      	/* cm->t[2][2] ROOT_IR (2) -> MATL_D  (4) */

  double t3t5;      	/* cm->t[3][0] MATL_ML (3) -> MATL_IL (5) */
  double t3t6;      	/* cm->t[3][1] MATL_ML (3) -> MATL_ML (6) */
  double t3t7;      	/* cm->t[3][2] MATL_ML (3) -> MATL_D  (7) */

  double t4t5;      	/* cm->t[4][0] MATL_D  (4) -> MATL_IL (5) */
  double t4t6;      	/* cm->t[4][1] MATL_D  (4) -> MATL_ML (6) */
  double t4t7;      	/* cm->t[4][2] MATL_D  (4) -> MATL_D  (7) */

  double t5t5;      	/* cm->t[5][0] MATL_IL (5) -> MATL_IL (5) */
  double t5t6;      	/* cm->t[5][1] MATL_IL (5) -> MATL_ML (6) */
  double t5t7;      	/* cm->t[5][2] MATL_IL (5) -> MATL_D  (7) */

  /* remainder of transitions are fixed, due to MATL_IL (8) being a detached insert, to remove ambiguity */
  /* cm->t[6][0] MATL_ML (6) -> MATL_IL (8)     will be 0.0, b/c state 8 is detached */
  /* cm->t[6][1] MATL_ML (6) -> END_E   (9)     will be 0.0, b/c state 8 is detached */

  /* cm->t[7][0] MATL_D  (7) -> MATL_IL (8)     will be 0.0, b/c state 8 is detached */
  /* cm->t[7][1] MATL_D  (7) -> END_E   (9)     will be 0.0, b/c state 8 is detached */

  /* cm->t[8][0] MATL_IL (8) -> MATL_IL (8)     is irrelevant, b/c state 8 is detached */
  /* cm->t[8][1] MATL_IL (8) -> END_E   (9)     is irrelevant, b/c state 8 is detached */

  double alpha;  	/* cm->e[v][A] emission for both match states (MATL_ML (v=3) and MATL_ML (v=6)) */
  double beta;  	/* cm->e[v][A] emission for all insert states (v = 1, 2, 5, 8) */

  double begin[4];	/* constructed from transitions when brute profile is configured. */
  double end;		/* internal ends, set when profile is configured */
};

static void        set_brute_matl_params(struct cm_brute_matl_param_s *prm);
static void        sample_zeropeppered_probvector(ESL_RANDOMNESS *r, double *p, int n);
static void        sample_brute_matl_params(ESL_RANDOMNESS *r, struct cm_brute_matl_param_s *prm);
static CM_t       *create_brute_matl_cm(ESL_ALPHABET *abc, char *errbuf, struct cm_brute_matl_param_s *prm);
static double      score_brute_matl_cm(struct cm_brute_matl_param_s *prm, double nullA, int do_cyk, double sc[3]);

int
main(int argc, char **argv)
{
  struct cm_brute_matl_param_s prm;
  int             status;
  ESL_GETOPTS    *go       = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_ALPHABET   *abc      = esl_alphabet_Create(eslRNA);
  ESL_RANDOMNESS *r        = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  CM_t           *cm       = NULL;
  char           *cmfile   = esl_opt_GetString (go, "-o");
  int             N        = esl_opt_GetInteger(go, "-N");
  int             do_local;
  double          brute_ins[5];	/* lod Inside scores for seqs L=0..4 calculated by brute force path enumeration */
  double          brute_cyk[5];	/* lod CYK    scores for seqs L=0..4 calculated by brute force path enumeration */
  float           ins_sc[5];	/* lod scores for seqs L=0..4 calculated by cm_InsideAlign() DP */
  float           cyk_sc[5];	/* lod scores for seqs L=0..4 calculated by cm_CYKAlign() DP */
  ESL_DSQ         dsq[6];
  int             L;
  int             i,j;
  float           cyk_precision, ins_precision; /* expected bound on absolute accuracy for CYK, Inside */
  char            errbuf[cmERRBUFSIZE];
  CM_MX          *mx   = NULL;       /* alpha DP matrix for non-banded CYK/Inside() */
  CM_SHADOW_MX   *shmx = NULL;       /* shadow matrix for non-banded tracebacks */

  /* setup logsum lookups */
  FLogsumInit();

  for (do_local = 0; do_local <= 0; do_local++) /* run tests in both glocal and local mode   */
    for (j = 0; j <= N; j++)	                /* #0 = fixed params; #1..N = sampled params */
      {
	if (esl_opt_GetBoolean(go, "-v")) 
	  printf("%s\n", do_local ? "Local mode" : "Glocal mode");
	  
	if (j == 0)      set_brute_matl_params(&prm);
	else             sample_brute_matl_params(r, &prm);

	cm     = create_brute_matl_cm(abc, errbuf, &prm);
	mx     = cm_mx_Create(cm);
	shmx   = cm_shadow_mx_Create(cm);

	score_brute_matl_cm(&prm, cm->null[0], TRUE,  brute_cyk);
	score_brute_matl_cm(&prm, cm->null[0], FALSE, brute_ins);
  
	if (cmfile)
	  {
	    FILE *ofp = fopen(cmfile, "w");
	    cm_file_WriteASCII(ofp, -1, cm);
	    fclose(ofp);
	  }

	for (L = 0; L <= 2; L++)
	  {
	    dsq[0] = dsq[L+1] = eslDSQ_SENTINEL;       /* Initialize dsq of length L at 0000... (all A) */
	    for (i = 1; i <= L; i++) dsq[i] = 0;
	    
	    if ((status = cm_CYKAlign    (cm, errbuf, dsq, L, 0, cm->M-1, 1, L, FALSE, 128., shmx, NULL, NULL, mx, &(cyk_sc[L])))  != eslOK) esl_fatal("CYK failed: %s", errbuf);
	    if (esl_opt_GetBoolean(go, "--vv")) cm_mx_Dump(stdout, mx);

	    if ((status = cm_InsideAlign(cm, errbuf, dsq, 1, L, 128., mx, &(ins_sc[L]))) != eslOK)  esl_fatal("Inside failed: %s", errbuf);
	    if (esl_opt_GetBoolean(go, "--vv")) cm_mx_Dump(stdout, mx);

	    cyk_precision = 1e-4;    /* default impl uses fp, should be accurate within machine precision      */
	    ins_precision = 0.015;   /* default impl uses FLogsum, tolerate 2^0.015 ~= 1% error in Forward probs */

	    if (esl_opt_GetBoolean(go, "-v")) 
	      printf("%d %-6s %6s %1d %8.4f %8.4f %8.4f %8.4f\n",
		     j,
		     do_local ? "local" : "glocal",
		     (j > 0)  ? "random": "fixed",
		     L, 
		     brute_ins[L], ins_sc[L], 
		     brute_cyk[L], cyk_sc[L]);

	    if (fabs(cyk_sc[L] - brute_cyk[L]) > cyk_precision)
	      esl_fatal("CYK    scores mismatched: %-6s %s  L=%1d brute=%8.4f cm_CYKAlign()=%8.4f (difference %g)",
			do_local ? "local" : "glocal",
			(j > 0)  ? "random": "fixed",
			L, 
			brute_cyk[L], cyk_sc[L], fabs(brute_cyk[L] - cyk_sc[L]));

	    /* verify that Inside scores match closely (within error introduced by FLogsum() */
	    if (fabs(ins_sc[L] - brute_ins[L]) > ins_precision) 
	      esl_fatal("Inside scores mismatched: %-6s %s L=%1d brute=%8.4f cm_InsideAlign()=%8.4f",
			do_local ? "local" : "glocal",
			(j > 0)  ? "random": "fixed",
			L, 
			brute_ins[L], ins_sc[L]);
	  }
	FreeCM(cm);
	cm_mx_Destroy(mx);
	cm_shadow_mx_Destroy(shmx);
      }

  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);

  printf("ok\n");
  return 0;
}

static void
set_brute_matl_params(struct cm_brute_matl_param_s *prm)
{
  prm->t0t1 = 0.15;      	/* cm->t[0][0] ROOT_S  (0) -> ROOT_IL (1) */
  prm->t0t2 = 0.05;      	/* cm->t[0][1] ROOT_S  (0) -> ROOT_IR (2) */
  prm->t0t3 = 0.70;      	/* cm->t[0][2] ROOT_S  (0) -> MATL_ML (3) */
  prm->t0t4 = 0.10;      	/* cm->t[0][3] ROOT_S  (0) -> MATL_D  (4) */

  prm->t1t1 = 0.10;      	/* cm->t[1][0] ROOT_IL (1) -> ROOT_IL (1) */
  prm->t1t2 = 0.10;      	/* cm->t[1][1] ROOT_IL (1) -> ROOT_IR (2) */
  prm->t1t3 = 0.60;      	/* cm->t[1][2] ROOT_IL (1) -> MATL_ML (3) */
  prm->t1t4 = 0.20;      	/* cm->t[1][3] ROOT_IL (1) -> MATL_D  (4) */

  prm->t2t2 = 0.25;      	/* cm->t[2][0] ROOT_IR (2) -> ROOT_IR (2) */
  prm->t2t3 = 0.70;      	/* cm->t[2][1] ROOT_IR (2) -> MATL_ML (3) */
  prm->t2t4 = 0.05;      	/* cm->t[2][2] ROOT_IR (2) -> MATL_D  (4) */

  prm->t3t5 = 0.08;      	/* cm->t[3][0] MATL_ML (3) -> MATL_IL (5) */
  prm->t3t6 = 0.80;      	/* cm->t[3][1] MATL_ML (3) -> MATL_ML (6) */
  prm->t3t7 = 0.12;      	/* cm->t[3][2] MATL_ML (3) -> MATL_D  (7) */

  prm->t4t5 = 0.03;      	/* cm->t[4][0] MATL_D  (4) -> MATL_IL (5) */
  prm->t4t6 = 0.89;      	/* cm->t[4][1] MATL_D  (4) -> MATL_ML (6) */
  prm->t4t7 = 0.08;      	/* cm->t[4][2] MATL_D  (4) -> MATL_D  (7) */

  prm->t5t5 = 0.07;      	/* cm->t[5][0] MATL_IL (5) -> MATL_IL (5) */
  prm->t5t6 = 0.82;      	/* cm->t[5][1] MATL_IL (5) -> MATL_ML (6) */
  prm->t5t7 = 0.11;      	/* cm->t[5][2] MATL_IL (5) -> MATL_D  (7) */

  /* remainder of transitions are fixed, due to detached insert */
  /* cm->t[6][0] MATL_ML (6) -> MATL_IL (8)     will be 0.0, b/c state 8 is detached */
  /* cm->t[6][1] MATL_ML (6) -> END_E   (9)     will be 1.0, b/c state 8 is detached */

  /* cm->t[7][0] MATL_D  (7) -> MATL_IL (8)     will be 0.0, b/c state 8 is detached */
  /* cm->t[7][1] MATL_D  (7) -> END_E   (9)     will be 1.0, b/c state 8 is detached */ 

  /* cm->t[8][0] MATL_IL (8) -> MATL_IL (8)     will be 0.0, b/c state 8 is detached */
  /* cm->t[8][1] MATL_IL (8) -> END_E   (9)     will be 1.0, b/c state 8 is detached */

  prm->alpha = 0.7;  	/* cm->e[v][A] emission for both match states (MATL_ML (v=3) and MATL_ML (v=6)) */
  prm->beta  = 0.25;  	/* must be 0.25, insert score 0, cm->e[v][A] emission for all insert states (v = 1, 2, 5, 8) */

  return;
}


static void
sample_zeropeppered_probvector(ESL_RANDOMNESS *r, double *p, int n)
{
  esl_dirichlet_DSampleUniform(r, n, p);
  if (esl_rnd_Roll(r, 2))	/* coin flip */
    {
      p[esl_rnd_Roll(r, n)] = 0.0;
      esl_vec_DNorm(p, n);
    }
}


static void
sample_brute_matl_params(ESL_RANDOMNESS *r, struct cm_brute_matl_param_s *prm)
{
  double tmp[4];

  /* make sure we can get into main match or S state of next node w/ nonzero prob, but pepper zeros elsewhere 
   * also make sure zero length path through all deletes is non-zero prob, (prm->t0t4 != 0 && prm-t4t7 != 0) */
  do { sample_zeropeppered_probvector(r, tmp, 4);  prm->t0t1 = tmp[0]; prm->t0t2 = tmp[1]; prm->t0t3 = tmp[2]; prm->t0t4 = tmp[3]; } while (prm->t0t1 == 0.0 || prm->t0t4 == 0.0);
  do { sample_zeropeppered_probvector(r, tmp, 4);  prm->t1t1 = tmp[0]; prm->t1t2 = tmp[1]; prm->t1t3 = tmp[2]; prm->t1t4 = tmp[3]; } while (prm->t1t3 == 0.0);
  do { sample_zeropeppered_probvector(r, tmp, 3);  prm->t2t2 = tmp[0]; prm->t2t3 = tmp[1]; prm->t2t4 = tmp[2];                     } while (prm->t2t3 == 0.0);
  do { sample_zeropeppered_probvector(r, tmp, 3);  prm->t3t5 = tmp[0]; prm->t3t6 = tmp[1]; prm->t3t7 = tmp[2];                     } while (prm->t3t6 == 0.0);
  do { sample_zeropeppered_probvector(r, tmp, 3);  prm->t4t5 = tmp[0]; prm->t4t6 = tmp[1]; prm->t4t7 = tmp[2];                     } while (prm->t4t6 == 0.0 || prm->t4t7 == 0.0);
  do { sample_zeropeppered_probvector(r, tmp, 3);  prm->t5t5 = tmp[0]; prm->t5t6 = tmp[1]; prm->t5t7 = tmp[2];                     } while (prm->t5t6 == 0.0);
  /* transition probabilities are fixed for the remainder of the states, b/c they have a detached insert as a child */

  /* make sure x=A emissions for match, insert are nonzero */
  prm->alpha = esl_rnd_UniformPositive(r);
  prm->beta  = esl_rnd_UniformPositive(r);	/* doesn't have to match background, although by default inserts score 0 */

  return;
}

static CM_t *
create_brute_matl_cm(ESL_ALPHABET *abc, char *errbuf, struct cm_brute_matl_param_s *prm)
{
  int           status;
  CM_t         *cm     = NULL;
  int           nnodes = 4;  /* ROOT, MATL, MATL, END */
  int           M      = 10; /* ROOT_S, ROOT_IL, ROOT_IR, MATL_ML, MATL_D, MATL_IL, MATL_ML, MATL_D, MATL_IL, END_E */
  int           clen   = 2;
  int           v, i, j;
  Parsetree_t  *gtr;		/* guide tree for alignment                   */
  CMConsensus_t *cons = NULL;
  if (abc->type != eslRNA) esl_fatal("brute CM uses RNA alphabet");

  cm = CreateCM(nnodes, M, clen, abc);
  CMZero(cm);

  gtr = CreateParsetree(4); /* the guide tree we'll construct the CM from */
  i = 1;
  j = 2;
  v = -1;
  /* add ROOT_nd */
  v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, ROOT_nd);
  /* add MATL_nd */
  v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, MATL_nd);
  i++;
  /* add MATL_nd */
  v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, MATL_nd);
  i++;
  /* add END_nd  */
  v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, END_nd);

  if((status = cm_from_guide(cm, errbuf, gtr, FALSE)) != eslOK) cm_Fail("Failed to create CM: %s", errbuf);
  FreeParsetree(gtr);

  cm->t[0][0] = prm->t0t1;      	/* cm->t[0][0] ROOT_S  (0) -> ROOT_IL (1) */
  cm->t[0][1] = prm->t0t2;      	/* cm->t[0][1] ROOT_S  (0) -> ROOT_IR (2) */
  cm->t[0][2] = prm->t0t3;      	/* cm->t[0][2] ROOT_S  (0) -> MATL_ML (3) */
  cm->t[0][3] = prm->t0t4;      	/* cm->t[0][3] ROOT_S  (0) -> MATL_D  (4) */

  cm->t[1][0] = prm->t1t1;      	/* cm->t[1][0] ROOT_IL (1) -> ROOT_IL (1) */
  cm->t[1][1] = prm->t1t2;      	/* cm->t[1][1] ROOT_IL (1) -> ROOT_IR (2) */
  cm->t[1][2] = prm->t1t3;      	/* cm->t[1][2] ROOT_IL (1) -> MATL_ML (3) */
  cm->t[1][3] = prm->t1t4;      	/* cm->t[1][3] ROOT_IL (1) -> MATL_D  (4) */

  cm->t[2][0] = prm->t2t2;      	/* cm->t[2][0] ROOT_IR (2) -> ROOT_IR (2) */
  cm->t[2][1] = prm->t2t3;      	/* cm->t[2][1] ROOT_IR (2) -> MATL_ML (3) */
  cm->t[2][2] = prm->t2t4;      	/* cm->t[2][2] ROOT_IR (2) -> MATL_D  (4) */

  cm->t[3][0] = prm->t3t5;      	/* cm->t[3][0] MATL_ML (3) -> MATL_IL (5) */
  cm->t[3][1] = prm->t3t6;      	/* cm->t[3][1] MATL_ML (3) -> MATL_ML (6) */
  cm->t[3][2] = prm->t3t7;      	/* cm->t[3][2] MATL_ML (3) -> MATL_D  (7) */

  cm->t[4][0] = prm->t4t5;      	/* cm->t[4][0] MATL_D  (4) -> MATL_IL (5) */
  cm->t[4][1] = prm->t4t6;      	/* cm->t[4][1] MATL_D  (4) -> MATL_ML (6) */
  cm->t[4][2] = prm->t4t7;      	/* cm->t[4][2] MATL_D  (4) -> MATL_D  (7) */

  cm->t[5][0] = prm->t5t5;      	/* cm->t[5][0] MATL_IL (5) -> MATL_IL (5) */
  cm->t[5][1] = prm->t5t6;      	/* cm->t[5][1] MATL_IL (5) -> MATL_ML (6) */
  cm->t[5][2] = prm->t5t7;      	/* cm->t[5][2] MATL_IL (5) -> MATL_D  (7) */

  cm->t[6][0] = 0.0;             	/* cm->t[6][0] MATL_ML (6) -> MATL_IL (8)     will be 0.0, b/c state 8 is detached */
  cm->t[6][1] = 1.0;            	/* cm->t[6][1] MATL_ML (6) -> END_E   (9) */

  cm->t[7][0] = 0.0;             	/* cm->t[7][0] MATL_D  (7) -> MATL_IL (8)     will be 0.0, b/c state 8 is detached */
  cm->t[7][1] = 1.0;            	/* cm->t[7][1] MATL_D  (7) -> END_E   (9) */

  /* cm->t[8] is irrelevant and not reset, because that state cannot be reached */

  for(v = 0; v < M; v++) { 
    if(cm->sttype[v] == ML_st) { 
      esl_vec_FSet(cm->e[v], abc->K, (1.0-prm->alpha)/(float)(abc->K-1));
      cm->e[v][0] = prm->alpha;
    }
    if(cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) { 
      esl_vec_FSet(cm->e[v], abc->K, (1.0-prm->beta)/(float)(abc->K-1));
      cm->e[v][0] = prm->beta;
    }
  }
  
  CMRenormalize(cm);
  CMLogoddsify(cm);

  /* Add mandatory annotation */
  cm_SetName(cm, "itest-brute");
  cm->W = cm->clen;
  CreateCMConsensus(cm, cm->abc, 3.0, 1.0, &(cons));
  if ((status = cm_SetConsensus(cm, cons, NULL)) != eslOK) cm_Fail("Failed to calculate consensus sequence");
  cm->nseq     = 1;
  cm->eff_nseq = 1;
  cm->checksum = 0;
  

  if(cons != NULL) FreeCMConsensus(cons);

  return cm;
}


/* score_brute_matl_cm() enumerates all paths combinatorially, and
 * calculates their Inside or CYK probabilities either by summing
 * or by max, for A* (polyA) sequences of lengths 0..2.
 */
static double
score_brute_matl_cm(struct cm_brute_matl_param_s *prm, double nullA, int do_cyk, double sc[3])
{
  double msc = prm->alpha / nullA;
  double isc = prm->beta  / nullA;
  double pL[3];			   /* summed odds of all paths of length 0..2  */
  double cp[19];		   /* odds of 19 possible paths through the CM */


  /* 1. There are 19 possible paths that up to L=2 residues can align
     to the core model. 
  */

  /* 0 residues */
  cp[0]  = prm->t0t4 * prm->t4t7 * 1.0; 	                 /* ROOT_S(0) MATL_D (4) MATL_D (7) END_E  (9)           (L=0) */ 

  /* 1 residue */
  cp[1]  = isc * prm->t0t1 * prm->t1t4 * prm->t4t7 * 1.0;        /* ROOT_S(0) ROOT_IL(1) MATL_D (4) MATL_D (7) END_E (9) (L=1) */
  cp[2]  = isc * prm->t0t2 * prm->t2t4 * prm->t4t7 * 1.0;        /* ROOT_S(0) ROOT_IR(2) MATL_D (4) MATL_D (7) END_E (9) (L=1) */
  cp[3]  = msc * prm->t0t3 * prm->t3t7 * 1.0;                    /* ROOT_S(0) MATL_ML(3) MATL_D (7) END_E  (9)           (L=1) */
  cp[4]  = msc * prm->t0t4 * prm->t4t6 * 1.0;                    /* ROOT_S(0) MATL_D (4) MATL_ML(6) END_E  (9)           (L=1) */
  cp[5]  = isc * prm->t0t4 * prm->t4t5 * prm->t5t7 * 1.0;        /* ROOT_S(0) MATL_D (4) MATL_IL(5) MATL_D (7) END_E (9) (L=1) */
  /* skipped because MATL_IL(8) is detached:                        ROOT_S(0) MATL_D (4) MATL_D (7) MATL_IL(8) END_E (9) (L=1) */

  /* 2 residues */
  /* one possibility with 2 ML emits */
  cp[6]  = msc * msc * prm->t0t3 * prm->t3t6 * 1.0;              /* ROOT_S(0) MATL_ML(3) MATL_ML(6) END_E  (9)           (L=2) */

  /* six possibilities with 1 ML emit, 1 insert emit (and two others that are only impossible because MATL_IL(8) is detached */
  cp[7]  = isc * msc * prm->t0t1 * prm->t1t3 * prm->t3t7 * 1.0;  /* ROOT_S(0) ROOT_IL(1) MATL_ML(3) MATL_D (7) END_E (9) (L=2) */
  cp[8]  = isc * msc * prm->t0t1 * prm->t1t4 * prm->t4t6 * 1.0;  /* ROOT_S(0) ROOT_IL(1) MATL_D (4) MATL_ML(6) END_E (9) (L=2) */
  cp[9]  = isc * msc * prm->t0t2 * prm->t2t3 * prm->t3t7 * 1.0;  /* ROOT_S(0) ROOT_IR(2) MATL_ML(3) MATL_D (7) END_E (9) (L=2) */
  cp[10] = isc * msc * prm->t0t2 * prm->t2t4 * prm->t4t6 * 1.0;  /* ROOT_S(0) ROOT_IR(2) MATL_D (4) MATL_ML(6) END_E (9) (L=2) */
  cp[11] = msc * isc * prm->t0t3 * prm->t3t5 * prm->t5t7 * 1.0;  /* ROOT_S(0) MATL_ML(3) MATL_IL(5) MATL_D (7) END_E (9) (L=2) */
  /* skipped because MATL_IL(8) is detached:                        ROOT_S(0) MATL_ML(3) MATL_D (7) MATL_IL(8) END_E (9) (L=2) */
  /* skipped because MATL_IL(8) is detached:                        ROOT_S(0) MATL_D (4) MATL_ML(6) MATL_IL(8) END_E (9) (L=2) */
  cp[12] = isc * msc * prm->t0t4 * prm->t4t5 * prm->t5t6 * 1.0;  /* ROOT_S(0) MATL_D (4) MATL_IL(5) MATL_ML(6) END_E (9) (L=2) */

  /* six possibilities with 2 inserts (and four others that are only impossible because MATL_IL(8) is detached */ 
  cp[13] = isc * isc * prm->t0t1 * prm->t1t1 * prm->t1t4 * prm->t4t7 * 1.0; /* ROOT_S(0) ROOT_IL(1) ROOT_IL(1) MATL_D (4) MATL_D (7) END_E  (9)           (L=2) */
  cp[14] = isc * isc * prm->t0t1 * prm->t1t2 * prm->t2t4 * prm->t4t7 * 1.0; /* ROOT_S(0) ROOT_IL(1) ROOT_IR(2) MATL_D (4) MATL_D (7) END_E  (9)           (L=2) */
  cp[15] = isc * isc * prm->t0t1 * prm->t1t4 * prm->t4t5 * prm->t5t7 * 1.0; /* ROOT_S(0) ROOT_IL(1) MATL_D (4) MATL_IL(5) MATL_D (7) END_E  (9)           (L=2) */
  /* skipped because MATL_IL(8) is detached:                                   ROOT_S(0) ROOT_IL(1) MATL_D (4) MATL_D (7) MATL_IL(8) END_E  (9)           (L=2) */

  cp[16] = isc * isc * prm->t0t2 * prm->t2t2 * prm->t2t4 * prm->t4t7 * 1.0; /* ROOT_S(0) ROOT_IR(2) ROOT_IR(2) MATL_D (4) MATL_D (7) END_E  (9)           (L=2) */
  cp[17] = isc * isc * prm->t0t2 * prm->t2t4 * prm->t4t5 * prm->t5t7 * 1.0; /* ROOT_S(0) ROOT_IR(2) MATL_D (4) MATL_IL(5) MATL_D (7) END_E  (9)           (L=2) */
  /* skipped because MATL_IL(8) is detached:                                   ROOT_S(0) ROOT_IR(2) MATL_D (4) MATL_D (7) MATL_IL(8) END_E  (9)           (L=2) */

  cp[18] = isc * isc * prm->t0t4 * prm->t4t5 * prm->t5t5 * prm->t5t7 * 1.0; /* ROOT_S(0) MATL_D (4) MATL_IL(5) MATL_IL(5) MATL_D (7) END_E  (9)           (L=2) */
  /* skipped because MATL_IL(8) is detached:                                   ROOT_S(0) MATL_D (4) MATL_IL(5) MATL_D (7) MATL_IL(8) END_E  (9)           (L=2) */

  /* skipped because MATL_IL(8) is detached:                                   ROOT_S(0) MATL_D (4) MATL_D (5) MATL_D (7) MATL_IL(8) MATL_IL(8) END_E (9) (L=2) */
  
  int i;
  for(i = 0; i <= 18; i++) { 
    printf("cp[%2d]: %.4f\n", i, cp[i]);
  }

  /* 2. Sum or max the total probability of L={0..2} aligned to one pass
        through the core model
   */
  if (do_cyk) 
    {
      pL[0] = cp[0];
      pL[1] = esl_vec_DMax(cp+1,  5);
      pL[2] = esl_vec_DMax(cp+6,  13);
      sc[0] = log(pL[0]) / log(2.);
      sc[1] = log(pL[1]) / log(2.);
      sc[2] = log(pL[2]) / log(2.);
      printf("CYK: L: %d prob: %10g sc: %10g\n", 0, pL[0], sc[0]);
      printf("CYK: L: %d prob: %10g sc: %10g\n", 1, pL[1], sc[1]);
      printf("CYK: L: %d prob: %10g sc: %10g\n", 2, pL[2], sc[2]);
    }
  else 
    {
      pL[0] = cp[0];
      pL[1] = esl_vec_DSum(cp+1,  5);
      pL[2] = esl_vec_DSum(cp+6,  13);
      sc[0] = log(pL[0]) / log(2.);
      sc[1] = log(pL[1]) / log(2.);
      sc[2] = log(pL[2]) / log(2.);
      printf("Ins: L: %d prob: %10g sc: %10g\n", 0, pL[0], sc[0]);
      printf("Ins: L: %d prob: %10g sc: %10g\n", 1, pL[1], sc[1]);
      printf("Ins: L: %d prob: %10g sc: %10g\n", 2, pL[2], sc[2]);
    }
  
  return eslOK;
}


  




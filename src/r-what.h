/* r-what.h
 *
 * Header support for r-what project
 *
 */

#ifndef R_WHAT_H_INCLUDED
#define R_WHAT_H_INCLUDED

/* Structure: PA_t (Partial Alignment)
 *
 * A minimal structure to track a partial SCFG alignment
 *
 * Does NOT maintain a full parsetree
 */
typedef struct partialalignment_s {
  int init_v;		/* Initial model state */
  int init_j;		/* Initial right sequence position */
  int init_d;		/* Initial sequence length */
  int init_i;		/* Initial left sequence position */
			/* While i is typically calculated from j and d,
                         * a separate variable is needed for marginal
                         * extension, where j and d are unknown         */

  int cur_v;
  int cur_j;
  int cur_d;
  int cur_i;

  int temp_v;		/* holder vars for provisional extension          */
  int temp_j;		/* Use with extreme caution - which ones have     */
  int temp_d;		/* defined values and whether they should replace */
  int temp_i;		/* init or cur equivalents depends on the         */
                        /* direction of extension.                        */
  int need_commit;	/* 0/1 - are there values in temp which need to   */
			/* be permanently recorded in init/cur            */

  float current_sc;
  float upper_bound_sc;

  int terminated;

} PA_t;

typedef struct branchedpartialalignment_s {
  PA_t  *chunk;
  struct branchedpartialalignment_s *left_child;
  struct branchedpartialalignment_s *right_child;

  float current_sc;
  float upper_bound_sc;
  int terminated;
} BPA_t;

PA_t*  PA_Copy(PA_t *orig);
BPA_t* BPA_Copy(BPA_t *orig);
float  BPA_Current_Score(BPA_t *root);
float  BPA_Upper_Bound(BPA_t *root);

void ConsensusD(CM_t *cm, int v, int *consensus_d);

void MaxSubsequenceScore(CM_t *cm, int W, float ***ret_max_sc);

PA_t* AstarExtension(CM_t *cm, char *dsq, int init_v, int init_j, int lower_d, int upper_d,
    	float init_sc, float **max_sc, float cutoff);

float LeftMarginalScore(float *esc, char syml);
float RightMarginalScore(float *esc, char symr);
void MarginalLeftInsideExtend(CM_t *cm, char *dsq, BPA_t *root, int rbound, float dropoff_sc, float *total_sc, float *delta_sc, int *commit, int *complete);

#endif /* R_WHAT_H_INCLUDED */

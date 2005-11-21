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

  int cur_v;
  int cur_j;
  int cur_d;

  float current_sc;
  float upper_bound_sc;

} PA_t;

PA_t* PA_Copy(PA_t *orig);

void MaxSubsequenceScore(CM_t *cm, int W, float ***ret_max_sc);

PA_t* AstarExtension(CM_t *cm, char *dsq, int init_v, int init_j, int lower_d, int upper_d,
    	float init_sc, float **max_sc);

#endif /* R_WHAT_H_INCLUDED */

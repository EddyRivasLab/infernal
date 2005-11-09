/* r-what.h
 *
 * Header support for r-what project
 *
 */

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

PA_t*
PA_Copy(PA_t *orig)
{
  PA_t *dup;
  dup = malloc(sizeof(PA_t));

  dup->init_v = orig->init_v;
  dup->init_j = orig->init_j;
  dup->init_d = orig->init_d;

  dup->cur_v = orig->cur_v;
  dup->cur_j = orig->cur_j;
  dup->cur_d = orig->cur_d;

  dup->current_sc = orig->current_sc;
  dup->upper_bound_sc = orig->upper_bound_sc;

  return dup;
}

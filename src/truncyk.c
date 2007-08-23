#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "easel.h"
#include "esl_stack.h"

#include "structs.h"
#include "funcs.h"
#include "cm_postprob.h"

#define BE_EFFICIENT 0
#define BE_PARANOID  1

#define USED_LOCAL_BEGIN 101
#define USED_EL          102

/* Structure: AlphaMats_t */
typedef struct alphamats_s {
   float ***J;
   float ***L;
   float ***R;
   float ***T;
} AlphaMats_t;

/* structure: BetaMats_t */
typedef struct betamats_s {
   float ***J;
   float  **L;
   float  **R;
   /* no T because T only applies at bifurcations, and beta/outside is only calculated on unbifurcated subgraphs */
} BetaMats_t;

/* Structure: ShadowMats_t */
typedef struct shadowmats_s {
   void ***J;
   void ***L;
   void ***Lmode;
   void ***R;
   void ***Rmode;
   void ***T;
} ShadowMats_t;

/* External API */
// CHANGEME - needed!

/* Divide and conquer */
float tr_generic_splitter(CM_t *cm, char *dsq, int L, Parsetree_t *tr,
                          int r, int vend, int i0, int j0, int allow_LM, int allow_RM);
float   tr_wedge_splitter(CM_t *cm, char *dsq, int L, Parsetree_t *tr,
                          int r, int z,    int i0, int j0, int allow_LM, int aloow_RM);
void        tr_v_splitter(CM_t *cm, char *dsq, int L, Parsetree_t *tr,
                          int r, int z,    int i0, int i1, int j1, int j0,
                          int useEL, int force_LM, int force_RM);

/* Alignment engine */
/* trinside is legacy, aviod use! */
float trinside (CM_t *cm, char *dsq, int L, int vroot, int vend, int i0, int j0, int do_full,
                void ****ret_shadow, void ****ret_L_shadow, void ****ret_R_shadow,
                void ****ret_T_shadow, void ****ret_Lmode_shadow, void ****ret_Rmode_shadow,
                int *ret_mode, int *ret_v, int *ret_i, int *ret_j);
float tr_inside(CM_t *cm, char *dsq, int L, int vroot, int vend, int i0, int j0, int do_full,
                AlphaMats_t *arg_alpha, AlphaMats_t *ret_alpha, 
                struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
                ShadowMats_t *ret_shadow, int *ret_mode, int *ret_v, int *ret_i, int *ret_j);
float tr_outside(CM_t *cm, char *dsq, int L, int vroot, int vend, int i0, int j0, int do_full,
                 BetaMats_t *arg_beta, BetaMats_t *ret_beta,
                 struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
                 int *ret_mode, int *ret_v, int *ret_j);

/* Traceback routine */
float trinsideT(CM_t *cm, char *dsq, int L, Parsetree_t *tr, int r, int z,
                int i0, int j0, int allow_begin, int *dmin, int *dmax);

/* Function: TrCYKInside()
 * Author:   DLK
 *
 * Purpose:  Full CYK alignment for truncated sequences
 *           with traceback
 * 
 *           Based on CYKInside()
 *
 * Args:     cm      - the covariance model
 *           dsq     - the sequence, 1..L
 *           L       - length of the sequence
 *           r       - root of subgraph to align to target subseq (usually 0, the model's root)
 *           i0      - start of target subsequence (usually 1, beginning of dsq)
 *           j0      - end of target subsequence (usually L, end of dsq)
 *           ret_tr  - RETURN: traceback (pass NULL if trace isn't wanted)
 *           dmin    - minimum d bound for each state v (NULL if non-banded)
 *           dmax    - maximum d bound for each state v (NULL if non-banded)
 *
 * Returns;  score of the alignment in bits
 */
float
TrCYKInside(CM_t *cm, char *dsq, int L, int r, int i0, int j0, Parsetree_t **ret_tr,
            int *dmin, int *dmax)
{
   Parsetree_t *tr;
   int          z;
   float        sc;

   /* Check input parameters */
   if ( cm->stid[r] != ROOT_S )
   {
      if (! (cm->flags & CM_LOCAL_BEGIN)) Die("internal error: we're not in local mode, but r is not root");
      if ( (cm->stid[r] != MATP_MP) &&
           (cm->stid[r] != MATL_ML) &&
           (cm->stid[r] != MATR_MR) &&
           (cm->stid[r] != BIF_B  )    )  Die("internal error: trying to do a local begin at a non-mainline start");
   }

   /* Create parse tree and initialize */
   tr = CreateParsetree();
   InsertTraceNode(tr, -1, TRACE_LEFT_CHILD, i0, j0, r);
   z = cm->M-1;
   sc = 0.0;

   /* If local begin is known */
   if ( r != 0 )
   {
      InsertTraceNode(tr, 0, TRACE_LEFT_CHILD, i0, j0, r);
      z = CMSubtreeFindEnd(cm, r);
      sc = cm->beginsc[r];
   }

   /* Solve by calling trinsideT() */
   sc += trinsideT(cm, dsq, L, tr, r, z, i0, j0, (r==0), dmin, dmax);

   if ( ret_tr != NULL ) *ret_tr = tr; else FreeParsetree(tr);

   return sc;
}

/* Function: tr_generic_splitter()
 * Author:   DLK
 *
 * Purpose:  Generic problem for divide-and-conquer
 *           Based closely on generic_splitter()
 *
 * Args:     
 *
 * Returns:
 */
float
tr_generic_splitter(CM_t *cm, char *dsq, int L, Parsetree_t *tr,
                    int r, int z, int i0, int j0, int allow_LM, int allow_RM)
{
   AlphaMats_t alpha;
   BetaMats_t  beta;
   struct deckpool_s *pool;
   int        v,w,y;
   int        wend, yend;
   int        jp;
   int        W;
   float      sc;
   int        j,d,k;
   float      best_sc;
   int        best_j, best_d, best_k;
   int        v_mode, w_mode, y_mode;
   int        b1_mode, b2mode;
   int        b1_v, b1_i, b1_j;
   int        b2_v, b2_i, b2_j;
   float      b1_sc, b2_sc;

   /* Case 1: problem size is small; solve with trinsideT()
    * size calculation is heuristic based on size of insideT() */
   if (5*insideT_size(cm, L, r, z, i0, j0) < RAMLIMIT)
   {
      sc = trinsideT(cm, dsq, L, tr, r, z, i0, j0, (r==0), NULL, NULL);
      return sc;
   }

   /* Case 2: find a bifurcation */
   for (v = r; v <= z-5; v++)
   {  if (cm->sttype[v] == B_st) break; }

   /* Case 3: no bifurcations -> wedge problem */
   if (cm->sttype[v] != B_st)
   {
      if (cm->sttype[z] != E_st) Die("z in tr_generic_splitter not E_st - that ain't right");
      sc = tr_wedge_splitter(cm, dsq, L, tr, r, z, i0, j0, allow_LM, allow_RM);
      return sc;
   }

   /* Unusual cases dispatched, back to case 2 (bifurcation) */
   w = cm->cfirst[v];
   y = cm->cnum[v];
   if (w < y) { wend = y-1; yend = z; }
   else       { yend = w-1; wend = z; }

   /* Calculate alphas for w and y
    * also pick up best local begins in each subtree */
   b1_sc = tr_inside(cm, dsq, L, w, wend, i0, j0, BE_EFFICIENT, NULL, &alpha, NULL, &pool, NULL, &b1_mode, &b1_v, &b1_i, &b1_j);
   if (!allow_begin) b1_sc = IMPOSSIBLE;
   b2_sc = tr_inside(cm, dsq, L, y, yend, i0, j0, BE_EFFICIENT,alpha, &alpha, pool,  NULL, NULL, &b2_mode, &b2_v, &b2_i, &b2_j);
   if (!allow_begin) b2_sc = IMPOSSIBLE;

   /* Calculate beta; release pool */
   b3_sc = tr_outside(cm, dsq, L, r, v, i0, j0, BE_EFFICIENT, NULL, &beta, NULL, NULL, &b3_mode, &b3_v, &b3_j);

   /* OK, to the point of actually finding the best split
    * We have a lot more types of splits than the non-truncated
    * version, so we need a better way to keep track of them    */
   W = j0 - i0 + 1;
   best_sc = IMPOSSIBLE;
   for (jp = 0; jp <= W; jp ++)
   {
      j = i0 - 1 + jp;
      for (d = 0; d <= jp; d++)
      {
         for (k = 0; k <= d; k++)
         {
            if ( (sc = alpha->J[w][j-k][d-k] + alpha->J[y][j][k] + beta->J[v][j][d]) > best_sc )
            {
               best_sc = sc;
               best_k  = k;
               best_j  = j;
               best_d  = d;
               v_mode = 3; w_mode = 3; y_mode = 3;
            }
            if ( allow_LM )
            if ( (sc = alpha->J[w][j-k][d-k] + alpha->L[y][j][k] + beta->L[v][j-d+1]) > best_sc )
            {
               best_sc = sc;
               best_k  = k;
               best_j  = j;
               best_d  = d;
               v_mode = 2; w_mode = 3; y_mode = 2;
            }
            if ( allow_RM )
            if ( (sc = alpha->R[w][j-k][d-k] + alpha->J[y][j][k] + beta->R[v][j]) > best_sc )
            {
               best_sc = sc;
               best_k  = k;
               best_j  = j;
               best_d  = d;
               v_mode = 1; w_mode = 1; y_mode = 3;
            }
            if ( allow_LM && allow_RM )
            if ( (sc = alpha->R[w][j-k][d-k] + alpha->L[y][j][k]) > best_sc )
            {
               best_sc = sc;
               best_k  = k;
               best_j  = j;
               best_d  = d;
               v_mode = 0; w_mode = 1; y_mode = 2;
            }
         }

         if ( allow_LM )
         if ( (sc = alpha->L[w][j][d] + beta->L[v][j-d+1]) > best_sc )
         {
            best_sc = sc;
            best_k  = 0;
            best_j  = j;
            best_d  = d;
            v_mode = 2; w_mode = 2; y_mode = 0;
         }
         if ( allow_RM )
         if ( (sc = alpha->R[y][j][d] + beta->R[v][j]) > best_sc )
         {
            best_sc = sc;
            best_k  = 0;
            best_j  = j;
            best_d  = d;
            v_mode = 1; w_mode = 0; y_mode = 1;
         }
         if ( (sc = beta->J[cm->M][j][d]) > best_sc) /* Joint parent to EL */
         {
            best_sc = sc;
            best_k  = 0;
            best_j  = j;
            best_d  = d;
            v_mode = 3; w_mode = 0; y_mode = 0;
         }
      }
   }

   /* Check for local entry in one of the child sub-trees */
   if (r == 0)
   {
      if (b1_sc > best_sc)
      {
         best_sc = b1_sc;
         best_k  = b1_v;
         best_j  = b1_j;
         best_d  = b1_j - b1_i + 1;
         v_mode = 0; w_mode = b1_mode; y_mode = 0;
      }
      if (b2_sc > best_sc)
      {
         best_sc = b2_sc;
         best_k  = b2_v;
         best_j  = b2_j;
         best_d  = b2_j - b2_i + 1;
         v_mode = 0; w_mode = 0; y_mode = b2_mode;
      }
   }

   /* local hit in parent (must be marginal) */
   if (b3_sc > best_sc)
   {
      best_sc = b3_sc;
      best_k  = b3_v;
      best_j  = b3_j;
      v_mode = b3_mode; w_mode = 0; y_mode = 0;
   }

   /* Free alphas */
   free_vjd_matrix(alpha->J, cm->M, i0, j0);
   free_vjd_matrix(alpha->L, cm->M, i0, j0);
   free_vjd_matrix(alpha->R, cm->M, i0, j0);
   free_vjd_matrix(alpha->T, cm->M, i0, j0);
   free_vjd_matrix( beta->J, cm->M, i0, j0);
// need to free the other stuff in beta too... */

   /* Found the best path, now to interpret and sub-divide */
   if ( v_mode ) /* parent graph is non-empty */
   {
      if ( w_mode == 0 && y_mode == 0 ) /* local hit in parent (marginal) */
      {
         tr_v_splitter(cm, dsq, L, tr, r, b3_v, i0, best_j+1, best_j, j0, (best_k == -1), (v_mode == 2), (v_mode == 1));
         return best_sc;
      }
      else
      {
         tr_v_splitter(cm, dsq, L, tr, r, v, i0, best_j-best_d+1, best_j, j0, (best_k == -1), (v_mode == 2), (v_mode == 1));
      }
   }
   else if ( w_mode == 0 || y_mode == 0 ) /* local entry to one of the children */
   {
      if ( w_mode )
      {
         InsertTraceNodeWithMode(tr, tr->n-1, TRACE_LEFT_CHILD, b1_i, b1_j, b1_v, b1_mode);
         z = CMSubtreeFindEnd(cm, b1_v);
         tr_generic_splitter(cm, dsq, L, tr, b1_v, z, b1_i, b1_j, (b1_mode == 2), (b1_mode == 1));
         return best_sc;
      }
      else if ( y_mode )
      {
         InsertTraceNodeWithMode(tr, tr->n-1, TRACE_LEFT_CHILD, b2_i, b2_j, b2_v, b2_mode);
         z = CMSubtreeFindEnd(cm, b2_v);
         tr_generic_splitter(cm, dsq, L, tr, b2_v, z, b2_i, b2_j, (b2_mode == 2), (b2_mode == 1));
         return best_sc;
      }
      else Die("Danger, danger!\n");
   }
   else /* case T: parent is empty, but both children are non-empty */
   {
      InsertTraceNodeWithMode(tr, tr->n-1, TRACE_LEFT_CHILD, best_j - best_d + 1, best_j, v, 0);
   }

   tv = tr->n - 1;
   if ( w_mode )
   {
      InsertTraceNodeWithMode(tr, tv, TRACE_LEFT_CHILD, best_j - best_d + 1, best_j - best_k, w, w_mode);
      tr_generic_splitter(cm, dsq, L, tr, w, wend, best_j - best_d + 1, best_j - best_k, (w_mode == 2), (w_mode == 1));
   }
   else
   {
      InsertTraceNodeWithMode(tr, tr->n-1, TRACE_LEFT_CHILD, best_j - best_d + 1, best_j - best_d, w, w_mode);
      InsertTraceNodeWithMode(tr, tr->n-1, TRACE_LEFT_CHILD, best_j - best_d + 1, best_j - best_d, cm->M, 3);
   }

   if ( y_mode )
   {
      InsertTraceNodeWithMode(tr, tv, TRACE_RIGHT_CHILD, best_j - best_k + 1, best_j, y, y_mode);
      tr_generic_splitter(cm, dsq, L, tr, y, yend, best_j - best_k + 1, best_j, (y_mode == 2), (y_mode == 1));
   }
   else 
   {
      InsertTraceNodeWithMode(tr, tv, TRACE_RIGHT_CHILD, best_j + 1, best_j, y, y_mode);
      InsertTraceNodeWithMode(tr, tv, TRACE_RIGHT_CHILD, best_j + 1, best_j, cm->M, 3);
   }

   return best_sc;
}

/* Function: tr_wedge_splitter()
 * Author:   DLK
 *
 * Purpose:  Wedge problem for divide-and-conquer
 *           Based closely on wedge_splitter()
 *
 * Args:     
 *
 * Returns:
 */
float
tr_wedge_splitter(CM_t *cm, char *dsq, int L, Parsetree_t *tr,
                  int r, int z, int i0, int j0, int allow_LM, int allow_RM)
{
   Alphamats_t *alpha;
   Betamats_t  *beta;
   float        sc;
   float        best_sc;
   int          v,w,y;
   int          W;
   int          jp, j, d;
   int          best_v, best_j, best_d;
   int          p_mode, c_mode;
   int          midnode;
   float        b1_sc, b2_sc;
   int          b1_mode, b1_v, b1_i, b1_j;
   int          b2_mode, b2_v, b2_j;

   /* Special case: problem is small enough to be solved with traceback */
   if ( (cm->ndidx[z] == cm->ndidx[r] + 1) || 
        (5 * insideT_size(cm, L, r, z, i0, j0) < RAMLIMIT) )
   {
      sc = trinsideT(cm, dsq, L, tr, r, z, i0, j0, (r==0), NULL, NULL);
      return sc;
   }

   /* Calculate a midpoint to split at */
   midnode = cm->ndidx[r] + ((cm->ndidx[z] - cm->ndidx[r])/2);
   w = cm->nodemap[midnode];
   y = cm->cfirst[w] - 1;

   /* Get alphas and betas */
   b1_sc = tr_inside(cm, dsq, L, w, z, i0, j0, BE_EFFICIENT, NULL, &alpha, NULL, NULL, NULL,
             &b1_mode, &b1_v, &b1_i, &b1_j);
   if (r != 0) b1_sc = IMPOSSIBLE;
   b2_sc = tr_outside(cm, dsq, L, r, y, i0, j0, BE_EFFICIENT, NULL, &beta, NULL, NULL,
             &b2_mode, &b2_v, &b2_j);
   if ( b2_mode == 2 && !allow_LM ) b2_sc = IMPOSSIBLE;
   if ( b2_mode == 1 && !allow_RM ) b2_sc = IMPOSSIBLE;

   /* Find the split */
   W = j0 - i0 + 1;
   best_sc = IMPOSSIBLE;
   for (v = w; v <= y; v++)
   {
      for (jp = 0; jp <= W; jp++)
      {
         j = i0 - 1 + jp;
         for (d = 0; d <= jp; d++)
         {
            if ( (sc = alpha->J[v][j][d] + beta->J[v][j][d]) > best_sc )
            {
               best_sc = sc;
               best_v  = v;
               best_d  = d;
               best_j  = j;
               p_mode = 3; c_mode = 3;
            }
            if ( allow_LM )
            if ( (sc = alpha->J[v][j][d] + beta->L[v][j-d+1]) > best_sc )
            {
               best_sc = sc;
               best_v  = v;
               best_d  = d;
               best_j  = j;
               p_mode = 2; c_mode = 3;
            }
            if ( allow_LM )
            if ( (sc = alpha->L[v][j][d] + beta->L[v][j-d+1]) > best_sc )
            {
               best_sc = sc;
               best_v  = v;
               best_d  = d;
               best_j  = j;
               p_mode = 2; c_mode = 2;
            }
            if ( allow_RM )
            if ( (sc = alpha->J[v][j][d] + beta->R[v][j]) > best_sc )
            {
               best_sc = sc;
               best_v  = v;
               best_d  = d;
               best_j  = j;
               p_mode = 1; c_mode = 3;
            }
            if ( allow_RM )
            if ( (sc = alpha->R[v][j][d] + beta->R[v][j]) > best_sc )
            {
               best_sc = sc;
               best_v  = v;
               best_d  = d;
               best_j  = j;
               p_mode = 1; c_mode = 1;
            }
         }
      }
   }

   /* Special case: joint parent to EL */
   for (jp = 0; jp <= W; jp++)
   {
      j = i0 - 1 + jp;
      for (d = 0; d <= jp; d++)
      {
         if ( (sc = beta->J[cm->M][j][d]) > best_sc )
         {
            best_sc = sc;
            best_v  = -1;
            best_j  = j;
            best_d  = d;
            p_mode = 3; c_mode = 0;
         }
      }
   }

   /* Special case: parent empty, child has local hit */
   if (b1_sc > best_sc)
   {
      best_sc = b1_sc;
      best_v  = b1_v;
      best_j  = b1_j;
      best_d  = b1_d;
      p_mode = 0; c_mode = b1_mode;
   }

   /* Special case: child empty, parent has local hit */
   /* 1 and 2 are the only appropriate values for b2_mode */
   if (b2_sc > best_sc)
   {
      best_sc = b2_sc;
      best_v  = b2_v;
      best_j  = b2_j;
      best_d  = 0;
      p_mode = b2_mode; c_mode = 0;
   }

   /* Free alpha and beta */
   free_vjd_matrix(alpha->J, cm->M, i0, j0);
   free_vjd_matrix(alpha->L, cm->M, i0, j0);
   free_vjd_matrix(alpha->R, cm->M, i0, j0);
   free_vjd_matrix(alpha->T, cm->M, i0, j0);
   free_vjd_matrix( beta->J, cm->M, i0, j0);
/* again, the remainder of beta... */
   
   if ( p_mode )
   {
      if ( c_mode == 0 ) /* child empty */
      {
         tr_v_splitter(cm, dsq, L, tr, r, (p_mode == 3 ? w : b2_v), i0, best_j - best_d + 1, best_j, j0, (p_mode == 3), (p_mode == 2), (p_mode == 1));
         return best_sc;
      }
      else
      {
         tr_v_splitter(cm, dsq, L, tr, r, best_v, i0, best_j - best_d + 1, best_j, j0, FALSE, (p_mode == 2), (p_mode == 1));
      }
   }

   if ( c_mode )
   {
      tr_wedge_splitter(cm, dsq, L, tr, best_v, z, best_j - best_d + 1, best_j, (c_mode == 2), (c_mode == 1));
   }
   else /* parent and child both empty */
      Die("Danger, danger!\n");

   return best_sc;
}

/* Function: tr_v_splitter()
 * Author:   DLK
 *
 * Purpose:  'V problem' - closely based on v_splitter()
 *
 * Args:     
 *
 * Returns;  none
 */
void
tr_v_splitter(CM_t *cm, char *dsq, int L, Parsetree_t *tr, int r, int z, int i0, int i1,
              int j1, int j0, int useEL, int force_LM, int force_RM)
{
   Alphamats_t *alpha;
   Betamats_t  *beta;
   float        sc, best_sc;
   int          v, w, y;
   int          best_v, best_i, best_j;
   int          midnode;

// Recommend a special handler for the fully marginal cases (linear alg.)
   /*
   if ( force_LM)
   {

   }
   else if ( force_RM )
   {

   }
   */

   /* Special case: solve without splitting for small problems and boundary conditions */
   if (cm->ndidx[z] == cm->ndidx[r] + 1 || r == z ||
       5*vinsideT_size(cm, r, z, i0, i1, j1, j0) < RAMLIMIT)
   {
      tr_vinsideT(cm, dsq, L, tr, r, z, i0, i1, j1, j0, useEL, force_LM, force_RM);
      return;
   }

   /* Find split set */
   midnode = cm->ndidx[r] + ((cm->ndidx[z] - cm->ndidx[r])/2);
   w = cm->nodemap[midnode];
   y = cm->cfirst[w] - 1;

   /* Calculate alphas and betas */
   b_sc =  vinside(cm, dsq, L, w, z, i0, i1, j1, j0, useEL, force_LM, force_RM,
                   BE_EFFICIENT, NULL, &alpha, NULL, NULL, &b_mode, &b_v, &b_i, &b_j);
   if (r != 0) b1_sc = IMPOSSIBLE;
   voutside(cm, dsq, L, r, y, i0, i1, j1, j0, useEL, force_LM, force_RM, BE_EFFICIENT, NULL, &beta, NULL, NULL);

   /* Find our best split */
   best_sc = IMPOSSIBLE;
   for (v = w; v <= y; v++)
   {
      for (ip = 0; ip <= i1-i0; ip++)
      {
         for (jp = 0; jp <= j0-j1, jp++)
         {
            if ( ! force_LM && ! force_RM )
            if ( (sc = alpha->J[v][jp][ip] + beta->J[v][jp][ip]) > best_sc )
            {
               best_sc = sc;
               best_v  = v;
               best_i  = ip + i0;
               best_j  = jp + j1;
               p_mode = 3; c_mode = 3;
            }
            if ( ! force_LM && ! force_RM )
            if ( (sc = alpha->J[v][jp][ip] + beta->L[v][ip]) > best_sc )
            {
               best_sc = sc;
               best_v  = v;
               best_i  = ip + i0;
               best_j  = jp + j1;
               p_mode = 2; c_mode = 3;
            }
            if ( ! force_RM )
            if ( (sc = alpha->L[v][jp][ip] + beta->L[v][ip]) > best_sc )
            {
               best_sc = sc;
               best_v  = v;
               best_i  = ip + i0;
               best_j  = jp + j1;
               p_mode = 2; c_mode = 2;
            }
            if ( ! force_LM && ! force_RM )
            if ( (sc = alpha->J[v][jp][ip] + beta->R[v][jp]) > best_sc )
            {
               best_sc = sc;
               best_v  = v;
               best_i  = ip + i0;
               best_j  = jp + j1;
               p_mode = 1; c_mode = 3;
            }
            if ( ! force_LM )
            if ( (sc = alpha->R[v][jp][ip] + beta->R[v][jp]) > best_sc )
            {
               best_sc = sc;
               best_v  = v;
               best_i  = ip + i0;
               best_j  = jp + j1;
               p_mode = 1; c_mode = 1;
            }
         }
      }
   }

   /* check EL */
   if ( useEL )
   {
      for (ip = 0; ip <= i1-i0; ip++)
      {
         for (jp = 0; jp <= j0-j1; jp++)
         {
            if ( (sc = beta->J[cm->M][jp][ip]) > best_sc )
            {
               best_sc = sc;
               best_v  = cm->M;
               best_i  = ip + i0;
               best_j  = jp + j1;
               p_mode = 3; c_mode = 0;
            }
         }
      }
   }

   /* check local begin */
   if (b_sc > best_sc)
   {
      best_sc = b_sc;
      best_v  = b_v;
      best_i  = b_i;
      best_j  = b_j;
      p_mode = 0; c_mode = b_mode;
   }

   /* Free memory */
   free_vji_matrix(alpha->J, cm->M, j1, j0);
   free_vji_matrix(alpha->L, cm->M, j1, j0);
   free_vji_matrix(alpha->R, cm->M, j1, j0);
   free_vji_matrix(alpha->T, cm->M, j1, j0);
   free_vji_matrix( beat->J, cm->M, j1, j0);
// Free the other betas....

   /* Interpret and subdivide */
   if ( p_mode )
   {
      if ( c_mode )
      {
         tr_v_splitter(cm, dsq, L, tr, r, best_v, i0, best_i, best_j, j0, FALSE, force_LM, force_RM);
         tr_v_splitter(cm, dsq, L, tr, best_v, z, best_i, i1, j1, best_j, useEL, force_LM, force_RM);  
      }
      else
      {
         tr_v_splitter(cm, dsq, L, tr, r, w, i0, best_i, best_j, j0, TRUE, FALSE, FALSE);
      }
   }
   else
   {
      if (best_v != z)
      {
         InsertTraceNodeWithMode(tr, tr->n-1, TRACE_LEFT_CHILD, best_i, best_j, best_v, best_mode);
      }
      tr_v_splitter(cm, dsq, L, tr, best_v, z, best_i, i1, j1, best_j, useEL, force_LM, force_RM);
   }

   return;
}

/* Legacy wrapper */
float
trinside (CM_t *cm, char *dsq, int L, int vroot, int vend, int i0, int j0, int do_full,
          void ****ret_shadow, void ****ret_L_shadow, void ****ret_R_shadow,
          void ****ret_T_shadow, void ****ret_Lmode_shadow, void ****ret_Rmode_shadow,
          int *ret_mode, int *ret_v, int *ret_i, int *ret_j)
{
   float sc;
   ShadowMats_t shadow;

   shadow.J = *ret_shadow;
   shadow.L = *ret_L_shadow;
   shadow.R = *ret_R_shadow;
   shadow.T = *ret_T_shadow;
   shadow.Lmode = *ret_Lmode_shadow;
   shadow.Rmode = *ret_Rmode_shadow;

   sc = tr_inside(cm, dsq, L, vroot, vend, i0, j0, do_full,
                  NULL, NULL, NULL, NULL, &shadow,
                  ret_mode, ret_v, ret_i, ret_j);
   return sc;
}

/* Function: tr_inside()
 * Author:   DLK
 *
 * Purpose:  inside phase of CYK on truncated sequence
 *           based on inside()
 *
 *           Score matrices and deckpool are only managed
 *           within this func, not available to caller
 *
 * Args:     cm      - the covariance model
 *           dsq     - the sequence, 1..L
 *           L       - length of the sequence
 *           r       - root of subgraph to align to target subseq (usually 0, the model's root)
 *           z       - last state of the subgraph
 *           i0      - start of target subsequence (usually 1, beginning of dsq)
 *           j0      - end of target subsequence (usually L, end of dsq)
 *           do_full - if TRUE, save all decks rather than re-using
 * 
 * Returns:  Score of the optimal alignment
 */
float
tr_inside(CM_t *cm, char *dsq, int L, int vroot, int vend, int i0, int j0, int do_full,
                AlphaMats_t *arg_alpha, AlphaMats_t *ret_alpha, 
                struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
                ShadowMats_t *ret_shadow, int *ret_mode, int *ret_v, int *ret_i, int *ret_j)
{
   float  **end;
   int      nends;
   int     *touch;
   int      v,y,z;
   int      j,d,i,k;
   float    sc,tsc;
   int      yoffset;
   int      W;
   int      jp;
   void  ***shadow;
   void  ***L_shadow;
   void  ***R_shadow;
   void  ***T_shadow;
   int   ***Lmode_shadow;
   int   ***Rmode_shadow;
   int    **kshad;
   char   **yshad;
   int      r_v, r_i, r_j, r_mode;
   float    r_sc;
   int      model_len;
   float    bsc;

   float ***alpha;
   float ***L_alpha;
   float ***R_alpha;
   float ***T_alpha;

   if ( arg_alpha == NULL )
   {
      alpha = NULL;
      L_alpha = NULL;
      R_alpha = NULL;
      T_alpha = NULL;
   }
   else
   {
        alpha = arg_alpha->J;
      L_alpha = arg_alpha->L;
      R_alpha = arg_alpha->R;
      T_alpha = arg_alpha->T;
   }

   /*Initialization */
   r_v = -1;
   r_i = i0;
   r_j = j0;
   r_mode = 3;
   r_sc = IMPOSSIBLE;
   W = j0-i0+1;
   model_len = 0;
   for ( v = vend; v >= vroot; v-- )
   {
      if      ( cm->stid[v] == MATP_MP ) model_len += 2;
      else if ( cm->stid[v] == MATL_ML ) model_len += 1;
      else if ( cm->stid[v] == MATR_MR ) model_len += 1;
   }
// CHANGEME!
// This correction doesn't belong here - it should be out at the top level API
// because in D and C, inside may be called more that once
// Existing drivers that call directly into the alignment engine (like the
// one for parameterizing mu and lambda) will need to be adjusted so that they
// get the correction as well
   /* 2.0 instead of 2 to force floating point division, not integer division */
   bsc = sreLOG2(2.0/(model_len*(model_len+1)));

   /* Make a deckpool */
   if ( dpool == NULL ) dpool = deckpool_create();
   if (! deckpool_pop(dpool, &end) )
   {  end = alloc_vjd_deck(L, i0, j0); }
   nends = CMSubtreeCountStatetype(cm, vroot, E_st);
   for ( jp=0; jp<=W; jp++ )
   {
      j = i0+jp-1;
      end[j][0] = 0.0;
      for ( d=1; d<=jp; d++ ) { end[j][d] = IMPOSSIBLE; }
   }

   /* Create score matrices */
   if ( alpha == NULL )
   {
      alpha = MallocOrDie(sizeof(float **) * (cm->M+1));
      for ( v=0; v<=cm->M; v++ ) { alpha[v] = NULL; }
   }
   if ( L_alpha == NULL )
   {
      L_alpha = MallocOrDie(sizeof(float **) * (cm->M+1));
      for ( v=0; v<=cm->M; v++ ) { L_alpha[v] = NULL; }
   }
   if ( R_alpha == NULL )
   {
      R_alpha = MallocOrDie(sizeof(float **) * (cm->M+1));
      for ( v=0; v<=cm->M; v++ ) { R_alpha[v] = NULL; }
   }
   if ( T_alpha == NULL )
   {
      T_alpha = MallocOrDie(sizeof(float **) * (cm->M+1));
      for ( v=0; v<=cm->M; v++ ) { T_alpha[v] = NULL; }
   }

   touch = MallocOrDie(sizeof(int) *cm->M);
   for ( v=0;      v<vroot; v++ ) { touch[v] = 0; }
   for ( v=vroot;  v<=vend; v++ ) { touch[v] = cm->pnum[v]; }
   for ( v=vend+1; v<cm->M; v++ ) {touch[v] = 0; }

   /* Create shadow matrices */
   if ( ret_shadow != NULL ) 
   {
      shadow = (void ***) MallocOrDie(sizeof(void **) * cm->M);
      for ( v=0; v<cm->M; v++ ) { shadow[v] = NULL; }

      L_shadow = (void ***) MallocOrDie(sizeof(void **) * cm->M);
      for ( v=0; v<cm->M; v++ ) { L_shadow[v] = NULL; }

      R_shadow = (void ***) MallocOrDie(sizeof(void **) * cm->M);
      for ( v=0; v<cm->M; v++ ) { R_shadow[v] = NULL; }

      T_shadow = (void ***) MallocOrDie(sizeof(void **) * cm->M);
      for ( v=0; v<cm->M; v++ ) { T_shadow[v] = NULL; }

      Lmode_shadow = (int ***) MallocOrDie(sizeof(int **) * cm->M);
      for ( v=0; v<cm->M; v++ ) { Lmode_shadow[v] = NULL; }

      Rmode_shadow = (int ***) MallocOrDie(sizeof(int **) * cm->M);
      for ( v=0; v<cm->M; v++ ) { Rmode_shadow[v] = NULL; }
   }

   /* Main recursion */
   for ( v = vend; v >= vroot; v-- )
   {
      if ( cm->sttype[v] == E_st )
      {
         alpha[v] = end;
         L_alpha[v] = end;
         R_alpha[v] = end;
         continue;
      }
      /* Assign alpha decks */
      if (! deckpool_pop(dpool, &(alpha[v])) )
      {  alpha[v] = alloc_vjd_deck(L, i0, j0); }
      if (! deckpool_pop(dpool, &(L_alpha[v])) )
      {  L_alpha[v] = alloc_vjd_deck(L, i0, j0); }
      if (! deckpool_pop(dpool, &(R_alpha[v])) )
      {  R_alpha[v] = alloc_vjd_deck(L, i0, j0); }
      if ( (cm->sttype[v] == B_st) && (! deckpool_pop(dpool, &(T_alpha[v])) ) )
      {  T_alpha[v] = alloc_vjd_deck(L, i0, j0); }

      /* Assign shadow decks */
      if ( ret_shadow != NULL )
      {
         if ( cm->sttype[v] == B_st )
         {
            kshad = alloc_vjd_kshadow_deck(L, i0, j0);
            shadow[v] = (void **) kshad;
         }
         else
         {
            yshad = alloc_vjd_yshadow_deck(L, i0, j0);
            shadow[v] = (void **) yshad;
         }

         if ( cm->sttype[v] == B_st )
         {
            kshad = alloc_vjd_kshadow_deck(L, i0, j0);
            L_shadow[v] = (void **) kshad;
         }
         else
         {
            yshad = alloc_vjd_yshadow_deck(L, i0, j0);
            L_shadow[v] = (void **) yshad;
         }

         if ( cm->sttype[v] == B_st )
         {
            kshad = alloc_vjd_kshadow_deck(L, i0, j0);
            R_shadow[v] = (void **) kshad;
         }
         else
         {
            yshad = alloc_vjd_yshadow_deck(L, i0, j0);
            R_shadow[v] = (void **) yshad;
         }

         if ( cm->sttype[v] == B_st )
         {
            kshad = alloc_vjd_kshadow_deck(L, i0, j0);
            T_shadow[v] = (void **) kshad;
         }

         kshad = alloc_vjd_kshadow_deck(L, i0, j0);
         Lmode_shadow[v] = (int **) kshad;

         kshad = alloc_vjd_kshadow_deck(L, i0, j0);
         Rmode_shadow[v] = (int **) kshad;
      }

      if ( cm->sttype[v] == D_st || cm->sttype[v] == S_st )
      {
         for ( jp=0; jp<=W; jp++ )
         {
            j = i0-1+jp;
            for ( d=0; d<=jp; d++ )
            {
               y = cm->cfirst[v];
               alpha[v][j][d]   = cm->endsc[v] + (cm->el_selfsc * (d));
/*
               if (cm->sttype[v] == S_st) alpha[v][j][d]   = cm->endsc[v] + (cm->el_selfsc * (d));
               else                       alpha[v][j][d]   = IMPOSSIBLE;
*/
               L_alpha[v][j][d] = IMPOSSIBLE;
               R_alpha[v][j][d] = IMPOSSIBLE;
               if ( ret_shadow   != NULL )
               {
                  ((char **)  shadow[v])[j][d] = USED_EL;
                  /* Set USED_EL to prevent a traceback bug */
                  ((char **)L_shadow[v])[j][d] = USED_EL;
                  ((char **)R_shadow[v])[j][d] = USED_EL;
               }

               for ( yoffset=0; yoffset<cm->cnum[v]; yoffset++ )
               {
                  if ( (sc = alpha[y+yoffset][j][d] + cm->tsc[v][yoffset]) > alpha[v][j][d] )
                  {
                     alpha[v][j][d] = sc;
                     if ( ret_shadow != NULL ) ((char **)shadow[v])[j][d] = yoffset;
                  }
                  if ( (sc = L_alpha[y+yoffset][j][d] + cm->tsc[v][yoffset]) > L_alpha[v][j][d] )
                  {
                     L_alpha[v][j][d] = sc;
                     if ( ret_shadow != NULL ) ((char **)L_shadow[v])[j][d] = yoffset;
                     if ( ret_shadow != NULL ) Lmode_shadow[v][j][d] = 2;
                  }
                  if ( (sc = R_alpha[y+yoffset][j][d] + cm->tsc[v][yoffset]) > R_alpha[v][j][d] )
                  {
                     R_alpha[v][j][d] = sc;
                     if ( ret_shadow != NULL ) ((char **)R_shadow[v])[j][d] = yoffset;
                     if ( ret_shadow != NULL ) Rmode_shadow[v][j][d] = 1;
                  }
               }

               if ( d == 0 )
               {
                  L_alpha[v][j][d] = IMPOSSIBLE;
                  R_alpha[v][j][d] = IMPOSSIBLE;
               }

               if (   alpha[v][j][d] < IMPOSSIBLE ) {   alpha[v][j][d] = IMPOSSIBLE; }
               if ( L_alpha[v][j][d] < IMPOSSIBLE ) { L_alpha[v][j][d] = IMPOSSIBLE; }
               if ( R_alpha[v][j][d] < IMPOSSIBLE ) { R_alpha[v][j][d] = IMPOSSIBLE; }
            }
         }
      }
      else if ( cm->sttype[v] == B_st )
      {
         for ( jp=0; jp<=W; jp++ )
         {
            j = i0-1+jp;
            for ( d=0; d<=jp; d++ )
            {
               int allow_L_exit = 0;
               int allow_R_exit = 0;
               int allow_J_exit = 0;
               y = cm->cfirst[v];
               z = cm->cnum[v];

               alpha[v][j][d]   =   alpha[y][j][d] +   alpha[z][j][0];
               L_alpha[v][j][d] = L_alpha[y][j][d]                   ;
               R_alpha[v][j][d] =                  + R_alpha[z][j][d];
               if ( ret_shadow != NULL )
               {
                  ((int **)  shadow[v])[j][d] = 0;
                  ((int **)L_shadow[v])[j][d] = 0;
                  ((int **)R_shadow[v])[j][d] = d;
                  Lmode_shadow[v][j][d] = 2;
                  Rmode_shadow[v][j][d] = 1;
               }

               if  ( (sc = alpha[y][j][d]                   ) > L_alpha[v][j][d] )
               {
                  L_alpha[v][j][d] = sc;
                  if ( ret_shadow != NULL ) { ((int **)L_shadow[v])[j][d] = 0; }
                  if ( ret_shadow != NULL ) { Lmode_shadow[v][j][d] = 3; }
               }
 
               for ( k=1; k<=d; k++ )
               {
                  if ( (sc = alpha[y][j-k][d-k] + alpha[z][j][k]) > alpha[v][j][d] )
                  {
                     alpha[v][j][d] = sc;
                     if (ret_shadow != NULL ) { ((int **)shadow[v])[j][d] = k; }
                     if (k == d) allow_J_exit = 0;
                     else        allow_J_exit = 1;
                  }
                  if ( (sc = alpha[y][j-k][d-k] + L_alpha[z][j][k]) > L_alpha[v][j][d] )
                  {
                     L_alpha[v][j][d] = sc;
                     if ( ret_shadow != NULL ) { ((int **)L_shadow[v])[j][d] = k; }
                     if ( ret_shadow != NULL ) { Lmode_shadow[v][j][d] = 3; }
                     allow_L_exit = 1;
                  }
               }
               if ( (sc =                        alpha[z][j][d]) > R_alpha[v][j][d] )
               {
                  R_alpha[v][j][d] = sc;
                  if ( ret_shadow != NULL ) { ((int **)R_shadow[v])[j][d] = d; }
                  if ( ret_shadow != NULL ) { Rmode_shadow[v][j][d] = 3; }
               }
               for ( k=0; k< d; k++ )
               {
                  if ( (sc = R_alpha[y][j-k][d-k] + alpha[z][j][k]) > R_alpha[v][j][d] )
                  {
                     R_alpha[v][j][d] = sc;
                     if ( ret_shadow != NULL ) { ((int **)R_shadow[v])[j][d] = k; }
                     if ( ret_shadow != NULL ) { Rmode_shadow[v][j][d] = 3; }
                     allow_R_exit = 1;
                  }
               }

               if (d == 0) {
                  L_alpha[v][j][d] = IMPOSSIBLE;
                  R_alpha[v][j][d] = IMPOSSIBLE;
               }

               if (d >= 2) {
                 T_alpha[v][j][d] = R_alpha[y][j-1][d-1] + L_alpha[z][j][1];
                 if ( ret_shadow != NULL ) { ((int **)T_shadow[v])[j][d] = 1; }
                 for ( k=2; k<d; k++ )
                 {
                    if ( (sc = R_alpha[y][j-k][d-k] + L_alpha[z][j][k]) > T_alpha[v][j][d] )
                    {
                       T_alpha[v][j][d] = sc;
                       if ( ret_shadow != NULL) { ((int **)T_shadow[v])[j][d] = k; }
                    }
                 }
               }
               else {
                 T_alpha[v][j][d] = IMPOSSIBLE;
               }

               if (   alpha[v][j][d] < IMPOSSIBLE ) {   alpha[v][j][d] = IMPOSSIBLE; }
               if ( L_alpha[v][j][d] < IMPOSSIBLE ) { L_alpha[v][j][d] = IMPOSSIBLE; }
               if ( R_alpha[v][j][d] < IMPOSSIBLE ) { R_alpha[v][j][d] = IMPOSSIBLE; }

               /* Shouldn't allow exit from marginal B if one of the children is NULL, sinee that is covered by the */
               /* root of the other child, and we haven't added anything above the bifurcation */
               if ((  alpha[v][j][d] + bsc > r_sc) && (allow_J_exit) )
               { r_mode = 3; r_v = v; r_j = j; r_i = j-d+1; r_sc = bsc +   alpha[v][j][d]; }
               if ((L_alpha[v][j][d] + bsc > r_sc) && (allow_L_exit) )
               { r_mode = 2; r_v = v; r_j = j; r_i = j-d+1; r_sc = bsc + L_alpha[v][j][d]; }
               if ((R_alpha[v][j][d] + bsc > r_sc) && (allow_R_exit) )
               { r_mode = 1; r_v = v; r_j = j; r_i = j-d+1; r_sc = bsc + R_alpha[v][j][d]; }
               if ( T_alpha[v][j][d] + bsc > r_sc )
               { r_mode = 0; r_v = v; r_j = j; r_i = j-d+1; r_sc = bsc + T_alpha[v][j][d]; }
            }
         }
      }
      else if ( cm->sttype[v] == MP_st )
      {
         for ( jp=0; jp<=W; jp++ )
         {
            j = i0-1+jp;
            y = cm->cfirst[v];
            alpha[v][j][0] = IMPOSSIBLE;
            L_alpha[v][j][0] = IMPOSSIBLE;
            R_alpha[v][j][0] = IMPOSSIBLE;;
            if ( jp > 0 ) {
               alpha[v][j][1] = IMPOSSIBLE;
               if ( dsq[j] < Alphabet_size ) { L_alpha[v][j][1] =  LeftMarginalScore(cm->esc[v],dsq[j]); }
               else { Die("Still can't deal with marginalizing degenerate residues! x %d dsq[x] %d",j,dsq[j]); }
               if ( dsq[j] < Alphabet_size ) { R_alpha[v][j][1] = RightMarginalScore(cm->esc[v],dsq[j]); }
               else { Die("Still can't deal with marginalizing degenerate residues! x %d dsq[x] %d",j,dsq[j]); }
            if ( ret_shadow != NULL ) { ((char **)L_shadow[v])[j][1] = USED_EL; }
            if ( ret_shadow != NULL ) { ((char **)R_shadow[v])[j][1] = USED_EL; }
            }
            for ( d=2; d<=jp; d++)
            {
               alpha[v][j][d] = cm->endsc[v] + (cm->el_selfsc * (d-2));
               L_alpha[v][j][d] = IMPOSSIBLE;
               R_alpha[v][j][d] = IMPOSSIBLE;
               if ( ret_shadow != NULL ) { ((char **)shadow[v])[j][d] = USED_EL; }

               for ( yoffset=0; yoffset<cm->cnum[v]; yoffset++ )
               {
                  if ( (sc = alpha[y+yoffset][j-1][d-2] + cm->tsc[v][yoffset]) > alpha[v][j][d] )
                  {
                     alpha[v][j][d] = sc;
                     if ( ret_shadow != NULL ) { ((char **)shadow[v])[j][d] = yoffset; }
                  }

                  if ( (sc = L_alpha[y+yoffset][j][d-1] + cm->tsc[v][yoffset]) > L_alpha[v][j][d])
                  {
                     L_alpha[v][j][d] = sc;
                     if ( ret_shadow != NULL ) { ((char **)L_shadow[v])[j][d] = yoffset; }
                     if ( ret_shadow != NULL ) { Lmode_shadow[v][j][d] = 2; }
                  }
                  if ( (sc = alpha[y+yoffset][j][d-1] + cm->tsc[v][yoffset]) > L_alpha[v][j][d])
                  {
                     L_alpha[v][j][d] = sc;
                     if ( ret_shadow != NULL ) { ((char **)L_shadow[v])[j][d] = yoffset; }
                     if ( ret_shadow != NULL ) { Lmode_shadow[v][j][d] = 3; }
                  }

                  if ( (sc = R_alpha[y+yoffset][j-1][d-1] + cm->tsc[v][yoffset]) > R_alpha[v][j][d])
                  {
                     R_alpha[v][j][d] = sc;
                     if ( ret_shadow != NULL ) { ((char **)R_shadow[v])[j][d] = yoffset; }
                     if ( ret_shadow != NULL ) { Rmode_shadow[v][j][d] = 1; }
                  }
                  if ( (sc = alpha[y+yoffset][j-1][d-1] + cm->tsc[v][yoffset]) > R_alpha[v][j][d])
                  {
                     R_alpha[v][j][d] = sc;
                     if ( ret_shadow != NULL ) { ((char **)R_shadow[v])[j][d] = yoffset; }
                     if ( ret_shadow != NULL ) { Rmode_shadow[v][j][d] = 3; }
                  }
               }

               i = j-d+1;
               if ( dsq[i] < Alphabet_size && dsq[j] < Alphabet_size )
               {  alpha[v][j][d] += cm->esc[v][(int) (dsq[i]*Alphabet_size+dsq[j])]; }
               else
               {  alpha[v][j][d] += DegeneratePairScore(cm->esc[v], dsq[i], dsq[j]); }
               if ( dsq[i] < Alphabet_size )
               { L_alpha[v][j][d] +=  LeftMarginalScore(cm->esc[v],dsq[i]); }
               else
               { Die("Still can't deal with marginalizing degenerate residues! x %d dsq[x] %d",i,dsq[i]); }
               if ( dsq[j] < Alphabet_size )
               { R_alpha[v][j][d] += RightMarginalScore(cm->esc[v],dsq[j]); }
               else
               { Die("Still can't deal with marginalizing degenerate residues! x %d dsq[x] %d",j,dsq[j]); }

               if (   alpha[v][j][d] < IMPOSSIBLE ) {   alpha[v][j][d] = IMPOSSIBLE; }
               if ( L_alpha[v][j][d] < IMPOSSIBLE ) { L_alpha[v][j][d] = IMPOSSIBLE; }
               if ( R_alpha[v][j][d] < IMPOSSIBLE ) { R_alpha[v][j][d] = IMPOSSIBLE; }
            }

            for ( d = 1; d <= jp; d++ )
            {
               if (   alpha[v][j][d] + bsc > r_sc ) { r_mode = 3; r_v = v; r_j = j; r_i = j-d+1; r_sc = bsc +   alpha[v][j][d]; }
               if ( L_alpha[v][j][d] + bsc > r_sc ) { r_mode = 2; r_v = v; r_j = j; r_i = j-d+1; r_sc = bsc + L_alpha[v][j][d]; }
               if ( R_alpha[v][j][d] + bsc > r_sc ) { r_mode = 1; r_v = v; r_j = j; r_i = j-d+1; r_sc = bsc + R_alpha[v][j][d]; }
            }
         }
      }
      else if ( cm->sttype[v] == IL_st || cm->sttype[v] == ML_st )
      {
         for ( jp = 0; jp <= W; jp++ )
         {
            j = i0-1+jp;
            y = cm->cfirst[v];

            alpha[v][j][0] = IMPOSSIBLE;
            L_alpha[v][j][0] = IMPOSSIBLE;
            R_alpha[v][j][0] = IMPOSSIBLE;

            for ( d = 1; d <= jp; d++ )
            {
               alpha[v][j][d]   = cm->endsc[v] + (cm->el_selfsc * (d-1));
               if (d == 1) L_alpha[v][j][d] = 0.0;
               else        L_alpha[v][j][d] = IMPOSSIBLE;
               R_alpha[v][j][d] = IMPOSSIBLE;
               if ( ret_shadow   != NULL ) { ((char **)  shadow[v])[j][d] = USED_EL; }
               if ( ret_shadow != NULL ) { ((char **)L_shadow[v])[j][d] = USED_EL; }
               if ( ret_shadow != NULL ) { Lmode_shadow[v][j][d] = 3; }

               for ( yoffset = 0; yoffset < cm->cnum[v]; yoffset++ )
               {
                  if  ( (sc = alpha[y+yoffset][j][d-1] + cm->tsc[v][yoffset]) > alpha[v][j][d] )
                  {
                     alpha[v][j][d] = sc;
                     if ( ret_shadow != NULL ) { ((char **)shadow[v])[j][d] = yoffset; }
                  }

                  if  ( d > 1 )
                  if  ( (sc = L_alpha[y+yoffset][j][d-1] + cm->tsc[v][yoffset]) > L_alpha[v][j][d] )
                  {
                     L_alpha[v][j][d] = sc;
                     if ( ret_shadow != NULL ) { ((char **)L_shadow[v])[j][d] = yoffset; }
                     if ( ret_shadow != NULL ) { Lmode_shadow[v][j][d] = 2; }
                  }

                  if  ( (sc = R_alpha[y+yoffset][j][d] + cm->tsc[v][yoffset]) > R_alpha[v][j][d] )
                  {
                     R_alpha[v][j][d] = sc;
                     if ( ret_shadow != NULL ) { ((char **)R_shadow[v])[j][d] = yoffset; }
                     if ( ret_shadow != NULL ) { Rmode_shadow[v][j][d] = 1; }
                  }

                  if  ( (sc = alpha[y+yoffset][j][d] + cm->tsc[v][yoffset]) > R_alpha[v][j][d] )
                  {
                     R_alpha[v][j][d] = sc;
                     if ( ret_shadow != NULL ) { ((char **)R_shadow[v])[j][d] = yoffset; }
                     if ( ret_shadow != NULL ) { Rmode_shadow[v][j][d] = 3; }
                  }
               }

               i = j-d+1;
               if ( dsq[i] < Alphabet_size ) 
               {
                    alpha[v][j][d] += cm->esc[v][(int) dsq[i]];
                  L_alpha[v][j][d] += cm->esc[v][(int) dsq[i]];
               }
               else
               {
                    alpha[v][j][d] += DegenerateSingletScore(cm->esc[v], dsq[i]);
                  L_alpha[v][j][d] += DegenerateSingletScore(cm->esc[v], dsq[i]);
               }

               if (   alpha[v][j][d] < IMPOSSIBLE ) {   alpha[v][j][d] = IMPOSSIBLE; }
               if ( L_alpha[v][j][d] < IMPOSSIBLE ) { L_alpha[v][j][d] = IMPOSSIBLE; }
               if ( R_alpha[v][j][d] < IMPOSSIBLE ) { R_alpha[v][j][d] = IMPOSSIBLE; }
            }

            for ( d = 1; d <= jp; d++ )
            {
               if (   alpha[v][j][d] + bsc > r_sc ) { r_mode = 3; r_v = v; r_j = j; r_i = j-d+1; r_sc = bsc +   alpha[v][j][d]; }
               if ( L_alpha[v][j][d] + bsc > r_sc ) { r_mode = 2; r_v = v; r_j = j; r_i = j-d+1; r_sc = bsc + L_alpha[v][j][d]; }
            }
         }
      }
      else if ( cm->sttype[v] == IR_st || cm->sttype[v] == MR_st )
      {
         for ( jp = 0; jp <= W; jp++ )
         {
            j = i0-1+jp;
            y = cm->cfirst[v];

            alpha[v][j][0] = IMPOSSIBLE;
            L_alpha[v][j][0] = IMPOSSIBLE;
            R_alpha[v][j][0] = IMPOSSIBLE;

            for ( d = 1; d <= jp; d++ )
            {
               alpha[v][j][d]   = cm->endsc[v] + (cm->el_selfsc * (d-1));
               L_alpha[v][j][d] = IMPOSSIBLE;
               if (d == 1) R_alpha[v][j][d] = 0.0;
               else        R_alpha[v][j][d] = IMPOSSIBLE;
               if ( ret_shadow   != NULL ) { ((char **)  shadow[v])[j][d] = USED_EL; }
               if ( ret_shadow != NULL ) { ((char **)R_shadow[v])[j][d] = USED_EL; }
               if ( ret_shadow != NULL ) { Rmode_shadow[v][j][d] = 3; }

               for ( yoffset = 0; yoffset < cm->cnum[v]; yoffset++ )
               {
                  if  ( (sc = alpha[y+yoffset][j-1][d-1] + cm->tsc[v][yoffset]) > alpha[v][j][d] )
                  {
                     alpha[v][j][d] = sc;
                     if ( ret_shadow != NULL ) { ((char **)shadow[v])[j][d] = yoffset; }
                  }

                  if  ( (sc = L_alpha[y+yoffset][j][d] + cm->tsc[v][yoffset]) > L_alpha[v][j][d] )
                  {
                     L_alpha[v][j][d] = sc;
                     if ( ret_shadow != NULL ) { ((char **)L_shadow[v])[j][d] = yoffset; }
                     if ( ret_shadow != NULL ) { Lmode_shadow[v][j][d] = 2; }
                  }

                  if  ( (sc = alpha[y+yoffset][j][d] + cm->tsc[v][yoffset]) > L_alpha[v][j][d] )
                  {
                     L_alpha[v][j][d] = sc;
                     if ( ret_shadow != NULL ) { ((char **)L_shadow[v])[j][d] = yoffset; }
                     if ( ret_shadow != NULL ) { Lmode_shadow[v][j][d] = 3; }
                  }

                  if  ( d > 1 )
                  if  ( (sc = R_alpha[y+yoffset][j-1][d-1] + cm->tsc[v][yoffset]) > R_alpha[v][j][d] )
                  {
                     R_alpha[v][j][d] = sc;
                     if ( ret_shadow != NULL ) { ((char **)R_shadow[v])[j][d] = yoffset; }
                     if ( ret_shadow != NULL ) { Rmode_shadow[v][j][d] = 1; }
                  }
               }

               if ( dsq[j] < Alphabet_size ) 
               {
                    alpha[v][j][d] += cm->esc[v][(int) dsq[j]];
                  R_alpha[v][j][d] += cm->esc[v][(int) dsq[j]];
               }
               else
               {
                    alpha[v][j][d] += DegenerateSingletScore(cm->esc[v], dsq[j]);
                  R_alpha[v][j][d] += DegenerateSingletScore(cm->esc[v], dsq[j]);
               }

               if (   alpha[v][j][d] < IMPOSSIBLE ) {   alpha[v][j][d] = IMPOSSIBLE; }
               if ( L_alpha[v][j][d] < IMPOSSIBLE ) { L_alpha[v][j][d] = IMPOSSIBLE; }
               if ( R_alpha[v][j][d] < IMPOSSIBLE ) { R_alpha[v][j][d] = IMPOSSIBLE; }
            }

            for ( d = 1; d <= jp; d++ )
            {
               if (   alpha[v][j][d] + bsc > r_sc ) { r_mode = 3; r_v = v; r_j = j; r_i = j-d+1; r_sc = bsc +   alpha[v][j][d]; }
               if ( R_alpha[v][j][d] + bsc > r_sc ) { r_mode = 1; r_v = v; r_j = j; r_i = j-d+1; r_sc = bsc + R_alpha[v][j][d]; }
            }
         }
      }
      else
      {
         Die("'Inconceivable!'\n'You keep using that word...'");
      }

      if ( v==0 )
      {
           alpha[0][j0][W] = r_sc;
         L_alpha[0][j0][W] = r_sc;
         R_alpha[0][j0][W] = r_sc;
         if ( ret_shadow   != NULL ) { ((char **)  shadow[0])[j0][W] = USED_LOCAL_BEGIN; }
         if ( ret_shadow != NULL ) { ((char **)L_shadow[0])[j0][W] = USED_LOCAL_BEGIN; }
         if ( ret_shadow != NULL ) { ((char **)R_shadow[0])[j0][W] = USED_LOCAL_BEGIN; }
         if ( ret_shadow != NULL ) { Lmode_shadow[0][j0][W] = r_mode; }
         if ( ret_shadow != NULL ) { Rmode_shadow[0][j0][W] = r_mode; }
      }

      if (! do_full)
      {
         if ( cm->sttype[v] == B_st )
         {
            y = cm->cfirst[v];
            deckpool_push(dpool,   alpha[y]);   alpha[y] = NULL;
            deckpool_push(dpool, L_alpha[y]); L_alpha[y] = NULL;
            deckpool_push(dpool, R_alpha[y]); R_alpha[y] = NULL;
            z = cm->cnum[v];
            deckpool_push(dpool,   alpha[z]);   alpha[z] = NULL;
            deckpool_push(dpool, L_alpha[z]); L_alpha[z] = NULL;
            deckpool_push(dpool, R_alpha[z]); R_alpha[z] = NULL;
         }
         else
         {
            for ( y = cm->cfirst[v]; y < cm->cfirst[v]+cm->cnum[v]; y++ )
            {
               touch[y]--;
               if ( touch[y] == 0 )
               {
                  if ( cm->sttype[y] == E_st )
                  {
                     nends--;
                     if ( nends == 0 ) { deckpool_push(dpool, end); end = NULL; }
                  }
                  else
                  {
		     deckpool_push(dpool,   alpha[y]);
                     deckpool_push(dpool, L_alpha[y]);
                     deckpool_push(dpool, R_alpha[y]);
                     if ( cm->sttype[y] == B_st ) { deckpool_push(dpool, T_alpha[y]); }
                  }

                    alpha[y] = NULL;
                  L_alpha[y] = NULL;
                  R_alpha[y] = NULL;
                  T_alpha[y] = NULL;
               }
            }
         }
      }
   } /* end loop over all v */

   sc = r_sc;
   if ( ret_v     != NULL ) { *ret_v     = r_v; } 
   if ( ret_i     != NULL ) { *ret_i     = r_i; } 
   if ( ret_j     != NULL ) { *ret_j     = r_j; } 
   if ( ret_mode  != NULL ) { *ret_mode  = r_mode; }

   /* Free or return score matrices */
   if ( ret_alpha == NULL )
   {
      for ( v = vroot; v <= vend; v++ )
      {
         if ( alpha[v] != NULL )
         {
            if ( cm->sttype[v] != E_st )
            {
               deckpool_push(dpool,   alpha[v]);   alpha[v] = NULL;
               deckpool_push(dpool, L_alpha[v]); L_alpha[v] = NULL;
               deckpool_push(dpool, R_alpha[v]); R_alpha[v] = NULL;
               if ( T_alpha[v] != NULL )
               {  deckpool_push(dpool, T_alpha[v]); T_alpha[v] = NULL; }
            }
            else
            {  end = alpha[v]; }
         }
      }
      if ( end != NULL) {deckpool_push(dpool, end); end = NULL; }
      free(  alpha);
      free(L_alpha);
      free(R_alpha);
      free(T_alpha);
   }
   else
   {
      ret_alpha->J = alpha;
      ret_alpha->L = alpha->L;
      ret_alpha->R = alpha->R;
      ret_alpha->T = alpha->T;
   }

   /* Free or return deckpool */
   if ( ret_dpool == NULL )
   {
      while ( deckpool_pop(dpool, &end)) free_vjd_deck(end, i0, j0);
      deckpool_free(dpool);
   }
   else
   {
      *ret_dpool = dpool;
   }

   free(touch);
   if ( ret_shadow != NULL ) ret_shadow->J = shadow;
   if ( ret_shadow != NULL ) ret_shadow->L = L_shadow;
   if ( ret_shadow != NULL ) ret_shadow->R = R_shadow;
   if ( ret_shadow != NULL ) ret_shadow->T = T_shadow;
   if ( ret_shadow != NULL ) ret_shadow->Lmode = (void ***)Lmode_shadow;
   if ( ret_shadow != NULL ) ret_shadow->Rmode = (void ***)Rmode_shadow;
   return sc;
}

/* Function: troutside()
 * Author:   DLK
 * 
 * Purpose:  outside version of truncated CYK run on
 *           an unbifurcated model segment vroot..vend
 *           Closely based on outside()
 * Args;
 *
 * Returns:  Score of best local hit (not extending to vend)
 */
float
tr_outside(CM_t *cm, char *dsq, int L, int vroot, int vend, int i0, int j0, int do_full,
           BetaMats_t *arg_beta, BetaMats_t *ret_beta,
           struct deckpool_s *dpool, struct deckpool_s **ret_dpool,
           int *ret_mode, int *ret_v, int *ret_j)
{
   int    v,y;
   int    j,d,i;
   float  sc;
   int   *touch;
   float  esc;
   int    W;
   int    jp;
   int    voffset;
   int    w1, w2;

   float  b_sc;
   int    b_mode, b_v, b_j;

   W = j0 - i0 + 1;
   if ( dpool == NULL ) dpool = deckpool_create();

   if ( beta == NULL )
   {
      beta->J = MallocOrDie(sizeof(float **) * (cm->M+1));
      beta->L[0] = MallocOrDie(sizeof(float) * (cm->M+1)*W);
      beta->R[0] = MallocOrDie(sizeof(float) * (cm->M+1)*W);
      for ( v = 0; v < cm->M+1; v++ )
      {
         beta->J[v] = NULL;
         beta->L[v] = beta->L[0] + (v * W);
         beta->R[v] = beta->R[0] + (v * W);
      }
   }

   /* Initialize the root deck, and its split set if applicable */
   w1 = cm->nodemap[cm->ndidx[vroot]];
   if (cm->sttype[vroot] == B_st)
   {
      w2 = w1;
      if (vend != vroot) Die("vroot B but vroot != vend!\n");
   }
   else
      w2 = cm->cfirst[w1] - 1;

   for (v = w1; v <= w2; v++)
   {
      if (! deckpool_pop(dpool, &(beta->J[v])) )
         beta->J[v] = alloc_vjd_deck(L, i0, j0);
      for (jp = 0; jp <= W; jp++)
      {
         j = i0 + jp - 1;
         for (d = 0; d <= jp; d++)
            beta->J[v][j][d] = IMPOSSIBLE;
         if ( vroot == 0 && v >= vroot )
         {
         beta->L[v][j] = 0.0;
         beta->R[v][j] = 0.0;
         }
         else
         {
         beta->L[v][j] = IMPOSSIBLE;
         beta->R[v][j] = IMPOSSIBLE;
         }
      }
   }
   beta->J[vroot][j0][W] = 0.0;
   beta->L[vroot][j0][W] = 0.0;
   beta->R[vroot][j0][W] = 0.0;

   /* Initialize EL */
   if (! deckpool_pop(dpool &(beta->J[cm->M])) )
      beta->J[cm->M] = alloc_vjd_deck(L, i0, j0);
   for (jp = 0; jp <= W; jp++)
   {
      j = i0 + jp - 1;
      for (d = 0; d <= jp; d++)
         beta->J[cm->M][j][d] = IMPOSSIBLE;
      beta->L[cm->M][j] = IMPOSSIBLE;
      beta->R[cm->M][j] = IMPOSSIBLE;
   }

   /* deal with vroot->EL */
   /* Marginal modes don't transition to EL,
    * so beta->L and beta->R remain at their
    * initialization values of IMPOSSIBLE */
   if (NOT_IMPOSSIBLE(cm->endsc[vroot]))
   {
      switch (cm->sttype[vroot])
      {
         case MP_st:
            if (W < 2) break;
            if (dsq[i0] < Alphabet_size && dsq[j0] < Alphabet_size)
               esc = cm->esc[vroot][(int) (dsq[i0]*Alphabet_size+dsq[j0])];
            else
               esc = DegeneratePairScore(cm->esc[vroot], dsq[i0], dsq[j0]);
            beta->J[cm->M][j0-1][W-2] = cm->endsc[vroot] + (cm->el_selfsc * (W-2)) + esc;
            if (beta->J[cm->M][j0-1][W-2] < IMPOSSIBLE) beta->J[cm->M][j0-1][W-2] = IMPOSSIBLE;
            break;
         case ML_st:
         case IL_st:
            if (W < 1) break;
            if (dsq[i0] < Alphabet_size)
               esc = cm->esc[vroot][(int) dsq[i0]];
            else
               esc = DegenerateSingletScore(cm->esc[vroot], dsq[i0]);
            beta->J[cm->M][j0][W-1] = cm->endsc[vroot] + (cm->el_selfsc * (W-1)) + esc;
            if (beta->J[cm->M][j0][W-1] < IMPOSSIBLE) beta[cm->M][j0][W-1] = IMPOSSIBLE;
            break;
         case MR_st:
         case IR_st:
            if (W < 1) break;
            if (dsq[j0] < Alphabet_size)
               esc = cm->esc[vroot][(int) dsq[j0]];
            else
               esc = DegenerateSingletScore(cm->esc[vroot], dsq[j0]);
            beta->J[cm->M][j0-1][W-1] = cm->endsc[vroot] + (cm->el_selfsc * (W-1)) + esc;
            if (beta->J[cm->M][j0-1][W-1] < IMPOSSIBLE) beta[cm->M][j0][W-1] = IMPOSSIBLE;
            break;
         case  S_st:
         case  D_st:
            beta->J[cm->M][j0][W] = cm->endsc[vroot] + (cm->el_selfsc * W);
            if (beta->J[cm->M][j0][W] < IMPOSSIBLE) beta->J[cm->M][j0][W] = IMPOSSIBLE;
            break;
         case  B_st: /* B_st can't go to EL? */
         default:
            Die("bogus parent state %d\n",cm->sttype[vroot]);
      }
   }

   /* Initialize touch vector for controlling deck de-allocation */
   touch = MallocOrDie(sizeof(int) * cm->M);
   for (v = 0;      v < w1;    v++) touch[v] = 0;
   for (v = vend+1; v < cm->M; v++) touch[v] = 0;
   for (v = w1;     v <= vend; v++)
   {
      if (cm->sttype[v] == B_st) touch[v] = 2;
      else                       touch[v] = cm->cnum[v];
   }

   b_sc = IMPOSSIBLE;
   b_v  = -1;
   b_j  = -1;
   b_mode = -1;

   /* Main loop through decks */
   for (v = w2+1; v <= vend; v++)
   {
      /* Get a deck */
      if (! deckpool_pop(dpool, &(beta->J[v])) )
         beta->J[v] = alloc_vjd_deck(L, i0, j0);
      for (jp = W; jp >= 0; jp--)
      {
         j = i0 + jp - 1;
         for (d = jp; d >= 0; d--)
         {
            if (vroot == 0)
               beta->J[v][j][d] = 0.0;
            else
               beta->J[v][j][d] = IMPOSSIBLE;
         }
         if (vroot == 0)
         {
            beta->L[v][j] = 0.0;
            beta->R[v][j] = 0.0;
         }
         else
         {
            beta->L[v][j] = IMPOSSIBLE;
            beta->R[v][j] = IMPOSSIBLE;
         }
      }

      /* mini-recursion for beta->L */
      for (j = i0; j <= j0; j++)
      {
         for (y = cm->plast[v]; y > cm->plast[v] - cm->pnum[v]; y--)
         {
            if (y < vroot) continue;
            voffset = v - cm->cfirst[y];

            switch (cm->sttype[y])
            {
               case MP_st:
               case ML_st:
               case IL_st:
                  if (j > i0)
                  {
                     if (dsq[j-1] < Alphabet_size)
                        esc = cm->esc[y][(int) dsq[j-1]];
                     else
                        esc = DegenerateSingletScore(cm->esc[y], dsq[j-1]);
                     if ( (sc = beta->L[y][j-1] + cm->tsc[y][voffset] + esc) > beta->L[v][j] )
                        beta->L[v][j] = sc;
                  }
                  break;
               case MR_st:
               case IR_st:
               case  S_st:
               case  E_st:
               case  D_st:
                  if ( (sc = beta->L[y][j] + cm->tsc[y][voffset]) > beta->L[v][j] )
                     beta->L[v][j] = sc;
                  break;
               default:
                  Die("Bogus parent type %d for y = %d, v = %d\n",cm->sttype[y],y,v);
            }
         }
         if (beta->L[v][j] > b_sc)
         {
            b_sc = beta->L[v][j];
            b_v  = v;
            b_j  = j;
            b_mode = 2;
         }
      }

      /* mini-recursion for beta->R */
      for (j = j0; j >= i0; j--)
      {
         for (y = cm->plast[v]; y > cm->plast[v] - cm->pnum[v]; y--)
         {
            if (y < vroot) continue;
            voffset = v - cm->cfirst[y];

            switch (cm->sttype[y])
            {
               case MP_st:
               case MR_st:
               case IR_st:
                  if (j < j0)
                  {
                     if (dsq[j+1] < Alphabet_size)
                        esc = cm->esc[y][(int) dsq[j+1]];
                     else
                        esc = DegenerateSingletScore(cm->esc[y], dsq[j+1]);
                     if ( (sc = beta->R[y][j+1] + cm->tsc[y][voffset] + esc) > beta->R[v][j] )
                        beta->R[v][j] = sc;
                  }
                  break;
               case ML_st:
               case IL_st:
               case  S_st:
               case  E_st:
               case  D_st:
                  if ( (sc = beta->R[y][j] + cm->tsc[y][voffset]) > beta->L[v][j] )
                     beta->R[v][j] = sc;
                  break;
               default:
                  Die("Bogus parent type %d for y = %d, v = %d\n",cm->sttype[y],y,v);
            }
         }
         if (beta->R[v][j] > b_sc)
         {
            b_sc = beta->L[v][j];
            b_v  = v;
            b_j  = j;
            b_mode = 1;
         }
      }

      /* main recursion */
      for (jp = W; jp >= 0; jp--)
      {
         j = i0 + jp - 1;
         for (d = jp; d >= 0; d--)
         {
            i = j - d + 1;
            for (y = cm->plast[v]; y > cm->plast[v] - cm->pnum[v]; y--)
            {
               if (y < vroot) continue;
               voffset = v - cm->cfirst[y];

               switch (cm->sttype[y])
               {
                  case MP_st:
                     if (j != j0 && d != jp)
                     {
                        if (dsq[i-1] < Alphabet_size && dsq[j+1] < Alphabet_size)
                           esc = cm->esc[y][(int) (dsq[i-1]*Alphabet_size + dsq[j+1])];
                        else
                           esc = DegeneratePairScore(cm->exc[y], dsq[i-1], dsq[j+1]);
                        if ( (sc = beta->J[y][j+1][d+2] + cm->tsc[y][voffset] + esc) > beta->J[v][j][d] )
                           beta->J[v][j][d] = sc;
                        if ( (sc = beta->L[y][i-1]      + cm->tsc[y][voffset] + esc) > beta->J[v][j][d] )
                           beta->J[v][j][d] = sc;
                        if ( (sc = beta->R[y][j+1]      + cm->tsc[y][voffset] + esc) > beta->J[v][j][d] )
                           beta->J[v][j][d] = sc;
                     }
                     break;
                  case ML_st:
                  case IL_st:
                     if (d != jp)
                     {
                        if (dsq[i-1] < Alphabet_size)
                           esc = cm->esc[y][(int) dsq[i-1]];
                        else
                           esc = DegenerateSingletScore(cm->esc[y], dsq[i-1]);
                        if ( (sc = beta->J[y][j][d+1] + cm->tsc[y][voffset] + esc) > beta->J[v][j][d] )
                           beta->J[v][j][d] = sc;
                        if ( (sc = beta->R[y][j]      + cm->tsc[y][voffset] + esc) > beta->J[v][j][d] )
                           beta->J[v][j][d] = sc;
                     }
                     break;
                  case MR_st:
                  case IR_st:
                     if (j != j0)
                     {
                        if (dsq[j+1] < Alphabet_size)
                           esc = cm->esc[y][(int) dsq[j+1]];
                        else
                           esc = DegenerateSingletScore(cm->esc[y], dsq[j+1]);
                        if ( (sc = beta->J[y][j+1][d+1] + cm->tsc[y][voffset] + esc) > beta->J[v][j][d] )
                           beta->J[v][j][d] = sc;
                        if ( (sc = beta->L[y][i]        + cm->tsc[y][voffset] + esc) > beta->J[v][j][d] )
                           beta->J[v][j][d] = sc;
                     }
                     break;
                  case  S_st:
                  case  E_st:
                  case  D_st:
                     if ( (sc = beta->J[y][j][d] + cm->tsc[y][voffset]) > beta->J[v][j][d] )
                        beta->J[v][j][d] = sc;
                     break;
                  default:
                     Die("Bogus parent type %d for y = %d, v = %d\n",cm->sttype[y],y,v);
               } /* End switch over state types */
            } /* End loop over parent states  - cell done */
         } /* End loop over d - row done */
      } /* End loop over jp - deck done */

      /* v->EL transitions (beta->J only */
      if (NOT_IMPOSSIBLE(cm->endsc[v]))
      {
         for (jp = 0; jp <= W; jp++)
         {
            j = i0-1+jp;
            for (d = 0; d <= jp; d++)
            {
               i = j-d+1;
               switch (cm->sttype[v])
               {
                  case MP_st:
                     if (j == j0 || d == jp) continue; /* boundary condition */
                     if (dsq[i-1] < Alphabet_size && dsq[j+1] < Alphabet_size)
                        esc = cm->esc[v][(int) (dsq[i-1]*Alphabet_size+dsq[j+1])];
                     else
                        esc = DegeneratePairScore(cm->esc[v], dsq[i-1], dsq[j+1]);
                     if ((sc = beta[v][j+1][d+2] + cm->endsc[v] + (cm->el_selfsc * d) + esc) > beta[cm->M][j][d])
                        beta[cm->M][j][d] = sc;
                     break;
                  case ML_st:
                  case IL_st:
                     if (d == jp) continue;
                     if (dsq[i-1] < Alphabet_size)
                        esc = cm->esc[v][(int) dsq[i-1]];
                     else
                        esc = DegenerateSingletScore(cm->esc[v], dsq[i-1]);
                     if ((sc = beta[v][j][d+1] + cm->endsc[v] + (cm->el_selfsc * d) + esc) > beta[cm->M][j][d])
                        /*(cm->el_selfsc * (d+1)) + esc) > beta[cm->M][j][d])*/
                        beta[cm->M][j][d] = sc;
                     break;
                  case MR_st:
                  case IR_st:
                     if (j == j0) continue;
                     if (dsq[j+1] < Alphabet_size)
                        esc = cm->esc[v][(int) dsq[j+1]];
                     else
                        esc = DegenerateSingletScore(cm->esc[v], dsq[j+1]);
                     if ((sc = beta[v][j+1][d+1] + cm->endsc[v] + (cm->el_selfsc * d) + esc) > beta[cm->M][j][d])
                        /*(cm->el_selfsc * (d+1)) + esc) > beta[cm->M][j][d])*/
                        beta[cm->M][j][d] = sc;
                     break;
                  case S_st:
                  case D_st:
                  case E_st:
                     if ((sc = beta[v][j][d] + cm->endsc[v] + (cm->el_selfsc * d)) > beta[cm->M][j][d])
                        beta[cm->M][j][d] = sc;
                     break;
                  case B_st:
                  default: Die("bogus parent state %d\n", cm->sttype[v]);
                /* note that although B is a valid vend for a segment we'd do
                   outside on, B->EL is set to be impossible, by the local alignment
                   config. There's no point in having a B->EL because B is a nonemitter
                   (indeed, it would introduce an alignment ambiguity). The same
                   alignment case is handled by the X->EL transition where X is the
                   parent consensus state (S, MP, ML, or MR) above the B. Thus,
                   this code is relying on the NOT_IMPOSSIBLE() test, above,
                   to make sure the sttype[vend]=B case gets into this switch.
                */
               } /* end switch over parent state type v */
            } /* end inner loop over d */
         } /* end outer loop over jp */
      } /* end conditional section for dealing w/ v->EL local end transitions */

      /* Recycle memory */
      if (! do_full)
      {
         for (y = cm->plast[v]; y > cm->plast[v] - cm->pnum[v]; y--)
         {
            touch[y]--;
            if (touch[y] == 0) { deckpool_push(dpool, beta->J[y]); beta->J[y] = NULL; }
         }
      }
   } /* end loop over decks v */

   /* Clean-up */
   if (ret_beta == NULL)
   {
      for (v = w1; v <= vend; v++)
         if (beta->J[v] != NULL) { deckpool_push(dpool, beta->J[v]); beta->J[v] = NULL; }
      deckpool_push(dpool, beta->J[cm->M]); beta->J[cm->M] = NULL;
      free(beta->L);
      free(beta->R);
   }
   else
   {
      *ret_beta = beta;
   }

   if (ret_dpool == NULL)
   {
      float **a;
      while (deckpool_pop(dpool,&a)) free_vjd_deck(a, i0, j0);
      deckpool_free(dpool);
   }
   else
   {
      *ret_dpool = dpool;
   }
   free(touch);

   if (ret_mode != NULL) *ret_mode = b_mode;
   if (ret_v    != NULL) *ret_v    = b_v;
   if (ret_j    != NULL) *ret_j    = b_j;

   return b_sc;
}

/* Function: trinsideT()
 * Author:   DLK
 * 
 * Purpose:  inside with traceback on truncated sequence
 *           based on insideT()
 *
 * Args:     cm      - the covariance model
 *           dsq     - the sequence, 1..L
 *           L       - length of the sequence
 *           tr      - Parsetree for traceback
 *           r       - root of subgraph to align to target subseq (usually 0, the model's root)
 *           z       - last state of the subgraph
 *           i0      - start of target subsequence (usually 1, beginning of dsq)
 *           j0      - end of target subsequence (usually L, end of dsq)
 *           allow_begin - allow local begin (true/false)
 *           dmin    - minimum d bound for each state v (NULL if non-banded)
 *           dmax    - maximum d bound for each state v (NULL if non-banded)
 *                     In current implementation, dmin and dmax are ignored!
 *                     placeholders for future compatibility if needed.
 *
 * Returns:  score of optimal alignment (float)
 */
float
trinsideT(CM_t *cm, char *dsq, int L, Parsetree_t *tr, int r, int z,
          int i0, int j0, int allow_begin, int *dmin, int *dmax)
{
   void    ***shadow;		/* standard shadow matrix with state information */
   void    ***L_shadow;		/* left marginal shadow matrix with state information */
   void    ***R_shadow;		/* right marginal shadow matrix with state information */
   void    ***T_shadow;		/* terminal shadow matrix with state information */
   void    ***Lmode_shadow;	/* left marginal shadow matrix with mode information */
   void    ***Rmode_shadow;	/* right marginal shadow matrix with mode information */
   float      sc;		/* score of the CYK alignment */
   ESL_STACK *pda;		/* stack for storing info of 2nd child at B_st */
   int        v,i,j,d;		/* indices for state, position, & distance */
   int        mode,nxtmode;
   int        k;
   int        y, yoffset;
   int        bifparent;

   sc = trinside(cm, dsq, L, r, z, i0, j0,
                 BE_EFFICIENT,
                 &shadow,
                 &L_shadow, &R_shadow, &T_shadow,
                 &Lmode_shadow, &Rmode_shadow,
                 &mode, &v, &i, &j );
   pda = esl_stack_ICreate();
   d = j-i+1;

   while (1)
   {
      if ( cm->sttype[v] == B_st )
      {
         if      ( mode == 3 )
         {
            k = ((int **) shadow[v])[j][d];

            esl_stack_IPush(pda, j);
            esl_stack_IPush(pda, k);
            esl_stack_IPush(pda, mode);
            esl_stack_IPush(pda, tr->n-1);
         }
         else if ( mode == 2 )
         {
            k = ((int **) L_shadow[v])[j][d];

            /* In left marginal mode, right child is always left marginal */
            esl_stack_IPush(pda, j);
            esl_stack_IPush(pda, k);
            esl_stack_IPush(pda, mode);
            esl_stack_IPush(pda, tr->n-1);

            /* Retrieve mode of left child (should be 3 or 2) */
            mode = ((int **)Lmode_shadow[v])[j][d];
         }
         else if ( mode == 1 )
         {
            k = ((int **) R_shadow[v])[j][d];

            /* Retrieve mode of right child (should be 3 or 1) */
            mode = ((int **)Rmode_shadow[v])[j][d];
            esl_stack_IPush(pda, j);
            esl_stack_IPush(pda, k);
            esl_stack_IPush(pda, mode);
            esl_stack_IPush(pda, tr->n-1);

            /* In right marginal mode, left child is always right marginal */
            mode = 1;
         }
         else if ( mode == 0 )
         {
            k = ((int **) T_shadow[v])[j][d];

             mode = 2;
             esl_stack_IPush(pda, j);
             esl_stack_IPush(pda, k);
             esl_stack_IPush(pda, mode);
             esl_stack_IPush(pda, tr->n-1);

             mode = 1;
         }
         else { Die("Unknown mode in traceback!"); }

         j = j-k;
         d = d-k;
         i = j-d+1;
         v = cm->cfirst[v];
         InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, v, mode);
      }
      else if ( (cm->sttype[v] == E_st) || (cm->sttype[v] == EL_st) )
      {
         if (esl_stack_IPop(pda, &bifparent) == eslEOD) break;
         esl_stack_IPop(pda, &mode);
         esl_stack_IPop(pda, &d);
         esl_stack_IPop(pda, &j);
         v = tr->state[bifparent];
         y = cm->cnum[v];
         i = j-d+1;

         v = y;
         InsertTraceNodewithMode(tr, bifparent, TRACE_RIGHT_CHILD, i, j, v, mode);
      }
      else
      {
         if      ( mode == 3 )
         {
            yoffset = ((char **)   shadow[v])[j][d];
            nxtmode = 3;
         }
         else if ( mode == 2 )
         {
            yoffset = ((char **) L_shadow[v])[j][d];
            nxtmode = ((int **)Lmode_shadow[v])[j][d];
         }
         else if ( mode == 1 )
         {
            yoffset = ((char **) R_shadow[v])[j][d];
            nxtmode = ((int **)Rmode_shadow[v])[j][d];
         }
         else { Die("Unknown mode in traceback!"); }

         switch (cm->sttype[v])
         {
            case  D_st:
               break;
            case MP_st:
               if ( mode == 3 )          i++;
               if ( mode == 2 && d > 0 ) i++;
               if ( mode == 3 )          j--;
               if ( mode == 1 && d > 0 ) j--;
               break;
            case ML_st:
               if ( mode == 3 )          i++;
               if ( mode == 2 && d > 0 ) i++;
               break;
            case MR_st:
               if ( mode == 3 )          j--;
               if ( mode == 1 && d > 0 ) j--;
               break;
            case IL_st:
               if ( mode == 3 )          i++;
               if ( mode == 2 && d > 0 ) i++;
               break;
            case IR_st:
               if ( mode == 3 )          j--;
               if ( mode == 1 && d > 0 ) j--;
               break;
            case  S_st:
               break;
            default:
               Die("'Inconceivable!'\n'You keep using that word...'");
         }
         d = j-i+1;
         mode = nxtmode;

         if ( yoffset == USED_EL )
         {
            v = cm->M;
            InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, v, mode);
         }
         else if ( yoffset == USED_LOCAL_BEGIN )
         {  /* local begin, can only happen once, from root */
            /* However, all hits from truncyk() are local hits, and this should have
               been dealt with immediately after return from the DP function.
               If we've reached this point, there's a major problem */
            Die("Impossible local begin in traceback\n");
         }
         else
         {
            v = cm->cfirst[v] + yoffset;
            InsertTraceNodewithMode(tr, tr->n-1, TRACE_LEFT_CHILD, i, j, v, mode);
         }
      }
   }

   esl_stack_Destroy(pda);
   free_vjd_shadow_matrix(shadow, cm, i0, j0);
   free_vjd_shadow_matrix(L_shadow, cm, i0, j0);
   free_vjd_shadow_matrix(R_shadow, cm, i0, j0);
   free_vjd_shadow_matrix(T_shadow, cm, i0, j0);
   free_vjd_shadow_matrix(Lmode_shadow, cm, i0, j0);
   free_vjd_shadow_matrix(Rmode_shadow, cm, i0, j0);

   return sc;
}

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

/* Alignment engine */
float trinside (CM_t *cm, ESL_SQ *sq, int vroot, int vend, int i0, int j0, int do_full,
                void ****ret_shadow, void ****ret_L_shadow, void ****ret_R_shadow,
                void ****ret_T_shadow, void ****ret_Lmode_shadow, void ****ret_Rmode_shadow,
                int *ret_mode, int *ret_v, int *ret_i, int *ret_j);


/* Traceback routine */
float trinsideT(CM_t *cm, ESL_SQ *sq, Parsetree_t *tr, int r, int z,
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
TrCYKInside(CM_t *cm, ESL_SQ *sq, int r, int i0, int j0, Parsetree_t **ret_tr,
            int *dmin, int *dmax)
{
  if(! (sq->flags & eslSQ_DIGITAL))
    esl_fatal("ERROR in TrCYKInside(), sq is not digitized.\n");

   Parsetree_t *tr;
   int          z;
   float        sc;
  
   /* Check input parameters */
   if ( cm->stid[r] != ROOT_S )
   {
      if (! (cm->flags & CM_LOCAL_BEGIN)) esl_fatal("internal error: we're not in local mode, but r is not root");
      if ( (cm->stid[r] != MATP_MP) &&
           (cm->stid[r] != MATL_ML) &&
           (cm->stid[r] != MATR_MR) &&
           (cm->stid[r] != BIF_B  )    )  esl_fatal("internal error: trying to do a local begin at a non-mainline start");
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
   sc += trinsideT(cm, sq, tr, r, z, i0, j0, (r==0), dmin, dmax);

   if ( ret_tr != NULL ) *ret_tr = tr; else FreeParsetree(tr);

   return sc;
}

/* Function: trinside()
 * Author:   DLK
 *
 * Purpose:  inside phase of CYK on truncated sequence
 *           based on inside()
 *
 *           Score matrices and deckpool are only managed
 *           within this func, not available to caller
 *
 * Args:     cm      - the covariance model
 *           sq      - the sequence, 1..sq->n
 *           r       - root of subgraph to align to target subseq (usually 0, the model's root)
 *           z       - last state of the subgraph
 *           i0      - start of target subsequence (usually 1, beginning of sq)
 *           j0      - end of target subsequence (usually L, end of sq)
 *           do_full - if TRUE, save all decks rather than re-using
 * 
 *           allow_begin - TRUE to allow 0->b local alignment begin transitions
 *           ret_b       - state for best local begin, NULL if unwanted
 *           ret_bsc     - score for best local begin, NULL if unwanted
 *           ret_bmode   - mode  for best local begin, NULL if unwanted
 *
 * Returns:  Score of the optimal alignment
 */
float
trinside (CM_t *cm, ESL_SQ *sq, int vroot, int vend, int i0, int j0, int do_full,
          void ****ret_shadow, void ****ret_L_shadow, void ****ret_R_shadow,
          void ****ret_T_shadow, void ****ret_Lmode_shadow, void ****ret_Rmode_shadow,
          int *ret_mode, int *ret_v, int *ret_i, int *ret_j)
{
   float  **end;
   int      nends;
   int     *touch;
   int      v,y,z;
   int      j,d,i,k;
   float    sc,tsc;
   /* float esc; */  /* removed to silence a compiler warning: it's unused at present (SRE) */
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

   struct deckpool_s *dpool = NULL;
   float ***alpha   = NULL;
   float ***L_alpha = NULL;
   float ***R_alpha = NULL;
   float ***T_alpha = NULL;

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
   /* 2.0 instead of 2 to force floating point division, not integer division */
   bsc = sreLOG2(2.0/(model_len*(model_len+1)));

   /* Make a deckpool */
   if ( dpool == NULL ) dpool = deckpool_create();
   if (! deckpool_pop(dpool, &end) )
   {  end = alloc_vjd_deck(sq->n, i0, j0); }
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
   }
   if ( ret_L_shadow != NULL ) 
   {
      L_shadow = (void ***) MallocOrDie(sizeof(void **) * cm->M);
      for ( v=0; v<cm->M; v++ ) { L_shadow[v] = NULL; }
   }
   if ( ret_R_shadow != NULL ) 
   {
      R_shadow = (void ***) MallocOrDie(sizeof(void **) * cm->M);
      for ( v=0; v<cm->M; v++ ) { R_shadow[v] = NULL; }
   }
   if ( ret_T_shadow != NULL ) 
   {
      T_shadow = (void ***) MallocOrDie(sizeof(void **) * cm->M);
      for ( v=0; v<cm->M; v++ ) { T_shadow[v] = NULL; }
   }
   if ( ret_Lmode_shadow != NULL ) 
   {
      Lmode_shadow = (int ***) MallocOrDie(sizeof(int **) * cm->M);
      for ( v=0; v<cm->M; v++ ) { Lmode_shadow[v] = NULL; }
   }
   if ( ret_Rmode_shadow != NULL ) 
   {
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
      {  alpha[v] = alloc_vjd_deck(sq->n, i0, j0); }
      if (! deckpool_pop(dpool, &(L_alpha[v])) )
      {  L_alpha[v] = alloc_vjd_deck(sq->n, i0, j0); }
      if (! deckpool_pop(dpool, &(R_alpha[v])) )
      {  R_alpha[v] = alloc_vjd_deck(sq->n, i0, j0); }
      if ( (cm->sttype[v] == B_st) && (! deckpool_pop(dpool, &(T_alpha[v])) ) )
      {  T_alpha[v] = alloc_vjd_deck(sq->n, i0, j0); }

      /* Assign shadow decks */
      if ( ret_shadow != NULL )
      {
         if ( cm->sttype[v] == B_st )
         {
            kshad = alloc_vjd_kshadow_deck(sq->n, i0, j0);
            shadow[v] = (void **) kshad;
         }
         else
         {
            yshad = alloc_vjd_yshadow_deck(sq->n, i0, j0);
            shadow[v] = (void **) yshad;
         }
      }
      if ( ret_L_shadow != NULL )
      {
         if ( cm->sttype[v] == B_st )
         {
            kshad = alloc_vjd_kshadow_deck(sq->n, i0, j0);
            L_shadow[v] = (void **) kshad;
         }
         else
         {
            yshad = alloc_vjd_yshadow_deck(sq->n, i0, j0);
            L_shadow[v] = (void **) yshad;
         }
      }
      if ( ret_R_shadow != NULL )
      {
         if ( cm->sttype[v] == B_st )
         {
            kshad = alloc_vjd_kshadow_deck(sq->n, i0, j0);
            R_shadow[v] = (void **) kshad;
         }
         else
         {
            yshad = alloc_vjd_yshadow_deck(sq->n, i0, j0);
            R_shadow[v] = (void **) yshad;
         }
      }
      if ( ret_T_shadow != NULL )
      {
         if ( cm->sttype[v] == B_st )
         {
            kshad = alloc_vjd_kshadow_deck(sq->n, i0, j0);
            T_shadow[v] = (void **) kshad;
         }
      }
      if ( ret_Lmode_shadow != NULL )
      {
         kshad = alloc_vjd_kshadow_deck(sq->n, i0, j0);
         Lmode_shadow[v] = (int **) kshad;
      }
      if ( ret_Rmode_shadow != NULL )
      {
         kshad = alloc_vjd_kshadow_deck(sq->n, i0, j0);
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
               if (cm->sttype[v] == S_st) alpha[v][j][d]   = cm->endsc[v] + (cm->el_selfsc * (d));
               else                       alpha[v][j][d]   = IMPOSSIBLE;
               L_alpha[v][j][d] = IMPOSSIBLE;
               R_alpha[v][j][d] = IMPOSSIBLE;
               if ( ret_shadow   != NULL ) { ((char **)  shadow[v])[j][d] = USED_EL; }

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
                     if ( ret_L_shadow != NULL ) ((char **)L_shadow[v])[j][d] = yoffset;
                     if ( ret_Lmode_shadow != NULL ) Lmode_shadow[v][j][d] = 2;
                  }
                  if ( (sc = R_alpha[y+yoffset][j][d] + cm->tsc[v][yoffset]) > R_alpha[v][j][d] )
                  {
                     R_alpha[v][j][d] = sc;
                     if ( ret_R_shadow != NULL ) ((char **)R_shadow[v])[j][d] = yoffset;
                     if ( ret_Rmode_shadow != NULL ) Rmode_shadow[v][j][d] = 1;
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
               if ( ret_shadow   != NULL ) { ((int **)  shadow[v])[j][d] = 0; }
               if ( ret_L_shadow != NULL ) { ((int **)L_shadow[v])[j][d] = 0; }
               if ( ret_R_shadow != NULL ) { ((int **)R_shadow[v])[j][d] = d; }
               if ( ret_Lmode_shadow != NULL) { Lmode_shadow[v][j][d] = 2; }
               if ( ret_Rmode_shadow != NULL) { Rmode_shadow[v][j][d] = 1; }

               if  ( (sc = alpha[y][j][d]                   ) > L_alpha[v][j][d] )
               {
                  L_alpha[v][j][d] = sc;
                  if ( ret_L_shadow != NULL ) { ((int **)L_shadow[v])[j][d] = 0; }
                  if ( ret_Lmode_shadow != NULL ) { Lmode_shadow[v][j][d] = 3; }
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
                     if ( ret_L_shadow != NULL ) { ((int **)L_shadow[v])[j][d] = k; }
                     if ( ret_Lmode_shadow != NULL ) { Lmode_shadow[v][j][d] = 3; }
                     allow_L_exit = 1;
                  }
               }
               if ( (sc =                        alpha[z][j][d]) > R_alpha[v][j][d] )
               {
                  R_alpha[v][j][d] = sc;
                  if ( ret_R_shadow != NULL ) { ((int **)R_shadow[v])[j][d] = d; }
                  if ( ret_Rmode_shadow != NULL ) { Rmode_shadow[v][j][d] = 3; }
               }
               for ( k=0; k< d; k++ )
               {
                  if ( (sc = R_alpha[y][j-k][d-k] + alpha[z][j][k]) > R_alpha[v][j][d] )
                  {
                     R_alpha[v][j][d] = sc;
                     if ( ret_R_shadow != NULL ) { ((int **)R_shadow[v])[j][d] = k; }
                     if ( ret_Rmode_shadow != NULL ) { Rmode_shadow[v][j][d] = 3; }
                     allow_R_exit = 1;
                  }
               }

               if (d == 0) {
                  L_alpha[v][j][d] = IMPOSSIBLE;
                  R_alpha[v][j][d] = IMPOSSIBLE;
               }

               if (d >= 2) {
                 T_alpha[v][j][d] = R_alpha[y][j-1][d-1] + L_alpha[z][j][1];
                 if ( ret_T_shadow != NULL ) { ((int **)T_shadow[v])[j][d] = 1; }
                 for ( k=2; k<d; k++ )
                 {
                    if ( (sc = R_alpha[y][j-k][d-k] + L_alpha[z][j][k]) > T_alpha[v][j][d] )
                    {
                       T_alpha[v][j][d] = sc;
                       if ( ret_T_shadow != NULL) { ((int **)T_shadow[v])[j][d] = k; }
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
            if ( jp > 0 ) { alpha[v][j][1] = IMPOSSIBLE; }
            if ( sq->dsq[j] < cm->abc->K ) { L_alpha[v][j][1] =  LeftMarginalScore(cm->abc, cm->esc[v],sq->dsq[j]); }
            else { esl_fatal("Still can't deal with marginalizing degenerate residues!"); }
            if ( sq->dsq[j] < cm->abc->K ) { R_alpha[v][j][1] = RightMarginalScore(cm->abc, cm->esc[v],sq->dsq[j]); }
            else { esl_fatal("Still can't deal with marginalizing degenerate residues!"); }
            if ( ret_L_shadow != NULL ) { ((char **)L_shadow[v])[j][d] = USED_EL; }
            if ( ret_R_shadow != NULL ) { ((char **)R_shadow[v])[j][d] = USED_EL; }
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
                     if ( ret_L_shadow != NULL ) { ((char **)L_shadow[v])[j][d] = yoffset; }
                     if ( ret_Lmode_shadow != NULL ) { Lmode_shadow[v][j][d] = 2; }
                  }
                  if ( (sc = alpha[y+yoffset][j][d-1] + cm->tsc[v][yoffset]) > L_alpha[v][j][d])
                  {
                     L_alpha[v][j][d] = sc;
                     if ( ret_L_shadow != NULL ) { ((char **)L_shadow[v])[j][d] = yoffset; }
                     if ( ret_Lmode_shadow != NULL ) { Lmode_shadow[v][j][d] = 3; }
                  }

                  if ( (sc = R_alpha[y+yoffset][j-1][d-1] + cm->tsc[v][yoffset]) > R_alpha[v][j][d])
                  {
                     R_alpha[v][j][d] = sc;
                     if ( ret_R_shadow != NULL ) { ((char **)R_shadow[v])[j][d] = yoffset; }
                     if ( ret_Rmode_shadow != NULL ) { Rmode_shadow[v][j][d] = 1; }
                  }
                  if ( (sc = alpha[y+yoffset][j-1][d-1] + cm->tsc[v][yoffset]) > R_alpha[v][j][d])
                  {
                     R_alpha[v][j][d] = sc;
                     if ( ret_R_shadow != NULL ) { ((char **)R_shadow[v])[j][d] = yoffset; }
                     if ( ret_Rmode_shadow != NULL ) { Rmode_shadow[v][j][d] = 3; }
                  }
               }

               i = j-d+1;
               if ( sq->dsq[i] < cm->abc->K && sq->dsq[j] < cm->abc->K )
               {  alpha[v][j][d] += cm->esc[v][(int) (sq->dsq[i]*cm->abc->K+sq->dsq[j])]; }
               else
               {  alpha[v][j][d] += DegeneratePairScore(cm->abc, cm->esc[v], sq->dsq[i], sq->dsq[j]); }
               if ( sq->dsq[i] < cm->abc->K )
               { L_alpha[v][j][d] +=  LeftMarginalScore(cm->abc, cm->esc[v],sq->dsq[i]); }
               else
               { esl_fatal("Still can't deal with marginalizing degenerate residues!"); }
               if ( sq->dsq[j] < cm->abc->K )
               { R_alpha[v][j][d] += RightMarginalScore(cm->abc, cm->esc[v],sq->dsq[j]); }
               else
               { esl_fatal("Still can't deal with marginalizing degenerate residues!"); }

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
               if ( ret_L_shadow != NULL ) { ((char **)L_shadow[v])[j][d] = USED_EL; }
               if ( ret_Lmode_shadow != NULL ) { Lmode_shadow[v][j][d] = 3; }

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
                     if ( ret_L_shadow != NULL ) { ((char **)L_shadow[v])[j][d] = yoffset; }
                     if ( ret_Lmode_shadow != NULL ) { Lmode_shadow[v][j][d] = 2; }
                  }

                  if  ( (sc = R_alpha[y+yoffset][j][d] + cm->tsc[v][yoffset]) > R_alpha[v][j][d] )
                  {
                     R_alpha[v][j][d] = sc;
                     if ( ret_R_shadow != NULL ) { ((char **)R_shadow[v])[j][d] = yoffset; }
                     if ( ret_Rmode_shadow != NULL ) { Rmode_shadow[v][j][d] = 1; }
                  }

                  if  ( (sc = alpha[y+yoffset][j][d] + cm->tsc[v][yoffset]) > R_alpha[v][j][d] )
                  {
                     R_alpha[v][j][d] = sc;
                     if ( ret_R_shadow != NULL ) { ((char **)R_shadow[v])[j][d] = yoffset; }
                     if ( ret_Rmode_shadow != NULL ) { Rmode_shadow[v][j][d] = 3; }
                  }
               }

               i = j-d+1;
               if ( sq->dsq[i] < cm->abc->K ) 
               {
                    alpha[v][j][d] += cm->esc[v][(int) sq->dsq[i]];
                  L_alpha[v][j][d] += cm->esc[v][(int) sq->dsq[i]];
               }
               else
               {
                    alpha[v][j][d] += esl_abc_FAvgScore(cm->abc, sq->dsq[i], cm->esc[v]);
                  L_alpha[v][j][d] += esl_abc_FAvgScore(cm->abc, sq->dsq[i], cm->esc[v]);
		  
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
               if ( ret_R_shadow != NULL ) { ((char **)R_shadow[v])[j][d] = USED_EL; }
               if ( ret_Rmode_shadow != NULL ) { Rmode_shadow[v][j][d] = 3; }

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
                     if ( ret_L_shadow != NULL ) { ((char **)L_shadow[v])[j][d] = yoffset; }
                     if ( ret_Lmode_shadow != NULL ) { Lmode_shadow[v][j][d] = 2; }
                  }

                  if  ( (sc = alpha[y+yoffset][j][d] + cm->tsc[v][yoffset]) > L_alpha[v][j][d] )
                  {
                     L_alpha[v][j][d] = sc;
                     if ( ret_L_shadow != NULL ) { ((char **)L_shadow[v])[j][d] = yoffset; }
                     if ( ret_Lmode_shadow != NULL ) { Lmode_shadow[v][j][d] = 3; }
                  }

                  if  ( d > 1 )
                  if  ( (sc = R_alpha[y+yoffset][j-1][d-1] + cm->tsc[v][yoffset]) > R_alpha[v][j][d] )
                  {
                     R_alpha[v][j][d] = sc;
                     if ( ret_R_shadow != NULL ) { ((char **)R_shadow[v])[j][d] = yoffset; }
                     if ( ret_Rmode_shadow != NULL ) { Rmode_shadow[v][j][d] = 1; }
                  }
               }

               if ( sq->dsq[j] < cm->abc->K ) 
               {
                    alpha[v][j][d] += cm->esc[v][(int) sq->dsq[j]];
                  R_alpha[v][j][d] += cm->esc[v][(int) sq->dsq[j]];
               }
               else
               {
                    alpha[v][j][d] += esl_abc_FAvgScore(cm->abc, sq->dsq[j], cm->esc[v]);
                  R_alpha[v][j][d] += esl_abc_FAvgScore(cm->abc, sq->dsq[j], cm->esc[v]);
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
         esl_fatal("'Inconceivable!'\n'You keep using that word...'");
      }

      if ( v==0 )
      {
           alpha[0][j0][W] = r_sc;
         L_alpha[0][j0][W] = r_sc;
         R_alpha[0][j0][W] = r_sc;
         if ( ret_shadow   != NULL ) { ((char **)  shadow[0])[j0][W] = USED_LOCAL_BEGIN; }
         if ( ret_L_shadow != NULL ) { ((char **)L_shadow[0])[j0][W] = USED_LOCAL_BEGIN; }
         if ( ret_R_shadow != NULL ) { ((char **)R_shadow[0])[j0][W] = USED_LOCAL_BEGIN; }
         if ( ret_Lmode_shadow != NULL ) { Lmode_shadow[0][j0][W] = r_mode; }
         if ( ret_Rmode_shadow != NULL ) { Rmode_shadow[0][j0][W] = r_mode; }
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

   /* No option for returning score matrices - delete them all */
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

   /* No option for returning deckpool - delete it */
   while ( deckpool_pop(dpool, &end)) free_vjd_deck(end, i0, j0);
   deckpool_free(dpool);

   free(touch);
   if ( ret_shadow != NULL ) *ret_shadow = shadow;
   if ( ret_L_shadow != NULL ) *ret_L_shadow = L_shadow;
   if ( ret_R_shadow != NULL ) *ret_R_shadow = R_shadow;
   if ( ret_T_shadow != NULL ) *ret_T_shadow = T_shadow;
   if ( ret_Lmode_shadow != NULL ) *ret_Lmode_shadow = (void ***)Lmode_shadow;
   if ( ret_Rmode_shadow != NULL ) *ret_Rmode_shadow = (void ***)Rmode_shadow;
   return sc;
}

/* Function: trinsideT()
 * Author:   DLK
 * 
 * Purpose:  inside with traceback on truncated sequence
 *           based on insideT()
 *
 * Args:     cm      - the covariance model
 *           sq      - the sequence, 1..sq->n
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
trinsideT(CM_t *cm, ESL_SQ *sq, Parsetree_t *tr, int r, int z,
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

   sc = trinside(cm, sq, r, z, i0, j0,
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
         else { esl_fatal("Unknown mode in traceback!"); }

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
         else { esl_fatal("Unknown mode in traceback!"); }

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
               esl_fatal("'Inconceivable!'\n'You keep using that word...'");
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
            esl_fatal("Impossible local begin in traceback\n");
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

#include "squid.h"
#include "msa.h"
#include "structs.h"

/* from alphabet.c
 */
extern char SymbolIndex(char sym);
extern void SingletCount(float *counters, char symidx, float wt);
extern void PairCount(float *counters, char syml, char symr, float wt);

/* from cm.c
 */
extern CM_t *CreateCM(int nnodes, int nstates);
extern void  CMZero(CM_t *cm);
extern void  FreeCM(CM_t *cm);
extern void  CMSetDefaultNullModel(CM_t *cm);
extern void  CMSimpleProbify(CM_t *cm);
extern void  CMLogoddsify(CM_t *cm);
extern int   CalculateStateIndex(CM_t *cm, int node, char utype);
extern void  PrintCM(FILE *fp, CM_t *cm);
extern void  SummarizeCM(FILE *fp, CM_t *cm);
extern char *Statetype(int type);
extern char *Nodetype(int type);
extern char *UniqueStatetype(int type);

/* from modelmaker.c
 */
extern void HandModelmaker(MSA *msa, int use_rf, float gapthresh, CM_t **ret_cm, Parsetree_t **ret_mtr);

/* from parsetree.c
 */
extern Parsetree_t *CreateParsetree(void);
extern void         GrowParsetree(Parsetree_t *tr);
extern void         FreeParsetree(Parsetree_t *tr);
extern int          InsertTraceNode(Parsetree_t *tr, int y, int whichway, 
				    int emitl, int emitr, int state);
extern void         ParsetreeCount(CM_t *cm, Parsetree_t *tr, char *seq, float wgt);
extern void         PrintParsetree(FILE *fp, Parsetree_t *tr);
extern void         SummarizeMasterTrace(FILE *fp, Parsetree_t *tr);

/* from rna_ops.c
 */
extern int KHS2ct(char *ss, int len, int allow_pseudoknots, int **ret_ct);

#include "squid.h"
#include "msa.h"
#include "structs.h"

/* from alphabet.c
 */
extern char   SymbolIndex(char sym);
extern void   SingletCount(float *counters, char symidx, float wt);
extern void   PairCount(float *counters, char syml, char symr, float wt);
extern float  DegeneratePairScore(float *esc, char syml, char symr);
extern float  DegenerateSingletScore(float *esc, char sym);
extern char  *DigitizeSequence(char *seq, int L);
extern char **DigitizeAlignment(char **aseq, int nseq, int alen);

/* from cm.c
 */
extern CM_t *CreateCM(int nnodes, int nstates);
extern void  CMZero(CM_t *cm);
extern void  FreeCM(CM_t *cm);
extern void  CMSetDefaultNullModel(CM_t *cm);
extern void  CMSimpleProbify(CM_t *cm);
extern void  CMLogoddsify(CM_t *cm);
extern int   CMCountStatetype(CM_t *cm, char type);
extern int   CMSegmentCountStatetype(CM_t *cm, int r, int z, char type);
extern int   CMSubtreeCountStatetype(CM_t *cm, int v, char type);
extern int   CMSubtreeFindEnd(CM_t *cm, int v);
extern int   CalculateStateIndex(CM_t *cm, int node, char utype);
extern void  PrintCM(FILE *fp, CM_t *cm);
extern void  SummarizeCM(FILE *fp, CM_t *cm);
extern char *Statetype(int type);
extern char *Nodetype(int type);
extern char *UniqueStatetype(int type);
extern CM_t *CMRebalance(CM_t *cm);

/* from cmio.c
 */
extern void WriteBinaryCM(FILE *fp, CM_t *cm);
extern int  ReadBinaryCM(FILE *fp, CM_t **ret_cm);

/* from display.c
 */
extern Fancyali_t    *CreateFancyAli(Parsetree_t *tr, CM_t *cm, CMConsensus_t *cons, char *dsq);
extern void           PrintFancyAli(FILE *fp, Fancyali_t *ali);
extern void           FreeFancyAli(Fancyali_t *ali);
extern CMConsensus_t *CreateCMConsensus(CM_t *cm, float pthresh, float sthresh);
extern void           FreeCMConsensus(CMConsensus_t *con);

/* from modelconfig.c
 */
extern void ConfigLocal(CM_t *cm, float p_internal_start, float p_internal_exit);

/* from modelmaker.c
 */
extern void HandModelmaker(MSA *msa, char **dsq, int use_rf, float gapthresh, 
			   CM_t **ret_cm, Parsetree_t **ret_mtr, Parsetree_t ***ret_tr);

/* from parsetree.c
 */
extern Parsetree_t *CreateParsetree(void);
extern void         GrowParsetree(Parsetree_t *tr);
extern void         FreeParsetree(Parsetree_t *tr);
extern int          InsertTraceNode(Parsetree_t *tr, int y, int whichway, 
				    int emitl, int emitr, int state);
extern void         ParsetreeCount(CM_t *cm, Parsetree_t *tr, char *seq, float wgt);
extern float        ParsetreeScore(CM_t *cm, Parsetree_t *tr, char *dsq);
extern void         PrintParsetree(FILE *fp, Parsetree_t *tr);
extern void         ParsetreeDump(FILE *fp, Parsetree_t *tr, CM_t *cm, char *dsq);
extern void         ParsetreeCompare(Parsetree_t *t1, Parsetree_t *t2);
extern void         SummarizeMasterTrace(FILE *fp, Parsetree_t *tr);
extern void         MasterTraceDisplay(FILE *fp, Parsetree_t *mtr, CM_t *cm);


/* from rna_ops.c
 */
extern int KHS2ct(char *ss, int len, int allow_pseudoknots, int **ret_ct);

/* from scancyk.c
 */
extern void  CYKScan(CM_t *cm, char *dsq, int L, int W, 
		     int *ret_nhits, int **ret_hitr, 
		     int **ret_hiti, int **ret_hitj, float **ret_hitsc);
extern float CYKScanRequires(CM_t *cm, int L, int W);

/* from smallcyk.c
 */
extern float CYKDivideAndConquer(CM_t *cm, char *dsq, int L, Parsetree_t **ret_tr);
extern float CYKLocalDivideAndConquer(CM_t *cm, char *dsq, int L, int r0, int i0, int j0,
				      Parsetree_t **ret_tr);
extern float CYKInside(CM_t *cm, char *dsq, int L, Parsetree_t **ret_tr);
extern float CYKInsideScore(CM_t *cm, char *dsq, int L);
extern void  CYKDemands(CM_t *cm, int L);
extern int   CYKDeckCount(CM_t *cm);


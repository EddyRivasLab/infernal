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

/* from bandcyk.c
 */
extern void     BandExperiment(CM_t *cm);
extern double **BandDistribution(CM_t *cm, int W);
extern void     BandBounds(double **gamma, int M, int W, double p, 
			   int **ret_min, int **ret_max);
extern void     PrintBandGraph(FILE *fp, double **gamma, int *min, int *max, int v, int W);

extern void     PrintDPCellsSaved(CM_t *cm, int *min, int *max, int W);
extern void     CYKBandedScan(CM_t *cm, char *dsq, int *dmin, int *dmax, int L, int W, 
			      int *ret_nhits, int **ret_hitr, 
			      int **ret_hiti, int **ret_hitj, float **ret_hitsc);
extern void     BandedParsetreeDump(FILE *fp, Parsetree_t *tr, CM_t *cm, char *dsq, 
				    double **gamma, int W, int *dmin, int *dmax);

/* from cm.c
 */
extern CM_t *CreateCM(int nnodes, int nstates);
extern CM_t *CreateCMShell(void);
extern void  CreateCMBody(CM_t *cm, int nnodes, int nstates);
extern void  CMZero(CM_t *cm);
extern void  CMRenormalize(CM_t *cm);
extern void  FreeCM(CM_t *cm);
extern void  CMSetDefaultNullModel(CM_t *cm);
extern void  CMSimpleProbify(CM_t *cm);
extern void  CMLogoddsify(CM_t *cm);
extern void  CMHackInsertScores(CM_t *cm);
extern int   CMCountStatetype(CM_t *cm, char type);
extern int   CMSegmentCountStatetype(CM_t *cm, int r, int z, char type);
extern int   CMSubtreeCountStatetype(CM_t *cm, int v, char type);
extern int   CMSubtreeFindEnd(CM_t *cm, int v);
extern int   CalculateStateIndex(CM_t *cm, int node, char utype);
extern int   TotalStatesInNode(int ndtype);
extern int   SplitStatesInNode(int ndtype);
extern int   InsertStatesInNode(int ndtype);
extern int   StateDelta(int sttype);
extern int   StateLeftDelta(int sttype);
extern int   StateRightDelta(int sttype);
extern void  PrintCM(FILE *fp, CM_t *cm);
extern void  SummarizeCM(FILE *fp, CM_t *cm);
extern char *Statetype(int type);
extern int   StateCode(char *s);
extern char *Nodetype(int type);
extern int   NodeCode(char *s);
extern char *UniqueStatetype(int type);
extern int   UniqueStateCode(char *s);
extern int   DeriveUniqueStateCode(int ndtype, int sttype);
extern CM_t *CMRebalance(CM_t *cm);

/* from cmio.c
 */
extern CMFILE *CMFileOpen(char *cmfile, char *env);
extern int     CMFileRead(CMFILE *cmf, CM_t **ret_cm);
extern void    CMFileClose(CMFILE *cmf);
extern void    CMFileRewind(CMFILE *cmf);
extern int     CMFilePositionByIndex(CMFILE *cmf, int idx);
extern int     CMFilePositionByKey(CMFILE *cmf, char *key);
extern void    CMFileWrite(FILE *fp, CM_t *cm, int do_binary);


/* from display.c
 */
extern Fancyali_t    *CreateFancyAli(Parsetree_t *tr, CM_t *cm, CMConsensus_t *cons, char *dsq);
extern void           PrintFancyAli(FILE *fp, Fancyali_t *ali);
extern void           FreeFancyAli(Fancyali_t *ali);
extern CMConsensus_t *CreateCMConsensus(CM_t *cm, float pthresh, float sthresh);
extern void           FreeCMConsensus(CMConsensus_t *con);
extern void           MainBanner(FILE *fp, char *banner); 

/* in emit.c
 */
extern void EmitParsetree(CM_t *cm, Parsetree_t **ret_tr, char **ret_seq, 
			  char **ret_dsq, int *ret_N);

/* in emitmap.c
 */
extern CMEmitMap_t *CreateEmitMap(CM_t *cm); 
extern void         DumpEmitMap(FILE *fp, CMEmitMap_t *map, CM_t *cm);
extern void         FreeEmitMap(CMEmitMap_t *map);

/* from modelconfig.c
 */
extern void ConfigLocal(CM_t *cm, float p_internal_start, float p_internal_exit);

/* from modelmaker.c
 */
extern void HandModelmaker(MSA *msa, char **dsq, int use_rf, float gapthresh, 
			   CM_t **ret_cm, Parsetree_t **ret_mtr);
extern Parsetree_t *Transmogrify(CM_t *cm, Parsetree_t *gtr, 
				 char *dsq, char *aseq, int alen);

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
extern int          ParsetreeCompare(Parsetree_t *t1, Parsetree_t *t2);
extern void         SummarizeMasterTrace(FILE *fp, Parsetree_t *tr);
extern void         MasterTraceDisplay(FILE *fp, Parsetree_t *mtr, CM_t *cm);


/* from rna_ops.c
 */
extern int KHS2ct(char *ss, int len, int allow_pseudoknots, int **ret_ct);
extern int IsCompensatory(float *pij, int symi, int symj);

/* from scancyk.c
 */
extern void  CYKScan(CM_t *cm, char *dsq, int L, int W, 
		     int *ret_nhits, int **ret_hitr, 
		     int **ret_hiti, int **ret_hitj, float **ret_hitsc);
extern float CYKScanRequires(CM_t *cm, int L, int W);

/* from smallcyk.c
 */
extern float CYKDivideAndConquer(CM_t *cm, char *dsq, int L,
				 int r, int i0, int j0, Parsetree_t **ret_tr);
extern float CYKInside(CM_t *cm, char *dsq, int L,
		       int r, int i0, int j0, Parsetree_t **ret_tr);
extern float CYKInsideScore(CM_t *cm, char *dsq, int L,
			    int r, int i0, int j0);
extern void  CYKDemands(CM_t *cm, int L);


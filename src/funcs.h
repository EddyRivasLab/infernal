#include "squid.h"
#include "msa.h"
#include "structs.h"
#include "hmmer_structs.h"
#include "cplan9.h"

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
extern double **BandDistribution(CM_t *cm, int W, int do_local);
extern int      BandCalculationEngine(CM_t *cm, int W, double p_thresh, 
				      int save_densities,
				      int **ret_dmin, int **ret_dmax, 
				      double ***ret_gamma, int do_local);
extern int      BandTruncationNegligible(double *density, int b, int W, double *ret_beta);
extern int      BandMonteCarlo(CM_t *cm, int nsample, int W, double ***ret_gamma);
extern void     FreeBandDensities(CM_t *cm, double **gamma);


extern void     BandBounds(double **gamma, int M, int W, double p, 
			   int **ret_min, int **ret_max);
extern void     PrintBandGraph(FILE *fp, double **gamma, int *min, int *max, int v, int W);

extern void     PrintDPCellsSaved(CM_t *cm, int *min, int *max, int W);
extern void     CYKBandedScan(CM_t *cm, char *dsq, int *dmin, int *dmax, int L, int W, 
			      int *ret_nhits, int **ret_hitr, 
			      int **ret_hiti, int **ret_hitj, float **ret_hitsc,
			      float min_thresh);
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

/*EPN 10.19.05*/
extern void  CMDefaultNullModel(float *null);
extern void  CMSetNullModel(CM_t *cm, float null[MAXABET]);
extern void  CMReadNullModel(char *rndfile, float *null);

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
extern int            IsCompensatory(float *pij, int symi, int symj);

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


/* from scancyk.c
 */
extern void  CYKScan(CM_t *cm, char *dsq, int L, int W, 
		     int *ret_nhits, int **ret_hitr, 
		     int **ret_hiti, int **ret_hitj, float **ret_hitsc,
		     float min_thresh);
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
extern float CYKDivideAndConquer_b(CM_t *cm, char *dsq, int L, int r, int i0, int j0, Parsetree_t **ret_tr,
				   int *dmin, int *dmax);
extern float CYKInside_b(CM_t *cm, char *dsq, int L, int r, int i0, int j0, Parsetree_t **ret_tr, 
			 int *dmin, int *dmax);

/* The memory management routines.
 */
extern struct  deckpool_s *deckpool_create(void);
extern void    deckpool_push(struct deckpool_s *dpool, float **deck);
extern int     deckpool_pop(struct deckpool_s *d, float ***ret_deck);
extern void    deckpool_free(struct deckpool_s *d);
extern float **alloc_vjd_deck(int L, int i, int j);
extern float   size_vjd_deck(int L, int i, int j);
extern void    free_vjd_deck(float **a, int i, int j);
extern void    free_vjd_matrix(float ***a, int M, int i, int j);
extern char  **alloc_vjd_yshadow_deck(int L, int i, int j);
extern float   size_vjd_yshadow_deck(int L, int i, int j);
extern void    free_vjd_yshadow_deck(char **a, int i, int j);
extern int   **alloc_vjd_kshadow_deck(int L, int i, int j);
extern float   size_vjd_kshadow_deck(int L, int i, int j);
extern void    free_vjd_kshadow_deck(int **a, int i, int j);
extern void    free_vjd_shadow_matrix(void ***shadow, CM_t *cm, int i, int j);
extern float **alloc_vji_deck(int i0, int i1, int j1, int j0);
extern float   size_vji_deck(int i0, int i1, int j1, int j0);
extern void    free_vji_deck(float **a, int j1, int j0);
extern void    free_vji_matrix(float ***a, int M, int j1, int j0);
extern char  **alloc_vji_shadow_deck(int i0, int i1, int j1, int j0);
extern float   size_vji_shadow_deck(int i0, int i1, int j1, int j0);
extern void    free_vji_shadow_deck(char **a, int j1, int j0);
extern void    free_vji_shadow_matrix(char ***a, int M, int j1, int j0);

extern float **alloc_banded_vjd_deck(int L, int i, int j, int min, int max);
extern char  **alloc_banded_vjd_yshadow_deck(int L, int i, int j, int min, int max);
extern int   **alloc_banded_vjd_kshadow_deck(int L, int i, int j, int min, int max);

extern void debug_print_alpha(float ***alpha, CM_t *cm, int L);
extern void debug_print_alpha_banded(float ***alpha, CM_t *cm, int L, int *dmin, int *dmax);
extern void debug_print_alpha_deck(int v, float **deck, CM_t *cm, int L);
extern void debug_print_shadow(void ***shadow, CM_t *cm, int L);
extern void debug_print_shadow_banded(void ***shadow, CM_t *cm, int L, int *dmin, int *dmax);
extern void debug_print_shadow_banded_deck(int v, void ***shadow, CM_t *cm, int L, int *dmin, int *dmax);

/* from cplan9.c 
 * CM Plan9 HMM structure support
 */
extern struct cplan9_s *AllocCPlan9(int M);
extern struct cplan9_s *AllocCPlan9Shell(void);
extern void AllocCPlan9Body(struct cplan9_s *hmm, int M);
extern void FreeCPlan9(struct cplan9_s *hmm);
extern void ZeroCPlan9(struct cplan9_s *hmm);
extern void CPlan9SetName(struct cplan9_s *hmm, char *name);
extern void CPlan9SetAccession(struct cplan9_s *hmm, char *acc);
extern void CPlan9SetDescription(struct cplan9_s *hmm, char *desc);
extern void CPlan9ComlogAppend(struct cplan9_s *hmm, int argc, char **argv);
extern void CPlan9SetCtime(struct cplan9_s *hmm);
extern void CPlan9SetNullModel(struct cplan9_s *hmm, float null[HMMERMAXABET], float p1);
extern void CP9Logoddsify(struct cplan9_s *hmm);
extern void CPlan9Rescale(struct cplan9_s *hmm, float scale);
extern void CPlan9Renormalize(struct cplan9_s *hmm);

extern struct cp9_dpmatrix_s *AllocCPlan9Matrix(int rows, int M, int ***mmx, 
						int ***imx, int ***dmx, int ***emx);
extern float SizeCPlan9Matrix(int rows, int M);
extern void FreeCPlan9Matrix(struct cp9_dpmatrix_s *mx);
extern struct cp9_dpmatrix_s *CreateCPlan9Matrix(int N, int M, int padN, int padM);
extern void ResizeCPlan9Matrix(struct cp9_dpmatrix_s *mx, int N, int M, 
			       int ***mmx, int ***imx, int ***dmx, int ***emx);

/* from cp9_hmmio.c 
 * CM Plan9 HMM Input/output (saving/reading)
 */
extern CP9HMMFILE *CP9_HMMFileOpen(char *hmmfile, char *env);
extern int      CP9_HMMFileRead(CP9HMMFILE *hmmfp, struct cplan9_s **ret_hmm);
extern void     CP9_HMMFileClose(CP9HMMFILE *hmmfp);
extern int      CP9_HMMFileFormat(CP9HMMFILE *hmmfp);
extern void     CP9_HMMFileRewind(CP9HMMFILE *hmmfp);
extern int      CP9_HMMFilePositionByName(CP9HMMFILE *hmmfp, char *name);
extern int      CP9_HMMFilePositionByIndex(CP9HMMFILE *hmmfp, int idx);
extern void     CP9_WriteAscHMM(FILE *fp, struct cplan9_s *hmm);
extern void     CP9_WriteBinHMM(FILE *fp, struct cplan9_s *hmm);

/* from hbandcyk.c
 */
extern float CYKInside_b_jd(CM_t *cm, char *dsq, int L, int r, int i0, int j0, 
			    Parsetree_t **ret_tr, int *jmin, int *jmax, 
			    int **hdmin, int **hdmax, int *dmin, int *dmax);
extern void PrintDPCellsSaved_jd(CM_t *cm, int *jmin, int *jmax, int **hdmin, int **hdmax,
		     int W);
extern void ij2d_bands(CM_t *cm, int L, int *imin, int *imax, int *jmin, int *jmax,
		       int **hdmin, int **hdmax, int debug_level);
extern void hd2safe_hd_bands(int M, int *jmin, int *jmax, int **hdmin, int **hdmax,
			     int *safe_hdmin, int *safe_hdmax);
extern void debug_print_hd_bands(CM_t *cm, int **hdmin, int **hdmax, int *jmin, int *jmax);

/* from CP9_scan.c */
extern float CP9ForwardScan(unsigned char *dsq, int L, int W, struct cplan9_s *hmm, 
			    struct cp9_dpmatrix_s **ret_mx, int *ret_nhits, int **ret_hitr,
			    int **ret_hiti, int **ret_hitj,  
			    float **ret_hitsc, float min_thresh);

extern float CP9ForwardScanRequires(struct cplan9_s *hmm, int L, int W);

/* from CP9_cm2wrhmm.c */
extern int CP9_cm2wrhmm(CM_t *cm, struct cplan9_s *hmm, int *node_cc_left, int *node_cc_right, 
			int *cc_node_map, int **cs2hn_map, int **cs2hs_map, int ***hns2cs_map, 
			int debug_level);
extern int CP9_check_wrhmm(CM_t *cm, struct cplan9_s *hmm, int ***hns2cs_map, int *cc_node_map,
			   int debug_level);

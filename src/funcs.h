#include "squid.h"
#include "easel.h"
#include "msa.h"
#include "structs.h"
#include "cplan9.h"
#include "prior.h"

#ifdef USE_MPI
#include "mpi.h"
#endif


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
extern float    CYKBandedScan(CM_t *cm, char *dsq, int *dmin, int *dmax, int i0, int j0, int W, 
			      float cutoff, float score_boost, scan_results_t *results);
extern void     BandedParsetreeDump(FILE *fp, Parsetree_t *tr, CM_t *cm, char *dsq, 
				    double **gamma, int W, int *dmin, int *dmax);
extern void     ExpandBands(CM_t *cm, int qlen, int *dmin, int *dmax);
extern void     qdb_trace_info_dump(CM_t *cm, Parsetree_t *tr, int *dmin, 
				    int *dmax, int bdump_level);

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
extern void ConfigNoLocalEnds(CM_t *cm);
extern void ConfigLocalEnds(CM_t *cm, float p_internal_exit);
extern void ConfigLocal_fullsub(CM_t *cm, float p_internal_start, 
				float p_internal_exit, int sstruct_nd,
				int estruct_nd);
extern void ConfigLocal_DisallowELEmissions(CM_t *cm);

/* from modelmaker.c
 */
extern void HandModelmaker(MSA *msa, char **dsq, int use_rf, float gapthresh, 
			   CM_t **ret_cm, Parsetree_t **ret_mtr);
extern void ConsensusModelmaker(char *ss_cons, int clen, CM_t **ret_cm, 
				Parsetree_t **ret_gtr);
extern Parsetree_t *Transmogrify(CM_t *cm, Parsetree_t *gtr, 
				 char *dsq, char *aseq, int alen);
extern void cm_from_guide(CM_t *cm, Parsetree_t *gtr);
extern int  cm_find_and_detach_dual_inserts(CM_t *cm, int do_check, int do_detach);
extern int  cm_check_before_detaching(CM_t *cm, int insert1, int insert2);
extern int  cm_detach_state(CM_t *cm, int insert1, int insert2);
				 
/* from parsetree.c
 */
extern Parsetree_t *CreateParsetree(void);
extern void         GrowParsetree(Parsetree_t *tr);
extern void         FreeParsetree(Parsetree_t *tr);
extern int          InsertTraceNode(Parsetree_t *tr, int y, int whichway, 
				    int emitl, int emitr, int state);
extern void         ParsetreeCount(CM_t *cm, Parsetree_t *tr, char *seq, float wgt);
extern float        ParsetreeScore(CM_t *cm, Parsetree_t *tr, char *dsq, int do_null2);
extern void         PrintParsetree(FILE *fp, Parsetree_t *tr);
extern void         ParsetreeDump(FILE *fp, Parsetree_t *tr, CM_t *cm, char *dsq);
extern int          ParsetreeCompare(Parsetree_t *t1, Parsetree_t *t2);
extern void         SummarizeMasterTrace(FILE *fp, Parsetree_t *tr);
extern void         MasterTraceDisplay(FILE *fp, Parsetree_t *mtr, CM_t *cm);
extern MSA *Parsetrees2Alignment(CM_t *cm, char **dsq, SQINFO *sqinfo, float *wgt, 
				 Parsetree_t **tr, int nseq, int do_full);

/* from scancyk.c
 */
extern scan_results_t *CreateResults (int size);
extern void ExpandResults (scan_results_t *r, int additional);
extern void FreeResults (scan_results_t *r);
extern int  compare_results (const void *a_void, const void *b_void);
extern void sort_results (scan_results_t *results);
extern void report_hit (int i, int j, int bestr, float score, scan_results_t *results);
extern void remove_overlapping_hits (scan_results_t *results, int L);
extern float  CYKScan(CM_t *cm, char *dsq, int i0, int j0, int W, 
		      float cutoff, float score_boost, scan_results_t *results);
extern float CYKScanRequires(CM_t *cm, int L, int W);

/* from smallcyk.c
 */
extern float CYKDivideAndConquer(CM_t *cm, char *dsq, int L, int r, int i0, int j0, 
				 Parsetree_t **ret_tr, int *dmin, int *dmax);
extern float CYKInside(CM_t *cm, char *dsq, int L, int r, int i0, int j0, 
		       Parsetree_t **ret_tr, int *dmin, int *dmax);
extern float CYKInsideScore(CM_t *cm, char *dsq, int L, int r, int i0, 
			    int j0, int *dmin, int *dmax);
extern void  CYKDemands(CM_t *cm, int L, int *dmin, int *dmax);
extern void  debug_print_bands(CM_t *cm, int *dmin, int *dmax);

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

/* from cm_eweight.c
 * Entropy-based sequence weighting */
extern double CM_Eweight(CM_t *cm,  Prior_t *pri, 
			 float numb_seqs, float targetent);
extern void ModelContent(float *ent1, float *ent2, int M);
extern void CMRescale(CM_t *hmm, float scale);
extern double CM_Eweight_RE(CM_t *cm, Prior_t *pri, float numb_seqs, 
			    float target_relent, float *randomseq);
extern double DRelEntropy(double *p, double *f, int n);


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
extern void CPlan9SetNullModel(struct cplan9_s *hmm, float null[MAXABET], float p1);
extern void CP9Logoddsify(struct cplan9_s *hmm);
extern void CPlan9Rescale(struct cplan9_s *hmm, float scale);
extern void CPlan9Renormalize(struct cplan9_s *hmm);
extern struct cp9_dpmatrix_s *AllocCPlan9Matrix(int rows, int M, int ***mmx, 
						int ***imx, int ***dmx, int ***emx);
extern float SizeCPlan9Matrix(int rows, int M);
extern void  FreeCPlan9Matrix(struct cp9_dpmatrix_s *mx);
extern struct cp9_dpmatrix_s *CreateCPlan9Matrix(int N, int M, int padN, int padM);
extern void  ResizeCPlan9Matrix(struct cp9_dpmatrix_s *mx, int N, int M, 
			       int ***mmx, int ***imx, int ***dmx, int ***emx);
extern void CPlan9SWConfig(struct cplan9_s *hmm, float pentry, float pexit);
extern void CPlan9GlobalConfig(struct cplan9_s *hmm);
extern void sub_CPlan9GlobalConfig(struct cplan9_s *hmm, int spos, int epos, double **phi);

extern void CPlan9RenormalizeExits(struct cplan9_s *hmm, int spos);
extern void CP9AllocTrace(int tlen, struct cp9trace_s **ret_tr);
extern void CP9ReallocTrace(struct cp9trace_s *tr, int tlen);
extern void CP9FreeTrace(struct cp9trace_s *tr);

extern void CP9_2sub_cp9(struct cplan9_s *orig_hmm, struct cplan9_s **ret_sub_hmm, int spos, int epos, double **orig_phi);
extern void CP9_reconfig2sub(struct cplan9_s *hmm, int spos, int epos, int spos_nd, int epos_nd, double **orig_phi);

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
extern void debug_print_alpha_banded_jd(float ***alpha, CM_t *cm, int L, int *jmin, int *jmax, 
					int **hdmin, int **hdmax);
extern float ** alloc_jdbanded_vjd_deck(int L, int i, int j, int jmin, int jmax, int *hdmin, int *hdmax);
extern void CYKBandedScan_jd(CM_t *cm, char *dsq, int *jmin, int *jmax, int **hdmin, int **hdmax, int i0, 
			     int j0, int W, int *ret_nhits, int **ret_hitr, int **ret_hiti, 
			     int **ret_hitj, float **ret_hitsc, float min_thresh);


/* from CP9_scan.c */
extern float CP9ForwardScan(char *dsq, int i0, int j0, int W, struct cplan9_s *hmm, 
			    struct cp9_dpmatrix_s **ret_mx, int *ret_nhits, int **ret_hitr,
			    int **ret_hiti, int **ret_hitj,  
			    float **ret_hitsc, float min_thresh);
extern float CP9BackwardScan(char *dsq, int i0, int j0, int W, struct cplan9_s *hmm, 
			     struct cp9_dpmatrix_s **ret_mx, int *ret_nhits, int **ret_hitr,
			     int **ret_hiti, int **ret_hitj,  
			     float **ret_hitsc, float min_thresh);
extern float CP9ForwardScanRequires(struct cplan9_s *hmm, int L, int W);
extern float CP9ForwardBackwardScan(char *dsq, int i0, int j0, int W, 
				    struct cplan9_s *hmm, struct cp9_dpmatrix_s **ret_fmx,
				    struct cp9_dpmatrix_s **ret_bmx, int *ret_nhits, 
				    int **ret_hitr, int **ret_hiti, 
				    int **ret_hitj, float **ret_hitsc, float min_thresh, int pad);
extern void CP9ScanFullPosterior(char *dsq, int L,
				 struct cplan9_s *hmm,
				 struct cp9_dpmatrix_s *fmx,
				 struct cp9_dpmatrix_s *bmx,
				 struct cp9_dpmatrix_s *mx);
extern void CP9_combine_FBscan_hits(int i0, int j0, int W, int fwd_nhits, int *fwd_hitr, int *fwd_hiti, 
				    int *fwd_hitj, float *fwd_hitsc, int bck_nhits, 
				    int *bck_hitr, int *bck_hiti, int *bck_hitj, float *bck_hitsc, 
				    int *ret_nhits, int **ret_hitr, int **ret_hiti, int **ret_hitj, 
				    float **ret_hitsc, int pad);

/* from CP9_cm2wrhmm.c */
extern int  build_cp9_hmm(CM_t *cm, struct cplan9_s **ret_hmm, CP9Map_t **ret_cp9map, 
			  int do_psi_test, float psi_vs_phi_threshold, int debug_level);
extern void CP9_map_cm2hmm(CM_t *cm, CP9Map_t *cp9map, int debug_level);
extern void map_helper(CM_t *cm, CP9Map_t *cp9map, int k, int ks, int v);
extern int  CP9_check_wrhmm(CM_t *cm, struct cplan9_s *hmm, CP9Map_t *cp9map, int debug_level);
extern void fill_psi(CM_t *cm, double *psi, char ***tmap);
extern void fill_phi_cp9(struct cplan9_s *hmm, double ***ret_phi, int spos);
extern void make_tmap(char ****ret_tmap);

extern int  CP9_check_by_sampling(CM_t *cm, struct cplan9_s *hmm, CMSubInfo_t *subinfo, int spos, int epos, 
				  float chi_thresh, int nsamples, int print_flag);
extern void CP9_fake_tracebacks(char **aseq, int nseq, int alen, int *matassign, struct cp9trace_s ***ret_tr);

extern void CP9TraceCount(struct cplan9_s *hmm, char *dsq, float wt, struct cp9trace_s *tr);
extern void debug_print_cp9_params(struct cplan9_s *hmm);
extern void debug_print_phi_cp9(struct cplan9_s *hmm, double **phi);
extern CP9Map_t *AllocCP9Map(CM_t *cm);
extern void FreeCP9Map(CP9Map_t *cp9map);



/* from scaninside.c */
extern void  InsideScan(CM_t *cm, char *dsq, int i0, int j0, int W, 
			int *ret_nhits, int **ret_hitr, 
			int **ret_hiti, int **ret_hitj, float **ret_hitsc,
			float min_thresh);
extern void  InsideBandedScan(CM_t *cm, char *dsq, int *dmin, int *dmax, int i0, int j0, int W, 
			      int *ret_nhits, int **ret_hitr, 
			      int **ret_hiti, int **ret_hitj, float **ret_hitsc,
			      float min_thresh);
extern void  InsideBandedScan_jd(CM_t *cm, char *dsq, int *jmin, int *jmax, int **hdmin, int **hdmax,
				 int i0, int j0, int W, 
				 int *ret_nhits, int **ret_hitr, 
				 int **ret_hiti, int **ret_hitj, float **ret_hitsc,
				 float min_thresh);
extern float LogSum2(float p1, float p2);

/* from cm_masks.c */
extern float CM_TraceScoreCorrection(CM_t *cm, Parsetree_t *tr, char *dsq);

/* from sub_cm.c */
extern int  build_sub_cm(CM_t *orig_cm, CM_t **ret_cm, int sstruct, int estruct, CMSubMap_t **ret_submap, 
			 int do_fullsub, int print_flag);
extern void CP9NodeForPosn(struct cplan9_s *hmm, int i0, int j0, int x, 
			   struct cp9_dpmatrix_s *post, int *ret_node, int *ret_type, int print_flag);
extern void StripWUSSGivenCC(MSA *msa, char **dsq, float gapthresh, int first_match, int last_match);
extern int  check_orig_psi_vs_sub_psi(CM_t *orig_cm, CM_t *sub_cm, CMSubMap_t *submap, double threshold, 
				       int print_flag);
extern int  check_sub_cm(CM_t *orig_cm, CM_t *sub_cm, CMSubMap_t *submap, CMSubInfo_t *subinfo, float pthresh,
			 int do_fullsub, int print_flag);
extern int check_sub_cm_by_sampling(CM_t *orig_cm, CM_t *sub_cm, CMSubMap_t *submap, CMSubInfo_t *subinfo,
				    float chi_thresh, int nsamples, int do_fullsub, int print_flag);
extern int  check_sub_cm_by_sampling2(CM_t *orig_cm, CM_t *sub_cm, int spos, int epos, int nseq);
extern int  sub_cm2cm_parsetree(CM_t *orig_cm, CM_t *sub_cm, Parsetree_t **ret_orig_tr, Parsetree_t *sub_tr, 
				CMSubMap_t *submap, int do_fullsub, int print_flag);
extern CMSubMap_t  *AllocSubMap(CM_t *sub_cm, CM_t *orig_cm, int sstruct, int estruct, int do_fullsub);
extern void         FreeSubMap(CMSubMap_t *submap);
extern CMSubInfo_t *AllocSubInfo(int clen);
extern void         FreeSubInfo(CMSubInfo_t *subinfo);
extern void  debug_print_cm_params(CM_t *cm);

/* from cm_wrappers.c
 */
extern void
AlignSeqsWrapper(CM_t *cm, char **dsq, SQINFO *sqinfo, int nseq, Parsetree_t ***ret_tr, 
		 int do_local, int do_small, int do_qdb, double qdb_beta, int do_hbanded, 
		 int use_sums, double hbandp, int do_sub, int do_fullsub, int do_hmmonly, 
		 int do_inside, int do_outside, int do_check, int do_post, char ***ret_postcode, 
		 int do_timings, int bdump_level, int debug_level, int silent_mode);
extern void 
serial_search_database (ESL_SQFILE *dbfp, CM_t *cm, CMConsensus_t *cons, int W, int cutoff_type, 
			float cutoff, int do_revcomp, int do_align, int do_stats, double *mu, 
			double *lambda, int *dmin, int *dmax);
extern void parallel_search_database (ESL_SQFILE *dbfp, CM_t *cm, CMConsensus_t *cons,
				      int W, int cutoff_type, float cutoff, 
				      int do_revcomp, int do_align, int do_stats,
				      double *mu, double *lambda, int *dmin, int *dmax,
				      int mpi_my_rank, int mpi_master_rank, 
				      int mpi_num_procs) ;




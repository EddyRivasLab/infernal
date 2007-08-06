#include "esl_config.h"
#include "config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_msa.h"

#include "structs.h"
#include "cplan9.h"
#include "prior.h"
#include "rnamat.h"

#ifdef USE_MPI
#include "mpi.h"
#endif


/* from alphabet.c
 */
extern void   PairCount(const ESL_ALPHABET *abc, float *counters, char syml, char symr, float wt);
extern float  DegeneratePairScore(const ESL_ALPHABET *abc, float *esc, char syml, char symr);
extern int    iDegeneratePairScore(const ESL_ALPHABET *abc, int *esc, char syml, char symr);
float         LeftMarginalScore(const ESL_ALPHABET *abc, float *esc, int dres);
float         RightMarginalScore(const ESL_ALPHABET *abc, float *esc, int dres);

/* from bandcyk.c
 */
extern void     BandExperiment(CM_t *cm);
extern double **BandDistribution(CM_t *cm, int W, int do_local);
extern int      BandCalculationEngine(CM_t *cm, int W, double p_thresh, 
				      int save_densities,
				      int **ret_dmin, int **ret_dmax, 
				      double ***ret_gamma);
extern int      BandTruncationNegligible(double *density, int b, int W, double *ret_beta);
extern int      BandMonteCarlo(CM_t *cm, int nsample, int W, double ***ret_gamma);
extern void     FreeBandDensities(CM_t *cm, double **gamma);


extern void     BandBounds(double **gamma, int M, int W, double p, 
			   int **ret_min, int **ret_max);
extern void     PrintBandGraph(FILE *fp, double **gamma, int *min, int *max, int v, int W);

extern void     PrintDPCellsSaved(CM_t *cm, int *min, int *max, int W);
extern float    CYKBandedScan(CM_t *cm, ESL_SQ *sq, int *dmin, int *dmax, int i0, int j0, int W, 
			      float cutoff, scan_results_t *results);
extern void     ExpandBands(CM_t *cm, int qlen, int *dmin, int *dmax);
extern void     qdb_trace_info_dump(CM_t *cm, Parsetree_t *tr, int *dmin, 
				    int *dmax, int bdump_level);

/* from cm.c
 */
extern CM_t *CreateCM(int nnodes, int nstates, const ESL_ALPHABET *abc);
extern CM_t *CreateCMShell(void);
extern void  CreateCMBody(CM_t *cm, int nnodes, int nstates, const ESL_ALPHABET *abc);
extern void  CMZero(CM_t *cm);
extern void  CMRenormalize(CM_t *cm);
extern void  FreeCM(CM_t *cm);
extern void  CMSimpleProbify(CM_t *cm);
extern int   rsearch_CMProbifyEmissions(CM_t *cm, fullmat_t *fullmat);
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
extern int **IMX2Alloc(int rows, int cols);
extern void  IMX2Free(int **mx);
extern float rsearch_calculate_gap_penalty (char from_state, char to_state, 
					    int from_node, int to_node, 
					    float input_alpha, float input_beta, 
					    float input_alphap, float input_betap);
extern int   ExponentiateCM(CM_t *cm, double z);
extern CM_t *DuplicateCM(CM_t *cm);
extern void  cm_banner(FILE *fp, char *progname, char *banner);
extern int   cm_Validate(CM_t *cm, float tol, char *errbuf);
extern char *CMStatetype(char st);
extern char *CMNodetype(char nd);
extern char *CMStateid(char st);
extern int   cm_SetName(CM_t *cm, char *name);
extern int   cm_SetAccession(CM_t *cm, char *acc);
extern int   cm_SetDescription(CM_t *cm, char *desc);
extern int   cm_AppendComlog(CM_t *cm, int argc, char **argv);
extern int   cm_SetCtime(CM_t *cm);

extern void  CMCreateNullModel(CM_t *cm);
extern void  CMSetNullModel(CM_t *cm, float *null);
extern void  CMReadNullModel(CM_t *cm, char *rndfile);


/* from cmio.c
 */
extern CMFILE *CMFileOpen(char *cmfile, char *env);
extern int     CMFileRead(CMFILE *cmf, ESL_ALPHABET **ret_abc, CM_t **ret_cm);
extern void    CMFileClose(CMFILE *cmf);
extern void    CMFileRewind(CMFILE *cmf);
extern int     CMFilePositionByIndex(CMFILE *cmf, int idx);
extern int     CMFilePositionByKey(CMFILE *cmf, char *key);
extern int     CMFileWrite(FILE *fp, CM_t *cm, int do_binary);


/* from display.c
 */
extern Fancyali_t    *CreateFancyAli(Parsetree_t *tr, CM_t *cm, CMConsensus_t *cons, ESL_SQ *sq);
extern void           PrintFancyAli(FILE *fp, Fancyali_t *ali, int offset, int in_revcomp);
extern void           FreeFancyAli(Fancyali_t *ali);
extern CMConsensus_t *CreateCMConsensus(CM_t *cm, float pthresh, float sthresh);
extern void           FreeCMConsensus(CMConsensus_t *con);
extern void           MainBanner(FILE *fp, char *banner); 
extern int            IsCompensatory(const ESL_ALPHABET *abc, float *pij, int symi, int symj);

/* in emit.c
 */
extern int EmitParsetree(CM_t *cm, ESL_RANDOMNESS *r, char *name, int do_digital, Parsetree_t **ret_tr, ESL_SQ **ret_sq, int *ret_N);

/* in emitmap.c
 */
extern CMEmitMap_t *CreateEmitMap(CM_t *cm); 
extern void         DumpEmitMap(FILE *fp, CMEmitMap_t *map, CM_t *cm);
extern void         FreeEmitMap(CMEmitMap_t *map);

/* from modelconfig.c
 */
extern void  ConfigCM(CM_t *cm, int *preset_dmin, int *preset_dmax);
extern void  ConfigCMEnforce(CM_t *cm);
extern void  ConfigLocal(CM_t *cm, float p_internal_start, float p_internal_exit);
extern void  ConfigGlobal(CM_t *cm);
extern void  ConfigNoLocalEnds(CM_t *cm);
extern void  ConfigLocalEnds(CM_t *cm, float p_internal_exit);
extern void  ConfigLocal_fullsub(CM_t *cm, float p_internal_start, 
				float p_internal_exit, int sstruct_nd,
				int estruct_nd);
extern void  ConfigLocal_DisallowELEmissions(CM_t *cm);
extern void  ConfigLocal_fullsub_post(CM_t *sub_cm, CM_t *orig_cm, CP9Map_t *orig_cp9map, CMSubMap_t *submap,
				     struct cp9_dpmatrix_s *post, int L);
extern void  ConfigLocalEnforce(CM_t *cm, float p_internal_start, float p_internal_exit);
extern int   EnforceSubsequence(CM_t *cm);
extern float EnforceScore(CM_t *cm);
extern int   EnforceFindEnfStart(CM_t *cm, int enf_cc_start);
extern int   ConfigForGumbelMode(CM_t *cm, int statmode);
extern int   ConfigQDB(CM_t *cm);

/* from modelmaker.c
 */
extern void HandModelmaker(ESL_MSA *msa, int use_rf, float gapthresh, 
			   CM_t **ret_cm, Parsetree_t **ret_mtr);
extern void ConsensusModelmaker(ESL_ALPHABET *abc, char *ss_cons, int clen, CM_t **ret_cm, 
				Parsetree_t **ret_gtr);
extern Parsetree_t *Transmogrify(CM_t *cm, Parsetree_t *gtr, 
				 ESL_DSQ *dsq, char *aseq, int alen);
extern void cm_from_guide(CM_t *cm, Parsetree_t *gtr);
extern int  cm_find_and_detach_dual_inserts(CM_t *cm, int do_check, int do_detach);
extern int  cm_check_before_detaching(CM_t *cm, int insert1, int insert2);
extern int  cm_detach_state(CM_t *cm, int insert1, int insert2);
extern int  clean_cs(char *cs, int alen);
				 
/* from parsetree.c
 */
extern Parsetree_t *CreateParsetree(void);
extern void         GrowParsetree(Parsetree_t *tr);
extern void         FreeParsetree(Parsetree_t *tr);
extern int          InsertTraceNode(Parsetree_t *tr, int y, int whichway, 
				    int emitl, int emitr, int state);
extern int          InsertTraceNodewithMode(Parsetree_t *tr, int y, int whichway, 
				    int emitl, int emitr, int state, int mode);
extern void         ParsetreeCount(CM_t *cm, Parsetree_t *tr, ESL_SQ *sq, float wgt);
extern float        ParsetreeScore(CM_t *cm, Parsetree_t *tr, ESL_SQ *sq, int do_null2);
extern void         PrintParsetree(FILE *fp, Parsetree_t *tr);
extern void         ParsetreeDump(FILE *fp, Parsetree_t *tr, CM_t *cm, ESL_SQ *sq, int *dmin, int *dmax);
extern int          ParsetreeCompare(Parsetree_t *t1, Parsetree_t *t2);
extern void         SummarizeMasterTrace(FILE *fp, Parsetree_t *tr);
extern void         MasterTraceDisplay(FILE *fp, Parsetree_t *mtr, CM_t *cm);
extern int          Parsetrees2Alignment(CM_t *cm, ESL_SQ **sq, float *wgt, 
					 Parsetree_t **tr, int nseq, int do_full, ESL_MSA **ret_msa);
extern float        ParsetreeScore_Global2Local(CM_t *cm, Parsetree_t *tr, ESL_SQ *sq, int print_flag);
extern int          Parsetree2CP9trace(CM_t *cm, Parsetree_t *tr, CP9trace_t **ret_cp9_tr);

/* from scancyk.c
 */
extern scan_results_t *CreateResults (int size);
extern void ExpandResults (scan_results_t *r, int additional);
extern void FreeResults (scan_results_t *r);
extern int  compare_results (const void *a_void, const void *b_void);
extern void sort_results (scan_results_t *results);
extern void report_hit (int i, int j, int bestr, float score, scan_results_t *results);
extern void remove_overlapping_hits (scan_results_t *results, int i0, int j0);
extern float CYKScan(CM_t *cm, ESL_SQ *sq, int i0, int j0, int W, 
		      float cutoff, scan_results_t *results);
extern float CYKScanRequires(CM_t *cm, int L, int W);

/* from smallcyk.c
 */
extern float CYKDivideAndConquer(CM_t *cm, ESL_SQ *sq, int r, int i0, int j0, 
				 Parsetree_t **ret_tr, int *dmin, int *dmax);
extern float CYKInside(CM_t *cm, ESL_SQ *sq, int r, int i0, int j0, 
		       Parsetree_t **ret_tr, int *dmin, int *dmax);
extern float CYKInsideScore(CM_t *cm, ESL_SQ *sq, int r, int i0, 
			    int j0, int *dmin, int *dmax);
extern float CYKDemands(CM_t *cm, int L, int *dmin, int *dmax, int be_quiet);
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

/* from truncyk.c
 */
float TrCYKInside(CM_t *cm, ESL_SQ *sq, int r, int i0, int j0, Parsetree_t **ret_tr, int *dmin, int *dmax);

/* from cm_eweight.c
 * Entropy-based sequence weighting */
extern double CM_Eweight(CM_t *cm,  const Prior_t *pri, 
			 float numb_seqs, float targetent);
extern void ModelContent(float *ent1, float *ent2, int M);
extern void CMRescale(CM_t *hmm, float scale);
extern double CM_Eweight_RE(CM_t *cm, Prior_t *pri, float numb_seqs, 
			    float target_relent, float *randomseq);
extern double DRelEntropy(double *p, double *f, int n);
extern float CMAverageMatchEntropy(CM_t *cm);
extern float CP9AverageMatchEntropy(CP9_t *cp9);

/* from cplan9.c 
 * CM Plan9 HMM structure support
 */
extern CP9_t *AllocCPlan9(int M, const ESL_ALPHABET *abc);
extern CP9_t *AllocCPlan9Shell(void);
extern void AllocCPlan9Body(CP9_t *hmm, int M, const ESL_ALPHABET *abc);
extern void FreeCPlan9(CP9_t *hmm);
extern void ZeroCPlan9(CP9_t *hmm);
extern void CPlan9SetName(CP9_t *hmm, char *name);
extern void CPlan9SetAccession(CP9_t *hmm, char *acc);
extern void CPlan9SetDescription(CP9_t *hmm, char *desc);
extern void CPlan9ComlogAppend(CP9_t *hmm, int argc, char **argv);
extern void CPlan9SetCtime(CP9_t *hmm);
extern void CPlan9SetNullModel(CP9_t *hmm, float null[MAXABET], float p1);
extern void CP9Logoddsify(CP9_t *hmm);
extern void CPlan9Rescale(CP9_t *hmm, float scale);
extern void CPlan9Renormalize(CP9_t *hmm);
extern struct cp9_dpmatrix_s *AllocCPlan9Matrix(int rows, int M, int ***mmx, 
						int ***imx, int ***dmx, int ***elmx, int **erow);
extern float SizeCPlan9Matrix(int rows, int M);
extern void  FreeCPlan9Matrix(struct cp9_dpmatrix_s *mx);
extern struct cp9_dpmatrix_s *CreateCPlan9Matrix(int N, int M, int padN, int padM);
extern void  ResizeCPlan9Matrix(struct cp9_dpmatrix_s *mx, int N, int M, 
			       int ***mmx, int ***imx, int ***dmx, int ***elmx, int **erow);
extern void CPlan9SWConfig(CP9_t *hmm, float pentry, float pexit);
extern void CPlan9SWConfigEnforce(CP9_t *hmm, float pentry, float pexit, 
				  int enf_start_pos, int enf_end_pos);
extern void CPlan9ELConfig(CM_t *cm);
extern void CPlan9NoEL(CM_t *cm);
extern void CPlan9InitEL(CM_t *cm, CP9_t *cp9);
extern void CPlan9GlobalConfig(CP9_t *hmm);
extern void sub_CPlan9GlobalConfig(CP9_t *hmm, int spos, int epos, double **phi);

extern void CPlan9RenormalizeExits(CP9_t *hmm, int spos);
extern void CP9AllocTrace(int tlen, CP9trace_t **ret_tr);
extern void CP9ReallocTrace(CP9trace_t *tr, int tlen);
extern void CP9FreeTrace(CP9trace_t *tr);

extern void CP9_2sub_cp9(CP9_t *orig_hmm, CP9_t **ret_sub_hmm, int spos, int epos, double **orig_phi);
extern void CP9_reconfig2sub(CP9_t *hmm, int spos, int epos, int spos_nd, int epos_nd, double **orig_phi);
extern void CP9HackInsertScores(CP9_t *cp9);
extern void CP9EnforceHackMatchScores(CP9_t *cp9, int enf_start_pos, int enf_end_pos);
extern void CP9_fake_tracebacks(ESL_MSA *msa, int *matassign, CP9trace_t ***ret_tr);
extern void  CP9TraceCount(CP9_t *hmm, ESL_SQ *sq, float wt, CP9trace_t *tr);
extern float CP9TraceScore(CP9_t *hmm, ESL_SQ *sq, CP9trace_t *tr);
extern void  CP9PrintTrace(FILE *fp, CP9trace_t *tr, CP9_t *hmm, ESL_SQ *sq);
extern char *CP9Statetype(char st);
extern int   CP9TransitionScoreLookup(struct cplan9_s *hmm, char st1, int k1, 
				    char st2, int k2);
extern void  CP9ViterbiTrace(struct cplan9_s *hmm, ESL_SQ *sq, int i0, int j0,
			     struct cp9_dpmatrix_s *mx, CP9trace_t **ret_tr);
extern void  CP9ReverseTrace(CP9trace_t *tr);
extern int   CP9Traces2Alignment(CM_t *cm, ESL_SQ **sq, float *wgt, int nseq, CP9trace_t **tr, 
				     int do_full, ESL_MSA **ret_msa);
extern void  DuplicateCP9(CM_t *src_cm, CM_t *dest_cm);

/* from hbandcyk.c
 */
extern float CYKInside_b_jd(CM_t *cm, ESL_SQ *sq, int r, int i0, int j0, 
			    Parsetree_t **ret_tr, int *jmin, int *jmax, 
			    int **hdmin, int **hdmax, int *dmin, int *dmax);
extern void PrintDPCellsSaved_jd(CM_t *cm, int *jmin, int *jmax, int **hdmin, int **hdmax,
		     int W);
extern void ij2d_bands(CM_t *cm, int L, int *imin, int *imax, int *jmin, int *jmax,
		       int **hdmin, int **hdmax, int debug_level);
extern void combine_qdb_hmm_d_bands(CM_t *cm, int *jmin, int *jmax, int **hdmin, int **hdmax);
extern void hd2safe_hd_bands(int M, int *jmin, int *jmax, int **hdmin, int **hdmax,
			     int *safe_hdmin, int *safe_hdmax);
extern void debug_print_hd_bands(CM_t *cm, int **hdmin, int **hdmax, int *jmin, int *jmax);
extern void debug_print_alpha_banded_jd(float ***alpha, CM_t *cm, int L, int *jmin, int *jmax, 
					int **hdmin, int **hdmax);
extern float ** alloc_jdbanded_vjd_deck(int L, int i, int j, int jmin, int jmax, int *hdmin, int *hdmax);
extern float CYKBandedScan_jd(CM_t *cm, ESL_SQ *sq, int *jmin, int *jmax, int **hdmin, int **hdmax, int i0, 
			      int j0, int W, float cutoff, scan_results_t *results);
extern float iInsideBandedScan_jd(CM_t *cm, ESL_SQ *sq, int *jmin, int *jmax, int **hdmin, int **hdmax, int i0, 
				  int j0, int W, float cutoff, scan_results_t *results);



/* from CP9_scan.c */
extern float CP9Forward(CM_t *cm, ESL_SQ *sq, int i0, int j0, int W, float cutoff, int **ret_isc, 
			int *ret_maxres, scan_results_t *results, int do_scan, int doing_align, 
			int doing_rescan, int be_efficient, CP9_dpmatrix_t **ret_mx);
extern float CP9Backward(CM_t *cm, ESL_SQ *sq, int i0, int j0, int W, float cutoff, int **ret_isc, 
			 int *ret_maxres, scan_results_t *results, int do_scan, int doing_align, 
			 int doing_rescan, int be_efficient, CP9_dpmatrix_t **ret_mx);
extern float CP9Scan_dispatch(CM_t *cm, ESL_SQ *sq, int i0, int j0, int W, float cm_cutoff, 
			      float cp9_cutoff, scan_results_t *results, int doing_cp9_stats, int *ret_flen);
extern float RescanFilterSurvivors(CM_t *cm, ESL_SQ *sq, scan_results_t *hmm_results, int i0, 
				   int j0, int W, int padmode, int ipad, int jpad, int do_collapse,
				   float cm_cutoff, float cp9_cutoff, scan_results_t *results, 
				   int *ret_flen);
extern void CP9ScanPosterior(ESL_SQ *sq, int i0, int j0, CP9_t *hmm, CP9_dpmatrix_t *fmx, 
			     CP9_dpmatrix_t *bmx, CP9_dpmatrix_t *mx);
extern float CP9ForwardScanDemands(CP9_t *cp9, int L);
extern float FindCP9FilterThreshold(CM_t *cm, CMStats_t *cmstats, ESL_RANDOMNESS *r, 
				    float Fmin, float Smin, float Starget, float Spad, int N, 
				    int use_cm_cutoff, float cm_ecutoff, int db_size, 
				    int emit_mode, int fthr_mode, int hmm_gum_mode, 
				    int do_fastfil, int do_Fstep, int my_rank, int nproc, 
				    int do_mpi, char *histfile, FILE *Rpts_fp, float *ret_F);
extern float FindExpFactor(CM_t *cm, CMStats_t *cmstats, ESL_RANDOMNESS *r, 
			   int use_cm_cutoff, float cm_ecutoff, int db_size, 
			   int emit_mode, int fthr_mode, int do_fastfil,
			   int ntrials, float fp_min, float fp_max);
extern float FindExpFactor_minmax(CM_t *cm, CMStats_t *cmstats, ESL_RANDOMNESS *r, 
				  float cm_min_ecutoff, float cm_max_ecutoff, int db_size, 
				  int emit_mode, int fthr_mode, int do_fastfil,
				  int ntrials, float fp_min, float fp_max);
extern float FindExpFactor_min(CM_t *cm, CMStats_t *cmstats, ESL_RANDOMNESS *r, 
			       float cm_min_ecutoff, int db_size, 
			       int emit_mode, int fthr_mode, int do_fastfil,
			       int ntrials, float fp_min);
extern float Filter_XFTableLookup(float X, float F, int emit_mode, int fthr_mode);
     
/* from CP9_cm2wrhmm.c */
extern int build_cp9_hmm(CM_t *cm, CP9_t **ret_hmm, CP9Map_t **ret_cp9map, int do_psi_test,
			 float psi_vs_phi_threshold, int debug_level);
extern void CP9_map_cm2hmm(CM_t *cm, CP9Map_t *cp9map, int debug_level);
extern void map_helper(CM_t *cm, CP9Map_t *cp9map, int k, int ks, int v);
extern int  CP9_check_wrhmm(CM_t *cm, CP9_t *hmm, CP9Map_t *cp9map, int debug_level);
extern void fill_psi(CM_t *cm, double *psi, char ***tmap);
extern void fill_phi_cp9(CP9_t *hmm, double ***ret_phi, int spos);
extern void make_tmap(char ****ret_tmap);

extern int  CP9_check_by_sampling(CM_t *cm, CP9_t *hmm, CMSubInfo_t *subinfo, int spos, int epos, 
				  float chi_thresh, int nsamples, int print_flag);
extern void debug_print_cp9_params(FILE *fp, CP9_t *hmm);
extern void debug_print_phi_cp9(CP9_t *hmm, double **phi);
extern CP9Map_t *AllocCP9Map(CM_t *cm);
extern void FreeCP9Map(CP9Map_t *cp9map);



/* from scaninside.c */
extern float  InsideScan(CM_t *cm, ESL_SQ *sq, int i0, int j0, int W, 
			 float cutoff, scan_results_t *results);
extern float  InsideBandedScan(CM_t *cm, ESL_SQ *sq, int *dmin, int *dmax, int i0, int j0, int W, 
			       float cutoff, scan_results_t *results);
extern void  InsideBandedScan_jd(CM_t *cm, ESL_SQ *sq, int *jmin, int *jmax, int **hdmin, int **hdmax,
				 int i0, int j0, int W, 
				 int *ret_nhits, int **ret_hitr, 
				 int **ret_hiti, int **ret_hitj, float **ret_hitsc,
				 float min_thresh);
extern float iInsideScan(CM_t *cm, ESL_SQ *sq, int i0, int j0, int W, 
			 float cutoff, scan_results_t *results);
extern float iInsideBandedScan(CM_t *cm, ESL_SQ *sq, int *dmin, int *dmax, int i0, int j0, int W, 
			       float cutoff, scan_results_t *results);
extern float LogSum2(float p1, float p2);

/* from cm_masks.c */
extern float CM_TraceScoreCorrection(CM_t *cm, Parsetree_t *tr, ESL_SQ *sq);

/* from sub_cm.c */
extern int  build_sub_cm(CM_t *orig_cm, CM_t **ret_cm, int sstruct, int estruct, CMSubMap_t **ret_submap, 
			 int do_fullsub, int print_flag);
extern void CP9NodeForPosn(CP9_t *hmm, int i0, int j0, int x, 
			   struct cp9_dpmatrix_s *post, int *ret_node, int *ret_type,
			   int do_fullsub, float pmass, int is_start, int print_flag);
extern void StripWUSSGivenCC(ESL_MSA *msa, char **dsq, float gapthresh, int first_match, int last_match);
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
extern void  debug_print_cm_params(FILE *fp, CM_t *cm);

/* Reading/writing of CP9 HMMs no longer supported. */
#if 0
/* from cp9_hmmio.c 
 * CM Plan9 HMM Input/output (saving/reading)
 */
extern CP9HMMFILE *CP9_HMMFileOpen(char *hmmfile, char *env);
extern int      CP9_HMMFileRead(CP9HMMFILE *hmmfp, CP9_t **ret_hmm);
extern void     CP9_HMMFileClose(CP9HMMFILE *hmmfp);
extern int      CP9_HMMFileFormat(CP9HMMFILE *hmmfp);
extern void     CP9_HMMFileRewind(CP9HMMFILE *hmmfp);
extern int      CP9_HMMFilePositionByName(CP9HMMFILE *hmmfp, char *name);
extern int      CP9_HMMFilePositionByIndex(CP9HMMFILE *hmmfp, int idx);
extern void     CP9_WriteAscHMM(FILE *fp, CP9_t *hmm);
extern void     CP9_WriteBinHMM(FILE *fp, CP9_t *hmm);
#endif

/* from rsearch_buildcm.c
 */
extern CM_t *build_cm (ESL_MSA *msa, fullmat_t *fullmat, int *querylen,
		       float alpha, float beta, float alphap, float betap,
		       float beginsc, float endsc);
extern CM_t *read_cm (char *queryfile);

/* from cm_cluster.c
 */
int MSADivide(ESL_MSA *mmsa, int do_all, int target_nc, float mindiff, int do_corig, 
	      int *ret_num_msa, ESL_MSA ***ret_cmsa);

/* cm_errors.c */
extern void cm_Die (char *format, ...);
extern void cm_Fail(char *format, ...);

#include "esl_config.h"
#include "config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_sqio.h"
#include "esl_msa.h"

#include "hmmer.h"

#include "structs.h"

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#define USE_NEWLOGSUM 1
#define USE_OLDLOGSUM 0

/* from cm.c */
extern CM_t *CreateCM(int nnodes, int nstates, const ESL_ALPHABET *abc);
extern CM_t *CreateCMShell(void);
extern void  CreateCMBody(CM_t *cm, int nnodes, int nstates, const ESL_ALPHABET *abc);
extern void  CMZero(CM_t *cm);
extern void  CMRenormalize(CM_t *cm);
extern void  FreeCM(CM_t *cm);
extern void  CMSimpleProbify(CM_t *cm);
extern int   rsearch_CMProbifyEmissions(CM_t *cm, fullmat_t *fullmat);
extern void  CMLogoddsify(CM_t *cm);
extern int   CMCountStatetype(CM_t *cm, char type);
extern int   CMCountNodetype(CM_t *cm, char type);
extern int   CMSegmentCountStatetype(CM_t *cm, int r, int z, char type);
extern int   CMSubtreeCountStatetype(CM_t *cm, int v, char type);
extern int   CMSubtreeCountNodetype(CM_t *cm, int v, char type);
extern int   CMSubtreeFindEnd(CM_t *cm, int v);
extern int   CalculateStateIndex(CM_t *cm, int node, char utype);
extern int   TotalStatesInNode(int ndtype);
extern int   SplitStatesInNode(int ndtype);
extern int   InsertStatesInNode(int ndtype);
extern int   StateDelta(int sttype);
extern int   StateLeftDelta(int sttype);
extern int   StateRightDelta(int sttype);
extern int   Emitmode(int sttype);
extern void  PrintCM(FILE *fp, CM_t *cm);
extern void  SummarizeCM(FILE *fp, CM_t *cm);
extern char *Statetype(int type);
extern int   StateCode(char *s);
extern char *Nodetype(int type);
extern int   NodeCode(char *s);
extern char *UniqueStatetype(int type);
extern int   UniqueStateCode(char *s);
extern int   DeriveUniqueStateCode(int ndtype, int sttype);
extern int   StateMapsLeft(char st);
extern int   StateMapsRight(char st);
extern int   StateMapsMatch(char st);
extern int   StateMapsInsert(char st);
extern int   StateMapsDelete(char st);
extern int   NodeMapsLeft(char ndtype);
extern int   NodeMapsRight(char ndtype);
extern int   StateIsDetached(CM_t *cm, int v);
extern CM_t *CMRebalance(CM_t *cm);
extern int **IMX2Alloc(int rows, int cols);
extern void  IMX2Free(int **mx);
extern float rsearch_calculate_gap_penalty (char from_state, char to_state, int from_node, int to_node, float input_alpha, float input_beta, float input_alphap, float input_betap);
extern int   ExponentiateCM(CM_t *cm, double z);
extern void  cm_banner(FILE *fp, char *progname, char *banner);
extern void  cm_CalcExpSc(CM_t *cm, float **ret_expsc, float **ret_expsc_noss);
extern int   cm_Validate(CM_t *cm, float tol, char *errbuf);
extern char *CMStatetype(char st);
extern char *CMNodetype(char nd);
extern char *CMStateid(char st);
extern int   cm_SetName(CM_t *cm, char *name);
extern int   cm_SetAccession(CM_t *cm, char *acc);
extern int   cm_SetDescription(CM_t *cm, char *desc);
extern int   cm_AppendComlog(CM_t *cm, int argc, char **argv);
extern int   cm_SetCtime(CM_t *cm);
extern int   DefaultNullModel(const ESL_ALPHABET *abc, float **ret_null);
extern int   CMAllocNullModel(CM_t *cm);
extern void  CMSetNullModel(CM_t *cm, float *null);
extern int   CMReadNullModel(const ESL_ALPHABET *abc, char *nullfile, float **ret_null);
extern int   IntMaxDigits();
extern int   IntDigits(int i);
extern ComLog_t * CreateComLog();
extern void       FreeComLog(ComLog_t *clog);
extern int        CopyComLog(const ComLog_t *src, ComLog_t *dest);
extern int        cm_GetAvgHitLen(CM_t *cm, char *errbuf, float *ret_avg_hit_len);
extern int        CompareCMGuideTrees(CM_t *cm1, CM_t *cm2);

/* from dispatch.c */
extern int DispatchSearch    (CM_t *cm, char *errbuf, int fround, ESL_DSQ *dsq, int i0, int j0, 
			      search_results_t **results, float size_limit, int *ret_flen, float *ret_sc);
extern int DispatchAlignments(CM_t *cm, char *errbuf, seqs_to_aln_t *seqs_to_aln, ESL_DSQ *dsq, search_results_t *results, 
			      int first_result, int bdump_level, int debug_level, int silent_mode, int do_null3, ESL_RANDOMNESS *r, float size_limit, FILE *ofp,
			      int pad7, int len7, float sc7, int end7, float mprob7, float mcprob7, float iprob7, float ilprob7);

/* from cm_dpalign.c */
extern int FastAlignHB        (CM_t *cm, char *errbuf, ESL_RANDOMNESS *r, ESL_DSQ *dsq, int L, int i0, int j0, float size_limit, CM_HB_MX *mx, CM_HB_SHADOW_MX *shmx, int do_optacc, int do_sample, CM_HB_MX *post_mx, Parsetree_t **ret_tr, char **ret_pcode1, char **ret_pcode2, float *ret_sc, float *ret_ins_sc);
extern int FastAlign          (CM_t *cm, char *errbuf, ESL_RANDOMNESS *r, ESL_DSQ *dsq, int L, int i0, int j0, float size_limit, float ****ret_mx, int do_optacc, int do_sample, float ****ret_post_mx, Parsetree_t **ret_tr, char **ret_pcode1, char **ret_pcode2, float *ret_sc, float *ret_ins_sc);
extern int FastInsideAlignHB  (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, float size_limit, CM_HB_MX *mx,     float *ret_sc);
extern int FastInsideAlign    (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, float size_limit, float ****ret_mx, float *ret_sc);
extern int FastOutsideAlignHB (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, float size_limit, CM_HB_MX *mx,    CM_HB_MX *ins_mx, int do_check, float *ret_sc);
extern int FastOutsideAlign   (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, float size_limit, float ****ret_mx, float ***ins_mx, int do_check, float *ret_sc);
extern int SampleFromInsideHB (ESL_RANDOMNESS *r, CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, CM_HB_MX *mx, Parsetree_t **ret_tr, float *ret_sc);
extern int SampleFromInside   (ESL_RANDOMNESS *r, CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float ***mx,  Parsetree_t **ret_tr, float *ret_sc);
extern int   ** alloc_jdbanded_vjd_kshadow_deck(int L, int i, int j, int jmin, int jmax, int *hdmin, int *hdmax);
extern char  ** alloc_jdbanded_vjd_yshadow_deck(int L, int i, int j, int jmin, int jmax, int *hdmin, int *hdmax);
extern int  CMPostalCode(CM_t *cm, char *errbuf, int i0, int j0, float ***post, Parsetree_t *tr, int do_marginalize, char **ret_pcode1, char **ret_pcode2, float *ret_avgp);
extern int  CMPostalCodeHB(CM_t *cm, char *errbuf, int i0, int j0, CM_HB_MX *post_mx, Parsetree_t *tr, int do_marginalize, char **ret_pcode1, char **ret_pcode2, float *ret_avgp);
extern float FScore2Prob(float sc, float null);
extern int   Fscore2postcode(float sc);
extern int CMPosteriorHB      (CM_t *cm, char *errbuf, int i0, int j0, float size_limit, CM_HB_MX *ins_mx, CM_HB_MX *out_mx, CM_HB_MX *post_mx);
extern int CMPosterior        (CM_t *cm, char *errbuf, int i0, int j0, float size_limit, float ***ins_mx, float ***out_mx, float ***post_mx);
extern int CMCheckPosteriorHB (CM_t *cm, char *errbuf, int i0, int j0, CM_HB_MX *post);
extern int CMCheckPosterior   (CM_t *cm, char *errbuf, int i0, int j0, float ***post);

/* from cm_dpsearch.c */
extern int  FastCYKScan      (CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq, int i0, int j0, float cutoff, search_results_t *results, int do_null3, float **ret_vsc, float *ret_sc);
extern int  FastIInsideScan  (CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq, int i0, int j0, float cutoff, search_results_t *results, int do_null3, float **ret_vsc, float *ret_sc);
extern int  XFastIInsideScan (CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq, int i0, int j0, float cutoff, search_results_t *results, int do_null3, float **ret_vsc, float *ret_sc);
extern int  X2FastIInsideScan(CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq, int i0, int j0, float cutoff, search_results_t *results, int do_null3, float **ret_vsc, float *ret_sc);
extern int  FastFInsideScan  (CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq, int i0, int j0, float cutoff, search_results_t *results, int do_null3, float **ret_vsc, float *ret_sc);
extern int  RefCYKScan       (CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq, int i0, int j0, float cutoff, search_results_t *results, int do_null3, float **ret_vsc, float *ret_sc);
extern int  RefIInsideScan   (CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq, int i0, int j0, float cutoff, search_results_t *results, int do_null3, float **ret_vsc, float *ret_sc);
extern int  XRefIInsideScan  (CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq, int i0, int j0, float cutoff, search_results_t *results, int do_null3, float **ret_vsc, float *ret_sc);
extern int  RefFInsideScan   (CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq, int i0, int j0, float cutoff, search_results_t *results, int do_null3, float **ret_vsc, float *ret_sc);
extern int  rsearch_CYKScan  (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float cutoff, int D, search_results_t *results, float *ret_sc);
extern int  FastCYKScanHB    (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, float cutoff, search_results_t *results, int do_null3, CM_HB_MX *mx, float size_limit, float *ret_sc);
extern int  FastFInsideScanHB(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, float cutoff, search_results_t *results, int do_null3, CM_HB_MX *mx, float size_limit, float *ret_sc);
extern int  cm_CountSearchDPCalcs(CM_t *cm, char *errbuf, int L, int *dmin, int *dmax, int W, int correct_for_first_W, float **ret_vcalcs, float *ret_calcs);
extern int  ProcessSearchWorkunit(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, search_results_t **ret_results, float mxsize_limit, int my_rank, float **ret_survfractA, int **ret_nhitsA);
extern int  DetermineSeqChunksize(int nproc, int L, int W);

/* from cm_dpsmall.c */
extern float CYKDivideAndConquer(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, int j0, 
				 Parsetree_t **ret_tr, int *dmin, int *dmax);
extern float CYKInside(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, int j0, 
		       Parsetree_t **ret_tr, int *dmin, int *dmax);
extern float CYKInsideScore(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, 
			    int j0, int *dmin, int *dmax);
extern float CYKDemands(CM_t *cm, int L, int *dmin, int *dmax, int be_quiet);
extern void  debug_print_bands(FILE *fp, CM_t *cm, int *dmin, int *dmax);
/* cm_dpsmall.c: size calculators - not normally part of external API, but truncyk.c currently uses them */
extern float insideT_size(CM_t *cm, int L, int r, int z, int i0, int j0);
extern float vinsideT_size(CM_t *cm, int r, int z, int i0, int i1, int j1, int j0);
/* cm_dpsmall.c: the memory management routines. */
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

/* from cm_io.c */
extern CMFILE *CMFileOpen(char *cmfile, char *env);
extern int     CMFileRead(CMFILE *cmf, char *errbuf, ESL_ALPHABET **ret_abc, CM_t **ret_cm);
extern void    CMFileClose(CMFILE *cmf);
extern void    CMFileRewind(CMFILE *cmf);
extern int     CMFilePositionByIndex(CMFILE *cmf, int64_t idx);
extern int     CMFilePositionByKey(CMFILE *cmf, char *key);
extern int     CMFileWrite(FILE *fp, CM_t *cm, int do_binary, char *errbuf);

/* from cm_modelconfig.c */
extern int   ConfigCM(CM_t *cm, char *errbuf, int always_calc_W, CM_t *mother_cm, CMSubMap_t *mother_map);
extern void  ConfigLocal(CM_t *cm, float p_internal_start, float p_internal_exit);
extern void  ConfigGlobal(CM_t *cm);
extern void  ConfigNoLocalEnds(CM_t *cm);
extern void  ConfigLocalEnds(CM_t *cm, float p_internal_exit);
extern void  ConfigLocal_fullsub(CM_t *cm, float p_internal_start, 
				float p_internal_exit, int sstruct_nd,
				int estruct_nd);
extern void  ConfigLocal_DisallowELEmissions(CM_t *cm);
extern int   ConfigQDBAndW(CM_t *cm, int do_calc_qdb);

/* from cm_modelmaker.c */
extern int  HandModelmaker(ESL_MSA *msa, char *errbuf, int use_rf, float gapthresh, CM_t **ret_cm, Parsetree_t **ret_mtr);
extern int  ConsensusModelmaker(const ESL_ALPHABET *abc, char *errbuf, char *ss_cons, int clen, int building_sub_model, CM_t **ret_cm, Parsetree_t **ret_gtr);
extern Parsetree_t *Transmogrify(CM_t *cm, Parsetree_t *gtr, 
				 ESL_DSQ *dsq, char *aseq, int alen);
extern int  cm_from_guide(CM_t *cm, char *errbuf, Parsetree_t *gtr, int will_never_localize);
extern int  cm_find_and_detach_dual_inserts(CM_t *cm, int do_check, int do_detach);
extern int  cm_check_before_detaching(CM_t *cm, int insert1, int insert2);
extern int  cm_detach_state(CM_t *cm, int insert1, int insert2);
extern int  clean_cs(char *cs, int alen, int be_quiet);

/* from cm_mx.c */
extern CM_HB_MX *       cm_hb_mx_Create            (int M);
extern int              cm_hb_mx_GrowTo            (CM_t *cm, CM_HB_MX *mx, char *errbuf, CP9Bands_t *cp9b, int L, float size_limit);
extern int              cm_hb_mx_Dump              (FILE *ofp, CM_HB_MX *mx);
extern void             cm_hb_mx_Destroy           (CM_HB_MX *mx);
extern CM_HB_SHADOW_MX *cm_hb_shadow_mx_Create     (CM_t *cm, int M);
extern int              cm_hb_shadow_mx_GrowTo     (CM_t *cm, CM_HB_SHADOW_MX *mx, char *errbuf, CP9Bands_t *cp9b, int L);
extern int              cm_hb_shadow_mx_Dump       (FILE *ofp, CM_t *cm, CM_HB_SHADOW_MX *mx);
extern void             cm_hb_shadow_mx_Destroy    (CM_HB_SHADOW_MX *mx);
extern ScanMatrix_t *   cm_CreateScanMatrix        (CM_t *cm, int W, int *dmin, int *dmax, double beta_W, double beta_qdb, int do_banded, int do_float, int do_int);
extern int              cm_CreateScanMatrixForCM   (CM_t *cm, int do_float, int do_int);           
extern int              cm_FloatizeScanMatrix      (CM_t *cm, ScanMatrix_t *smx);
extern int              cm_IntizeScanMatrix        (CM_t *cm, ScanMatrix_t *smx);
extern int              cm_UpdateScanMatrixForCM   (CM_t *cm);
extern int              cm_FreeFloatsFromScanMatrix(CM_t *cm, ScanMatrix_t *smx);
extern int              cm_FreeIntsFromScanMatrix  (CM_t *cm, ScanMatrix_t *smx);
extern void             cm_FreeScanMatrix          (CM_t *cm, ScanMatrix_t *smx);
extern void             cm_FreeScanMatrixForCM     (CM_t *cm);
extern void             cm_DumpScanMatrixAlpha     (CM_t *cm, int j, int i0, int doing_float);
extern float **         FCalcOptimizedEmitScores   (CM_t *cm);
extern int **           ICalcOptimizedEmitScores   (CM_t *cm);
extern int **           ICopyOptimizedEmitScoresFromFloats(CM_t *cm, float **oesc);
extern void             DumpOptimizedEmitScores    (CM_t *cm, FILE *fp);
extern void             FreeOptimizedEmitScores    (float **fesc_vAA, int **iesc_vAA, int M);
extern float **         FCalcInitDPScores          (CM_t *cm);
extern int **           ICalcInitDPScores          (CM_t *cm);
extern GammaHitMx_t    *CreateGammaHitMx           (int L, int i0, int be_greedy, float cutoff, int do_backward);
extern void             FreeGammaHitMx             (GammaHitMx_t *gamma);
extern int              UpdateGammaHitMxCM         (CM_t *cm, char *errbuf, GammaHitMx_t *gamma, int j, float *alpha_row, int dn, int dx, int using_hmm_bands, int *bestr, search_results_t *results, int W, double **act);
extern int              UpdateGammaHitMxCP9Forward (CP9_t *cp9, char *errbuf, GammaHitMx_t *gamma, int i, int j, float hit_sc, search_results_t *results, int W, double **act, int clen);
extern int              UpdateGammaHitMxCP9Backward(CP9_t *cp9, char *errbuf, GammaHitMx_t *gamma, int i, int j, float hit_sc, search_results_t *results, int W, double **act);
extern void             TBackGammaHitMxForward     (GammaHitMx_t *gamma, search_results_t *results, int i0, int j0);
extern void             TBackGammaHitMxBackward    (GammaHitMx_t *gamma, search_results_t *results, int i0, int j0);

/* from cm_parsetree.c */
extern Parsetree_t *CreateParsetree(int size);
extern void         GrowParsetree(Parsetree_t *tr);
extern void         FreeParsetree(Parsetree_t *tr);
extern int          InsertTraceNode(Parsetree_t *tr, int y, int whichway, int emitl, int emitr, int state);
extern int          InsertTraceNodewithMode(Parsetree_t *tr, int y, int whichway, int emitl, int emitr, int state, int mode);
extern void         PrintParsetree(FILE *fp, Parsetree_t *tr);
extern void         ParsetreeCount(CM_t *cm, Parsetree_t *tr, ESL_DSQ *dsq, float wgt);
extern int          ParsetreeScore(CM_t *cm, CMEmitMap_t *emap, char *errbuf, Parsetree_t *tr, ESL_DSQ *dsq, int do_null2, float *ret_sc, float *ret_struct_sc, float *ret_primary_sc, int *ret_spos, int *ret_epos);
extern void         ParsetreeDump(FILE *fp, Parsetree_t *tr, CM_t *cm, ESL_DSQ *dsq, int *dmin, int *dmax);
extern int          ParsetreeCompare(Parsetree_t *t1, Parsetree_t *t2);
extern void         SummarizeMasterTrace(FILE *fp, Parsetree_t *tr);
extern void         MasterTraceDisplay(FILE *fp, Parsetree_t *mtr, CM_t *cm);
extern int          Parsetrees2Alignment(CM_t *cm, const ESL_ALPHABET *abc, ESL_SQ **sq, float *wgt, 
					 Parsetree_t **tr, int nseq, int do_full, int do_matchonly, ESL_MSA **ret_msa);
extern int          Alignment2Parsetrees(ESL_MSA *msa, CM_t *cm, Parsetree_t *mtr, char *errbuf, ESL_SQ ***ret_sq, Parsetree_t ***ret_tr);
extern float        ParsetreeScore_Global2Local(CM_t *cm, Parsetree_t *tr, ESL_DSQ *dsq, int print_flag);
extern int          Parsetree2CP9trace(CM_t *cm, Parsetree_t *tr, CP9trace_t **ret_cp9_tr);
extern void         rightjustify(const ESL_ALPHABET *abc, char *s, int n);
extern void         leftjustify(const ESL_ALPHABET *abc, char *s, int n);
extern int          EmitParsetree(CM_t *cm, char *errbuf, ESL_RANDOMNESS *r, char *name, int do_digital, Parsetree_t **ret_tr, ESL_SQ **ret_sq, int *ret_N);
extern int          ParsetreeScoreCorrectionNull2(CM_t *cm, char *errbuf, Parsetree_t *tr, ESL_DSQ *dsq, int start, float *ret_sc);
extern int          ParsetreeScoreCorrectionNull3(CM_t *cm, char *errbuf, Parsetree_t *tr, ESL_DSQ *dsq, int start, float *ret_sc);
extern int          ParsetreeCountMPEmissions(CM_t *cm, Parsetree_t *tr);
extern void         ScoreCorrectionNull3(const ESL_ALPHABET *abc, float *null0, float *comp, int len, float *ret_sc);
extern void         ScoreCorrectionNull3CompUnknown(const ESL_ALPHABET *abc, float *null0, ESL_DSQ *dsq, int start, int stop, float *ret_sc);

/* from cm_qdband.c */
extern void     BandExperiment(CM_t *cm);
extern double **BandDistribution(CM_t *cm, int W, int do_local);
extern int      BandCalculationEngine(CM_t *cm, int W, double p_thresh, 
				      int save_densities,
				      int **ret_dmin, int **ret_dmax, 
				      double ***ret_gamma, float **ret_seqlen);
extern int      BandTruncationNegligible(double *density, int b, int W, double *ret_beta);
extern int      BandMonteCarlo(CM_t *cm, int nsample, int W, double ***ret_gamma);
extern void     FreeBandDensities(CM_t *cm, double **gamma);
extern void     BandBounds(double **gamma, int M, int W, double p, 
			   int **ret_min, int **ret_max);
extern void     PrintBandGraph(FILE *fp, double **gamma, int *min, int *max, int v, int W);

extern void     PrintDPCellsSaved(CM_t *cm, int *min, int *max, int W);
extern void     ExpandBands(CM_t *cm, int qlen, int *dmin, int *dmax);
extern void     qdb_trace_info_dump(CM_t *cm, Parsetree_t *tr, int *dmin, int *dmax, int bdump_level);
extern int      cm_GetNCalcsPerResidueForGivenBeta(CM_t *cm, char *errbuf, int no_qdb, double beta, float *ret_cm_ncalcs_per_res, int *ret_W);

/* from cm_submodel.c */
extern int  build_sub_cm(CM_t *orig_cm, char *errbuf, CM_t **ret_cm, int sstruct, int estruct, CMSubMap_t **ret_submap, int print_flag);
extern void CP9NodeForPosn(CP9_t *hmm, int i0, int j0, int x, CP9_MX *post, int *ret_node, int *ret_type, float pmass, int is_start, int print_flag);
extern void StripWUSSGivenCC(ESL_MSA *msa, float gapthresh, int first_match, int last_match);
extern int  check_orig_psi_vs_sub_psi(CM_t *orig_cm, CM_t *sub_cm, CMSubMap_t *submap, double threshold, 
				       int print_flag);
extern int  check_sub_cm(CM_t *orig_cm, CM_t *sub_cm, CMSubMap_t *submap, CMSubInfo_t *subinfo, float pthresh, int print_flag);
extern int  check_sub_cm_by_sampling(CM_t *orig_cm, CM_t *sub_cm, ESL_RANDOMNESS *r, CMSubMap_t *submap, CMSubInfo_t *subinfo,
				     float chi_thresh, int nsamples, int print_flag);
extern int  check_sub_cm_by_sampling2(CM_t *orig_cm, CM_t *sub_cm, ESL_RANDOMNESS *r, int spos, int epos, int nseq);
extern int  sub_cm2cm_parsetree(CM_t *orig_cm, CM_t *sub_cm, Parsetree_t **ret_orig_tr, Parsetree_t *sub_tr, 
				CMSubMap_t *submap, int print_flag);
extern CMSubMap_t  *AllocSubMap(CM_t *sub_cm, CM_t *orig_cm, int sstruct, int estruct);
extern void         FreeSubMap(CMSubMap_t *submap);
extern CMSubInfo_t *AllocSubInfo(int clen);
extern void         FreeSubInfo(CMSubInfo_t *subinfo);
extern void  debug_print_cm_params(FILE *fp, CM_t *cm);
extern int   SubCMLogoddsify(CM_t *cm, char *errbuf, CM_t *mother_cm, CMSubMap_t *mother_map);
extern float ** SubFCalcAndCopyOptimizedEmitScoresFromMother(CM_t *cm, CM_t *mother_cm, CMSubMap_t *mother_map);

/* from cp9.c */
extern CP9_t *AllocCPlan9(int M, const ESL_ALPHABET *abc);
extern CP9_t *AllocCPlan9Shell(void);
extern void   AllocCPlan9Body(CP9_t *hmm, int M, const ESL_ALPHABET *abc);
extern void   FreeCPlan9(CP9_t *hmm);
extern void   ZeroCPlan9(CP9_t *hmm);
extern void   CPlan9SetNullModel(CP9_t *hmm, float *null, float p1);
extern void   DuplicateCP9(CM_t *src_cm, CM_t *dest_cm);
extern int    cp9_GetNCalcsPerResidue(CP9_t *cp9, char *errbuf, float *ret_cp9_ncalcs_per_res);

/* from cp9_dp.c */
extern int cp9_Viterbi(CM_t *cm, char *errbuf, CP9_MX *mx, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, search_results_t *results, 
		       int do_scan, int doing_align, int be_efficient, int do_null3, int **ret_psc, int *ret_maxres, 
		       CP9trace_t **ret_tr, float *ret_sc);
extern int cp9_ViterbiBackward(CM_t *cm, char *errbuf, CP9_MX *mx, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, search_results_t *results, 
			       int do_scan, int doing_align, int be_efficient, int do_null3, int **ret_psc, int *ret_maxres, 
			       CP9trace_t **ret_tr, float *ret_sc);
extern int cp9_Forward(CM_t *cm, char *errbuf, CP9_MX *mx, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, search_results_t *results, 
		       int do_scan, int doing_align, int be_efficient, int do_null3, int **ret_psc, int *ret_maxres, float *ret_sc);
extern int cp9_FastForward(CM_t *cm, char *errbuf, CP9_MX *mx, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, search_results_t *results, 
			   int do_scan, int doing_align, int be_efficient, int be_safe, int do_null3, int **ret_psc, int *ret_maxres, float *ret_sc);
extern int cp9_Backward(CM_t *cm, char *errbuf, CP9_MX *mx, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, search_results_t *results, 
			int do_scan, int doing_align, int be_efficient, int do_null3, int **ret_psc, int *ret_maxres, 
			float *ret_sc);
extern int cp9_CheckFB(CP9_MX *fmx, CP9_MX *bmx, CP9_t *hmm, char *errbuf, float sc, int i0, int j0, ESL_DSQ *dsq);
extern int cp9_WorstForward(CM_t *cm, char *errbuf, CP9_MX *mx, int thresh, int doing_scan, int doing_align, int *ret_L);
extern int cp9_CheckTransitionGuarantees(CP9_t *cp9, char *errbuf);
extern int cp9_GetLocalityMode(CP9_t *cp9, char *errbuf, int *ret_mode);

/* from cp9_modelconfig.c */
extern void  CP9Logoddsify(CP9_t *hmm);
extern void  CPlan9Renormalize(CP9_t *hmm);
extern void  CPlan9SWConfig(CP9_t *hmm, float pentry, float pexit, int do_match_local_cm, int first_cm_ndtype);
extern void  CPlan9ELConfig(CM_t *cm);
extern void  CPlan9NoEL(CM_t *cm);
extern void  CPlan9InitEL(CM_t *cm, CP9_t *cp9);
extern void  CPlan9RenormalizeExits(CP9_t *hmm, int spos);
extern int   Prob2Score(float p, float null);
extern float Score2Prob(int sc, float null);
extern float Scorify(int sc);
extern void  CPlan9CMLocalBeginConfig(CM_t *cm);
extern void  CP9_reconfig2sub(CP9_t *hmm, int spos, int epos, int spos_nd, int epos_nd, double **orig_phi);

/* from cp9_modelmaker.c */
extern CP9Map_t *AllocCP9Map(CM_t *cm);
extern void FreeCP9Map(CP9Map_t *cp9map);
extern int  build_cp9_hmm(CM_t *cm, CP9_t **ret_hmm, CP9Map_t **ret_cp9map, int do_psi_test,
			  float psi_vs_phi_threshold, int debug_level);
extern void CP9_map_cm2hmm(CM_t *cm, CP9Map_t *cp9map, int debug_level);
extern void fill_psi(CM_t *cm, double *psi, char ***tmap);
extern void fill_phi_cp9(CP9_t *hmm, double ***ret_phi, int spos, int entered_only);
extern void map_helper(CM_t *cm, CP9Map_t *cp9map, int k, int ks, int v);
extern void make_tmap(char ****ret_tmap);
extern int  CP9_check_by_sampling(CM_t *cm, CP9_t *hmm, ESL_RANDOMNESS *r, CMSubInfo_t *subinfo, int spos, int epos, 
				  float chi_thresh, int nsamples, int print_flag);
extern void debug_print_cp9_params(FILE *fp, CP9_t *hmm, int print_scores);
extern void debug_print_phi_cp9(CP9_t *hmm, double **phi);
extern int  MakeDealignedString(const ESL_ALPHABET *abc, char *aseq, int alen, char *ss, char **ret_s);
extern int  sub_build_cp9_hmm_from_mother(CM_t *cm, char *errbuf, CM_t *mother_cm, CMSubMap_t *mother_map, CP9_t **ret_hmm, CP9Map_t **ret_cp9map, int do_psi_test,
					  float psi_vs_phi_threshold, int debug_level);

/* from cp9_mx.c */
extern CP9_MX *CreateCP9Matrix(int N, int M);
extern void    FreeCP9Matrix  (CP9_MX *mx);
extern int     GrowCP9Matrix  (CP9_MX *mx, char *errbuf, int N, int M, int *kmin, int *kmax, int ***mmx, int ***imx, int ***dmx, int ***elmx, int **erow);
extern void    InitializeCP9Matrix(CP9_MX *mx);

/* from cp9_trace.c */
extern void  CP9AllocTrace(int tlen, CP9trace_t **ret_tr);
extern void  CP9ReallocTrace(CP9trace_t *tr, int tlen);
extern void  CP9FreeTrace(CP9trace_t *tr);
extern void  CP9_fake_tracebacks(ESL_MSA *msa, int *matassign, CP9trace_t ***ret_tr);
extern void  CP9TraceCount(CP9_t *hmm, ESL_DSQ *dsq, float wt, CP9trace_t *tr);
extern float CP9TraceScore(CP9_t *hmm, ESL_DSQ *dsq, CP9trace_t *tr);
extern char *CP9Statetype(char st);
extern void  CP9PrintTrace(FILE *fp, CP9trace_t *tr, CP9_t *hmm, ESL_DSQ *dsq);
extern int   CP9TransitionScoreLookup(CP9_t *hmm, char st1, int k1, 
				    char st2, int k2);
extern void  CP9ViterbiTrace(CP9_t *hmm, ESL_DSQ *dsq, int i0, int j0,
			     CP9_MX *mx, CP9trace_t **ret_tr);
extern void  CP9ReverseTrace(CP9trace_t *tr);
extern int   CP9Traces2Alignment(CM_t *cm, const ESL_ALPHABET *abc, ESL_SQ **sq, float *wgt, 
				 int nseq, CP9trace_t **tr, int do_full, int do_matchonly, ESL_MSA **ret_msa);
extern int   CP9TraceScoreCorrectionNull2(CP9_t *hmm, char *errbuf, CP9trace_t *tr, ESL_DSQ *dsq, int start, float *ret_sc);

/* from alphabet.c */
extern void   PairCount(const ESL_ALPHABET *abc, float *counters, ESL_DSQ syml, ESL_DSQ symr, float wt);
extern float  DegeneratePairScore(const ESL_ALPHABET *abc, float *esc, ESL_DSQ syml, ESL_DSQ symr);
extern int    iDegeneratePairScore(const ESL_ALPHABET *abc, int *esc, ESL_DSQ syml, ESL_DSQ symr);
extern char   resolve_degenerate (ESL_RANDOMNESS *r, char c);
extern int    revcomp(const ESL_ALPHABET *abc, ESL_SQ *comp, ESL_SQ *sq);
float  LeftMarginalScore(const ESL_ALPHABET *abc, float *esc, ESL_DSQ dres);
float  RightMarginalScore(const ESL_ALPHABET *abc, float *esc, ESL_DSQ dres);
extern float  FastPairScoreBothDegenerate(int K, float *esc, float *left, float *right);
extern int  iFastPairScoreBothDegenerate(int K, int *esc, float *left, float *right);
extern float FastPairScoreLeftOnlyDegenerate(int K, float *esc, float *left, ESL_DSQ symr);
extern int  iFastPairScoreLeftOnlyDegenerate(int K, int *iesc, float *left, ESL_DSQ symr);
extern float FastPairScoreRightOnlyDegenerate(int K, float *esc, float *right, ESL_DSQ syml);
extern float iFastPairScoreRightOnlyDegenerate(int K, int *iesc, float *right, ESL_DSQ syml);
extern float  FastPairScoreBothDegenerate(int K, float *esc, float *left, float *right);
extern int  iFastPairScoreBothDegenerate(int K, int *esc, float *left, float *right);
extern float FastPairScoreLeftOnlyDegenerate(int K, float *esc, float *left, ESL_DSQ symr);
extern int  iFastPairScoreLeftOnlyDegenerate(int K, int *iesc, float *left, ESL_DSQ symr);
extern float FastPairScoreRightOnlyDegenerate(int K, float *esc, float *right, ESL_DSQ syml);
extern float iFastPairScoreRightOnlyDegenerate(int K, int *iesc, float *right, ESL_DSQ syml);


/* from display.c */
extern Fancyali_t    *CreateFancyAli(const ESL_ALPHABET *abc, Parsetree_t *tr, CM_t *cm, CMConsensus_t *cons, ESL_DSQ *dsq, int do_noncanonical, char *pcode1, char *pcode2);
extern void           PrintFancyAli(FILE *fp, Fancyali_t *ali, int offset, int in_revcomp, int do_top);
extern void           FreeFancyAli(Fancyali_t *ali);
extern int            CreateCMConsensus(CM_t *cm, const ESL_ALPHABET *abc, float pthresh, float sthresh, CMConsensus_t **ret_cons);
extern void           FreeCMConsensus(CMConsensus_t *con);
extern void           MainBanner(FILE *fp, char *banner); 
extern int            IsCompensatory(const ESL_ALPHABET *abc, float *pij, int symi, int symj);
extern CMEmitMap_t   *CreateEmitMap(CM_t *cm); 
extern void           DumpEmitMap(FILE *fp, CMEmitMap_t *map, CM_t *cm);
extern void           FreeEmitMap(CMEmitMap_t *map);
extern void           FormatTimeString(char *buf, double sec, int do_frac);
extern int            GetDate(char *errbuf, char **ret_date);

/* from errors.c */
extern void cm_Die (char *format, ...);
extern void cm_Fail(char *format, ...);

/* from eweight.c */
extern int    cm_EntropyWeight(CM_t *cm, const Prior_t *pri, double etarget, int pretend_cm_is_hmm, double *ret_hmm_re, double *ret_Neff);
extern void   cm_Rescale(CM_t *hmm, float scale);
extern void   cp9_Rescale(CP9_t *hmm, float scale);
extern double cm_MeanMatchInfo(const CM_t *cm);
extern double cm_MeanMatchEntropy(const CM_t *cm);
extern double cm_MeanMatchRelativeEntropy(const CM_t *cm);
extern double cm_MeanMatchInfoHMM(const CM_t *cm);
extern double cm_MeanMatchEntropyHMM(const CM_t *cm);
extern double cm_MeanMatchRelativeEntropyHMM(const CM_t *cm);
extern double cp9_MeanMatchInfo(const CM_t *cm);
extern double cp9_MeanMatchEntropy(const CM_t *cm);
extern double cp9_MeanMatchRelativeEntropy(const CM_t *cm);

/* from hmmband.c */
extern int          cp9_HMM2ijBands(CM_t *cm, char *errbuf, CP9Bands_t *cp9b, CP9Map_t *cp9map, int i0, int j0, int doing_search, int debug_level);
extern int          cp9_HMM2ijBands_OLD(CM_t *cm, char *errbuf, CP9Bands_t *cp9b, CP9Map_t *cp9map, int i0, int j0, int doing_search, int debug_level);
extern CP9Bands_t * AllocCP9Bands(CM_t *cm, CP9_t *hmm);
extern void         FreeCP9Bands(CP9Bands_t *cp9bands);
extern int          cp9_Seq2Bands     (CM_t *cm, char *errbuf, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *pmx, ESL_DSQ *dsq, int i0, int j0, CP9Bands_t *cp9b, int doing_search, int debug_level);
extern int          cp9_Seq2Posteriors(CM_t *cm, char *errbuf, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *pmx, ESL_DSQ *dsq, int i0, int j0, int debug_level);
extern void         cp9_Posterior(ESL_DSQ *dsq, int i0, int j0, CP9_t *hmm, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *mx, int did_scan);
extern void         cp9_IFillPostSums(CP9_MX *post, CP9Bands_t *cp9, int i0, int j0);
extern double       DScore2Prob(int sc, float null);
extern int          cp9_FB2HMMBands        (CP9_t *hmm, char *errbuf, ESL_DSQ *dsq, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *pmx, CP9Bands_t *cp9b, 
				            int i0, int j0, int M, double p_thresh, int did_scan, int do_old_hmm2ij, int debug_level);
extern int          cp9_FB2HMMBandsWithSums(CP9_t *hmm, char *errbuf, ESL_DSQ *dsq, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *pmx, CP9Bands_t *cp9b, 
					    int i0, int j0, int M, double p_thresh, int did_scan, int do_old_hmm2ij, int debug_level);
extern int          HMMBandsEnforceValidParse(CM_t *cm, CP9Bands_t *cp9b, CP9Map_t *cp9map, char *errbuf, int i0, int j0, int doing_search, int *ret_did_expand, 
					      int **ret_r_mn, int **ret_r_mx, int **ret_r_in,  int **ret_r_ix, int **ret_r_dn, int **ret_r_dx,
					      int **ret_r_nn_i, int **ret_r_nx_i, int **ret_r_nn_j, int **ret_r_nx_j);
extern int          HMMBandsFixUnreachable(CP9Bands_t *cp9b, char *errbuf, int k, int r_prv_min, int r_prv_max, int r_insert_prv_min);
extern int          HMMBandsFillGap(CP9Bands_t *cp9b, char *errbuf, int k, int min1, int max1, int min2, int max2, int prv_nd_r_mn, int prv_nd_r_dn);
extern int          CMBandsCheckValidParse(CM_t *cm, CP9Bands_t *cp9b, char *errbuf, int i0, int j0, int doing_search);
/*extern void         cp9_RelaxRootBandsForSearch(CM_t *cm, int i0, int j0, int *imin, int *imax, int *jmin, int *jmax);*/
extern void         cp9_DebugPrintHMMBands(FILE *ofp, int L, CP9Bands_t *cp9b, double hmm_bandp, int debug_level);
extern void         cp9_CompareBands(CP9Bands_t *cp9b1, CP9Bands_t *cp9b2);
extern int          cp9_GrowHDBands(CP9Bands_t *cp9b, char *errbuf);
extern void         ijBandedTraceInfoDump(CM_t *cm, Parsetree_t *tr, int *imin, int *imax, 
					  int *jmin, int *jmax, int debug_level);
extern void         ijdBandedTraceInfoDump(CM_t *cm, Parsetree_t *tr, int *imin, int *imax, 
					   int *jmin, int *jmax, int **hdmin, int **hdmax, 
					   int debug_level);
extern int          cp9_ValidateBands(CM_t *cm, char *errbuf, CP9Bands_t *cp9b, int i0, int j0);
extern void         ij2d_bands(CM_t *cm, int L, int *imin, int *imax, int *jmin, int *jmax,
			       int **hdmin, int **hdmax, int debug_level);
extern void         combine_qdb_hmm_d_bands(CM_t *cm, int *jmin, int *jmax, int **hdmin, int **hdmax);
extern void         hd2safe_hd_bands(int M, int *jmin, int *jmax, int **hdmin, int **hdmax,
				     int *safe_hdmin, int *safe_hdmax);
extern void         debug_print_hd_bands(CM_t *cm, int **hdmin, int **hdmax, int *jmin, int *jmax);
extern void         PrintDPCellsSaved_jd(CM_t *cm, int *jmin, int *jmax, int **hdmin, int **hdmax, int W);
extern void         debug_print_ij_bands(CM_t *cm);
extern void         debug_print_parsetree_and_ij_bands(FILE *fp, Parsetree_t *tr, CM_t *cm, ESL_DSQ *dsq, CP9Bands_t *cp9b);

/* from p7_modelmaker.c */
#if 0
extern int          BuildP7HMM_MatchEmitsOnly(CM_t *cm, P7_HMM **ret_p7, P7_PROFILE **ret_gm, P7_OPROFILE **ret_om);
#endif
extern int          BuildP7HMM_MatchEmitsOnly(CM_t *cm, P7_HMM **ret_p7, P7_PROFILE **ret_gm);

/* from my_p7_msvfilter.c */
#if 0
extern int          my_p7_MSVFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, P7_GMX *gx, float *ret_sc);
extern int          p7_omx_CopyMSVRow2gmx(P7_OMX *ox, const P7_OPROFILE *om, P7_GMX *gx, int rowi, uint8_t xE, uint8_t xN, uint8_t xJ, uint8_t xB, uint8_t xC);
#endif
extern int          p7_gmx_Match2DMatrix(P7_GMX *gx, int do_diff, ESL_DMATRIX **ret_D, double *ret_min, double *ret_max);
extern int          my_dmx_Visualize(FILE *fp, ESL_DMATRIX *D, double min, double max, double min2fill);
extern int          my_p7_GTraceMSV(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_GMX *gx, P7_TRACE *tr, int **ret_i2k, int **ret_k2i, float **ret_sc, int **ret_iconflict);
extern int          Parsetree2i_to_k(CM_t *cm, CMEmitMap_t *emap, int L, char *errbuf, Parsetree_t *tr, int **ret_i2k);
extern int          prune_i2k(int *i2k, int *iconflict, float *isc, int L, double **phi, float min_sc, int min_len, int min_end, float min_mprob, float min_mcprob, float max_iprob, float max_ilprob);
extern int          p7_pins2bands(int *i2k, char *errbuf, int L, int M, int pad, int **ret_imin, int **ret_imax, int *ret_ncells);
extern int          DumpP7Bands(FILE *fp, int *i2k, int *kmin, int *kmax, int L);
extern int          cp9_ForwardP7B(CM_t *cm, char *errbuf, CP9_MX *mx, ESL_DSQ *dsq, int L, int *kmin, int *kmax, float *ret_sc);
extern int          cp9_ForwardP7B_OLD_WITH_EL(CM_t *cm, char *errbuf, CP9_MX *mx, ESL_DSQ *dsq, int L, int *kmin, int *kmax, float *ret_sc);
extern int          cp9_BackwardP7B(CM_t *cm, char *errbuf, CP9_MX *mx, ESL_DSQ *dsq, int L, int *kmin, int *kmax, float *ret_sc);
extern int          cp9_CheckFBP7B(CP9_MX *fmx, CP9_MX *bmx, CP9_t *hmm, char *errbuf, float sc, int i0, int j0, ESL_DSQ *dsq, int *kmin, int *kmax);
extern int          cp9_Seq2BandsP7B     (CM_t *cm, char *errbuf, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *pmx, ESL_DSQ *dsq, int L, CP9Bands_t *cp9b, int *kmin, int *kmax, int debug_level);
extern int          cp9_Seq2PosteriorsP7B(CM_t *cm, char *errbuf, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *pmx, ESL_DSQ *dsq, int L, int *kmin, int *kmax, int debug_level);
extern int          cp9_PosteriorP7B(ESL_DSQ *dsq, char *errbuf, int L, CP9_t *hmm, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *pmx, int *kmin, int *kmax);
extern int          cp9_FB2HMMBandsP7B(CP9_t *hmm, char *errbuf, ESL_DSQ *dsq, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *pmx, CP9Bands_t *cp9b, int L, int M, double p_thresh, int do_old_hmm2ij, int *kmin, int *kmax, int debug_level);
extern int          p7_Seq2Bands(CM_t *cm, char *errbuf, P7_GMX *gx, P7_BG *bg, P7_TRACE *p7_tr, ESL_DSQ *dsq, int L, 
				 double **phi, float sc7, int len7, int end7, float mprob7, float mcprob7, float iprob7, float ilprob7, int pad7,
				 int **ret_i2k, int **ret_kmin, int **ret_kmax, int *ret_ncells);

extern int          CP9NodeForPosnP7B(CP9_t *hmm, char *errbuf, int x, CP9_MX *post, int kn, int kx, int *ret_node, int *ret_type, int print_flag);
extern int          P7BandsAdjustForSubCM(int *kmin, int *kmax, int L, int spos, int epos);

/* from hybridsearch.c */
extern int                cm_cp9_HybridScan(CM_t *cm, char *errbuf, CP9_MX *mx, ESL_DSQ *dsq, HybridScanInfo_t *hsi, int i0, int j0, int W, 
					    float cutoff, search_results_t *results, int **ret_psc, int *ret_maxres, float *ret_sc);
extern int                predict_xsub(CM_t *cm, float *cm_vcalcs, float *cm_expsc, float *cp9_expsc);
extern void               cm_CalcAvgHitLength(CM_t *cm, double beta, float **ret_hitlen);
extern HybridScanInfo_t * cm_CreateHybridScanInfo(CM_t *cm, double hsi_beta, float full_cm_ncalcs);
extern int                cm_AddRootToHybridScanInfo(CM_t *cm, HybridScanInfo_t *hsi, int vroot_to_add);
extern int                cm_CheckCompatibleWithHybridScanInfo(CM_t *cm, HybridScanInfo_t *hsi, int v_root_to_add);
extern int                cm_ValidateHybridScanInfo(CM_t *cm, HybridScanInfo_t *hsi);
extern void               cm_FreeHybridScanInfo(HybridScanInfo_t *hsi, CM_t *cm);

/* from logsum.c */
extern void  init_ilogsum(void);
extern int   ILogsum(int s1, int s2);
extern int   ILogsumNI(int s1, int s2);
extern int   ILogsumNI_diff(int s1a, int s1b, int s2a, int s2b, int db);
extern void  FLogsumInit(void);
extern float LogSum2(float p1, float p2);
extern float FLogsum(float p1, float p2);

/* from mpisupport.c */
#if HAVE_MPI
extern int cm_master_MPIBcast(CM_t *cm, int tag, MPI_Comm comm, char **buf, int *nalloc);
extern int cm_worker_MPIBcast(int tag, MPI_Comm comm, char **buf, int *nalloc, ESL_ALPHABET **abc, CM_t **ret_cm);
extern int cm_MPIUnpack(ESL_ALPHABET **abc, char *buf, int n, int *pos, MPI_Comm comm, CM_t **ret_cm);
extern int cm_MPIPack(CM_t *cm, char *buf, int n, int *pos, MPI_Comm comm);
extern int cm_MPIPackSize(CM_t *cm, MPI_Comm comm, int *ret_n);
extern int cm_justread_MPIUnpack(ESL_ALPHABET **abc, char *buf, int n, int *pos, MPI_Comm comm, CM_t **ret_cm);
extern int cm_justread_MPIPack(CM_t *cm, char *buf, int n, int *pos, MPI_Comm comm);
extern int cm_justread_MPIPackSize(CM_t *cm, MPI_Comm comm, int *ret_n);
extern int cm_dsq_MPISend(ESL_DSQ *dsq, int L, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
extern int cm_dsq_MPIRecv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, ESL_DSQ **ret_dsq, int *ret_L);
extern int cm_search_results_MPISend(search_results_t *results, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
extern int cm_search_results_MPIRecv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, search_results_t  **ret_results);
extern int cm_search_results_MPIPackSize(const search_results_t *results, MPI_Comm comm, int *ret_n);
extern int cm_search_results_MPIPack(const search_results_t *results, char *buf, int n, int *position, MPI_Comm comm);
extern int cm_search_results_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, search_results_t **ret_results);
extern int cm_search_result_node_MPIPackSize(const search_result_node_t *rnode, MPI_Comm comm, int *ret_n) ;
extern int cm_search_result_node_MPIPack(const search_result_node_t *rnode, char *buf, int n, int *position, MPI_Comm comm);
extern int cm_search_result_node_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, search_result_node_t *ret_rnode);
extern int cm_seqs_to_aln_MPISend(seqs_to_aln_t *seqs_to_aln, int offset, int nseq_to_send, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
extern int cm_seqs_to_aln_MPIRecv(const ESL_ALPHABET *abc, int source, int tag, MPI_Comm comm, char **buf, int *nalloc, seqs_to_aln_t **ret_seqs_to_aln);
extern int cm_seqs_to_aln_MPIPackSize(const seqs_to_aln_t *results, int offset, int nseq_to_pack, MPI_Comm comm, int *ret_n);
extern int cm_seqs_to_aln_MPIPack(const seqs_to_aln_t *seqs_to_aln, int offset, int nseq_to_pack, char *buf, int n, int *position, MPI_Comm comm);
extern int cm_seqs_to_aln_MPIUnpack(const ESL_ALPHABET *abc, char *buf, int n, int *pos, MPI_Comm comm, seqs_to_aln_t **ret_seqs_to_aln);
extern int cm_parsetree_MPIPackSize(const Parsetree_t *tr, MPI_Comm comm, int *ret_n);
extern int cm_parsetree_MPIPack(const Parsetree_t *tr, char *buf, int n, int *position, MPI_Comm comm);
extern int cm_parsetree_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, Parsetree_t **ret_tr);
extern int cm_cp9trace_MPIPackSize(const CP9trace_t *cp9_tr, MPI_Comm comm, int *ret_n);
extern int cm_cp9trace_MPIPack(const CP9trace_t *cp9_tr, char *buf, int n, int *position, MPI_Comm comm);
extern int cm_cp9trace_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, CP9trace_t **ret_cp9_tr);
extern int cm_digitized_sq_MPIPackSize(const ESL_SQ *sq, MPI_Comm comm, int *ret_n);
extern int cm_digitized_sq_MPIPack(const ESL_SQ *sq, char *buf, int n, int *position, MPI_Comm comm);
extern int cm_digitized_sq_MPIUnpack(const ESL_ALPHABET *abc, char *buf, int n, int *pos, MPI_Comm comm, ESL_SQ **ret_sq);
extern int cmstats_MPIPackSize(CMStats_t *cmstats, MPI_Comm comm, int *ret_n);
extern int cmstats_MPIPack(CMStats_t *cmstats, char *buf, int n, int *position, MPI_Comm comm);
extern int cmstats_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, CMStats_t **ret_cmstats);
extern int exp_info_MPIPackSize(ExpInfo_t *exp, MPI_Comm comm, int *ret_n);
extern int exp_info_MPIPack(ExpInfo_t *exp, char *buf, int n, int *position, MPI_Comm comm);
extern int exp_info_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, ExpInfo_t **ret_exp);
extern int hmm_filter_info_MPIPackSize(HMMFilterInfo_t *hfi, MPI_Comm comm, int *ret_n);
extern int hmm_filter_info_MPIPack(HMMFilterInfo_t *hfi, char *buf, int n, int *position, MPI_Comm comm);
extern int hmm_filter_info_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, HMMFilterInfo_t **ret_hfi);
extern int best_filter_info_MPIPackSize(BestFilterInfo_t *bf, MPI_Comm comm, int *ret_n);
extern int best_filter_info_MPIPack(BestFilterInfo_t *bf, char *buf, int n, int *position, MPI_Comm comm);
extern int best_filter_info_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, BestFilterInfo_t **ret_bf);
extern int cmcalibrate_exp_results_MPIPackSize(float *scA, int nseq, MPI_Comm comm, int *ret_n);
extern int cmcalibrate_exp_results_MPIPack(float *scA, int nseq, char *buf, int n, int *position, MPI_Comm comm);
extern int cmcalibrate_exp_results_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, float **ret_scA, int *ret_nseq);
extern int cmcalibrate_filter_results_MPIPackSize(int nseq, MPI_Comm comm, int *ret_n);
extern int cmcalibrate_filter_results_MPIPack(float *cyk_scA, float *ins_scA, float *fwd_scA, int *partA, int nseq, char *buf, int n, int *position, MPI_Comm comm);
extern int cmcalibrate_filter_results_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, float **ret_cyk_scA, float **ret_ins_scA, float **ret_fwd_scA, int **ret_partA, int *ret_nseq);
extern int comlog_MPIPackSize(ComLog_t *comlog, MPI_Comm comm, int *ret_n);
extern int comlog_MPIPack    (ComLog_t *comlog, char *buf, int n, int *position, MPI_Comm comm);
extern int comlog_MPIUnpack  (char *buf, int n, int *pos, MPI_Comm comm, ComLog_t **ret_comlog);

#endif

/* from prior.c */
extern Prior_t *Prior_Create(void);
extern void     Prior_Destroy(Prior_t *pri);
extern Prior_t *Prior_Read(FILE *fp);
extern void     PriorifyCM(CM_t *cm, const Prior_t *pri);
extern Prior_t *Prior_Default(void);
extern struct p7prior_s *P7DefaultInfernalPrior(void);

/* from rnamat.c */
extern int numbered_nucleotide (char c);
extern int numbered_basepair (char c, char d);
extern FILE *MatFileOpen (char *matfile);
extern fullmat_t *ReadMatrix(const ESL_ALPHABET *abc, FILE *matfp);
extern int ribosum_MSA_resolve_degeneracies(fullmat_t *fullmat, ESL_MSA *msa);
extern int ribosum_calc_targets(fullmat_t *fullmat);
extern void FreeMat(fullmat_t *fullmat);

/* from searchinfo.c */
extern int  CreateSearchInfo(CM_t *cm, int cutoff_type, float sc_cutoff, float e_cutoff);
extern int  AddFilterToSearchInfo(CM_t *cm, int cyk_filter, int inside_filter, int viterbi_filter, int forward_filter, int hybrid_filter, 
				  ScanMatrix_t *smx, HybridScanInfo_t *hsi, int cutoff_type, float sc_cutoff, float e_cutoff, int do_null3);
extern void FreeSearchInfo(SearchInfo_t *si, CM_t *cm);
extern void DumpSearchInfo(SearchInfo_t *si);
extern void DumpSearchOpts(int search_opts);
extern void ValidateSearchInfo(CM_t *cm, SearchInfo_t *fi);
extern void UpdateSearchInfoCutoff(CM_t *cm, int nround, int cutoff_type, float sc_cutoff, float e_cutoff);
extern void UpdateSearchInfoForExpMode(CM_t *cm, int round, int exp_mode);
extern void UpdateSearchInfoForNewSMX(CM_t *cm);

extern search_results_t *CreateResults (int size);
extern void ExpandResults              (search_results_t *r, int additional);
extern void AppendResults              (search_results_t *src_results, search_results_t *dest_results, int i0);
extern void FreeResults                (search_results_t *r);
extern int  CompareResultsByScore      (const void *a_void, const void *b_void);
extern int  CompareResultsByEndPoint   (const void *a_void, const void *b_void);
extern void SortResultsByScore         (search_results_t *results);
extern void SortResultsByEndPoint      (search_results_t *results);
extern void PrintResults               (CM_t *cm, FILE *fp, FILE *tabfp, SearchInfo_t *si, const ESL_ALPHABET *abc, CMConsensus_t *cons, dbseq_t *dbseq, int do_top, int do_bottom, int do_noncompensatory, int do_noncanonical, int namewidth);
extern void ReportHit                  (int i, int j, int bestr, float score, search_results_t *results);
extern int  UpdateHitScoresWithNull2Or3(CM_t *cm, char *errbuf, SearchInfo_t *si, search_results_t *results, ESL_DSQ *dsq, int first_result, float sc_cutoff, int do_null2, int do_null3, int sort_by_score, int sort_by_endpoint);
extern void RemoveOverlappingHits      (search_results_t *results, int i0, int j0);
extern int  RemoveHitsOverECutoff      (CM_t *cm, char *errbuf, SearchInfo_t *si, int sround, search_results_t *results, ESL_DSQ *dsq, int first_result, int sort_by_score, int sort_by_endpoint);
extern int  ScoresFromResults          (search_results_t *results, char *errbuf, float **ret_scA, int *ret_scN); 
extern float CountScanDPCalcs          (CM_t *cm, int L, int use_qdb);
extern BestFilterInfo_t *CreateBestFilterInfo();
extern int  SetBestFilterInfoHMM(BestFilterInfo_t *bf, char *errbuf, int cm_M, float cm_eval, float F, int N, int db_size, float full_cm_ncalcs, int ftype, float e_cutoff, float fil_ncalcs, float fil_plus_surv_ncalcs);
extern int  SetBestFilterInfoHybrid(BestFilterInfo_t *bf, char *errbuf, int cm_M, float cm_eval, float F, int N, int db_size, float full_cm_ncalcs, float e_cutoff, float fil_ncalcs, float fil_plus_surv_ncalcs, HybridScanInfo_t *hsi, int np, ExpInfo_t **hexpA);
extern void FreeBestFilterInfo(BestFilterInfo_t *bf);
extern void DumpBestFilterInfo(BestFilterInfo_t *bf);
extern HMMFilterInfo_t *CreateHMMFilterInfo();
extern int  SetHMMFilterInfoHMM(HMMFilterInfo_t *hfi, char *errbuf, float F, int N, int dbsize, int ncut, float *cm_E_cut, float *fwd_E_cut, int always_better_than_Smax);
extern void FreeHMMFilterInfo(HMMFilterInfo_t *hfi);
extern int  DumpHMMFilterInfo(FILE *fp, HMMFilterInfo_t *hfi, char *errbuf, CM_t *cm, int cm_mode, int hmm_mode, long dbsize, int ncm, int namewidth, char *namedashes);
extern int  DumpHMMFilterInfoForCME(FILE *fp, HMMFilterInfo_t *hfi, char *errbuf, CM_t *cm, int cm_mode, int hmm_mode, long dbsize, int cmi, float cm_E, int do_header, int namewidth, char *namedashes,
				    float *ret_cm_bit_sc, float *ret_hmm_E, float *ret_hmm_bit_sc, float *ret_S, float *ret_xhmm, float *ret_spdup, float *ret_cm_ncalcs_per_res, float *ret_hmm_ncalcs_per_res, int *ret_do_filter);
extern int  DumpHMMFilterInfoForCMBitScore(FILE *fp, HMMFilterInfo_t *hfi, char *errbuf, CM_t *cm, int cm_mode, int hmm_mode, long dbsize, int cmi, float cm_bit_sc, int do_header, int namewidth, char *namedashes,
					   float *ret_cm_E, float *ret_hmm_E, float *ret_hmm_bit_sc, float *ret_S, float *ret_xhmm, float *ret_spdup, float *ret_cm_ncalcs_per_res, float *ret_hmm_ncalcs_per_res, int *ret_do_filter);
extern int   PlotHMMFilterInfo(FILE *fp, HMMFilterInfo_t *hfi, char *errbuf, CM_t *cm, int cm_mode, int hmm_mode, long dbsize, int mode);
extern float GetHMMFilterS(HMMFilterInfo_t *hfi, int cut, int W, float avg_hit_len);
extern float GetHMMFilterTotalCalcs(HMMFilterInfo_t *hfi, int cut, int W, float avg_hit_len, float cm_ncalcs_per_res, float hmm_ncalcs_per_res);
extern float GetHMMFilterXHMM(HMMFilterInfo_t *hfi, int cut, int W, float avg_hit_len, float cm_ncalcs_per_res, float hmm_ncalcs_per_res);
extern float GetHMMFilterSpeedup(HMMFilterInfo_t *hfi, int cut, int W, float avg_hit_len, float cm_ncalcs_per_res, float hmm_ncalcs_per_res);
extern int   GetHMMFilterFwdECutGivenCME(HMMFilterInfo_t *hfi, char *errbuf, float cm_E, long dbsize, int *ret_cut_pt);
extern int   GetHMMFilterFwdECutGivenCMBitScore(HMMFilterInfo_t *hfi, char *errbuf, float cm_bit_sc, long dbsize, int *ret_cut_pt, CM_t *cm, int cm_mode);
extern float SurvFract2E(float S, int W, float avg_hit_len, long dbsize);
extern float E2SurvFract(float E, int W, float avg_hit_len, long dbsize, int do_pad);
extern int   Results2SurvFract(CM_t *cm, char *errbuf, int i0, int j0, search_results_t *results, int do_pad, int do_collapse, float *ret_survfract);

/* from seqstoaln.c */
extern seqs_to_aln_t *CreateSeqsToAln(int size, int i_am_mpi_master);
extern seqs_to_aln_t *CreateSeqsToAlnFromSq(ESL_SQ **sq, int size, int i_am_mpi_master);
extern int            GrowSeqsToAln(seqs_to_aln_t *seqs_to_aln, int new_alloc, int i_am_mpi_master); 
extern void           FreeSeqsToAln(seqs_to_aln_t *seqs_to_aln);
extern void           FreePartialSeqsToAln(seqs_to_aln_t *s, int do_free_sq, int do_free_tr, int do_free_cp9_tr, int do_free_post, int do_free_sc, int do_free_pp, int do_free_struct_sc);
extern int            ReadSeqsToAln(const ESL_ALPHABET *abc, ESL_SQFILE *seqfp, int nseq, int do_read_all, seqs_to_aln_t *seqs_to_aln, int i_am_mpi_master); 
extern seqs_to_aln_t *CMEmitSeqsToAln(ESL_RANDOMNESS *r, CM_t *cm, int ncm, int nseq, int padW, double *pdist, int i_am_mpi_master);
extern seqs_to_aln_t *RandomEmitSeqsToAln(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, double *pdist, int extranum, int nseq, double *L_distro, int Lmax, int i_am_mpi_master); 

/* from stats.c */
extern CMStats_t *AllocCMStats(int np);
extern void       FreeCMStats(CMStats_t *cmstats);
extern int        debug_print_cmstats(CM_t *cm, char *errbuf, CMStats_t *cmstats, int has_fthr);
extern int        debug_print_expinfo(ExpInfo_t *exp);
extern int        get_gc_comp(const ESL_ALPHABET *abc, ESL_DSQ *dsq, int start, int stop);
extern int        get_alphabet_comp(const ESL_ALPHABET *abc, ESL_DSQ *dsq, int start, int stop, float **ret_freq); 
extern int        GetDBSize (ESL_SQFILE *sqfp, char *errbuf, long *ret_N, int *ret_nseq, int *ret_namewidth);
extern int        GetDBInfo(const ESL_ALPHABET *abc, ESL_SQFILE *sqfp, char *errbuf, long *ret_N, int *ret_nseq, double **ret_gc_ct);
extern int        E2MinScore(CM_t *cm, char *errbuf, int exp_mode, float E,  float *ret_sc);
extern int        E2ScoreGivenExpInfo(ExpInfo_t *exp, char *errbuf, float E, float *ret_sc);
extern int        Score2MaxE(CM_t *cm, char *errbuf, int exp_mode, float sc, float *ret_E);
extern double     Score2E(float x, double mu, double lambda, long eff_dbsize);
extern int        CM2ExpMode(CM_t *cm, int search_opts, int *ret_cm_exp_mode, int *ret_cp9_exp_mode);
extern int        CM2FthrMode(CM_t *cm, char *errbuf, int search_opts, int *ret_fthr_mode);
extern int        ExpModeIsLocal(int exp_mode);
extern int        ExpModeIsForCM(int exp_mode);
extern int        ExpModeToFthrMode(int exp_mode);
extern ExpInfo_t *CreateExpInfo();
extern void       SetExpInfo(ExpInfo_t *exp, double lambda, double mu_orig, long dbsize, int nrandhits, double tailp);
extern ExpInfo_t *DuplicateExpInfo(ExpInfo_t *src);
extern char      *DescribeExpMode(int exp_mode);
extern char      *DescribeFthrMode(int fthr_mode);
extern int        UpdateExpsForDBSize(CM_t *cm, char *errbuf, long dbsize);

/* from truncyk.c */
void  SetMarginalScores(CM_t *cm);
float TrCYK_DnC(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, int j0, Parsetree_t **ret_tr);
float TrCYK_Inside(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, int j0, int lenCORREX, Parsetree_t **ret_tr);
/* legacy, avoid use: */
float trinside (CM_t *cm, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0, int do_full,
                void ****ret_shadow, void ****ret_L_shadow, void ****ret_R_shadow,
                void ****ret_T_shadow, void ****ret_Lmode_shadow, void ****ret_Rmode_shadow,
                int *ret_mode, int *ret_v, int *ret_i, int *ret_j);




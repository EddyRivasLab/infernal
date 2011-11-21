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
extern CM_t *CreateCM(int nnodes, int nstates, int clen, const ESL_ALPHABET *abc);
extern CM_t *CreateCMShell(void);
extern void  CreateCMBody(CM_t *cm, int nnodes, int nstates, int clen, const ESL_ALPHABET *abc);
extern void  CMZero(CM_t *cm);
extern void  CMRenormalize(CM_t *cm);
extern void  FreeCM(CM_t *cm);
extern void  CMSimpleProbify(CM_t *cm);
extern int   rsearch_CMProbifyEmissions(CM_t *cm, fullmat_t *fullmat);
extern int   CMLogoddsify(CM_t *cm);
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
extern int   NumReachableInserts(int stid);
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
extern char *MarginalMode(char mode); 
extern int   cm_SetName(CM_t *cm, char *name);
extern int   cm_SetAccession(CM_t *cm, char *acc);
extern int   cm_SetDescription(CM_t *cm, char *desc);
extern int   cm_SetConsensus(CM_t *cm, CMConsensus_t *cons, ESL_SQ *sq);
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
extern int        CloneCMJustReadFromFile(CM_t *cm, char *errbuf, CM_t **ret_cm);
extern void       DumpCMFlags(FILE *fp, CM_t *cm);
extern ESL_GETOPTS *cm_CreateDefaultApp(ESL_OPTIONS *options, int nargs, int argc, char **argv, char *banner, char *usage);
extern CM_P7_OM_BLOCK *cm_p7_oprofile_CreateBlock(int size);
extern void            cm_p7_oprofile_DestroyBlock(CM_P7_OM_BLOCK *block);
extern float **FCalcOptimizedEmitScores      (CM_t *cm);
extern int   **ICalcOptimizedEmitScores      (CM_t *cm);
extern int   **ICopyOptimizedEmitScoresFromFloats(CM_t *cm, float **oesc);
extern void    DumpOptimizedEmitScores       (CM_t *cm, FILE *fp);
extern void    FreeOptimizedEmitScores       (float **fesc_vAA, int **iesc_vAA, int M);
extern float **FCalcInitDPScores             (CM_t *cm);
extern int   **ICalcInitDPScores             (CM_t *cm);

/* from dispatch.c */
extern int DispatchSearch    (CM_t *cm, char *errbuf, int fround, ESL_DSQ *dsq, int i0, int j0, int hit_len_guess, 
			      CM_TOPHITS **hitlistA, float size_limit, int *ret_flen, float *ret_sc);
extern int DispatchAlignments(CM_t *cm, char *errbuf, seqs_to_aln_t *seqs_to_aln, 
			      int bdump_level, int debug_level, int silent_mode, int do_null3, TrScanInfo_t *trsi, ESL_RANDOMNESS *r, float size_limit, FILE *ofp, FILE *sfp, int iidx,
			      int pad7, int len7, float sc7, int end7, float mprob7, float mcprob7, float iprob7, float ilprob7);

/* from cm_dpalign.c */
extern int   cm_Align             (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, int do_optacc, int do_sample, CM_MX *mx,    CM_SHADOW_MX    *shmx, CM_MX    *post_mx, CM_EMIT_MX *emit_mx, ESL_RANDOMNESS *r, char **ret_ppstr, float *ret_ins_sc, Parsetree_t **ret_tr, float *ret_sc);
extern int   cm_AlignHB           (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, int do_optacc, int do_sample, CM_HB_MX *mx, CM_HB_SHADOW_MX *shmx, CM_HB_MX *post_mx, CM_HB_EMIT_MX *emit_mx, ESL_RANDOMNESS *r, char **ret_ppstr, float *ret_ins_sc, Parsetree_t **ret_tr, float *ret_sc);
extern int   cm_CYKInsideAlign    (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit,               CM_MX    *mx, CM_SHADOW_MX    *shmx, int *ret_b, float *ret_sc);
extern int   cm_CYKInsideAlignHB  (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit,               CM_HB_MX *mx, CM_HB_SHADOW_MX *shmx, int *ret_b, float *ret_sc);
extern int   cm_InsideAlign       (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit,               CM_MX    *mx, float *ret_sc);
extern int   cm_InsideAlignHB     (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit,               CM_HB_MX *mx, float *ret_sc);
extern int   cm_OptAccAlign       (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit,               CM_MX    *mx, CM_SHADOW_MX    *shmx, CM_EMIT_MX *emit_mx,    int *ret_b, float *ret_sc);
extern int   cm_OptAccAlignHB     (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit,               CM_HB_MX *mx, CM_HB_SHADOW_MX *shmx, CM_HB_EMIT_MX *emit_mx, int *ret_b, float *ret_sc);
extern int   cm_CYKOutsideAlign   (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, int do_check, CM_MX    *mx, CM_MX *inscyk_mx, float *ret_sc);
extern int   cm_CYKOutsideAlignHB (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, int do_check, CM_HB_MX *mx, CM_HB_MX *ins_mx, float *ret_sc);
extern int   cm_OutsideAlign      (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, int do_check, CM_MX    *mx, CM_MX    *ins_mx, float *ret_sc);
extern int   cm_OutsideAlignHB    (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, int do_check, CM_HB_MX *mx, CM_HB_MX *ins_mx, float *ret_sc);
extern int   cm_PosteriorHB       (CM_t *cm, char *errbuf,               int L, float size_limit,               CM_HB_MX *ins_mx, CM_HB_MX *out_mx, CM_HB_MX *post_mx);
extern int   cm_Posterior         (CM_t *cm, char *errbuf,               int L, float size_limit,               CM_MX    *ins_mx, CM_MX    *out_mx, CM_MX    *post_mx);
extern int   cm_EmitterPosterior  (CM_t *cm, char *errbuf,               int L, float size_limit, CM_MX    *post, CM_EMIT_MX    *emit_mx, int do_check);
extern int   cm_EmitterPosteriorHB(CM_t *cm, char *errbuf,               int L, float size_limit, CM_HB_MX *post, CM_HB_EMIT_MX *emit_mx, int do_check);
extern int   cm_SampleParsetree  (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, CM_MX    *mx, ESL_RANDOMNESS *r, Parsetree_t **ret_tr, float *ret_sc);
extern int   cm_SampleParsetreeHB(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, CM_HB_MX *mx, ESL_RANDOMNESS *r, Parsetree_t **ret_tr, float *ret_sc);
extern int   cm_PostCode  (CM_t *cm, char *errbuf, int L, CM_EMIT_MX    *emit_mx, Parsetree_t *tr, char **ret_ppstr, float *ret_avgp);
extern int   cm_PostCodeHB(CM_t *cm, char *errbuf, int L, CM_HB_EMIT_MX *emit_mx, Parsetree_t *tr, char **ret_ppstr, float *ret_avgp);
extern int   cm_InitializeOptAccShadowDZero  (CM_t *cm, char *errbuf, char ***yshadow, int L);
extern int   cm_InitializeOptAccShadowDZeroHB(CM_t *cm, CP9Bands_t *cp9b, char *errbuf, char ***yshadow, int L);
extern float FScore2Prob(float sc, float null);
extern char  Fscore2postcode(float sc);

/* from cm_dpalign_trunc.c */
extern int  cm_TrAlign              (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, char opt_mode, int do_optacc, int do_sample, CM_TR_MX    *mx, CM_TR_SHADOW_MX    *shmx, CM_TR_MX    *post_mx, CM_TR_EMIT_MX    *emit_mx, ESL_RANDOMNESS *r, char **ret_ppstr, float *ret_ins_sc, Parsetree_t **ret_tr, float *ret_sc);
extern int  cm_TrAlignHB            (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, char opt_mode, int do_optacc, int do_sample, CM_TR_HB_MX *mx, CM_TR_HB_SHADOW_MX *shmx, CM_TR_HB_MX *post_mx, CM_TR_HB_EMIT_MX *emit_mx, ESL_RANDOMNESS *r, char **ret_ppstr, float *ret_ins_sc, Parsetree_t **ret_tr, float *ret_sc);
extern int  cm_TrCYKInsideAlign     (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, char opt_mode, CM_TR_MX    *mx, CM_TR_SHADOW_MX    *shmx, int *ret_Jb, int *ret_Lb, int *ret_Rb, int *ret_Tb, char *ret_mode, float *ret_sc);
extern int  cm_TrCYKInsideAlignHB   (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, char opt_mode, CM_TR_HB_MX *mx, CM_TR_HB_SHADOW_MX *shmx, int *ret_Jb, int *ret_Lb, int *ret_Rb, int *ret_Tb, char *ret_mode, float *ret_sc);
extern int  cm_TrInsideAlign        (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, char opt_mode, CM_TR_MX    *mx, char *ret_mode, float *ret_sc);
extern int  cm_TrInsideAlignHB      (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, char opt_mode, CM_TR_HB_MX *mx, char *ret_mode, float *ret_sc);
extern int  cm_TrOptAccAlign        (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, char opt_mode, CM_TR_MX    *mx, CM_TR_SHADOW_MX    *shmx, CM_TR_EMIT_MX    *emit_mx, int *ret_b, float *ret_sc);
extern int  cm_TrOptAccAlignHB      (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, char opt_mode, CM_TR_HB_MX *mx, CM_TR_HB_SHADOW_MX *shmx, CM_TR_HB_EMIT_MX *emit_mx, int *ret_b, float *ret_sc);
extern int  cm_TrCYKOutsideAlign    (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, char opt_mode, int do_check, CM_TR_MX    *mx, CM_TR_MX    *inscyk_mx);
extern int  cm_TrCYKOutsideAlignHB  (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, char opt_mode, int do_check, CM_TR_HB_MX *mx, CM_TR_HB_MX *inscyk_mx);
extern int  cm_TrOutsideAlign       (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, char opt_mode, int do_check, CM_TR_MX *mx,    CM_TR_MX    *ins_mx);
extern int  cm_TrOutsideAlignHB     (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float size_limit, char opt_mode, int do_check, CM_TR_HB_MX *mx, CM_TR_HB_MX *ins_mx);
extern int  cm_TrPosterior          (CM_t *cm, char *errbuf,               int L, float size_limit, char opt_mode, CM_TR_MX    *ins_mx, CM_TR_MX    *out_mx, CM_TR_MX    *post_mx);
extern int  cm_TrPosteriorHB        (CM_t *cm, char *errbuf,               int L, float size_limit, char opt_mode, CM_TR_HB_MX *ins_mx, CM_TR_HB_MX *out_mx, CM_TR_HB_MX *post_mx);
extern int  cm_TrEmitterPosterior   (CM_t *cm, char *errbuf,               int L, float size_limit, char opt_mode, int do_check, CM_TR_MX    *post, CM_TR_EMIT_MX    *emit_mx);
extern int  cm_TrEmitterPosteriorHB (CM_t *cm, char *errbuf,               int L, float size_limit, char opt_mode, int do_check, CM_TR_HB_MX *post, CM_TR_HB_EMIT_MX *emit_mx);
extern int  cm_TrPostCode           (CM_t *cm, char *errbuf,               int L, CM_TR_EMIT_MX    *emit_mx, Parsetree_t *tr, char **ret_ppstr, float *ret_avgp);
extern int  cm_TrPostCodeHB         (CM_t *cm, char *errbuf,               int L, CM_TR_HB_EMIT_MX *emit_mx, Parsetree_t *tr, char **ret_ppstr, float *ret_avgp);
extern int  cm_TrFillFromMode       (char mode, int *ret_fill_L, int *ret_fill_R, int *ret_fill_T);

/* from cm_dpsearch.c */
extern int  FastCYKScan      (CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq, int64_t i0, int64_t j0, float cutoff, CM_TOPHITS *hitlist, int do_null3, float env_cutoff, int64_t *ret_envi, int64_t *ret_envj, float **ret_vsc, float *ret_sc);
extern int  RefCYKScan       (CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq, int64_t i0, int64_t j0, float cutoff, CM_TOPHITS *hitlist, int do_null3, float env_cutoff, int64_t *ret_envi, int64_t *ret_envj, float **ret_vsc, float *ret_sc);
extern int  FastIInsideScan  (CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq, int64_t i0, int64_t j0, float cutoff, CM_TOPHITS *hitlist, int do_null3, float **ret_vsc, float *ret_sc);
extern int  RefIInsideScan   (CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq, int64_t i0, int64_t j0, float cutoff, CM_TOPHITS *hitlist, int do_null3, float **ret_vsc, float *ret_sc);
extern int  FastFInsideScan  (CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq, int64_t i0, int64_t j0, float cutoff, CM_TOPHITS *hitlist, int do_null3, float **ret_vsc, float *ret_sc);
extern int  RefFInsideScan   (CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq, int64_t i0, int64_t j0, float cutoff, CM_TOPHITS *hitlist, int do_null3, float **ret_vsc, float *ret_sc);
extern int  FastCYKScanHB    (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int64_t i0, int64_t j0, float cutoff, CM_TOPHITS *hitlist, int do_null3, CM_HB_MX *mx, float size_limit, float env_cutoff, int64_t *ret_envi, int64_t *ret_envj, float *ret_sc);
extern int  FastFInsideScanHB(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int64_t i0, int64_t j0, float cutoff, CM_TOPHITS *hitlist, int do_null3, CM_HB_MX *mx, float size_limit, float env_cutoff, int64_t *ret_envi, int64_t *ret_envj, float *ret_sc);
extern int  cm_CountSearchDPCalcs(CM_t *cm, char *errbuf, int L, int *dmin, int *dmax, int W, int correct_for_first_W, float **ret_vcalcs, float *ret_calcs);
extern int  DetermineSeqChunksize(int nproc, int L, int W);

/* from cm_dpsearch_trunc.c */
extern int  RefTrCYKScan    (CM_t *cm, char *errbuf, TrScanMatrix_t *trsmx, TrScanInfo_t *trsi, ESL_DSQ *dsq, int64_t i0, int64_t j0, float cutoff, CM_TOPHITS *hitlist, 
			     int do_null3, float env_cutoff, int64_t *ret_envi, int64_t *ret_envj, float **ret_vsc, char *ret_mode, float *ret_sc);
extern int  RefITrInsideScan(CM_t *cm, char *errbuf, TrScanMatrix_t *trsmx, TrScanInfo_t *trsi, ESL_DSQ *dsq, int64_t i0, int64_t j0, float cutoff, CM_TOPHITS *hitlist,
			     int do_null3, float env_cutoff, int64_t *ret_envi, int64_t *ret_envj, float **ret_vsc, char *ret_mode, float *ret_sc);
extern int  TrCYKScanHB(CM_t *cm, char *errbuf, TrScanInfo_t *trsi, ESL_DSQ *dsq, int64_t i0, int64_t j0, float cutoff, CM_TOPHITS *hitlist, int do_null3, 
			CM_TR_HB_MX *mx, float size_limit, float env_cutoff, int64_t *ret_envi, int64_t *ret_envj, char *ret_mode, float *ret_sc);
extern int  FTrInsideScanHB(CM_t *cm, char *errbuf, TrScanInfo_t *trsi, ESL_DSQ *dsq, int64_t i0, int64_t j0, float cutoff, CM_TOPHITS *hitlist, int do_null3, 
			    CM_TR_HB_MX *mx, float size_limit, float env_cutoff, int64_t *ret_envi, int64_t *ret_envj, char *ret_mode, float *ret_sc);

/* from cm_dpsmall.c */
extern float CYKDivideAndConquer(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, int j0, 
				 Parsetree_t **ret_tr, int *dmin, int *dmax);
extern float CYKInside(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, int j0, 
		       Parsetree_t **ret_tr, int *dmin, int *dmax);
extern float CYKInsideScore(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, 
			    int j0, int *dmin, int *dmax);
extern float CYKDemands(CM_t *cm, int L, int *dmin, int *dmax, int be_quiet);
extern float CYKNonQDBSmallMbNeeded(CM_t *cm, int L);
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

/* from cm_file.c */
extern int     cm_file_Open(char *filename, char *env, int allow_1p0, CM_FILE **ret_cmfp, char *errbuf);
extern int     cm_file_OpenNoDB(char *filename, char *env, int allow_1p0, CM_FILE **ret_cmfp, char *errbuf);
extern int     cm_file_OpenBuffer(char *buffer, int size, int allow_1p0, CM_FILE **ret_cmfp);
extern void    cm_file_Close(CM_FILE *cmfp);
extern int     cm_file_CreateLock(CM_FILE *cmfp);
extern int     cm_file_WriteASCII(FILE *fp, int format, CM_t *cm);
extern int     cm_file_WriteBinary(FILE *fp, int format, CM_t *cm, off_t *opt_fp7_offset);
extern int     cm_file_Read(CM_FILE *cmfp, int read_fp7, ESL_ALPHABET **ret_abc, CM_t **opt_cm);
extern int     cm_file_PositionByKey(CM_FILE *cmfp, const char *key);
extern int     cm_file_Position(CM_FILE *cmfp, const off_t offset);
extern int     cm_p7_hmmfile_Read(CM_FILE *cmfp, ESL_ALPHABET *abc, off_t offset, P7_HMM **ret_hmm);
extern int     cm_p7_oprofile_Write(FILE *ffp, FILE *pfp, off_t cm_offset, int cm_len, int cm_W, float gfmu, float gflambda, P7_OPROFILE *om);
extern int     cm_p7_oprofile_ReadMSV(CM_FILE *cmfp, int read_scores, ESL_ALPHABET **byp_abc, off_t *ret_cm_offset, int *ret_cm_clen, int *ret_cm_W, float *ret_gfmu, float *ret_gflambda, P7_OPROFILE **ret_om);
extern int     cm_p7_oprofile_ReadBlockMSV(CM_FILE *cmfp, ESL_ALPHABET **byp_abc, CM_P7_OM_BLOCK *hmmBlock);
extern int     cm_p7_oprofile_Position(CM_FILE *cmfp, off_t offset);

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
extern int  HandModelmaker(ESL_MSA *msa, char *errbuf, int use_rf, int use_wts, float gapthresh, CM_t **ret_cm, Parsetree_t **ret_mtr);
extern int  ConsensusModelmaker(const ESL_ALPHABET *abc, char *errbuf, char *ss_cons, int clen, int building_sub_model, CM_t **ret_cm, Parsetree_t **ret_gtr);
extern Parsetree_t *Transmogrify(CM_t *cm, Parsetree_t *gtr, 
				 ESL_DSQ *dsq, char *aseq, int alen);
extern int  cm_from_guide(CM_t *cm, char *errbuf, Parsetree_t *gtr, int will_never_localize);
extern int  cm_find_and_detach_dual_inserts(CM_t *cm, int do_check, int do_detach);
extern int  cm_check_before_detaching(CM_t *cm, int insert1, int insert2);
extern int  cm_detach_state(CM_t *cm, int insert1, int insert2);
extern int  clean_cs(char *cs, int alen, int be_quiet);

/* from cm_mx.c */
extern CM_MX           *cm_mx_Create                  (CM_t *cm);
extern int              cm_mx_GrowTo                  (CM_t *cm, CM_MX *mx, char *errbuf, int L, float size_limit);
extern int              cm_mx_Dump                    (FILE *ofp, CM_MX *mx);
extern void             cm_mx_Destroy                 (CM_MX *mx);
extern int              cm_mx_SizeNeeded              (CM_t *cm, char *errbuf, int L, int64_t *ret_ncells, float *ret_Mb);

extern CM_TR_MX        *cm_tr_mx_Create               (CM_t *cm);
extern int              cm_tr_mx_GrowTo               (CM_t *cm, CM_TR_MX *mx, char *errbuf, int L, float size_limit);
extern int              cm_tr_mx_Dump                 (FILE *ofp, CM_TR_MX *mx, char opt_mode);
extern void             cm_tr_mx_Destroy              (CM_TR_MX *mx);
extern int              cm_tr_mx_SizeNeeded           (CM_t *cm, char *errbuf, int L, int64_t *ret_Jncells, int64_t *ret_Lncells, int64_t *ret_Rncells, int64_t *ret_Tncells, float *ret_Mb);

extern CM_HB_MX        *cm_hb_mx_Create               (int M);
extern int              cm_hb_mx_GrowTo               (CM_t *cm, CM_HB_MX *mx, char *errbuf, CP9Bands_t *cp9b, int L, float size_limit);
extern int              cm_hb_mx_Dump                 (FILE *ofp, CM_HB_MX *mx);
extern void             cm_hb_mx_Destroy              (CM_HB_MX *mx);
extern int              cm_hb_mx_SizeNeeded           (CM_t *cm, char *errbuf, CP9Bands_t *cp9b, int L, int64_t *ret_ncells, float *ret_Mb);

extern CM_TR_HB_MX     *cm_tr_hb_mx_Create            (CM_t *cm);
extern int              cm_tr_hb_mx_GrowTo            (CM_t *cm, CM_TR_HB_MX *mx, char *errbuf, CP9Bands_t *cp9b, int L, float size_limit);
extern int              cm_tr_hb_mx_Dump              (FILE *ofp, CM_TR_HB_MX *mx, char opt_mode);
extern void             cm_tr_hb_mx_Destroy           (CM_TR_HB_MX *mx);
extern int              cm_tr_hb_mx_SizeNeeded        (CM_t *cm, char *errbuf, CP9Bands_t *cp9b, int L, int64_t *ret_Jncells, int64_t *ret_Lncells, 
						       int64_t *ret_Rncells, int64_t *ret_Tncells, float *ret_Mb);

extern CM_SHADOW_MX    *cm_shadow_mx_Create           (CM_t *cm);
extern int              cm_shadow_mx_GrowTo           (CM_t *cm, CM_SHADOW_MX *mx, char *errbuf, int L, float size_limit);
extern int              cm_shadow_mx_Dump             (FILE *ofp, CM_t *cm, CM_SHADOW_MX *mx);
extern void             cm_shadow_mx_Destroy          (CM_SHADOW_MX *mx);
extern int              cm_shadow_mx_SizeNeeded       (CM_t *cm, char *errbuf, int L, int64_t *ret_ny_cells, int64_t *ret_nk_cells, float *ret_Mb);

extern CM_TR_SHADOW_MX *cm_tr_shadow_mx_Create        (CM_t *cm);
extern int              cm_tr_shadow_mx_GrowTo        (CM_t *cm, CM_TR_SHADOW_MX *mx, char *errbuf, int L, float size_limit);
extern int              cm_tr_shadow_mx_Dump          (FILE *ofp, CM_t *cm, CM_TR_SHADOW_MX *mx, char opt_mode);
extern void             cm_tr_shadow_mx_Destroy       (CM_TR_SHADOW_MX *mx);
extern int              cm_tr_shadow_mx_SizeNeeded    (CM_t *cm, char *errbuf, int L, int64_t *ret_Jny_cells, int64_t *ret_Lny_cells, int64_t *ret_Rny_cells, 
						       int64_t *ret_Jnk_cells, int64_t *ret_Lnk_cells, int64_t *ret_Rnk_cells, int64_t *ret_Tnk_cells, float *ret_Mb);

extern CM_HB_SHADOW_MX *cm_hb_shadow_mx_Create        (CM_t *cm);
extern int              cm_hb_shadow_mx_GrowTo        (CM_t *cm, CM_HB_SHADOW_MX *mx, char *errbuf, CP9Bands_t *cp9b, int L, float size_limit);
extern int              cm_hb_shadow_mx_Dump          (FILE *ofp, CM_t *cm, CM_HB_SHADOW_MX *mx);
extern void             cm_hb_shadow_mx_Destroy       (CM_HB_SHADOW_MX *mx);
extern int              cm_hb_shadow_mx_SizeNeeded    (CM_t *cm, char *errbuf, CP9Bands_t *cp9b, int64_t *ret_ny_cells, int64_t *ret_nk_cells, float *ret_Mb);

extern CM_TR_HB_SHADOW_MX *cm_tr_hb_shadow_mx_Create  (CM_t *cm);
extern int              cm_tr_hb_shadow_mx_GrowTo     (CM_t *cm, CM_TR_HB_SHADOW_MX *mx, char *errbuf, CP9Bands_t *cp9b, int L, float size_limit);
extern int              cm_tr_hb_shadow_mx_Dump       (FILE *ofp, CM_t *cm, CM_TR_HB_SHADOW_MX *mx, char opt_mode);
extern void             cm_tr_hb_shadow_mx_Destroy    (CM_TR_HB_SHADOW_MX *mx);
extern int              cm_tr_hb_shadow_mx_SizeNeeded (CM_t *cm, char *errbuf, CP9Bands_t *cp9b, int64_t *ret_Jny_cells, int64_t *ret_Lny_cells, int64_t *ret_Rny_cells, 
						       int64_t *ret_Jnk_cells, int64_t *ret_Lnk_cells, int64_t *ret_Rnk_cells, int64_t *ret_Tnk_cells, float *ret_Mb);

extern CM_EMIT_MX      *cm_emit_mx_Create     (CM_t *cm);
extern int              cm_emit_mx_GrowTo     (CM_t *cm, CM_EMIT_MX *mx, char *errbuf, int L, float size_limit);
extern int              cm_emit_mx_Dump       (FILE *ofp, CM_t *cm, CM_EMIT_MX *mx);
extern void             cm_emit_mx_Destroy    (CM_EMIT_MX *mx);
extern int              cm_emit_mx_SizeNeeded (CM_t *cm, char *errbuf, int L, int64_t *ret_l_ncells, int64_t *ret_r_ncells, float *ret_Mb);

extern CM_TR_EMIT_MX   *cm_tr_emit_mx_Create     (CM_t *cm);
extern int              cm_tr_emit_mx_GrowTo     (CM_t *cm, CM_TR_EMIT_MX *mx, char *errbuf, int L, float size_limit);
extern int              cm_tr_emit_mx_Dump       (FILE *ofp, CM_t *cm, CM_TR_EMIT_MX *mx, char opt_mode);
extern void             cm_tr_emit_mx_Destroy    (CM_TR_EMIT_MX *mx);
extern int              cm_tr_emit_mx_SizeNeeded (CM_t *cm, char *errbuf, int L, int64_t *ret_l_ncells, int64_t *ret_r_ncells, float *ret_Mb);

extern CM_HB_EMIT_MX   *cm_hb_emit_mx_Create     (CM_t *cm);
extern int              cm_hb_emit_mx_GrowTo     (CM_t *cm, CM_HB_EMIT_MX *mx, char *errbuf, CP9Bands_t *cp9b, int L, float size_limit);
extern int              cm_hb_emit_mx_Dump       (FILE *ofp, CM_t *cm, CM_HB_EMIT_MX *mx);
extern void             cm_hb_emit_mx_Destroy    (CM_HB_EMIT_MX *mx);
extern int              cm_hb_emit_mx_SizeNeeded (CM_t *cm, char *errbuf, CP9Bands_t *cp9b, int L, int64_t *ret_l_ncells, int64_t *ret_r_ncells, float *ret_Mb);

extern CM_TR_HB_EMIT_MX *cm_tr_hb_emit_mx_Create     (CM_t *cm);
extern int               cm_tr_hb_emit_mx_GrowTo     (CM_t *cm, CM_TR_HB_EMIT_MX *mx, char *errbuf, CP9Bands_t *cp9b, int L, float size_limit);
extern int               cm_tr_hb_emit_mx_Dump       (FILE *ofp, CM_t *cm, CM_TR_HB_EMIT_MX *mx, char opt_mode);
extern void              cm_tr_hb_emit_mx_Destroy    (CM_TR_HB_EMIT_MX *mx);
extern int               cm_tr_hb_emit_mx_SizeNeeded (CM_t *cm, char *errbuf, CP9Bands_t *cp9b, int L, int64_t *ret_l_ncells, int64_t *ret_r_ncells, float *ret_Mb);

extern ScanMatrix_t    *cm_CreateScanMatrix           (CM_t *cm, int W, int *dmin, int *dmax, double beta_W, double beta_qdb, int do_banded, int do_float, int do_int);
extern int              cm_CreateScanMatrixForCM      (CM_t *cm, int do_float, int do_int);           
extern int              cm_FloatizeScanMatrix         (CM_t *cm, ScanMatrix_t *smx);
extern int              cm_IntizeScanMatrix           (CM_t *cm, ScanMatrix_t *smx);
extern int              cm_UpdateScanMatrixForCM      (CM_t *cm);
extern int              cm_FreeFloatsFromScanMatrix   (CM_t *cm, ScanMatrix_t *smx);
extern int              cm_FreeIntsFromScanMatrix     (CM_t *cm, ScanMatrix_t *smx);
extern void             cm_FreeScanMatrix             (CM_t *cm, ScanMatrix_t *smx);
extern void             cm_FreeScanMatrixForCM        (CM_t *cm);
extern void             cm_DumpScanMatrixAlpha        (CM_t *cm, int j, int i0, int doing_float);

extern TrScanMatrix_t  *cm_CreateTrScanMatrix         (CM_t *cm, int W, int *dmax, double beta_W, double beta_qdb, int do_banded, int do_float, int do_int);
extern int              cm_FloatizeTrScanMatrix       (CM_t *cm, TrScanMatrix_t *trsmx);
extern int              cm_IntizeTrScanMatrix         (CM_t *cm, TrScanMatrix_t *trsmx);
extern int              cm_FreeFloatsFromTrScanMatrix (CM_t *cm, TrScanMatrix_t *trsmx);
extern int              cm_FreeIntsFromTrScanMatrix   (CM_t *cm, TrScanMatrix_t *trsmx);
extern void             cm_FreeTrScanMatrix           (CM_t *cm, TrScanMatrix_t *trsmx);
extern void             cm_DumpTrScanMatrixAlpha      (CM_t *cm, TrScanMatrix_t *trsmx, int j, int i0, int doing_float);

extern TrScanInfo_t *   CreateTrScanInfo();

extern GammaHitMx_t    *CreateGammaHitMx              (int L, int64_t i0, float cutoff);
extern void             FreeGammaHitMx                (GammaHitMx_t *gamma);
extern int              UpdateGammaHitMx              (CM_t *cm, char *errbuf, GammaHitMx_t *gamma, int j, int dmin, int dmax, float *bestsc, int *bestr, char *bestmode, TrScanInfo_t *trsi, int W, double **act);
extern int              ReportHitsGreedily            (CM_t *cm, char *errbuf, int j, int dmin, int dmax, float *bestsc, int *bestr, char *bestmode, TrScanInfo_t *trsi, int W, double **act, int64_t i0, int64_t j0, float cutoff, CM_TOPHITS *hitlist);
extern void             TBackGammaHitMx               (GammaHitMx_t *gamma, CM_TOPHITS *hitlist, int64_t i0, int64_t j0);



/* from cm_parsetree.c */
extern Parsetree_t *CreateParsetree(int size);
extern void         GrowParsetree(Parsetree_t *tr);
extern void         FreeParsetree(Parsetree_t *tr);
extern float        SizeofParsetree(Parsetree_t *tr);
extern int          InsertTraceNode(Parsetree_t *tr, int y, int whichway, int emitl, int emitr, int state);
extern int          InsertTraceNodewithMode(Parsetree_t *tr, int y, int whichway, int emitl, int emitr, int state, char mode);
extern void         PrintParsetree(FILE *fp, Parsetree_t *tr);
extern void         ParsetreeCount(CM_t *cm, Parsetree_t *tr, ESL_DSQ *dsq, float wgt);
extern int          ParsetreeScore(CM_t *cm, CMEmitMap_t *emap, char *errbuf, Parsetree_t *tr, ESL_DSQ *dsq, int do_null2, float *ret_sc, float *ret_struct_sc, float *ret_primary_sc, int *ret_spos, int *ret_epos);
extern void         ParsetreeDump(FILE *fp, Parsetree_t *tr, CM_t *cm, ESL_DSQ *dsq, int *dmin, int *dmax);
extern int          ParsetreeCompare(Parsetree_t *t1, Parsetree_t *t2);
extern void         SummarizeMasterTrace(FILE *fp, Parsetree_t *tr);
extern void         MasterTraceDisplay(FILE *fp, Parsetree_t *mtr, CM_t *cm);
extern int          Parsetrees2Alignment(CM_t *cm, char *errbuf, const ESL_ALPHABET *abc, ESL_SQ **sq, float *wgt, Parsetree_t **tr, 
					 char **postcode, int nseq, FILE *insertfp, FILE *elfp, 
					 int do_full, int do_matchonly, int be_efficient, ESL_MSA **ret_msa);
extern int          Alignment2Parsetrees(ESL_MSA *msa, CM_t *cm, Parsetree_t *mtr, char *errbuf, ESL_SQ ***ret_sq, Parsetree_t ***ret_tr);
extern float        ParsetreeScore_Global2Local(CM_t *cm, Parsetree_t *tr, ESL_DSQ *dsq, int print_flag);
extern int          Parsetree2CP9trace(CM_t *cm, Parsetree_t *tr, CP9trace_t **ret_cp9_tr);
extern void         rightjustify(const ESL_ALPHABET *abc, char *s, int n);
extern void         leftjustify(const ESL_ALPHABET *abc, char *s, int n);
extern int          EmitParsetree(CM_t *cm, char *errbuf, ESL_RANDOMNESS *r, char *name, int do_digital, Parsetree_t **ret_tr, ESL_SQ **ret_sq, int *ret_N);
extern int          ParsetreeScoreCorrectionNull2(CM_t *cm, char *errbuf, Parsetree_t *tr, ESL_DSQ *dsq, int start, float omega, float *ret_sc);
extern int          ParsetreeScoreCorrectionNull3(CM_t *cm, char *errbuf, Parsetree_t *tr, ESL_DSQ *dsq, int start, float omega, float *ret_sc);
extern int          ParsetreeCountMPEmissions(CM_t *cm, Parsetree_t *tr);
extern void         ScoreCorrectionNull3(const ESL_ALPHABET *abc, float *null0, float *comp, int len, float omega, float *ret_sc);
extern void         ScoreCorrectionNull3CompUnknown(const ESL_ALPHABET *abc, float *null0, ESL_DSQ *dsq, int start, int stop, float omega, float *ret_sc);
extern char         ParsetreeMode(Parsetree_t *tr);

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
extern int cp9_Viterbi(CM_t *cm, char *errbuf, CP9_MX *mx, ESL_DSQ *dsq, int i0, int j0, int W, int hit_len_guess, float cutoff, 
		       int do_scan, int doing_align, int be_efficient, int do_null3, int **ret_psc, int *ret_maxres, 
		       CP9trace_t **ret_tr, float *ret_sc);
extern int cp9_ViterbiBackward(CM_t *cm, char *errbuf, CP9_MX *mx, ESL_DSQ *dsq, int i0, int j0, int W, int hit_len_guess, float cutoff,
			       int do_scan, int doing_align, int j_is_fixed, int be_efficient, int do_null3, int **ret_psc, int *ret_maxres, 
			       CP9trace_t **ret_tr, float *ret_sc);
extern int cp9_Forward(CM_t *cm, char *errbuf, CP9_MX *mx, ESL_DSQ *dsq, int i0, int j0, int W, int hit_len_guess, float cutoff,
		       int do_scan, int doing_align, int be_efficient, int do_null3, int **ret_psc, int *ret_maxres, float *ret_sc);
extern int cp9_FastForward(CM_t *cm, char *errbuf, CP9_MX *mx, ESL_DSQ *dsq, int i0, int j0, int W, int hit_len_guess, float cutoff,
			   int do_scan, int doing_align, int be_efficient, int be_safe, int do_null3, int **ret_psc, int *ret_maxres, float *ret_sc);
extern int cp9_Backward(CM_t *cm, char *errbuf, CP9_MX *mx, ESL_DSQ *dsq, int i0, int j0, int W, int hit_len_guess, float cutoff,
			int do_scan, int doing_align, int j_is_fixed, int be_efficient, int do_null3, int **ret_psc, int *ret_maxres, 
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
extern int   CP9TraceScoreCorrectionNull2(CP9_t *hmm, char *errbuf, CP9trace_t *tr, ESL_DSQ *dsq, int start, float omega, float *ret_sc);

/* from alphabet.c */
extern void   PairCount(const ESL_ALPHABET *abc, float *counters, ESL_DSQ syml, ESL_DSQ symr, float wt);
extern float  DegeneratePairScore(const ESL_ALPHABET *abc, float *esc, ESL_DSQ syml, ESL_DSQ symr);
extern int    iDegeneratePairScore(const ESL_ALPHABET *abc, int *esc, ESL_DSQ syml, ESL_DSQ symr);
extern char   resolve_degenerate (ESL_RANDOMNESS *r, char c);
extern int    revcomp(const ESL_ALPHABET *abc, ESL_SQ *comp, ESL_SQ *sq);
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
extern Fancyali_t    *CreateFancyAli(const ESL_ALPHABET *abc, Parsetree_t *tr, CM_t *cm, CMConsensus_t *cons, ESL_DSQ *dsq, int do_noncanonical, char *pcode);
extern void           PrintFancyAli(FILE *fp, Fancyali_t *ali, int64_t offset, int in_revcomp, int do_top, int linewidth);
extern void           FreeFancyAli(Fancyali_t *ali);
extern int            CreateCMConsensus(CM_t *cm, const ESL_ALPHABET *abc, float pthresh, float sthresh, CMConsensus_t **ret_cons);
extern void           FreeCMConsensus(CMConsensus_t *con);
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
extern int    cm_EntropyWeight(CM_t *cm, const Prior_t *pri, double etarget, double min_Neff, int pretend_cm_is_hmm, double *ret_hmm_re, double *ret_Neff);
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
extern int          cp9_HMM2ijBands(CM_t *cm, char *errbuf, CP9Bands_t *cp9b, CP9Map_t *cp9map, int i0, int j0, int doing_search, int do_trunc, int debug_level);
extern int          cp9_HMM2ijBands_OLD(CM_t *cm, char *errbuf, CP9Bands_t *cp9b, CP9Map_t *cp9map, int i0, int j0, int doing_search, int debug_level);
extern CP9Bands_t  *AllocCP9Bands(int cm_M, int hmm_M);
extern void         FreeCP9Bands(CP9Bands_t *cp9bands);
extern int          cp9_Seq2Bands     (CM_t *cm, char *errbuf, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *pmx, ESL_DSQ *dsq, int i0, int j0, CP9Bands_t *cp9b, int doing_search, TrScanInfo_t *trsi, int debug_level);
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
			       int **hdmin, int **hdmax, int do_trunc, int debug_level);
extern void         combine_qdb_hmm_d_bands(CM_t *cm, int *jmin, int *jmax, int **hdmin, int **hdmax);
extern void         hd2safe_hd_bands(int M, int *jmin, int *jmax, int **hdmin, int **hdmax,
				     int *safe_hdmin, int *safe_hdmax);
extern void         debug_print_hd_bands(CM_t *cm, int **hdmin, int **hdmax, int *jmin, int *jmax);
extern void         PrintDPCellsSaved_jd(CM_t *cm, int *jmin, int *jmax, int **hdmin, int **hdmax, int W);
extern void         debug_print_ij_bands(CM_t *cm);
extern void         debug_print_parsetree_and_ij_bands(FILE *fp, Parsetree_t *tr, CM_t *cm, ESL_DSQ *dsq, CP9Bands_t *cp9b);
extern void         cp9_ShiftCMBands(CM_t *cm, int i, int j, int do_trunc);
extern CP9Bands_t  *cp9_CloneBands(CP9Bands_t *src_cp9b, char *errbuf);
extern void         cp9_PredictStartAndEndPositions(CP9_MX *pmx, CP9Bands_t *cp9b, int i0, int j0);
extern void         cp9_MarginalCandidatesFromStartEndPositions(CM_t *cm, CP9Bands_t *cp9b, TrScanInfo_t *trsi);

/* from cm_p7_modelconfig_trunc.c */
extern int p7_ProfileConfig5PrimeTrunc(P7_PROFILE *gm, int L);
extern int p7_ProfileConfig3PrimeTrunc(const P7_HMM *hmm, P7_PROFILE *gm, int L);
extern int p7_ProfileConfig5PrimeAnd3PrimeTrunc(P7_PROFILE *gm, int L);
extern int p7_ReconfigLength5PrimeTrunc(P7_PROFILE *gm, int L);
extern int p7_ReconfigLength3PrimeTrunc(P7_PROFILE *gm, int L);

/* from cm_p7_modelmaker.c */
#if 0
extern int          BuildP7HMM_MatchEmitsOnly(CM_t *cm, P7_HMM **ret_p7);
#endif
extern int          BuildP7HMM_MatchEmitsOnly(CM_t *cm, P7_HMM **ret_p7);
extern int          cm_cp9_to_p7(CM_t *cm);
extern int          cm_p7_Calibrate(P7_HMM *hmm, char *errbuf, int ElmL, int ElvL, int ElfL, int EgfL, int ElmN, int ElvN, int ElfN, int EgfN, double ElfT, double EgfT, double *ret_gfmu, double *ret_gflambda);
extern int          cm_p7_Tau(ESL_RANDOMNESS *r, char *errbuf, P7_OPROFILE *om, P7_PROFILE *gm, P7_BG *bg, int L, int N, double lambda, double tailp, double *ret_tau);
extern int          cm_SetFilterHMM(CM_t *cm, P7_HMM *hmm, double gfmu, double gflambda);
extern int          dump_p7(P7_HMM *hmm, FILE *fp);

/* from cm_p7_band.c */
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
extern int          p7_Seq2Bands(CM_t *cm, char *errbuf, P7_PROFILE *gm, P7_GMX *gx, P7_BG *bg, P7_TRACE *p7_tr, ESL_DSQ *dsq, int L, 
				 double **phi, float sc7, int len7, int end7, float mprob7, float mcprob7, float iprob7, float ilprob7, int pad7,
				 int **ret_i2k, int **ret_kmin, int **ret_kmax, int *ret_ncells);

extern int          CP9NodeForPosnP7B(CP9_t *hmm, char *errbuf, int x, CP9_MX *post, int kn, int kx, int *ret_node, int *ret_type, int print_flag);
extern int          P7BandsAdjustForSubCM(int *kmin, int *kmax, int L, int spos, int epos);

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
extern int cm_justread_MPIUnpack(ESL_ALPHABET **abc, char *buf, int n, int *pos, MPI_Comm comm, CM_t **ret_cm);
extern int cm_justread_MPIPack(CM_t *cm, char *buf, int n, int *pos, MPI_Comm comm);
extern int cm_justread_MPIPackSize(CM_t *cm, MPI_Comm comm, int *ret_n);
extern int cm_dsq_MPISend(ESL_DSQ *dsq, int L, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
extern int cm_dsq_MPIRecv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, ESL_DSQ **ret_dsq, int *ret_L);
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
extern int exp_info_MPIPackSize(ExpInfo_t *exp, MPI_Comm comm, int *ret_n);
extern int exp_info_MPIPack(ExpInfo_t *exp, char *buf, int n, int *position, MPI_Comm comm);
extern int exp_info_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, ExpInfo_t **ret_exp);
extern int cmcalibrate_exp_results_MPIPackSize(float *scA, int nseq, MPI_Comm comm, int *ret_n);
extern int cmcalibrate_exp_results_MPIPack(float *scA, int nseq, char *buf, int n, int *position, MPI_Comm comm);
extern int cmcalibrate_exp_results_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, float **ret_scA, int *ret_nseq);
extern int cmcalibrate_filter_results_MPIPackSize(int nseq, MPI_Comm comm, int *ret_n);
extern int cmcalibrate_filter_results_MPIPack(float *cyk_scA, float *ins_scA, float *fwd_scA, int nseq, int seq_offset, char *buf, int n, int *position, MPI_Comm comm);
extern int cmcalibrate_filter_results_MPIUnpack(char *buf, int n, int *pos, MPI_Comm comm, float **ret_cyk_scA, float **ret_ins_scA, float **ret_fwd_scA, int *ret_nseq, int *ret_seq_offset);
extern int comlog_MPIPackSize(ComLog_t *comlog, MPI_Comm comm, int *ret_n);
extern int comlog_MPIPack    (ComLog_t *comlog, char *buf, int n, int *position, MPI_Comm comm);
extern int comlog_MPIUnpack  (char *buf, int n, int *pos, MPI_Comm comm, ComLog_t **ret_comlog);
extern int cm_pipeline_MPISend(CM_PIPELINE *pli, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
extern int cm_pipeline_MPIRecv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, ESL_GETOPTS *go, CM_PIPELINE **ret_pli);
extern int cm_tophits_MPISend(CM_TOPHITS *th, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
extern int cm_tophits_MPIRecv(int source, int tag, MPI_Comm comm, char **buf, int *nalloc, CM_TOPHITS **ret_th);
#endif

/* from prior.c */
extern Prior_t *Prior_Create(void);
extern void     Prior_Destroy(Prior_t *pri);
extern Prior_t *Prior_Read(FILE *fp);
extern void     PriorifyCM(CM_t *cm, const Prior_t *pri);
extern Prior_t *Prior_Default(void);
extern Prior_t *Prior_Default_v0p56_through_v1p02(void);

/* from rnamat.c */
extern int numbered_nucleotide (char c);
extern int numbered_basepair (char c, char d);
extern FILE *MatFileOpen (char *matfile);
extern fullmat_t *ReadMatrix(const ESL_ALPHABET *abc, FILE *matfp);
extern int ribosum_MSA_resolve_degeneracies(fullmat_t *fullmat, ESL_MSA *msa);
extern int ribosum_calc_targets(fullmat_t *fullmat);
extern void FreeMat(fullmat_t *fullmat);

/* from seqstoaln.c */
extern seqs_to_aln_t *CreateSeqsToAln(int size, int i_am_mpi_master);
extern seqs_to_aln_t *CreateSeqsToAlnFromSq(ESL_SQ **sq, int size, int i_am_mpi_master);
extern int            GrowSeqsToAln(seqs_to_aln_t *seqs_to_aln, int new_alloc, int i_am_mpi_master); 
extern void           FreeSeqsToAln(seqs_to_aln_t *seqs_to_aln);
extern void           FreePartialSeqsToAln(seqs_to_aln_t *s, int do_free_sq, int do_free_tr, int do_free_cp9_tr, int do_free_post, int do_free_sc, int do_free_pp, int do_free_struct_sc);
extern int            ReadSeqsToAln(const ESL_ALPHABET *abc, ESL_SQFILE *seqfp, int nseq, seqs_to_aln_t *seqs_to_aln, int i_am_mpi_master); 
extern seqs_to_aln_t *CMEmitSeqsToAln(ESL_RANDOMNESS *r, CM_t *cm, int ncm, int nseq, int padW, double *pdist, int i_am_mpi_master);
extern seqs_to_aln_t *RandomEmitSeqsToAln(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, double *pdist, int extranum, int nseq, double *L_distro, int Lmax, int i_am_mpi_master); 

/* from stats.c */
extern int        debug_print_expinfo_array(CM_t *cm, char *errbuf, ExpInfo_t **expA);
extern int        debug_print_expinfo(ExpInfo_t *exp);
extern int        get_gc_comp(const ESL_ALPHABET *abc, ESL_DSQ *dsq, int start, int stop);
extern int        get_alphabet_comp(const ESL_ALPHABET *abc, ESL_DSQ *dsq, int start, int stop, float **ret_freq); 
extern int        GetDBSize (ESL_SQFILE *sqfp, char *errbuf, long *ret_N, int *ret_nseq, int *ret_namewidth);
extern int        GetDBInfo(const ESL_ALPHABET *abc, ESL_SQFILE *sqfp, char *errbuf, long *ret_N, int *ret_nseq, double **ret_gc_ct);
extern int        E2ScoreGivenExpInfo(ExpInfo_t *exp, char *errbuf, float E, float *ret_sc);
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
extern int        CreateGenomicHMM(const ESL_ALPHABET *abc, char *errbuf, double **ret_sA, double ***ret_tAA, double ***ret_eAA, int *ret_nstates);
extern int        SampleGenomicSequenceFromHMM(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, char *errbuf, double *sA, double **tAA, double **eAA, int nstates, int L, ESL_DSQ **ret_dsq);
extern int        CopyExpInfo(ExpInfo_t *src, ExpInfo_t *dest);

/* from truncyk.c */
void  SetMarginalScores_reproduce_bug_i27(CM_t *cm);
float LeftMarginalScore_reproduce_bug_i27(const ESL_ALPHABET *abc, float *esc, ESL_DSQ dres);
float RightMarginalScore_reproduce_bug_i27(const ESL_ALPHABET *abc, float *esc, ESL_DSQ dres);

float TrCYK_DnC(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, int j0, Parsetree_t **ret_tr);
float TrCYK_Inside(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, int j0, int lenCORREX, Parsetree_t **ret_tr);
/* legacy, avoid use: */
float trinside (CM_t *cm, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0, int do_full,
                void ****ret_shadow, void ****ret_L_shadow, void ****ret_R_shadow,
                void ****ret_T_shadow, void ****ret_Lmode_shadow, void ****ret_Rmode_shadow,
                int *ret_mode, int *ret_v, int *ret_i, int *ret_j);


/* from cm_pipeline.c */
extern CM_PIPELINE *cm_pipeline_Create (ESL_GETOPTS *go, ESL_ALPHABET *abc, int clen_hint, int L_hint, int64_t Z, enum cm_zsetby_e Z_setby, enum cm_pipemodes_e mode);
extern int          cm_pipeline_Reuse  (CM_PIPELINE *pli);
extern void         cm_pipeline_Destroy(CM_PIPELINE *pli, CM_t *cm);
extern int          cm_pipeline_Merge  (CM_PIPELINE *p1, CM_PIPELINE *p2);

extern int cm_pli_TargetReportable  (CM_PIPELINE *pli, float score,     double Eval);
extern int cm_pli_TargetIncludable  (CM_PIPELINE *pli, float score,     double Eval);
extern int cm_pli_NewModel          (CM_PIPELINE *pli, int modmode, CM_t *cm, int cm_clen, int cm_W, int need_fsmx, int need_smx, int *fcyk_dmin, int *fcyk_dmax, int *final_dmin, int *final_dmax, P7_OPROFILE *om, P7_BG *bg, int64_t cur_cm_idx);
extern int cm_pli_NewModelThresholds(CM_PIPELINE *pli, CM_t *cm);
extern int cm_pli_NewSeq            (CM_PIPELINE *pli, const ESL_SQ *sq, int64_t cur_seq_idx);
extern int cm_pli_p7Filter          (CM_PIPELINE *pli, P7_OPROFILE *om, P7_BG *bg, float *p7_evparam, const ESL_SQ *sq, int64_t **ret_ws, int64_t **ret_we, int *ret_nwin);
extern int cm_pli_p7EnvelopeDef     (CM_PIPELINE *pli, P7_OPROFILE *om, P7_BG *bg, float *p7_evparam, const ESL_SQ *sq, int64_t *ws, int64_t *we, int nwin, P7_PROFILE **opt_gm, P7_PROFILE **opt_Rgm, P7_PROFILE **opt_Lgm, P7_PROFILE **opt_Tgm, int64_t **ret_es, int64_t **ret_ee, int *ret_nenv);
extern int cm_pli_CMStage           (CM_PIPELINE *pli, off_t cm_offset, int cm_config_opts, const ESL_SQ *sq, int64_t *es, int64_t *ee, int nenv, CM_TOPHITS *hitlist, CM_t **opt_cm, CMConsensus_t **opt_cmcons);
extern int cm_pli_AlignHit          (CM_PIPELINE *pli, CM_t *cm, CMConsensus_t *cmcons, const ESL_SQ *sq, int do_trunc, CM_HIT *hit, int first_hit, CP9Bands_t *scan_cp9b);
extern int cm_Pipeline              (CM_PIPELINE *pli, off_t cm_offset, int cm_config_opts, P7_OPROFILE *om, P7_BG *bg, float *p7_evparam, const ESL_SQ *sq, CM_TOPHITS *hitlist, P7_PROFILE **opt_gm, P7_PROFILE **opt_Rgm, P7_PROFILE **opt_Lgm, P7_PROFILE **opt_Tgm, CM_t **opt_cm, CMConsensus_t **opt_cmcons);
extern int cm_pli_Statistics(FILE *ofp, CM_PIPELINE *pli, ESL_STOPWATCH *w);

/* from cm_p7_domaindef.c */
extern int p7_domaindef_GlocalByPosteriorHeuristics(const ESL_SQ *sq, P7_PROFILE *gm, P7_GMX *gxf, P7_GMX *gxb,
						    P7_GMX *fwd, P7_GMX *bck, P7_DOMAINDEF *ddef, int do_null2);
/* from cm_tophits.c */
extern CM_TOPHITS *cm_tophits_Create(void);
extern int         cm_tophits_Grow(CM_TOPHITS *h);
extern int         cm_tophits_CreateNextHit(CM_TOPHITS *h, CM_HIT **ret_hit);
extern int         cm_tophits_SortByScore(CM_TOPHITS *h);
extern int         cm_tophits_SortByPosition(CM_TOPHITS *h);
extern int         cm_tophits_Merge(CM_TOPHITS *h1, CM_TOPHITS *h2);
extern int         cm_tophits_GetMaxPositionLength(CM_TOPHITS *h);
extern int         cm_tophits_GetMaxNameLength(CM_TOPHITS *h);
extern int         cm_tophits_GetMaxAccessionLength(CM_TOPHITS *h);
extern int         cm_tophits_GetMaxShownLength(CM_TOPHITS *h);
extern int         cm_tophits_Reuse(CM_TOPHITS *h);
extern void        cm_tophits_Destroy(CM_TOPHITS *h);
extern int         cm_tophits_CloneHitMostly(CM_TOPHITS *src_th, int h, CM_TOPHITS *dest_th);
extern int         cm_tophits_ComputeEvalues(CM_TOPHITS *th, double eZ, int istart);
extern int         cm_tophits_RemoveDuplicates(CM_TOPHITS *th);
extern int         cm_tophits_UpdateHitPositions(CM_TOPHITS *th, int hit_start, int64_t seq_start, int in_revcomp);

extern int cm_tophits_Threshold(CM_TOPHITS *th, CM_PIPELINE *pli);
extern int cm_tophits_Targets(FILE *ofp, CM_TOPHITS *th, CM_PIPELINE *pli, int textw);
extern int cm_tophits_HitAlignments(FILE *ofp, CM_TOPHITS *th, CM_PIPELINE *pli, int textw);
extern int cm_tophits_HitAlignmentStatistics(FILE *ofp, CM_TOPHITS *th, int used_cyk);
extern int cm_tophits_TabularTargets(FILE *ofp, char *qname, char *qacc, CM_TOPHITS *th, CM_PIPELINE *pli, int show_header);
extern int cm_tophits_TabularTail(FILE *ofp, const char *progname, enum cm_pipemodes_e pipemode, 
				  const char *qfile, const char *tfile, const ESL_GETOPTS *go);

extern int cm_tophits_Dump(FILE *fp, const CM_TOPHITS *th);
extern int cm_hit_Dump    (FILE *fp, const CM_HIT *h);

/* cm_alidisplay.c */
extern CM_ALIDISPLAY *cm_alidisplay_Create(const ESL_ALPHABET *abc, Parsetree_t *tr, CM_t *cm, CMConsensus_t *cons, const ESL_SQ *sq, 
					   int64_t seqoffset, char *pcode, float aln_sc, int used_optacc, int used_hbands, float matrix_Mb, double elapsed_secs);
extern CM_ALIDISPLAY *cm_alidisplay_Clone(const CM_ALIDISPLAY *ad);
extern size_t         cm_alidisplay_Sizeof(const CM_ALIDISPLAY *ad);
extern int            cm_alidisplay_Serialize(CM_ALIDISPLAY *ad);
extern int            cm_alidisplay_Deserialize(CM_ALIDISPLAY *ad);
extern void           cm_alidisplay_Destroy(CM_ALIDISPLAY *ad);
extern char           cm_alidisplay_EncodePostProb(float p);
extern float          cm_alidisplay_DecodePostProb(char pc);
extern int            cm_alidisplay_Print(FILE *fp, CM_ALIDISPLAY *ad, int min_aliwidth, int linewidth, int show_accessions, int do_noncanonicals);
extern int            cm_alidisplay_Dump(FILE *fp, const CM_ALIDISPLAY *ad);
extern int            cm_alidisplay_Compare(const CM_ALIDISPLAY *ad1, const CM_ALIDISPLAY *ad2);


#include "esl_config.h"
#include "config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_sqio.h"
#include "esl_msa.h"

#include "structs.h"

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#define USE_NEWLOGSUM 1
#define USE_OLDLOGSUM 0

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
				      double ***ret_gamma, float **ret_seqlen);
extern int      BandTruncationNegligible(double *density, int b, int W, double *ret_beta);
extern int      BandMonteCarlo(CM_t *cm, int nsample, int W, double ***ret_gamma);
extern void     FreeBandDensities(CM_t *cm, double **gamma);


extern void     BandBounds(double **gamma, int M, int W, double p, 
			   int **ret_min, int **ret_max);
extern void     PrintBandGraph(FILE *fp, double **gamma, int *min, int *max, int v, int W);

extern void     PrintDPCellsSaved(CM_t *cm, int *min, int *max, int W);
extern float    CYKBandedScan(CM_t *cm, ESL_DSQ *dsq, int *dmin, int *dmax, int i0, int j0, int W, 
			      float cutoff, search_results_t *results);
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
extern Fancyali_t    *CreateFancyAli(const ESL_ALPHABET *abc, Parsetree_t *tr, CM_t *cm, CMConsensus_t *cons, ESL_DSQ *dsq, char *pcode1, char *pcode2);
extern void           PrintFancyAli(FILE *fp, Fancyali_t *ali, int offset, int in_revcomp);
extern void           FreeFancyAli(Fancyali_t *ali);
extern int            CreateCMConsensus(CM_t *cm, const ESL_ALPHABET *abc, float pthresh, float sthresh, CMConsensus_t **ret_cons);
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
extern int   ConfigCM(CM_t *cm, int *preset_dmin, int *preset_dmax);
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
				     CP9_MX *post, int L);
extern void  ConfigLocalEnforce(CM_t *cm, float p_internal_start, float p_internal_exit);
extern int   EnforceSubsequence(CM_t *cm);
extern float EnforceScore(CM_t *cm);
extern int   EnforceFindEnfStart(CM_t *cm, int enf_cc_start);
extern int   ConfigForGumbelMode(CM_t *cm, int statmode);
extern int   ConfigQDB(CM_t *cm);
extern void  CMHackInsertScores(CM_t *cm);
extern void  CP9HackInsertScores(CP9_t *cp9);

/* from modelmaker.c
 */
extern void HandModelmaker(ESL_MSA *msa, int use_rf, float gapthresh, 
			   CM_t **ret_cm, Parsetree_t **ret_mtr);
extern void ConsensusModelmaker(const ESL_ALPHABET *abc, char *ss_cons, int clen, CM_t **ret_cm, 
				Parsetree_t **ret_gtr);
extern Parsetree_t *Transmogrify(CM_t *cm, Parsetree_t *gtr, 
				 ESL_DSQ *dsq, char *aseq, int alen);
extern void cm_from_guide(CM_t *cm, Parsetree_t *gtr);
extern int  cm_find_and_detach_dual_inserts(CM_t *cm, int do_check, int do_detach);
extern int  cm_check_before_detaching(CM_t *cm, int insert1, int insert2);
extern int  cm_detach_state(CM_t *cm, int insert1, int insert2);
extern int  clean_cs(char *cs, int alen, int be_quiet);
				 
/* from parsetree.c
 */
extern Parsetree_t *CreateParsetree(int size);
extern void         GrowParsetree(Parsetree_t *tr);
extern void         FreeParsetree(Parsetree_t *tr);
extern int          InsertTraceNode(Parsetree_t *tr, int y, int whichway, 
				    int emitl, int emitr, int state);
extern int          InsertTraceNodewithMode(Parsetree_t *tr, int y, int whichway, 
				    int emitl, int emitr, int state, int mode);
extern void         ParsetreeCount(CM_t *cm, Parsetree_t *tr, ESL_DSQ *dsq, float wgt);
extern float        ParsetreeScore(CM_t *cm, Parsetree_t *tr, ESL_DSQ *dsq, int do_null2);
extern void         PrintParsetree(FILE *fp, Parsetree_t *tr);
extern void         ParsetreeDump(FILE *fp, Parsetree_t *tr, CM_t *cm, ESL_DSQ *dsq, int *dmin, int *dmax);
extern int          ParsetreeCompare(Parsetree_t *t1, Parsetree_t *t2);
extern void         SummarizeMasterTrace(FILE *fp, Parsetree_t *tr);
extern void         MasterTraceDisplay(FILE *fp, Parsetree_t *mtr, CM_t *cm);
extern int          Parsetrees2Alignment(CM_t *cm, const ESL_ALPHABET *abc, ESL_SQ **sq, float *wgt, 
					 Parsetree_t **tr, int nseq, int do_full, int do_matchonly, ESL_MSA **ret_msa);
extern float        ParsetreeScore_Global2Local(CM_t *cm, Parsetree_t *tr, ESL_DSQ *dsq, int print_flag);
extern int          Parsetree2CP9trace(CM_t *cm, Parsetree_t *tr, CP9trace_t **ret_cp9_tr);
extern void         rightjustify(const ESL_ALPHABET *abc, char *s, int n);
extern void         leftjustify(const ESL_ALPHABET *abc, char *s, int n);

/* from smallcyk.c
 */
extern float CYKDivideAndConquer(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, int j0, 
				 Parsetree_t **ret_tr, int *dmin, int *dmax);
extern float CYKInside(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, int j0, 
		       Parsetree_t **ret_tr, int *dmin, int *dmax);
extern float CYKInsideScore(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, 
			    int j0, int *dmin, int *dmax);
extern float CYKDemands(CM_t *cm, int L, int *dmin, int *dmax, int be_quiet);
extern void  debug_print_bands(FILE *fp, CM_t *cm, int *dmin, int *dmax);

/* size calculators - not normally part of external API, but truncyk.c currently uses them */
extern float insideT_size(CM_t *cm, int L, int r, int z, int i0, int j0);
extern float vinsideT_size(CM_t *cm, int r, int z, int i0, int i1, int j1, int j0);

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
float TrCYK_DnC(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, int j0, Parsetree_t **ret_tr);
float TrCYK_Inside(CM_t *cm, ESL_DSQ *dsq, int L, int r, int i0, int j0, Parsetree_t **ret_tr);
/* legacy, avoid use: */
float trinside (CM_t *cm, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0, int do_full,
                void ****ret_shadow, void ****ret_L_shadow, void ****ret_R_shadow,
                void ****ret_T_shadow, void ****ret_Lmode_shadow, void ****ret_Rmode_shadow,
                int *ret_mode, int *ret_v, int *ret_i, int *ret_j);
/*EPN, Thu Sep 20 13:40:01 2007, what I had in my 1.0 branch prior to merge: 
  float TrCYKInside(CM_t *cm, ESL_DSQ *sq, int L, int r, int i0, int j0, Parsetree_t **ret_tr, int *dmin, int *dmax);
*/

/* from cm_eweight.c
 * Entropy-based sequence weighting */
extern double CM_Eweight(CM_t *cm,  const Prior_t *pri, 
			 float numb_seqs, float targetent);
extern void ModelContent(float *ent1, float *ent2, int M);
extern void CMRescale(CM_t *hmm, float scale);
extern double CM_Eweight_RE(CM_t *cm, const Prior_t *pri, float numb_seqs, 
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
extern void CP9EnforceHackMatchScores(CP9_t *cp9, int enf_start_pos, int enf_end_pos);
extern void CP9_fake_tracebacks(ESL_MSA *msa, int *matassign, CP9trace_t ***ret_tr);
extern void  CP9TraceCount(CP9_t *hmm, ESL_DSQ *dsq, float wt, CP9trace_t *tr);
extern float CP9TraceScore(CP9_t *hmm, ESL_DSQ *dsq, CP9trace_t *tr);
extern void  CP9PrintTrace(FILE *fp, CP9trace_t *tr, CP9_t *hmm, ESL_DSQ *dsq);
extern char *CP9Statetype(char st);
extern int   CP9TransitionScoreLookup(CP9_t *hmm, char st1, int k1, 
				    char st2, int k2);
extern void  CP9ViterbiTrace(CP9_t *hmm, ESL_DSQ *dsq, int i0, int j0,
			     CP9_MX *mx, CP9trace_t **ret_tr);
extern void  CP9ReverseTrace(CP9trace_t *tr);
extern int   CP9Traces2Alignment(CM_t *cm, const ESL_ALPHABET *abc, ESL_SQ **sq, float *wgt, 
				 int nseq, CP9trace_t **tr, int do_full, int do_matchonly, ESL_MSA **ret_msa);
extern void  DuplicateCP9(CM_t *src_cm, CM_t *dest_cm);
/* more from cplan9.c: functions stolen from HMMER-2.4::mathsupport.c */
extern int   Prob2Score(float p, float null);
extern float Score2Prob(int sc, float null);
extern float Scorify(int sc);
extern int   DegenerateSymbolScore(float *p, float *null, int ambig);

/* from cm_fastsearch.c */
extern int  FastCYKScan      (CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq, int i0, int j0, float cutoff, search_results_t *results, float **ret_vsc, float *ret_sc);
extern int  FastIInsideScan  (CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq, int i0, int j0, float cutoff, search_results_t *results, float **ret_vsc, float *ret_sc);
extern int  XFastIInsideScan (CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq, int i0, int j0, float cutoff, search_results_t *results, float **ret_vsc, float *ret_sc);
extern int  X2FastIInsideScan(CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq, int i0, int j0, float cutoff, search_results_t *results, float **ret_vsc, float *ret_sc);
extern int  FastFInsideScan  (CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq, int i0, int j0, float cutoff, search_results_t *results, float **ret_vsc, float *ret_sc);
extern int  RefCYKScan       (CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq, int i0, int j0, float cutoff, search_results_t *results, float **ret_vsc, float *ret_sc);
extern int  RefIInsideScan   (CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq, int i0, int j0, float cutoff, search_results_t *results, float **ret_vsc, float *ret_sc);
extern int  XRefIInsideScan  (CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq, int i0, int j0, float cutoff, search_results_t *results, float **ret_vsc, float *ret_sc);
extern int  RefFInsideScan   (CM_t *cm, char *errbuf, ScanMatrix_t *smx, ESL_DSQ *dsq, int i0, int j0, float cutoff, search_results_t *results, float **ret_vsc, float *ret_sc);
extern int  rsearch_CYKScan  (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float cutoff, int D, search_results_t *results, float *ret_sc);
extern int  FastCYKScanHB    (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, float cutoff, search_results_t *results, CM_HB_MX *mx, float *ret_sc);
extern int  FastFInsideScanHB(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, float cutoff, search_results_t *results, CM_HB_MX *mx, float *ret_sc);
extern int  cm_CountSearchDPCalcs(CM_t *cm, char *errbuf, int L, int *dmin, int *dmax, int W, float **ret_vcalcs, float *ret_calcs);

/* from CP9_scan.c */
extern float CP9Viterbi(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, int **ret_sc, 
			int *ret_bestpos, search_results_t *results, int do_scan, int doing_align, 
			int be_efficient, CP9_MX **ret_mx, CP9trace_t **ret_tr);
extern float CP9Forward(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, int **ret_isc, 
			int *ret_maxres, search_results_t *results, int do_scan, int doing_align, 
			int be_efficient, CP9_MX **ret_mx);
extern float CP9Backward(CM_t *cm, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, int **ret_isc, 
			 int *ret_maxres, search_results_t *results, int do_scan, int doing_align, 
			 int be_efficient, CP9_MX **ret_mx);
extern float CP9ForwardScanDemands(CP9_t *cp9, int L);
/*
  extern float CP9Scan_dispatch(CM_t *cm,ESL_DSQ *dsq, int i0, int j0, int W, float cm_cutoff, 
  float cp9_cutoff, search_results_t *results, int doing_cp9_stats, int *ret_flen);
extern float RescanFilterSurvivors(CM_t *cm, ESL_DSQ *dsq, search_results_t *hmm_results, int i0, 
				   int j0, int W, int padmode, int ipad, int jpad, int do_collapse,
				   float cm_cutoff, float cp9_cutoff, search_results_t *results, 
				   int *ret_flen);
extern float FindCP9FilterThreshold(CM_t *cm, CMStats_t *cmstats, ESL_RANDOMNESS *r, 
				    float Fmin, float Smin, float Starget, float Spad, int N, 
				    int use_cm_cutoff, float cm_ecutoff, int db_size, 
				    int emit_mode, int fthr_mode, int hmm_gum_mode, 
				    int do_fastfil, int do_Fstep, int my_rank, int nproc, 
				    int do_mpi, char *histfile, FILE *Rpts_fp, float *ret_F);
extern float FindSubFilterThreshold(CM_t *cm, CMStats_t *cmstats, ESL_RANDOMNESS *r, 
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
*/     
/* from CP9_cm2wrhmm.c */
extern int build_cp9_hmm(CM_t *cm, CP9_t **ret_hmm, CP9Map_t **ret_cp9map, int do_psi_test,
			 float psi_vs_phi_threshold, int debug_level);
extern void CP9_map_cm2hmm(CM_t *cm, CP9Map_t *cp9map, int debug_level);
extern void map_helper(CM_t *cm, CP9Map_t *cp9map, int k, int ks, int v);
extern int  CP9_check_wrhmm(CM_t *cm, CP9_t *hmm, CP9Map_t *cp9map, int debug_level);
extern void fill_psi(CM_t *cm, double *psi, char ***tmap);
extern void fill_phi_cp9(CP9_t *hmm, double ***ret_phi, int spos);
extern void make_tmap(char ****ret_tmap);

extern int  CP9_check_by_sampling(CM_t *cm, CP9_t *hmm, ESL_RANDOMNESS *r, CMSubInfo_t *subinfo, int spos, int epos, 
				  float chi_thresh, int nsamples, int print_flag);
extern void debug_print_cp9_params(FILE *fp, CP9_t *hmm, int print_scores);
extern void debug_print_phi_cp9(CP9_t *hmm, double **phi);
extern CP9Map_t *AllocCP9Map(CM_t *cm);
extern void FreeCP9Map(CP9Map_t *cp9map);
extern int  MakeDealignedString(const ESL_ALPHABET *abc, char *aseq, int alen, char *ss, char **ret_s);




/* from cm_masks.c */
extern float CM_TraceScoreCorrection(CM_t *cm, Parsetree_t *tr, ESL_DSQ *dsq);

/* from sub_cm.c */
extern int  build_sub_cm(CM_t *orig_cm, CM_t **ret_cm, int sstruct, int estruct, CMSubMap_t **ret_submap, int print_flag);
extern void CP9NodeForPosn(CP9_t *hmm, int i0, int j0, int x, CP9_MX *post, int *ret_node, int *ret_type, float pmass, int is_start, int print_flag);
extern void StripWUSSGivenCC(ESL_MSA *msa, float gapthresh, int first_match, int last_match);
extern int  check_orig_psi_vs_sub_psi(CM_t *orig_cm, CM_t *sub_cm, CMSubMap_t *submap, double threshold, 
				       int print_flag);
extern int  check_sub_cm(CM_t *orig_cm, CM_t *sub_cm, CMSubMap_t *submap, CMSubInfo_t *subinfo, float pthresh, int print_flag);
extern int check_sub_cm_by_sampling(CM_t *orig_cm, CM_t *sub_cm, ESL_RANDOMNESS *r, CMSubMap_t *submap, CMSubInfo_t *subinfo,
				    float chi_thresh, int nsamples, int print_flag);
extern int  check_sub_cm_by_sampling2(CM_t *orig_cm, CM_t *sub_cm, ESL_RANDOMNESS *r, int spos, int epos, int nseq);
extern int  sub_cm2cm_parsetree(CM_t *orig_cm, CM_t *sub_cm, Parsetree_t **ret_orig_tr, Parsetree_t *sub_tr, 
				CMSubMap_t *submap, int print_flag);
extern CMSubMap_t  *AllocSubMap(CM_t *sub_cm, CM_t *orig_cm, int sstruct, int estruct);
extern void         FreeSubMap(CMSubMap_t *submap);
extern CMSubInfo_t *AllocSubInfo(int clen);
extern void         FreeSubInfo(CMSubInfo_t *subinfo);
extern void  debug_print_cm_params(FILE *fp, CM_t *cm);

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

/* from cm_errors.c */
extern void cm_Die (char *format, ...);
extern void cm_Fail(char *format, ...);

/* from rnamat.c */
extern int numbered_nucleotide (char c);
extern int numbered_basepair (char c, char d);
extern FILE *MatFileOpen (char *matfile);
extern fullmat_t *ReadMatrix(const ESL_ALPHABET *abc, FILE *matfp);
extern int ribosum_calc_targets(fullmat_t *fullmat);

/* from cm_dispatch.c */
extern int ActuallySearchTarget(CM_t *cm, char *errbuf, int fround, ESL_DSQ *dsq, int i0, int j0, 
				search_results_t **results, int *ret_flen, float *ret_sc);
extern int ActuallyAlignTargets(CM_t *cm, char *errbuf, seqs_to_aln_t *seqs_to_aln, ESL_DSQ *dsq, search_results_t *results, 
				int first_result, int bdump_level, int debug_level, int silent_mode, ESL_RANDOMNESS *r);
extern int  revcomp(const ESL_ALPHABET *abc, ESL_SQ *comp, ESL_SQ *sq);
extern int  read_next_search_seq(const ESL_ALPHABET *abc, ESL_SQFILE *seqfp, int do_revcomp, dbseq_t **ret_dbseq);
extern void print_results (CM_t *cm, SearchInfo_t *si, const ESL_ALPHABET *abc, CMConsensus_t *cons, dbseq_t *dbseq, int do_complement);
extern void remove_hits_over_e_cutoff (CM_t *cm, SearchInfo_t *si, search_results_t *results, ESL_SQ *sq);
extern seqs_to_aln_t *CreateSeqsToAln(int size, int i_am_mpi_master);
extern seqs_to_aln_t *CreateSeqsToAlnFromSq(ESL_SQ **sq, int size, int i_am_mpi_master);
extern int GrowSeqsToAln(seqs_to_aln_t *seqs_to_aln, int new_alloc, int i_am_mpi_master); 
extern void FreeSeqsToAln(seqs_to_aln_t *seqs_to_aln);
extern void FreePartialSeqsToAln(seqs_to_aln_t *s, int do_free_sq, int do_free_tr, int do_free_cp9_tr, int do_free_post, int do_free_sc);
extern int  ReadSeqsToAln(const ESL_ALPHABET *abc, ESL_SQFILE *seqfp, int nseq, int do_read_all, seqs_to_aln_t *seqs_to_aln, int i_am_mpi_master); 
extern seqs_to_aln_t *CMEmitSeqsToAln(ESL_RANDOMNESS *r, CM_t *cm, int ncm, int nseq, int i_am_mpi_master);
extern seqs_to_aln_t *RandomEmitSeqsToAln(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, double *pdist, int extranum, int nseq, int L, int i_am_mpi_master); 

/* from logsum.c: (stolen from HMMER3 dev code) EPN, Fri Sep  7 16:56:45 2007 */

/* from prior.c */
extern Prior_t *Prior_Create(void);
extern void     Prior_Destroy(Prior_t *pri);
extern Prior_t *Prior_Read(FILE *fp);
extern void     PriorifyCM(CM_t *cm, const Prior_t *pri);
extern Prior_t *Prior_Default(void);
extern struct p7prior_s *P7DefaultInfernalPrior(void);

/* from stats.c */
extern CMStats_t *AllocCMStats(int np);
extern void FreeCMStats(CMStats_t *cmstats);
extern int debug_print_cmstats(CMStats_t *cmstats, int has_fthr);
extern int debug_print_gumbelinfo(GumbelInfo_t *evd);
extern int debug_print_filterthrinfo(CMStats_t *cmstats, CP9FilterThr_t *fthr);

extern int  get_gc_comp(ESL_SQ *sq, int start, int stop);
extern void OLD_serial_make_histogram (int *gc_count, int *partitions, int num_partitions,
				       CM_t *cm, int num_samples, 
				       int sample_length, int doing_cp9_stats,
				       int use_easel);
extern void GetDBInfo(const ESL_ALPHABET *abc, ESL_SQFILE *sqfp, long *ret_N, double **ret_gc_ct);

extern float e_to_score (float E, double *mu, double *lambda);

extern double RJK_ExtremeValueE (float x, double mu, double lambda);

extern char resolve_degenerate (ESL_RANDOMNESS *r, char c);

extern float MinScCutoff (CM_t *cm, SearchInfo_t *si, int n);
extern int   CM2Gumbel_mode(CM_t *cm, int search_opts, int *ret_cm_gum_mode, int *ret_cp9_gum_mode);
extern int   CopyFThrInfo(CP9FilterThr_t *src, CP9FilterThr_t *dest);
extern int   CopyCMStatsGumbel(CMStats_t *src, CMStats_t *dest);
extern int   CopyCMStats(CMStats_t *src, CMStats_t *dest);
/*extern int SetCMCutoff(CM_t *cm, int cm_cutoff_type, float cm_sc_cutoff, float cm_e_cutoff);
  extern int SetCP9Cutoff(CM_t *cm, int cp9_cutoff_type, float cp9_sc_cutoff, float cp9_e_cutoff,
  float cm_e_cutoff);
extern int PrintSearchInfo(FILE *fp, CM_t *cm, int cm_mode, int cp9_mode, long N);
extern float MinCP9ScCutoff (CM_t *cm);
extern void parallel_make_histogram (int *gc_count, int *partitions, int num_partitions, 
			      CM_t *cm, int num_samples, int sample_length,
			      int doing_cp9_stats,
			      int mpi_my_rank, int mpi_num_procs, 
			      int mpi_master_rank);
*/

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



extern void mpi_worker_search_target(CM_t *cm, int my_rank);
extern void mpi_worker_cm_and_cp9_search(CM_t *cm, int do_fast, int my_rank);
extern void mpi_worker_cm_and_cp9_search_maxsc(CM_t *cm, int do_fast, int do_minmax, int my_rank);
extern int dsq_MPISend(ESL_DSQ *dsq, int L, int dest);
extern int dsq_MPIRecv(ESL_DSQ **ret_dsq, int *ret_L);
extern int dsq_maxsc_MPISend(char *dsq, int L, float maxsc, int dest);
extern int dsq_maxsc_MPIRecv(char **ret_dsq, int *ret_L, float *ret_maxsc);


#endif

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

/* from logsum.c */
extern void  init_ilogsum(void);
extern int   ILogsum(int s1, int s2);
extern int   ILogsumNI(int s1, int s2);
extern int   ILogsumNI_diff(int s1a, int s1b, int s2a, int s2b, int db);

extern void  FLogsumInit(void);
extern float LogSum2(float p1, float p2);
extern float FLogsum(float p1, float p2);

/* from cp9_fastsearch.c */
extern int cp9_FastViterbi(CM_t *cm, char *errbuf, CP9_MX *mx, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, search_results_t *results, 
			   int do_scan, int doing_align, int be_efficient, int **ret_psc, int *ret_maxres, 
			   CP9trace_t **ret_tr, float *ret_sc);
extern int cp9_FastViterbiBackward(CM_t *cm, char *errbuf, CP9_MX *mx, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, search_results_t *results, 
				   int do_scan, int doing_align, int be_efficient, int **ret_psc, int *ret_maxres, 
				   CP9trace_t **ret_tr, float *ret_sc);
extern int cp9_FastForward(CM_t *cm, char *errbuf, CP9_MX *mx, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, search_results_t *results, 
			   int do_scan, int doing_align, int be_efficient, int **ret_psc, int *ret_maxres, float *ret_sc);
extern int Xcp9_FastForward(CM_t *cm, char *errbuf, CP9_MX *mx, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, search_results_t *results, 
			    int do_scan, int doing_align, int be_efficient, int be_safe, int **ret_psc, int *ret_maxres, float *ret_sc);
extern int Xcp9_FastBackward(CM_t *cm, char *errbuf, CP9_MX *mx, ESL_DSQ *dsq, int i0, int j0, int W, float cutoff, search_results_t *results, 
			     int do_scan, int doing_align, int be_efficient, int **ret_psc, int *ret_maxres, 
			     float *ret_sc);
extern int cp9_WorstForward(CM_t *cm, char *errbuf, CP9_MX *mx, int thresh, int doing_scan, int doing_align, int *ret_L);
extern int cp9_CheckTransitionGuarantees(CP9_t *cp9, char *errbuf);
extern int cp9_GetLocalityMode(CP9_t *cp9, char *errbuf, int *ret_mode);

/* from cm_fastalign.c */
extern int fast_cyk_align_hb (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0, void ****ret_shadow, 
				int allow_begin, int *ret_b, float *ret_bsc, CM_HB_MX *mx, float *ret_sc);
extern int fast_cyk_align    (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, int vroot, int vend, int i0, int j0, void ****ret_shadow,  
			        int allow_begin, int *ret_b, float *ret_bsc, float ****ret_mx, float *ret_sc);
extern int optimal_accuracy_align_hb(CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, int i0, int j0, void ****ret_shadow,  
				     int *ret_b, float *ret_bsc, CM_HB_MX *mx, CM_HB_MX *post_mx, float *ret_pp);
extern int optimal_accuracy_align   (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, int i0, int j0, void ****ret_shadow,  
				     int *ret_b, float *ret_bsc, float ****ret_mx, float ***post_mx, float *ret_pp);
extern int fast_alignT_hb    (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, Parsetree_t *tr, int r, int z, int i0, int j0, 
			        int allow_begin, CM_HB_MX *mx, int do_optacc, CM_HB_MX *post_mx, float *ret_sc);
extern int fast_alignT       (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, Parsetree_t *tr, int r, int z, int i0, int j0, 
			        int allow_begin, float ****ret_mx, int do_optacc, float ***post_mx, float *ret_sc);
extern int FastAlignHB        (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, int i0, int j0, CM_HB_MX *mx,     int do_optacc, CM_HB_MX *post_mx, Parsetree_t **ret_tr, char **ret_pcode1, char **ret_pcode2, float *ret_sc);
extern int FastAlign          (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, int i0, int j0, float ****ret_mx, int do_optacc, float ****ret_post_mx, Parsetree_t **ret_tr, char **ret_pcode1, char **ret_pcode2, float *ret_sc);
extern int FastInsideAlignHB  (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, CM_HB_MX *mx,     float *ret_sc);
extern int FastInsideAlign    (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, float ****ret_mx, float *ret_sc);
extern int FastOutsideAlignHB (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, CM_HB_MX *mx,    CM_HB_MX *ins_mx, int do_check, float *ret_sc);
extern int FastOutsideAlign   (CM_t *cm, char *errbuf, ESL_DSQ *dsq, int i0, int j0, float ****ret_mx, float ***ins_mx, int do_check, float *ret_sc);
extern int SampleFromInsideHB (ESL_RANDOMNESS *r, CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, CM_HB_MX *mx, Parsetree_t **ret_tr, float *ret_sc);
extern int SampleFromInside   (ESL_RANDOMNESS *r, CM_t *cm, char *errbuf, ESL_DSQ *dsq, int L, float ***mx,  Parsetree_t **ret_tr, float *ret_sc);
extern int   ** alloc_jdbanded_vjd_kshadow_deck(int L, int i, int j, int jmin, int jmax, int *hdmin, int *hdmax);
extern char  ** alloc_jdbanded_vjd_yshadow_deck(int L, int i, int j, int jmin, int jmax, int *hdmin, int *hdmax);

/* TO BE cm_postprob.c */
extern void CMPostalCode(CM_t *cm, int L, float ***post, Parsetree_t *tr, char **ret_pcode1, char **ret_pcode2);
extern void CMPostalCodeHB(CM_t *cm, int L, CM_HB_MX *post_mx, Parsetree_t *tr, char **ret_pcode1, char **ret_pcode2);
extern float FScore2Prob(float sc, float null);
extern int   Fscore2postcode(float sc);
extern int CMPosteriorHB      (CM_t *cm, char *errbuf, int i0, int j0, CM_HB_MX *ins_mx, CM_HB_MX *out_mx, CM_HB_MX *post_mx);
extern int CMPosterior        (CM_t *cm, char *errbuf, int i0, int j0, float ***ins_mx, float ***out_mx, float ***post_mx);
extern int CMCheckPosteriorHB (CM_t *cm, char *errbuf, int i0, int j0, CM_HB_MX *post);
extern int CMCheckPosterior   (CM_t *cm, char *errbuf, int i0, int j0, float ***post);


/* from hmmband.c */
extern CP9Bands_t * AllocCP9Bands(CM_t *cm, CP9_t *hmm);
extern void         FreeCP9Bands(CP9Bands_t *cp9bands);
extern double       DScore2Prob(int sc, float null);
extern int          cp9_Seq2Bands     (CM_t *cm, char *errbuf, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *pmx, ESL_DSQ *dsq, int i0, int j0, CP9Bands_t *cp9b, int doing_search, int debug_level);
extern int          cp9_Seq2Posteriors(CM_t *cm, char *errbuf, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *pmx, ESL_DSQ *dsq, int i0, int j0, int debug_level);
extern int          cp9_FB2HMMBands        (CP9_t *hmm, char *errbuf, ESL_DSQ *dsq, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *pmx, CP9Bands_t *cp9b, 
				            int i0, int j0, int M, double p_thresh, int did_scan, int debug_level);
extern int          cp9_FB2HMMBandsWithSums(CP9_t *hmm, char *errbuf, ESL_DSQ *dsq, CP9_MX *fmx, CP9_MX *bmx, CP9_MX *pmx, CP9Bands_t *cp9b, 
					    int i0, int j0, int M, double p_thresh, int did_scan, int debug_level);
extern void         cp9_Posterior(ESL_DSQ *dsq, int i0, int j0, CP9_t *hmm, CP9_MX *fmx, CP9_MX *bmx,
				  CP9_MX *mx, int did_scan);
extern void         cp9_IFillPostSums(CP9_MX *post, CP9Bands_t *cp9, int i0, int j0);
extern int          cp9_HMM2ijBands(CM_t *cm, char *errbuf, CP9Bands_t *cp9b, CP9Map_t *cp9map, int i0, int j0, int debug_level);
extern void         cp9_RelaxRootBandsForSearch(CM_t *cm, int *imin, int *imax, int *jmin, int *jmax);
extern void         cp9_DebugPrintHMMBands(FILE *ofp, int L, CP9Bands_t *cp9b, double hmm_bandp, int debug_level);
extern int          cp9_CheckFB(CP9_MX *fmx, CP9_MX *bmx, CP9_t *hmm, char *errbuf, float sc, int i0, int j0, ESL_DSQ *dsq);
extern void         cp9_CompareBands(CP9Bands_t *cp9b1, CP9Bands_t *cp9b2);
extern int          cp9_GrowHDBands(CP9Bands_t *cp9b, char *errbuf);
extern void         ijBandedTraceInfoDump(CM_t *cm, Parsetree_t *tr, int *imin, int *imax, 
					  int *jmin, int *jmax, int debug_level);
extern void         ijdBandedTraceInfoDump(CM_t *cm, Parsetree_t *tr, int *imin, int *imax, 
					   int *jmin, int *jmax, int **hdmin, int **hdmax, 
					   int debug_level);
extern int          cp9_ValidateBands(CM_t *cm, char *errbuf, CP9Bands_t *cp9b, int i0, int j0);
extern void ij2d_bands(CM_t *cm, int L, int *imin, int *imax, int *jmin, int *jmax,
		       int **hdmin, int **hdmax, int debug_level);
extern void combine_qdb_hmm_d_bands(CM_t *cm, int *jmin, int *jmax, int **hdmin, int **hdmax);
extern void hd2safe_hd_bands(int M, int *jmin, int *jmax, int **hdmin, int **hdmax,
			     int *safe_hdmin, int *safe_hdmax);
extern void debug_print_hd_bands(CM_t *cm, int **hdmin, int **hdmax, int *jmin, int *jmax);
extern void PrintDPCellsSaved_jd(CM_t *cm, int *jmin, int *jmax, int **hdmin, int **hdmax, int W);

/* old functions (get rid of them ?) */
extern float CP9ViterbiAlign (ESL_DSQ *dsq, int i0, int j0, CP9_t *hmm, CP9_MX *mx, CP9trace_t **ret_tr);
extern float CP9ForwardAlign (ESL_DSQ *dsq, int i0, int j0, CP9_t *hmm, CP9_MX *mx);
extern float CP9BackwardAlign(ESL_DSQ *dsq, int i0, int j0, CP9_t *hmm, CP9_MX *mx);

/* from cm_mx.c */
extern CM_HB_MX *       cm_hb_mx_Create            (int M);
extern int              cm_hb_mx_GrowTo            (CM_t *cm, CM_HB_MX *mx, char *errbuf, CP9Bands_t *cp9b, int L);
extern int              cm_hb_mx_Dump              (FILE *ofp, CM_HB_MX *mx);
extern void             cm_hb_mx_Destroy           (CM_HB_MX *mx);
extern ScanMatrix_t *   cm_CreateScanMatrix        (CM_t *cm, int W, int *dmin, int *dmax, double beta, int do_banded, int do_float, int do_int);
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
extern void             FreeOptimizedEmitScores    (float **fesc_vAA, int **iesc_vAA, int M);
extern float **         FCalcInitDPScores          (CM_t *cm);
extern int **           ICalcInitDPScores          (CM_t *cm);
extern GammaHitMx_t    *CreateGammaHitMx           (int L, int i0, int be_greedy, float cutoff, int do_backward);
extern void             FreeGammaHitMx             (GammaHitMx_t *gamma);
extern void             UpdateGammaHitMxCM         (GammaHitMx_t *gamma, int j, float *alpha_row, int dn, int dx, int using_hmm_bands, int *bestr, int doing_inside, search_results_t *results);
extern void             UpdateGammaHitMxCP9Forward (GammaHitMx_t *gamma, int i, int j, float hit_sc, search_results_t *results);
extern void             UpdateGammaHitMxCP9Backward(GammaHitMx_t *gamma, int i, int j, float hit_sc, search_results_t *results);
extern void             TBackGammaHitMxForward     (GammaHitMx_t *gamma, search_results_t *results, int i0, int j0);
extern void             TBackGammaHitMxBackward    (GammaHitMx_t *gamma, search_results_t *results, int i0, int j0);

/* from cp9_mx.c */
extern CP9_MX *CreateCP9Matrix(int N, int M);
extern void    FreeCP9Matrix  (CP9_MX *mx);
extern int     GrowCP9Matrix  (CP9_MX *mx, char *errbuf, int N, int M, int ***mmx, int ***imx, int ***dmx, int ***elmx, int **erow);

/* from cm_cp9_hybridsearch.c */
extern int cm_cp9_HybridScan(CM_t *cm, char *errbuf, CP9_MX *mx, ESL_DSQ *dsq, HybridScanInfo_t *hsi, int i0, int j0, int W, 
			     float cutoff, search_results_t *results, int **ret_psc, int *ret_maxres, float *ret_sc);
extern int predict_xsub(CM_t *cm, float *cm_vcalcs, float *cm_expsc, float *cp9_expsc);
extern void cm_CalcAvgHitLength(CM_t *cm, double beta, float **ret_hitlen);
extern HybridScanInfo_t * cm_CreateHybridScanInfo(CM_t *cm, double hsi_beta, float full_cm_ncalcs);
extern int cm_AddRootToHybridScanInfo(CM_t *cm, HybridScanInfo_t *hsi, int vroot_to_add);
extern int cm_ValidateHybridScanInfo(CM_t *cm, HybridScanInfo_t *hsi);
extern void cm_FreeHybridScanInfo(HybridScanInfo_t *hsi, CM_t *cm);

/* from cm_searchinfo.c */
extern int  cm_CreateSearchInfo(CM_t *cm, int cutoff_type, float cutoff);
extern int  cm_AddFilterToSearchInfo(CM_t *cm, int cyk_filter, int inside_filter, int viterbi_filter, int forward_filter,
				     int hybrid_filter, ScanMatrix_t *smx, HybridScanInfo_t *hsi, int cutoff_type, float cutoff);
extern void cm_FreeSearchInfo(SearchInfo_t *si, CM_t *cm);
extern void cm_DumpSearchInfo(SearchInfo_t *si);
extern void DumpSearchOpts(int search_opts);
extern void cm_ValidateSearchInfo(CM_t *cm, SearchInfo_t *fi);
extern void cm_UpdateSearchInfoCutoff(CM_t *cm, int nround, int cutoff_type, float cutoff);
extern search_results_t *CreateResults (int size);
extern void ExpandResults (search_results_t *r, int additional);
extern void AppendResults (search_results_t *src_results, search_results_t *dest_results, int i0);
extern void FreeResults   (search_results_t *r);
extern int  compare_results (const void *a_void, const void *b_void);
extern void sort_results (search_results_t *results);
extern void report_hit (int i, int j, int bestr, float score, search_results_t *results);
extern void remove_overlapping_hits (search_results_t *results, int i0, int j0);
extern float CountScanDPCalcs(CM_t *cm, int L, int use_qdb);


/* from cm_theta.c */
extern int cm_CalcMaxSc(CM_t *cm, double **ret_maxsc, double **ret_maxsc_noss);
extern Theta_t *cm_CalcTheta(CM_t *cm, Theta_t **ret_theta, float stepsize);

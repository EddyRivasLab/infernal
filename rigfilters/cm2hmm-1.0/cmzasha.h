/*
This file copyright (c) 2004, Zasha Weinberg
All rights reserved.

Redistribution and use in source and binary forms, 
with or without modification, are permitted 
provided that the following conditions are met:

- Redistributions of source code must retain the 
above copyright notice, this list of conditions 
and the following disclaimer. 
- Redistributions in binary form must reproduce 
the above copyright notice, this list of 
conditions and the following disclaimer in the 
documentation and/or other materials provided 
with the distribution. 
- Neither the name of the University of Washington 
nor the names of its contributors may be used to 
endorse or promote products derived from this 
software without specific prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS 
AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS 
OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN 
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN 
IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
// generic #includes for stuff in my project


// some other defines in cmzasha.cpp,ScancykBranchAndBound_Heuristic.cpp
/* these were for branch&bound, which doesn't work well.  Okay, at all.
#define ENABLE_PRUNING
#define USE_CELL_VISITED
#define MAX_CELL_VISITED 2
*/

// forward decls
class TransitionCounter;
class BoolBasePairMatrix;
class HmmWithPenalties;
class LeftNucEventList;
class HmmOptimizedForRememberNucScanning;
class SequenceSet;
// end forward decls

extern bool dumpHmmCommitteePruningStats;
extern const char *overrideHmmCacheFileName; // NULL means don't override
extern const char *dumpHmmScoresFileName; // NULL means don't dump it.
extern const char *dumpFracLetsThruByScoreThresholdFileName; // NULL means don't dump it
extern FILE *dumpHmmScoresFile;
extern const char *dumpLastHeurScoreInSequenceFileName; // NULL means don't dump
extern FILE *dumpLastHeurScoreInSequenceFile;
extern bool useThresholdOnCm;
extern bool cmpWithInfernal;
extern bool cmpWithPureHmm;
extern char *nucs;
extern bool disableHmmBuildInfoDump; // for testing
extern bool analScoreDumping; // I swear, it's not what it sounds like
extern bool dumpAverageScore;
extern FILE *dumpFilterResultsFile;
extern bool enableProgressFile;

#ifndef DISABLE_ZRAND
extern zrand::ZRandom *GetRander (void);
extern int RandInt (int N);
extern double rand01 (void);
extern void InitRand (void);
extern void DestroyRand (void);
extern void ReinitRandWithFixedKnuth (void);
#endif

extern void AnotherParam (void);
extern bool HasAnotherParam (int a,int argc);
extern void AnotherParam (int a,int argc);
extern std::string DumpProgramParams(int argc,char **argv,bool doPreamble);

inline double log2 (double x) {
	return log(x)/log(2.0);
}
inline double pow2 (double x) {
	return pow(2.0,x);
}

extern std::string LineBreaksToSpaces (std::string inStr);
extern std::string LineBreaksToTabs (std::string inStr);
extern std::string GetFileNameFromFullPath (std::string path);

#include <NaryCounter.h>
#include <stl_extra.h>
#include "MarkovModelStats.h"


extern std::string GetSubsequence (const char *rnaSequence,int first,int last);

struct TopLevelMatch { // matches corresponding to starting at state=0 (start state)
	int windowLast,windowLen;
	float score;

	// redundant
	int windowFirst;

	bool operator < (const TopLevelMatch& t) const {
		// sort first by windowFirst
		if (windowFirst!=t.windowFirst) {
			return windowFirst<t.windowFirst;
		}
		// then by windowLen
		return windowLen<t.windowLen;
	}
};
struct CykscanStats {
	MultiplyArray3d<float> scores;
	vector<float> scoresPerWindowLast;
	std::string programParams;

	vector<float> hmmScoresPerWindowLast; // for things with 2nd struct whose scores shouldn't ever be higher than a pure HMM
	vector<vector<float> > fullHmmDynProgTable; // first dimension is windowLast, next is hmm state

	bool isValid;
	bool collectScores;
	bool collectAboveThreshold;
	bool collectScoresPerWindowLast; // for comparing with HMM, just collects the max score from Start node over all possible windowLen
	bool collectHmmScoresPerWindowLast;
	bool collectHmmFullDynProgTable;

	typedef std::list<TopLevelMatch> TopLevelMatchList;
	TopLevelMatchList topLevelMatchesAboveThreshold;
	double threshold;

	CykscanStats (void) {
		isValid=false;
		collectScores=false;
		collectAboveThreshold=false;
		collectScoresPerWindowLast=false;
		collectHmmScoresPerWindowLast=false;
		collectHmmFullDynProgTable=false;
	}
};

extern void CYKScanZasha(CM_t *cm, char *dsq, int L, int W, 
			 int *ret_nhits, int **ret_hitr, int **ret_hiti, int **ret_hitj, float **ret_hitsc,
			 CykscanStats& cykscanStats);
extern void CYKScan_OptionalBug(CM_t *cm, char *dsq, int L, int W, 
			 int *ret_nhits, int **ret_hitr, int **ret_hiti, int **ret_hitj, float **ret_hitsc,
			 bool fixBug=true);

#include "CovarianceModel.h"
#include "InfernalHmm.h"

	// (taken out of HmmWithPenalties to avoid a cyclic compilation dependency, which was difficult to solve)
	// this structure describes a fact that we can remember about what happened with the left nuc, in terms of a set of possible, mutually exclusive events.  The potential elements of the set: either (1) the left nuc was in the state corresponding to an MP (this requires type 1/Expanded-type/separated-MP HMMs) and was one of A,C,G,U, or (2) the left nuc didn't use the MP state, so was an aberrant case relative to the consensus alignment
	// and now I'm too lazy to change the name of this class even though it also does right-nucs
	class LeftNucEvent {
	protected:
		bool member[MAXABET+1];
	public:
		LeftNucEvent (void) { SetEmpty(); }
		~LeftNucEvent () {}
		bool IsEmpty (void) const { for (int i=0; i<MAXABET+1; i++) { if (member[i]) { return false; } } return true; }
		void SetEmpty (void) { for (int i=0; i<MAXABET+1; i++) { member[i]=false; isLeft=true; } }
		void SetLeftNucInMP (int nuc,bool isInSet) { assert(nuc>=0 && nuc<MAXABET); member[nuc]=isInSet; }
		bool HasLeftNucInMP (int nuc) const { assert(nuc>=0 && nuc<MAXABET); return member[nuc]; }
		void SetLeftNotInMP (bool isInSet) { member[MAXABET]=isInSet; }
		bool HasLeftNotInMP (void) const { return member[MAXABET]; }
		void Dump (FILE *out) const;
		void DumpFilterSpec(FILE *out) const;
		void DumpFilterSpec(char *out) const;
		bool operator == (const LeftNucEvent& t) const;
		bool operator != (const LeftNucEvent& t) const;

		bool isLeft;
	};
	class LeftNucEventList : public std::list<LeftNucEvent> {
	public:
		void Dump(FILE *out) const;
		void DumpFilterSpec(FILE *out) const;
		void DumpFilterSpec(char *out) const;
	};

// I think this is a reasonable upper bound, & I think it'll be worth it in memory usage
#define MAX_CHILDREN_WITH_BIFURICATIONS 512


// the meta-state that is represented in the Alpha dyn prog table
struct AlphaState {
	CovarianceModel::State state;
	int windowLast,windowLen;
};


// code to iterate in a really, really generic way -- too generic for some purposes (like I think confusing for prioritization), but nice for others
struct ChildGenericTransition {

	typedef FixedArrayWithSize<AlphaState,2> ChildrenToSum; // size is 1 for most state types, 2 for bifurication nodes
	ChildrenToSum childrenToSum;

	float tsc; // cost of this transition (tsc=0 for bifurication nodes);
};
typedef FixedArrayWithSize<ChildGenericTransition,MAX_CHILDREN_WITH_BIFURICATIONS> ChildGenericTransitionVector;
// another way of doing things that doesn't require the expensive vector
inline int GetNumChildrenAlphaStates(const CovarianceModel& cm,const AlphaState& alphaState)
{
	if (cm.IsBifurication(alphaState.state)) {
		return alphaState.windowLen+1;
	}
	else {
		return cm.GetNumChildren(alphaState.state);
	}
}
inline void GetNthChildAlphaState (ChildGenericTransition& childTransition,int childNum,const CovarianceModel& cm,const AlphaState& alphaState)
{
	if (cm.IsBifurication(alphaState.state)) {
		assert(childNum>=0 && childNum<GetNumChildrenAlphaStates(cm,alphaState));

		int leftWindowLen=childNum;
		childTransition.childrenToSum.resize(2);
		childTransition.childrenToSum[0].state=cm.GetLeftBifurifactionChild(alphaState.state);
		childTransition.childrenToSum[0].windowLast=alphaState.windowLast - alphaState.windowLen + leftWindowLen;
		childTransition.childrenToSum[0].windowLen=leftWindowLen;
		childTransition.childrenToSum[1].state=cm.GetRightBifurifactionChild(alphaState.state);
		childTransition.childrenToSum[1].windowLast=alphaState.windowLast;
		childTransition.childrenToSum[1].windowLen=alphaState.windowLen - leftWindowLen;
		childTransition.tsc=0;
	}
	else {
		childTransition.childrenToSum.resize(1);
		childTransition.childrenToSum[0]=alphaState;
		childTransition.childrenToSum[0].state=cm.GetNthChildState(alphaState.state,childNum);
		childTransition.tsc=cm.GetNthChildTsc(alphaState.state,childNum);
	}
}
inline void GetChildren (ChildGenericTransitionVector& children,const CovarianceModel& cm,AlphaState alphaState)
{
	int numChildren=GetNumChildrenAlphaStates(cm,alphaState);
	children.resize(numChildren);
	for (int i=0; i<numChildren; i++) {
		GetNthChildAlphaState(children[i],i,cm,alphaState);
	}
}

inline bool IsBaseCase (float& get_score,const CovarianceModel& cm,CovarianceModel::State state,int windowLen) {
	if (windowLen<cm.GetNumSymbolsEmitted(state)) {
		// this is invalid; there's no space left
		get_score=(float)IMPOSSIBLE;
		return true;
	}
	if (cm.IsEndState(state)) {
		if (windowLen==0) {
			// we're done
			get_score=0;
		}
		else {
			// we're not done, but we think we are
			get_score=(float)IMPOSSIBLE;
		}
		return true;
	}
	return false;
}
inline bool IsBaseCase (float& get_score,const CovarianceModel& cm,AlphaState alphaState)
{
	return IsBaseCase(get_score,cm,alphaState.state,alphaState.windowLen);
}

struct SequenceInfo {
	const char *sequenceName;
	int sequenceNum; // this is actually unique for each seq (we don't get seq #0 forward, seq #0 reversed ; rather seq #0 forward, seq #1 reversed.)
	bool isReversed;
};

// forward declarations
class BranchAndBoundHeuristicStuff;
class BranchAndBound_UpperBoundHeuristic;
class BranchAndBound_SearchHeuristic;

#ifndef CM2HMM_ONLY
#include "AlphaDynProgTable.h"
#include "ScancykBranchAndBound.h"
#include "ScancykBranchAndBound_Heuristics.h"
#endif

class HitList : public std::list<std::pair<int,int> > {
public:
	void Init (int length); // one interval from [0,length)
	void Init (int first,int second); // one interval
	void InitEmpty (void);
	void Dump (FILE *file) const;
	int GetOverallLast (void) const;

	int TotalSize (void) const;
	__int64 SizeIn2D (int windowLen) const; // what part of the dynamic programming table must we look at; it's basically TotalSize() * windowLen, except that at the beginning of each interval in the hit list, we only have to worry about a triangular part of the dynamic programming table.
};

extern void GetNucNumsFromHalfOpenInterval(int& startNuc,int& endNuc,int first,int last,int sequenceLen,bool isReversed);
extern void GetNucNumsFromHalfOpenInterval(int& startNuc,int& endNuc,int first,int last,SequenceSet& sequenceSet);

class FracLetsThruCounter {
protected:
	__int64 nucsInAllSeqs,nucsLetsThru;
	__int64 nucsSinceLastProgressReport,nucsLetThruSinceLastProgressReport;
	__int64 size2dOfAllSeqs,size2dLetThru;
	__int64 size2dSinceLastProgressReport,size2dLetThruSinceLastProgressReport;
public:
	FracLetsThruCounter ();
	~FracLetsThruCounter ();

	void DumpFracLetsThru (FILE *out,const char *messagePrefix,bool sinceLastProgressReport);
	__int64 GetNucsInAllSeqs (void);

	void ProcessPruning (const HitList& inputHitList,const HitList& outputHitList,int windowLen);

	double GetFilteringFraction (void) const;
	void ResetCounts (void);
};

#include "Cm2HMM.h"
#ifndef CM2HMM_ONLY
#include "SubCM.h"
#endif
#include "HmmType1.h"
#ifndef CM2HMM_ONLY
#include "ScanHMM.h"
#endif
#include "SymbolicMath.h"
#include "Cm2HmmOptimize.h"
#ifndef CM2HMM_ONLY
#include "Cm2HmmOptimize_Deprecated.h"
#endif

extern int DumpRandomCMWalks (char *cmFileName,int numDumps,double biasingExponent,bool useHackInsertScores);
extern void SearchCMSansMP (int windowLen,float minLodScoreForHit,char *cmFileName,bool doLocalAlignment,char *fastaFileToSearch,bool hackCM,const std::string& programParams);

extern int GetPid (void);
extern const char *GetSqinfoDesc (const SQINFO& sqinfo);
extern void AnnounceFastaFile (FILE *out,const char *fastaFileName,const char *virtualFastaFileName="");
extern void AnnounceCmFile(FILE *out,const char *cmFileName,bool doLocalAlignment);

// assumes that thisWindow is sliding to the right
void AddWindowToList(HitList& hitList,const std::pair<int,int>& thisWindow);

struct RfamSeq {
	std::string name;
	std::string seq;
	char *alloc_digitizedSeq,*digitizedSeq; // place-holder for the client
	bool isFound;
};
typedef std::list<RfamSeq> RfamSeqList;
void ParseRfamMembers(RfamSeqList& rfamSeqList,const char *rfamID,const char *rfamFullCsvFileName);

#include "SequenceSet.h"
#ifndef CM2HMM_ONLY
#include "SearchPruner.h"
#include "SmithWater.h"

// Cm2hmmOptimize.cpp
extern void GeneticHmmSet(int windowLen,float minLodScoreForHit,char *cmFileName,bool doLocalAlignment,char *fastaFileToSearch,const std::string& programParams,char *savedHmmsFileName,int numHmmsToConsider,double adjacentLinkageProb,int numGeneticIters,bool onlyForwardStrand);
extern void TrainHmmEMish(int windowLen,float minLodScoreForHit,char *cmFileName,bool doLocalAlignment,char *fastaFileToSearch,const std::string& programParams,int numIterations,bool onlyForwardStrand);
extern void DumpHmmInflationsPerInequality (char *cmFileName,const std::string& programParams,bool doLocalAlignment);

// Cm2HmmOptimizeCorrel.cpp
extern void GlobalHmmOptimizerCorrel (char *cmFileName,bool doLocalAlignment,const std::string& programParams,bool useProbRatioForEmissions,const char *distanceMeasure);
extern void BuildHmm_MaxLikeliPathCorrespondence(char *cmFileName,bool doLocalAlignment,const char *hmmBinFileName,const std::string& programParams,Cm2Hmm_HmmBuildType hmmType);
extern void BuildHmm_MaxLikeliPathCorrespondence(char *cmFileName,bool doLocalAlignment,InfernalHmm& createdInfernalHmm,const std::string& programParams,Cm2Hmm_HmmBuildType hmmType);

// cmzashaSearch.cpp
extern void SetHeuristics(BranchAndBound_UpperBoundHeuristic *(&upperBoundHeuristic),BranchAndBound_SearchHeuristic *(&searchHeuristic),const char *searchTypeName,const char *BB);
extern void SearchWithPruning(int windowLen,float minLodScoreForHit,const char *cmFileName,bool doLocalAlignment,SequenceSet& sequenceSet,const std::string& programParams,bool runCM,const CovarianceModel& cm,SearchPruner& searchPruner);
extern void SearchWithPruning(double& filteringFraction,double& runTimeInSeconds,int windowLen,float minLodScoreForHit,const char *cmFileName,bool doLocalAlignment,SequenceSet& sequenceSet,const std::string& programParams,bool runCM,const CovarianceModel& cm,SearchPruner& searchPruner,bool enablePrintfs=true);
extern void SearchWithHmmCommittee(const int hmmCommitteeSize,const int windowLen,const float minLodScoreForHit,char *cmFileName,bool doLocalAlignment,SequenceSet& sequenceSet,const std::string& programParams,const bool runCM,const int effectiveCommitteeSize);
extern void SearchWithHmmCommittee_EvaluateEveryone(const int hmmCommitteeSize,const int windowLen,const float minLodScoreForHit,char *cmFileName,bool doLocalAlignment,SequenceSet& sequenceSet,const std::string& programParams,const bool runCM,const char *inputSequenceFileName);
extern void SearchWithHmm(int windowLen,float minLodScoreForHit,char *cmFileName,bool doLocalAlignment,SequenceSet& sequenceSet,const std::string& programParams,bool runCM);
extern void SearchWithHmmWithSubCm(const char *subCmFileName,int windowLen,float minLodScoreForHit,char *cmFileName,bool doLocalAlignment,SequenceSet& sequenceSet,const std::string& programParams,bool runCM);
extern void SearchWithBlastish(int valueOf11,const char *rfamCsvFileName,const char *rfamId,int windowLen,float minLodScoreForHit,char *cmFileName,bool doLocalAlignment,SequenceSet& sequenceSet,const std::string& programParams,bool runCM);
extern void SearchFromBlast(int windowLen,float minLodScoreForHit,char *cmFileName,bool doLocalAlignment,char *fastaFileToSearch,SequenceSet& sequenceSet,char *blastOutputFile,const std::string& programParams,bool extendHitsToWindowLen);
extern void SearchForExactRfamMembers(const char *rfamID,const char *rfamFullCsvFileName,SequenceSet& sequenceSet);
extern void ZashaSearch(const char *searchTypeName,int windowLen,float minLodScoreForHit,char *cmFileName,bool doLocalAlignment,char *fastaFileToSearch);
extern void SearchToLearnMarkov(SequenceSet& sequenceSet,const std::string& programParams,int markovModelOrder,const char *markovSaveFile);
extern void DumpHmmAndCmScores (int windowLen,float minLodScoreForHit,char *cmFileName,bool doLocalAlignment,SequenceSet& sequenceSet,const std::string& programParams);
extern void DumpHmmBlockHeuristic (int windowLen,float minLodScoreForHit,char *cmFileName,bool doLocalAlignment,SequenceSet& sequenceSet,const std::string& programParams,int blockSize,const char *dumpFileName);
extern void CountNucs(char *seqfile);
extern void PartitionSequenceFiles (SequenceSet& sequenceSet,const char *targetDirectory,__int64 maxNucsPerFile);
extern void HashTestSequences (SequenceSet& sequenceSet,const char *targetFileName);

// in FakeCmbuild.cpp
extern void FakeCmbuild (char *alifile,char *cmfile,const char *hmmBinFileName,const char *strategy,const std::string& programParams);

// in CmSimplePathEmit.cpp
extern void CreateRandomStockholmMSA (const char *msaFileName,char *cmFileName,bool doLocalAlignment,int numSeqs);
#endif

// in cmzasha.cpp
extern void MarkovGenerate(int length,int maxSizePerSequence,const char *outputFileName);

extern MarkovModelStats *defaultTrainingMarkovModelStats;
extern void SetNew_defaultTrainingMarkovModelStats (MarkovModelStats *newMarkov);

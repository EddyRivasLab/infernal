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

// classes to do float<-->type conversions for using 32-bit ints instead of floats.  In practice, looks like uses floats works best, and it's certainly easiest.
class FloatScoreHelper {
public:
	typedef float ScoreType;

	inline static ScoreType GetZeroScore (void) {
		return 0;
	}
	inline static ScoreType GetImpossibleScore (void) {
		return (float)IMPOSSIBLE;
	}
	inline static ScoreType GetNegativeInfinityScore (void) {
		return -FLT_MAX;
	}
	inline static ScoreType FloatToScoreType (float t) {
		return t;
	}
	inline static float ScoreTypeToFloat (ScoreType t) {
		return t;
	}
	inline static ScoreType DoubleToScoreType (double t) {
		return (float)t;
	}
};
class DoubleScoreHelper {
public:
	typedef double ScoreType;

	inline static ScoreType GetZeroScore (void) {
		return 0;
	}
	inline static ScoreType GetImpossibleScore (void) {
		return IMPOSSIBLE;
	}
	inline static ScoreType GetNegativeInfinityScore (void) {
		return -DBL_MAX;
	}
	inline static ScoreType FloatToScoreType (float t) {
		return t;
	}
	inline static float ScoreTypeToFloat (ScoreType t) {
		return (float)t;
	}
};
// for 32-bit int scores  --   HACKY for now
class IntScoreHelper {
public:
	typedef int ScoreType;

	inline static int GetBitsOfDecimalAccuracy (void) {
		return 23; // give 23 bits of decimal, leaving 8 bits of integer range (0-255), and 1 sign bit
	}

	inline static ScoreType FloatToScoreType (float t) {
		double scoreAsDouble=ldexp(t,GetBitsOfDecimalAccuracy());
		double rounded=ceil(scoreAsDouble); // make sure we round up, because we want upper bounds
		return (int)rounded;
	}
	inline static float ScoreTypeToFloat (ScoreType t) {
		return (float)(ldexp((double)t,-GetBitsOfDecimalAccuracy()));
	}
	inline static ScoreType GetZeroScore (void) {
		return 0;
	}
	inline static ScoreType GetImpossibleScore (void) {
		return -INT_MAX + (20<<GetBitsOfDecimalAccuracy());
	}
	inline static ScoreType GetNegativeInfinityScore (void) {
		return -INT_MAX;
	}
};

// I re-templated this, because I wanted to be able to use doubles.  So, ignore: I have un-templated HmmType1, because templating seems to be making g++ produce slower code.  Originally, I wanted to be able to use 32-bit ints, but this doesn't give enough of a performance improvement to justify the complexity.  So, I've decided to just back away from templating.
//template <class ScoreHelperType>
class HmmType1 {
public:
    typedef FloatScoreHelper ScoreHelperType; // this is no longer a template param
	typedef ScoreHelperType::ScoreType ScoreType;
	typedef int State;
protected:
	int numStates;
	vector<vector<ScoreType> > singletEmissionScores; // singleEmissionScores[nuc#][state]
	enum {MAX_CHILDREN=5};
	struct StateInfo {
		int numChildren;
		int firstChild;
		ScoreType tsc[MAX_CHILDREN];
		bool isEmitting;
		bool isRightState;
	};
	vector<StateInfo> stateInfoVector;

	// for doLocal
	bool doLocal;
	enum {MAX_ENDSC_LINK_TO_LEFTS=2};
	struct LinkToLeftInfo {
		State hmmLeftState;
		float endsc;
	};
	typedef FixedArrayWithSize<LinkToLeftInfo,MAX_ENDSC_LINK_TO_LEFTS> EndscLinksToLeft;
	struct LocalStateInfo {
		float leftwardBeginsc; // beginsc for HMM left-side states, corresponding to the start of the HMM
		float rightwardBeginsc;

		EndscLinksToLeft endscLinksToLeft;
	};
	vector<LocalStateInfo> localStateInfoVector;

	// probabilities, which are used in the Forward Alg
	struct ForwardAlgInfo {
		float transitionProbs[MAX_CHILDREN];
		float singletEmissionProbs[MAXABET];
	};
	vector<ForwardAlgInfo> forwardAlgInfoVector;
	
	Hmm2CmStateVector hmm2CmStateVector;
	Cm2HmmStateVector cm2HmmStateVector;
public:

	inline HmmType1 (void) {
	}
	inline ~HmmType1 () {
	}

	void Dump (FILE *out) const;
	void Dump (const char *fileName) const;

	void Init (const InfernalHmm& sourceHmm); // construct from the CovarianceModel representation of an HMM

	inline const std::list<CovarianceModel::State>& GetCmStateListOfHmmState (State hmmState) const {
		return hmm2CmStateVector[InfernalHmm::IntToState(hmmState)];
	}
	inline State GetHmmLeftState (CovarianceModel::State cmState) const {
		return InfernalHmm::StateToInt(cm2HmmStateVector[cmState].hmmLeftState);
	}
	inline State GetHmmRightState (CovarianceModel::State cmState) const {
		return InfernalHmm::StateToInt(cm2HmmStateVector[cmState].hmmRightState);
	}

	inline int GetNumStates (void) const {
		return numStates;
	}
	inline State GetFirstState (void) const {
		return 0;
	}
	inline State GetLastState (void) const {
		return GetNumStates();
	}
	inline State GetActualLastState (void) const {
		return GetLastState()-1;
	}
	inline static State GetInvalidState (void) {
		return -1;
	}
	// all states are either left-emitting, or non-emitting.  No other types are necessary for HMMs
	inline bool IsEmittingState (State state) const {
		assert(state>=0 && state<GetNumStates());
		return stateInfoVector[state].isEmitting;
	}
	inline ScoreType GetSingletEmissionScore (State state,int nuc) const {
		assert(state>=0 && state<GetNumStates());
		assert(IsEmittingState(state));
		assert(nuc>=0 && nuc<Alphabet_size);
		return singletEmissionScores[nuc][state];
	}
	inline int GetNumChildren (State state) const {
		assert(state>=0 && state<GetNumStates());
		assert(stateInfoVector[state].numChildren<=MAX_CHILDREN); // seems a reasonable place to validate this property (although, I suppose I should check closer to the actual setting of stateInfoVector[state].numChildren)
		return stateInfoVector[state].numChildren;
	}
	inline State GetNthChildState (State state,int childNum) const {
		assert(state>=0 && state<GetNumStates());
		assert(childNum>=0 && childNum<GetNumChildren(state));
		return stateInfoVector[state].firstChild + childNum;
	}
	inline ScoreType GetNthChildTsc (State state,int childNum) const {
		assert(state>=0 && state<GetNumStates());
		assert(childNum>=0 && childNum<GetNumChildren(state));
		return stateInfoVector[state].tsc[childNum];
	}
	inline int GetChildNum_Slow(State fromState,State toState) const {
		assert(fromState>=0 && fromState<GetNumStates());
		assert(toState>=0 && toState<GetNumStates());
		for (int childNum=0; childNum<GetNumChildren(fromState); childNum++) {
			if (GetNthChildState(fromState,childNum)==toState) {
				return childNum;
			}
		}
		assert(false); // wasn't a valid transition -- toState is not a direct child of fromState
		return -1;
	}

	// actual probs
	inline ScoreType GetNthChildTransitionProb (State state,int childNum) const {
		return forwardAlgInfoVector[state].transitionProbs[childNum];
	}
	inline ScoreType GetSingletEmissionProb (State state,int nuc) const {
		return forwardAlgInfoVector[state].singletEmissionProbs[nuc];
	}
	// setting -- hey it's easy, & convenient
	inline void SetNthChildTransitionProb (State state,int childNum,ScoreType score) {
		forwardAlgInfoVector[state].transitionProbs[childNum]=score;
	}
	inline void SetSingletEmissionProb (State state,int nuc,ScoreType score) {
		forwardAlgInfoVector[state].singletEmissionProbs[nuc]=score;
	}
	// just so I don't put the explicit double->float conversions elsewhere & forget about them if I want to change it (though there are undoubtably many littered around the code)
	inline void SetNthChildTransitionProb (State state,int childNum,double score) {
		forwardAlgInfoVector[state].transitionProbs[childNum]=ScoreHelperType::DoubleToScoreType(score);
	}
	inline void SetSingletEmissionProb (State state,int nuc,double score) {
		forwardAlgInfoVector[state].singletEmissionProbs[nuc]=ScoreHelperType::DoubleToScoreType(score);
	}

	inline const std::list<CovarianceModel::State>& GetCmState (State state) const {
		return hmm2CmStateVector[InfernalHmm::IntToState(state)];
	}

	// accessor functions for doLocal
	inline bool DoLocal (void) const {
		return doLocal;
	}
	inline bool IsRightState (State state) const {
		return stateInfoVector[state].isRightState;
	}
	inline float GetLeftwardBeginsc (State state) const {
		return localStateInfoVector[state].leftwardBeginsc;
	}
	inline float GetRightwardBeginsc (State state) const {
		return localStateInfoVector[state].rightwardBeginsc;
	}
	inline int GetNumEndscLinksToLeft (State state) const {
		return localStateInfoVector[state].endscLinksToLeft.size();
	}
	inline State GetEndscLinksToLeft_State (State state,int linkNum) const {
		return localStateInfoVector[state].endscLinksToLeft[linkNum].hmmLeftState;
	}
	inline float GetEndscLinksToLeft_Endsc (State state,int linkNum) const {
		return localStateInfoVector[state].endscLinksToLeft[linkNum].endsc;
	}
};

// attempt to make scanning faster by using un-C++'d version of HmmType1.  Didn't affect performance.
struct HmmType1_OldSchool {
	typedef int State;

	int numStates;
	float **esc; // [nuc][state]
	enum {MAX_CHILDREN=4};
	struct StateInfo {
		int numChildren;
		int firstChild;
		float tsc[MAX_CHILDREN];
		bool isEmitting;
	};
	StateInfo *stateInfo;

	HmmType1_OldSchool (void);
	~HmmType1_OldSchool ();
	void Init (const HmmType1& hmm);
};


// Count # of times transitions were used in a series of Viterbi parses.  Was used for the deprecated EM-ish alg.
class TransitionCounter {
protected:
	vector<vector<int> > transitionCounts;
	vector<vector<int> > emitCounts;
	int numSamples;
	const HmmType1& hmm;
public:
	TransitionCounter (const HmmType1& _hmm);
	~TransitionCounter ();

	// get actual counts
	// "_Unreversed" refers to the fact that the hmm gets reversed because that made
	// more sense for the HmmScan function.  This means that the HMM-building code uses one direction, but the 
	// scanning code uses another.  These functions are used by the HMM-building code, so
	// they use the un-reversed HMM.  It's irrelevant to emission, tho (since it's the same state in either
	// direction), so doesn't matter.
	int GetTransitionCount_Unreversed (int unreversed_fromState,int unreversed_toState) const;
	int GetTransitionCount_Unreversed (CovarianceModel::State unreversed_fromState,CovarianceModel::State unreversed_toState) const;
	int GetEmitCount (int state,int nuc) const;
	int GetEmitCount (CovarianceModel::State state,int nuc) const;

	// Okay, I'm stopping with the pseudocount stuff -- I'll add something into the final #s
	// get frequencies, with optional pseudocounts (e.g. pseudocount=1).  pseudocounts are applied to each edge (which I think is reasonable, but I'm not sure).  Note that with pseudocounts, the forward & reverse direction can be different, but I'm not sure this really matters.  Perhaps a better idea with pseudocounts is to add say 0.0001 to each inequality's weight once we've got everything (since the inequality weights are kind of probabilities)
	// GetTraversalProbability_Unreversed gets Pr(going into toState as opposed to states [alternateToStateFirst,alternateToStateeLast)), starting at the beginning of the HMM.  In other words, if [alternateToStateFirst,alternateToStateLast) describe the "Normal" states within some node, then what is the probability of an optimal parse going thru toState?   BTW, requires that toState is in the alternateToState list
	double GetEntryProbability_Unreversed (int unreversed_toState,int unreversed_alternateToStateFirst,int unreversed_alternateToStateLast) const;
	double GetEntryProbability_Unreversed (CovarianceModel::State unreversed_toState,CovarianceModel::State unreversed_alternateToStateFirst,CovarianceModel::State unreversed_alternateToStateLast) const;
	// gets empirical Pr(transitioning to this childNum | start in state)
	double GetTransitionFrequency_Unreversed (int unreversed_fromState,int unreversed_toState,double pseudocount=0) const;
	double GetTransitionFrequency_Unreversed (CovarianceModel::State unreversed_fromState,CovarianceModel::State unreversed_toState,double pseudocount=0) const;
	// gets empirical Pr(emitting this nuc | in state)
	double GetEmitFrequency (int state,int nuc,double pseudocount=0) const;
	double GetEmitFrequency (CovarianceModel::State state,int nuc,double pseudocount=0) const;
	// a sample is one full alignment
	int GetNumSamples (void) const;

	void AddTransition (int state,int childNum);
	void AddEmit (int state,int nuc);
	// a sample is one full alignment
	void AddSample (void);

	void Dump (FILE *out);
};

// in ScanHMM_NonTemplated.cpp
extern void ScanHmm_HmmType1Float_NonTemplated (HitList& hmmHitList,const HitList& inputHitList,const HmmType1& hmm,float minLodScoreForHit,CykscanStats& cykscanStats,const char *rnaSequence, int windowLen);
extern void ScanHmm_HmmType1Float_NonTemplated (HitList& hmmHitList,const HitList& inputHitList,const HmmType1& hmmType1,const HmmType1_OldSchool& hmm,float minLodScoreForHit,CykscanStats& cykscanStats,const char *rnaSequence, int windowLen);

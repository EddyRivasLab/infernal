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

struct SolveScoresPath {
	CovarianceModel::Node cmFirstNode,cmEndingNode;
};
typedef std::list<SolveScoresPath> SolveScoresPathList;

// collects the info we use as we build an HMM model
struct HmmAndBuildInfo {
	InfernalHmm hmm;
	Cm2HmmStateVector cm2HmmState;
	Hmm2CmStateVector hmm2CmStateVector;
	CovarianceModel::State leftToRightPassthruState; // the old end of the CM block, which is now in the middle of the HMM

	SolveScoresPathList solveScoresPathList; // remember this info so we can solve scores easily at the end
};

typedef VectorByCmState<std::vector<int> > TransitionToVariableNumVector; // first dim is HMM state, then transition #, and we get the variable number for the Linear Program.  If variable number is -1, then we don't use LP to find it (and just set it 0)
struct EmissionInfo {
	InfernalHmm::State state;
	int nuc;
};
struct TransitionOrEmissionInfo {
	bool isEmission; // else is transition
	bool isUsed; // sanity check

	InfernalHmm::EdgeInfo edgeInfo;
	EmissionInfo emissionInfo;
};
typedef std::vector<TransitionOrEmissionInfo> TransitionOrEmissionInfoVector;

extern void GetGlobalVarsFromInfernalHmm (std::vector<double>& globalVars,const InfernalHmm& infernalHmm,const TransitionOrEmissionInfoVector& transitionOrEmissionInfoVector);
extern void SetGlobalVarsIntoInfernalHmm (InfernalHmm& infernalHmm,const std::vector<double>& globalVars,const TransitionOrEmissionInfoVector& transitionOrEmissionInfoVector);

struct ScoreVariablesInfo {
	TransitionToVariableNumVector transitionToVariableNumVector,emissionToVariableNumVector;
	TransitionOrEmissionInfoVector globalVariableToTransitionOrEmissionVector;
	int numVariables;

};
struct InequalityTerm {
	// no need for coefficient, since they're all +1.
	int variableNum;
	bool operator < (const InequalityTerm& t) const {
		return variableNum<t.variableNum;
	}
	bool operator == (const InequalityTerm& t) const {
		return variableNum==t.variableNum;
	}
	bool operator != (const InequalityTerm& t) const {
		return variableNum!=t.variableNum;
	}
};
enum InequalityType {
	IneqType_GE,IneqType_Less
};
struct InequalityBase {
	// all inequalities are lhs >= rhs
	std::list<InequalityTerm> lhs;
	float rhs;

	// this code is used by the old implementation of the inf-len forward alg; now I evaluate it fully symbolically
	CovarianceModel::State pathStartState,pathEndState;
	float sumOfConstantsInHmm; // for convenience if we're tracing paths in the HMM.  The inequalities themselves mention any variables in the path, but some transitions don't have associated variables, so this field stores the sum of these constants (and remember that the constants represent logs)
	InfernalHmm::StateList hmmInsertStatesInPath; // again, for convenience of adding scores in paths; in this case, we technically have to consider unbounded-length paths thru insert states

	// # of nucs of each type emitted on this path
	// since insert states don't emit on the path, they're not reflected here
	int nucEmitCount[MAXABET];

	// 'weight' is used by a defunct method to optimize HMMs ("EMish"), which is not as good as the inf-len forward alg
	double weight; // 'weight' is the weight of this inequality's slack variable in the objective function.  weight==1.0 if we're not weighting the inequalities

	InequalityType inequalityType; // for non-standard inequalities

	// convenience functions
	bool ContainsLocalVar (int localVar) const {
		InequalityTerm term;
		term.variableNum=localVar;
		return std::find(lhs.begin(),lhs.end(),term)!=lhs.end();
	}
};
struct Inequality : public InequalityBase {
	Inequality (void) {
		inequalityType=IneqType_GE; // what we normally use
	}
	~Inequality () {
	}
	void operator = (const Inequality& t) {
		InequalityBase::operator =(t);
	}
	Inequality (const Inequality&t) {
		*this=t;
	}
};
typedef std::list<Inequality> InequalityList;
typedef std::list<InequalityList> InequalityListList;

struct InequalitiesAndLocalVariables {
	std::vector<int> globalToLocalVariables;
	std::vector<int> localToGlobalVariables;
	int numLocalVariables;
	InequalityList inequalityList;
};
typedef VectorByCmNode<InequalitiesAndLocalVariables> InequalitiesAndLocalVariablesByCmNode;

struct ExtraCm2HmmInfo {
	// outputs of Cm2Hmm
	InequalitiesAndLocalVariablesByCmNode inequalitiesAndLocalVariables;
	ScoreVariablesInfo scoreVariablesInfo;

	// inputs
	bool actuallySolveScores; // set this to false to save time when we only want the inequalities, not the actual solutions

	ExtraCm2HmmInfo (void) {
		actuallySolveScores=true;
	}
};

extern void AddInequalitiesInTermsOfGlobalVars(InequalityList& inequalityList,const ExtraCm2HmmInfo& extraCm2HmmInfo,const CovarianceModel::Node cmNode);
extern void AddCrossProductOfInequalityLists(InequalityList& inequalityList,InequalityListList& inequalityListList);
extern void ConvertInequalityListVarNums (InequalityList& inequalityList,const std::vector<int>& varNumMapping);

extern void GetGlobalVariableNumsForNode (std::list<int>& globalVarNums,const CovarianceModel& cm,InfernalHmm& infernalHmm,const CovarianceModel::Node cmNode,const ExtraCm2HmmInfo& extraCm2HmmInfo);
extern void GetLocalVariablesValueForNode(std::vector<float>& localVariablesToValue,const CovarianceModel& cm,InfernalHmm& infernalHmm,const CovarianceModel::Node cmNode,const ExtraCm2HmmInfo& extraCm2HmmInfo);
extern void SetLocalVariablesValueForNode(std::vector<float>& localVariablesToValue,const CovarianceModel& cm,InfernalHmm& infernalHmm,const CovarianceModel::Node cmNode,const ExtraCm2HmmInfo& extraCm2HmmInfo);


class TemporarilyModifyInequality {
protected:
	Inequality& inequalitySoFar;
	int numVariablesAdded;
	float startingScore;
	double startingWeight;
	float starting_sumOfConstantsInHmm;
	size_t starting_hmmInsertStatesInPath_size;
	int startingNucEmitCount[MAXABET];
public:
	TemporarilyModifyInequality (Inequality& _inequalitySoFar);
	~TemporarilyModifyInequality ();

	void AddScore (float addToScore);
	void AddVariable (int globalVariableNum,std::vector<int>& globalToLocalVariables,int& numlocalVariables);
	void MultiplyWeight (double mult);
	void PushInsertState (InfernalHmm::State insertState);
};

struct WeightedInequalitiesInfo {
	const TransitionCounter *transitionCounter;
	double pseudocount; // should be 1, but maybe 0.5 makes sense.  Setting this to 0 is very dangerous, since then you can get divide by 0 for states that were never visited.  While in the ensuing debugging, you might miss your bus, and be stranded at school forever.  (Not really; it's only a 30-minute walk.)
};

// HMM committees (a series of HMMs to apply) are deprecated, since it doesn't help much when the HMMs are of the same type (i.e. all compact or all expanded) and are optimized with the inf-len forward alg.  Anyway, the same functionality can be specified on the command line.
class HmmCommittee {
protected:
	std::string description;
	int recommendedCommitteeSize; // is computed by evaluate_everyone, otherwise just unknown (0)
public:
	int GetRecommendedCommitteeSize (void) const; // will 'Die' if unknown

	typedef std::list<InfernalHmm> HmmList;
	HmmList hmmList;

	HmmCommittee (void);
	~HmmCommittee ();

	bool /*success*/ LoadInBinary (const char *fileName,int effectiveCommitteeSize);
	bool SaveInBinary (const char *fileName,const char *additionalDescription);

	void DumpMemberBuildInfo (FILE *out);

	// information stored on a linear program, so when we want to solve it again, we can generate a new solution.  For now, I'm just planning on finding all solutions at once, and then picking from them.  Another simplification, is that I assume the only distinction between solutions is which constraints they satisfy perfectly (i.e. which slack variables are 0).  In practice, there seem to be multiple solutions in this sense, but clearly you should be able to generate intermediates.
	struct LinearProgramStoredInfo {

		unsigned int numTimesGeneratedSolution; // our pattern is to first go thru all solutions, and then to pick randomly
		std::vector<std::vector<float> > localVariablesToValuePerSolution; // first dimension is solution #, 2nd is the local variables' values in that solution
		std::vector<float> avgInflationPerSolution;

		/*
		   std::vector<double> sumOfSlack; // each constraint has a slack variable.  This is the sum of the slack variables in the optimal solutions generated so far
		_Bvector slackHasBeenZero; // each constraint may have been solved perfectly in one of the solutions, if so slack has zero is set to true

		double optimalValueOfObjectiveFunc; // when we come up with alternate solutions, we want them to all be just as optimal, i.e. have the same objective func value

		bool isFirstSolutionFound;
		LinearProgramStoredInfo (void) {
			isFirstSolutionFound=false;
		}
		*/
	};
	// allows me to add ideas later.  Usually passed around as a pointer, where NULL implies no need to do anything about committee.  Stores input parameters, and work-in-progress
	struct CommitteeBuildInfo {
		VectorByCmNode<LinearProgramStoredInfo> cmStartNodeToLinearProgramInfo;
		int currCommitteeMemberNum; // first is 0, then 1, ...
	};
};

extern void Cm2Hmm (InfernalHmm& hmm,Cm2Hmm_HmmBuildType hmmBuildType,const CovarianceModel& sourceCM,const char *cmFileName,const std::string& programParams,bool forceCreate=false); // cmFileName just used for caching the translation to an HMM

extern void Cm2Hmm (HmmCommittee& hmmCommittee,const CovarianceModel& sourceCM,const char *cmFileName,int committeeSize,int effectiveCommitteeSize);

// slightly lower-level function
void Cm2Hmm_WithWeighting_NoCaching (InfernalHmm& hmm,Cm2Hmm_HmmBuildType hmmBuildType,const CovarianceModel& sourceCM,const char *cmFileName,WeightedInequalitiesInfo *weightedInequalitiesInfo,ExtraCm2HmmInfo *extraCm2HmmInfo=NULL,const int nodesToSpanWhileSolvingScores=1);

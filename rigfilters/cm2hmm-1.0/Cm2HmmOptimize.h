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

void VectorFloatToDouble (vector<double>& varsDouble,const vector<float>& varsFloat);
void VectorDoubleToFloat (vector<float>& varsFloat,const vector<double>& varsDouble);
void ExponentiateVariables(vector<double>& varsExp,const vector<double>& vars);

// abstract class to represent an objective func, with gradient & hessian, so I can re-use code
class ObjectiveFunc {
public:
	virtual double EvalActualObjectiveFuncLog2 (const vector<double>& problemVars);
	// convenience
	double EvalValueOnly (const vector<double>& problemVars);

	virtual ~ObjectiveFunc ();
	virtual void Eval (double& f,vector<double>& gradient,vector2d<double>& hessian,const vector<double>& problemVars,bool calculateHessian=true,bool calculateGradient=true) = 0;
	virtual int GetNumProblemVars (void) = 0;
	virtual void LocalToProblemVars (vector<double>& problemVars,const vector<float>& localVars) = 0;
	virtual void ProblemToLocalVars (vector<float>& localVars,const vector<double>& problemVars) = 0;
	virtual const InequalityList& GetInequalityList (void) = 0;
};

// adapts an ObjectiveFunc object to return the log(obj-func)
class LnObjectiveFuncAdaptor : public ObjectiveFunc {
protected:
	ObjectiveFunc *adaptee;
public:
	// adaptee is owned by this class & freed in destructor
	LnObjectiveFuncAdaptor (ObjectiveFunc *adaptee_);
	~LnObjectiveFuncAdaptor ();

	void Eval (double& f,vector<double>& gradient,vector2d<double>& hessian,const vector<double>& problemVars,bool calculateHessian,bool calculateGradient);
	int GetNumProblemVars (void);
	void LocalToProblemVars (vector<double>& problemVars,const vector<float>& localVars);
	void ProblemToLocalVars (vector<float>& localVars,const vector<double>& problemVars);
	const InequalityList& GetInequalityList (void);
};

// class to capture the objective function of a node by running an HMM forward alg on a "number" class that symbolically represents the vars.  Currently, it's only set up to work with the infinite-length forward algorithm
class SymbolicProbVariableMath : public SymbolicMath {
	friend class Hmm;
protected:
	vector<int> problemToGlobalVars;
	void CreateGlobalVars(vector<double>& globalVars,const vector<double>& problemVars);
public:
	SymbolicProbVariableMath (const vector<int>& problemToGlobalVars_,const vector<int>& globalToProblemVars,Expression expression); // for more general functionality
	SymbolicProbVariableMath (const vector<int>& problemToGlobalVars_,const vector<int>& globalToProblemVars,const CovarianceModel& cm_,const InfernalHmm& infernalHmm_,const ScoreVariablesInfo& scoreVariablesInfo_,bool applyLog2=false);
	~SymbolicProbVariableMath ();

	int GetNumProblemVars (void);

	// a version of HMM that returns "numbers" of type Expression
	class Hmm {
	protected:
		vector<int> globalToProblemVars; // vector mapping each global variable to a problem variable.  If entry is -1, then that global variable is not used.
		const HmmType1& scanningHmm;
		const InfernalHmm& infernalHmm;
		const ScoreVariablesInfo& scoreVariablesInfo;

		SymbolicProbVariableMath::Expression GetVarIfUsedElseConst (int globalVar,double constVal) const;
	public:
		Hmm (const vector<int>& globalToProblemVars_,const HmmType1& scanningHmm_,const InfernalHmm& infernalHmm_,const ScoreVariablesInfo& scoreVariablesInfo_);
		~Hmm ();

		typedef int State;
		int GetNumStates (void) const;
		State GetFirstState (void) const;
		State GetLastState (void) const;
		State GetActualLastState (void) const;
		static State GetInvalidState (void);
		bool IsEmittingState (State state) const;
		int GetNumChildren (State state) const;
		State GetNthChildState (State state,int childNum) const;
		SymbolicProbVariableMath::Expression GetNthChildTransitionProb (State state,int childNum) const;
		SymbolicProbVariableMath::Expression GetSingletEmissionProb (State state,int nuc) const;
	};

	class GenericObjectiveFunc : public ::ObjectiveFunc {
	protected:
		SymbolicProbVariableMath& master;
		InequalityList inequalityList;
	public:
		GenericObjectiveFunc (SymbolicProbVariableMath& master_);
		GenericObjectiveFunc (SymbolicProbVariableMath& master_,const InequalityList& inequalityList_);
		~GenericObjectiveFunc ();

		void Eval (double& f,vector<double>& gradient,vector2d<double>& hessian,const vector<double>& problemVars,bool calculateHessian,bool calculateGradient);
		int GetNumProblemVars (void);
		void LocalToProblemVars (vector<double>& problemVars,const vector<float>& localVars); // function isn't implemented, since this class doesn't know about local vars
		void ProblemToLocalVars (vector<float>& localVars,const vector<double>& problemVars);
		const InequalityList& GetInequalityList (void);
	};
};

class GenericSymbolicObjectiveFunc : public ::ObjectiveFunc {
protected:
	SymbolicMath& master;
	InequalityList inequalityList;
	int numProblemVars;
public:
	GenericSymbolicObjectiveFunc (SymbolicMath& master_,const InequalityList& inequalityList_,int numProblemVars_);
	~GenericSymbolicObjectiveFunc ();

	void Eval (double& f,vector<double>& gradient,vector2d<double>& hessian,const vector<double>& problemVars,bool calculateHessian,bool calculateGradient);
	int GetNumProblemVars (void);
	void LocalToProblemVars (vector<double>& problemVars,const vector<float>& localVars); // function isn't implemented, since this class doesn't know about local vars
	void ProblemToLocalVars (vector<float>& localVars,const vector<double>& problemVars);
	const InequalityList& GetInequalityList (void);
};

class ProblemDefinition {
public:
	virtual ~ProblemDefinition ();
	virtual ObjectiveFunc *NewObjectiveFunc (void) = 0;
};
class ProblemDefinitionInstantiator {
public:
	virtual ~ProblemDefinitionInstantiator ();
	virtual ProblemDefinition *NewProblemDefinition (const InfernalHmm& sourceHmm,const CovarianceModel& cm,int windowLen,CovarianceModel::Node cmNode,const ExtraCm2HmmInfo& extraCm2HmmInfo,char *fastaFileToSearch,int numLocalVariables,const vector<float>& localVars) = 0;
};

// uses symbolic math to do objective function for a single node
class SymbolicObjective_OneNode_ForwardAlg : public ProblemDefinition {
protected:

	const ExtraCm2HmmInfo& extraCm2HmmInfo;
	CovarianceModel::Node cmNode;

	SymbolicProbVariableMath *symbolic;
public:
	SymbolicObjective_OneNode_ForwardAlg (const InfernalHmm& sourceHmm,const CovarianceModel& cm,int windowLen,const CovarianceModel::Node cmNode,const ExtraCm2HmmInfo& extraCm2HmmInfo,const char *fastaFileToSearch,int _numLocalVariables);
	~SymbolicObjective_OneNode_ForwardAlg ();

	class ObjectiveFunc : public ::ObjectiveFunc {
		SymbolicProbVariableMath& symbolic;
		SymbolicObjective_OneNode_ForwardAlg& master;
		int numLocalVariables;
		const InequalityList& inequalityList;
	public:
		ObjectiveFunc (SymbolicObjective_OneNode_ForwardAlg& _master,SymbolicProbVariableMath& symbolic_,int numLocalVariables_,const InequalityList& inequalityList_);
		void Eval (double& f,vector<double>& gradient,vector2d<double>& hessian,const vector<double>& problemVars,bool calculateHessian,bool calculateGradient);
		int GetNumProblemVars (void);
		void LocalToProblemVars (vector<double>& problemVars,const vector<float>& localVars);
		void ProblemToLocalVars (vector<float>& localVars,const vector<double>& problemVars);
		const InequalityList& GetInequalityList (void);
	};
	::ObjectiveFunc *NewObjectiveFunc (void);

	class ProblemDefinitionInstantiator : public ::ProblemDefinitionInstantiator {
	public:
		ProblemDefinitionInstantiator (void);
		~ProblemDefinitionInstantiator ();
		ProblemDefinition *NewProblemDefinition (const InfernalHmm& sourceHmm,const CovarianceModel& cm,int windowLen,CovarianceModel::Node cmNode,const ExtraCm2HmmInfo& extraCm2HmmInfo,char *fastaFileToSearch,int numLocalVariables,const vector<float>& localVars);
	};
};

// uses the symbolic math stuff to define the global problem (rather than a per-node problem)
class GlobalForwardInfSymbolicObjectiveFunc : public ObjectiveFunc {
protected:
	SymbolicProbVariableMath *symbolic;
	InequalityList inequalityList; // in terms of globalVars
	int numGlobalVars;
public:
	GlobalForwardInfSymbolicObjectiveFunc (const InfernalHmm& sourceHmm,const CovarianceModel& cm,const ExtraCm2HmmInfo& extraCm2HmmInfo,int numAdjacentNodesToMerge);
	~GlobalForwardInfSymbolicObjectiveFunc ();
	void Eval (double& f,vector<double>& gradient,vector2d<double>& hessian,const vector<double>& problemVars,bool calculateHessian,bool calculateGradient);
	int GetNumProblemVars (void);
	void LocalToProblemVars (vector<double>& problemVars,const vector<float>& localVars);
	void ProblemToLocalVars (vector<float>& localVars,const vector<double>& problemVars);
	const InequalityList& GetInequalityList (void);
};

class NodeCombinerForwardInfSymbolicObjectiveFunc : public ObjectiveFunc {
protected:
	SymbolicProbVariableMath *symbolic;
	InequalityList inequalityList;

	int numGlobalVars,numProblemVars;
	vector<int> problemToGlobalVars,globalToProblemVars;
public:
	NodeCombinerForwardInfSymbolicObjectiveFunc (InfernalHmm& sourceHmm,const CovarianceModel& cm,const ExtraCm2HmmInfo& extraCm2HmmInfo,int numAdjacentNodesToMerge,int maxNodesAtATime,const CovarianceModel::Node cmStartNode,bool applyLog2=false);
	~NodeCombinerForwardInfSymbolicObjectiveFunc ();

	void Eval (double& f,vector<double>& gradient,vector2d<double>& hessian,const vector<double>& problemVars,bool calculateHessian,bool calculateGradient);
	int GetNumProblemVars (void);
	void LocalToProblemVars (vector<double>& problemVars,const vector<float>& localVars);
	void ProblemToLocalVars (vector<float>& localVars,const vector<double>& problemVars);
	const InequalityList& GetInequalityList (void);

	void GlobalToProblemVars (vector<double>& problemVars,const vector<double>& globalVars);
	void UpdateGlobalVarsFromProblemVars (vector<double>& globalVars,const vector<double>& problemVars);
};

// wraps a solver library, such as Opt++ or CFSQP
class SolverWrapper {
public:
	SolverWrapper (void);
	virtual ~SolverWrapper ();

	class MessageReceiver {
	public:
		virtual ~MessageReceiver ();
		virtual void EvaluatedObjectiveFunc (double functionValue,const vector<double>& problemVars); // default: do nothing
	};

	virtual vector<double> /* optimal problem vars */ Solve (ObjectiveFunc *objectiveFunc,const vector<double>& inputProblemVars,double maxVariableMagnitudeForUpperLowerBounds,bool importantBoundsAreSet=false,double importantLowerBoundAllVars=0,double importantUppderBoundAllVars=0,MessageReceiver *messageReceiver=NULL) = 0;
};

SolverWrapper *NewSolverWrapper_cfsqp (int B,int C);
SolverWrapper *NewSolverWrapper_OptNIPS (void);
SolverWrapper *MakeSolverWrapper (int& a,const int argc,char **argv);

extern float ComputeLhs (const vector<float>& localVariablesToValue,const Inequality& ineq);

extern void GlobalHmmOptimizer (char *cmFileName,bool doLocalAlignment,const std::string& programParams,int numAdjacentNodesToMerge);

extern void HmmOptimizer_NodeCombiner (char *cmFileName,bool doLocalAlignment,const std::string& programParams,int numAdjacentNodesToMerge,int maxNodesAtATime,int numIters,Cm2Hmm_HmmBuildType hmmType,SolverWrapper *solverWrapper,bool saveHmmsInProgress=true,const char *hmmSaveFileName="combiner-hmm.bin",HmmFileFormat hmmFileFormat=HmmFileFormat_Binary); // defaults are for benefit of legacy code


// actually defined in ForwardHMM.cpp
extern SymbolicProbVariableMath::Expression InfiniteLengthForwardAlg_Symbolic (const SymbolicProbVariableMath::Hmm& hmm,MarkovModelStats& markovModelStats);

extern void MakeCmHmmForOptimization (CovarianceModel& cm,InfernalHmm& infernalHmm,ExtraCm2HmmInfo& extraCm2HmmInfo,char *cmFileName,bool doLocalAlignment,const std::string& programParams,Cm2Hmm_HmmBuildType hmmType=HmmBuildType_Original);
extern float ComputeLhs (const vector<float>& localVariablesToValue,const Inequality& ineq);

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
#include "stdafx.h"

#include "UseDebugNew.h"

#include "cmzasha.h"

#include "NoUnderflowDouble.h"

template <class Real1,class Real2>
void Vector_ConvertReals (std::vector<Real1>& v1,const std::vector<Real2>& v2)
{
	v1.resize(v2.size());
	for (size_t i=0; i<v2.size(); i++) {
		v1[i]=(Real1)v2[i];
	}
}

void VectorFloatToDouble (std::vector<double>& varsDouble,const std::vector<float>& varsFloat)
{
	Vector_ConvertReals(varsDouble,varsFloat);
	/*
	varsDouble.resize(varsFloat.size());
	for (size_t i=0; i<varsFloat.size(); i++) {
		varsDouble[i]=varsFloat[i];
	}
	*/
}
void VectorDoubleToFloat (std::vector<float>& varsFloat,const std::vector<double>& varsDouble)
{
	Vector_ConvertReals(varsFloat,varsDouble);
	/*
	varsFloat.resize(varsDouble.size());
	for (size_t i=0; i<varsDouble.size(); i++) {
		varsFloat[i]=(float)(varsDouble[i]);
	}
	*/
}
template <class Real>
void ExponentiateVariables_templ(std::vector<Real>& varsExp,const std::vector<Real>& vars)
{
	varsExp.resize(vars.size());
	Real log_of_2=log((Real)(2));
	for (size_t i=0; i<vars.size(); i++) {
		varsExp[i]=exp(vars[i]*log_of_2);
	}
}
void ExponentiateVariables(std::vector<double>& varsExp,const std::vector<double>& vars)
{
	return ExponentiateVariables_templ<double>(varsExp,vars);
}
void ExponentiateVariables(std::vector<NoUnderflowDouble>& varsExp,const std::vector<NoUnderflowDouble>& vars)
{
	return ExponentiateVariables_templ<NoUnderflowDouble>(varsExp,vars);
}

////////////////////////
// ObjectiveFunc

ObjectiveFunc::~ObjectiveFunc ()
{
}
double ObjectiveFunc::EvalActualObjectiveFuncLog2 (const std::vector<double>& problemVars)
{
	double fx;
	std::vector<double> gradient;
	vector2d<double> hessian;
	Eval(fx,gradient,hessian,problemVars,false);
	return fx/log(2.0);
}
double ObjectiveFunc::EvalValueOnly (const std::vector<double>& problemVars)
{
	double fx;
	std::vector<double> gradient;
	vector2d<double> hessian;
	Eval(fx,gradient,hessian,problemVars,false,false);
	return fx;
}

////////////////////////
// LnObjectiveFuncAdaptor

LnObjectiveFuncAdaptor::LnObjectiveFuncAdaptor (ObjectiveFunc *adaptee_)
{
	adaptee=adaptee_;
}
LnObjectiveFuncAdaptor::~LnObjectiveFuncAdaptor ()
{
	delete adaptee;
}
void LnObjectiveFuncAdaptor::Eval (double& f,std::vector<double>& gradient,vector2d<double>& hessian,const std::vector<double>& problemVars,bool calculateHessian,bool calculateGradient)
{
	double fX;
	adaptee->Eval(fX,gradient,hessian,problemVars,calculateHessian,calculateGradient);

	f=log(fX);

	if (calculateHessian) {
		// do hessian first, before we clobber gradient
		// d ((df(X)/dx_i)/f(X) / dx_j   [formula from gradient, now differentiating wrt x_j]
		// = (d(df(X)/dx_i)/dx_j)/f(X)   +    -(df(X)/dx_i)*(df(X)/dx_j)/(f(X)^2)
		// so hessian = hessian(f(X))/f(X)  - gradient(f(X))*gradient(f(X))' /f(X)^2
		for (int i=0; i<(int)gradient.size(); i++) {
			for (int j=0; j<(int)gradient.size(); j++) {
				hessian[i][j] /= fX;
				hessian[i][j] -= gradient[i]*gradient[j]/(fX*fX);
			}
		}
	}

	// f(X) = objective function
	// dln(f(X))/dx_i 
	// = dln(f(x))/df(x) * df(X)/dx_i = (df(X)/dx_i) / f(X)
	// so gradient = gradient(f(X))/f(X)
	for (size_t i=0; i<gradient.size(); i++) {
		gradient[i] /= fX;
	}
}
int LnObjectiveFuncAdaptor::GetNumProblemVars (void)
{
	return adaptee->GetNumProblemVars();
}
void LnObjectiveFuncAdaptor::LocalToProblemVars (std::vector<double>& problemVars,const std::vector<float>& localVars)
{
	return adaptee->LocalToProblemVars(problemVars,localVars);
}
void LnObjectiveFuncAdaptor::ProblemToLocalVars (std::vector<float>& localVars,const std::vector<double>& problemVars)
{
	return adaptee->ProblemToLocalVars (localVars,problemVars);
}
const InequalityList& LnObjectiveFuncAdaptor::GetInequalityList (void)
{
	return adaptee->GetInequalityList();
}

////////////////////////
// ProblemDefinition

ProblemDefinition::~ProblemDefinition ()
{
}

ProblemDefinitionInstantiator::~ProblemDefinitionInstantiator ()
{
}

///////////////////////
// SolverWrapper

SolverWrapper::SolverWrapper (void)
{
}
SolverWrapper::~SolverWrapper ()
{
}

SolverWrapper::MessageReceiver::~MessageReceiver ()
{
}
void SolverWrapper::MessageReceiver::EvaluatedObjectiveFunc (double functionValue,const std::vector<double>& problemVars)
{
}

SolverWrapper *MakeSolverWrapper (int& a,const int argc,char **argv)
{
	AnotherParam(a,argc);
	const char *method=argv[a];
	a++;
	if (strcmp(method,"OptNIPS")==0) {
	  /*
		printf("WARNING: OptNIPS can be very slow compared to CFSQP, I believe since it has to manipulate the hessian matrix, and since the code has to find the O(n^2) double derivatives for the hessian.\n");
		return NewSolverWrapper_OptNIPS();
	  */
	  	printf("ERROR: Depracating OptNIPS - this solver is not availabe!\n");
	}
	if (strcmp(method,"cfsqp")==0) {
		AnotherParam(a,argc);
		int B=atoi(argv[a]);
		a++;
		AnotherParam(a,argc);
		int C=atoi(argv[a]);
		a++;
		return NewSolverWrapper_cfsqp(B,C);
	}
	throw SimpleStringException("unknown solver: \"%s\"",method);
}


//////////////////////////
// SymbolicProbVariableMath

SymbolicProbVariableMath::SymbolicProbVariableMath (const std::vector<int>& problemToGlobalVars_,const std::vector<int>& globalToProblemVars,Expression expression)
{
	problemToGlobalVars=problemToGlobalVars_;
	SetRootExpression(expression);
}
SymbolicProbVariableMath::SymbolicProbVariableMath (const std::vector<int>& problemToGlobalVars_,const std::vector<int>& globalToProblemVars,const CovarianceModel& cm,const InfernalHmm& infernalHmm,const ScoreVariablesInfo& scoreVariablesInfo,bool applyLog2)
{
	problemToGlobalVars=problemToGlobalVars_;

	HmmType1 scanningHmm;
	scanningHmm.Init(infernalHmm);

	Hmm symbolicHmm(globalToProblemVars,scanningHmm,infernalHmm,scoreVariablesInfo);
	SymbolicProbVariableMath::Expression rootExpression=InfiniteLengthForwardAlg_Symbolic (symbolicHmm,*defaultTrainingMarkovModelStats);
	if (applyLog2) {
		rootExpression=Expression::Log2(rootExpression);
	}

	SetRootExpression(rootExpression);
}
SymbolicProbVariableMath::~SymbolicProbVariableMath ()
{
}
void SymbolicProbVariableMath::CreateGlobalVars(std::vector<double>& globalVars,const std::vector<double>& problemVars)
{
	for (size_t i=0; i<problemVars.size(); i++) {
		globalVars[problemToGlobalVars[i]]=problemVars[i];
	}
}
int SymbolicProbVariableMath::GetNumProblemVars (void)
{
	return (int)(problemToGlobalVars.size());
}

SymbolicProbVariableMath::GenericObjectiveFunc::GenericObjectiveFunc (SymbolicProbVariableMath& master_)
: master(master_)
{
}
SymbolicProbVariableMath::GenericObjectiveFunc::GenericObjectiveFunc (SymbolicProbVariableMath& master_,const InequalityList& inequalityList_)
: master(master_)
, inequalityList(inequalityList_)
{
}
SymbolicProbVariableMath::GenericObjectiveFunc::~GenericObjectiveFunc ()
{
}
void SymbolicProbVariableMath::GenericObjectiveFunc::Eval (double& f,std::vector<double>& gradient,vector2d<double>& hessian,const std::vector<double>& problemVars,bool calculateHessian,bool calculateGradient)
{
	//printf("problemVars: "); for (size_t i=0; i<problemVars.size(); i++) {	printf("%d=%lg,  ",i,problemVars[i]);	} printf("\n");
	master.Eval(GetNumProblemVars(),f,gradient,hessian,problemVars,calculateHessian,calculateGradient);
	fprintf(stderr,"obj func = %.15lg\n",f);
	//if (calculateGradient) {		printf("gradient: "); for (size_t i=0; i<problemVars.size(); i++) {	printf("d %d=%lg,  ",i,gradient[i]);	} printf("\n");  	}
}
int SymbolicProbVariableMath::GenericObjectiveFunc::GetNumProblemVars (void)
{
	return master.GetNumProblemVars();
}
void SymbolicProbVariableMath::GenericObjectiveFunc::LocalToProblemVars (std::vector<double>& problemVars,const std::vector<float>& localVars)
{
	throw SimpleStringException("SymbolicProbVariableMath::GenericObjectiveFunc::LocalToProblemVars not implemented, %s:%d",__FILE__,__LINE__);
}
void SymbolicProbVariableMath::GenericObjectiveFunc::ProblemToLocalVars (std::vector<float>& localVars,const std::vector<double>& problemVars)
{
	throw SimpleStringException("SymbolicProbVariableMath::GenericObjectiveFunc::ProblemToLocalVars not implemented, %s:%d",__FILE__,__LINE__);
}
const InequalityList& SymbolicProbVariableMath::GenericObjectiveFunc::GetInequalityList (void)
{
	return inequalityList;
}

SymbolicProbVariableMath::Hmm::Hmm (const std::vector<int>& globalToProblemVars_,const HmmType1& scanningHmm_,const InfernalHmm& infernalHmm_,const ScoreVariablesInfo& scoreVariablesInfo_)
: globalToProblemVars(globalToProblemVars_)
, scanningHmm(scanningHmm_)
, infernalHmm(infernalHmm_)
, scoreVariablesInfo(scoreVariablesInfo_)
{
}
SymbolicProbVariableMath::Hmm::~Hmm ()
{
}
int SymbolicProbVariableMath::Hmm::GetNumStates (void) const
{
	return scanningHmm.GetNumStates();
}
SymbolicProbVariableMath::Hmm::State SymbolicProbVariableMath::Hmm::GetFirstState (void) const
{
	return scanningHmm.GetFirstState();
}
SymbolicProbVariableMath::Hmm::State SymbolicProbVariableMath::Hmm::GetLastState (void) const
{
	return scanningHmm.GetLastState();
}
SymbolicProbVariableMath::Hmm::State SymbolicProbVariableMath::Hmm::GetActualLastState (void) const
{
	return scanningHmm.GetActualLastState();
}
SymbolicProbVariableMath::Hmm::State SymbolicProbVariableMath::Hmm::GetInvalidState (void)
{
	return HmmType1::GetInvalidState();
}
bool SymbolicProbVariableMath::Hmm::IsEmittingState (State state) const
{
	return scanningHmm.IsEmittingState(state);
}
int SymbolicProbVariableMath::Hmm::GetNumChildren (State state) const
{
	return scanningHmm.GetNumChildren(state);
}
SymbolicProbVariableMath::Hmm::State SymbolicProbVariableMath::Hmm::GetNthChildState (State state,int childNum) const
{
	return scanningHmm.GetNthChildState(state,childNum);
}
SymbolicMath::Expression SymbolicProbVariableMath::Hmm::GetVarIfUsedElseConst (int globalVar,double constVal) const
{
	if (globalVar==-1) {
		return Expression(constVal);
	}
	else {
		int problemVarNum=globalToProblemVars[globalVar];
		if (problemVarNum==-1) {
			// no problem var -- so, just return a constant
			return Expression(constVal);
		}
		else {
			return Expression::ExpressionOfVarPow2(problemVarNum);
		}
	}
}
SymbolicMath::Expression SymbolicProbVariableMath::Hmm::GetNthChildTransitionProb (State state,int childNum) const
{
	const HmmType1::State toState=scanningHmm.GetNthChildState(state,childNum);
	const int infernalChildNum=infernalHmm.GetChildNum_Slow(InfernalHmm::IntToState(toState),InfernalHmm::IntToState(state));
	int globalVarNum=-1;
	if (!scoreVariablesInfo.transitionToVariableNumVector[InfernalHmm::IntToState(toState)].empty()) {
		globalVarNum=scoreVariablesInfo.transitionToVariableNumVector[InfernalHmm::IntToState(toState)][infernalChildNum];
	}
	return GetVarIfUsedElseConst(globalVarNum,scanningHmm.GetNthChildTransitionProb(state,childNum));
}
SymbolicMath::Expression SymbolicProbVariableMath::Hmm::GetSingletEmissionProb (State state,int nuc) const
{
	int globalVarNum=-1;
	if (!scoreVariablesInfo.emissionToVariableNumVector[InfernalHmm::IntToState(state)].empty()) {
		globalVarNum=scoreVariablesInfo.emissionToVariableNumVector[InfernalHmm::IntToState(state)][nuc];
	}
	return GetVarIfUsedElseConst(globalVarNum,scanningHmm.GetSingletEmissionProb(state,nuc));
}

////////////////////
// GenericSymbolicObjectiveFunc

GenericSymbolicObjectiveFunc::GenericSymbolicObjectiveFunc (SymbolicMath& master_,const InequalityList& inequalityList_,int numProblemVars_)
: master(master_)
, inequalityList(inequalityList_)
, numProblemVars(numProblemVars_)
{
}
GenericSymbolicObjectiveFunc::~GenericSymbolicObjectiveFunc ()
{
}
void GenericSymbolicObjectiveFunc::Eval (double& f,std::vector<double>& gradient,vector2d<double>& hessian,const std::vector<double>& problemVars,bool calculateHessian,bool calculateGradient)
{
	master.Eval(numProblemVars,f,gradient,hessian,problemVars,calculateHessian,calculateGradient);
}
int GenericSymbolicObjectiveFunc::GetNumProblemVars (void)
{
	return numProblemVars;
}
void GenericSymbolicObjectiveFunc::LocalToProblemVars (std::vector<double>& problemVars,const std::vector<float>& localVars)
{
	throw SimpleStringException("Not implemented %s:%d",__FILE__,__LINE__);
}
void GenericSymbolicObjectiveFunc::ProblemToLocalVars (std::vector<float>& localVars,const std::vector<double>& problemVars)
{
	throw SimpleStringException("Not implemented %s:%d",__FILE__,__LINE__);
}
const InequalityList& GenericSymbolicObjectiveFunc::GetInequalityList (void)
{
	return inequalityList;
}


////////////////////
// SymbolicObjective_OneNode_ForwardAlg
SymbolicObjective_OneNode_ForwardAlg::SymbolicObjective_OneNode_ForwardAlg (const InfernalHmm& sourceHmm,const CovarianceModel& cm,int windowLen,const CovarianceModel::Node cmNode_,const ExtraCm2HmmInfo& extraCm2HmmInfo_,const char *fastaFileToSearch,int _numLocalVariables)
: extraCm2HmmInfo(extraCm2HmmInfo_)
{
	cmNode=cmNode_;

	symbolic=new SymbolicProbVariableMath (extraCm2HmmInfo.inequalitiesAndLocalVariables[cmNode].localToGlobalVariables,extraCm2HmmInfo.inequalitiesAndLocalVariables[cmNode].globalToLocalVariables,cm,sourceHmm,extraCm2HmmInfo.scoreVariablesInfo);
}
SymbolicObjective_OneNode_ForwardAlg::~SymbolicObjective_OneNode_ForwardAlg ()
{
	delete symbolic;
}
ObjectiveFunc *SymbolicObjective_OneNode_ForwardAlg::NewObjectiveFunc (void)
{
	return new LnObjectiveFuncAdaptor(new ObjectiveFunc(*this,*symbolic,extraCm2HmmInfo.inequalitiesAndLocalVariables[cmNode].numLocalVariables,extraCm2HmmInfo.inequalitiesAndLocalVariables[cmNode].inequalityList));
}

SymbolicObjective_OneNode_ForwardAlg::ObjectiveFunc::ObjectiveFunc (SymbolicObjective_OneNode_ForwardAlg& _master,SymbolicProbVariableMath& symbolic_,int numLocalVariables_,const InequalityList& inequalityList_)
: master(_master)
, symbolic(symbolic_)
, inequalityList(inequalityList_)
{
	numLocalVariables=numLocalVariables_;
}
void SymbolicObjective_OneNode_ForwardAlg::ObjectiveFunc::Eval (double& f,std::vector<double>& gradient,vector2d<double>& hessian,const std::vector<double>& problemVars,bool calculateHessian,bool calculateGradient)
{
	symbolic.Eval(numLocalVariables,f,gradient,hessian,problemVars,calculateHessian);
}
int SymbolicObjective_OneNode_ForwardAlg::ObjectiveFunc::GetNumProblemVars (void)
{
	return numLocalVariables;
}
void SymbolicObjective_OneNode_ForwardAlg::ObjectiveFunc::LocalToProblemVars (std::vector<double>& problemVars,const std::vector<float>& localVars)
{
	VectorFloatToDouble(problemVars,localVars);
}
void SymbolicObjective_OneNode_ForwardAlg::ObjectiveFunc::ProblemToLocalVars (std::vector<float>& localVars,const std::vector<double>& problemVars)
{
	VectorDoubleToFloat(localVars,problemVars);
}
const InequalityList& SymbolicObjective_OneNode_ForwardAlg::ObjectiveFunc::GetInequalityList (void)
{
	return inequalityList;
}

SymbolicObjective_OneNode_ForwardAlg::ProblemDefinitionInstantiator::ProblemDefinitionInstantiator (void)
{
}
SymbolicObjective_OneNode_ForwardAlg::ProblemDefinitionInstantiator::~ProblemDefinitionInstantiator ()
{
}
ProblemDefinition *SymbolicObjective_OneNode_ForwardAlg::ProblemDefinitionInstantiator::NewProblemDefinition (const InfernalHmm& sourceHmm,const CovarianceModel& cm,int windowLen,CovarianceModel::Node cmNode,const ExtraCm2HmmInfo& extraCm2HmmInfo,char *fastaFileToSearch,int numLocalVariables,const std::vector<float>& localVars)
{
	return new SymbolicObjective_OneNode_ForwardAlg(sourceHmm,cm,windowLen,cmNode,extraCm2HmmInfo,fastaFileToSearch,numLocalVariables);
}


/////////////////////////
// GlobalForwardInfSymbolicObjectiveFunc

GlobalForwardInfSymbolicObjectiveFunc::GlobalForwardInfSymbolicObjectiveFunc (const InfernalHmm& sourceHmm,const CovarianceModel& cm,const ExtraCm2HmmInfo& extraCm2HmmInfo,int numAdjacentNodesToMerge)
{
	numGlobalVars=extraCm2HmmInfo.scoreVariablesInfo.numVariables;

	std::vector<int> identityMap;
	identityMap.resize(numGlobalVars);
	for (int i=0; i<numGlobalVars; i++) {
		identityMap[i]=i;
	}
	symbolic=new SymbolicProbVariableMath (identityMap,identityMap,cm,sourceHmm,extraCm2HmmInfo.scoreVariablesInfo);

	inequalityList.clear();
	for (CovarianceModel::Node cmNode=cm.GetFirstNode(); cmNode!=cm.GetLastNode(); ) {
		InequalityListList inequalityListList;

		int numNodesIncluded=0;
		while (cmNode!=cm.GetLastNode() && numNodesIncluded<numAdjacentNodesToMerge) {
			InequalityList dummy;
			inequalityListList.push_back(dummy);
			AddInequalitiesInTermsOfGlobalVars(inequalityListList.back(),extraCm2HmmInfo,cmNode);

			cmNode++;
			numNodesIncluded++;
		}

		AddCrossProductOfInequalityLists(inequalityList,inequalityListList);
	}
}
GlobalForwardInfSymbolicObjectiveFunc::~GlobalForwardInfSymbolicObjectiveFunc ()
{
	delete symbolic;
}
void GlobalForwardInfSymbolicObjectiveFunc::Eval (double& f,std::vector<double>& gradient,vector2d<double>& hessian,const std::vector<double>& problemVars,bool calculateHessian,bool calculateGradient)
{
	symbolic->Eval(numGlobalVars,f,gradient,hessian,problemVars,calculateHessian);
}
int GlobalForwardInfSymbolicObjectiveFunc::GetNumProblemVars (void)
{
	return numGlobalVars;
}
void GlobalForwardInfSymbolicObjectiveFunc::LocalToProblemVars (std::vector<double>& problemVars,const std::vector<float>& localVars)
{
	Die("GlobalForwardInfSymbolicObjectiveFunc::LocalToProblemVars is meaningless");
}
void GlobalForwardInfSymbolicObjectiveFunc::ProblemToLocalVars (std::vector<float>& localVars,const std::vector<double>& problemVars)
{
	Die("GlobalForwardInfSymbolicObjectiveFunc::LocalToProblemVars is meaningless");
}
const InequalityList& GlobalForwardInfSymbolicObjectiveFunc::GetInequalityList (void)
{
	return inequalityList;
}

//////////////////////////////
// NodeCombinerForwardInfSymbolicObjectiveFunc

NodeCombinerForwardInfSymbolicObjectiveFunc::NodeCombinerForwardInfSymbolicObjectiveFunc (InfernalHmm& sourceHmm,const CovarianceModel& cm,const ExtraCm2HmmInfo& extraCm2HmmInfo,int numAdjacentNodesToMerge,int maxNodesAtATime,const CovarianceModel::Node cmStartNode,bool applyLog2)
{
	numGlobalVars=extraCm2HmmInfo.scoreVariablesInfo.numVariables;

	globalToProblemVars.assign(numGlobalVars,-1);

	inequalityList.clear();
	for (CovarianceModel::Node cmNode=cmStartNode; cmNode<cm.GetLastNode() && cmNode<cmStartNode.PlusInt(maxNodesAtATime); ) {
		InequalityListList inequalityListList;

		int numNodesIncluded=0;
		while (cmNode!=cm.GetLastNode() && numNodesIncluded<numAdjacentNodesToMerge) {
			InequalityList dummy;
			inequalityListList.push_back(dummy);
			AddInequalitiesInTermsOfGlobalVars(inequalityListList.back(),extraCm2HmmInfo,cmNode);

			std::list<int> globalVarNums;
			if (extraCm2HmmInfo.inequalitiesAndLocalVariables[cmNode].numLocalVariables>0) {
				GetGlobalVariableNumsForNode (globalVarNums,cm,sourceHmm,cmNode,extraCm2HmmInfo);
			}
			for (std::list<int>::iterator i=globalVarNums.begin(); i!=globalVarNums.end(); i++) {
				//printf("%d\n",*i);
				globalToProblemVars[*i]=0;
			}

			cmNode++;
			numNodesIncluded++;
		}

		AddCrossProductOfInequalityLists(inequalityList,inequalityListList);
	}

	numProblemVars=0;
	for (int var=0; var<numGlobalVars; var++) {
		if (globalToProblemVars[var]!=-1) {
			globalToProblemVars[var]=numProblemVars;
			numProblemVars++;
		}
	}
	problemToGlobalVars.resize(numProblemVars);
	for (int var=0; var<numGlobalVars; var++) {
		if (globalToProblemVars[var]!=-1) {
			problemToGlobalVars[globalToProblemVars[var]]=var;
		}
	}

	ConvertInequalityListVarNums (inequalityList,globalToProblemVars);

	symbolic=new SymbolicProbVariableMath (problemToGlobalVars,globalToProblemVars,cm,sourceHmm,extraCm2HmmInfo.scoreVariablesInfo,applyLog2);
}
NodeCombinerForwardInfSymbolicObjectiveFunc::~NodeCombinerForwardInfSymbolicObjectiveFunc ()
{
	delete symbolic;
}
void NodeCombinerForwardInfSymbolicObjectiveFunc::Eval (double& f,std::vector<double>& gradient,vector2d<double>& hessian,const std::vector<double>& problemVars,bool calculateHessian,bool calculateGradient)
{
	symbolic->Eval(numProblemVars,f,gradient,hessian,problemVars,calculateHessian);
}
int NodeCombinerForwardInfSymbolicObjectiveFunc::GetNumProblemVars (void)
{
	return numProblemVars;
}
void NodeCombinerForwardInfSymbolicObjectiveFunc::LocalToProblemVars (std::vector<double>& problemVars,const std::vector<float>& localVars)
{
	Die("NodeCombinerForwardInfSymbolicObjectiveFunc::LocalToProblemVars is meaningless");
}
void NodeCombinerForwardInfSymbolicObjectiveFunc::ProblemToLocalVars (std::vector<float>& localVars,const std::vector<double>& problemVars)
{
	Die("NodeCombinerForwardInfSymbolicObjectiveFunc::ProblemToLocalVars is meaningless");
}
const InequalityList& NodeCombinerForwardInfSymbolicObjectiveFunc::GetInequalityList (void)
{
	return inequalityList;
}
void NodeCombinerForwardInfSymbolicObjectiveFunc::GlobalToProblemVars (std::vector<double>& problemVars,const std::vector<double>& globalVars)
{
	problemVars.resize(numProblemVars);
	assert(numGlobalVars==(int)(globalVars.size()));
	for (size_t i=0; i<globalVars.size(); i++) {
		if (globalToProblemVars[i]!=-1) {
			problemVars[globalToProblemVars[i]]=globalVars[i];
		}
	}
}
void NodeCombinerForwardInfSymbolicObjectiveFunc::UpdateGlobalVarsFromProblemVars (std::vector<double>& globalVars,const std::vector<double>& problemVars)
{
	assert(numProblemVars==(int)(problemVars.size()));
	assert(numGlobalVars==(int)(globalVars.size()));

	for (size_t i=0; i<globalVars.size(); i++) {
		if (globalToProblemVars[i]!=-1) {
			globalVars[i]=problemVars[globalToProblemVars[i]];
		}
	}
}


float ComputeLhs (const std::vector<float>& localVariablesToValue,const Inequality& ineq)
{
	float lhs=0;
	for (std::list<InequalityTerm>::const_iterator termIter=ineq.lhs.begin(); termIter!=ineq.lhs.end(); termIter++) {
		int localVar=termIter->variableNum;
		lhs += localVariablesToValue[localVar];
	}
	return lhs;
}

void DumpHmmInflationsPerInequality (char *cmFileName,const std::string& programParams,bool doLocalAlignment)
{
	CovarianceModel cm;
	cm.Load(cmFileName,doLocalAlignment);

	InfernalHmm infernalHmm;
	Cm2Hmm (infernalHmm,HmmBuildType_Original,cm,cmFileName,programParams);

	// get inequalities
	ExtraCm2HmmInfo extraCm2HmmInfo;
	extraCm2HmmInfo.actuallySolveScores=false;
	InfernalHmm dummyInfernalHmm;
	Cm2Hmm_WithWeighting_NoCaching (dummyInfernalHmm,HmmBuildType_Original,cm,cmFileName,NULL,&extraCm2HmmInfo);

	float minInflation=+FLT_MAX;

	for (CovarianceModel::Node cmNode=cm.GetFirstNode(); cmNode!=cm.GetLastNode(); cmNode++) {

		if (extraCm2HmmInfo.inequalitiesAndLocalVariables[cmNode].numLocalVariables>0) {

			// get scores' curr values
			std::vector<float> localVariablesToValue;
			GetLocalVariablesValueForNode(localVariablesToValue,cm,infernalHmm,cmNode,extraCm2HmmInfo);

			const InequalityList& inequalityList=extraCm2HmmInfo.inequalitiesAndLocalVariables[cmNode].inequalityList;
			InequalityList::const_iterator ineqIter;
			for (ineqIter=inequalityList.begin(); ineqIter!=inequalityList.end(); ineqIter++) {

				float lhs=ComputeLhs(localVariablesToValue,*ineqIter);
				float thisInflation=lhs - ineqIter->rhs;
				minInflation=std::min(minInflation,thisInflation);

				printf("node #%d,ineq,%f\n",CovarianceModel::NodeToInt(cmNode),thisInflation);
			}
		}
	}

	printf("\n\nlowest inflation = %f (there's a problem if this is much below -1e-6)\n",minInflation);
}

void MakeCmHmmForOptimization (CovarianceModel& cm,InfernalHmm& infernalHmm,ExtraCm2HmmInfo& extraCm2HmmInfo,char *cmFileName,bool doLocalAlignment,const std::string& programParams,Cm2Hmm_HmmBuildType hmmType)
{
	cm.Load(cmFileName,doLocalAlignment);

	Cm2Hmm (infernalHmm,hmmType,cm,cmFileName,programParams);
	infernalHmm.BuildNonSavedInfoIfNecessary();
	infernalHmm.AddBuildDescription(programParams);
	infernalHmm.SetHmmBuildType(hmmType);

	// get inequalities
	extraCm2HmmInfo.actuallySolveScores=false;
	InfernalHmm dummyInfernalHmm;
	Cm2Hmm_WithWeighting_NoCaching (dummyInfernalHmm,hmmType,cm,cmFileName,NULL,&extraCm2HmmInfo);
}

void GlobalHmmOptimizer (char *cmFileName,bool doLocalAlignment,const std::string& programParams,int numAdjacentNodesToMerge)
{
	CovarianceModel cm;
	InfernalHmm infernalHmm;
	ExtraCm2HmmInfo extraCm2HmmInfo;
	MakeCmHmmForOptimization (cm,infernalHmm,extraCm2HmmInfo,cmFileName,doLocalAlignment,programParams);

	int B=0;
	int C=1;
	SolverWrapper *solverWrapper=NewSolverWrapper_cfsqp (B,C);

	std::vector<double> globalVars;
	GetGlobalVarsFromInfernalHmm (globalVars,infernalHmm,extraCm2HmmInfo.scoreVariablesInfo.globalVariableToTransitionOrEmissionVector);

	printf("merging %d adjacent nodes\n",numAdjacentNodesToMerge);
	GlobalForwardInfSymbolicObjectiveFunc *globalObjectiveFunc=new GlobalForwardInfSymbolicObjectiveFunc(infernalHmm,cm,extraCm2HmmInfo,numAdjacentNodesToMerge);
	LnObjectiveFuncAdaptor lnObjectiveFunc(globalObjectiveFunc);
	std::vector<double> gradient;
	vector2d<double> hessian;
	double currLogVal;
	lnObjectiveFunc.Eval(currLogVal,gradient,hessian,globalVars,false,false);

	printf("NOTE: this function won't actually save the optimized HMM, since global optimization doesn't seem to perform any better than node-at-a-time optimization; it doesn't run faster (at least not for reasonably-sized HMMs), and doesn't come up with a better HMM.\n");

	printf("starting log val = %lg  (log_2 = %lg)\n",currLogVal,currLogVal/log(2.0));

	solverWrapper->Solve(&lnObjectiveFunc,globalVars,100.0);

	delete solverWrapper;
}

void HmmOptimizer_NodeCombiner (char *cmFileName,bool doLocalAlignment,const std::string& programParams,int numAdjacentNodesToMerge,int maxNodesAtATime,int numIters,Cm2Hmm_HmmBuildType hmmType,SolverWrapper *solverWrapper,bool saveHmmsInProgress,const char *hmmSaveFileName,HmmFileFormat hmmFileFormat)
{
	CovarianceModel cm;
	InfernalHmm infernalHmm;
	ExtraCm2HmmInfo extraCm2HmmInfo;
	MakeCmHmmForOptimization (cm,infernalHmm,extraCm2HmmInfo,cmFileName,doLocalAlignment,programParams,hmmType);

	printf("Optimizing HMM using infinite-length forward algorithm\n");
	printf("merging %d adjacent nodes, %d nodes at a time,\n",numAdjacentNodesToMerge,maxNodesAtATime);

	CovarianceModel::Node cmStartNode=cm.GetFirstNode();
	double lastScoreAtFirstNode=+DBL_MAX;
	for (int iter=0; iter<numIters; iter++) {

		std::vector<double> globalVars;
		GetGlobalVarsFromInfernalHmm (globalVars,infernalHmm,extraCm2HmmInfo.scoreVariablesInfo.globalVariableToTransitionOrEmissionVector);

		NodeCombinerForwardInfSymbolicObjectiveFunc *objectiveFunc=new NodeCombinerForwardInfSymbolicObjectiveFunc(infernalHmm,cm,extraCm2HmmInfo,numAdjacentNodesToMerge,maxNodesAtATime,cmStartNode,true);
		//LnObjectiveFuncAdaptor lnObjectiveFunc(objectiveFunc);

		if (objectiveFunc->GetNumProblemVars()==0) {

			printf("iter #%d (start node=%d):  skipping, since no vars\n",iter,InfernalHmm::NodeToInt(cmStartNode));
		}
		else {

			std::vector<double> problemVars;
			objectiveFunc->GlobalToProblemVars(problemVars,globalVars);

			std::vector<double> gradient;
			vector2d<double> hessian;
			double currLogVal;
			objectiveFunc->Eval(currLogVal,gradient,hessian,problemVars,false,false);
			currLogVal *= log(2.0);

			printf("iter #%d (start node=%d, # vars=%d):  starting log val = %lg  (log_2 = %lg)\n",iter,InfernalHmm::NodeToInt(cmStartNode),objectiveFunc->GetNumProblemVars(),currLogVal,currLogVal/log(2.0));

			if (cmStartNode==cm.GetFirstNode()) {
				if (lastScoreAtFirstNode!=+DBL_MAX) {
					printf("At first node.  Prev score @ first node = %f .   Curr score = %f .\n",lastScoreAtFirstNode,currLogVal);
				}
				if (fabs(lastScoreAtFirstNode-currLogVal)<1e-5) {
					printf("We haven't made enough improvement since last run thru the nodes.  Stopping.\n");
					break;
				}
				lastScoreAtFirstNode=currLogVal;

				// seems a good time to clear any unused stuff out of the cache
				SymbolicMath::Expression::ClearConstCache();
			}

			problemVars=solverWrapper->Solve(objectiveFunc,problemVars,100.0);

			objectiveFunc->UpdateGlobalVarsFromProblemVars(globalVars,problemVars);
			SetGlobalVarsIntoInfernalHmm(infernalHmm,globalVars,extraCm2HmmInfo.scoreVariablesInfo.globalVariableToTransitionOrEmissionVector);
		}

		cmStartNode += maxNodesAtATime;
		if (cmStartNode>=cm.GetLastNode()) {
			cmStartNode=cm.GetFirstNode();
		}

		if (saveHmmsInProgress) {
			infernalHmm.SaveInFormat(hmmSaveFileName,hmmFileFormat);
		}

		delete objectiveFunc;
	}
	infernalHmm.SaveInFormat(hmmSaveFileName,hmmFileFormat);

	delete solverWrapper;
}

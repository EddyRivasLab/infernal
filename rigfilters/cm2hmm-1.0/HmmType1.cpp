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
#include <UseDebugNew.h>
#include "cmzasha.h"


///////////////////////
// HmmType1

void HmmType1::Init (const InfernalHmm& unreversedSourceHmm)
{
	InfernalHmm sourceHmm;
	sourceHmm.CopyReverseOf(unreversedSourceHmm);

	numStates=sourceHmm.GetNumStates();

	stateInfoVector.resize(numStates);
	singletEmissionScores.resize(Alphabet_size);
	for (int nuc=0; nuc<Alphabet_size; nuc++) {
		singletEmissionScores[nuc].resize(numStates);
	}

	doLocal=false;
	if (sourceHmm.DoLocal()) {
		doLocal=true;
		localStateInfoVector.resize(numStates);
	}

	hmm2CmStateVector=sourceHmm.GetCmStateVector();
	cm2HmmStateVector=sourceHmm.GetHmmStateVector();

	CovarianceModel::State state;
	for (state=sourceHmm.GetFirstState(); state!=sourceHmm.GetLastState(); state++) {
		int i_state=CovarianceModel::StateToInt(state);

		int stateType=sourceHmm.GetStateType(state);
		assert(stateType==ML_st || stateType==IL_st || stateType==E_st || stateType==D_st || stateType==PASSTHRU_st || stateType==S_st); // these are the only types that HMM should have
		if (stateType==ML_st || stateType==IL_st) {
			stateInfoVector[i_state].isEmitting=true;
		}
		else {
			if (stateType==E_st || stateType==D_st || stateType==S_st || stateType==PASSTHRU_st) {
				stateInfoVector[i_state].isEmitting=false;
			}
			else {
				Die("Internal error: Unexpected state type in HMM represented as a CovarianceModel");
			}
		}

		for (int nuc=0; nuc<Alphabet_size; nuc++) {
			float floatScore;
			if (stateInfoVector[i_state].isEmitting) {
				floatScore=sourceHmm.GetSingletEmissionScore(state,nuc);
			}
			else {
				floatScore=0;
			}
			singletEmissionScores[nuc][i_state]=ScoreHelperType::FloatToScoreType(floatScore);
		}

		stateInfoVector[i_state].numChildren=sourceHmm.GetNumChildren(state);
		assert(stateInfoVector[i_state].numChildren<=MAX_CHILDREN);
		stateInfoVector[i_state].firstChild=-1;
		stateInfoVector[i_state].isRightState=sourceHmm.IsRightState(state);

		if (stateInfoVector[i_state].numChildren>0) {
			stateInfoVector[i_state].firstChild=CovarianceModel::StateToInt(sourceHmm.GetNthChildState(state,0));

			for (int childNum=0; childNum<sourceHmm.GetNumChildren(state); childNum++) {
				if (GetNthChildState(i_state,childNum)!=CovarianceModel::StateToInt(sourceHmm.GetNthChildState(state,childNum))) {
					Die("Internal error: CovarianceModel apparently doesn't have consecutively-numbered states any more.  I hate it when things change on me, so I'm going to die now.");
				}
				stateInfoVector[i_state].tsc[childNum]=ScoreHelperType::FloatToScoreType(sourceHmm.GetNthChildTsc(state,childNum));
			}
		}

		if (sourceHmm.DoLocal()) {

			localStateInfoVector[i_state].leftwardBeginsc=sourceHmm.GetLeftwardBeginsc(state);
			localStateInfoVector[i_state].rightwardBeginsc=sourceHmm.GetRightwardBeginsc(state);

			int numLinks=sourceHmm.GetNumEndscLinksToLeft(state);
			localStateInfoVector[i_state].endscLinksToLeft.resize(numLinks);
			for (int link=0; link<numLinks; link++) {
				localStateInfoVector[i_state].endscLinksToLeft[link].hmmLeftState=InfernalHmm::StateToInt(sourceHmm.GetEndscLinkToLeft_State(state,link));
				localStateInfoVector[i_state].endscLinksToLeft[link].endsc=sourceHmm.GetEndscLinkToLeft_Endsc(state,link);
			}
		}
	}

	// Forward Alg info
	forwardAlgInfoVector.resize(GetNumStates());
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		for (int i=0; i<GetNumChildren(state); i++) {
			forwardAlgInfoVector[state].transitionProbs[i]=(float)(pow2(GetNthChildTsc(state,i)));
		}
		assert(Alphabet_size==MAXABET); // I used MAXABET to define the bounds of singletEmissionProbs in ScanHMM.h
		for (int nuc=0; nuc<Alphabet_size; nuc++) {
			if (IsEmittingState(state)) {
				forwardAlgInfoVector[state].singletEmissionProbs[nuc]=(float)(pow2(GetSingletEmissionScore(state,nuc)));
			}
		}
	}
}
void HmmType1::Dump (FILE *out) const
{
	fprintf(out,"# states: %d\n",GetNumStates());
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		fprintf(out,"state #%d : %s\n",state,IsEmittingState(state) ? "EMITTING" : "silent");
		if (IsEmittingState(state)) {
			fprintf(out,"\t");
			for (int nuc=0; nuc<MAXABET; nuc++) {
				fprintf(out,"%c: %f , ",nucs[nuc],GetSingletEmissionScore(state,nuc));
			}
			fprintf(out,"\n");
		}
		fprintf(out,"\t");
		for (int child=0; child<GetNumChildren(state); child++) {
			fprintf(out,"#%d: %f , ",GetNthChildState(state,child),GetNthChildTsc(state,child));
		}
		fprintf(out,"\n");
	}
}
void HmmType1::Dump (const char *fileName) const
{
	FILE *out=ThrowingFopen(fileName,"wt");
	Dump(out);
	fclose(out);
}



////////////////////
// TransitionCounter

TransitionCounter::TransitionCounter (const HmmType1& _hmm)
: hmm(_hmm)
{
	transitionCounts.resize(hmm.GetNumStates());
	emitCounts.resize(hmm.GetNumStates());

	HmmType1::State state;
	for (state=hmm.GetFirstState(); state!=hmm.GetLastState(); state++) {
		transitionCounts[state].assign(hmm.GetNumChildren(state),0);
		if (hmm.IsEmittingState(state)) {
			emitCounts[state].assign(Alphabet_size,0);
		}
	}

	numSamples=0;
}
TransitionCounter::~TransitionCounter ()
{
}
void TransitionCounter::AddSample (void)
{
	numSamples++;
}
int TransitionCounter::GetNumSamples (void) const
{
	return numSamples;
}
int TransitionCounter::GetTransitionCount_Unreversed (int unreversed_fromState,int unreversed_toState) const
{
	int fromState=unreversed_toState;
	int toState=unreversed_fromState;
	return transitionCounts[fromState][hmm.GetChildNum_Slow(fromState,toState)];
}
int TransitionCounter::GetTransitionCount_Unreversed (CovarianceModel::State unreversed_fromState,CovarianceModel::State unreversed_toState) const
{
	return GetTransitionCount_Unreversed(CovarianceModel::StateToInt(unreversed_fromState),CovarianceModel::StateToInt(unreversed_toState));
}
int TransitionCounter::GetEmitCount (int state,int nuc) const
{
	return emitCounts[state][nuc];
}
int TransitionCounter::GetEmitCount (CovarianceModel::State state,int nuc) const
{
	return GetEmitCount(CovarianceModel::StateToInt(state),nuc);
}
double TransitionCounter::GetEntryProbability_Unreversed (int unreversed_toState,int unreversed_alternateToStateFirst,int unreversed_alternateToStateLast) const
{
	assert(unreversed_alternateToStateLast>=unreversed_alternateToStateFirst);
	assert(unreversed_toState>=unreversed_alternateToStateFirst && unreversed_toState<unreversed_alternateToStateLast);

	int fromState=unreversed_toState;
	int alternateFromStateFirst=unreversed_alternateToStateFirst;
	int alternateFromStateLast=unreversed_alternateToStateLast;

	int total=0;
	int numerator=-1;
	for (int currFromState=alternateFromStateFirst; currFromState<alternateFromStateLast; currFromState++) {

		int totalForState=0;
		for (int childNum=0; childNum<hmm.GetNumChildren(currFromState); childNum++) {
			if (currFromState!=hmm.GetNthChildState(currFromState,childNum)) {
				totalForState += transitionCounts[currFromState][childNum];
			}
		}
		if (currFromState==fromState) {
			numerator=totalForState;
		}

		total += totalForState;
	}
	assert(numerator>=0);

	if (total==0) {
		// special case for start node
		assert(numerator==0);
		return 0; // path is never taken at all, so prob is 0
	}
	else {
		double result=(double)(numerator)/(double)(total);
		assert(result>=0 && result<=1);
		return result;
	}
} 
double TransitionCounter::GetEntryProbability_Unreversed (CovarianceModel::State unreversed_toState,CovarianceModel::State unreversed_alternateToStateFirst,CovarianceModel::State unreversed_alternateToStateLast) const
{
	return GetEntryProbability_Unreversed(CovarianceModel::StateToInt(unreversed_toState),CovarianceModel::StateToInt(unreversed_alternateToStateFirst),CovarianceModel::StateToInt(unreversed_alternateToStateLast));
}
double TransitionCounter::GetTransitionFrequency_Unreversed (int unreversed_fromState,int unreversed_toState,double pseudocount) const
{
	int fromState=unreversed_toState;
	int toState=unreversed_fromState;

	if (fromState==toState) {
		assert(false); // we shouldn't need to explore self-loops
		return 1;
	}

	double prChild=-1; // we'll actually set this later
	double prState=0;
	for (size_t i=0; i<transitionCounts[fromState].size(); i++) {
		HmmType1::State childState=hmm.GetNthChildState(fromState,(int)i);
		if (childState!=fromState) { // avoid self-loops
			prState += (double)(transitionCounts[fromState][i]) + pseudocount;
			if (childState==toState) {
				assert(prChild==-1); // should only set this once
				prChild=(double)(transitionCounts[fromState][i]) + pseudocount;
			}
		}
	}
	assert(prChild!=-1); // should set this

	if (prState==0) {
		assert(prChild==0);
		return 0;
	}
	else {
		double result=prChild/prState;
		assert(result>=0 && result<=1);
		return result;
	}
}
double TransitionCounter::GetTransitionFrequency_Unreversed (CovarianceModel::State unreversed_fromState,CovarianceModel::State unreversed_toState,double pseudocount) const
{
	return GetTransitionFrequency_Unreversed (CovarianceModel::StateToInt(unreversed_fromState),CovarianceModel::StateToInt(unreversed_toState),pseudocount);
}
double TransitionCounter::GetEmitFrequency (int state,int nuc,double pseudocount) const
{
	double prNuc=(double)(emitCounts[state][nuc]) + pseudocount;

	double prState=0;
	for (size_t i=0; i<emitCounts[state].size(); i++) {
		prState += (double)(emitCounts[state][i]) + pseudocount;
	}

	if (prState==0) {
		assert(prNuc==0);
		return 0;
	}
	else {
		double result=prNuc/prState;
		assert(result>=0 && result<=1);
		return result;
	}
}
double TransitionCounter::GetEmitFrequency (CovarianceModel::State state,int nuc,double pseudocount) const
{
	return GetEmitFrequency(CovarianceModel::StateToInt(state),nuc,pseudocount);
}
void TransitionCounter::AddTransition (int state,int childNum)
{
	transitionCounts[state][childNum]++;
}
void TransitionCounter::AddEmit (int state,int nuc)
{
	emitCounts[state][nuc]++;
}
void TransitionCounter::Dump (FILE *out)
{
	HmmType1::State state,childState;
	for (state=hmm.GetFirstState(); state!=hmm.GetLastState(); state++) {
		fprintf(out,"state #%d:\n",state);
		for (int childNum=0; childNum<hmm.GetNumChildren(state); childNum++) {
			childState=hmm.GetNthChildState(state,childNum);
			fprintf(out,"\t--> state #%d: %d\n",childState,transitionCounts[state][childNum]);
		}
		if (hmm.IsEmittingState(state)) {
			for (int nuc=0; nuc<Alphabet_size; nuc++) {
				fprintf(out,"\t\temit %c: %d\n",nucs[nuc],emitCounts[state][nuc]);
			}
		}
	}
}

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


// helper function for InfiniteLengthForwardAlg_LastStateDoesntEmit (see below)
template <class Hmm,class Real>
Real CalcExpectedEmitProb_0order (const Hmm& hmm,const typename Hmm::State state,MarkovModelStats& markovModelStats)
{
	if (hmm.IsEmittingState(state)) {
		Real emitProb=0.0;
		for (int nuc=0; nuc<Alphabet_size; nuc++) {
			emitProb += markovModelStats.GetProbOfNuc_0order(nuc) * hmm.GetSingletEmissionProb(state,nuc);
		}
		return emitProb;
	}
	else {
		return 1;
	}
}

// helper funcs for below
double RealToConcreteDouble (double selfLoopProb)
{
	return selfLoopProb;
}
double RealToConcreteDouble (SymbolicProbVariableMath::Expression selfLoopProb)
{
	return selfLoopProb.ToConstDouble();
}


// returns the sum of the probabilities of paths of any length in the HMM starting at startState and ending at endState.  If startState or endState emits, the emission probability is counted.  If you don't want their emission probs, you'll have to divide it out.  It is illegal for startState or endState to be of type IL_st (really, it's illegal for them to self loop).
// idea: use a dyn prog table that has only 1 dim, corresponding to states.  The dyn prog stores the probability of generating a path that ends at that state, including any emission(s) done at this state (could be more than 1 emission for insert states).  For last state, divide out emitScore, if any.  To simplify, we disallow start/end states being insert states; this is fine because for the node-probs, we don't have to stitch together paths at insert states.
template <class Hmm,class Real>
Real /* expected score */ InfiniteLengthForwardAlg_LastStateDoesntEmit (const Hmm& hmm,const typename Hmm::State startState,const typename Hmm::State endState,MarkovModelStats& markovModelStats)
{
	assert(markovModelStats.GetOrder()==0);
	if (markovModelStats.GetOrder()!=0) {
		throw SimpleStringException("%s:%d",__FILE__,__LINE__);
	}

	vector<Real> table;
	table.assign(hmm.GetNumStates(),0.0); // clear to 0, so we don't have to do any explicit initialization

	typename Hmm::State state=startState;
	table[state]=CalcExpectedEmitProb_0order<Hmm,Real>(hmm,startState,markovModelStats); // by the definition of the table, it must include emissions
	while (1) {
		state++;

		if (state>endState) {
			break;
		}

		Real emitProb=CalcExpectedEmitProb_0order<Hmm,Real>(hmm,state,markovModelStats);

		Real totalTransitionProb=0.0;
		Real selfLoopMultiplier=1.0;
		for (int childNum=0; childNum<hmm.GetNumChildren(state); childNum++) {
			const typename Hmm::State childState=hmm.GetNthChildState(state,childNum);
			const Real tsc(hmm.GetNthChildTransitionProb(state,childNum));
			if (childState==state) {
				// this is a self loop
				// compute the prob, which is a sum of a geometric series.  This prob will be multiplied into other things
				Real selfLoopProb_Real=emitProb*tsc;
				double selfLoopProb=RealToConcreteDouble (selfLoopProb_Real);
				if (selfLoopProb>=1) {
					assert(false);
					throw SimpleStringException("Self-loop (presumably in insert state) has >1 probability, which will converge to infinity, and make it impossible to apply the forward alg with infinite-length strings.  One non-hacky solution: model the probabilities of the input strings s.t. longer strings are less likely, but this requires a big re-think of the infinite-length forward alg.  A quick-fix solution is to adjust the 0-order Markov model you're using (with the --set-input-mm-file flag) such that the insert state is okay.  If you set the input Markov model to be exactly the same as the null (background) model used for the CM, then insert self-loops can never have >1 probability.  Hmm state #%d",state);
				}
				selfLoopMultiplier=1.0/(1.0-selfLoopProb); // sum of infinite geometric series
			}
			else {
				totalTransitionProb += tsc*table[childState];
			}
		}

		// justification of this formula for insert state (I think the tricky case): probability of being at one of our child states, times probability of transitioning here, which is given by totalTransitionProb.  This is our probability of getting into the state at all.  At this point, we have to do an emission (if we're an insert state), so multiply in emitProb.  Now, if we're an insert state, we could just go to the next state, with whatever probability that has (but that's taken care of by the other state).  Or we could loop, with probability emitProb*self_loop_tsc.  Since we can loop an unbounded amount of times, we should sum the infinite geometric series, as done for selfLoopMultiplier.
		Real thisProb=emitProb*totalTransitionProb*selfLoopMultiplier;
		table[state]=thisProb;
	}

	return table[endState];
}
template <class Hmm,class Real>
Real /* expected score */ InfiniteLengthForwardAlg_LastStateDoesntEmit_StartWithInfernalHmm (const InfernalHmm& infernalHmm,const InfernalHmm::State infernalStartState,const InfernalHmm::State infernalEndState,MarkovModelStats& markovModelStats)
{
	assert(infernalHmm.GetStateType(infernalStartState)!=IL_st && infernalHmm.GetStateType(infernalEndState)!=IL_st);

	Hmm hmm;
	hmm.Init(infernalHmm);
	const typename Hmm::State startState=InfernalHmm::StateToInt(infernalStartState);
	const typename Hmm::State endState=InfernalHmm::StateToInt(infernalEndState);

	return InfiniteLengthForwardAlg_LastStateDoesntEmit<Hmm,Real>(hmm,startState,endState,markovModelStats);
}

// I didn't really implement the following function.  It looks like higher-order Markov models aren't going to help the HMMs all that much, and it's tricky to implement.
// generalizes InfiniteLengthForwardAlg_LastStateDoesntEmit to arbitrary-order Markov models.  I've re-written the function, since much has to change.  Also, for this function, I'm planning to generate the objective functions to optimize HMMs using a fake 'Real' class that actually does symbolic expressions.  Therefore, there's no need to be able to set a start/end state, or make sure that we don't over/under-count emissions (in InfiniteLengthForwardAlg_LastStateDoesntEmit, I had to make sure it was as if endState didn't emit).  This is really convenient, since emissions will be tricky here.
// issues: IL self-loops with a general MM -- the infinite geometric series could get weird...
template <class Hmm,class Real>
Real /* expected score */ InfiniteLengthForwardAlg_GeneralMarkov (const Hmm& hmm,MarkovModelStats& markovModelStats)
{
	typename Hmm::State state;
	const typename Hmm::State startState=hmm.GetStartState();
	const typename Hmm::State endState=hmm.GetEndState();

	int order=markovModelStats.GetOrder();

	// table is first indexed by a state.  Within the state, it is indexed by the Markov-model context, with 4^N entries for an N-order model.  Within each of these entries is a Real.  The entry for state V and context C is Score(V^C) the total score for paths ending in state V that produce strings ending in C.
	vector<VariableDimVector<Real> > table;
	table.resize(hmm.GetNumStates());
	// initialize each state-cell
	for (state=startState; state<=endState; state++) {
		table[state].resize(order,Alphabet_size);
	}

	// initialize the start state.  The total probability of the start state sums to 1.  However, it's marginalized on the probability of each possible context of size N (for N-order model)
	VariableDimVector<double> contextDistribution; // type 'double', since it's constants
	markovModelStats.GetContextDistribution(contextDistribution);
	// since contextDistribution is of type 'double', and we want to set it into type 'Real' (which might be abstract), we must do it manually
	NaryCounter counter(order,Alphabet_size);
	bool counting=true;
	while (counting) {
		table[startState].GetRef(counter) = contextDistribution.Get(counter);
		counting=counter.Next();
	}

	// temp variable, declared here to avoid reallocation
	vector<int> childContext;
	childContext.reserve(order+1);

	// now start computing values of states
	state=startState;
	while (1) {
		state++;

		if (state>endState) {
			break;
		}

		// check if the state self-loops.  if so, we'll have to do the infinite geometric series, only this time we need to solve a system of linear equations to find the limit
		// the trick: for 0-order model case, we need to solve the usual sum of infinite geometric series.  Here's the usual solution.  Let S be the sum we wish to solve, i.e. S=1+x+x^2+x^3+x^4+... .  Then, S=1+x*S.  Therefore, S=1/(1-x)
		// now, here's the 1-order case.  We can write out
		// S_a = 1 + p_aa*S_a + p_ac*S_c + p_ag*S_g + p_au*S_u
		// S_c = 1 + ...
		// ...
		//
		// we then solve this system for S_{context}, and each S_{context} is our self-loop multiplier.
		// note that in this case, p_xy includes the Markov model probability Pr(x|y), the self-loop tsc, and the HMM's emission score.
		assert(false); // not implemented

		// first loop over contexts for curr state
		counter.Init();
		counting=true;
		while (counting) {

			// now loop over children -- we need to do this now, since the emission probs are dependant on which child (it's simpler if the state doesn't emit, but no point in dealing with that now)
			Real probOfThisContext=0;
			for (int childNum=0; childNum<hmm.GetNumChildren(state); childNum++) {
				const typename Hmm::State childState=hmm.GetNthChildState(state,childNum);
				const Real tsc=hmm.GetNthChildTransitionProb(state,childNum);

				if (childState==state) {
					assert(false); // implement logic to work out self-loop prob... well, use the self-prob that was calculated above
				}
				else {

					Real emitProb;
					if (hmm.IsEmittingState(state)) {

						// handle case with emission
						childContext.clear();
						childContext.push_back(0); // dummy to use the cell -- this is the nucleotide that goes out of the context when the current state emits its nucleotide
						childContext.insert(childContext.end(),counter.begin(),counter.end());
						int emittedNuc=childContext.back(); // this is the nucleotide the current state emits
						childContext.pop_back();

						emitProb=0;
						for (int outOfContextNuc=0; outOfContextNuc<Alphabet_size; outOfContextNuc++) {
							childContext[0]=outOfContextNuc;
							// Pr(child at this context)*Pr(emittedNuc | childContext)*Score(emittedNuc | state)
							emitProb += table[childState].Get(childContext)*markovModelStats.GetProbOfNuc(emittedNuc,childContext) * hmm.GetSingletEmissionProb(state,emittedNuc);
						}
					}
					else {
						emitProb=table[childState].Get(counter);
					}
					probOfThisContext += tsc*emitProb;
				}

				table[state].GetRef(counter) = probOfThisContext*selfLoopSomething;
			}

			counting=counter.Next();
		}
	}

	// finally, add up the marginalized probabilities in the last state
	Real total=0;
	counter.Init();
	counting=true;
	while (counting) {

		total += table[endState].Get(counter);
		counting=counter.Next();
	}

	return total;
}

double CalcExpectedEmitProb_0order_InfernalHmmDouble (const InfernalHmm& infernalHmm,InfernalHmm::State state,MarkovModelStats& markovModelStats)
{
	return CalcExpectedEmitProb_0order<InfernalHmm,double>(infernalHmm,state,markovModelStats);
}
double InfiniteLengthForwardAlg_LastStateDoesntEmit_HmmType1Double (const InfernalHmm& infernalHmm,const InfernalHmm::State infernalStartState,const InfernalHmm::State infernalEndState,MarkovModelStats& markovModelStats)
{
	return InfiniteLengthForwardAlg_LastStateDoesntEmit_StartWithInfernalHmm<HmmType1,double>(infernalHmm,infernalStartState,infernalEndState,markovModelStats);
}
double InfiniteLengthForwardAlg_LastStateDoesntEmit_HmmType1Double (const InfernalHmm& infernalHmm,MarkovModelStats& markovModelStats)
{
	return InfiniteLengthForwardAlg_LastStateDoesntEmit_HmmType1Double(infernalHmm,infernalHmm.GetFirstState(),infernalHmm.GetActualLastState(),markovModelStats);
}
SymbolicProbVariableMath::Expression InfiniteLengthForwardAlg_Symbolic (const SymbolicProbVariableMath::Hmm& hmm,MarkovModelStats& markovModelStats)
{
	return InfiniteLengthForwardAlg_LastStateDoesntEmit<SymbolicProbVariableMath::Hmm,SymbolicMath::Expression>(hmm,hmm.GetFirstState(),hmm.GetActualLastState(),markovModelStats);
}
double InfiniteLengthForwardAlg_LastStateDoesntEmit_HmmType1Double (const HmmType1& hmm,MarkovModelStats& markovModelStats)
{
	return InfiniteLengthForwardAlg_LastStateDoesntEmit<HmmType1,double>(hmm,hmm.GetFirstState(),hmm.GetActualLastState(),markovModelStats);
}

#ifndef CM2HMM_ONLY
SymbolicProbVariableMath::Expression InfiniteLengthForwardAlg_Symbolic (const SymbolicMath_HmmWithPenalties::Hmm& hmm,MarkovModelStats& markovModelStats)
{
	return InfiniteLengthForwardAlg_LastStateDoesntEmit<SymbolicMath_HmmWithPenalties::Hmm,SymbolicMath::Expression>(hmm,hmm.GetFirstState(),hmm.GetActualLastState(),markovModelStats);
}
double InfiniteLengthForwardAlg_LastStateDoesntEmit_HmmWithPenaltiesDouble (const HmmWithPenalties& hmm,MarkovModelStats& markovModelStats)
{
	return InfiniteLengthForwardAlg_LastStateDoesntEmit<HmmWithPenalties,double>(hmm,hmm.GetFirstState(),hmm.GetActualLastState(),markovModelStats);
}
#endif

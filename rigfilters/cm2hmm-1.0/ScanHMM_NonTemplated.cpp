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

/*
This has related functionality to ScanHMM.cpp, but it's un-templated code, which I believe may run faster.
*/

void ScanHmm_HmmType1Float_NonTemplated_Window(HitList& hmmHitList,int startPosToScan,int endPosToScan,const HmmType1& hmm,float minLodScoreForHit,CykscanStats& cykscanStats,const char *rnaSequence, int windowLen,float *prevTable,float *currTable)
{
	// init table
	HmmType1::State state=hmm.GetFirstState();
	prevTable[state]=0;
	while (1) {

		state++;
		if (state==hmm.GetLastState()) {
			break;
		}

		if (hmm.IsEmittingState(state)) {
			// during initialization, nothing can emit
			prevTable[state]=(float)IMPOSSIBLE;
		}
		else {
			float bestTransitionScore=(float)IMPOSSIBLE;
			for (int child=0; child<hmm.GetNumChildren(state); child++) {
				const float thisTransitionScore=prevTable[hmm.GetNthChildState(state,child)] + hmm.GetNthChildTsc(state,child);
				bestTransitionScore=std::max(bestTransitionScore,thisTransitionScore);
			}
			prevTable[state]=bestTransitionScore;
		}

		if (cykscanStats.collectHmmFullDynProgTable) {
			const float templatedHmmScore=cykscanStats.fullHmmDynProgTable[0][state];
			if (fabs(templatedHmmScore - prevTable[state])>2e-6) {
				assert(false);
				throw SimpleStringException("Sanity check failed: HMM scores don't match.  I am insane.  Initializing table.  templatedHmmScore=%g, prevTable[state]=%g.  state=%d",templatedHmmScore,prevTable[state],state);
			}
		}
	}

	// scan
	for (int windowLast=startPosToScan+1; windowLast<=endPosToScan; windowLast++) {
		const char actualNucLetter=rnaSequence[windowLast-1];

		// compute this position
		HmmType1::State state=hmm.GetFirstState();
		currTable[state]=0;

		while (1) {

			state++;
			if (state==hmm.GetLastState()) {
				break;
			}

			const float *table;
			float emitScore;
			if (hmm.IsEmittingState(state)) {
				if (actualNucLetter<MAXABET) {
					emitScore=hmm.GetSingletEmissionScore(state,actualNucLetter);
				}
				else {
					emitScore=0;
					assert(Alphabet_size==4);
					for (int concreteNuc = 0; concreteNuc <Alphabet_size; concreteNuc++) {

						if (Degenerate[actualNucLetter][concreteNuc]) {
							emitScore += hmm.GetSingletEmissionScore(state,concreteNuc) / (float)(DegenCount[actualNucLetter]);
						}
					}
				}
				table=prevTable;
			}
			else {
				emitScore=0;
				table=currTable;
			}

			float bestTransitionScore=(float)IMPOSSIBLE;
			for (int child=0; child<hmm.GetNumChildren(state); child++) {
				const float thisTransitionScore=table[hmm.GetNthChildState(state,child)] + hmm.GetNthChildTsc(state,child);
				bestTransitionScore=std::max(bestTransitionScore,thisTransitionScore);
			}
			currTable[state]=emitScore + bestTransitionScore;

			if (cykscanStats.collectHmmFullDynProgTable) {
				const float templatedHmmScore=cykscanStats.fullHmmDynProgTable[windowLast][state];
				if (fabs(templatedHmmScore - currTable[state])>1e-4) {
					assert(false);
					for (HmmType1::State i=hmm.GetFirstState(); i!=hmm.GetLastState(); i++) {
						printf("table[%d]=%g\n",i,table[i]);
					}
					throw SimpleStringException("Sanity check failed: HMM scores don't match.  I am insane.  windowLast=%d.  templatedHmmScore=%g, currTable[state]=%g, state=%d.  emitScore=%g",windowLast,templatedHmmScore,currTable[state],state,emitScore);
				}
			}
		}

		const float bestScoreAtPos=currTable[hmm.GetActualLastState()];
		if (bestScoreAtPos >= minLodScoreForHit) {
			//printf("bestScoreAtPos >= minLodScoreForHit\n");
			std::pair<int,int> thisWindow;
			thisWindow.first=std::max(startPosToScan,windowLast-windowLen);
			thisWindow.second=windowLast;

			AddWindowToList(hmmHitList,thisWindow);
		}

		std::swap(currTable,prevTable);
	}
}

void ScanHmm_HmmType1Float_NonTemplated (HitList& hmmHitList,const HitList& inputHitList,const HmmType1& hmm,float minLodScoreForHit,CykscanStats& cykscanStats,const char *rnaSequence, int windowLen)
{
	float *currTable=new float[hmm.GetNumStates()];
	float *prevTable=new float[hmm.GetNumStates()];

	hmmHitList.clear();

#ifndef CM2HMM_ONLY
	if (cmpWithPureHmm) {
		cykscanStats.collectHmmFullDynProgTable=true;
		HitList dummyHitList;
		ScanHmm_HmmType1Float(dummyHitList,inputHitList,hmm,minLodScoreForHit,cykscanStats,rnaSequence,windowLen);
	}
#endif

	HitList::const_iterator inputHitListIterator;
	for (inputHitListIterator=inputHitList.begin(); inputHitListIterator!=inputHitList.end(); inputHitListIterator++) {

		int startPosToScan=inputHitListIterator->first;
		int endPosToScan=inputHitListIterator->second;

		ScanHmm_HmmType1Float_NonTemplated_Window(hmmHitList,startPosToScan,endPosToScan,hmm,minLodScoreForHit,cykscanStats,rnaSequence,windowLen,prevTable,currTable);
	}


	delete [] currTable;
	delete [] prevTable;
}


///////////////////////
// HmmType1_OldSchool

HmmType1_OldSchool::HmmType1_OldSchool (void)
{
}
HmmType1_OldSchool::~HmmType1_OldSchool ()
{
	for (int nuc=0; nuc<MAXABET; nuc++) {
		delete [] esc[nuc];
	}
	delete [] esc;
	delete [] stateInfo;
}
void HmmType1_OldSchool::Init (const HmmType1& hmm)
{
	numStates=hmm.GetNumStates();
	stateInfo=new StateInfo [numStates];
	esc=new float *[MAXABET];
	for (int nuc=0; nuc<MAXABET; nuc++) {
		esc[nuc]=new float[numStates];
	}

	for (State state=0; state<numStates; state++) {
		if (hmm.IsEmittingState(state)) {
			for (int nuc=0; nuc<MAXABET; nuc++) {
				esc[nuc][state]=hmm.GetSingletEmissionScore(state,nuc);
			}
		}
		stateInfo[state].isEmitting=hmm.IsEmittingState(state);
		stateInfo[state].numChildren=hmm.GetNumChildren(state);
		if (stateInfo[state].numChildren>0) {
			stateInfo[state].firstChild=hmm.GetNthChildState(state,0);
		}
		for (int child=0; child<stateInfo[state].numChildren; child++) {
			stateInfo[state].tsc[child]=hmm.GetNthChildTsc(state,child);
		}
	}
}




void ScanHmm_HmmType1Float_NonTemplated_Window(HitList& hmmHitList,int startPosToScan,int endPosToScan,const HmmType1_OldSchool& hmm,float minLodScoreForHit,CykscanStats& cykscanStats,const char *rnaSequence, int windowLen,float *prevTable,float *currTable)
{
	// init table
	HmmType1_OldSchool::State state=0;
	prevTable[state]=0;
	while (1) {

		state++;
		if (state==hmm.numStates) {
			break;
		}

		if (hmm.stateInfo[state].isEmitting) {
			// during initialization, nothing can emit
			prevTable[state]=(float)IMPOSSIBLE;
		}
		else {
			float bestTransitionScore=(float)IMPOSSIBLE;
			for (int child=0; child<hmm.stateInfo[state].numChildren; child++) {
				const float thisTransitionScore=prevTable[hmm.stateInfo[state].firstChild + child] + hmm.stateInfo[state].tsc[child];
				bestTransitionScore=std::max(bestTransitionScore,thisTransitionScore);
			}
			prevTable[state]=bestTransitionScore;
		}

		if (cykscanStats.collectHmmFullDynProgTable) {
			const float templatedHmmScore=cykscanStats.fullHmmDynProgTable[0][state];
			if (fabs(templatedHmmScore - prevTable[state])>2e-6) {
				assert(false);
				throw SimpleStringException("Sanity check failed: HMM scores don't match.  I am insane.");
			}
		}
	}

	// scan
	for (int windowLast=startPosToScan+1; windowLast<=endPosToScan; windowLast++) {
		const char actualNucLetter=rnaSequence[windowLast-1];

		// compute this position
		HmmType1::State state=0;
		currTable[state]=0;

		while (1) {

			state++;
			if (state==hmm.numStates) {
				break;
			}

			const float *table;
			float emitScore;
			if (hmm.stateInfo[state].isEmitting) {
				if (actualNucLetter<MAXABET) {
					emitScore=hmm.esc[actualNucLetter][state];
				}
				else {
					emitScore=0;
					assert(Alphabet_size==4);
					for (int concreteNuc = 0; concreteNuc <Alphabet_size; concreteNuc++) {

						if (Degenerate[actualNucLetter][concreteNuc]) {
							emitScore += hmm.esc[concreteNuc][state] / (float)(DegenCount[actualNucLetter]);
						}
					}
				}
				table=prevTable;
			}
			else {
				emitScore=0;
				table=currTable;
			}

			float bestTransitionScore=(float)IMPOSSIBLE;
			for (int child=0; child<hmm.stateInfo[state].numChildren; child++) {
				const float thisTransitionScore=table[hmm.stateInfo[state].firstChild + child] + hmm.stateInfo[state].tsc[child];
				bestTransitionScore=std::max(bestTransitionScore,thisTransitionScore);
			}
			currTable[state]=emitScore + bestTransitionScore;

			if (cykscanStats.collectHmmFullDynProgTable) {
				const float templatedHmmScore=cykscanStats.fullHmmDynProgTable[windowLast][state];
				if (fabs(templatedHmmScore - currTable[state])>2e-6) {
					assert(false);
					throw SimpleStringException("Sanity check failed: HMM scores don't match.  I am insane.");
				}
			}
		}

		const float bestScoreAtPos=currTable[state];
		if (bestScoreAtPos >= minLodScoreForHit) {
			std::pair<int,int> thisWindow;
			thisWindow.first=std::max(startPosToScan,windowLast-windowLen);
			thisWindow.second=windowLast;

			AddWindowToList(hmmHitList,thisWindow);
		}

		std::swap(currTable,prevTable);
	}
}

void ScanHmm_HmmType1Float_NonTemplated (HitList& hmmHitList,const HitList& inputHitList,const HmmType1& hmmType1,const HmmType1_OldSchool& hmm,float minLodScoreForHit,CykscanStats& cykscanStats,const char *rnaSequence, int windowLen)
{
	float *currTable=new float[hmm.numStates];
	float *prevTable=new float[hmm.numStates];

	hmmHitList.clear();

#ifndef CM2HMM_ONLY
	if (cmpWithPureHmm) {
		cykscanStats.collectHmmFullDynProgTable=true;
		HitList dummyHitList;
		ScanHmm_HmmType1Float(dummyHitList,inputHitList,hmmType1,minLodScoreForHit,cykscanStats,rnaSequence,windowLen);
	}
#endif

	HitList::const_iterator inputHitListIterator;
	for (inputHitListIterator=inputHitList.begin(); inputHitListIterator!=inputHitList.end(); inputHitListIterator++) {

		int startPosToScan=inputHitListIterator->first;
		int endPosToScan=inputHitListIterator->second;

		ScanHmm_HmmType1Float_NonTemplated_Window(hmmHitList,startPosToScan,endPosToScan,hmm,minLodScoreForHit,cykscanStats,rnaSequence,windowLen,prevTable,currTable);
	}


	delete [] currTable;
	delete [] prevTable;
}

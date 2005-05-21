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

/*
namespace lpsolve {
#include <lpkit.h>
}
*/

// note: what I call "normal states" and "(post) insert states", Sean Eddy calls "split set" and "insert set" in Eddy, S.R. (2002), BMC Bioinformatics, 3:18 (see Table 3 and text on that page).  I guess I should ideally have read that paper before doing this...

// my LP problems were unbounded, so just in case the problem recurs (this was from my very early work, but the lower bound seems pretty conservative)
#define HACK_LOWERBOUND -8000.0
#define PARAMETER_GUESS 0.0

#define ENABLE_CACHING // building the HMM can take kind of a while.  Caching is somewhat deprecated (it seems better to deal with loading/saving HMMs in specific files), and I don't expect it to make it into a production version of this code.

#ifdef CM2HMM_ONLY
// this caching thing doesn't make sense for distributed code -- it's just a complication to explain
#undef ENABLE_CACHING
#endif

//#define DEBUG_DUMP // enable dumping of HMM build work, into hmmdump.txt
FILE *dumpFile=NULL;

// THE FOLLOWING EQUIVOCATION IS NO LONGER NECESSARY - THIS VARIABLE WILL ALWAYS BE SET TO 'TRUE'
// I'm not sure if this is necessary/good.  currently I think it is
// if true, then every CM node has at least 1 left & at least 1 right node in the HMM
// e.g. a MATL node will have a dummy PASSTHRU node on the right
// if false, MATL nodes wouldn't have anything on the right
// advantage of true: it's well-defined how to move thru the HMM based on the CM, whereas
// without the passthru nodes, it's weird.  Also, without the passthru node, we get more info
// in the HMM, but this isn't reflected in the CM, so it doesn't help with an upper bound
// of the CM (e.g. a CM has ROOT,MATL,MATL,MATL,MATL,MATR.  if everyCmNodeHasLeftAndRightHmmNodes==false,
// we can hook up the ROOT directly to MATR & make better decisions on what state to enter in the MATR,
// but the CM doesn't do this anyway).
#define everyCmNodeHasLeftAndRightHmmNodes true

bool IsFirstHmmBuilt(HmmCommittee::CommitteeBuildInfo *committeeBuildInfo)
{
	if (committeeBuildInfo==NULL) {
		return true;
	}
	else {
		return committeeBuildInfo->currCommitteeMemberNum==0;
	}
}

// the HMM consists of Left-emitting stuff from the CM, and Right-emitting stuff from the CM
enum HMM_Direction {
	HMM_Left,HMM_Right
};

struct CmNodeTypeData {
	// total # of states in the Covariance Model
	int numCmStates;
	// all nodes start with something normal, like a S_st, MP, ML, MR, D, that has to do with their function.  Then they can have an IL and/or IR state.  the normal states of a node always connect to the insert states of that node and the normal states of the next node.  The insert states of a node connect only to the normal states of the next node.
	int numNormalStates,numInsertPostStates;
	// state types, for verifying that they're the way I expect
	const int *stateTypes;
	
	// starts/ends segment
	bool isStartingNode,isEndingNode;

	// stuff about HMM encoding
	// # of HMM states to the left, # to the right
	int numHmmLeftStates,numHmmRightStates;
	// how many of those states are IL/IR
	int numHmmLeftInsertStates,numHmmRightInsertStates; // these are either 0 or 1
};
static int BEGR_stateTypes[]={S_st,IL_st,INVALID_st},
	BEGL_stateTypes[]={S_st,INVALID_st},
	ROOT_stateTypes[]={S_st,IL_st,IR_st,INVALID_st},
	BIF_stateTypes[]={B_st,INVALID_st},
	END_stateTypes[]={E_st,INVALID_st},
	MATP_stateTypes[]={MP_st,ML_st,MR_st,D_st,IL_st,IR_st,INVALID_st},
	MATL_stateTypes[]={ML_st,D_st,IL_st,INVALID_st},
	MATR_stateTypes[]={MR_st,D_st,IR_st,INVALID_st};
static CmNodeTypeData BEGR_data,BEGL_data,ROOT_data,BIF_data,END_data,MATP_data,MATL_data,MATR_data;
static bool isInit_CmNodeTypeData=false;
void InferStuffFromStateTypes(CmNodeTypeData& data,const int *stateTypes)
{
	data.stateTypes=stateTypes; // that was easy :-)

	data.isStartingNode=stateTypes[0]==S_st;
	data.isEndingNode=stateTypes[0]==E_st || stateTypes[0]==EL_st || stateTypes[0]==B_st;

	data.numCmStates=0;
	data.numNormalStates=0;
	data.numInsertPostStates=0;
	data.numHmmLeftStates=0;
	data.numHmmRightStates=0;
	data.numHmmLeftInsertStates=0;
	data.numHmmRightInsertStates=0;
	for (int i=0; stateTypes[i]!=INVALID_st; i++) {
		bool isInsertState=stateTypes[i]==IL_st || stateTypes[i]==IR_st;
		bool isStartingOrEndingState=stateTypes[i]==S_st || stateTypes[i]==E_st || stateTypes[i]==EL_st || stateTypes[i]==B_st;
		bool isLeft=stateTypes[i]==IL_st || stateTypes[i]==ML_st || stateTypes[i]==MP_st || stateTypes[i]==D_st || isStartingOrEndingState;
		bool isRight=stateTypes[i]==IR_st || stateTypes[i]==MR_st || stateTypes[i]==MP_st || stateTypes[i]==D_st || isStartingOrEndingState;

		data.numCmStates++;
		if (isInsertState) {
			data.numInsertPostStates++;
		}
		else {
			data.numNormalStates++;
		}
		if (stateTypes[i]==MP_st) {
			// don't put anything extra for HMM for MP_st - ML_st and MR_st should put it in
		}
		else {
			if (isLeft) {
				data.numHmmLeftStates++;
				if (isInsertState) {
					data.numHmmLeftInsertStates++;
				}
			}
			if (isRight) {
				data.numHmmRightStates++;
				if (isInsertState) {
					data.numHmmRightInsertStates++;
				}
			}
		}
	}

	if (stateTypes[0]==ML_st) {
		data.numHmmRightStates--; // generic alg overcounts because it thinks the D_st is on both sides
	}
	if (stateTypes[0]==MR_st) {
		data.numHmmLeftStates--; // generic alg overcounts because it thinks the D_st is on both sides
	}

	if (everyCmNodeHasLeftAndRightHmmNodes) {
		if (data.numHmmLeftStates==0) {
			data.numHmmLeftStates++;
		}
		if (data.numHmmRightStates==0) {
			data.numHmmRightStates++;
		}
	}
}
void InitCmNodeTypeData(void)
{
	InferStuffFromStateTypes(BEGR_data,BEGR_stateTypes);
	InferStuffFromStateTypes(BEGL_data,BEGL_stateTypes);
	InferStuffFromStateTypes(ROOT_data,ROOT_stateTypes);
	InferStuffFromStateTypes(BIF_data,BIF_stateTypes);
	InferStuffFromStateTypes(END_data,END_stateTypes);
	InferStuffFromStateTypes(MATP_data,MATP_stateTypes);
	InferStuffFromStateTypes(MATL_data,MATL_stateTypes);
	InferStuffFromStateTypes(MATR_data,MATR_stateTypes);

	isInit_CmNodeTypeData=true;
}
CmNodeTypeData GetCmNodeTypeData (int nodeType)
{
	if (!isInit_CmNodeTypeData) {
		InitCmNodeTypeData();
	}

	switch (nodeType) {
case BIF_nd:
	return BIF_data;
case END_nd:
	return END_data;
case MATP_nd:
	return MATP_data;
case MATL_nd:
	return MATL_data;
case MATR_nd:
	return MATR_data;
case BEGL_nd:
	return BEGL_data;
case BEGR_nd:
	return BEGR_data;
case ROOT_nd:
	return ROOT_data;
case DUMMY_nd:
	assert(false); // I thought DUMMY nodes get removed once the CM is made
	break;
default:
	assert(false);
	Die("unknown node type when converting CM to HMM\n");
	break;
	}

	// keep compiler happy
	CmNodeTypeData data;
	data.isStartingNode=false;
	return data;
}

// make sure we got everything about CovarianceModels right
// don't use asserts -- I want to run this at run time
bool ValidateCmNodeTypeData (const CovarianceModel& cm)
{
	CovarianceModel::Node node;
	for (node=cm.GetFirstNode(); node!=cm.GetLastNode(); node++) {
		int thisNodeType=cm.GetNodeType(node);
		CmNodeTypeData data=GetCmNodeTypeData(thisNodeType);
		CovarianceModel::State state;
		int numCmStates=0;
		int numInsertPostStates=0;
		for (state=cm.GetFirstStateOfNode(node); state!=cm.GetLastStateOfNode(node); state++) {
			if (cm.GetStateType(state)==IL_st || cm.GetStateType(state)==IR_st) {
				numInsertPostStates++;
			}

			CovarianceModel::Node nextNode(node);
			nextNode++;

			// the curr node is not an ending node, then the next node should be physically next
			if (nextNode!=cm.GetLastNode() && !cm.IsBifuricationNode(node) && !cm.IsEndNode(node)) {
				CovarianceModel::State state=cm.GetFirstStateOfNode(node); // I don't trust GetLastStateOfNode, since I might change the alg
				while (cm.GetNode(state)==node) {
					state++;
				}
				// now that we're past 'node', we should be at 'nextNode', since nextNode should be physically next, as well as numerically
				if (cm.GetNode(state)!=nextNode) {
					assert(false);
					return false;
				}
			}

			if (nextNode!=cm.GetLastNode() && !cm.IsBifuricationNode(node)) { // last node doesn't link to anything & I'll ignore bifurication nodes
				CmNodeTypeData nextNodeData=GetCmNodeTypeData(cm.GetNodeType(nextNode));
				CovarianceModel::State expectedFirstChild;
				int numTransitions=cm.GetNumChildren(state);
				int expectedNumTransitions;

				if (cm.GetNodeType(node)==END_nd) {
					// special case for end nodes
					expectedNumTransitions=0;
				}
				else {
					if (numCmStates<data.numNormalStates) {
						// should link to IL/IR states + all normal states in the next node
						expectedNumTransitions=data.numInsertPostStates + nextNodeData.numNormalStates;
						expectedFirstChild=cm.GetFirstStateOfNode(node).PlusInt(data.numNormalStates);
					}
					else {
						// isa IL/IR state.  should link to itself and any remaining IL/IR states & all normal states im next node
						expectedNumTransitions=(data.numCmStates-numCmStates) + nextNodeData.numNormalStates;
						expectedFirstChild=state; // I should be my first child if I'm an IL/IR
					}
				}

				// verify we have the exact child states we expect
				if (numTransitions!=expectedNumTransitions) {
					assert(false);
					return false;
				}
				for (int childNum=0; childNum<numTransitions; childNum++) {
					if (cm.GetNthChildState(state,childNum)!=expectedFirstChild.PlusInt(childNum)) {
						assert(false);
						return false;
					}
				}
			}

			numCmStates++;
		}
		if (cm.GetNodeType(node)!=BIF_nd) { // ignore bifurications
			if (data.numCmStates!=numCmStates || data.numInsertPostStates!=numInsertPostStates || data.numNormalStates+data.numInsertPostStates!=data.numCmStates) {
				assert(false);
				return false;
			}
		}
	}
	return true;
}

int GetNumExtraMATPStatesVersusOriginalHmmBuildType (Cm2Hmm_HmmBuildType hmmType)
{
	switch (hmmType) {
		case HmmBuildType_Original:
			return 0;
		case HmmBuildType_separateMPandMLMR:
			return 1;
		case HmmBuildType_separateMPMLMRD:
			return 2;
	}
	assert(false);
	throw SimpleStringException("internal error %s:%d",__FILE__,__LINE__);
}

int GetNumNormalHmmStatesInNextLeftNode(const CovarianceModel& cm,CovarianceModel::Node currNode,Cm2Hmm_HmmBuildType hmmBuildType)
{
	// check that it's not an ending node, in which case there is no next node
	{
		int currNodeType=cm.GetNodeType(currNode);
		CmNodeTypeData data=GetCmNodeTypeData(currNodeType);
		if (data.isEndingNode) {
			return 0;
		}
	}

	// else, see what's next
	currNode++;
	while (1) {
		int currNodeType=cm.GetNodeType(currNode);
		CmNodeTypeData data=GetCmNodeTypeData(currNodeType);
		if (data.isEndingNode) {
			return 1;
		}

		if (data.numHmmLeftStates>0) {
			int num=data.numHmmLeftStates - data.numHmmLeftInsertStates;
			if (currNodeType==MATP_nd) {
				num += GetNumExtraMATPStatesVersusOriginalHmmBuildType(hmmBuildType);
			}
			return num;
		}

		currNode++;
	}
}

CovarianceModel::Node FindPreviousRightNode_NotAtStart(const CovarianceModel& cm,CovarianceModel::Node currNode)
{
	// check that it's not a starting node, in which case there is no previous node, and this violates our precondition
	assert(!GetCmNodeTypeData(cm.GetNodeType(currNode)).isStartingNode);

	// else see what's next
	currNode--;
	while (1) {
		int currNodeType=cm.GetNodeType(currNode);
		CmNodeTypeData data=GetCmNodeTypeData(currNodeType);
		if (data.isStartingNode) {
			return currNode;
		}

		if (data.numHmmLeftStates>0) {
			return currNode;
		}

		currNode++;
	}
}

// find the first (lowest-numbered) state in the given CM node that has a state on the Right side of the HMM
CovarianceModel::State GetFirstHmmRightState(const HmmAndBuildInfo& newHMM,const CovarianceModel& cm,CovarianceModel::Node prevRightNodeInCM)
{
	CovarianceModel::State firstCmState=cm.GetFirstStateOfNode(prevRightNodeInCM);
	CovarianceModel::State lastCmState=cm.GetLastStateOfNode(prevRightNodeInCM);

	CovarianceModel::State firstState=CovarianceModel::GetInvalidState();
	CovarianceModel::State cmState;
	for (cmState=firstCmState; cmState!=lastCmState; cmState++) {
		CovarianceModel::State hmmRightState=newHMM.cm2HmmState[cmState].hmmRightState;
		if (hmmRightState!=CovarianceModel::GetInvalidState()) {
			if (firstState==CovarianceModel::GetInvalidState()) {
				firstState=hmmRightState;
			}
			else {
				if (CovarianceModel::StateToInt(hmmRightState) < CovarianceModel::StateToInt(firstState)) {
					firstState=hmmRightState;
				}
			}
		}
	}

	assert(firstState!=CovarianceModel::GetInvalidState()); // should have got something
	return firstState;
}

// does all the work to set up the children of an HMM state, based on a CM source state
// the tricky thing about this function is that the children are relative to a start-to-end traversal of the HMM.  Thus, on the HMM that corresponds to Right-side emissions in the CM, we have to set up the children backwards, since CM nodes appear in the reverse order on the right-emitting side.
void SetupChildren (HmmAndBuildInfo& newHMM,const CovarianceModel& cm,CovarianceModel::State cmState,CovarianceModel::Node cmNode,HMM_Direction hmmDirection,CovarianceModel::Node cmFirstNode,CovarianceModel::Node cmLastNode,Cm2Hmm_HmmBuildType hmmBuildType)
{
	assert(hmmDirection==HMM_Left || hmmDirection==HMM_Right);

	CovarianceModel::State firstCmStateOfNode=cm.GetFirstStateOfNode(cmNode);

	CovarianceModel::State hmmState;
	if (hmmDirection==HMM_Left) {
		hmmState=newHMM.cm2HmmState[cmState].hmmLeftState;
	}
	else {
		hmmState=newHMM.cm2HmmState[cmState].hmmRightState;
	}

	int thisNodeType=cm.GetNodeType(cmNode);
	CmNodeTypeData thisNodeData=GetCmNodeTypeData(thisNodeType);

	assert(CovarianceModel::IsStateInRange(cmState,firstCmStateOfNode,firstCmStateOfNode.PlusInt(thisNodeData.numCmStates)));

	int thisNode_separate_extraStates=0;
	if (thisNodeType==MATP_nd) {
		thisNode_separate_extraStates=GetNumExtraMATPStatesVersusOriginalHmmBuildType(hmmBuildType);
	}

	if (hmmDirection==HMM_Left) {

		int numNormalStatesInNextLeftNode=GetNumNormalHmmStatesInNextLeftNode(cm,cmNode,hmmBuildType);
		if (CovarianceModel::IsStateInRange(cmState,firstCmStateOfNode,firstCmStateOfNode.PlusInt(thisNodeData.numNormalStates))) {
			// normal state on left side is connected to insert states in this node on the left side & normal states in next node in CM

			newHMM.hmm.SetFirstChild(hmmState,newHMM.cm2HmmState[firstCmStateOfNode].hmmLeftState.PlusInt(thisNodeData.numHmmLeftStates+thisNode_separate_extraStates-thisNodeData.numHmmLeftInsertStates)); // first insert state in this node
			newHMM.hmm.SetNumChildren(hmmState,thisNodeData.numHmmLeftInsertStates + numNormalStatesInNextLeftNode);
		}
		else {
			assert(cm.GetStateType(cmState)!=MP_st); // this shouldn't happen for MP state, since separateMPandMLMR logic is not applied here

			// insert state on left side is connected to itself and normal states in the next node in the CM
			newHMM.hmm.SetFirstChild(hmmState,newHMM.cm2HmmState[cmState].hmmLeftState);
			assert(thisNodeData.numHmmLeftInsertStates==1); // I'm an insert state, and I can only have <=1 IL_st, and <=1 IR_st
			newHMM.hmm.SetNumChildren(hmmState,1 + numNormalStatesInNextLeftNode);
		}
	}
	else {
		assert(hmmDirection==HMM_Right);
		if (CovarianceModel::IsStateInRange(cmState,firstCmStateOfNode,firstCmStateOfNode.PlusInt(thisNodeData.numNormalStates))) {
			// normal state on right side is connected to all states (insert and normal) of the _previous_ CM node that emits on the right
			if (thisNodeData.isStartingNode || cmNode==cmFirstNode) { // if it's an actual start node, or if it's acting like a start node (presumably because it's in DoLocal mode)
				// start node on right is connected to nothing (it's really an end node);
				newHMM.hmm.SetNoChildren(hmmState);
			}
			else {
				CovarianceModel::Node prevRightNode=FindPreviousRightNode_NotAtStart(cm,cmNode);
				//CovarianceModel::State firstStateOfPrevRightNode=cm.GetFirstStateOfNode(prevRightNode);
				int prevRightNodeType=cm.GetNodeType(prevRightNode);
				CmNodeTypeData prevRightNodeData=GetCmNodeTypeData(prevRightNodeType);
				newHMM.hmm.SetFirstChild(hmmState,GetFirstHmmRightState(newHMM,cm,prevRightNode));
				int extraStateForPrevNode=0;
				if (prevRightNodeType==MATP_nd) {
					extraStateForPrevNode=GetNumExtraMATPStatesVersusOriginalHmmBuildType(hmmBuildType);
				}
				newHMM.hmm.SetNumChildren(hmmState,prevRightNodeData.numHmmRightStates + extraStateForPrevNode);
			}
		}
		else {
			// insert state on right side is connected to itself & all normal states of the _current_ CM node that emit on the right
			newHMM.hmm.SetFirstChild(hmmState,GetFirstHmmRightState(newHMM,cm,cmNode));
			newHMM.hmm.SetNumChildren(hmmState,thisNodeData.numHmmRightStates+thisNode_separate_extraStates);
		}
	}
}

void ReverseMapCm2HmmState(HmmAndBuildInfo& newHMM,const CovarianceModel& cm)
{
	// set up reverse mapping 'hmm2CmStateVector'
	// WARNING: might be ambiguous since some HMM states are used by multiple CM states, but for now it doesn't matter.  So, I'm using a list of states
	newHMM.hmm2CmStateVector.clear(); // forget anything that's here already
	newHMM.hmm2CmStateVector.resize(newHMM.hmm.GetNumStates());
	CovarianceModel::State cmState;
	for (cmState=cm.GetFirstState(); cmState!=cm.GetLastState(); cmState++) {

		CovarianceModel::State hmmLeftState=newHMM.cm2HmmState[cmState].hmmLeftState;
		if (hmmLeftState!=CovarianceModel::GetInvalidState()) {
			newHMM.hmm2CmStateVector[hmmLeftState].push_back(cmState);
		}

		CovarianceModel::State hmmRightState=newHMM.cm2HmmState[cmState].hmmRightState;
		if (hmmRightState!=CovarianceModel::GetInvalidState()) {
			newHMM.hmm2CmStateVector[hmmRightState].push_back(cmState);
		}
	}
}

/*
Converts the block of the CM in the closed interval [firstNode,lastNode]
to an HMM structure (i.e. not filling in the transition/emission weights)

firstNode must be a ROOT or BEGIN node
lastNode must be END, local END or BIFURICATION NODE
there cannot be any bifurication node in the range firstNode..(lastNode-1)
*/
void Cm2Hmm_Structurally_Block (HmmAndBuildInfo& newHMM,const CovarianceModel& sourceCM,CovarianceModel::Node firstNode,CovarianceModel::Node lastNode,ExtraCm2HmmInfo *extraCm2HmmInfo,Cm2Hmm_HmmBuildType hmmBuildType)
{
	assert(sourceCM.IsRootOrBeginNode(firstNode) || sourceCM.DoLocal()); // firstNode must be BEGL, BEGR or ROOT, unless we're doing local, in which case we skip the root
	assert(sourceCM.IsEndNode(lastNode) || sourceCM.IsBifuricationNode(lastNode)); // last state must be END or BIF

	CovarianceModel::Node lastNodePlus1(lastNode);
	lastNodePlus1++;

	// other sanity checks
	CovarianceModel::Node cmNode;
	for (cmNode=firstNode; cmNode!=lastNode; cmNode++) {
		assert(cmNode==lastNode || !sourceCM.IsBifuricationNode(cmNode)); // can't have bifurication node, except at end
		assert(!(sourceCM.GetNodeType(cmNode)==DUMMY_nd)); // I thought those were removed by the time the CM was built
	}

	// work out how many states we'll need in the HMM
	int numHmmStates=0;
	for (CovarianceModel::Node node=firstNode; node!=lastNodePlus1; node++) {

		const CmNodeTypeData nodeTypeData=GetCmNodeTypeData(sourceCM.GetNodeType(node));
		numHmmStates += nodeTypeData.numHmmLeftStates + nodeTypeData.numHmmRightStates;

		if (sourceCM.GetNodeType(node)==MATP_nd) {
			numHmmStates += 2*GetNumExtraMATPStatesVersusOriginalHmmBuildType(hmmBuildType);
		}
	}
	numHmmStates--; // the end states are over counted since they're both left & right

	newHMM.hmm.Init(numHmmStates);
	Cm2HmmState nullCm2HmmState;
	nullCm2HmmState.hmmLeftState=CovarianceModel::GetInvalidState();
	nullCm2HmmState.hmmRightState=CovarianceModel::GetInvalidState();
	newHMM.cm2HmmState.assign(sourceCM.GetNumStates(),nullCm2HmmState);
	CovarianceModel::State firstHmmState=newHMM.hmm.GetFirstState(),lastHmmState=newHMM.hmm.GetLastState();

	// first set up the new HMM states & which CM states map to them (this info will be used in setting up the children of those states)
	for (cmNode=firstNode; cmNode!=lastNodePlus1; cmNode++) {

		assert(firstHmmState!=lastHmmState); // else we miscalculated somewhere

		int thisNodeType=sourceCM.GetNodeType(cmNode);
		CmNodeTypeData thisNodeData=GetCmNodeTypeData(thisNodeType);
		CovarianceModel::State firstCmState=sourceCM.GetFirstStateOfNode(cmNode);
		CovarianceModel::State lastCmState=sourceCM.GetLastStateOfNode(cmNode);

		if (thisNodeData.isEndingNode) {
			lastHmmState--;  assert(firstHmmState==lastHmmState); // end in the middle
			newHMM.hmm.SetStateType(firstHmmState,PASSTHRU_st);
			newHMM.cm2HmmState[firstCmState].hmmLeftState=firstHmmState;
			newHMM.cm2HmmState[firstCmState].hmmRightState=firstHmmState;
			assert(thisNodeData.numCmStates==1); // should have checked this already, but whatever
			newHMM.leftToRightPassthruState=firstHmmState;
			firstHmmState++;
		}
		else {
			// do normal states of node
			switch (sourceCM.GetNodeType(cmNode)) {
				case ROOT_nd:
				case BEGL_nd:
				case BEGR_nd:
					// left side start
					newHMM.hmm.SetStateType(firstHmmState,S_st);
					newHMM.cm2HmmState[firstCmState].hmmLeftState=firstHmmState;
					firstHmmState++;
					// right side start
					lastHmmState--;
					newHMM.hmm.SetStateType(lastHmmState,E_st);
					newHMM.cm2HmmState[firstCmState].hmmRightState=lastHmmState;
					break;
				case MATL_nd:
					// ML_st
					newHMM.hmm.SetStateType(firstHmmState,ML_st);
					newHMM.cm2HmmState[firstCmState].hmmLeftState=firstHmmState;
					newHMM.cm2HmmState[firstCmState].hmmRightState=CovarianceModel::GetInvalidState();
					firstHmmState++;
					// D_st
					newHMM.hmm.SetStateType(firstHmmState,D_st);
					newHMM.cm2HmmState[firstCmState.PlusInt(1)].hmmLeftState=firstHmmState;
					newHMM.cm2HmmState[firstCmState.PlusInt(1)].hmmRightState=CovarianceModel::GetInvalidState();
					firstHmmState++;
					if (everyCmNodeHasLeftAndRightHmmNodes) {
						lastHmmState--;
						newHMM.hmm.SetStateType(lastHmmState,PASSTHRU_st);
						newHMM.cm2HmmState[firstCmState].hmmRightState=lastHmmState;
						newHMM.cm2HmmState[firstCmState.PlusInt(1)].hmmRightState=lastHmmState;
					}
					break;
				case MATR_nd:
					// MR_st
					lastHmmState--;
					newHMM.hmm.SetStateType(lastHmmState,ML_st);
					newHMM.cm2HmmState[firstCmState].hmmLeftState=CovarianceModel::GetInvalidState();
					newHMM.cm2HmmState[firstCmState].hmmRightState=lastHmmState;
					// D_st
					lastHmmState--;
					newHMM.hmm.SetStateType(lastHmmState,D_st);
					newHMM.cm2HmmState[firstCmState.PlusInt(1)].hmmLeftState=CovarianceModel::GetInvalidState();
					newHMM.cm2HmmState[firstCmState.PlusInt(1)].hmmRightState=lastHmmState;
					if (everyCmNodeHasLeftAndRightHmmNodes) {
						newHMM.hmm.SetStateType(firstHmmState,PASSTHRU_st);
						newHMM.cm2HmmState[firstCmState].hmmLeftState=firstHmmState;
						newHMM.cm2HmmState[firstCmState.PlusInt(1)].hmmLeftState=firstHmmState;
						firstHmmState++;
					}
					break;
				case MATP_nd:
					switch (hmmBuildType) {
						case HmmBuildType_Original: // also called "Compacted" or "type 0"
							// normal case

							// left match
							newHMM.hmm.SetStateType(firstHmmState,ML_st);
							// left delete
							newHMM.hmm.SetStateType(firstHmmState.PlusInt(1),D_st);
							// right stuff
							lastHmmState--;
							lastHmmState--;
							// right match
							newHMM.hmm.SetStateType(lastHmmState,ML_st);
							// right delete
							newHMM.hmm.SetStateType(lastHmmState.PlusInt(1),D_st);

							// hook up CM to equivalent HMM states
							// MP_st
							newHMM.cm2HmmState[firstCmState.PlusInt(0)].hmmLeftState=firstHmmState;
							newHMM.cm2HmmState[firstCmState.PlusInt(0)].hmmRightState=lastHmmState;
							// ML_st
							newHMM.cm2HmmState[firstCmState.PlusInt(1)].hmmLeftState=firstHmmState;
							newHMM.cm2HmmState[firstCmState.PlusInt(1)].hmmRightState=lastHmmState.PlusInt(1);
							// MR_st
							newHMM.cm2HmmState[firstCmState.PlusInt(2)].hmmLeftState=firstHmmState.PlusInt(1);
							newHMM.cm2HmmState[firstCmState.PlusInt(2)].hmmRightState=lastHmmState;
							// D_st
							newHMM.cm2HmmState[firstCmState.PlusInt(3)].hmmLeftState=firstHmmState.PlusInt(1);
							newHMM.cm2HmmState[firstCmState.PlusInt(3)].hmmRightState=lastHmmState.PlusInt(1);

							firstHmmState++;
							firstHmmState++;
							break;
						case HmmBuildType_separateMPandMLMR: // also called "expanded" or "type 1"

							// with extra states

							// left match MP
							newHMM.hmm.SetStateType(firstHmmState,ML_st);
							// left match ML
							newHMM.hmm.SetStateType(firstHmmState.PlusInt(1),ML_st);
							// left delete
							newHMM.hmm.SetStateType(firstHmmState.PlusInt(2),D_st);
							// right stuff
							lastHmmState--;
							lastHmmState--;
							lastHmmState--;
							// right match MP
							newHMM.hmm.SetStateType(lastHmmState,ML_st);
							// right match MR
							newHMM.hmm.SetStateType(lastHmmState.PlusInt(1),ML_st);
							// right delete
							newHMM.hmm.SetStateType(lastHmmState.PlusInt(2),D_st);

							// hook up CM to equivalent HMM states
							// MP_st
							newHMM.cm2HmmState[firstCmState.PlusInt(0)].hmmLeftState=firstHmmState;
							newHMM.cm2HmmState[firstCmState.PlusInt(0)].hmmRightState=lastHmmState;
							// ML_st
							newHMM.cm2HmmState[firstCmState.PlusInt(1)].hmmLeftState=firstHmmState.PlusInt(1);
							newHMM.cm2HmmState[firstCmState.PlusInt(1)].hmmRightState=lastHmmState.PlusInt(2);
							// MR_st
							newHMM.cm2HmmState[firstCmState.PlusInt(2)].hmmLeftState=firstHmmState.PlusInt(2);
							newHMM.cm2HmmState[firstCmState.PlusInt(2)].hmmRightState=lastHmmState.PlusInt(1);
							// D_st
							newHMM.cm2HmmState[firstCmState.PlusInt(3)].hmmLeftState=firstHmmState.PlusInt(2);
							newHMM.cm2HmmState[firstCmState.PlusInt(3)].hmmRightState=lastHmmState.PlusInt(2);

							firstHmmState++;
							firstHmmState++;
							firstHmmState++;
							break;
						case HmmBuildType_separateMPMLMRD:

							// with more extra states

							// left match MP
							newHMM.hmm.SetStateType(firstHmmState,ML_st);
							// left match ML
							newHMM.hmm.SetStateType(firstHmmState.PlusInt(1),ML_st);
							// left delete MR
							newHMM.hmm.SetStateType(firstHmmState.PlusInt(2),D_st);
							// left delete D
							newHMM.hmm.SetStateType(firstHmmState.PlusInt(3),D_st);
							// right stuff
							lastHmmState--;
							lastHmmState--;
							lastHmmState--;
							lastHmmState--;
							// right match MP
							newHMM.hmm.SetStateType(lastHmmState,ML_st);
							// right match MR
							newHMM.hmm.SetStateType(lastHmmState.PlusInt(1),ML_st);
							// right delete ML
							newHMM.hmm.SetStateType(lastHmmState.PlusInt(2),D_st);
							// right delete D
							newHMM.hmm.SetStateType(lastHmmState.PlusInt(3),D_st);

							// hook up CM to equivalent HMM states
							// MP_st
							newHMM.cm2HmmState[firstCmState.PlusInt(0)].hmmLeftState=firstHmmState;
							newHMM.cm2HmmState[firstCmState.PlusInt(0)].hmmRightState=lastHmmState;
							// ML_st
							newHMM.cm2HmmState[firstCmState.PlusInt(1)].hmmLeftState=firstHmmState.PlusInt(1);
							newHMM.cm2HmmState[firstCmState.PlusInt(1)].hmmRightState=lastHmmState.PlusInt(2);
							// MR_st
							newHMM.cm2HmmState[firstCmState.PlusInt(2)].hmmLeftState=firstHmmState.PlusInt(2);
							newHMM.cm2HmmState[firstCmState.PlusInt(2)].hmmRightState=lastHmmState.PlusInt(1);
							// D_st
							newHMM.cm2HmmState[firstCmState.PlusInt(3)].hmmLeftState=firstHmmState.PlusInt(3);
							newHMM.cm2HmmState[firstCmState.PlusInt(3)].hmmRightState=lastHmmState.PlusInt(3);

							firstHmmState++;
							firstHmmState++;
							firstHmmState++;
							firstHmmState++;
							break;
					}
					break;
				default:
					assert(false);
					break;
			}

			// do insert states
			CovarianceModel::State firstCmInsertState=firstCmState.PlusInt(thisNodeData.numNormalStates);
			CovarianceModel::State lastCmInsertState=lastCmState;
			CovarianceModel::State cmState;
			for (cmState=firstCmInsertState; cmState!=lastCmInsertState; cmState++) {
				if (sourceCM.GetStateType(cmState)==IL_st) {
					newHMM.hmm.SetStateType(firstHmmState,IL_st);
					newHMM.cm2HmmState[cmState].hmmLeftState=firstHmmState;
					newHMM.cm2HmmState[cmState].hmmRightState=CovarianceModel::GetInvalidState();
					firstHmmState++;
				}
				else {
					assert(sourceCM.GetStateType(cmState)==IR_st);

					lastHmmState--;
					newHMM.hmm.SetStateType(lastHmmState,IL_st);
					newHMM.cm2HmmState[cmState].hmmLeftState=CovarianceModel::GetInvalidState();
					newHMM.cm2HmmState[cmState].hmmRightState=lastHmmState;
				}
			}
		}
	}

	// now figure out child structure
	for (cmNode=firstNode; cmNode!=lastNodePlus1; cmNode++) {

		CovarianceModel::State cmState;
		for (cmState=sourceCM.GetFirstStateOfNode(cmNode); cmState!=sourceCM.GetLastStateOfNode(cmNode); cmState++) {

			if (newHMM.cm2HmmState[cmState].hmmLeftState!=CovarianceModel::GetInvalidState()) {
				SetupChildren (newHMM,sourceCM,cmState,cmNode,HMM_Left,firstNode,lastNode,hmmBuildType);
			}
			if (newHMM.cm2HmmState[cmState].hmmRightState!=CovarianceModel::GetInvalidState()) {
				SetupChildren (newHMM,sourceCM,cmState,cmNode,HMM_Right,firstNode,lastNode,hmmBuildType);
			}
		}
	}

	ReverseMapCm2HmmState(newHMM,sourceCM);
}

void DumpHmm (FILE *file,const HmmAndBuildInfo& newHMM,const CovarianceModel& cm)
{
	newHMM.hmm.DumpInfernalHmm(file,cm);
}

void DumpHmmAndBuildInfo (FILE *file,const HmmAndBuildInfo& newHMM,const CovarianceModel& cm)
{
	DumpHmm(file,newHMM,cm);

	CovarianceModel::State state;
	fprintf(file,"CM had %d states\n",newHMM.cm2HmmState.size());
	for (state=CovarianceModel::IntToState(0); state!=CovarianceModel::IntToState((int)newHMM.cm2HmmState.size()); state++) {
		fprintf(file,"CM state %d maps to HMM: left = %d, right = %d\n",CovarianceModel::StateToInt(state),CovarianceModel::StateToInt(newHMM.cm2HmmState[state].hmmLeftState),CovarianceModel::StateToInt(newHMM.cm2HmmState[state].hmmRightState));
	}
}

CovarianceModel::Node FindFirstEndingNode (const CovarianceModel& sourceCM,CovarianceModel::Node firstNode)
{
	CovarianceModel::Node lastNode=sourceCM.GetLastNode();
	CovarianceModel::Node currNode;
	for (currNode=firstNode; currNode!=lastNode; currNode++) {
		int thisNodeType=sourceCM.GetNodeType(currNode);
		CmNodeTypeData thisNodeData=GetCmNodeTypeData(thisNodeType);

		if (thisNodeData.isEndingNode) {
			return currNode;
		}
	}
	return lastNode;
}

void SetupReverseMapping(ScoreVariablesInfo& scoreVariablesInfo,const HmmAndBuildInfo& newHMM,const CovarianceModel& sourceCM)
{
	TransitionOrEmissionInfo dummy;
	dummy.isUsed=false;
	scoreVariablesInfo.globalVariableToTransitionOrEmissionVector.assign(scoreVariablesInfo.numVariables,dummy);

	CovarianceModel::State hmmState;
	for (hmmState=newHMM.hmm.GetFirstState(); hmmState!=newHMM.hmm.GetLastState(); hmmState++) {

		for (unsigned int i=0;	i<scoreVariablesInfo.transitionToVariableNumVector[hmmState].size(); i++) {
			int globalVar=scoreVariablesInfo.transitionToVariableNumVector[hmmState][i];
			if (globalVar!=-1) {
				CovarianceModel::State toState=newHMM.hmm.GetNthChildState(hmmState,i);

				assert(!scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].isUsed);
				scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].isUsed=true;
				scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].isEmission=false;
				scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].edgeInfo.fromState=hmmState;
				scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].edgeInfo.toState=toState;
				scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].edgeInfo.childNum=(int)i;
			}
		}
		for (unsigned int i=0;	i<scoreVariablesInfo.emissionToVariableNumVector[hmmState].size(); i++) {
			int globalVar=scoreVariablesInfo.emissionToVariableNumVector[hmmState][i];
			if (globalVar!=-1) {

				assert(!scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].isUsed);
				scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].isUsed=true;
				scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].isEmission=true;
				scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].emissionInfo.state=hmmState;
				scoreVariablesInfo.globalVariableToTransitionOrEmissionVector[globalVar].emissionInfo.nuc=(int)i;
			}
		}
	}
}

void SetupTransitionAndEmissionVariables(ScoreVariablesInfo& scoreVariablesInfo,HmmAndBuildInfo& newHMM,const CovarianceModel& sourceCM)
{
	scoreVariablesInfo.transitionToVariableNumVector.resize(newHMM.hmm.GetNumStates());
	scoreVariablesInfo.emissionToVariableNumVector.resize(newHMM.hmm.GetNumStates());

	int nextVariableNum=0;
	CovarianceModel::State hmmState;
	for (hmmState=newHMM.hmm.GetFirstState(); hmmState!=newHMM.hmm.GetLastState(); hmmState++) {

		switch (newHMM.hmm.GetNumChildren(hmmState)) {
			case 0:
				// nothing to do
				break;
			case 1:
				// with only 1 transition, it doesn't matter -- just set it to 0 (this case can arise with PASSTHRU nodes)
				newHMM.hmm.SetTransitionLogScore(hmmState,0,0.0);
				// resize to 1, for convenience
				scoreVariablesInfo.transitionToVariableNumVector[hmmState].resize(1);
				scoreVariablesInfo.transitionToVariableNumVector[hmmState][0]=-1;
				break;
			default:
				// the interesting case -- multiple transitions
				scoreVariablesInfo.transitionToVariableNumVector[hmmState].resize(newHMM.hmm.GetNumChildren(hmmState));
				for (int childNum=0; childNum<newHMM.hmm.GetNumChildren(hmmState); childNum++) {
					if (newHMM.hmm.GetNthChildState(hmmState,childNum)==hmmState) {

						// self-loops are easy -- they must have the same score as the corresponding loop in the CM

						assert(newHMM.hmm.IsInsertState(hmmState)); // this is the only self loops I'm expecting
						assert(newHMM.hmm2CmStateVector[hmmState].size()==1); // insert state should only map to one state
						CovarianceModel::State cmState=newHMM.hmm2CmStateVector[hmmState].front();
						assert(sourceCM.GetNthChildState(cmState,0)==cmState && sourceCM.IsInsertState(cmState)); // this is what I expect based on the stereotyped structure of CMs in the infernal code -- otherwise, I'd have to search for the self-loop

						// NOTE: the following 2 lines are a workaround to a scary (possible) g++ optimizer bug that makes this code crash in release mode, but do fine in debug.  By using the variable 'tsc' (rather than just putting 'sourceCM.GetNthChildTsc(cmState,0)' directly as a param), it doesn't crash, & appears to work.  This happens for RF00032 and RF00016, and I'd guess everything else.  The code appears to work properly with this fix.  I've checked hmm-dump.txt for RF00032 and it achieves the same inflations; the actual weights are different by around 1e-6, but I already know the exact solution is different on different platforms.  Pretty scary, because I'm really not sure what caused this bug, and I can't think of an obvious, benign reason for it to suddenly show up when I changed to code to do structure entirely first, and then solve the scores.  Oh well.
						float tsc=sourceCM.GetNthChildTsc(cmState,0);
						newHMM.hmm.SetTransitionLogScore(hmmState,childNum,tsc);
						scoreVariablesInfo.transitionToVariableNumVector[hmmState][childNum]=-1;

						// while I'm at it, might as well set the emission scores, which should also be the same
						scoreVariablesInfo.emissionToVariableNumVector[hmmState].resize(Alphabet_size);
						for (int nuc=0; nuc<Alphabet_size; nuc++) {
							scoreVariablesInfo.emissionToVariableNumVector[hmmState][nuc]=-1;
							newHMM.hmm.SetSingletEmissionLogScore(hmmState,nuc,sourceCM.GetSingletEmissionScore(cmState,nuc));
						}
					}
					else {
						// add a new variable
						scoreVariablesInfo.transitionToVariableNumVector[hmmState][childNum]=nextVariableNum;
						nextVariableNum++;
					}
				}
				break;
		}

		// we already did IL_st, IR_st, since they're self-looping, and we can just copy the probs from the CM.  Now, we have to do ML_st and MR_st.
		assert(newHMM.hmm.GetNumSymbolsEmitted(hmmState)==0 || newHMM.hmm.GetNumSymbolsEmitted(hmmState)==1); // HMMs don't have MP_st
		if (newHMM.hmm.GetNumSymbolsEmitted(hmmState)==1 && !newHMM.hmm.IsInsertState(hmmState)) {
			scoreVariablesInfo.emissionToVariableNumVector[hmmState].resize(Alphabet_size);
			for (int nuc=0; nuc<Alphabet_size; nuc++) {
				scoreVariablesInfo.emissionToVariableNumVector[hmmState][nuc]=nextVariableNum;
				nextVariableNum++;
			}
		}
	}

	scoreVariablesInfo.numVariables=nextVariableNum;
}

void DumpVariables(FILE *file,const ScoreVariablesInfo& scoreVariablesInfo,const HmmAndBuildInfo& newHMM)
{
	fprintf(file,"\nHMM variables: \n");
	CovarianceModel::State hmmState;
	for (hmmState=newHMM.hmm.GetFirstState(); hmmState!=newHMM.hmm.GetLastState(); hmmState++) {

		fprintf(file,"state %d:\n",CovarianceModel::StateToInt(hmmState));

		fprintf(file,"\tTransitions:\n");
		int childNum;
		for (childNum=0; childNum<(int)scoreVariablesInfo.transitionToVariableNumVector[hmmState].size(); childNum++) {

			if (scoreVariablesInfo.transitionToVariableNumVector[hmmState][childNum]==-1) {
				fprintf(file,"\t\tchild #%d (to state %d): fixed at %f\n",childNum,CovarianceModel::StateToInt(newHMM.hmm.GetNthChildState(hmmState,childNum)),newHMM.hmm.GetNthChildTsc(hmmState,childNum));
			}
			else {
				fprintf(file,"\t\tchild #%d (to state %d): var #%d\n",childNum,CovarianceModel::StateToInt(newHMM.hmm.GetNthChildState(hmmState,childNum)),scoreVariablesInfo.transitionToVariableNumVector[hmmState][childNum]);
			}
		}

		if (newHMM.hmm.GetNumSymbolsEmitted(hmmState)==1) {
			fprintf(file,"\tEmissions:\n");
			for (int nuc=0; nuc<(int)scoreVariablesInfo.emissionToVariableNumVector[hmmState].size(); nuc++) {
				if (scoreVariablesInfo.emissionToVariableNumVector[hmmState][nuc]==-1) {
					fprintf(file,"\t\temit %c : fixed at %f\n",nucs[nuc],newHMM.hmm.GetSingletEmissionScore(hmmState,nuc));
				}
				else {
					fprintf(file,"\t\temit %c : var #%d\n",nucs[nuc],scoreVariablesInfo.emissionToVariableNumVector[hmmState][nuc]);
				}
			}
		}
	}
}

TemporarilyModifyInequality::TemporarilyModifyInequality (Inequality& _inequalitySoFar)
: inequalitySoFar(_inequalitySoFar)
{
	numVariablesAdded=0;
	startingScore=inequalitySoFar.rhs;
	startingWeight=inequalitySoFar.weight;
	starting_sumOfConstantsInHmm=inequalitySoFar.sumOfConstantsInHmm;
	starting_hmmInsertStatesInPath_size=inequalitySoFar.hmmInsertStatesInPath.size();
	for (int i=0; i<MAXABET; i++) {
		startingNucEmitCount[i]=inequalitySoFar.nucEmitCount[i];
	}
}
TemporarilyModifyInequality::~TemporarilyModifyInequality ()
{
	for (int i=0; i<numVariablesAdded; i++) {
		assert(!inequalitySoFar.lhs.empty());
		inequalitySoFar.lhs.pop_back();
	}

	inequalitySoFar.rhs=startingScore;
	inequalitySoFar.weight=startingWeight;
	inequalitySoFar.sumOfConstantsInHmm=starting_sumOfConstantsInHmm;

	assert(inequalitySoFar.hmmInsertStatesInPath.size()>=starting_hmmInsertStatesInPath_size); // how could it get smaller as we do more of the path?!
	while (inequalitySoFar.hmmInsertStatesInPath.size()>starting_hmmInsertStatesInPath_size) {
		inequalitySoFar.hmmInsertStatesInPath.pop_back();
	}

	for (int i=0; i<MAXABET; i++) {
		inequalitySoFar.nucEmitCount[i]=startingNucEmitCount[i];
	}
}
void TemporarilyModifyInequality::PushInsertState (InfernalHmm::State insertState)
{
	inequalitySoFar.hmmInsertStatesInPath.push_back(insertState);
	assert(inequalitySoFar.hmmInsertStatesInPath.size()<=2);
}
void TemporarilyModifyInequality::AddScore (float addToScore)
{
	inequalitySoFar.rhs += addToScore;
}
void TemporarilyModifyInequality::MultiplyWeight (double mult)
{
	inequalitySoFar.weight *= mult;
}
void TemporarilyModifyInequality::AddVariable (int globalVariableNum,std::vector<int>& globalToLocalVariables,int& numLocalVariables)
{
	InequalityTerm newTerm;
	int localVariableNum=globalToLocalVariables[globalVariableNum];
	if (localVariableNum==-1) {
		localVariableNum=numLocalVariables;
		globalToLocalVariables[globalVariableNum]=localVariableNum;
		numLocalVariables++;
	}
	newTerm.variableNum=localVariableNum;
	inequalitySoFar.lhs.push_back(newTerm);
	numVariablesAdded++;
}

void Cm2Hmm_MakeInequalitiesForPath_HmmTransition(TemporarilyModifyInequality& temporarilyModifyInequality,const CovarianceModel::State hmmState,CovarianceModel::State& hmmNextState,const HmmAndBuildInfo& newHMM,const ScoreVariablesInfo& scoreVariablesInfo,std::vector<int>& globalToLocalVariables,int& numLocalVariables,WeightedInequalitiesInfo *weightedInequalitiesInfo)
{
	if (hmmNextState==CovarianceModel::GetInvalidState()) {
		// stays the same, and no transition cost
		hmmNextState=hmmState;
	}
	else {
		if (hmmNextState==hmmState) {
			// nothing to do -- this HMM state didn't change on that CM transition
		}
		else {
			// check if we're going to an insert state
			if (newHMM.hmm.IsInsertState(hmmNextState)) {
				temporarilyModifyInequality.PushInsertState(hmmNextState);
			}

			// update inequality to take into account this transition
			CovarianceModel::State fromState=hmmState,toState=hmmNextState;
			if (fromState>=toState) {
				// we must be on the right side of the HMM, so we're going backwards
				std::swap(fromState,toState);
			}
			int hmmChildNum=newHMM.hmm.GetChildNum_Slow(fromState,toState);
			int globalVariableNum=scoreVariablesInfo.transitionToVariableNumVector[fromState][hmmChildNum];
			if (globalVariableNum==-1) {
				// oh, don't bother
			}
			else {
				// okay, add the variable for this transition
				temporarilyModifyInequality.AddVariable(globalVariableNum,globalToLocalVariables,numLocalVariables);

				// and weight the inequality 
				if (weightedInequalitiesInfo!=NULL) {
					temporarilyModifyInequality.MultiplyWeight(weightedInequalitiesInfo->transitionCounter->GetTransitionFrequency_Unreversed(fromState,toState));
				}
			}
		}
	}
}


void DumblyWorkOutHmmStates(CovarianceModel::State& hmmFirstNormalState,CovarianceModel::State& hmmLastNormalState,const HmmAndBuildInfo& newHMM,CovarianceModel::Node cmNode,const CovarianceModel& sourceCM,const CovarianceModel::State Cm2HmmState::*leftOrRight)
{
	// do this the dumb way that works;

	CmNodeTypeData nodeData=GetCmNodeTypeData(sourceCM.GetNodeType(cmNode));
	CovarianceModel::State cmFirstNormalState=sourceCM.GetFirstStateOfNode(cmNode);
	CovarianceModel::State cmLastNormalState=cmFirstNormalState.PlusInt(nodeData.numNormalStates);

	hmmFirstNormalState=newHMM.cm2HmmState[cmFirstNormalState].*leftOrRight;
	hmmLastNormalState=newHMM.cm2HmmState[cmFirstNormalState].*leftOrRight;

	CovarianceModel::State cmRootState;
	for (cmRootState=cmFirstNormalState; sourceCM.IsStateInRange(cmRootState,cmFirstNormalState,cmLastNormalState); cmRootState++) {
		CovarianceModel::State thisHmmState=newHMM.cm2HmmState[cmRootState].*leftOrRight;
		hmmFirstNormalState=std::min(hmmFirstNormalState,thisHmmState);
		hmmLastNormalState=std::max(hmmLastNormalState,thisHmmState);
	}
	// make half-open
	hmmLastNormalState++;
	assert(hmmFirstNormalState!=CovarianceModel::GetInvalidState() && hmmLastNormalState!=CovarianceModel::GetInvalidState() && hmmFirstNormalState!=CovarianceModel::GetInvalidState() && hmmLastNormalState!=CovarianceModel::GetInvalidState());
	assert(hmmFirstNormalState<hmmLastNormalState); // every node should have at least 1 "normal" state
}

// recursively make the inequalities for each node along a path of consecutive nodes.
void Cm2Hmm_MakeInequalitiesForPath(HmmAndBuildInfo& newHMM,const CovarianceModel& sourceCM,const ScoreVariablesInfo& scoreVariablesInfo,const CovarianceModel::Node cmStartPathNode,const CovarianceModel::Node cmEndPathNode,const CovarianceModel::State cmRootState,const CovarianceModel::State hmmLeftState,const CovarianceModel::State hmmRightState,const bool doneEmission,InequalityList& inequalityList,std::vector<int>& globalToLocalVariables,int& numLocalVariables,Inequality& inequalitySoFar,WeightedInequalitiesInfo *weightedInequalitiesInfo)
{
	//fprintf(stderr,"MakeInequalitiesForPath: root state=%d,done emit=%c\n",CovarianceModel::StateToInt(cmRootState),doneEmission?'T':'F');

	// check for base case
	if (sourceCM.GetNode(cmRootState)==cmEndPathNode) {

		// it's the base case, end of path
		inequalitySoFar.pathEndState=cmRootState;

		/* // nope, not for now at least
		// add to totalWeight (for later normalization)
		totalWeight += inequalitySoFar.weight;
		*/

		if (weightedInequalitiesInfo!=NULL) {
			// do the left side of the HMM in-edges at the end of the walk (search for comment with "JHKOOIUYOIUY" for details on why)
			CovarianceModel::State hmmLeftFirstNormalState,hmmLeftLastNormalState;
			DumblyWorkOutHmmStates(hmmLeftFirstNormalState,hmmLeftLastNormalState,newHMM,cmEndPathNode,sourceCM,&Cm2HmmState::hmmLeftState);
			double leftProb=weightedInequalitiesInfo->transitionCounter->GetEntryProbability_Unreversed(hmmLeftState,hmmLeftFirstNormalState,hmmLeftLastNormalState);
			inequalitySoFar.weight *= leftProb;
		}

		// add this equation to the list
		inequalityList.push_back(inequalitySoFar);
	}
	else {

		// not base case

		// take another step
		if (sourceCM.IsEmitting(cmRootState) && !doneEmission) {

			if (sourceCM.GetNumSymbolsEmitted(cmRootState)==1) {

				// case 1: only 1 symbol is emitted

				// NOTE: we could be a position where cmRootState emits 1 symbol, but both hmmLeftState AND hmmRightState are emitting.  This happens, e.g., in ROOT node in a transition from the IL state to the IR state, where hmmLeftState is still on the IL state.
				CovarianceModel::State hmmEmitState;
				if (sourceCM.EmitsLeft(cmRootState)) {
					assert(newHMM.hmm.IsEmitting(hmmLeftState));
					hmmEmitState=hmmLeftState;
				}
				else {
					assert(sourceCM.EmitsRight(cmRootState));
					assert(newHMM.hmm.IsEmitting(hmmRightState));
					hmmEmitState=hmmRightState;
				}

				// special case: given the current implementation of infernal-0.54, and the fact that we don't need to solve much for IL/IR nodes, we'll have many IL/IR nodes where the emissions don't have variables assigned for any symbol.  In this case, we only need to do the recursion for one of the vars.  This saves time, and may allow us to build more sophisticated sets of equations
				if (sourceCM.IsInsertState(cmRootState)) {
					// is IL/IR: activate special case, but verify it's true

					// sanity checking code
					for (int nuc=0; nuc<Alphabet_size; nuc++) {
						float emitScore=sourceCM.GetSingletEmissionScore(cmRootState,nuc);
						assert (scoreVariablesInfo.emissionToVariableNumVector[hmmEmitState][nuc]==-1); // IL/IR states should never need vars for emissions
						assert(newHMM.hmm.GetSingletEmissionScore(hmmEmitState,nuc)==emitScore); // I should have set it
					}

					// NOTE: don't add anything to inequalitySoFar.sumOfConstantsInHmm or inequalitySoFar.nucEmitCount, since this is an insert state, and insert states (since their emits are not vars, and are equal) get treated weirdly in other code

					// we got here -- it's safe, so recurse
					Cm2Hmm_MakeInequalitiesForPath(newHMM,sourceCM,scoreVariablesInfo,cmStartPathNode,cmEndPathNode,cmRootState,hmmLeftState,hmmRightState,true,inequalityList,globalToLocalVariables,numLocalVariables,inequalitySoFar,weightedInequalitiesInfo);
				}
				else {

					// normal case
					// do the emission
					for (int nuc=0; nuc<Alphabet_size; nuc++) {

						TemporarilyModifyInequality temporarilyModifyInequality(inequalitySoFar);

						// do weighting.  Technically, weighting shouldn't care about whether it's an insert or match (i.e. whether or not we need a variable for the emission cost).  However, we can put this code here: for the insert case, we're really exploring all nucs at once (since we don't use a variable), and the total probability for the possible emissions is 1
						if (weightedInequalitiesInfo!=NULL) {
							temporarilyModifyInequality.MultiplyWeight(weightedInequalitiesInfo->transitionCounter->GetEmitFrequency(hmmEmitState,nuc));
						}

						inequalitySoFar.nucEmitCount[nuc]++;

						float emitScore=sourceCM.GetSingletEmissionScore(cmRootState,nuc);
						if (scoreVariablesInfo.emissionToVariableNumVector[hmmEmitState][nuc]==-1) {
							assert(newHMM.hmm.GetSingletEmissionScore(hmmEmitState,nuc)==emitScore); // I should have set it this way, if no variable was assigned, in which case, they cancel out
							assert(false); // now that I'm doing inequalitySoFar.sumOfConstantsInHmm, all emits except insert states should have variables.  Otherwise, the code has to change to accomodate this.
						}
						else {
							// put it in
							int globalVariableNum=scoreVariablesInfo.emissionToVariableNumVector[hmmEmitState][nuc];
							temporarilyModifyInequality.AddVariable (globalVariableNum,globalToLocalVariables,numLocalVariables);
							temporarilyModifyInequality.AddScore(emitScore);
						}

						// and recurse
						Cm2Hmm_MakeInequalitiesForPath(newHMM,sourceCM,scoreVariablesInfo,cmStartPathNode,cmEndPathNode,cmRootState,hmmLeftState,hmmRightState,true,inequalityList,globalToLocalVariables,numLocalVariables,inequalitySoFar,weightedInequalitiesInfo);
					}
				}
			}
			else {
				assert(sourceCM.GetNumSymbolsEmitted(cmRootState)==2); // only remaining possibility

				for (int leftNuc=0; leftNuc<Alphabet_size; leftNuc++) {
					for (int rightNuc=0; rightNuc<Alphabet_size; rightNuc++) {

						TemporarilyModifyInequality temporarilyModifyInequality(inequalitySoFar);
						float emitScore=sourceCM.GetPairEmissionScore(cmRootState,leftNuc,rightNuc);

						inequalitySoFar.nucEmitCount[leftNuc]++;
						inequalitySoFar.nucEmitCount[rightNuc]++;

						// weight by join prob
						if (weightedInequalitiesInfo!=NULL) {
							temporarilyModifyInequality.MultiplyWeight(weightedInequalitiesInfo->transitionCounter->GetEmitFrequency(hmmLeftState,leftNuc));
							temporarilyModifyInequality.MultiplyWeight(weightedInequalitiesInfo->transitionCounter->GetEmitFrequency(hmmRightState,rightNuc));
						}

						int globalLeftVariableNum=scoreVariablesInfo.emissionToVariableNumVector[hmmLeftState][leftNuc];
						int globalRightVariableNum=scoreVariablesInfo.emissionToVariableNumVector[hmmRightState][rightNuc];
						assert(globalLeftVariableNum!=-1 && globalRightVariableNum!=-1); // MP_st should map to two ML_st in the HMM, and since they're not insert states, we should have variables for both of their emissions

						temporarilyModifyInequality.AddVariable (globalLeftVariableNum,globalToLocalVariables,numLocalVariables);
						temporarilyModifyInequality.AddVariable (globalRightVariableNum,globalToLocalVariables,numLocalVariables);
						temporarilyModifyInequality.AddScore (emitScore);

						// and recurse
						Cm2Hmm_MakeInequalitiesForPath(newHMM,sourceCM,scoreVariablesInfo,cmStartPathNode,cmEndPathNode,cmRootState,hmmLeftState,hmmRightState,true,inequalityList,globalToLocalVariables,numLocalVariables,inequalitySoFar,weightedInequalitiesInfo);
					}
				}
			}
		}
		else {

			// take a transition to another state
			for (int childNum=0; childNum<sourceCM.GetNumChildren(cmRootState); childNum++) {
				CovarianceModel::State cmNextState=sourceCM.GetNthChildState(cmRootState,childNum);
				if (cmNextState!=cmRootState) { // no need to try self loops

					TemporarilyModifyInequality temporarilyModifyInequality(inequalitySoFar);

					float cmTsc=sourceCM.GetNthChildTsc(cmRootState,childNum);
					temporarilyModifyInequality.AddScore(cmTsc);

					CovarianceModel::State hmmNextLeftState=newHMM.cm2HmmState[cmNextState].hmmLeftState;
					CovarianceModel::State hmmNextRightState=newHMM.cm2HmmState[cmNextState].hmmRightState;

					Cm2Hmm_MakeInequalitiesForPath_HmmTransition(temporarilyModifyInequality,hmmLeftState,hmmNextLeftState,newHMM,scoreVariablesInfo,globalToLocalVariables,numLocalVariables,weightedInequalitiesInfo);
					Cm2Hmm_MakeInequalitiesForPath_HmmTransition(temporarilyModifyInequality,hmmRightState,hmmNextRightState,newHMM,scoreVariablesInfo,globalToLocalVariables,numLocalVariables,weightedInequalitiesInfo);

					// and recurse
					Cm2Hmm_MakeInequalitiesForPath(newHMM,sourceCM,scoreVariablesInfo,cmStartPathNode,cmEndPathNode,cmNextState,hmmNextLeftState,hmmNextRightState,false,inequalityList,globalToLocalVariables,numLocalVariables,inequalitySoFar,weightedInequalitiesInfo);
				}
			}
		}
	}
}

void DumpInequalities(FILE *file,const HmmAndBuildInfo& newHMM,const CovarianceModel::Node cmStartPathNode,const CovarianceModel::Node cmEndPathNode,const InequalityList& inequalityList,const std::vector<int>& globalToLocalVariables,int numLocalVariables,const ScoreVariablesInfo& scoreVariablesInfo,const std::vector<float>& localVariableToValue,HmmCommittee::CommitteeBuildInfo *committeeBuildInfo)
{
	// build reverse-mapping info
	// global variables to string description
	std::vector<std::string> globalVarToString;
	globalVarToString.resize(scoreVariablesInfo.numVariables);
	CovarianceModel::State hmmState;
	for (hmmState=newHMM.hmm.GetFirstState(); hmmState!=newHMM.hmm.GetLastState(); hmmState++) {
		for (unsigned int i=0;	i<scoreVariablesInfo.transitionToVariableNumVector[hmmState].size(); i++) {
			int globalVar=scoreVariablesInfo.transitionToVariableNumVector[hmmState][i];
			if (globalVar!=-1) {
				char buf[256];
				CovarianceModel::State toState=newHMM.hmm.GetNthChildState(hmmState,i);
				sprintf(buf,"trans(%s %d -> %s %d)",newHMM.hmm.GetStateTypeName(hmmState),CovarianceModel::StateToInt(hmmState),newHMM.hmm.GetStateTypeName(toState),CovarianceModel::StateToInt(toState));
				globalVarToString[globalVar]=buf;
			}
		}
		for (unsigned int i=0;	i<scoreVariablesInfo.emissionToVariableNumVector[hmmState].size(); i++) {
			int globalVar=scoreVariablesInfo.emissionToVariableNumVector[hmmState][i];
			if (globalVar!=-1) {
				char buf[256];
				CovarianceModel::State toState;
				sprintf(buf,"emit(%s %d , %c)",newHMM.hmm.GetStateTypeName(hmmState),CovarianceModel::StateToInt(hmmState),nucs[i]);
				globalVarToString[globalVar]=buf;
			}
		}
	}
	// local vars to globals
	std::vector<int> localToGlobalVar;
	localToGlobalVar.resize(numLocalVariables);
	for (unsigned int i=0; i<globalToLocalVariables.size(); i++) {
		int localVar=globalToLocalVariables[i];
		if (localVar!=-1) {
			localToGlobalVar[localVar]=i;
		}
	}

	fprintf(file,"\nInequalities for node %d --> node %d\n",CovarianceModel::NodeToInt(cmStartPathNode),CovarianceModel::NodeToInt(cmEndPathNode));

	fprintf(file,"\nglobal to local var mapping:\n");
	for (unsigned int i=0; i<globalToLocalVariables.size(); i++) {
		if (globalToLocalVariables[i]!=-1) {
			fprintf(file,"\tglobal #%d --> local #%d\n",i,globalToLocalVariables[i]);
		}
	}

	fprintf(file,"\ninequalities:\n");
	InequalityList::const_iterator ineqIter;
	for (ineqIter=inequalityList.begin(); ineqIter!=inequalityList.end(); ineqIter++) {

		float totalValue=0;

		// low level
		fprintf(file,"\t");
		std::list<InequalityTerm>::const_iterator lhsIter;
		for (lhsIter=ineqIter->lhs.begin(); lhsIter!=ineqIter->lhs.end(); lhsIter++) {
			if (lhsIter!=ineqIter->lhs.begin()) {
				fprintf(file," + ");
			}
			int variableNum=lhsIter->variableNum;
			fprintf(file,"%d",variableNum);
			float value=localVariableToValue[variableNum];
			totalValue += value;
		}
		fprintf(file," >= %f",ineqIter->rhs);
		//fprintf(file,"  (weight=%lf)",ineqIter->weight);
		fprintf(file,"\n");

		// easier to read
		fprintf(file,"\t");
		for (lhsIter=ineqIter->lhs.begin(); lhsIter!=ineqIter->lhs.end(); lhsIter++) {
			if (lhsIter!=ineqIter->lhs.begin()) {
				fprintf(file," + ");
			}
			int localVar=lhsIter->variableNum;
			int globalVar=localToGlobalVar[localVar];
			std::string description=globalVarToString[globalVar];
			fprintf(file,"%s",description.c_str());
		}
		fprintf(file,"\n");
		fprintf(file,"\t(weight=%lf)\n",ineqIter->weight);
		fprintf(file,"\tachieved value=%f\ninflation=%f\n\n",totalValue,totalValue - ineqIter->rhs);
	}

	if (IsFirstHmmBuilt(committeeBuildInfo) && committeeBuildInfo!=NULL) {
		fprintf(file,"\nSolutions within committee (tab-delimited, for cut&paste into Excel):\n");
		for (int localVar=0; localVar<numLocalVariables; localVar++) {
			int globalVar=localToGlobalVar[localVar];
			std::string localVarDescription=globalVarToString[globalVar];
			fprintf(file,"%s\t",localVarDescription.c_str());

			const HmmCommittee::LinearProgramStoredInfo& storedInfo=committeeBuildInfo->cmStartNodeToLinearProgramInfo
[cmStartPathNode];
			for (size_t solutionNum=0; solutionNum<storedInfo.localVariablesToValuePerSolution.size(); solutionNum++) {
				fprintf(file,"%lg\t",storedInfo.localVariablesToValuePerSolution[solutionNum][localVar]);
			}
			fprintf(file,"\n");
		}
	}

	fprintf(file,"\n");
}

/*
void AddConstraints_set_mat(lpsolve::lprec *linearProgram,const InequalityList& inequalityList,std::vector<double>& weightedUsesOfVariable,const int numLocalVariables)
{
	int result;
	int constraintNum=0+1;
	InequalityList::const_iterator ineqIter;
	for (ineqIter=inequalityList.begin(); ineqIter!=inequalityList.end(); ineqIter++) {
		// set everything to 0 first
		for (int i=0; i<numLocalVariables; i++) {
			result=(int) lpsolve::set_mat(linearProgram,constraintNum,i+1,0.0);
			if (!result) {
				Die("LPSOLVE is unhappy: set_mat failed");
			}
		}
		// now go thru variables in the constraint, & set them to 1
		std::list<InequalityTerm>::const_iterator termIter;
		for (termIter=ineqIter->lhs.begin(); termIter!=ineqIter->lhs.end(); termIter++) {
			weightedUsesOfVariable[termIter->variableNum] += ineqIter->weight;
			result=(int) lpsolve::set_mat(linearProgram,constraintNum,termIter->variableNum+1,1.0);
			if (!result) {
				Die("LPSOLVE is unhappy: set_mat failed");
			}
		}
        result=(int) lpsolve::set_constr_type(linearProgram, constraintNum, GE); // all >= constraints
		if (!result) {
			Die("LPSOLVE is unhappy: set_constr_type failed");
		}
		result=(int) lpsolve::set_rh(linearProgram,constraintNum,ineqIter->rhs);
		if (!result) {
			Die("LPSOLVE is unhappy: set_rh failed");
		}

		constraintNum++;
	}
}
*/

/*
// I wonder if using the add_constraint function is faster
void AddConstraints_add_constraint(lpsolve::lprec *linearProgram,const InequalityList& inequalityList,std::vector<double>& weightedUsesOfVariable,const int numLocalVariables)
{
	std::vector<REAL> coeffs;
	coeffs.resize(numLocalVariables+1);

	int result;
	InequalityList::const_iterator ineqIter;
	for (ineqIter=inequalityList.begin(); ineqIter!=inequalityList.end(); ineqIter++) {
		// set everything to 0 first
		for (int i=0; i<numLocalVariables; i++) {
			coeffs[i+1]=0;
		}
		// now go thru variables in the constraint, & set them to 1
		std::list<InequalityTerm>::const_iterator termIter;
		for (termIter=ineqIter->lhs.begin(); termIter!=ineqIter->lhs.end(); termIter++) {
			weightedUsesOfVariable[termIter->variableNum] += ineqIter->weight;
			coeffs[termIter->variableNum+1]=1;
		}
if (ineqIter->rhs > 0) fprintf(stdout,"WARN: rhs = %f\n",ineqIter->rhs);

		result=(int) lpsolve::add_constraint(linearProgram,&coeffs.front(),GE,ineqIter->rhs);
		if (!result) {
			Die("LPSOLVE is unhappy: add_constraint failed");
		}
	}
}
*/

/*
void SolveInequalities_SetVariableLowerBound(lpsolve::lprec *linearProgram,const int numLocalVariables)
{
	// by default, all variables are constrained to be >0 (I have no idea why, but I'm not about to write my own LP solver...).  Take this off.
	for (int i=0; i<numLocalVariables; i++) {
		// crap -- it doesn't like it when I set it to -infinity
		int result=(int) lpsolve::set_lowbo(linearProgram,i+1,HACK_LOWERBOUND);
		if (!result) {
			Die("LPSOLVE is unhappy: set_lowbo failed");
		}
	}
}
*/

/*
void SolveInequalities_SetObjectiveFunction(lpsolve::lprec *linearProgram,int numLocalVariables,const std::vector<double>& weightedUsesOfVariable)
{
	// now the objective function
	lpsolve::set_minim(linearProgram); // we want to minimize the objective function

	// I think we want the obj func to be the sum of the lhs of all the inequalities -- we really want to minimize the sum of the slack variables -- except now I allow weighting
	for (int i=0; i<numLocalVariables; i++) {
		int result=(int) lpsolve::set_mat(linearProgram,0,i+1,weightedUsesOfVariable[i]);
		if (!result) {
			Die("LPSOLVE is unhappy: set_mat failed");
		}
	}
}
*/

/*
void SolveInequalities_CalcInflation(float& get_avgInflation,lpsolve::lprec *linearProgram,const int numLocalVariables,const std::vector<float>& localVariableToValue,const InequalityList& inequalityList)
{
	REAL *constraints;
	int result=(int) lpsolve::get_ptr_constraints(linearProgram,&constraints);
	if (!result) {
		Die("LPSOLVE is unhappy: 'get_ptr_constraints' failed\n");
	}

	int numRows=(int)(inequalityList.size());
	double totalInflation=0;
	for (int r=0; r<numRows; r++) {
		float lhs=0;
		for (int i=0; i<numLocalVariables; i++) {
			lhs += (float)(lpsolve::get_mat(linearProgram,r+1,i+1) * localVariableToValue[i]);
		}
		float rhs=(float)(lpsolve::get_rh(linearProgram,r+1));
		float inflation=lhs-rhs;
		assert(inflation>=-0.00001);

		double putativeSlack=constraints[r] - lpsolve::get_rh(linearProgram,r+1);
		assert(putativeSlack>=-0.00001);
		double deltaToPutativeSlack=fabs((double)(inflation)-putativeSlack);
		assert(deltaToPutativeSlack/inflation<=0.00001 || deltaToPutativeSlack<=0.00001); // plausible equal to 0, modulo precision errors

		totalInflation += putativeSlack;
	}
	get_avgInflation = (float)(totalInflation)/(float)(numRows);
}
*/

/*
void SolveInequalities_FindVariablesValues(lpsolve::lprec *linearProgram,const int numLocalVariables,std::vector<float>& localVariableToValue)
{
	REAL *variables;
	int result=(int) lpsolve::get_ptr_variables(linearProgram,&variables);
	if (!result) {
		Die("LPSOLVE is unhappy: 'get_ptr_variables' failed");
	}
	localVariableToValue.resize(numLocalVariables);
	for (int i=0; i<numLocalVariables; i++) {
		localVariableToValue[i]=(float)(variables[i]);
	}
}
*/

/*
void SolveInequalities_DumpSolution(lpsolve::lprec *linearProgram,const int numLocalVariables,std::vector<float>& localVariableToValue,HmmCommittee::CommitteeBuildInfo *committeeBuildInfo)
{
#if 0
	lpsolve::print_lp(linearProgram);
	lpsolve::print_solution(linearProgram);
#endif
#ifdef DEBUG_DUMP
	if (IsFirstHmmBuilt(committeeBuildInfo)) {
		fprintf(dumpFile,"\nSolution:\n");
		for (int i=0; i<numLocalVariables; i++) {
			fprintf(dumpFile,"local var #%d = %f\n",i,localVariableToValue[i]);
		}
		fprintf(dumpFile,"Objective func value with solution: %lf\n",lpsolve::get_objective(linearProgram));
		fprintf(dumpFile,"\n");
	}
#endif
}
*/

/*
void SolveInequalities_DumpProgram(lpsolve::lprec *linearProgram,const CovarianceModel::Node cmStartPathNode,const CovarianceModel::Node cmEndPathNode,HmmCommittee::CommitteeBuildInfo *committeeBuildInfo)
{
#ifdef DEBUG_DUMP
	if (IsFirstHmmBuilt(committeeBuildInfo)) {
		fprintf(dumpFile,"\ndumping linear program (node %d --> %d)\n",CovarianceModel::NodeToInt(cmStartPathNode),CovarianceModel::NodeToInt(cmEndPathNode));
		lpsolve::write_LP(linearProgram,dumpFile);
	}
#endif
}
*/

/*
void SolveInequalities_ConstructLP(lpsolve::lprec *(&linearProgram),const InequalityList& inequalityList,const int numLocalVariables,const CovarianceModel::Node cmStartPathNode,const CovarianceModel::Node cmEndPathNode,const int nodesToSpanWhileSolvingScores,HmmCommittee::CommitteeBuildInfo *committeeBuildInfo)
{
	const bool use_add_constraint=true;

	int numRows=(int)(inequalityList.size());
	linearProgram=lpsolve::make_lp(use_add_constraint?0:numRows,numLocalVariables);
	if (linearProgram==NULL) {
		Die("lpsolve couldn't make the LP");
	}

	SolveInequalities_SetVariableLowerBound(linearProgram,numLocalVariables);

	//fprintf(stderr,"Setting up %d constraints on %d variables...\n",numRows,numLocalVariables);
	std::vector<double> weightedUsesOfVariable;
	weightedUsesOfVariable.assign(numLocalVariables,0);
	if (use_add_constraint) {
		AddConstraints_add_constraint(linearProgram,inequalityList,weightedUsesOfVariable,numLocalVariables);
	}
	else {
		AddConstraints_set_mat(linearProgram,inequalityList,weightedUsesOfVariable,numLocalVariables);
	}
	SolveInequalities_SetObjectiveFunction(linearProgram,numLocalVariables,weightedUsesOfVariable);

	SolveInequalities_DumpProgram(linearProgram,cmStartPathNode,cmEndPathNode,committeeBuildInfo);
}
*/

void SolveInequalities(std::vector<float>& localVariableToValue,const InequalityList& inequalityList,const int numLocalVariables,CovarianceModel::Node cmStartPathNode,CovarianceModel::Node cmEndPathNode,float& get_avgInflation,const int nodesToSpanWhileSolvingScores,HmmCommittee::CommitteeBuildInfo *committeeBuildInfo)
{
	/* This first resize is very important, as it actually initializes the vector at this point.
	 * localVariableToValue is a NULL list when received by SolveInequalities()
	 */
	localVariableToValue.resize(numLocalVariables);
	for (int i=0; i<numLocalVariables; i++) {
		localVariableToValue[i] = PARAMETER_GUESS;
	}
}

/*
void SolveInequalities_BeforeCommittees(std::vector<float>& localVariableToValue,const InequalityList& inequalityList,const int numLocalVariables,CovarianceModel::Node cmStartPathNode,CovarianceModel::Node cmEndPathNode,float& get_avgInflation,const int nodesToSpanWhileSolvingScores)
{
	lpsolve::lprec *linearProgram;
	SolveInequalities_ConstructLP(linearProgram,inequalityList,numLocalVariables,cmStartPathNode,cmEndPathNode,nodesToSpanWhileSolvingScores,NULL);

	int result=lpsolve::solve(linearProgram);
	if (result!=OPTIMAL) {
		assert(false);
		Die("LPSOLVE is unhappy: 'solve' failed");
	}

	//fprintf(stderr,"Finding optimal local variables' values.\n");
	SolveInequalities_FindVariablesValues(linearProgram,numLocalVariables,localVariableToValue);

	SolveInequalities_DumpSolution(linearProgram,numLocalVariables,localVariableToValue,NULL);

	SolveInequalities_CalcInflation(get_avgInflation,linearProgram,numLocalVariables,localVariableToValue,inequalityList);

	lpsolve::delete_lp(linearProgram);
	//fprintf(stderr,"Done LP.\n");
}
*/

/*
void SolveInequality_WithCommittee_RegurgitateSolution(std::vector<float>& localVariableToValue,const int numLocalVariables,const CovarianceModel::Node cmStartPathNode,const CovarianceModel::Node cmEndPathNode,float& get_avgInflation,HmmCommittee::CommitteeBuildInfo *committeeBuildInfo)
{
	assert(committeeBuildInfo!=NULL);
	HmmCommittee::LinearProgramStoredInfo& storedInfo=committeeBuildInfo->cmStartNodeToLinearProgramInfo[cmStartPathNode];

	assert(!storedInfo.localVariablesToValuePerSolution.empty());

	unsigned int solutionToUse;
	if (storedInfo.numTimesGeneratedSolution < storedInfo.localVariablesToValuePerSolution.size()) {
		// just use all solutions in sequence
		solutionToUse=storedInfo.numTimesGeneratedSolution;
	}
	else {
		// and then pick a random solution
#ifdef DISABLE_ZRAND
		throw SimpleStringException("Random numbers (zrandlib) are disabled, so this won't work %s:%d",__FILE__,__LINE__);
#else
		solutionToUse=(unsigned int)(RandInt((int)(storedInfo.localVariablesToValuePerSolution.size())));
#endif
	}
	storedInfo.numTimesGeneratedSolution++;

#ifdef DEBUG_DUMP
	fprintf(dumpFile,"For start path node #%d: regurgitate solution #%u\n",CovarianceModel::NodeToInt(cmStartPathNode),solutionToUse);
#endif
	localVariableToValue=storedInfo.localVariablesToValuePerSolution[solutionToUse];
	get_avgInflation=storedInfo.avgInflationPerSolution[solutionToUse];
}
*/

/*
void FindWhichSlackVarsAreZero(_Bvector& slackVariablesFoundWithZero,lpsolve::lprec *linearProgram,size_t numRows,int numLocalVariables)
{
	REAL *constraints;
	int result=(int) lpsolve::get_ptr_constraints(linearProgram,&constraints);
	if (!result) {
		Die("LPSOLVE is unhappy: 'get_ptr_constraints' failed\n");
	}

	for (unsigned int r=0; r<numRows; r++) {
		REAL rhs=lpsolve::get_rh(linearProgram,r+1);
		REAL lhs=constraints[r];
		double slackVar=lhs-rhs;
		assert(slackVar>=-0.000001); // can't be negative

		if (slackVar < 0.000001) { // close enough to zero
			if (!(slackVariablesFoundWithZero[r])) {
#ifdef DEBUG_DUMP
				fprintf(dumpFile,"Slack variable/inequality #%d found to be 0\n",r);
#endif
				slackVariablesFoundWithZero[r]=true;
			}
		}
	}
}
*/

void StoreSolution(HmmCommittee::LinearProgramStoredInfo& storedInfo,const std::vector<float>& localVariableToValue,const float get_avgInflation)
{
	size_t numSolutionsSoFar=storedInfo.localVariablesToValuePerSolution.size();
	assert(numSolutionsSoFar==storedInfo.avgInflationPerSolution.size());

	storedInfo.localVariablesToValuePerSolution.resize(numSolutionsSoFar+1);
	storedInfo.localVariablesToValuePerSolution[numSolutionsSoFar]=localVariableToValue;
	storedInfo.avgInflationPerSolution.resize(numSolutionsSoFar+1);
	storedInfo.avgInflationPerSolution[numSolutionsSoFar]=get_avgInflation;
}

/*
Note: lpsolve does weird things if you re-use the linearProgram (i.e. I was creating one linearProgram instance, and then just modifying individual constraints, and I got weird results for the slack vars).  So, I create a new linearProgram instance from scratch every time I want to call lpsolve::solve.
*/
/*
void SolveInequality_WithCommittee_FindSolutions(const InequalityList& inequalityList,const int numLocalVariables,CovarianceModel::Node cmStartPathNode,CovarianceModel::Node cmEndPathNode,const int nodesToSpanWhileSolvingScores,HmmCommittee::CommitteeBuildInfo *committeeBuildInfo)
{
	assert(committeeBuildInfo!=NULL);
	HmmCommittee::LinearProgramStoredInfo& storedInfo=committeeBuildInfo->cmStartNodeToLinearProgramInfo[cmStartPathNode];
	assert(storedInfo.localVariablesToValuePerSolution.empty()); // we shouldn't have done anything yet

	std::vector<float> localVariableToValue;
	localVariableToValue.resize(numLocalVariables);

	// make the LP
	lpsolve::lprec *linearProgram;
	SolveInequalities_ConstructLP(linearProgram,inequalityList,numLocalVariables,cmStartPathNode,cmEndPathNode,nodesToSpanWhileSolvingScores,committeeBuildInfo);

	// find our first solution
	int result=lpsolve::solve(linearProgram);
	if (result!=OPTIMAL) {
		assert(false);
		Die("LPSOLVE is unhappy: 'solve' failed");
	}
#ifdef DEBUG_DUMP
	fprintf(dumpFile,"First solution for start path node #%d\n",CovarianceModel::NodeToInt(cmStartPathNode));
#endif
	SolveInequalities_FindVariablesValues(linearProgram,numLocalVariables,localVariableToValue);
	SolveInequalities_DumpSolution(linearProgram,numLocalVariables,localVariableToValue,committeeBuildInfo);
	float unconstrainedAvgInflation;
	SolveInequalities_CalcInflation(unconstrainedAvgInflation,linearProgram,numLocalVariables,localVariableToValue,inequalityList);
#ifdef DEBUG_DUMP
	fprintf(dumpFile,"Most general LP gives avg inflation=%f\n",unconstrainedAvgInflation);
#endif

	storedInfo.localVariablesToValuePerSolution.reserve(10); // try to minimize reallocs
	storedInfo.avgInflationPerSolution.reserve(10);

	StoreSolution(storedInfo,localVariableToValue,unconstrainedAvgInflation);

	// now the fun starts, as we try to find other solutions

	// initialize work variables
	_Bvector slackVariablesFoundWithZero; // we stop when we've found a solution with each slack variable 0, or given up
	slackVariablesFoundWithZero.assign(inequalityList.size(),false);
	FindWhichSlackVarsAreZero(slackVariablesFoundWithZero,linearProgram,inequalityList.size(),numLocalVariables);

	lpsolve::delete_lp(linearProgram);

	// go thru the rows (aka inequalities), and see which ones can be zero
	for (size_t currRow=0; currRow<inequalityList.size(); currRow++) {

		if (slackVariablesFoundWithZero[currRow]) {
			// already found a solution with this slack variable being zero (aka this inequality being solved perfectly), don't do anything
		}
		else {
			// haven't solved this inequality perfectly yet -- see if it can be done
			SolveInequalities_ConstructLP(linearProgram,inequalityList,numLocalVariables,cmStartPathNode,cmEndPathNode,nodesToSpanWhileSolvingScores,committeeBuildInfo);

			size_t constraintNum=currRow+1;
			result=(int) lpsolve::set_constr_type(linearProgram,(int)constraintNum,EQ); // require it to be perfect
			if (!result) {
				Die("LPSOLVE is unhappy: set_constr_type failed");
			}

			// and attempt to solve
			result=lpsolve::solve(linearProgram);
			if (result!=OPTIMAL && result!=INFEASIBLE) {
				assert(false);
				Die("LPSOLVE is unhappy: 'solve' failed");
			}

			if (result==OPTIMAL) {

				// it's solvable -- see what the cost is
#ifdef DEBUG_DUMP
				fprintf(dumpFile,"Found solution with ineq #%d perfect\n",currRow);
#endif
				SolveInequalities_FindVariablesValues(linearProgram,numLocalVariables,localVariableToValue);
				SolveInequalities_DumpSolution(linearProgram,numLocalVariables,localVariableToValue,committeeBuildInfo);
				float avgInflation;
				SolveInequalities_CalcInflation(avgInflation,linearProgram,numLocalVariables,localVariableToValue,inequalityList);
#ifdef DEBUG_DUMP
				fprintf(dumpFile,"avg inflation = %f\n",avgInflation);
#endif

				if (fabs(unconstrainedAvgInflation-avgInflation) < 0.000001) {

					// this solution is good
					StoreSolution(storedInfo,localVariableToValue,avgInflation);
					FindWhichSlackVarsAreZero(slackVariablesFoundWithZero,linearProgram,inequalityList.size(),numLocalVariables);
				}
				else {
					// this solution is less optimal, so I'm going to leave it
#ifdef DEBUG_DUMP
					fprintf(dumpFile,"inflation got worse, so not using this alternate solution\n");
#endif
				}
			}
			else {
				// that constraint can't be satisfied perfectly, so just give up on it
#ifdef DEBUG_DUMP
				fprintf(dumpFile,"ineq #%d can't be solved perfectly\n",currRow);
#endif
			}

			// set constraint back to normal
			result=(int) lpsolve::set_constr_type(linearProgram,(int)constraintNum,GE);
			if (!result) {
				Die("LPSOLVE is unhappy: set_constr_type failed");
			}

			lpsolve::delete_lp(linearProgram);
		}
	}

#ifdef DEBUG_DUMP
	fprintf(dumpFile,"Found %u distinct optimal solutions\n",storedInfo.localVariablesToValuePerSolution.size());
#endif
}
*/

/*
void SolveInequality_WithCommittee(std::vector<float>& localVariableToValue,const InequalityList& inequalityList,const int numLocalVariables,CovarianceModel::Node cmStartPathNode,CovarianceModel::Node cmEndPathNode,float& get_avgInflation,const int nodesToSpanWhileSolvingScores,HmmCommittee::CommitteeBuildInfo *committeeBuildInfo)
{
	assert(committeeBuildInfo!=NULL);
	HmmCommittee::LinearProgramStoredInfo& storedInfo=committeeBuildInfo->cmStartNodeToLinearProgramInfo[cmStartPathNode];

	// check if we've already done this
	if (storedInfo.localVariablesToValuePerSolution.empty()) {

		// nope -- we'll have to solve them now
		SolveInequality_WithCommittee_FindSolutions(inequalityList,numLocalVariables,cmStartPathNode,cmEndPathNode,nodesToSpanWhileSolvingScores,committeeBuildInfo);

		// and finally, regurgitate a solution
		storedInfo.numTimesGeneratedSolution=0;
		SolveInequality_WithCommittee_RegurgitateSolution(localVariableToValue,numLocalVariables,cmStartPathNode,cmEndPathNode,get_avgInflation,committeeBuildInfo);
	}
	else {
		// yup -- just find a solution to re-use
		SolveInequality_WithCommittee_RegurgitateSolution(localVariableToValue,numLocalVariables,cmStartPathNode,cmEndPathNode,get_avgInflation,committeeBuildInfo);
	}

}
*/

/* The new SolveInequalities() above replaces this one
void SolveInequalities(std::vector<float>& localVariableToValue,const InequalityList& inequalityList,const int numLocalVariables,CovarianceModel::Node cmStartPathNode,CovarianceModel::Node cmEndPathNode,float& get_avgInflation,const int nodesToSpanWhileSolvingScores,HmmCommittee::CommitteeBuildInfo *committeeBuildInfo)
{
	if (committeeBuildInfo==NULL) {
		SolveInequalities_BeforeCommittees(localVariableToValue,inequalityList,numLocalVariables,cmStartPathNode,cmEndPathNode,get_avgInflation,nodesToSpanWhileSolvingScores);
	}
	else {
		SolveInequality_WithCommittee(localVariableToValue,inequalityList,numLocalVariables,cmStartPathNode,cmEndPathNode,get_avgInflation,nodesToSpanWhileSolvingScores,committeeBuildInfo);
	}
}
*/

void SetHmmScores(HmmAndBuildInfo& newHMM,const ScoreVariablesInfo& scoreVariablesInfo,const std::vector<int>& globalToLocalVariables,const std::vector<float>& localVariablesToValue)
{
	CovarianceModel::State hmmState;
	for (hmmState=newHMM.hmm.GetFirstState(); hmmState!=newHMM.hmm.GetLastState(); hmmState++) {
		for (unsigned int i=0;	i<scoreVariablesInfo.transitionToVariableNumVector[hmmState].size(); i++) {
			int globalVar=scoreVariablesInfo.transitionToVariableNumVector[hmmState][i];
			if (globalVar!=-1) {
				int localVar=globalToLocalVariables[globalVar];
				if (localVar!=-1) {
					int childNum=(int)i;
					newHMM.hmm.SetTransitionLogScore(hmmState,childNum,localVariablesToValue[localVar]);
				}
			}
		}
		for (unsigned int i=0;	i<scoreVariablesInfo.emissionToVariableNumVector[hmmState].size(); i++) {
			int globalVar=scoreVariablesInfo.emissionToVariableNumVector[hmmState][i];
			if (globalVar!=-1) {
				int localVar=globalToLocalVariables[globalVar];
				if (localVar!=-1) {
					int nuc=i;
					newHMM.hmm.SetSingletEmissionLogScore(hmmState,nuc,localVariablesToValue[localVar]);
				}
			}
		}
	}
}

// warning: I think this function takes some short cuts that only work when walking delete state (well, obviously it doesn't do any emissions)
void WalkHmmDeleteStates(const HmmAndBuildInfo& newHMM,const CovarianceModel& cm)
{
	float cmCost=0,hmmCost=0;
	CovarianceModel::State cmState,hmmLeftState,hmmRightState;
	cmState=cm.GetFirstState();
	hmmLeftState=newHMM.cm2HmmState[cmState].hmmLeftState;
	hmmRightState=newHMM.cm2HmmState[cmState].hmmRightState;
	while (cm.GetStateType(cmState)!=E_st && cm.GetStateType(cmState)!=B_st) {

		fprintf(dumpFile,"CM=(%d) $%f   /   HMM=(%d,%d) $%f\n",CovarianceModel::StateToInt(cmState),cmCost,CovarianceModel::StateToInt(hmmLeftState),CovarianceModel::StateToInt(hmmRightState),hmmCost);

		CovarianceModel::State cmOkayChildState=CovarianceModel::GetInvalidState();
		float tscCmCost;
		for (int childNum=0; childNum<cm.GetNumChildren(cmState); childNum++) {
			CovarianceModel::State cmChildState=cm.GetNthChildState(cmState,childNum);
			int cmStateType=cm.GetStateType(cmChildState);
			if (cmStateType==D_st || cmStateType==PASSTHRU_st || cmStateType==E_st || cmStateType==B_st) {
				assert(cmOkayChildState==CovarianceModel::GetInvalidState()); // I think there should be only one of these
				cmOkayChildState=cmChildState;
				tscCmCost=cm.GetNthChildTsc(cmState,childNum);
			}
		}
		assert(cmOkayChildState!=CovarianceModel::GetInvalidState()); // should have got one

		cmState=cmOkayChildState;
		cmCost += tscCmCost;

		CovarianceModel::State hmmNextLeftState=newHMM.cm2HmmState[cmState].hmmLeftState;
		CovarianceModel::State hmmNextRightState=newHMM.cm2HmmState[cmState].hmmRightState;
		if (hmmNextLeftState==CovarianceModel::GetInvalidState()) {
			hmmNextLeftState=hmmLeftState;
		}
		if (hmmNextRightState==CovarianceModel::GetInvalidState()) {
			hmmNextRightState=hmmRightState;
		}
		if (hmmNextLeftState!=hmmLeftState) {
			CovarianceModel::State fromState=hmmLeftState, toState=hmmNextLeftState;
			if (fromState>=toState) {
				std::swap(fromState,toState);
			}
			int childNum=newHMM.hmm.GetChildNum_Slow(fromState,toState);
			hmmCost += newHMM.hmm.GetNthChildTsc(fromState,childNum);
		}
		if (hmmNextRightState!=hmmRightState) {
			CovarianceModel::State fromState=hmmRightState, toState=hmmNextRightState;
			if (fromState>=toState) {
				std::swap(fromState,toState);
			}
			int childNum=newHMM.hmm.GetChildNum_Slow(fromState,toState);
			hmmCost += newHMM.hmm.GetNthChildTsc(fromState,childNum);
		}
		hmmLeftState=hmmNextLeftState;
		hmmRightState=hmmNextRightState;
	}
}

/* // this isn't the way I want to do this
void NormalizeInequalityWeights(InequalityList& inequalityList,WeightedInequalitiesInfo *weightedInequalitiesInfo,double totalWeight)
{
	double pseudocountTotal=weightedInequalitiesInfo->pseudocount * (double)(inequalityList.size());

	InequalityList::iterator ineqIter;
	for (ineqIter=inequalityList.begin(); ineqIter!=inequalityList.end(); ineqIter++) {
		ineqIter->weight = (ineqIter->weight+weightedInequalitiesInfo->pseudocount)/(totalWeight + pseudocountTotal);
	}
}
*/

std::vector<int> ReverseGlobalToLocalVariables(const std::vector<int>& globalToLocalVariables,int numLocalVariables)
{
	std::vector<int> localToGlobalVariables;
	localToGlobalVariables.resize(numLocalVariables);

	int globalVar;
	for (globalVar=0; globalVar<(int)(globalToLocalVariables.size()); globalVar++) {
		int localVar=globalToLocalVariables[globalVar];
		if (localVar!=-1) {
			localToGlobalVariables[localVar]=globalVar;
		}
	}

	return localToGlobalVariables;
}

void Cm2Hmm_SolveScoresForPath(HmmAndBuildInfo& newHMM,const CovarianceModel& sourceCM,const ScoreVariablesInfo& scoreVariablesInfo,CovarianceModel::Node cmStartPathNode,CovarianceModel::Node cmEndPathNode,CovarianceModel::Node firstNode,CovarianceModel::Node lastNode,float& get_avgInflation,const int nodesToSpanWhileSolvingScores,HmmCommittee::CommitteeBuildInfo *committeeBuildInfo,WeightedInequalitiesInfo *weightedInequalitiesInfo,ExtraCm2HmmInfo *extraCm2HmmInfo)
{
	// start in all normal nodes in cmStartPathNode, and end in any normal node in cmEndPathNode
	assert(cmStartPathNode!=cmEndPathNode); // else path is 0 len
	CovarianceModel::State cmFirstStateOfStartPathNode=sourceCM.GetFirstStateOfNode(cmStartPathNode);

	// first initialize mapping of overall variables, to variables we'll use in this Linear Program
	std::vector<int> globalToLocalVariables;
	globalToLocalVariables.assign(scoreVariablesInfo.numVariables,-1); // they all don't have a local mapping
	int numLocalVariables=0;

	// we'll build up a set of inequalities
	InequalityList inequalityList;

	double totalWeight=0;
	CmNodeTypeData startNodeData=GetCmNodeTypeData(sourceCM.GetNodeType(cmStartPathNode));
	CovarianceModel::State cmRootState; // the root of this path thru the CM - try all normal nodes
	CovarianceModel::State cmFirstNormalState=cmFirstStateOfStartPathNode;
	CovarianceModel::State cmLastNormalState=cmFirstStateOfStartPathNode.PlusInt(startNodeData.numNormalStates);
	CovarianceModel::State cmActualLastNormalState=cmLastNormalState.PlusInt(-1);

	// NOTE (JHKOOIUYOIUY): the left side of the HMM goes up in state#, the right side goes down.  To handle this, we check the in-edges for the left side at cmEndPathNode (at the end of walking a path), and the in-edges for the right side at the beginning using cmStartPathNode (i.e. in this function).  Things may also be backwards because TransitionCounter (ScanHMM.h) uses the reversed HMM used for scanning, i.e. HmmType1 (ScanHMM.h), rather than the one that we build in this file using CovarianceModel.
	CovarianceModel::State hmmRightFirstNormalState,hmmRightLastNormalState;
	DumblyWorkOutHmmStates(hmmRightFirstNormalState,hmmRightLastNormalState,newHMM,cmStartPathNode,sourceCM,&Cm2HmmState::hmmRightState);

	for (cmRootState=sourceCM.GetFirstStateOfNode(cmStartPathNode); sourceCM.IsStateInRange(cmRootState,cmFirstNormalState,cmLastNormalState); cmRootState++) {

		// find the HMM states
		CovarianceModel::State hmmLeftState=newHMM.cm2HmmState[cmRootState].hmmLeftState;
		CovarianceModel::State hmmRightState=newHMM.cm2HmmState[cmRootState].hmmRightState;
		assert(hmmLeftState!=CovarianceModel::GetInvalidState() && hmmRightState!=CovarianceModel::GetInvalidState()); // my current code requires that we always get a state, at least for normal states (as opposed to post insert states)

		// prepare root state for search
		Inequality inequalitySoFar;
		inequalitySoFar.lhs.clear();
		inequalitySoFar.rhs=0.0;
		inequalitySoFar.weight=1.0;
		inequalitySoFar.pathStartState=cmRootState;
		inequalitySoFar.sumOfConstantsInHmm=0;
		for (int i=0; i<MAXABET; i++) {
			inequalitySoFar.nucEmitCount[i]=0;
		}
		if (weightedInequalitiesInfo==NULL) {
			// no need to do anything with the weight
		}
		else {
			// start off with in-edge probabilities, on right side only (search for "JHKOOIUYOIUY" for explanation)
			double rightProb=weightedInequalitiesInfo->transitionCounter->GetEntryProbability_Unreversed(hmmRightState,hmmRightFirstNormalState,hmmRightLastNormalState);
			inequalitySoFar.weight *= rightProb;
		}

		// and start exploring paths
		Cm2Hmm_MakeInequalitiesForPath(newHMM,sourceCM,scoreVariablesInfo,cmStartPathNode,cmEndPathNode,cmRootState,hmmLeftState,hmmRightState,false,inequalityList,globalToLocalVariables,numLocalVariables,inequalitySoFar,weightedInequalitiesInfo);
		assert(inequalitySoFar.lhs.empty());
	}

	if (weightedInequalitiesInfo!=NULL) {
		// adjust weights by pseudocounts
		double numSamples=weightedInequalitiesInfo->transitionCounter->GetNumSamples();
		double pseudocount=weightedInequalitiesInfo->pseudocount;
		for (InequalityList::iterator ineqIter=inequalityList.begin(); ineqIter!=inequalityList.end(); ineqIter++) {
			assert(ineqIter->weight>=0 && ineqIter->weight<=1);  // should be a probability of the traversal, and not have any pseudocount-type adjustment
			ineqIter->weight += pseudocount/numSamples;
		}
	}

	if (extraCm2HmmInfo!=NULL) {
		extraCm2HmmInfo->inequalitiesAndLocalVariables[cmStartPathNode].numLocalVariables=numLocalVariables;
		extraCm2HmmInfo->inequalitiesAndLocalVariables[cmStartPathNode].inequalityList=inequalityList;
		extraCm2HmmInfo->inequalitiesAndLocalVariables[cmStartPathNode].globalToLocalVariables=globalToLocalVariables;
		extraCm2HmmInfo->inequalitiesAndLocalVariables[cmStartPathNode].localToGlobalVariables=ReverseGlobalToLocalVariables(globalToLocalVariables,numLocalVariables);
	}

	std::vector<float> localVariableToValue;
	bool actuallySolveScores=true;
	if (extraCm2HmmInfo!=NULL) {
		actuallySolveScores=extraCm2HmmInfo->actuallySolveScores;
	}
	if (actuallySolveScores) {
		//printf("# inequalities: %u, # local vars: %d\n",inequalityList.size(),numLocalVariables);
		SolveInequalities(localVariableToValue,inequalityList,numLocalVariables,cmStartPathNode,cmEndPathNode,get_avgInflation,nodesToSpanWhileSolvingScores,committeeBuildInfo);
	}

#ifdef DEBUG_DUMP
	if (actuallySolveScores) {
		if (IsFirstHmmBuilt(committeeBuildInfo)) {
			DumpInequalities(dumpFile,newHMM,cmStartPathNode,cmEndPathNode,inequalityList,globalToLocalVariables,numLocalVariables,scoreVariablesInfo,localVariableToValue,committeeBuildInfo);
			fprintf(dumpFile,"Avg inflation here = %f\n",get_avgInflation);
		}
	}
#endif

	if (actuallySolveScores) {
		SetHmmScores(newHMM,scoreVariablesInfo,globalToLocalVariables,localVariableToValue);
	}
}

void Cm2Hmm_FindScores (HmmAndBuildInfo& newHMM,const CovarianceModel& sourceCM,const ScoreVariablesInfo& scoreVariablesInfo,CovarianceModel::Node firstNode,CovarianceModel::Node lastNode,const int nodesToSpanWhileSolvingScores,HmmCommittee::CommitteeBuildInfo *committeeBuildInfo,WeightedInequalitiesInfo *weightedInequalitiesInfo,ExtraCm2HmmInfo *extraCm2HmmInfo)
{
	float totalAvgInflation=0;

	CovarianceModel::Node cmStartPathNode,cmEndPathNode;
	cmStartPathNode=firstNode;
	bool keepGoing=true;
	while (keepGoing && cmStartPathNode!=lastNode) {

		//fprintf(stderr,"cmStartPathNode=%d\n",CovarianceModel::NodeToInt(cmStartPathNode));

		cmEndPathNode=cmStartPathNode;
		for (int i=0; i<nodesToSpanWhileSolvingScores; i++) {
			if (cmEndPathNode==lastNode) {
				keepGoing=false; // we're at the end, just do this last segment, and don't re-loop
				break;
			}
			cmEndPathNode++;
		}

		float avgInflation;
		Cm2Hmm_SolveScoresForPath(newHMM,sourceCM,scoreVariablesInfo,cmStartPathNode,cmEndPathNode,firstNode,lastNode,avgInflation,nodesToSpanWhileSolvingScores,committeeBuildInfo,weightedInequalitiesInfo,extraCm2HmmInfo);
		totalAvgInflation += avgInflation;

		cmStartPathNode=cmEndPathNode;
	}

#ifdef DEBUG_DUMP
	fprintf(dumpFile,"total avg inflation for this sub-CM: %f\n",totalAvgInflation);
#endif

	//WalkHmmDeleteStates(newHMM,sourceCM);
}

void IncCm2HmmState(const CovarianceModel& cm,Cm2HmmStateVector& cm2HmmState,CovarianceModel::State firstHmmStateToInc,int incrementBy)
{
	CovarianceModel::State cmState;
	for (cmState=cm.GetFirstState(); cmState!=cm.GetLastState(); cmState++) {
		if (cm2HmmState[cmState].hmmLeftState!=CovarianceModel::GetInvalidState()) {
			if (cm2HmmState[cmState].hmmLeftState>firstHmmStateToInc || firstHmmStateToInc==CovarianceModel::GetInvalidState()) { // sic -- left side only inc's if it's _greater_
				cm2HmmState[cmState].hmmLeftState += incrementBy;
			}
		}
		if (cm2HmmState[cmState].hmmRightState!=CovarianceModel::GetInvalidState()) {
			if (cm2HmmState[cmState].hmmRightState>=firstHmmStateToInc || firstHmmStateToInc==CovarianceModel::GetInvalidState()) {
				cm2HmmState[cmState].hmmRightState += incrementBy;
			}
		}
	}
}

void MergeCm2HmmState(Cm2HmmStateVector& cm2HmmState,const Cm2HmmStateVector& cm2HmmStateToAdd,CovarianceModel::State state,CovarianceModel::State Cm2HmmState::*hmmState)
{
	bool has=cm2HmmState[state].*hmmState!=CovarianceModel::GetInvalidState();
	bool hasToAdd=cm2HmmStateToAdd[state].*hmmState!=CovarianceModel::GetInvalidState();
	if (has && hasToAdd) {
		assert(false); // they both do the same state??
		Die("Internal error %s:%d",__FILE__,__LINE__);
	}
	if (!has && hasToAdd) {
		// get it from ToAdd
		cm2HmmState[state].*hmmState=cm2HmmStateToAdd[state].*hmmState;
	}
}
void MergeCm2HmmState(const CovarianceModel& cm,Cm2HmmStateVector& cm2HmmState,const Cm2HmmStateVector& cm2HmmStateToAdd)
{
	assert(cm2HmmState.size()==cm2HmmStateToAdd.size());

	CovarianceModel::State state;
	for (state=cm.GetFirstState(); state!=cm.GetLastState(); state++) {
		MergeCm2HmmState(cm2HmmState,cm2HmmStateToAdd,state,&Cm2HmmState::hmmLeftState);
		MergeCm2HmmState(cm2HmmState,cm2HmmStateToAdd,state,&Cm2HmmState::hmmRightState);
	}
}

/*
Purpose: build an HMM structurally-only based on a CM, or a part of a CM.  The function calls itself recursively to build
parts of the CM to deal with bifurication nodes.

Params:
firstNode - the CM node of type BEGL, BEGR or ROOT that will correspond to the start state of the HMM
justBuildToBIFNode - if justBuildToBIFNode==false, and we get to a BIF node, the function should recurse and build each component, and
stitch them together.  if justBuildToBIFNode==true, and we get to a BIF node, the we're already recursing, so just build up to the BIF node as if it were an END node.  Regarless of the value of justBuildToBIFNode, if we get to an END node, that's just the end.
*/
void Cm2Hmm_Structurally (HmmAndBuildInfo& newHMM,const CovarianceModel& sourceCM,const CovarianceModel::Node firstNode,ExtraCm2HmmInfo *extraCm2HmmInfo,Cm2Hmm_HmmBuildType hmmBuildType,const bool justBuildToBIFNode=false)
{
	CovarianceModel::Node endingNode=FindFirstEndingNode(sourceCM,firstNode);

	int endingNodeType=sourceCM.GetNodeType(endingNode);
	assert(endingNodeType==END_nd || endingNodeType==BIF_nd);
	if (justBuildToBIFNode || endingNodeType==END_nd) {
		// base case
		Cm2Hmm_Structurally_Block(newHMM,sourceCM,firstNode,endingNode,extraCm2HmmInfo,hmmBuildType);

		SolveScoresPath thisPath;
		thisPath.cmFirstNode=firstNode;
		thisPath.cmEndingNode=endingNode;
		newHMM.solveScoresPathList.push_back(thisPath);

#ifdef DEBUG_DUMP
		DumpHmmAndBuildInfo (dumpFile,newHMM,sourceCM);
#endif
	}
	else {
		// BIF case
		// Build 3 sub-HMMs: one rooted at the current firstNode, and one at each of the child states of the BIF
		// then hook these HMMs together.

		// work out where the subtrees are
		CovarianceModel::State bifState=sourceCM.GetFirstStateOfNode(endingNode);
		if (!sourceCM.IsBifurication(bifState)) {
			Die("BIF nodes should have only a single B_st.  What kind of a world do we live in?");
		}
		CovarianceModel::State leftChildState=sourceCM.GetLeftBifurifactionChild(bifState);
		CovarianceModel::Node leftChildNode=sourceCM.GetNode(leftChildState);
		CovarianceModel::State rightChildState=sourceCM.GetRightBifurifactionChild(bifState);
		CovarianceModel::Node rightChildNode=sourceCM.GetNode(rightChildState);

		// and recurse
		HmmAndBuildInfo rootHmm,leftChildHmm,rightChildHmm;
		Cm2Hmm_Structurally(rootHmm,sourceCM,firstNode,extraCm2HmmInfo,hmmBuildType,true);
		assert(rootHmm.leftToRightPassthruState!=CovarianceModel::GetInvalidState());
		Cm2Hmm_Structurally(leftChildHmm,sourceCM,leftChildNode,extraCm2HmmInfo,hmmBuildType,false);
		Cm2Hmm_Structurally(rightChildHmm,sourceCM,rightChildNode,extraCm2HmmInfo,hmmBuildType,false);

#ifdef DEBUG_DUMP
		fprintf(dumpFile,"Recursed HMM construction.  root @%d,  left child @%d,  right child @%d\n",CovarianceModel::NodeToInt(firstNode),CovarianceModel::NodeToInt(leftChildNode),CovarianceModel::NodeToInt(rightChildNode));
		fprintf(dumpFile,"Starting HMM merge.\n");
		fprintf(dumpFile,"\nroot hmm: \n");
		DumpHmm(dumpFile,rootHmm,sourceCM);
		fprintf(dumpFile,"\nleft child hmm: \n");
		DumpHmm(dumpFile,leftChildHmm,sourceCM);
		fprintf(dumpFile,"\nright child hmm: \n");
		DumpHmm(dumpFile,rightChildHmm,sourceCM);
#endif

		// order is: left part of rootHmm, left child, right child, right part of rootHmm

		// first splice rightChildHmm onto leftChildHmm s.t. the start state of rightChildHmm takes over the end state of leftChildHmm
		// in the cm2HmmState vectors
		IncCm2HmmState(sourceCM,rightChildHmm.cm2HmmState,CovarianceModel::GetInvalidState(),leftChildHmm.hmm.GetNumStates()-1);
		MergeCm2HmmState(sourceCM,leftChildHmm.cm2HmmState,rightChildHmm.cm2HmmState);
		// and in the HMMs
		CovarianceModel::State originalLeftChildEndState=leftChildHmm.hmm.GetActualLastState();
		assert(leftChildHmm.hmm.GetStateType(originalLeftChildEndState)==E_st); // this is what we're clobbering
		leftChildHmm.hmm.AddStates(rightChildHmm.hmm.GetNumStates()-1); // -1 for the E_st we're clobbering
		rightChildHmm.hmm.MoveStatesHigher(rightChildHmm.hmm.GetFirstState(),originalLeftChildEndState);
		leftChildHmm.hmm.CopyStatesVerbatimFrom(rightChildHmm.hmm,originalLeftChildEndState,leftChildHmm.hmm.GetLastState());

		// now splice the curr leftChildHmm (which is a concatenation of the original leftChildHmm and rightChildHmm) into the (middle-ish) part of rootHmm that corresponded to the original E_st in the CM that rootHmm was built from (i.e. rootHmm.leftToRightPassthruState).
		// move the right-side of rootHmm up to make room, and also move the left-to-right passthru state up, because it knows what its children are.
		int numLeftChildStatesThatAreMoving=leftChildHmm.hmm.GetNumStates()-1; // -1 since its E_st is going to be tragically lost in the move
		// in cm2HmmState
		IncCm2HmmState(sourceCM,rootHmm.cm2HmmState,rootHmm.leftToRightPassthruState,numLeftChildStatesThatAreMoving);
		// in the actual HMM
		rootHmm.hmm.MoveStatesHigher(rootHmm.leftToRightPassthruState,rootHmm.leftToRightPassthruState.PlusInt(numLeftChildStatesThatAreMoving)); 
		// prepare leftChildHmm to be spliced in
		// in cm2HmmState
		IncCm2HmmState(sourceCM,leftChildHmm.cm2HmmState,CovarianceModel::GetInvalidState(),CovarianceModel::StateToInt(rootHmm.leftToRightPassthruState));
		// in actual HMM
		leftChildHmm.hmm.MoveStatesHigher(leftChildHmm.hmm.GetFirstState(),rootHmm.leftToRightPassthruState);
		// and copy it in
		rootHmm.hmm.CopyStatesVerbatimFrom(leftChildHmm.hmm,rootHmm.leftToRightPassthruState,rootHmm.leftToRightPassthruState.PlusInt(numLeftChildStatesThatAreMoving));
		MergeCm2HmmState(sourceCM,rootHmm.cm2HmmState,leftChildHmm.cm2HmmState);

		// re-build HMM-to-CM state mapping
		ReverseMapCm2HmmState(rootHmm,sourceCM);

		// merge solveScoresPathList
		rootHmm.solveScoresPathList.insert(rootHmm.solveScoresPathList.end(),leftChildHmm.solveScoresPathList.begin(),leftChildHmm.solveScoresPathList.end());
		rootHmm.solveScoresPathList.insert(rootHmm.solveScoresPathList.end(),rightChildHmm.solveScoresPathList.begin(),rightChildHmm.solveScoresPathList.end());

		// et voila -- copy it back
		newHMM.hmm.CopyFrom(rootHmm.hmm);
		newHMM.cm2HmmState=rootHmm.cm2HmmState;
		newHMM.hmm2CmStateVector=rootHmm.hmm2CmStateVector;
		newHMM.solveScoresPathList=rootHmm.solveScoresPathList;
		newHMM.leftToRightPassthruState=CovarianceModel::GetInvalidState(); // we shouldn't need this any more, and don't have a meaningful answer anyway
#ifdef DEBUG_DUMP
		fprintf(dumpFile,"\nResulting hmm:\n");
		DumpHmm(dumpFile,newHMM,sourceCM);
#endif
	}
}

void SetupLeftRightness(HmmAndBuildInfo& newHMM,const CovarianceModel& sourceCM)
{
	// set left/right-ness of HMM states
	InfernalHmm::State hmmState;
	for (hmmState=newHMM.hmm.GetFirstState(); hmmState!=newHMM.hmm.GetLastState(); hmmState++) {
		newHMM.hmm.SetLeftState(hmmState,false);
		newHMM.hmm.SetRightState(hmmState,false);
	}
	CovarianceModel::State cmState;
	for (cmState=sourceCM.GetFirstState(); cmState!=sourceCM.GetLastState(); cmState++) {

		InfernalHmm::State hmmLeftState=newHMM.cm2HmmState[cmState].hmmLeftState;
		if (hmmLeftState!=InfernalHmm::GetInvalidState()) {
			newHMM.hmm.SetLeftState(hmmLeftState,true);
		}

		InfernalHmm::State hmmRightState=newHMM.cm2HmmState[cmState].hmmRightState;
		if (hmmRightState!=InfernalHmm::GetInvalidState()) {
			newHMM.hmm.SetRightState(hmmRightState,true);
		}
	}
}

void Cm2Hmm_Structurally(HmmAndBuildInfo& newHMM,const CovarianceModel& sourceCM,ExtraCm2HmmInfo *extraCm2HmmInfo,Cm2Hmm_HmmBuildType hmmBuildType)
{
	CovarianceModel::Node cmFirstNodeForConversion=sourceCM.GetFirstNode();
	if (sourceCM.DoLocal()) {
		cmFirstNodeForConversion++;  // don't do anything for the whole first node, since that gets clobbered by the ConfigLocal function in modelconfig.c.
	}
	Cm2Hmm_Structurally(newHMM,sourceCM,cmFirstNodeForConversion,extraCm2HmmInfo,hmmBuildType);

	newHMM.hmm.AllocHmmData();

	SetupLeftRightness(newHMM,sourceCM);
}

void SetLocal(HmmAndBuildInfo& newHMM,const CovarianceModel& sourceCM)
{
	assert(sourceCM.DoLocal()); // otherwise there's no reason to call this function

	newHMM.hmm.SetDoLocal(true);

	InfernalHmm::State hmmState;
	for (hmmState=newHMM.hmm.GetFirstState(); hmmState!=newHMM.hmm.GetLastState(); hmmState++) {
		newHMM.hmm.SetLeftwardBeginsc(hmmState,(float)IMPOSSIBLE);
		newHMM.hmm.SetRightwardBeginsc(hmmState,(float)IMPOSSIBLE);
		newHMM.hmm.SetNumEndscLinksToLeft(hmmState,0);
	}

	CovarianceModel::State cmState;
	for (cmState=sourceCM.GetFirstState(); cmState!=sourceCM.GetLastState(); cmState++) {

		if (sourceCM.GetNode(cmState)!=sourceCM.GetFirstNode()) { // don't do anything for the whole first node, since that gets clobbered by the ConfigLocal function in modelconfig.c. 
			InfernalHmm::State hmmLeftState=newHMM.cm2HmmState[cmState].hmmLeftState;
			InfernalHmm::State hmmRightState=newHMM.cm2HmmState[cmState].hmmRightState;

			float beginsc=sourceCM.GetBeginsc(cmState);
			if (beginsc!=(float)IMPOSSIBLE) {
				

				if (hmmLeftState==InfernalHmm::GetInvalidState() || hmmRightState==InfernalHmm::GetInvalidState()) {
					assert(false);
					Die("internal error %s:%d",__FILE__,__LINE__);
				}

				float hmmOldLeftBeginsc=newHMM.hmm.GetLeftwardBeginsc(hmmLeftState);
				float hmmOldRightBeginsc=newHMM.hmm.GetRightwardBeginsc(hmmRightState);

				// scores must sum to original beginsc (or more), for now I'll just divide the score evenly into 2.
				float newBeginsc=beginsc/(float)2.0;

				// sometimes multiple CM states can use the same HMM left or right state; so, we must take the max of the beginsc
				newHMM.hmm.SetLeftwardBeginsc(hmmLeftState,std::max(hmmOldLeftBeginsc,newBeginsc));
				newHMM.hmm.SetRightwardBeginsc(hmmRightState,std::max(hmmOldRightBeginsc,newBeginsc));
			}

			if (sourceCM.GetEndsc(cmState)!=(float)IMPOSSIBLE) {

				// link right state to this left state, for simulating endsc
				if (hmmLeftState==InfernalHmm::GetInvalidState() || hmmRightState==InfernalHmm::GetInvalidState()) {
					assert(false);
					Die("For making an HMM local, it's pretty tricky if the insert states (IL,IR) in the CM have non-IMPOSSIBLE endsc scores, because in my mapping, IL,IR states don't have both a left&right hmmState, and it's not easy to see what it should be.  So, I just assumed it wouldn't happen, making an ass out of me and you, the anonymous user of this program.  Sorry.");
				}

				int numLinks=newHMM.hmm.GetNumEndscLinksToLeft(hmmRightState);
				newHMM.hmm.SetNumEndscLinksToLeft(hmmRightState,numLinks+1);
				newHMM.hmm.SetEndscLinkToLeft_State(hmmRightState,numLinks,hmmLeftState);
				newHMM.hmm.SetEndscLinkToLeft_Endsc(hmmRightState,numLinks,sourceCM.GetEndsc(cmState));
			}
		}
	}
}

void Cm2Hmm (HmmAndBuildInfo& newHMM,Cm2Hmm_HmmBuildType hmmBuildType,const CovarianceModel& sourceCM,const int nodesToSpanWhileSolvingScores,HmmCommittee::CommitteeBuildInfo *committeeBuildInfo=NULL,WeightedInequalitiesInfo *weightedInequalitiesInfo=NULL,ExtraCm2HmmInfo *extraCm2HmmInfo=NULL)
{
	Cm2Hmm_Structurally(newHMM,sourceCM,extraCm2HmmInfo,hmmBuildType);
	newHMM.hmm.SetHmm2CmState(newHMM.hmm2CmStateVector);
	newHMM.hmm.SetCm2HmmState(newHMM.cm2HmmState);

	ScoreVariablesInfo scoreVariablesInfo;
	SetupTransitionAndEmissionVariables(scoreVariablesInfo,newHMM,sourceCM);
	SetupReverseMapping(scoreVariablesInfo,newHMM,sourceCM);
#ifdef DEBUG_DUMP
	DumpVariables(dumpFile,scoreVariablesInfo,newHMM);
#endif
	if (extraCm2HmmInfo!=NULL) {
		extraCm2HmmInfo->scoreVariablesInfo=scoreVariablesInfo;
		InequalitiesAndLocalVariables dummyInequalitiesAndLocalVariables;
		dummyInequalitiesAndLocalVariables.numLocalVariables=0;
		extraCm2HmmInfo->inequalitiesAndLocalVariables.assign(sourceCM.GetNumNodes(),dummyInequalitiesAndLocalVariables);
	}

	SolveScoresPathList::const_iterator pathIter;
	for (pathIter=newHMM.solveScoresPathList.begin(); pathIter!=newHMM.solveScoresPathList.end(); pathIter++) {
		Cm2Hmm_FindScores(newHMM,sourceCM,scoreVariablesInfo,pathIter->cmFirstNode,pathIter->cmEndingNode,nodesToSpanWhileSolvingScores,committeeBuildInfo,weightedInequalitiesInfo,extraCm2HmmInfo);
	}

	if (sourceCM.DoLocal()) {
		SetLocal(newHMM,sourceCM);
	}
}


std::string EscapeCmFileName(const std::string& cmFileName)
{
	std::string s(cmFileName);
	for (unsigned int i=0; i<s.size(); i++) {
		if (s[i]=='/') {
			s[i]='_';
		}
	}

	return s;
}

std::string MakeCm2HmmCacheFileName(const char *cmFileName,const bool doLocalAlignment,const int nodesToSpanWhileSolvingScores,const int committeeSize,Cm2Hmm_HmmBuildType hmmBuildType)
{
	std::string cacheFileName("cache/cm2hmm_");
	std::string escapeCmFileName=EscapeCmFileName(cmFileName);
	cacheFileName += escapeCmFileName;
	if (doLocalAlignment) {
		cacheFileName += "_local";
	}
	if (nodesToSpanWhileSolvingScores>1) {
		char buf[16];
		sprintf(buf,"_spl%d",nodesToSpanWhileSolvingScores); // "spl"="Sub-Path Length"
		cacheFileName += buf;
	}
	if (committeeSize>1) {
		char buf[16];
		sprintf(buf,"_cmtee%d",committeeSize);
		cacheFileName += buf;
	}
	{
		char buf[32];
		sprintf(buf,"_buildtype%d",(int)hmmBuildType);
		cacheFileName += buf;
	}
	return cacheFileName;
}

void Cm2Hmm_WithWeighting_NoCaching (InfernalHmm& hmm,Cm2Hmm_HmmBuildType hmmBuildType,const CovarianceModel& sourceCM,const char *cmFileName,WeightedInequalitiesInfo *weightedInequalitiesInfo,ExtraCm2HmmInfo *extraCm2HmmInfo,const int nodesToSpanWhileSolvingScores)
{
	VerifyValidHmmBuildType(hmmBuildType);

	// verify that the structure of the CM is what we expect, so we can just assume it for the rest of this function
	if (!ValidateCmNodeTypeData(sourceCM)) {
		Die("(while converting Covariance Model to HMM) CM doesn't have the structure we expect");
	}

	HmmAndBuildInfo newHMM;

	Cm2Hmm(newHMM,hmmBuildType,sourceCM,nodesToSpanWhileSolvingScores,NULL,weightedInequalitiesInfo,extraCm2HmmInfo);
	newHMM.hmm.SetFromCmFileName(cmFileName);

#ifdef DEBUG_DUMP
	DumpHmmAndBuildInfo (dumpFile,newHMM,sourceCM);
#endif

	hmm.CopyFrom(newHMM.hmm);
}

void Cm2Hmm (InfernalHmm& hmm,Cm2Hmm_HmmBuildType hmmBuildType,const CovarianceModel& sourceCM,const char *cmFileName,const std::string& programParams,bool forceCreate)
{
	const int nodesToSpanWhileSolvingScores=1; // making it higher didn't seem to lead to any more optimal scores (using RF00032==Histone3, which is small), and it takes forever to build the HMM -- it spends most of its time entering the constraints to lp_solve, since lp_solve's sparse matrix format requires it to do a lot of copying & reallocation.  I could probably optimize this by eliminating redundant constraints myself, but it doesn't look like we can improve things this way anyway.

	if (dumpFile==NULL) {
		dumpFile=ThrowingFopen("hmm-dump.txt","wt");
	}

	InfernalHmm creatingHmm;
#ifdef ENABLE_CACHING
	std::string cacheFileName;
	if (overrideHmmCacheFileName!=NULL) {
		cacheFileName=overrideHmmCacheFileName;
	}
	else {
		cacheFileName=MakeCm2HmmCacheFileName(cmFileName,sourceCM.DoLocal(),nodesToSpanWhileSolvingScores,1,hmmBuildType);
	}
	fprintf(stderr,"Trying to load HMM from cached file '%s'\n",cacheFileName.c_str());
	printf("Trying to load HMM from cached file '%s'\n",cacheFileName.c_str());

	bool loadHmmFromCache=!forceCreate;
	if (loadHmmFromCache) {
		loadHmmFromCache=creatingHmm.LoadInBinary(cacheFileName.c_str());
	}
	if (loadHmmFromCache) {
#ifdef DEBUG_DUMP
		fprintf(dumpFile,"Loaded HMM from cache.\n");
		fprintf(stderr,"Loaded HMM from cache.\n");
#ifdef DEBUG_DUMP
	creatingHmm.DumpInfernalHmm(dumpFile,sourceCM);
#endif
#endif
	}
	else {
#endif
	fprintf(stderr,"Building HMM based on CM...\n");

	Cm2Hmm_WithWeighting_NoCaching(creatingHmm,hmmBuildType,sourceCM,cmFileName,NULL,NULL,nodesToSpanWhileSolvingScores);
	std::string build("(default construction of HMM) ");
	build += programParams;
	creatingHmm.AddBuildDescription(build);
	creatingHmm.SetHmmBuildType(hmmBuildType);

	fprintf(stderr,"Built HMM.\n");

#ifdef ENABLE_CACHING
		
		creatingHmm.SaveInBinary(cacheFileName.c_str());
	}
#endif

	hmm.CopyFrom(creatingHmm);

	fflush(dumpFile);
}


HmmCommittee::HmmCommittee ()
{
	recommendedCommitteeSize=0;
}
HmmCommittee::~HmmCommittee ()
{
}
int HmmCommittee::GetRecommendedCommitteeSize (void) const
{
	if (recommendedCommitteeSize==0) {
		Die("HmmCommittee::GetRecommendedCommitteeSize called, but we don't have a good effective committee size -- you must explicitly give one.");
	}
	return recommendedCommitteeSize;
}
bool HmmCommittee::LoadInBinary (const char *fileName,int effectiveCommitteeSize)
{
	FILE *file=fopen(fileName,"rb");
	if (file==NULL) {
		return false;
	}

	size_t committeeSize;
	fread(&committeeSize,sizeof(committeeSize),1,file);

	if (committeeSize==0) {

		int format;
		fread(&format,sizeof(format),1,file);
		if (format!=0) {
			Die("Expecting HMM committee to be saved in file format #0\n");
		}

		size_t descriptionLen;
		fread(&descriptionLen,sizeof(descriptionLen),1,file);
		char *descriptionTmp=new char [descriptionLen+1];
		fread(descriptionTmp,1,descriptionLen,file);
		descriptionTmp[descriptionLen]=0;
		description=descriptionTmp;
		delete [] descriptionTmp;

		fprintf(stderr,"----hmmCommitteeOrigin: %s\n",description.c_str());
		printf("----hmmCommitteeOrigin: %s\n",description.c_str());

		fread(&recommendedCommitteeSize,sizeof(recommendedCommitteeSize),1,file);

		// remember to read the committee size
		fread(&committeeSize,sizeof(committeeSize),1,file);
	}
	else {
		fprintf(stderr,"----hmmCommitteeOrigin: Old file format -- unknown origin\n");
		printf("----hmmCommitteeOrigin: Old file format -- unknown origin\n");
		recommendedCommitteeSize=0;
	}
	if (effectiveCommitteeSize==0) {
		effectiveCommitteeSize=GetRecommendedCommitteeSize();
	}

	printf("Loading %d HMMs (effective committee size) out of max %u\n",effectiveCommitteeSize,committeeSize);
	committeeSize=std::min(committeeSize,(size_t)effectiveCommitteeSize);
	for (size_t i=0; i<committeeSize; i++) {
		InfernalHmm dummy;
		hmmList.push_back(dummy);
		hmmList.back().LoadInBinary(file);
	}
	return true;
}
bool HmmCommittee::SaveInBinary (const char *fileName,const char *additionalDescription)
{
	FILE *file=fopen(fileName,"wb");
	if (file==NULL) {
		return false;
	}

	description += additionalDescription;

	size_t zero=0; // signify new file format
	fwrite(&zero,sizeof(zero),1,file);

	int format=0;
	fwrite(&format,sizeof(format),1,file);

	// write description string
	size_t descriptionLen=description.size();
	fwrite(&descriptionLen,sizeof(descriptionLen),1,file);
	fwrite(description.c_str(),1,description.size(),file);

	// and the actual committee

	fwrite(&recommendedCommitteeSize,sizeof(recommendedCommitteeSize),1,file);
	size_t committeeSize=hmmList.size();
	fwrite(&committeeSize,sizeof(committeeSize),1,file);

	for (HmmList::iterator i=hmmList.begin(); i!=hmmList.end(); i++) {
		i->SaveInBinary(file);
	}
	return true;
}
void HmmCommittee::DumpMemberBuildInfo (FILE *out)
{
	if (disableHmmBuildInfoDump) {
		// don't dump anything
		return;
	}

	fprintf(out,"Dumping build info of committee:\n");
	int n=0;
	for (HmmList::iterator i=hmmList.begin(); i!=hmmList.end(); i++) {
		fprintf(out,"\tHMM #%d: %s\n",n,i->GetBuildDescription().c_str());
		n++;
	}
}

void Cm2Hmm (HmmCommittee& hmmCommittee,const CovarianceModel& sourceCM,const char *cmFileName,int committeeSize,int effectiveCommitteeSize)
{
	dumpFile=fopen("hmm-dump-committee.txt","wt");

	const int nodesToSpanWhileSolvingScores=1; // see single-building HMM function for details

#ifdef ENABLE_CACHING
	std::string cacheFileName;
	if (overrideHmmCacheFileName!=NULL) {
		cacheFileName=overrideHmmCacheFileName;
	}
	else {
		cacheFileName=MakeCm2HmmCacheFileName(cmFileName,sourceCM.DoLocal(),nodesToSpanWhileSolvingScores,committeeSize,HmmBuildType_Original);
	}
	fprintf(stderr,"Trying to load HMM Committee from cached file '%s'\n",cacheFileName.c_str());
	printf("Trying to load HMM Committee from cached file '%s'\n",cacheFileName.c_str());

	if (hmmCommittee.LoadInBinary(cacheFileName.c_str(),effectiveCommitteeSize)) {
#ifdef DEBUG_DUMP
		fprintf(dumpFile,"Loaded HMM Committee from cache.\n");
		fprintf(stderr,"Loaded HMM Committee from cache.\n");
#endif
	}
	else {
#endif
	fprintf(stderr,"Building HMM Committee of size %d based on CM...\n",committeeSize);

	// verify that the structure of the CM is what we expect, so we can just assume it for the rest of this function
	if (!ValidateCmNodeTypeData(sourceCM)) {
		Die("(while converting Covariance Model to HMM) CM doesn't have the structure we expect");
	}

	HmmCommittee::CommitteeBuildInfo committeeBuildInfo;
	committeeBuildInfo.cmStartNodeToLinearProgramInfo.resize(sourceCM.GetNumNodes());
	for (int i=0; i<committeeSize; i++) {
		committeeBuildInfo.currCommitteeMemberNum=i;

		HmmAndBuildInfo newHMM;
		Cm2Hmm(newHMM,HmmBuildType_Original,sourceCM,nodesToSpanWhileSolvingScores,&committeeBuildInfo);

		InfernalHmm dummy;
		hmmCommittee.hmmList.push_back(dummy);
		hmmCommittee.hmmList.back().CopyFrom(newHMM.hmm);
	}

	fprintf(stderr,"Built HMM Committee.\n");

#ifdef ENABLE_CACHING
		
		hmmCommittee.SaveInBinary(cacheFileName.c_str(),cmFileName);
	}
#endif

#ifdef DEBUG_DUMP
	fprintf(dumpFile,"\n\nDUMPING HMM COMMITTEE\n\n");
	int c=0;
	for (HmmCommittee::HmmList::const_iterator i=hmmCommittee.hmmList.begin(); i!=hmmCommittee.hmmList.end(); i++) {
		fprintf(dumpFile,"Commitee member #%d\n",c);
		HmmAndBuildInfo hmmAndBuildInfo;
		hmmAndBuildInfo.hmm.CopyFrom(*i);
		DumpHmmAndBuildInfo (dumpFile,hmmAndBuildInfo,sourceCM);
		fprintf(dumpFile,"\n");
		c++;
	}
#endif

#ifdef DEBUG_DUMP
	fclose(dumpFile);
#endif
}

void GetOrSet_LocalVariablesValueForNode(std::vector<float>& localVariablesToValue,const CovarianceModel& cm,InfernalHmm& infernalHmm,const CovarianceModel::Node cmNode,const ExtraCm2HmmInfo& extraCm2HmmInfo,std::list<int>& globalVarNums,bool isSetting)
{
	if (!isSetting) {
		localVariablesToValue.assign(extraCm2HmmInfo.inequalitiesAndLocalVariables[cmNode].numLocalVariables,FLT_MAX); // MAX_FLT doesn't seem reasonable -- protect against getting something twice
	}
	
	assert(extraCm2HmmInfo.inequalitiesAndLocalVariables[cmNode].numLocalVariables>0); // else you probably shouldn't be calling

	// first emissions
	InfernalHmm::StateFromNodeList hmmStates;
	infernalHmm.GetHmmStatesByCmNode (hmmStates,cm,cmNode);

	for (InfernalHmm::StateFromNodeList::iterator i=hmmStates.begin(); i!=hmmStates.end(); i++) {
		InfernalHmm::State state=i->state;
		if (infernalHmm.IsEmitting(state)) {
			size_t s=extraCm2HmmInfo.scoreVariablesInfo.emissionToVariableNumVector.size(); // for aid in debugger
			if (extraCm2HmmInfo.scoreVariablesInfo.emissionToVariableNumVector[state].size() == 0) {
				// although this state emits, its emission probs weren't configured as variables, because they could reasonably be solved directly.  This is the case for insert states.
			}
			else {
				for (int nuc=0; nuc<Alphabet_size; nuc++) {
					int globalVar=extraCm2HmmInfo.scoreVariablesInfo.emissionToVariableNumVector[state][nuc];
					if (globalVar==-1) {
						// wasn't configured as variable
					}
					else {
						globalVarNums.push_back(globalVar);

						int localVar=extraCm2HmmInfo.inequalitiesAndLocalVariables[cmNode].globalToLocalVariables[globalVar];
						assert(localVar!=-1);

						if (isSetting) {
							infernalHmm.SetSingletEmissionLogScore(state,nuc,localVariablesToValue[localVar]);
						}
						else {
							assert(localVariablesToValue[localVar]==FLT_MAX); // hasn't been set before
							localVariablesToValue[localVar]=infernalHmm.GetSingletEmissionScore(state,nuc);
						}
					}
				}
			}
		}
	}

	// transitions
	InfernalHmm::EdgeInfoList edgeInfoList;
	infernalHmm.GetHmmEdgesByCmNode (edgeInfoList,hmmStates,cm,cmNode);

	for (InfernalHmm::EdgeInfoList::iterator i=edgeInfoList.begin(); i!=edgeInfoList.end(); i++) {
		InfernalHmm::State fromState=i->fromState;
		InfernalHmm::State toState=i->toState;
		int childNum=i->childNum;
		if (extraCm2HmmInfo.scoreVariablesInfo.transitionToVariableNumVector[i->fromState].size() == 0) {
			// this transition wasn't configured as a variable, which is done for the dummy PASSTHRU states
		}
		else {
			int globalVar=extraCm2HmmInfo.scoreVariablesInfo.transitionToVariableNumVector[i->fromState][i->childNum];
			if (globalVar==-1) {
				// wasn't configured as variable
			}
			else {

				globalVarNums.push_back(globalVar);

				/* // was for aid in debugger
				CovarianceModel::Node otherNode;
				for (otherNode=cm.GetFirstNode(); otherNode!=cm.GetLastNode(); otherNode++) {
					if (extraCm2HmmInfo.inequalitiesAndLocalVariables[cmNode].globalToLocalVariables[globalVar]!=-1) {
						break;
					}
				}
				*/
				int localVar=extraCm2HmmInfo.inequalitiesAndLocalVariables[cmNode].globalToLocalVariables[globalVar];
				assert(localVar!=-1);

				if (isSetting) {
					infernalHmm.SetTransitionLogScore(i->fromState,i->childNum,localVariablesToValue[localVar]);
				}
				else {
					assert(localVariablesToValue[localVar]==FLT_MAX); // hasn't been set before
					localVariablesToValue[localVar]=infernalHmm.GetNthChildTsc(i->fromState,i->childNum);
				}
			}
		}
	}
}

void GetGlobalVariableNumsForNode (std::list<int>& globalVarNums,const CovarianceModel& cm,InfernalHmm& infernalHmm,const CovarianceModel::Node cmNode,const ExtraCm2HmmInfo& extraCm2HmmInfo)
{
	std::vector<float> localVariablesToValue;
	GetOrSet_LocalVariablesValueForNode(localVariablesToValue,cm,infernalHmm,cmNode,extraCm2HmmInfo,globalVarNums,false);
}

void GetLocalVariablesValueForNode(std::vector<float>& localVariablesToValue,const CovarianceModel& cm,InfernalHmm& infernalHmm,const CovarianceModel::Node cmNode,const ExtraCm2HmmInfo& extraCm2HmmInfo)
{
	std::list<int> dummy;
	GetOrSet_LocalVariablesValueForNode(localVariablesToValue,cm,infernalHmm,cmNode,extraCm2HmmInfo,dummy,false);
}

void SetLocalVariablesValueForNode(std::vector<float>& localVariablesToValue,const CovarianceModel& cm,InfernalHmm& infernalHmm,const CovarianceModel::Node cmNode,const ExtraCm2HmmInfo& extraCm2HmmInfo)
{
	std::list<int> dummy;
	GetOrSet_LocalVariablesValueForNode(localVariablesToValue,cm,infernalHmm,cmNode,extraCm2HmmInfo,dummy,true);
}


void GetGlobalVarsFromInfernalHmm (std::vector<double>& globalVars,const InfernalHmm& infernalHmm,const TransitionOrEmissionInfoVector& transitionOrEmissionInfoVector)
{
	size_t numGlobalVars=transitionOrEmissionInfoVector.size();
	globalVars.resize(numGlobalVars);

	for (size_t i=0; i<numGlobalVars; i++) {
		assert(transitionOrEmissionInfoVector[i].isUsed);
		if (transitionOrEmissionInfoVector[i].isEmission) {
			globalVars[i]=infernalHmm.GetSingletEmissionScore(transitionOrEmissionInfoVector[i].emissionInfo.state,transitionOrEmissionInfoVector[i].emissionInfo.nuc);
		}
		else {
			globalVars[i]=infernalHmm.GetNthChildTsc(transitionOrEmissionInfoVector[i].edgeInfo.fromState,transitionOrEmissionInfoVector[i].edgeInfo.childNum);
		}
	}
}
void SetGlobalVarsIntoInfernalHmm (InfernalHmm& infernalHmm,const std::vector<double>& globalVars,const TransitionOrEmissionInfoVector& transitionOrEmissionInfoVector)
{
	size_t numGlobalVars=transitionOrEmissionInfoVector.size();
	assert(numGlobalVars==globalVars.size());

	for (size_t i=0; i<numGlobalVars; i++) {
		assert(transitionOrEmissionInfoVector[i].isUsed);
		if (transitionOrEmissionInfoVector[i].isEmission) {
			infernalHmm.SetSingletEmissionLogScore(transitionOrEmissionInfoVector[i].emissionInfo.state,transitionOrEmissionInfoVector[i].emissionInfo.nuc,(float)(globalVars[i]));
		}
		else {
			infernalHmm.SetTransitionLogScore(transitionOrEmissionInfoVector[i].edgeInfo.fromState,transitionOrEmissionInfoVector[i].edgeInfo.childNum,(float)(globalVars[i]));
		}
	}
}


void AddCrossProductOfInequalityLists(InequalityList& inequalityList,InequalityListList& inequalityListList)
{
	// if one of the list is empty, we should just preserve the list (I guess it's not quite a proper cross product).  In other words, remove any empty inequality lists from the list of lists.
	std::list<InequalityListList::iterator> killList;
	for (InequalityListList::iterator i=inequalityListList.begin(); i!=inequalityListList.end(); i++) {
		if (i->empty()) {
			killList.push_back(i);
		}
	}
	for (std::list<InequalityListList::iterator>::iterator k=killList.begin(); k!=killList.end(); k++) {
		inequalityListList.erase(*k);
	}

	if (inequalityListList.size()==0) {
		return;
	}
	if (inequalityListList.size()==1) {
		inequalityList.insert(inequalityList.end(),inequalityListList.front().begin(),inequalityListList.front().end());
		return;
	}
	if (inequalityListList.size()==2) {


		InequalityListList::const_iterator lli=inequalityListList.begin();
		const InequalityList& l1=*lli;
		lli++;
		const InequalityList& l2=*lli;

		for (InequalityList::const_iterator i1=l1.begin(); i1!=l1.end(); i1++) {
			for (InequalityList::const_iterator i2=l2.begin(); i2!=l2.end(); i2++) {

				Inequality newIneq;
				newIneq=*i1;
				newIneq.rhs=i1->rhs + i2->rhs;
				newIneq.lhs.clear();
				newIneq.lhs.insert(newIneq.lhs.end(),i1->lhs.begin(),i1->lhs.end());
				newIneq.lhs.insert(newIneq.lhs.end(),i2->lhs.begin(),i2->lhs.end());

				inequalityList.push_back(newIneq);
			}
		}

		return;
	}

	assert(false); // full generality is not implemented
	Die("something bad @ %s:%d",__FILE__,__LINE__);
}
void AddInequalitiesInTermsOfGlobalVars(InequalityList& inequalityList,const ExtraCm2HmmInfo& extraCm2HmmInfo,const CovarianceModel::Node cmNode)
{
	const std::vector<int>& localToGlobal=extraCm2HmmInfo.inequalitiesAndLocalVariables[cmNode].localToGlobalVariables;
	const InequalityList& nodeInequalityList=extraCm2HmmInfo.inequalitiesAndLocalVariables[cmNode].inequalityList;
	for (InequalityList::const_iterator ineqIter=nodeInequalityList.begin(); ineqIter!=nodeInequalityList.end(); ineqIter++) {
		Inequality newIneq=*ineqIter;
		newIneq.lhs.clear();
		for (std::list<InequalityTerm>::const_iterator termIter=ineqIter->lhs.begin(); termIter!=ineqIter->lhs.end(); termIter++) {
			InequalityTerm newTerm;
			newTerm.variableNum=localToGlobal[termIter->variableNum];
			newIneq.lhs.push_back(newTerm);
		}

		inequalityList.push_back(newIneq);
	}
}
void ConvertInequalityListVarNums (InequalityList& inequalityList,const std::vector<int>& varNumMapping)
{
	for (InequalityList::iterator ineqIter=inequalityList.begin(); ineqIter!=inequalityList.end(); ineqIter++) {
		for (std::list<InequalityTerm>::iterator termIter=ineqIter->lhs.begin(); termIter!=ineqIter->lhs.end(); termIter++) {
			int inputVarNum=termIter->variableNum;
			//("%d   ",inputVarNum);
			int outputVarNum=varNumMapping[inputVarNum];
			assert(outputVarNum!=-1);
			termIter->variableNum=outputVarNum;
		}
	}
}

void VerifyValidHmmBuildType (Cm2Hmm_HmmBuildType hmmType)
{
	switch (hmmType) {
	case HmmBuildType_Original:
	case HmmBuildType_separateMPandMLMR:
	case HmmBuildType_separateMPMLMRD:
		// okay
		break;
	default:
		throw SimpleStringException("The given HMM build type (%d) doesn't appear to be valid",(int)hmmType);
	}
}

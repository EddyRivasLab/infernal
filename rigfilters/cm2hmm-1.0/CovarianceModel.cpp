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

int CovarianceModelBase::nextUniqueId=0;
bool CovarianceModelBase::enableInsertHack=true; // by default, we should be using this for compatibility with Infernal, and also because as I write this, my code is dependent on it
void CovarianceModelBase::DisableInsertHack (void)
{
	enableInsertHack=false;
}
CovarianceModelBase::CovarianceModelBase (void)
{
	cmMemoryOwned=true;
	cm=NULL;
	cmfp=NULL;

	uniqueId=nextUniqueId;
	nextUniqueId++;
}
void CovarianceModelBase::Destruct (void)
{
	if (cm!=NULL && cmMemoryOwned) {
		FreeCM(cm);
		cm=NULL;
	}
	if (cmfp!=NULL) {
		CMFileClose(cmfp);
		cmfp=NULL;
	}
	cmMemoryOwned=true;
}
CovarianceModelBase::~CovarianceModelBase ()
{
	Destruct();
}
int CovarianceModelBase::GetUniqueId (void) const
{
	return uniqueId;
}
void CovarianceModelBase::Save (char *cmFileName)
{
	FILE *out=ThrowingFopen(cmFileName,"wb");
	CMFileWrite(out,cm,FALSE);
	fclose(out);
}
bool CovarianceModelBase::Load (char *cmFileName,bool doLocalAlignment)
{
	// if it looks like it came from tRNAscan-SE, make sure enableInsertHack==true
	if (strlen(cmFileName)>=4) {

		const char *cursor=cmFileName+strlen(cmFileName);
		while (cursor!=cmFileName) {
			if (*cursor=='/' || *cursor=='\\') {
				cursor++;
				break;
			}
			cursor--;
		}

		if ((cursor[0]=='t' || cursor[0]=='T')
		&& (cursor[1]=='r' || cursor[1]=='R')
		&& (cursor[2]=='n' || cursor[2]=='N')
		&& (cursor[3]=='a' || cursor[3]=='A')) {
			if (enableInsertHack) {
				throw SimpleStringException("CM file name begins with 'trna', which suggests it might have come from tRNAscan-SE.  Are you sure you didn't mean to use the --no-insert-hack flag?  To get around this error (__IF__ it's a false alarm), change the file name.\n");
			}
		}
	}

	bool isRsearchCM=false;
	char RSEARCHCM[]=".rsearchcm";
	if (strlen(cmFileName)>strlen(RSEARCHCM)) {
		const char *cursor=cmFileName+strlen(cmFileName)-strlen(RSEARCHCM);
		if (strcmp(cursor,RSEARCHCM)==0) {
			isRsearchCM=true;
		}
	}

	if ((cmfp = CMFileOpen(cmFileName, NULL)) == NULL)
		Die("Failed to open covariance model save file %s\n", cmFileName);

	if (isRsearchCM) {
#ifdef CM2HMM_ONLY
		throw SimpleStringException("Can't handle RSEARCH-format CMs (which use fake probabilities) in this version of the program.");
#else
		fprintf(stderr,"Detected CM from RSEARCH\n");
		printf("Detected CM from RSEARCH\n");
		if (! read_ascii_cm_nonrenormalizing(cmfp,&cm)) {
			Die("Failed to read RSEARCH CM from %s\n",cmFileName);
		}
#endif
	}
	else {
		if (! CMFileRead(cmfp, &cm))
			Die("Failed to read a CM from %s -- file corrupt?\n", cmFileName);
	}
	if (cm == NULL) 
		Die("%s empty?\n", cmFileName);

	if (doLocalAlignment) {

		Die("Local alignments aren't properly handled, since local-ends turned out to be more tricky than I expected.");

		ConfigLocal(cm, 0.5, 0.5);
	}

	CMLogoddsify(cm);
	if (enableInsertHack) {
		CMHackInsertScores(cm);	/* make insert emissions score zero. "TEMPORARY" FIX. */
	}

	// for debugging local stuff -- remember, tho, that endsc can't really be debugged easily
	/*
	int state;
	for (state=0; state!=GetNumStates(); state++) {
		cm->endsc[state]=(float)IMPOSSIBLE;
	}
	for (state=0; state!=363; state++) {
		//cm->beginsc[state]=(float)IMPOSSIBLE;
	}
	for (state=364; state!=GetNumStates(); state++) {
		//cm->beginsc[state]=(float)IMPOSSIBLE;
	}
	//cm->beginsc[3]=0;
	//cm->beginsc[363]=(float)-7.72;

	for (state=363; state!=GetNumStates(); state++) {
		//printf("%d,%f,%f\n",state,cm->beginsc[state],cm->endsc[state]);
	}
	//cm->endsc[363]=(float)IMPOSSIBLE;
	*/

	/*
	FILE *f=fopen("local.cm","wt");
	CMFileWrite(f,cm,0);
	fclose(f);
	*/

	return true;
}

namespace tRNAscanSE { // stuff related to reverse engineering of tRNAscan-SE by Todd Lowe & Sean Eddy

	// this is from tRNAscan-SE-1.23/structs.h
	enum NodeType {
		BIFURC_NODE=0,
		MATP_NODE=1,
		MATL_NODE=2,
		MATR_NODE=3,
		BEGINL_NODE=4,
		BEGINR_NODE=5,
		ROOT_NODE=6,
		END_NODE=BIFURC_NODE
	};
	enum StateType {
		DEL_ST=0,
		MATP_ST=1,
		MATL_ST=2,
		MATR_ST=3,
		INSL_ST=4,
		INSR_ST=5,
		STATETYPES=6, /* MATP nodes contain 6 states */
		BEGIN_ST=DEL_ST,
		BIFURC_ST=DEL_ST,
		END_ST=DEL_ST,
	};

	struct NodeTypeInfo {
		NodeType nodeType;
		int infernalNodeType;

		int numSplitSetStates;
		int numInsertStates;
		StateType *states;
		bool firstStateIsStart;
	};
	StateType matpStates[]={MATP_ST,MATL_ST,MATR_ST,DEL_ST,INSL_ST,INSR_ST};
	StateType matlStates[]={MATL_ST,DEL_ST,INSL_ST};
	StateType matrStates[]={MATR_ST,DEL_ST,INSR_ST};
	StateType rootStates[]={DEL_ST,INSL_ST,INSR_ST};
	StateType beglStates[]={DEL_ST};
	StateType begrStates[]={DEL_ST,INSL_ST};
	StateType bifStates[]={BIFURC_ST};
	NodeTypeInfo nodeTypeInfo[]={
		{BIFURC_NODE,BIF_nd,1,0,bifStates,false}, // BIF handled specially, though
		{MATP_NODE,MATP_nd,4,2,matpStates,false},
		{MATL_NODE,MATL_nd,2,1,matlStates,false},
		{MATR_NODE,MATR_nd,2,1,matrStates,false},
		{BEGINL_NODE,BEGL_nd,1,0,beglStates,true},
		{BEGINR_NODE,BEGR_nd,1,1,begrStates,true},
		{ROOT_NODE,ROOT_nd,1,2,rootStates,true}
	};
	NodeTypeInfo GetNodeTypeInfo (NodeType nodeType) {
		for (unsigned int i=0; i<sizeof(nodeTypeInfo)/sizeof(NodeTypeInfo); i++) {
			if (nodeTypeInfo[i].nodeType==nodeType) {
				return nodeTypeInfo[i];
			}
		}
		throw SimpleStringException("GetNodeTypeInfo failed for node type %d",(int)nodeType);
	}

	struct StateTypeInfo {
		StateType stateType;
		int infernalStateType;

		int numEmits;
	};
	StateTypeInfo stateTypeInfo[]={
		{DEL_ST,D_st,0},
		{MATP_ST,MP_st,2},
		{MATL_ST,ML_st,1},
		{MATR_ST,MR_st,1},
		{INSL_ST,IL_st,1},
		{INSR_ST,IR_st,1},
		{BIFURC_ST,B_st,0},
		{BEGIN_ST,S_st,0},
		{END_ST,E_st,0},
		{END_ST,EL_st,0}
	};
	StateTypeInfo GetStateTypeInfo (StateType stateType) {
		for (unsigned int i=0; i<sizeof(stateTypeInfo)/sizeof(StateTypeInfo); i++) {
			if (stateTypeInfo[i].stateType==stateType) {
				return stateTypeInfo[i];
			}
		}
		throw SimpleStringException("GetStateTypeInfo failed for state type %d",(int)stateType);
	}
	StateType InfernalToCoveStateType (int infernalStateType) {
		for (unsigned int i=0; i<sizeof(stateTypeInfo)/sizeof(StateTypeInfo); i++) {
			if (stateTypeInfo[i].infernalStateType==infernalStateType) {
				return stateTypeInfo[i].stateType;
			}
		}
		throw SimpleStringException("InfernalToCoveStateType failed for state type %d",infernalStateType);
	}

	typedef float qTransitionMatrix[STATETYPES][STATETYPES];
	struct TransitionMatrix {qTransitionMatrix array;};
	typedef float qSingletEmit[4];
	struct SingletEmit {qSingletEmit array;} ;
	typedef float qPairEmit[4][4];
	struct PairEmit {qPairEmit array;} ;
};
void CopySingletEmitFromCove(CM_t *cm,int state,const tRNAscanSE::SingletEmit& coveEmit)
{
	for (int nuc=0; nuc<4; nuc++) {
		cm->e[state][nuc]=coveEmit.array[nuc];
		cm->esc[state][nuc]=(float)(sreLOG2(coveEmit.array[nuc]));
	}
}
void CovarianceModelBase::LoadCove2 (const char *cmFileName)
{
	Destruct();

	bool useLocalEndStates=false; // I wonder what EL_st states are for, looks like it's always E_st

	FILE *cmFile=ThrowingFopen(cmFileName,"rt");
	char version[256]="";
	fscanf(cmFile,"### cove %s\n",version);
	if (strcmp(version,"V2")!=0) {
		throw SimpleStringException("Cove-format file didn't begin with \"### cove V2\" --> it's not in Cove V2 format.");
	}
	int numNodes;
	fscanf(cmFile, "%d \tnodes\n", &numNodes);

	// can allocate Infernal structure now, possibly using more memory than we need (because we're using an upper bound for required # of states), but I don't think that matters
	int upperBoundOnNumStates=numNodes*6;
	int upperBoundOnNumNodes=numNodes*2;
	cm=CreateCM(upperBoundOnNumNodes,upperBoundOnNumStates);
	assert(cm->nodes==upperBoundOnNumNodes); // just verifying my expectation of the CreateCM function
	assert(cm->flags==0); // shouldn't have any set, right?

	cm->name=(char *)MallocOrDie(strlen(cmFileName)+1);
	strcpy(cm->name,cmFileName);

	cm->null[0]=cm->null[1]=cm->null[2]=cm->null[3]=0.25; // default NULL model

	vector<int> coveLeftChild(numNodes),coveRightChild(numNodes);
	vector<tRNAscanSE::NodeType> coveNodeType(numNodes);
	vector<tRNAscanSE::SingletEmit> coveEmitIL(numNodes),coveEmitIR(numNodes),coveEmitML(numNodes),coveEmitMR(numNodes);
	vector<tRNAscanSE::TransitionMatrix> coveTransitionMatrix(numNodes);
	vector<tRNAscanSE::PairEmit> coveEmitMP(numNodes);
	vector<int> infernalNodeToCoveNode(upperBoundOnNumNodes); // upper bound on # of nodes, if every cove node were an END node

	// first, read in the file
	for (int currNode=0; currNode<numNodes; currNode++) {
		int currNodeFromFile;
		fscanf(cmFile, "### node %d",&currNodeFromFile);
		if (currNode!=currNodeFromFile) {
			throw SimpleStringException("Cove V2 file's nodes are not in sequential order, but this dinky code assumes they are.\n");
		}

		int thisCoveNodeType;
		fscanf(cmFile, " type %d\n", &thisCoveNodeType);
		coveNodeType[currNode]=(tRNAscanSE::NodeType)thisCoveNodeType;

		int leftCoveChildNode,rightCoveChildNode;
		fscanf(cmFile, "%d  %d\n", &leftCoveChildNode,&rightCoveChildNode);
		if (leftCoveChildNode!=-1 && leftCoveChildNode!=currNode+1) {
			throw SimpleStringException("Cove V2 file's left child (of node #%d) is not the next node or -1, which I was assuming it would be",currNode);
		}
		coveLeftChild[currNode]=leftCoveChildNode;
		coveRightChild[currNode]=rightCoveChildNode;

		for (int from=0; from<tRNAscanSE::STATETYPES; from++) {
			for (int to=0; to<tRNAscanSE::STATETYPES; to++) {
				fscanf(cmFile,"%f ",&(coveTransitionMatrix[currNode].array[from][to]));
			}
			fscanf(cmFile, "\n");
		}

		for (int i = 0; i < 4; i++) {
			fscanf(cmFile, "%f ", &(coveEmitIL[currNode].array[i]));
		}
		fscanf(cmFile, "# INSL\n");
		for (int i = 0; i < 4; i++) {
			fscanf(cmFile, "%f ", &(coveEmitIR[currNode].array[i]));
		}
		fscanf(cmFile, "# INSR\n");
		for (int i = 0; i < 4; i++) {
			for (int j=0; j<4; j++) {
				fscanf(cmFile, "%f ", &(coveEmitMP[currNode].array[i][j]));
			}
			fscanf(cmFile, "# MATP\n");
		}
		for (int i = 0; i < 4; i++) {
			fscanf(cmFile, "%f ", &(coveEmitML[currNode].array[i]));
		}
		fscanf(cmFile, "# MATL\n");
		for (int i = 0; i < 4; i++) {
			fscanf(cmFile, "%f ", &(coveEmitMR[currNode].array[i]));
		}
		fscanf(cmFile, "# MATR\n");
	}

	fclose(cmFile);

	// now, linearize the states, as in Infernal
	int nextInfernalState=0;
	int nextInfernalNode=0;
	for (int currNode=0; currNode<numNodes; currNode++) {

		infernalNodeToCoveNode[nextInfernalNode]=currNode;

		cm->nodemap[nextInfernalNode]=nextInfernalState;
		tRNAscanSE::NodeTypeInfo nodeTypeInfo=tRNAscanSE::GetNodeTypeInfo (coveNodeType[currNode]);

		cm->ndtype[nextInfernalNode]=nodeTypeInfo.infernalNodeType;
		if (coveNodeType[currNode]==tRNAscanSE::BIFURC_NODE) {
			if (coveLeftChild[currNode]==-1 && coveRightChild[currNode]==-1) {
				// it's really an end node
				cm->ndtype[nextInfernalNode]=END_nd;
				cm->sttype[nextInfernalState]=EL_st;
				cm->ndidx[nextInfernalState]=nextInfernalNode;
				cm->cfirst[nextInfernalState]=-1;
				cm->cnum[nextInfernalState]=0;
				if (currNode+1==numNodes || !useLocalEndStates) {
					// final end state
					cm->sttype[nextInfernalState]=E_st;
				}
			}
			else {
				if (coveLeftChild[currNode]==-1 || coveRightChild[currNode]==-1) {
					throw SimpleStringException("Cove V2 had Bif/End node (#%d) that had exactly one child.  I agree that's senseless.",currNode);
				}

				// bif node
				cm->sttype[nextInfernalState]=B_st;
				cm->ndidx[nextInfernalState]=nextInfernalNode;
				// set up the children later
			}
			nextInfernalState++;
		}
		else {

			int numSplitSetStates=nodeTypeInfo.numSplitSetStates;
			int numInsertStates=nodeTypeInfo.numInsertStates;

			int nextSplitSetStates;
			if (coveLeftChild[currNode]==-1 && coveRightChild[currNode]==-1) {
				// next node is a phantom END node
				nextSplitSetStates=1;
			}
			else {
				assert(currNode+1<numNodes); // else this should have been an END node
				tRNAscanSE::NodeTypeInfo nextNodeTypeInfo=tRNAscanSE::GetNodeTypeInfo (coveNodeType[currNode+1]);
				nextSplitSetStates=nextNodeTypeInfo.numSplitSetStates;
			}

			int firstChildState=nextInfernalState+numSplitSetStates;
			for (int currState=0; currState<numSplitSetStates+numInsertStates; currState++) {
				tRNAscanSE::StateType coveStateType=nodeTypeInfo.states[currState];
				tRNAscanSE::StateTypeInfo stateTypeInfo=GetStateTypeInfo(coveStateType);
				cm->sttype[nextInfernalState]=stateTypeInfo.infernalStateType;
				if (nodeTypeInfo.firstStateIsStart && currState==0) {
					cm->sttype[nextInfernalState]=S_st;
				}
				cm->ndidx[nextInfernalState]=nextInfernalNode;
				cm->cfirst[nextInfernalState]=firstChildState;
				cm->cnum[nextInfernalState]=nextSplitSetStates+nodeTypeInfo.numInsertStates;
				if (cm->sttype[nextInfernalState]==IR_st && numInsertStates==2) {
					// special case -- IR doesn't go to IL
					cm->cfirst[nextInfernalState]++;
					cm->cnum[nextInfernalState]--;
				}
				
				// emits
				switch (stateTypeInfo.numEmits) {
					case 0:
						// nothing to do for emits
						break;
					case 1:
						if (coveStateType==tRNAscanSE::INSL_ST) { // MSVC++ screws up tabbing if I use switch/case
							CopySingletEmitFromCove(cm,nextInfernalState,coveEmitIL[currNode]);
						}
						if (coveStateType==tRNAscanSE::INSR_ST) { // MSVC++ screws up tabbing if I use switch/case
							CopySingletEmitFromCove(cm,nextInfernalState,coveEmitIR[currNode]);
						}
						if (coveStateType==tRNAscanSE::MATL_ST) { // MSVC++ screws up tabbing if I use switch/case
							CopySingletEmitFromCove(cm,nextInfernalState,coveEmitML[currNode]);
						}
						if (coveStateType==tRNAscanSE::MATR_ST) { // MSVC++ screws up tabbing if I use switch/case
							CopySingletEmitFromCove(cm,nextInfernalState,coveEmitMR[currNode]);
						}
						break;
					case 2:
						{
							for (int leftNuc=0; leftNuc<4; leftNuc++) {
								for (int rightNuc=0; rightNuc<4; rightNuc++) {
									float emitProb=coveEmitMP[currNode].array[leftNuc][rightNuc];
									cm->e[nextInfernalState][GetPairIndex(leftNuc,rightNuc)]=emitProb;
									cm->esc[nextInfernalState][GetPairIndex(leftNuc,rightNuc)]=(float)(sreLOG2(emitProb));
								}
							}
						}
						break;
					default:
						assert(false);
						break;
				}

				nextInfernalState++;
			}

			// (annoying) special case: a node of regular type (like MATL) can also be an END node.  In this case, we have to put in the INFERNAL END node afterwards.
			if (coveLeftChild[currNode]==-1 && coveRightChild[currNode]==-1) {
				nextInfernalNode++;
				infernalNodeToCoveNode[nextInfernalNode]=currNode;
				cm->ndtype[nextInfernalNode]=END_nd;
				cm->nodemap[nextInfernalNode]=nextInfernalState;
				cm->sttype[nextInfernalState]=EL_st;
				cm->ndidx[nextInfernalState]=nextInfernalNode;
				cm->cfirst[nextInfernalState]=-1;
				cm->cnum[nextInfernalState]=0;
				if (currNode+1==numNodes || !useLocalEndStates) {
					// final end state
					cm->sttype[nextInfernalState]=E_st;
				}
				nextInfernalState++;
			}
		}

		nextInfernalNode++;
	}
	cm->M=nextInfernalState;
	cm->nodes=nextInfernalNode;

	// do transitions, now that every state has its type
	for (int state=0; state<cm->M; state++) {
		int currInfernalNode=cm->ndidx[state];
		int currNode=infernalNodeToCoveNode[currInfernalNode];

		if (cm->sttype[state]!=B_st) {
			tRNAscanSE::StateType coveFromStateType=tRNAscanSE::InfernalToCoveStateType(cm->sttype[state]);
			for (int i=0; i<cm->cnum[state]; i++) {
				int childState=cm->cfirst[state]+i;
				tRNAscanSE::StateType coveToStateType=tRNAscanSE::InfernalToCoveStateType(cm->sttype[childState]);
				float p=coveTransitionMatrix[currNode].array[coveFromStateType][coveToStateType];
				cm->t[state][i]=p;
				cm->tsc[state][i]=(float)(sreLOG2(p));
			}
		}
	}

	// now set up right children of BIF nodes
	// first reverse-map nodes
	vector<int> coveNodeToInfernalNode(numNodes);
	for (int infernalNode=0; infernalNode<cm->nodes; infernalNode++) {
		coveNodeToInfernalNode[infernalNodeToCoveNode[infernalNode]]=infernalNode;
	}
	// now do it
	for (int currNode=0; currNode<numNodes; currNode++) {
		if (coveNodeType[currNode]==tRNAscanSE::BIFURC_NODE) {
			if (coveLeftChild[currNode]!=-1 && coveRightChild[currNode]!=-1) {
				int bifState=cm->nodemap[coveNodeToInfernalNode[currNode]];
				int leftChildState=cm->nodemap[coveNodeToInfernalNode[coveLeftChild[currNode]]];
				int rightChildState=cm->nodemap[coveNodeToInfernalNode[coveRightChild[currNode]]];
				cm->cfirst[bifState]=leftChildState;
				cm->cnum[bifState]=rightChildState;
			}
		}
	}

	// now set up parent pointers & full stid
	// first clear parents
	for (int state=0; state<cm->M; state++) {
		cm->plast[state]=-1;
		cm->pnum[state]=0;
	}
	// now set everything
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		cm->stid[StateToInt(state)]=DeriveUniqueStateCode(GetNodeType(GetNode(state)), GetStateType(state));
		if (IsBifurcation(state)) {
			State childState=GetLeftBifurcationChild(state);
			cm->plast[StateToInt(childState)]=StateToInt(state);
			cm->pnum[StateToInt(childState)]++;
			childState=GetRightBifurcationChild(state);
			cm->plast[StateToInt(childState)]=StateToInt(state);
			cm->pnum[StateToInt(childState)]++;
		}
		else {
			for (int childNum=0; childNum<GetNumChildren(state); childNum++) {
				State childState=GetNthChildState(state,childNum);
				cm->plast[StateToInt(childState)]=StateToInt(state);
				cm->pnum[StateToInt(childState)]++;
			}
		}
	}
}
void CovarianceModelBase::CopyFrom (const CovarianceModelBase& t)
{
	Destruct();

	cm=CreateCM(t.GetNumNodes(),t.GetNumStates());

	cm->name=(char *)MallocOrDie(strlen(t.GetName())+1);
	strcpy(cm->name,t.GetName());

	COPY_ARRAY(sttype,t.GetNumStates()+1);
	COPY_STATE_ARRAY(ndidx);
	COPY_ARRAY(stid,t.GetNumStates()+1);
	COPY_STATE_ARRAY(cfirst);
	COPY_STATE_ARRAY(cnum);
	COPY_STATE_ARRAY(plast);
	COPY_STATE_ARRAY(pnum);

	COPY_ARRAY(nodemap,t.GetNumNodes());
	COPY_ARRAY(ndtype,t.GetNumNodes());

	COPY_ARRAY_2D(t,t.GetNumStates(),MAXCONNECT);
	COPY_ARRAY_2D(e,t.GetNumStates(),Alphabet_size*Alphabet_size);
	COPY_STATE_ARRAY(begin);
	COPY_STATE_ARRAY(end);
	COPY_ARRAY_2D(tsc,t.GetNumStates(),MAXCONNECT);
	COPY_ARRAY_2D(esc,t.GetNumStates(),Alphabet_size*Alphabet_size);
	COPY_STATE_ARRAY(beginsc);
	COPY_STATE_ARRAY(endsc);

	cm->flags  = t.cm->flags;
}
void CovarianceModelBase::MirrorFrom (CM_t *t)
{
	Destruct();

	cmMemoryOwned=false;
	cm=t;
	cmfp=NULL;
}

void CovarianceModelBase::Realloc2d(float **(&array),int oldNumStates,int newNumStates,int sizeDim2)
{
	float **old=array;
	array=FMX2Alloc(newNumStates, sizeDim2);
	for (int i=0; i<oldNumStates; i++) {
		for (int j=0; j<sizeDim2; j++) {
			array[i][j]=old[i][j];
		}
	}
	FMX2Free(old);
}
// WARNING: only for HMMs
void CovarianceModelBase::AddStates(int numNewStates)
{
	int newNumStates=numNewStates+GetNumStates();

	cm->sttype=(char *)ReallocOrDie(cm->sttype,sizeof(cm->sttype[0])*(newNumStates+1));
	cm->stid=(char *)ReallocOrDie(cm->stid,sizeof(cm->stid[0])*(newNumStates+1));
	cm->ndidx=(int *)ReallocOrDie(cm->ndidx,sizeof(cm->ndidx[0])*newNumStates);
	cm->cfirst=(int *)ReallocOrDie(cm->cfirst,sizeof(cm->cfirst[0])*newNumStates);
	cm->cnum=(int *)ReallocOrDie(cm->cnum,sizeof(cm->cnum[0])*newNumStates);
	cm->plast=(int *)ReallocOrDie(cm->plast,sizeof(cm->plast[0])*newNumStates);
	cm->pnum=(int *)ReallocOrDie(cm->pnum,sizeof(cm->pnum[0])*newNumStates);

	cm->begin=(float *)ReallocOrDie(cm->begin,sizeof(cm->begin[0])*newNumStates);
	cm->end=(float *)ReallocOrDie(cm->end,sizeof(cm->end[0])*newNumStates);
	cm->beginsc=(float *)ReallocOrDie(cm->beginsc,sizeof(cm->beginsc[0])*newNumStates);
	cm->endsc=(float *)ReallocOrDie(cm->endsc,sizeof(cm->endsc[0])*newNumStates);

	Realloc2d(cm->tsc,GetNumStates(),newNumStates,MAXCONNECT);
	Realloc2d(cm->esc,GetNumStates(),newNumStates,Alphabet_size*Alphabet_size);

	Realloc2d(cm->t,GetNumStates(),newNumStates,MAXCONNECT);
	Realloc2d(cm->e,GetNumStates(),newNumStates,MAXCONNECT);

	cm->M=newNumStates;
}
// WARNING: only for HMMs
void CovarianceModelBase::MoveStatesHigher (State st_firstState,State st_destOfFirstState)
{
	int firstState=StateToInt(st_firstState);
	int destOfFirstState=StateToInt(st_destOfFirstState);

	assert(destOfFirstState>=firstState);
	int increase=destOfFirstState-firstState;
	int newNumStates=GetNumStates() + increase;

	AddStates(increase);

	// do the move, careful about the overlap
	for (int state=newNumStates-1; state>=destOfFirstState; state--) {
		cm->sttype[state]=cm->sttype[state-increase];
		cm->cfirst[state]=cm->cfirst[state-increase];
		cm->cnum[state]=cm->cnum[state-increase];
		for (int i=0; i<MAXCONNECT; i++) {
			cm->tsc[state][i]=cm->tsc[state-increase][i];
		}
		for (int i=0; i<Alphabet_size*Alphabet_size; i++) {
			cm->esc[state][i]=cm->esc[state-increase][i];
		}
	}

	// now adjust 'cfirst' members of the moved thingies
	for (int state=destOfFirstState; state<newNumStates; state++) {
		if (cm->cfirst[state]!=-1) {
			assert(cm->cfirst[state]>=firstState);
			cm->cfirst[state] += increase;
		}
	}
}
void CovarianceModelBase::CopyStatesVerbatimFrom(const CovarianceModelBase& t,State st_firstState,State st_lastState)
{
	int firstState=StateToInt(st_firstState);
	int lastState=StateToInt(st_lastState);
	int state;
	for (state=firstState; state<lastState; state++) {

		cm->sttype[state]=t.cm->sttype[state];
		cm->cfirst[state]=t.cm->cfirst[state];
		cm->cnum[state]=t.cm->cnum[state];
		for (int i=0; i<MAXCONNECT; i++) {
			cm->tsc[state][i]=t.cm->tsc[state][i];
		}
		for (int i=0; i<Alphabet_size*Alphabet_size; i++) {
			cm->esc[state][i]=t.cm->esc[state][i];
		}
	}
}

void CovarianceModelBase::Init (int numStates)
{
	cm=CreateCM(0,numStates);

	char dummyName[]="unnamed";
	cm->name = sre_strdup(dummyName, (int)(strlen(dummyName)));
}
void CovarianceModelBase::DivideOutUniformNullModel(void)
{
	int v, x, y;
	for (v = 0; v < cm->M; v++)
	{
		if (cm->sttype[v] == MP_st)
			for (x = 0; x < Alphabet_size; x++)
				for (y = 0; y < Alphabet_size; y++) {
					cm->e[v][x*Alphabet_size+y] *= 16.0;
					cm->esc[v][x*Alphabet_size+y] = (float)sreLOG2(cm->e[v][x*Alphabet_size+y]);
				}

		if (cm->sttype[v] == ML_st || cm->sttype[v] == MR_st ||
			cm->sttype[v] == IL_st || cm->sttype[v] == IR_st)
			for (x = 0; x < Alphabet_size; x++) {
				cm->e[v][x] *= 4.0;
				cm->esc[v][x] = (float)sreLOG2(cm->e[v][x]);
			}
	}
}
void CovarianceModelBase::CMRenormalize (void)
{
	::CMRenormalize(cm);

	// copied from ::CMLogoddsify, but in this case we just want to Log
	int v, x, y;
	for (v = 0; v < cm->M; v++)
	{
		if (cm->sttype[v] != B_st && cm->sttype[v] != E_st)
			for (x = 0; x < cm->cnum[v]; x++)
				cm->tsc[v][x] = (float)sreLOG2(cm->t[v][x]);

		if (cm->sttype[v] == MP_st)
			for (x = 0; x < Alphabet_size; x++)
				for (y = 0; y < Alphabet_size; y++)
					cm->esc[v][x*Alphabet_size+y] = (float)sreLOG2(cm->e[v][x*Alphabet_size+y]);

		if (cm->sttype[v] == ML_st || cm->sttype[v] == MR_st ||
			cm->sttype[v] == IL_st || cm->sttype[v] == IR_st)
			for (x = 0; x < Alphabet_size; x++)
				cm->esc[v][x] = (float)sreLOG2(cm->e[v][x]);

		/* These work even if begin/end distributions are inactive 0's,
		* sreLOG2 will set beginsc, endsc to -infinity.
		*/
		cm->beginsc[v] = (float)sreLOG2(cm->begin[v]);
		cm->endsc[v]   = (float)sreLOG2(cm->end[v]);
	}

	if (enableInsertHack) {
		::CMHackInsertScores(cm);
	}
}
void CovarianceModelBase::HackInsertScoresToStrictProbs (void)
{
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		if (IsInsertState(state)) {
			for (int nuc=0; nuc<Alphabet_size; nuc++) {
				SetEmissionProbDirectly(state,nuc,0.25);
				SetSingletEmissionLogScore(state,nuc,-2.0);
			}
		}
	}
}
void CovarianceModelBase::SetTransitionProbDirectly (State state,int child,float prob)
{
	cm->t[StateToInt(state)][child]=prob;
}
float CovarianceModelBase::GetTransitionProbDirectly (State state,int child) const
{
	return cm->t[StateToInt(state)][child];
}
void CovarianceModelBase::SetEmissionProbDirectly (State state,int nuc,float prob)
{
	cm->e[StateToInt(state)][nuc]=prob;
}
float CovarianceModelBase::GetEmissionProbDirectly (State state,int nuc) const
{
	return cm->e[StateToInt(state)][nuc];
}
float CovarianceModelBase::GetPairEmissionProbDirectly (State state,int leftNuc,int rightNuc) const
{
	return cm->e[StateToInt(state)][leftNuc*Alphabet_size + rightNuc];
}
void CovarianceModelBase::SetPairEmissionProbDirectly(State state,int leftNuc,int rightNuc,float prob)
{
	cm->e[StateToInt(state)][leftNuc*Alphabet_size + rightNuc]=prob;
}
void CovarianceModelBase::NormalizeEmissionsToStrictProbabilitiesViaProbabilities (void)
{
	for (State state=GetFirstState(); state!=GetLastState(); state++) {

		switch (GetNumSymbolsEmitted(state)) {
			case 1:
				{
					double totalProb=0;
					for (int nuc=0; nuc<Alphabet_size; nuc++) {
						totalProb += GetEmissionProbDirectly(state,nuc);
					}
					for (int nuc=0; nuc<Alphabet_size; nuc++) {
						float t=GetEmissionProbDirectly(state,nuc);
						float prob=(float)(t/totalProb);
						if (t==0 && totalProb==0) {
							prob=1;
						}
						SetEmissionProbDirectly(state,nuc,prob);
						SetSingletEmissionLogScore(state,nuc,(float)(sreLOG2(prob)));
					}
				}
				break;
			case 2:
				{
					double totalProb=0;
					for (int leftNuc=0; leftNuc<Alphabet_size; leftNuc++) {
						for (int rightNuc=0; rightNuc<Alphabet_size; rightNuc++) {
							totalProb += GetPairEmissionProbDirectly(state,leftNuc,rightNuc);
						}
					}
					for (int leftNuc=0; leftNuc<Alphabet_size; leftNuc++) {
						for (int rightNuc=0; rightNuc<Alphabet_size; rightNuc++) {
							float t=GetPairEmissionProbDirectly(state,leftNuc,rightNuc);
							float prob=(float)(t/totalProb);
							if (t==0 && totalProb==0) {
								prob=1;
							}
							SetPairEmissionProbDirectly(state,leftNuc,rightNuc,prob);
							SetPairEmissionScore(state,leftNuc,rightNuc,(float)(sreLOG2(prob)));
						}
					}
				}
				break;
		}
	}
}
void CovarianceModelBase::NormalizeTransitionsToStrictProbabilitiesViaProbabilities (void)
{
	for (State state=GetFirstState(); state!=GetLastState(); state++) {

		double totalProb=0;
		for (int child=0; child<GetNumChildren(state); child++) {
			totalProb += GetTransitionProbDirectly(state,child);
		}
		for (int child=0; child<GetNumChildren(state); child++) {
			float t=GetTransitionProbDirectly(state,child);
			float prob=(float)(t/totalProb);
			if (t==0 && totalProb==0) {
				prob=1;
			}
			SetTransitionProbDirectly(state,child,prob);
			SetTransitionLogScore(state,child,(float)(sreLOG2(prob)));
		}
	}
}
void CovarianceModelBase::ZeroAllEmitProbs (void)
{
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		assert(GetNumSymbolsEmitted(state)!=2); // only meant for HMMs
		if (IsEmitting(state)) {
			for (int nuc=0; nuc<Alphabet_size; nuc++) {
				SetEmissionProbDirectly(state,nuc,(float)0);
			}
		}
	}
}
void CovarianceModelBase::ZeroAllTransitionProbsExceptSelfLoops (void)
{
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		for (int child=0; child<GetNumChildren(state); child++) {
			State toState=GetNthChildState(state,child);
			if (state!=toState) { // no a self-loop
				SetTransitionProbDirectly(state,child,(float)0);
			}
		}
	}
}
void CovarianceModelBase::ZeroAllTransitionProbs (void)
{
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		for (int child=0; child<GetNumChildren(state); child++) {
			SetTransitionProbDirectly(state,child,(float)0);
		}
	}
}
void CovarianceModelBase::SetFirstChild (State state,State firstChild)
{
	cm->cfirst[StateToInt(state)]=StateToInt(firstChild);
}
void CovarianceModelBase::SetNumChildren (State state,int numChildren)
{
	assert(numChildren>=0 && numChildren<MAXCONNECT);
	cm->cnum[StateToInt(state)]=numChildren;
}
void CovarianceModelBase::SetNoChildren (State state)
{
	cm->cfirst[StateToInt(state)]=-1;
	cm->cnum[StateToInt(state)]=0;
}
void CovarianceModelBase::SetStateType (State state,int stateType)
{
	cm->sttype[StateToInt(state)]=stateType;
}
void CovarianceModelBase::SetTransitionLogScore (State state,int childNum,float score)
{
	assert(childNum>=0 && childNum<GetNumChildren(state));
	cm->tsc[StateToInt(state)][childNum]=score;
}
void CovarianceModelBase::SetSingletEmissionLogScore (State state,int symbol,float score)
{
	assert(symbol>=0 && symbol<Alphabet_size);
	cm->esc[StateToInt(state)][symbol]=score;
}
void CovarianceModelBase::SetBeginsc (State state,float score)
{
	cm->beginsc[StateToInt(state)]=score;
}
float CovarianceModelBase::CYKDivideAndConquer(char *seq,int L,State state,int i,int j,Parsetree_t **retr) const {
	CM_t *nonConstCM=(CM_t *)cm;
	return ::CYKDivideAndConquer(nonConstCM,seq,L,StateToInt(state),i,j,retr);
}
float CovarianceModelBase::CYKInsideScore (char *seq,int L,State state,int i,int j) const {
	CM_t *nonConstCM=(CM_t *)cm;
	return ::CYKInsideScore(nonConstCM,seq,L,StateToInt(state),i,j);
}
CMConsensus_t *CovarianceModelBase::CreateCMConsensus (float x,float y) const
{
	CM_t *nonConstCM=(CM_t *)cm;
	return ::CreateCMConsensus(nonConstCM,x,y);
}
void CovarianceModelBase::CYKScanZasha (char *dsq, int L, int W,int *ret_nhits, int **ret_hitr, int **ret_hiti, int **ret_hitj, float **ret_hitsc,CykscanStats& cykscanStats) const {
#ifndef CM2HMM_ONLY
	CM_t *nonConstCM=(CM_t *)cm;
	::CYKScanZasha(nonConstCM,dsq,L,W,ret_nhits,ret_hitr,ret_hiti,ret_hitj,ret_hitsc,cykscanStats);
#endif
}
void CovarianceModelBase::CYKScan_OptionalBug (char *dsq, int L, int W,int *ret_nhits, int **ret_hitr, int **ret_hiti, int **ret_hitj, float **ret_hitsc,bool fixBug) const
{
#ifndef CM2HMM_ONLY
	CM_t *nonConstCM=(CM_t *)cm;
	::CYKScan_OptionalBug(nonConstCM,dsq,L,W,ret_nhits,ret_hitr,ret_hiti,ret_hitj,ret_hitsc,fixBug);
#endif
}
Fancyali_t *CovarianceModelBase::CreateFancyAli(Parsetree_t *tr, CMConsensus_t *cons, char *rnaSequence) const
{
	CM_t *nonConstCM=(CM_t *)cm;
	return ::CreateFancyAli(tr,nonConstCM,cons,rnaSequence);
}
float CovarianceModelBase::ParsetreeScore(Parsetree_t *tr,char *rnaSequence) const
{
#ifndef CM2HMM_ONLY
	CM_t *nonConstCM=(CM_t *)cm;
	return ::ParsetreeScore(nonConstCM,tr,rnaSequence);
#else
	throw SimpleStringException("Not implemented %s:%d",__FILE__,__LINE__);
#endif
}
void CovarianceModelBase::ParsetreeDump(FILE *out, Parsetree_t *tr, char *rnaSequence) const
{
#ifndef CM2HMM_ONLY
	CM_t *nonConstCM=(CM_t *)cm;
	::ParsetreeDump(out, tr, nonConstCM, rnaSequence);
#endif
}
const char *CovarianceModelBase::GetName (void) const {
	return cm->name;
}
void CovarianceModelBase::CYKScan (char *dsq, int L, int W,int *ret_nhits, int **ret_hitr, int **ret_hiti, int **ret_hitj, float **ret_hitsc) const
{
	CM_t *nonConstCM=(CM_t *)cm;
	::CYKScan(nonConstCM,dsq,L,W,ret_nhits,ret_hitr,ret_hiti,ret_hitj,ret_hitsc);
}

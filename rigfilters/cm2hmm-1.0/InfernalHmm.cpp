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


///////////////////
// InfernalHmm

InfernalHmm::InfernalHmm (void)
{
	loadedFileFormat=GetCurrFileFormatNum();
	hmmBuildType=HmmBuildType_Original; // assume this for now
}
void InfernalHmm::SetHmmBuildType (Cm2Hmm_HmmBuildType hmmBuildType_)
{
	hmmBuildType=hmmBuildType_;
}
void InfernalHmm::SetFromCmFileName (const std::string& _fromCmFileName)
{
	fromCmFileName=_fromCmFileName;
}
void InfernalHmm::AddBuildDescription (const std::string addedBuildDescription)
{
	fullBuildDescription += ";";
	fullBuildDescription += addedBuildDescription;
}
std::string InfernalHmm::GetBuildDescription (void) const
{
	return LineBreaksToTabs(fullBuildDescription);
}
const std::string& InfernalHmm::GetFromCmFileName (void) const
{
	return fromCmFileName;
}
int InfernalHmm::GetFileFormat (void) const
{
	return loadedFileFormat;
}
void InfernalHmm::Init (int numStates)
{
	CovarianceModelBase::Init(numStates);
	endscLinksToLeftVector.resize(numStates);
	otherStateInfoVector.resize(numStates);
	hmm2CmStateVector.resize(numStates);
	cm2HmmStateVector.clear();
	nonSavedInfoVector.clear();
}
void InfernalHmm::BuildNonSavedInfoIfNecessary (void)
{
	if (nonSavedInfoVector.empty()) {

		// must build
		nonSavedInfoVector.resize(GetNumStates());

		State parentState;
		for (parentState=GetFirstState(); parentState!=GetLastState(); parentState++) {
			for (int childNum=0; childNum<GetNumChildren(parentState); childNum++) {

				ParentAndMyChildNum parentAndMyChildNum;
				parentAndMyChildNum.parentState=parentState;
				parentAndMyChildNum.myChildNum=childNum;
				State myState=GetNthChildState(parentState,childNum);

				nonSavedInfoVector[myState].parentAndMyChildNumVector.push_back(parentAndMyChildNum);
			}
		}

		State state;
		for (state=GetFirstState(); state!=GetLastState(); state++) {
			nonSavedInfoVector[state].childNumOfSelfLoop=-1; // assume nothing
			for (int childNum=0; childNum<GetNumChildren(state); childNum++) {
				State childState=GetNthChildState(state,childNum);
				if (state==childState) {
					nonSavedInfoVector[state].childNumOfSelfLoop=childNum;
				}
			}
		}
	}
}
int InfernalHmm::GetChildNumOfSelfLoop (State state) const
{
	assert(!nonSavedInfoVector.empty()); // didn't call BuildNonSavedInfoIfNecessary
	return nonSavedInfoVector[state].childNumOfSelfLoop;
}
void InfernalHmm::SetHmm2CmState (const Hmm2CmStateVector& _hmm2CmStateVector)
{
	hmm2CmStateVector=_hmm2CmStateVector;
}
void InfernalHmm::SetCm2HmmState (const Cm2HmmStateVector& _cm2HmmStateVector)
{
	cm2HmmStateVector=_cm2HmmStateVector;
}
void InfernalHmm::CopyFrom (const InfernalHmm& t)
{
	CovarianceModelBase::CopyFrom(t);
	endscLinksToLeftVector=t.endscLinksToLeftVector;
	otherStateInfoVector=t.otherStateInfoVector;
	hmm2CmStateVector=t.hmm2CmStateVector;
	cm2HmmStateVector=t.cm2HmmStateVector;
	nonSavedInfoVector=t.nonSavedInfoVector;
	loadedFileFormat=t.loadedFileFormat;
	fullBuildDescription=t.fullBuildDescription;
	fromCmFileName=t.fromCmFileName;
	hmmBuildType=t.hmmBuildType;
}
void InfernalHmm::MirrorFromWithHackedExtraInfo(CM_t *cm)
{
	MirrorFrom(cm);
	endscLinksToLeftVector.resize(GetNumStates());
	otherStateInfoVector.resize(GetNumStates());
	hmm2CmStateVector.resize(GetNumStates());
	cm2HmmStateVector.clear();
}
void InfernalHmm::ClobberIR (void)
{
	State state=GetFirstState();
	state++;
	state++;
	if (GetStateType(state)!=IR_st) {
		Die("InfernalHmm::ClobberIR: expected IR to be 3rd state.");
	}

	// change IR to a D state, with the same transitions.  I think this is a reasonable hack, and whatever we do the IR state will have low probability and there's only 1 IR state, so it won't make that much of a difference.
	SetStateType(state,D_st);

	// now verify there's nothing else right-emitting.
	for (state=GetFirstState(); state!=GetLastState(); state++) {
		if (EmitsRight(state)) {
			Die("InfernalHmm::ClobberIR: found unexpected right-emitting state.");
		}
	}
}
void InfernalHmm::CopyFromWithEscHack(const InfernalHmm& t)
{
	CopyFrom(t);
	
	for (State state=GetFirstState(); state!=GetLastState(); state++) {
		if (IsEmitting(state)) {
			for (int nuc=0; nuc<Alphabet_size; nuc++) {
				float esc=GetSingletEmissionScore(state,nuc);
				esc += (float)2.0;
				SetSingletEmissionLogScore (state,nuc,esc);
			}
		}
	}
}
void InfernalHmm::SetLeftState (State state,bool isLeftState)
{
	otherStateInfoVector[state].isLeftState=isLeftState;
}
void InfernalHmm::SetRightState (State state,bool isRightState)
{
	otherStateInfoVector[state].isRightState=isRightState;
}
void InfernalHmm::SetLeftwardBeginsc (State state,float sc)
{
	otherStateInfoVector[state].leftwardBeginsc=sc;
}
void InfernalHmm::SetRightwardBeginsc (State state,float sc)
{
	otherStateInfoVector[state].rightwardBeginsc=sc;
}
void InfernalHmm::DumpInfernalHmm (FILE *file,const CovarianceModel& cm) const
{
	if (!disableHmmBuildInfoDump) {
		fprintf(file,"----hmm-from-cmFileName: %s\n",GetFromCmFileName().c_str());
		fprintf(file,"----hmm-build-info: %s\n",GetBuildDescription().c_str());
		fprintf(file,"----fileFormatNum: %d\n",loadedFileFormat);
	}

	int numStates=GetNumStates();
	fprintf(file,"# states=%d\n",numStates);

	State state;
	for (state=GetFirstState(); state!=GetLastState(); state++) {
		fprintf(file,"\tState #%d: type=%s.  (left=%s,right=%s)",StateToInt(state),GetStateTypeName(state),IsLeftState(state)?"T":"f",IsRightState(state)?"T":"f");
		if (DoLocal()) {
			fprintf(file,"   beginsc left=%f, right=%f, links-to-left: ",GetLeftwardBeginsc(state),GetRightwardBeginsc(state));
			for (int linkNum=0; linkNum<GetNumEndscLinksToLeft(state); linkNum++) {
				fprintf(file,"(%d,%f) ",StateToInt(GetEndscLinkToLeft_State(state,linkNum)),GetEndscLinkToLeft_Endsc(state,linkNum));
			}
		}
		fprintf(file,"\n");
		if (HasChildren(state)) {
			int numChildren=GetNumChildren(state);
			fprintf(file,"\t\tChildren = ",StateToInt(GetNthChildState(state,0)),StateToInt(GetNthChildState(state,numChildren-1)));
			for (int childNum=0; childNum<numChildren; childNum++) {
				fprintf(file,"state %d (score=%f)  ",StateToInt(GetNthChildState(state,childNum)),GetNthChildTsc(state,childNum));
			}
			fprintf(file,"\n");
		}
		else {
			fprintf(file,"\t\tChildren = none.\n");
		}
		if (IsEmitting(state)) {
			fprintf(file,"\t\tEmits: ");
			for (int nuc=0; nuc<Alphabet_size; nuc++) {
				fprintf(file,"%c (score=%f)  ",nucs[nuc],GetSingletEmissionScore(state,nuc));
			}
			fprintf(file,"\n");
		}
	}

	fprintf(file,"\nCM to HMM mappings\n");
	if (cm2HmmStateVector.empty()) {
		fprintf(file,"\tMAPPINGS HAVE NOT BE CALCULATED & RECORDED YET.\n");
	}
	else {
		fprintf(file,"\n\tCM state --> HMM left state , HMM right state\n");
		for (CovarianceModel::State cmState=cm.GetFirstState(); cmState!=cm.GetLastState(); cmState++) {
			State hmmLeftState=cm2HmmStateVector[cmState].hmmLeftState;
			State hmmRightState=cm2HmmStateVector[cmState].hmmRightState;
			char leftStr[32],rightStr[32];
			if (hmmLeftState==GetInvalidState()) {
				sprintf(leftStr,"---");
			}
			else {
				sprintf(leftStr,"%d",StateToInt(hmmLeftState));
			}
			if (hmmRightState==GetInvalidState()) {
				sprintf(rightStr,"---");
			}
			else {
				sprintf(rightStr,"%d",StateToInt(hmmRightState));
			}
			fprintf(file,"\t%d --> %s , %s\n",CovarianceModel::StateToInt(cmState),leftStr,rightStr);
		}
		fprintf(file,"\n\tCM node --> {HMM state set}\n");
		for (CovarianceModel::Node cmNode=cm.GetFirstNode(); cmNode!=cm.GetLastNode(); cmNode++) {

			fprintf(file,"\t%d --> {",CovarianceModel::NodeToInt(cmNode));

			StateFromNodeList hmmStateList;
			GetHmmStatesByCmNode (hmmStateList,cm,cmNode);
			for (StateFromNodeList::iterator i=hmmStateList.begin(); i!=hmmStateList.end(); i++) {
				fprintf(file," %d",StateToInt(i->state));
			}
			fprintf(file,"}\n");
		}
	}
}
void InfernalHmm::AllocHmmData (void)
{
	endscLinksToLeftVector.resize(GetNumStates());
	otherStateInfoVector.resize(GetNumStates());
}
int InfernalHmm::GetCurrFileFormatNum (void)
{
	return 4;
}
bool InfernalHmm::LoadInBinary (FILE *file)
{
	char firstFourChars[5];
	fread(firstFourChars,4,1,file);
	firstFourChars[4]=0;
	if (strcmp(firstFourChars,"Prof")==0) {

		// new text format

		fscanf(file,"ileHMMFilter\n");

		if (fscanf(file,"FileFormatNum:\t%d\n",&loadedFileFormat)!=1) {
			throw SimpleStringException("Input profile HMM missing 'FileFormatNum'");
		}

		int numStates;
		if (fscanf(file,"states:\t%d\n",&numStates)!=1) {
			throw SimpleStringException("Input profile HMM missing 'states'");
		}
		Init(numStates);

		char cmFileNameTemp[4096];
		if (fscanf(file,"cmFileName:\t%s\n",cmFileNameTemp)!=1) {
			throw SimpleStringException("Input profile HMM missing 'cmFileName'");
		}
		fromCmFileName=cmFileNameTemp;

		fullBuildDescription.clear();
		fscanf(file,"hmmBuildCommands: ");
		while (true) {
			int ch=fgetc(file);
			if (ch=='\r' || ch=='\n' || ch==EOF) {
				break;
			}
			if (ch=='\t') {
				ch='\n'; // I'm finicky about restoring the '\n's because by making things exactly the same, it helps me to show that this text format works like the binary format by loading a binary-format file, saving in text format, loading that, saving in binary, and verifying that the file is the same as the original.  Other than that, it doesn't matter.
			}
			fullBuildDescription += ch;
		}

		char hmmTypeText[256];
		if (fscanf(file,"hmmType:\t%s\n",hmmTypeText)!=1) {
			throw SimpleStringException("Input profile HMM missing 'hmmType'");
		}
		hmmBuildType=GetHmmBuildTypeByText(hmmTypeText);

		int doLocalInt;
		if (fscanf(file,"isLocal:\t%d\n",&doLocalInt)!=1) {
			throw SimpleStringException("Input profile HMM missing 'isLocal'");
		}
		SetDoLocal(doLocalInt!=0);
		if (DoLocal()) {
			throw SimpleStringException("Hey, wait!  Input profile HMM specifies that it uses local alignments, but I haven't implemented that yet.  Are you sure about that?");
		}

		for (State state=GetFirstState(); state!=GetLastState(); state++) {
			int i_state=StateToInt(state);

			char sttypeText[256];
			int stid;
			if (fscanf(file,"%d\t%s\t%d\t%d\t",&(stid),sttypeText,&(cm->cfirst[i_state]),&(cm->cnum[i_state]))!=4) {
				throw SimpleStringException("Input profile HMM state #%d didn't have the expected first 4 fields (unique id,state type,first child,# children)",i_state);
			}
			cm->stid[i_state]=(char)stid;
			if (strcmp(sttypeText,"PASSTHRU")==0) {
				cm->sttype[i_state]=PASSTHRU_st;
			}
			else {
				cm->sttype[i_state]=StateCode(sttypeText);
			}

			int isLeftStateInt,isRightStateInt;
			if (fscanf(file,"%d\t%d\t%g\t%g\t",&isLeftStateInt,&isRightStateInt,&(otherStateInfoVector[state].leftwardBeginsc),&(otherStateInfoVector[state].rightwardBeginsc))!=4) {
				throw SimpleStringException("Input profile HMM state #%d didn't have the expected vestigal isLocal fields (isLeft,isRight,leftwardBeginsc,rightwardBegins)",i_state);
			}
			otherStateInfoVector[state].isLeftState=isLeftStateInt!=0;
			otherStateInfoVector[state].isRightState=isRightStateInt!=0;

			size_t numRelatedCmStates;
			if (fscanf(file,"%u\t",&numRelatedCmStates)!=1) {
				throw SimpleStringException("Input profile HMM state #%d didn't have the specified of the # of related CM states",i_state);
			}
			for (size_t i=0; i<numRelatedCmStates; i++) {
				int cmState;
				if (fscanf(file,"%d\t",&cmState)!=1) {
					throw SimpleStringException("Input profile HMM state #%d supposedly had %u related CM states, but they weren't all there (I could only read %d)",i_state,cmState,i);
				}
				hmm2CmStateVector[state].push_back(CovarianceModel::IntToState(cmState));
			}

			if (IsEmittingState(state)) {
				for (int nuc=0; nuc<Alphabet_size; nuc++) {
					if (fscanf(file,"%g\t",&(cm->esc[i_state][nuc]))!=1) {
						throw SimpleStringException("Input profile HMM state #%d was supposedly emitting, but I couldn't read all its 4 emit scores",i_state);
					}
				}
			}

			// tsc
			for (int child=0; child<GetNumChildren(state); child++) {
				float tsc;
				if (fscanf(file,"%g\t",&tsc)!=1) {;
					throw SimpleStringException("Input profile HMM state #%d supposedly has %d children, but I couldn't reaed the transition score for its %dth child",i_state,GetNumChildren(state),child);
				}
				cm->tsc[i_state][child]=tsc;
			}
			fscanf(file,"\n");
		}

		size_t numCmStates;
		if (fscanf(file,"%u\t",&numCmStates)!=1) {
			throw SimpleStringException("Input profile HMM didn't specify the # of CM states in the underlying CM");
		}
		cm2HmmStateVector.resize(numCmStates);
		for (size_t i=0; i<numCmStates; i++) {
			CovarianceModel::State cmState=CovarianceModel::IntToState((int)i);
			fscanf(file,"%d\t%d\n",&(cm2HmmStateVector[cmState].hmmLeftState),&(cm2HmmStateVector[cmState].hmmRightState));
		}
	}
	else {
		// old binary format -- and BTW those 4 bytes we read we the # states in the original format
		int numStates=*(int *)(firstFourChars);

		if (numStates==0) {

			// new format
			int format;
			fread(&format,sizeof(format),1,file);
			loadedFileFormat=format;

			//fprintf(stderr,"Reading InfernalHmm file format #%d\n",format);
			if (format<0 || format>=GetCurrFileFormatNum()+1) {
				Die("Unknown file format for InfernalHmm");
			}

			if (format>=2) {

				// read cmFileName and fullBuildDescription
				size_t s;
				fread(&s,sizeof(s),1,file);
				char *p;
				p=new char[s+1];
				fread(p,1,s,file);
				p[s]=0;
				fromCmFileName=p;
				delete [] p;

				fread(&s,sizeof(s),1,file);
				p=new char[s+1];
				fread(p,1,s,file);
				p[s]=0;
				fullBuildDescription=p;
				delete [] p;
			}

			if (format>=3) {
				fread(&hmmBuildType,sizeof(hmmBuildType),1,file);
			}
			else {
				hmmBuildType=HmmBuildType_Original;
			}

			// convenient to know # of states before format-specific stuff
			fread(&numStates,sizeof(numStates),1,file);
			Init(numStates);

			int int_doLocal;
			fread(&int_doLocal,sizeof(int_doLocal),1,file);
			SetDoLocal(int_doLocal!=0);

			if (DoLocal()) {
				fread(&(endscLinksToLeftVector.front()),sizeof(EndscLinksToLeft),numStates,file);
			}

			fread(&(otherStateInfoVector.front()),sizeof(OtherStateInfo),numStates,file);

			hmm2CmStateVector.resize(numStates);
			for (State state=IntToState(0); state<IntToState(numStates); state++) {
				int n;
				fread(&n,sizeof(n),1,file);
				for (int i=0; i<n; i++) {
					int i_state;
					fread(&i_state,sizeof(i_state),1,file);
					hmm2CmStateVector[state].push_back(IntToState(i_state));
				}
			}

			if (format>=1) {
				size_t numCmStates;
				fread(&numCmStates,sizeof(numCmStates),1,file);
				cm2HmmStateVector.resize(numCmStates);

				fread(&(cm2HmmStateVector.front()),sizeof(Cm2HmmState),numCmStates,file);
			}
		}
		else {
			// original format
			Init(numStates);
		}

		fread(cm->stid,sizeof(cm->stid[0]),numStates+1,file);
		fread(cm->sttype,sizeof(cm->sttype[0]),numStates+1,file);
		fread(cm->cfirst,sizeof(cm->cfirst[0]),numStates,file);
		fread(cm->cnum,sizeof(cm->cnum[0]),numStates,file);

		for (int state=0; state<numStates; state++) {
			fread(cm->esc[state],sizeof(cm->esc[0][0]),Alphabet_size*Alphabet_size,file);
			fread(cm->tsc[state],sizeof(cm->tsc[0][0]),MAXCONNECT,file);
		}
	}

	return true;
}
bool InfernalHmm::LoadInBinary (const char *cmFileName)
{
	// WARNING: only does what HMM needs -- no nodes, null model

	FILE *file=fopen(cmFileName,"rb");
	if (file==NULL) {
		return false;
	}

	bool result=LoadInBinary(file);
	
	fclose(file);

	return result;
}
void InfernalHmm::SaveInDeprecatedBinaryFormat (const char *fileName)
{
	FILE *file=ThrowingFopen(fileName,"wb");
	SaveInDeprecatedBinaryFormat (file);
	fclose(file);
}
void InfernalHmm::SaveInDeprecatedBinaryFormat (FILE *file)
{
	int zero=0;
	fwrite(&zero,sizeof(zero),1,file);

	int format=3;
	fwrite(&format,sizeof(format),1,file);

	size_t s;
	s=fromCmFileName.size();
	fwrite(&s,sizeof(s),1,file);
	fwrite(fromCmFileName.c_str(),1,s,file);
	s=fullBuildDescription.size();
	fwrite(&s,sizeof(s),1,file);
	fwrite(fullBuildDescription.c_str(),1,s,file);

	fwrite(&hmmBuildType,sizeof(hmmBuildType),1,file);

	fwrite(&(cm->M),sizeof(cm->M),1,file);

	int int_doLocal=DoLocal()?1:0;
	fwrite(&int_doLocal,sizeof(int_doLocal),1,file);
	if (DoLocal()) {
		fwrite(&(endscLinksToLeftVector.front()),sizeof(EndscLinksToLeft),GetNumStates(),file);
	}

	fwrite(&(otherStateInfoVector.front()),sizeof(OtherStateInfo),GetNumStates(),file);
	for (State state=IntToState(0); state<IntToState(GetNumStates()); state++) {
		size_t n=hmm2CmStateVector[state].size();
		fwrite(&n,sizeof(n),1,file);
		for (std::list<State>::const_iterator i=hmm2CmStateVector[state].begin(); i!=hmm2CmStateVector[state].end(); i++) {
			int i_state=StateToInt(*i);
			fwrite(&i_state,sizeof(i_state),1,file);
		}
	}

	size_t numCmStates=cm2HmmStateVector.size();
	fwrite(&numCmStates,sizeof(numCmStates),1,file);
	fwrite(&(cm2HmmStateVector.front()),sizeof(Cm2HmmState),numCmStates,file);

	fwrite(cm->stid,sizeof(cm->stid[0]),GetNumStates()+1,file);
	fwrite(cm->sttype,sizeof(cm->sttype[0]),GetNumStates()+1,file);
	fwrite(cm->cfirst,sizeof(cm->cfirst[0]),GetNumStates(),file);
	fwrite(cm->cnum,sizeof(cm->cnum[0]),GetNumStates(),file);

	for (int state=0; state<GetNumStates(); state++) {
		fwrite(cm->esc[state],sizeof(cm->esc[0][0]),Alphabet_size*Alphabet_size,file);
		fwrite(cm->tsc[state],sizeof(cm->tsc[0][0]),MAXCONNECT,file);
	}
}
void InfernalHmm::SaveInBinary (FILE *file)
{
	fprintf(file,"ProfileHMMFilter\n");
	fprintf(file,"FileFormatNum:\t%d\n",GetCurrFileFormatNum());
	fprintf(file,"states:\t%d\n",GetNumStates());
	fprintf(file,"cmFileName:\t%s\n",fromCmFileName.c_str());
	fprintf(file,"hmmBuildCommands:\t%s\n",GetBuildDescription().c_str());
	fprintf(file,"hmmType:\t%s\n",GetHmmBuildTypeDescription(hmmBuildType));
	fprintf(file,"isLocal:\t%d\n",DoLocal()?1:0);

	for (State state=GetFirstState(); state!=GetLastState(); state++) {

		int i_state=StateToInt(state);

		// basic info
		const char *sttypeText;
		if (cm->sttype[i_state]==PASSTHRU_st) {
			sttypeText="PASSTHRU";
		}
		else {
			sttypeText=Statetype(cm->sttype[i_state]);
		}
		fprintf(file,"%d\t%s\t%d\t%d\t",(int)(cm->stid[i_state]),sttypeText,cm->cfirst[i_state],cm->cnum[i_state]);

		// vestigal info for local alignments
		fprintf(file,"%d\t%d\t%g\t%g\t",otherStateInfoVector[state].isLeftState?1:0,otherStateInfoVector[state].isRightState?1:0,otherStateInfoVector[state].leftwardBeginsc,otherStateInfoVector[state].rightwardBeginsc);

		// hmm2cm state mapping (list of CM states that relate to this HMM state)
		size_t numRelatedCmStates=hmm2CmStateVector[state].size();
		fprintf(file,"%u\t",numRelatedCmStates);
		for (std::list<State>::const_iterator i=hmm2CmStateVector[state].begin(); i!=hmm2CmStateVector[state].end(); i++) {
			int i_cmState=StateToInt(*i);
			fprintf(file,"%d\t",i_cmState);
		}

		// esc
		if (IsEmittingState(state)) {
			for (int nuc=0; nuc<Alphabet_size; nuc++) {
				fprintf(file,"%.9g\t",GetSingletEmissionScore(state,nuc));
			}
		}
		// tsc
		for (int child=0; child<GetNumChildren(state); child++) {
			fprintf(file,"%.9g\t",GetNthChildTsc(state,child));
		}
		fprintf(file,"\n");
	}

	// mapping of CM states to HMM left & right states
	fprintf(file,"%u\n",cm2HmmStateVector.size());
	for (size_t i=0; i<cm2HmmStateVector.size(); i++) {
		CovarianceModel::State cmState=CovarianceModel::IntToState((int)i);
		fprintf(file,"%d\t%d\n",InfernalHmm::StateToInt(cm2HmmStateVector[cmState].hmmLeftState),InfernalHmm::StateToInt(cm2HmmStateVector[cmState].hmmRightState));
	}
}
void InfernalHmm::SaveInBinary (const char *cmFileName)
{
	FILE *file=fopen(cmFileName,"wb");
	if (file==NULL) {
		Die("Could not open '%s' for writing.  %s:%d",cmFileName,__FILE__,__LINE__);;
	}

	SaveInBinary(file);

	fclose(file);
}
void InfernalHmm::CopyReverseOf(const InfernalHmm& t)
{
	// WARNING: only for HMMs
	cm=CreateCM(t.GetNumNodes(),t.GetNumStates());

	cm->name=(char *)MallocOrDie(strlen(t.GetName())+1);
	strcpy(cm->name,t.GetName());

	COPY_ARRAY_2D(esc,t.GetNumStates(),Alphabet_size*Alphabet_size);

	COPY_STATE_ARRAY(beginsc);
	endscLinksToLeftVector=t.endscLinksToLeftVector;
	otherStateInfoVector=t.otherStateInfoVector;
	hmm2CmStateVector=t.hmm2CmStateVector;
	cm2HmmStateVector=t.cm2HmmStateVector;
	cm->flags=t.cm->flags;

	// initially no children
	State myState;
	for (myState=GetFirstState(); myState!=GetLastState(); myState++) {
		SetStateType(myState,t.GetStateType(myState));
		SetNoChildren(myState);
	}

	State srcState;
	for (srcState=t.GetFirstState(); srcState!=t.GetLastState(); srcState++) {
		for (int childNum=0; childNum<t.GetNumChildren(srcState); childNum++) {
			State srcToState=t.GetNthChildState(srcState,childNum);
			float tsc=t.GetNthChildTsc(srcState,childNum);
			if (GetNumChildren(srcToState)==0) {
				SetFirstChild(srcToState,srcState);
				SetNumChildren(srcToState,1);
				SetTransitionLogScore(srcToState,0,tsc);
			}
			else {
				int int_srcToState=StateToInt(srcToState);
				int int_srcState=cm->cfirst[int_srcToState] + cm->cnum[int_srcToState];
				assert(srcState==IntToState(int_srcState));
				cm->tsc[int_srcToState][cm->cnum[int_srcToState]]=tsc;
				cm->cnum[int_srcToState]++;
			}
		}
	}
}
int InfernalHmm::GetNumStatesForCmNode (const CovarianceModel& cm,CovarianceModel::Node cmNode) const
{
	// the cheezy, slow implementation
	StateFromNodeList hmmStates;
	GetHmmStatesByCmNode(hmmStates,cm,cmNode);
	return (int)(hmmStates.size());
}
void InfernalHmm::GetHmmStatesByCmNode (StateFromNodeList& get_hmmStates,const CovarianceModel& cm,CovarianceModel::Node cmNode) const
{
	get_hmmStates.clear();

	BoolVectorByCmState hmmStateUsed;
	hmmStateUsed.assign(GetNumStates(),false);

	CovarianceModel::State cmState;
	for (cmState=cm.GetFirstStateOfNode(cmNode); cmState!=cm.GetLastStateOfNode(cmNode); cmState++) {


		InfernalHmm::State hmmLeftState=cm2HmmStateVector[cmState].hmmLeftState;
		if (hmmLeftState!=GetInvalidState()) {
			if (!hmmStateUsed[hmmLeftState]) {
				hmmStateUsed[hmmLeftState]=true;
				StateFromNode add;
				add.state=hmmLeftState;
				add.isRight=false;
				add.isNormal=!cm.IsInsertState(cmState);
				get_hmmStates.push_back(add);
			}
		}

		InfernalHmm::State hmmRightState=cm2HmmStateVector[cmState].hmmRightState;
		if (hmmRightState!=GetInvalidState()) {
			if (!hmmStateUsed[hmmRightState]) {
				hmmStateUsed[hmmRightState]=true;
				StateFromNode add;
				add.state=hmmRightState;
				add.isRight=true;
				add.isNormal=!cm.IsInsertState(cmState);
				get_hmmStates.push_back(add);
			}
		}
	}
}
void InfernalHmm::GetNormalHmmStatesOfLeftOrRightByCmNode (StateList& get_stateList,const CovarianceModel& cm,CovarianceModel::Node cmNode,bool getIsRight) const
{
	get_stateList.clear();

	StateFromNodeList stateFromNodeList;
	GetHmmStatesByCmNode (stateFromNodeList,cm,cmNode);
	for (StateFromNodeList::iterator i=stateFromNodeList.begin(); i!=stateFromNodeList.end(); i++) {

		if (i->isNormal && i->isRight==getIsRight) {
			get_stateList.push_back(i->state);
		}
	}
}
void InfernalHmm::GetHmmEdgesByCmNode (std::list<EdgeInfo>& get_edgeInfoList,const StateFromNodeList& hmmStateList,const CovarianceModel& cm,CovarianceModel::Node cmNode)
{
	BuildNonSavedInfoIfNecessary();

	get_edgeInfoList.clear();

	for (StateFromNodeList::const_iterator i=hmmStateList.begin(); i!=hmmStateList.end(); i++) {

		InfernalHmm::State fromState=i->state;

		bool isRight=i->isRight;

		if (!isRight) {
			assert(IsLeftState(fromState)); // could also be right, but I don't care
			for (int childNum=0; childNum<GetNumChildren(fromState); childNum++) {
				InfernalHmm::State toState=GetNthChildState(fromState,childNum);
				EdgeInfo edge;
				edge.fromState=fromState;
				edge.toState=toState;
				edge.childNum=childNum;
				edge.isRightSide=false;
				get_edgeInfoList.push_back(edge);
			}
		}
		else {
			assert(IsRightState(fromState));

			for (ParentAndMyChildNumVector::const_iterator parentIter=nonSavedInfoVector[fromState].parentAndMyChildNumVector.begin(); parentIter!=nonSavedInfoVector[fromState].parentAndMyChildNumVector.end(); parentIter++) {
				EdgeInfo edge;
				edge.fromState=parentIter->parentState;
				edge.toState=fromState;
				edge.childNum=parentIter->myChildNum;
				edge.isRightSide=true;
				get_edgeInfoList.push_back(edge);
			}
		}
	}
}
void InfernalHmm::CopyProbabilitiesForNode(const InfernalHmm& sourceHmm,const CovarianceModel& cm,const CovarianceModel::Node cmNode)
{
	StateFromNodeList hmmStateList;
	GetHmmStatesByCmNode(hmmStateList,cm,cmNode);
	EdgeInfoList edgeInfoList;
	GetHmmEdgesByCmNode (edgeInfoList,hmmStateList,cm,cmNode);

	// emissions
	for (StateFromNodeList::iterator i=hmmStateList.begin(); i!=hmmStateList.end(); i++) {

		InfernalHmm::State fromState=i->state;

		if (IsEmitting(fromState)) {
			for (int nuc=0; nuc<Alphabet_size; nuc++) {
				SetSingletEmissionLogScore(fromState,nuc,sourceHmm.GetSingletEmissionScore(fromState,nuc));
			}
		}
	}

	// transitions (i.e. edges)
	for (EdgeInfoList::iterator i=edgeInfoList.begin(); i!=edgeInfoList.end(); i++) {
		State fromState=i->fromState;
		State toState=i->toState;
		int childNum=i->childNum;

		SetTransitionLogScore(fromState,childNum,sourceHmm.GetNthChildTsc(fromState,childNum));
	}
}
void InfernalHmm::DumpProbabilitiesCsvByNode(FILE *nodeDump,const CovarianceModel& cm,bool doValues)
{
	BuildNonSavedInfoIfNecessary();

	CovarianceModel::Node cmNode;
	StateFromNodeList hmmStateList;
	EdgeInfoList edgeInfoList;
	for (cmNode=cm.GetFirstNode(); cmNode!=cm.GetLastNode(); cmNode++) {

		if (doValues) {
			fprintf(nodeDump,",,");
		}
		else {
			fprintf(nodeDump,",node #%d,",CovarianceModel::NodeToInt(cmNode));
		}

		GetHmmStatesByCmNode(hmmStateList,cm,cmNode);
		GetHmmEdgesByCmNode (edgeInfoList,hmmStateList,cm,cmNode);

		// emissions
		for (StateFromNodeList::iterator i=hmmStateList.begin(); i!=hmmStateList.end(); i++) {

			InfernalHmm::State fromState=i->state;
			if (IsEmitting(fromState)) {
				for (int nuc=0; nuc<Alphabet_size; nuc++) {
					if (doValues) {
						fprintf(nodeDump,"%f,",GetSingletEmissionScore(fromState,nuc));
					}
					else {
						fprintf(nodeDump,"%d @ %c,",StateToInt(fromState),nucs[nuc]);
					}
				}
			}
		}

		// transitions (i.e. edges)
		for (EdgeInfoList::iterator i=edgeInfoList.begin(); i!=edgeInfoList.end(); i++) {
			State fromState=i->fromState;
			State toState=i->toState;
			int childNum=i->childNum;

			if (doValues) {
				fprintf(nodeDump,"%f,",GetNthChildTsc(fromState,childNum));
			}
			else {
				if (i->isRightSide) {
					std::swap(fromState,toState); // just for printing's sake
				}
				fprintf(nodeDump,"%d --> %d,",StateToInt(fromState),StateToInt(toState));
			}
		}
	}
	fprintf(nodeDump,"\n");
}
bool InfernalHmm::AreSingletEmissionScoresSameForAllNucs(const State& state) const
{
	assert(IsEmitting(state));
	for (int nuc=0; nuc<Alphabet_size; nuc++) {
		if (GetSingletEmissionScore(state,nuc)!=GetSingletEmissionScore(state,0)) {
			return false;
		}
	}
	return true;
}
void InfernalHmm::HmmNormalize (void)
{
	for (State state=GetFirstState(); state!=GetLastState(); state++) {

		if (IsEmitting(state)) {
			double sum=0;
			for (int nuc=0; nuc<Alphabet_size; nuc++) {
				double e=GetSingletEmissionScore(state,nuc);
				sum += pow2(e);
			}
			double logSum=sreLOG2(sum);
			for (int nuc=0; nuc<Alphabet_size; nuc++) {
				SetSingletEmissionLogScore(state,nuc,(float)(GetSingletEmissionScore(state,nuc) - logSum));
			}
		}

		double sum=0;
		for (int child=0; child<GetNumChildren(state); child++) {
			double e=GetNthChildTsc(state,child);
			sum += pow2(e);
		}
		double logSum=sreLOG2(sum);
		for (int child=0; child<GetNumChildren(state); child++) {
			SetTransitionLogScore(state,child,(float)(GetNthChildTsc(state,child)-logSum));
		}
	}
}
void InfernalHmm::ConvertToInfernalSavableFormat (const CovarianceModel& sourceCm)
{
	throw SimpleStringException("Sorry InfernalHmm::ConvertToInfernalSavableFormat doesn't work; Infernal is highly dependent on node type, and it'd take a fair amount of work to make my HMMs conform to the nodes defined by Infernal, & this info is used in Infernal's alignment code.");

	// alg: go thru HMM states in order.  Whenever the corresponding CM node or the left/right side of the HMM state changes, that's a new node

	// first alloc node stuff.  GetNumStates() is an upper bound on # of nodes we'll need
	cm->ndtype = (char *)MallocOrDie(GetNumStates()  * sizeof(char));
	cm->nodemap= (int *)MallocOrDie(GetNumStates()  * sizeof(int));

	int nextNode=-1; // first HMM state will cause an increment, so we'll start at 0
	CovarianceModel::Node prevCmNode=CovarianceModel::GetInvalidNode();
	bool prevIsLeft=false; // just init to avoid compiler warnings
	for (InfernalHmm::State hmmState=GetFirstState(); hmmState!=GetLastState(); hmmState++) {

		// PASSTHRU_st -> S_st
		if (GetStateType(hmmState)==PASSTHRU_st) {
			cm->sttype[InfernalHmm::StateToInt(hmmState)]=S_st;
		}

		// what node are we in
		CovarianceModel::Node cmNode=CovarianceModel::GetInvalidNode();
		for (std::list<CovarianceModel::State>::iterator cmStateIter=hmm2CmStateVector[hmmState].begin(); cmStateIter!=hmm2CmStateVector[hmmState].end(); cmStateIter++) {
			CovarianceModel::State cmState=*cmStateIter;
			CovarianceModel::Node thisNode=sourceCm.GetNode(cmState);
			if (cmNode!=CovarianceModel::GetInvalidNode()) {
				//assert(cmNode==thisNode); // they should all map to the same node, or something's weird
				// this happens with states that are from splicing after bifurcation.  Just take the lowest-# state; it odesn't really matter
				if (thisNode<cmNode) {
					cmNode=thisNode;
				}
			}
			cmNode=thisNode;
		}
		assert(cmNode!=CovarianceModel::GetInvalidNode());
		//printf("%d\n",CovarianceModel::NodeToInt(cmNode));

		// did it change?
		if (prevCmNode!=cmNode || prevIsLeft!=IsLeftState(hmmState)) {
			// a change
			nextNode++;
			cm->nodemap[nextNode]=InfernalHmm::StateToInt(hmmState);
			if (GetStateType(hmmState)==S_st) {
				cm->ndtype[nextNode]=MATL_nd;
			}
			else {
				if (IsEndState(hmmState)) {
					cm->ndtype[nextNode]=END_nd;
				}
				else {
					//crap: on the right side, the nodes are in the reverse order.  Maybe I can just call it MATL_nd... assert(GetStateType(hmmState)==ML_st); // what other type of HMM node is there?
					cm->ndtype[nextNode]=MATL_nd;
				}
			}
		}
		cm->ndidx[InfernalHmm::StateToInt(hmmState)]=nextNode;

		// remember previous
		prevCmNode=cmNode;
		prevIsLeft=IsLeftState(hmmState);
	}

	cm->nodes=nextNode;
	cm->null[0]=cm->null[1]=cm->null[2]=cm->null[3]=1;

	for (int s=0; s<GetNumStates(); s++) {
		if (IsEmitting(IntToState(s))) {
			for (int nuc=0; nuc<Alphabet_size; nuc++) {
				cm->e[s][nuc]=(float)(pow2(cm->esc[s][nuc]));
			}
		}
		for (int child=0; child<cm->cnum[s]; child++) {
			cm->t[s][child]=(float)(pow2(cm->tsc[s][child]));
		}
	}

	// set plast/pnum
	for (InfernalHmm::State hmmState=GetFirstState(); hmmState!=GetLastState(); hmmState++) {
		int i_hmmState=InfernalHmm::StateToInt(hmmState);
		cm->pnum[i_hmmState]=0;
		cm->plast[i_hmmState]=-1;
	}
	for (InfernalHmm::State hmmState=GetFirstState(); hmmState!=GetLastState(); hmmState++) {
		for (int child=0; child<GetNumChildren(hmmState); child++) {
			int i_hmmState=InfernalHmm::StateToInt(hmmState);
			int i_childState=InfernalHmm::StateToInt(GetNthChildState(hmmState,child));
			cm->plast[i_childState]=i_hmmState;
			cm->pnum[i_childState]++;
		}
	}
}
void InfernalHmm::DumpCmAndHmmCorrespondingScores (const char *outFileName,const CovarianceModel& cm)
{
	FILE *out=ThrowingFopen(outFileName,"wt");

	fprintf(out,"CM node,node type,CM state,state type,emit/transition?,emit left nuc OR transition to state#,emit right nuc OR transition to state type,CM score,HMM left score,HMM right score\n");


	for (CovarianceModel::State cmState=cm.GetFirstState(); cmState!=cm.GetLastState(); cmState++) {

		if (cm.IsBifurcation(cmState)) {
			continue;
		}

		fprintf(out,"\n");

		InfernalHmm::State leftHmmState=GetHmmLeftStateOfCmState(cmState);
		InfernalHmm::State rightHmmState=GetHmmRightStateOfCmState(cmState);

		if (cm.IsEmittingState(cmState)) {
			int numLeftNucs=cm.EmitsLeft(cmState) ? Alphabet_size : 1;
			int numRightNucs=cm.EmitsRight(cmState) ? Alphabet_size : 1;
			for (int leftNuc=0; leftNuc<numLeftNucs; leftNuc++) {
				for (int rightNuc=0; rightNuc<numRightNucs; rightNuc++) {
					fprintf(out,"%d,%s,%d,%s,emit,",CovarianceModel::NodeToInt(cm.GetNode(cmState)),Nodetype(cm.GetNodeType(cm.GetNode(cmState))),CovarianceModel::StateToInt(cmState),Statetype(cm.GetStateType(cmState)));
					if (cm.EmitsLeft(cmState)) {
						fprintf(out,"%c",nucs[leftNuc]);
					}
					fprintf(out,",");
					if (cm.EmitsRight(cmState)) {
						fprintf(out,"%c",nucs[rightNuc]);
					}
					fprintf(out,",");
					if (cm.EmitsLeftAndRight(cmState)) {
						fprintf(out,"%f,",cm.GetPairEmissionScore(cmState,leftNuc,rightNuc));
					}
					else {
						const int nuc=cm.EmitsLeft(cmState) ? leftNuc : rightNuc;
						fprintf(out,"%f,",cm.GetSingletEmissionScore(cmState,nuc));
					}
					if (cm.EmitsLeft(cmState)) {
						fprintf(out,"%f",GetSingletEmissionScore(leftHmmState,leftNuc));
					}
					fprintf(out,",");
					if (cm.EmitsRight(cmState)) {
						fprintf(out,"%f",GetSingletEmissionScore(rightHmmState,rightNuc));
					}
					fprintf(out,",");
					fprintf(out,"\n");
				}
			}
		}

		for (int child=0; child<cm.GetNumChildren(cmState); child++) {
			CovarianceModel::State cmToState=cm.GetNthChildState(cmState,child);
			fprintf(out,"%d,%s,%d,%s,transition,%d,%s,%f,",CovarianceModel::NodeToInt(cm.GetNode(cmState)),Nodetype(cm.GetNodeType(cm.GetNode(cmState))),CovarianceModel::StateToInt(cmState),Statetype(cm.GetStateType(cmState)),CovarianceModel::StateToInt(cmToState),Statetype(cm.GetStateType(cmToState)),cm.GetNthChildTsc(cmState,child));

			InfernalHmm::State leftToHmmState=GetHmmLeftStateOfCmState(cmToState);
			if (leftToHmmState!=InfernalHmm::GetInvalidState() && leftHmmState!=InfernalHmm::GetInvalidState()) {
				int hmmChild=GetChildNum_Slow(leftHmmState,leftToHmmState);
				float tsc=GetNthChildTsc(leftHmmState,hmmChild);
				fprintf(out,"%f",tsc);
			}
			fprintf(out,",");

			InfernalHmm::State rightToHmmState=GetHmmRightStateOfCmState(cmToState);
			if (rightToHmmState!=InfernalHmm::GetInvalidState() && rightHmmState!=InfernalHmm::GetInvalidState()) {
				int hmmChild=GetChildNum_Slow(rightToHmmState,rightHmmState);
				float tsc=GetNthChildTsc(rightToHmmState,hmmChild);
				fprintf(out,"%f",tsc);
			}
			fprintf(out,",");
			fprintf(out,"\n");
		}
	}

	fclose(out);
}
CovarianceModel::NodeList InfernalHmm::GetCmNode (const CovarianceModel& cm,State hmmState) const
{
	CovarianceModel::NodeList nodeList;
	const CovarianceModel::StateList& stateList=GetCmState(hmmState);
	for (CovarianceModel::StateList::const_iterator i=stateList.begin(); i!=stateList.end(); i++) {
		CovarianceModel::Node node=cm.GetNode(*i);
		if (std::find(nodeList.begin(),nodeList.end(),node)==nodeList.end()) {
			nodeList.push_back(node);
		}
	}

	return nodeList;
}
void InfernalHmm::ComputePairInflationMatrix (PairInflationMatrix& inflationPerPair,const CovarianceModel& cm,CovarianceModel::State cmConsensusState) const
{
	assert(cm.GetStateType(cmConsensusState)==MP_st); // otherwise it's weird to call this function

	const InfernalHmm::State leftHmmConsensusState=GetHmmLeftStateOfCmState(cmConsensusState);
	const InfernalHmm::State rightHmmConsensusState=GetHmmRightStateOfCmState(cmConsensusState);

	float minInflation=+FLT_MAX; // okay really, the min inflation should be 0, but I'm cheating by calculating this under the assumption that the minimum "inflation" (HMM score minus CM score) is like inflation of 0.  So, later I subtract this out
	for (int leftNuc=0; leftNuc<Alphabet_size; leftNuc++) {
		for (int rightNuc=0; rightNuc<Alphabet_size; rightNuc++) {
			float hmmScore=GetSingletEmissionScore(leftHmmConsensusState,leftNuc) + GetSingletEmissionScore(rightHmmConsensusState,rightNuc);
			inflationPerPair[leftNuc][rightNuc] = hmmScore - cm.GetPairEmissionScore(cmConsensusState,leftNuc,rightNuc);
			minInflation=std::min(minInflation,inflationPerPair[leftNuc][rightNuc]);
		}
	}
	for (int leftNuc=0; leftNuc<Alphabet_size; leftNuc++) {
		for (int rightNuc=0; rightNuc<Alphabet_size; rightNuc++) {
			inflationPerPair[leftNuc][rightNuc] -= minInflation;
		}
	}
}
#ifndef CM2HMM_ONLY
void InfernalHmm::AdjustLeftAndRightEscForRewardedPairsOnly (SingletEmissionScores& leftNucScores,SingletEmissionScores& rightNucScores,const CovarianceModel& cm,const CovarianceModel::State cmConsensusState,const BoolBasePairMatrix& boolBasePairMatrix) const
{
	assert(cm.GetStateType(cmConsensusState)==MP_st); // otherwise it's weird to call this function

	const InfernalHmm::State leftHmmConsensusState=GetHmmLeftStateOfCmState(cmConsensusState);
	const InfernalHmm::State rightHmmConsensusState=GetHmmRightStateOfCmState(cmConsensusState);

	// find inflation matrix
	PairInflationMatrix inflationPerPair;
	ComputePairInflationMatrix(inflationPerPair,cm,cmConsensusState);

	// now actually compute extra penalties
	SingletEmissionScores minInflationForLeftNuc={+FLT_MAX,+FLT_MAX,+FLT_MAX,+FLT_MAX},minInflationForRightNuc={+FLT_MAX,+FLT_MAX,+FLT_MAX,+FLT_MAX};
	for (int leftNuc=0; leftNuc<Alphabet_size; leftNuc++) {
		for (int rightNuc=0; rightNuc<Alphabet_size; rightNuc++) {
			if (boolBasePairMatrix.IsPenalized(leftNuc,rightNuc)) {
				// this is a penalized pair
				minInflationForLeftNuc[leftNuc]=std::min(minInflationForLeftNuc[leftNuc],inflationPerPair[leftNuc][rightNuc]);
				minInflationForRightNuc[rightNuc]=std::min(minInflationForRightNuc[rightNuc],inflationPerPair[leftNuc][rightNuc]);
			}
		}
	}

	for (int nuc=0; nuc<Alphabet_size; nuc++) {

		float score;

		// left
		score=GetSingletEmissionScore(leftHmmConsensusState,nuc);
		score -= minInflationForLeftNuc[nuc]; // if there are no penalized pairs for this leftNuc (kind of weird, but whatever), minInflationForLeftNuc will be +FLT_MAX/2.0, so score will be some large negative #
		if (score<(float)IMPOSSIBLE) {
			score=(float)IMPOSSIBLE;
		}
		leftNucScores[nuc]=score;

		// right
		score=GetSingletEmissionScore(rightHmmConsensusState,nuc);
		score -= minInflationForRightNuc[nuc]; // if there are no penalized pairs for this leftNuc (kind of weird, but whatever), minInflationForLeftNuc will be +FLT_MAX/2.0, so score will be some large negative #
		if (score<(float)IMPOSSIBLE) {
			score=(float)IMPOSSIBLE;
		}
		rightNucScores[nuc]=score;
	}
}
#endif
Cm2Hmm_HmmBuildType InfernalHmm::GetHmmBuildType (void) const
{
	return hmmBuildType;
}
void InfernalHmm::SaveInFormat (const char *hmmFileName,HmmFileFormat hmmFileFormat)
{
	switch (hmmFileFormat) {
		case HmmFileFormat_Binary:
			SaveInDeprecatedBinaryFormat(hmmFileName); // really save in binary
			break;
		case HmmFileFormat_Text:
			SaveInBinary(hmmFileName);
			break;
		case HmmFileFormat_HMMER:
			throw SimpleStringException("Save format HMMER not implemented yet.");
			if (hmmBuildType!=HmmBuildType_Original) {
				throw SimpleStringException("This HMM is not a compact-type profile HMM, so it cannot be saved in HMMER format.");
			}
			break;
		default:
			throw SimpleStringException("Invalid save file format set, or it wasn't implemented.  Sorry, it's undoubtably my fault.");
	}
}

HmmFileFormat InfernalHmm::GetDefaultHmmFileFormat (void)
{
	return HmmFileFormat_Binary;
}

char *GetHmmBuildTypeDescription (Cm2Hmm_HmmBuildType hmmType)
{
	switch (hmmType) {
	case HmmBuildType_Original:
		return "compact";
	case HmmBuildType_separateMPandMLMR:
		return "expanded";
	case HmmBuildType_separateMPMLMRD:
		return "overexpanded"; // deprecated too
	default:
		assertr(false);
		return "unknown type";
	}
}
Cm2Hmm_HmmBuildType GetHmmBuildTypeByText (const char *text)
{
	if (strcmp(text,"compact")==0) {
		return HmmBuildType_Original;
	}
	if (strcmp(text,"expanded")==0) {
		return HmmBuildType_separateMPandMLMR;
	}
	throw SimpleStringException("HMM type text string (e.g. compact or expanded) was not recognized -- is some input file in the wrong format?");
}

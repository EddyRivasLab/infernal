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


void AnnounceSequence(FILE *out,const SQINFO& sqinfo,int sequenceNum)
{
	fprintf(out,"----sequence: #%d,%s %s\n",sequenceNum,sqinfo.name,GetSqinfoDesc(sqinfo));
}
void AnnounceSequenceDirection (FILE *out,bool isReversed)
{
	fprintf(out,"----isReversed: %d\n",isReversed?1:0);
}

////////////////////////
// SequenceSet

SequenceSet::~SequenceSet ()
{
}
std::string SequenceSet::GetEmblId (void)
{
	const char *name=GetSeqName();
	const char EMBL[]="emb|";
	const char SP[]="sp|";
	const int emblSize=sizeof(EMBL)-1;
	const int spSize=sizeof(SP)-1;
	const char *startId=NULL;
	if (strncmp(name,EMBL,emblSize)==0) {
		startId=name+emblSize;
	}
	if (strncmp(name,SP,spSize)==0) {
		startId=name+spSize;
	}
	if (startId!=NULL) {
		// seems right format
		size_t len=strcspn(startId,"|");
		if (len<20) {

			return std::string(startId,len);
		}
		else {
			// that's definitely way too big
		}
	}
	else {
		// not EMBL sequence (well, not an EMBL sequence that was converted to FASTA using the sp2fasta command that's a part of the WU-BLAST package)
		// check if it's NCBI
		const char REF[]="ref|";
		const char *cursor=strstr(name,REF);
		if (cursor!=NULL) {
			cursor += strlen(REF);
			size_t len=strcspn(cursor,"|");
			if (len<20) {
				return std::string(cursor,len);
			}
			else {
				// that's way too big
			}
		}
		else {
			// wrong format -- fall thru
		}
	}

	throw SimpleStringException("(SequenceSet::GetEmblId %s:%d) requested EMBL Id of something that's either not in EMBL format, or I don't know how to parse.  seq description='%s'",__FILE__,__LINE__,name);
}
const char *SequenceSet::GetCurrVirtualFileName (void)
{
	return GetCurrFileName();
}

int SequenceSet_UsingSQINFO::GetLength (void)
{
	return GetSQINFO().len;
}
const char *SequenceSet_UsingSQINFO::GetSqinfoDesc (void)
{
	return ::GetSqinfoDesc(GetSQINFO());
}
const char *SequenceSet_UsingSQINFO::GetSeqName (void)
{
	return GetSQINFO().name;
}


SequenceSet_Reversing::SequenceSet_Reversing (bool onlyForwardStrand_)
: onlyForwardStrand(onlyForwardStrand_)
{
	isReversed=true; // forces getting the next seq
	seq=NULL;
	digitizedSeq=NULL;
}
SequenceSet_Reversing::~SequenceSet_Reversing ()
{
	if (digitizedSeq!=NULL) {
		free(digitizedSeq);
	}
}
char *SequenceSet_Reversing::GetDigitizedSeq  (void)
{
	return digitizedSeq;
}
const char *SequenceSet_Reversing::GetTextSeq (void)
{
	return seq;
}
bool SequenceSet_Reversing::Next (void)
{
	if (digitizedSeq!=NULL) {
		free(digitizedSeq);
	}
	digitizedSeq=NULL;

	if (isReversed || onlyForwardStrand) { // onlyForwardStrand --> skip reversed strand
		isReversed=false;
		if (NextSeq(seq,sqinfo)) {
			/*
			static int loopy=0;
			loopy++;
			if (loopy==47) {
				DebugBreak();
			}
			digitizedSeq=DigitizeSequence(seq, sqinfo.len);
			free(digitizedSeq);
			*/
			digitizedSeq=DigitizeSequence(seq, sqinfo.len);
			AnnounceSequenceDirection(stderr,isReversed);
			return true;
		}
		else {
			return false;
		}
	}
	else {
		isReversed=true;
		revcomp(seq,seq);
		digitizedSeq=DigitizeSequence(seq, sqinfo.len);
		AnnounceSequenceDirection(stderr,isReversed);
		return true;
	}
}
const SQINFO& SequenceSet_Reversing::GetSQINFO (void)
{
	return sqinfo;
}
bool SequenceSet_Reversing::IsReversed (void)
{
	return isReversed;
}
void SequenceSet_Reversing::DumpContext (FILE *out)
{
	AnnounceSequenceDirection(out,isReversed);
}


int SequenceSet_OneFile::squidFormat=SQFILE_FASTA; // default to FASTA -- by forcing this, we can read .gz stuff without mucking around with squid
SequenceSet_OneFile::SequenceSet_OneFile (const char *_fileName,bool onlyForwardStrand)
: SequenceSet_Reversing(onlyForwardStrand)
{
	fileName=_fileName;
	char *nonConstFileName=(char *)_fileName;

	int format = squidFormat;
	if ((sqfp = SeqfileOpen(nonConstFileName, format, NULL)) == NULL) {
		throw SimpleStringException("Failed to open sequence database file %s\n", fileName.c_str());
	}

	AnnounceFastaFile(stderr,fileName.c_str());
	AnnounceFastaFile(stdout,fileName.c_str()); // doesn't change that often

	sequenceNum=-1;
	seq=NULL;
}
SequenceSet_OneFile::~SequenceSet_OneFile ()
{
	if (seq!=NULL) {
		FreeSequence(seq,&sqinfo);
	}
}
bool SequenceSet_OneFile::NextSeq (char *(&get_seq),SQINFO& get_sqinfo)
{
	if (sqfp==NULL) {
		return false;
	}

	if (seq!=NULL) {
		FreeSequence(seq,&sqinfo);
		seq=NULL;
	}

	if (ReadSeq(sqfp, sqfp->format, &seq, &sqinfo)) {
		get_seq=seq;
		get_sqinfo=sqinfo;
		sequenceNum++;
		AnnounceSequence(stderr,sqinfo,sequenceNum);

		char NOTIFY_SOURCE_FASTA[]="NOTIFY_SOURCE_FASTA=";
		if (strncmp(sqinfo.name,NOTIFY_SOURCE_FASTA,strlen(NOTIFY_SOURCE_FASTA))==0) {
			virtualFileName=sqinfo.name+strlen(NOTIFY_SOURCE_FASTA);
		}
		return true;
	}
	else {
		SeqfileClose(sqfp);
		sqfp=NULL;
		return false;
	}
}
const char *SequenceSet_OneFile::GetCurrFileName (void)
{
	return fileName.c_str();
}
const char *SequenceSet_OneFile::GetCurrVirtualFileName (void)
{
	if (virtualFileName.empty()) {
		return GetCurrFileName();
	}
	else {
		return virtualFileName.c_str();
	}
}
void SequenceSet_OneFile::DumpContext (FILE *out)
{
	AnnounceFastaFile(out,fileName.c_str(),virtualFileName.c_str());
	AnnounceSequence(out,sqinfo,sequenceNum);
	SequenceSet_Reversing::DumpContext(out);
}


SequenceSet_FileList::SequenceSet_FileList (const char *fileListFileName,bool onlyForwardStrand_)
: fileListFile(fileListFileName,-1) // -1 is unlikely to be found in a file name...
, onlyForwardStrand(onlyForwardStrand_)
{
	currFile=NULL;
}
SequenceSet_FileList::~SequenceSet_FileList ()
{
	if (currFile!=NULL) {
		delete currFile;
		currFile=NULL;
	}
}
const SQINFO& SequenceSet_FileList::GetSQINFO (void)
{
	return currFile->GetSQINFO();
}
char *SequenceSet_FileList::GetDigitizedSeq (void)
{
	return currFile->GetDigitizedSeq();
}
bool SequenceSet_FileList::IsReversed (void)
{
	return currFile->IsReversed();
}
bool SequenceSet_FileList::Next (void)
{
	while (1) {
		if (currFile!=NULL) {
			bool hasNext=currFile->Next();
			if (hasNext) {
				return true;
			}
			else {
				delete currFile;
				currFile=NULL;
			}
		}

		// see if we can load another file
		assert(currFile==NULL);
		if (!fileListFile.ReadLine()) {
			return false;
		}

		const char *nextFileName=fileListFile.GetField(0);
		if (strlen(nextFileName)==0) {
			printf("file list has 0-length name\n");
			fprintf(stderr,"file list has 0-length name\n");
		}
		else {
			try {
				currFile=new SequenceSet_OneFile(nextFileName,onlyForwardStrand);
			}
			catch (const std::exception& e) {
				printf("ERROR: (SequenceSet_FileList::NextSeq) eating problem opening file: %s\n",e.what());
				fprintf(stderr,"ERROR: (SequenceSet_FileList::NextSeq) eating problem opening file: %s\n",e.what());
			}
		}
	}
}
const char *SequenceSet_FileList::GetCurrFileName (void)
{
	assert(currFile!=NULL);
	return currFile->GetCurrFileName();
}
void SequenceSet_FileList::DumpContext (FILE *out)
{
	currFile->DumpContext(out);
}
const char *SequenceSet_FileList::GetTextSeq (void)
{
	return currFile->GetTextSeq();
}


SequenceSet *MakeSequenceSet(const char *fastaFileString,bool onlyForwardStrand)
{
	if (fastaFileString[0]=='@') {
		fastaFileString++;
		return new SequenceSet_FileList(fastaFileString,onlyForwardStrand);
	}

	return new SequenceSet_OneFile(fastaFileString,onlyForwardStrand);
}

InMemorySequenceSet::InMemorySequenceSet (SequenceSet& sourceSequenceSet)
{
	Seq dummy;
	while (sourceSequenceSet.Next()) {
		seqList.push_back(dummy);
		Seq& seq=seqList.back();
		seq.name=sourceSequenceSet.GetSeqName();
		seq.desc=sourceSequenceSet.GetSqinfoDesc();
		seq.fileName=sourceSequenceSet.GetCurrFileName();
		seq.virtualFileName=sourceSequenceSet.GetCurrVirtualFileName();
		seq.rnaSequence.resize(sourceSequenceSet.GetLength()+1);
		std::copy(sourceSequenceSet.GetDigitizedSeq(),sourceSequenceSet.GetDigitizedSeq()+sourceSequenceSet.GetLength()+1,seq.rnaSequence.begin());
		seq.rnaSequenceLen=sourceSequenceSet.GetLength();
		seq.isReversed=sourceSequenceSet.IsReversed();
	}

	Rewind();
}
InMemorySequenceSet::~InMemorySequenceSet ()
{
}
void InMemorySequenceSet::Rewind (void)
{
	currIter=seqList.end();
}
bool InMemorySequenceSet::Next (void)
{
	if (currIter==seqList.end()) {
		// sneaky, special value for the extra beginning
		currIter=seqList.begin();
		return true;
	}
	else {
		currIter++;
		return currIter!=seqList.end();
	}
}
int InMemorySequenceSet::GetLength (void)
{
	return currIter->rnaSequenceLen;
}
char *InMemorySequenceSet::GetDigitizedSeq (void)
{
	return (char *)(&(currIter->rnaSequence.front()));
}
const char *InMemorySequenceSet::GetTextSeq (void)
{
	throw SimpleStringException("Not implemented %s:%d",__FILE__,__LINE__); // takes lots of RAM & probably won't be used
}
bool InMemorySequenceSet::IsReversed (void)
{
	return currIter->isReversed;
}
const char *InMemorySequenceSet::GetCurrFileName (void)
{
	return currIter->fileName.c_str();
}
const char *InMemorySequenceSet::GetCurrVirtualFileName (void)
{
	return currIter->virtualFileName.c_str();
}
const char *InMemorySequenceSet::GetSqinfoDesc (void)
{
	return currIter->desc.c_str();
}
const char *InMemorySequenceSet::GetSeqName (void)
{
	return currIter->name.c_str();
}
void InMemorySequenceSet::DumpContext (FILE *out)
{
	// do nothing
}

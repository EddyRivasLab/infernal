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

// global vars
char *nucs=Alphabet; // defined in globals.c in infernal-0.54/src
bool cmpWithInfernal=false;
bool cmpWithPureHmm=false;
bool useThresholdOnCm=false;
bool dumpHmmCommitteePruningStats=false;
bool analScoreDumping=false;
const char *overrideHmmCacheFileName=NULL;
const char *dumpHmmScoresFileName=NULL;
const char *dumpLastHeurScoreInSequenceFileName=NULL;
FILE *dumpLastHeurScoreInSequenceFile=NULL;
FILE *dumpHmmScoresFile=NULL;
bool disableHmmBuildInfoDump=false;
const char *dumpFracLetsThruByScoreThresholdFileName=NULL;
bool dumpAverageScore=false;
FILE *dumpFilterResultsFile=NULL;
bool enableProgressFile=true;

/* miscellaneous utility-type functions & classes */

std::string GetFileNameFromFullPath (std::string path)
{
	size_t i=path.find_last_of("/\\");
	if (i==std::string::npos) {
		return path;
	}
	else {
		return std::string(path,i+1,std::string::npos);
	}
}

std::string DumpProgramParams(int argc,char **argv,bool doPreamble)
{
	std::string programParams;
	programParams += "Params: ";
#ifdef _MSC_VER
	int firstArgv=1;
#else
	int firstArgv=0;
#endif
	for (int i=firstArgv; i<argc; i++) {
		programParams += argv[i];
		programParams += " ";
	}
	programParams += "\n";
#ifdef _DEBUG
	programParams += "Build: debug\n";
#else
	programParams += "Build: release\n";
#endif

#ifdef _MSC_VER
	WORD wVersionRequested = MAKEWORD( 2, 2 );
	WSADATA wsaData;
	WSAStartup(wVersionRequested,&wsaData);
#endif
	char hostname[4096]="???";
	gethostname(hostname,4096);
	programParams += "Host: ";
	programParams += hostname;
	programParams += "\n";

	if (doPreamble) {
		printf("%s",programParams.c_str());
#ifndef _WIN32
		pid_t pid=getpid();
		double d_pid=(double)pid;
		printf("pid: %lg\n",d_pid);
#endif
	}

	return programParams;
}

void AnotherParam (void)
{
	Die("insufficient params.  Run \"cmzasha --help\" for help.");
}
bool HasAnotherParam (int a,int argc)
{
	return a<argc;
}
void AnotherParam (int a,int argc)
{
	if (!HasAnotherParam(a,argc)) {
		AnotherParam();
	}
}

std::string LineBreaksToTabs (std::string inStr)
{
	char *s=new char[inStr.size()+1];
	strcpy(s,inStr.c_str());
	char *cursor=s;
	while (true) {
		size_t offset=strcspn(cursor,"\r\n");
		if (offset==strlen(s)) {
			break;
		}
		s[offset]='\t';
	}
	std::string str(s);
	delete [] s;
	return str;
}

std::string LineBreaksToSpaces (std::string inStr)
{
	char *s=new char[inStr.size()+1];
	strcpy(s,inStr.c_str());
	char *cursor=s;
	while (true) {
		size_t offset=strcspn(cursor,"\r\n");
		if (offset==strlen(s)) {
			break;
		}
		s[offset]=' ';
	}
	std::string str(s);
	delete [] s;
	return str;
}


///////////////////////////
// HitList

int HitList::TotalSize (void) const
{
	int total=0;
	const_iterator i;
	for (i=begin(); i!=end(); i++) {
		total += (i->second - i->first);
	}
	return total;
}
cm_int64 HitList::SizeIn2D (int windowLen) const
{
	cm_int64 total=0;
	const_iterator i;
	for (i=begin(); i!=end(); i++) {
		int len=i->second - i->first;
		if (len<=windowLen) {
			total += len*len/2; // it's just a triangular region
		}
		else {
			total += windowLen*windowLen/2; // the triangular starting part
			total += (len-windowLen)*windowLen; // the rectangular part
		}
		total += (i->second - i->first);
	}
	return total;
}
void HitList::Dump (FILE *file) const
{
	const_iterator i;
	for (i=begin(); i!=end(); i++) {
		if (i!=begin()) {
			fprintf(file,", ");
		}
		fprintf(file,"[%d-%d)",i->first,i->second);
	}
	fprintf(file,"\n");
}
int HitList::GetOverallLast (void) const
{
	if (empty()) {
		return 0;
	}
	else {
		return back().second;
	}
}
void HitList::Init (int first,int second)
{
	clear();
	std::pair<int,int> fullWindow;
	fullWindow.first=first;
	fullWindow.second=second;
	push_back(fullWindow);
}
void HitList::Init (int length)
{
	Init(0,length);
}
void HitList::InitEmpty (void)
{
	clear();
}

void GetNucNumsFromHalfOpenInterval(int& startNuc,int& endNuc,int first,int last,int sequenceLen,bool isReversed)
{
	int start=first,end=last;
	// emblcsv format is a fully closed interval, whereas we used half-open intervals
	end--;
	// if reverse strand, we must reflect the #s
	if (isReversed) {
		start=sequenceLen-1-start;
		end=sequenceLen-1-end;
	}
	// emblcsv format is 1-based, but we're using 0-based
	start++;
	end++;

	startNuc=start;
	endNuc=end;
}
void GetNucNumsFromHalfOpenInterval(int& startNuc,int& endNuc,int first,int last,SequenceSet& sequenceSet)
{
	GetNucNumsFromHalfOpenInterval(startNuc,endNuc,first,last,sequenceSet.GetLength(),sequenceSet.IsReversed());
}


////////////////////////////
// FracLetsThruCounter

FracLetsThruCounter::FracLetsThruCounter ()
{
	ResetCounts();
}
FracLetsThruCounter::~FracLetsThruCounter ()
{
}
void FracLetsThruCounter::ResetCounts (void)
{
	nucsInAllSeqs=0;
	nucsLetsThru=0;
	nucsSinceLastProgressReport=0;
	nucsLetThruSinceLastProgressReport=0;

	size2dOfAllSeqs=0;
	size2dLetThru=0;
	size2dSinceLastProgressReport=0;
	size2dLetThruSinceLastProgressReport=0;
}
cm_int64 FracLetsThruCounter::GetNucsInAllSeqs (void)
{
	return nucsInAllSeqs;
}
double FracLetsThruCounter::GetFilteringFraction (void) const
{
	return (double)(nucsLetsThru)/(double)(nucsInAllSeqs);
}
void FracLetsThruCounter::DumpFracLetsThru (FILE *out,const char *messagePrefix,bool sinceLastProgressReport)
{
	fprintf(out,"%sdone %lg nucs (frac lets thru so far=%lg, 2d-fracLetsThru=%lg).",messagePrefix,(double)nucsInAllSeqs,(double)(nucsLetsThru)/(double)(nucsInAllSeqs),(double)(size2dLetThru)/(double)(size2dOfAllSeqs));;
	if (sinceLastProgressReport) {
		fprintf(out,"  (since last report did %lg nucs, fracLetThru=%lg, 2d-fracLetsThru=%lg)\n",(double)(nucsSinceLastProgressReport),(double)(nucsLetThruSinceLastProgressReport)/(double)(nucsSinceLastProgressReport),(double)(size2dLetThruSinceLastProgressReport)/(double)(size2dSinceLastProgressReport));

		nucsSinceLastProgressReport=0;
		nucsLetThruSinceLastProgressReport=0;
		size2dSinceLastProgressReport=0;
		size2dLetThruSinceLastProgressReport=0;
	}
	else {
		fprintf(out,"\n");
	}
}
void FracLetsThruCounter::ProcessPruning (const HitList& inputHitList,const HitList& outputHitList,int windowLen)
{
	int inputNucs=inputHitList.TotalSize();
	int outputNucs=outputHitList.TotalSize();
	nucsInAllSeqs += inputNucs;
	nucsLetsThru += outputNucs;
	nucsSinceLastProgressReport += inputNucs;
	nucsLetThruSinceLastProgressReport += outputNucs;

	cm_int64 input2d=inputHitList.SizeIn2D(windowLen);
	cm_int64 output2d=outputHitList.SizeIn2D(windowLen);
	size2dOfAllSeqs += input2d;
	size2dLetThru += output2d;
	size2dSinceLastProgressReport += input2d;
	size2dLetThruSinceLastProgressReport += output2d;
}



#ifndef DISABLE_ZRAND

zrand::ZRandom *rander=NULL;
zrand::ZRandom *GetRander (void)
{
	return rander;
}
void InitRand (void)
{
	bool initWithRandomSeed=true;
#ifdef _DEBUG
	initWithRandomSeed=false;
#endif
	rander=new zrand::KnuthRngDouble(initWithRandomSeed);
}
void DestroyRand (void)
{
	delete rander;
	rander=NULL;
}
void ReinitRandWithFixedKnuth (void)
{
	delete rander;
	bool initWithRandomSeed=false;
	rander=new zrand::KnuthRngDouble(initWithRandomSeed);
}
double rand01 (void)
{
	/*
	double ri=rand();
	double maxRand=RAND_MAX;
	double r=ri/maxRand;
	assert(r>=0.0 && r<1.0);
	return r;
	*/
	return rander->Get0To1();
}
int RandInt (int N)
{
	int r=(int)(rand01()*(double)(N));
	assert(r>=0 && r<N);
	return r;
}

#endif


std::string GetSubsequence (const char *rnaSequence,int first,int last)
{
	std::string s;
	for (int i=first; i<last; i++) {
		s += nucs[rnaSequence[i]];
	}
	return s;
}

const char *GetSqinfoDesc (const SQINFO& sqinfo)
{
	if ((sqinfo.flags & SQINFO_DESC)!=0) {
		return sqinfo.desc;
	}
	else {
		return "";
	}
}
void AnnounceFastaFile (FILE *out,const char *fastaFileName,const char *virtualFastaFileName)
{
	fprintf(out,"----fastaFile: %s\n",fastaFileName);
	if (strcmp(virtualFastaFileName,"")!=0) {
		fprintf(out,"----virtualFastaFile: %s\n",virtualFastaFileName);
	}
	fflush(out); // doesn't happen that often, since FASTA files are usually reasonably large, & flushing means that the user can see where we got to
}
void AnnounceCmFile (FILE *out,const char *cmFileName,bool doLocalAlignment)
{
	fprintf(out,"----cmFile: %s (%s)\n",cmFileName,doLocalAlignment?"local":"global");
}



void AddWindowToList(HitList& hitList,const std::pair<int,int>& thisWindow)
{
	assert(thisWindow.second>=thisWindow.first);

	if (hitList.empty()) {
		hitList.push_back(thisWindow);
	}
	else {
		std::pair<int,int>& prev=hitList.back();
		if (prev.second >= thisWindow.first) {
			// merge them
			assert(prev.first <= thisWindow.first); // else we weren't calculating the left extent correctly, or it's not sliding properly
			assert(thisWindow.second >= prev.second);
			prev.second=thisWindow.second;
		}
		else {
			// this is a new one
			hitList.push_back(thisWindow);
		}
	}
}


void ParseRfamMembers(RfamSeqList& rfamSeqList,const char *rfamID,const char *rfamFullCsvFileName)
{
	rfamSeqList.clear();

	CommaSepFileReader rfamFullCsv(rfamFullCsvFileName,',');
	while (rfamFullCsv.ReadLine()) {
		if (rfamFullCsv.GetNumFields()!=4) {
			throw SimpleStringException("Rfam.full.csv file (line #%d) has a line with # fields !=4\n",rfamFullCsv.GetLineNum());
		}
		if (strcmp(rfamID,rfamFullCsv.GetField(0))==0 || strcmp(rfamID,rfamFullCsv.GetField(1))==0) {
			RfamSeq thisSeq;
			thisSeq.name=rfamFullCsv.GetField(2);
			thisSeq.seq=rfamFullCsv.GetField(3);
			thisSeq.isFound=false;

			if (thisSeq.seq.size()>=5) { // there's a 1-length member of RF00001...
                rfamSeqList.push_back(thisSeq);
			}
		}
	}
}


void SetNew_defaultTrainingMarkovModelStats (MarkovModelStats *newMarkov)
{
	if (defaultTrainingMarkovModelStats!=NULL) {
		delete defaultTrainingMarkovModelStats;
	}
	defaultTrainingMarkovModelStats=newMarkov;
}


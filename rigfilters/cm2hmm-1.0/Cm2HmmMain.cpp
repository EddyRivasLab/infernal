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
// Cm2HmmMain.cpp: the 'main' function for the cm2hmm command, intended to become a part of Infernal

#include "stdafx.h"
#include <UseDebugNew.h>
#include "cmzasha.h"

void Cm2Hmm_Create(char *cmFileName,const char *hmmFileName,Cm2Hmm_HmmBuildType hmmBuildType,HmmFileFormat hmmFileFormat,SolverWrapper *solverWrapper,const std::string& programParams)
{
	// do one node at a time, as described in paper (doing more doesn't appear to help anyway)
	int numAdjacentNodesToMerge=1;
	int numNodesAtATime=1;
	// generous # of iterations
	int numIters=100000;
	// can't do local alignments
	bool doLocalAlignment=false;

	printf("Using the following 0th-order Markov model:\n");
	defaultTrainingMarkovModelStats->Dump(stdout);


	// run the infinite-length forward alg to create the HMM
	HmmOptimizer_NodeCombiner (cmFileName,doLocalAlignment,programParams,numAdjacentNodesToMerge,numNodesAtATime,numIters,hmmBuildType,solverWrapper,false,hmmFileName,hmmFileFormat);
}

void Cm2Hmm_Create(int a,int argc,char **argv,const std::string& programParams)
{
	bool okay;

	// CM file name
	AnotherParam (a,argc);
	char *cmFileName=argv[a];
	a++;

	// HMM file name
	AnotherParam(a,argc);
	const char *hmmFileName=argv[a];
	a++;

	// markov model specification
	okay=false;
	AnotherParam (a,argc);
	int a_markov=a;
	const char *markovType=argv[a];
	a++;
	if (strcmp(markovType,"uniform")==0) {
		SetNew_defaultTrainingMarkovModelStats(MarkovModelStats::NewUniformMarkov());
		okay=true;
	}
	if (strcmp(markovType,"gc")==0) {
		okay=true;
		AnotherParam (a,argc);
		double gcFrac=atof(argv[a]);
		a++;
		if (gcFrac<0 || gcFrac>1) {
			throw SimpleStringException("G+C fraction used to specify 0th-order Markov model needs to be between 0 and 1, but isn't.");
		}
		double gFrac=gcFrac/2;
		double aFrac=0.5-gFrac;
		double nucProbs[4]={aFrac,gFrac,gFrac,aFrac};
		SetNew_defaultTrainingMarkovModelStats(MarkovModelStats::NewMarkov0 (nucProbs));
	}
	if (strcmp(markovType,"file")==0) {
		okay=true;
		AnotherParam(a,argc);
		const char *fileName=argv[a];
		a++;
		FILE *file=ThrowingFopen(fileName,"rt");
		MarkovModelStats *mm=new MarkovModelStats(file);
		fclose(file);
		SetNew_defaultTrainingMarkovModelStats(mm);
	}
	if (!okay) {
		throw SimpleStringException("Unrecognized specification of 0th-order Markov model: %s",argv[a_markov]);
	}

	// HMM type & output format
	Cm2Hmm_HmmBuildType hmmBuildType;
	HmmFileFormat hmmFileFormat;
	okay=false;
	AnotherParam (a,argc);
	int a_hmm=a;
	const char *hmmType=argv[a];
	a++;
	if (strcmp(hmmType,"compact")==0) {
		okay=true;
		hmmBuildType=HmmBuildType_Original;
		hmmFileFormat=HmmFileFormat_Text;
	}
	if (strcmp(hmmType,"expanded")==0) {
		okay=true;
		hmmBuildType=HmmBuildType_separateMPandMLMR;
		hmmFileFormat=HmmFileFormat_Text;
	}
	if (strcmp(hmmType,"compact-HMMER")==0) {
		okay=true;
		hmmBuildType=HmmBuildType_Original;
		hmmFileFormat=HmmFileFormat_HMMER;
	}
	if (!okay) {
		throw SimpleStringException("Unrecognized specification of HMM output type & format: %s",argv[a_hmm]);
	}

	// solver spec
	SolverWrapper *solverWrapper=MakeSolverWrapper (a,argc,argv);

	Cm2Hmm_Create(cmFileName,hmmFileName,hmmBuildType,hmmFileFormat,solverWrapper,programParams);
}

int try_main(int argc, char **argv)
{
#if defined(_DEBUG) && defined(_MSC_VER)
	// enable MSVC++ debug heap
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF|_CRTDBG_LEAK_CHECK_DF);
	// _CRTDBG_CHECK_ALWAYS_DF|_CRTDBG_DELAY_FREE_MEM_DF
#endif
#ifdef WIN32
	// for this process to allow me to work while it's running
	// set it to low priority
	SetPriorityClass(GetCurrentProcess(),BELOW_NORMAL_PRIORITY_CLASS); // only works on Win 2K/XP
#endif

	bool doHelp=false;

	if (argc>=2) {
		if (strcmp(argv[1],"--help")==0) {
			doHelp=true;
		}
	}
	if (argc<2 || doHelp) {

		fprintf(stderr,"Make a profile HMM for rigorous filtering of a CM:\n");

		fprintf(stderr,"cm2hmm <input CM file name> <output HMM file name> <0th-order Markov model specification> <HMM type & output format> <solver specification>\n");

		fprintf(stderr,"\t<input CM file name> : file name of a CM in Infernal-0.55 format.\n");
		fprintf(stderr,"\t<output HMM file name> : file name of HMM to create.\n");

		fprintf(stderr,"\t<0th-order Markov model specification> : one of the following:\n");
		fprintf(stderr,"\t\tuniform : use a uniform 0th-order model (all nucleotides have probability 0.25)\n");
		fprintf(stderr,"\t\tgc <fraction> : the G+C content is <fraction>, a number from 0 to 1.\n");
		fprintf(stderr,"\t\tfile <file name> : load it from a file (logic to create these files from an input sequence may or may not be implemented in distribution.\n");

		fprintf(stderr,"\t<HMM type & output format> : one of the following: \n");
		fprintf(stderr,"\t\tcompact : create a compact-type profile HMM in the default text format.\n");
		fprintf(stderr,"\t\texpanded : create an expanded-type profile HMM in the default text format.\n");
		//fprintf(stderr,"\t\tcompact-HMMER : create a compact-type profile HMM in HMMER format.\n");

		fprintf(stderr,"\tsolver-specification : one option currently:\n");
		fprintf(stderr,"\t\tcfsqp <B> <C> : use CFSQP.  <B>=0, <C>=1 are reasonable parameters.  Refer to the CFSQP manual for details.\n");
		fprintf(stderr,"\n");

		return 0;
	}

	std::string programParams=DumpProgramParams(argc,argv,true);


	int a=1;
	Cm2Hmm_Create(a,argc,argv,programParams);
	return 0;
}
int main(int argc, char **argv)
{
	int result=1;
	try {
		result=try_main(argc,argv);
	}
	catch (const std::exception& e) {
		printf("FATAL: %s\n",e.what());
	}

	delete defaultTrainingMarkovModelStats;

	return result;
}

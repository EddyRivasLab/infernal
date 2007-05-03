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
// Cm2HmmMainSearch.cpp: the 'main' function for the cm2hmmsearch command, intended to (indirectly) become a part of Infernal

#include "stdafx.h"
#include "UseDebugNew.h"
#include "cmzasha.h"

void AddHmmFilter(int& numHmmFilters,HmmType1 **hmmFilterArray,const char *hmmFileName)
{
	if (strcmp(hmmFileName,"-")==0) {
		// ignore
	}
	else {
		InfernalHmm infernalHmm;
		infernalHmm.LoadInBinary(hmmFileName);
		hmmFilterArray[numHmmFilters]=new HmmType1;
		hmmFilterArray[numHmmFilters]->Init(infernalHmm);
		numHmmFilters++;
	}
}

void ApplyFilters(HitList& hitList,const HitList& inputHitList,float scoreThreshold,const char *rnaSequence,int rnaSequenceLen,int windowLen,int numHmmFilters,HmmType1 **hmmFilterArray)
{
	hitList=inputHitList;
	for (int currHmmFilter=0; currHmmFilter<numHmmFilters; currHmmFilter++) {
		HmmType1& hmm=*(hmmFilterArray[currHmmFilter]);
		//hmm.Dump(stdout);
		//for (int i=0; i<inputHitList.front().second; i++) { printf("%c",nucs[rnaSequence[i]]); }
		CykscanStats scancykStatsDummy;
		HitList currHitList;
		ScanHmm_HmmType1Float_NonTemplated (currHitList,hitList,hmm,scoreThreshold,scancykStatsDummy,rnaSequence,windowLen);
		//hitList.Dump(stdout); currHitList.Dump(stdout);
		hitList=currHitList;
	}
}

void DumpHitHeader(char *rnaSequence,int rnaSequenceLen,int start,int end,int reversed,float score,int& currHitNum)
{
	printf("----hitSequence: ");
	for (int i=start; i<=end; i++) {
		printf("%c",nucs[rnaSequence[i]]);
	}
	printf("\n");

	printf("hit %-4d: %6d %6d (%s) %8.2f bits\n", currHitNum, 
		reversed ? rnaSequenceLen - start + 1 : start, 
		reversed ? rnaSequenceLen - end + 1 : end,
		reversed ? "rev" : "fwd",
		score);
	if (analScoreDumping) {
		unsigned int *ui=(unsigned int *)&score;
		printf("anal score: %.8f (in hex of IEEE format: %x)\n",score,*ui);
	}
	currHitNum++;
}

// for the HMM filtering, we run CYKScan on a subsequence of the overall sequence.  this function translates the sequence position info relative to the subsequence (what we ran CYKScan on), to the overall sequence, so we can report hits
void TranslateHitsFromSubsequence (int nhits,int *hitr,int *hiti,int *hitj,float *hitsc,int subsequenceFirst,int subsequenceLast)
{
	for (int i=0; i<nhits; i++) {

		hiti[i] += subsequenceFirst;
		hitj[i] += subsequenceFirst;
	}
}

void FixupAlignmentByShiftingDatabasePosition(Fancyali_t *ali,int startOfLocalWithinLargerSequence)
{
	ali->sqfrom += startOfLocalWithinLargerSequence;
	ali->sqto += startOfLocalWithinLargerSequence;

	for (int pos=0; pos<ali->len; pos++) {
		if (ali->scoord[pos] != 0) {
			ali->scoord[pos] += startOfLocalWithinLargerSequence;
		}
	}
}
// version of DumpHitsAndAlignments to handle larger sequences (like human chromosomes)
// localRnaSequence,localRnaSequenceLength is within the current interval in which a hit was found; this is a semi-hack that will allow us to handle huge sequences (like human chromosomes) without major changes that relate to Infernal, which would add lots of risk a week before the ISMB deadline (the problem is that the routines to print alignments allocate too much RAM).  This trick won't work if filtering fraction=1, since then we'll be in the same position we started with.  A perhaps more general solution would be to locally break up long sequences into smaller chunks, being careful to overlap chunks to windowLen nucs; the overlap thing seems tricky though.  Also, we still have to do this trick of adjusting the nucleotide positions in our dump to indicate that they came from somewhere else.
void DumpHitsAndAlignments_LocalCoords(int nhits,int *hitr,int *hiti,int *hitj,float *hitsc,char *localRnaSequence,int localRnaSequenceLen,int startOfLocalWithinLargerSequence,int reversed,int rnaSequenceLen,float scoreThreshold,const CovarianceModel& cm,CMConsensus_t *cons,int& currHitNum,SequenceSet *sequenceSet)
{
	Parsetree_t     *tr;		/* parse of an individual hit */
	Fancyali_t      *ali;         /* alignment, formatted for display */

	bool gotHitAboveThreshold=false;

	for (int i = 0; i < nhits; i++)
	{
		bool isAboveThreshold= (hitsc[i] - scoreThreshold) >= -2e-6; // allow some extra fudge
		if (isAboveThreshold) {

			printf ("\n");

			if (!gotHitAboveThreshold) {
				if (sequenceSet!=NULL) {
					sequenceSet->DumpContext(stdout);
				}
				gotHitAboveThreshold=true;
			}

			// this function can know about the full sequence, since it doesn't allocate memory
			DumpHitHeader(localRnaSequence-startOfLocalWithinLargerSequence,rnaSequenceLen,hiti[i],hitj[i],reversed,hitsc[i],currHitNum);

			cm.CYKDivideAndConquer(localRnaSequence, localRnaSequenceLen, CovarianceModel::IntToState(hitr[i]), hiti[i]-startOfLocalWithinLargerSequence, hitj[i]-startOfLocalWithinLargerSequence, &tr);

			ali = cm.CreateFancyAli(tr, cons, localRnaSequence);
			FixupAlignmentByShiftingDatabasePosition(ali,startOfLocalWithinLargerSequence);
			PrintFancyAli(stdout, ali, 
				      0,  /* offset in seq index */
				      reversed); /*TRUE if reverse complement*/

			printf("----endhit\n");

			FreeFancyAli(ali);
			FreeParsetree(tr);
		}

	}

	if (gotHitAboveThreshold) {
		fflush(stdout);
	}

	free(hitr);
	free(hiti);
	free(hitj);
	free(hitsc);
}

void Cm2Hmm_Search(int windowLen,float scoreThreshold,const CovarianceModel& cm,int numHmmFilters,HmmType1 **hmmFilterArray,SequenceSet& sequenceSet,bool runCM)
{
	int maxlen;
	int currHitNum=0;

	int    nhits=0;			/* number of hits in a seq */
	int   *hitr=NULL;			/* initial states for hits */
	int   *hiti=NULL;                  /* start positions of hits */
	int   *hitj=NULL;                  /* end positions of hits */
	float *hitsc=NULL;			/* scores of hits */
	CMConsensus_t   *cons=NULL;	/* precalculated consensus info for display */
	FracLetsThruCounter fracLetsThru;

	Stopwatch_t     *watch;
	watch = StopwatchCreate();

	if (runCM) {
		cons = cm.CreateCMConsensus(3.0, 1.0); 
	}

	StopwatchZero(watch);
	StopwatchStart(watch);

	maxlen   = 0;
	while (sequenceSet.Next()) // iterate over all sequences
	{
		if (sequenceSet.GetLength() == 0) continue; 	/* silently skip len 0 seqs */
		if (sequenceSet.GetLength() > maxlen) maxlen = sequenceSet.GetLength();
		char *rnaSequence = sequenceSet.GetDigitizedSeq();
		const char *rnaSequenceFixed=rnaSequence+1; // my code uses 0-based offsets

		HitList hitList;
		HitList inputHitList;
		inputHitList.Init(sequenceSet.GetLength()); // one interval spanning the whole sequence
		ApplyFilters(hitList,inputHitList,scoreThreshold,rnaSequenceFixed,sequenceSet.GetLength(),windowLen,numHmmFilters,hmmFilterArray);
		fracLetsThru.ProcessPruning(inputHitList,hitList,windowLen); // update filtering fraction statistics

		// The HitList now is a list of intervals that the profile HMM filter was unable to eliminate; scan each of these intervals with the CM
		HitList::iterator i;
		for (i=hitList.begin(); i!=hitList.end(); i++) {
			assert(i->second >= i->first);
			//fprintf(stderr,"[%d,%d) ",i->first,i->second);

			int numNucs=i->second - i->first;
			// do the searches in "local" nucleotide coordinates, i.e. where the start of the interval (i->first) is numbered zero.  This means that for large human chromosomes, we don't have to allocate so much RAM in CYKScan.
			char *localRnaSequence=rnaSequence + i->first;
			int localRnaSequenceLen=i->second - i->first;
			if (runCM) {
				// run the CM
				cm.CYKScanRFWrapper(localRnaSequence, localRnaSequenceLen, windowLen, 
					&nhits, &hitr, &hiti, &hitj, &hitsc);
				// translate the local nucleotide coordinates into coordinates relative to the start of the whole sequence
				TranslateHitsFromSubsequence (nhits,hitr,hiti,hitj,hitsc,i->first,i->second);
				// print out info on any hits that the CM found
				DumpHitsAndAlignments_LocalCoords(nhits,hitr,hiti,hitj,hitsc,localRnaSequence,localRnaSequenceLen,i->first,sequenceSet.IsReversed(),sequenceSet.GetLength(),scoreThreshold,cm,cons,currHitNum,&sequenceSet);
			}
		}
	}

	StopwatchStop(watch);

	printf("\nFound %d hits.\n\n",currHitNum);
	fracLetsThru.DumpFracLetsThru(stdout,"",false);
	StopwatchDisplay(stdout, "\nCPU time: ", watch);

	if (runCM) {
		FreeCMConsensus(cons);
	}
	StopwatchFree(watch);
	SqdClean();
}

void Cm2Hmm_Search(int windowLen,float scoreThreshold,char *cmFileName,const char *hmmFileName1,const char *hmmFileName2,char *sequenceFileName,bool runCM)
{
	CovarianceModel cm;
	cm.Load(cmFileName,false);

	int numHmmFilters=0;
	HmmType1 *hmmFilterArray[2];
	AddHmmFilter(numHmmFilters,hmmFilterArray,hmmFileName1);
	AddHmmFilter(numHmmFilters,hmmFilterArray,hmmFileName2);

	SequenceSet *sequenceSet=MakeSequenceSet(sequenceFileName,false);

	Cm2Hmm_Search(windowLen,scoreThreshold,cm,numHmmFilters,hmmFilterArray,*sequenceSet,runCM);

	delete sequenceSet;
}

void Cm2Hmm_Search(int a,int argc,char **argv,const std::string& programParams)
{
	bool okay;

	// window len
	AnotherParam (a,argc);
	int windowLen=atoi(argv[a]);
	a++;

	// score threshold
	AnotherParam (a,argc);
	float scoreThreshold=(float)(atof(argv[a]));
	a++;

	// CM file name
	AnotherParam (a,argc);
	char *cmFileName=argv[a];
	a++;

	// compact HMM file name
	AnotherParam(a,argc);
	const char *hmmFileName1=argv[a];
	a++;

	// expanded HMM file name
	AnotherParam(a,argc);
	const char *hmmFileName2=argv[a];
	a++;

	// sequence file
	AnotherParam(a,argc);
	char *sequenceFileName=argv[a];
	a++;

	// run CM?
	AnotherParam(a,argc);
	bool runCM=atoi(argv[a])!=0;
	a++;

	Cm2Hmm_Search(windowLen,scoreThreshold,cmFileName,hmmFileName1,hmmFileName2,sequenceFileName,runCM);
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

		fprintf(stderr,"Do a CM search using profile HMM filter(s):\n");

		fprintf(stderr,"cm2hmmsearch <window len> <score threshold> <CM file name> <compact profile HMM file name> <expanded profile HMM file name> <sequence file> <run CM?>\n");

		fprintf(stderr,"\t<window len> : window length parameter for CM scan.\n");

		fprintf(stderr,"\t<score threshold> : hits below this threshold will be ignored (and likely filtered out by the profile HMMs).\n");

		fprintf(stderr,"\t<CM file name> : file name of a CM in Infernal-0.55 format.\n");

		fprintf(stderr,"\t<compact profile HMM file name> : name of a profile HMM to do filtering, or \"-\" (a single dash) to not use this HMM.  Although this HMM is presumed to be compact type, this is not enforced.\n");
		fprintf(stderr,"\t<expanded profile HMM file name> : same idea as previous field.\n");

		fprintf(stderr,"\t<sequence file> : name of a sequence file, which is presumed to be in FASTA format.\n");
		fprintf(stderr,"\t<run CM?> : if \"0\" do NOT actually run the CM, just do the filtering and report the filtering fraction.  If \"1\", run the CM to find hits.\n");

		return 0;
	}

	std::string programParams=DumpProgramParams(argc,argv,true);

	int a=1;
	Cm2Hmm_Search(a,argc,argv,programParams);
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

	return result;
}

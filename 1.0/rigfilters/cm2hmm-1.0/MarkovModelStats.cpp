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


///////////////////////////
// MarkovModelStats

// unnecessary: #include <gsl/gsl_linalg.h>

MarkovModelStats *defaultTrainingMarkovModelStats=NULL;

MarkovModelStats::MarkovModelStats (void)
{
	order=0;
	alphabetSize=0;
}
MarkovModelStats::MarkovModelStats (FILE *file)
{
	LoadFromDump(file);
}
MarkovModelStats::MarkovModelStats (int _order,int alphabetSize_) 
: data(_order+1,alphabetSize_)
{
	alphabetSize=alphabetSize_;
	order=_order;
	path.assign(order+1,0);
}
MarkovModelStats::~MarkovModelStats () 
{
}
void MarkovModelStats::ClearData (void)
{
	data.SetAll(0);
}
void MarkovModelStats::PushOnPath (int nuc) 
{
	assert(nuc>=0 && nuc<alphabetSize);
	path.pop_front();
	path.push_back(nuc);
}
void MarkovModelStats::Inc (void)
{
	VariableDimVector<double>::LowLevelOffset offset=data.GetOffset(path);
	data.Set(offset,data.Get(offset)+1.0);
}
double MarkovModelStats::Get (void) {
	return data.Get(path);
}
void MarkovModelStats::Set (double t) {
	data.Set(path,t);
}
template <class Container>
double MarkovModelStats::GetProbOfNuc_templ (int nucForPr,const Container& container)
{
	assert(nucForPr>=0 && nucForPr<alphabetSize);

	int context[16];
	assert(order<=8);
	if (order>16) {
		Die("%s:%d",__FILE__,__LINE__);
	}

	double totalProb=0;
	double prOfNuc=-1;
	for (int nuc=0; nuc<alphabetSize; nuc++) {

		int *first=context;
		int *lastOfContext=std::copy(container.begin(),container.end(),first);
		*lastOfContext=nuc;
		assert(lastOfContext-first==order); // that's how long the context is

		int *last=lastOfContext+1;
		double p=data.Get(first,last);
		totalProb += p;
		if (nuc==nucForPr) {
			prOfNuc=p;
		}
	}
	assert(prOfNuc>=0);
	if (totalProb==0) {
		return 0.25; // hedge in this case
	}
	else {
		return prOfNuc/totalProb;
	}
}
double MarkovModelStats::GetProbOfNuc (int nucForPr,const NaryCounter& counter)
{
	return GetProbOfNuc_templ(nucForPr,counter);
}
double MarkovModelStats::GetProbOfNuc (int nucForPr,const std::list<int>& context)
{
	assert(context.size()==(size_t)order); // that's how much context an 'order'-order Markov model needs
	return GetProbOfNuc_templ(nucForPr,context);
}
double MarkovModelStats::GetProbOfNuc (int nucForPr,const std::vector<int>& context)
{
	assert(context.size()==(size_t)order); // that's how much context an 'order'-order Markov model needs
	return GetProbOfNuc_templ(nucForPr,context);
}
int MarkovModelStats::GetOrder (void)
{
	return order;
}
double MarkovModelStats::GetProbOfNuc_0order (int nuc)
{
	assert(order==0); // this function is only appropriate for 0-order models, since it doesn't take any context
	std::list<int> context;
	return GetProbOfNuc(nuc,context);
}
void MarkovModelStats::GetContextDistribution (VariableDimVector<double>& probabilityOfEachContext)
{
	probabilityOfEachContext.resize(order,alphabetSize);
	if (order==0) {
		// nothing much to do
		return;
	}

#if 0
	// CRAP -- I was right all along & this whole linear equation thing (which requires the GSL) is unnecessary.  I was bothered by the fact that he order-(N-1) Markov model learned from a sequence is different from the order-(N-1) Markov model inferred from the order-N Markov model learned from the sequence.  It's not the same because when I start with the order-N model, I ignore the first N-1 characters -- this explains why the numbers are close, but not quite the same.  I now believe the non linear algebra code is correct.  (Plus this GSL stuff had some problem, which I don't care to solve.)


	// here's the idea.  For a 1-order MM, where p_xy = Pr(Y is next char | X is prev), we can write the equation
	// p_a*p_aa + p_c*p_ca + p_g*p_ga + p_u*p_ua = p_a
	// writing this out for all cases leads to a system of linear equations
	// we'll also require p_a + p_c + p_g + p_u = 1, otherwise p_x=0 is okay.  To be able to introduce this equation, while keeping the matrix square, we add a dummy variable that's 0.  To make it 0, we place that dummy var in every equation with coefficient 0, and make the rhs of the equation 1.

	// A*x=b, where 'x' is the entries of probabilityOfEachContext
	std::vector<double> A_data,b_data;
	// init A=0, since many entries will be 0
	A_data.assign((size_t)((probabilityOfEachContext.GetLinearSize()+1)*(probabilityOfEachContext.GetLinearSize()+1)),0);
	// b=0 is the final value, except for the last equation
	b_data.assign((size_t)probabilityOfEachContext.GetLinearSize()+1,0);
	b_data.back()=1.0; // result of last equation

	// last row of A is all 1s (sum of vars is 1)
	for (int i=0; i<probabilityOfEachContext.GetLinearSize()+1; i++) {
		A_data[(int)(A_data.size())-1-i]=1;
	}
	// last column of A is all 1s (make the dummy variable 0)
	for (int i=0; i<probabilityOfEachContext.GetLinearSize(); i++) {
		A_data[i*(probabilityOfEachContext.GetLinearSize()+1) + probabilityOfEachContext.GetLinearSize()]=1;
	}
	// but lower,right corner of A is 0 (don't want the last variable influencing)
	A_data.back()=0;
	
	gsl_vector_view b = gsl_vector_view_array(&(*b_data.begin()), probabilityOfEachContext.GetLinearSize()+1);
	gsl_vector *x = gsl_vector_alloc (probabilityOfEachContext.GetLinearSize()+1);

	std::vector<int> currContextInput,currContextOutput;
	currContextInput.reserve(order);
	currContextOutput.reserve(order);

	{
		NaryCounter counter(order-1,alphabetSize);
		bool counting=true;
		while (counting) {

			for (int preNuc=0; preNuc<alphabetSize; preNuc++) {
				for (int postNuc=0; postNuc<alphabetSize; postNuc++) {

					currContextInput.assign(1,preNuc);
					currContextInput.insert(currContextInput.end(),counter.begin(),counter.end());

					currContextOutput.assign(counter.begin(),counter.end());
					currContextOutput.push_back(postNuc);

					std::list<int> currContextInputList;
					currContextInputList.insert(currContextInputList.end(),currContextInput.begin(),currContextInput.end());
					double prob=GetProbOfNuc(postNuc,currContextInputList);
					if (currContextInput==currContextOutput) {
						// same variable, so subtract the occurrence on the rhs
						prob -= 1.0;
					}

					int col=probabilityOfEachContext.LowLevelOffset2LinearOffset(probabilityOfEachContext.GetOffset(currContextInput));
					int row=probabilityOfEachContext.LowLevelOffset2LinearOffset(probabilityOfEachContext.GetOffset(currContextOutput));

					A_data[(size_t)(row*(probabilityOfEachContext.GetLinearSize()+1)+col)]=prob;
				}
			}
			counting=counter.Next();
		}
	}

	gsl_matrix_view A = gsl_matrix_view_array(&(*A_data.begin()), probabilityOfEachContext.GetLinearSize()+1, probabilityOfEachContext.GetLinearSize()+1);

	int s;
	gsl_permutation * p = gsl_permutation_alloc (probabilityOfEachContext.GetLinearSize()+1);
	int status=gsl_linalg_LU_decomp (&A.matrix, p, &s);
	if (status!=GSL_SUCCESS) {
		throw SimpleStringException("gsl_linalg_LU_decomp didn't work: %s",gsl_strerror (status));
	}
	status=gsl_linalg_LU_solve (&A.matrix, p, &b.vector, x);
	if (status!=GSL_SUCCESS) {
		throw SimpleStringException("gsl_linalg_LU_solve didn't work: %s",gsl_strerror (status));
	}

	for (int i=0; i<probabilityOfEachContext.GetLinearSize(); i++) {
		probabilityOfEachContext.Set(probabilityOfEachContext.LinearOffset2LowLevelOffset(i),gsl_vector_get(x,(size_t)i));
	}

	gsl_permutation_free (p);
#endif

#if 1
	probabilityOfEachContext.SetAll(0);

	if (order==0) {
		// nothing more to do
		return;
	}

	// we have the counts of strings of length 'order+1', so each string of length (order+1) has two strings of length (order)
	double totalCount=0;
	{
		NaryCounter counter(order+1,alphabetSize);
		bool counting=true;
		while (counting) {

			double thisCount=data.Get(counter);
			totalCount += thisCount;

			NaryCounter::const_iterator first,last;
			first=counter.begin();
			last=counter.end();

			last--;
			probabilityOfEachContext.GetRef(first,last) += thisCount;

			first++;
			last++;
			probabilityOfEachContext.GetRef(first,last) += thisCount;

			counting=counter.Next();
		}
	}

	// divide total probability out, to make it a probability distribution
	{
		NaryCounter counter(order,alphabetSize);
		bool counting=true;
		while (counting) {

			probabilityOfEachContext.GetRef(counter) /= (totalCount*2.0);

			counting=counter.Next();
		}
	}
#endif
}
void MarkovModelStats::LoadFromDump (const char *fileName)
{
	FILE *file=ThrowingFopen(fileName,"rt");
	LoadFromDump(file);
	fclose(file);
}
void MarkovModelStats::LoadFromDump (FILE *underlyingFile)
{
	CommaSepFileReader file(underlyingFile,',');
	while (file.ReadLine()) {
		if (file.GetNumFields()>0) {
			if (strcmp(file.GetField(0),"order & count-dump list: ")==0) {

				// found the line
				if (file.GetNumFields()<2) {
					throw SimpleStringException("MarkovModelStats::LoadFromDump: found line with data, but it didn't actually have the data.");
				}
				int fieldNum=1;
				order=file.GetFieldAsInt(fieldNum);
				fieldNum++;
				if (order==-1) {
					// new format, which includes alphabetSize
					order=file.GetFieldAsInt(fieldNum);
					fieldNum++;
					alphabetSize=file.GetFieldAsInt(fieldNum);
					fieldNum++;
				}
				else {
					alphabetSize=4; // default, from old file format
				}
				if (alphabetSize==4) { // I don't feel like implementing this for other #s
					if (file.GetNumFields()!=2+(1<<(2*(order+1)))) {
						throw SimpleStringException("MarkovModelStats::LoadFromDump: found line with data, but it didn't have the right amount of data.");
					}
				}

				data.resize(order+1,alphabetSize);

				NaryCounter counter(order+1,alphabetSize);
				bool counting=true;
				while (counting) {

					data.Set(counter,file.GetFieldAsDouble(fieldNum));
					counting=counter.Next();
					fieldNum++;
				}

				// doesn't really matter what the path is
				path.assign(order+1,0);

				return;
			}
		}
	}
	throw SimpleStringException("MarkovModelStats::LoadFromDump: couldn't find data.");
}
void MarkovModelStats::Dump (FILE *out) {
	fprintf(out,"%d-order Markov model:\n",order);

	{
		if (alphabetSize==4) {
			// old format
			fprintf(out,"order & count-dump list: ,%d",order);
		}
		else {
			// new format, specifying alphabetSize
			fprintf(out,"order & count-dump list: ,-1,%d,%d",order,alphabetSize);
		}
		NaryCounter counter(order+1,alphabetSize);
		bool counting=true;
		while (counting) {

			fprintf(out,",%lg",data.Get(counter));

			counting=counter.Next();
		}
		fprintf(out,"\n");
	}
	{
		fprintf(out,"conditional probs:\n");
		NaryCounter counter(order,alphabetSize);
		bool counting=true;
		while (counting) {

			std::list<int> countingPath;
			for (int i=0; i<order; i++) {
				countingPath.push_back(counter[i]);
			}
			for (int nuc=0; nuc<alphabetSize; nuc++) {
				fprintf(out,"\t");
				for (int i=0; i<order; i++) {
					if (alphabetSize==4) {
						fprintf(out,"%c ",nucs[counter[i]]);
					}
					else {
						fprintf(out,"%d ",counter[i]);
					}
				}
				if (alphabetSize==4) {
					fprintf(out,"%c  = %lg\n",nucs[nuc],GetProbOfNuc(nuc,countingPath));
				}
				else {
					fprintf(out,"%d  = %lg\n",nuc,GetProbOfNuc(nuc,countingPath));
				}
			}

			counting=counter.Next();
		}
	}
}
MarkovModelStats *MarkovModelStats::NewMarkov0 (double *nucProbs)
{
	MarkovModelStats *markov=new MarkovModelStats(0);
	markov->alphabetSize=4;
	assert(markov->data.GetDim()==1 && markov->data.GetSizeOfDim(0)==4);
	std::vector<int> offset;
	offset.resize(1);
	for (int i=0; i<4; i++) {
		offset[0]=i;
		markov->data.Set(offset,nucProbs[i]);
	}
	return markov;
}
MarkovModelStats *MarkovModelStats::NewUniformMarkov (void)
{
	double nucProbs[4]={0.25,0.25,0.25,0.25};
	return NewMarkov0(nucProbs);
}
MarkovModelStats *MarkovModelStats::NewDecrementedOrderMarkov (MarkovModelStats *inputMM)
{
	MarkovModelStats *newMM=new MarkovModelStats;
	inputMM->GetContextDistribution (newMM->data);
	newMM->alphabetSize=inputMM->alphabetSize;
	newMM->order=inputMM->order - 1;
	newMM->path.assign(newMM->order+1,0);
	return newMM;
}

#ifndef DISABLE_ZRAND
void MarkovModelStats::GenerateSeq (std::vector<char>& seq,MarkovModelStats& markovModelStats,zrand::ZRandom *rander,int seqLen)
{
	seq.clear();
	seq.reserve(seqLen);

	std::list<int> context;
	markovModelStats.GetRandomContext(context,rander);

	for (int i=0; i<seqLen; i++) {
		int nuc=markovModelStats.GenerateNuc(context,rander);
		seq.push_back((char)nuc);
	}
}
int MarkovModelStats::GenerateNuc (std::list<int>& context,zrand::ZRandom *rander)
{
	double r=rander->Get0To1();
	double probSumSoFar=0;
	int nucToPick=-1;
	for (int nuc=0; nuc<alphabetSize; nuc++) {

		double thisProb=GetProbOfNuc(nuc,context);
		probSumSoFar += thisProb;
		if (r<probSumSoFar) {
			nucToPick=nuc;
			break;
		}
	}
	assert(nucToPick!=-1);
	if (nucToPick==-1) {
		// defensive
		nucToPick=alphabetSize-1; // most likely, we were just the victim of a rounding error, and this is what it should be
	}

	// adjust context (note the order of operation is important for a 0th-order Model, in which case, 'context' must end up empty
	context.push_back(nucToPick);
	context.pop_front();

	return nucToPick;
}
void MarkovModelStats::GetRandomContext (std::list<int>& context,zrand::ZRandom *rander)
{
	if (order==0) {
		// nothing needs to be done -- there's no context
		context.clear();
	}
	else {
		VariableDimVector<double> probabilityOfEachContext;
		GetContextDistribution (probabilityOfEachContext);

		double r=rander->Get0To1();
		double probSumSoFar=0;

		NaryCounter counter(order,alphabetSize);
		bool counting=true;
		while (counting) {

			double thisProb=probabilityOfEachContext.Get(counter);
			probSumSoFar += thisProb;
			if (r<probSumSoFar) {
				context.clear();
				context.insert(context.end(),counter.begin(),counter.end());
				return;
			}

			counting=counter.Next();
		}

		assert(false); // shouldn't get here
		context.assign(order,alphabetSize-1); // defensive programming -- and perhaps it's just rounding error, in which case this is most likely what we should have selected
	}
}
#endif // DISABLE_ZRAND

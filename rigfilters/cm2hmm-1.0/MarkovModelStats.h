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

class MarkovModelStats {
	int order;
	int alphabetSize;
	VariableDimVector<double> data;
	std::list<int> path;

	MarkovModelStats (void);

	template <class Container>
	double GetProbOfNuc_templ (int nucForPr,const Container& container);
public:
	MarkovModelStats (FILE *file); // calls LoadFromDump
	MarkovModelStats (int order_,int alphabetSize_=4);
	~MarkovModelStats ();
	void ClearData (void);
	void PushOnPath (int nuc);
	void Dump (FILE *out);
	void LoadFromDump (FILE *file);
	void LoadFromDump (const char *fileName);
	// get raw number
	double Get (void);
	void Set (double t);
	void Inc (void); // equiv to Set(Get()+1)
	int GetOrder (void);
	// get Pr(nuc|context)
	double GetProbOfNuc (int nucForPr,const std::list<int>& context); // context is a list of 'order' nucs, with context.back() being the most recent
	double GetProbOfNuc (int nucForPr,const NaryCounter& context);
	double GetProbOfNuc (int nucForPr,const vector<int>& context);
	// convenience function for 0-order models
	double GetProbOfNuc_0order (int nucForPr);

	// a Markov model assigns Pr(x|Y) where x is a nuc, and Y is a string of nucs of length 'order'.  I'm calling Y the 'context'.  In some cases, we'd like to find the probability of these contexts according to the model, e.g. to initialize the context for generating a random string according to the model.
	void GetContextDistribution (VariableDimVector<double>& probabilityOfEachContext);

	// makes 0-order Markov model with uniform stats
	static MarkovModelStats *NewUniformMarkov (void);
	// makes 0-order Markov model with given nucleotide probs (nucProbs is array of size 4)
	static MarkovModelStats *NewMarkov0 (double *nucProbs);
	// makes a Markov model of order (order-1), where 'order' is the order of the given MM.  This is pretty much just for testing my code
	static MarkovModelStats *NewDecrementedOrderMarkov (MarkovModelStats *inputMM);

#ifndef DISABLE_ZRAND // I don't want to have to also include this library
	// sets a random context, from the distribution
	void GetRandomContext (std::list<int>& context,zrand::ZRandom *rander);
	// using context, generates random nuc, and also updates context
	int GenerateNuc (std::list<int>& context,zrand::ZRandom *rander);
	// generate a full sequence from the model
	static void GenerateSeq (vector<char>& seq,MarkovModelStats& markovModelStats,zrand::ZRandom *rander,int seqLen);
#endif // DISABLE_ZRAND
};

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
/*
NaryCounter:
Counts in base-N with D digits, i.e. exploring all combinations.
Assumes that N<256, since otherwise going thru all combinations would take too long

  See also 'BooleanCounter.h'

  example code:

  NaryCounter counter(numDigits,base);
  counter.Init();
  bool counting=true;
  while (counting) {

	// use counter[0...(numDigits-1)]

	counting=counter.Next();
  }
*/

class NaryCounter {
protected:
	std::vector<int> array;
	int numDigits;
	int base;
public:
	NaryCounter (int _numDigits,int _base);
	~NaryCounter ();

	// re-initialize to all 0s
	void Init (void);

	// increment by 1
	bool /* has next */ Next (void);

	inline int operator [] (int i) const {
		return array[i];
	}

	// for convenience, make this look like a container
	typedef std::vector<int>::const_iterator const_iterator;
	const_iterator begin (void) const;
	const_iterator end (void) const;
};

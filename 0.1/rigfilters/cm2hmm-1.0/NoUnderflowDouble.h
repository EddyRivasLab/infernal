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
NoUnderflowDouble:
an attempt to generically deal with underflow problems caused
by multiplying <1 probabilities together many times.

  HACK WARNING: this class should work with the kinds of situations that
  I've been having trouble with, but it's not general enough to solve all
  related problems.
*/

class NoUnderflowDouble {
protected:
	double value;
	int extraExp;

	inline static int TrapExpsBelow (void) { return 256; }

	inline void Normalize (void) {
		int e;
		frexp(value,&e); // throw away the mantissa (the return value of frexp)
		if (e>+TrapExpsBelow()) {
			int fudge = ((e/TrapExpsBelow()))*TrapExpsBelow();
			value=ldexp(value,-fudge);
			extraExp += fudge;
		}
		if (e<=-TrapExpsBelow()) {
			int fudge = (-e/TrapExpsBelow())*TrapExpsBelow();
			value=ldexp(value,+fudge);
			extraExp -= fudge;
		}
		int remainder=extraExp%TrapExpsBelow();
		if (remainder!=0) {
			extraExp -= remainder;
			assert((extraExp%TrapExpsBelow())==0);
			value=ldexp(value,+remainder);
		}
	}

	inline NoUnderflowDouble (double _value,int _extraExp) {
		value=_value;
		extraExp=_extraExp;
	}
public:

	// sometimes my implementation is too dinky -- I'd like to know
	class NoUnderflowDoubleOverflow : public SimpleStringException {
	public:
		NoUnderflowDoubleOverflow (const char *file,int line)
			: SimpleStringException("NoUnderflowDouble suffered underflow/overflow at %s:%d.  Excuse: I'm an engineer, not a perfectionist.")
		{
		}
	};

	inline NoUnderflowDouble (void) { }

	inline void operator = (const NoUnderflowDouble& t) { 
		value=t.value; 
		extraExp=t.extraExp; 
	}
	inline void operator = (double t) {
		value=t;
		extraExp=0;
		Normalize();
	}

	inline NoUnderflowDouble (const NoUnderflowDouble& t) { 
		*this=t;
	}
	inline NoUnderflowDouble (double t) {
		*this=t;
	}
	inline NoUnderflowDouble (int t) {
		double d=t;
		*this=d;
	}

	inline operator double () const {
		// this may underflow, but there's nothing else to do
		return ldexp(value,extraExp);
	}

	inline bool IsOverUnderFlowedForDoubles (void) const {
		return !(extraExp>=-1000 && extraExp<=+1000); // not in safe range
	}

	// this function is basically 'operator double',
	// except the caller is saying it's okay to set it to 0.0 on an underflow
	inline double ToDouble_ZeroOnUnderflow (void) const
	{
		if (extraExp<=-768) { // 768 is a bit conservative, but it doesn't really matter
			return 0.0;
		}
		else {
			return (double)*this;
		}
	}

	inline double Log2 (void) const
	{
		double result=log2(value);
		result += extraExp;
		return result;
	}

	// for a weird bug I saw, where the value became negative.  I think this may have been a compiler
	// problem (when I re-built all it went away, but I'm not sure if that's why it went away).
	inline bool IsInClosed0To1 (void) const {
		if (value<0.0) {
			return false;
		}
		int e;
		frexp(value,&e); // throw away the mantissa (the return value of frexp)
		if (e<=0) {
			return true;
		}
		if (e>1) {
			return false;
		}
		// e==1
		const double doubleVal=*this;
		return doubleVal==1.0;
	}

	inline void operator - (void) {
		value=-value;
	}

	inline void operator += (const NoUnderflowDouble& t) {
		if (value==0.0) {
			*this=t;
			return;
		}
		if (t.value==0.0) {
			return;
		}
		if (extraExp==t.extraExp) {
			value += t.value;
		}
		else {
			if (extraExp-t.extraExp==TrapExpsBelow()) {
				value += ldexp(t.value,-TrapExpsBelow());
			}
			else {
				if (extraExp-t.extraExp==-TrapExpsBelow()) {
					extraExp=t.extraExp;
					value=t.value + ldexp(value,-TrapExpsBelow());
				}
				else {
					// pick the more significant one - the other one'll get lost
					if (t.extraExp-extraExp>0) {
						value=t.value;
						extraExp=t.extraExp;
					}
				}
			}
		}
	}

	inline void operator -= (const NoUnderflowDouble& t) {
		// re-use +=, paying a minor performance penalty
		*this += NoUnderflowDouble(-t.value,t.extraExp);
	}

	inline void operator *= (const NoUnderflowDouble& t) {
		extraExp += t.extraExp;
		value *= t.value;
		Normalize();
	}

	inline void operator /= (const NoUnderflowDouble& t) {
		extraExp -= t.extraExp;
		value /= t.value;
		Normalize();
	}

	inline NoUnderflowDouble operator * (const NoUnderflowDouble& t) const {
		NoUnderflowDouble temp(*this);
		temp *= t;
		return temp;
	}

	inline NoUnderflowDouble operator / (const NoUnderflowDouble& t) const {
		NoUnderflowDouble temp(*this);
		temp /= t;
		return temp;
	}

	inline NoUnderflowDouble operator * (double t) const {
		NoUnderflowDouble temp(*this);
		temp *= t;
		return temp;
	}

	inline NoUnderflowDouble operator / (double t) const {
		NoUnderflowDouble temp(*this);
		temp /= t;
		return temp;
	}
	NoUnderflowDouble exp (void) const {
		assert(!IsOverUnderFlowedForDoubles()); // this function is extra hard if the input value could be overflowed (or underflowed), and in the programs I'm interested in, this case doesn't arise.  So, I ignore it, and just assert it's not happening.
		if (IsOverUnderFlowedForDoubles()) {
			throw NoUnderflowDoubleOverflow(__FILE__,__LINE__);
		}
		double td=*this;

		// okay, this is tricky.  First, we agressively decompose the input number, because it's easy to overflow when you're exponentiating something
		int inputExp;
		double inputMantissa=frexp(value,&inputExp);

		if (inputExp<0) {
			// believe it or not, but this case seems tricky, and anyway, there's no possibility of over or under flow, so I'll just special case it
			return NoUnderflowDouble(::exp(td));
		}

		// input now in form: inputMantissa*2^{inputExp}
		// now, e^{input} = e^{inputMantissa*2^{inputExp}} = {e^{inputMantissa}}^{2^{inputExp}}, by some law of exponents

		NoUnderflowDouble expOfMantissa=::exp(inputMantissa); // I'd be very surprised if this overflowed or underflowed, since input mantissa is supposed to be around 1

		// decompose expOfMantissa into its mantissa and exp
		int expOfMantissaExp;
		double expOfMantissaMantissa=frexp(expOfMantissa,&expOfMantissaExp);

		// now then, {e^{inputMantissa}}^{2^{inputExp}} = {expOfMantissa}^{2^{inputExp}}
		//  = {{expOfMantissaMantissa}*2^{expOfMantissaExp}}^{2^{inputExp}}   , substituting our decomposition of expOfMantissa
		//  = {expOfMantissaMantissa}^{2^{inputExp}} * {2^{expOfMantissaExp}}^{2^{inputExp}}  , distributing the exponentiation over the multiplication
		//  Now, {expOfMantissaMantissa}^{2^{inputExp}} should be within range, since I'm assuming the input's exponent wasn't ridiculous, and since expOfMantissaMantissa is a mantissa, so it should be close to 1
		//  And, {2^{expOfMantissaExp}}^{2^{inputExp}} = 2^{ expOfMantissaExp * 2^inputExp }  , by the law of exponents I used earlier, but in reverse
		double expOfMantissaMantissa_part=::exp(::log(expOfMantissaMantissa) * ldexp(1.0,inputExp));
		assert(inputExp>=0); // else this doesn't work, & I'm not sure what to do (which is why I special case it)
		int newExp=expOfMantissaExp * (1<<inputExp);

		NoUnderflowDouble result(expOfMantissaMantissa_part,newExp);
		result.Normalize();
#ifdef _DEBUG
		double direct=::exp(td); // for comparison, to see if my code's right (at least in the case where there's trivially no overflow & we don't really need all this logic)
		double resultAsDouble=result;
		if (resultAsDouble>100) {
			int q=9;
		}
		//printf("%lf,%lf,%lf\n",direct-resultAsDouble,direct,resultAsDouble);
#endif
		return result;
	}
	NoUnderflowDouble log (void) const {
		double valueLog=::log((double)value);
		double extraExpLog=(double)(extraExp)*::log(2.0);
		NoUnderflowDouble result(valueLog+extraExpLog);
#ifdef _DEBUG
		// for comparison, at least when the input is not too high
		double direct=::log(ldexp(value,extraExp));
		double resultAsDouble=result;
#endif
		return result;
	}
};
inline NoUnderflowDouble exp (NoUnderflowDouble t)
{
	//return ::exp((double)t);
	return t.exp();
}
inline NoUnderflowDouble log (NoUnderflowDouble t)
{
	//return ::log((double)t);
	return t.log();
}

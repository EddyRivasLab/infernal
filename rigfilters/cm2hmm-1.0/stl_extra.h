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

//////////////
// miscellaneous things that add to the STL

template <class SourceIter,class Key,class Pred,class Alloc>
void insert (SourceIter first,SourceIter last,std::set<Key,Pred,Alloc>& destSet)
{
	SourceIter i;
	for (i=first; i!=last; i++) {
		destSet.insert(*i);
	}
}

// looks like a vector, but it's at a fixed position -- really more like an array
template <class T>
class FixedPositionVector {
protected:
	T *first,*last;
public:
	FixedPositionVector (void) {
		first=last=NULL;
	}
	FixedPositionVector (T *_first,T *_last) {
		first=_first;
		last=_last;
	}
	FixedPositionVector (vector<T>& t) {
		first=t.begin();
		last=t.end();
	}
	void Set (T *_first,size_t size) {
		first=_first;
		last=first+size;
	}
	FixedPositionVector (T *first,size_t size) {
		Set(first,size);
	}
	void operator = (const FixedPositionVector& t) {
		first=t.first;
		last=t.last;
	}
	FixedPositionVector (const FixedPositionVector& t) {
		*this=t;
	}
	size_t size (void) const {
		return last-first;
	}
	T& operator [] (size_t i) {
		assert(i<size());
		return first[i];
	}
	const T& operator [] (size_t i) const {
		assert(i<size());
		return first[i];
	}
	typedef T * iterator;
	typedef const T * const_iterator;
	iterator begin (void) {
		return first;
	}
	iterator end (void) {
		return last;
	}
	const_iterator begin (void) const {
		return first;
	}
	const_iterator end (void) const {
		return last;
	}
};

// same as above, but data inside is constant
template <class T>
class ConstFixedPositionVector {
protected:
	const T *first,*last;
public:
	ConstFixedPositionVector (void) {
		first=last=NULL;
	}
	ConstFixedPositionVector (const T *_first,const T *_last) {
		first=_first;
		last=_last;
	}
	ConstFixedPositionVector (const vector<T>& t) {
		first=t.begin();
		last=t.end();
	}
	void Set (const T *_first,size_t size) {
		first=_first;
		last=first+size;
	}
	ConstFixedPositionVector (const T *first,size_t size) {
		Set(first,size);
	}
	void operator = (const ConstFixedPositionVector& t) {
		first=t.first;
		last=t.last;
	}
	ConstFixedPositionVector (const ConstFixedPositionVector& t) {
		*this=t;
	}
	void operator = (const FixedPositionVector<T>& t) {
		first=t.begin();
		last=t.begin();
	}
	ConstFixedPositionVector (const FixedPositionVector<T>& t) {
		*this=t;
	}
	size_t size (void) const {
		return last-first;
	}
	const T& operator [] (size_t i) const {
		assert(i<size());
		return first[i];
	}
	typedef const T * const_iterator;
	const_iterator begin (void) const {
		return first;
	}
	const_iterator end (void) const {
		return last;
	}
};

// NOT TESTED -- I decided it wasn't worth it, since the performance problems I'm having with std::list aren't that important right now
// circular array on top of some kind of array type, e.g. ArrayType=vector<T>.  ArrayType must have resize and operator []
// only has some functions so far, because I'm just using it for paths thru a Markov model
template <class T,class ArrayType>
class CircularArrayWithMaxSize {
	friend class iterator;
protected:
	ArrayType array;
	int maxSize;
	int first,currSize;
	// redundant
	int last; // (first+currSize)%maxSize.  Note that we need currSize, since currSize could be 0 or maxSize, which are both 0 mod maxSize

	inline void IncWithinArray (int& offset) {
		offset++;
		if (offset==maxSize) {
			offset=0;
		}
	}
public:
	class OverflowException : public SimpleStringException {
	public:
		OverflowException()
			: SimpleStringException("Overflow in CircularArrayWithMaxSize") {
		}
	};

	CircularArrayWithMaxSize (int _maxSize) {
		maxSize=_maxSize;
		array.resize(maxSize);
		currHead=0;
		currSize=0;
	}
	~CircularArrayWithMaxSize () {
	}

	void push_back (const T& t) {
		assert(currSize<maxSize);
		if (currSize>=maxSize) {
			throw OverflowException();
		}
		array[last]=t;
		IncWithinArray(last);
		currSize++;
	}
	void pop_front (void) {
		assert(currSize>0);
		IncWithinArray(first);
	}

	class const_iterator {
		friend class CircularArray;
	protected:
		CircularArray *circularArray;
		int position,sizeLeft;

		const_iterator (CircularArray *_circularArray,int _position,int _sizeLeft) {
			circularArray=_circularArray;
			position=_position;
			sizeLeft=_sizeLeft;
		}
	public:
		const_iterator (void) {
			circularArray=NULL;
		}
		void operator = (const const_iterator& t) {
			circularArray=t.circularArray;
			position=t.position;
			sizeLeft=t.sizeLeft;
		}
		const_iterator (const const_iterator& t) {
			*this=t;
		}
		~const_iterator () {
		}

		bool operator == (const const_iterator& t) const {
			return circularArray==t.circularArray && position==t.position && sizeLeft=t.sizeLeft;
		}

		void operator ++ (void) {
			circularArray->IncWithinArray(position);
			sizeLeft--;
		}
		void operator ++ (int) {
			++*this;
		}

		const T& operator * (void) const {
			return circularArray->array[position];
		}
	};

	class iterator : public const_iterator {
	};

	const_iterator begin (void) const {
		return const_iterator(*this,first,currSize);
	}
	const_iterator end (void) const {
		return const_iterator(*this,last,0);
	}
};

template <class T>
class CircularArrayWithMaxSizeOverVector : public CircularArrayWithMaxSize<T,vector<T> > {
};

// variable-dimension array, i.e. an array whose dimension is not fixed at compile time
template <class T>
class VariableDimVector {
public:
	typedef int LowLevelOffset; // low-level offset into array for quick lookups
protected:
	vector<T> array;
	int numDim;
	vector<int> sizes;
	static int ProductOfSizes (const vector<int>& sizes) {
		int p=1;
		for (int i=0; i<(int)(sizes.size()); i++) {
			p *= sizes[i];
		}
		return p;
	}

public:

	void resize (const vector<int>& _sizes) {
		sizes=_sizes;
		numDim=(int)(sizes.size());
		array.resize(ProductOfSizes(sizes));
	}
	void resize (int _numDim,int sizeOfEveryDim) {
		numDim=_numDim;
		sizes.assign(numDim,sizeOfEveryDim);
		array.resize(ProductOfSizes(sizes));
	}

	VariableDimVector (void) {
		numDim=0;
	}
	VariableDimVector (const vector<int>& _sizes) {
		resize(_sizes);
	}
	VariableDimVector (int _numDim,int sizeOfEveryDim) {
		resize(_numDim,sizeOfEveryDim);
	}
	~VariableDimVector () {
	}

	int GetDim (void) const {
		return numDim;
	}
	int GetSizeOfDim (int dim) const {
		return sizes[dim];
	}

	int GetLinearSize (void) const {
		return (int)(array.size());
	}
	int LowLevelOffset2LinearOffset (LowLevelOffset i) const {
		return i;
	}
	LowLevelOffset LinearOffset2LowLevelOffset (int i) const {
		return i;
	}

	// sets everything to 't', but doesn't change the dimension or sizes
	void SetAll (const T& t) {
		array.assign(array.size(),t);
	}

	// quick functions if you know the offset
	const T& Get (LowLevelOffset i) const {
		assert(i>=0 && i<(int)(array.size()));
		return array[i];
	}
	T& GetRef (LowLevelOffset i) {
		assert(i>=0 && i<(int)(array.size()));
		return array[i];
	}
	void Set (LowLevelOffset i,const T& t) {
		assert(i>=0 && i<(int)(array.size()));
		array[i]=t;
	}

	// Iter is a forward iterator.  [first,last) controls numDim elements, each of which is an index
	template <class Iter>
	LowLevelOffset GetOffset (const Iter& first,const Iter& last) const {
		int offset=0;
		int dim=0;
		for (Iter i=first; i!=last; i++) {
			assert(dim<numDim); // else [first,last) has too many elements
			offset *= sizes[dim];
			offset += *i;
			dim++;
		}
		assert(dim==numDim); // else [first,last) has too few elements
		return offset;
	}
	// convenience functions, which wrap GetOffset and the low-level Get,Set
	template <class Iter>
	const T& Get (const Iter& first,const Iter& last) const {
		return Get(GetOffset(first,last));
	}
	template <class Iter>
	T& GetRef (const Iter& first,const Iter& last) {
		return GetRef(GetOffset(first,last));
	}
	template <class Iter>
	void Set (const Iter& first,const Iter& last,const T& t) {
		Set(GetOffset(first,last),t);
	}
	// other versions that take a container, rather than two iterators
	template <class Container>
	LowLevelOffset GetOffset (const Container& container) const {
		return GetOffset(container.begin(),container.end());
	}
	template <class Container>
	const T& Get (const Container& container) const {
		return Get(container.begin(),container.end());
	}
	template <class Container>
	T& GetRef (const Container& container) {
		return GetRef(container.begin(),container.end());
	}
	template <class Container>
	void Set (const Container& container,const T& t) {
		Set(container.begin(),container.end(),t);
	}
};

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
Convenience functions for multi-dimensional vectors.
Uses 'vector' from vectorPlus.h (which is basically the same
as the STL vector class)

  Since I can't think of a good generic way to do this, I'm just going to make
  small fixed dimensionalities that I'm likely to use.
*/

template <class T>
class vector2d {
protected:
	vector<vector<T> > vec;
	int sizes[2];
public:
	vector2d () {
	}
	~vector2d () {
	}

	void resize (int firstDim,int secondDim) {
		vec.resize(firstDim);
		int i;
		for (i=0; i<firstDim; i++) {
			vec[i].resize(secondDim);
		}
		sizes[0]=firstDim;
		sizes[1]=secondDim;
	}

	void assign (int firstDim,int secondDim,const T& t) {
		resize(firstDim,secondDim);
		int i,j;
		for (i=0; i<size(0); i++) {
			for (j=0; j<size(1); j++) {
				(*this)[i][j]=t;
			}
		}
	}

	int size (int dimension) const {
		assert(dimension>=0 && dimension<2);
		return sizes[dimension];
	}

	class Dim2Ref {
	protected:
		vector2d<T>& owner;
		int firstIndex;
	public:
		Dim2Ref (vector2d<T>& _owner,int _firstIndex)
			: owner(_owner),firstIndex(_firstIndex)
		{
		}
		T& operator [] (int secondIndex) {
			return owner.vec[firstIndex][secondIndex];
		}
	};
	class const_Dim2Ref {
	protected:
		const vector2d<T>& owner;
		int firstIndex;
	public:
		const_Dim2Ref (const vector2d<T>& _owner,int _firstIndex)
			: owner(_owner),firstIndex(_firstIndex)
		{
		}
		const T& operator [] (int secondIndex) const {
			return owner.vec[firstIndex][secondIndex];
		}
	};

	Dim2Ref operator [] (int firstIndex) {
		return Dim2Ref(*this,firstIndex);
	}
	const_Dim2Ref operator [] (int firstIndex) const {
		return const_Dim2Ref(*this,firstIndex);
	}
	const T& Get (int firstIndex,int secondIndex) const {
		return (*this)[firstIndex][secondIndex];
	}
	void Set (int firstIndex,int secondIndex,const T& t) const {
		(*this)[firstIndex][secondIndex]=t;
	}

	friend class Dim2Ref;
	friend class const_Dim2Ref;
};



template <class T>
class MultiplyArray3d {
protected:
	int s1,s2,s3;
	vector<T> array;
	inline int GetOffset (int i1,int i2,int i3) const {
		return i1*s2*s3 + i2*s3 + i3;
	}
public:
	MultiplyArray3d (void) {
		s1=s2=s3=0;
	}
	void Init (int _s1,int _s2,int _s3) {
		s1=_s1; s2=_s2; s3=_s3;
		array.resize(s1*s2*s3);
	}
	MultiplyArray3d (int _s1,int _s2,int _s3) {
		Load(_s1,_s2,_s3);
	}
	~MultiplyArray3d () {
	}
	const T& Get (int i1,int i2,int i3) const {
		assert(i1>=0 && i1<s1);
		assert(i2>=0 && i2<s2);
		assert(i3>=0 && i3<s3);
		return array[GetOffset(i1,i2,i3)];
	}
	void Set (int i1,int i2,int i3,const T& t) {
		assert(i1>=0 && i1<s1);
		assert(i2>=0 && i2<s2);
		assert(i3>=0 && i3<s3);
		array[GetOffset(i1,i2,i3)]=t;
	}
};

template <class T>
class MultiplyArray4d {
protected:
	int s1,s2,s3,s4;
	vector<T> array;
	inline int GetOffset (int i1,int i2,int i3,int i4) const {
		return i1*s2*s3*s4 + i2*s3*s4 + i3*s4 + i4;
	}
public:
	MultiplyArray4d (void) {
		s1=s2=s3=s4=0;
	}
	void Init (int _s1,int _s2,int _s3,int _s4) {
		s1=_s1; s2=_s2; s3=_s3; s4=_s4;
		array.resize(s1*s2*s3*s4);
	}
	MultiplyArray4d (int _s1,int _s2,int _s3,int _s4) {
		Init(_s1,_s2,_s3,_s4);
	}
	~MultiplyArray4d () {
	}
	const T& Get (int i1,int i2,int i3,int i4) const {
		assert(i1>=0 && i1<s1);
		assert(i2>=0 && i2<s2);
		assert(i3>=0 && i3<s3);
		assert(i4>=0 && i4<s4);
		return array[GetOffset(i1,i2,i3,i4)];
	}
	void Set (int i1,int i2,int i3,int i4,const T& t) {
		assert(i1>=0 && i1<s1);
		assert(i2>=0 && i2<s2);
		assert(i3>=0 && i3<s3);
		assert(i4>=0 && i4<s4);
		array[GetOffset(i1,i2,i3,i4)]=t;
	}

	void Load (FILE *file) {
		fread(&s1,sizeof(s1),1,file);
		fread(&s2,sizeof(s1),1,file);
		fread(&s3,sizeof(s1),1,file);
		fread(&s4,sizeof(s1),1,file);
		Init(s1,s2,s3,s4);

		fread(&*(array.begin()),sizeof(T),array.size(),file);
	}
	void Save (FILE *file) {
		fwrite(&s1,sizeof(s1),1,file);
		fwrite(&s2,sizeof(s1),1,file);
		fwrite(&s3,sizeof(s1),1,file);
		fwrite(&s4,sizeof(s1),1,file);

		fwrite(&*(array.begin()),sizeof(T),array.size(),file);
	}
};

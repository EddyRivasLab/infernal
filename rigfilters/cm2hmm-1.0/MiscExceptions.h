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
Miscellaneous exception classes derived from std::exception

  Requires:
  <exception>,<string>  (STD C++)
*/

#ifndef MISC_EXCEPTIONS_INCLUDED
#define MISC_EXCEPTIONS_INCLUDED

// stores the message as a string
class SimpleStringException : public std::exception {
protected:
	std::string msg;

	// NOTE: you're not supposed to throw pointers to these exceptions (unlike my previous
	// exception classes), so use of operator new probably means that I'm being absent-minded
	// thus, I'm making this non-public
protected: void * operator new (size_t bytes);

public:
	SimpleStringException (const char *format,...);
	SimpleStringException (const std::string& s);
	SimpleStringException (const SimpleStringException& t);
	~SimpleStringException () throw ();

    const char *what() const throw();
};

// 'fopen' call failed
class FopenException : public SimpleStringException {
protected:
	static std::string BuildErrorMessage(const char *fileName);
public:
	FopenException (const char *fileName);
	~FopenException () throw ();
};
// convenience: throws exception on failure
extern FILE *ThrowingFopen (const char *fileName,const char *mode);

// some ANSI C call failed
class ANSICLibException : public SimpleStringException {
protected:
	static std::string BuildErrorMessage(const char *description,const char *failedFunctionName);
public:
	ANSICLibException (const char *description,const char *failedFunctionName);
	~ANSICLibException () throw ();
};

// convenience function for ANSI-C errors, using the errno/strerror interface
std::string GetAnsiCErrorMessage (void);
// same, for WIN32 errors, using GetLastError
std::string GetWin32ErrorMessage (void);

// assert, even in release mode; throw exception on failure
#ifdef _DEBUG
// in debug mode, use regular assert
#define assertr(exp) assert(exp)
#else
#define assertr(exp) if (!(exp)) { throw SimpleStringException("Internal error (release mode assertion failed \"%s\") %s:%d",#exp,__FILE__,__LINE__); }
#endif


#endif


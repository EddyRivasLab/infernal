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
  (c) Zasha Weinberg, 1999-2002.
  Released into the public domain.
 */

#include "stdafx.h"

#include "MiscExceptions.h"

#include <errno.h>
#include <stdlib.h>
#include <stdarg.h>

#ifdef _MSC_VER
#define vsnprintf _vsnprintf
#include <malloc.h>
#endif

#ifdef WIN32
#include <windows.h>
#endif

////////////////////
// SimpleStringException

SimpleStringException::SimpleStringException (const char *format,...)
{
    va_list(arglist);
    va_start(arglist, format);

	vector<char> buf;
	buf.resize(128); // starting value.  Actually, I'm not too concerned about having too many
		// reallocations, because this is an exception, so it shouldn't happen too often
		// and as for memory allocation errors, I've completely punted on that.

	while (1) {
		int snprintfReturn=vsnprintf(&*(buf.begin()),buf.size(),format,arglist);

		if (snprintfReturn>=0 && (size_t)(snprintfReturn)+1+1<=buf.size()) { // extra +1 is for snprintf on g++ 2.95.3/SGI, which seems to return the # of bytes that fit when not the whole string fit (very silly).  So, for its benefit, we conservatively try again in that case
			// vsnprintf says it fit in buf, so we're done
			break;
		}

		// buffer wasn't long enough

		// MSVC++ returns <0 for too long condition
		if (snprintfReturn<0 || snprintfReturn+1==buf.size()) { // again, last case is for g++ 2.95.3/SGI
			// double buffer & try again
			buf.resize(buf.size()*2);
		}
		else {
			assert(snprintfReturn>=buf.size());
			// C99 standard and relatively recent gcc's return the # of bytes that are needed in buffer (much more useful!)
			// alloc buf to required size, & try again
			buf.resize(snprintfReturn+1); // not sure if +1 is necessary, but can't hurt much...
		}
	}
	msg=&*(buf.begin());

	va_end(arglist);
}
SimpleStringException::SimpleStringException (const std::string& s)
: msg(s)
{
}
SimpleStringException::SimpleStringException (const SimpleStringException& t)
{
	msg=t.msg;
}
SimpleStringException::~SimpleStringException () throw ()
{
}
const char *SimpleStringException::what() const throw()
{
	return msg.c_str();
}

std::string GetAnsiCErrorMessage (void)
{
	return std::string(strerror(errno));
}
std::string GetWin32ErrorMessage (void)
{
#ifdef WIN32
	DWORD error=GetLastError();
	char *buffer=NULL;
	if (FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER|FORMAT_MESSAGE_FROM_SYSTEM|FORMAT_MESSAGE_IGNORE_INSERTS,
		NULL,error,0,(char *)&buffer,1,NULL)==0) {
		char buf[1024];
		sprintf(buf,"Error code %d (sorry FormatMessage failed with error code %d)",
			error,GetLastError());
		return std::string(buf);
	}
	else {
		std::string msg(buffer);
		LocalFree(buffer);
		return msg;
	}
#else
	return std::string("GetWin32ErrorMessage ("__FILE__") called, but we're not on WIN32");
#endif
}


//////////////////////
// FileNotFoundException

FopenException::FopenException (const char *fileName)
: SimpleStringException(BuildErrorMessage(fileName))
{
}
FopenException::~FopenException () throw ()
{
}
std::string FopenException::BuildErrorMessage(const char *fileName)
{
	static const char format[]="Cannot open file '%s': '%s'";
	std::string errorType=GetAnsiCErrorMessage();

	char *msg=(char *)alloca(strlen(format)+strlen(errorType.c_str())+strlen(fileName));
	sprintf(msg,format,fileName,errorType.c_str());
	return std::string(msg);
}

FILE *ThrowingFopen (const char *fileName,const char *mode)
{
	FILE *file;
	file=fopen(fileName,mode);
	if (file==NULL) {
		throw FopenException(fileName);
	}
	return file;
}

////////////////////////////
// ANSICLibException
ANSICLibException::ANSICLibException (const char *description,const char *failedFunctionName)
: SimpleStringException(BuildErrorMessage(description,failedFunctionName))
{
}
ANSICLibException::~ANSICLibException () throw ()
{
}
std::string ANSICLibException::BuildErrorMessage(const char *description,const char *failedFunctionName)
{
	static const char format[]="%s: '%s' failed, '%s'";
	std::string errorType=GetAnsiCErrorMessage();

	char *msg=(char *)alloca(strlen(format)+strlen(errorType.c_str())+strlen(description)+strlen(failedFunctionName));
	sprintf(msg,format,description,failedFunctionName,errorType.c_str());
	return std::string(msg);
}

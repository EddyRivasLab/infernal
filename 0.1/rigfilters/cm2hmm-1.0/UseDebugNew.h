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
#include this file in .cpp files immediately after the main #includes
Works also in template implementation files & will change file name

NOTE: code to enable debug head is:
#if defined(_DEBUG) && defined(WIN32)
	// do this after initialization - most of the above
	// has the lifetime of the program, so it's not a big deal
	// if it's cleaned, or not (and it's a pain to sift thru it all)
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF|_CRTDBG_LEAK_CHECK_DF);
#endif
*/

#if defined(_DEBUG) && defined(WIN32)
// Memory Leak detection stuff.  Thanks VC++.
#undef DEBUG_NEW
#define DEBUG_NEW new(_NORMAL_BLOCK,__FILE__,__LINE__)
#define new DEBUG_NEW
#define malloc(X) _malloc_dbg(X,_NORMAL_BLOCK,__FILE__,__LINE__)
#else
#define DEBUG_NEW new
#endif

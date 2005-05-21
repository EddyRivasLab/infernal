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
CommaSepFileReader:
Reads comma-, tab-, or any-character separated files
for slightly extra convenience.
*/
#ifndef __COMMASEPFILEREADER_H
#define __COMMASEPFILEREADER_H

// a nice class for SGD files, which have internal fields that are pipe-delimited
// this class is also used by the actual CommaSepFileReader
// This class will split up given strings (in a vector<char>) by the delimiter, and
// allow the user to query the fields
class CommaSepSeparator {
protected:
	const char delimiterChar;
	std::vector<const char *> fieldsInCurrLine;
	std::vector<char> currLine;

	virtual std::string GetAdditionalInformationForException (void) const;
	void SeparateCurrLine(void);
public:
	CommaSepSeparator (char _delimiterChar);
	~CommaSepSeparator ();

	CommaSepSeparator (const CommaSepSeparator& t);

	// copies line to an internal vector, and separates the fields in it
	// this replaces any line that was previoualy separated (i.e. all previous fields are lost)
	void SeparateLine (const std::vector<char>& line);
	// same thing, but here line is a 0-terminated string
	void SeparateLine (const char *line);
	void SeparateLine (const std::string& line);

	int GetNumFields (void) const;
	const char *GetField (int fieldNum) const;

	// for convenience: returns # of fields, not counting blank fields on the right (when exporting
	// to .csv format, Excel will include blank fields to pad lines that don't have so many fields)
	int GetNumFieldsExcludingBlankPadding (void) const;
	// for convenience: returns true iff original line was blank
	bool IsLineBlank (void) const;
	// for convenience, returns true iff field is the empty string
	bool IsFieldBlank (int fieldNum) const;
	// convenience, returns true iff at least 1 field has the value of valueStr
	bool FieldsContainValue (const char *value) const;

	// convenience: getting field as other data types
	int GetFieldAsInt (int fieldNum) const; // note: it's an error for field to be blank, or contain any non-int characters
	double GetFieldAsDouble (int fieldNum) const; // note: it's an error for field to be blank, or contain any non-double characters
};

class CommaSepAbstractFile : public CommaSepSeparator {
public:
	CommaSepAbstractFile (char delimiterChar);
	virtual ~CommaSepAbstractFile ();

	virtual int GetLineNum (void) const = 0;
	virtual bool /* has next line */ ReadLine (void) = 0;
	// gets the 0-based field in the current line
	// for convenience, throws an exception if there aren't that many fields in the curr
	// line (as opposed to asserting, which wouldn't work in release mode & thus isn't appropriate
	// for possible errors with input files)
};

class CommaSepFileReader : public CommaSepAbstractFile {
protected:
	FILE *inFile;
	bool deleteFileOnDestructor; // does this class own the file
	int lineNum;

	std::vector<char> currLine;

	// implement this to put in line #s
	std::string GetAdditionalInformationForException (void) const;
	// make these inherited members protected
	inline void SeparateLine (const std::vector<char>& line) { CommaSepSeparator::SeparateLine(line); }
	void SeparateLine (const char *line);
public:
	CommaSepFileReader (const char *fileName,char _delimiterChar);
	CommaSepFileReader (FILE *_inFile,char _delimiterChar,int currLineNum=0); // will not close file on exit if this constructor used
	~CommaSepFileReader ();

	bool ReadLine (void);
	int GetLineNum (void) const;

	// gets a CommaSepSeparator class that corresponds to the current line, so the caller
	// can store this data conveniently
	const CommaSepSeparator& GetCommaSepSeparator (void) const;
};

// similar to reading a comma-separated file, but the entire file is just one string.  One delimeter separates lines, another separates fields within the lines.  This is useful for code that reads a comma-sep file, when I want to optionally be able to put the file as a command-line parameter, so I don't have to create a whole file.
class CommaSepMetaSep : public CommaSepAbstractFile {
protected:
	CommaSepSeparator lines;
	int lineNum;
public:
	CommaSepMetaSep (const char *fullString,char lineDelimiterChar='/',char fieldDelimiterChar=',');
	~CommaSepMetaSep ();
	bool ReadLine (void);
	int GetLineNum (void) const;
};

#endif // __COMMASEPFILEREADER_H


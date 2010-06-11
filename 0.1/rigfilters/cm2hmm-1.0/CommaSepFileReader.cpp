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

#include "CommaSepFileReader.h"
#include "MiscExceptions.h"

////////////////////
// CommaSepSeparator

CommaSepSeparator::CommaSepSeparator (char _delimiterChar)
: delimiterChar(_delimiterChar)
{
}
CommaSepSeparator::CommaSepSeparator (const CommaSepSeparator& t)
: delimiterChar(t.delimiterChar)
{
	currLine=t.currLine;
	fieldsInCurrLine=t.fieldsInCurrLine;
	size_t i;
	for (i=0; i<fieldsInCurrLine.size(); i++) {
		fieldsInCurrLine[i] -= (long)&(*t.currLine.begin());
		fieldsInCurrLine[i] += (long)&(*currLine.begin());
	}
}
CommaSepSeparator::~CommaSepSeparator ()
{
}
void CommaSepSeparator::SeparateCurrLine(void)
{
	fieldsInCurrLine.clear();
	char *cursor=&*currLine.begin();
	while (1) {
		char *nextTab=strchr(cursor,delimiterChar);
		if (nextTab!=NULL) {
			*nextTab=0;
		}
		fieldsInCurrLine.push_back(cursor);
		if (nextTab==NULL) {
			break;
		}
		cursor=nextTab+1;
	}
}
void CommaSepSeparator::SeparateLine (const char *line)
{
	const char *first(line);
	const char *last(line+strlen(line)+1);
	currLine.reserve(last-first);
	currLine.clear();
	const char *i;
	for (i=first; i!=last; i++) {
		currLine.push_back(*i);
	}
	SeparateCurrLine();
}
void CommaSepSeparator::SeparateLine (const std::vector<char>& line)
{
	currLine=line;
	SeparateCurrLine();
}
void CommaSepSeparator::SeparateLine (const std::string& line)
{
	SeparateLine(line.c_str());
}
const char *CommaSepSeparator::GetField (int fieldNum) const
{
	if (fieldNum<0 || fieldNum>=GetNumFields()) {
		throw SimpleStringException("Reading delimited data file: a required field was missing (0-based field #%d, %s) .",
			fieldNum,GetAdditionalInformationForException().c_str());
	}
	return fieldsInCurrLine[fieldNum];
}
int CommaSepSeparator::GetNumFields (void) const
{
	return fieldsInCurrLine.size();
}
bool CommaSepSeparator::IsFieldBlank (int fieldNum) const
{
	return GetField(fieldNum)[0]==0;
}
int CommaSepSeparator::GetNumFieldsExcludingBlankPadding (void) const
{
	int numFields=GetNumFields();
	while (numFields>0) {
		if (!IsFieldBlank(numFields-1)) {
			break;
		}
		numFields--;
	}
	return numFields;
}
bool CommaSepSeparator::IsLineBlank (void) const
{
	switch (GetNumFields()) {
	case 0:
		return true;
	case 1:
		return GetField(0)[0]==0;
	default:
		return false;
	}
}
std::string CommaSepSeparator::GetAdditionalInformationForException (void) const
{
	return std::string("");
}
int CommaSepSeparator::GetFieldAsInt (int fieldNum) const
{
	const char *field=GetField(fieldNum);
	char *endptr;
	int result=strtol(field,&endptr,10);
	if (*endptr!=0) {
		throw SimpleStringException("Int field had some non-numeric content, field text='%s', %s",
			field,GetAdditionalInformationForException().c_str());
	}
	return result;
}
double CommaSepSeparator::GetFieldAsDouble (int fieldNum) const
{
	const char *field=GetField(fieldNum);
	char *endptr;
	double result=strtod(field,&endptr);
	if (*endptr!=0) {
		throw SimpleStringException("Double field had some non-numeric content, field text='%s', %s",
			field,GetAdditionalInformationForException().c_str());
	}
	return result;
}
bool CommaSepSeparator::FieldsContainValue (const char *value) const
{
	int f;
	for (f=0; f<GetNumFields(); f++) {
		if (strcmp(GetField(f),value)==0) {
			return true;
		}
	}
	return false;
}


////////////////////////////
// CommaSepAbstractFile

CommaSepAbstractFile::CommaSepAbstractFile (char delimiterChar)
: CommaSepSeparator(delimiterChar)
{
}
CommaSepAbstractFile::~CommaSepAbstractFile ()
{
}


////////////////////
// CommaSepFileReader

CommaSepFileReader::CommaSepFileReader (const char *fileName,char _delimiterChar)
: CommaSepAbstractFile(_delimiterChar)
{
	inFile=ThrowingFopen(fileName,"rt");
	deleteFileOnDestructor=true;
	lineNum=0;

	currLine.resize(128);
}
CommaSepFileReader::CommaSepFileReader (FILE *_inFile,char _delimiterChar,int currLineNum)
: CommaSepAbstractFile(_delimiterChar)
{
	inFile=_inFile;
	deleteFileOnDestructor=false;
	lineNum=currLineNum;

	currLine.resize(128);
}
CommaSepFileReader::~CommaSepFileReader ()
{
	if (deleteFileOnDestructor && inFile!=NULL) {
		fclose(inFile);
	}
}
int CommaSepFileReader::GetLineNum (void) const
{
	return lineNum;
}
bool CommaSepFileReader::ReadLine (void)
{
	lineNum++;

	// read whole line
	int bufferPos=0;
	while (1) {
		if (fgets(&(currLine[bufferPos]),currLine.size()-bufferPos,inFile)==NULL) {
			if (feof(inFile)) {
				return false;
			}
			throw ANSICLibException("Couldn't read next line","fgets");
		}
		const char *s=&*currLine.begin();
		if (strchr(s,'\n')!=NULL) {
			// already read entire line
			break;
		}
		else {
			// line was too long - extend buffer & try again
			bufferPos=strlen(s);
			currLine.resize(currLine.size()*2);
		}
	}

	const char *s=&*currLine.begin();
	currLine[strcspn(s,"\r\n")]=0;

	SeparateLine(currLine);
	return true;
}
std::string CommaSepFileReader::GetAdditionalInformationForException (void) const
{
	char buf[256];
	sprintf(buf,"line #%d",lineNum);
	return std::string(buf);
}
const CommaSepSeparator& CommaSepFileReader::GetCommaSepSeparator (void) const
{
	return (const CommaSepSeparator&) (*this);
}

///////////////////
// CommaSepMetaSep

CommaSepMetaSep::CommaSepMetaSep (const char *fullString,char lineDelimiterChar,char fieldDelimiterChar)
: lines(lineDelimiterChar)
, CommaSepAbstractFile(fieldDelimiterChar)
{
	lines.SeparateLine(fullString);
	lineNum=0;
}
CommaSepMetaSep::~CommaSepMetaSep ()
{
}
bool CommaSepMetaSep::ReadLine (void)
{
	if (lineNum==lines.GetNumFields()) {
		return false;
	}
	SeparateLine(lines.GetField(lineNum));
	lineNum++;
	return true;
}
int CommaSepMetaSep::GetLineNum (void) const
{
	return lineNum;
}

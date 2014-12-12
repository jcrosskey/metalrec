/*
 * Reads.h
 *
 * Created on: Fri Nov  7 11:51:37 EST 2014
 * Author: JJ Chai
 */

#include "Read.h"


/**********************************************************************************************************************
	Default constructor
**********************************************************************************************************************/
Read::Read(void)
{
	// Initialize the variables.
	ID = 0;
	frequency = 0;
	superReadID = 0;
	numInEdges = 0;
	numOutEdges = 0;

	listOfEdgesForward = new vector<Edge *>;
	listOfEdgesForward->resize(listOfEdgesForward->size());			// Resize to 0 to reduce space.
	locationOnEdgesForward = new vector<UINT64>;
	locationOnEdgesForward->resize(locationOnEdgesForward->size());	// Resize to 0 to reduce space.

	containedReadIDs = new vector<UINT64>;
	containedReadIDs->resize(containedReadIDs->size());	// Resize to 0 to reduce space.

	overlapReadIDs = new vector<UINT64>;
	overlapReadIDs->resize(overlapReadIDs->size());	// Resize to 0 to reduce space.

	overlapReadOffsets = new vector<UINT64>;
	overlapReadOffsets->resize(overlapReadOffsets->size());	// Resize to 0 to reduce space.
}

/**********************************************************************************************************************
	Another constructor, from sam align record in string format
**********************************************************************************************************************/
Read::Read(const string & s)
{
	// Initialize the variables.
	ID = 0;
	frequency = 0;
	superReadID = 0;
	numInEdges = 0;
	numOutEdges = 0;

	listOfEdgesForward = new vector<Edge *>;
	listOfEdgesForward->resize(listOfEdgesForward->size());			// Resize to 0 to reduce space.
	locationOnEdgesForward = new vector<UINT64>;
	locationOnEdgesForward->resize(locationOnEdgesForward->size());	// Resize to 0 to reduce space.

	containedReadIDs = new vector<UINT64>;
	containedReadIDs->resize(containedReadIDs->size());	// Resize to 0 to reduce space.

	overlapReadIDs = new vector<UINT64>;
	overlapReadIDs->resize(overlapReadIDs->size());	// Resize to 0 to reduce space.

	overlapReadOffsets = new vector<UINT64>;
	overlapReadOffsets->resize(overlapReadOffsets->size());	// Resize to 0 to reduce space.

	setRead(s);
}

/**********************************************************************************************************************
	Another constructor, from sam align record in BamAlignmentRecord format
**********************************************************************************************************************/
Read::Read(const seqan::BamAlignmentRecord & record)
{
	// Initialize the variables.
	ID = 0;
	frequency = 0;
	superReadID = 0;
	numInEdges = 0;
	numOutEdges = 0;

	listOfEdgesForward = new vector<Edge *>;
	listOfEdgesForward->resize(listOfEdgesForward->size());			// Resize to 0 to reduce space.
	locationOnEdgesForward = new vector<UINT64>;
	locationOnEdgesForward->resize(locationOnEdgesForward->size());	// Resize to 0 to reduce space.

	containedReadIDs = new vector<UINT64>;
	containedReadIDs->resize(containedReadIDs->size());	// Resize to 0 to reduce space.

	overlapReadIDs = new vector<UINT64>;
	overlapReadIDs->resize(overlapReadIDs->size());	// Resize to 0 to reduce space.

	overlapReadOffsets = new vector<UINT64>;
	overlapReadOffsets->resize(overlapReadOffsets->size());	// Resize to 0 to reduce space.

	setRead(record);
}

/**********************************************************************************************************************
	Default destructor
**********************************************************************************************************************/
Read::~Read(void)
{
	// delete all the pointers.
	delete listOfEdgesForward;
	delete locationOnEdgesForward;
	delete containedReadIDs;
	delete overlapReadIDs;
	delete overlapReadOffsets;
}

/**********************************************************************************************************************
	Function to store the read from an alignment record in sam file
**********************************************************************************************************************/
bool Read::setRead(const string & s)
{
	//FILE_LOG(logDEBUG3) << "Use generic sam record parsing\n";
	setFrequency(1);	// Set the frequency to 1.
	stringstream ss(s);
	int flag;
	if (ss.good())
	{
		vector<string> fields = Utils::StringToVector(s, '\t');
		for (size_t i = 0 ; i < fields.size(); i++)
			FILE_LOG(logDEBUG4) << "fileds[" << i << "]: " << fields[i]; 
		readName = seqan::CharString(fields.at(0));	// 1. Query/Read name
		flag = Utils::stringToInt(fields.at(1));	// 2. bitwise Flag
		// 3. Reference Name
		rStart = Utils::stringToInt(fields.at(3)) - 1;	// 4. 1-based leftmost mapping position, changes to 0-based 
		mapQV = Utils::stringToInt(fields.at(4));	// 5. Mapping quality value
		isReverseComplement = ( (flag & 0x10) == 0x10);	// Is SEQ reverse complemented
		cigarString = fields.at(5);	// 6. Cigar string
		setClip(cigarString); // Set the clipping lengths at both ends
		// 7. Ref name of the mate/next read
		// 8. Position of the mate/next read
		// 9. observed template length
		readDnaString = seqan::DnaString(fields.at(9));	// 10. Read string
		FILE_LOG(logDEBUG4) << "read String: " << readDnaString;
		// more fields? XL, XQ, AS, etc
		alignScore = getTag("AS", s);	// 13. Alignment Score generated by BLASR, negative value (smaller score means better alignment)
		FILE_LOG(logDEBUG4) << "AS: " << alignScore;
	}
	readReverseDnaString = readDnaString;	// Initialize the reverse complement string to be the same as forward string
	seqan::reverseComplement(readReverseDnaString);	// Reverse and complement in place
	startCoord = rStart - leftClip;
	return true;
}

/**********************************************************************************************************************
	Function to store the read from an alignment record in sam file, using seqAn BamStream class
**********************************************************************************************************************/
bool Read::setRead(const seqan::BamAlignmentRecord & record) 		// Set the read from alignment record string, (with seqAn BamStream class)
{
	setFrequency(1);	// Set the frequency to 1.
	readName = record.qName;	// query/read name
	readDnaString = record.seq;	// fragment sequence
	readReverseDnaString = readDnaString;	// Initialize the reverse complement string to be the same as forward string
	seqan::reverseComplement(readReverseDnaString);	// Reverse and complement in place
	rStart = record.beginPos;	// Begin position of the alignment on the PacBio read, 0-based (instead of 1 as in the sam file)
	mapQV = record.mapQ;	// mapping quality, 255 for * (not available)
	isReverseComplement = seqan::hasFlagRC(record);	// Is SEQ reverse complemented
	seqan::CharString tagString = record.tags;	// CharString, raw tag string
	seqan::BamTagsDict tagsDict(tagString);	// BamTagsDict, wrapper around the tags member
	seqan::CharString AS_tag = "AS";
	unsigned int AS_pos = 0;	// id of the found tag
	bool AS_Found = seqan::findTagKey(AS_pos, tagsDict, AS_tag);	// See if AS tag exists in the record
	if (!AS_Found)	// If the tag does not exist, set it to 0 (negative is better score)
	{
		FILE_LOG(logWARNING) << "Warning: AS key does not exist in alignment record";
		alignScore = 0;
	}
	else
		seqan::extractTagValue(alignScore, tagsDict, AS_pos);

	seqan::String<seqan::CigarElement<> > cigar_String = record.cigar;	// cigar string of type String<CigarElement<> >
	// cigar string as string, instead of a seqan class, may not be useful at all later.
	stringstream ss;
	for ( size_t i = 0; i < seqan::length(cigar_String); i++)
	{
		ss << cigar_String[i].count;
		ss << cigar_String[i].operation;
	}
	cigarString = ss.str();
	leftClip = 0;	// Set them to 0 first
	rightClip = 0;
	if ( seqan::length(cigar_String) > 0 )	// If cigar_String is not empty, check whether there are clippings at the ends
	{
		if (cigar_String[0].operation == 'S')
			leftClip = cigar_String[0].count;

		if (cigar_String[(int) length(cigar_String) - 1].operation == 'S')
			rightClip = cigar_String[(int) length(cigar_String) - 1].count;
	}
	startCoord = rStart - leftClip;
	return true;
}

/**********************************************************************************************************************
  Get the clipped length on both ends from cigar string.
Note: rightClip is not used currently.
**********************************************************************************************************************/
bool Read::setClip(const string & cigarString)
{
	size_t pos = 0;
	size_t S_pos = cigarString.find('S');
	if ( S_pos != string::npos )
	{
		leftClip = Utils::stringToUnsignedInt( cigarString.substr(pos, S_pos) );
		pos = S_pos + 1;
		S_pos = cigarString.find('S', pos);
		if ( S_pos != string::npos )
		{
			size_t last_S_pos = S_pos;
			while( isdigit ( cigarString.at(S_pos - 1)) ) // if the previous digit is a number
			{
				S_pos = S_pos - 1 ;
			}
			rightClip = Utils::stringToUnsignedInt( cigarString.substr(S_pos, last_S_pos) );

		}
		else
			rightClip = 0;
	}
	else
	{
		leftClip = 0;
		rightClip = 0;
	}
	return true;
}

/**********************************************************************************************************************
	This function assigns an ID to the read.
**********************************************************************************************************************/
bool Read::setReadID(UINT64 id)
{
	ID = id;	// Set the read number.
	return true;
}

/**********************************************************************************************************************
	This function sets the frequency of the read.
**********************************************************************************************************************/
bool Read::setFrequency(UINT32 freq)
{
	if(freq < 1) 
		MYEXIT("Frequency less than 1.");
	frequency = freq;	// Set the frequency of the read.
	return true;
}

/**********************************************************************************************************************
	Return the ending coordinates of the read, disregarding the alignment to PacBio read, 
	assume that all the mismatches are substitutions in Illumina reads.
**********************************************************************************************************************/
INT32 Read::getEndCoord(void)
{
	return ( startCoord + length(readDnaString) ); // 0-based, position after the last bp of the read
}

/**********************************************************************************************************************
	Return the value for a certain required tag from the align record
**********************************************************************************************************************/
INT32 Read::getTag(const string & tagName, const string & alignRecord)
{
	size_t tag_position = alignRecord.find(tagName + ":i:"); // Look for the string in the alignment record
	if (tag_position != string::npos)	// The tag is found 
	{
		size_t val_start_position = tag_position + tagName.size() + 3;	// Shift to the right, past the tagname, and the :i: (format indicator)
		size_t val_end_position = alignRecord.find('\t', val_start_position);	// Find the end of the field, which is one char to the right of the value
		string val;	// Get the value
		val.assign(alignRecord.begin() + val_start_position , alignRecord.begin() + val_end_position);	// Save the value for the tag
		return Utils::stringToInt(val);
	}
	else
	{
		FILE_LOG(logWARNING) << "Tag " << tagName << " does not exist in the record " << alignRecord;
		return -1;	// if the tag is not found, set the value to -1 ( another value might be better)
	}

}

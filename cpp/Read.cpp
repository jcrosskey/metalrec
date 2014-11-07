/*
 * Reads.h
 *
 * Created on: Fri Nov  7 11:51:37 EST 2014
 * Author: JJ Chai
 */

#include "Common.h"
#include "Read.h"


/**********************************************************************************************************************
	Default constructor
**********************************************************************************************************************/
Read::Read(void)
{
	// Initialize the variables.
	ID = 0;
	frequency = 0;
	isContainedRead = false;
	superReadID = 0;

	listOfEdgesForward = new vector<Edge *>;
	listOfEdgesForward->resize(listOfEdgesForward->size());			// Resize to 0 to reduce space.
	listOfEdgesReverse = new vector<Edge *>;
	listOfEdgesReverse->resize(listOfEdgesReverse->size());			// Resize to 0 to reduce space.
	locationOnEdgesForward = new vector<UINT64>;
	locationOnEdgesForward->resize(locationOnEdgeForward->size());	// Resize to 0 to reduce space.
	locationOnEdgesReverse = new vector<UINT64>;
	locationOnEdgesReverse->resize(locationOnEdgeReverse->size());	// Resize to 0 to reduce space.

}

/**********************************************************************************************************************
	Another constructor
**********************************************************************************************************************/
Read::Read(const string & s)
{
	// Initialize the variables.
	ID = 0;
	frequency = 0;
	isContainedRead = false;
	superReadID = 0;

	listOfEdgesForward = new vector<Edge *>;
	listOfEdgesForward->resize(listOfEdgesForward->size());			// Resize to 0 to reduce space.
	listOfEdgesReverse = new vector<Edge *>;
	listOfEdgesReverse->resize(listOfEdgesReverse->size());			// Resize to 0 to reduce space.
	locationOnEdgesForward = new vector<UINT64>;
	locationOnEdgesForward->resize(locationOnEdgeForward->size());	// Resize to 0 to reduce space.
	locationOnEdgesReverse = new vector<UINT64>;
	locationOnEdgesReverse->resize(locationOnEdgeReverse->size());	// Resize to 0 to reduce space.
	setRead(s);
}

/**********************************************************************************************************************
	Default destructor
**********************************************************************************************************************/
Read::~Read(void)
{
	// delete all the pointers.
	delete listOfEdgesForward;
	delete listOfEdgesReverse;
	delete locationOnEdgesForward;
	delete locationOnEdgesReverse;
}

/**********************************************************************************************************************
	Function to store the read from an alignment record in sam file
**********************************************************************************************************************/
bool Read::setRead(const string & s)
{
	setFrequency(1);	// Set the frequency to 1.
	stringstream ss(s);
	int flag;
	if (ss.good())
	{
		vector<string> fields = Utils::StringToVector(s, '\t')
		readName = fields.at(0);	// 1. Query/Read name
		flag = Utils::stringToInt(fields.at(1));	// 2. bitwise Flag
		// 3. Reference Name
		rStart = Utils::stringToInt(fields.at(3)) - 1;	// 4. 1-based leftmost mapping position 
		mapQV = Utils::stringToInt(fields.at(4));	// 5. Mapping quality value
		cigarString = fields.at(5);	// 6. Cigar string
		setClip(cigarString); // Set the clipping lengths at both ends
		// 7. Ref name of the mate/next read
		// 8. Position of the mate/next read
		readString = fields.at(8);	// 9. Read string
		// more fields? XL, XQ, AS, etc
	}
	readReverseString = reverseComplement(readString);
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
	if ( S_pos != npos )
	{
		leftClip = Utils::stringToUnsignedInt( cigarString.substr(pos, S_pos) );
		pos = S_pos + 1;
		S_pos = cigarString.find('S', pos);
		if ( S_pos != npos )
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
bool Read::setReadNumber(UINT64 id)
{
	if(id <= 0) 
		MYEXIT("ID less than 1.");
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
	Returns the reverse complement of a read.
**********************************************************************************************************************/
string Read::reverseComplement(const string & read)
{
	UINT64 stringLength = read.length();
	string reverse(stringLength,'0');
	for(UINT64 i = 0;i < stringLength; i++)						// Then complement the string. Change A to T, C to G, G to C and T to A.
	{
		if( read[i] & 0X02 ) // C or G
			reverse.at(stringLength -  i - 1 ) = read[i] ^ 0X04;
		else // A <==> T
			reverse.at(stringLength -  i - 1 ) = read[i] ^ 0X15;
	}
	return reverse; // return the reverse complement as a string
}

/**********************************************************************************************************************
	Return the starting and ending coordinates of the read, disregarding the alignment to PacBio read, 
	assume that all the mismatches are substitutions in Illumina reads.
**********************************************************************************************************************/
INT32 Read::getStartCoord(void)
{
	return (rStart - leftClip);	// 0-based, position of the first bp
}

INT32 Read::getEndCoord(void)
{
	return ( rStart - leftClip + readString.length() ); // 0-based, position after the last bp of the read
}

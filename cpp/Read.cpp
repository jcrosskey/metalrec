/*
 * Reads.cpp
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
	readName = "";
	readDnaString = "";
	cigarString = "";
	rStart = 0;
	leftClip = 0;
	rightClip = 0;
	numOfInsertions = 0;
	flag = 0;
	numOfDeletions = 0;
	numOfMatches = 0;
	numOfSubstitutions = 0;
	editDistance = 0;
	startCoord = 0;
	refName = "";
	superReadID = 0;
	numInEdges = 0;
	numOutEdges = 0;

	listOfEdgesForward = new vector<Edge *>;
	listOfEdgesForward->resize(listOfEdgesForward->size());			// Resize to 0 to reduce space.
	locationOnEdgesForward = new vector<UINT64>;
	locationOnEdgesForward->resize(locationOnEdgesForward->size());	// Resize to 0 to reduce space.

	containedReadIDs = new vector<UINT64>;
	containedReadIDs->resize(containedReadIDs->size());	// Resize to 0 to reduce space.

	overlapReadOffsets = new vector<UINT64>;
	overlapReadOffsets->resize(overlapReadOffsets->size());	// Resize to 0 to reduce space.
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  Read
 *  Description:  Copy constructor, make deep copy of known member, instead of default
 *  	          shallow copy
 * =====================================================================================
 */
Read::Read(const Read & R)
{
	ID = R.ID;
	readName = R.readName;
	readDnaString = R.readDnaString;
	cigarString = R.cigarString;
	rStart = R.rStart;
	leftClip = R.leftClip;
	rightClip = R.rightClip;
	numOfInsertions = R.numOfInsertions;
	flag = R.flag;
	numOfDeletions = R.numOfDeletions;
	numOfMatches = R.numOfMatches;
	numOfSubstitutions = R.numOfSubstitutions;
	editDistance = R.editDistance;
	startCoord = R.startCoord;
	refName = R.refName;
	superReadID = R.superReadID;
	numInEdges = R.numInEdges;
	numOutEdges = R.numOutEdges;
	size_t k = 0;

	listOfEdgesForward = new vector<Edge *>;
	listOfEdgesForward->resize((R.listOfEdgesForward)->size());
	for(k = 0; k < listOfEdgesForward->size(); k++){
		Edge * EdgeForward = new Edge();
		*EdgeForward = *((R.listOfEdgesForward)->at(k));
		listOfEdgesForward->at(k) = EdgeForward;
	}

	locationOnEdgesForward = new vector<UINT64>;
	locationOnEdgesForward->resize(R.locationOnEdgesForward->size());
	for(k = 0; k < locationOnEdgesForward->size(); k++){
		locationOnEdgesForward->at(k) = R.locationOnEdgesForward->at(k);
	}

	containedReadIDs = new vector<UINT64>;
	containedReadIDs->resize((R.containedReadIDs)->size());
	for(k = 0; k < containedReadIDs->size(); k++){
		containedReadIDs->at(k) = (R.containedReadIDs)->at(k);
	}

	overlapReadOffsets = new vector<UINT64>;
	overlapReadOffsets->resize((R.overlapReadOffsets)->size());
	for(k = 0; k < overlapReadOffsets->size(); k++){
		overlapReadOffsets->at(k) = (R.overlapReadOffsets)->at(k);
	}

}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  operator=
 *  Description:  Copy assignment, make deep copy of known member, instead of default
 *  	          shallow copy
 * =====================================================================================
 */
Read & Read::operator= (const Read & R)
{
	ID = R.ID;
	readName = R.readName;
	readDnaString = R.readDnaString;
	cigarString = R.cigarString;
	rStart = R.rStart;
	leftClip = R.leftClip;
	rightClip = R.rightClip;
	//mapQV = R.mapQV;
	numOfInsertions = R.numOfInsertions;
	flag = R.flag;
	numOfDeletions = R.numOfDeletions;
	numOfMatches = R.numOfMatches;
	numOfSubstitutions = R.numOfSubstitutions;
	editDistance = R.editDistance;
	//alignScore = R.alignScore;
	startCoord = R.startCoord;
	refName = R.refName;
	superReadID = R.superReadID;
	numInEdges = R.numInEdges;
	numOutEdges = R.numOutEdges;
	size_t k = 0;

	for(k = 0; k < listOfEdgesForward->size(); k++){
		delete listOfEdgesForward->at(k);
	}
	delete listOfEdgesForward;
	listOfEdgesForward = new vector<Edge *>;
	size_t numEdgesForward = (R.listOfEdgesForward)->size();
	listOfEdgesForward->resize(numEdgesForward);
	for(k = 0; k < numEdgesForward; k++){
		Edge * EdgeForward = new Edge();
		*EdgeForward = *((R.listOfEdgesForward)->at(k));
		listOfEdgesForward->at(k) = EdgeForward;
	}

	delete locationOnEdgesForward;
	locationOnEdgesForward = new vector<UINT64>;
	locationOnEdgesForward->resize(numEdgesForward);
	for(k = 0; k < numEdgesForward; k++){
		locationOnEdgesForward->at(k) = R.locationOnEdgesForward->at(k);
	}

	delete containedReadIDs;
	containedReadIDs = new vector<UINT64>;
	containedReadIDs->resize((R.containedReadIDs)->size());
	for(k = 0; k < containedReadIDs->size(); k++){
		containedReadIDs->at(k) = (R.containedReadIDs)->at(k);
	}


	delete overlapReadOffsets;
	overlapReadOffsets = new vector<UINT64>;
	overlapReadOffsets->resize((R.overlapReadOffsets)->size());
	for(k = 0; k < overlapReadOffsets->size(); k++){
		overlapReadOffsets->at(k) = (R.overlapReadOffsets)->at(k);
	}
	return *this;
}


/**********************************************************************************************************************
	Another constructor, from sam align record in string format
**********************************************************************************************************************/
Read::Read(const string & s)
{
	// Initialize the variables.
	ID = 0;
	readName = "";
	readDnaString = "";
	cigarString = "";
	rStart = 0;
	leftClip = 0;
	rightClip = 0;
	//mapQV = 0;
	numOfInsertions = 0;
	flag = 0;
	numOfDeletions = 0;
	numOfMatches = 0;
	numOfSubstitutions = 0;
	editDistance = 0;
//	alignScore = 0;
	startCoord = 0;
	refName = "";
	superReadID = 0;
	numInEdges = 0;
	numOutEdges = 0;

	listOfEdgesForward = new vector<Edge *>;
	listOfEdgesForward->resize(listOfEdgesForward->size());			// Resize to 0 to reduce space.

	locationOnEdgesForward = new vector<UINT64>;
	locationOnEdgesForward->resize(locationOnEdgesForward->size());	// Resize to 0 to reduce space.

	containedReadIDs = new vector<UINT64>;
	containedReadIDs->resize(containedReadIDs->size());	// Resize to 0 to reduce space.

	overlapReadOffsets = new vector<UINT64>;
	overlapReadOffsets->resize(overlapReadOffsets->size());	// Resize to 0 to reduce space.

	setRead(s);
}


/**********************************************************************************************************************
	Default destructor
**********************************************************************************************************************/
Read::~Read(void)
{
	delete listOfEdgesForward;
	delete locationOnEdgesForward;
	delete containedReadIDs;
	delete overlapReadOffsets;
}

/**********************************************************************************************************************
	Function to store the read from an alignment record in sam file
**********************************************************************************************************************/
bool Read::setRead(const string & s)
{
	//FILE_LOG(logDEBUG3) << "Use generic sam record parsing\n";
	stringstream ss(s);
	if (ss.good())
	{
		vector<string> fields = Utils::StringToVector(s, '\t');
		if (fields.size() < 11)
			FILE_LOG(logERROR) << "Align record does not have 11 mandatory fields, quit...";
		else
		{
			refName = fields.at(2);
			for (size_t i = 0 ; i < fields.size(); i++)
				FILE_LOG(logDEBUG4) << "fileds[" << i << "]: " << fields[i]; 
			readName = fields.at(0);	// 1. Query/Read name
			flag = stoi(fields.at(1));	// 2. bitwise Flag
			if(!isReadUnmapped())     /* If read is mapped, set other numbers */
			{
				// 3. Reference Name
				rStart = stoi(fields.at(3)) - 1;	// 4. 1-based leftmost mapping position, changes to 0-based 
				//mapQV = stoi(fields.at(4));	// 5. Mapping quality value
				//isReverseComplement = ( (flag & 0x10) == 0x10);	// Is SEQ reverse complemented
				cigarString = fields.at(5);	// 6. Cigar string
				//setClip(cigarString); // Set the clipping lengths at both ends
				setEdits(cigarString);
				// 7. Ref name of the mate/next read
				// 8. Position of the mate/next read
				// 9. observed template length
				readDnaString = fields.at(9);	// 10. Read string
				FILE_LOG(logDEBUG4) << "read String: " << readDnaString;
				// more fields? XL, XQ, AS, etc
				//alignScore = getTag("AS", s);	// 13. Alignment Score generated by BLASR, negative value (smaller score means better alignment)
				if(numOfMatches == 0 && numOfSubstitutions == 0)
				{
					editDistance = getTag("NM", s); // edit distance from the PacBio read (NM)
					if (editDistance == -1)
						flag = flag | 0x4;
					else if (editDistance >= (int)(numOfInsertions + numOfDeletions))
						numOfSubstitutions = editDistance - numOfInsertions - numOfDeletions;
					else
						FILE_LOG(logWARNING) << "In read " << readName << " editDistance " << editDistance << " is smaller than the number of indels " << numOfDeletions+numOfInsertions << " , NM is not reliable and number of substitutions is set to 0.";
				}
				//FILE_LOG(logDEBUG4) << "AS: " << alignScore << "; NM: " << editDistance;
			}
		}
		startCoord = rStart - leftClip;
	}
	return true;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  setEdits
 *  Description:  Set the number of indels and substitutions of the alignment, from cigar string
 * =====================================================================================
 */
bool Read::setEdits(const string & s)
{
//	FILE_LOG(logINFO) << s;
	size_t pos, pos1 = 0, length = 0;
	unsigned int numOfMatchMismatches = 0 ;
	int num;
	bool firstOp = true;                    /* if the operation is the first one */
	for(pos = 0; pos < s.length(); pos++)
	{
		if (isdigit(s.at(pos)))         /* if this position is digit */
		{
			length++;
		}
		else    /* if this position is alphabetic, meaning a new operation is found */
		{
			num = stoi(s.substr(pos1, length));
			switch(s.at(pos))
			{
				case 'D' :
					numOfDeletions += num;
					break;
				case 'I':
					numOfInsertions += num;
					break;
				case 'M':
					numOfMatchMismatches += num;
					break;
				case '=':
					numOfMatches += num;
					break;
				case 'X':
					numOfSubstitutions += num;
					break;
				case 'S':
					if(firstOp)
					{
						leftClip = num;
//						FILE_LOG(logINFO) << readName << " leftClip: " << leftClip;
					}
					else{
						rightClip = num;
//						FILE_LOG(logINFO) << readName << " rightClip: " << rightClip;
					}
				case 'N':
					break;
				case 'H':
					break;
				case 'P':
					break;
				default:
					FILE_LOG(logWARNING) << "Unknown operation string in CIGAR string " << s.at(pos);
			}
			pos1 = pos + 1;         /* update the positions to find the next number */
			length = 0;
			if(firstOp)
				firstOp = false;
		}
	}
	return true;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  isReadGood
 *  Description:  Test and see if read is good (whether it should be mapped to this PacBio read)
 *  		1. clipped length should be less than half of the whole read
 *  		2. whole read should be contained in the LR
 *  		3. sub, ins, del errors should be within the thresholds
 * =====================================================================================
 */
bool Read::isReadGood(const float & indelRate, const float & insRate, const float & delRate, const float & subRate, const UINT16 & PacBioReadLength, const float & percentInLR)
{
	if (isReadUnmapped())
		return false;
	UINT64 readLength = readDnaString.length();
	UINT64 mappedLength = readDnaString.length() - leftClip - rightClip; /* length of the aligned segment, not counting the left and right clipping */
	FILE_LOG(logDEBUG3) << numOfSubstitutions << " substitutions in " << mappedLength << " mapped bps";
	FILE_LOG(logDEBUG3) << numOfDeletions << " deletions and " << numOfInsertions << " insertion in " << mappedLength << " mapped bps";
	FILE_LOG(logDEBUG3) << "leftClip: " << leftClip << " rightClip: " << rightClip << " in " << readLength << " bps";
	/* Too many clipped bases TODO: need to put these in the config file */
	if (float(leftClip + rightClip) > readLength * 0.5 || leftClip > 100 || rightClip > 100) 
		return false;
	/* If the read is not totally included in the PacBio read */
	if (startCoord < 0 || getEndCoord() > PacBioReadLength)               
	{
		UINT64 basesOutOfLR = (startCoord > 0 ? 0 : -startCoord);
		basesOutOfLR += ((getEndCoord() > PacBioReadLength)?(getEndCoord()-PacBioReadLength):0);
		if (basesOutOfLR > readLength * (1-percentInLR))
			return false;
	}
	if (numOfSubstitutions >  mappedLength * subRate) /* too many substitution errors */
		return false;
	if ( numOfInsertions  > mappedLength * insRate) /* too many insertion erros */
		return false;
	if ( numOfDeletions  > mappedLength * delRate) /* too many deletion erros */
		return false;
	if ( (numOfDeletions + numOfInsertions) > mappedLength * indelRate) /* too many indel erros */
		return false;

	return true;
}


/**********************************************************************************************************************
	Return the ending coordinates of the read, disregarding the alignment to PacBio read, 
	assume that all the mismatches are substitutions in Illumina reads.
**********************************************************************************************************************/
INT32 Read::getEndCoord(void)
{
	return ( startCoord + readDnaString.length()) ; // 0-based, position after the last bp of the read
}

INT32 Read::getEndCoordInLR(void)
{
	return ( startCoord + readDnaString.length() - rightClip) ; // 0-based, position after the last bp of the read, in the long read
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
		return stoi(val);
	}
	else
	{
		FILE_LOG(logWARNING) << "Tag " << tagName << " does not exist in the record " << alignRecord;
		return -1;	// if the tag is not found, set the value to -1 ( another value might be better)
	}

}

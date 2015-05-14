/*
 * Reads.h
 *
 * Created on: Fri Nov  7 11:51:37 EST 2014
 * Author: JJ Chai
 */


#ifndef READ_H_
#define READ_H_

#include "Common.h"
#include "Edge.h"

class Edge;

// TODO: Add paired end reads support, add mate as a variable for read so that later paired end reads can be used.

/**********************************************************************************************************************
	Class to store a read.
**********************************************************************************************************************/
class Read
{
	private:
		UINT64 ID;	// Unique Identification of the read.
		string readName;	// Query/Read name (original name in the input file)
		int flag;
		string readDnaString;	// Same as below, but with seqAn's String<Dna> implementation
		string cigarString;	// Cigar string, might be useful in the future. 
		UINT32 rStart;	// Alignment start position on the PacBio read
		UINT32 leftClip;	// Clipped length at the left end
		UINT32 rightClip;	// Clipped length at the right end
		UINT64 numOfInsertions;         /* number of insertions into the PacBio read */
		UINT64 numOfDeletions;          /* number of deletions from the PacBio read */
		UINT64 numOfMatches;
		UINT64 numOfSubstitutions;        /* number of substitutions */
		int editDistance;            /* edit distance of the read from the PacBio read */
		INT32 startCoord;	// Coordinate of the first bp in the read
		string refName;                 /* reference sequence name */

		/*  pointer members, need to be careful when destroy and copy object */
		vector<Edge *> *listOfEdgesForward;   	// List of edges that contain the forward string of this read.
		vector<UINT64> *locationOnEdgesForward;	// List of locations on the edges that contain the forward string of the current read.
		vector<UINT64> *containedReadIDs;	// Pointer to vector of IDs of reads contained in this read
		vector<UINT64> *overlapReadOffsets;	// Pointer to vector of offsets for the corresponding overlapping reads

		bool setEdits(const string & s);
		INT32 getTag(const string & tagName, const string & alignRecord);

	public:
		Read(void);	// Default constructor.
		Read(const string & s);	// Another constructor, from alignment record in string format.
		~Read(void);	// Destructor.
		Read(const Read & R); // copy constructor
		Read & operator= (const Read & R); // Copy assignment
		/*  TODO: should we include move constructor and assignment functions? */

		UINT64 superReadID;	// ID of the (longest) read containing this read. 
					// 0 = not a contained read, otherwise superReadID contains the ID of the uniqe super read.
		UINT16 numInEdges;	// Number of edges going into the read
		UINT16 numOutEdges;	// Number of edges going out of the read

		/* mutators */
		bool setRead(const string & s); 		// Set the read from alignment record string, (with generic string parsing)
		void setReadID(UINT64 id){ID = id;} 			// Set the read ID.
		void setStartCoord(INT32 start_coord){startCoord = start_coord;}	// Set the starting coordinate

		/* accessors */
		bool isReadUnmapped(void) {return ((flag & 0x4) == 0x4);}        /* from tag determine if the read is mapped */
		UINT64 getNumOfSubstitutionsInRead(void){return numOfSubstitutions;}
		UINT64 getNumOfInsertionsInRead(void){return numOfInsertions;}
		UINT64 getNumOfDeletionsInRead(void){return numOfDeletions;}
		bool isReadGood(const float & indelRate, const float & insRate, const float & delRate, const float & subRate, const UINT16 & PacBioReadLength, const float & percentInLR=0.80); /* Test if read is good enough */
		bool isContainedRead(void){return (superReadID == 0 ? false : true);};	// Whether the read is contained in another read
		string getDnaStringForward(void){return readDnaString;}	// Get the forward string of the current read.
		string getReadName(void){return readName;}	// Get the name of the current read.
		INT32 getStartCoord(void){return startCoord;}	// Get the starting coordinate of the read, with clipped part added back
		UINT32 getLeftClip(void) {return leftClip;}
		UINT32 getRightClip(void) {return rightClip;}
		INT32 getEndCoord(void);	// Get the ending coordinate of the read, with the clipped part considered
		INT32 getEndCoordInLR(void);    /* Get the ending coordinate of the read with clipped part removed */
		UINT32 getrStart(void) {return rStart;}	// Get the leftmost alignment position on the reference sequence.
		size_t getReadLength(void){return readDnaString.length();}	// Get the length of the string in the current read.
		size_t getAlignedLength(void){return (readDnaString.length()-leftClip-rightClip);}	// Get the length of the string in the current read.
		UINT64 getClippedLength(void){return (leftClip+rightClip);}	// Get the length of the string in the current read.
		UINT64 getID(void) {return ID;}	// Get the read number of the current read.
		bool isReverseComplemented(void) {return  (flag & 0x10) == 0x10;}	// Get the reverse complement status
		string getRefName(void){return refName;} /* Get reference sequence name */

		vector<Edge *> * getListOfEdgesForward(void){return listOfEdgesForward;}	// Get the list of edges that contain the forward string of the current read.
		vector<UINT64> * getLocationOnEdgesForward(void){return locationOnEdgesForward;}	// Get the list of locations on the edges that contain the forward string of the current read.
		vector<UINT64> * getContainedReadIDs(void){return containedReadIDs;}	// Get list of IDs of reads contained in this read.
		vector<UINT64> * getOverlapReadOffsets(void){return overlapReadOffsets;}	// Get list of offsets for the corresponding overlapping reads.
};

#endif /* READS_H_ */

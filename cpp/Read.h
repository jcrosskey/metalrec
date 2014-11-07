/*
 * Reads.h
 *
 * Created on: Fri Nov  7 11:51:37 EST 2014
 * Author: JJ Chai
 */


#ifndef READ_H_
#define READ_H_

#include "Common.h"
#include "Utils.h"

class Edge;

/* No mate pair information will be considered for the current implementation of error correction */

/**********************************************************************************************************************
	Class to store a read.
**********************************************************************************************************************/

class Read
{
	private:
		UINT64 ID;	// Unique Identification of the read.
		string readName;
		string readString;	// Sequence of the read.
		string readReverseString;	// Reverse complement string of the read ( need for bidirectional graph )
		string cigarString;	// Cigar string, might be useful in the future
		UINT32 frequency;	// Frequency of the read. Number of times this read is present in the dataset. Used in some statistical analysis.
		UINT32 rStart;	// Alignment start position on the PacBio read
		UINT32 rEnd;	// Alignment end position on the PacBio read
		UINT32 leftClip;	// Clipped length at the left end
		UINT32 rightClip;	// Clipped length at the right end
		INT32 mapQV;	// Mapping quality value, smaller the better, 255 means quality is not available

		/* Only the forward edges will be considered later on, but the bidirectional graph structure is kept */
		vector<Edge *> *listOfEdgesForward;   	// List of edges that contain the forward string of this read.
		vector<UINT64> *locationOnEdgesForward;	// List of locations on the edges that contain the forward string of the current read.
		vector<Edge *> *listOfEdgesReverse;	// List of edges that contain the reverse string of this read.
		vector<UINT64> *locationOnEdgesReverse;	// List of locations on the edges that contain the reverse string of the current read.

	public:
		Read(void);	// Default constructor.
		Read(const string & s);	// Another constructor.
		~Read(void);	// Destructor.

		bool isContainedRead;	// Whether the read is contained in another read
		UINT64 superReadID;	// ID of the (longest) read containing this read. \
					0 = not a contained read, otherwise superReadID contains the ID of the uniqe super read.

		/* mutators */
		bool setRead(const string & s); 		// Set the read sequence string.
		bool setReadID(UINT64 id); 			// Set the read ID.
		bool setFrequency(UINT32 freq);	// Set the ferquency of the read.

		/* accessors */
		string getStringForward(void){return readString;}	// Get the forward string of the current read.
		string getStringReverse(void){return readReverseString;}	// Get the reverse string of the current read.
		INT32 getStartCoord(void);	// Get the starting coordinate of the read, with clipped part added back
		INT32 getEndCoord(void);	// Get the ending coordinate of the read
		UINT16 getReadLength(void){return readString.length();}	// Get the length of the string in the current read.
		UINT64 getID(void) {return ID;}	// Get the read number of the current read.
		UINT32 getFrequency(void) {return frequency;}	// Get the frequency of the current read.

		vector<Edge *> * getListOfEdgesForward(void){return listOfEdgesForward;}	// Get the list of edges that contain the forward string of the current read.
		vector<UINT64> * getLocationOnEdgesForward(void){return locationOnEdgesForward;}	// Get the list of locations on the edges that contain the forward string of the current read.
		vector<Edge *> * getListOfEdgesReverse(void){return listOfEdgesReverse;}	// Get the list of edges that contain the reverse string of the current read.
		vector<UINT64> * getLocationOnEdgesReverse(void){return locationOnEdgesReverse;}	// Get the list of locations on the edges that contain the reverse string of the current read.

};

#endif /* READS_H_ */

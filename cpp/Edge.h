/*
 * Edge.h
 *
 * Created on: Fri Nov  7 15:13:39 EST 2014
 * Author: JJ Chai
 */

#ifndef EDGE_H_
#define EDGE_H_
#include "Common.h"
#include "Read.h"

/**********************************************************************************************************************
  Class to store an edge.
 **********************************************************************************************************************/
class Edge{
	private:
		Read *source; 	// Source read u
		Read *destination; 	// Destination read v
		UINT64 overlapOffset;	// Length of the overlap.
		// overlap offset in the following example is 6
		// 	  012345678901234567890123456789
		//	u ACTTACGGGATTATACCATCGAGA
		//	v       GGGATTATACCATCGAGATTCAAT
		UINT16 numOfSubstitions;             /* number of substitutions in the edge */
		vector<UINT64> * listOfReads; 	// List of ordered reads in the current edge.
		vector<UINT16> * listOfOverlapOffsets; 	// List of overlap offsets of the ordered reads in the current edge.
		vector<UINT64> * listOfSubstitutionPoses; /* List of substitution positions */

	public:
		bool transitiveRemovalFlag;	// Used to mark transitive edges.
		UINT16 flow;	// Store the flow in the current edge.
		UINT64 coverageDepth;	// Estimated depth of coverage.

		Edge(void);	// Default constructor.
		Edge(Read *from, Read *to, UINT64 length, UINT16 numSub, vector<UINT64> *listSubs); 	// Another constructor.
		Edge(Read *from, Read *to, UINT64 length, UINT16 numSub, vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets, vector<UINT64> *listSubs);
		~Edge();	// Destructor.

		bool makeEdge(Read *from, Read *to, UINT64 length, UINT16 numSub, vector<UINT64> *listSubs);
		bool makeEdge(Read *from, Read *to, UINT64 length,  UINT16 numSub, vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets, vector<UINT64> *listSubs);

		Read * getSourceRead() {return source;}	// Get the read object of the source node.
		Read * getDestinationRead() {return destination; }	// Get the read object of the destination node.
		UINT64 getOverlapOffset() {return overlapOffset;}	// Return the overlap offset.
		UINT16 getNumOfSubstitutions() {return numOfSubstitions;}	// Return the number of substitutions in the edge.

		vector<UINT64> * getListOfReads() {return listOfReads;}	// Get the ordered list of reads in the current edge.
		vector<UINT16> * getListOfOverlapOffsets() {return listOfOverlapOffsets;} // Get the list of ordered offset.
		vector<UINT64> * getListOfSubstitutionPoses() {return listOfSubstitutionPoses;} // Get the list of substitution positions.
};

#endif /* EDGE_H_ */

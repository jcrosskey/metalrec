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

class Read;

/**********************************************************************************************************************
  Class to store an edge.
 **********************************************************************************************************************/
class Edge{
	private:
		Read *source; 	// Source read u
		Read *destination; 	// Destination read v
		UINT64 overlapOffset;	// Length of the overlap.
		UINT64 edgeID;                  /* ID of the edge, used to track different edges, especially the ones connecting same nodes (in bubbles) in the graph */
		// overlap offset in the following example is 6
		// 	  012345678901234567890123456789
		//	u ACTTACGGGATTATACCATCGAGA
		//	v       GGGATTATACCATCGAGATTCAAT
		UINT16 numOfSubstitions;             /* number of substitutions in the edge */
		UINT32 endCoordinateLimit;      /* Limit of the ending coordinate */
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
		Edge(const Edge & E);
		Edge & operator= (const Edge & E);

		bool makeEdge(Read *from, Read *to, UINT64 length, UINT16 numSub, vector<UINT64> *listSubs);
		bool makeEdge(Read *from, Read *to, UINT64 length,  UINT16 numSub, vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets, vector<UINT64> *listSubs);
		bool setEndCorrdinateLimit(UINT64 endCoordLimit){endCoordinateLimit = endCoordLimit; return true;} /* set the end coordinate limit */
		bool operator<(Edge & anotherEdge);
		bool setID(UINT64 edge_id){edgeID=edge_id; return true;} /* set ID for an edge */

		UINT64 getStringLengthInRange(); /* Get the string length in the edge that is also in the range limit */
		Read * getSourceRead() {return source;}	// Get the read object of the source node.
		Read * getDestinationRead() {return destination; }	// Get the read object of the destination node.
		UINT64 getOverlapOffset() {return overlapOffset;}	// Return the overlap offset.
		UINT16 getNumOfSubstitutions() {return numOfSubstitions;}	// Return the number of substitutions in the edge.
		UINT64 getEdgeID(void){return edgeID;} /* Return edgeID */


		vector<UINT64> * getListOfReads() {return listOfReads;}	// Get the ordered list of reads in the current edge.
		vector<UINT16> * getListOfOverlapOffsets() {return listOfOverlapOffsets;} // Get the list of ordered offset.
		vector<UINT64> * getListOfSubstitutionPoses() {return listOfSubstitutionPoses;} // Get the list of substitution positions.
};

#endif /* EDGE_H_ */

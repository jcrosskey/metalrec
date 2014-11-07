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
		vector<UINT64> * listOfReads; 	// List of ordered reads in the current edge.
		vector<UINT16> * listOfOverlapOffsets; 	// List of overlap offsets of the ordered reads in the current edge.
		vector<UINT8> * listOfOrientations;	// List of orientations of the ordered reads in the current edge.
		Edge *reverseEdge;	// This points to the reverse edge in the overlap graph.
		// Each edge (u,v) is present twice in the graph.
		// Edge (u,v) is present in the list of edges of node u
		// Edge (v,u) is present in the list of edges of node v.
		// Edge (u,v) and (v,u) are called twin edge and are reverse of one another.


	public:
		bool transitiveRemovalFlag;	// Used to mark transitive edges.
		UINT16 flow;	// Store the flow in the current edge.
		UINT64 coverageDepth;	// Estimated depth of coverage.
		UINT64 SD;
		Edge(void);	// Default constructor.
		Edge(Read *from, Read *to, UINT64 orient, UINT64 length); 	// Another constructor.
		Edge(Read *from, Read *to, UINT64 orient, UINT64 length, vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets, vector<UINT8> * listOrientations);
		~Edge();	// Destructor.
		bool makeEdge(Read *from, Read *to, UINT64 orient, UINT64 length);
		bool makeEdge(Read *from, Read *to, UINT64 orient, UINT64 length,  vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets, vector<UINT8> * listOrientations);
		string getStringInEdge(void); 	// return the string in the edge
		bool setReverseEdge(Edge * edge);	// Set the pointer to the reverse edge.
		Read * getSourceRead() {return source;}	// Get the read object of the source node.
		Read * getDestinationRead() {return destination; }	// Get the read object of the destination node.
		UINT8 getOrientation() {return overlapOrientation;}	// Return the orientation of the edge.
		UINT64 getOverlapOffset() {return overlapOffset;}	// Return the overlap offset.
		vector<UINT64> * getListOfReads() {return listOfReads;}	// Get the ordered list of reads in the current edge.
		vector<UINT16> * getListOfOverlapOffsets() {return listOfOverlapOffsets;} // Get the list of ordered offset.
		vector<UINT8> * getListOfOrientations() {return listOfOrientations;}	// Get the ordered orientation of the reads. 1 means forward. 0 means reverse.
		Edge * getReverseEdge() {return reverseEdge;}	// Return the pointer to the reverse edge.
};

#endif /* EDGE_H_ */

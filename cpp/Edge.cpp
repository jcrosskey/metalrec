/*
 * Edge.cpp
 *
 * Created on: Fri Nov  7 15:13:39 EST 2014
 * Author: JJ Chai
 */

#include "Common.h"
#include "Edge.h"


/**********************************************************************************************************************
	Default Constructor
**********************************************************************************************************************/
Edge::Edge(void)
{

	// Initialize the variables.
	overlapOffset = 0;
	overlapOrientation = 10;
	transitiveRemovalFlag = false;
	flow = 0;
	coverageDepth = 0;
/*
	listOfReads = new vector<UINT64>;
	listOfReads->resize(listOfReads->size());						// Resize to reduce space.

	listOfOverlapLengths = new vector<UINT16>;
	listOfOverlapLengths->resize(listOfOverlapLengths->size());		// Resize to reduce space.

	listOfOrientations = new vector<UINT8>;
	listOfOrientations->resize(listOfOrientations->size());			// Resize to reduce space.
*/
}

/**********************************************************************************************************************
	Another Constructor
**********************************************************************************************************************/
Edge::Edge(Read *from, Read *to, UINT64 orient, UINT64 length)
{
	makeEdge(from, to, orient, length);
}

/**********************************************************************************************************************
 	 Another Constructor
**********************************************************************************************************************/
Edge::Edge(Read *from, Read *to, UINT64 orient, UINT64 length, vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets, vector<UINT8> * listOrientations)
{
	makeEdge(from, to, orient, length, listReads, listOverlapOffsets, listOrientations);
}

/**********************************************************************************************************************
	Default destructor
**********************************************************************************************************************/
Edge::~Edge()
{
	// Free the memory used by the current edge.
	delete listOfReads;
	delete listOfOverlapOffsets;
	delete listOfOrientations;
}

/**********************************************************************************************************************
	Function to insert a simple edge
**********************************************************************************************************************/
bool Edge::makeEdge(Read *from, Read *to, UINT64 orient, UINT64 length)
{
	source = from;
	destination = to;
	overlapOrientation = orient;
	overlapOffset = length;

	// Initialize variables.
	transitiveRemovalFlag = false;
	flow = 0;
	coverageDepth = 0;

	listOfReads = new vector<UINT64>;
	listOfReads->resize(listOfReads->size());						// Resize to reduce space.

	listOfOverlapOffsets = new vector<UINT16>;
	listOfOverlapOffsets->resize(listOfOverlapOffsets->size());		// Resize to reduce space.

	listOfOrientations = new vector<UINT8>;
	listOfOrientations->resize(listOfOrientations->size());			// Resize to reduce space.

	return true;
}

/**********************************************************************************************************************
	Function to add a composite edge
**********************************************************************************************************************/
bool Edge::makeEdge(Read *from, Read *to, UINT64 orient, UINT64 length,  vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets, vector<UINT8> * listOrientations)
{
	source = from;
	destination = to;
	overlapOrientation = orient;
	overlapOffset = length;
	transitiveRemovalFlag = false;
	flow = 0;
	coverageDepth = 0;
	listOfReads = listReads;
	listOfOverlapOffsets = listOverlapOffsets;
	listOfOrientations = listOrientations;
	listOfReads->resize(listOfReads->size());					// Resize to reduce space.
	listOfOverlapOffsets->resize(listOfOverlapOffsets->size());	// Resize to reduce space.
	listOfOrientations->resize(listOfOrientations->size());		// Resize to reduce space.

	return true;
}

/**********************************************************************************************************************
	Function to set the pointer to the reverse edge;
**********************************************************************************************************************/
bool Edge::setReverseEdge(Edge * edge)
{
	reverseEdge = edge;
	return true;
}

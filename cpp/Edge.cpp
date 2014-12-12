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
	transitiveRemovalFlag = false;
	flow = 0;
	coverageDepth = 0;
}

/**********************************************************************************************************************
	Another Constructor
**********************************************************************************************************************/
Edge::Edge(Read *from, Read *to, UINT64 length)
{
	makeEdge(from, to, length);
}

/**********************************************************************************************************************
 	 Another Constructor
**********************************************************************************************************************/
Edge::Edge(Read *from, Read *to, UINT64 length, vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets)
{
	makeEdge(from, to, length, listReads, listOverlapOffsets);
}

/**********************************************************************************************************************
	Default destructor
**********************************************************************************************************************/
Edge::~Edge()
{
	// Free the memory used by the current edge.
	delete listOfReads;
	delete listOfOverlapOffsets;
}

/**********************************************************************************************************************
	Function to insert a simple edge
**********************************************************************************************************************/
bool Edge::makeEdge(Read *from, Read *to, UINT64 length)
{
	source = from;
	destination = to;
	overlapOffset = length;

	// Initialize variables.
	transitiveRemovalFlag = false;
	flow = 0;
	coverageDepth = 0;

	listOfReads = new vector<UINT64>;
	listOfReads->resize(listOfReads->size());	// Resize to reduce space.

	listOfOverlapOffsets = new vector<UINT16>;
	listOfOverlapOffsets->resize(listOfOverlapOffsets->size());	// Resize to reduce space.

	return true;
}

/**********************************************************************************************************************
	Function to add a composite edge
**********************************************************************************************************************/
bool Edge::makeEdge(Read *from, Read *to, UINT64 length,  vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets)
{
	source = from;
	destination = to;
	overlapOffset = length;

	// Initialize variables.
	transitiveRemovalFlag = false;
	flow = 0;
	coverageDepth = 0;

	listOfReads = listReads;
	listOfReads->resize(listOfReads->size());	// Resize to reduce space.

	listOfOverlapOffsets = listOverlapOffsets;
	listOfOverlapOffsets->resize(listOfOverlapOffsets->size());	// Resize to reduce space.

	return true;
}

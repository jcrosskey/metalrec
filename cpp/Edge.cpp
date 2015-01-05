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
	endCoordinateLimit = 0;
}

/**********************************************************************************************************************
	Another Constructor
**********************************************************************************************************************/
Edge::Edge(Read *from, Read *to, UINT64 length, UINT16 numSub, vector<UINT64> *listSubs) 	// Another constructor.
{
	makeEdge(from, to, length, numSub, listSubs);
	endCoordinateLimit = 0;
}

/**********************************************************************************************************************
 	 Another Constructor
**********************************************************************************************************************/
Edge::Edge(Read *from, Read *to, UINT64 length, UINT16 numSub, vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets, vector<UINT64> *listSubs)
{
	makeEdge(from, to, length, numSub, listReads, listOverlapOffsets, listSubs);
	endCoordinateLimit = 0;
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
bool Edge::makeEdge(Read *from, Read *to, UINT64 length, UINT16 numSub, vector<UINT64> *listSubs)
{
	source = from;
	destination = to;
	overlapOffset = length;
	numOfSubstitions = numSub;

	// Initialize variables.
	transitiveRemovalFlag = false;
	flow = 0;
	coverageDepth = 0;

	listOfReads = new vector<UINT64>;
	listOfReads->resize(listOfReads->size());	// Resize to reduce space.

	listOfOverlapOffsets = new vector<UINT16>;
	listOfOverlapOffsets->resize(listOfOverlapOffsets->size());	// Resize to reduce space.

	listOfSubstitutionPoses = listSubs;
	listOfSubstitutionPoses->resize(listOfSubstitutionPoses->size());	// Resize to reduce space.

	return true;
}

/**********************************************************************************************************************
	Function to add a composite edge
**********************************************************************************************************************/
bool Edge::makeEdge(Read *from, Read *to, UINT64 length,  UINT16 numSub, vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets, vector<UINT64> *listSubs)
{
	source = from;
	destination = to;
	overlapOffset = length;
	numOfSubstitions = numSub;

	// Initialize variables.
	transitiveRemovalFlag = false;
	flow = 0;
	coverageDepth = 0;

	listOfReads = listReads;
	listOfReads->resize(listOfReads->size());	// Resize to reduce space.

	listOfOverlapOffsets = listOverlapOffsets;
	listOfOverlapOffsets->resize(listOfOverlapOffsets->size());	// Resize to reduce space.

	listOfSubstitutionPoses = listSubs;
	listOfSubstitutionPoses->resize(listOfSubstitutionPoses->size());	// Resize to reduce space.

	return true;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getStringLengthInRange
 *  Description:  Get the string length in the edge that is also in the range limit
 * =====================================================================================
 */
UINT64 Edge::getStringLengthInRange() 
{
	INT32 totalLength = getOverlapOffset() + getDestinationRead()->getReadLength();
	if (getSourceRead()->getStartCoord() < 0 )
	{
		totalLength = totalLength + getSourceRead()->getStartCoord();
	}
	if ((UINT32) (getDestinationRead()->getEndCoord()) > endCoordinateLimit )
	{
		totalLength = totalLength - (getDestinationRead()->getEndCoord() - endCoordinateLimit);
	}
	if (totalLength > 0)
		return (INT64) totalLength;
	else
		return 0;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  <
 *  Description:  Overloading comparing operator for edge object
 * =====================================================================================
 */
bool Edge::operator<(Edge & anotherEdge)
{
	return (getStringLengthInRange() < anotherEdge.getStringLengthInRange());
}

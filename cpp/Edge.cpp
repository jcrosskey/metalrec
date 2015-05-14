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
}

/**********************************************************************************************************************
	copy constructor
**********************************************************************************************************************/
Edge::Edge(const Edge & E)
{

	overlapOffset = E.overlapOffset;
	numOfSubstitions = E.numOfSubstitions;
	transitiveRemovalFlag = E.transitiveRemovalFlag;

	source = new Read;
	*source = *(E.source);

	destination = new Read;
	*destination = *(E.destination);

	size_t k = 0;
	listOfReads = new vector<UINT64>;
	listOfReads->resize((E.listOfReads)->size());
	for(k = 0; k < listOfReads->size(); k++){
		listOfReads->at(k) = (E.listOfReads)->at(k);
	}
	listOfOverlapOffsets = new vector<UINT16>;
	listOfOverlapOffsets->resize((E.listOfOverlapOffsets)->size());
	for(k = 0; k < listOfOverlapOffsets->size(); k++){
		listOfOverlapOffsets->at(k) = (E.listOfOverlapOffsets)->at(k);
	}
	listOfSubstitutionPoses = new vector<UINT64>;
	listOfSubstitutionPoses->resize((E.listOfSubstitutionPoses)->size());
	for(k = 0; k < listOfSubstitutionPoses->size(); k++){
		listOfSubstitutionPoses->at(k) = (E.listOfSubstitutionPoses)->at(k);
	}
}

/**********************************************************************************************************************
	copy assignment
**********************************************************************************************************************/
Edge & Edge::operator= (const Edge & E)
{

	overlapOffset = E.overlapOffset;
	numOfSubstitions = E.numOfSubstitions;
	transitiveRemovalFlag = E.transitiveRemovalFlag;

	delete this->source;
	source = new Read;
	*source = *(E.source);

	delete this->destination;
	destination = new Read;
	*destination = *(E.destination);

	size_t k = 0;
	delete this->listOfReads;
	listOfReads = new vector<UINT64>;
	listOfReads->resize((E.listOfReads)->size());
	for(k = 0; k < listOfReads->size(); k++){
		listOfReads->at(k) = (E.listOfReads)->at(k);
	}

	delete this->listOfOverlapOffsets;
	listOfOverlapOffsets = new vector<UINT16>;
	listOfOverlapOffsets->resize((E.listOfOverlapOffsets)->size());
	for(k = 0; k < listOfOverlapOffsets->size(); k++){
		listOfOverlapOffsets->at(k) = (E.listOfOverlapOffsets)->at(k);
	}

	delete this->listOfSubstitutionPoses;
	listOfSubstitutionPoses = new vector<UINT64>;
	listOfSubstitutionPoses->resize((E.listOfSubstitutionPoses)->size());
	for(k = 0; k < listOfSubstitutionPoses->size(); k++){
		listOfSubstitutionPoses->at(k) = (E.listOfSubstitutionPoses)->at(k);
	}
	return *this;
}

/**********************************************************************************************************************
	Another Constructor
**********************************************************************************************************************/
Edge::Edge(Read *from, Read *to, UINT64 length, UINT16 numSub, vector<UINT64> *listSubs) 	// Another constructor.
{
	source = from;
	destination = to;
	overlapOffset = length;
	numOfSubstitions = numSub;

	// Initialize variables.
	transitiveRemovalFlag = false;

	listOfReads = new vector<UINT64>;
	listOfReads->resize(listOfReads->size());	// Resize to reduce space.

	listOfOverlapOffsets = new vector<UINT16>;
	listOfOverlapOffsets->resize(listOfOverlapOffsets->size());	// Resize to reduce space.

	listOfSubstitutionPoses = new vector<UINT64>;
	listOfSubstitutionPoses->resize(listSubs->size());	// Resize to reduce space.
	for(size_t k = 0; k < listOfSubstitutionPoses->size(); k++){
		listOfSubstitutionPoses->at(k) = listSubs->at(k);
	}
}

/**********************************************************************************************************************
 	 Another Constructor
**********************************************************************************************************************/
Edge::Edge(Read *from, Read *to, UINT64 length, UINT16 numSub, vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets, vector<UINT64> *listSubs)
{
	source = from;
	destination = to;
	overlapOffset = length;
	numOfSubstitions = numSub;

	// Initialize variables.
	transitiveRemovalFlag = false;

	size_t k = 0;
	/* Copy the pointers and the pointed values to the Edge, the passed parameters can be deleted later. */
	listOfReads = new vector<UINT64>;
	listOfReads->resize(listReads->size());	// Resize to reduce space.
	for(k = 0; k < listOfReads->size(); k++){
		listOfReads->at(k) = listReads->at(k);
	}

	listOfOverlapOffsets = new vector<UINT16>;
	listOfOverlapOffsets->resize(listOverlapOffsets->size());	// Resize to reduce space.
	for(k = 0; k < listOfOverlapOffsets->size(); k++){
		listOfOverlapOffsets->at(k) = listOverlapOffsets->at(k);
	}
	listOfOverlapOffsets->resize(listOfOverlapOffsets->size());	// Resize to reduce space.

	listOfSubstitutionPoses = new vector<UINT64>;
	listOfSubstitutionPoses->resize(listSubs->size());	// Resize to reduce space.
	for(k = 0; k < listOfSubstitutionPoses->size(); k++){
		listOfSubstitutionPoses->at(k) = listSubs->at(k);
	}
}

/**********************************************************************************************************************
	Default destructor
**********************************************************************************************************************/
Edge::~Edge()
{
	// Free the memory used by the current edge.
	delete listOfReads;
	delete listOfOverlapOffsets;
	delete listOfSubstitutionPoses;
}


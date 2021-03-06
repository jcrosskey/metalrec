/*
 * OverlapGraph.cpp
 *
 * Created on: Tue Nov 11 12:45:42 EST 2014
 * Author: JJ Chai
 */

#include "OverlapGraph.h"

bool compareEdgeByStringLength (Edge *edge1, Edge *edge2)
{
	return (edge1->getOverlapOffset() + edge1->getDestinationRead()->getReadLength() > edge2->getOverlapOffset() + edge2->getDestinationRead()->getReadLength());
}

/*******************************************************
 * Function to compare two edges. Used for sorting.
 * Edges are sorted by the overlap offset.
 *******************************************************/
bool compareEdgeByOverlapOffset (Edge *edge1, Edge* edge2)
{
	return (edge1->getOverlapOffset() < edge2->getOverlapOffset());
}


/**********************************************************************************************************************
  Default Constructor
 **********************************************************************************************************************/
OverlapGraph::OverlapGraph(void)
{
}


/**********************************************************************************************************************
  Copy Constructor
 **********************************************************************************************************************/
OverlapGraph::OverlapGraph(const OverlapGraph & O)
{
	minimumOverlapLength = O.minimumOverlapLength;
	maxError = O.maxError;
	maxErrorRate = O.maxErrorRate;
	rubberPos = O.rubberPos;
	numberOfNodes = O.numberOfNodes;
	numberOfEdges = O.numberOfEdges;
	dataSet = new Dataset;
	*dataSet = *(O.dataSet);
	hashTable = new HashTable;
	*hashTable = *(O.hashTable);
}


/**********************************************************************************************************************
  Copy assignment
 **********************************************************************************************************************/
OverlapGraph & OverlapGraph::operator= (const OverlapGraph & O)
{
	minimumOverlapLength = O.minimumOverlapLength;
	maxError = O.maxError;
	maxErrorRate = O.maxErrorRate;
	rubberPos = O.rubberPos;
	numberOfNodes = O.numberOfNodes;
	numberOfEdges = O.numberOfEdges;
	delete dataSet;
	dataSet = new Dataset;
	*dataSet = *(O.dataSet);
	delete hashTable;
	hashTable = new HashTable;
	*hashTable = *(O.hashTable);
	return *this;
}


/**********************************************************************************************************************
  Another Constructor. Build the overlap grpah from the data_Set, and specified parameter values
 **********************************************************************************************************************/
OverlapGraph::OverlapGraph(HashTable *ht, const UINT64 & minOverlap, const UINT32 & max_Error, const float & max_ErrorRate, const UINT32 & rubber_pos)
{
	minimumOverlapLength = minOverlap;
	maxError = max_Error;
	maxErrorRate = max_ErrorRate;
	numberOfNodes = 0;
	numberOfEdges = 0;
	rubberPos = rubber_pos;
	buildOverlapGraphFromHashTable(ht);
}


/**********************************************************************************************************************
  Default destructor.
 **********************************************************************************************************************/
OverlapGraph::~OverlapGraph()
{
	// Free the memory used by the overlap graph.
	for(UINT64 i = 0; i < graph->size(); i++)
	{
		for(UINT64 j = 0; j< graph->at(i)->size(); j++)
		{
			delete graph->at(i)->at(j);
		}
		delete graph->at(i);
	}
	delete graph;
}


/**********************************************************************************************************************
 * Build the overlap graph from data set
 **********************************************************************************************************************/
bool OverlapGraph::buildOverlapGraphFromHashTable(HashTable *ht)
{
	CLOCKSTART;
	numberOfNodes = 0;
	numberOfEdges = 0;
	numberOfEdgesEverMade = 0;
	hashTable = ht;                         /* set hashtable */
	dataSet = ht->getDataset();             /* set dataset */
	UINT64 counter = 0;	// Number of reads explored so far (maybe)

	// initialize and reserve space for the type vectors
	vector<nodeType> *exploredReads = new vector<nodeType>;
	exploredReads->reserve(dataSet->getNumberOfUniqueReads()+1);

	vector<UINT64> * queue = new vector<UINT64>;
	queue->reserve(dataSet->getNumberOfUniqueReads()+1);

	vector<markType> *markedNodes = new vector<markType>;
	markedNodes->reserve(dataSet->getNumberOfUniqueReads()+1);

	graph = new vector< vector<Edge *> * >; /* initialize the graph */
	graph->reserve(dataSet->getNumberOfUniqueReads()+1);

	for(UINT64 i = 0; i <= dataSet->getNumberOfUniqueReads(); i++) // Initialization, one edge list for each unique read
	{
		vector<Edge *> *newList = new vector<Edge *>;	// Vector of edges
		graph->push_back(newList);	// Insert this vector (list) of edges into the graph
		exploredReads->push_back(UNEXPLORED);	// Initialize all reads to be UNEXPLORED
		queue->push_back(0); /* Initialize the queue that manages reads to explore */
		markedNodes->push_back(VACANT); /* Initialize all reads to be VACANT */
	}

	markContainedReads();
	buildReverseGraph();

	for(UINT64 i = 1; i <= dataSet->getNumberOfUniqueReads(); i++)	// i starts at 1, because read number starts with 1
	{
		if(exploredReads->at(i) == UNEXPLORED && !dataSet->getReadFromID(i)->isContainedRead()) /* explore the unexplored, and non-contained reads */
		{
			UINT64 start = 0, end = 0; 	// Initialize queue start and end.
			queue->at(end++) = i;
			while(start < end) 	// This loop will explore all connected component starting from read i.
			{
				counter++;	// Number of reads explored so far
				UINT64 read1 = queue->at(start++);
				/*  Explore node if not explored already */
				if(exploredReads->at(read1) == UNEXPLORED)
				{
					insertAllEdgesOfRead(read1, exploredReads);	// Explore current node, insert all its adjacent edges
					exploredReads->at(read1) = EXPLORED;	// This read is marked as explored
				}
				/* If the just explored node has some neighbors, continue to find transitive edges and remove them */
				if(graph->at(read1)->size() != 0) 	// Read has some edges (required only for the first read when a new queue starts).
				{
					if(exploredReads->at(read1) == EXPLORED) 	// Explore unexplored neighbors, so that the transitive edges adjacent to read1 can be marked.
					{
						for(UINT64 index1 = 0; index1 < graph->at(read1)->size(); index1++ )
						{
							UINT64 read2 = graph->at(read1)->at(index1)->getDestinationRead()->getID();
							if(exploredReads->at(read2) == UNEXPLORED) 	// Not explored.
							{
								queue->at(end++) = read2; 	// Put in the queue.
								insertAllEdgesOfRead(read2, exploredReads);	// Explore this unexplored neighbor
								exploredReads->at(read2) = EXPLORED;
							}
						}
						markTransitiveEdges(read1, markedNodes); // Mark transitive edges
						exploredReads->at(read1) = EXPLORED_AND_TRANSITIVE_EDGES_MARKED;
					}
					if(exploredReads->at(read1) == EXPLORED_AND_TRANSITIVE_EDGES_MARKED)	// Read has transitive edges marked, explore neighbors' neighbors so that transitive edges adjacent to read1 can be removed.
					{
						for(UINT64 index1 = 0; index1 < graph->at(read1)->size(); index1++) 	// Then explore all neighbors' neighbors
						{
							UINT64 read2 = graph->at(read1)->at(index1)->getDestinationRead()->getID();
							if(exploredReads->at(read2) == EXPLORED)
							{
								for(UINT64 index2 = 0; index2 < graph->at(read2)->size(); index2++) 	// Explore all neighbors neighbors
								{
									UINT64 read3 = graph->at(read2)->at(index2)->getDestinationRead()->getID();
									if(exploredReads->at(read3) == UNEXPLORED) 	// Not explored
									{
										queue->at(end++) = read3; 	// Put in the queue
										insertAllEdgesOfRead(read3, exploredReads);
										exploredReads->at(read3) = EXPLORED;
									}
								}
								markTransitiveEdges(read2, markedNodes); // Mark transitive edge
								exploredReads->at(read2) = EXPLORED_AND_TRANSITIVE_EDGES_MARKED;
							}
						}
						removeTransitiveEdges(read1); // Remove the transitive edges, function need to rewrite
					}
				}
				if(counter%100000==0)	// Show the progress.
				{
					FILE_LOG(logDEBUG)<<"counter: " << counter << " Nodes: " << numberOfNodes << " Edges: " << numberOfEdges;
				}
			}
		}
	}
	// report total number
	FILE_LOG(logINFO)<< "counter: " << counter << " Nodes: " << numberOfNodes << " Edges: " << numberOfEdges;
	FILE_LOG(logINFO)<< "number of edges made so far: " << numberOfEdgesEverMade;

	delete exploredReads;
	delete queue;
	delete markedNodes;


/* 	vector<Edge *> contigEdges;
 * 
 * 	counter = 0;
 * 	// contracting composite edges, do not remove dead end nodes for now 
 * 	if (numberOfEdges > 0)
 * 	{
 * 		if (loglevel > 3){
 * 			getEdges(contigEdges);
 * 			printGraph("debug.gdl",contigEdges);
 * 		}
 * 		unsigned int iteration = 0;
 * 		do
 * 		{
 * 
 * 			string prefix = "test" + Utils::intToString(iteration);
 * 			counter = contractCompositePaths();	// need to rewrite contractCompositePaths function
 * 			FILE_LOG(logINFO)<< "number of edges made so far: " << numberOfEdgesEverMade;
 * 			getEdges(contigEdges);
 * 			printGraph(prefix+".comp.gdl",contigEdges);
 * 
 * 			counter += removeDeadEndNodes();
 * 			getEdges(contigEdges);
 * 			printGraph(prefix+".dead.gdl",contigEdges);
 * 
 * 			counter += popBubbles();
 * 			iteration++;
 * 		} while (counter > 0);
 * 
 * 	}
 */
	CLOCKSTOP;
	return true;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  buildReverseGraph
 *  Description:  Store the connection (only) between nodes with edges in the reverse direction
 * =====================================================================================
 */
void OverlapGraph::buildReverseGraph(void) /* build the reverse graph */
{
	CLOCKSTART;
	FILE_LOG(logINFO) << "size of reverse graph is : " << reverseGraph.size();
	for(UINT64 k = 0; k <= dataSet->getNumberOfUniqueReads(); k++){
		vector<UINT64> edgeVector;
		reverseGraph.push_back(edgeVector);
	}
	FILE_LOG(logINFO) << "size of reverse graph is : " << reverseGraph.size();
	for(UINT64 k = 0; k < graph->size(); k++){
		if (graph->at(k)->size() > 0){
			UINT64 readNumber1 = k; /* source read ID */
			for (UINT64 i = 0; i < graph->at(k)->size(); i++){
				UINT64 readNumber2 = graph->at(k)->at(i)->getDestinationRead()->getID(); /* destination read ID */
				(reverseGraph.at(readNumber2)).push_back(readNumber1);
			}
		}
	}
	CLOCKSTOP;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  markContainedReads
 *  Description:  This function check if a read contains other small reads. 
 *  		  If a read is contained in more than one super read then it is assigned to the longest such super read.
 * =====================================================================================
 */
void OverlapGraph::markContainedReads(void)
{
	CLOCKSTART;
	if(dataSet->longestReadLength == dataSet->shortestReadLength) // If all reads are of same length, then no need to do look for contained reads.
	{
		FILE_LOG(logINFO) << "All reads are of same length " << dataSet->longestReadLength << ". No contained reads.";
		CLOCKSTOP;
		return;
	}
	UINT64 counter = 0;
	for(UINT64 i = 1; i < graph->size(); i++) // For each read, except the last one
	{
		Read *read1 = dataSet->getReadFromID(i); // Get the read
		if (read1->superReadID == 0)
		{
			FILE_LOG(logDEBUG3) << "Check contained reads with read " << i;
			string readString = read1->getDnaStringForward(); // Get the forward of the read
			string subString; /* sub string of read1 with length hashStringLength */
			UINT64 readNumber = read1->getID();
			for(UINT64 j = 0; j < read1->getReadLength() - hashTable->getHashStringLength(); j++) // for each substring of read1 of length getHashStringLength
			{
				subString = readString.substr(j,hashTable->getHashStringLength()); // Get the substring from read1
				vector<UINT64> * listOfReads=hashTable->getListOfReads(subString); // Search the substring in the hash table
				//FILE_LOG(logDEBUG3) << "\tNumber of reads with this string hashed: " << listOfReads->size() << " at position " << j;
				if(!listOfReads->empty()) // If other reads contain the substring as prefix or suffix
				{
					for(UINT64 k = 0; k < listOfReads->size(); k++) // For each read in the list.
					{
						UINT64 data = listOfReads->at(k); // We used bit operation in the hash table to store read ID and orientation
						Read *read2 = dataSet->getReadFromID(data & 0X3FFFFFFFFFFFFFFF); 	// Least significant 62 bits store the read number.
						if (read2->superReadID == 0) /* Only check the reads that are not marked as contained already */
						{
							UINT64 readNumber2 = read2->getID();
							// Most significant 1 bits store the orientation.
							// Orientation 0 means prefix of forward string of the read
							// Orientation 1 means suffix of forward string of the read

							bool contain_possible = false; /* check if overlap between the 2 reads is possible or not, based on the coordinates */
							if ( readNumber < readNumber2 ) {
								contain_possible = ((read1->getEndCoord() - read2->getStartCoord() + rubberPos) > 0);
							}
							else if (readNumber > readNumber2)
							{
								contain_possible = ((read2->getEndCoord() - read1->getStartCoord() + rubberPos) > 0);
							}
							if( readString.length() > read2->getDnaStringForward().length() /* read1 has to be longer than read2, also eliminate the possibility of comparing a read to itself */
									&& contain_possible /* they have to have the possibility to contain first */
									&& checkOverlapForContainedRead(read1,read2,(data >> 62),j)) // Check if the remaining of the strings also match
							{       /* What if read1 is contained in some other read? */
								if(read2->superReadID == 0) // This is the first super read found. we store the ID of the super read.
								{
									FILE_LOG(logDEBUG3) << "\tFound read " << readNumber2 << " contained in read " << readNumber;
									read2->superReadID = i;
									counter ++;
								}
								else // Already found some previous super read. Now we have another super read.
								{
									if(readString.length() > dataSet->getReadFromID(read2->superReadID)->getReadLength()) // This super read is longer than the previous super read. Update the super read ID.
										read2->superReadID = i;
								}
							}
						}
					}
				}
			}
		}
		if(i % 1000000 == 0)
		{
			FILE_LOG(logDEBUG) << counter << " contained reads in " << i << " super reads.";
		}
	}

	// Get some statistics
	UINT64 containedReads = 0, nonContainedReads = 0;
	for(UINT64 i = 1; i < graph->size(); i++)
	{
		Read *rr = dataSet->getReadFromID(i);
		if(rr->superReadID == 0) // Count the number of reads that are not contained by some other reads.
			nonContainedReads++;
		else					// Count the number of reads that are contained by some other read.
			containedReads++;
	}
	FILE_LOG(logINFO) << nonContainedReads << " Non-contained reads.";
	FILE_LOG(logINFO)<< containedReads << " contained reads.";
	CLOCKSTOP;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  checkOverlapForContainedRead
 *  Description:  Hash table search found that a proper substring of read1 is a prefix or suffix of read2 
 *  		need to also check the starting coordinate so that if two reads are not around the same region, don't attempt to do this 
 *  		orient 0 means prefix of forward of the read2 
 *  		orient 1 means suffix of forward of the read2 
 *  		We need to check if the remaining of the stings match to see if read2 is contained in read1.
 * =====================================================================================
 */
bool OverlapGraph::checkOverlapForContainedRead(Read *read1, Read *read2, UINT64 orient, UINT64 start)
{
	// start: index of the position where the prefix/suffix is found in read1
	string string1=read1->getDnaStringForward(); // Get the forward of read1
	UINT64 hashStringLength = hashTable->getHashStringLength(); 
	UINT64 lengthRemaining1, lengthRemaining2;
	string string2 = read2->getDnaStringForward(); // Get the string in read2 based on the orientation.
	UINT16 numSub;
	vector<UINT64> * listSubs = new vector<UINT64>;
	bool contained = false;
	if(orient == 0)
		// orient 0 read1 has read2's prefix
		//   >--------MMMMMMMMMMMMMMM*******------> read1      M means match found by hash table
		//            MMMMMMMMMMMMMMM*******>       read2      * means we need to check these characters for match
	{
		lengthRemaining1 = string1.length() - start - hashStringLength; 	// This is the remaining of read1
		lengthRemaining2 = string2.length() - hashStringLength; 	// This is the remaining of read2
		if(lengthRemaining1 >= lengthRemaining2)
		{
			contained = checkOverlapWithSub(string1.substr(start + hashStringLength, lengthRemaining2), string2.substr(hashStringLength, lengthRemaining2), numSub, listSubs, orient); // If the remaining of the string match, then read2 is contained in read1
		}
	}
	else                                 // orient 1, read1 has read2's suffix
		//   >---*****MMMMMMMMMMMMMMM-------------> read1      M means match found by hash table
		//      >*****MMMMMMMMMMMMMMM       		read2      * means we need to check these characters for match
	{
		lengthRemaining1 = start;
		lengthRemaining2 = string2.length() - hashStringLength;
		if(lengthRemaining1 >= lengthRemaining2)
		{
			contained = checkOverlapWithSub(string1.substr(start - lengthRemaining2, lengthRemaining2),  string2.substr(0, lengthRemaining2), numSub, listSubs, orient); // If the remaining of the string match, then read2 is contained in read1
		}
	}
	delete listSubs;
	return contained;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  checkOverlap
 *  Description:  Checks if two read overlaps.  
 *  		Hash table search found that a proper substring of read1 is a prefix or suffix of read2 
 *  		Orientation 0 means prefix of forward of the read2 
 *  		Orientation 1 means suffix of forward of the read2
 * =====================================================================================
 */
bool OverlapGraph::checkOverlap(Read *read1, Read *read2, UINT64 orient, UINT64 start, UINT16 & numSub, vector<UINT64> * listSubs)
{
	string string1=read1->getDnaStringForward(); // Get the forward string of read1
	UINT64 hashStringLength = hashTable->getHashStringLength();
	string string2 = read2->getDnaStringForward(); // Get the string from read2. 
	if(orient == 0)		// orient 0
		//   >--------MMMMMMMMMMMMMMM*************> 			read1      M means match found by hash table
		//            MMMMMMMMMMMMMMM*************------->      read2      * means we need to check these characters for match
	{
		if(string1.length()- start - hashStringLength >= string2.length() - hashStringLength || string1.length() - start < minimumOverlapLength) // The overlap must continue till the end, also overlap length has to be as long as the minimumOverlapLength, assuming that they overlap
			return false;
		return checkOverlapWithSub(string1.substr(start + hashStringLength, string1.length()-(start + hashStringLength)), string2.substr(hashStringLength,  string1.length()-(start + hashStringLength)), numSub, listSubs, orient); // If the remaining strings match.
	}
	else								// orient 1
		//   	>********MMMMMMMMMMMMMMM-------------> 			read1      M means match found by hash table
		//  >----********MMMMMMMMMMMMMMM       		    		read2      * means we need to check these characters for match
	{
		if(string2.length()-hashStringLength < start || start + hashStringLength < minimumOverlapLength)
			return false;
		return checkOverlapWithSub(string1.substr(0, start), string2.substr(string2.length()-hashStringLength-start, start), numSub, listSubs, orient); // If the remaining strings match.
	}
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  checkOverlapWithSub
 *  Description:  check if two string overlap with number of substitutions allowed, perfectly matched prefix/suffix is not included in these strings. 
 *  		  In other words, str1 and str2 are the string to compare before/after the perfectly matched hash string
 * =====================================================================================
 */
bool OverlapGraph::checkOverlapWithSub(const string & str1, const string & str2, UINT16 & numSub, vector<UINT64> * listSubs, UINT64 orient)
{
	UINT64 match_len = hashTable->getHashStringLength(); /* number of perfect matched bases */
	UINT64 comp_len = min(str1.length(), str2.length() ); /* number of characters need to compare at the most */
	UINT64 maxMismatch = max(UINT32 (maxErrorRate * ( match_len + comp_len)), maxError ); /* Total overlap length between 2 reads, including the hash string */
	UINT64 mismatchPosStart = (orient == 0 ? hashTable->getHashStringLength() : 0 ); /* substitution positions count from after the overlap offset */
	/* For example: (MMMMs are the matched hash string, m is a match, and here there is a mismatch between A and G
	 * coord:               012345678901234567
	 * Read1 >--------------MMMMMMMMMmmAmmmmmm
	 * Read2                MMMMMMMMMmmGmmmmmm------------->
	 * In this case, the substitution position will be 11
	 * */

	numSub = 0;                        /* number of substitutions so far */
	listSubs->clear(); /* Clear the content of the list of substitution positions */

	bool overlap = true;
	for ( UINT64 i=0; i < comp_len; i++ ) {
		if (str1.at(i) != str2.at(i))
		{
			numSub++;
			listSubs->push_back(mismatchPosStart + i); /* list of substitution positions: relative to the beginning of the overlap */
			if (numSub > maxMismatch)
			{
				overlap = false;
				break;
			}
		}
	}
	
	if ( overlap && numSub > 0 ) {
		FILE_LOG(logDEBUG3) << "Found overlap with " << numSub << " substitutions and overlap length " << match_len + comp_len;
	}
	return overlap;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  markTransitiveEdges
 *  Description:  Mark all the transitive edges of a read.  
 *  For Details: E.W. Myers. The fragment assembly string graph. Bioinformatics, 21(suppl 2):ii79-ii85, 2005.
 * =====================================================================================
 */
bool OverlapGraph::markTransitiveEdges(UINT64 readNumber, vector<markType> * markedNodes)
{

	for(UINT64 i = 0; i < graph->at(readNumber)->size(); i++) // Mark all the neighbours of the current read as INPLAY
		markedNodes->at(graph->at(readNumber)->at(i)->getDestinationRead()->getID()) = INPLAY; // Inplay

	for(UINT64 i = 0; i < graph->at(readNumber)->size(); i++) // Traverse through the list of edges according to their overlap offset.
	{
		UINT64 read2 = graph->at(readNumber)->at(i)->getDestinationRead()->getID(); // For each neighbor
		if(markedNodes->at(read2) == INPLAY) 	// If the neighbor is marked as INPLAY
		{
			for(UINT64 j = 0; j < graph->at(read2)->size(); j++)
			{
				UINT64 read3 = graph->at(read2)->at(j)->getDestinationRead()->getID(); // Get the neighbors neighbors
				if(markedNodes->at(read3) == INPLAY)
				{
					markedNodes->at(read3) = ELIMINATED; 	// Mark as ELIMINATED
				}
			}
		}
	}
	for(UINT64 i = 0;i < graph->at(readNumber)->size(); i++)
	{
		if(markedNodes->at(graph->at(readNumber)->at(i)->getDestinationRead()->getID()) == ELIMINATED) // Current read to a node marked as ELIMINATED
		{
			graph->at(readNumber)->at(i)->transitiveRemovalFlag = true; 	// Mark this edge as transitive edge. Will remove this edge later.
		}
	}

	for(UINT64 i = 0; i < graph->at(readNumber)->size(); i++)
		markedNodes->at(graph->at(readNumber)->at(i)->getDestinationRead()->getID()) = VACANT; 	// Change back all the variables modified in this function to VACANT

	markedNodes->at(readNumber) = VACANT; 	// Mark as vacant.
	return true;
}


/**********************************************************************************************************************
  Insert an edge in the overlap graph.
 **********************************************************************************************************************/
bool OverlapGraph::insertEdge(Edge * edge)
{

	UINT64 SourceID = edge->getSourceRead()->getID(); // This is the source read.
	UINT64 DestinationID = edge->getDestinationRead()->getID(); // This is the destination read.
	if((dataSet->getReadFromID(SourceID)->numInEdges + dataSet->getReadFromID(SourceID)->numOutEdges) == 0)	// If there is no edge incident to the node
		numberOfNodes++;	// Then a new node is inserted in the graph. Number of nodes increased.
	if((dataSet->getReadFromID(DestinationID)->numOutEdges + dataSet->getReadFromID(DestinationID)->numInEdges) == 0)
		numberOfNodes++;

	// update the number of in edges and out edges incident to the corresponding node
	dataSet->getReadFromID(SourceID)->numOutEdges++;
	dataSet->getReadFromID(DestinationID)->numInEdges++;

	graph->at(SourceID)->push_back(edge);	// Insert the edge in the list of edges of ID
	numberOfEdges++;	// Increase the number of edges.
	updateReadLocations(edge);	// If the current edge contains some reads, then we need to update their location information.
	reverseGraph.at(DestinationID).push_back(SourceID);
	numberOfEdgesEverMade++;
	edge->setID(numberOfEdgesEverMade);
	return true;
}


/**********************************************************************************************************************
  Insert an edge in the graph.
 **********************************************************************************************************************/
bool OverlapGraph::insertEdge(Read *read1, Read *read2, UINT16 overlapOffset, UINT16 numSub, vector<UINT64> *listSubs)
{
	FILE_LOG(logDEBUG4) << " Insert edge from " << read1->getID() << " to " << read2->getID() << " with overlap offset " << overlapOffset;
	Edge * edge1 = new Edge(read1,read2,overlapOffset, numSub, listSubs);	// Create a new edge in the graph to insert.
	insertEdge(edge1);	// Insert the edge in the overlap graph.
	return true;
}


/**********************************************************************************************************************
  Contract composite paths in the overlap graph.
 **********************************************************************************************************************/
UINT64 OverlapGraph::contractCompositePaths(void)
{
	CLOCKSTART;
	UINT64 counter = 0;
	UINT64 index = 0;
	bool found_composite_read = false;
	while (index < graph->size())
	{
		if(graph->at(index)->size() > 0)	// If this read has some edges adjacent to it
		{
			found_composite_read = false;
			for(UINT64 i = 0; i < graph->at(index)->size(); i++)	// Check all the edges 
			{
				Read *read = graph->at(index)->at(i)->getDestinationRead(); /* Destination read */

				if (read->numInEdges == 1 && read->numOutEdges == 1)	// If the destination read only has 1 edge in and 1 edge out, there are a pair of composite edges
				{
					insertEdge(mergeEdges(graph->at(index)->at(i), graph->at(read->getID())->at(0)));	// merge the two edges and insert into the graph, update the related reads' information, need to check
					removeEdge(graph->at(index)->at(i));	// delete the two edges that were contracted
					removeEdge(graph->at(read->getID())->at(0));
					counter++;
					found_composite_read = true;
				}
			}
			if (!found_composite_read)
				index++;
		}
		else
			index++;
	}
	FILE_LOG(logINFO) << counter << " composite Edges merged.";
	CLOCKSTOP;
	return counter;
}


/**********************************************************************************************************************
  remove an edge from the overlap graph.
 **********************************************************************************************************************/
bool OverlapGraph::removeEdge(Edge *edge)
{
	removeReadLocations(edge);	// If the current edge contains some reads. We have to update their location formation.
	UINT64 ID1 = edge->getSourceRead()->getID();  // Get the source read ID.
	UINT64 ID2 = edge->getDestinationRead()->getID();  // Get the destation read ID.

	for(UINT64 i = 0; i< graph->at(ID1)->size(); i++) // Delete the edge.
	{
		if(graph->at(ID1)->at(i) == edge)
		{
			delete graph->at(ID1)->at(i);
			graph->at(ID1)->at(i) = graph->at(ID1)->at(graph->at(ID1)->size()-1);	// Move the last edge to this entry
			graph->at(ID1)->pop_back();	// delete the last entry
			dataSet->getReadFromID(ID1)->numOutEdges--;
			dataSet->getReadFromID(ID2)->numInEdges--;
			if((dataSet->getReadFromID(ID1)->numInEdges + dataSet->getReadFromID(ID1)->numOutEdges) == 0)
				numberOfNodes--;
			if((dataSet->getReadFromID(ID2)->numInEdges + dataSet->getReadFromID(ID2)->numOutEdges) == 0)
				numberOfNodes--;
			numberOfEdges--;
			break;
		}
	}
	for(UINT64 i = 0; i < reverseGraph.at(ID2).size(); i++){
		if(reverseGraph.at(ID2).at(i) == ID1){
			reverseGraph.at(ID2).at(i) = reverseGraph.at(ID2).at(reverseGraph.at(ID2).size()-1);
			reverseGraph.at(ID2).pop_back();
			break;
		}
			
	}
	return true;
}



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getEdges
 *  Description:  save all the edges in the graph in a vector of pointers to these edges
 * =====================================================================================
 */
bool OverlapGraph::getEdges(vector<Edge *> & contigEdges)
{
	CLOCKSTART;
	contigEdges.clear();
	for(UINT64 i = 1; i<= dataSet->getNumberOfUniqueReads(); i++)
	{
		if(!graph->at(i)->empty()) // if this read has some edge(s) going out of it (since now the graph is directed)
		{
			for(UINT64 j = 0; j < graph->at(i)->size(); j++)
			{
				Edge * e = graph->at(i)->at(j);
				if( e->getDestinationRead()->superReadID==0 || e->getSourceRead()->superReadID==0 ) /* Only consider the edges between non-contained reads */
				{
					contigEdges.push_back(e); // List of contigs.
					//e->setEndCorrdinateLimit(dataSet->getPacBioReadLength());
				}
			}
		}
	}
	CLOCKSTOP;
	return true;
}

/**********************************************************************************************************************
 * This function prints the overlap graph in overlap_graph->gdl file. 
 * The graph can be viewed by aisee (free software available at http://www.aisee.com/) 
 * It also stores the contigs in a file.
 **********************************************************************************************************************/
bool OverlapGraph::printGraph(string graphFileName, const vector<Edge *> & contigEdges)
{
	CLOCKSTART;
	string edgeColor;
	ofstream graphFilePointer; 

	/************************* Store the graph in a file. ************************/
	graphFilePointer.open(graphFileName.c_str());
	if(!graphFilePointer.is_open())
		MYEXIT("Unable to open file: "+graphFileName);

    graphFilePointer << "source\ttarget\tnumReads\toffset\n";
    for (UINT64 i = 0; i < contigEdges.size(); i++)
    {
        Edge *e = contigEdges.at(i);
        int num_reads = 0;
        if(e->getListOfReads())
            num_reads = e->getListOfReads()->size();
        UINT64 source = e->getSourceRead()->getID(), destination = e->getDestinationRead()->getID();
        // Edge label: (first overlap length, edge length, number of reads, overlap offset, last overlap length)
        graphFilePointer << source << "\t" << destination << "\t" 
            << num_reads << "\t" << e->getOverlapOffset() << "\n";
    }

	graphFilePointer.close();
	FILE_LOG(logINFO) << "Graph written." << endl;
	CLOCKSTOP;
	return true;
	/************************* Store the graph in a file done. ************************/
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  printContigs
 *  Description:  Print contigs(strings spelled by edges) to file
 * =====================================================================================
 */
bool OverlapGraph::printContigs(string outputFastaName, vector<Edge *> & contigEdges, bool longestOnly)
{
	/************************* Store the contigs in a file. ************************/
	if (contigEdges.size() > 0) // Sort the contigs by their length if there are edges in the graph.
	{
		sort(contigEdges.begin(),contigEdges.end(),compareEdgeByStringLength);	// Sort the contigs by their total string length, in decreasing order
		/* Output the longest one for result */
		ofstream outputContigFilePointer;
		outputContigFilePointer.open(outputFastaName.c_str());
		if(!outputContigFilePointer.is_open())
			MYEXIT("Unable to open file: " + outputFastaName);
		if (longestOnly || contigEdges.size() == 1) /* Only the longest one is written or there is only 1 edge in the graph */
		{
			string s = getStringInEdge(contigEdges.at(0)); // get the string in the longest edge. This function need to be rewritten too.
			outputContigFilePointer << ">" << dataSet->getPacBioReadName() << " Edge ("  << contigEdges.at(0)->getSourceRead()->getID() << ":"<< contigEdges.at(0)->getSourceRead()->getStartCoord() << ", " << contigEdges.at(0)->getDestinationRead()->getID() << ":" << contigEdges.at(0)->getDestinationRead()->getEndCoord() << ") String Length: " << s.length() << endl; /* If only the longest one is printed, then it's named the same as the PacBio read name */
			UINT32 start=0;
			do
			{
				outputContigFilePointer << s.substr(start, 100) << endl;  // save 100 BP in each line.
				start+=100;
			} while (start < s.length());
		}
		else 
		{
			for (size_t i = 0; i< contigEdges.size(); i++)
			{
				string s = getStringInEdge(contigEdges.at(i)); // get the string in the longest edge. This function need to be rewritten too.
				outputContigFilePointer << ">contig_" << i+1  << " Edge ("  << contigEdges.at(i)->getSourceRead()->getID() << ":"<< contigEdges.at(i)->getSourceRead()->getStartCoord() << ", " << contigEdges.at(i)->getDestinationRead()->getID() << ":" << contigEdges.at(i)->getDestinationRead()->getEndCoord() << ") String Length: " << s.length() << ". Contains " << contigEdges.at(i)->getListOfOverlapOffsets()->size() << " reads" << endl;

				UINT32 start=0;
				do
				{
					outputContigFilePointer << s.substr(start, 100) << endl;  // save 100 BP in each line.
					start+=100;
				} while (start < s.length());
			}
			
		}
		outputContigFilePointer.close();
	}
	else
	{
		FILE_LOG(logWARNING) << "No contigs (edges) left in the graph, nothing to write";
	}
	return true;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  insertAllEdgesOfRead
 *  Description:  Insert all edges of a read in the overlap graph 
 *  		  If a read is already explored, it won't be explored against for another read again.
 * =====================================================================================
 */
bool OverlapGraph::insertAllEdgesOfRead(UINT64 readNumber, vector<nodeType> * exploredReads)
{
	Read *read1 = dataSet->getReadFromID(readNumber); 	// Get the current read read1.
	string readString = read1->getDnaStringForward(); 		// Get the forward string of read1.
	string subString;
	/* Shoud the scan of the string start from the second base, as in original OMEGA? QQQ 
	 * For Now, we'll start from the beginning. Think about it, will this result in duplicates? */
	for(UINT64 j = 1; j < read1->getReadLength()-hashTable->getHashStringLength() - 1; j++) // For each proper substring of length getHashStringLength of read1
	{
		subString = readString.substr(j,hashTable->getHashStringLength());  // Get the proper substring s of read1.
		vector<UINT64> * listOfReads=hashTable->getListOfReads(subString); // Search the string in the hash table.
		if(!listOfReads->empty()) // If there are some reads that contain s as prefix or suffix of the read or their reverse complement
		{
			for(UINT64 k = 0; k < listOfReads->size(); k++) // For each such reads.
			{
				UINT64 data = listOfReads->at(k);			// We used bit operations in the hash table. Most significant 2 bits store orientation and least significant 62 bits store read ID.
				UINT16 overlapOffset; 
				UINT64 coordOffset;
				UINT8 orientation;
				Read *read2 = dataSet->getReadFromID(data & 0X3FFFFFFFFFFFFFFF); 	// Least significant 62 bits store the read number.
				UINT64 readNumber2 = read2->getID();

				/* See if the edge between the two reads has already been found, if so, do not try to find the same edge again.
				 * For the overlap between two reads, if there is mismatch, there are two chances to find them. 
				 * Cannot simply depend on if one node has been explored or not. */
//				if(isEdgePresent(readNumber, readNumber2) || isEdgePresent(readNumber2, readNumber))
//					continue;

				/* If no errors are allowed, then every pair of reads only need to be checked once */
				if ( maxError == 0 && maxErrorRate == 0 && exploredReads->at(readNumber2)!= UNEXPLORED )
					continue;

				UINT16 numSub;
				vector<UINT64> *listSubs = new vector<UINT64>;

				bool overlap_possible = false; /* check if overlap between the 2 reads is possible or not, based on the coordinates */
				/* Do not consider the overlap between a read and itself */
				if ( readNumber < readNumber2 ) {
					overlap_possible = ((read1->getEndCoord() - read2->getStartCoord() + rubberPos) > 0);
				}
				else
				{
					overlap_possible = ((read2->getEndCoord() - read1->getStartCoord() + rubberPos) > 0);
				}
				if(read1->superReadID == 0 && read2->superReadID == 0 /* Both read need to be non-contained. */
						&& readNumber != readNumber2 /* Do not look for overlap with itself */
						&& overlap_possible /* They have to have the possibility to contain first */
						&& checkOverlap(read1, read2, (data >> 62), j, numSub, listSubs)) /* Check if they have valid overlap */
				{
					orientation = data >> 62;
					
					/* overlap found must be within +- 2*rubberPos distance from the coordinate offset */
					// PacBio ++++++++++++++++++++++++++++++++++++++++++++++
					// read1 	---------MMMMMMMMMMMM
					// read2		 MMMMMMMMMMMM-----------
					if ( orientation == 0 && readNumber < readNumber2 ) {
						coordOffset = read2->getStartCoord() - read1->getStartCoord();
						overlapOffset =  j;
						/* make sure the distance between overlap offset and coordinate offset is smaller than 2*rubberpos */
						if ((coordOffset + 2*rubberPos > overlapOffset) && (overlapOffset > ((coordOffset > 2*rubberPos)?(coordOffset - 2*rubberPos):(2*rubberPos-coordOffset))))
						{
							insertEdge(read1, read2, overlapOffset, numSub, listSubs);
							FILE_LOG(logDEBUG4) << "coord offset: " << coordOffset << "\t overlapOffset: " << overlapOffset;
							FILE_LOG(logDEBUG4) << read1->getDnaStringForward();
							FILE_LOG(logDEBUG4) << string(overlapOffset, '-') + read2->getDnaStringForward();
							FILE_LOG(logDEBUG4) << " ";
						}
					}
					// PacBio ++++++++++++++++++++++++++++++++++++++++++++++
					// read1		 MMMMMMMMMMMM-----------
					// read2 	---------MMMMMMMMMMMM
					else if ( orientation == 1 && readNumber > readNumber2)
					{
						coordOffset = read1->getStartCoord() - read2->getStartCoord();
						overlapOffset = read2->getReadLength() - hashTable->getHashStringLength() - j;
						if ((coordOffset + 2*rubberPos > overlapOffset) && (overlapOffset > ((coordOffset > 2*rubberPos)?(coordOffset - 2*rubberPos):(2*rubberPos-coordOffset))))
						{
							insertEdge(read2, read1, overlapOffset, numSub, listSubs);
							FILE_LOG(logDEBUG4) << "coord offset: " << coordOffset << "\t overlapOffset: " << overlapOffset;
							FILE_LOG(logDEBUG4) << string(overlapOffset, '-') + read1->getDnaStringForward();
							FILE_LOG(logDEBUG4) << read2->getDnaStringForward();
							FILE_LOG(logDEBUG4) << " ";
						}
					}
					// PacBio ++++++++++++++++++++++++++++++++++++++++++++++
					// read1		 MMMMMMMMMMMM-----------
					// read2 	          ---------MMMMMMMMMMMM
					else if (orientation == 1 && readNumber < readNumber2)
					{
						coordOffset = read2->getStartCoord() - read1->getStartCoord();
						overlapOffset = read2->getReadLength() - hashTable->getHashStringLength() - j;
						if (2*rubberPos > coordOffset &&  2*rubberPos-coordOffset > overlapOffset)
						{
							insertEdge(read2, read1, overlapOffset, numSub, listSubs);
							FILE_LOG(logDEBUG4) << "coord offset: " << coordOffset << "\t overlapOffset: " << overlapOffset;
							FILE_LOG(logDEBUG4) << string(overlapOffset, '-') + read1->getDnaStringForward();
							FILE_LOG(logDEBUG4) << read2->getDnaStringForward();
							FILE_LOG(logDEBUG4) << " ";
						}
					}
					// PacBio ++++++++++++++++++++++++++++++++++++++++++++++
					// read1 	          ---------MMMMMMMMMMMM
					// read2		 MMMMMMMMMMMM-----------
					// readNumber > readNumber2 and orientation == 0
					else {
						coordOffset = read1->getStartCoord() - read2->getStartCoord();
						overlapOffset = j;
						if (2*rubberPos > coordOffset &&  2*rubberPos-coordOffset > overlapOffset)
						{
							insertEdge(read1, read2, overlapOffset, numSub, listSubs);
							FILE_LOG(logDEBUG4) << "coord offset: " << coordOffset << "\t overlapOffset: " << overlapOffset;
							FILE_LOG(logDEBUG4) << read1->getDnaStringForward();
							FILE_LOG(logDEBUG4) << string(overlapOffset, '-') + read2->getDnaStringForward();
							FILE_LOG(logDEBUG4) << " ";
						}
						
					}
				}
				delete listSubs;
			}
		}
	}
	if(graph->at(readNumber)->size() != 0)
		sort(graph->at(readNumber)->begin(),graph->at(readNumber)->end(), compareEdgeByOverlapOffset); // Sort the list of edges of the current node according to the overlap offset (ascending).
	return true;
}


/**********************************************************************************************************************
  Remove all transitive edges of a given read.
  For Details: E.W. Myers. The fragment assembly string graph. Bioinformatics, 21(suppl 2):ii79-ii85, 2005.
 **********************************************************************************************************************/
bool OverlapGraph::removeTransitiveEdges(UINT64 readNumber)
{
	UINT64 j=0;
	for(UINT64 index=0; index < graph->at(readNumber)->size(); index++) // We will remove all the transitive edges of the current read.
	{
		if(graph->at(readNumber)->at(index)->transitiveRemovalFlag == false)	// We move all the non-transitive edges at the beginning of the list
			graph->at(readNumber)->at(j++) = graph->at(readNumber)->at(index);
		else	// Free the transitive edge
		{
			numberOfEdges--;
			UINT64 ID2 = graph->at(readNumber)->at(index)->getDestinationRead()->getID();
			dataSet->getReadFromID(readNumber)->numOutEdges--;
			dataSet->getReadFromID(ID2)->numInEdges--;
			if ( (dataSet->getReadFromID(ID2)->numInEdges + dataSet->getReadFromID(ID2)->numOutEdges)==0)
				numberOfNodes--;
			delete graph->at(readNumber)->at(index);
			for(UINT64 i = 0; i < reverseGraph.at(ID2).size(); i++){
				if(reverseGraph.at(ID2).at(i) == readNumber){
					reverseGraph.at(ID2).at(i) = reverseGraph.at(ID2).at(reverseGraph.at(ID2).size()-1);
					reverseGraph.at(ID2).pop_back();
					break;
				}

			}
		}
	}
	graph->at(readNumber)->resize(j);
	if ( (dataSet->getReadFromID(readNumber)->numInEdges + dataSet->getReadFromID(readNumber)->numOutEdges)==0)
		numberOfNodes--;
	return true;
}


/**********************************************************************************************************************
  Remove all edges whose source or destination read has read number between (endpoints included) the two given numbers
 **********************************************************************************************************************/
bool OverlapGraph::removeEdgesBetweenReadNumbers(UINT64 readNumber1, UINT64 readNumber2)
{
	FILE_LOG(logDEBUG3) << "Remove edges between read numebers " << readNumber1 << " and " << readNumber2;
	for (UINT64 index = 0; index < graph->size(); index++)	// Going through all the edges
	{
		if (graph->at(index)->size() > 0)
		{
			for( UINT64 j = 0; j < graph->at(index)->size(); j++)
			{
				Edge *edge = graph->at(index)->at(j);
				UINT64 sourceID = edge->getSourceRead()->getID();
				UINT64 destID = edge->getDestinationRead()->getID();
				/* If either of the source or destination is in the given range */
				if ((sourceID >= readNumber1 && sourceID <= readNumber2) || (destID >= readNumber1 && destID <= readNumber2))
				{
					FILE_LOG(logDEBUG3) << "Remove edges between read numebers " << sourceID << " and " << destID;
					removeEdge(edge);
					j--;    /* this is because when an edge is removed, the last edge incident to the node is moved to this position */
				}
			}
		}
	}
	return true;
}


/**********************************************************************************************************************
  Merge two edges in the overlap graph.
 **********************************************************************************************************************/
Edge * OverlapGraph::mergeEdges(Edge *edge1, Edge *edge2)
{
	//Edge *newEdge = new Edge();
	Read *read1 = edge1->getSourceRead(), *read2 = edge2->getDestinationRead();

	vector<UINT64> * listReadsForward = new vector<UINT64>;	// List of reads in the forward edge.
	vector<UINT16> * listOverlapOffsetsForward= new vector<UINT16>;	// List of Overlaps in the forward edge.
	vector<UINT64> * listOfSubstitutionPoses= new vector<UINT64>;	// List of Overlaps in the forward edge.

	mergeList(edge1, edge2, listReadsForward, listOverlapOffsetsForward, listOfSubstitutionPoses); // Merge the lists from the two edges.

	Edge * newEdge = new Edge(read1,read2, edge1->getOverlapOffset() + edge2->getOverlapOffset(), edge1->getNumOfSubstitutions() + edge2->getNumOfSubstitutions(), listReadsForward, listOverlapOffsetsForward, listOfSubstitutionPoses); // Make the forward edge TODO: number of substitutions and list of substitution positions are not handled correctly now.

	//insertEdge(newEdge);	// Insert the new forward edge in the graph.

	// Delete edges that were merged into the new edge 
	// Q: These were only deleted when there is no flow left in the edges. Does it matter?
	// A: There could be other cases where edges are merged, other than composite edge contraction, so the reads should not be removed in this function)
	//removeEdge(edge1);
	//removeEdge(edge2);
	delete listReadsForward;
	delete listOverlapOffsetsForward;
	delete listOfSubstitutionPoses;

	return newEdge;
	//return true;

}


/**********************************************************************************************************************
  Merge the list of reads, list of overlap offsets, list of substitutions, and list of orientations of two edges.
 **********************************************************************************************************************/
bool OverlapGraph::mergeList(Edge *edge1, Edge *edge2, vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets, vector<UINT64> *listOfSubstitutionPoses)
{
	UINT64 sum = 0;	// Overlap offset sum
	for(UINT64 i = 0; i < edge1->getListOfReads()->size(); i++) 	// Take the list from edge1.
	{
		listReads->push_back(edge1->getListOfReads()->at(i));
		listOverlapOffsets->push_back(edge1->getListOfOverlapOffsets()->at(i));
		sum += edge1->getListOfOverlapOffsets()->at(i);
	}
	/* Update substitution positions */
	for(UINT64 i = 0; i < edge1->getListOfSubstitutionPoses()->size();i++)
		listOfSubstitutionPoses->push_back(edge1->getListOfSubstitutionPoses()->at(i));

	listReads->push_back(edge1->getDestinationRead()->getID()); 	// Insert the common node of the two edges

	/* when merge two edges, they share a common node, no more substitutions will be added, just the ones from the 2 edges */

	listOverlapOffsets->push_back(edge1->getOverlapOffset() - sum);	// Get the overlap offset.

	for(UINT64 i = 0; i < edge2->getListOfReads()->size(); i++)	// take the list from edge2.
	{
		listReads->push_back(edge2->getListOfReads()->at(i));
		listOverlapOffsets->push_back(edge2->getListOfOverlapOffsets()->at(i));
	}
	for(UINT64 i = 0; i < edge2->getListOfSubstitutionPoses()->size();i++)
		listOfSubstitutionPoses->push_back(edge2->getListOfSubstitutionPoses()->at(i));
	return true;
}


/**********************************************************************************************************************
  Update location of reads for the new edge.
  The new edge many contain some reads. We need to update the location information of all such read.
 **********************************************************************************************************************/
bool OverlapGraph::updateReadLocations(Edge *edge)
{
	UINT64 distance = 0;
	for(UINT64 i = 0; i < edge->getListOfReads()->size(); i++) 	// For each read in the edge
	{
		distance += edge->getListOfOverlapOffsets()->at(i);	// Distance of the read in the edge.
		Read *read = dataSet->getReadFromID(edge->getListOfReads()->at(i));	// Get the read object.

		read->getListOfEdgesForward()->push_back(edge);	// Insert the current edge in the list of forward edges.
		read->getLocationOnEdgesForward()->push_back(distance);	// Also insert the distance within the edge.
		read->getListOfEdgesForward()->resize(read->getListOfEdgesForward()->size());
		read->getLocationOnEdgesForward()->resize(read->getLocationOnEdgesForward()->size());
	}
	return true;
}


/**********************************************************************************************************************
  Remove read mapping information of the current edge. This function is called when an edge is removed.
 **********************************************************************************************************************/
bool OverlapGraph::removeReadLocations(Edge *edge)
{
	for(UINT64 i = 0; i < edge->getListOfReads()->size(); i++) 	// For each read in this edge.
	{
		Read *read = dataSet->getReadFromID(edge->getListOfReads()->at(i)); 	// Get the read object.
		for(UINT64 j = 0; j < read->getListOfEdgesForward()->size(); j++)  	// Check the list of edges that contain forward of this read.
		{
			if(read->getListOfEdgesForward()->at(j) == edge)	// Remove the current edge from the list.
			{
				read->getListOfEdgesForward()->at(j) = read->getListOfEdgesForward()->at(read->getListOfEdgesForward()->size()-1);
				read->getLocationOnEdgesForward()->at(j) = read->getLocationOnEdgesForward()->at(read->getLocationOnEdgesForward()->size()-1);

				read->getListOfEdgesForward()->pop_back();
				read->getLocationOnEdgesForward()->pop_back();

				read->getListOfEdgesForward()->resize(read->getListOfEdgesForward()->size());
				read->getLocationOnEdgesForward()->resize(read->getLocationOnEdgesForward()->size());
			}
		}

	}
	return true;
}


/**********************************************************************************************************************
  return edge between source and destination
 **********************************************************************************************************************/
vector<Edge *> OverlapGraph::findEdge(UINT64 source, UINT64 destination)
{
	vector<Edge *> edgesFound;
	for(UINT64 i = 0; i < graph->at(source)->size(); i++) // For the list of edges of the source node.
	{
		if(graph->at(source)->at(i)->getDestinationRead()->getID() == destination)	// check if there is an edge to destination
			edgesFound.push_back(graph->at(source)->at(i));	// return the edge.
	}
	if(edgesFound.empty()){
		FILE_LOG(logINFO) << "Cannot find edge from " << source << " to " << destination;
	}
	return edgesFound;
}


Edge * OverlapGraph::findEdge(UINT64 source, UINT64 destination, UINT64 edgeID)
{
	for(UINT64 i = 0; i < graph->at(source)->size(); i++) // For the list of edges of the source node.
	{
		if(graph->at(source)->at(i)->getDestinationRead()->getID() == destination && graph->at(source)->at(i)->getEdgeID() == edgeID)	// check if there is an edge to destination
			return graph->at(source)->at(i);	// return the edge.
	}
	FILE_LOG(logINFO) << "Check for error " << source << " to " << destination << " with edgeID " << edgeID;
	MYEXIT("Unable to find edge");
}


/**********************************************************************************************************************
  Checks if there is an edge (source, destination)
 **********************************************************************************************************************/
bool OverlapGraph::isEdgePresent(UINT64 source, UINT64 destination)
{
	for(UINT64 i = 0; i < reverseGraph.at(destination).size(); i++){
		if(reverseGraph.at(destination).at(i) == source)
			return true;
	}
//	for(UINT64 i = 0; i < graph->at(source)->size(); i++)	// The list of edges of the source node
//	{
//		if(graph->at(source)->at(i)->getDestinationRead()->getID() == destination)	// check if there is an edge to destination
//			return true;	// return true if there is an edge (source,destination)
//	}
	return false;	// edge not found between source and destination
}


string OverlapGraph::getStringInEdge(Edge *edge, bool includeLast, bool leftClip, bool rightClip)
{
	string source_string, dest_string, readTemp, returnString;
	// strings of the source and destination reads
	source_string =  edge->getSourceRead()->getDnaStringForward();
	dest_string =  edge->getDestinationRead()->getDnaStringForward();

	if(rightClip){
		dest_string = dest_string.substr(0, edge->getDestinationRead()->getReadLength() - edge->getDestinationRead()->getRightClip());
	}
	if(edge->getListOfReads()->empty()){
		returnString = source_string.substr(0,edge->getOverlapOffset());
	}
	else{
		returnString = source_string.substr(0,edge->getListOfOverlapOffsets()->at(0));
		UINT64 lastOffset = edge->getOverlapOffset() - edge->getListOfOverlapOffsets()->at(0);
		UINT64 i  = 0;
		for(i = 0; i < (edge->getListOfReads()->size() - 1); i++)
		{
			readTemp = dataSet->getReadFromID(edge->getListOfReads()->at(i))->getDnaStringForward();
			returnString = returnString + readTemp.substr(0, edge->getListOfOverlapOffsets()->at(i+1));
			lastOffset -= edge->getListOfOverlapOffsets()->at(i+1);
		}
		readTemp = dataSet->getReadFromID(edge->getListOfReads()->at(i))->getDnaStringForward();
		returnString = returnString + readTemp.substr(0, lastOffset);
	}
	if(leftClip){
		UINT64 left_clip = edge->getSourceRead()->getLeftClip();
		if (returnString.length() < left_clip)
			returnString = "";
		else
			returnString = returnString.substr(edge->getSourceRead()->getLeftClip(), string::npos);
	}
	return includeLast?(returnString+dest_string):returnString;
}


/**********************************************************************************************************************
  removes composite path. TODO: any other operation needed? remove transitive edges?? should already be removed
 **********************************************************************************************************************/
bool OverlapGraph::simplifyGraph(void)
{
	//cout << "call simplify graph function" << endl;
	UINT64 counter = 0;
	do
	{
		counter = contractCompositePaths();	// Contract composite paths
	} while (counter > 0);
	return true;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  DFS
 *  Description:  Depth first search of the graph, do topological sort of the nodes that have incident edges. 
 *  Also find back edges that make the graph cyclic, and delete them.
 *  Note: Discover time and finish time is not really needed, can be removed if desired.
 * =====================================================================================
 */
bool OverlapGraph::DFS(vector< UINT64 > * topoSortedNodes) /* Depth first search of the graph, can also be used to determine when the graph is cyclic */
{
	CLOCKSTART;
	vector<Edge *> * backEdges = new vector<Edge*>;
	/* Topological sort of the nodes */
	topoSortedNodes->clear();
	/* Color the node differently to represent the status of the DFS at the nodes */
	vector<int> * searchStatus = new vector<int>;
	/* Predecessors for each node */
	vector<UINT64> *predecessors = new vector<UINT64>;
	/* Discover and finish time for each node in DFS */
	vector<int> * discoverTime = new vector<int>;
	vector<int> * finishTime = new vector<int>;

	UINT64 sTime = 0; /* Searching time */
	for(UINT64 i = 0; i <= dataSet->getNumberOfUniqueReads(); i++){ /* initialize the record vectors */
		searchStatus->push_back(0);
		predecessors->push_back(0);
		discoverTime->push_back(0);
		finishTime->push_back(0);
	}
	for(UINT64 i = 1; i <= dataSet->getNumberOfUniqueReads(); i++){
		if (searchStatus->at(i) == 0)
		{
			DFS_visit(i, searchStatus, predecessors, sTime, discoverTime, finishTime, backEdges, topoSortedNodes);
		}
	}
	/* reverse the topologically sorted list (vector) */
	reverse(topoSortedNodes->begin(), topoSortedNodes->end()); 

	/* Debug, print the topological sorted nodes, and check the result */
	if (loglevel > 5){
		FILE_LOG(logDEBUG2) << "Nodes in topological sorted order:";
		for(UINT64 k = 0; k < topoSortedNodes->size(); k++)
		{
			cout << topoSortedNodes->at(k) << " ";
		}
		cout << endl;
	}
	/* End Debug */

	/* Remove back edges and therefore cycles in the graph */
	if (backEdges->size() > 0){
		FILE_LOG(logINFO) << "Number of back edges found in DFS: " << backEdges->size();
		for(UINT64 k = 0; k < backEdges->size(); k++){
			removeEdge(backEdges->at(k));
		}
	}
	delete backEdges;
	delete searchStatus;
	delete predecessors;
	delete discoverTime;
	delete finishTime;
	CLOCKSTOP;
	return true;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  DFS_visit
 *  Description:  visit all the neighbors of a node, and update the DFS vectors
 * =====================================================================================
 */
bool OverlapGraph::DFS_visit(UINT64 u, vector<int> * searchStatus, vector<UINT64> * predecessors, UINT64 & sTime, vector<int> * discoverTime, vector<int> * finishTime, vector<Edge *> * backEdges, vector<UINT64> * topoSortedNodes)
{
	sTime = sTime + 1;                      /* node i has just been discovered */
	discoverTime->at(u) = sTime;            /* set the discover time */
	searchStatus->at(u) = 1;             /* First time a node is discovered, color it 1 */
	if (graph->at(u)->size() >0)
	{
		for(UINT64 j = 0; j < graph->at(u)->size(); j++){ /* explore all the edges incident to this node */
			UINT64 v = graph->at(u)->at(j)->getDestinationRead()->getID(); /* read ID for the destination read */
			if (searchStatus->at(v) == 0)
			{
				predecessors->at(v) = u; /* v's predecessor is u */
				DFS_visit(v, searchStatus, predecessors, sTime, discoverTime, finishTime, backEdges, topoSortedNodes); /* visit v's neighbors */
			}
			else if (searchStatus->at(v) == 1 ) /* back edge */
			{
				backEdges->push_back(graph->at(u)->at(j));
			}
		}
	}
	searchStatus->at(u) = 2;
	finishTime->at(u) = ++sTime;
	if(graph->at(u)->size()>0 || reverseGraph.at(u).size() > 0)
		topoSortedNodes->push_back(u);          /* insert finished node into the topologically sorted list, only for the non-isolated nodes */
	return true;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  FindLongestPath
 *  Description:  Find longest path in the graph
 * =====================================================================================
 */
bool OverlapGraph::FindLongestPath(vector<UINT64> * topoSortedNodes, string & finalString)         /* Find longest path in the graph */
{
	CLOCKSTART;
	UINT64 numOfNodes = topoSortedNodes->size(); /*  number of nodes that have edges in the graph */
	vector<double> *lengthUntilNodes = new vector<double>; /* maximum length until a node */
	/* bool vector to indicate if the longest path till a node is already calculated */
	vector<bool> *calculated = new vector<bool>;
	/* maximum length of path found so far */
	double maxWeight = 0;

	/* longest path until each node */
	vector<vector<UINT64> *> *longestPathsUntilNodes = new vector<vector<UINT64>*>;
	/* EdgeIDs on the longest path for each node */
	vector<vector<UINT64> *> *edgeIDsUntilNodes = new vector<vector<UINT64>*>;

	UINT64 nodeWithLongestPath; /* node with the longest path found so far */
	/* Initialization of the vectors */
	for(UINT64 i = 0; i <= dataSet->getNumberOfUniqueReads(); i ++){
		lengthUntilNodes->push_back(0); /* Every node has starting weight of 0, so that any can be the start or end of a path */
		vector<UINT64> *longestPathAtnode = new vector<UINT64>;
		vector<UINT64> *edgeIDsAtnode = new vector<UINT64>;
		longestPathsUntilNodes->push_back(longestPathAtnode);
		edgeIDsUntilNodes->push_back(edgeIDsAtnode);
		calculated->push_back(false);
	}

	for(UINT64 i = 0; i < numOfNodes; i++)  /* Initialize the max length until each node */
	{
		UINT64 readID = topoSortedNodes->at(i);
		longestPathsUntilNodes->at(readID)->push_back(readID);
		if(dataSet->getReadFromID(readID)->numInEdges == 0) /* source nodes have initial length 0 */
		{
			lengthUntilNodes->at(readID) = 0;
			calculated->at(readID) = true;
		}
	}
	for(UINT64 i = 0; i < numOfNodes; i++){ /* update max length until each node in topological sorted order */
		UINT64 readID = topoSortedNodes->at(i);
		if(!calculated->at(readID))
		{
			UINT64 maxSource = readID;
			UINT64 maxEdgeID;
			for(UINT64 j = 0; j < reverseGraph.at(readID).size(); j++)
			{
				UINT64 source = reverseGraph.at(readID).at(j); /* source of the edge going into readID */
				vector<Edge *> connectEdges = findEdge(source, readID);
				for(size_t i = 0; i < connectEdges.size(); i++)
				{
					double lengthFromSource = dataSet->getWeight(connectEdges.at(i)) + lengthUntilNodes->at(source);
					if (lengthFromSource > lengthUntilNodes->at(readID)){
						lengthUntilNodes->at(readID) = lengthFromSource;
						maxSource = source;
						maxEdgeID = connectEdges.at(i)->getEdgeID();
					}
				}
			}
			/* Update the path */
			for(UINT64 k = 0; k < longestPathsUntilNodes->at(maxSource)->size(); k++)
				longestPathsUntilNodes->at(readID)->push_back(longestPathsUntilNodes->at(maxSource)->at(k));
			for(UINT64 k = 0; k < edgeIDsUntilNodes->at(maxSource)->size(); k++)
				edgeIDsUntilNodes->at(readID)->push_back(edgeIDsUntilNodes->at(maxSource)->at(k));
			edgeIDsUntilNodes->at(readID)->push_back(maxEdgeID);
			if (lengthUntilNodes->at(readID) > maxWeight)
			{
				maxWeight = lengthUntilNodes->at(readID);
				nodeWithLongestPath = readID;
			}
			calculated->at(readID) = true;
		}
		
	}
	FILE_LOG(logINFO) << "node where the longest path ends is " << nodeWithLongestPath;
	vector<UINT64>::reverse_iterator rit = longestPathsUntilNodes->at(nodeWithLongestPath)->rbegin();
	vector<UINT64>::iterator it = edgeIDsUntilNodes->at(nodeWithLongestPath)->begin();
	UINT64 beginNodeID = *rit;
	//Read * beginR = dataSet->getReadFromID(beginNodeID);
	//FILE_LOG(logINFO) << "\nFirst read in the path has length " << beginR->getReadLength() <<  " bps and " << beginR->getNumOfSubstitutionsInRead() << " subs.";
	if(loglevel > 2){
		for(rit=longestPathsUntilNodes->at(nodeWithLongestPath)->rbegin(); rit!=(longestPathsUntilNodes->at(nodeWithLongestPath)->rend()-1);rit++,it++)
			cout << *rit << " (--" << *it << ", " << findEdge(*rit, *(rit+1), *(it))->getOverlapOffset() << "->) ";
		cout << *rit;
		cout << endl;
	}
	FILE_LOG(logINFO) << "The longest path has weight " << maxWeight;

	/* Print the longest path */
	finalString = "";
	it = edgeIDsUntilNodes->at(nodeWithLongestPath)->begin();
	bool first_edge = true;
	for(rit=longestPathsUntilNodes->at(nodeWithLongestPath)->rbegin(); rit!=(longestPathsUntilNodes->at(nodeWithLongestPath)->rend()-2);rit++,it++)
	{
		Edge * e = findEdge(*rit, *(rit+1), *(it));
		if(first_edge)
		{
			finalString += getStringInEdge(e, false, true, false);
			first_edge = false;
		}
		else
			finalString += getStringInEdge(e, false, false, false);

	}
	/* Last edge in the path need to include the last read string */
	Edge *e = findEdge(*rit, *(rit+1), *(it));
	finalString += getStringInEdge(e, true, false, true);

	/* Remove the edges covered by the longest path */
	removeEdgesBetweenReadNumbers(beginNodeID, *(rit+1));

	delete lengthUntilNodes;
	delete calculated;
	for(UINT64 i = 0; i < longestPathsUntilNodes->size(); i++){
		delete longestPathsUntilNodes->at(i);
	}
	delete longestPathsUntilNodes;
	for(UINT64 i = 0; i < edgeIDsUntilNodes->size(); i++){
		delete edgeIDsUntilNodes->at(i);
	}
	delete edgeIDsUntilNodes;

	CLOCKSTOP;
	return true;
}

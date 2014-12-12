/*
 * OverlapGraph.cpp
 *
 * Created on: Tue Nov 11 12:45:42 EST 2014
 * Author: JJ Chai
 */

#include "OverlapGraph.h"
#include "CS2/cs2.h"

/**************************************************
 * Function to compare two edges. Used for sorting.
 * Edges are sorted by the destination read number.
 **************************************************/
bool compareEdgeByDestination (Edge *edge1, Edge* edge2)
{
	return (edge1->getDestinationRead()->getID() < edge2->getDestinationRead()->getID());
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
	// Initialize the variables.
	minimumOverlapLength = 0;
	maxError = 0;
	maxErrorRate = 0.0;
	numberOfNodes = 0;
	numberOfEdges = 0;
	rubberPos = 0;
}


/**********************************************************************************************************************
  Another Constructor. Build the overlap grpah from the data_Set, and specified parameter values
 **********************************************************************************************************************/
OverlapGraph::OverlapGraph(Dataset *data_Set, const UINT64 & minOverlap, const UINT32 & max_Error, const float & max_ErrorRate, const INT32 & rubber_pos)
{
	// Initialize the variables.
	minimumOverlapLength = minOverlap;
	maxError = max_Error;
	maxErrorRate = max_ErrorRate;
	numberOfNodes = 0;
	numberOfEdges = 0;
	rubberPos = rubber_pos;
	//rubberPos = 10;	// Need to put this in the argument list TODO
	buildOverlapGraphFromDataSet(data_Set);
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
bool OverlapGraph::buildOverlapGraphFromDataSet(Dataset *data_Set)
{
	CLOCKSTART;
	dataSet = data_Set;	// Corresponding data set for this overlap graph
	numberOfNodes = 0;
	numberOfEdges = 0;
	UINT64 counter = 0;	// Number of reads explored so far (maybe)

	// initialize and reserve space for the type vectors
	// They all have size 1 bigger than the number of UniqueReads, why??
	vector<nodeType> *exploredReads = new vector<nodeType>;
	exploredReads->reserve(dataSet->getNumberOfUniqueReads()+1);

	vector<UINT64> * queue = new vector<UINT64>;
	queue->reserve(dataSet->getNumberOfUniqueReads()+1);

	vector<markType> *markedNodes = new vector<markType>;
	markedNodes->reserve(dataSet->getNumberOfUniqueReads()+1);

	graph = new vector< vector<Edge *> * >;
	graph->reserve(dataSet->getNumberOfUniqueReads()+1);

	for(UINT64 i = 0; i <= dataSet->getNumberOfUniqueReads(); i++) // Initialization, one edge list for each unique read
	{
		vector<Edge *> *newList = new vector<Edge *>;	// Vector of edges
		graph->push_back(newList);	// Insert this vector (list) of edges into the graph
		exploredReads->push_back(UNEXPLORED);	// Initialize all reads to be UNEXPLORED
		queue->push_back(0);
		markedNodes->push_back(VACANT);
	}

	UINT64 total_overlap = getAllOverlaps();
	for(UINT64 i = 1; i <= dataSet->getNumberOfUniqueReads(); i++)	// i starts at 1, why?? (exploredReads' first entry for what?)
	{
		if(exploredReads->at(i) == UNEXPLORED && !dataSet->getReadFromID(i)->isContainedRead())
		{
			UINT64 start = 0, end = 0; 	// Initialize queue start and end.
			queue->at(end++) = i;
			while(start < end) 	// This loop will explore all connected component starting from read i.
			{
				counter++;	// Number of reads explored so far
				UINT64 read1 = queue->at(start++);
				if(exploredReads->at(read1) == UNEXPLORED)
				{
					insertAllEdgesOfRead(read1, exploredReads);	// Explore current node, insert all its adjacent edges
					exploredReads->at(read1) = EXPLORED;	// This read is marked as explored
				}
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
					FILE_LOG(logDEBUG)<<"counter: " << setw(10) << counter << " Nodes: " << setw(10) << numberOfNodes << " Edges: " << setw(10) << numberOfEdges/2;
				}
			}
		}
	}

	// report total number
	FILE_LOG(logINFO)<<"counter: " << setw(10) << counter << " Nodes: " << setw(10) << numberOfNodes << " Edges: " << setw(10) << numberOfEdges;

	delete exploredReads;
	delete queue;
	delete markedNodes;

	counter = 0;
	/* contracting composite edges, do not remove dead end nodes for now */
	if (numberOfEdges > 0)
	{
		int iteration = 0;
		do
		{
			if (loglevel > 4 )
			{
				FILE_LOG(logDEBUG1) << iteration << " iteration, after contracting " << counter <<  " composite edges";
				stringstream ss;
				ss << iteration;
				this->printGraph("start_" + ss.str() + ".gdl","start_" + ss.str() + ".fasta");	// printGraph function needs to be rewritten 
				iteration += 1;
			}
			counter = contractCompositePaths();	// need to rewrite contractCompositePaths function
		} while (counter > 0);
	}
	CLOCKSTOP;
	return true;
}


/**********************************************************************************************************************
  Mark all the transitive edges of a read.
  For Details: E.W. Myers. The fragment assembly string graph. Bioinformatics, 21(suppl 2):ii79-ii85, 2005.
 **********************************************************************************************************************/
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
	return true;
}


/**********************************************************************************************************************
  Insert an edge in the graph.
 **********************************************************************************************************************/
bool OverlapGraph::insertEdge(Read *read1, Read *read2, UINT16 overlapOffset)
{
	FILE_LOG(logDEBUG4) << " Insert edge from " << read1->getID() << " to " << read2->getID() << " with overlap offset " << overlapOffset;
	Edge * edge1 = new Edge(read1,read2,overlapOffset);	// Create a new edge in the graph to insert.
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
				UINT64 read = graph->at(index)->at(i)->getDestinationRead()->getID();
				if (graph->at(read)->size() == 1)	// If the destination read only has 1 edge adjacent to it, there are a pair of composite edges
				{
					mergeEdges(graph->at(index)->at(i), graph->at(read)->at(0));	// merge the two edges, need to check
					removeEdge(graph->at(index)->at(i));	// delete the two edges that were contracted
					removeEdge(graph->at(read)->at(0));
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
	
	//for(UINT64 index = 0; index < graph->size(); index++)	// Go through all the nodes in the graph
	//{
	//	if(graph->at(index)->size() > 0)	// If this read has some edges adjacent to it
	//	{
	//		for(UINT64 i = 0; i < graph->at(index)->size(); i++)	// Check all the edges 
	//		{
	//			UINT64 read = graph->at(index)->at(i)->getDestinationRead()->getID();
	//			if (graph->at(read)->size() == 1)	// If the destination read only has 1 edge adjacent to it, there are a pair of composite edges
	//			{
	//				mergeEdges(graph->at(index)->at(i), graph->at(read)->at(0));	// merge the two edges, need to check
	//				removeEdge(graph->at(index)->at(i));	// delete the two edges that were contracted
	//				removeEdge(graph->at(read)->at(0));
	//				counter++;
	//			}
	//		}
	//	}
	//}
	FILE_LOG(logDEBUG) << setw(10) << counter << " composite Edges merged.";
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
	return true;
}



/**********************************************************************************************************************
  Remove an all simple edge in the overlap graph that does not have any flow.
  Definition of simple edge: just a simple overlap, no tiling or intermediate reads in the edge.
 **********************************************************************************************************************/
UINT64 OverlapGraph::removeAllSimpleEdgesWithoutFlow()
{
	CLOCKSTART;
	vector <Edge *> listOfEdges;
	for(UINT64 i = 1; i < graph->size(); i++) // For each read.
	{
		if(!graph->at(i)->empty())	// If the read has some edges.
		{
			for(UINT64 j=0; j < graph->at(i)->size(); j++) // For each edge
			{
				Edge * edge = graph->at(i)->at(j);
				// The edge is simple edge with no flow, and has string shorter than deadEndBp. JJ
				//if(edge->getSourceRead()->getReadNumber() < edge->getDestinationRead()->getReadNumber() && edge->getListOfReads()->empty() && edge->flow == 0 && getStringInEdge(edge).size() < deadEndBp ) 
				if(edge->getSourceRead()->getID() < edge->getDestinationRead()->getID() && edge->getListOfReads()->empty() && edge->flow == 0 ) 
				{
					listOfEdges.push_back(edge); // Put in the list of edges to be removed.
					FILE_LOG(logDEBUG1)  << "removing edge ("<< edge->getSourceRead()->getID()<<","  << edge->getDestinationRead()->getID()<<") OverlapOffset : " << setw(10) << edge->getOverlapOffset(); 
				}
			}
		}
	}
	for(UINT64 i = 0 ; i < listOfEdges.size(); i++)
		removeEdge(listOfEdges.at(i));		// remove the edges from the list.
	FILE_LOG(logDEBUG) << "Simple edges without flow removed: " << listOfEdges.size();
	CLOCKSTOP;
	return listOfEdges.size();
}


/**********************************************************************************************************************
 * This function prints the overlap graph in overlap_graph->gdl file. 
 * The graph can be viewed by aisee (free software available at http://www.aisee.com/) 
 * It also stores the contigs in a file.
 **********************************************************************************************************************/
bool OverlapGraph::printGraph(string graphFileName, string contigFileName)
{
	CLOCKSTART;
	vector <Edge *> contigEdges;	// vector of edges, each of which will form a contig
	vector <Read *> isolateReads;
	UINT64 thickness, highestDegree = 0, highestDegreeNode;
	ofstream graphFilePointer; 

	/************************* Store the graph in a file. ************************/
	graphFilePointer.open(graphFileName.c_str());
	if(!graphFilePointer.is_open())
		MYEXIT("Unable to open file: "+graphFileName);

	// Graph specification before the nodes and edges
	graphFilePointer << "graph: {" << endl <<  "layoutalgorithm :forcedir" << endl <<  "fdmax:704" << endl <<  "tempmax:254" << endl <<  "tempmin:0" << endl <<  "temptreshold:3" << endl <<  "tempscheme:3" << endl <<  "tempfactor:1.08" << endl <<  "randomfactor:100" << endl <<  "gravity:0.0" << endl <<  "repulsion:161" << endl <<  "attraction:43" << endl <<  "ignore_singles:yes" << endl <<  "node.fontname:\"helvB10\"" << endl << "edge.fontname:\"helvB10\"" << endl <<  "node.shape:box" << endl <<  "node.width:80" << endl <<  "node.height:20" << endl <<  "node.borderwidth:1" << endl <<  "node.bordercolor:31" << endl;

	// All the nodes, title and label
	for(UINT64 i = 1; i<= dataSet->getNumberOfUniqueReads(); i++)
	{
		if(dataSet->getReadFromID(i)->superReadID ==0)
			graphFilePointer << "node: { title:\""<< i <<"\" label: \"" << i << "\" }" << endl;	// Print nodes even if there are no edge connected to it
	}

	// All the edges
	for(UINT64 i = 1; i<= dataSet->getNumberOfUniqueReads(); i++)
	{
		if(dataSet->getReadFromID(i)->superReadID ==0 && !graph->at(i)->empty()) // if this read has some edge(s) going out of it (since now the graph is directed)
		{
			// check the degree of the node
			if(graph->at(i)->size() > highestDegree)
			{
				highestDegree = graph->at(i)->size();
				highestDegreeNode = i;

			}
			for(UINT64 j=0; j < graph->at(i)->size(); j++)
			{
				Edge * e = graph->at(i)->at(j);
				UINT64 source = e->getSourceRead()->getID(), destination = e->getDestinationRead()->getID();
				if( e->getDestinationRead()->superReadID==0 && source < destination )	// Print edge only if source node number is smaller than destination node number
				{
					contigEdges.push_back(e); // List of contigs.
					thickness = e->getListOfReads()->empty() ? 1 : 3;	// Thicker edges if composite
					// Now there is only 1 kind of edge since the graph is directed, color doesn't matter either.
					graphFilePointer << "edge: { source:\"" << source << "\" target:\"" << destination << "\" thickness: " << thickness << " arrowstyle: solid backarrowstyle: none color: red label: \"(" <<  e->coverageDepth << "x," << e->getOverlapOffset() << "," << e->getListOfReads()->size() << ")\" }" << endl;
				}
			}
		}
	}
	graphFilePointer << "}";
	graphFilePointer.close();
	FILE_LOG(logINFO) << "Aisee graph written." << endl;
	/************************* Store the graph in a file done. ************************/


	/************************* Store the contigs in a file. ************************/
	if (contigEdges.size() > 0) // Sort the contigs by their length if there are edges in the graph.
	{
		sort(contigEdges.begin(),contigEdges.end(),compareEdgeByOverlapOffset);	// But this is sorted by overlap offset..?? (in a sense, if the offset is longer, the spelled edge is longer, given that the read lengths are more or less comparable.)
		reverse(contigEdges.begin(), contigEdges.end());	// Reverse the order of edges in contigEdges, so that the edges are ordered increasingly by length
	}
	UINT64 sum = 0;	// sum of contig lengths
	ofstream contigFilePointer;
	contigFilePointer.open(contigFileName.c_str());
	if(!contigFilePointer.is_open())
		MYEXIT("Unable to open file: " + contigFileName);
	// All contigs read from the edges
	if (contigEdges.size() > 0)
	{
		for(UINT64 i = 0; i < contigEdges.size(); i++) 
		{
			seqan::DnaString s = getStringInEdge(contigEdges.at(i)); // get the string in the edge. This function need to be rewritten too.
			contigFilePointer << ">contig_"<< i+1 << " Edge  (" << contigEdges.at(i)->getSourceRead()->getID() << ", " << contigEdges.at(i)->getDestinationRead()->getID() << ") String Length: " << seqan::length(s) << " Coverage: " << contigEdges.at(i)->coverageDepth << endl;
			sum += seqan::length(s);
			//contigFilePointer << s << endl;  // save 100 BP in each line.
			UINT32 start=0;
			do
			{
				contigFilePointer << seqan::infix(s,start,min(start+100, seqan::length(s))) << endl;  // save 100 BP in each line.
				start+=100;
			} while (start < seqan::length(s));
		}
	}
	contigFilePointer.close();
	/************************* Store the contigs in a file, done. ************************/

	// Print some statistics about the graph.
	FILE_LOG(logINFO)<< "Total contig length: " << sum <<" bps";
	FILE_LOG(logINFO)<< "Number of Nodes in the graph: " << getNumberOfNodes();
	FILE_LOG(logINFO)<< "Number of Edges in the graph: " << getNumberOfEdges();

	if (contigEdges.size() > 0)
	{
		UINT64 simEdges = 0, comEdges = 0;
		for(UINT64 i=0; i < graph->at(highestDegreeNode)->size(); i++)
		{
			if(graph->at(highestDegreeNode)->at(i)->getListOfReads()->empty())
				simEdges++;
			else
				comEdges++;
		}
		// Print some more statistics on the node with highest degree.
		FILE_LOG(logINFO)<< "Highest Degree Read " << highestDegreeNode << " has " << highestDegree << " neighbors.";
		FILE_LOG(logDEBUG2)<< "Read string: " << dataSet->getReadFromID(highestDegreeNode)->getDnaStringForward();
		FILE_LOG(logINFO) << "It has " << simEdges << " simple edges," << "  and " << comEdges << " composite edges."; 
	}

	CLOCKSTOP;
	return true;
}



/**********************************************************************************************************************
 * Find all the overlaps between reads in the graph, and mark the contained reads in the process.
 **********************************************************************************************************************/
UINT64 OverlapGraph::getAllOverlaps(void)
{
	UINT64 total_overlap = 0;
	for (UINT64 readNumber = 1; readNumber <= dataSet->getNumberOfUniqueReads(); readNumber++)
	{
		FILE_LOG(logDEBUG3) << "find overlap reads of read " << readNumber;
		Read *read1 = dataSet->getReadFromID(readNumber); 	// Get the current read read1.
		INT64 endCoord = read1->getEndCoord() - minimumOverlapLength - rubberPos;
		FILE_LOG(logDEBUG3) << "End coordinate is " << endCoord;
		seqan::DnaString readDnaString = read1->getDnaStringForward(); 	// Get the forward string of read1.
		for(UINT64 j = readNumber + 1 ; j <= dataSet->getNumberOfUniqueReads(); j++) // For each read with start coordinate bigger than the current read, see if they overlap
		{
			if (dataSet->getReadFromID(j)->getStartCoord() > endCoord)	// Too far away to overlap with the current read
				break;
			// Otherwise, realign the two reads and see if they overlap
			Read *read2 = dataSet->getReadFromID(j);
			INT16 overlapOffset;
			if (read1->superReadID == 0 && read2->superReadID == 0 && DoReadsOverlap(read1, read2, overlapOffset))	// If the two reads do overlap, and both non-contained.
			{
				if (overlapOffset > 0 )
				{
					FILE_LOG(logDEBUG4)<< "read" << readNumber << " and read" << j << " overlap";
				}
				else if (overlapOffset < 0 )
				{
					FILE_LOG(logDEBUG4)<< "read" << j << " and read" << readNumber << " overlap";
				}
				total_overlap++;
			}
			else if (read1->superReadID != 0)
			{
				FILE_LOG(logDEBUG4)<< "read" << readNumber << " is contained in read" << read1->superReadID;
				break;
			}
			else if (read2->superReadID != 0)
				FILE_LOG(logDEBUG4)<< "read" << j << " is contained in read" << read2->superReadID;
			else
				FILE_LOG(logDEBUG4)<< "read" << j << "  and read" << readNumber  << "\'s Overlap is not valid!";

		}
	}
	FILE_LOG(logDEBUG)<< "Total number of overlaps found is: " << total_overlap;
	return total_overlap;
}


/**********************************************************************************************************************
 * Insert all edges of a read in the overlap graph.  
 * If a read is already explored, it won't be explored against for another read again.
 **********************************************************************************************************************/
bool OverlapGraph::insertAllEdgesOfRead(UINT64 readNumber, vector<nodeType> * exploredReads)
{
	FILE_LOG(logDEBUG4) << "Insert all edges of read " << readNumber;
	Read *read1 = dataSet->getReadFromID(readNumber); 	// Get the current read read1.
	for(size_t i = 0; i < read1->getOverlapReadIDs()->size(); i++)
	{
		Read *r1 = dataSet->getReadFromID(read1->getOverlapReadIDs()->at(i));
		if(!r1->isContainedRead())
			insertEdge(read1, r1, read1->getOverlapReadOffsets()->at(i));
	}
	if(graph->at(readNumber)->size() != 0)	// Sort the list of edges of the current node according to the overlap offset (ascending) if there are multiple edges
		sort(graph->at(readNumber)->begin(),graph->at(readNumber)->end(), compareEdgeByOverlapOffset); 
	return true;
}


/**********************************************************************************************************************
 * Check if two reads overlap with each other, and if so, find the overlapOffset between them
 **********************************************************************************************************************/
bool OverlapGraph::DoReadsOverlap(Read * read1, Read * read2, INT16 & OverlapOffset)
{
	if (read1->superReadID != 0 || read2->superReadID != 0)	// Do not consider overlap for contained reads
	{
		FILE_LOG(logDEBUG4)<< "One of the reads is contained in other read(s), do not consider overlap between them";
		return false;
	}
		
	UINT64 readNumber1 = read1->getID(), readNumber2 = read2->getID();
	typedef seqan::DnaString TDna;	// Sequence type is DNA sequence
	typedef seqan::Align<TDna, seqan::ArrayGaps> TAlign;	// Align type
	typedef seqan::Row<TAlign>::Type TRow;	// gapped sequence type (GAPS class)
	typedef seqan::Iterator<TRow>::Type TRowIterator;	// Iterator over gapped sequence

	TDna seq1 = read1->getDnaStringForward();	// Get the DnaString of the reads
	TDna seq2 = read2->getDnaStringForward();

	TAlign align;	// Align object of type TAlign
	seqan::resize(rows(align),2);	// Pairwise alignment, so the size is set to 2
	seqan::assignSource(row(align,0), seq1);	// Put the sequences in the align object
	seqan::assignSource(row(align,1), seq2);

	int score = seqan::globalAlignment(align, seqan::Score<int, seqan::Simple>(1,-1,-1), seqan::AlignConfig<true, true, true, true>());	// Get the best score, for overlap alignment, with AlignConfig<t,t,t,t>()

	/* Find opening and ending gaps in each sequence of the alignment */
	TRow row1 = row(align,0);
	TRow row2 = row(align,1);
	TRowIterator it1 = seqan::begin(row1);
	TRowIterator it1End = seqan::end(row1);
	TRowIterator it2 = seqan::begin(row2);
	TRowIterator it2End = seqan::end(row2);

	size_t openGap1 = 0, endGap1 = 0, openGap2 = 0, endGap2 = 0;	
	for(TRowIterator it = it1; it != it1End; ++it)
	{
		if(seqan::isGap(it))
		{
			openGap1++;
			//cout << *it;
		}
		else
			break;
	}
	for(TRowIterator it = --it1End; it != it1; --it)
	{
		if(seqan::isGap(it))
		{
			endGap1++;
			//cout << *it;
		}
		else
			break;
	}
	for(TRowIterator it = it2; it != it2End; ++it)
	{
		if(seqan::isGap(it))
		{
			openGap2++;
			//cout << *it;
		}
		else
			break;
	}
	for(TRowIterator it = --it2End; it != it2; --it)
	{
		if(seqan::isGap(it))
		{
			endGap2++;
			//cout << *it;
		}
		else
			break;
	}

	size_t alignLength = seqan::length(row1);	// Number of base pairs in the alignment, could also use length(row2)
	OverlapOffset = openGap2 - openGap1;	// Overlapoffset, if negative, second sequence has left overhang
	UINT32 misMatches = (alignLength - openGap1 - openGap2 - endGap1 - endGap2 - score)/2;	// Number of mismatches, including gaps
	UINT32 overlapLength = alignLength - openGap1 - openGap2 - endGap1 - endGap2;
	//int ref_overlapOffset = read2->getStartCoord() - read1->getStartCoord();

	FILE_LOG(logDEBUG4) << "Total align length: " << alignLength << "; OverlapOffset: " << OverlapOffset << "; misMatches: " << misMatches << "; overlapLength: " << overlapLength;
	/* Check and see if one read is contained in the other, they are not considered as overlapping in this case */
	if (misMatches <= maxError || misMatches <= overlapLength * maxErrorRate) 	// Only if two sequences are similar enough
	{
		/* Change the starting coordinate of the second read, 
		 * depending on the overlap offset with the first read, whose coordinate is already set
		 * DO NOT DO THIS NOW, SINCE STARTING COORDINATES ARE NOT USED ANYWHERE
		 */
		//if (read2->getStartCoord() != (read1->getStartCoord() + OverlapOffset))
		//{
		//	bool change_coord = read2->setStartCoord(read1->getStartCoord() + OverlapOffset);
		//	//FILE_LOG(logDEBUG2) << "Depending on read" << readNumber1 << "\'s coord " << read1->getStartCoord() << ", Change read" << readNumber2 << "\'s starting coordinate to " << read2->getStartCoord();
		//}
		if (misMatches != 0)
		{
			FILE_LOG(logDEBUG2) << "Total align length: " << alignLength << "; OverlapOffset: " << OverlapOffset << "; misMatches: " << misMatches << "; overlapLength: " << overlapLength;
			FILE_LOG(logDEBUG4) << align;
		}

		if (openGap1 + endGap1 == 0)	// Seq1 does not contain any gap (contains seq2)
		{
			// Neither of the two reads are contained in other reads, so currently this will be the read's first contained read
			read2->superReadID = read1->getID();	// Assign read2's superReadID to read1's ID
			read1->getContainedReadIDs()->push_back(read2->getID());	// Add read2 to read1's containedReadIDs
			if (!read2->getContainedReadIDs()->empty())	// read2 also contains other reads
			{
				//reads contained in read2 now should also be contained in read1
				for(size_t i = 0; i < read2->getContainedReadIDs()->size(); i++)
				{
					read1->getContainedReadIDs()->push_back(read2->getContainedReadIDs()->at(i));	// Put reads contained in Read2 as Read1's contained reads
					dataSet->getReadFromID(read2->getContainedReadIDs()->at(i))->superReadID = read1->getID();	// Set those contained reads' super read ID to read1's ID
				}
			}
			dataSet->numberOfNonContainedReads = dataSet->numberOfNonContainedReads - 1;
			FILE_LOG(logDEBUG4) << "read" << readNumber1 << " contains read" << readNumber2;
			return false;
		}
		else if (openGap2 + endGap2 == 0)	// Seq2 does not contain any gap (contains seq1)
		{
			read1->superReadID = read2->getID();
			read2->getContainedReadIDs()->push_back(read1->getID());
			if (!read1->getContainedReadIDs()->empty())
			{
				for(size_t i = 0; i < read1->getContainedReadIDs()->size(); i++)
				{
					read2->getContainedReadIDs()->push_back(read1->getContainedReadIDs()->at(i));
					dataSet->getReadFromID(read1->getContainedReadIDs()->at(i))->superReadID = read2->getID();
				}
			}
			//removeEdgesOfRead(read1);	// Remove all the edges adjacent to read2 because it's contained in other reads. TODO: write this function
			dataSet->numberOfNonContainedReads = dataSet->numberOfNonContainedReads - 1;
			FILE_LOG(logDEBUG4) << "Read" << readNumber2 << " contains Read" << readNumber1;
			return false;
		}
		/* If neither is contained in the other, check if the overlap length passes the threshold */
		else if (overlapLength >= minimumOverlapLength && -OverlapOffset < rubberPos)
		//else if (overlapLength >= minimumOverlapLength)
		{
			FILE_LOG(logDEBUG4) << "Found overlap between " << readNumber1 << " and " << readNumber2;
			if (OverlapOffset > 0)
			{
				read1->getOverlapReadIDs()->push_back(readNumber2);
				read1->getOverlapReadOffsets()->push_back(OverlapOffset);
			}
			else
			{
				read2->getOverlapReadIDs()->push_back(readNumber1);
				read2->getOverlapReadOffsets()->push_back(-OverlapOffset);
			}
			return true;
		}
		else if (overlapLength < minimumOverlapLength)
		{
			FILE_LOG(logDEBUG4) << readNumber1 << " and " << readNumber2 << " do not contain each other and do not overlap long enough ";
			return false;
		}
		else
		{
			FILE_LOG(logDEBUG2) << readNumber1 << " and " << readNumber2 << " overlap offset " << OverlapOffset << " does not fit the coordiates " << read1->getStartCoord() << ", " << read2->getStartCoord();
			return false;
		}
	}
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
		}
	}
	graph->at(readNumber)->resize(j);
	if ( (dataSet->getReadFromID(readNumber)->numInEdges + dataSet->getReadFromID(readNumber)->numOutEdges)==0)
		numberOfNodes--;
	return true;
}


/**********************************************************************************************************************
  Remove all edges in and out of a given read.
 **********************************************************************************************************************/
bool OverlapGraph::removeEdgesOfRead(Read * read)
{
	UINT64 readNumber = read->getID();
	for (UINT64 index = 0; index < graph->size(); index++)	// Going through all the edges
	{
		if (graph->at(index)->size() > 0)
		{
			for( UINT64 j = 0; j < graph->at(index)->size(); j++)
			{
				Edge *edge = graph->at(index)->at(j);
				if (edge->getSourceRead()->getID() == readNumber || edge->getDestinationRead()->getID() == readNumber)
					removeEdge(edge);
			}
		}
	}
	return true;
}
/**********************************************************************************************************************
  Merge two edges in the overlap graph.
 **********************************************************************************************************************/
bool OverlapGraph::mergeEdges(Edge *edge1, Edge *edge2)
{
	Edge *newEdge = new Edge();
	Read *read1 = edge1->getSourceRead(), *read2 = edge2->getDestinationRead();

	vector<UINT64> * listReadsForward = new vector<UINT64>;	// List of reads in the forward edge.
	vector<UINT16> * listOverlapOffsetsForward= new vector<UINT16>;	// List of Overlaps in the forward edge.

	mergeList(edge1, edge2, listReadsForward, listOverlapOffsetsForward); // Merge the lists from the two edges.

	newEdge->makeEdge(read1,read2, edge1->getOverlapOffset() + edge2->getOverlapOffset(), listReadsForward, listOverlapOffsetsForward); // Make the forward edge

	insertEdge(newEdge);	// Insert the new forward edge in the graph.

	// Delete edges that were merged into the new edge 
	// Q: These were only deleted when there is no flow left in the edges. Does it matter?
	// A: There could be other cases where edges are merged, other than composite edge contraction, so the reads should not be removed in this function)
	//removeEdge(edge1);
	//removeEdge(edge2);

	return true;

}


/**********************************************************************************************************************
  Merge the list of reads, list of overlap offsets and list of orientations of two edges.
 **********************************************************************************************************************/
bool OverlapGraph::mergeList(Edge *edge1, Edge *edge2, vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets)
{
	UINT64 sum = 0;	// Overlap offset sum
	for(UINT64 i = 0; i < edge1->getListOfReads()->size(); i++) 	// Take the list from edge1.
	{
		listReads->push_back(edge1->getListOfReads()->at(i));
		listOverlapOffsets->push_back(edge1->getListOfOverlapOffsets()->at(i));
		sum += edge1->getListOfOverlapOffsets()->at(i);
	}
	listReads->push_back(edge1->getDestinationRead()->getID()); 	// Insert the common node of the two edges

	listOverlapOffsets->push_back(edge1->getOverlapOffset() - sum);	// Get the overlap offset.

	for(UINT64 i = 0; i < edge2->getListOfReads()->size(); i++)	// take the list from edge2.
	{
		listReads->push_back(edge2->getListOfReads()->at(i));
		listOverlapOffsets->push_back(edge2->getListOfOverlapOffsets()->at(i));
	}
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
  This function calculates the cost and bounds for an edge in the overlap graph.
  This function is very sensitive to the assembled contigs.
 **********************************************************************************************************************/
bool OverlapGraph::calculateBoundAndCost(Edge *edge, INT64* FLOWLB, INT64* FLOWUB, INT64* COST)
{
	for(UINT64 i = 0; i < 3; i++)		// For the simple edges we put very high cost
	{
		FLOWLB[i] = 0; FLOWUB[i] = 10; COST[i] = 500000;
	}

	if(!edge->getListOfReads()->empty()) // Composite Edge
	{
		if(edge->getListOfReads()->size() > 5 || edge->getOverlapOffset() > 1000 ) // Composite Edge of at least 20 reads, or with length at least 1000. Must have at least one unit of flow.
		{
			FLOWLB[0] = 1; FLOWUB[0] = 1; COST[0] = 1;
			FLOWLB[1] = 0; FLOWUB[1] = 1; COST[1] = 50000;
			FLOWLB[2] = 0; FLOWUB[2] = 8; COST[2] = 100000;
		}
		else // Short composite edge containing less than 20 reads or with overlap offset less than 1000. May have zero flow.
		{
			FLOWLB[0] = 0; FLOWUB[0] = 1; COST[0] = 1;
			FLOWLB[1] = 0; FLOWUB[1] = 1; COST[1] = 50000;
			FLOWLB[2] = 0; FLOWUB[2] = 8; COST[2] = 100000;
		}
	}

	return true;
}


/********************************
 * Calculate minimum cost flow
 ********************************/
bool OverlapGraph::calculateFlow(string inputFileName, string outputFileName)
{
	CLOCKSTART;
	// Add super source and super sink nodes, add edge from super sink to super source with very big cost
	// Add edge from super source to every node in the graph, also from every node in the graph to the super sink
	// Every edge will be assigned a lower bound and an upper bound of the flow (capacity), and the cost associated with the edge
	UINT64 V = numberOfNodes + 2, E = numberOfEdges * 3 + numberOfNodes * 2 + 1 , SUPERSOURCE = 1, SUPERSINK = V;
	INT64 FLOWLB[3], FLOWUB[3], COST[3];			// Flow bounds and cost of the edges, cost function originally is a piecewise function with 3 segments
	ofstream outputFile;
	outputFile.open(inputFileName.c_str());
	if(!outputFile.is_open())
		MYEXIT("Unable to open file: "+inputFileName);
	stringstream ss;
	ss << "p min " << setw(10) << V << " " << setw(10) << E << endl;  	// Problem description: Number of nodes and edges in the graph.
	ss << "n " << setw(10) << SUPERSOURCE << setw(10) << " 0" << endl;	// Flow in the super source
	ss << "n " << setw(10) << SUPERSINK << setw(10) << " 0" << endl;	// Flow in the super sink.

	// Flow lower bound and upper bound, and the cost for the first segment in the piecewise cost function 
	FLOWLB[0] = 1; FLOWUB[0] = 1000000; COST[0] = 1000000;
	ss << "a " << setw(10) << SUPERSINK << " " << setw(10) << SUPERSOURCE << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl; // Add an edge from super sink to super source with very high cost (almost infinity), also at most can be used once


	// If the ID of a node in the original graph is 100 and directed graph is 5
	// Then listOfNodes->at(100) is equal to 5
	// and ListOfNodesReverse->at(5) is equal to 100.
	vector<UINT64> *listOfNodes = new vector<UINT64>;
	vector<UINT64> *listOfNodesReverse = new vector<UINT64>;

	for(UINT64 i = 0; i <= graph->size(); i++)		// For n nodes in the graph, CS2 requires that the nodes are numbered from 1 to n. In the overlap graph, the nodes does not have sequencinal ID. We need to convert them to 1 - n
	{
		listOfNodes->push_back(0);
		listOfNodesReverse->push_back(0);
	}

	// This loop set lower bound and upper bound from super source to each node, and from each node to super sink. All costs are 0.
	UINT64 currentIndex = 1;
	for(UINT64 i = 1; i < graph->size(); i++)
	{
		if( (dataSet->getReadFromID(i)->numInEdges + dataSet->getReadFromID(i)->numOutEdges) != 0 ) // edges to and from the super source and super sink
		{
			FILE_LOG(logDEBUG2) << "Found node " << i << " corresponding to index " << currentIndex;
			listOfNodes->at(i) = currentIndex;					// Mapping between original node ID and cs2 node ID
			listOfNodesReverse->at(currentIndex) = i;			// Mapping between original node ID and cs2 node ID
			FLOWLB[0] = 0; FLOWUB[0] = 1000000; COST[0] = 0;
			ss << "a " << setw(10) << SUPERSOURCE << " " << setw(10) << currentIndex + 1 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
			ss << "a " << setw(10) << currentIndex + 1 << " " << setw(10) << SUPERSINK << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
			currentIndex++;
		}
	}

	// This loop converts the original bi-directed edges to directed edges (1 becomes 6).
	// This loop set the lower and upper bounds of the flow in each edge, and the cost
	for(UINT64 i = 1; i < graph->size(); i++)
	{
		if(!graph->at(i)->empty()) // edges to and from the super source and super sink
		{
			for(UINT64 j = 0; j < graph->at(i)->size(); j++)
			{
				Edge *edge = graph->at(i)->at(j);
				UINT64 u = listOfNodes->at(edge->getSourceRead()->getID());	// Node number of source read in the new graph
				UINT64 v = listOfNodes->at(edge->getDestinationRead()->getID());	// Node number of destination read in the new graph

				calculateBoundAndCost(edge, FLOWLB, FLOWUB, COST);	// Calculate bounds and cost depending on the edge property (number of reads in edge, or string length)

				// Here for each edge we add three edges with different values of cost and bounds, 3 pieces in the piecewise function
				ss << "a " << setw(10) << u + 1 << " " << setw(10) << v + 1 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
				ss << "a " << setw(10) << u + 1 << " " << setw(10) << v + 1 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;
				ss << "a " << setw(10) << u + 1 << " " << setw(10) << v + 1 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;
			}
		}
	}
	outputFile << ss.str();		// Write the string in a file for CS2
	outputFile.close();
	ss.str(std::string());

	char * inFile = new char[inputFileName.size() + 1];
	std::copy(inputFileName.begin(), inputFileName.end(), inFile);
	inFile[inputFileName.size()] = '\0';

	char * outFile = new char[outputFileName.size() + 1];
	std::copy(outputFileName.begin(), outputFileName.end(), outFile);
	outFile[outputFileName.size()] = '\0';


	FILE_LOG(logINFO) << "Calling CS2";
	main_cs2(inFile,outFile);			// Call CS2
	FILE_LOG(logINFO) << "CS2 finished";

	delete[] inFile;
	delete[] outFile;

	ifstream inputFile;
	inputFile.open(outputFileName.c_str());
	if(!inputFile.is_open())
		MYEXIT("Unable to open file: "+outputFileName);


	string s, d, f;
	UINT64 lineNum = 0;
	while(!inputFile.eof())
	{
		lineNum ++;
		UINT64 source, destination, flow;
		inputFile >> source >> destination >> flow;		// get the flow from CS2

		if(source != SUPERSINK && source != SUPERSOURCE && destination != SUPERSOURCE && destination != SUPERSINK && flow!=0)
		{
			UINT64 mySource = listOfNodesReverse->at(source-1);				// Map the source to the original graph
			UINT64 myDestination = listOfNodesReverse->at(destination-1);	// Map the destination in the original graph
			Edge *edge = findEdge(mySource, myDestination);	// Find the edge in the original graph.
			edge->flow += flow;												// Add the flow in the original graph.
			FILE_LOG(logDEBUG2) << "Edge from " << mySource << " to " << myDestination << " has flow " << edge->flow;
		}
	}
	inputFile.close();
	delete listOfNodes;
	delete listOfNodesReverse;
	this->flowComputed = true;
	CLOCKSTOP;
	return true;
}


/**********************************************************************************************************************
  return edge between source and destination
 **********************************************************************************************************************/
Edge * OverlapGraph::findEdge(UINT64 source, UINT64 destination)
{
	for(UINT64 i = 0; i < graph->at(source)->size(); i++) // For the list of edges of the source node.
	{
		if(graph->at(source)->at(i)->getDestinationRead()->getID() == destination)	// check if there is an edge to destination
			return graph->at(source)->at(i);	// return the edge.
	}
	FILE_LOG(logINFO) << "Check for error " << source << " to " << destination;
	MYEXIT("Unable to find edge");
}


/**********************************************************************************************************************
  Checks if there is an edge (source, destination)
 **********************************************************************************************************************/
bool OverlapGraph::isEdgePresent(UINT64 source, UINT64 destination)
{
	for(UINT64 i = 0; i < graph->at(source)->size(); i++)	// The list of edges of the source node
	{
		if(graph->at(source)->at(i)->getDestinationRead()->getID() == destination)	// check if there is an edge to destination
			return true;	// return true if there is an edge (source,destination)
	}
	return false;	// edge not found between source and destination
}


/**********************************************************************************************************************
  This function returns the string spelled out by overlapping the reads in an edge in the overlap graph
 **********************************************************************************************************************/
seqan::DnaString OverlapGraph::getStringInEdge(Edge *edge)
{
	seqan::DnaString read1_string, read2_string, readTemp, returnString;
	// strings of the source and destination reads
	read1_string =  edge->getSourceRead()->getDnaStringForward();
	read2_string =  edge->getDestinationRead()->getDnaStringForward();
	returnString = read1_string;

	UINT64 previousLength = seqan::length(read1_string);
	UINT64 substringLength;

	// Going through all the reads on the edge (if this edge is composite)
	for(UINT64 i = 0; i < edge->getListOfReads()->size(); i++)
	{
		readTemp = dataSet->getReadFromID(edge->getListOfReads()->at(i))->getDnaStringForward();

		substringLength =  seqan::length(readTemp) + edge->getListOfOverlapOffsets()->at(i) - previousLength;	// Length of the added substring
		if( edge->getListOfOverlapOffsets()->at(i) ==  previousLength)	// Overlap offset is equal to the length of the previous read (not really overlap, two reads just touch each other at the ends)
			seqan::append(returnString, "N");
		seqan::append(returnString, seqan::infix(readTemp, seqan::length(readTemp) - substringLength, seqan::length(readTemp)) );
		previousLength = seqan::length(readTemp);
	}

	if(edge->getListOfReads()->empty()) // Simple edge
	{
		substringLength =  seqan::length(read2_string) + edge->getOverlapOffset() - seqan::length(read1_string);
		seqan::append(returnString, seqan::infix(read2_string, seqan::length(read2_string) - substringLength, seqan::length(read2_string)));
	}
	else
	{
		substringLength = edge->getListOfOverlapOffsets()->at(0);
		seqan::append(returnString, seqan::infix(read2_string, seqan::length(read2_string) - substringLength, seqan::length(read2_string)));
	}
	return returnString;
}


/**********************************************************************************************************************
  removes composite path. TODO: any other operation needed? remove transitive edges?? should already be removed
 **********************************************************************************************************************/
bool OverlapGraph::simplifyGraph(void)
{
	UINT64 counter = 0;
	do
	{
		counter = contractCompositePaths();	// Contract composite paths
	} while (counter > 0);
	return true;
}


/**********************************************************************************************************************
  Calculate the coverage depth of an edge for every basepair and then update the Mean and SD of coverage depth in
  the edge. Only consider reads that are unique to the edge.
 **********************************************************************************************************************/
void OverlapGraph::getBaseByBaseCoverage(Edge *edge)
{
	vector<UINT64> * coverageBaseByBase = new vector<UINT64>;
	UINT64 length = edge->getOverlapOffset() + edge->getDestinationRead()->getReadLength();	// Array lenght same as the string length in the edge.
	for(UINT64 i = 0; i <=length; i++)
	{
		coverageBaseByBase->push_back(0);	// At first all the bases are covered 0 times.
	}
	UINT64 overlapOffset = 0;
	for(UINT64 i = 0; i < edge->getListOfReads()->size(); i++)	// For each read in the edge.
	{
		Read *read = dataSet->getReadFromID(edge->getListOfReads()->at(i));
		overlapOffset += edge->getListOfOverlapOffsets()->at(i);	// Where the current read starts in the string.
		UINT64 readLength = read->getReadLength();
		for(UINT64 j = overlapOffset; j < overlapOffset + readLength; j++)
		{
			coverageBaseByBase->at(j) += read->getFrequency();	// Increase the coverage of all bases by the frequency of the read.
		}
	}

	overlapOffset = 0;
	for(UINT64 i = 0; i < edge->getListOfReads()->size(); i++) // Scan the reads again.
	{
		Read *read = dataSet->getReadFromID(edge->getListOfReads()->at(i));
		overlapOffset += edge->getListOfOverlapOffsets()->at(i);
		UINT64 readLength = read->getReadLength();
		if(read->getListOfEdgesForward()->size() > 1)	// Clear the bases that are covered by reads apperaing in multiple places in the graph.
		{
			for(UINT64 j = overlapOffset; j < overlapOffset + readLength; j++)
			{
				coverageBaseByBase->at(j) = 0;
			}
		}
	}
	for(UINT64 i = 0; i < edge->getSourceRead()->getReadLength(); i++)	// For the source read, clear the bases. Because this read is present in multiple places or not all reads are considered for these bases.
	{
		coverageBaseByBase->at(i) = 0;
	}

	for(UINT64 i = 0; i < edge->getDestinationRead()->getReadLength(); i++)	// Similarly clear the bases covered by the destination read.
	{
		coverageBaseByBase->at(coverageBaseByBase->size() - 1 - i) = 0;
	}

	UINT64 sum = 0, variance=0, count = 0, mean = 0, sd = 0;

	for(UINT64 i = 0; i < coverageBaseByBase->size(); i++)
	{
		if(coverageBaseByBase->at(i) != 0)	// Count only the non-zero values.
		{
			sum += coverageBaseByBase->at(i);
			count++;
		}
	}
	if ( count != 0 )
	{
		mean = sum/count;	// Calculate the mean.

		for(UINT64 i = 0; i < coverageBaseByBase->size(); i++)
		{
			if(coverageBaseByBase->at(i) != 0)	// Calculate the variance.
			{
				variance += (mean - coverageBaseByBase->at(i)) * (mean - coverageBaseByBase->at(i));
			}
		}
		sd = sqrt(variance/count);	// Calculate the standard deviation.
	}
	edge->coverageDepth = mean;	// Update the mean of the current edge.
	edge->SD = sd;	// Update the standard deviation of the current edge.
	delete coverageBaseByBase;
}


/**********************************************************************************************************************
  For each node in the graph, sort all its incident edges according to destination read ID.
 **********************************************************************************************************************/
void OverlapGraph::sortEdges()
{
	for(UINT64 i = 1; i < graph->size(); i++)
	{
		if(!graph->at(i)->empty())
		{
			sort(graph->at(i)->begin(), graph->at(i)->end(), compareEdgeByDestination);
		}
	}
}


/**********************************************************************************************************************
  This function returns the edit distance between two strings.
  Code downloaded from http://rosettacode.org/wiki/Levenshtein_distance#C.2B.2B
 **********************************************************************************************************************/
UINT64 OverlapGraph::calculateEditDistance(const std::string &s1, const std::string &s2)
{
	const UINT64 m(s1.size());
	const UINT64 n(s2.size());
	if( m==0 )
		return n;
	if( n==0 )
		return m;
	UINT64 *costs = new UINT64[n + 1];
	for( UINT64 k=0; k<=n; k++ )
		costs[k] = k;

	UINT64 i = 0;
	for ( std::string::const_iterator it1 = s1.begin(); it1 != s1.end(); ++it1, ++i )
	{
		costs[0] = i+1;
		UINT64 corner = i;
		UINT64 j = 0;
		for ( std::string::const_iterator it2 = s2.begin(); it2 != s2.end(); ++it2, ++j )
		{
			UINT64 upper = costs[j+1];
			if( *it1 == *it2 )
			{
				costs[j+1] = corner;
			}
			else
			{
				UINT64 t(upper<corner?upper:corner);
				costs[j+1] = (costs[j]<t?costs[j]:t)+1;
			}
			corner = upper;
		}
	}
	UINT64 result = costs[n];
	delete [] costs;
	//cout << s1 << endl << s2 << endl << result<< endl;
	return result;
}

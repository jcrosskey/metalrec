/*
 * OverlapGraph.cpp
 *
 * Created on: Tue Nov 11 12:45:42 EST 2014
 * Author: JJ Chai
 */

#include "OverlapGraph.h"

/**************************************************
 * Function to compare two edges. Used for sorting.
 * Edges are sorted by the destination read number.
**************************************************/
bool compareEdgeByDestination (Edge *edge1, Edge* edge2)
{
	return (edge1->getDestinationRead()->getReadNumber() < edge2->getDestinationRead()->getReadNumber());
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
	numberOfNodes = 0;
	numberOfEdges = 0;
}

/**********************************************************************************************************************
	Another Constructor. Build the overlap grpah from the data_Set
**********************************************************************************************************************/
OverlapGraph::OverlapGraph(DataSet *data_Set, UINT64 minOverlap, UINT32 max_Error, float max_ErrorRate)
{
	// Initialize the variables.
	minimumOverlapLength = minOverlap;
	maxError = max_Error;
	maxErrorRate = max_ErrorRate;
	numberOfNodes = 0;
	numberOfEdges = 0;
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
	Build the overlap graph from data set
**********************************************************************************************************************/
bool OverlapGraph::buildOverlapGraphFromDataSet(DataSet *data_Set)
{
	CLOCKSTART;
	numberOfNodes = 0;
	numberOfEdges = 0;
	dataSet = data_Set;	// Corresponding data set for this overlap graph
	UINT64 counter = 0;
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

	//markContainedReads();	// Do not mark contained reads here, instead, do it at the pairwise alignment stage

	for(UINT64 i = 1; i <= dataSet->getNumberOfUniqueReads(); i++)
	{
		if(exploredReads->at(i) == UNEXPLORED)
		{
			UINT64 start = 0, end = 0; 	// Initialize queue start and end.
			queue->at(end++) = i;
			while(start < end) 	// This loop will explore all connected component starting from read i.
			{
				counter++;	// Number of nodes explored in the graph ??
				UINT64 read1 = queue->at(start++);
				if(exploredReads->at(read1) == UNEXPLORED)
				{
					insertAllEdgesOfRead(read1, exploredReads);	// Explore current node, insert all its adjacent edges
					exploredReads->at(read1) = EXPLORED;	// This read is marked as explored
				}
				if(graph->at(read1)->size() != 0) 	// Read has some edges (required only for the first read when a new queue starts.
				{
					if(exploredReads->at(read1) == EXPLORED) 	// Explore unexplored neighbors first.
					{
						for(UINT64 index1 = 0; index1 < graph->at(read1)->size(); index1++ )
						{
							UINT64 read2 = graph->at(read1)->at(index1)->getDestinationRead()->getReadNumber();
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
					if(exploredReads->at(read1) == EXPLORED_AND_TRANSITIVE_EDGES_MARKED)
					{
						for(UINT64 index1 = 0; index1 < graph->at(read1)->size(); index1++) 	// Then explore all neighbour's neighbors
						{
							UINT64 read2 = graph->at(read1)->at(index1)->getDestinationRead()->getReadNumber();
							if(exploredReads->at(read2) == EXPLORED)
							{
								for(UINT64 index2 = 0; index2 < graph->at(read2)->size(); index2++) 	// Explore all neighbors neighbors
								{
									UINT64 read3 = graph->at(read2)->at(index2)->getDestinationRead()->getReadNumber();
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
						removeTransitiveEdges(read1); // Remove the transitive edges
					}
				}
				if(counter%100000==0)	// Show the progress.
				{
					FILE_LOG(logDEBUG)<<"counter: " << setw(10) << counter << " Nodes: " << setw(10) << numberOfNodes << " Edges: " << setw(10) << numberOfEdges/2;
				}
			}
		}
	}
	FILE_LOG(logINFO)<<"counter: " << setw(10) << counter << " Nodes: " << setw(10) << numberOfNodes << " Edges: " << setw(10) << numberOfEdges/2;

	delete exploredReads;
	delete queue;
	delete markedNodes;
	delete hashTable;	// Do not need the hash table any more.
	//JJ : print graph and contigs before doing composite edge contraction and deadend removal
	if (loglevel > 2)
	{
		FILE_LOG(logINFO) << "Print graph before contracting composite edges and removing dead ends to start_0.*";
		this->printGraph("start_0.gdl","start_0.fasta");
	}
	if (numberOfEdges > 0)
	{
		int iteration = 0;
		do
		{
			if (loglevel > 2 )
			{
				FILE_LOG(logDEBUG1) << iteration << " iteration, after contracting " << counter <<  " composite edges";
				stringstream ss;
				ss << iteration;
				this->printGraph("start_" + ss.str() + ".gdl","start_" + ss.str() + ".fasta");
				iteration += 1;
			}
			counter = contractCompositePaths();
	//counter += removeDeadEndNodes();
		} while (counter > 0);
	}
	CLOCKSTOP;
	return true;
}

/**********************************************************************************************************************
	This function check if a read contains other small reads. If a read is contained in more than one super read
	then it is assigned to the longest such super read.
**********************************************************************************************************************/
//void OverlapGraph::markContainedReads(void)
//{
//	CLOCKSTART;
//	if(dataSet->longestReadLength == dataSet->shortestReadLength) // If all reads are of same length, then no need to do look for contained reads.
//	{
//		FILE_LOG(logINFO) << "All reads are of same length. No contained reads.";
//		CLOCKSTOP;
//		return;
//	}
//	UINT64 counter = 0;
//	for(UINT64 i = 1; i < graph->size(); i++) // For each read
//	{
//		Read *read1 = dataSet->getReadFromID(i); // Get the read
//		seqan::DnaString readString = read1->getDnaStringForward(); // Get the forward string of the read
//		string subString;
//		for(UINT64 j = 1; j < read1->getReadLength() - hashTable->getHashStringLength(); j++) // For each substring of read1 of length getHashStringLength
//		{
//			subString = readString.substr(j,hashTable->getHashStringLength()); // Get the substring from read1
//			vector<UINT64> * listOfReads=hashTable->getListOfReads(subString); // Search the substring in the hash table
//			if(!listOfReads->empty()) // If other reads contain the substring as prefix or suffix
//			{
//				for(UINT64 k = 0; k < listOfReads->size(); k++) // For each read in the list.
//				{
//					UINT64 data = listOfReads->at(k); // We used bit operation in the hash table to store read ID and orientation
//					Read *read2 = dataSet->getReadFromID(data & 0X3FFFFFFFFFFFFFFF); 	// Least significant 62 bits store the read number.
//	// Most significant 2 bits store the orientation.
//	// Orientation 0 means prefix of forward of the read
//	// Orientation 1 means suffix of forward of the read
//	// Orientation 2 means prefix of reverse of the read
//	// Orientation 3 means prefix of reverse of the read
//
//					if(readString.length() > read2->getStringForward().length() && checkOverlapForContainedRead(read1,read2,(data >> 62),j)) // read1 need to be longer than read2 in order to contain read2
//																																			 // Check if the remaining of the strings also match
//					{
//						if(read2->superReadID == 0) // This is the rist super read found. we store the ID of the super read.
//						{
//							read2->superReadID = i;
//							counter ++;
//						}
//						else // Already found some previous super read. Now we have another super read.
//						{
//							if(readString.length() > dataSet->getReadFromID(read2->superReadID)->getReadLength()) // This super read is longer than the previous super read. Update the super read ID.
//								read2->superReadID = i;
//						}
//					}
//				}
//			}
//		}
//		if(i%1000000 == 0)
//		{
//			FILE_LOG(logDEBUG) << setw(10) << counter << " contained reads in " << setw(10) << i << " super reads.";
//		}
//	}
//
//	// Get some statistics
//	UINT64 containedReads = 0, nonContainedReads = 0;
//	for(UINT64 i = 1; i < graph->size(); i++)
//	{
//		Read *rr = dataSet->getReadFromID(i);
//		if(rr->superReadID == 0) // Count the number of reads that are not contained by some other reads.
//			nonContainedReads++;
//		else	// Count the number of reads that are contained by some other read.
//			containedReads++;
//	}
//	FILE_LOG(logINFO) << setw(10) << nonContainedReads << " Non-contained reads. (Keep as is)";
//	FILE_LOG(logINFO)<< setw(10) << containedReads << " contained reads. (Need to change their mate-pair information)";
//	CLOCKSTOP;
//}

/**********************************************************************************************************************
	Hash table search found that a proper substring of read1 is a prefix or suffix of read2 or reverse complement of
	read2 (depents on the orient).
	orient 0 means prefix of forward of the read2
	orient 1 means suffix of forward of the read2
	orient 2 means prefix of reverse of the read2
	orient 3 means prefix of reverse of the read2
	We need to check if the remaining of the stings match to see if read2 is contained in read1.
**********************************************************************************************************************/
bool OverlapGraph::checkOverlapForContainedRead(Read *read1, Read *read2, UINT64 orient, UINT64 start)
{
	string string1=read1->getStringForward(); // Get the forward of read1
	UINT64 hashStringLength = hashTable->getHashStringLength(), lengthRemaining1, lengthRemaining2;
	string string2 = (orient == 0 || orient== 1) ? read2->getStringForward() : read2->getStringReverse(); // Get the string in read2 based on the orientation.
	if(orient == 0 || orient == 2)
	// orient 0
	//   >--------MMMMMMMMMMMMMMM*******------> read1      M means match found by hash table
	//            MMMMMMMMMMMMMMM*******>       read2      * means we need to check these characters for match
	//				OR
	// orient 2
	//	 >---*****MMMMMMMMMMMMMMM*******------> read1
	//		      MMMMMMMMMMMMMMM*******<	    Reverese complement of read2
	{
		lengthRemaining1 = string1.length() - start - hashStringLength; 	// This is the remaining of read1
		lengthRemaining2 = string2.length() - hashStringLength; 	// This is the remaining of read2
		if(lengthRemaining1 >= lengthRemaining2)
		{
			return string1.substr(start + hashStringLength, lengthRemaining2) == string2.substr(hashStringLength, lengthRemaining2); // If the remaining of the string match, then read2 is contained in read1
		}
	}
	else	// orient 1
	//   >---*****MMMMMMMMMMMMMMM-------------> read1      M means match found by hash table
	//      >*****MMMMMMMMMMMMMMM       		read2      * means we need to check these characters for match
	//				OR
	// orient 3
	//	 >---*****MMMMMMMMMMMMMMM-------------> read1
	//		<*****MMMMMMMMMMMMMMM				Reverse Complement of Read2
	{
		lengthRemaining1 = start;
		lengthRemaining2 = string2.length() - hashStringLength;
		if(lengthRemaining1 >= lengthRemaining2)
		{
			return string1.substr(start - lengthRemaining2, lengthRemaining2) == string2.substr(0, lengthRemaining2); // If the remaining of the string match, then read2 is contained in read1
		}
	}
	return false;

}

/**********************************************************************************************************************
	Checks if two read overlaps.
	Hash table search found that a proper substring of read1 is a prefix or suffix of read2 or reverse complement of
	read2 (depents on the orient).
	Orientation 0 means prefix of forward of the read2
	Orientation 1 means suffix of forward of the read2
	Orientation 2 means prefix of reverse of the read2
	Orientation 3 means prefix of reverse of the read2
	We need to check if the remaining of the stings match to see if read1 and read2 overlap.
**********************************************************************************************************************/
bool OverlapGraph::checkOverlap(Read *read1, Read *read2, UINT64 orient, UINT64 start)
{
	string string1=read1->getStringForward(); // Get the forward string of read1
	UINT64 hashStringLength = hashTable->getHashStringLength();
	string string2 = (orient == 0 || orient== 1) ? read2->getStringForward() : read2->getStringReverse(); // Get the string from read2 according to orient.
	if(orient == 0 || orient == 2)	// orient 0
	//   >--------MMMMMMMMMMMMMMM*************> 			read1      M means match found by hash table
	//            MMMMMMMMMMMMMMM*************------->      read2      * means we need to check these characters for match
	//				OR
	// orient 2
	//	 >---*****MMMMMMMMMMMMMMM*************> 			read1
	//		      MMMMMMMMMMMMMMM*************-------<	    Reverese complement of read2
	{
		if(string1.length()- start - hashStringLength >= string2.length() - hashStringLength) // The overlap must continue till the end.
			return false;
		return string1.substr(start + hashStringLength, string1.length()-(start + hashStringLength)) == string2.substr(hashStringLength,  string1.length()-(start + hashStringLength)); // If the remaining strings match.
	}
	else	// orient 1
	//   	>********MMMMMMMMMMMMMMM-------------> 			read1      M means match found by hash table
	//  >----********MMMMMMMMMMMMMMM       		    		read2      * means we need to check these characters for match
	//				OR
	// orient 3
	//	 	>********MMMMMMMMMMMMMMM-------------> 			read1
	//	<----********MMMMMMMMMMMMMMM						Reverse Complement of Read2
	{
		if(string2.length()-hashStringLength < start)
			return false;
		return string1.substr(0, start) == string2.substr(string2.length()-hashStringLength-start, start); // If the remaining strings match.
	}
}

/**********************************************************************************************************************
	Insert an edge in the overlap graph.
**********************************************************************************************************************/
bool OverlapGraph::insertEdge(Edge * edge)
{

	UINT64 ID = edge->getSourceRead()->getReadNumber(); // This is the source read.
	if(graph->at(ID)->empty()) 	// If there is no edge incident to the node
		numberOfNodes++;	// Then a new node is inserted in the graph. Number of nodes increased.
	graph->at(ID)->push_back(edge);	// Insert the edge in the list of edges of ID
		numberOfEdges++;	// Increase the number of edges.
	updateReadLocations(edge);	// If the current edge contains some reads, then we need to update their location information.
	return true;
}

/**********************************************************************************************************************
	Insert an edge in the graph.
**********************************************************************************************************************/
bool OverlapGraph::insertEdge(Read *read1, Read *read2, UINT8 orient, UINT16 overlapOffset)
{
	Edge * edge1 = new Edge(read1,read2,orient,overlapOffset);	// Create a new edge in the graph to insert.
	UINT16 overlapOffsetReverse = read2->getReadLength() + overlapOffset - read1->getReadLength();	// Set the overlap offset accordingly for the reverse edge. Note that read lengths are different.
	// If read lengths are the same. Then the reverse edge has the same overlap offset.
	Edge * edge2 = new Edge(read2,read1,twinEdgeOrientation(orient),overlapOffsetReverse);	// Create a new edge for the reverses string.

	edge1->setReverseEdge(edge2);	// Set the reverse edge pointer.
	edge2->setReverseEdge(edge1);	// Set the reverse edge pinter.
	insertEdge(edge1);	// Insert the edge in the overlap graph.
	insertEdge(edge2);	// Insert the edge in the overlap graph.
	return true;
}

/**********************************************************************************************************************
	This function prints the overlap graph in overlap_graph->gdl file. The graph can be viewed by
	aisee (free software available at http://www.aisee.com/)
	It also stores the contigs in a file.

**********************************************************************************************************************/
bool OverlapGraph::printGraph(string graphFileName, string contigFileName)
{
	CLOCKSTART;
	vector <Edge *> contigEdges;
	vector <Read *> isolateReads;
	UINT64 thickness, highestDegree = 0, highestDegreeNode;
	ofstream graphFilePointer; 
	//cout << "logging level is " << loglevel << endl;

	// Store the graph in a file.
	if (loglevel > 2 ) // only print graph in DEBUG mode
	{
		graphFilePointer.open(graphFileName.c_str());
		if(!graphFilePointer.is_open())
			MYEXIT("Unable to open file: "+graphFileName);

		graphFilePointer << "graph: {" << endl <<  "layoutalgorithm :forcedir" << endl <<  "fdmax:704" << endl <<  "tempmax:254" << endl <<  "tempmin:0" << endl <<  "temptreshold:3" << endl <<  "tempscheme:3" << endl <<  "tempfactor:1.08" << endl <<  "randomfactor:100" << endl <<  "gravity:0.0" << endl <<  "repulsion:161" << endl <<  "attraction:43" << endl <<  "ignore_singles:yes" << endl <<  "node.fontname:\"helvB10\"" << endl << "edge.fontname:\"helvB10\"" << endl <<  "node.shape:box" << endl <<  "node.width:80" << endl <<  "node.height:20" << endl <<  "node.borderwidth:1" << endl <<  "node.bordercolor:31" << endl;
		for(UINT64 i = 1; i<= dataSet->getNumberOfUniqueReads(); i++)
			if(!graph->at(i)->empty())
				graphFilePointer << "node: { title:\""<<i<<"\" label: \"" << i << "\" }" << endl;
	}

	for(UINT64 i = 1; i<= dataSet->getNumberOfUniqueReads(); i++)
	{
		if(!graph->at(i)->empty()) // if this read has some edge(s) connected to it
		{
			if(graph->at(i)->size() > highestDegree)
			{
				highestDegree = graph->at(i)->size();
				highestDegreeNode = i;

			}
			for(UINT64 j=0; j < graph->at(i)->size(); j++)
			{
				Edge * e = graph->at(i)->at(j);
				UINT64 source = e->getSourceRead()->getReadNumber(), destination = e->getDestinationRead()->getReadNumber();
				if(source < destination || (source == destination && e < e->getReverseEdge()) )
				{
					contigEdges.push_back(e); // List of contigs.
					if (loglevel > 2)
					{
						thickness = e->getListOfReads()->empty() ? 1: 3;
						if(e->getOrientation() == 0)
							graphFilePointer << "edge: { source:\"" << source << "\" target:\"" << destination << "\" thickness: " << thickness << " arrowstyle: none backarrowstyle: solid color: red label: \"(" << e->flow << "," << e->coverageDepth << "x," << e->getOverlapOffset() << "," << e->getListOfReads()->size() << ")\" }" << endl;
						else if(e->getOrientation() == 1)
							graphFilePointer << "edge: { source:\"" << source << "\" target:\"" << destination << "\" thickness: " << thickness << " backarrowstyle:solid arrowstyle:solid color: green label: \"(" << e->flow << "," << e->coverageDepth <<  "x," << e->getOverlapOffset() << "," << e->getListOfReads()->size() << ")\" }" << endl;
						else if(e->getOrientation() == 2)
							graphFilePointer << "edge: { source:\"" << source << "\" target:\"" << destination << "\" thickness: " << thickness << " arrowstyle: none color: blue label: \"(" << e->flow << "," << e->coverageDepth << "x," << e->getOverlapOffset() << "," << e->getListOfReads()->size() << ")\" }" << endl;
						else if(e->getOrientation() == 3)
							graphFilePointer << "edge: { source:\"" << source << "\" target:\"" << destination << "\" thickness: " << thickness << " arrowstyle:solid color: red label: \"(" << e->flow << "," << e->coverageDepth <<  "x," << e->getOverlapOffset() << "," << e->getListOfReads()->size() << ")\" }" << endl;
					}
				}
			}
		}
	// if this read does not have any  edge connected to it and is not on any composite edge
	//else if ( dataSet->getReadFromID(i)->superReadID == 0 && dataSet->getReadFromID(i)->getListOfEdgesForward()->size() == 0 && dataSet->getReadFromID(i)->getListOfEdgesReverse()->size() == 0 )
	//{
	//	// NEED TO FILL IN THE OPERATIONS HERE TODO FOR THE READS THAT DO NOT OVERLAP WITH ANYBODY ELSE	
	//	// How to print isolate nodes in graph??
	//	isolateReads.push_back(dataSet->getReadFromID(i));
	//}
	}
	if (loglevel > 2 )
	{
		graphFilePointer << "}";
		graphFilePointer.close();
		FILE_LOG(logINFO) << "Aisee graph written." << endl;
	}
	

	// Store the contigs in a file.
	if (contigEdges.size() > 0) // Sort the contigs by their length if there are edges in the graph.
	{
		sort(contigEdges.begin(),contigEdges.end(),compareEdgeByOverlapOffset); 
		reverse(contigEdges.begin(), contigEdges.end());
	}
	UINT64 sum = 0;
	//UINT64 isolate_sum = 0;
	ofstream contigFilePointer;
	contigFilePointer.open(contigFileName.c_str());
	if(!contigFilePointer.is_open())
		MYEXIT("Unable to open file: "+contigFileName);
	// All the edges
	if (contigEdges.size() > 0)
	{
		for(UINT64 i = 0; i < contigEdges.size(); i++) 
		{
			string s = getStringInEdge(contigEdges.at(i)); // get the string in the edge.
			contigFilePointer << ">contig_"<< i+1 << " Flow: " << setw(10) << contigEdges.at(i)->flow <<  " Edge  (" << setw(10) << contigEdges.at(i)->getSourceRead()->getReadNumber() << ", " << setw(10) << contigEdges.at(i)->getDestinationRead()->getReadNumber() << ") String Length: " << setw(10) << s.length() << " Coverage: " << setw(10) << contigEdges.at(i)->coverageDepth << endl;
			sum += s.length();
			UINT64 start=0;
			do
			{
				contigFilePointer << s.substr(start,100) << endl;  // save 100 BP in each line.
				start+=100;
			} while (start < s.length());
		}
	}
	// All the reads that do not have Edges? JJ
	//if (isolateReads.size() > 0)
	//{
	//	for(UINT64 i = 0; i < isolateReads.size(); i++)
	//	{
	//		string s = isolateReads.at(i)->getStringForward(); // get the string in the edge.
	//		contigFilePointer << ">isolate_"<< i+1 << " Read: " << setw(10) << isolateReads.at(i)->getReadNumber() << " String Length: " << setw(10) << s.length() << endl;
	//		isolate_sum += s.length();
	//		UINT64 start=0;
	//		do
	//		{
	//			contigFilePointer << s.substr(start,100) << endl;  // save 100 BP in each line.
	//			start+=100;
	//		} while (start < s.length());

	//	}
	//	FILE_LOG(logINFO)<< "Number of isolated Unique Reads in the graph: " << isolateReads.size();
	//	FILE_LOG(logINFO)<< "Total isolate read length: " << isolate_sum <<" BP";
	//}
	contigFilePointer.close();

	// Print some statistics about the graph.
	FILE_LOG(logINFO)<< "Total contig length: " << sum <<" BP";
	FILE_LOG(logINFO)<< "Number of Nodes in the graph: " << getNumberOfNodes();
	FILE_LOG(logINFO)<< "Number of Edges in the graph: " << getNumberOfEdges()/2;

	if (contigEdges.size() > 0)
	{
	UINT64 simEdges = 0, comEdges = 0, inEdges = 0, outEdges = 0;
	for(UINT64 i=0; i < graph->at(highestDegreeNode)->size(); i++)
	{
		if(graph->at(highestDegreeNode)->at(i)->getListOfReads()->empty())
			simEdges++;
		else
			comEdges++;
		if(graph->at(highestDegreeNode)->at(i)->getOrientation() == 0 || graph->at(highestDegreeNode)->at(i)->getOrientation() == 1)
			inEdges++;
		else
			outEdges++;
	}
	// Print some more statistics on the node with highest degree.
	FILE_LOG(logINFO)<< "Highest Degree Read " << highestDegreeNode << " has " << highestDegree << " neighbors.";
	FILE_LOG(logINFO) << "In Edges: " << inEdges << " Out Edges: " << outEdges << " Simple Edges: " << simEdges << " Composite Edges: " << comEdges;
	FILE_LOG(logDEBUG2)<< "Read forward string: " << dataSet->getReadFromID(highestDegreeNode)->getStringForward() << endl;
	}

	CLOCKSTOP;
	return true;
}

/**********************************************************************************************************************
 * Insert all edges of a read in the overlap graph.  
 * If a read is already explored, it won't be explored against for another read again.
**********************************************************************************************************************/
bool OverlapGraph::insertAllEdgesOfRead(UINT64 readNumber, vector<nodeType> * exploredReads)
{
	Read *read1 = dataSet->getReadFromID(readNumber); 	// Get the current read read1.
	seqan::DnaString readDnaString = read1->getDnaStringForward(); 	// Get the forward string of read1.
	seqan::DnaString subString;
	for(UINT64 j = readNumber + 1 ; j < dataSet->getNumberOfUniqueReads(); j++) // For each read with start coordinate bigger than the current read, see if they overlap
	{
		subString = readString.substr(j,hashTable->getHashStringLength());  // Get the proper substring s of read1.
		vector<UINT64> * listOfReads=hashTable->getListOfReads(subString); // Search the string in the hash table.
		if(!listOfReads->empty()) // If there are some reads that contain s as prefix or suffix of the read or their reverse complement
		{
			for(UINT64 k = 0; k < listOfReads->size(); k++) // For each such reads.
			{
				UINT64 data = listOfReads->at(k);	// We used bit operations in the hash table. Most significant 2 bits store orientation and least significant 62 bits store read ID.
				UINT16 overlapOffset;
				UINT8 orientation;
				Read *read2 = dataSet->getReadFromID(data & 0X3FFFFFFFFFFFFFFF); 	// Least significant 62 bits store the read number.
				if(exploredReads->at(read2->getReadNumber())!= UNEXPLORED) 	// No need to discover the same edge again. All edges of read2 is already inserted in the graph.
						continue;
				if(read1->superReadID == 0 && read2->superReadID == 0 && checkOverlap(read1,read2,(data >> 62),j)) // Both read need to be non contained.
				{
					switch (data >> 62) // Most significant 2 bit represents  00 - prefix forward, 01 - suffix forward, 10 -  prefix reverse, 11 -  suffix reverse.
					{
						case 0: orientation = 3; overlapOffset = read1->getReadLength() - j; break; 	// 3 = r1>------->r2
						case 1: orientation = 0; overlapOffset = hashTable->getHashStringLength() + j; break; 	// 0 = r1<-------<r2
						case 2: orientation = 2; overlapOffset = read1->getReadLength() - j; break; 	// 2 = r1>-------<r2
						case 3: orientation = 1; overlapOffset = hashTable->getHashStringLength() + j; break; 	// 1 = r2<------->r2
					}
					insertEdge(read1,read2,orientation,read1->getStringForward().length()-overlapOffset); 	// Insert the edge in the graph.
				}
			}
		}
	}
	if(graph->at(readNumber)->size() != 0)
		sort(graph->at(readNumber)->begin(),graph->at(readNumber)->end(), compareEdgeByOverlapOffset); // Sort the list of edges of the current node according to the overlap offset (ascending).
	return true;
}

/**********************************************************************************************************************
	Mark all the transitive edges of a read.
	For Details: E.W. Myers. The fragment assembly string graph. Bioinformatics, 21(suppl 2):ii79-ii85, 2005.
**********************************************************************************************************************/
bool OverlapGraph::markTransitiveEdges(UINT64 readNumber, vector<markType> * markedNodes)
{

	for(UINT64 i = 0; i < graph->at(readNumber)->size(); i++) // Mark all the neighbours of the current read as INPLAY
		markedNodes->at(graph->at(readNumber)->at(i)->getDestinationRead()->getReadNumber()) = INPLAY; // Inplay

	for(UINT64 i = 0; i < graph->at(readNumber)->size(); i++) // Traverse through the list of edges according to their overlap offset.
	{
		UINT64 read2 = graph->at(readNumber)->at(i)->getDestinationRead()->getReadNumber(); // For each neighbor
		if(markedNodes->at(read2) == INPLAY) 	// If the neighbor is marked as INPLAY
		{
			for(UINT64 j = 0; j < graph->at(read2)->size(); j++)
			{
				UINT64 read3 = graph->at(read2)->at(j)->getDestinationRead()->getReadNumber(); // Get the neighbors neighbors
				if(markedNodes->at(read3) == INPLAY)
				{

					UINT8 type1 = graph->at(readNumber)->at(i)->getOrientation();
					UINT8 type2 = graph->at(read2)->at(j)->getOrientation();
					if((type1 == 0 ||  type1 == 2) && (type2==0 || type2==1)) 	// Check edge orientation
						markedNodes->at(read3) = ELIMINATED; 	// Mark as ELIMINATED
					else if((type1==1||type1==3) && (type2==2 || type2==3)) 	// Check edge orientation
						markedNodes->at(read3) = ELIMINATED; 	// Mark as ELIMINATED
				}
			}
		}
	}
	for(UINT64 i = 0;i < graph->at(readNumber)->size(); i++)
	{
		if(markedNodes->at(graph->at(readNumber)->at(i)->getDestinationRead()->getReadNumber()) == ELIMINATED) // Current read to a node marked as ELIMINATED
		{
			graph->at(readNumber)->at(i)->transitiveRemovalFlag = true; 	// Mark this edge as transitive edge. Will remove this edge later.
			graph->at(readNumber)->at(i)->getReverseEdge()->transitiveRemovalFlag = true;	// Mark also the reverse edge. Will remove this edge later.
		}
	}

	for(UINT64 i = 0; i < graph->at(readNumber)->size(); i++)
		markedNodes->at(graph->at(readNumber)->at(i)->getDestinationRead()->getReadNumber()) = VACANT; 	// Change back all the variables modified in this function to VACANT

	markedNodes->at(readNumber) = VACANT; 	// Mark as vacant.
	return true;
}

/**********************************************************************************************************************
	Remove all transitive edges of a given read.
	For Details: E.W. Myers. The fragment assembly string graph. Bioinformatics, 21(suppl 2):ii79-ii85, 2005.
**********************************************************************************************************************/
bool OverlapGraph::removeTransitiveEdges(UINT64 readNumber)
{
	for(UINT64 index = 0; index < graph->at(readNumber)->size(); index++)  	// Go through the list of edges of the current read.
	{
		if(graph->at(readNumber)->at(index)->transitiveRemovalFlag == true)	// This edge is marked as transitive. We will first remove the reverese edge.
		{
			Edge *twinEdge = graph->at(readNumber)->at(index)->getReverseEdge();
			UINT64 ID = twinEdge->getSourceRead()->getReadNumber();
			for(UINT64 index1 = 0; index1 < graph->at(ID)->size(); index1++) 	// Get the reverse edge first
			{
				if(graph->at(ID)->at(index1) == twinEdge)
				{
					delete twinEdge;
					graph->at(ID)->at(index1) = graph->at(ID)->at(graph->at(ID)->size()-1); // Move the transitive edges at the back of the list and remove.
					graph->at(ID)->pop_back();
					if(graph->at(ID)->empty())
						numberOfNodes--;
					numberOfEdges--;
					break;
				}
			}
		}
	}
	UINT64 j=0;
	for(UINT64 index=0; index < graph->at(readNumber)->size(); index++) // Then we will remove all the transitive edges of the current read.
	{
		if(graph->at(readNumber)->at(index)->transitiveRemovalFlag == false)	// We move all the non-transitive edges at the beginning of the list
			graph->at(readNumber)->at(j++) = graph->at(readNumber)->at(index);
		else	// Free the transitive edge
		{
			numberOfEdges--;
			delete graph->at(readNumber)->at(index);
		}
	}
	graph->at(readNumber)->resize(j);
	if(graph->at(readNumber)->empty())
		numberOfNodes--;
	return true;
}

/**********************************************************************************************************************
	Contract composite paths in the overlap graph.
	u*-------->v>---------*w  => u*--------------------*w
	u*--------<v<---------*w  => u*--------------------*w
**********************************************************************************************************************/
UINT64 OverlapGraph::contractCompositePaths(void)
{
	CLOCKSTART;
	UINT64 counter = 0;
	for(UINT64 index = 1 ; index < graph->size(); index++)
	{
		if(graph->at(index)->size() == 2) // Check if the node has only two edges.
		{
			Edge *edge1 = graph->at(index)->at(0);  // First edge.
			Edge *edge2 = graph->at(index)->at(1);  // Second Edge.
			if(flowComputed == true || !isEdgePresent(edge1->getDestinationRead()->getReadNumber(), edge2->getDestinationRead()->getReadNumber()))
	// Before flow is computed we do not insert multiple edges between the same endpoints.
			{
				if( matchEdgeType(edge1->getReverseEdge(), edge2) && edge1->getSourceRead() != edge1->getDestinationRead()) // One incoming edge and one outgoing edge.
				{
					mergeEdges(edge1->getReverseEdge(),edge2);	// Merge the edges.
					counter++;	// Counter how many edges merged.
				}
			}
		}

	}
	FILE_LOG(logDEBUG) << setw(10) << counter << " composite Edges merged.";
	CLOCKSTOP;
	return counter;
}

/**********************************************************************************************************************
	Merge two edges in the overlap graph.
**********************************************************************************************************************/
bool OverlapGraph::mergeEdges(Edge *edge1, Edge *edge2)
{
	Edge *edgeForward = new Edge(); // New forward edge.
	Read *read1 = edge1->getSourceRead(), *read2 = edge2->getDestinationRead();
	Edge *edgeReverse = new Edge(); // New reverse edge.

	UINT8 orientationForward = mergedEdgeOrientation(edge1,edge2);	// Orientation of the forward edge.
	UINT8 orientationReverse = twinEdgeOrientation(orientationForward);	// Orientation of the reverse edge.

	vector<UINT64> * listReadsForward = new vector<UINT64>;	// List of reads in the forward edge.
	vector<UINT16> * listOverlapOffsetsForward= new vector<UINT16>;	// List of Overlaps in the forward edge.
	vector<UINT8> * listOrientationsForward = new vector<UINT8>;	// List of orientations in the forward edge.

	mergeList(edge1, edge2, listReadsForward, listOverlapOffsetsForward, listOrientationsForward); // Merge the lists from the two edges.

	edgeForward->makeEdge(read1,read2,orientationForward, edge1->getOverlapOffset() + edge2->getOverlapOffset(), listReadsForward, listOverlapOffsetsForward, listOrientationsForward); // Make the forward edge

	vector<UINT64> * listReadsReverse = new vector<UINT64>;	// List of reads in the reverse edge.
	vector<UINT16> * listOverlapOffsetsReverse= new vector<UINT16>;	// List of overlaps in the reverse edge.
	vector<UINT8> * listOrientationsReverse = new vector<UINT8>;	// List of orientations in the reverse edge.

	mergeList(edge2->getReverseEdge(),edge1->getReverseEdge(), listReadsReverse, listOverlapOffsetsReverse,listOrientationsReverse);	// Merge the lists from the two reverse edges.

	edgeReverse->makeEdge(read2, read1, orientationReverse, edge2->getReverseEdge()->getOverlapOffset() + edge1->getReverseEdge()->getOverlapOffset(), listReadsReverse, listOverlapOffsetsReverse, listOrientationsReverse); // Make the reverse edge

	edgeForward->setReverseEdge(edgeReverse);	// Update the reverse edge pointer.
	edgeReverse->setReverseEdge(edgeForward);	// Update the reverse edge pointer.


	UINT16 flow = min(edge1->flow,edge2->flow);	// We take the minimum of the flow in the new edge.

	edgeForward->flow = flow;	// Modify the flow in the new forward edge.
	edgeReverse->flow = flow;	// Modify the flow in the new reverse edge.


	insertEdge(edgeForward);	// Insert the new forward edge in the graph.
	insertEdge(edgeReverse);	// Insert the new reverse edge in the graph.

	edge1->flow = edge1->flow - flow;	// Remove the used flow from edge1.
	edge1->getReverseEdge()->flow = edge1->flow;	// Revmod the used flow from the reverse of edge1.

	edge2->flow = edge2->flow - flow;	// Remove the used flow from edge2.
	edge2->getReverseEdge()->flow = edge2->flow;	// Remove the used flow from teh reverse of edge2.

	if(edge1->flow == 0 || flow == 0)	// If no flow left in edge1
		removeEdge(edge1);	// edge1 is deleted from the graph.
	if(edge2->flow == 0 || flow == 0)	// If now flow left in edge2
		removeEdge(edge2);	// edge 2 is deleted from the graph.

	return true;

}

/**********************************************************************************************************************
	Merge the list of reads, list of overlap offsets and list of orientations of two edges.
**********************************************************************************************************************/

bool OverlapGraph::mergeList(Edge *edge1, Edge *edge2, vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets, vector<UINT8> * listOrientations)
{
	UINT64 sum = 0;
	for(UINT64 i = 0; i < edge1->getListOfOrientations()->size(); i++) 	// Take the list from edge1.
	{
		listReads->push_back(edge1->getListOfReads()->at(i));
		listOverlapOffsets->push_back(edge1->getListOfOverlapOffsets()->at(i));
		listOrientations->push_back(edge1->getListOfOrientations()->at(i));
		sum += edge1->getListOfOverlapOffsets()->at(i);
	}
	listReads->push_back(edge1->getDestinationRead()->getReadNumber()); 	// Insert the common node of the two edges

	listOverlapOffsets->push_back(edge1->getOverlapOffset() - sum);	// Get the overlap offset.

	if(edge1->getOrientation() == 1 || edge1->getOrientation() == 3)	// Orientation of the common node. Depends on the orientation of the edges.
		listOrientations->push_back(1);
	else
		listOrientations->push_back(0);
	for(UINT64 i = 0; i < edge2->getListOfOrientations()->size(); i++)	// take the list from edge2.
	{
		listReads->push_back(edge2->getListOfReads()->at(i));
		listOverlapOffsets->push_back(edge2->getListOfOverlapOffsets()->at(i));
		listOrientations->push_back(edge2->getListOfOrientations()->at(i));
	}
	return true;
}




/**********************************************************************************************************************
	Orientation of the merged edge.
	e1(u,v) e2(v,w)
	Orientation e1 0 = u<-------<v		Orientation e2 0 = v<-------<w
	Orientation e1 1 = u<------->v		Orientation e2 1 = v<------->w
	Orientation e1 2 = u>-------<v		Orientation e2 2 = v>-------<w
	Orientation e1 3 = u>------->v		Orientation e2 3 = v>------->w

	0 + 0 = u<-------<w
	0 + 1 = u<------->w
	...................

**********************************************************************************************************************/
UINT8 OverlapGraph::mergedEdgeOrientation(Edge *edge1, Edge *edge2)
{
	UINT8 or1 = edge1->getOrientation(), or2 = edge2->getOrientation(),returnValue;
	if(or1 == 0 && or2 == 0)
		returnValue = 0;
	else if(or1 == 0 && or2 == 1)
		returnValue = 1;
	else if(or1 == 1 && or2 == 2)
		returnValue = 0;
	else if(or1 == 1 && or2 == 3)
		returnValue = 1;
	else if(or1 == 2 && or2 == 0)
		returnValue = 2;
	else if(or1 == 2 && or2 == 1)
		returnValue = 3;
	else if(or1 == 3 && or2 == 2)
			returnValue = 2;
	else if(or1 == 3 && or2 == 3)
			returnValue = 3;
	else
	{
		FILE_LOG(logDEBUG)<<(int)or1<<" "<<(int)or2;
		MYEXIT("Unable to merge.")
	}
	return returnValue;
}





/**********************************************************************************************************************
	Orientation of a reverse edge;
	Twin edge of Orientation 0 = <-------< is Orientation 3 = >------->
	Twin edge of Orientation 1 = <-------> is Orientation 1 = <------->
	Twin edge of Orientation 2 = >-------< is Orientation 2 = >-------<
	Twin edge of Orientation 3 = >-------> is Orientation 0 = <-------<
**********************************************************************************************************************/
UINT8 OverlapGraph::twinEdgeOrientation(UINT8 orientation)
{
	UINT8 returnValue;
	if(orientation == 0)
		returnValue = 3;
	else if(orientation == 1)
		returnValue = 1;
	else if(orientation == 2)
		returnValue = 2;
	else if(orientation == 3)
		returnValue = 0;
	else
		MYEXIT("Unsupported edge orientation.")
	return returnValue;
}




/**********************************************************************************************************************
	remove an edge from the overlap graph.
**********************************************************************************************************************/
bool OverlapGraph::removeEdge(Edge *edge)
{
	removeReadLocations(edge);	// If the current edge contains some reads. We have to update their location formation.
	removeReadLocations(edge->getReverseEdge());	// Same for the reverse edge.
	Edge *twinEdge = edge->getReverseEdge();
	UINT64 ID1 = edge->getSourceRead()->getReadNumber(), ID2 = edge->getDestinationRead()->getReadNumber();  // Get the source and destation read IDs.
	for(UINT64 i = 0; i< graph->at(ID2)->size(); i++) // Delete the twin edge first.
	{
		if(graph->at(ID2)->at(i) == twinEdge)
		{
			delete graph->at(ID2)->at(i);
			graph->at(ID2)->at(i) = graph->at(ID2)->at(graph->at(ID2)->size()-1);
			graph->at(ID2)->pop_back();
			if(graph->at(ID2)->empty())
				numberOfNodes--;
			numberOfEdges--;
			break;
		}
	}
	for(UINT64 i = 0; i< graph->at(ID1)->size(); i++) // Delete the edge then.
	{
		if(graph->at(ID1)->at(i) == edge)
		{
			delete graph->at(ID1)->at(i);
			graph->at(ID1)->at(i) = graph->at(ID1)->at(graph->at(ID1)->size()-1);
			graph->at(ID1)->pop_back();
			if(graph->at(ID1)->empty())
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
	Sometimes these edges could have a lot of basepairs and probably should be retained, instead of being deleted. JJ
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
				if(edge->getSourceRead()->getReadNumber() < edge->getDestinationRead()->getReadNumber() && edge->getListOfReads()->empty() && edge->flow == 0 ) 
				{
					listOfEdges.push_back(edge); // Put in the list of edges to be removed.
					FILE_LOG(logDEBUG1)  << "removing edge ("<< edge->getSourceRead()->getReadNumber()<<","  << edge->getDestinationRead()->getReadNumber()<<") OverlapOffset : " << setw(10) << edge->getOverlapOffset(); 
				}
			}
		}
	}
	for(UINT64 i = 0 ; i < listOfEdges.size(); i++)
		removeEdge(listOfEdges.at(i));	// remove the edges from the list.
	FILE_LOG(logDEBUG) << "Simple edges without flow removed: " << listOfEdges.size();
	CLOCKSTOP;
	return listOfEdges.size();
}


/**********************************************************************************************************************
	Remove nodes with all simple edges and all same arrow type
**********************************************************************************************************************/
UINT64 OverlapGraph::removeDeadEndNodes(void)
{
	CLOCKSTART;
	vector <UINT64> listOfNodes; // for saving nodes that should be deleted
	
	if ( getNumberOfNodes() < 3){
		FILE_LOG(logINFO) << "Number of nodes in graph is " << getNumberOfNodes() << ", cannot reduce any more";
		return listOfNodes.size();
	}
	else{
		for(UINT64 i = 1; i < graph->size(); i++) // For each read.
		{
			if(!graph->at(i)->empty())	// If the read has some edges.
			{
				UINT64 flag = 0, inEdge = 0 , outEdge = 0;
				for(UINT64 j=0; j < graph->at(i)->size(); j++) // For each edge
				{
					Edge * edge = graph->at(i)->at(j);
	// Break case:
	// 1. if composite edge is more than deadEndLength reads
	// 2. if the edge is loop for the current node
	// 3. if the string formed in this edge is longer than deadEndBp (500) (JJ)
	// Then flag=1 and exit the loop
	//if(edge->getListOfReads()->size() > deadEndLength || getStringInEdge(edge).size() > deadEndBp || edge->getSourceRead()->getReadNumber() == edge->getDestinationRead()->getReadNumber())
					if(edge->getListOfReads()->size() > deadEndLength || edge->getSourceRead()->getReadNumber() == edge->getDestinationRead()->getReadNumber())
					{
						flag = 1;
						break;
					}

					if(edge->getOrientation() == 0 || edge->getOrientation() == 1)
						inEdge++;
					else
						outEdge++;
				}
				if(flag == 0) // If not break case
				{
					if( (inEdge > 0 && outEdge == 0) || (inEdge == 0 && outEdge > 0)) // only one type of simple edges
					{
						listOfNodes.push_back(i);
					}
				}
			}
		}
	
		UINT64 edgesRemoved = 0;
		vector <Edge *> listOfEdges;
		for(UINT64 i = 0 ; i < listOfNodes.size(); i++)
		{
			listOfEdges.clear();
			if(!graph->at(listOfNodes.at(i))->empty())	// If the read has some edges.
			{
				edgesRemoved += graph->at(listOfNodes.at(i))->size();
	//FILE_LOG(logDEBUG1) << "Removing dead-end node " << listOfNodes.at(i) ;
				for(UINT64 j=0; j < graph->at(listOfNodes.at(i))->size(); j++) // For each edge
				{
					listOfEdges.push_back(graph->at(listOfNodes.at(i))->at(j));
					FILE_LOG(logDEBUG1) << "Removing dead-end node edge ( " <<  graph->at(listOfNodes.at(i))->at(j)->getSourceRead()->getReadNumber() << ", " << graph->at(listOfNodes.at(i))->at(j)->getDestinationRead()->getReadNumber() << ") Length: " << getStringInEdge(graph->at(listOfNodes.at(i))->at(j)).length();
				}
				for(UINT64 j = 0; j< listOfEdges.size(); j++)
				{
					removeEdge(listOfEdges.at(j));	// Remove all the edges of the current node.
				}
			}
		}
		FILE_LOG(logDEBUG)<< "Dead-end nodes removed: " << listOfNodes.size();
		FILE_LOG(logDEBUG)<<  "Total Edges removed: " << edgesRemoved;
		CLOCKSTOP;
		return listOfNodes.size();
	}
}



/**********************************************************************************************************************
	Estimate the genome size. The estimated genome size is not used anywere in the program.
**********************************************************************************************************************/
bool OverlapGraph::estimateGenomeSize(void)
{
	CLOCKSTART;
	UINT64 delta, deltaSum, freqSum, currentGenomeSize = 0, previousGenomeSize = 0, freq, counter = 0;
	do
	{
		counter++; 	// To avoid infinite loop.
		deltaSum = 0;
		freqSum = 0;
		float a_statistic;
		for(UINT64 i = 1; i < graph->size(); i++)
		{
			for(UINT64 j = 0; j < graph->at(i)->size(); j++)
			{
				Edge * currentEdge = graph->at(i)->at(j);
				if(currentEdge->getSourceRead()->getReadNumber() < currentEdge->getDestinationRead()->getReadNumber())
				{
					delta = currentEdge->getOverlapOffset();
					freq = 0;
					for(UINT64 k = 0; k < currentEdge->getListOfReads()->size(); k++)
						freq +=  dataSet->getReadFromID( currentEdge->getListOfReads()->at(k) )->getFrequency();
					if(previousGenomeSize != 0) 	// If not the first estimation
					{
						a_statistic=(float)((float)delta*(float)((float)dataSet->getNumberOfReads() / (float)previousGenomeSize) - (float)((float)freq * log((float)2)));
						if(a_statistic >= aStatisticsThreshold && delta >= minDelta)
						{
							deltaSum += delta;
							freqSum += freq;
						}
					}
					else if(currentEdge->getOverlapOffset() > 500)  // Only take edges longer than 500 for the first estimation
					{
						deltaSum += delta;
						freqSum += freq;
					}
				}
			}
		}
		previousGenomeSize = currentGenomeSize;
		currentGenomeSize = (int)((float)dataSet->getNumberOfReads()/(float)freqSum*(float)deltaSum);
		FILE_LOG(logINFO)<<"Current estimated genome size: " << currentGenomeSize;
	}while (currentGenomeSize != previousGenomeSize && counter < 10);
	estimatedGenomeSize=currentGenomeSize;
	FILE_LOG(logINFO)<<"Final estimated genome size: " << getEstimatedGenomeSize();
	CLOCKSTOP;
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
		if(edge->getListOfOrientations()->at(i) == 1)	// Orientation 1 means that the forward string of the read is contained in this edge.
		{
			read->getListOfEdgesForward()->push_back(edge);	// Insert the current edge in the list of forward edges.
			read->getLocationOnEdgeForward()->push_back(distance);	// Also insert the distance within the edge.
			read->getListOfEdgesForward()->resize(read->getListOfEdgesForward()->size());
			read->getLocationOnEdgeForward()->resize(read->getLocationOnEdgeForward()->size());
		}
		else	// Orientation 0 means that the reverser string of the read is contained in this edge.
		{
			read->getListOfEdgesReverse()->push_back(edge);	// Insert the current edge in the list of reverser edges.
			read->getLocationOnEdgeReverse()->push_back(distance);	// Also insert the distance withing the edge.
			read->getListOfEdgesReverse()->resize(read->getListOfEdgesReverse()->size());
			read->getLocationOnEdgeReverse()->resize(read->getLocationOnEdgeReverse()->size());
		}
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
				read->getLocationOnEdgeForward()->at(j) = read->getLocationOnEdgeForward()->at(read->getLocationOnEdgeForward()->size()-1);

				read->getListOfEdgesForward()->pop_back();
				read->getLocationOnEdgeForward()->pop_back();

				read->getListOfEdgesForward()->resize(read->getListOfEdgesForward()->size());
				read->getLocationOnEdgeForward()->resize(read->getLocationOnEdgeForward()->size());
			}
		}

		for(UINT64 j = 0; j < read->getListOfEdgesReverse()->size(); j++) 	// Check the list of edges that contain reverse of this read.
		{
			if(read->getListOfEdgesReverse()->at(j) == edge) 	// Remove the current edge from the list.
			{
				read->getListOfEdgesReverse()->at(j) = read->getListOfEdgesReverse()->at(read->getListOfEdgesReverse()->size()-1);
				read->getLocationOnEdgeReverse()->at(j) = read->getLocationOnEdgeReverse()->at(read->getLocationOnEdgeReverse()->size()-1);

				read->getListOfEdgesReverse()->pop_back();
				read->getLocationOnEdgeReverse()->pop_back();

				read->getListOfEdgesReverse()->resize(read->getListOfEdgesReverse()->size());
				read->getLocationOnEdgeReverse()->resize(read->getLocationOnEdgeReverse()->size());
			}
		}
	}
	return true;
}

/**********************************************************************************************************************
	Calculate Mean and Standard Deviation of insert size of each dataset.
**********************************************************************************************************************/
bool OverlapGraph::calculateMeanAndSdOfInsertSize(void)
{
	CLOCKSTART;
	if(dataSet->pairedEndDatasetFileNames.size() == 0)	// If no paired-edge dataset
		return true;	// No need to proceed.

	vector<UINT64> *insertSizes = new vector<UINT64>;
	vector<Edge *> listOfEdgesRead1, listOfEdgesRead2;
	vector<UINT64> locationOnEdgeRead1, locationOnEdgeRead2;

	for(UINT64 d = 0; d < dataSet->pairedEndDatasetFileNames.size(); d++)	// For each dataset.
	{

		FILE_LOG(logINFO) << "Calculating mean and SD of dataset: " << d;
		insertSizes->clear();

		for(UINT64 i = 1; i < graph->size(); i++)	// for each read in the dataset
		{
			Read *read1 = dataSet->getReadFromID(i), *read2;	// Get a read read1 in the graph.

			for(UINT64 j = 0; j < read1->getMatePairList()->size(); j++) 	// For each mate-pair read2
			{
				if(read1->getMatePairList()->at(j).datasetNumber == d)	// if it is in the current dataset.
				{
					listOfEdgesRead1.clear();  locationOnEdgeRead1.clear();
					listOfEdgesRead1.insert(listOfEdgesRead1.end(), read1->getListOfEdgesForward()->begin(), read1->getListOfEdgesForward()->end());	// All the edges that contain forward string of read1
					listOfEdgesRead1.insert(listOfEdgesRead1.end(), read1->getListOfEdgesReverse()->begin(), read1->getListOfEdgesReverse()->end());	// All the edges that contain reverse string of read1

					locationOnEdgeRead1.insert(locationOnEdgeRead1.end(), read1->getLocationOnEdgeForward()->begin(), read1->getLocationOnEdgeForward()->end());  	// Location on the corresponding edge.
					locationOnEdgeRead1.insert(locationOnEdgeRead1.end(), read1->getLocationOnEdgeReverse()->begin(), read1->getLocationOnEdgeReverse()->end());	// Location on the corresponding edge.

					read2 = dataSet->getReadFromID(read1->getMatePairList()->at(j).matePairID);	// Get the read2 object.
					listOfEdgesRead2.clear();  locationOnEdgeRead2.clear();

					listOfEdgesRead2.insert(listOfEdgesRead2.end(), read2->getListOfEdgesForward()->begin(), read2->getListOfEdgesForward()->end());	// All the edges that contain forward string of read2
					listOfEdgesRead2.insert(listOfEdgesRead2.end(), read2->getListOfEdgesReverse()->begin(), read2->getListOfEdgesReverse()->end());	// All the edges that contain reverse string of read1

					locationOnEdgeRead2.insert(locationOnEdgeRead2.end(), read2->getLocationOnEdgeForward()->begin(), read2->getLocationOnEdgeForward()->end());	// Location on the corresponding edge.
					locationOnEdgeRead2.insert(locationOnEdgeRead2.end(), read2->getLocationOnEdgeReverse()->begin(), read2->getLocationOnEdgeReverse()->end());	// Location on the corresponding edge.

					for(UINT64 k = 0; k < listOfEdgesRead1.size(); k++)	// For each edge containing read1
					{
						for(UINT64 l = 0; l < listOfEdgesRead2.size(); l++)	// For each edge containing read2
						{
							if(listOfEdgesRead1.at(k) == listOfEdgesRead2.at(l) && locationOnEdgeRead1.at(k) > locationOnEdgeRead2.at(l))	// Both reads are on the same edge
							{
									if(locationOnEdgeRead1.at(k) - locationOnEdgeRead2.at(l) < 1000)	// Distance between the two edges is less than 1000. Some times some mate pairs are far off the actual value. This upper bound is used to get only good mate pairs. We may need to change the threshold for datasets with longer insert size.
									{
										insertSizes->push_back(locationOnEdgeRead1.at(k) - locationOnEdgeRead2.at(l));	// Insert the distance between the  reads in the list.
									}
							}
						}
					}
				}
			}
		}

		UINT64 sum = 0, variance=0;
		if(insertSizes->size() == 0) // If no insert size found
		{
			FILE_LOG(logINFO) << "No insert-size found for dataset: " << d;
			meanOfInsertSizes.push_back(0);
			sdOfInsertSizes.push_back(0);
			continue;
		}

		for(UINT64 i = 0; i < insertSizes->size(); i++)
			sum += insertSizes->at(i);	// Calculate the sum of all the insert sizes of the current dataset.

		meanOfInsertSizes.push_back(sum/insertSizes->size());	// Calculate the mean. In push it in the variable of the OverlapGraph object.

		for(UINT64 i = 0; i < insertSizes->size(); i++)	// Calculate the variance of the current dataset.
			variance += (meanOfInsertSizes.at(d) - insertSizes->at(i)) * (meanOfInsertSizes.at(d) - insertSizes->at(i));

		sdOfInsertSizes.push_back(sqrt(variance/insertSizes->size()));  // Calculate and insert the standard deviation.


	// Print the values of the current dataset.
		FILE_LOG(logINFO) << "Mean set to: " << meanOfInsertSizes.at(d);
		FILE_LOG(logINFO) << "SD set to: " << sdOfInsertSizes.at(d);
		FILE_LOG(logINFO) << "Reads on same edge: " << insertSizes->size();

	}

	delete insertSizes;
	CLOCKSTOP;
	return true;
}

/**********************************************************************************************************************
	Save the unitig graph in a text file
**********************************************************************************************************************/
bool OverlapGraph::saveGraphToFile(string fileName)
{
	CLOCKSTART;
	ofstream filePointer;
	filePointer.open(fileName.c_str());
	if(!filePointer.is_open())
		MYEXIT("Unable to open file: "+fileName);

	vector<UINT64> *list = new vector<UINT64>;
	for(UINT64 i = 1; i < graph->size(); i++) // for each node
	{
		if(!graph->at(i)->empty())
		{
			for(UINT64 j = 0; j < graph->at(i)->size(); j++)	// for each edge of the node
			{
				Edge * e = graph->at(i)->at(j);
				UINT64 source = e->getSourceRead()->getReadNumber();
				UINT64 destination = e->getDestinationRead()->getReadNumber();
				if(source < destination || (source == destination && e < e->getReverseEdge()))
				{
					list->push_back(source);	// store the edge information first
					list->push_back(destination);
					list->push_back(e->getOrientation());
					list->push_back(e->getOverlapOffset());
					list->push_back(e->getListOfReads()->size());	// store the size of the vector of elements in the edge.
					for(UINT64 k = 0; k < e->getListOfReads()->size(); k++)	// store information about the reads in the edge.
					{
						list->push_back(e->getListOfReads()->at(k));
						list->push_back(e->getListOfOverlapOffsets()->at(k));
						list->push_back(e->getListOfOrientations()->at(k));
					}
				}
			}
		}
	}

	for(UINT64 i = 0; i < list->size(); i++)	// store in a file for future use.
		filePointer<<list->at(i)<<endl;
	filePointer.close();
	delete list;
	CLOCKSTOP;
	return true;
}





/**********************************************************************************************************************
	Read the unitig graph from a binary file
**********************************************************************************************************************/
bool OverlapGraph::readGraphFromFile(string fileName)
{
	CLOCKSTART;
	ifstream filePointer;
	filePointer.open(fileName.c_str());
	if(!filePointer.is_open())
		MYEXIT("Unable to open file: "+fileName);

	graph = new vector< vector<Edge *> * >;
	graph->reserve(dataSet->getNumberOfUniqueReads()+1);
	for(UINT64 i = 0; i <= dataSet->getNumberOfUniqueReads(); i++) // initialize the list
	{
		vector<Edge *> *newList = new vector<Edge *>;
		graph->push_back(newList);
	}

	vector<UINT64> *list =  new vector<UINT64>;

	while(filePointer.good())	// read from the file.
	{
		UINT64 temp;
		filePointer>>temp;
		list->push_back(temp);
	};
	filePointer.close();
	for(UINT64 i = 0; i < list->size() -1;) // parse the numbers.
	{

		UINT64 source = list->at(i++);	// first number is the source read id.
		UINT64 destination = list->at(i++);	// destination read id.
		UINT64 orientation = list->at(i++);	// orientation of the edge.
		UINT64 overlapOffset = list->at(i++);	// overlap offset of the edge
		UINT64 numberOfReadsInEdge = list->at(i++);	// Array size
		vector<UINT64> *listReads = new vector<UINT64>;
		vector<UINT16> *listOverlapOffsets = new vector<UINT16>;
		vector<UINT8> *listOrientations = new vector<UINT8>;
		UINT64 length = 0;
		for(UINT64 j = 0; j < 3 * numberOfReadsInEdge; j += 3)	// Read all the three arrays. list of reads, overlap offsets and orientation of the reads.
		{
			listReads->push_back(list->at(i + j));
			listOverlapOffsets->push_back(list->at(i + j + 1));
			listOrientations->push_back(list->at(i + j + 2));
			length += list->at(i + j + 1);
		}

		Edge *edgeForward = new Edge();	// create an edge
		Edge *edgeReverse = new Edge();
		Read *read1, *read2;

		read1 = dataSet->getReadFromID(source);
		read2 = dataSet->getReadFromID(destination);


		vector<UINT64> *listReadsReverse = new vector<UINT64>;
		vector<UINT16> *listOverlapOffsetsReverse = new vector<UINT16>;
		vector<UINT8> *listOrientationsReverse = new vector<UINT8>;
		UINT64 size = listReads->size();

		UINT64 length1, length2, overlapOffsetForward, revereseOverlap;

		for(UINT64 j = 0; j < size; j++)	// creat the list of reads, overlap offsets and orientations for the reverse edge based on the number in the forward edge.
		{
			listReadsReverse->push_back(listReads->at(size-j-1));

			if(j == 0) // last/first read
			{
				length1 = read2->getReadLength();
				overlapOffsetForward = overlapOffset - length;
			}
			else // any read within the edge
			{
				length1 = dataSet->getReadFromID(listReads->at(size-j))->getReadLength();
				overlapOffsetForward = listOverlapOffsets->at(size-j);
			}
			length2 = dataSet->getReadFromID(listReads->at(size-j-1))->getReadLength();
			revereseOverlap = length1 + overlapOffsetForward - length2;
			listOverlapOffsetsReverse->push_back(revereseOverlap);

			listOrientationsReverse->push_back(!(listOrientations->at(size-j-1)));
		}
		UINT64 revereseOverlapOffset =  overlapOffset + read2->getReadLength() - read1->getReadLength();

		edgeForward->makeEdge(read1, read2, orientation, overlapOffset, listReads, listOverlapOffsets, listOrientations); // make the forward edge.
		edgeReverse->makeEdge(read2, read1, twinEdgeOrientation(orientation), revereseOverlapOffset, listReadsReverse, listOverlapOffsetsReverse, listOrientationsReverse);	// make the reverse edge.

		edgeForward->setReverseEdge(edgeReverse); // set the reverse edge pointer
		edgeReverse->setReverseEdge(edgeForward); // set the reverse edge pinter.

		insertEdge(edgeForward);	// insert the forward edge
		insertEdge(edgeReverse);	// insert the reverse edge.


		i += numberOfReadsInEdge * 3;
	}
	delete list;
	CLOCKSTOP;
	return true;
}




/**********************************************************************************************************************
	Remove dead ends containing less than threshold number of reads.
**********************************************************************************************************************/
/*UINT64 OverlapGraph::removeDeadEnds(void)
{
	CLOCKSTART;
	UINT64 counter = 0;
	for(UINT64 i = 0; i < graph->size(); i++)	// For each node in the graph.
	{
		if(graph->at(i)->size() == 1) 	// If the node has only one edge.
		{
			Edge *edge = graph->at(i)->at(0); // The only edge of the current node.
			if(edge->getListOfReads()->size() < deadEndLength) // If it has fewer reads in it.
			{
				removeEdge(edge);	// We then remove the edge/
				counter++;	// Count the number of such dead-end removed.
			}
		}
	}
	cout << setw(10) << counter << " Dead ends removed." << endl;
	CLOCKSTOP;
	return counter;
}*/




/**********************************************************************************************************************
	Calculate min cost flow
**********************************************************************************************************************/
bool OverlapGraph::calculateFlow(string inputFileName, string outputFileName)
{
	CLOCKSTART;
	UINT64 V = numberOfNodes * 2 + 2, E = numberOfEdges * 3 + numberOfNodes * 4 + 1 , SUPERSOURCE = 1, SUPERSINK = V;
	// For each edge in the overlap graph, we add 3 edges in the directed graph. Two nodes are created for each node in the original graph.
	// A super source and a super sink is added. Each node is connected to the super source and super sink.
	INT64 FLOWLB[3], FLOWUB[3], COST[3];	// Flow bounds and cost of the edges.
	ofstream outputFile;
	outputFile.open(inputFileName.c_str());
	if(!outputFile.is_open())
		MYEXIT("Unable to open file: "+inputFileName);
	stringstream ss;
	ss << "p min " << setw(10) << V << " " << setw(10) << E << endl;  	// Number of nodes and edges in the new graph.
	ss << "n " << setw(10) << SUPERSOURCE << setw(10) << " 0" << endl;	// Flow in the super source
	ss << "n " << setw(10) << SUPERSINK << setw(10) << " 0" << endl;	// Flow in the super sink.
	FLOWLB[0] = 1; FLOWUB[0] = 1000000; COST[0] = 1000000;
	ss << "a " << setw(10) << SUPERSINK << " " << setw(10) << SUPERSOURCE << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl; // Add an edge from super sink to super source with very high cost.
	vector<UINT64> *listOfNodes = new vector<UINT64>;
	vector<UINT64> *listOfNodesReverse = new vector<UINT64>;

	for(UINT64 i = 0; i <= graph->size(); i++)	// For n nodes in the graph, CS2 requires that the nodes are numbered from 1 to n. In the overlap graph, the nodes does not have sequencinal ID. We need to convert them to 1 - n
	{
		listOfNodes->push_back(0);
		listOfNodesReverse->push_back(0);
	// If the ID of a node in the original graph is 100 and directed graph is 5
	// Then listOfNodes->at(100) is equal to 5
	// and ListOfNodesReverse->at(5) is equal to 100.
	}

	// This loop set lower bound and upper bound from super source to super sink. All costs are 0.
	UINT64 currentIndex = 1;
	for(UINT64 i = 1; i < graph->size(); i++)
	{
		if(!graph->at(i)->empty()) // edges to and from the super source and super sink
		{
			listOfNodes->at(i) = currentIndex;	// Mapping between original node ID and cs2 node ID
			listOfNodesReverse->at(currentIndex) = i;	// Mapping between original node ID and cs2 node ID
			FLOWLB[0] = 0; FLOWUB[0] = 1000000; COST[0] = 0;
			ss << "a " << setw(10) << SUPERSOURCE << " " << setw(10) << 2 * currentIndex << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
			ss << "a " << setw(10) << SUPERSOURCE << " " << setw(10) << 2 * currentIndex + 1 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
			ss << "a " << setw(10) << 2 * currentIndex << " " << setw(10) << SUPERSINK << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
			ss << "a " << setw(10) << 2 * currentIndex + 1 << " " << setw(10) << SUPERSINK << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
			currentIndex++;
		}
	}

	// This loop converts the original bi-directed edges to directed edges (1 becomes 6).
	for(UINT64 i = 1; i < graph->size(); i++)
	{
		if(!graph->at(i)->empty()) // edges to and from the super source and super sink
		{
			for(UINT64 j = 0; j < graph->at(i)->size(); j++)
			{
				Edge *edge = graph->at(i)->at(j);
				UINT64 u = listOfNodes->at(edge->getSourceRead()->getReadNumber());
				UINT64 v = listOfNodes->at(edge->getDestinationRead()->getReadNumber());

	// set the bound and cost here
	// if edge has more than 20 reads:
	//   FLOWLB[0] = 1; FLOWUB[0] = 1; COST[0] = 1;
	//   FLOWLB[1] = 0; FLOWUB[1] = 1; COST[1] = 50000;
	//   FLOWLB[2] = 0; FLOWUB[2] = 8; COST[2] = 100000;
	// else:
	//   FLOWLB[0] = 0; FLOWUB[0] = 1; COST[0] = 1;
	//   FLOWLB[1] = 0; FLOWUB[1] = 1; COST[1] = 50000;
	//   FLOWLB[2] = 0; FLOWUB[2] = 8; COST[2] = 100000;
				calculateBoundAndCost(edge, FLOWLB, FLOWUB, COST);

				if(u < v || (u == v && edge < edge->getReverseEdge()))
				{
	// Here for each edge we add three edges with different values of cost and bounds.
	// Total 6 edges considering the reverse edge too.
	// For details on how to convert the edges off different types please see my thesis.
					UINT64 u1 = 2 * u, u2 = 2 * u + 1, v1 =  2 * v, v2 = 2 * v + 1;
					if(edge->getOrientation() == 0)
					{
						ss << "a " << setw(10) << v1 << " " << setw(10) << u1 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
						ss << "a " << setw(10) << u2 << " " << setw(10) << v2 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;

						ss << "a " << setw(10) << v1 << " " << setw(10) << u1 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;
						ss << "a " << setw(10) << u2 << " " << setw(10) << v2 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;

						ss << "a " << setw(10) << v1 << " " << setw(10) << u1 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;
						ss << "a " << setw(10) << u2 << " " << setw(10) << v2 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;
					}
					else if(edge->getOrientation() == 1)
					{
						ss << "a " << setw(10) << v2 << " " << setw(10) << u1 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
						ss << "a " << setw(10) << u2 << " " << setw(10) << v1 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;

						ss << "a " << setw(10) << v2 << " " << setw(10) << u1 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;
						ss << "a " << setw(10) << u2 << " " << setw(10) << v1 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;

						ss << "a " << setw(10) << v2 << " " << setw(10) << u1 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;
						ss << "a " << setw(10) << u2 << " " << setw(10) << v1 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;

					}
					else if(edge->getOrientation() == 2)
					{
						ss << "a " << setw(10) << u1 << " " << setw(10) << v2 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
						ss << "a " << setw(10) << v1 << " " << setw(10) << u2 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;

						ss << "a " << setw(10) << u1 << " " << setw(10) << v2 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;
						ss << "a " << setw(10) << v1 << " " << setw(10) << u2 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;

						ss << "a " << setw(10) << u1 << " " << setw(10) << v2 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;
						ss << "a " << setw(10) << v1 << " " << setw(10) << u2 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;

					}
					else if(edge->getOrientation() == 3)
					{
						ss << "a " << setw(10) << u1 << " " << setw(10) << v1 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;
						ss << "a " << setw(10) << v2 << " " << setw(10) << u2 << " " << setw(10) << FLOWLB[0] << " " << setw(10) << FLOWUB[0] << " " << setw(10) << COST[0] << endl;

						ss << "a " << setw(10) << u1 << " " << setw(10) << v1 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;
						ss << "a " << setw(10) << v2 << " " << setw(10) << u2 << " " << setw(10) << FLOWLB[1] << " " << setw(10) << FLOWUB[1] << " " << setw(10) << COST[1] << endl;

						ss << "a " << setw(10) << u1 << " " << setw(10) << v1 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;
						ss << "a " << setw(10) << v2 << " " << setw(10) << u2 << " " << setw(10) << FLOWLB[2] << " " << setw(10) << FLOWUB[2] << " " << setw(10) << COST[2] << endl;

					}
				}
			}
		}
	}
	outputFile << ss.str();	// Write the string in a file for CS2
	outputFile.close();
	ss.str(std::string());

	char * inFile = new char[inputFileName.size() + 1];
	std::copy(inputFileName.begin(), inputFileName.end(), inFile);
	inFile[inputFileName.size()] = '\0';

	char * outFile = new char[outputFileName.size() + 1];
	std::copy(outputFileName.begin(), outputFileName.end(), outFile);
	outFile[outputFileName.size()] = '\0';


	FILE_LOG(logINFO) << "Calling CS2";
	main_cs2(inFile,outFile);	// Call CS2
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
		inputFile >> source >> destination >> flow;	// get the flow from CS2

		if(source != SUPERSINK && source != SUPERSOURCE && destination != SUPERSOURCE && destination != SUPERSINK && flow!=0)
		{
			UINT64 mySource = listOfNodesReverse->at(source/2);	// Map the source to the original graph
			UINT64 myDestination = listOfNodesReverse->at(destination/2);	// Map the destination in the original graph
			Edge *edge = findEdge(mySource, myDestination);	// Find the edge in the original graph.
			edge->flow += flow;	// Add the flow in the original graph.
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
		if(graph->at(source)->at(i)->getDestinationRead()->getReadNumber() == destination)	// check if there is an edge to destination
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
	for(UINT64 i = 0; i < graph->at(source)->size(); i++)	// flro the list of edges of the source node
	{
		if(graph->at(source)->at(i)->getDestinationRead()->getReadNumber() == destination)	// check if there is an edge to destination
			return true;	// return true if there is an edge (source,destination)
	}
	return false;	// edge not found between source and destination
}


/**********************************************************************************************************************
	This function calculates the cost and bounds for an edge in the overlap graph.
	This function is very sensitive to the assembled contigs.
**********************************************************************************************************************/
bool OverlapGraph::calculateBoundAndCost(Edge *edge, INT64* FLOWLB, INT64* FLOWUB, INT64* COST)
{
	for(UINT64 i = 0; i < 3; i++)	// For the simple edges we put very high cost
	{
		FLOWLB[i] = 0; FLOWUB[i] = 10; COST[i] = 500000;
	}

	if(!edge->getListOfReads()->empty()) // Composite Edge
	{
		if(edge->getListOfReads()->size() > 20) // Composite Edge of at least 20 reads. Must have at least one unit of flow.
		{
			FLOWLB[0] = 1; FLOWUB[0] = 1; COST[0] = 1;
			FLOWLB[1] = 0; FLOWUB[1] = 1; COST[1] = 50000;
			FLOWLB[2] = 0; FLOWUB[2] = 8; COST[2] = 100000;
		}
		else // Short composite edge containing less than 20 reads. May have zero flow.
		{
			FLOWLB[0] = 0; FLOWUB[0] = 1; COST[0] = 1;
			FLOWLB[1] = 0; FLOWUB[1] = 1; COST[1] = 50000;
			FLOWLB[2] = 0; FLOWUB[2] = 8; COST[2] = 100000;
		}
	}

	return true;
}


/**********************************************************************************************************************
	This function finds all the paths between two reads in the overlap graph.
**********************************************************************************************************************/

bool OverlapGraph::findPathBetweenMatepairs(Read * read1, Read * read2, UINT8 orient, UINT8 datasetNumber, vector <Edge *> &copyOfPath, vector <UINT64> &copyOfFlags)
{
	UINT64 pathFound = 0;
	vector <Edge *> firstPath;
	vector <UINT64> flags;
	copyOfPath.clear();
	copyOfFlags.clear();
	//UINT64 flag = 0;

	vector<Edge *> *listRead1, *listRead2;
	vector<UINT64> *locationOnEdgeRead1, *locationOnEdgeRead2;

	listRead1 = (orient == 2 || orient == 3) ? read1->getListOfEdgesForward() : read1->getListOfEdgesReverse();
	locationOnEdgeRead1 = (orient == 2 || orient == 3) ? read1->getLocationOnEdgeForward() : read1->getLocationOnEdgeReverse();

	//listRead2 = (orient == 1 || orient == 3) ? read2->getListOfEdgesForward() : read2->getListOfEdgesReverse(); // if the matepairs are forward - forward
	//locationOnEdgeRead2 = (orient == 1 || orient == 3) ? read2->getLocationOnEdgeForward() : read2->getLocationOnEdgeReverse(); // if the matepairs are forward - forward


	listRead2 = (orient == 0 || orient == 2) ? read2->getListOfEdgesForward() : read2->getListOfEdgesReverse(); // if the matepairs are forward - reverse
	locationOnEdgeRead2 = (orient == 0 || orient == 2) ? read2->getLocationOnEdgeForward() : read2->getLocationOnEdgeReverse(); // if the matepairs are forward - reverse

	if(listRead1->size()==0 || listRead2->size()==0)	// Not  present in any edges.
	{
		return false;
	}
	else
	{
		for(UINT64 i = 0; i < listRead1->size(); i++)
		{
			for(UINT64 j = 0; j < listRead2->size(); j++)
			{
				Edge * firstEdge = listRead1->at(i);
				Edge * lastEdge = listRead2->at(j);
				if(firstEdge ==lastEdge || firstEdge == lastEdge->getReverseEdge()) // Both are on the same edge
				{
					return false;
				}
			}
		}
	}

	for(UINT64 i = 0; i < listRead1->size(); i++)	// If not on the same edge.
	{
		for(UINT64 j = 0; j < listRead2->size(); j++)
		{
			Edge * firstEdge = listRead1->at(i);
			Edge * lastEdge = listRead2->at(j);
			if(firstEdge !=lastEdge && firstEdge != lastEdge->getReverseEdge()) //not on the same edge
			{
				UINT64 distanceOnFirstEdge = firstEdge->getOverlapOffset() - locationOnEdgeRead1->at(i);
				UINT64 distanceOnLastEdge = locationOnEdgeRead2->at(j);
				if(distanceOnFirstEdge + distanceOnLastEdge < getMean(datasetNumber) + 3 * getSD(datasetNumber))
				{
					UINT64 newPaths= exploreGraph(firstEdge, lastEdge, distanceOnFirstEdge, distanceOnLastEdge, datasetNumber, 0, firstPath, flags);	// from firs edge  try to find a path to the last edge.
					if(newPaths > 0)	// If some path found.
					{
						pathFound+=newPaths;	// How many paths found.
						if(copyOfPath.empty())	// Path found for the first time.
						{
							for(UINT64 k = 0; k < firstPath.size(); k++)	// Add the paths in the list.
								copyOfPath.push_back(firstPath.at(k));
							for(UINT64 k = 0; k < firstPath.size() - 1; k++)	// Also add the flag.
								copyOfFlags.push_back(flags.at(k));
						}
						else	// Not the first path
						{
							UINT64 k , l;
							for( k = 0; k< copyOfPath.size() - 1; k++)	// Check if the previous path contains the same pair of edges adjacent to the new path and has flag 1.
							{
								for( l = 0; l < firstPath.size() - 1; l++)
								{
									if(copyOfPath.at(k) == firstPath.at(l) &&  copyOfPath.at(k+1) == firstPath.at(l+1) && flags.at(l) == 1)
										break;
								}
								if(l == firstPath.size() - 1)	// This pair is not supportd
									copyOfFlags.at(k) = 0;
							}
						}
					}
				}
			}
		}
	}
	return true;
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


/**********************************************************************************************************************
	Explore all paths starting from firstEdge and tries to reach lastEdge using depth first search.
	Depth is limited to 10 to reduce running time
**********************************************************************************************************************/

UINT64 OverlapGraph::exploreGraph(Edge* firstEdge, Edge * lastEdge, UINT64 distanceOnFirstEdge, UINT64 distanceOnLastEdge, UINT64 datasetNumber, UINT64 level, vector <Edge *> &firstPath, vector <UINT64> &flags)
{
	static UINT64 pathFound;
	static vector <Edge *> listOfEdges;
	static vector <UINT64> pathLengths;
	if(level == 0)
	{
		pathFound = 0;
		firstPath.resize(0);
		flags.resize(0);
		listOfEdges.resize(0);
		pathLengths.resize(0);
	}
	else
	{
		listOfEdges.resize(level);
		pathLengths.resize(level);
	}

	if(level > 100) return 0; // Do not go very deep.


	if(level == 0)
	{
		listOfEdges.push_back(firstEdge);
		pathLengths.push_back(distanceOnFirstEdge);
	}
	else
	{
		if(firstEdge == lastEdge) // Destination found.
		{
			if(distanceOnLastEdge + pathLengths.at(level - 1) >= getMean(datasetNumber) - 3 * getSD(datasetNumber) && distanceOnLastEdge + pathLengths.at(level - 1) <= getMean(datasetNumber) + 3 * getSD(datasetNumber))
			{
				listOfEdges.push_back(firstEdge);
				pathLengths.push_back(distanceOnLastEdge + pathLengths.at(level - 1));
				pathFound++;
				if(pathFound == 1)	// First path
				{

					for(UINT64 i = 0; i <listOfEdges.size() ; i++) // Copy the path.
						firstPath.push_back(listOfEdges.at(i));

					for(UINT64 i = 0; i <listOfEdges.size() - 1 ; i++) // All adjacent edges in the path are supported.
						flags.push_back(1);
				}
				else	// Not the first path.
				{
					UINT64 i , j;
					for( i = 0; i< firstPath.size() - 1; i++)	// Compare to the first path.
					{
						for( j = 0; j < listOfEdges.size() - 1; j++)
						{
							if(firstPath.at(i) == listOfEdges.at(j) &&  firstPath.at(i+1) == listOfEdges.at(j+1)) // A pair of edges adjacent in the first path is also adjacent in this path. We keep the support.
								break;
						}
						if(j == listOfEdges.size() - 1)	// A pair of adjacent edge in the first path is not adjacent in this path. So this pair is not supported anymore. So we clear the flag.
							flags.at(i) = 0;
					}
				}
				/*for(UINT64 i = 0; i <firstPath.size() ; i++) // Copy the path.
				{
					if(i <listOfEdges.size() - 1 )
						cout<< " (" << firstPath.at(i)->getSourceRead()->getReadNumber() <<", " <<firstPath.at(i)->getDestinationRead()->getReadNumber()<<") " <<flags.at(i);
					else
						cout<< " (" << firstPath.at(i)->getSourceRead()->getReadNumber() <<", " <<firstPath.at(i)->getDestinationRead()->getReadNumber()<<") ";
				}
				cout << endl;*/
				return 1;
			}
			else	// add the new edge in the path
			{
				listOfEdges.push_back(firstEdge);
				pathLengths.push_back( distanceOnFirstEdge + pathLengths.at(level - 1) );
			}
		}
		else	// add the new edge in the path.
		{
			listOfEdges.push_back(firstEdge);
			pathLengths.push_back( distanceOnFirstEdge + pathLengths.at(level - 1) );
		}
	}

	for(UINT64 i = 0 ; i < graph->at(firstEdge->getDestinationRead()->getReadNumber())->size(); i++ )	// go deepter in the graph to reach destination edge
	{
		Edge * nextEdge =  graph->at(firstEdge->getDestinationRead()->getReadNumber())->at(i);
		if(matchEdgeType(firstEdge, nextEdge) && pathLengths.at(level) < getMean(datasetNumber) + 3 * getSD(datasetNumber))
			exploreGraph(nextEdge, lastEdge, nextEdge->getOverlapOffset(), distanceOnLastEdge, datasetNumber, level+1, firstPath, flags);
	}
	return pathFound;
}


// Please move this into header file
struct pairedEdges
{
		Edge * edge1;
		Edge * edge2;
		UINT64 support;
		UINT64 distance;
		bool isFreed;
		bool operator < (const pairedEdges& rhs) const
		{
		       return support > rhs.support;
		}

};


/**********************************************************************************************************************
	Check for paths in the overlap graph between each matepair and calculate the support by using the paths.
**********************************************************************************************************************/
UINT64 OverlapGraph::findSupportByMatepairsAndMerge(void)
{
	CLOCKSTART;

	// if the file set is not mate-pair, then just skip
	if(meanOfInsertSizes.size() == 0) // no mate-pair
		return 0;
	vector <Edge *> copyOfPath;
	vector <UINT64> copyOfFlags;
	UINT64 noPathsFound = 0, pathsFound = 0, mpOnSameEdge=0;

	vector <pairedEdges> listOfPairedEdges;
	for(UINT64 i = 1; i < graph->size(); i++)
	{
		Read * read1 = dataSet->getReadFromID(i);
	// initial mate-pair information is saved before building the graph using storeMatePairInfo()
		for(UINT64 j = 0; j < read1->getMatePairList()->size(); j++)
		{
			Read * read2 = dataSet->getReadFromID(read1->getMatePairList()->at(j).matePairID);
			if(read1->getReadNumber() > read2->getReadNumber()) 	// To avoid finding path twice
				continue;
			if(meanOfInsertSizes.at(read1->getMatePairList()->at(j).datasetNumber) == 0)
				continue;

			if( findPathBetweenMatepairs(read1, read2, read1->getMatePairList()->at(j).matePairOrientation, read1->getMatePairList()->at(j).datasetNumber, copyOfPath, copyOfFlags) == true)
			{
	// Matepair on different edge
				if(copyOfPath.size() == 0)
					noPathsFound++;
				else
					pathsFound++;
			}
			else // Mate pair on the same edge
			{
				mpOnSameEdge++;
			}

			if(copyOfPath.size() > 1)	// Path found
			{
				for(UINT64 k = 0; k < copyOfFlags.size(); k++)
				{
					if(copyOfFlags.at(k) == 1)	// edge at k and k+1 is supported. We need to add it in the list if not present. If already the pair of edges present then increase the support
					{
						UINT64 l;
						for(l = 0; l < listOfPairedEdges.size(); l++)
						{
							if(listOfPairedEdges.at(l).edge1 == copyOfPath.at(k) && listOfPairedEdges.at(l).edge2 == copyOfPath.at(k+1)) // already present in the list
							{
								listOfPairedEdges.at(l).support = listOfPairedEdges.at(l).support + 1;	// only increase the support
								break;
							}
							else if(listOfPairedEdges.at(l).edge2->getReverseEdge() == copyOfPath.at(k) && listOfPairedEdges.at(l).edge1->getReverseEdge() == copyOfPath.at(k+1)) // already present in the list
							{
								listOfPairedEdges.at(l).support = listOfPairedEdges.at(l).support + 1;	// only increase the support
								break;
							}
						}
						if(l == listOfPairedEdges.size()) // not present in the list
						{
							if(copyOfPath.at(k)->getSourceRead()->getReadNumber() != copyOfPath.at(k)->getDestinationRead()->getReadNumber() || copyOfPath.at(k+1)->getSourceRead()->getReadNumber() != copyOfPath.at(k+1)->getDestinationRead()->getReadNumber()) // add in the list with support 1
	//if(copyOfPath.at(k)!=copyOfPath.at(k+1) && copyOfPath.at(k)!=copyOfPath.at(k+1)->getReverseEdge()) // do not want to add support between edge (a,a) and (a,a)
							{
								pairedEdges newPair;
								newPair.edge1 = copyOfPath.at(k);
								newPair.edge2 = copyOfPath.at(k+1);
								newPair.support = 1;
								newPair.isFreed = false;
								listOfPairedEdges.push_back(newPair);
							}
						}
					}
				}
			}
		}
	}

	sort(listOfPairedEdges.begin(), listOfPairedEdges.end());

	UINT64 pairsOfEdgesMerged = 0;

	for(UINT64 i = 0; i<listOfPairedEdges.size(); i++)
	{
		if(listOfPairedEdges.at(i).isFreed == false && listOfPairedEdges.at(i).support >= minimumSupport)// && listOfPairedEdges.at(i).edge1->flow > 0 && listOfPairedEdges.at(i).edge2->flow > 0)
		{
			pairsOfEdgesMerged++;
			FILE_LOG(logDEBUG) << setw(4) << i + 1 << " Merging (" << setw(10) << listOfPairedEdges.at(i).edge1->getSourceRead()->getReadNumber() << "," << setw(10) <<listOfPairedEdges.at(i).edge1->getDestinationRead()->getReadNumber() << ") Length: " << setw(8) << listOfPairedEdges.at(i).edge1->getOverlapOffset() << " Flow: " << setw(3) << listOfPairedEdges.at(i).edge1->flow << " and (" << setw(10) << listOfPairedEdges.at(i).edge2->getSourceRead()->getReadNumber() << "," << setw(10) << listOfPairedEdges.at(i).edge2->getDestinationRead()->getReadNumber() << ") Length: " << setw(8) << listOfPairedEdges.at(i).edge2->getOverlapOffset() << " Flow: " << setw(3) << listOfPairedEdges.at(i).edge2->flow << " are supported " << setw(4) << listOfPairedEdges.at(i).support<<" times";

			Edge * e1f = listOfPairedEdges.at(i).edge1, *e1r = listOfPairedEdges.at(i).edge1->getReverseEdge();
			Edge * e2f = listOfPairedEdges.at(i).edge2, *e2r = listOfPairedEdges.at(i).edge2->getReverseEdge();

			mergeEdges(listOfPairedEdges.at(i).edge1, listOfPairedEdges.at(i).edge2);

			for(UINT64 j = i + 1; j<listOfPairedEdges.size(); j++)
			{
				if( listOfPairedEdges.at(j).edge1 == e1f || listOfPairedEdges.at(j).edge1 == e1r || listOfPairedEdges.at(j).edge1 == e2f || listOfPairedEdges.at(j).edge1 == e2r )
					listOfPairedEdges.at(j).isFreed = true;
				if( listOfPairedEdges.at(j).edge2 == e1f || listOfPairedEdges.at(j).edge2 == e1r || listOfPairedEdges.at(j).edge2 == e2f || listOfPairedEdges.at(j).edge2 == e2r )
					listOfPairedEdges.at(j).isFreed = true;
			}
		}
	}
	FILE_LOG(logINFO) << pairsOfEdgesMerged <<" Pairs of Edges merged out of " << listOfPairedEdges.size() << " supported pairs of edges" ;
	FILE_LOG(logINFO) << "No paths found between " << noPathsFound << " matepairs that are on different edge.";
	FILE_LOG(logINFO) << "Paths found between " << pathsFound << " matepairs that are on different edge.";
	FILE_LOG(logINFO) << "Total matepairs on different edges " << pathsFound+ noPathsFound;
	FILE_LOG(logINFO) << "Total matepairs on same edge " << mpOnSameEdge;
	FILE_LOG(logINFO) << "Total matepairs " << pathsFound+noPathsFound+mpOnSameEdge;
	CLOCKSTOP;
	return pairsOfEdgesMerged;

}


/**********************************************************************************************************************
	This function returns the string by overlapping the reads in an edge in the overlap graph
**********************************************************************************************************************/

string OverlapGraph::getStringInEdge(Edge *edge)
{
	string read1, read2, readTemp, returnString;
	read1 = (edge->getOrientation() == 2 || edge->getOrientation() == 3) ?  edge->getSourceRead()->getStringForward() : edge->getSourceRead()->getStringReverse();
	read2 = (edge->getOrientation() == 1 || edge->getOrientation() == 3) ?  edge->getDestinationRead()->getStringForward() : edge->getDestinationRead()->getStringReverse();
	returnString = read1;
	UINT64 previousLength = read1.length(), substringLength;
	for(UINT64 i = 0; i < edge->getListOfReads()->size(); i++)
	{
		readTemp = (edge->getListOfOrientations()->at(i) == 1) ? dataSet->getReadFromID(edge->getListOfReads()->at(i))->getStringForward(): dataSet->getReadFromID(edge->getListOfReads()->at(i))->getStringReverse();

		substringLength =  readTemp.length() + edge->getListOfOverlapOffsets()->at(i) - previousLength;
		if( edge->getListOfOverlapOffsets()->at(i) ==  previousLength)
			returnString = returnString + 'N' + readTemp.substr(readTemp.length() - substringLength, substringLength);
		else
			returnString = returnString + readTemp.substr(readTemp.length() - substringLength, substringLength);
		previousLength = readTemp.length();


	}

	if(edge->getListOfReads()->empty()) // Simple edge
	{
		substringLength =  read2.length() + edge->getOverlapOffset() - read1.length();
		returnString = returnString + read2.substr(read2.length() - substringLength, substringLength);
	}
	else
	{
		substringLength = edge->getReverseEdge()->getListOfOverlapOffsets()->at(0);
		returnString = returnString + read2.substr(read2.length() - substringLength, substringLength);
	}
	return returnString;
}


/**********************************************************************************************************************
	This function removes in-trees and out-trees.
**********************************************************************************************************************/

UINT64 OverlapGraph::reduceTrees(void)
{
	CLOCKSTART;
	UINT64 NumOfInEdges, NumOfOutEdges, inFlow, outFlow, nodeMerged = 0;
	vector <Edge *> listOfInEdges, listOfOutEdges;
	for(UINT64 i = 0; i < graph->size(); i++)	// For each node in the graph
	{

		NumOfInEdges = 0; NumOfOutEdges = 0; inFlow = 0; outFlow = 0;
		listOfInEdges.clear(); listOfOutEdges.clear();
		for(UINT64 j = 0; j< graph->at(i)->size(); j++)
		{
			if(graph->at(i)->at(j)->flow == 0 || graph->at(i)->at(j)->flow != graph->at(i)->at(j)->getReverseEdge()->flow || graph->at(i)->at(j)->getSourceRead()->getReadNumber() == graph->at(i)->at(j)->getDestinationRead()->getReadNumber() )  // Some conditions for not being considered as a tree
					break;
			if(graph->at(i)->at(j)->getOrientation() == 0 || graph->at(i)->at(j)->getOrientation() == 1 )	// It is an in-edge
			{
				NumOfInEdges++;
				inFlow += graph->at(i)->at(j)->flow;	// Count the in flow
				listOfInEdges.push_back(graph->at(i)->at(j));
			}
			else	// It is an out-edge
			{
				NumOfOutEdges++;
				outFlow += graph->at(i)->at(j)->flow;	// Count the out flow
				listOfOutEdges.push_back(graph->at(i)->at(j));
			}

			if(inFlow == outFlow && ( (NumOfInEdges == 1 && NumOfOutEdges > 1) || (NumOfInEdges > 1 && NumOfOutEdges == 1) ) )	// Either an in tree or an out tree
			{
				nodeMerged++;
				for(UINT64 k = 0; k < listOfInEdges.size(); k++)
				{
					for(UINT64 l = 0; l < listOfOutEdges.size(); l++)
					{
						mergeEdges(listOfInEdges.at(k)->getReverseEdge(), listOfOutEdges.at(l));
					}
				}
			}
		}
	}
	FILE_LOG(logINFO) << setw(10) << nodeMerged << " trees removed.";
	CLOCKSTOP;
	return nodeMerged;
}



/**********************************************************************************************************************
	removes composite path, dead-ends and trees
**********************************************************************************************************************/
bool OverlapGraph::simplifyGraph(void)
{
	UINT64 counter = 0;
	do
	{
		 //counter  = removeDeadEndNodes();	// Remove dead-ends
		 counter = contractCompositePaths();	// Contract composite paths
		 counter += removeSimilarEdges();
		 //counter += reduceTrees();	// Remove trees.
		 //counter += reduceLoops();

	} while (counter > 0);
	return true;
}


/**********************************************************************************************************************
	Original Scaffolder function.
	*******************************************************************************************************************/
UINT64 OverlapGraph::scaffolder(void)
{
	CLOCKSTART;
	UINT64 pairsOfEdgesMerged = 0, dist;
	UINT8 orient, datasetNumber;
	Read *read1, *read2;
	vector <pairedEdges> listOfPairedEdges;
	vector<Edge *> *listRead1, *listRead2;
	vector<UINT64> *locationOnEdgeRead1, *locationOnEdgeRead2;

	for(UINT64 i = 1; i < graph->size(); i++)	// For each read
	{
		read1 = dataSet->getReadFromID(i);	// Get the read object
		for(UINT64 j = 0; j < read1->getMatePairList()->size(); j++)// For each matepair
		{
			read2 = dataSet->getReadFromID(read1->getMatePairList()->at(j).matePairID);	// Get the read object of the matepair.

			if(read1->getReadNumber() > read2->getReadNumber()) // To avoid duplicate computation
				continue;

			orient = read1->getMatePairList()->at(j).matePairOrientation;	// Get the matepair orientation
			datasetNumber = read1->getMatePairList()->at(j).datasetNumber;	// Get the dataset number

	// 0 = 00 means the reverse of r1 and the reverse of r2 are matepairs.
	// 1 = 01 means the reverse of r1 and the forward of r2 are matepairs.
	// 2 = 10 means the forward of r1 and the reverse of r2 are matepairs.
	// 3 = 11 means the forward of r1 and the forward of r2 are matepairs.
	// To calculate distance of forward read, flip the read and get the location of the offset.

			listRead1 = (orient == 0 || orient == 1) ? read1->getListOfEdgesForward() : read1->getListOfEdgesReverse();
			locationOnEdgeRead1 = (orient == 0 || orient == 1) ? read1->getLocationOnEdgeForward() : read1->getLocationOnEdgeReverse();

	// To calculate distance of reverse read, flip the read and get the location of the offset.
			listRead2 = (orient == 0 || orient == 2) ? read2->getListOfEdgesForward() : read2->getListOfEdgesReverse();
			locationOnEdgeRead2 = (orient == 0 || orient == 2) ? read2->getLocationOnEdgeForward() : read2->getLocationOnEdgeReverse();

	// Only consider uniquely mapped reads and the distance is less than mean+3*SD
			if( listRead1->size() == 1 && listRead2->size() == 1 && locationOnEdgeRead1->at(0) + locationOnEdgeRead2->at(0) < (getMean(datasetNumber) + 3 * getSD(datasetNumber)) )  // Both the reads are present on only on edge and the distance is less that mean+3*sd
			{
				dist = locationOnEdgeRead1->at(0) + locationOnEdgeRead2->at(0);
	// if there are already in the same edge, don't do anything
				if(listRead1->at(0) == listRead2->at(0) ||  listRead1->at(0) == listRead2->at(0)->getReverseEdge()) // Not on the same edge
					continue;

				UINT64 k;

				for(k = 0; k < listOfPairedEdges.size(); k++)
				{
					if(listOfPairedEdges.at(k).edge1 == listRead1->at(0)->getReverseEdge() && listOfPairedEdges.at(k).edge2 == listRead2->at(0))	// This pair of edge was supported previously by another mate-pair. Only increase the support.
					{
						listOfPairedEdges.at(k).support += 1;
						listOfPairedEdges.at(k).distance += dist;
						break;
					}
					if(listOfPairedEdges.at(k).edge1 == listRead2->at(0)->getReverseEdge() && listOfPairedEdges.at(k).edge2 == listRead1->at(0))	// This pair of edge was supported previously by another mate-pair. Only increase the support.
					{
						listOfPairedEdges.at(k).support += 1;
						listOfPairedEdges.at(k).distance += dist;
						break;
					}
				}

				if(k == listOfPairedEdges.size())	// This mate pair is the first reads to support the pair of edges. Add the new pair of edges in the list with support 1.
				{
					pairedEdges newPair;
					newPair.edge1 = listRead1->at(0)->getReverseEdge();
					newPair.edge2 = listRead2->at(0);
					newPair.support = 1;
					newPair.distance = dist;
					newPair.isFreed = false;
					listOfPairedEdges.push_back(newPair);
				}
			}

		}
	}

	sort(listOfPairedEdges.begin(), listOfPairedEdges.end());	// Sort the list according to support.

	for(UINT64 i = 0; i < listOfPairedEdges.size(); i++)
	{
		if(listOfPairedEdges.at(i).isFreed == false && listOfPairedEdges.at(i).support >= minimumSupport)
		{
			pairsOfEdgesMerged++;
			listOfPairedEdges.at(i).distance = listOfPairedEdges.at(i).distance / listOfPairedEdges.at(i).support;
			FILE_LOG(logDEBUG) << setw(4) << i + 1 << " (" << setw(10) << listOfPairedEdges.at(i).edge1->getSourceRead()->getReadNumber() << "," << setw(10) <<listOfPairedEdges.at(i).edge1->getDestinationRead()->getReadNumber() << ") Length: " << setw(8) << listOfPairedEdges.at(i).edge1->getOverlapOffset() << " Flow: " << setw(3) << listOfPairedEdges.at(i).edge1->flow << " and (" << setw(10) << listOfPairedEdges.at(i).edge2->getSourceRead()->getReadNumber() << "," << setw(10) << listOfPairedEdges.at(i).edge2->getDestinationRead()->getReadNumber() << ") Length: " << setw(8) << listOfPairedEdges.at(i).edge2->getOverlapOffset() << " Flow: " << setw(3) << listOfPairedEdges.at(i).edge2->flow << " are supported " << setw(4) << listOfPairedEdges.at(i).support << " times. Average distance: "<< setw(4) << listOfPairedEdges.at(i).distance;
			Edge * e1f = listOfPairedEdges.at(i).edge1, *e1r = listOfPairedEdges.at(i).edge1->getReverseEdge();
			Edge * e2f = listOfPairedEdges.at(i).edge2, *e2r = listOfPairedEdges.at(i).edge2->getReverseEdge();

			mergeEdgesDisconnected(listOfPairedEdges.at(i).edge1, listOfPairedEdges.at(i).edge2,listOfPairedEdges.at(i).distance);	// Merge the edges.

			for(UINT64 j = i + 1; j<listOfPairedEdges.size(); j++)
			{
				if( listOfPairedEdges.at(j).edge1 == e1f || listOfPairedEdges.at(j).edge1 == e1r || listOfPairedEdges.at(j).edge1 == e2f || listOfPairedEdges.at(j).edge1 == e2r )
					listOfPairedEdges.at(j).isFreed = true;
				if( listOfPairedEdges.at(j).edge2 == e1f || listOfPairedEdges.at(j).edge2 == e1r || listOfPairedEdges.at(j).edge2 == e2f || listOfPairedEdges.at(j).edge2 == e2r )
					listOfPairedEdges.at(j).isFreed = true;
			}
		}
	}

	CLOCKSTOP;
	return pairsOfEdgesMerged;
}

/**********************************************************************************************************************
	New Scaffolder function.
**********************************************************************************************************************/
/*
struct keyInfo
{
    Edge *edge1;   // key1 -> *edge1
    Edge *edge2;   // key2 -> *edge2

    keyInfo(Edge A, Edge B) :
        edge1(&A), edge2(&B) {}

    bool operator<(const keyInfo& A) const
    {
        return ((edge1<A.edge1) || (edge2<A.edge2)); 
    }

};

struct valueInfo
{
    UINT64 support;
    UINT64 distance;
    bool isFreed;

    valueInfo(UINT64 A,UINT64 B,bool C) : 
    support(A),distance(B),isFreed(C) {}
};

typedef map<keyInfo, valueInfo> MapType;


UINT64 OverlapGraph::scaffolder(void)
{
	CLOCKSTART;
	UINT64 pairsOfEdgesMerged = 0, dist;
	UINT8 orient, datasetNumber;
	Read *read1, *read2;
	vector <pairedEdges> listOfPairedEdges;
    // new map for struct key and val instead of vector
	MapType mapOfPairedEdges;
	vector<Edge *> *listRead1, *listRead2;
	vector<UINT64> *locationOnEdgeRead1, *locationOnEdgeRead2;

	for(UINT64 i = 1; i < graph->size(); i++)	// For each read
	{
		read1 = dataSet->getReadFromID(i);	// Get the read object
		for(UINT64 j = 0; j < read1->getMatePairList()->size(); j++)// For each matepair
		{
			read2 = dataSet->getReadFromID(read1->getMatePairList()->at(j).matePairID);	// Get the read object of the matepair.

			if(read1->getReadNumber() > read2->getReadNumber()) // To avoid duplicate computation
				continue;

			orient = read1->getMatePairList()->at(j).matePairOrientation;	// Get the matepair orientation
			datasetNumber = read1->getMatePairList()->at(j).datasetNumber;	// Get the dataset number

	// 0 = 00 means the reverse of r1 and the reverse of r2 are matepairs.
	// 1 = 01 means the reverse of r1 and the forward of r2 are matepairs.
	// 2 = 10 means the forward of r1 and the reverse of r2 are matepairs.
	// 3 = 11 means the forward of r1 and the forward of r2 are matepairs.
	// To calculate distance of forward read, flip the read and get the location of the offset.
			listRead1 = (orient == 0 || orient == 1) ? read1->getListOfEdgesForward() : read1->getListOfEdgesReverse();
			locationOnEdgeRead1 = (orient == 0 || orient == 1) ? read1->getLocationOnEdgeForward() : read1->getLocationOnEdgeReverse();
	// To calculate distance of reverse read, flip the read and get the location of the offset.
			listRead2 = (orient == 0 || orient == 2) ? read2->getListOfEdgesForward() : read2->getListOfEdgesReverse();
			locationOnEdgeRead2 = (orient == 0 || orient == 2) ? read2->getLocationOnEdgeForward() : read2->getLocationOnEdgeReverse();

	// get distance
			dist = locationOnEdgeRead1->at(0) + locationOnEdgeRead2->at(0);

	// Only consider uniquely mapped reads and the distance is less than mean+3*SD
			if( listRead1->size() == 1 && listRead2->size() == 1 && dist < (getMean(datasetNumber) + 3 * getSD(datasetNumber)) )  // Both the reads are present on only on edge and the distance is less that mean+3*sd
			{
	// if there are already in the same edge, don't do anything
				if(listRead1->at(0) == listRead2->at(0) ||  listRead1->at(0) == listRead2->at(0)->getReverseEdge()) // Not on the same edge
					continue;

                // get key set (*edge1, *edge2)
                keyInfo keySet1(*listRead1->at(0)->getReverseEdge(), *listRead2->at(0));
                keyInfo keySet2(*listRead2->at(0)->getReverseEdge(), *listRead1->at(0));

                // find the keySet1 and keySet2
                MapType::iterator iter1=mapOfPairedEdges.find(keySet1);
                MapType::iterator iter2=mapOfPairedEdges.find(keySet2);

                // check keySet1
                if (iter1 != mapOfPairedEdges.end())        // if the keySet1 exists
                {
                    iter1->second.support += 1;             // add +1 the support
                    iter1->second.distance += dist;         // add +dist the distance
                }
                else if (iter2 != mapOfPairedEdges.end())   // if the keySet2 exists
                {
                    iter2->second.support += 1;             // add +1 the support
                    iter2->second.distance += dist;         // add +dist the distance
                }
                else                                        // if the keySet1 and keySet2 do not exist
                {
                    valueInfo valueSet(1, dist, false);     // get valueSet (support, distance, isFreed)
                    mapOfPairedEdges.insert(MapType::value_type(keySet1,valueSet)); // insert valueSet to the map
                }
			}
		}
	}

    


	sort(listOfPairedEdges.begin(), listOfPairedEdges.end());	// Sort the list according to support.

	for(UINT64 i = 0; i < listOfPairedEdges.size(); i++)
	{
		if(listOfPairedEdges.at(i).isFreed == false && listOfPairedEdges.at(i).support >= minimumSupport)
		{
			pairsOfEdgesMerged++;
			listOfPairedEdges.at(i).distance = listOfPairedEdges.at(i).distance / listOfPairedEdges.at(i).support;
			cout << setw(4) << i + 1 << " (" << setw(10) << listOfPairedEdges.at(i).edge1->getSourceRead()->getReadNumber() << "," << setw(10) <<listOfPairedEdges.at(i).edge1->getDestinationRead()->getReadNumber() << ") Length: " << setw(8) << listOfPairedEdges.at(i).edge1->getOverlapOffset() << " Flow: " << setw(3) << listOfPairedEdges.at(i).edge1->flow << " and (" << setw(10) << listOfPairedEdges.at(i).edge2->getSourceRead()->getReadNumber() << "," << setw(10) << listOfPairedEdges.at(i).edge2->getDestinationRead()->getReadNumber() << ") Length: " << setw(8) << listOfPairedEdges.at(i).edge2->getOverlapOffset() << " Flow: " << setw(3) << listOfPairedEdges.at(i).edge2->flow << " are supported " << setw(4) << listOfPairedEdges.at(i).support << " times. Average distance: "<< setw(4) << listOfPairedEdges.at(i).distance << endl;
			Edge * e1f = listOfPairedEdges.at(i).edge1, *e1r = listOfPairedEdges.at(i).edge1->getReverseEdge();
			Edge * e2f = listOfPairedEdges.at(i).edge2, *e2r = listOfPairedEdges.at(i).edge2->getReverseEdge();

			mergeEdgesDisconnected(listOfPairedEdges.at(i).edge1, listOfPairedEdges.at(i).edge2,listOfPairedEdges.at(i).distance);	// Merge the edges.

			for(UINT64 j = i + 1; j<listOfPairedEdges.size(); j++)
			{
				if( listOfPairedEdges.at(j).edge1 == e1f || listOfPairedEdges.at(j).edge1 == e1r || listOfPairedEdges.at(j).edge1 == e2f || listOfPairedEdges.at(j).edge1 == e2r )
					listOfPairedEdges.at(j).isFreed = true;
				if( listOfPairedEdges.at(j).edge2 == e1f || listOfPairedEdges.at(j).edge2 == e1r || listOfPairedEdges.at(j).edge2 == e2f || listOfPairedEdges.at(j).edge2 == e2r )
					listOfPairedEdges.at(j).isFreed = true;
			}
		}
	}


	CLOCKSTOP;
	return pairsOfEdgesMerged;
}
*/

/**********************************************************************************************************************
	Check if two strings overlap. At least 10 bp must overlap.
**********************************************************************************************************************/

UINT64 OverlapGraph::findOverlap(string string1, string string2)
{
	UINT64 minimum = min(string1.length(), string2.length());
	for(UINT64 i = minimum - 1; i >= 10; i--)
	{
		if(string1.substr(string1.length() - i, i) == string2.substr(0,i))
		{
			return i;
		}
	}
	return 0;
}



/**********************************************************************************************************************
	Merge two edges that do not share any node.
**********************************************************************************************************************/
bool OverlapGraph::mergeEdgesDisconnected(Edge *edge1, Edge *edge2, UINT64 gapLength)
{
	if(edge1->getDestinationRead()->getReadNumber() == edge2->getSourceRead()->getReadNumber() && matchEdgeType (edge1, edge2)) // If the two edges are already connected. A --->B and B---->C. They share a common node
	{
		mergeEdges(edge1,edge2); // Merge the edges.
		return true;
	}

	// A------>B and C------>D. We need to check if the nodes B and C overlaps or not
	string string1, string2;
	string1 = ( edge1->getOrientation() == 1 || edge1->getOrientation() == 3 ) ? edge1->getDestinationRead()->getStringForward() : edge1->getDestinationRead()->getStringReverse(); // We will check if the two nodes overlap or not
	string2 = ( edge2->getOrientation() == 2 || edge2->getOrientation() == 3 ) ? edge2->getSourceRead()->getStringForward() : edge2->getSourceRead()->getStringReverse();

	UINT64 overlapLength = findOverlap(string1,string2); // Find the overlap between B and C. If they do not overlap then the return will be zero. We check for at least 10 bp overlap

	UINT64 overlapOffset1, overlapOffset2;
	Edge *edgeForward = new Edge();	// Forward edge
	Read *read1 = edge1->getSourceRead(), *read2 = edge2->getDestinationRead(); // Get the reads.
	Edge *edgeReverse = new Edge();	// Reverse edge.
	UINT8 orientationForward = mergedEdgeOrientationDisconnected(edge1,edge2); // Orientation of the forward edge based on the orientations of edge1 and edge2
	UINT8 orientationReverse = twinEdgeOrientation(orientationForward);	// Orientation of the reverse edge.

	vector<UINT64> * listReadsForward = new vector<UINT64>;	// List of reads in the forward edge.
	vector<UINT16> * listOverlapOffsetsForward= new vector<UINT16>;	// List of overlap offsets in the reads of the forward edge.
	vector<UINT8> * listOrientationsForward = new vector<UINT8>;	// List of orientations of the reads of the forward edge. 1 means forward string of the reads, 0 means reverse string of the read

	if(overlapLength == 0) // Strings in the read B and C do not overlap
	{
		overlapOffset1 = edge1->getDestinationRead()->getReadLength();	// In this case we concatenate the strings in the edges. So the offset is the length of the read B
		overlapOffset2 = edge2->getSourceRead()->getReadLength();	// This is the overlap offset of the reverse edge.
	}
	else	// Strings in the read B and C do overlap
	{
		overlapOffset1 = edge1->getDestinationRead()->getReadLength() - overlapLength; // Overlap offset of the forward edge is taken according to the length of the string in B
		overlapOffset2 = edge2->getSourceRead()->getReadLength() - overlapLength; // overlap offset of the reverse edge is taken according to the length of the string in C
	}


	mergeListDisconnected(edge1, edge2, overlapOffset1, gapLength, listReadsForward, listOverlapOffsetsForward, listOrientationsForward);	// Merge the list of reads, overlaps etc for the forward edge
	edgeForward->makeEdge(read1,read2,orientationForward, edge1->getOverlapOffset() + edge2->getOverlapOffset() + overlapOffset1, listReadsForward, listOverlapOffsetsForward, listOrientationsForward);

	vector<UINT64> * listReadsReverse = new vector<UINT64>;	// List of reads in the reverse edge.
	vector<UINT16> * listOverlapOffsetsReverse= new vector<UINT16>;	// List of overlap offsets in the reads of the reverse edge.
	vector<UINT8> * listOrientationsReverse = new vector<UINT8>;	// List of orientations of the reads of the reverse edge. 1 means forward string of the reads, 0 means reverse string of the read

	mergeListDisconnected(edge2->getReverseEdge(),edge1->getReverseEdge(), overlapOffset2, gapLength, listReadsReverse, listOverlapOffsetsReverse,listOrientationsReverse); // Merge the list of reads, overlaps etc for the reverse edge
	edgeReverse->makeEdge(read2, read1, orientationReverse, edge1->getReverseEdge()->getOverlapOffset() + edge2->getReverseEdge()->getOverlapOffset() + overlapOffset2, listReadsReverse, listOverlapOffsetsReverse, listOrientationsReverse);

	edgeForward->setReverseEdge(edgeReverse);	// set the pointer of reverse edge
	edgeReverse->setReverseEdge(edgeForward);	// set the pointer of reverse edge

	UINT16 flow = min(edge1->flow,edge2->flow);	// Take the minimum of the flow from the two original edges.
	UINT64 coverage = min(edge1->coverageDepth, edge2->coverageDepth);	// not used
	edgeForward->flow = flow;	// set the flow in the forward edge.
	edgeForward->coverageDepth = coverage;

	edgeReverse->flow = flow;	// set the flow in the reverse edge.
	edgeReverse->coverageDepth = coverage;

	//if(flowComputed && flow == 0 && edgeForward->getOverlapOffset() > 1000)
	//{
	//	cout << "Check for error inserting edge between " << edgeForward->getSourceRead()->getReadNumber() << " and " << edgeForward->getDestinationRead()->getReadNumber() << " Length: " << edgeForward->getOverlapOffset() << endl;
	//}
	insertEdge(edgeForward); // insert forward the edge in the graph.
	insertEdge(edgeReverse); // insert the reverse edge in the graph.

	edge1->flow = edge1->flow - flow;	// decrease the flow in the original edge.
	edge1->getReverseEdge()->flow = edge1->getReverseEdge()->flow - flow; // decrease the flow in the original edge.
	edge1->coverageDepth = edge1->coverageDepth - coverage;
	edge1->getReverseEdge()->coverageDepth = edge1->getReverseEdge()->coverageDepth - coverage;

	edge2->flow = edge2->flow - flow; // decrease the flow in the original edge.
	edge2->getReverseEdge()->flow = edge2->getReverseEdge()->flow - flow; // decrease the flow in the original edge.
	edge2->coverageDepth = edge2->coverageDepth - coverage;
	edge2->getReverseEdge()->coverageDepth = edge2->getReverseEdge()->coverageDepth - coverage;

	if(edge1->flow == 0 || flow == 0)	// Remove the edge1 if the flow is used.
		removeEdge(edge1);
	if(edge2->flow == 0 || flow == 0)	// Remove the edge2 if the flow is used.
		removeEdge(edge2);

	return true;
}



/**********************************************************************************************************************
	Merge the list of reads, list of overlap offsets and list of orientations of two edges.
**********************************************************************************************************************/

bool OverlapGraph::mergeListDisconnected(Edge *edge1, Edge *edge2, UINT64 overlapOffset, UINT64 gapLength, vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets, vector<UINT8> * listOrientations)
{
	UINT64 sum = 0;
	for(UINT64 i = 0; i < edge1->getListOfOrientations()->size(); i++)	// Get the list from the first edge.
	{
		listReads->push_back(edge1->getListOfReads()->at(i));
		listOverlapOffsets->push_back(edge1->getListOfOverlapOffsets()->at(i));
		listOrientations->push_back(edge1->getListOfOrientations()->at(i));
		sum += edge1->getListOfOverlapOffsets()->at(i);
	}
	listReads->push_back(edge1->getDestinationRead()->getReadNumber());	// Add the destination read of the first edge.

		listOverlapOffsets->push_back(edge1->getOverlapOffset() - sum);

	if(edge1->getOrientation() == 1 || edge1->getOrientation() == 3)
		listOrientations->push_back(1);
	else
		listOrientations->push_back(0);

	listReads->push_back(edge2->getSourceRead()->getReadNumber());	// Add the source read of the second edge.
	listOverlapOffsets->push_back(overlapOffset);

	if(edge2->getOrientation() == 2 || edge2->getOrientation() == 3)
		listOrientations->push_back(1);
	else
		listOrientations->push_back(0);



	for(UINT64 i = 0; i < edge2->getListOfOrientations()->size(); i++)	// Get the list from the second edge.
	{
		listReads->push_back(edge2->getListOfReads()->at(i));
		listOverlapOffsets->push_back(edge2->getListOfOverlapOffsets()->at(i));
		listOrientations->push_back(edge2->getListOfOrientations()->at(i));
	}
	return true;
}




/**********************************************************************************************************************
	Orientation of the merged edge
**********************************************************************************************************************/
UINT8 OverlapGraph::mergedEdgeOrientationDisconnected(Edge *edge1, Edge *edge2)
{
	UINT8 or1 = edge1->getOrientation(), or2 = edge2->getOrientation(),returnValue;
	     if( (or1 == 1 || or1 == 0) && (or2 == 0 || or2 == 2) )	// <---------*  and *-----------<
		returnValue = 0;
	else if( (or1 == 1 || or1 == 0) && (or2 == 1 || or2 == 3) )	// <---------*  and *----------->
		returnValue = 1;
	else if( (or1 == 2 || or1 == 3) && (or2 == 0 || or2 == 2) )	// >---------*  and *-----------<
		returnValue = 2;
	else if( (or1 == 2 || or1 == 3) && (or2 == 1 || or2 == 3) )	// >---------*  and *----------->
		returnValue = 3;
	else
	{
		FILE_LOG(logDEBUG)<<(int)or1<<" "<<(int)or2;
		MYEXIT("Unable to merge.")
	}
	return returnValue;
}



/**********************************************************************************************************************
	Remove edges with similar endpoint in the overlap graph
**********************************************************************************************************************/

UINT64 OverlapGraph::removeSimilarEdges(void)
{
	CLOCKSTART;
	UINT64 counter = 0;
	vector <Edge *> listOfEdges1, listOfEdges2;
	vector <UINT64> listOfEditDistance;
	for(UINT64 i = 1; i < graph->size(); i++)	// For each node.
	{
		if(!graph->at(i)->empty())	// The node has some edges in the graph.
		{
			for(UINT64 j = 0; j < graph->at(i)->size(); j++)	// For all the edges.
			{
				Edge * e1 = graph->at(i)->at(j);
				UINT64 source1 = e1->getSourceRead()->getReadNumber(), destination1 = e1->getDestinationRead()->getReadNumber();
				if( source1 < destination1)
				{
					for(UINT64 k = j + 1; k < graph->at(i)->size(); k++)
					{
						Edge * e2 = graph->at(i)->at(k);
						UINT64 source2 = e2->getSourceRead()->getReadNumber(), destination2 = e2->getDestinationRead()->getReadNumber();
						if( source1 == source2 && destination1 == destination2)	// Check both edges starts and ends at same end points
						{
							if( abs((int)(e1->getOverlapOffset() - e2->getOverlapOffset())) <  abs((int)(e2->getOverlapOffset()/20)) ) // The lengths are more than 95% similar
							{
								string s1 = getStringInEdge(e1);	// Get the string on the first edge.
								string s2 = getStringInEdge(e2);	// Get the string on the second edge.
								UINT64 editDistance = calculateEditDistance(s1,s2);	// Calculate the edit distance between the strings.
								if( editDistance < min(e1->getOverlapOffset(), e2->getOverlapOffset())/20 )	// If the edit distance is less than 5% of the length of the shortest string.
								{
									UINT64 l;
									for(l= 0; l <  listOfEdges1.size(); l++)	// Already in the list.
									{
										if(listOfEdges2.at(l) ==  e2 || listOfEdges2.at(l) == e1) // check if the edges are already used.
											break;
									}
									if(l ==  listOfEdges1.size())	// Not in the list. Add it in the list.
									{
										getBaseByBaseCoverage(e1);	// JJ: Keep the edge with higher coverage depth
										getBaseByBaseCoverage(e2);
										if (e1->coverageDepth > e2->coverageDepth)
										{
											listOfEdges1.push_back(e1);	// We will keep this edge. Just keep the first edge?? how does one know this is the right one?? -JJ
											listOfEdges2.push_back(e2);	// This edge will be deleted and the flow will be moved to the first edge.
										}
										else
										{
											listOfEdges1.push_back(e2);	// We will keep this edge. Just keep the first edge?? how does one know this is the right one?? -JJ
											listOfEdges2.push_back(e1);	// This edge will be deleted and the flow will be moved to the first edge.
										}
										listOfEditDistance.push_back(editDistance);	// Also store the edit distance.
									}
								}
							}
						}
					}
				}
			}
		}
	}
	FILE_LOG(logINFO) << listOfEdges2.size()<< " edges to remove";
	for(UINT64 i = 0; i < listOfEdges2.size(); i++)
	{
		FILE_LOG(logDEBUG) <<  ++ counter << ": removing edge ("<<  listOfEdges2.at(i)->getSourceRead()->getReadNumber()<<"," <<  listOfEdges2.at(i)->getDestinationRead()->getReadNumber()<<")" ;
		FILE_LOG(logDEBUG) << "---> " << "Lengths : " <<  listOfEdges1.at(i)->getOverlapOffset() << " and " <<  listOfEdges2.at(i)->getOverlapOffset() ; 
		FILE_LOG(logDEBUG) << "---> " << "Flows: " << setw(3) << listOfEdges1.at(i)->flow << " and " << setw(3) << listOfEdges2.at(i)->flow ;
		FILE_LOG(logDEBUG) << "---> " << "Edit Distance: " << setw(5) << listOfEditDistance.at(i) ;
		FILE_LOG(logDEBUG) << "---> " << "Reads on the edges: " << listOfEdges1.at(i)->getListOfReads()->size() << " and " << listOfEdges2.at(i)->getListOfReads()->size() ;
		FILE_LOG(logDEBUG) << "---> " << "Coverage Depths: " << listOfEdges1.at(i)->coverageDepth << " and " << listOfEdges2.at(i)->coverageDepth << endl;
		listOfEdges1.at(i)->flow += listOfEdges2.at(i)->flow;	// Move the flow of the delete edge to this edge.
		listOfEdges1.at(i)->getReverseEdge()->flow += listOfEdges2.at(i)->getReverseEdge()->flow;	// Same for the reverse edge.
		removeEdge(listOfEdges2.at(i));	// Then remove the edge with similar string.
	}
	FILE_LOG(logINFO) << counter << " edges removed.";
	CLOCKSTOP;
	return counter;
}



/**********************************************************************************************************************
	Resolve ambiguous nodes according to coverage.
**********************************************************************************************************************/
UINT64 OverlapGraph::resolveNodes(void)
{
	CLOCKSTART;
	UINT64 counter = 0;
	vector <Edge *> listOfInEdges, listOfOutEdges;
	for(UINT64 i = 1; i < graph->size(); i++)
	{
		listOfInEdges.clear();
		listOfOutEdges.clear();
		if(graph->at(i)->size()==4)	// Nodes with 4 edges.
		{
			for(UINT64 j=0; j < graph->at(i)->size(); j++)
			{
				if(graph->at(i)->at(j)->getSourceRead() == graph->at(i)->at(j)->getDestinationRead())	// Will not consider the nodes that has loop
				{
					listOfInEdges.clear();
					listOfOutEdges.clear();
					break;
				}
				if(graph->at(i)->at(j)->getOrientation() == 0 || graph->at(i)->at(j)->getOrientation() == 1)	// In-edges.
				{
					listOfInEdges.push_back(graph->at(i)->at(j)->getReverseEdge());
				}
				else	// Out-edges
				{
					listOfOutEdges.push_back(graph->at(i)->at(j));
				}
			}
			if(listOfInEdges.size() == 2 && listOfOutEdges.size() == 2)	// If two in-edges and two out-edges.
			{

				Edge *inEdge1, *inEdge2, *outEdge1, *outEdge2;
	// Calculate the mean and SD of coverage depth of the unique reads in these edges.
				getBaseByBaseCoverage(listOfInEdges.at(0));
				getBaseByBaseCoverage(listOfInEdges.at(1));
				getBaseByBaseCoverage(listOfOutEdges.at(0));
				getBaseByBaseCoverage(listOfOutEdges.at(1));
				if(listOfInEdges.at(0)->coverageDepth > listOfInEdges.at(1)->coverageDepth)	// Sort the in-Edges.
				{
					inEdge1 = listOfInEdges.at(0);
					inEdge2 = listOfInEdges.at(1);
				}
				else
				{
					inEdge1 = listOfInEdges.at(1);
					inEdge2 = listOfInEdges.at(0);
				}

				if(listOfOutEdges.at(0)->coverageDepth > listOfOutEdges.at(1)->coverageDepth)	// Sort the out-edges.
				{
					outEdge1 = listOfOutEdges.at(0);
					outEdge2 = listOfOutEdges.at(1);
				}
				else
				{
					outEdge1 = listOfOutEdges.at(1);
					outEdge2 = listOfOutEdges.at(0);
				}

				UINT64 flag1 = 0, flag2 = 0;
				if(isOverlappintInterval(inEdge1->coverageDepth, inEdge1->SD, outEdge1->coverageDepth, outEdge1->SD) && !isOverlappintInterval(inEdge1->coverageDepth, inEdge1->SD, outEdge2->coverageDepth, outEdge2->SD) && !isOverlappintInterval(inEdge2->coverageDepth, inEdge2->SD, outEdge1->coverageDepth, outEdge1->SD))
				{
					flag1 = 1;
				}
				if(isOverlappintInterval(inEdge2->coverageDepth, inEdge2->SD, outEdge2->coverageDepth, outEdge2->SD) && !isOverlappintInterval(inEdge2->coverageDepth, inEdge2->SD, outEdge1->coverageDepth, outEdge1->SD) && !isOverlappintInterval(inEdge1->coverageDepth, inEdge1->SD, outEdge2->coverageDepth, outEdge2->SD))
				{
					flag2 = 1;
				}
				if(flag1 == 1 )
				{
					counter++;
					FILE_LOG(logDEBUG) << setw(10) << counter << " Merging edges (" << setw(10) << inEdge1->getSourceRead()->getReadNumber() << "," << setw(10) <<inEdge1->getDestinationRead()->getReadNumber() << ") Length: " << setw(6) << inEdge1->getOverlapOffset() << " Flow: " << setw(3) << inEdge1->flow << " Coverage: " << setw(4) << inEdge1->coverageDepth << " SD: " << setw(3) << inEdge1->SD << " and (" << setw(10) << outEdge1->getSourceRead()->getReadNumber() << "," << setw(10) << outEdge1->getDestinationRead()->getReadNumber() << ") Length: " << setw(6) << outEdge1->getOverlapOffset() << " Flow: " << setw(3) << outEdge1->flow << " Coverage: " << setw(4) << outEdge1->coverageDepth << " SD: " << setw(3) << outEdge1->SD;
					mergeEdges(inEdge1,outEdge1);	// Merge the edges.
				}
				if(flag2 == 1)
				{
					counter++;
					FILE_LOG(logDEBUG) << setw(10) << counter << " Merging edges (" << setw(10) << inEdge2->getSourceRead()->getReadNumber() << "," << setw(10) <<inEdge2->getDestinationRead()->getReadNumber() << ") Length: " << setw(6) << inEdge2->getOverlapOffset() << " Flow: " << setw(3) << inEdge2->flow << " Coverage: " << setw(4) << inEdge2->coverageDepth << " SD: " << setw(3) << inEdge2->SD << " and (" << setw(10) << outEdge2->getSourceRead()->getReadNumber() << "," << setw(10) << outEdge2->getDestinationRead()->getReadNumber() << ") Length: " << setw(6) << outEdge2->getOverlapOffset() << " Flow: " << setw(3) << outEdge2->flow << " Coverage: " << setw(4) << outEdge2->coverageDepth << " SD: " << setw(3) << outEdge2->SD;
					mergeEdges(inEdge2,outEdge2);	// Merge the edges.
				}
			}
		}
	}
	FILE_LOG(logINFO)<< counter << " edges merged.";
	CLOCKSTOP;
	return counter;
}




struct stackElement
{
	Edge * edge;
	UINT64 distance;
};

struct overlappingReads
{
	Edge *edge;
	UINT64 readID;
	UINT64 distance;
};



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
	This function remove loops that can be traversed in only one way.
	a--->b--->b--->c
**********************************************************************************************************************/
UINT64 OverlapGraph::reduceLoops(void)
{
	CLOCKSTART;
	UINT64 counter = 0;
	Edge *ab,*bb,*bc;
	for(UINT64 i = 1; i < graph->size(); i++)
	{
		if(graph->at(i)->size() == 4) // only four edges. The loop is counted twice.
		{
			UINT64 loopCount = 0, incomingEdgeCount = 0, outgoingEdgeCount = 0;
			for(UINT64 j = 0; j< graph->at(i)->size(); j++)
			{
				if(graph->at(i)->at(j)->getDestinationRead()->getReadNumber() == i) // This is a loop
				{
					loopCount++;
					bb = graph->at(i)->at(j);
				}
				else if(graph->at(i)->at(j)->getOrientation() == 0 || graph->at(i)->at(j)->getOrientation() == 1) // incoming edge
				{
					incomingEdgeCount++;
					ab = graph->at(i)->at(j)->getReverseEdge();
				}
				else if(graph->at(i)->at(j)->getOrientation() == 2 || graph->at(i)->at(j)->getOrientation() == 3) // outgoing edge
				{
					outgoingEdgeCount++;
					bc = graph->at(i)->at(j);
				}
			}
			if(loopCount==2 && incomingEdgeCount == 1 && outgoingEdgeCount == 1)  // two in the loop and one incoming and one outgoing
			{
				FILE_LOG(logDEBUG)<<"Loop found at node: " << i <<  " loop edge length: " << bb->getOverlapOffset() << " flow: " << bb->flow << " Other edge lengths: " << ab->getOverlapOffset() << " and " << bc->getOverlapOffset();
	//cout << getStringInEdge(bb) << endl;
				if(bb->getOrientation() == 0)
				{
					counter++;
					mergeEdges(ab,bb->getReverseEdge());
				}
				else if(bb->getOrientation() == 3)
				{
					counter++;
					mergeEdges(ab,bb);
				}
				else	// Arrow in the edge is >---< or <---->. In this case it is not possible to decide which way to go.
				{
					FILE_LOG(logDEBUG) << "Unable to reduce loop because of the edge type.";
				}
			}
		}
	}
	FILE_LOG(logINFO) <<" Loops removed: " << counter; // Number of loop we were able to reduce
	CLOCKSTOP;
	return counter;
}

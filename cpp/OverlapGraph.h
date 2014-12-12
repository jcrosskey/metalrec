/*
 * OverlapGraph.h
 *
 * Created on: Mon Dec  1 14:25:04 EST 2014
 * Author: JJ Chai
 */


#ifndef OVERLAPGRAPH_H_
#define OVERLAPGRAPH_H_
#include "Common.h"
#include "Dataset.h"
#include "Edge.h"
#include <seqan/align.h>

/**********************************************************************************************************************
 * Class to store the overlap graph, hash table is no longer used here. 
 * Overlap graph is a directed graph, built from dataSet including all reads' information and pairwise alignment
 * between the reads.
 **********************************************************************************************************************/

enum nodeType {
	UNEXPLORED = 0,	// Current node u is not explored yet. Meaning that there is no edge (u,v) in the graph.
	EXPLORED = 1,	//  Current node u is explored. Meaning that all edges (u,v) are inserted in the dataset.
	EXPLORED_AND_TRANSITIVE_EDGES_MARKED = 2	// Meaning that all transitive edges (u,v) of current node u is marked and its neighbors transitive edges are also marked. Now it is safe to remove any transitive edge (u,v) from node u.
};

/* mark types for the nodes */
enum markType{
	VACANT = 0,
	INPLAY = 1,
	ELIMINATED = 2
};

class OverlapGraph
{
	private:
		UINT64 minimumOverlapLength;	// Length of the minimum overlap to consider two reads as overlapping
		UINT32 maxError;	// maximum number of substitutions allowed in overlap between Illumina reads
		float maxErrorRate;	// maximum error rate
		UINT64 numberOfNodes;	// Number of nodes in the overlap graph.
		UINT64 numberOfEdges;	// Number of edges in the overlap graph.
		INT32 rubberPos;	// Number of base pairs allowed for the errors of the starting and ending coordinates.

		Dataset * dataSet; 	// Pointer to the dataset containing all the reads.
		vector< vector<Edge *> * > *graph;	// Adjacency list of the graph.

		bool mergeList(Edge *edge1, Edge *edge2, vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets);
		bool DoReadsOverlap(Read * read1, Read * read2, INT16 & OverlapOffset);

	public:
		bool flowComputed;	// Flag to check wheather the flow is computed or not.
		OverlapGraph(void);	// Default constructor.
		OverlapGraph(Dataset *data_Set, const UINT64 & minOverlap, const UINT32 & max_Error, const float & max_ErrorRate, const INT32 & rubber_pos);	// Another constructor, from dataSet including all the reads, with specified values for other parameters. Default values are also included
		~OverlapGraph();	// Destructor.

		bool buildOverlapGraphFromDataSet(Dataset *data_Set);	// Build the overlap graph using dataSet.
		bool markTransitiveEdges(UINT64 readNumber, vector<markType> * markedNodes); // Mark transitive edges of a read.
		bool insertEdge(Edge * edge); 	// Insert an edge in the overlap graph.
		bool insertEdge(Read *read1, Read *read2, UINT16 overlapOffset); // Insert an edge in the overlap graph.
		UINT64 contractCompositePaths(void); 	// Contract composite paths in the overlap graph.
		bool removeEdge(Edge *edge); 	// Remove an edge from the overlap graph.
		bool printGraph(string graphFileName, string contigFileName);	// Store the overlap graph for visual display and also store the contigs/scaffods in a file.
		UINT64 getAllOverlaps(void);	// Find all the overlaps between reads in the graph, and mark the contained reads in the process.
		bool insertAllEdgesOfRead(UINT64 readNumber, vector<nodeType> * exploredReads);	// Insert into the overlap graph all edges of a read.
		bool removeTransitiveEdges(UINT64 readNumber);	// Remove all transitive edges from the overlap graph incident to a given read.
		bool removeEdgesOfRead(Read * read);	// Remove all edges adjacent to a read
		bool mergeEdges(Edge *edge1, Edge *edge2);	// Merge two edges in the  overlap graph.
		UINT64 removeAllSimpleEdgesWithoutFlow();	// Remove simple edges without flow
		bool updateReadLocations(Edge *edge);	// Update the location of all the reads in the current edge. This function is called when a new edge is inserted.
		bool removeReadLocations(Edge *edge);	// Remove the location of all the reads from the current edge. This function is called when an edge is removed.
		UINT64 getNumberOfEdges(void){return numberOfEdges;}	// Get the number of edges in the overlap graph.
		UINT64 getNumberOfNodes(void){return numberOfNodes;}	// Get the number of nodes in the overlap graph.
		bool setDataset(Dataset *data_Set){dataSet=data_Set; return true;}	// Set the dataset pointer.
		bool setRubberPos(const INT32 rubber_pos){rubberPos = rubber_pos; return true;}	// Set the rubber base pairs.
		Edge *findEdge(UINT64 source, UINT64 destination);	// Find an edge from source to destination in the overlap graph.
		bool isEdgePresent(UINT64 source, UINT64 destination);	// Check if an edge is present in the overlap graph between source and destination.
		bool calculateBoundAndCost(Edge *edge, INT64* FLOWLB, INT64* FLOWUB, INT64* COST); // Calculate bounds and costs of flow for minimum cost flow in the overlap graph.
		bool calculateFlow(string inputFileName, string outputFileName);									// Calculate the minimum cost flow of the overlap graph.
		seqan::DnaString getStringInEdge(Edge *edge);	// Get the string in an edge by overlapping the ordered reads in the edge.
		bool simplifyGraph(void);	// Some simple simplification.
		void getBaseByBaseCoverage(Edge *edge);	// Get the coverage Mean and SD of an edge. Only considering the unique reads.
		void sortEdges();	// Sort edges of each read based on ID of the destination read.
		UINT64 calculateEditDistance(const std::string &s1, const std::string &s2);	// Find the edit distance between two strings.
};



#endif /* OVERLAPGRAPH_H_ */

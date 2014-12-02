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
		UINT64 minimumOverlapLength;	// Length of the shortest read in the dataset.
		UINT32 maxError;	// maximum number of substitutions allowed in overlap between Illumina reads
		float maxErrorRate;	// maximum error rate
		Dataset * dataSet; 	// Pointer to the dataset containing all the reads.
		vector< vector<Edge *> * > *graph;	// Adjacency list of the graph.
		UINT64 numberOfNodes;	// Number of nodes in the overlap graph.
		UINT64 numberOfEdges;	// Number of edges in the overlap graph.
		bool mergeList(Edge *edge1, Edge *edge2, vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets);
		UINT64 exploreGraph(Edge* firstEdge, Edge * lastEdge, UINT64 distanceOnFirstEdge, UINT64 distanceOnLastEdge, UINT64 datasetNumber, UINT64 level, vector <Edge *> &firstPath, vector <UINT64> &flags);

	public:
		OverlapGraph(void);	// Default constructor.
		OverlapGraph(DataSet *data_Set, UINT64 minOverlap, UINT32 max_Error, float max_ErrorRate);	// Another constructor, from dataSet including all the reads
		~OverlapGraph();	// Destructor.
		bool markTransitiveEdges(UINT64 readNumber, vector<markType> * markedNodes); // Mark transitive edges of a read.
		bool buildOverlapGraphFromDataSet(DataSet *data_Set);	// Build the overlap graph using dataSet.
		bool insertEdge(Edge * edge); 	// Insert an edge in the overlap graph.
		bool insertEdge(Read *read1, Read *read2, UINT16 overlapOffset); // Insert an edge in the overlap graph.
		UINT64 contractCompositePaths(void); 	// Contract composite paths in the overlap graph.
		UINT64 removeDeadEndNodes(void); 					        // Remove dead-ends from the overlap graph.
		bool removeEdge(Edge *edge); 	// Remove an edge from the overlap graph.
		//bool checkOverlap(Read *read1, Read *read2, UINT64 orient, UINT64 start); // Check overlap between two reads after a match is found using the hash table.
		//bool checkOverlapForContainedRead(Read *read1, Read *read2, UINT64 orient, UINT64 start);
		bool printGraph(string graphFileName, string contigFileName);	// Store the overlap graph for visual display and also store the contigs/scaffods in a file.
		bool insertAllEdgesOfRead(UINT64 readNumber, vector<nodeType> * exploredReads);	// Insert into the overlap graph all edges of a read.
		bool removeTransitiveEdges(UINT64 readNumber);	// Remove all transitive edges from the overlap graph incident to a given read.
		bool mergeEdges(Edge *edge1, Edge *edge2);	// Merge two edges in the  overlap graph.
		bool updateReadLocations(Edge *edge);	// Update the location of all the reads in the current edge. This function is called when a new edge is inserted.
		bool removeReadLocations(Edge *edge);	// Remove the location of all the reads from the current edge. This function is called when an edge is removed.
		UINT64 getNumberOfEdges(void){return numberOfEdges;}	// Get the number of edges in the overlap graph.
		UINT64 getNumberOfNodes(void){return numberOfNodes;}	// Get the number of nodes in the overlap graph.
		bool saveGraphToFile(string fileName);	// Save the unitig graph to a file. This file can be used to reproduced the graph. Useful to avoid recomputing the graph.
		bool readGraphFromFile(string fileName);	// Read the graph from a file stored before.
		bool setDataset(Dataset *data_Set){dataSet=data_Set; return true;}	// Set the dataset pointer.
		Edge *findEdge(UINT64 source, UINT64 destination);	// Find an edge from source to destination in the overlap graph.
		bool isEdgePresent(UINT64 source, UINT64 destination);	// Check if an edge is present in the overlap graph between source and destination.
		string getStringInEdge(Edge *edge);	// Get the string in an edge by overlapping the ordered reads in the edge.
		//UINT64 reduceTrees(void);	// Remove trees in the overlap graph.
		bool simplifyGraph(void);	// Some simple simplification.
		//UINT64 scaffolder(void);	// Construct scaffolds using matepair information.
		UINT64 removeSimilarEdges(void);	// remove multi-edges with similar strings
		//UINT64 resolveNodes(void);	// Merge edges based on coverage depth
		void markContainedReads(void);	// Find superReads for each read and mark them as contained read.
		void getBaseByBaseCoverage(Edge *edge);	// Get the coverage Mean and SD of an edge. Only considering the unique reads.
		void sortEdges();	// Sort edges of each read based on ID of the destination read.
		UINT64 findOverlap(string string1, string string2);	// Find overlap length between two strings.
		UINT64 calculateEditDistance(const std::string &s1, const std::string &s2);	// Find the edit distance between two strings.
		UINT64 reduceLoops(void);	// loops that can be traversed only one way
};



#endif /* OVERLAPGRAPH_H_ */

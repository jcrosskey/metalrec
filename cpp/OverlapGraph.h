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
#include "HashTable.h"
#include "Edge.h"

/**********************************************************************************************************************
 * Class to store the overlap graph, hash table is used to find overlaps, no more pairwise alignment is involved
 * Overlap graph is a directed graph, since the strand is already decided by the reads' alignments to the PacBio read.
 * Mismatches of substitution type is allowed for now
 **********************************************************************************************************************/

enum nodeType {
	UNEXPLORED = 0,	// Current node u is not explored yet. Meaning that there is no edge (u,v) in the graph.
	EXPLORED = 1,	//  Current node u is explored. Meaning that all edges (u,v) are inserted in the dataset.
	EXPLORED_AND_TRANSITIVE_EDGES_MARKED = 2	// Meaning that all transitive edges (u,v) of current node u is marked and its neighbors transitive edges are also marked. Now it is safe to remove any transitive edge (u,v) from node u.
};

enum nodeColor{                                 /* Used for depth first search */
	WHITE = 0, /* unexplored at all */
	GRAY = 1, /* First time the node is discovered */
	BLACK = 2 /* Node is finished, i.e. all the children have been explored */
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
		UINT32 rubberPos;	// Number of base pairs allowed for the errors of the starting and ending coordinates.
		UINT64 numberOfNodes;	// Number of nodes in the overlap graph.
		UINT64 numberOfEdges;	// Number of edges in the overlap graph.
		UINT64 numberOfEdgesEverMade;	// Number of edges created, simple or composite, during the layout process

		Dataset * dataSet; 	// Pointer to the dataset containing all the reads.
		HashTable *hashTable;           /* pointer to the hash table */
		vector< vector<Edge *> * > *graph;	// Adjacency list of the graph.
		vector<vector<UINT64> > reverseGraph; /* graph with edges in reverse direction, used for easy track of edges from nodes. */

		bool mergeList(Edge *edge1, Edge *edge2, vector<UINT64> *listReads, vector<UINT16> *listOverlapOffsets, vector<UINT64> *listOfSubstitutionPoses); /* Merge edges into composite edge */

	public:
		OverlapGraph(void);	// Default constructor.
		OverlapGraph(HashTable *ht, const UINT64 & minOverlap, const UINT32 & max_Error, const float & max_ErrorRate, const UINT32 & rubber_pos);	// Another constructor, from dataSet including all the reads, with specified values for other parameters. Default values are also included
		~OverlapGraph();	// Destructor.
		OverlapGraph(const OverlapGraph & O);	// Copy constructor
		OverlapGraph & operator= (const OverlapGraph & O);	// Copy assignment

		bool buildOverlapGraphFromHashTable(HashTable *ht);	// Build the overlap graph using dataSet.
		void buildReverseGraph(void); /* build the reverse graph */
		void markContainedReads(void);	// Find superReads for each read and mark them as contained read.
		bool checkOverlapForContainedRead(Read *read1, Read *read2, UINT64 orient, UINT64 start);
		bool checkOverlap(Read *read1, Read *read2, UINT64 orient, UINT64 start, UINT16 & numSub, vector<UINT64> * listSubs);
		bool checkOverlapWithSub(const string & str1, const string & str2, UINT16 & numSub, vector<UINT64> * listSubs, UINT64 orient); /* check if two strings overlap with substitutions within specified range */
		bool markTransitiveEdges(UINT64 readNumber, vector<markType> * markedNodes); // Mark transitive edges of a read.
		bool insertEdge(Edge * edge); 	// Insert an edge in the overlap graph.
		bool insertEdge(Read *read1, Read *read2, UINT16 overlapOffset, UINT16 numSub, vector<UINT64> *listSubs); // Insert an edge in the overlap graph.
		UINT64 contractCompositePaths(void); 	// Contract composite paths in the overlap graph.
		bool removeEdge(Edge *edge); 	// Remove an edge from the overlap graph.
		bool getEdges(vector<Edge *> & contigEdges); /* save all the edges in the graph in a vector of pointers to these edges */
		bool printGraph(string graphFileName, const vector<Edge *> & contigEdges);	// Store the overlap graph for visual display and also store the contigs/scaffods in a file.
		bool printContigs(string outputFastaName, vector<Edge *> & contigEdges, bool longestOnly);
		bool insertAllEdgesOfRead(UINT64 readNumber, vector<nodeType> * exploredReads);	// Insert into the overlap graph all edges of a read.
		bool removeTransitiveEdges(UINT64 readNumber);	// Remove all transitive edges from the overlap graph incident to a given read.
		bool removeEdgesBetweenReadNumbers(UINT64 readNumber1, UINT64 readNumber2);
		Edge * mergeEdges(Edge *edge1, Edge *edge2);	// Merge two edges in the  overlap graph.
		bool updateReadLocations(Edge *edge);	// Update the location of all the reads in the current edge. This function is called when a new edge is inserted.
		bool removeReadLocations(Edge *edge);	// Remove the location of all the reads from the current edge. This function is called when an edge is removed.
		UINT64 getNumberOfEdges(void){return numberOfEdges;}	// Get the number of edges in the overlap graph.
		UINT64 getNumberOfNodes(void){return numberOfNodes;}	// Get the number of nodes in the overlap graph.
		bool setDataset(Dataset *data_Set){dataSet=data_Set; return true;}	// Set the dataset pointer.
		bool setRubberPos(const UINT32 rubber_pos){rubberPos = rubber_pos; return true;}	// Set the rubber base pairs.
		vector<Edge *> findEdge(UINT64 source, UINT64 destination);	// Find an edge from source to destination in the overlap graph.
		Edge * findEdge(UINT64 source, UINT64 destination, UINT64 edgeID);
		bool isEdgePresent(UINT64 source, UINT64 destination);	// Check if an edge is present in the overlap graph between source and destination.
		string getStringInEdge(Edge *edge, bool includeLast=true, bool leftClip=false, bool rightClip=false);	// Get the string in an edge by overlapping the ordered reads in the edge.
		bool simplifyGraph(void);	// Some simple simplification.
		bool DFS(vector< UINT64 > * topoSortedNodes); /* Depth first search of the graph, can also be used to determine when the graph is cyclic */
		bool DFS_visit(UINT64 i, vector<int> *searchStatus, vector<UINT64> *predecessors, UINT64 & sTime, vector<int> * discoverTime, vector<int> * finishTime, vector<Edge *> * backEdges, vector<UINT64> * topoSortedNodes); /* visit all the neighbors of a node, and update the DFS vectors */
		bool FindLongestPath(vector<UINT64> * topoSortedNodes, string & finalString, UINT64 & beginCoord, UINT64 & endCoord, double & maxWeight, UINT64 & totalOverlapLength, UINT64 & totalReads, UINT64 & totalSubs, UINT64 & totalIns, UINT64 & totalDel,UINT64 & totalClip, UINT64 & totalReadLength);         /* Find longest path in the graph */
};

#endif /* OVERLAPGRAPH_H_ */

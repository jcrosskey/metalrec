/*
 * main_ec.h
 *
 * Created on: Tue Nov 11 10:58:50 EST 2014
 * Author: JJ Chai
 */

/*** utility functions ***/
#include "Common.h"
/*** Read class ***/
#include "Read.h"
/*** Dataset class ***/
#include "Dataset.h"
/*** Edge class ***/
#include "Edge.h"
/*** HashTable class ***/
#include "HashTable.h"
/*** OverlapGraph class ***/
#include "OverlapGraph.h"
/* check file existence and permissions */
#include <unistd.h>
int loglevel; // logging level in integer format, for different levels of verbosity

/************************
 * main function
 * **********************/
int main_ec(FILE * samStream, const  string & outputFastaName) /* From sam file in stream, do error correction, and write output to outputFastaName */
{
	CLOCKSTART;

	/** Read sam file and store all the reads **/
	Dataset * dataSet = new Dataset(samStream, 40, 0.25, 0.05);// read from the stdin stream
	FILE_LOG(logINFO) << "Length of the PacBio read is " << dataSet->getPacBioReadLength(); /* print PacBio read length */
	if (dataSet->getNumberOfReads() <= 1)
		FILE_LOG(logERROR) << "Data set has no more than 1 read in it, quitting...";

	else /* If the dataset has some reads in it, build hash table and overlap graph next */
	{
		FILE_LOG(logINFO) << "number of unique reads in dataset is " << dataSet->getNumberOfUniqueReads();
		HashTable *ht = new HashTable();
		ht->insertDataset(dataSet, 10);
		OverlapGraph *graph = new OverlapGraph(ht, 40, 0, 0, 10);
		delete ht;                      /* delete hash table after overlap graph is built */
		if (graph->getNumberOfEdges() == 0)
			FILE_LOG(logERROR) << "Data set  has no edge in it, quitting...";

		else /* If there is at least 1 edge in the data set, try to calculate flow and output contigs */
		{
			graph->calculateFlow();
			FILE_LOG(logINFO) << "nodes: " << graph->getNumberOfNodes() << " edges: " << graph->getNumberOfEdges() << endl;
			graph->removeAllSimpleEdgesWithoutFlow();
			graph->simplifyGraph();
			vector<Edge *> contigEdges;
			graph->getEdges(contigEdges);
			graph->printContigs(outputFastaName, contigEdges,true);

		}
		delete graph;
	}

	delete dataSet;

	CLOCKSTOP;
	return 0;
}

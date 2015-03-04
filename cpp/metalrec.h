#include <stdexcept> // for standard exceptions out_or_range
#include "Utils.h"
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

/********************************************************
 * metalrec.h

 * metalrec that correct 1 PacBio read from a number of bam files.
 * Input file is BAM file(s) with Illumina reads aligned to PacBio reads.
 * samtools is required
 *
 * Created by JJ Chai  Thu Sep 18 09:35:28 EDT 2014. Last modified Feb 09, 2015. 
 * Copyright (c) 2014 JJ Chai (ORNL). Allrights reserved.

 ********************************************************/

// error correction job
void metalrec(const vector<string> & bamFiles, const string & PacBioName, const string & allFileName,
		const string & samtools_path, const string & outDir,
		const UINT64 & minimumOverlapLength, const UINT64 & hashStringLength,
		const UINT32 & maxError, const UINT32 &rubberPos,
		const float & indelRate, const float & subRate, const float & maxErrorRate)
{
	CLOCKSTART;
	if(PacBioName.length() > 0)
	{
		/* output file */
		string outFile = outDir + "/" + allFileName + ".fasta"; // output corrected fasta file
		if (Utils::isFileExist(outFile))
		{
			FILE_LOG(logINFO) << allFileName << " is already done. skip";
		}
		else{
			Utils::mkdirIfNonExist(outDir);

			/* Initiate Dataset object, and set the indel and substitution rates allowed in the alignment */
			Dataset * dataSet = new Dataset();
			dataSet->setIndelRate(indelRate);
			dataSet->setSubRate(subRate);

			/* Read all bam file to collect Illumina reads aligned to this PacBio read */
			for (size_t i = 0; i < bamFiles.size(); i++)
			{
				string getSamCmd = samtools_path + " view " + bamFiles[i] + " " + PacBioName;
				FILE_LOG(logINFO) << "Reading bam file " << bamFiles[i];
				/* popen returns NULL if fork or pipe calls fail, or if it cannot allocate memory */
				FILE * sam_pipe = popen(getSamCmd.c_str(), "r"); /* stream to capture the output of "samtools view bamfile refname" */
				if(!sam_pipe){
					cerr << "   *** Failed command " << getSamCmd << endl;
					perror("Error encountered in popen()"); /* If fork or pipe fails, errno is set */
					return;
				}

				else{
					/** Read sam file and store all the reads **/
					dataSet->AddDataset(sam_pipe);
					int exit_status = pclose(sam_pipe); /* pclose() waits for the associated process to terminate and returns the exit status of the command as returned by wait4; 
									       returns -1 if wait4 returns an error, or some other error is detected. */
					if (exit_status == -1)
					{
						perror("Error encountered in pclose()");
						return;
					}
				}
			}
			dataSet->finalize();
			dataSet->setPacBioReadName(PacBioName);

			if (dataSet->getNumberOfReads() <= 1)
				FILE_LOG(logWARNING) << "Data set has no more than 1 read in it, quitting...";

			else /* If the dataset has some reads in it, build hash table and overlap graph next */
				/* Work flow: 
				 * 1. build hash table 
				 * 2. build overlap graph based on the hash table
				 * 3. contract composite edges
				 * 4. flow analysis
				 * 5. remove edges without flow
				 * 6. simplify graph (contract composite edges, and pop bubbles)
				 * */
			{
				FILE_LOG(logINFO) << "number of unique reads in dataset is " << dataSet->getNumberOfUniqueReads();
				HashTable *ht = new HashTable();
				ht->insertDataset(dataSet, hashStringLength);
				OverlapGraph *graph = new OverlapGraph(ht, minimumOverlapLength, maxError, maxErrorRate, rubberPos);
				if(loglevel > 4)
				{
					vector<Edge *> contigEdges;
					graph->getEdges(contigEdges);
					graph->printGraph(outDir + "/" + allFileName + ".init.gdl", contigEdges);
				}
				delete ht;                      /* delete hash table after overlap graph is built */
				if (graph->getNumberOfEdges() == 0)
					FILE_LOG(logWARNING) << "Data set  has no edge in it, quitting...";

				else /* If there is at least 1 edge in the data set, try to calculate flow and output contigs */
				{
					graph->calculateFlow();
					FILE_LOG(logINFO) << "nodes: " << graph->getNumberOfNodes() << " edges: " << graph->getNumberOfEdges();
					graph->removeAllSimpleEdgesWithoutFlow();
					if(loglevel > 4)
					{
						vector<Edge *> contigEdges;
						graph->getEdges(contigEdges);
						graph->printGraph(outDir + "/" + allFileName + ".afterFlow.gdl", contigEdges);
					}
					graph->simplifyGraph();
					vector<Edge *> contigEdges;
					graph->getEdges(contigEdges);
					if(loglevel > 4)
					{
						graph->printGraph(outDir + "/" + allFileName + ".final.gdl", contigEdges);
						graph->printContigs(outDir + "/" + allFileName + ".final.fasta",contigEdges,false);
					}
					graph->printContigs(outFile, contigEdges,true);
					/* Use the path with most reads in it instead of the longest contigs */
					//				vector< Edge *>  allPaths;
					//				graph->findPaths(allPaths);
					//				if(loglevel > 4)
					//				{
					//					graph->printPaths(outDir + "/" + allFileName + ".paths.fasta",allPaths,false);
					//				}
					//				graph->printPaths(outFile, allPaths, true); /* Use the path with longest span as final output */

				}
				delete graph;
			}
			delete dataSet;
		}
	}
	else{
		cerr << "PacBio read name is empty\n";
	}

	CLOCKSTOP;
}

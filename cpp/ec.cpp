/*
 * =====================================================================================
 *
 *       Filename:  ec.cpp
 *
 *    Description:  ec corrects 1 PacBio read from a number of bam files.  
 *    Input file is BAM file(s) with Illumina reads aligned to PacBio reads.  
 *    samtools is required
 *
 *        Version:  1.0
 *        Created:  03/20/2015 13:04:47
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  JJ Chai (jjchai), jjchai01@gmail.com
 *   Organization:  ORNL
 *
 * =====================================================================================
 */

#include "ec.h"
/* Version where the output stream for error corrected contigs is not given */
void ec(const vector<string> & bamFiles, const string & PacBioName,  
		const UINT16 & PacBioLength, const string & allFileName,
		const string & samtools_path, const string & outDir,
		const UINT64 & minimumOverlapLength, const UINT64 & hashStringLength,
		const UINT32 & maxError, const UINT32 &rubberPos,
		const float & indelRate, const float & insRate, const float & delRate, const float & subRate, const float & maxErrorRate, const UINT16 minPacBioLength)
{
	string outFile = outDir + "/" + allFileName + ".fasta"; // output corrected fasta file
	ofstream outFastaStream;
	outFastaStream.open(outFile.c_str(), ofstream::out); /* append to the output file */
	if(!outFastaStream.is_open())
	{
		MYEXIT("Unable to open file: " + outFile);
	}
	else
	{
		ec_stream(bamFiles, PacBioName, PacBioLength, allFileName, outFastaStream, samtools_path, outDir, minimumOverlapLength, hashStringLength, maxError, rubberPos,indelRate, insRate, delRate, subRate, maxErrorRate, minPacBioLength);
	}
	outFastaStream.close();
}

/* Version where the output stream for error corrected contigs is given */
void ec_stream(const vector<string> & bamFiles, const string & PacBioName,  
		const UINT16 & PacBioLength, const string & allFileName, ofstream & outFastaStream,
		const string & samtools_path, const string & outDir,
		const UINT64 & minimumOverlapLength, const UINT64 & hashStringLength,
		const UINT32 & maxError, const UINT32 &rubberPos,
		const float & indelRate, const float & insRate, const float & delRate, const float & subRate, const float & maxErrorRate, const UINT16 minPacBioLength)
{
	CLOCKSTART;
	if(PacBioName.length() > 0)
	{
		if (PacBioLength < minPacBioLength){
			cerr << "PacBio read length " << PacBioLength << " is smaller than " << minPacBioLength << endl;
		}
		else{
			/* output file */
			string outFile = outDir + "/" + allFileName + ".fasta"; // output corrected fasta file
			Utils::mkdirIfNonExist(outDir);

			/* Initiate Dataset object, and set the indel and substitution rates allowed in the alignment */
			Dataset * dataSet = new Dataset();
			dataSet->setIndelRate(indelRate);
			dataSet->setInsRate(insRate);
			dataSet->setDelRate(delRate);
			dataSet->setSubRate(subRate);
			dataSet->setLRLength(PacBioLength);

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
			if(loglevel > 3)
			{
				dataSet->saveReads(outDir + "/" + allFileName + ".reads");
			}
			FILE_LOG(logINFO) << "Average coverage depth is: " << dataSet->getAvgCoverage();

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
				if(loglevel > 3)
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
					UINT64 stringLen, beginCoord, endCoord;
					int iter = 0;
					do{
						vector<UINT64> * topoSortedNodes = new vector<UINT64>;
						string finalString;
						graph->DFS(topoSortedNodes);
						graph->FindLongestPath(topoSortedNodes, finalString, beginCoord, endCoord);
						delete topoSortedNodes;
						stringLen = finalString.length();
						iter++;

						if(loglevel > 3)
						{
							vector<Edge *> contigEdges;
							graph->getEdges(contigEdges);
							graph->printGraph(outDir + "/" + allFileName + to_string(iter) + ".gdl", contigEdges);
						}
						FILE_LOG(logINFO) << "After longest path number " << iter << " is printed, number of edges left is: " << graph->getNumberOfEdges();
						// Print the sequence to output file
						outFastaStream << ">" << dataSet->getPacBioReadName() << "#" << iter << " Length_" << stringLen \
							<< " from_" << beginCoord << " to_" << endCoord \
							<< " numReads_" << dataSet->getNumberOfUniqueReads() \
							<< " covDepth_" << dataSet->getAvgCoverage() << " origLen_" << PacBioLength  \
							<< endl;

						UINT32 start=0;
						do
						{
							outFastaStream << finalString.substr(start, 100) << endl;  // save 100 BP in each line.
							start+=100;
						} while (start < stringLen);
					}while( stringLen > 1000 && graph->getNumberOfEdges() > 0);
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

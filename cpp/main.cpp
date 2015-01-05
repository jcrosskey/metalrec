/*
 * main.cpp
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
int loglevel; // logging level in integer format, for different levels of verbosity

void usage();
void parseArguments(int argc, char **argv, string & inputSamFile, string & allFileName, UINT64 & minimumOverlapLength, UINT64 & hashStringLength, UINT32 & maxError, UINT32 & rubPos, float & maxErrorRate, bool & useCoverageDepth, string & outputDir, string & outputFastaName, float & indelRate, float & subRate);

/************************
 * Help usage
 * **********************/
void usage()
{
	std::cout << std::endl
		<< "  Usage:" << std::endl
		<< "    metalrec [OPTION]...<PARAM>..." << std::endl
		<< std::endl
		<< "  <PARAM>" << std::endl
		<< "    -s\tinputSamFile" << std::endl  // Input sam file, generated by BLASR, with soft clipping
		<< "    -l\tminimum overlap length" << std::endl  // minimum overlap length
		<< std::endl
		<< "  [OPTION]" << std::endl
		<< "    -h/--help" << std::endl
		<< "    -f\tOutput file name prefix (default: metalrec)" << std::endl    // all output file with have this name with different extensions.
		<< "    -od\tOutput directory (default: ./ current working directory)" << std::endl    // output directory
		<< "    -o\tFasta files full path name with corrected sequence (default: outputDir + inputSam.fasta)" << std::endl    // output directory
		<< "    -k\thash string length (kmer length, default: 31)" << std::endl  // hash string length
		<< "    -r\trubber length (default: 10)" << std::endl  // hash string length
		<< "    -e\tMaximum (substitution) errors allowed in the overlap (default: 0)" << std::endl    // Maximum (substitution) errors allowed in the overlap between Illumina reads
		<< "    -er\tMaximum (substitution) error rate allowed in the overlap (default: 0)" << std::endl    // Maximum (substitution) error rate allowed in the overlap
		<< "    -indelRate\tMaximum indel error rate allowed in the alignment to the PacBio read (default: 1)" << std::endl
		<< "    -subRate\tMaximum substitution error rate allowed in the alignment to the PacBio read (default: 1)" << std::endl
		<< "    -noCV\tdo not use coverage depth information, specify if the reads are not random or normalized" << std::endl 
		<< "    -log\tSpecify log/output verbosity level, from ERROR, WARNING, INFO, DEBUG (default: INFO)" << std::endl 
		<< std::endl;
}

/**********************************************************************************************************************
  Parse the input arguments
 **********************************************************************************************************************/
void parseArguments(int argc, char **argv, string & inputSamFile, string & allFileName, UINT64 & minimumOverlapLength, UINT64 & hashStringLength, UINT32 & maxError, UINT32 & rubPos, float & maxErrorRate, bool & useCoverageDepth, string & outputDir, string & outputFastaName, float & indelRate, float & subRate)
{
	allFileName = "metalrec";
	minimumOverlapLength = 0;
	hashStringLength = 31;
	maxError = 0;
	maxErrorRate = 0.0;
	indelRate = 1.0;
	subRate = 1.0;
	rubPos = 10;
	useCoverageDepth = true;
	inputSamFile = "";
	outputDir = "./";
	outputFastaName = "";
	FILELog::ReportingLevel();	// Initialize the log level to the default (INFO)
	vector<string> argumentsList;
	cout << endl;
	while(argc--)
		argumentsList.push_back(*argv++);

	if(argumentsList.size() == 1)
	{
		usage();
		exit(0);
	}

	for(UINT64 i = 1; i <= argumentsList.size()-1; i++)
	{
		if (argumentsList[i] == "-h" || argumentsList[i] == "--help")
		{
			usage();
			exit(0);
		}
		else if (argumentsList[i] == "-s")
		{
			inputSamFile = argumentsList[++i];
		}
		else if (argumentsList[i] == "-f")
			allFileName = argumentsList[++i];
		else if (argumentsList[i] == "-l")
			minimumOverlapLength = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-k")
			hashStringLength = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-e")
			maxError = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-er")
			maxErrorRate = atof(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-indelRate")
			indelRate = atof(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-subRate")
			subRate = atof(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-r")
			rubPos = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-od")
			outputDir = argumentsList[++i];
		else if (argumentsList[i] == "-o")
			outputFastaName = argumentsList[++i];
		else if (argumentsList[i] == "-log"){
			try{
				FILELog::ReportingLevel() = FILELog::FromString(argumentsList[++i]);
			}
			catch(const std::exception& e)
			{
				FILE_LOG(logERROR) << e.what();
				FILELog::FromString("INFO");
			}
		}
		else if (argumentsList[i] == "-noCV")
		{
			useCoverageDepth = false;
		}
		else
		{
			usage();
			exit(0);
		}
	}

	if(minimumOverlapLength == 0)
	{
		FILE_LOG(logERROR) << "missed -l option!" << std::endl;
		usage();
		exit(0);
	}

	if ( hashStringLength >= minimumOverlapLength ) {
		FILE_LOG(logWARNING) << "hash string length is greater than the minimum overlap length, reset to minOverlap - 1" << std::endl;
		hashStringLength = minimumOverlapLength - 1;
	}
	if(maxErrorRate < 0 || maxErrorRate > 1)
	{
		FILE_LOG(logERROR) << "maximum error rate has to be between 0 and 1!" << std::endl;
		usage();
		exit(0);
	}

	if(inputSamFile.size() == 0)
	{
		FILE_LOG(logERROR) << "missed input sam file!" << std::endl;
		usage();
		exit(0);
	}
	Utils::mkdirIfNonExist(outputDir);
	if (outputFastaName.size() == 0)
		outputFastaName = outputDir + "/" + Utils::getFilename(inputSamFile) + ".fasta";

	FILE_LOG(logDEBUG) << "input sam file: " << inputSamFile;
	FILE_LOG(logDEBUG) << "output file name prefix: " << allFileName;
	FILE_LOG(logDEBUG) << "minimum overlap length: " << minimumOverlapLength;
	FILE_LOG(logDEBUG) << "hash string length: " << hashStringLength;
	FILE_LOG(logDEBUG) << "maximum error allowed in the overlap: " << maxError;
	FILE_LOG(logDEBUG) << "maximum error rate allowed in the overlap: " << maxErrorRate;
	FILE_LOG(logDEBUG) << "maximum indel error rate allowed in the alignment: " << indelRate;
	FILE_LOG(logDEBUG) << "maximum substitution error rate allowed in the alignment: " << subRate;
	FILE_LOG(logDEBUG) << "rubber length: " << rubPos;
	FILE_LOG(logDEBUG) << "use coverage depth: " << (useCoverageDepth ? "True":"False");
}

/************************
 * main function
 * **********************/
int main(int argc, char **argv)
{
	CLOCKSTART;
	/** Parse command line arguments **/
	string inputSamFile, allFileName, outputDir, outputFastaName;
	UINT64 minimumOverlapLength, hashStringLength;
	UINT32 maxError;
	float maxErrorRate, indelRate, subRate;
	UINT32 rubberPos;
	bool useCoverageDepth;
	parseArguments(argc, argv, inputSamFile, allFileName, minimumOverlapLength, hashStringLength, maxError, rubberPos, maxErrorRate, useCoverageDepth, outputDir, outputFastaName,indelRate,subRate);
	loglevel = FILELog::ReportingLevel(); // logging level in integer
	FILE_LOG(logDEBUG1) << "Parsing argument list finished";
	//cout << "logging level is " << loglevel << endl;	// For debugging only

	/** Declare variables **/
	//UINT64 counter;
	//UINT64 iteration = 0;
	/** Read sam file and store all the reads **/
	Dataset *dataSet = new Dataset(inputSamFile, minimumOverlapLength, indelRate, subRate);	// now reads the .sam file, later should be able to take the string stream TODO**
	FILE_LOG(logINFO) << "Length of the PacBio read is " << dataSet->getPacBioReadLength();
	if (dataSet->getNumberOfReads() == 0)
		FILE_LOG(logERROR) << "Data set " << inputSamFile << " has no read in it, quitting...";
	else
	{
		//UINT64 mostLikelyID = dataSet->findMostLikelyReadID(); /* find the most likely read/contig, for testing here */
		//cout << "most likely contig has ID " << mostLikelyID << " and the string is " << dataSet->getReadFromID(mostLikelyID)->getDnaStringForward() << endl;
		//dataSet->printReadsTiling(allFileName + "0_reads.tiling");
		FILE_LOG(logINFO) << "number of unique reads in dataset is " << dataSet->getNumberOfUniqueReads();
		HashTable *ht = new HashTable();
		ht->insertDataset(dataSet, hashStringLength);
		OverlapGraph *graph = new OverlapGraph(ht, minimumOverlapLength, maxError, maxErrorRate, rubberPos);
		//dataSet->saveReads(allFileName + "_reads.fasta");
		if (graph->getNumberOfEdges() == 0)
			FILE_LOG(logERROR) << "Data set " << inputSamFile << " has no edge in it, quitting...";
		else
		{
			graph->calculateFlow(outputDir + "/" + allFileName+"_flow.input", outputDir + "/" + allFileName+"_flow.output");
			FILE_LOG(logINFO) << "nodes: " << graph->getNumberOfNodes() << " edges: " << graph->getNumberOfEdges() << endl;
			graph->removeAllSimpleEdgesWithoutFlow();
			graph->simplifyGraph();
			if (loglevel > 3)
				graph->printGraph(outputDir + "/" + allFileName+"_graph.gdl", outputDir + "/" + allFileName+"_contigs.fasta");
			//graph->printGraph(outputFastaName);

			/* Use the paths instead of the unitigs for the output. First find all the paths, then pick the longest one.
			 * At this point, they are not chosen based on the likelihood, but on the lengths only. */

			vector<Edge *> paths;
			graph->findPaths(paths);
			// Save all the strings contained in the paths in fasta file
			if (loglevel > 3)
				graph->printContigs(outputDir + "/" + allFileName + "_paths.fasta", paths, false);
			// Only the longest path, as the final output
			graph->printContigs(outputFastaName, paths, true);

			delete ht;
			delete graph;
		}
	}

	delete dataSet;
	/*** For debugging, print the read containing information and the overlap information ***/
	//dataSet->printReadsTiling(allFileName + "_reads.tiling");

	CLOCKSTOP;
}

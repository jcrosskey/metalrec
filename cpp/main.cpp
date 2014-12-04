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
/*** OverlapGraph class ***/
#include "OverlapGraph.h"
int loglevel; // logging level in integer format, for different levels of verbosity

void usage();
void parseArguments(int argc, char **argv, string & inputSamFile, string & allFileName, UINT64 & minimumOverlapLength, UINT32 & maxError, float & maxErrorRate, bool & useCoverageDepth);

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
		<< "    -f\tOutput file name prefix (default: output)" << std::endl    // all output file with have this name with different extensions.
		<< "    -e\tMaximum (substitution) errors allowed in the overlap (default: 0)" << std::endl    // Maximum (substitution) errors allowed in the overlap between Illumina reads
		<< "    -er\tMaximum (substitution) error rate allowed in the overlap (default: 0)" << std::endl    // Maximum (substitution) error rate allowed in the overlap
		<< "    -noCV\tdo not use coverage depth information, specify if the reads are not random or normalized" << std::endl 
		<< "    -log\tSpecify log/output verbosity level, from ERROR, WARNING, INFO, DEBUG (default: INFO)" << std::endl 
		<< std::endl;
}

/**********************************************************************************************************************
  Parse the input arguments
 **********************************************************************************************************************/
void parseArguments(int argc, char **argv, string & inputSamFile, string & allFileName, UINT64 & minimumOverlapLength, UINT32 & maxError, float & maxErrorRate, bool & useCoverageDepth)
{
	allFileName = "output";
	minimumOverlapLength = 0;
	maxError = 0;
	maxErrorRate = 0.0;
	useCoverageDepth = true;
	inputSamFile = "";
	FILELog::ReportingLevel();	// Initialize the log level to the default (INFO)
	vector<string> argumentsList;
	//cout << "PRINTING ARGUMENTS" << endl;
	//for(UINT64 i = 0; i < argc; i++)
	//{
	//	cout << argv[i] << ' ';
	//}
	cout << endl << endl;
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
		else if (argumentsList[i] == "-e")
			maxError = atoi(argumentsList[++i].c_str());
		else if (argumentsList[i] == "-er")
			maxErrorRate = atof(argumentsList[++i].c_str());
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
	FILE_LOG(logDEBUG) << "input sam file: " << inputSamFile;
	FILE_LOG(logDEBUG) << "output file name prefix: " << allFileName;
	FILE_LOG(logDEBUG) << "minimum overlap length: " << minimumOverlapLength;
	FILE_LOG(logDEBUG) << "maximum error allowed in the overlap: " << maxError;
	FILE_LOG(logDEBUG) << "maximum error rate allowed in the overlap: " << maxErrorRate;
	FILE_LOG(logDEBUG) << "use coverage depth: " << (useCoverageDepth ? "True":"False");
}


/************************
 * main function
 * **********************/
int main(int argc, char **argv)
{
	CLOCKSTART;
	/** Parse command line arguments **/
	string inputSamFile, allFileName;
	UINT64 minimumOverlapLength;
	UINT32 maxError;
	float maxErrorRate;
	UINT32 rubberPos = 10;
	bool useCoverageDepth;
	parseArguments(argc, argv, inputSamFile, allFileName, minimumOverlapLength, maxError, maxErrorRate, useCoverageDepth);
	loglevel = FILELog::ReportingLevel(); // logging level in integer
	FILE_LOG(logDEBUG1) << "Parsing argument list finished";
	//cout << "logging level is " << loglevel << endl;	// For debugging only

	/** Declare variables **/
	//UINT64 counter;
	//UINT64 iteration = 0;
	/** Read sam file and store all the reads **/
	Dataset *dataSet = new Dataset(inputSamFile, minimumOverlapLength, maxError, maxErrorRate);	// now reads the .sam file, later should be able to take the string stream TODO**
	dataSet->saveReads(allFileName + "_reads.fasta");
	dataSet->printReadsTiling(allFileName + "_reads.tiling");
	FILE_LOG(logDEBUG4) << "number of unique reads in dataset is " << dataSet->getNumberOfUniqueReads();
	//OverlapGraph *graph = new OverlapGraph();
	OverlapGraph *graph = new OverlapGraph(dataSet, minimumOverlapLength, maxError, maxErrorRate, rubberPos);
	
	CLOCKSTOP;
}

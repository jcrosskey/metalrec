#include <mpi.h> // for mpi
#include <stdexcept> // for standard exceptions out_or_range
#include <iostream>
#include <stdio.h>
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
 * mpi_metalrec.cpp

 * MPI version of metalrec, correct reads in parallel.
 * Input file is BAM file with Illumina reads aligned to PacBio reads.
 * samtools is required
 * Naive parallelization, start the next job once there is a job done. 
 *
 * Created by JJ Chai  Thu Sep 18 09:35:28 EDT 2014. Last modified Sep 18, 2014. 
 * Copyright (c) 2014 JJ Chai (ORNL). Allrights reserved.

 ********************************************************/

#define WORKTAG    1
#define DIETAG     2
#define	MAX_FILE_NUM 20000			/*  */
using namespace std;

/* print usage */
void usage()
{
	cout << " [Usage]" << endl
		<< "  metalrec_mpi [options] -b <bamFile> -od <outDir> -samtools <samtools_path>" << endl
		<< endl

		<< " [Inputs]" << endl
		<< "   -b <bamFile> BAM file with alignment from Illumina reads to PacBio reads" << endl
		<< "   -samtools <samtools_path> path to samtools executable" << endl
		<< endl

		<< " [Outputs]" << endl
		<< "   -od <outDir> directory for fasta files for corrected PacBio sequences, maybe some log files too" << endl
		<< endl

		<< " [options]" << endl
		<< "   -h/--help Display this help message" << endl
		<< endl;
}

/* Parse command line arguments */

int initializeArguments(int argc, char ** argv,
		string & bamFile,               /* BAM file */
		string & outDir,       /* output directory */
		string & samtools_path)         /* path to samtools executable */
{
	vector<string> Arguments;
	while(argc--)
		Arguments.push_back(*argv++);

	outDir = "";
	bamFile  = "";
	samtools_path = "samtools";

	for(int i = 1; i < (int)Arguments.size(); i++)
	{
		if (Arguments[i] == "-b") {
			bamFile = Arguments.at(++i);
		}

		else if (Arguments[i] == "-od"){
			outDir = Arguments.at(++i);
		}

		else if (Arguments[i] == "-samtools"){
			samtools_path = Arguments.at(++i);
		}

		else if (Arguments[i] == "-log"){
			try{
				FILELog::ReportingLevel() = FILELog::FromString(Arguments[++i]);
			}
			catch(const std::exception& e)
			{
				FILE_LOG(logERROR) << e.what();
				FILELog::FromString("INFO");
			}
		}
		// help, print Usage
		else if (Arguments[i] == "-h" || Arguments[i] == "--help")
		{
			usage();
			return 1;
		}
		else
		{
			cerr << "Unknown option " << Arguments[i] << endl << endl;
			cerr << "Use -h/--help for usage. \n";
			return 1;
		}
	}

	if ( bamFile == "")
	{
		cerr << "missing bamFile input"<<endl;
		cerr << "use -h/--help for help" << endl;
		return 1;
	}

	if (outDir == "")
		outDir = Utils::get_cwd(); /* output directory, default: current directory */

	if (!Utils::isDirectory(outDir)) {
		mkdir(outDir.c_str(), 0777);
	}

	return 0;
}

// master job
void MasterProcess(const size_t & fileNum)
{
	int i, num_proc_spawn; // total number of processes to be spawned

	int currentWorkID, num_proc, num_slave, result;

	cout << "number of PacBio sequences is " << fileNum << endl;

	// MPI calls and initiation
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD,&num_proc);

	num_slave = num_proc - 1; /* number of slave processes */

	num_proc_spawn = (((int)fileNum <= num_slave) ? fileNum : num_slave); /* number of processes to really spawn */

	for (i = 1; i <= num_proc_spawn; i++) {
		// index of the file working on now
		currentWorkID = i - 1;
		// send the index to the slave process
		MPI_Send(&currentWorkID,1,MPI_INT, i, WORKTAG, MPI_COMM_WORLD);
	}

	if ( (int)fileNum > num_slave ) {
		// next file after each process gets a job
		currentWorkID = num_slave ;

		while (currentWorkID < (int)fileNum) { /* currentWorkID starts at 0, fileNum is the total number */
			// receive message from any process that finishes its job
			MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			//cout << "Master: received signal from " << status.MPI_SOURCE << endl;
			// send it a new job
			MPI_Send(&currentWorkID, 1, MPI_INT, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);
			//cout << "Master: send job id " << currentWorkID << " to " << status.MPI_SOURCE << endl;
			// another job is being done!
			currentWorkID++;
		}
	}

	int finish = 0;
	// when all jobs are done, send signal everybody to quit (BROADCAST?)
	for (i = 1; i <= num_slave; i++) {
		MPI_Send(&finish, 1, MPI_INT, i, DIETAG, MPI_COMM_WORLD);
	}
	//cout << "Master process is done." << endl;
}

// slave job
void SlaveProcess(const string & bamFile, const vector<string> & PacBioNames,
		const string & samtools_path, const string & outDir)
{
	MPI_Status status;
	int currentWorkID, myid;
	int finish = 0;

	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	// do jobs until master tells to stop
	while (true) {
		MPI_Recv(&currentWorkID, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		cout << myid << ": receive workID " << currentWorkID << endl;

		if (status.MPI_TAG == DIETAG) {
			break;
		}

		/* output file */
		stringstream workID, folderID;
		workID << currentWorkID;
		int folderNum = currentWorkID / MAX_FILE_NUM;
		folderID << folderNum;
		string outFile = outDir + "/" + folderID.str() + "/" +  workID.str() + ".fasta"; // output corrected fasta file
		Utils::mkdirIfNonExist(outDir + "/" + folderID.str());

		string refName = PacBioNames[currentWorkID];
		cout << myid << ": working on " << refName << endl;

		// if the corresponding sam file exists, try to correct sequence
		//string exec_cmd = samtools_path + " view " + bamFile + " " + refName + " | /chongle/shared/software/metalrec/cpp/metalrec -o " + outFile;
		//cout << exec_cmd << endl;
		string getSamCmd = samtools_path + " view " + bamFile + " " + refName;
		FILE * sam_pipe = popen(getSamCmd.c_str(), "r"); /* stream to capture the output of "samtools view bamfile refname" */
		if(!sam_pipe){
			cerr << "   *** Failed command " << getSamCmd << endl;
			perror("Error encountered"); /* If fork or pipe fails, errno is set */
		}

		/** Read sam file and store all the reads **/
		Dataset * dataSet = new Dataset(sam_pipe, 40, 0.25, 0.05);// read from the stdin stream
		//FILE_LOG(logINFO) << "Length of the PacBio read is " << dataSet->getPacBioReadLength(); /* print PacBio read length */
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
				graph->printContigs(outFile, contigEdges,true);

			}
			delete graph;
		}
		delete dataSet;
		pclose(sam_pipe);
		// signal master when done
		MPI_Send(&finish, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

	}
}


// main function
// function hmm_file prop count_file fastafile
int main(int argc, char ** argv){

	MPI_Init(&argc, &argv);

	string bamFile, samtools_path, outDir;
	vector<string> PacBioNames;
	int myid;

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	/* initialize command line arguments */
	int init = initializeArguments(argc, argv, bamFile, outDir, samtools_path);

	if(init != 0){
		MPI_Finalize();
		return 0;
	}
	else{
		double start_time, finish_time;
		/* Use samtools to get the header lines, and then get the names of the PacBio names */
		string getRefNameCmd = samtools_path + " view -H " + bamFile; 
		FILE * pipe = popen(getRefNameCmd.c_str(), "r");
		if (!pipe) 
		{
			Utils::exitWithError(" *** Failed command " + getRefNameCmd);
		}
		Utils::getRefNames(pipe, PacBioNames);
		pclose(pipe);

		size_t fileNum = PacBioNames.size();

		if (myid == 0) {
			cout << "output directory is: " << outDir << endl;
			cout << "input bam file is: " << bamFile << endl;
			start_time = MPI_Wtime();
			// display work start and time record
			cout << "Initialization succeeded " << endl;
			cout << endl
				<< "============================================================================"
				<< endl << Utils::currentDateTime() << endl
				<< " Beginning Error Correction" << endl;
			cout << " [Step 1] Error correction: Running -> " << std::flush;
			MasterProcess(fileNum);
		}

		else
		{
			SlaveProcess(bamFile, PacBioNames, samtools_path, outDir);
		}

		if( myid == 0 )
		{
			cout << " Done!" << endl;
			finish_time = MPI_Wtime();

			// display work end and time record
			cout << Utils::currentDateTime() << " Ending Corrections "<< endl
				<< "Total Elapsed Time =  "
				<< double(finish_time - start_time) << " [seconds]" << endl
				<< "============================================================================"
				<< std::endl << std::endl;
		}

		MPI_Finalize();

		return 0;
	}
}

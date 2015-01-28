#include <mpi.h> // for mpi
#include <stdexcept> // for standard exceptions out_or_range
#include <iostream>
#include <map>
#include <stdio.h>
/* check file existence and permissions */
#include <unistd.h>
#include "directoryStructure.hpp"
#include "Utils.h"

/********************************************************
 * mpi_align.cpp

 * Mapping Illumina reads to PacBio reads in parallel, using the program 
 * and parameters specified by user.
 *
 * Created by JJ Chai  Tue Jan 28 09:35:28 EDT 2015. Last modified Tue 28, 2015. 
 * Copyright (c) 2015 JJ Chai (ORNL). Allrights reserved.

 ********************************************************/

#define WORKTAG    1
#define DIETAG     2
#define	MAX_FILE_NUM 20000			/*  */
using namespace std;

/* Function declaration */
void usage();
int initializeArguments(int argc, char ** argv, string & IlluminaFiles, string & PacBioFiles, string & outDir, string & configFile, string & blasrCmd, string & samtoolsPath);
void parseConfig(const string & configFile, map<string, string> & param_map);
/* print usage */
void usage()
{
	cout << " [Usage]" << endl
		<< "  mpi_align [-od <outDir> ] -c <config.cfg> " << endl
		<< endl

		<< " [Inputs]" << endl
		<< "   -c configuration file" << endl
		<< endl

		<< " [Outputs]" << endl
		<< "   -od <outDir> directory for fasta files for corrected PacBio sequences (default: ./)" << endl
		<< endl

		<< " [options]" << endl
		<< "   -h/--help Display this help message" << endl
		<< endl;
}

/* Parse command line arguments */

int initializeArguments(int argc, char ** argv, vector<string> & IlluminaFiles, vector<string> & PacBioFiles, string & outDir, string & configFile,
		string & blasrCmd, string & samtoolsPath)
{
	vector<string> Arguments;
	while(argc--)
		Arguments.push_back(*argv++);

	outDir = "";
	configFile = "";

	for(int i = 1; i < (int)Arguments.size(); i++)
	{
		if (Arguments[i] == "-od"){
			outDir = Arguments.at(++i);
		}

		else if (Arguments[i] == "-c"){
			configFile = Arguments.at(++i);
		}


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

	if (outDir == "")
		outDir = Utils::get_cwd(); /* output directory, default: current directory */

	if (!Utils::isDirectory(outDir)) {
		mkdir(outDir.c_str(), 0777);
	}

	if (configFile.length() == 0)
	{
		cerr << "Missing config file\n";
		return 1;
	}
	else
		parseConfig(configFile, IlluminaFiles, PacBioFiles, blasrCmd, samtoolsPath);
	return 0;
}

void parseConfig(const string & configFile, vector<string> & IlluminaFiles, vector<string> & PacBioFiles, string & blasrCmd, string & samtoolsPath)
{
	map<string, string> param_map;

	ifstream configFilePointer;
	configFilePointer.open(configFile.c_str());
	if(!configFilePointer.is_open())
		cerr << "Unable to open file: " << configFile << ", Use default parameter value instead.";
	else
	{
		string line;
		while(!configFilePointer.eof())
		{
			getline(configFilePointer, line);
			if (line.at(0)!='#' && line.length() != 0)
			{
				istringstream is_line(line);
				string key;
				if (getline(is_line, key, '='))
				{
					string value;
					if (getline(is_line, value))
					{
						try
						{ param_map.at(key) = value;}
						catch(const exception & e)
						{FILE_LOG(logERROR) << e.what();}
					}
				}
			}
		}
	}
	configFilePointer.close();
	blasrCmd = param_map["blasrCmd"];
	samtoolsPath = param_map["samtoolsPath"];
	IlluminaDir = param_map["IlluminaDir"];
	PacBioDir = param_map["PacBioDir"];

	DirectoryStructure illumina_dir(IlluminaDir);
	illumina_dir.setPattern(".fasta");
	illumina_dir.getFiles(IlluminaFiles);
	DirectoryStructure pacbio_dir(PacBioDir);
	pacbio_dir.setPattern(".fasta");
	pacbio_dir.getFiles(PacBioFiles);
}
// master job
void MasterProcess(const size_t & IlluminaFileNum, const size_t & PacBioFileNum)
{
	int i, num_proc_spawn; // total number of processes to be spawned

	int currentWorkID, num_proc, num_slave, result;

	int fileNum = IlluminaFileNum * PacBioFileNum; /* Total number of BLASR jobs to run */
	cout << "number of mapping jobs is " << fileNum << endl;

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
void SlaveProcess( const vector<string> & IlluminaFiles,const vector<string> & PacBioFiles,
		const string & blasrCmd, const string & samtoolsPath, const string & outDir)
{
	MPI_Status status;
	int currentWorkID, myid;
	int finish = 0;

	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	// do jobs until master tells to stop
	while (true) {
		MPI_Recv(&currentWorkID, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		//cout << myid << ": receive workID " << currentWorkID << endl;

		if (status.MPI_TAG == DIETAG) {
			cout << myid << ": receive workID " << currentWorkID << endl;
			break;
		}

		/* output file */
		stringstream IlluminaID, PacBioID;
		IlluminaID << (currentWorkID % PacBioFiles.size());
		PacBioID << (currentWorkID / PacBioFiles.size());

		string PacBioFile = PacBioFiles.at(currentWorkID / PacBioFiles.size()); // PacBio File name for this job
		string IlluminaFile = IlluminaFiles.at(currentWorkID % PacBioFiles.size()); // Illumina file name for this job

		string outFile = outDir + "/" + IlluminaID.str() + "_" +  PacBioID.str() + ".bam"; // output corrected fasta file

		cout << myid << ": working on " << IlluminaID.str() + "_" + PacBioID.str() << endl;

		string wholeCmd = blasrCmd + " " + IlluminaFile + " " + PacBioFile  + " 2> /dev/null | " + \
				  samtoolsPath + " view -@ 16 -bT " + PacBioFile + " | " + \
				  samtoolsPath + " sort -@ 16 -o " + outFile;

		res = system(wholeCmd.c_str());
		if (res != 0)
		{
			cerr << "Fail command " << wholeCmd << endl;
		}
		// signal master when done
		MPI_Send(&finish, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
}


// main function
// function hmm_file prop count_file fastafile
int main(int argc, char ** argv){

	MPI_Init(&argc, &argv);

	string blasrCmd, samtoolsPath, configFile, outDir;

	vector<string> PacBioFiles, IlluminaFiles;
	int myid;

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	/* initialize command line arguments */
	int init = initializeArguments(argc, argv, IlluminaFiles, PacBioFiles, outDir, configFile, blasrCmd, samtoolsPath);

	if(init != 0){
		MPI_Finalize();
		return 0;
	}
	else{
		size_t IlluminaFileNum = IlluminaFiles.size();
		size_t PacBioFileNum = PacBioFiles.size();

		if (myid == 0) {
			cout << endl
				<< "============================================================================"
				<< endl << Utils::currentDateTime() << endl
				<< " Beginning Error Correction" << endl;
			cout << " [Step 1] Error correction: Running -> " << std::flush;
			MasterProcess(IlluminaFileNum, PacBioFileNum);
		}

		else
		{
			SlaveProcess(IlluminaFiles,PacBioFiles, blasrCmd, samtoolsPath, outDir);
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

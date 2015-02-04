#include <mpi.h> // for mpi
#include <stdexcept> // for standard exceptions out_or_range
#include <iostream>
#include <string>
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
#define	MAX_FILE_NUM 20000
using namespace std;

/* Function declaration */
void usage();
int initializeArguments(int argc, char ** argv, string & IlluminaFiles, string & PacBioFiles, string & outDir, string & configFile, string & blasrCmd, string & samtoolsPath);
void parseConfig(const string & configFile, string & IlluminaDir, string & PacBioDir, vector<string> & IlluminaFiles, vector<string> & PacBioFiles, string & blasrCmd, string & samtoolsPath, string & threads);
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

int initializeArguments(int argc, char ** argv, string & IlluminaDir, string & PacBioDir, vector<string> & IlluminaFiles, vector<string> & PacBioFiles, 
		string & outDir, string & configFile,
		string & blasrCmd, string & samtoolsPath, string & threads)
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
		parseConfig(configFile, IlluminaDir, PacBioDir, IlluminaFiles, PacBioFiles, blasrCmd, samtoolsPath, threads);
	return 0;
}

void parseConfig(const string & configFile, string & IlluminaDir, string & PacBioDir, vector<string> & IlluminaFiles, vector<string> & PacBioFiles, string & blasrCmd, string & samtoolsPath, string & threads)
{
	threads = "1";
	ifstream configFilePointer;
	configFilePointer.open(configFile.c_str());
	if(!configFilePointer.is_open())
		cerr << "Unable to open file: " << configFile << ", Use default parameter value instead.";
	else
	{
		string line, key, value;
		while(!configFilePointer.eof())
		{
			getline(configFilePointer, line);
			//cout << "line: " <<  line << " length: " << line.length() << endl;
			if (line.length() != 0 && line.at(0)!='#')
			{
				istringstream is_line(line);
				if (getline(is_line, key, '='))
				{
					key = key.substr(key.find_first_not_of(" "), key.find_last_not_of(" ") - key.find_first_not_of(" ")+1);
					//cout << key << "=";
					if (getline(is_line, value))
					{
						value = value.substr(value.find_first_not_of(" "), value.find_last_not_of(" ") - value.find_first_not_of(" ") + 1);
						//cout << value << endl;
						if (key.compare("blasrCmd")==0)
							blasrCmd = value;
						else if (key.compare("samtoolsPath")==0)
							samtoolsPath = value;
						else if (key.compare("IlluminaDir")==0)
							IlluminaDir = value;
						else if (key.compare("PacBioDir")==0)
							PacBioDir = value;
						else if (key.compare("Threads") ==0)
							threads = value;
					}
				}
			}
		}
	}
	configFilePointer.close();

	DirectoryStructure illumina_dir(IlluminaDir);
	illumina_dir.setPattern(".fasta");
	illumina_dir.getFiles(IlluminaFiles);
	illumina_dir.setPattern(".fastq");
	illumina_dir.getFiles(IlluminaFiles);
	DirectoryStructure pacbio_dir(PacBioDir);
	pacbio_dir.setPattern(".fasta");
	pacbio_dir.getFiles(PacBioFiles);
	pacbio_dir.setPattern(".fastq");
	pacbio_dir.getFiles(PacBioFiles);
	//cout << "Illumina directory is " << IlluminaDir << endl;
	//cout << "PacBio directory is " << PacBioDir << endl;
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
void SlaveProcess(const string & IlluminaDir, const string & PacBioDir, const vector<string> & IlluminaFiles,const vector<string> & PacBioFiles,
		const string & blasrCmd, const string & samtoolsPath, const string & outDir, const string & threads)
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
		IlluminaID << (currentWorkID % IlluminaFiles.size());
		PacBioID << (currentWorkID / IlluminaFiles.size());

		string PacBioFile = PacBioDir + "/" + PacBioFiles.at(currentWorkID / IlluminaFiles.size()); // PacBio File name for this job
		string IlluminaFile = IlluminaDir + "/" + IlluminaFiles.at(currentWorkID % IlluminaFiles.size()); // Illumina file name for this job

		string outFile = outDir + "/" + IlluminaID.str() + "_" +  PacBioID.str(); // output corrected fasta file
		string samFile = outFile + ".sam";
		string bamFile = outFile + ".bam";

		cout << myid << ": working on " << IlluminaID.str() + "_" + PacBioID.str() << endl;

		string wholeCmd = blasrCmd + " -nproc " + threads + " -out " + samFile +  " " + IlluminaFile + " " + PacBioFile + " 2> /dev/null && " + \
				  samtoolsPath + " view -F 4 -@ " + threads + " -bT " + PacBioFile + " " + samFile + " | " + \
				  samtoolsPath + " sort -@ " + threads + " -o " + bamFile + " -T " + outFile + "_tmp" + " && " + \
				  samtoolsPath + " index " + bamFile;

		//cout << myid << ": whole command is \"" << wholeCmd <<  "\" " << endl;
		cout << myid << ": " << IlluminaFiles.at(currentWorkID % IlluminaFiles.size()) << "\t" << PacBioFiles.at(currentWorkID / IlluminaFiles.size()) <<  "\t" << IlluminaID.str() + "_" +  PacBioID.str() + ".bam" <<  endl;
		int res = system(wholeCmd.c_str());
		if (res != 0)
		{
			cerr << "Fail command for " << IlluminaID.str() + "_" +  PacBioID.str() << endl;
		}
		// signal master when done
		MPI_Send(&finish, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
}


// main function
// function hmm_file prop count_file fastafile
int main(int argc, char ** argv){

	MPI_Init(&argc, &argv);

	string IlluminaDir, PacBioDir, blasrCmd, samtoolsPath, configFile, outDir, threads;

	vector<string> PacBioFiles, IlluminaFiles;
	int myid;

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	/* initialize command line arguments */
	int init = initializeArguments(argc, argv, IlluminaDir, PacBioDir, IlluminaFiles, PacBioFiles, outDir, configFile, blasrCmd, samtoolsPath, threads);

	if(init != 0){
		MPI_Finalize();
		return 0;
	}
	else{
		double start_time, finish_time;
		size_t IlluminaFileNum = IlluminaFiles.size();
		size_t PacBioFileNum = PacBioFiles.size();

		if (myid == 0) {
			cout << "Number of illumina files: " << IlluminaFiles.size() << endl;
			cout << "Number of PacBio files: " << PacBioFiles.size() << endl;
			cout << "blasr command is " << blasrCmd << endl;
			cout << "samtools path is " << samtoolsPath << endl;
			cout << endl
				<< "============================================================================"
				<< endl << Utils::currentDateTime() << endl
				<< " Beginning Error Correction" << endl;
			start_time = MPI_Wtime();
			cout << " [Step 1] Error correction: Running -> " << std::flush;
			MasterProcess(IlluminaFileNum, PacBioFileNum);
		}

		else
		{
			SlaveProcess(IlluminaDir, PacBioDir,IlluminaFiles,PacBioFiles, blasrCmd, samtoolsPath, outDir, threads);
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

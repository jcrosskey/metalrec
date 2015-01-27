#include <mpi.h> // for mpi
#include <stdexcept> // for standard exceptions out_or_range
#include <iostream>
#include <map>
#include <stdio.h>
/* check file existence and permissions */
#include <unistd.h>

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
//		DirectoryStructure sam_dir(samDir);
//		sam_dir.setPattern(".sam");
//		sam_dir.getFiles(samFilenames);
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
	DirectoryStructure pacbio_dir(IlluminaDir);
	illumina_dir.setPattern(".fasta");
	illumina_dir.getFiles(IlluminaFiles);
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
		const string & samtools_path, const string & outDir,
		const UINT64 & minimumOverlapLength, const UINT64 & hashStringLength,
		const UINT32 & maxError, const UINT32 &rubberPos,
		const float & indelRate, const float & subRate, const float & maxErrorRate)
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
		/* popen returns NULL if fork or pipe calls fail, or if it cannot allocate memory */
		FILE * sam_pipe = popen(getSamCmd.c_str(), "r"); /* stream to capture the output of "samtools view bamfile refname" */
		if(!sam_pipe){
			cerr << "   *** Failed command " << getSamCmd << endl;
			perror("Error encountered in popen()"); /* If fork or pipe fails, errno is set */
		}

		else{
			/** Read sam file and store all the reads **/
			Dataset * dataSet = new Dataset(sam_pipe, minimumOverlapLength, indelRate, subRate);// read from the stdin stream
			int exit_status = pclose(sam_pipe); /* pclose() waits for the associated process to terminate and returns the exit status of the command as returned by wait4; 
							       returns -1 if wait4 returns an error, or some other error is detected. */
			if (exit_status == -1)
			{
				perror("Error encountered in pclose()");
			}
			//FILE_LOG(logINFO) << "Length of the PacBio read is " << dataSet->getPacBioReadLength(); /* print PacBio read length */
			if (dataSet->getNumberOfReads() <= 1)
				FILE_LOG(logERROR) << "Data set has no more than 1 read in it, quitting...";

			else /* If the dataset has some reads in it, build hash table and overlap graph next */
			{
				FILE_LOG(logINFO) << "number of unique reads in dataset is " << dataSet->getNumberOfUniqueReads();
				HashTable *ht = new HashTable();
				ht->insertDataset(dataSet, hashStringLength);
				OverlapGraph *graph = new OverlapGraph(ht, minimumOverlapLength, maxError, maxErrorRate, rubberPos);
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
		}
		// signal master when done
		MPI_Send(&finish, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
}


// main function
// function hmm_file prop count_file fastafile
int main(int argc, char ** argv){

	MPI_Init(&argc, &argv);

	string bamFile, samtools_path, configFile, outDir, allFileName;
	UINT64 minimumOverlapLength, hashStringLength;
	UINT32 maxError, rubberPos;
	float indelRate, subRate, maxErrorRate;
	map<string, string> param_map;          /* mapping from argument key to arg value, initialization */
	param_map["allFileName"]= "metalrec";
	param_map["minimumOverlapLength"] = "40"; 
	param_map["hashStringLength"] = "10";
	param_map["maxError"]="0";
	param_map["maxErrorRate"] = "0.0";
	param_map["rubberPos"] ="10";
	param_map["indelRate"] = "0.25"; 
	param_map["subRate"] = "0.05";

	vector<string> PacBioNames;
	int myid;

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	/* initialize command line arguments */
	int init = initializeArguments(argc, argv, bamFile, outDir, configFile, samtools_path);

	if(init != 0){
		MPI_Finalize();
		return 0;
	}
	else{
		param_map["outDir"] = outDir;
		param_map["samtools_path"] = samtools_path;
		if (configFile.length() != 0)
		{
			FILE_LOG(logINFO) << "config file is specified, parsing the file to get new parameter values";
			parseConfig(configFile,param_map);
		}

		minimumOverlapLength = (UINT64) Utils::stringToUnsignedInt(param_map["minimumOverlapLength"]);
		hashStringLength = (UINT64) Utils::stringToUnsignedInt(param_map["hashStringLength"]);
		maxError = (UINT32) Utils::stringToUnsignedInt(param_map["maxError"]);
		rubberPos = (UINT32) Utils::stringToUnsignedInt(param_map["rubberPos"]);
		indelRate = Utils::stringToFloat(param_map["indelRate"]);
		maxErrorRate = Utils::stringToFloat(param_map["maxErrorRate"]);
		subRate = Utils::stringToFloat(param_map["subRate"]);
		allFileName = param_map["allFileName"];
		outDir = param_map["outDir"];
		samtools_path = param_map["samtools_path"];

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
			FILE_LOG(logINFO) << "samtools path: " << samtools_path;
			FILE_LOG(logINFO) << "bam file: " << bamFile;
			FILE_LOG(logINFO) << "output directory: " << outDir;
			FILE_LOG(logINFO) << "allFileName: " << allFileName;
			FILE_LOG(logINFO) << "minimum overlap length: " << minimumOverlapLength;
			FILE_LOG(logINFO) << "hash string length: " << hashStringLength;
			FILE_LOG(logINFO) << "maxError: " << maxError;
			FILE_LOG(logINFO) << "max error rate : " << maxErrorRate;
			FILE_LOG(logINFO) << "indel rate: " << indelRate;
			FILE_LOG(logINFO) << "substitution rate: " << subRate;
			cout << endl
				<< "============================================================================"
				<< endl << Utils::currentDateTime() << endl
				<< " Beginning Error Correction" << endl;
			cout << " [Step 1] Error correction: Running -> " << std::flush;
			MasterProcess(fileNum);
		}

		else
		{
			SlaveProcess(bamFile, PacBioNames, samtools_path, outDir, minimumOverlapLength, hashStringLength, maxError, rubberPos, 
					indelRate, subRate, maxErrorRate);
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

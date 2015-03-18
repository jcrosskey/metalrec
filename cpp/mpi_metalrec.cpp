#include <mpi.h>
#include "metalrec.h"
int loglevel;                                   /* level of log information(verbosity) */

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
#define	MAX_FILE_NUM 20000
using namespace std;

/* Function declaration */
void usage();
int initializeArguments(int argc, char ** argv, string & outDir, string & configFile, string & PacBioReadName, string & PacBio_file);
void parseConfig(const string & configFile, map<string, string> & param_map);
/* print usage */
void usage()
{
	cout << " [Usage]" << endl
		<< "  metalrec [-od <outDir> ] -c <config.cfg>] " << endl
		<< endl

		<< " [Inputs]" << endl
		<< "   -c configuration file with non-default parameter values" << endl
		<< endl

		<< " [Outputs]" << endl
		<< "   -od <outDir> directory for fasta files for corrected PacBio sequences (default: ./)" << endl
		<< endl

		<< " [options]" << endl
		<< "   -name <pacbio_name> PacBio read to error correct (default: Try to correct all)" << endl
		<< "   -f <PacBio_file> File including names of PacBio reads to correct, one name on each line (default: none)" << endl
		<< "   -h/--help Display this help message" << endl
		<< endl;
}

/* Parse command line arguments */

int initializeArguments(int argc, char ** argv,               /* BAM file */
		string & outDir,       /* output directory */
		string & configFile,            /* Configuration file, required */
		string & PacBioReadName, string & PacBio_file)        /* A particular PacBio read to correct, or a file including all the PacBio read names */ 
{
	vector<string> Arguments;
	while(argc--)
		Arguments.push_back(*argv++);

	outDir = "";
	configFile = "";
	PacBioReadName = "";
	PacBio_file = "";

	for(int i = 1; i < (int)Arguments.size(); i++)
	{
		if (Arguments[i] == "-od"){
			outDir = Arguments.at(++i);
		}

		else if (Arguments[i] == "-c"){
			configFile = Arguments.at(++i);
		}

		else if (Arguments[i] == "-name"){
			PacBioReadName = Arguments.at(++i);
		}

		else if (Arguments[i] == "-f"){
			PacBio_file = Arguments.at(++i);
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

		else if (Arguments[i] == "-h" || Arguments[i] == "--help")
		{
			usage();
			return 1;
		}
		else
		{
			cerr << "Unknown option " << Arguments[i] << endl;
			cerr << "Use -h/--help for usage. \n";
			return 1;
		}
	}

	if (PacBioReadName.length()>0 && PacBio_file.length() > 0)
	{
		FILE_LOG(logWARNING) << "Both long read name and file with long reads specified, correct only 1 by the name";
		PacBio_file = "";

	}
	if (configFile.length() == 0)
	{
		cerr << "Missing config file...\nUse option -h/--h to see help.\n";
		return 1;
	}
	if (outDir == "")
		outDir = Utils::get_cwd(); /* output directory, default: current directory */

	if (!Utils::isDirectory(outDir)) {
		mkdir(outDir.c_str(), 0777);
	}

	return 0;
}

void parseConfig(const string & configFile, map<string, string> & param_map)
{
	CLOCKSTART;
	ifstream configFilePointer;
	configFilePointer.open(configFile.c_str());
	if(!configFilePointer.is_open())
		FILE_LOG(logERROR) << "Unable to open file: " << configFile << ", Use default parameter value instead.";
	else
	{
		string line;
		while(!configFilePointer.eof())
		{
			getline(configFilePointer, line);
			if (line.length() != 0 && line.at(0)!='#')
			{
				istringstream is_line(line);
				string key;
				if (getline(is_line, key, '='))
				{
					key = key.substr(key.find_first_not_of(" "), key.find_last_not_of(" ") - key.find_first_not_of(" ")+1);
					string value;
					if (getline(is_line, value))
					{
						value = value.substr(value.find_first_not_of(" "), value.find_last_not_of(" ") - value.find_first_not_of(" ") + 1);
						param_map[key] = value;
					}
				}
			}
		}
	}

	configFilePointer.close();
	CLOCKSTOP;
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
void SlaveProcess(const vector<string> & bamFiles, const vector<string> & PacBioNames, const vector<UINT16> & PacBioLengths,
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
		int folderNum = currentWorkID / MAX_FILE_NUM; /* Corrected sequences are in different folders since number of files in a folder is limited */
		folderID << folderNum;          /* folder name */
		string prefixName = workID.str(); /* prefix for all output files */
		string outputDir = outDir + "/" + folderID.str() + "/"; /* output directory, including the numbered folder */
		Utils::mkdirIfNonExist(outputDir);

		string refName = PacBioNames[currentWorkID];
		UINT16 PacBioLength = PacBioLengths[currentWorkID];
		cout << myid << ": working on " << refName << endl;
		metalrec(bamFiles, refName, PacBioLength, prefixName, samtools_path, outputDir, minimumOverlapLength, hashStringLength, maxError, rubberPos, 
				indelRate, subRate, maxErrorRate);

		// signal master when done
		MPI_Send(&finish, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
}


// main function
// function hmm_file prop count_file fastafile
int main(int argc, char ** argv){

	MPI_Init(&argc, &argv);

	string samtools_path, configFile, outDir, allFileName, PacBioName, PacBio_file;
	int myid, num_proc;                               /* worker ID in mpi setting */
	UINT64 minimumOverlapLength, hashStringLength;
	UINT32 maxError, rubberPos;
	vector<string> PacBioNames;
	vector<UINT16> PacBioLengths;
	float indelRate, subRate, maxErrorRate;
	double start_time, finish_time;
	map<string, string> param_map;          /* mapping from argument key to arg value, initialization */

	/* Set default values */
	param_map["allFileName"]= "metalrec";
	param_map["minimumOverlapLength"] = "40"; 
	param_map["hashStringLength"] = "10";
	param_map["maxError"]="0";
	param_map["maxErrorRate"] = "0.0";
	param_map["rubberPos"] ="10";
	param_map["indelRate"] = "0.25"; 
	param_map["subRate"] = "0.05";
	param_map["samtools_path"] = "samtools";


	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&num_proc);

	/* initialize command line arguments */
	int init = initializeArguments(argc, argv, outDir, configFile, PacBioName, PacBio_file);
	loglevel = FILELog::ReportingLevel(); // logging level in integer

	if(init != 0){
		MPI_Finalize();
		return 0;
	}
	else{
		FILE_LOG(logINFO) << "Parsing the config file to get parameter values";
		parseConfig(configFile,param_map);

		vector<string> bamFiles = Utils::StringToVector(param_map["bamfiles"],' ');

		if (param_map["bamDir"].length()!=0) /* Attach the bam directory to the bam files to get full path */
		{
			FILE_LOG(logINFO) << "bam directory is: " << param_map["bamDir"];
			for (size_t i = 0 ; i < bamFiles.size(); i++)
				bamFiles.at(i) = param_map["bamDir"] + "/" + bamFiles.at(i);
		}
		if ( bamFiles.size() == 0 )
		{
			Utils::exitWithError("No input bam files!\n");
		}

		param_map["outDir"] = outDir;
		minimumOverlapLength = (UINT64) Utils::stringToUnsignedInt(param_map["minimumOverlapLength"]);
		hashStringLength = (UINT64) Utils::stringToUnsignedInt(param_map["hashStringLength"]);
		maxError = (UINT32) Utils::stringToUnsignedInt(param_map["maxError"]);
		rubberPos = (UINT32) Utils::stringToUnsignedInt(param_map["rubberPos"]);
		indelRate = stof(param_map["indelRate"]);
		maxErrorRate = stof(param_map["maxErrorRate"]);
		subRate = stof(param_map["subRate"]);
		allFileName = param_map["allFileName"];
		outDir = param_map["outDir"];
		samtools_path = param_map["samtools_path"];

		/* If no specific PacBio read(s) is requested from command line, 
		 * Use samtools to get the header lines, and then get the names of all PacBio names */
		if ( PacBioName.length() == 0 && PacBio_file.length() == 0)
		{
			string getRefNameCmd = samtools_path + " view -H " + bamFiles[0]; 
			FILE * pipe = popen(getRefNameCmd.c_str(), "r");
			if (!pipe) 
			{
				Utils::exitWithError(" *** Failed command " + getRefNameCmd);
			}
			Utils::getRefNames(pipe, PacBioNames, PacBioLengths);
			pclose(pipe);
		}
		else if (PacBio_file.length() > 0)
		{
			Utils::saveLinesToVec(PacBio_file, PacBioNames);
		}
		else                            /* Only 1 read to correct */
			PacBioNames.push_back(PacBioName);

		size_t fileNum = PacBioNames.size();

		if (myid == 0) {
			cout << "Number of PacBio reads is: " << PacBioNames.size() << endl;
			cout << "File with PacBio names is: " << PacBio_file << endl;
			cout << "samtools path: " << samtools_path << endl;
			cout << "number of bam files: " << bamFiles.size() << endl;
			cout << "bamfiles: " << param_map["bamfiles"] << endl;
			cout << "output directory: " << outDir << endl;
			cout << "number of pacbio reads to work on: " << PacBioNames.size() << endl;
			cout << "allFileName: " << allFileName << endl;
			cout << "minimum overlap length: " << minimumOverlapLength << endl;
			cout << "hash string length: " << hashStringLength << endl;
			cout << "maxError: " << maxError << endl;
			cout << "max error rate : " << maxErrorRate << endl;
			cout << "indel rate: " << indelRate << endl;
			cout << "substitution rate: " << subRate << endl;
			cout << "============================================================================" << endl;
			cout << "Beginning Error Correction" << Utils::currentDateTime() << endl;

			start_time = MPI_Wtime();
			if (PacBioNames.size() == 1)
			{
				metalrec(bamFiles, PacBioNames.at(0), 65000, allFileName, samtools_path, outDir, minimumOverlapLength, hashStringLength, maxError, rubberPos, 
						indelRate, subRate, maxErrorRate);

			}
			if(num_proc > 1)        /* If there is more than 1 mpi processes */
			{
				MasterProcess(fileNum);
			}
			else                    /* If there is only 1 mpi process, correct sequences in serial instead */
			{
				for ( size_t j = 0; j < PacBioNames.size(); j++)
				{
					FILE_LOG(logINFO) << "Read " << PacBioNames.at(j);
					string prefixName = Utils::intToString(j);
					metalrec(bamFiles, PacBioNames.at(j), PacBioLengths.at(j),prefixName, samtools_path, outDir, minimumOverlapLength, hashStringLength, maxError, rubberPos, 
							indelRate, subRate, maxErrorRate);
				}

				cout << " Done!"<< Utils::currentDateTime() << endl;
				cout << "============================================================================" << endl;
			}

		}

		else
		{
			SlaveProcess(bamFiles, PacBioNames, PacBioLengths, samtools_path, outDir, minimumOverlapLength, hashStringLength, maxError, rubberPos, 
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

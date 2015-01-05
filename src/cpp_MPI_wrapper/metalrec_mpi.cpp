#include <mpi.h> // for mpi
#include <stdexcept> // for standard exceptions out_or_range
#include <iostream>
#include <stdio.h>
#include "utils.hpp" // for useful utilities - namespace
#include "directoryStructure.hpp" // for searching files in directory - class

/********************************************************
 * metalrec_mpi.cpp
 
 * MPI version of metalrec.py, correct reads in parallel.
 * All the PacBio reads are stored in a directory(fasta_dir) with one fasta file for each read.
 * All the sam files are stored in another directory (sam_dir) with one sam file for each read, unless there was no read mapped to the PacBio read.
 * Naive parallelization, start the next job once there is a job done. Call python system command.
 *
 * Created by JJ Chai  Thu Sep 18 09:35:28 EDT 2014. Last modified Sep 18, 2014. 
 * Copyright (c) 2014 JJ Chai (ORNL). Allrights reserved.
 
 ********************************************************/

#define WORKTAG    1
#define DIETAG     2

using namespace std;

/* print usage */
void usage()
{
    cout << " [Usage]" << endl
    << "  metalrec_mpi [options] -f <readNameFile> -fd <fastaDir> -sd <samDir> -od <outDir> -path <py_path> " << endl
    << endl
    
    << " [Inputs]" << endl
    << "   -f <readNameFile> File including base names of read/sam files" << endl
    << "   -sd <samDir> directory with sam files including mapping results for each PacBio sequence" << endl
    << "   -p <exec_path> full path of error correction executable" << endl
    << endl
    
    << " [Outputs]" << endl
    << "   -od <outDir> directory for fasta files for corrected PacBio sequences, maybe some log files too" << endl
    << endl
    
    << " [options]" << endl
    << "   -h/--help Display this help message" << endl
    << endl;
}

/* Parse command line arguments */
/* align_fasta [options] -d <fastaDir> -o <outDir> -s <py_path> */

int initializeArguments(int argc, char ** argv,
                         vector<string> & samFilenames, // input fasta files
                         vector<string> & outFilenames, //output fasta files for corrected sequences
			 string & samDir, // output directory name
                         string & exec_path)  // error correction executable file's path
{
	vector<string> Arguments;
	while(argc--)
		Arguments.push_back(*argv++);

	string baseNameFile = ""; /* File including base names of read/sam files */
	string outDir = ""; /* output directory */

	samDir  = ""; /* input directory for the sam files of PacBio sequences */
	exec_path = ""; /* python commmand's prefix (including path and some options, no input and output) */

	for(int i = 1; i < (int)Arguments.size(); i++)
	{
		// File including base names of read/sam files, one name on each row
		if (Arguments[i] == "-f") {
			try{
			    baseNameFile = Arguments.at(++i);
			}
			catch(const out_of_range& oor){
				baseNameFile = "";
			}
		}

		// input directory for the sam files of PacBio sequences
		else if (Arguments[i] == "-sd"){
			try{
			    samDir = Arguments.at(++i);
			}
			catch(const out_of_range& oor){
				samDir = "";
			}
		}

		// output directory for the fasta files of corrected PacBio sequences
		else if (Arguments[i] == "-od"){
			try{
			    outDir = Arguments.at(++i);
			}
			catch(const out_of_range& oor){
				outDir = "";
			}
		}

		// python commmand's prefix (including path and some options, no input and output)
		else if (Arguments[i] == "-p"){
			try{
			    exec_path = Arguments.at(++i);
			}
			catch(const out_of_range& oor){
				exec_path = "";
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
	/* check required arguments
	* 1. fastaDir; 2. samDir 3. error correction executable path */
	//if ( ((baseNameFile == "") && (fastaDir== "")) || (exec_path == "") || (samDir == ""))
	if ( ((baseNameFile == "") && (samDir== "")) || (exec_path == "") )
	{
		cerr << "missing necessary parameter(s)"<<endl;
		cerr << "use -h/--help for help" << endl;
		return 1;
	}

	if (outDir == "")
	    outDir = Utils::get_cwd(); /* output directory, default: current directory */

	//cout << "python command's path is: " << py_cmd_path << endl;
	//cout << "omega executable file's path is: " << omega_cmd_path << endl;
	//cout << "output directory is: " << outDir << endl;
	//cout << "input fasta directory is: " << fastaDir << endl;
	//cout << "input sam file directory is: " << samDir << endl;

	// first check if starting from a file, instead of scanning a directory
	if (baseNameFile != "") {
		ifstream baseIn;
		baseIn.open(baseNameFile.c_str());
		string line;
		while (getline(baseIn, line)){
			samFilenames.push_back(line + ".sam");
		}
		baseIn.close();
	}
	// find all the sam files in the input directory, ends with .sam
	else if (samDir != "") {
		DirectoryStructure sam_dir(samDir);

		sam_dir.setPattern(".sam");
		sam_dir.getFiles(samFilenames);

		//cout << "total number of sam files is " << samFilenames.size() << endl;
	}

	/* number of fasta files found in the input directory
	* if no fasta files are found in the directories, exit
	*/
	int fileNum = (int) samFilenames.size();
	if (fileNum == 0)
	{
		cerr << "no sam file in the input directory"<<endl;
		return 1;
	}

	// full path of the input files and the output files
	for (int i = 0; i < fileNum; i++)
	{
		string outFilename = outDir + "/" + Utils::getFilename2(samFilenames.at(i)) + ".fasta"; // output file name for the corrected PacBio sequence
		samFilenames.at(i) = samDir + "/" + samFilenames.at(i); // full path for the input sam files
		outFilenames.push_back(outFilename);
	}

	// create the output directory if it does not exist already
	if (!Utils::isDirectory(outDir)) {
		//cout << " create directory: " << outDir << endl;
		mkdir(outDir.c_str(), 0777);
	}

	return 0;
}

// master job
void MasterProcess(const int & fileNum)
{
    int i, num_proc_spawn; // total number of processes to be spawned
    
    int currentWorkID, num_proc, num_slave, result;

    cout << "number of PacBio sequences is " << fileNum << endl;
    
    // MPI calls and initiation
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD,&num_proc);
    
    num_slave = num_proc - 1; /* number of slave processes */
    
    num_proc_spawn = ((fileNum <= num_slave) ? fileNum : num_slave); /* number of processes to really spawn */
    
    for (i = 1; i <= num_proc_spawn; i++) {
        // index of the file working on now
        currentWorkID = i - 1;
        // send the index to the slave process
        MPI_Send(&currentWorkID,1,MPI_INT, i, WORKTAG, MPI_COMM_WORLD);
    }
    
    if ( fileNum > num_slave ) {
        // next file after each process gets a job
        currentWorkID = num_slave ;
        
        while (currentWorkID < fileNum) { /* currentWorkID starts at 0, fileNum is the total number */
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
void SlaveProcess(const vector<string> & samFilenames, const vector<string> & outFilenames,
                  const string & exec_path, const string & samDir)
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
        
	string samFile = samFilenames[currentWorkID];	// sam file to work on
	string basename = Utils::getFilename2(samFile); // base name of the fasta file, same for sam file
	string outFile = outFilenames[currentWorkID]; // output corrected fasta file

	cout << myid << ": working on " << basename << endl;
	string output_prefix = basename + "_metalrec"; // prefix for omega's output files

	// if the corresponding sam file exists, try to correct sequence
	if(Utils::isFileExist(samFile)) {
		string metalrec_cmd =  exec_path + " -s " + samFile + " -o " + outFile + " -od " + samDir + " -f " + output_prefix + " -l 40 -k 15 -er 0.000 -log ERROR -indelRate 0.25 -subRate 0.1";
		//cout << metalrec_cmd << endl;
		int res = system(metalrec_cmd.c_str());
		if(res != 0){
			cout << "   *** Failed command " << metalrec_cmd << endl;
		}
	}
	// else just print message
	else{
		cout << "   *** samFile " << samFile << " does not exist" << endl;
	}
        
        // signal master when done
        MPI_Send(&finish, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        
    }
    return;
    
}


// main function
// function hmm_file prop count_file fastafile
int main(int argc, char ** argv){
    
    MPI_Init(&argc, &argv);

    string exec_path, samDir;
    vector<string> samFilenames, outFilenames;
    int myid;
    
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);

    /* some simple test for MPI
    * int world_size;
    * MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    * char processor_name[MPI_MAX_PROCESSOR_NAME];
    * int name_len;
    * MPI_Get_processor_name(processor_name, &name_len);
    * cout << "Hello world from processor " << processor_name << ", " << myid << " out of " << world_size << " processors\n";
    * cout << argc << " arguments on the command line in total:\n";
    * for (int i = 0; i < argc; i++){
    *         cout << argv[i] << "\t";
    * }
    * cout << endl;
    */
	
    /* initialize command line arguments */
    int init = initializeArguments(argc, argv, samFilenames, outFilenames, samDir, exec_path);

    if(init != 0){
	    MPI_Finalize();
	    return 0;
    }
    else{
	    double start_time, finish_time;
	    int fileNum = (int) samFilenames.size();
	    
	    if (myid == 0) {
		start_time = MPI_Wtime();
		// display work start and time record
		cout << "Initialization succeeded " << endl;
		cout << endl
		<< "============================================================================"
		<< endl << Utils::currentDateTime() << endl
		<< " Beginning Error Correction" << endl
		<< " [Step 1] Looking for sam files: Running -> " << std::flush;
	    }
	    
	    if ( myid == 0 ) {
		cout << " Done!" << endl;
		cout << " [Step 2] Error correction: Running -> " << std::flush;
		MasterProcess(fileNum);
	    }
	    
	    else
	    {
		SlaveProcess(samFilenames, outFilenames, exec_path, samDir);
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

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
    << "   -fd <samDir> directory with fasta files for each PacBio sequence" << endl
    << "   -ec <ec_path> full path of error correction executable file" << endl
    << "   -pick <pick_path> full path of contig picking executable file" << endl
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
                         vector<string> & samFilenames, // input fasta files
			 string & outDir,       /* output directory */
			 string & samDir, // PacBio sam file directory name
			 string & fastaDir, // PacBio fasta directory name
			 string & pick_path, // contig picking executable file's path
                         string & ec_path)  // error correction executable file's path
{
	vector<string> Arguments;
	while(argc--)
		Arguments.push_back(*argv++);

	string baseNameFile = ""; /* File including base names of read/sam files */
	outDir = ""; /* output directory */

	samDir  = ""; /* input directory for the sam files of PacBio sequences */
	fastaDir = ""; /* directory with fasta files for each PacBio sequence */
	ec_path = ""; /* full path of error correction executable file */
	pick_path = ""; /* full path of contig picking executable file */

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
		else if (Arguments[i] == "-fd"){
			try{
			    fastaDir = Arguments.at(++i);
			}
			catch(const out_of_range& oor){
				fastaDir = "";
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

		// contig picking executable file's path
		else if (Arguments[i] == "-pick"){
			try{
			    pick_path = Arguments.at(++i);
			}
			catch(const out_of_range& oor){
				pick_path = "";
			}
		}
		
		// error correction executable file's path
		else if (Arguments[i] == "-ec"){
			try{
			    ec_path = Arguments.at(++i);
			}
			catch(const out_of_range& oor){
				ec_path = "";
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
	* 1. fastaDir; 2. samDir 3. error correction executable path 4. contig picking executable path */
	//if ( ((baseNameFile == "") && (fastaDir== "")) || (exec_path == "") || (samDir == ""))
	if ( ((baseNameFile == "") && (samDir== "")) || (ec_path == "") || (fastaDir == "") || (pick_path == "") )
	{
		cerr << "missing necessary parameter(s)"<<endl;
		cerr << "use -h/--help for help" << endl;
		return 1;
	}

	if (outDir == "")
	    outDir = Utils::get_cwd(); /* output directory, default: current directory */

//	cout << "error correction executable's path is: " << ec_path << endl;
//	cout << "contig picking executable's path is: " << pick_path << endl;
//	cout << "output directory is: " << outDir << endl;
//	cout << "input fasta directory is: " << fastaDir << endl;
//	cout << "input sam file directory is: " << samDir << endl;

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
		samFilenames.at(i) = samDir + "/" + samFilenames.at(i); // full path for the input sam files
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
void SlaveProcess(const vector<string> & samFilenames,
                  const string & ec_path, const string & pick_path, const string & samDir, const string & fastaDir, const string & outDir)
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
        
	/* input files */
	string samFile = samFilenames[currentWorkID];	// sam file to work on
	string basename = Utils::getFilename2(samFile); // base name of the fasta file, same for sam file
	string fastaFile = fastaDir + "/" + basename + ".fasta"; // fasta file

	/* output file */
	string outFile = outDir + "/" + basename + ".fasta"; // output corrected fasta file

	/* intermediate files */
	string contigFile = outDir + "/" + basename + ".contigs.fasta"; // contigs from miniassembly
	string contigSamFile = outDir + "/" + basename + ".contigs.sam"; /* alignment of contigs to PacBio read */

	cout << myid << ": working on " << basename << endl;
	string output_prefix = basename + "_metalrec"; // prefix for omega's output files

	// if the corresponding sam file exists, try to correct sequence
	if(Utils::isFileExist(samFile)) {
		string metalrec_cmd =  ec_path + " -s " + samFile + " -o " + contigFile + " -od " + samDir + " -f " + output_prefix + " -l 40 -k 15 -log ERROR";
		//cout << metalrec_cmd << endl;
		int res = system(metalrec_cmd.c_str());
		if(res != 0){
			cout << "   *** Failed command " << metalrec_cmd << endl;
		}
		else if (Utils::isFileExist(contigFile))
		{
			/* First check how many contigs are assembled from the Illumina reads */
			ifstream contig_in;
			contig_in.open(contigFile.c_str());
			if(!contig_in.is_open())
			{
				cout << "Unable to open file: " << contigFile << endl;
			}
			else
			{
				int num_reads = 0;
				string s;
				while(getline(contig_in, s))
				{
					if ( s[0] == ">"[0])
						num_reads = num_reads + 1;
					if (num_reads > 1)
						break;
				}
				contig_in.close();
				if(num_reads == 1)
				{
					cout << "mini assembly generated only 1 contig, No need to choose.\nDone.\n";
					rename(contigFile.c_str(), outFile.c_str()); // rename ec's output to output File
				}
				else
				{
					cout << "mini assembly generated more than 1 contig\n";

					// Step 4: if there is more than 1 output contig, use blasr to align them to the PacBio sequence, and choose one
					string blasr_cmd = "blasr " + contigFile + " " + fastaFile + " -noSplitSubreads -sam -clipping soft -out " + contigSamFile + " 2> /dev/null"; /* stderr not output to stream */
					blasr_cmd += " | " + pick_path + " -log ERROR -o " + outFile;

					res = system(blasr_cmd.c_str());
					if(res != 0){
						cout << "   *** Failed command " << blasr_cmd << endl;
					}
				}
			}
		}
	}
	// else just print message
	else{
		cout << "   *** samFile " << samFile << " does not exist" << endl;
	}
        
        // signal master when done
        MPI_Send(&finish, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        
    }
}


// main function
// function hmm_file prop count_file fastafile
int main(int argc, char ** argv){
    
    MPI_Init(&argc, &argv);

    string samDir, fastaDir, outDir, pick_path, ec_path;
    vector<string> samFilenames;
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
    int init = initializeArguments(argc, argv, samFilenames, outDir, samDir, fastaDir, pick_path, ec_path);

    if(init != 0){
	    MPI_Finalize();
	    return 0;
    }
    else{
	    double start_time, finish_time;
	    int fileNum = (int) samFilenames.size();
	    
	    if (myid == 0) {
		    cout << "error correction executable's path is: " << ec_path << endl;
		    cout << "contig picking executable's path is: " << pick_path << endl;
		    cout << "output directory is: " << outDir << endl;
		    cout << "input fasta directory is: " << fastaDir << endl;
		    cout << "input sam file directory is: " << samDir << endl;
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
		SlaveProcess(samFilenames, ec_path, pick_path, samDir, fastaDir, outDir);
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

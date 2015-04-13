#include "ec.h"
int loglevel;                                   /* level of log information(verbosity) */

/********************************************************
 * metalrec.cpp

 * serial version of metalrec.
 * Input file is BAM file(s) with Illumina reads aligned to PacBio reads.
 * samtools is required
 *
 * Created by JJ Chai  Thu Sep 18 09:35:28 EDT 2014. Last modified Sep 18, 2014. 
 * Copyright (c) 2014 JJ Chai (ORNL). Allrights reserved.

 ********************************************************/

#define WORKTAG    1
#define DIETAG     2
#define	MAX_FILE_NUM 20000	/* Limit of number of files in each folder (on some old file system) */
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
		Utils::exitWithError("Missing config file...\nUse option -h/--h to see help.\n");
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

// main function
// function hmm_file prop count_file fastafile
int main(int argc, char ** argv){
	string samtools_path, configFile, outDir, allFileName, PacBioName, PacBio_file;
	UINT64 minimumOverlapLength, hashStringLength;
	UINT32 maxError, rubberPos;
	UINT16 minPacBioLength;
	vector<string> PacBioNames;
	vector<UINT16> PacBioLengths;
	float indelRate, insRate, delRate, subRate, minPercentInLR, maxErrorRate;
	map<string, string> param_map;          /* mapping from argument key to arg value, initialization */

	/* Set default values */
	param_map["allFileName"]= "metalrec";
	param_map["minimumOverlapLength"] = "40"; 
	param_map["hashStringLength"] = "10";
	param_map["maxError"]="0";
	param_map["maxErrorRate"] = "0.0";
	param_map["rubberPos"] ="10";
	param_map["indelRate"] = "0.25"; 
	param_map["insRate"] = "0.15"; 
	param_map["delRate"] = "0.10"; 
	param_map["subRate"] = "0.05";
	param_map["minPercentInLR"] = "0.80";
	param_map["samtools_path"] = "samtools";
	param_map["minPacBioLength"] = "1000";

	/* initialize command line arguments */
	int init = initializeArguments(argc, argv, outDir, configFile, PacBioName, PacBio_file);
	loglevel = FILELog::ReportingLevel(); // logging level in integer

	if(init != 0){
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
		insRate = stof(param_map["insRate"]);
		delRate = stof(param_map["delRate"]);
		maxErrorRate = stof(param_map["maxErrorRate"]);
		subRate = stof(param_map["subRate"]);
		minPercentInLR = stof(param_map["minPercentInLR"]);
		allFileName = param_map["allFileName"];
		outDir = param_map["outDir"];
		samtools_path = param_map["samtools_path"];
		minPacBioLength = stoi(param_map["minPacBioLength"]);


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
			for(size_t i = 0; i < PacBioNames.size(); i++){
				PacBioLengths.push_back(Utils::getPacBioLength(PacBioNames.at(i)));
			}
		}
		else                            /* Only 1 read to correct */
		{
			PacBioNames.push_back(PacBioName);
			PacBioLengths.push_back(Utils::getPacBioLength(PacBioName));
			FILE_LOG(logINFO) << "PacBio Read Length is: " << PacBioLengths.at(0);
		}

		FILE_LOG(logINFO) << "Number of PacBio reads in first bam file is: " << PacBioNames.size();
		FILE_LOG(logINFO) << "samtools path: " << samtools_path;
		FILE_LOG(logINFO) << "number of bam files: " << bamFiles.size();
		FILE_LOG(logINFO) << "bamfiles: " << param_map["bamfiles"];
		FILE_LOG(logINFO) << "output directory: " << outDir;
		FILE_LOG(logINFO) << "number of pacbio reads to work on: " << PacBioNames.size();
		FILE_LOG(logINFO) << "allFileName: " << allFileName;
		FILE_LOG(logINFO) << "minimum overlap length: " << minimumOverlapLength;
		FILE_LOG(logINFO) << "hash string length: " << hashStringLength;
		FILE_LOG(logINFO) << "maxError: " << maxError;
		FILE_LOG(logINFO) << "max error rate : " << maxErrorRate;
		FILE_LOG(logINFO) << "indel rate: " << indelRate;
		FILE_LOG(logINFO) << "insertion rate: " << insRate;
		FILE_LOG(logINFO) << "deletion rate: " << delRate;
		FILE_LOG(logINFO) << "substitution rate: " << subRate;
		FILE_LOG(logINFO) << "minimum length percent of short read in long read: " << minPercentInLR;
		FILE_LOG(logINFO) << "rubber length for alignment buffer: " << rubberPos;
		FILE_LOG(logINFO) << "============================================================================";
		FILE_LOG(logINFO) << "Beginning Error Correction";

		if (PacBioNames.size() == 1)    /* When there is only 1 PacBio read to correct, use the specified prefix name */
		{
			FILE_LOG(logINFO) << "Read " << PacBioNames.at(0);
			ec(bamFiles, PacBioNames.at(0), PacBioLengths.at(0), allFileName, samtools_path, outDir, minimumOverlapLength, hashStringLength, maxError, rubberPos, 
					indelRate, insRate, delRate, subRate, maxErrorRate, minPercentInLR, minPacBioLength);
		}
		else{
			for ( size_t j = 0; j < PacBioNames.size(); j++) /* Otherwise the prefix has to be changed so that each PacBio read has different prefix */
			{
				FILE_LOG(logINFO) << "Read " << PacBioNames.at(j);
				string prefixName = Utils::intToString(j);
				ec(bamFiles, PacBioNames.at(j), PacBioLengths.at(j), prefixName, samtools_path, outDir, minimumOverlapLength, hashStringLength, maxError, rubberPos, 
						indelRate, insRate, delRate, subRate, maxErrorRate, minPercentInLR, minPacBioLength);
			}
		}

		FILE_LOG(logINFO) << " Done!";
		FILE_LOG(logINFO) << "============================================================================";
		return 0;
	}
}

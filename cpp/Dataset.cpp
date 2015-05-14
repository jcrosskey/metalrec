/*
 * Dataset.cpp
 *
 * Created on: Tue Nov 11 12:45:42 EST 2014
 * Author: JJ Chai
 */

#include "Dataset.h"
#include <limits>
#include <unistd.h>

/**********************************************************************************************************************
 * Function to compare two reads by their start mapping coordinate. Used for sorting.
 * If two reads have same start mapping coordinate, sort by their lengths (i.e. by their ending coordinates).
 * i.e. longer read will be placed before the shorter read when they have the same starting coordinates.
 **********************************************************************************************************************/
bool compareReads (Read *read1, Read *read2)
{
	if (read1 -> getStartCoord() != read2 -> getStartCoord())
		return read1->getStartCoord() < read2->getStartCoord();
	else
		return read1->getReadLength() > read2->getReadLength();	// When the staring coordinates are the same, the longer one comes first
}

/**********************************************************************************************************************
  Default constructor
 **********************************************************************************************************************/
Dataset::Dataset(void)
{
	// Initialize the variable values.
	totalBps = 0;
	numberOfReads = 0;
	numberOfUniqueReads = 0;
	numberOfNonContainedReads = 0;
	minimumOverlapLength = 0;
	shortestReadLength = 0XFFFFFFFFFFFFFFFF;
	longestReadLength = 0X0000000000000000;
	PacBioReadName = "";
	PacBioReadLength = 0;
	subRate = 0.01;
	indelRate = 0.15;
	insRate = 0.10;
	delRate = 0.05;
	percentInLR = 0.80;
	reads = new vector<Read *>;
}

/**********************************************************************************************************************
  Copy constructor
 **********************************************************************************************************************/
Dataset::Dataset(const Dataset & D)
{
	// Initialize the variable values.
	totalBps = D.totalBps;
	numberOfReads = D.numberOfReads;
	numberOfUniqueReads = D.numberOfUniqueReads;
	numberOfNonContainedReads = D.numberOfNonContainedReads;
	minimumOverlapLength = D.minimumOverlapLength;
	shortestReadLength = D.shortestReadLength;
	longestReadLength = D.longestReadLength;
	PacBioReadName = D.PacBioReadName;
	PacBioReadLength = D.PacBioReadLength;
	subRate = D.subRate;
	indelRate = D.indelRate;
	insRate = D.insRate;
	delRate = D.delRate;
	percentInLR = D.percentInLR;
	reads = new vector<Read *>;
	reads->resize((D.reads)->size());
	for(size_t k = 0; k < reads->size(); k++){
		Read * read1 = new Read;
		*read1 = (*((D.reads)->at(k)));
		reads->at(k) = read1;
	}
}

/**********************************************************************************************************************
  Copy assignment
 **********************************************************************************************************************/
Dataset & Dataset::operator= (const Dataset & D)
{
	// Initialize the variable values.
	totalBps = D.totalBps;
	numberOfReads = D.numberOfReads;
	numberOfUniqueReads = D.numberOfUniqueReads;
	numberOfNonContainedReads = D.numberOfNonContainedReads;
	minimumOverlapLength = D.minimumOverlapLength;
	shortestReadLength = D.shortestReadLength;
	longestReadLength = D.longestReadLength;
	PacBioReadName = D.PacBioReadName;
	PacBioReadLength = D.PacBioReadLength;
	subRate = D.subRate;
	indelRate = D.indelRate;
	insRate = D.insRate;
	delRate = D.delRate;
	percentInLR = D.percentInLR;

	delete this->reads;
	reads = new vector<Read *>;
	reads->resize((D.reads)->size());
	for(size_t k = 0; k < reads->size(); k++){
		Read * read1 = new Read;
		*read1 = (*((D.reads)->at(k)));
		reads->at(k) = read1;
	}
	return *this;
}

/**********************************************************************************************************************
 * Another constructor, from BLASR generated sam file, do not use BamAlignmentRecord class
 **********************************************************************************************************************/
Dataset::Dataset(const string & inputSamFile, UINT64 minOverlap, const float & indel_rate, const float & sub_rate, const float & ins_rate, const float & del_rate, const float & percent_inLR)
{
	/** get reads from sam file , instead of FILE object that identifies and contains information to control a stream**/
	FILE * samIn;
	samIn = fopen(inputSamFile.c_str(), "r");
	if (samIn == NULL)
		perror("Error opening file");
	else
	{
		Dataset(samIn, minOverlap, indel_rate, sub_rate, ins_rate, del_rate, percent_inLR);
		fclose(samIn);// Close the sam file
	}
}

/**********************************************************************************************************************
 * Another constructor, from BLASR generated sam file, do not use BamAlignmentRecord class
 **********************************************************************************************************************/
Dataset::Dataset(FILE * inputSamStream, UINT64 minOverlap, const float & indel_rate, const float & sub_rate, const float & ins_rate, const float & del_rate, const float & percent_inLR)
{
	CLOCKSTART;
	// Initialize the variables.
	numberOfUniqueReads = 0;
	numberOfReads = 0;                      /* number of good reads in the dataset */
	numberOfNonContainedReads = 0;
	shortestReadLength = 0XFFFFFFFFFFFFFFFF;
	longestReadLength = 0X0000000000000000;
	reads = new vector<Read *>;
	minimumOverlapLength = minOverlap;
	indelRate = indel_rate;
	subRate = sub_rate;
	insRate = ins_rate;
	delRate = del_rate;
	percentInLR = percent_inLR;
	PacBioReadName = "";

	AddDataset(inputSamStream);
	CLOCKSTOP;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  AddDataset
 *  Description:  Add more aligned reads to this LR from another input stream
 * =====================================================================================
 */
bool Dataset::AddDataset(FILE * inputSamStream)
{
	CLOCKSTART;
	UINT64 goodReads = 0, badReads = 0;
	/** get reads from sam file **/
	string line = "";
	char buffer[BUFFER_SIZE];
	while(1)
	{
		if(fgets(buffer, BUFFER_SIZE, inputSamStream)==NULL)
		{
			//perror("In function Dataset, Encountered in function fgets");
			if (ferror(inputSamStream))
			{FILE_LOG(logERROR) << "error reading dataset";}
			break;
		}
		line += buffer;
		if (line.at(line.size() - 1) == '\n')        /* A whole line has been read, processing it */
		{
			line.erase(line.size()-1, 1);        /* Delete last character from line, which is EOL */
			FILE_LOG(logDEBUG4) << "Line: " << line;
			if (line[0] == '@')
			{
				FILE_LOG(logDEBUG3) << "header line";
				if ( line.substr(0,3).compare("@SQ") ==0) /* SQ line */
				{
					size_t LNPos = line.find("LN:");
					size_t nextTabPos = line.find("\t", LNPos+3);
					PacBioReadLength = Utils::stringToUnsignedInt(line.substr(LNPos+3,nextTabPos-LNPos-3)); /* length of the PacBio read */
					LNPos = line.find("SN:");
					nextTabPos = line.find("\t", LNPos+3);
					PacBioReadName = line.substr(LNPos+3,nextTabPos-LNPos-3);
				}
			}
			else	// not header line, alignment record
			{
				FILE_LOG(logDEBUG4) << "Read line with length: " << line.length();
				if (line.length() !=0)
				{
					Read *r = new Read(line);
					if ( PacBioReadName.length() == 0 )
						PacBioReadName = r->getRefName();
					FILE_LOG(logDEBUG4) << "Scanned read: " << r->getReadName();
					if ( r->isReadGood(indelRate,insRate,delRate, subRate, PacBioReadLength, percentInLR) && testRead(r->getDnaStringForward()))
					{
						UINT32 len = r->getReadLength();
						if (len > longestReadLength)
							longestReadLength = len;
						if (len < shortestReadLength)
							shortestReadLength = len;
						Read *goodRead = new Read(*r);
						reads->push_back(goodRead);
						numberOfReads++;
						goodReads++;
						totalBps += r->getAlignedLength();
					}
					else
					{
						badReads++;
						if (!(r->isReadGood(indelRate,insRate,delRate, subRate, PacBioReadLength, percentInLR)))
							FILE_LOG(logDEBUG3) << "Read " << r->getReadName() << " has too many errors";
						else
							FILE_LOG(logDEBUG3) << "Read is not valid";
					}
					delete r;
				}
			}
			line = "";              /* Set line to empty again */
		}
	}
	FILE_LOG(logINFO) << "In this data set: ";
	FILE_LOG(logINFO) << "Shortest read length: " << setw(5) << shortestReadLength;
	FILE_LOG(logINFO) << "Longest read length: " << setw(5) << longestReadLength;
	FILE_LOG(logINFO) << "Number of good reads: " << setw(5) << goodReads;
	FILE_LOG(logINFO) << "Number of bad reads: " << setw(5) << badReads;
	FILE_LOG(logINFO) << "Total number of reads: " << setw(5) << badReads + goodReads;
	FILE_LOG(logINFO) << "Total number of reads in the dataset so far: " << setw(5) << reads->size();
	FILE_LOG(logINFO) << "Total number of aligned bps in the dataset so far: " << setw(5) << totalBps;
	CLOCKSTOP;
	return true;
}

bool Dataset::finalize(void)
{
	sortReads();
	removeDupicateReads();	// Remove duplicated reads for the dataset.
	numberOfNonContainedReads = numberOfUniqueReads;
	return true;
}

/**********************************************************************************************************************
  Default destructor
 **********************************************************************************************************************/
Dataset::~Dataset(void)
{
	// Free the memory used by the dataset.
	for(UINT64 i = 0; i < reads->size(); i++)
	{
		delete reads->at(i);
		reads->at(i) = NULL;
	}
	delete reads;
	reads = NULL;
}

/****************************************************
 * Sort reads by their start mapping coordinate
 * **************************************************/
void Dataset::sortReads(void)
{
	CLOCKSTART;
	sort(reads->begin(),reads->end(), compareReads);	// Sort the reads lexicographically.
	CLOCKSTOP;
}

/**********************************************************************************************************************
 * This function returns the number of unique reads and assign ID to the reads.
 * Reads are duplicated iff their mapping start coordinates are the same and their strings are the same.
 * Same reads mapped to different positions are considered as different reads.
 **********************************************************************************************************************/
bool Dataset::removeDupicateReads(void)
{
	CLOCKSTART;
	UINT64 j = 0;	// Unique read counter
	Read *temp;
	if ( reads->size() == 0 )
		numberOfUniqueReads = 0;
	else
	{
		for(UINT64 i = 0; i < reads->size(); i++)	// Move the unique reads to the top of the sorted list. Store the frequency of the duplicated reads.
		{
			if(reads->at(j)->getStartCoord()!= reads->at(i)->getStartCoord() \
					|| reads->at(j)->getDnaStringForward() != reads->at(i)->getDnaStringForward())	// Two reads have different start mapping coordinates, or different strings
			{
				j++;	// Number of unique reads increases by 1
				temp = reads->at(j);	// Save the read that was just checked
				reads->at(j) = reads->at(i);	// Switch reads.at(i) and reads.at(j)
				reads->at(i) = temp;
			}
		}
		numberOfUniqueReads = j+1;
		for(UINT64 i = 0 ; i < reads->size(); i++) 	// Assing ID to the reads.
		{
			if(i < getNumberOfUniqueReads())
				reads->at(i)->setReadID(i+1);	// Read's ID, 1 based, 0 (superReadID default for non-contained reads)
			else
				delete reads->at(i); 	// Free the unused reads.
		}
		reads->resize(numberOfUniqueReads);	//Resize the vector of reads.
	}
	FILE_LOG(logINFO) <<"Number of unique reads: " << numberOfUniqueReads;
	CLOCKSTOP;
	return true;
}

/**********************************************************************************************************************
  Returns true if the read is longer than the overlap offset and contains only {A,C,G,T} 
 **********************************************************************************************************************/
bool Dataset::testRead(const string & readDnaString)
{
	if ( readDnaString.length() < minimumOverlapLength)
	{
		FILE_LOG(logDEBUG3) << "Read's length is smaller than the minimumOverlapLength " << readDnaString;
		return false;
	}
	/* make sure there is no characters other than ACGT */
	if(readDnaString.find_first_not_of("ACGT") != string::npos){
		FILE_LOG(logDEBUG3) << "Read has characters other than ACGT " << readDnaString;
		return false;
	}
	/* We might not want to check the composition of the read any more, since this might not be what the user wants 
	 * However, homopolymer run is more problematic than simply high composition of a single nucleotide.*/
	UINT64 cnt[4] = {0,0,0,0};              /* count of the 4 nucleotides */
	UINT64 homopolymerRuns[4] = {0,0,0,0};  /* longest homopolymer runs for each nucleotide */
	UINT64 longestRunPos[4] = {0,0,0,0};    /* position where the longest homopolymer runs end */
	if(readDnaString.length() > 0){
		UINT64 i = 0;
		while(i < readDnaString.length()){
			char base = readDnaString[i];
			cnt[(base >> 1) & 0X03]++;
			UINT64 currentRun = 1;
			i++;
			while(readDnaString[i] == base && i < readDnaString.length()){
				cnt[(base >> 1) & 0X03]++;
				currentRun++;
				i++;
			}
			if (homopolymerRuns[(base >> 1) & 0X03] < currentRun)
			{
				homopolymerRuns[(base >> 1) & 0X03] = currentRun;
				longestRunPos[(base>>1) & 0X03] = i-1;
			}
		}
	}
/* 	for(UINT64 i = 0; i < readDnaString.length(); i++) // Count the number of A's, C's , G's and T's in the string.
 * 	{
 * 		if(readDnaString[i]!= 'A' && readDnaString[i] != 'C' && readDnaString[i] != 'G' && readDnaString[i] != 'T')
 * 		{
 * 			FILE_LOG(logDEBUG3) << "Read has characters other than ACGT " << readDnaString;
 * 			return false;
 * 		}
 * 		cnt[(readDnaString[i] >> 1) & 0X03]++; // Least significant 2nd and 3rd bits of ASCII value used here
 * 	}
 */
	/* TODO: put the percentage and homopolymer run length in the config file */
	UINT64 threshold = readDnaString.length()*.8;	// 80% of the length.
	if(cnt[0] >= threshold || cnt[1] >= threshold || cnt[2] >= threshold || cnt[3] >= threshold)
	{
		FILE_LOG(logDEBUG3) << "More than 80%% of the read string has the same character " << readDnaString;
		return false;	// If 80% bases are the same base.
	}
	/* If a read has homopolymer run longer than 15 bps, consider it as bad */
	if(homopolymerRuns[0] >= 15 || homopolymerRuns[1] >= 15 || homopolymerRuns[2] >= 15 || homopolymerRuns[3] >= 15)
	{
		FILE_LOG(logDEBUG3) << "Long homopolymer run(s) > 15 bps exists in the read " << readDnaString;
		return false;	// If 80% bases are the same base.
	}

	return true;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getNumberOfNonContainedReads
 *  Description:  as function name suggested
 * =====================================================================================
 */
UINT64 Dataset::getNumberOfNonContainedReads(void)
{
	UINT64 N = 0;
	for(UINT64 i = 1; i <= numberOfUniqueReads; i++)
	{
		if(getReadFromID(i)->superReadID == 0)
			N++;
	}
	return N;
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getBpsOfNonContainedReads
 *  Description:  Find total number of bases of the nonContained reads
 * =====================================================================================
 */
UINT64 Dataset::getBpsOfNonContainedReads(void)
{
	UINT64 N = 0;
	for(UINT64 i = 1; i <= numberOfUniqueReads; i++)
	{
		if(getReadFromID(i)->superReadID == 0)
			N += getReadFromID(i)->getReadLength();
	}
	return N;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getCoveredInfo
 *  Description:  Find number of bases covered by Illumina reads, and number of regions covered
 * =====================================================================================
 */
bool Dataset::getCoveredInfo(INT64 & coveredBases, INT64 & coveredRegions)
{
	CLOCKSTART;
	coveredBases = 0;
	coveredRegions = 0;
	vector<INT64> startCoords, endCoords;
	INT64 istart = -(numeric_limits<INT64>::max()-10);
	INT64 iend = -(numeric_limits<INT64>::max()-10);
	for(UINT64 i = 1; i <= numberOfUniqueReads; i++){
		Read * read = getReadFromID(i);
		if (read->superReadID == 0){
			if(read->getStartCoord() > iend) /* this read does not overlap with intervals already found */
			{
				if(iend > 0){         /* an interval ends here */
					startCoords.push_back(istart<0?0:istart);
					endCoords.push_back(iend<PacBioReadLength?iend:PacBioReadLength);
				}
				istart = read->getStartCoord();
				iend = read->getEndCoord();
			}
			else if (read->getEndCoord() > iend){ /* only if this read has a chance to cover more bases */
				iend = read->getEndCoord();
			}
		}
	}
	startCoords.push_back(istart<0?0:istart);
	endCoords.push_back(iend<PacBioReadLength?iend:PacBioReadLength);
	for(size_t i = 0; i < startCoords.size(); i++){
		coveredBases += (endCoords.at(i) - startCoords.at(i));
	}
	coveredRegions = startCoords.size();
	FILE_LOG(logINFO) << coveredBases << " bps are covered by short reads, split into " << coveredRegions << " intervals.";
	FILE_LOG(logDEBUG2) << "Coordinates of these intervals are: ";
	if(loglevel >= 3){
		for(size_t i = 0; i < startCoords.size(); i++)
			cout << "(" << startCoords.at(i) << ", " << endCoords.at(i) << ") ";
		cout << endl;
	}
	CLOCKSTOP;
	return true;
}

/**********************************************************************************************************************
  Search a read in the dataset using binary search
 **********************************************************************************************************************/
vector<Read *> Dataset::getReadsFromCoord(const INT32 & coord)	// Find a read in the database given the start mapping coord. Uses binary search in the list of reads.
{
	vector<Read *> readsFound;
	UINT64 min = 0, max = getNumberOfUniqueReads()-1, mid, midLast;
	int comparator;
	while (max >= min) 	// At first search for the forward string.
	{
		mid = (min + max) / 2; 	// Determine which subarray to search.
		comparator = reads->at(mid)->getStartCoord() - coord;
		if(comparator == 0)
			break;                  /* Found a read with the specified starting coordinate */
		else if (comparator < 0) 	// Change min index to search upper subarray.
			min = mid + 1;
		else if (comparator > 0) 	// Change max index to search lower subarray.
			max = mid - 1;
	}
	if (comparator == 0)
	{
		midLast = mid;
		while (mid > 0){
			if(reads->at(mid-1)->getStartCoord() == coord)
				mid--;
			else
				break;
		}
		while (midLast < getNumberOfUniqueReads()-1){
			if(reads->at(midLast+1)->getStartCoord() == coord)
				midLast++;
			else
				break;
		}
		for(UINT64 i = mid; i <= midLast; i++)
		{
			readsFound.push_back(reads->at(i));
		}
	}
	return readsFound;
}

/**********************************************************************************************************************
  This function returns a read from its ID.
 **********************************************************************************************************************/
Read * Dataset::getReadFromID(UINT64 ID)
{
	if( ID < 1 || ID > numberOfUniqueReads)	// ID outside the range.
	{
		stringstream ss;
		ss << "ID " << ID << " out of bound.";
		string s = ss.str();
		MYEXIT(s);
	}
	return reads->at(ID - 1);
}

/**********************************************************************************************************************
  This function saves the starting coordinate sorted unique reads in a file. Used for debugging only.
 **********************************************************************************************************************/
void Dataset::saveReads(string fileName, bool writeContained)
{
	CLOCKSTART;
	ofstream outputFile;
	outputFile.open(fileName.c_str());
	if(!outputFile.is_open())
		MYEXIT("Unable to open file: " + fileName);
	for(UINT64 i = 1; i <= numberOfUniqueReads; i++)
	{
		Read * read1 = getReadFromID(i);
		if(read1->superReadID!=0) 
		{
			if (writeContained)
			{
				outputFile << ">" << i <<  " # " << read1->getReadName() \
					<< " # Start coord: " << read1->getStartCoord() \
					<< " # End coord: " << read1->getEndCoord() \
					<< " # Length: " << read1->getReadLength() \
					<< " # Subs: " << read1->getNumOfSubstitutionsInRead() \
					<< " # Ins: " << read1->getNumOfInsertionsInRead() \
					<< " # Del: " << read1->getNumOfDeletionsInRead() \
					<< " # leftClip: " << read1->getLeftClip() \
					<< " # rightClip: " << read1->getRightClip() \
					<< " # Contained in "  << read1->superReadID \
					<< "\n" << read1->getDnaStringForward() << endl;
			}
		}
		else
		{
			outputFile << ">" << i <<  " # " << read1->getReadName() \
				<< " # Start coord: " << read1->getStartCoord() \
				<< " # End coord: " << read1->getEndCoord() \
				<< " # Length: " << read1->getReadLength() \
				<< " # Subs: " << read1->getNumOfSubstitutionsInRead() \
				<< " # Ins: " << read1->getNumOfInsertionsInRead() \
				<< " # Del: " << read1->getNumOfDeletionsInRead() \
				<< " # leftClip: " << read1->getLeftClip() \
				<< " # rightClip: " << read1->getRightClip() \
				<< " # Noncontained, contains: ";
			for(size_t i = 0; i < read1->getContainedReadIDs()->size(); i++)
				outputFile << read1->getContainedReadIDs()->at(i) << ", ";
			outputFile << "\n" << read1->getDnaStringForward() << endl;
		}
	}
	outputFile.close();
	CLOCKSTOP;
}

/**********************************************************************************************************************
 * This function prints the reads in tiling format, for checking the overlap alignment, and debugging.
 **********************************************************************************************************************/
void Dataset::printReadsTiling(string fileName)	// Print all the reads in tiling format. Used for checking the overlap (debugging)
{
	CLOCKSTART;
	ofstream outputFile;
	outputFile.open(fileName.c_str());
	if(!outputFile.is_open())
		MYEXIT("Unable to open file: " + fileName);
	Read * read0 = getReadFromID(1);	// Read with leftmost mapping start coordinate
	int LeftMostCoord = read0->getStartCoord();
	if (LeftMostCoord < 0 )
	{
		string overHang(0-LeftMostCoord, '.');
		outputFile << overHang;
	}
	string PacBioSeq(1000,'*');
	outputFile << PacBioSeq << endl;
	for(UINT64 i = 1; i <= numberOfUniqueReads; i++)
	{
		Read * read1 = getReadFromID(i);
		if (read1->getStartCoord() > LeftMostCoord)
		{
			string OffsetString(read1->getStartCoord() - LeftMostCoord, '.');
			outputFile << OffsetString;
			//cout << i << "\t" << read1->getStartCoord() << endl;
		}
		outputFile << read1->getDnaStringForward() << endl;
	}
	outputFile.close();
	CLOCKSTOP;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getWeight
 *  Description:  get weight of an edge
 * =====================================================================================
 */
double Dataset::getWeight(Edge * e)
{
	/* coverage depth */
	UINT64 numOfOverlapOffsets = e->getListOfOverlapOffsets()->size()+1; 

	return ((double)e->getOverlapOffset() + (double)numOfOverlapOffsets*PacBioReadLength/(double)numberOfUniqueReads - getSubsOnEdge(e)*20);
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  getSubsOnEdge
 *  Description:  Find expected number of substitutions on an edge
 * =====================================================================================
 */
double Dataset::getSubsOnEdge(Edge *e)
{
	double substitutionsInEdge = 0.0;
	Read *sourceRead = e->getSourceRead();
	if(e->getListOfReads()->size() == 0){   /* If this edge is a simple edge, then just 1 overlap offset */
		/* #(substitution) * overlapoffset / alignedLength = expected # of subs in the overlapoffset part */
		substitutionsInEdge += (double)(sourceRead->getNumOfSubstitutionsInRead() * e->getOverlapOffset())/(double)(sourceRead->getAlignedLength());
	}
	else{                                   /* composite edge, add up substitutions for each overlap offset */
		UINT64 lastOverlap = e->getOverlapOffset(); /* last overlap offset, not stored in the list */
		/* first overlap offset, on the source read */
		substitutionsInEdge += (double)(sourceRead->getNumOfSubstitutionsInRead() * e->getListOfOverlapOffsets()->at(0))/(double)(sourceRead->getAlignedLength());
		lastOverlap -= e->getListOfOverlapOffsets()->at(0);
		/* The following overlaps except the last one, since the offset length is not stored in the list */
		for(UINT64 i = 0; i < (e->getListOfReads()->size()-1); i++){
			Read *r = getReadFromID(e->getListOfReads()->at(i));
			UINT16 ovl = e->getListOfOverlapOffsets()->at(i+1);
			substitutionsInEdge += (double)(r->getNumOfSubstitutionsInRead() * ovl)/(double)r->getAlignedLength();
			lastOverlap -= ovl;
		}
		Read *r = getReadFromID(e->getListOfReads()->at(e->getListOfReads()->size()-1));
		substitutionsInEdge += (double)(r->getNumOfSubstitutionsInRead() * lastOverlap)/(double)r->getAlignedLength();
	}
	return substitutionsInEdge;
}

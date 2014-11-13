/*
 * Dataset.h
 *
 * Created on: Tue Nov 11 12:45:42 EST 2014
 * Author: JJ Chai
 */

#include "Dataset.h"

/**********************************************************************************************************************
  Function to compare two reads by their start mapping coordinate. Used for sorting.
 **********************************************************************************************************************/
bool compareReads (Read *read1, Read *read2)
{
	return read1->getStartCoord() < read2->getStartCoord();
}

/**********************************************************************************************************************
  Default constructor
 **********************************************************************************************************************/
Dataset::Dataset(void)
{
	// Initialize the variables.
	numberOfUniqueReads = 0;
	numberOfReads = 0;
	minimumOverlapLength = 0;
	maxError = 0;
	maxErrorRate = 0.0;
	shortestReadLength = 0XFFFFFFFFFFFFFFFF;
	longestReadLength = 0X0000000000000000;
}


/**********************************************************************************************************************
 * Another constructor, from BLASR generated sam file
 **********************************************************************************************************************/
Dataset::Dataset(const string & inputSamFile, UINT64 minOverlap, UINT32 maxError, float maxErrorRate)
{
	CLOCKSTART;
	// Initialize the variables.
	numberOfUniqueReads = 0;
	numberOfReads = 0;
	shortestReadLength = 0XFFFFFFFFFFFFFFFF;
	longestReadLength = 0X0000000000000000;
	reads = new vector<Read *>;
	minimumOverlapLength = minOverlap;
	maxError = 0;
	maxErrorRate = 0.0;
	UINT64 counter = 0;
	PacBioReadName = Utils::getFilename(inputSamFile);
	UINT64 goodReads = 0, badReads = 0;

	/** get reads from sam file **/
	ifstream samIn;
	samIn.open(inputSamFile.c_str());
	if(!samIn.is_open())
	{
		MYEXIT( "Input sam file " + inputSamFile + " cannot be opened." );
	}
	else
	{
		string line;
		while(samIn.good())
		{
			char c = samIn.peek();	// See if the line starts with @, if so, it's header line
			if (c == '@')
			{
				FILE_LOG(logDEBUG3) << "header line";
				samIn.ignore(maxCharInLine, '\n');
			}
			else	// not header line, alignment record
			{
				getline(samIn, line);
				Read *r = new Read(line);
				FILE_LOG(logDEBUG4) << "Scanned read: " << r->getReadName();
				if (testRead(r->getDnaStringForward()))
				{
					UINT32 len = r->getReadLength();
					if (len > longestReadLength)
						longestReadLength = len;
					if (len < shortestReadLength)
						shortestReadLength = len;
					reads->push_back(r);
					numberOfReads++;
					goodReads++;
				}
				else
					badReads++;
			}
		}
	}
	samIn.close();	// Close the sam file

	FILE_LOG(logINFO) << "Shortest read length: " << setw(5) << shortestReadLength;
	FILE_LOG(logINFO) << "Longest read length: " << setw(5) << longestReadLength;
	FILE_LOG(logINFO) << "Number of good reads: " << setw(5) << goodReads;
	FILE_LOG(logINFO) << "Number of bad reads: " << setw(5) << badReads;
	FILE_LOG(logINFO) << "Total number of reads: " << setw(5) << badReads + goodReads;
	sortReads();
	removeDupicateReads();	// Remove duplicated reads for the dataset.
	CLOCKSTOP;
}

/**********************************************************************************************************************
 * Another constructor, from BLASR generated sam file in piped stream instead of reading the file
 **********************************************************************************************************************/
Dataset::Dataset(stringstream * inputSamStream, UINT64 minOverlap, UINT32 maxError, float maxErrorRate)
{
	CLOCKSTART;
	// Initialize the variables.
	numberOfUniqueReads = 0;
	numberOfReads = 0;
	shortestReadLength = 0XFFFFFFFFFFFFFFFF;
	longestReadLength = 0X0000000000000000;
	reads = new vector<Read *>;
	minimumOverlapLength = minOverlap;
	maxError = 0;
	maxErrorRate = 0.0;
	UINT64 counter = 0;
	PacBioReadName = Utils::getFilename(inputSamFile);
	UINT64 goodReads = 0, badReads = 0;

	/** get reads from sam file **/
	//string line;
	//while(inputSamStream->good())
	//{
	//	char c = inputSamStream->peek();	// See if the line starts with @, if so, it's header line, ignore until EOL
	//	if (c == '@')
	//	{
	//		FILE_LOG(logDEBUG2) << "header line";
	//		inputSamStream->ignore(maxCharInLine, '\n');
	//	}
	//	else	// not header line, alignment record
	//	{
	//		stringstream record;
	//		inputSamStream >> record;
	//		/* this function has some problem... need to be fixed TODO */
	//		getline(inputSamStream, line);
	//		Read *r = new Read(line);
	//		FILE_LOG(logDEBUG2) << "Read in read: " << r->readName;
	//		if (testRead(r->getDnaStringForward))
	//		{
	//			UINT32 len = r->getReadLength();
	//			if (len > longestReadLength)
	//				longestReadLength = len;
	//			if (len < shortestReadLength)
	//				shortestReadLength = len;
	//			reads.push_back(r);
	//			numberOfReads++;
	//			goodReads++;
	//		}
	//		else
	//			badReads++;
	//	}
	//}

	//FILE_LOG(logINFO) << "Shortest read length: " << setw(5) << shortestReadLength;
	//FILE_LOG(logINFO) << "Longest read length: " << setw(5) << longestReadLength;
	//FILE_LOG(logINFO) << "Number of good reads: " << setw(5) << goodReads;
	//FILE_LOG(logINFO) << "Number of bad reads: " << setw(5) << badReads;
	//FILE_LOG(logINFO) << "Total number of reads: " << setw(5) << badReads + goodReads;
	//sortReads();
	//removeDupicateReads();	// Remove duplicated reads for the dataset.
	CLOCKSTOP;
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
	}
	delete reads;

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
	for(UINT64 i = 0; i < reads->size(); i++)	// Move the unique reads in the top of the sorted list. Store the frequencey of the duplicated reads.
	{
		if(reads->at(j)->getStartCoord()!= reads->at(i)->getStartCoord() \
				|| reads->at(j)->getDnaStringForward() != reads->at(i)->getDnaStringForward())	// Two reads have different start mapping coordinates, or different strings
		{
			j++;	// Number of unique reads increases by 1
			temp = reads->at(j);	// Save the read that was just checked
			reads->at(j) = reads->at(i);	// Switch reads.at(i) and reads.at(j)
			reads->at(i) = temp;
		}
		else if(i!=j)	// No more same reads as the current one, set the frequencies for this read
			reads->at(j)->setFrequency(reads->at(j)->getFrequency() + 1);
	}
	numberOfUniqueReads = j+1;
	FILE_LOG(logINFO) <<"Number of unique reads: " << numberOfUniqueReads;
	for(UINT64 i = 0 ; i < reads->size(); i++) 	// Assing ID's to the reads.
	{
		if(i < getNumberOfUniqueReads())
			reads->at(i)->setReadID(i);	// Read's ID, 0 based
		else
			delete reads->at(i); 	// Free the unused reads.
	}
	reads->resize(numberOfUniqueReads);	//Resize the vector of reads.
	CLOCKSTOP;
	return true;
}

/**********************************************************************************************************************
  Returns true if the read contains only {A,C,G,T} and does not contain more than 80% of the same nucleotide
 **********************************************************************************************************************/
bool Dataset::testRead(const seqan::DnaString & readDnaString)
{

	UINT64 cnt[4] = {0,0,0,0};
	size_t readLength = seqan::length(readDnaString);
	if ( readLength < minimumOverlapLength)
	{
		FILE_LOG(logDEBUG3) << "Read's length is smaller than the minimumOverlapLength " << readDnaString;
		return false;
	}
	seqan::Dna a= 'A', c = 'C', g = 'G', t = 'T';
	for(UINT64 i = 0; i < readLength; i++) // Count the number of A's, C's , G's and T's in the string.
	{
		if(readDnaString[i]!= a && readDnaString[i] != c && readDnaString[i] != g && readDnaString[i] != t)
		{
			FILE_LOG(logDEBUG3) << "Read has characters other than ACGT " << readDnaString;
			return false;
		}
		cnt[(seqan::ordValue(readDnaString[i]) >> 1) & 0X03]++; // Least significant 2nd and 3rd bits of ASCII value used here
	}
	UINT64 threshold = seqan::length(readDnaString)*.8;	// 80% of the length.
	if(cnt[0] >= threshold || cnt[1] >= threshold || cnt[2] >= threshold || cnt[3] >= threshold)
	{
		FILE_LOG(logDEBUG3) << "More than 80\% of the read string has the same character " << readDnaString;
		return false;	// If 80% bases are the same base.
	}
	return true;
}

/**********************************************************************************************************************
  Search a read in the dataset using binary search
 **********************************************************************************************************************/
Read * Dataset::getReadFromCoord(const INT32 & coord)	// Find a read in the database given the start mapping coord. Uses binary search in the list of reads.
{
	UINT64 min = 0, max = getNumberOfUniqueReads()-1;
	int comparator;
	while (max >= min) 	// At first search for the forward string.
	{
		UINT64 mid = (min + max) / 2; 	// Determine which subarray to search.
		comparator = reads->at(mid)->getStartCoord() - coord;
		if(comparator == 0)
			return reads->at(mid);
		else if (comparator < 0) 	// Change min index to search upper subarray.
			min = mid + 1;
		else if (comparator > 0) 	// Change max index to search lower subarray.
			max = mid - 1;
	}
	MYEXIT( "No reads were found in Dataset starting at coordinate: " + Utils::intToString(coord));
}

/**********************************************************************************************************************
  This function returns a read from its ID.
 **********************************************************************************************************************/
Read * Dataset::getReadFromID(UINT64 ID)
{
	if( ID > numberOfUniqueReads)	// ID outside the range.
	{
		stringstream ss;
		ss << "ID " << ID << " out of bound.";
		string s = ss.str();
		MYEXIT(s);
	}
	return reads->at(ID);
}

/**********************************************************************************************************************
  This function saves the starting coordinate sorted unique reads in a file. Used for debugging only.
 **********************************************************************************************************************/
void Dataset::saveReads(string fileName)
{
	ofstream outputFile;
	outputFile.open(fileName.c_str());
	if(!outputFile.is_open())
		MYEXIT("Unable to open file: " + fileName);
	for(UINT64 i = 0; i < numberOfUniqueReads; i++)
	{
		Read * read1 = getReadFromID(i);
		if(read1->superReadID!=0)
		{
			outputFile << ">" << i <<  " # " << read1->getReadName() \
				<< " # Start coord: " << read1->getStartCoord() \
				<< " # End coord: " << read1->getEndCoord() \
				<< " # Length: " << read1->getReadLength() \
				<< " # Contained in "  << read1->superReadID \
				<< "\n" << read1->getDnaStringForward() << endl;
		}
		else
		{
			outputFile << ">" << i <<  " # " << read1->getReadName() \
				<< " # Start coord: " << read1->getStartCoord() \
				<< " # End coord: " << read1->getEndCoord() \
				<< " # Length: " << read1->getReadLength() \
				<< " # Noncontained\n" << read1->getDnaStringForward() << endl;
		}
	}
	outputFile.close();
}

/**********************************************************************************************************************
 * This function prints the reads in tiling format, for checking the overlap alignment, and debugging.
 **********************************************************************************************************************/
void Dataset::printReadsTiling(string fileName)	// Print all the reads in tiling format. Used for checking the overlap (debugging)
{
	ofstream outputFile;
	outputFile.open(fileName.c_str());
	if(!outputFile.is_open())
		MYEXIT("Unable to open file: " + fileName);
	Read * read0 = reads->at(0);	// Read with leftmost mapping start coordinate
	int LeftMostCoord = read0->getStartCoord();
	if (LeftMostCoord < 0 )
	{
		string overHang(0-LeftMostCoord, '.');
		outputFile << overHang;
	}
	string PacBioSeq(1000,'*');
	outputFile << PacBioSeq << endl;
	for(UINT64 i = 0; i < numberOfUniqueReads; i++)
	{
		Read * read1 = getReadFromID(i);
		if (read1->getStartCoord() > LeftMostCoord)
		{
			string OffsetString(read1->getStartCoord() - LeftMostCoord, '.');
			outputFile << OffsetString;
		}
		outputFile << read1->getDnaStringForward() << endl;
	}
	outputFile.close();
}

/*
 * Dataset.h
 *
 * Created on: Tue Nov 11 12:45:42 EST 2014
 * Author: JJ Chai
 */


#ifndef DATASET_H_
#define DATASET_H_
#include "Common.h"
#include "Read.h"


/**********************************************************************************************************************
  Class to store the dataset
 **********************************************************************************************************************/
class Dataset
{
	private:
		UINT64 numberOfReads;	// Number of total reads present in the dataset.
		UINT64 numberOfUniqueReads;	// number of unique reads in the dataset.
		UINT64 numberOfNonContainedReads;	// number of unique reads in the dataset.
		UINT64 minimumOverlapLength;	// Length of the shortest read in the dataset.
		UINT32 maxError;	// maximum number of substitutions allowed in overlap between Illumina reads
		float maxErrorRate;	// maximum error rate
		string PacBioReadName;	// Name of the PacBio read

	//map< int, vector<Reads *> > *readMaps;	// Use map instead of vector to store all the reads
		vector<Read *> *reads;	// List of reads in the dataset.
		bool testRead(const string & readString);	// Test if the read is good. Contains only {A,C,G,T} and does not contain more than 80% of same base. 
	// The dataset contains only good quality reads.
		bool removeDupicateReads(void);	// Remove duplicate reads from the dataset. Frequency is stored for each of the reads.
		void sortReads(void);	// Sort the reads by their start mapping coordinates

	public:
		string inputSamFile;
		UINT64 shortestReadLength;
		UINT64 longestReadLength;

		Dataset(void);	// Default constructor.
		Dataset(const string & inputSamFile, UINT64 minOverlap, UINT32 maxError, float maxErrorRate);// another constructor, from a BLASR generated sam file
		Dataset(stringstream * inputSamStream, UINT64 minOverlap, UINT32 maxError, float maxErrorRate);// anotherconstructor, uses string stream directly instead of reading the file
		~Dataset(void);	// Default destructor.

		bool setPacBioReadName(const string & name){PacBioReadName = name; return true;}	// Set the name for the PacBio filtered subread
		UINT64 getNumberOfReads(void){return numberOfReads;}	// Get the number of total reads in the database.
		UINT64 getNumberOfUniqueReads(void){return numberOfUniqueReads;}	// Get the number of unique reads in the database.

		Read * getReadFromCoord(const INT32 & coord);	// Find a read in the database given the start mapping coord. Uses binary search in the list of reads.
		Read * getReadFromID(UINT64 ID);	// Find a read in the database given the ID in constant time.
		void saveReads(string fileName);	// Save all the sorted unique reads in a text file. Used for debugging.
		void printReadsTiling(string fileName);	// Print all the reads in tiling format. Used for checking the overlap (debugging)

};


#endif /* DATASET_H_ */

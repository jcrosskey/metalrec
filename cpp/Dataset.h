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
		UINT64 totalBps; /* Total number of aligned bases in the dataset, can be used to calculate the coverage depth */
		UINT64 numberOfReads;	// Number of total reads present in the dataset.
		UINT64 numberOfUniqueReads;	// number of unique reads in the dataset.
		UINT64 minimumOverlapLength;	// Length of the shortest read in the dataset.
		string PacBioReadName;	// Name of the PacBio read
		UINT16 PacBioReadLength; // Length of the PacBio read
		float subRate; /* substitution rate allowed in the Illumina read */
		float indelRate; /* indel rate allowed in the Illumina read */
		float insRate; /* insertion rate allowed */
		float delRate; /* deletion rate allowed. */
		float percentInLR; /* percentage of short read length contained in long read */

		vector<Read *> *reads;	// List of reads in the dataset.
		bool testRead(const string & readDnaString);	// Test if the read is good
		bool removeDupicateReads(void);	// Remove duplicate reads from the dataset. Frequency is stored for each of the reads.
		void sortReads(void);	// Sort the reads by their start mapping coordinates

	public:
		/*  Members */
		string inputSamFile;
		UINT64 shortestReadLength;
		UINT64 longestReadLength;
		UINT64 numberOfNonContainedReads;	// number of unique reads in the dataset. 0 in initialization stage, contained reads will be marked in overlapgraph building stage

		/* Constructors and destructor */
		Dataset(void);	// Default constructor.
		Dataset(const string & inputSamFile, UINT64 minOverlap, const float & indelRate, const float & subRate, const float & insRate, const float & delRate, const float & percent_inLR);// another constructor, from sam file with given file name
		Dataset(FILE * inputSamStream, UINT64 minOverlap, const float & indelRate, const float & subRate, const float & insRate, const float & delRate, const float & percent_inLR);// another constructor, uses input stream directly instead of reading the file
		Dataset(const Dataset & D);
		Dataset & operator= (const Dataset & D);
		~Dataset(void);	// Default destructor.

		/* mutators */
		bool AddDataset(FILE * inputSamStream);
		void setPacBioReadName(const string & name){PacBioReadName = name;}	// Set the name for the PacBio filtered subread
		void setIndelRate(const float & indel_rate){indelRate = indel_rate;}
		void setInsRate(const float & ins_rate){insRate = ins_rate; }
		void setDelRate(const float & del_rate){delRate = del_rate; }
		void setSubRate(const float & sub_rate){subRate = sub_rate; }
		void setPercentInLR(const float & percent_inLR){percentInLR = percent_inLR; }
		void setLRLength(const UINT16 & length){PacBioReadLength = length; }
		bool finalize(void);            /* Sort reads and remove duplicated reads */

		/* accessors */
		string getPacBioReadName(void){return PacBioReadName;}	// Get the name for the PacBio filtered subread
		UINT64 getPacBioReadLength(void){return PacBioReadLength;}	// Get the length of the PacBio filtered subread
		UINT64 getNumberOfReads(void){return numberOfReads;}	// Get the number of total reads in the database.
		UINT64 getNumberOfUniqueReads(void){return numberOfUniqueReads;}	// Get the number of unique reads in the database.
		UINT64 getNumberOfNonContainedReads(void);
		UINT64 getBpsOfNonContainedReads(void);
		UINT64 getTotalBps(void){return totalBps;}
		bool getCoveredInfo(INT64 & coveredBases, INT64 & coveredRegions);

		vector<Read *> getReadsFromCoord(const INT32 & coord);	// Find a read in the database given the start mapping coord. Uses binary search in the list of reads.
		Read * getReadFromID(UINT64 ID);	// Find a read in the database given the ID in constant time.
		void saveReads(string fileName, bool writeContained=false);	// Save all the sorted unique reads in a text file. Used for debugging.
		void printReadsTiling(string fileName);	// Print all the reads in tiling format. Used for checking the overlap (debugging)
		double getWeight(Edge * e);
		double getSubsOnEdge(Edge *e);
};


#endif /* DATASET_H_ */

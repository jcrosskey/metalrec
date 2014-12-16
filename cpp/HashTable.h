/*
 * HashTable.h
 *
 *  Created on: Apr 22, 2013
 *      Author: b72
 */

#ifndef HASHTABLE_H_
#define HASHTABLE_H_
#include "Common.h"
#include "Dataset.h"

/**********************************************************************************************************************
	Class to store hashtable.
**********************************************************************************************************************/
class HashTable{
	private:
		Dataset *dataSet;							// Pointer of the dataset.
		UINT64 hashTableSize; 						// Ted: Size of the hash table. This is the prime number of mod operation.
		vector < vector<UINT64> *> *hashTable; 		// Ted: List of hash keys of read number (ID) and orientation.
		UINT16 hashStringLength;					// Ted: Length of prefix and suffix of the reads to hash. This is equal to the minumum overlap length.
		UINT64 numberOfHashCollision;				// Counter to count the number of hash collisions. For debugging only.
		bool insertIntoTable(Read *read, string substring, UINT64 orientation);	// Insert a string in the hash table.
		bool hashRead(Read *read); 					// Ted: Hash prefix and suffix of the read and its reverse complement in the hash table. Turn over to the constant
		void setHashTableSize(UINT64 size); 		// Set the size of the hash table.
	public:
		HashTable(void);							// Default constructor.
		HashTable(Dataset *d);						// Another constructor.
		~HashTable();								// Destructor.
		bool insertDataset(Dataset *d, UINT64 minOverlapLength);	// Insert the dataset in the hash table.
		vector<UINT64> * getListOfReads(string subString); 			// Get the list of reads that contain subString as prefix or suffix.
		UINT64 hashFunction(string subString); 						// Hash function.
		UINT64 getHashTableSize(void){return hashTableSize;}		// Get the size of the hash table.
		UINT64 getHashStringLength() {return hashStringLength;}		// Get the hash string length.
		Dataset * getDataset(void){return dataSet;}					// Get the pointer to the dataset.
};


#endif /* HASHTABLE_H_ */

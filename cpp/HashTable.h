/*
 * HashTable.h
 * Hash table used to build directed graph for assembly of reads, while allow sustitution errors in the overlap instead of only perfect matches.
 *
 *  Created on: Tue Dec 16 10:22:58 EST 2014
 *      Author: JJ Chai
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
		Dataset *dataSet;	// Pointer to the dataset from which the hash table is built
		UINT64 hashTableSize; 	// Size of the hash table. This is the prime number of mod operation.
		vector < vector<UINT64> *> *hashTable; 	// List of hash keys of read number (ID) and orientation.
		UINT16 hashStringLength;	// Length of prefix and suffix of the reads to hash. This is equal to the minumum overlap length.
		UINT64 numberOfHashCollision;	// Counter to count the number of hash collisions. For debugging only.
		bool insertIntoTable(Read *read, string substring, UINT64 orientation);	// Insert a string in the hash table.
		bool hashRead(Read *read); 	// Hash prefix and suffix of the read and its reverse complement in the hash table. Turn over to the constant
		void setHashTableSize(UINT64 size); 	// Set the size of the hash table.
	public:
		HashTable(void);	// Default constructor.
		~HashTable();	// Destructor.
		bool insertDataset(Dataset *d, UINT64 hashLength);	// Insert the dataset in the hash table.
		vector<UINT64> * getListOfReads(string subString); 	// Get the list of reads that contain subString as prefix or suffix.
		UINT64 hashFunction(string subString); 	// Hash function.
		UINT64 getHashTableSize(void){return hashTableSize;}	// Get the size of the hash table.
		UINT64 getHashStringLength() {return hashStringLength;}	// Get the hash string length.
		Dataset * getDataset(void){return dataSet;}	// Get the pointer to the dataset.
};


#endif /* HASHTABLE_H_ */

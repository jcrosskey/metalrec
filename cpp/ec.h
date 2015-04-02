#ifndef EC_H
#define EC_H

#include <stdexcept> // for standard exceptions out_or_range
#include "Utils.h"
/*** utility functions ***/
#include "Common.h"
/*** Read class ***/
#include "Read.h"
/*** Dataset class ***/
#include "Dataset.h"
/*** Edge class ***/
#include "Edge.h"
/*** HashTable class ***/
#include "HashTable.h"
/*** OverlapGraph class ***/
#include "OverlapGraph.h"
/* check file existence and permissions */
#include <unistd.h>

/********************************************************
 * ec.h
 *
 * Created by JJ Chai  Thu Sep 18 09:35:28 EDT 2014. Last modified Fri Mar 20 EDT 2015
 * Copyright (c) 2014 JJ Chai (ORNL). Allrights reserved.

 ********************************************************/

// error correction job
void ec(const vector<string> & bamFiles, const string & PacBioName,  
		const UINT16 & PacBioLength, const string & allFileName,
		const string & samtools_path, const string & outDir,
		const UINT64 & minimumOverlapLength, const UINT64 & hashStringLength,
		const UINT32 & maxError, const UINT32 &rubberPos,
		const float & indelRate, const float & insRate, const float & delRate, const float & subRate, const float & maxErrorRate, const UINT16 minPacBioLength);

void ec_stream(const vector<string> & bamFiles, const string & PacBioName,  
		const UINT16 & PacBioLength, const string & allFileName,ofstream & outFastaStream,
		const string & samtools_path, const string & outDir,
		const UINT64 & minimumOverlapLength, const UINT64 & hashStringLength,
		const UINT32 & maxError, const UINT32 &rubberPos,
		const float & indelRate, const float & insRate, const float & delRate, const float & subRate, const float & maxErrorRate, const UINT16 minPacBioLength);
#endif

#ifndef DIRECTORYSTRUCTURE_HPP
#define DIRECTORYSTRUCTURE_HPP

/********************************************************************
 * directoryStructure.hpp
 
 * Class used to manage files in a given directory.
 * Reads a particular directory and gets all the *.chro files from that directory.
 * Has utility functions to provide file count and individual files.
 *
 * Created by JJ Chai  on 10/17/2013 based on Yingfeng's directorystructure.hpp.
 * Copyright (c) 2013 JJ Chai (ORNL). Allrights reserved.
********************************************************************/


#ifndef _WIN32
#include <unistd.h>
#endif


#include <cstring>
#include <iostream>
#include <vector>
#include <string>
#include <sys/types.h>
#include <dirent.h>

using namespace std;

class DirectoryStructure
{
	public:
		/** Constructors. Create an instance by providing a directory name. **/
    
		DirectoryStructure();
		DirectoryStructure( const string & );

        /** destructor **/
		virtual ~DirectoryStructure();

		/* 
		 * Set the directory explicitly, use it when you want to use
		 * the same instance for a different directory.
		 */
		void setDirectory( const string & );
		
		/*
		 * Set the tail pattern for the files to be detected.
		 * Eg: To detect *.chro files, set the pattern as
		 * ".chro".
		 *
		 * BEWARE: This class doesn't support full fledged pattern
		 * matching
		 */
		virtual void setPattern( string pattern );

		/*
		 * Get the number of matched files 
		 */
		size_t getFileCount();

		/*
		 * Get all the files at a time.
		 */
		void getFiles( vector<string> & fileList );

		/*
		 * Re-entrant function, gives you serial access to all the files.
		 */
		// start from the first file 
		void resetListIterator();
    
		// return the next filename,
		// return empty string "", if there is no more file
		string getNextFile();


	private:

		// the total number of files in this directory
        size_t fileCount;

		// the name of the directory
		string 	Directory;

		// the list of file names in this directory
		vector<string> fileList;

		// iterator used for serial access
		vector<string>::iterator fileListItr;


		// determine if this filename matches this pattern
		bool matchPattern(  const string & pattern, const string & filename );

};

#endif //DIRECTORYSTRUCTURE_HPP


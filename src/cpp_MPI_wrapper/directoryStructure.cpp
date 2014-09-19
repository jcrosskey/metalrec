/********************************************************************
 * directoryStructure.cpp
 
 * Class used to manage files in a given directory.
 * Reads a particular directory and gets all the *.chro files from that directory.
 * Has utility functions to provide file count and individual files.
 *
 * Created by JJ Chai  on 10/17/2013 based on Yingfeng's directorystructure.hpp.
 * Copyright (c) 2013 JJ Chai (ORNL). Allrights reserved.
 ********************************************************************/
#include "directoryStructure.hpp"

DirectoryStructure::DirectoryStructure()
{
	fileCount = 0;
    Directory = "";
}

DirectoryStructure::DirectoryStructure( const string & dir )
{
	fileCount = 0;
	Directory = dir;
}


DirectoryStructure::~DirectoryStructure()
{

}

void DirectoryStructure::setDirectory( const string & dir )
{
	Directory = dir;
}

void DirectoryStructure::setPattern( string pattern )
{

	fileList.clear();
	
	DIR *d;
	struct dirent *dir;
	d = opendir(Directory.c_str());
	if (d)
	{
	  while ((dir = readdir(d)) != NULL)
	  {
	    string filename( dir->d_name );
	    if( matchPattern( pattern, filename) )
		    fileList.push_back( filename );
	  }

	  closedir(d);
	}
	
	resetListIterator();
	fileCount = fileList.size();
}

size_t DirectoryStructure::getFileCount()
{
	return fileCount;
}

void DirectoryStructure::getFiles( vector<string> & files )
{

	for( unsigned int i = 0; i < fileList.size(); ++i)
		files.push_back( fileList[i] );
}

void DirectoryStructure::resetListIterator()
{
	fileListItr = fileList.begin();
}

string DirectoryStructure::getNextFile()
{
	if ( fileListItr != fileList.end() )
	{
		string nextfile = *fileListItr;
		fileListItr++;
		return nextfile;
	}
	else
	{
		return "";
	}
}

bool DirectoryStructure::matchPattern( const string & pattern, const string & filename )
{
	// if this filename is ".", ".." or "lost+found"
	// it is not considered
	if( ( filename.compare( "." ) == 0 )||
		( filename.compare( ".." ) == 0 ) ||
		( filename.compare("lost+found") == 0 ) ||
		( filename.length() <= pattern.size() ) )
	{
		return false;
	}
	
	// if the pattern is empty or "*", 
	// then all filenames are matched
	if( pattern == "" || pattern == "*" )
		return true;

	// if the filename ends with pattern
	// then it is matched
	if ( ( filename.compare( filename.length() - pattern.length(), pattern.length(), pattern ) == 0 ) )
		return true;
	else
		return false;

}


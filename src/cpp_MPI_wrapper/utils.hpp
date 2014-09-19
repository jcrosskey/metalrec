#ifndef UTILS_HPP
#define UTILS_HPP

/********************************************************************
 * utils.hpp
 *
 * utils.cpp includes useful common uility functions.
 * 
 * Created by JJ Chai  on 10/17/2013 based on Ted's utils.hpp.
 * Copyright (c) 2013 JJ Chai (ORNL). Allrights reserved.
********************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <sys/stat.h>
#include "directoryStructure.hpp"
#include <stdio.h>
// for getting current working directory
#ifdef WINDOWS
	#include <direct.h>
	#define GetCurrentDir _getcwd
#else
	#include <unistd.h>
	#define GetCurrentDir getcwd
#endif


/**
 * utils: utility functions for file (directory) path, filename,
 *        base filename without extension, and format changing.
 */
namespace Utils
{
    /** get path separator **/
    std::string getPathSeparator();

    /** check file exist **/
    bool isFileExist(const std::string & filename);

    /** if file exists, then remove the file **/
    void ifFileExistRemove(const std::string & filepath);
    
    /** create directory if it did not exist **/
    void mkdirIfNonExist(const std::string & dirname);
    
    /** check if file is empty **/
    bool isFileEmpty(const std::string & filename);
    
    /** if file is empty, remove the file **/
    void ifFileEmptyRemove(const std::string & filename);

    /** get program name only from program path **/
    std::string getProgramName(const std::string & program_path);

    /** get program directory from program path **/
    std::string getProgramDir(const std::string & program_path);

    /** if error, then print error and exit program **/
    void exitWithError(const std::string & error);

    /** method to convert string to int **/
    int stringToInt(const std::string & input_string);

    /** method to convert string to unsigned int **/
    unsigned int stringToUnsignedInt(const std::string & input_string);

     /** method to convert string to float **/
    float stringToFloat( const std::string& input_string );

    /** method to convert string to double **/
    double stringToDouble(const std::string & input_string);

    /** method to convert int to string **/
    std::string intToString(const int & input_value);

    /** method to convert unsigned int to string **/
    std::string unsignedIntToString(const unsigned int & input_value);

    /** method to convert double to string **/
    std::string doubleToString(const double & input_value);

    /** check directory **/
    bool isDirectory(const std::string & directory);

    /** check fasta format **/
    bool isFastaFormat(const std::string & input_string);

    /** check fastq format **/
    bool isFastqFormat(const std::string & input_string);

    /** get file name only from file path **/
    std::string getFilename(const std::string & filepath);
    
    /** get filebase from filename without extension and path **/
    std::string getFilename2(const std::string & filename);


    /** get filebase from filename without extension **/
    std::string getFilebase(const std::string & filename);

    /** get filebase from filename without extension **/
    std::string getFilebase2(const std::string & filename);

    /** Get current date/time, format is YYYY-MM-DD.HH:mm:ss **/
    std::string currentDateTime();

    /** string vector to comma separated string **/
    std::string vectorToCommaString(const std::vector<std::string> & input_vector);
    
    /** comma (or other delimiter) separated string to string vector **/
    std::vector<std::string> commaStringToVector(std::string & comma_string, const char & delimit);
    
    /** determine if this filename matches this pattern **/
    bool matchPattern(  const std::string & pattern, const std::string & filename );
    
    /** find all files with specified extensions in a directory **/
    void getFiles(const std::string & directory, const std::string & extension,
                  std::vector<std::string> & filenames);
    
    /** find the last non-empty line of a file **/
    std::string last_line(const std::string & filename);

    /** get current working directory, where the program is invoked **/
    std::string get_cwd();
};

#endif //UTILS.HPP

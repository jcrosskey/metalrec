/********************************************************************
 * utils.cpp
 *
 * utils.cpp includes useful common uility functions.
 * 
 * Created by JJ Chai on 10/17/2013 based on Ted's script.
 * Copyright (c) 2013 JJ Chai (ORNL). Allrights reserved.
********************************************************************/

#include "utils.hpp"

using namespace std;

namespace Utils
{

    /** get path separator **/
    std::string getPathSeparator()
    {
    #if _WIN32
        return "\\" ;
    #else
        return "/" ;
    #endif
    }

    /** check file exist **/
    bool isFileExist(const std::string & filename)
    {
        // infile
        std::ifstream infile(filename.c_str());

        // return true if none of the stream's error state flags (eofbit, failbit and badbit) is set
        return infile.good();
    }

   
    /** if file exists, then remove the file **/
    void ifFileExistRemove(const std::string & filename)
    {
        // if file exists,
        if (isFileExist(filename))
        {
            // then remove if there is no error
            if (remove(filename.c_str()) != 0)
                exitWithError("*** Error: Failed to remove: " + filename);
        }
    }

    /** create directory if it did not exist **/
    void mkdirIfNonExist(const std::string & dirname)
    {
        if (!Utils::isDirectory(dirname)) {
            string cmd = "mkdir -p "+ dirname;
            const int res = system(cmd.c_str());
            if (res != 0 )
                Utils::exitWithError("*** Error: Failed command: " + cmd);
            
        }
    }
    /** check if file is empty **/
    bool isFileEmpty(const std::string & filename)
    {

        
        //check if file exists first
        if (!isFileExist(filename)) {
            //std::cout << "*** file " << filename << "does not exist .... " << std::endl ;
            return true;
        }
        else
        {
            //infile
            std::ifstream infile(filename.c_str());
            
            //check file size
            infile.seekg(0,std::ios::end);
            size_t size = infile.tellg();
            infile.close();
            
            if( size == 0 )
                return true;
            else
                return false;
        }
        
    }
    
    /** if file is empty, remove the file **/
    void ifFileEmptyRemove(const std::string & filename)
    {
        // if file is empty
        if (isFileEmpty(filename))
        {
            // then remove if there is no error
            if (remove(filename.c_str()) != 0)
                exitWithError("*** Error: Failed to remove: " + filename);
        }
    }
    
    
    
    /** get program name only from program path **/
    std::string getProgramName(const std::string & program_path)
    {
        std::string program_name;

        // position
        size_t program_position = program_path.find_last_of(getPathSeparator());

        // find path separator
        if (program_position != std::string::npos)
            program_name.assign(program_path.begin()+program_position+1,program_path.end());
        
        // couldn't find path separator
        else
            program_name = program_path;

        // return 
        return program_name;
    }

    /** get program directory from program path **/
    std::string getProgramDir(const std::string & program_path)
    {
        std::string program_dir;

        // position
        size_t program_position = program_path.find_last_of(getPathSeparator());

        // find path separator
        if (program_position != std::string::npos)
            program_dir.assign(program_path.begin(),program_path.begin()+program_position+1);
        // couldn't find path separator
        else
            program_dir = "." + getPathSeparator();

        // return 
        return program_dir;
    }

    /** exit program when error happens **/
    void exitWithError(const std::string &error) 
    {
        std::cerr << error;
        exit(EXIT_FAILURE);
    }

    /** method to convert string to unsigned int **/
    int stringToInt( const std::string &input_string ) 
    {
        std::istringstream input_string_stream(input_string);
        int x;
        if (!(input_string_stream >> x)) {
            std::cerr << input_string;
            exitWithError("*** Error: string was not converted to int correctly!\n");
        }   

        // return 
        return x;
    } 

    /** method to convert string to unsigned int **/
    unsigned int stringToUnsignedInt( const std::string &input_string ) 
    {
        std::istringstream input_string_stream(input_string);
        unsigned int x;
        if (!(input_string_stream >> x)) {
            exitWithError("*** Error: string was not converted to unsigned int correctly!\n");
        }   

        // return 
        return x;
    } 

    /** method to convert string to float **/
    float stringToFloat( const std::string& input_string ) 
    {
        std::istringstream input_string_stream(input_string);
        float x;
        if (!(input_string_stream >> x)) {
            std::cerr << "input string is " << input_string << "\n";
            exitWithError("*** Error: string was not converted to float correctly!");
        }   

        // return 
        return x;
    } 
    /** method to convert string to double **/
    double stringToDouble( const std::string& input_string ) 
    {
        std::istringstream input_string_stream(input_string);
        double x;
        if (!(input_string_stream >> x)) {
            std::cerr << "input string is " << input_string << "\n";
            exitWithError("*** Error: string was not converted to double correctly!");
        }   

        // return 
        return x;
    } 

    /** method to convert int to string **/
    std::string intToString( const int & input_value )
    {
        std::string input_value_str;
        std::ostringstream convert;
        if (!(convert << input_value)) {
            exitWithError("*** Error: unsigned int was not converted to string correctly!");
        }   
        input_value_str = convert.str();

        // return 
        return input_value_str;
    }

    /** method to convert unsigned int to string **/
    std::string unsignedIntToString( const unsigned int & input_value )
    {
        std::string input_value_str;
        std::ostringstream convert;
        if (!(convert << input_value)) {
            exitWithError("*** Error: unsigned int was not converted to string correctly!");
        }   
        input_value_str = convert.str();

        // return 
        return input_value_str;
    }

    /** method to convert double to string **/
    std::string doubleToString( const double & input_value )
    {
        std::string input_value_str;
        std::ostringstream convert;
        if (!(convert << input_value)) {
            exitWithError("*** Error: unsigned int was not converted to string correctly!");
        }   
        input_value_str = convert.str();

        // return 
        return input_value_str;
    }

    /** check directory **/
    bool isDirectory( const std::string & directory )
    {
        struct stat sb;
        if (stat(directory.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode) )
            return true;
        else
            return false;
    }

    /** check fasta format **/
    bool isFastaFormat(const std::string & input_string)
    {

        // list of fasta format
        std::vector<std::string> file_format;
        file_format.push_back("fasta");
        file_format.push_back("fa");
        file_format.push_back("fas");
        file_format.push_back("fna");
        file_format.push_back("ffn");

        // get file extension
        std::string file_ext = input_string.substr(input_string.find_last_of(".") + 1);

        // transform to lower characters
        std::transform(file_ext.begin(), file_ext.end(), file_ext.begin(), ::tolower);

        // loop list
        for(unsigned int i=0; i < file_format.size(); i++) {
            if (file_ext == file_format[i]) 
                return true;
        }
        // return 
        return false;
    }

    /** check fastq format **/
    bool isFastqFormat(const std::string & input_string)
    {

        // list of fasta format
        std::vector<std::string> file_format;
        file_format.push_back("fastq");
        file_format.push_back("fq");
        file_format.push_back("faq");

        // get file extension
        std::string file_ext = input_string.substr(input_string.find_last_of(".") + 1);
        std::transform(file_ext.begin(), file_ext.end(), file_ext.begin(), ::tolower);

        // loop list
        for(unsigned int i=0; i < file_format.size(); i++) {
            if (file_ext == file_format[i]) 
                // return 
                return true;
        }
        // return 
        return false;
    }

    /** get file name only from file path **/
    std::string getFilename(const std::string & filepath)
    {
        std::string filename;

        size_t pos = filepath.find_last_of(getPathSeparator());
        if (pos != std::string::npos)
            filename.assign(filepath.begin()+pos+1,filepath.end());
        else
            filename = filepath;

        // return 
        return filename;
    }
    
    /** get filebase from filename without extension and path **/
    std::string getFilename2(const std::string & filename)
    {
        std::string filebase;
        filebase = getFilename(filename);
        filebase = getFilebase(filebase);
        return filebase;
    }

    /** get filebase from filename without extension **/
    std::string getFilebase(const std::string & filename)
    {
        std::string filebase;

        size_t pos = filename.find_last_of(".");
        if (pos != std::string::npos)
            filebase.assign(filename.begin(),filename.begin()+pos);
        else
            filebase = filename;

        // return 
        return filebase;
    }

    /** get filebase from filename i.e., AAA.CC.X -> AAA **/
    std::string getFilebase2(const std::string & filename)
    {
        std::string filebase;

        // parse base filename
        std::istringstream filename_stream(filename);
        std::vector<std::string> fields_vector;
        for (std::string field; getline(filename_stream, field, '.'); ) { 
           fields_vector.push_back(field);
        }   
        if (fields_vector.size() > 2) {
            filebase = ""; 
            for(unsigned int i = 0; i < fields_vector.size() - 2; i++) {
                filebase += fields_vector[i] + '.';
            }   
            filebase = filebase.substr(0, filebase.size()-1);
        }   
        else {
            filebase = filename;
        }

        // return 
        return filebase;
    }

    /** Get current date/time, format is YYYY-MM-DD.HH:mm:ss **/
    std::string currentDateTime() 
    {
        time_t     now = time(0);
        struct tm  tstruct;
        char       buf[80];
        tstruct = *localtime(&now);
        strftime(buf, sizeof(buf), "[%Y-%m-%d.%X]", &tstruct);

        // return 
        return buf;
    }

    /** string vector to comma separated string **/
    std::string vectorToCommaString(const std::vector<std::string> & input_vector)
    {
        // variables
        std::stringstream ss;
        std::string comma_string;

        // loop list
        for (unsigned int i=0; i<input_vector.size(); i++ ) {
            if (i != 0)
                ss << ",";
            ss << input_vector[i];
        }
  
        // comma_string
        comma_string = ss.str();

        // return 
        return comma_string;
    }
    
    /** comma separated string to string vector **/
    std::vector<std::string> commaStringToVector(std::string & comma_string, const char & delimit)
    {
        std::vector<std::string> strVec;
        std::size_t pos = 0;
        std::size_t comma_pos = comma_string.find(delimit,pos);
        string str;
        while (comma_pos != std::string::npos) {
            str = comma_string.substr(pos, comma_pos-pos);
            strVec.push_back(str);
            pos = comma_pos + 1;
            comma_pos = comma_string.find(delimit,pos);
        }
        str.assign(comma_string.begin()+pos,comma_string.end());
        strVec.push_back(str);
        return strVec;
    }
    
    /** determine if this filename matches this pattern **/
    bool matchPattern(  const std::string & pattern, const std::string & filename )
    {
        // if filename is ".", "..", or "lost+found", return false
        if ( (filename.compare(".") == 0) ||
             (filename.compare("..") == 0) ||
             (filename.compare("lost+found") == 0 ) ||
             (filename.length() <= pattern.size())
            ){
            return false;
        }
        
        // if pattern is empty or "*", match all files
        if ((pattern.compare("") == 0) || (pattern.compare("*") == 0)) {
            return true;
        }
        // get file extension
        std::string file_ext = filename.substr(filename.find_last_of(".") + 1);
        
        // compare extension with the pattern
        if (file_ext.compare(pattern) == 0)
            return true;
        else
            return false;
    }
    
    /** find all files with specified extensions in a directory **/
    void getFiles(const std::string & directory, const std::string & extension,
                  std::vector<std::string> & filenames)
    {
        DirectoryStructure searchDir(directory);
        
        searchDir.setPattern(extension);
        searchDir.getFiles(filenames);
        
        int i, fileNum;
        fileNum = (int) searchDir.getFileCount();
        vector<string>::reverse_iterator rit = filenames.rbegin();
        // full path of the input files and the output files
        for (i = 0; i < fileNum; rit++,i++) {
            //cout << *rit << " -> ";
            *rit = directory + getPathSeparator() + *rit;
            //cout << *rit << endl;
        }

    }
    
    /** find the last non-empty line of a file **/
    std::string last_line(const std::string & filename)
    {
        std::ifstream fin(filename.c_str());
        std::string lastline;
        
        if(fin)
        {
            char ch;
            
            // go to one spot before the EOF
            fin.seekg(-1,std::ios_base::end);
            
            fin.get(ch);
            
            // skip the last EOF
            if (ch =='\n')
            {
                //std::cout << "file ended with EOF \n";
                fin.seekg(-2,std::ios_base::cur);
            }
            else
            {
                fin.unget();
                fin.seekg(-2,std::ios_base::cur);
            }
            
            // indicator of end of looping
            bool keepLooping = true;
            
            while (keepLooping) {
                
                // get current byte's data
                fin.get(ch);

                // if the data was at or before the beginning
                // first line is the last line
                // stop here
                if ((int)fin.tellg() <= 1) {
                    //std::cout << "fin.tellg = " << (int)fin.tellg() << std::endl;
                    fin.seekg(0);
                    keepLooping = false;
                    //std::cout << "now at beginning \n";
                }
                // if the data is a new line, stop
                else if (ch == '\n')
                {
                    keepLooping = false;
                    //std::cout << "now at a new line \n";
                }
                
                // keep moving until the first 2 cases
                else
                {
                    fin.seekg(-2,std::ios_base::cur);
                }
            }
            
            //std::cout << "looping ended" << std::endl;
            
            getline(fin,lastline);
            
            fin.close();
        }
        
        
        return lastline;
    }

    std::string get_cwd()
    {
	char path[FILENAME_MAX];
	return ( GetCurrentDir(path, sizeof(path)) ? std::string(path) : std::string("") );
    }

}



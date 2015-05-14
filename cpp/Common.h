/*
 * Common.h
 *
 * Created on: April 22, 2013
 * Author: Md. Bahlul Haider
 */


#ifndef COMMON_H_
#define COMMON_H_

#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <time.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <math.h>
#include <streambuf>
#include <sys/stat.h>
#include <map>
#include <ctype.h>
#include "Utils.h"

#ifdef LOG1
#include "logcpp/log1.h"
#else
#include "logcpp/log.h"
#endif
using namespace std;

typedef unsigned char UINT8;
typedef unsigned short UINT16;
typedef short INT16;
typedef unsigned long UINT32;
typedef long INT32;
typedef unsigned long long UINT64;
typedef long long INT64;


extern int loglevel;
#define maxCharInLine 10000
// Exit code that displays the place of exit and message.
#define MYEXIT(a) { FILE_LOG(logERROR) << "Exit from File: " << __FILE__ << " Line: " << __LINE__ << " Function: " << __FUNCTION__ << "()" << endl << "Message: " << a; exit(0);}
// Print which function is currently executing. Only for functions that take long time


// To keep time information of functions.
#define CLOCKSTART clock_t begin = clock(); FILE_LOG(logDEBUG2) <<"Currently in file: " << __FILE__ << " Function: "<< __FUNCTION__ << "()";
#define CLOCKSTOP clock_t end = clock(); FILE_LOG(logDEBUG2)  << "Function " << __FUNCTION__ << "() finished in " << double(end - begin) / CLOCKS_PER_SEC<< " Seconds." << endl;

// To change the log level

#endif /* COMMON_H_ */

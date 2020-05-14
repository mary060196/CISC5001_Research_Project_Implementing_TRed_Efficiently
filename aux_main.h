//   |=================================================================|
//   |   TRED: a tool for detecting Tandem Repeats within sequences,   |
//   |               using the Edit Distance metric.                   |
//   |=================================================================|

// Copyright © 2007-2009 Dina Sokol, Justin Tojeira
// Distributed under the Aladdin Free Public License

// THIS SOFTWARE SHOULD BE ACCOMPANIED BY readme.txt AND license.html
// WE STRONGLY ENCOURAGE YOU TO READ BOTH BEFORE PROCEEDING
//------------------------------------------------------------------------------

#ifndef AUX_MAIN_H
#define AUX_MAIN_H

#include <stdio.h>       // For FILE type
#include "parameters.h"  // For MAX_PERIOD
#include <string>        // For string
// Library mingw-std-threads by Mega, downloaded at:
//   https://github.com/meganz/mingw-std-threads
// and added on 04.23.2020 by Miriam Briskman
// The reason is that, on win32 threading, the Thread C++11 library is
//   not fully supported, which prevents the program from using the
//   #include <thread> preprocessor directive.
#include "mingw.thread.h"

using namespace std;

/*          'thread_info' struct
A structure containing information passed to
functions that in which threads work.     */
// Miriam Briskman, 04.24.2020
struct thread_info {
    char *beg_str;
    string filename;
    long long times;
    unsigned long long offset;
};

void read_file (char*, unsigned long long&, FILE*); // Miriam Briskman, 02.23.2020
void threads_func (thread_info *); // Miriam Briskman, 04.24.2020
void threads_func_last (thread_info *); // Miriam Briskman, 04.24.2020

// Miriam Briskman, 03.02.2020
const long DOUBLE_PERIOD = 2*MAX_PERIOD,   // Twice the maximum period (see parameters.h for definition)
           QUAD_PERIOD = 4*MAX_PERIOD,     // Four times the maximum period (see parameters.h for definition)
           D_ERRORS_PLUS_1 = 2*MAX_ERRORS + 1; // Twice the maximum allowed errors + 1

// Number of threads used in the program (besides the main thread):
const unsigned short THREAD_NUM = thread::hardware_concurrency();

// External variables
extern int *numarray;
extern unsigned long long wholestr_len;

#endif

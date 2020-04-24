//   |=================================================================|
//   |   TRED: a tool for detecting Tandem Repeats within sequences,   |
//   |               using the Edit Distance metric.                   |
//   |=================================================================|

// Copyright © 2007-2009 Dina Sokol, Justin Tojeira
// Distributed under the Aladdin Free Public License

// THIS SOFTWARE SHOULD BE ACCOMPANIED BY readme.txt AND license.html
// WE STRONGLY ENCOURAGE YOU TO READ BOTH BEFORE PROCEEDING
//------------------------------------------------------------------------------



// This program will take an input string and process it completely using the
// framework of the Main-Lorentz algorithm. This file contains the function
// MAIN, which together with the function ONEITERATION constitutes the full
// implementation of the Main-Lorentz algorithm, adapted to work with edit
// distance matrices. In the MAIN function, the majority of the space that
// will be dynamically allocated is declared. This will be the "workspace" for
// our program, and most of the data we need throughout the program will be
// stored in this space.

// The MAIN function will handle the recursive portion of the Main-Lorentz
// algorithm, but it will do this iteratively. Its primary purpose is to pass
// the appropriate string to the function ONEITERATION.

// As per the Main-Lorentz algorithm, the MAIN function will pass the entire
// string, and each half of the entire string, and each half of each of those
// halves, and so on.

// ONEITERATION will handle the remainder of the Main-Lorentz algorithm. It will
// find the halfway point of the string it receives from MAIN, and compare the
// substrings originating at that point to those originating from each other
// point on the string. It will do this from the center out, alternating left
// then right. The variable "j" in the loop that comprises most of the function
// represents how far from the center of the string the points that we are
// comparing are. For each j, this loop is run twice, once on the left and once
// on the right.

#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include "parameters.h"
#include "aux_main.h"
#include <iostream>
#include <string>
#include <thread>
// Library mingw-std-threads by Mega, downloaded at:
//   https://github.com/meganz/mingw-std-threads
// and added on 04.23.2020 by Miriam Briskman
// The reason is that, on win32 threading, the Thread C++11 library is
//   not fully supported, which prevents the program from using the
//   #include <thread> preprocessor directive.
#include "mingw.thread.h"

using namespace std;

// Globals:
int *init_increments;      // External variable declared in 'buildk.h'
int *numarray;             // External variable declared in 'aux_main.h'
// File's sile minus comments and whitespaces:
unsigned long long wholestr_len = 0;  // External variable declared in 'aux_main.h'

int main (int argc, char *argv[])
{
   // Capture the beginning time of the program:
   auto start_time = chrono::high_resolution_clock::now();

   thread threads[THREAD_NUM]; // Threads array
   thread_info info[THREAD_NUM]; // Create an array of thread_info structs

   FILE *infile = NULL, *outfile = NULL, *tempfile = NULL;
   string filename[2];
   char *wholestr = NULL, byte;
   int i, j, levels = 0;
   unsigned long long offset = START_POS, // Miriam Briskman, 02.23.2020:
                      fileSize;           // Variable to hold the size of a file in terms of bytes.                      
   long long to_be_read,         // The number of characters remaining until end of wholestr.
             res; // Number of iterations of loop each thread has to complete.

   /*********** FOR INPUT/OUTPUT FILES AS COMMAND LINE ARGUMENTS ************/

   if (argc != 3)
   {
       cout << "Please manually enter the following 2 filenames below:" << endl << endl
            << "1) Enter the name of the sequence file : >> ";
       cin  >> filename[0];
       cout << endl << "2) Enter the name of the intermediary file : >> ";
       cin  >> filename[1];
       cout << endl << endl << "Thank you! Information is being processed..." << endl << endl;

       if ((infile = fopen(filename[0].c_str(), "a+")) == NULL){ // Changed to "append + read" mode
           cout << "Error opening sequence file.\n";
           exit (EXIT_FAILURE); 
       }
       if ((outfile = fopen(filename[1].c_str(), "wb")) == NULL){ 
           cout << "Error opening intermediary file.\n";
           exit (EXIT_FAILURE);
       }
   }
   else
   {
       if ((infile = fopen(argv[1], "a+")) == NULL){ // Changed to "append + read" mode
           cout << "Error opening sequence file.\n";
           exit (EXIT_FAILURE); 
       }
       if ((outfile = fopen(argv[2], "w")) == NULL){ 
           cout << "Error opening intermediary file.\n";
           exit (EXIT_FAILURE);
       }
       filename[0] = argv[1];
       filename[1] = argv[2];
   }   

   /********************** READING THE SEQUENCE ***************************/

   // Miriam Briskman, 02.23.2020:
   // Put the position of reading at the end of the file:
   fseek (infile, 0L, SEEK_END);

   // Thereby finding the size of the file in bytes:
   fileSize = ftell(infile);

   // Place a newline character at the end of the file (to ease reading):
   fseek (infile, -1, SEEK_CUR);
   if (getc(infile) != '\n')
   {
      fseek (infile, 1, SEEK_CUR);
      fputc ('\n', infile);
   }

   // Re-open the file in reading-only mode:
   if (argc != 3)
   {
       if ((infile = freopen(filename[0].c_str(), "r", infile)) == NULL){
       cout << "Error opening sequence file.\n";
       exit (EXIT_FAILURE);
       }
   }
   else if ((infile = freopen(argv[1], "r", infile)) == NULL){
       cout << "Error opening sequence file.\n";
       exit (EXIT_FAILURE); 
   }

   // Allocate only the memory needed for the file:
   wholestr = (char*) malloc((fileSize + 1)*sizeof(char));
   if (wholestr == NULL){
       cout << "Not enough memory in main.\n";
       exit (EXIT_FAILURE);
   }

   // Read-in the entire file (without comments or whitespaces) into the memory:
   read_file (wholestr, wholestr_len, infile);

   /************************* MEMORY ALLOCATION ******************************/

   // Miriam Briskman, 03.03.2020
   // Allocating memory for the initialization increment array [1][2][3]...[n]
   //   used to build the kxk matrices. The pointer to the array will be passed
   //   to 'oneiteration' and then to 'buildmatrixforward' and 'buildmatrixbackward'.
   // Also: (Copied by Miriam Briskman, on 03.03.2020, from 'buildk.cpp')
   // lower holds the information of the row currently being calculated, and
   // upper holds the information of the row above that, on which the
   // calculations are based. temp is used to swap pointers so that after a row
   // is finished upper points to it and lower is free for the next calculations
   init_increments = (int*) malloc ((D_ERRORS_PLUS_1) * sizeof(int));

   if (init_increments == NULL)
   {
       cout << "Memory allocation error in main.\n";
       exit (EXIT_FAILURE);
   }

   // Initialize the increments array:
   for (i = MAX_ERRORS; i < D_ERRORS_PLUS_1; i++)
     init_increments[i] = i - MAX_ERRORS;

   // Allocating memory for 'numarray': (Copied from the old 'createarray.cpp' by Miriam Briskman, 02.23.2020)
   j = (QUAD_PERIOD - 1)/MIN_LENGTH + 1;
   for (i = 1; i < j; i *= 2)
       levels++;

   numarray = (int*) malloc ((i*2 - 1)*sizeof(int)); // Size sufficient for the entire program
   if (numarray == NULL)
   {
       cout << "Memory allocation error in main.\n";
       exit (EXIT_FAILURE);
   }

   numarray[0] = DOUBLE_PERIOD;
   i = i/2 - 1;
   levels--;

   for (j = 0; j < i; j++) // Copied from old 'createarray.cpp' by Miriam Briskman, 02.23.2020
   {
       numarray[j*2 + 1] = numarray[j]/2;
       numarray[j*2 + 2] = numarray[j] - numarray[j*2 + 1];
   }

   /************************** CREATING THREADS *******************************/

   // Miriam Briskman, 04.24.2020

   // How many portions of DOUBLE_PERIOD will each thread analyze?
   res = (wholestr_len / DOUBLE_PERIOD) / THREAD_NUM;

   for (i = 0; i < THREAD_NUM - 1; i++)
   {
      // Setting relevant thread_info:
      (info[i]).beg_str = wholestr + res*i*DOUBLE_PERIOD;
      (info[i]).filename = filename[1] + to_string(i + 1) + ".data";
      (info[i]).times = res;
      (info[i]).offset = offset + res*i*DOUBLE_PERIOD;
  
      threads[i] = thread (threads_func, &info[i]); // Call to 'thread' constructor
   }

   // Setting relevant thread_info for last thread:
   (info[i]).beg_str = wholestr + res*(THREAD_NUM - 1)*DOUBLE_PERIOD;
   (info[i]).filename = filename[1] + to_string(THREAD_NUM) + ".data";
   (info[i]).times = res;
   (info[i]).offset = offset + res*(THREAD_NUM - 1)*DOUBLE_PERIOD;
  
   threads[i] = thread (threads_func_last, &info[i]); // Call to 'thread' constructor

   /************************** JOINING THREADS *******************************/

   for (i = 0; i < THREAD_NUM; i++)
      threads[i].join();

   /*************************** APPENDING ALL FILES **************************/

   for (i = 0; i < THREAD_NUM; i++)
   {
      if ((tempfile = fopen(((info[i]).filename).c_str(), "r")) == NULL)
      {
          cout << "Error opening outfile files.\n";
          exit (EXIT_FAILURE); 
      }

      if (!feof(tempfile))
          byte = fgetc (tempfile);
      while (!feof(tempfile)) 
      {
          fputc (byte, outfile);
          byte = fgetc (tempfile);
      }
      fclose(tempfile);
      remove (((info[i]).filename).c_str()); // Remove temporarily-created files from directory
   }

   /************************ MEMORY DEALLOCATION *****************************/

   free (init_increments);
   free (numarray); 
   free (wholestr);

   /*************************** CLOSING FILES ********************************/

   fclose (infile);
   fclose (outfile);

   // Capture the ending time of the program:
   auto end_time = std::chrono::high_resolution_clock::now();
   auto mtime = end_time - start_time;

   cout << "Total execution time: " << mtime/std::chrono::milliseconds(1) << " milliseconds." << endl;

   return EXIT_SUCCESS; // END MAIN
}

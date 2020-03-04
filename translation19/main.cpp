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
#include "buildk.h"
#include "oneiteration.h"
#include "parameters.h"
#include "tracek.h"
#include <iostream>

// Function Prototypes:
void read_file (char*, unsigned long long&, FILE*); // Miriam Briskman, 02.23.2020

using namespace std;

// Miriam Briskman, 03.02.2020
const long DOUBLE_PERIOD = 2*MAX_PERIOD;   // Twice the maximum period (see parameters.h for definition)
const long QUAD_PERIOD = 4*MAX_PERIOD;     // Four times the maximum period (see parameters.h for definition)
const long D_ERRORS_PLUS_1 
           = 2*MAX_ERRORS + sizeof(char); // Twice the maximum allowed errors + sizeof(char)

int *upper, *lower, *init_increments;      // External variables declared in 'buildk.h'
char *forward_pattern1, *forward_pattern2; // External variables declared in 'tracek.h'

int main (int argc, char *argv[])
{
   // Capture the beginning time of the program:
   auto start_time = chrono::high_resolution_clock::now();

   FILE *infile = NULL, *outfile = NULL;
   char filename[2][40];
   char *partialstr = NULL, *wholestr = NULL, *strpnt = NULL;
   int i, j, k, levels = 0,
       *numarray = NULL,
       **matrix_forward = NULL, **matrix_reverse = NULL,
       **left_matrix_forward = NULL, **left_matrix_reverse = NULL,
       **right_matrix_forward = NULL, **right_matrix_reverse = NULL;
   int num, temp;
   unsigned long long offset = START_POS, // Miriam Briskman, 02.23.2020:
                      fileSize,           // Variable to hold the size of a file in terms of bytes.
                      wholestr_len = 0,   // The size of file minus comments and whitespaces.
                      to_be_read;         // The number of characters remaining until end of wholestr.
   Entry *down_errors = NULL, *right_errors = NULL, *left_down_errors = NULL,
         *left_right_errors = NULL, *right_down_errors = NULL, *right_right_errors = NULL;

   /*********** FOR INPUT/OUTPUT FILES AS COMMAND LINE ARGUMENTS ************/

   if (argc != 3)
   {
       cout << "Please manually enter the following 2 filenames below:" << endl << endl
            << "1) Enter the name of the sequence file : >> ";
       cin  >> filename[0];
       cout << endl << "2) Enter the name of the intermediary file : >> ";
       cin  >> filename[1];
       cout << endl << endl << "Thank you! Information is being processed..." << endl << endl;

       if ((infile = fopen(filename[0], "a+")) == NULL){ // Changed to "append + read" mode
           cout << "Error opening sequence file.\n";
           exit (EXIT_FAILURE); 
       }
       if ((outfile = fopen(filename[1], "w")) == NULL){ 
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
   }   

   /************************* MEMORY ALLOCATION ******************************/

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
   if ((infile = freopen(argv[1], "r", infile)) == NULL){
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
   
   // Allocating memory for the matrices:
   matrix_forward = (int **) malloc ((D_ERRORS_PLUS_1) * sizeof(int *));
   matrix_reverse = (int **) malloc ((D_ERRORS_PLUS_1) * sizeof(int *));
   left_matrix_forward = (int **) malloc ((D_ERRORS_PLUS_1) * sizeof(int *));
   left_matrix_reverse = (int **) malloc ((D_ERRORS_PLUS_1) * sizeof(int *));
   right_matrix_forward = (int **) malloc ((D_ERRORS_PLUS_1) * sizeof(int *));
   right_matrix_reverse = (int **) malloc ((D_ERRORS_PLUS_1) * sizeof(int *));

   if (matrix_forward == NULL || matrix_reverse == NULL || left_matrix_forward == NULL ||
      left_matrix_reverse == NULL || right_matrix_forward == NULL || right_matrix_reverse == NULL) 
   { 
       cout << "Memory allocation failed.\n"; 
       exit (EXIT_FAILURE); 
   }
   
   for (i = 0; i < D_ERRORS_PLUS_1; i++)
   {
      temp = (i <= MAX_ERRORS) ? i + 1 : D_ERRORS_PLUS_1 - i + 1;
      matrix_forward[i] = (int *) malloc (temp * sizeof(int));
      matrix_reverse[i] = (int *) malloc (temp * sizeof(int));
      left_matrix_forward[i] = (int *) malloc (temp * sizeof(int));
      left_matrix_reverse[i] = (int *) malloc (temp * sizeof(int));
      right_matrix_forward[i] = (int *) malloc (temp * sizeof(int));
      right_matrix_reverse[i] = (int *) malloc (temp * sizeof(int));
      if (matrix_forward[i] == NULL || matrix_reverse[i] == NULL || left_matrix_forward[i] == NULL ||
          left_matrix_reverse[i] == NULL || right_matrix_forward[i] == NULL || right_matrix_reverse[i] == NULL)
      { 
          cout << "Memory allocation failed.\n"; 
          exit (EXIT_FAILURE);
      }
   }

   // Allocating memory for the entries:
   down_errors = (Entry*) malloc((MAX_ERRORS + 1)*sizeof(Entry));
   right_errors = (Entry*) malloc((MAX_ERRORS + 1)*sizeof(Entry));
   left_down_errors = (Entry*) malloc((MAX_ERRORS + 1)*sizeof(Entry));
   left_right_errors = (Entry*) malloc((MAX_ERRORS + 1)*sizeof(Entry));
   right_down_errors = (Entry*) malloc((MAX_ERRORS + 1)*sizeof(Entry));
   right_right_errors = (Entry*) malloc((MAX_ERRORS + 1)*sizeof(Entry));

   if (down_errors == NULL || right_errors == NULL || left_down_errors == NULL 
       || left_right_errors == NULL || right_down_errors == NULL || right_right_errors == NULL)
   {  cout << "Error allocating errors array.\n" << endl;  exit (EXIT_FAILURE);  }

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
   lower = (int *) malloc ((D_ERRORS_PLUS_1) * sizeof(int));
   upper = (int *) malloc ((D_ERRORS_PLUS_1) * sizeof(int));

   if (init_increments == NULL || upper == NULL || lower == NULL)
   {
       cout << "Memory allocation error in main.\n";
       exit (EXIT_FAILURE);
   }

   // Initialize the increments array:
   for (i = MAX_ERRORS; i < D_ERRORS_PLUS_1; i++)
     init_increments[i] = i - MAX_ERRORS;

   // Miriam Briskman, 03.02.2020:
   forward_pattern1 = (char*) malloc ((QUAD_PERIOD + 1)*sizeof(char)), 
   forward_pattern2 = (char*) malloc ((QUAD_PERIOD + 1)*sizeof(char));
   if (forward_pattern1 == NULL || forward_pattern2 == NULL)
   {
       printf ("Memory allocation error in main.\n");
       exit (EXIT_FAILURE);
   }

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

   /********************DECLARE AND INITIALIZE STRUCTS USING INITIALIZATION LISTS***********************/
    
   //the parameters are {j, left, right, length, rating, left_errs, right_errs, total_errs, 
   //                    matrix_forward, matrix_reverse, offset, down_errors,right_errors,newhigh}
    
   LastReportedRepeat left = {0, 0, 0, 0, 0, 0, 0, 0, left_matrix_forward, left_matrix_reverse, 
                              0, left_down_errors, left_right_errors, false};
    
   LastReportedRepeat right = {0, 0, 0, 0, 0, 0, 0, 0, right_matrix_forward, right_matrix_reverse,
                               0, right_down_errors, right_right_errors, false};

   /****************** ALGORITHMIC PORTION OF MAIN FUNCTION ********************/

   // Performs the (usually) recursive part of the Main-Lorentz algorithm
   // iteratively. The algorithm is modified so as to work with overlapping
   // parts, allowing it to process longer input strings.

   to_be_read = wholestr_len;
   partialstr = wholestr;
   strpnt = partialstr;

   numarray[0] = DOUBLE_PERIOD;
   i = i/2 - 1;
   levels--;
   for (j = 0; j < i; j++) // Copied from old 'createarray.cpp' by Miriam Briskman, 02.23.2020
   {
       numarray[j*2 + 1] = numarray[j]/2;
       numarray[j*2 + 2] = numarray[j] - numarray[j*2 + 1];
   }

   to_be_read -= DOUBLE_PERIOD;

   while (to_be_read >= DOUBLE_PERIOD)
   {
        strpnt = partialstr;
        for (i = 0, j = 1; i < levels; i++)
        {
            if (PROCESSING)
                printf(".");
            for (k = 0; k < j; k++)
            {
                num = numarray[j + k - 1];
                OneIteration (strpnt, num, matrix_forward, matrix_reverse, down_errors, right_errors,
                              left, right, offset + (long)(strpnt - partialstr), 0, outfile);
                strpnt += num;
            }
            j *= 2;
            strpnt = partialstr;
        }
        offset += DOUBLE_PERIOD;
        
        OneIteration(partialstr, QUAD_PERIOD, matrix_forward, matrix_reverse, down_errors, 
                     right_errors, left, right, offset - DOUBLE_PERIOD, 0, outfile);
        fprintf (outfile, "-1 %ld\n\n", offset);

        partialstr += DOUBLE_PERIOD;
        to_be_read -= DOUBLE_PERIOD;
   } // End of while loop

   OneIteration (partialstr, DOUBLE_PERIOD + to_be_read, matrix_forward, matrix_reverse, down_errors, 
                 right_errors, left, right, offset, 1, outfile);

   strpnt = partialstr;

   // Recompute the arrays of numarray:
   levels = 0;
   j = (DOUBLE_PERIOD + to_be_read - 1)/MIN_LENGTH + 1;
   for (i = 1; i < j; i *= 2)
       levels++;

   numarray[0] = DOUBLE_PERIOD + to_be_read;
   for (j = 0; j < i - 1; j++)
   {
       numarray[j*2 + 1] = numarray[j]/2;
       numarray[j*2 + 2] = numarray[j] - numarray[j*2 + 1];
   }

   // Final loop:
   for (i = 1, j = 2; i < levels; i++)
   {
       if (PROCESSING)
           printf(".");
       for (k = 0; k < j; k++){
           num = numarray[j + k - 1];
           OneIteration(strpnt, num, matrix_forward, matrix_reverse, down_errors, right_errors,
                        left, right, offset + (long)(strpnt - partialstr), k==j-1?1:0, outfile);
           strpnt += num;
       }
       j *= 2;
       strpnt = partialstr;
   }
   fprintf (outfile,"-1 -1\n\n");

   /************************ MEMORY DEALLOCATION *****************************/

   free (init_increments); free (upper); free (lower);
   free (numarray); free (wholestr);
   free(forward_pattern1); free(forward_pattern2);
   for (i = 0; i < D_ERRORS_PLUS_1; i++)
   {
      free (matrix_forward[i]); free (matrix_reverse[i]);
      free (left_matrix_forward[i]); free (left_matrix_reverse[i]);
      free (right_matrix_forward[i]); free (right_matrix_reverse[i]);
   }
   free (matrix_forward); free (matrix_reverse);
   free (left_matrix_forward); free (left_matrix_reverse);
   free (right_matrix_forward); free (right_matrix_reverse);
   free (down_errors); free (right_errors);
   free (left_down_errors); free (left_right_errors);
   free (right_down_errors); free (right_right_errors);

   fclose (infile);
   fclose (outfile);

   // Capture the ending time of the program:
   auto end_time = std::chrono::high_resolution_clock::now();
   auto mtime = end_time - start_time;

   cout << "Total execution time: " << mtime/std::chrono::milliseconds(1) << " milliseconds." << endl;

   return EXIT_SUCCESS; // END MAIN
}

/************************** READ_FILE FUNCTION + FASTA ******************************/
/************************* Miriam Briskman, 02.23.2020 ******************************/

// Version of fgets that removes all white space and FASTA comments.
// Reads the entire file. Does not stop at a new line.
// A null byte is appended to s to mark the end of the string.
// Assumptions:
// 1) There is always a newline character immediately after a FASTA comment.
// 2) A base {A, C, G, T, N} is followed immediately only by a base or a newline character.
void read_file (char* s, unsigned long long &numRead, FILE* stream)
{
    char next;

    while (next != EOF)
    {
        next = getc(stream);
        switch (next)
        {
            case '>':  while (getc(stream) != '\n');
            case '\n':
            case '\t':
            case '\r':
            case '\f':
            case '\v':
            case ' ':
            case EOF:  break;
            default:   s[numRead++] = next;
                       while ((s[numRead++] = getc(stream)) != '\n');
                       numRead--; // Discard the '\n' character.
                       break;
        }
    }
    s[numRead] = '\0';
}

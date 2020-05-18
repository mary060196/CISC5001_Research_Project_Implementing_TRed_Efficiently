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
#include "buildk.h"
#include "oneiteration.h"
#include "parameters.h"
#include "NewMatrixL.h"
#include "tracek.h"
#include "aux_main.h"
#include <iostream>
#include <string>

using namespace std;

// Function called by the different threads in the program:
void threads_func (thread_info *my_info)
{
   FILE *outfile = NULL;
   char *partialstr = NULL, *strpnt = NULL, *forward_pattern1, *forward_pattern2;
   int d, i, j, k, l, levels = 0,
       **matrix_forward = NULL, **matrix_reverse = NULL,
       **left_matrix_forward = NULL, **left_matrix_reverse = NULL,
       **right_matrix_forward = NULL, **right_matrix_reverse = NULL,
       *upper, *lower, num, temp;
   unsigned long long offset; // Miriam Briskman, 02.23.2020
   long long times; 
   Entry *down_errors = NULL, *right_errors = NULL, *left_down_errors = NULL,
         *left_right_errors = NULL, *right_down_errors = NULL, *right_right_errors = NULL;

   char *joinedstr[2];
   suffix **suffixes[2];
   int **L, rows = 2*MAX_ERRORS + 3, cols = MAX_ERRORS + 3,
       *lcp[4], *sizeToIndex[4], *sizeToIndexFixed[4], *smallest01[4][2];
   bool *isSuffixConsidered[4];

   /* As a reminder,

            struct thread_info {
                char *beg_str;
                string filename;
                long long times;
                unsigned long long offset;
            }; 
   */

   /***************** USING my_info TO ASSIGN VARIABLES **********************/

   partialstr = my_info->beg_str;
   times = my_info->times;
   offset = my_info->offset;

   /****************** OPENING THE THREAD'S OUTPUT FILE **********************/

   if ((outfile = fopen((my_info->filename).c_str(), "w")) == NULL)
   { 
           cout << "Error opening intermediary file.\n";
           exit (EXIT_FAILURE);
   }

   /************************* MEMORY ALLOCATION ******************************/

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

   lower = (int *) malloc ((D_ERRORS_PLUS_1) * sizeof(int));
   upper = (int *) malloc ((D_ERRORS_PLUS_1) * sizeof(int));

   if (upper == NULL || lower == NULL)
   {
       cout << "Memory allocation error in main.\n";
       exit (EXIT_FAILURE);
   }

   // Miriam Briskman, 03.02.2020:
   forward_pattern1 = (char*) malloc ((QUAD_PERIOD + 1)*sizeof(char)), 
   forward_pattern2 = (char*) malloc ((QUAD_PERIOD + 1)*sizeof(char));
   if (forward_pattern1 == NULL || forward_pattern2 == NULL)
   {
       printf ("Memory allocation error in main.\n");
       exit (EXIT_FAILURE);
   }

   // Allocating memory for the L matrix:
   L = (int**) malloc (rows * sizeof(int*));
   if (L == NULL)
   {
      printf ("Memory allocation error in 'aux_main'.\n");
      exit (EXIT_FAILURE);
   }
   for (i = 0; i < rows; ++i)
   {
       L[i] = (int*) calloc (cols, sizeof(int));
       if (L[i] == NULL)
       {
         printf ("Memory allocation error in 'aux_main'.\n");
         exit (EXIT_FAILURE);
       }
   }

   // Initializing the matrix according to the Landau-Vishkin Paper:
   for (d = 0; d <= MAX_ERRORS; d++)
   {
      L[d][MAX_ERRORS - d + 2] = MAX_ERRORS - d;
      L[d][MAX_ERRORS - d + 1] = -2;
   }
   for (d = MAX_ERRORS + 1; d <= MAX_ERRORS*2 + 2; d++)
   {
       L[d][d - MAX_ERRORS] = -1;
       L[d][d - MAX_ERRORS - 1] = -2;
   }

   joinedstr[0] = (char*) malloc (JOINED_LEN * sizeof(char));
   joinedstr[1] = (char*) malloc (JOINED_LEN * sizeof(char));

   suffixes[0] = (suffix**) malloc (JOINED_LEN * sizeof(suffix*));
   suffixes[1] = (suffix**) malloc (JOINED_LEN * sizeof(suffix*));

   if (suffixes[0] == NULL || suffixes[1] == NULL 
       || joinedstr[0] == NULL || joinedstr[1] == NULL)
   {
      printf ("Memory allocation issues inside 'prepareSuffixArrays'.\n");
      exit (EXIT_FAILURE);
   }

   for (i = 0; i < JOINED_LEN; i++)
   {
       suffixes[0][i] = (suffix*) malloc (sizeof(suffix)); // One suffix instance
       suffixes[1][i] = (suffix*) malloc (sizeof(suffix)); // One suffix instance
       if (suffixes[0][i] == NULL || suffixes[1][i] == NULL)
       {
          printf ("Memory allocation issues inside 'prepareSuffixArrays'.\n");
          exit (EXIT_FAILURE);
       }
   }

   for (i = 0; i < 4; i++)
   {
      lcp[i] = (int*) malloc (JOINED_LEN * sizeof(int));
      sizeToIndex[i] = (int*) malloc (JOINED_LEN * sizeof(int));
      sizeToIndexFixed[i] = (int*) malloc (JOINED_LEN * sizeof(int));
      isSuffixConsidered[i] = (bool*) malloc (JOINED_LEN * sizeof(bool));
      smallest01[i][0] = (int*) malloc (JOINED_LEN * sizeof(int));
      smallest01[i][1] = (int*) malloc (JOINED_LEN * sizeof(int));
      if (lcp[i] == NULL || sizeToIndex[i] == NULL || sizeToIndexFixed[i] == NULL 
          || isSuffixConsidered[i] == NULL || smallest01[i][0] == NULL || smallest01[i][1] == NULL)
      {
          printf ("Memory allocation issues inside 'aux_main'.\n");
          exit (EXIT_FAILURE);
      }
   }

   /******************* DECLARE AND INITIALIZE STRUCTS USING INITIALIZATION LISTS **********************/
    
   // Miriam Briskman, 04.17.2020 (Expanding initialization list)
   //the parameters are {j, left, right, length, rating, NoMinusOnesUpperRowForward, NoMinusOnesLowerRowBackward,
   //                    left_errs, right_errs, total_errs, 
   //                    matrix_forward, matrix_reverse, offset, down_errors,right_errors,newhigh}
    
   LastReportedRepeat left = {0, 0, 0, 0, 0, D_ERRORS_PLUS_1, 0,
                              0, 0, 0, 
                              left_matrix_forward, left_matrix_reverse, 
                              0, left_down_errors, left_right_errors, false};
    
   LastReportedRepeat right = {0, 0, 0, 0, 0, D_ERRORS_PLUS_1, 0,
                               0, 0, 0, 
                               right_matrix_forward, right_matrix_reverse,
                               0, right_down_errors, right_right_errors, false};

   /****************** ALGORITHMIC PORTION OF MAIN FUNCTION ********************/

   // Performs the (usually) recursive part of the Main-Lorentz algorithm
   // iteratively. The algorithm is modified so as to work with overlapping
   // parts, allowing it to process longer input strings.

   j = (QUAD_PERIOD - 1)/MIN_LENGTH + 1;
   for (i = 1; i < j; i *= 2)
       levels++;
   levels--;

   // Skip the regions full with UNKNOWNs, as these never form
   //    tandem repeats. In a case of a DNA sequence with frequent
   //    unknown regions, this would save a great amount of work.
   // The regions will be skipped only if they are multiples of the
   //    length of DOUBLE_PERIOD.
   // Miriam Briskman, 04.17.2020
   while (*partialstr == UNKNOWN)
   {
       for (l = 1; l < DOUBLE_PERIOD && partialstr[l] == UNKNOWN; l++);
       if (l < DOUBLE_PERIOD)
          break;
       else
       {
           offset += DOUBLE_PERIOD;
           fprintf (outfile, "-1 %llu\n\n", offset);
           partialstr += DOUBLE_PERIOD;
           --times;
       }
   }

   for (int g = 0; g < times; g++)
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
                              left, right, offset + (long)(strpnt - partialstr), upper, lower,
                              forward_pattern1, forward_pattern2, outfile, joinedstr, suffixes,
                              L, lcp, sizeToIndex, sizeToIndexFixed, smallest01, isSuffixConsidered);
                strpnt += num;
            }
            j *= 2;
            strpnt = partialstr;
        }
        offset += DOUBLE_PERIOD;
        
        OneIteration(partialstr, QUAD_PERIOD, matrix_forward, matrix_reverse, down_errors, 
                     right_errors, left, right, offset - DOUBLE_PERIOD, upper, lower,
                     forward_pattern1, forward_pattern2, outfile, joinedstr, suffixes,
                     L, lcp, sizeToIndex, sizeToIndexFixed, smallest01, isSuffixConsidered);
        fprintf (outfile, "-1 %llu\n\n", offset);

        partialstr += DOUBLE_PERIOD;

        // Skip the regions full with UNKNOWNs, as these never form
        //    tandem repeats. In a case of a DNA sequence with frequent
        //    unknown regions, this would save a great amount of work.
        // The regions will be skipped only if they are of multiples of the
        //    length of DOUBLE_PERIOD.
        // Miriam Briskman, 04.17.2020
        while (*partialstr == UNKNOWN)
        {
            for (l = 1; l < DOUBLE_PERIOD && partialstr[l] == UNKNOWN; l++);

            if (l < DOUBLE_PERIOD)
               break;
            else
            {
                offset += DOUBLE_PERIOD;
                fprintf (outfile, "-1 %llu\n\n", offset);
                partialstr += DOUBLE_PERIOD;
                --times;
            }
        }
   } // End of while loop

   /************************ MEMORY DEALLOCATION *****************************/

   free (upper); free (lower);
   free (forward_pattern1); free (forward_pattern2);
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

   for (i = 0; i < rows; ++i)
     free (L[i]);
   free (L);
   free (joinedstr[0]); free (joinedstr[1]);
   for (i = 0; i < JOINED_LEN; i++)
   {
      free (suffixes[0][i]); free (suffixes[1][i]);
   }
   free (suffixes[0]); free (suffixes[1]);
   for (i = 0; i < 4; i++)
   {
      free (lcp[i]); free (sizeToIndex[i]); 
      free (sizeToIndexFixed[i]); free (isSuffixConsidered[i]); 
      free (smallest01[i][0]); free (smallest01[i][1]);
   }
   

   fclose (outfile);
}

// Function called by the last thread in the program:
void threads_func_last (thread_info *my_info)
{
   FILE *outfile = NULL;
   char *partialstr = NULL, *strpnt = NULL, *forward_pattern1, *forward_pattern2;
   int d, i, j, k, l, levels = 0,
       *numarrayBig = NULL,
       **matrix_forward = NULL, **matrix_reverse = NULL,
       **left_matrix_forward = NULL, **left_matrix_reverse = NULL,
       **right_matrix_forward = NULL, **right_matrix_reverse = NULL,
       *upper, *lower, num, temp;
   unsigned long long offset; // Miriam Briskman, 02.23.2020
   long long to_be_read,  // The number of characters remaining until end of wholestr.
             times; 
   Entry *down_errors = NULL, *right_errors = NULL, *left_down_errors = NULL,
         *left_right_errors = NULL, *right_down_errors = NULL, *right_right_errors = NULL;

   char *joinedstr[2];
   suffix **suffixes[2];
   int **L, rows = 2*MAX_ERRORS + 3, cols = MAX_ERRORS + 3,
       *lcp[4], *sizeToIndex[4], *sizeToIndexFixed[4], *smallest01[4][2];
   bool *isSuffixConsidered[4];

   /* As a reminder,

            struct thread_info {
                char *beg_str;
                string filename;
                long long times;
                unsigned long long offset;
            }; 
   */

   /***************** USING my_info TO ASSIGN VARIABLES **********************/

   partialstr = my_info->beg_str;
   times = my_info->times;
   offset = my_info->offset;

   /****************** OPENING THE THREAD'S OUTPUT FILE **********************/

   if ((outfile = fopen((my_info->filename).c_str(), "w")) == NULL)
   { 
           cout << "Error opening intermediary file.\n";
           exit (EXIT_FAILURE);
   }

   /************************* MEMORY ALLOCATION ******************************/

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

   lower = (int *) malloc ((D_ERRORS_PLUS_1) * sizeof(int));
   upper = (int *) malloc ((D_ERRORS_PLUS_1) * sizeof(int));

   if (upper == NULL || lower == NULL)
   {
       cout << "Memory allocation error in main.\n";
       exit (EXIT_FAILURE);
   }

   // Miriam Briskman, 03.02.2020:
   forward_pattern1 = (char*) malloc ((QUAD_PERIOD + 1)*sizeof(char)), 
   forward_pattern2 = (char*) malloc ((QUAD_PERIOD + 1)*sizeof(char));
   if (forward_pattern1 == NULL || forward_pattern2 == NULL)
   {
       printf ("Memory allocation error in main.\n");
       exit (EXIT_FAILURE);
   }

   // Allocating memory for the L matrix:
   L = (int**) malloc (rows * sizeof(int*));
   if (L == NULL)
   {
      printf ("Memory allocation error in 'aux_main'.\n");
      exit (EXIT_FAILURE);
   }
   for (i = 0; i < rows; ++i)
   {
       L[i] = (int*) calloc (cols, sizeof(int));
       if (L[i] == NULL)
       {
         printf ("Memory allocation error in 'aux_main'.\n");
         exit (EXIT_FAILURE);
       }
   }

   // Initializing the matrix according to the Landau-Vishkin Paper:
   for (d = 0; d <= MAX_ERRORS; d++)
   {
      L[d][MAX_ERRORS - d + 2] = MAX_ERRORS - d;
      L[d][MAX_ERRORS - d + 1] = -2;
   }
   for (d = MAX_ERRORS + 1; d <= MAX_ERRORS*2 + 2; d++)
   {
       L[d][d - MAX_ERRORS] = -1;
       L[d][d - MAX_ERRORS - 1] = -2;
   }

   joinedstr[0] = (char*) malloc (JOINED_LEN * sizeof(char));
   joinedstr[1] = (char*) malloc (JOINED_LEN * sizeof(char));

   suffixes[0] = (suffix**) malloc (JOINED_LEN * sizeof(suffix*));
   suffixes[1] = (suffix**) malloc (JOINED_LEN * sizeof(suffix*));

   if (suffixes[0] == NULL || suffixes[1] == NULL 
       || joinedstr[0] == NULL || joinedstr[1] == NULL)
   {
      printf ("Memory allocation issues inside 'prepareSuffixArrays'.\n");
      exit (EXIT_FAILURE);
   }

   for (i = 0; i < JOINED_LEN; i++)
   {
       suffixes[0][i] = (suffix*) malloc (sizeof(suffix)); // One suffix instance
       suffixes[1][i] = (suffix*) malloc (sizeof(suffix)); // One suffix instance
       if (suffixes[0][i] == NULL || suffixes[1][i] == NULL)
       {
          printf ("Memory allocation issues inside 'prepareSuffixArrays'.\n");
          exit (EXIT_FAILURE);
       }
   }

   for (i = 0; i < 4; i++)
   {
      lcp[i] = (int*) malloc (JOINED_LEN * sizeof(int));
      sizeToIndex[i] = (int*) malloc (JOINED_LEN * sizeof(int));
      sizeToIndexFixed[i] = (int*) malloc (JOINED_LEN * sizeof(int));
      isSuffixConsidered[i] = (bool*) malloc (JOINED_LEN * sizeof(bool));
      smallest01[i][0] = (int*) malloc (JOINED_LEN * sizeof(int));
      smallest01[i][1] = (int*) malloc (JOINED_LEN * sizeof(int));
      if (lcp[i] == NULL || sizeToIndex[i] == NULL || sizeToIndexFixed[i] == NULL 
          || isSuffixConsidered[i] == NULL || smallest01[i][0] == NULL || smallest01[i][1] == NULL)
      {
          printf ("Memory allocation issues inside 'aux_main'.\n");
          exit (EXIT_FAILURE);
      }
   }

   /******************* DECLARE AND INITIALIZE STRUCTS USING INITIALIZATION LISTS **********************/
    
   // Miriam Briskman, 04.17.2020 (Expanding initialization list)
   //the parameters are {j, left, right, length, rating, NoMinusOnesUpperRowForward, NoMinusOnesLowerRowBackward,
   //                    left_errs, right_errs, total_errs, 
   //                    matrix_forward, matrix_reverse, offset, down_errors,right_errors,newhigh}
    
   LastReportedRepeat left = {0, 0, 0, 0, 0, D_ERRORS_PLUS_1, 0,
                              0, 0, 0, 
                              left_matrix_forward, left_matrix_reverse, 
                              0, left_down_errors, left_right_errors, false};
    
   LastReportedRepeat right = {0, 0, 0, 0, 0, D_ERRORS_PLUS_1, 0,
                               0, 0, 0, 
                               right_matrix_forward, right_matrix_reverse,
                               0, right_down_errors, right_right_errors, false};

   /****************** ALGORITHMIC PORTION OF MAIN FUNCTION ********************/

   // Performs the (usually) recursive part of the Main-Lorentz algorithm
   // iteratively. The algorithm is modified so as to work with overlapping
   // parts, allowing it to process longer input strings.

   j = (QUAD_PERIOD - 1)/MIN_LENGTH + 1;
   for (i = 1; i < j; i *= 2)
       levels++;
   levels--;

   to_be_read = wholestr_len - times*(THREAD_NUM-1)*DOUBLE_PERIOD;

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
                              left, right, offset + (long)(strpnt - partialstr), upper, lower,
                              forward_pattern1, forward_pattern2, outfile, joinedstr, suffixes,
                              L, lcp, sizeToIndex, sizeToIndexFixed, smallest01, isSuffixConsidered);
                strpnt += num;
            }
            j *= 2;
            strpnt = partialstr;
        }
        offset += DOUBLE_PERIOD;
        
        OneIteration(partialstr, QUAD_PERIOD, matrix_forward, matrix_reverse, down_errors, 
                     right_errors, left, right, offset - DOUBLE_PERIOD, upper, lower,
                     forward_pattern1, forward_pattern2, outfile, joinedstr, suffixes,
                     L, lcp, sizeToIndex, sizeToIndexFixed, smallest01, isSuffixConsidered);

        fprintf (outfile, "-1 %llu\n\n", offset);

        partialstr += DOUBLE_PERIOD;
        to_be_read -= DOUBLE_PERIOD;

        // Skip the regions full with UNKNOWNs, as these never form
        //    tandem repeats. In a case of a DNA sequence with frequent
        //    unknown regions, this would save a great amount of work.
        // The regions will be skipped only if they are of multiples of the
        //    length of DOUBLE_PERIOD.
        // Miriam Briskman, 04.17.2020
        while (*partialstr == UNKNOWN)
        {
            for (l = 1; l < DOUBLE_PERIOD && partialstr[l] == UNKNOWN; l++);

            if (l < DOUBLE_PERIOD)
               break;
            else
            {
                offset += DOUBLE_PERIOD;
                fprintf (outfile, "-1 %llu\n\n", offset);
                partialstr += DOUBLE_PERIOD;
                to_be_read -= DOUBLE_PERIOD; // Might reduce to a negative number
            }
        }
   } // End of while loop

   if (to_be_read >= 0) // Check if there were no UNKNOWNs at the end of the file
   {
      OneIteration (partialstr, DOUBLE_PERIOD + to_be_read, matrix_forward, matrix_reverse, down_errors, 
                    right_errors, left, right, offset, upper, lower, forward_pattern1, forward_pattern2,
                    outfile, joinedstr, suffixes, L, lcp, sizeToIndex, sizeToIndexFixed, smallest01, 
                    isSuffixConsidered);

      strpnt = partialstr;

      // Recompute the arrays of numarray:
      levels = 0;
      j = (DOUBLE_PERIOD + to_be_read - 1)/MIN_LENGTH + 1;
      for (i = 1; i < j; i *= 2)
          levels++;

      numarrayBig = (int*) malloc ((i*2 - 1)*sizeof(int)); // Size sufficient for the entire program
      if (numarrayBig == NULL)
      {
          cout << "Memory allocation error in main.\n";
          exit (EXIT_FAILURE);
      }

      numarrayBig[0] = DOUBLE_PERIOD + to_be_read;
      for (j = 0; j < i - 1; j++)
      {
          numarrayBig[j*2 + 1] = numarrayBig[j]/2;
          numarrayBig[j*2 + 2] = numarrayBig[j] - numarrayBig[j*2 + 1];
      }

      // Final loop:
      for (i = 1, j = 2; i < levels; i++)
      {
          if (PROCESSING)
              printf(".");
          for (k = 0; k < j; k++){
              num = numarrayBig[j + k - 1];
              OneIteration(strpnt, num, matrix_forward, matrix_reverse, down_errors, right_errors,
                           left, right, offset + (long)(strpnt - partialstr), upper, lower,
                           forward_pattern1, forward_pattern2, outfile, joinedstr, suffixes,
                           L, lcp, sizeToIndex, sizeToIndexFixed, smallest01, isSuffixConsidered);
              strpnt += num;
          }
          j *= 2;
          strpnt = partialstr;
      }
   }
   fprintf (outfile,"-1 -1\n\n");

   /************************ MEMORY DEALLOCATION *****************************/

   free (upper); free (lower);
   free (numarrayBig);
   free (forward_pattern1); free (forward_pattern2);
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

   for (i = 0; i < rows; ++i)
     free (L[i]);
   free (L);
   free (joinedstr[0]); free (joinedstr[1]);
   for (i = 0; i < JOINED_LEN; i++)
   {
      free (suffixes[0][i]); free (suffixes[1][i]);
   }
   free (suffixes[0]); free (suffixes[1]);
   for (i = 0; i < 4; i++)
   {
      free (lcp[i]); free (sizeToIndex[i]); 
      free (sizeToIndexFixed[i]); free (isSuffixConsidered[i]); 
      free (smallest01[i][0]); free (smallest01[i][1]);
   }

   fclose (outfile);
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
    char next = 'a'; // Assign a char not equal to EOF to let the loop execute.

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
                       break;
        }
    }
    s[numRead] = '\0';
}

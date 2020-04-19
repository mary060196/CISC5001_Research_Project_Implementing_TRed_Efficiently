//   |=================================================================|
//   |   TRED: a tool for detecting Tandem Repeats within sequences,   |
//   |               using the Edit Distance metric.                   |
//   |=================================================================|

// Copyright © 2007-2009 Dina Sokol, Justin Tojeira
// Distributed under the Aladdin Free Public License

// THIS SOFTWARE SHOULD BE ACCOMPANIED BY readme.txt AND license.html
// WE STRONGLY ENCOURAGE YOU TO READ BOTH BEFORE PROCEEDING
//------------------------------------------------------------------------------



// This function performs a traceback of a KxK edit distance matrix created by
// (and described in the comments of) buildk.cpp

//           ~~~~~ COMPUTER REPRESENTATION OF THE KxK MATRIX ~~~~~

// It is important to note that the KxK matrix is represented in the computer
// differently from how we think of it.

//   In the computer                      How we think of it
//       0  1  2  3                            0  1  2  3
//   0   6                                 -3           6
//   1   5  6                              -2        5  6
//   2   4  5  6                           -1     4  5  6
//   3   1  2  5  6                         0  1  2  5  6
//   4   1  2  5                            1     1  2  5
//   5   2  3                               2        2  3
//   6   2                                  3           2

// We like to think of the main diagonal as "O", but the computer can't do this
// because it would result in negative array subscripting. Also, I shifted all
// the data to the beginning of each array to save space. I will describe the
// traceback algorithm in terms of our representation of the KxK matrix. In
// the code, the numbers may be a bit different but the algorithm it follows
// is the same.

//                   ~~~~~ TRACEBACK ALGORITHM ~~~~~

// During the traceback, we do not recreate any part of the NxN matrix. We do,
// however, deal with a "virtual" NxN matrix. We do this by tracking our
// position on each of the strings used to create the NxN matrix, and we use
// the position as a set of "coordinates".

//         A  C  D  B  E  H                  0  1  2  3
//      ?  ?  ?  ?  ?  ?  ?              -3           6
//   A  ?  ?  ?  ?  ?  ?  ?              -2        5  6
//   B  ?  ?  ?  ?  ?  ?  ?              -1     4  5  6
//   C  ?  ?  ?  ?  *  ?  ?               0  1  2  5  6
//   D  ?  ?  ?  ?  ?  ?  ?               1     1  2  5
//   E  ?  ?  ?  ?  ?  ?  ?               2        2  3
//   F  ?  ?  ?  ?  ?  ?  ?               3           2

// At the position indicated by the *, position 2 on S1 and 3 on S2, we can
// tell that we're on diagonal 1. Furthermore, our row in the "virtual" NxN
// matrix is 3 (remember that the NxN matrix begins 1 before the strings do).
// Were we to start at this point, we would begin by checking row 1 in the KxK
// matrix to find the column containing the lowest number which is greater than
// or equal to 3. In this case, it is the column 3, which holds the number 5.
// So we begin in the KxK matrix in row 1 and column 3.

// From this position in the KxK matrix, if the characters we're up to in S1
// and S2 match up, it signals that we can move diagonally in the "virtual"
// NxN matrix without changing the number of errors. To the KxK matrix, this
// means no change, as the row is the same (same diagonal on the NxN) and the
// errors are the same. If the characters we're up to in S1 and S2 don't match
// up, a move must be made in the KxK matrix. This move can be either left
// (representing a diagonal move in the NxN, and 1 less error), up and left
// (a left move in the NxN) or down and left (an up move in the NxN). In any
// case, a test must be made that the position we're moving to in the virtual
// NxN matrix has the correct number of errors. We do this by making sure that
// the value of the position we're moving to in the KxK matrix is greater than
// or equal to the row we're moving to on the virtual NxN matrix. Note that the
// row in the virtual NxN matrix is equal to the position on S1, plus 1.

//                    ~~~~~ FUNCTION SPECIFICATIONS ~~~~~
//    THE PARAMETERS HAVE BEEN MODIFIED IN THIS VERSION TO FIT WITH MAIN

// Traceback takes: a KxK matrix (m), the number of errors used to create the
// matrix (max_errors, this needs to be known as it determines the size of the
// matrix), the two strings used to create the matrix (s1 and s2), and pointers
// to the beginning of two strings, o1 and o2, to contain the aligned strings.

// It is assumed that the traceback will begin at the last letter of each input
// string. If this is not the case, then the user must also enter the position
// in each string to start at (s1pos and s2pos).

// This function depends on its parameters matching up and does not have
// safeguards for all possible invalid cases. The parameters max_errors, s1 and
// s2 should be the same parameters used to create the matrix m. Furthermore,
// the points of origin (s1pos and s2pos), whether specified or left to default,
// should be points that that the function buildmatrix was able to reach in the
// allowed amount of errors or fewer.

//                      ~~~~~ POINTS OF INTEREST ~~~~~

// Traceback does not actually recreate any part of the NxN matrix.
// Traceback outputs the aligned strings, output1 and output2, backwards (I'm
// pretty sure there was some reason to do it this way but I can't remember
// what it was).

//                   ~~~~~ VERSION 2 MODIFICATIONS ~~~~~

// A macro "UNKNOWN" has been added for usage with biological sequences.
// UNKNOWN should be defined as the character that represents an unknown element
// of the sequence (For example, the default is 'N' which represents an unknown
// base pair in a DNA sequence). This character will always be considered a
// mismatch, even against itself. UNKNOWN is defined in buildk.h

#include "tracek.h"
#include "buildk.h"
#include "parameters.h" 
#include <stdio.h>

#define max_errors MAX_ERRORS // Miriam Briskman, 03.03.2020
#define ABS_ROW_ERRORS ((row > max_errors) ? row - max_errors : max_errors - row) // Miriam Briskman, 02.23.2020

/************************** TRACEALIGN *****************************/
/***************** Miriam Briskman, 03.02.2020 *********************/
void tracealign (int** m, char* s1, char* s2, int lowerRow, int s1pos, int s2pos, 
                 long &pos1, long &pos2, bool &state, FILE* outfile)
{

/************************* VARIABLE INITIALIZATION ****************************/
   // row and col refer to the computer representation of the KxK matrix.
   int row = s2pos - s1pos + max_errors,
       col = max_errors - ABS_ROW_ERRORS,
       init_s1pos = s1pos, init_s2pos = s2pos;

/******************** CHECK TO MAKE SURE ROW IS IN BOUNDS *********************/

   if (row < 0 || row > 2*max_errors)
   {
        printf("\n\n          *!*!*!*!*!*!*!*!* ERROR *!*!*!*!*!*!*!*!*!*\n\n");
        printf("          Traceback cannot be started from that point\n");
        printf("          *!*!*!*!*!*!*!*!* ERROR *!*!*!*!*!*!*!*!*!*\n\n");
        return;
   }

   fprintf (outfile, "%ld %ld ", pos1, pos2);

/************** DECREMENTS COL WHEN THE ERRORS WE'RE ALLOWING ARE *************/
/************** LESS THAN THE ACTUAL ERRORS AT OUR STARTING POINT *************/

   while (col > 0 && m[row][col - 1] > s1pos)
        col--;

/*************************** MAIN TRACEBACK LOOP ******************************/

                // LOOPS UNTIL WE ARE ON THE LEFTMOST COLUMN
                // OR TOP ROW OF THE VIRTUAL NxN MATRIX
   while (s1pos >= 0 && s2pos >= 0)
   {
      if (row > max_errors && m[row - 1][col] > s1pos)
      {
          row--;                         // left move in the virtual NxN
          s2pos--;                       // up and left move in our KxK
          state = 1;                     // up move in computer's KxK
      }
      else if (row < max_errors && m[row + 1][col] >= s1pos)
      {
          row++;                     // up move in the virtual NxN
          s1pos--;                   // down and left move in our KxK
          state = 1;                 // down move in computer's KxK
      }
      else if (row <= max_errors && row >= lowerRow && col > 1 && m[row - 1][col - 2] > s1pos)
      {
          row--;                         // left move in the virtual NxN
          col -= 2;                      // up and left move in our KxK
          s2pos--;                       // up and 2 left move in computer's KxK
          state = 1;
      }
      else if (row >= max_errors && col > 1 && m[row + 1][col - 2] >= s1pos)
      {
          row++;                       // up move in the virtual NxN
          col -= 2;                    // down and left move in our KxK
          s1pos--;                     // down and 2 left move in computer's KxK
          state = 1;
      }
      else if (col > 0 && m[row][col - 1] >= s1pos)
      {
          if (state == 1)
          {
               fprintf (outfile, "%ld %ld ", pos1 + init_s2pos - s2pos, pos2 + init_s1pos - s1pos);
               state = 0;
          }
          col--;                          // diagonal move in the virtual NxN
          s1pos--;                  	  // left move in our KxK
          s2pos--;                        // left move in the computer's KxK
      }
      else
      {                                  
          if (state == 1)
          {
               fprintf (outfile, "%ld %ld ", pos1 + init_s2pos - s2pos, pos2 + init_s1pos - s1pos);
               state = 0;
          }                              // in the case of the characters we're up to in S1 and 
          s1pos--;                       // S2 matching, we make a diagonal move in the virtual NxN and
          s2pos--;                       // no move in either representation of the KxK
      }
   }

   /************************ FINISHES UP THE TRACEBACK ***************************/

   pos1 += init_s2pos - s2pos;
   pos2 += init_s1pos - s1pos;

   if (s1pos >= 0 || s2pos >= 0)
       state = 1;
   if (s1pos >= 0)
       pos2 += s1pos + 1;
   if (s2pos >= 0)
       pos1 += s2pos + 1;

   return;
} // End of tracealign

/************************* TRACEFORWARD ****************************/
/***************** Miriam Briskman, 02.23.2020 *********************/
void traceforward (int** m, char* s1, char* s2, int upperRow, int s1pos, int s2pos, 
                   long pos1, long pos2, bool state, FILE* outfile)
{

/************************* VARIABLE INITIALIZATION ****************************/
   // row and col refer to the computer representation of the KxK matrix.
   // o1pos and o2pos are used to keep position on output1 and output2.
   int row = s2pos - s1pos + max_errors,
       col = max_errors - ABS_ROW_ERRORS,
       oPos = 0;

/******************** CHECK TO MAKE SURE ROW IS IN BOUNDS *********************/

   if (row < 0 || row > 2*max_errors)
   {
        printf("\n\n          *!*!*!*!*!*!*!*!* ERROR *!*!*!*!*!*!*!*!*!*\n\n");
        printf("          Traceback cannot be started from that point\n");
        printf("          *!*!*!*!*!*!*!*!* ERROR *!*!*!*!*!*!*!*!*!*\n\n");
        return;
   }

/************** DECREMENTS COL WHEN THE ERRORS WE'RE ALLOWING ARE *************/
/************** LESS THAN THE ACTUAL ERRORS AT OUR STARTING POINT *************/

   while (col > 0 && m[row][col - 1] > s1pos)
       col--;

/*************************** MAIN TRACEBACK LOOP ******************************/

                // LOOPS UNTIL WE ARE ON THE LEFTMOST COLUMN
                // OR TOP ROW OF THE VIRTUAL NxN MATRIX
   while (s1pos >= 0 && s2pos >= 0)
   {
      if (row > max_errors && m[row - 1][col] > s1pos)
      {
          forward_pattern1[oPos] = '-';
          forward_pattern2[oPos] = '\0'; // left move in the virtual NxN
          row--;                         // up and left move in our KxK
          s2pos--;                       // up move in computer's KxK
      }
      else if (row < max_errors && m[row + 1][col] >= s1pos)
      {
          forward_pattern1[oPos] = '\0';
          forward_pattern2[oPos] = '-'; // up move in the virtual NxN
          row++;                        // down and left move in our KxK
          s1pos--;                      // down move in computer's KxK
      }
      else if (row <= max_errors && col > 1 && m[row - 1][col - 2] > s1pos)
      {
          forward_pattern1[oPos] = '-';
          forward_pattern2[oPos] = '\0'; // left move in the virtual NxN
          row--;                         // up and left move in our KxK
          col -= 2;                      // up and 2 left move in computer's KxK
          s2pos--;
      }
      else if (row >= max_errors && row <= upperRow && col > 1 && m[row + 1][col - 2] >= s1pos)
      {
          forward_pattern1[oPos] = '\0';
          forward_pattern2[oPos] = '-'; // up move in the virtual NxN
          row++;                        // down and left move in our KxK
          col -= 2;                     // down and 2 left move in computer's KxK
          s1pos--;
      }
      else if (col > 0 && m[row][col - 1] >= s1pos)
      {
          forward_pattern1[oPos] = '\0';
          forward_pattern2[oPos] = '\0';
          col--;                          // diagonal move in the virtual NxN
          s1pos--;                  	  // left move in our KxK
          s2pos--;                        // left move in the computer's KxK
      }
   	  else
      {                                  
          forward_pattern1[oPos] = '\0';
          forward_pattern2[oPos] = '\0'; // in the case of the characters we're up
          s1pos--;                       // to in S1 and S2 matching, we make a
          s2pos--;                       // diagonal move in the virtual NxN and
                                         // no move in either representation of the KxK
      }
      oPos++;
   }                 

/************************ FINISHES UP THE TRACEBACK ***************************/

           // TRACES UP ALONG THE LEFTMOST COLUMN OF THE NxN MATRIX IF NEEDED
   while (s1pos >= 0)
   {
      forward_pattern1[oPos] = '\0';
      forward_pattern2[oPos] = '-';
      s1pos--;
      oPos++;
   }
           // TRACES LEFT ALONG THE TOP ROW OF THE NxN MATRIX IF NEEDED
   while (s2pos >= 0)
   {
      forward_pattern1[oPos] = '-';
      forward_pattern2[oPos] = '\0';
      s2pos--;
      oPos++;
   }

   oPos--;

           // PRINT TO OUTPUT FILE IN SPECIAL FORMAT FOR POST-PROCESSING
   while (oPos >= 0)
   {
       if (forward_pattern2[oPos] == '-')
       {
           pos2++;
           state = 1;
       }
       else if (forward_pattern1[oPos] == '-')
       {
           pos1++;
           state = 1;
       }
       else 
       {
           if (state == 1)
           {
               fprintf (outfile, "%ld %ld ", pos1, pos2);
               state = 0;
           }
           pos1++;
           pos2++;
       }
       oPos--;
   }

   fprintf (outfile, "%ld 0\n\n", pos1);
   return;
} // End of traceforward

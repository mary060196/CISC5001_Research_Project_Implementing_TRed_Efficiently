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
// ONEITERATION, which together with the function MAIN constitutes the full
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
#include "tracek.h"
#include "oneiteration.h"
#include "parameters.h"
#include "NewMatrixL.h"

// MAX_COLUMNS_ROWS (Miriam Briskman, 02.23.2020)
#define MAX_COLUMNS_ROWS down_errors[left_errs].col + right_errors[right_errs].col \
                         > down_errors[left_errs].row + right_errors[right_errs].row \
                         ? down_errors[left_errs].col + right_errors[right_errs].col \
                         : down_errors[left_errs].row+right_errors[right_errs].row;
// COLUMN (Copied by Miriam Briskman, 02.23.2020, from 'errorsarray.h')
#define COLUMN ( i - abs(k - MAX_ERRORS) )

void OneIteration(char* str, long length, int ** matrix_forward, int ** matrix_reverse,
                  Entry* down_errors, Entry* right_errors, struct LastReportedRepeat left, 
                  struct LastReportedRepeat right, long int offset, int *upper, int *lower,
                  char *forward_pattern1, char *forward_pattern2, FILE* outfile, char *joinedstr[2],
                  suffix **suffixes[2], int **L, int *lcp[4], int *sizeToIndex[4], 
                  int *sizeToIndexFixed[4], int *smallest01[4][2], bool *isSuffixConsidered[4])
{

   /***************************** Variable Initialization *********************************/

   char *mid_forward = NULL, *mid_backward = NULL, *right_forward = NULL, *right_backward = NULL, 
        *left_forward = NULL, *left_backward = NULL;
   int i, j, k,
       begpos, endpos,
       n = length/2,
       curr_left, curr_right, left_errs, right_errs, curr_errs, curr_length, curr_rating, curr_off,
       comparison, diff, direction, longerstring,
       right_rating = 0, right_left = n, right_right = n-1, left_rating = 0, left_left = n, left_right = n-1,
       trunc_right, trunc_left,
       mid_forward_len = 0,  mid_backward_len = 0,
       right_forward_len = 0, right_backward_len = 0,
       left_forward_len = 0, left_backward_len = 0,
       rev_pattern_len = 0, 
       rev_text_len = 0, 
       mid_text_len = 0, 
       mid_pattern_len = 0,
       **temp = NULL,
       tempNum,
       k1,
       best_left, best_right, best_total,
       rowMinusOnes, // Miriam Briskman, 04.17.2020
       rowMinusOnesNeg, // Miriam Briskman, 04.17.2020
       rowMinusOnesLimit; // Miriam Briskman, 04.17.2020
   bool exitLoops = false, state;
   long pos1, pos2;
   struct LastReportedRepeat *prev = NULL;
   Entry * etemp = NULL;
   long long curr_text_len[4];
 
   /***************************** Setting Pointers *********************************/
    
   mid_forward = str + n;
   mid_forward_len = length - n;

   mid_backward = str + n - 1;
   mid_backward_len = n;

   right_forward = mid_forward + MIN_PERIOD - 1;
   right_forward_len = mid_forward_len - MIN_PERIOD + 1;

   left_forward = mid_forward - MIN_PERIOD + 1;
   left_forward_len = mid_forward_len + MIN_PERIOD - 1;

   right_backward = mid_backward + MIN_PERIOD - 1;
   right_backward_len = mid_backward_len + MIN_PERIOD - 1;

   left_backward = mid_backward - MIN_PERIOD + 1;
   left_backward_len = mid_backward_len - MIN_PERIOD + 1;

   left.left = n;
   left.right = n-1;
   right.left = n;
   right.right = n-1;

   /************************ PREPARE SUFFIX ARRAYS ***************************/

   curr_text_len[0] = 2*mid_forward_len + 2;

   prepareSuffixArrays (joinedstr[0], length + mid_forward_len + 1, str, length, mid_forward, 
                        mid_forward_len, curr_text_len[0], lcp[0], 
                        sizeToIndexFixed[0], sizeToIndex[0], suffixes[0], isSuffixConsidered[0], 
                        smallest01[0], 0);
 
   /***************************** OUTER LOOP *********************************/

   for (j = MIN_PERIOD; j <= n - 1 + length%2; j++)
   {
    for (direction = 0 ; direction < 2 ; direction++)
    {
    	if (j == n && direction == 0)
      	    continue;
        else if (direction == 0)
        {
            left_forward--;
            left_forward_len++;
            left_backward--;
            left_backward_len--;
            prev = &left;
            rev_text_len = left_backward_len;
            rev_pattern_len = mid_backward_len;
            mid_text_len = left_forward_len;
            mid_pattern_len = mid_forward_len;
            // FILLS UP MATRIX_FORWARD & MATRIX_REVERSE

            // Preparation to building the KxK matrix with the suffix array (complex.: O(n))
            // Miriam Briskman, 05.17.2020
            stringIncrement (joinedstr[0], length + mid_forward_len + 1, curr_text_len[0] + j, 
                             lcp[0], sizeToIndexFixed[0], sizeToIndex[0], suffixes[0], 
                             isSuffixConsidered[0], smallest01[0]);

            // Creating KxK with new algorithm (complex.: O(k^2))
            // Miriam Briskman, 05.17.2020
            buildMatrixL (matrix_forward, L, curr_text_len[0] + j, lcp[0], sizeToIndex[0], smallest01[0], 
                          mid_forward_len + 1 + j, mid_forward_len, MAX_ERRORS + j, 
                          MAX_ERRORS - mid_forward_len + 1);

            buildmatrixbackward (matrix_reverse, upper, lower, mid_backward, mid_backward_len, 
                                 left_backward, left_backward_len, -j);
        }
        else {
            right_forward++;
            right_forward_len--;
            right_backward++;
            right_backward_len++;
            prev = &right;
            rev_text_len = mid_backward_len;
            rev_pattern_len = right_backward_len;
            mid_text_len = mid_forward_len;
            mid_pattern_len = right_forward_len;
            // FILLS UP MATRIX_FORWARD & MATRIX_REVERSE
            buildmatrixforward (matrix_forward, upper, lower, right_forward, right_forward_len, 
                                mid_forward, mid_forward_len, j);

            buildmatrixbackward (matrix_reverse, upper, lower, right_backward, right_backward_len, 
                                 mid_backward, mid_backward_len, -j);
        }

        diff = j - prev->j;

        // Efficiently implementing the code from 'errorsarray.cpp' (Copied by Miriam Briskman, 02.20.2020):
        // FILLS UP RIGHT_ERRORS & DOWN_ERRORS
        // Miriam Briskman, 04.17.2020
        rowMinusOnes = (j <= mid_text_len) ? j - 1: mid_text_len;
        rowMinusOnesNeg = (j <= rev_pattern_len) ? j - 1: rev_pattern_len;

        // Fill errors arrays:
        for (i = 0; i <= MAX_ERRORS; i++)
        {
      	    // Zero out errors arrays:
            right_errors[i].row = 0;
            right_errors[i].col = 0;
            down_errors[i].row = 0;
            down_errors[i].col = 0;

            rowMinusOnesLimit = MAX_ERRORS + ((i <= rowMinusOnes) ? i : rowMinusOnes);
            for (k = MAX_ERRORS - ((i <= mid_pattern_len) ? i : mid_pattern_len); k <= rowMinusOnesLimit; k++)
      	    {
                tempNum = matrix_forward[k][COLUMN] + k - MAX_ERRORS;
                if (tempNum > right_errors[i].col)
                {
                    right_errors[i].row = matrix_forward[k][COLUMN];
                    right_errors[i].col = tempNum;
                }
            }
            // Special case bottom row of -1s applies:
            if (rowMinusOnes < i && i - 1 > right_errors[i].col)
            {
                right_errors[i].row = -1;
                right_errors[i].col = i - 1;
            }

            rowMinusOnesLimit = MAX_ERRORS - ((i <= rowMinusOnesNeg) ? i : rowMinusOnesNeg);
            for (k = MAX_ERRORS + ((i <= rev_text_len) ? i : rev_text_len); k >= rowMinusOnesLimit; k--)
                if (matrix_reverse[k][COLUMN] > down_errors[i].row)
                {
            	    down_errors[i].row = matrix_reverse[k][COLUMN];
                    down_errors[i].col = matrix_reverse[k][COLUMN] + k - MAX_ERRORS;
                }
        }

        // FIGURES OUT THE RANGE OF THE REPEATS IN THE CURRENT J
        curr_left = n;
        curr_right = n - 1;
        left_errs = 0;
        right_errs = 0;
        for (i = MAX_ERRORS; i >= 0; i--) 
        {
            k1 = MAX_ERRORS - i;
            if (i > 0 && right_errors[i - 1].row >= right_errors[i].row - 1)
                continue;
            while (k1 > 0 && down_errors[k1 - 1].col >= down_errors[k1].col - 1)
                k1--;
            if (right_errors[i].col + down_errors[k1].row >= j)
            {
                if (n + j*direction + right_errors[i].row - 1 > curr_right)
                {
                    curr_right = n + j*direction + right_errors[i].row - 1;
                    right_errs = i;
                }
                if (n + j*(direction - 1) - down_errors[k1].col < curr_left)
                {
                    curr_left = n + j*(direction - 1) - down_errors[k1].col;
                    left_errs = k1;
                }
            }
        }

        // IF NO REPEAT IS DETECTED, GO STRAIGHT TO CHECKING IF A PREVIOUS
        // REPEAT SHOULD BE REPORTED
        if (curr_left > curr_right)
        {
            comparison = -1;
            exitLoops = true;
        }

        if (!exitLoops)
        {
            // REMOVES POORLY-MATCHING PORTIONS (ACCORDING TO ERROR_VAL) FROM EACH END
            // OF THE POTENTIAL REPEAT
            for (i = right_errs - 1, trunc_right = right_errs; i >= 0; i--)
                if (right_errors[trunc_right].row - right_errors[i].row < (ERROR_VAL + 1)*(trunc_right - i))
                    trunc_right = i;
            for (i = left_errs - 1, trunc_left = left_errs; i >= 0; i--)
                if (down_errors[trunc_left].col - down_errors[i].col <= (ERROR_VAL + 1)*(trunc_left - i))
                    trunc_left = i;

            // IF THE REPEAT WAS REDUCED TO UNDER 2 FULL PERIODS BY THE PREVIOUS
            // FILTER, IT IS EXTENDED TO AT LEAST 2 FULL PERIODS WITH THE LEAST
            // POSSIBLE ERRORS
            if (right_errors[trunc_right].col + down_errors[trunc_left].row < j)
            {
                 for (i = trunc_right; i <= right_errs; i++)
             	     if (right_errors[i].col + down_errors[trunc_left].row >= j) 
                         break;
                 for (k = trunc_left; k <= left_errs; k++)
                     if (right_errors[i].col + down_errors[k].row >= j) 
                         break;
                 best_left = k;
                 best_right = i;
                 best_total = k+i;
                 for (i--; i > trunc_right; i--)
                 {
             	     while (k < left_errs && right_errors[i].col + down_errors[k].row < j)
                         k++;
                     if (right_errors[i].col + down_errors[k].row < j) 
                         break;
                     if (k+i < best_total)
                     {
                         best_left = k;
                         best_right = i;
                         best_total = k+i;
                     }
                 }
                 trunc_right = best_right;
                 trunc_left = best_left;
            }
            right_errs = trunc_right;
            left_errs = trunc_left;
            curr_right = n + j*direction + right_errors[right_errs].row - 1;
            curr_left = n + j*(direction - 1) - down_errors[left_errs].col;
            curr_off = n + j*direction - down_errors[left_errs].row;

            // GIVES THE REPEAT A RATING
            curr_errs = left_errs + right_errs;
            longerstring = MAX_COLUMNS_ROWS;
            curr_rating = longerstring - curr_errs * (ERROR_VAL + 1);
            curr_length = curr_right - curr_left + 1;

            // Comparisons: -1 is <, 0 is =, 1 is >
            if (curr_rating > prev->rating)
                comparison = 1;
            else if (prev->rating > curr_rating)
                comparison = -1;
            else
                comparison = 0;
        }
        exitLoops = false; // Reset condition for next iteration.
        
        // DETERMINES WHEN A (PREVIOUSLY SAVED) REPEAT SHOULD BE REPORTED AND
        // REPORTS IT
        if ( prev->newhigh == 1 && comparison == -1 )
        {
            prev->newhigh = 0;
            if ( (direction == 0 && (prev->rating > right_rating
                                     ||prev->left <= right_left - SHIFT
                                     ||prev->right >= right_right + SHIFT)) ||
                 (direction == 1 && (prev->rating > left_rating
                                     ||prev->left <= left_left - SHIFT
                                     ||prev->right >= left_right + SHIFT)) )
            {
              if (direction == 0)
              {
                  left_rating = prev->rating;
                  left_left = prev->left;
                  left_right = prev->right;
              }
              else
              {
                  right_rating = prev->rating;
                  right_left = prev->left;
                  right_right = prev->right;
              }

              begpos = prev->left;
              endpos = prev->right;

              fprintf (outfile, "%ld %ld %d\n", begpos + offset, 
                       endpos + offset, prev->right_errs + prev->left_errs);

              state = 0; // state 1 means the previous position had an insertion or deletion
              pos1 = (long)offset + begpos;
              pos2 = (long)offset + prev->offset;

                  // first, traceback the reverse matrix and print the 2-sequence alignment.
                  //    then, print repeat as an alignment. align1 begins at n/2
                  //    and align2 begins at n/2+j:
              if (direction) // Function calls placed in one conditional branch (Miriam Briskman, 03.03.2020)
              {
                    tracealign (prev->matrix_reverse, right_backward-diff, mid_backward,
                                prev->NoMinusOnesLowerRowBackward,
                                prev->down_errors[prev->left_errs].row-1,
                                prev->down_errors[prev->left_errs].col-1,
                                pos1, pos2, state, outfile);
                    traceforward (prev->matrix_forward, right_forward-diff, mid_forward,
                                  prev->NoMinusOnesUpperRowForward,
                                  prev->right_errors[prev->right_errs].row-1,
                                  prev->right_errors[prev->right_errs].col-1,
                                  pos1, pos2, state, forward_pattern1, forward_pattern2,
                                  outfile);
              }
              else
              {
                    tracealign (prev->matrix_reverse, mid_backward, left_backward+diff,
                                prev->NoMinusOnesLowerRowBackward,
                                prev->down_errors[prev->left_errs].row-1,
                                prev->down_errors[prev->left_errs].col-1,
                                pos1, pos2, state, outfile);
                    traceforward (prev->matrix_forward, mid_forward, left_forward+diff,
                                  prev->NoMinusOnesUpperRowForward,
                                  prev->right_errors[prev->right_errs].row-1,
                                  prev->right_errors[prev->right_errs].col-1,
                                  pos1, pos2, state, forward_pattern1, forward_pattern2,
                                  outfile);
              }

              // In case the repeat does not reach the midpoint n, we can take a
              // part of only one side

              // PSEUDOCODE

              // IF SIDE A > MINLENGTH AND SIDE B < MINLENGTH
              // FROM BEGINNING OF SIDE A CLOSEST TO CENTER TO END OF SIDE A
              // IF DELETING FROM END OF SIDE B TO CURRENT POINT WOULD RAISE
              // SCORE AND REPEAT WOULD REMAIN OVER MIN LENGTH, AND REPEAT WOULD REMAIN 2 FULL PERIODS, DO IT
            }
        }

        // SAVES CURRENT VALUES IN "PREV" IF THIS J IS A CANDIDATE TO BE REPORTED
        if (comparison > 0 && curr_length >= MIN_LENGTH && curr_rating >= MIN_RATING)
        {
            prev->newhigh        = 1;
            prev->left           = curr_left;
            prev->right          = curr_right;
            prev->right_errs     = right_errs;
            prev->left_errs      = left_errs;
            prev->length         = curr_length;
            prev->rating         = curr_rating;
            prev->j              = j;
            prev->offset         = curr_off;
            temp                 = matrix_forward;
            matrix_forward       = prev->matrix_forward;
            prev->matrix_forward = temp;
            temp                 = matrix_reverse;
            matrix_reverse       = prev->matrix_reverse;
            prev->matrix_reverse = temp;
            temp                 = NULL;
            etemp                = down_errors;
            down_errors          = prev->down_errors;
            prev->down_errors    = etemp;
            etemp                = right_errors;
            right_errors         = prev->right_errors;
            prev->right_errors   = etemp;
            etemp                = NULL;
            // Miriam Briskman, 04.17.2020:
            prev->NoMinusOnesUpperRowForward = MAX_ERRORS + rowMinusOnes;
            prev->NoMinusOnesLowerRowBackward = MAX_ERRORS - rowMinusOnesNeg;
        }  // endif
    } // end for one direction
   } // end for j

   /***************************** OUTER LOOP END *****************************/
   return;
}

/*
                          NewMatrixL.cpp

  The new L (arrow shaped) matrix is of the form (MAX_ERRORS = 20):

 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -2 20
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -2 19  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -2 18  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -2 17  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -2 16  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -2 15  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -2 14  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0 -2 13  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0 -2 12  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0 -2 11  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0 -2 10  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0 -2  9  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0 -2  8  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0 -2  7  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0 -2  6  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
 0  0  0  0  0 -2  4  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5
 0  0  0  0 -2  3  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5
 0  0  0 -2  2  3  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5
 0  0 -2  1  2  3  4  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5
 0 -2  0  1  2  3  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5
-2 -1  0  1  2  4  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5
 0 -2 -1  0  1  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5
 0  0 -2 -1  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5
 0  0  0 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
 0  0  0  0 -2 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0 -2 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0 -2 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0 -2 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0 -2 -1  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0 -2 -1  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0 -2 -1  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0 -2 -1  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0 -2 -1  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0 -2 -1  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0 -2 -1  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -2 -1  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -2 -1  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -2 -1  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -2 -1  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -2 -1  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -2 -1  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -2 -1

The elements to the right of the arrow and between the -1 rows are
the elements that are being copied into the KxK matrix.

Function Description:

prepareSuffixArrays: 1) Creates the suffix array (sort the suffixes) 
                        based on the concatenated text and pattern.
                     2) Creates the longest common prefix (lcp) array. 
                     3) Creates the 'sizeToIndexFixed' array, which is
                        a map from the sizes of suffixes to their
                        position in the array.
                     4) Creates the 'sizeToIndex' array, which is a map
                        from the sizes of suffixes of current size to 
                        their position in the array.
                     5) Create the 'isSuffixConsidered' array, in which
                        each element is 0 is the suffix at that index is
                        longer than the group of suffixes we consider,
                        or 1 if it is included in the group of suffixes
                        that we consider.
                     6) Populates the 'smallest01' arrays.

stringIncrement: 1) Updates the 'isSuffixConsidered' array by placing '1'
                    in the slot correponding to the suffix that is added
                    in the current increment.
                 2) Updates the 'lcp' array.
                 3) Updates the 'sizeToIndex' array.
                 4) Updates the 'smallest01' arrays.

buildMatrixL: 1) Initializes L matrix based on the region of the KxK
                 matrix that we are filling.
              1) Queries the suffix and lcp array.
              2) Populates the L matrix, and copies corresponding elements
                 to the KxK matrix.

printSuffixes: Prints the suffix array at its current state.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// Header files within the same folder:
#include "NewMatrixL.h"
#include "parameters.h"


// Creating a new suffix array:
void prepareSuffixArrays (char *str, const long long str_len, const char* text, const int text_len, 
                          const char* pattern, const int pattern_len, const int curr_text_len,
                          int *lcp, int *sizeToIndexFixed, int *sizeToIndex, suffix **suffixes, 
                          bool *isSuffixConsidered, int *smallest01[], bool reversed)
{
        suffix *tempSuf, *highSuf;
        long long min, compareValue;
        int i, j, k, p, temp, count,
            l = 0, h = str_len - 1,
        // Create an auxiliary stack for the iterative quicksort of suffixes:
            *stack = (int*) malloc ((h - l + 1) * sizeof(int)),
        // initialize top of stack 
            top = -1;
        // Set each of the indices to be a value outside of the indices of the suffixes array!
        // Index 0 is for the array at smallest01[0], etc.
        long long tempIndex[3] = {str_len, str_len, str_len};

        // Create a new char array, called str, that has: <text>#<pattern>\0.
        if (reversed)
        {
            for (i = 0; i < text_len; i++)
               str[i] = text[-i];
            str[text_len] = '#';
            for (i = 0; i < pattern_len; i++)
               str[text_len + 1 + i] = pattern[-i];
        }
        else
        {
            for (i = 0; i < text_len; i++)
               str[i] = text[i];
            str[text_len] = '#';
            for (i = 0; i < pattern_len; i++)
               str[text_len + 1 + i] = pattern[i];
        }
        str[text_len + pattern_len + 1] = '\0';
           
        
        // Initializing the suffix arrays:
        for (i = 0; i < str_len; i++)
        {
            suffixes[i][0].index = i;
            suffixes[i][0].text_len = str_len - i;
        }

        // Quicksort to sort the suffixes array:
        // The code is based on: https://www.geeksforgeeks.org/iterative-quick-sort/

        // push initial values of l and h to stack 
        stack[++top] = l; 
        stack[++top] = h;
  
        // Keep popping from stack while is not empty 
        while (top >= 0) 
        { 
           // Pop h and l 
           h = stack[top--]; 
           l = stack[top--]; 
  
           // Set pivot element at its correct position 
           // in sorted array 
           // The following is from the 'partition' method:
           /**/ highSuf = suffixes[h];
           /**/ i = l - 1; 
           /**/
           /**/ for (j = l; j <= h - 1; j++)
           /**/ {
           /**/    // We want to see if highSuf is 'larger than' (*suffixes)[j].
           /**/    // The code from 'compareTo':
           /**/    min = (highSuf[0].text_len < suffixes[j][0].text_len) 
           /**/            ? highSuf[0].text_len : suffixes[j][0].text_len;
           /**/    k = 0;
           /**/    while (k < min && str[highSuf[0].index + k] == str[suffixes[j][0].index + k])
           /**/        k++;
           /**/    if (k < min)
           /**/         if (str[highSuf[0].index + k] < str[suffixes[j][0].index + k])
           /**/             compareValue = -1;
           /**/         else
           /**/             compareValue = 1;
           /**/    else
           /**/        compareValue = highSuf[0].text_len - suffixes[j][0].text_len;
           /**/
           /**/    if (compareValue > 0) // If 'larger than':
           /**/    { 
           /**/        i++;
           /**/        // Swap:
           /**/        tempSuf = suffixes[i];
           /**/        suffixes[i] = suffixes[j];
           /**/        suffixes[j] = tempSuf;
           /**/    }
           /**/ }
           /**/
           /**/ // Swap:
           /**/ tempSuf = suffixes[i + 1];
           /**/ suffixes[i + 1] = suffixes[h];
           /**/ suffixes[h] = tempSuf;

           p = i + 1;
  
           // If there are elements on left side of pivot, 
           // then push left side to stack 
           if (p - 1 > l) 
           { 
               stack[++top] = l; 
               stack[++top] = p - 1; 
           } 
  
           // If there are elements on right side of pivot, 
           // then push right side to stack 
           if (p + 1 < h) 
           { 
               stack[++top] = p + 1; 
               stack[++top] = h; 
           } 
        } // End of Quicksort

        free (stack);

        j = 0;
        // Initializing the 'sizeToIndex' array:
        for (i = 0; i < str_len; i++)
        {
           temp = suffixes[i][0].text_len;
           sizeToIndexFixed[temp] = i;
           if (temp <= curr_text_len)
           {
              isSuffixConsidered[i] = true;
              sizeToIndex[temp] = j++;
           }
           else
             isSuffixConsidered[i] = false;
        }

        j = str_len - 1;
        // Finding a relevant suffix from the end:
        while (!isSuffixConsidered[j])
          j--;

        // Initializing the 'lcp' array:
        for (count = curr_text_len - 1; count >= 1; count--)
        {
            i = j;
            j--;
            while (!isSuffixConsidered[j] && j >= 0)
               j--;
            if (j < -1)
              break;
            
            min = (suffixes[i][0].text_len < suffixes[j][0].text_len) 
                  ? suffixes[i][0].text_len : suffixes[j][0].text_len;
            // Optimization: we have 'pattern_len' instances of the followng form:
            // suffixes[i - 1][0] = 123
            // suffixes[i][0] = 123#89123
            // The lcp value is len(123) = 3, which is 'min'.
            // Hence, instead of letting the 'while' loop run for 'min' iterations,
            //   we simply assign 'min' to the corresponding 'lcp' slot.
            if ((suffixes[i][0].text_len - suffixes[j][0].text_len == pattern_len + 1) &&
                suffixes[j][0].text_len <= pattern_len)
                  lcp[count] = min;
            else
            {
                k = 0;
                while (k < min && str[suffixes[i][0].index + k] == str[suffixes[j][0].index + k])
                    k++;
                lcp[count] = k;
            }
            if (lcp[count] < 2)
                tempIndex[lcp[count]] = count;
            smallest01[0][count] = tempIndex[0];
            smallest01[1][count] = tempIndex[1];
    	}
} // End of 'prepareSuffixArrays'

void stringIncrement (char *str, const long long str_len, const int curr_text_len, int *lcp, 
                      int *sizeToIndexFixed, int *sizeToIndex, suffix **suffixes, bool *isSuffixConsidered, 
                      int *smallest01[])
{
    int i, j, k, l, n, m, newIndex, leftSufIndex, min, prevLCP;

    newIndex = sizeToIndexFixed[curr_text_len];
    isSuffixConsidered[newIndex] = true;
    
    // We break the update of the lcp array into 2 cases:
    //   1. The newly added suffix (of length curr_text_len) is at the very bottom of the suffix array.
    //   2. The newly added suffix is anywhere else in the array.

    i = newIndex - 1;
    while (!isSuffixConsidered[i])
       i--;
    // Now i points to an already existing suffix in the scope of lcp.
    // Next, let's extract the index of the suffix of lesser size to its left:
    leftSufIndex = sizeToIndex[suffixes[i][0].text_len];
    // We check if we fall into case 1 above:
    if (leftSufIndex == curr_text_len - 2)
    {
        // Update the sizeToIndex array:
        sizeToIndex[curr_text_len] = leftSufIndex + 1;

        // Update the lcp array:
        min = suffixes[i][0].text_len;
        k = 0;
        while (k < min && str[suffixes[i][0].index + k] == str[suffixes[newIndex][0].index + k])
           k++;
        lcp[curr_text_len - 1] = k;

        // Update the smallest01 arrays:
        smallest01[0][curr_text_len - 1] = str_len;
        smallest01[1][curr_text_len - 1] = str_len;
        if (k < 2)
        {
           i = curr_text_len - 1;
           smallest01[k][i] = i;
           j = i - 1;
           while (j > 0 && lcp[j] != k)
             smallest01[k][j--] = i;
        }
    }
    else
    {
        // Update the sizeToIndex array:
        //   1. Increase by 1 the indices of the elements to the right of the insertion:
        l = newIndex + 1;
        while (!isSuffixConsidered[l])
               l++;
        m = l;
        prevLCP = lcp[leftSufIndex + 1];
        while (l < str_len)
        {
            while (!isSuffixConsidered[l] && l < str_len)
               l++;
            if (l < str_len)
            {
               sizeToIndex[suffixes[l][0].text_len]++;
               l++;
            }
        }

        l = m;

        //   2. Update the index of the newly added suffix:
        sizeToIndex[curr_text_len] = leftSufIndex + 1;

        // Update the lcp array:
        //   1. Move all elements to the right of the new one 1 place to the right:
        for (k = curr_text_len - 1; k > leftSufIndex + 2; k--)
           lcp[k] = lcp[k - 1];

        //   2. Compare the new suffix with the suffix to the left (whose index is i):
        min = suffixes[i][0].text_len;
        k = 0;
        while (k < min && str[suffixes[i][0].index + k] == str[suffixes[newIndex][0].index + k])
           k++;
        lcp[leftSufIndex + 1] = k;

        //   3. Compare the new suffix with the suffix to the right (whose index is l):
        n = 0;
        while (n < min && str[suffixes[l][0].index + n] == str[suffixes[newIndex][0].index + n])
           n++;
        lcp[leftSufIndex + 2] = n;

        // Update the smallest01 arrays:
        for (m = 0; m < 2; m++)
        {
            i = str_len;
            for (j = curr_text_len - 1; j > 0 && smallest01[m][j] != j - 1; j--)
            {
               if (lcp[j] == m)
                  i = j;
               smallest01[m][j] = i;
            }
            if (smallest01[m][j] == j - 1)
            {
               if (lcp[j] == m)
                  i = j;
               smallest01[m][j--] = i;
               while (j > 0 && lcp[j] != m)
                  smallest01[m][j--] = i;
            }
        }
    }
}

//MatrixL constructor
void buildMatrixL (int **matrix, int **L, const long long str_len, int *lcp, int *sizeToIndex, 
                   int *smallest01[], const int text_len, const int pattern_len, int upLim, int lowLim)
{
    // Variable initialization:
    int rows = 2*MAX_ERRORS + 3,
        cols = MAX_ERRORS + 3,
        i, d, e, temp, temp2, row, rank1, rank2;

    // Initialization of L matrix is required only once, when
    //    the memory for the matrix is allocated, for the loop
    //    below (1) does not touch the 'border' elements that
    //    were initialized, and (2) overwrites all the other
    //    elements anyway.

    // We do need, however, to set the upper and lower boundaries 
    //    of the L matrix to -1s (this is O(k), linear in the
    //    number of errors) to demarcate the region in the matrix 
    //    where element are to be computed.
    // Set lower boundary:
    for (i = MAX_ERRORS + 3 - lowLim; i < MAX_ERRORS + 3; i++)
       L[lowLim - 1][i] = -1;
    // Set upper boundary:
    for (i = upLim - MAX_ERRORS + 2; i < MAX_ERRORS + 3; i++)
       L[upLim + 1][i] = -1;

    // Main Nested Loop:
    //iterate over each column:
    for (e = 0; e <= MAX_ERRORS; e++) 
    {
        temp2 = (e + 1 + MAX_ERRORS < upLim) ? e + 1 + MAX_ERRORS : upLim;
        //iterate over select rows:
        for (d = ((-e + 1 + MAX_ERRORS > lowLim) ? -e + 1 + MAX_ERRORS : lowLim) ; d <= temp2; d++) 
        {
             row = (L[d + 1][e + 1] + 1 > L[d][e + 1] + 1) 
                   ? L[d + 1][e + 1] + 1 : L[d][e + 1] + 1;
             row = (row > L[d - 1][e + 1]) ? row : L[d - 1][e + 1];

             row = (row < pattern_len) ? row : pattern_len;

             // Below is the code from the original calculateLCP function.
             // Author: Robin Cohen Date: 11/9/2014

             //if either query is beyond the bounds of its respective string within str, return a 0
         	 if (row + d - MAX_ERRORS >= text_len || row >= pattern_len)
                  L[d][e + 2] = row;
         	 else
             {
                 // compute the rank in the lcp array of the substrings starting at query1 and query2
                 // New! We made the computation constant. Miriam Briskman, 05.03.2020.
                 rank1 = sizeToIndex[str_len - row - d + MAX_ERRORS];
                 rank2 = sizeToIndex[str_len - row - text_len - 1];

             	 //the lower rank is not the lower bound in the range minimum query of the LCP array;
             	 //rather the lower rank+1 is the lower bound
                 if (rank1 < rank2)
                 	rank1++;
                 else
                 {
                     // Exchanging the ranks for easier work afterwards:
                     temp = rank1;
                     rank1 = rank2 + 1;
                     rank2 = temp;
                 }

                 // Most of the cases fall into the first 2 cases of the if-else-if
                 //    branch below, which executes at constant time.
                 // Miriam Briskman, 05.06.2020.
                 if (smallest01[0][rank1] <= rank2)
                    L[d][e + 2] = row;
                 else if (smallest01[1][rank1] <= rank2)
                    L[d][e + 2] = row + 1;
                 else
                 {
                     // The very small number of cases that does not fall into
                     //   the branch above requires a sequential approach. However,
                     //   since all of 0s and 1s are very common in the lcp
                     //   array, only very short subsequences of lcp have non of 
                     //   0 or 1. The loop should execute for about lg(n) time
                     //   (less than 10 elements on average.) Since these cases
                     //   are rare, we preserve an expected constant work in
                     //   'calculateLCP'.
                     // Miriam Briskman, 05.06.2020.
                     temp = rank1;
                     for (i = rank1 + 1; i <= rank2; ++i)
                     {
                         if (lcp[i] < lcp[temp])
                             temp = i;
                         if (lcp[temp] == 2)
                             break;
                     }
                     L[d][e + 2] = row + lcp[temp];
                 }
             }
 
             // Assign the computed element into the corresponding location in the KxK matrix:
             matrix[d - 1][e - ((MAX_ERRORS >= d - 1) ? MAX_ERRORS - d + 1 : d - MAX_ERRORS - 1)] = L[d][e + 2];

        }//end for

    } //end for

    // end of filling L matrix
}

void printSuffixes (char *str, const long long str_len, suffix **suffixes)
{
    printf ("Suffix Array:\n");
    for (int i = 0; i < str_len; i++)
       printf ("%s\n", (str + suffixes[i][0].index));
}

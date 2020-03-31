/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\

* A new algorithm for using Edit Distance matrices of the form AxB to compute Edit  *
     Distance matrices of the form aAxB, where 'A' and 'B' are strings, 'B' is a 
*    suffix of 'A', and 'a' is a character. This algorithm reduces the constant of  *
     the Edit Distance running time by at least a quarter. That is, given that the 
*    running time of the conventional algorithm for computing an Edit Distance      *
     matrix from scratch uses about 7xAB computations, the new algorithm uses about 
*    5.5xAB computations.                                                           *
  Written by Miriam Briskman, for Brooklyn College CISC 5001 course of Spring 2020.
* Supervised by Professor Dina Sokol.                                               *

\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>       // for 'time'
#include <math.h>       // for 'floor'

#include <iostream>

using namespace std;

// Global Constants:
const int DOUBLE_PERIOD = 500,
          TRIALS = 10000;

// Function Prototypes:
void displayMatrix (char const *argv[]);
void expendNbyNMatrix (char const *);
bool testingEditDistance (int **, char const *, char const *, int);
bool algo (char const *argv[]);

// Global Variables:
int text_len,
    **matrix, // pointer to int matrix (double array)
    tempVal, // Multi-purpose value
    i, j, k, // Loop Indices
    whatRow,
    half_text,
    pure_half;
FILE *outfile; // Output file pointer

int main (int argc, char *argv[])
{
   // Check for correct number of arguments:
   if (argc != 2)
   {
      cout << "Program requires 1 command line arguments: the name of output file!" << endl;
      exit (EXIT_FAILURE);
   }

   // Open the output file.
   if ((outfile = fopen(argv[1], "w")) == NULL)
   { 
       cout << "Error opening output file." << endl;
       exit (EXIT_FAILURE);
   }

   // Allocate memory for a c-string of size DOUBLE_PERIOD:
   char * ranString = (char*) malloc ((DOUBLE_PERIOD + 2) * sizeof (char)); 
   char const * arr[2];
   int e;
   short temp;

   srand (time(NULL));

   arr[0] = NULL;
   arr[1] = ranString;

   printf ("Testing began ... Please Wait ...\n");
   
   for (int l = 0; l < TRIALS; l++)
   {
       text_len = rand() % DOUBLE_PERIOD + 2;

       // Create the random string:
       for (e = 0; e < text_len; e++)
       {
           temp = floor(rand()/(double)RAND_MAX * 4);
           switch (temp)
           {
              case 0: ranString[e] = 'A';
                      break;
              case 1: ranString[e] = 'C';
                      break;
              case 2: ranString[e] = 'G';
                      break;
              case 3: ranString[e] = 'T';
                      break;
              default: ranString[e] = 'N';
                       break;
           }
       }

       ranString[e] = '\0';
       
       if (!algo (arr))
       {
          // Cleaning:
          free (ranString); 
          fclose (outfile);
          return EXIT_FAILURE;
       }
   }

   printf ("Testing Successful!\n");

   // Cleaning:
   free (ranString); 
   fclose (outfile);

   return EXIT_SUCCESS;
}

bool algo (char const *argv[])
{
   // Find the size of the string:
   text_len = strlen (argv[1]);
   half_text = text_len/2 + text_len%2;
   pure_half = text_len/2;

   // Allocate memory for the matrix:
   matrix = (int**) malloc ((text_len + 1) * sizeof (int*));
   if (matrix == NULL)
   {
      cout << "Memory allocation issues in main!" << endl;
      return 0;
   }
   for (i = 0; i < text_len + 1; i++)
   {
      matrix[i] = (int*) malloc ((half_text + 1) * sizeof (int));
      if (matrix[i] == NULL)
      {
         cout << "Memory allocation issues in main!" << endl;
         return 0;
      }
   }

   whatRow = pure_half;

   // Compute the symmetric edit distance. 
   // By symmetric edit distance, we mean the Edit Distance matrix where both the 
   //   'pattern' on the top and the 'text' on the left are the second half
   //   of the text. For example, if the text is "CGACCTGAGC", the corresponding
   //   symmetric Edit Distance matrix is:
   //                 T  G  A  G  C
   //              0  1  2  3  4  5
   //           T  1  0  1  2  3  4
   //           G  2  1  0  1  2  3
   //           A  3  2  1  0  1  2
   //           G  4  3  2  1  0  1
   //           C  5  4  3  2  1  0
   // which is clearly symmetric across its main diagonal.
   // Here, each value is given
   //   by the formula:
   //       matrix[i][j] = j - i + whatRow
   // where i begins as whatRow , and j begins at i - whatRow:
   for (i = whatRow; i < text_len + 1; i++)
       for (j = i - whatRow; j < half_text + 1; j++)
       {
           matrix[i][j] = j - i + whatRow;
           matrix[j + whatRow][i - whatRow] = matrix[i][j];
       }

   // Print the symmetric matrix:
   //displayMatrix (argv);

   // Expand the symmetric matrix by adding a prefix to the text on the left
   //   from the rest of the characters at a time. For example, the first expansion
   //   of the symmetric matrix from the previous example of "CGACCTGAGC" is:
   //           T  G  A  G  C
   //        0  1  2  3  4  5
   //     C  1  1  2  3  4  4
   //     T  2  1  2  3  4  5
   //     G  3  2  1  2  3  4
   //     A  4  3  2  1  2  3
   //     G  5  4  3  2  1  2
   //     C  6  5  4  3  2  1
   //   and the next expansion is:
   //           T  G  A  G  C
   //        0  1  2  3  4  5
   //     C  1  1  2  3  4  4
   //     C  2  2  2  3  4  4
   //     T  3  2  3  3  4  5
   //     G  4  3  2  3  3  4
   //     A  5  4  3  2  3  4
   //     G  6  5  4  3  2  3
   //     C  7  6  5  4  3  2

   for (k = 0; k < pure_half; k++)
   {
       // Expand the matrix once:
       expendNbyNMatrix (argv[1]);

       // Print the expended matrix:
       //displayMatrix (argv);

       if (!testingEditDistance (matrix, argv[1] + whatRow, argv[1] + pure_half, whatRow))
          return 0;
   }

   for (i = 0; i < text_len + 1; i++)
      free (matrix[i]);
   free (matrix);

   return 1;
}


/* displayMatrix
   Prints a matrix of the form
              T  G  A  G  C
           0  1  2  3  4  5
        C  1  1  2  3  4  4
        T  2  1  2  3  4  5
        G  3  2  1  2  3  4
        A  4  3  2  1  2  3
        G  5  4  3  2  1  2
        C  6  5  4  3  2  1
   to the 'output' file.
   Referring to this example, the first
   loop prints
              T  G  A  G  C
   The second loop prints
           0  1  2  3  4  5
   And the 3rd nested loop prints
        C  1  1  2  3  4  4
        T  2  1  2  3  4  5
        G  3  2  1  2  3  4
        A  4  3  2  1  2  3
        G  5  4  3  2  1  2
        C  6  5  4  3  2  1
*/
void displayMatrix (char const* argv[])
{
   // Printing the pattern to the file:
   fprintf (outfile, "      ");
   for (int i = 0; i < pure_half + 1; i++)
      fprintf (outfile, "%3c", argv[1][i + pure_half]);
   fprintf (outfile, "\n");
   
   // Printing the first row of the matrix:
   fprintf (outfile, "   ");
   for (int i = 0; i < half_text + 1; i++)
      fprintf (outfile, "%3d", matrix[whatRow][i]);
   fprintf (outfile, "\n");

   // Printing the rest of the matrix:
   for (int i = whatRow + 1; i < text_len + 1; i++)
   {
      fprintf (outfile, "%3c", argv[1][i-1]);
      for (int j = 0; j < half_text + 1; j++)
         fprintf (outfile, "%3d", matrix[i][j]);
      fprintf (outfile, "\n");
   }
   fprintf (outfile, "\n");
}

/* expendNbyNMatrix
   Computes the Edit Distance matrix for aAxB from an existing
   Edit Distance matrix for AxB.
*/
void expendNbyNMatrix (char const *theText)
{
   // Compute the other matrices by using the incremental string 
   //   comparison method.
   // Increase every element in the current matrix by 1:
   // Here we use 2xAB computations (one for the incrementing 
   //    index j and one for the incrementing matrix element:)
   for (i = whatRow; i < text_len + 1; ++i)
       for (j = 0; j < half_text + 1; ++j)
           ++matrix[i][j];

   --whatRow;
 
   // Populate the new 0th row with 0..half_text:
   for (j = 0; j < half_text + 1; ++j)
      matrix[whatRow][j] = j;

   // For the rest of the columns, find the top element, and then
   //   check if the elements below are the same as the ones
   //   currently in place:
   for (j = 0; j < half_text; j++) // Loop controlling the columns
   {
       // Check if the element on a row below on the same column is the same
       //   as the current one. If yes, go to the next column. If not, keep
       //   inceasing / decreasing until the same value is met
       if (theText[whatRow] != theText[pure_half + j])
           tempVal = matrix[whatRow][j] + 1;
       else
           tempVal = matrix[whatRow][j];

       if (tempVal > matrix[whatRow + 1][j] + 1)
           tempVal = matrix[whatRow + 1][j] + 1;

       i = whatRow;

       do
       {
          // Here we use 3.5xAB computations (7 steps done for
          //    about 1/2 of the matrix:)
          matrix[i + 1][j + 1] = tempVal;
          ++i;

          if (theText[i] == theText[pure_half + j])
             tempVal = matrix[i][j]; // On the same diagonal, up-left place.
          else
             tempVal = matrix[i][j] + 1;

          if (tempVal > matrix[i][j + 1] + 1) // Element above.
             tempVal = matrix[i][j + 1] + 1;
          if (tempVal > matrix[i + 1][j] + 1) // Element to the left.
             tempVal = matrix[i + 1][j] + 1;

       } while (matrix[i + 1][j + 1] != tempVal);

   }
}

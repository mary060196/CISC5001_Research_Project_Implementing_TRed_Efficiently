/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\

* A new algorithm for using Edit Distance matrices of the form AxbB to compute an
     Edit Distance matrice of the form AxB, where 'A' and 'B' are strings, 'B' is a 
*    suffix of 'A', and 'b' is a character. This algorithm reduces the constant of 
     the Edit Distance running time by at least a quarter. That is, given that the 
*    running time of the conventional algorithm for computing an Edit Distance 
     matrix from scratch uses about 7xAB computations, the new algorithm uses about 
*    4xAB computations.
  Written by Miriam Briskman, for Brooklyn College CISC 5001 course of Spring 2020.
* Supervised by Professor Dina Sokol.

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
void shrinkNbyNMatrix (char const *);
bool testingEditDistance (int **, char const *, char const *, int);
bool algo (char const *argv[]);

// Global Variables:
int text_len,
    **matrix, // pointer to int matrix (double array)
    tempVal, // Multi-purpose value
    i, j, k, // Loop Indices
    whatCol;
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

   // Allocate memory for a square matrix:
   matrix = (int**) malloc ((text_len + 1) * sizeof (int*));
   if (matrix == NULL)
   {
      cout << "Memory allocation issues in main!" << endl;
      return 0;
   }
   for (i = 0; i < text_len + 1; i++)
   {
      matrix[i] = (int*) malloc ((text_len + 1) * sizeof (int));
      if (matrix[i] == NULL)
      {
         cout << "Memory allocation issues in main!" << endl;
         return 0;
      }
   }

   whatCol = 0;

   // Compute the symmetric edit distance. 
   // By symmetric edit distance, we mean the Edit Distance matrix where both the 
   //   'pattern' on the top and the 'text' on the left are the second half
   //   of the text. For example, if the text is "TGAGC", the corresponding
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
   //       matrix[i][j] = j - i

   for (i = 0; i < text_len + 1; i++)
       for (j = i; j < text_len + 1; j++)
           matrix[i][j] = matrix[j][i] = j - i;

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

   for (k = 1; k < text_len; k++)
   {
       // Expand the matrix once:
       shrinkNbyNMatrix (argv[1]);

       // Print the expended matrix:
       //displayMatrix (argv);

       if (!testingEditDistance (matrix, argv[1], argv[1] + whatCol, whatCol))
       {
          for (i = 0; i < text_len + 1; i++)
             free (matrix[i]);
          free (matrix);
          fclose (outfile);
          return 0;
       }
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
   for (int i = whatCol; i < text_len; i++)
      fprintf (outfile, "%3c", argv[1][i]);
   fprintf (outfile, "\n");
   
   // Printing the first row of the matrix:
   fprintf (outfile, "   ");
   for (int i = whatCol; i < text_len + 1; i++)
      fprintf (outfile, "%3d", matrix[0][i]);
   fprintf (outfile, "\n");

   // Printing the rest of the matrix:
   for (int i = 1; i < text_len + 1; i++)
   {
      fprintf (outfile, "%3c", argv[1][i-1]);
      for (int j = whatCol; j < text_len + 1; j++)
         fprintf (outfile, "%3d", matrix[i][j]);
      fprintf (outfile, "\n");
   }
   fprintf (outfile, "\n");
}

/* shrinkNbyNMatrix
   Computes the Edit Distance matrix for AxB from an existing
   Edit Distance matrix for AxbB.
*/
void shrinkNbyNMatrix (char const *theText)
{
   // Copying elements under the main digonal 1 place to the
   //   right.
   // Here we use 1xAB computations (1 'i' decrement and 1 copy,
   //    done on about 1/2 of the matrix, which is (1+1)/2 = 
   //    1xAB computations.)
   for (j = text_len; j > whatCol - 1; --j)
       for (i = text_len; i > j; --i)
           matrix[i][j] = matrix[i][j - 1];

   // Initializing the leftmost column (Linear running time):
   for (i = 1; i < text_len + 1; ++i)
       matrix[i][whatCol + 1] = i;   

   // Initializing the top row (Linear running time):
   for (j = text_len; j > whatCol; --j)
       matrix[0][j] = matrix[0][j - 1];

   ++whatCol;

   // Here, we use 3xAB computations (1 'i' increment, 2 for the 'elseif'
   //    statement, and 1.5*2 = 3 steps for both 'if' statments, each of which
   //    is done on about 1/2 of the matrix, making it (1+2+3)/2 = 3xAB 
   //    computations.)
   for (j = whatCol + 1; j < text_len + 1; j++) // Loop controlling the columns
   {
        for (i = 1; i <= j; i++)
        {
          if (theText[i - 1] == theText[j - 1])
             matrix[i][j] = matrix[i - 1][j - 1]; // On the same diagonal, up-left place.
          else
             matrix[i][j] = matrix[i - 1][j - 1] + 1;

          if (matrix[i][j] > matrix[i - 1][j] + 1) // Element above.
             matrix[i][j] = matrix[i - 1][j] + 1;
          if (matrix[i][j] > matrix[i][j - 1] + 1) // Element to the left.
             matrix[i][j] = matrix[i][j - 1] + 1;
        }
   }
}

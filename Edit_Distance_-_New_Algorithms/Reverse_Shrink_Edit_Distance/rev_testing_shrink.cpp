#include <stdio.h>
#include <stdlib.h>
#include <string.h>

bool testingEditDistance (int **matrixToTest, char const *text, int text_len, 
                                              char const *pattern, int pattern_len, int whatCol)
{
   // First, create the Edit Distance matrix using the
   //   straightforward algorithm:
   int **matrix, // pointer to int matrix (double array)
       tempVal, // Multi-purpose value
       i, j; // Loop Indices

   // Allocate memory for matrix:
   matrix = (int**) malloc ((text_len + 1) * sizeof (int*));
   if (matrix == NULL)
   {
      printf ("Memory allocation issues in main!");
      exit (EXIT_FAILURE);
   }
   for (i = 0; i < text_len + 1; i++)
   {
      matrix[i] = (int*) malloc ((pattern_len + 1) * sizeof (int));
      if (matrix[i] == NULL)
      {
         printf ("Memory allocation issues in main!");
         exit (EXIT_FAILURE);
      }
   }
   
   // Compute Edit Distance matrix:
   // Initialization:
   // Left (1st) Column:
   for (i = 0; i < text_len + 1; i++)
      matrix[i][0] = i;
   // Top (1st) Row:
   for (i = 0; i < pattern_len + 1; i++)
      matrix[0][i] = i;

   // Filling the Rest of the Matrix:
   for (i = 1; i < text_len + 1; i++)
       for (j = 1; j < pattern_len + 1; j++)
       {
          if (text[-i+1] == pattern[-j+1])
             tempVal = matrix[i - 1][j - 1]; // On the same diagonal, up-left place.
          else
             tempVal = matrix[i - 1][j - 1] + 1;
          if (tempVal > matrix[i - 1][j] + 1) // Element above.
             tempVal = matrix[i - 1][j] + 1;
          if (tempVal > matrix[i][j - 1] + 1) // Element to the left.
             tempVal = matrix[i][j - 1] + 1;
          matrix[i][j] = tempVal;
       }

   // Second, compare the matrices:
   for (i = 1; i < text_len + 1; i++)
       for (j = 1; j < pattern_len + 1; j++)
       {
           if (matrix[i][j] != matrixToTest[i][j + whatCol])
           {
               printf ("The matrices are not equal! The CORRECT matrix is:\n");
               // Printing the pattern to the file:
               printf ("      ");
               for (i = 0; i < pattern_len; i++)
                  printf ("%3c", pattern[-i]);
               printf ("\n");
   
               // Printing the first row of the matrix:
               printf ("   ");
               for (i = 0; i < pattern_len + 1; i++)
                  printf ("%3d", matrix[0][i]);
               printf ("\n");

               // Printing the rest of the matrix:
               for (i = 0; i < text_len; i++)
               {
                  printf ("%3c", text[-i]);
                  for (j = 0; j < pattern_len + 1; j++)
                     printf ("%3d", matrix[i+1][j]);
                  printf ("\n");
               }

               printf ("\nAnd the WRONG matrix is:\n");
               // Printing the pattern to the file:
               printf ("      ");
               for (i = 0; i < pattern_len; i++)
                  printf ("%3c", pattern[-i]);
               printf ("\n");
   
               // Printing the first row of the matrix:
               printf ("   ");
               for (i = 0; i < pattern_len + 1; i++)
                  printf ("%3d", matrixToTest[0][i + whatCol]);
               printf ("\n");

               // Printing the rest of the matrix:
               for (i = 0; i < text_len; i++)
               {
                  printf ("%3c", text[-i]);
                  for (j = 0; j < pattern_len + 1; j++)
                     printf ("%3d", matrixToTest[i + 1][j + whatCol]);
                  printf ("\n");
               }

               return false;
           }
       }

   // Cleanup: freeing memory: 
   for (i = 0; i < text_len + 1; i++)
      free (matrix[i]);
   free (matrix);
   return true;
}

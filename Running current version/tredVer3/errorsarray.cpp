//   |=================================================================|
//   |   TRED: a tool for detecting Tandem Repeats within sequences,   |
//   |               using the Edit Distance metric.                   |
//   |=================================================================|

// Copyright © 2007-2009 Dina Sokol, Justin Tojeira
// Distributed under the Aladdin Free Public License

// THIS SOFTWARE SHOULD BE ACCOMPANIED BY readme.txt AND license.html
// WE STRONGLY ENCOURAGE YOU TO READ BOTH BEFORE PROCEEDING
//------------------------------------------------------------------------------



/*/ Errors array
Returns an array whose values[i] are the maximum distance you can go in the strings with i errors
Parameters:
	int **matrix	-	the edit distance matrix calculated by build matrix
	int rows		-	the length of s1+1, (or number of rows in the matrix)
	int cols		-	the length of s2+2, (or number of cols in the matrix)
	int max_err		-	The maximum number of errors you will accept.  return array goes from 0 to max_err (inclusive)
	int direction	-	Indicates the direction that we are processing
							RIGHT (or any positive number not 0) to decide using the rightmost lowest value
							DOWN  (or any negative number) to decide using the lowest, rightmost value
							(Difference is precedence,  RIGHT cares more about lateral position in the matrix while
							DOWN cares more about the vertical direction in the matrix)
(Entry is a row,col entry in the edit distance matrix.)
Return Value:
	Entry *			-	Array comtaining edit distance values.  If the total edit distance of the array is less than
						max_errs, then all values from k(actual edit distance) to max_val will equal k.
/*/

#include "errorsarray.h"

#define COLUMN (i-(((j>max_err)?j-max_err:max_err-j)))

void errors_array(Entry* errors, int **m, int max_err, int direction)
{
	int i,j;
	//Entry *errors;

/*	errors= (Entry *)calloc(max_err+1,sizeof(Entry));
	if (errors == NULL)
	  {printf("error allocating errors array\n");  exit(1);  }    */

   // Zero out errors array
   for (i=0;i<=max_err;i++) {
		errors[i].row=0;
		errors[i].col=0;
	}

	if (direction>0) {
   	for (i=0; i<=max_err; i++)
      	for (j=max_err-i; j<=max_err+i; j++)
         	if (m[j][COLUMN]+j-max_err > errors[i].col){
               errors[i].col=m[j][COLUMN]+j-max_err;
               errors[i].row=m[j][COLUMN];
            }
	}
   else {
   	for (i=0; i<=max_err; i++)
      	for (j=max_err+i; j>=max_err-i; j--)
         	if (m[j][COLUMN] > errors[i].row){
            	errors[i].row=m[j][COLUMN];
               errors[i].col=m[j][COLUMN]+j-max_err;
            }
	}

	//return errors;
}


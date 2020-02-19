//   |=================================================================|
//   |   TRED: a tool for detecting Tandem Repeats within sequences,   |
//   |               using the Edit Distance metric.                   |
//   |=================================================================|

// Copyright © 2007-2009 Dina Sokol, Justin Tojeira
// Distributed under the Aladdin Free Public License

// THIS SOFTWARE SHOULD BE ACCOMPANIED BY readme.txt AND license.html
// WE STRONGLY ENCOURAGE YOU TO READ BOTH BEFORE PROCEEDING
//------------------------------------------------------------------------------



// This function builds a KxK edit distance matrix, which is a compression of a
// full (NxN) edit distance matrix.

//						~~~~~ EXPLAINING THE NxN MATRIX  ~~~~~

// This program is based on a matrix created from two strings, S1 and S2, of
// approximate length N, henceforth to be called the NxN matrix. The NxN matrix
// is created by putting S1 down the left side, S2 along the top, and filling
// in each element of the matrix to the number of mismatches, insertions or
// deletions that separate the two strings up to the point corresponding to
// the element in question.

// For example, if S1 were "ABCDEF" and S2 were "ACDBEH", the NxN matrix would
// be:

//         A  C  D  B  E  H
//      0  1  2  3  4  5  6
//   A  1  0  1  2  3  4  5
//   B  2  1  1  2  2  3  4
//   C  3  2 *1* 2  3  3  4
//   D  4  3  2  1  2  3  4
//   E  5  4  3  2  2  2  3
//   F  6  5  4  3  3  3  3

// In this example, the 1 surrounded by *s represents the number of mismatches
// between the substrings "ABC" and "AC". In this case, that mismatch is the
// omission of the letter "B" from the string "AC".

// This produces two certain results regarding the diagonals of the matrix:
// 1) Each diagonal has a minimum number of errors equal to its distance from
//		the main diagonal.
// 2) As each diagonal progresses the number may only stay the same, or increase
//		by 1 over the previous element on the same diagonal.

// Because of these properties, given a maximum number of errors allowable, it
// is possible to disregard a large portion of the NxN matrix. Any diagonal
// farther from the main than the number of allowable errors may be immediately
// disregarded, as may any diagonal which goes over the number of allowable
// errors. Then, we may store this matrix in a smaller matrix, to be called the
// KxK matrix, deriving its name from the amount of allowable errors, K.


//						~~~~~ EXPLAINING THE KxK MATRIX  ~~~~~

// The KxK matrix does not record the number of mismatches in each element, but
// rather records when that number increases on each diagonal in the NxN matrix.
// The diagonals which interest us are the main diagonal, and K diagonals in
// either direction. Each gets a row in the KxK matrix, giving the matrix 2K+1
// rows. In each row, it is recorded how far down, in the NxN matrix, the
// corresponding diagonal reaches for each possible number of errors, up to and
// including the maximum allowable errors (K). However, because each diagonal
// has a different number of minimum errors possible, each row in the KxK matrix
// has a different number of elements. Thus, the KxK matrix is formatted as
// follows:

//          NxN matrix           Corresponding KxK matrix allowing 3 errors
//
//         A  C  D  B  E  H                  0  1  2  3
//      0  1  2  3  4  5  6              -3           6
//   A  1  0  1  2  3  4  5              -2       *5* 6
//   B  2  1  1  2  2  3  4              -1     4  5  6
//   C  3  2  1  2  3  3  4               0  1  2  5  6
//   D  4  3  2  1  2  3  4               1     1  2  5
//   E  5  4  3  2  2  2  3               2        2  3
//   F  6  5  4  3  3  3  3               3           2

// The diagonals are numbered along the left, with 0 being the main diagonal,
// negative numbers being diagonals left of the main, and positive to the right.
// The number of errors are listed along the top, and as noted, diagonals can't
// have fewer errors than their distance from the main diagonal. The elements of
// the matrix represent the lowest row in the NxN matrix along the corresponding
// diagonal that has the corresponding number of errors. To calculate the actual
// position in the NxN matrix from this information, simply use the element as
// the row, and add the value of the diagonal to that to find the column.
// Thus, in this example, the 5 surrounded by *s represents the value in the
// NxN matrix at row 5 and column 3, which as we see is the number 2, and is
// in fact the last 2 on the -1 diagonal, which is what the KxK matrix states.


//				  	  	 ~~~~~ NOTES ABOUT THE PROGRAM ~~~~~

// Buildmatrix does not allocate space for the entire NxN matrix, but rather
// calculates it a row at a time, needing to know only the row it's working on
// and the row above it. Thus, it stores the row it's working on in the array
// "lower" and the row above in the array "upper". After it finishes a row, it
// switches the pointers to upper and lower so that the row it just finished is
// in "upper", and it can calculate the next row in "lower". Buildmatrix
// attempts to calculate as little of the NxN matrix as it needs to. It begins
// tracking only the diagonals which start with the allowed number of errors
// or less, and stops tracking diagonals when they exceed the allowed errors.
// To do this it uses the variables "low", "high", and "indent", processing
// each row from (indent + low) to (indent + high). Indent increases by 1 for
// each subsequent row because the area of interest to us moves diagonally.
// Meanwhile, low and high are used as boundaries for the area we need to
// process.

// While it processes the NxN matrix in this fashion, buildmatrix checks each
// result against the previous result on the same diagonal, in order to note
// when the number of errors on a given diagonal increases. When it notes an
// increase, it invokes the macro FILL_KxK_MATRIX, which puts the value of the
// previous row (which is the last row to have a certain number of errors, since
// the current row has one more) into the appropriate cell of the KxK matrix.

//				  	  	   ~~~~~ INPUT AND OUTPUT ~~~~~
// THIS VERSION HAS BEEN MODIFIED TO FILL UP AN ALREADY ALLOCATED KxK MATRIX

// Input are the character strings S1 and S2, and "errors", which is the number
// of allowable errors. In creating the NxN matrix, S1 goes on the left and S2
// on top.

//				  	  	  ~~~~~ POINTS OF INTEREST ~~~~~

// In the KxK matrix (see diagram in "explaining the KxK matrix", above), note
// that not all rows possess the same number of elements. Note also that the
// first element in the top row contains the maximum number of errors, the
// second row starts at one less error, and so on. Because of this, the matrix
// cannot simply be indexed by (diagonal,errors). Rather, calculations similar
// to that used in the macro FILL_KxK_MATRIX must be used.

// S1 goes down the left side of the NxN matrix, S2 along the top.

// Because of the diagonal qualities of the matrix, the area of interest shifts
// one to the right each time we go down a row. Thus the column in which "upper"
// and "lower" start increases by one each time we go down a row. Thus,
// lower[i] is not directly under upper[i] in the NxN matrix, but rather they
// are on the same diagonal, with lower[i] being the element right after
// upper[i] on the same diagonal.

//				  	  	 ~~~~~ VERSION 2 MODIFICATIONS ~~~~~

// The parameter "limit" has been added to suit the case of one string being
// a suffix of the other string. In such a case, we will not want the string to
// "catch up" with its suffix. To prevent this, limit, when positive, will
// denote that S1 is a suffix of S2. When negative, S2 is a suffix of S1. In
// either case, the difference between the strings will be equal to the absolute
// value of limit.

// A macro "UNKNOWN" has been added for usage with biological sequences.
// UNKNOWN should be defined as the character that represents an unknown element
// of the sequence (For example, the default is 'N' which represents an unknown
// base pair in a DNA sequence). This character will always be considered a
// mismatch, even against itself. UNKNOWN is defined in buildk.h


#include "buildk.h"                  
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define FILL_KxK_MATRIX m[i][upper[i]-((errors>i)?errors-i:i-errors)]=j-1
	// fills in one cell of the KxK matrix, generally whenever a cell in the NxN
	// matrix is calculated which is greater than the previous cell on the same
	// diagonal

void buildmatrix (int **m, char * s1, char * s2, int errors, int limit)
{
//	int   ** m;                   // Pointer to the KxK matrix

   int   rows = 2 * errors + 1;	// Number of diagonals of interest in the NxN
											// matrix. Each diagonal of interest in the NxN
									      // matrix is represented by 1 row in the KxK

   int   cols = errors + 1;      // Number of columns in the KxK matrix
   int   reps = strlen(s1) + 1;  // Number of rows in the NxN matrix
   int	length = strlen(s2);    // Number of columns in the NxN matrix

   int   indent = 0-errors,   low = errors,   high = (length<errors)?errors+length:2*errors;
   // The diagonals of interest in the NxN matrix are calculated by calculating
   // the elements of each row that fall between columns indent + low and
   // indent + high. indent increases by 1 for each subsequent row of the
   // matrix. low and high are adjusted based on any outermost diagonals that
   // are determined to no longer be of interest.

   if (limit>0)
   	high = (high<errors+limit-1)? high:errors+limit-1;

   int   *upper, *lower, *temp;
   // lower holds the information of the row currently being calculated, and
   // upper holds the information of the row above that, on which the
   // calculations are based. temp is used to swap pointers so that after a row
   // is finished upper points to it and lower is free for the next calculations

   int   i,j,k;		// loop control variables

/******************************MEMORY ALLOCATION*******************************/
   upper = (int *) malloc (rows * sizeof(int));
   if (upper == NULL){
   	printf("Memory allocation failed");
      exit(1);
   }
   lower = (int *) malloc (rows * sizeof(int));
   if (lower == NULL){
   	printf("Memory allocation failed");
      exit(1);
   }
/*   m = (int **) malloc (rows * sizeof(int *));
   if (m == NULL){
   	printf("Memory allocation failed");
      exit(1);
   }
	for (i=0;i<rows;i++){
      m[i] = (int *) malloc (((i<cols)?i+1:rows-i+1) * sizeof(int));
      if (m[i] == NULL){
   		printf("Memory allocation failed");
      	exit(1);
      }
   }

/**********************BEGINS CALCULATIONS OF NxN MATRIX***********************/

						// PRINTS TITLE, "THE NxN MATRIX"
//	printf("THE NxN MATRIX:\n");

						// PRINTS STRING S2 ALONG TOP OF NxN MATRIX
/*	printf("      ");
   for (i=0;i<length;i++)
   	printf("%c  ",s2[i]);
   printf("\n ");                     */

					// PUTS THE TOP ROW OF THE NxN MATRIX INTO UPPER
   for (i=low;i<=high;i++){
   	upper[i]=i-low;
//      printf("%3d",upper[i]);          // printing
   }
//   for (k=indent+i-1;k<length;k++)     // printing
//   	printf("  +");

   if (errors<0){                    // errors must be a non-negative integer
   	printf("\n\nNegative error limit is invalid");
   	exit(1);
   }
   if (errors==0){           		      // when errors == 0 only the main
      j=1;                             // diagonal is of interest
      i=0;
      goto onediagonal;
   }

/******************CALCULATES THE BULK OF THE NxN MATRIX***********************/
/*******************FILLS IN THE KxK MATRIX AS IT GOES*************************/

   for (j=1;j<reps;j++){	// J IS THE ROW BEING WORKED ON

      					// FIRST ELEMENT OF THE ROW
      indent++;

      if (errors>=j && !(limit<0 && low-errors-1<=limit)) {
                                       // for the first few rows where the area
      	i=--low;                      // of interest begins along the leftmost
         lower[i]=j;                   // side of the matrix

//         printf("\n%c",s1[j-1]);       // printing
//         printf("%3d",lower[i]);
      }

      else {                             // processes the first element of the
      	i=low;                          // row when there are no special cases
         if (s1[j-1]==s2[indent+i-1]&&s1[j-1]!=UNKNOWN)
         	lower[i]=upper[i];
         else lower[i]=upper[i]+1;
         if (upper[i+1]+1<lower[i])
         	lower[i]=upper[i+1]+1;
         if (lower[i]>upper[i]) {
         	FILL_KxK_MATRIX;
            
         	while (lower[i]>errors) {      // increases the lower boundary (low)
            	i=++low;                    // if the leftmost diagonal has gone
               if (low==high) {            // over the error limit. Continues
               	indent--;                // doing this as long as the leftmost
                  goto onediagonal;         // diagonal is over the error limit
               }                             // and does not meet the rightmost
         		if (s1[j-1]==s2[indent+i-1]&&s1[j-1]!=UNKNOWN)   // diagonal.
         			lower[i]=upper[i];
         		else lower[i]=upper[i]+1;
         		if (upper[i+1]+1<lower[i])
         			lower[i]=upper[i+1]+1;
               if (lower[i]>upper[i]&&upper[i]<=errors)
         			FILL_KxK_MATRIX;
            }
         }

/*      	printf("\n%c",s1[j-1]);         // printing
      	for (k=0;k<indent+low;k++)
      		printf("  +");
   		printf("%3d",lower[i]);       */

      }
      i++;
      
      					// MIDDLE ELEMENTS OF THE ROW
   	// processes the middle elements of each row. Because the boundaries of
      // the part of the matrix we're dealing with is adjusted when we process
      // the first and last elements of each row, the middle elements have no
      // special cases.
      for (;i<high;i++){
      	if (s1[j-1]==s2[indent+i-1]&&s1[j-1]!=UNKNOWN)
         	lower[i]=upper[i];
         else lower[i]=upper[i]+1;
         if (upper[i+1]+1<lower[i])
         	lower[i]=upper[i+1]+1;
         if (lower[i-1]+1<lower[i])
         	lower[i]=lower[i-1]+1;
         if (lower[i]>upper[i]&&upper[i]<=errors)
         	FILL_KxK_MATRIX;

/*         if (lower[i]<=errors)         // printing
         	printf("%3d",lower[i]);
         else
         	printf("  x");                    */
      }

      					// LAST ELEMENT OF THE ROW

      // when the rightmost diagonal goes past the right edge of the NxN matrix,
      // that diagonal is finished and the corresponding row on the KxK matrix
      // must be finished off as well
      if (indent+high>length){
      	for (k=upper[i];k<=errors;k++)
         	m[i][k-((errors>i)?errors-i:i-errors)]=j-1;
         high--;
      }

      // processes the last element of the row when there are no special cases
      else {
      	if (s1[j-1]==s2[indent+i-1]&&s1[j-1]!=UNKNOWN)
      		lower[i]=upper[i];
      	else lower[i]=upper[i]+1;
      	if (lower[i-1]+1<lower[i])
      		lower[i]=lower[i-1]+1;
         if (lower[i]>upper[i])
         	FILL_KxK_MATRIX;

/*         if (lower[i]<=errors)           // printing
      		printf("%3d",lower[i]);
         else
         	printf("  y");               */
      }

      while (lower[high]>errors&&high>low)    // decreases the upper limit if
      	high--;                              // the rightmost diagonal has gone
                                              // over the error limit

//      for (k=indent+i;k<length;k++)      // printing
//      	printf("  +");

      temp = upper;         // switch upper and lower so that the next row may
      upper = lower;        // be calculated off this one, which is being put
      lower = temp;         // into upper

      if (low==high){
      	i=high;
         j++;
         break;
      }
	}

/*************WHEN THE AREA OF INTEREST IS DOWN TO ONE DIAGONAL ***************/
onediagonal:

   for (;j<reps;j++){
   	indent++;                          // processes the one element in the row
      if (s1[j-1]==s2[indent+i-1]&&s1[j-1]!=UNKNOWN)
      	lower[i]=upper[i];
      else lower[i]=upper[i]+1;
      if (lower[i]>upper[i]) {
      	FILL_KxK_MATRIX;
      	if (lower[i]>errors) break;
      }
      
/*      printf("\n%c",s1[j-1]);            // printing
      for (k=0;k<indent+i;k++)
      	printf("  +");
      printf("%3d",lower[i]);
      for (;k<length;k++)
      	printf("  +");                */

      temp = upper;
      upper = lower;
      lower = temp;
      }

/********* PUTS THE BOTTOM ROW OF THE NxN MATRIX IN THE KxK MATRIX ************/

	// If a diagonal reaches the bottom row of the NxN matrix with the allowed
   // number of errors or fewer, then the bottom row is counted as the last row
   // with that number of errors. Additionally, if the number of errors is fewer
   // than the allowed number, then the bottom row of the NxN matrix will be
   // counted as the last row with the number of errors from the current up to
   // and including the maximum allowed. For instance, if we allow 3 errors
   // and reach the bottom row with 2 in diagonal -1, then the bottom row will
   // be considered the last row that has both 2 errors and 3 errors for
   // diagonal -1, and be entered into the KxK matrix accordingly.

	for (i=low;i<=high;i++){
   	for (k=upper[i];k<=errors;k++)
      	m[i][k-((errors>i)?errors-i:i-errors)]=j-1;
   }

/********* CASE IF ALLOWED ERRORS > STRING LENGTH OR LIMIT IS REACHED**********/
  	for (i=reps-1; i<errors; i++)
   	for (j=0; j<errors-i; j++)
      	m[errors-i-1][j]=-1;

   for (i=length; i<errors; i++)
   	for (j=0; j<errors-i; j++)
      	m[errors+i+1][j]=-1;

   if (limit>0){
   	for (i=limit-1; i<errors; i++)
   		for (j=0; j<errors-i; j++)
      		m[errors+i+1][j]=-1;
   }
   else{
   	for (i=-limit-1; i<errors; i++)
   		for (j=0; j<errors-i; j++)
      		m[errors-i-1][j]=-1;
   }

/*************************** PRINTS THE KxK MATRIX ****************************/
/*	printf("\n\nTHE KxK MATRIX:\n   ");
   for (i=0;i<cols;i++)
   	printf("%3d",i);

 	for (j=0;j<rows;j++){
   	printf("\n%3d",j-errors);
      for (i=0;i<((j>errors)?j-errors:errors-j);i++)
      	printf("  +");
   	for (i=0;i<((j<errors)?j+1:rows-j);i++)
      	printf("%3d",m[j][i]);
   }

   printf("\n\n");

/************************** DEALLOCATES AND RETURNS ***************************/

   free (upper);
   free (lower);
//   return m;
}


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

//	   		 ~~~~~ COMPUTER REPRESENTATION OF THE KxK MATRIX ~~~~~

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

//					   	~~~~~ TRACEBACK ALGORITHM ~~~~~

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

//					   	~~~~~ FUNCTION SPECIFICATIONS ~~~~~
//		THE PARAMETERS HAVE BEEN MODIFIED IN THIS VERSION TO FIT WITH MAIN

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

//				  	  	  ~~~~~ POINTS OF INTEREST ~~~~~

// Traceback does not actually recreate any part of the NxN matrix.
// Traceback outputs the aligned strings, output1 and output2, backwards (I'm
// pretty sure there was some reason to do it this way but I can't remember
// what it was).

//                ~~~~~ VERSION 2 MODIFICATIONS ~~~~~

// A macro "UNKNOWN" has been added for usage with biological sequences.
// UNKNOWN should be defined as the character that represents an unknown element
// of the sequence (For example, the default is 'N' which represents an unknown
// base pair in a DNA sequence). This character will always be considered a
// mismatch, even against itself. UNKNOWN is defined in buildk.h

#include "tracek.h"
#include "buildk.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define PRINT_INFO printf("row=%d (%d) col=%d (%d) s1pos=%d s2pos=%d\n",row-max_errors,row,col+(row>max_errors?row-max_errors:max_errors-row),col,s1pos,s2pos)
	// This macro prints the row, column, s1pos and s2pos we're up to during
   // the course of the function. For row and column it prints the values of
   // our representation of the KxK matrix, with the computer's internal values
   // in parentheses.

void traceback (int** m, int max_errors, char* s1, char* s2, char* o1, char* o2) //, int* matches)
{
	traceback (m,max_errors,s1,s2,o1,o2,strlen(s1)-1,strlen(s2)-1);
}

void traceback (int** m, int max_errors, char* s1, char* s2, char* o1, char* o2, int s1pos, int s2pos) //, int* matches, int s1pos, int s2pos)
{

/******************************MEMORY ALLOCATION*******************************/
/*	int max_string_length = ((s1pos<s2pos)?s1pos:s2pos) + max_errors + 2;

	char * o1 = (char *) malloc (max_string_length * sizeof(char));
   if (o1 == NULL){
   	printf("Memory allocation failed");
      exit(1);
   }
	char * o2 = (char *) malloc (max_string_length * sizeof(char));
   if (o2 == NULL){
   	printf("Memory allocation failed");
      exit(1);
   }                                                                     */

/************************* VARIABLE INITIALIZATION ****************************/
	// row and col refer to the computer representation of the KxK matrix.
   // o1pos and o2pos are used to keep position on output1 and output2.
	int row = s2pos - s1pos + max_errors;
   int col = max_errors-((row>max_errors)?row-max_errors:max_errors-row);
   int o1pos = 0, o2pos = 0;
//   *matches=0;

/******************** CHECK TO MAKE SURE ROW IS IN BOUNDS *********************/
   if (row<0||row>2*max_errors){
   	printf("\n\n          *!*!*!*!*!*!*!*!* ERROR *!*!*!*!*!*!*!*!*!*\n\n");
   	printf("          Traceback cannot be started from that point\n");
      //PRINT_INFO;
      printf("          *!*!*!*!*!*!*!*!* ERROR *!*!*!*!*!*!*!*!*!*\n\n");
      //printf("Press <enter> twice\n");
		//getchar();
      //getchar();
      o1[0]='\0';
      o2[0]='\0';
      //*output1=o1;
      //*output2=o2;
      return;
   }

/************** DECREMENTS COL WHEN THE ERRORS WE'RE ALLOWING ARE *************/
/************** LESS THAN THE ACTUAL ERRORS AT OUR STARTING POINT *************/
//   printf("Traceback starting\n");
//   PRINT_INFO;

   while (col>0&&m[row][col-1]>s1pos){
   	col--;
      //PRINT_INFO;
   }

/*************************** MAIN TRACEBACK LOOP ******************************/

                // LOOPS UNTIL WE ARE ON THE LEFTMOST COLUMN
                // OR TOP ROW OF THE VIRTUAL NxN MATRIX
   while (s1pos>=0&&s2pos>=0){

      if (row>max_errors&&m[row-1][col]>s1pos){
      	o1[o1pos]='-';
         o2[o2pos]=s2[s2pos];      // left move in the virtual NxN
         row--;                    // up and left move in our KxK
         s2pos--;                  // up move in computer's KxK
         o1pos++;
         o2pos++;
      }

      else if (row<max_errors&&m[row+1][col]>=s1pos){
      	o1[o1pos]=s1[s1pos];
         o2[o2pos]='-';             // up move in the virtual NxN
         row++;                     // down and left move in our KxK
         s1pos--;                   // down move in computer's KxK
         o1pos++;
         o2pos++;
      }

      else if (row<=max_errors&&col>1&&m[row-1][col-2]>s1pos){
      	o1[o1pos]='-';
         o2[o2pos]=s2[s2pos];         // left move in the virtual NxN
         row--;                       // up and left move in our KxK
         col = col - 2;               // up and 2 left move in computer's KxK
         s2pos--;
         o1pos++;
         o2pos++;
      }

      else if (row>=max_errors&&col>1&&m[row+1][col-2]>=s1pos){
      	o1[o1pos]=s1[s1pos];
         o2[o2pos]='-';             // up move in the virtual NxN
         row++;                     // down and left move in our KxK
         col = col - 2;             // down and 2 left move in computer's KxK
         s1pos--;
         o1pos++;
         o2pos++;
      }

      else if (col>0&&m[row][col-1]>=s1pos){
      	o1[o1pos]=s1[s1pos];
         o2[o2pos]=s2[s2pos];         // diagonal move in the virtual NxN
         col--;                       // left move in our KxK
         s1pos--;                  	  // left move in the computer's KxK
         s2pos--;
         o1pos++;
         o2pos++;
      }

   	else if (s1[s1pos]==s2[s2pos]&&s1[s1pos]!=UNKNOWN){
      	o1[o1pos]=s1[s1pos];
         o2[o2pos]=s2[s2pos];         // in the case of the characters we're up
         s1pos--;                     // to in S1 and S2 matching, we make a
         s2pos--;                     // diagonal move in the virtual NxN and
         o1pos++;                     // no move in either representation of the
         o2pos++;                     // KxK
//         *matches++;
      }

      else{
      	printf("\n\n          *!*!*!*!*!*!*!*!* ERROR *!*!*!*!*!*!*!*!*!*\n\n");
         printf("                     Error in traceback\n\n");
         printf("          *!*!*!*!*!*!*!*!* ERROR *!*!*!*!*!*!*!*!*!*\n\n");
         //printf("Press <enter> twice\n");
         //getchar();
         //getchar();
         o1[o1pos]='\0';
         o2[o2pos]='\0';            // if none of the above cases are true,
         //*output1=o1;               // there is an error in the traceback,
         //*output2=o2;               // probably due to inconsistent
         return;                    // passed parameters
      }
      //PRINT_INFO;
   }                

/************************ FINISHES UP THE TRACEBACK ***************************/

           // TRACES UP ALONG THE LEFTMOST COLUMN OF THE NxN MATRIX IF NEEDED
   while (s1pos>=0){
   	o1[o1pos]=s1[s1pos];
      o2[o2pos]='-';
      s1pos--;
      o1pos++;
      o2pos++;
      row++;
      //PRINT_INFO;
   }

           // TRACES LEFT ALONG THE TOP ROW OF THE NxN MATRIX IF NEEDED
   while (s2pos>=0){
   	o1[o1pos]='-';
      o2[o2pos]=s2[s2pos];
      s2pos--;
      o1pos++;
      o2pos++;
      row--;
      //PRINT_INFO;
   }



           // NULL-TERMINATES O1 AND O2
   o1[o1pos]='\0';
   o2[o2pos]='\0';

//   strrev(o1);
//   strrev(o2);

   //printf("Traceback finished\n\n");

   //*output1=o1;
   //*output2=o2;

}

/*char *strrev(char *string) {

	char *start = string;
	char *left = string;
	char ch;

	while(*string++)
		;
	string-=2;
	while(left<string) {
		ch = *left;
		*left++ = *string;
		*string-- = ch;
	}
	return start;
}*/


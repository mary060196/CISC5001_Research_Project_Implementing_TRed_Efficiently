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


#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "errorsarray.h"
#include "createarray.h"
#include "oneiteration.h"
#include "parameters.h"

char* myfgets(char* s, int n, FILE* stream);

int main(int argc, char *argv[])
{
   char filename[40];
   FILE *infile, *outfile;
   int i,j,k;
   struct LastReportedRepeat left, right;

   /*********** FOR INPUT/OUTPUT FILES AS COMMAND LINE ARGUMENTS ************/

   if (argc!=3){
   	printf("File parameters required. See readme for details\n");    
      exit(1);
   }
	if ((infile=fopen(argv[1],"r"))==NULL)
   		{ printf("Error opening sequence file.\n"); exit(1); }
   if ((outfile=fopen(argv[2],"w"))==NULL)
   	{ printf("Error opening intermediary file.\n"); exit(1); }

   /**************** TO MANUALLY ENTER INPUT/OUTPUT FILES *******************/

/*   printf("Enter sequence file: ");
   cin >> filename;
	if ((infile=fopen(filename,"r"))==NULL)
   		{ printf("Error opening sequence file.\n"); exit(1); }

	printf("Enter intermediary file: ");
   cin >> filename;
   if ((outfile=fopen(filename,"w"))==NULL)
   	{ printf("Error opening intermediary file.\n"); exit(1); }

  /************************* MEMORY ALLOCATION ******************************/

	char *str=(char*) malloc((2*MAX_PERIOD+1)*sizeof(char));
	if (str==NULL)
	  { printf("not enough memory in main.\n"); exit(1);  }
	char *doublestr=(char*) malloc((4*MAX_PERIOD+1)*sizeof(char));
	if (doublestr==NULL)
	  { printf("not enough memory in main.\n"); exit(1);  }
	char *aheadstr=(char*) malloc((2*MAX_PERIOD+1)*sizeof(char));
	if (aheadstr==NULL)
	  { printf("not enough memory in main.\n"); exit(1);  }
   char *substring = (char*) malloc((2*MAX_PERIOD+1)*sizeof(char));
   if (substring==NULL)
     { printf("not enough memory in main.\n"); exit(1);  }

   int ** matrix_forward = (int **) malloc ((2*MAX_ERRORS+1) * sizeof(int *));
   if (matrix_forward == NULL) { printf("Memory allocation failed"); exit(1); }
	for (i=0;i<2*MAX_ERRORS+1;i++){
      matrix_forward[i] = (int *) malloc (((i<=MAX_ERRORS)?i+1:2*MAX_ERRORS-i+2) * sizeof(int));
      if (matrix_forward[i] == NULL){ printf("Memory allocation failed"); exit(1); }
   }
   int ** matrix_reverse = (int **) malloc ((2*MAX_ERRORS+1) * sizeof(int *));
   if (matrix_reverse == NULL){ printf("Memory allocation failed"); exit(1); }
	for (i=0;i<2*MAX_ERRORS+1;i++){
      matrix_reverse[i] = (int *) malloc (((i<=MAX_ERRORS)?i+1:2*MAX_ERRORS-i+2) * sizeof(int));
      if (matrix_reverse[i] == NULL){ printf("Memory allocation failed"); exit(1); }
   }
   int ** left_matrix_forward = (int **) malloc ((2*MAX_ERRORS+1) * sizeof(int *));
   if (left_matrix_forward == NULL){ printf("Memory allocation failed"); exit(1); }
	for (i=0;i<2*MAX_ERRORS+1;i++){
      left_matrix_forward[i] = (int *) malloc (((i<=MAX_ERRORS)?i+1:2*MAX_ERRORS-i+2) * sizeof(int));
      if (left_matrix_forward[i] == NULL){ printf("Memory allocation failed"); exit(1); }
   }
   int ** left_matrix_reverse = (int **) malloc ((2*MAX_ERRORS+1) * sizeof(int *));
   if (left_matrix_reverse == NULL) { printf("Memory allocation failed"); exit(1); }
	for (i=0;i<2*MAX_ERRORS+1;i++){
      left_matrix_reverse[i] = (int *) malloc (((i<=MAX_ERRORS)?i+1:2*MAX_ERRORS-i+2) * sizeof(int));
      if (left_matrix_reverse[i] == NULL){ printf("Memory allocation failed"); exit(1); }
   }
   int ** right_matrix_forward = (int **) malloc ((2*MAX_ERRORS+1) * sizeof(int *));
   if (right_matrix_forward == NULL) { printf("Memory allocation failed"); exit(1); }
	for (i=0;i<2*MAX_ERRORS+1;i++){
      right_matrix_forward[i] = (int *) malloc (((i<=MAX_ERRORS)?i+1:2*MAX_ERRORS-i+2) * sizeof(int));
      if (right_matrix_forward[i] == NULL){ printf("Memory allocation failed"); exit(1); }
   }
   int ** right_matrix_reverse = (int **) malloc ((2*MAX_ERRORS+1) * sizeof(int *));
   if (right_matrix_reverse == NULL) { printf("Memory allocation failed"); exit(1); }
	for (i=0;i<2*MAX_ERRORS+1;i++){
      right_matrix_reverse[i] = (int *) malloc (((i<=MAX_ERRORS)?i+1:2*MAX_ERRORS-i+2) * sizeof(int));
      if (right_matrix_reverse[i] == NULL){ printf("Memory allocation failed"); exit(1); }
   }

 	Entry *down_errors= (Entry *)malloc((MAX_ERRORS+1)*sizeof(Entry));
 	Entry *right_errors= (Entry *)malloc((MAX_ERRORS+1)*sizeof(Entry));
 	Entry *left_down_errors= (Entry *)malloc((MAX_ERRORS+1)*sizeof(Entry));
 	Entry *left_right_errors= (Entry *)malloc((MAX_ERRORS+1)*sizeof(Entry));
 	Entry *right_down_errors= (Entry *)malloc((MAX_ERRORS+1)*sizeof(Entry));
 	Entry *right_right_errors= (Entry *)malloc((MAX_ERRORS+1)*sizeof(Entry));
   if (down_errors==NULL || right_errors==NULL || left_down_errors==NULL || left_right_errors==NULL || right_down_errors==NULL || right_right_errors==NULL)
	   {  printf("error allocating errors array\n");  exit(1);  }

  /********************** VARIABLE INITIALIZATION ***************************/

   left.j = 0; left.length = 0; left.rating = 0;
   left.left_errs = 0; left.right_errs = 0; left.total_errs = 0;
   left.matrix_forward = left_matrix_forward;
   left.matrix_reverse = left_matrix_reverse;
   left.down_errors = left_down_errors;
   left.right_errors = left_right_errors;
   left.newhigh=0;
   right.j = 0; right.length = 0; right.rating = 0;
   right.left_errs = 0; right.right_errs = 0; right.total_errs = 0;
   right.matrix_forward = right_matrix_forward;
   right.matrix_reverse = right_matrix_reverse;
   right.down_errors = right_down_errors;
   right.right_errors = right_right_errors;
   right.newhigh=0;
   doublestr[0] = '\0';

  /****************** ALGORITHMIC PORTION OF MAIN FUNCTION ********************/

   // Performs the (usually) recursive part of the Main-Lorentz algorithm
   // iteratively. The algorithm is modified so as to work with overlapping
   // parts, allowing it to process longer input strings.

   long offset = START_POS;
   int num, levels;

   myfgets(str,2*MAX_PERIOD+1,infile);
   while (true) {
      if (PROCESSING)
   		printf("\nProcessing");
   	if (offset!=START_POS) {
      	strcat(doublestr,str);
         OneIteration(doublestr, matrix_forward, matrix_reverse, down_errors, right_errors, left, right, offset-2*MAX_PERIOD, 0, outfile, 0);
         fprintf(outfile,"-1 %ld\n\n",offset);
      }

      if (!myfgets(aheadstr,2*MAX_PERIOD+1,infile)){
      	aheadstr[0]='\0';
         break;
      }
      if (strlen(aheadstr)!=2*MAX_PERIOD) break;

   	char *strpnt = str;
      levels = 0;
      int *numarray = createarray(str, MIN_LENGTH, &levels);
      for (i=0,j=1;i<levels;i++){
         if (PROCESSING)
      		printf(".");
      	for (k=0;k<j;k++){
         	num = numarray[j+k-1];
            strncpy(substring,strpnt,num);
            substring[num]='\0';
            OneIteration(substring, matrix_forward, matrix_reverse, down_errors, right_errors, left, right, offset+(long int)(strpnt-str), 0, outfile, i+1);
            strpnt+=num;
         }
         j=j*2;
         strpnt=str;
      }
      offset += 2*MAX_PERIOD;
      free (numarray);
      strcpy(doublestr,str);
      strcpy(str,aheadstr);
   }

   strcpy(doublestr,str);
   strcat(doublestr,aheadstr);
   OneIteration(doublestr, matrix_forward, matrix_reverse, down_errors, right_errors, left, right, offset, 1, outfile, 0);
   char *strpnt = doublestr;
   num = 0; levels = 0;
   int *numarray = createarray(doublestr, MIN_LENGTH, &levels);
   for (i=1,j=2;i<levels;i++){
   	if (PROCESSING)
      	printf(".");
   	for (k=0;k<j;k++){
      	num = numarray[j+k-1];
         strncpy(substring,strpnt,num);
         substring[num]='\0';
         OneIteration(substring, matrix_forward, matrix_reverse, down_errors, right_errors, left, right, offset+(long int)(strpnt-doublestr), k==j-1?1:0, outfile, i+1);
         strpnt+=num;
      }
   	j=j*2;
   	strpnt=doublestr;
   }
   fprintf(outfile,"-1 -1\n\n");
   free (numarray);

  /************************ MEMORY DEALLOCATION *****************************/

   free (str); free (doublestr); free (aheadstr); free (substring);
   for (i=0;i<2*MAX_ERRORS+1;i++){
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

   fflush(NULL);
   fclose (infile);
   fclose (outfile);

	return 0; // END MAIN
}


  /************************** MYFGETS FUNCTION ******************************/

// Version of fgets that removes all white space.
// Reads n-1 characters or until end of file. Does not stop at a new line.
// A null byte is appended to s to mark the end of the string.
char* myfgets(char* s, int n, FILE* stream)
{
	char* pos = s;
   char next;

   do {
      next=getc(stream);
      if (next==EOF) return NULL;
   } while (next==' '||next=='\n'||next=='\t'||next=='\r');
   *pos++=next;
   for (int i=1; i<n-1; ){
   	next=getc(stream);
      if (next==EOF) break;
      if (next!=' '&&next!='\n'&&next!='\t'&&next!='\r'){
      	*pos++=next;
         i++;
      }
   }
   *pos='\0';
   return s;
}


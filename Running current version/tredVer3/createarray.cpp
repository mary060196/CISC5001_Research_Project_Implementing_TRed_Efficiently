//   |=================================================================|
//   |   TRED: a tool for detecting Tandem Repeats within sequences,   |
//   |               using the Edit Distance metric.                   |
//   |=================================================================|

// Copyright © 2007-2009 Dina Sokol, Justin Tojeira
// Distributed under the Aladdin Free Public License

// THIS SOFTWARE SHOULD BE ACCOMPANIED BY readme.txt AND license.html
// WE STRONGLY ENCOURAGE YOU TO READ BOTH BEFORE PROCEEDING
//------------------------------------------------------------------------------



#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>

int* createarray (char* str, int min, int *levels){
	int length = strlen(str);
   int i,j;

   j = (length-1)/min +1;
   for (i=1;i<j;i*=2)
   	(*levels)++;

   int *numarray = (int*) malloc ((i*2-1)*sizeof(int));
   if (numarray == NULL){
   	printf("Memory allocation error in createarray");
      exit(1);
   }

   numarray[0]=length;
   for (j=0;j<i-1;j++){
   	numarray[j*2+1]=numarray[j]/2;
      numarray[j*2+2]=numarray[j]-numarray[j*2+1];
   }

   return numarray;
}


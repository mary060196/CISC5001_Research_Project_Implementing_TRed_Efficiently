//   |=================================================================|
//   |   TRED: a tool for detecting Tandem Repeats within sequences,   |
//   |               using the Edit Distance metric.                   |
//   |=================================================================|

// Copyright © 2007-2009 Dina Sokol, Justin Tojeira
// Distributed under the Aladdin Free Public License

// THIS SOFTWARE SHOULD BE ACCOMPANIED BY readme.txt AND license.html
// WE STRONGLY ENCOURAGE YOU TO READ BOTH BEFORE PROCEEDING
//------------------------------------------------------------------------------



// Prints out a condensed form of a 2-sequence alignment. At the beginning, and
// after each insertion or deletion, the position of each sequence is printed.
// This is sufficient information to represent the layout of any optimal,
// maximal 2-sequence alignment (which are the only kind we are interested in).
// This condensed form is used in a data file  which is then read in by the
// post-processor.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int printcompact(char *seq1, char *seq2, long int pos1, long int pos2, FILE* outfile)
{
	bool state = 0; // state 1 means the previous position had an insertion or deletion

	fprintf(outfile,"%ld %ld ", pos1, pos2);
   while (*seq1!='\0' && *seq2!='\0'){
      if (*seq1=='-'){
      	pos2++;
         state=1;
      }
      else if (*seq2=='-'){
      	pos1++;
         state=1;
      }
   	else {
      	if (state==1){
         	fprintf(outfile,"%ld %ld ", pos1, pos2);
         	state=0;
         }
         pos1++;
         pos2++;
      }
      seq1++;
      seq2++;
   }
   fprintf(outfile,"%ld 0\n\n", pos1);
   return 0;
}


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
// ONEITERATION, which together with the function MAIN constitutes the full
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
#include "buildk.h"
#include "tracek.h"
#include "errorsarray.h"
#include "oneiteration.h"
#include "printcompact.h"
#include "parameters.h"

int OneIteration(char* str, int ** matrix_forward, int ** matrix_reverse, Entry* down_errors, Entry* right_errors, struct LastReportedRepeat left, struct LastReportedRepeat right, long int offset, bool isend, FILE* outfile, int level)

{
	char *mid_forward, *mid_backward, *right_forward, *right_backward, *left_forward, *left_backward;
   char *rev_pattern, *rev_text, *mid_text, *mid_pattern;
	int i,j,k;
	int begpos, endpos;
   int length = strlen(str), n = length/2;
   int curr_left, curr_right, left_errs, right_errs, curr_errs, curr_length, curr_rating, curr_off;
   int comparison, diff, direction, longerstring;
   int right_rating=0, right_left=n, right_right=n-1, left_rating=0, left_left=n, left_right=n-1;
   int trunc_right, trunc_left;
   struct LastReportedRepeat *prev;

	char *forw_align1=(char*) malloc (length+1);    // need this much?
	char *align1=(char*) malloc (length*2+1);       // possibly length of string
	char *forw_align2=(char*) malloc (length+1);    // *6 total!
	char *align2=(char*) malloc (length*2+1);
	if (align1==NULL || align2==NULL || forw_align1==NULL || forw_align2==NULL)
   	{ printf("memory allocation error in OneIteration.");  exit(1); }


   char *rev_str = (char *)malloc(length+1);

   if (rev_str==NULL)

   	{ printf("memory allocation error in OneIteration.");  exit(1); }


   strcpy (rev_str, str);

   strrev (rev_str);


   mid_forward = str + n;

   mid_backward = rev_str+length-n;

   right_forward = mid_forward+MIN_PERIOD-1;

   left_forward = mid_forward-MIN_PERIOD+1;

   right_backward = mid_backward-MIN_PERIOD+1;

   left_backward = mid_backward+MIN_PERIOD-1;

   left.left=n;

   left.right=n-1;

   right.left=n;

   right.right=n-1;



  /***************************** OUTER LOOP *********************************/

	for(j=MIN_PERIOD;j<=n-1+length%2;j++) {
    for (direction=0;direction<2;direction++){
    	if (j==n&&direction==0)
      	continue;
    	else if (direction==0){
      	left_forward--;
    		left_backward++;
    		prev = &left;
    		mid_text = left_forward;
    		rev_text = left_backward;
    		mid_pattern = mid_forward;
    		rev_pattern = mid_backward;
      }
      else {
    		right_forward++;
    		right_backward--;
    		prev = &right;
         mid_text = mid_forward;
         rev_text = mid_backward;
         mid_pattern = right_forward;
         rev_pattern = right_backward;
      }

      diff = j - prev->j;

      //shift = MAX_ERRORS<diff?MAX_ERRORS:diff;

      // FILLS UP MATRIX_FORWARD, MATRIX_REVERSE, RIGHT_ERRORS, DOWN_ERRORS
      buildmatrix(matrix_forward,mid_pattern,mid_text,MAX_ERRORS,j);
      buildmatrix(matrix_reverse,rev_pattern,rev_text,MAX_ERRORS,-j);
		errors_array(right_errors,matrix_forward,MAX_ERRORS,1);
		errors_array(down_errors,matrix_reverse,MAX_ERRORS,-1);

      // BETWEEN ITERATIONS FILTER
//      if (right_errors[MAX_ERRORS].row + down_errors[MAX_ERRORS].col + j >= length && length<=2*MAX_PERIOD){
//      	if (offset!=START_POS && !isend) {/*fprintf(outfile,"Between Iterations Filter\n"); */return 0;}
//         if (offset!=START_POS || !isend){
//         	for (i=0;i<=MAX_ERRORS;i++){
//         		if (right_errors[i].row + down_errors[MAX_ERRORS-i].col + j >= length) {/*fprintf(outfile,"Between Iterations Filter\n"); */return 0;}
//         	}
//         }
//      }

      // FIGURES OUT THE RANGE OF THE REPEATS IN THE CURRENT J
      curr_left=n;
      curr_right=n-1;
      left_errs=0;
      right_errs=0;
 //     int prev_report=-10;
      for(i=MAX_ERRORS;i>=0;i--) {
      	int k1=MAX_ERRORS-i;
         if(i>0 && right_errors[i-1].row >= right_errors[i].row-1)
         	continue;
		  	while(k1>0 && down_errors[k1-1].col >= down_errors[k1].col-1)
         	k1--;
 //		  	if (prev_report >= k1)
 //				continue;
         if(right_errors[i].col + down_errors[k1].row >= j){
         	if (n + j*direction + right_errors[i].row-1 > curr_right){
            	curr_right = n + j*direction + right_errors[i].row-1;
               right_errs=i;
            }
            if (n + j*(direction-1) - down_errors[k1].col < curr_left){
            	curr_left = n + j*(direction-1) - down_errors[k1].col;
               left_errs=k1;
            }
         }
      }

      // IF NO REPEAT IS DETECTED, GO STRAIGHT TO CHECKING IF A PREVIOUS
      // REPEAT SHOULD BE REPORTED
      if (curr_left > curr_right){
      	comparison = -1;
         goto REPORT_PREVIOUS;
      }

      // REMOVES POORLY-MATCHING PORTIONS (ACCORDING TO ERROR_VAL) FROM EACH END
      // OF THE POTENTIAL REPEAT
      for (i=right_errs-1,trunc_right=right_errs; i>=0; i--){
      	if ( right_errors[trunc_right].row-right_errors[i].row < (ERROR_VAL+1)*(trunc_right-i) )
         	trunc_right=i;
      }
      for (i=left_errs-1,trunc_left=left_errs; i>=0; i--){
      	if ( down_errors[trunc_left].col-down_errors[i].col <= (ERROR_VAL+1)*(trunc_left-i) )
         	trunc_left=i;
      }

      // IF THE REPEAT WAS REDUCED TO UNDER 2 FULL PERIODS BY THE PREVIOUS
      // FILTER, IT IS EXTENDED TO AT LEAST 2 FULL PERIODS WITH THE LEAST
      // POSSIBLE ERRORS
      if (right_errors[trunc_right].col + down_errors[trunc_left].row < j){
      	int best_left, best_right, best_total;
         for (i=trunc_right;i<=right_errs;i++){
         	if (right_errors[i].col + down_errors[trunc_left].row >= j) break;
         }
         for (k=trunc_left;k<=left_errs;k++){
            if (right_errors[i].col + down_errors[k].row >= j) break;
         }
         best_left = k;
         best_right = i;
         best_total = k+i;
      	for (i--;i>trunc_right;i--){
         	while (k<left_errs && right_errors[i].col+down_errors[k].row<j)
            	k++;
            if (right_errors[i].col+down_errors[k].row<j) break;
            if (k+i < best_total){
            	best_left = k;
               best_right = i;
               best_total = k+i;
            }
         }

			trunc_right = best_right;
         trunc_left = best_left;
      }
                                            
      right_errs = trunc_right;
      left_errs = trunc_left;
      curr_right = n + j*direction + right_errors[right_errs].row-1;
      curr_left = n + j*(direction-1) - down_errors[left_errs].col;
      curr_off = n + j*direction - down_errors[left_errs].row;


      // GIVES THE REPEAT A RATING
      curr_errs = left_errs + right_errs;
      longerstring = down_errors[left_errs].col+right_errors[right_errs].col > down_errors[left_errs].row+right_errors[right_errs].row ? down_errors[left_errs].col+right_errors[right_errs].col : down_errors[left_errs].row+right_errors[right_errs].row;
      curr_rating = longerstring - curr_errs * (ERROR_VAL+1);
      curr_length = curr_right - curr_left + 1;

      // Comparisons: -1 is <, 0 is =, 1 is >
      if (curr_rating > prev->rating)
      	comparison=1;
      else if (prev->rating > curr_rating)
      	comparison=-1;
      else
      	comparison=0;


REPORT_PREVIOUS:

      // DETERMINES WHEN A (PREVIOUSLY SAVED) REPEAT SHOULD BE REPORTED AND
      // REPORTS IT
   	if ( (prev->newhigh==1) && (comparison==-1) ){
      	prev->newhigh=0;
         if ( (direction==0&&(prev->rating>right_rating||prev->left<=right_left-SHIFT||prev->right>=right_right+SHIFT)) || (direction==1&&(prev->rating>left_rating||prev->left<=left_left-SHIFT||prev->right>=left_right+SHIFT)) ){
            if (direction==0){
            	left_rating=prev->rating;
               left_left=prev->left;
               left_right=prev->right;
            }
            else{
            	right_rating=prev->rating;
               right_left=prev->left;
               right_right=prev->right;
            }

         	begpos = prev->left;
         	endpos = prev->right;

				// now print repeat as an alignment. align1 begins at n/2
				// and align2 begins at n/2+j
            if (direction)
		  			traceback(prev->matrix_forward,MAX_ERRORS,mid_pattern-diff,mid_text,forw_align2,forw_align1,prev->right_errors[prev->right_errs].row-1,prev->right_errors[prev->right_errs].col-1);
            else
               traceback(prev->matrix_forward,MAX_ERRORS,mid_pattern,mid_text+diff,forw_align2,forw_align1,prev->right_errors[prev->right_errs].row-1,prev->right_errors[prev->right_errs].col-1);
				strrev(forw_align1);
				strrev(forw_align2);

				// then, traceback the reverse matrix, concatenate alignment strings
				// and call printrepeat on the 2-sequence alignment.
            if (direction)
					traceback(prev->matrix_reverse,MAX_ERRORS,rev_pattern+diff,rev_text,align2,align1,prev->down_errors[prev->left_errs].row-1,prev->down_errors[prev->left_errs].col-1);
            else
            	traceback(prev->matrix_reverse,MAX_ERRORS,rev_pattern,rev_text-diff,align2,align1,prev->down_errors[prev->left_errs].row-1,prev->down_errors[prev->left_errs].col-1);

            // In case the repeat does not reach the midpoint n, we can take a
            // part of only one side

            // PSEUDOCODE

            // IF SIDE A > MINLENGTH AND SIDE B <MINLENGTH
            // FROM BEGINNING OF SIDE A CLOSEST TO CENTER TO END OF SIDE A
            // IF DELETING FROM END OF SIDE B TO CURRENT POINT WOULD RAISE
            // SCORE AND REPEAT WOULD REMAIN OVER MIN LENGTH, AND REPEAT WOULD REMAIN 2 FULL PERIODS, DO IT

				strcat(align1,forw_align1);
            strcat(align2,forw_align2);


				// PRINTS TO OUTPUT FILE IN SPECIAL FORMAT FOR POST-PROCESSING

            fprintf(outfile,"%ld %ld %d\n",begpos+offset,endpos+offset,prev->right_errs+prev->left_errs);
            printcompact(align1, align2, (long int)offset+begpos, (long int)offset+prev->offset, outfile);

 	 	  	}
		}


      // SAVES CURRENT VALUES IN "PREV" IF THIS J IS A CANDIDATE TO BE REPORTED
      if (comparison>0 && curr_length>=MIN_LENGTH && curr_rating>=MIN_RATING){
      	prev->newhigh=1;
			prev->left = curr_left;
      	prev->right = curr_right;
      	prev->right_errs = right_errs;
      	prev->left_errs = left_errs;
      	prev->length = curr_length;
         prev->rating = curr_rating;
      	prev->j=j;
         prev->offset = curr_off;
      	int** temp = matrix_forward;          // switches current with previous
      	matrix_forward = prev->matrix_forward;
      	prev->matrix_forward = temp;
      	temp = matrix_reverse;
      	matrix_reverse = prev->matrix_reverse;
      	prev->matrix_reverse = temp;
      	Entry * etemp = down_errors;          // switches current with previous
      	down_errors = prev->down_errors;
      	prev->down_errors = etemp;
      	etemp = right_errors;
      	right_errors = prev->right_errors;
      	prev->right_errors = etemp;
      }  // endif
    } // end for one direction
	} // end for j


  /***************************** OUTER LOOP END *****************************/

	free(align1);
	free(align2);
	free(forw_align1);
	free(forw_align2);
	free(rev_str);
	return 0;
}


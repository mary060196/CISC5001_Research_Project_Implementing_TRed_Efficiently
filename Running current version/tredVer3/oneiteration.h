//   |=================================================================|
//   |   TRED: a tool for detecting Tandem Repeats within sequences,   |
//   |               using the Edit Distance metric.                   |
//   |=================================================================|

// Copyright © 2007-2009 Dina Sokol, Justin Tojeira
// Distributed under the Aladdin Free Public License

// THIS SOFTWARE SHOULD BE ACCOMPANIED BY readme.txt AND license.html
// WE STRONGLY ENCOURAGE YOU TO READ BOTH BEFORE PROCEEDING
//------------------------------------------------------------------------------



#ifndef ONEITERATION_H
#define ONEITERATION_H

int OneIteration(char*, int**, int**, Entry*, Entry*, struct LastReportedRepeat, struct LastReportedRepeat, long int, bool, FILE*, int);

struct LastReportedRepeat{
	int j, left, right, length, rating;
   int left_errs, right_errs, total_errs;
   int ** matrix_forward, ** matrix_reverse;
   int offset;
   Entry *down_errors, *right_errors;
   bool newhigh;
};


#endif


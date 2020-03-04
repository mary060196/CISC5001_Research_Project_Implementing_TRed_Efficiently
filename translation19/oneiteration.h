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

// Copied from old 'errorsarray.h' by Miriam Briskman, 02.23.2020
typedef struct entry {
	int row;
	int col;
} Entry;

void OneIteration(char*, long, int**, int**, Entry*, Entry*, struct LastReportedRepeat, 
                  struct LastReportedRepeat, long int, bool, FILE*);

struct LastReportedRepeat{
    int   j, left, right, length, rating,
          left_errs, right_errs, total_errs,
          ** matrix_forward, ** matrix_reverse,
          offset;
    Entry *down_errors, *right_errors;
    bool  newhigh;
};

#endif

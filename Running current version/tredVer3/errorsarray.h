//   |=================================================================|
//   |   TRED: a tool for detecting Tandem Repeats within sequences,   |
//   |               using the Edit Distance metric.                   |
//   |=================================================================|

// Copyright © 2007-2009 Dina Sokol, Justin Tojeira
// Distributed under the Aladdin Free Public License

// THIS SOFTWARE SHOULD BE ACCOMPANIED BY readme.txt AND license.html
// WE STRONGLY ENCOURAGE YOU TO READ BOTH BEFORE PROCEEDING
//------------------------------------------------------------------------------



#ifndef ERRORSARRAY_H
#define ERRORSARRAY_H

#include <stdio.h>
#include <stdlib.h>

typedef struct entry {
	int row;
	int col;
} Entry;

void errors_array(Entry*, int**, int, int);

#endif


//   |=================================================================|
//   |   TRED: a tool for detecting Tandem Repeats within sequences,   |
//   |               using the Edit Distance metric.                   |
//   |=================================================================|

// Copyright © 2007-2009 Dina Sokol, Justin Tojeira
// Distributed under the Aladdin Free Public License

// THIS SOFTWARE SHOULD BE ACCOMPANIED BY readme.txt AND license.html
// WE STRONGLY ENCOURAGE YOU TO READ BOTH BEFORE PROCEEDING
//------------------------------------------------------------------------------

#ifndef BUILDK_H
#define BUILDK_H

void buildmatrixforward (int**, const char*, const int, const char*, const int, const int);
void buildmatrixbackward (int**, const char*, const int, const char*, const int, const int);

// External variables used in 'main.cpp' and 'buildk.cpp' (Miriam Briskman, 03.03.2020)
extern int *init_increments, *upper, *lower;

#define UNKNOWN 'N'

#endif

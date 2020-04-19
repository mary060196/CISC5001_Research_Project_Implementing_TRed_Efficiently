//   |=================================================================|
//   |   TRED: a tool for detecting Tandem Repeats within sequences,   |
//   |               using the Edit Distance metric.                   |
//   |=================================================================|

// Copyright © 2007-2009 Dina Sokol, Justin Tojeira
// Distributed under the Aladdin Free Public License

// THIS SOFTWARE SHOULD BE ACCOMPANIED BY readme.txt AND license.html
// WE STRONGLY ENCOURAGE YOU TO READ BOTH BEFORE PROCEEDING
//------------------------------------------------------------------------------

#ifndef TRACEK_H
#define TRACEK_H

#include <stdio.h>

void tracealign (int**, char*, char*, int, int, int, long&, long&, bool&, FILE*); // Miriam Briskman, 02.20.2020
void traceforward (int**, char*, char*, int, int, int, long, long, bool, FILE*); // Miriam Briskman, 02.20.2020

// External variables used in 'main.cpp' and 'tracek.cpp' (Miriam Briskman, 03.04.2020)
extern char *forward_pattern1, *forward_pattern2;

#endif

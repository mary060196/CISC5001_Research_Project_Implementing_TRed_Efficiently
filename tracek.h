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

// Miriam Briskman, 02.20.2020:
void tracealign (int**, char*, char*, int, int, int, long&, long&, bool&, FILE*); 
void traceforward (int**, char*, char*, int, int, int, long, long, bool, char*, char*, FILE*);

#endif

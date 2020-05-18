#ifndef NEW_MATRIXL_H
#define NEW_MATRIXL_H

#include "parameters.h"

// Constants
const long long JOINED_LEN = MAX_PERIOD*6 + 10;

// Structures:
struct suffix {
    long long text_len, index;
};

// Function Prototypes:
void prepareSuffixArrays (char*, const long long, const char*, const int, const char*, const int, const int,
                          int*, int*, int*, suffix**, bool*, int*[], bool);
void stringIncrement (char*, const long long, const int, int*, int*, int*, suffix**, bool*, int*[]);
void buildMatrixL (int**, int**, const long long, int*, int*, int*[], const int, const int, int, int);
void printSuffixes (char*, const long long, suffix**);

#endif // NEW_MATRIXL_H

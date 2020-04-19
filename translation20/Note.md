# Translation 20 (April 19, 2020)

## To compile the main program, enter:

    g++ -std=c++11 -o tred main.cpp oneiteration.cpp buildk.cpp tracek.cpp

## Changes from previous translation (Translation 19):

- A new loop was added; the loop skips substrings of length `2*MAX_PERIOD` filled entirely with `UNKNOWN` (`N`) characters. The idea behind this is that substrings with `UNKNOWN` characters do no produce any tandem repeats, so we can skip large areas thereof without loss.

- In the `buildmatrixforward` and `buildmatrixbackward` functions, no element is initialized to `-1` anymore, as the code of the program checks when those special regions are entered, without the need to check if any of the elements is actually equal to `-1`. Thus, the code in `buildk.cpp` focuses only on building the meaningful regions of the `KxK` matrix.

- In the `tracealign` function in `tracek.cpp`:

  - `pos1` and `pos2` were prevented from incrementing in every iteration of the loop. Instead, the position values are computed when needed and at the end of the loop.

  - The unused variable `oPos` was removed.

  - The function accepts a new argument `int lowerRow` to check if we do not enter into a `-1`s region in the `KxK` matrix.

  - Only one extra condition `row >= lowerRow` was added to the function to eliminate entrance into regions in the `KxK` matrix that once contained `-1`s.

  - The condition
  
        else if (*(s1 - s1pos) == *(s2 - s2pos) && *(s1 - s1pos) != UNKNOWN)
        
    was changed to
    
        else
        
    since the cases in this `if-else-if` branch are all mutually exclusive, hence preventing the need for a trailing `else` statement after it.

- In the `traceforward` function in `tracek.cpp`:

  - The unneeded `int size` argument accepted by the function was removed.

  - The function accepts a new argument `int upperRow` to check if we do not enter into a `-1`s region in the `KxK` matrix.

  - Only one extra condition `row <= upperRow` was added to the function to eliminate entrance into regions in the `KxK` matrix that once contained `-1`s.

  - The condition
  
        else if (*(s1 - s1pos) == *(s2 - s2pos) && *(s1 - s1pos) != UNKNOWN)
        
    was changed to
    
        else
        
    since the cases in this `if-else-if` branch are all mutually exclusive, thus preventing the need for a trailing `else` statement after it.

- In the `OneIteration` function:

  - The unused argument (bool isend) was removed.

  - The calls to `buildmatrixforward` and `buildmatrixbackward` function were placed inside the `if-else` branch preceding them. This allowed to eliminate the use of the variables `rev_pattern`, `rev_text`, `mid_text,` and `mid_pattern`.

  - The algorithm for filling up the `right_errors` and `down_errors` arrays was changed such that it does not refer to the regions filled with `-1`s inside the KxK matrices.

  - The `LastReportedRepeat` structure defined in `oneiteration.h` was added two new fields: `NoMinusOnesUpperRowForward` and `NoMinusOnesLowerRowBackward`, which are being updated in `OneIteration` when the information about a candidate tandem repeat is saved (near the end of the function's definition.)

- Bug fixing: Users can enter the sequence and intermediary filenames *after* the program launches if the filenames were not entered as command line arguments.

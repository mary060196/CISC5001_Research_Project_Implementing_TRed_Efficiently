# Edit Distance - New Algorithms

## Description

The new algorithms for computing the Edit Distance matrix (also known as the `D` or `NxN` matrix) aspire to reduce the amount of computation by constant time (averagely, by 2.5 / 7 =~ 36% of the original work.)

## Directories

You would find 4 directories here:

1. `Expand_Edit_Distance`: computing expanding edit distance matrices from one another, with the string given in a forward direction.

2. `Reverse_Expand_Edit_Distance`: computing expanding edit distance matrices from one another, with the string given in a backward (reverse) direction.

3. `Shrink_Edit_Distance`: computing shrinking edit distance matrices from one another, with the string given in a forward direction.

2. `Reverse_Shrink_Edit_Distance`: computing shrinking edit distance matrices from one another, with the string given in a backward (reverse) direction.

The algorithms in `Expand_Edit_Distance` and `Reverse_Expand_Edit_Distance`, and the algorithms in `Shrink_Edit_Distance` and `Reverse_Shrink_Edit_Distance`, are the same except for the direction of the sequence.

## Other Notes

Please view the README files of the individual algorithms inside the respective directories for further information.
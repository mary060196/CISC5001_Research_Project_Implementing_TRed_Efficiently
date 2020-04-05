# Shrink Edit Distance

## Goal

Make the construction of the edit distance matrix (known also as `NxN` or `D`) faster, at least by constant time.

## How the Program Works

The program implements the new algorithm and tests for its correctness by comparing the obtained matrix to a matrix built using the conventional Edit Distance algorithm.

It builts 10,000 random strings of random sizes between 1 and 500 that serve as the 'text' and conducts the test upon the resulting matrices.

If an inaccuracy is detected, both the wrong and the correct matrix are displayed and the program haults. Otherwise, each of the constructed matrices is correct, and the statement

    Testing Successful!

is output to `stdout`. If the testing is successful, there is very little chance that an obtained matrix by the algorithm is incorrect.

## Getting Started

### Printing all the Created Matrices to a file

In the `shrink_edit_distance` C++ file, lines 161 and 191 are commented out. If those lines are uncommented, the obtained matrices during the testing will be printed to a file specified as an argument to the program (see below).

### To Compile

Use

    g++ -std=c++11 -o shrink_edit_distance shrink_edit_distance.cpp testing_shrink.cpp

to compile the program.

### To Run

Use

    shrink_edit_distance.exe <output data file>

on Windows and
    ./shrink_edit_distance <output data file>

on Mac / Linux,

where **<output data file>** should be replaced by a name of a file to which the constructed matrices should be printed to run the program. This command line argument must be passed to the program even if the printing-to-file lines remain commented out.

## Troubleshooting

If the program runs for too long, and no error messages nor the `Testing Successful!` line is printed, you might need to press ENTER in the terminal / command prompt to finish the program and have the `Testing Successful!` line printed out. This does not mean that a matrix obtained in the program is incorrect.


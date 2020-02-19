02.19.2020
			string_functions_testing

This program attempts to measure the time it takes to
1) Read a string of length n
2) Reverse it,
3) Concatenate it with "bla-bla",
4) Find its length, and
5) Change its content using iteration and dereferencing
using both C-string and C++ string class functions.

What the program does:
1) It asks the user to enter an integer n > 0 representing the length of a desired string. 
Note that the user is not expected to provide the actual string, since the program creates a demo sequence in a new file and operates with it.
2) As mentioned, the program opens a new file named 'newSeq.fsa' and writes a demo DNA sequence (with random bases A, C, G, T, and N) of the specified length n into it. The created file has a FASTA format and syntax.
3) The program then opens the file for reading, sets a time counter, and reads the sequence. Afterwards, the timer is stopped and the time it took for the reading process is computed and displayed on the screen (in milliseconds.)
4) Step 3 is repeated, but now the program reverses the string.
5) Step 3 is repeated, but now the program concatenates it with "bla-bla".
6) Step 3 is repeated, but now the program finds the size of the string.
7) Step 3 is repeated, but now the program iterates through the string and changes the values.
Steps 3-7 are conducted both on c-strings (char arrays) and C++ string objects, separately. Both C++ files (<fstream>) and C files (<stdio.h>) are used, separately, for comparison.

Hypothesis / Expectation:
It is anticipated that each of the char-array operations would take a time proportional to the length of the string; that is, each of the operations is of O(n). Moreover, since the running time is a function of n, the larger the integer n, the longer time it takes for the program to complete the operations. The reason is that the complexity of each of the three C functions used in the program is O(n).
On the other hand, the length measurement of a C++ function should take a constant time, while the other functions will be of linear running time.

Environment:
The program was run on a desktop computer having the following features:
1) Processor:      Intel(R) Core(TM) i5-2300 CPU @ 2.80GHz  2.80GHz
2) Installed RAM:  6.00 GB
3) System Type:    64-bit operating system, x64-based processor

Results:
The program was run 6 times on input of various sizes, and output the following results:
for n = 1000:
It took 0 milliseconds to create a FASTA file with 1000 bases.
It took 0 milliseconds to read the sequence in using 'fgetc' and "+=" into a string object.
It took 0 milliseconds to read the sequence in using 'fgetc' into a c-string.
It took 0 milliseconds to convert a c-string into a string object.
It took 0 milliseconds to read the sequence in using fstream 'get' and "+=" into a string object.
It took 0 milliseconds to read the sequence in using fstream 'get' into a c-string.
It took 0 milliseconds to convert a c-string into a string object.
It took 0 milliseconds to reverse the sequence using 'reverse' (<algorithm>).
It took 0 milliseconds to reverse the sequence using 'strrev'.
It took 0 milliseconds to concatenate the sequence with 'bla-bla' using '+='.
It took 0 milliseconds to concatenate the sequence with 'bla-bla' using 'strcat'.
The length of the string is: 1007 (including 7 chars of 'bla-bla'.)
It took 0 milliseconds to find the sequence's length using 'length'.
The length of the string is: 1007 (including 7 chars of 'bla-bla'.)
It took 0 milliseconds to find the sequence's length using 'strlen'.
It took 0 milliseconds to change the content of a string object to all 'A's.
It took 0 milliseconds to change the content of a c-string to all 'A's.

for n = 10000:
It took 1 milliseconds to create a FASTA file with 10000 bases.
It took 1 milliseconds to read the sequence in using 'fgetc' and "+=" into a string object.
It took 0 milliseconds to read the sequence in using 'fgetc' into a c-string.
It took 0 milliseconds to convert a c-string into a string object.
It took 0 milliseconds to read the sequence in using fstream 'get' and "+=" into a string object.
It took 0 milliseconds to read the sequence in using fstream 'get' into a c-string.
It took 0 milliseconds to convert a c-string into a string object.
It took 0 milliseconds to reverse the sequence using 'reverse' (<algorithm>).
It took 0 milliseconds to reverse the sequence using 'strrev'.
It took 0 milliseconds to concatenate the sequence with 'bla-bla' using '+='.
It took 0 milliseconds to concatenate the sequence with 'bla-bla' using 'strcat'.
The length of the string is: 10007 (including 7 chars of 'bla-bla'.)
It took 0 milliseconds to find the sequence's length using 'length'.
The length of the string is: 10007 (including 7 chars of 'bla-bla'.)
It took 0 milliseconds to find the sequence's length using 'strlen'.
It took 0 milliseconds to change the content of a string object to all 'A's.
It took 0 milliseconds to change the content of a c-string to all 'A's.

for n = 100000:
It took 12 milliseconds to create a FASTA file with 100000 bases.
It took 11 milliseconds to read the sequence in using 'fgetc' and "+=" into a string object.
It took 5 milliseconds to read the sequence in using 'fgetc' into a c-string.
It took 0 milliseconds to convert a c-string into a string object.
It took 12 milliseconds to read the sequence in using fstream 'get' and "+=" into a string object.
It took 6 milliseconds to read the sequence in using fstream 'get' into a c-string.
It took 0 milliseconds to convert a c-string into a string object.
It took 1 milliseconds to reverse the sequence using 'reverse' (<algorithm>).
It took 0 milliseconds to reverse the sequence using 'strrev'.
It took 0 milliseconds to concatenate the sequence with 'bla-bla' using '+='.
It took 0 milliseconds to concatenate the sequence with 'bla-bla' using 'strcat'.
The length of the string is: 100007 (including 7 chars of 'bla-bla'.)
It took 0 milliseconds to find the sequence's length using 'length'.
The length of the string is: 100007 (including 7 chars of 'bla-bla'.)
It took 0 milliseconds to find the sequence's length using 'strlen'.
It took 1 milliseconds to change the content of a string object to all 'A's.
It took 0 milliseconds to change the content of a c-string to all 'A's.

for n = 1000000:
It took 137 milliseconds to create a FASTA file with 1000000 bases.
It took 114 milliseconds to read the sequence in using 'fgetc' and "+=" into a string object.
It took 58 milliseconds to read the sequence in using 'fgetc' into a c-string.
It took 0 milliseconds to convert a c-string into a string object.
It took 131 milliseconds to read the sequence in using fstream 'get' and "+=" into a string object.
It took 67 milliseconds to read the sequence in using fstream 'get' into a c-string.
It took 0 milliseconds to convert a c-string into a string object.
It took 13 milliseconds to reverse the sequence using 'reverse' (<algorithm>).
It took 1 milliseconds to reverse the sequence using 'strrev'.
It took 0 milliseconds to concatenate the sequence with 'bla-bla' using '+='.
It took 0 milliseconds to concatenate the sequence with 'bla-bla' using 'strcat'.
The length of the string is: 1000007 (including 7 chars of 'bla-bla'.)
It took 0 milliseconds to find the sequence's length using 'length'.
The length of the string is: 1000007 (including 7 chars of 'bla-bla'.)
It took 0 milliseconds to find the sequence's length using 'strlen'.
It took 8 milliseconds to change the content of a string object to all 'A's.
It took 3 milliseconds to change the content of a c-string to all 'A's.

for n = 10000000:
It took 1322 milliseconds to create a FASTA file with 10000000 bases.
It took 1172 milliseconds to read the sequence in using 'fgetc' and "+=" into a string object.
It took 645 milliseconds to read the sequence in using 'fgetc' into a c-string.
It took 11 milliseconds to convert a c-string into a string object.
It took 1395 milliseconds to read the sequence in using fstream 'get' and "+=" into a string object.
It took 679 milliseconds to read the sequence in using fstream 'get' into a c-string.
It took 11 milliseconds to convert a c-string into a string object.
It took 140 milliseconds to reverse the sequence using 'reverse' (<algorithm>).
It took 15 milliseconds to reverse the sequence using 'strrev'.
It took 0 milliseconds to concatenate the sequence with 'bla-bla' using '+='.
It took 15 milliseconds to concatenate the sequence with 'bla-bla' using 'strcat'.
The length of the string is: 10000007 (including 7 chars of 'bla-bla'.)
It took 0 milliseconds to find the sequence's length using 'length'.
The length of the string is: 10000007 (including 7 chars of 'bla-bla'.)
It took 0 milliseconds to find the sequence's length using 'strlen'.
It took 93 milliseconds to change the content of a string object to all 'A's.
It took 38 milliseconds to change the content of a c-string to all 'A's.

for n = 100000000:
It took 13243 milliseconds to create a FASTA file with 100000000 bases.
It took 11631 milliseconds to read the sequence in using 'fgetc' and "+=" into a string object.
It took 6194 milliseconds to read the sequence in using 'fgetc' into a c-string.
It took 120 milliseconds to convert a c-string into a string object.
It took 13265 milliseconds to read the sequence in using fstream 'get' and "+=" into a string object.
It took 6815 milliseconds to read the sequence in using fstream 'get' into a c-string.
It took 120 milliseconds to convert a c-string into a string object.
It took 1451 milliseconds to reverse the sequence using 'reverse' (<algorithm>).
It took 118 milliseconds to reverse the sequence using 'strrev'.
It took 62 milliseconds to concatenate the sequence with 'bla-bla' using '+='.
It took 82 milliseconds to concatenate the sequence with 'bla-bla' using 'strcat'.
The length of the string is: 100000007 (including 7 chars of 'bla-bla'.)
It took 0 milliseconds to find the sequence's length using 'length'.
The length of the string is: 100000007 (including 7 chars of 'bla-bla'.)
It took 23 milliseconds to find the sequence's length using 'strlen'.
It took 891 milliseconds to change the content of a string object to all 'A's.
It took 299 milliseconds to change the content of a c-string to all 'A's.

Discussion:
1) As we conjectured, the running time of each function, besides 'length' of the C++ string, depended on the input size. 
2) To approximate, 1 ms is needed to read an input of size 10000 into a string object, either using 'fgetc' of <stdio.h> or 'get' of <iostream>. On the other hand, 1 ms is sufficient to read an input of size 20000, double than before, using C-strings. 
3) It takes 10 time longer to revert a string object than it does to revert a C-string. Concatenating C++ strings takes shorter time than it does concatenating C-string does.
4) 'strlen' for C-string is of complexity O(n), while 'length' of string object is of constant complexity.
5) Iterating through a string object is 4 times more expensive than doing so to a C-string.

Conclusions:
1) In terms of running time in general, C-strings defeat string objects.
2) It is better to use a string object to accomplish a work involving the measurement of the length of a string.
3) For relatively small input sizes, the difference in the running time is insignificant.
4) No major difference exists between reading from C-style files and from C++-style files.
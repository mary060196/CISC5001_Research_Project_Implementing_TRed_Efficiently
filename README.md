# TRed Version 4.0 
## A tool for detecting Tandem Repeats within sequences, using the Edit Distance metric.

May 17, 2020
							      
Copyright © 2007-2020 Dina Sokol, Justin Tojeira, Miriam Briskman.

Distributed under the Aladdin Free Public License.

See [`license.htm`](./License.htm) for details.

<hr>

# Documentation + Usage Instructions

## VERSION NOTES

Version 4.0 includes a post-processing combiner/filter which must be run after the main program. For details on this, see [RUNNING TRed 4.0](#running-tred-40).

<hr>

## TABLE OF CONTENTS

1) [RUNNING TRed 4.0](#running-tred-40)

2) [INPUT/OUTPUT](#inputoutput)

3) [PARAMETERS](#parameters)

4) [SOFTWARE AND HARDWARE REQUIREMENTS](#software-and-hardware-requirements)

5) [EXAMPLE OF RUNNING PROGRAMS](#example-of-running-programs)

<hr>

## RUNNING TRed 4.0:

### Running the `main` Program

To run TRed, you must first run the main program, then run one of two post-processing programs. The main program is headed by `main.cpp` and includes all `.cpp` and `.h` files with the exception of `filter.cpp` and `nofilter.cpp`. If you use UNIX, this may be compiled via the makefile. The files `filter.cpp` and `nofilter.cpp` must be compiled and run separately.

To compile the main program on a non-UNIX or Linux computer, enter:

    g++ -g -std=c++11 -o tred main.cpp aux_main.cpp oneiteration.cpp buildk.cpp tracek.cpp

inside the terminal / command line application.

Using any compilation, the executable's name (excluding the extension) will be `tred`.

The main program requires two filenames to be typed in at the command line. The first we will call the "sequence file". This will contain the sequence (nucleotide, protein, etc.) that you wish to analyze. 

**NEW!** TRed 4.0 Supports FASTA formatted sequences. Please refer to the [INPUT/OUTPUT](#inputoutput) section below for details.

The second filename that the program requires is called the "intermediary file". This file will be the output of the main program, and will be formatted to be read by the post-processing programs.

Hence, to run the main program, enter:

    tred.exe <sequence file's name> <intermediary file's name>

on non-UNIX/Linux systems, or:

    ./tred <sequence file's name> <intermediary file's name>

on UNIX/Linux systems.

**NOTE:** If the "intermediary file" whose name is being passed to the program does not already exist in the directory where the executable resides, the program will create a new file in the same directory with that name. If the file already exists in that directory, all its content **WILL BE ERASED PERMANENTLY** and the file will be overwritten.

### Running the Post-Processing Programs (`filter` and `nofilter`)

There are two post-processing programs, `filter` and `nofilter`. `filter` will combine the repeat fragments found in the main program, and discard any repeats judged to be extraneous. `nofilter` will report all repeat fragments found in the main program, without combining. 

To compile the `filter`, enter:

    g++ -g -std=c++11 -o filter filter.cpp

inside the terminal / command line application.

To compile the `nofilter`, enter:

    g++ -g -std=c++11 -o nofilter nofilter.cpp

inside the terminal / command line application.

After compilation, the executable's name (excluding the extension) will be `filter` or `nofilter`, respectively.

Each of `filter` or `nofilter` takes 4 command line parameters:

1. The sequence file (the same one you used for the main program)
2. The intermediate file (the output from the main program)
3. The alignment file, and
4. The output table file. The last two will contain the final output from the TRed program. See [INPUT/OUTPUT](#inputoutput) for details.

Hence, to run `filter`, enter:

    filter.exe <sequence file's name> <intermediary file's name> <alignment file's name> <table file's name>

on non-UNIX/Linux systems, or:

    ./filter <sequence file's name> <intermediary file's name> <alignment file's name> <table file's name>

on UNIX/Linux systems. To run `nofilter`, replace the instances of `filter` above by `nofilter`.

**NOTE:** `filter` and `nofilter` are now programmed to post-process output generated based on FASTA files.

<hr>

## INPUT/OUTPUT

Input (the sequence file) must be one file containing one sequence.

**NEW!** TRed 4.0 Supports FASTA formatted sequences. For the program to recognize the sequence, the following two conditions about the format must hold:

1) Each FASTA comment is assumed to begin with a `>` character.
2) Each FASTA comment is assumed to be followed by at least one `\n` (newline) character.

The program also supports sequence files containing whitespace characters from the list below:

1) `\n` (newline)
2) `\t` (tab)
3) `\r` (carriage return)
4) `\f` (formfeed page break)
5) `\v` (vertical tab), and
6) ` `  (space).

These whitespaces can appear at any location in the sequence file. The program is able to recognize those and filter them out, so that they do not become a part of a tandem repeat. Clearly, the program is able to process a sequence with no whitespaces or FASTA comments.

Output consists of two files. The first, which we will call "output table", is a summary of the detected repeats in table form. The second, which we will call "alignment file", shows the full alignment of the repeats. 

**NOTE:** If any of the "output table" or "alignment file" whose names are being passed to the program does not already exist in the directory where the executable `filter` or `nofilter` resides, the program will create a new file in the same directory with that name. If the file already exists in that directory, all its content **WILL BE ERASED PERMANENTLY**, and the file will be overwritten.

The two output files show the following information:

Output Table:
 - `StartIndx`: The index of the first character in the repeat
 - `EndIndex`: The index of the last character in the repeat
 - `Length`: The length of the repeat
 - `Period`: The average length of a period of the repeat, as the period of an evolutive tandem repeat may vary in length throughout the repeat.
 - `Reps`: The number of periods in the repeat
 - `Errs`: The number of errors (insertions, deletions and mismatches) between consecutive periods of the repeat
 - `Score`: This is the rating of the repeat, which is an approximation of the difference between the number of matches in the repeat, and the number of errors in the repeat multiplied by the parameter `ERROR_VAL`.

Alignment File:

The repeats in the alignment file appear in the same order as they do in the output table. Each alignment is preceded by a line stating the length of the repeat, the start and end indices of the repeat, and the number of errors in the repeat, and is followed by its rating. Each period of the repeat is shown as its own line in the alignment, with the first and last indices of that period. Periods are aligned vertically, with each character matching against the character above it (from the previous period) and below it (from the next period). Lines of the alignment which have no indices are the same period as the line above it, shown to align with the line below it.

**Note:** alignment file must be in monospace font in order to properly line up. A change in the number of digits in the indices (for example, from 999 to 1000) may also disrupt the proper alignment.

<hr>

## PARAMETERS

Understanding the following parameters is crucial in order to be able to use TRed effectively. More restrictive parameters will result in the program running faster. These parameters are being used in all of the main, `filter` and `nofilter` programs.

 - `MAX_ERRORS`: Repeats with more than this amount of errors between periods may not be detected, meaning `MAX_ERRORS` may have to be raised to detect divergent repeats with large periods. Also, segments of the repeat must meet the minimum length (`MIN_LENGTH`) and minimum rating (`MIN_RATING`) requirements in the allowed amount of errors, meaning that even when searching for repeats of very short period, you may have to allow a few errors. Raising this value will make the program run much slower, lowering it will make the program run much faster. If this value is too low, repeats may appear fragmented, with only the better-matching portions of a repeat being reported. Default value is 10. Must be 0 or greater. Huge effect on speed of the program.

 - `MIN_LENGTH`: Repeats with length less than `MIN_LENGTH` will not be reported, even if they are detected. Repeats with length greater than or equal to `MIN_LENGTH` are guaranteed to be detected, and will be reported. Increasing this value will make the program run a little faster, decreasing it will make the program run a little slower. Default value is 20. Must be 1 or greater. Recommended that this is less than `MAX_PERIOD`.

 - `MIN_RATING`: Repeats with a rating (also called score) less than this value will not be reported. The rating is an approximation of the number of matches in the repeat minus the number of errors in the repeat multiplied by `ERROR_VAL`. In other words, each match adds 1 to the rating, while each error subtracts from the rating an amount equal to `ERROR_VAL`. Thus, in order for a repeat to be reported it needs `ERROR_VAL` matches for each error, plus an additional `MIN_RATING` matches. Aside from influencing the amount of output the program gives, this value has a negligible effect on the speed of the program. Default value is 15.

 - `MIN_PERIOD`: The period of an evolutive tandem repeat may change throughout the repeat. Thus, only repeats whose periods are larger than or equal to this value throughout the repeat are guaranteed to be detected. Repeats whose periods are always smaller than `MIN_PERIOD` will not be detected. Repeats whose periods are sometimes larger and sometimes smaller may or may not be detected. Increasing this value will make the program run a little faster. Default value is 1. Must be 1 or greater.

 - `MAX_PERIOD`: Large inputs are processed by breaking them up into pieces of size `2*MAX_PERIOD`. Repeats with periods up to the size of `MAX_PERIOD` will be detected. Repeats with period sizes up to `2*MAX_PERIOD` may or may not be detected. Repeats with period sizes greater than `2*MAX_PERIOD` will not be detected. Increasing this value will make the program run a little slower, decreasing it will make the program run a little faster. Default value is 250. Must be at least 1/4 the value of `MIN_LENGTH`. Recommended that this is greater than `MIN_LENGTH`.

 - `ERROR_VAL`: This is how much an error counts against the rating. See `MIN_RATING` for details. When modifying this parameter, keep in mind that for a repeat to be reported, it will need to have a number of matches equal to the value of `ERROR_VAL` for each error it has, making it `(1-100/(ERROR_VAL+1))` percent matching. So a value of 1 is `50%` matching, 2 is `66%` matching, 3 is `75%` matching, etc. Aside from influencing the amount of output the program gives, this value has a negligible effect on the speed of the program. Default value is 3. Recommended value is 1 or greater.

 - `SHIFT`: Repeats found in the same iteration will be reported if they either have a better rating than previously reported repeats or if they extend `SHIFT` characters to the left or right. Raising this value may help get rid of redundant output (slightly). My recommendation is not to bother. Default value is 3. Strongly recommended value is 1 or greater.

 - `START_POS`: This is the index of the first character of your input string. Normally the first index will be 1, but in certain cases it may not be. For example, if you wish to examine only the middle part of a sequence, say from 1500 to 2500, you would paste the part you wish to run the program on into a new file, and set `START_POS` to 1500. Default value is 1.

 - `PROCESSING`: During execution, the program will repeatedly print `Processing.....` to the screen show the program's progress. Each line represents `2*MAX_PERIOD` base pairs, so this can be used to gauge the speed of the program using different parameters. To shut off this feature, set `PROCESSING` to 0. To turn it on, set `PROCESSING` to 1. Default value is 0.

 - `UNKNOWN`: This is found in `buildk.h`. It designates a character to represent an unknown in the sequence. This character, when detected in the input string, will not be considered a match to any other character, not even itself. For example, in DNA strings, the character `N` is used to denote an unknown base in the sequence. Because we don't know what it is, we do not consider it to match anything. When changing this value, the unknown character must be enclosed in single quotes. For example, `'N'` or `'n'` or `'X'` (uppercase or lowercase matters). If there is no unknown character in your sequence, you may set this to a character which does not appear in your sequence, or you may set it to `-1`. If you set it to `-1`, do not use quotes. Default character is `'N'`.

<hr>

## SOFTWARE AND HARDWARE REQUIREMENTS

### Software Requirements

- The software was tested both on Windows 10 and Ubuntu 18.04 (Linux). It is expected to operate correctly on any Windows or UNIX/Linux operating system.

- A `C++` compiler of version `C++11` and up is required to compile the programs.

### Hardware Requirements

- To run the main program, a disk space of two (2) times the input sequence file is required.

- It is **UNKNOWN** ahead of time how much space the 'alignment' output file, whose content is produced by either of `filter` or `nofilter` programs, will require. The reason is that the size of the 'alignment' output file depends on the detected repeats. For example, when running the `nofilter` program on Chromosome Y (27 MB size,) which is the shortest chromosome, we manually terminated the `nofilter` program when about 1/8 of the job was completed. At the point after the termination, the 'alignment' file reached the size of 8GB. Hence, it is unknown what size the 'alignment' file for a larger sequence will reach when the `nofilter` or `filter` programs complete their jobs. It is suggested, hence, to run the programs on portions of the sequences once at a time.

- The 'table' file will reach a size not larger than the size of the sequence.

<hr>

## EXAMPLE OF RUNNING PROGRAMS

In the directory, you will find the files:

- `newSeq.txt`,
- `intermed.txt`,
- `align.txt`, and
- `table.txt`.

The file `newSeq.txt` contains the first 33K bases (or so) of Chromosome Y, including the first FASTA comment. We ran the main program:

    tred.exe newSeq.txt intermed.txt

and then the `nofilter` program:

    nofilter.exe newSeq.txt intermed.txt align.txt table.txt

and obtained what the `intermed.txt`, `align.txt`, and `table.txt` files contain.

<hr>

For **questions and comments** regarding the software and accompanying documentation, contact [jtojeira@gmail.com](mailto:jtojeira@gmail.com)

For updates, [http://www.sci.brooklyn.cuny.edu/~sokol/trepeats/index.html](http://www.sci.brooklyn.cuny.edu/~sokol/trepeats/index.html)
For all other queries, contact [sokol@sci.brooklyn.cuny.edu](mailto:sokol@sci.brooklyn.cuny.edu)

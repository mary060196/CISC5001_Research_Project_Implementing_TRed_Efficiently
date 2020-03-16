# Study's Progress Diary
This file keeps track of how the research project advances. Mentions of milestones and important notes are introduced.

### Sunday, March 15, 2020
- The original TRed Version 3 and the translated code in [translation19](https://github.com/mary060196/CISC5001_Research_Project_Implementing_TRed_Efficiently/tree/master/translation19) were run, separately, on the entire Chromosome Y fasta file.
- The fasta file, together with all the other Chromosome sequence files, was earlier downloaded from

      [ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/ARCHIVE/ANNOTATION_RELEASE.109/](ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/ARCHIVE/ANNOTATION_RELEASE.109/)

- The original program ran for 8420086 miliseconds (= 2.34 hours,) while the translated code ran for 6800100 miliseconds (= 1.89 hours.)
- This is a decrease of **19.24%** in the total running time.
- The programs were run on a desktop computer having the following specifications:

          Processor:      Intel(R) Core(TM) i5-2300 CPU @ 2.80GHz  2.80GHz
          Installed RAM:  6.00 GB
          System Type:    64-bit operating system, x64-based processor
          
- It is expected that the decrease when running the programs on larger chromosome sequences (e.g., Chromosome 1 or 2,) will be greater than **20%**.

### Wednesday, March 4, 2020
- Just a periodic update of how the semester advances:
- All the logical errors in the program detected by now were eliminated. This is confirmed by running the program against a relatively long biological sequence and checking if the output of the program with modifications has **the exact same output** as the one obtained by the original TRed version 3.
- No c-string functions are used in the program. Instead, strings are accessed using `char` pointers and integer variables keeping track of the strings' length.
- Repetitive code with similar structure was combined. For example, two consecutive `for` loops, both iterating from some fixed `i` to some fixed `i+k`, were combined into a single loop running on the same indices, as long as their functionality is not harmed.
- Dates were added across the program next to new lines of code.
- Indentation across the program was fixed.
- The `makefile` has been updated.
- Memory of a **fixed and known ahead of time** size is now allocated **only once** inside `main.cpp`. Usually, the size is given by one or more of the macros in `parameters.h`. External variables are used instead of function arguments to access these fixed data fields.
- Declaration of variables takes place at the beginning of functions rather than scatteredly in the middle.
- The total time it takes the program to execute is measured and displayed to `stdout`.
- The latest update to the program is located in the directory `translation19` which can be found on the tandem repeats server and in this GitHub repository.
- Lastly, the overall running time of the main program on a 33,000 byte size file, with a parameter specification of 

                     #define MAX_ERRORS 20
                     #define MIN_LENGTH 5
                     #define MIN_RATING 15
                     #define MIN_PERIOD 1
                     #define MAX_PERIOD 250
                     #define ERROR_VAL 3
                     #define SHIFT 3
                     #define START_POS 1
                     #define PROCESSING 0
                     
takes now 9.5 seconds instead of 10.5 seconds it takes for the original TRed version 3. The file on which the program was run, `newSeq.txt`, is located in the `translation19` directory in this repository.
- At this point, we move into a new phase in this semester, during which we work on suffix arrays and their inclusion into the program.

### Wednesday, February 26, 2020
- 4th meeting with Professor Sokol in the Spring 2020 Semester.
- By now, the program reads the entire bio-sequence into a single buffer, whose size was computed using the standard library `fseek` and `ftell`.
- The concatenation is, therefore, not necessary anymore, so the corresponding `strcat` functions were eliminated.
- An attempt was made to replace string reveral by pointer manipulation. However, the program produces an error related to this change.
- During the meeting, Professor Sokol approved the current implementation of `main.cpp`.
- Nonetheless, the naming of variables in `main.cpp` and in other files of the program is quite unclear. Variables should carry mnemonic names.
- Professor Sokol explained the idea behind the generation of the `kxk` matrix and how it is achieved.
- One important aspect that deserves focus is how many times the macro `FILL_KXK_MATRIX` is invoked during the running of the program on some fixed sequence. This will allow us decide how helpful the integration of suffix arrays into the program will.
- It would also be useful to compare the running time of the current version of TRed with other, similar tandem repeat software that works similar to TRed to understand whether there might exist a faster way to implement the program.
- Finally, memory in the program should be allocated dynamically to eliminate stack overflow scenarios and allow memory allocation of sizes unknown to the compiler before runtime.
- Below is the section from `tasks.txt` pertaining to today:

@@@@@

Wednesday, 02.26.2020 (Meeting w/ Professor Sokol)
Short Term:
- Keep on debugging the program to allow pointer manipulation to replace string reversal methods. (Done)
- Change the names of variables so that they imply about the purpose of the variables. (Done)
- Keep on working on proper indentation of the code to ease readability. (Done)
- After all the key short-term goals (including those above) are completed, begin reading about suffix arrays and thinking how to integrate them into the program, hoping to achieve greater efficiency.

@@@@@

### Wednesday, February 19, 2020
- 3rd meeting with Professor Sokol in the Spring 2020 Semester.
- The essence of the `main` source file was discussed; in particular, the concern was about how the seqeunce is read into the file.
- It is imperative to remove all the "linear time" functions from the program, including the two string concatenation functions residing in `main.cpp` and the string reversal function in `oneiteration.cpp`.
- Since the structure of the program changes, and some source files, such as `createarray.cpp`, `errorsarray.cpp` and `printcompact.cpp` are being replaced, it is necessary to modify the `makefile` accordingly to allow Linux users to compile easily.
- What might cause additional running time in the program are function calls which could be replaced by macros. For example, a procedure in `buildk.cpp` compares two characters in the two strings that were provided, and changes the value of an array entry according to the result of the conditional. Since the body of this method is very short, but it requires passing several arguments, it would be more efficient to implement it as a macro.
- Finally, while changes are made to the program regularly, it is important to keep track of any new aspect, so comments, names and dates should be added to new lines of code.
- Following is the section from `tasks.txt` pertaining to today:

@@@@@

Wednesday, 02.19.2020 (Meeting w/ Professor Sokol)
Short Term:
- Change `main.cpp` such that is reads the entire sequence at once and stores it in a large buffer. This will eliminate the need to use `strcat` functions. (Done)
- Replace the `strrev` function by pointer operation: to do so, make pointers decrement to search throughout the string. (Done)
- Keep track of the sizes of the processed string segments. (Done)
- Replace short and frequently used functions or methods by macros. (Done)
- Update the `makefile` according to the setup of the program. (Done)
- Add comments and dates to new lines of code. (Done)

@@@@@

### Monday, February 10, 2020
- 2nd meeting with Professor Sokol in the Spring 2020 Semester.
- The focus was driven towards assessing the correctness of Professor Mermelstein's version but keeping the iterative (rather than recursive) approach when calling the `OneIteration` function.
- Following is the section from `tasks.txt` pertaining to today:

@@@@@

Monday, 02.10.2020 (Meeting w/ Professor Sokol)
Short Term:
- Compile Professor Mermelstein's version of the program as well.
- Make the program read input from FASTA-formatted files. (Done)
- Fetch the most recent Homo Sapiens DNA sequences as FASTA files (Done).
  - The new files are in the directory miriam-new/fasta-new/ .
- Adopt the iterative version of the Main-Lorentz approach, as in Professor Tojeira's code. (Done)

@@@@@

### Wednesday, February 5, 2020
- First meeting with Professor Sokol in the Spring 2020 Semester.
- Code efficiency issues at the technical level and their potential solutions were discussed.
- The work for the following few weeks was outlined (see Tasks below.)
- A Linux account pertaining to the TRedD website was created.
- Professor Sokol included the work of Professor Mermelstein (`tredAri` folder) in the account.
- The folder `miriam-new` was created. It will contain the program files created during this research.
- `tasks.txt`, which will include short and long term tasks as part of the research, was created.

### Tasks
The textfile `tasks.txt` (`/home/mbriskman/miriam-new/tasks.txt`) was added the following targets for today (02.05.2020):

@@@@@

Wednesday, 02.05.2020 (Meeting w/ Professor Sokol)

Short Term:
- Compile and run the TRed (Version 3) program in its current form to check if no errors found. (Done)
- Understand the program's flow, and what role each of its files plays. (Done)
- Translate the TRed (Version 3) program into C++. (Done)
- Pay attention to replacement of:
  - unrecommended techniques, such as goto statements (Done)
  - strlen, strcat  and other O(n) functions (Done)
  - difficultly readable lines, such as ones without spacing and parenthesization between intricate expressions (Done)
  - short functions that could be specified as inline in C++ for optimization (Done)
- Create a short program reading in, reversing, and finding the lenght of a chromosome's nucleotide sequence, and measure how much time each of these 3 tasks took to run (in milliseconds): (Done)
  - Use functions from <ctime> to measure the time.
  - Specify the features of the computer in which the program was run (RAM, # of cores, etc.)
    - Here is the information from the Settings -> System -> About page of my Windows 10 computer:
    
          Processor:      Intel(R) Core(TM) i5-2300 CPU @ 2.80GHz  2.80GHz
          Installed RAM:  6.00 GB
          System Type:    64-bit operating system, x64-based processor
    - This will be at least one of the computers to run this program.

Long Term:
- Replace string reversal and concatenation methods with pointer operations. (Done)
- Implant the suffix arrays implementation into the program.
- Attempt to detect other remediable aspects of the program and address them.
- Test the correctness of the program with the made changes. 
- Run the program on the human chromosomes and note if an improvement in the running time was made.
  - To be capable of comparing the total running time, we must run the TRed Version 3 (current) on the same device(s). 

@@@@@

In addition, the `tasks.txt` file mentions the C++ documentation website

http://www.cplusplus.com/reference/

which can provide a great deal of useful information about C++ procedures, including their complexity and thrown exceptions.

<hr>

### Monday, January 6, 2020
- Registration permission processed, and registration for class completed.
- Spring 2020 Semester begins on January 27.

<hr>

### Friday, January 3, 2020
- Draft approved by Professor Sokol and sent to the Chairperson for approval and class registration permission.
- Proposal approved by the Chairperson.

<hr>

### Thursday, January 2, 2020
- 2nd draft of the project proposal sent for verification.

<hr>

### Wednesday, January 1, 2020
- 1st draft of the project proposal sent to Professor Sokol for verification.

<hr>

### Thursday, December 9, 2019
- What the topic of the study should be is discussed with Professor Sokol.
- The agreed topic theme: **"Efficiently Implementing + Applying the Landau-Vishkin Algorithm to TRed"**.

<hr>

**You have reached the bottom of the file; no older records exist!**

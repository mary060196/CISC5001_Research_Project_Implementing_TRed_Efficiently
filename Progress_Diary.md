# Study's Progress Diary
This file keeps track of how the research project advances. Mentions of milestones and important notes are introduced.

**Author**: Miriam Briskman

<hr>

### Friday, May 22, 2020

- End of Spring 2020 semester.

<hr>

### Wednesday, May 20, 2020

- Releases `v4.0.1` and `v4.1.1-alpha` were created. See:

[https://github.com/mary060196/CISC5001_Research_Project_Implementing_TRed_Efficiently/releases/tag/v4.0.1](https://github.com/mary060196/CISC5001_Research_Project_Implementing_TRed_Efficiently/releases/tag/v4.0.1)

and

[https://github.com/mary060196/CISC5001_Research_Project_Implementing_TRed_Efficiently/releases/tag/v4.1.1-alpha](https://github.com/mary060196/CISC5001_Research_Project_Implementing_TRed_Efficiently/releases/tag/v4.1.1-alpha)

On this pages, one can download a zipped folder (either `zip` or `tar.gz`) containing the source code files and a few text files with an example. 

- The releases are based on the `multi-threaded` and `suffix-array` branches created early in the repository.

- **Changes**:

   - Reverting all macros in `parameters.h` to their default values (see README.md for a listing of the macros.)
   - Adding a note in README.md stating that the example output files provided in the directory were created after changing some of the default macros in `parameters.h` (see README.md for exact values.)
   
<hr>

### Monday, May 18, 2020

- Releases `v4.0` and `v4.1-alpha` were created. See:

[https://github.com/mary060196/CISC5001_Research_Project_Implementing_TRed_Efficiently/releases/tag/v4.0](https://github.com/mary060196/CISC5001_Research_Project_Implementing_TRed_Efficiently/releases/tag/v4.0)

and

[https://github.com/mary060196/CISC5001_Research_Project_Implementing_TRed_Efficiently/releases/tag/v4.1-alpha](https://github.com/mary060196/CISC5001_Research_Project_Implementing_TRed_Efficiently/releases/tag/v4.1-alpha)

On this pages, one can download a zipped folder (either `zip` or `tar.gz`) containing the source code files and a few text files with an example. 

- The releases are based on the `multi-threaded` and `suffix-array` branches created early in the repository.

- The `translation21` directory was updated, and `translation22` directory was created. They contain the programs corresponding to `v4.0` and `v4.1-alpha`, respectively.

- The `README` files in each directory were updated to reflect the changes.

The `README` for version `v4.0` can be found at:

[https://github.com/mary060196/CISC5001_Research_Project_Implementing_TRed_Efficiently/blob/multi-threaded/README.md](https://github.com/mary060196/CISC5001_Research_Project_Implementing_TRed_Efficiently/blob/multi-threaded/README.md)

and this for version `v4.1-alpha` can be found at:

[https://github.com/mary060196/CISC5001_Research_Project_Implementing_TRed_Efficiently/blob/suffix-array/README.md](https://github.com/mary060196/CISC5001_Research_Project_Implementing_TRed_Efficiently/blob/suffix-array/README.md)

- The TRed website that was built for the repository at

[https://mary060196.github.io/CISC5001_Research_Project_Implementing_TRed_Efficiently/](https://mary060196.github.io/CISC5001_Research_Project_Implementing_TRed_Efficiently/)

was updated to contain the links to the releases described above. The website is active and can be accessed from any internet-connected device.

- The Final Report for the CISC 5001 class is ready, and can be found at

[https://github.com/mary060196/CISC5001_Research_Project_Implementing_TRed_Efficiently/blob/master/CISC-5001-Spring-2020-Efficiently-Implementing-TRed-Final-Report-Miriam-Briskman.pdf](https://github.com/mary060196/CISC5001_Research_Project_Implementing_TRed_Efficiently/blob/master/CISC-5001-Spring-2020-Efficiently-Implementing-TRed-Final-Report-Miriam-Briskman.pdf)

<hr>

### Sunday, May 17, 2020

- The `filter` and `nofilter` programs were updated to work with FASTA formatted input files.

<hr>

### Sunday, May 10, 2020

- It was discovered that by changing the initialization of the `L` matrix, it is possible to generate elements of the `KxK` matrix. The initialization changes consist of creating arrow-like boundaries to the `L` matrix. The rest of the querying job is the same as in the original implementation of the algorithm by Brani Cohen at

[https://github.com/Brani/Implementation-of-the-Landau-Vishkin-Algorithm](https://github.com/Brani/Implementation-of-the-Landau-Vishkin-Algorithmhttps://github.com/Brani/Implementation-of-the-Landau-Vishkin-Algorithm) 

- The `NewMatrixL.cpp` and `NewMatrixL.h` source code files, which contain the updated Landau-Vishkin suffix arrays algorithm, were created and tested against the current edit distance algorithm.

- A new `translation22` will be created soon to reflect this progress and allow global access to the TRed program operating with the Landau-Vishkin suffix arrays algorithm.

- A GitHub Pages website for this repository was created at

[https://mary060196.github.io/CISC5001_Research_Project_Implementing_TRed_Efficiently/](https://mary060196.github.io/CISC5001_Research_Project_Implementing_TRed_Efficiently/)

The website will contain links to the programs that will be ready soon. It is designed to show the date when the website was recently modified and a unique visitors hit counter.

<hr>

### Wednesday, April 29, 2020

- 6th meeting with Professor Sokol in the Spring 2020 Semester.
- The features of the program in `translation21` were discussed.
- The idea behind changing the suffixes array algorithm was discussed.
- The structure of the final report that is going to be written for this CISC 5001 class was discussed.

<hr>

### Tuesday, April 28, 2020

- We began working on merging the suffixes array (Landau-Vishkin '98) algorithm into the program. 
- This requires changing the programming techniques in the current C++ version of the algorithm. For example, we replace all C++ strings by C strings (we proved earlier this semester that the code runs faster this way.) To speed up the algorithm further, we extract the code from the various class definitions and place it into routines (this way, less function calls, which consume time, are made.)
- More algorithmic simplifications can be made to the program, such as the usage of suffix properties, to replace `O(n)` or `O(lg n)` operations by `O(1)` ones.
- Three goals are pursued while changing the suffixes array algorithm:
  1. Making the code fit into the current programming approach in TRed (as explained above.) (Done)
  2. Finding out how to extract the `KxK` matrix out of the `L` matrix, which is what the algorithm currently builds. (Done)
  3. Having the algorithm run faster than the currently used algorithm for computing the `KxK` matrix.
     - The current suffix-arrays algorithm runs slower than the original `KxK` matrix algorithm.
  
<hr>

### Friday, April 24, 2020

- It was recalled that, before this Spring 2020 semester began, Professor Sokol indicated that, in one of the previous projects, Professor Ari Mermelstein of Brooklyn College wanted to add a multithreading feature to the program to make it utilize all the existing CPU cores to decrease the total running time of the program.
- Since almost every modern computer is multicore, there is a very high chance that the time it take the program to run will decrease by at least 50% (for a minimum of 2 cores). In general, when disregarding the synchronization overhead, the running time should decrease by `N` times, where `N` is the number of present cores.
- Having this idea in mind, a new version of the program, `translation21`, which creates `N` threads, where `N` is the number of present cores, was created, and can be found in this repository at [https://github.com/mary060196/CISC5001_Research_Project_Implementing_TRed_Efficiently/tree/master/translation21](https://github.com/mary060196/CISC5001_Research_Project_Implementing_TRed_Efficiently/tree/master/translation21).
- The premise that allows using multithreading in this program is that the tandem repeats in one section with `DOUBLE_PERIOD` characters do not depend on the tandem repeats found in another, non-overlapping section of that size. Therefore, we can divide all such sections among the threads created in the program to do an individual work that will not interfere with the results of other threads.
- When testing the program with different number of threads activated, it was noticed that *the program runs the fastest when the number of threads is equal to the number of cores on the machine*.
- Therefore, the program in `translation21` retrieves the number of cores in the computer and activates threads of that quantity.
- The overhead in using threading in this program constitutes of:
  - Creating and joining the threads.
  - Allocation and deallocation of heap memory for each thread. Each of the threads uses its own memory storage for computations to prevent race conditions. This includes the memory allocated for various arrays and matrices used in the program. Mutexing is not used in the program, as doing so in the context of *this* program would cancel the benefits of having each thread working on its own without blocking and context switching.
  - Appending all the output files into a single output file. Each of the threads writes the detected tandem repeats locations within the sequence's region on which it works into its own output file created in the same directory as the program's binary to prevent mixing the output. We want to have a sequential list of the detected repeats, so it is important to keep them in the order they were found. After all threads are joined, the main thread copies each of the output files into one output file (whose name the user entered as a command line argument) and destroys the temporary output files.
- The computer with specifications

          Processor:      Intel(R) Core(TM) i5-2300 CPU @ 2.80GHz  2.80GHz
          Installed RAM:  6.00 GB
          System Type:    64-bit operating system, x64-based processor
          
  has 4 cores, so the gaining in the running time of `translation21` compared to `translation20` is expected to be of **75%**.
- When running the `translation21` program on Chromosome Y, it took `1656296` milliseconds for the program to run, which is about `0.46` hours (slightly less than 1/2 hour.) Compared to the `translation20` program, which ran `6436886` milliseconds, the `translation21` program ran **74.2687%** faster, which is quite close to the gain we hypothesized.
- Compared to the original TRed version 3, which ran `8420086` milliseconds, the `translation21` program ran **80.3292%** faster.

<hr>

### Monday, April 20, 2020

- Both the original TRed Version 3 program and the program in `translation20` were run on Chromosome 2, which is currently the chromosome with the longest provided sequence (242,694 KB).
- The TRed Version 3 program ran for `66246724` milliseconds, with is about 18.4 hours.
- The `translation20` program ran for `57406157` milliseconds, with is about 15.95 hours.
- Therefore, `translation20` ran **13.35%** faster than TRed Version 3 did.
- This rejects the hypothesis that program in `translation20` will run about **25%** faster than TRed Version 3 would.
- It is possible, though, that the percentage gaining in time varies from one chromosome sequence to another, even if their length is similar, because their sequences are different.
- Thus, we might notice a different result for current Chromosome 1 sequence, whose length is a bit smaller than this of Chromosome 2.

<hr>

### Sunday, April 19, 2020

- An attempt to use the new algorithm for reducing the computations of the `NxN` matrices was made, but the program ran 10 times longer than before, not resulting in any speedup of the program.
- The reason for that is that the original code in `buildk.cpp` computes only about `2*K` rows of the `NxN` matrix, where `K` is the number of allowed errors.
- The region in the `NxN` matrix that potentially allows the reduction in work is usually below those rows, which prevents us from introducting the new algorithm into the program.
- However, the newly added `translation20` into the repository completely prevents the dependency on the `-1`s regions in the `KxK` matrices, letting us stop initializing the elements in those regions to `-1` every time the `buildmatrixforward` and `buildmatrixbackward` functions are called.
- Please refer to the directory at [translation20](https://github.com/mary060196/CISC5001_Research_Project_Implementing_TRed_Efficiently/tree/master/translation20) to view the new version of the program. The `Note.md` file describes the new features in `translation20` and differences between the programs in `translation20` and `translation19`.
- When running the program in `translation20` on Chromosome Y, it took `6436886` milliseconds to run, which is about `1.79` hours. This is about **23.5532%** less running time compared to this of the original TRed Version 3 program.
- The expectation is that the program will run about **25%** time less than the original TRed Version 3 program does on the large chromosomes.
- From now until the end of the semester, the only focused-upon topic is going to be the integration of the suffix arrays algorithm based on the Landau-Vishkin '98 suffix trees algorithm into the program and measuring how this addition impacts the running time of the program.
- Short-term tasks:

  - Setting up the suffix arrays algorithm to fit the context (char arrays rather than strings, structs rather than classes, etc) (Done)
  - Recognition of the `KxK` matrix inside the `L` matrix constructed in the algorithm (Done)
  - Implantation of the algorithm into the program (Done)
  - Running the program on Chromosome Y and other chromosomes to check how fast it runs
    - The new program was ran on a portion of the Y chromosome, and showed to run slower than the lastest, non-suffix-arrays version of the program.
  
<hr>

### Sunday, April 5, 2020

- During the last week, a way to reduce the work done while computing the Edit Distance matrix was discovered.
- Corresponding algorithms are added to the repository in the `Edit_Distance_-_New_Algorithms` directory.
- The work is reduced by constant time. Since on average, the completed work is 4.5xAB out of the 7xAB, where A is the `text` sequence and B is the `pattern` sequence, the savings are (7-4.5)/7 =~ 36% of the job.
- The algorithms were tested. The resulting matrices were matched against matrices created via the traditional Edit Distance algorithm.
- Please view the individual directories in the `Edit_Distance_-_New_Algorithms` directory for further information and instructions on how to run the algorithms.
- Short-term tasks:

  - Devising an algorithm to construct a `KxK` matrix from an already given Edit Distance matrix. (Done)
  - Testing the validity of the algorithm against the currently used algorithm in `buildk.cpp` file of the TRed version 3 program. (Done)
  - Integration of the algorithms into a new "translated" code, which will be called `translation20`, as a continuation to the `translation19` program, whose code is located in the `translation19` directory in this repository. (Done)
    - The program in `translation22` and the one in the `suffix-array` branch exhibit this implementation.
  
<hr>

### Monday, March 16, 2020

- 5th meeting with Professor Sokol in the Spring 2020 Semester.
- Transition from the `L` matrix to the `KxK` matrix discussed.
- The algorithm of the suffix arrays initializes the `L` matrix differently than the `KxK` algorithm does.
- Changing the initialization could result in obtaining a portion of a `KxK` matrix within the corresponding `L` matrix for the same text and pattern.
- The slides at [https://u.cs.biu.ac.il/~amir/PMslides/IndexingProblem.pdf](https://u.cs.biu.ac.il/~amir/PMslides/IndexingProblem.pdf) discussed and explained.
- Short-term tasks (for the remaining weeks):

  - Integration of the suffix arrays algorithm into the TRed version 3 program by changing the way the initialization takes place. (Done)
  - Running the new program and finding out how much time, compared to the TRed version 3 program, it runs on the same sequence. (Done)
  - Adoption of the *incremental string comparison* approach? This task might not take place since the article where the approach is described is difficult to understand in a practical way.
    - From one `for` iteration inside `oneiteration.cpp`, only `O(n)` work is needed to update the `lcp` (longest common prefix) array. Is this considered an implementation of the *incremental string comparison* approach?
  
<hr>

### Sunday, March 15, 2020
- The original TRed Version 3 and the translated code in [translation19](https://github.com/mary060196/CISC5001_Research_Project_Implementing_TRed_Efficiently/tree/master/translation19) were run, separately, on the entire Chromosome Y fasta file.
- The fasta file, together with all the other Chromosome sequence files, was earlier downloaded from

      [ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/ARCHIVE/ANNOTATION_RELEASE.109/](ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/ARCHIVE/ANNOTATION_RELEASE.109/)

- The original program ran for `8420086` miliseconds (= 2.34 hours,) while the translated code ran for `6800100` miliseconds (= 1.89 hours.)
- This is a decrease of **19.24%** in the total running time.
- The programs were run on a desktop computer having the following specifications:

          Processor:      Intel(R) Core(TM) i5-2300 CPU @ 2.80GHz  2.80GHz
          Installed RAM:  6.00 GB
          System Type:    64-bit operating system, x64-based processor
          
- It is expected that the decrease when running the programs on larger chromosome sequences (e.g., Chromosome 1 or 2,) will be greater than **20%**.

<hr>

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

<hr>

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
- After all the key short-term goals (including those above) are completed, begin reading about suffix arrays and thinking how to integrate them into the program, hoping to achieve greater efficiency. (Done)

@@@@@

<hr>

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

<hr>

### Monday, February 10, 2020
- 2nd meeting with Professor Sokol in the Spring 2020 Semester.
- The focus was driven towards assessing the correctness of Professor Mermelstein's version but keeping the iterative (rather than recursive) approach when calling the `OneIteration` function.
- Following is the section from `tasks.txt` pertaining to today:

@@@@@

Monday, 02.10.2020 (Meeting w/ Professor Sokol)
Short Term:
- Compile Professor Mermelstein's version of the program as well.
  - Not completed, as Mermelstein's version is recursive -- a possible source for stack overflow.
- Make the program read input from FASTA-formatted files. (Done)
- Fetch the most recent Homo Sapiens DNA sequences as FASTA files (Done).
  - The new files are in the directory `miriam-new/fasta-new/`.
- Adopt the iterative version of the Main-Lorentz approach, as in Professor Tojeira's code. (Done)

@@@@@

<hr>

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
- Compile and run the TRed (Version 3) program in its current form to check if no errors are found. (Done)
- Understand the program's flow, and what role each of its files plays. (Done)
- Translate the TRed (Version 3) program into C++. (Done)
- Pay attention to replacement of:
  - unrecommended techniques, such as goto statements (Done)
  - strlen, strcat  and other O(n) functions (Done)
  - difficultly readable lines, such as ones without spacing and parenthesization between intricate expressions (Done)
  - short functions that could be specified as inline in C++ for optimization (Done)
- Create a short program reading in, reversing, and finding the length of a chromosome's nucleotide sequence, and measure how much time each of these 3 tasks takes to run (in milliseconds): (Done)
  - Use functions from `<ctime>` to measure the time. (Done)
    - We use more precise time measuring options from the `chrono` C++ library.
  - Specify the features of the computer in which the program was run (RAM, # of cores, etc.) (Done)
    - Here is the information from the Settings -> System -> About page of my Windows 10 computer:
    
          Processor:      Intel(R) Core(TM) i5-2300 CPU @ 2.80GHz  2.80GHz
          Installed RAM:  6.00 GB
          System Type:    64-bit operating system, x64-based processor
    - This will be at least one of the computers to run this program.

Long Term:
- Replace string reversal and concatenation methods with pointer operations. (Done)
- Implant the suffix arrays implementation into the program. (Done)
- Attempt to detect other remediable aspects of the program and address them. (Done)
- Test the correctness of the program with the made changes. (Done)
- Run the program on the human chromosomes and note if an improvement in the running time was made. (Done)
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

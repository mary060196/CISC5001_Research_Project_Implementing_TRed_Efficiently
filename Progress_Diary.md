# Study's Progress Diary
This file tracks how the research project advances. Mentions of milestones and important notes are introduced.

### Monday, February 10, 2020
- 2nd meeting with Professor Sokol in the Spring 2020 Semester.
- The focus was driven towards assessing the correctness of Professor Mermelstein's version but keeping the iterative (rather than recursive) approach when calling the `OneIteration` function.
- Following is the section from `tasks.txt` pertaining to today:

@@@@@

Monday, 02.10.2020 (Meeting w/ Professor Sokol)
Short Term:
- Compile Professor Mermelstein's version of the program as well.
- Make the program read input from FASTA-formatted files.
- Fetch the most recent Homo Sapiens DNA sequences as FASTA files (Done).
  - The new files are in the directory miriam-new/fasta-new/ .
- Adopt the iterative version of the Main-Lorentz approach, as in Professor Tojeira's code.

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
- Compile and run the TRed (Version 3) program in its current form to check if no errors found.
- Understand the program's flow, and what role each of its files plays.
- Translate the TRed (Version 3) program into C++.
- Pay attention to replacement of:
  - unrecommended techniques, such as goto statements
  - strlen, strcat  and other O(n) functions
  - difficultly readable lines, such as ones without spacing and parenthesization between intricate expressions
  - short functions that could be specified as inline in C++ for optimization.
- Create a short program reading in, reversing, and finding the lenght of a chromosome's nucleotide sequence, and measure how much time each of these 3 tasks took to run (in milliseconds):
  - Use functions from <ctime> to measure the time.
  - Specify the features of the computer in which the program was run (RAM, # of cores, etc.)
    - Here is the information from the Settings -> System -> About page of my Windows 10 computer:
    
          Processor:      Intel(R) Core(TM) i5-2300 CPU @ 2.80GHz  2.80GHz
          Installed RAM:  6.00 GB
          System Type:    64-bit operating system, x64-based processor
    - This will be at least one of the computers to run this program.

Long Term:
- Replace string reversal and concatenation methods with pointer operations.
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

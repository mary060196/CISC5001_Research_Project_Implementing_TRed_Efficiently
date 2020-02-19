/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\ 
*                Comparing Running Time of C-Strings and C++ String Object                       *
* Author:              Supervising Professor:     For:                                           *
*    Miriam Briskman       Dina Sokol                CISC 5001 Independent Study                 *
*                                                    Spring 2020                                 *
*                                                    Brooklyn College.                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cmath>
#include <string>
#include <chrono>
#include <time.h>
#include <iostream>
#include <fstream>
#include <algorithm>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Function Prototypes:
void user (unsigned long long&);
void create_seq (unsigned long long);
void process_seq (unsigned long long);

using namespace std;

const char* filename  = "newSeq.fsa"; // The sequence file

int main (int argc, char *argv[])
{
    // Variable Declarations:
    unsigned long long n;   // To-be-retrieved length of the sequence 
	
    // Instantiating n or user - prompting:
    if (argc > 1)
    {
        n = strtoull (argv[1], NULL, 10);
        cout << "\nThank you! The integer you entered is " 
             << n << "." << endl << endl << "Processing..." << endl;
    }
    else
        user (n);

    // Create demo sequence in FASTA format:
    create_seq (n);

    // Process the sequence:
    process_seq (n);

    return 0;
}

void user (unsigned long long &n)
{
    cout << "Please enter a positive integer corresponding to the length of the sequence: >> ";
    cin >> n;
    cout << "\nThank you! The integer you entered is " 
         << n << "." << endl << endl << "Processing..." << endl;
    return;   
}

void create_seq (unsigned long long n)
{
    FILE *outfile = fopen (filename, "w"); // Open a new file to write into.
    srand((long)time(NULL));               // Seed rand().

    unsigned long long i, c;
    short temp;

    // Capture the beginning time of the creation of the sequence:
    auto start_time = chrono::high_resolution_clock::now();

    fprintf (outfile, "> Bla Bla FASTA Format\n");
    for (i = 1, c = 1; i <= n; i++, c++)
    {
        temp = floor(rand()/(double)RAND_MAX * 4);
        switch (temp)
        {
           case 0: fprintf (outfile, "A");
                   break;
           case 1: fprintf (outfile, "C");
                   break;
           case 2: fprintf (outfile, "G");
                   break;
           case 3: fprintf (outfile, "T");
                   break;
           case 4: fprintf (outfile, "N\n> Bla Bla FASTA Format\n");
                   c--;
                   break;
        }
        if (c % 80 == 0) 
           fprintf (outfile, "\n");
    }
    // Capture the ending time of the creation:
    auto end_time = std::chrono::high_resolution_clock::now();
    auto mtime = end_time - start_time;

    cout << "It took " << mtime/std::chrono::milliseconds(1) 
         << " milliseconds to create a FASTA file with " << n << " bases." << endl;

    fclose (outfile);
    return;
}

void process_seq (unsigned long long n)
{
    FILE *infile = fopen (filename, "r"); // Open a new file to read from.
    ifstream infile_cpp (filename); // Open another file to read from.
    string temporary;
    char *cstr;
    int c;
    unsigned long long i = 0;

    // Read the sequence in:

    // Capture the beginning time of the string's reading:
    auto start_time = chrono::high_resolution_clock::now();
    
    while ((c = fgetc(infile)) != EOF)
    {
        switch (c)
        {
            case '>':  while ((c = fgetc(infile)) != '\n');
                       break;
            case '\n': break;
            default:   temporary += (char) c;
                       break;
        }
    }

    // Capture the ending time of the reading:
    auto end_time = std::chrono::high_resolution_clock::now();
    auto mtime = end_time - start_time;

    cout << "It took " << mtime/std::chrono::milliseconds(1) 
         << " milliseconds to read the sequence in using 'fgetc' and \"+=\" into a string object." << endl;

    fseek (infile, 0, SEEK_SET); // Place file position to beginning
    
    i = 0;

    // Capture the beginning time of the string's reading:
    start_time = chrono::high_resolution_clock::now();
    
    cstr = (char*) malloc ((n + 10) * sizeof(char));
    
    while ((c = fgetc(infile)) != EOF)
    {
        switch (c)
        {
            case '>':  while ((c = fgetc(infile)) != '\n');
                       break;
            case '\n': break;
            default:   cstr[i++] = (char) c;
                       break;
        }
    }
    cstr[i] = '\0';
    // Capture the ending time of the reading:
    end_time = std::chrono::high_resolution_clock::now();
    mtime = end_time - start_time;

    cout << "It took " << mtime/std::chrono::milliseconds(1) 
         << " milliseconds to read the sequence in using 'fgetc' into a c-string." << endl;

    fclose (infile); // Close the file

    // Capture the beginning time of the string's conversion:
    start_time = chrono::high_resolution_clock::now();
    string str (cstr);
    // Capture the ending time of the reading:
    end_time = std::chrono::high_resolution_clock::now();
    mtime = end_time - start_time;

    cout << "It took " << mtime/std::chrono::milliseconds(1) 
         << " milliseconds to convert a c-string into a string object." << endl;

    i = 0;
    temporary.clear();
    free (cstr);

    // Capture the beginning time of the string's reading:
    start_time = chrono::high_resolution_clock::now();
    
    while ((c = infile_cpp.get()) != EOF)
    {
        switch (c)
        {
            case '>':  while ((c = infile_cpp.get()) != '\n');
                       break;
            case '\n': break;
            default:   temporary += (char) c;
                       break;
        }
    }

    // Capture the ending time of the reading:
    end_time = std::chrono::high_resolution_clock::now();
    mtime = end_time - start_time;

    cout << "It took " << mtime/std::chrono::milliseconds(1) 
         << " milliseconds to read the sequence in using fstream 'get' and \"+=\" into a string object." << endl;

    infile_cpp.close(); // Close the file
    infile_cpp.open (filename); // Reopen the file
    i = 0;

    // Capture the beginning time of the string's reading:
    start_time = chrono::high_resolution_clock::now();
    cstr = (char*) malloc ((n + 10) * sizeof(char));
    
    while ((c = infile_cpp.get()) != EOF)
    {
        switch (c)
        {
            case '>':  while ((c = infile_cpp.get()) != '\n');
                       break;
            case '\n': break;
            default:   cstr[i++] = (char) c;
                       break;
        }
    }
    cstr[i] = '\0';
    // Capture the ending time of the reading:
    end_time = std::chrono::high_resolution_clock::now();
    mtime = end_time - start_time;

    cout << "It took " << mtime/std::chrono::milliseconds(1) 
         << " milliseconds to read the sequence in using fstream 'get' into a c-string." << endl;

    infile_cpp.close(); // Close the file

    // Capture the beginning time of the string's conversion:
    start_time = chrono::high_resolution_clock::now();
    string str2 (cstr);
    // Capture the ending time of the reading:
    end_time = std::chrono::high_resolution_clock::now();
    mtime = end_time - start_time;

    cout << "It took " << mtime/std::chrono::milliseconds(1) 
         << " milliseconds to convert a c-string into a string object." << endl;

    // Reverse it:

    // Capture the beginning time of the string's reversal:
    start_time = chrono::high_resolution_clock::now();
    reverse (str.begin(), str.end());
    // Capture the ending time of the reversal:
    end_time = std::chrono::high_resolution_clock::now();
    mtime = end_time - start_time;

    cout << "It took " << mtime/std::chrono::milliseconds(1) 
         << " milliseconds to reverse the sequence using 'reverse' (<algorithm>)." << endl;

    // Capture the beginning time of the string's reversal:
    start_time = chrono::high_resolution_clock::now();
    strrev (cstr);
    // Capture the ending time of the reversal:
    end_time = std::chrono::high_resolution_clock::now();
    mtime = end_time - start_time;

    cout << "It took " << mtime/std::chrono::milliseconds(1) 
         << " milliseconds to reverse the sequence using 'strrev'." << endl;

    // Concatenate it:

    // Capture the beginning time of the string's concatenation:
    start_time = chrono::high_resolution_clock::now();
    str += "bla-bla";
    // Capture the ending time of the concatenation:
    end_time = std::chrono::high_resolution_clock::now();
    mtime = end_time - start_time;

    cout << "It took " << mtime/std::chrono::milliseconds(1) 
         << " milliseconds to concatenate the sequence with 'bla-bla' using '+='." << endl;

    // Capture the beginning time of the string's concatenation:
    start_time = chrono::high_resolution_clock::now();
    strcat (cstr, "bla-bla");
    // Capture the ending time of the concatenation:
    end_time = std::chrono::high_resolution_clock::now();
    mtime = end_time - start_time;

    cout << "It took " << mtime/std::chrono::milliseconds(1) 
         << " milliseconds to concatenate the sequence with 'bla-bla' using 'strcat'." << endl;

    // Find the string's length:

    // Capture the beginning time of the string's length-measuring:
    start_time = chrono::high_resolution_clock::now();
	cout << "The length of the string is: " << str.length() << " (including 7 chars of 'bla-bla'.)" << endl;
    // Capture the ending time of the measurement:
    end_time = std::chrono::high_resolution_clock::now();
    mtime = end_time - start_time;

    cout << "It took " << mtime/std::chrono::milliseconds(1) 
         << " milliseconds to find the sequence's length using 'length'." << endl;

    // Capture the beginning time of the string's length-measuring:
    start_time = chrono::high_resolution_clock::now();
	cout << "The length of the string is: " << strlen(cstr) << " (including 7 chars of 'bla-bla'.)" << endl;
    // Capture the ending time of the measurement:
    end_time = std::chrono::high_resolution_clock::now();
    mtime = end_time - start_time;

    cout << "It took " << mtime/std::chrono::milliseconds(1) 
         << " milliseconds to find the sequence's length using 'strlen'." << endl;

    // Change content of string:

    // Capture the beginning time of the string's length-measuring:
    start_time = chrono::high_resolution_clock::now();
	for (int j = 0; j < n; j++)
        str[j] = 'A';
    // Capture the ending time of the measurement:
    end_time = std::chrono::high_resolution_clock::now();
    mtime = end_time - start_time;

    cout << "It took " << mtime/std::chrono::milliseconds(1) 
         << " milliseconds to change the content of a string object to all 'A's." << endl;

    // Capture the beginning time of the string's length-measuring:
    start_time = chrono::high_resolution_clock::now();
	for (int j = 0; j < n; j++)
        cstr[j] = 'A';
    // Capture the ending time of the measurement:
    end_time = std::chrono::high_resolution_clock::now();
    mtime = end_time - start_time;

    cout << "It took " << mtime/std::chrono::milliseconds(1) 
         << " milliseconds to change the content of a c-string to all 'A's." << endl;

    free (cstr);
}

//   |=================================================================|
//   |   TRED: a tool for detecting Tandem Repeats within sequences,   |
//   |               using the Edit Distance metric.                   |
//   |=================================================================|

// Copyright © 2007-2009 Dina Sokol, Justin Tojeira
// Distributed under the Aladdin Free Public License

// THIS SOFTWARE SHOULD BE ACCOMPANIED BY readme.txt AND license.html
// WE STRONGLY ENCOURAGE YOU TO READ BOTH BEFORE PROCEEDING
//------------------------------------------------------------------------------

#include <stdio.h>
#include <vector>
#include <iostream>
#include <string>
#include "parameters.h"

using namespace std;

// Structure keeping info on a repeat:
typedef struct{
   long start, end, actualStart, actualEnd, errs, length,
        actualPosition1, actualPosition2;
   vector <long> ssa;
} rep;

// Function Prototypes:
void filter (vector<rep>*, FILE*, FILE*, FILE*, FILE*, FILE*, FILE*, int);
void computeActualPos (rep*, FILE*);
long moves (FILE*);
void convert (rep*, FILE*, FILE*, FILE*, FILE*, FILE*, FILE*);

int main (int argc, char *argv[])
{
   FILE *infile,       // Points to the intermediary file.
        *outfile,      // Points to the alignments file.
        *tabfile,      // Points to the table file.
        *seqfile,      // Points to the CLEANED sequence file.
        *seqfile2,     // Points to the CLEANED sequence file as well.
        *temporary,    // Points to original sequence file.
        *temporary2;   // Points to original sequence file as well.
   string filename[5]; // C++ Strings to hold names of files.
   long temp;          // Temporary
   vector <rep> rv;    // Vector holding info on repeats.

   /*********** FOR INPUT/OUTPUT FILES AS COMMAND LINE ARGUMENTS ************/

   if (argc != 5)
   {
       cout << "Please manually enter the following 4 filenames below:" << endl << endl
            << "1) Enter the name of the sequence file : >> ";
       cin  >> filename[0];
       cout << endl << "2) Enter the name of the intermediary file : >> ";
       cin  >> filename[1];
       cout << endl << "3) Enter the name of the output alignment file : >> ";
       cin  >> filename[2];
       cout << endl << "4) Enter the name of the output table file : >> ";
       cin  >> filename[3];
       cout << endl << endl << "Thank you! Information is being processed..." << endl << endl;

       filename[5] = filename[0] + "temporary.DATA";

       // Open the sequence file and copy relevant content to 'temporary' (no comments, 
       //   no whitespaces):
       if ((seqfile = fopen(filename[0].c_str(), "r")) == NULL){
           cout << "Error opening sequence file.\n";
           exit (EXIT_FAILURE); 
       }
       if ((temporary = fopen(filename[5].c_str(), "w")) == NULL){
           cout << "Error opening temporary file.\n";
           exit (EXIT_FAILURE); 
       }

       char next = 'a'; // A char not equal to EOF to let 'while' execute.

       // Clean the original file from whitespaces and FASTA comments:
       while (next != EOF)
       {
           next = getc(seqfile);
           switch (next)
           {
               case '>':  while (getc(seqfile) != '\n');
               case '\n':
               case '\t':
               case '\r':
               case '\f':
               case '\v':
               case ' ':
               case EOF:  break;
               default:   fputc (next, temporary);
                          break;
           }
       }
       fputc ('#', temporary); // Place '#' at end of cleaned file.
       fflush (temporary);
       fclose (temporary);
       fclose (seqfile);

       // Re-open the file in reading-only mode:
       if ((seqfile = fopen(filename[5].c_str(), "r")) == NULL){
          cout << "Error opening sequence file.\n";
          exit (EXIT_FAILURE);
       }

       if ((seqfile2 = fopen(filename[5].c_str(), "r")) == NULL){ 
           cout << "Error opening sequence file.\n";
           exit (EXIT_FAILURE);
       }
       if ((infile = fopen(filename[1].c_str(), "r")) == NULL){
           cout << "Error opening intermediary file.\n";
           exit (EXIT_FAILURE); 
       }
       if ((outfile = fopen(filename[2].c_str(), "w")) == NULL){ 
           cout << "Error opening alignment file.\n";
           exit (EXIT_FAILURE);
       }
       if ((tabfile = fopen(filename[3].c_str(), "w")) == NULL){ 
           cout << "Error opening output table file.\n";
           exit (EXIT_FAILURE);
       }
       if ((temporary = fopen(filename[0].c_str(), "r")) == NULL){
           cout << "Error opening temporary file.\n";
           exit (EXIT_FAILURE); 
       }
       if ((temporary2 = fopen(filename[0].c_str(), "r")) == NULL){
           cout << "Error opening temporary file.\n";
           exit (EXIT_FAILURE); 
       }
   }
   else
   {
       if ((seqfile = fopen(argv[1], "r")) == NULL)
           { printf("Error opening sequence file.\n"); exit(EXIT_FAILURE); }

       filename[5] = argv[1];
       filename[5] += "temporary.DATA";

       if ((temporary = fopen(filename[5].c_str(), "w")) == NULL){
           cout << "Error opening temporary file.\n";
           exit (EXIT_FAILURE); 
       }

       char next = 'a'; // A char not equal to EOF to let 'while' execute.

       // Clean the original file from whitespaces and FASTA comments:
       while (next != EOF)
       {
           next = getc(seqfile);
           switch (next)
           {
               case '>':  while (getc(seqfile) != '\n');
               case '\n':
               case '\t':
               case '\r':
               case '\f':
               case '\v':
               case ' ':
               case EOF:  break;
               default:   fputc (next, temporary);
                          break;
           }
       }
       fputc ('#', temporary); // Place '#' at end of cleaned file.
       fflush (temporary);
       fclose (temporary);
       fclose (seqfile);

       // Re-open the file in reading-only mode:
       if ((seqfile = fopen(filename[5].c_str(), "r")) == NULL){
          cout << "Error opening sequence file.\n";
          exit (EXIT_FAILURE);
       }

       if ((seqfile2 = fopen(filename[5].c_str(),"r")) == NULL)
           { printf("Error opening sequence file.\n"); exit(EXIT_FAILURE); }
       if ((infile = fopen(argv[2],"r")) == NULL)
           { printf("Error opening intermediary file.\n"); exit(EXIT_FAILURE); }
       if ((outfile = fopen(argv[3],"w")) == NULL)
           { printf("Error opening alignment file.\n"); exit(EXIT_FAILURE); }
       if ((tabfile = fopen(argv[4],"w")) == NULL)
           { printf("Error opening output table file.\n"); exit(EXIT_FAILURE); }
       if ((temporary = fopen(argv[1], "r")) == NULL)
           { printf("Error opening temporary file.\n"); exit(EXIT_FAILURE); }
       if ((temporary2 = fopen(argv[1], "r")) == NULL)
           { printf("Error opening temporary file.\n"); exit(EXIT_FAILURE); }
       printf ("\nInformation is being processed...\n\n");
   }

   /************************************************************************/

   // Print header of table file:
   fprintf (tabfile,"StartIndx  EndIndex Length Period   Reps Errs Score\n");

   while (fscanf(infile, "%ld", &temp) != EOF)
   {
      if (temp == -1) 
      {
          fscanf(infile, "%ld", &temp);
          // Call the 'filter' function, which calls 'convert', which prints the repeats:
          filter(&rv, seqfile, seqfile2, outfile, tabfile, temporary, temporary2, temp);
      }
      else 
      {
          rv.push_back(rep());
          rv.back().start = temp;
          fscanf (infile, "%ld %ld", &(rv.back().end), &(rv.back().errs));
          // Populate vector with position values:
          do {
               fscanf (infile, "%ld", &temp);
               rv.back().ssa.push_back(temp);
               fscanf (infile, "%ld", &temp);
               rv.back().ssa.push_back(temp);
          } while (rv.back().ssa.back() != 0);
      }
   }

   // Closing files and removing temporaries:
   fclose (seqfile);
   fclose (seqfile2);
   fclose (infile);
   fclose (outfile);
   fclose (tabfile);
   fclose (temporary);
   fclose (temporary2);
   remove (filename[5].c_str());
   printf ("Program finished executing!\n");
   return 0;
}

// 'filter'
// Input: rv, r2 (rep vector pointer), seqfile, seqfile2, outfile, tabfile, 
//        temporary1, temporary2 (file pointers), limit (integer).
// Procedure:
// 1) Filter repeats, if requested,
// 2) Send repeat info for printing.
// Author: Justin, Tojeira, 04.29.2009
void filter (vector<rep> *rv, FILE *seqfile, FILE *seqfile2, FILE *outfile, 
             FILE *tabfile, FILE *temporary, FILE *temporary2, int limit)
{
   long pos1, pos2, i, j;

   if (limit == -1)
      for (i = 0; i < rv->size();)
      {
         rep newrepeat = *(rv->begin() + i);
         convert (&newrepeat, seqfile, seqfile2, outfile, tabfile, temporary, temporary2);
         rv->erase (rv->begin() + i);
      }
   else 
   	   for(i = 0; i < rv->size();)
   		   if (rv->at(i).end < limit)
           {
               convert (&(*(rv->begin() + i)), seqfile, seqfile2, outfile, tabfile, temporary, temporary2);
               rv->erase (rv->begin() + i);
      	   }
      	   else 
               i++;
   return; 
}

// 'computeActualPos'
// Input: r1 (rep struct pointer), temporary (file pointer).
// Output: none.
// Procedure:
// 1) Compute the length of repeat
// 2) Compute the actualStart, actualEnd, actualPosition1,
//    and actualPosition2, which are fields of r1.
// 3) The actual position is the position of a particular
//    character in the original (crude) sequence file based
//    on the cleaned sequence file. It is necessary when the
//    original sequence file has FASTA comments or whitespaces.
// Author: Miriam Briskman, 05.17.2020
void computeActualPos (rep *r1, FILE *temporary)
{
    long temp;
    char next = 'a';

    // Compute the length of the repeat and place in 'length':
    r1 -> length = r1->end - r1->start + 1;
    // Compute the starting position of the repeat (which is
    //   also the first position index in the ssa vector):
    // 1) Rewinding the sequence file to the beginning:
    rewind (temporary);
    temp = 1;
    
    r1->actualStart = 1;
    // 2) Finding the actual start position (omitting FASTA comments
    //    and whitespaces.)
    while (temp != r1->start)
    {
       next = getc(temporary);
       switch (next)
       {
            case '>':  do
                         r1->actualStart++;
                       while (getc(temporary) != '\n');
            case '\n':
            case '\t':
            case '\r':
            case '\f':
            case '\v':
            case ' ':
            case EOF: break;
            default: temp++;
                     break;
       }
       r1->actualStart++;
    }

    //printf ("Finished while\n");
    // Repeat the same steps as above for the end position:
    //   Recall that temp == r1->start
    r1->actualEnd = r1->actualStart;
    while (temp != r1->end)
    {
       next = getc(temporary);
       switch (next)
       {
            case '>':  do
                          r1->actualEnd++;
                       while (getc(temporary) != '\n');
            case '\n':
            case '\t':
            case '\r':
            case '\f':
            case '\v':
            case ' ':
            case EOF: break;
            default: temp++;
                     break;
       }
       r1->actualEnd++;
    }
    
    // Compute the actual position in the sequence textfile of
    //   the 1st and 2nd in the ssa vector:
    // 1) The first vector element is always the start position, so
    r1->actualPosition1 = r1->actualStart;
    // 2) Set temp, file position and r1->actualPosition2 to correct initial
    //    values.
    temp = r1->start;
    fseek (temporary, (r1->actualStart - r1->actualEnd)*sizeof(char), SEEK_CUR);
    r1->actualPosition2 = r1->actualPosition1;
    // 3) Find the position of the second element:
    while (temp != r1->ssa.at(1))
    {
        next = getc(temporary);
        switch (next)
        {
           case '>':  do
                         r1->actualPosition2++;
                      while (getc(temporary) != '\n');
           case '\n':
           case '\t':
           case '\r':
           case '\f':
           case '\v':
           case ' ':
           case EOF:  break;
           default: temp++;
                    break;
        }
        r1->actualPosition2++;
    }
    return;
}


// 'moves'
// Input: temporary (file pointer).
// Output: movements (long int).
// Procedure:
// 1) Let 'movements' be 1.
// 2) As long as the file is currently pointing
//    at a whitespace or FASTA comment, increment
//    'movements'.
// 3) If the file is pointing at a relevant character
//    (say, a base), return 'movements'.
// Author: Miriam Briskman, 05.17.2020
long moves (FILE *temporary)
{
   char next = 'a'; // A char not equal to EOF to let 'while' execute.
   int movements = 1;

   while (next != EOF)
   {
       next = getc(temporary);
       switch (next)
       {
           case '>':  do
                        movements++;
                      while (getc(temporary) != '\n');
           case '\n':
           case '\t':
           case '\r':
           case '\f':
           case '\v':
           case ' ':
           case EOF: break;
           default: return movements;
       }
       movements++;
   }
   return movements;
}

// 'convert'
// Input: r1 (rep struct pointer), seqfile, seqfile2, outfile, tabfile, 
//        temporary1, temporary2 (file pointers).
// Output: none.
// Procedure:
//    Print the information about a repeat into the 'alignment' and
//    table file, preserving the correct positions according to the
//    crude sequence.
// Author: Justin, Tojeira, 04.29.2009
//         Miriam Briskman, 05.17.2020
void convert (rep *r1, FILE *seqfile, FILE *seqfile2, 
              FILE *outfile, FILE *tabfile, FILE *temporary1, FILE *temporary2)
{
   // declare and initialize
   long pos1, pos1a, pos2, pos2a, next1, next1a, next2, next2a,
        actualPos1 = 0, actualPos2 = 0,
        nextper, perstart, lastper,
        errors = 0,
        j, i = 4, i2 = 0,
        score,
        ts1, ts2, num;
   int report1 = 2;
   char bp;
   bool done = false, last = false;
   float reps = 1, period;

   computeActualPos (r1, temporary1);

   // Printing the header of the repeat:
   fprintf (outfile, "Start: %ld End: %ld Length: %ld\n\n", 
            r1->actualStart, r1->actualEnd, r1->length);

   pos1 = r1->ssa.at(0);
   pos2 = r1->ssa.at(1);
   actualPos1 = r1->actualPosition1; // Miriam Briskman, 05.17.2020
   actualPos2 = r1->actualPosition2; // Miriam Briskman, 05.17.2020
   
   fseek (seqfile, (pos1 - START_POS)*sizeof(char), 0);
   fseek (seqfile2, (pos2 - START_POS)*sizeof(char), 0);
   fseek (temporary1, (actualPos1 - START_POS)*sizeof(char), 0); // Miriam Briskman, 05.17.2020
   fseek (temporary2, (actualPos2 - START_POS)*sizeof(char), 0); // Miriam Briskman, 05.17.2020

   next1 = r1->ssa.at(2);
   next2 = r1->ssa.at(3);

   if (r1->ssa.size() == 4)
   {
        last = true;
        next2 = next1 + pos2 - pos1;
   }
   lastper = pos2 - pos1;

   // print repeat
   while (!done) // outer loop
   {   
       if (report1)
       {
      	   if (report1 == 2)
         	   fprintf (outfile, "%ld  ", actualPos1);
           else 
           {
         	   for (j = 1; j <= actualPos1; j *= 10)
         		   fprintf (outfile, " ");
         	   fprintf (outfile, "  ");
           }
           pos1a = pos1;
           pos2a = pos2;
           i2 = i - 4;
           while (r1->ssa.at(i2) > pos1) 
              i2 -= 2;
           i2 += 2;

           next1a = r1->ssa.at(i2++);
           next2a = r1->ssa.at(i2++);
           if (i2 >= r1->ssa.size())
              next2a = r1->end + 1;

           while (pos1a < pos2 && pos2a <= r1->end)
           {
         	   if (pos1a < next1a && pos2a < next2a)
               {
         		  fprintf (outfile, "%c", fgetc(seqfile));
                  pos1a++;
                  pos2a++;
               }
               else if (pos1a < next1a)
               {
            	  fprintf (outfile, "%c", fgetc(seqfile));
                  pos1a++;
               }
               else if (pos2a < next2a)
               {
            	  fprintf (outfile, "-");
                  pos2a++;
               }
               else 
               {
                  next1a = r1->ssa.at(i2++);
                  next2a = r1->ssa.at(i2++);
                  if (i2 >= r1->ssa.size())
                     next2a = r1->end + 1;
               }
           }
           fseek (seqfile, (pos1 - START_POS)*sizeof(char), 0);
           fseek (temporary1, (actualPos1 - START_POS)*sizeof(char), 0); // Miriam Briskman, 05.17.2020
           if (report1 == 2)
         	  fprintf (outfile,"  %ld", actualPos2 - 1);
           fprintf (outfile, "\n");
           report1 = 0;
       } // end if(report1)

       perstart = pos1;
       nextper = pos2;
       fprintf (outfile, "%ld  ", actualPos2);

   	   while (pos1 < nextper) // inner loop
       {     
      	   if (pos1 < next1 && pos2 < next2)
           {
               num = moves (temporary2); // Miriam Briskman, 05.17.2020
               bp = fgetc (seqfile2);
               if (bp != '#') // Miriam Briskman, 05.17.2020
                  fprintf (outfile, "%c", bp);
               if (bp != fgetc (seqfile) && bp != '#') 
                  errors++;
               else if (bp == '#') // Miriam Briskman, 05.17.2020
                  actualPos2 -= num;
               pos1++;
               pos2++;
               actualPos1 += moves (temporary1); // Miriam Briskman, 05.17.2020
               actualPos2 += num; // Miriam Briskman, 05.17.2020

         	   if (!report1 && pos2a <= r1->end)
               {
            	   if (pos1a == next1a && pos2a == next2a)
                   {
                       next1a = r1->ssa.at(i2++);
               	       next2a = r1->ssa.at(i2++);
               	       if (i2 >= r1->ssa.size()) 
                          next2a = r1->end + 1;
                       pos1a++; 
                       pos2a++;
                   }
            	   else if (pos1a == next1a) 
                       report1 = 1;
                   else if (pos2a == next2a) 
                       pos1a++;
                   else 
                   {
                       pos1a++; 
                       pos2a++;
                   }
               }
           }
           else 
           {
              if (last) 
              {
                  done = true;
                  break;
              }
         	  while (pos2 < next2) 
              {
                  fprintf (outfile, "%c", fgetc (seqfile2));
                  pos2++;
                  actualPos2 += moves (temporary2); // Miriam Briskman, 05.17.2020
                  errors++;
                  if (!report1 && pos2a <= r1->end)
                  {
                      if (pos1a == next1a && pos2a == next2a)
                      {
                          next1a = r1->ssa.at(i2++);
                          next2a = r1->ssa.at(i2++);
                          if (i2 >= r1->ssa.size()) 
                             next2a = r1->end + 1;
                          pos1a++; 
                          pos2a++;
                      }
                      else if (pos1a == next1a) 
                          report1 = 1;
                      else if (pos2a == next2a) 
                          pos1a++;
                      else 
                      {
                          pos1a++; 
                          pos2a++;
                      }
            	  }
         	  }
         	  while (pos1 < next1) 
              {
                  fprintf (outfile, "-");
                  // Discard next character:
                  fgetc (seqfile);
                  pos1++;
                  actualPos1 += moves (temporary1); // Miriam Briskman, 05.17.2020
                  errors++;
                  if (!report1 && pos2a <= r1->end)
                  {
                      if (pos1a == next1a && pos2a < next2a) 
                          pos2a++;
                      else 
                          report1 = 1;
                  }
              }
              next1 = r1->ssa.at(i++);
              next2 = r1->ssa.at(i++);
              if (i >= r1->ssa.size())
              {
                  last = true;
                  next2 = next1 + pos2 - pos1;
              }
           }
       }// end inner while loop
       fprintf (outfile, "  %ld\n", actualPos2 - 1);
       if (last && pos1 >= next1)
       {
          reps += float(pos1 - perstart)/lastper;
          done = true;
       }
       else
       {
          lastper = nextper - perstart;
          reps++;
       }
   }// end outer while loop

   period = (r1->length)/reps;
   ts1 = r1->ssa.at(r1->ssa.size() - 2) - r1->start;
   ts2 = r1->end - r1->ssa.at(1) + 1;
   score = ((ts1 > ts2) ? ts1: ts2) - errors * (ERROR_VAL + 1);
   fprintf (outfile, "\nErrors: %ld   Score: %ld\n\n\n", errors, score);
   fprintf (tabfile, "%9ld %9ld %6ld %6.1f %6.1f %4ld %5ld\n", r1->actualStart, r1->actualEnd,
            r1->length, period, reps, errors, score);
}

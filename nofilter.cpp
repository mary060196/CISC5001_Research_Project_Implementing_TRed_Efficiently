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
#include <algorithm>
#include "parameters.h"


using namespace std;
typedef struct{
	long start, end, errs;
   vector <long> ssa;
} rep;

const int STRSIZE=40;
void filter (vector<rep> *rv, FILE *seqfile, FILE *seqfile2, FILE *outfile, FILE *tabfile, int limit);bool matchpoint (rep r1, rep r2, int *pos1, int *pos2);
bool repcomp (rep r1, rep r2);
void convert (rep *r1, FILE *seqfile, FILE *seqfile2, FILE *outfile, FILE *tabfile);

int main (int argc, char *argv[])
{
   FILE *infile, *outfile, *tabfile, *seqfile, *seqfile2;
   char filename[STRSIZE];
   long temp;   rep r;
   vector <rep> rv;

   /*********** FOR INPUT/OUTPUT FILES AS COMMAND LINE ARGUMENTS ************/

   if (argc!=5){
   	printf("File parameters required. See readme for details\n");
      exit(1);
   }
	if ((seqfile=fopen(argv[1],"r"))==NULL)
   		{ printf("Error opening sequence file.\n"); exit(1); }
   if ((seqfile2=fopen(argv[1],"r"))==NULL)
   	{ printf("Error opening sequence file.\n"); exit(1); }
	if ((infile=fopen(argv[2],"r"))==NULL)
   		{ printf("Error opening intermediary file.\n"); exit(1); }
   if ((outfile=fopen(argv[3],"w"))==NULL)
   	{ printf("Error opening alignment file.\n"); exit(1); }
   if ((tabfile=fopen(argv[4],"w"))==NULL)
   	{ printf("Error opening output table file.\n"); exit(1); }

   /**************** TO MANUALLY ENTER INPUT/OUTPUT FILES *******************/

/*   printf("Enter sequence file: ");
   cin >> filename;
	if ((seqfile=fopen(filename,"r"))==NULL)
   		{ printf("Error opening sequence file.\n"); exit(1); }   if ((seqfile2=fopen(filename,"r"))==NULL)
   		{ printf("Error opening sequence file.\n"); exit(1); }

   printf("Enter intermediary file: ");
   cin >> filename;
	if ((infile=fopen(filename,"r"))==NULL)
   		{ printf("Error opening intermediary file.\n"); exit(1); }
	printf("Enter output alignment file: ");   cin >> filename;
   if ((outfile=fopen(filename,"w"))==NULL)   	{ printf("Error opening alignment file.\n"); exit(1); }

	printf("Enter output table file: ");
   cin >> filename;
   if ((tabfile=fopen(filename,"w"))==NULL)   	{ printf("Error opening output table file.\n"); exit(1); }

   /************************************************************************/

   fprintf(tabfile,"StartIndx  EndIndex Length Period   Reps Errs Score\n");

   while (fscanf(infile, "%ld", &temp)!=EOF){
   	if (temp==-1) {
         fscanf(infile, "%ld", &temp);
      	filter(&rv,seqfile,seqfile2,outfile,tabfile,temp);
      }
   	else {
      	rv.push_back(rep());
			rv.back().start = temp;
      	fscanf(infile, "%ld %ld", &(rv.back().end), &(rv.back().errs));
      	do {
            fscanf(infile, "%ld", &temp);
            rv.back().ssa.push_back(temp);
            fscanf(infile, "%ld", &temp);
            rv.back().ssa.push_back(temp);
      	} while (rv.back().ssa.back()!=0);
   	}
   };

	return 0;
}

void filter (vector<rep> *rv, FILE *seqfile, FILE *seqfile2, FILE *outfile, FILE *tabfile, int limit)
{
/*	int pos1, pos2;

   sort(rv->begin(),rv->end(),repcomp);

   for(int i=0; i<rv->size(); i++){
   	for(int j=i+1; j<rv->size(); ){
      	if (matchpoint(rv->at(i),rv->at(j),&pos1,&pos2)){
         	if (rv->at(i).end < rv->at(j).end){
            	rv->at(i).ssa.erase(rv->at(i).ssa.begin()+pos1,rv->at(i).ssa.end());
               while (pos2 < rv->at(j).ssa.size()){
               	rv->at(i).ssa.push_back(rv->at(j).ssa.at(pos2));
                  pos2++;
               }
               rv->at(i).end=rv->at(j).end;
            }
            rv->erase(rv->begin()+j);
         }
         else j++;
      }
   }                         */

   if (limit==-1){
      for(int i=0; i<rv->size(); ){
      	//fprintf(outfile,"%ld %ld %d\n",rv->at(i).start,rv->at(i).end,rv->at(i).errs);
         //for(int j=0; j<rv->at(i).ssa.size(); j++)
         //	fprintf(outfile,"%ld ",rv->at(i).ssa.at(j));
         //fprintf(outfile,"\n\n");
         rep newrepeat = *(rv->begin()+i);
         convert(&newrepeat, seqfile, seqfile2, outfile, tabfile);         rv->erase(rv->begin()+i);
      }
   }
   else {
   	for(int i=0; i<rv->size(); ){
   		if (rv->at(i).end<limit){
      		//fprintf(outfile,"%ld %ld %d\n",rv->at(i).start,rv->at(i).end,rv->at(i).errs);
            //for(int j=0; j<rv->at(i).ssa.size(); j++)
            //	fprintf(outfile,"%ld ",rv->at(i).ssa.at(j));
            //fprintf(outfile,"\n\n");
            convert(&(*(rv->begin()+i)), seqfile, seqfile2, outfile, tabfile);
         	rv->erase(rv->begin()+i);
      	}
      	else i++;
      }
   }
}

bool matchpoint (rep r1, rep r2, int *pos1, int *pos2)
{
   // Need to iterate through vectors in r1 and r2. For alignment at position in
   // r1, if position in r2 has same difference and falls before next pair in r1,
   // then there is a match. Must save and return positions also for use in main
   // program.
   int i=0, j=0;
   int i1, i2, j1, j2, k1, k2, l1, l2;

   while(i<r1.ssa.size()/2-1 && j<r2.ssa.size()/2-1){
      i1 = r1.ssa.at(i*2);
      i2 = r1.ssa.at(i*2+1);
      j1 = r2.ssa.at(j*2);
      j2 = r2.ssa.at(j*2+1);
      k1 = r1.ssa.at(i*2+2);
      k2 = r1.ssa.at(i*2+3);
      l1 = r2.ssa.at(j*2+2);
      l2 = r2.ssa.at(j*2+3);

      if (k2==0){
      	k1++;
         k2=k1+(i2-i1);
      }

      if (l2==0){
      	l1++;
         l2=l1+(j2-j1);
      }

      if (j1>=k1 || j2>=k2) i++;
      else if (j2-j1 == i2-i1){
      	*pos1 = i*2+2;
         *pos2 = j*2+2;
         return 1;
      }
      else {
         if (l1>k1 && l2>k2 && (j2-j1 == k2-k1)){
         	*pos1 = i*2+4;
            *pos2 = j*2+2;
            return 1;
         }
      	j++;
      }
   }
	return 0;
}

bool repcomp (rep r1, rep r2)
{
	if (r1.start<r2.start) return 1;
   if (r1.start>r2.start) return 0;
   if (r1.end>r2.end) return 1;
   return 0;
}

void convert (rep *r1, FILE *seqfile, FILE *seqfile2, FILE *outfile, FILE *tabfile)
{
	// print header
	fprintf(outfile, "Start: %ld End: %ld Length: %ld\n\n",r1->start,r1->end,r1->end-r1->start+1);

   // declare and initialize
	long pos1, pos1a, pos2, pos2a, next1, next1a, next2, next2a;
   long nextper, perstart, lastper;
   long errors = 0;
   long i = 4, i2 = 0;
   char bp;
   bool done = 0, last = 0;
   int report1 = 2;
   float reps = 1, period;
   long score;
   pos1 = r1->ssa.at(0);
   pos2 = r1->ssa.at(1);
   fseek(seqfile, (pos1-START_POS)*sizeof(char) ,0);
   fseek(seqfile2, (pos2-START_POS)*sizeof(char), 0);
   next1 = r1->ssa.at(2);
   next2 = r1->ssa.at(3);
   if (r1->ssa.size()==4){
   	last = 1;
		next2 = next1+pos2-pos1;
   }
   lastper = pos2-pos1;

   // print repeat
   while (!done){   // outer loop
   	if (report1){
      	if (report1==2)
         	fprintf(outfile,"%ld  ",pos1);
         else {
         	for (long j=1; j<=pos1; j*=10)
         		fprintf(outfile," ");
         	fprintf(outfile,"  ");
         }
         pos1a = pos1;
         pos2a = pos2;
         i2 = i-4;
         while (r1->ssa.at(i2)>pos1) i2 -= 2;
         i2 += 2;
         next1a = r1->ssa.at(i2++);
         next2a = r1->ssa.at(i2++);
         if (i2>=r1->ssa.size()) next2a = r1->end+1;

         while (pos1a<pos2 && pos2a<=r1->end){
         	if (pos1a<next1a && pos2a<next2a){
         		fprintf(outfile,"%c",fgetc(seqfile));
               pos1a++;
               pos2a++;
            }
            else if (pos1a<next1a){
            	fprintf(outfile, "%c", fgetc(seqfile));
               pos1a++;
            }
            else if (pos2a<next2a){
            	fprintf(outfile, "-");
               pos2a++;
            }
            else {
               next1a = r1->ssa.at(i2++);
               next2a = r1->ssa.at(i2++);
               if (i2>=r1->ssa.size()) next2a = r1->end+1;
            }
         }
         fseek(seqfile, (pos1-START_POS)*sizeof(char) ,0);
         if (report1==2)
         	fprintf(outfile,"  %ld",pos2-1);
         fprintf(outfile,"\n");
         report1 = 0;
      } // end if(report1)

      perstart = pos1;
      nextper = pos2;
      fprintf(outfile, "%ld  ",pos2);

   	while (pos1<nextper){     // inner loop
      	if (pos1<next1 && pos2<next2){
         	bp = getc(seqfile2);
            fprintf(outfile, "%c", bp);
            if (bp != fgetc(seqfile)) errors++;
            pos1++;
            pos2++;
         	if (!report1 && pos2a <= r1->end ){
            	if (pos1a==next1a && pos2a==next2a){
               	next1a = r1->ssa.at(i2++);
               	next2a = r1->ssa.at(i2++);
               	if (i2>=r1->ssa.size()) next2a = r1->end+1;
                  pos1a++; pos2a++;
               }
            	else if (pos1a==next1a) report1=1;
               else if (pos2a==next2a) pos1a++;
               else {pos1a++; pos2a++;}
            }
         }
         else {
         	if (last) {
            	done = 1;
               break;
            }
         	while (pos2<next2) {
         		fprintf(outfile, "%c", fgetc(seqfile2));
            	pos2++;
               errors++;
               if (!report1 && pos2a <= r1->end ){
               	if (pos1a==next1a && pos2a==next2a){
               		next1a = r1->ssa.at(i2++);
               		next2a = r1->ssa.at(i2++);
               		if (i2>=r1->ssa.size()) next2a = r1->end+1;
                     pos1a++; pos2a++;
                  }
                  else if (pos1a==next1a) report1=1;
                  else if (pos2a==next2a) pos1a++;
               	else {pos1a++; pos2a++;}
            	}
         	}
         	while (pos1<next1) {
            	fprintf(outfile, "-");
               fgetc(seqfile);
            	pos1++;
               errors++;
               if (!report1 && pos2a <= r1->end ){
                  if (pos1a==next1a && pos2a<next2a) pos2a++;
                  else report1=1;
            	}
            }
            next1 = r1->ssa.at(i++);
      		next2 = r1->ssa.at(i++);
      		if (i>=r1->ssa.size()){
            	last = 1;
               next2 = next1+pos2-pos1;
            }
         }
      }// end inner while loop
      fprintf(outfile, "  %ld\n",pos2-1);
      if (last && pos1>=next1){
         reps += float(pos1-perstart)/lastper;
      	done = 1;
      }
      else{
      	lastper = nextper-perstart;
      	reps++;
      }
   }// end outer while loop

   period = (r1->end-r1->start+1)/reps;
   long ts1 = r1->ssa.at(r1->ssa.size()-2) - r1->start;
   long ts2 = r1->end - r1->ssa.at(1) +1;
   score = ((ts1>ts2)? ts1:ts2) - errors * (ERROR_VAL+1);
   fprintf(outfile, "\nErrors: %ld   Score: %ld\n\n\n",errors,score);
   fprintf(tabfile,"%9ld %9ld %6ld %6.1f %6.1f %4ld %5ld\n",r1->start,r1->end,r1->end-r1->start+1,period,reps,errors,score);
}


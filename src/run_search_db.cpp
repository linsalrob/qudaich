/***********************************************************************************

    Copyright (C) 2013 by Sajia Akhter, Edwards Lab, San Diego State University

    This file is part of Qudaich.

    Qudaich is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

**********************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h> /* required for MAX_PATH */


void version(char *);
char *execpath(char *);
void inputError(const char *);
void help(char *);

int TOTAL_DB_SEQ = 0, TOTAL_QRY_SEQ = 0, LINE_LEN, SIZE_SEQ = 0, SIZE_SEQ_NAME = 0;
char mypath[PATH_MAX];


void inputError(const char *msg) {
  fprintf(stderr, "%s", msg);
  exit(-1);
}


/* return the path of the current executable so we can append bin to it
 * arg: exec is argv[0]
 * returns: a string representing the current path
 *
 * Uses realpath in stdlib.h and MAX_PATH in limits.h
 */
char *execpath(char *exec) {
	char *mp = realpath(exec, mypath); /*  this is the complete path */
	char *ptr = strrchr(mp, '/'); /* this is from the last / forwards - should be path separator! */
	int posn = strlen(mp) - strlen(ptr); /* the position of the last / */
	mypath[posn] = '\0';
	return mypath;
}

/* print out the version number and then the help menu. Note that we require argv[0] so we can pass it forward */
void version(char *p) {
	FILE * fp;
	char v[256];
	fp = fopen("VERSION", "r");
	if (fp != NULL) {
		char *check_fgets = fgets(v, 256, fp);
	}
	printf("\n%s %s\n", p, v);
	help(p);
}



/* help() prints the help menu and exits. Please provide argv[0] as an option so we can include that in the help menu */
void help(char * p) {
      printf("\n%s -q [query] -d [database] -p [program] or\n%s -query [query file] -db [database file] -prog [program]\n\n", p, p);
      printf("-q|-query                      Name of the query file (Required)\n");
      printf("-d|-db                         Name of database file (Required)\n");
      printf("-p|-prog                       Alignment options (Required):  n (nucleotide),\n"); 
      printf("                                                           p (protein),\n");
      printf("                                                           trn (translated nucleotide)\n");
      printf("                                                           trnx (translated nucleotide vs protein database)\n");
      printf("-t|-top                        Number of alignments per query sequence (default 1)\n");
      printf("-f|-freqFile                   Frequency file Name\n");
      printf("-c|-heuristic                  Heuristic options: 1 (default) or 2 (See the README or paper for more information)\n");
      printf("-h|-help                       Show command line options and exit\n");
      printf("-v|-version                    Print the version and help options and exit\n");
      exit(0);
}



int input_size(FILE *f, int flag)
{
  int L = 2000, SIZE = 0, t, max = 0, total_seq = 0, temp = 0, maxsize=0 ;
  char name[L];
  
  while(fgets(name,L,f)!=NULL)
  {
        if (name[0] == '>'){
		total_seq++;
		t = strlen(name);
		if(t>max) max = t; 
		SIZE += temp;
		if (temp >  maxsize) maxsize = temp;
		temp = 0;
                continue;
	}
        else
	{	
		t = strlen(name)-1;
		if((name[t]) != '\n')
            		t++;
                temp = temp + t;
	}
  }
  SIZE += temp; 
  if (temp >  maxsize) maxsize = temp;

  if (flag == 0){//DB
	TOTAL_DB_SEQ = total_seq;
        if (SIZE_SEQ<maxsize) SIZE_SEQ = maxsize;
	if (SIZE_SEQ_NAME<max) SIZE_SEQ_NAME = max;
  }
  else//query
  {
 	TOTAL_QRY_SEQ = total_seq;
        if (SIZE_SEQ<maxsize) SIZE_SEQ = maxsize;
        if (SIZE_SEQ_NAME<max) SIZE_SEQ_NAME = max;
  }
  return SIZE;
}

int main(int argc, char **argv) {
  int i;
  char queryFile[1024], refFile[1024], outputfile[1024], program[16];
  int top = 1, hypothesis = 1;
  FILE *f;

  strcpy(outputfile, "freq.txt");
  
  for(i=1;i<argc;i++) {
    if(!strcmp(argv[i], "-query") || !strcmp(argv[i], "-q")) {
      if(i+1 >= argc) {
	inputError("Query File name required\n");
      }
      strcpy(queryFile, argv[i+1]);
      i++;
    }
    else if(!strcmp(argv[i], "-ref") || !strcmp(argv[i], "-db") || !strcmp(argv[i], "-d")) {
      if(i+1 >= argc) {
	inputError("Database File name required\n");
      }
      strcpy(refFile, argv[i+1]);
      i++;
    }
    else if(!strcmp(argv[i], "-prog") || !strcmp(argv[i], "-p")) {
      if(i+1 >= argc) {
        inputError("Alignment option required\n");
      }
      strcpy(program, argv[i+1]);
      if (strcmp(program,"n")!=0 and strcmp(program,"p")!=0 and strcmp(program,"trn")!=0 and strcmp(program,"trnx")!=0)
	inputError("Valid alignment option required\nPlease use -h for details information\n");
      i++;
    }

    else if(!strcmp(argv[i], "-freqFile") || !strcmp(argv[i], "-f")) {
      if(i+1 >= argc) {
        inputError("Frequency File name required\n");
      }
      strcpy(outputfile, argv[i+1]);
      i++;
    }
    else if(!strcmp(argv[i], "-hypo") || !strcmp(argv[i], "-heuristic") || !strcmp(argv[i], "-c")) {
      if(i+1 >= argc) {
        inputError("Argument required\n");
      }
      sscanf(argv[i+1], "%d", &hypothesis);
      if(hypothesis <= 0 || hypothesis >=3) {
        inputError("Value of the heuristic must be either 1 (default) or 2\n");
      }
      i++;
    }
    else if(!strcmp(argv[i], "-top") || !strcmp(argv[i], "-t")) {
      if(i+1 >= argc) {
	inputError("Argument required\n");
      }
      sscanf(argv[i+1], "%d", &top);
      if(top <= 0) {
	inputError("Value of top must be positive\n");
      }
      i++;
    }
    else if(!strcmp(argv[i], "-h") || !strcmp(argv[i], "-help"))
	    help(argv[0]);
    else if (!strcmp(argv[i], "-version") || !strcmp(argv[i], "-v"))
	    version(argv[0]);
    else {
      inputError("Invalid arguments\n");
      inputError(argv[i]);	
    }
  }
    
  if(!queryFile[0] || !refFile[0] || !program[0]) {
    help(argv[0]);
  }

  // CHECK whether these are valid files
  unsigned int SIZEr = 0, SIZEq = 0, SIZE = 0;
  
  //read file:
  f = fopen(refFile,"r");
  SIZEr = input_size(f,0);
  fclose(f);
  f = fopen(queryFile,"r");
  SIZEq = input_size(f,1);
  fclose(f);
	
  char progBuf[1024],write_buff[1024];
  LINE_LEN = 2000;

  SIZE_SEQ++;
  SIZE_SEQ_NAME++;
  sprintf(write_buff,"%d %d %d %d %d", TOTAL_DB_SEQ, TOTAL_QRY_SEQ, LINE_LEN, SIZE_SEQ, SIZE_SEQ_NAME);

  char *app_path;
  if ((app_path = execpath(argv[0])) == NULL) {
	  fprintf(stderr, "ERROR: We can not get an application path from %s\n", argv[0]);
	  exit(-1);
  }
  

  
  if (strcmp(program,"n")==0)
  {
	SIZE = SIZEr * 2 + SIZEq + 2;
        if (SIZE >= 2000000000){
		 printf("Current version of Quadich can handle maximum 2 billion bp.\nThe input database and query file contain %u bps [= total database bp (%u)  * 2 + total query bp (%u)]\n",SIZE,SIZEr,SIZEq);
		inputError("");
	}

	if(hypothesis == 1)
	{ 
		if(top == 1)
			sprintf(progBuf, "%s/bin/search_dna_hypo1 %s %s %s %d %d %s", app_path, refFile, queryFile, outputfile, SIZE, top, write_buff);
		else	
			sprintf(progBuf, "%s/bin/search_dna_hypo1_top %s %s %s %d %d %s", app_path, refFile, queryFile, outputfile, SIZE, top, write_buff);
	}
	else
	{
	        if(top == 1)
                        sprintf(progBuf, "%s/bin/search_dna_hypo2 %s %s %s %d %d %s", app_path, refFile, queryFile, outputfile, SIZE, top, write_buff);
                else
                        sprintf(progBuf, "%s/bin/search_dna_hypo2_top %s %s %s %d %d %s", app_path, refFile, queryFile, outputfile, SIZE, top, write_buff);
	}
  }
  else if (strcmp(program,"p")==0)
  {
        SIZE = SIZEr + SIZEq + 2;
        if (SIZE >= 2000000000){ 
                          printf("Current version of Quadich can handle maximum 2 billion bp.\nThe input database and query file contain %u bps [= total database bp (%u) + total query bp (%u)]\n",SIZE,SIZEr,SIZEq);
                inputError("");
        }
	if(hypothesis == 1)
        { 
                if(top == 1)
                        sprintf(progBuf, "%s/bin/search_pro_hypo1 %s %s %s %d %d %s", app_path, refFile, queryFile, outputfile, SIZE, top, write_buff);
                else
                        sprintf(progBuf, "%s/bin/search_pro_hypo1_top %s %s %s %d %d %s", app_path, refFile, queryFile, outputfile, SIZE, top, write_buff);
        }
        else
        {
                if(top == 1)
                        sprintf(progBuf, "%s/bin/search_pro_hypo2 %s %s %s %d %d %s", app_path, refFile, queryFile, outputfile, SIZE, top, write_buff);
                else
                        sprintf(progBuf, "%s/bin/search_pro_hypo2_top %s %s %s %d %d %s", app_path, refFile, queryFile, outputfile, SIZE, top, write_buff);
        }
  }
  else if (strcmp(program,"trn")==0)
  {
        SIZE = SIZEr * 2 + SIZEq + 2;
        if (SIZE >= 2000000000) {
                          printf("Current version of Quadich can handle maximum 2 billion bp.\nThe input database and query file contain %u bps [= total database bp (%u)  * 2 + total query bp (%u)]\n",SIZE,SIZEr,SIZEq);
                inputError("");
        }
	if(hypothesis == 1)
        {
                if(top == 1)
                        sprintf(progBuf, "%s/bin/search_trn_hypo1 %s %s %s %d %d %s", app_path, refFile, queryFile, outputfile, SIZE, top, write_buff);
                else
                        sprintf(progBuf, "%s/bin/search_trn_hypo1_top %s %s %s %d %d %s", app_path, refFile, queryFile, outputfile, SIZE, top, write_buff);
        }
        else
        {
                if(top == 1)
                        sprintf(progBuf, "%s/bin/search_trn_hypo2 %s %s %s %d %d %s", app_path, refFile, queryFile, outputfile, SIZE, top, write_buff);
                else
                        sprintf(progBuf, "%s/bin/search_trn_hypo2_top %s %s %s %d %d %s", app_path, refFile, queryFile, outputfile, SIZE, top, write_buff);
        }
  }

//////////////////
  else if (strcmp(program,"trnx")==0)
  {
        SIZE = SIZEr + SIZEq*2 + 2;
        if (SIZE >= 2000000000) {
                          printf("Current version of Quadich can handle maximum 2 billion bp.\nThe input database and query file contain %u bps [= total database bp (%u) + total query bp (%u) * 2]\n",SIZE,SIZEr,SIZEq);
                inputError("");
        }
        if(hypothesis == 1)
        {
                        sprintf(progBuf, "%s/bin/search_trnx_hypo1_top %s %s %s %d %d %s", app_path, refFile, queryFile, outputfile, SIZE, top, write_buff);
        }
        else
        {
		printf("Current version of Quadich supports only hypothesis I for translated nucleotide vs protein database alignmnet.\n\n");
                inputError("");
                /*if(top == 1)
                        sprintf(progBuf, "%s/bin/search_trn_hypo2 %s %s %s %d %d %s", app_path, refFile, queryFile, outputfile, SIZE, top, write_buff);
                else
                        sprintf(progBuf, "%s/bin/search_trn_hypo2_top %s %s %s %d %d %s", app_path, refFile, queryFile, outputfile, SIZE, top, write_buff);
      
	*/
        }
  }
  else{
        inputError("Please specify the correct program option\n");
      }


  printf("%s\n",progBuf);
  int returnval = system(progBuf);

  return returnval;
}


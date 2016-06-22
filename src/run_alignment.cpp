/*****************************************************************************

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

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

void inputError(const char *msg) {
  fprintf(stderr, "%s", msg);
  exit(-1);
}

int main(int argc, char **argv) {

  char dummy[1024], queryFile[1024] = "", outputfile[1024] = "output_qudaich.txt", freqFile[1024]= "freq.txt", program[16] = "", fvalue[64] = "avg";
  int top = 1, freq;
  int GAP_OPENING_COST = -1, GAP_EXTENSION_COST = -2, MATCH_WT = 1, MISMATCH_WT = -3;

  FILE *f;
  

  for(int i=1;i<argc;i++) {
   
    if(!strcmp(argv[i], "-f")) {
      if(i+1 >= argc) {
        inputError("Accepted argument for -f is: avg, all or a INT value\n");
      }
      strcpy(fvalue, argv[i+1]);
      if( strcmp(fvalue,"avg") != 0 and strcmp(fvalue,"all") != 0)
	{
		for (unsigned j = 0; j< strlen(fvalue);j++)
		{
			if(isdigit(fvalue[j]))	continue;
			else inputError("Accepted argument for -f is: avg, all or a INT value\n");
		}
	}
      i++;
    }

    else if(!strcmp(argv[i], "-freqFile")) {
      if(i+1 >= argc) {
	inputError("Frequency File name required\n");
      }
      strcpy(freqFile, argv[i+1]);
      i++;
    }
    else if(!strcmp(argv[i], "-output")) {
      if(i+1 >= argc) {
        inputError("Output File name required\n");
      }
      strcpy(outputfile, argv[i+1]);
      i++;
    }
    else if(!strcmp(argv[i], "-match")) {
      if(i+1 >= argc) {
	inputError("Argument required\n");
      }
      sscanf(argv[i+1], "%d", &MATCH_WT);
      i++;
    }
    else if(!strcmp(argv[i], "-mismatch")) {
      if(i+1 >= argc) {
        inputError("Argument required\n");
      }
      sscanf(argv[i+1], "%d", &MISMATCH_WT);
      i++;
    }
    else if(!strcmp(argv[i], "-gap_open")) {
      if(i+1 >= argc) {
        inputError("Argument required\n");
      }
      sscanf(argv[i+1], "%d", &GAP_OPENING_COST);
      i++;
    }
    else if(!strcmp(argv[i], "-gap_ext")) {
      if(i+1 >= argc) {
        inputError("Argument required\n");
      }
      sscanf(argv[i+1], "%d", &GAP_EXTENSION_COST);
      i++;
    }

    else if(!strcmp(argv[i], "-h")) {
      printf("-f                           Options: all = generate alignments for all query sequences\n");
      printf("                                      avg (default) = generate alignments for those query sequences whose frequency or sum(lcp) >= average of all query sequences\n");
      printf("                                      an integer value = generate alignments for those query sequences whose frequency or sum(lcp) >= given integer value\n");
      printf("-freqFile                    Name of the frequency file\n");
      printf("-output                      Name of output file\n");
      printf("-match                       Match weight (default 1)\n");
      printf("-mismatch                    Mismatch penalty (default -3)\n");
      printf("-gap_open                    Gap opening penalty (default -1)\n");
      printf("-gap_ext                     Gap extension penalty (default -2)\n");
      printf("-h                           Show command line options\n");
      exit(0);
    }
    else {
      inputError("Invalid arguments\n");
    }
  }
  
  f = fopen(freqFile,"r");
  if( f == NULL ) {
    inputError("Frequency File not found\n");
  }
  else{ // Read input from fre file
	int ignore;

        ignore = fscanf(f,"%d %d",&freq, &top);
        ignore = fscanf(f,"%s",dummy);
	ignore = fscanf(f,"%s %s %s", queryFile, dummy, program);
	ignore++; /* this is to suppress a compiler warning :) */

	fclose(f);
  }
  
  char progBuf[1024];
  
  if (strcmp(program,"n")==0)
  {
	if(top == 1)
		sprintf(progBuf, "bin/align_dna %s %s %d %d %d %d %s", freqFile, outputfile, MATCH_WT, MISMATCH_WT, GAP_OPENING_COST, GAP_EXTENSION_COST, fvalue);
	else	
		sprintf(progBuf, "bin/align_dna_top %s %s %d %d %d %d %s", freqFile, outputfile, MATCH_WT, MISMATCH_WT, GAP_OPENING_COST, GAP_EXTENSION_COST, fvalue);
  }
  else if (strcmp(program,"p")==0)
  {
        if(top == 1)
                sprintf(progBuf, "bin/align_pro %s %s %d %d %d %d %s", freqFile, outputfile, MATCH_WT, MISMATCH_WT, GAP_OPENING_COST, GAP_EXTENSION_COST, fvalue);
        else
                sprintf(progBuf, "bin/align_pro_top %s %s %d %d %d %d %s", freqFile, outputfile, MATCH_WT, MISMATCH_WT, GAP_OPENING_COST, GAP_EXTENSION_COST, fvalue);
  }
  else if (strcmp(program,"trn")==0)
  {
                if(top == 1)
                        sprintf(progBuf, "bin/align_trn %s %s %d %d %d %d %s", freqFile, outputfile, MATCH_WT, MISMATCH_WT, GAP_OPENING_COST, GAP_EXTENSION_COST, fvalue);
                else
                        sprintf(progBuf, "bin/align_trn_top %s %s %d %d %d %d %s", freqFile, outputfile, MATCH_WT, MISMATCH_WT, GAP_OPENING_COST, GAP_EXTENSION_COST, fvalue);
  }
  else if (strcmp(program,"trnx")==0)
  {
                sprintf(progBuf, "bin/align_trnx_top %s %s %d %d %d %d %s", freqFile, outputfile, MATCH_WT, MISMATCH_WT, GAP_OPENING_COST, GAP_EXTENSION_COST, fvalue);
  }
  else{
	inputError("Something happend in previous step\n");
      }
	

  printf("%s\n",progBuf);
  int returnval = system(progBuf);
  return returnval;
}


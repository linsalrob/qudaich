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
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <sys/time.h>
#include "match_matrix.h"
using namespace std;

#define INF 1000000000
#define EPS       1e-6

#define GO_DIAG  0     // maps [i-1][j-1]
#define GO_UP    1     // maps [i-1][j]
#define GO_LEFT  2     // maps [i][j-1]
#define GAP_OPEN 1
#define GAP_EXT  2
#define ssSIZE_SEQ 1000

#define INDEX(X,Y) ((X)*TOP_FREQ+Y) 


float MATCH_WT = 1;
float MISMATCH_WT = -3;

int TOTAL_DB_SEQ;
int TOTAL_QRY_SEQ;
int SIZE_SEQ_NAME;
int SIZE_SEQ;
int LINE_LEN = 2000;
int ignore;

int THRESHOLD_FREQ = 0;

float GAP_OPENING_COST = -1;
float GAP_EXTENSION_COST = -2;
float NEW_GAP_COST = GAP_OPENING_COST; // + GAP_EXTENSION_COST;

int TOP_FREQ;

char QRY[ssSIZE_SEQ];
char DB[ssSIZE_SEQ];
char PROTEIN[ssSIZE_SEQ];
char output_str[1000]; //to write on file

float sMat[ssSIZE_SEQ+1][ssSIZE_SEQ+1];
float qrMat[ssSIZE_SEQ+1][ssSIZE_SEQ+1];
float dbMat[ssSIZE_SEQ+1][ssSIZE_SEQ+1];

char sMatPath[ssSIZE_SEQ+1][ssSIZE_SEQ+1];
char qrMatPath[ssSIZE_SEQ+1][ssSIZE_SEQ+1];
char dbMatPath[ssSIZE_SEQ+1][ssSIZE_SEQ+1];

char **db_seq_name;
char **db_seq;

char complement[128];
char convertStr[128];
char converter[128];
#include "sequence_processing.h"

int ind;

char inputfile1[1024];
char inputfile2[1024];
char db_id_file[1024];
char output_file[1024];

void initialization(void)
{
  int i;
  //initialization
  db_seq_name = new char* [TOTAL_DB_SEQ];
  db_seq = new char* [TOTAL_DB_SEQ];

  for ( i = 0; i< TOTAL_DB_SEQ; i++)
  {
        db_seq_name[i] = new char [SIZE_SEQ_NAME];
        db_seq[i] = new char [SIZE_SEQ];
  }

  memset(converter,0,128);
  converter['A'] = 0;
  converter['C'] = 1;
  converter['G'] = 2;
  converter['T'] = 3;
  converter['a'] = 0;
  converter['c'] = 1;
  converter['g'] = 2;
  converter['t'] = 3;

  memset(complement,'a',128);
  complement['A'] = 't';
  complement['C'] = 'g';
  complement['G'] = 'c';
  complement['T'] = 'a';
  complement['a'] = 't';
  complement['c'] = 'g';
  complement['g'] = 'c';
  complement['t'] = 'a';

  init_codonList();
  initialize_blosum_62();
}

void make_reverse_complement(int dlen)
{
  int start, end = dlen-1;
  char ch;

  for(start = 0; start<= end ; end--,start++)
    {
      ch = QRY[start];
      QRY[start] = complement[QRY[end]];
      QRY[end] = complement[ch];
    }
}

void translate_dna(int qi, int l)
{
  int i, j, k = strlen(PROTEIN);

  strcpy(QRY,PROTEIN);

  if (qi > 2)
    {
      make_reverse_complement(k);
      qi = qi - 3;
    }
  for(i = qi, j = 0; i< k -2; i= i +3, j++)
          QRY[j] = codonList[converter[QRY[i]]][converter[QRY[i+1]]][converter[QRY[i+2]]];
  QRY[j] = '\0';
}

int input_db_protein(FILE *f)
{
  int i, total_seq = 0, temp_length, n = 0; //DB ID
  char *check_fgets, name[LINE_LEN];

  if(f == NULL){
    printf("No DB file\n");
    exit(0);
  }

  check_fgets = fgets(name,LINE_LEN,f);

  temp_length = strlen(name)-1;
  if((name[temp_length]) == '\n')
       name[temp_length] = 0;
   strcpy(db_seq_name[total_seq],name);

  while ((check_fgets = fgets(name,LINE_LEN,f)) != NULL)
    {
      if (name[0] == '>')
        {
          total_seq++;
          temp_length = strlen(name)-1;
          if((name[temp_length]) == '\n')
                 name[temp_length] = 0;
          strcpy(db_seq_name[total_seq],name);

          n = 0;
        }
      else
        {
          int t = strlen(name)-1;
          if((name[t]) == '\n')
            name[t] = 0;
          strcpy(db_seq[total_seq]+ n,name);
          n += t;
        }
    }

  total_seq++;
  printf("DONE this function: %d %d\n", n, total_seq);
  fflush(stdout);
  return total_seq;
}


int input_db_not_used(FILE *f)
{
  int i, total_seq = 0, temp_length, n = 0; //DB ID
  char *check_fgets, name[LINE_LEN];

  if(f == NULL){
    printf("No DB file\n");
    exit(0);
  }

  check_fgets = fgets(name,LINE_LEN,f);

  temp_length = strlen(name)-1;
  if((name[temp_length]) == '\n')
       name[temp_length] = 0;
   strcpy(db_seq_name[total_seq],name);

  while ((check_fgets = fgets(name,LINE_LEN,f)) !=NULL)
    {
      if (name[0] == '>')
        {
          total_seq++;
	  temp_length = strlen(name)-1;
  	  if((name[temp_length]) == '\n')
      		 name[temp_length] = 0;
   	  strcpy(db_seq_name[total_seq],name);

          n = 0;
        }
      else
        {
          int t = strlen(name)-1;
          if((name[t]) == '\n')
            name[t] = 0;
          strcpy(db_seq[total_seq]+ n,name);
          n += t;
        }
    }

  total_seq++;
  printf("DONE this function: %d %d\n", n, total_seq);
  fflush(stdout);
  return total_seq;
}


int isEqual(double a, double b) {
  if(fabs(a-b) < EPS) return 1;
  return 0;
}

float _max(double x, double y) {
  return x > y ? x : y;
}

float _max(double x, double y, double z) {
  return x > y ? _max(x, z) : _max(y, z);
}


float similarity_score(char a, char b) {
  return blosum62[convert_blosum62[a]][convert_blosum62[b]];
}

void SW(int flag_for_rev) {
  
  //strcpy(DB,"aatggagaggagtggattggaatggaatggaatggaatcgaatcgaatggaatgctagaaaagaatatgtgttgagattgagcca");
  //strcpy(QRY,"gaacaacaatagttgatggaacaacaatagaatggaatggaatggaatgaaatatccggagtatttggtttgttggatggaacaacgaataggttgatggaacaacaatag");
  //convert_similar_seq();
  int n = strlen(QRY) + 1;  // QRY is 'i' in sMat, DB is 'j' in sMat
  int m = strlen(DB)+ 1;

  //printf("QRY= %s len= %d DB= %s len= %d \n", QRY, n, DB, m);

  int i, j, myi, myj;
  float mymax = -INF, t;

  for(i=1;i<n;i++) {
    sMat[i][0] = -INF;    // since no character exist in DB(-1)
    //    qrMat[i][0]= GAP_OPENING_COST + (i-1) * GAP_EXTENSION_COST;
    qrMat[i][0] = 0;      // Leading gaps not counted
    dbMat[i][0]= -INF;    // Impossible case
    qrMatPath[i][0] = GAP_EXT;
    if(i==1)
      qrMatPath[i][0] = GAP_OPEN;
  }
  for(j=1;j<m;j++) {
    sMat[0][j] = -INF;    // since no character exist in QRY(-1)
    //    dbMat[0][j]= GAP_OPENING_COST + (j-1) * GAP_EXTENSION_COST; 
    dbMat[0][j] = 0;      // Leading gaps not counted
    qrMat[0][j]= -INF;    // Impossible case
    dbMatPath[0][j] = GAP_EXT;
    if(j==1)
      dbMatPath[0][j] = GAP_OPEN;
  }

  for (i = 1; i < n; i++) {
    for (j = 1; j < m; j++) {
      sMat[i][j] = sMat[i-1][j-1] + similarity_score(QRY[i-1], DB[j-1]);
      sMatPath[i][j] = GO_DIAG;

      // QR_i matches DB_j, but QR_(i-1) matches a gap
      t = qrMat[i-1][j-1] + similarity_score(QRY[i-1], DB[j-1]);
      if(t > sMat[i][j]) {
	sMat[i][j] = t;
	sMatPath[i][j] = GO_UP;
      }

      // QR_i matches DB_j, but DB_(j-1) matches a gap
      t = dbMat[i-1][j-1] + similarity_score(QRY[i-1], DB[j-1]);
      if(t > sMat[i][j]) {
	sMat[i][j] = t;
	sMatPath[i][j] = GO_LEFT;
      }

      // We are considering local alignment, 
      // if we match QR_i and DB_j, best score is negative
      if(sMat[i][j] < 0) {
      	sMat[i][j] = 0;
      	sMatPath[i][j] = -1;
      }

      // Update qrMat[i][j], QR_i matches with a gap
      // two cases - QR_(i-1) can match with a gap (gap_extension), or QR_(i-1) matches
      // with DB_j and QR_i matchs with a gap opening
      qrMat[i][j] = qrMat[i-1][j] + GAP_EXTENSION_COST;
      qrMatPath[i][j] = GAP_EXT;
      t = sMat[i-1][j] + GAP_OPENING_COST;
      if(t > qrMat[i][j]) {
	qrMat[i][j] = t;
	qrMatPath[i][j] = GAP_OPEN;
      }

      // Same for dbMat
      dbMat[i][j] = dbMat[i][j-1] + GAP_EXTENSION_COST;
      dbMatPath[i][j] = GAP_EXT;
      t = sMat[i][j-1] + GAP_OPENING_COST;
      if(t > dbMat[i][j]) {
	dbMat[i][j] = t;
	dbMatPath[i][j] = GAP_OPEN;
      }
      
      if(mymax < sMat[i][j]) {
	mymax = sMat[i][j];
	myi = i;
	myj = j;
      }
    }
  }
  // printf("\nscore = %.2lf pos (index starts from 0) %d %d\n", mymax, myi - 1, myj - 1);

  // CAUTION: GAP_OPEN+GAP_EXT > LOWEST_MISMATCH
  //  Find the alignment
  string sQRY = QRY, sDB = DB;
  int tempgap = 0;

  i = myi, j = myj;
  while(sMat[i][j]) {
    assert(sMatPath[i][j] != -1);

    // Change to qrMat
    if(sMatPath[i][j] == GO_UP) {
      // QR_(i-1) matches a gap
      i--; j--;
      while(qrMatPath[i][j] == GAP_EXT) {
  	sDB.insert(j, 1, '-');
	tempgap++;
  	i--;
      }
      assert(qrMatPath[i][j] == GAP_OPEN);
      sDB.insert(j, 1, '-'); 
      tempgap++;
      i--;
    }
    else if(sMatPath[i][j] == GO_LEFT) {
      // DB_(j-1) matches a gap
      i--; j--;
      while(dbMatPath[i][j] == GAP_EXT) {
  	sQRY.insert(i, 1, '-');
  	j--;
      }
      assert(dbMatPath[i][j] == GAP_OPEN);
      sQRY.insert(i, 1, '-');
      j--;
    }
    else {
      i--;
      j--;
    }
  }

  /////////////////////////my part//////////////////////////////////////
  
  int length, gap=0, similarity=0, identity=0, tempi,tempj, qstart = i, dbstart = j;


  tempj = j; tempi = i;

  while(sQRY[tempi] == '-' or sDB[tempj] == '-')
    {
      tempi++;tempj++;
    }

  if(sQRY[i] == '-')
    dbstart = tempj;
  else
    qstart = tempi;
  
  length = myj - tempj +tempgap; 
  


  for (; tempj < myj + tempgap; tempj++, tempi++)
    {
      if(sQRY[tempi] == '-' or sDB[tempj] == '-')
	gap++;

      else if (sQRY[tempi] == sDB[tempj])
	identity++;
    }

  //printf("score = %.2lf\nlength = %d\nidentity =  %d\ngap = %d\nqstart = %d, qend = %d\ndbstart = %d, dbend = %d\n", mymax, length, identity, gap, qstart+1, myi, dbstart+1, myj);

  if (flag_for_rev>=0)
        sprintf(output_str,"%lu\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2lf" ,strlen(QRY),identity,length,gap,qstart+1,myi,dbstart+1,myj,mymax);
  else
        sprintf(output_str,"%lu\t%d\t%d\t%d\t%lu\t%lu\t%d\t%d\t%.2lf" ,strlen(QRY),identity,length,gap,strlen(QRY)-(qstart+1)+1,strlen(QRY)-myi+1,dbstart+1,myj,mymax);

  ////////// end of my part ///////////////////////////////////////
  
  if(i < j) {
    sQRY.insert(0, j-i, ' ');
  }
  else {
    sDB.insert(0, i-j, ' ');
  }
}

int main(int argc, char **argv)
{
  FILE *f, *fw;
  int total_db,i,j, n,k, l;
  int *match, *freq, temp_length, fscanf_result;
  char name[LINE_LEN], qry_name[LINE_LEN], dummy[1024], temp_freq_threshold[64];
  char *check_fgets;

  ////////////////////// input from top program //////////////////////////
  
  assert(argc == 8);
  strcpy(db_id_file, argv[1]);
  strcpy(output_file,argv[2]);
  sscanf(argv[3], "%f", &MATCH_WT);
  sscanf(argv[4], "%f", &MISMATCH_WT);
  sscanf(argv[5], "%f", &GAP_OPENING_COST);
  sscanf(argv[6], "%f", &GAP_EXTENSION_COST);
  sscanf(argv[7], "%s", temp_freq_threshold);
  //////////////////////  /////////////////////////
  
  NEW_GAP_COST = GAP_OPENING_COST;
  //printf("%f %f %f %f\n", MATCH_WT, MISMATCH_WT, GAP_OPENING_COST, GAP_EXTENSION_COST);

  ////////////////////////////////////////////////////
  f = fopen(db_id_file, "r");
  fscanf_result = fscanf(f,"%d %d",&THRESHOLD_FREQ, &TOP_FREQ);
  fscanf_result = fscanf(f,"%s",dummy);
  fscanf_result = fscanf(f,"%s %s %s", inputfile2, inputfile1, dummy);
  fscanf_result = fscanf(f,"%d %d %d %d %d", &TOTAL_DB_SEQ, &TOTAL_QRY_SEQ, &LINE_LEN, &SIZE_SEQ, &SIZE_SEQ_NAME);
  fscanf_result++; /* suppressing a compiler error */
  fclose(f);

  ////////////////////////////////////////////////////////
  if(SIZE_SEQ > 1999)
        {
                printf("Current version of qudaich can not process sequence length > 2000\n");
                exit(0);
        }
  ///////////////////////////////////////////////////////

  if(temp_freq_threshold[0] != 'a')
        THRESHOLD_FREQ = atoi(temp_freq_threshold);
  else
  {
        if (strcmp(temp_freq_threshold,"all")==0)
                THRESHOLD_FREQ = 0;
  }

  //////////////////////////////////////////////////
  initialization();

  printf("Start menory initialization\n");
  //Initializing memory
  match = new int[TOTAL_QRY_SEQ * 6 * TOP_FREQ];
  freq = new int[TOTAL_QRY_SEQ * 6 * TOP_FREQ];

  printf("Done menory initialization\n");

  //Input for DB
  f = fopen(inputfile1, "r");
  total_db = input_db_protein(f);
  fclose(f);

  printf("Done DB Input\n");
  fflush(stdout);

  //Input for match

  f = fopen(db_id_file, "r");
  assert(f);
  check_fgets = fgets(name,LINE_LEN,f);
  check_fgets = fgets(name,LINE_LEN,f);
  check_fgets = fgets(name,LINE_LEN,f);
  check_fgets = fgets(name,LINE_LEN,f);

  while ((check_fgets = fgets(name,LINE_LEN,f)) != NULL)
    {
        l = 0;
        while(l<TOP_FREQ){
 	     sscanf(name, "%d\t%d\t%d",&i,&j, &k);
      	     match[INDEX(i,l)] = j;
             freq[INDEX(i,l)] = k;
             l++;
             if (l<TOP_FREQ) 
		     check_fgets = fgets(name,LINE_LEN,f);
        }
    }

  fclose(f);

  printf("Done match input\n");

  //Input for query
  f = fopen(inputfile2, "r");
  fw = fopen(output_file,"w");

  check_fgets = fgets(name,LINE_LEN,f);

  temp_length = strlen(name)-1;
  if((name[temp_length]) == '\n')
        name[temp_length] = 0;
   strcpy(qry_name,name);

  i = -1;
  n = 0;
  strcpy(output_str,"dd");
  
  while ((check_fgets = fgets(name,LINE_LEN,f)) !=NULL)
  {
          if (name[0] == '>')
            {
                /////////////////////// CALL SW ///////////////////////////
		//frame 1
		i++;
  		for (l = 0; l < TOP_FREQ; l++)
		{               
		if (freq[INDEX(i,l)] >= THRESHOLD_FREQ){
		        strcpy(DB,db_seq[match[INDEX(i,l)]]);
			if(!l) strcpy(PROTEIN,QRY);
        		translate_dna(i%6, l);
                                SW(1);
			fprintf(fw,"%s\t%s\t%s\n",qry_name, db_seq_name[match[INDEX(i,l)]], output_str);
		}
		}
		//frame 2
		i++;
		for (l = 0; l < TOP_FREQ; l++)
		{
		if (freq[INDEX(i,l)] >= THRESHOLD_FREQ){
			strcpy(DB,db_seq[match[INDEX(i,l)]]);
        		translate_dna(i%6, l);
                                SW(1);
			fprintf(fw,"%s\t%s\t%s\n",qry_name, db_seq_name[match[INDEX(i,l)]], output_str);
		}
		}
		
		//frame 3
		i++;
		for (l = 0; l < TOP_FREQ; l++)
		{
                if (freq[INDEX(i,l)] >= THRESHOLD_FREQ){
			strcpy(DB,db_seq[match[INDEX(i,l)]]);
		        translate_dna(i%6, l);
                                SW(1);
        		fprintf(fw,"%s\t%s\t%s\n",qry_name, db_seq_name[match[INDEX(i,l)]], output_str);
		}
		}

		///

                //frame 4
                i++;
                for (l = 0; l < TOP_FREQ; l++)
                {
                if (freq[INDEX(i,l)] >= THRESHOLD_FREQ){
                        strcpy(DB,db_seq[match[INDEX(i,l)]]);
                        //if(!l) strcpy(PROTEIN,QRY);
                        translate_dna(i%6, l);
                                SW(-1);
			 fprintf(fw,"%s\t%s\t%s\n",qry_name, db_seq_name[match[INDEX(i,l)]], output_str);
                }
                }
                //frame 5
                i++;
                for (l = 0; l < TOP_FREQ; l++)
                {
                if (freq[INDEX(i,l)] >= THRESHOLD_FREQ){
                        strcpy(DB,db_seq[match[INDEX(i,l)]]);
                        translate_dna(i%6, l);
                                SW(-1);
			fprintf(fw,"%s\t%s\t%s\n",qry_name, db_seq_name[match[INDEX(i,l)]], output_str);
                }
                }

                //frame 6
                i++;
                for (l = 0; l < TOP_FREQ; l++)
                {
                if (freq[INDEX(i,l)] >= THRESHOLD_FREQ){
                        strcpy(DB,db_seq[match[INDEX(i,l)]]);
                        translate_dna(i%6, l);
                                SW(-1);

                        fprintf(fw,"%s\t%s\t%s\n",qry_name, db_seq_name[match[INDEX(i,l)]], output_str);
                }
                }

		/////
  		//if(i == 100000 ) exit(0);//printf("%d\n",i);	
       	        temp_length = strlen(name)-1;
  	        if((name[temp_length]) == '\n')
        		name[temp_length] = 0;
   	        strcpy(qry_name,name);
	
                n = 0;
            }
          else
            {
              int t = strlen(name)-1;
              if((name[t]) == '\n')
                name[t] = 0;
              strcpy((char*)QRY+ n,name);
              n += t;
            }
  }

  //for last sequence
   //frame 1
   i++;
   for (l = 0; l < TOP_FREQ; l++)
                {
                if (freq[INDEX(i,l)] >= THRESHOLD_FREQ){
                        strcpy(DB,db_seq[match[INDEX(i,l)]]);
                        if(!l) strcpy(PROTEIN,QRY);
                        translate_dna(i%6, l);
                                SW(1);
                        fprintf(fw,"%s\t%s\t%s\n",qry_name, db_seq_name[match[INDEX(i,l)]], output_str);
                }
   }
   //frame 2
   i++;
   for (l = 0; l < TOP_FREQ; l++)
   {
                if (freq[INDEX(i,l)] >= THRESHOLD_FREQ){
                        strcpy(DB,db_seq[match[INDEX(i,l)]]);
                        translate_dna(i%6, l);
                                SW(1);
                        fprintf(fw,"%s\t%s\t%s\n",qry_name, db_seq_name[match[INDEX(i,l)]], output_str);
                }
   }

   //frame 3
   i++;
   for (l = 0; l < TOP_FREQ; l++)
   {
                if (freq[INDEX(i,l)] >= THRESHOLD_FREQ){
                        strcpy(DB,db_seq[match[INDEX(i,l)]]);
                        translate_dna(i%6, l);
                                SW(1);

                        fprintf(fw,"%s\t%s\t%s\n",qry_name, db_seq_name[match[INDEX(i,l)]], output_str);
                }
   }

    ///
   //frame 4
   i++;
   for (l = 0; l < TOP_FREQ; l++)
   {
                if (freq[INDEX(i,l)] >= THRESHOLD_FREQ){
                        strcpy(DB,db_seq[match[INDEX(i,l)]]);
                        //if(!l) strcpy(PROTEIN,QRY);
                        translate_dna(i%6, l);
                                SW(-1);
                        fprintf(fw,"%s\t%s\t%s\n",qry_name, db_seq_name[match[INDEX(i,l)]], output_str);
                }
   }
   //frame 5
   i++;
   for (l = 0; l < TOP_FREQ; l++)
   {
         if (freq[INDEX(i,l)] >= THRESHOLD_FREQ){
                        strcpy(DB,db_seq[match[INDEX(i,l)]]);
                        translate_dna(i%6, l);
                                SW(-1);
                        fprintf(fw,"%s\t%s\t%s\n",qry_name, db_seq_name[match[INDEX(i,l)]], output_str);
                }
   }

   //frame 6
   i++;
   for (l = 0; l < TOP_FREQ; l++)
   {
                if (freq[INDEX(i,l)] >= THRESHOLD_FREQ){
                        strcpy(DB,db_seq[match[INDEX(i,l)]]);
                        translate_dna(i%6, l);
                                SW(-1);

                        fprintf(fw,"%s\t%s\t%s\n",qry_name, db_seq_name[match[INDEX(i,l)]], output_str);
                }
   }

   /////

  
  
  
  //  CS//////mithWatermanGotoh gotoh;
  //gotoh.Align();

  fclose(f);
  fclose(fw);
  return 0;
}
                                       

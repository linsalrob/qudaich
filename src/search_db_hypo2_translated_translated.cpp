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

/*This program make 1 suffix array using both query and database sequences. 
Then find the best match by traversing the suffix array file once.   
 */

/*
sequence starts at 0, seq id starts at 0, rev complements are odd DB id

To Do - 
-----
 */


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sais_lcp.h"
#include <sys/time.h>
#include <omp.h>
#include <map>
#include <vector>
#include <algorithm>

using namespace std;

#define LINE_LEN 2000
#define SINGLE_DNA_LEN 5000

char inputfile1[1024]; 
char inputfile2[1024];
char outputfile[1024];


unsigned int SIZE;

int TOP_FREQ;

unsigned char *dna;
int *SA, *LCP;
int (*map_input_seq)[2];
int total_input_seq = 0;
int start_qry_seq;
vector<int> input_query_length;
vector<int> input_DB_length;

char query_seq[LINE_LEN];
char name[LINE_LEN];

char converter[128];
#include "sequence_processing.h"
char convertStr[128];
char complement[128];

char protein[SINGLE_DNA_LEN];


struct frequency_for_db {
  int db_id;
  short freq;
};

struct lcp_db_vec {
  int db_id;
  unsigned short score;
};

double getTime() {
  struct timeval TV;
  struct timezone TZ;
  
  const double kMicro = 1e6;
  const int RC = gettimeofday(&TV, &TZ);
  if(RC == -1) {
    printf("ERROR: Bad call to gettimeofday\n");
    return(-1);
  }
  return( ((double)TV.tv_sec)*kMicro + ((double)TV.tv_usec) );
}  

void make_reverse_complement(int qlen)
{
  int start, end = qlen-1;
  char ch;

  for(start = 0; start<= end ; end--,start++)
    {
      ch = protein[start];
      protein[start] = complement[protein[end]];
      protein[end] = complement[ch];
    }
}

int input(FILE *f, int n, int type) // type = -1 means reverse complemet; otherwise type = 1
{
  int pren = n, i,length_each_seq,j,k,l;
  char *check_fgets;

  if(f == NULL){
    printf("No input file\n");
    exit(0);
  }

  check_fgets = fgets(name,LINE_LEN,f);

  length_each_seq = 0;
  while(fgets(name,LINE_LEN,f)!=NULL)
    {
      if (name[0] == '>')
        {
         for (i = 0, k = pren; i < 3; i++) // i = frame shift, k = actual position of translated protein seq
            {
              for(j = i, l = 0; j< length_each_seq - 2; j = j+3, k++, l++)
                {
                  dna[k] = codonList[converter[protein[j]]][converter[protein[j+1]]][converter[protein[j+2]]];
                  map_input_seq[k][0] = total_input_seq;
                  map_input_seq[k][1] = l;
                }
              total_input_seq++;
              n = n + l;
              if (type == 1) {
                  input_query_length.push_back(l);
              }
	      else
		  input_DB_length.push_back(l);

            }


         if (type == -1)
            {
              make_reverse_complement(length_each_seq);
              for (i = 0, k = n; i < 3; i++) // i = frame shift, k = actual position of translated protein seq
                {
                  for(j = i, l = 0; j< length_each_seq - 2; j = j+3, k++, l++)
                    {
                      dna[k] = codonList[converter[protein[j]]][converter[protein[j+1]]][converter[protein[j+2]]];
                      map_input_seq[k][0] = total_input_seq;
		      map_input_seq[k][1] = l;
                    }
                n = n + l;
                total_input_seq++;
                input_DB_length.push_back(l);
                }
            }

          //next seq
          pren=n;
          length_each_seq = 0;
        }
      else
        {
          int t = strlen(name)-1;
          if((name[t]) == '\n')
            name[t] = 0;
          strcpy((char*)protein + length_each_seq, name); //not protein input 
          length_each_seq += t;
        }
    }

  //for last sequence
  for (i = 0, k = pren; i < 3; i++) // i = frame shift, k = actual position of translated protein seq
     {
          for(j = i, l = 0; j< length_each_seq - 2; j = j+3, k++, l++)
          {
             dna[k] = codonList[converter[protein[j]]][converter[protein[j+1]]][converter[protein[j+2]]];
             map_input_seq[k][0] = total_input_seq;
	     map_input_seq[k][1] = l;
          }
          total_input_seq++;
          n = n + l;
          if (type == 1) {
             input_query_length.push_back(l);
          }
	  else
              input_DB_length.push_back(l);
      }

  if (type == -1)
      {
           make_reverse_complement(length_each_seq);
           for (i = 0, k = n; i < 3; i++) // i = frame shift, k = actual position of translated protein seq
             {
               for(j = i, l = 0; j< length_each_seq - 2; j = j+3, k++, l++)
                 {
                   dna[k] = codonList[converter[protein[j]]][converter[protein[j+1]]][converter[protein[j+2]]];
                   map_input_seq[k][0] = total_input_seq;
		   map_input_seq[k][1] = l;
                 }
               n = n + l;
               total_input_seq++;
               input_DB_length.push_back(l);
             }
      }

  return n;
}

// int SA_sort_freq(frequency_for_db a, frequency_for_db b) {
//   return b.freq - a.freq;
// }

int SA_sort_freq(const void *a, const void *b) {
  return ((frequency_for_db *)b)->freq - ((frequency_for_db *)a)->freq;
}

bool SA_sort_db(struct lcp_db_vec a, struct lcp_db_vec b) {
  return (b.db_id > a.db_id);
}

// sort by each qry
void best_exact_match_dna(int n, FILE *fw2)  {
  int i, j, qry_id, SA_top_db, SA_bottom_db, k;
  unsigned total_freq;
  int total_qry = total_input_seq - start_qry_seq;

  unsigned short  score_up, score_bottom, score_considering_length_up, score_considering_length_bottom; //lcp
  int start_i, stop_i; //lcp

  assert(total_qry == (int)input_query_length.size());
  int *SA_qry_start_index, *SA_qry_cur_index;

  SA_qry_start_index = new int[total_qry];
  SA_qry_cur_index   = new int[total_qry];
  assert(SA_qry_start_index && SA_qry_cur_index);

  SA_qry_start_index[0] = SA_qry_cur_index[0] = 0;
  for(i=1;i<total_qry;i++) { 
    SA_qry_start_index[i] = SA_qry_start_index[i-1] + 2*input_query_length[i-1];
    SA_qry_cur_index[i] = SA_qry_start_index[i];
  }

  lcp_db_vec *SA_db_vec = new lcp_db_vec[SA_qry_start_index[total_qry-1] + 2 * input_query_length[total_qry-1]]; //lcp
  assert(SA_db_vec);
  
  char *temp_check = new char[n+1];
  double tm = 0;

  printf("input parameter n = %d\n",n);
  //tm  = -getTime();

  for(i = 0; i<n;i++)
    temp_check[i] = map_input_seq[SA[i]][0] >= start_qry_seq;
  temp_check[n] = 0;

  //tm += getTime();
  //printf("SA tm1 = %lf\n", tm/1e6);
  
  // storing the immediate upper DB match
  vector<int> SA_cur_qry_grp;
  SA_top_db = -1, SA_bottom_db = -1;
  start_i = -1, stop_i = -1;
 
  //tm = -getTime(); 
  for(i=1;i<n;i++) {
    if(temp_check[i]) {  // Start of a query sequence

      if(i > 1){ 
	SA_top_db = map_input_seq[SA[i-1]][0];
        start_i = i; //lcp
      }
    
      SA_cur_qry_grp.clear();
      while(temp_check[i]) {
	qry_id = map_input_seq[SA[i]][0] - start_qry_seq;
	SA_cur_qry_grp.push_back(qry_id);
	i++;
      }
      
      if(i != n){         // Doesn't end with the query group
	SA_bottom_db = map_input_seq[SA[i]][0];
        stop_i = i-1; //lcp
      }
      
      if(SA_top_db!=-1 && SA_bottom_db!=-1) {        
        /////lcp start
        score_up = LCP[start_i]; 
        if (input_DB_length[SA_top_db] - map_input_seq[SA[start_i-1]][1] < score_up)
          score_up = input_DB_length[SA_top_db] - map_input_seq[SA[start_i-1]][1];
        score_considering_length_up = score_up;

        for(j=0 ;j<(int)SA_cur_qry_grp.size();j++) {
          qry_id = SA_cur_qry_grp[j];

          if (input_query_length[qry_id] - map_input_seq[SA[start_i+j]][1] <score_up)
              score_considering_length_up = input_query_length[qry_id] - map_input_seq[SA[start_i+j]][1];

          SA_db_vec[SA_qry_cur_index[qry_id]].db_id = SA_top_db;
          SA_db_vec[SA_qry_cur_index[qry_id]].score = score_considering_length_up;
          SA_qry_cur_index[qry_id]++;

          if(LCP[start_i+j+1] < score_up)
             score_up = LCP[start_i+j+1];
          score_considering_length_up = score_up;
        }     
        ////////////////////////////////////////////

	score_bottom = LCP[stop_i+1];
        if (input_DB_length[SA_bottom_db] - map_input_seq[SA[stop_i+1]][1] < score_bottom)
          score_bottom = input_DB_length[SA_bottom_db] - map_input_seq[SA[stop_i+1]][1];
        score_considering_length_bottom = score_bottom;

	for(j=(int)SA_cur_qry_grp.size()-1 , k = 0;j>-1;j--, k++) {
	  qry_id = SA_cur_qry_grp[j];
	  
         if (input_query_length[qry_id] - map_input_seq[SA[stop_i-k]][1] < score_bottom)
              score_considering_length_bottom = input_query_length[qry_id] - map_input_seq[SA[stop_i-k]][1];

	  SA_db_vec[SA_qry_cur_index[qry_id]].db_id = SA_bottom_db;
          SA_db_vec[SA_qry_cur_index[qry_id]].score = score_considering_length_bottom;
          SA_qry_cur_index[qry_id]++;
	  
          if(LCP[stop_i-k] < score_bottom)
             score_bottom = LCP[stop_i-k];
          score_considering_length_bottom = score_bottom;
	}
      }
      else if(SA_top_db!=-1) {        // Qry group not followed by a DB
	score_up = LCP[start_i];
	if (input_DB_length[SA_top_db] - map_input_seq[SA[start_i-1]][1]  < score_up)
          score_up = input_DB_length[SA_top_db] - map_input_seq[SA[start_i-1]][1];
        score_considering_length_up = score_up;


      	for(j=0;j<(int)SA_cur_qry_grp.size();j++) {
          qry_id = SA_cur_qry_grp[j];
      	 
      	  if (input_query_length[qry_id] - map_input_seq[SA[start_i+j]][1] <score_up)
              score_considering_length_up = input_query_length[qry_id] - map_input_seq[SA[start_i+j]][1];

          SA_db_vec[SA_qry_cur_index[qry_id]].db_id = SA_top_db;
          SA_db_vec[SA_qry_cur_index[qry_id]].score = score_considering_length_up;
          SA_qry_cur_index[qry_id]++;
     
          if(LCP[start_i+j+1] < score_up)
             score_up = LCP[start_i+j+1];
          score_considering_length_up = score_up;
	  
      	}
      }
      else if(SA_bottom_db!=-1) {     // Qry group not preceded by a DB
	score_bottom = LCP[stop_i+1];
        if (input_DB_length[SA_bottom_db] - map_input_seq[SA[stop_i+1]][1] < score_bottom)
          score_bottom = input_DB_length[SA_bottom_db] - map_input_seq[SA[stop_i+1]][1];
        score_considering_length_bottom = score_bottom;

      	for(j=(int)SA_cur_qry_grp.size()-1, k = 0;j>-1;j--, k++) {
      	  qry_id = SA_cur_qry_grp[j];
      	  
          if (input_query_length[qry_id] - map_input_seq[SA[stop_i-k]][1] <score_bottom)
              score_considering_length_bottom = input_query_length[qry_id] - map_input_seq[SA[stop_i-k]][1];

          SA_db_vec[SA_qry_cur_index[qry_id]].db_id = SA_bottom_db;
          SA_db_vec[SA_qry_cur_index[qry_id]].score = score_considering_length_bottom;
          SA_qry_cur_index[qry_id]++;

	  if(LCP[stop_i-k] < score_bottom)
             score_bottom = LCP[stop_i-k];
          score_considering_length_bottom = score_bottom;
      	}
      }	
      i--;
    }
  }
  //tm += getTime();
  //printf("SA tm2= %lf\n", tm/1e6);

  // Sanity check, can be removed later. For DEBUG purpose only
  //for(i=1;i<total_qry;i++)
    //assert(SA_qry_start_index[i] >= SA_qry_cur_index[i-1]);
  // Sanity check, can be removed later. For DEBUG purpose only

  int SA_max_sz = 0;
  //tm = -getTime();
  for(i=0;i<total_qry;i++) {
    sort(SA_db_vec + SA_qry_start_index[i], SA_db_vec + SA_qry_cur_index[i],SA_sort_db);
    if(SA_qry_cur_index[i] - SA_qry_start_index[i] > SA_max_sz)
      SA_max_sz = SA_qry_cur_index[i] - SA_qry_start_index[i];
  }
  //tm += getTime();
  //printf("SA tm4= %lf, SA_max_sz = %d\n", tm/1e6, SA_max_sz);

  //tm = -getTime();
  // now count frequency
  frequency_for_db most_freq_db;
  //frequency_for_db *freq_db = new frequency_for_db[SA_max_sz];
  assert(SA_max_sz);

  total_freq = 0;
  int prev_db, count, SA_n;
  for ( i=0; i<total_qry; i++) {
    most_freq_db.db_id = -1, most_freq_db.freq = 0;
    prev_db = SA_db_vec[SA_qry_start_index[i]].db_id; count = SA_db_vec[SA_qry_start_index[i]].score;
    SA_n = 0;
    for(j=SA_qry_start_index[i]+1;j<SA_qry_cur_index[i];j++) {
      if (prev_db == SA_db_vec[j].db_id)
  	count = count + SA_db_vec[j].score; 
      else {
	if(count > most_freq_db.freq) {
          most_freq_db.freq = count;
          most_freq_db.db_id = prev_db;
        }
        count = SA_db_vec[j].score;
        prev_db = SA_db_vec[j].db_id;
      }    
    }
    if(count > most_freq_db.freq) {
      most_freq_db.freq = count;
      most_freq_db.db_id = prev_db;
    }

    // Do the processing.
    fprintf(fw2,"%d\t%d\t%d\n", i , most_freq_db.db_id, most_freq_db.freq);
    total_freq += most_freq_db.freq;
  }
  //tm += getTime();
  //printf("SA tm5= %lf\n", tm/1e6);
  printf("total frequency = %d, total query =  %d, avg frequency = %f\n",total_freq, total_qry, float(total_freq)/ total_qry);

  fseek(fw2,0,SEEK_SET);
  fprintf(fw2,"%.0f %d\n",float(total_freq)/ total_qry, TOP_FREQ);
	
  //delete []freq_db;
  return;
}


void initialization(void)
{
  //initialization
  memset(converter,0,128);
  converter['A'] = 0;
  converter['C'] = 1;
  converter['G'] = 2;
  converter['T'] = 3;
  converter['a'] = 0;
  converter['c'] = 1;
  converter['g'] = 2;
  converter['t'] = 3;

  memset(convertStr,'a',128);
  convertStr['A'] = 'a';
  convertStr['C'] = 'c';
  convertStr['G'] = 'g';
  convertStr['T'] = 't';
  convertStr['a'] = 'a';
  convertStr['c'] = 'c';
  convertStr['g'] = 'g';
  convertStr['t'] = 't';
  
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
  
  dna = new unsigned char[SIZE];
  SA = new int[SIZE];
  LCP = new int[SIZE];
  map_input_seq = new int[SIZE][2];
}

int get_input_file()
{
  FILE *f;
  int n;
 
  f = fopen(inputfile1, "r");
  n = input(f,0,-1);
  fclose(f);
  
  start_qry_seq = total_input_seq; //query seq starts from total_input_seq; that is <total_input_seq are db seq  

  //query file 
  f = fopen(inputfile2, "r");
  n = input(f,n,1);
  fclose(f);
  
  dna[n] = '$';
  n++;
  return n;
}

int main(int argc, char **argv) 
{
  int flag, n, total_qry, total_db;
  double my_time = 0;
  FILE *fw2;
  int sTOTAL_DB_SEQ, sTOTAL_QRY_SEQ, sLINE_LEN, sSIZE_SEQ, sSIZE_SEQ_NAME;
  
  ////////////////////// input from top program //////////////////////////
  assert(argc == 11);
  strcpy(inputfile1, argv[1]);
  strcpy(inputfile2, argv[2]);
  strcpy(outputfile, argv[3]);
  sscanf(argv[4], "%d", &SIZE);
  sscanf(argv[5], "%d", &TOP_FREQ);
  sscanf(argv[6], "%d", &sTOTAL_DB_SEQ);
  sscanf(argv[7], "%d", &sTOTAL_QRY_SEQ);
  sscanf(argv[8], "%d", &sLINE_LEN);
  sscanf(argv[9], "%d", &sSIZE_SEQ);
  sscanf(argv[10], "%d", &sSIZE_SEQ_NAME);
  //            //////////////////////  //////////////////////////

  //initialization
  my_time -= getTime();
  
  initialization();

  my_time += getTime();
  printf("Initial memory declearation and initialization %lfs\n", my_time/1e6);
  my_time = 0;
  my_time -= getTime();
  // done initialization

  //db file
  n = get_input_file();

  my_time += getTime();
  printf("Reading and input dna %lfs\n", my_time/1e6);
  my_time = 0;
  my_time -= getTime();
  // input done
    
  
  total_qry = total_input_seq-start_qry_seq;
  total_db = start_qry_seq;

  my_time += getTime();
  printf("Memory declearation and initialization for 2nd step %lfs\n", my_time/1e6);
  my_time = 0;
  my_time -= getTime();
  

  //call SA construction
  flag = sais(dna, SA, LCP, n);
  if (flag != 0)
    {
      printf("SA_IS did not work.\n");
      exit(0);
    }
  else
    printf("SA_IS worked.\n");
  

  my_time += getTime();
  printf("Suffix array construction %lfs\n", my_time/1e6);
  my_time = 0;
  my_time -= getTime();


  fw2 = fopen(outputfile, "w");
  if(fw2 == NULL)
  {
    printf("ERROR: Can't write output file\n");
    exit(0);
  }

  //call hypothesis I
  fprintf(fw2,"111111111111111111111111111111111111111111111111111\n");
  fprintf(fw2,"%s %s trn\n%d %d %d %d %d\n", inputfile2, inputfile1, sTOTAL_DB_SEQ, sTOTAL_QRY_SEQ, sLINE_LEN, sSIZE_SEQ, sSIZE_SEQ_NAME);
  best_exact_match_dna(n, fw2);
  
  printf("Done!!\n");
  return 0;
}
